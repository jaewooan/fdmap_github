module hydrofrac

  use fd, only : limits

  implicit none
  ! -- Second derivative finite difference type
  ! Contains:
  ! NI = number of points in interior stencil
  ! MI = minimum index in interior stencil
  ! PI = maximum index in interior stencil
  ! NBND = number of points in boundary stencil
  ! NBST = number of boundary stencils
  ! H00I = SBP H-norm weight at boundary grid point
  ! DI = interior derivative coefficients
  ! DL/DR = left/right boundary derivative coefficients
  ! D1L/D1R = first derivative boundary coefficients
  ! D2 = full second derivative matrix including boundary treatment (used by implicit solve)
  ! H1Li/H1Ri = inverted diagonal norm coefficients near boundaries (needed by penalty terms of the
  ! form D^T)
  ! HL/HR = diagonal norm entries near boundary (normalized so H=1 in interior)
  ! FDmethod = finite difference method for second derivative
  type :: fd2_type
    integer :: nI,mI,pI,nbnd,nbst,nD1
    real :: H00i
    real,dimension(:),allocatable :: DI,DmI,DpI,HL,HR,D1L,D1R,HLi,HRi
    real,dimension(:,:),allocatable :: DL,DR,D2
    character(4) :: FDmethod
  end type


  ! USE_HF = flag indicating hydraulic fracture
  ! BNDM,BNDP = process is responsible for minus/plus boundary 
  ! DIRECTION = direction of fracture
  ! COUPLED = full fluid-solid coupling (F=neglect crack opening in mass balance)
  ! LINEARIZED_WALLS = linearize wall opening term in mass balance
  ! SOURCE_TERM = add source term to mass balance equation
  ! OPERATOR_SPLITTING = use operator splitting to handle stiffness
  ! RHO0 = (initial) density
  ! K0 = bulk modulus
  ! C0 = sound speed
  ! mu = dynamic viscosity
  ! N = number of points resolving width
  ! WM,WP = position of minus (lower) and upper walls (wp-wm=width)
  ! WM0,WP0 = initial values of wm,wp
  ! U = width-averaged velocity
  ! V = velocity
  ! P = pressure
  ! DU = rate of change of width-averaged velocity
  ! DV = rate of change of velocity
  ! DP = rate of change of pressure
  ! DWM,DWP = rate of wall position
  ! Hw = width-averaging operator 
  ! SATM,SATP = SAT penalty weights (~wave speed / grid spacing)
  ! SATYM, SATYP = SAT penalty terms used by viscous terms 
  ! L = limits type (indices for FD operations)
  ! XSOURCE,YSOURCE = source position
  ! TSOURCE = source duration
  ! ASOURCE = source amplitude
  ! WSOURCE = source width (in space)
  ! fd2_type  = finite difference second derivative operator

  type :: hf_type
     logical :: use_HF,coupled,linearized_walls,source_term,operator_splitting,bndm,bndp,inviscid,slope
     integer :: n
     character(1) :: direction
     real :: rho0,K0,mu,c0,h,SATm,SATp,xsource,ysource,tsource,Asource,wsource
     real,dimension(:),allocatable :: wm,wp,wm0,wp0,dwm0dx,dwp0dx,u,p,Du,Dp,Dwm,Dwp,Hw,hy,SATyp,SATym
     real,dimension(:,:),allocatable :: v,Dv
     type(limits) :: L
     type(fd2_type) :: fd2
  end type hf_type


contains


  subroutine init_hydrofrac(iface,HF,m,p,x,y,hmin,hmax,refine,input,echo,skip,direction,mg,pg, &
       process_m,process_p,comm_m,comm_p,array)

    use fd_coeff, only : H00i
    use mpi_routines, only : is_master
    use io, only : warning,error,write_matlab,seek_to_string  

    implicit none

    integer,intent(in) :: iface,m,p,input,echo,mg,pg,comm_m,comm_p,array
    type(hf_type),intent(out) :: HF
    real,intent(in) :: hmin,hmax,refine
    real,dimension(m:p),intent(in) :: x,y
    logical,intent(in) :: skip,process_m,process_p
    character(1),intent(in) :: direction

    integer :: stat
    character(256) :: HFstr
    character(256) :: str

    logical :: coupled,linearized_walls,source_term,operator_splitting,inviscid
    integer :: n,i,j
    real :: rho0,K0,mu,w0,xsource,ysource,tsource,Asource,wsource
    character(256) :: geom_file,initial_conds_file
    character(4) :: FDmethod
    integer,parameter :: nb=3 ! number of additional boundary (ghost) points needed for FD operators

    namelist /hydrofrac_list/ rho0,K0,mu,w0,xsource,ysource,tsource,Asource,wsource, &
         n,coupled,linearized_walls,source_term,operator_splitting,geom_file,initial_conds_file, &
         inviscid,FDmethod

    ! defaults
       
    rho0 = 0d0
    K0 = 1d40
    w0 = 0d0
    mu = 0d0
    xsource = 0d0
    ysource = 0d0
    tsource = 1d0
    Asource = 1d0
    wsource = 1d0
    n = 1
    coupled = .true.
    linearized_walls = .false.
    source_term = .true.
    operator_splitting = .true.
    inviscid = .true.
    geom_file = ''
    initial_conds_file = ''
    FDmethod = 'SBP6'

    ! read in hydraulic fracture parameters

    write(str,'(a,i0,a)') '!---IFACE',iface,'---'
    call seek_to_string(input,str)
    read(input,nml=hydrofrac_list,iostat=stat)
    if (stat>0) call error('Error in hydrofrac_list','init_hydrofrac')

    HF%use_HF = (rho0/=0d0) ! flag: T if using hydraulic fracture routines

    if (.not.HF%use_HF) return

    select case(direction)
    case('x')
       HF%direction = 'y'
    case('y')
       HF%direction = 'x'
    end select

    ! input parameters

    HF%rho0 = rho0
    HF%K0 = K0
    HF%mu = mu
    HF%coupled = coupled
    HF%linearized_walls = linearized_walls
    HF%operator_splitting = operator_splitting
    HF%inviscid =  inviscid
    HF%source_term = source_term
    HF%xsource = xsource
    HF%ysource = ysource
    HF%tsource = tsource
    HF%Asource = Asource
    HF%wsource = wsource
    HF%n = ceiling(dble(n-1)*refine)+1
    HF%h = 0.5d0*(hmin+hmax)
    if (abs(hmin-hmax)>1d-12.and.is_master) &
         call warning('Possible nonuniform grid spacing: hydrofrac routines may not work','init_hydrofrac')

    ! other parameters

    HF%c0 = sqrt(HF%K0/HF%rho0)
    HF%SATm = H00i*HF%c0/HF%h ! SAT penalty weight, minus side
    HF%SATp = H00i*HF%c0/HF%h ! SAT penalty weight, plus  side



    ! output hydraulic fracture parameters
    
    if (is_master) then
       write(HFstr,'(a,i0,a)') 'HF{',iface,'}'
       call write_matlab(echo,'rho0',HF%rho0,HFstr)
       call write_matlab(echo,'K0',HF%K0,HFstr)
       call write_matlab(echo,'c0',HF%c0,HFstr)
       call write_matlab(echo,'w0',w0,HFstr)
       call write_matlab(echo,'coupled',HF%coupled,HFstr)
       call write_matlab(echo,'linearized_walls',HF%linearized_walls,HFstr)
       call write_matlab(echo,'operator_splitting',HF%operator_splitting,HFstr)
       call write_matlab(echo,'inviscid',HF%inviscid,HFstr)
       call write_matlab(echo,'source_term',HF%source_term,HFstr)
       if (HF%source_term) then
          call write_matlab(echo,'xsource',HF%xsource,HFstr)
          call write_matlab(echo,'ysource',HF%ysource,HFstr)
          call write_matlab(echo,'tsource',HF%tsource,HFstr)
          call write_matlab(echo,'Asource',HF%Asource,HFstr)
          call write_matlab(echo,'wsource',HF%wsource,HFstr)
       end if
    end if

    ! return if not needed

    if (skip) return

    ! initialize grid and fields

    HF%bndm = (m==mg) ! check if process handles minus boundary
    HF%bndp = (p==pg) ! check if process handles plus  boundary
    HF%L = limits(nb,m,p,mg,pg,m-nb,p+nb,HF%bndm,HF%bndp) ! derived type used for FD
    ! where L%m:L%p  = range of indices handled by process
    ! and  L%mb:L%pb = range of indices including ghost points


    ! allocate arrays and initialize with unreasonable values to facilitate debugging

    allocate( &
         HF%wm (HF%L%mb :HF%L%pb ),HF%wp (HF%L%mb :HF%L%pb ), &
         HF%wm0(HF%L%mb :HF%L%pb ),HF%wp0(HF%L%mb :HF%L%pb ), &
         HF%Dwm(HF%L%m :HF%L%p ),HF%Dwp(HF%L%m :HF%L%p ), &
         HF%u  (HF%L%mb:HF%L%pb),HF%p  (HF%L%mb:HF%L%pb), &
         HF%Du (HF%L%m :HF%L%p ),HF%Dp (HF%L%m :HF%L%p ), &
         HF%v (HF%n,HF%L%m:HF%L%p), &
         HF%Dv(HF%n,HF%L%m:HF%L%p), &
         HF%Hw(HF%n), &
         HF%hy(HF%L%mb:HF%L%pb), &
         HF%SATyp(4),HF%SATym(4) &
            )

    HF%wm0 = 1d40
    HF%wp0 = 1d40
    HF%wm  = 1d40
    HF%wp  = 1d40
    HF%u   = 1d40
    HF%v   = 1d40
    HF%p   = 1d40
    HF%Dwm = 1d40
    HF%Dwp = 1d40
    HF%Du  = 1d40
    HF%Dv  = 1d40
    HF%Dp  = 1d40

    ! initial conditions on wall position, velocity, and pressure perturbation
    ! (spatially uniform values can be set from input file parameters)

    HF%wm0 = -w0
    HF%wp0 =  w0

    ! Slope profile (only used if profiles are read from file)
    HF%slope = .false.

    HF%u = 0d0
    HF%v = 0d0
    HF%p = 0d0

     ! Load non planar profiles

     if (geom_file /= '' ) then

        allocate(HF%dwm0dx(HF%L%m:HF%L%p),HF%dwp0dx(HF%L%m:HF%L%p))

        ! both sides read file (so process may read file twice)
        if (process_m) call read_hydrofrac_geom(HF,geom_file,comm_m,array)
        if (process_p) call read_hydrofrac_geom(HF,geom_file,comm_p,array)

        HF%slope = .true.

    end if


     
    ! override uniform initial conditions with those from file, if desired

     if (initial_conds_file /= '' ) then
        ! both sides read file (so process may read file twice)
        if (process_m) call read_hydrofrac_initial_conds(HF,initial_conds_file,comm_m,array)
        if (process_p) call read_hydrofrac_initial_conds(HF,initial_conds_file,comm_p,array)
     end if

    ! current (=initial) wall positions

    HF%wm = HF%wm0
    HF%wp = HF%wp0

    ! Compute Width-averaging norm operator
    call Hnorm(HF%Hw)

    if(HF%inviscid) return

    ! Compute grid spacings in the cross-section direction (y)
    HF%hy = (HF%wp0 - HF%wm0)/(HF%n - 1)
    

    ! Set SAT penalty weights for viscous terms
    ! S(1)*(v - g) + S(2)*(Dv - g) + S(3)*D^T*(v - g) + S(4)*D^T*(Dv - g)

    ! Impose velocities only
    HF%SATym(:) =   HF%mu/HF%rho0*(/  0d0, 0d0, 1d0, 0d0  /)
    HF%SATyp(:) = - HF%mu/HF%rho0*(/  0d0, 0d0, 1d0, 0d0  /) 
    
    ! Initialize second derivative SBP operator
    if( .not. inviscid) then
        HF%fd2%FDmethod = FDmethod
        call init_fd2(FDmethod,HF%fd2)
        if(operator_splitting) call init_fd2_full(HF,HF%n)
    end if 
    

  end subroutine init_hydrofrac


   subroutine read_hydrofrac_geom(HF,geom_file,comm,array)

     use io, only : file_distributed,open_file_distributed, &
          read_file_distributed,close_file_distributed
     use mpi_routines, only : pw

     implicit none

     type(hf_type),intent(inout) :: HF
     character(*),intent(in) :: geom_file
     integer,intent(in) :: comm,array

     integer :: j
     type(file_distributed) :: fh

     call open_file_distributed(fh,geom_file,'read',comm,array,pw)

     call read_file_distributed(fh,HF%wm0(HF%L%m:HF%L%p))
     call read_file_distributed(fh,HF%wp0(HF%L%m:HF%L%p))
     call read_file_distributed(fh,HF%dwm0dx(HF%L%m:HF%L%p))
     call read_file_distributed(fh,HF%dwp0dx(HF%L%m:HF%L%p))

     call close_file_distributed(fh)

   end subroutine read_hydrofrac_geom
  
  subroutine read_hydrofrac_initial_conds(HF,initial_conds_file,comm,array)

    use io, only : file_distributed,open_file_distributed, &
         read_file_distributed,close_file_distributed
    use mpi_routines, only : pw

    implicit none

    type(hf_type),intent(inout) :: HF
    character(*),intent(in) :: initial_conds_file
    integer,intent(in) :: comm,array

    integer :: j
    type(file_distributed) :: fh

    call open_file_distributed(fh,initial_conds_file,'read',comm,array,pw)

     do j = 1,HF%n
        call read_file_distributed(fh,HF%v(j,HF%L%m:HF%L%p))
     end do

    call read_file_distributed(fh,HF%u  (HF%L%m:HF%L%p))
    call read_file_distributed(fh,HF%p  (HF%L%m:HF%L%p))

    call close_file_distributed(fh)

  end subroutine read_hydrofrac_initial_conds


  subroutine destroy_hydrofrac(HF)

    implicit none

    type(hf_type),intent(inout) :: HF

    if (allocated(HF%wm0)) deallocate(HF%wm0)
    if (allocated(HF%wp0)) deallocate(HF%wp0)
    if (allocated(HF%dwm0dx)) deallocate(HF%dwm0dx)
    if (allocated(HF%dwp0dx)) deallocate(HF%dwp0dx)
    if (allocated(HF%wm )) deallocate(HF%wm )
    if (allocated(HF%wp )) deallocate(HF%wp )
    if (allocated(HF%u  )) deallocate(HF%u  )
    if (allocated(HF%v  )) deallocate(HF%v  )
    if (allocated(HF%p  )) deallocate(HF%p  )
    if (allocated(HF%Dwm)) deallocate(HF%Dwm)
    if (allocated(HF%Dwp)) deallocate(HF%Dwp)
    if (allocated(HF%Du )) deallocate(HF%Du )
    if (allocated(HF%Dv )) deallocate(HF%Dv )
    if (allocated(HF%Dp )) deallocate(HF%Dp )
    if (allocated(HF%Hw )) deallocate(HF%Hw )
    if (allocated(HF%hy )) deallocate(HF%hy )
    if (allocated(HF%SATyp )) deallocate(HF%SATyp )
    if (allocated(HF%SATym )) deallocate(HF%SATym )
    if (allocated(HF%fd2%D2 )) deallocate(HF%fd2%D2)
    if (allocated(HF%fd2%DI )) deallocate(HF%fd2%DI)
    if (allocated(HF%fd2%DmI )) deallocate(HF%fd2%DmI)
    if (allocated(HF%fd2%DpI )) deallocate(HF%fd2%DpI)
    if (allocated(HF%fd2%D1L )) deallocate(HF%fd2%D1L)
    if (allocated(HF%fd2%D1R )) deallocate(HF%fd2%D1R)
    if (allocated(HF%fd2%HL )) deallocate(HF%fd2%HL)
    if (allocated(HF%fd2%HR )) deallocate(HF%fd2%HR)
    if (allocated(HF%fd2%HLi )) deallocate(HF%fd2%HLi)
    if (allocated(HF%fd2%HRi )) deallocate(HF%fd2%HRi)
    if (allocated(HF%fd2%DL )) deallocate(HF%fd2%DL)
    if (allocated(HF%fd2%DR )) deallocate(HF%fd2%DR)

  end subroutine destroy_hydrofrac


  subroutine checkpoint_hydrofrac(fh,operation,HF)

    use io, only : file_distributed, &
         read_file_distributed,write_file_distributed

    implicit none

    type(file_distributed),intent(in) :: fh
    character(*),intent(in) :: operation
    type(hf_type),intent(inout) :: HF

    integer :: j,m,p

    ! fields read/written to same file as iface fields,
    ! routine is only called by io_process

    if (.not.HF%use_HF) return

    m = HF%L%m
    p = HF%L%p

    ! save all fields (but not rates) when checkpointing

    select case(operation)
    case('read')
       call read_file_distributed(fh,HF%wm(m:p))
       call read_file_distributed(fh,HF%wp(m:p))
       call read_file_distributed(fh,HF%u (m:p))
       call read_file_distributed(fh,HF%p (m:p))
       do j = 1,HF%n
          call read_file_distributed(fh,HF%v(j,m:p))
       end do
    case('write')
       call write_file_distributed(fh,HF%wm(m:p))
       call write_file_distributed(fh,HF%wp(m:p))
       call write_file_distributed(fh,HF%u (m:p))
       call write_file_distributed(fh,HF%p (m:p))
       do j = 1,HF%n
          call write_file_distributed(fh,HF%v(j,m:p))
       end do
    end select

  end subroutine checkpoint_hydrofrac


  subroutine scale_rates_hydrofrac(HF,A)

    implicit none

    type(hf_type),intent(inout) :: HF
    real,intent(in) :: A

    if (.not.HF%use_HF) return

    HF%Dwm = A*HF%Dwm
    HF%Dwp = A*HF%Dwp
    HF%Du  = A*HF%Du
    HF%Dv  = A*HF%Dv
    HF%Dp  = A*HF%Dp

  end subroutine scale_rates_hydrofrac


  subroutine update_fields_hydrofrac(HF,dt)

    implicit none

    type(hf_type),intent(inout) :: HF
    real,intent(in) :: dt

    integer :: i

    if (.not.HF%use_HF) return

    do i = HF%L%m,HF%L%p
       HF%wm(i) = HF%wm(i)+dt*HF%Dwm(i)
       HF%wp(i) = HF%wp(i)+dt*HF%Dwp(i)
       HF%p (i) = HF%p (i)+dt*HF%Dp (i)
       HF%v(:,i) = HF%v(:,i)+dt*HF%Dv(:,i)
    end do

  end subroutine update_fields_hydrofrac


  subroutine set_rates_hydrofrac(HF,C,m,p,phip,vnm,vnp,vtm,vtp,sntm,sntp,x,y,t)

    use mpi_routines2d, only : cartesian
    use fd, only : diff

    implicit none

    type(hf_type),intent(inout) :: HF
    type(cartesian),intent(in) :: C
    integer,intent(in) :: m,p
    real,intent(in) :: phip(m:p),vnm(m:p),vnp(m:p),vtm(m:p),vtp(m:p), &
                       sntm(m:p),sntp(m:p),x(m:p),y(m:p),t

    integer :: i
    real :: uhat,phat, Z
    real,dimension(:),allocatable :: dudx,dpdx,b

    if (.not.HF%use_HF) return
    
    ! Width-average v to obtain u
    call width_average(HF)

    ! populate ghost points of neighboring processors

    call share_hydrofrac(HF,C)

    ! allocate auxiliary arrays

    allocate(dudx(HF%L%m:HF%L%p),dpdx(HF%L%m:HF%L%p),b(HF%L%mb:HF%L%pb))


    ! SBP differentiation of velocity and pressure
    ! (this could certainly be improved for compatibility with external mesh)

    b = HF%wp0 - HF%wm0
    call diff(HF%L,HF%u*b,dudx)
    call diff(HF%L,HF%p,dpdx)
    dudx = dudx/HF%h
    dpdx = dpdx/HF%h

    ! set rates from mass and momentum balance equations,
    ! starting with linearized acoustics (rigid walls)


    do i = HF%L%m,HF%L%p
       HF%Dv(:,i) = HF%Dv(:,i)-dpdx(i)/HF%rho0
       HF%Dp(i) = HF%Dp(i)-dudx(i)*HF%K0/b(i)
    end do
    
    do i = HF%L%m,HF%L%p
       HF%Dwp(i) = HF%Dwp(i) + vnp(i)
       HF%Dwm(i) = HF%Dwm(i) + vnm(i)
    end do


    ! and source terms (explosion source only appears in mass balance)

    if (HF%source_term) then
       do i = HF%L%m,HF%L%p
          HF%Dp(i) = HF%Dp(i)+HF%Asource* &
               exp(-0.5d0*((x(i)-HF%xsource)**2+(y(i)-HF%ysource)**2)/HF%wsource**2)* &
               (t/HF%tsource)*exp(-t/HF%tsource)
       end do
    end if


    call set_rates_viscous(HF,m,p,vtm,vtp,sntm,sntp)

    ! add penalty terms to enforce BC with SAT method

    if (HF%bndm) then ! check if process handles minus boundary
       i = HF%L%mg
       uhat = 0d0 ! zero fluid velocity at crack tip
       phat = HF%p(i)-HF%rho0*HF%c0*HF%u(i) ! set by preserving characteristic variable into fluid
       HF%Du(i) = HF%Du(i)-HF%SATm*(HF%u(i)-uhat)
       HF%Dp(i) = HF%Dp(i)-HF%SATm*(HF%p(i)-phat)
    end if

    if (HF%bndp) then ! check if process handles plus  boundary
       i = HF%L%pg
       uhat = 0d0 ! zero fluid velocity at crack tip
       phat = HF%p(i)+HF%rho0*HF%c0*HF%u(i) ! set by preserving characteristic variable into fluid
       HF%Du(i) = HF%Du(i)-HF%SATp*(HF%u(i)-uhat)
       HF%Dp(i) = HF%Dp(i)-HF%SATp*(HF%p(i)-phat)
    end if
    
    if (.not. HF%coupled) return

    ! and then adding crack opening/closing term in mass balance

    if (HF%linearized_walls .and. .not. HF%operator_splitting) then
        ! Couple to solid without geometric correction
        do i = HF%L%m,HF%L%p
           HF%Dp(i) = HF%Dp(i)-HF%K0*(vnp(i)-vnm(i))/(HF%wp0(i)-HF%wm0(i))
        end do

        ! Add geometric correction
        if(HF%slope) then
           do i = HF%L%m,HF%L%p
              HF%Dp(i) = HF%Dp(i)+HF%K0*(vtp(i)*HF%dwp0dx(i) - vtm(i)*HF%dwm0dx(i))/(HF%wp0(i)-HF%wm0(i))
           end do
        end if
    end if

    ! TODO: Implement this properly
    if (.not. HF%linearized_walls) then
       if (HF%operator_splitting) then
          do i = HF%L%m,HF%L%p
             HF%Dp(i) = HF%Dp(i)-HF%K0*phip(i)/(HF%wp(i)-HF%wm(i)) ! *INCORRECT*
          end do
       else
          do i = HF%L%m,HF%L%p
             HF%Dp(i) = HF%Dp(i)-HF%K0*(vnp(i)-vnm(i))/(HF%wp(i)-HF%wm(i))
          end do
       end if
    end if

    ! deallocate temporary arrays

    deallocate(dudx,dpdx,b)

  end subroutine set_rates_hydrofrac


  subroutine fluid_stresses(HF,i,p,taum,taup,vtm,vtp)

    implicit none

    type(HF_type),intent(in) :: HF
    integer,intent(in) :: i
    real,intent(out) :: p,taum,taup,vtm,vtp

    ! fluid pressure

    p = HF%p(i)

    ! Tangential velocities
    vtm = HF%v(1,i)
    vtp = HF%v(HF%n,i)

    ! shear stress on bottom (taum) and top (taup) walls,
    ! evaluate these as mu*dv/dy for viscous fluid with fluid velocity v


    ! Inviscid case
    taum = 0d0
    taup = 0d0

    if(HF%inviscid) return

    ! Viscous case
    taum = HF%mu*diff_bnd_m(HF%v(:,i),HF%fd2)/HF%hy(i)
    taup = HF%mu*diff_bnd_p(HF%v(:,i),HF%fd2)/HF%hy(i)

    ! Add geometric correction
    if(.not. HF%slope) return

    taum = taum + HF%p(i)*HF%dwm0dx(i)
    taup = taup + HF%p(i)*HF%dwp0dx(i)

  end subroutine fluid_stresses

  subroutine update_fields_hydrofrac_implicit(HF,m,p,Zm,Zp,phip,vtp,vtm,dt)
  
    implicit none

    type(HF_type),intent(inout) :: HF
    integer,intent(in) :: m,p
    real,dimension(m:p),intent(in) :: Zm,Zp,phip, & ! P-wave impedances on minus,plus sides
                                      vtp,vtm ! tangential solid velocities on
                                              ! plus, minus sides
    real,intent(in) :: dt

    integer :: i,j,info
    real :: Z,Zi
    real,dimension(HF%n) :: b ! Right hand side vector containing boundary terms 
                              ! and solution at previous time step
    real,dimension(HF%n,HF%n) :: A ! matrix to LU decompose
    integer :: ipiv(HF%n)

    if (.not.HF%operator_splitting) return ! return if fully explicit

    ! implicit update for pressure
    
    if (HF%linearized_walls .and. .not. HF%slope) then
       do i = m,p
          Zi = (Zp(i) + Zm(i) )/(Zp(i)*Zm(i)) ! combine impedances of two sides
          !HF%p(i) = HF%p(i)*(1d0-dt*HF%K0*Z/(HF%wp0(i)-HF%wm0(i))) ! forward Euler
          HF%p(i) = (HF%p(i) -  dt*HF%K0*phip(i)/(HF%wp0(i)-HF%wm0(i)) ) &
                    /(1d0+(dt*HF%K0*Zi)/(HF%wp0(i)-HF%wm0(i))) ! backward Euler
              ! HF%Dp(i) = HF%Dp(i)+HF%K0*(vtp(i)*HF%dwp0dx(i) - vtm(i)*HF%dwm0dx(i))/(HF%wp0(i)-HF%wm0(i))
       end do
    end if
    ! Implicit update with geometric correction

    if (HF%linearized_walls .and. HF%slope) then
       do i = m,p
          Zi = (Zp(i) + Zm(i) )/(Zp(i)*Zm(i)) ! combine impedances of two sides
          !HF%p(i) = HF%p(i)*(1d0-dt*HF%K0*Z/(HF%wp0(i)-HF%wm0(i))) ! forward Euler
          HF%p(i) = (HF%p(i) -  dt*HF%K0*& 
                        (       &
                            phip(i) - vtp(i)*HF%dwp0dx(i) + vtm(i)*HF%dwm0dx(i) &
                        )       &              
                        /(HF%wp0(i)-HF%wm0(i)) &
                    ) &
                    /(1d0+(dt*HF%K0*Zi)/(HF%wp0(i)-HF%wm0(i))) ! backward Euler
               ! HF%Dp(i) = HF%Dp(i)+HF%K0*()/(HF%wp0(i)-HF%wm0(i))
       end do
    end if
    
    if (.not. HF%linearized_walls) then
       do i = m,p
          Z = 1d0/(1d0/Zm(i)+1d0/Zp(i)) ! combine impedances of two sides
          ! at what time should we use wp-wm?
          HF%p(i) = HF%p(i)/(1d0+dt*HF%K0/(Z*(HF%wp (i)-HF%wm (i)))) ! backward Euler
       end do
    end if

    ! implicit update for viscous diffusion terms

    if(HF%inviscid) return ! Return if inviscid

    ! Solve (I - dt/h^2*M)v^N+1 = v^N + dt*b 
    ! M contains D2 operator and SAT penalty terms

    do i = m,p
       b = HF%v(:,i)
       call diff_T_bnd_p(b,-dt*HF%SATyp(3)*vtp(i)/HF%hy(i)**2,HF%fd2)
       call diff_T_bnd_m(b,-dt*HF%SATym(3)*vtm(i)/HF%hy(i)**2,HF%fd2)
       A =  - dt*HF%fd2%D2/HF%hy(i)**2
       ! Add identity matrix
       do j=1,HF%n
       A(j,j) = A(j,j) + 1
       end do
       call dgesv(HF%n,1,A,HF%n,ipiv,b,HF%n,info)
       HF%v(:,i) = b
    end do

  end subroutine update_fields_hydrofrac_implicit


  subroutine share_hydrofrac(HF,C)

    use mpi_routines2d, only : cartesian,populate_ghost_cells

    implicit none

    type(HF_type),intent(inout) :: HF
    type(cartesian),intent(in) :: C

    ! communicate fields to ghost points

    call populate_ghost_cells(C,HF%L%m,HF%L%p,HF%L%nb, &
         HF%bndm,HF%bndp,HF%u,HF%direction)
    call populate_ghost_cells(C,HF%L%m,HF%L%p,HF%L%nb, &
         HF%bndm,HF%bndp,HF%p,HF%direction)

  end subroutine share_hydrofrac
  
  subroutine Hnorm(Hw)


      use fd_coeff, only : HL, HR, nbst 

      implicit none

      real,intent(out),dimension(:) :: Hw
      integer :: n


      Hw = 1d0 ! interior values
      n = size(Hw,1)

      !todo: add assertions to ensure that n is large enough to hold the norm

      Hw(1:nbst) = HL
      Hw((n-nbst+1):n) = HR

  end subroutine Hnorm

  ! Apply the width-averaging operator along each cross-section of the hydrofrac
  ! layer
  ! The velocity v is width-averaged to u
  subroutine width_average(HF)

      implicit none

      type(hf_type), intent(inout) :: HF
      integer :: i

       do i = HF%L%m,HF%L%p
          HF%u(i) = HF%u(i) + dot_product(HF%v(:,i),HF%Hw)
       end do
       HF%u = HF%u/(HF%n - 1)

  end subroutine
  
  subroutine set_rates_viscous(HF,m,p,vtm,vtp,sntm,sntp)

      implicit none

      type(hf_type),intent(inout) :: HF
      integer,intent(in) :: m,p
      real,intent(in) :: vtm(m:p),vtp(m:p),sntm(m:p),sntp(m:p)

      real,dimension(:),allocatable :: v_yy
      real :: gp,gm ! Boundary conditions of the form: gp = SAT_WEIGHT*b.c/h 
      integer :: i

      if(HF%inviscid) return
      if(HF%operator_splitting) return

      allocate(v_yy(HF%n))
      v_yy = 0d0
    
      ! and viscous terms (explicit time-stepping)
      ! SAT terms impose velocities in the fluid

      do i = HF%L%m,HF%L%p
         call second_derivative(HF%fd2,HF%v(:,i),v_yy)
         v_yy = v_yy/(HF%hy(i)**2)
         HF%Dv(:,i) = HF%Dv(:,i)+HF%mu/HF%rho0*v_yy

         ! SAT terms, minus boundary
         ! Penalty terms not used
         !HF%Dv(1,i) = HF%Dv(1,i) + HF%SATym(1,i)*HF%fd2%H00i*HF%v(1,i)
         !HF%Dv(1,i) = HF%Dv(1,i) + HF%SATym(2,i)*HF%fd2%H00i*(diff_bnd_m(HF%v(:,i),HF%fd2)/HF%hy(i)) 
         gm = HF%SATym(3)*( HF%v(1,i)  - vtp(i) )/HF%hy(i)**2
         call diff_T_bnd_m(HF%Dv(:,i),gm,HF%fd2)

         ! SAT terms, plus boundary
         ! Penalty terms not used
         !HF%Dv(HF%n,i) = HF%Dv(HF%n,i) + HF%SATyp(1,i)*HF%fd2%H00i*HF%v(i,HF%n)
         !HF%Dv(HF%n,i) = HF%Dv(HF%n,i) + HF%SATyp(2,i)*HF%fd2%H00i*(diff_bnd_p(HF%v(:,i),HF%fd2)/HF%hy(i))
         gp = HF%SATyp(3)*( HF%v(HF%n,i) - vtm(i) )/HF%hy(i)**2
         call diff_T_bnd_p(HF%Dv(:,i),gp,HF%fd2)
      end do

      deallocate(v_yy)

  end subroutine 
  
  ! Compute the second derivative of v and store the result in v_yy
  subroutine second_derivative(fd2,v,v_yy)

      implicit none

      type(fd2_type),intent(in) :: fd2
      real,dimension(:),intent(in) :: v
      real,dimension(:),intent(out) :: v_yy

      integer :: i,n

      n = size(v,1)
      
      ! Interior
      do i = fd2%nbst+1,n-fd2%nbst
         v_yy(i) = dot_product(fd2%DI,v(i+fd2%mI:i+fd2%pI))
      end do

      ! Left boundary
      do i = 0,fd2%nbst-1
         v_yy(1 + i) = dot_product( fd2%DL(:, i),v(1:fd2%nbnd) ) 
      end do
    
      ! Right boundary
      do i = 0,fd2%nbst-1         
         v_yy(n-i) = dot_product( fd2%DR(:,-i),v((n-fd2%nbnd+1):n) )
      end do
    

  end subroutine

  ! Differentiate a vector u on minus boundary (used by penalty terms)
  pure function diff_bnd_m(u,fd2) result(v)

      implicit none

      real,dimension(:),intent(in) :: u
      type(fd2_type),intent(in) :: fd2
      real :: v
      integer :: n

      n = size(u,1)

      v = dot_product(u(1:fd2%nD1),fd2%D1L)

  end function

  ! Differentiate a vector u on plus boundary (used by penalty terms)
  pure function diff_bnd_p(u,fd2) result(v)


      implicit none

      real,dimension(:),intent(in) :: u
      type(fd2_type),intent(in) :: fd2
      real :: v
      integer :: n

      n = size(u,1)

      v = dot_product(u((n-fd2%nD1+1):n),fd2%D1R)

  end function
  
  
  ! Compute Dv = Dv + D^T*g on minus boundary (used by penalty terms)
  subroutine diff_T_bnd_m(Dv,g,fd2)

      implicit none

      real,dimension(:),intent(inout) :: Dv
      real,intent(in) :: g
      type(fd2_type),intent(in) :: fd2
      
      integer :: n

      n = size(Dv,1)
      
      Dv(1:fd2%nD1) = Dv(1:fd2%nD1) + fd2%HLi*fd2%D1L*g

  end subroutine

  ! Compute Dv = Dv + D^T*g on plus boundary (used by penalty terms)
  subroutine diff_T_bnd_p(Dv,g,fd2)
      

      implicit none

      real,dimension(:),intent(inout) :: Dv
      real,intent(in) :: g
      type(fd2_type),intent(in) :: fd2
      
      integer :: n

      n = size(Dv,1)
      
      Dv((n-fd2%nD1+1):n) = Dv((n-fd2%nD1+1):n) + fd2%HRi*fd2%D1R*g

  end subroutine
  
  ! Initializes stencils for compact, second derivatives and 
  ! compatible first derivative boundary stencils
  subroutine init_fd2(FDmethod,fd2)

      use io, only : error

      implicit none

      character(*),intent(in) :: FDmethod

      integer :: i
      type(fd2_type), intent(out) :: fd2

      select case(FDmethod)
      case default
        call error('Invalid second derivative FD method','init_fd2')
      case('SBP2')
        fd2%nI = 3
        fd2%nbnd = 3
        fd2%nbst = 1
        fd2%H00i = 2d0
        fd2%nD1  = 3
      case('SBP4')
        fd2%nI = 5
        fd2%nbnd = 6
        fd2%nbst = 4
        fd2%H00i = 2.823529411764706d0
        fd2%nD1  = 4
      case('SBP6')
        fd2%nI = 7
        fd2%nbnd = 9
        fd2%nbst = 6
        fd2%H00i = 3.165067037878233d0
        fd2%nD1  = 5
      end select

    fd2%mI = -(fd2%nI-1)/2
    fd2%pI =  (fd2%nI-1)/2

    allocate(fd2%DI(fd2%mI:fd2%pI))
    allocate(fd2%DL (0:fd2%nbnd-1,0:fd2%nbst-1),fd2%DR (-(fd2%nbnd-1):0,-(fd2%nbst-1):0))
    allocate(fd2%D1L(fd2%nD1),fd2%D1R(fd2%nD1))
    allocate(fd2%HLi(fd2%nD1),fd2%HRi(fd2%nD1))

    select case(FDmethod)

    case('SBP2')

        ! interior

        fd2%DI = (/ 1.0d0, -2.0d0, 1.0d0 /)

        ! left boundary

        fd2%DL(:,0) = fd2%DI

        ! First derivative for the left boundary

        fd2%D1L = (/ -3/2d0, 2d0, -1/2d0 /)  

        ! Norm for the left boundary

        fd2%HLi = (/ 1/2d0, 1d0, 1d0 /)

    case('SBP4')

        ! interior 

        fd2%DI = (/ -0.083333333333333d0,  1.333333333333333d0, -2.500000000000000d0,       &
                     1.333333333333333d0, -0.083333333333333d0                              /)

        ! left boundary

        fd2%DL(:,0) = (/ 2d0, -5d0, 4d0, -1d0, 0d0, 0d0 /)
        fd2%DL(:,1) = (/ 1d0, -2d0, 1d0,  0d0, 0d0, 0d0 /)
        fd2%DL(:,2) = (/ -0.093023255813953d0,  1.372093023255814d0, -2.558139534883721d0,  &
                          1.372093023255814d0, -0.093023255813953d0,  0.000000000000000d0   /)
        fd2%DL(:,3) = (/ -0.020408163265306d0,  0.000000000000000d0,  1.204081632653061d0,  &
                         -2.408163265306122d0,  1.306122448979592d0, -0.081632653061224d0   /)

        ! First derivative for the left boundary

        fd2%D1L = (/ -11/6d0, 3d0, -3/2d0, 1/3d0 /)  

        ! Norm for the left boundary

        fd2%HLi = (/ 0.354166666666667d0,1.229166666666667d0,0.895833333333333d0,1.020833333333333d0 /)

    case('SBP6')
        
        ! interior

        fd2%DI = (/  0.011111111111111d0,  -0.150000000000000d0,   1.500000000000000d0,     & 
                    -2.722222222222222d0,   1.500000000000000d0,  -0.150000000000000d0,     &  
                     0.011111111111111d0                                                    /)

        fd2%DL(:,0) = (/  2.788238454587638d0,  -8.024525606271522d0,   8.215717879209711d0, &  
                         -3.382384545876377d0,   0.274525606271522d0,   0.128428212079029d0, &
                          0.000000000000000d0,   0.000000000000000d0,   0.000000000000000d0  /)
        fd2%DL(:,1) = (/  1.053412969283277d0,  -2.350398179749716d0,   1.867463026166098d0, &  
                         -1.034129692832765d0,   0.600398179749716d0,  -0.136746302616610d0, &
                          0.000000000000000d0,   0.000000000000000d0,   0.000000000000000d0  /)
        fd2%DL(:,2) = (/ -0.644178040083610d0,   4.137556867084717d0,  -8.108447067502766d0, &   
                          6.941780400836100d0,  -2.887556867084716d0,   0.560844706750277d0, &
                          0.000000000000000d0,   0.000000000000000d0,   0.000000000000000d0  /)
        fd2%DL(:,3) = (/  0.213357591590471d0,  -1.159078186228774d0,   3.511693723953474d0, &  
                         -4.723144865335573d0,   2.489690240716552d0,  -0.341475399639236d0, &   
                          0.008956894943086d0,   0.000000000000000d0,   0.000000000000000d0  /)
        fd2%DL(:,4) = (/ -0.179078329313190d0,   0.915651051584783d0,  -1.987601032542000d0, &   
                          3.387647581566586d0,  -3.719859506580339d0,   1.735582497566756d0, &  
                         -0.164529643265203d0,   0.012187380982608d0,   0.000000000000000d0  /)
        fd2%DL(:,5) = (/  0.040020014763742d0,  -0.187522354892963d0,   0.347126777927445d0, &  
                         -0.417791070219097d0,   1.560601736642238d0,  -2.684870208442730d0, &   
                          1.479418278121504d0,  -0.147941827812150d0,   0.010958653912011d0  /) 

        ! First derivative for the left boundary
        
        fd2%D1L = (/ -25/12d0, 4d0, -3d0, 4/3d0, -1/4d0 /)
         
        ! Norm for the left boundary

        fd2%HLi = (/ 0.315949074074074d0,1.390393518518518d0,0.627546296296296d0, &
        1.240509259259259d0,0.911689814814815d0 /)
    end select
        
        ! right boundary

        do i = 0,fd2%nbst-1
            fd2%DR(:,-i) = fd2%DL(fd2%nbnd-1:0:-1,i)
        end do

        ! First derivative for the right boundary

        fd2%D1R = -fd2%D1L(fd2%nD1:1:-1)

        ! Invert norm for the left boundary
        fd2%HLi = 1/fd2%HLi

        ! Norm for the right boundary
        fd2%HRi = fd2%HLi(fd2%nD1:1:-1)


  end subroutine

  ! Construct full, second derivative matrix of size n x n (used for linear system solve)
  subroutine init_fd2_full(HF,n)

      implicit none

      type(hf_type),intent(inout) :: HF
      integer,intent(in) :: n

      integer :: i
      real :: DT_BL, DT_BR


      allocate(HF%fd2%D2(n,n))
      HF%fd2%D2 = 0d0

      ! Put second derivative stencils in matrix
      
      ! Interior
      do i = HF%fd2%nbst+1,n-HF%fd2%nbst
         HF%fd2%D2(i,i+HF%fd2%mI:i+HF%fd2%pI) = HF%fd2%DI
      end do

      ! Left boundary
      do i = 0,HF%fd2%nbst-1
         HF%fd2%D2(1 + i,:) = HF%fd2%DL(:, i)
      end do
    
      ! Right boundary
      do i = 0,HF%fd2%nbst-1         
         HF%fd2%D2(n - i,n - HF%fd2%nbnd+1:n ) = HF%fd2%DR(:,-i)
      end do

      ! SAT penalty matrices
      ! D2*v + BLm*(v - g) + BLD*(Dv - g) + D^TBL*(v - g) + D^TBLD*(v - g) 

      ! Starting by implementing velocities in the fluid only 
      !gm = *( HF%v(1,i)  - vtp(i) )/HF%hy(i)
      HF%fd2%D2 = HF%mu/HF%rho0*HF%fd2%D2
      call diff_T_bnd_m(HF%fd2%D2(:,1)   ,HF%SATym(3),HF%fd2)
      call diff_T_bnd_p(HF%fd2%D2(:,HF%n),HF%SATyp(3),HF%fd2)

      ! Used for debugging
      !call print_matrix(HF%fd2%D2)

  end subroutine


subroutine print_matrix(A)
      implicit none

      real(kind=8),dimension(:,:),intent(in) :: A
      integer :: i

      do i=1,size(A,1) 
        print *,A(i,:)
        ! write (*,'(F10.2)') A(i,:)
        ! print *,""
      end do
end subroutine 
end module hydrofrac
