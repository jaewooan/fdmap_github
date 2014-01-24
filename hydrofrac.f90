module hydrofrac

  use fd, only : limits

  implicit none

  ! USE_HF = flag indicating hydraulic fracture
  ! BNDM,BNDP = process is resposible for minus/plus boundary 
  ! DIRECTION = direction of fracture
  ! COUPLED = full fluid-solid coupling (F=neglect crack opening in mass balance)
  ! LINEARIZED_WALLS = linearize wall opening term in mass balance
  ! SOURCE_TERM = add source term to mass balance equation
  ! OPERATOR_SPLITTING = use operator splitting to handle stiffness
  ! RHO0 = (initial) density
  ! K0 = bulk modulus
  ! C0 = sound speed
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
  ! SATM,SATP = SAT penalty weights (~wave speed / grid spacing)
  ! L = limits type (indices for FD operations)
  ! XSOURCE,YSOURCE = source position
  ! TSOURCE = source duration
  ! ASOURCE = source amplitude
  ! WSOURCE = source width (in space)

  type :: hf_type
     logical :: use_HF,coupled,linearized_walls,source_term,operator_splitting,bndm,bndp
     integer :: n
     character(1) :: direction
     real :: rho0,K0,c0,h,SATm,SATp,xsource,ysource,tsource,Asource,wsource
     real,dimension(:),allocatable :: wm,wp,wm0,wp0,u,p,Du,Dp,Dwm,Dwp
     real,dimension(:,:),allocatable :: v,Dv
     type(limits) :: L
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

    logical :: hydrofrac_file,coupled,linearized_walls,source_term,operator_splitting
    integer :: n
    real :: rho0,K0,w0,xsource,ysource,tsource,Asource,wsource
    character(256) :: filename
    integer,parameter :: nb=3 ! number of additional boundary (ghost) points needed for FD operators

    namelist /hydrofrac_list/ rho0,K0,w0,xsource,ysource,tsource,Asource,wsource, &
         n,coupled,linearized_walls,source_term,operator_splitting,hydrofrac_file,filename

    ! defaults
       
    rho0 = 0d0
    K0 = 1d40
    w0 = 0d0
    xsource = 0d0
    ysource = 0d0
    tsource = 1d0
    Asource = 1d0
    wsource = 1d0
    n = 1
    coupled = .true.
    linearized_walls = .false.
    hydrofrac_file = .false.
    source_term = .true.
    operator_splitting = .true.
    filename = ''

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
    HF%coupled = coupled
    HF%linearized_walls = linearized_walls
    HF%operator_splitting = operator_splitting
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
         HF%wm (HF%L%m :HF%L%p ),HF%wp (HF%L%m :HF%L%p ), &
         HF%wm0(HF%L%m :HF%L%p ),HF%wp0(HF%L%m :HF%L%p ), &
         HF%Dwm(HF%L%m :HF%L%p ),HF%Dwp(HF%L%m :HF%L%p ), &
         HF%u  (HF%L%mb:HF%L%pb),HF%p  (HF%L%mb:HF%L%pb), &
         HF%Du (HF%L%m :HF%L%p ),HF%Dp (HF%L%m :HF%L%p ), &
         HF%v (HF%n,HF%L%m:HF%L%p), &
         HF%Dv(HF%n,HF%L%m:HF%L%p))

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

    HF%u = 0d0
    HF%v = 0d0
    HF%p = 0d0
    
    !HF%p(HF%L%m:HF%L%p) = exp(-x**2) ! gaussian initial condition (hard-coded)

    ! override uniform initial conditions with those from file, if desired

    if (hydrofrac_file) then
       ! both sides read file (so process may read file twice)
       if (process_m) call read_hydrofrac(HF,filename,comm_m,array)
       if (process_p) call read_hydrofrac(HF,filename,comm_p,array)
    end if

    ! current (=initial) wall positions

    HF%wm = HF%wm0
    HF%wp = HF%wp0

  end subroutine init_hydrofrac


  subroutine read_hydrofrac(HF,filename,comm,array)

    use io, only : file_distributed,open_file_distributed, &
         read_file_distributed,close_file_distributed
    use mpi_routines, only : pw

    implicit none

    type(hf_type),intent(inout) :: HF
    character(*),intent(in) :: filename
    integer,intent(in) :: comm,array

    integer :: j
    type(file_distributed) :: fh

    call open_file_distributed(fh,filename,'read',comm,array,pw)

    call read_file_distributed(fh,HF%wm0(HF%L%m:HF%L%p))
    call read_file_distributed(fh,HF%wp0(HF%L%m:HF%L%p))
    call read_file_distributed(fh,HF%u  (HF%L%m:HF%L%p))
    call read_file_distributed(fh,HF%p  (HF%L%m:HF%L%p))
    do j = 1,HF%n
       call read_file_distributed(fh,HF%v(j,HF%L%m:HF%L%p))
    end do

    call close_file_distributed(fh)

  end subroutine read_hydrofrac


  subroutine destroy_hydrofrac(HF)

    implicit none

    type(hf_type),intent(inout) :: HF

    if (allocated(HF%wm0)) deallocate(HF%wm0)
    if (allocated(HF%wp0)) deallocate(HF%wp0)
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
       HF%u (i) = HF%u (i)+dt*HF%Du (i)
       HF%p (i) = HF%p (i)+dt*HF%Dp (i)
       HF%v(:,i) = HF%v(:,i)+dt*HF%Dv(:,i)
    end do

  end subroutine update_fields_hydrofrac


  subroutine set_rates_hydrofrac(HF,C,m,p,phip,vnm,vnp,x,y,t)

    use mpi_routines2d, only : cartesian
    use fd, only : diff

    implicit none

    type(hf_type),intent(inout) :: HF
    type(cartesian),intent(in) :: C
    integer,intent(in) :: m,p
    real,intent(in) :: phip(m:p),vnm(m:p),vnp(m:p),x(m:p),y(m:p),t

    integer :: i
    real :: uhat,phat
    real,dimension(:),allocatable :: dudx,dpdx

    if (.not.HF%use_HF) return

    ! populate ghost points of neighboring processors

    call share_hydrofrac(HF,C)

    ! allocate auxiliary arrays

    allocate(dudx(HF%L%m:HF%L%p),dpdx(HF%L%m:HF%L%p))

    ! SBP differentiation of velocity and pressure
    ! (this could certainly be improved for compatibility with external mesh)

    call diff(HF%L,HF%u,dudx)
    call diff(HF%L,HF%p,dpdx)
    dudx = dudx/HF%h
    dpdx = dpdx/HF%h

    ! set rates from mass and momentum balance equations,
    ! starting with linearized acoustics (rigid walls)

    do i = HF%L%m,HF%L%p
       HF%Du(i) = HF%Du(i)-dpdx(i)/HF%rho0
       !HF%Dv(:,i) = HF%Dv(:,i)-dpdx(i)/HF%rho0
       HF%Dp(i) = HF%Dp(i)-dudx(i)*HF%K0
    end do
    
    do i = HF%L%m,HF%L%p
       HF%Dwp(i) = HF%Dwp(i) + vnp(i)
       HF%Dwm(i) = HF%Dwm(i) + vnm(i)
    end do

    ! and then adding crack opening/closing term in mass balance

    if (HF%coupled) then
       if (HF%linearized_walls) then
          if (HF%operator_splitting) then
             do i = HF%L%m,HF%L%p
                HF%Dp(i) = HF%Dp(i)-HF%K0*phip(i)/(HF%wp0(i)-HF%wm0(i)) ! *INCORRECT*
             end do
          else
             do i = HF%L%m,HF%L%p
                HF%Dp(i) = HF%Dp(i)-HF%K0*(vnp(i)-vnm(i))/(HF%wp0(i)-HF%wm0(i))
             end do
          end if
       else
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
    end if

    ! and source terms (explosion source only appears in mass balance)

    if (HF%source_term) then
       do i = HF%L%m,HF%L%p
          HF%Dp(i) = HF%Dp(i)+HF%Asource* &
               exp(-0.5d0*((x(i)-HF%xsource)**2+(y(i)-HF%ysource)**2)/HF%wsource**2)* &
               (t/HF%tsource)*exp(-t/HF%tsource)
       end do
    end if

    ! and viscous terms (explicit time-stepping)

    !if (.not.HF%operator_splitting) then
    !   do i = HF%L%m,HF%L%p
    !      call second_derivative(v(:,i),v_yy)
    !      HF%Dv(:,i) = HF%Dv(:,i)+HF%mu*v_yy+SAT_terms
    !   end do
    !end if

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

    ! deallocate temporary arrays

    deallocate(dudx,dpdx)

  end subroutine set_rates_hydrofrac


  subroutine fluid_stresses(HF,i,p,taum,taup)

    implicit none

    type(HF_type),intent(in) :: HF
    integer,intent(in) :: i
    real,intent(out) :: p,taum,taup

    ! fluid pressure

    p = HF%p(i)

    ! shear stress on bottom (taum) and top (taup) walls,
    ! evaluate these as mu*dv/dy for viscous fluid with fluid velocity v

    taum = 0d0
    taup = 0d0

  end subroutine fluid_stresses

  subroutine update_fields_hydrofrac_implicit(HF,m,p,Zm,Zp,dt)
  
    implicit none

    type(HF_type),intent(inout) :: HF
    integer,intent(in) :: m,p
    real,dimension(m:p),intent(in) :: Zm,Zp ! P-wave impedances on minus,plus sides
    real,intent(in) :: dt

    integer :: i
    real :: Z

    if (.not.HF%operator_splitting) return ! return if fully explicit

    ! implicit update for pressure
    
    if (HF%linearized_walls) then
       do i = m,p
          Z = 1d0/(1d0/Zm(i)+1d0/Zp(i)) ! combine impedances of two sides
          HF%p(i) = HF%p(i)*(1d0-dt*HF%K0/(Z*(HF%wp0(i)-HF%wm0(i)))) ! forward Euler
          !HF%p(i) = HF%p(i)/(1d0+dt*HF%K0/(Z*(HF%wp0(i)-HF%wm0(i)))) ! backward Euler
       end do
    else
       do i = m,p
          Z = 1d0/(1d0/Zm(i)+1d0/Zp(i)) ! combine impedances of two sides
          ! at what time should we use wp-wm?
          HF%p(i) = HF%p(i)/(1d0+dt*HF%K0/(Z*(HF%wp (i)-HF%wm (i)))) ! backward Euler
       end do
    end if

    ! implicit update for viscous diffusion terms

    !do i = m,p
       !HF%v(:,i) = HF%v(:,i)+dt*mu*v_yy ! forward Euler
       !call solve_linear_system(HF%v(:,i),...) ! backward Euler
    !end do

    ! use LAPACK for linear system solve
    ! call dgesv(...)

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


end module hydrofrac
