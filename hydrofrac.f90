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
  ! D2 = full second derivative matrix including boundary terms with scaling h^2 (used by implicit solve)
  ! bD2 = Same as D2 but stored in LAPACK's banded matrix format
  ! H1Li/H1Ri = inverted diagonal norm coefficients near boundaries (needed by penalty terms of the
  ! form D^T)
  ! HL/HR = diagonal norm entries near boundary (normalized so H=1 in interior)
  ! FDmethod = finite difference method for second derivative
  type :: fd2_type
    integer :: nI,mI,pI,nbnd,nbst,nD1
    real :: H00i
    real,dimension(:),allocatable :: DI,DmI,DpI,HL,HR,D1L,D1R,HLi,HRi
    real,dimension(:,:),allocatable :: DL,DR,D2,bD2
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
  ! DYDETA = Gradient for variable gridspacing in the y-direction
  ! L = limits type (indices for FD operations)
  ! nsource = number of sources
  ! XSOURCE,YSOURCE = source position
  ! TSOURCE = source duration
  ! T0SOURCE = source start time
  ! TW0SOURCE = source start time interval
  ! ASOURCE = source amplitude
  ! WSOURCE = source width (in space)
  ! Lsource = length which multiple sources occupy
  ! bcL/bcR = boundary conditions on left and right boundary (sin_p for
  ! oscillatory pressure input) (zero_u for no velocity (default))
  ! bcLA = Amplitude
  ! bcLomega = Angular frequency
  ! bcLphi = phase
  ! bcLtend = shut off time. When t > tend the oscillatory pressure input is
  ! disabled and the zero_u b.c. is used instead. 
  ! fd2_type  = finite difference second derivative operator
  ! USE_MMS = enable MMS if true

  type :: hf_type
     logical :: use_HF,coupled,linearized_walls,source_term,operator_splitting,bndm,bndp,inviscid,slope
     logical :: banded_storage,variable_grid_spacing,use_mms
     integer :: n,nsource
     character(1) :: direction
     real :: &
     rho0,K0,mu,c0,&
     h,SATm,SATp,xsource,ysource,tsource,Asource,wsource,Lsource,tw0source,t0source, &
     bcLA,bcLomega,bcLphi,bcRA,bcRomega,bcRphi,bcLtend,bcRtend
     character(256) :: bcL,bcR,distsource
     real,dimension(:),allocatable :: &
     wm,wp,wm0,wp0,dwm0dx,dwp0dx,u,p,Du,Dp,Dwm,Dwp,Hw,hy,SATyp,SATym,dydetai,y,&
     rxsource,rtsource
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
    character(256) :: str,bcL,bcR

    logical :: coupled,linearized_walls,source_term,operator_splitting, &
               inviscid,banded_storage,variable_grid_spacing,use_mms
    integer :: n,nsource
    real :: rho0,K0,mu,w0,&
            xsource,ysource,tsource,Asource,wsource,Lsource,tw0source,t0source, &
            bcLA,bcLomega,bcLphi,bcRA,bcRomega,bcRphi,bcLtend,bcRtend
    character(256) :: geom_file,initial_conds_file,var_file,distsource
    character(4) :: FDmethod
    integer,parameter :: nb=3 ! number of additional boundary (ghost) points needed for FD operators

    namelist /hydrofrac_list/ rho0,K0,mu,w0,&
         nsource,xsource,ysource,tsource,Asource,wsource,Lsource,tw0source,t0source,distsource, &
         n,coupled,linearized_walls,source_term,operator_splitting,geom_file,initial_conds_file, &
         inviscid,FDmethod,bcL,bcR,bcLA,bcLomega,bcLphi,bcRA,bcRomega,bcRphi,bcLtend,bcRtend,&
         banded_storage,var_file, &
         use_mms

    ! defaults
       
    rho0 = 0d0
    K0 = 1d40
    w0 = 0d0
    mu = 0d0
    distsource = 'equidistant'
    nsource = 1
    xsource = 0d0
    ysource = 0d0
    tsource = 1d0
    t0source = 0d0
    tw0source = 0d0
    Asource = 1d0
    Lsource = 0d0
    wsource = 1d0
    n = 1
    coupled = .true.
    linearized_walls = .false.
    source_term = .true.
    operator_splitting = .true.
    inviscid = .true.
    variable_grid_spacing = .false.
    geom_file = ''
    initial_conds_file = ''
    var_file = ''
    FDmethod = 'SBP6'
    bcL = 'zero_u'
    bcLtend = -1d0
    bcLA = 1d0
    bcLomega = 1d0
    bcLphi = 0d0
    bcR = 'zero_u'
    bcRtend = -1d0
    bcRA = 1d0
    bcRomega = 1d0
    bcRphi = 0d0
    banded_storage = .true.
    use_mms = .false.

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
    HF%nsource = nsource
    HF%ysource = ysource
    HF%Lsource = Lsource
    HF%tsource = tsource
    HF%t0source = t0source
    HF%tw0source = tw0source
    HF%distsource = distsource
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

    ! Boundary conditions 
    HF%bcL = bcL
    HF%bcLtend = bcLtend
    HF%bcLA = bcLA
    HF%bcLomega = bcLomega
    HF%bcR = bcR
    HF%bcRA = bcRA
    HF%bcRomega = bcRomega
    HF%bcRtend = bcRtend

    HF%banded_storage = banded_storage

    ! MMS
    HF%use_mms = use_mms

    ! Geometry
    HF%slope = .false.
    HF%variable_grid_spacing = .false.

    ! Geometry file 
    if (geom_file /= '' ) then
         HF%slope = .true.
    end if

    ! Variable grid spacing file
    if (var_file /= '' ) then
      HF%variable_grid_spacing = .true.
    end if

    ! Initialize source terms that use random distributions
    call init_source_terms(HF)

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
          call write_matlab(echo,'nsource',HF%nsource,HFstr)
          call write_matlab(echo,'xsource',HF%xsource,HFstr)
          call write_matlab(echo,'ysource',HF%ysource,HFstr)
          call write_matlab(echo,'tsource',HF%tsource,HFstr)
          call write_matlab(echo,'t0source',HF%t0source,HFstr)
          call write_matlab(echo,'tw0source',HF%tw0source,HFstr)
          call write_matlab(echo,'Lsource',HF%Lsource,HFstr)
          call write_matlab(echo,'Asource',HF%Asource,HFstr)
          call write_matlab(echo,'wsource',HF%wsource,HFstr)
       end if

       ! Boundary conditions
       call write_matlab(echo,'bcL',HF%bcL,HFstr)
       call write_matlab(echo,'bcLA',HF%bcLA,HFstr)
       call write_matlab(echo,'bcLomega',HF%bcLomega,HFstr)
       call write_matlab(echo,'bcR',HF%bcR,HFstr)
       call write_matlab(echo,'bcRA',HF%bcRA,HFstr)
       call write_matlab(echo,'bcRomega',HF%bcRomega,HFstr)

       ! Geometry
       call write_matlab(echo,'custom_hf_profile',HF%slope,HFstr)
       call write_matlab(echo,'variable_grid_spacing',HF%variable_grid_spacing,HFstr)

       ! MMS
       call write_matlab(echo,'use_mms',HF%use_mms,HFstr)

       ! Geometry file 
       if (geom_file /= '' ) then
         call write_matlab(echo,'geom_file',trim(geom_file),HFstr)
            HF%slope = .true.
       end if

       ! Variable grid spacing file
       if (var_file /= '' ) then
         call write_matlab(echo,'var_file',trim(var_file),HFstr)
            HF%slope = .true.
         HF%variable_grid_spacing = .true.
       end if

       ! Initial conditions file
       if (initial_conds_file /= '' ) then
         call write_matlab(echo,'initial_conds_file',initial_conds_file,HFstr)
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
         HF%SATyp(4),HF%SATym(4), &
         HF%dydetai(HF%n), &
         HF%y(HF%n) )

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
    HF%dydetai = 1d40
    HF%y = 1d40

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

        HF%slope = .true.
        ! both sides read file (so process may read file twice)
        if (process_m) call read_hydrofrac_geom(HF,geom_file,comm_m,array)
        if (process_p) call read_hydrofrac_geom(HF,geom_file,comm_p,array)


    end if

    ! Load grid points with variable grid spacing
    if (var_file /= '' ) then

        HF%variable_grid_spacing = .true.
        ! both sides read file (so process may read file twice)
        if (process_m) call read_hydrofrac_var_grid(HF,var_file,comm_m)
        if (process_p) call read_hydrofrac_var_grid(HF,var_file,comm_p)

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

    HF%SATym(1) = 0d0
    HF%SATym(2) = 0d0
    HF%SATym(3) = HF%mu/HF%rho0
    HF%SATym(4) = 0d0
    
    HF%SATyp(1) = 0d0
    HF%SATyp(2) = 0d0
    HF%SATyp(3) = - HF%mu/HF%rho0
    HF%SATyp(4) = 0d0
    
    
    ! Initialize second derivative SBP operator
    if( .not. inviscid) then
        HF%fd2%FDmethod = FDmethod
        call init_fd2(FDmethod,HF%fd2)
        if(operator_splitting) then
            ! Either use variable grid spacing or 
            ! uniform grid spacing
            ! the D2 operator can be either full or
            ! packed using banded matrix format
            if(HF%variable_grid_spacing) then
              call init_fd2var_full(HF,HF%n)
            else
              call init_fd2_full(HF,HF%n)
            end if
        end if
    end if 
    
    ! Initialize solution using MMS
    if(HF%use_mms) call mms_init(HF,m,p,x,y)
    

  end subroutine init_hydrofrac


   subroutine read_hydrofrac_geom(HF,geom_file,comm,array)

     use io, only : file_distributed,open_file_distributed, &
          read_file_distributed,close_file_distributed
     use mpi_routines, only : pw, myid

     implicit none

     type(hf_type),intent(inout) :: HF
     character(*),intent(in) :: geom_file
     integer,intent(in) :: comm,array

     type(file_distributed) :: fh

     call open_file_distributed(fh,geom_file,'read',comm,array,pw)

     call read_file_distributed(fh,HF%wm0(HF%L%m:HF%L%p))
     call read_file_distributed(fh,HF%wp0(HF%L%m:HF%L%p))
     call read_file_distributed(fh,HF%dwm0dx(HF%L%m:HF%L%p))
     call read_file_distributed(fh,HF%dwp0dx(HF%L%m:HF%L%p))

     call close_file_distributed(fh)

   end subroutine read_hydrofrac_geom
   
   subroutine read_hydrofrac_var_grid(HF,var_file,comm)


     use io, only : file_distributed,open_file_distributed, &
          read_file_distributed,close_file_distributed
     use mpi_routines, only : pw, myid
     use mpi, only: MPI_MODE_RDONLY,MPI_INFO_NULL

     implicit none

     type(hf_type),intent(inout) :: HF
     character(*),intent(in) :: var_file
     integer,intent(in) :: comm

     type(file_distributed) :: fh
     integer :: ierr


     call MPI_File_open(comm,var_file,MPI_MODE_RDONLY, &
              MPI_INFO_NULL,fh%fh,ierr)
     call read_file_distributed(fh,HF%dydetai(1:HF%n))
     call read_file_distributed(fh,HF%y(1:HF%n))
     call close_file_distributed(fh)



   end subroutine read_hydrofrac_var_grid
  
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
    if (allocated(HF%rxsource  )) deallocate(HF%rxsource  )
    if (allocated(HF%rtsource  )) deallocate(HF%rtsource  )
    if (allocated(HF%Dwm)) deallocate(HF%Dwm)
    if (allocated(HF%Dwp)) deallocate(HF%Dwp)
    if (allocated(HF%Du )) deallocate(HF%Du )
    if (allocated(HF%Dv )) deallocate(HF%Dv )
    if (allocated(HF%Dp )) deallocate(HF%Dp )
    if (allocated(HF%Hw )) deallocate(HF%Hw )
    if (allocated(HF%hy )) deallocate(HF%hy )
    if (allocated(HF%SATyp )) deallocate(HF%SATyp )
    if (allocated(HF%SATym )) deallocate(HF%SATym )
    if (allocated(HF%dydetai )) deallocate(HF%dydetai)
    if (allocated(HF%y )) deallocate(HF%y)
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
    if (allocated(HF%fd2%bD2 )) deallocate(HF%fd2%bD2)

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


  subroutine set_rates_hydrofrac(HF,C,m,p,phip,vnm,vnp,vtm,vtp,x,y,t)

    use mpi_routines2d, only : cartesian
    use fd, only : diff
    use mms, only : mms_hydrofrac

    implicit none

    type(hf_type),intent(inout) :: HF
    type(cartesian),intent(in) :: C
    integer,intent(in) :: m,p
    real,intent(in) :: phip(m:p),vnm(m:p),vnp(m:p),vtm(m:p),vtp(m:p), &
                       x(m:p),y(m:p),t

    integer :: i
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
    
    !if(HF%use_mms) then
    !    do i = HF%L%mb,HF%L%pb
    !        HF%p(i) = mms_hydrofrac(x(i),0d0,t,0,'p')
    !    end do
    !end if

    b = HF%wp0 - HF%wm0
    call diff(HF%L,HF%u*b,dudx)
    call diff(HF%L,HF%p,dpdx)
    dudx = dudx/HF%h
    dpdx = dpdx/HF%h
    
    !if(HF%use_mms) then
    !    do i = HF%L%m,HF%L%p
    !        dpdx(i) = mms_hydrofrac(x(i),0d0,t,0,'dpdx')
    !    end do
    !end if

    ! set rates from mass and momentum balance equations,
    ! starting with linearized acoustics (rigid walls)

    if(HF%use_mms) then
        do i = HF%L%m,HF%L%p
           HF%Dv(:,i) = HF%Dv(:,i)- dpdx(i)/HF%rho0
           HF%Dp(i)   = HF%Dp(i)-  dudx(i)*HF%K0/b(i)
        end do
    else
        do i = HF%L%m,HF%L%p
           HF%Dv(:,i) = HF%Dv(:,i)-dpdx(i)/HF%rho0
           HF%Dp(i) = HF%Dp(i)-dudx(i)*HF%K0/b(i)
        end do
    end if
    
    ! Apply boundary conditions at the left and right boundary
    call set_bc(HF,t)

    ! and source terms (explosion source only appears in mass balance)
    call set_source_terms(HF,i,x,y,t)




    call set_rates_viscous(HF,m,p,vtm,vtp)
    
    ! Use mms
    if(HF%use_mms) call mms_source_terms(HF,m,p,x,y,t)

    if (.not. HF%coupled) return


    do i = HF%L%m,HF%L%p
       HF%Dwp(i) = HF%Dwp(i) + vnp(i)
       HF%Dwm(i) = HF%Dwm(i) + vnm(i)
    end do

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
    
  subroutine init_source_terms(HF)

      implicit none

      type(hf_type),intent(inout) :: HF
      
      if(.not. HF%source_term) return

      select case(HF%distsource)
      case('uniform')
      allocate(HF%rxsource(HF%nsource),HF%rtsource(HF%nsource))
      !todo: each process gets a different seed which is bad
      ! for now let the random numbers be the same each time
      !call init_random_seed()
      call random_number(HF%rxsource)
      call random_number(HF%rtsource)
      end select

  end subroutine init_source_terms
  
  subroutine init_random_seed()
            use iso_fortran_env, only: int64
            implicit none
            integer, allocatable :: seed(:)
            integer :: i, n, un, istat, dt(8), pid
            integer(int64) :: t
          
            call random_seed(size = n)
            allocate(seed(n))
            ! First try if the OS provides a random number generator
            open(newunit=un, file="/dev/urandom", access="stream", &
                 form="unformatted", action="read", status="old", iostat=istat)
            if (istat == 0) then
               read(un) seed
               close(un)
            else
               ! Fallback to XOR:ing the current time and pid. The PID is
               ! useful in case one launches multiple instances of the same
               ! program in parallel.
               call system_clock(t)
               if (t == 0) then
                  call date_and_time(values=dt)
                  t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
                       + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
                       + dt(3) * 24_int64 * 60 * 60 * 1000 &
                       + dt(5) * 60 * 60 * 1000 &
                       + dt(6) * 60 * 1000 + dt(7) * 1000 &
                       + dt(8)
               end if
               pid = getpid()
               t = ieor(t, int(pid, kind(t)))
               do i = 1, n
                  seed(i) = lcg(t)
               end do
            end if
            call random_seed(put=seed)
          contains
            ! This simple PRNG might not be good enough for real work, but is
            ! sufficient for seeding a better PRNG.
            function lcg(s)
              integer :: lcg
              integer(int64) :: s
              if (s == 0) then
                 s = 104729
              else
                 s = mod(s, 4294967296_int64)
              end if
              s = mod(s * 279470273_int64, 4294967291_int64)
              lcg = int(mod(s, int(huge(0), int64)), kind(0))
            end function lcg
  end subroutine init_random_seed

  subroutine set_source_terms(HF,i,x,y,t)

      use io, only: error
      implicit none

      type(hf_type),intent(inout) :: HF
      integer,intent(in) :: i
      real,dimension(HF%L%m:HF%L%p),intent(in) :: x,y
      real,intent(in) :: t

      if(.not. HF%source_term) return

      select case(HF%distsource)
      case('equidistant')
          call set_equidistant_source_terms(HF,i,x,y,t)
      case('uniform')
          call set_uniform_source_terms(HF,i,x,y,t)
      case default
          call error('Invalid type of source term distribution', &
                    'hydrofrac:: set_source_terms')
      end select

  end subroutine set_source_terms

  subroutine set_equidistant_source_terms(HF,i,x,y,t)

      implicit none
      
      type(hf_type),intent(inout) :: HF
      integer,intent(in) :: i
      real,dimension(HF%L%m:HF%L%p),intent(in) :: x,y
      real,intent(in) :: t

      real :: x0,y0,dx,dy,x_,y_
      integer :: j

      ! For now, we always pretend that crack is horizontal
      ! todo: handle the general case when (x,y) is used

      ! Special case, only one source term
      if(HF%nsource == 1) then
      HF%Dp = HF%Dp+HF%Asource*&
              exp(-0.5d0*((x-HF%xsource)**2+(y-HF%ysource)**2)/HF%wsource**2) &
              *exp(-0.5d0*(t-HF%t0source)**2/HF%tsource**2)
          return
      end if

      ! More than one source term
      
      x0 = - HF%Lsource*0.5 + HF%xsource
      y0 = - HF%Lsource*0.5 + HF%ysource
      !xn = + HF%Lsource*0.5 + HF%xsource
      dx = HF%Lsource/(HF%nsource - 1)
      dy = 0d0

      do j=1,HF%nsource
        x_ = x0 + dx*(j-1)
        y_ = 0d0*y0 + dy*(j-1)
        HF%Dp = HF%Dp+HF%Asource*&
              exp(-0.5d0*((x-x_)**2+(y-y_)**2)/HF%wsource**2)* &
              exp(-0.5d0*(t-HF%t0source)**2/HF%tsource**2)
      end do

  end subroutine set_equidistant_source_terms
  
  subroutine set_uniform_source_terms(HF,i,x,y,t)

      implicit none
      
      type(hf_type),intent(inout) :: HF
      integer,intent(in) :: i
      real,dimension(HF%L%m:HF%L%p),intent(in) :: x,y
      real,intent(in) :: t

      real :: x0,y0,dx,dy,x_,y_,t0_
      integer :: j
      
      ! Special case, only one source term
      if(HF%nsource == 1) then
      HF%Dp = HF%Dp+HF%Asource*&
              exp(-0.5d0*((x-HF%xsource)**2+(y-HF%ysource)**2)/HF%wsource**2) &
              *exp(-0.5d0*(t-HF%t0source)**2/HF%tsource**2)
          return
      end if

      x0 = - HF%Lsource*0.5 + HF%xsource
      y0 = - HF%Lsource*0.5 + HF%ysource

      dx = HF%Lsource/(HF%nsource - 1)
      dy = 0d0
      

      do j=1,HF%nsource
        x_ = x0 + HF%Lsource*HF%rxsource(j)
        t0_ = HF%tw0source*HF%rtsource(j)
        y_ = 0d0*y0 + dy*(j-1)
        HF%Dp = HF%Dp+HF%Asource*&
              exp(-0.5d0*((x-x_)**2+(y-y_)**2)/HF%wsource**2)* &
              exp(-0.5d0*(t-t0_)**2/HF%tsource**2)
      end do



  end subroutine set_uniform_source_terms




  subroutine set_bc(HF,t)

      use io, only : error

      type(hf_type),intent(inout) :: HF
      real,intent(in) :: t

      integer :: i
    
      ! add penalty terms to enforce BC with SAT method

      ! Left boundary
      select case(HF%bcL)
        case('zero_u')
            call bcL_zero_u(HF)
        case('sin_p')
            call bcL_sin_p(HF,t)
        case('sin_u')
            call bcL_sin_u(HF,t)
        case('mms-hydrofrac-w')
            call bcL_mms(HF,t)
        case default
            call error("Invalid boundary condition at left boundary", &
                "hydrofrac::set_bc")
      end select
      
      ! Right boundary
      select case(HF%bcR)
        case('zero_u')
            call bcR_zero_u(HF)
        case('sin_p')
            call bcR_sin_p(HF,t)
        case('sin_u')
            call bcR_sin_u(HF,t)
        case('mms-hydrofrac-w')
            call bcR_mms(HF,t)
        case default
            call error("Invalid boundary condition at right boundary", &
                "hydrofrac::set_bc")
      end select


  end subroutine

  subroutine bcL_sin_u(HF,t)
    
    type(hf_type),intent(inout) :: HF
    real,intent(in) :: t

    integer :: i
    real :: uhat,phat

    if (HF%bndm) then ! check if process handles minus boundary
       i = HF%L%mg
       uhat = hf%bcLA*sin(hf%bcLomega*(t - hf%bcLphi)) ! oscillatory velocity input
       phat = (hf%rho0*hf%c0*uhat + (hf%p(i) - hf%rho0*hf%c0*hf%u(i)) )

       HF%Dv(:,i) = HF%Dv(:,i)-HF%SATm*(HF%u(i)-uhat)
       HF%Dp(i) = HF%Dp(i)-HF%SATm*(HF%p(i)-phat)
    end if

  end subroutine
  
  subroutine bcR_sin_u(HF,t)
    
    type(hf_type),intent(inout) :: HF
    real,intent(in) :: t

    integer :: i
    real :: uhat,phat
    
    if (HF%bndp) then ! check if process handles minus boundary
       i = HF%L%pg
       uhat = HF%bcRA*sin(HF%bcRomega*(t - HF%bcRphi)) ! Oscillatory pressure input
       phat = ( - hf%rho0*hf%c0*uhat + (hf%p(i) + hf%rho0*hf%c0*hf%u(i)) )
       HF%Du(i) = HF%Du(i)-HF%SATp*(HF%u(i)-uhat)
       HF%Dp(i) = HF%Dp(i)-HF%SATp*(HF%p(i)-phat)
    end if


  end subroutine
  
  subroutine bcL_sin_p(HF,t)
    
    type(hf_type),intent(inout) :: HF
    real,intent(in) :: t

    integer :: i
    real :: uhat,phat


    if (HF%bndm) then ! check if process handles minus boundary
       i = HF%L%mg
       phat = hf%bcLA*sin(hf%bcLomega*(t - hf%bcLphi))! oscillatory pressure input
       if(HF%bcLtend < t .and. HF%bcLtend > 0) then
           !call bcL_zero_u(HF)
           !return
           phat = 0d0
       end if
       uhat = (phat - (hf%p(i) - hf%rho0*hf%c0*hf%u(i)) )/(hf%rho0*hf%c0)

       HF%Dv(:,i) = HF%Dv(:,i)-HF%SATm*(HF%u(i)-uhat)
       !HF%Dv(:,i) = HF%Dv(:,i)-HF%SATm*(HF%v(:,i)-uhat)
       HF%Dp(i) = HF%Dp(i)-HF%SATm*(HF%p(i)-phat)
    end if

  end subroutine
  
  subroutine bcR_sin_p(HF,t)
    
    type(hf_type),intent(inout) :: HF
    real,intent(in) :: t

    integer :: i
    real :: uhat,phat
    
    if(HF%bcRtend < t .and. HF%bcRtend > 0) then
        call bcR_zero_u(HF)
        return
    end if
    
    if (HF%bndp) then ! check if process handles minus boundary
       i = HF%L%pg
       phat = HF%bcRA*(sin(HF%bcRomega*(t - HF%bcRphi))) ! Oscillatory pressure input
       uhat = (HF%p(i)+HF%rho0*HF%c0*HF%u(i) - phat)/(HF%rho0*HF%c0)
       HF%Du(i) = HF%Du(i)-HF%SATp*(HF%u(i)-uhat)
       HF%Dp(i) = HF%Dp(i)-HF%SATp*(HF%p(i)-phat)
    end if


  end subroutine

  subroutine bcL_zero_u(HF)
      
    type(hf_type),intent(inout) :: HF
    integer :: i
    real :: uhat,phat
    if (HF%bndm) then ! check if process handles minus boundary
       i = HF%L%mg
       uhat = 0d0 ! zero fluid velocity at crack tip
       phat = HF%p(i) - HF%rho0*HF%c0*HF%u(i)

       HF%Dv(:,i) = HF%Dv(:,i)- HF%SATm*(HF%u(i)-uhat)
       HF%Dp(i)   = HF%Dp(i)  - HF%SATm*(HF%p(i)-phat)
    end if

  end subroutine

  subroutine bcR_zero_u(HF)
    
    type(hf_type),intent(inout) :: HF
    integer :: i
    real :: uhat,phat

    if (HF%bndp) then ! check if process handles plus  boundary
       i = HF%L%pg
       uhat = 0d0 ! zero fluid velocity at crack tip
       phat = HF%p(i) + HF%rho0*HF%c0*HF%u(i) ! set by preserving characteristic variable into fluid

       HF%Dv(:,i) = HF%Dv(:,i)- HF%SATp*(HF%u(i)-uhat)
       HF%Dp(i)   = HF%Dp(i)  - HF%SATp*(HF%p(i)-phat)
    end if

  end subroutine
  
  subroutine bcL_mms(HF,t)
      use mms, only : mms_hydrofrac
      
    type(hf_type),intent(inout) :: HF
    real,intent(in) :: t
    integer :: i
    real :: uhat,phat
    if (HF%bndm) then ! check if process handles minus boundary
       i = HF%L%mg
       phat = mms_hydrofrac(0d0,0d0,t,0,'bcL_phat')
       uhat = (phat - (hf%p(i) - hf%rho0*hf%c0*hf%u(i)) )/(hf%rho0*hf%c0)

       HF%Dv(:,i) = HF%Dv(:,i)- HF%SATm*(HF%u(i)-uhat)
       HF%Dp(i)   = HF%Dp(i)  - HF%SATm*(HF%p(i)-phat)
    end if

  end subroutine

  subroutine bcR_mms(HF,t)
    use mms, only : mms_hydrofrac
    
    type(hf_type),intent(inout) :: HF
    real,intent(in) :: t
    integer :: i
    real :: uhat,phat

    if (HF%bndp) then ! check if process handles plus  boundary
       i = HF%L%pg
       phat = mms_hydrofrac(0d0,0d0,t,0,'bcR_phat')
       uhat = (HF%p(i)+HF%rho0*HF%c0*HF%u(i) - phat)/(HF%rho0*HF%c0)
       !phat = mms_hydrofrac(0d0,0d0,t,0,'bcR_phat')
       !uhat = (HF%p(i)+HF%rho0*HF%c0*HF%u(i) - phat)/(HF%rho0*HF%c0)

       HF%Dv(:,i) = HF%Dv(:,i)- HF%SATp*(HF%u(i)-uhat)
       HF%Dp(i)   = HF%Dp(i)  - HF%SATp*(HF%p(i)-phat)
    end if

  end subroutine

  subroutine mms_init(HF,m,p,x,y)
      use mms, only: mms_hydrofrac

      implicit none

      type(HF_type),intent(inout) :: HF
      integer,intent(in) :: m,p
      real,dimension(m:p),intent(in) :: x,y

      integer :: i,j
      real :: t,y_



      if(.not. HF%use_mms) return


      t = 0d0

      do i=m,p
            HF%p(i) = mms_hydrofrac(x(i),0d0,t,0,'p')
      end do
      ! Init MMS solution on equidistant grid
      if(.not. HF%variable_grid_spacing) then
        do i=m,p
          do j=1,HF%n
              y_ = HF%wm0(i) + HF%hy(i)*(j-1)
              HF%v(j,i) = mms_hydrofrac(x(i),y_,t,0,'v')
          end do
        end do
      end if

      ! Init MMS solution on curvilinear grid
      if(HF%variable_grid_spacing) then
        do i=m,p
          do j=1,HF%n
              y_ = HF%wm0(i) + (HF%wp0(i) - HF%wm0(i) )*HF%y(j)
              HF%v(j,i) = mms_hydrofrac(x(i),y_,t,0,'v')
          end do
        end do
      end if


  end subroutine

  subroutine mms_stresses(x,y,t,taum,taup)
      use mms, only: mms_hydrofrac

      implicit none

      real,intent(in) :: x,y,t
      real,intent(inout) :: taum,taup

      taum = taum + mms_hydrofrac(x,y,t,0,'I_taum')
      taup = taup + mms_hydrofrac(x,y,t,0,'I_taup')

  end subroutine
  
  subroutine mms_velocities(x,y,t,vtm,vtp)
      use mms, only: mms_hydrofrac

      implicit none

      real,intent(in) :: x,y,t
      real,intent(inout) :: vtm,vtp

      vtm = vtm + mms_hydrofrac(x,y,t,0,'I_vtm')
      vtp = vtp + mms_hydrofrac(x,y,t,0,'I_vtp')

  end subroutine

  subroutine mms_source_terms(HF,m,p,x,y,t)
      use mms, only: mms_hydrofrac

      implicit none

      type(HF_type),intent(inout) :: HF
      integer,intent(in) :: m,p
      real,intent(in) :: x(m:p),y(m:p),t

      integer :: i,j
      real :: y_
      
      do i=m,p
        do j=1,HF%n
            y_ = HF%wm0(i) + HF%hy(i)*(j-1)
            HF%Dv(j,i) = HF%Dv(j,i) + mms_hydrofrac(x(i),y_,t,0,'s_v')
        end do
            HF%Dp(i) = HF%Dp(i) + mms_hydrofrac(x(i),0d0,t,0,'s_p')
      end do

  end subroutine

  subroutine fluid_stresses(HF,i,x,y,t,p,taum,taup)
    use mms, only: mms_hydrofrac
    implicit none

    type(HF_type),intent(in) :: HF
    integer,intent(in) :: i
    real,intent(in) :: x,y,t
    real,intent(out) :: p,taum,taup

    ! fluid pressure

    p = HF%p(i)


    ! shear stress on bottom (taum) and top (taup) walls,
    ! evaluate these as mu*dv/dy for viscous fluid with fluid velocity v


    ! Inviscid case
    taum = 0d0
    taup = 0d0
    
    ! Add geometric correction
    if(HF%slope) then
      taum = taum + HF%p(i)*HF%dwm0dx(i)
      taup = taup + HF%p(i)*HF%dwp0dx(i)
    end if

    if(HF%inviscid) return

    ! Viscous case
    if(HF%variable_grid_spacing) then
        taum = taum + HF%mu*HF%dydetai(1)*diff_bnd_m(HF%v(:,i),HF%fd2)/HF%hy(i)
        taup = taup + HF%mu*HF%dydetai(HF%n)*diff_bnd_p(HF%v(:,i),HF%fd2)/HF%hy(i)
    end if

    if(.not. HF%variable_grid_spacing) then
        taum = taum + HF%mu*diff_bnd_m(HF%v(:,i),HF%fd2)/HF%hy(i)
        taup = taup + HF%mu*diff_bnd_p(HF%v(:,i),HF%fd2)/HF%hy(i)
    end if

    if(HF%use_mms) call mms_stresses(x,y,t,taum,taup)
    !if(HF%use_mms) p = mms_hydrofrac(x,y,t,0,'p')
    !if(HF%use_mms) then
    ! taum = mms_hydrofrac(x,y,t,0,'sxy')
    ! taup = mms_hydrofrac(x,y,t,0,'sxy')
    !end if


  end subroutine fluid_stresses

  subroutine update_fields_hydrofrac_implicit(HF,m,p,Zm,Zp,phip,vtp,vtm,dt)
  
    implicit none

    type(HF_type),intent(inout) :: HF
    integer,intent(in) :: m,p
    real,dimension(m:p),intent(in) :: Zm,Zp,phip, & ! P-wave impedances on minus,plus sides
                                      vtp,vtm ! tangential solid velocities on
                                              ! plus, minus sides
    real,intent(in) :: dt

    integer :: i
    real :: Z,Zi

    if (.not.HF%operator_splitting) return ! return if fully explicit

    ! implicit update for pressure
    
    if (HF%coupled .and. HF%linearized_walls .and. .not. HF%slope) then
       do i = m,p
          Zi = (Zp(i) + Zm(i) )/(Zp(i)*Zm(i)) ! combine impedances of two sides
          HF%p(i) = (HF%p(i) -  dt*HF%K0*phip(i)/(HF%wp0(i)-HF%wm0(i)) ) &
                    /(1d0+(dt*HF%K0*Zi)/(HF%wp0(i)-HF%wm0(i))) ! backward Euler
       end do
    end if

    ! Implicit update with geometric correction

    if (HF%coupled .and. HF%linearized_walls .and. HF%slope) then
       do i = m,p
          Zi = (Zp(i) + Zm(i) )/(Zp(i)*Zm(i)) ! combine impedances of two sides
          HF%p(i) = (HF%p(i) -  dt*HF%K0*& 
                        (       &
                            phip(i) - vtp(i)*HF%dwp0dx(i) + vtm(i)*HF%dwm0dx(i) &
                        )       &              
                        /(HF%wp0(i)-HF%wm0(i)) &
                    ) &
                    /(1d0+(dt*HF%K0*Zi)/(HF%wp0(i)-HF%wm0(i))) ! backward Euler
       end do
    end if
    
    if (HF%coupled .and. .not. HF%linearized_walls) then
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

    call update_fields_full_matrix(HF,m,p,dt,vtp,vtm)
    call update_fields_banded_matrix(HF,m,p,dt,vtp,vtm)

  end subroutine update_fields_hydrofrac_implicit

  subroutine update_fields_full_matrix(HF,m,p,dt,vtp,vtm)

    implicit none

    type(hf_type),intent(inout) :: HF
    integer,intent(in) :: m,p
    real,intent(in) :: dt
    real,dimension(m:p),intent(in) :: vtp,vtm 
    
    real,dimension(HF%n) :: b ! Right hand side vector containing boundary terms 
                              ! and solution at previous time step
    real,dimension(HF%n,HF%n) :: A ! matrix to LU decompose
    integer :: i,j,info,ipiv(HF%n),n

    n = HF%n

    if(HF%banded_storage) return

    if(HF%variable_grid_spacing) then
    do i = m,p
       b = HF%v(:,i)
       b(1) = b(1)       - dt*HF%SATym(2)*HF%dydetai(1)*HF%fd2%H00i*vtm(i)/HF%hy(i)
       b(HF%n) = b(HF%n) - dt*HF%SATyp(2)*HF%dydetai(HF%n)*HF%fd2%H00i*vtp(i)/Hf%hy(i)
       call diff_T_bnd_p(b,-dt*HF%dydetai(HF%n)**2*HF%SATyp(3)*vtp(i)/HF%hy(i)**2,HF%fd2)
       call diff_T_bnd_m(b,-dt*HF%dydetai(1)**2   *HF%SATym(3)*vtm(i)/HF%hy(i)**2,HF%fd2)
       A =  - dt*HF%fd2%D2/HF%hy(i)**2
       ! Add identity matrix
       do j=1,HF%n
       A(j,j) = A(j,j) + 1
       end do
       call dgesv(HF%n,1,A,HF%n,ipiv,b,HF%n,info)
       HF%v(:,i) = b
    end do
    end if
    
    if(.not. HF%variable_grid_spacing) then
    do i = m,p
       b = HF%v(:,i)
       b(1) = b(1)  - dt*HF%SATym(2)*HF%fd2%H00i*vtm(i)/HF%hy(i)
       b(n) = b(n)  - dt*HF%SATyp(2)*HF%fd2%H00i*vtp(i)/HF%hy(i)
       
       call diff_T_bnd_m(b,-dt*HF%SATym(3)*vtm(i)/HF%hy(i)**2,HF%fd2)
       call diff_T_bnd_p(b,-dt*HF%SATyp(3)*vtp(i)/HF%hy(i)**2,HF%fd2)

       A =  - dt*(HF%fd2%D2/HF%hy(i)**2)
       !A(1,1) = A(1,1) -dt*HF%SATym(1)*HF%fd2%H00i/HF%hy(i)**2
       !A(n,n) = A(n,n) -dt*HF%SATyp(1)*HF%fd2%H00i/HF%hy(i)**2
       ! Add identity matrix
       do j=1,HF%n
       A(j,j) = A(j,j) + 1
       end do
       call dgesv(HF%n,1,A,HF%n,ipiv,b,HF%n,info)
       HF%v(:,i) = b
    end do
    end if

  end subroutine

  subroutine update_fields_banded_matrix(HF,m,p,dt,vtp,vtm)

    implicit none

    type(hf_type),intent(inout) :: HF
    integer,intent(in) :: m,p
    real,intent(in) :: dt
    real,dimension(m:p),intent(in) :: vtp,vtm 
    
    real,dimension(HF%n) :: b ! Right hand side vector containing boundary terms 
                              ! and solution at previous time step
    real,dimension(:,:),allocatable :: A ! matrix to LU decompose
    integer :: i,j,info,ipiv(HF%n),kl,ku,n,ldb

    if(.not. HF%banded_storage) return

    kl = HF%fd2%nbst-1  ! Number of subdiagonals
    ku = HF%fd2%nbnd-1  ! Number of superdiagonals

    n = HF%n
    ldb = n ! Number of equations
    
    
    allocate(A(2*kl + ku +1,n))
    
    if(HF%variable_grid_spacing) then
       do i = m,p
          b = HF%v(:,i)
          b(1) = b(1) - dt*HF%SATym(2)*HF%fd2%H00i*HF%dydetai(1)*vtm(i)/HF%hy(i)
          b(n) = b(n) - dt*HF%SATyp(2)*HF%fd2%H00i*HF%dydetai(n)*vtp(i)/Hf%hy(i)
          call diff_T_bnd_m(b,-dt*HF%SATym(3)*HF%dydetai(1)**2*vtm(i)/HF%hy(i)**2,HF%fd2)
          call diff_T_bnd_p(b,-dt*HF%SATyp(3)*HF%dydetai(n)**2*vtp(i)/HF%hy(i)**2,HF%fd2)
          A =  - dt*HF%fd2%bD2/HF%hy(i)**2
          ! Add identity matrix to banded matrix
          do j=1,n
          A(kl+ku+1,j) = A(kl+ku+1,j) + 1
          end do
          call dgbsv(n,kl,ku,1,A,2*kl+ku+1,ipiv,b,ldb,info)
          HF%v(:,i) = b
       end do
    end if
    
    if(.not. HF%variable_grid_spacing) then
       do i = m,p
          b = HF%v(:,i)
          b(1) = b(1)       - dt*HF%SATym(2)*HF%fd2%H00i*vtm(i)/HF%hy(i)
          b(HF%n) = b(HF%n) - dt*HF%SATyp(2)*HF%fd2%H00i*vtp(i)/Hf%hy(i)
          call diff_T_bnd_m(b,-dt*HF%SATym(3)*vtm(i)/HF%hy(i)**2,HF%fd2)
          call diff_T_bnd_p(b,-dt*HF%SATyp(3)*vtp(i)/HF%hy(i)**2,HF%fd2)
          A =  - dt*HF%fd2%bD2/HF%hy(i)**2
          ! Add identity matrix to banded matrix
          do j=1,n
          A(kl+ku+1,j) = A(kl+ku+1,j) + 1
          end do
          call dgbsv(n,kl,ku,1,A,2*kl+ku+1,ipiv,b,ldb,info)
          HF%v(:,i) = b
       end do
    end if

    deallocate(A)

  end subroutine


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

    ! Communicate initial wall positions
    ! TODO: This communication only needs to happen once 
    ! (unless we use time dependent walls)
    call populate_ghost_cells(C,HF%L%m,HF%L%p,HF%L%nb, &
         HF%bndm,HF%bndp,HF%wm0,HF%direction)
    call populate_ghost_cells(C,HF%L%m,HF%L%p,HF%L%nb, &
         HF%bndm,HF%bndp,HF%wp0,HF%direction)

  end subroutine share_hydrofrac
  
  subroutine Hnorm(Hw)


      use fd_coeff, only : HL, HR, nbst 

      implicit none

      real,intent(inout),dimension(:) :: Hw
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

      ! Width average with variable grid spacing
      if(HF%variable_grid_spacing) then
       do i = HF%L%m,HF%L%p
          HF%u(i) = HF%u(i) + dot_product(HF%v(:,i),1/HF%dydetai*HF%Hw)
       end do
      end if
       
      ! Width average with uniform grid spacing
      if(.not. HF%variable_grid_spacing) then
       do i = HF%L%m,HF%L%p
          HF%u(i) = HF%u(i) + dot_product(HF%v(:,i),HF%Hw)
       end do
      end if

      HF%u = HF%u/(HF%n - 1)

  end subroutine
  
  subroutine set_rates_viscous(HF,m,p,vtm,vtp)

      implicit none

      type(hf_type),intent(inout) :: HF
      integer,intent(in) :: m,p
      real,intent(in) :: vtm(m:p),vtp(m:p)

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
         ! Penalize v using v-hat = u_dot-hat
         ! Two penalties are applied:  
         ! sigma1*(v - vhat) + sigma2*D^T*(v - vhat)
         HF%Dv(1,i) = HF%Dv(1,i) + HF%SATym(2)*HF%fd2%H00i*(HF%v(1,i) - vtm(i))/HF%hy(i)
         gm = HF%SATym(3)*(HF%v(1,i) - vtm(i))/HF%hy(i)**2
         call diff_T_bnd_m(HF%Dv(:,i),gm,HF%fd2)

         ! SAT terms, plus boundary
         HF%Dv(HF%n,i) = HF%Dv(HF%n,i) + HF%SATyp(2)*HF%fd2%H00i*(HF%v(HF%n,i) - vtp(i))/HF%hy(i) 
         gp = HF%SATyp(3)*(HF%v(HF%n,i) - vtp(i))/HF%hy(i)**2
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
  function diff_bnd_p(u,fd2) result(v)


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
      
      Dv(1:fd2%nD1) = Dv(1:fd2%nD1) + fd2%HLi(1:fd2%nD1)*fd2%D1L(1:fd2%nD1)*g

  end subroutine

  ! Compute Dv = Dv + D^T*g on plus boundary (used by penalty terms)
  subroutine diff_T_bnd_p(Dv,g,fd2)
      

      implicit none

      real,dimension(:),intent(inout) :: Dv
      real,intent(in) :: g
      type(fd2_type),intent(in) :: fd2
      
      integer :: n

      n = size(Dv,1)
      
      Dv((n-fd2%nD1+1):n) = Dv((n-fd2%nD1+1):n) + &
      fd2%HRi(1:fd2%nD1)*fd2%D1R(1:fd2%nD1)*g

  end subroutine
  
  ! Initializes stencils for compact, second derivatives and 
  ! compatible first derivative boundary stencils
  subroutine init_fd2(FDmethod,fd2)

      use io, only : error

      implicit none

      character(*),intent(in) :: FDmethod

      integer :: i
      type(fd2_type), intent(inout) :: fd2


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
    allocate(fd2%HLi(fd2%nbst),fd2%HRi(fd2%nbst))

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
                     1.240509259259259d0,0.911689814814815d0,1.013912037037037d0 /)
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
        fd2%HRi = fd2%HLi(fd2%nbst:1:-1)


  end subroutine

  ! Construct full, second derivative matrix of size n x n (used for linear system solve)
  subroutine init_fd2_full(HF,n)

      implicit none

      type(hf_type),intent(inout) :: HF
      integer,intent(in) :: n

      integer :: i


      allocate(HF%fd2%D2(n,n))
      HF%fd2%D2 = 0d0

      ! Put second derivative stencils in matrix
      
      ! Interior
      do i = HF%fd2%nbst+1,n-HF%fd2%nbst
         HF%fd2%D2(i,i+HF%fd2%mI:i+HF%fd2%pI) = HF%fd2%DI
      end do

      ! Left boundary
      do i = 0,HF%fd2%nbst-1
         HF%fd2%D2(1 + i,1:HF%fd2%nbnd) = HF%fd2%DL(0:HF%fd2%nbnd-1, i)
      end do
    
      ! Right boundary
      do i = 0,HF%fd2%nbst-1         
         HF%fd2%D2(n - i,n - HF%fd2%nbnd+1:n ) = HF%fd2%DR(:,-i)
      end do

      ! SAT penalty matrices
      ! D2*v + BLm*(v - g) + BLD*(Dv - g) + D^TBL*(v - g) + D^TBLD*(v - g) 
      
      ! Note that one penalty term is not included in D2 since it does not scale
      ! as 1/h^2
      HF%fd2%D2 = HF%mu/HF%rho0*HF%fd2%D2

      HF%fd2%D2(1,1) = HF%fd2%D2(1,1) + HF%SATym(2)*HF%fd2%H00i
      HF%fd2%D2(n,n) = HF%fd2%D2(n,n) + HF%SATyp(2)*HF%fd2%H00i

      call diff_T_bnd_m(HF%fd2%D2(:,1)   ,HF%SATym(3),HF%fd2)
      call diff_T_bnd_p(HF%fd2%D2(:,HF%n),HF%SATyp(3),HF%fd2)

      ! Pack the matrix into a banded matrix if requested
      if(HF%banded_storage) then
          call banded_matrix(HF%fd2%D2,HF%fd2%nbnd-1,HF%fd2%nbst-1,HF%fd2%bD2)
          !deallocate(HF%fd2%D2)
          print *,"Banded storage", size(HF%fd2%bD2,1), size(HF%fd2%bD2,2) 
      end if

      ! Used for debugging
      !call print_matrix(HF%fd2%bD2)

  end subroutine

  ! This routine initializes the full matrix for the 
  ! second derivative SBP operator with variable coefficients
  subroutine init_fd2var_full(HF,n)

      use io, only: error
      implicit none

      type(hf_type),intent(inout) :: HF
      integer,intent(in) :: n

      integer :: i,j

      real,dimension(:,:),allocatable :: DS
      logical,parameter :: debug = .false.


      allocate(HF%fd2%D2(n,n),DS(n,n))
      HF%fd2%D2 = 0d0             

      ! The D2 matrix with variable coefficients is
      ! D2 = HI*(-M + DS)
      select case(HF%fd2%FDmethod)
      case('SBP4')
      call init_fd2var4_R(HF%fd2%D2,HF%dydetai,n)
      print *,"FD2VAR: SBP4"
      case('SBP6')
      call init_fd2var6_R(HF%fd2%D2,HF%dydetai,n)
      print *,"FD2VAR: SBP6"
      case('SBP2')
          call error('SBP2 not yet implemented!','init_fd2var_full') 
      end select

      ! Computes: -Hi*M
      HF%fd2%D2 = -HF%fd2%D2           
      do i=1,HF%fd2%nbst
        HF%fd2%D2(i,:) = HF%fd2%HLi(i)*HF%fd2%D2(i,:)
      end do
      j = 1
      do i=n-HF%fd2%nbst+1,n
        HF%fd2%D2(i,:) = HF%fd2%HRi(j)*HF%fd2%D2(i,:)
        j = j + 1
      end do
      
      if(debug) call write_matrix(HF%fd2%D2,'debug/Mi.txt')

      ! Computes: HI*DS
      DS = 0d0
      DS(1,1:HF%fd2%nD1) = -HF%fd2%D1L(1:HF%fd2%nD1)*HF%dydetai(1)*HF%fd2%HLi(1)
      DS(n,n-HF%fd2%nD1+1:n) = &
      HF%fd2%D1R(1:HF%fd2%nD1)*HF%dydetai(n)*HF%fd2%HLi(1)
      if(debug) call write_matrix(DS,'debug/DS.txt')

      HF%fd2%D2 = HF%fd2%D2 + DS

      if(debug) call write_matrix(HF%fd2%D2,'debug/D2.txt')


      ! Multiply by dydetai
       do i=1,HF%n
         HF%fd2%D2(i,:) = HF%fd2%D2(i,:)*HF%dydetai(i)
       end do

      ! SAT penalty matrices
      ! D2*v + BLm*(v - g) + BLD*(Dv - g) + D^TBL*(v - g) + D^TBLD*(v - g) 

      ! Starting by implementing velocities in the fluid only 
      !gm = *( HF%v(1,i)  - vtp(i) )/HF%hy(i)
      HF%fd2%D2 = HF%mu/HF%rho0*HF%fd2%D2
      ! Note that one penalty term is not included in D2 since it does not scale
      ! as 1/h^2  (looks like it should scale like that, why?)
      HF%fd2%D2(1,1) = HF%fd2%D2(1,1) + HF%SATym(2)*HF%fd2%H00i*HF%dydetai(1)
      HF%fd2%D2(n,n) = HF%fd2%D2(n,n) + HF%SATyp(2)*HF%fd2%H00i*HF%dydetai(n)

      call diff_T_bnd_m(HF%fd2%D2(:,1)   ,HF%dydetai(1)**2*HF%SATym(3),HF%fd2)
      call diff_T_bnd_p(HF%fd2%D2(:,HF%n),HF%dydetai(n)**2*HF%SATyp(3),HF%fd2)

      ! Pack the matrix into a banded matrix if requested
      if(HF%banded_storage) then
          call banded_matrix(HF%fd2%D2,HF%fd2%nbnd-1,HF%fd2%nbst-1,HF%fd2%bD2)
          !deallocate(HF%fd2%D2)
          deallocate(DS)
          print *,"Banded storage", size(HF%fd2%bD2,1), size(HF%fd2%bD2,2) 
      end if

  end subroutine

  ! Store the full matrix D2 constructed in "init_fd2_full" using
  ! LAPACK's banded matrix format

subroutine banded_matrix(A,ku,kl,B)

     implicit none

     real(kind=8),dimension(:,:),intent(in) :: A
     real(kind=8),dimension(:,:),allocatable,intent(out) :: B

     real(kind=8),dimension(:),allocatable :: d
     integer :: nbnd,nbst,n,i,ku,kl

    !Extract upper diagonals from full matrix
    nbnd = ku+1
    nbst = kl+1
    n = size(A,1)
    
    allocate(B(2*kl + ku + 1,n))
    allocate(d(n))
    B = 0d0

    
    do i=1,nbnd
      call extract_diagonal(A,1,nbnd-i+1,d)
      B(kl+i,:) = d
    end do
    
    !Extract Lower diagonals from full matrix
    do i=1,nbst-1
      call extract_diagonal(A,i+1,1,d)
      B(kl+i+nbnd,:) = d
    end do  

    deallocate(d)

  end subroutine

  ! Extract a diagonal from a matrix
  subroutine extract_diagonal(A,i,j,d)

   implicit none

   real(kind=8),dimension(:,:),intent(in) :: A
   integer,intent(in) :: i,j
   real(kind=8),dimension(:),intent(out) :: d

   integer :: nx,ny,k

   ny = size(A,1)
   nx = size(A,2)
   
   d = 0d0
   do k=0,ny-1
       if(j+k <= nx .and. i+k <= ny) then
         d(k+j) = a(i+k,j+k)
     end if
   end do
 
  end subroutine

subroutine print_matrix(A)
      implicit none

      real,dimension(:,:),intent(in) :: A
      integer :: i

      do i=1,size(A,1) 
       print '(20G12.4)' , A(i,:)
      end do
end subroutine 

subroutine write_matrix(A,file_name)
    implicit none

      real,dimension(:,:),intent(in) :: A
      character(*),intent(in) :: file_name
      integer :: i,j,N
      N = size(A,1)
open(unit=2, file=file_name, ACTION="write", STATUS="replace")
do i=1,N
    write(2,"(100g15.5)") ( A(i,j), j=1,n )
end do
close(2)

end subroutine

subroutine init_fd2var4_R(R,c,m)

      implicit none

      real,dimension(:,:),intent(inout) :: R
      real,dimension(:),intent(in),allocatable :: c
      integer,intent(in) :: m

      integer :: i

      do i=5,m-4
 R(i,i-2:i+2)=(/-c(i-1)/0.6d1+c(i-2)/0.8d1+c(i)/0.8d1,&
                -c(i-2)/0.6d1-c(i+1)/0.6d1-c(i-1)/0.2d1-c(i)/0.2d1,&
                 c(i-2)/0.24d2+0.5d1/0.6d1*c(i-1)+0.5d1/0.6d1*c(i+1)+c(i+2)/0.24d2+0.3d1/0.4d1*c(i),&
                -c(i-1)/0.6d1- c(i+2)/0.6d1-c(i)/0.2d1-c(i+1)/0.2d1,&
                -c(i+1)/0.6d1+c(i)/0.8d1+c(i+2)/0.8d1 /) 
      end do


 R(1,1:6)=(/0.12d2/0.17d2*&
c(1)+&
0.59d2/0.192d3*&
c(2)+&
0.27010400129d11/0.345067064608d12*&
c(3)+&
0.69462376031d11/0.2070402387648d13*&
c(4),&
-0.59d2/0.68d2*&
c(1)-0.6025413881d10/0.21126554976d11*&
c(3)-&
0.537416663d9/0.7042184992d10*&
c(4),&
0.2d1/0.17d2*&
c(1)-0.59d2/0.192d3*&
c(2)+&
0.213318005d9/0.16049630912d11*&
c(4)+&
0.2083938599d10/0.8024815456d10*&
c(3),&
0.3d1/0.68d2*&
c(1)-0.1244724001d10/0.21126554976d11*&
c(3)+&
0.752806667d9/0.21126554976d11*&
c(4),&
0.49579087d8/0.10149031312d11*&
c(3)-0.49579087d8/0.10149031312d11*&
c(4),&
-c(4)/0.784d3+&
c(3)/0.784d3/)

 R(1+1,1:6) = (/ &

    -0.59d2/0.68d2*&
c(1)-0.6025413881d10/0.21126554976d11*&
c(3)-0.537416663d9/0.7042184992d10*&
c(4),&
0.3481d4/0.3264d4*&
c(1)+&
0.9258282831623875d16/0.7669235228057664d16*&
c(3)+&
0.236024329996203d15/0.1278205871342944d16*&
c(4),&
-0.59d2/0.408d3*&
c(1)-0.29294615794607d14/0.29725717938208d14*&
c(3)-0.2944673881023d13/0.29725717938208d14*&
c(4),&
-0.59d2/0.1088d4*&
c(1)+&
0.260297319232891d15/0.2556411742685888d16*&
c(3)-0.60834186813841d14/0.1278205871342944d16*&
c(4),&
-0.1328188692663d13/0.37594290333616d14*&
c(3)+&
0.1328188692663d13/0.37594290333616d14*&
c(4),&
-0.8673d4/0.2904112d7*&
c(3)+&
0.8673d4/0.2904112d7*&
c(4)/)

 R(1+2,1:6) = (/ &

    0.2d1/0.17d2*&
c(1)-0.59d2/0.192d3*&
c(2)+&
0.213318005d9/0.16049630912d11*&
c(4)+&
0.2083938599d10/0.8024815456d10*&
c(3),&
-0.59d2/0.408d3*&
c(1)-0.29294615794607d14/0.29725717938208d14*&
c(3)-0.2944673881023d13/0.29725717938208d14*&
c(4),&
c(1)/0.51d2+&
0.59d2/0.192d3*&
c(2)+&
0.13777050223300597d17/0.26218083221499456d17*&
c(4)+&
0.564461d6/0.13384296d8*&
c(5)+&
0.378288882302546512209d21/0.270764341349677687456d21*&
c(3),&
c(1)/0.136d3-0.125059d6/0.743572d6*&
c(5)-0.4836340090442187227d19/0.5525802884687299744d19*&
c(3)-0.17220493277981d14/0.89177153814624d14*&
c(4),&
-0.10532412077335d14/0.42840005263888d14*&
c(4)+&
0.1613976761032884305d19/0.7963657098519931984d19*&
c(3)+&
0.564461d6/0.4461432d7*&
c(5),&
-0.960119d6/0.1280713392d10*&
c(4)-0.3391d4/0.6692148d7*&
c(5)+&
0.33235054191d11/0.26452850508784d14*&
c(3)/)
 R(1+3,1:6) = (/ &

    0.3d1/0.68d2*&
c(1)-0.1244724001d10/0.21126554976d11*&
c(3)+&
0.752806667d9/0.21126554976d11*&
c(4),&
-0.59d2/0.1088d4*&
c(1)+&
0.260297319232891d15/0.2556411742685888d16*&
c(3)-0.60834186813841d14/0.1278205871342944d16*&
c(4),&
c(1)/0.136d3-0.125059d6/0.743572d6*&
c(5)-0.4836340090442187227d19/0.5525802884687299744d19*&
c(3)-0.17220493277981d14/0.89177153814624d14*&
c(4),&
0.3d1/0.1088d4*&
c(1)+&
0.507284006600757858213d21/0.475219048083107777984d21*&
c(3)+&
0.1869103d7/0.2230716d7*&
c(5)+&
c(6)/0.24d2+&
0.1950062198436997d16/0.3834617614028832d16*&
c(4),&
-0.4959271814984644613d19/0.20965546238960637264d20*&
c(3)-c(6)/0.6d1-0.15998714909649d14/0.37594290333616d14*&
c(4)-0.375177d6/0.743572d6*&
c(5),&
-0.368395d6/0.2230716d7*&
c(5)+&
0.752806667d9/0.539854092016d12*&
c(3)+&
0.1063649d7/0.8712336d7*&
c(4)+&
c(6)/0.8d1/)
 R(1+4,1:6) = (/ &

    0.49579087d8 / 0.10149031312d11*&
c(3)-0.49579087d8/0.10149031312d11*&
c(4),&
-0.1328188692663d13/0.37594290333616d14*&
c(3)&
+&
 0.1328188692663d13/0.37594290333616d14*&
c(4),&
-0.10532412077335d14/0.42840005263888d14*&
c(4)+&
0.1613976761032884305d19/0.7963657098519931984d19*&
c(3)+&
0.564461d6/0.4461432d7*&
c(5),&
-0.4959271814984644613d19/0.20965546238960637264d20&
*&
 c(3) - c(6) / 0.6d1 - 0.15998714909649d14 / 0.37594290333616d14 *&
 c(4) - 0.375177d6 / 0.743572d6 *&
 c(5) ,& 
 0.8386761355510099813d19 / 0.128413970713633903242d21 *&
 c(3) +&
 0.2224717261773437d16 / 0.2763180339520776d16 *&
 c(4) +&
 0.5d1&
/ 0.6d1 *&
 c(6) +&
 c(7) / 0.24d2 +&
 0.280535d6 / 0.371786d6 *&
 c(5) ,& 
 -0.35039615d8 / 0.213452232d9 *&
 c(4) - c(7) / 0.6d1 - 0.13091810925d11 / 0.13226425254392d14 *&
 c(3) - 0.1118749d7 / 0.2230716d7 *&
 c(5) - c(6) / 0.2d1/)
 R(1+5,1:6) = (/ & 
    -c(4)/0.784d3+&
c(3)/0.784d3,&
-0.8673d4/0.2904112d7*&
c(3)+&
0.8673d4/0.2904112d7*&
c(4),&
-0.960119d6/0.1280713392d10*&
c(4)-0.3391d4/0.6692148d7*&
c(5)+&
0.33235054191d11/0.26452850508784d14*&
c(3),&
-0.368395d6/0.2230716d7*&
c(5)+&
0.752806667d9/0.539854092016d12*&
c(3)+&
0.1063649d7/0.8712336d7*&
c(4)+&
c(6)/0.8d1,&
-0.35039615d8&
/ 0.213452232d9 *&
 c(4) - c(7) / 0.6d1 - 0.13091810925d11 / 0.13226425254392d14 *&
 c(3) - 0.1118749d7 / 0.2230716d7 *&
 c(5) - c(6) / 0.2d1 ,& 
 0.3290636d7 / 0.80044587d8 *&
 c(4) +&
 0.5580181d7 / 0.6692148d7 *&
 c(5) +&
 0.5d1 / 0.6d1 *&
 c(7) +&
 c(8) / 0.24d2 +&
 0.660204843d9 / 0.13226425254392d14 *&
 c(3) +&
 0.3d1 / 0.4d1 *&
 c(6)/)


R(m-5,m-5:m)=(/c(m-7)/0.24d2+&
0.5d1/0.6d1 *&
 c(m-6) +&
 0.5580181d7 / 0.6692148d7&
*&
 c(m-4) +&
 0.4887707739997d13 / 0.119037827289528d15 *&
 c(m-3) +&
 0.3d1 / 0.4d1&
*&
 c(m-5) +&
 0.660204843d9 / 0.13226425254392d14 *&
 c(m-2) +&
 0.660204843d9 / 0.13226425254392d14&
*&
 c(m-1) ,& 
 -c(m-6) / 0.6d1 - 0.1618585929605d13 / 0.9919818940794d13 *&
 c(m-3) - c(m-5) / 0.2d1&
- 0.1118749d7 / 0.2230716d7 *&
 c(m-4) - 0.13091810925d11 / 0.13226425254392d14 *&
 c(m-2) - 0.13091810925d11&
/ 0.13226425254392d14 *&
 c(m-1) ,& 
 -0.368395d6 / 0.2230716d7 *&
 c(m-4) +&
 c(m-5) / 0.8d1 +&
 0.48866620889d11&
/ 0.404890569012d12 *&
 c(m-3) +&
 0.752806667d9 / 0.539854092016d12 *&
 c(m-2) +&
 0.752806667d9 / 0.539854092016d12&
*&
 c(m-1) ,& 
 -0.3391d4 / 0.6692148d7 *&
 c(m-4) - 0.238797444493d12 / 0.119037827289528d15 *&
 c(m-3) +&
 0.33235054191d11&
/ 0.26452850508784d14 *&
 c(m-2) +&
 0.33235054191d11 / 0.26452850508784d14 *&
 c(m-1) ,& 
 -0.8673d4 / 0.2904112d7 *&
 c(m-2)&
- 0.8673d4 / 0.2904112d7 *&
 c(m-1) +&
 0.8673d4 / 0.1452056d7 *&
 c(m-3) ,& 
 -c(m-3) / 0.392d3 +&
 c(m-2) / 0.784d3 +&
 c(m-1)&
/ 0.784d3/)
 R(m-5+1,m-5:m) = (/ & 
    -c(m-6) / 0.6d1 - 0.1618585929605d13 / 0.9919818940794d13 *&
 c(m-3) - c(m-5)&
/ 0.2d1 - 0.1118749d7 / 0.2230716d7 *&
 c(m-4) - 0.13091810925d11 / 0.13226425254392d14&
*&
 c(m-2) - 0.13091810925d11 / 0.13226425254392d14 *&
 c(m-1) ,& 
 c(m-6) / 0.24d2 +&
0.5d1 / 0.6d1 *&
 c(m-5) +&
 0.3896014498639d13 / 0.4959909470397d13 *&
 c(m-3) +&
 0.8386761355510099813d19&
/ 0.128413970713633903242d21 *&
 c(m-2) +&
 0.280535d6 / 0.371786d6 *&
 c(m-4) +&
 0.3360696339136261875d19 /&
0.171218627618178537656d21 *&
 c(m-1) ,& 
 -c(m-5) / 0.6d1 - 0.4959271814984644613d19 / 0.20965546238960637264d20&
*&
 c(m-2) - 0.375177d6 / 0.743572d6 *&
 c(m-4) - 0.13425842714d11 / 0.33740880751d11 *&
 c(m-3) - 0.193247108773400725d18&
/ 0.6988515412986879088d19 *&
 c(m-1) ,& 
 -0.365281640980d12 / 0.1653303156799d13 *&
 c(m-3) +&
 0.564461d6 / 0.4461432d7&
*&
 c(m-4) +&
 0.1613976761032884305d19 / 0.7963657098519931984d19 *&
 c(m-2) - 0.198407225513315475d18 / 0.7963657098519931984d19&
*&
 c(m-1) ,& 
 -0.1328188692663d13 / 0.37594290333616d14 *&
 c(m-2) +&
 0.2226377963775d13 / 0.37594290333616d14 *&
 c(m-1) - 0.8673d4&
/ 0.363014d6 *&
 c(m-3) ,& 
 c(m-3) / 0.49d2 +&
 0.49579087d8 / 0.10149031312d11 *&
 c(m-2) - 0.256702175d9 / 0.10149031312d11 *&
 c(m-1)/)
 R(m-5+2,m-5:m) = (/ &

    -0.368395d6 / 0.2230716d7 *&
 c(m-4) +&
 c(m-5) / 0.8d1 +&
 0.48866620889d11&
/ 0.404890569012d12 *&
 c(m-3) +&
 0.752806667d9 / 0.539854092016d12 *&
 c(m-2) +&
0.752806667d9 / 0.539854092016d12 *&
 c(m-1) ,& 
 -c(m-5) / 0.6d1 - 0.4959271814984644613d19&
/ 0.20965546238960637264d20 *&
 c(m-2) - 0.375177d6 / 0.743572d6 *&
 c(m-4) - 0.13425842714d11&
/ 0.33740880751d11 *&
 c(m-3) - 0.193247108773400725d18 / 0.6988515412986879088d19 *&
 c(m-1) ,&

 c(m-5) / 0.24d2 +&
 0.1869103d7 / 0.2230716d7 *&
 c(m-4) +&
 0.507284006600757858213d21 / 0.475219048083107777984d21&
*&
 c(m-2) +&
 0.3d1 / 0.1088d4 *&
 c(m) +&
 0.31688435395d11 / 0.67481761502d11 *&
 c(m-3) +&
 0.27769176016102795561d20&
/ 0.712828572124661666976d21 *&
 c(m-1) ,& 
 -0.125059d6 / 0.743572d6 *&
 c(m-4) +&
 c(m) / 0.136d3 - 0.23099342648d11 / 0.101222642253d12&
*&
 c(m-3) - 0.4836340090442187227d19 / 0.5525802884687299744d19 *&
 c(m-2) +&
 0.193950157930938693d18 / 0.5525802884687299744d19 *&
 c(m-1)&
,& 
 0.260297319232891d15 / 0.2556411742685888d16 *&
 c(m-2) - 0.59d2 / 0.1088d4 *&
 c(m) - 0.106641839640553d15 / 0.1278205871342944d16 *&
 c(m-1)&
+&
 0.26019d5 / 0.726028d6 *&
 c(m-3) ,& 
 -0.1244724001d10 / 0.21126554976d11 *&
 c(m-2) +&
 0.3d1 / 0.68d2 *&
 c(m) +&
 0.752806667d9 / 0.21126554976d11&
*&
 c(m-1)/)
 R(m-5+3,m-5:m) = (/ & 
    -0.3391d4 / 0.6692148d7 *&
 c(m-4) - 0.238797444493d12 / 0.119037827289528d15&
*&
 c(m-3) +&
 0.33235054191d11 / 0.26452850508784d14 *&
 c(m-2) +&
 0.33235054191d11&
/ 0.26452850508784d14 *&
 c(m-1) ,& 
 -0.365281640980d12 / 0.1653303156799d13 *&
 c(m-3)&
+&
 0.564461d6 / 0.4461432d7 *&
 c(m-4) +&
 0.1613976761032884305d19 / 0.7963657098519931984d19&
*&
 c(m-2) - 0.198407225513315475d18 / 0.7963657098519931984d19 *&
 c(m-1) ,& 
 -0.125059d6 / 0.743572d6&
*&
 c(m-4) +&
 c(m) / 0.136d3 - 0.23099342648d11 / 0.101222642253d12 *&
 c(m-3) - 0.4836340090442187227d19&
/ 0.5525802884687299744d19 *&
 c(m-2) +&
 0.193950157930938693d18 / 0.5525802884687299744d19 *&
 c(m-1) ,&

 0.564461d6 / 0.13384296d8 *&
 c(m-4) +&
 0.470299699916357d15 / 0.952302618316224d15 *&
 c(m-3) +&
 0.550597048646198778781d21&
/ 0.1624586048098066124736d22 *&
 c(m-1) +&
 c(m) / 0.51d2 +&
 0.378288882302546512209d21 / 0.270764341349677687456d21 *&
 c(m-2)&
,& 
 -0.59d2 / 0.408d3 *&
 c(m) - 0.29294615794607d14 / 0.29725717938208d14 *&
 c(m-2) - 0.2234477713167d13 / 0.29725717938208d14&
*&
 c(m-1) - 0.8673d4 / 0.363014d6 *&
 c(m-3) ,& 
 -0.59d2 / 0.3136d4 *&
 c(m-3) - 0.13249937023d11 / 0.48148892736d11 *&
 c(m-1)&
+&
 0.2d1 / 0.17d2 *&
 c(m) +&
 0.2083938599d10 / 0.8024815456d10 *&
 c(m-2)/)
 R(m-5+4,m-5:m) = (/ & 
    -0.8673d4 / 0.2904112d7 *&
 c(m-2) - 0.8673d4 / 0.2904112d7 *&
 c(m-1) +&
 0.8673d4&
/ 0.1452056d7 *&
 c(m-3) ,& 
 -0.1328188692663d13 / 0.37594290333616d14 *&
 c(m-2) +&
 0.2226377963775d13&
/ 0.37594290333616d14 *&
 c(m-1) - 0.8673d4 / 0.363014d6 *&
 c(m-3) ,& 
 0.260297319232891d15 / 0.2556411742685888d16&
*&
 c(m-2) - 0.59d2 / 0.1088d4 *&
 c(m) - 0.106641839640553d15 / 0.1278205871342944d16 *&
 c(m-1) +&
 0.26019d5 / 0.726028d6&
*&
 c(m-3) ,& 
 -0.59d2 / 0.408d3 *&
 c(m) - 0.29294615794607d14 / 0.29725717938208d14 *&
 c(m-2) - 0.2234477713167d13 / 0.29725717938208d14&
*&
 c(m-1) - 0.8673d4 / 0.363014d6 *&
 c(m-3) ,& 
 0.9258282831623875d16 / 0.7669235228057664d16 *&
 c(m-2) +&
 0.3481d4 / 0.3264d4 *&
 c(m)&
+&
 0.228389721191751d15 / 0.1278205871342944d16 *&
 c(m-1) +&
 0.8673d4 / 0.1452056d7 *&
 c(m-3) ,& 
 -0.6025413881d10 / 0.21126554976d11 *&
c(m-2) - 0.59d2 / 0.68d2 *&
 c(m) - 0.537416663d9 / 0.7042184992d10 *&
 c(m-1)/)
 R(m-5+5,m-5:m) = (/ &
    -c(m-3) / 0.392d3 +&
 c(m-2) / 0.784d3 +&
 c(m-1) / 0.784d3 ,& 
 c(m-3) / 0.49d2&
+&
 0.49579087d8 / 0.10149031312d11 *&
 c(m-2) - 0.256702175d9 / 0.10149031312d11 *&
c(m-1) ,& 
 -0.1244724001d10 / 0.21126554976d11 *&
 c(m-2) +&
 0.3d1 / 0.68d2 *&
 c(m)&
+&
 0.752806667d9 / 0.21126554976d11 *&
 c(m-1) ,& 
 -0.59d2 / 0.3136d4 *&
 c(m-3) - 0.13249937023d11&
/ 0.48148892736d11 *&
 c(m-1) +&
 0.2d1 / 0.17d2 *&
 c(m) +&
 0.2083938599d10 / 0.8024815456d10 *&
 c(m-2)&
,& 
 -0.6025413881d10 / 0.21126554976d11 *&
 c(m-2) - 0.59d2 / 0.68d2 *&
 c(m) - 0.537416663d9 / 0.7042184992d10&
*&
 c(m-1) ,& 
 0.3d1 / 0.3136d4 *&
 c(m-3) +&
 0.27010400129d11 / 0.345067064608d12 *&
 c(m-2) +&
 0.234566387291d12&
/ 0.690134129216d12 *&
 c(m-1) +&
 0.12d2 / 0.17d2 *&
 c(m)/)

  end subroutine !init_fd2var6_R

  subroutine init_fd2var6_R(R,c,m)

      implicit none

      real,dimension(:,:),intent(inout) :: R
      real,dimension(:),intent(in),allocatable :: c
      integer,intent(in) :: m

      integer :: i


      do i=6,m-5
R(i,i-3:i+3)=(/c(i-2) / 0.40d2 + c(i-1) / 0.40d2 - 0.11d2 / 0.360d3 *&
 c(i-3)- 0.11d2 / 0.360d3 *c(i) ,& 
 c(i-3) / 0.20d2 - 0.3d1 / 0.10d2 *c(i-1) +&
 c(i+1)/ 0.20d2 + 0.7d1 / 0.40d2 *c(i) + 0.7d1 / 0.40d2 * c(i-2) ,& 
 -c(i-3) / 0.40d2 -&
0.3d1 / 0.10d2 *&
 c(i-2) - 0.3d1 / 0.10d2 *&
 c(i+&
1) - c(i+&
2) / 0.40d2 - 0.17d2 / 0.40d2&
*&
 c(i) - 0.17d2 / 0.40d2 *&
 c(i-1) ,& 
 c(i-3) / 0.180d3 +&
 c(i-2) / 0.8d1 +&
 0.19d2 / 0.20d2&
*&
 c(i-1) +&
 0.19d2 / 0.20d2 *&
 c(i+&
1) +&
 c(i+&
2) / 0.8d1 +&
 c(i+&
3) / 0.180d3 +&
 0.101d3&
/ 0.180d3 *&
 c(i) ,& 
 -c(i-2) / 0.40d2 - 0.3d1 / 0.10d2 *&
 c(i-1) - 0.3d1 / 0.10d2 *&
 c(i+&
2) -&
c(i+&
3) / 0.40d2 - 0.17d2 / 0.40d2 *&
 c(i) - 0.17d2 / 0.40d2 *&
 c(i+&
1) ,& 
 c(i-1) / 0.20d2 - 0.3d1&
/ 0.10d2 *&
 c(i+&
1) +&
 c(i+&
3) / 0.20d2 +&
 0.7d1 / 0.40d2 *&
 c(i) +&
 0.7d1 / 0.40d2 *&
 c(i+&
2) ,&

 c(i+&
1) / 0.40d2 +&
 c(i+&
2) / 0.40d2 - 0.11d2 / 0.360d3 *&
 c(i) - 0.11d2 / 0.360d3 *&
 c(i+&
3)/)
     end do

R(1,1:9)=(/0.7912667594695582093926295d0 *&
 c(1) +&
 0.2968472090638000742888467d0&
*&
 c(2) +&
 0.3185519088796429015220016d-2 *&
 c(3) +&
 0.1632404042590951953384672d-1&
*&
 c(4) +&
 0.3160302244094415087693968d-1 *&
 c(5) +&
 0.3167964748016105299646518d-1&
*&
 c(6) +&
 0.3148577733947253920469418d-1 *&
 c(7) ,& 
 -0.1016689339350338144430605d1&
*&
 c(1) - 0.2845627370491611369031341d-1 *&
 c(3) - 0.4128029838349298819825156d-1 *&
 c(4)&
- 0.1392281451620140507549866d0 *&
 c(5) - 0.1195777325611201766551392d0 *&
 c(6) - 0.1194267756529333410855186d0&
*&
 c(7) ,& 
 0.7075642937243715046279337d-1 *&
 c(1) - 0.1845476106024151050283847d0 *&
 c(2) - 0.4364163147111892346990101d-1&
*&
 c(4) +&
 0.2432367907207732460874765d0 *&
 c(5) +&
 0.1582127073537215443965653d0 *&
 c(6) +&
 0.1602348578364786307613271d0&
*&
 c(7) ,& 
 0.2251991532891353212689574d0 *&
 c(1) - 0.1662748711097054895317080d0 *&
 c(2) +&
 0.2710530961648671297733465d-1&
*&
 c(3) - 0.1916646185968439909125616d0 *&
 c(5) - 0.7684117160199014594442072d-1 *&
 c(6) - 0.8219586949831697575883635d-1 *&
c(7) ,& 
 -0.5224403464202056316702078d-1 *&
 c(1) +&
 0.4440063948509876221050939d-1 *&
 c(2) - 0.1023976547309387874453988d-2 *&
c(3) +&
 0.7403484645316174090533193d-1 *&
 c(4) +&
 0.1241625568998496895352046d-1 *&
 c(6) +&
 0.7188652847892601282652228d-1 *&
c(5) +&
 0.1379362997104735503447960d-1 *&
 c(7) ,& 
 -0.1828896813877197352675410d-1 *&
 c(1) +&
 0.9574633163221758060736592d-2 *&
c(2) - 0.8105784530576404277872603d-3 *&
 c(3) - 0.7348845587775519698437916d-2 *&
 c(4) +&
 0.1063601949723906997026904d-1 *&
 c(5) -&
0.1315967038382618382356495d-1 *&
 c(6) - 0.2117936478838753524581943d-1 *&
 c(7) ,& 
 0.1911888563316170927411831d-2 *&
 c(4) - 0.4068130355529149936100229d-1&
*&
 c(5) +&
 0.1319674981073749167009902d-1 *&
 c(6) +&
 0.2557266518123783676349144d-1 *&
 c(7) ,& 
 0.1559652871136785763960685d-1 *&
 c(5) - 0.6486184157331537899459796d-2&
*&
 c(6) - 0.9110344554036319740147054d-2 *&
 c(7) ,& 
 0.5593983696629863059347067d-3 *&
 c(6) - 0.1384822535100796372263822d-2 *&
 c(5) +&
 0.8254241654378100663291154d-3 *&
c(7)/)


 R(1+1,1:9) = (/ -0.1016689339350338144430605d1 *&
 c(1) - 0.2845627370491611369031341d-1 *&
 c(3) - 0.4128029838349298819825156d-1 *&
 c(4) - 0.1392281451620140507549866d0 *&
 c(5) - 0.1195777325611201766551392d0&
*&
 c(6) - 0.1194267756529333410855186d0 *&
 c(7) ,& 
 0.1306332157111667628555907d1 *&
 c(1) +&
 0.2542001760457345743492403d0 *&
 c(3) +&
 0.1043897828092562609502636d0 *&
 c(4) +&
 0.6672328021032112950919876d0&
*&
 c(5) +&
 0.4681819359722749441073885d0 *&
 c(6) +&
 0.4676415410195836920069412d0 *&
 c(7) ,& 
 -0.9091410269992464604926176d-1 *&
 c(1) +&
 0.1103611313171476425250639d0 *&
 c(4) - 0.1290397544997518887000350d1&
*&
 c(5) - 0.6639605248735044787146222d0 *&
 c(6) - 0.6615974464005206184151509d0 *&
 c(7) ,& 
 -0.2893557395653431666593814d0 *&
 c(1) - 0.2421320004064592721552708d0 *&
 c(3) +&
 0.1187670255028031027693374d1 *&
 c(5)&
+&
 0.3956598149904136332753521d0 *&
 c(6) +&
 0.3860048921755800000681479d0 *&
 c(7) ,& 
 0.6712774475803763988977355d-1 *&
 c(1) +&
 0.9147192682075630179962131d-2 *&
 c(3) - 0.1872196143003808021730728d0 *&
 c(4) - 0.1319358558853174530078498d0&
*&
 c(6) - 0.4871575736811911887376923d0 *&
 c(5) - 0.1047516312275448138054418d0 *&
 c(7) ,& 
 0.2349927974590068869356781d-1 *&
 c(1) +&
 0.7240905383565181316381731d-2 *&
 c(3) +&
 0.1858378996391679448655070d-1 *&
 c(4) - 0.9289616133938676174345208d-1&
*&
 c(5) +&
 0.1223513270418807666970488d0 *&
 c(6) +&
 0.1113520320436295033894092d0 *&
 c(7) ,& 
 -0.4834791406446907590553793d-2 *&
 c(4) +&
 0.2310683832687820403062715d0 *&
 c(5) - 0.1080774142196007991746827d0 *&
 c(6) - 0.1181561776427343335410349d0&
*&
 c(7) ,& 
 -0.8368141434403455353724691d-1 *&
 c(5) +&
 0.4093499466767054661591066d-1 *&
 c(6) +&
 0.4274641967636400692133626d-1 *&
 c(7) ,& 
 -0.3576545132696983143406173d-2 *&
 c(6) +&
 0.7389399124121078682094445d-2 *&
 c(5) - 0.3812853991424095538688273d-2&
*&
 c(7)/)
 R(1+2,1:9) = (/ 0.7075642937243715046279337d-1 *&
 c(1) - 0.1845476106024151050283847d0 *&
 c(2) - 0.4364163147111892346990101d-1 *&
 c(4) +&
 0.2432367907207732460874765d0 *&
 c(5) +&
 0.1582127073537215443965653d0 *&
 c(6) +&
 0.1602348578364786307613271d0 *&
 c(7) ,&

 -0.9091410269992464604926176d-1 *&
 c(1) +&
 0.1103611313171476425250639d0 *&
 c(4) - 0.1290397544997518887000350d1 *&
 c(5) - 0.6639605248735044787146222d0 *&
 c(6) - 0.6615974464005206184151509d0 *&
 c(7) ,& 
 0.6327161147136873807796515d-2 *&
 c(1) +&
 0.1147318200715868527529827d0&
*&
 c(2) +&
 0.1166740554279680007487795d0 *&
 c(4) +&
 0.2766610808285444037240703d1 *&
 c(5) +&
 0.1070920689960817104203947d1 *&
 c(6) +&
 0.1013161391032973057171717d1 *&
 c(7) ,& 
 0.2013769413884797246646959d-1 *&
 c(1) +&
 0.1033717994630886401730470d0 *&
 c(2) - 0.2913221621151742724258117d1&
*&
 c(5) - 0.8755807343482262259774782d0 *&
 c(6) - 0.6909957183488812426508351d0 *&
 c(7) ,& 
 -0.4671751091575462868310238d-2 *&
 c(1) - 0.2760353365637712827793337d-1 *&
 c(2) - 0.1979290298620869974478871d0 *&
 c(4) +&
 0.5402985338373433052255418d0 *&
 c(6) +&
 0.1239177593031911077924537d1 *&
c(5) +&
 0.2628038050247358227280031d0 *&
 c(7) ,& 
 -0.1635430866921887819487473d-2 *&
 c(1) - 0.5952475275883259619711594d-2 *&
 c(2) +&
 0.1964682777744275219350831d-1 *&
 c(4) +&
 0.3236640012639046600590714d0 *&
 c(5) - 0.4659516693228870973898560d0 *&
 c(6) - 0.2217272720941736859420432d0 *&
 c(7)&
,& 
 -0.5111353189352474549563559d-2 *&
 c(4) - 0.5355878163774754346032096d0 *&
 c(5) +&
 0.3328335104489738933610597d0 *&
 c(6) +&
 0.2078656591178540157917135d0 *&
 c(7) ,& 
 0.1824328174134289562208038d0 *&
 c(5) - 0.1059816030196818445908057d0 *&
 c(6) - 0.7645121439374711162999809d-1 *&
 c(7) ,& 
0.9209089963443799485648361d-2 *&
 c(6) - 0.1591502818872493167091475d-1 *&
 c(5) +&
 0.6705938225281132185266388d-2 *&
 c(7)/)
 R(1+3,1:9) = (/ 0.2251991532891353212689574d0 *&
 c(1) - 0.1662748711097054895317080d0 *&
 c(2) +&
 0.2710530961648671297733465d-1 *&
 c(3) - 0.1916646185968439909125616d0 *&
 c(5) - 0.7684117160199014594442072d-1&
*&
 c(6) - 0.8219586949831697575883635d-1 *&
 c(7) ,& 
 -0.2893557395653431666593814d0 *&
 c(1) - 0.2421320004064592721552708d0 *&
 c(3) +&
 0.1187670255028031027693374d1 *&
 c(5) +&
 0.3956598149904136332753521d0 *&
 c(6) +&
 0.3860048921755800000681479d0 *&
 c(7) ,& 
 0.2013769413884797246646959d-1 *&
 c(1) +&
 0.1033717994630886401730470d0&
*&
 c(2) - 0.2913221621151742724258117d1 *&
 c(5) - 0.8755807343482262259774782d0 *&
 c(6) - 0.6909957183488812426508351d0 *&
 c(7) ,& 
 0.6409299775987186986730499d-1 *&
 c(1) +&
 0.9313657638804699948929701d-1 *&
 c(2) +&
 0.2306367624634749229113646d0 *&
 c(3) +&
 0.3689440308283716620260816d1 *&
 c(5) +&
 0.1190550338687608873798462d1 *&
c(6) +&
 0.5912479546888856519443605d0 *&
 c(7) ,& 
 -0.1486895819265604128572498d-1 *&
 c(1) - 0.2487040599390160764166412d-1 *&
 c(2) - 0.8712928907711754187084757d-2 *&
 c(3) - 0.1263507837371824205693950d1 *&
 c(6) - 0.3058317397843997326920898d0 *&
 c(7) - 0.1470691926045802954795783d1 *&
 c(5) ,& 
 -0.5205147429855955657625694d-2 *&
 c(1)&
- 0.5363098747528542488971874d-2 *&
 c(2) - 0.6897142765790609546343709d-2 *&
 c(3) - 0.7857524521667450101721993d0 *&
 c(5) +&
 0.2291148005423734600066709d0 *&
 c(7) +&
 0.9977064356292750529201981d0 *&
 c(6) ,& 
 0.6697297488067662265210608d0 *&
 c(5) - 0.5013247356072127938999311d0 *&
 c(6) - 0.1795161243106645437322408d0 *&
 c(7) ,& 
 -0.2022909060111751565150958d0&
*&
 c(5) +&
 0.1453421858063658498587377d0 *&
 c(6) +&
 0.5694872020480930665635812d-1 *&
 c(7) ,& 
 -0.1200429618441003833696998d-1 *&
 c(6) - 0.4776915669385923841535432d-2 *&
 c(7) +&
 0.1678121185379596217850541d-1 *&
 c(5)/)
 R(1+4,1:9) = (/ -0.5224403464202056316702078d-1 *&
 c(1) +&
 0.4440063948509876221050939d-1 *&
 c(2) - 0.1023976547309387874453988d-2 *&
 c(3) +&
 0.7403484645316174090533193d-1&
*&
 c(4) +&
 0.1241625568998496895352046d-1 *&
 c(6) +&
 0.7188652847892601282652228d-1 *&
 c(5) +&
 0.1379362997104735503447960d-1 *&
 c(7) ,& 
 0.6712774475803763988977355d-1 *&
 c(1) +&
 0.9147192682075630179962131d-2 *&
 c(3) - 0.1872196143003808021730728d0 *&
 c(4) - 0.1319358558853174530078498d0 *&
 c(6) - 0.4871575736811911887376923d0 *&
 c(5) - 0.1047516312275448138054418d0 *&
 c(7)&
,& 
 -0.4671751091575462868310238d-2 *&
 c(1) - 0.2760353365637712827793337d-1 *&
 c(2) - 0.1979290298620869974478871d0 *&
 c(4) +&
 0.5402985338373433052255418d0 *&
 c(6) +&
 0.1239177593031911077924537d1 *&
 c(5) +&
 0.2628038050247358227280031d0 *&
 c(7) ,& 
 -0.1486895819265604128572498d-1 *&
 c(1) - 0.2487040599390160764166412d-1 *&
 c(2) - 0.8712928907711754187084757d-2 *&
 c(3) - 0.1263507837371824205693950d1&
*&
 c(6) - 0.3058317397843997326920898d0 *&
 c(7) - 0.1470691926045802954795783d1 *&
 c(5) ,& 
 0.3449455095910233625229891d-2 *&
 c(1) +&
 0.6641183499427826101618457d-2 *&
 c(2) +&
 0.3291545083271862858501887d-3 *&
 c(3) +&
 0.3357721707576477199985656d0 *&
 c(4) +&
 0.2096413329579026439044119d1 *&
 c(6) +&
 0.2317323204183126854954203d0 *&
 c(7) +&
 0.6107825764368264576481962d-2 *&
 c(8) +&
 0.7109125850683376695640722d0&
*&
 c(5) ,& 
 0.1207544072304193806052558d-2 *&
 c(1) +&
 0.1432116665752147607469646d-2 *&
 c(2) +&
 0.2605582646183255957264249d-3 *&
 c(3) - 0.3332941113251635390801278d-1 *&
 c(4) - 0.2808241697385532683612407d0 *&
 c(7) - 0.2720908083525083608370563d-1 *&
 c(8) +&
 0.1045865435682921987447929d0 *&
 c(5) - 0.1348436986667115543203552d1 *&
 c(6) ,& 
 0.8671038084174692625075159d-2 *&
 c(4) +&
 0.1736073411355428563685818d0 *&
c(6) +&
 0.5331362125287625412555844d-1 *&
 c(8) - 0.2424935262404526301801157d0 *&
 c(5) +&
 0.1569015257678588270609004d0 *&
 c(7) ,& 
 -0.8631683980217122275970376d-1 *&
 c(6) +&
 0.2698842360470999243492629d-1 *&
 c(7) +&
 0.8098194147715651085292754d-1 *&
 c(5) - 0.3276463639080639163926118d-1 *&
 c(8) ,& 
 0.7462059484530855073291365d-2 *&
 c(6) - 0.8121640361668678949573496d-3 *&
 c(7) +&
 0.5522702088127090209264064d-3 *&
c(8) - 0.7202165657176696199260422d-2 *&
 c(5)/)
 R(1+5,1:9) = (/ -0.1828896813877197352675410d-1 *&
 c(1) +&
 0.9574633163221758060736592d-2 *&
 c(2) - 0.8105784530576404277872603d-3 *&
 c(3) - 0.7348845587775519698437916d-2 *&
 c(4) +&
 0.1063601949723906997026904d-1 *&
 c(5) - 0.1315967038382618382356495d-1 *&
 c(6) - 0.2117936478838753524581943d-1 *&
 c(7) ,& 
 0.2349927974590068869356781d-1 *&
 c(1) +&
 0.7240905383565181316381731d-2 *&
 c(3)&
+&
 0.1858378996391679448655070d-1 *&
 c(4) - 0.9289616133938676174345208d-1 *&
 c(5) +&
 0.1223513270418807666970488d0 *&
 c(6) +&
 0.1113520320436295033894092d0 *&
 c(7) ,& 
 -0.1635430866921887819487473d-2 *&
 c(1) - 0.5952475275883259619711594d-2 *&
 c(2) +&
 0.1964682777744275219350831d-1 *&
 c(4) +&
 0.3236640012639046600590714d0 *&
 c(5) - 0.4659516693228870973898560d0 *&
 c(6) - 0.2217272720941736859420432d0 *&
 c(7) ,& 
 -0.5205147429855955657625694d-2&
*&
 c(1) - 0.5363098747528542488971874d-2 *&
 c(2) - 0.6897142765790609546343709d-2 *&
 c(3) - 0.7857524521667450101721993d0 *&
 c(5) +&
 0.2291148005423734600066709d0 *&
 c(7) +&
 0.9977064356292750529201981d0 *&
 c(6) ,& 
 0.1207544072304193806052558d-2 *&
 c(1) +&
 0.1432116665752147607469646d-2 *&
 c(2) +&
 0.2605582646183255957264249d-3 *&
 c(3) - 0.3332941113251635390801278d-1 *&
 c(4) - 0.2808241697385532683612407d0 *&
 c(7) - 0.2720908083525083608370563d-1&
*&
 c(8) +&
 0.1045865435682921987447929d0 *&
 c(5) - 0.1348436986667115543203552d1 *&
 c(6) ,& 
 0.4227226173449345042468960d-3 *&
 c(1) +&
 0.3088241944378964404772302d-3 *&
 c(2) +&
 0.2062575706647430620228133d-3 *&
 c(3) +&
 0.3308343404200968256656458d-2 *&
 c(4) +&
 0.5828047016405001815804837d0 *&
 c(5) +&
 0.8054174220366215473556835d0 *&
 c(7) +&
 0.1338363233410033443348225d0 *&
 c(8) +&
 0.5555555555555555555555556d-2 *&
 c(9) +&
 0.1190362071861893051132274d1&
*&
 c(6) ,& 
 -0.8607044252686413302647675d-3 *&
 c(4) - 0.1748074708673904989293256d0 *&
 c(5) - 0.3132544850115050165022338d0 *&
 c(8) - 0.2500000000000000000000000d-1 *&
 c(9) - 0.3169166305310429271303167d0 *&
 c(7) - 0.6691607091647929161078591d0 *&
 c(6) ,& 
 0.3354661791693352108660900d-1 *&
 c(5) - 0.3343620022386971405018586d0 *&
 c(7) +&
 0.5000000000000000000000000d-1 *&
 c(9) +&
 0.2169790609807602750804271d0 *&
 c(6) +&
 0.1838363233410033443348225d0 *&
 c(8)&
,& 
 0.2912518476823004642951502d-1 *&
 c(7) +&
 0.2279091916474916391629437d-1 *&
 c(8) - 0.3068985997518740530511593d-1 *&
 c(6) - 0.1781799513347360596249022d-2 *&
 c(5) - 0.3055555555555555555555556d-1 *&
 c(9)/)
 R(1+6,1:9) = (/ 0.1911888563316170927411831d-2 *&
 c(4) - 0.4068130355529149936100229d-1 *&
 c(5) +&
 0.1319674981073749167009902d-1 *&
 c(6) +&
 0.2557266518123783676349144d-1 *&
 c(7) ,& 
 -0.4834791406446907590553793d-2 *&
 c(4) +&
 0.2310683832687820403062715d0 *&
 c(5)&
- 0.1080774142196007991746827d0 *&
 c(6) - 0.1181561776427343335410349d0 *&
 c(7) ,& 
 -0.5111353189352474549563559d-2 *&
 c(4) - 0.5355878163774754346032096d0 *&
 c(5) +&
 0.3328335104489738933610597d0 *&
 c(6) +&
 0.2078656591178540157917135d0 *&
 c(7) ,& 
 0.6697297488067662265210608d0 *&
 c(5) - 0.5013247356072127938999311d0 *&
 c(6) - 0.1795161243106645437322408d0 *&
 c(7) ,& 
 0.8671038084174692625075159d-2 *&
 c(4) +&
 0.1736073411355428563685818d0 *&
 c(6) +&
 0.5331362125287625412555844d-1&
*&
 c(8) - 0.2424935262404526301801157d0 *&
 c(5) +&
 0.1569015257678588270609004d0 *&
 c(7) ,& 
 -0.8607044252686413302647675d-3 *&
 c(4) - 0.1748074708673904989293256d0 *&
 c(5) - 0.3132544850115050165022338d0 *&
 c(8) - 0.2500000000000000000000000d-1 *&
 c(9) - 0.3169166305310429271303167d0 *&
 c(7) - 0.6691607091647929161078591d0 *&
 c(6) ,& 
 0.2239223735771599178951297d-3 *&
 c(4) +&
 0.1275437785430956673825710d0 *&
 c(5) +&
 0.1011699483929608164601067d1 *&
 c(6) +&
 0.9698817275172575247533506d0&
*&
 c(8) +&
 0.1250000000000000000000000d0 *&
 c(9) +&
 0.5555555555555555555555556d-2 *&
 c(10) +&
 0.4823177543031281500117826d0 *&
 c(7) ,& 
 -0.3784113973033012949863031d-1 *&
 c(5) - 0.2997556885134827361576001d0 *&
 c(6) - 0.3000000000000000000000000d0 *&
 c(9) - 0.2500000000000000000000000d-1 *&
 c(10) - 0.3991486867446821178415359d0 *&
 c(7) - 0.4382544850115050165022338d0 *&
 c(8) ,& 
 0.4698146218022683933926520d-1 *&
 c(6) - 0.2966863787471237458744416d0 *&
 c(8) +&
 0.5000000000000000000000000d-1&
*&
 c(10) +&
 0.1716355704146006481727960d0 *&
 c(7) +&
 0.3069346152296258362380356d-2 *&
 c(5) +&
 0.1750000000000000000000000d0 *&
 c(9)/)
 R(1+7,1:9) = (/ 0.1559652871136785763960685d-1 *&
 c(5) - 0.6486184157331537899459796d-2 *&
 c(6) - 0.9110344554036319740147054d-2 *&
 c(7) ,& 
 -0.8368141434403455353724691d-1 *&
 c(5) +&
 0.4093499466767054661591066d-1 *&
 c(6) +&
 0.4274641967636400692133626d-1 *&
 c(7) ,& 
 0.1824328174134289562208038d0 *&
 c(5) - 0.1059816030196818445908057d0 *&
 c(6) - 0.7645121439374711162999809d-1&
*&
 c(7) ,& 
 -0.2022909060111751565150958d0 *&
 c(5) +&
 0.1453421858063658498587377d0 *&
 c(6) +&
 0.5694872020480930665635812d-1 *&
 c(7) ,& 
 -0.8631683980217122275970376d-1 *&
 c(6) +&
 0.2698842360470999243492629d-1 *&
 c(7) +&
 0.8098194147715651085292754d-1 *&
 c(5) - 0.3276463639080639163926118d-1 *&
 c(8) ,& 
 0.3354661791693352108660900d-1 *&
 c(5) - 0.3343620022386971405018586d0 *&
 c(7) +&
 0.5000000000000000000000000d-1 *&
 c(9) +&
 0.2169790609807602750804271d0 *&
 c(6) +&
 0.1838363233410033443348225d0&
*&
 c(8) ,& 
 -0.3784113973033012949863031d-1 *&
 c(5) - 0.2997556885134827361576001d0 *&
 c(6) - 0.3000000000000000000000000d0 *&
 c(9) - 0.2500000000000000000000000d-1 *&
 c(10) - 0.3991486867446821178415359d0 *&
 c(7) - 0.4382544850115050165022338d0 *&
 c(8) ,& 
 0.1230328942716804455358698d-1 *&
 c(5) +&
 0.1183647529645898332481833d0 *&
 c(6) +&
 0.9410511898227943334189628d0 *&
 c(7) +&
 0.9500000000000000000000000d0 *&
 c(9) +&
 0.1250000000000000000000000d0 *&
 c(10) +&
 0.5555555555555555555555556d-2 *&
 c(11)&
+&
 0.5699474344521144554459336d0 *&
 c(8) ,& 
 -0.2308067892671916339568942d-1 *&
 c(6) - 0.2986625053775149497180439d0 *&
 c(7) - 0.3000000000000000000000000d0 *&
 c(10) - 0.2500000000000000000000000d-1 *&
 c(11) - 0.1047734860515050802561078d-2 *&
 c(5) - 0.4272090808352508360837056d0 *&
 c(8) - 0.4250000000000000000000000d0 *&
 c(9)/)
 R(1+8,1:9) = (/ 0.5593983696629863059347067d-3 *&
 c(6) - 0.1384822535100796372263822d-2 *&
 c(5) +&
 0.8254241654378100663291154d-3 *&
 c(7) ,& 
 -0.3576545132696983143406173d-2 *&
 c(6) +&
 0.7389399124121078682094445d-2&
*&
 c(5) - 0.3812853991424095538688273d-2 *&
 c(7) ,& 
 0.9209089963443799485648361d-2 *&
 c(6) - 0.1591502818872493167091475d-1 *&
 c(5) +&
 0.6705938225281132185266388d-2 *&
 c(7) ,& 
 -0.1200429618441003833696998d-1 *&
 c(6) - 0.4776915669385923841535432d-2 *&
 c(7) +&
 0.1678121185379596217850541d-1 *&
 c(5) ,& 
 0.7462059484530855073291365d-2 *&
 c(6) - 0.8121640361668678949573496d-3 *&
 c(7) +&
 0.5522702088127090209264064d-3 *&
 c(8) - 0.7202165657176696199260422d-2 *&
 c(5) ,& 
 0.2912518476823004642951502d-1 *&
 c(7) +&
 0.2279091916474916391629437d-1&
*&
 c(8) - 0.3068985997518740530511593d-1 *&
 c(6) - 0.1781799513347360596249022d-2 *&
 c(5) - 0.3055555555555555555555556d-1 *&
 c(9) ,& 
 0.4698146218022683933926520d-1 *&
 c(6) - 0.2966863787471237458744416d0 *&
 c(8) +&
 0.5000000000000000000000000d-1 *&
 c(10) +&
 0.1716355704146006481727960d0 *&
 c(7) +&
 0.3069346152296258362380356d-2 *&
 c(5) +&
 0.1750000000000000000000000d0 *&
 c(9) ,& 
 -0.2308067892671916339568942d-1 *&
 c(6) - 0.2986625053775149497180439d0 *&
 c(7) - 0.3000000000000000000000000d0 *&
 c(10) - 0.2500000000000000000000000d-1 *&
 c(11)&
- 0.1047734860515050802561078d-2 *&
 c(5) - 0.4272090808352508360837056d0 *&
 c(8) - 0.4250000000000000000000000d0 *&
 c(9) ,& 
 0.5139370221149109977041877d-2 *&
 c(6) +&
 0.1247723215009422001393184d0 *&
 c(7) +&
 0.9505522702088127090209264d0 *&
 c(8) +&
 0.9500000000000000000000000d0 *&
 c(10) +&
 0.1250000000000000000000000d0 *&
 c(11) +&
 0.5555555555555555555555556d-2 *&
 c(12) +&
 0.9159362465153641826887659d-4 *&
 c(5) +&
 0.5611111111111111111111111d0 *&
 c(9)/)

R(m-8,m-8:m)=(/0.5555555555555555555555556d-2 *&
 c(m-11) +&
 0.1250000000000000000000000d0&
*&
 c(m-10) +&
 0.9500000000000000000000000d0 *&
 c(m-9) +&
 0.9505522702088127090209264d0 *&
c(m-7) +&
 0.1247205076844361998744053d0 *&
 c(m-6) +&
 0.5139370221149109977041877d-2 *&
 c(m-5)&
+&
 0.5611111111111111111111111d0 *&
 c(m-8) +&
 0.1434074411575366831819799d-3 *&
 c(m-4) ,& 
 -0.2500000000000000000000000d-1&
*&
 c(m-10) - 0.3000000000000000000000000d0 *&
 c(m-9) - 0.2980649679116425253322056d0 *&
 c(m-6) - 0.2308067892671916339568942d-1&
*&
 c(m-5) - 0.4250000000000000000000000d0 *&
 c(m-8) - 0.4272090808352508360837056d0 *&
 c(m-7) - 0.1645272326387475188399322d-2 *&
c(m-4) ,& 
 0.5000000000000000000000000d-1 *&
 c(m-9) - 0.2966863787471237458744416d0 *&
 c(m-7) +&
 0.4698146218022683933926520d-1 *&
c(m-5) +&
 0.1750000000000000000000000d0 *&
 c(m-8) +&
 0.1700291833903489463825077d0 *&
 c(m-6) +&
 0.4675733176547960152668626d-2 *&
c(m-4) ,& 
 0.2279091916474916391629437d-1 *&
 c(m-7) +&
 0.3097763128598982561225538d-1 *&
 c(m-6) - 0.3055555555555555555555556d-1 *&
c(m-8) - 0.3068985997518740530511593d-1 *&
 c(m-5) - 0.3634246031107139778989373d-2 *&
 c(m-4) ,& 
 0.5522702088127090209264064d-3 *&
 c(m-7)&
- 0.3265435411305071914756373d-2 *&
 c(m-6) +&
 0.7462059484530855073291365d-2 *&
 c(m-5) - 0.4748894282038492179461399d-2 *&
 c(m-4) ,& 
 0.6272075574042975468177820d-3&
*&
 c(m-6) - 0.1200429618441003833696998d-1 *&
 c(m-5) +&
 0.1137708862700574079015220d-1 *&
 c(m-4) ,& 
 0.9209089963443799485648361d-2 *&
 c(m-5) - 0.3129629392354775191148163d-3&
*&
 c(m-6) - 0.8896127024208321966533544d-2 *&
 c(m-4) ,& 
 -0.3576545132696983143406173d-2 *&
 c(m-5) +&
 0.4335019854436220306755673d-3 *&
 c(m-6) +&
 0.3143043147253361112730605d-2&
*&
 c(m-4) ,& 
 0.5593983696629863059347067d-3 *&
 c(m-5) - 0.1446656414398166805849327d-3 *&
 c(m-6) - 0.4147327282231696253497740d-3 *&
 c(m-4)/)
 R(m-8+1,m-8:m) = (/ -0.2500000000000000000000000d-1 *&
c(m-10) - 0.3000000000000000000000000d0 *&
 c(m-9) - 0.2980649679116425253322056d0 *&
 c(m-6) - 0.2308067892671916339568942d-1 *&
 c(m-5) - 0.4250000000000000000000000d0 *&
 c(m-8) - 0.4272090808352508360837056d0&
*&
 c(m-7) - 0.1645272326387475188399322d-2 *&
 c(m-4) ,& 
 0.5555555555555555555555556d-2 *&
 c(m-10) +&
 0.1250000000000000000000000d0 *&
 c(m-9) +&
 0.9500000000000000000000000d0 *&
 c(m-8) +&
 0.9341601509609901526962449d0&
*&
 c(m-6) +&
 0.1183647529645898332481833d0 *&
 c(m-5) +&
 0.1919432828897222527630486d-1 *&
 c(m-4) +&
 0.5699474344521144554459336d0 *&
 c(m-7) ,& 
 -0.2500000000000000000000000d-1 *&
 c(m-9) - 0.3000000000000000000000000d0&
*&
 c(m-8) - 0.2997556885134827361576001d0 *&
 c(m-5) - 0.5636663150858098975790317d-1 *&
 c(m-4) - 0.4382544850115050165022338d0 *&
 c(m-7) - 0.3806231949664312575822630d0 *&
 c(m-6) ,& 
 0.5000000000000000000000000d-1 *&
 c(m-8)&
- 0.3557251496099816106154206d0 *&
 c(m-6) +&
 0.5490976528821799120017102d-1 *&
 c(m-4) +&
 0.1838363233410033443348225d0 *&
 c(m-7) +&
 0.2169790609807602750804271d0 *&
 c(m-5) ,& 
 0.5528052133944605740009217d-1 *&
 c(m-6) - 0.8631683980217122275970376d-1&
*&
 c(m-5) - 0.3276463639080639163926118d-1 *&
 c(m-7) +&
 0.5268984374242044588776166d-1 *&
 c(m-4) ,& 
 -0.5373770512016897565958305d-2 *&
 c(m-6) +&
 0.1453421858063658498587377d0 *&
 c(m-5) - 0.1399684152943489522927794d0 *&
 c(m-4) ,& 
 -0.1059816030196818445908057d0&
*&
 c(m-5) +&
 0.1014880675788250237247178d0 *&
 c(m-4) +&
 0.4493535440856820866087846d-2 *&
 c(m-6) ,& 
 0.4093499466767054661591066d-1 *&
 c(m-5) - 0.3471075437892810033585296d-1 *&
 c(m-4) - 0.6224240288742446280057699d-2 *&
 c(m-6) ,& 
 -0.6486184157331537899459796d-2&
*&
 c(m-5) +&
 0.4409068609809831485979484d-2 *&
 c(m-4) +&
 0.2077115547521706413480312d-2 *&
 c(m-6)/)
 R(m-8+2,m-8:m) = (/ 0.5000000000000000000000000d-1 *&
 c(m-9) - 0.2966863787471237458744416d0 *&
 c(m-7) +&
 0.4698146218022683933926520d-1 *&
 c(m-5) +&
 0.1750000000000000000000000d0 *&
c(m-8) +&
 0.1700291833903489463825077d0 *&
 c(m-6) +&
 0.4675733176547960152668626d-2 *&
 c(m-4) ,& 
 -0.2500000000000000000000000d-1 *&
 c(m-9) - 0.3000000000000000000000000d0 *&
 c(m-8) - 0.2997556885134827361576001d0 *&
 c(m-5) - 0.5636663150858098975790317d-1 *&
 c(m-4)&
- 0.4382544850115050165022338d0 *&
 c(m-7) - 0.3806231949664312575822630d0 *&
 c(m-6) ,& 
 0.5555555555555555555555556d-2 *&
 c(m-9) +&
 0.1250000000000000000000000d0 *&
 c(m-8) +&
 0.9698817275172575247533506d0 *&
 c(m-7) +&
 0.1011699483929608164601067d1 *&
 c(m-5) +&
 0.1773466968705924819112984d0&
*&
 c(m-4) +&
 0.2239223735771599178951297d-3 *&
 c(m-3) +&
 0.4325148359756313354830552d0 *&
 c(m-6) ,& 
 -0.2500000000000000000000000d-1 *&
 c(m-8) - 0.3132544850115050165022338d0 *&
 c(m-7) - 0.2322389872063761557916742d0 *&
 c(m-4) - 0.8607044252686413302647675d-3 *&
 c(m-3) - 0.2594851141920572702679681d0&
*&
 c(m-6) - 0.6691607091647929161078591d0 *&
 c(m-5) ,& 
 0.5331362125287625412555844d-1 *&
 c(m-7) +&
 0.1736073411355428563685818d0 *&
 c(m-5) +&
 0.8671038084174692625075159d-2 *&
 c(m-3) +&
 0.8084259844422177692569663d-1 *&
 c(m-6) - 0.1664345989168155800449120d0 *&
 c(m-4) ,& 
 -0.5013247356072127938999311d0&
*&
 c(m-5) +&
 0.5021853752328231128475915d0 *&
 c(m-4) - 0.1197175073672143005877150d-1 *&
 c(m-6) ,& 
 0.3328335104489738933610597d0 *&
 c(m-5) - 0.3179803804558436283847901d0 *&
 c(m-4) - 0.5111353189352474549563559d-2 *&
 c(m-3) - 0.9741776803777790426705996d-2 *&
 c(m-6) ,& 
 -0.1080774142196007991746827d0 *&
c(m-5) +&
 0.9941834083648937298100811d-1 *&
 c(m-4) - 0.4834791406446907590553793d-2 *&
 c(m-3) +&
 0.1349386478955833378422842d-1 *&
 c(m-6) ,& 
 0.1319674981073749167009902d-1 *&
 c(m-5) - 0.1060554802883657391328704d-1 *&
 c(m-4) +&
 0.1911888563316170927411831d-2 *&
 c(m-3) - 0.4503090345217088684223814d-2 *&
c(m-6)/)
 R(m-8+3,m-8:m) = (/ 0.2279091916474916391629437d-1 *&
 c(m-7) +&
 0.3097763128598982561225538d-1 *&
 c(m-6) - 0.3055555555555555555555556d-1 *&
 c(m-8) - 0.3068985997518740530511593d-1 *&
 c(m-5) - 0.3634246031107139778989373d-2 *&
 c(m-4) ,& 
 0.5000000000000000000000000d-1 *&
 c(m-8) - 0.3557251496099816106154206d0 *&
 c(m-6)&
+&
 0.5490976528821799120017102d-1 *&
 c(m-4) +&
 0.1838363233410033443348225d0 *&
 c(m-7) +&
 0.2169790609807602750804271d0 *&
 c(m-5) ,& 
 -0.2500000000000000000000000d-1 *&
 c(m-8) - 0.3132544850115050165022338d0 *&
 c(m-7) - 0.2322389872063761557916742d0 *&
 c(m-4) - 0.8607044252686413302647675d-3 *&
 c(m-3) - 0.2594851141920572702679681d0&
*&
 c(m-6) - 0.6691607091647929161078591d0 *&
 c(m-5) ,& 
 0.5555555555555555555555556d-2 *&
 c(m-8) +&
 0.1338363233410033443348225d0 *&
 c(m-7) +&
 0.7391887916719206077121040d0 *&
 c(m-6) +&
 0.6490333320052011212240632d0 *&
 c(m-4) +&
 0.3308343404200968256656458d-2 *&
 c(m-3) +&
 0.2062575706647430620228133d-3 *&
 c(m-2) +&
 0.3088241944378964404772302d-3&
*&
 c(m-1) +&
 0.4227226173449345042468960d-3 *&
 c(m) +&
 0.1190362071861893051132274d1 *&
 c(m-5) ,& 
 -0.2720908083525083608370563d-1 *&
 c(m-7) - 0.1931148612480615118957263d0 *&
 c(m-6) - 0.3332941113251635390801278d-1 *&
 c(m-3) +&
 0.2605582646183255957264249d-3 *&
 c(m-2) +&
 0.1432116665752147607469646d-2 *&
 c(m-1) +&
 0.1207544072304193806052558d-2 *&
c(m) - 0.1348436986667115543203552d1 *&
 c(m-5) +&
 0.1687723507780044227927853d-1 *&
 c(m-4) ,& 
 0.3590669644811151307464697d-1 *&
 c(m-6) - 0.5925443480724830632401754d0 *&
 c(m-4) - 0.6897142765790609546343709d-2 *&
 c(m-2) - 0.5363098747528542488971874d-2 *&
 c(m-1) - 0.5205147429855955657625694d-2 *&
 c(m) +&
 0.9977064356292750529201981d0 *&
 c(m-5) ,& 
 0.7272438906214475928744770d-1&
*&
 c(m-4) +&
 0.1964682777744275219350831d-1 *&
 c(m-3) - 0.5952475275883259619711594d-2 *&
 c(m-1) - 0.1635430866921887819487473d-2 *&
 c(m) +&
 0.2921234010758621482958052d-1 *&
 c(m-6) - 0.4659516693228870973898560d0 *&
 c(m-5) ,& 
 0.5891947149681041048896399d-1 *&
 c(m-4) +&
 0.1858378996391679448655070d-1 *&
 c(m-3) +&
 0.7240905383565181316381731d-2 *&
 c(m-2) +&
 0.2349927974590068869356781d-1&
*&
 c(m) - 0.4046360079256766884300687d-1 *&
 c(m-6) +&
 0.1223513270418807666970488d0 *&
 c(m-5) ,& 
 -0.2404661162020836566908542d-1 *&
 c(m-4) - 0.7348845587775519698437916d-2 *&
 c(m-3) - 0.8105784530576404277872603d-3 *&
 c(m-2) +&
 0.9574633163221758060736592d-2 *&
 c(m-1) - 0.1828896813877197352675410d-1 *&
 c(m) +&
 0.1350326632905990039353503d-1 *&
 c(m-6) - 0.1315967038382618382356495d-1 *&
 c(m-5)/)
 R(m-8+4,m-8:m) = (/ &
0.5522702088127090209264064d-3 *&
 c(m-7) - 0.3265435411305071914756373d-2 *&
 c(m-6) +&
 0.7462059484530855073291365d-2 *&
 c(m-5) - 0.4748894282038492179461399d-2 *&
 c(m-4) ,& 
 0.5528052133944605740009217d-1 *&
 c(m-6) - 0.8631683980217122275970376d-1 *&
 c(m-5) - 0.3276463639080639163926118d-1 *&
 c(m-7) +&
 0.5268984374242044588776166d-1 *&
 c(m-4) ,& 
 0.5331362125287625412555844d-1 *&
 c(m-7) +&
 0.1736073411355428563685818d0&
*&
 c(m-5) +&
 0.8671038084174692625075159d-2 *&
 c(m-3) +&
 0.8084259844422177692569663d-1 *&
 c(m-6) - 0.1664345989168155800449120d0 *&
 c(m-4) ,& 
 -0.2720908083525083608370563d-1 *&
 c(m-7) - 0.1931148612480615118957263d0 *&
 c(m-6) - 0.3332941113251635390801278d-1 *&
 c(m-3) +&
 0.2605582646183255957264249d-3 *&
 c(m-2) +&
 0.1432116665752147607469646d-2 *&
 c(m-1) +&
 0.1207544072304193806052558d-2 *&
 c(m) - 0.1348436986667115543203552d1&
*&
 c(m-5) +&
 0.1687723507780044227927853d-1 *&
 c(m-4) ,& 
 0.6107825764368264576481962d-2 *&
 c(m-7) +&
 0.1155752633643216628010304d0 *&
 c(m-6) +&
 0.2096413329579026439044119d1 *&
 c(m-5) +&
 0.3357721707576477199985656d0 *&
 c(m-3) +&
 0.3291545083271862858501887d-3 *&
 c(m-2) +&
 0.6641183499427826101618457d-2 *&
 c(m-1) +&
 0.3449455095910233625229891d-2 *&
 c(m) +&
 0.8270696421223286922584620d0 *&
 c(m-4) ,& 
 -0.4995827370863505253765970d-1&
*&
 c(m-6) - 0.1263507837371824205693950d1 *&
 c(m-5) - 0.8712928907711754187084757d-2 *&
 c(m-2) - 0.2487040599390160764166412d-1 *&
 c(m-1) - 0.1486895819265604128572498d-1 *&
 c(m) - 0.1726565392121567634950213d1 *&
 c(m-4) ,& 
 0.5402985338373433052255418d0 *&
 c(m-5) - 0.1979290298620869974478871d0 *&
 c(m-3) - 0.2760353365637712827793337d-1 *&
 c(m-1) - 0.4671751091575462868310238d-2 *&
 c(m) - 0.6952587985456154591014641d-1 *&
 c(m-6) +&
 0.1571507277911208446562686d1&
*&
 c(m-4) ,& 
 -0.1319358558853174530078498d0 *&
 c(m-5) - 0.1872196143003808021730728d0 *&
 c(m-3) +&
 0.9147192682075630179962131d-2 *&
 c(m-2) +&
 0.6712774475803763988977355d-1 *&
 c(m) +&
 0.9630407686703666967100804d-1 *&
 c(m-6) - 0.6882132817757726722141421d0 *&
 c(m-4) ,& 
 0.1241625568998496895352046d-1 *&
 c(m-5) +&
 0.7403484645316174090533193d-1 *&
 c(m-3) - 0.1023976547309387874453988d-2 *&
 c(m-2) +&
 0.4440063948509876221050939d-1 *&
 c(m-1) - 0.5224403464202056316702078d-1&
*&
 c(m) - 0.3213800979246298453953842d-1 *&
 c(m-6) +&
 0.1178181682424363524005403d0 *&
 c(m-4)/)
 R(m-8+5,m-8:m) = (/ 0.6272075574042975468177820d-3 *&
 c(m-6) - 0.1200429618441003833696998d-1 *&
 c(m-5) +&
 0.1137708862700574079015220d-1 *&
 c(m-4) ,& 
 -0.5373770512016897565958305d-2 *&
 c(m-6) +&
 0.1453421858063658498587377d0 *&
 c(m-5) - 0.1399684152943489522927794d0 *&
 c(m-4) ,& 
 -0.5013247356072127938999311d0 *&
 c(m-5) +&
 0.5021853752328231128475915d0 *&
 c(m-4) - 0.1197175073672143005877150d-1 *&
c(m-6) ,& 
 0.3590669644811151307464697d-1 *&
 c(m-6) - 0.5925443480724830632401754d0 *&
 c(m-4) - 0.6897142765790609546343709d-2 *&
 c(m-2) - 0.5363098747528542488971874d-2 *&
 c(m-1) - 0.5205147429855955657625694d-2 *&
 c(m) +&
 0.9977064356292750529201981d0 *&
 c(m-5) ,& 
 -0.4995827370863505253765970d-1 *&
 c(m-6) - 0.1263507837371824205693950d1 *&
 c(m-5) - 0.8712928907711754187084757d-2 *&
 c(m-2) - 0.2487040599390160764166412d-1 *&
 c(m-1) - 0.1486895819265604128572498d-1 *&
 c(m) - 0.1726565392121567634950213d1&
*&
 c(m-4) ,& 
 0.2760393423824887721078848d-1 *&
 c(m-6) +&
 0.1190550338687608873798462d1 *&
 c(m-5) +&
 0.4253084328734353394994388d1 *&
 c(m-4) +&
 0.2306367624634749229113646d0 *&
 c(m-2) +&
 0.9313657638804699948929701d-1 *&
 c(m-1) +&
 0.6409299775987186986730499d-1 *&
 c(m) ,& 
 -0.8755807343482262259774782d0 *&
 c(m-5) - 0.3645285178085761821545207d1 *&
 c(m-4) +&
 0.1033717994630886401730470d0 *&
 c(m-1) +&
 0.2013769413884797246646959d-1 *&
 c(m) +&
 0.4106783858513785463625543d-1 *&
 c(m-6) ,& 
 0.3956598149904136332753521d0&
*&
 c(m-5) +&
 0.1630560443616104907615866d1 *&
 c(m-4) - 0.2421320004064592721552708d0 *&
 c(m-2) - 0.2893557395653431666593814d0 *&
 c(m) - 0.5688529641249387985434413d-1 *&
 c(m-6) ,& 
 -0.7684117160199014594442072d-1 *&
 c(m-5) - 0.2928439026361256842196229d0 *&
 c(m-4) +&
 0.2710530961648671297733465d-1 *&
 c(m-2) - 0.1662748711097054895317080d0 *&
 c(m-1) +&
 0.2251991532891353212689574d0 *&
 c(m) +&
 0.1898341454096471754822498d-1 *&
 c(m-6)/)
 R(m-8+6,m-8:m) = (/ 0.9209089963443799485648361d-2 *&
 c(m-5) - 0.3129629392354775191148163d-3 *&
 c(m-6) -&
0.8896127024208321966533544d-2 *&
 c(m-4) ,& 
 -0.1059816030196818445908057d0 *&
 c(m-5) +&
 0.1014880675788250237247178d0 *&
 c(m-4) +&
 0.4493535440856820866087846d-2 *&
 c(m-6) ,& 
 0.3328335104489738933610597d0 *&
 c(m-5) - 0.3179803804558436283847901d0 *&
 c(m-4) - 0.5111353189352474549563559d-2 *&
 c(m-3) - 0.9741776803777790426705996d-2 *&
 c(m-6) ,& 
 0.7272438906214475928744770d-1 *&
 c(m-4) +&
 0.1964682777744275219350831d-1 *&
 c(m-3) - 0.5952475275883259619711594d-2 *&
 c(m-1) - 0.1635430866921887819487473d-2 *&
 c(m) +&
 0.2921234010758621482958052d-1&
*&
 c(m-6) - 0.4659516693228870973898560d0 *&
 c(m-5) ,& 
 0.5402985338373433052255418d0 *&
 c(m-5) - 0.1979290298620869974478871d0 *&
 c(m-3) - 0.2760353365637712827793337d-1 *&
 c(m-1) - 0.4671751091575462868310238d-2 *&
 c(m) - 0.6952587985456154591014641d-1 *&
 c(m-6) +&
 0.1571507277911208446562686d1 *&
 c(m-4) ,& 
 -0.8755807343482262259774782d0 *&
 c(m-5) - 0.3645285178085761821545207d1 *&
 c(m-4) +&
 0.1033717994630886401730470d0 *&
 c(m-1) +&
 0.2013769413884797246646959d-1 *&
 c(m) +&
 0.4106783858513785463625543d-1 *&
 c(m-6) ,& 
 0.1070920689960817104203947d1&
*&
 c(m-5) +&
 0.3717418466925056542408153d1 *&
 c(m-4) +&
 0.1166740554279680007487795d0 *&
 c(m-3) +&
 0.1147318200715868527529827d0 *&
 c(m-1) +&
 0.6327161147136873807796515d-2 *&
 c(m) +&
 0.6235373239336055200426697d-1 *&
 c(m-6) ,& 
 -0.6639605248735044787146222d0 *&
 c(m-5) - 0.1865625445986772763641423d1 *&
 c(m-4) +&
 0.1103611313171476425250639d0 *&
 c(m-3) - 0.9091410269992464604926176d-1 *&
 c(m) - 0.8636954541126674177407762d-1 *&
 c(m-6) ,& 
 0.1582127073537215443965653d0 *&
 c(m-5) +&
 0.3746489300753517635549495d0 *&
 c(m-4) - 0.4364163147111892346990101d-1&
*&
 c(m-3) - 0.1845476106024151050283847d0 *&
 c(m-1) +&
 0.7075642937243715046279337d-1 *&
 c(m) +&
 0.2882271848190011329385407d-1 *&
 c(m-6)/)
 R(m-8+7,m-8:m) = (/ -0.3576545132696983143406173d-2 *&
 c(m-5) +&
 0.4335019854436220306755673d-3 *&
 c(m-6) +&
 0.3143043147253361112730605d-2 *&
 c(m-4) ,& 
 0.4093499466767054661591066d-1 *&
 c(m-5) - 0.3471075437892810033585296d-1 *&
 c(m-4) - 0.6224240288742446280057699d-2 *&
 c(m-6) ,& 
 -0.1080774142196007991746827d0 *&
 c(m-5) +&
 0.9941834083648937298100811d-1 *&
 c(m-4) - 0.4834791406446907590553793d-2 *&
 c(m-3) +&
 0.1349386478955833378422842d-1&
*&
 c(m-6) ,& 
 0.5891947149681041048896399d-1 *&
 c(m-4) +&
 0.1858378996391679448655070d-1 *&
 c(m-3) +&
 0.7240905383565181316381731d-2 *&
 c(m-2) +&
 0.2349927974590068869356781d-1 *&
 c(m) - 0.4046360079256766884300687d-1 *&
 c(m-6) +&
 0.1223513270418807666970488d0 *&
 c(m-5) ,& 
 -0.1319358558853174530078498d0 *&
 c(m-5) - 0.1872196143003808021730728d0 *&
 c(m-3) +&
 0.9147192682075630179962131d-2 *&
 c(m-2) +&
 0.6712774475803763988977355d-1 *&
 c(m) +&
 0.9630407686703666967100804d-1 *&
 c(m-6) - 0.6882132817757726722141421d0 *&
 c(m-4) ,& 
 0.3956598149904136332753521d0 *&
c(m-5) +&
 0.1630560443616104907615866d1 *&
 c(m-4) - 0.2421320004064592721552708d0 *&
 c(m-2) - 0.2893557395653431666593814d0 *&
 c(m) - 0.5688529641249387985434413d-1 *&
 c(m-6) ,& 
 -0.6639605248735044787146222d0 *&
 c(m-5) - 0.1865625445986772763641423d1 *&
 c(m-4) +&
 0.1103611313171476425250639d0 *&
 c(m-3) - 0.9091410269992464604926176d-1 *&
 c(m) - 0.8636954541126674177407762d-1 *&
 c(m-6) ,& 
 0.4681819359722749441073885d0 *&
 c(m-5) +&
 0.1015239189167790053447110d1 *&
 c(m-4) +&
 0.1043897828092562609502636d0 *&
 c(m-3) +&
 0.2542001760457345743492403d0 *&
 c(m-2) +&
 0.1306332157111667628555907d1&
*&
 c(m) +&
 0.1196351539550049336518187d0 *&
 c(m-6) ,& 
 -0.1195777325611201766551392d0 *&
 c(m-5) - 0.2187310061229745694542609d0 *&
 c(m-4) - 0.4128029838349298819825156d-1 *&
 c(m-3) - 0.2845627370491611369031341d-1 *&
 c(m-2) - 0.1016689339350338144430605d1 *&
 c(m) - 0.3992391469197282238624438d-1 *&
 c(m-6)/)
 R(m-8+8,m-8:m) = (/ 0.5593983696629863059347067d-3 *&
 c(m-5) - 0.1446656414398166805849327d-3 *&
 c(m-6) - 0.4147327282231696253497740d-3 *&
 c(m-4) ,& 
 -0.6486184157331537899459796d-2 *&
 c(m-5) +&
 0.4409068609809831485979484d-2 *&
 c(m-4) +&
 0.2077115547521706413480312d-2 *&
 c(m-6) ,& 
 0.1319674981073749167009902d-1&
*&
 c(m-5) - 0.1060554802883657391328704d-1 *&
 c(m-4) +&
 0.1911888563316170927411831d-2 *&
 c(m-3) - 0.4503090345217088684223814d-2 *&
 c(m-6) ,& 
 -0.2404661162020836566908542d-1 *&
 c(m-4) - 0.7348845587775519698437916d-2 *&
 c(m-3) - 0.8105784530576404277872603d-3 *&
 c(m-2) +&
 0.9574633163221758060736592d-2 *&
 c(m-1) - 0.1828896813877197352675410d-1 *&
 c(m) +&
 0.1350326632905990039353503d-1 *&
 c(m-6) - 0.1315967038382618382356495d-1 *&
 c(m-5) ,& 
 0.1241625568998496895352046d-1 *&
 c(m-5) +&
 0.7403484645316174090533193d-1 *&
 c(m-3) - 0.1023976547309387874453988d-2 *&
 c(m-2) +&
 0.4440063948509876221050939d-1&
*&
 c(m-1) - 0.5224403464202056316702078d-1 *&
 c(m) - 0.3213800979246298453953842d-1 *&
 c(m-6) +&
 0.1178181682424363524005403d0 *&
 c(m-4) ,& 
 -0.7684117160199014594442072d-1 *&
 c(m-5) - 0.2928439026361256842196229d0 *&
 c(m-4) +&
 0.2710530961648671297733465d-1 *&
 c(m-2) - 0.1662748711097054895317080d0 *&
 c(m-1) +&
 0.2251991532891353212689574d0 *&
 c(m) +&
 0.1898341454096471754822498d-1 *&
 c(m-6) ,& 
 0.1582127073537215443965653d0 *&
 c(m-5) +&
 0.3746489300753517635549495d0 *&
 c(m-4) - 0.4364163147111892346990101d-1 *&
 c(m-3) - 0.1845476106024151050283847d0 *&
 c(m-1) +&
 0.7075642937243715046279337d-1 *&
 c(m)&
+&
 0.2882271848190011329385407d-1 *&
 c(m-6) ,& 
 -0.1195777325611201766551392d0 *&
 c(m-5) - 0.2187310061229745694542609d0 *&
 c(m-4) - 0.4128029838349298819825156d-1 *&
 c(m-3) - 0.2845627370491611369031341d-1 *&
 c(m-2) - 0.1016689339350338144430605d1 *&
 c(m) - 0.3992391469197282238624438d-1 *&
 c(m-6) ,& 
 0.3167964748016105299646518d-1 *&
 c(m-5) +&
 0.4976563420877041544013670d-1 *&
 c(m-4) +&
 0.1632404042590951953384672d-1 *&
 c(m-3) +&
 0.3185519088796429015220016d-2 *&
 c(m-2) +&
 0.2968472090638000742888467d0 *&
 c(m-1) +&
 0.7912667594695582093926295d0 *&
 c(m) +&
 0.1332316557164627464149716d-1 *&
 c(m-6)/)
  end subroutine
end module hydrofrac
