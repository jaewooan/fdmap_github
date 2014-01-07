module hydrofrac

  use fd, only : limits

  implicit none

  ! USE_HF = flag indicating hydraulic fracture
  ! BNDM,BNDP = process is resposible for minus/plus boundary 
  ! DIRECTION = direction of fracture
  ! COUPLED = full fluid-solid coupling (F=neglect crack opening in mass balance)
  ! RHO0 = (initial) density
  ! K0 = bulk modulus
  ! C0 = sound speed
  ! A,B = position of lower and upper walls (b-a=width)
  ! A0,B0 = initial values of a,b
  ! U = velocity
  ! P = pressure
  ! DU = rate of change of velocity
  ! DP = rate of change of pressure
  ! DA,DB = rate of wall position
  ! SATM,SATP = SAT penalty weights (~wave speed / grid spacing)
  ! L = limits type (indices for FD operations)
  ! H = grid spacing (*HARD-CODED NOW, LATER UPDATE FOR COMPABILITY WITH GRID*)

  type :: hf_type
     logical :: use_HF,coupled,bndm,bndp
     character(1) :: direction
     real :: rho0,K0,c0,h,SATm,SATp
     real,dimension(:),allocatable :: a,b,a0,b0,u,p,Du,Dp,Da,Db
     type(limits) :: L
  end type hf_type

contains


  subroutine init_hydrofrac(iface,HF,m,p,x,y,input,echo,skip,direction,mg,pg, &
       process_m,process_p,comm_m,comm_p,array)

    use fd_coeff, only : H00i
    use mpi_routines, only : is_master
    use io, only : error,write_matlab,seek_to_string  

    implicit none

    integer,intent(in) :: iface,m,p,input,echo,mg,pg,comm_m,comm_p,array
    type(hf_type),intent(out) :: HF
    real,dimension(m:p),intent(in) :: x,y
    logical,intent(in) :: skip,process_m,process_p
    character(1),intent(in) :: direction

    integer :: stat
    character(256) :: HFstr
    character(256) :: str

    logical :: hydrofrac_file,coupled
    real :: rho0,K0,w0,h
    character(256) :: filename
    integer,parameter :: nb=3 ! number of additional boundary (ghost) points needed for FD operators

    namelist /hydrofrac_list/ rho0,K0,w0,h,coupled,hydrofrac_file,filename

    ! defaults
       
    rho0 = 0d0
    K0 = 1d40
    w0 = 0d0
    h = 0d0
    coupled = .true.
    hydrofrac_file = .false.
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
    HF%h = h
    HF%coupled = coupled

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
    end if

    ! return if not needed

    if (skip) return

    ! initialize grid and fields

    HF%bndm = (m==mg)
    HF%bndp = (p==pg)
    HF%L = limits(nb,m,p,mg,pg,m-nb,p+nb,HF%bndm,HF%bndp)

    ! allocate arrays and initialize with unreasonable values to facilitate debugging

    allocate( &
         HF%a (HF%L%mb:HF%L%pb),HF%b (HF%L%mb:HF%L%pb), &
         HF%a0(HF%L%mb:HF%L%pb),HF%b0(HF%L%mb:HF%L%pb), &
         HF%u (HF%L%mb:HF%L%pb),HF%p (HF%L%mb:HF%L%pb), &
         HF%Da(HF%L%mb:HF%L%pb),HF%Db(HF%L%mb:HF%L%pb), &
         HF%Du(HF%L%mb:HF%L%pb),HF%Dp(HF%L%mb:HF%L%pb))

    HF%a0 = 1d40
    HF%b0 = 1d40
    HF%a  = 1d40
    HF%b  = 1d40
    HF%u  = 1d40
    HF%p  = 1d40
    HF%Da = 1d40
    HF%Db = 1d40
    HF%Du = 1d40
    HF%Dp = 1d40

    ! initial conditions on wall position, velocity, and pressure perturbation

    HF%a0 = -w0
    HF%b0 =  w0

    HF%u = 0d0
    HF%p = 0d0

    HF%p(m:p) = exp(-x**2)

    ! override uniform initial conditions with those from file, if needed

    if (hydrofrac_file) then
       ! both sides read file (so process may read file twice)
       if (process_m) call read_hydrofrac(HF,filename,comm_m,array)
       if (process_p) call read_hydrofrac(HF,filename,comm_p,array)
    end if

    ! current (=initial, for now) wall positions

    HF%a = HF%a0
    HF%b = HF%b0

  end subroutine init_hydrofrac


  subroutine read_hydrofrac(HF,filename,comm,array)

    use io, only : file_distributed,open_file_distributed, &
         read_file_distributed,close_file_distributed
    use mpi_routines, only : pw

    implicit none

    type(hf_type),intent(inout) :: HF
    character(*),intent(in) :: filename
    integer,intent(in) :: comm,array

    type(file_distributed) :: fh

    call open_file_distributed(fh,filename,'read',comm,array,pw)

    call read_file_distributed(fh,HF%a0(HF%L%m:HF%L%p))
    call read_file_distributed(fh,HF%b0(HF%L%m:HF%L%p))
    call read_file_distributed(fh,HF%u (HF%L%m:HF%L%p))
    call read_file_distributed(fh,HF%p (HF%L%m:HF%L%p))

    call close_file_distributed(fh)

  end subroutine read_hydrofrac


  subroutine destroy_hydrofrac(HF)

    implicit none

    type(hf_type),intent(inout) :: HF

    if (allocated(HF%a0)) deallocate(HF%a0)
    if (allocated(HF%b0)) deallocate(HF%b0)
    if (allocated(HF%a )) deallocate(HF%a )
    if (allocated(HF%b )) deallocate(HF%b )
    if (allocated(HF%u )) deallocate(HF%u )
    if (allocated(HF%p )) deallocate(HF%p )
    if (allocated(HF%Da)) deallocate(HF%Da)
    if (allocated(HF%Db)) deallocate(HF%Db)
    if (allocated(HF%Du)) deallocate(HF%Du)
    if (allocated(HF%Dp)) deallocate(HF%Dp)

  end subroutine destroy_hydrofrac


  subroutine checkpoint_hydrofrac(fh,operation,HF)

    use io, only : file_distributed, &
         read_file_distributed,write_file_distributed

    implicit none

    type(file_distributed),intent(in) :: fh
    character(*),intent(in) :: operation
    type(hf_type),intent(inout) :: HF

    integer :: m,p

    ! fields read/written to same file as iface fields,
    ! routine is only called by io_process

    if (.not.HF%use_HF) return

    m = HF%L%m
    p = HF%L%p

    select case(operation)
    case('read')
       call read_file_distributed(fh,HF%a(m:p))
       call read_file_distributed(fh,HF%b(m:p))
       call read_file_distributed(fh,HF%u(m:p))
       call read_file_distributed(fh,HF%p(m:p))
    case('write')
       call write_file_distributed(fh,HF%a(m:p))
       call write_file_distributed(fh,HF%b(m:p))
       call write_file_distributed(fh,HF%u(m:p))
       call write_file_distributed(fh,HF%p(m:p))
    end select

  end subroutine checkpoint_hydrofrac


  subroutine scale_rates_hydrofrac(HF,A)

    implicit none

    type(hf_type),intent(inout) :: HF
    real,intent(in) :: A

    if (.not.HF%use_HF) return

    HF%Da = A*HF%Da
    HF%Db = A*HF%Db
    HF%Du = A*HF%Du
    HF%Dp = A*HF%Dp

  end subroutine scale_rates_hydrofrac


  subroutine update_fields_hydrofrac(HF,dt)

    implicit none

    type(hf_type),intent(inout) :: HF
    real,intent(in) :: dt

    integer :: i

    if (.not.HF%use_HF) return

    do i = HF%L%m,HF%L%p
       HF%a(i) = HF%a(i)+dt*HF%Da(i)
       HF%b(i) = HF%b(i)+dt*HF%Db(i)
       HF%u(i) = HF%u(i)+dt*HF%Du(i)
       HF%p(i) = HF%p(i)+dt*HF%Dp(i)
    end do

  end subroutine update_fields_hydrofrac


  subroutine set_rates_hydrofrac(HF,C,comm,m,p,x,y,t)

    use mpi_routines2d, only : cartesian
    use fd, only : diff

    implicit none

    type(hf_type),intent(inout) :: HF
    type(cartesian),intent(in) :: C
    integer,intent(in) :: comm,m,p
    real,intent(in) :: x(m:p),y(m:p),t

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
       HF%Dp(i) = HF%Dp(i)-dudx(i)*HF%K0
    end do

    ! and then adding crack opening/closing term in mass balance

    if (HF%coupled) then
       do i = HF%L%m,HF%L%p
          HF%Dp(i) = HF%Dp(i)-HF%K0*(HF%Db(i)-HF%Da(i))/(HF%b0(i)-HF%a0(i))
       end do
    end if

    ! and source terms

    !do i = HF%L%m,HF%L%p
    !   HF%Dp(i) = HF%Dp(i)+0.1d0*exp(-0.5d0*(x(i)/0.1d0)**2)*exp(-0.5d0*((t-1d0)/0.1d0)**2)
    !end do

    ! add penalty terms to enforce BC with SAT method

    if (HF%bndm) then
       i = HF%L%mg
       uhat = 0d0
       phat = HF%p(i)-HF%rho0*HF%c0*HF%u(i)
       HF%Du(i) = HF%Du(i)-HF%SATm*(HF%u(i)-uhat)
       HF%Dp(i) = HF%Dp(i)-HF%SATm*(HF%p(i)-phat)
    end if

    if (HF%bndp) then
       i = HF%L%pg
       uhat = 0d0
       phat = HF%p(i)+HF%rho0*HF%c0*HF%u(i)
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
    ! evaluate these as mu*d(v-vwall)/dy for viscous fluid with fluid velocity v
    ! and wall velocity vwall *may require simultaneously solving elasticity 
    ! with wall shear stress relation for both tau and vwall

    taum = 0d0
    taup = 0d0

  end subroutine fluid_stresses


  subroutine wall_velocities(HF,i,vnm,vnp)

    implicit none

    type(HF_type),intent(inout) :: HF
    integer,intent(in) :: i
    real,intent(in) :: vnm,vnp

    HF%Da(i) = vnm
    HF%Db(i) = vnp

  end subroutine wall_velocities
  

  subroutine share_hydrofrac(HF,C)

    use mpi_routines2d, only : cartesian,populate_ghost_cells

    implicit none

    type(HF_type),intent(inout) :: HF
    type(cartesian),intent(in) :: C

    call populate_ghost_cells(C,HF%L%m,HF%L%p,HF%L%nb, &
         HF%bndm,HF%bndp,HF%u,HF%direction)
    call populate_ghost_cells(C,HF%L%m,HF%L%p,HF%L%nb, &
         HF%bndm,HF%bndp,HF%p,HF%direction)

  end subroutine share_hydrofrac


end module hydrofrac
