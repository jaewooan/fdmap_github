module erupt

  use fd, only : limits

  implicit none

  ! FLUID_TYPE = type of fluid (liquid, gas, mixture)
  ! P0 = reference pressure for eos
  ! RHO0 = density at reference pressure
  ! K0 = compressibility
  ! MU = viscosity
  ! RT0 = gas constant * Temperature
  ! exp = relevant exponent for henry's law
  ! s = solubility constant
  ! n0 = gas fraction at exsolution point
  ! TEX = exsolution time

  type :: fluid
     character(16) :: fluid_type
     real :: p0,rho0,K0,mu,RT0,exp,s,n0,Tex
  end type fluid

  ! BCND = boundary condition
  ! U0 = velocity (linearized BC)
  ! C0 = sound speed (linearized BC)
  ! M0 = Mach number at boundary (linearized BC)
  ! Z0 = wave impedance at boundary (linearized BC)
  ! U = velocity
  ! P = pressure
  ! K = SAT penalty parameter
  ! Q = conservation variables (BC value)

  type :: bc_type
     character(16) :: bcnd
     real :: u0,c0,M0,Z0,u,p,K
     real,dimension(:),allocatable :: q
  end type bc_type

  ! USE_ER = use eruption model
  ! BNDM,BNDP = process is resposible for minus/plus boundary 
  ! COUPLED = allow conduit width to evolve
  ! DIRECTION = direction of flow
  ! METHOD = numerical method for solving balance equations
  ! R = order of accuracy
  ! EQUILIBRATE = solve for w given excess pressure p-Smin
  ! MI,PI = indices of interior points
  ! ME,PE = indices of cell edges surrounding interior points
  ! FL = fluid parameters
  ! BCM,BCP = boundary conditions on minus/plus boundary
  ! G = gravitational acceleration
  ! W0 = initial conduit half-width
  ! PATM = atmospheric pressure
  ! RHO = density
  ! U = velocity
  ! P = pressure
  ! C = sound speed
  ! W = conduit half-width
  ! N = mass fraction of exsolved gas
  ! DL = grid spacing
  ! Q = conservation variables
  ! DQ = rate of change of q
  ! L = limits (indices)

  type :: er_type
     logical :: use_ER,bndm,bndp,coupled,equilibrate
     character(1) :: direction
     character(16) :: method
     integer :: r,mi,pi,me,pe
     type(fluid) :: FL
     type(bc_type) :: BCm,BCp
     real :: g,patm
     real,dimension(:),allocatable :: rho,u,p,c,w,dl,Smin,n,Dn
     real,dimension(:,:),allocatable :: q,Dq
     type(limits) :: L
  end type er_type

  ! NQ = number of conservation variables

  integer,parameter :: nq = 2

contains


  subroutine init_erupt(iface,ER,m,p,x,y,input,echo,skip,direction,mg,pg, &
       process_m,process_p,comm_m,comm_p,array,dl)

    use fd_coeff, only : nbst,init_weno
    use mpi_routines, only : is_master
    use io, only : error,write_matlab,seek_to_string  

    implicit none

    integer,intent(in) :: iface,m,p,input,echo,mg,pg,comm_m,comm_p,array
    type(er_type),intent(out) :: ER
    real,dimension(m:p),intent(in) :: x,y
    logical,intent(in) :: skip,process_m,process_p
    character(1),intent(in) :: direction
    real,dimension(m:p),intent(in),optional :: dl

    integer :: stat,i,r
    character(256) :: ERstr
    character(256) :: str
    character(16) :: method,BCm,BCp,fluid_type

    logical :: coupled,erupt_file,equilibrate
    real :: p0,rho0,K0,mu,RT0,exp,s,n0,Tex, &
         g,w0,BCm_p,BCm_u,BCp_p,BCp_u,patm,fraction,rhos
    character(256) :: filename
    integer,parameter :: nb=3

    namelist /erupt_list/ method,r,coupled,equilibrate, &
         fluid_type,p0,rho0,K0,mu,RT0,exp,s,n0,Tex, &
         g,w0,BCm,BCp,BCm_p,BCm_u,BCp_p,BCp_u,patm,fraction,rhos, &
         erupt_file,filename

    ! defaults
       
    method = 'SBP'
    r = 1
    coupled = .true.
    equilibrate = .false.
    fluid_type = 'mixture' 
    p0 = 0d0
    rho0 = 0d0
    K0 = 1d40
    mu = 0d0
    RT0 = 0d0
    exp = 0d0
    s = 0d0
    n0 = 0d0
    Tex = 0d0
    g = 0d0
    w0 = 0d0
    BCm = 'inflow-p'
    BCp = 'outflow-p'
    BCm_p = 1d40
    BCm_u = 1d40
    BCp_p = 1d40
    BCp_u = 1d40
    patm = 0d0
    fraction = 0d0
    rhos = 0d0
    erupt_file = .false.
    filename = ''

    ! read in eruption parameters

    write(str,'(a,i0,a)') '!---IFACE',iface,'---'
    call seek_to_string(input,str)
    read(input,nml=erupt_list,iostat=stat)
    if (stat>0) call error('Error in erupt_list','init_erupt')

    ER%use_ER = (rho0/=0d0)

    if (.not.ER%use_ER) return

    call init_weno

    select case(direction)
    case('x')
       ER%direction = 'y'
    case('y')
       ER%direction = 'x'
    end select

    ! input parameters

    ER%method = method
    ER%r = r
    ER%coupled = coupled
    ER%equilibrate = equilibrate
    ER%FL%fluid_type = fluid_type 
    ER%FL%p0 = p0
    ER%FL%rho0 = rho0
    ER%FL%K0 = K0
    ER%FL%mu = mu
    ER%FL%RT0 = RT0
    ER%FL%exp = exp
    ER%FL%s = s
    ER%FL%n0 = n0
    ER%FL%Tex = Tex
    ER%g = g

    ! output eruption parameters
    
    if (is_master) then
       write(ERstr,'(a,i0,a)') 'ER{',iface,'}'
       call write_matlab(echo,'method',ER%method,ERstr)
       call write_matlab(echo,'r' ,ER%r ,ERstr)
       call write_matlab(echo,'coupled' ,ER%coupled ,ERstr)
       call write_matlab(echo,'fluid_type',ER%FL%fluid_type,ERstr)
       call write_matlab(echo,'fraction',fraction,ERstr)
       call write_matlab(echo,'p0',ER%FL%p0,ERstr)
       call write_matlab(echo,'rho0',ER%FL%rho0,ERstr)
       call write_matlab(echo,'K0',ER%FL%K0,ERstr)
       call write_matlab(echo,'RT0',ER%FL%RT0,ERstr)
       call write_matlab(echo,'s',ER%FL%s,ERstr)
       call write_matlab(echo,'n0',ER%FL%n0,ERstr)
       call write_matlab(echo,'Tex',ER%FL%Tex,ERstr)
       call write_matlab(echo,'m',ER%FL%exp,ERstr)
       call write_matlab(echo,'mu',ER%FL%mu,ERstr)
       call write_matlab(echo,'g',ER%g,ERstr)
       call write_matlab(echo,'w0',w0,ERstr)
       call write_matlab(echo,'BCm',BCm,ERstr)
       call write_matlab(echo,'BCm_p',BCm_p,ERstr)
       call write_matlab(echo,'BCm_u',BCm_u,ERstr)
       call write_matlab(echo,'BCp',BCp,ERstr)
       call write_matlab(echo,'BCp_p',BCp_p,ERstr)
       call write_matlab(echo,'BCp_u',BCp_u,ERstr)
       call write_matlab(echo,'patm',patm,ERstr)
       call write_matlab(echo,'rhos',rhos,ERstr)
    end if

    ! return if not needed

    if (skip) return

    ! initialize

    ER%bndm = (m==mg)
    ER%bndp = (p==pg)
    ER%L = limits(nb,m,p,mg,pg,m-nb,p+nb,ER%bndm,ER%bndp)

    ! interior grid points at which weno will be used
    ER%mi = max(m,mg+nbst) ! special boundary FD stencils used at mg,...,mg+nbst-1
    ER%pi = min(p,pg-nbst) ! special boundary FD stencils used at pg-nbst+1,...,pg
    ! edges surrounding interior grid points
    ER%me = ER%mi-1
    ER%pe = ER%pi

    allocate(ER%q(ER%L%mb:ER%L%pb,nq),ER%Dq(ER%L%m:ER%L%p,nq))
    allocate(ER%rho(ER%L%mb:ER%L%pb),ER%u(ER%L%mb:ER%L%pb), &
         ER%p(ER%L%mb:ER%L%pb),ER%c(ER%L%mb:ER%L%pb), &
         ER%w(ER%L%mb:ER%L%pb),ER%dl(ER%L%m:ER%L%p), &
         ER%Smin(ER%L%m:ER%L%p),ER%n(ER%L%mb:ER%L%pb),ER%Dn(ER%L%m:ER%L%p))

    ER%dl = dl

    ER%rho = 1d40
    ER%u = 1d40
    ER%p = 1d40
    ER%c = 1d40
    ER%w = 1d40
    ER%q = 1d40
    ER%Dq = 1d40
    ER%n = 1d40
    ER%Dn = 1d40

    do i = ER%L%m,ER%L%p
       ! remote compression is fraction of lithostatic load
       ER%Smin(i) = fraction*(patm-rhos*ER%g*y(i))
       ! initial pressure is lithospheric load
       ER%p(i) = patm-rhos*ER%g*y(i)
       ! uniform width
       ER%w(i) = w0
    end do

    if (erupt_file) then
       ! both sides read file (so process may read file twice)
       if (process_m) call read_erupt(ER,filename,comm_m,array)
       if (process_p) call read_erupt(ER,filename,comm_p,array)
    end if

    do i = ER%L%m,ER%L%p
       ER%n(i) = n_eos(ER%FL,ER%p(i))
       ER%rho(i) = rho_eos(ER%FL,ER%p(i),ER%n(i))
       ER%c(i) = c_eos(ER%FL,ER%p(i),ER%n(i),ER%rho(i))
       ER%u(i) = 0d0
    end do

    call set_conservation_vars(ER)

    if (ER%bndm) then
       if (BCm(1:6)/='inflow') &
            call error('Need inflow BC at inlet','init_erupt')
       call init_bc(ER%BCm,BCm,nq,BCm_u,BCm_p)
    end if

    if (ER%bndp) then
       if (BCp(1:7)/='outflow') &
            call error('Need outflow BC at outlet','init_erupt')
       call init_bc(ER%BCp,BCp,nq,BCp_u,BCp_p)
    end if

  end subroutine init_erupt


  subroutine read_erupt(ER,filename,comm,array)

    use io, only : file_distributed,open_file_distributed, &
         read_file_distributed,close_file_distributed
    use mpi_routines, only : pw

    implicit none

    type(er_type),intent(inout) :: ER
    character(*),intent(in) :: filename
    integer,intent(in) :: comm,array

    type(file_distributed) :: fh

    call open_file_distributed(fh,filename,'read',comm,array,pw)

    call read_file_distributed(fh,ER%Smin(ER%L%m:ER%L%p))
    call read_file_distributed(fh,ER%p   (ER%L%m:ER%L%p))
    call read_file_distributed(fh,ER%w   (ER%L%m:ER%L%p))

    call close_file_distributed(fh)

  end subroutine read_erupt


  subroutine set_conservation_vars(ER)

    implicit none

    type(er_type),intent(inout) :: ER

    integer :: i

    do i = ER%L%m,ER%L%p
       ER%q(i,1) = ER%rho(i)*ER%w(i)
       ER%q(i,2) = ER%rho(i)*ER%u(i)*ER%w(i)
    end do

  end subroutine set_conservation_vars


  subroutine destroy_erupt(ER)

    implicit none

    type(er_type),intent(inout) :: ER

    if (allocated(ER%q   )) deallocate(ER%q)
    if (allocated(ER%Dq  )) deallocate(ER%Dq)
    if (allocated(ER%rho )) deallocate(ER%rho)
    if (allocated(ER%u   )) deallocate(ER%u)
    if (allocated(ER%p   )) deallocate(ER%p)
    if (allocated(ER%c   )) deallocate(ER%c)
    if (allocated(ER%w   )) deallocate(ER%w)
    if (allocated(ER%dl  )) deallocate(ER%dl)
    if (allocated(ER%Smin)) deallocate(ER%Smin)
    if (allocated(ER%n   )) deallocate(ER%n)
    if (allocated(ER%Dn  )) deallocate(ER%Dn)

    if (allocated(ER%BCm%q)) deallocate(ER%BCm%q)
    if (allocated(ER%BCp%q)) deallocate(ER%BCp%q)

  end subroutine destroy_erupt


  subroutine checkpoint_erupt(fh,operation,ER)

    use io, only : file_distributed, &
         read_file_distributed,write_file_distributed

    implicit none

    type(file_distributed),intent(in) :: fh
    character(*),intent(in) :: operation
    type(er_type),intent(inout) :: ER

    integer :: m,p

    ! fields read/written to same file as iface fields,
    ! routine is only called by io_process

    if (.not.ER%use_ER) return

    m = ER%L%m
    p = ER%L%p

    select case(operation)
    case('read')
       call read_file_distributed(fh,ER%Smin(m:p))
       call read_file_distributed(fh,ER%rho(m:p))
       call read_file_distributed(fh,ER%u(m:p))
       call read_file_distributed(fh,ER%p(m:p))
       call read_file_distributed(fh,ER%c(m:p))
       call read_file_distributed(fh,ER%w(m:p))
       call read_file_distributed(fh,ER%n(m:p))
       call set_conservation_vars(ER)
    case('write')
       call write_file_distributed(fh,ER%Smin(m:p))
       call write_file_distributed(fh,ER%rho(m:p))
       call write_file_distributed(fh,ER%u(m:p))
       call write_file_distributed(fh,ER%p(m:p))
       call write_file_distributed(fh,ER%c(m:p))
       call write_file_distributed(fh,ER%w(m:p))
       call write_file_distributed(fh,ER%n(m:p))
    end select

  end subroutine checkpoint_erupt


  subroutine init_bc(BC,bcnd,nq,u,p)

    implicit none

    type(bc_type),intent(out) :: BC
    character(*),intent(in) :: bcnd
    integer,intent(in) :: nq
    real,intent(in) :: u,p

    BC%bcnd = bcnd

    BC%u = u
    BC%p = p

    BC%u0 = 1d40
    BC%c0 = 1d40
    BC%M0 = 1d40
    BC%Z0 = 1d40
    BC%K  = 1d40

    allocate(BC%q(nq))
    BC%q = 1d40

  end subroutine init_bc


  subroutine scale_rates_erupt(ER,A)

    implicit none

    type(er_type),intent(inout) :: ER
    real,intent(in) :: A

    if (.not.ER%use_ER) return

    ER%Dq = A*ER%Dq
    ER%Dn = A*ER%Dn

  end subroutine scale_rates_erupt


  subroutine update_fields_erupt(ER,dt,DDn)

    implicit none

    type(er_type),intent(inout) :: ER
    real,intent(in) :: dt,DDn(:)

    integer :: i,m,p

    if (.not.ER%use_ER) return

    m = ER%L%m
    p = ER%L%p

    ! update width

    if (ER%coupled) ER%w(m:p) = ER%w(m:p)+dt*DDn

    ! return without updating q, etc., if solving for w given p-Smin

    if (ER%equilibrate) return

    ! update q and n, then convert to primitive variables
    
    do i = m,p
       ER%n(i) = ER%n(i)+dt*ER%Dn(i)
       ER%q(i,:) = ER%q(i,:)+dt*ER%Dq(i,:)
       ER%rho(i) = ER%q(i,1)/ER%w(i)
       ER%u(i) = ER%q(i,2)/ER%q(i,1)
       ER%p(i) = p_eos(ER%FL,ER%n(i),ER%rho(i))
       ER%c(i) = c_eos(ER%FL,ER%p(i),ER%n(i),ER%rho(i))
    end do

  end subroutine update_fields_erupt


  subroutine set_rates_erupt(ER,C,comm,t,m,p,y)

    use mpi_routines2d, only : cartesian
    use fd, only : diff

    implicit none

    type(er_type),intent(inout) :: ER
    type(cartesian),intent(in) :: C
    integer,intent(in) :: comm,m,p
    real,intent(in) :: t,y(m:p)

    integer :: i,j
    real :: tau
    real,dimension(:),allocatable :: Df,Dw,Dn
    real,dimension(:,:),allocatable :: f

    if (.not.ER%use_ER) return

    ! populate ghost cells of neighboring processors

    call share_erupt(ER,C)

    ! allocate auxiliary arrays

    allocate(f(ER%L%mb:ER%L%pb,2),Df(ER%L%m:ER%L%p), &
         Dw(ER%L%m:ER%L%p),Dn(ER%L%m:ER%L%p))

    ! evaluate flux at all grid points, including ghost points

    do i = ER%L%mb,ER%L%pb
       f(i,1) = ER%rho(i)*ER%u(i)*ER%w(i)
       f(i,2) = (ER%rho(i)*ER%u(i)**2+ER%p(i))*ER%w(i)
    end do

    ! SBP differentiation of flux

    do j = 1,2
       call diff(ER%L,f(:,j),Df)
       do i = ER%L%m,ER%L%p
          ER%Dq(i,j) = ER%Dq(i,j)-Df(i)/ER%dl(i)
       end do
    end do

    ! source term for variable width

    call diff(ER%L,ER%w,Dw)
    do i = ER%L%m,ER%L%p
       ER%Dq(i,2) = ER%Dq(i,2)+ER%p(i)*Dw(i)/ER%dl(i)
    end do

    ! wall shear stress and weight

    do i = ER%L%m,ER%L%p
       tau = wall_shear(ER%rho(i),ER%u(i),ER%w(i),ER%FL)
       ER%Dq(i,2) = ER%Dq(i,2)-tau-ER%rho(i)*ER%g*ER%w(i)
    end do

    ! add perturbation in momentum balance

    do i = ER%L%m,ER%L%p
       ER%Dq(i,2) = ER%Dq(i,2)+1d6*exp(-0.5d0*((y(i)+750d0)/10d0)**2)*exp(-0.5d0*((t-20d0)/0.01d0)**2)
    end do

    ! gas exsolution

    call diff(ER%L,ER%n,Dn)
    do i = ER%L%m,ER%L%p
       ER%Dn(i) = ER%Dn(i)-ER%u(i)*Dn(i)/ER%dl(i)- &
            (ER%n(i)-n_eos(ER%FL,ER%p(i)))/ER%FL%Tex
    end do

    ! add penalty terms to enforce BC with SAT method

    if (ER%bndm) then
       i = ER%L%mg
       call linearize_bc(ER%rho(i),ER%u(i),ER%c(i),ER%BCm)
       call bc_erupt(ER%rho(i),ER%u(i),ER%p(i),ER%w(i),ER%n(i),ER%FL,ER%BCm)
       ER%Dq(i,:) = ER%Dq(i,:)-ER%BCm%K/ER%dl(i)*(ER%q(i,:)-ER%BCm%q)
    end if

    if (ER%bndp) then
       i = ER%L%pg
       call linearize_bc(ER%rho(i),ER%u(i),ER%c(i),ER%BCp)
       call bc_erupt(ER%rho(i),ER%u(i),ER%p(i),ER%w(i),ER%n(i),ER%FL,ER%BCp)
       ER%Dq(i,:) = ER%Dq(i,:)-ER%BCp%K/ER%dl(i)*(ER%q(i,:)-ER%BCp%q)
    end if

    ! deallocate temporary arrays

    deallocate(f,Df,Dw,Dn)

  end subroutine set_rates_erupt


  function n_eos(FL,p) result(n)

    use io, only : error

    implicit none

    real,intent(in) :: p
    type(fluid),intent(in) :: FL
    real :: n

    select case(FL%fluid_type)
    case default
       call error('Invalid fluid_type','n_eos')
    case('liquid')
       n = 0d0
    case('gas')
       n = 1d0
    case('mixture')
       n = max(0d0,FL%n0-FL%s*P**FL%exp)
    end select

  end function n_eos


  function rho_eos(FL,p,n) result(rho)

    use io, only : error

    implicit none

    type(fluid),intent(in) :: FL
    real,intent(in) :: p,n
    real :: rho

    real :: rhoL,rhoG

    select case(FL%fluid_type)
    case default
       call error('Invalid fluid_type','rho_eos')
    case('liquid')
       rho = FL%rho0*(1d0+(p-FL%p0)/FL%K0)
    case('gas')
       rho = p/FL%RT0
    case('mixture')
       rhoL = FL%rho0*(1d0+(p-FL%p0)/FL%K0)
       rhoG = p/FL%RT0
       rho = 1d0/(n/rhoG+(1d0-n)/rhoL)
    end select

  end function rho_eos


  function p_eos(FL,n,rho) result(p)

    use io, only : error

    implicit none

    type(fluid),intent(in) :: FL
    real,intent(in) :: n,rho
    real :: p

    real :: a,b,c

    select case(FL%fluid_type)
    case default
       call error('Invalid fluid_type','p_eos')
    case('liquid')
       p = FL%p0+FL%K0*(rho/FL%rho0-1d0)
    case('gas')
       p = rho*FL%RT0
    case('mixture')
       if (rho<0d0) then
          print *, 'rho = ',rho
          call error('Negative density','p_eos')
       end if
       a = (FL%rho0/rho)*(FL%p0/FL%K0)
       b = (FL%rho0/rho)*(1d0-FL%p0/FL%K0)-n*FL%RT0*FL%rho0/FL%K0-(1d0-n)
       c = -n*FL%RT0*FL%rho0/FL%p0*(1d0-FL%p0/FL%K0)
       p = FL%p0*(-b+sqrt(b**2-4d0*a*c))/(2d0*a)
    end select

  end function p_eos


  function c_eos(FL,p,n,rho) result(c)
    
    use io, only : error
    
    implicit none

    type(fluid),intent(in) :: FL
    real,intent(in) :: p,n,rho
    real :: c

    real :: rhoG,rhoL,drhoGdp,drhoLdp,drhodp

    select case(FL%fluid_type)
    case default
       call error('Invalid fluid_type','drhodp_eos')
    case('liquid')
       c = sqrt(FL%K0/FL%rho0)
    case('gas')
       c = sqrt(FL%RT0)
    case('mixture')
       rhoG = p/FL%RT0
       drhoGdp = 1d0/FL%RT0
       rhoL = FL%rho0*(1d0+(p-FL%p0)/FL%K0)
       drhoLdp = FL%rho0/FL%K0
       drhodp  = rho**2*(n*drhoGdp/rhoG**2+(1d0-n)*drhoLdp/rhoL**2)
       c = 1d0/sqrt(drhodp)
    end select

  end function c_eos


  function wall_shear(rho,u,w,FL) result(tau)

    implicit none

    real,intent(in) :: rho,u,w
    type(fluid),intent(in) :: FL
    real :: tau

    tau = 3d0*FL%mu*u/w

  end function wall_shear


  subroutine linearize_bc(rho,u,c,BC)

    use fd_coeff, only : H00i
    use io, only : error

    implicit none

    real,intent(in) :: rho,u,c
    type(bc_type),intent(inout) :: BC

    BC%u0 = u
    BC%c0 = c
    BC%M0 = u/c
    BC%Z0 = rho*c

    ! SAT relaxation parameters

    select case(BC%bcnd)
    case('inflow-p','inflow-u')
       if     (BC%M0<=-1d0) then ! supersonic outflow (no BC)
          BC%K = 0d0 ! to set penalty term to zero
       elseif (BC%M0<  1d0) then ! subsonic (one BC)
          BC%K = (c+u)*H00i
       elseif (BC%M0>= 1d0) then ! supersonic inflow (two BC)
          BC%K = u*H00i
       end if
    case('outflow-p','outflow-u')
       if     (BC%M0<=-1d0) then ! supersonic inflow (two BC)
          BC%K = -u*H00i
       elseif (BC%M0<  1d0) then ! subsonic (one BC)
          BC%K = (c-u)*H00i
       elseif (BC%M0>= 1d0) then ! supersonic outflow (no BC)
          BC%K = 0d0 ! to set penalty term to zero
       end if
    case default
       call error('Invalid eruption BC','linearize_bc')
    end select

  end subroutine linearize_bc


  subroutine bc_erupt(rhoFD,uFD,pFD,w,n,FL,BC)

    use io, only : error

    implicit none

    real,intent(in) :: rhoFD,uFD,pFD,w,n
    type(fluid),intent(in) :: FL
    type(bc_type),intent(inout) :: BC

    real :: rho,u,p

    rho = 1d40
    u = 1d40
    p = 1d40

    select case(BC%bcnd)

    case default

       call error('Invalid BC','bc_erupt')

    case('inflow-p')

       if     (BC%M0<=-1d0) then ! supersonic outflow (no BC)
          u = uFD
          p = pFD
       elseif (BC%M0<  1d0) then ! subsonic (one BC)
          p = BC%p ! set p
          u = uFD+(p-pFD)/BC%Z0
       elseif (BC%M0>= 1d0) then ! supersonic inflow (two BC)
          p = BC%p ! set p
          u = BC%u ! set u
       end if

    case('inflow-u')

       if     (BC%M0<=-1d0) then ! supersonic outflow (no BC)
          u = uFD
          p = pFD
       elseif (BC%M0<  1d0) then ! subsonic (one BC)
          u = BC%u ! set u
          p = pFD+BC%Z0*(u-uFD)
       elseif (BC%M0>= 1d0) then ! supersonic inflow (two BC)
          p = BC%p ! set p
          u = BC%u ! set u
       end if

    case('outflow-p')

       if     (BC%M0<=-1d0) then ! supersonic inflow (two BC)
          p = BC%p ! set p
          u = BC%u ! set u
       elseif (BC%M0<  1d0) then ! subsonic (one BC)
          p = BC%p ! set p
          u = uFD-(p-pFD)/BC%Z0
       elseif (BC%M0>= 1d0) then ! supersonic outflow (no BC)
          u = uFD
          p = pFD
       end if

    case('outflow-u')

       if     (BC%M0<=-1d0) then ! supersonic inflow (two BC)
          p = BC%p ! set p
          u = BC%u ! set u
       elseif (BC%M0<  1d0) then ! subsonic (one BC)
          u = BC%u ! set u
          p = pFD-BC%Z0*(u-uFD)
       elseif (BC%M0>= 1d0) then ! supersonic outflow (no BC)
          u = uFD
          p = pFD
       end if

    end select

    rho = rho_eos(FL,p,n)

    BC%q = (/ rho*w,rho*u*w /)

  end subroutine bc_erupt


  subroutine pressure_erupt(ER,i,p)

    implicit none

    type(er_type),intent(in) :: ER
    integer,intent(in) :: i
    real,intent(inout) :: p

    ! solid only feels excess pressure (and solves for stress changes)

    if (ER%use_ER) p = p+ER%p(i)-ER%Smin(i)

  end subroutine pressure_erupt


  subroutine share_erupt(ER,C)

    use mpi_routines2d, only : cartesian,populate_ghost_cells

    implicit none

    type(ER_type),intent(inout) :: ER
    type(cartesian),intent(in) :: C

    integer :: j

    call populate_ghost_cells(C,ER%L%m,ER%L%p,ER%L%nb, &
         ER%bndm,ER%bndp,ER%rho,ER%direction)
    call populate_ghost_cells(C,ER%L%m,ER%L%p,ER%L%nb, &
         ER%bndm,ER%bndp,ER%u  ,ER%direction)
    call populate_ghost_cells(C,ER%L%m,ER%L%p,ER%L%nb, &
         ER%bndm,ER%bndp,ER%p  ,ER%direction)
    call populate_ghost_cells(C,ER%L%m,ER%L%p,ER%L%nb, &
         ER%bndm,ER%bndp,ER%c  ,ER%direction)
    call populate_ghost_cells(C,ER%L%m,ER%L%p,ER%L%nb, &
         ER%bndm,ER%bndp,ER%w  ,ER%direction)
    call populate_ghost_cells(C,ER%L%m,ER%L%p,ER%L%nb, &
         ER%bndm,ER%bndp,ER%n  ,ER%direction)

    do j = 1,nq
       call populate_ghost_cells(C,ER%L%m,ER%L%p,ER%L%nb, &
            ER%bndm,ER%bndp,ER%q(:,j),ER%direction)
    end do

  end subroutine share_erupt


  subroutine set_rates_erupt_old(ER,C,comm,t)

    use mpi_routines2d, only : cartesian
    use fd_coeff, only : DI,mI,pI,nbnd,nbst,reconstruct3,reconstruct5,reconstructW
    use fd, only : diff,diff_bndm,diff_bndp,diff_point
    use utilities, only : within
    use mpi_routines, only : MPI_REAL_PW
    use mpi

    implicit none

    ! .  |  .  |  .  |  .  |  .  |  .  |  .
    ! mg         i-1    i    i+1          pg (grid points)
    !    mg         i-1    i    i+1   pg-1   (edges)

    ! each process responsible for rates at
    ! grid points (cell centers) m:p with surrounding edges m-1:p

    type(er_type),intent(inout) :: ER
    type(cartesian),intent(in) :: C
    integer,intent(in) :: comm
    real,intent(in) :: t

    integer :: i,j,im,ip,stencil,ierr
    real :: uav,cav,alpha(nq),alpha_send(nq),tau
    real,dimension(:,:),allocatable :: f,fe,lambda,fm,fp
    real,dimension(:,:,:),allocatable :: R,Ri
    real,dimension(:),allocatable :: Df,Dw
    character(1) :: location

    if (.not.ER%use_ER) return

    ! populate ghost cells of neighboring processors

    call share_erupt(ER,C)

    ! allocate auxiliary arrays

    allocate(f(ER%L%mb:ER%L%pb,nq),Df(ER%L%m:ER%L%p),Dw(ER%L%m:ER%L%p))

    allocate(fe(ER%me:ER%pe,nq),lambda(ER%me:ER%pe,nq))
    allocate(R(ER%me:ER%pe,nq,nq),Ri(ER%me:ER%pe,nq,nq))
    allocate(fm(1-ER%r:ER%r,nq),fp(1-ER%r:ER%r,nq))

    ! evaluate flux at all grid points, including ghost points

    do i = ER%L%mb,ER%L%pb
       f(i,:) = flux(ER%rho(i),ER%u(i),ER%p(i),ER%w(i))
    end do

    ! SBP differentiation

    select case(ER%method)

    case('SBP') ! no upwinding

       do j = 1,nq
          call diff(ER%L,f(:,j),Df)
          ER%Dq(:,j) = ER%Dq(:,j)-Df/ER%dl
       end do

    case('SBPp') ! no upwinding, point-by-point

       ! minus boundary
       if (ER%bndm) then
          do j = 1,nq
             call diff_bndm(ER%L,f(ER%L%mg:ER%L%mg+nbnd-1,j), &
                  Df(ER%L%mg:ER%L%mg+nbst-1))
             ER%Dq(ER%L%mg:ER%L%mg+nbst-1,j) = ER%Dq(ER%L%mg:ER%L%mg+nbst-1,j)- &
                  Df(ER%L%mg:ER%L%mg+nbst-1)/ER%dl(ER%L%mg:ER%L%mg+nbst-1)
          end do
       end if

       ! plus  boundary
       if (ER%bndp) then
          do j = 1,nq
             call diff_bndp(ER%L,f(ER%L%pg-nbnd+1:ER%L%pg,j), &
                  Df(ER%L%pg-nbst+1:ER%L%pg))
             ER%Dq(ER%L%pg-nbst+1:ER%L%pg,j) = ER%Dq(ER%L%pg-nbst+1:ER%L%pg,j)- &
                  Df(ER%L%pg-nbst+1:ER%L%pg)/ER%dl(ER%L%pg-nbst+1:ER%L%pg)
          end do
       end if

       ! interior
       do i = max(ER%L%m,ER%L%mg+nbst),min(ER%L%p,ER%L%pg-nbst)
          do j = 1,nq
             Df(i) = dot_product(DI,f(i+mI:i+pI,j))
             ER%Dq(i,j) = ER%Dq(i,j)-Df(i)/ER%dl(i)
          end do
       end do

    case('SBPp2') ! no upwinding, point-by-point

       ! at each point to be updated,

       do i = ER%L%m,ER%L%p

          ! (1) define indices of points in stencil
          
          if     (within(ER%L%mg,i,ER%L%mg+nbst-1)) then
             im = ER%L%mg
             ip = ER%L%mg+nbnd-1
             location = 'm'
             stencil = i-ER%L%mg
          elseif (within(ER%L%pg-nbst+1,i,ER%L%pg)) then
             im = ER%L%pg-nbnd+1
             ip = ER%L%pg
             location = 'p'
             stencil = i-ER%L%pg
          else
             im = i+mI
             ip = i+pI
             location = 'i'
             stencil = 0
          end if

          ! (2) approximate derivative
          
          do j = 1,nq
             Df(i) = diff_point(f(im:ip,j),location,stencil)
             ER%Dq(i,j) = ER%Dq(i,j)-Df(i)/ER%dl(i)
          end do

       end do

    case('upwind') ! point-by-point, upwinding in interior

       ! preparation for upwinding in interior:
       ! note: edge i lies between grid points i and i+1

       ! at all cell edges of interior grid points,

       ! (1) average fields and perform local characteristic decomposition

       do i = ER%me,ER%pe
          call average_fields(ER%u(i:i+1),ER%c(i:i+1),uav,cav)
          call characteristics(uav,cav,lambda(i,:),R(i,:,:),Ri(i,:,:))
       end do
       
       ! (2) determine global Lax-Friedrichs splitting parameter for each characteristic variable
       
       alpha = maxval(abs(lambda),1)
       ! MPI-1 version
       alpha_send = alpha
       call MPI_Allreduce(alpha_send,alpha,nq,MPI_REAL_PW,MPI_MAX,comm,ierr)
       ! MPI-2 version
       !call MPI_Allreduce(MPI_IN_PLACE,alpha,nq,MPI_REAL_PW,MPI_MAX,comm,ierr)

       ! (3) split flux at cell edge

       do i = ER%me,ER%pe

          ! split flux (characteristic variables) at neighboring grid points
             
          do j = 1-ER%r,ER%r
             call split_flux(ER%q(i+j,:),f(i+j,:),Ri(i,:,:),alpha,fm(j,:),fp(j,:))
          end do
          
          ! reconstruct flux (characteristic variables) at cell edge
          
          select case(ER%r)
          case(1)
             fe(i,:) = fp(0,:)+fm(1,:)
          case(2)
             fe(i,1) = reconstruct3(1,fp(-1:1,1))+reconstruct3(0,fm(0:2,1))
             fe(i,2) = reconstruct3(1,fp(-1:1,2))+reconstruct3(0,fm(0:2,2))
          case(3)
             fe(i,1) = reconstruct5(2,fp(-2:2,1))+reconstruct5(1,fm(-1:3,1))
             fe(i,2) = reconstruct5(2,fp(-2:2,2))+reconstruct5(1,fm(-1:3,2))
             !fe(i,1) = reconstructW(2,fp(-2:2,1))+reconstructW(1,fm(-1:3,1))
             !fe(i,2) = reconstructW(2,fp(-2:2,2))+reconstructW(1,fm(-1:3,2))
          end select
          
          ! transform flux from characteristic to physical variables
          
          fe(i,:) = matmul(R(i,:,:),fe(i,:))
          
       end do

       ! at each point to be updated,

       do i = ER%L%m,ER%L%p

          ! (1) define indices of points in stencil
          
          if     (within(ER%L%mg,i,ER%L%mg+nbst-1)) then
             im = ER%L%mg
             ip = ER%L%mg+nbnd-1
             location = 'm'
             stencil = i-ER%L%mg
          elseif (within(ER%L%pg-nbst+1,i,ER%L%pg)) then
             im = ER%L%pg-nbnd+1
             ip = ER%L%pg
             location = 'p'
             stencil = i-ER%L%pg
          else
             im = i+mI
             ip = i+pI
             location = 'i'
             stencil = 0
          end if

          ! (2) approximate derivative

          select case(location)

          case default ! boundaries, so use linear SBP operators

             do j = 1,nq
                Df(i) = diff_point(f(im:ip,j),location,stencil)
                ER%Dq(i,j) = ER%Dq(i,j)-Df(i)/ER%dl(i)
             end do

          case('i') ! interior, so upwind

             ER%Dq(i,:) = ER%Dq(i,:)-(fe(i,:)-fe(i-1,:))/ER%dl(i)

          end select
             
       end do

    end select

    ! add source term for variable width

    call diff(ER%L,ER%w,Dw)
    ER%Dq(:,2) = ER%Dq(:,2)+ER%p(ER%L%m:ER%L%p)*Dw/ER%dl

    ! add source term for drag

    do i = ER%L%m,ER%L%p
       tau = wall_shear(ER%rho(i),ER%u(i),ER%w(i),ER%FL)
       ER%Dq(i,2) = ER%Dq(i,2)-tau
    end do

    ! add source term for gravity
    
    if (ER%g/=0d0) ER%Dq(:,2) = ER%Dq(:,2)- &
         ER%g*ER%rho(ER%L%m:ER%L%p)*ER%w(ER%L%m:ER%L%p)

    ! add perturbation in momentum balance

    !do i = ER%L%m,ER%L%p
    !   ER%Dq(i,2) = ER%Dq(i,2)+1d5*exp(-0.25d0*dble(i-20)**2)*exp(-0.04d0*(t-100d0)**2)
    !end do

    ! add penalty terms to enforce BC with SAT method

    if (ER%bndm) then
       i = ER%L%mg
       call linearize_bc(ER%rho(i),ER%u(i),ER%c(i),ER%BCm)
       call bc_erupt(ER%rho(i),ER%u(i),ER%p(i),ER%w(i),ER%n(i),ER%FL,ER%BCm)
       ER%Dq(i,:) = ER%Dq(i,:)-ER%BCm%K/ER%dl(i)*(ER%q(i,:)-ER%BCm%q)
    end if

    if (ER%bndp) then
       i = ER%L%pg
       call linearize_bc(ER%rho(i),ER%u(i),ER%c(i),ER%BCp)
       call bc_erupt(ER%rho(i),ER%u(i),ER%p(i),ER%w(i),ER%n(i),ER%FL,ER%BCp)
       ER%Dq(i,:) = ER%Dq(i,:)-ER%BCp%K/ER%dl(i)*(ER%q(i,:)-ER%BCp%q)
    end if

    ! deallocate temporary arrays

    deallocate(fe,lambda,R,Ri,fm,fp)
    deallocate(f,Df,Dw)

  end subroutine set_rates_erupt_old


  function flux(rho,u,p,w) result(f)
    
    implicit none

    real,intent(in) :: rho,u,p,w
    real :: f(nq)

    f(1) = rho*u*w
    f(2) = (rho*u**2+p)*w
    
  end function flux

  
  subroutine average_fields(u,c,uav,cav)

    implicit none

    real,intent(in) :: u(2),c(2)
    real,intent(out) :: uav,cav

    uav = 0.5d0*sum(u)
    cav = 0.5d0*sum(c)

  end subroutine average_fields


  subroutine characteristics(u,c,lambda,R,Ri)

    implicit none

    real,intent(in) :: u,c
    real,intent(out) :: lambda(nq),R(nq,nq),Ri(nq,nq)

    lambda = (/ u-c,u+c /)

    R(1,1) = 1d0
    R(1,2) = 1d0
    R(2,1) = u-c
    R(2,2) = u+c
    
    Ri(1,1) = c+u
    Ri(1,2) = -1d0
    Ri(2,1) = c-u
    Ri(2,2) = 1d0
    Ri = Ri/(2d0*c)

  end subroutine characteristics


  subroutine split_flux(q,f,Ri,alpha,fm,fp)

    implicit none

    real,intent(in) :: q(:),f(:),Ri(:,:),alpha(:)
    real,intent(out) :: fm(:),fp(:)

    real :: qhat(size(q)),fhat(size(f))

    ! transform to local characteristic variables
  
    qhat = matmul(Ri,q)
    fhat = matmul(Ri,f)
  
    ! split flux (scalar Lax-Friedrichs on each characteristic variable)
    
    fp = 0.5d0*(fhat+alpha*qhat)
    fm = 0.5d0*(fhat-alpha*qhat)
  
  end subroutine split_flux


end module erupt
