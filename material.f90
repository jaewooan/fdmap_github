module material

  implicit none

  type :: block_material
     character(256) :: response,pmltype
     real :: G,cs,cp, &
          rho,nu,lambda,K,M,gamma,m2,Zs,Zsi,Zp,Zpi, &
          beta,mu,b,h,eta,pmlsx,pmlsy,pmlax,pmlay,pmllx,pmlly
     logical :: pmlx,pmly
  end type block_material


contains


  subroutine init_material(iblock,M,input,echo)

    use mpi_routines, only : is_master
    use io, only : error,write_matlab,seek_to_string

    implicit none

    integer,intent(in) :: iblock,input,echo
    type(block_material),intent(out) :: M

    integer :: stat
    character(256) :: Mstr
    character(256) :: str
    character(256) :: response = 'elastic'
    real :: rho,G,cs,cp,beta,mu,b,h,eta,eta_factor

    namelist /material_list/ response,rho,G,cs,cp,beta,mu,b,h,eta,eta_factor

    ! defaults

    rho = 0d0
    G = 0d0
    cs = 0d0
    cp = 0d0
    beta = 0d0
    mu = 1d0
    b = 0d0
    h = 0d0
    eta = 0d0
    eta_factor = 0d0

    ! read in material parameters

    write(str,'(a,i0,a)') '!---BLOCK',iblock,'---'
    call seek_to_string(input,str)
    read(input,nml=material_list,iostat=stat)
    if (stat>0) call error('Error in material_list','init_material')

    select case(response)
    case('elastic','plastic')
    case default
       call error('Invalid material response: ' // trim(response),'init_material')
    end select

    if (cs==0d0) call error('Must specify S-wave speed','init_material')

    ! set cp assuming Poisson solid if not otherwise set on input

    if (cp==0d0) cp = sqrt(3d0)*cs

    ! must specify either G or rho in addition to cs and cp

    if (rho==0d0) rho = G/cs**2
    if (G==0d0) G = rho*cs**2

    ! store parameters

    M%response = response
    M%rho = rho
    M%G = G
    M%cs = cs
    M%cp = cp
    M%beta = beta
    M%mu = mu
    M%b = b
    M%h = h

    if (eta_factor==0d0) then
       M%eta = eta
    else
       ! special for BSSA papers
       M%eta = 0.277482678983834d0*eta_factor
    end if

    ! other parameters, calculated as a function of input parameters

    M%M = M%rho*M%cp**2
    M%lambda = M%M-2d0*M%G
    M%nu = 0.5d0*M%lambda/(M%lambda+M%G)
    M%K = M%lambda+2d0*M%G/3d0
    M%m2 = (M%cs/M%cp)**2
    M%gamma = 1d0-2d0*M%m2
    M%Zs = M%rho*M%cs
    M%Zsi = 1d0/M%Zs ! can be NaN for fluid
    M%Zp = M%rho*M%cp
    M%Zpi = 1d0/M%Zp

    ! output material parameters

    if (is_master) then
       write(Mstr,'(a,i0,a)') 'M{',iblock,'}'
       call write_matlab(echo,'response',M%response,Mstr)
       call write_matlab(echo,'rho',M%rho,Mstr)
       call write_matlab(echo,'G',M%G,Mstr)
       call write_matlab(echo,'lambda',M%lambda,Mstr)
       call write_matlab(echo,'K',M%K,Mstr)
       call write_matlab(echo,'M',M%M,Mstr)
       call write_matlab(echo,'cs',M%cs,Mstr)
       call write_matlab(echo,'cp',M%cp,Mstr)
       call write_matlab(echo,'nu',M%nu,Mstr)
       call write_matlab(echo,'Zs',M%Zs,Mstr)
       call write_matlab(echo,'Zp',M%Zp,Mstr)
       call write_matlab(echo,'beta',M%beta,Mstr)
       call write_matlab(echo,'mu',M%mu,Mstr)
       call write_matlab(echo,'b',M%b,Mstr)
       call write_matlab(echo,'h',M%h,Mstr)
       call write_matlab(echo,'eta',M%eta,Mstr)
    end if

  end subroutine init_material

  subroutine init_pml(iblock,M,input,echo)

    use mpi_routines, only : is_master
    use io, only : error,write_matlab,seek_to_string

    implicit none

    integer,intent(in) :: iblock,input,echo
    type(block_material),intent(inout) :: M

    integer :: stat
    character(256) :: Mstr
    character(256) :: str
    character(256) :: pmltype = 'constant'
    logical :: pmlx, pmly
    real :: pmlsx,pmlsy,pmlax,pmlay,pmllx,pmlly

    namelist /pml_list/ pmltype,pmlx,pmly,pmlsx,pmlsy,pmlax,pmlay,pmllx,pmlly

    ! defaults

    pmlx  = .false.
    pmlsx = 0d0
    pmlax = 0d0
    pmllx = 0d0

    pmly  = .false.
    pmlsy = 0d0
    pmlay = 0d0
    pmlly = 0d0

    ! read in material parameters

    write(str,'(a,i0,a)') '!---BLOCK',iblock,'---'
    call seek_to_string(input,str)
    read(input,nml=pml_list,iostat=stat)
    if (stat>0) call error('Error in pml_list','init_pml')

    ! store parameters
    M%pmltype  = pmltype

    M%pmlx  = pmlx
    M%pmlsx = pmlsx
    M%pmlax = pmlax
    M%pmllx = pmllx

    M%pmly  = pmly
    M%pmlsy = pmlsy
    M%pmlay = pmlay
    M%pmlly = pmlly

    ! output material parameters

    if (is_master) then
       write(Mstr,'(a,i0,a)') 'M{',iblock,'}'
       call write_matlab(echo,'pmlx',M%pmlx,Mstr)
       call write_matlab(echo,'pmly',M%pmly,Mstr)
       call write_matlab(echo,'pmlsx',M%pmlsx,Mstr)
       call write_matlab(echo,'pmlsy',M%pmlsy,Mstr)
       call write_matlab(echo,'pmlax',M%pmlax,Mstr)
       call write_matlab(echo,'pmlay',M%pmlay,Mstr)
    end if

  end subroutine init_pml


end module material
