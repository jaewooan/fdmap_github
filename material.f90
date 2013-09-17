module material

  implicit none

  type :: elastic_type
     real,dimension(:,:),allocatable :: rho,cs,cp,G,M,gamma,nu
  end type elastic_type

  type :: pml_type
    logical :: pmlx,pmly
    real :: n ! Exponent
    real :: sigma ! pml sigma weight
    real :: kappa ! complex frequency shift
    real :: mx, px, my, py ! start location of the pml in minus x, plus x, etc.., direction
  end type pml_type

  type :: block_material
     character(256) :: response
     real :: rho,cs,cp, &
          G,nu,lambda,K,M,gamma,Zs,Zsi,Zp,Zpi, &
          beta,mu,b,h,eta
     type(pml_type) :: pml
     logical :: heterogeneous,plastic_strain_tensor
  end type block_material

contains


  subroutine init_material(iblock,M,E,C,input,echo)

    use mpi_routines, only : is_master
    use mpi_routines2d, only : cartesian
    use io, only : error,write_matlab,seek_to_string

    implicit none

    integer,intent(in) :: iblock,input,echo
    type(block_material),intent(out) :: M
    type(elastic_type),intent(inout) :: E
    type(cartesian),intent(in) :: C

    integer :: stat
    character(256) :: Mstr,str,response,filename
    real :: rho,G,cs,cp,beta,mu,b,h,eta
    logical :: heterogeneous,plastic_strain_tensor

    namelist /material_list/ response,heterogeneous,filename,rho,G,cs,cp,beta,mu,b,h,eta,plastic_strain_tensor

    ! defaults
    
    response = 'elastic'
    heterogeneous = .false.
    plastic_strain_tensor = .false.
    filename = ''
    rho = 0d0
    G = 0d0
    cs = 0d0
    cp = 0d0
    beta = 0d0
    mu = 1d0
    b = 0d0
    h = 0d0
    eta = 0d0

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

    if (heterogeneous.and.(.not.allocated(E%rho))) &
         call init_medium(filename,C,E)

    if (cs==0d0.and.(.not.heterogeneous)) call error('Must specify S-wave speed','init_material')
    
    ! set cp assuming Poisson solid if not otherwise set on input
    
    if (cp==0d0) cp = sqrt(3d0)*cs
    
    ! must specify either G or rho in addition to cs and cp
    
    if (rho==0d0) rho = G/cs**2
    if (G==0d0) G = rho*cs**2

    ! store parameters

    M%response = response
    M%heterogeneous = heterogeneous
    M%plastic_strain_tensor = plastic_strain_tensor
    M%rho = rho
    M%G = G
    M%cs = cs
    M%cp = cp
    M%beta = beta
    M%mu = mu
    M%b = b
    M%h = h
    M%eta = eta

    ! other parameters, calculated as a function of input parameters

    M%M = M%rho*M%cp**2
    M%lambda = M%M-2d0*M%G
    M%nu = 0.5d0*M%lambda/(M%lambda+M%G)
    M%K = M%lambda+2d0*M%G/3d0
    M%gamma = 1d0-2d0*(M%cs/M%cp)**2
    M%Zs = M%rho*M%cs
    M%Zsi = 1d0/M%Zs ! can be NaN for fluid
    M%Zp = M%rho*M%cp
    M%Zpi = 1d0/M%Zp

    ! output material parameters

    if (is_master) then
       write(Mstr,'(a,i0,a)') 'M{',iblock,'}'
       call write_matlab(echo,'response',M%response,Mstr)
       call write_matlab(echo,'heterogeneous',M%heterogeneous,Mstr)
       if (.not.M%heterogeneous) then
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
       end if
       if (M%response=='plastic') then
          call write_matlab(echo,'beta',M%beta,Mstr)
          call write_matlab(echo,'mu',M%mu,Mstr)
          call write_matlab(echo,'b',M%b,Mstr)
          call write_matlab(echo,'h',M%h,Mstr)
          call write_matlab(echo,'eta',M%eta,Mstr)
       end if
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
    real :: sigma, kappa, n, mx, px, my, py

    namelist /pml_list/ pmlx, pmly, n, sigma, kappa, mx, px, my, py

    ! defaults

    pmlx  = .false.
    pmly  = .false.

    sigma = 0d0
    kappa = 0d0

    mx = -huge(0d0)
    px =  huge(0d0)

    my = -huge(0d0)
    py =  huge(0d0)

    ! read in material parameters
    write(str,'(a,i0,a)') '!---BLOCK',iblock,'---'
    call seek_to_string(input,str)
    read(input,nml=pml_list,iostat=stat)
    if (stat>0) call error('Error in pml_list','init_pml')

    ! store parameters
    M%pml%n = n

    M%pml%pmlx = pmlx
    M%pml%pmly = pmly

    M%pml%sigma = sigma
    M%pml%kappa = kappa

    M%pml%mx = mx
    M%pml%px = px

    M%pml%my = my
    M%pml%py = py

    ! output material parameters

    if (is_master) then
       write(Mstr,'(a,i0,a)') 'M{',iblock,'}'
       call write_matlab(echo,'pmlx',M%pml%pmlx,Mstr)
       call write_matlab(echo,'pmly',M%pml%pmly,Mstr)

       call write_matlab(echo,'pml_n',M%pml%n,Mstr)
       call write_matlab(echo,'pml_sigma',M%pml%sigma,Mstr)
       call write_matlab(echo,'pml_kappa',M%pml%kappa,Mstr)

       call write_matlab(echo,'pml_mx',M%pml%mx,Mstr)
       call write_matlab(echo,'pml_px',M%pml%px,Mstr)
       call write_matlab(echo,'pml_my',M%pml%my,Mstr)
       call write_matlab(echo,'pml_py',M%pml%py,Mstr)
    end if

  end subroutine init_pml


  subroutine init_medium(filename,C,E)

    use mpi_routines2d, only : cartesian,exchange_all_neighbors

    implicit none

    character(*),intent(in) :: filename
    type(cartesian),intent(in) :: C
    type(elastic_type),intent(inout) :: E

    integer :: i,j

    ! initialize heterogeneous medium

    allocate(E%rho  (C%mbx:C%pbx,C%mby:C%pby))
    allocate(E%cs   (C%mbx:C%pbx,C%mby:C%pby))
    allocate(E%cp   (C%mbx:C%pbx,C%mby:C%pby))
    allocate(E%G    (C%mbx:C%pbx,C%mby:C%pby))
    allocate(E%M    (C%mbx:C%pbx,C%mby:C%pby))
    allocate(E%gamma(C%mbx:C%pbx,C%mby:C%pby))
    allocate(E%nu   (C%mbx:C%pbx,C%mby:C%pby))

    E%rho = 1d40
    E%cs = 1d40
    E%cp = 1d40
    E%G = 1d40
    E%M = 1d40
    E%gamma = 1d40
    E%nu = 1d40

    call mediumIO('read',filename,C,E)

    do j = C%my,C%py
       do i = C%mx,C%px
          E%G(i,j) = E%rho(i,j)*E%cs(i,j)**2
          E%M(i,j) = E%rho(i,j)*E%cp(i,j)**2
          E%gamma(i,j) = 1d0-2d0*(E%cs(i,j)/E%cp(i,j))**2
          E%nu(i,j) = 0.5d0*(E%M(i,j)-2d0*E%G(i,j))/(E%M(i,j)-E%G(i,j))
       end do
    end do

    call exchange_all_neighbors(C,E%rho)
    call exchange_all_neighbors(C,E%cs)
    call exchange_all_neighbors(C,E%cp)
    call exchange_all_neighbors(C,E%G)
    call exchange_all_neighbors(C,E%M)
    call exchange_all_neighbors(C,E%gamma)
    call exchange_all_neighbors(C,E%nu)

  end subroutine init_medium


  subroutine mediumIO(operation,filename,C,E)

    use io, only : file_distributed,open_file_distributed, &
         read_file_distributed,write_file_distributed,close_file_distributed
    use mpi_routines2d, only : cartesian
    use mpi_routines, only : pw,is_master
    use mpi

    implicit none

    character(*),intent(in) :: operation,filename
    type(cartesian),intent(in) :: C
    type(elastic_type),intent(inout) :: E

    type(file_distributed) :: fh
    integer :: ierr

    if (operation=='delete') then
       if (is_master) call MPI_file_delete(filename,MPI_INFO_NULL,ierr)
       return
    end if

    call open_file_distributed(fh,filename,operation,C%c2d%comm,C%c2d%array_w,pw)

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    select case(operation)
    case('read')
       call  read_file_distributed(fh,E%rho(C%mx:C%px,C%my:C%py))
       call  read_file_distributed(fh,E%cs (C%mx:C%px,C%my:C%py))
       call  read_file_distributed(fh,E%cp (C%mx:C%px,C%my:C%py))
    case('write')
       call write_file_distributed(fh,E%rho(C%mx:C%px,C%my:C%py))
       call write_file_distributed(fh,E%cs (C%mx:C%px,C%my:C%py))
       call write_file_distributed(fh,E%cp (C%mx:C%px,C%my:C%py))
    end select

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    call close_file_distributed(fh)

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

  end subroutine mediumIO


  subroutine destroy_elastic(E)

    implicit none

    type(elastic_type),intent(inout) :: E

    if (allocated(E%rho  )) deallocate(E%rho  )
    if (allocated(E%cs   )) deallocate(E%cs   )
    if (allocated(E%cp   )) deallocate(E%cp   )
    if (allocated(E%G    )) deallocate(E%G    )
    if (allocated(E%M    )) deallocate(E%M    )
    if (allocated(E%gamma)) deallocate(E%gamma)
    if (allocated(E%nu   )) deallocate(E%nu   )

  end subroutine destroy_elastic


end module material
