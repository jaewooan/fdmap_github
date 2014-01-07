module thermpres

  implicit none

  type :: tp_type
     logical :: use_TP
     integer :: n,nsubsteps
     real :: h,alphath,alphahy,Lambda,rhoc,w,theta
     real,dimension(:),allocatable :: z
     real,dimension(:,:),allocatable :: g0,T,p
     real,dimension(:,:,:),allocatable :: M0,Rthi,Rhyi,Rthe,Rhye
  end type tp_type


contains
  

  subroutine init_thermpres(iface,TP,m,p,input,echo,skip,refine,dt)
    
    use mpi_routines, only : is_master
    use io, only : error,write_matlab,seek_to_string  

    implicit none

    integer,intent(in) :: iface,m,p,input,echo
    type(tp_type),intent(out) :: TP
    logical,intent(in) :: skip
    real,intent(in) :: refine,dt(:)

    integer :: stat
    character(256) :: TPstr
    character(256) :: str

    integer :: n,k,nsubsteps
    real :: h,alphath,alphahy,Lambda,rhoc,w,T0,p0,theta

    namelist /thermpres_list/ n,h,alphath,alphahy,Lambda,rhoc,w,T0,p0,theta,nsubsteps

    ! defaults
       
    n = 0
    nsubsteps = 1
    h = 1d0
    alphath = 0d0
    alphahy = 0d0
    Lambda = 0d0
    rhoc = 1d0
    w = 1d0
    T0 = 0d0
    p0 = 0d0
    theta = 1d0
    ! theta=0   (explicit Euler)
    ! theta=0.5 (Crank-Nicolson)
    ! theta=1   (implicit Euler)

    ! read in thermal pressurization parameters

    write(str,'(a,i0,a)') '!---IFACE',iface,'---'
    call seek_to_string(input,str)
    read(input,nml=thermpres_list,iostat=stat)
    if (stat>0) call error('Error in thermpres_list','init_thermpres')

    TP%use_TP = (n/=0)

    if (.not.TP%use_TP) return

    TP%n = n
    if (TP%n/=1) TP%n = ceiling(dble(TP%n-1)*refine)+1
    TP%h = h/refine
    TP%nsubsteps = nsubsteps
    TP%alphath = alphath
    TP%alphahy = alphahy
    TP%Lambda = Lambda
    TP%rhoc = rhoc
    TP%w = w
    TP%theta = theta

    ! output thermal pressurization parameters
    
    if (is_master) then
       write(TPstr,'(a,i0,a)') 'TP{',iface,'}'
       call write_matlab(echo,'n',TP%n,TPstr)
       call write_matlab(echo,'nsubsteps',TP%nsubsteps,TPstr)
       call write_matlab(echo,'h',TP%h,TPstr)
       call write_matlab(echo,'alphath',TP%alphath,TPstr)
       call write_matlab(echo,'alphahy',TP%alphahy,TPstr)
       call write_matlab(echo,'Lambda',TP%Lambda,TPstr)
       call write_matlab(echo,'rhoc',TP%rhoc,TPstr)
       call write_matlab(echo,'w',TP%w,TPstr)
       call write_matlab(echo,'theta',TP%theta,TPstr)
    end if

    ! return if not needed

    if (skip) return
    
    ! allocate and initialize arrays

    allocate(TP%z(TP%n))
    allocate(TP%T(TP%n,m:p),TP%p(TP%n,m:p))

    TP%z = dble( (/ (k, k=0,TP%n-1) /) )*TP%h
    TP%T = T0
    TP%p = p0

    call init_solver(TP,dt/dble(nsubsteps))

  end subroutine init_thermpres


  subroutine init_solver(TP,dt)

    implicit none

    type(tp_type),intent(inout) :: TP
    real,dimension(:),intent(in) :: dt

    integer :: j,s,nstage
    real :: r
    real,dimension(:,:),allocatable :: I,D2,Rthi,Rhyi,Rthe,Rhye
    real,parameter :: sqrt2pi = 2.506628274631d0

    ! number of stages

    nstage = size(dt)

    ! arrays for linear system solver

    allocate( &
         TP%Rthi(3,TP%n,nstage),TP%Rhyi(3,TP%n,nstage), &
         TP%Rthe(3,TP%n,nstage),TP%Rhye(3,TP%n,nstage))
    allocate(TP%M0(TP%n,TP%n,nstage),TP%g0(TP%n,nstage))

    ! temporary arrays used when setting up

    allocate(I(TP%n,TP%n),D2(TP%n,TP%n))
    allocate( &
         Rthi(TP%n,TP%n),Rhyi(TP%n,TP%n), &
         Rthe(TP%n,TP%n),Rhye(TP%n,TP%n))

    ! identity matrix

    I = 0d0
    forall (j=1:TP%n) I(j,j) = 1d0

    ! second derivative matrix

    call D2matrix(TP%n,TP%h,D2)

    ! set up matrices for all stages

    do s = 1,nstage

       ! diffusion operators (i=implicit, e=explicit)
       
       r = dt(s)*TP%alphahy/TP%h**2
       Rhyi = I-TP%theta*r*D2
       Rhye = I+(1d0-TP%theta)*r*D2
       
       r = dt(s)*TP%alphath/TP%h**2
       Rthi = I-TP%theta*r*D2
       Rthe = I+(1d0-TP%theta)*r*D2

       ! product of diffusion operators (pentadiagonal)
       
       TP%M0(:,:,s) = matmul(Rthi,Rhyi)
       
       ! shear heating source profile
       
       if (TP%w==0) then
          ! slip on plane model, source term from discretization of BC
          TP%g0(:,s) = 0d0
          TP%g0(1,s) = dt(s)/(TP%h*TP%rhoc)
       else
          TP%g0(:,s) = dt(s)/(sqrt2pi*TP%w*TP%rhoc)*exp(-0.5d0*(TP%z/TP%w)**2)
       end if

       ! store tridiagonal matrices in banded format
       
       call pack_banded_matrix(TP%n,1,1,Rthi,TP%Rthi(:,:,s))
       call pack_banded_matrix(TP%n,1,1,Rthe,TP%Rthe(:,:,s))
       call pack_banded_matrix(TP%n,1,1,Rhyi,TP%Rhyi(:,:,s))
       call pack_banded_matrix(TP%n,1,1,Rhye,TP%Rhye(:,:,s))
       
    end do

    ! deallocate temporary arrays

    deallocate(I,D2)
    deallocate(Rthi,Rhyi,Rthe,Rhye)

  end subroutine init_solver


  subroutine destroy_thermpres(TP)
    
    implicit none

    type(tp_type),intent(inout) :: TP

    if (allocated(TP%z)) deallocate(TP%z)
    if (allocated(TP%T)) deallocate(TP%T)
    if (allocated(TP%p)) deallocate(TP%p)
    if (allocated(TP%g0)) deallocate(TP%g0)
    if (allocated(TP%M0)) deallocate(TP%M0)
    if (allocated(TP%Rthi)) deallocate(TP%Rthi)
    if (allocated(TP%Rhyi)) deallocate(TP%Rhyi)
    if (allocated(TP%Rthe)) deallocate(TP%Rthe)
    if (allocated(TP%Rhye)) deallocate(TP%Rhye)

  end subroutine destroy_thermpres


  subroutine checkpoint_thermpres(operation,name,checkpoint_number, &
       iface,side,TP,comm,io_process)

    use io, only : file_distributed,open_file_distributed, &
         read_file_distributed,write_file_distributed,close_file_distributed
    use mpi_routines, only : MPI_REAL_PW,pw,subarray,is_master
    use mpi

    implicit none

    character(*),intent(in) :: operation,name
    integer,intent(in) :: checkpoint_number,iface,side
    type(tp_type),intent(inout) :: TP
    integer,intent(in),optional :: comm
    logical,intent(in),optional :: io_process

    type(file_distributed) :: fh
    character(256) :: filename
    integer :: array,ierr

    ! must be called by all processes in MPI_COMM_WORLD,
    ! as master process deletes files

    if (.not.TP%use_TP) return

    write(filename,'(a,i0,a,i0,a,i0)') trim(adjustl(name)) // &
         'iface',iface,'side',side,'TP.ckpt',checkpoint_number

    if (operation=='delete') then
       if (is_master) call MPI_File_delete(filename,MPI_INFO_NULL,ierr)
       return
    end if

    if (.not.io_process) return

    stop 'checkpoint_thermpres not functional yet'
    !call subarray(TP%n,pg-mg+1,1,TP%n,m-mg+1,p-mg+1,MPI_REAL_PW,array)

    call open_file_distributed(fh,filename,operation,comm,array,pw)

    select case(operation)
    case('read')
       call read_file_distributed(fh,TP%T)
       call read_file_distributed(fh,TP%p)
    case('write')
       call write_file_distributed(fh,TP%T)
       call write_file_distributed(fh,TP%p)
    end select

    call close_file_distributed(fh)

  end subroutine checkpoint_thermpres


  subroutine update_fields_thermpres(TP,m,p,f,V,sigma,stage)

    implicit none

    type(tp_type),intent(inout) :: TP
    integer,intent(in) :: m,p,stage
    real,dimension(m:p),intent(in) :: f,V,sigma
    
    integer :: i,n

    if (.not.TP%use_TP) return

    ! DEBUG
    !i = m
    !TP%p(:,i) = 100d0
    !TP%T(:,i) = 100d0
    !call solve_thermpres_explicit(TP,1d0,1d0,200d0,TP%p(:,i),TP%T(:,i),stage)
    !call solve_thermpres_implicit(TP,1d0,1d0,200d0,TP%p(:,i),TP%T(:,i),stage)
    !call write_matrix(TP%n,2,(/TP%p(:,i),TP%T(:,i)/))
    !stop
    ! END DEBUG

    do n = 1,TP%nsubsteps
       do i = m,p
          if (TP%theta==0d0) then
             call solve_thermpres_explicit(TP,f(i),V(i),sigma(i),TP%p(:,i),TP%T(:,i),stage)
          else
             call solve_thermpres_implicit(TP,f(i),V(i),sigma(i),TP%p(:,i),TP%T(:,i),stage)
          end if
       end do
    end do

  end subroutine update_fields_thermpres


  subroutine solve_thermpres_explicit(TP,f,V,sigma,p,T,s)

    use io, only : error

    implicit none

    type(tp_type),intent(inout) :: TP
    real,intent(in) :: f,V,sigma
    real,dimension(:),intent(inout) :: p,T
    integer,intent(in) :: s

    real :: a
    real,dimension(:),allocatable :: g,qth,qhy

    ! special case for Lambda=0

    if (TP%Lambda==0d0) call error('Lambda=0 case not yet implemented','solve_thermpres')
    
    ! allocate arrays

    allocate(g(TP%n),qth(TP%n),qhy(TP%n))

    ! shear heating

    g = f*V*TP%g0(:,s)

    ! update T (and temporarily store Tnew in qth)

    ! qth = matmul(Rthe,T)+(sigma-p(1))*g
    qth = g
    a = sigma-p(1)
    call dgbmv('n',TP%n,TP%n,1,1,1d0,TP%Rthe(:,:,s),3,T,1,a,qth,1)

    ! update p (and temporarily store pnew in qhy)

    ! qhy = matmul(Rhye,p)+Lambda*(Tnew-Told)
    qhy = qth-T ! Tnew-Told
    a = TP%Lambda
    call dgbmv('n',TP%n,TP%n,1,1,1d0,TP%Rhye(:,:,s),3,p,1,a,qhy,1)

    ! set T and p

    T = qth
    p = qhy

    ! deallocate arrays
    
    deallocate(g,qth,qhy)

  end subroutine solve_thermpres_explicit


  subroutine solve_thermpres_implicit(TP,f,V,sigma,p,T,s)

    use io, only : error

    implicit none

    type(tp_type),intent(inout) :: TP
    real,intent(in) :: f,V,sigma
    real,dimension(:),intent(inout) :: p,T
    integer,intent(in) :: s

    integer :: info
    integer,dimension(:),allocatable :: ipvt
    real :: a
    real,dimension(:),allocatable :: g,qth,qhy
    real,dimension(:,:),allocatable :: M

    ! special case for Lambda=0

    if (TP%Lambda==0d0) call error('Lambda=0 case not yet implemented','solve_thermpres')
    
    ! allocate arrays

    allocate(g(TP%n),M(TP%n,TP%n),ipvt(TP%n),qth(TP%n),qhy(TP%n))

    ! shear heating

    g = f*V*TP%g0(:,s)

    ! matrix in linear system for p

    M = TP%M0(:,:,s)
    M(:,1) = M(:,1)+TP%Lambda*TP%theta*g

    ! right-hand sides

    ! qhy = matmul(Rhye,p)-TP%Lambda*T
    qhy = T
    a = -TP%Lambda
    call dgbmv('n',TP%n,TP%n,1,1,1d0,TP%Rhye(:,:,s),3,p,1,a,qhy,1)

    ! qth = matmul(Rthe,T)+(sigma-(1d0-TP%theta)*p(1))*g
    qth = g
    a = sigma-(1d0-TP%theta)*p(1)
    call dgbmv('n',TP%n,TP%n,1,1,1d0,TP%Rthe(:,:,s),3,T,1,a,qth,1)

    ! solve linear system for p

    !p = matmul(Rth,qhy)+TP%Lambda*qth
    p = qth
    a = TP%Lambda
    call dgbmv('n',TP%n,TP%n,1,1,1d0,TP%Rthi(:,:,s),3,qhy,1,a,p,1)

    call dgesv(TP%n,1,M,TP%n,ipvt,p,TP%n,info)

    ! evaluate T

    !T = (1d0/TP%Lambda)*matmul(Rhy,p)-(1d0/TP%Lambda)*qhy
    T = -qhy
    a = 1d0/TP%Lambda
    call dgbmv('n',TP%n,TP%n,1,1,a,TP%Rhyi(:,:,s),3,p,1,a,T,1)

    ! deallocate arrays
    
    deallocate(g,qth,qhy,M,ipvt)

  end subroutine solve_thermpres_implicit


  subroutine D2matrix(n,h,D2)

    implicit none

    integer,intent(in) :: n
    real,intent(in) :: h
    real,intent(out) :: D2(n,n)

    integer :: i,j

    ! second derivative matrix

    D2 = 0d0

    do i = 1,n
       D2(i,i) = -2d0
       do j = i-1,i+1,2
          if (1<=j.and.j<=n) D2(i,j) = 1d0
       end do
    end do

    ! homogeneous Neumann boundary conditions at both boundaries

    D2(1,2  ) = D2(1,2  )+1d0
    D2(n,n-1) = D2(n,n-1)+1d0

  end subroutine D2matrix


  subroutine write_matrix(m,n,A)

    implicit none

    integer,intent(in) :: m,n
    real,intent(in) :: A(m,n)

    integer :: i
    character(16) :: str

    !write(str,'(a,i0,a)') '(',n,'e18.10)'
    write(str,'(a,i0,a)') '(',n,'f10.4)'

    print *

    do i = 1,m
       write(6,str) A(i,:)
    end do

    print *

  end subroutine write_matrix


  subroutine pressure_thermpres(TP,i,p)

    implicit none

    type(tp_type),intent(in) :: TP
    integer,intent(in) :: i
    real,intent(inout) :: p

    if (TP%use_TP) p = p+TP%p(1,i)

  end subroutine pressure_thermpres


  subroutine pack_banded_matrix(n,kl,ku,A,M)

    implicit none

    integer,intent(in) :: n,kl,ku
    real,intent(in) :: A(n,n)
    real,intent(out) :: M(kl+ku+1,n)

    integer :: i,j,k 

    M = 1d40 ! to prevent accidental use of unused entries

    do j = 1,n
       k = ku+1-j
       do i = max(1,j-ku),min(n,j+kl)
          M(k+i,j) = A(i,j)
       end do
    end do

  end subroutine pack_banded_matrix


end module thermpres
