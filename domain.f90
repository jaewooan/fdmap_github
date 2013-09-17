module domain

  use material, only : block_material
  use grid, only : block_grid,grid_type
  use fields, only : block_fields,fields_type
  use boundaries, only : block_boundaries,iface_type
  use mpi_routines2d, only : cartesian
  use gravity, only : gravity_type

  implicit none

  type :: block_type
     logical :: dissipation
     real :: c_dissipation
     type(block_material) :: M
     type(block_grid) :: G
     type(block_fields) :: F
     type(block_boundaries) :: B
  end type block_type

  type :: domain_type
     logical :: operator_split,exact_metric
     integer :: mode,nblocks_x,nblocks_y,nblocks,nifaces,n
     character(10) :: FDmethod
     character(16) :: method
     real :: t,tD
     type(block_type),allocatable :: B(:)
     type(iface_type),allocatable :: I(:)
     type(grid_type) :: G
     type(fields_type) :: F
     type(cartesian) :: C
     type(gravity_type) :: grav
  end type domain_type


contains


  subroutine init_domain(D,t,refine,CFL,dt,input,echo,cRK)

    use material, only : init_material, init_pml
    use grid, only : read_grid,block_limits,init_grid,init_grid_partials,grid_spacing
    use fields, only : init_fields,exchange_fields,init_fields_side_rate
    use boundaries, only : init_iface,init_iface_blocks,init_iface_fields, &
         init_boundaries,apply_bc
    use fault, only : couple_blocks
    use fd_coeff, only : init_fd
    use gravity , only : init_gravity
    use mpi_routines2d, only : decompose2d,exchange_all_neighbors
    use mpi_routines, only : is_master
    use io, only : error,write_matlab,message,messages,seek_to_string

    implicit none

    type(domain_type),intent(out) :: D
    real,intent(in) :: t,refine,cRK(:)
    real,intent(inout) :: CFL,dt
    integer,intent(in) :: input,echo

    integer :: i,im,ip,mode,nblocks_x,nblocks_y,nblocks,nifaces, &
         nx,ny,nF,nprocs_x,nprocs_y,stat
    integer,dimension(:),allocatable :: blkxm,blkym
    real :: cs,c_dissipation
    real,allocatable :: dtRK(:)
    character(10) :: FDmethod
    character(6) :: mpi_method
    character(16) :: method
    character(256) :: str1,str2,str3,str
    logical :: operator_split,decomposition_info,energy_balance,displacement,exact_metric
    logical,parameter :: periodic_x=.false.,periodic_y=.false.

    namelist /domain_list/ mode,FDmethod,nblocks_x,nblocks_y,nblocks,nifaces, &
         nx,ny,mpi_method,nprocs_x,nprocs_y,decomposition_info,operator_split, &
         energy_balance,displacement,method,exact_metric

    namelist /operator_list/ c_dissipation

    ! defaults

    mode = 2
    FDmethod = 'SBP6' ! SBP2,SBP4,SBP6
    nblocks_x = 0
    nblocks_y = 0
    nblocks = 0
    nifaces = 0
    exact_metric = .false.

    nx = 1
    ny = 1

    operator_split = .false.

    energy_balance = .false.
    displacement   = .false.

    method = 'corrected'

    mpi_method = '2d'
    nprocs_x = 1
    nprocs_y = 1
    decomposition_info = .false.

    ! read in domain parameters

    rewind(input)
    read(input,nml=domain_list,iostat=stat)
    if (stat>0) call error('Error in domain_list','init_domain')

    if (nblocks_x==0) then
       nblocks_x = 1
       if (is_master) call message('Assuming that nblocks_x = 1')
    end if
    if (nblocks_y==0) then
       nblocks_y = 1
       if (is_master) call message('Assuming that nblocks_y = 1')
    end if
    if (nblocks==0) then
       nblocks = 1
       if (is_master) call message('Assuming that nblocks = 1')
    end if

    D%mode = mode
    D%FDmethod = FDmethod
    D%nblocks_x = nblocks_x
    D%nblocks_y = nblocks_y
    D%nblocks = nblocks
    D%nifaces = nifaces
    D%exact_metric = exact_metric

    D%operator_split = operator_split

    D%method = method

    ! initial time

    D%t = t

    ! refine number of cells (ncells = npoints-nblocks)

    if (nx/=1) nx = ceiling(dble(nx-D%nblocks_x)*refine)+D%nblocks_x
    if (ny/=1) ny = ceiling(dble(ny-D%nblocks_y)*refine)+D%nblocks_y

    D%C%nx = nx
    D%C%ny = ny

    ! allocate memory for each block and interface

    allocate(D%B(D%nblocks))
    allocate(D%I(D%nifaces))

    ! initialize FD coefficients

    call init_fd(D%FDmethod)

    ! read operator information for each block

    do i = 1,D%nblocks

       c_dissipation = 0d0 ! default

       write(str,'(a,i0,a)') '!---BLOCK',i,'---'
       call seek_to_string(input,str)
       read(input,nml=operator_list,iostat=stat)
       if (stat>0) call error('Error in operator_list','init_domain')

       D%B(i)%c_dissipation = c_dissipation
       D%B(i)%dissipation = (c_dissipation>0d0)

       if (is_master) then
          write(str,'(a,i0,a)') 'B{',i,'}'
          call write_matlab(echo,'c_dissipation',D%B(i)%c_dissipation,str)
       end if

    end do

    ! read grid information for each block;
    ! determine starting indices of each block (used in domain decomposition routine,
    ! to prevent process from having only one point in block and such things)

    allocate(blkxm(D%nblocks),blkym(D%nblocks))

    do i = 1,D%nblocks
       call read_grid(i,D%B(i)%G,refine,input,echo)
       call block_limits(D%B(i)%G,blkxm(i),blkym(i))
    end do

    ! MPI decomposition of global 2D domain

    select case(mode)
    case(2)
       nF = 6
    case(3)
       nF = 3
    end select

    call decompose2d(D%C,nF,periodic_x,periodic_y, &
         mpi_method,nprocs_x,nprocs_y,blkxm,blkym)

    deallocate(blkxm,blkym)

    ! decomposition information

    if (decomposition_info) then
       write(str1,'(i12)') D%C%c2d%coord(1)
       write(str2,'(i12)') D%C%c2d%coord(2)
       write(str3,'(a)') 'coord = (' // &
            trim(adjustl(str1)) // ',' // trim(adjustl(str2)) // ')'
       write(str1,'(a,i6,a,i6,a,i6,a,i6)') &
            'x:',D%C%mx,'<=i<=',D%C%px,' lnx=',D%C%lnx,' nprocs_x=',D%C%c2d%nprocs_x
       write(str2,'(a,i6,a,i6,a,i6,a,i6)') &
            'y:',D%C%my,'<=j<=',D%C%py,' lny=',D%C%lny,' nprocs_y=',D%C%c2d%nprocs_y
       call messages(str3,str1,str2)
    end if

    ! initialize grid blocks and build global grid

    do i = 1,D%nblocks
       call init_grid(i,D%B(i)%G,D%G,D%C,D%FDmethod,D%exact_metric)
    end do

    call exchange_all_neighbors(D%C,D%G%x)
    call exchange_all_neighbors(D%C,D%G%y)

    do i = 1,D%nblocks
       call init_grid_partials(D%B(i)%G,D%G,D%exact_metric)
    end do

    call grid_spacing(D%C,D%G)

    call exchange_all_neighbors(D%C,D%G%xq)
    call exchange_all_neighbors(D%C,D%G%xr)
    call exchange_all_neighbors(D%C,D%G%yq)
    call exchange_all_neighbors(D%C,D%G%yr)
    call exchange_all_neighbors(D%C,D%G%J )

    ! initialize interfaces (grid and processor information)

    do i = 1,D%nifaces
       call init_iface(i,D%I(i),refine,energy_balance,input)
       im = D%I(i)%iblockm
       ip = D%I(i)%iblockp
       call init_iface_blocks(i,D%I(i),D%B(im)%G,D%B(ip)%G,D%C, &
            refine,input,echo)
       !call exchange_grid_edges(D%I(i),D%C,D%B(im)%G,D%B(ip)%G)
    end do

    ! initialize blocks (material properties, initial fields, etc.)

    do i = 1,D%nblocks
       call init_material(i,D%B(i)%M,input,echo)
       call init_pml(i,D%B(i)%M,input,echo)
       call init_boundaries(i,D%B(i)%B,input,echo)
       call init_fields(D%mode,i,D%C,D%B(i)%G,D%G,D%B(i)%F,D%F,D%t, &
            input,D%B(i)%M%response,energy_balance,displacement, &
            D%B(i)%M,D%B(i)%M%pmlx,D%B(i)%M%pmly)
       call interior_to_edges(D%B(i)%G,D%F,D%B(i)%F)
       call apply_bc(D%B(i)%G,D%B(i)%F,D%B(i)%B,D%B(i)%M,D%mode,D%t,i)
    end do

    call exchange_fields(D%C,D%F)

    ! gravity parameters

    call init_gravity(D%grav,input,echo)

    ! time step and CFL parameter

    cs = 0d0
    do i = 1,D%nblocks
       cs = max(cs,D%B(i)%M%cs)
    end do

    dt = dt/refine

    if (CFL==0d0) then
       if (dt==0d0) then ! default values
          CFL = 0.5d0
          dt = CFL*D%G%hmin/cs
       else ! set CFL using dt
          CFL = cs*dt/D%G%hmin
       end if
    else
       if (dt==0d0) then ! set dt using CFL
          dt = CFL*D%G%hmin/cs
       else ! both dt and CFL input, default to dt
          call message('Both dt and CFL specified; defaulting to dt','init')
          CFL = cs*dt/D%G%hmin
       end if
    end if
    if (CFL>1d0) then
      write(str,'(a,e,a)') 'Using CFL=',CFL,'>1; numerical instability likely'
      call message(str,'init')
    end if

    ! Runge-Kutta stage lengths

    if (D%operator_split) then
       allocate(dtRK(1))
       dtRK = dt
    else
       allocate(dtRK(size(cRK)-1))
       do i = 1,size(cRK)-1
          dtRK(i) = dt*(cRK(i+1)-cRK(i))
       end do
    end if

    ! initialize interfaces (fields)

    do i = 1,D%nifaces
       im = D%I(i)%iblockm
       ip = D%I(i)%iblockp
       call init_iface_fields(i,D%I(i),refine,input,echo,dtRK)
       call exchange_fields_edges(D%I(i),D%C,D%B(im)%F,D%B(ip)%F)
       ! JK: Maybe not the best place for this!
       if(D%I(i)%bndr_rates) then
         select case(D%I(i)%direction)
         case('x')
           call init_fields_side_rate(D%B(im)%F%bndFR)
           call init_fields_side_rate(D%B(ip)%F%bndFL)
         case('y')
           call init_fields_side_rate(D%B(im)%F%bndFT)
           call init_fields_side_rate(D%B(ip)%F%bndFB)
         end select
       end if
       call couple_blocks(D%I(i),D%B(im)%F,D%B(ip)%F,D%B(im)%M,D%B(ip)%M, &
            D%C,D%mode,D%t,initialize=.true.)
    end do

    deallocate(dtRK)

    ! output domain parameters

    if (is_master) then
       call write_matlab(echo,'mode',D%mode)
       call write_matlab(echo,'FDmethod',D%FDmethod)
       call write_matlab(echo,'method',D%method)
       call write_matlab(echo,'c_dissipation',c_dissipation)
       call write_matlab(echo,'nblocks',D%nblocks)
       call write_matlab(echo,'nifaces',D%nifaces)
       call write_matlab(echo,'nx',D%C%nx)
       call write_matlab(echo,'ny',D%C%ny)
       call write_matlab(echo,'dt',dt)
       call write_matlab(echo,'CFL',CFL)
       call write_matlab(echo,'hmin',D%G%hmin)
       call write_matlab(echo,'hmax',D%G%hmax)
       call write_matlab(echo,'h',0.5d0*(D%G%hmin+D%G%hmax))
       call write_matlab(echo,'operator_split',operator_split)
    end if

  end subroutine init_domain


  subroutine finish_domain(D)

    use fd_coeff, only : destroy_fd
    use grid, only : destroy_grid
    use fields, only : destroy_block_fields,destroy_fields
    use boundaries, only : destroy_iface

    implicit none

    type(domain_type),intent(inout) :: D

    integer :: i

    call destroy_fd

    do i = 1,D%nblocks
       call destroy_grid(D%G,D%B(i)%G)
       call destroy_block_fields(D%B(i)%F)
    end do
    deallocate(D%B)

    do i = 1,D%nifaces
       call destroy_iface(D%I(i))
    end do
    deallocate(D%I)

    call destroy_fields(D%F)

  end subroutine finish_domain


  subroutine check_nucleation(D,minV,slipping)

    use boundaries, only : maxV
    use mpi_routines, only : MPI_REAL_PW
    use mpi

    implicit none

    type(domain_type),intent(in) :: D
    real,intent(in) :: minV
    logical,intent(out) :: slipping

    real :: V,Vmax
    integer :: i,ierr

    V = 0d0

    do i = 1,D%nifaces
       call maxV(D%I(i),V)
    end do

    call MPI_Allreduce(V,Vmax,1,MPI_REAL_PW,MPI_MAX,MPI_COMM_WORLD,ierr)

    slipping = (Vmax>=minV)

  end subroutine check_nucleation


  subroutine write_field(fh,D,field,location,mx,px,sx,my,py,sy,i)

    use io, only : file_distributed,write_file_distributed,error

    implicit none

    type(file_distributed),intent(in) :: fh
    type(domain_type),intent(in) :: D
    character(*),intent(in) :: field,location
    integer,intent(in) :: mx,px,sx,my,py,sy
    integer,intent(in),optional :: i

    select case(location)

    case default

       select case(field)
       case default
          call error('Field ' // trim(field) // ' not defined','get_field')
       case('Etot')
          call write_file_distributed(fh,D%F%Etot)
       case('x')
          call write_file_distributed(fh,D%G%x(mx:px:sx,my:py:sy))
       case('y')
          call write_file_distributed(fh,D%G%y(mx:px:sx,my:py:sy))
       case('xq')
          call write_file_distributed(fh,D%G%xq(mx:px:sx,my:py:sy))
       case('xr')
          call write_file_distributed(fh,D%G%xr(mx:px:sx,my:py:sy))
       case('yq')
          call write_file_distributed(fh,D%G%yq(mx:px:sx,my:py:sy))
       case('yr')
          call write_file_distributed(fh,D%G%yr(mx:px:sx,my:py:sy))
       case('J')
          call write_file_distributed(fh,D%G%J(mx:px:sx,my:py:sy))
       case('ux','uz')
          call write_file_distributed(fh,D%F%U(mx:px:sx,my:py:sy,1))
       case('uy')
          call write_file_distributed(fh,D%F%U(mx:px:sx,my:py:sy,2))
       case('vx','vz')
          call write_file_distributed(fh,D%F%F(mx:px:sx,my:py:sy,1))
       case('vy','sxz')
          call write_file_distributed(fh,D%F%F(mx:px:sx,my:py:sy,2))
       case('sxx','syz')
          call write_file_distributed(fh,D%F%F(mx:px:sx,my:py:sy,3))
       case('sxy')
          call write_file_distributed(fh,D%F%F(mx:px:sx,my:py:sy,4))
       case('syy')
          call write_file_distributed(fh,D%F%F(mx:px:sx,my:py:sy,5))
       case('szz')
          call write_file_distributed(fh,D%F%F(mx:px:sx,my:py:sy,6))
       case('lambda')
          call write_file_distributed(fh,D%F%lambda(mx:px:sx,my:py:sy))
       case('gammap')
          call write_file_distributed(fh,D%F%gammap(mx:px:sx,my:py:sy))
       case('F')
          call write_file_distributed(fh,D%F%F(mx:px:sx,my:py:sy,:))
       end select

    case('ifacex','point_ifacex')

       select case(field)
       case default
          call error('Field ' // trim(field) // ' not defined','get_field')
       case('x')
          call write_file_distributed(fh,D%B(D%I(i)%iblockm)%G%bndR%x(my:py:sy))
       case('y')
          call write_file_distributed(fh,D%B(D%I(i)%iblockm)%G%bndR%y(my:py:sy))
       case('V')
          call write_file_distributed(fh,D%I(i)%FR%V(my:py:sy))
       case('O')
          call write_file_distributed(fh,D%I(i)%FR%O(my:py:sy))
       case('S')
          call write_file_distributed(fh,D%I(i)%FR%S(my:py:sy))
       case('N')
          call write_file_distributed(fh,D%I(i)%FR%N(my:py:sy))
       case('a')
          call write_file_distributed(fh,D%I(i)%FR%rs%a(my:py:sy))
       case('b')
          call write_file_distributed(fh,D%I(i)%FR%rs%b(my:py:sy))
       case('L')
          call write_file_distributed(fh,D%I(i)%FR%rs%L(my:py:sy))
       case('V0')
          call write_file_distributed(fh,D%I(i)%FR%rs%V0(my:py:sy))
       case('Vw')
          call write_file_distributed(fh,D%I(i)%FR%rs%Vw(my:py:sy))
       case('f0')
          call write_file_distributed(fh,D%I(i)%FR%rs%f0(my:py:sy))
       case('fw')
          call write_file_distributed(fh,D%I(i)%FR%rs%fw(my:py:sy))
       case('S0')
          call write_file_distributed(fh,D%I(i)%FR%S0(my:py:sy))
       case('N0')
          call write_file_distributed(fh,D%I(i)%FR%N0(my:py:sy))
       case('Ds')
          call write_file_distributed(fh,D%I(i)%FR%Ds(my:py:sy))
       case('Dn')
          call write_file_distributed(fh,D%I(i)%FR%Dn(my:py:sy))
       case('E')
          call write_file_distributed(fh,D%I(i)%E(my:py:sy))
       case('D')
          call write_file_distributed(fh,D%I(i)%FR%D(my:py:sy))
       case('Psi')
          call write_file_distributed(fh,D%I(i)%FR%Psi(my:py:sy))
       case('trup')
          call write_file_distributed(fh,D%I(i)%FR%trup(my:py:sy))
       case('z')
          call write_file_distributed(fh,D%I(i)%TP%z(mx:px:sx))
       case('T')
          call write_file_distributed(fh,D%I(i)%TP%T(mx:px:sx,my:py:sy))
       case('P')
          call write_file_distributed(fh,D%I(i)%TP%p(mx:px:sx,my:py:sy))
       case('p')
          call write_file_distributed(fh,D%I(i)%ER%p(my:py:sy))
       case('rho')
          call write_file_distributed(fh,D%I(i)%ER%rho(my:py:sy))
       case('u')
          call write_file_distributed(fh,D%I(i)%ER%u(my:py:sy))
       case('c')
          call write_file_distributed(fh,D%I(i)%ER%c(my:py:sy))
       case('w')
          call write_file_distributed(fh,D%I(i)%ER%w(my:py:sy))
       case('Smin')
          call write_file_distributed(fh,D%I(i)%ER%Smin(my:py:sy))
       case('n')
          call write_file_distributed(fh,D%I(i)%ER%n(my:py:sy))
       case('q1')
          call write_file_distributed(fh,D%I(i)%ER%q(my:py:sy,1))
       case('q2')
          call write_file_distributed(fh,D%I(i)%ER%q(my:py:sy,2))
       end select

    case('ifacey','point_ifacey')

       select case(field)
       case default
          call error('Field ' // trim(field) // ' not defined','get_field')
       case('x')
          call write_file_distributed(fh,D%B(D%I(i)%iblockm)%G%bndT%x(mx:px:sx))
       case('y')
          call write_file_distributed(fh,D%B(D%I(i)%iblockm)%G%bndT%y(mx:px:sx))
       case('V')
          call write_file_distributed(fh,D%I(i)%FR%V(mx:px:sx))
       case('O')
          call write_file_distributed(fh,D%I(i)%FR%O(mx:px:sx))
       case('S')
          call write_file_distributed(fh,D%I(i)%FR%S(mx:px:sx))
       case('N')
          call write_file_distributed(fh,D%I(i)%FR%N(mx:px:sx))
       case('a')
          call write_file_distributed(fh,D%I(i)%FR%rs%a(mx:px:sx))
       case('b')
          call write_file_distributed(fh,D%I(i)%FR%rs%b(mx:px:sx))
       case('L')
          call write_file_distributed(fh,D%I(i)%FR%rs%L(mx:px:sx))
       case('V0')
          call write_file_distributed(fh,D%I(i)%FR%rs%V0(mx:px:sx))
       case('Vw')
          call write_file_distributed(fh,D%I(i)%FR%rs%Vw(mx:px:sx))
       case('f0')
          call write_file_distributed(fh,D%I(i)%FR%rs%f0(mx:px:sx))
       case('fw')
          call write_file_distributed(fh,D%I(i)%FR%rs%fw(mx:px:sx))
       case('S0')
          call write_file_distributed(fh,D%I(i)%FR%S0(mx:px:sx))
       case('N0')
          call write_file_distributed(fh,D%I(i)%FR%N0(mx:px:sx))
       case('Ds')
          call write_file_distributed(fh,D%I(i)%FR%Ds(mx:px:sx))
       case('Dn')
          call write_file_distributed(fh,D%I(i)%FR%Dn(mx:px:sx))
       case('E')
          call write_file_distributed(fh,D%I(i)%E(mx:px:sx))
       case('D')
          call write_file_distributed(fh,D%I(i)%FR%D(mx:px:sx))
       case('Psi')
          call write_file_distributed(fh,D%I(i)%FR%Psi(mx:px:sx))
       case('trup')
          call write_file_distributed(fh,D%I(i)%FR%trup(mx:px:sx))
       case('z')
          call write_file_distributed(fh,D%I(i)%TP%z(my:py:sy))
       case('T')
          call write_file_distributed(fh,transpose(D%I(i)%TP%T(my:py:sy,mx:px:sx)))
       case('P')
          call write_file_distributed(fh,transpose(D%I(i)%TP%p(my:py:sy,mx:px:sx)))
       case('p')
          call write_file_distributed(fh,D%I(i)%ER%p(mx:px:sx))
       case('rho')
          call write_file_distributed(fh,D%I(i)%ER%rho(mx:px:sx))
       case('u')
          call write_file_distributed(fh,D%I(i)%ER%u(mx:px:sx))
       case('c')
          call write_file_distributed(fh,D%I(i)%ER%c(mx:px:sx))
       case('w')
          call write_file_distributed(fh,D%I(i)%ER%w(mx:px:sx))
       case('Smin')
          call write_file_distributed(fh,D%I(i)%ER%Smin(mx:px:sx))
       case('n')
          call write_file_distributed(fh,D%I(i)%ER%n(mx:px:sx))
       case('q1')
          call write_file_distributed(fh,D%I(i)%ER%q(mx:px:sx,1))
       case('q2')
          call write_file_distributed(fh,D%I(i)%ER%q(mx:px:sx,2))
       end select

    case('Eblock')

       select case(field)
       case default
          call error('Field ' // trim(field) // ' not defined','get_field')
       case('Etot')
          call write_file_distributed(fh,D%B(i)%F%Etot)
       case('E')
          call write_file_distributed(fh,sum(D%B(i)%F%E(:)))
       case('EL')
          call write_file_distributed(fh,D%B(i)%F%E(1))
       case('ER')
          call write_file_distributed(fh,D%B(i)%F%E(2))
       case('EB')
          call write_file_distributed(fh,D%B(i)%F%E(3))
       case('ET')
          call write_file_distributed(fh,D%B(i)%F%E(4))
       case('DE')
          call write_file_distributed(fh,sum(D%B(i)%F%DE(:)))
       case('DEL')
          call write_file_distributed(fh,D%B(i)%F%DE(1))
       case('DER')
          call write_file_distributed(fh,D%B(i)%F%DE(2))
       case('DEB')
          call write_file_distributed(fh,D%B(i)%F%DE(3))
       case('DET')
          call write_file_distributed(fh,D%B(i)%F%DE(4))
       end select

    end select

  end subroutine write_field


  subroutine check_field(D,field,location,i)

    use io, only : error
    use utilities, only : within

    implicit none

    type(domain_type),intent(in) :: D
    character(*),intent(in) :: field,location
    integer,intent(in),optional :: i

    logical :: ok

    select case(location)

    case default

       select case(field)
       case default
          ok = .false.
       case('Etot')
          ok = .true.
       case('x')
          ok = allocated(D%G%x)
       case('y')
          ok = allocated(D%G%y)
       case('xq')
          ok = allocated(D%G%xq)
       case('xr')
          ok = allocated(D%G%xr)
       case('yq')
          ok = allocated(D%G%yq)
       case('yr')
          ok = allocated(D%G%yr)
       case('J')
          ok = allocated(D%G%J)
       case('ux','uz')
          ok = allocated(D%F%U)
          if (ok) ok = within(lbound(D%F%U,3),1,ubound(D%F%U,3))
       case('uy')
          ok = allocated(D%F%U)
          if (ok) ok = within(lbound(D%F%U,3),2,ubound(D%F%U,3))
       case('vx','vz')
          ok = allocated(D%F%F)
          if (ok) ok = within(lbound(D%F%F,3),1,ubound(D%F%F,3))
       case('vy','sxz')
          ok = allocated(D%F%F)
          if (ok) ok = within(lbound(D%F%F,3),2,ubound(D%F%F,3))
       case('sxx','syz')
          ok = allocated(D%F%F)
          if (ok) ok = within(lbound(D%F%F,3),3,ubound(D%F%F,3))
       case('sxy')
          ok = allocated(D%F%F)
          if (ok) ok = within(lbound(D%F%F,3),4,ubound(D%F%F,3))
       case('syy')
          ok = allocated(D%F%F)
          if (ok) ok = within(lbound(D%F%F,3),5,ubound(D%F%F,3))
       case('szz')
          ok = allocated(D%F%F)
          if (ok) ok = within(lbound(D%F%F,3),6,ubound(D%F%F,3))
       case('lambda')
          ok = allocated(D%F%lambda)
       case('gammap')
          ok = allocated(D%F%gammap)
       case('F')
          ok = allocated(D%F%F)
       end select

    case('ifacex','point_ifacex','ifacey','point_ifacey')

       select case(field)
       case default
          ok = .false.
       case('x')
          if (index(location,'ifacex')/=0) ok = allocated(D%B(D%I(i)%iblockm)%G%bndR%x)
          if (index(location,'ifacey')/=0) ok = allocated(D%B(D%I(i)%iblockm)%G%bndT%x)
       case('y')
          if (index(location,'ifacex')/=0) ok = allocated(D%B(D%I(i)%iblockm)%G%bndR%y)
          if (index(location,'ifacey')/=0) ok = allocated(D%B(D%I(i)%iblockm)%G%bndT%y)
       case('V')
          ok = allocated(D%I(i)%FR%V)
       case('O')
          ok = allocated(D%I(i)%FR%O)
       case('S')
          ok = allocated(D%I(i)%FR%S)
       case('N')
          ok = allocated(D%I(i)%FR%N)
       case('a')
          ok = allocated(D%I(i)%FR%rs%a)
       case('b')
          ok = allocated(D%I(i)%FR%rs%b)
       case('f0')
          ok = allocated(D%I(i)%FR%rs%f0)
       case('L')
          ok = allocated(D%I(i)%FR%rs%L)
       case('V0')
          ok = allocated(D%I(i)%FR%rs%V0)
       case('fw')
          ok = allocated(D%I(i)%FR%rs%fw)
       case('Vw')
          ok = allocated(D%I(i)%FR%rs%Vw)
       case('S0')
          ok = allocated(D%I(i)%FR%S0)
       case('N0')
          ok = allocated(D%I(i)%FR%N0)
       case('Ds')
          ok = allocated(D%I(i)%FR%Ds)
       case('Dn')
          ok = allocated(D%I(i)%FR%Dn)
       case('E')
          ok = allocated(D%I(i)%E)
       case('D')
          ok = allocated(D%I(i)%FR%D)
       case('Psi')
          ok = allocated(D%I(i)%FR%Psi)
       case('trup')
          ok = allocated(D%I(i)%FR%trup)
       case('z')
          ok = allocated(D%I(i)%TP%z)
       case('T')
          ok = allocated(D%I(i)%TP%T)
       case('P')
          ok = allocated(D%I(i)%TP%p)
       case('p')
          ok = allocated(D%I(i)%ER%p)
       case('rho')
          ok = allocated(D%I(i)%ER%rho)
       case('u')
          ok = allocated(D%I(i)%ER%u)
       case('c')
          ok = allocated(D%I(i)%ER%c)
       case('w')
          ok = allocated(D%I(i)%ER%w)
       case('Smin')
          ok = allocated(D%I(i)%ER%Smin)
       case('n')
          ok = allocated(D%I(i)%ER%n)
       case('q1','q2')
          ok = allocated(D%I(i)%ER%q)
       end select

    case('Eblock')

       select case(field)
       case default
          ok = .false.
       case('Etot','E','EL','ER','EB','ET','DE','DEL','DER','DEB','DET')
          ok = .true.
       end select

    end select

    if (.not.ok) call error('Problem with output field ' // trim(field) // &
         ' at location ' // trim(location),'check_field')

  end subroutine check_field


  subroutine exchange_rates_field_edge(I,C,Fm,Fp)

    use fields, only : block_fields
    use mpi_routines2d,only : cartesian,exchange_edge

    implicit none

    type(iface_type),intent(in) :: I
    type(cartesian),intent(in) :: C
    type(block_fields),intent(inout) :: Fm,Fp

    integer :: l,nF

    if (.not.I%share) return ! no need to share fields

    ! initialize arrays in blocks not handled by process
    ! with values from neighboring process

    select case(I%direction)
    case('x')
       if (.not.allocated(Fm%bndFR%DF)) return
       nF = size(Fm%bndFR%DF,2)
       do l = 1,nF
          call exchange_edge(C,I%m,I%p,Fm%bndFR%DF(I%m:I%p,l),Fp%bndFL%DF(I%m:I%p,l), &
               I%rank_m,I%rank_p,'x')
       end do
    case('y')
       if (.not.allocated(Fm%bndFT%DF)) return
       nF = size(Fm%bndFT%DF,2)
       do l = 1,nF
          call exchange_edge(C,I%m,I%p,Fm%bndFT%DF(I%m:I%p,l),Fp%bndFB%DF(I%m:I%p,l), &
               I%rank_m,I%rank_p,'y')
       end do
    end select

  end subroutine exchange_rates_field_edge

  subroutine exchange_fields_edges(I,C,Fm,Fp)

    use fields, only : block_fields
    use mpi_routines2d,only : cartesian,exchange_edge

    implicit none

    type(iface_type),intent(in) :: I
    type(cartesian),intent(in) :: C
    type(block_fields),intent(inout) :: Fm,Fp

    integer :: l,nF

    if (.not.I%share) return ! no need to share fields

    ! initialize arrays in blocks not handled by process
    ! with values from neighboring process

    select case(I%direction)
    case('x')
       nF = size(Fm%bndFR%F,2)
       do l = 1,nF
          call exchange_edge(C,I%m,I%p,Fm%bndFR%F(I%m:I%p,l),Fp%bndFL%F(I%m:I%p,l), &
               I%rank_m,I%rank_p,'x')
       end do
    case('y')
       nF = size(Fm%bndFT%F,2)
       do l = 1,nF
          call exchange_edge(C,I%m,I%p,Fm%bndFT%F(I%m:I%p,l),Fp%bndFB%F(I%m:I%p,l), &
               I%rank_m,I%rank_p,'y')
       end do
    end select

  end subroutine exchange_fields_edges


  subroutine exchange_grid_edges(I,C,Gm,Gp)

    use grid, only : block_grid
    use mpi_routines2d,only : cartesian,exchange_edge

    implicit none

    type(iface_type),intent(in) :: I
    type(cartesian),intent(in) :: C
    type(block_grid),intent(inout) :: Gm,Gp

    if (.not.I%share) return ! no need to share fields

    ! initialize arrays in blocks not handled by process
    ! with values from neighboring process

    select case(I%direction)
    case('x')
       call exchange_edge(C,I%m,I%p,Gm%bndR%x(I%m:I%p),Gp%bndL%x(I%m:I%p), &
            I%rank_m,I%rank_p,'x')
       call exchange_edge(C,I%m,I%p,Gm%bndR%y(I%m:I%p),Gp%bndL%y(I%m:I%p), &
            I%rank_m,I%rank_p,'x')
       call exchange_edge(C,I%m,I%p,Gm%bndR%n(I%m:I%p,1),Gp%bndL%n(I%m:I%p,1), &
            I%rank_m,I%rank_p,'x')
       call exchange_edge(C,I%m,I%p,Gm%bndR%n(I%m:I%p,2),Gp%bndL%n(I%m:I%p,2), &
            I%rank_m,I%rank_p,'x')
    case('y')
       call exchange_edge(C,I%m,I%p,Gm%bndT%x(I%m:I%p),Gp%bndB%x(I%m:I%p), &
            I%rank_m,I%rank_p,'y')
       call exchange_edge(C,I%m,I%p,Gm%bndT%y(I%m:I%p),Gp%bndB%y(I%m:I%p), &
            I%rank_m,I%rank_p,'y')
       call exchange_edge(C,I%m,I%p,Gm%bndT%n(I%m:I%p,1),Gp%bndB%n(I%m:I%p,1), &
            I%rank_m,I%rank_p,'y')
       call exchange_edge(C,I%m,I%p,Gm%bndT%n(I%m:I%p,2),Gp%bndB%n(I%m:I%p,2), &
            I%rank_m,I%rank_p,'y')
    end select

  end subroutine exchange_grid_edges


  subroutine interior_to_edges(B,F,BF)

    use grid, only : block_grid
    use fields, only : fields_type,block_fields
    use io, only : error

    implicit none

    type(block_grid),intent(in) :: B
    type(fields_type),intent(in) :: F
    type(block_fields),intent(inout) :: BF

    ! reconstruct fields on edges of block, which are later
    ! used when setting boundary and interface conditions

    if (B%skip) return ! process has no cells in this block

    if (B%sideL.and.B%nx/=1) BF%bndFL%F = F%F(B%mgx,B%my:B%py,:)
    if (B%sideR.and.B%nx/=1) BF%bndFR%F = F%F(B%pgx,B%my:B%py,:)
    if (B%sideB.and.B%ny/=1) BF%bndFB%F = F%F(B%mx:B%px,B%mgy,:)
    if (B%sideT.and.B%ny/=1) BF%bndFT%F = F%F(B%mx:B%px,B%pgy,:)

  end subroutine interior_to_edges

  subroutine interior_rates_to_edges(B,F,BF,A)

    use grid, only : block_grid
    use fields, only : fields_type,block_fields
    use io, only : error
    use mpi_routines, only : myid

    implicit none

    type(block_grid),intent(in) :: B
    type(fields_type),intent(in) :: F
    type(block_fields),intent(inout) :: BF
    real,intent(in) :: A

    ! if (B%skip) return ! process has no cells in this block

    ! Store last scaled rate and copy over new (unscaled) rate
    if (allocated(BF%bndFL%DF)) then
      BF%bndFL%DF2 = A*BF%bndFL%DF
      if(B%sideL.and.B%nx/=1) BF%bndFL%DF  = F%DF(B%mgx,B%my:B%py,:)
    end if
    if (allocated(BF%bndFR%DF)) then
      BF%bndFR%DF2 = A*BF%bndFR%DF
      if(B%sideR.and.B%nx/=1) BF%bndFR%DF  = F%DF(B%pgx,B%my:B%py,:)
    end if
    if (allocated(BF%bndFB%DF)) then
      BF%bndFB%DF2 = A*BF%bndFB%DF
      if(B%sideB.and.B%ny/=1) BF%bndFB%DF  = F%DF(B%mx:B%px,B%mgy,:)
    end if
    if (allocated(BF%bndFT%DF)) then
      BF%bndFT%DF2 = A*BF%bndFT%DF
      if(B%sideT.and.B%ny/=1) BF%bndFT%DF  = F%DF(B%mx:B%px,B%pgy,:)
    end if

  end subroutine interior_rates_to_edges


end module domain
