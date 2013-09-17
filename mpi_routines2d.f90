module mpi_routines2d

  implicit none

  type :: cartesian2d
     integer :: comm,myid,nprocs_x,nprocs_y, &
          coord(2),rank_mx,rank_px,rank_my,rank_py, &
          array_w,array_s,arrayF_w,arrayF_s
  end type cartesian2d
  
  type :: cartesian1d
     integer :: comm,myid,nprocs,coord(1),rank_m,rank_p, &
          array_w,array_s
  end type cartesian1d

  type :: cartesian
     type(cartesian2d) :: c2d
     type(cartesian1d) :: c1dx,c1dy
     integer :: nb=3, &
          nx,mx,px,mbx,pbx,lnx, &
          ny,my,py,mby,pby,lny
     integer :: line_x,line_y,block_x,block_y
  end type cartesian

  interface allocate_array_body
     module procedure allocate_array_body_2d,allocate_array_body_3d
  end interface
  
contains


  subroutine exchange_all_neighbors(C,F)

    implicit none

    type(cartesian),intent(in) :: C
    real,dimension(C%mbx:C%pbx,C%mby:C%pby),intent(inout) :: F

    call exchange_neighbors(C,F,'xm')
    call exchange_neighbors(C,F,'xp')
    call exchange_neighbors(C,F,'ym')
    call exchange_neighbors(C,F,'yp')

  end subroutine exchange_all_neighbors
  

  subroutine exchange_neighbors(C,F,side,method_in)

    use mpi_routines, only : MPI_REAL_PW
    use mpi

    implicit none

    type(cartesian),intent(in) :: C
    real,dimension(C%mbx:C%pbx,C%mby:C%pby),intent(inout) :: F
    character(*),intent(in) :: side
    integer,intent(in),optional :: method_in

    integer :: ierr,method
    integer,parameter :: tagxm=1, tagxp=2, tagym=3, tagyp=4
    real,dimension(:,:),allocatable :: S,R ! contiguous send/receive buffers

    ! default method

    if (present(method_in)) then
       method = method_in
    else
       method = 1
    end if

    ! share data between processes

    select case(method)
       
    case(1) ! blocking sendrecv with 1D communicators
       
       select case(side)
       case('xm') ! px --> mbx
          call MPI_SendRecv( &
               F(C%px-C%nb+1,C%my),1,C%block_x,C%c1dx%rank_p,tagxm, &
               F(C%mx-C%nb  ,C%my),1,C%block_x,C%c1dx%rank_m,tagxm, &
               C%c1dx%comm,MPI_STATUS_IGNORE,ierr)
       case('xp') ! pbx <-- mx
          call MPI_SendRecv( &
               F(C%mx  ,C%my),1,C%block_x,C%c1dx%rank_m,tagxp, &
               F(C%px+1,C%my),1,C%block_x,C%c1dx%rank_p,tagxp, &
               C%c1dx%comm,MPI_STATUS_IGNORE,ierr)
       case('ym') ! py --> mby
          call MPI_SendRecv( &
               F(C%mx,C%py-C%nb+1),1,C%block_y,C%c1dy%rank_p,tagym, &
               F(C%mx,C%my-C%nb  ),1,C%block_y,C%c1dy%rank_m,tagym, &
               C%c1dy%comm,MPI_STATUS_IGNORE,ierr)
       case('yp') ! pby <-- my
          call MPI_SendRecv( &
               F(C%mx,C%my  ),1,C%block_y,C%c1dy%rank_m,tagyp, &
               F(C%mx,C%py+1),1,C%block_y,C%c1dy%rank_p,tagyp, &
               C%c1dy%comm,MPI_STATUS_IGNORE,ierr)
       end select

    case(2) ! blocking sendrecv with 2D communicator
       
       select case(side)
       case('xm') ! px --> mbx
          call MPI_SendRecv( &
               F(C%px-C%nb+1,C%my),1,C%block_x,C%c2d%rank_px,tagxm, &
               F(C%mx-C%nb  ,C%my),1,C%block_x,C%c2d%rank_mx,tagxm, &
               C%c2d%comm,MPI_STATUS_IGNORE,ierr)
       case('xp') ! pbx <-- mx
          call MPI_SendRecv( &
               F(C%mx  ,C%my),1,C%block_x,C%c2d%rank_mx,tagxp, &
               F(C%px+1,C%my),1,C%block_x,C%c2d%rank_px,tagxp, &
               C%c2d%comm,MPI_STATUS_IGNORE,ierr)
       case('ym') ! py --> mby
          call MPI_SendRecv( &
               F(C%mx,C%py-C%nb+1),1,C%block_y,C%c2d%rank_py,tagym, &
               F(C%mx,C%my-C%nb  ),1,C%block_y,C%c2d%rank_my,tagym, &
               C%c2d%comm,MPI_STATUS_IGNORE,ierr)
       case('yp') ! pby <-- my
          call MPI_SendRecv( &
               F(C%mx,C%my  ),1,C%block_y,C%c2d%rank_my,tagyp, &
               F(C%mx,C%py+1),1,C%block_y,C%c2d%rank_py,tagyp, &
               C%c2d%comm,MPI_STATUS_IGNORE,ierr)
       end select

    case(3) ! blocking sendrecv with 1D communicators, pack into contiguous array
       
       select case(side)
       case('xm') ! px --> mbx
          allocate(S(C%nb,C%lny),R(C%nb,C%lny))
          S = F(C%px-C%nb+1:C%px,C%my:C%py)
          call MPI_SendRecv( &
               S(1,1),C%nb*C%lny,MPI_REAL_PW,C%c1dx%rank_p,tagxm, &
               R(1,1),C%nb*C%lny,MPI_REAL_PW,C%c1dx%rank_m,tagxm, &
               C%c1dx%comm,MPI_STATUS_IGNORE,ierr)
          F(C%mx-C%nb:C%mx-1,C%my:C%py) = R
       case('xp') ! pbx <-- mx
          allocate(S(C%nb,C%lny),R(C%nb,C%lny))
          S = F(C%mx:C%mx+C%nb-1,C%my:C%py)
          call MPI_SendRecv( &
               S(1,1),C%nb*C%lny,MPI_REAL_PW,C%c1dx%rank_m,tagxp, &
               R(1,1),C%nb*C%lny,MPI_REAL_PW,C%c1dx%rank_p,tagxp, &
               C%c1dx%comm,MPI_STATUS_IGNORE,ierr)
          F(C%px+1:C%px+C%nb,C%my:C%py) = R
       case('ym') ! py --> mby
          allocate(S(C%lnx,C%nb),R(C%lnx,C%nb))
          S = F(C%mx:C%px,C%py-C%nb+1:C%py)
          call MPI_SendRecv( &
               S(1,1),C%nb*C%lnx,MPI_REAL_PW,C%c1dy%rank_p,tagym, &
               R(1,1),C%nb*C%lnx,MPI_REAL_PW,C%c1dy%rank_m,tagym, &
               C%c1dy%comm,MPI_STATUS_IGNORE,ierr)
          F(C%mx:C%px,C%my-C%nb:C%my-1) = R
       case('yp') ! pby <-- my
          allocate(S(C%lnx,C%nb),R(C%lnx,C%nb))
          S = F(C%mx:C%px,C%my:C%my+C%nb-1)
          call MPI_SendRecv( &
               S(1,1),C%nb*C%lnx,MPI_REAL_PW,C%c1dy%rank_m,tagyp, &
               R(1,1),C%nb*C%lnx,MPI_REAL_PW,C%c1dy%rank_p,tagyp, &
               C%c1dy%comm,MPI_STATUS_IGNORE,ierr)
          F(C%mx:C%px,C%py+1:C%py+C%nb) = R
       end select
       deallocate(S,R)

    end select

  end subroutine exchange_neighbors
  

  subroutine exchange_edge(C,m,p,Fm,Fp,rank_m,rank_p,direction)

    use mpi_routines, only : MPI_REAL_PW
    use mpi

    implicit none

    type(cartesian),intent(in) :: C
    integer,intent(in) :: m,p,rank_m,rank_p
    real,dimension(m:p),intent(inout) :: Fm,Fp
    character(*),intent(in) :: direction

    integer :: n,myid,comm,ierr
    integer,parameter :: tagm=5, tagp=6

    n = p-m+1

    select case(direction)
    case('x')
       myid = C%c1dx%myid
       comm = C%c1dx%comm
    case('y')
       myid = C%c1dy%myid
       comm = C%c1dy%comm
    end select

    ! send m, receive p
    
    if (myid==rank_m) call MPI_SendRecv( &
         Fm(m),n,MPI_REAL_PW,rank_p,tagm, &
         Fp(m),n,MPI_REAL_PW,rank_p,tagp, &
         comm,MPI_STATUS_IGNORE,ierr)
       
    ! send p, receive m
    
    if (myid==rank_p) call MPI_SendRecv( &
         Fp(m),n,MPI_REAL_PW,rank_m,tagp, &
         Fm(m),n,MPI_REAL_PW,rank_m,tagm, &
         comm,MPI_STATUS_IGNORE,ierr)

  end subroutine exchange_edge


  subroutine populate_ghost_cells(C,m,p,nb,endm,endp,F,direction)

    use mpi_routines, only : MPI_REAL_PW
    use mpi

    implicit none

    type(cartesian),intent(in) :: C
    integer,intent(in) :: m,p,nb
    logical,intent(in) :: endm,endp
    real,dimension(m-nb:p+nb),intent(inout) :: F
    character(*),intent(in) :: direction

    integer :: comm,rank_m,rank_p,ierr
    integer,parameter :: tagm=7, tagp=8

    select case(direction)
    case('x')
       comm = C%c1dx%comm
       rank_m = C%c1dx%rank_m
       rank_p = C%c1dx%rank_p
    case('y')
       comm = C%c1dy%comm
       rank_m = C%c1dy%rank_m
       rank_p = C%c1dy%rank_p
    end select

    if (endm) rank_m = MPI_PROC_NULL
    if (endp) rank_p = MPI_PROC_NULL

    ! p --> m
    
    call MPI_SendRecv( &
         F(p-nb+1),nb,MPI_REAL_PW,rank_p,tagm, &
         F(m-nb  ),nb,MPI_REAL_PW,rank_m,tagm, &
         comm,MPI_STATUS_IGNORE,ierr)

    ! p <-- m
    
    call MPI_SendRecv( &
         F(m     ),nb,MPI_REAL_PW,rank_m,tagp, &
         F(p+1   ),nb,MPI_REAL_PW,rank_p,tagp, &
         comm,MPI_STATUS_IGNORE,ierr)

  end subroutine populate_ghost_cells


  subroutine share_1d_array(C,source,array)

    use mpi_routines, only : MPI_REAL_PW
    use mpi

    implicit none

    type(cartesian1d),intent(in) :: C
    logical,intent(in) :: source
    real,dimension(:),intent(inout) :: array

    integer :: rank_source,ierr
    integer :: rank_source_send ! MPI-1 only

    ! share array on process having source=T with all other processes
    ! in x or y direction (i.e., all processes in 1D communicator C)

    if (C%nprocs==1) return ! no other process in communicator

    ! all processes determine rank of source process

    rank_source = 0
    if (source) rank_source = C%myid

    ! MPI-2 version
    !call MPI_Allreduce(MPI_IN_PLACE,rank_source,1,MPI_INTEGER,MPI_MAX,C%comm,ierr)

    ! MPI-1 version
    rank_source_send = rank_source
    call MPI_Allreduce(rank_source_send,rank_source,1,MPI_INTEGER,MPI_MAX,C%comm,ierr)

    ! broadcast array from source

    call MPI_Bcast(array(1),size(array),MPI_REAL_PW,rank_source,C%comm,ierr)

  end subroutine share_1d_array
  

  subroutine test_exchange(C)

    implicit none

    type(cartesian),intent(in) :: C

    integer :: i,j
    real :: G
    real,dimension(C%mbx:C%pbx,C%mby:C%pby) :: F

    ! initialize array

    F = -1d0

    do j = C%my,C%py
       do i = C%mx,C%px
          F(i,j) = real(i)+1000d0*real(j)
       end do
    end do

    ! exchange

    call exchange_all_neighbors(C,F)

    ! check

    ! left
    if (C%mx/=1) then
       do j = C%my,C%py
          do i = C%mbx,C%mx-1
             G = real(i)+1000d0*real(j)
             if (F(i,j)/=G.and.j==4) print *, i,j,G,F(i,j)
          end do
       end do
    end if

    ! right
    if (C%px/=C%nx) then
       do j = C%my,C%py
          do i = C%px+1,C%pbx
             G = real(i)+1000d0*real(j)
             if (F(i,j)/=G.and.j==4) print *, i,j,G,F(i,j)
          end do
       end do
    end if

    ! bottom
    if (C%my/=1) then
       do j = C%mby,C%my-1
          do i = C%mx,C%px
             G = real(i)+1000d0*real(j)
             if (F(i,j)/=G) print *, i,j,G,F(i,j)
          end do
       end do
    end if

    ! top
    if (C%py/=C%ny) then
       do j = C%py+1,C%pby
          do i = C%mx,C%px
             G = real(i)+1000d0*real(j)
             if (F(i,j)/=G) print *, i,j,G,F(i,j)
          end do
       end do
    end if

  end subroutine test_exchange


  subroutine decompose2d(C,nF,periodic_x,periodic_y, &
       method,nprocs_x_in,nprocs_y_in,blkxm,blkym)

    use io, only : error
    use mpi_routines, only : nprocs,MPI_REAL_PW,MPI_REAL_PS,pw,ps,decompose1d,is_master,subarray
    use mpi

    implicit none

    type(cartesian),intent(inout) :: C
    integer,intent(in) :: nF,nprocs_x_in,nprocs_y_in
    logical,intent(in) :: periodic_x,periodic_y
    character(*),intent(inout) :: method
    integer,dimension(:),intent(in),optional :: blkxm,blkym

    integer,parameter :: dim1=1,dim2=2
    integer :: ierr,index,shift,nprocs_xy(2),l
    logical :: periodic(2),reorder,in_comm(2)
    integer,dimension(:),allocatable :: blocklengths,types
    integer(MPI_ADDRESS_KIND),dimension(:),allocatable :: displacements
    integer,dimension(:),allocatable :: displacements_mpi1

    ! subroutine requires that C%nx, C%ny, and C%nb be set on input

    ! processor layout in Cartesian topology

    if (method=='manual') then

       if (nprocs/=nprocs_x_in*nprocs_y_in) then
          if (is_master) &
               call error('Error: Incorrect number of processors for manual decomposition', &
               'decompose2d')
       end if

       nprocs_xy = (/ nprocs_x_in,nprocs_y_in /)

    else

       select case(method)
       case default ! 2D
          nprocs_xy = 0
       case('1Dx') ! 1D, x
          nprocs_xy(1) = 0
          nprocs_xy(2) = 1
       case('1Dy') ! 1D, y
          nprocs_xy(1) = 1
          nprocs_xy(2) = 0
       end select

       if (C%nx==1) nprocs_xy(1) = 1
       if (C%ny==1) nprocs_xy(2) = 1
       
       call MPI_Dims_create(nprocs,dim2,nprocs_xy,ierr)

    end if

    C%c2d%nprocs_x = nprocs_xy(1)
    C%c2d%nprocs_y = nprocs_xy(2)

    ! 2D Cartesian communicator and coordinates

    periodic = (/ periodic_x,periodic_y /)
    reorder = .true.
    call MPI_Cart_create(MPI_COMM_WORLD,dim2,nprocs_xy,periodic,reorder, &
         C%c2d%comm,ierr)
    call MPI_Comm_rank(C%c2d%comm,C%c2d%myid,ierr)
    call MPI_Cart_coords(C%c2d%comm,C%c2d%myid,dim2,C%c2d%coord,ierr)

    ! nearest neighbors in x-direction

    index = 0; shift = 1
    call MPI_Cart_shift(C%c2d%comm,index,shift,C%c2d%rank_mx,C%c2d%rank_px,ierr)

    ! nearest neighbors in y-direction

    index = 1; shift = 1
    call MPI_Cart_shift(C%c2d%comm,index,shift,C%c2d%rank_my,C%c2d%rank_py,ierr)

    ! 1D communicators (subgrids of 2D Cartesian grid)

    ! x-direction communicator
    in_comm  = (/ .true.,.false. /)
    call MPI_Cart_sub(C%c2d%comm,in_comm,C%c1dx%comm,ierr)
    call MPI_Comm_size(C%c1dx%comm,C%c1dx%nprocs,ierr)
    call MPI_Comm_rank(C%c1dx%comm,C%c1dx%myid,ierr)
    call MPI_Cart_coords(C%c1dx%comm,C%c1dx%myid,dim1,C%c1dx%coord,ierr)
    index = 0; shift = 1
    call MPI_Cart_shift(C%c1dx%comm,index,shift,C%c1dx%rank_m,C%c1dx%rank_p,ierr)

    ! y-direction communicator
    in_comm  = (/ .false.,.true. /)
    call MPI_Cart_sub(C%c2d%comm,in_comm,C%c1dy%comm,ierr)
    call MPI_Comm_size(C%c1dy%comm,C%c1dy%nprocs,ierr)
    call MPI_Comm_rank(C%c1dy%comm,C%c1dy%myid,ierr)
    call MPI_Cart_coords(C%c1dy%comm,C%c1dy%myid,dim1,C%c1dy%coord,ierr)
    index = 0; shift = 1
    call MPI_Cart_shift(C%c1dy%comm,index,shift,C%c1dy%rank_m,C%c1dy%rank_p,ierr)

    ! initial data distribution on processors

    if (present(blkxm)) then
       call decompose1d(C%nx,C%c2d%nprocs_x,C%c2d%coord(1),C%mx,C%px,C%lnx,blkxm)
    else
       call decompose1d(C%nx,C%c2d%nprocs_x,C%c2d%coord(1),C%mx,C%px,C%lnx)
    end if
    if (present(blkym)) then
       call decompose1d(C%ny,C%c2d%nprocs_y,C%c2d%coord(2),C%my,C%py,C%lny,blkym)
    else
       call decompose1d(C%ny,C%c2d%nprocs_y,C%c2d%coord(2),C%my,C%py,C%lny)
    end if
    C%mbx = C%mx-C%nb
    C%pbx = C%px+C%nb
    C%mby = C%my-C%nb
    C%pby = C%py+C%nb

    ! MPI types containing lines of constant x and y

    ! variable x, constant y
    call MPI_Type_vector(C%lnx,1,1,MPI_REAL_PW,C%line_x,ierr)
    call MPI_Type_commit(C%line_x,ierr)

    ! variable y, constant x
    call MPI_Type_vector(C%lny,1,C%lnx+2*C%nb,MPI_REAL_PW,C%line_y,ierr)
    call MPI_Type_commit(C%line_y,ierr)

    ! MPI types containing boundary blocks

    allocate(blocklengths(C%nb),displacements(C%nb),types(C%nb))
    allocate(displacements_mpi1(C%nb))
    blocklengths = 1

    !displacements = (/ (l,l=0,C%nb-1) /)
    displacements_mpi1 = (/ (l,l=0,C%nb-1) /)
    types = C%line_y
    ! MPI-2
    !call MPI_Type_create_struct(C%nb,blocklengths,displacements*pw,types,C%block_x,ierr)
    ! MPI-1
    call MPI_Type_struct(C%nb,blocklengths,displacements_mpi1*pw,types,C%block_x,ierr)
    call MPI_Type_commit(C%block_x,ierr)

    displacements = (/ (l,l=0,C%nb-1) /)*(C%lnx+2*C%nb)
    displacements_mpi1 = (/ (l,l=0,C%nb-1) /)*(C%lnx+2*C%nb)
    types = C%line_x
    ! MPI-2
    !call MPI_Type_create_struct(C%nb,blocklengths,displacements*pw,types,C%block_y,ierr)
    ! MPI-1
    call MPI_Type_struct(C%nb,blocklengths,displacements_mpi1*pw,types,C%block_y,ierr)
    call MPI_Type_commit(C%block_y,ierr)

    deallocate(blocklengths,displacements,types)

    ! MPI types used to set file views for data I/O;
    ! set processor-specific offset and extent

    ! 1D array section along x-axis

    call subarray(C%nx,C%mx,C%px,MPI_REAL_PW,C%c1dx%array_w)
    call subarray(C%nx,C%mx,C%px,MPI_REAL_PS,C%c1dx%array_s)

    ! 1D array section along y-axis

    call subarray(C%ny,C%my,C%py,MPI_REAL_PW,C%c1dy%array_w)
    call subarray(C%ny,C%my,C%py,MPI_REAL_PS,C%c1dy%array_s)

    ! 2D array section

    call subarray(C%nx,C%ny,C%mx,C%px,C%my,C%py,MPI_REAL_PW,C%c2d%array_w)
    call subarray(C%nx,C%ny,C%mx,C%px,C%my,C%py,MPI_REAL_PS,C%c2d%array_s)

    ! 3D array section (2D array section with third index for fields)

    call subarray(C%nx,C%ny,nF,C%mx,C%px,C%my,C%py,1,nF,MPI_REAL_PW,C%c2d%arrayF_w)
    call subarray(C%nx,C%ny,nF,C%mx,C%px,C%my,C%py,1,nF,MPI_REAL_PS,C%c2d%arrayF_s)

  end subroutine decompose2d


  subroutine allocate_array_body_2d(F,C,ghost_nodes,Fval)

    implicit none

    real,dimension(:,:),allocatable,intent(inout) :: F
    type(cartesian),intent(in) :: C
    logical,intent(in) :: ghost_nodes
    real,intent(in),optional :: Fval

    if (.not.allocated(F)) then
       if (ghost_nodes) then
          allocate(F(C%mbx:C%pbx,C%mby:C%pby))
       else
          allocate(F(C%mx :C%px ,C%my :C%py ))
       end if
       if (present(Fval)) then
          F = Fval
       else
          F = 1d40
       end if
    end if

  end subroutine allocate_array_body_2d


  subroutine allocate_array_body_3d(F,C,nF,ghost_nodes,Fval)

    implicit none

    real,dimension(:,:,:),allocatable,intent(inout) :: F
    type(cartesian),intent(in) :: C
    integer,intent(in) :: nF
    logical,intent(in) :: ghost_nodes
    real,intent(in),optional :: Fval

    if (.not.allocated(F)) then
       if (ghost_nodes) then
          allocate(F(C%mbx:C%pbx,C%mby:C%pby,nF))
       else
          allocate(F(C%mx :C%px ,C%my :C%py ,nF))
       end if
       if (present(Fval)) then
          F = Fval
       else
          F = 1d40
       end if
    end if

  end subroutine allocate_array_body_3d


end module mpi_routines2d
