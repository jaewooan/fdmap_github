module fields

  implicit none

  type :: bnd_fields
     ! F = current fields
     ! DF = field rates
     ! F0 = initial fields
     ! DE = rate of energy flow through boundary
     ! (integrate over boundary to get energy dissipation rate)
     ! E = total energy lost through boundary (DE integrated over time)
     logical :: holds_bnd
     integer :: comm
     real,dimension(:,:),allocatable :: U,DU,F,F0,DF,DF2 ! DF is the previous stages DF
     real,dimension(:),allocatable :: E,DE
  end type bnd_fields

  type :: block_fields
     logical :: energy_balance,displacement
     real :: F0(9),DE(4),E(4),Etot
     real,dimension(:),allocatable :: Hx,Hy
     type(bnd_fields) :: bndFL,bndFR,bndFB,bndFT
     real,dimension(:,:,:),allocatable :: Wx,Wy ! PML fields
     real,dimension(:,:,:),allocatable :: DWx,DWy ! PML fields
  end type block_fields

  type :: fields_type
     character(256) :: problem
     logical :: displacement,energy_balance
     integer :: nF,nC,nU,ns
     real :: Etot
     real,dimension(:,:,:),allocatable :: F,U,DU,DF
     real,dimension(:,:),allocatable :: gammap,Dgammap,lambda
     real,dimension(:),allocatable :: Hx,Hy
  end type fields_type

  type :: fields_perturb
    character(256) :: shape
    real :: x0,y0,Lx,Ly,vx,vy,vz,sxx,sxy,sxz,syy,syz,szz
  end type fields_perturb

  real,save :: besselA

  interface rotate_fields_xy2nt
     module procedure rotate_fields_xy2nt_mode3,rotate_fields_xy2nt_mode2
  end interface

  interface rotate_fields_nt2xy
     module procedure rotate_fields_nt2xy_mode3,rotate_fields_nt2xy_mode2
  end interface


contains


  subroutine init_fields(mode,iblock,C,B,G,BF,F,t,input,response, &
      energy_balance,displacement,M,pmlx,pmly)

    use mpi_routines, only : new_communicator
    use mpi_routines2d, only : cartesian
    use grid, only : block_grid,grid_type,set_block_limits
    use io, only : error,seek_to_string
    use utilities, only : deg2rad
    use fd, only : limits,Hnorm
    use material, only : block_material

    implicit none

    integer,intent(in) :: mode,iblock,input
    type(cartesian),intent(in) :: C
    type(block_grid),intent(in) :: B
    type(block_material),intent(in) :: M
    type(grid_type),intent(in) :: G
    type(block_fields),intent(out) :: BF
    type(fields_type),intent(inout) :: F
    real,intent(in) :: t
    character(*),intent(in) :: response
    logical,intent(in) :: energy_balance,displacement,pmlx,pmly

    integer :: stat
    real :: V0,Psi,vx0,vy0,vz0,sxx0,sxy0,sxz0,syy0,syz0,szz0
    character(256) :: problem
    character(256) :: str
    type(fields_perturb) :: P1,P2,P3,P4
    real,dimension(:),allocatable :: F0
    type(limits) :: lim

    namelist /fields_list/ problem,V0,Psi,vx0,vy0,vz0, &
         sxx0,sxy0,sxz0,syy0,syz0,szz0, &
         P1,P2,P3,P4,besselA

    ! defaults

    problem = ''

    V0 = 0d0
    Psi = 0d0

    vx0 = 0d0
    vy0 = 0d0
    vz0 = 0d0

    sxx0 = 0d0
    sxy0 = 0d0
    sxz0 = 0d0
    syy0 = 0d0
    syz0 = 0d0
    szz0 = 0d0

    P1 = fields_perturb('',0d0,0d0,1d0,1d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0)
    P2 = fields_perturb('',0d0,0d0,1d0,1d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0)
    P3 = fields_perturb('',0d0,0d0,1d0,1d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0)
    P4 = fields_perturb('',0d0,0d0,1d0,1d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0)

    besselA = 0d0

    ! read in field parameters

    write(str,'(a,i0,a)') '!---BLOCK',iblock,'---'
    call seek_to_string(input,str)
    read(input,nml=fields_list,iostat=stat)
    if (stat>0) call error('Error in fields_list','init_fields')

    F%problem = problem

    F%energy_balance = energy_balance
    BF%energy_balance = energy_balance
    BF%displacement = displacement
    F%displacement = displacement

    select case(F%problem)
    case default
    case('uniform','rough') ! initial velocity = 0 => uniform state variable (set elsewhere)
       ! syy0 and sxy0 given, set sxx0 using Psi
       if (Psi==45d0) then
          sxx0 = syy0 ! special case, avoid tan(pi/2)=Inf
       else
          sxx0 = (1d0-2d0*sxy0/(syy0*tan(2d0*deg2rad(Psi))))*syy0
       end if
       ! set szz0 as average of sxx0 and syy0
       szz0 = 0.5d0*(sxx0+syy0)
    case('rough-old') ! initial velocity nonzero => heterogeneous state variable for rough fault
       ! syy0 and sxy0 given, set sxx0 using Psi
       if (Psi==45d0) then
          sxx0 = syy0 ! special case, avoid tan(pi/2)=Inf
       else
          sxx0 = (1d0-2d0*sxy0/(syy0*tan(2d0*deg2rad(Psi))))*syy0
       end if
       ! set szz0 as average of sxx0 and syy0
       szz0 = 0.5d0*(sxx0+syy0)
       V0 = 5.4857d-12*exp((sxy0-30.605904d0)/2.016d0)
       if (mod(iblock,2)==0) then
          vx0 =  0.5d0*V0 ! above fault
       else
          vx0 = -0.5d0*V0 ! below fault
       end if
    end select

    F%Etot = 1d40

    ! store initial fields

    BF%F0 = (/ vx0,vy0,vz0,sxx0,sxy0,sxz0,syy0,syz0,szz0 /)

    select case(mode)
    case(2)
       F%nF = 6
       F%nC = 5
       F%nU = 2
       F%ns = 4
       allocate(F0(F%nF))
       F0 = (/ vx0,vy0,sxx0,sxy0,syy0,szz0 /)
    case(3)
       F%nF = 3
       F%nC = 3
       F%nU = 1
       F%ns = 2
       allocate(F0(F%nF))
       F0 = (/ vz0,sxz0,syz0 /)
    end select

    ! allocate memory for fields

    if (.not.allocated(F%F)) then
       allocate(F%F(C%mbx:C%pbx,C%mby:C%pby,F%nF))
       if(F%displacement) allocate(F%U(C%mx:C%px,C%my:C%py,F%nU))
       if(F%displacement) allocate(F%DU(C%mx:C%px,C%my:C%py,F%nU))
       allocate(F%DF(C%mx:C%px,C%my:C%py,F%nF))
       allocate(F%Hx(C%mx:C%px),F%Hy(C%my:C%py))
       F%F = 1d40
       if(F%displacement) F%U = 0d0
       if(F%displacement) F%DU = 1d40
       F%DF = 1d40
       F%Hx = 1d40
       F%Hy = 1d40
       if (index(response,'plastic')/=0) then
          allocate(F%gammap(C%mx:C%px,C%my:C%py),F%Dgammap(C%mx:C%px,C%my:C%py))
          allocate(F%lambda(C%mx:C%px,C%my:C%py))
          F%gammap = 0d0
          F%Dgammap = 0d0
          F%lambda = 0d0
       end if
    end if
    if (pmlx) then
      allocate(BF%Wx(B%mx:B%px,B%my:B%py,F%nF))
      BF%Wx = 0d0
      allocate(BF%DWx(B%mx:B%px,B%my:B%py,F%nF))
      BF%DWx = 1d40
    end if
    if (pmly) then
      allocate(BF%Wy(B%mx:B%px,B%my:B%py,F%nF))
      BF%Wy = 0d0
      allocate(BF%DWy(B%mx:B%px,B%my:B%py,F%nF))
      BF%DWy = 1d40
    end if

    ! initialize fields on sides of this block

    call init_fields_side(B%bndL,BF%bndFL,B%skip,B%my,B%py,F%nF,F%nU,F0,mode,problem, &
         t,iblock,P1,P2,P3,P4,F%energy_balance,F%displacement,M)
    call init_fields_side(B%bndR,BF%bndFR,B%skip,B%my,B%py,F%nF,F%nU,F0,mode,problem, &
         t,iblock,P1,P2,P3,P4,F%energy_balance,F%displacement,M)
    call init_fields_side(B%bndB,BF%bndFB,B%skip,B%mx,B%px,F%nF,F%nU,F0,mode,problem, &
         t,iblock,P1,P2,P3,P4,F%energy_balance,F%displacement,M)
    call init_fields_side(B%bndT,BF%bndFT,B%skip,B%mx,B%px,F%nF,F%nU,F0,mode,problem, &
         t,iblock,P1,P2,P3,P4,F%energy_balance,F%displacement,M)

    ! initialize fields in interior of this block

    call init_fields_interior(B,G,F,F0,mode,problem,t,iblock,P1,P2,P3,P4,M)

    deallocate(F0)

    ! initialize communicators for each side

    call new_communicator(B%sideL,BF%bndFL%comm)
    call new_communicator(B%sideR,BF%bndFR%comm)
    call new_communicator(B%sideB,BF%bndFB%comm)
    call new_communicator(B%sideT,BF%bndFT%comm)

    ! return if process not responsible for block

    if (B%skip) return

    if (F%energy_balance) then

       ! energy and energy dissipation rate

       BF%Etot = 1d40 ! total mechanical energy
       BF%DE = 1d40 ! energy dissipation rate
       BF%E = 0d0 ! cumulative dissipated energy

       ! initialize diagonal norm for this block (assuming unit grid spacing)

       allocate(BF%Hx(B%mx:B%px),BF%Hy(B%my:B%py))

       call set_block_limits(B,lim,'x')
       call Hnorm(lim,BF%Hx)

       call set_block_limits(B,lim,'y')
       call Hnorm(lim,BF%Hy)

       ! initialize norm globally

       F%Hx(B%mx:B%px) = BF%Hx
       F%Hy(B%my:B%py) = BF%Hy

    end if

  end subroutine init_fields


  subroutine destroy_fields(F)

    implicit none

    type(fields_type),intent(inout) :: F

    if (allocated(F%F      )) deallocate(F%F      )
    if (allocated(F%DF     )) deallocate(F%DF     )
    if (allocated(F%U      )) deallocate(F%U      )
    if (allocated(F%DU     )) deallocate(F%DU     )
    if (allocated(F%gammap )) deallocate(F%gammap )
    if (allocated(F%Dgammap)) deallocate(F%Dgammap)
    if (allocated(F%lambda )) deallocate(F%lambda )
    if (allocated(F%Hx     )) deallocate(F%Hx     )
    if (allocated(F%Hy     )) deallocate(F%Hy     )

  end subroutine destroy_fields


  subroutine destroy_block_fields(BF)

    implicit none

    type(block_fields),intent(inout) :: BF

    call destroy_bnd_fields(BF%bndFL)
    call destroy_bnd_fields(BF%bndFR)
    call destroy_bnd_fields(BF%bndFB)
    call destroy_bnd_fields(BF%bndFT)

  end subroutine destroy_block_fields


  subroutine destroy_bnd_fields(bndF)

    implicit none

    type(bnd_fields),intent(inout) :: bndF

    if (allocated(bndF%F  )) deallocate(bndF%F  )
    if (allocated(bndF%U  )) deallocate(bndF%U  )
    if (allocated(bndF%DU )) deallocate(bndF%DU )
    if (allocated(bndF%DF )) deallocate(bndF%DF )
    if (allocated(bndF%DF2)) deallocate(bndF%DF2)
    if (allocated(bndF%F0 )) deallocate(bndF%F0 )
    if (allocated(bndF%E  )) deallocate(bndF%E  )
    if (allocated(bndF%DE )) deallocate(bndF%DE )

  end subroutine destroy_bnd_fields


  subroutine checkpoint_fields(operation,name,checkpoint_number,C,F)

    use io, only : file_distributed,open_file_distributed, &
         read_file_distributed,write_file_distributed,close_file_distributed
    use mpi_routines2d, only : cartesian
    use mpi_routines, only : pw,is_master
    use mpi

    implicit none

    character(*),intent(in) :: operation,name
    integer,intent(in) :: checkpoint_number
    type(cartesian),intent(in) :: C
    type(fields_type),intent(inout) :: F

    type(file_distributed) :: fh
    character(256) :: filename
    integer :: l,ierr

    write(filename,'(a,i0)') trim(adjustl(name)) // 'F.ckpt',checkpoint_number

    if (operation=='delete') then
       if (is_master) call MPI_file_delete(filename,MPI_INFO_NULL,ierr)
       return
    end if

    call open_file_distributed(fh,filename,operation,C%c2d%comm,C%c2d%array_w,pw)

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    select case(operation)
    case('read')
       do l = 1,F%nF
          call read_file_distributed(fh,F%F(C%mx:C%px,C%my:C%py,l))
       end do
       do l = 1,F%nF
          call read_file_distributed(fh,F%DF(:,:,l))
       end do
       do l = 1,F%nU
         if(F%displacement) call read_file_distributed(fh,F%U(C%mx:C%px,C%my:C%py,l))
         if(F%displacement) call read_file_distributed(fh,F%DU(C%mx:C%px,C%my:C%py,l))
       end do
       if (allocated(F%gammap )) call read_file_distributed(fh,F%gammap )
       if (allocated(F%Dgammap)) call read_file_distributed(fh,F%Dgammap)
       if (allocated(F%lambda )) call read_file_distributed(fh,F%lambda )
    case('write')
       do l = 1,F%nF
          call write_file_distributed(fh,F%F(C%mx:C%px,C%my:C%py,l))
       end do
       do l = 1,F%nF
          call write_file_distributed(fh,F%DF(:,:,l))
       end do
       do l = 1,F%nU
         if(F%displacement) call write_file_distributed(fh,F%U(:,:,l))
         if(F%displacement) call write_file_distributed(fh,F%DU(:,:,l))
       end do
       if (allocated(F%gammap )) call write_file_distributed(fh,F%gammap )
       if (allocated(F%Dgammap)) call write_file_distributed(fh,F%Dgammap)
       if (allocated(F%lambda )) call write_file_distributed(fh,F%lambda )
    end select

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    call close_file_distributed(fh)

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

  end subroutine checkpoint_fields


  subroutine checkpoint_block_fields(operation,name,checkpoint_number,iblock,B,BF,F)

    use grid, only : block_grid
    use io, only : file_distributed,open_file_distributed, &
         read_file_distributed,write_file_distributed,close_file_distributed
    use mpi_routines, only : MPI_REAL_PW,pw,subarray,new_communicator,is_master
    use mpi

    implicit none

    character(*),intent(in) :: operation,name
    integer,intent(in) :: checkpoint_number,iblock
    type(block_grid),intent(in) :: B
    type(block_fields),intent(inout) :: BF
    type(fields_type),intent(in) :: F

    type(file_distributed) :: fh
    character(1) :: side_str
    character(256) :: filename,filenameE,filenameWx,filenameWy
    integer :: comm,array,side,l,ierr,rank
    logical :: io_process

    do side = 1,4

       select case(side)
       case(1) ! L
          side_str = 'L'
          io_process = B%sideL
       case(2) ! R
          side_str = 'R'
          io_process = B%sideR
       case(3) ! B
          side_str = 'B'
          io_process = B%sideB
       case(4) ! T
          side_str = 'T'
          io_process = B%sideT
       end select

       write(filename ,'(a,i0,a,i0)') trim(adjustl(name)) // &
            'block',iblock,side_str // '.ckpt',checkpoint_number
       write(filenameE,'(a,i0,a,i0)') trim(adjustl(name)) // &
            'Eblock',iblock,'.ckpt',checkpoint_number
       write(filenameWx,'(a,i0,a,i0)') trim(adjustl(name)) // &
            'pmlWx',iblock,'.ckpt',checkpoint_number
       write(filenameWy,'(a,i0,a,i0)') trim(adjustl(name)) // &
            'pmlWy',iblock,'.ckpt',checkpoint_number

       if (operation=='delete') then
          if (is_master) call MPI_file_delete(filename ,MPI_INFO_NULL,ierr)
          if (is_master) call MPI_file_delete(filenameE,MPI_INFO_NULL,ierr)
          if (side==4) return
          cycle
       end if

       call new_communicator(io_process,comm)

       if (io_process) then
          select case(side)
          case(1,2)
             call subarray(B%pgy-B%mgy+1,B%my-B%mgy+1,B%py-B%mgy+1,MPI_REAL_PW,array)
          case(3,4)
             call subarray(B%pgx-B%mgx+1,B%mx-B%mgx+1,B%px-B%mgx+1,MPI_REAL_PW,array)
          end select
          call open_file_distributed(fh,filename,operation,comm,array,pw)
       end if

       call MPI_Barrier(MPI_COMM_WORLD,ierr)

       if (io_process) then
          select case(operation)
          case('read')
             select case(side)
             case(1) ! L
                do l = 1,F%nF
                   call read_file_distributed(fh,BF%bndFL%F (:,l))
                   call read_file_distributed(fh,BF%bndFL%F0(:,l))
                end do
                do l = 1,F%nU
                  if(F%displacement) call read_file_distributed(fh,BF%bndFL%U(:,l))
                  if(F%displacement) call read_file_distributed(fh,BF%bndFL%DU(:,l))
                end do
                if (BF%energy_balance) call read_file_distributed(fh,BF%bndFL%E)
                if (BF%energy_balance) call read_file_distributed(fh,BF%bndFL%DE)
             case(2) ! R
                do l = 1,F%nF
                   call read_file_distributed(fh,BF%bndFR%F (:,l))
                   call read_file_distributed(fh,BF%bndFR%F0(:,l))
                end do
                do l = 1,F%nU
                  if(F%displacement) call read_file_distributed(fh,BF%bndFR%U(:,l))
                  if(F%displacement) call read_file_distributed(fh,BF%bndFR%DU(:,l))
                end do
                if (BF%energy_balance) call read_file_distributed(fh,BF%bndFR%E)
                if (BF%energy_balance) call read_file_distributed(fh,BF%bndFR%DE)
             case(3) ! B
                do l = 1,F%nF
                   call read_file_distributed(fh,BF%bndFB%F (:,l))
                   call read_file_distributed(fh,BF%bndFB%F0(:,l))
                end do
                do l = 1,F%nU
                  if(F%displacement) call read_file_distributed(fh,BF%bndFB%U(:,l))
                  if(F%displacement) call read_file_distributed(fh,BF%bndFB%DU(:,l))
                end do
                if (BF%energy_balance) call read_file_distributed(fh,BF%bndFB%E)
                if (BF%energy_balance) call read_file_distributed(fh,BF%bndFB%DE)
             case(4) ! T
                do l = 1,F%nF
                   call read_file_distributed(fh,BF%bndFT%F (:,l))
                   call read_file_distributed(fh,BF%bndFT%F0(:,l))
                end do
                do l = 1,F%nU
                  if(F%displacement) call read_file_distributed(fh,BF%bndFT%U(:,l))
                  if(F%displacement) call read_file_distributed(fh,BF%bndFT%DU(:,l))
                end do
                if (BF%energy_balance) call read_file_distributed(fh,BF%bndFT%E)
                if (BF%energy_balance) call read_file_distributed(fh,BF%bndFT%DE)
             end select
          case('write')
             select case(side)
             case(1) ! L
                do l = 1,F%nF
                   call write_file_distributed(fh,BF%bndFL%F (:,l))
                   call write_file_distributed(fh,BF%bndFL%F0(:,l))
                end do
                do l = 1,F%nU
                  if(F%displacement) call write_file_distributed(fh,BF%bndFL%U(:,l))
                  if(F%displacement) call write_file_distributed(fh,BF%bndFL%DU(:,l))
                end do
                if (BF%energy_balance) call write_file_distributed(fh,BF%bndFL%E)
                if (BF%energy_balance) call write_file_distributed(fh,BF%bndFL%DE)
             case(2) ! R
                do l = 1,F%nF
                   call write_file_distributed(fh,BF%bndFR%F (:,l))
                   call write_file_distributed(fh,BF%bndFR%F0(:,l))
                end do
                do l = 1,F%nU
                  if(F%displacement) call write_file_distributed(fh,BF%bndFR%U(:,l))
                  if(F%displacement) call write_file_distributed(fh,BF%bndFR%DU(:,l))
                end do
                if (BF%energy_balance) call write_file_distributed(fh,BF%bndFR%E)
                if (BF%energy_balance) call write_file_distributed(fh,BF%bndFR%DE)
             case(3) ! B
                do l = 1,F%nF
                   call write_file_distributed(fh,BF%bndFB%F (:,l))
                   call write_file_distributed(fh,BF%bndFB%F0(:,l))
                end do
                do l = 1,F%nU
                  if(F%displacement) call write_file_distributed(fh,BF%bndFB%U(:,l))
                  if(F%displacement) call write_file_distributed(fh,BF%bndFB%DU(:,l))
                end do
                if (BF%energy_balance) call write_file_distributed(fh,BF%bndFB%E)
                if (BF%energy_balance) call write_file_distributed(fh,BF%bndFB%DE)
             case(4) ! T
                do l = 1,F%nF
                   call write_file_distributed(fh,BF%bndFT%F (:,l))
                   call write_file_distributed(fh,BF%bndFT%F0(:,l))
                end do
                do l = 1,F%nU
                  if(F%displacement) call write_file_distributed(fh,BF%bndFT%U(:,l))
                  if(F%displacement) call write_file_distributed(fh,BF%bndFT%DU(:,l))
                end do
                if (BF%energy_balance) call write_file_distributed(fh,BF%bndFT%E)
                if (BF%energy_balance) call write_file_distributed(fh,BF%bndFT%DE)
             end select
          end select
       end if

       call MPI_Barrier(MPI_COMM_WORLD,ierr)

       if (io_process) then
          call close_file_distributed(fh)
          call MPI_Type_free(array,ierr)
          call MPI_Comm_free(comm,ierr)
       end if

       call MPI_Barrier(MPI_COMM_WORLD,ierr)

    end do

    if (.not.B%skip.and.BF%energy_balance) then
       call subarray(4,1,4,MPI_REAL_PW,array)
       select case(operation)
       case('read') ! all block processes read same file
          call open_file_distributed(fh,filenameE,operation,MPI_COMM_SELF,array,pw)
          call read_file_distributed(fh,BF%E)
          call read_file_distributed(fh,BF%DE)
          call close_file_distributed(fh)
       case('write') ! only block master process writes file
          call MPI_Comm_rank(B%comm,rank,ierr)
          if (rank==0) then
             call open_file_distributed(fh,filenameE,operation,MPI_COMM_SELF,array,pw)
             call write_file_distributed(fh,BF%E)
             call write_file_distributed(fh,BF%DE)
             call close_file_distributed(fh)
          end if
       end select
       call MPI_Type_free(array,ierr)
    end if

    ! Handle the pml variables
    io_process = .not.B%skip
    call new_communicator(io_process,comm)

    if (.not.B%skip.and.allocated(BF%Wx)) then
       call subarray(B%pgx,B%pgy,F%nf,B%mx,B%px,B%my,B%py,1,F%nf,MPI_REAL_PW,array)
       call open_file_distributed(fh,filenameWx,operation,comm,array,pw)
       select case(operation)
       case('read') ! all block processes read same file
          call read_file_distributed(fh,BF%Wx)
       case('write') ! only block master process writes file
         call write_file_distributed(fh,BF%Wx)
       end select
       call close_file_distributed(fh)
       call MPI_Type_free(array,ierr)
    end if

    if (io_process.and.allocated(BF%Wy)) then
       call subarray(B%pgx,B%pgy,F%nf,B%mx,B%px,B%my,B%py,1,F%nf,MPI_REAL_PW,array)
       call open_file_distributed(fh,filenameWy,operation,comm,array,pw)
       select case(operation)
       case('read') ! all block processes read same file
          call read_file_distributed(fh,BF%Wy)
       case('write') ! only block master process writes file
         call write_file_distributed(fh,BF%Wy)
       end select
       call close_file_distributed(fh)
       call MPI_Type_free(array,ierr)
    end if

    if(io_process) call MPI_Comm_free(comm,ierr)

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

  end subroutine checkpoint_block_fields


  subroutine cycle_stress_fields(F,G,BF,BG,BM,mode,z,H,vp,dt)

    use grid, only : block_grid, grid_type
    use material, only : block_material
    implicit none

    type(fields_type),intent(inout) :: F
    type(grid_type),intent(inout) :: G
    type(block_fields),intent(inout) :: BF
    type(block_grid),intent(in) :: BG
    type(block_material),intent(inout) :: BM

    integer,intent(in) :: mode
    real,intent(in) :: z,H,vp,dt

    integer :: l,lv,i,j
    real :: stress


    if(BG%skip) return

    ! set particle velocities to zero

    select case(mode)
    case(2)
       lv = 2
    case(3)
       lv = 1
    end select

    do l = 1,lv

       do j = BG%my,BG%py
          do i = BG%mx,BG%px
             F%F(i,j,l) = 0d0
          end do
       end do

       if(BG%sideL) then
       BF%bndFL%F (:,l) = 0d0
       BF%bndFL%F0(:,l) = 0d0
       endif

       if(BG%sideR) then
       BF%bndFR%F (:,l) = 0d0
       BF%bndFR%F0(:,l) = 0d0
       endif

       if(BG%sideB) then
       BF%bndFB%F (:,l) = 0d0
       BF%bndFB%F0(:,l) = 0d0
       end if

       if(BG%sideT) then
       BF%bndFT%F (:,l) = 0d0
       BF%bndFT%F0(:,l) = 0d0
       endif

    end do

    select case(mode)
    case(2)
       l = 4 ! sxy
    case(3)
       l = 3
    end select


    do j = BG%my,BG%py
       do i = BG%mx,BG%px
       call interseismic_stress(G%x(i,j),G%y(i,j),z,H,BM%G,vp,dt,stress)
          F%F(i,j,l) = F%F(i,j,l)+stress
       end do
    end do



    if(BG%sideL) then
       do j = BG%my,BG%py

          call interseismic_stress(BG%bndL%x(j),BG%bndL%y(j),z,H,BM%G,vp,dt,stress)

          BF%bndFL%F (j,l) = BF%bndFL%F (j,l)+stress
          BF%bndFL%F0(j,l) = BF%bndFL%F0(j,l)+stress
       enddo
    endif

    if(BG%sideR) then
       do j = BG%my,BG%py
          call interseismic_stress(BG%bndR%x(j),BG%bndR%y(j),z,H,BM%G,vp,dt,stress)

          BF%bndFR%F (j,l) = BF%bndFR%F (j,l)+stress
          BF%bndFR%F0(j,l) = BF%bndFR%F0(j,l)+stress

       end do
    endif

    if(BG%sideB) then
       do i = BG%mx,BG%px
          call interseismic_stress(BG%bndB%x(i),BG%bndB%y(i),z,H,BM%G,vp,dt,stress)

          BF%bndFB%F (i,l) = BF%bndFB%F (i,l)+stress
          BF%bndFB%F0(i,l) = BF%bndFB%F0(i,l)+stress
       enddo
    endif

    if(BG%sideT) then
       do i = BG%mx,BG%px
          call interseismic_stress(BG%bndT%x(i),BG%bndT%y(i),z,H,BM%G,vp,dt,stress)

          BF%bndFT%F (i,l) = BF%bndFT%F (i,l)+stress
          BF%bndFT%F0(i,l) = BF%bndFT%F0(i,l)+stress

       end do
    endif


  end subroutine cycle_stress_fields


  subroutine interseismic_stress(x,y,z,H,shear_mod,vp,dt,stress)

    implicit none

    real,intent(in) :: x,y,z,H,shear_mod,vp,dt  !Dislocation at z = H,  z = h is where cross-section is taken
    real,intent(out) :: stress

    real,parameter :: pi = 3.141592653589793d0


    stress = (shear_mod*vp*dt/(2d0*pi))*(   (z+H)/((z+H)**2 + y**2) + (H-z)/((H-z)**2+y**2)  )

  end subroutine interseismic_stress


  subroutine init_fields_interior(B,G,F,F0,mode,problem,t,iblock,P1,P2,P3,P4,M)

    use grid, only : grid_type,block_grid
    use material, only : block_material

    implicit none

    type(block_grid),intent(in) :: B
    type(grid_type),intent(in) :: G
    type(block_material), intent(in) :: M
    type(fields_type),intent(inout) :: F
    real,intent(in) :: F0(:),t
    integer,intent(in) :: mode,iblock
    character(*),intent(in) :: problem
    type(fields_perturb),intent(in) :: P1,P2,P3,P4

    integer :: i,j
    real :: A1,A2,A3,A4
    real,dimension(F%nF) :: F1,F2,F3,F4

    if (B%skip) return ! no points in this block

    select case(mode)
    case(2)
       F1 = (/ P1%vx,P1%vy,P1%sxx,P1%sxy,P1%syy,P1%szz /)
       F2 = (/ P2%vx,P2%vy,P2%sxx,P2%sxy,P2%syy,P2%szz /)
       F3 = (/ P3%vx,P3%vy,P3%sxx,P3%sxy,P3%syy,P3%szz /)
       F4 = (/ P4%vx,P4%vy,P4%sxx,P4%sxy,P4%syy,P4%szz /)
    case(3)
       F1 = (/ P1%vz,P1%sxz,P1%syz /)
       F2 = (/ P2%vz,P2%sxz,P2%syz /)
       F3 = (/ P3%vz,P3%sxz,P3%syz /)
       F4 = (/ P4%vz,P4%sxz,P4%syz /)
    end select

    do j = B%my,B%py
       do i = B%mx,B%px
          A1 = perturb_fields(G%x(i,j),G%y(i,j),P1)
          A2 = perturb_fields(G%x(i,j),G%y(i,j),P2)
          A3 = perturb_fields(G%x(i,j),G%y(i,j),P3)
          A4 = perturb_fields(G%x(i,j),G%y(i,j),P4)
          F%F(i,j,:) = F0+A1*F1+A2*F2+A3*F3+A4*F4
          call initial_fields(F%F(i,j,:),G%x(i,j),G%y(i,j),t,problem,iblock,M)
       end do
    end do

  end subroutine init_fields_interior


  subroutine init_fields_side(bnd,bndF,Bskip,m,p,nF,nU,F0,mode,problem,t,iblock,&
      P1,P2,P3,P4,energy_balance,displacement,Mat)

    use geometry, only : curve
    use material, only : block_material

    implicit none

    type(curve),intent(in) :: bnd
    type(bnd_fields),intent(inout) :: bndF
    logical,intent(in) :: Bskip,energy_balance,displacement
    integer,intent(in) :: m,p,nF,nU
    real,intent(in) :: F0(:),t
    integer,intent(in) :: mode,iblock
    character(*),intent(in) :: problem
    type(block_material),intent(in) :: Mat
    type(fields_perturb),intent(in) :: P1,P2,P3,P4

    integer :: i
    real :: A1,A2,A3,A4
    real,dimension(nF) :: F1,F2,F3,F4

    select case(mode)
    case(2)
       F1 = (/ P1%vx,P1%vy,P1%sxx,P1%sxy,P1%syy,P1%szz /)
       F2 = (/ P2%vx,P2%vy,P2%sxx,P2%sxy,P2%syy,P2%szz /)
       F3 = (/ P3%vx,P3%vy,P3%sxx,P3%sxy,P3%syy,P3%szz /)
       F4 = (/ P4%vx,P4%vy,P4%sxx,P4%sxy,P4%syy,P4%szz /)
    case(3)
       F1 = (/ P1%vz,P1%sxz,P1%syz /)
       F2 = (/ P2%vz,P2%sxz,P2%syz /)
       F3 = (/ P3%vz,P3%sxz,P3%syz /)
       F4 = (/ P4%vz,P4%sxz,P4%syz /)
    end select

    allocate(bndF%F(m:p,nF),bndF%F0(m:p,nF))
    bndF%F  = 1d40
    bndF%F0 = 1d40

    if(displacement) then
      allocate(bndF%U(m:p,nU))
      allocate(bndF%DU(m:p,nU))
      bndF%U  = 0d0
      bndF%DU  = 1d40
    end if

    if (energy_balance) then
       allocate(bndF%E(m:p),bndF%DE(m:p))
       bndF%E  = 0d0
       bndF%DE = 1d40
    end if

    ! note that boundary fields are allocated (but not initialized with proper values)
    ! even if process is not responsible for this block -- this is because fields on
    ! opposite side of interface (for which another process is responsible) may be needed
    ! to enforce interface conditions (values are set during exchange with other process)

    if (Bskip) return

    do i = m,p
       A1 = perturb_fields(bnd%x(i),bnd%y(i),P1)
       A2 = perturb_fields(bnd%x(i),bnd%y(i),P2)
       A3 = perturb_fields(bnd%x(i),bnd%y(i),P3)
       A4 = perturb_fields(bnd%x(i),bnd%y(i),P4)
       bndF%F0(i,:) = F0+A1*F1+A2*F2+A3*F3+A4*F4
       call initial_fields(bndF%F0(i,:),bnd%x(i),bnd%y(i),t,problem,iblock,Mat,.true.)
    end do

  end subroutine init_fields_side

  subroutine init_fields_side_rate(bndF)

    implicit none

    type(bnd_fields),intent(inout) :: bndF
    integer :: m,p,nF

    if(.not.allocated(bndF%F)) return

    nF = size(bndF%F,2)
    m = lbound(bndF%F,1)
    p = ubound(bndF%F,1)
    allocate(bndF%DF(m:p,nF), bndF%DF2(m:p,nF))
    bndF%DF = 1d40
    bndF%DF2 = 1d40

  end subroutine init_fields_side_rate


  function perturb_fields(x,y,P) result(A)

    use utilities, only : boxcar,gaussian,smooth,triangle
    use io, only : error

    implicit none

    real,intent(in) :: x,y
    type(fields_perturb),intent(in) :: P
    real :: A

    real :: r

    r = sqrt(((x-P%x0)/P%Lx)**2+((y-P%y0)/P%Ly)**2)

    select case(P%shape)
    case default
       call error('Invalid perturbation shape','perturb_fields')
    case('')
       A = 0d0
    case('linear')
       A = (x-P%x0)/P%Lx+(y-P%y0)/P%Ly
    case('boxcar')
       A = boxcar  (r,1d0,1d0,0d0)
    case('triangle')
       A = triangle(r,1d0,1d0,0d0)
    case('smooth')
       A = smooth  (r,1d0,1d0,0d0)
    case('gaussian')
       A = gaussian(r,1d0,1d0,0d0)
    end select

  end function perturb_fields


  subroutine initial_fields(F,x,y,t,problem,iblock,M,bnd_call)

    use utilities, only : step
    use material, only : block_material

    implicit none

    real,dimension(:),intent(inout) :: F
    real,intent(in) :: x,y,t
    character(*),intent(in) :: problem
    integer,intent(in) :: iblock
    type(block_material), intent(in) :: M
    logical,intent(in),optional :: bnd_call
    logical :: b_call

    real :: dip

    if(present(bnd_call)) then
        b_call = bnd_call
    else
        b_call = .false.
    endif

    select case(problem)
    case('break')
        if(b_call) then
            F(1) = F(1)+0d0
            F(2) = F(2)+0d0
            F(3) = F(3)+0d0
        else
            F(1) = F(1)+1d0
            F(2) = F(2)+1d0
            F(3) = F(3)+1d0
        endif
    case('bessel')
       ! verification using method of manufactured solutions
       ! (with bessel function solution in interior)
       ! mode III only
       F(1) = F(1)+bessel(x,y,t,iblock,'vz')
       F(2) = F(2)+bessel(x,y,t,iblock,'sxz')
       F(3) = F(3)+bessel(x,y,t,iblock,'syz')
    case('inplane-bessel')
       F(1) = F(1)+inplane_bessel(x,y,t,M,'vx')
       F(2) = F(2)+inplane_bessel(x,y,t,M,'vy')
       F(3) = F(3)+inplane_bessel(x,y,t,M,'sxx')
       F(4) = F(4)+inplane_bessel(x,y,t,M,'sxy')
       F(5) = F(5)+inplane_bessel(x,y,t,M,'syy')
       F(6) = F(6)+inplane_bessel(x,y,t,M,'szz')
    case('inplane-fault-mms','inplane-fault-mms-nostate')
       F(1) = F(1)+inplane_fault_mms(x,y,t,iblock,'vx')
       F(2) = F(2)+inplane_fault_mms(x,y,t,iblock,'vy')
       F(3) = F(3)+inplane_fault_mms(x,y,t,iblock,'sxx')
       F(4) = F(4)+inplane_fault_mms(x,y,t,iblock,'sxy')
       F(5) = F(5)+inplane_fault_mms(x,y,t,iblock,'syy')
       F(6) = F(6)+inplane_fault_mms(x,y,t,iblock,'szz')
    case('kenneth')
       F(1) = exp(-log(2d0)*(x**2 + (y-0.5d0)**2)/0.01d0)
       F(2) = exp(-log(2d0)*(x**2 + (y-0.5d0)**2)/0.01d0)
       F(3) = 0d0
       F(4) = 0d0
       F(5) = 0d0
       F(6) = 0d0
    case('gaussian')
       ! verification using method of manufactured solutions
       ! (with bessel function solution in interior)
       ! mode III only
       F(1) = F(1)+100d0*exp(-sqrt((x-0d0)**2 + (y+150d0)**2)/2.5d0)
       F(2) = F(2)+100d0*exp(-sqrt((x-0d0)**2 + (y+150d0)**2)/2.5d0)
       ! F(1) = F(1)+100d0*exp(-sqrt((x-0d0)**2 + (y-0d0)**2)/2.5d0)
       ! F(2) = F(2)+100d0*exp(-sqrt((x-0d0)**2 + (y-0d0)**2)/2.5d0)
    case('mms-sin','mms-sin-nostate')
       ! verification using method of manufactured solutions
       ! for curvilinear
       ! mode II only
       F(1) = F(1)+mms_sin(x,y,t,iblock,'vx')
       F(2) = F(2)+mms_sin(x,y,t,iblock,'vy')
       F(3) = F(3)+mms_sin(x,y,t,iblock,'sxx')
       F(4) = F(4)+mms_sin(x,y,t,iblock,'sxy')
       F(5) = F(5)+mms_sin(x,y,t,iblock,'syy')
       F(6) = F(6)+mms_sin(x,y,t,iblock,'szz')
    case('TPV12','TPV13')
       dip = 2d0*x
       F(3) = F(3)+step(dip,13.8d0, 5.8243d0*y,16.6600d0*y) ! sxx
       F(5) = F(5)+step(dip,13.8d0,16.6600d0*y,16.6600d0*y) ! syy
       F(6) = F(6)+step(dip,13.8d0,11.2422d0*y,16.6600d0*y) ! szz
    end select

  end subroutine initial_fields


  subroutine initial_stress(s,x,y,s0,problem,mode)

    ! used only by plasticity routines

    use utilities, only : step

    implicit none

    real,intent(out) :: s(6)
    real,intent(in) :: x,y,s0(6)
    character(*),intent(in) :: problem
    integer,intent(in) :: mode

    real :: dip

    ! sxx  sxy  sxz  syy  syz  szz
    ! s(1) s(2) s(3) s(4) s(5) s(6)

    select case(problem)
    case default ! only stress components unchanged by slip
       select case(mode)
       case(2)
          s(1) = 0d0 ! sxx
          s(2) = 0d0 ! sxy
          s(3) = s0(3)
          s(4) = 0d0 ! syy
          s(5) = s0(5)
          s(6) = 0d0 ! szz
       case(3)
          s(1) = s0(1)
          s(2) = s0(2)
          s(3) = 0d0 ! sxz
          s(4) = s0(4)
          s(5) = 0d0 ! syz
          s(6) = s0(6)
       end select
    case('TPV12alt','TPV13alt')
       s = s0
       dip = 2d0*x
       s(1) = s(1)+step(dip,13.8d0, 5.8243d0*y,16.6600d0*y) ! sxx
       s(4) = s(4)+step(dip,13.8d0,16.6600d0*y,16.6600d0*y) ! syy
       s(6) = s(6)+step(dip,13.8d0,11.2422d0*y,16.6600d0*y) ! szz
    case('YM')
       s = s0
       s(1) = s(1)+ 5.8243d0*y ! sxx
       s(4) = s(4)+16.6600d0*y ! syy
       s(6) = s(6)+11.2422d0*y ! szz
    case('BU')
       s = s0
       s(1) = s(1)+ 8.6d0*y ! sxx
       s(4) = s(4)+15.2d0*y ! syy
       s(6) = s(6)+10.4d0*y ! szz
    case('SB')
       s = s0
       s(1) = s(1)-75.46d0+34.3d0*y ! sxx
       s(4) = s(4)-24.2d0 +11.0d0*y ! syy
       s(6) = s(6)-28.6d0 +13.0d0*y ! szz
    end select

  end subroutine initial_stress


  subroutine scale_rates_interior(C,F,A)

    use mpi_routines2d, only : cartesian

    implicit none

    type(cartesian),intent(in) :: C
    type(fields_type),intent(inout) :: F
    real,intent(in) :: A

    call sri3d(C%mx,C%my,1,F%DF,A,C%mx,C%px,C%my,C%py,1,F%nF)
    if (allocated(F%Dgammap)) &
      call sri2d(C%mx,C%my,F%Dgammap,A,C%mx,C%px,C%my,C%py)
    if(allocated(F%DU)) &
      call sri3d(C%mx,C%my,1,F%DU,A,C%mx,C%px,C%my,C%py,1,F%nU)

  end subroutine scale_rates_interior

  subroutine scale_rates_pml(BF,A,B,nf)
    use grid, only : block_grid

    implicit none

    type(block_fields),intent(inout) :: BF
    type(block_grid),intent(in) :: B
    real,intent(in) :: A
    integer,intent(in) :: nf

    if(allocated(BF%DWx))&
      call sri3d(B%mx,B%my,1,BF%DWx,A,B%mx,B%px,B%my,B%py,1,nF)
      ! BF%DWx = A*BF%DWx
    if(allocated(BF%DWy))&
      call sri3d(B%mx,B%my,1,BF%DWy,A,B%mx,B%px,B%my,B%py,1,nF)
      ! BF%DWy = A*BF%DWy

  end subroutine scale_rates_pml


  subroutine sri2d(lbndx,lbndy,F,A,mx,px,my,py)

    integer,intent(in) :: lbndx,lbndy,mx,px,my,py
    real,dimension(lbndx:,lbndy:),intent(inout) :: F
    real,intent(in) :: A

    integer :: i,j

    do j = my,py
       do i = mx,px
          F(i,j) = A*F(i,j)
       end do
    end do

  end subroutine sri2d


  subroutine sri3d(lbndx,lbndy,lbndz,F,A,mx,px,my,py,mz,pz)

    integer,intent(in) :: lbndx,lbndy,lbndz,mx,px,my,py,mz,pz
    real,dimension(lbndx:,lbndy:,lbndz:),intent(inout) :: F
    real,intent(in) :: A

    integer :: i,j,k

    do k = mz,pz
       do j = my,py
          do i = mx,px
             F(i,j,k) = A*F(i,j,k)
          end do
       end do
    end do

  end subroutine sri3d


  subroutine update_pml(B,BF,dt,nf)

    use grid, only : block_grid

    implicit none

    type(block_fields),intent(inout) :: BF
    type(block_grid),intent(in) :: B
    real,intent(in) :: dt
    integer,intent(in) :: nf

    if(allocated(BF%Wx))&
      call ufi3d(B%mx,B%my,1,B%mx,B%my,1,BF%Wx,BF%DWx,dt,B%mx,B%px,B%my,B%py,1,nF)
      ! BF%Wx = BF%Wx + dt*BF%DWx
    if(allocated(BF%Wy))&
      call ufi3d(B%mx,B%my,1,B%mx,B%my,1,BF%Wy,BF%DWy,dt,B%mx,B%px,B%my,B%py,1,nF)
      ! BF%Wy = BF%Wy + dt*BF%DWy

  end subroutine update_pml

  subroutine update_fields_interior(C,F,dt)

    use mpi_routines2d, only : cartesian

    implicit none

    type(cartesian),intent(in) :: C
    type(fields_type),intent(inout) :: F
    real,intent(in) :: dt

    call ufi3d(C%mbx,C%mby,1,C%mx,C%my,1,F%F,F%DF,dt,C%mx,C%px,C%my,C%py,1,F%nF)
    if (allocated(F%gammap)) &
      call ufi2d(C%mx,C%my,C%mx,C%my,F%gammap,F%Dgammap,dt,C%mx,C%px,C%my,C%py)
    if(F%displacement) then
      ! First add the velocity to the rate
      call ufi3d(C%mx,C%my,1,C%mbx,C%mby,1,F%DU,F%F,1d0,C%mx,C%px,C%my,C%py,1,F%nU)

      ! Then update the displacement
      call ufi3d(C%mx,C%my,1,C%mx,C%my,1,F%U,F%DU,dt,C%mx,C%px,C%my,C%py,1,F%nU)
    end if

  end subroutine update_fields_interior

  subroutine update_displacement(B,BF,dt,nU)

    use grid, only : block_grid

    implicit none

    type(block_grid),intent(in) :: B
    type(block_fields),intent(inout) :: BF
    real,intent(in) :: dt
    integer,intent(in) :: nU

    if (.not.BF%displacement) return

    if (B%sideL) then
      BF%bndFL%U = BF%bndFL%U+dt*BF%bndFL%F(:,1:nU)
      BF%bndFL%U = BF%bndFL%U+dt*BF%bndFL%DU
    end if

    if (B%sideR) then
      BF%bndFR%U = BF%bndFR%U+dt*BF%bndFR%F(:,1:nU)
      BF%bndFR%U = BF%bndFR%U+dt*BF%bndFR%DU
    end if

    if (B%sideB) then
      BF%bndFB%U = BF%bndFB%U+dt*BF%bndFB%F(:,1:nU)
      BF%bndFB%U = BF%bndFB%U+dt*BF%bndFB%DU
    end if

    if (B%sideT) then
      BF%bndFT%U = BF%bndFT%U+dt*BF%bndFT%F(:,1:nU)
      BF%bndFT%U = BF%bndFT%U+dt*BF%bndFT%DU
    end if

  end subroutine update_displacement

  subroutine ufi2d(lbndFx,lbndFy,lbndDFx,lbndDFy,F,DF,dt,mx,px,my,py)

    integer,intent(in) :: lbndFx,lbndFy,lbndDFx,lbndDFy,mx,px,my,py
    real,dimension(lbndFx :,lbndFy :),intent(inout) :: F
    real,dimension(lbndDFx:,lbndDFy:),intent(in) :: DF
    real,intent(in) :: dt

    integer :: i,j

    do j = my,py
       do i = mx,px
          F(i,j) = F(i,j)+dt*DF(i,j)
       end do
    end do

  end subroutine ufi2d


  subroutine ufi3d(lbndFx,lbndFy,lbndFz,lbndDFx,lbndDFy,lbndDFz,F,DF,dt,mx,px,my,py,mz,pz)

    integer,intent(in) :: lbndFx,lbndFy,lbndFz,lbndDFx,lbndDFy,lbndDFz,mx,px,my,py,mz,pz
    real,dimension(lbndFx :,lbndFy :,lbndFz :),intent(inout) :: F
    real,dimension(lbndDFx:,lbndDFy:,lbndDFz:),intent(in) :: DF
    real,intent(in) :: dt

    integer :: i,j,k

    do k = mz,pz
       do j = my,py
          do i = mx,px
             F(i,j,k) = F(i,j,k)+dt*DF(i,j,k)
          end do
       end do
    end do

  end subroutine ufi3d


  subroutine exchange_fields(C,F)

    use mpi_routines2d, only : cartesian,exchange_all_neighbors

    implicit none

    type(cartesian),intent(in) :: C
    type(fields_type),intent(inout) :: F

    integer :: l

    do l = 1,F%nF
       call exchange_all_neighbors(C,F%F(:,:,l))
    end do

  end subroutine exchange_fields


  subroutine rotate_fields_xy2nt_mode3(F,normal,vz,stz,snz)

    use geometry, only : rotate_xy2nt

    implicit none

    real,intent(in) :: F(3),normal(2)
    real,intent(out) :: vz,stz,snz

    vz  = F(1)
    call rotate_xy2nt(F(2),F(3),stz,snz,normal)

  end subroutine rotate_fields_xy2nt_mode3


  subroutine rotate_fields_xy2nt_mode2(F,normal,vt,vn,stt,snt,snn,szz)

    use geometry, only : rotate_xy2nt

    implicit none

    real,intent(in) :: F(6),normal(2)
    real,intent(out) :: vt,vn,stt,snt,snn,szz

    call rotate_xy2nt(F(1),F(2),vt,vn,normal)
    call rotate_xy2nt(F(3),F(4),F(5),stt,snt,snn,normal)
    szz = F(6)

  end subroutine rotate_fields_xy2nt_mode2


  subroutine rotate_fields_nt2xy_mode3(F,normal,vz,stz,snz)

    use geometry, only : rotate_nt2xy

    implicit none

    real,intent(out) :: F(3)
    real,intent(in) :: normal(2),vz,stz,snz

    F(1) = vz
    call rotate_nt2xy(F(2),F(3),stz,snz,normal)

  end subroutine rotate_fields_nt2xy_mode3


  subroutine rotate_fields_nt2xy_mode2(F,normal,vt,vn,stt,snt,snn,szz)

    use geometry, only : rotate_nt2xy

    implicit none

    real,intent(out) :: F(6)
    real,intent(in) :: normal(2),vt,vn,stt,snt,snn,szz

    call rotate_nt2xy(F(1),F(2),vt,vn,normal)
    call rotate_nt2xy(F(3),F(4),F(5),stt,snt,snn,normal)
    F(6) = szz

  end subroutine rotate_fields_nt2xy_mode2


  function bessel(x,y,t,side,field) result(F)

    use ifport ! for Intel Fortran compiler

    implicit none

    real,intent(in) :: x,y,t
    integer,intent(in) :: side
    character(*),intent(in) :: field
    real :: F

    real :: r,J,dJdr,drdx,drdy,F0
    real,parameter :: vz0=1d0,syz0=20d0

    r = sqrt(x**2+y**2)

    J = dbesj0(r)
    dJdr = -dbesj1(r)

    if (r<epsilon(r)) then

       select case(field)
       case('vz')
          F0 = J*cos(t)
          F = (vz0+F0)-0.5d0*besselA*F0
          if (side==1) F = -F
       case('V')
          F0 = J*cos(t)
          F = 2d0*(vz0+F0)-besselA*F0
       case('sxz')
          F = 0d0
       case('syz','S')
          F = syz0
       end select

    else

       select case(field)
       case('vz')
          F0 = J*cos(t)
          F = (vz0+F0)-0.5d0*besselA*F0
          if (side==1) F = -F
       case('V')
          F0 = J*cos(t)
          F = 2d0*(vz0+F0)-besselA*F0
       case('sxz')
          drdx = x/r
          F0 = drdx*dJdr*sin(t)
          F = (1d0-0.5d0*besselA)*F0
          if (side==1) F = -F
       case('syz')
          drdy = y/r
          F0 = drdy*dJdr*sin(t)
          F = (1d0-0.5d0*besselA)*F0
          if (side==1) F = -F
          F = F+syz0
       case('S')
          F = syz0
       end select

    end if

  end function bessel


  function inplane_bessel(x,y,t,M,field) result(F)

    use ifport ! for Intel Fortran compiler
    use material, only : block_material
    use io, only : error

    implicit none

    real,intent(in) :: x,y,t
    character(*),intent(in) :: field
    type(block_material),intent(in) :: M
    real :: F

    real :: r,r_x,r_y,r_xy
    real :: rp,rp_x,rp_y
    real :: J0p,J1p,J2p
    real :: ftp,ftp_t
    real :: kp,kcp
    real :: rs,rs_x,rs_y
    real :: J0s,J1s,J2s
    real :: fts,fts_t
    real :: ks,kcs
    real :: ux,uy,ux_x,ux_y,uy_x,uy_y
    real :: pi

    pi= 4d0 * datan(1d0)
    kp = 20d0*pi
    ks = 10d0*pi

    kcp = kp/M%cp
    kcs = ks/M%cs

    F = 1d40

    r = sqrt(x**2+y**2)
    if (r<epsilon(r)) then
      r = 0d0

      r_x = 0d0
      r_y = 0d0
    else
      r_x = x / r
      r_y = y / r
    end if

    ! P-wave
    rp     = kcp*r
    rp_x   = kcp*r_x
    rp_y   = kcp*r_y

    J0p   = dbesjn(0,rp)
    J1p   = dbesjn(1,rp)
    J2p   = dbesjn(2,rp)
    ftp    = sin(kp*t)
    ftp_t  = kp*cos(kp*t)

    ! S-wave
    rs     = kcs*r
    rs_x   = kcs*r_x
    rs_y   = kcs*r_y

    J0s   = dbesjn(0,rs)
    J1s   = dbesjn(1,rs)
    J2s   = dbesjn(2,rs)
    fts    = sin(ks*t)
    fts_t  = ks*cos(ks*t)

    select case(field)
    case('vx')
      F = -ftp_t * rp_x * J1p - fts_t * rs_y * J1s
    case('vy')
      F = -ftp_t * rp_y * J1p + fts_t * rs_x * J1s
    case('sxx','syy')
      ux_x = -ftp*0.5d0*kcp**2
      uy_y = -ftp*0.5d0*kcp**2
      if (r>epsilon(r)) then
        ux_x = ftp*(kcp*(x**2-y**2)*J1p - kcp**2*x**2*r*J0p)/r**3 + fts*kcs**2*x*y*J2s/r**2
        uy_y = ftp*(kcp*(y**2-x**2)*J1p - kcp**2*y**2*r*J0p)/r**3 - fts*kcs**2*x*y*J2s/r**2
      end if
      select case(field)
      case('sxx')
        F = (M%lambda+2d0*M%G)*ux_x + M%lambda*uy_y
      case('syy')
        F = M%lambda*ux_x + (M%lambda+2d0*M%G)*uy_y
      end select
    case('sxy')
      ux_x = -fts*0.5d0*kcs**2
      uy_y =  fts*0.5d0*kcs**2
      if (r>epsilon(r)) then
        ux_y = ftp*kcp**2*x*y*J2p/r**2 + fts*(kcs*(y**2-x**2)*J1s - kcs**2*y**2*r*J0s)/r**3
        uy_x = ftp*kcp**2*x*y*J2p/r**2 - fts*(kcs*(x**2-y**2)*J1s - kcs**2*x**2*r*J0s)/r**3
      end if

      F = M%G*(ux_y+uy_x)
    case('szz')
      F = 0d0
    end select

  end function inplane_bessel

  function inplane_fault_mms(x,y,t,iblock,field) result(F)

    use ifport ! for Intel Fortran compiler
    use material, only : block_material
    use io, only : error

    implicit none

    real,intent(in) :: x,y,t
    integer,intent(in) :: iblock
    character(*),intent(in) :: field
    real :: F

    real :: ftp,ftp_t,ftp_tt
    real :: kp,kcp

    real :: vx  ,vy
    real :: vx_t,vy_t
    real :: vx_x,vy_x
    real :: vx_y,vy_y

    real :: v, v_x, v_y, v_t

    real :: nx  ,ny
    real :: nx_x,ny_x
    real :: nx_y,ny_y

    real :: mx  ,my
    real :: mx_x,my_x
    real :: mx_y,my_y

    real :: sxx  ,syy  ,sxy  ,szz
    real :: sxx_t,syy_t,sxy_t,szz_t
    real :: sxx_x,syy_x,sxy_x
    real :: sxx_y,syy_y,sxy_y

    real :: pi

    real :: G = 32d0, cs = 3d0, cp = 5d0, rho, lambda
    rho = G/cs**2
    lambda = rho*cp**2 - 2d0*G

    F = mms_sin(x,y,t,iblock,field)
    return

    pi = 4d0 * datan(1d0)


    kp = 10d0*pi
    kcp = kp/cp

    ftp    =       sin(kp*t)/kp
    ftp_t  =    kp*cos(kp*t)/kp
    ftp_tt =-kp**2*sin(kp*t)/kp

    szz   = 0d0
    szz_t = 0d0

    sxx    =     ftp  *sin(kcp*x)*sin(kcp*y) - 126d0
    sxx_t  =     ftp_t*sin(kcp*x)*sin(kcp*y)
    sxx_x  = kcp*ftp  *cos(kcp*x)*sin(kcp*y)
    sxx_y  = kcp*ftp  *sin(kcp*x)*cos(kcp*y)

    syy    =     ftp  *sin(kcp*x)*sin(kcp*y) - 126d0
    syy_t  =     ftp_t*sin(kcp*x)*sin(kcp*y)
    syy_x  = kcp*ftp  *cos(kcp*x)*sin(kcp*y)
    syy_y  = kcp*ftp  *sin(kcp*x)*cos(kcp*y)

    sxy    =     ftp  *sin(kcp*x)*sin(kcp*y) + 0.6*126d0
    sxy_t  =     ftp_t*sin(kcp*x)*sin(kcp*y)
    sxy_x  = kcp*ftp  *cos(kcp*x)*sin(kcp*y)
    sxy_y  = kcp*ftp  *sin(kcp*x)*cos(kcp*y)

    ! WARNING: This must match the coordinate transform!!!!!
    ! nx =  pi*cos(pi*x)
    ! ny = -1d0

    ! nx_x = -pi**2*sin(pi*x)
    ! ny_x = 0d0

    ! nx_y = 0d0
    ! ny_y = 0d0
    nx =  pi*cos(pi*x) / sqrt(10d0**2 + (pi*cos(pi*x))**2)
    ny =         -10d0 / sqrt(10d0**2 + (pi*cos(pi*x))**2)

    nx_x = -10d0**2*pi**2*sin(    pi*x) / sqrt(10d0**2 + (pi*cos(pi*x))**2)**3
    ny_x =     -5d0*pi**3*sin(2d0*pi*x) / sqrt(10d0**2 + (pi*cos(pi*x))**2)**3

    nx_y = 0d0
    ny_y = 0d0

    mx =  ny
    my = -nx

    mx_x =  ny_x
    my_x = -nx_x

    mx_y =  ny_y
    my_y = -nx_y

    v   = (3d0-2d0*iblock)*    (ftp_t *sin(kcp*x)*sin(kcp*y) - 1.1d0)
    v_t = (3d0-2d0*iblock)*    (ftp_tt*sin(kcp*x)*sin(kcp*y))
    v_x = (3d0-2d0*iblock)*kcp*(ftp_t *cos(kcp*x)*sin(kcp*y))
    v_y = (3d0-2d0*iblock)*kcp*(ftp_t *sin(kcp*x)*cos(kcp*y))

    vx   = -  v*ny
    vx_t = -v_t*ny
    vx_x = -v_x*ny - v*ny_x
    vx_y = -v_y*ny - v*ny_y

    vy   = v  *nx
    vy_t = v_t*nx
    vy_x = v_x*nx + v*nx_x
    vy_y = v_y*nx + v*nx_y

    select case(field)
    case default
       call error('Invalid field (' // trim(field) // ')')
    case('V')
      F = 2d0*(mx*vx + my*vy)
    case('N')
      F = -(nx*sxx*nx + ny*syy*ny + 2d0*nx*ny*sxy)
    case('S')
      F = nx*sxx*mx + ny*syy*my + (nx*my + ny*mx)*sxy
    case('Vt')
      F = 2d0*(mx*vx_t + my*vy_t)
    case('Nt')
      F = -(nx*sxx_t*nx + ny*syy_t*ny + 2d0*nx*ny*sxy_t)
    case('St')
      F = nx*sxx_t*mx + ny*syy_t*my + (nx*my + ny*mx)*sxy_t
    case('vx')
      F = vx
    case('vy')
      F = vy
    case('s_vx')
      F = vx_t-(sxx_x+sxy_y)/rho
    case('s_vy')
      F = vy_t-(sxy_x+syy_y)/rho
    case('sxx')
      F = sxx
    case('syy')
      F = syy
    case('sxy')
      F = sxy
    case('szz')
      F = szz
    case('s_sxx')
      F = sxx_t-((lambda+2d0*G)*vx_x + lambda*vy_y)
    case('s_syy')
      F = syy_t-(lambda*vx_x + (lambda+2d0*G)*vy_y)
    case('s_sxy')
      F = sxy_t-(G*(vy_x+vx_y))
    case('s_szz')
      F = szz_t-(lambda*(vx_x+vy_y))
    end select

  end function inplane_fault_mms


  function inplane_fault_mms_old(x,y,t,M,field) result(F)

    use ifport ! for Intel Fortran compiler
    use material, only : block_material
    use io, only : error

    implicit none

    real,intent(in) :: x,y,t
    character(*),intent(in) :: field
    type(block_material),intent(in) :: M
    real :: F

    real :: r,r_x,r_y,r_xy,r_xx,r_yy
    real :: rp,rp_x,rp_y
    real :: ftp,ftp_t,ftp_tt
    real :: gp,gp_r,gp_rr
    real :: gp_x,gp_y,gp_xx,gp_yy,gp_xy
    real :: kp,kcp
    real :: ux,uy
    real :: ux_t,uy_t
    real :: ux_x,uy_x
    real :: ux_y,uy_y
    real :: ux_xx,uy_xx
    real :: ux_xy,uy_xy
    real :: ux_yy,uy_yy
    real :: ux_tt,uy_tt
    real :: pi

    pi = 4d0 * datan(1d0)
    r = sqrt(x**2+y**2)

    r_x = x / r
    r_y = y / r

    r_xx = y**2 / r**3
    r_yy = x**2 / r**3
    r_xy = x*y  / r**3

    kp = 10d0*pi
    kcp = kp/M%cp

    ftp    =       sin(kp*t)
    ftp_t  =    kp*cos(kp*t)
    ftp_tt =-kp**2*sin(kp*t)

    gp    =        cos(kcp*r)
    gp_r  =   -kcp*sin(kcp*r)
    gp_rr =-kcp**2*cos(kcp*r)

    gp_x = gp_r*r_x
    gp_y = gp_r*r_y

    gp_xx = gp_rr*r_x**2  + gp_r*r_xx
    gp_yy = gp_rr*r_y**2  + gp_r*r_yy
    gp_xy = gp_rr*r_y*r_x + gp_r*r_xy

    if(r < epsilon(r)) then
      gp_x = 0d0
      gp_y = 0d0

      gp_xx = -kcp**2
      gp_yy = -kcp**2
      gp_xy = 0d0
    end if
    
    gp = x**2 + y**2
    gp_x = 2d0*x
    gp_xx = 2d0
    gp_y = 2d0*y
    gp_yy = 2d0
    gp_xy = 0d0
    
    gp    =        sin(kcp*x)*sin(kcp*y)
    gp_x  =    kcp*cos(kcp*x)*sin(kcp*y)
    gp_xx =-kcp**2*sin(kcp*x)*sin(kcp*y)
    gp_y  =    kcp*sin(kcp*x)*cos(kcp*y)
    gp_yy =-kcp**2*sin(kcp*x)*sin(kcp*y)
    gp_xy = kcp**2*cos(kcp*x)*cos(kcp*y)
    
    ux    = ftp   *gp
    ux_t  = ftp_t *gp
    ux_tt = ftp_tt*gp
    ux_x  = ftp   *gp_x
    ux_y  = ftp   *gp_y
    ux_xx = ftp   *gp_xx
    ux_xy = ftp   *gp_xy
    ux_yy = ftp   *gp_yy
    
    uy    = ftp   *gp
    uy_t  = ftp_t *gp
    uy_tt = ftp_tt*gp
    uy_x  = ftp   *gp_x
    uy_y  = ftp   *gp_y
    uy_xx = ftp   *gp_xx
    uy_xy = ftp   *gp_xy
    uy_yy = ftp   *gp_yy
    
    select case(field)
    case('vx')
      F = ux_t
    case('vy')
      F = uy_t
    case('s_vx')
      F = ux_tt-(&
        (M%lambda+2d0*M%G)*ux_xx + M%lambda*uy_xy&
        +M%G*(ux_yy+uy_xy)&
        )/M%rho
    case('s_vy')
      F = uy_tt-(&
        M%G*(ux_xy+uy_xx)&
        +M%lambda*ux_xy + (M%lambda+2d0*M%G)*uy_yy&
        )/M%rho
    case('sxx')
      F = (M%lambda+2d0*M%G)*ux_x + M%lambda*uy_y
    case('syy')
      F = M%lambda*ux_x + (M%lambda+2d0*M%G)*uy_y
    case('sxy')
      F = M%G*(ux_y+uy_x)
    case('szz','s_sxx','s_syy','s_sxy','s_szz')
      F = 0d0
    end select

  end function inplane_fault_mms_old

  function mms_sin_old(x,y,t,side,field,M) result(F)

    use ifport ! for Intel Fortran compiler
    use material, only : block_material
    use io, only : error

    implicit none

    real,intent(in) :: x,y,t
    integer,intent(in) :: side
    character(*),intent(in) :: field
    type(block_material),intent(in),optional :: M
    real :: F, F2

    real :: syy,sxy,sxx,vx,vy,g,gp,pi,mv
    real :: syy0,sxy0,sxx0,vx0,vy0
    real :: vx_t,vx_x,vx_y
    real :: vy_t,vy_x,vy_y
    real :: syy_t,syy_x,syy_y
    real :: sxy_t,sxy_x,sxy_y
    real :: sxx_t,sxx_x,sxx_y
    real :: szz_t
    real :: nx,ny,mx,my
    real :: ft,ft_t
    real :: nx_x,ny_x,mx_x,my_x,gp_x
    real :: nx_y,ny_y,mx_y,my_y
    real :: w,kx,ky,mf,kf

    pi = 4d0 * datan(1d0)
    kx  = 2d0*pi
    ky  = 2d0*pi
    mv  = 2d0
    w  = (2d0/0.46d0)*pi
    ft = cos(w*t)
    ft_t = -w*sin(w*t)

    ! WARNING: This must match the coordinate transform!!!!!
    mf   = 0.1d0
    kf   = 1d0*pi
    g    = mf*sin(kf*x)
    gp   = mf*kf*cos(kf*x)
    gp_x =-mf*(kf**2)*sin(kf*x)
    ! g = 0d0
    ! gp = 0d0
    ! gp_x = 0d0

    nx =  gp/sqrt(1+gp**2)
    ny = -1 /sqrt(1+gp**2)
    mx =  ny
    my = -nx

    ! need to be careful since slip and sxy must have the same sign at the
    ! fault for the friction law
    vx   = mv*(3d0-2d0*side)*   1d0*cos(kx*x)*cos(ky*y)
    vx_x =-mv*(3d0-2d0*side)*kx*1d0*sin(kx*x)*cos(ky*y)
    vx_y =-mv*(3d0-2d0*side)*ky*1d0*cos(kx*x)*sin(ky*y)

    ! vx   =   cos(kx*x)*cos(ky*y)
    ! vx_x =-kx*sin(kx*x)*cos(ky*y)
    ! vx_y =-ky*cos(kx*x)*sin(ky*y)

    vy   =    1d0*cos(kx*x)*cos(ky*y)
    vy_x =-kx*1d0*sin(kx*x)*cos(ky*y)
    vy_y =-ky*1d0*cos(kx*x)*sin(ky*y)

    syy   =   -50d0*cos(kx*x)*cos(ky*y)
    syy_x = kx*50d0*sin(kx*x)*cos(ky*y)
    syy_y = ky*50d0*cos(kx*x)*sin(ky*y)

    sxy   =    20d0*cos(kx*x)*cos(ky*y)
    sxy_x =-kx*20d0*sin(kx*x)*cos(ky*y)
    sxy_y =-ky*20d0*cos(kx*x)*sin(ky*y)

    sxx   = (2d0*side-3d0)*   50d0*cos(kx*x)*cos(ky*y)
    sxx_x =-(2d0*side-3d0)*kx*50d0*sin(kx*x)*cos(ky*y)
    sxx_y =-(2d0*side-3d0)*ky*50d0*cos(kx*x)*sin(ky*y)

    ! sxx   =   cos(kx*x)*cos(ky*y)
    ! sxx_x =-kx*sin(kx*x)*cos(ky*y)
    ! sxx_y =-ky*cos(kx*x)*sin(ky*y)

    vx0  = (mv+1)*(3d0-2d0*side)
    vy0  = 0d0
    syy0 = -126d0
    sxy0 = -syy0*0.6d0
    sxx0 = 0d0

    ! sxy = 0d0
    ! sxx = 0d0
    ! syy = 0d0
    ! vx = 0d0
    ! vy = 0d0

    select case(field)
    case('szz')
        F = 0d0
        return
    case('vx','vy','V')
        vx = ft*vx + vx0
        vy = ft*vy + vy0
        select case(field)
        case('V')
            F = 2d0*vx
            return
        case('vx')
            F = ny*vx-my*vy
            return
        case('vy')
            F = mx*vy-nx*vx
            return
        end select
    case('sxx','syy','sxy','N','S')
        syy = ft*syy + syy0
        sxx = ft*sxx + sxx0
        sxy = ft*sxy + sxy0
        select case(field)
        case('S')
            F = sxy
            return
        case('N')
            F =-syy
            return
        case('sxx')
            F = syy*my**2+sxx*ny**2-2d0*sxy*my*ny
            return
        case('syy')
            F = syy*mx**2+sxx*nx**2-2d0*sxy*mx*nx
            return
        case('sxy')
            F =-syy*my*mx-sxx*ny*nx+sxy*(ny*mx+nx*my)
            return
        end select
    case('Vt')
         F = 2d0*vx*ft_t
         return
    case('St')
         F = sxy*ft_t
         return
    case('Nt')
         F =  -syy*ft_t
         return
    case('s_vx','s_vy')
        ! need to be careful since slip and sxy must have the same sign at the
        ! fault for the friction law

        ! Need the time derivatives of velocity
        vx_t = ft_t*vx
        vy_t = ft_t*vy

        ! For the stresses we need both parts
        syy = ft*syy + syy0
        sxx = ft*sxx + sxx0
        sxy = ft*sxy + sxy0

        syy_x = syy_x*ft
        sxy_x = sxy_x*ft
        sxx_x = sxx_x*ft

        syy_y = syy_y*ft
        sxy_y = sxy_y*ft
        sxx_y = sxx_y*ft

        ! We will also need the derivatives of the normal vectors
        nx_x =     gp_x/sqrt((1+gp**2)**3)
        ny_x =  gp*gp_x/sqrt((1+gp**2)**3)
        mx_x =  ny_x
        my_x = -nx_x

        nx_y = 0d0
        ny_y = 0d0
        mx_y =  ny_y
        my_y = -nx_y

        select case(field)
        case('s_vx')
            ! vx_t - (1/rho)*(sxx_x+sxy_y)
            vx_t  = ny*vx_t-my*vy_t
            sxx_x = syy_x*my**2+sxx_x*ny**2-2*sxy_x*my*ny &
                  + 2d0*syy*my*my_x+2d0*sxx*ny*ny_x-2*sxy*(my_x*ny+my*ny_x)
            sxy_y = -syy_y*my*mx-sxx_y*ny*nx+sxy_y*(ny*mx+nx*my) &
                  -syy*(my_y*mx+my*mx_y)-sxx*(ny_y*nx+ny*nx_y)+sxy*(ny_y*mx+ny*mx_y+nx_y*my+nx*my_y)
            F = vx_t - (1d0/M%rho)*(sxx_x+sxy_y)
            return
        case('s_vy')
            ! vy_t - (1/rho)*(sxy_x+syy_y)
            vy_t = mx*vy_t-nx*vx_t
            sxy_x = -syy_x*my*mx-sxx_x*ny*nx+sxy_x*(ny*mx+nx*my) &
                  -syy*(my_x*mx+my*mx_x)-sxx*(ny_x*nx+ny*nx_x)+sxy*(ny_x*mx+ny*mx_x+nx_x*my+nx*my_x)
            syy_y = syy_y*mx**2+sxx_y*nx**2-2*sxy_y*mx*nx &
                  + 2d0*syy*mx*mx_y+2d0*sxx*nx*nx_y-2*sxy*(mx_y*nx+mx*nx_y)
            F = vy_t - (1d0/M%rho)*(sxy_x+syy_y)
            return
        end select
    case('s_sxx','s_syy','s_sxy','s_szz')
        ! Need the time derivatives of stresses
        syy_t = ft_t*syy
        sxy_t = ft_t*sxy
        sxx_t = ft_t*sxx

        ! For the velocities we need both parts
        vx = ft*vx + vx0
        vy = ft*vy + vy0

        vx_x = vx_x*ft
        vy_x = vy_x*ft

        vx_y = vx_y*ft
        vy_y = vy_y*ft

        ! We will also need the derivatives of the normal vectors
        nx_x =     gp_x/sqrt((1+gp**2)**3)
        ny_x =  gp*gp_x/sqrt((1+gp**2)**3)
        mx_x =  ny_x
        my_x = -nx_x

        nx_y = 0d0
        ny_y = 0d0
        mx_y =  ny_y
        my_y = -nx_y

        select case(field)
        case('s_sxx')
            ! sxx_t - (lam+2*G) vx_x - lam * vy_y
            sxx_t = syy_t*my**2+sxx_t*ny**2-2*sxy_t*my*ny
            vx_x = ny_x*vx-my_x*vy + ny*vx_x-my*vy_x
            vy_y = mx_y*vy-nx_y*vx + mx*vy_y-nx*vx_y
            F = sxx_t - (M%lambda+2d0*M%G)*vx_x - M%lambda*vy_y
            return
        case('s_syy')
            ! syy_t - lam * vx_x - (lam+2*G) vy_y
            vx_x = ny_x*vx-my_x*vy + ny*vx_x-my*vy_x
            vy_y = mx_y*vy-nx_y*vx + mx*vy_y-nx*vx_y
            syy_t = syy_t*mx**2+sxx_t*nx**2-2*sxy_t*mx*nx
            F = syy_t - M%lambda*vx_x - (M%lambda+2d0*M%G)*vy_y
            return
        case('s_sxy')
            ! sxy_t - G * vy_x + G * vx_y
            vx_y = ny_y*vx-my_y*vy + ny*vx_y-my*vy_y
            vy_x = mx_x*vy-nx_x*vx + mx*vy_x-nx*vx_x
            sxy_t =-syy_t*my*mx-sxx_t*ny*nx+sxy_t*(ny*mx+nx*my)
            F = sxy_t - M%G*(vy_x+vx_y)
            return
        case('s_szz')
            ! sxx_t - (lam+2*G) vx_x - lam * vy_y
            szz_t = 0d0
            vx_x = ny_x*vx-my_x*vy + ny*vx_x-my*vy_x
            vy_y = mx_y*vy-nx_y*vx + mx*vy_y-nx*vx_y
            F = szz_t - M%lambda*vx_x - M%lambda*vy_y
            return
        end select
    end select

       call error('Invalid field (' // trim(field) // ')')

    F = 0

  end function mms_sin_old

  function mms_sin(x,y,t,side,field) result(F)

    use ifport ! for Intel Fortran compiler
    use material, only : block_material
    use io, only : error

    implicit none

    real,intent(in) :: x,y,t
    integer,intent(in) :: side
    character(*),intent(in) :: field
    real :: F, F2

    real :: syy,sxy,sxx,vx,vy,ftp,ftpp,pi,mv
    real :: vx_t,vx_x,vx_y
    real :: vy_t,vy_x,vy_y
    real :: syy_t,syy_x,syy_y
    real :: sxy_t,sxy_x,sxy_y
    real :: sxx_t,sxx_x,sxx_y
    real :: szz_t
    real :: nx,ny,mx,my
    real :: ft,ft_t
    real :: nx_x,ny_x,mx_x,my_x,ftpp_x
    real :: nx_y,ny_y,mx_y,my_y
    real :: w,kx,ky,mf,kf

    real :: G = 32d0, cs = 3d0, cp = 5d0, rho, lambda, Zs
    rho = G/cs**2
    lambda = rho*cp**2 - 2d0*G
    Zs = rho*cs


    pi = 4d0 * datan(1d0)
    kx  = 2d0*pi
    ky  = 2d0*pi
    mv  = 2d0
    w  = (4d0/0.46d0)*pi
    ! w = 10d0*pi
    ft = cos(w*t)
    ft_t = -w*sin(w*t)

    ! WARNING: This must match the coordinate transform!!!!!
    mf   = 0.1d0
    kf   = 1d0*pi
    ftp    = mf*sin(kf*x)
    ftpp   = mf*kf*cos(kf*x)
    ftpp_x =-mf*(kf**2)*sin(kf*x)
    ! ftp = 0d0
    ! ftpp = 0d0
    ! ftpp_x = 0d0

    nx =  ftpp/sqrt(1d0+ftpp**2)
    ny = -1d0 /sqrt(1d0+ftpp**2)
    mx =  ny
    my = -nx

    ! We will also need the derivatives of the normal vectors
    nx_x =     ftpp_x/sqrt((1+ftpp**2)**3)
    ny_x =  ftpp*ftpp_x/sqrt((1+ftpp**2)**3)
    mx_x =  ny_x
    my_x = -nx_x

    nx_y = 0d0
    ny_y = 0d0
    mx_y =  ny_y
    my_y = -nx_y


    ! need to be careful since slip and sxy must have the same sign at the
    ! fault for the friction law
    vx   = mv*(3d0-2d0*side)*   (2d0*12.6/Zs)*cos(kx*x)*cos(ky*y)*ft
    vx_t = mv*(3d0-2d0*side)*   (2d0*12.6/Zs)*cos(kx*x)*cos(ky*y)*ft_t
    vx_x =-mv*(3d0-2d0*side)*kx*(2d0*12.6/Zs)*sin(kx*x)*cos(ky*y)*ft
    vx_y =-mv*(3d0-2d0*side)*ky*(2d0*12.6/Zs)*cos(kx*x)*sin(ky*y)*ft

    vy   =    (2d0*12.6/Zs)*cos(kx*x)*cos(ky*y)*ft
    vy_t =    (2d0*12.6/Zs)*cos(kx*x)*cos(ky*y)*ft_t
    vy_x =-kx*(2d0*12.6/Zs)*sin(kx*x)*cos(ky*y)*ft
    vy_y =-ky*(2d0*12.6/Zs)*cos(kx*x)*sin(ky*y)*ft

    syy   =   -63d0*cos(kx*x)*cos(ky*y)*ft - 126d0
    syy_t =   -63d0*cos(kx*x)*cos(ky*y)*ft_t
    syy_x = kx*63d0*sin(kx*x)*cos(ky*y)*ft
    syy_y = ky*63d0*cos(kx*x)*sin(ky*y)*ft

    sxy   =    0.6d0*126d0*cos(kx*x)*cos(ky*y)*ft
    sxy_t =    0.6d0*126d0*cos(kx*x)*cos(ky*y)*ft_t
    sxy_x =-kx*0.6d0*126d0*sin(kx*x)*cos(ky*y)*ft
    sxy_y =-ky*0.6d0*126d0*cos(kx*x)*sin(ky*y)*ft

    sxx   = (2d0*side-3d0)*   0.6d0*1260*cos(kx*x)*cos(ky*y)*ft
    sxx_t = (2d0*side-3d0)*   0.6d0*1260*cos(kx*x)*cos(ky*y)*ft_t
    sxx_x =-(2d0*side-3d0)*kx*0.6d0*1260*sin(kx*x)*cos(ky*y)*ft
    sxx_y =-(2d0*side-3d0)*ky*0.6d0*1260*cos(kx*x)*sin(ky*y)*ft

    select case(field)
    case('szz')
        F = 0d0
        return
    case('V')
        F = 2d0*vx
        return
    case('vx')
        F = ny*vx-my*vy
        return
    case('vy')
        F = mx*vy-nx*vx
        return
    case('S')
      F = sxy
      return
    case('N')
      F =-syy
      return
    case('sxx')
      F = syy*my**2+sxx*ny**2-2d0*sxy*my*ny
      return
    case('syy')
      F = syy*mx**2+sxx*nx**2-2d0*sxy*mx*nx
      return
    case('sxy')
      F =-syy*my*mx-sxx*ny*nx+sxy*(ny*mx+nx*my)
      return
    case('Vt')
         F = 2d0*vx*ft_t
         return
    case('St')
         F = sxy*ft_t
         return
    case('Nt')
         F =  -syy*ft_t
         return
    case('s_vx')
        ! vx_t - (1/rho)*(sxx_x+sxy_y)
        vx_t  = ny*vx_t-my*vy_t
        sxx_x = syy_x*my**2+sxx_x*ny**2-2*sxy_x*my*ny &
              + 2d0*syy*my*my_x+2d0*sxx*ny*ny_x-2*sxy*(my_x*ny+my*ny_x)
        sxy_y = -syy_y*my*mx-sxx_y*ny*nx+sxy_y*(ny*mx+nx*my) &
              -syy*(my_y*mx+my*mx_y)-sxx*(ny_y*nx+ny*nx_y)+sxy*(ny_y*mx+ny*mx_y+nx_y*my+nx*my_y)
        F = vx_t - (1d0/rho)*(sxx_x+sxy_y)
        return
    case('s_vy')
        ! vy_t - (1/rho)*(sxy_x+syy_y)
        vy_t = mx*vy_t-nx*vx_t
        sxy_x = -syy_x*my*mx-sxx_x*ny*nx+sxy_x*(ny*mx+nx*my) &
              -syy*(my_x*mx+my*mx_x)-sxx*(ny_x*nx+ny*nx_x)+sxy*(ny_x*mx+ny*mx_x+nx_x*my+nx*my_x)
        syy_y = syy_y*mx**2+sxx_y*nx**2-2*sxy_y*mx*nx &
              + 2d0*syy*mx*mx_y+2d0*sxx*nx*nx_y-2*sxy*(mx_y*nx+mx*nx_y)
        F = vy_t - (1d0/rho)*(sxy_x+syy_y)
        return
    case('s_sxx')
        ! sxx_t - (lam+2*G) vx_x - lam * vy_y
        sxx_t = syy_t*my**2+sxx_t*ny**2-2*sxy_t*my*ny
        vx_x = ny_x*vx-my_x*vy + ny*vx_x-my*vy_x
        vy_y = mx_y*vy-nx_y*vx + mx*vy_y-nx*vx_y
        F = sxx_t - (lambda+2d0*G)*vx_x - lambda*vy_y
        return
    case('s_syy')
        ! syy_t - lam * vx_x - (lam+2*G) vy_y
        vx_x = ny_x*vx-my_x*vy + ny*vx_x-my*vy_x
        vy_y = mx_y*vy-nx_y*vx + mx*vy_y-nx*vx_y
        syy_t = syy_t*mx**2+sxx_t*nx**2-2*sxy_t*mx*nx
        F = syy_t - lambda*vx_x - (lambda+2d0*G)*vy_y
        return
    case('s_sxy')
        ! sxy_t - G * vy_x + G * vx_y
        vx_y = ny_y*vx-my_y*vy + ny*vx_y-my*vy_y
        vy_x = mx_x*vy-nx_x*vx + mx*vy_x-nx*vx_x
        sxy_t =-syy_t*my*mx-sxx_t*ny*nx+sxy_t*(ny*mx+nx*my)
        F = sxy_t - G*(vy_x+vx_y)
        return
    case('s_szz')
        ! sxx_t - (lam+2*G) vx_x - lam * vy_y
        szz_t = 0d0
        vx_x = ny_x*vx-my_x*vy + ny*vx_x-my*vy_x
        vy_y = mx_y*vy-nx_y*vx + mx*vy_y-nx*vx_y
        F = szz_t - lambda*vx_x - lambda*vy_y
        return
    end select

    call error('Invalid field (' // trim(field) // ')')

    F = 0d0

  end function mms_sin

  function mms_simple(x,y,t,side,field,M) result(F)

    use ifport ! for Intel Fortran compiler
    use material, only : block_material

    implicit none

    real,intent(in) :: x,y,t
    integer,intent(in) :: side
    character(*),intent(in) :: field
    type(block_material),intent(in),optional :: M
    real :: F

    real :: r2,syy,sxy,sxx,vx,vy,g,gp,pi
    real :: vx_t,vx_x,vx_y
    real :: vy_t,vy_x,vy_y
    real :: syy_t,syy_x,syy_y
    real :: sxy_t,sxy_x,sxy_y
    real :: sxx_t,sxx_x,sxx_y
    real :: w,k

    ! WARNING: This must match the coordinate transform!!!!!
    pi = 4d0 * datan(1d0)
    w = 2d0*pi
    k = 4d0*pi

    select case(field)
    case('vx','vy','V')
        ! need to be careful since slip and sxy must have the same sign at the
        ! fault for the friction law
        vx = cos(w*t)*cos(k*x)*cos(k*y)
        vy = cos(w*t)*cos(k*x)*cos(k*y)
        select case(field)
        case('V')
            F = 2*vx
            return
        case('vx')
            F = vx
            return
        case('vy')
            F = vy
            return
        end select
    case('sxx','syy','sxy','N','S')
        syy = cos(w*t)*cos(k*x)*cos(k*y)
        sxy = cos(w*t)*cos(k*x)*cos(k*y)
        sxx = cos(w*t)*cos(k*x)*cos(k*y)
        select case(field)
        case('S')
            F = sxy
            return
        case('N')
            F = syy
            return
        case('sxx')
            F = sxx
            return
        case('syy')
            F = syy
            return
        case('sxy')
            F = sxy
            return
        end select
    case('s_vx','s_vy')
        ! need to be careful since slip and sxy must have the same sign at the
        ! fault for the friction law

        ! Need the time derivatives of velocity
        vx_t = -w*sin(w*t)*cos(k*x)*cos(k*y)
        vy_t = -w*sin(w*t)*cos(k*x)*cos(k*y)

        ! For the stresses we need both parts
        syy_x = -k*cos(w*t)*sin(k*x)*cos(k*y)
        sxy_x = -k*cos(w*t)*sin(k*x)*cos(k*y)
        sxx_x = -k*cos(w*t)*sin(k*x)*cos(k*y)

        syy_y = -k*cos(w*t)*cos(k*x)*sin(k*y)
        sxy_y = -k*cos(w*t)*cos(k*x)*sin(k*y)
        sxx_y = -k*cos(w*t)*cos(k*x)*sin(k*y)

        select case(field)
        case('s_vx')
            ! vx_t - (1/rho)*(sxx_x+sxy_y)
            F = vx_t - (1/M%rho)*(sxx_x+sxy_y)
            return
        case('s_vy')
            ! vy_t - (1/rho)*(sxy_x+syy_y)
            F = vy_t - (1/M%rho)*(sxy_x+syy_y)
            return
        end select

    case('s_sxx','s_syy','s_sxy')
        ! Need the time derivatives of stresses
        syy_t = -w*sin(w*t)*cos(k*x)*cos(k*y)
        sxy_t = -w*sin(w*t)*cos(k*x)*cos(k*y)
        sxx_t = -w*sin(w*t)*cos(k*x)*cos(k*y)

        ! For the velocities we need both parts
        vx_x = -k*cos(w*t)*sin(k*x)*cos(k*y)
        vy_x = -k*cos(w*t)*sin(k*x)*cos(k*y)

        vx_y = -k*cos(w*t)*cos(k*x)*sin(k*y)
        vy_y = -k*cos(w*t)*cos(k*x)*sin(k*y)

        select case(field)
        ! For the velocities we will use the chain rule and need all of these
        case('s_sxx')
            ! sxx_t - (lam+2*G) vx_x - lam * vy_y
            F = sxx_t - (M%lambda+2*M%G)*vx_x - M%lambda*vy_y
            return
        case('s_syy')
            ! syy_t - lam * vx_x - (lam+2*G) vy_y
            F = syy_t - (M%lambda+2*M%G)*vy_y - M%lambda*vx_x
            return
        case('s_sxy')
            ! sxy_t - G * vy_x + G * vx_y
            F = sxy_t - M%G*(vy_x+vx_y)
            return
        end select
    end select

    F = 0

  end function mms_simple

  subroutine scale_rates_displacement(B,BF,A)

    use grid, only : block_grid

    type(block_grid),intent(in) :: B
    type(block_fields),intent(inout) :: BF
    real,intent(in) :: A

    if (.not.BF%displacement) return

    if (B%sideL) BF%bndFL%DU = A*BF%bndFL%DU
    if (B%sideR) BF%bndFR%DU = A*BF%bndFR%DU
    if (B%sideB) BF%bndFB%DU = A*BF%bndFB%DU
    if (B%sideT) BF%bndFT%DU = A*BF%bndFT%DU

  end subroutine scale_rates_displacement

end module fields
