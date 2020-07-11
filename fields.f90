module fields

  implicit none

  type :: bnd_fields
     ! F = fields (grid data)
     ! FHAT = fields (hat variables)
     ! DF = field rates
     ! F0 = initial fields
     ! M = material properties
     ! DE = rate of energy flow through boundary
     ! (integrate over boundary to get energy dissipation rate)
     ! E = total energy lost through boundary (DE integrated over time)
     ! U = boundary displacement
     ! DU = boundary displacement rate
     logical :: holds_bnd
     integer :: comm
     real,dimension(:,:),allocatable :: U,DU,F,Fhat,F0,DF,M
     real,dimension(:),allocatable :: E,DE
  end type bnd_fields

  type :: block_fields
     logical :: energy_balance
     real :: F0(9),DE(4),E(4),Etot
     real,dimension(:),allocatable :: Hx,Hy
     type(bnd_fields) :: bndFL,bndFR,bndFB,bndFT
  end type block_fields

  type :: fields_type
     character(256) :: problem,prestress_filename
     logical :: displacement,energy_balance,peak,prestress_from_file
     integer :: nF,nC,nU,ns,nEP
     real :: Etot
     real,dimension(:,:,:),allocatable :: F,U,DU,DF,EP,pgv,pga,S0
     real,dimension(:,:,:),allocatable :: Wx,Wy,DWx,DWy ! PML fields
     real,dimension(:,:),allocatable :: gammap,lambda,Wp
     real,dimension(:),allocatable :: Hx,Hy
  end type fields_type

  type :: fields_perturb
    character(256) :: shape
    real :: x0,y0,Lx,Ly,vx,vy,vz,sxx,sxy,sxz,syy,syz,szz
  end type fields_perturb

  interface rotate_fields_xy2nt
     module procedure rotate_fields_xy2nt_mode3,rotate_fields_xy2nt_mode2
  end interface

  interface rotate_fields_nt2xy
     module procedure rotate_fields_nt2xy_mode3,rotate_fields_nt2xy_mode2
  end interface


contains


  subroutine init_fields(mode,iblock,C,B,G,BF,F,t,input, &
      energy_balance,displacement,peak,M,pmlx,pmly)

    use mpi_routines, only : new_communicator
    use mpi_routines2d, only : cartesian,allocate_array_body
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
    logical,intent(in) :: energy_balance,displacement,pmlx,pmly,peak

    integer :: stat
    logical :: prestress_from_file
    real :: Psi,vx0,vy0,vz0,sxx0,sxy0,sxz0,syy0,syz0,szz0
    character(256) :: problem,prestress_filename
    character(256) :: str
    ! P5 does not influence boundary conditions only interior
    type(fields_perturb) :: P1,P2,P3,P4,P5
    real,dimension(:),allocatable :: F0
    type(limits) :: lim

    namelist /fields_list/ problem,Psi,vx0,vy0,vz0, &
         sxx0,sxy0,sxz0,syy0,syz0,szz0, &
         P1,P2,P3,P4,P5, &
         prestress_from_file,prestress_filename

    ! defaults

    problem = ''

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
    P5 = fields_perturb('',0d0,0d0,1d0,1d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0)

    prestress_from_file = .false.
    prestress_filename = ''
    
    ! read in field parameters

    write(str,'(a,i0,a)') '!---BLOCK',iblock,'---'
    call seek_to_string(input,str)
    read(input,nml=fields_list,iostat=stat)
    if (stat>0) call error('Error in fields_list','init_fields')

    F%problem = problem

    F%energy_balance = energy_balance
    BF%energy_balance = energy_balance
    F%displacement = displacement
    F%peak = peak
    
    select case(F%problem)
    case default
    case('uniform') ! initial velocity = 0 => uniform state variable (set elsewhere)
       ! syy0 and sxy0 given, set sxx0 using Psi
       if (Psi==45d0) then
          sxx0 = syy0 ! special case, avoid tan(pi/2)=Inf
       else
          sxx0 = (1d0-2d0*sxy0/(syy0*tan(2d0*deg2rad(Psi))))*syy0
       end if
       ! set szz0 as average of sxx0 and syy0
       szz0 = 0.5d0*(sxx0+syy0)
    end select

    F%Etot = 1d40

    ! store initial fields

    BF%F0 = (/ vx0,vy0,vz0,sxx0,sxy0,sxz0,syy0,syz0,szz0 /)

    select case(mode)
    case(2)
       F%nF = 6
       F%nEP = 6
       F%nC = 5
       F%nU = 2
       F%ns = 4
       allocate(F0(F%nF))
       F0 = (/ vx0,vy0,sxx0,sxy0,syy0,szz0 /)
    case(3)
       F%nF = 3
       F%nEP = 6
       F%nC = 3
       F%nU = 1
       F%ns = 2
       allocate(F0(F%nF))
       F0 = (/ vz0,sxz0,syz0 /)
    end select

    ! allocate memory for fields

    if (.not.allocated(F%F)) then
       call allocate_array_body(F%F ,C,F%nF,ghost_nodes=.true.)
       call allocate_array_body(F%DF,C,F%nF,ghost_nodes=.false.)
       if(F%displacement) then
          call allocate_array_body(F%U ,C,F%nU,ghost_nodes=.false.,Fval=0d0)
          call allocate_array_body(F%DU,C,F%nU,ghost_nodes=.false.)
       end if
       if(F%peak) then
          call allocate_array_body(F%pgv,C,1,ghost_nodes=.false.,Fval=0d0)
          call allocate_array_body(F%pga,C,1,ghost_nodes=.false.,Fval=0d0)
       end if
       allocate(F%Hx(C%mx:C%px),F%Hy(C%my:C%py))
       F%Hx = 1d40
       F%Hy = 1d40
       if (M%response=='plastic') then
          call allocate_array_body(F%gammap ,C,ghost_nodes=.false.,Fval=0d0)
          call allocate_array_body(F%lambda ,C,ghost_nodes=.false.,Fval=0d0)
          if(M%plastic_strain_tensor) &
               call allocate_array_body(F%Ep,C,F%nEP,ghost_nodes=.false.,Fval=0d0)
          if(M%plastic_work) &
               call allocate_array_body(F%Wp,C,ghost_nodes=.false.,Fval=0d0)
       end if
    end if
    if (pmlx .and. .not.allocated(F%Wx)) then
      call allocate_array_body(F%Wx, C,F%nF,ghost_nodes=.false.,Fval=0d0)
      call allocate_array_body(F%DWx,C,F%nF,ghost_nodes=.false.,Fval=0d0)
    end if
    if (pmly .and. .not.allocated(F%Wy)) then
      call allocate_array_body(F%Wy, C,F%nF,ghost_nodes=.false.,Fval=0d0)
      call allocate_array_body(F%DWy,C,F%nF,ghost_nodes=.false.,Fval=0d0)
    end if

    ! initialize prestress, if input from file
    ! This prestress is stored in F%S0 and is only used by plasticity routines;
    ! it is assumed (but never verified) that S0 satisfies equilbrium.
    ! Note that S0 is NOT used by friction routines; those routines require prestress
    ! to be input as separate (fault prestress) file.
    
    if (prestress_from_file) call init_prestress_from_file(F,C,prestress_filename)
    
    ! initialize fields on sides of this block

    call init_fields_side(B%bndL,BF%bndFL,B%skip,B%my,B%py,F%nF,F%nU,F0,mode,problem, &
         t,iblock,P1,P2,P3,P4,M)
    call init_fields_side(B%bndR,BF%bndFR,B%skip,B%my,B%py,F%nF,F%nU,F0,mode,problem, &
         t,iblock,P1,P2,P3,P4,M)
    call init_fields_side(B%bndB,BF%bndFB,B%skip,B%mx,B%px,F%nF,F%nU,F0,mode,problem, &
         t,iblock,P1,P2,P3,P4,M)
    call init_fields_side(B%bndT,BF%bndFT,B%skip,B%mx,B%px,F%nF,F%nU,F0,mode,problem, &
         t,iblock,P1,P2,P3,P4,M)

    ! initialize fields in interior of this block

    call init_fields_interior(B,G,F,F0,mode,problem,t,iblock,P1,P2,P3,P4,P5,M)

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
    if (allocated(F%lambda )) deallocate(F%lambda )
    if (allocated(F%EP     )) deallocate(F%EP     )
    if (allocated(F%Wp     )) deallocate(F%Wp     )
    if (allocated(F%Hx     )) deallocate(F%Hx     )
    if (allocated(F%Hy     )) deallocate(F%Hy     )
    if (allocated(F%S0     )) deallocate(F%S0     )

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

    if (allocated(bndF%F   )) deallocate(bndF%F   )
    if (allocated(bndF%Fhat)) deallocate(bndF%Fhat)
    if (allocated(bndF%DF  )) deallocate(bndF%DF  )
    if (allocated(bndF%F0  )) deallocate(bndF%F0  )
    if (allocated(bndF%M   )) deallocate(bndF%M   )
    if (allocated(bndF%U   )) deallocate(bndF%U   )
    if (allocated(bndF%DU  )) deallocate(bndF%DU  )
    if (allocated(bndF%E   )) deallocate(bndF%E   )
    if (allocated(bndF%DE  )) deallocate(bndF%DE  )

  end subroutine destroy_bnd_fields


  subroutine checkpoint_fields(operation,name,checkpoint_number,C,F,adjoint,mode)

    use io, only : file_distributed,open_file_distributed, &
         read_file_distributed,write_file_distributed,close_file_distributed
    use mpi_routines2d, only : cartesian
    use mpi_routines, only : pw,is_master
    use mpi

    implicit none

    character(*),intent(in) :: operation,name
    integer,intent(in) :: checkpoint_number,mode
    type(cartesian),intent(in) :: C
    type(fields_type),intent(inout) :: F
    logical,intent(in) :: adjoint

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
         if(F%displacement) call read_file_distributed(fh,F%U( C%mx:C%px,C%my:C%py,l))
         if(F%displacement) call read_file_distributed(fh,F%DU(C%mx:C%px,C%my:C%py,l))
       end do
       if(F%peak) call read_file_distributed(fh,F%pgv(C%mx:C%px,C%my:C%py,1))
       if(F%peak) call read_file_distributed(fh,F%pga(C%mx:C%px,C%my:C%py,1))
       if (allocated(F%gammap )) call read_file_distributed(fh,F%gammap )
       if (allocated(F%lambda )) call read_file_distributed(fh,F%lambda )
       if (allocated(F%Wp     )) call read_file_distributed(fh,F%Wp     )
       if (allocated(F%Wx     )) call read_file_distributed(fh,F%Wx     )
       if (allocated(F%Wy     )) call read_file_distributed(fh,F%Wy     )
       if (allocated(F%EP     )) then
         do l = 1,F%nEP
           call read_file_distributed(fh,F%EP(C%mx:C%px,C%my:C%py,l))
         end do
       end if
    case('write')
       do l = 1,F%nF
          call write_file_distributed(fh,F%F(C%mx:C%px,C%my:C%py,l))
       end do
       do l = 1,F%nF
          call write_file_distributed(fh,F%DF(:,:,l))
       end do
       do l = 1,F%nU
         if(F%displacement) call write_file_distributed(fh,F%U( :,:,l))
         if(F%displacement) call write_file_distributed(fh,F%DU(:,:,l))
       end do
       if(F%peak) call write_file_distributed(fh,F%pgv(:,:,1))
       if(F%peak) call write_file_distributed(fh,F%pga(:,:,1))
       if (allocated(F%gammap )) call write_file_distributed(fh,F%gammap )
       if (allocated(F%lambda )) call write_file_distributed(fh,F%lambda )
       if (allocated(F%Wp     )) call write_file_distributed(fh,F%Wp     )
       if (allocated(F%Wx     )) call write_file_distributed(fh,F%Wx     )
       if (allocated(F%Wy     )) call write_file_distributed(fh,F%Wy     )
       if (allocated(F%EP     )) then
         do l = 1,F%nEP
           call write_file_distributed(fh,F%EP(C%mx:C%px,C%my:C%py,l))
         end do
       end if

      if (adjoint) then

       select case(mode)
          case(2)
             F%F(:,:,1:2) = -F%F(:,:,1:2)
          case(3)
             F%F(:,:,1) = -F%F(:,:,1)
          end select

       end if
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
    character(256) :: filename,filenameE
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
                   call read_file_distributed(fh,BF%bndFL%F   (:,l))
                   call read_file_distributed(fh,BF%bndFL%Fhat(:,l))
                   call read_file_distributed(fh,BF%bndFL%F0  (:,l))
                end do
                do l = 1,F%nU
                  call read_file_distributed(fh,BF%bndFL%U( :,l))
                  call read_file_distributed(fh,BF%bndFL%DU(:,l))
                end do
                call read_file_distributed(fh,BF%bndFL%E)
                call read_file_distributed(fh,BF%bndFL%DE)
             case(2) ! R
                do l = 1,F%nF
                   call read_file_distributed(fh,BF%bndFR%F   (:,l))
                   call read_file_distributed(fh,BF%bndFR%Fhat(:,l))
                   call read_file_distributed(fh,BF%bndFR%F0  (:,l))
                end do
                do l = 1,F%nU
                  call read_file_distributed(fh,BF%bndFR%U( :,l))
                  call read_file_distributed(fh,BF%bndFR%DU(:,l))
                end do
                 call read_file_distributed(fh,BF%bndFR%E)
                 call read_file_distributed(fh,BF%bndFR%DE)
             case(3) ! B
                do l = 1,F%nF
                   call read_file_distributed(fh,BF%bndFB%F   (:,l))
                   call read_file_distributed(fh,BF%bndFB%Fhat(:,l))
                   call read_file_distributed(fh,BF%bndFB%F0  (:,l))
                end do
                do l = 1,F%nU
                  call read_file_distributed(fh,BF%bndFB%U( :,l))
                  call read_file_distributed(fh,BF%bndFB%DU(:,l))
                end do
                 call read_file_distributed(fh,BF%bndFB%E)
                 call read_file_distributed(fh,BF%bndFB%DE)
             case(4) ! T
                do l = 1,F%nF
                   call read_file_distributed(fh,BF%bndFT%F   (:,l))
                   call read_file_distributed(fh,BF%bndFT%Fhat(:,l))
                   call read_file_distributed(fh,BF%bndFT%F0  (:,l))
                end do
                do l = 1,F%nU
                  call read_file_distributed(fh,BF%bndFT%U( :,l))
                  call read_file_distributed(fh,BF%bndFT%DU(:,l))
                end do
                 call read_file_distributed(fh,BF%bndFT%E)
                 call read_file_distributed(fh,BF%bndFT%DE)
             end select
          case('write')
             select case(side)
             case(1) ! L
                do l = 1,F%nF
                   call write_file_distributed(fh,BF%bndFL%F   (:,l))
                   call write_file_distributed(fh,BF%bndFL%Fhat(:,l))
                   call write_file_distributed(fh,BF%bndFL%F0  (:,l))
                end do
                do l = 1,F%nU
                  call write_file_distributed(fh,BF%bndFL%U( :,l))
                  call write_file_distributed(fh,BF%bndFL%DU(:,l))
                end do
                 call write_file_distributed(fh,BF%bndFL%E)
                 call write_file_distributed(fh,BF%bndFL%DE)
             case(2) ! R
                do l = 1,F%nF
                   call write_file_distributed(fh,BF%bndFR%F   (:,l))
                   call write_file_distributed(fh,BF%bndFR%Fhat(:,l))
                   call write_file_distributed(fh,BF%bndFR%F0  (:,l))
                end do
                do l = 1,F%nU
                  call write_file_distributed(fh,BF%bndFR%U( :,l))
                  call write_file_distributed(fh,BF%bndFR%DU(:,l))
                end do
                 call write_file_distributed(fh,BF%bndFR%E)
                 call write_file_distributed(fh,BF%bndFR%DE)
             case(3) ! B
                do l = 1,F%nF
                   call write_file_distributed(fh,BF%bndFB%F   (:,l))
                   call write_file_distributed(fh,BF%bndFB%Fhat(:,l))
                   call write_file_distributed(fh,BF%bndFB%F0  (:,l))
                end do
                do l = 1,F%nU
                  call write_file_distributed(fh,BF%bndFB%U( :,l))
                  call write_file_distributed(fh,BF%bndFB%DU(:,l))
                end do
                 call write_file_distributed(fh,BF%bndFB%E)
                 call write_file_distributed(fh,BF%bndFB%DE)
             case(4) ! T
                do l = 1,F%nF
                   call write_file_distributed(fh,BF%bndFT%F   (:,l))
                   call write_file_distributed(fh,BF%bndFT%Fhat(:,l))
                   call write_file_distributed(fh,BF%bndFT%F0  (:,l))
                end do
                do l = 1,F%nU
                  call write_file_distributed(fh,BF%bndFT%U( :,l))
                  call write_file_distributed(fh,BF%bndFT%DU(:,l))
                end do
                 call write_file_distributed(fh,BF%bndFT%E)
                 call write_file_distributed(fh,BF%bndFT%DE)
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

    if(io_process) call MPI_Comm_free(comm,ierr)

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

  end subroutine checkpoint_block_fields


  subroutine init_prestress_from_file(F,C,prestress_filename)

    use mpi_routines2d, only : cartesian,allocate_array_body

    implicit none

    type(fields_type),intent(inout) :: F
    type(cartesian),intent(in) :: C
    character(*) :: prestress_filename

    if (.not.allocated(F%S0)) then
       ! allocate array
       call allocate_array_body(F%S0,C,F%ns,ghost_nodes=.false.)
       ! read values from file
       call prestressIO('read',prestress_filename,C,F)
    end if
       
  end subroutine init_prestress_from_file

  subroutine init_force_term_from_file(F,C,prestress_filename)
    ! Abrahams Edit Edit
  end subroutine init_force_term_from_file



  subroutine prestressIO(operation,filename,C,F)

    use io, only : file_distributed,open_file_distributed, &
         read_file_distributed,write_file_distributed,close_file_distributed
    use mpi_routines2d, only : cartesian
    use mpi_routines, only : pw,is_master
    use mpi

    implicit none

    character(*),intent(in) :: operation,filename
    type(cartesian),intent(in) :: C
    type(fields_type),intent(inout) :: F

    type(file_distributed) :: fh
    integer :: l,ierr

    if (operation=='delete') then
       if (is_master) call MPI_file_delete(filename,MPI_INFO_NULL,ierr)
       return
    end if

    call open_file_distributed(fh,filename,operation,C%c2d%comm,C%c2d%array_w,pw)

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    select case(operation)
    case('read')
       do l = 1,F%ns
          call  read_file_distributed(fh,F%S0(C%mx:C%px,C%my:C%py,l))
       end do
    case('write')
       do l = 1,F%ns
          call write_file_distributed(fh,F%S0(C%mx:C%px,C%my:C%py,l))
       end do
    end select

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    call close_file_distributed(fh)

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

  end subroutine prestressIO
  
  
  subroutine init_fields_interior(B,G,F,F0,mode,problem,t,iblock,P1,P2,P3,P4,P5,M)

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
    type(fields_perturb),intent(in) :: P1,P2,P3,P4,P5

    integer :: i,j
    real :: A1,A2,A3,A4,A5
    real,dimension(F%nF) :: F1,F2,F3,F4,F5

    if (B%skip) return ! no points in this block

    select case(mode)
    case(2)
       F1 = (/ P1%vx,P1%vy,P1%sxx,P1%sxy,P1%syy,P1%szz /)
       F2 = (/ P2%vx,P2%vy,P2%sxx,P2%sxy,P2%syy,P2%szz /)
       F3 = (/ P3%vx,P3%vy,P3%sxx,P3%sxy,P3%syy,P3%szz /)
       F4 = (/ P4%vx,P4%vy,P4%sxx,P4%sxy,P4%syy,P4%szz /)
       F5 = (/ P5%vx,P5%vy,P5%sxx,P5%sxy,P5%syy,P5%szz /)
    case(3)
       F1 = (/ P1%vz,P1%sxz,P1%syz /)
       F2 = (/ P2%vz,P2%sxz,P2%syz /)
       F3 = (/ P3%vz,P3%sxz,P3%syz /)
       F4 = (/ P4%vz,P4%sxz,P4%syz /)
       F5 = (/ P5%vz,P5%sxz,P5%syz /)
    end select

    do j = B%my,B%py
       do i = B%mx,B%px
          A1 = perturb_fields(G%x(i,j),G%y(i,j),P1)
          A2 = perturb_fields(G%x(i,j),G%y(i,j),P2)
          A3 = perturb_fields(G%x(i,j),G%y(i,j),P3)
          A4 = perturb_fields(G%x(i,j),G%y(i,j),P4)
          A5 = perturb_fields(G%x(i,j),G%y(i,j),P5)
          F%F(i,j,:) = F0+A1*F1+A2*F2+A3*F3+A4*F4+A5*F5
          call initial_fields(F%F(i,j,:),G%x(i,j),G%y(i,j),t,problem,iblock,M)
       end do
    end do

  end subroutine init_fields_interior


  subroutine init_fields_side(bnd,bndF,Bskip,m,p,nF,nU,F0,mode,problem,t,iblock,&
      P1,P2,P3,P4,Mat)

    use geometry, only : curve
    use material, only : block_material

    implicit none

    type(curve),intent(in) :: bnd
    type(bnd_fields),intent(inout) :: bndF
    logical,intent(in) :: Bskip
    integer,intent(in) :: m,p,nF,nU
    real,intent(in) :: F0(:),t
    integer,intent(in) :: mode,iblock
    character(*),intent(in) :: problem
    type(block_material),intent(in) :: Mat
    type(fields_perturb),intent(in) :: P1,P2,P3,P4

    integer :: i
    real :: A1,A2,A3,A4
    real,dimension(nF) :: F1,F2,F3,F4

    allocate(bndF%F(m:p,nF),bndF%Fhat(m:p,nF),bndF%F0(m:p,nF),bndF%M(m:p,5))
    bndF%F    = 1d40
    bndF%Fhat = 1d40
    bndF%F0   = 1d40
    bndF%M    = 1d40

    allocate(bndF%U( m:p,nU),bndF%DU(m:p,nU))
    bndF%U  = 0d0
    bndF%DU  = 1d40

    allocate(bndF%E(m:p),bndF%DE(m:p))
    bndF%E  = 0d0
    bndF%DE = 1d40

    ! note that boundary fields are allocated (but not initialized with proper values)
    ! even if process is not responsible for this block -- this is because fields on
    ! opposite side of interface (for which another process is responsible) may be needed
    ! to enforce interface conditions (values are set during exchange with other process)

    if (Bskip) return

    ! initial fields on boundaries

    ! store amplitudes for field perturbations

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

    ! then set possibly spatially variable fields

    do i = m,p
       A1 = perturb_fields(bnd%x(i),bnd%y(i),P1)
       A2 = perturb_fields(bnd%x(i),bnd%y(i),P2)
       A3 = perturb_fields(bnd%x(i),bnd%y(i),P3)
       A4 = perturb_fields(bnd%x(i),bnd%y(i),P4)
       bndF%F0(i,:) = F0+A1*F1+A2*F2+A3*F3+A4*F4
       call initial_fields(bndF%F0(i,:),bnd%x(i),bnd%y(i),t,problem,iblock,Mat)
    end do

  end subroutine init_fields_side


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
    case('dy_gaussian')
       A = -(y-P%y0)/(P%Ly**2)*exp(-0.5*(y - P%y0)**2/(P%Ly**2))
    end select

  end function perturb_fields


  subroutine initial_fields(F,x,y,t,problem,iblock,M)

    use utilities, only : step
    use material, only : block_material
    use mms, only : inplane_fault_mms,mms_sin,mms_simple,mms_hydrofrac

    implicit none

    real,dimension(:),intent(inout) :: F
    real,intent(in) :: x,y,t
    character(*),intent(in) :: problem
    integer,intent(in) :: iblock
    type(block_material), intent(in) :: M

    select case(problem)
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
    case('mms-simple')
       ! verification using method of manufactured solutions
       ! for curvilinear
       ! mode II only
       F(1) = F(1)+mms_simple(x,y,t,iblock,'vx')
       F(2) = F(2)+mms_simple(x,y,t,iblock,'vy')
       F(3) = F(3)+mms_simple(x,y,t,iblock,'sxx')
       F(4) = F(4)+mms_simple(x,y,t,iblock,'sxy')
       F(5) = F(5)+mms_simple(x,y,t,iblock,'syy')
       F(6) = F(6)+mms_simple(x,y,t,iblock,'szz')
    case('mms-hydrofrac')
       ! verification using method of manufactured solutions
       ! for curvilinear
       ! mode II only
       F(1) = F(1)+mms_hydrofrac(x,y,t,iblock,'vx')
       F(2) = F(2)+mms_hydrofrac(x,y,t,iblock,'vy')
       F(3) = F(3)+mms_hydrofrac(x,y,t,iblock,'sxx')
       F(4) = F(4)+mms_hydrofrac(x,y,t,iblock,'sxy')
       F(5) = F(5)+mms_hydrofrac(x,y,t,iblock,'syy')
       F(6) = F(6)+mms_hydrofrac(x,y,t,iblock,'szz')
    end select

  end subroutine initial_fields


  subroutine initial_stress(s,x,y,s0,problem,mode,s0_file)

    ! used only by plasticity routines

    use utilities, only : step

    implicit none

    real,intent(out) :: s(6)
    real,intent(in) :: x,y,s0(6)
    character(*),intent(in) :: problem
    integer,intent(in) :: mode
    real,intent(in),optional :: s0_file(:)

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
    end select

    ! add prestress that was input in file
    if (present(s0_file)) then
       select case(mode)
       case(2)
          s(1) = s(1)+s0_file(1) ! sxx
          s(2) = s(2)+s0_file(2) ! sxy
          s(4) = s(4)+s0_file(3) ! syy
          s(6) = s(6)+s0_file(4) ! szz
       case(3)
          s(3) = s(3)+s0_file(1) ! sxz
          s(5) = s(5)+s0_file(2) ! syz
       end select
    end if
    
  end subroutine initial_stress


  subroutine scale_rates_interior(C,F,A)

    use mpi_routines2d, only : cartesian

    implicit none

    type(cartesian),intent(in) :: C
    type(fields_type),intent(inout) :: F
    real,intent(in) :: A

    call   sri3d(C%mx,C%my,1,F%DF ,A,C%mx,C%px,C%my,C%py,1,F%nF)
    if(allocated(F%DU)) &
      call sri3d(C%mx,C%my,1,F%DU ,A,C%mx,C%px,C%my,C%py,1,F%nU)
    if(allocated(F%Wx)) &
      call sri3d(C%mx,C%my,1,F%DWx,A,C%mx,C%px,C%my,C%py,1,F%nF)
    if(allocated(F%Wy)) &
      call sri3d(C%mx,C%my,1,F%DWy,A,C%mx,C%px,C%my,C%py,1,F%nF)

  end subroutine scale_rates_interior


  subroutine scale_rates_boundary(B,BF,A)

    use grid, only : block_grid

    implicit none

    type(block_grid),intent(in) :: B
    type(block_fields),intent(inout) :: BF
    real,intent(in) :: A

    if (B%sideL) BF%bndFL%DU = A*BF%bndFL%DU
    if (B%sideR) BF%bndFR%DU = A*BF%bndFR%DU
    if (B%sideB) BF%bndFB%DU = A*BF%bndFB%DU
    if (B%sideT) BF%bndFT%DU = A*BF%bndFT%DU

  end subroutine scale_rates_boundary


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


  subroutine update_fields_peak(C,F,mode)

    use mpi_routines2d, only : cartesian

    implicit none

    type(cartesian),intent(in) :: C
    type(fields_type),intent(inout) :: F
    integer,intent(in) :: mode
    integer :: i,j

    select case(mode)
    case(2)
      do j = C%my,C%py
        do i = C%mx,C%px
          F%pgv(i,j,1) = max(F%pgv(i,j,1),sqrt(F%F( i,j,1)**2 + F%F( i,j,2)**2))
          F%pga(i,j,1) = max(F%pga(i,j,1),sqrt(F%DF(i,j,1)**2 + F%DF(i,j,2)**2))
        end do
      end do
    case(3)
      do j = C%my,C%py
        do i = C%mx,C%px
          F%pgv(i,j,1) = max(F%pgv(i,j,1),abs(F%F( i,j,1)))
          F%pga(i,j,1) = max(F%pga(i,j,1),abs(F%DF(i,j,1)))
        end do
      end do
    end select

  end subroutine update_fields_peak


  subroutine update_fields_interior(C,F,dt)

    use mpi_routines2d, only : cartesian

    implicit none

    type(cartesian),intent(in) :: C
    type(fields_type),intent(inout) :: F
    real,intent(in) :: dt

    call ufi3d(C%mbx,C%mby,1,C%mx,C%my,1,F%F,F%DF,dt,C%mx,C%px,C%my,C%py,1,F%nF)

    if (F%displacement) &
         call ufi3d(C%mx,C%my,1,C%mx,C%my,1,F%U ,F%DU ,dt,C%mx,C%px,C%my,C%py,1,F%nU)
    if (allocated(F%Wx)) &
         call ufi3d(C%mx,C%my,1,C%mx,C%my,1,F%Wx,F%DWx,dt,C%mx,C%px,C%my,C%py,1,F%nF)
    if (allocated(F%Wy)) &
         call ufi3d(C%mx,C%my,1,C%mx,C%my,1,F%Wy,F%DWy,dt,C%mx,C%px,C%my,C%py,1,F%nF)
    
  end subroutine update_fields_interior


  subroutine update_fields_boundary(B,BF,dt)

    use grid, only : block_grid

    implicit none

    type(block_grid),intent(in) :: B
    type(block_fields),intent(inout) :: BF
    real,intent(in) :: dt

    if (B%sideL) BF%bndFL%U = BF%bndFL%U+dt*BF%bndFL%DU
    if (B%sideR) BF%bndFR%U = BF%bndFR%U+dt*BF%bndFR%DU
    if (B%sideB) BF%bndFB%U = BF%bndFB%U+dt*BF%bndFB%DU
    if (B%sideT) BF%bndFT%U = BF%bndFT%U+dt*BF%bndFT%DU

  end subroutine update_fields_boundary


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


  subroutine set_rates_displacement(C,F)

    use mpi_routines2d, only : cartesian

    implicit none

    type(cartesian),intent(in) :: C
    type(fields_type),intent(inout) :: F

    ! DU = DU+F(1:nU)
    if (F%displacement) &
         call ufi3d(C%mx,C%my,1,C%mbx,C%mby,1,F%DU,F%F,1d0,C%mx,C%px,C%my,C%py,1,F%nU)

  end subroutine set_rates_displacement


  subroutine set_rates_boundary(B,BF,nU)

    use grid, only : block_grid

    implicit none

    type(block_grid),intent(in) :: B
    type(block_fields),intent(inout) :: BF
    integer,intent(in) :: nU

    if (B%sideL) BF%bndFL%DU = BF%bndFL%DU+BF%bndFL%Fhat(:,1:nU)
    if (B%sideR) BF%bndFR%DU = BF%bndFR%DU+BF%bndFR%Fhat(:,1:nU)
    if (B%sideB) BF%bndFB%DU = BF%bndFB%DU+BF%bndFB%Fhat(:,1:nU)
    if (B%sideT) BF%bndFT%DU = BF%bndFT%DU+BF%bndFT%Fhat(:,1:nU)

  end subroutine set_rates_boundary


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


  subroutine scale_rates_displacement(B,BF,A)

    use grid, only : block_grid

    type(block_grid),intent(in) :: B
    type(block_fields),intent(inout) :: BF
    real,intent(in) :: A

    if (B%sideL) BF%bndFL%DU = A*BF%bndFL%DU
    if (B%sideR) BF%bndFR%DU = A*BF%bndFR%DU
    if (B%sideB) BF%bndFB%DU = A*BF%bndFB%DU
    if (B%sideT) BF%bndFT%DU = A*BF%bndFT%DU

  end subroutine scale_rates_displacement


end module fields
