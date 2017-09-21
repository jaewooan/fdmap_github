module checkpoint

  implicit none

  ! CKPT = enable checkpointing
  ! RESUME = 
  ! INTERVAL = time step interval between successive checkpoints

  type :: checkpoint_type
     logical :: ckpt,resume,del_final
     character(256) :: filename
     integer :: unit,begin,interval,current,previous
     real :: max_time
  end type checkpoint_type

  ! Checkpoints are written to files NAME*.ckptNNN, where NNN is time step.
  ! After writing new checkpoint, previous checkpoint is deleted.
  ! In addition to writing checkpoints at regular intervals, a checkpoint is
  ! written when aborting. There are two ways to initiate an abort:
  ! (1) Exceed specified wall clock limit (MAX_TIME);
  ! (2) Abort file NAME.abort contains "T" (logical true).
  ! There are two ways to start from a checkpoint:
  ! (1) Resume from last checkpoint (whose number is stored in file NAME.ckpt)
  ! (2) Explicitly specify the checkpoint number
  ! Method (1) supersedes (2) if both are requested.

contains


  subroutine init_checkpoint(name,C,input,nstart,adjoint)

    use io, only : message,error,new_io_unit
    use mpi_routines, only : is_master,master,MPI_REAL_PW
    use mpi

    implicit none
    
    character(*),intent(in) :: name
    integer,intent(in) :: input
    type(checkpoint_type),intent(out) :: C
    integer,intent(out) :: nstart
    logical,intent(out) :: adjoint

    logical :: ckpt,resume,del_final
    integer :: stat,begin,interval,ierr
    real :: max_time

    namelist /checkpoint_list/ ckpt,resume,begin,interval,max_time,del_final,adjoint

    ! defaults

    ckpt = .false.
    resume = .false.
    del_final = .true.
    adjoint = .false.
    begin = 0
    interval = 1
    max_time = 1d40

    rewind(input)
    read(input,nml=checkpoint_list,iostat=stat)
    if (stat>0) call error('Error in checkpoint_list','init_checkpoint')

    C%ckpt = ckpt
    C%resume = resume
    C%del_final = del_final
    C%begin = begin
    C%interval = interval
    C%max_time = max_time
!    C%adjoint = adjoint

    ! open checkpoint number file (master only)
    ! if resuming from last checkpoint, determine its number (supersedes value from input file)

    if (is_master) then

       C%unit = new_io_unit()
       C%filename = trim(name) // '.ckpt'

       if (C%resume) then
          open(C%unit,file=C%filename,iostat=stat,status='old')
          if (stat==0) then
             ! file exists, read checkpoint number
             read(C%unit,*) C%begin
             if (C%begin==-1) then
                call message('Problem successfully completed, no need to restart.')
                call error('Delete ' // trim(C%filename) // ' to solve problem again')
             end if
          else
             ! file does not exist, so create file and set contents to 0
             open(C%unit,file=C%filename,status='new') 
             write(C%unit,*) 0
             rewind(C%unit)
          end if
       else
          ! create file, overwriting if needed
          open(C%unit,file=C%filename,iostat=stat,status='replace')
          if (stat/=0) call error('Error opening ' // trim(C%filename),'init_checkpoint')
       end if
       
    end if
    
    call MPI_Bcast(C%begin,1,MPI_REAL_PW,master,MPI_COMM_WORLD,ierr)
        
    C%begin = max(C%begin,0)
    C%previous = C%begin
    C%current = C%begin

    ! set nstart (initial time step, starting from 1 instead of 0)
    
    nstart = C%begin+1

  end subroutine init_checkpoint


  function abort_now(name,C) result(abort)

    use io, only : new_io_unit
    use mpi_routines, only : time_elapsed

    implicit none

    character(*),intent(in) :: name
    type(checkpoint_type),intent(in) :: C
    logical :: abort

    integer :: unit,stat

    abort = .false.

    ! examine contents of abort file (is one exists)

    unit = new_io_unit()
    open(unit,file=trim(name) // '.abort',iostat=stat,status='old')
    if (stat==0) then
       read(unit,'(l1)',end=100) abort
       close(unit)
    end if
100 continue

    ! check if elapsed time exceeds limit (with max_time in minutes)
    
    abort = (abort.or.(time_elapsed()>=C%max_time*60d0))

  end function abort_now


  subroutine checkpoint_read(name,C,D,adjoint,mode)

    use domain, only : domain_type
    use io, only : message
    use mpi_routines, only : is_master

    implicit none

    character(*),intent(in) :: name
    type(checkpoint_type),intent(in) :: C
    type(domain_type),intent(inout) :: D
    logical,intent(in) :: adjoint
    integer,intent(in) :: mode

    character(128) :: str

    if (.not.C%ckpt) return
    if (C%begin==0) return

    ! read checkpoint at C%begin

    write(str,'(a,i0)') '...reading  checkpoint at time step ',C%begin
    if (is_master) call message(str)

    call checkpoint_op('read',name,C%begin,D,adjoint,mode)

  end subroutine checkpoint_read


  subroutine checkpoint_write_delete(name,number,C,D,abort,adjoint,mode)

    use domain, only : domain_type
    use io, only : message
    use mpi_routines, only : is_master

    implicit none

    character(*),intent(in) :: name
    integer,intent(in) :: number,mode
    type(checkpoint_type),intent(inout) :: C
    type(domain_type),intent(inout) :: D
    logical,intent(in) :: abort,adjoint

    character(256) :: str

    if (abort) then
       C%ckpt = .true.
    else
       if (mod(number,C%interval)/=0) return
    end if

    if (.not.C%ckpt) return

    write(str,'(a,i0)') '...writing  checkpoint at time step ',number
    if (is_master) call message(str)
    
    call checkpoint_op('write',name,number,D,adjoint,mode)

    if (is_master) then
       rewind(C%unit)
       write(C%unit,*) number
    end if

    call checkpoint_delete(name,C,D)

    C%previous = number

  end subroutine checkpoint_write_delete


  subroutine checkpoint_delete(name,C,D,final,adjoint,mode)

    use domain, only : domain_type
    use io, only : message
    use mpi_routines, only : is_master

    implicit none

    character(*),intent(in) :: name
    type(checkpoint_type),intent(in) :: C
    type(domain_type),intent(inout) :: D
    logical,intent(in),optional :: final,adjoint
    integer,intent(in),optional :: mode

    character(256) :: str

    if (.not.C%ckpt) return
    if (C%previous==0) return

    ! delete checkpoint at C%previous

    write(str,'(a,i0)') '...deleting checkpoint at time step ',C%previous
    if (is_master) call message(str)

    call checkpoint_op('delete',name,C%previous,D,adjoint,mode)

    if (present(final)) then
       if (final.and.is_master) then
          rewind(C%unit)
          write(C%unit,*) -1
       end if
    end if

  end subroutine checkpoint_delete


  subroutine checkpoint_op(operation,name,number,D,adjoint,mode)

    use domain, only : domain_type
    use interfaces, only : checkpoint_iface
    use fields, only : checkpoint_fields,checkpoint_block_fields,exchange_fields
    use io, only : error

    implicit none

    character(*),intent(in) :: operation,name
    integer,intent(in) :: number,mode
    type(domain_type),intent(inout) :: D
    logical,intent(in) :: adjoint

    integer :: i

    select case(operation)
    case('read','write','delete')
       call checkpoint_fields(operation,name,number,D%C,D%F,adjoint,mode)
       do i = 1,D%nblocks
          call checkpoint_block_fields(operation,name,number,i,D%B(i)%G,D%B(i)%F,D%F)
       end do
       do i = 1,D%nifaces
          call checkpoint_iface(operation,name,number,i,D%I(i))
       end do
    case default
       call error('Invalid checkpoint operation','checkpoint_op')
    end select

    select case(operation)
    case('read')
       call exchange_fields(D%C,D%F)
    end select

  end subroutine checkpoint_op


  subroutine finish_checkpoint(C)
    
    use mpi_routines, only : is_master

    implicit none

    type(checkpoint_type),intent(inout) :: C

    if (is_master) close(C%unit)

  end subroutine finish_checkpoint


end module checkpoint
