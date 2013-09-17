module io

  implicit none

  logical,parameter :: collective=.true.

  type :: file_distributed
     character(256) :: name
     integer :: fh,array,pr,MPI_REAL_PR
  end type file_distributed

  interface write_file_distributed
     module procedure write_file_distributed_0d,write_file_distributed_1d, &
          write_file_distributed_2d,write_file_distributed_3d
  end interface

  interface read_file_distributed
     module procedure read_file_distributed_0d,read_file_distributed_1d, &
          read_file_distributed_2d,read_file_distributed_3d
  end interface

  interface write_matlab
     module procedure write_matlab_character
     module procedure write_matlab_logical
     module procedure write_matlab_integer
     module procedure write_matlab_real
  end interface

contains


  function new_io_unit() result(next)
    ! NEW_IO_UNIT returns the next available file unit

    implicit none

    ! NEXT = unit number of new file

    integer :: next

    ! MIN = minimum unit number
    ! MAX = maximum unit number
    ! LAST = unit number that was last returned

    integer,parameter :: min=10,max=999
    integer,save :: last=10
    logical :: open,top=.false.

    next = min-1
    if (last>0) then
       next = last+1
       inquire(unit=next,opened=open)
       if (.not.open) last = next
       return
    else
       do
          next = next+1
          inquire(unit=next,opened=open)
          if (.not.open) then
             last = next
             exit
          end if
          if (next==max) then
             if (top) stop 'No available file unit'
             top = .true.
             next = min-1
          end if
       end do
    end if

  end function new_io_unit


  subroutine message(str,routine)
    ! MESSAGE writes informational message

    use mpi_routines, only : myid,stdout

    implicit none

    ! STR = string to be written to standard out
    ! ROUTINE =  subroutine name in which message originated

    character(*),intent(in) :: str
    character(*),intent(in),optional :: routine
   
    character(12) :: id
    character(256) :: str0,str1

    ! write message and subroutine name (if given)

    write(id,'(i10)') myid
    write(id,'(a)') '(' // trim(adjustl(id)) // ') '
    write(str1,'(a)') id // trim(adjustl(str))

    if (present(routine)) then
       write(str0,'(a)') id // 'Message from subroutine: ' &
            // trim(adjustl(routine))
       write(stdout,'(a,/,a)') trim(str0),trim(str1)
    else
       write(stdout,'(a)') trim(str1)
    end if

  end subroutine message


  subroutine messages(str1,str2,str3)
    ! MESSAGES writes informational messages

    use mpi_routines, only : myid,stdout

    implicit none

    ! STR = string to be written to standard out

    character(*),intent(in) :: str1,str2,str3
   
    character(12) :: id
    character(256) :: str(3)

    ! write messages

    write(id,'(i10)') myid
    write(id,'(a)') '(' // trim(adjustl(id)) // ') '
    write(str(1),'(a)') id // trim(adjustl(str1))
    write(str(2),'(a)') id // trim(adjustl(str2))
    write(str(3),'(a)') id // trim(adjustl(str3))

    write(stdout,'(a,/,a,/,a)') trim(str(1)),trim(str(2)),trim(str(3))

  end subroutine messages


  subroutine warning(str,routine)
    ! WARNING writes error message, but does not terminate program

    use mpi_routines, only : myid,stderr

    implicit none

    ! STR = string to be written to standard error
    ! ROUTINE =  subroutine name in which error occurred

    character(*),intent(in) :: str
    character(*),intent(in),optional :: routine

    character(12) :: id
    character(256) :: str0,str1

    ! write error message and subroutine name (if given)
    
    write(id,'(i10)') myid
    write(id,'(a)') '(' // trim(adjustl(id)) // ') '
    write(str1,'(a)') id // trim(adjustl(str))

    if (present(routine)) then
       write(str0,'(a)') id // 'Warning in subroutine: ' &
            // trim(adjustl(routine))
       write(stderr,'(a,/,a)') trim(str0),trim(str1)
    else
       write(stderr,'(a)') trim(str1)
    end if

  end subroutine warning


  subroutine error(str,routine)
    ! ERROR writes error message and terminates program

    use mpi_routines, only : myid,stderr
    use mpi

    implicit none

    ! STR = string to be written to standard error
    ! ROUTINE =  subroutine name in which error occurred

    character(*),intent(in) :: str
    character(*),intent(in),optional :: routine
   
    ! IERR = MPI error flag

    integer :: ierr
    character(12) :: id
    character(256) :: str0,str1,str2

    ! write error message and subroutine name (if given)
    
    write(id,'(i10)') myid
    write(id,'(a)') '(' // trim(adjustl(id)) // ') '
    write(str1,'(a)') id // trim(adjustl(str))
    write(str2,'(a)') id // 'Terminating program'

    if (present(routine)) then
       write(str0,'(a)') id // 'Error in subroutine: ' &
            // trim(adjustl(routine))
       write(stderr,'(a,/,a,/,a)') trim(str0),trim(str1),trim(str2)
    else
       write(stderr,'(a,/,a)') trim(str1),trim(str2)
    end if

    ! terminate program

    call MPI_Abort(MPI_COMM_WORLD,0,ierr)

  end subroutine error


  subroutine seek_to_string(unit,str)

    implicit none

    integer,intent(in) :: unit
    character(*),intent(in) :: str

    character(256) :: temp

    rewind(unit)

    do
       read(unit,*,end=100) temp
       if (temp==str) return
    end do

100 call error("String '" // trim(adjustl(str)) // &
         "' not found in file",'seek_to_string')

  end subroutine seek_to_string


  subroutine copy_text_file(name,input)
    
    implicit none
    
    character(*),intent(in) :: name
    integer,intent(in) :: input

    integer :: echo,stat
    character(256) :: txt

    echo = new_io_unit()
    open(echo,file=trim(name) // '.in',form='formatted',iostat=stat, &
         status='replace')
    if (stat/=0) call error('Error opening ' // trim(name) // '.in','copy_text_file')
    rewind(input)
    do
       read(input,'(a)',end=100) txt
       write(echo,'(a)') trim(txt)
    end do
100 continue
    rewind(input)
    close(echo,iostat=stat)
    if (stat/=0) call error('Error closing ' // trim(name) // '.in','copy_text_file')

  end subroutine copy_text_file


  subroutine open_file_distributed(fid,name,operation,comm,array,precision,offset_in)

    use mpi_routines, only : pw,ps,MPI_REAL_PW,MPI_REAL_PS
    use mpi

    implicit none

    type(file_distributed),intent(out) :: fid
    character(*),intent(in) :: name,operation
    integer,intent(in) :: comm,array,precision
    integer(MPI_OFFSET_KIND),intent(in),optional :: offset_in

    integer :: ierr,try,ntry=5
    integer(MPI_OFFSET_KIND) :: offset
    integer(MPI_OFFSET_KIND),parameter :: zero=0
    character(256) :: str

    ! offset (displacement) at beginning of file

    if (present(offset_in)) then
       offset = offset_in
    else
       offset = zero
    end if

    ! file name

    fid%name = name

    ! precision

    if (precision==pw) then
       fid%pr = pw
       fid%MPI_REAL_PR = MPI_REAL_PW
    else
       fid%pr = ps
       fid%MPI_REAL_PR = MPI_REAL_PS
    end if

    ! open file
    
    do try = 1,ntry,1
      select case(operation)
      case('read')
         call MPI_File_open(comm,fid%name,MPI_MODE_RDONLY, &
              MPI_INFO_NULL,fid%fh,ierr)
      case('write','append')
         call MPI_File_open(comm,fid%name,MPI_MODE_CREATE+MPI_MODE_WRONLY, &
              MPI_INFO_NULL,fid%fh,ierr)
      case default
         call error('Invalid file operation','open_file_distributed') 
      end select
      if (ierr==MPI_SUCCESS) exit
    end do
    if (ierr/=MPI_SUCCESS) then
        call io_error_message('MPI_File_open (on file ' // trim(name) // ')',ierr,str)
        call error(str,'open_file_distributed') 
    end if

    ! if writing, delete contents of file (by setting size to zero)

    if (operation=='write') then
       call MPI_File_set_size(fid%fh,zero,ierr)
       if (ierr/=MPI_SUCCESS) then
          call io_error_message('MPI_File_set_size',ierr,str)
          call error(str,'open_file_distributed') 
       end if
    end if

    ! set view

    call MPI_File_set_view(fid%fh,offset,fid%MPI_REAL_PR,array,'native',MPI_INFO_NULL,ierr)
    if (ierr/=MPI_SUCCESS) then
       call io_error_message('MPI_File_set_view',ierr,str)
       call error(str,'open_file_distributed') 
    end if

  end subroutine open_file_distributed

  
  subroutine write_file_distributed_0d(fid,data)

    use mpi_routines, only : pw,ps
    use mpi

    implicit none
    
    type(file_distributed),intent(in) :: fid
    real,intent(in) :: data

    integer :: ierr
    character(256) :: str

    if (fid%pr==pw) then
       if (collective) then
          call MPI_File_write_all(fid%fh,data,1, &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       else
          call MPI_File_write    (fid%fh,data,1, &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       end if
    else
       if (collective) then
          call MPI_File_write_all(fid%fh,real(data,ps),1, &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       else
          call MPI_File_write    (fid%fh,real(data,ps),1, &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       end if
    end if

    if (ierr/=MPI_SUCCESS) then
       call io_error_message('MPI_File_write_all',ierr,str)
       call error(str,'write_file_distributed_0d')
    end if

  end subroutine write_file_distributed_0d


  subroutine write_file_distributed_1d(fid,data)

    use mpi_routines, only : pw,ps
    use mpi

    implicit none
    
    type(file_distributed),intent(in) :: fid
    real,dimension(:),intent(in) :: data

    integer :: ierr
    character(256) :: str

    if (fid%pr==pw) then
       if (collective) then
          call MPI_File_write_all(fid%fh,data,size(data), &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       else
          call MPI_File_write    (fid%fh,data,size(data), &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       end if
    else
       if (collective) then
          call MPI_File_write_all(fid%fh,real(data,ps),size(data), &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       else
          call MPI_File_write    (fid%fh,real(data,ps),size(data), &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       end if
    end if

    if (ierr/=MPI_SUCCESS) then
       call io_error_message('MPI_File_write_all',ierr,str)
       call error(str,'write_file_distributed_1d')
    end if

  end subroutine write_file_distributed_1d


  subroutine write_file_distributed_2d(fid,data)

    use mpi_routines, only : pw,ps
    use mpi

    implicit none
    
    type(file_distributed),intent(in) :: fid
    real,dimension(:,:),intent(in) :: data

    integer :: ierr
    character(256) :: str

    if (fid%pr==pw) then
       if (collective) then
          call MPI_File_write_all(fid%fh,data,size(data), &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       else
          call MPI_File_write    (fid%fh,data,size(data), &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       end if
    else
       if (collective) then
          call MPI_File_write_all(fid%fh,real(data,ps),size(data), &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       else
          call MPI_File_write    (fid%fh,real(data,ps),size(data), &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       end if
    end if

    if (ierr/=MPI_SUCCESS) then
       call io_error_message('MPI_File_write_all',ierr,str)
       call error(str,'write_file_distributed_2d')
    end if

  end subroutine write_file_distributed_2d


  subroutine write_file_distributed_3d(fid,data)

    use mpi_routines, only : pw,ps
    use mpi

    implicit none
    
    type(file_distributed),intent(in) :: fid
    real,dimension(:,:,:),intent(in) :: data

    integer :: ierr
    character(256) :: str

    if (fid%pr==pw) then
       if (collective) then
          call MPI_File_write_all(fid%fh,data,size(data), &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       else
          call MPI_File_write    (fid%fh,data,size(data), &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       end if
    else
       if (collective) then
          call MPI_File_write_all(fid%fh,real(data,ps),size(data), &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       else
          call MPI_File_write    (fid%fh,real(data,ps),size(data), &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       end if
    end if

    if (ierr/=MPI_SUCCESS) then
       call io_error_message('MPI_File_write_all',ierr,str)
       call error(str,'write_file_distributed_3d')
    end if

  end subroutine write_file_distributed_3d


  subroutine read_file_distributed_0d(fid,data)

    use mpi_routines, only : MPI_REAL_PW
    use mpi

    implicit none
    
    type(file_distributed),intent(in) :: fid
    real,intent(out) :: data

    integer :: ierr
    character(256) :: str

    if (collective) then
       call MPI_File_read_all(fid%fh,data,1,MPI_REAL_PW,MPI_STATUS_IGNORE,ierr)
    else
       call MPI_File_read    (fid%fh,data,1,MPI_REAL_PW,MPI_STATUS_IGNORE,ierr)
    end if
    if (ierr/=MPI_SUCCESS) then
       call io_error_message('MPI_File_read_all',ierr,str)
       call error(str,'read_file_distributed_0d')
    end if

  end subroutine read_file_distributed_0d


  subroutine read_file_distributed_1d(fid,data)

    use mpi_routines, only : MPI_REAL_PW
    use mpi

    implicit none
    
    type(file_distributed),intent(in) :: fid
    real,dimension(:),intent(out) :: data

    integer :: ierr
    character(256) :: str

    if (collective) then
       call MPI_File_read_all(fid%fh,data,size(data),MPI_REAL_PW,MPI_STATUS_IGNORE,ierr)
    else
       call MPI_File_read    (fid%fh,data,size(data),MPI_REAL_PW,MPI_STATUS_IGNORE,ierr)
    end if
    if (ierr/=MPI_SUCCESS) then
       call io_error_message('MPI_File_read_all',ierr,str)
       call error(str,'read_file_distributed_1d')
    end if

  end subroutine read_file_distributed_1d


  subroutine read_file_distributed_2d(fid,data)

    use mpi_routines, only : MPI_REAL_PW
    use mpi

    implicit none
    
    type(file_distributed),intent(in) :: fid
    real,dimension(:,:),intent(out) :: data

    integer :: ierr
    character(256) :: str

    if (collective) then
       call MPI_File_read_all(fid%fh,data,size(data),MPI_REAL_PW,MPI_STATUS_IGNORE,ierr)
    else
       call MPI_File_read    (fid%fh,data,size(data),MPI_REAL_PW,MPI_STATUS_IGNORE,ierr)
    end if
    if (ierr/=MPI_SUCCESS) then
       call io_error_message('MPI_File_read_all',ierr,str)
       call error(str,'read_file_distributed_2d')
    end if

  end subroutine read_file_distributed_2d


  subroutine read_file_distributed_3d(fid,data)

    use mpi_routines, only : MPI_REAL_PW
    use mpi

    implicit none
    
    type(file_distributed),intent(in) :: fid
    real,dimension(:,:,:),intent(out) :: data

    integer :: ierr
    character(256) :: str

    if (collective) then
       call MPI_File_read_all(fid%fh,data,size(data),MPI_REAL_PW,MPI_STATUS_IGNORE,ierr)
    else
       call MPI_File_read    (fid%fh,data,size(data),MPI_REAL_PW,MPI_STATUS_IGNORE,ierr)
    end if
    if (ierr/=MPI_SUCCESS) then
       call io_error_message('MPI_File_read_all',ierr,str)
       call error(str,'read_file_distributed_3d')
    end if

  end subroutine read_file_distributed_3d


  subroutine close_file_distributed(fid)

    use mpi

    implicit none

    type(file_distributed),intent(inout) :: fid

    integer :: ierr
    character(256) :: str

    call MPI_File_close(fid%fh,ierr)
    if (ierr/=MPI_SUCCESS) then
       call io_error_message('MPI_File_close',ierr,str)
       call error(str,'close_file_distributed')
    end if

  end subroutine close_file_distributed


  subroutine write_matlab_character(nunit,name,str,struct)
    ! WRITE_MATLAB_CHARACTER writes character string to matlab file

    implicit none

    ! NUNIT = unit number of output file
    ! NAME = name of variable
    ! STR = value of variable
    ! STRUCT = name of structure in which to write variable

    integer,intent(in) :: nunit
    character(*),intent(in) :: name,str
    character(*),intent(in),optional :: struct

    ! STAT = I/O error flag

    integer :: stat

    ! write character string to file

    if (present(struct)) then
       write(nunit,'(a)',iostat=stat) &
            trim(adjustl(struct)) // '.' // &
            trim(adjustl(name)) // " = '" // trim(adjustl(str)) // "';"
    else
       write(nunit,'(a)',iostat=stat) &
            trim(adjustl(name)) // " = '" // trim(adjustl(str)) // "';"
    end if

    if (stat/=0) call error('Error writing ' // trim(adjustl(name)) // &
         ' to matlab file','write_matlab_character')

  end subroutine write_matlab_character


  subroutine write_matlab_integer(nunit,name,data,struct)
    ! WRITE_MATLAB_INTEGER writes integer to matlab file

    implicit none

    ! NUNIT = unit number of output file
    ! NAME = name of variable
    ! DATA = value of variable
    ! STRUCT = name of structure in which to write variable

    integer,intent(in) :: nunit,data
    character(*),intent(in) :: name
    character(*),intent(in),optional :: struct

    ! STAT = I/O error flag

    integer :: stat
    
    ! write integer to file

    if (present(struct)) then
       write(nunit,'(a,i12,a)',iostat=stat) &
            trim(adjustl(struct)) // '.' // &
            trim(adjustl(name)) // ' = ',data,';'
    else
       write(nunit,'(a,i12,a)',iostat=stat) &
            trim(adjustl(name)) // ' = ',data,';'
    end if
    
    if (stat/=0) call error('Error writing ' // trim(adjustl(name)) // &
         ' to matlab file','write_matlab_integer')

  end subroutine write_matlab_integer


  subroutine write_matlab_real(nunit,name,data,struct)
    ! WRITE_MATLAB_REAL writes real to matlab file

    implicit none

    ! NUNIT = unit number of output file
    ! NAME = name of variable
    ! DATA = value of variable
    ! STRUCT = name of structure in which to write variable

    integer,intent(in) :: nunit
    character(*),intent(in) :: name
    real,intent(in) :: data
    character(*),intent(in),optional :: struct

    ! STAT = I/O error flag

    integer :: stat

    ! write real to file

    if (present(struct)) then
       write(nunit,'(a,e30.15,a)',iostat=stat) &
            trim(adjustl(struct)) // '.' // &
            trim(adjustl(name)) // ' = ',data,';'
    else
       write(nunit,'(a,e30.15,a)',iostat=stat) &
            trim(adjustl(name)) // ' = ',data,';'
    end if

    if (stat/=0) call error('Error writing ' // trim(adjustl(name)) // &
         ' to matlab file','write_matlab_real')

  end subroutine write_matlab_real


  subroutine write_matlab_logical(nunit,name,data,struct)
    ! WRITE_MATLAB_LOGICAL writes logical to matlab file

    implicit none

    ! NUNIT = unit number of output file
    ! NAME = name of variable
    ! DATA = value of variable
    ! STRUCT = name of structure in which to write variable

    integer,intent(in) :: nunit
    character(*),intent(in) :: name
    logical,intent(in) :: data
    character(*),intent(in),optional :: struct

    ! STR = string to hold 'T' or 'F'

    character(1) :: str

    ! assign logical value to string

    if (data) then
       str = 'T'
    else
       str = 'F'
    end if

    ! write string to matlab file

    if (present(struct)) then
       call write_matlab_character(nunit,name,str,struct)
    else
       call write_matlab_character(nunit,name,str)
    end if

  end subroutine write_matlab_logical


  subroutine io_error_message(routine,ierr,str)

    use mpi

    implicit none

    character(*),intent(in) :: routine
    integer,intent(in) :: ierr
    character(*),intent(out) :: str

    !select case(ierr)
    !case(MPI_ERR_IO)
    !   write(str,'(a)') 'Problem with ' // trim(routine) // ': MPI_ERR_IO'
    !case(MPI_ERR_NO_SPACE)
    !   write(str,'(a)') 'Problem with ' // trim(routine) // ': MPI_ERR_NO_SPACE'
    !case(MPI_ERR_NO_SUCH_FILE)
    !   write(str,'(a)') 'Problem with ' // trim(routine) // ': MPI_ERR_NO_SUCH_FILE'
    !case default
       write(str,'(a,i6)') 'Problem with ' // trim(routine) // ': ierr=',ierr
    !end select

  end subroutine io_error_message


  subroutine get_endian(endian)

    implicit none

    character(1),intent(out) :: endian

    integer,parameter :: int=selected_int_kind(9), byte=selected_int_kind(2), &
         n=bit_size(1_int)/bit_size(1_byte)
    integer(int) :: input = 1+256+256**3
    integer(byte) :: bytes(n)

    bytes = transfer (source=input, mold=bytes)

    if (all(bytes==(/ 1,1,0,1 /))) then
       endian = 'l'
    elseif (all(bytes==(/ 1,0,1,1 /))) then
       endian = 'b'
    else
       call error('Endian test failed','get_endian')
    end if

  end subroutine get_endian


end module io
