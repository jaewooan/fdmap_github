module output

  use io, only : file_distributed

  implicit none

  ! OUTPUT_ITEM is a derived type containing output variables
  !
  ! NAME = name of item
  ! FIELD = field components to output
  ! FH = file handle
  ! FHT = file handle for time
  ! FHTD = file handle for dynamic time
  ! ST = stride in time
  ! TMIN = minimum value of time output window
  ! TMAX = maximum value of t
  ! NMIN = minimum time step in time output window
  ! NMAX = maximum time step in time output window
  ! MX = minimum index in x direction (point/grid output)
  ! PX = maximum index in x direction (point/grid output)
  ! SX = stride in x direction (point/grid output)
  ! MY = minimum index in y direction (point/grid output)
  ! PY = maximum index in y direction (point/grid output)
  ! SY = stride in y direction (point/grid output)
  ! X = value of x (point output)
  ! Y = value of y (point output)
  ! LOCATION = location of field (e.g., body, block, iface, bndX)
  ! IFACE = iface index (only used if location=iface)
  ! BLOCK = block index (only used if location=block)
  ! OUTPUT_PROCESS = flag indicating if process holds data to be output
  ! OUTPUT_PROCESS_T = flag indicating if process holds time data to be output
  ! ARRAY = MPI distributed array type
  ! COMM = MPI communicator

  type :: output_item
     character(256) :: name,location
     character(6) :: field
     logical :: output_process,output_process_t
     integer :: nmin,nmax,st,mx,px,sx,my,py,sy,nx,ny,comm,array,iface,block
     real :: x,y,tmin,tmax
     type(file_distributed) :: fh,fht
  end type output_item

  ! OUTPUT_NODE is a derived type for basic node in singly-linked list 
  ! (single output item variables and pointer to next node)
  !
  ! ITEM = output item
  ! NEXT = pointer to next node

  type :: output_node
     private
     type(output_item),pointer :: item=>null()
     type(output_node),pointer :: next=>null()
  end type output_node

  ! OUTPUT_LIST is a derived type containing all output variables
  !
  ! ROOT = start of singly-linked list of output items

  type :: output_list
     type(output_node),pointer :: root=>null()
  end type output_list

contains


  subroutine read_output(list,D,input)
    ! READ_OUTPUT reads in output variables from file

    use domain, only : domain_type
    use utilities, only : word_count,extract_words,within
    use io, only : error

    implicit none

    ! LIST = output list
    ! D = domain variables
    ! INPUT = file unit for input file

    type(output_list),intent(out) :: list
    type(domain_type),intent(in) :: D
    integer,intent(in) :: input

    ! NAME = name of output item
    ! TEMP = temporary character string
    ! IFLD = index of field
    ! NFLD = number of fields
    ! FLDS = field names
    ! ITEM = output item

    integer :: ifld,nfld,mx,sx,px,my,sy,py,i,j
    character(256) :: name,temp,loci,locj
    character(6),dimension(:),allocatable :: flds
    type(output_item) :: item

    ! move to start of output list

    rewind(input)
    do
       read(input,'(a)',end=100) temp
       if (temp=='!---begin:output_list---') exit
    end do
       
    ! read output data and add nodes to list

    do

       ! read item data

       ! name

       read(input,'(a)',end=100) name

       ! return if at end of list

       if (name=='!---end:output_list---') return

       ! location

       item%iface = -1 ! default
       item%block = -1 ! default

       read(input,'(a)') item%location
       if (index(item%location,'iface')/=0) then
          if (item%location(1:5)=='point') then
             read(item%location(12:),'(i6)') item%iface
             item%location = 'point_iface' // D%I(item%iface)%direction
          elseif (item%location(1:4)=='grid') then
             read(item%location(11:),'(i6)') item%iface
             item%location = 'grid_iface' // D%I(item%iface)%direction
          else
             read(item%location(6:),'(i6)') item%iface
             item%location = 'iface' // D%I(item%iface)%direction
          end if
          if (.not.within(1,item%iface,D%nifaces)) &
               call error('Invalid iface number','read_output')
       end if
       if (item%location(1:5)=='block') then
          read(item%location(6:),'(i6)') item%block
          if (.not.within(1,item%block,D%nblocks)) &
               call error('Invalid block number','read_output')
          item%location = 'block'
       end if
       if (item%location(1:3)=='bnd') then
          read(item%location(5:),'(i6)') item%block
          if (.not.within(1,item%block,D%nblocks)) &
               call error('Invalid block number','read_output')
          item%location = item%location(1:4)
       end if
       if (item%location(1:6)=='Eblock') then
          read(item%location(7:),'(i6)') item%block
          if (.not.within(1,item%block,D%nblocks)) &
               call error('block number exceeds nblocks','read_output')
          item%location = 'Eblock'
       end if

       ! fields

       read(input,'(a)') temp
       nfld = word_count(temp)
       allocate(flds(nfld))
       call extract_words(temp,flds)

       ! time

       read(input,*) item%tmin,item%tmax,item%st

       ! indices (dependent on location)

       item%sx = 1 ! default
       item%sy = 1 ! default

       select case(item%location)
       case default
          call error('Invalid output location: ' // trim(adjustl(item%location)) , &
               'read_output')
       case('body','Ebody','block','Eblock','ifacex','ifacey','bndL','bndR','bndB','bndT')
          ! nothing to be done
       case('slicex')
          read(input,*) item%my
          if (.not.(within(1,item%my,D%C%ny))) &
               call error('Indices out of range for output item ' // &
               trim(adjustl(name)) ,'read_output')
       case('slicey')
          read(input,*) item%mx
          if (.not.(within(1,item%mx,D%C%nx))) &
               call error('Indices out of range for output item ' // &
               trim(adjustl(name)) ,'read_output')
       case('point_ifacey')
          read(input,*) item%mx
          if (.not.(within(D%I(item%iface)%mg,item%mx,D%I(item%iface)%pg))) &
               call error('Indices out of range for output item ' // &
               trim(adjustl(name)) ,'read_output')
       case('point_ifacex')
          read(input,*) item%my
          if (.not.(within(D%I(item%iface)%mg,item%my,D%I(item%iface)%pg))) &
               call error('Indices out of range for output item ' // &
               trim(adjustl(name)) ,'read_output')
       case('grid_ifacey')
          read(input,*) mx,sx,px
          if (.not.( &
               within(D%I(item%iface)%mg,mx,D%I(item%iface)%pg).and. &
               within(D%I(item%iface)%mg,px,D%I(item%iface)%pg).and.(mx<=px))) &
               call error('Indices out of range for output item ' // &
               trim(adjustl(name)) ,'read_output')
       case('grid_ifacex')
          read(input,*) my,sy,py
          if (.not.( &
               within(D%I(item%iface)%mg,my,D%I(item%iface)%pg).and. &
               within(D%I(item%iface)%mg,py,D%I(item%iface)%pg).and.(my<=py))) &
               call error('Indices out of range for output item ' // &
               trim(adjustl(name)) ,'read_output')
       case('point_body')
          read(input,*) item%mx,item%my
          if (.not.( &
               within(1,item%mx,D%C%nx).and. &
               within(1,item%my,D%C%ny))) &
               call error('Indices out of range for output item ' // &
               trim(adjustl(name)) ,'read_output')
       case('grid_body')
          read(input,*) mx,sx,px,my,sy,py
          if (.not.( &
               within(1,mx,D%C%nx).and.within(1,px,D%C%nx).and. &
               within(1,my,D%C%ny).and.within(1,py,D%C%ny).and. &
               (mx<=px).and.(my<=py))) &
               call error('Indices out of range for output item ' // &
               trim(adjustl(name)) ,'read_output')
       case('body_grid')
          read(input,*) item%mx,item%sx,item%px,item%my,item%sy,item%py
          if (.not.( &
               within(1,item%mx,D%C%nx).and.within(1,item%px,D%C%nx).and. &
               within(1,item%my,D%C%ny).and.within(1,item%py,D%C%ny).and. &
               (item%mx<=item%px).and.(item%my<=item%py))) &
               call error('Indices out of range for output item ' // &
               trim(adjustl(name)) ,'read_output')
       end select

       ! store as nodes in list

       select case(item%location)
       case default
          do ifld = 1,nfld
             item%name = trim(name) // '_' // flds(ifld)
             item%field = flds(ifld)
             call create_output_node(list%root,item)
          end do
       case('grid_ifacey')
          item%location = 'point_ifacey'
          do i = mx,px,sx
             item%mx = i
             do ifld = 1,nfld
                write(loci,'(i64)') item%mx
                item%name = trim(name) // '_' // trim(adjustl(loci)) &
                     // '_' // flds(ifld)
                item%field = flds(ifld)
                call create_output_node(list%root,item)
             end do
          end do
       case('grid_ifacex')
          item%location = 'point_ifacex'
          do j = my,py,sy
             item%my = j
             do ifld = 1,nfld
                write(loci,'(i64)') item%my
                item%name = trim(name) // '_' // trim(adjustl(loci)) &
                     // '_' // flds(ifld)
                item%field = flds(ifld)
                call create_output_node(list%root,item)
             end do
          end do
       case('grid_body')
          item%location = 'point_body'
          do j = my,py,sy
             item%my = j
             do i = mx,px,sx
                item%mx = i
                do ifld = 1,nfld
                   write(loci,'(i64)') item%mx
                   write(locj,'(i64)') item%my
                   item%name = trim(name) // '_' // trim(adjustl(loci)) &
                        // '_' // trim(adjustl(locj)) // '_' // flds(ifld)
                   item%field = flds(ifld)
                   call create_output_node(list%root,item)
                end do
             end do
          end do
       end select
       deallocate(flds)

    end do
    
100 call error("Error reading list 'output' in input file",'read_output')

  end subroutine read_output


  recursive subroutine create_output_node(node,item)
    ! CREATE_OUTPUT_NODE recursively traverses list and adds node at end

    implicit none

    ! I/O Parameters
    ! NODE = output node
    ! ITEM = output item

    type(output_node),pointer :: node
    type(output_item),intent(in) :: item
    
    if (.not.associated(node)) then
       allocate(node)
       nullify(node%next)
       allocate(node%item)
       node%item = item
    else
        call create_output_node(node%next,item)
    endif

  end subroutine create_output_node


  subroutine init_output(echo,pb_name,list,D,checkpoint_number,dt)
    ! INIT_OUTPUT initializes output variables and creates output items

    use domain, only : domain_type

    implicit none
    
    ! I/O Parameters:
    ! ECHO = unit number for output file
    ! PB_NAME = name of problem
    ! LIST = output variables

    integer,intent(in) :: echo
    character(*),intent(in) :: pb_name
    type(output_list),intent(inout) :: list
    type(domain_type),intent(in) :: D
    integer,intent(in) :: checkpoint_number
    real,intent(in) :: dt

    ! traverse list and initialize all nodes

    call init_output_node(list%root,echo,D,pb_name,checkpoint_number,dt)

  end subroutine init_output

  
  recursive subroutine init_output_node(node,echo,D,pb_name,checkpoint_number,dt)
    ! INIT_OUTPUT_NODE recursively initializes all output nodes

    use domain, only : domain_type,check_field

    implicit none
    
    ! I/O Parameters:
    ! NODE = output node
    ! ECHO = unit number for output file
    ! PB_NAME = name of problem

    type(output_node),pointer :: node
    integer,intent(in) :: echo
    type(domain_type),intent(in) :: D
    character(*),intent(in) :: pb_name
    integer,intent(in) :: checkpoint_number
    real,intent(in) :: dt

    ! exit recursion if at end of list

    if (.not.associated(node)) return

    ! allocate and initialize output item

    call init_output_item(echo,node%item,D,dt)

    ! check if field can be output (i.e., check allocation status, etc.)

    if (node%item%output_process) &
         call check_field(D,node%item%field,node%item%location,node%item%iface)

    ! open file for output
    
    call open_output_file(node%item,pb_name,checkpoint_number)

    ! repeat for next node

    call init_output_node(node%next,echo,D,pb_name,checkpoint_number,dt)    

  end subroutine init_output_node


  subroutine init_output_item(echo,item,D,dt)
    ! INIT_OUTPUT_ITEM creates and initializes a output item

    use domain, only : domain_type
    use io, only : write_matlab,error
    use utilities, only : within
    use mpi_routines, only : is_master,master,MPI_REAL_PW,MPI_REAL_PS,subarray,new_communicator
    use mpi

    implicit none
    
    ! ECHO = unit number for output file
    ! ITEM = output item

    integer,intent(in) :: echo
    type(output_item),intent(inout) :: item
    type(domain_type),intent(in) :: D
    real,intent(in) :: dt

    integer :: ierr,mx,px,my,py,mgx,mgy,nx,ny,myid
    integer,parameter :: tagx=11,tagy=12

    ! output time step window

    item%nmin = ceiling((item%tmin-D%t)/dt)
    if ((item%tmax-D%t)/dt>dble(huge(1))) then
       item%nmax = huge(1)
    else
       item%nmax = floor((item%tmax-D%t)/dt)
    end if

    ! spatial indices, etc.

    item%nx = 1
    item%ny = 1

    select case(item%location)

    case default

       call error('Invalid item location','init_output_item')

    case('point_body')

       item%px = item%mx
       item%nx = 1

       item%py = item%my
       item%ny = 1

       item%output_process = &
            within(D%C%mx,item%mx,D%C%px).and. &
            within(D%C%my,item%my,D%C%py)
  
       if (item%output_process) then

          item%comm = MPI_COMM_SELF
          item%array = MPI_REAL_PS

          item%x = D%G%x(item%mx,item%my)
          item%y = D%G%y(item%mx,item%my)
          call MPI_Send(item%x,1,MPI_REAL_PW,master,tagx,MPI_COMM_WORLD,ierr)
          call MPI_Send(item%y,1,MPI_REAL_PW,master,tagy,MPI_COMM_WORLD,ierr)

       end if

    case('body')

       item%mx = D%C%mx
       item%px = D%C%px
       item%nx = D%C%nx

       item%my = D%C%my
       item%py = D%C%py
       item%ny = D%C%ny

       item%output_process = .true.
       call MPI_Comm_dup(D%C%c2d%comm,item%comm,ierr)

       select case(item%field)
       case default
          item%array = D%C%c2d%array_s
       case('F')
          item%array = D%C%c2d%arrayF_s
       end select

    case('body_grid')

       ! for now, item%m and item%p are global limits of section of original array

       ! item%m is exact starting index within original array;
       ! ensure that item%p is proper ending index within original array
       item%px = floor(real(item%px-item%mx)/real(item%sx))*item%sx+item%mx
       item%py = floor(real(item%py-item%my)/real(item%sy))*item%sy+item%my

       item%output_process = &
            (item%mx<=D%C%px.and.D%C%mx<=item%px).and. &
            (item%my<=D%C%py.and.D%C%my<=item%py)

       call new_communicator(item%output_process,item%comm,D%C%c2d%comm)

       if (item%output_process) then

          ! size of strided array
          item%nx = (item%px-item%mx)/item%sx+1
          item%ny = (item%py-item%my)/item%sy+1

          !print *, '(',myid,')  ',item%nx,item%mx,item%px

          ! global limits of original array
          mgx = item%mx
          mgy = item%my

          ! processor-specific limits of original array
          item%mx = max(item%mx,ceiling(real(D%C%mx-mgx)/real(item%sx))*item%sx+mgx)
          item%px = min(item%px,floor  (real(D%C%px-mgx)/real(item%sx))*item%sx+mgx)
          item%my = max(item%my,ceiling(real(D%C%my-mgy)/real(item%sy))*item%sy+mgy)
          item%py = min(item%py,floor  (real(D%C%py-mgy)/real(item%sy))*item%sy+mgy)

          ! processor-specific limits of strided array (section of 1:item%n)
          mx = (item%mx-mgx)/item%sx+1
          px = (item%px-mgx)/item%sx+1
          my = (item%my-mgy)/item%sy+1
          py = (item%py-mgy)/item%sy+1
          
          !print *, '(',myid,')  ',mx,px,item%mx,item%px
          
          call subarray(item%nx,item%ny,mx,px,my,py,MPI_REAL_PS,item%array)

       end if

    case('slicex')

       item%output_process = within(D%C%my,item%my,D%C%py)

       if (item%output_process) then

          item%mx = D%C%mx
          item%px = D%C%px
          item%nx = D%C%nx
          
          item%py = item%my
          item%ny = 1
          
          call MPI_Comm_dup(D%C%c1dx%comm,item%comm,ierr)
          item%array = D%C%c1dx%array_s

       end if

    case('slicey')

       item%output_process = within(D%C%mx,item%mx,D%C%px)

       if (item%output_process) then

          item%px = item%mx
          item%nx = 1
          
          item%my = D%C%my
          item%py = D%C%py
          item%ny = D%C%ny
          
          call MPI_Comm_dup(D%C%c1dy%comm,item%comm,ierr)
          item%array = D%C%c1dy%array_s
          
       end if

    case('point_ifacex')

       item%output_process = D%I(item%iface)%output_process &
            .and.within(D%I(item%iface)%m,item%my,D%I(item%iface)%p)

       if (item%output_process) then

          item%comm = MPI_COMM_SELF
          
          item%py = item%my
          item%ny = 1
          
          select case(item%field)
          case default
             item%mx = 1
             item%px = 1
             item%nx = 1
             item%array = MPI_REAL_PS
          case('T','P')
             item%mx = 1
             item%px = D%I(item%iface)%TP%n
             item%nx = D%I(item%iface)%TP%n
             call subarray(item%nx,item%mx,item%px,MPI_REAL_PS,item%array)
          case('z')
             print *, 'add ability to output z'
             stop
          end select
          
          item%x = D%B(D%I(item%iface)%iblockm)%G%bndR%x(item%my)
          item%y = D%B(D%I(item%iface)%iblockm)%G%bndR%y(item%my)
          call MPI_Send(item%x,1,MPI_REAL_PW,master,tagx,MPI_COMM_WORLD,ierr)
          call MPI_Send(item%y,1,MPI_REAL_PW,master,tagy,MPI_COMM_WORLD,ierr)

       end if

    case('point_ifacey')

       item%output_process = D%I(item%iface)%output_process &
            .and.within(D%I(item%iface)%m,item%mx,D%I(item%iface)%p)

       if (item%output_process) then

          item%comm = MPI_COMM_SELF
          
          item%px = item%mx
          item%nx = 1
          
          select case(item%field)
          case default
             item%my = 1
             item%py = 1
             item%ny = 1
             item%array = MPI_REAL_PS
          case('T','P')
             item%my = 1
             item%py = D%I(item%iface)%TP%n
             item%ny = D%I(item%iface)%TP%n
             call subarray(item%ny,item%my,item%py,MPI_REAL_PS,item%array)
          case('z')
             print *, 'add ability to output z'
             stop
          end select
          
          item%x = D%B(D%I(item%iface)%iblockm)%G%bndT%x(item%mx)
          item%y = D%B(D%I(item%iface)%iblockm)%G%bndT%y(item%mx)
          call MPI_Send(item%x,1,MPI_REAL_PW,master,tagx,MPI_COMM_WORLD,ierr)
          call MPI_Send(item%y,1,MPI_REAL_PW,master,tagy,MPI_COMM_WORLD,ierr)

       end if

    case('ifacex')

       item%output_process = D%I(item%iface)%output_process

       if (item%output_process) then
          
          call MPI_Comm_dup(D%I(item%iface)%comm,item%comm,ierr)
          
          item%my = D%I(item%iface)%m
          item%py = D%I(item%iface)%p
          item%ny = D%I(item%iface)%pg-D%I(item%iface)%mg+1
          
          select case(item%field)
          case default
             item%mx = 1
             item%px = 1
             item%nx = 1
             item%array = D%I(item%iface)%array_s
          case('T','P')
             item%mx = 1
             item%px = D%I(item%iface)%TP%n
             item%nx = D%I(item%iface)%TP%n
             call subarray(item%nx,item%ny, &
                  item%mx,item%px, &
                  item%my-D%I(item%iface)%mg+1,item%py-D%I(item%iface)%mg+1, &
                  MPI_REAL_PS,item%array)
          case('z')
             print *, 'add ability to output z'
             stop
          end select

       end if

    case('ifacey')

       item%output_process = D%I(item%iface)%output_process

       if (item%output_process) then

          call MPI_Comm_dup(D%I(item%iface)%comm,item%comm,ierr)
          
          item%mx = D%I(item%iface)%m
          item%px = D%I(item%iface)%p
          item%nx = D%I(item%iface)%pg-D%I(item%iface)%mg+1
          
          select case(item%field)
          case default
             item%my = 1
             item%py = 1
             item%ny = 1
             item%array = D%I(item%iface)%array_s
          case('T','P')
             item%my = 1
             item%py = D%I(item%iface)%TP%n
             item%ny = D%I(item%iface)%TP%n
             call subarray(item%nx,item%ny, &
                  item%mx-D%I(item%iface)%mg+1,item%px-D%I(item%iface)%mg+1, &
                  item%my,item%py, &
                  MPI_REAL_PS,item%array)
          case('z')
             print *, 'add ability to output z'
             stop
          end select

       end if

    case('block')

       item%output_process = .not.D%B(item%block)%G%skip

       if (item%output_process) then
          
          item%mx = D%B(item%block)%G%mx
          item%px = D%B(item%block)%G%px
          item%nx = D%B(item%block)%G%nx

          item%my = D%B(item%block)%G%my
          item%py = D%B(item%block)%G%py
          item%ny = D%B(item%block)%G%ny

          call MPI_Comm_dup(D%B(item%block)%G%comm,item%comm,ierr)
          call subarray(item%nx,item%ny, &
               item%mx-D%B(item%block)%G%mgx+1,item%px-D%B(item%block)%G%mgx+1, &
               item%my-D%B(item%block)%G%mgy+1,item%py-D%B(item%block)%G%mgy+1, &
               MPI_REAL_PS,item%array)

       end if

    case('Eblock')

       if (D%B(item%block)%G%skip) then
          item%output_process = .false.
       else
          item%output_process = (D%B(item%block)%G%myid==0)
       end if

       if (item%output_process) then
          
          item%mx = 1
          item%px = 1
          item%nx = 1

          item%my = 1
          item%py = 1
          item%ny = 1

          item%comm = MPI_COMM_SELF
          call subarray(item%nx,item%mx,item%px,MPI_REAL_PS,item%array)

       end if

    case('Ebody')

       item%output_process = is_master

       if (item%output_process) then
          
          item%mx = 1
          item%px = 1
          item%nx = 1

          item%my = 1
          item%py = 1
          item%ny = 1

          item%comm = MPI_COMM_SELF
          call subarray(item%nx,item%mx,item%px,MPI_REAL_PS,item%array)

       end if

    case('bndL')

       item%output_process = D%B(item%block)%G%sideL

       if (item%output_process) then

          item%mx = D%B(item%block)%G%mx
          item%px = D%B(item%block)%G%mx
          item%nx = 1

          item%my = D%B(item%block)%G%my
          item%py = D%B(item%block)%G%py
          item%ny = D%B(item%block)%G%ny

          call MPI_Comm_dup(D%B(item%block)%F%bndFL%comm,item%comm,ierr)
          call subarray(item%ny, &
               item%my-D%B(item%block)%G%mgy+1,item%py-D%B(item%block)%G%mgy+1, &
               MPI_REAL_PS,item%array)

       end if

    case('bndR')

       item%output_process = D%B(item%block)%G%sideR
       
       if (item%output_process) then
          
          item%mx = D%B(item%block)%G%px
          item%px = D%B(item%block)%G%px
          item%nx = 1

          item%my = D%B(item%block)%G%my
          item%py = D%B(item%block)%G%py
          item%ny = D%B(item%block)%G%ny

          call MPI_Comm_dup(D%B(item%block)%F%bndFR%comm,item%comm,ierr)
          call subarray(item%ny, &
               item%my-D%B(item%block)%G%mgy+1,item%py-D%B(item%block)%G%mgy+1, &
               MPI_REAL_PS,item%array)

       end if

    case('bndB')

       item%output_process = D%B(item%block)%G%sideB

       if (item%output_process) then

          item%mx = D%B(item%block)%G%mx
          item%px = D%B(item%block)%G%px
          item%nx = D%B(item%block)%G%nx

          item%my = D%B(item%block)%G%my
          item%py = D%B(item%block)%G%my
          item%ny = 1

          call MPI_Comm_dup(D%B(item%block)%F%bndFB%comm,item%comm,ierr)
          call subarray(item%nx, &
               item%mx-D%B(item%block)%G%mgx+1,item%px-D%B(item%block)%G%mgx+1, &
               MPI_REAL_PS,item%array)

       end if

    case('bndT')

       item%output_process = D%B(item%block)%G%sideT

       if (item%output_process) then

          item%mx = D%B(item%block)%G%mx
          item%px = D%B(item%block)%G%px
          item%nx = D%B(item%block)%G%nx

          item%my = D%B(item%block)%G%py
          item%py = D%B(item%block)%G%py
          item%ny = 1

          call MPI_Comm_dup(D%B(item%block)%F%bndFT%comm,item%comm,ierr)
          call subarray(item%nx, &
               item%mx-D%B(item%block)%G%mgx+1,item%px-D%B(item%block)%G%mgx+1, &
               MPI_REAL_PS,item%array)

       end if

    end select

    ! time output handled by master process in item%comm

    if (item%output_process) then
       call MPI_Comm_rank(item%comm,myid,ierr)
       item%output_process_t = (myid==0)
    else
       item%output_process_t = .false.
    end if

    ! make sure master has proper sizes
    call MPI_Reduce(item%nx,nx,1,MPI_INTEGER,MPI_MAX,master,MPI_COMM_WORLD,ierr)
    call MPI_Reduce(item%ny,ny,1,MPI_INTEGER,MPI_MAX,master,MPI_COMM_WORLD,ierr)

    if (is_master) then

       call write_matlab(echo,'name',item%name,item%name)
       call write_matlab(echo,'field',item%field,item%name)
       call write_matlab(echo,'location',item%location,item%name)
       
       call write_matlab(echo,'tmin',item%tmin,item%name)
       call write_matlab(echo,'tmax',item%tmax,item%name)
       call write_matlab(echo,'stride_t',item%st,item%name)

       call write_matlab(echo,'nx',nx,item%name)
       call write_matlab(echo,'ny',ny,item%name)

       if (index(item%location,'point')/=0) then
          call MPI_Recv(item%x,1,MPI_REAL_PW,MPI_ANY_SOURCE,tagx, &
               MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
          call MPI_Recv(item%y,1,MPI_REAL_PW,MPI_ANY_SOURCE,tagy, &
               MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
          call write_matlab(echo,'x',item%x,item%name)
          call write_matlab(echo,'y',item%y,item%name)
       end if

       if (index(item%location,'block')/=0.or.index(item%location,'bnd')/=0) &
            call write_matlab(echo,'block',item%block,item%name)
       if (index(item%location,'iface')/=0) &
            call write_matlab(echo,'iface',item%iface,item%name)

    end if

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

  end subroutine init_output_item


  subroutine open_output_file(item,pb_name,checkpoint_number)
    ! OPEN_OUTPUT_FILE opens a output file

    use io, only : open_file_distributed
    use mpi_routines, only : ps,MPI_REAL_PS
    use mpi

    implicit none

    ! I/O Parameters:
    ! ITEM = output item
    ! PB_NAME = problem name
    
    type(output_item),intent(inout) :: item
    character(*),intent(in) :: pb_name
    integer,intent(in) :: checkpoint_number

    ! Internal Parameters:
    ! NAME = file name
    ! OPERATION = action to perform on file
    ! IERR = MPI error flag

    character(256) :: name,operation
    integer :: noffset,ierr
    integer(MPI_OFFSET_KIND) :: offset

    if (checkpoint_number>0) then
       operation = 'append'
    else
       operation = 'write'
    end if

    ! number of previously output time steps to skip, excluding present step if output before

    if (checkpoint_number<item%nmin) then ! nothing output
       noffset = 0
    else
       if (mod(checkpoint_number-item%nmin,item%st)==0) then
          ! current step was output (and will now be overwritten)
          noffset = (checkpoint_number-item%nmin)/item%st
       else ! current step was not output
          noffset = (checkpoint_number-item%nmin)/item%st+1
       end if
    end if

    ! special case of checkpoint_number>item%nmax, no further output

    if (noffset>item%nmax) noffset = 0

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    ! open data file with appropriate offset if restarting from checkpoint
       
    if (item%output_process) then
       name = trim(adjustl(pb_name)) // '_' // trim(adjustl(item%name)) // '.dat'
       offset = int(noffset,kind(offset))*int(item%nx,kind(offset))*int(item%ny,kind(offset))*int(ps,kind(offset))
       call open_file_distributed(item%fh,name,operation,item%comm,item%array,ps,offset)
    end if

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    ! open output time file with appropriate offset if restarting from checkpoint

    if (item%output_process_t) then
       name = trim(adjustl(pb_name)) // '_' // trim(adjustl(item%name)) // '.t'
       offset = noffset*ps
       call open_file_distributed(item%fht ,name,operation,MPI_COMM_SELF,MPI_REAL_PS,ps,offset)
    end if

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

  end subroutine open_output_file


  subroutine write_output(list,D)
    ! WRITE_OUTPUT obtains array of field values, appropriately strided, and calls
    ! a routine to write them to a data file

    use domain, only : domain_type

    implicit none

    ! LIST = output variables
    ! T = current time

    type(output_list),intent(inout) :: list
    type(domain_type),intent(in) :: D

    ! recursively traverse list to output all output items

    call write_output_node(list%root,D)

  end subroutine write_output


  recursive subroutine write_output_node(node,D)
    ! WRITE_OUTPUT_NODE recursively outputs data for nodes in list

    use domain, only : domain_type,write_field
    use io, only : write_file_distributed
    use utilities, only : within

    implicit none

    ! NODE = output node
    ! D = domain variables

    type(output_node),pointer :: node
    type(domain_type),intent(in) :: D

    ! return if node does not exist (end of list)

    if (.not.associated(node)) return

    ! proceed only if current time is within output time window
          
    if (within(node%item%nmin,D%n,node%item%nmax)) then

       ! proceed only if stride of potential time output index is acceptable
       
       if (mod(D%n-node%item%nmin,node%item%st)==0) then
          
          ! write output time
          
          if (node%item%output_process_t) &
               call write_file_distributed(node%item%fht ,(/D%t /))

          ! write data
          
          if (node%item%output_process) then
             if     (node%item%iface/=-1) then
                call write_field(node%item%fh,D, &
                     node%item%field,node%item%location, &
                     node%item%mx,node%item%px,node%item%sx, &
                     node%item%my,node%item%py,node%item%sy, &
                     node%item%iface)
             elseif (node%item%block/=-1) then
                call write_field(node%item%fh,D, &
                     node%item%field,node%item%location, &
                     node%item%mx,node%item%px,node%item%sx, &
                     node%item%my,node%item%py,node%item%sy, &
                     node%item%block)
             else
                call write_field(node%item%fh,D, &
                     node%item%field,node%item%location, &
                     node%item%mx,node%item%px,node%item%sx, &
                     node%item%my,node%item%py,node%item%sy)
             end if
          end if

       end if

    end if

    ! move to next output node
    
    call write_output_node(node%next,D)
       
  end subroutine write_output_node


  subroutine destroy_output(list)
    ! DESTROY_OUTPUT destroys derived type list

    implicit none

    ! LIST = output variables

    type(output_list),intent(inout) :: list

    ! recursively traverse list, destroying each node (slow implementation)

    do while(associated(list%root)) 
       call destroy_output_node(list%root)
    end do

  end subroutine destroy_output


  recursive subroutine destroy_output_node(node)
    ! DESTROY_OUTPUT_NODE recursively destroys all nodes in list

    use io, only : close_file_distributed

    implicit none

    ! NODE = output node

    type(output_node),pointer :: node

    ! move to end of list

    if (associated(node%next)) then

       ! not at end, move to next node

       call destroy_output_node(node%next)

    else

       ! end of list, destroy node

       if (associated(node%item)) then
          
          ! close time file
          
          if (node%item%output_process_t) &
               call close_file_distributed(node%item%fht )

          ! close data file
          
          if (node%item%output_process) call close_file_distributed(node%item%fh)
          
          ! deallocate and nullify output item
          
          deallocate(node%item)
          node%item => null()
          
       end if
       
       ! deallocate and nullify output node

       deallocate(node)
       node => null()
       
       return

    end if

  end subroutine destroy_output_node


end module output

