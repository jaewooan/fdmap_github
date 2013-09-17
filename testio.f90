program testio

  use mpi_routines, only : ps,start_mpi,finish_mpi,is_master,MPI_REAL_PS,master,nprocs,new_communicator,subarray
  use mpi_routines2d, only : cartesian,decompose2d
  use io, only : file_distributed,open_file_distributed, &
       write_file_distributed,close_file_distributed
  use mpi

  implicit none

  ! test IO routines in FDMAP

  integer :: nF,nx,ny,nb,nprocs_x,nprocs_y,n,nt,nfiles,stat, &
       i,comm_fault,array_fault,ierr
  integer,dimension(2) :: blkxm,blkym
  character(256) :: strin,mpi_method,fname
  real :: starttime,totaltime,averagetime,filesize,bandwidth
  real,dimension(:,:,:),allocatable :: F
  real,dimension(:,:),allocatable :: G
  type(cartesian) :: C
  type(file_distributed),pointer :: fh(:),fhf(:)
  logical :: fault
  logical,parameter :: periodic_x=.false.,periodic_y=.false.

  call start_mpi

  if (is_master) then

     ! number of files
     
     call get_command_argument(1,strin,status=stat)
     if (stat/=0) strin = '1'
     read(strin,*) nfiles
     
     ! number of time steps
     
     call get_command_argument(2,strin,status=stat)
     if (stat/=0) strin = '1'
     read(strin,*) nt
     
     ! body, fault, or both
     
     call get_command_argument(3,strin,status=stat)
     if (stat/=0) strin = 'both'

  end if

  call MPI_Bcast(nfiles,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(nt,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(strin,len(strin),MPI_CHARACTER,master,MPI_COMM_WORLD,ierr)

  nF = 1
  nx = 4801
  ny = 2402
  nb = 3
  blkxm = (/ 1,1 /)
  blkym = (/ 1,ny/2+1 /)
  mpi_method = '2D'
  nprocs_x = 1
  nprocs_y = 1

  C%nx = nx
  C%ny = ny

  call decompose2d(C,nF,periodic_x,periodic_y, &
       mpi_method,nprocs_x,nprocs_y,blkxm,blkym)

  if (is_master) print *,'Number of processes = ',nprocs
     
  ! body

  if (strin=='body'.or.strin=='both') then

     if (is_master) print *, 'Testing body file output'

     allocate(F(C%mx-nb:C%px+nb,C%my-nb:C%py+nb,nF))
     F = 0d0

     allocate(fh(nfiles))

     do i = 1,nfiles
        write(fname,'(a,i0,a)') 'data/iotest',i,'.body'
        call open_file_distributed(fh(i),trim(adjustl(fname)),'write',C%c2d%comm,C%c2d%array_s,ps)
     end do

     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     starttime = MPI_WTime()
     do n = 1,nt
        do i = 1,nfiles
           call write_file_distributed(fh(i),F(C%mx:C%px,C%my:C%py,nF))
        end do
     end do
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     totaltime = MPI_WTime()-starttime

     if (is_master) then
        
        averagetime = totaltime/dble(nt)
        filesize = dble(nx*ny*4)/dble(1048576)
        bandwidth = dble(nfiles)*filesize/averagetime

        print *, 'Number of files = ',nfiles
        print *, 'File size = ',filesize,' MB'
        print *, 'Total output time for ',nt,' time steps (all files) = ',totaltime,' s'
        print *, 'Average output time per time step (all files) = ',averagetime,' s'
        print *, 'Bandwidth = ',bandwidth,'MB/s'

     end if

     do i = 1,nfiles
        call close_file_distributed(fh(i))
     end do

     deallocate(F)
     deallocate(fh)
     nullify(fh)
     
  end if

  ! fault

  if (strin=='fault'.or.strin=='both') then

     if (is_master) then
        print *
        print *, 'Testing fault file output'
     end if

     fault = (C%my<=ny/2).and.(ny/2<=C%py)
          
     call new_communicator(fault,comm_fault)
     
     call subarray(nx,C%mx,C%px,MPI_REAL_PS,array_fault)
     
     if (fault) then

        allocate(G(C%mx:C%px,nF))
        G = 0d0

        allocate(fhf(nfiles))

        do i = 1,nfiles
           write(fname,'(a,i0,a)') 'data/iotest',i,'.fault'
           call open_file_distributed(fhf(i),trim(adjustl(fname)),'write',comm_fault,array_fault,ps)
        end do
        
     end if

     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     starttime = MPI_WTime()
     do n = 1,nt
        do i = 1,nfiles
           if (fault) call write_file_distributed(fhf(i),G(C%mx:C%px,nF))
        end do
     end do
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     totaltime = MPI_WTime()-starttime

     if (is_master) then

        averagetime = totaltime/dble(nt)
        filesize = dble(nx*4)/dble(1048576)
        bandwidth = dble(nfiles)*filesize/averagetime
        
        print *, 'Number of files = ',nfiles
        print *, 'File size = ',filesize,' MB'
        print *, 'Total output time for ',nt,' time steps (all files) = ',totaltime,' s'
        print *, 'Average output time per time step (all files) = ',averagetime,' s'
        print *, 'Bandwidth = ',bandwidth,'MB/s'
        
     end if
     
     do i = 1,nfiles
        if (fault) call close_file_distributed(fhf(i))
     end do

     if (fault) then
        deallocate(G)
        deallocate(fhf)
        nullify(fhf)
     end if
     
  end if

  call finish_mpi

end program testio
