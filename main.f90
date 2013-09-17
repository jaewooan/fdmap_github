program main

  use checkpoint, only : checkpoint_type,init_checkpoint,checkpoint_read, &
       checkpoint_write_delete,checkpoint_delete,finish_checkpoint,abort_now
  use domain, only : domain_type,init_domain,finish_domain,check_nucleation
  use time_step, only : RK_type,init_RK,time_step_LS
  use io, only : new_io_unit,error,message,write_matlab,get_endian,copy_text_file
  use output, only : output_list,read_output,init_output,destroy_output
  use utilities, only : convert_time,within
  use mpi_routines, only : start_mpi,nprocs,is_master,finish_mpi,time_elapsed
  use mpi

  implicit none

  character(256) :: name,RKmethod
  character(256) :: str,input_file
  character(1) :: endian
  integer :: n,nt,ninfo,nstart,input,echo,stat,RKorder,hr,mn
  real :: CFL,t,dt,refine,tend,sc, &
       tnuc_min,tnuc_max,minV, &
       walltime_start,walltime_total,walltime_per_step_ave, &
       walltime1,walltime2,walltime_per_step
  logical :: abort,slipping,constant_dt,final_step,solid
  type(RK_type) :: RK
  type(domain_type) :: D
  type(output_list) :: outlist
  type(checkpoint_type) :: C

  namelist /problem_list/ name,t,nt,CFL,dt,refine,ninfo,Rkmethod,RKorder, &
       tnuc_min,tnuc_max,minV,constant_dt,solid,tend

  ! get problem name

  call get_command_argument(1,input_file,status=stat)
  if (stat/=0) input_file = 'default.in'

  ! start MPI

  call start_mpi

  ! print number of processors
  
  if (is_master) then
     call message('Reading input file ' // trim(adjustl(input_file)))
     write(str,'(i12)') nprocs
     call message('Number of processors = ' // trim(str))
  end if

  ! open input file

  input = new_io_unit()
  open(input,file=input_file,iostat=stat,status='old')
  if (stat/=0) call error('Error opening ' // trim(input_file),'main')

  ! defaults

  name = 'default'
  t = 0d0
  nt = 0
  CFL = 0d0
  dt = 0d0
  refine = 1d0
  tend = 0d0
  ninfo = 1
  constant_dt = .true.

  Rkmethod ='LS'
  RKorder = 3

  solid = .true.

  tnuc_min = huge(tnuc_min)
  tnuc_max = huge(tnuc_max)
  minV = 0d0

  ! read in problem parameters

  read(input,nml=problem_list,iostat=stat)
  if (stat>0) call error('Error in problem_list','main')

  if (is_master) &
     call message('Starting problem: ' // trim(adjustl(name)))

  ! initialize checkpointing

  call init_checkpoint(name,C,input,nstart)

  ! copy input file

  if (is_master) call copy_text_file(name,input)

  ! open parameter file

  if (is_master) then
     echo = new_io_unit()
     open(echo,file=trim(name) // '.m',form='formatted',iostat=stat, &
          status='replace')
     if (stat/=0) call error('Error opening ' // trim(name) // '.m','main')
  end if

  ! initialize Runge-Kutta time stepping

  RK%method = RKmethod
  RK%order = RKorder
  call init_RK(RK)

  ! initialize domain and set time step corresponding to CFL condition

  call init_domain(D,t,refine,CFL,dt,input,echo,RK%c)

  D%n = 0
  if(nt .eq. 0) then
    nt = ceiling(tend/dt)
  else
    nt = ceiling(dble(nt-1)*refine)+1
  end if

  ! write parameters

  if (is_master)  then
     call get_endian(endian)
     call write_matlab(echo,'endian',endian)
     call write_matlab(echo,'nt',nt)
     call write_matlab(echo,'refine',refine)
     call write_matlab(echo,'solid',solid)
     call write_matlab(echo,'method',RK%method,'RK')
     call write_matlab(echo,'order',RK%order,'RK')
  end if

  ! read output list

  call read_output(outlist,D,input)

  ! intialize output

  call init_output(echo,name,outlist,D,C%begin,dt)

  ! close input file

  close(input,iostat=stat)
  if (stat/=0) call error('Error closing ' // trim(input_file),'main')

  ! close parameter file

  if (is_master) then
     close(echo,iostat=stat)
     if (stat/=0) call error('Error closing ' // trim(name) // '.m','main')
  end if

  ! read checkpoint

  D%t = t+dble(C%begin)*dt
  D%n = C%begin
  call checkpoint_read(name,C,D)

  ! initialization complete

  if (is_master) then
     write(str,'(a,f0.6,a)') 'Initialization complete (',time_elapsed(),' s)'
     call message(str)
  end if

  ! loop over all time steps

  abort = .false.
  walltime_start = time_elapsed() ! start timing for average time/step

  do n = nstart,nt

     ! time each step

     walltime1 = time_elapsed()

     ! manually set time (to minimize round-off errors from successively adding dt)
     
     if (constant_dt) D%t = t+dble(n-1)*dt
     
     final_step = (n==nt)
     
     ! advance by one time step

     select case(RK%method)
     case('LS')
        call time_step_LS(D,dt,RK,outlist,final_step,solid)
     case default
        call error('Invalid RK method','main')
     end select
     
     D%n = n

     ! determine if abort (immediate checkpoint and termination) is requested

     abort = abort_now(name,C)
     if (abort.and.is_master) call message('Aborting now...')

     ! write new checkpoint and delete previous checkpoint after successful write

     call checkpoint_write_delete(name,n,C,D,abort)

     ! status update
     
     if (mod(n-1,ninfo)==0) then
        walltime2 = time_elapsed()
        walltime_per_step = walltime2-walltime1
        if (is_master) then
           write(str,'(a,i0,a,i0,a,f0.6,a,f0.6,a,f0.6,a)') &
                'time step ',n,' of ',nt,': t = ',D%t,' (wall clock: ', &
                walltime2,' s, ',walltime_per_step,' s/step)'
           call message(str)
        end if
     end if
  
     ! check if current slip velocity is too small

     slipping = .true.
     if (within(tnuc_min,D%t,tnuc_max)) call check_nucleation(D,minV,slipping)
     if (.not.slipping) then
        if (is_master) call message( &
             'Fault slip velocity too small (nucleation failed) -- terminating run...')
        exit
     end if
     
     ! exit loop if aborting

     if (abort) exit

     ! calculate average wall clock time per time step

     if (n==nt-1) walltime_per_step_ave = (time_elapsed()-walltime_start)/dble(nt-nstart)

  end do

  ! delete final checkpoint (unless aborting)

  if (.not.abort.and.C%del_final) call checkpoint_delete(name,C,D,final=.true.)

  ! collect timing information

  walltime_total = time_elapsed()

  ! display timing information

  if (is_master) then
     call convert_time(walltime_total,hr,mn,sc)
     write(str,'(a,f0.8,a,i0,a,f0.4,a)') &
          'total: ',walltime_total,' s for ',nt,' time steps (',walltime_per_step_ave,' s/step)'
     call message(str)
     write(str,'(a,i0,a,i0,a,f0.8,a)') &
          '= ',hr,' hr ',mn,' min ',sc, ' s'
     call message(str)
     write(str,'(a,i0,a,f0.8)') &
          'timing: ',nprocs,'  ',walltime_per_step_ave
     call message(str)
     if (n-1==nt) then
        call message('Finished problem: ' // trim(adjustl(name)))
     else
        call message('Early termination of problem: ' // trim(adjustl(name)))
     end if
  end if

  ! finish
  
  call finish_checkpoint(C)
  call finish_domain(D)
  call destroy_output(outlist)

  ! finish MPI

  call finish_mpi

end program main
