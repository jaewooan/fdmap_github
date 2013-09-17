program main

  use checkpoint, only : checkpoint_type,init_checkpoint,checkpoint_read, &
       checkpoint_write_delete,checkpoint_delete,finish_checkpoint,abort_now
  use domain, only : domain_type,init_domain,finish_domain,check_nucleation
  use time_step, only : RK_type,init_RK,time_step_LS,time_step_interseismic
  use io, only : new_io_unit,error,message,write_matlab,get_endian,copy_text_file
  use output, only : output_list,read_output,init_output,destroy_output
  use utilities, only : convert_time,within
  use mpi_routines, only : start_mpi,nprocs,is_master,finish_mpi,clock_split
  use mpi

  implicit none

  character(256) :: name,RKmethod
  character(256) :: str,input_file
  character(1) :: endian
  integer :: n,nt,ninfo,nstart,input,echo,stat,RKorder,hr,mn
  real :: CFL,t,dt,refine,time_per_step, &
       start_time,end_time,total_time,initial_time,sc, &
       tnuc_min,tnuc_max,minV,cyc_dt,cyc_z,cyc_H,cyc_plate_rate, tend
  logical :: cyc,abort,slipping,constant_dt,final_step,interseismic_step,solid
  type(RK_type) :: RK
  type(domain_type) :: D
  type(output_list) :: outlist
  type(checkpoint_type) :: C

  namelist /problem_list/ name,t,nt,CFL,dt,refine,ninfo,Rkmethod,RKorder, &
       tnuc_min,tnuc_max,minV,constant_dt,cyc,cyc_dt,cyc_z,cyc_H,cyc_plate_rate,solid,&
       tend

  ! get problem name

  call get_command_argument(1,input_file,status=stat)
  if (stat/=0) input_file = 'default.in'

  ! start MPI

  call start_mpi
  call clock_split(current_time=start_time)
  initial_time = start_time

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

  cyc = .false.
  cyc_dt = 0d0
  cyc_z = 0d0
  cyc_H = 0d0
  cyc_plate_rate = 0d0

  ! read in problem parameters

  read(input,nml=problem_list,iostat=stat)
  if (stat>0) call error('Error in problem_list','main')

  if (cyc) constant_dt = .false.

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
  D%tD = D%t
  if (cyc.and.C%begin/=0) &
       call message('Must set initial D%t and D%tD properly when starting cycles from checkpoint')
  D%n = C%begin
  call checkpoint_read(name,C,D)

  ! initialization complete

  call clock_split(current_time=end_time)
  total_time = end_time-start_time
  start_time = end_time
  if (is_master) then
     write(str,'(a,f0.6,a)') 'Initialization complete (',total_time,' s)'
     call message(str)
  end if

  ! loop over all time steps

  abort = .false.
  interseismic_step = .false.

  do n = nstart,nt

     ! manually set time (to minimize round-off errors from successively adding dt)
     
     if (constant_dt) D%t = t+dble(n-1)*dt
     
     final_step = (n==nt)
     
     ! advance by one time step

     if (interseismic_step) then
        
        ! increase shear stress and evolve state for cycle simulations
        
        call time_step_interseismic(D,cyc_z,cyc_H,cyc_plate_rate,cyc_dt,outlist)

     else

        ! elastodynamic time step

        select case(RK%method)
        case('LS')
           call time_step_LS(D,dt,CFL,RK,outlist,final_step,solid)
           D%tD = D%tD+dt ! only updated during dynamic steps
        case default
           call error('Invalid RK method','main')
        end select
        
     end if

     D%n = n

     ! determine if abort (immediate checkpoint and termination) is requested

     abort = abort_now(name,C,initial_time)
     if (abort.and.is_master) call message('Aborting now...')

     ! write new checkpoint and delete previous checkpoint after successful write

     call checkpoint_write_delete(name,n,C,D,abort)

     ! status update
     
     if (mod(n-1,ninfo)==0) then
        call clock_split(time_per_step,end_time)
        time_per_step = time_per_step/dble(min(ninfo,n-nstart+1))
        if (is_master) then
           if (cyc) then
              if (interseismic_step) then
                 write(str,'(a,i0,a,i0,a,f0.6,a,f0.6,a,f0.6,a,f0.6,a)') &
                      'time step (QS ) ',n,' of ',nt,': t = ',D%t,', tD = ',D%tD, &
                      ' (wall clock: ',end_time-start_time,' s, ',time_per_step,' s/step)'
              else
                 write(str,'(a,i0,a,i0,a,f0.6,a,f0.6,a,f0.6,a,f0.6,a)') &
                      'time step (DYN) ',n,' of ',nt,': t = ',D%t,', tD = ',D%tD, &
                      ' (wall clock: ',end_time-start_time,' s, ',time_per_step,' s/step)'
              end if
           else
              write(str,'(a,i0,a,i0,a,f0.6,a,f0.6,a,f0.6,a)') &
                   'time step ',n,' of ',nt,': t = ',D%t,' (wall clock: ', &
                   end_time-start_time,' s, ',time_per_step,' s/step)'
           end if
           call message(str)
        end if
     end if
  
     ! check if current slip velocity is too small

     if (cyc) then
        call check_nucleation(D,minV,slipping)
        interseismic_step = .not.slipping
     else
        slipping = .true.
        if (within(tnuc_min,D%t,tnuc_max)) call check_nucleation(D,minV,slipping)
        if (.not.slipping) then
           if (is_master) call message( &
                'Fault slip velocity too small (nucleation failed) -- terminating run...')
           exit
        end if
     end if
     
     ! exit loop if aborting

     if (abort) exit
   
  end do

  ! collect timing information

  call clock_split(current_time=end_time)
  total_time = end_time-start_time

  ! delete final checkpoint (unless aborting)

  if (.not.abort.and.C%del_final) call checkpoint_delete(name,C,D,final=.true.)

  ! display timing information

  if (is_master) then
     call convert_time(total_time,hr,mn,sc)
     write(str,'(a,f0.8,a,i0,a,f0.4,a)') &
          'total: ',total_time,' s for ',n,' time steps (', &
          total_time/dble(n-nstart+1),' s/step)'
     call message(str)
     write(str,'(a,i0,a,i0,a,f0.8,a)') &
          '= ',hr,' hr ',mn,' min ',sc, ' s'
     call message(str)
     write(str,'(i0,a,f0.8,a,f0.4)') &
          nprocs,'  ',total_time,'  ',total_time/dble(n-nstart+1)
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
