module acoustic_gravity

  ! module for additional physics (equation) for acoustic-gravity waves

  implicit none

  ! RHOPRIME = density perturbation
  ! DRHOPRIME = d/dt of rhoprime
  ! G = acceleration due to gravity
  ! DRHO0DY = d/dy of background density
  
  type :: ag_type
     real :: g
     real,dimension(:,:),allocatable :: rhoprime,Drhoprime,drho0dy
  end type ag_type

contains


  subroutine init_acoustic_gravity(AG,C,input,echo)

    use mpi_routines, only : is_master
    use mpi_routines2d, only : cartesian,allocate_array_body
    use io, only : error,write_matlab,seek_to_string  
    
    implicit none

    type(ag_type),intent(out) :: AG
    type(cartesian),intent(in) :: C
    integer,intent(in) :: input,echo

    character(256) :: filename
    integer :: stat
    real :: g,drho0dy

    namelist /acoustic_gravity_list/ filename,g,drho0dy
    
    ! default parameters (spatially uniform)

    filename = ''
    g = 0d0
    drho0dy = 0d0
    
    ! read in acoustic gravity wave parameters

    rewind(input)
    read(input,nml=acoustic_gravity_list,iostat=stat)
    if (stat>0) call error('Error in acoustic_gravity_list','init_acoustic_gravity')

    AG%g = g
    
    ! output acoustic-gravity-wave parameters
    ! (but only if uniform <==> filename=='')
    
    if (is_master) then
       if (filename=='') then
          call write_matlab(echo,'AG.g',g)
          call write_matlab(echo,'AG.drho0dy',drho0dy)
       end if
    end if

    ! allocate memory for fields, set default value to uniform value

    call allocate_array_body(AG%drho0dy  ,C,ghost_nodes=.false.,Fval=drho0dy)
    call allocate_array_body(AG%rhoprime ,C,ghost_nodes=.false.,Fval=0d0)
    call allocate_array_body(AG%Drhoprime,C,ghost_nodes=.false.)

    ! initialize fields from file

    if (filename/='') call acoustic_gravity_IO('read',filename,C,AG)
    
  end subroutine init_acoustic_gravity


  subroutine acoustic_gravity_IO(operation,filename,C,AG)

    use io, only : file_distributed,open_file_distributed, &
         read_file_distributed,write_file_distributed,close_file_distributed
    use mpi_routines2d, only : cartesian
    use mpi_routines, only : pw,is_master
    use mpi

    implicit none

    character(*),intent(in) :: operation,filename
    type(cartesian),intent(in) :: C
    type(ag_type),intent(inout) :: AG

    type(file_distributed) :: fh
    integer :: ierr

    if (operation=='delete') then
       if (is_master) call MPI_file_delete(filename,MPI_INFO_NULL,ierr)
       return
    end if

    call open_file_distributed(fh,filename,operation,C%c2d%comm,C%c2d%array_w,pw)

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    select case(operation)
    case('read')
       call  read_file_distributed(fh,AG%drho0dy(C%mx:C%px,C%my:C%py))
    case('write')
       call write_file_distributed(fh,AG%drho0dy(C%mx:C%px,C%my:C%py))
    end select

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    call close_file_distributed(fh)

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

  end subroutine acoustic_gravity_IO


  subroutine destroy_acoustic_gravity(AG)

    implicit none

    type(ag_type),intent(inout) :: AG

    if (allocated(AG%drho0dy  )) deallocate(AG%drho0dy)
    if (allocated(AG%rhoprime )) deallocate(AG%rhoprime)
    if (allocated(AG%Drhoprime)) deallocate(AG%Drhoprime)

  end subroutine destroy_acoustic_gravity


  subroutine set_rates_acoustic_gravity(C,F,E,AG,mode)
    
    use mpi_routines2d, only : cartesian
    use fields, only : fields_type
    use material, only : elastic_type
    use io, only : error
    
    implicit none

    type(cartesian),intent(in) :: C
    type(fields_type),intent(inout) :: F
    type(ag_type),intent(inout) :: AG
    type(elastic_type),intent(in) :: E
    integer,intent(in) :: mode

    integer :: i,j
    
    select case(mode)
    case(2)
       ! F(:,:,1) = vx
       ! F(:,:,2) = vy
       ! F(:,:,3) = sxx
       ! F(:,:,4) = sxy
       ! F(:,:,5) = syy
       ! DF = d/dt of F
       ! low storage RK method requires DF = DF + new rate, not DF = new rate
       do j = C%my,C%py
          do i = C%mx,C%px
             ! vertical momentum equation (update vy)
             F%DF(i,j,2) = F%DF(i,j,2)-AG%rhoprime(i,j)*AG%g/E%rho(i,j)
             ! normal stress (= -pressure) equation
             !  F(:,:,3) = sxx stress = -pressure
             ! DF(:,:,3) = d/dt of sxx stress = d/dt of -pressure
             !  F(:,:,5) = syy stress = -pressure
             ! DF(:,:,5) = d/dt of syy stress = d/dt of -pressure
             F%DF(i,j,3) = F%DF(i,j,3)-E%rho(i,j)*AG%g*F%F(i,j,2)
             F%DF(i,j,5) = F%DF(i,j,5)-E%rho(i,j)*AG%g*F%F(i,j,2)
             ! density perturbation equation
             AG%Drhoprime(i,j) = AG%Drhoprime(i,j)- &
                  0.5d0/E%cp(i,j)**2*(F%DF(i,j,3)+F%DF(i,j,5))- &
                  (E%rho(i,j)*AG%g/E%cp(i,j)**2+AG%drho0dy(i,j))*F%F(i,j,2)
          end do
       end do
    case(3)
       call error('acoustic gravity not valid in mode 3','set_rates_acoustic_gravity')
    end select

  end subroutine set_rates_acoustic_gravity
  

  subroutine scale_rates_acoustic_gravity(AG,A)

    implicit none

    type(ag_type),intent(inout) :: AG
    real,intent(in) :: A
    
    AG%Drhoprime = A*AG%Drhoprime
    
  end subroutine scale_rates_acoustic_gravity


  subroutine update_fields_acoustic_gravity(AG,dt)

    implicit none

    type(ag_type),intent(inout) :: AG
    real,intent(in) :: dt
    
    AG%rhoprime = AG%rhoprime+dt*AG%Drhoprime
    
  end subroutine update_fields_acoustic_gravity

  
end module acoustic_gravity
