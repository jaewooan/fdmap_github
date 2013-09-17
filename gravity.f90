module gravity

  implicit none

  type :: gravity_type
     logical :: use_gravity
     real :: gx,gy,gz
  end type gravity_type

contains

  
  subroutine init_gravity(g,input,echo)

    use mpi_routines, only : is_master
    use io, only : error,write_matlab

    implicit none
    
    type(gravity_type),intent(out) :: g
    integer,intent(in) :: input,echo

    real :: gx,gy,gz
    integer :: stat

    namelist /gravity_list/ gx,gy,gz

    ! defaults

    gx = 0d0
    gy = 0d0
    gz = 0d0

    ! read in material parameters
        
    rewind(input)
    read(input,nml=gravity_list,iostat=stat)
    if (stat>0) call error('Error in gravity_list','init_gravity')

    g%use_gravity = (gx/=0d0.or.gy/=0d0.or.gz/=0d0)

    if (.not.g%use_gravity) return

    g%gx = gx
    g%gy = gy
    g%gz = gz

    ! output gravity parameters

    if (is_master) then
       call write_matlab(echo,'gx',g%gx,'g')
       call write_matlab(echo,'gy',g%gy,'g')
       call write_matlab(echo,'gz',g%gz,'g')
    end if

  end subroutine init_gravity


  subroutine set_rates_gravity(F,g,mode)

    use fields, only : fields_type

    implicit none

    type(fields_type),intent(inout) :: F
    type(gravity_type),intent(in) :: g
    integer,intent(in) :: mode

    if (.not.g%use_gravity) return
    
    select case(mode)
    case(2)
       F%DF(:,:,1) = F%DF(:,:,1)+g%gx
       F%DF(:,:,2) = F%DF(:,:,2)+g%gy
    case(3)
       F%DF(:,:,1) = F%DF(:,:,1)+g%gz
    end select

  end subroutine set_rates_gravity


end module gravity
