module tsunami

  implicit none

contains


  subroutine initial_sea_surface(x,eta,problem)

    implicit none

    real,intent(in) :: x
    real,intent(out) :: eta
    character(*),intent(in) :: problem

    select case(problem)
    case default
       eta = 0d0
    case('sin')
       eta = sin(x)
    end select

  end subroutine initial_sea_surface


  subroutine seafloor_velocity(x,y,t,vx,vy,problem)

    implicit none

    real,intent(in) :: x,y,t
    real,intent(out) :: vx,vy
    character(*),intent(in) :: problem

    select case(problem)
    case default
       vx = 0d0
       vy = 0d0
    case('slope')
       vx = t*exp(-t)
       vy = 1d0
    end select

  end subroutine seafloor_velocity


end module tsunami
