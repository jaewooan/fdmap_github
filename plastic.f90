module plastic

  implicit none

contains


  subroutine set_rates_plastic(B,G,F,BF,M,mode,dt)

    use grid, only : grid_type,block_grid
    use fields, only : fields_type,block_fields,initial_stress
    use material, only : block_material

    implicit none

    type(block_grid),intent(in) :: B
    type(grid_type),intent(in) :: G
    type(fields_type),intent(inout) :: F
    type(block_fields),intent(in) :: BF
    type(block_material),intent(in) :: M
    integer,intent(in) :: mode
    real,intent(in) :: dt

    integer :: i,j
    real :: s0(6)

    if (M%response=='plastic') then

       do j = B%my,B%py
          do i = B%mx,B%px

             call initial_stress(s0,G%x(i,j),G%y(i,j),BF%F0(4:9),F%problem,mode)
             call plastic_flow(F%DF(i,j,F%nU+1:F%nF),F%Dgammap(i,j),F%lambda(i,j), &
                  M,mode,F%F(i,j,F%nU+1:F%nF),F%gammap(i,j),s0,dt)

          end do
       end do

    end if

  end subroutine set_rates_plastic


  subroutine plastic_flow(Ds,Dgammap,lambda,M,mode,s,gammap,s0,dt)

    use material, only : block_material
    use io, only : error

    implicit none

    ! DS = stress rate in current time step (or RK stage)
    ! LAMBDA = plastic strain rate in current time step (or RK stage)
    ! S0 = initial stress
    ! S = stress change (added to S0 to get absolute stress)
    ! GAMMAP = plastic strain

    real,intent(inout) :: Ds(:),Dgammap,lambda,s(:),gammap
    type(block_material),intent(in) :: M
    integer,intent(in) :: mode
    real,intent(in) :: s0(6),dt

    ! SA = absolute stress tensor
    ! SD = deviatoric stress tensor
    ! W = flow direction in stress space

    real :: tau,sigma,Y
    real,dimension(6) :: sa,sd,W
    real,dimension(6),parameter :: delta = (/ 1d0,0d0,0d0,1d0,0d0,1d0 /) ! Kronecker delta

    ! check if trial update violates yield condition

    call stress(s,s0,sa,tau,sigma,mode)
    call yield(tau,sigma,gammap,Y,M)

    if (Y<=0d0) then ! no plastic flow
       lambda = 0d0
       return
    end if

    ! deviatoric stress (trial)

    sd = sa-sigma*delta

    ! plastic flow:
    ! implicit solution for lambda,s,gammap to
    ! s = s_trial-dt*lambda*W(s)
    ! gammap = gammap_trial+dt*lambda
    ! eta*lambda = Y(s,gammap)
    ! rate-independent limit is eta = 0
    
    ! solution for lambda is closed-form:
    ! regular case for correction onto pressure-dependent yield surface
    ! special case for correction onto tau = 0 yield surface
    ! (selection of special case may not be correct with hardening)

    if (tau*M%mu*M%beta*M%K>(M%mu*sigma-M%b)*M%G) then ! usual case
       !print *, 'regular',Y,tau,sigma
       lambda = (tau+M%mu*sigma-(M%b+M%h*gammap))/(M%eta+dt*(M%h+M%G+M%mu*M%beta*M%K))
    else ! special case, occurs for TPV13
       !print *, 'special',Y,tau,sigma
       lambda = tau/(M%eta+dt*M%G)
    end if
    
    ! update gammap

    gammap = gammap+dt*lambda

    ! update tau and sigma

    tau = tau-dt*lambda*M%G
    sigma = sigma-dt*lambda*M%beta*M%K

    ! check yield using new tau,sigma (UNNECESSARY, DEBUGGING ONLY)
    !call yield(tau,sigma,gammap,Y,M)
    !print *, 'CHECK: Y(tau,sigma,gammap) = ',Y
    
    ! calculate flow direction, using old sd and new tau

    W = M%K*M%beta*delta+M%G*sd/(tau+dt*lambda*M%G)

    ! update deviatoric stress

    sd = sd*tau/(tau+dt*lambda*M%G)

    ! calculate flow direction, using new sd and new tau

    !W = M%K*M%beta*delta+M%G*sd/tau

    ! determine stress components (both methods give identical results)

    sa = sigma*delta+sd
    !sa = sa-dt*lambda*W

    ! extract s from sa = s+s0

    select case(mode)
    case(2)
       s(1) = sa(1)-s0(1) ! sxx
       s(2) = sa(2)-s0(2) ! sxy
       s(3) = sa(4)-s0(4) ! syy
       s(4) = sa(6)-s0(6) ! szz
    case(3)
       s(1) = sa(3)-s0(3) ! sxz
       s(2) = sa(5)-s0(5) ! syz
    end select

    ! check yield using new stress components (UNNECESSARY, DEBUGGING ONLY)
    !call stress(s,s0,sa,tau,sigma,mode)
    !call yield(tau,sigma,gammap,Y,M)
    !print *, 'CHECK: Y(s_ij,gammap) = ',Y
        
    ! adjust rates

    Dgammap = Dgammap+lambda

    select case(mode)
    case(2)
       Ds = Ds-lambda* (/ W(1),W(2),W(4),W(6) /)
    case(3)
       Ds = Ds-lambda* (/ W(3),W(5) /)
    end select
       
  end subroutine plastic_flow


  subroutine invariants(sxx,sxy,sxz,syy,syz,szz,tau,sigma)

    implicit none

    real,intent(in) :: sxx,sxy,sxz,syy,syz,szz
    real,intent(out) :: tau,sigma

    real :: J2

    sigma = (sxx+syy+szz)/3d0
    J2 = ((sxx-syy)**2+(syy-szz)**2+(szz-sxx)**2)/6d0+ &
         sxy**2+sxz**2+syz**2
    tau = sqrt(J2)

  end subroutine invariants


  subroutine stress(s,s0,sa,tau,sigma,mode)

    implicit none

    real,intent(in) :: s(:),s0(6)
    real,intent(out) :: sa(6),tau,sigma
    integer,intent(in) :: mode

    ! s0 = (/ sxx0,sxy0,sxz0,syy0,syz0,szz0 /)
    !           1    2    3    4    5    6

    select case(mode)
    case(2)
       sa(1) = s0(1)+s(1)
       sa(2) = s0(2)+s(2)
       sa(3) = s0(3)
       sa(4) = s0(4)+s(3)
       sa(5) = s0(5)
       sa(6) = s0(6)+s(4)
    case(3)
       sa(1) = s0(1)
       sa(2) = s0(2)
       sa(3) = s0(3)+s(1)
       sa(4) = s0(4)
       sa(5) = s0(5)+s(2)
       sa(6) = s0(6)
    end select

    call invariants(sa(1),sa(2),sa(3),sa(4),sa(5),sa(6),tau,sigma)

  end subroutine stress


  subroutine yield(tau,sigma,gammap,Y,M)

    use material, only : block_material

    implicit none

    real,intent(in) :: tau,sigma,gammap
    real,intent(out) :: Y
    type(block_material),intent(in) :: M

    real :: b

    b = M%b+gammap*M%h ! with hardening
    !if (-M%mu*sigma+b<0d0) print *, 'beyond corner'
    Y = tau-max(-M%mu*sigma+b,0d0)

  end subroutine yield


end module plastic
