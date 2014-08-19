module plastic

  implicit none

contains


  subroutine update_fields_plastic(B,G,F,BF,M,E,mode,dt)

    use grid, only : grid_type,block_grid
    use fields, only : fields_type,block_fields,initial_stress
    use material, only : block_material,elastic_type

    implicit none

    type(block_grid),intent(in) :: B
    type(grid_type),intent(in) :: G
    type(fields_type),intent(inout) :: F
    type(block_fields),intent(in) :: BF
    type(block_material),intent(in) :: M
    type(elastic_type),intent(in) :: E
    integer,intent(in) :: mode
    real,intent(in) :: dt

    integer :: i,j
    real :: Gel,Kel,s0(6),ep(6),Wp

    logical :: heterogeneous

    if (M%response/='plastic') return

    heterogeneous = allocated(E%rho)
    
    do j = B%my,B%py
       do i = B%mx,B%px
          
          if (heterogeneous) then
             Gel = E%G(i,j)
             Kel = E%M(i,j)-1.333333333333333d0*E%G(i,j)
          else
             Gel = M%G
             Kel = M%K
          end if

          call initial_stress(s0,G%x(i,j),G%y(i,j),BF%F0(4:9),F%problem,mode)
          call plastic_flow(F%DF(i,j,F%nU+1:F%nF),F%lambda(i,j), &
               M,Gel,Kel,mode,F%F(i,j,F%nU+1:F%nF),F%gammap(i,j),s0,dt,ep,Wp)

          if (allocated(F%Ep)) F%ep(i,j,:) = F%ep(i,j,:)+ep
          if (allocated(F%Wp)) F%Wp(i,j,:) = F%Wp(i,j,:)+Wp
          
       end do
    end do

  end subroutine update_fields_plastic


  subroutine plastic_flow(Ds,lambda,M,G,K,mode,s,gammap,s0,dt,ep,Wp)

    use material, only : block_material
    use io, only : error

    implicit none

    ! DS = stress rate in current time step (or RK stage)
    ! LAMBDA = plastic strain rate in current time step (or RK stage)
    ! S0 = initial stress
    ! S = stress change (added to S0 to get absolute stress)
    ! GAMMAP = plastic strain
    ! EP = plastic strain tensor increment
    ! WP = plastic work increment

    real,intent(inout) :: Ds(:),lambda,s(:),gammap
    real,intent(out) :: ep(:),Wp
    type(block_material),intent(in) :: M
    integer,intent(in) :: mode
    real,intent(in) :: G,K,s0(6),dt

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

    if (tau*M%mu*M%beta*K>(M%mu*sigma-M%b)*G) then ! usual case
       !print *, 'regular',Y,tau,sigma
       lambda = (tau+M%mu*sigma-(M%b+M%h*gammap))/(M%eta+dt*(M%h+G+M%mu*M%beta*K))
    else ! special case, occurs for TPV13
       !print *, 'special',Y,tau,sigma
       lambda = tau/(M%eta+dt*G)
    end if
    
    ! update gammap

    gammap = gammap+dt*lambda

    ! update tau and sigma

    tau = tau-dt*lambda*G
    sigma = sigma-dt*lambda*M%beta*K

    ! check yield using new tau,sigma (UNNECESSARY, DEBUGGING ONLY)
    !call yield(tau,sigma,gammap,Y,M)
    !print *, 'CHECK: Y(tau,sigma,gammap) = ',Y
    
    ! calculate flow direction, using old sd and new tau

    W = K*M%beta*delta+G*sd/(tau+dt*lambda*G)

    ! update deviatoric stress

    sd = sd*tau/(tau+dt*lambda*G)

    ! calculate flow direction, using new sd and new tau

    !W = K*M%beta*delta+G*sd/tau

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

    gammap = gammap+dt*lambda

    select case(mode)
    case(2)
       Ds = Ds-lambda* (/ W(1),W(2),W(4),W(6) /)
    case(3)
       Ds = Ds-lambda* (/ W(3),W(5) /)
    end select
    
    ! plastic strain tensor increment
    ep = dt*lambda*(sd/(2d0*tau)+(M%beta/3d0)*delta)

    ! plastic work rate increment
    Wp = dt*lambda*(tau+(M%beta/3d0)*sigma)
       
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
