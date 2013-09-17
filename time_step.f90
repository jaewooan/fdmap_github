module time_step

  implicit none

  type :: RK_type
     character(3) :: method
     integer :: order,nstage
     real,dimension(:),allocatable :: A,B,C
  end type RK_type

contains

  
  subroutine init_RK(RK)

    use io, only : error

    implicit none

    type(RK_type),intent(inout) :: RK

    ! assumes RK%method and RK%order have been set
    ! (or will use default values)

    select case(RK%method)

    case default

       call error('Invalid Runge-Kutta method','init_RK')

    case('LS')

       select case(RK%order)
       case default
          call error('Invalid Runge-Kutta order','init_RK')
       case(1) ! explicit Euler
          RK%nstage = 1
          allocate(RK%A(RK%nstage),RK%B(RK%nstage),RK%C(0:RK%nstage))
          RK%A = (/ 0d0 /)
          RK%B = (/ 1d0 /)
          RK%C = (/ 0d0,1d0 /)
       case(2) ! Heun
          RK%nstage = 2
          allocate(RK%A(RK%nstage),RK%B(RK%nstage),RK%C(0:RK%nstage))
          RK%A = (/ 0d0,-1d0 /)
          RK%B = (/ 1d0,1d0/2d0 /)
          RK%C = (/ 0d0,1d0,1d0 /)
       case(3) ! Williamson (3,3)
          RK%nstage = 3
          allocate(RK%A(RK%nstage),RK%B(RK%nstage),RK%C(0:RK%nstage))
          RK%A = (/ 0d0, -5d0/9d0, -153d0/128d0 /)
          RK%B = (/ 1d0/3d0, 15d0/16d0, 8d0/15d0 /)
          RK%C = (/ 0d0, 1d0/3d0, 3d0/4d0, 1d0 /)
       case(4) ! Kennedy-Carpenter (5,4)
          RK%nstage = 5
          allocate(RK%A(RK%nstage),RK%B(RK%nstage),RK%C(0:RK%nstage))
          RK%A = (/ 0d0, -567301805773d0/1357537059087d0, -2404267990393d0/2016746695238d0 &
               , -3550918686646d0/2091501179385d0, -1275806237668d0/842570457699d0 /)
          RK%B = (/ 1432997174477d0/9575080441755d0, 5161836677717d0/13612068292357d0 &
               , 1720146321549d0/2090206949498d0, 3134564353537d0/4481467310338d0 &
               , 2277821191437d0/14882151754819d0 /)
          RK%C = (/ 0d0, 1432997174477d0/9575080441755d0, 2526269341429d0/6820363962896d0 &
               , 2006345519317d0/3224310063776d0, 2802321613138d0/2924317926251d0, 1d0 /)
       end select

    end select

  end subroutine init_RK


  subroutine time_step_LS(D,dt,CFL,RK,outlist,final_step,solid)
    ! low storage Runge-Kutta methods

    use domain, only : domain_type,exchange_fields_edges,interior_to_edges,interior_rates_to_edges,exchange_rates_field_edge
    use fields, only : scale_rates_interior,update_fields_interior, &
         exchange_fields,scale_rates_pml,update_pml,scale_rates_displacement,update_displacement
    use boundaries, only : apply_bc,scale_rates_iface,update_fields_iface, &
         update_fields_iface_implicit
    use fault, only : couple_blocks
    use energy, only : energy_interior,energy_block, &
         scale_rates_energy,set_rates_energy,update_energy
    use output, only : output_list,write_output
    use plastic, only : set_rates_plastic
    use gravity, only : set_rates_gravity
    use fd_coeff, only : H00i
    use io, only : error

    implicit none

    type(domain_type),intent(inout) :: D
    real,intent(in) :: dt,CFL
    type(RK_type),intent(in) :: RK
    type(output_list),intent(inout) :: outlist
    logical,intent(in) :: final_step,solid

    integer :: i,im,ip,stage
    real :: t0,Ks,Kp

    ! initial time
    
    t0 = D%t

    ! loop over RK stages

    do stage = 1,RK%nstage

       ! multiply interior rates by RK coefficient A

       if (solid) call scale_rates_interior(D%C,D%F,RK%A(stage))

       ! multiply interface rates by RK coefficient A

       do i = 1,D%nifaces
          call scale_rates_iface(D%I(i),RK%A(stage))
       end do

       ! multiply energy dissipation and pml rates by RK coefficient A

       do i = 1,D%nblocks
          call scale_rates_energy(D%B(i)%G,D%B(i)%F,RK%A(stage))
          call scale_rates_displacement(D%B(i)%G,D%B(i)%F,RK%A(stage))
          call scale_rates_pml(D%B(i)%F,RK%A(stage),D%B(i)%G,D%F%nf)
       end do

       ! adjust fields on boundaries to satisfy bc
       
       do i = 1,D%nblocks
         !JK: done at end of step now
         ! call interior_to_edges(D%B(i)%G,D%F,D%B(i)%F)
          call apply_bc(D%B(i)%G,D%B(i)%F,D%B(i)%B,D%B(i)%M,D%mode,D%t,i)
       end do
       
       ! adjust fields on interfaces to satisfy jump conditions
       ! (and set state rate)
       
       do i = 1,D%nifaces
          im = D%I(i)%iblockm
          ip = D%I(i)%iblockp
          !JK: done at end of step now call
          ! exchange_fields_edges(D%I(i),D%C,D%B(im)%F,D%B(ip)%F)
          call couple_blocks(D%I(i),D%B(im)%F,D%B(ip)%F,D%B(im)%M,D%B(ip)%M, &
               D%C,D%mode,D%t,initialize=.false.)
       end do
       
       ! enforce boundary conditions with SAT method: add additional source term to rates

       do i = 1,D%nblocks
          Ks = D%B(i)%M%cs*H00i ! SAT penalty parameter, S-wave
          Kp = D%B(i)%M%cp*H00i ! SAT penalty parameter, P-wave
          call set_rates_SAT(D%B(i)%G,D%G,D%F,D%B(i)%F,D%mode,Ks,Kp)
       end do
       
       ! set rates

       if (solid) then ! can turn off this time-consuming step

          do i = 1,D%nblocks
             if (D%B(i)%dissipation) then
                call set_rates_dissipation(D%B(i)%G,D%G,D%F,D%B(i)%M,D%C,D%mode, &
                     D%B(i)%c_dissipation*CFL/dt)
                !call set_rates_dissipation_xonly(D%B(i)%G,D%G,D%F,D%B(i)%M,D%C,D%mode, &
                !     D%B(i)%c_dissipation*CFL/dt) ! enable manually if desired
             else
                call set_rates(D%B(i)%G,D%G,D%F,D%B(i)%F,D%B(i)%M,D%C,D%mode)
             end if
          end do
       
       end if

       call set_rates_gravity(D%F,D%grav,D%mode)

       ! energy dissipation rate

       do i = 1,D%nblocks
          call set_rates_energy(D%B(i)%G,D%B(i)%F)
       end do

       ! compute energy

       if (stage==1) then
          call energy_interior(D)
       end if

       ! output fields
       
       if (stage==1) call write_output(outlist,D)

       ! return if final time step

       if (final_step) return

       ! set the source terms if this is an mms problem

       do i = 1,D%nblocks
         call set_source(D%B(i)%G,D%G,D%F,D%B(i)%M,D%C,D%t,D%mode,i)
       end do

       ! update interior fields
       
       D%t = t0+RK%C(stage)*dt
       if (solid) call update_fields_interior(D%C,D%F,RK%B(stage)*dt)
       
       ! adjust rates and fields for plastic response (implicit-explicit approach)
       
       !do i = 1,D%nblocks
       !   call set_rates_plastic(D%B(i)%G,D%F,D%B(i)%F,D%B(i)%M,D%mode,RK%B(stage)*dt)
       !end do
       
       ! update iface fields, except thermal pressurization
       

       do i = 1,D%nblocks
         call interior_to_edges(D%B(i)%G,D%F,D%B(i)%F)
         call interior_rates_to_edges(D%B(i)%G,D%F,D%B(i)%F,RK%A(stage))
       end do
       do i = 1,D%nifaces
         im = D%I(i)%iblockm
         ip = D%I(i)%iblockp
         call exchange_fields_edges(D%I(i),D%C,D%B(im)%F,D%B(ip)%F)
         call exchange_rates_field_edge(D%I(i),D%C,D%B(im)%F,D%B(ip)%F)
       end do
       do i = 1,D%nifaces
         im = D%I(i)%iblockm
         ip = D%I(i)%iblockp
         ! we need both the current and future time
         call update_fields_iface(D%I(i),RK%B(stage)*dt,t0+RK%C(stage-1)*dt,D%t,&
           D%B(im)%F,D%B(ip)%F,D%B(im)%M,D%B(ip)%M,RK%A(stage))
       end do

       ! update iface fields (implicit only, thermal pressurization)
       
       if (.not.D%operator_split) then
          do i = 1,D%nifaces
             call update_fields_iface_implicit(D%I(i),stage)
          end do
       end if

       ! update dissipated energy and pml
       
       do i = 1,D%nblocks
          call update_displacement(D%B(i)%G,D%B(i)%F,RK%B(stage)*dt,D%F%nU)
          call update_energy(D%B(i)%G,D%B(i)%F,RK%B(stage)*dt)
          call update_pml(D%B(i)%G,D%B(i)%F,RK%B(stage)*dt,D%F%nf)
          ! if(i .eq. 1) then
          !   print*,'W',D%B(i)%F%Wy(21,21,:)
          !   print*,'F',D%F%F(21,21,:)
          ! end if
       end do
       
       if (stage/=RK%nstage) call exchange_fields(D%C,D%F)

    end do

    ! operator splitting: implicit Euler integration of plasticity and thermal pressurization

    ! update iface fields (thermal pressurization)
    ! note this requires that f,V,N be solved for at end of time step
    ! (otherwise values would be at final RK stage)

    if (D%operator_split) then
       
       ! adjust fields on interfaces to satisfy jump conditions
       
       do i = 1,D%nifaces
          if (D%I(i)%TP%use_TP) then
             im = D%I(i)%iblockm
             ip = D%I(i)%iblockp
             call exchange_fields_edges(D%I(i),D%C,D%B(im)%F,D%B(ip)%F)
             call couple_blocks(D%I(i),D%B(im)%F,D%B(ip)%F,D%B(im)%M,D%B(ip)%M, &
                  D%C,D%mode,D%t,initialize=.false.)
          end if
       end do
       
       do i = 1,D%nifaces
          if (D%I(i)%TP%use_TP) call update_fields_iface_implicit(D%I(i),1)
       end do

    end if

    ! adjust rates and fields for plastic response
    
    do i = 1,D%nblocks
       call set_rates_plastic(D%B(i)%G,D%G,D%F,D%B(i)%F,D%B(i)%M,D%mode,dt)
    end do
              
    call exchange_fields(D%C,D%F)

  end subroutine time_step_LS


  subroutine set_rates(B,G,F,BF,M,C,mode)

    use grid, only : grid_type,block_grid
    use fields, only : fields_type, block_fields
    use rates, only : &
         rates_interior_mode2,rates_boundary_mode2, &
         rates_interior_mode3,rates_boundary_mode3, &
         rates_szz,&
         rates_interior_mode2_pml,rates_boundary_mode2_pml
    use material, only : block_material
    use mpi_routines2d, only : cartesian

    implicit none

    type(block_grid),intent(in) :: B
    type(grid_type),intent(in) :: G
    type(fields_type),intent(inout) :: F
    type(block_fields),intent(inout) :: BF
    type(block_material),intent(in) :: M
    type(cartesian),intent(in) :: C
    integer,intent(in) :: mode

    if (B%skip) return ! process has no cells in this block

    select case(mode)

    case(2)

      if(M%pmlx .or. M%pmly) then
        ! calculate rates at interior points

        call rates_interior_mode2_pml(C%mbx,C%mby,C%mx,C%my,B%mx,B%my,&
          F%F,F%DF,BF%Wx,BF%DWx,BF%Wy,BF%DWy,&
          M%rho,M%G,M%M,M%gamma,&
          M%pmltype,M%pmlax,M%pmlsx,M%pmllx,M%pmlay,M%pmlsy,M%pmlly,&
          G%x,G%y,&
          G%xq,G%xr,G%yq,G%yr,G%Ji,B%mix,B%pix,B%miy,B%piy,M%pmlx,M%pmly)

        ! calculate rates at boundary points

        if (B%sideL) call rates_boundary_mode2_pml(C%mbx,C%mby,C%mx,C%my,B%mx,B%my,&
          F%F,F%DF,BF%Wx,BF%DWx,BF%Wy,BF%DWy,&
          M%rho,M%G,M%M,M%gamma,&
          M%pmltype,M%pmlax,M%pmlsx,M%pmllx,M%pmlay,M%pmlsy,M%pmlly,&
          G%x,G%y,&
          G%xq,G%xr,G%yq,G%yr,G%Ji,B%mgx  ,B%mix-1,B%miy  ,B%piy,M%pmlx,M%pmly,&
          'L')
        if (B%sideR) call rates_boundary_mode2_pml(C%mbx,C%mby,C%mx,C%my,B%mx,B%my,&
          F%F,F%DF,BF%Wx,BF%DWx,BF%Wy,BF%DWy,&
          M%rho,M%G,M%M,M%gamma,&
          M%pmltype,M%pmlax,M%pmlsx,M%pmllx,M%pmlay,M%pmlsy,M%pmlly,&
          G%x,G%y,&
          G%xq,G%xr,G%yq,G%yr,G%Ji,B%pix+1,B%pgx  ,B%miy  ,B%piy,M%pmlx,M%pmly,&
          'R')
        if (B%sideB) call rates_boundary_mode2_pml(C%mbx,C%mby,C%mx,C%my,B%mx,B%my,&
          F%F,F%DF,BF%Wx,BF%DWx,BF%Wy,BF%DWy,&
          M%rho,M%G,M%M,M%gamma,&
          M%pmltype,M%pmlax,M%pmlsx,M%pmllx,M%pmlay,M%pmlsy,M%pmlly,&
          G%x,G%y,&
          G%xq,G%xr,G%yq,G%yr,G%Ji,B%mix  ,B%pix  ,B%mgy,B%miy-1,M%pmlx,M%pmly,&
          'B')
        if (B%sideT) call rates_boundary_mode2_pml(C%mbx,C%mby,C%mx,C%my,B%mx,B%my,&
          F%F,F%DF,BF%Wx,BF%DWx,BF%Wy,BF%DWy,&
          M%rho,M%G,M%M,M%gamma,&
          M%pmltype,M%pmlax,M%pmlsx,M%pmllx,M%pmlay,M%pmlsy,M%pmlly,&
          G%x,G%y,&
          G%xq,G%xr,G%yq,G%yr,G%Ji,B%mix  ,B%pix  ,B%piy+1,B%pgy,M%pmlx,M%pmly,&
          'T')

        if (B%sideL.and.B%sideB) call rates_boundary_mode2_pml(C%mbx,C%mby,C%mx,C%my,B%mx,B%my,&
          F%F,F%DF,BF%Wx,BF%DWx,BF%Wy,BF%DWy,&
          M%rho,M%G,M%M,M%gamma,&
          M%pmltype,M%pmlax,M%pmlsx,M%pmllx,M%pmlay,M%pmlsy,M%pmlly,&
          G%x,G%y,&
          G%xq,G%xr,G%yq,G%yr,G%Ji,B%mgx  ,B%mix-1,B%mgy,B%miy-1,M%pmlx,M%pmly,&
          'LB')
        if (B%sideR.and.B%sideB) call rates_boundary_mode2_pml(C%mbx,C%mby,C%mx,C%my,B%mx,B%my,&
          F%F,F%DF,BF%Wx,BF%DWx,BF%Wy,BF%DWy,&
          M%rho,M%G,M%M,M%gamma,&
          M%pmltype,M%pmlax,M%pmlsx,M%pmllx,M%pmlay,M%pmlsy,M%pmlly,&
          G%x,G%y,&
          G%xq,G%xr,G%yq,G%yr,G%Ji,B%pix+1,B%pgx  ,B%mgy,B%miy-1,M%pmlx,M%pmly,&
          'RB')
        if (B%sideL.and.B%sideT) call rates_boundary_mode2_pml(C%mbx,C%mby,C%mx,C%my,B%mx,B%my,&
          F%F,F%DF,BF%Wx,BF%DWx,BF%Wy,BF%DWy,&
          M%rho,M%G,M%M,M%gamma,&
          M%pmltype,M%pmlax,M%pmlsx,M%pmllx,M%pmlay,M%pmlsy,M%pmlly,&
          G%x,G%y,&
          G%xq,G%xr,G%yq,G%yr,G%Ji,B%mgx  ,B%mix-1,B%piy+1,B%pgy,M%pmlx,M%pmly,&
          'LT')
        if (B%sideR.and.B%sideT) call rates_boundary_mode2_pml(C%mbx,C%mby,C%mx,C%my,B%mx,B%my,&
          F%F,F%DF,BF%Wx,BF%DWx,BF%Wy,BF%DWy,&
          M%rho,M%G,M%M,M%gamma,&
          M%pmltype,M%pmlax,M%pmlsx,M%pmllx,M%pmlay,M%pmlsy,M%pmlly,&
          G%x,G%y,&
          G%xq,G%xr,G%yq,G%yr,G%Ji,B%pix+1,B%pgx  ,B%piy+1,B%pgy,M%pmlx,M%pmly,&
          'RT')
      else
        ! calculate rates at interior points

        call rates_interior_mode2(C%mbx,C%mby,C%mx,C%my, &
          F%F,F%DF,M%rho,M%G,M%M,M%gamma,G%xq,G%xr,G%yq,G%yr,G%Ji,B%mix,B%pix,B%miy,B%piy)

        ! calculate rates at boundary points

        if (B%sideL) call rates_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
          F%F,F%DF,M%rho,M%G,M%M,M%gamma,G%xq,G%xr,G%yq,G%yr,G%Ji,B%mgx  ,B%mix-1,B%miy  ,B%piy  ,'L')
        if (B%sideR) call rates_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
          F%F,F%DF,M%rho,M%G,M%M,M%gamma,G%xq,G%xr,G%yq,G%yr,G%Ji,B%pix+1,B%pgx  ,B%miy  ,B%piy  ,'R')
        if (B%sideB) call rates_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
          F%F,F%DF,M%rho,M%G,M%M,M%gamma,G%xq,G%xr,G%yq,G%yr,G%Ji,B%mix  ,B%pix  ,B%mgy  ,B%miy-1,'B')
        if (B%sideT) call rates_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
          F%F,F%DF,M%rho,M%G,M%M,M%gamma,G%xq,G%xr,G%yq,G%yr,G%Ji,B%mix  ,B%pix  ,B%piy+1,B%pgy  ,'T')

        if (B%sideL.and.B%sideB) call rates_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
          F%F,F%DF,M%rho,M%G,M%M,M%gamma,G%xq,G%xr,G%yq,G%yr,G%Ji,B%mgx  ,B%mix-1,B%mgy,B%miy-1,'LB')
        if (B%sideR.and.B%sideB) call rates_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
          F%F,F%DF,M%rho,M%G,M%M,M%gamma,G%xq,G%xr,G%yq,G%yr,G%Ji,B%pix+1,B%pgx  ,B%mgy,B%miy-1,'RB')
        if (B%sideL.and.B%sideT) call rates_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
          F%F,F%DF,M%rho,M%G,M%M,M%gamma,G%xq,G%xr,G%yq,G%yr,G%Ji,B%mgx  ,B%mix-1,B%piy+1,B%pgy,'LT')
        if (B%sideR.and.B%sideT) call rates_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
          F%F,F%DF,M%rho,M%G,M%M,M%gamma,G%xq,G%xr,G%yq,G%yr,G%Ji,B%pix+1,B%pgx  ,B%piy+1,B%pgy,'RT')
      end if

       ! use plane strain relation to set d(szz)/dt
       
       call rates_szz(C%mx,C%my,F%DF,M%nu,B%mx,B%px,B%my,B%py)

    case(3)

       ! calculate rates at interior points
       
       call rates_interior_mode3(C%mbx,C%mby,C%mx,C%my, &
            F%F,F%DF,M%rho,M%G,G%xq,G%xr,G%yq,G%yr,G%Ji,B%mix,B%pix,B%miy,B%piy)

       ! calculate rates at boundary points
       
       if (B%sideL) call rates_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
            F%F,F%DF,M%rho,M%G,G%xq,G%xr,G%yq,G%yr,G%Ji,B%mgx  ,B%mix-1,B%miy  ,B%piy  ,'L')
       if (B%sideR) call rates_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
            F%F,F%DF,M%rho,M%G,G%xq,G%xr,G%yq,G%yr,G%Ji,B%pix+1,B%pgx  ,B%miy  ,B%piy  ,'R')
       if (B%sideB) call rates_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
            F%F,F%DF,M%rho,M%G,G%xq,G%xr,G%yq,G%yr,G%Ji,B%mix  ,B%pix  ,B%mgy  ,B%miy-1,'B')
       if (B%sideT) call rates_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
            F%F,F%DF,M%rho,M%G,G%xq,G%xr,G%yq,G%yr,G%Ji,B%mix  ,B%pix  ,B%piy+1,B%pgy  ,'T')

       if (B%sideL.and.B%sideB) call rates_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
            F%F,F%DF,M%rho,M%G,G%xq,G%xr,G%yq,G%yr,G%Ji,B%mgx  ,B%mix-1,B%mgy,B%miy-1,'LB')
       if (B%sideR.and.B%sideB) call rates_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
            F%F,F%DF,M%rho,M%G,G%xq,G%xr,G%yq,G%yr,G%Ji,B%pix+1,B%pgx  ,B%mgy,B%miy-1,'RB')
       if (B%sideL.and.B%sideT) call rates_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
            F%F,F%DF,M%rho,M%G,G%xq,G%xr,G%yq,G%yr,G%Ji,B%mgx  ,B%mix-1,B%piy+1,B%pgy,'LT')
       if (B%sideR.and.B%sideT) call rates_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
            F%F,F%DF,M%rho,M%G,G%xq,G%xr,G%yq,G%yr,G%Ji,B%pix+1,B%pgx  ,B%piy+1,B%pgy,'RT')

    end select
  end subroutine set_rates

  subroutine set_rates_dissipation(B,G,F,M,C,mode,CFL_dt)

    use grid, only : grid_type,block_grid
    use fields, only : fields_type
    use rates, only : &
         rates_dissipation_interior_mode2,rates_dissipation_boundary_mode2, &
         rates_dissipation_interior_mode3,rates_dissipation_boundary_mode3, &
         rates_szz
    use material, only : block_material
    use mpi_routines2d, only : cartesian

    implicit none

    type(block_grid),intent(in) :: B
    type(grid_type),intent(in) :: G
    type(fields_type),intent(inout) :: F
    type(block_material),intent(in) :: M
    type(cartesian),intent(in) :: C
    integer,intent(in) :: mode
    real,intent(in) :: CFL_dt

    if (B%skip) return ! process has no cells in this block

    select case(mode)

    case(2)

       ! calculate rates at interior points
       
       call rates_dissipation_interior_mode2(C%mbx,C%mby,C%mx,C%my, &
            F%F,F%DF,M%rho,M%G,M%M,M%gamma,CFL_dt, &
            G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mix,B%pix,B%miy,B%piy)

       ! calculate rates at boundary points
       
       if (B%sideL) call rates_dissipation_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
            F%F,F%DF,M%rho,M%G,M%M,M%gamma,CFL_dt, &
            G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mgx  ,B%mix-1,B%miy  ,B%piy  ,'L')
       if (B%sideR) call rates_dissipation_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
            F%F,F%DF,M%rho,M%G,M%M,M%gamma,CFL_dt, &
            G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%pix+1,B%pgx  ,B%miy  ,B%piy  ,'R')
       if (B%sideB) call rates_dissipation_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
            F%F,F%DF,M%rho,M%G,M%M,M%gamma,CFL_dt, &
            G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mix  ,B%pix  ,B%mgy  ,B%miy-1,'B')
       if (B%sideT) call rates_dissipation_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
            F%F,F%DF,M%rho,M%G,M%M,M%gamma,CFL_dt, &
            G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mix  ,B%pix  ,B%piy+1,B%pgy  ,'T')

       if (B%sideL.and.B%sideB) call rates_dissipation_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
            F%F,F%DF,M%rho,M%G,M%M,M%gamma,CFL_dt, &
            G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mgx  ,B%mix-1,B%mgy,B%miy-1,'LB')
       if (B%sideR.and.B%sideB) call rates_dissipation_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
            F%F,F%DF,M%rho,M%G,M%M,M%gamma,CFL_dt, &
            G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%pix+1,B%pgx  ,B%mgy,B%miy-1,'RB')
       if (B%sideL.and.B%sideT) call rates_dissipation_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
            F%F,F%DF,M%rho,M%G,M%M,M%gamma,CFL_dt, &
            G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mgx  ,B%mix-1,B%piy+1,B%pgy,'LT')
       if (B%sideR.and.B%sideT) call rates_dissipation_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
            F%F,F%DF,M%rho,M%G,M%M,M%gamma,CFL_dt, &
            G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%pix+1,B%pgx  ,B%piy+1,B%pgy,'RT')

       ! use plane strain relation to set d(szz)/dt
       
       call rates_szz(C%mx,C%my,F%DF,M%nu,B%mx,B%px,B%my,B%py)

    case(3)

       ! calculate rates at interior points
       
       call rates_dissipation_interior_mode3(C%mbx,C%mby,C%mx,C%my, &
            F%F,F%DF,M%rho,M%G,CFL_dt, &
            G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mix,B%pix,B%miy,B%piy)

       ! calculate rates at boundary points
       
       if (B%sideL) call rates_dissipation_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
            F%F,F%DF,M%rho,M%G,CFL_dt, &
            G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mgx  ,B%mix-1,B%miy  ,B%piy  ,'L')
       if (B%sideR) call rates_dissipation_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
            F%F,F%DF,M%rho,M%G,CFL_dt, &
            G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%pix+1,B%pgx  ,B%miy  ,B%piy  ,'R')
       if (B%sideB) call rates_dissipation_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
            F%F,F%DF,M%rho,M%G,CFL_dt, &
            G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mix  ,B%pix  ,B%mgy  ,B%miy-1,'B')
       if (B%sideT) call rates_dissipation_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
            F%F,F%DF,M%rho,M%G,CFL_dt, &
            G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mix  ,B%pix  ,B%piy+1,B%pgy  ,'T')

       if (B%sideL.and.B%sideB) call rates_dissipation_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
            F%F,F%DF,M%rho,M%G,CFL_dt, &
            G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mgx  ,B%mix-1,B%mgy,B%miy-1,'LB')
       if (B%sideR.and.B%sideB) call rates_dissipation_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
            F%F,F%DF,M%rho,M%G,CFL_dt, &
            G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%pix+1,B%pgx  ,B%mgy,B%miy-1,'RB')
       if (B%sideL.and.B%sideT) call rates_dissipation_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
            F%F,F%DF,M%rho,M%G,CFL_dt, &
            G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mgx  ,B%mix-1,B%piy+1,B%pgy,'LT')
       if (B%sideR.and.B%sideT) call rates_dissipation_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
            F%F,F%DF,M%rho,M%G,CFL_dt, &
            G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%pix+1,B%pgx  ,B%piy+1,B%pgy,'RT')

    end select

  end subroutine set_rates_dissipation


  subroutine set_rates_dissipation_xonly(B,G,F,M,C,mode,CFL_dt)

    use grid, only : grid_type,block_grid
    use fields, only : fields_type
    use rates, only : &
         rates_dissipation_xonly_interior_mode3,rates_dissipation_xonly_boundary_mode3
    use material, only : block_material
    use mpi_routines2d, only : cartesian
    use io, only : error

    implicit none

    type(block_grid),intent(in) :: B
    type(grid_type),intent(in) :: G
    type(fields_type),intent(inout) :: F
    type(block_material),intent(in) :: M
    type(cartesian),intent(in) :: C
    integer,intent(in) :: mode
    real,intent(in) :: CFL_dt

    if (B%skip) return ! process has no cells in this block

    select case(mode)

    case(2)

       call error('No routines for x-only artificial dissipation','set_rates_dissipation_xonly')

    case(3)

       ! calculate rates at interior points
       
       call rates_dissipation_xonly_interior_mode3(C%mbx,C%mby,C%mx,C%my, &
            F%F,F%DF,M%rho,M%G,CFL_dt, &
            G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mix,B%pix,B%miy,B%piy)

       ! calculate rates at boundary points
       
       if (B%sideL) call rates_dissipation_xonly_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
            F%F,F%DF,M%rho,M%G,CFL_dt, &
            G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mgx  ,B%mix-1,B%miy  ,B%piy  ,'L')
       if (B%sideR) call rates_dissipation_xonly_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
            F%F,F%DF,M%rho,M%G,CFL_dt, &
            G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%pix+1,B%pgx  ,B%miy  ,B%piy  ,'R')
       if (B%sideB) call rates_dissipation_xonly_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
            F%F,F%DF,M%rho,M%G,CFL_dt, &
            G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mix  ,B%pix  ,B%mgy  ,B%miy-1,'B')
       if (B%sideT) call rates_dissipation_xonly_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
            F%F,F%DF,M%rho,M%G,CFL_dt, &
            G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mix  ,B%pix  ,B%piy+1,B%pgy  ,'T')

       if (B%sideL.and.B%sideB) call rates_dissipation_xonly_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
            F%F,F%DF,M%rho,M%G,CFL_dt, &
            G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mgx  ,B%mix-1,B%mgy,B%miy-1,'LB')
       if (B%sideR.and.B%sideB) call rates_dissipation_xonly_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
            F%F,F%DF,M%rho,M%G,CFL_dt, &
            G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%pix+1,B%pgx  ,B%mgy,B%miy-1,'RB')
       if (B%sideL.and.B%sideT) call rates_dissipation_xonly_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
            F%F,F%DF,M%rho,M%G,CFL_dt, &
            G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mgx  ,B%mix-1,B%piy+1,B%pgy,'LT')
       if (B%sideR.and.B%sideT) call rates_dissipation_xonly_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
            F%F,F%DF,M%rho,M%G,CFL_dt, &
            G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%pix+1,B%pgx  ,B%piy+1,B%pgy,'RT')

    end select

  end subroutine set_rates_dissipation_xonly


  subroutine set_rates_SAT(B,G,F,BF,mode,Ks,Kp)

    use grid, only : block_grid,grid_type
    use fields, only : fields_type,block_fields

    implicit none

    type(block_grid),intent(in) :: B
    type(grid_type),intent(in) :: G
    type(fields_type),intent(inout) :: F
    type(block_fields),intent(in) :: BF
    integer,intent(in) :: mode
    real,intent(in) :: Ks,Kp

    integer :: i,j
    real :: dx,dy

    ! use simultaneous approximation term to enforce boundary conditions

    if (B%skip) return ! process has no cells in this block

    ! zero-dimensional case
    
    if (B%nx==1.and.B%ny==1) return
    
    ! one- and two-dimensional cases

    if (B%sideL.and.B%nx/=1) then
       i = B%mgx
       do j = B%my,B%py
          dx = abs(G%J(i,j)/sqrt(G%xr(i,j)**2+G%yr(i,j)**2))
          call SAT_term(F%DF(i,j,:),BF%bndFL%F(j,:),F%F(i,j,:), &
               mode,dx,B%bndL%n(j,:),Ks,Kp)
       end do
    end if
    
    if (B%sideR.and.B%nx/=1) then
       i = B%pgx
       do j = B%my,B%py
          dx = abs(G%J(i,j)/sqrt(G%xr(i,j)**2+G%yr(i,j)**2))
          call SAT_term(F%DF(i,j,:),BF%bndFR%F(j,:),F%F(i,j,:), &
               mode,dx,B%bndR%n(j,:),Ks,Kp)
       end do
    end if
    
    if (B%sideB.and.B%ny/=1) then
       j = B%mgy
       do i = B%mx,B%px
          dy = abs(G%J(i,j)/sqrt(G%xq(i,j)**2+G%yq(i,j)**2))
          call SAT_term(F%DF(i,j,:),BF%bndFB%F(i,:),F%F(i,j,:), &
               mode,dy,B%bndB%n(i,:),Ks,Kp)
       end do
    end if
    
    if (B%sideT.and.B%ny/=1) then
       j = B%pgy
       do i = B%mx,B%px
          dy = abs(G%J(i,j)/sqrt(G%xq(i,j)**2+G%yq(i,j)**2))
          call SAT_term(F%DF(i,j,:),BF%bndFT%F(i,:),F%F(i,j,:), &
               mode,dy,B%bndT%n(i,:),Ks,Kp)
       end do
    end if
 
  end subroutine set_rates_SAT


  subroutine SAT_term(DF,Fbnd,F,mode,h,normal,Ks,Kp)

    implicit none

    real,dimension(:),intent(inout) :: DF
    real,dimension(:),intent(in) :: Fbnd,F,normal
    integer,intent(in) :: mode
    real,intent(in) :: h,Ks,Kp

    real,dimension(size(F)) :: Fs,Fp

    select case(mode)
    case(2)
       call split_sp(Fbnd-F,normal,Fs,Fp)
       DF = DF+Ks/h*Fs
       DF = DF+Kp/h*Fp
    case(3)
       DF = DF+Ks/h*(Fbnd-F)
    end select

  end subroutine SAT_term


  subroutine split_sp(F,normal,Fs,Fp)

    use fields, only : rotate_fields_xy2nt,rotate_fields_nt2xy

    implicit none

    real,dimension(:),intent(in) :: F,normal
    real,dimension(:),intent(out) :: Fs,Fp

    real :: vn,vt,snn,snt,stt,szz

    call rotate_fields_xy2nt(F,normal,vt,vn,stt,snt,snn,szz)
    call rotate_fields_nt2xy(Fs,normal,vt ,0d0,0d0,snt,0d0,0d0)
    call rotate_fields_nt2xy(Fp,normal,0d0,vn ,stt,0d0,snn,szz)

  end subroutine split_sp


  subroutine set_source(B,G,F,M,C,t,mode,iblock)

    use io, only : error
    use grid, only : grid_type,block_grid
    use fields, only : fields_type, mms_sin, inplane_fault_mms
    use material, only : block_material
    use mpi_routines2d, only : cartesian

    implicit none

    type(block_grid),intent(in) :: B
    type(grid_type),intent(in) :: G
    type(fields_type),intent(inout) :: F
    type(block_material),intent(in) :: M
    type(cartesian),intent(in) :: C
    integer,intent(in) :: mode,iblock
    real,intent(in) :: t

    real :: x,y
    integer :: i,j

    if (B%skip) return ! process has no cells in this block

    select case(F%problem)
    case('mms-sin')
      do j = B%my,B%py
        do i = B%mx,B%px
          x = G%x(i,j)
          y = G%y(i,j)
          F%DF(i,j,1) = F%DF(i,j,1) + mms_sin(x,y,t,iblock,'s_vx')
          F%DF(i,j,2) = F%DF(i,j,2) + mms_sin(x,y,t,iblock,'s_vy')
          F%DF(i,j,3) = F%DF(i,j,3) + mms_sin(x,y,t,iblock,'s_sxx')
          F%DF(i,j,4) = F%DF(i,j,4) + mms_sin(x,y,t,iblock,'s_sxy')
          F%DF(i,j,5) = F%DF(i,j,5) + mms_sin(x,y,t,iblock,'s_syy')
          F%DF(i,j,6) = F%DF(i,j,6) + mms_sin(x,y,t,iblock,'s_szz')
        end do
      end do
    case('inplane-fault-mms')
      do j = B%my,B%py
        do i = B%mx,B%px
          x = G%x(i,j)
          y = G%y(i,j)
          F%DF(i,j,1) = F%DF(i,j,1) + inplane_fault_mms(x,y,t,iblock,'s_vx')
          F%DF(i,j,2) = F%DF(i,j,2) + inplane_fault_mms(x,y,t,iblock,'s_vy')
          F%DF(i,j,3) = F%DF(i,j,3) + inplane_fault_mms(x,y,t,iblock,'s_sxx')
          F%DF(i,j,4) = F%DF(i,j,4) + inplane_fault_mms(x,y,t,iblock,'s_sxy')
          F%DF(i,j,5) = F%DF(i,j,5) + inplane_fault_mms(x,y,t,iblock,'s_syy')
          F%DF(i,j,6) = F%DF(i,j,6) + inplane_fault_mms(x,y,t,iblock,'s_szz')
        end do
      end do
    end select

  end subroutine set_source


  subroutine time_step_interseismic(D,z,H,vp,dt,outlist)

    use domain, only : domain_type,exchange_fields_edges,interior_to_edges
    use fields, only : cycle_stress_fields,exchange_fields
    use boundaries, only : apply_bc
    use fault, only : couple_blocks, cycle_stress_fault
    use output, only : output_list,write_output

    implicit none

    type(domain_type),intent(inout) :: D
    real,intent(in) :: z,H,vp,dt
    type(output_list),intent(inout) :: outlist

    integer :: i,im,ip

    call write_output(outlist,D)

    D%t = D%t+dt

    do i = 1,D%nblocks
       call cycle_stress_fields(D%F,D%G,D%B(i)%F,D%B(i)%G,D%B(i)%M,D%mode,z,H,vp,dt)
       call interior_to_edges(D%B(i)%G,D%F,D%B(i)%F)
       call apply_bc(D%B(i)%G,D%B(i)%F,D%B(i)%B,D%B(i)%M,D%mode,D%t,i)
    end do

    call exchange_fields(D%C,D%F)
   
    do i = 1,D%nifaces
       im = D%I(i)%iblockm
       ip = D%I(i)%iblockp
       
       call cycle_stress_fault(D%I(i)%FR)
       call exchange_fields_edges(D%I(i),D%C,D%B(im)%F,D%B(ip)%F)
       call couple_blocks(D%I(i),D%B(im)%F,D%B(ip)%F,D%B(im)%M,D%B(ip)%M, &
            D%C,D%mode,D%t,initialize=.true.,dt_in=dt)
    end do

  end subroutine time_step_interseismic

end module time_step
