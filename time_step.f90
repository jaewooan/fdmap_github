module time_step

  implicit none

  type :: RK_type
     character(256) :: method
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


  subroutine time_step_LS(D,dt,RK,outlist,final_step,solid)
    ! low storage Runge-Kutta methods

    use domain, only : domain_type,prepare_edges,enforce_edge_conditions
    use fields, only : exchange_fields,update_fields_peak
    use energy, only : energy_interior,energy_block
    use plastic, only : update_fields_plastic
    use interfaces, only : update_fields_iface_implicit
    use output, only : output_list,write_output
    use io, only : error

    implicit none

    type(domain_type),intent(inout) :: D
    real,intent(in) :: dt
    type(RK_type),intent(in) :: RK
    type(output_list),intent(inout) :: outlist
    logical,intent(in) :: final_step,solid

    integer :: i,im,ip,stage
    real :: t0

    ! initial time
    
    t0 = D%t

    ! loop over RK stages

    do stage = 1,RK%nstage

       ! multiply rates by RK coefficient A

       call scale_rates_all(D,RK%A(stage),solid)

       ! enforce boundary and interface conditions (set hat variables)
       ! also, set rates on interfaces

       call prepare_edges(D)
       call enforce_edge_conditions(D)
       
       ! set rates (except interfaces)

       call set_rates_all(D,solid)

       ! compute energy

       if (stage==1) call energy_interior(D)

       ! output fields (needs to be done here after setting rates, 
       ! so slip velocity and related fields are correct)
       
       if (stage==1 .and. D%F%peak) call update_fields_peak(D%C,D%F,D%mode)

       if (stage==1) call write_output(outlist,D)

       ! return if final time step

       if (final_step) return

       ! update fields

       D%t = t0+RK%C(stage)*dt
       call update_fields_all(D,RK%B(stage)*dt,solid)
       
       ! exchange interior fields between processes
       ! (stage==RK%nstage case is handled below after operator splitting update)

       if (stage/=RK%nstage) call exchange_fields(D%C,D%F)

    end do

    ! operator splitting: implicit updates of stiff terms or DAEs

    ! implicit Euler integration of plasticity
       
    do i = 1,D%nblocks
       call update_fields_plastic(D%B(i)%G,D%G,D%F,D%B(i)%F,D%B(i)%M,D%E,D%mode,dt)
    end do

    ! implicit Euler integration of interface physics
 
    do i = 1,D%nifaces
       im = D%I(i)%iblockm
       ip = D%I(i)%iblockp
       call update_fields_iface_implicit(D%I(i),D%B(im)%F,D%B(ip)%F,D%t,dt)
    end do

    ! exchange fields between processes

    call exchange_fields(D%C,D%F)

  end subroutine time_step_LS


  subroutine scale_rates_all(D,A,solid)

    use domain, only : domain_type
    use fields, only : scale_rates_interior,scale_rates_boundary
    use interfaces, only : scale_rates_iface
    use energy, only : scale_rates_energy

    implicit none

    type(domain_type),intent(inout) :: D
    real,intent(in) :: A
    logical,intent(in) :: solid

    integer :: i

    ! multiply interior rates by RK coefficient A
    
    if (solid) call scale_rates_interior(D%C,D%F,A)

    ! multiply interface rates by RK coefficient A
    
    do i = 1,D%nifaces
       call scale_rates_iface(D%I(i),A)
    end do
    
    ! multiply energy dissipation and boundary displacement rates by RK coefficient A
    
    do i = 1,D%nblocks
       call scale_rates_energy(  D%B(i)%G,D%B(i)%F,A)
       call scale_rates_boundary(D%B(i)%G,D%B(i)%F,A)
    end do

  end subroutine scale_rates_all


  subroutine set_rates_all(D,solid)

    use domain, only : domain_type
    use energy, only : set_rates_energy
    use fields, only : set_rates_boundary,set_rates_displacement
    use source, only : set_source

    implicit none

    type(domain_type),intent(inout) :: D
    logical,intent(in) :: solid

    integer :: i

    ! SAT penalty term
    
    do i = 1,D%nblocks
       call set_rates_SAT(D%B(i)%G,D%F,D%B(i)%F,D%mode)
    end do
    
    ! interior rates
    
    if (solid) then ! can turn off this time-consuming step
       
       do i = 1,D%nblocks
          if (D%B(i)%dissipation) then
             call set_rates_dissipation(D%B(i)%G,D%G,D%F,D%B(i)%M,D%E,D%C,D%mode,D%B(i)%Cdiss)
          else
             call set_rates(            D%B(i)%G,D%G,D%F,D%B(i)%M,D%E,D%C,D%mode)
          end if
       end do

       call set_rates_displacement(D%C,D%F)
       
    end if
    
    ! source terms
    
    do i = 1,D%nblocks
       call set_source(D%B(i)%G,D%G,D%F,D%B(i)%M,D%S,D%t,D%mode,i)
    end do
    
    ! energy dissipation and boundary displacement rates
    
    do i = 1,D%nblocks
       call set_rates_energy(  D%B(i)%G,D%B(i)%F,D%mode)
       call set_rates_boundary(D%B(i)%G,D%B(i)%F,D%F%nU)
    end do
    
  end subroutine set_rates_all


  subroutine update_fields_all(D,Bdt,solid)

    use domain, only : domain_type
    use fields, only : update_fields_interior,update_fields_boundary
    use interfaces, only : update_fields_iface
    use energy, only : update_energy

    implicit none

    type(domain_type),intent(inout) :: D
    real,intent(in) :: Bdt
    logical,intent(in) :: solid

    integer :: i

    ! update interior fields
    
    if (solid) call update_fields_interior(D%C,D%F,Bdt)
    
    ! update iface fields

    do i = 1,D%nifaces
       call update_fields_iface(D%I(i),Bdt)
    end do
    
    ! update energy and boundary displacement
    
    do i = 1,D%nblocks
       call update_energy(         D%B(i)%G,D%B(i)%F,Bdt)
       call update_fields_boundary(D%B(i)%G,D%B(i)%F,Bdt)
    end do
    
  end subroutine update_fields_all


  subroutine set_rates(B,G,F,M,E,C,mode)

    use grid, only : grid_type,block_grid
    use fields, only : fields_type
    use rates, only : &
         rates_interior_mode2,rates_boundary_mode2, &
         rates_interior_mode3,rates_boundary_mode3, &
         rates_szz,&
         rates_interior_mode2_pml,rates_boundary_mode2_pml
    use rates_heterogeneous, only : &
         rates_h_interior_mode2,rates_h_boundary_mode2, &
         rates_h_interior_mode3,rates_h_boundary_mode3, &
         rates_h_szz
    use material, only : block_material,elastic_type
    use mpi_routines2d, only : cartesian

    implicit none

    type(block_grid),intent(in) :: B
    type(grid_type),intent(in) :: G
    type(fields_type),intent(inout) :: F
    type(block_material),intent(in) :: M
    type(elastic_type),intent(in) :: E
    type(cartesian),intent(in) :: C
    integer,intent(in) :: mode

    logical :: heterogeneous

    heterogeneous = allocated(E%rho)

    if (B%skip) return ! process has no cells in this block

    select case(mode)

    case(2)

      if(M%pml%pmlx .or. M%pml%pmly) then
        ! calculate rates at interior points

        call rates_interior_mode2_pml(C%mbx,C%mby,C%mx,C%my,C%mx,C%my,&
          F%F,F%DF,F%Wx,F%DWx,F%Wy,F%DWy,&
          M%rho,M%G,M%M,M%gamma,&
          M%pml,&
          G%x,G%y,&
          G%xq,G%xr,G%yq,G%yr,G%Ji,B%mix,B%pix,B%miy,B%piy)

        ! calculate rates at boundary points

        if (B%sideL) call rates_boundary_mode2_pml(C%mbx,C%mby,C%mx,C%my,C%mx,C%my,&
          F%F,F%DF,F%Wx,F%DWx,F%Wy,F%DWy,&
          M%rho,M%G,M%M,M%gamma,&
          M%pml,&
          G%x,G%y,&
          G%xq,G%xr,G%yq,G%yr,G%Ji,B%mgx  ,B%mix-1,B%miy  ,B%piy,&
          'L')
        if (B%sideR) call rates_boundary_mode2_pml(C%mbx,C%mby,C%mx,C%my,C%mx,C%my,&
          F%F,F%DF,F%Wx,F%DWx,F%Wy,F%DWy,&
          M%rho,M%G,M%M,M%gamma,&
          M%pml,&
          G%x,G%y,&
          G%xq,G%xr,G%yq,G%yr,G%Ji,B%pix+1,B%pgx  ,B%miy  ,B%piy,&
          'R')
        if (B%sideB) call rates_boundary_mode2_pml(C%mbx,C%mby,C%mx,C%my,C%mx,C%my,&
          F%F,F%DF,F%Wx,F%DWx,F%Wy,F%DWy,&
          M%rho,M%G,M%M,M%gamma,&
          M%pml,&
          G%x,G%y,&
          G%xq,G%xr,G%yq,G%yr,G%Ji,B%mix  ,B%pix  ,B%mgy,B%miy-1,&
          'B')
        if (B%sideT) call rates_boundary_mode2_pml(C%mbx,C%mby,C%mx,C%my,C%mx,C%my,&
          F%F,F%DF,F%Wx,F%DWx,F%Wy,F%DWy,&
          M%rho,M%G,M%M,M%gamma,&
          M%pml,&
          G%x,G%y,&
          G%xq,G%xr,G%yq,G%yr,G%Ji,B%mix  ,B%pix  ,B%piy+1,B%pgy,&
          'T')

        if (B%sideL.and.B%sideB) call rates_boundary_mode2_pml(C%mbx,C%mby,C%mx,C%my,C%mx,C%my,&
          F%F,F%DF,F%Wx,F%DWx,F%Wy,F%DWy,&
          M%rho,M%G,M%M,M%gamma,&
          M%pml,&
          G%x,G%y,&
          G%xq,G%xr,G%yq,G%yr,G%Ji,B%mgx  ,B%mix-1,B%mgy,B%miy-1,&
          'LB')
        if (B%sideR.and.B%sideB) call rates_boundary_mode2_pml(C%mbx,C%mby,C%mx,C%my,C%mx,C%my,&
          F%F,F%DF,F%Wx,F%DWx,F%Wy,F%DWy,&
          M%rho,M%G,M%M,M%gamma,&
          M%pml,&
          G%x,G%y,&
          G%xq,G%xr,G%yq,G%yr,G%Ji,B%pix+1,B%pgx  ,B%mgy,B%miy-1,&
          'RB')
        if (B%sideL.and.B%sideT) call rates_boundary_mode2_pml(C%mbx,C%mby,C%mx,C%my,C%mx,C%my,&
          F%F,F%DF,F%Wx,F%DWx,F%Wy,F%DWy,&
          M%rho,M%G,M%M,M%gamma,&
          M%pml,&
          G%x,G%y,&
          G%xq,G%xr,G%yq,G%yr,G%Ji,B%mgx  ,B%mix-1,B%piy+1,B%pgy,&
          'LT')
        if (B%sideR.and.B%sideT) call rates_boundary_mode2_pml(C%mbx,C%mby,C%mx,C%my,C%mx,C%my,&
          F%F,F%DF,F%Wx,F%DWx,F%Wy,F%DWy,&
          M%rho,M%G,M%M,M%gamma,&
          M%pml,&
          G%x,G%y,&
          G%xq,G%xr,G%yq,G%yr,G%Ji,B%pix+1,B%pgx  ,B%piy+1,B%pgy,&
          'RT')

      else

         if (heterogeneous) then

            ! calculate rates at interior points
            
            call rates_h_interior_mode2(C%mbx,C%mby,C%mx,C%my, &
                 F%F,F%DF,E%rho,E%G,E%M,E%gamma,G%xq,G%xr,G%yq,G%yr,G%Ji,B%mix,B%pix,B%miy,B%piy)
            
            ! calculate rates at boundary points
            
            if (B%sideL) call rates_h_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
                 F%F,F%DF,E%rho,E%G,E%M,E%gamma,G%xq,G%xr,G%yq,G%yr,G%Ji,B%mgx  ,B%mix-1,B%miy  ,B%piy  ,'L')
            if (B%sideR) call rates_h_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
                 F%F,F%DF,E%rho,E%G,E%M,E%gamma,G%xq,G%xr,G%yq,G%yr,G%Ji,B%pix+1,B%pgx  ,B%miy  ,B%piy  ,'R')
            if (B%sideB) call rates_h_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
                 F%F,F%DF,E%rho,E%G,E%M,E%gamma,G%xq,G%xr,G%yq,G%yr,G%Ji,B%mix  ,B%pix  ,B%mgy  ,B%miy-1,'B')
            if (B%sideT) call rates_h_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
                 F%F,F%DF,E%rho,E%G,E%M,E%gamma,G%xq,G%xr,G%yq,G%yr,G%Ji,B%mix  ,B%pix  ,B%piy+1,B%pgy  ,'T')
            
            if (B%sideL.and.B%sideB) call rates_h_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
                 F%F,F%DF,E%rho,E%G,E%M,E%gamma,G%xq,G%xr,G%yq,G%yr,G%Ji,B%mgx  ,B%mix-1,B%mgy,B%miy-1,'LB')
            if (B%sideR.and.B%sideB) call rates_h_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
                 F%F,F%DF,E%rho,E%G,E%M,E%gamma,G%xq,G%xr,G%yq,G%yr,G%Ji,B%pix+1,B%pgx  ,B%mgy,B%miy-1,'RB')
            if (B%sideL.and.B%sideT) call rates_h_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
                 F%F,F%DF,E%rho,E%G,E%M,E%gamma,G%xq,G%xr,G%yq,G%yr,G%Ji,B%mgx  ,B%mix-1,B%piy+1,B%pgy,'LT')
            if (B%sideR.and.B%sideT) call rates_h_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
                 F%F,F%DF,E%rho,E%G,E%M,E%gamma,G%xq,G%xr,G%yq,G%yr,G%Ji,B%pix+1,B%pgx  ,B%piy+1,B%pgy,'RT')

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

      end if

       ! use plane strain relation to set d(szz)/dt
       
      if (heterogeneous) then
         call rates_h_szz(C%mbx,C%mby,C%mx,C%my,F%DF,E%nu,B%mx,B%px,B%my,B%py)
      else
         call rates_szz(C%mx,C%my,F%DF,M%nu,B%mx,B%px,B%my,B%py)
      end if

    case(3)

       if (heterogeneous) then

          ! calculate rates at interior points
          
          call rates_h_interior_mode3(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,E%rho,E%G,G%xq,G%xr,G%yq,G%yr,G%Ji,B%mix,B%pix,B%miy,B%piy)
          
          ! calculate rates_h at boundary points
          
          if (B%sideL) call rates_h_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,E%rho,E%G,G%xq,G%xr,G%yq,G%yr,G%Ji,B%mgx  ,B%mix-1,B%miy  ,B%piy  ,'L')
          if (B%sideR) call rates_h_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,E%rho,E%G,G%xq,G%xr,G%yq,G%yr,G%Ji,B%pix+1,B%pgx  ,B%miy  ,B%piy  ,'R')
          if (B%sideB) call rates_h_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,E%rho,E%G,G%xq,G%xr,G%yq,G%yr,G%Ji,B%mix  ,B%pix  ,B%mgy  ,B%miy-1,'B')
          if (B%sideT) call rates_h_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,E%rho,E%G,G%xq,G%xr,G%yq,G%yr,G%Ji,B%mix  ,B%pix  ,B%piy+1,B%pgy  ,'T')
          
          if (B%sideL.and.B%sideB) call rates_h_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,E%rho,E%G,G%xq,G%xr,G%yq,G%yr,G%Ji,B%mgx  ,B%mix-1,B%mgy,B%miy-1,'LB')
          if (B%sideR.and.B%sideB) call rates_h_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,E%rho,E%G,G%xq,G%xr,G%yq,G%yr,G%Ji,B%pix+1,B%pgx  ,B%mgy,B%miy-1,'RB')
          if (B%sideL.and.B%sideT) call rates_h_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,E%rho,E%G,G%xq,G%xr,G%yq,G%yr,G%Ji,B%mgx  ,B%mix-1,B%piy+1,B%pgy,'LT')
          if (B%sideR.and.B%sideT) call rates_h_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,E%rho,E%G,G%xq,G%xr,G%yq,G%yr,G%Ji,B%pix+1,B%pgx  ,B%piy+1,B%pgy,'RT')

       else

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

       end if

    end select

  end subroutine set_rates


  subroutine set_rates_dissipation(B,G,F,M,E,C,mode,Cdiss)

    use grid, only : grid_type,block_grid
    use fields, only : fields_type
    use rates, only : &
         rates_dissipation_interior_mode2,rates_dissipation_boundary_mode2, &
         rates_dissipation_interior_mode3,rates_dissipation_boundary_mode3, &
         rates_szz
    use rates_heterogeneous, only : &
         rates_h_dissipation_interior_mode2,rates_h_dissipation_boundary_mode2, &
         rates_h_dissipation_interior_mode3,rates_h_dissipation_boundary_mode3, &
         rates_h_szz
    use material, only : block_material,elastic_type
    use mpi_routines2d, only : cartesian

    implicit none

    type(block_grid),intent(in) :: B
    type(grid_type),intent(in) :: G
    type(fields_type),intent(inout) :: F
    type(block_material),intent(in) :: M
    type(elastic_type),intent(in) :: E
    type(cartesian),intent(in) :: C
    integer,intent(in) :: mode
    real,intent(in) :: Cdiss

    logical :: heterogeneous

    heterogeneous = allocated(E%rho)

    if (B%skip) return ! process has no cells in this block

    select case(mode)

    case(2)

       if (heterogeneous) then

          ! calculate rates at interior points
          
          call rates_h_dissipation_interior_mode2(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,E%rho,E%G,E%M,E%gamma,Cdiss, &
               G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mix,B%pix,B%miy,B%piy)
          
          ! calculate rates at boundary points
          
          if (B%sideL) call rates_h_dissipation_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,E%rho,E%G,E%M,E%gamma,Cdiss, &
               G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mgx  ,B%mix-1,B%miy  ,B%piy  ,'L')
          if (B%sideR) call rates_h_dissipation_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,E%rho,E%G,E%M,E%gamma,Cdiss, &
               G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%pix+1,B%pgx  ,B%miy  ,B%piy  ,'R')
          if (B%sideB) call rates_h_dissipation_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,E%rho,E%G,E%M,E%gamma,Cdiss, &
               G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mix  ,B%pix  ,B%mgy  ,B%miy-1,'B')
          if (B%sideT) call rates_h_dissipation_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,E%rho,E%G,E%M,E%gamma,Cdiss, &
               G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mix  ,B%pix  ,B%piy+1,B%pgy  ,'T')
          
          if (B%sideL.and.B%sideB) call rates_h_dissipation_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,E%rho,E%G,E%M,E%gamma,Cdiss, &
               G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mgx  ,B%mix-1,B%mgy,B%miy-1,'LB')
          if (B%sideR.and.B%sideB) call rates_h_dissipation_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,E%rho,E%G,E%M,E%gamma,Cdiss, &
               G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%pix+1,B%pgx  ,B%mgy,B%miy-1,'RB')
          if (B%sideL.and.B%sideT) call rates_h_dissipation_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,E%rho,E%G,E%M,E%gamma,Cdiss, &
               G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mgx  ,B%mix-1,B%piy+1,B%pgy,'LT')
          if (B%sideR.and.B%sideT) call rates_h_dissipation_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,E%rho,E%G,E%M,E%gamma,Cdiss, &
               G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%pix+1,B%pgx  ,B%piy+1,B%pgy,'RT')
          
          ! use plane strain relation to set d(szz)/dt
          
          call rates_h_szz(C%mbx,C%mby,C%mx,C%my,F%DF,E%nu,B%mx,B%px,B%my,B%py)

       else

          ! calculate rates at interior points
          
          call rates_dissipation_interior_mode2(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,M%rho,M%G,M%M,M%gamma,Cdiss, &
               G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mix,B%pix,B%miy,B%piy)
          
          ! calculate rates at boundary points
          
          if (B%sideL) call rates_dissipation_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,M%rho,M%G,M%M,M%gamma,Cdiss, &
               G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mgx  ,B%mix-1,B%miy  ,B%piy  ,'L')
          if (B%sideR) call rates_dissipation_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,M%rho,M%G,M%M,M%gamma,Cdiss, &
               G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%pix+1,B%pgx  ,B%miy  ,B%piy  ,'R')
          if (B%sideB) call rates_dissipation_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,M%rho,M%G,M%M,M%gamma,Cdiss, &
               G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mix  ,B%pix  ,B%mgy  ,B%miy-1,'B')
          if (B%sideT) call rates_dissipation_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,M%rho,M%G,M%M,M%gamma,Cdiss, &
               G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mix  ,B%pix  ,B%piy+1,B%pgy  ,'T')
          
          if (B%sideL.and.B%sideB) call rates_dissipation_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,M%rho,M%G,M%M,M%gamma,Cdiss, &
               G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mgx  ,B%mix-1,B%mgy,B%miy-1,'LB')
          if (B%sideR.and.B%sideB) call rates_dissipation_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,M%rho,M%G,M%M,M%gamma,Cdiss, &
               G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%pix+1,B%pgx  ,B%mgy,B%miy-1,'RB')
          if (B%sideL.and.B%sideT) call rates_dissipation_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,M%rho,M%G,M%M,M%gamma,Cdiss, &
               G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mgx  ,B%mix-1,B%piy+1,B%pgy,'LT')
          if (B%sideR.and.B%sideT) call rates_dissipation_boundary_mode2(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,M%rho,M%G,M%M,M%gamma,Cdiss, &
               G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%pix+1,B%pgx  ,B%piy+1,B%pgy,'RT')
          
          ! use plane strain relation to set d(szz)/dt
          
          call rates_szz(C%mx,C%my,F%DF,M%nu,B%mx,B%px,B%my,B%py)

       end if

    case(3)

       if (heterogeneous) then

          ! calculate rates at interior points
          
          call rates_h_dissipation_interior_mode3(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,E%rho,E%G,Cdiss, &
               G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mix,B%pix,B%miy,B%piy)
          
          ! calculate rates at boundary points
          
          if (B%sideL) call rates_h_dissipation_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,E%rho,E%G,Cdiss, &
               G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mgx  ,B%mix-1,B%miy  ,B%piy  ,'L')
          if (B%sideR) call rates_h_dissipation_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,E%rho,E%G,Cdiss, &
               G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%pix+1,B%pgx  ,B%miy  ,B%piy  ,'R')
          if (B%sideB) call rates_h_dissipation_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,E%rho,E%G,Cdiss, &
               G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mix  ,B%pix  ,B%mgy  ,B%miy-1,'B')
          if (B%sideT) call rates_h_dissipation_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,E%rho,E%G,Cdiss, &
               G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mix  ,B%pix  ,B%piy+1,B%pgy  ,'T')
          
          if (B%sideL.and.B%sideB) call rates_h_dissipation_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,E%rho,E%G,Cdiss, &
               G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mgx  ,B%mix-1,B%mgy,B%miy-1,'LB')
          if (B%sideR.and.B%sideB) call rates_h_dissipation_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,E%rho,E%G,Cdiss, &
               G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%pix+1,B%pgx  ,B%mgy,B%miy-1,'RB')
          if (B%sideL.and.B%sideT) call rates_h_dissipation_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,E%rho,E%G,Cdiss, &
               G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mgx  ,B%mix-1,B%piy+1,B%pgy,'LT')
          if (B%sideR.and.B%sideT) call rates_h_dissipation_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,E%rho,E%G,Cdiss, &
               G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%pix+1,B%pgx  ,B%piy+1,B%pgy,'RT')

       else
          
          ! calculate rates at interior points
          
          call rates_dissipation_interior_mode3(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,M%rho,M%G,Cdiss, &
               G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mix,B%pix,B%miy,B%piy)
          
          ! calculate rates at boundary points
          
          if (B%sideL) call rates_dissipation_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,M%rho,M%G,Cdiss, &
               G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mgx  ,B%mix-1,B%miy  ,B%piy  ,'L')
          if (B%sideR) call rates_dissipation_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,M%rho,M%G,Cdiss, &
               G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%pix+1,B%pgx  ,B%miy  ,B%piy  ,'R')
          if (B%sideB) call rates_dissipation_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,M%rho,M%G,Cdiss, &
               G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mix  ,B%pix  ,B%mgy  ,B%miy-1,'B')
          if (B%sideT) call rates_dissipation_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,M%rho,M%G,Cdiss, &
               G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mix  ,B%pix  ,B%piy+1,B%pgy  ,'T')
          
          if (B%sideL.and.B%sideB) call rates_dissipation_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,M%rho,M%G,Cdiss, &
               G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mgx  ,B%mix-1,B%mgy,B%miy-1,'LB')
          if (B%sideR.and.B%sideB) call rates_dissipation_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,M%rho,M%G,Cdiss, &
               G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%pix+1,B%pgx  ,B%mgy,B%miy-1,'RB')
          if (B%sideL.and.B%sideT) call rates_dissipation_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,M%rho,M%G,Cdiss, &
               G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%mgx  ,B%mix-1,B%piy+1,B%pgy,'LT')
          if (B%sideR.and.B%sideT) call rates_dissipation_boundary_mode3(C%mbx,C%mby,C%mx,C%my, &
               F%F,F%DF,M%rho,M%G,Cdiss, &
               G%xq,G%xr,G%yq,G%yr,G%J,G%Ji,B%pix+1,B%pgx  ,B%piy+1,B%pgy,'RT')

       end if

    end select

  end subroutine set_rates_dissipation


  subroutine set_rates_SAT(B,F,BF,mode)

    use grid, only : block_grid
    use fields, only : fields_type,block_fields

    implicit none

    type(block_grid),intent(in) :: B
    type(fields_type),intent(inout) :: F
    type(block_fields),intent(in) :: BF
    integer,intent(in) :: mode

    integer :: i,j

    ! use simultaneous approximation term to enforce boundary conditions

    if (B%skip) return ! process has no cells in this block

    ! zero-dimensional case
    
    if (B%nx==1.and.B%ny==1) return
    
    ! one- and two-dimensional cases

    if (B%sideL.and.B%nx/=1) then
       i = B%mgx
       do j = B%my,B%py
          call SAT_term(F%DF(i,j,:),BF%bndFL%Fhat(j,:),BF%bndFL%F(j,:), &
               mode,B%bndL%n(j,:),BF%bndFL%M(j,4),BF%bndFL%M(j,5))
       end do
    end if
    
    if (B%sideR.and.B%nx/=1) then
       i = B%pgx
       do j = B%my,B%py
          call SAT_term(F%DF(i,j,:),BF%bndFR%Fhat(j,:),BF%bndFR%F(j,:), &
               mode,B%bndR%n(j,:),BF%bndFR%M(j,4),BF%bndFR%M(j,5))
       end do
    end if
    
    if (B%sideB.and.B%ny/=1) then
       j = B%mgy
       do i = B%mx,B%px
          call SAT_term(F%DF(i,j,:),BF%bndFB%Fhat(i,:),BF%bndFB%F(i,:), &
               mode,B%bndB%n(i,:),BF%bndFB%M(i,4),BF%bndFB%M(i,5))
       end do
    end if
    
    if (B%sideT.and.B%ny/=1) then
       j = B%pgy
       do i = B%mx,B%px
          call SAT_term(F%DF(i,j,:),BF%bndFT%Fhat(i,:),BF%bndFT%F(i,:), &
               mode,B%bndT%n(i,:),BF%bndFT%M(i,4),BF%bndFT%M(i,5))
       end do
    end if
 
  end subroutine set_rates_SAT


  subroutine SAT_term(DF,Fhat,F,mode,normal,Ks,Kp)

    implicit none

    real,dimension(:),intent(inout) :: DF
    real,dimension(:),intent(in) :: Fhat,F,normal
    integer,intent(in) :: mode
    real,intent(in) :: Ks,Kp

    real,dimension(size(F)) :: Fs,Fp

    select case(mode)
    case(2)
       call split_sp(Fhat-F,normal,Fs,Fp)
       DF = DF+Ks*Fs
       DF = DF+Kp*Fp
    case(3)
       DF = DF+Ks*(Fhat-F)
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


end module time_step
