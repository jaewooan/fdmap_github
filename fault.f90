module fault

  implicit none

contains


  subroutine couple_blocks(I,Fm,Fp,Mm,Mp,C,mode,t,initialize,dt_in)

    use boundaries, only : iface_type
    use fields, only : block_fields
    use material, only : block_material
    use mpi_routines2d, only : cartesian
    use thermpres, only : pressure_thermpres
    use erupt, only : pressure_erupt,set_rates_erupt
    use friction, only : state_rate,rupture_front

    implicit none

    type(iface_type),intent(inout) :: I
    type(block_fields),intent(inout) :: Fm,Fp
    type(block_material),intent(in) :: Mm,Mp
    type(cartesian),intent(in) :: C
    integer,intent(in) :: mode
    real,intent(in) :: t
    logical,intent(in) :: initialize
    real,intent(in),optional :: dt_in

    integer :: j
    real :: p,dt

    if (present(dt_in)) then
       dt = dt_in
    else
       dt = 0d0
    end if

    if (I%skip) return ! process has no cells adjacent to interface

    ! couple blocks, point by point

    do j = I%m,I%p

       ! set fluid pressure

       p = 0d0
       select case(I%coupling)
       case('friction')
          call pressure_thermpres(I%TP,j,p)
       case('eruption')
          call pressure_erupt(I%ER,j,p)
       end select

       ! enforce elasticity relations (and possibly friction)
       
       if (allocated(I%DE)) then
          select case(I%direction)
          case('x')
             call couple_points(I%FR,I%FR%V(j),I%FR%O(j),I%FR%S(j),I%FR%N(j),I%FR%Neff(j), &
                  I%FR%S0(j),I%FR%N0(j),I%FR%Dn(j),I%FR%D(j),I%FR%Psi(j),p, &
                  Fm%bndFR%F(j,:),Fp%bndFL%F(j,:),Fp%F0, &
                  I%nhat(j,:),I%coupling, &
                  j,I%x(j),I%y(j),t,Mm,Mp,mode,initialize,dt, &
                  I%DE(j),Fm%bndFR%DE(j),Fp%bndFL%DE(j))
          case('y')
             call couple_points(I%FR,I%FR%V(j),I%FR%O(j),I%FR%S(j),I%FR%N(j),I%FR%Neff(j), &
                  I%FR%S0(j),I%FR%N0(j),I%FR%Dn(j),I%FR%D(j),I%FR%Psi(j),p, &
                  Fm%bndFT%F(j,:),Fp%bndFB%F(j,:),Fp%F0, &
                  I%nhat(j,:),I%coupling, &
                  j,I%x(j),I%y(j),t,Mm,Mp,mode,initialize,dt, &
                  I%DE(j),Fm%bndFT%DE(j),Fp%bndFB%DE(j))
          end select
       else
          select case(I%direction)
          case('x')
             call couple_points(I%FR,I%FR%V(j),I%FR%O(j),I%FR%S(j),I%FR%N(j),I%FR%Neff(j), &
                  I%FR%S0(j),I%FR%N0(j),I%FR%Dn(j),I%FR%D(j),I%FR%Psi(j),p, &
                  Fm%bndFR%F(j,:),Fp%bndFL%F(j,:),Fp%F0, &
                  I%nhat(j,:),I%coupling, &
                  j,I%x(j),I%y(j),t,Mm,Mp,mode,initialize,dt)
          case('y')
             call couple_points(I%FR,I%FR%V(j),I%FR%O(j),I%FR%S(j),I%FR%N(j),I%FR%Neff(j), &
                  I%FR%S0(j),I%FR%N0(j),I%FR%Dn(j),I%FR%D(j),I%FR%Psi(j),p, &
                  Fm%bndFT%F(j,:),Fp%bndFB%F(j,:),Fp%F0, &
                  I%nhat(j,:),I%coupling, &
                  j,I%x(j),I%y(j),t,Mm,Mp,mode,initialize,dt)
          end select
       end if

       ! store rates for integration

       I%FR%DDs(j) = I%FR%DDs(j)+I%FR%V(j)
       I%FR%DDn(j) = I%FR%DDn(j)+I%FR%O(j)

       ! set rates for constitutive models (at point)
       
       select case(I%coupling)
       case('friction')
          I%FR%DPsi(j) = I%FR%DPsi(j)+state_rate(I%FR,I%FR%V(j),I%FR%Psi(j),j,I%x(j),I%y(j),t)
          call rupture_front(I%FR,I%FR%D(j),I%FR%V(j),t,I%FR%trup(j))
       end select

    end do

    ! set rates for constitutive models (at set of points)
    
    select case(I%coupling)
    case('eruption')
       if (.not.initialize) call set_rates_erupt(I%ER,C,I%comm_mp,t,I%m,I%p,I%y)
    end select

  end subroutine couple_blocks


  subroutine couple_points(FR,V,O,S,N,Neff,S0,N0,Dn,D,Psi,p,Fm,Fp,F0,&
       normal,coupling,i,x,y,t,Mm,Mp,mode,initialize,dt,DE,DEm,DEp)

    use friction, only : fr_type,load_stress
    use material, only : block_material
    use geometry, only : rotate_xy2nt

    implicit none

    type(fr_type),intent(in) :: FR
    real,intent(inout) :: V,O,S,N,S0,N0,Psi,Neff
    real,intent(in) :: Dn,D,p,dt
    real,dimension(:),intent(inout) :: Fm,Fp
    real,intent(in) :: F0(9),normal(2)
    character(*),intent(in) :: coupling
    integer,intent(in) :: i
    real,intent(in) :: x,y,t
    type(block_material),intent(in) :: Mm,Mp
    integer,intent(in) :: mode
    logical,intent(in) :: initialize
    real,intent(inout),optional :: DE,DEm,DEp

    real :: stt,snt,snn

    ! set load
    call load_stress(FR,x,y,t,S0,N0,normal,Mm,Mp,FR%bS0(i),FR%bN0(i))

    ! couple points

    if (present(DE)) then
       select case(mode)
       case(2)
          call couple_mode2(FR,V,O,S,N,Neff,S0,N0,Dn,D,Psi,p, &
               Fm,Fp,normal,coupling,i,x,y,t,Mm,Mp,initialize,dt,DE,DEm,DEp)
       case(3)
          ! normal stress contribution from resolving in-plane stress field onto fault
          call rotate_xy2nt(F0(4),F0(5),F0(7),stt,snt,snn,normal)
          N0 = N0-snn
          call couple_mode3(FR,V,O,S,N,S0,N0,D,Psi,p, &
               Fm,Fp,normal,coupling,i,x,y,t,Mm,Mp,initialize,dt,DE,DEm,DEp)
       end select
    else
       select case(mode)
       case(2)
          call couple_mode2(FR,V,O,S,N,Neff,S0,N0,Dn,D,Psi,p, &
               Fm,Fp,normal,coupling,i,x,y,t,Mm,Mp,initialize,dt)
       case(3)
          ! normal stress contribution from resolving in-plane stress field onto fault
          call rotate_xy2nt(F0(4),F0(5),F0(7),stt,snt,snn,normal)
          N0 = N0-snn
          call couple_mode3(FR,V,O,S,N,S0,N0,D,Psi,p, &
               Fm,Fp,normal,coupling,i,x,y,t,Mm,Mp,initialize,dt)
       end select
    end if

  end subroutine couple_points


  subroutine couple_mode3(FR,V,O,S,N,S0,N0,D,Psi,p,Fm,Fp, &
       normal,coupling,i,x,y,t,Mm,Mp,initialize,dt,DE,DEm,DEp)

    use friction, only : fr_type,initial_state,solve_friction
    use material, only : block_material
    use fields, only : rotate_fields_xy2nt,rotate_fields_nt2xy
    use io, only : error

    implicit none

    type(fr_type),intent(in) :: FR
    real,intent(inout) :: V,O,S,N,Psi
    real,intent(in) :: S0,N0,D,p
    real,dimension(3),intent(inout) :: Fm,Fp
    real,dimension(2),intent(in) :: normal
    character(*),intent(in) :: coupling
    integer,intent(in) :: i
    real,intent(in) :: x,y,t,dt
    type(block_material),intent(in) :: Mm,Mp
    logical,intent(in) :: initialize
    real,intent(inout),optional :: DE,DEm,DEp

    logical :: compressive
    real :: Nlock,Slock,phis,etas, &
         vzp,snzp,vzFDp,snzFDp,stzFDp, &
         vzm,snzm,vzFDm,snzFDm,stzFDm
    
    ! 1. snzp+Zsp*vzp = snzFDp+Zsp*vzFDp 
    !    (S-wave into interface from plus  side)
    ! 2. snzm-Zsm*vzm = snzFDm-Zsm*vzFDm 
    !    (S-wave into interface from minus side)
    ! 3. snzp = snzm
    !    (continuity of shear stress)
    !
    ! 1-3 yield equation of form snz = phis-etas*V with
    ! phis = etas*(snzFDp/Zsp+snzFDm/Zsm+vzFDp-vzFDm),
    ! V = vzp-vzm, and 1/etas = 1/Zsp+1/Zsm
    
    ! rotate into local normal and tangential coordinates

    call rotate_fields_xy2nt(Fp,normal,vzFDp,stzFDp,snzFDp)
    call rotate_fields_xy2nt(Fm,normal,vzFDm,stzFDm,snzFDm)

    ! apply interface conditions in local coordinates

    etas = 1d0/(Mp%Zsi+Mm%Zsi)
    phis = etas*(snzFDp*Mp%Zsi+snzFDm*Mm%Zsi+vzFDp-vzFDm)

    ! set V,O,S,N

    select case(coupling)

    case default

       call error('Invalid coupling (' // trim(coupling) &
            // ')','couple_mode3')

    case('friction')

       ! no friction law to solve in the case of RFSL, values come in correct!
       if(FR%friction_law .eq. 'RFSL') then
         if(initialize) then
           S = S0+phis
           N = N0-p
           V = 0d0
           O = 0d0
         end if
       else
         Slock = S0+phis
         Nlock = N0-p
         compressive = (Nlock>0d0)
         N = max(Nlock,0d0)
         O = 0d0 ! no opening in mode III, even for tensile normal stresses

         if (initialize) then
           V = vzFDp-vzFDm
           S = Slock-etas*V
           call initial_state(FR,Psi,V,S,N,i,x,y,dt)
         end if

         if (compressive) then ! apply friction law
           call solve_friction(FR,V,S,N,Slock,etas,D,Psi,i,x,y,t,.false.)
         else ! no shear resistance to slip (cannot allow opening in mode III)
           S = 0d0
           V = Slock/etas
         end if
       end if

    case('locked')

       V = 0d0
       O = 0d0

       S = S0+phis
       N = N0-p

    end select

    snzp = S-S0
    snzm = S-S0
    vzp = vzFDp+Mp%Zsi*(snzFDp-snzp)
    vzm = vzFDm-Mm%Zsi*(snzFDm-snzm)

    ! energy flux (positive when energy flows out of medium into fault)

    if (present(DE)) then
       DEm = DEm+snzm*vzm-0.25d0*Mm%Zsi*((snzFDm-snzm)+Mm%Zs*(vzFDm-vzm))**2
       DEp = DEp-snzp*vzp-0.25d0*Mp%Zsi*((snzFDp-snzp)-Mp%Zs*(vzFDp-vzp))**2
       DE = DE-S*V- &
            0.25d0*Mm%Zsi*((snzFDm-snzm)+Mm%Zs*(vzFDm-vzm))**2- &
            0.25d0*Mp%Zsi*((snzFDp-snzp)-Mp%Zs*(vzFDp-vzp))**2
    end if

    ! rotate back to x-y coordinates

    call rotate_fields_nt2xy(Fp,normal,vzp,stzFDp,snzp)
    call rotate_fields_nt2xy(Fm,normal,vzm,stzFDm,snzm)

  end subroutine couple_mode3


  subroutine couple_mode2(FR,V,O,S,N,Neff,S0,N0,Dn,D,Psi,p,Fm,Fp, &
       normal,coupling,i,x,y,t,Mm,Mp,initialize,dt,DE,DEm,DEp)

    use friction, only : fr_type,initial_state,solve_friction
    use material, only : block_material
    use fields, only : rotate_fields_xy2nt,rotate_fields_nt2xy, inplane_fault_mms
    use io, only : error

    implicit none

    type(fr_type),intent(in) :: FR
    real,intent(inout) :: V,O,S,N,Psi,Neff
    real,intent(in) :: S0,N0,Dn,D,p
    real,dimension(6),intent(inout) :: Fm,Fp
    real,dimension(2),intent(in) :: normal
    character(*),intent(in) :: coupling
    integer,intent(in) :: i
    real,intent(in) :: x,y,t,dt
    type(block_material),intent(in) :: Mm,Mp
    logical,intent(in) :: initialize
    real,intent(inout),optional :: DE,DEm,DEp

    logical :: previously_closed,compressive
    real :: Slock,Nlock,phis,phip,etas,etap, &
         vnp,vtp,snnp,sntp,sttp,szzp,vnFDp,vtFDp,snnFDp,sntFDp,sttFDp,szzFDp, &
         vnm,vtm,snnm,sntm,sttm,szzm,vnFDm,vtFDm,snnFDm,sntFDm,sttFDm,szzFDm

    ! 1. sntp+Zsp*vtp = sntFDp+Zsp*vtFDp 
    !    (S-wave into interface from plus  side)
    ! 2. sntm-Zsm*vtm = sntFDm-Zsm*vtFDm 
    !    (S-wave into interface from minus side)
    ! 3. sntp = sntm
    !    (continuity of shear stress)
    !
    ! 1-3 yield equation of form snt = phis-etas*V with
    ! phis = etas*(sntFDp/Zsp+sntFDm/Zsm+vtFDp-vtFDm),
    ! V = vtp-vtm, and 1/etas = 1/Zsp+1/Zsm
    !
    ! 4. snnp+Zpp*vnp = snnFDp+Zpp*vnFDp 
    !    (P-wave into interface from plus  side)
    ! 5. snnm-Zpm*vnm = snnFDm-Zpm*vnFDm 
    !    (P-wave into interface from minus side)
    ! 6. snnp = snnm
    !    (continuity of normal stress)
    !
    ! 1-3 yield equation of form snn = phip-Zp*O with
    ! phip = etap*(snnFDp/Zpp+snnFDm/Zpm+vnFDp-vnFDm),
    ! O = vnp-vnm, and 1/etap = 1/Zpp+1/Zpm

    ! rotate into local normal and tangential coordinates

    call rotate_fields_xy2nt(Fp,normal,vtFDp,vnFDp,sttFDp,sntFDp,snnFDp,szzFDp)
    call rotate_fields_xy2nt(Fm,normal,vtFDm,vnFDm,sttFDm,sntFDm,snnFDm,szzFDm)

    ! apply interface conditions in local coordinates

    etas = 1d0/(Mp%Zsi+Mm%Zsi)
    phis = etas*(sntFDp*Mp%Zsi+sntFDm*Mm%Zsi+vtFDp-vtFDm)

    etap = 1d0/(Mp%Zpi+Mm%Zpi)
    phip = etap*(snnFDp*Mp%Zpi+snnFDm*Mm%Zpi+vnFDp-vnFDm)

    ! stresses assuming no additional opening or slip
    
    Slock = S0+phis
    Nlock = N0-p-phip

    ! set V,O,S,N

    select case(coupling)

    case default

       call error('Invalid coupling (' // trim(coupling) &
            // ')','couple_mode2')

    case('locked') ! points are locked (no additional slip or opening)

       V = 0d0
       O = 0d0

       S = Slock
       N = Nlock

       snnp = -(N-N0+p)
       snnm = -(N-N0+p)

    case('friction')
       ! no friction law to solve in the case of RFSL, values come in correct!
       if(FR%friction_law .eq. 'RFSL') then
         if(initialize) then
           N = Nlock
           Neff = Nlock
           V = FR%Psi0
           S = Slock-etas*V
           O = 0d0
         end if
       else

         ! check for fault opening:
         ! Nlock = effective normal stress if fault is constrained against opening/closing
         ! (a) if previously closed (Dn=0) and Nlock>0 (compressive), then fault remains closed => friction
         ! (b) if previously closed (Dn=0) and Nlock<=0 (tensile), then fault opens => free surface
         ! (c) if previously opened (Dn>0) and Nlock<=0 (tensile), then fault remains open => free surface
         ! (d) if previously opened (Dn>0) and Nlock>0 (compressive), then fault is closing => free surface
         ! Note that with explicit time stepping methods, (d) may result in interpenetration

         !Dn = max(Dn,0d0) ! SHOULD ADD TO CORRECT INTERPENETRATION (BUT MAY BE MORE APPROPRIATE BEFORE OUTPUT)
         previously_closed = (Dn<=0d0) ! SHOULD BE STRICT EQUALITY, BUT NOT WITH INTERPENETRATION
         compressive = (Nlock>0d0)

         if (initialize) then
           O = 0d0
           N = Nlock
           V = vtFDp-vtFDm
           S = Slock-etas*V
           call initial_state(FR,Psi,V,S,N,i,x,y,dt)
         end if

         if (previously_closed) then
           if (compressive) then ! fault closed => friction
             O = 0d0
             N = Nlock
             Neff = Nlock + FR%subduction_B*phip
             if(FR%friction_law == 'RSL-mms') then
               N = inplane_fault_mms(x,y,t,1,'N')
             end if
             call solve_friction(FR,V,S,Neff,Slock,etas,D,Psi,i,x,y,t,.false.)
           else
             if (FR%opening) then ! fault opens => free surface
               N = 0d0
               O = -Nlock/etap
               S = 0d0
               V = Slock/etas
             else
               ! fault constrained against opening, so tensile normal stress,
               ! but offers no shear resistance to slip
               O = 0d0
               N = Nlock
               S = 0d0
               V = Slock/etas
             end if
           end if
         else ! fault open => free surface
           if (FR%opening) then ! fault opening permitted
             N = 0d0
             O = -Nlock/etap
             S = 0d0
             V = Slock/etas
           else
             ! fault constrained against opening, so tensile normal stress,
             ! but offers no shear resistance to slip
             O = 0d0
             N = Nlock
             S = 0d0
             V = Slock/etas
           end if
         end if
       end if

       snnp = -(N-N0+p)
       snnm = -(N-N0+p)

    case('eruption')
       
       ! fluid pressure balances normal stress

       Nlock = N0-phip ! instead of Nlock = N0-p-phip, as used for friction
       N = p
       O = (N-Nlock)/etap

       snnp = -(N-N0)
       snnm = -(N-N0)

       ! solve simultaneous equations for wall shear stress and interface-parallel velocities

       ! no shear resistance
       S = 0d0
       V = Slock/etas

       !tauFDp = sntFDp+Mp%Zs*vtFDp
       !tauFDm = sntFDm-Mm%Zs*vtFDm
       
       !call solve_wall_shear(ER,tauFDp,tauFDm,Mp%Zs,Mm%Zs,vtp,vtm,sntp,sntm,i)

    end select

    vnp = vnFDp+Mp%Zpi*(snnFDp-snnp)
    vnm = vnFDm-Mm%Zpi*(snnFDm-snnm)
    
    sttp = sttFDp-Mp%gamma*(snnFDp-snnp)
    sttm = sttFDm-Mm%gamma*(snnFDm-snnm)
    
    szzp = szzFDp-Mp%gamma*(snnFDp-snnp)
    szzm = szzFDm-Mm%gamma*(snnFDm-snnm)
       
    sntp = S-S0 ! not valid for eruption
    sntm = S-S0 ! not valid for eruption
    vtp = vtFDp+Mp%Zsi*(sntFDp-sntp)
    vtm = vtFDm-Mm%Zsi*(sntFDm-sntm)
    
    ! energy flux (positive when energy flows out of medium into fault)

    if (present(DE)) then
       DEm = DEm+sntm*vtm-0.25d0*Mm%Zsi*((sntFDm-sntm)+Mm%Zs*(vtFDm-vtm))**2
       DEp = DEp-sntp*vtp-0.25d0*Mp%Zsi*((sntFDp-sntp)-Mp%Zs*(vtFDp-vtp))**2
       DEm = DEm+snnm*vnm-0.25d0*Mm%Zpi*((snnFDm-snnm)+Mm%Zs*(vnFDm-vnm))**2
       DEp = DEp-snnp*vnp-0.25d0*Mp%Zpi*((snnFDp-snnp)-Mp%Zs*(vnFDp-vnp))**2
       ! DE = DE + DEp + DEm
       DE = DE-S*V
       ! DE = DE-S*V- &
       !      0.25d0*Mm%Zsi*((snzFDm-snzm)+Mm%Zs*(vzFDm-vzm))**2- &
       !      0.25d0*Mp%Zsi*((snzFDp-snzp)-Mp%Zs*(vzFDp-vzp))**2
    end if

    ! rotate back to x-y coordinates

    call rotate_fields_nt2xy(Fp,normal,vtp,vnp,sttp,sntp,snnp,szzp)
    call rotate_fields_nt2xy(Fm,normal,vtm,vnm,sttm,sntm,snnm,szzm)

  end subroutine couple_mode2


  subroutine cycle_stress_fault(FR)
    use friction, only : fr_type
    implicit none
    type(fr_type),intent(inout) :: FR
    
    if(allocated(FR%V)) FR%V = 0d0

  end subroutine cycle_stress_fault

end module fault
