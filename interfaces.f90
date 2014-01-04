module interfaces

  implicit none


contains


  subroutine enforce_interface_conditions(I,Fm,Fp,C,mode,t)

    use boundaries, only : iface_type
    use fields, only : block_fields
    use mpi_routines2d, only : cartesian
    use friction, only : set_rates_friction

    implicit none

    type(iface_type),intent(inout) :: I
    type(block_fields),intent(inout) :: Fm,Fp
    type(cartesian),intent(in) :: C
    integer,intent(in) :: mode
    real,intent(in) :: t

    if (I%skip) return ! process has no grid points adjacent to interface

    ! set hat variables

    select case(I%coupling)

    case('locked')

       select case(I%direction)
       case('x')
          call enforce_locked_interface(I%m,I%p, &
               Fm%bndFR%Fhat(I%m:I%p, : ),Fp%bndFL%Fhat(I%m:I%p, : ), &
               Fm%bndFR%F   (I%m:I%p, : ),Fp%bndFL%F   (I%m:I%p, : ), &
               Fm%bndFR%M   (I%m:I%p,1:3),Fp%bndFL%M   (I%m:I%p,1:3), &
               I%nhat(I%m:I%p,:),mode)
       case('y')
          call enforce_locked_interface(I%m,I%p, &
               Fm%bndFT%Fhat(I%m:I%p, : ),Fp%bndFB%Fhat(I%m:I%p, : ), &
               Fm%bndFT%F   (I%m:I%p, : ),Fp%bndFB%F   (I%m:I%p, : ), &
               Fm%bndFT%M   (I%m:I%p,1:3),Fp%bndFB%M   (I%m:I%p,1:3), &
               I%nhat(I%m:I%p,:),mode)
       end select

    case('friction')

       select case(I%direction)
       case('x')
          call enforce_frictional_interface(I%m,I%p, &
               Fm%bndFR%Fhat(I%m:I%p, : ),Fp%bndFL%Fhat(I%m:I%p, : ), &
               Fm%bndFR%F   (I%m:I%p, : ),Fp%bndFL%F   (I%m:I%p, : ), &
               Fm%bndFR%M   (I%m:I%p,1:3),Fp%bndFL%M   (I%m:I%p,1:3), &
               I%nhat(I%m:I%p,:),mode, &
               I%x(I%m:I%p),I%y(I%m:I%p),t,I%FR,I%TP,Fp%F0)
       case('y')
          call enforce_frictional_interface(I%m,I%p, &
               Fm%bndFT%Fhat(I%m:I%p, : ),Fp%bndFB%Fhat(I%m:I%p, : ), &
               Fm%bndFT%F   (I%m:I%p, : ),Fp%bndFB%F   (I%m:I%p, : ), &
               Fm%bndFT%M   (I%m:I%p,1:3),Fp%bndFB%M   (I%m:I%p,1:3), &
               I%nhat(I%m:I%p,:),mode, &
               I%x(I%m:I%p),I%y(I%m:I%p),t,I%FR,I%TP,Fp%F0)
       end select

       ! set rates for auxiliary fields (like slip and state variables),
       ! also calculate rupture time
          
       call set_rates_friction(I%FR,I%m,I%p,I%x(I%m:I%p),I%y(I%m:I%p),t)

    end select

  end subroutine enforce_interface_conditions


  subroutine enforce_locked_interface(m,p,Fhatm,Fhatp,Fm,Fp,Mm,Mp,nhat,mode)

    implicit none

    integer,intent(in) :: m,p,mode
    real,dimension(m:,:),intent(out) :: Fhatm,Fhatp
    real,dimension(m:,:),intent(in) :: Fm,Fp,Mm,Mp,nhat

    integer :: i

    do i = m,p
       select case(mode)
       case(2)
          call locked_interface_mode2(Fhatm(i,:),Fhatp(i,:),Fm(i,:),Fp(i,:),nhat(i,:), &
               Mm(i,1),Mp(i,1),Mm(i,2),Mp(i,2),Mm(i,3),Mp(i,3))
       case(3)
          call locked_interface_mode3(Fhatm(i,:),Fhatp(i,:),Fm(i,:),Fp(i,:),nhat(i,:), &
               Mm(i,1),Mp(i,1))
       end select
    end do
    
  end subroutine enforce_locked_interface


  subroutine locked_interface_mode3(Fhatm,Fhatp,Fm,Fp,normal,Zsm,Zsp)

    use fields, only : rotate_fields_xy2nt,rotate_fields_nt2xy

    implicit none

    real,intent(out) :: Fhatm(3),Fhatp(3)
    real,intent(in) :: Fm(3),Fp(3),normal(2),Zsm,Zsp

    real :: Zsim,Zsip,phis,etas, &
         vzp,snzp,vzFDp,snzFDp,stzFDp, &
         vzm,snzm,vzFDm,snzFDm,stzFDm

    ! rotate into local normal and tangential coordinates

    call rotate_fields_xy2nt(Fp,normal,vzFDp,stzFDp,snzFDp)
    call rotate_fields_xy2nt(Fm,normal,vzFDm,stzFDm,snzFDm)

    ! enforce interface conditions

    Zsip = 1d0/Zsp
    Zsim = 1d0/Zsm

    etas = 1d0/(Zsip+Zsim)
    phis = etas*(snzFDp*Zsip+snzFDm*Zsim+vzFDp-vzFDm)

    snzp = phis
    snzm = phis
    vzp = vzFDp+Zsip*(snzFDp-snzp)
    vzm = vzFDm-Zsim*(snzFDm-snzm)

    ! rotate back to x-y coordinates

    call rotate_fields_nt2xy(Fhatp,normal,vzp,stzFDp,snzp)
    call rotate_fields_nt2xy(Fhatm,normal,vzm,stzFDm,snzm)
    
  end subroutine locked_interface_mode3


  subroutine locked_interface_mode2(Fhatm,Fhatp,Fm,Fp,normal,Zsm,Zsp,Zpm,Zpp,gammam,gammap)

    use fields, only : rotate_fields_xy2nt,rotate_fields_nt2xy

    implicit none

    real,intent(out) :: Fhatm(6),Fhatp(6)
    real,intent(in) :: Fm(6),Fp(6),normal(2),Zsm,Zsp,Zpm,Zpp,gammam,gammap

    real :: Zsim,Zsip,Zpim,Zpip,phis,phip,etas,etap, &
         vnp,vtp,snnp,sntp,sttp,szzp,vnFDp,vtFDp,snnFDp,sntFDp,sttFDp,szzFDp, &
         vnm,vtm,snnm,sntm,sttm,szzm,vnFDm,vtFDm,snnFDm,sntFDm,sttFDm,szzFDm

    ! rotate into local normal and tangential coordinates

    call rotate_fields_xy2nt(Fp,normal,vtFDp,vnFDp,sttFDp,sntFDp,snnFDp,szzFDp)
    call rotate_fields_xy2nt(Fm,normal,vtFDm,vnFDm,sttFDm,sntFDm,snnFDm,szzFDm)

    ! enforce interface conditions

    Zsip = 1d0/Zsp
    Zsim = 1d0/Zsm
    Zpip = 1d0/Zpp
    Zpim = 1d0/Zpm

    ! normal components

    etap = 1d0/(Zpip+Zpim)
    phip = etap*(snnFDp*Zpip+snnFDm*Zpim+vnFDp-vnFDm)

    snnp = phip
    snnm = phip

    vnp = vnFDp+Zpip*(snnFDp-snnp)
    vnm = vnFDm-Zpim*(snnFDm-snnm)
    
    sttp = sttFDp-gammap*(snnFDp-snnp)
    sttm = sttFDm-gammam*(snnFDm-snnm)
    
    szzp = szzFDp-gammap*(snnFDp-snnp)
    szzm = szzFDm-gammam*(snnFDm-snnm)

    ! shear components

    etas = 1d0/(Zsip+Zsim)
    phis = etas*(sntFDp*Zsip+sntFDm*Zsim+vtFDp-vtFDm)
       
    sntp = phis
    sntm = phis
    vtp = vtFDp+Zsip*(sntFDp-sntp)
    vtm = vtFDm-Zsim*(sntFDm-sntm)
    
    ! rotate back to x-y coordinates

    call rotate_fields_nt2xy(Fhatp,normal,vtp,vnp,sttp,sntp,snnp,szzp)
    call rotate_fields_nt2xy(Fhatm,normal,vtm,vnm,sttm,sntm,snnm,szzm)
    
  end subroutine locked_interface_mode2


  subroutine enforce_frictional_interface(m,p,Fhatm,Fhatp,Fm,Fp,Mm,Mp,nhat,mode,x,y,t,FR,TP,F0)

    use friction, only : fr_type
    use thermpres, only : tp_type

    implicit none

    integer,intent(in) :: m,p,mode
    real,dimension(m:,:),intent(out) :: Fhatm,Fhatp
    real,dimension(m:,:),intent(in) :: Fm,Fp,Mm,Mp,nhat
    real,dimension(m:),intent(in) :: x,y
    real,intent(in) :: t,F0(9)
    type(fr_type),intent(inout) :: FR
    type(tp_type),intent(in) :: TP

    integer :: i

    do i = m,p
       select case(mode)
       case(2)
          call frictional_interface_mode2(Fhatm(i,:),Fhatp(i,:),Fm(i,:),Fp(i,:),nhat(i,:), &
               Mm(i,1),Mp(i,1),Mm(i,2),Mp(i,2),Mm(i,3),Mp(i,3),x(i),y(i),t,i,FR,TP,F0)
       case(3)
          call frictional_interface_mode3(Fhatm(i,:),Fhatp(i,:),Fm(i,:),Fp(i,:),nhat(i,:), &
               Mm(i,1),Mp(i,1),x(i),y(i),t,i,FR,TP,F0)
       end select
    end do
    
  end subroutine enforce_frictional_interface


  subroutine frictional_interface_mode3(Fhatm,Fhatp,Fm,Fp,normal,Zsm,Zsp,x,y,t,i,FR,TP,F0)

    use geometry, only : rotate_xy2nt
    use fields, only : rotate_fields_xy2nt,rotate_fields_nt2xy
    use friction, only : fr_type,solve_friction
    use thermpres, only : tp_type,pressure_thermpres

    implicit none

    real,intent(out) :: Fhatm(3),Fhatp(3)
    real,intent(in) :: Fm(3),Fp(3),normal(2),Zsm,Zsp,x,y,t,F0(9)
    integer,intent(in) :: i
    type(fr_type),intent(inout) :: FR
    type(tp_type),intent(in) :: TP

    real :: Zsim,Zsip,phis,etas, &
         vzp,snzp,vzFDp,snzFDp,stzFDp, &
         vzm,snzm,vzFDm,snzFDm,stzFDm
    real :: stt,snt,snn,p,phip

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

    ! balance tractions and define stress transfer term phis

    Zsip = 1d0/Zsp
    Zsim = 1d0/Zsm

    etas = 1d0/(Zsip+Zsim)
    phis = etas*(snzFDp*Zsip+snzFDm*Zsim+vzFDp-vzFDm)

    ! normal stress contribution from resolving in-plane stress field onto fault

    call rotate_xy2nt(F0(4),F0(5),F0(7),stt,snt,snn,normal)

    ! reduction of effective normal stress from pore pressure

    p = 0d0
    call pressure_thermpres(TP,i,p)
    phip = snn+p

    ! solve for V,S,O,N by enforcing friction law and no opening
        
    call solve_friction(FR,FR%V(i),FR%S(i),FR%O(i),FR%N(i),phip,phis,etas,FR%D(i),FR%Psi(i),i,x,y,t,.false.)

    snzp = phis
    snzm = phis
    vzp = vzFDp+Zsip*(snzFDp-snzp)
    vzm = vzFDm-Zsim*(snzFDm-snzm)

    ! rotate back to x-y coordinates

    call rotate_fields_nt2xy(Fhatp,normal,vzp,stzFDp,snzp)
    call rotate_fields_nt2xy(Fhatm,normal,vzm,stzFDm,snzm)
    
  end subroutine frictional_interface_mode3


  subroutine frictional_interface_mode2(Fhatm,Fhatp,Fm,Fp,normal,Zsm,Zsp,Zpm,Zpp,gammam,gammap,x,y,t,i,FR,TP,F0)

    use fields, only : rotate_fields_xy2nt,rotate_fields_nt2xy
    use friction, only : fr_type,solve_friction
    use thermpres, only : tp_type,pressure_thermpres

    implicit none

    real,intent(out) :: Fhatm(6),Fhatp(6)
    real,intent(in) :: Fm(6),Fp(6),normal(2),Zsm,Zsp,Zpm,Zpp,gammam,gammap,x,y,t,F0(9)
    integer,intent(in) :: i
    type(fr_type),intent(inout) :: FR
    type(tp_type),intent(in) :: TP

    real :: Zsim,Zsip,Zpim,Zpip,phis,phip,etas,etap, &
         vnp,vtp,snnp,sntp,sttp,szzp,vnFDp,vtFDp,snnFDp,sntFDp,sttFDp,szzFDp, &
         vnm,vtm,snnm,sntm,sttm,szzm,vnFDm,vtFDm,snnFDm,sntFDm,sttFDm,szzFDm
    real :: p

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

    ! balance tractions and define stress transfer terms phis and phip

    Zsip = 1d0/Zsp
    Zsim = 1d0/Zsm
    Zpip = 1d0/Zpp
    Zpim = 1d0/Zpm

    etas = 1d0/(Zsip+Zsim)
    phis = etas*(sntFDp*Zsip+sntFDm*Zsim+vtFDp-vtFDm)

    etap = 1d0/(Zpip+Zpim)
    phip = etap*(snnFDp*Zpip+snnFDm*Zpim+vnFDp-vnFDm)

    ! reduction of effective normal stress from pore pressure

    p = 0d0
    call pressure_thermpres(TP,i,p)
    phip = phip+p

    ! solve for V,S,O,N by enforcing friction law and no opening

    call solve_friction(FR,FR%V(i),FR%S(i),FR%O(i),FR%N(i),phip,phis,etas,FR%D(i),FR%Psi(i),i,x,y,t,.false.)

    snnp = -(FR%N(i)-FR%N0(i))
    snnm = -(FR%N(i)-FR%N0(i))

    vnp = vnFDp+Zpip*(snnFDp-snnp)
    vnm = vnFDm-Zpim*(snnFDm-snnm)
    
    sttp = sttFDp-gammap*(snnFDp-snnp)
    sttm = sttFDm-gammam*(snnFDm-snnm)
    
    szzp = szzFDp-gammap*(snnFDp-snnp)
    szzm = szzFDm-gammam*(snnFDm-snnm)
    
    sntp = FR%S(i)-FR%S0(i)
    sntm = FR%S(i)-FR%S0(i)

    vtp = vtFDp+Zsip*(sntFDp-sntp)
    vtm = vtFDm-Zsim*(sntFDm-sntm)

    ! rotate back to x-y coordinates

    call rotate_fields_nt2xy(Fhatp,normal,vtp,vnp,sttp,sntp,snnp,szzp)
    call rotate_fields_nt2xy(Fhatm,normal,vtm,vnm,sttm,sntm,snnm,szzm)
    
  end subroutine frictional_interface_mode2


end module interfaces
