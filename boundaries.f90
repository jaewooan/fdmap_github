module boundaries

  implicit none

  type :: block_boundaries
     character(24) :: bcL,bcR,bcB,bcT
     real :: rhog ! for tsunami calculations, not ideal place to store this, but easy
     character(256) :: problem ! might want to change variable name
  end type block_boundaries


contains


  subroutine init_boundaries(iblock,B,input,echo)

    use mpi_routines, only : is_master
    use io, only : error,write_matlab,seek_to_string

    implicit none
    
    integer,intent(in) :: iblock,input,echo
    type(block_boundaries),intent(out) :: B

    integer :: stat
    character(24) :: bcL,bcR,bcB,bcT
    real :: rhog
    character(256) :: str
    character(256) :: Bstr
    character(256) :: problem

    namelist /boundaries_list/ bcL,bcR,bcB,bcT,rhog,problem

    ! defaults

    bcL = 'none'
    bcR = 'none'
    bcB = 'none'
    bcT = 'none'
    rhog = 0d0
    problem = ''

    ! read in boundary condition parameters
    
    write(str,'(a,i0,a)') '!---BLOCK',iblock,'---'
    call seek_to_string(input,str)
    read(input,nml=boundaries_list,iostat=stat)
    if (stat>0) call error(trim(str) // ' :: Error in boundaries_list','init')

    ! store

    B%bcL = bcL
    B%bcR = bcR
    B%bcB = bcB
    B%bcT = bcT
    B%rhog = rhog
    B%problem = problem

    ! output boundary condition parameters

    if (is_master) then
       write(Bstr,'(a,i0,a)') 'B{',iblock,'}'
       call write_matlab(echo,'bcL',B%bcL,Bstr)
       call write_matlab(echo,'bcR',B%bcR,Bstr)
       call write_matlab(echo,'bcB',B%bcB,Bstr)
       call write_matlab(echo,'bcT',B%bcT,Bstr)
       if (B%rhog/=0d0) call write_matlab(echo,'rhog',B%rhog,Bstr)
       if (B%problem/='') call write_matlab(echo,'problem',B%problem,Bstr)
    end if

  end subroutine init_boundaries


  subroutine enforce_boundary_conditions(G,F,B,mode,t,iblock)

    use grid, only : block_grid
    use fields, only : block_fields

    implicit none

    type(block_grid),intent(in) :: G
    type(block_fields),intent(inout) :: F
    type(block_boundaries),intent(in) :: B
    integer,intent(in) :: mode,iblock
    real,intent(in) :: t

    if (G%skip) return ! process has no points in this block

    ! q-direction

    if (G%nx/=1) then
       if (G%sideL) call apply_bc_side(G%bndL,F%bndFL,G%my,G%py,B%bcL,B%rhog,mode,t,iblock,B%problem)
       if (G%sideR) call apply_bc_side(G%bndR,F%bndFR,G%my,G%py,B%bcR,B%rhog,mode,t,iblock,B%problem)
    end if

    ! r-direction

    if (G%ny/=1) then
       if (G%sideB) call apply_bc_side(G%bndB,F%bndFB,G%mx,G%px,B%bcB,B%rhog,mode,t,iblock,B%problem)
       if (G%sideT) call apply_bc_side(G%bndT,F%bndFT,G%mx,G%px,B%bcT,B%rhog,mode,t,iblock,B%problem)
    end if

  end subroutine enforce_boundary_conditions


  subroutine apply_bc_side(bndC,bndF,m,p,bc,rhog,mode,t,iblock,problem)
    
    use geometry, only : curve
    use fields, only : bnd_fields

    implicit none
    
    type(curve),intent(in) :: bndC
    type(bnd_fields),intent(inout) :: bndF
    integer,intent(in) :: m,p
    character(*),intent(in) :: bc,problem
    real,intent(in) :: rhog,t
    integer,intent(in) :: mode,iblock

    integer :: i

    if (bc=='none') return

    do i = m,p
       call set_bc(bndF%Fhat(i,:),bndF%F(i,:),bndF%F0(i,:),bndF%U(i,:),rhog, &
            bndC%x(i),bndC%y(i),bndC%n(i,:),bc,bndF%M(i,1:3),mode,t,iblock,problem)
    end do
    
  end subroutine apply_bc_side


  subroutine set_bc(Fhat,F,F0,U,rhog,x,y,normal,bc,M,mode,t,iblock,problem)

    implicit none

    ! F  = fields (grid data)
    ! FHAT = fields (hat variables)
    ! F0 = initial fields

    real,intent(out) :: Fhat(:)
    real,intent(in) :: F(:),F0(:),U(:),rhog,x,y,normal(2),M(3),t
    character(*),intent(in) :: bc,problem
    integer,intent(in) :: mode,iblock

    select case(mode)
    case(2)
       call set_bc_mode2(Fhat,F,F0,U,rhog,normal,bc,M(1),M(2),M(3),x,y,t,iblock,problem)
    case(3)
       call set_bc_mode3(Fhat,F,F0,normal,bc,M(1),x,y,t,iblock)
    end select

  end subroutine set_bc


  subroutine set_bc_mode3(Fhat,F,F0,normal,bc,Zs,x,y,t,iblock)

    use fields, only : rotate_fields_xy2nt,rotate_fields_nt2xy
    use io, only : error

    implicit none

    real,intent(out) :: Fhat(3)
    real,intent(in) :: F(3),F0(3),normal(2),Zs,x,y,t
    character(*),intent(in) :: bc
    integer,intent(in) :: iblock

    real :: Zsi,vz,snz,vzFD,snzFD,stzFD,vzEX,snzEX,stzEX,Fe(3)

    if (Zs==0d0) then
       ! fluid, no boundary condition
       Fhat = F
       return
    end if
       
    Zsi = 1d0/Zs

    ! subtract initial fields (if needed for BC) and
    ! rotate into local normal and tangential coordinates
    
    select case(bc)
    case default
       call rotate_fields_xy2nt(F-F0,normal,vzFD,stzFD,snzFD)
    end select

    ! then enforce BC

    select case(bc)
    case('absorbing')
       vz   = 0.5d0*( vzFD-Zsi*snzFD)
       snz  = 0.5d0*(snzFD-Zs * vzFD)
    case('free')
       vz   = vzFD-Zsi*snzFD
       snz  = 0d0
    case('rigid')
       vz   = 0d0
       snz  = snzFD-Zs*vzFD
    case default
       call error('Invalid boundary condition (' // &
            trim(bc) // ')','set_bc_mode3')
    end select

    ! rotate back to x-y coordinates

    call rotate_fields_nt2xy(Fhat,normal,vz,stzFD,snz)

    ! and add initial fields (if needed)

    select case(bc)
    case default
       Fhat = Fhat+F0
    end select

  end subroutine set_bc_mode3


  subroutine set_bc_mode2(Fhat,F,F0,U,rhog,normal,bc,Zs,Zp,gamma,x,y,t,iblock,problem)

    use fields, only : rotate_fields_xy2nt,rotate_fields_nt2xy
    use mms, only : mms_sin,inplane_fault_mms,mms_hydrofrac
    use tsunami, only : seafloor_velocity
    use geometry, only : rotate_xy2nt
    use io, only : error

    implicit none

    real,intent(out) :: Fhat(6)
    real,intent(in) :: F(6),F0(6),U(2),rhog,normal(2),Zs,Zp,gamma,x,y,t
    character(*),intent(in) :: bc,problem
    integer,intent(in) :: iblock

    real :: Zsi,Zpi,vn,vt,snn,snt,stt,szz,vnFD,vtFD,snnFD,sntFD,sttFD,szzFD
    real :: vnEX,vtEX,snnEX,sntEX,sttEX,szzEX,FEX(6),Ut,Un,vx,vy

    if (Zs==0d0) then
       Zsi = 0d0 ! avoid division by zero for fluid case
    else
       Zsi = 1d0/Zs
    end if

    Zpi = 1d0/Zp

    ! subtract initial fields (if needed for BC) and
    ! rotate into local normal and tangential coordinates

    select case(bc)
    case default
       call rotate_fields_xy2nt(F-F0,normal,vtFD,vnFD,sttFD,sntFD,snnFD,szzFD)
    case('mms-sin-w','mms-hydrofrac-w','inplane-fault-mms-w','rigid-0','absorbing-0','free-0')
       call rotate_fields_xy2nt(F  ,normal,vtFD,vnFD,sttFD,sntFD,snnFD,szzFD)
    end select

    ! then enforce BC

    select case(bc)
       
    case('mms-sin-w','inplane-fault-mms-w','mms-hydrofrac-w') ! MMS
       
       select case(bc)
       case('inplane-fault-mms-w')
          FEX(1) = inplane_fault_mms(x,y,t,iblock,'vx')
          FEX(2) = inplane_fault_mms(x,y,t,iblock,'vy')
          FEX(3) = inplane_fault_mms(x,y,t,iblock,'sxx')
          FEX(4) = inplane_fault_mms(x,y,t,iblock,'sxy')
          FEX(5) = inplane_fault_mms(x,y,t,iblock,'syy')
          FEX(6) = inplane_fault_mms(x,y,t,iblock,'szz')
       case('mms-sin-w')
          FEX(1) = mms_sin(x,y,t,iblock,'vx')
          FEX(2) = mms_sin(x,y,t,iblock,'vy')
          FEX(3) = mms_sin(x,y,t,iblock,'sxx')
          FEX(4) = mms_sin(x,y,t,iblock,'sxy')
          FEX(5) = mms_sin(x,y,t,iblock,'syy')
          FEX(6) = 0d0
       case('mms-hydrofrac-w')
          FEX(1) = mms_hydrofrac(x,y,t,iblock,'vx')
          FEX(2) = mms_hydrofrac(x,y,t,iblock,'vy')
          FEX(3) = mms_hydrofrac(x,y,t,iblock,'sxx')
          FEX(4) = mms_hydrofrac(x,y,t,iblock,'sxy')
          FEX(5) = mms_hydrofrac(x,y,t,iblock,'syy')
          FEX(6) = 0d0
       end select
       call rotate_fields_xy2nt(FEX,normal,vtEX,vnEX,sttEX,sntEX,snnEX,szzEX)
       vt  = 0.5d0*(vtEX+vtFD+Zsi*(sntEX-sntFD))
       snt = 0.5d0*(sntEX+sntFD+Zs*(vtEX-vtFD))
       vn  = 0.5d0*(vnEX+vnFD+Zpi*(snnEX-snnFD))
       snn = 0.5d0*(snnEX+snnFD+Zp*(vnEX-vnFD))
       stt = sttFD-gamma*(snnFD-snn)
       szz = szzFD-gamma*(snnFD-snn)
       
    case('rigid-0','absorbing-0','free-0') 
       
       ! apply boundary conditions in local coordinates
       
       select case(bc)
       case('absorbing-0')
          vt  = 0.5d0*( vtFD-Zsi*sntFD)
          snt = 0.5d0*(sntFD-Zs * vtFD)
          vn  = 0.5d0*( vnFD-Zpi*snnFD)
          snn = 0.5d0*(snnFD-Zp * vnFD)
       case('free-0')
          vt  = vtFD-Zsi*sntFD
          snt = 0d0
          vn  = vnFD-Zpi*snnFD
          snn = 0d0
       case('rigid-0')
          vt  = 0d0
          snt = sntFD-Zs*vtFD
          vn  = 0d0
          snn = snnFD-Zp*vnFD
       end select
       stt = sttFD-gamma*(snnFD-snn)
       szz = szzFD-gamma*(snnFD-snn)

    case('rigid','absorbing','free','absorbing-velocity','free-velocity','tsunami','seafloor') 
       
       ! apply boundary conditions in local coordinates
       
       select case(bc)
       case('absorbing')
          vt  = 0.5d0*( vtFD-Zsi*sntFD)
          snt = 0.5d0*(sntFD-Zs * vtFD)
          vn  = 0.5d0*( vnFD-Zpi*snnFD)
          snn = 0.5d0*(snnFD-Zp * vnFD)
       case('free')
          vt  = vtFD-Zsi*sntFD
          snt = 0d0
          vn  = vnFD-Zpi*snnFD
          snn = 0d0
       case('rigid')
          vt  = 0d0
          snt = sntFD-Zs*vtFD
          vn  = 0d0
          snn = snnFD-Zp*vnFD
       case('absorbing-velocity')
          vt  = -Zsi*sntFD
          snt = sntFD
          vn  = -Zpi*snnFD
          snn = snnFD
       case('free-velocity')
          vt  = vtFD-Zsi*sntFD
          snt = sntFD
          vn  = vnFD-Zpi*snnFD
          snn = snnFD
       case('tsunami')
          vt  = vtFD-Zsi*sntFD
          snt = 0d0
          ! rotate displacement vector
          call rotate_xy2nt(U(1),U(2),Ut,Un,normal)
          snn = -rhog*Un
          vn  = vnFD-Zpi*(snnFD-snn)
       case('seafloor')
          call seafloor_velocity(x,y,t,vx,vy,problem)
          call rotate_xy2nt(vx,vy,vt,vn,normal)
          snt = sntFD-Zs*(vtFD-vt)
          snn = snnFD-Zp*(vnFD-vn)
       end select
       stt = sttFD-gamma*(snnFD-snn)
       szz = szzFD-gamma*(snnFD-snn)
       
    case default
       call error('Invalid boundary condition (' // &
            trim(bc) // ')','set_bc_mode2')
    end select

    ! rotate back to x-y coordinates
    
    call rotate_fields_nt2xy(Fhat,normal,vt,vn,stt,snt,snn,szz)
    
    ! and add initial fields (if needed)

    select case(bc)
    case default
       Fhat = Fhat+F0
    case('mms-sin-w','mms-hydrofrac-w','inplane-fault-mms-w','rigid-0','absorbing-0','free-0')
       ! initial fields not required 
    end select
    
  end subroutine set_bc_mode2


end module boundaries
