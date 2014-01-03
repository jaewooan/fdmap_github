module boundaries

  use friction, only : fr_type
  use thermpres, only : tp_type
  use erupt, only : er_type

  implicit none

  type :: block_boundaries
     character(24) :: bcL,bcR,bcB,bcT
     real :: rhog ! for tsunami calculations, not ideal place to store this, but easy
  end type block_boundaries

  type :: iface_type
     logical :: skip,share,process_m,process_p,output_process
     character(24) :: coupling
     character(1) :: direction
     integer :: iblockm,iblockp,m,p,mg,pg,mb,pb, &
          rank_m,rank_p,comm_m,comm_p,comm_mp,comm,array_w,array_s
     type(fr_type) :: FR
     type(tp_type) :: TP
     type(er_type) :: ER
     real,dimension(:),allocatable :: x,y,dl
     real,dimension(:,:),allocatable :: nhat
  end type iface_type


contains


  subroutine init_iface(iface,I,input)

    use io, only : error,seek_to_string

    implicit none

    integer,intent(in) :: iface,input
    type(iface_type),intent(out) :: I

    integer :: stat,iblockm,iblockp,mg,pg
    character(24) :: coupling
    character(1) :: direction
    character(256) :: str
    
    namelist /interface_list/ coupling,iblockm,iblockp,direction,mg,pg
    
    ! defaults

    coupling = 'locked'

    ! read in interface parameters
        
    write(str,'(a,i0,a)') '!---IFACE',iface,'---'
    call seek_to_string(input,str)
    read(input,nml=interface_list,iostat=stat)
    if (stat>0) call error('Error in interface_list','init_iface')

    I%coupling = coupling
    I%iblockm = iblockm
    I%iblockp = iblockp
    I%direction = direction
    I%mg = mg
    I%pg = pg

  end subroutine init_iface


  subroutine init_iface_blocks(iface,I,Gm,Gp,C,refine)

    use grid, only : block_grid
    use mpi_routines2d, only : cartesian
    use mpi_routines, only : MPI_REAL_PW,MPI_REAL_PS,is_master,subarray,new_communicator
    use io, only : message,error

    implicit none

    integer,intent(in) :: iface
    type(iface_type),intent(inout) :: I
    type(block_grid),intent(in) :: Gm,Gp
    type(cartesian),intent(in) :: C
    real,intent(in) :: refine

    character(74) :: str
    integer :: mg,pg,m,p,iblock

    ! global indices prior to refinement

    mg = I%mg
    pg = I%pg

    ! adjust indices for refinement

    select case(I%direction)
    case default
       call error('Invalid direction','init_iface_blocks')
    case('x')
       iblock = Gm%iblock_y
    case('y')
       iblock = Gm%iblock_x
    end select

    if (mg==pg) then
       I%mg = ceiling(dble(mg-iblock)*refine)+iblock
       I%pg = I%mg
    else
       I%mg = ceiling(dble(mg-iblock)*refine)+iblock
       I%pg = ceiling(dble(pg-iblock)*refine)+iblock
    end if

    if (is_master) then
       write(str,'(a,i0,2a,2(a,i2),2(a,i6))') 'iface',iface, &
            ':  dir=',I%direction, &
            ',  block_m=',I%iblockm,',  block_p=',I%iblockp, &
            ',  mg=',I%mg,',  pg=',I%pg
       call message(str)
    end if

    ! determine if process has cells adjacent to interface (skip=.false.);
    ! if so, and if blocks on each side of interface are held by 
    ! different processes, share fields between processes (share=.true.)

    select case(I%direction)
    case default
       call error('Invalid direction','init_iface_blocks')
    case('x')
       m = max(Gm%my,Gp%my)
       p = min(Gm%py,Gp%py)
       I%skip = (.not.(Gm%sideR.or.Gp%sideL)).or.I%pg<m.or.p<I%mg
       I%share = (Gm%sideR.neqv.Gp%sideL).and.(.not.I%skip)
       if (I%share.and.Gm%sideR) then
          I%rank_m = C%c1dx%myid
          I%rank_p = C%c1dx%rank_p
       end if
       if (I%share.and.Gp%sideL) then
          I%rank_p = C%c1dx%myid
          I%rank_m = C%c1dx%rank_m
       end if
       I%process_m = Gm%sideR.and.(.not.I%skip)
       I%process_p = Gp%sideL.and.(.not.I%skip)
    case('y')
       m = max(Gm%mx,Gp%mx)
       p = min(Gm%px,Gp%px)
       I%skip = (.not.(Gm%sideT.or.Gp%sideB)).or.I%pg<m.or.p<I%mg
       I%share = (Gm%sideT.neqv.Gp%sideB).and.(.not.I%skip)
       if (I%share.and.Gm%sideT) then
          I%rank_m = C%c1dy%myid
          I%rank_p = C%c1dy%rank_p
       end if
       if (I%share.and.Gp%sideB) then
          I%rank_p = C%c1dy%myid
          I%rank_m = C%c1dy%rank_m
       end if
       I%process_m = Gm%sideT.and.(.not.I%skip)
       I%process_p = Gp%sideB.and.(.not.I%skip)
    end select

    I%output_process = I%process_m ! minus-side processes used for output

    ! limits of interface between blocks for this process

    select case(I%direction)
    case('x')
       I%m = max(I%mg,m)
       I%p = min(I%pg,p)
    case('y')
       I%m = max(I%mg,m)
       I%p = min(I%pg,p)
    end select

    ! create distributed array data type for I/O

    if (I%process_m.or.I%process_p.or.I%output_process) then
       call subarray(I%pg-I%mg+1,I%m-I%mg+1,I%p-I%mg+1,MPI_REAL_PW,I%array_w)
       call subarray(I%pg-I%mg+1,I%m-I%mg+1,I%p-I%mg+1,MPI_REAL_PS,I%array_s)
    end if

    ! communicators
    
    call new_communicator(I%process_m,I%comm_m)
    call new_communicator(I%process_p,I%comm_p)
    call new_communicator(.not.I%skip,I%comm_mp)
    call new_communicator(I%output_process,I%comm)

    ! return if needed (except master, which must 
    ! continue for proper output to echo file)

    if (I%skip.and.(.not.is_master)) return

    if (.not.I%skip) then

       ! allocate memory for interface fields 
       
       allocate(I%x   (I%m:I%p),I%y   (I%m:I%p))
       allocate(I%dl  (I%m:I%p),I%nhat(I%m:I%p,2))

       I%x    = 1d40
       I%y    = 1d40
       I%nhat = 1d40
       I%dl   = 1d40

       ! initialize coordinates and unit normal
       
       select case(I%direction)
       case('x')
          if (.not.Gm%skip) then
             I%x    =  Gm%bndR%x(I%m:I%p)
             I%y    =  Gm%bndR%y(I%m:I%p)
             I%nhat =  Gm%bndR%n(I%m:I%p,:)
             I%dl   =  sqrt(Gm%xrR(I%m:I%p)**2+Gm%yrR(I%m:I%p)**2)
          else
             I%x    =  Gp%bndL%x(I%m:I%p)
             I%y    =  Gp%bndL%y(I%m:I%p)
             I%nhat = -Gp%bndL%n(I%m:I%p,:) ! note minus sign
             I%dl   =  sqrt(Gp%xrL(I%m:I%p)**2+Gp%yrL(I%m:I%p)**2)
          end if
       case('y')
          if (.not.Gm%skip) then
             I%x    =  Gm%bndT%x(I%m:I%p)
             I%y    =  Gm%bndT%y(I%m:I%p)
             I%nhat =  Gm%bndT%n(I%m:I%p,:)
             I%dl   =  sqrt(Gm%xqT(I%m:I%p)**2+Gm%yqT(I%m:I%p)**2)
          else
             I%x    =  Gp%bndB%x(I%m:I%p)
             I%y    =  Gp%bndB%y(I%m:I%p)
             I%nhat = -Gp%bndB%n(I%m:I%p,:) ! note minus sign
             I%dl   =  sqrt(Gp%xqB(I%m:I%p)**2+Gp%yqB(I%m:I%p)**2)
          end if
       end select
       
    end if

  end subroutine init_iface_blocks


  subroutine init_iface_fields(iface,I,refine,input,echo,dt)

    use friction, only : init_friction
    use thermpres, only : init_thermpres
    use erupt, only : init_erupt

    implicit none

    integer,intent(in) :: iface,input,echo
    type(iface_type),intent(inout) :: I
    real,intent(in) :: refine,dt(:)

    ! initialize various interface models

    call init_friction( iface,I%FR,I%m,I%p,input,echo,I%skip, &
         I%process_m,I%process_p,I%comm_m,I%comm_p,I%array_w)
    call init_thermpres(iface,I%TP,I%m,I%p,input,echo,I%skip,refine,dt)
    if (I%skip) then
       ! commented this because I%x was not always allocated (DEBUG THIS)
       ! could use allocatable attribute when passing x,y,dl to init_erupt
       ! another idea is to always allocate I%x, but with zero length when I%skip=true
       !call init_erupt( iface,I%ER,I%m,I%p,I%x,I%y,input,echo,I%skip,I%direction,I%mg,I%pg, &
       !     I%process_m,I%process_p,I%comm_m,I%comm_p,I%array_w) ! otherwise error since dl not allocated
    else
       call init_erupt( iface,I%ER,I%m,I%p,I%x,I%y,input,echo,I%skip,I%direction,I%mg,I%pg, &
            I%process_m,I%process_p,I%comm_m,I%comm_p,I%array_w,I%dl)
    end if

  end subroutine init_iface_fields


  subroutine destroy_iface(I)

    use friction, only : destroy_friction
    use thermpres, only : destroy_thermpres
    use erupt, only : destroy_erupt

    implicit none

    type(iface_type),intent(inout) :: I

    if (allocated(I%x   )) deallocate(I%x   )
    if (allocated(I%y   )) deallocate(I%y   )
    if (allocated(I%nhat)) deallocate(I%nhat)
    if (allocated(I%dl  )) deallocate(I%dl  )

    call destroy_friction( I%FR)
    call destroy_thermpres(I%TP)
    call destroy_erupt(    I%ER)

  end subroutine destroy_iface
  

  subroutine checkpoint_iface(operation,name,checkpoint_number,iface,I)

    use friction, only : checkpoint_friction
    use thermpres, only : checkpoint_thermpres
    use erupt, only : checkpoint_erupt
    use io, only : file_distributed,open_file_distributed, &
         read_file_distributed,write_file_distributed,close_file_distributed
    use mpi_routines, only : pw,is_master
    use mpi

    implicit none

    character(*),intent(in) :: operation,name
    integer,intent(in) :: checkpoint_number,iface
    type(iface_type),intent(inout) :: I

    type(file_distributed) :: fh
    character(256) :: filename
    integer :: side,comm,ierr
    logical :: io_process

    ! must be called by all processes in MPI_COMM_WORLD,
    ! as master process deletes files

    do side = 1,2

       write(filename,'(a,i0,a,i0,a,i0)') trim(adjustl(name)) // &
            'iface',iface,'side',side,'.ckpt',checkpoint_number

       if (operation=='delete') then
          if (is_master) call MPI_File_delete(filename,MPI_INFO_NULL,ierr)
          if (I%coupling=='friction') call checkpoint_thermpres(operation,name, &
               checkpoint_number,iface,side,I%TP)
          if (side==2) return
          cycle
       end if

       select case(side)
       case(1)
          io_process = I%process_m
          if (io_process) call MPI_Comm_dup(I%comm_m,comm,ierr)
       case(2)
          io_process = I%process_p
          if (io_process) call MPI_Comm_dup(I%comm_p,comm,ierr)
       end select

       if (io_process) call open_file_distributed(fh,filename,operation,comm,I%array_w,pw)

       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       
       if (io_process) then
          if (I%coupling=='friction') call checkpoint_friction(fh,operation,I%FR)
          call checkpoint_erupt(fh,operation,I%ER)
       end if

       if (I%coupling=='friction') call checkpoint_thermpres(operation,name, &
            checkpoint_number,iface,side,I%TP,comm,io_process)

       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       
       if (io_process) then
          call close_file_distributed(fh)
          call MPI_Comm_free(comm,ierr)
       end if
       
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       
    end do

  end subroutine checkpoint_iface


  subroutine scale_rates_iface(I,A)

    use friction, only : scale_rates_friction
    use erupt, only : scale_rates_erupt

    implicit none

    type(iface_type),intent(inout) :: I
    real,intent(in) :: A

    if (I%skip) return

    select case(I%coupling)
    case('friction')
       call scale_rates_friction(I%FR,A)
    case('eruption')
       call scale_rates_erupt(   I%ER,A)
    end select

  end subroutine scale_rates_iface


  subroutine update_fields_iface(I,dt)

    use friction, only : update_fields_friction

    implicit none

    type(iface_type),intent(inout) :: I
    real,intent(in) :: dt

    if (I%skip) return

    select case(I%coupling)
    case('friction')
       call update_fields_friction(I%FR,dt)
    case('eruption')
       !JK call update_fields_erupt(I%ER,dt,I%FR%DDn(I%m:I%p))
    end select

  end subroutine update_fields_iface


  subroutine update_fields_iface_implicit(I,stage)

    use thermpres, only : update_fields_thermpres

    implicit none

    type(iface_type),intent(inout) :: I
    integer,intent(in) :: stage

    real,dimension(:),allocatable :: f,sigma

    if (I%skip) return

    !! JK: NEED TO FIX!
    if (I%TP%use_TP) then
       allocate(f(I%m:I%p),sigma(I%m:I%p))
       where (I%FR%N==0d0)
          f = 0d0
       elsewhere
          f = I%FR%S/I%FR%N
       end where
       sigma = I%FR%N+I%TP%p(1,:)
       call update_fields_thermpres(I%TP,I%m,I%p,f,I%FR%V,sigma,stage)
       deallocate(f,sigma)
    end if

  end subroutine update_fields_iface_implicit


  subroutine maxV(I,V)

    implicit none

    type(iface_type),intent(in) :: I
    real,intent(inout) :: V

    if (I%skip) return
    !! JK: NEED TO FIX!
    V = max(V,maxval(I%FR%V))

  end subroutine maxV


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

    namelist /boundaries_list/ bcL,bcR,bcB,bcT,rhog

    ! defaults

    bcL = 'none'
    bcR = 'none'
    bcB = 'none'
    bcT = 'none'
    rhog = 0d0

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

    ! output boundary condition parameters

    if (is_master) then
       write(Bstr,'(a,i0,a)') 'B{',iblock,'}'
       call write_matlab(echo,'bcL',B%bcL,Bstr)
       call write_matlab(echo,'bcR',B%bcR,Bstr)
       call write_matlab(echo,'bcB',B%bcB,Bstr)
       call write_matlab(echo,'bcT',B%bcT,Bstr)
       if (B%rhog/=0d0) call write_matlab(echo,'rhog',B%rhog,Bstr)
    end if

  end subroutine init_boundaries


  subroutine apply_bc(G,F,B,mode,t,iblock)

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
       if (G%sideL) call apply_bc_side(G%bndL,F%bndFL,G%my,G%py,B%bcL,B%rhog,mode,t,iblock)
       if (G%sideR) call apply_bc_side(G%bndR,F%bndFR,G%my,G%py,B%bcR,B%rhog,mode,t,iblock)
    end if

    ! r-direction

    if (G%ny/=1) then
       if (G%sideB) call apply_bc_side(G%bndB,F%bndFB,G%mx,G%px,B%bcB,B%rhog,mode,t,iblock)
       if (G%sideT) call apply_bc_side(G%bndT,F%bndFT,G%mx,G%px,B%bcT,B%rhog,mode,t,iblock)
    end if

  end subroutine apply_bc


  subroutine apply_bc_side(bndC,bndF,m,p,bc,rhog,mode,t,iblock)
    
    use geometry, only : curve
    use fields, only : bnd_fields

    implicit none
    
    type(curve),intent(in) :: bndC
    type(bnd_fields),intent(inout) :: bndF
    integer,intent(in) :: m,p
    character(*),intent(in) :: bc
    real,intent(in) :: rhog,t
    integer,intent(in) :: mode,iblock

    integer :: i

    if (bc=='none') return

    do i = m,p
       call set_bc(bndF%Fhat(i,:),bndF%F(i,:),bndF%F0(i,:),bndF%U(i,:),rhog, &
            bndC%x(i),bndC%y(i),bndC%n(i,:),bc,bndF%M(i,1:3),mode,t,iblock)
    end do
    
  end subroutine apply_bc_side


  subroutine set_bc(Fhat,F,F0,U,rhog,x,y,normal,bc,M,mode,t,iblock)

    implicit none

    ! F  = fields (grid data)
    ! FHAT = fields (hat variables)
    ! F0 = initial fields

    real,intent(out) :: Fhat(:)
    real,intent(in) :: F(:),F0(:),U(:),rhog,x,y,normal(2),M(3),t
    character(*),intent(in) :: bc
    integer,intent(in) :: mode,iblock

    select case(mode)
    case(2)
       call set_bc_mode2(Fhat,F,F0,U,rhog,normal,bc,M(1),M(2),M(3),x,y,t,iblock)
    case(3)
       call set_bc_mode3(Fhat,F,F0,normal,bc,M(1),x,y,t,iblock)
    end select

  end subroutine set_bc


  subroutine set_bc_mode3(Fhat,F,F0,normal,bc,Zs,x,y,t,iblock)

    use fields, only : rotate_fields_xy2nt,rotate_fields_nt2xy
    use mms, only : bessel
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
    case('bessel-v','bessel-w')
       call rotate_fields_xy2nt(F   ,normal,vzFD,stzFD,snzFD)
    end select

    ! then enforce BC

    select case(bc)
    case('bessel-v')
       vz = bessel(x,y,t,iblock,'vz')
       snz = snzFD-Zs*(vzFD-vz)
    case('bessel-w')
       Fe(1) = bessel(x,y,t,iblock,'vz')
       Fe(2) = bessel(x,y,t,iblock,'sxz')
       Fe(3) = bessel(x,y,t,iblock,'syz')
       call rotate_fields_xy2nt(Fe,normal,vzEX,stzEX,snzEX)
       vz = 0.5d0*(vzEX+vzFD+Zsi*(snzEX-snzFD))
       snz = 0.5d0*(snzEX+snzFD+Zs*(vzEX-vzFD))
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
    case('bessel-v','bessel-w')
       ! initial fields not required 
    end select

  end subroutine set_bc_mode3


  subroutine set_bc_mode2(Fhat,F,F0,U,rhog,normal,bc,Zs,Zp,gamma,x,y,t,iblock)

    use fields, only : rotate_fields_xy2nt,rotate_fields_nt2xy
    use mms, only : mms_sin,inplane_bessel,inplane_fault_mms
    use geometry, only : rotate_xy2nt
    use io, only : error

    implicit none

    real,intent(out) :: Fhat(6)
    real,intent(in) :: F(6),F0(6),U(2),rhog,normal(2),Zs,Zp,gamma,x,y,t
    character(*),intent(in) :: bc
    integer,intent(in) :: iblock

    real :: Zsi,Zpi,vn,vt,snn,snt,stt,szz,vnFD,vtFD,snnFD,sntFD,sttFD,szzFD
    real :: vnEX,vtEX,snnEX,sntEX,sttEX,szzEX,FEX(6),Ut,Un

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
    case('inplane-bessel-w','mms-sin-w','inplane-fault-mms-w','rigid-0','absorbing-0','free-0')
       call rotate_fields_xy2nt(F  ,normal,vtFD,vnFD,sttFD,sntFD,snnFD,szzFD)
    end select

    ! then enforce BC

    select case(bc)
       
    case('inplane-bessel-w','mms-sin-w','inplane-fault-mms-w') ! MMS
       
       select case(bc)
       case('inplane-bessel-w')
          FEX(1) = inplane_bessel(x,y,t,'vx')
          FEX(2) = inplane_bessel(x,y,t,'vy')
          FEX(3) = inplane_bessel(x,y,t,'sxx')
          FEX(4) = inplane_bessel(x,y,t,'sxy')
          FEX(5) = inplane_bessel(x,y,t,'syy')
          FEX(6) = 0d0
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
       end select
       call rotate_fields_xy2nt(FEX,normal,vtEX,vnEX,sttEX,sntEX,snnEX,szzEX)
       vt  = 0.5d0*(vtEX+vtFD+Zsi*(sntEX-sntFD))
       snt = 0.5d0*(sntEX+sntFD+Zs*(vtEX-vtFD))
       vn  = 0.5d0*(vnEX+vnFD+Zpi*(snnEX-snnFD))
       snn = 0.5d0*(snnEX+snnFD+Zp*(vnEX-vnFD))
       stt = sttFD-gamma*(snnFD-snn)
       szz = szzFD-gamma*(snnFD-snn)
       call rotate_fields_nt2xy(Fhat,normal,vt,vn,stt,snt,snn,szz)
       
    case('rigid-0','absorbing-0','free-0') 
       
       ! apply boundary conditions in local coordinates
       
       select case(bc)
       case('absorbing','absorbing-0')
          vt  = 0.5d0*( vtFD-Zsi*sntFD)
          snt = 0.5d0*(sntFD-Zs * vtFD)
          vn  = 0.5d0*( vnFD-Zpi*snnFD)
          snn = 0.5d0*(snnFD-Zp * vnFD)
       case('free','free-0')
          vt  = vtFD-Zsi*sntFD
          snt = 0d0
          vn  = vnFD-Zpi*snnFD
          snn = 0d0
       case('rigid','rigid-0')
          vt  = 0d0
          snt = sntFD-Zs*vtFD
          vn  = 0d0
          snn = snnFD-Zp*vnFD
       end select
       stt = sttFD-gamma*(snnFD-snn)
       szz = szzFD-gamma*(snnFD-snn)

    case('rigid','absorbing','free','absorbing-velocity','free-velocity','tsunami') 
       
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
    case('inplane-bessel-w','mms-sin-w','inplane-fault-mms-w','rigid-0','absorbing-0','free-0')
       ! initial fields not required 
    end select
    
  end subroutine set_bc_mode2


end module boundaries
