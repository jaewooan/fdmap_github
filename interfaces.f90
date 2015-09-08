module interfaces

  use friction, only : fr_type
  use thermpres, only : tp_type
  use hydrofrac, only : hf_type

  implicit none

  type :: iface_type
     logical :: skip,share,process_m,process_p,output_process
     character(24) :: coupling
     character(1) :: direction
     integer :: iblockm,iblockp,m,p,mg,pg,mb,pb, &
          rank_m,rank_p,comm_m,comm_p,comm_mp,comm,array_w,array_s
     type(fr_type) :: FR
     type(tp_type) :: TP
     type(hf_type) :: HF
     real,dimension(:),allocatable :: x,y,dl
     real,dimension(:,:),allocatable :: nhat
  end type iface_type


contains


  subroutine init_iface(iface,I,im,ip,Gm,Gp,dir,input,refine)

    use io, only : error,seek_to_string
    use grid, only : block_grid

    implicit none

    integer,intent(in) :: iface,input,im,ip
    real,intent(in) :: refine
    type(iface_type),intent(out) :: I
    type(block_grid),intent(in) :: Gm,Gp
    character(1),intent(in) :: dir
    logical :: read_params,no_refine
    integer :: stat,iblockm,iblockp,iblock,mg,pg
    character(24) :: coupling
    character(1) :: direction
    character(256) :: str,temp
    
    namelist /interface_list/ coupling,iblockm,iblockp,direction,mg,pg
    
    ! defaults

    coupling = 'locked'
    iblockm = im
    iblockp = ip
    direction = dir
    mg = 0
    pg = 0
    no_refine = .false.


    ! read in interface parameters if the interface
    ! is explicitly stated in the input file
        
    write(str,'(a,i0,a)') '!---IFACE',iface,'---'
    call seek_to_string_no_error(input,str,read_params)
    if(read_params) then
      read(input,nml=interface_list,iostat=stat)
      if (stat>0) call error('Error in interface_list','init_iface')
    end if
    
    ! Automatically determine mg and pg accounting for grid refinement
    if(mg == 0 .or. pg == 0) then
      select case(dir)
        case('x')
          mg = Gm%mgy
          pg = Gm%pgy
        case('y')
          mg = Gm%mgx
          pg = Gm%pgx
      end select
      no_refine = .true.
    end if


    I%coupling = coupling
    I%iblockm = iblockm
    I%iblockp = iblockp
    I%direction = direction
    I%mg = mg
    I%pg = pg

    if(no_refine) return

    ! adjust indices for refinement

    select case(I%direction)
    case default
       call error('Invalid direction','init_iface_blocks')
    case('x')
       iblock = Gm%iblock_y
    case('y')
       iblock = Gm%iblock_x
    end select

    ! Refine mg and pg that has been manually inputted
    if (mg==pg) then
       I%mg = ceiling(dble(mg-iblock)*refine)+iblock
       I%pg = I%mg
    else
       I%mg = ceiling(dble(mg-iblock)*refine)+iblock
       I%pg = ceiling(dble(pg-iblock)*refine)+iblock
    end if

  end subroutine init_iface
  
  subroutine seek_to_string_no_error(unit,str,stat)

    implicit none

    integer,intent(in) :: unit
    character(*),intent(in) :: str
    logical,intent(out) :: stat

    character(256) :: temp

    stat = .true.

    rewind(unit)

    do
       read(unit,*,end=100) temp
       if (temp==str) return
    end do

100 stat = .false.

  end subroutine seek_to_string_no_error


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


  subroutine init_iface_fields(iface,I,hmin,hmax,refine,input,echo,dt)

    use friction, only : init_friction
    use thermpres, only : init_thermpres
    use hydrofrac, only : init_hydrofrac

    implicit none

    integer,intent(in) :: iface,input,echo
    type(iface_type),intent(inout) :: I
    real,intent(in) :: hmin,hmax,refine,dt(:)

    ! initialize various interface models
    select case(I%coupling)
    case('friction')
      call init_friction( iface,I%FR,I%m,I%p,input,echo,I%skip, &
           I%process_m,I%process_p,I%comm_m,I%comm_p,I%array_w)
      call init_thermpres(iface,I%TP,I%m,I%p,input,echo,I%skip,refine,dt)
    case('hydrofrac')
     call init_hydrofrac(iface,I%HF,I%m,I%p,I%x,I%y,hmin,hmax,refine,input,echo,I%skip,I%direction,I%mg,I%pg, &
          I%process_m,I%process_p,I%comm_m,I%comm_p,I%array_w)
    end select

  end subroutine init_iface_fields


  subroutine destroy_iface(I)

    use friction, only : destroy_friction
    use thermpres, only : destroy_thermpres
    use hydrofrac, only : destroy_hydrofrac

    implicit none

    type(iface_type),intent(inout) :: I

    if (allocated(I%x   )) deallocate(I%x   )
    if (allocated(I%y   )) deallocate(I%y   )
    if (allocated(I%nhat)) deallocate(I%nhat)
    if (allocated(I%dl  )) deallocate(I%dl  )

    call destroy_friction( I%FR)
    call destroy_thermpres(I%TP)
    call destroy_hydrofrac(I%HF)

  end subroutine destroy_iface
  

  subroutine checkpoint_iface(operation,name,checkpoint_number,iface,I)

    use friction, only : checkpoint_friction
    use thermpres, only : checkpoint_thermpres
    use hydrofrac, only : checkpoint_hydrofrac
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
          select case(I%coupling)
          case('friction') 
             call checkpoint_friction(fh,operation,I%FR)
             !call checkpoint_thermpres(fh,operation,I%TP) ! NOT IMPLEMENTED YET
          case('hydrofrac')
             call checkpoint_hydrofrac(fh,operation,I%HF)
          end select
       end if

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
    use hydrofrac, only : scale_rates_hydrofrac

    implicit none

    type(iface_type),intent(inout) :: I
    real,intent(in) :: A

    if (I%skip) return

    select case(I%coupling)
    case('friction')
       call scale_rates_friction (I%FR,A)
    case('hydrofrac')
       call scale_rates_hydrofrac(I%HF,A)
    end select

  end subroutine scale_rates_iface


  subroutine update_fields_iface(I,dt)

    use friction, only : update_fields_friction
    use hydrofrac, only : update_fields_hydrofrac

    implicit none

    type(iface_type),intent(inout) :: I
    real,intent(in) :: dt

    if (I%skip) return

    select case(I%coupling)
    case('friction')
       call update_fields_friction (I%FR,dt)
    case('hydrofrac')
       call update_fields_hydrofrac(I%HF,dt)
    end select

  end subroutine update_fields_iface


  subroutine update_fields_iface_implicit(I,Fm,Fp,t,dt)

    use fields, only : block_fields
    !use thermpres, only : update_fields_thermpres
    use hydrofrac, only : update_fields_hydrofrac_implicit

    implicit none

    type(iface_type),intent(inout) :: I
    type(block_fields),intent(in) :: Fm,Fp
    real,intent(in) :: t,dt
    real,dimension(I%m:I%p) :: phip,vtp,vtm

    integer,parameter :: mode = 2

    !real,dimension(:),allocatable :: f,sigma

    if (I%skip) return

    select case(I%coupling)
    case('friction')
       !if (I%TP%use_TP) then
       !   allocate(f(I%m:I%p),sigma(I%m:I%p))
       !   where (I%FR%N==0d0)
       !      f = 0d0
       !   elsewhere
       !      f = I%FR%S/I%FR%N
       !   end where
       !   sigma = I%FR%N+I%TP%p(1,:)
       !   call update_fields_thermpres(I%TP,I%m,I%p,f,I%FR%V,sigma,stage)
       !   deallocate(f,sigma)
       !end if
    case('hydrofrac')
       select case(I%direction)
       case('x')
          ! x 
          call get_hydrofrac_interface(I%HF,I%m,I%p, &
               Fm%bndFR%F   (I%m:I%p, : ),Fp%bndFL%F   (I%m:I%p, : ), &
               Fm%bndFR%M   (I%m:I%p,1:3),Fp%bndFL%M   (I%m:I%p,1:3), &
               I%nhat(I%m:I%p,:),mode, &
               I%x(I%m:I%p),I%y(I%m:I%p),t,phip,vtp,vtm)
          call update_fields_hydrofrac_implicit(I%HF,I%m,I%p, &
               Fm%bndFR%M(I%m:I%p,2),Fp%bndFL%M(I%m:I%p,2), &
               phip,vtp,vtm,dt)
       case('y')
          ! y
          call get_hydrofrac_interface(I%HF,I%m,I%p, &
               Fm%bndFT%F   (I%m:I%p, : ),Fp%bndFB%F   (I%m:I%p, : ), &
               Fm%bndFT%M   (I%m:I%p,1:3),Fp%bndFB%M   (I%m:I%p,1:3), &
               I%nhat(I%m:I%p,:),mode, &
               I%x(I%m:I%p),I%y(I%m:I%p),t,phip,vtp,vtm)
          call update_fields_hydrofrac_implicit(I%HF,I%m,I%p, &
               Fm%bndFT%M(I%m:I%p,2),Fp%bndFB%M(I%m:I%p,2), &
               phip,vtp,vtm,dt)
       end select

    end select

  end subroutine update_fields_iface_implicit


  subroutine maxV(I,V)

    implicit none

    type(iface_type),intent(in) :: I
    real,intent(inout) :: V

    if (I%skip) return
    !! JK: NEED TO FIX!
    V = max(V,maxval(I%FR%V))

  end subroutine maxV


  subroutine enforce_interface_conditions(I,Fm,Fp,C,mode,t,dt)

    use fields, only : block_fields
    use mpi_routines2d, only : cartesian
    use friction, only : set_rates_friction
    use hydrofrac, only : set_rates_hydrofrac

    implicit none

    type(iface_type),intent(inout) :: I
    type(block_fields),intent(inout) :: Fm,Fp
    type(cartesian),intent(in) :: C
    integer,intent(in) :: mode
    real,intent(in) :: t,dt

    real,dimension(:),allocatable :: phip,vnm,vnp,vtm,vtp,sntm,sntp

    if (I%skip) return ! process has no grid points adjacent to interface

    ! set hat variables Fhat for different interface conditions

    select case(I%coupling)

    case('locked')

       ! locked (or welded) interface: no opening or slip

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

       ! enforce friction law by solving elasticity+friction for 
       ! slip velocity, shear stress, opening velocity (assumed zero for no opening),
       ! and normal stress; those interface fields are stored within FR (friction type)

       select case(I%direction)
       case('x')
          call enforce_frictional_interface(I%m,I%p, &
               Fm%bndFR%Fhat(I%m:I%p, : ),Fp%bndFL%Fhat(I%m:I%p, : ), &
               Fm%bndFR%F   (I%m:I%p, : ),Fp%bndFL%F   (I%m:I%p, : ), &
               Fm%bndFR%M   (I%m:I%p,1:3),Fp%bndFL%M   (I%m:I%p,1:3), &
               I%nhat(I%m:I%p,:),mode, &
               I%x(I%m:I%p),I%y(I%m:I%p),t,I%FR,I%TP,Fp%F0,dt)
       case('y')
          call enforce_frictional_interface(I%m,I%p, &
               Fm%bndFT%Fhat(I%m:I%p, : ),Fp%bndFB%Fhat(I%m:I%p, : ), &
               Fm%bndFT%F   (I%m:I%p, : ),Fp%bndFB%F   (I%m:I%p, : ), &
               Fm%bndFT%M   (I%m:I%p,1:3),Fp%bndFB%M   (I%m:I%p,1:3), &
               I%nhat(I%m:I%p,:),mode, &
               I%x(I%m:I%p),I%y(I%m:I%p),t,I%FR,I%TP,Fp%F0,dt)
       end select

       ! set rates for auxiliary fields (like slip and state variables),
       ! also calculate rupture time
          
       call set_rates_friction(I%FR,I%m,I%p,I%x(I%m:I%p),I%y(I%m:I%p),t)

    case('hydrofrac')

       allocate(phip(I%m:I%p))
       allocate(vnm(I%m:I%p))
       allocate(vnp(I%m:I%p))
       allocate(vtm(I%m:I%p))
       allocate(vtp(I%m:I%p))
       allocate(sntm(I%m:I%p))
       allocate(sntp(I%m:I%p))

       ! set hat variables (Fhat) by balancing fluid and solid tractions
       ! and matching fluid and solid velocities

       select case(I%direction)
       case('x')
          call enforce_hydrofrac_interface(I%m,I%p, &
               Fm%bndFR%Fhat(I%m:I%p, : ),Fp%bndFL%Fhat(I%m:I%p, : ), &
               Fm%bndFR%F   (I%m:I%p, : ),Fp%bndFL%F   (I%m:I%p, : ), &
               Fm%bndFR%M   (I%m:I%p,1:3),Fp%bndFL%M   (I%m:I%p,1:3), &
               I%nhat(I%m:I%p,:),mode, &
               I%x(I%m:I%p),I%y(I%m:I%p),t,I%HF,phip,vnm,vnp,vtm,vtp,sntm,sntp)
       case('y')

          ! Fm = block fields on minus side (use bndFT = top    boundary of this block)
          ! Fp = block fields on plus  side (use bndFB = bottom boundary of this block)
          ! Fhat = hat variables: 
          !        first index = I%m:I%p = range of indices along boundary, no ghost points)
          !        second index = fields index (vx,vy,sxx,sxy,syy,szz for plane strain)
          ! F = grid values of fields on boundaries, same indexing as Fhat
          ! M = elastic properties of material adjacent to boundary,
          !     second index: 1=Zs, 2=Zp, 3=gamma
          ! nhat = unit normal, second index=components of nhat in x and y directions
          ! mode = 2 for plane strain
          ! x,y = coordinates along boundary
          ! t = time
          ! HF = hydraulic fracture type

          call enforce_hydrofrac_interface(I%m,I%p, &
               Fm%bndFT%Fhat(I%m:I%p, : ),Fp%bndFB%Fhat(I%m:I%p, : ), &
               Fm%bndFT%F   (I%m:I%p, : ),Fp%bndFB%F   (I%m:I%p, : ), &
               Fm%bndFT%M   (I%m:I%p,1:3),Fp%bndFB%M   (I%m:I%p,1:3), &
               I%nhat(I%m:I%p,:),mode, &
               I%x(I%m:I%p),I%y(I%m:I%p),t,I%HF,phip,vnm,vnp,vtm,vtp,sntm,sntp)
          ! output P-wave stress transfer phip from above subroutine

       end select

       ! set rates for auxiliary fields (fluid velocity, pressure),
       ! passing in P-wave stress transfer phip, vnm, vnp

       call set_rates_hydrofrac(I%HF,C,I%m,I%p, &
           phip,vnm,vnp,vtm,vtp,I%x(I%m:I%p),I%y(I%m:I%p),t)

       deallocate(phip)

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


  subroutine enforce_frictional_interface(m,p,Fhatm,Fhatp,Fm,Fp,Mm,Mp,nhat,mode,x,y,t,FR,TP,F0,dt)

    use friction, only : fr_type
    use thermpres, only : tp_type

    implicit none

    integer,intent(in) :: m,p,mode
    real,dimension(m:,:),intent(out) :: Fhatm,Fhatp
    real,dimension(m:,:),intent(in) :: Fm,Fp,Mm,Mp,nhat
    real,dimension(m:p),intent(in) :: x,y
    real,intent(in) :: t,F0(9),dt
    type(fr_type),intent(inout) :: FR
    type(tp_type),intent(in) :: TP

    integer :: i

    do i = m,p
       select case(mode)
       case(2)
          call frictional_interface_mode2(Fhatm(i,:),Fhatp(i,:),Fm(i,:),Fp(i,:),nhat(i,:), &
               Mm(i,1),Mp(i,1),Mm(i,2),Mp(i,2),Mm(i,3),Mp(i,3),x(i),y(i),t,i,FR,TP,F0,dt)
       case(3)
          call frictional_interface_mode3(Fhatm(i,:),Fhatp(i,:),Fm(i,:),Fp(i,:),nhat(i,:), &
               Mm(i,1),Mp(i,1),x(i),y(i),t,i,FR,TP,F0,dt)
       end select
    end do
    
  end subroutine enforce_frictional_interface


  subroutine frictional_interface_mode3(Fhatm,Fhatp,Fm,Fp,normal,Zsm,Zsp,x,y,t,i,FR,TP,F0,dt)

    use geometry, only : rotate_xy2nt
    use fields, only : rotate_fields_xy2nt,rotate_fields_nt2xy
    use friction, only : fr_type,solve_friction
    use thermpres, only : tp_type,pressure_thermpres

    implicit none

    real,intent(out) :: Fhatm(3),Fhatp(3)
    real,intent(in) :: Fm(3),Fp(3),normal(2),Zsm,Zsp,x,y,t,F0(9),dt
    integer,intent(in) :: i
    type(fr_type),intent(inout) :: FR
    type(tp_type),intent(in) :: TP

    real :: Zsim,Zsip,phis,etas, &
         vzp,snzp,vzFDp,snzFDp,stzFDp, &
         vzm,snzm,vzFDm,snzFDm,stzFDm
    real :: stt,snt,snn,p,phip,N0

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
    ! note that as defined, phis includes prestress since snzFD does not have
    ! prestress subtracted out -- this is fine provided the traction components
    ! of prestress are continuous across the interface
    
    ! normal stress contribution from resolving in-plane stress field onto fault

    call rotate_xy2nt(F0(4),F0(5),F0(7),stt,snt,snn,normal)

    ! reduction of effective normal stress from pore pressure

    p = 0d0
    call pressure_thermpres(TP,i,p)
    N0 = -snn-p
    phip = 0d0
    ! here prestress-pore pressure is N0 and phip is only the stress transfer (0 for mode 3)

    ! solve for V,S,O,N by enforcing friction law and no opening
        
    call solve_friction(FR,FR%V(i),FR%S(i),FR%O(i),FR%N(i),N0,phip,phis,etas,FR%D(i),FR%Psi(i),i,x,y,t,.false.,dt)

    snzp = FR%S(i)-FR%S0(i)
    snzm = FR%S(i)-FR%S0(i)
    vzp = vzFDp+Zsip*(snzFDp-snzp)
    vzm = vzFDm-Zsim*(snzFDm-snzm)

    ! rotate back to x-y coordinates

    call rotate_fields_nt2xy(Fhatp,normal,vzp,stzFDp,snzp)
    call rotate_fields_nt2xy(Fhatm,normal,vzm,stzFDm,snzm)
    
  end subroutine frictional_interface_mode3


  subroutine frictional_interface_mode2(Fhatm,Fhatp,Fm,Fp,normal,Zsm,Zsp,Zpm,Zpp,gammam,gammap,x,y,t,i,FR,TP,F0,dt)

    use geometry, only : rotate_xy2nt
    use fields, only : rotate_fields_xy2nt,rotate_fields_nt2xy
    use friction, only : fr_type,solve_friction
    use thermpres, only : tp_type,pressure_thermpres

    implicit none

    real,intent(out) :: Fhatm(6),Fhatp(6)
    real,intent(in) :: Fm(6),Fp(6),normal(2),Zsm,Zsp,Zpm,Zpp,gammam,gammap,x,y,t,F0(9),dt
    integer,intent(in) :: i
    type(fr_type),intent(inout) :: FR
    type(tp_type),intent(in) :: TP

    real :: Zsim,Zsip,Zpim,Zpip,phis,phip,etas,etap, &
         vnp,vtp,snnp,sntp,sttp,szzp,vnFDp,vtFDp,snnFDp,sntFDp,sttFDp,szzFDp, &
         vnm,vtm,snnm,sntm,sttm,szzm,vnFDm,vtFDm,snnFDm,sntFDm,sttFDm,szzFDm
    real :: stt0,snt0,snn0,p,N0

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
    ! 4-6 yield equation of form snn = phip-Zp*O with
    ! phip = etap*(snnFDp/Zpp+snnFDm/Zpm+vnFDp-vnFDm),
    ! O = vnp-vnm, and 1/etap = 1/Zpp+1/Zpm

    ! rotate into local normal and tangential coordinates

    call rotate_fields_xy2nt(Fp,normal,vtFDp,vnFDp,sttFDp,sntFDp,snnFDp,szzFDp)
    call rotate_fields_xy2nt(Fm,normal,vtFDm,vnFDm,sttFDm,sntFDm,snnFDm,szzFDm)

    ! normal stress contribution from prestress

    call rotate_xy2nt(F0(4),F0(5),F0(7),stt0,snt0,snn0,normal)
    
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

    ! separate effective normal stress into prestress and stress transfer

    N0 = -snn0-p ! include pore pressure change here (so won't be subject to Skempton reduction)
    phip = phip-snn0 ! subtract out prestress

    ! solve for V,S,O,N by enforcing friction law and no opening

    call solve_friction(FR,FR%V(i),FR%S(i),FR%O(i),FR%N(i),N0,phip,phis,etas,FR%D(i),FR%Psi(i),i,x,y,t,.false.,dt)

    snnp = -(FR%N(i)-FR%N0(i)-N0)
    snnm = -(FR%N(i)-FR%N0(i)-N0)

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


  subroutine enforce_hydrofrac_interface(m,p,Fhatm,Fhatp,Fm,Fp,Mm,Mp,nhat,mode,x,y,t,HF,phip,vnm,vnp,vtm,vtp,sntm,sntp)

    use hydrofrac, only : hf_type
    use io, only : error

    implicit none

    integer,intent(in) :: m,p,mode
    real,dimension(m:,:),intent(out) :: Fhatm,Fhatp
    real,dimension(m:,:),intent(in) :: Fm,Fp,Mm,Mp,nhat
    real,dimension(m:p),intent(in) :: x,y
    real,intent(in) :: t
    type(hf_type),intent(inout) :: HF
    real,dimension(m:p),intent(out) :: phip,vnm,vnp,vtm,vtp,sntm,sntp

    integer :: i

    do i = m,p
       select case(mode)
       case(2)
          call hydrofrac_interface_mode2(Fhatm(i,:),Fhatp(i,:),Fm(i,:),Fp(i,:),nhat(i,:), &
               Mm(i,1),Mp(i,1),Mm(i,2),Mp(i,2),Mm(i,3),Mp(i,3),x(i),y(i),t,i,HF,phip(i),vnm(i),vnp(i),vtp(i),vtm(i),sntm(i),sntp(i))
       case(3)
          call error('No hydraulic fractures in antiplane shear','enforce_hydrofrac_interface')
       end select
    end do
    
  end subroutine enforce_hydrofrac_interface


  subroutine hydrofrac_interface_mode2(Fhatm,Fhatp,Fm,Fp,normal,Zsm,Zsp,Zpm,Zpp,&
                                       gammam,gammap,x,y,t,i,HF,phip,vnm,vnp,vtm,vtp,sntm,sntp)

    use fields, only : rotate_fields_xy2nt,rotate_fields_nt2xy
    use hydrofrac, only : hf_type,fluid_stresses,mms_velocities
    use mms, only : mms_hydrofrac

    implicit none

    real,intent(out) :: Fhatm(6),Fhatp(6),phip,vnm,vnp,vtm,vtp,sntm,sntp
    real,intent(in) :: Fm(6),Fp(6),normal(2),Zsm,Zsp,Zpm,Zpp,gammam,gammap,x,y,t
    integer,intent(in) :: i
    type(hf_type),intent(inout) :: HF

    real :: Zsim,Zsip,Zpim,Zpip, &
            snnp,sttp,szzp,vnFDp,vtFDp,snnFDp,sntFDp,sttFDp,szzFDp, &
            snnm,sttm,szzm,vnFDm,vtFDm,snnFDm,sntFDm,sttFDm,szzFDm
    real :: etap,taum,taup,p,wm,wp

    ! naming convention:
    ! vtFDp = tangential component of velocity, FD = grid value, p = plus side
    ! vtp = hat variable tangential component of velocity
    ! etc.

    ! 1. sntp+Zsp*vtp = sntFDp+Zsp*vtFDp 
    !    (S-wave into interface from plus  side)
    ! 2. sntm-Zsm*vtm = sntFDm-Zsm*vtFDm 
    !    (S-wave into interface from minus side)
    ! 3. sntp = taup (wall shear stress, plus  side)
    ! 4. sntm = taum (wall shear stress, minus side)
    !
    ! 5. snnp+Zpp*vnp = snnFDp+Zpp*vnFDp 
    !    (P-wave into interface from plus  side)
    ! 6. snnm-Zpm*vnm = snnFDm-Zpm*vnFDm 
    !    (P-wave into interface from minus side)
    ! 7. snnp = -p (fluid pressure)
    ! 8. snnm = -p (fluid pressure)
    !    (continuity of normal stress)

    ! rotate into local normal and tangential coordinates

    call rotate_fields_xy2nt(Fp,normal,vtFDp,vnFDp,sttFDp,sntFDp,snnFDp,szzFDp)
    call rotate_fields_xy2nt(Fm,normal,vtFDm,vnFDm,sttFDm,sntFDm,snnFDm,szzFDm)

    ! reciprocal impedances

    Zsip = 1d0/Zsp
    Zsim = 1d0/Zsm
    Zpip = 1d0/Zpp
    Zpim = 1d0/Zpm

    ! shear and normal tractions exerted by the fluid on the solid walls

    call fluid_stresses(HF,i,x,y,t,p,taum,taup)

    ! balance normal tractions

    snnp = -p
    snnm = -p

    vnp = vnFDp+Zpip*(snnFDp-snnp)
    vnm = vnFDm-Zpim*(snnFDm-snnm)
    
    sttp = sttFDp-gammap*(snnFDp-snnp)
    sttm = sttFDm-gammam*(snnFDm-snnm)
    
    szzp = szzFDp-gammap*(snnFDp-snnp)
    szzm = szzFDm-gammam*(snnFDm-snnm)
    
    ! balance shear tractions

    sntp = taup
    sntm = taum

    vtp = vtFDp+Zsip*(sntFDp-sntp)
    vtm = vtFDm-Zsim*(sntFDm-sntm)


    ! rotate back to x-y coordinates

    call rotate_fields_nt2xy(Fhatp,normal,vtp,vnp,sttp,sntp,snnp,szzp)
    call rotate_fields_nt2xy(Fhatm,normal,vtm,vnm,sttm,sntm,snnm,szzm)

    ! Send grid-values to fluid
    vtp = vtFDp
    vtm = vtFDm
    if(HF%use_mms) call mms_velocities(x,y,t,vtm,vtp)

    ! Exact tangential velocity used by mms
    !if(HF%use_mms) then
    !    vtm = mms_hydrofrac(x,0d0,t,0,'vtm')
    !    vtp = mms_hydrofrac(x,0d0,t,0,'vtp')
    !end if

    ! calculate P-wave stress transfer for use in implicit-explicit time-stepping

    !etap = 1d0/(Zpip+Zpim)
    !phip = etap*(snnFDp*Zpip+snnFDm*Zpim+vnFDp-vnFDm)
    
    phip =  snnFDp*Zpip+snnFDm*Zpim+vnFDp-vnFDm

  end subroutine hydrofrac_interface_mode2

  subroutine get_hydrofrac_interface(HF,m,p,Fm,Fp,Mm,Mp,nhat,mode,x,y,t, &
                                     phip,vtp,vtm)

    use hydrofrac, only : hf_type
    use io, only : error

    implicit none

    type(hf_type),intent(in) :: HF
    integer,intent(in) :: m,p,mode
    real,dimension(m:,:),intent(in) :: Fm,Fp,Mm,Mp,nhat
    real,dimension(m:p),intent(in) :: x,y
    real,intent(in) :: t
    real,dimension(m:p),intent(out) :: phip,vtp,vtm

    integer :: i

    do i = m,p
       select case(mode)
       case(2)
          call get_hydrofrac_interface_mode2(HF,Fm(i,:),Fp(i,:),nhat(i,:), &
               Mm(i,1),Mp(i,1),Mm(i,2),Mp(i,2),Mm(i,3),Mp(i,3),x(i),y(i),t, &
               phip(i),vtp(i),vtm(i))
       case(3)
          call error('No hydraulic fractures in antiplane shear','enforce_hydrofrac_interface')
       end select
    end do
    
 end subroutine get_hydrofrac_interface

 subroutine get_hydrofrac_interface_mode2(HF,Fm,Fp,normal,Zsm,Zsp,Zpm,Zpp,&
                                       gammam,gammap,x,y,t,phip,vtp,vtm)

    use fields, only : rotate_fields_xy2nt,rotate_fields_nt2xy
    use hydrofrac, only : hf_type,fluid_stresses,mms_velocities
    use mms, only : mms_hydrofrac

    implicit none

    type(hf_type),intent(in) :: HF
    real,intent(in) :: Fm(6),Fp(6),normal(2),Zsm,Zsp,Zpm,Zpp,gammam,gammap,x,y,t

    real :: Zpim,Zpip, &
            snnp,sttp,szzp,vnFDp,snnFDp,sntFDp,sttFDp,szzFDp, &
            snnm,sttm,szzm,vnFDm,snnFDm,sntFDm,sttFDm,szzFDm,vtFDp,vtFDm
    real,intent(out) :: phip,vtp, &
                             vtm

    ! rotate into local normal and tangential coordinates

    call rotate_fields_xy2nt(Fp,normal,vtFDp,vnFDp,sttFDp,sntFDp,snnFDp,szzFDp)
    call rotate_fields_xy2nt(Fm,normal,vtFDm,vnFDm,sttFDm,sntFDm,snnFDm,szzFDm)

    ! reciprocal impedances

    Zpip = 1d0/Zpp
    Zpim = 1d0/Zpm

    vtm = vtFDm
    vtp = vtFDp
    
    if(HF%use_mms) call mms_velocities(x,y,t,vtm,vtp)

    ! Exact tangential velocity used by mms
    !if(HF%use_mms) then
    !    vtm = mms_hydrofrac(x,0d0,t,0,'vtm')
    !    vtp = mms_hydrofrac(x,0d0,t,0,'vtp')
    !end if

    ! calculate P-wave stress transfer for use in implicit-explicit time-stepping
    
    phip =  snnFDp*Zpip+snnFDm*Zpim+vnFDp-vnFDm

 end subroutine get_hydrofrac_interface_mode2 


end module interfaces
