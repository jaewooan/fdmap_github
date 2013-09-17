module grid

  use geometry, only : point,curve

  implicit none

  type :: profile
     character(256) :: topo,filename
     real :: gamma,lambda
  end type profile

  type :: block_grid
     ! mgx,pgx are global limits; mx,px are process-specific limits
     integer :: nb=3, &
          iblock_x,mgx,pgx,nx,mx,px,mbx,pbx,mix,pix, &
          iblock_y,mgy,pgy,ny,my,py,mby,pby,miy,piy, &
          comm,myid,commx,commy
     logical :: skip,sideL,sideR,sideB,sideT
     character(256) :: bndLfile,bndRfile,bndBfile,bndTfile
     real :: bndLmap,bndRmap,bndBmap,bndTmap
     real,dimension(:),allocatable :: xrL,yrL,xrR,yrR,xqB,yqB,xqT,yqT,xyrL,xyrR,xyqB,xyqT
     type(point) :: LB,LT,RB,RT
     type(curve) :: bndL,bndR,bndB,bndT
  end type block_grid

  type :: grid_type
     real :: hmin=1d40,hmax=0d0
     real,dimension(:,:),allocatable :: x,y,xq,xr,yq,yr,J,Ji
  end type grid_type
  

contains


  subroutine read_grid(iblock,B,refine,input,echo)

    use io, only : seek_to_string,error,write_matlab,message
    use mpi_routines, only : is_master

    implicit none

    integer,intent(in) :: iblock,input,echo
    type(block_grid),intent(out) :: B ! block-specific
    real,intent(in) :: refine

    character(256) :: bndLfile,bndRfile,bndBfile,bndTfile
    character(256) :: str
    character(256) :: Bstr
    integer :: stat,iblock_x,iblock_y,nx,ny,mgx,mgy
    real :: bndLmap,bndRmap,bndBmap,bndTmap
    type(point) :: LB,LT,RB,RT

    namelist /grid_list/ iblock_x,iblock_y,nx,ny,mgx,mgy, &
         LB,LT,RB,RT, &
         bndLfile,bndRfile,bndBfile,bndTfile, &
         bndLmap,bndRmap,bndBmap,bndTmap

    ! defaults

    iblock_x = 0
    iblock_y = 0
    nx = 1
    ny = 1
    mgx = 1
    mgy = 1

    LB = point(huge(1d0),huge(1d0))
    LT = point(huge(1d0),huge(1d0))
    RB = point(huge(1d0),huge(1d0))
    RT = point(huge(1d0),huge(1d0))

    bndLfile = ''
    bndRfile = ''
    bndBfile = ''
    bndTfile = ''

    bndLmap = 0d0
    bndRmap = 0d0
    bndBmap = 0d0
    bndTmap = 0d0

    ! read in block-specific grid parameters

    write(str,'(a,i0,a)') '!---BLOCK',iblock,'---'
    call seek_to_string(input,str)
    read(input,nml=grid_list,iostat=stat)
    if (stat>0) call error('Error in grid_list','init_grid_block')

    write(str,'(i13)') iblock
    if (iblock_x==0) then
       iblock_x = 1
       if (is_master) call message('Block ' // trim(adjustl(str)) // &
            ': Assuming that iblock_x = 1')
    end if
    if (iblock_y==0) then
       if (is_master) call message('Block ' // trim(adjustl(str)) // &
            ': Assuming that iblock_y = 1')
       iblock_y = 1
    end if

    B%iblock_x = iblock_x
    B%iblock_y = iblock_y

    B%bndLfile = bndLfile
    B%bndRfile = bndRfile
    B%bndBfile = bndBfile
    B%bndTfile = bndTfile

    B%bndLmap = bndLmap
    B%bndRmap = bndRmap
    B%bndBmap = bndBmap
    B%bndTmap = bndTmap

    ! refine number of cells and indices of first point in block
    ! (ncells = npoints-1 for node-centered FD method)

    if (nx/=1) nx = ceiling(dble(nx-1)*refine)+1
    if (ny/=1) ny = ceiling(dble(ny-1)*refine)+1
    mgx = ceiling(dble(mgx-B%iblock_x)*refine)+B%iblock_x
    mgy = ceiling(dble(mgy-B%iblock_y)*refine)+B%iblock_y

    ! global (not process specific) block indices

    ! total number of points in block
    B%nx = nx
    B%ny = ny

    ! min/max indices
    B%mgx = mgx
    B%mgy = mgy
    B%pgx = mgx+nx-1
    B%pgy = mgy+ny-1

    ! corners of block

    B%LB = LB
    B%LT = LT
    B%RB = RB
    B%RT = RT
    
    ! write block information

    if (is_master) then
       write(Bstr,'(a,i0,a)') 'B{',iblock,'}'
       call write_matlab(echo,'LB.x',B%LB%x,Bstr)
       call write_matlab(echo,'LB.y',B%LB%y,Bstr)
       call write_matlab(echo,'LT.x',B%LT%x,Bstr)
       call write_matlab(echo,'LT.y',B%LT%y,Bstr)
       call write_matlab(echo,'RB.x',B%RB%x,Bstr)
       call write_matlab(echo,'RB.y',B%RB%y,Bstr)
       call write_matlab(echo,'RT.x',B%RT%x,Bstr)
       call write_matlab(echo,'RT.y',B%RT%y,Bstr)
       call write_matlab(echo,'mx',B%mgx,Bstr)
       call write_matlab(echo,'px',B%pgx,Bstr)
       call write_matlab(echo,'my',B%mgy,Bstr)
       call write_matlab(echo,'py',B%pgy,Bstr)
    end if

  end subroutine read_grid
  

  subroutine init_grid(iblock,B,G,C,FDmethod,exact_metric)

    use io, only : error,message
    use fd_coeff, only : nbst
    use mpi_routines, only : is_master,new_communicator,myid
    use mpi_routines2d, only : cartesian
    use mpi

    implicit none

    integer,intent(in) :: iblock
    type(block_grid),intent(inout) :: B ! block-specific
    type(grid_type),intent(inout) :: G ! global grid coordinates
    type(cartesian),intent(in) :: C ! global grid indices
    character(*),intent(in) :: FDmethod
    logical, intent(in) :: exact_metric

    character(256) :: str,str2
    integer :: ierr

    ! check block limits

    if (B%mgx<1.or.C%nx<B%pgx.or.B%mgy<1.or.C%ny<B%pgy) then
       write(str2,'(a,i6,a)') 'Block ',iblock,' indices too small or large'
       call error(str2,'init_grid')
    end if

    if (is_master) then
       write(str2,'(a,i0,a,2i6,a,2i6)') &
            'block',iblock,':  i=',B%mgx,B%pgx,',  j=',B%mgy,B%pgy
       call message(str2)
    end if

    ! process-specific block indices

    B%mx = max(B%mgx,C%mx)
    B%px = min(B%pgx,C%px)
    B%my = max(B%mgy,C%my)
    B%py = min(B%pgy,C%py)

    B%mbx = B%mx-B%nb
    B%pbx = B%px+B%nb
    B%mby = B%my-B%nb
    B%pby = B%py+B%nb

    ! interior indices (excluding boundary regions using one-sided FD stencils)
    B%mix = max(B%mx,B%mgx+nbst)
    B%pix = min(B%px,B%pgx-nbst)
    B%miy = max(B%my,B%mgy+nbst)
    B%piy = min(B%py,B%pgy-nbst)

    ! determine if process has any data for this block (if not, set skip=T)
    
    B%skip = .not.(B%mx<=B%px.and.B%my<=B%py)

    ! check that process has enough points in each block for FD method

    if (.not.B%skip) then
       write(str,'(i13)') iblock
       select case(FDmethod)
       case('SBP2','UPW1')
          if (B%px-B%mx+1<3 .and.B%nx/=1) &
               call error('Process has too few x points in block ' // trim(str),'init_grid')
          if (B%py-B%my+1<3 .and.B%ny/=1) &
               call error('Process has too few y points in block ' // trim(str),'init_grid')
       case('SBP4','UPW3')
          if (B%px-B%mx+1<8 .and.B%nx/=1) &
               call error('Process has too few x points in block ' // trim(str),'init_grid')
          if (B%py-B%my+1<8 .and.B%ny/=1) &
               call error('Process has too few y points in block ' // trim(str),'init_grid')
       case('SBP6','UPW5')
          if (B%px-B%mx+1<12.and.B%nx/=1) &
               call error('Process has too few x points in block ' // trim(str),'init_grid')
          if (B%py-B%my+1<12.and.B%ny/=1) &
               call error('Process has too few y points in block ' // trim(str),'init_grid')
       end select
    end if

    ! determine if process contains boundaries and/or fault

    B%sideL = (B%mx==B%mgx).and.(.not.B%skip)
    B%sideR = (B%px==B%pgx).and.(.not.B%skip)
    B%sideB = (B%my==B%mgy).and.(.not.B%skip)
    B%sideT = (B%py==B%pgy).and.(.not.B%skip)

    ! allocate memory for distributed coordinate arrays

    if (.not.allocated(G%x)) then
       allocate(G%x(C%mbx:C%pbx,C%mby:C%pby))
       allocate(G%y(C%mbx:C%pbx,C%mby:C%pby))
       G%x = 1d40 ! for debugging
       G%y = 1d40 ! for debugging
    end if

    if (.not.allocated(G%xq)) then
       allocate(G%xq(C%mbx:C%pbx,C%mby:C%pby))
       allocate(G%xr(C%mbx:C%pbx,C%mby:C%pby))
       allocate(G%yq(C%mbx:C%pbx,C%mby:C%pby))
       allocate(G%yr(C%mbx:C%pbx,C%mby:C%pby))
       allocate(G%J (C%mbx:C%pbx,C%mby:C%pby))
       allocate(G%Ji(C%mbx:C%pbx,C%mby:C%pby))
       G%xq = 1d0
       G%xr = 0d0
       G%yq = 0d0
       G%yr = 1d0
       G%J  = 1d0
       G%Ji = 1d0
    end if

    ! allocate memory to edges (and initialize to huge values to prevent inadvertent use)

    if (.not.B%skip) then

       allocate(B%bndL%x(B%mby:B%pby),B%bndL%y(B%mby:B%pby),B%bndL%xt(B%mby:B%pby),B%bndL%yt(B%mby:B%pby),B%bndL%n(B%my:B%py,2))
       allocate(B%bndR%x(B%mby:B%pby),B%bndR%y(B%mby:B%pby),B%bndR%xt(B%mby:B%pby),B%bndR%yt(B%mby:B%pby),B%bndR%n(B%my:B%py,2))
       allocate(B%bndB%x(B%mbx:B%pbx),B%bndB%y(B%mbx:B%pbx),B%bndB%xt(B%mbx:B%pbx),B%bndB%yt(B%mbx:B%pbx),B%bndB%n(B%mx:B%px,2))
       allocate(B%bndT%x(B%mbx:B%pbx),B%bndT%y(B%mbx:B%pbx),B%bndT%xt(B%mbx:B%pbx),B%bndT%yt(B%mbx:B%pbx),B%bndT%n(B%mx:B%px,2))

       B%bndL%x = 1d40
       B%bndL%y = 1d40
       B%bndL%xt = 1d40
       B%bndL%yt = 1d40
       B%bndL%n = 1d40
       B%bndR%x = 1d40
       B%bndR%y = 1d40
       B%bndR%xt = 1d40
       B%bndR%yt = 1d40
       B%bndR%n = 1d40
       B%bndB%x = 1d40
       B%bndB%y = 1d40
       B%bndB%xt = 1d40
       B%bndB%yt = 1d40
       B%bndB%n = 1d40
       B%bndT%x = 1d40
       B%bndT%y = 1d40
       B%bndT%xt = 1d40
       B%bndT%yt = 1d40
       B%bndT%n = 1d40

    end if

    if (B%sideL) then
       allocate(B%xrL(B%my:B%py),B%yrL(B%my:B%py),B%xyrL(B%my:B%py))
       B%xrL = 1d40
       B%yrL = 1d40
       B%xyrL = 1d40
    end if
    if (B%sideR) then
       allocate(B%xrR(B%my:B%py),B%yrR(B%my:B%py),B%xyrR(B%my:B%py))
       B%xrR = 1d40
       B%yrR = 1d40
       B%xyrR = 1d40
    end if
    if (B%sideB) then
       allocate(B%xqB(B%mx:B%px),B%yqB(B%mx:B%px),B%xyqB(B%mx:B%px))
       B%xqB = 1d40
       B%yqB = 1d40
       B%xyqB = 1d40
    end if
    if (B%sideT) then
       allocate(B%xqT(B%mx:B%px),B%yqT(B%mx:B%px),B%xyqT(B%mx:B%px))
       B%xqT = 1d40
       B%yqT = 1d40
       B%xyqT = 1d40
    end if

    ! block communicator (2D)

    call new_communicator(.not.B%skip,B%comm,C%c2d%comm)
    if (.not.B%skip) call MPI_Comm_rank(B%comm,B%myid,ierr)

    ! return if process does not manipulate interior of this block
    
    if (B%skip) return

    ! x and y direction (1D) communicators within block
    ! same value of coord(2) means same y-position => same commx
    ! same value of coord(1) means same x-position => same commy

    call MPI_Comm_split(B%comm,C%c2d%coord(2),B%myid,B%commx,ierr)
    call MPI_Comm_split(B%comm,C%c2d%coord(1),B%myid,B%commy,ierr)

    ! initialize edges and interpolate interior points

    call interp_block(B,G,exact_metric)

  end subroutine init_grid
  

  subroutine block_limits(B,mgx,mgy,pgx,pgy)

    type(block_grid),intent(in) :: B
    integer,intent(out) :: mgx,mgy
    integer,intent(out),optional :: pgx,pgy

    mgx = B%mgx
    mgy = B%mgy
    if (present(pgx)) pgx = B%pgx
    if (present(pgy)) pgy = B%pgy

  end subroutine block_limits


  subroutine interp_block(B,G,exact_metric)

    use geometry, only : line_map,transfinite_interp

    implicit none

    type(block_grid),intent(inout) :: B
    type(grid_type),intent(inout) :: G
    logical, intent(in) :: exact_metric

    integer :: i,j
    real :: qmin,qmax,rmin,rmax
    real :: q(B%mx:B%px),r(B%my:B%py)
    real,dimension(2),parameter :: xhat = (/1d0,0d0/), yhat = (/0d0,1d0/)

    ! return if process does not manipulate interior of this block

    if (B%skip) return

    ! q,r coordinate vectors with unit spacing

    qmin = dble(B%mgx)
    qmax = dble(B%pgx)
    rmin = dble(B%mgy)
    rmax = dble(B%pgy)

    q = dble( (/ (i, i=B%mx,B%px) /) )
    r = dble( (/ (j, j=B%my,B%py) /) )
    
    ! edges

    if (B%ny==1) then
       B%bndL = curve((/B%LB%x/),(/B%LB%y/),(/1d40/),(/1d40/),reshape(-xhat,(/1,2/)))
       B%bndR = curve((/B%RB%x/),(/B%RB%y/),(/1d40/),(/1d40/),reshape( xhat,(/1,2/)))
    else
       if (B%bndLfile=='') then
          call line_map(B%bndL,B%LB,B%LT,B%my,B%py,r,rmin,rmax,B%bndLmap)
       else
          call read_bnd(B%bndL,B,trim(B%bndLfile),'x',B%commy,exact_metric)
       end if
       B%bndL%n = -B%bndL%n ! fix sign for outward normal
       if (B%bndRfile=='') then
          call line_map(B%bndR,B%RB,B%RT,B%my,B%py,r,rmin,rmax,B%bndRmap)
       else
          call read_bnd(B%bndR,B,trim(B%bndRfile),'x',B%commy,exact_metric)
       end if
    end if

    if (B%nx==1) then
       B%bndB = curve((/B%LB%x/),(/B%LB%y/),(/1d40/),(/1d40/),reshape(-yhat,(/1,2/)))
       B%bndT = curve((/B%LT%x/),(/B%LT%y/),(/1d40/),(/1d40/),reshape( yhat,(/1,2/)))
    else
       if (B%bndBfile=='') then
          call line_map(B%bndB,B%LB,B%RB,B%mx,B%px,q,qmin,qmax,B%bndBmap)
       else
          call read_bnd(B%bndB,B,trim(B%bndBfile),'y',B%commx,exact_metric)
       end if
       if (B%bndTfile=='') then
          call line_map(B%bndT,B%LT,B%RT,B%mx,B%px,q,qmin,qmax,B%bndTmap)
       else
          call read_bnd(B%bndT,B,trim(B%bndTfile),'y',B%commx,exact_metric)
       end if
       B%bndT%n = -B%bndT%n ! fix sign for outward normal
    end if

    ! transfinite interpolation for interior region

    call transfinite_interp(B%mx,B%px,B%my,B%py,qmin,qmax,rmin,rmax, &
       q,r,B%bndL,B%bndR,B%bndB,B%bndT,B%LB,B%LT,B%RB,B%RT, &
       G%x(B%mx:B%px,B%my:B%py),G%y(B%mx:B%px,B%my:B%py))

  end subroutine interp_block


  subroutine init_grid_partials(B,G,exact_metric)

    use fd, only : limits,diff
    use geometry, only : transfinite_interp_derivative

    implicit none

    type(block_grid),intent(inout) :: B ! block-specific
    type(grid_type),intent(inout) :: G ! global grid coordinates
    logical,intent(in) :: exact_metric

    integer :: i,j
    real :: dx,dy
    type(limits) :: L

    real :: qmin,qmax,rmin,rmax
    real :: q(B%mx:B%px),r(B%my:B%py)


    ! return if process does not manipulate interior of this block

    if (B%skip) return

    if(exact_metric) then

      qmin = dble(B%mgx)
      qmax = dble(B%pgx)
      rmin = dble(B%mgy)
      rmax = dble(B%pgy)

      q = dble( (/ (i, i=B%mx,B%px) /) )
      r = dble( (/ (j, j=B%my,B%py) /) )

      call transfinite_interp_derivative(B%mx,B%px,B%my,B%py,qmin,qmax,rmin,rmax, &
        q,r,B%bndL,B%bndR,B%bndB,B%bndT,B%LB,B%LT,B%RB,B%RT, &
        G%xq(B%mx:B%px,B%my:B%py),G%xr(B%mx:B%px,B%my:B%py),&
        G%yq(B%mx:B%px,B%my:B%py),G%yr(B%mx:B%px,B%my:B%py))

    else
      ! q-derivatives

      if (B%nx==1) then
        G%xq(B%mx:B%px,B%my:B%py) = 1d0
        G%yq(B%mx:B%px,B%my:B%py) = 0d0
      else
        call set_block_limits(B,L,'x')
        do j = B%my,B%py
          call diff(L,G%x(B%mbx:B%pbx,j),G%xq(B%mx:B%px,j))
          call diff(L,G%y(B%mbx:B%pbx,j),G%yq(B%mx:B%px,j))
        end do
      end if

      ! r-derivatives

      if (B%ny==1) then
        G%xr(B%mx:B%px,B%my:B%py) = 0d0
        G%yr(B%mx:B%px,B%my:B%py) = 1d0
      else
        call set_block_limits(B,L,'y')
        do i = B%mx,B%px
          call diff(L,G%x(i,B%mby:B%pby),G%xr(i,B%my:B%py))
          call diff(L,G%y(i,B%mby:B%pby),G%yr(i,B%my:B%py))
        end do
      end if
    end if

    G%J(B%mx:B%px,B%my:B%py) = &
         G%xq(B%mx:B%px,B%my:B%py)*G%yr(B%mx:B%px,B%my:B%py)- &
         G%xr(B%mx:B%px,B%my:B%py)*G%yq(B%mx:B%px,B%my:B%py)
    G%Ji(B%mx:B%px,B%my:B%py) = 1d0/G%J(B%mx:B%px,B%my:B%py)

    ! effective grid spacings, used to set time step with CFL condition

    if     (B%nx==1) then
       
       do j = B%my,B%py
          dy = abs(G%J(B%mx,j)/G%xq(B%mx,j))
          G%hmin = min(G%hmin,dy)
          G%hmax = max(G%hmax,dy)
       end do

    elseif (B%ny==1) then

       do i = B%mx,B%px
          dx = abs(G%J(i,B%my)/G%yr(i,B%my))
          G%hmin = min(G%hmin,dx)
          G%hmax = max(G%hmax,dx)
       end do

    else
       
       do j = B%my,B%py
          do i = B%mx,B%px
             dx = abs(G%J(i,j)/sqrt(G%xr(i,j)**2+G%yr(i,j)**2))
             dy = abs(G%J(i,j)/sqrt(G%xq(i,j)**2+G%yq(i,j)**2))
             G%hmin = min(G%hmin,dx,dy)
             ! if(G%hmin < 4d-6) then
             !   print*,i,j
             ! end if
             G%hmax = max(G%hmax,dx,dy)
         end do
       end do

    end if

    ! partial derivatives on boundaries (after filling ghost cells)

    if (B%sideL) then
       if (B%ny==1) then
          B%xrL = 0d0
          B%yrL = 1d0
          B%xyrL = 1d0
       else
          B%xrL = G%xr(B%mgx,B%my:B%py)
          B%yrL = G%yr(B%mgx,B%my:B%py)
          B%xyrL = sqrt(B%xrL**2 + B%yrL**2)
       end if

    end if

    if (B%sideR) then
       if (B%ny==1) then
          B%xrR = 0d0
          B%yrR = 1d0
          B%xyrR = 1d0
       else
          B%xrR = G%xr(B%pgx,B%my:B%py)
          B%yrR = G%yr(B%pgx,B%my:B%py)
          B%xyrR = sqrt(B%xrR**2 + B%yrR**2)
       end if
    end if

    if (B%sideB) then
       if (B%nx==1) then
          B%xqB = 1d0
          B%yqB = 0d0
          B%xyqB = 1d0
       else
          B%xqB = G%xq(B%mx:B%px,B%mgy)
          B%yqB = G%yq(B%mx:B%px,B%mgy)
          B%xyqB = sqrt(B%xqB**2 + B%yqB**2)
       end if
    end if

    if (B%sideT) then
       if (B%nx==1) then
          B%xqT = 1d0
          B%yqT = 0d0
          B%xyqT = 1d0
       else
          B%xqT = G%xq(B%mx:B%px,B%pgy)
          B%yqT = G%yq(B%mx:B%px,B%pgy)
          B%xyqT = sqrt(B%xqT**2 + B%yqT**2)
       end if
    end if

  end subroutine init_grid_partials


  subroutine grid_spacing(C,G)

    use io, only : error
    use mpi_routines2d, only : cartesian
    use mpi_routines, only : MPI_REAL_PW
    use mpi

    implicit none

    type(cartesian),intent(in) :: C
    type(grid_type),intent(inout) :: G

    integer :: ierr
    real :: hmin,hmax
    
    ! minimum/maximum grid spacing over all processors

    !write(6,*) 'before',G%hmin,G%hmax

    ! MPI_IN_PLACE broken on OpenMPI for Mac
    !call MPI_Allreduce(MPI_IN_PLACE,G%hmin,1,MPI_REAL_PW,MPI_MIN,C%c2d%comm,ierr)
    !call MPI_Allreduce(MPI_IN_PLACE,G%hmax,1,MPI_REAL_PW,MPI_MAX,C%c2d%comm,ierr)

    hmin = G%hmin
    hmax = G%hmax
    call MPI_Allreduce(hmin,G%hmin,1,MPI_REAL_PW,MPI_MIN,C%c2d%comm,ierr)
    call MPI_Allreduce(hmax,G%hmax,1,MPI_REAL_PW,MPI_MAX,C%c2d%comm,ierr)

    if (G%hmin<=0d0) call error('Minimum grid spacing must be positive','grid_spacing')

    !write(6,*) 'after ',G%hmin,G%hmax

  end subroutine grid_spacing


  subroutine destroy_grid(G,B)

    use geometry, only : destroy_curve

    implicit none

    type(grid_type),intent(inout) :: G
    type(block_grid),intent(inout) :: B

    if (allocated(G%x)) deallocate(G%x,G%y)
    if (allocated(G%xq)) deallocate(G%xq,G%xr,G%yq,G%yr,G%J,G%Ji)

    call destroy_curve(B%bndL)
    call destroy_curve(B%bndR)
    call destroy_curve(B%bndB)
    call destroy_curve(B%bndT)

    if (allocated(B%xrL)) deallocate(B%xrL,B%yrL,B%xyrL)
    if (allocated(B%xrR)) deallocate(B%xrR,B%yrR,B%xyrR)
    if (allocated(B%xqB)) deallocate(B%xqB,B%yqB,B%xyqB)
    if (allocated(B%xqT)) deallocate(B%xqT,B%yqT,B%xyqT)

  end subroutine destroy_grid


  subroutine read_bnd(C,B,filename,direction,comm,exact_metric)

    use io, only : file_distributed,open_file_distributed, &
         read_file_distributed,close_file_distributed,error
    use mpi_routines, only : pw,MPI_REAL_PW,subarray
    use mpi

    implicit none

    type(curve),intent(inout) :: C
    type(block_grid),intent(in) :: B
    character(*),intent(in) :: filename,direction
    integer,intent(in) :: comm
    logical,intent(in) :: exact_metric

    integer :: m,p,n,mg,array,ierr
    type(file_distributed) :: fh

    select case(direction)
    case default
       call error('Invalid direction','read_bnd')
    case('x')
       m = B%my
       p = B%py
       n = B%ny
       mg = B%mgy
    case('y')
       m = B%mx
       p = B%px
       n = B%nx
       mg = B%mgx
    end select

    if (m>p) call error('CHECK WHY THIS CASE ARISES','read_bnd')

    ! distributed read of binary file

    call subarray(n,m-mg+1,p-mg+1,MPI_REAL_PW,array)

    call open_file_distributed(fh,filename,'read',comm,array,pw)
    call read_file_distributed(fh,C%x(m:p))
    call read_file_distributed(fh,C%y(m:p))
    call read_file_distributed(fh,C%n(m:p,1))
    call read_file_distributed(fh,C%n(m:p,2))
    if(exact_metric) then
      call read_file_distributed(fh,C%xt(m:p))
      call read_file_distributed(fh,C%yt(m:p))
      C%xt(m:p) = C%xt(m:p)/(dble(n-1))
      C%yt(m:p) = C%yt(m:p)/(dble(n-1))
    end if
    call close_file_distributed(fh)

    ! clean up

    call MPI_Type_free(array,ierr)

  end subroutine read_bnd


  subroutine set_block_limits(B,L,direction)

    use fd, only : limits
    use utilities, only : within
    use io, only : error

    implicit none

    type(block_grid),intent(in) :: B
    type(limits),intent(out) :: L
    character(*),intent(in) :: direction

    L%nb = B%nb

    select case(direction)
    case default
       call error('Invalid direction','set_block_limits')
    case('x')
       L%m = B%mx
       L%p = B%px
       L%mg = B%mgx
       L%pg = B%pgx
    case('y')
       L%m = B%my
       L%p = B%py
       L%mg = B%mgy
       L%pg = B%pgy
    end select

    L%mb = L%m-L%nb
    L%pb = L%p+L%nb

    L%endm = within(L%m,L%mg,L%p)
    L%endp = within(L%m,L%pg,L%p)

  end subroutine set_block_limits


end module grid
