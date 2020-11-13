module source

  use io, only : file_distributed
  
  implicit none

  type :: source_type
     logical :: use_gravity,use_singular_source,forcing_from_file
     character(256) :: source_time_function,problem,forcing_filename
     real :: gx,gy,gz ! components of gravitational acceleration
     real :: x0,y0,t_duration,t_center,Fx,Fy,Fz,Mxx,Mxy,Myy,Mzx,Mzy ! point force and moment tensor
     real :: gauss_amp,gauss_sigma,gauss_sigmat ! tsunami forcing, Gaussian in space and time
     real,dimension(:,:,:),allocatable :: forcing_old,forcing_new ! time-dependent forcing from file
     real :: t_old,t_new ! time-dependent forcing from file
     type(file_distributed) :: fh
  end type source_type

contains


  subroutine init_source(F,C,checkpoint_number,S,input,echo)

    use mpi_routines, only : is_master
    use mpi_routines2d, only : cartesian
    use io, only : error,write_matlab
    use fields, only : fields_type

    implicit none

    type(fields_type),intent(in) :: F
    type(cartesian),intent(in) :: C
    integer,intent(in) :: checkpoint_number
    type(source_type),intent(out) :: S
    integer,intent(in) :: input,echo

    logical :: forcing_from_file
    real :: gx,gy,gz,x0,y0,t_duration,t_center,Fx,Fy,Fz,Mxx,Mxy,Myy,Mzx,Mzy,gauss_amp,gauss_sigma,gauss_sigmat
    character(256) :: source_time_function,problem,forcing_filename
    integer :: stat

    namelist /source_list/ problem,gauss_amp,gauss_sigma,gauss_sigmat, &
         gx,gy,gz,x0,y0,t_duration,t_center,Fx,Fy,Fz,Mxx,Mxy,Myy,Mzx,Mzy,source_time_function, &
         forcing_from_file,forcing_filename
         
    ! defaults

    gx = 0d0
    gy = 0d0
    gz = 0d0

    x0 = 0d0
    y0 = 0d0
    t_duration = 1d0
    t_center = 0d0

    Fx = 0d0
    Fy = 0d0
    Fz = 0d0

    Mxx = 0d0
    Mxy = 0d0
    Myy = 0d0
    Mzx = 0d0
    Mzy = 0d0

    source_time_function = 'gaussian'

    problem = ''
    gauss_amp = 0d0
    gauss_sigma = 0d0
    gauss_sigmat = 0d0

    forcing_from_file = .false.
    forcing_filename = ''
    
    ! read in source parameters
        
    rewind(input)
    read(input,nml=source_list,iostat=stat)
    if (stat>0) call error('Error in source_list','init_source')

    S%use_gravity = (gx/=0d0.or.gy/=0d0.or.gz/=0d0)
    S%use_singular_source = (Fx/=0d0.or.Fy/=0d0.or.Fz/=0d0).or. &
         (Mxx/=0d0.or.Mxy/=0d0.or.Myy/=0d0.or.Mzx/=0d0.or.Mzy/=0d0)

    S%gx = gx
    S%gy = gy
    S%gz = gz

    S%x0 = x0
    S%y0 = y0
    S%t_duration = t_duration
    S%t_center = t_center
    
    S%Fx = Fx
    S%Fy = Fy
    S%Fz = Fz

    S%Mxx = Mxx
    S%Mxy = Mxy
    S%Myy = Myy
    S%Mzx = Mzx
    S%Mzy = Mzy

    S%source_time_function = source_time_function

    S%problem = problem
    S%gauss_amp = gauss_amp
    S%gauss_sigma = gauss_sigma
    S%gauss_sigmat = gauss_sigmat

    S%forcing_from_file = forcing_from_file
    S%forcing_filename = forcing_filename

    if (S%forcing_from_file) call init_forcing_from_file(F,C,S,checkpoint_number)
    
    ! output source parameters

    if (is_master) then
       call write_matlab(echo,'gx',S%gx,'S')
       call write_matlab(echo,'gy',S%gy,'S')
       call write_matlab(echo,'gz',S%gz,'S')
       call write_matlab(echo,'x0',S%x0,'S')
       call write_matlab(echo,'y0',S%y0,'S')
       call write_matlab(echo,'t_duration',S%t_duration,'S')
       call write_matlab(echo,'t_center',S%t_center,'S')
       call write_matlab(echo,'Fx',S%Fx,'S')
       call write_matlab(echo,'Fy',S%Fy,'S')
       call write_matlab(echo,'Fz',S%Fz,'S')
       call write_matlab(echo,'Mxx',S%Mxx,'S')
       call write_matlab(echo,'Mxy',S%Mxy,'S')
       call write_matlab(echo,'Myy',S%Myy,'S')
       call write_matlab(echo,'Mzx',S%Mzx,'S')
       call write_matlab(echo,'Mzy',S%Mzy,'S')
       call write_matlab(echo,'source_time_function',S%source_time_function,'S')
       call write_matlab(echo,'problem',S%problem,'S')
       call write_matlab(echo,'gauss_amp',S%gauss_amp,'S')
       call write_matlab(echo,'gauss_sigmat',S%gauss_sigmat,'S')
       call write_matlab(echo,'gauss_sigma',S%gauss_sigma,'S')
       call write_matlab(echo,'forcing_from_file',S%forcing_from_file,'S')
       call write_matlab(echo,'forcing_filename',S%forcing_filename,'S')
    end if

  end subroutine init_source


  subroutine set_source(B,G,F,M,S,t,mode,iblock)
    
    use io, only : error
    use grid, only : grid_type,block_grid
    use fields, only : fields_type
    use mms, only : mms_sin,inplane_fault_mms,mms_hydrofrac
    use material, only : block_material
    use utilities, only : search_binary
    
    implicit none

    type(block_grid),intent(in) :: B
    type(grid_type),intent(in) :: G
    type(fields_type),intent(inout) :: F
    type(block_material),intent(in) :: M
    type(source_type),intent(in) :: S
    integer,intent(in) :: mode,iblock
    real,intent(in) :: t

    real,parameter :: pi = 3.14159265358979323846264338327950288419716939937510d0

    real :: x,y
    integer :: i,j

    integer :: i0,j0,mm,pp
    logical :: outx,outy
    real :: hx,hy,ax,ay,dx,dy,wx,wy,A
    real :: weight_old,weight_new
    
    if (B%skip) return ! process has no cells in this block

    if (S%use_gravity) then
       select case(mode)
       case(2)
          do j = B%my,B%py
             do i = B%mx,B%px
                F%DF(i,j,1) = F%DF(i,j,1)+S%gx
                F%DF(i,j,2) = F%DF(i,j,2)+S%gy
             end do
          end do
       case(3)
          do j = B%my,B%py
             do i = B%mx,B%px
                F%DF(i,j,1) = F%DF(i,j,1)+S%gz
             end do
          end do
       end select
    end if

    if (S%use_singular_source) then
       
       mm = max(B%mbx,B%mgx); pp = min(B%pbx,B%pgx); 
       i0 = search_binary(G%x(mm:pp,B%my),S%x0,'floor',outx)+mm-1
       if (outx) return

       mm = max(B%mby,B%mgy); pp = min(B%pby,B%pgy); 
       j0 = search_binary(G%y(B%mx,mm:pp),S%y0,'floor',outy)+mm-1
       if (outy) return

       hx = G%x(B%mx+1,B%my)-G%x(B%mx,B%my)
       hy = G%y(B%mx,B%my+1)-G%y(B%mx,B%my)

       !print *, myid,i0,j0,G%x(i0,j0),G%y(i0,j0)

       ax = (S%x0-G%x(i0,B%my))/hx
       ay = (S%y0-G%y(B%mx,j0))/hy

       select case(S%source_time_function)
       case('ramp')
          A = max(min((t-S%t_center)/S%t_duration,1d0),0d0) ! ramp in time
       case('sinusoid')
          A = sin(2d0*pi*(t-S%t_center)/S%t_duration) ! sinusoid
       case('gaussian')
          A = exp(-0.5d0*(t-S%t_center)**2/S%t_duration**2) ! Gaussian
       case('Ricker')
          A = (1d0-2d0*pi**2*(t-S%t_center)**2/S%t_duration**2)*exp(-pi**2*(t-S%t_center)**2/S%t_duration**2) ! Ricker wavelet
       case default
          call error('Invalid source time function','set_source')
       end select
       
       do j = B%my,B%py
          do i = B%mx,B%px
             
             call singular_source(i,i0,hx,ax,dx,wx)
             call singular_source(j,j0,hy,ay,dy,wy)

             select case(mode)
             case(2)
                ! effective body force approach
                !F%DF(i,j,1) = F%DF(i,j,1)+A*(S%Fx*dx*dy-S%Mxx*wx*dy-S%Mxy*dx*wy)/M%rho
                !F%DF(i,j,2) = F%DF(i,j,2)+A*(S%Fy*dx*dy-S%Mxy*wx*dy-S%Myy*dx*wy)/M%rho
                F%DF(i,j,1) = F%DF(i,j,1)+A*S%Fx*dx*dy/M%rho
                F%DF(i,j,2) = F%DF(i,j,2)+A*S%Fy*dx*dy/M%rho
                ! transformation strain approach
                F%DF(i,j,3) = F%DF(i,j,3)-A*S%Mxx*dx*dy
                F%DF(i,j,4) = F%DF(i,j,4)-A*S%Mxy*dx*dy
                F%DF(i,j,5) = F%DF(i,j,5)-A*S%Myy*dx*dy
             case(3)
                ! effective body force approach
                !F%DF(i,j,1) = F%DF(i,j,1)+A*(S%Fz*dx*dy-S%Mzx*wx*dy-S%Mzy*dx*wy)/M%rho
                F%DF(i,j,1) = F%DF(i,j,1)+A*S%Fz*dx*dy/M%rho
                ! transformation strain approach
                F%DF(i,j,2) = F%DF(i,j,2)-A*S%Mzx*dx*dy
                F%DF(i,j,3) = F%DF(i,j,3)-A*S%Mzy*dx*dy
             end select

          end do
       end do

    end if

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
    case('mms-hydrofrac')
       do j = B%my,B%py
          do i = B%mx,B%px
             x = G%x(i,j)
             y = G%y(i,j)
             F%DF(i,j,1) = F%DF(i,j,1) + mms_hydrofrac(x,y,t,iblock,'s_vx')
             F%DF(i,j,2) = F%DF(i,j,2) + mms_hydrofrac(x,y,t,iblock,'s_vy')
             F%DF(i,j,3) = F%DF(i,j,3) + mms_hydrofrac(x,y,t,iblock,'s_sxx')
             F%DF(i,j,4) = F%DF(i,j,4) + mms_hydrofrac(x,y,t,iblock,'s_sxy')
             F%DF(i,j,5) = F%DF(i,j,5) + mms_hydrofrac(x,y,t,iblock,'s_syy')
             F%DF(i,j,6) = F%DF(i,j,6) + mms_hydrofrac(x,y,t,iblock,'s_szz')
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

    select case(S%problem)
    case('Forcing') ! tsunami forcing
       select case(mode)
       case(2)
       case(3)
          do j = B%my,B%py
             do i = B%mx,B%px
                x = G%x(i,j)
                y = G%y(i,j)
                F%DF(i,j,1) = F%DF(i,j,1) + &            
                     S%gauss_amp*exp(-0.5d0*(x/S%gauss_sigma)**2 - 0.5d0*(y/S%gauss_sigma)**2) * &
                     exp(-0.5d0*((t-4d0*S%gauss_sigmat)/S%gauss_sigmat)**2) / (S%gauss_sigmat*sqrt(2d0*pi))
             end do
          end do
       end select
    end select

    if (S%forcing_from_file) then
       ! linearly interpolate in time
       weight_old = (t-S%t_old)/(S%t_new-S%t_old)
       weight_new = 1d0-weight_old
       F%DF(B%mx:B%px,B%my:B%py,:) = F%DF(B%mx:B%px,B%my:B%py,:) + &
            weight_old*S%forcing_old(B%mx:B%px,B%my:B%py,:) + &
            weight_new*S%forcing_new(B%mx:B%px,B%my:B%py,:)
    end if
       
  end subroutine set_source


  subroutine singular_source(i,i0,h,a,d,w)
    
    implicit none

    integer,intent(in) :: i,i0
    real,intent(in) :: h,a
    real,intent(out) :: d,w

    real :: beta,gam,gap
    
    ! d = delta function
    ! w = derivative of delta function

    beta = (1d0-2d0*a)/4d0 + 0.25d0
    gam = (2d0*a-1d0)/16d0
    gap = (2d0*a+1d0)/16d0

    if     (i==i0-2) then
       d = (-gam)/h
    elseif (i==i0-1) then
       d = (beta+3d0*gam-gap)/h
    elseif (i==i0  ) then
       d = (1d0-a-2d0*beta-3d0*gam+3d0*gap)/h
    elseif (i==i0+1) then
       d = (a+beta+gam-3d0*gap)/h
    elseif (i==i0+2) then
       d = (gap)/h
    else
       d = 0d0
    end if

    w = 0d0
    
!!$    ! two point discretization
!!$    if     (i==i0  ) then
!!$       d = (1d0-a)/h
!!$       w =  1d0/h**2
!!$    elseif (i==i0+1) then
!!$       d =      a /h
!!$       w = -1d0/h**2
!!$    else
!!$       d = 0d0
!!$       w = 0d0
!!$    end if

    ! centered three point discretization, correct for w only with a/=0
    !if     (i==i0-1) then
    !   d = 0d0
    !   w = ( 1d0-2d0*a)/(2d0*h**2)
    !elseif (i==i0  ) then
    !   d = 1d0/h
    !   w = 2d0*a/h**2
    !elseif (i==i0+1) then
    !   d = 0d0
    !   w = (-1d0-2d0*a)/(2d0*h**2)
    !else
    !   d = 0d0
    !   w = 0d0
    !end if


!!$    ! four point discretization (orthogonal to pi mode)
!!$    if     (i==i0-1) then
!!$       d = ( 1d0+2d0*a*(a-2d0))/(8d0*h)
!!$       w = ( 1d0-a)/(2d0*h**2)
!!$    elseif (i==i0  ) then
!!$       d = ( 5d0-2d0*a**2    )/(8d0*h)
!!$       w = (   a)/(2d0*h**2)
!!$    elseif (i==i0+1) then
!!$       d = ( 3d0-2d0*a*(a-2d0))/(8d0*h)
!!$       w = (-1d0+a)/(2d0*h**2)
!!$    elseif (i==i0+2) then
!!$       d = (-1d0+2d0*a**2    )/(8d0*h)
!!$       w = (  -a)/(2d0*h**2)
!!$    else
!!$       d = 0d0
!!$       w = 0d0
!!$    end if
    
  end subroutine singular_source


  subroutine init_forcing_from_file(F,C,S,checkpoint_number)

    use mpi_routines2d, only : cartesian,allocate_array_body
    use mpi_routines, only : pw
    use mpi
    use io, only : read_file_distributed,open_file_distributed
    use fields, only : fields_type
    
    implicit none

    type(fields_type),intent(in) :: F
    type(cartesian),intent(in) :: C
    type(source_type),intent(inout) :: S
    integer,intent(in) :: checkpoint_number
    
    integer(MPI_OFFSET_KIND) :: offset
    integer :: l,ierr

    if (.not.S%forcing_from_file) return
    
    if (.not.allocated(S%forcing_old)) then
       ! allocate arrays
       call allocate_array_body(S%forcing_old,C,F%nF,ghost_nodes=.false.)
       call allocate_array_body(S%forcing_new,C,F%nF,ghost_nodes=.false.)
    end if

    ! read values from file (at initial time, adjusting if restarting from checkpoint
    
    offset = int(checkpoint_number,kind(offset))*int(C%nx,kind(offset))*int(C%ny,kind(offset))*int(pw,kind(offset))
    print *, offset,checkpoint_number
    call open_file_distributed(S%fh,S%forcing_filename,'read',C%c2d%comm,C%c2d%array_w,pw,offset)

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    do l = 1,F%nF
       call read_file_distributed(S%fh,S%forcing_new(C%mx:C%px,C%my:C%py,l))
    end do

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    
  end subroutine init_forcing_from_file


  subroutine destroy_forcing(S)

    use io, only : close_file_distributed
    
    implicit none

    type(source_type),intent(inout) :: S

    if (allocated(S%forcing_old)) deallocate(S%forcing_old)
    if (allocated(S%forcing_new)) deallocate(S%forcing_new)

    call close_file_distributed(S%fh)
    
  end subroutine destroy_forcing

  
  subroutine load_forcing(F,C,S)

    use mpi_routines2d, only : cartesian
    use fields, only : fields_type
    use io, only : read_file_distributed
    use mpi
    
    implicit none

    type(fields_type),intent(in) :: F
    type(cartesian),intent(in) :: C
    type(source_type),intent(inout) :: S

    integer :: l,ierr

    if (.not.S%forcing_from_file) return
    
    ! overwrite old rate with previous new rate
    S%forcing_old = S%forcing_new

    ! read new rate from file

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    do l = 1,F%nF
       call read_file_distributed(S%fh,S%forcing_new(C%mx:C%px,C%my:C%py,l))
    end do

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

  end subroutine load_forcing

  
end module source
