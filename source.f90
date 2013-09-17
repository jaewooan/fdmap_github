module source

  implicit none

  type :: source_type
     logical :: use_gravity,use_singular_source
     real :: gx,gy,gz ! components of gravitational acceleration
     real :: x0,y0,T,Fx,Fy,Fz,Mxx,Mxy,Myy,Mzx,Mzy ! point force and moment tensor
  end type source_type

contains


  subroutine init_source(S,input,echo)

    use mpi_routines, only : is_master
    use io, only : error,write_matlab

    implicit none
    
    type(source_type),intent(out) :: S
    integer,intent(in) :: input,echo

    real :: gx,gy,gz,x0,y0,T,Fx,Fy,Fz,Mxx,Mxy,Myy,Mzx,Mzy
    integer :: stat

    namelist /source_list/ gx,gy,gz,x0,y0,T,Fx,Fy,Fz,Mxx,Mxy,Myy,Mzx,Mzy

    ! defaults

    gx = 0d0
    gy = 0d0
    gz = 0d0

    x0 = 0d0
    y0 = 0d0
    T = 1d0

    Fx = 0d0
    Fy = 0d0
    Fz = 0d0

    Mxx = 0d0
    Mxy = 0d0
    Myy = 0d0
    Mzx = 0d0
    Mzy = 0d0

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
    S%T = T

    S%Fx = Fx
    S%Fy = Fy
    S%Fz = Fz

    S%Mxx = Mxx
    S%Mxy = Mxy
    S%Myy = Myy
    S%Mzx = Mzx
    S%Mzy = Mzy

    ! output source parameters

    if (is_master) then
       call write_matlab(echo,'gx',S%gx,'S')
       call write_matlab(echo,'gy',S%gy,'S')
       call write_matlab(echo,'gz',S%gz,'S')
       call write_matlab(echo,'x0',S%x0,'S')
       call write_matlab(echo,'y0',S%y0,'S')
       call write_matlab(echo,'T',S%T,'S')
       call write_matlab(echo,'Fx',S%Fx,'S')
       call write_matlab(echo,'Fy',S%Fy,'S')
       call write_matlab(echo,'Fz',S%Fz,'S')
       call write_matlab(echo,'Mxx',S%Mxx,'S')
       call write_matlab(echo,'Mxy',S%Mxy,'S')
       call write_matlab(echo,'Myy',S%Myy,'S')
       call write_matlab(echo,'Mzx',S%Mzx,'S')
       call write_matlab(echo,'Mzy',S%Mzy,'S')
    end if

  end subroutine init_source


  subroutine set_source(B,G,F,M,S,t,mode,iblock)
    
    use io, only : error
    use grid, only : grid_type,block_grid
    use fields, only : fields_type
    use mms, only : mms_sin,inplane_fault_mms
    use material, only : block_material
    use utilities, only : search_binary
    
    !USE MPI_ROUTINES, ONLY : MYID

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

    if (B%skip) return ! process has no cells in this block

    if (S%use_gravity) then
       select case(mode)
       case(2)
          do j = B%my,B%py
             do i = B%mx,B%px
                F%DF(:,:,1) = F%DF(:,:,1)+S%gx
                F%DF(:,:,2) = F%DF(:,:,2)+S%gy
             end do
          end do
       case(3)
          do j = B%my,B%py
             do i = B%mx,B%px
                F%DF(:,:,1) = F%DF(:,:,1)+S%gz
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

       !A = min(t/S%T,1d0) ! ramp in time
       !A = sqrt(t)
       A = sin(2d0*pi*t/S%T) ! sinusoid
             
       do j = B%my,B%py
          do i = B%mx,B%px
             
             call singular_source(i,i0,hx,ax,dx,wx)
             call singular_source(j,j0,hy,ay,dy,wy)

             select case(mode)
             case(2)
                F%DF(i,j,1) = F%DF(i,j,1)+A*(S%Fx*dx*dy-S%Mxx*wx*dy-S%Mxy*dx*wy)/M%rho
                F%DF(i,j,2) = F%DF(i,j,2)+A*(S%Fy*dx*dy-S%Mxy*wx*dy-S%Myy*dx*wy)/M%rho
             case(3)
                F%DF(i,j,1) = F%DF(i,j,1)+A*(S%Fz*dx*dy-S%Mzx*wx*dy-S%Mzy*dx*wy)/M%rho
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


  subroutine singular_source(i,i0,h,a,d,w)
    
    implicit none

    integer,intent(in) :: i,i0
    real,intent(in) :: h,a
    real,intent(out) :: d,w

    ! d = delta function
    ! w = derivative of delta function

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
    if     (i==i0-1) then
       d = 0d0
       w = ( 1d0-2d0*a)/(2d0*h**2)
    elseif (i==i0  ) then
       d = 1d0/h
       w = 2d0*a/h**2
    elseif (i==i0+1) then
       d = 0d0
       w = (-1d0-2d0*a)/(2d0*h**2)
    else
       d = 0d0
       w = 0d0
    end if


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


end module source
