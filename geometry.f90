module geometry

  implicit none

  type :: point
     real :: x,y
  end type point

  type :: curve
     real,dimension(:),allocatable :: x,y,xt,yt
     real,dimension(:,:),allocatable :: n
  end type curve

  interface rotate_xy2nt
     module procedure rotate_xy2nt_vec,rotate_xy2nt_mat
  end interface

  interface rotate_nt2xy
     module procedure rotate_nt2xy_vec,rotate_nt2xy_mat
  end interface


contains


  subroutine destroy_curve(C)

    type(curve),intent(inout) :: C

    if (allocated(C%x)) deallocate(C%x)
    if (allocated(C%y)) deallocate(C%y)
    if (allocated(C%xt)) deallocate(C%xt)
    if (allocated(C%yt)) deallocate(C%yt)
    if (allocated(C%n)) deallocate(C%n)

  end subroutine destroy_curve


  subroutine line(C,P1,P2,m,p,t,t1,t2)

    implicit none

    type(curve),intent(inout) :: C
    type(point),intent(in) :: P1,P2
    integer,intent(in) :: m,p
    real,intent(in) :: t(m:p),t1,t2
    
    real :: L,w(m:p),wt(m:p)

    if (m>p) return

    w = (t-t1)/(t2-t1) ! 0<=w<=1
    wt = 1d0/(t2-t1)

    C%x(m:p) = P1%x+(P2%x-P1%x)*w
    C%y(m:p) = P1%y+(P2%y-P1%y)*w

    C%xt(m:p) = (P2%x-P1%x)*wt
    C%yt(m:p) = (P2%y-P1%y)*wt

    L = sqrt((P2%x-P1%x)**2+(P2%y-P1%y)**2)
    if (L/=0) then
       C%n(m:p,1) = (P2%y-P1%y)/L
       C%n(m:p,2) = (P1%x-P2%x)/L
    end if

  end subroutine line


  subroutine line_map(C,P1,P2,m,p,t,t1,t2,B)

    implicit none

    type(curve),intent(inout) :: C
    type(point),intent(in) :: P1,P2
    integer,intent(in) :: m,p
    real,intent(in) :: t(m:p),t1,t2,B
    
    real :: L,w(m:p),wt(m:p)

    if (m>p) return

    w  = (t-t1)/(t2-t1) ! 0<=w<=1
    wt = 1d0/(t2-t1)

    if     (B>0d0) then ! concentrate points around P1
       w  = 1d0+tanh(B*(w-1d0))/tanh(B)
       wt = wt*B/tanh(B)/cosh(B*(w-1d0))**2 ! NOT CHECKED!
    elseif (B<0d0) then ! concentrate points around P2
       w  = tanh(B*w)/tanh(B)
       wt = wt*B/tanh(B)/cosh(B*w)**2 ! NOT CHECKED!
    end if

    C%x(m:p) = P1%x+(P2%x-P1%x)*w
    C%y(m:p) = P1%y+(P2%y-P1%y)*w

    C%xt(m:p) = (P2%x-P1%x)*wt
    C%yt(m:p) = (P2%y-P1%y)*wt

    L = sqrt((P2%x-P1%x)**2+(P2%y-P1%y)**2)
    if (L/=0) then
       C%n(m:p,1) = (P2%y-P1%y)/L
       C%n(m:p,2) = (P1%x-P2%x)/L
    end if

  end subroutine line_map


  subroutine read_curve(C,name,m,p)

    use io, only : new_io_unit,error

    implicit none

    type(curve),intent(inout) :: C
    character(*),intent(in) :: name
    integer,intent(in) :: m,p

    integer :: i,fid,stat
    real :: x,y,nx,ny

    if (m>p) return

    fid = new_io_unit()
    
    open(fid,file=name,iostat=stat,status='old')
    if (stat/=0) call error('Error opening curve file ' // trim(name),'read_curve')

    do

       read(fid,*,end=100) i,x,y,nx,ny

       if     (i<m) then
          cycle
       elseif (i<=p) then
          C%x (i) = x
          C%y (i) = y
          C%n(i,1) = nx
          C%n(i,2) = ny
          if (i==p) then
             close(fid)
             return
          end if
       end if

    end do

100 call error('Index of curve not in file ' //  trim(name),'read_curve')

  end subroutine read_curve


  subroutine transfinite_interp(mx,px,my,py,qmin,qmax,rmin,rmax, &
       q,r,bndL,bndR,bndB,bndT,LB,LT,RB,RT,x,y)

    implicit none

    integer,intent(in) :: mx,px,my,py
    real,intent(in) :: qmin,qmax,rmin,rmax,q(mx:px),r(my:py)
    type(point),intent(in) :: LB,LT,RB,RT
    type(curve),intent(in) :: bndL,bndR,bndB,bndT
    real,dimension(mx:,my:),intent(out) :: x,y

    integer :: i,j
    real :: Wq,Wr

    ! transfinite interpolation (see Boyd, Chebyshev and Fourier Spectral
    ! Methods, second edition, 2001, pp. 114-115)

    Wq = qmax-qmin
    Wr = rmax-rmin

    ! prevent division by 0 in 1D cases
    if (Wq==0d0) Wq = 1d0
    if (Wr==0d0) Wr = 1d0

    do j = my,py
       do i = mx,px

          x(i,j) = &
               (qmax-q(i))*bndL%x(j)/Wq+(q(i)-qmin)*bndR%x(j)/Wq+ &
               (rmax-r(j))*bndB%x(i)/Wr+(r(j)-rmin)*bndT%x(i)/Wr- &
               (qmax-q(i))*(rmax-r(j))*LB%x/(Wq*Wr)- &
               (qmax-q(i))*(r(j)-rmin)*LT%x/(Wq*Wr)- &
               (q(i)-qmin)*(rmax-r(j))*RB%x/(Wq*Wr)- &
               (q(i)-qmin)*(r(j)-rmin)*RT%x/(Wq*Wr)

          y(i,j) = &
               (qmax-q(i))*bndL%y(j)/Wq+(q(i)-qmin)*bndR%y(j)/Wq+ &
               (rmax-r(j))*bndB%y(i)/Wr+(r(j)-rmin)*bndT%y(i)/Wr- &
               (qmax-q(i))*(rmax-r(j))*LB%y/(Wq*Wr)- &
               (qmax-q(i))*(r(j)-rmin)*LT%y/(Wq*Wr)- &
               (q(i)-qmin)*(rmax-r(j))*RB%y/(Wq*Wr)- &
               (q(i)-qmin)*(r(j)-rmin)*RT%y/(Wq*Wr)

       end do
    end do


  end subroutine transfinite_interp


  subroutine transfinite_interp_derivative(mx,px,my,py,qmin,qmax,rmin,rmax, &
       q,r,bndL,bndR,bndB,bndT,LB,LT,RB,RT,xq,xr,yq,yr)

    implicit none

    integer,intent(in) :: mx,px,my,py
    real,intent(in) :: qmin,qmax,rmin,rmax,q(mx:px),r(my:py)
    type(point),intent(in) :: LB,LT,RB,RT
    type(curve),intent(in) :: bndL,bndR,bndB,bndT
    real,dimension(mx:,my:),intent(out) :: xq,xr,yq,yr

    integer :: i,j
    real :: Wq,Wr

    ! transfinite interpolation (see Boyd, Chebyshev and Fourier Spectral
    ! Methods, second edition, 2001, pp. 114-115)

    Wq = qmax-qmin
    Wr = rmax-rmin

    ! prevent division by 0 in 1D cases
    if (Wq==0d0) Wq = 1d0
    if (Wr==0d0) Wr = 1d0

    do j = my,py
      do i = mx,px

        !x(i,j) = &
        !  (qmax-q(i))*bndL%x(j)/Wq+(q(i)-qmin)*bndR%x(j)/Wq &
        !  +(rmax-r(j))*bndB%x(i)/Wr+(r(j)-rmin)*bndT%x(i)/Wr &
        !  -(qmax-q(i))*(rmax-r(j))*LB%x/(Wq*Wr) &
        !  -(qmax-q(i))*(r(j)-rmin)*LT%x/(Wq*Wr) &
        !  -(q(i)-qmin)*(rmax-r(j))*RB%x/(Wq*Wr) &
        !  -(q(i)-qmin)*(r(j)-rmin)*RT%x/(Wq*Wr)

        xq(i,j) = &
          -bndL%x(j)/Wq+bndR%x(j)/Wq &
          +(rmax-r(j))*bndB%xt(i)/Wr+(r(j)-rmin)*bndT%xt(i)/Wr &
          +(rmax-r(j))*LB%x/(Wq*Wr) &
          +(r(j)-rmin)*LT%x/(Wq*Wr) &
          -(rmax-r(j))*RB%x/(Wq*Wr) &
          -(r(j)-rmin)*RT%x/(Wq*Wr)

        xr(i,j) = &
          (qmax-q(i))*bndL%xt(j)/Wq+(q(i)-qmin)*bndR%xt(j)/Wq &
          -bndB%x(i)/Wr+bndT%x(i)/Wr &
          +(qmax-q(i))*LB%x/(Wq*Wr) &
          -(qmax-q(i))*LT%x/(Wq*Wr) &
          +(q(i)-qmin)*RB%x/(Wq*Wr) &
          -(q(i)-qmin)*RT%x/(Wq*Wr)

        !y(i,j) = &
        !  (qmax-q(i))*bndL%y(j)/Wq+(q(i)-qmin)*bndR%y(j)/Wq &
        !  +(rmax-r(j))*bndB%y(i)/Wr+(r(j)-rmin)*bndT%y(i)/Wr &
        !  -(qmax-q(i))*(rmax-r(j))*LB%y/(Wq*Wr) &
        !  -(qmax-q(i))*(r(j)-rmin)*LT%y/(Wq*Wr) &
        !  -(q(i)-qmin)*(rmax-r(j))*RB%y/(Wq*Wr) &
        !  -(q(i)-qmin)*(r(j)-rmin)*RT%y/(Wq*Wr)

        yq(i,j) = &
          -bndL%y(j)/Wq+bndR%y(j)/Wq &
          +(rmax-r(j))*bndB%yt(i)/Wr+(r(j)-rmin)*bndT%yt(i)/Wr &
          +(rmax-r(j))*LB%y/(Wq*Wr) &
          +(r(j)-rmin)*LT%y/(Wq*Wr) &
          -(rmax-r(j))*RB%y/(Wq*Wr) &
          -(r(j)-rmin)*RT%y/(Wq*Wr)

        yr(i,j) = &
          (qmax-q(i))*bndL%yt(j)/Wq+(q(i)-qmin)*bndR%yt(j)/Wq &
          -bndB%y(i)/Wr+bndT%y(i)/Wr &
          +(qmax-q(i))*LB%y/(Wq*Wr) &
          -(qmax-q(i))*LT%y/(Wq*Wr) &
          +(q(i)-qmin)*RB%y/(Wq*Wr) &
          -(q(i)-qmin)*RT%y/(Wq*Wr)

      end do
    end do


  end subroutine transfinite_interp_derivative


  subroutine rotate_xy2nt_vec(wx,wy,wt,wn,n)

    implicit none

    real,intent(in) :: wx,wy,n(2)
    real,intent(out) :: wt,wn

    wt =  n(2)*wx-n(1)*wy
    wn =  n(1)*wx+n(2)*wy

  end subroutine rotate_xy2nt_vec


  subroutine rotate_nt2xy_vec(wx,wy,wt,wn,n)

    implicit none

    real,intent(in) :: wt,wn,n(2)
    real,intent(out) :: wx,wy

    wx =  n(2)*wt+n(1)*wn
    wy = -n(1)*wt+n(2)*wn

  end subroutine rotate_nt2xy_vec


  subroutine rotate_xy2nt_mat(wxx,wxy,wyy,wtt,wnt,wnn,n)

    implicit none

    real,intent(in) :: wxx,wxy,wyy,n(2)
    real,intent(out) :: wtt,wnt,wnn

    wtt = -2d0*n(1)*n(2)*wxy+n(2)**2*wxx+n(1)**2*wyy
    wnt = n(1)*n(2)*(wxx-wyy)+(n(2)**2-n(1)**2)*wxy
    wnn = 2d0*n(1)*n(2)*wxy+n(1)**2*wxx+n(2)**2*wyy

  end subroutine rotate_xy2nt_mat


  subroutine rotate_nt2xy_mat(wxx,wxy,wyy,wtt,wnt,wnn,n)

    implicit none

    real,intent(in) :: wtt,wnt,wnn,n(2)
    real,intent(out) :: wxx,wxy,wyy

    wxx = 2d0*n(1)*n(2)*wnt+n(2)**2*wtt+n(1)**2*wnn
    wxy = -n(1)*n(2)*(wtt-wnn)+(n(2)**2-n(1)**2)*wnt
    wyy = -2d0*n(1)*n(2)*wnt+n(1)**2*wtt+n(2)**2*wnn

  end subroutine rotate_nt2xy_mat


end module geometry
