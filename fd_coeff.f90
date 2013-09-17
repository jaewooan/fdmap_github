module fd_coeff

  implicit none

  save

  ! cx(j,r) is coefficient of jth point in x-pt stencil with left shift of r

  ! dw combines three 3-pt stencils to get one 5-pt stencil
  ! dw(s,r) corresponds to linear weight for stencil with left-shift of r-s
  ! and with furthest left-shifted stencil having left shift of r

  real,private :: c6(0:5,-1:5),c5(0:4,-1:4),c4(0:3,-1:3),c3(0:2,-3:4), &
       dw(0:2,-1:4),gw(1:2,-3:4),g1(0:2,-3:4),g2(0:2,-3:4)
  !real,parameter :: eps=1d-40
  real,parameter :: eps=1d-6
  real :: H00i,c5p(-3:2),c5m(-2:3),c5u(-3:3), &
       c5pL( 0:6, 0:3),c5mL( 0:6, 0:3),c5uL( 0:6, 0:3), &
       c5pR(-6:0,-3:0),c5mR(-6:0,-3:0),c5uR(-6:0,-3:0)

  ! NI = number of points in interior stencil
  ! MI = minimum index in interior stencil
  ! PI = maximum index in interior stencil
  ! NBND = number of points in boundary stencil
  ! NBST = number of boundary stencils
  ! DI = interior derivative coefficients
  ! DL/DR = left/right boundary derivative coefficients
  ! AI = interior artificial dissipation coefficients
  ! AL/AR = left/right boundary artificial dissipation coefficients
  ! DMI/DPI = minus/plus flux interior derivative coefficients
  ! DML/DPL,DMR/DPR = minus/plus flux left/right derivative coefficients
  ! HL/HR = diagonal norm entries near boundary (normalized so H=1 in interior)

  integer :: nI,mI,pI,nbnd,nbst
  logical :: break = .false.
  real,dimension(:),allocatable :: DI,AI,DmI,DpI,HL,HR
  real,dimension(:,:),allocatable :: DL,DR,AL,AR,DmL,DmR,DpL,DpR

contains


  subroutine init_fd(FDmethod)

    use io, only : error

    implicit none

    character(*),intent(in) :: FDmethod

    integer :: i

    select case(FDmethod)
    case default
       call error('Invalid FD method','init_fd')
    case('WENO')
       call init_weno
    case('SBP2','SBP2_break')
       nI = 3
       nbnd = 2
       nbst = 1
       H00i = 2d0
    case('SBP4','SBP4_break')
       nI = 5
       nbnd = 6
       nbst = 4
       H00i = 2.823529411764706d0
    case('SBP6','SBP6_break')
       nI = 7
       nbnd = 9
       nbst = 6
       H00i = 3.165067037878233d0
    end select

    mI = -(nI-1)/2
    pI =  (nI-1)/2

    allocate(DI(mI:pI),AI(mI:pI))
    allocate(DL (0:nbnd-1,0:nbst-1),DR (-(nbnd-1):0,-(nbst-1):0))
    allocate(AL (0:nbnd-1,0:nbst-1),AR (-(nbnd-1):0,-(nbst-1):0))
    allocate(HL (         0:nbst-1),HR (            -(nbst-1):0))

    select case(FDmethod)

    case('SBP2','SBP2_break')
       
       ! interior

       DI = (/ -0.5d0, 0d0, 0.5d0 /)
       AI = (/  0.5d0,-1d0, 0.5d0 /)

       ! left boundary

       DL(:,0) = (/ -1d0, 1d0 /)
       AL(:,0) = (/ -1d0, 1d0 /)

       HL = (/ 0.5d0 /)

    case('SBP4','SBP4_break')
       
       ! interior

       DI = (/  0.083333333333333d0, -0.666666666666667d0,  0.000000000000000d0, &
            0.666666666666667d0, -0.083333333333333d0 /)
       AI = (/ -0.083333333333333d0,  0.333333333333333d0, -0.500000000000000d0, &
            0.333333333333333d0, -0.083333333333333d0 /)

       ! left boundary

       DL(:,0) = (/ -1.411764705882353d0,  1.735294117647059d0, -0.235294117647059d0, &
            -0.088235294117647d0,  0.000000000000000d0,  0.000000000000000d0 /)
       DL(:,1) = (/ -0.500000000000000d0,  0.000000000000000d0,  0.500000000000000d0, &
            0.000000000000000d0,  0.000000000000000d0,  0.000000000000000d0 /)
       DL(:,2) = (/ 0.093023255813953d0, -0.686046511627907d0,  0.000000000000000d0, &
            0.686046511627907d0, -0.093023255813953d0,  0.000000000000000d0 /)
       DL(:,3) = (/ 0.030612244897959d0,  0.000000000000000d0, -0.602040816326531d0, &
            0.000000000000000d0,  0.653061224489796d0, -0.081632653061224d0 /)

       AL(:,0) = (/ -0.470588235294118d0, 0.941176470588235d0, -0.470588235294118d0, &
            0.000000000000000d0,  0.000000000000000d0,  0.000000000000000d0 /)
       AL(:,1) = (/ 0.271186440677966d0, -0.610169491525424d0,  0.406779661016949d0, &
            -0.067796610169492d0, 0.000000000000000d0,  0.000000000000000d0 /)
       AL(:,2) = (/ -0.186046511627907d0, 0.558139534883721d0, -0.651162790697674d0, &
            0.372093023255814d0, -0.093023255813953d0,  0.000000000000000d0 /)
       AL(:,3) = (/ 0.000000000000000d0, -0.081632653061224d0,  0.326530612244898d0, &
            -0.489795918367347d0, 0.326530612244898d0, -0.081632653061224d0 /)

       HL = (/ 0.354166666666667d0,1.229166666666667d0,0.895833333333333d0,1.020833333333333d0 /)

    case('SBP6','SBP6_break')

      ! interior

      DI = (/ -0.016666666666667d0,  0.150000000000000d0, -0.750000000000000d0,  &
        0.000000000000000d0,  0.750000000000000d0, -0.150000000000000d0, 0.016666666666667d0 /)
      AI = (/ 0.016666666666667d0,  -0.100000000000000d0,  0.250000000000000d0,  &
        -0.333333333333333d0, 0.250000000000000d0, -0.100000000000000d0, 0.016666666666667d0 /)

      ! left boundary

      DL(:,0) = (/ -1.582533518939116d0,  1.905066305223826d0,  0.371736635162527d0, &
        -1.220272547439373d0,  0.617737563191443d0, -0.091734437199306d0, &
        0.000000000000000d0,  0.000000000000000d0,  0.000000000000000d0 /)
      DL(:,1) = (/ -0.432901856322317d0,  0.000000000000000d0, -0.004314770110158d0, &
        0.841962873553650d0, -0.506472155165238d0,  0.101725908044063d0, &
        0.000000000000000d0,  0.000000000000000d0,  0.000000000000000d0 /)
      DL(:,2) = (/ -0.187157260543465d0,  0.009559818025329d0,  0.000000000000000d0, &
        -0.685786302717324d0,  1.269119636050658d0, -0.405735890815197d0, &
        0.000000000000000d0,  0.000000000000000d0,  0.000000000000000d0 /)
      DL(:,3) = (/ 0.310794924426199d0, -0.943692853144243d0,  0.346924177396280d0, &
        0.000000000000000d0,  0.193459600671767d0,  0.079078808235367d0, &
        0.013435342414630d0,  0.000000000000000d0,  0.000000000000000d0 /)
      DL(:,4) = (/ -0.214078964072617d0,  0.772407007744065d0, -0.873577080952986d0, &
        -0.263234734035800d0,  0.000000000000000d0,  0.724732343108629d0, &
        -0.164529643265203d0,  0.018281071473911d0,  0.000000000000000d0 /)
      DL(:,5) = (/ 0.028585724831244d0, -0.139498337176472d0,  0.251124403552430d0, &
        -0.096751976743301d0, -0.651665106580520d0,  0.000000000000000d0, &
        0.739709139060752d0, -0.147941827812150d0,  0.016437980868017d0 /)

      AL(:,0) = (/ -0.105502234595941d0, 0.316506703787823d0,   -0.316506703787823d0, &
        0.105502234595941d0,   0.000000000000000d0,  0.000000000000000d0, &
        0.000000000000000d0,   0.000000000000000d0,  0.000000000000000d0 /)
      AL(:,1) = (/ 0.071922084408557d0, -0.227753267293765d0,  0.251727295429951d0, &
        -0.107883126612836d0,  0.011987014068093d0,  0.000000000000000d0, &
        0.000000000000000d0,   0.000000000000000d0,  0.000000000000000d0 /)
      AL(:,2) = (/ -0.159350793065290d0, 0.557727775728513d0, -0.743637034304685d0, &
        0.478052379195869d0,  -0.159350793065290d0,  0.026558465510882d0, &
        0.000000000000000d0,   0.000000000000000d0,  0.000000000000000d0 /)
      AL(:,3) = (/ 0.026870684829259d0, -0.120918081731666d0,  0.241836163463333d0, &
        -0.282142190707221d0,  0.201530136219444d0, -0.080612054487778d0, &
        0.013435342414630d0,   0.000000000000000d0,  0.000000000000000d0 /)
      AL(:,4) = (/ 0.000000000000000d0,  0.018281071473911d0, -0.109686428843468d0, &
        0.274216072108671d0,  -0.365621429478228d0,  0.274216072108671d0, &
        -0.109686428843468d0,  0.018281071473911d0,  0.000000000000000d0 /)
      AL(:,5) = (/ 0.000000000000000d0,  0.000000000000000d0,  0.016437980868017d0, &
        -0.098627885208100d0,  0.246569713020251d0, -0.328759617360334d0, &
        0.246569713020251d0,  -0.098627885208100d0,  0.016437980868017d0 /)

      HL = (/ 0.315949074074074d0,1.390393518518518d0,0.627546296296296d0, &
        1.240509259259259d0,0.911689814814815d0,1.013912037037037d0 /)

      ! DI = (/-1d0/60d0,3d0/20d0,-3d0/4d0,0d0,3d0/4d0,-3d0/20d0,1d0/60d0/)

      ! DL(:,0) = (/-21600d0/13649d0,104009d0/54596d0,30443d0/81894d0,&
      !   -33311d0/27298d0,16863d0/27298d0,-15025d0/163788d0,0d0,0d0,0d0/)
      ! DL(:,1) = (/-104009d0/240260d0,0d0,-311d0/72078d0,20229d0/24026d0,&
      !   -24337d0/48052d0,36661d0/360390d0,0d0,0d0,0d0/)
      ! DL(:,2) = (/-30443d0/162660d0,311d0/32532d0,0d0,-11155d0/16266d0,&
      !   41287d0/32532d0,-21999d0/54220d0,0d0,0d0,0d0/)
      ! DL(:,3) = (/33311d0/107180d0,-20229d0/21436d0,485d0/1398d0,0d0,4147d0/21436d0,&
      !   25427d0/321540d0,72d0/5359d0,0d0,0d0/)
      ! DL(:,4) = (/-16863d0/78770d0,24337d0/31508d0,-41287d0/47262d0,-4147d0/15754d0,&
      !   0d0,342523d0/472620d0,-1296d0/7877d0,144d0/7877d0,0d0/)
      ! DL(:,5) = (/15025d0/525612d0,-36661d0/262806d0,21999d0/87602d0,-25427d0/262806d0,&
      !   -342523d0/525612d0,0d0,32400d0/43801d0,-6480d0/43801d0,720d0/43801d0/)

      ! HL = (/13649d0/43200d0,12013d0/8640d0,2711d0/4320d0,&
      !   5359d0/4320d0,7877d0/8640d0,43801d0/43200d0/)

    end select

    ! right boundary, set using symmetry conditions

    do i = 0,nbst-1
      DR (:,-i) = -DL(nbnd-1:0:-1,i)
       AR (:,-i) =  AL(nbnd-1:0:-1,i)
       HR (  -i) =  HL(            i)
    end do

    select case(FDmethod)
    case('SBP2_break','SBP4_break','SBP6_break')
      break = .true.
    case default
      break = .false.
    end select
    
  end subroutine init_fd


  subroutine check_fd

    implicit none

    integer :: i
    real :: zero,one

    ! check consistency of FD coefficients

    ! interior

    one = dot_product(DI,dble( (/ (i, i=mI,pI) /) ))
    print *, 'DI:',one
    zero = sum(AI)
    print *, 'AI:',zero

    print *

    ! left boundary

    do i = 0,nbst-1
       one = dot_product(DL(:,i),dble( (/ (i, i=1,nbnd) /) ))
       print *, 'DL:',one
       zero = sum(AL(:,i))
       print *, 'AL:',zero
    end do

  end subroutine check_fd


  subroutine init_weno

    implicit none

    integer :: r,j

    ! calculate reconstruction coefficients, c 

    do r = 5,-1,-1
       do j = 0,5
          c6(j,r) = coeff_c(r,j,6)
       end do
    end do

    do r = 4,-1,-1
       do j = 0,4
          c5(j,r) = coeff_c(r,j,5)
       end do
    end do

    do r = 3,-1,-1
       do j = 0,3
          c4(j,r) = coeff_c(r,j,4)
       end do
    end do

    do r = 4,-3,-1
       do j = 0,2
          c3(j,r) = coeff_c(r,j,3)
       end do
    end do

    ! calculate linear weights, d

    do r = 4,-1,-1
       dw(0,r) = coeff_c(r,0,5)/coeff_c(r  ,0,3)
       dw(2,r) = coeff_c(r,4,5)/coeff_c(r-2,2,3)
       dw(1,r) = 1d0-dw(0,r)-dw(2,r)
    end do

    ! calculate weights for smoothness indicators

    gw = 1d40
    g1 = 1d40
    g2 = 1d40

    do r = 2,0,-1
       gw(:,r) = (/ 1d0/4d0, 13d0/12d0 /)
       g2(:,r) = (/ 1d0, -2d0, 1d0 /)
    end do
    
    g1(:, 2) = (/ 1d0, -4d0, 3d0 /)
    g1(:, 1) = (/ 1d0, 0d0,-1d0 /)
    g1(:, 0) = (/ 3d0, -4d0, 1d0 /)

  end subroutine init_weno

  
  subroutine destroy_fd

    implicit none

    deallocate(DI,AI)
    deallocate(DL,DR)
    deallocate(AL,AR)
    deallocate(HL,HR)

  end subroutine destroy_fd


  function coeff_c(r,j,k) result(c)

    implicit none

    integer,intent(in) :: r,j,k
    real :: c

    integer :: m,l,q
    real :: A,B,D

    c = 0d0
    
    do m = j+1,k
       
       ! compute A
       
       A = 0d0
       do l = 0,k

          if (l==m) cycle
          
          ! compute B

          B = 1d0
          do q = 0,k
             if (q==m.or.q==l) cycle
             B = B*dble(r-q+1)
          end do
          
          A = A+B

       end do
       
       ! compute D

       D = 1d0
       do l = 0,k
          if (l==m) cycle
          D = D*dble(m-l)
       end do
       
       c = c+A/D

    end do
    
  end function coeff_c


  function reconstruct6(r,f) result(fweno)

    implicit none

    integer,intent(in) :: r
    real,intent(in) :: f(-r:-r+5)
    real :: fweno

    fweno = dot_product(c6(:,r  ),f)

  end function reconstruct6


  function reconstruct5(r,f) result(fweno)

    implicit none

    integer,intent(in) :: r
    real,intent(in) :: f(-r:-r+4)
    real :: fweno

    fweno = dot_product(c5(:,r  ),f)

  end function reconstruct5


  function reconstruct4(r,f) result(fweno)

    implicit none

    integer,intent(in) :: r
    real,intent(in) :: f(-r:-r+3)
    real :: fweno

    fweno = dot_product(c4(:,r  ),f)

  end function reconstruct4


  function reconstruct3(r,f) result(fweno)

    implicit none

    integer,intent(in) :: r
    real,intent(in) :: f(-r:-r+2)
    real :: fweno

    fweno = dot_product(c3(:,r  ),f)

  end function reconstruct3


  function reconstructW(r,f) result(fweno)

    implicit none

    integer,intent(in) :: r
    real,intent(in) :: f(-r:-r+4)
    real :: fweno

    real,dimension(0:2) :: fe,B,a,w
    !real :: tau5

    ! evaluate alternative reconstructions

    fe(0) = dot_product(c3(:,r  ),f(-r  :-r+2))
    fe(1) = dot_product(c3(:,r-1),f(-r+1:-r+3))
    fe(2) = dot_product(c3(:,r-2),f(-r+2:-r+4))

    !fweno = dot_product(fe(:),dw(:,r)) ! linear 5th order
    !return

    ! evaluate smoothness indicators

    ! Jiang and Shu (1996) smoothness indicators
    B(0) = gw(1,r  )*dot_product(g1(:,r  ),f(-r  :-r+2))**2+ &
         gw(2,r  )*dot_product(g2(:,r  ),f(-r  :-r+2))**2
    B(1) = gw(1,r-1)*dot_product(g1(:,r-1),f(-r+1:-r+3))**2+ &
         gw(2,r-1)*dot_product(g2(:,r-1),f(-r+1:-r+3))**2
    B(2) = gw(1,r-2)*dot_product(g1(:,r-2),f(-r+2:-r+4))**2+ &
         gw(2,r-2)*dot_product(g2(:,r-2),f(-r+2:-r+4))**2
    
    ! Zhang and Shu (2006) smoothness indicators
    !B(0) = dot_product(g2(:,r  ),f(-r  :-r+2))**2
    !B(1) = dot_product(g2(:,r-1),f(-r+1:-r+3))**2
    !B(2) = dot_product(g2(:,r-2),f(-r+2:-r+4))**2

    ! evaluate weights

    ! Jiang and Shu (1996)
    a = dw(:,r)/(B+eps)**2

    ! WENO-Z, Borges et al. (2008)
    !tau5 = abs(B(0)-B(2))
    !a = dw(:,r)*(1d0+(tau5/(B+eps))**2)

    w = a/sum(a)
    
    ! combine stencils for WENO approximation

    fweno = dot_product(w,fe)

  end function reconstructW


end module fd_coeff
