module utilities

  implicit none

  interface within
     module procedure within_integer,within_real
  end interface


contains


  function search_binary(x,y,mode,out_range,ig,z) result(mid)
    ! SEARCH_BINARY performs a binary search of ordered list x(1:n), 
    ! returning index i of x whose value is closest to y and satisfies either
    ! x(i)<=y<x(i+1) (floor mode) or x(i-1)<y<=x(i) (ceiling mode)
    ! as well as z = x(i).  The values of i must satisfy 1<=i<=n.  If
    ! the value y lies below (or equal to, for ceiling mode) x(1) then 
    ! i=1 is returned.  If the value y lies above (or equal to, for floor
    ! mode) x(n) then i=n is returned.  In these cases, out_range is
    ! set to T to indicate that y lies outside range of x.

    use io, only : error

    implicit none

    ! X = ordered list
    ! Y = value to be found in list
    ! MID = mean of upper and low bounds
    ! MODE = whether x(i) is greater or lesser than y
    !      'floor' = x(i)<=y<x(i+1)
    !      'ceiling' = x(i-1)<y<=x(i)
    ! OUT_RANGE = flag used to signify that y lies outside range of x
    !      T = outside range of x
    !      F = inside range of x
    ! IG = guess at index i
    ! Z = value of x(i) that satisfies inequality

    real,dimension(:),intent(in) :: x
    real,intent(in) :: y
    character(*),intent(in) :: mode
    logical,intent(out) :: out_range
    integer,intent(in),optional :: ig
    real,intent(out),optional :: z
    integer :: mid

    ! N = length of x
    ! LOW = lower bound on index
    ! HIGH = upper bound on index
    ! I = used in expanding binary search
    ! J = used in expanding binary search

    integer :: n,low,high,i,j

    n = size(x)

    ! special case of size(x)==1

    if (n==1) then
       mid = 1
       if (y==x(1)) then
          out_range = .false.
       else
          out_range = .true.
       end if
       if (present(z)) z = x(1)
       return
    end if

    ! use initial guess (if given) to isolate range containing index,
    ! otherwise start from full range

    if (present(ig)) then

       ! make sure initial guess is within bounds

       i = min(max(ig,1),n)
       
       ! determine if y is above or below guess x(i), set low and high

       if (y==x(i)) then
          mid = i
          if (present(z)) z = x(i)
          out_range = .false.
          return
       elseif (y<x(i)) then
          ! y less than x(i) => i is upper bound
          high = i
          ! find lower bound, moving out in powers of 2
          j = 0
          do
             low = max(1,high-2**j)
             if (x(low)<y.or.low==1) exit
             j = j+1
          end do
       else ! y>x(i)
          ! y greater than x(i) => i is lower bound
          low = i
          ! find upper bound, moving out in powers of 2
          j = 0
          do
             high = min(low+2**j,n)
             if (y<x(high).or.high==n) exit
             j = j+1
          end do
       end if

    else

       low = 1
       high = n

    end if

    ! perform binary search to locate index

    select case(mode)

    case('floor')

       if (y<x(low)) then
          ! out of bounds below
          mid = low
          out_range = .true.
       elseif (y>=x(high)) then
          ! out of bounds above
          mid = high-1
          out_range = .true.
       else
          mid = search_binary_floor(x,y,low,high)
          out_range = .false.
       end if

    case('ceiling')

       if (y<=x(low)) then
          ! out of bounds below
          mid = low+1
          out_range = .true.
       elseif (y>x(high)) then
          ! out of bounds above
          mid = high
          out_range = .true.
       else
          mid = search_binary_ceiling(x,y,low,high)
          out_range = .false.
       end if

    case default

       call error('Invalid mode in binary search','search_binary')

    end select

    if (present(z)) z = x(mid)

  end function search_binary


  recursive function search_binary_floor(x,y,low,high) result(mid)
    ! SEARCH_BINARY_FLOOR performs a binary search of ordered list x, 
    ! returning index i of x whose value is closest to y and satisfies x(i)<=y<x(i+1)

    implicit none

    ! X = ordered list
    ! Y = value to be found in list
    ! LOW = lower bound on index
    ! HIGH = upper bound on index
    ! MID = mean of upper and low bounds

    real,dimension(:),intent(in) :: x
    real,intent(in) :: y
    integer,intent(in) :: low,high
    integer :: mid

    ! find mean of low and high, rounding down
    mid = (low+high)/2

    if (x(mid)<=y) then
       ! value greater than or equal to mid
       if (y<x(mid+1)) then
          ! value less than mid+1, properly bracketed
          return
       else
          ! value greater than mid+1, make this new lower bound
          mid = search_binary_floor(x,y,mid+1,high)
       end if
    else
       ! value less than mid, make this new upper bound
       mid = search_binary_floor(x,y,low,mid)
    end if

  end function search_binary_floor


  recursive function search_binary_ceiling(x,y,low,high) result(mid)
    ! SEARCH_BINARY_CEILING performs a binary search of ordered list x, 
    ! returning index i of x whose value is closest to y and satisfies x(i-1)<y<=x(i)

    implicit none

    ! X = ordered list
    ! Y = value to be found in list
    ! LOW = lower bound on index
    ! HIGH = upper bound on index
    ! MID = mean of upper and low bounds

    real,dimension(:),intent(in) :: x
    real,intent(in) :: y
    integer,intent(in) :: low,high
    integer :: mid

    ! find mean of low and high, rounding up
    mid = ceiling(real(low+high)/2d0)

    if (x(mid)>=y) then
       ! value less than mid
       if (y>x(mid-1)) then
          ! value less than or equal to mid-1, properly bracketed
          return
       else
          ! value less than mid-1, make this new upper bound
          mid = search_binary_ceiling(x,y,low,mid-1)
       end if
    else
       ! value greater than or equal to mid, make this new lower bound
       mid = search_binary_ceiling(x,y,mid,high)
    end if

  end function search_binary_ceiling


  function word_count(str) result(nw)
    ! WORD_COUNT counts the number of words (delimited by spaces) in a string

    implicit none
    
    ! I/O Parameters:
    ! STR = string to be processed
    ! NW = number of words in string

    character(*),intent(in) :: str
    integer :: nw
    
    ! Internal Parameters:
    ! I = location of space in string
    ! C = current position in string

    integer :: i,c

    nw = 0
    c = 1
    do
       i = index(trim(str(c:)),' ')
       nw = nw+1
       if (i==0) return
       c = c+i
    end do

  end function word_count


  subroutine extract_words(str,words)
    ! EXTRACT_WORDS extracts space-delimited words from string
    
    implicit none

    ! I/O Parameters:
    ! STR = string to be processed
    ! WORDS = extracted words

    character(*),intent(in) :: str
    character(*),dimension(:),intent(out) :: words
    
    ! Internal Parameters:
    ! I = location of space in string
    ! N = length of trimmed string
    ! C = current position in string
    ! NW = number of words

    integer :: i,n,c,nw

    n = len_trim(str)

    nw = 0
    c = 1
    do
       i = index(trim(str(c:)),' ')
       nw = nw+1
       if (i==0) then
          words(nw) = str(c:n)
          return
       else
          words(nw) = str(c:c+i-2)
          c = c+i
       end if
    end do

  end subroutine extract_words
  

  elemental function step(r,r0,fin,fout) result(f)
    
    implicit none
    
    real,intent(in) :: r,r0,fin,fout
    real :: f

    real,parameter :: tol = 100d0*epsilon(r)

    if (r<r0-tol) then ! inside
       f = fin
    elseif (abs(r-r0)<=tol) then ! boundary
       f = 0.5d0*(fin+fout)
    else ! outside
       f = fout
    end if

  end function step


  elemental function boxcar(r,w,fin,fout) result(f)
    
    implicit none
    
    real,intent(in) :: r,w,fin,fout
    real :: f

    real,parameter :: tol = 100d0*epsilon(r)

    if (-w+tol<r.and.r<w-tol) then ! inside
       f = fin
    elseif (abs(-w-r)<=tol.or.abs(r-w)<=tol) then ! boundary
       f = 0.5d0*(fin+fout)
    else ! outside
       f = fout
    end if

  end function boxcar


  elemental function gaussian(r,w,fin,fout) result(f)
    
    implicit none
    
    real,intent(in) :: r,w,fin,fout
    real :: f

    f = fout+(fin-fout)*exp(-0.5d0*(r/w)**2)

  end function gaussian


  ! transition between f = a1*x+b1 and f = a2*x+b2
  elemental function smooth_wedge(a,b,r,c,m,x) result(f)
    
    implicit none
    
    real,intent(in) :: a,b,r,c,m,x
    real :: f

    ! f = @(a,b,m,r,c,x) (1/2).*x.*(a+b)+(a-b).*log(cosh(m.*(r-x)))./(2.*m)+c;
    f = (1/2.0d0)*x*(a+b)+(a-b)*log(cosh(m*(r-x)))/(2.0d0*m)+c

  end function smooth_wedge

  elemental function smooth(r,w,fin,fout) result(f)
    
    implicit none
    
    real,intent(in) :: r,w,fin,fout
    real :: f

    real,parameter ::  tol = 100d0*epsilon(r)

    if (r>=w-tol) then
       F = fout
    else
       F = fout+(fin-fout)*exp(r**2/(r**2-w**2))
    end if

  end function smooth


  elemental function triangle(r,w,fin,fout) result(f)
    
    implicit none
    
    real,intent(in) :: r,w,fin,fout
    real :: f
    
    if (-w<r.and.r<w) then ! inside
       f = fout+(fin-fout)*(1d0-r/w)
    else ! outside
       f = fout
    end if

  end function triangle



  elemental function smooth_boxcar(r,w,l,fin,fout) result(f)
    
    implicit none
    
    real,intent(in) :: r,w,l,fin,fout
    real :: f

    if (r<=w) then
       f = 1d0
    elseif (r<w+l) then
       f = 0.5d0*(1d0+tanh(l/(r-w-l)+l/(r-w)))
    else
       f = 0d0
    end if

    f = fout+(fin-fout)*f

  end function smooth_boxcar


  elemental function decaying_step(x,x0,w,fin,fout) result(f)
    
    implicit none
    
    real,intent(in) :: x,x0,w,fin,fout
    real :: f

    if (x==x0) then
       f = 0.5d0
    elseif (x<x0) then
       f = exp((x-x0)/w)
    else
       f = 0d0
    end if

    f = fout+(fin-fout)*f

  end function decaying_step


  function within_integer(m,i,p) result(within)

    integer,intent(in) :: m,i,p
    logical :: within

    within = (m<=i.and.i<=p)

  end function within_integer


  function within_real(m,i,p) result(within)

    real,intent(in) :: m,i,p
    logical :: within

    within = (m<=i.and.i<=p)

  end function within_real


  elemental function deg2rad(deg) result(rad)

    implicit none

    real,intent(in) :: deg
    real :: rad

    real,parameter :: pi_over_180 = 0.017453292519943d0

    rad = pi_over_180*deg

  end function deg2rad


  elemental subroutine convert_time(total_sc,hr,mn,sc)

    real,intent(in) :: total_sc
    integer,intent(out) :: hr,mn
    real,intent(out) :: sc

    hr = total_sc/3600
    mn = (total_sc-3600*hr)/60
    sc = total_sc-3600d0*dble(hr)-60d0*dble(mn)

  end subroutine convert_time


end module utilities
