module fd

  implicit none

  type :: limits
     integer :: nb,m,p,mg,pg,mb,pb
     logical :: endm,endp
  end type limits


contains


  subroutine Hnorm(L,H)

    use fd_coeff, only : HL,HR,nbst

    implicit none

    type(limits),intent(in) :: L
    real,dimension(L%m:L%p),intent(out) :: H

    integer :: i

    H = 1d0 ! interior values

    if (L%endm) then
       do i = 0,nbst-1
          H(L%mg+i) = HL(i)
       end do
    end if
    
    if (L%endp) then
       do i = 0,nbst-1
          H(L%pg-i) = HR(-i)
       end do
    end if
    
  end subroutine Hnorm


  subroutine diff(L,f,Df)

    use fd_coeff, only : DI,DL,DR,mI,pI,nbnd,nbst
    use io, only : error

    implicit none

    type(limits),intent(in) :: L
    real,dimension(L%mb:L%pb),intent(in) :: f
    real,dimension(L%m:L%p),intent(out) :: Df

    integer :: i

    ! Df is first derivative of f with unit grid spacing

    ! special case of single point
    
    if (L%m==L%p) then
       Df = 0d0
       return
    end if
       
    if (L%endm) then
       do i = 0,nbst-1
          Df(L%mg+i) = dot_product(DL(:, i),f(L%mg:L%mg+nbnd-1))
       end do
    end if
    
    if (L%endp) then
       do i = 0,nbst-1
          Df(L%pg-i) = dot_product(DR(:,-i),f(L%pg-nbnd+1:L%pg))
       end do
    end if
    
    do i = max(L%m,L%mg+nbst),min(L%p,L%pg-nbst)
       Df(i) = dot_product(DI,f(i+mI:i+pI))
    end do

  end subroutine diff


  subroutine diff_bndm(L,f,Df)

    use fd_coeff, only : DL,nbnd,nbst
    use io, only : error

    implicit none

    type(limits),intent(in) :: L
    real,dimension(L%mg:L%mg+nbnd-1),intent(in) :: f
    real,dimension(L%mg:L%mg+nbst-1),intent(out) :: Df

    integer :: i

    ! Df is first derivative of f with unit grid spacing
       
    if (L%endm) then
       do i = 0,nbst-1
          Df(L%mg+i) = dot_product(DL(:, i),f(L%mg:L%mg+nbnd-1))
       end do
    end if

  end subroutine diff_bndm


  subroutine diff_bndp(L,f,Df)

    use fd_coeff, only : DR,nbnd,nbst
    use io, only : error

    implicit none

    type(limits),intent(in) :: L
    real,dimension(L%pg-nbnd+1:L%pg),intent(in) :: f
    real,dimension(L%pg-nbst+1:L%pg),intent(out) :: Df

    integer :: i

    if (L%endp) then
       do i = 0,nbst-1
          Df(L%pg-i) = dot_product(DR(:,-i),f(L%pg-nbnd+1:L%pg))
       end do
    end if
    
  end subroutine diff_bndp


  function diff_point(f,location,stencil) result(Df)

    use fd_coeff, only : DI,DL,DR

    implicit none

    real,intent(in) :: f(:)
    character(*),intent(in) :: location
    integer,intent(in) :: stencil
    real :: Df

    select case(location)
    case('m')
       Df = dot_product(DL(:,stencil),f)
    case('p')
       Df = dot_product(DR(:,stencil),f)
    case default
       Df = dot_product(DI,f)
    end select

  end function diff_point


  subroutine diss(L,f,Af)

    use fd_coeff, only : AI,AL,AR,mI,pI,nbnd,nbst
    use io, only : error

    implicit none

    type(limits),intent(in) :: L
    real,dimension(L%mb:L%pb),intent(in) :: f
    real,dimension(L%m:L%p),intent(inout) :: Af

    integer :: i

    ! Af is A*f for artificial dissipation operator A with unit grid spacing
    ! NOTE: A*f is added to input Af (no overwriting of input), unlike diff() subroutine

    ! special case of single point
    
    if (L%m==L%p) return
       
    if (L%endm) then
       do i = 0,nbst-1
          Af(L%mg+i) = Af(L%mg+i)+dot_product(AL(:, i),f(L%mg:L%mg+nbnd-1))
       end do
    end if
    
    if (L%endp) then
       do i = 0,nbst-1
          Af(L%pg-i) = Af(L%pg-i)+dot_product(AR(:,-i),f(L%pg-nbnd+1:L%pg))
       end do
    end if
    
    do i = max(L%m,L%mg+nbst),min(L%p,L%pg-nbst)
       Af(i) = Af(i)+dot_product(AI,f(i+mI:i+pI))
    end do

  end subroutine diss


  subroutine reconstruct_flux(L,fp,fm,fhatp,fhatm,FDmethod)

    use fd_coeff, only : reconstruct3,reconstruct5,reconstructW
    use io, only : error

    implicit none

    type(limits),intent(in) :: L
    real,dimension(L%mb:L%pb),intent(in) :: fp,fm
    real,dimension(L%m-1:L%p),intent(out) :: fhatp,fhatm
    character(*),intent(in) :: FDmethod

    integer :: i ! edge index (NOT cell index)

    ! |  .  |  .  |  .  |
    !      j-1 j  j
    ! cell j has left edge j-1 and right edge j

    select case(FDmethod)

    case default

       call error('Invalid FDmethod','reconstruct_flux')

    case('UPW1') ! first order reconstruction
       
       ! plus flux (waves moving to right, upwind to left)
       
       ! plus flux at left edge (i_edge = mG-1) to be set by boundary condition
       if (L%endm) fhatp(L%mg-1) = 1d40
       
       ! cell edge i, plus flux upwinds to left => cell center i
       do i = max(L%m-1,L%mg),L%p ! mg <= i_edge <= pg
          fhatp(i) = fp(i)
       end do
       
       ! minus flux (waves moving to left, upwind to right)
       
       ! cell edge i, minus flux upwinds to right => cell center i+1
       do i = L%m-1,min(L%p,L%pg-1) ! mg-1 <= i_edge <= pg-1
          fhatm(i) = fm(i+1)
       end do
       
       ! minus flux at right edge (i_edge = pG) to be set by boundary condition
       if (L%endp) fhatm(L%pg) = 1d40
       
    case('UPW3') ! third order reconstruction

       ! plus flux (waves moving to right, upwind to left)
       
       if (L%endm) then
          fhatp(L%mg-1) = 1d40 ! set by boundary condition
          fhatp(L%mg  ) = reconstruct3( 0,fp(L%mg:L%mg+2))
       end if

       do i = max(L%m-1,L%mg+1),min(L%p,L%pg-1)
          fhatp(i) = reconstruct3(1,fp(i-1:i+1))
       end do

       if (L%endp) then
          fhatp(L%pg  ) = reconstruct3( 2,fp(L%pg-2:L%pg))
       end if

       ! minus flux (waves moving to left, upwind to right)
       
       if (L%endm) then
          fhatm(L%mg-1) = reconstruct3(-1,fm(L%mg:L%mg+2))
       end if

       do i = max(L%m-1,L%mg),min(L%p,L%pg-2)
          fhatm(i) = reconstruct3(0,fm(i:i+2))
       end do
       
       if (L%endp) then
          fhatm(L%pg-1) = reconstruct3( 1,fm(L%pg-2:L%pg))
          fhatm(L%pg  ) = 1d40 ! set by boundary condition
       end if
       
    case('UPW5') ! fifth-order interior, third-order boundaries

       ! plus flux (waves moving to right, upwind to left)
       
       if (L%endm) then
          fhatp(L%mg-1) = 1d40 ! set by boundary condition
          fhatp(L%mg  ) = reconstruct3( 0,fp(L%mg:L%mg+2))
          fhatp(L%mg+1) = reconstruct3( 1,fp(L%mg:L%mg+2))
       end if

       do i = max(L%m-1,L%mg+2),min(L%p,L%pg-2)
          fhatp(i) = reconstruct5(2,fp(i-2:i+2))
       end do

       if (L%endp) then
          fhatp(L%pg-1) = reconstruct3( 1,fp(L%pg-2:L%pg))
          fhatp(L%pg  ) = reconstruct3( 2,fp(L%pg-2:L%pg))
       end if

       ! minus flux (waves moving to left, upwind to right)
       
       if (L%endm) then
          fhatm(L%mg-1) = reconstruct3(-1,fm(L%mg:L%mg+2))
          fhatm(L%mg  ) = reconstruct3( 0,fm(L%mg:L%mg+2))
       end if

       do i = max(L%m-1,L%mg+1),min(L%p,L%pg-3)
          fhatm(i) = reconstruct5(1,fm(i-1:i+3))
       end do
       
       if (L%endp) then
          fhatm(L%pg-2) = reconstruct3( 0,fm(L%pg-2:L%pg))
          fhatm(L%pg-1) = reconstruct3( 1,fm(L%pg-2:L%pg))
          fhatm(L%pg  ) = 1d40 ! set by boundary condition
       end if

    case('WENO') ! fifth-order WENO interior, third-order boundaries

       ! plus flux (waves moving to right, upwind to left)
       
       if (L%endm) then
          fhatp(L%mg-1) = 1d40 ! set by boundary condition
          fhatp(L%mg  ) = reconstruct3( 0,fp(L%mg:L%mg+2))
          fhatp(L%mg+1) = reconstruct3( 1,fp(L%mg:L%mg+2))
       end if

       do i = max(L%m-1,L%mg+2),min(L%p,L%pg-2)
          fhatp(i) = reconstructW(2,fp(i-2:i+2))
       end do

       if (L%endp) then
          fhatp(L%pg-1) = reconstruct3( 1,fp(L%pg-2:L%pg))
          fhatp(L%pg  ) = reconstruct3( 2,fp(L%pg-2:L%pg))
       end if

       ! minus flux (waves moving to left, upwind to right)
       
       if (L%endm) then
          fhatm(L%mg-1) = reconstruct3(-1,fm(L%mg:L%mg+2))
          fhatm(L%mg  ) = reconstruct3( 0,fm(L%mg:L%mg+2))
       end if

       do i = max(L%m-1,L%mg+1),min(L%p,L%pg-3)
          !fhatm(i) = reconstructW(1,fm(i-1:i+3)) ! incorrect
          fhatm(i) = reconstructW(2,fm(i+3:i-1:-1))
       end do
       
       if (L%endp) then
          fhatm(L%pg-2) = reconstruct3( 0,fm(L%pg-2:L%pg))
          fhatm(L%pg-1) = reconstruct3( 1,fm(L%pg-2:L%pg))
          fhatm(L%pg  ) = 1d40 ! set by boundary condition
       end if

    end select

  end subroutine reconstruct_flux


end module fd
