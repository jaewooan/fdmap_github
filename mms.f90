module mms

  implicit none

  real,save :: besselA ! not a good way to store this


contains

  

  function bessel(x,y,t,side,field) result(F)

    !use ifport ! for Intel Fortran compiler

    implicit none

    real,intent(in) :: x,y,t
    integer,intent(in) :: side
    character(*),intent(in) :: field
    real :: F

    real :: r,J,dJdr,drdx,drdy,F0
    real,parameter :: vz0=1d0,syz0=20d0

    r = sqrt(x**2+y**2)

    J = dbesj0(r)
    dJdr = -dbesj1(r)

    if (r<epsilon(r)) then

       select case(field)
       case('vz')
          F0 = J*cos(t)
          F = (vz0+F0)-0.5d0*besselA*F0
          if (side==1) F = -F
       case('V')
          F0 = J*cos(t)
          F = 2d0*(vz0+F0)-besselA*F0
       case('sxz')
          F = 0d0
       case('syz','S')
          F = syz0
       end select

    else

       select case(field)
       case('vz')
          F0 = J*cos(t)
          F = (vz0+F0)-0.5d0*besselA*F0
          if (side==1) F = -F
       case('V')
          F0 = J*cos(t)
          F = 2d0*(vz0+F0)-besselA*F0
       case('sxz')
          drdx = x/r
          F0 = drdx*dJdr*sin(t)
          F = (1d0-0.5d0*besselA)*F0
          if (side==1) F = -F
       case('syz')
          drdy = y/r
          F0 = drdy*dJdr*sin(t)
          F = (1d0-0.5d0*besselA)*F0
          if (side==1) F = -F
          F = F+syz0
       case('S')
          F = syz0
       end select

    end if

  end function bessel


  function inplane_bessel(x,y,t,field) result(F)

    !use ifport ! for Intel Fortran compiler
    use material, only : block_material
    use io, only : error

    implicit none

    real,intent(in) :: x,y,t
    character(*),intent(in) :: field
    real :: F

    real :: r,r_x,r_y,r_xy
    real :: rp,rp_x,rp_y
    real :: J0p,J1p,J2p
    real :: ftp,ftp_t
    real :: kp,kcp
    real :: rs,rs_x,rs_y
    real :: J0s,J1s,J2s
    real :: fts,fts_t
    real :: ks,kcs
    real :: ux,uy,ux_x,ux_y,uy_x,uy_y
    real :: pi

    real,parameter :: rho=1d0,cs=1d0,cp=1.732050807568877d0,lambda=1d0,G=1d0

    pi= 4d0 * datan(1d0)
    kp = 20d0*pi
    ks = 10d0*pi

    kcp = kp/cp
    kcs = ks/cs

    F = 1d40

    r = sqrt(x**2+y**2)
    if (r<epsilon(r)) then
      r = 0d0

      r_x = 0d0
      r_y = 0d0
    else
      r_x = x / r
      r_y = y / r
    end if

    ! P-wave
    rp     = kcp*r
    rp_x   = kcp*r_x
    rp_y   = kcp*r_y

    J0p   = dbesjn(0,rp)
    J1p   = dbesjn(1,rp)
    J2p   = dbesjn(2,rp)
    ftp    = sin(kp*t)
    ftp_t  = kp*cos(kp*t)

    ! S-wave
    rs     = kcs*r
    rs_x   = kcs*r_x
    rs_y   = kcs*r_y

    J0s   = dbesjn(0,rs)
    J1s   = dbesjn(1,rs)
    J2s   = dbesjn(2,rs)
    fts    = sin(ks*t)
    fts_t  = ks*cos(ks*t)

    select case(field)
    case('vx')
      F = -ftp_t * rp_x * J1p - fts_t * rs_y * J1s
    case('vy')
      F = -ftp_t * rp_y * J1p + fts_t * rs_x * J1s
    case('sxx','syy')
      ux_x = -ftp*0.5d0*kcp**2
      uy_y = -ftp*0.5d0*kcp**2
      if (r>epsilon(r)) then
        ux_x = ftp*(kcp*(x**2-y**2)*J1p - kcp**2*x**2*r*J0p)/r**3 + fts*kcs**2*x*y*J2s/r**2
        uy_y = ftp*(kcp*(y**2-x**2)*J1p - kcp**2*y**2*r*J0p)/r**3 - fts*kcs**2*x*y*J2s/r**2
      end if
      select case(field)
      case('sxx')
        F = (lambda+2d0*G)*ux_x + lambda*uy_y
      case('syy')
        F = lambda*ux_x + (lambda+2d0*G)*uy_y
      end select
    case('sxy')
      ux_x = -fts*0.5d0*kcs**2
      uy_y =  fts*0.5d0*kcs**2
      if (r>epsilon(r)) then
        ux_y = ftp*kcp**2*x*y*J2p/r**2 + fts*(kcs*(y**2-x**2)*J1s - kcs**2*y**2*r*J0s)/r**3
        uy_x = ftp*kcp**2*x*y*J2p/r**2 - fts*(kcs*(x**2-y**2)*J1s - kcs**2*x**2*r*J0s)/r**3
      end if

      F = G*(ux_y+uy_x)
    case('szz')
      F = 0d0
    end select

  end function inplane_bessel

  function inplane_fault_mms(x,y,t,iblock,field) result(F)

    !use ifport ! for Intel Fortran compiler
    use material, only : block_material
    use io, only : error

    implicit none

    real,intent(in) :: x,y,t
    integer,intent(in) :: iblock
    character(*),intent(in) :: field
    real :: F

    real :: ftp,ftp_t,ftp_tt
    real :: kp,kcp

    real :: vx  ,vy
    real :: vx_t,vy_t
    real :: vx_x,vy_x
    real :: vx_y,vy_y

    real :: v, v_x, v_y, v_t

    real :: nx  ,ny
    real :: nx_x,ny_x
    real :: nx_y,ny_y

    real :: mx  ,my
    real :: mx_x,my_x
    real :: mx_y,my_y

    real :: sxx  ,syy  ,sxy  ,szz
    real :: sxx_t,syy_t,sxy_t,szz_t
    real :: sxx_x,syy_x,sxy_x
    real :: sxx_y,syy_y,sxy_y

    real :: pi

    real :: G = 32d0, cs = 3d0, cp = 5d0, rho, lambda
    rho = G/cs**2
    lambda = rho*cp**2 - 2d0*G

    F = mms_sin(x,y,t,iblock,field)
  end function inplane_fault_mms

  function mms_sin(x,y,t,side,field) result(F)

    !use ifport ! for Intel Fortran compiler
    use io, only : error

    implicit none

    real,intent(in) :: x,y,t
    integer,intent(in) :: side
    character(*),intent(in) :: field
    real :: F, F2

    real :: syy,sxy,sxx,vx,vy,ftp,ftpp,pi,mv
    real :: vx_t,vx_x,vx_y
    real :: vy_t,vy_x,vy_y
    real :: syy_t,syy_x,syy_y
    real :: sxy_t,sxy_x,sxy_y
    real :: sxx_t,sxx_x,sxx_y
    real :: szz_t
    real :: nx,ny,mx,my
    real :: ft,ft_t
    real :: nx_x,ny_x,mx_x,my_x,ftpp_x
    real :: nx_y,ny_y,mx_y,my_y
    real :: w,kx,ky,mf,kf

    real :: G = 32d0, cs = 3d0, cp = 5d0, rho, lambda, Zs
    rho = G/cs**2
    lambda = rho*cp**2 - 2d0*G
    Zs = rho*cs


    pi = 4d0 * datan(1d0)
    kx  = 2d0*pi
    ky  = 2d0*pi
    mv  = 2d0
    w  = (4d0/0.46d0)*pi
    ! w = 10d0*pi
    ft = cos(w*t)
    ft_t = -w*sin(w*t)

    ! WARNING: This must match the coordinate transform!!!!!
    mf   = 0.1d0
    kf   = 1d0*pi
    ftp    = mf*sin(kf*x)
    ftpp   = mf*kf*cos(kf*x)
    ftpp_x =-mf*(kf**2)*sin(kf*x)
    ! ftp = 0d0
    ! ftpp = 0d0
    ! ftpp_x = 0d0

    nx =  ftpp/sqrt(1d0+ftpp**2)
    ny = -1d0 /sqrt(1d0+ftpp**2)
    mx =  ny
    my = -nx

    ! We will also need the derivatives of the normal vectors
    nx_x =     ftpp_x/sqrt((1+ftpp**2)**3)
    ny_x =  ftpp*ftpp_x/sqrt((1+ftpp**2)**3)
    mx_x =  ny_x
    my_x = -nx_x

    nx_y = 0d0
    ny_y = 0d0
    mx_y =  ny_y
    my_y = -nx_y


    ! need to be careful since slip and sxy must have the same sign at the
    ! fault for the friction law
    vx   = mv*(3d0-2d0*side)*   (2d0*12.6/Zs)*cos(kx*x)*cos(ky*y)*ft
    vx_t = mv*(3d0-2d0*side)*   (2d0*12.6/Zs)*cos(kx*x)*cos(ky*y)*ft_t
    vx_x =-mv*(3d0-2d0*side)*kx*(2d0*12.6/Zs)*sin(kx*x)*cos(ky*y)*ft
    vx_y =-mv*(3d0-2d0*side)*ky*(2d0*12.6/Zs)*cos(kx*x)*sin(ky*y)*ft

    vy   =    (2d0*12.6/Zs)*cos(kx*x)*cos(ky*y)*ft
    vy_t =    (2d0*12.6/Zs)*cos(kx*x)*cos(ky*y)*ft_t
    vy_x =-kx*(2d0*12.6/Zs)*sin(kx*x)*cos(ky*y)*ft
    vy_y =-ky*(2d0*12.6/Zs)*cos(kx*x)*sin(ky*y)*ft

    syy   =   -63d0*cos(kx*x)*cos(ky*y)*ft - 126d0
    syy_t =   -63d0*cos(kx*x)*cos(ky*y)*ft_t
    syy_x = kx*63d0*sin(kx*x)*cos(ky*y)*ft
    syy_y = ky*63d0*cos(kx*x)*sin(ky*y)*ft

    sxy   =    0.6d0*126d0*cos(kx*x)*cos(ky*y)*ft
    sxy_t =    0.6d0*126d0*cos(kx*x)*cos(ky*y)*ft_t
    sxy_x =-kx*0.6d0*126d0*sin(kx*x)*cos(ky*y)*ft
    sxy_y =-ky*0.6d0*126d0*cos(kx*x)*sin(ky*y)*ft

    sxx   = (2d0*side-3d0)*   0.6d0*1260*cos(kx*x)*cos(ky*y)*ft
    sxx_t = (2d0*side-3d0)*   0.6d0*1260*cos(kx*x)*cos(ky*y)*ft_t
    sxx_x =-(2d0*side-3d0)*kx*0.6d0*1260*sin(kx*x)*cos(ky*y)*ft
    sxx_y =-(2d0*side-3d0)*ky*0.6d0*1260*cos(kx*x)*sin(ky*y)*ft

    select case(field)
    case('szz')
        F = 0d0
        return
    case('V')
        F = 2d0*vx
        return
    case('vx')
        F = ny*vx-my*vy
        return
    case('vy')
        F = mx*vy-nx*vx
        return
    case('S')
      F = sxy
      return
    case('N')
      F =-syy
      return
    case('sxx')
      F = syy*my**2+sxx*ny**2-2d0*sxy*my*ny
      return
    case('syy')
      F = syy*mx**2+sxx*nx**2-2d0*sxy*mx*nx
      return
    case('sxy')
      F =-syy*my*mx-sxx*ny*nx+sxy*(ny*mx+nx*my)
      return
    case('Vt')
         F = 2d0*vx*ft_t
         return
    case('St')
         F = sxy*ft_t
         return
    case('Nt')
         F =  -syy*ft_t
         return
    case('s_vx')
        ! vx_t - (1/rho)*(sxx_x+sxy_y)
        vx_t  = ny*vx_t-my*vy_t
        sxx_x = syy_x*my**2+sxx_x*ny**2-2*sxy_x*my*ny &
              + 2d0*syy*my*my_x+2d0*sxx*ny*ny_x-2*sxy*(my_x*ny+my*ny_x)
        sxy_y = -syy_y*my*mx-sxx_y*ny*nx+sxy_y*(ny*mx+nx*my) &
              -syy*(my_y*mx+my*mx_y)-sxx*(ny_y*nx+ny*nx_y)+sxy*(ny_y*mx+ny*mx_y+nx_y*my+nx*my_y)
        F = vx_t - (1d0/rho)*(sxx_x+sxy_y)
        return
    case('s_vy')
        ! vy_t - (1/rho)*(sxy_x+syy_y)
        vy_t = mx*vy_t-nx*vx_t
        sxy_x = -syy_x*my*mx-sxx_x*ny*nx+sxy_x*(ny*mx+nx*my) &
              -syy*(my_x*mx+my*mx_x)-sxx*(ny_x*nx+ny*nx_x)+sxy*(ny_x*mx+ny*mx_x+nx_x*my+nx*my_x)
        syy_y = syy_y*mx**2+sxx_y*nx**2-2*sxy_y*mx*nx &
              + 2d0*syy*mx*mx_y+2d0*sxx*nx*nx_y-2*sxy*(mx_y*nx+mx*nx_y)
        F = vy_t - (1d0/rho)*(sxy_x+syy_y)
        return
    case('s_sxx')
        ! sxx_t - (lam+2*G) vx_x - lam * vy_y
        sxx_t = syy_t*my**2+sxx_t*ny**2-2*sxy_t*my*ny
        vx_x = ny_x*vx-my_x*vy + ny*vx_x-my*vy_x
        vy_y = mx_y*vy-nx_y*vx + mx*vy_y-nx*vx_y
        F = sxx_t - (lambda+2d0*G)*vx_x - lambda*vy_y
        return
    case('s_syy')
        ! syy_t - lam * vx_x - (lam+2*G) vy_y
        vx_x = ny_x*vx-my_x*vy + ny*vx_x-my*vy_x
        vy_y = mx_y*vy-nx_y*vx + mx*vy_y-nx*vx_y
        syy_t = syy_t*mx**2+sxx_t*nx**2-2*sxy_t*mx*nx
        F = syy_t - lambda*vx_x - (lambda+2d0*G)*vy_y
        return
    case('s_sxy')
        ! sxy_t - G * vy_x + G * vx_y
        vx_y = ny_y*vx-my_y*vy + ny*vx_y-my*vy_y
        vy_x = mx_x*vy-nx_x*vx + mx*vy_x-nx*vx_x
        sxy_t =-syy_t*my*mx-sxx_t*ny*nx+sxy_t*(ny*mx+nx*my)
        F = sxy_t - G*(vy_x+vx_y)
        return
    case('s_szz')
        ! sxx_t - (lam+2*G) vx_x - lam * vy_y
        szz_t = 0d0
        vx_x = ny_x*vx-my_x*vy + ny*vx_x-my*vy_x
        vy_y = mx_y*vy-nx_y*vx + mx*vy_y-nx*vx_y
        F = szz_t - lambda*vx_x - lambda*vy_y
        return
    end select

    call error('Invalid field (' // trim(field) // ')')

    F = 0d0

  end function mms_sin

  function mms_simple(x,y,t,side,field) result(F)

    !use ifport ! for Intel Fortran compiler

    implicit none

    real,intent(in) :: x,y,t
    integer,intent(in) :: side
    character(*),intent(in) :: field
    real :: F

    real :: r2,syy,sxy,sxx,vx,vy,pi
    real :: vx_t,vx_x,vx_y
    real :: vy_t,vy_x,vy_y
    real :: syy_t,syy_x,syy_y
    real :: sxy_t,sxy_x,sxy_y
    real :: sxx_t,sxx_x,sxx_y
    real :: w,k

    real,parameter :: rho=1d0,lambda=1d0,G=1d0

    ! WARNING: This must match the coordinate transform!!!!!
    pi = 4d0 * datan(1d0)
    w = 2d0*pi
    k = 4d0*pi


    select case(field)
    case('vx','vy','V')
        ! need to be careful since slip and sxy must have the same sign at the
        ! fault for the friction law
        vx = cos(w*t)*cos(k*x)*cos(k*y)
        vy = cos(w*t)*cos(k*x)*cos(k*y)
        select case(field)
        case('V')
            F = 2*vx
            return
        case('vx')
            F = vx
            return
        case('vy')
            F = vy
            return
        end select
    case('sxx','syy','sxy','N','S')
        syy = cos(w*t)*cos(k*x)*cos(k*y)
        sxy = cos(w*t)*cos(k*x)*cos(k*y)
        sxx = cos(w*t)*cos(k*x)*cos(k*y)
        select case(field)
        case('S')
            F = sxy
            return
        case('N')
            F = syy
            return
        case('sxx')
            F = sxx
            return
        case('syy')
            F = syy
            return
        case('sxy')
            F = sxy
            return
        end select
    case('s_vx','s_vy')
        ! need to be careful since slip and sxy must have the same sign at the
        ! fault for the friction law

        ! Need the time derivatives of velocity
        vx_t = -w*sin(w*t)*cos(k*x)*cos(k*y)
        vy_t = -w*sin(w*t)*cos(k*x)*cos(k*y)

        ! For the stresses we need both parts
        syy_x = -k*cos(w*t)*sin(k*x)*cos(k*y)
        sxy_x = -k*cos(w*t)*sin(k*x)*cos(k*y)
        sxx_x = -k*cos(w*t)*sin(k*x)*cos(k*y)

        syy_y = -k*cos(w*t)*cos(k*x)*sin(k*y)
        sxy_y = -k*cos(w*t)*cos(k*x)*sin(k*y)
        sxx_y = -k*cos(w*t)*cos(k*x)*sin(k*y)

        select case(field)
        case('s_vx')
            ! vx_t - (1/rho)*(sxx_x+sxy_y)
            F = vx_t - (1/rho)*(sxx_x+sxy_y)
            return
        case('s_vy')
            ! vy_t - (1/rho)*(sxy_x+syy_y)
            F = vy_t - (1/rho)*(sxy_x+syy_y)
            return
        end select

    case('s_sxx','s_syy','s_sxy')
        ! Need the time derivatives of stresses
        syy_t = -w*sin(w*t)*cos(k*x)*cos(k*y)
        sxy_t = -w*sin(w*t)*cos(k*x)*cos(k*y)
        sxx_t = -w*sin(w*t)*cos(k*x)*cos(k*y)

        ! For the velocities we need both parts
        vx_x = -k*cos(w*t)*sin(k*x)*cos(k*y)
        vy_x = -k*cos(w*t)*sin(k*x)*cos(k*y)

        vx_y = -k*cos(w*t)*cos(k*x)*sin(k*y)
        vy_y = -k*cos(w*t)*cos(k*x)*sin(k*y)

        select case(field)
        ! For the velocities we will use the chain rule and need all of these
        case('s_sxx')
            ! sxx_t - (lam+2*G) vx_x - lam * vy_y
            F = sxx_t - (lambda+2*G)*vx_x - lambda*vy_y
            return
        case('s_syy')
            ! syy_t - lam * vx_x - (lam+2*G) vy_y
            F = syy_t - (lambda+2*G)*vy_y - lambda*vx_x
            return
        case('s_sxy')
            ! sxy_t - G * vy_x + G * vx_y
            F = sxy_t - G*(vy_x+vx_y)
            return
        end select
    end select

    F = 0

  end function mms_simple

  function mms_hydrofrac(x,y,t,side,field) result(F)

    !use ifport ! for Intel Fortran compiler

    implicit none

    real,intent(in) :: x,y,t
    integer,intent(in) :: side
    character(*),intent(in) :: field
    real :: F

    real :: r2,syy,sxy,sxx,vx,vy,pi
    real :: vx_t,vx_x,vx_y
    real :: vy_t,vy_x,vy_y
    real :: M_x, M_y
    real :: syy_t,syy_x,syy_y
    real :: sxy_t,sxy_x,sxy_y
    real :: sxx_t,sxx_x,sxx_y
    real :: w,k,A,B
    real :: p,v,u
    real :: I_w,I_vtp,I_vtm,I_taum,I_taup,s_p,s_v, &
            vtm,vtp,&
            bcL_uhat,bcL_phat,bcR_uhat,bcR_phat, &
            dpdx


    ! MMS using displacement field:
    ! ux =   Ax*cos(k*x + k*y + w*t)
    ! uy = - Ax*cos(k*x + k*y + w*t)
    ! uz = 0

    !WARNING: These parameters must match the parameters in the input file
    real :: G,cs = 3d3, cp = 5d3, rho = 2.6d0, lambda, Zs
    real :: K0 = 1d3, rho0 = 1d0,mu = 1d-6
    
    ! Fracture width and crack tips
    real :: wm0 = -1.0d0, wp0 = 1.0d0, xL = -5d0, xR = 5d0 
    
    cs = 3d0
    cp = 5d0
    rho = 2.6d0
    K0 = 1d0
    rho0 = 1d0
    mu = 1d0

    G = rho*cs**2
    lambda = rho*cp**2 - 2d0*G
    Zs = rho*cs

    pi = 4d0 * datan(1d0)
    w = 2d0*pi ! omega
    k = 2d0*pi/5d0 ! wave number
    A = 2d0*pi/(lambda)   ! amplitude
    B = 1d0 ! Fluid amplitude

! Solid

! Velocities:
vx =  A*w*cos(k*x)*cos(k*y)*cos(t*w)
vy =  -A*w*cos(k*x)*cos(k*y)*cos(t*w)
! Stresses:
sxx =  -2.0*A*G*k*sin(k*x)*sin(t*w)*cos(k*y) + lambda*(-1.0*A*k*sin(k*x)*sin(t*w)*cos(k*y)&
+ 1.0*A*k*sin(k*y)*sin(t*w)*cos(k*x))
syy =  2.0*A*G*k*sin(k*y)*sin(t*w)*cos(k*x) + lambda*(-1.0*A*k*sin(k*x)*sin(t*w)*cos(k*y)&
+ 1.0*A*k*sin(k*y)*sin(t*w)*cos(k*x))
sxy =  2*G*(0.5*A*k*sin(k*x)*sin(t*w)*cos(k*y) - 0.5*A*k*sin(k*y)*sin(t*w)*cos(k*x))
! Source function for momentum balance:
M_x =  -A*w**2*sin(t*w)*cos(k*x)*cos(k*y) - (-2.0*A*G*k**2*sin(t*w)*cos(k*x)*cos(k*y)&
+ 2*G*(-0.5*A*k**2*sin(k*x)*sin(k*y)*sin(t*w) - 0.5*A*k**2*sin(t*w)*cos(k*x)*cos(k*y))&
+ lambda*(-1.0*A*k**2*sin(k*x)*sin(k*y)*sin(t*w) - 1.0*A*k**2*sin(t*w)*cos(k*x)*cos(k*y)))/rho
M_y =  A*w**2*sin(t*w)*cos(k*x)*cos(k*y) - (2.0*A*G*k**2*sin(t*w)*cos(k*x)*cos(k*y)&
+ 2*G*(0.5*A*k**2*sin(k*x)*sin(k*y)*sin(t*w) + 0.5*A*k**2*sin(t*w)*cos(k*x)*cos(k*y))&
+ lambda*(1.0*A*k**2*sin(k*x)*sin(k*y)*sin(t*w) + 1.0*A*k**2*sin(t*w)*cos(k*x)*cos(k*y)))/rho

! Fluid

! Solutions
p =  1.0*A*k*lambda*sin(k*x)*sin(t*w)
v =  w*cos(k*x)*cos(t*w)
u =  w*cos(k*x)*cos(t*w)

! Forcing functions

! Boundary conditions

 bcL_uhat =  w*cos(k*xL)*cos(t*w)
 bcL_phat =  1.0*A*k*lambda*sin(k*xL)*sin(t*w)
 bcR_uhat =  w*cos(k*xR)*cos(t*w)
 bcR_phat =  1.0*A*k*lambda*sin(k*xR)*sin(t*w)

! Interface conditions

! Normal velocity
I_w =  0
! Tangential velocity
! v - vtm = I_vtm
I_vtm =  -A*w*cos(k*x)*cos(t*w) + w*cos(k*x)*cos(t*w)
I_vtp =  -A*w*cos(k*x)*cos(t*w) + w*cos(k*x)*cos(t*w)
! Tractions
I_taum =  1.0*A*G*k*sin(k*x)*sin(t*w)
I_taup =  1.0*A*G*k*sin(k*x)*sin(t*w)

! Governing equations

! Mass balance
s_p =  k*w*(1.0*A*lambda - K0)*sin(k*x)*cos(t*w)
! Momentum balance
s_v =  1.0*A*k**2*lambda*sin(t*w)*cos(k*x)/rho0 - w**2*sin(t*w)*cos(k*x)

! Exact interface conditions
vtm =  w*cos(k*x)*cos(t*w)
vtp =  w*cos(k*x)*cos(t*w)

! Other
dpdx =  1.0*A*k**2*lambda*sin(t*w)*cos(k*x)



    ! Fields
    select case(field)
    case('vx','vy')
        select case(field)
        case('vx')
            F = vx
            return
        case('vy')
            F = vy
            return
        end select
    case('sxx','syy','sxy')
        select case(field)
        case('sxx')
            F = sxx
            return
        case('syy')
            F = syy
            return
        case('sxy')
            F = sxy
            return
        end select
    ! Source (force function)
    case('s_vx','s_vy')
        select case(field)
        case('s_vx')
            F = M_x
            return
        case('s_vy')
            F = M_y
            return
        end select
    case('s_sxx','s_syy','s_sxy')
        select case(field)
        case('s_sxx')
            F = 0d0 
            return
        case('s_syy')
            F = 0d0
            return
        case('s_xy')
            F = 0d0
            return
        end select
    case('p','v')
        select case(field)
        case('p')
            F = p
            return
        case('v')
            F = v
            return
        end select
    ! Forcing function for fluid mass and momentum balance
    case('s_p','s_v')
        select case(field)
        case('s_v')
            F = s_v 
            return
        case('s_p')
            F = s_p
            return
        end select
    ! Forcing function for boundary conditions
    case('bcL_uhat','bcL_phat','bcR_uhat','bcR_phat')
        select case(field)
        case('bcL_uhat')
            F = bcL_uhat
            return
        case('bcL_phat')
            F = bcL_phat
            return
        case('bcR_uhat')
            F = bcR_uhat
            return
        case('bcR_phat')
            F = bcR_phat
            return
        end select
    ! Forcing function for interface conditions
    case('I_w','I_vm','I_vp','I_taum','I_taup')
        select case(field)
        case('I_w')
            F = I_w
            return
        case('I_vtm')
            F = I_vtm
            return
        case('I_vtp')
            F = I_vtp
            return
        case('I_taum')
            F = I_taum
            return
        case('I_taup')
            F = I_taup
            return
        end select
    case('vtm','vtp')
        select case(field)
        case('vtp')
            F = vtp
            return
        case('vtm')
            F = vtm
            return
        end select
    case('dpdx')
        F = dpdx
        return
    end select

    F = 0d0

  end function mms_hydrofrac

end module mms
