module friction

  implicit none

  type :: ratestate_constant
     real :: a,b,V0,f0,L,fw,Vw
  end type ratestate_constant

  type :: slipweak_constant
     real :: c0,fs,fd,Dc
  end type slipweak_constant

  type :: ratestate
     real,dimension(:),allocatable :: a,b,V0,f0,L,fw,Vw
  end type ratestate

  type :: slipweak
     real,dimension(:),allocatable :: c0,fs,fd,Dc
  end type slipweak

  type :: kinematic
     real :: c0,fd,dfdx,vL,vR,xL,xR,tmax
  end type kinematic

  type :: pseudodynamic
     real,dimension(:),allocatable :: Dmax,trup,Vpeak
  end type pseudodynamic

  type :: load
     character(256) :: shape
     real :: S,N,T,R,x0,y0
  end type load

  type :: fr_type
     logical :: opening,force
     character(256) :: friction_law,problem,rup_field
     real :: Psi0,f_S0,angle,rup_threshold,uni_x0,uni_dSdx,xlockm,xlockp,flock, &
     subduction_f0,subduction_N0,subduction_B
     type(ratestate) :: rs
     type(slipweak) :: sw
     type(kinematic) :: kn
     type(pseudodynamic) :: pd
     type(load) :: ld
     real,dimension(:),allocatable :: D,Psi,DPsi,trup,Ds,Dn,F,V,O,S,N,S0,N0,DDs,DDn
  end type fr_type


contains

  
  subroutine init_friction(iface,FR,m,p,input,echo,skip, &
       process_m,process_p,comm_m,comm_p,array)

    use mpi_routines, only : is_master
    use io, only : error,write_matlab,seek_to_string  

    implicit none

    integer,intent(in) :: iface,m,p,input,echo,comm_m,comm_p,array
    type(fr_type),intent(out) :: FR
    logical,intent(in) :: skip,process_m,process_p

    integer :: stat
    character(256) :: FRstr
    character(256) :: str

    logical :: opening,force,friction_file
    character(256) :: friction_law,problem,rup_field,filename
    real :: Psi0,S0,angle,rup_threshold,uni_x0,uni_dSdx,xlockm,xlockp,flock, &
    subduction_f0,subduction_N0,subduction_B
    type(ratestate_constant) :: rs
    type(slipweak_constant) :: sw
    type(kinematic) :: kn
    type(load) :: ld

    namelist /friction_list/ opening,force,friction_law,problem,friction_file,filename, &
         S0,Psi0,angle,rs,sw,kn,ld,rup_field,rup_threshold,uni_x0,uni_dSdx,xlockm,xlockp,flock, &
         subduction_f0,subduction_N0,subduction_B

    ! defaults

    opening = .true.
    force = .false.
    friction_law = 'frictionless'
    problem = ''
    friction_file = .false.
    filename = ''
    Psi0 = 0d0
    S0 = 0d0
    angle = 0d0
    rs = ratestate_constant(0.01d0,0.014d0,1d-6,0.6d0,0.4d0,0.2d0,1d20)
    sw = slipweak_constant(0d0,1d0,0d0,1d0)
    kn = kinematic(0d0,0d0,1d0,0d0,0d0,0d0,0d0,huge(1d0))
    ld = load('smooth',0d0,0d0,0d0,1d0,0d0,0d0)
    rup_field = 'V'
    rup_threshold = 0d0
    uni_x0 = 0d0
    uni_dSdx = 0d0
    xlockm = -1d10
    xlockp = 1d10
    flock = 1d10
    subduction_f0=0d0
    subduction_N0=0d0
    subduction_B=0d0

    ! read in friction parameters

    write(str,'(a,i0,a)') '!---IFACE',iface,'---'
    call seek_to_string(input,str)
    read(input,nml=friction_list,iostat=stat)
    if (stat>0) call error('Error in friction_list','init_friction')

    FR%opening = opening
    FR%force = force
    FR%friction_law = friction_law
    FR%problem = problem
    FR%Psi0 = Psi0
    FR%f_S0 = S0
    FR%angle = 0.017453292519943d0*angle ! convert deg to rad
    FR%kn = kn
    FR%ld = ld
    FR%rup_field = rup_field
    FR%rup_threshold = rup_threshold
    FR%uni_x0 = uni_x0
    FR%uni_dSdx = uni_dSdx
    FR%xlockm = xlockm
    FR%xlockp = xlockp
    FR%flock = flock
    FR%subduction_f0 = subduction_f0
    FR%subduction_N0 = subduction_N0
    FR%subduction_B = subduction_B

    ! output friction parameters
    
    if (is_master) then

       write(FRstr,'(a,i0,a)') 'FR{',iface,'}'

       call write_matlab(echo,'friction_law',FR%friction_law,FRstr)
       call write_matlab(echo,'problem',FR%problem,FRstr)
       call write_matlab(echo,'opening',FR%opening,FRstr)
       call write_matlab(echo,'force',FR%force,FRstr)
       call write_matlab(echo,'angle',FR%angle,FRstr)
       call write_matlab(echo,'S0',FR%f_S0,FRstr)
       call write_matlab(echo,'subduction_f0',FR%subduction_f0,FRstr)
       call write_matlab(echo,'subduction_N0',FR%subduction_N0,FRstr)
       call write_matlab(echo,'subduction_B',FR%subduction_B,FRstr)

       call write_matlab(echo,'xlockm',FR%xlockm,FRstr)
       call write_matlab(echo,'xlockp',FR%xlockp,FRstr)
       call write_matlab(echo,'flock',FR%flock,FRstr)

       select case(FR%friction_law)

       case('SL','FL','RSL','RFL','RSF','RSL-mms','RSL-mms-nostate')

          call write_matlab(echo,'Psi0',FR%Psi0,FRstr)
          call write_matlab(echo,'rs.a' ,rs%a ,FRstr)
          call write_matlab(echo,'rs.b' ,rs%b ,FRstr)
          call write_matlab(echo,'rs.V0',rs%V0,FRstr)
          call write_matlab(echo,'rs.f0',rs%f0,FRstr)
          call write_matlab(echo,'rs.L' ,rs%L ,FRstr)
          call write_matlab(echo,'rs.fw',rs%fw,FRstr)
          call write_matlab(echo,'rs.Vw',rs%Vw,FRstr)

       case('SW')

          call write_matlab(echo,'sw.c0',sw%c0,FRstr)
          call write_matlab(echo,'sw.fs',sw%fs,FRstr)
          call write_matlab(echo,'sw.fd',sw%fd,FRstr)
          call write_matlab(echo,'sw.Dc',sw%Dc,FRstr)

       end select

       call write_matlab(echo,'kn.c0'  ,FR%kn%c0  ,FRstr)
       call write_matlab(echo,'kn.fd'  ,FR%kn%fd  ,FRstr)
       call write_matlab(echo,'kn.dfdx',FR%kn%dfdx,FRstr)
       call write_matlab(echo,'kn.vL'  ,FR%kn%vL  ,FRstr)
       call write_matlab(echo,'kn.vR'  ,FR%kn%vR  ,FRstr)
       call write_matlab(echo,'kn.xL'  ,FR%kn%xL  ,FRstr)
       call write_matlab(echo,'kn.xR'  ,FR%kn%xR  ,FRstr)
       call write_matlab(echo,'kn.tmax',FR%kn%tmax,FRstr)

       call write_matlab(echo,'ld.shape',FR%ld%shape,FRstr)
       call write_matlab(echo,'ld.S'    ,FR%ld%S    ,FRstr)
       call write_matlab(echo,'ld.N'    ,FR%ld%N    ,FRstr)
       call write_matlab(echo,'ld.T'    ,FR%ld%T    ,FRstr)
       call write_matlab(echo,'ld.R'    ,FR%ld%R    ,FRstr)
       call write_matlab(echo,'ld.x0'   ,FR%ld%x0   ,FRstr)
       call write_matlab(echo,'ld.y0'   ,FR%ld%y0   ,FRstr)

       call write_matlab(echo,'rup.field'    ,FR%rup_field    ,FRstr)
       call write_matlab(echo,'rup.threshold',FR%rup_threshold,FRstr)

       call write_matlab(echo,'un.x0'  ,FR%uni_x0  ,FRstr)
       call write_matlab(echo,'un.dSdx',FR%uni_dSdx,FRstr)

    end if

    ! return if not needed

    if (skip) return

    ! allocate and initialize arrays

    allocate(FR%V   (m:p),FR%O   (m:p))
    allocate(FR%S   (m:p),FR%N   (m:p))
    allocate(FR%S0  (m:p),FR%N0  (m:p))
    allocate(FR%Ds  (m:p),FR%DDs (m:p))
    allocate(FR%Dn  (m:p),FR%DDn (m:p))

    FR%Ds  = 0d0
    FR%Dn  = 0d0
    FR%DDs = 0d0
    FR%DDn = 0d0
    FR%V    = 1d40
    FR%O    = 1d40
    FR%S    = 1d40
    FR%N    = 1d40
    FR%S0   = 1d40
    FR%N0   = 1d40
    
    allocate(FR%Psi(m:p),FR%DPsi(m:p),FR%D(m:p),FR%trup(m:p))

    FR%Psi  = 1d40
    FR%DPsi = 1d40
    FR%D = 0d0
    FR%trup = 1d10

    select case(FR%friction_law)
    case('SL','FL','RSL','RFL','RSF','RSL-mms','RSL-mms-nostate')
       allocate(FR%rs%a(m:p),FR%rs%b(m:p),FR%rs%V0(m:p),FR%rs%f0(m:p),FR%rs%L(m:p),FR%rs%fw(m:p),FR%rs%Vw(m:p))
       FR%rs%a  = rs%a
       FR%rs%b  = rs%b
       FR%rs%V0 = rs%V0
       FR%rs%f0 = rs%f0
       FR%rs%L  = rs%L
       FR%rs%fw = rs%fw
       FR%rs%Vw = rs%Vw
    case('SW')
       allocate(FR%sw%c0(m:p),FR%sw%fs(m:p),FR%sw%fd(m:p),FR%sw%Dc(m:p))
       FR%sw%c0 = sw%c0
       FR%sw%fs = sw%fs
       FR%sw%fd = sw%fd
       FR%sw%Dc = sw%Dc
    case('pseudodynamic')
       allocate(FR%pd%Dmax(m:p),FR%pd%trup(m:p),FR%pd%Vpeak(m:p))
       FR%pd%Dmax = 0d0
       FR%pd%trup = 1d40
       FR%pd%Vpeak = 0d0
    end select

    if (friction_file) then
       ! both sides read file (so process may read file twice)
       if (process_m) call read_friction(FR,filename,comm_m,array)
       if (process_p) call read_friction(FR,filename,comm_p,array)
    end if

  end subroutine init_friction


  subroutine destroy_friction(FR)

    implicit none

    type(fr_type),intent(inout) :: FR

    if (allocated(FR%Psi )) deallocate(FR%Psi )
    if (allocated(FR%DPsi)) deallocate(FR%DPsi)
    if (allocated(FR%D   )) deallocate(FR%D   )
    if (allocated(FR%trup)) deallocate(FR%trup)
    if (allocated(FR%V   )) deallocate(FR%V   )
    if (allocated(FR%O   )) deallocate(FR%O   )
    if (allocated(FR%S   )) deallocate(FR%S   )
    if (allocated(FR%N   )) deallocate(FR%N   )
    if (allocated(FR%S0  )) deallocate(FR%S0  )
    if (allocated(FR%N0  )) deallocate(FR%N0  )
    if (allocated(FR%Ds  )) deallocate(FR%Ds  )
    if (allocated(FR%Dn  )) deallocate(FR%Dn  )
    if (allocated(FR%DDs )) deallocate(FR%DDs )
    if (allocated(FR%DDn )) deallocate(FR%DDn )

    if (allocated(FR%rs%a )) deallocate(FR%rs%a )
    if (allocated(FR%rs%b )) deallocate(FR%rs%b )
    if (allocated(FR%rs%V0)) deallocate(FR%rs%V0)
    if (allocated(FR%rs%f0)) deallocate(FR%rs%f0)
    if (allocated(FR%rs%L )) deallocate(FR%rs%L )
    if (allocated(FR%rs%fw)) deallocate(FR%rs%fw)
    if (allocated(FR%rs%Vw)) deallocate(FR%rs%Vw)

    if (allocated(FR%sw%c0)) deallocate(FR%sw%c0)
    if (allocated(FR%sw%fs)) deallocate(FR%sw%fs)
    if (allocated(FR%sw%fd)) deallocate(FR%sw%fd)
    if (allocated(FR%sw%Dc)) deallocate(FR%sw%Dc)

    if (allocated(FR%pd%Dmax)) deallocate(FR%pd%Dmax)
    if (allocated(FR%pd%trup)) deallocate(FR%pd%trup)
    if (allocated(FR%pd%Vpeak)) deallocate(FR%pd%Vpeak)

  end subroutine destroy_friction


  subroutine read_friction(FR,filename,comm,array)

    use io, only : file_distributed,open_file_distributed, &
         read_file_distributed,close_file_distributed
    use mpi_routines, only : pw

    implicit none

    type(fr_type),intent(inout) :: FR
    character(*),intent(in) :: filename
    integer,intent(in) :: comm,array

    type(file_distributed) :: fh

    call open_file_distributed(fh,filename,'read',comm,array,pw)

    select case(FR%friction_law)

    case('SL','FL','RSL','RFL','RSF','RSL-mms','RSL-mms-nostate')
       
       call read_file_distributed(fh,FR%rs%a )
       call read_file_distributed(fh,FR%rs%b )
       call read_file_distributed(fh,FR%rs%V0)
       call read_file_distributed(fh,FR%rs%f0)
       call read_file_distributed(fh,FR%rs%L )
       call read_file_distributed(fh,FR%rs%fw)
       call read_file_distributed(fh,FR%rs%Vw)

    case('SW')

       call read_file_distributed(fh,FR%sw%c0)
       call read_file_distributed(fh,FR%sw%fs)
       call read_file_distributed(fh,FR%sw%fd)
       call read_file_distributed(fh,FR%sw%Dc)

    case('pseudodynamic')

       call read_file_distributed(fh,FR%pd%Dmax)
       call read_file_distributed(fh,FR%pd%trup)
       call read_file_distributed(fh,FR%pd%Vpeak)

    end select

    call close_file_distributed(fh)

  end subroutine read_friction


  subroutine checkpoint_friction(fh,operation,FR)

    use io, only : file_distributed, &
         read_file_distributed,write_file_distributed

    implicit none

    type(file_distributed),intent(in) :: fh
    character(*),intent(in) :: operation
    type(fr_type),intent(inout) :: FR

    ! fields read/written to same file as iface fields,
    ! routine is only called by io_process

    select case(operation)
    case('read')
       call read_file_distributed(fh,FR%D)
       call read_file_distributed(fh,FR%Psi)
       call read_file_distributed(fh,FR%DPsi)
       call read_file_distributed(fh,FR%trup)
       call read_file_distributed(fh,FR%V)
       call read_file_distributed(fh,FR%O)
       call read_file_distributed(fh,FR%S)
       call read_file_distributed(fh,FR%N)
       call read_file_distributed(fh,FR%S0)
       call read_file_distributed(fh,FR%N0)
       call read_file_distributed(fh,FR%Ds)
       call read_file_distributed(fh,FR%Dn)
       call read_file_distributed(fh,FR%DDs)
       call read_file_distributed(fh,FR%DDn)
    case('write')
       call write_file_distributed(fh,FR%D)
       call write_file_distributed(fh,FR%Psi)
       call write_file_distributed(fh,FR%DPsi)
       call write_file_distributed(fh,FR%trup)
       call write_file_distributed(fh,FR%V)
       call write_file_distributed(fh,FR%O)
       call write_file_distributed(fh,FR%S)
       call write_file_distributed(fh,FR%N)
       call write_file_distributed(fh,FR%S0)
       call write_file_distributed(fh,FR%N0)
       call write_file_distributed(fh,FR%Ds)
       call write_file_distributed(fh,FR%Dn)
       call write_file_distributed(fh,FR%DDs)
       call write_file_distributed(fh,FR%DDn)
    end select

  end subroutine checkpoint_friction


  subroutine scale_rates_friction(FR,A)

    implicit none

    type(fr_type),intent(inout) :: FR
    real,intent(in) :: A

    FR%DDs  = A*FR%DDs
    FR%DDn  = A*FR%DDn
    FR%DPsi = A*FR%DPsi

  end subroutine scale_rates_friction
  

  subroutine update_fields_friction(FR,dt)

    implicit none

    type(fr_type),intent(inout) :: FR
    real,intent(in) :: dt

    FR%Ds  = FR%Ds +dt*FR%DDs
    FR%Dn  = FR%Dn +dt*FR%DDn
    FR%D   = FR%D  +dt*abs(FR%DDs)
    FR%Psi = FR%Psi+dt*FR%DPsi

  end subroutine update_fields_friction


  subroutine initial_state(FR,Psi,V,S,N,i,x,y,dt)

    use io, only : error
    use mms, only: inplane_fault_mms

    implicit none

    type(fr_type),intent(in) :: FR
    real,intent(inout) :: Psi
    real,intent(in) :: V,S,N,x,y,dt
    integer,intent(in) :: i

    real :: Vex,Nex,Sex
    character(256) :: str
    type(ratestate_constant) :: rs
    real :: theta0

    select case(FR%friction_law)
    case('frictionless','pseudodynamic','SW')
       Psi = 0d0
       return
    end select

    call ratestate_param(i,x,y,FR,rs)

    if (dt>0d0) then
       theta0 = (rs%L/rs%V0)*exp((Psi-rs%f0)/rs%b)
       Psi = Psi+rs%b*log(1d0+dt/theta0)
       return
    end if

    ! set Psi to constant initial value if V = 0
    ! (which probably does not satisfy S = f(V,Psi)*N)

    select case(FR%problem)
    case default
       if (abs(V)==0d0) then
          Psi = FR%Psi0
          return
       end if
    case('tpv105')
       Psi = rs%a*log((2d0*rs%V0/1d-16)*sinh(FR%f_S0/(rs%a*N)))
       return
    end select

    ! otherwise set Psi such that S = f(V,Psi)*N

    select case(FR%friction_law)
    case('SL','FL')
       Psi = abs(S/N)-rs%a*log(abs(V)/rs%V0)
    case('RSL','RFL','RSF')
       if (V*S<0d0) then
          write(str,*) 'Incompatible signs: V = ',V,' S = ',S,' at x = ',x,', y = ',y
          call error(str,'initial_state')
       else
          Psi = rs%a*log(2d0*rs%V0/V*sinh(S/(rs%a*N)))
       end if
    case('RSL-mms','RSL-mms-nostate')
       if (V*S<0d0) then
          write(str,*) 'Incompatible signs: V = ',V,' S = ',S,' at x = ',x,', y = ',y
          call error(str,'initial_state')
       else
          Vex = inplane_fault_mms(x,y,0d0,1,'V')
          Sex = inplane_fault_mms(x,y,0d0,1,'S')
          Nex = inplane_fault_mms(x,y,0d0,1,'N')
          Psi = rs%a*log(2d0*rs%V0*sinh(Sex/(rs%a*Nex))/Vex)
          Psi = 0d0
       end if
    case default
       Psi = 0d0
    end select

  end subroutine initial_state


  subroutine solve_friction(FR,V,S,N,Slock,eta,D,Psi,i,x,y,t,info)

    use mms, only : bessel,inplane_fault_mms
    use io, only : error

    implicit none

    type(fr_type),intent(in) :: FR
    real,intent(inout) :: V,S
    real,intent(in) :: N,Slock,eta,D,Psi,x,y,t
    integer,intent(in) :: i
    logical,intent(in) :: info

    logical :: fail
    real :: Sk,xm,xp,fkm,fkp,xx,Vex,Sex,Nex,Psiex
    type(slipweak_constant) :: sw
    type(ratestate_constant) :: rs

    ! lock fault outside specified region, specified in rotated coordinates
    
    xx =  x*cos(FR%angle)+y*sin(FR%angle) ! distance along fault

    if (xx<FR%xlockm.or.xx>FR%xlockp) then
       S = FR%flock*N
       if (abs(Slock)>=S) then ! slipping
          S = sign(S,Slock)
          V = (Slock-S)/eta
       else ! locked
          V = 0d0
          S = Slock
       end if
       return
    end if

    ! special treatment of simple fault behaviors

    select case(FR%friction_law)
    case('frictionless')
       S = 0d0
       V = Slock/eta
       return
    case('locked')
       V = 0d0
       S = Slock
       return
    case('pseudodynamic')
       V = pseudodynamicV(i,t,FR%pd)
       S = Slock-eta*V
       return
    end select

    ! frictional strength from kinematic forcing

    if (FR%force) then

       ! work in rotated coordinates

       xx =  x*cos(FR%angle)+y*sin(FR%angle) ! distance along fault
       !yy = -x*sin(FR%angle)+y*cos(FR%angle)

       ! left tip
       xm = FR%kn%xL+FR%kn%vL*min(t,FR%kn%tmax)
       if (xx<xm) then
          fkm = FR%kn%fd+FR%kn%dfdx*(xm-xx)
       else
          fkm = FR%kn%fd
       end if

       ! right tip
       xp = FR%kn%xR+FR%kn%vR*min(t,FR%kn%tmax)
       if (xx>xp) then
          fkp = FR%kn%fd+FR%kn%dfdx*(xx-xp)
       else
          fkp = FR%kn%fd
       end if

       ! combine
       Sk = FR%kn%c0+N*max(fkm,fkp)

    end if

    ! frictional strength set by friction law

    select case(FR%friction_law)

    case default

       call error('Invalid friction law','solve_friction')

    case('SW') ! slip-weakening friction

       call slipweak_param(i,x,y,FR,sw)

       ! friction strength from slip-weakening friction (S non-negative)
       
       S = sw%c0+max(N,0d0)* &
            ((sw%fs-(sw%fs-sw%fd)*min(D,sw%Dc)/sw%Dc))

       ! if forcing, strength is minimum of frictional strength and forced strength

       if (FR%force) S = min(S,Sk)

       ! set stress and slip velocity (with appropriate signs)

       if (abs(Slock)>=S) then ! slipping
          S = sign(S,Slock)
          V = (Slock-S)/eta
       else ! locked
          V = 0d0
          S = Slock
       end if

    case('SL','FL','RSL','RFL','RSF','RSL-mms','RSL-mms-nostate') ! rate-and-state friction

       call ratestate_param(i,x,y,FR,rs)

       ! solve nonlinear equation with Newton's method

       select case(FR%problem)
       case default
          call newton_solver(FR%friction_law,rs,V,S,Slock,N,eta,Psi,x,y,t,info,fail)
          if (fail) call newton_solver(FR%friction_law,rs,V,S,Slock,N,eta,Psi,x,y,t,.true.,fail)
       case('bessel')
          Vex = bessel(x,y,t,1,'V')
          Sex = bessel(x,y,t,1,'S')
          Psiex = rs%a*log(2d0*rs%V0*sinh(Sex/(rs%a*N))/Vex)
          !if (x==0d0) print *, 'ex',x,t,Vex,Sex,Psiex
          call newton_solver(FR%friction_law,rs,V,S,Slock,N,eta,Psiex,x,y,t,info,fail)
          if (fail) call newton_solver(FR%friction_law,rs,V,S,Slock,N,eta,Psiex,x,y,t,.true.,fail)
       case('inplane-fault-mms-nostate')
          Vex = inplane_fault_mms(x,y,t,1,'V')
          Sex = inplane_fault_mms(x,y,t,1,'S')
          Nex = inplane_fault_mms(x,y,t,1,'N')
          Psiex = rs%a*log(2d0*rs%V0*sinh(Sex/(rs%a*Nex))/Vex)
          !if (x==0d0) print *, 'ex',x,t,Vex,Sex,Psiex
          call newton_solver(FR%friction_law,rs,V,S,Slock,N,eta,Psiex,x,y,t,info,fail)
          if (fail) call newton_solver(FR%friction_law,rs,V,S,Slock,N,eta,Psiex,x,y,t,.true.,fail)
       case('inplane-fault-mms')
          Vex = inplane_fault_mms(x,y,t,1,'V')
          Sex = inplane_fault_mms(x,y,t,1,'S')
          Nex = inplane_fault_mms(x,y,t,1,'N')
          Psiex = rs%a*log(2d0*rs%V0*sinh(Sex/(rs%a*Nex))/Vex)
          call newton_solver(FR%friction_law,rs,V,S,Slock,N,eta,Psiex+Psi,x,y,t,info,fail)
          if (fail) call newton_solver(FR%friction_law,rs,V,S,Slock,N,eta,Psiex+Psi,x,y,t,.true.,fail)
       end select

       ! if forcing, and forced strength is less than frictional strength,
       ! adjust stress and slip velocity (otherwise rate-and-state prevails)

       if (FR%force) then
          if (Sk<abs(S)) then
             S = sign(Sk,Slock)
             V = (Slock-S)/eta
          end if
       end if

    end select

  end subroutine solve_friction


  subroutine slipweak_param(i,x,y,FR,sw)

    use utilities, only : boxcar

    implicit none

    integer,intent(in) :: i
    real,intent(in) :: x,y
    type(fr_type),intent(in) :: FR
    type(slipweak_constant),intent(out) :: sw

    real :: dip
    real,parameter :: tol = 100d0*epsilon(x)

    ! unmodified parameters
    
    sw = slipweak_constant(FR%sw%c0(i),FR%sw%fs(i),FR%sw%fd(i),FR%sw%Dc(i))

    select case(FR%problem)
    case('TPV5','TPV205')
       ! barrier at ends
       if (abs(x)>15d0+tol) sw%c0 = 1d10
    case('TPV10','TPV11','TPV210')
       dip = 2d0*x
       ! barrier at depth
       if (dip>15d0+tol) sw%c0 = 1d10
    case('TPV12','TPV13','TPV12alt','TPV13alt')
       dip = 2d0*x
       ! reduced friction coefficient inside nucleation region
       sw%fs = boxcar(dip-12d0,1.5d0,0.54d0,0.7d0)
       ! barrier at depth
       if (dip>15d0+tol) sw%c0 = 1d10
    case('TPV14BRANCH','TPV15BRANCH','TPV14bBRANCH','TPV15bBRANCH')
       ! first node on the branch is locked
       if (x**2+y**2<tol) sw%c0 = 1d10
    end select

  end subroutine slipweak_param


  subroutine ratestate_param(i,x,y,FR,rs)

    use utilities, only : smooth_boxcar

    implicit none

    integer,intent(in) :: i
    real,intent(in) :: x,y
    type(fr_type),intent(in) :: FR
    type(ratestate_constant),intent(out) :: rs

    ! unmodified parameters

    rs = ratestate_constant(FR%rs%a(i),FR%rs%b(i),FR%rs%V0(i), &
         FR%rs%f0(i),FR%rs%L(i),FR%rs%fw(i),FR%rs%Vw(i))
    
    ! modification for specific problems

    select case(FR%problem)
    case('rough-old') ! increase in a and Vw with decreasing x
       rs%a = rs%a+0.008d0*max(0d0,min(1d0,-(x+6d0)/6d0))
       rs%Vw = rs%Vw+5d0*max(0d0,-(x+6d0)/6d0)
    case('tpv105')
       rs%a  = rs%a +smooth_boxcar(abs(x),15d0,3d0,0d0,0.01d0)
       rs%Vw = rs%Vw+smooth_boxcar(abs(x),15d0,3d0,0d0, 0.9d0)
    end select

  end subroutine ratestate_param


  subroutine newton_solver(friction_law,rs,V,S,Slock,N,eta,Psi,x,y,t,info,fail)

    use io, only : warning

    implicit none

    character(*),intent(in) :: friction_law
    type(ratestate_constant),intent(in) :: rs
    real,intent(inout) :: V,S
    real,intent(in) :: Slock,N,eta,Psi,x,y,t
    logical,intent(in) :: info
    logical,intent(out) :: fail
    
    character(256) :: str
    integer :: i
    integer,parameter :: imax=99
    real :: inV,inS,Vmin,Vmax,dSdV,Sf,dSfdV,dV,R,dRdV
    !real,parameter :: rtolV=0d-12, atolV=0d0, atolR=1d-14
    !real,parameter :: rtolV=1d-12, atolV=0d-20, atolR=1d-12
    !real,parameter :: rtolV=1d-12, atolV=0d0, atolR=epsilon(1d0)*1d4
    real,parameter :: rtolV=0d0, atolV=0d0, atolR=epsilon(1d0)*1d4
    !real,parameter :: rtolV=epsilon(1d0)*1d2, atolV=0d0, atolR=0d0

    ! assume Newton's method succeeds

    fail = .false.

    ! save input values
    
    inV = V
    inS = S
    
    ! bracket solution and ensure initial guess lies within
    
    Vmin = min(0d0,Slock/eta)
    Vmax = max(0d0,Slock/eta)

    ! modify initial guess for nearly locked conditions
    
    call estimate_V(friction_law,rs,V,Slock,N,Psi)

    ! ensure initial guess lies within bracketed region

    if (V<=Vmin.or.Vmax<=V) V = 0.5d0*(Vmin+Vmax)

    if (info) then
       call warning("Newton's method info:")
       write(str,*) 'x=',x,' y=',y,' t=',t
       call warning(trim(str))
       write(str,*) 'input Psi=',Psi,' Slock=',Slock,' N=',N
       call warning(trim(str))
       write(str,*) 'input   V=',inV,' S=',inS
       call warning(trim(str))
       write(str,'(a,10x,a,21x,a,21x,a,24x,a,21x,a)') &
            'iteration','R(V)','Vmin','V','Vmax','dV'
       call warning(trim(str))
    end if

    ! solve S = Slock-eta*V and S = Sf(V,Psi)

    do i = 1,imax
       
       S = Slock-eta*V
       dSdV = -eta
       call strength(friction_law,rs,V,Psi,N,Sf,dSfdV)
       R = -S+Sf
       
       if (R>0d0) Vmax = V
       if (R<0d0) Vmin = V
       
       dRdV = -dSdV+dSfdV
       dV = -R/dRdV

       if (info) then
          write (str,'(i0,5(2x,f0.15))') i,R,Vmin,V,Vmax,dV
          call warning(trim(str))
          write (str,'(3(2x,f0.15))') abs(dV), atolV+rtolV*(abs(V)+abs(dV))
          call warning(trim(str))
       end if
          
       if (abs(R)<=atolR) exit
       if (i==imax) exit
       
       V = V+dV
       if (V<=Vmin.or.Vmax<=V) then
          !V = 0.5d0*(Vmin+Vmax) ! bisection on V
          if     (Vmin==0d0) then
             V = 0.1d0*Vmax
          elseif (Vmax==0d0) then
             V = 0.1d0*Vmin
          else
             V = sign(sqrt(Vmin*Vmax),Slock) ! bisection on log(V)
          end if
          ! Define a new dV which is the bracket size
          dV = Vmax-Vmin
       end if
       
       if (abs(dV)<=atolV+rtolV*(abs(V)+abs(dV))) exit
       
    end do

    if (i==imax) then ! failure

       if (info) stop

       fail = .true.

       call warning("Newton's method info:")
       write(str,*) 'x=',x,' y=',y,' t=',t
       call warning(trim(str))
       write(str,*) 'input Psi=',Psi,' Slock=',Slock,' N=',N
       call warning(trim(str))
       write(str,*) 'input   V=',inV,' S=',inS
       call warning(trim(str))
       write(str,*) 'current V=',  V,' S=',S
       call warning(trim(str))
       write(str,*) '       dV=', dV,' R=',R
       call warning(trim(str))
       write(str,*) '       rtolV=',rtolV,' atolV=',atolV,' atolR=',atolR
       call warning(trim(str))
       write(str,*) '       Sf=',Sf,' dSfdV=',dSfdV
       call warning(trim(str))
       
       call warning('Friction-elasticity solver failed','newton_solver')
       
    end if

  end subroutine newton_solver


  subroutine strength(friction_law,rs,V,Psi,N,S,dSdV)

    implicit none

    character(*),intent(in) :: friction_law
    type(ratestate_constant),intent(in) :: rs
    real,intent(in) :: V,Psi,N
    real,intent(out) :: S,dSdV

    real :: absV,f,O

    select case(friction_law)

    case('SL','FL')

       if (V==0d0) then
          S = 0d0
          dSdV = 1d10 ! simply to help Newton solver to move away from V=0
       else
          absV = abs(V)
          f = rs%a*log(absV/rs%V0)+Psi
          S = N*sign(f,V)
          dSdV = N*rs%a/absV
       end if

    case('RSL','RFL','RSF','RSL-mms','RSL-mms-nostate')

       O = 0.5d0/rs%V0*exp(Psi/rs%a)
       S = rs%a*N*arcsinh(O*V)
       dSdV = rs%a*N*O/sqrt(1d0+(V*O)**2)

    end select

  end subroutine strength


  subroutine estimate_V(friction_law,rs,V,Slock,N,Psi)

    implicit none

    character(*),intent(in) :: friction_law
    type(ratestate_constant),intent(in) :: rs
    real,intent(inout) :: V
    real,intent(in) :: Slock,N,Psi

    real :: SaN

    SaN = Slock/(rs%a*N)
    if (abs(SaN)>10d0) return ! estimate not helpful

    select case(friction_law)
    case('SL','FL')
       V = sign(rs%V0*exp(-Psi/rs%a)*exp(abs(SaN)),SaN)
    case('RSL','RFL','RSL-mms','RSL-mms-nostate')
       V = 2d0*rs%V0*exp(-Psi/rs%a)*sinh(SaN)
    end select

  end subroutine estimate_V


  function state_rate(FR,V,Psi,i,x,y,t) result(DPsi)

    use mms, only : inplane_fault_mms

    implicit none

    type(fr_type),intent(in) :: FR
    integer,intent(in) :: i
    real,intent(in) :: V,Psi,x,y,t
    real :: DPsi
    real :: Vex,Sex,Nex,Psiex,Vtex,Stex,Ntex,fssex,frsex

    real :: absV,frs,fss,fLV,Psiss,Hmss
    type(ratestate_constant) :: rs

    if (V==0d0) then
       DPsi = 0d0
       return
    end if

    absV = abs(V)

    select case(FR%friction_law)
    case('frictionless','pseudodynamic','SW')
       DPsi = 0d0
       return
    end select

    call ratestate_param(i,x,y,FR,rs)

    select case(FR%friction_law)

    case('SL')

       frs = rs%a*log(absV/rs%V0)+Psi
       fss = rs%f0-(rs%b-rs%a)*log(absV/rs%V0)
       DPsi = -(absV/rs%L)*(frs-fss)

    case('FL')

       frs = rs%a*log(absV/rs%V0)+Psi
       
       fLV = rs%f0-(rs%b-rs%a)*log(absV/rs%V0)
       fss = rs%fw+(fLV-rs%fw)/(1d0+(V/rs%Vw)**8)**0.125d0
       !fss = rs%fw+(fLV-rs%fw)/sqrt(1d0+(V/rs%Vw)**2)

       !if (absV<rs%Vw) then
       !   fss = fLV
       !else
       !   fss = rs%fw+(fLV-rs%fw)*rs%Vw/absV
       !end if

       DPsi = -(absV/rs%L)*(frs-fss)


    case('RSL','RSL-mms-nostate')

       fss = rs%f0-(rs%b-rs%a)*log(absV/rs%V0)

       !Psiss = rs%a*log(2d0*rs%V0/absV*sinh(fss/rs%a))
       !DPsi = -(absV/rs%L)*(Psi-Psiss)

       frs = rs%a*arcsinh(0.5d0*absV/rs%V0*exp(Psi/rs%a))
       DPsi = -(absV/rs%L)*(frs-fss)

    case('RSL-mms')

       ! And now we do the mms portion
       Vex   = abs(inplane_fault_mms(x,y,t,1,'V'))
       Sex   = abs(inplane_fault_mms(x,y,t,1,'S'))
       Nex   = inplane_fault_mms(x,y,t,1,'N')
       Psiex = rs%a*log(2d0*rs%V0*sinh(Sex/(rs%a*Nex))/Vex)

       ! Vtex  = inplane_fault_mms(x,y,t,1,'Vt')
       ! Stex  = inplane_fault_mms(x,y,t,1,'St')
       ! Ntex  = inplane_fault_mms(x,y,t,1,'Nt')
       ! Hmss = (-(Ntex/Nex) * (Sex/Nex) + (Stex/Nex)) * coth( (Sex / Nex) / rs%a)
       ! Hmss = Hmss  - rs%a * (Vtex / Vex)
       Hmss = 0d0

       fssex = rs%f0-(rs%b-rs%a)*log(Vex/rs%V0)
       frsex = rs%a*arcsinh(0.5d0*Vex/rs%V0*exp(Psiex/rs%a))
       frsex = Sex / Nex

       ! first we do the discrete portion
       fss = rs%f0-(rs%b-rs%a)*log(absV/rs%V0)
       frs = rs%a*arcsinh(0.5d0*absV/rs%V0*exp((Psi+Psiex)/rs%a))

       DPsi =  (Vex/rs%L)*(frsex-fssex) - (absV/rs%L)*(frs-fss)
       ! DPsi =  (Vex/rs%L)*((frsex-(absV/Vex)*frs)-(fssex-(absV/Vex)*fss))

    case('RFL')

       fLV = rs%f0-(rs%b-rs%a)*log(absV/rs%V0)
       fss = rs%fw+(fLV-rs%fw)/(1d0+(V/rs%Vw)**8)**0.125d0

!!$       if (absV<rs%Vw) then
!!$          fss = fLV
!!$       else
!!$          fss = rs%fw+(fLV-rs%fw)*rs%Vw/absV
!!$       end if

       frs = rs%a*arcsinh(0.5d0*absV/rs%V0*exp(Psi/rs%a))
       DPsi = (absV/rs%L)*(fss-frs)

    case('RSF')

       fLV = rs%f0-(rs%b-rs%a)*log(absV/rs%V0)
       fss = rs%fw+(fLV-rs%fw)/(1d0+(V/rs%Vw)**8)**0.125d0
       Psiss = rs%a*log(2d0*rs%V0/absV*sinh(fss/rs%a))
       DPsi = (absV/rs%L)*(Psiss-Psi)

    case default

       DPsi = 0d0

    end select

  end function state_rate


  subroutine rupture_front(FR,D,V,t,trup)

    use io, only : error

    implicit none

    type(fr_type),intent(in) :: FR
    real,intent(in) :: D,V,t
    real,intent(inout) :: trup

    logical :: ruptured

    select case(FR%rup_field)
    case('D')
       ruptured = (D>FR%rup_threshold)
    case('V')
       ruptured = (abs(V)>FR%rup_threshold)
    case default
       call error('Invalid field to define rupture front','rupture_front')
    end select

    if (ruptured) trup = min(trup,t)

  end subroutine rupture_front


  elemental function arcsinh(x) result(f)

    real,intent(in) :: x
    real :: f

    real :: y

    y = abs(x)
    f = log(y+sqrt(y**2+1d0))
    f = sign(f,x)

  end function arcsinh


  elemental function coth(x) result(f)

    real,intent(in) :: x
    real :: f

    f = (exp(2d0*x)+1d0)/(exp(2d0*x)-1d0)

  end function coth


  subroutine load_stress(FR,x,y,t,S0,N0)

    use utilities, only : step,boxcar,gaussian,smooth,triangle,decaying_step
    use io, only : error
    use material, only : block_material
    use geometry, only : rotate_xy2nt

    implicit none

    type(fr_type),intent(in) :: FR
    real,intent(in) :: x,y,t
    real,intent(out) :: S0,N0

    real :: r,A,B,xx,dip

    select case(FR%problem)

    case default ! specified via input parameters

       ! work in rotated coordinates

       xx =  x*cos(FR%angle)+y*sin(FR%angle) ! distance along fault
       !yy = -x*sin(FR%angle)+y*cos(FR%angle)
       r = abs((xx-FR%ld%x0)/FR%ld%R)
       N0 = 0
       S0 = 0

       select case(FR%ld%shape)
       case default
          call error('Invalid load shape','load_stress')
       case('')
          A = 0d0
       case('uniform')
          A = 1d0
       case('linear')
          A = (xx-FR%ld%x0)/FR%ld%R
       case('boxcar')
          A = boxcar  (r,1d0,1d0,0d0)
       case('triangle')
          A = triangle(r,1d0,1d0,0d0)
       case('smooth')
          A = smooth  (r,1d0,1d0,0d0)
       case('gaussian')
          A = gaussian(r,1d0,1d0,0d0)
       case('decaying_step')
          A = decaying_step(x,FR%ld%x0,FR%ld%R,1d0,0d0)
       end select

       B = load_ramp(t,FR%ld)
       
       S0 = S0 + FR%ld%S*A*B
       N0 = N0 + FR%ld%N*A*B

    case('TPV5','TPV205')

       S0 = 0d0
       N0 = 0d0
       S0 = S0+boxcar(abs(x),1.5d0,11.6d0,0d0)
       S0 = S0+boxcar(abs(x-7.5d0),1.5d0,-8d0,0d0)
       S0 = S0+boxcar(abs(x+7.5d0),1.5d0, 8d0,0d0)

    case('TPV10','TPV210')

       dip = 2d0*x
       N0 = 7.378d0*dip
       S0 = boxcar(dip-12d0,1.5d0,0.2d0+(0.76d0+0.0057d0)*N0,0.55d0*N0)

    case('TPV11')

       dip = 2d0*x
       N0 = 7.378d0*dip
       S0 = boxcar(dip-12d0,1.5d0,0.2d0+(0.57d0+0.0057d0)*N0,0.55d0*N0)

    case('TPV12alt','TPV13alt')

       dip = 2d0*x
       S0 = step(dip,13.8d0,-4.6920d0*y,0d0)
       N0 = step(dip,13.8d0,-8.5332d0*y,-16.66d0*y)

    ! alt only matters for the branch
    case('TPV14','TPV14alt')

       N0 = 0d0
       S0 = 70d0
       S0 = S0 + boxcar(abs(x+8.0d0),1.5d0,11.6d0,0d0)

    case('TPV14altBRANCH','TPV14BRANCH')

       N0 = 0d0
       S0 = 70d0

    ! alt only matters for the branch
    case('TPV15','TPV15alt')

       N0 = 0d0
       S0 = -70d0
       S0 = S0 - boxcar(abs(x+8.0d0),1.5d0,11.6d0,0d0)

    case('TPV15altBRANCH','TPV15BRANCH')

       N0 = 0d0
       S0 = -78d0


    ! alt only matters for the branch
    case('TPV14b','TPV14balt')

       N0 = 0d0
       S0 = 0d0
       S0 = boxcar(abs(x+8.0d0),1.5d0,17.6d0,0d0)

    case('TPV14baltBRANCH','TPV14bBRANCH')

       N0 = 0d0
       S0 = 0d0

    ! alt only matters for the branch
    case('TPV15b','TPV15balt')

       N0 = 0d0
       S0 = 0d0
       S0 = S0 - boxcar(abs(x+8.0d0),1.5d0,17.6d0,0d0)

    case('TPV15baltBRANCH','TPV15bBRANCH')

       N0 = 0d0
       S0 = 0d0

    case('shaw_coon') ! linear decrease in stress outside region

       r = abs(x/FR%ld%R)
       A = smooth(r,1d0,1d0,0d0)
       B = load_ramp(t,FR%ld)

       S0 = FR%ld%S*A*B
       if (abs(x)>10d0) S0 = S0-1d0*(abs(x)-10d0)
       N0 = FR%ld%N*A*B

    case('unilateral') ! drop stress to left for unilateral pulse

       xx =  x*cos(FR%angle)+y*sin(FR%angle) ! distance along fault
       r = abs((xx-FR%ld%x0)/FR%ld%R)
       A = gaussian(r,1d0,1d0,0d0)
       B = load_ramp(t,FR%ld)

       S0 = FR%ld%S*A*B+min(0d0,FR%uni_dSdx*(xx+FR%uni_x0))
       N0 = FR%ld%N*A*B

    case('gaussian_time') ! shear stress carried by incident Gaussian pulse

       A = gaussian(3d0*t-10d0,1d0,1d0,0d0)
       S0 = FR%ld%S*A
       N0 = FR%ld%N*A

    case('dip7deg')

       S0 = 2.16d0
       N0 = 4.78d0
       S0 = S0+boxcar(abs(x+60.4781d0),15d0,0.72d0,0d0)

    case('BU')

       S0 = -2.1d0*y
       N0 = -7.2d0*y

       r = abs((y-FR%ld%y0)/FR%ld%R)
       A = gaussian(r,1d0,1d0,0d0)
       B = load_ramp(t,FR%ld)

       S0 = S0+FR%ld%S*A*B
       N0 = N0+FR%ld%N*A*B

    case('SB')
       
       S0 = 22.2d0-10.09d0*y
       N0 = 37.0d0-16.82d0*y
       
       xx =  x*cos(FR%angle)+y*sin(FR%angle) ! distance along fault       
       r = abs((xx-FR%ld%x0)/FR%ld%R)
       A = gaussian(r,1d0,1d0,0d0)
       B = load_ramp(t,FR%ld)
       
       S0 = S0+FR%ld%S*A*B
       N0 = N0+FR%ld%N*A*B

    end select

  end subroutine load_stress


  function load_ramp(t,ld) result(B)

    implicit none

    real,intent(in) :: t
    type(load),intent(in) :: ld
    real :: B
    
    if (ld%T<=0d0) then
       if (t<0d0) then
          B = 0d0
       else
          B = 1d0
       end if
    else
       if (t<0d0) then
          B = 0d0
       else
          B = 1d0-exp(-t/ld%T)
       end if
       !if (t<=0d0) then
       !   B = 0d0
       !elseif (t<ld%T) then
       !   B = exp((t-ld%T)**2/(t*(t-2d0*ld%T)))
       !else
       !   B = 1d0
       !end if
    end if

  end function load_ramp


  function pseudodynamicV(i,t,pd) result(V)

    implicit none

    integer,intent(in) :: i
    real,intent(in) :: t
    type(pseudodynamic),intent(in) :: pd

    real :: V

    real :: t0

    if (t<=pd%trup(i)) then
       V = 0d0
    else
       t0 = pd%Dmax(i)/(exp(1d0)*pd%Vpeak(i)) ! peak time/exponential decay time constant
       V = pd%Dmax(i)/t0*((t-pd%trup(i))/t0)*exp(-(t-pd%trup(i))/t0) ! slip velocity
    end if

  end function pseudodynamicV


end module friction
