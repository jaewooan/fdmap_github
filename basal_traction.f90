module basal_traction

  ! module for adding basal tractions to 2D plane stress model
  ! fields F%F in medium are velocity and stress perturbations about initial equilibrium state
  
  implicit none

  ! H = layer thickness
  ! N = effective normal stress, positive in compression
  ! S = (scalar) shear strength (absolute, not change relative to initial strength)
  ! Sx0 = initial shear stress in x direction
  ! Sy0 = initial shear stress in y direction
  ! a = rate-and-state direct effect parameter (set to zero for constant friction f0)
  ! b = rate-and-state evolution effect parameter
  ! V0 = rate-and-state reference velocity
  ! f0 = rate-and-state reference friction
  ! L = rate-and-state state evolution distance
  ! Psi = state variable
  ! DPsi = d/dt of Psi
  
  type :: bt_type
     real,dimension(:,:),allocatable :: H,N,S,Sx0,Sy0,a,b,V0,f0,L,Psi,DPsi
  end type bt_type

  
contains


  subroutine init_basal_traction(BT,C,input,echo)

    use mpi_routines, only : is_master
    use mpi_routines2d, only : cartesian,allocate_array_body
    use io, only : error,write_matlab,seek_to_string  
    
    implicit none

    type(bt_type),intent(out) :: BT
    type(cartesian),intent(in) :: C
    integer,intent(in) :: input,echo

    character(256) :: filename
    integer :: stat
    real :: H,N,Sx0,Sy0,a,b,V0,f0,L,Psi

    namelist /basal_traction_list/ filename,H,N,Sx0,Sy0,a,b,V0,f0,L,Psi
    
    ! default parameters (spatially uniform)

    filename = ''
    H = 1d0
    N = 1d0
    Sx0 = 0d0
    Sy0 = 0d0
    a = 0d0
    b = 0d0
    V0 = 1d-6
    f0 = 0d0
    L = 1d0
    Psi = 1d0
    
    ! read in basal traction parameters

    rewind(input)
    read(input,nml=basal_traction_list,iostat=stat)
    if (stat>0) call error('Error in basal_traction_list','init_basal_traction')

    ! output basal traction parameters
    ! (but only if uniform <==> filename=='')
    
    if (is_master) then
       if (filename=='') then
          call write_matlab(echo,'BT.H',H)
          call write_matlab(echo,'BT.N',N)
          call write_matlab(echo,'BT.Sx0',Sx0)
          call write_matlab(echo,'BT.Sy0',Sy0)
          call write_matlab(echo,'BT.a',a)
          call write_matlab(echo,'BT.b',b)
          call write_matlab(echo,'BT.V0',V0)
          call write_matlab(echo,'BT.f0',f0)
          call write_matlab(echo,'BT.L',L)
          call write_matlab(echo,'BT.Psi',Psi)
       end if
    end if

    ! allocate memory for fields, set default value to uniform value

    call allocate_array_body(BT%H   ,C,ghost_nodes=.false.,Fval=H)
    call allocate_array_body(BT%N   ,C,ghost_nodes=.false.,Fval=N)
    call allocate_array_body(BT%S   ,C,ghost_nodes=.false.)
    call allocate_array_body(BT%Sx0 ,C,ghost_nodes=.false.,Fval=Sx0)
    call allocate_array_body(BT%Sy0 ,C,ghost_nodes=.false.,Fval=Sy0)
    call allocate_array_body(BT%a   ,C,ghost_nodes=.false.,Fval=a)
    call allocate_array_body(BT%b   ,C,ghost_nodes=.false.,Fval=b)
    call allocate_array_body(BT%V0  ,C,ghost_nodes=.false.,Fval=V0)
    call allocate_array_body(BT%f0  ,C,ghost_nodes=.false.,Fval=f0)
    call allocate_array_body(BT%L   ,C,ghost_nodes=.false.,Fval=L)
    call allocate_array_body(BT%Psi ,C,ghost_nodes=.false.,Fval=Psi)
    call allocate_array_body(BT%DPsi,C,ghost_nodes=.false.)

    ! initialize fields from file

    if (filename/='') call basal_traction_IO('read',filename,C,BT)
    
  end subroutine init_basal_traction


  subroutine basal_traction_IO(operation,filename,C,BT)

    use io, only : file_distributed,open_file_distributed, &
         read_file_distributed,write_file_distributed,close_file_distributed
    use mpi_routines2d, only : cartesian
    use mpi_routines, only : pw,is_master
    use mpi

    implicit none

    character(*),intent(in) :: operation,filename
    type(cartesian),intent(in) :: C
    type(bt_type),intent(inout) :: BT

    type(file_distributed) :: fh
    integer :: ierr

    if (operation=='delete') then
       if (is_master) call MPI_file_delete(filename,MPI_INFO_NULL,ierr)
       return
    end if

    call open_file_distributed(fh,filename,operation,C%c2d%comm,C%c2d%array_w,pw)

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    select case(operation)
    case('read')
       call  read_file_distributed(fh,BT%H  (C%mx:C%px,C%my:C%py))
       call  read_file_distributed(fh,BT%N  (C%mx:C%px,C%my:C%py))
       call  read_file_distributed(fh,BT%Sx0(C%mx:C%px,C%my:C%py))
       call  read_file_distributed(fh,BT%Sy0(C%mx:C%px,C%my:C%py))
       call  read_file_distributed(fh,BT%a  (C%mx:C%px,C%my:C%py))
       call  read_file_distributed(fh,BT%b  (C%mx:C%px,C%my:C%py))
       call  read_file_distributed(fh,BT%V0 (C%mx:C%px,C%my:C%py))
       call  read_file_distributed(fh,BT%f0 (C%mx:C%px,C%my:C%py))
       call  read_file_distributed(fh,BT%L  (C%mx:C%px,C%my:C%py))
       call  read_file_distributed(fh,BT%Psi(C%mx:C%px,C%my:C%py))
    case('write')
       call write_file_distributed(fh,BT%H  (C%mx:C%px,C%my:C%py))
       call write_file_distributed(fh,BT%N  (C%mx:C%px,C%my:C%py))
       call write_file_distributed(fh,BT%Sx0(C%mx:C%px,C%my:C%py))
       call write_file_distributed(fh,BT%Sy0(C%mx:C%px,C%my:C%py))
       call write_file_distributed(fh,BT%a  (C%mx:C%px,C%my:C%py))
       call write_file_distributed(fh,BT%b  (C%mx:C%px,C%my:C%py))
       call write_file_distributed(fh,BT%V0 (C%mx:C%px,C%my:C%py))
       call write_file_distributed(fh,BT%f0 (C%mx:C%px,C%my:C%py))
       call write_file_distributed(fh,BT%L  (C%mx:C%px,C%my:C%py))
       call write_file_distributed(fh,BT%Psi(C%mx:C%px,C%my:C%py))
    end select

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    call close_file_distributed(fh)

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

  end subroutine basal_traction_IO


  subroutine destroy_basal_traction(BT)

    implicit none

    type(bt_type),intent(inout) :: BT

    if (allocated(BT%H   )) deallocate(BT%H   )
    if (allocated(BT%N   )) deallocate(BT%N   )
    if (allocated(BT%S   )) deallocate(BT%S   )
    if (allocated(BT%Sx0 )) deallocate(BT%Sx0 )
    if (allocated(BT%Sy0 )) deallocate(BT%Sy0 )
    if (allocated(BT%a   )) deallocate(BT%a   )
    if (allocated(BT%b   )) deallocate(BT%b   )
    if (allocated(BT%V0  )) deallocate(BT%V0  )
    if (allocated(BT%f0  )) deallocate(BT%f0  )
    if (allocated(BT%L   )) deallocate(BT%L   )
    if (allocated(BT%Psi )) deallocate(BT%Psi )
    if (allocated(BT%DPsi)) deallocate(BT%DPsi)

  end subroutine destroy_basal_traction


  subroutine set_basal_traction(C,F,BT,mode)
    
    use mpi_routines2d, only : cartesian
    use fields, only : fields_type
    use io, only : error
    
    implicit none

    type(cartesian),intent(in) :: C
    type(fields_type),intent(inout) :: F
    type(bt_type),intent(inout) :: BT
    integer,intent(in) :: mode

    integer :: i,j
    real :: frs,fss,Sx,Sy,V
    
    select case(mode)
    case(2)
       do j = C%my,C%py
          do i = C%mx,C%px
             V = sqrt(F%F(i,j,1)**2+F%F(i,j,2)**2) ! slip velocity
             if (BT%a(i,j)==0d0.or.V==0d0) then
                ! constant friction
                frs = BT%f0(i,j)
                fss = BT%f0(i,j)
             else
                ! rate-and-state friction
                frs = BT%a(i,j)*arcsinh(0.5d0*V/BT%V0(i,j)*exp(BT%Psi(i,j)/BT%a(i,j)))
                fss = BT%f0(i,j)-(BT%b(i,j)-BT%a(i,j))*log(V/BT%V0(i,j))
             end if
             BT%S(i,j) = BT%N(i,j)*frs ! shear strength
             ! align strength with slip velocity (Vx=F%F(i,j,1),Vy=F%F(i,j,2))
             if (V==0d0) then
                Sx = 0d0
                Sy = 0d0
             else
                Sx = BT%S(i,j)*F%F(i,j,1)/V
                Sy = BT%S(i,j)*F%F(i,j,2)/V
             end if
             ! add basal traction (but only perturbation from initial traction)
             F%DF(i,j,1) = F%DF(i,j,1)-(Sx-BT%Sx0(i,j))/BT%H(i,j)
             F%DF(i,j,2) = F%DF(i,j,2)-(Sy-BT%Sy0(i,j))/BT%H(i,j)
             ! state evolution
             BT%DPsi(i,j) = BT%DPsi(i,j)-V/BT%L(i,j)*(frs-fss)
          end do
       end do
    case(3)
       call error('basal traction not valid in mode 3','set_basal_traction')
    end select

  end subroutine set_basal_traction


  subroutine scale_rates_basal_traction(BT,A)

    implicit none

    type(bt_type),intent(inout) :: BT
    real,intent(in) :: A
    
    BT%DPsi = A*BT%DPsi
    
  end subroutine scale_rates_basal_traction


  subroutine update_fields_basal_traction(BT,dt)

    implicit none

    type(bt_type),intent(inout) :: BT
    real,intent(in) :: dt
    
    BT%Psi = BT%Psi+dt*BT%DPsi
    
  end subroutine update_fields_basal_traction

  
  elemental function arcsinh(x) result(f)

    real,intent(in) :: x
    real :: f

    real :: y

    y = abs(x)
    f = log(y+sqrt(y**2+1d0))
    f = sign(f,x)

  end function arcsinh

  
end module basal_traction
