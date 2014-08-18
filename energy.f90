module energy

  implicit none


contains


  subroutine energy_interior(D)
    
    use domain, only : domain_type

    implicit none

    type(domain_type),intent(inout) :: D

    real :: E
    integer :: i,j,ierr

    if (.not.D%F%energy_balance) return

    ! loop over x and y, compute inner product

    E = 0d0
    do i = 1,D%nblocks
      call energy_block(D%C,D%B(i)%G,D%B(i)%F,D%B(i)%M,D%F,D%mode,D%G)
      E = E + D%B(i)%F%Etot
    end do
    ! call energy_interior(D%C,D%F,D%mode,0.5d0*(D%G%hmin+D%G%hmax))
    D%F%Etot = E

  end subroutine energy_interior


  subroutine energy_block(C,B,BF,M,F,mode,G)
    
    use grid, only : block_grid, grid_type
    use mpi_routines2d, only : cartesian
    use fields, only : block_fields,fields_type
    use material, only : block_material
    use mpi_routines, only : MPI_REAL_PW
    use mpi

    implicit none

    type(cartesian),intent(in) :: C
    type(grid_type),intent(in) :: G
    type(block_grid),intent(in) :: B
    type(block_fields),intent(inout) :: BF
    type(block_material),intent(in) :: M
    type(fields_type),intent(in) :: F
    integer,intent(in) :: mode

    integer :: i,j,ierr
    real :: E
    real :: v1,v2,v3,v4

    ! If we aren't doing the calculation return
    if (.not.BF%energy_balance) return

    ! initialize the energy
    E = 0d0

    ! loop over the block
    if (.not.B%skip) then

      ! loop over x and y, compute inner product
      select case(mode)
      case(2)
         v1 = 0.5d0*M%rho
         v2 = 0.5d0*M%G
         v3 = 0.25d0*M%G
         v4 =-0.25d0*M%G*(M%lambda/(3.0d0*M%lambda + 2.0d0*M%G));
         do j = B%my,B%py
            do i = B%mx,B%px
               E = E+F%Hx(i)*F%Hy(j)*G%J(i,j)*( &
                      v1*(F%F(i,j,1)**2+F%F(i,j,2)**2) &
                    + v2*F%F(i,j,4)**2 &
                    + v3*(F%F(i,j,3)**2 + F%F(i,j,5)**2 + F%F(i,j,6)**2) &
                    + v4*(F%F(i,j,3)    + F%F(i,j,5)    + F%F(i,j,6)   )**2)
            end do
         end do
      case(3)
         do j = B%my,B%py
            do i = B%mx,B%px
               E = E+0.5d0*G%J(i,j)*(M%rho*F%F(i,j,1)**2+ &
                    sum(F%F(i,j,2:3)**2)/M%G)*F%Hx(i)*F%Hy(j)
            end do
         end do
      end select
    endif

    ! call MPI_Allreduce(E,BF%Etot,1,MPI_REAL_PW,MPI_SUM,B%comm,ierr)
    call MPI_Allreduce(E,BF%Etot,1,MPI_REAL_PW,MPI_SUM,C%c2d%comm,ierr)

  end subroutine energy_block


  subroutine scale_rates_energy(B,BF,A)

    use grid, only : block_grid
    use fields, only : block_fields

    type(block_grid),intent(in) :: B
    type(block_fields),intent(inout) :: BF
    real,intent(in) :: A

    if (.not.BF%energy_balance) return

    if (B%sideL) BF%bndFL%DE = A*BF%bndFL%DE
    if (B%sideR) BF%bndFR%DE = A*BF%bndFR%DE
    if (B%sideB) BF%bndFB%DE = A*BF%bndFB%DE
    if (B%sideT) BF%bndFT%DE = A*BF%bndFT%DE

  end subroutine scale_rates_energy


  subroutine update_energy(B,BF,dt)

    use grid, only : block_grid
    use fields, only : block_fields

    implicit none

    type(block_grid),intent(in) :: B
    type(block_fields),intent(inout) :: BF
    real,intent(in) :: dt

    if (.not.BF%energy_balance) return

    BF%E = BF%E+dt*BF%DE

    if (B%sideL) BF%bndFL%E = BF%bndFL%E+dt*BF%bndFL%DE
    if (B%sideR) BF%bndFR%E = BF%bndFR%E+dt*BF%bndFR%DE
    if (B%sideB) BF%bndFB%E = BF%bndFB%E+dt*BF%bndFB%DE
    if (B%sideT) BF%bndFT%E = BF%bndFT%E+dt*BF%bndFT%DE

  end subroutine update_energy


  subroutine set_rates_energy(B,BF,mode)

    use grid, only : block_grid
    use fields, only : block_fields
    use mpi_routines, only : MPI_REAL_PW
    use mpi

    implicit none
    
    type(block_grid),intent(in) :: B
    type(block_fields),intent(inout) :: BF
    integer,intent(in) :: mode

    real :: DE,DEL,DER,DEB,DET,DELRBT(4),DEglobal(4)
    integer :: i,j,ierr

    if (B%skip.or.(.not.BF%energy_balance)) return

    ! calculate energy rates on block boundaries

    if (B%sideL) call energy_rates_boundary(B%bndL,BF%bndFL,B%my,B%py,mode)
    if (B%sideR) call energy_rates_boundary(B%bndR,BF%bndFR,B%my,B%py,mode)
    if (B%sideB) call energy_rates_boundary(B%bndB,BF%bndFB,B%my,B%py,mode)
    if (B%sideT) call energy_rates_boundary(B%bndT,BF%bndFT,B%my,B%py,mode)

    ! integrate energy flux across boundary: multiply by norm
    ! sum all points held by process, global sum across all processes

    if (B%sideL) then
       DE = 0d0
       do j = B%my,B%py
          DE = DE+BF%bndFL%DE(j)*BF%Hy(j)*B%xyrL(j)
       end do
       call MPI_Allreduce(DE,DEL,1,MPI_REAL_PW,MPI_SUM,B%commy,ierr)
    else
       DEL = huge(DEL)
    end if
 
    if (B%sideR) then
       DE = 0d0
       do j = B%my,B%py
          DE = DE+BF%bndFR%DE(j)*BF%Hy(j)*B%xyrR(j)
       end do
       call MPI_Allreduce(DE,DER,1,MPI_REAL_PW,MPI_SUM,B%commy,ierr)
    else
       DER = huge(DER)
    end if

    if (B%sideB) then
       DE = 0d0
       do i = B%mx,B%px
          DE = DE+BF%bndFB%DE(i)*BF%Hx(i)*B%xyqB(i)
       end do
       call MPI_Allreduce(DE,DEB,1,MPI_REAL_PW,MPI_SUM,B%commx,ierr)
    else
       DEB = huge(DEB)
    end if

    if (B%sideT) then
       DE = 0d0
       do i = B%mx,B%px
          DE = DE+BF%bndFT%DE(i)*BF%Hx(i)*B%xyqT(i)
       end do
       call MPI_Allreduce(DE,DET,1,MPI_REAL_PW,MPI_SUM,B%commx,ierr)
    else
       DET = huge(DET)
    end if

    DELRBT = (/ DEL,DER,DEB,DET /)
    call MPI_Allreduce(DELRBT,DEglobal,4,MPI_REAL_PW,MPI_MIN,B%comm,ierr)

    BF%DE = DEglobal

  end subroutine set_rates_energy


  subroutine energy_rates_boundary(bndC,bndF,m,p,mode)

    use geometry, only : curve
    use fields, only : bnd_fields

    implicit none
    
    type(curve),intent(in) :: bndC
    type(bnd_fields),intent(inout) :: bndF
    integer,intent(in) :: m,p,mode

    integer :: i
    
    select case(mode)
    case(3)
       do i = m,p
          ! DE = DE+... (use subroutines below, but will need to rotate fields, etc.)
          ! use fields bndF%DE(i),bndF%Fhat(i,:),bndF%F(i,:),bndF%F0(i,:),bndC%n(i,:),bndF%M(i,1:3)
       end do
    case(2)
    end select

  end subroutine energy_rates_boundary


  subroutine energy_rate_mode3(DE,snz,vz,snzFD,vzFD,Zs,Zsi)

    implicit none

    real,intent(inout) :: DE
    real,intent(in) :: snz,vz,snzFD,vzFD,Zs,Zsi

    DE = DE+snz*vz-0.25d0*Zsi*((snzFD-snz)+Zs*(vzFD-vz))**2

  end subroutine energy_rate_mode3

  
  subroutine energy_rate_mode2(DE,snt,snn,vt,vn,sntFD,snnFD,vtFD,vnFD,Zs,Zsi,Zp,Zpi)
    
    implicit none
    
    real,intent(inout) :: DE
    real,intent(in) :: snt,snn,vt,vn,sntFD,snnFD,vtFD,vnFD,Zs,Zsi,Zp,Zpi
    
    DE = DE+snt*vt+vn*snn- &
         0.25d0*Zsi*((sntFD-snt)+Zs*(vtFD-vt))**2- &
         0.25d0*Zpi*((snnFD-snn)+Zp*(vnFD-vn))**2

  end subroutine energy_rate_mode2


end module energy
