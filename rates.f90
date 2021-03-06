module rates

  implicit none

contains


  subroutine rates_interior_mode3(lbndFx,lbndFy,lbndDFx,lbndDFy, &
       F,DF,rho,G,xq,xr,yq,yr,JacInv,mx,px,my,py)

    use fd_coeff, only : DI

    implicit none

    integer,intent(in) :: lbndFx,lbndFy,lbndDFx,lbndDFy,mx,px,my,py
    real,dimension(lbndFx :,lbndFy :,1:),intent(in) :: F
    real,dimension(lbndDFx:,lbndDFy:,1:),intent(inout) :: DF
    real,intent(in) :: rho,G
    real,dimension(lbndFx:,lbndFy:),intent(in) :: xq,xr,yq,yr,JacInv

    integer :: i,j,s
    real,dimension(3) :: a,b

    a=0d0; b=0d0
    do i = 1,ubound(DI,1)
       a(i) = DI(i)/rho
       b(i) = DI(i)*G
    end do

    do j = my,py
       do i = mx,px
          do s = 1,3

             DF(i,j,1) = DF(i,j,1)+( &
                  0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,2)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,3)- &
                  0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,2)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,3)- &
                  0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,2)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,3)+ &
                  0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,2)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,3) &
                  )*a(s)*JacInv(i,j)

             DF(i,j,2) = DF(i,j,2)+( &
                  0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,1)-0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,1)- &
                  0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,1)+0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,1) &
                  )*b(s)*JacInv(i,j)

             DF(i,j,3) = DF(i,j,3)+( &
                 -0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,1)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,1)+ &
                  0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,1)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,1) &
                  )*b(s)*JacInv(i,j)


          end do
       end do
    end do

  end subroutine rates_interior_mode3


  subroutine rates_boundary_mode3(lbndFx,lbndFy,lbndDFx,lbndDFy, &
       F,DF,rho,G,xq,xr,yq,yr,JacInv,mx,px,my,py,bnd)

    use fd_coeff, only : DL,DR,DI,pI,nbnd

    implicit none

    integer,intent(in) :: lbndFx,lbndFy,lbndDFx,lbndDFy,mx,px,my,py
    real,dimension(lbndFx :,lbndFy :,1:),intent(in) :: F
    real,dimension(lbndDFx:,lbndDFy:,1:),intent(inout) :: DF
    real,intent(in) :: rho,G
    real,dimension(lbndFx:,lbndFy:),intent(in) :: xq,xr,yq,yr,JacInv
    character(*),intent(in) :: bnd

    integer :: i,j,s

    select case(bnd)

    case('L')

       do j = my,py
          do i = mx,px

             do s = 0,nbnd-1 ! one-sided for q
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,2)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,3) &
                     )*DL(s,i-mx)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,1) &
                     )*DL(s,i-mx)*G*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     -0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,1) &
                     )*DL(s,i-mx)*G*JacInv(i,j)
                
             end do

             do s = 1,pI ! central for r
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,2)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,3)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,2)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,3) &
                     )*DI(s)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,1)+0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,1) &
                     )*DI(s)*G*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,1)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,1) &
                     )*DI(s)*G*JacInv(i,j)

             end do

          end do
       end do

    case('R')
       
       do j = my,py
          do i = mx,px

             do s = 0,nbnd-1 ! one-sided for q
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,2)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,3) &
                     )*DR(-s,i-px)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,1) &
                     )*DR(-s,i-px)*G*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     -0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,1) &
                     )*DR(-s,i-px)*G*JacInv(i,j)

             end do

             do s = 1,pI ! central for r
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,2)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,3)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,2)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,3) &
                     )*DI(s)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,1)+0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,1) &
                     )*DI(s)*G*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,1)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,1) &
                     )*DI(s)*G*JacInv(i,j)
                
             end do

          end do
       end do

    case('B')
       
       do j = my,py
          do i = mx,px

             do s = 0,nbnd-1 ! one-sided for r
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,2)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,3) &
                     )*DL(s,j-my)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,1) &
                     )*DL(s,j-my)*G*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,1) &
                     )*DL(s,j-my)*G*JacInv(i,j)
                
             end do

             do s = 1,pI ! central for q
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,2)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,3)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,2)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,3) &
                     )*DI(s)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,1)-0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,1) &
                     )*DI(s)*G*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     -0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,1)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,1) &
                     )*DI(s)*G*JacInv(i,j)
                
             end do

          end do
       end do
       
    case('T')
    
       do j = my,py
          do i = mx,px

             do s = 0,nbnd-1 ! one-sided for r
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,2)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,3) &
                     )*DR(-s,j-py)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,1) &
                     )*DR(-s,j-py)*G*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,1) &
                     )*DR(-s,j-py)*G*JacInv(i,j)
                
                
             end do

             do s = 1,pI ! central for q
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,2)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,3)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,2)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,3) &
                     )*DI(s)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,1)-0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,1) &
                     )*DI(s)*G*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     -0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,1)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,1) &
                     )*DI(s)*G*JacInv(i,j)
                
                
             end do

          end do
       end do
       
    case('LB')

       do j = my,py
          do i = mx,px
             do s = 0,nbnd-1 ! one-sided for q and r
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,2)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,3) &
                     )*DL(s,i-mx)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,1) &
                     )*DL(s,i-mx)*G*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     -0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,1) &
                     )*DL(s,i-mx)*G*JacInv(i,j)
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,2)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,3) &
                     )*DL(s,j-my)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,1) &
                     )*DL(s,j-my)*G*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,1) &
                     )*DL(s,j-my)*G*JacInv(i,j)
                
             end do
          end do
       end do

    case('RB')

       do j = my,py
          do i = mx,px
             do s = 0,nbnd-1 ! one-sided for q and r

                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,2)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,3) &
                     )*DR(-s,i-px)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,1) &
                     )*DR(-s,i-px)*G*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     -0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,1) &
                     )*DR(-s,i-px)*G*JacInv(i,j)
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,2)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,3) &
                     )*DL(s,j-my)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,1) &
                     )*DL(s,j-my)*G*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,1) &
                     )*DL(s,j-my)*G*JacInv(i,j)

             end do
          end do
       end do

    case('LT')

       do j = my,py
          do i = mx,px
             do s = 0,nbnd-1 ! one-sided for q and r
              
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,2)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,3) &
                     )*DL(s,i-mx)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,1) &
                     )*DL(s,i-mx)*G*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     -0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,1) &
                     )*DL(s,i-mx)*G*JacInv(i,j)

                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,2)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,3) &
                     )*DR(-s,j-py)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,1) &
                     )*DR(-s,j-py)*G*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,1) &
                     )*DR(-s,j-py)*G*JacInv(i,j)
                
             end do
          end do
       end do

    case('RT')
    
       do j = my,py
          do i = mx,px
             do s = 0,nbnd-1 ! one-sided for q and r
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,2)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,3) &
                     )*DR(-s,i-px)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,1) &
                     )*DR(-s,i-px)*G*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     -0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,1) &
                     )*DR(-s,i-px)*G*JacInv(i,j)
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,2)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,3) &
                     )*DR(-s,j-py)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,1) &
                     )*DR(-s,j-py)*G*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,1) &
                     )*DR(-s,j-py)*G*JacInv(i,j)
                
             end do
          end do
       end do

    end select

  end subroutine rates_boundary_mode3


  subroutine rates_interior_mode2(lbndFx,lbndFy,lbndDFx,lbndDFy, &
       F,DF,rho,G,M,gamma,xq,xr,yq,yr,JacInv,mx,px,my,py)

    use fd_coeff, only : DI

    implicit none

    integer,intent(in) :: lbndFx,lbndFy,lbndDFx,lbndDFy,mx,px,my,py
    real,dimension(lbndFx :,lbndFy :,1:),intent(in) :: F
    real,dimension(lbndDFx:,lbndDFy:,1:),intent(inout) :: DF
    real,intent(in) :: rho,G,M,gamma
    real,dimension(lbndFx:,lbndFy:),intent(in) :: xq,xr,yq,yr,JacInv

    integer :: i,j,s
    real,dimension(3) :: a,b,c

    a=0d0; b=0d0; c=0d0
    do i = 1,ubound(DI,1)
       a(i) = DI(i)/rho
       b(i) = DI(i)*G
       c(i) = DI(i)*M
    end do

    do j = my,py
       do i = mx,px
          do s = 1,3

             DF(i,j,1) = DF(i,j,1)+( &
                  0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,3)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,4)- &
                  0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,3)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,4)- &
                  0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,3)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,4)+ &
                  0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,3)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,4) &
                  )*a(s)*JacInv(i,j)

             DF(i,j,2) = DF(i,j,2)+( &
                  0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,4)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,5)- &
                  0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,4)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,5)- &
                  0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,4)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,5)+ &
                  0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,4)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,5) &
                  )*a(s)*JacInv(i,j)

             DF(i,j,3) = DF(i,j,3)+( &
                  0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,1)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,2)*gamma- &
                  0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,1)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,2)*gamma- &
                  0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,1)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,2)*gamma+ &
                  0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,1)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,2)*gamma &
                  )*c(s)*JacInv(i,j)

             DF(i,j,4) = DF(i,j,4)+( &
                  0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,2)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,1)- &
                  0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,2)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,1)- &
                  0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,2)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,1)+ &
                  0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,2)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,1) &
                  )*b(s)*JacInv(i,j)

             DF(i,j,5) = DF(i,j,5)+( &
                  0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,1)*gamma-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,2)- &
                  0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,1)*gamma+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,2)- &
                  0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,1)*gamma+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,2)+ &
                  0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,1)*gamma-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,2) &
                  )*c(s)*JacInv(i,j)

          end do
       end do
    end do

  end subroutine rates_interior_mode2


  subroutine rates_boundary_mode2(lbndFx,lbndFy,lbndDFx,lbndDFy, &
       F,DF,rho,G,M,gamma,xq,xr,yq,yr,JacInv,mx,px,my,py,bnd)

    use fd_coeff, only : DL,DR,DI,pI,nbnd

    implicit none

    integer,intent(in) :: lbndFx,lbndFy,lbndDFx,lbndDFy,mx,px,my,py
    real,dimension(lbndFx :,lbndFy :,1:),intent(in) :: F
    real,dimension(lbndDFx:,lbndDFy:,1:),intent(inout) :: DF
    real,intent(in) :: rho,G,M,gamma
    real,dimension(lbndFx:,lbndFy:),intent(in) :: xq,xr,yq,yr,JacInv
    character(*),intent(in) :: bnd

    integer :: i,j,s

    select case(bnd)

    case('L')

       do j = my,py
          do i = mx,px

             do s = 0,nbnd-1 ! one-sided for q
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,3)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,4) &
                     )*DL(s,i-mx)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,4)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,5) &
                     )*DL(s,i-mx)/rho*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,1)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,2)*gamma &
                     )*DL(s,i-mx)*M*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,2)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,1) &
                     )*DL(s,i-mx)*G*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,1)*gamma-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,2) &
                     )*DL(s,i-mx)*M*JacInv(i,j)

             end do

             do s = 1,pI ! central for r
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,3)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,4)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,3)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,4) &
                     )*DI(s)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,4)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,5)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,4)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,5) &
                     )*DI(s)/rho*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,1)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,2)*gamma+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,1)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,2)*gamma &
                     )*DI(s)*M*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,2)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,1)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,2)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,1) &
                     )*DI(s)*G*JacInv(i,j)

                DF(i,j,5) = DF(i,j,5)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,1)*gamma+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,2)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,1)*gamma-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,2) &
                     )*DI(s)*M*JacInv(i,j)

             end do

          end do
       end do

    case('R')
       
       do j = my,py
          do i = mx,px

             do s = 0,nbnd-1 ! one-sided for q
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,3)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,4) &
                     )*DR(-s,i-px)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,4)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,5) &
                     )*DR(-s,i-px)/rho*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,1)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,2)*gamma &
                     )*DR(-s,i-px)*M*JacInv(i,j)

                DF(i,j,4) = DF(i,j,4)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,2)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,1) &
                     )*DR(-s,i-px)*G*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,1)*gamma-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,2) &
                     )*DR(-s,i-px)*M*JacInv(i,j)

             end do

             do s = 1,pI ! central for r
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,3)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,4)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,3)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,4) &
                     )*DI(s)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,4)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,5)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,4)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,5) &
                     )*DI(s)/rho*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,1)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,2)*gamma+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,1)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,2)*gamma &
                     )*DI(s)*M*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,2)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,1)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,2)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,1) &
                     )*DI(s)*G*JacInv(i,j)

                DF(i,j,5) = DF(i,j,5)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,1)*gamma+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,2)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,1)*gamma-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,2) &
                     )*DI(s)*M*JacInv(i,j)

             end do

          end do
       end do

    case('B')
       
       do j = my,py
          do i = mx,px

             do s = 0,nbnd-1 ! one-sided for r
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,3)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,4) &
                     )*DL(s,j-my)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,4)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,5) &
                     )*DL(s,j-my)/rho*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,1)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,2)*gamma &
                     )*DL(s,j-my)*M*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,2)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,1) &
                     )*DL(s,j-my)*G*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,1)*gamma+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,2) &
                     )*DL(s,j-my)*M*JacInv(i,j)
                
             end do

             do s = 1,pI ! central for q
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,3)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,4)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,3)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,4) &
                     )*DI(s)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,4)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,5)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,4)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,5) &
                     )*DI(s)/rho*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,1)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,2)*gamma- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,1)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,2)*gamma &
                     )*DI(s)*M*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,2)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,1)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,2)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,1) &
                     )*DI(s)*G*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,1)*gamma-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,2)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,1)*gamma+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,2) &
                     )*DI(s)*M*JacInv(i,j)
                
             end do

          end do
       end do
       
    case('T')
    
       do j = my,py
          do i = mx,px

             do s = 0,nbnd-1 ! one-sided for r
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,3)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,4) &
                     )*DR(-s,j-py)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,4)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,5) &
                     )*DR(-s,j-py)/rho*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,1)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,2)*gamma &
                     )*DR(-s,j-py)*M*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,2)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,1) &
                     )*DR(-s,j-py)*G*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,1)*gamma+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,2) &
                     )*DR(-s,j-py)*M*JacInv(i,j)
                
             end do

             do s = 1,pI ! central for q
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,3)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,4)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,3)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,4) &
                     )*DI(s)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,4)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,5)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,4)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,5) &
                     )*DI(s)/rho*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,1)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,2)*gamma- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,1)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,2)*gamma &
                     )*DI(s)*M*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,2)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,1)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,2)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,1) &
                     )*DI(s)*G*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,1)*gamma-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,2)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,1)*gamma+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,2) &
                     )*DI(s)*M*JacInv(i,j)
                
             end do

          end do
       end do
       
    case('LB')

       do j = my,py
          do i = mx,px
             do s = 0,nbnd-1 ! one-sided for q and r
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,3)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,4) &
                     )*DL(s,i-mx)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,4)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,5) &
                     )*DL(s,i-mx)/rho*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,1)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,2)*gamma &
                     )*DL(s,i-mx)*M*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,2)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,1) &
                     )*DL(s,i-mx)*G*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,1)*gamma-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,2) &
                     )*DL(s,i-mx)*M*JacInv(i,j)
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,3)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,4) &
                     )*DL(s,j-my)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,4)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,5) &
                     )*DL(s,j-my)/rho*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,1)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,2)*gamma &
                     )*DL(s,j-my)*M*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,2)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,1) &
                     )*DL(s,j-my)*G*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,1)*gamma+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,2) &
                     )*DL(s,j-my)*M*JacInv(i,j)
                
             end do
          end do
       end do

    case('RB')

       do j = my,py
          do i = mx,px
             do s = 0,nbnd-1 ! one-sided for q and r

                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,3)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,4) &
                     )*DR(-s,i-px)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,4)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,5) &
                     )*DR(-s,i-px)/rho*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,1)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,2)*gamma &
                     )*DR(-s,i-px)*M*JacInv(i,j)

                DF(i,j,4) = DF(i,j,4)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,2)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,1) &
                     )*DR(-s,i-px)*G*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,1)*gamma-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,2) &
                     )*DR(-s,i-px)*M*JacInv(i,j)

                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,3)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,4) &
                     )*DL(s,j-my)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,4)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,5) &
                     )*DL(s,j-my)/rho*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,1)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,2)*gamma &
                     )*DL(s,j-my)*M*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,2)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,1) &
                     )*DL(s,j-my)*G*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,1)*gamma+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,2) &
                     )*DL(s,j-my)*M*JacInv(i,j)
                
             end do
          end do
       end do

    case('LT')

       do j = my,py
          do i = mx,px
             do s = 0,nbnd-1 ! one-sided for q and r
              
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,3)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,4) &
                     )*DL(s,i-mx)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,4)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,5) &
                     )*DL(s,i-mx)/rho*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,1)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,2)*gamma &
                     )*DL(s,i-mx)*M*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,2)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,1) &
                     )*DL(s,i-mx)*G*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,1)*gamma-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,2) &
                     )*DL(s,i-mx)*M*JacInv(i,j)
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,3)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,4) &
                     )*DR(-s,j-py)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,4)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,5) &
                     )*DR(-s,j-py)/rho*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,1)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,2)*gamma &
                     )*DR(-s,j-py)*M*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,2)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,1) &
                     )*DR(-s,j-py)*G*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,1)*gamma+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,2) &
                     )*DR(-s,j-py)*M*JacInv(i,j)
                
             end do
          end do
       end do

    case('RT')
    
       do j = my,py
          do i = mx,px
             do s = 0,nbnd-1 ! one-sided for q and r
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,3)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,4) &
                     )*DR(-s,i-px)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,4)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,5) &
                     )*DR(-s,i-px)/rho*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,1)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,2)*gamma &
                     )*DR(-s,i-px)*M*JacInv(i,j)

                DF(i,j,4) = DF(i,j,4)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,2)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,1) &
                     )*DR(-s,i-px)*G*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,1)*gamma-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,2) &
                     )*DR(-s,i-px)*M*JacInv(i,j)

                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,3)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,4) &
                     )*DR(-s,j-py)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,4)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,5) &
                     )*DR(-s,j-py)/rho*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,1)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,2)*gamma &
                     )*DR(-s,j-py)*M*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,2)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,1) &
                     )*DR(-s,j-py)*G*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,1)*gamma+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,2) &
                     )*DR(-s,j-py)*M*JacInv(i,j)
                
             end do
          end do
       end do

    end select

  end subroutine rates_boundary_mode2


  subroutine rates_szz(lbndDFx,lbndDFy,DF,nu,mx,px,my,py)

    implicit none

    integer,intent(in) :: lbndDFx,lbndDFy,mx,px,my,py
    real,dimension(lbndDFx:,lbndDFy:,1:),intent(inout) :: DF
    real,intent(in) :: nu

    integer :: i,j

    do j = my,py
       do i = mx,px
          DF(i,j,6) = nu*(DF(i,j,3)+DF(i,j,5))
       end do
    end do

  end subroutine rates_szz


  subroutine rates_dissipation_interior_mode3(lbndFx,lbndFy,lbndDFx,lbndDFy, &
       F,DF,rho,G,Cdiss,xq,xr,yq,yr,Jac,JacInv,mx,px,my,py)

    use fd_coeff, only : DI,AI

    implicit none

    integer,intent(in) :: lbndFx,lbndFy,lbndDFx,lbndDFy,mx,px,my,py
    real,dimension(lbndFx :,lbndFy :,1:),intent(in) :: F
    real,dimension(lbndDFx:,lbndDFy:,1:),intent(inout) :: DF
    real,intent(in) :: rho,G,Cdiss
    real,dimension(lbndFx:,lbndFy:),intent(in) :: xq,xr,yq,yr,Jac,JacInv

    integer :: i,j,s
    real :: a(1:3),b(1:3),c(0:3)

    a=0d0; b=0d0; c=0d0
    do i = 1,ubound(DI,1)
       a(i) = DI(i)/rho
       b(i) = DI(i)*G
       c(i) = AI(i)*Cdiss
    end do
    c(0) = 2d0*AI(0)*Cdiss ! coefficient used twice (both for q and r derivatives)

    do j = my,py
       do i = mx,px

          DF(i,j,1) = DF(i,j,1)+c(0)*F(i,j,1)
          DF(i,j,2) = DF(i,j,2)+c(0)*F(i,j,2)
          DF(i,j,3) = DF(i,j,3)+c(0)*F(i,j,3)

          do s = 1,3

             DF(i,j,1) = DF(i,j,1)+( &
                  0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,2)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,3)- &
                  0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,2)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,3)- &
                  0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,2)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,3)+ &
                  0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,2)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,3) &
                  )*a(s)*JacInv(i,j)+c(s)*JacInv(i,j)*( &
                  Jac(i+s,j)*F(i+s,j,1)+Jac(i-s,j)*F(i-s,j,1)+ &
                  Jac(i,j+s)*F(i,j+s,1)+Jac(i,j-s)*F(i,j-s,1))

             DF(i,j,2) = DF(i,j,2)+( &
                  0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,1)-0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,1)- &
                  0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,1)+0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,1) &
                  )*b(s)*JacInv(i,j)+c(s)*JacInv(i,j)*( &
                  Jac(i+s,j)*F(i+s,j,2)+Jac(i-s,j)*F(i-s,j,2)+ &
                  Jac(i,j+s)*F(i,j+s,2)+Jac(i,j-s)*F(i,j-s,2))

             DF(i,j,3) = DF(i,j,3)+( &
                 -0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,1)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,1)+ &
                  0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,1)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,1) &
                  )*b(s)*JacInv(i,j)+c(s)*JacInv(i,j)*( &
                  Jac(i+s,j)*F(i+s,j,3)+Jac(i-s,j)*F(i-s,j,3)+ &
                  Jac(i,j+s)*F(i,j+s,3)+Jac(i,j-s)*F(i,j-s,3))

          end do
       end do
    end do

  end subroutine rates_dissipation_interior_mode3


  subroutine rates_dissipation_boundary_mode3(lbndFx,lbndFy,lbndDFx,lbndDFy, &
       F,DF,rho,G,Cdiss,xq,xr,yq,yr,Jac,JacInv,mx,px,my,py,bnd)

    use fd_coeff, only : DL,DR,DI,AL,AR,AI,pI,nbnd

    implicit none

    integer,intent(in) :: lbndFx,lbndFy,lbndDFx,lbndDFy,mx,px,my,py
    real,dimension(lbndFx :,lbndFy :,1:),intent(in) :: F
    real,dimension(lbndDFx:,lbndDFy:,1:),intent(inout) :: DF
    real,intent(in) :: rho,G,Cdiss
    real,dimension(lbndFx:,lbndFy:),intent(in) :: xq,xr,yq,yr,Jac,JacInv
    character(*),intent(in) :: bnd

    integer :: i,j,s

    select case(bnd)

    case('L')

       do j = my,py
          do i = mx,px

             do s = 0,nbnd-1 ! one-sided for q
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,2)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,3) &
                     )*DL(s,i-mx)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,1) &
                     )*DL(s,i-mx)*G*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     -0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,1) &
                     )*DL(s,i-mx)*G*JacInv(i,j)
                
                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(mx+s,j)*F(mx+s,j,:) &
                     )*AL(s,i-mx)*JacInv(i,j)*Cdiss
                
             end do

             DF(i,j,:) = DF(i,j,:)+F(i,j,:)*AI(0)*Cdiss

             do s = 1,pI ! central for r
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,2)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,3)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,2)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,3) &
                     )*DI(s)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,1)+0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,1) &
                     )*DI(s)*G*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,1)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,1) &
                     )*DI(s)*G*JacInv(i,j)

                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(i,j+s)*F(i,j+s,:)+Jac(i,j-s)*F(i,j-s,:) &
                     )*AI(s)*JacInv(i,j)*Cdiss

             end do

          end do
       end do

    case('R')
       
       do j = my,py
          do i = mx,px

             do s = 0,nbnd-1 ! one-sided for q
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,2)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,3) &
                     )*DR(-s,i-px)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,1) &
                     )*DR(-s,i-px)*G*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     -0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,1) &
                     )*DR(-s,i-px)*G*JacInv(i,j)

                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(px-s,j)*F(px-s,j,:) &
                     )*AR(-s,i-px)*JacInv(i,j)*Cdiss

             end do

             DF(i,j,:) = DF(i,j,:)+F(i,j,:)*AI(0)*Cdiss

             do s = 1,pI ! central for r
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,2)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,3)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,2)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,3) &
                     )*DI(s)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,1)+0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,1) &
                     )*DI(s)*G*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,1)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,1) &
                     )*DI(s)*G*JacInv(i,j)
                
                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(i,j+s)*F(i,j+s,:)+Jac(i,j-s)*F(i,j-s,:) &
                     )*AI(s)*JacInv(i,j)*Cdiss

             end do

          end do
       end do

    case('B')
       
       do j = my,py
          do i = mx,px

             do s = 0,nbnd-1 ! one-sided for r
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,2)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,3) &
                     )*DL(s,j-my)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,1) &
                     )*DL(s,j-my)*G*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,1) &
                     )*DL(s,j-my)*G*JacInv(i,j)
                
                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(i,my+s)*F(i,my+s,:) &
                     )*AL(s,j-my)*JacInv(i,j)*Cdiss
                
             end do

             DF(i,j,:) = DF(i,j,:)+F(i,j,:)*AI(0)*Cdiss

             do s = 1,pI ! central for q
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,2)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,3)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,2)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,3) &
                     )*DI(s)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,1)-0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,1) &
                     )*DI(s)*G*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     -0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,1)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,1) &
                     )*DI(s)*G*JacInv(i,j)
                
                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(i+s,j)*F(i+s,j,:)+Jac(i-s,j)*F(i-s,j,:) &
                     )*AI(s)*JacInv(i,j)*Cdiss
                
             end do

          end do
       end do
       
    case('T')
    
       do j = my,py
          do i = mx,px

             do s = 0,nbnd-1 ! one-sided for r
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,2)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,3) &
                     )*DR(-s,j-py)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,1) &
                     )*DR(-s,j-py)*G*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,1) &
                     )*DR(-s,j-py)*G*JacInv(i,j)
                
                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(i,py-s)*F(i,py-s,:) &
                     )*AR(-s,j-py)*JacInv(i,j)*Cdiss
                
             end do

             DF(i,j,:) = DF(i,j,:)+F(i,j,:)*AI(0)*Cdiss

             do s = 1,pI ! central for q
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,2)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,3)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,2)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,3) &
                     )*DI(s)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,1)-0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,1) &
                     )*DI(s)*G*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     -0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,1)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,1) &
                     )*DI(s)*G*JacInv(i,j)
                
                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(i+s,j)*F(i+s,j,:)+Jac(i-s,j)*F(i-s,j,:) &
                     )*AI(s)*JacInv(i,j)*Cdiss
                
             end do

          end do
       end do
       
    case('LB')

       do j = my,py
          do i = mx,px
             do s = 0,nbnd-1 ! one-sided for q and r
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,2)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,3) &
                     )*DL(s,i-mx)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,1) &
                     )*DL(s,i-mx)*G*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     -0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,1) &
                     )*DL(s,i-mx)*G*JacInv(i,j)
                
                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(mx+s,j)*F(mx+s,j,:) &
                     )*AL(s,i-mx)*JacInv(i,j)*Cdiss
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,2)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,3) &
                     )*DL(s,j-my)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,1) &
                     )*DL(s,j-my)*G*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,1) &
                     )*DL(s,j-my)*G*JacInv(i,j)
                
                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(i,my+s)*F(i,my+s,:) &
                     )*AL(s,j-my)*JacInv(i,j)*Cdiss
                
             end do
          end do
       end do

    case('RB')

       do j = my,py
          do i = mx,px
             do s = 0,nbnd-1 ! one-sided for q and r

                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,2)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,3) &
                     )*DR(-s,i-px)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,1) &
                     )*DR(-s,i-px)*G*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     -0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,1) &
                     )*DR(-s,i-px)*G*JacInv(i,j)
                
                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(px-s,j)*F(px-s,j,:) &
                     )*AR(-s,i-px)*JacInv(i,j)*Cdiss
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,2)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,3) &
                     )*DL(s,j-my)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,1) &
                     )*DL(s,j-my)*G*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,1) &
                     )*DL(s,j-my)*G*JacInv(i,j)

                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(i,my+s)*F(i,my+s,:) &
                     )*AL(s,j-my)*JacInv(i,j)*Cdiss

             end do
          end do
       end do

    case('LT')

       do j = my,py
          do i = mx,px
             do s = 0,nbnd-1 ! one-sided for q and r
              
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,2)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,3) &
                     )*DL(s,i-mx)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,1) &
                     )*DL(s,i-mx)*G*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     -0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,1) &
                     )*DL(s,i-mx)*G*JacInv(i,j)

                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(mx+s,j)*F(mx+s,j,:) &
                     )*AL(s,i-mx)*JacInv(i,j)*Cdiss

                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,2)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,3) &
                     )*DR(-s,j-py)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,1) &
                     )*DR(-s,j-py)*G*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,1) &
                     )*DR(-s,j-py)*G*JacInv(i,j)
                
                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(i,py-s)*F(i,py-s,:) &
                     )*AR(-s,j-py)*JacInv(i,j)*Cdiss
                
             end do
          end do
       end do

    case('RT')
    
       do j = my,py
          do i = mx,px
             do s = 0,nbnd-1 ! one-sided for q and r
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,2)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,3) &
                     )*DR(-s,i-px)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,1) &
                     )*DR(-s,i-px)*G*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     -0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,1) &
                     )*DR(-s,i-px)*G*JacInv(i,j)
                
                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(px-s,j)*F(px-s,j,:) &
                     )*AR(-s,i-px)*JacInv(i,j)*Cdiss
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,2)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,3) &
                     )*DR(-s,j-py)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,1) &
                     )*DR(-s,j-py)*G*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,1) &
                     )*DR(-s,j-py)*G*JacInv(i,j)
                
                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(i,py-s)*F(i,py-s,:) &
                     )*AR(-s,j-py)*JacInv(i,j)*Cdiss
                
             end do
          end do
       end do

    end select

  end subroutine rates_dissipation_boundary_mode3


  subroutine rates_dissipation_interior_mode2(lbndFx,lbndFy,lbndDFx,lbndDFy, &
       F,DF,rho,G,M,gamma,Cdiss,xq,xr,yq,yr,Jac,JacInv,mx,px,my,py)

    use fd_coeff, only : DI,AI

    implicit none

    integer,intent(in) :: lbndFx,lbndFy,lbndDFx,lbndDFy,mx,px,my,py
    real,dimension(lbndFx :,lbndFy :,1:),intent(in) :: F
    real,dimension(lbndDFx:,lbndDFy:,1:),intent(inout) :: DF
    real,intent(in) :: rho,G,M,gamma,Cdiss
    real,dimension(lbndFx:,lbndFy:),intent(in) :: xq,xr,yq,yr,Jac,JacInv

    integer :: i,j,s
    real :: a(1:3),b(1:3),c(1:3),d(0:3)

    a=0d0; b=0d0; c=0d0; d=0d0
    do i = 1,ubound(DI,1)
       a(i) = DI(i)/rho
       b(i) = DI(i)*G
       c(i) = DI(i)*M
       d(i) = AI(i)*Cdiss
    end do
    d(0) = 2d0*AI(0)*Cdiss ! coefficient used twice (both for q and r derivatives)

    do j = my,py
       do i = mx,px

          DF(i,j,1) = DF(i,j,1)+d(0)*F(i,j,1)
          DF(i,j,2) = DF(i,j,2)+d(0)*F(i,j,2)
          DF(i,j,3) = DF(i,j,3)+d(0)*F(i,j,3)
          DF(i,j,4) = DF(i,j,4)+d(0)*F(i,j,4)
          DF(i,j,5) = DF(i,j,5)+d(0)*F(i,j,5)

          do s = 1,3

             DF(i,j,1) = DF(i,j,1)+(a(s)*( &
                  0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,3)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,4)- &
                  0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,3)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,4)- &
                  0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,3)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,4)+ &
                  0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,3)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,4) &
                  )+d(s)*( &
                  Jac(i+s,j)*F(i+s,j,1)+Jac(i-s,j)*F(i-s,j,1)+ &
                  Jac(i,j+s)*F(i,j+s,1)+Jac(i,j-s)*F(i,j-s,1)))*JacInv(i,j)

             DF(i,j,2) = DF(i,j,2)+(a(s)*( &
                  0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,4)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,5)- &
                  0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,4)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,5)- &
                  0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,4)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,5)+ &
                  0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,4)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,5) &
                  )+d(s)*( &
                  Jac(i+s,j)*F(i+s,j,2)+Jac(i-s,j)*F(i-s,j,2)+ &
                  Jac(i,j+s)*F(i,j+s,2)+Jac(i,j-s)*F(i,j-s,2)))*JacInv(i,j)

             DF(i,j,3) = DF(i,j,3)+(c(s)*( &
                  0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,1)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,2)*gamma- &
                  0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,1)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,2)*gamma- &
                  0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,1)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,2)*gamma+ &
                  0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,1)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,2)*gamma &
                  )+d(s)*( &
                  Jac(i+s,j)*F(i+s,j,3)+Jac(i-s,j)*F(i-s,j,3)+ &
                  Jac(i,j+s)*F(i,j+s,3)+Jac(i,j-s)*F(i,j-s,3)))*JacInv(i,j)

             DF(i,j,4) = DF(i,j,4)+(b(s)*( &
                  0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,2)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,1)- &
                  0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,2)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,1)- &
                  0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,2)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,1)+ &
                  0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,2)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,1) &
                  )+d(s)*( &
                  Jac(i+s,j)*F(i+s,j,4)+Jac(i-s,j)*F(i-s,j,4)+ &
                  Jac(i,j+s)*F(i,j+s,4)+Jac(i,j-s)*F(i,j-s,4)))*JacInv(i,j)

             DF(i,j,5) = DF(i,j,5)+(c(s)*( &
                  0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,1)*gamma-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,2)- &
                  0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,1)*gamma+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,2)- &
                  0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,1)*gamma+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,2)+ &
                  0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,1)*gamma-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,2) &
                  )+d(s)*( &
                  Jac(i+s,j)*F(i+s,j,5)+Jac(i-s,j)*F(i-s,j,5)+ &
                  Jac(i,j+s)*F(i,j+s,5)+Jac(i,j-s)*F(i,j-s,5)))*JacInv(i,j)
             
          end do
       end do
    end do

  end subroutine rates_dissipation_interior_mode2


  subroutine rates_dissipation_boundary_mode2(lbndFx,lbndFy,lbndDFx,lbndDFy, &
       F,DF,rho,G,M,gamma,Cdiss,xq,xr,yq,yr,Jac,JacInv,mx,px,my,py,bnd)

    use fd_coeff, only : DL,DR,DI,AL,AR,AI,pI,nbnd

    implicit none

    integer,intent(in) :: lbndFx,lbndFy,lbndDFx,lbndDFy,mx,px,my,py
    real,dimension(lbndFx :,lbndFy :,1:),intent(in) :: F
    real,dimension(lbndDFx:,lbndDFy:,1:),intent(inout) :: DF
    real,intent(in) :: rho,G,M,gamma,Cdiss
    real,dimension(lbndFx:,lbndFy:),intent(in) :: xq,xr,yq,yr,Jac,JacInv
    character(*),intent(in) :: bnd

    integer :: i,j,s

    select case(bnd)

    case('L')

       do j = my,py
          do i = mx,px

             do s = 0,nbnd-1 ! one-sided for q
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,3)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,4) &
                     )*DL(s,i-mx)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,4)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,5) &
                     )*DL(s,i-mx)/rho*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,1)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,2)*gamma &
                     )*DL(s,i-mx)*M*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,2)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,1) &
                     )*DL(s,i-mx)*G*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,1)*gamma-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,2) &
                     )*DL(s,i-mx)*M*JacInv(i,j)

                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(mx+s,j)*F(mx+s,j,:) &
                     )*AL(s,i-mx)*JacInv(i,j)*Cdiss
                
             end do

             DF(i,j,:) = DF(i,j,:)+F(i,j,:)*AI(0)*Cdiss

             do s = 1,pI ! central for r
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,3)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,4)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,3)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,4) &
                     )*DI(s)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,4)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,5)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,4)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,5) &
                     )*DI(s)/rho*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,1)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,2)*gamma+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,1)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,2)*gamma &
                     )*DI(s)*M*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,2)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,1)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,2)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,1) &
                     )*DI(s)*G*JacInv(i,j)

                DF(i,j,5) = DF(i,j,5)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,1)*gamma+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,2)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,1)*gamma-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,2) &
                     )*DI(s)*M*JacInv(i,j)

                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(i,j+s)*F(i,j+s,:)+Jac(i,j-s)*F(i,j-s,:) &
                     )*AI(s)*JacInv(i,j)*Cdiss

             end do

          end do
       end do

    case('R')
       
       do j = my,py
          do i = mx,px

             do s = 0,nbnd-1 ! one-sided for q
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,3)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,4) &
                     )*DR(-s,i-px)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,4)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,5) &
                     )*DR(-s,i-px)/rho*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,1)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,2)*gamma &
                     )*DR(-s,i-px)*M*JacInv(i,j)

                DF(i,j,4) = DF(i,j,4)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,2)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,1) &
                     )*DR(-s,i-px)*G*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,1)*gamma-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,2) &
                     )*DR(-s,i-px)*M*JacInv(i,j)

                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(px-s,j)*F(px-s,j,:) &
                     )*AR(-s,i-px)*JacInv(i,j)*Cdiss

             end do

             DF(i,j,:) = DF(i,j,:)+F(i,j,:)*AI(0)*Cdiss

             do s = 1,pI ! central for r
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,3)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,4)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,3)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,4) &
                     )*DI(s)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,4)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,5)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,4)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,5) &
                     )*DI(s)/rho*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,1)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,2)*gamma+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,1)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,2)*gamma &
                     )*DI(s)*M*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,2)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,1)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,2)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,1) &
                     )*DI(s)*G*JacInv(i,j)

                DF(i,j,5) = DF(i,j,5)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,1)*gamma+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,2)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,1)*gamma-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,2) &
                     )*DI(s)*M*JacInv(i,j)

                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(i,j+s)*F(i,j+s,:)+Jac(i,j-s)*F(i,j-s,:) &
                     )*AI(s)*JacInv(i,j)*Cdiss

             end do

          end do
       end do

    case('B')
       
       do j = my,py
          do i = mx,px

             do s = 0,nbnd-1 ! one-sided for r
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,3)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,4) &
                     )*DL(s,j-my)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,4)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,5) &
                     )*DL(s,j-my)/rho*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,1)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,2)*gamma &
                     )*DL(s,j-my)*M*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,2)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,1) &
                     )*DL(s,j-my)*G*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,1)*gamma+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,2) &
                     )*DL(s,j-my)*M*JacInv(i,j)
                
                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(i,my+s)*F(i,my+s,:) &
                     )*AL(s,j-my)*JacInv(i,j)*Cdiss
                
             end do

             DF(i,j,:) = DF(i,j,:)+F(i,j,:)*AI(0)*Cdiss

             do s = 1,pI ! central for q
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,3)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,4)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,3)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,4) &
                     )*DI(s)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,4)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,5)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,4)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,5) &
                     )*DI(s)/rho*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,1)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,2)*gamma- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,1)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,2)*gamma &
                     )*DI(s)*M*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,2)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,1)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,2)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,1) &
                     )*DI(s)*G*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,1)*gamma-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,2)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,1)*gamma+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,2) &
                     )*DI(s)*M*JacInv(i,j)
                
                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(i+s,j)*F(i+s,j,:)+Jac(i-s,j)*F(i-s,j,:) &
                     )*AI(s)*JacInv(i,j)*Cdiss
                
             end do

          end do
       end do
       
    case('T')
    
       do j = my,py
          do i = mx,px

             do s = 0,nbnd-1 ! one-sided for r
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,3)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,4) &
                     )*DR(-s,j-py)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,4)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,5) &
                     )*DR(-s,j-py)/rho*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,1)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,2)*gamma &
                     )*DR(-s,j-py)*M*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,2)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,1) &
                     )*DR(-s,j-py)*G*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,1)*gamma+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,2) &
                     )*DR(-s,j-py)*M*JacInv(i,j)
                
                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(i,py-s)*F(i,py-s,:) &
                     )*AR(-s,j-py)*JacInv(i,j)*Cdiss
                
             end do

             DF(i,j,:) = DF(i,j,:)+F(i,j,:)*AI(0)*Cdiss

             do s = 1,pI ! central for q
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,3)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,4)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,3)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,4) &
                     )*DI(s)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,4)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,5)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,4)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,5) &
                     )*DI(s)/rho*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,1)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,2)*gamma- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,1)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,2)*gamma &
                     )*DI(s)*M*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,2)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,1)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,2)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,1) &
                     )*DI(s)*G*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,1)*gamma-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,2)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,1)*gamma+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,2) &
                     )*DI(s)*M*JacInv(i,j)
                
                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(i+s,j)*F(i+s,j,:)+Jac(i-s,j)*F(i-s,j,:) &
                     )*AI(s)*JacInv(i,j)*Cdiss
                
             end do

          end do
       end do
       
    case('LB')

       do j = my,py
          do i = mx,px
             do s = 0,nbnd-1 ! one-sided for q and r
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,3)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,4) &
                     )*DL(s,i-mx)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,4)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,5) &
                     )*DL(s,i-mx)/rho*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,1)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,2)*gamma &
                     )*DL(s,i-mx)*M*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,2)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,1) &
                     )*DL(s,i-mx)*G*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,1)*gamma-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,2) &
                     )*DL(s,i-mx)*M*JacInv(i,j)
                
                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(mx+s,j)*F(mx+s,j,:) &
                     )*AL(s,i-mx)*JacInv(i,j)*Cdiss
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,3)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,4) &
                     )*DL(s,j-my)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,4)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,5) &
                     )*DL(s,j-my)/rho*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,1)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,2)*gamma &
                     )*DL(s,j-my)*M*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,2)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,1) &
                     )*DL(s,j-my)*G*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,1)*gamma+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,2) &
                     )*DL(s,j-my)*M*JacInv(i,j)
                
                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(i,my+s)*F(i,my+s,:) &
                     )*AL(s,j-my)*JacInv(i,j)*Cdiss
                
             end do
          end do
       end do

    case('RB')

       do j = my,py
          do i = mx,px
             do s = 0,nbnd-1 ! one-sided for q and r

                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,3)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,4) &
                     )*DR(-s,i-px)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,4)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,5) &
                     )*DR(-s,i-px)/rho*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,1)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,2)*gamma &
                     )*DR(-s,i-px)*M*JacInv(i,j)

                DF(i,j,4) = DF(i,j,4)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,2)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,1) &
                     )*DR(-s,i-px)*G*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,1)*gamma-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,2) &
                     )*DR(-s,i-px)*M*JacInv(i,j)

                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(px-s,j)*F(px-s,j,:) &
                     )*AR(-s,i-px)*JacInv(i,j)*Cdiss
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,3)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,4) &
                     )*DL(s,j-my)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,4)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,5) &
                     )*DL(s,j-my)/rho*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,1)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,2)*gamma &
                     )*DL(s,j-my)*M*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,2)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,1) &
                     )*DL(s,j-my)*G*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,1)*gamma+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,2) &
                     )*DL(s,j-my)*M*JacInv(i,j)
                
                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(i,my+s)*F(i,my+s,:) &
                     )*AL(s,j-my)*JacInv(i,j)*Cdiss

             end do
          end do
       end do

    case('LT')

       do j = my,py
          do i = mx,px
             do s = 0,nbnd-1 ! one-sided for q and r
              
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,3)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,4) &
                     )*DL(s,i-mx)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,4)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,5) &
                     )*DL(s,i-mx)/rho*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,1)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,2)*gamma &
                     )*DL(s,i-mx)*M*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,2)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,1) &
                     )*DL(s,i-mx)*G*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,1)*gamma-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,2) &
                     )*DL(s,i-mx)*M*JacInv(i,j)
                
                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(mx+s,j)*F(mx+s,j,:) &
                     )*AL(s,i-mx)*JacInv(i,j)*Cdiss

                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,3)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,4) &
                     )*DR(-s,j-py)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,4)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,5) &
                     )*DR(-s,j-py)/rho*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,1)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,2)*gamma &
                     )*DR(-s,j-py)*M*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,2)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,1) &
                     )*DR(-s,j-py)*G*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,1)*gamma+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,2) &
                     )*DR(-s,j-py)*M*JacInv(i,j)
                
                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(i,py-s)*F(i,py-s,:) &
                     )*AR(-s,j-py)*JacInv(i,j)*Cdiss
                
             end do
          end do
       end do

    case('RT')
    
       do j = my,py
          do i = mx,px
             do s = 0,nbnd-1 ! one-sided for q and r
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,3)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,4) &
                     )*DR(-s,i-px)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,4)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,5) &
                     )*DR(-s,i-px)/rho*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,1)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,2)*gamma &
                     )*DR(-s,i-px)*M*JacInv(i,j)

                DF(i,j,4) = DF(i,j,4)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,2)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,1) &
                     )*DR(-s,i-px)*G*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,1)*gamma-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,2) &
                     )*DR(-s,i-px)*M*JacInv(i,j)

                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(px-s,j)*F(px-s,j,:) &
                     )*AR(-s,i-px)*JacInv(i,j)*Cdiss
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,3)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,4) &
                     )*DR(-s,j-py)/rho*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,4)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,5) &
                     )*DR(-s,j-py)/rho*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,1)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,2)*gamma &
                     )*DR(-s,j-py)*M*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,2)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,1) &
                     )*DR(-s,j-py)*G*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,1)*gamma+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,2) &
                     )*DR(-s,j-py)*M*JacInv(i,j)
                
                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(i,py-s)*F(i,py-s,:) &
                     )*AR(-s,j-py)*JacInv(i,j)*Cdiss
                
             end do
          end do
       end do

    end select

  end subroutine rates_dissipation_boundary_mode2


  subroutine rates_interior_mode2_pml(lbndFx,lbndFy,lbndDFx,lbndDFy,lbndWx,lbndWy, &
      F,DF,Wx,DWx,Wy,DWy,&
      rho,G,M,gamma,&
      pml,&
      x,y,&
      xq,xr,yq,yr,JacInv,mx,px,my,py)

    use fd_coeff, only : DI
    use material, only : pml_type

    implicit none

    integer,intent(in) :: lbndFx,lbndFy,lbndDFx,lbndDFy,lbndWx,lbndWy,mx,px,my,py
    real,dimension(lbndFx :,lbndFy :,1:),intent(in) :: F
    real,dimension(lbndDFx:,lbndDFy:,1:),intent(inout) :: DF
    real,dimension(lbndWx:,lbndWy:,1:),intent(inout) :: DWx,DWy
    real,dimension(lbndWx:,lbndWy:,1:),intent(in) :: Wx,Wy
    real,intent(in) :: rho,G,M,gamma
    type(pml_type),intent(in) :: pml
    real,dimension(lbndFx:,lbndFy:),intent(in) :: xq,xr,yq,yr,JacInv,x,y

    integer :: i,j,s
    real :: hxi,hyi,pmlsx,pmlsy
    real,dimension(3) :: a,b,c,DIP

    pmlsx = pml%sigma
    pmlsy = pml%sigma

    a=0d0; b=0d0; c=0d0; DIP=0d0
    do i = 1,ubound(DI,1)
      a(i) = DI(i)/rho
      b(i) = DI(i)*G
      c(i) = DI(i)*M
      DIP(i) = DI(i)
    end do

    hxi = JacInv(mx,my)*yr(mx,my)
    hyi = JacInv(mx,my)*xq(mx,my)
    do j = my,py
      do i = mx,px

        if(pml%n > 0d0) then
          pmlsx = pml%sigma*abs(max(0d0,max(pml%mx-x(i,j),x(i,j)-pml%px)))**pml%n
          pmlsy = pml%sigma*abs(max(0d0,max(pml%my-y(i,j),y(i,j)-pml%py)))**pml%n
        end if

        if(pml%pmlx) then
          DF(i,j,1) = DF(i,j,1) - Wx(i,j,3)/rho
          DF(i,j,2) = DF(i,j,2) - Wx(i,j,4)/rho
          DF(i,j,3) = DF(i,j,3) - Wx(i,j,1)*M
          DF(i,j,4) = DF(i,j,4) - Wx(i,j,2)*G
          DF(i,j,5) = DF(i,j,5) - Wx(i,j,1)*gamma*M

          DWx(i,j,:) = DWx(i,j,:) - (pmlsx+pml%kappa)*Wx(i,j,:)
        end if

        if(pml%pmly) then
          DF(i,j,1) = DF(i,j,1) - Wy(i,j,4)/rho
          DF(i,j,2) = DF(i,j,2) - Wy(i,j,5)/rho
          DF(i,j,3) = DF(i,j,3) - Wy(i,j,2)*gamma*M
          DF(i,j,4) = DF(i,j,4) - Wy(i,j,1)*G
          DF(i,j,5) = DF(i,j,5) - Wy(i,j,2)*M

          DWy(i,j,:) = DWy(i,j,:) - (pmlsy+pml%kappa)*Wy(i,j,:)
        end if

        do s = 1,3

          if(pml%pmlx) then
            DWx(i,j,:) = DWx(i,j,:) + pmlsx*hxi*(F(i+s,j,:)-F(i-s,j,:))*DIP(s)
          end if

          if(pml%pmly) then
            DWy(i,j,:) = DWy(i,j,:) + pmlsy*hyi*(F(i,j+s,:)-F(i,j-s,:))*DIP(s)
          end if

          DF(i,j,1) = DF(i,j,1)+( &
            hxi*(F(i+s,j,3)-F(i-s,j,3))+&
            hyi*(F(i,j+s,4)-F(i,j-s,4)) &
            )*a(s)

          DF(i,j,2) = DF(i,j,2)+( &
            hxi*(F(i+s,j,4)-F(i-s,j,4))+&
            hyi*(F(i,j+s,5)-F(i,j-s,5)) &
            )*a(s)

          DF(i,j,3) = DF(i,j,3)+( &
            hxi*(F(i+s,j,1)-F(i-s,j,1))+&
            hyi*(F(i,j+s,2)-F(i,j-s,2))*gamma &
            )*c(s)

          DF(i,j,4) = DF(i,j,4)+( &
            hxi*(F(i+s,j,2)-F(i-s,j,2))+&
            hyi*(F(i,j+s,1)-F(i,j-s,1))&
            )*b(s)

          DF(i,j,5) = DF(i,j,5)+( &
            hxi*(F(i+s,j,1)-F(i-s,j,1))*gamma+&
            hyi*(F(i,j+s,2)-F(i,j-s,2))&
            )*c(s)

        end do
      end do
    end do

  end subroutine rates_interior_mode2_pml

  subroutine rates_boundary_mode2_pml(lbndFx,lbndFy,lbndDFx,lbndDFy,lbndWx,lbndWy, &
      F,DF,Wx,DWx,Wy,DWy,&
      rho,G,M,gamma,&
      pml,&
      x,y,&
      xq,xr,yq,yr,JacInv,mx,px,my,py,bnd)

    use fd_coeff, only : DL,DR,DI,pI,nbnd
    use material, only : pml_type

    implicit none

    integer,intent(in) :: lbndFx,lbndFy,lbndDFx,lbndDFy,lbndWx,lbndWy,mx,px,my,py
    real,dimension(lbndFx :,lbndFy :,1:),intent(in) :: F
    real,dimension(lbndDFx:,lbndDFy:,1:),intent(inout) :: DF
    real,dimension(lbndWx:,lbndWy:,1:),intent(inout) :: DWx,DWy
    real,dimension(lbndWx:,lbndWy:,1:),intent(in) :: Wx,Wy
    real,intent(in) :: rho,G,M,gamma
    real,dimension(lbndFx:,lbndFy:),intent(in) :: xq,xr,yq,yr,JacInv,x,y
    character(*),intent(in) :: bnd
    type(pml_type), intent(in) :: pml

    real,dimension(5) :: dx,dy
    integer :: i,j,s
    real :: hxi,hyi,pmlsx,pmlsy

    pmlsx = pml%sigma
    pmlsy = pml%sigma
    hxi = JacInv(mx,my)*yr(mx,my)
    hyi = JacInv(mx,my)*xq(mx,my)

    select case(bnd)

    case('L')

      do j = my,py
        do i = mx,px
          if(pml%n > 0d0) then
            pmlsx = pml%sigma*abs(max(0d0,max(pml%mx-x(i,j),x(i,j)-pml%px)))**pml%n
            pmlsy = pml%sigma*abs(max(0d0,max(pml%my-y(i,j),y(i,j)-pml%py)))**pml%n
          end if

          if(pml%pmlx) then
            DF(i,j,1) = DF(i,j,1) - Wx(i,j,3)/rho
            DF(i,j,2) = DF(i,j,2) - Wx(i,j,4)/rho
            DF(i,j,3) = DF(i,j,3) - Wx(i,j,1)*M
            DF(i,j,4) = DF(i,j,4) - Wx(i,j,2)*G
            DF(i,j,5) = DF(i,j,5) - Wx(i,j,1)*gamma*M

            DWx(i,j,:) = DWx(i,j,:) - (pmlsx+pml%kappa)*Wx(i,j,:)
          end if

          if(pml%pmly) then
            DF(i,j,1) = DF(i,j,1) - Wy(i,j,4)/rho
            DF(i,j,2) = DF(i,j,2) - Wy(i,j,5)/rho
            DF(i,j,3) = DF(i,j,3) - Wy(i,j,2)*gamma*M
            DF(i,j,4) = DF(i,j,4) - Wy(i,j,1)*G
            DF(i,j,5) = DF(i,j,5) - Wy(i,j,2)*M

            DWy(i,j,:) = DWy(i,j,:) - (pmlsy+pml%kappa)*Wy(i,j,:)
          end if

          do s = 0,nbnd-1 ! one-sided for x

            if(pml%pmlx) then
              DWx(i,j,:) = DWx(i,j,:) + pmlsx*hxi*F(mx+s,j,:)*DL(s,i-mx)
            end if

            DF(i,j,1) = DF(i,j,1)+hxi*F(mx+s,j,3)*DL(s,i-mx)/rho
            DF(i,j,2) = DF(i,j,2)+hxi*F(mx+s,j,4)*DL(s,i-mx)/rho
            DF(i,j,3) = DF(i,j,3)+hxi*F(mx+s,j,1)*DL(s,i-mx)*M
            DF(i,j,4) = DF(i,j,4)+hxi*F(mx+s,j,2)*DL(s,i-mx)*G
            DF(i,j,5) = DF(i,j,5)+hxi*F(mx+s,j,1)*DL(s,i-mx)*M*gamma

          end do

          do s = 1,pI ! central for y
            if(pml%pmly) then
              DWy(i,j,:) = DWy(i,j,:) + pmlsy*hyi*(F(i,j+s,:)-F(i,j-s,:))*DI(s)
            end if

            DF(i,j,1) = DF(i,j,1)+hyi*(F(i,j+s,4)-F(i,j-s,4))*DI(s)/rho
            DF(i,j,2) = DF(i,j,2)+hyi*(F(i,j+s,5)-F(i,j-s,5))*DI(s)/rho
            DF(i,j,3) = DF(i,j,3)+hyi*(F(i,j+s,2)-F(i,j-s,2))*DI(s)*M*gamma
            DF(i,j,4) = DF(i,j,4)+hyi*(F(i,j+s,1)-F(i,j-s,1))*DI(s)*G
            DF(i,j,5) = DF(i,j,5)+hyi*(F(i,j+s,2)-F(i,j-s,2))*DI(s)*M

          end do

        end do
      end do

    case('R')

      do j = my,py
        do i = mx,px
          if(pml%n > 0d0) then
            pmlsx = pml%sigma*abs(max(0d0,max(pml%mx-x(i,j),x(i,j)-pml%px)))**pml%n
            pmlsy = pml%sigma*abs(max(0d0,max(pml%my-y(i,j),y(i,j)-pml%py)))**pml%n
          end if

          if(pml%pmlx) then
            DF(i,j,1) = DF(i,j,1) - Wx(i,j,3)/rho
            DF(i,j,2) = DF(i,j,2) - Wx(i,j,4)/rho
            DF(i,j,3) = DF(i,j,3) - Wx(i,j,1)*M
            DF(i,j,4) = DF(i,j,4) - Wx(i,j,2)*G
            DF(i,j,5) = DF(i,j,5) - Wx(i,j,1)*gamma*M

            DWx(i,j,:) = DWx(i,j,:) - (pmlsx+pml%kappa)*Wx(i,j,:)
          end if

          if(pml%pmly) then
            DF(i,j,1) = DF(i,j,1) - Wy(i,j,4)/rho
            DF(i,j,2) = DF(i,j,2) - Wy(i,j,5)/rho
            DF(i,j,3) = DF(i,j,3) - Wy(i,j,2)*gamma*M
            DF(i,j,4) = DF(i,j,4) - Wy(i,j,1)*G
            DF(i,j,5) = DF(i,j,5) - Wy(i,j,2)*M

            DWy(i,j,:) = DWy(i,j,:) - (pmlsy+pml%kappa)*Wy(i,j,:)
          end if

          do s = 0,nbnd-1 ! one-sided for x
            if(pml%pmlx) then
              DWx(i,j,:) = DWx(i,j,:) + pmlsx*hxi*F(px-s,j,:)*DR(-s,i-px)
            end if

            DF(i,j,1) = DF(i,j,1)+hxi*F(px-s,j,3)*DR(-s,i-px)/rho
            DF(i,j,2) = DF(i,j,2)+hxi*F(px-s,j,4)*DR(-s,i-px)/rho
            DF(i,j,3) = DF(i,j,3)+hxi*F(px-s,j,1)*DR(-s,i-px)*M
            DF(i,j,4) = DF(i,j,4)+hxi*F(px-s,j,2)*DR(-s,i-px)*G
            DF(i,j,5) = DF(i,j,5)+hxi*F(px-s,j,1)*DR(-s,i-px)*M*gamma

          end do

          do s = 1,pI ! central for y
            if(pml%pmly) then
              DWy(i,j,:) = DWy(i,j,:) + pmlsy*hyi*(F(i,j+s,:)-F(i,j-s,:))*DI(s)
            end if

            DF(i,j,1) = DF(i,j,1)+hyi*(F(i,j+s,4)-F(i,j-s,4))*DI(s)/rho
            DF(i,j,2) = DF(i,j,2)+hyi*(F(i,j+s,5)-F(i,j-s,5))*DI(s)/rho
            DF(i,j,3) = DF(i,j,3)+hyi*(F(i,j+s,2)-F(i,j-s,2))*DI(s)*M*gamma
            DF(i,j,4) = DF(i,j,4)+hyi*(F(i,j+s,1)-F(i,j-s,1))*DI(s)*G
            DF(i,j,5) = DF(i,j,5)+hyi*(F(i,j+s,2)-F(i,j-s,2))*DI(s)*M

          end do

        end do
      end do

    case('B')

      do j = my,py
        do i = mx,px
          if(pml%n > 0d0) then
            pmlsx = pml%sigma*abs(max(0d0,max(pml%mx-x(i,j),x(i,j)-pml%px)))**pml%n
            pmlsy = pml%sigma*abs(max(0d0,max(pml%my-y(i,j),y(i,j)-pml%py)))**pml%n
          end if

          if(pml%pmlx) then
            DF(i,j,1) = DF(i,j,1) - Wx(i,j,3)/rho
            DF(i,j,2) = DF(i,j,2) - Wx(i,j,4)/rho
            DF(i,j,3) = DF(i,j,3) - Wx(i,j,1)*M
            DF(i,j,4) = DF(i,j,4) - Wx(i,j,2)*G
            DF(i,j,5) = DF(i,j,5) - Wx(i,j,1)*gamma*M

            DWx(i,j,:) = DWx(i,j,:) - (pmlsx+pml%kappa)*Wx(i,j,:)
          end if

          if(pml%pmly) then
            DF(i,j,1) = DF(i,j,1) - Wy(i,j,4)/rho
            DF(i,j,2) = DF(i,j,2) - Wy(i,j,5)/rho
            DF(i,j,3) = DF(i,j,3) - Wy(i,j,2)*gamma*M
            DF(i,j,4) = DF(i,j,4) - Wy(i,j,1)*G
            DF(i,j,5) = DF(i,j,5) - Wy(i,j,2)*M

            DWy(i,j,:) = DWy(i,j,:) - (pmlsy+pml%kappa)*Wy(i,j,:)
          end if

          do s = 1,pI ! central for x
            if(pml%pmlx) then
              DWx(i,j,:) = DWx(i,j,:) + pmlsx*hxi*(F(i+s,j,:)-F(i-s,j,:))*DI(s)
            end if

            DF(i,j,1) = DF(i,j,1)+hxi*(F(i+s,j,3)-F(i-s,j,3))*DI(s)/rho
            DF(i,j,2) = DF(i,j,2)+hxi*(F(i+s,j,4)-F(i-s,j,4))*DI(s)/rho
            DF(i,j,3) = DF(i,j,3)+hxi*(F(i+s,j,1)-F(i-s,j,1))*DI(s)*M
            DF(i,j,4) = DF(i,j,4)+hxi*(F(i+s,j,2)-F(i-s,j,2))*DI(s)*G
            DF(i,j,5) = DF(i,j,5)+hxi*(F(i+s,j,1)-F(i-s,j,1))*DI(s)*M*gamma

          end do

          do s = 0,nbnd-1 ! one-sided for y
            if(pml%pmly) then
              DWy(i,j,:) = DWy(i,j,:) + pmlsy*hyi*F(i,my+s,:)*DL(s,j-my)
            end if

            DF(i,j,1) = DF(i,j,1)+hyi*F(i,my+s,4)*DL(s,j-my)/rho
            DF(i,j,2) = DF(i,j,2)+hyi*F(i,my+s,5)*DL(s,j-my)/rho
            DF(i,j,3) = DF(i,j,3)+hyi*F(i,my+s,2)*DL(s,j-my)*M*gamma
            DF(i,j,4) = DF(i,j,4)+hyi*F(i,my+s,1)*DL(s,j-my)*G
            DF(i,j,5) = DF(i,j,5)+hyi*F(i,my+s,2)*DL(s,j-my)*M

          end do

        end do
      end do

    case('T')

      do j = my,py
        do i = mx,px
          if(pml%n > 0d0) then
            pmlsx = pml%sigma*abs(max(0d0,max(pml%mx-x(i,j),x(i,j)-pml%px)))**pml%n
            pmlsy = pml%sigma*abs(max(0d0,max(pml%my-y(i,j),y(i,j)-pml%py)))**pml%n
          end if

          if(pml%pmlx) then
            DF(i,j,1) = DF(i,j,1) - Wx(i,j,3)/rho
            DF(i,j,2) = DF(i,j,2) - Wx(i,j,4)/rho
            DF(i,j,3) = DF(i,j,3) - Wx(i,j,1)*M
            DF(i,j,4) = DF(i,j,4) - Wx(i,j,2)*G
            DF(i,j,5) = DF(i,j,5) - Wx(i,j,1)*gamma*M

            DWx(i,j,:) = DWx(i,j,:) - (pmlsx+pml%kappa)*Wx(i,j,:)
          end if

          if(pml%pmly) then
            DF(i,j,1) = DF(i,j,1) - Wy(i,j,4)/rho
            DF(i,j,2) = DF(i,j,2) - Wy(i,j,5)/rho
            DF(i,j,3) = DF(i,j,3) - Wy(i,j,2)*gamma*M
            DF(i,j,4) = DF(i,j,4) - Wy(i,j,1)*G
            DF(i,j,5) = DF(i,j,5) - Wy(i,j,2)*M

            DWy(i,j,:) = DWy(i,j,:) - (pmlsy+pml%kappa)*Wy(i,j,:)
          end if

          do s = 1,pI ! central for x
            if(pml%pmlx) then
              DWx(i,j,:) = DWx(i,j,:) + pmlsx*hxi*(F(i+s,j,:)-F(i-s,j,:))*DI(s)
            end if

            DF(i,j,1) = DF(i,j,1)+hxi*(F(i+s,j,3)-F(i-s,j,3))*DI(s)/rho
            DF(i,j,2) = DF(i,j,2)+hxi*(F(i+s,j,4)-F(i-s,j,4))*DI(s)/rho
            DF(i,j,3) = DF(i,j,3)+hxi*(F(i+s,j,1)-F(i-s,j,1))*DI(s)*M
            DF(i,j,4) = DF(i,j,4)+hxi*(F(i+s,j,2)-F(i-s,j,2))*DI(s)*G
            DF(i,j,5) = DF(i,j,5)+hxi*(F(i+s,j,1)-F(i-s,j,1))*DI(s)*M*gamma

          end do

          do s = 0,nbnd-1 ! one-sided for y
            if(pml%pmly) then
              DWy(i,j,:) = DWy(i,j,:) + pmlsy*hyi*F(i,py-s,:)*DR(-s,j-py)
            end if

            DF(i,j,1) = DF(i,j,1)+hyi*F(i,py-s,4)*DR(-s,j-py)/rho
            DF(i,j,2) = DF(i,j,2)+hyi*F(i,py-s,5)*DR(-s,j-py)/rho
            DF(i,j,3) = DF(i,j,3)+hyi*F(i,py-s,2)*DR(-s,j-py)*M*gamma
            DF(i,j,4) = DF(i,j,4)+hyi*F(i,py-s,1)*DR(-s,j-py)*G
            DF(i,j,5) = DF(i,j,5)+hyi*F(i,py-s,2)*DR(-s,j-py)*M

          end do

        end do
      end do

    case('LB')

      do j = my,py
        do i = mx,px
          if(pml%n > 0d0) then
            pmlsx = pml%sigma*abs(max(0d0,max(pml%mx-x(i,j),x(i,j)-pml%px)))**pml%n
            pmlsy = pml%sigma*abs(max(0d0,max(pml%my-y(i,j),y(i,j)-pml%py)))**pml%n
          end if

          if(pml%pmlx) then
            DF(i,j,1) = DF(i,j,1) - Wx(i,j,3)/rho
            DF(i,j,2) = DF(i,j,2) - Wx(i,j,4)/rho
            DF(i,j,3) = DF(i,j,3) - Wx(i,j,1)*M
            DF(i,j,4) = DF(i,j,4) - Wx(i,j,2)*G
            DF(i,j,5) = DF(i,j,5) - Wx(i,j,1)*gamma*M

            DWx(i,j,:) = DWx(i,j,:) - (pmlsx+pml%kappa)*Wx(i,j,:)
          end if

          if(pml%pmly) then
            DF(i,j,1) = DF(i,j,1) - Wy(i,j,4)/rho
            DF(i,j,2) = DF(i,j,2) - Wy(i,j,5)/rho
            DF(i,j,3) = DF(i,j,3) - Wy(i,j,2)*gamma*M
            DF(i,j,4) = DF(i,j,4) - Wy(i,j,1)*G
            DF(i,j,5) = DF(i,j,5) - Wy(i,j,2)*M

            DWy(i,j,:) = DWy(i,j,:) - (pmlsy+pml%kappa)*Wy(i,j,:)
          end if

          do s = 0,nbnd-1 ! one-sided for q and r

            if(pml%pmlx) then
              DWx(i,j,:) = DWx(i,j,:) + pmlsx*hxi*F(mx+s,j,:)*DL(s,i-mx)
            end if

            DF(i,j,1) = DF(i,j,1)+hxi*F(mx+s,j,3)*DL(s,i-mx)/rho
            DF(i,j,2) = DF(i,j,2)+hxi*F(mx+s,j,4)*DL(s,i-mx)/rho
            DF(i,j,3) = DF(i,j,3)+hxi*F(mx+s,j,1)*DL(s,i-mx)*M
            DF(i,j,4) = DF(i,j,4)+hxi*F(mx+s,j,2)*DL(s,i-mx)*G
            DF(i,j,5) = DF(i,j,5)+hxi*F(mx+s,j,1)*DL(s,i-mx)*M*gamma

            if(pml%pmly) then
              DWy(i,j,:) = DWy(i,j,:) + pmlsy*hyi*F(i,my+s,:)*DL(s,j-my)
            end if

            DF(i,j,1) = DF(i,j,1)+hyi*F(i,my+s,4)*DL(s,j-my)/rho
            DF(i,j,2) = DF(i,j,2)+hyi*F(i,my+s,5)*DL(s,j-my)/rho
            DF(i,j,3) = DF(i,j,3)+hyi*F(i,my+s,2)*DL(s,j-my)*M*gamma
            DF(i,j,4) = DF(i,j,4)+hyi*F(i,my+s,1)*DL(s,j-my)*G
            DF(i,j,5) = DF(i,j,5)+hyi*F(i,my+s,2)*DL(s,j-my)*M

          end do
        end do
      end do

    case('RB')

      do j = my,py
        do i = mx,px
          if(pml%n > 0d0) then
            pmlsx = pml%sigma*abs(max(0d0,max(pml%mx-x(i,j),x(i,j)-pml%px)))**pml%n
            pmlsy = pml%sigma*abs(max(0d0,max(pml%my-y(i,j),y(i,j)-pml%py)))**pml%n
          end if

          if(pml%pmlx) then
            DF(i,j,1) = DF(i,j,1) - Wx(i,j,3)/rho
            DF(i,j,2) = DF(i,j,2) - Wx(i,j,4)/rho
            DF(i,j,3) = DF(i,j,3) - Wx(i,j,1)*M
            DF(i,j,4) = DF(i,j,4) - Wx(i,j,2)*G
            DF(i,j,5) = DF(i,j,5) - Wx(i,j,1)*gamma*M

            DWx(i,j,:) = DWx(i,j,:) - (pmlsx+pml%kappa)*Wx(i,j,:)
          end if

          if(pml%pmly) then
            DF(i,j,1) = DF(i,j,1) - Wy(i,j,4)/rho
            DF(i,j,2) = DF(i,j,2) - Wy(i,j,5)/rho
            DF(i,j,3) = DF(i,j,3) - Wy(i,j,2)*gamma*M
            DF(i,j,4) = DF(i,j,4) - Wy(i,j,1)*G
            DF(i,j,5) = DF(i,j,5) - Wy(i,j,2)*M

            DWy(i,j,:) = DWy(i,j,:) - (pmlsy+pml%kappa)*Wy(i,j,:)
          end if

          do s = 0,nbnd-1 ! one-sided for q and r

            if(pml%pmlx) then
              DWx(i,j,:) = DWx(i,j,:) + pmlsx*hxi*F(px-s,j,:)*DR(-s,i-px)
            end if

            DF(i,j,1) = DF(i,j,1)+hxi*F(px-s,j,3)*DR(-s,i-px)/rho
            DF(i,j,2) = DF(i,j,2)+hxi*F(px-s,j,4)*DR(-s,i-px)/rho
            DF(i,j,3) = DF(i,j,3)+hxi*F(px-s,j,1)*DR(-s,i-px)*M
            DF(i,j,4) = DF(i,j,4)+hxi*F(px-s,j,2)*DR(-s,i-px)*G
            DF(i,j,5) = DF(i,j,5)+hxi*F(px-s,j,1)*DR(-s,i-px)*M*gamma

            if(pml%pmly) then
              DWy(i,j,:) = DWy(i,j,:) + pmlsy*hyi*F(i,my+s,:)*DL(s,j-my)
            end if

            DF(i,j,1) = DF(i,j,1)+hyi*F(i,my+s,4)*DL(s,j-my)/rho
            DF(i,j,2) = DF(i,j,2)+hyi*F(i,my+s,5)*DL(s,j-my)/rho
            DF(i,j,3) = DF(i,j,3)+hyi*F(i,my+s,2)*DL(s,j-my)*M*gamma
            DF(i,j,4) = DF(i,j,4)+hyi*F(i,my+s,1)*DL(s,j-my)*G
            DF(i,j,5) = DF(i,j,5)+hyi*F(i,my+s,2)*DL(s,j-my)*M

          end do
        end do
      end do

    case('LT')

      do j = my,py
        do i = mx,px
          if(pml%n > 0d0) then
            pmlsx = pml%sigma*abs(max(0d0,max(pml%mx-x(i,j),x(i,j)-pml%px)))**pml%n
            pmlsy = pml%sigma*abs(max(0d0,max(pml%my-y(i,j),y(i,j)-pml%py)))**pml%n
          end if

          if(pml%pmlx) then
            DF(i,j,1) = DF(i,j,1) - Wx(i,j,3)/rho
            DF(i,j,2) = DF(i,j,2) - Wx(i,j,4)/rho
            DF(i,j,3) = DF(i,j,3) - Wx(i,j,1)*M
            DF(i,j,4) = DF(i,j,4) - Wx(i,j,2)*G
            DF(i,j,5) = DF(i,j,5) - Wx(i,j,1)*gamma*M

            DWx(i,j,:) = DWx(i,j,:) - (pmlsx+pml%kappa)*Wx(i,j,:)
          end if

          if(pml%pmly) then
            DF(i,j,1) = DF(i,j,1) - Wy(i,j,4)/rho
            DF(i,j,2) = DF(i,j,2) - Wy(i,j,5)/rho
            DF(i,j,3) = DF(i,j,3) - Wy(i,j,2)*gamma*M
            DF(i,j,4) = DF(i,j,4) - Wy(i,j,1)*G
            DF(i,j,5) = DF(i,j,5) - Wy(i,j,2)*M

            DWy(i,j,:) = DWy(i,j,:) - (pmlsy+pml%kappa)*Wy(i,j,:)
          end if

          do s = 0,nbnd-1 ! one-sided for q and r

            if(pml%pmlx) then
              DWx(i,j,:) = DWx(i,j,:) + pmlsx*hxi*F(mx+s,j,:)*DL(s,i-mx)
            end if

            DF(i,j,1) = DF(i,j,1)+hxi*F(mx+s,j,3)*DL(s,i-mx)/rho
            DF(i,j,2) = DF(i,j,2)+hxi*F(mx+s,j,4)*DL(s,i-mx)/rho
            DF(i,j,3) = DF(i,j,3)+hxi*F(mx+s,j,1)*DL(s,i-mx)*M
            DF(i,j,4) = DF(i,j,4)+hxi*F(mx+s,j,2)*DL(s,i-mx)*G
            DF(i,j,5) = DF(i,j,5)+hxi*F(mx+s,j,1)*DL(s,i-mx)*M*gamma

            if(pml%pmly) then
              DWy(i,j,:) = DWy(i,j,:) + pmlsy*hyi*F(i,py-s,:)*DR(-s,j-py)
            end if

            DF(i,j,1) = DF(i,j,1)+hyi*F(i,py-s,4)*DR(-s,j-py)/rho
            DF(i,j,2) = DF(i,j,2)+hyi*F(i,py-s,5)*DR(-s,j-py)/rho
            DF(i,j,3) = DF(i,j,3)+hyi*F(i,py-s,2)*DR(-s,j-py)*M*gamma
            DF(i,j,4) = DF(i,j,4)+hyi*F(i,py-s,1)*DR(-s,j-py)*G
            DF(i,j,5) = DF(i,j,5)+hyi*F(i,py-s,2)*DR(-s,j-py)*M

          end do
        end do
      end do

    case('RT')

      do j = my,py
        do i = mx,px
          if(pml%n > 0d0) then
            pmlsx = pml%sigma*abs(max(0d0,max(pml%mx-x(i,j),x(i,j)-pml%px)))**pml%n
            pmlsy = pml%sigma*abs(max(0d0,max(pml%my-y(i,j),y(i,j)-pml%py)))**pml%n
          end if

          if(pml%pmlx) then
            DF(i,j,1) = DF(i,j,1) - Wx(i,j,3)/rho
            DF(i,j,2) = DF(i,j,2) - Wx(i,j,4)/rho
            DF(i,j,3) = DF(i,j,3) - Wx(i,j,1)*M
            DF(i,j,4) = DF(i,j,4) - Wx(i,j,2)*G
            DF(i,j,5) = DF(i,j,5) - Wx(i,j,1)*gamma*M

            DWx(i,j,:) = DWx(i,j,:) - (pmlsx+pml%kappa)*Wx(i,j,:)
          end if

          if(pml%pmly) then
            DF(i,j,1) = DF(i,j,1) - Wy(i,j,4)/rho
            DF(i,j,2) = DF(i,j,2) - Wy(i,j,5)/rho
            DF(i,j,3) = DF(i,j,3) - Wy(i,j,2)*gamma*M
            DF(i,j,4) = DF(i,j,4) - Wy(i,j,1)*G
            DF(i,j,5) = DF(i,j,5) - Wy(i,j,2)*M

            DWy(i,j,:) = DWy(i,j,:) - (pmlsy+pml%kappa)*Wy(i,j,:)
          end if

          do s = 0,nbnd-1 ! one-sided for q and r

            if(pml%pmlx) then
              DWx(i,j,:) = DWx(i,j,:) + pmlsx*hxi*F(px-s,j,:)*DR(-s,i-px)
            end if

            DF(i,j,1) = DF(i,j,1)+hxi*F(px-s,j,3)*DR(-s,i-px)/rho
            DF(i,j,2) = DF(i,j,2)+hxi*F(px-s,j,4)*DR(-s,i-px)/rho
            DF(i,j,3) = DF(i,j,3)+hxi*F(px-s,j,1)*DR(-s,i-px)*M
            DF(i,j,4) = DF(i,j,4)+hxi*F(px-s,j,2)*DR(-s,i-px)*G
            DF(i,j,5) = DF(i,j,5)+hxi*F(px-s,j,1)*DR(-s,i-px)*M*gamma

            if(pml%pmly) then
              DWy(i,j,:) = DWy(i,j,:) + pmlsy*hyi*F(i,py-s,:)*DR(-s,j-py)
            end if

            DF(i,j,1) = DF(i,j,1)+hyi*F(i,py-s,4)*DR(-s,j-py)/rho
            DF(i,j,2) = DF(i,j,2)+hyi*F(i,py-s,5)*DR(-s,j-py)/rho
            DF(i,j,3) = DF(i,j,3)+hyi*F(i,py-s,2)*DR(-s,j-py)*M*gamma
            DF(i,j,4) = DF(i,j,4)+hyi*F(i,py-s,1)*DR(-s,j-py)*G
            DF(i,j,5) = DF(i,j,5)+hyi*F(i,py-s,2)*DR(-s,j-py)*M

          end do
        end do
      end do

    end select

  end subroutine rates_boundary_mode2_pml


end module rates
