module rates_heterogeneous

  implicit none

contains


  subroutine rates_h_interior_mode3(lbndFx,lbndFy,lbndDFx,lbndDFy, &
       F,DF,rho,G,xq,xr,yq,yr,JacInv,mx,px,my,py)

    use fd_coeff, only : DI

    implicit none

    integer,intent(in) :: lbndFx,lbndFy,lbndDFx,lbndDFy,mx,px,my,py
    real,dimension(lbndFx :,lbndFy :,1:),intent(in) :: F
    real,dimension(lbndDFx:,lbndDFy:,1:),intent(inout) :: DF
    real,dimension(lbndFx:,lbndFy:),intent(in) :: xq,xr,yq,yr,JacInv,rho,G

    integer :: i,j,s

    do j = my,py
       do i = mx,px
          do s = 1,3

             DF(i,j,1) = DF(i,j,1)+( &
                  0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,2)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,3)- &
                  0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,2)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,3)- &
                  0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,2)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,3)+ &
                  0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,2)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,3) &
                  )*DI(s)*JacInv(i,j)/rho(i,j)

             DF(i,j,2) = DF(i,j,2)+( &
                  0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,1)-0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,1)- &
                  0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,1)+0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,1) &
                  )*DI(s)*G(i,j)*JacInv(i,j)

             DF(i,j,3) = DF(i,j,3)+( &
                 -0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,1)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,1)+ &
                  0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,1)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,1) &
                  )*DI(s)*G(i,j)*JacInv(i,j)


          end do
       end do
    end do

  end subroutine rates_h_interior_mode3


  subroutine rates_h_boundary_mode3(lbndFx,lbndFy,lbndDFx,lbndDFy, &
       F,DF,rho,G,xq,xr,yq,yr,JacInv,mx,px,my,py,bnd)

    use fd_coeff, only : DL,DR,DI,pI,nbnd

    implicit none

    integer,intent(in) :: lbndFx,lbndFy,lbndDFx,lbndDFy,mx,px,my,py
    real,dimension(lbndFx :,lbndFy :,1:),intent(in) :: F
    real,dimension(lbndDFx:,lbndDFy:,1:),intent(inout) :: DF
    real,dimension(lbndFx:,lbndFy:),intent(in) :: xq,xr,yq,yr,JacInv,rho,G
    character(*),intent(in) :: bnd

    integer :: i,j,s

    select case(bnd)

    case('L')

       do j = my,py
          do i = mx,px

             do s = 0,nbnd-1 ! one-sided for q
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,2)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,3) &
                     )*DL(s,i-mx)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,1) &
                     )*DL(s,i-mx)*G(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     -0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,1) &
                     )*DL(s,i-mx)*G(i,j)*JacInv(i,j)
                
             end do

             do s = 1,pI ! central for r
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,2)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,3)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,2)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,3) &
                     )*DI(s)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,1)+0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,1) &
                     )*DI(s)*G(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,1)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,1) &
                     )*DI(s)*G(i,j)*JacInv(i,j)

             end do

          end do
       end do

    case('R')
       
       do j = my,py
          do i = mx,px

             do s = 0,nbnd-1 ! one-sided for q
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,2)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,3) &
                     )*DR(-s,i-px)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,1) &
                     )*DR(-s,i-px)*G(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     -0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,1) &
                     )*DR(-s,i-px)*G(i,j)*JacInv(i,j)

             end do

             do s = 1,pI ! central for r
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,2)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,3)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,2)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,3) &
                     )*DI(s)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,1)+0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,1) &
                     )*DI(s)*G(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,1)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,1) &
                     )*DI(s)*G(i,j)*JacInv(i,j)
                
             end do

          end do
       end do

    case('B')
       
       do j = my,py
          do i = mx,px

             do s = 0,nbnd-1 ! one-sided for r
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,2)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,3) &
                     )*DL(s,j-my)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,1) &
                     )*DL(s,j-my)*G(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,1) &
                     )*DL(s,j-my)*G(i,j)*JacInv(i,j)
                
             end do

             do s = 1,pI ! central for q
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,2)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,3)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,2)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,3) &
                     )*DI(s)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,1)-0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,1) &
                     )*DI(s)*G(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     -0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,1)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,1) &
                     )*DI(s)*G(i,j)*JacInv(i,j)
                
             end do

          end do
       end do
       
    case('T')
    
       do j = my,py
          do i = mx,px

             do s = 0,nbnd-1 ! one-sided for r
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,2)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,3) &
                     )*DR(-s,j-py)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,1) &
                     )*DR(-s,j-py)*G(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,1) &
                     )*DR(-s,j-py)*G(i,j)*JacInv(i,j)
                
                
             end do

             do s = 1,pI ! central for q
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,2)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,3)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,2)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,3) &
                     )*DI(s)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,1)-0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,1) &
                     )*DI(s)*G(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     -0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,1)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,1) &
                     )*DI(s)*G(i,j)*JacInv(i,j)
                
                
             end do

          end do
       end do
       
    case('LB')

       do j = my,py
          do i = mx,px
             do s = 0,nbnd-1 ! one-sided for q and r
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,2)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,3) &
                     )*DL(s,i-mx)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,1) &
                     )*DL(s,i-mx)*G(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     -0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,1) &
                     )*DL(s,i-mx)*G(i,j)*JacInv(i,j)
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,2)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,3) &
                     )*DL(s,j-my)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,1) &
                     )*DL(s,j-my)*G(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,1) &
                     )*DL(s,j-my)*G(i,j)*JacInv(i,j)
                
             end do
          end do
       end do

    case('RB')

       do j = my,py
          do i = mx,px
             do s = 0,nbnd-1 ! one-sided for q and r

                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,2)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,3) &
                     )*DR(-s,i-px)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,1) &
                     )*DR(-s,i-px)*G(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     -0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,1) &
                     )*DR(-s,i-px)*G(i,j)*JacInv(i,j)
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,2)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,3) &
                     )*DL(s,j-my)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,1) &
                     )*DL(s,j-my)*G(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,1) &
                     )*DL(s,j-my)*G(i,j)*JacInv(i,j)

             end do
          end do
       end do

    case('LT')

       do j = my,py
          do i = mx,px
             do s = 0,nbnd-1 ! one-sided for q and r
              
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,2)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,3) &
                     )*DL(s,i-mx)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,1) &
                     )*DL(s,i-mx)*G(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     -0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,1) &
                     )*DL(s,i-mx)*G(i,j)*JacInv(i,j)

                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,2)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,3) &
                     )*DR(-s,j-py)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,1) &
                     )*DR(-s,j-py)*G(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,1) &
                     )*DR(-s,j-py)*G(i,j)*JacInv(i,j)
                
             end do
          end do
       end do

    case('RT')
    
       do j = my,py
          do i = mx,px
             do s = 0,nbnd-1 ! one-sided for q and r
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,2)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,3) &
                     )*DR(-s,i-px)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,1) &
                     )*DR(-s,i-px)*G(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     -0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,1) &
                     )*DR(-s,i-px)*G(i,j)*JacInv(i,j)
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,2)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,3) &
                     )*DR(-s,j-py)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,1) &
                     )*DR(-s,j-py)*G(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,1) &
                     )*DR(-s,j-py)*G(i,j)*JacInv(i,j)
                
             end do
          end do
       end do

    end select

  end subroutine rates_h_boundary_mode3


  subroutine rates_h_interior_mode2(lbndFx,lbndFy,lbndDFx,lbndDFy, &
       F,DF,rho,G,M,gamma,xq,xr,yq,yr,JacInv,mx,px,my,py)

    use fd_coeff, only : DI

    implicit none

    integer,intent(in) :: lbndFx,lbndFy,lbndDFx,lbndDFy,mx,px,my,py
    real,dimension(lbndFx :,lbndFy :,1:),intent(in) :: F
    real,dimension(lbndDFx:,lbndDFy:,1:),intent(inout) :: DF
    real,dimension(lbndFx:,lbndFy:),intent(in) :: xq,xr,yq,yr,JacInv,rho,G,M,gamma

    integer :: i,j,s

    do j = my,py
       do i = mx,px
          do s = 1,3

             DF(i,j,1) = DF(i,j,1)+( &
                  0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,3)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,4)- &
                  0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,3)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,4)- &
                  0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,3)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,4)+ &
                  0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,3)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,4) &
                  )*DI(s)*JacInv(i,j)/rho(i,j)

             DF(i,j,2) = DF(i,j,2)+( &
                  0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,4)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,5)- &
                  0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,4)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,5)- &
                  0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,4)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,5)+ &
                  0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,4)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,5) &
                  )*DI(s)*JacInv(i,j)/rho(i,j)

             DF(i,j,3) = DF(i,j,3)+( &
                  0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,1)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,2)*gamma(i,j)- &
                  0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,1)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,2)*gamma(i,j)- &
                  0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,1)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,2)*gamma(i,j)+ &
                  0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,1)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,2)*gamma(i,j) &
                  )*DI(s)*JacInv(i,j)*M(i,j)

             DF(i,j,4) = DF(i,j,4)+( &
                  0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,2)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,1)- &
                  0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,2)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,1)- &
                  0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,2)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,1)+ &
                  0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,2)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,1) &
                  )*DI(s)*JacInv(i,j)*G(i,j)

             DF(i,j,5) = DF(i,j,5)+( &
                  0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,1)*gamma(i,j)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,2)- &
                  0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,1)*gamma(i,j)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,2)- &
                  0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,1)*gamma(i,j)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,2)+ &
                  0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,1)*gamma(i,j)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,2) &
                  )*DI(s)*JacInv(i,j)*M(i,j)

          end do
       end do
    end do

  end subroutine rates_h_interior_mode2


  subroutine rates_h_boundary_mode2(lbndFx,lbndFy,lbndDFx,lbndDFy, &
       F,DF,rho,G,M,gamma,xq,xr,yq,yr,JacInv,mx,px,my,py,bnd)

    use fd_coeff, only : DL,DR,DI,pI,nbnd

    implicit none

    integer,intent(in) :: lbndFx,lbndFy,lbndDFx,lbndDFy,mx,px,my,py
    real,dimension(lbndFx :,lbndFy :,1:),intent(in) :: F
    real,dimension(lbndDFx:,lbndDFy:,1:),intent(inout) :: DF
    real,dimension(lbndFx:,lbndFy:),intent(in) :: xq,xr,yq,yr,JacInv,rho,G,M,gamma
    character(*),intent(in) :: bnd

    integer :: i,j,s

    select case(bnd)

    case('L')

       do j = my,py
          do i = mx,px

             do s = 0,nbnd-1 ! one-sided for q
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,3)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,4) &
                     )*DL(s,i-mx)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,4)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,5) &
                     )*DL(s,i-mx)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,1)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,2)*gamma(i,j) &
                     )*DL(s,i-mx)*M(i,j)*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,2)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,1) &
                     )*DL(s,i-mx)*G(i,j)*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,1)*gamma(i,j)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,2) &
                     )*DL(s,i-mx)*M(i,j)*JacInv(i,j)

             end do

             do s = 1,pI ! central for r
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,3)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,4)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,3)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,4) &
                     )*DI(s)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,4)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,5)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,4)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,5) &
                     )*DI(s)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,1)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,2)*gamma(i,j)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,1)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,2)*gamma(i,j) &
                     )*DI(s)*M(i,j)*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,2)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,1)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,2)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,1) &
                     )*DI(s)*G(i,j)*JacInv(i,j)

                DF(i,j,5) = DF(i,j,5)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,1)*gamma(i,j)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,2)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,1)*gamma(i,j)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,2) &
                     )*DI(s)*M(i,j)*JacInv(i,j)

             end do

          end do
       end do

    case('R')
       
       do j = my,py
          do i = mx,px

             do s = 0,nbnd-1 ! one-sided for q
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,3)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,4) &
                     )*DR(-s,i-px)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,4)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,5) &
                     )*DR(-s,i-px)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,1)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,2)*gamma(i,j) &
                     )*DR(-s,i-px)*M(i,j)*JacInv(i,j)

                DF(i,j,4) = DF(i,j,4)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,2)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,1) &
                     )*DR(-s,i-px)*G(i,j)*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,1)*gamma(i,j)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,2) &
                     )*DR(-s,i-px)*M(i,j)*JacInv(i,j)

             end do

             do s = 1,pI ! central for r
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,3)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,4)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,3)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,4) &
                     )*DI(s)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,4)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,5)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,4)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,5) &
                     )*DI(s)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,1)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,2)*gamma(i,j)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,1)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,2)*gamma(i,j) &
                     )*DI(s)*M(i,j)*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,2)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,1)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,2)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,1) &
                     )*DI(s)*G(i,j)*JacInv(i,j)

                DF(i,j,5) = DF(i,j,5)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,1)*gamma(i,j)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,2)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,1)*gamma(i,j)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,2) &
                     )*DI(s)*M(i,j)*JacInv(i,j)

             end do

          end do
       end do

    case('B')
       
       do j = my,py
          do i = mx,px

             do s = 0,nbnd-1 ! one-sided for r
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,3)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,4) &
                     )*DL(s,j-my)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,4)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,5) &
                     )*DL(s,j-my)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,1)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,2)*gamma(i,j) &
                     )*DL(s,j-my)*M(i,j)*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,2)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,1) &
                     )*DL(s,j-my)*G(i,j)*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,1)*gamma(i,j)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,2) &
                     )*DL(s,j-my)*M(i,j)*JacInv(i,j)
                
             end do

             do s = 1,pI ! central for q
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,3)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,4)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,3)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,4) &
                     )*DI(s)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,4)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,5)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,4)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,5) &
                     )*DI(s)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,1)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,2)*gamma(i,j)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,1)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,2)*gamma(i,j) &
                     )*DI(s)*M(i,j)*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,2)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,1)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,2)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,1) &
                     )*DI(s)*G(i,j)*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,1)*gamma(i,j)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,2)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,1)*gamma(i,j)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,2) &
                     )*DI(s)*M(i,j)*JacInv(i,j)
                
             end do

          end do
       end do
       
    case('T')
    
       do j = my,py
          do i = mx,px

             do s = 0,nbnd-1 ! one-sided for r
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,3)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,4) &
                     )*DR(-s,j-py)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,4)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,5) &
                     )*DR(-s,j-py)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,1)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,2)*gamma(i,j) &
                     )*DR(-s,j-py)*M(i,j)*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,2)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,1) &
                     )*DR(-s,j-py)*G(i,j)*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,1)*gamma(i,j)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,2) &
                     )*DR(-s,j-py)*M(i,j)*JacInv(i,j)
                
             end do

             do s = 1,pI ! central for q
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,3)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,4)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,3)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,4) &
                     )*DI(s)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,4)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,5)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,4)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,5) &
                     )*DI(s)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,1)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,2)*gamma(i,j)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,1)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,2)*gamma(i,j) &
                     )*DI(s)*M(i,j)*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,2)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,1)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,2)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,1) &
                     )*DI(s)*G(i,j)*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,1)*gamma(i,j)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,2)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,1)*gamma(i,j)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,2) &
                     )*DI(s)*M(i,j)*JacInv(i,j)
                
             end do

          end do
       end do
       
    case('LB')

       do j = my,py
          do i = mx,px
             do s = 0,nbnd-1 ! one-sided for q and r
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,3)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,4) &
                     )*DL(s,i-mx)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,4)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,5) &
                     )*DL(s,i-mx)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,1)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,2)*gamma(i,j) &
                     )*DL(s,i-mx)*M(i,j)*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,2)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,1) &
                     )*DL(s,i-mx)*G(i,j)*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,1)*gamma(i,j)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,2) &
                     )*DL(s,i-mx)*M(i,j)*JacInv(i,j)
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,3)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,4) &
                     )*DL(s,j-my)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,4)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,5) &
                     )*DL(s,j-my)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,1)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,2)*gamma(i,j) &
                     )*DL(s,j-my)*M(i,j)*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,2)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,1) &
                     )*DL(s,j-my)*G(i,j)*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,1)*gamma(i,j)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,2) &
                     )*DL(s,j-my)*M(i,j)*JacInv(i,j)
                
             end do
          end do
       end do

    case('RB')

       do j = my,py
          do i = mx,px
             do s = 0,nbnd-1 ! one-sided for q and r

                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,3)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,4) &
                     )*DR(-s,i-px)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,4)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,5) &
                     )*DR(-s,i-px)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,1)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,2)*gamma(i,j) &
                     )*DR(-s,i-px)*M(i,j)*JacInv(i,j)

                DF(i,j,4) = DF(i,j,4)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,2)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,1) &
                     )*DR(-s,i-px)*G(i,j)*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,1)*gamma(i,j)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,2) &
                     )*DR(-s,i-px)*M(i,j)*JacInv(i,j)

                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,3)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,4) &
                     )*DL(s,j-my)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,4)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,5) &
                     )*DL(s,j-my)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,1)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,2)*gamma(i,j) &
                     )*DL(s,j-my)*M(i,j)*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,2)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,1) &
                     )*DL(s,j-my)*G(i,j)*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,1)*gamma(i,j)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,2) &
                     )*DL(s,j-my)*M(i,j)*JacInv(i,j)
                
             end do
          end do
       end do

    case('LT')

       do j = my,py
          do i = mx,px
             do s = 0,nbnd-1 ! one-sided for q and r
              
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,3)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,4) &
                     )*DL(s,i-mx)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,4)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,5) &
                     )*DL(s,i-mx)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,1)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,2)*gamma(i,j) &
                     )*DL(s,i-mx)*M(i,j)*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,2)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,1) &
                     )*DL(s,i-mx)*G(i,j)*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,1)*gamma(i,j)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,2) &
                     )*DL(s,i-mx)*M(i,j)*JacInv(i,j)
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,3)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,4) &
                     )*DR(-s,j-py)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,4)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,5) &
                     )*DR(-s,j-py)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,1)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,2)*gamma(i,j) &
                     )*DR(-s,j-py)*M(i,j)*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,2)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,1) &
                     )*DR(-s,j-py)*G(i,j)*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,1)*gamma(i,j)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,2) &
                     )*DR(-s,j-py)*M(i,j)*JacInv(i,j)
                
             end do
          end do
       end do

    case('RT')
    
       do j = my,py
          do i = mx,px
             do s = 0,nbnd-1 ! one-sided for q and r
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,3)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,4) &
                     )*DR(-s,i-px)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,4)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,5) &
                     )*DR(-s,i-px)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,1)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,2)*gamma(i,j) &
                     )*DR(-s,i-px)*M(i,j)*JacInv(i,j)

                DF(i,j,4) = DF(i,j,4)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,2)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,1) &
                     )*DR(-s,i-px)*G(i,j)*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,1)*gamma(i,j)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,2) &
                     )*DR(-s,i-px)*M(i,j)*JacInv(i,j)

                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,3)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,4) &
                     )*DR(-s,j-py)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,4)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,5) &
                     )*DR(-s,j-py)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,1)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,2)*gamma(i,j) &
                     )*DR(-s,j-py)*M(i,j)*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,2)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,1) &
                     )*DR(-s,j-py)*G(i,j)*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,1)*gamma(i,j)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,2) &
                     )*DR(-s,j-py)*M(i,j)*JacInv(i,j)
                
             end do
          end do
       end do

    end select

  end subroutine rates_h_boundary_mode2


  subroutine rates_h_szz(lbndFx,lbndFy,lbndDFx,lbndDFy,DF,nu,mx,px,my,py)

    implicit none

    integer,intent(in) :: lbndFx,lbndFy,lbndDFx,lbndDFy,mx,px,my,py
    real,dimension(lbndDFx:,lbndDFy:,1:),intent(inout) :: DF
    real,dimension(lbndFx:,lbndFy:),intent(in) :: nu

    integer :: i,j

    do j = my,py
       do i = mx,px
          DF(i,j,6) = nu(i,j)*(DF(i,j,3)+DF(i,j,5))
       end do
    end do

  end subroutine rates_h_szz


  subroutine rates_h_dissipation_interior_mode3(lbndFx,lbndFy,lbndDFx,lbndDFy, &
       F,DF,rho,G,Cdiss,xq,xr,yq,yr,Jac,JacInv,mx,px,my,py)

    use fd_coeff, only : DI,AI

    implicit none

    integer,intent(in) :: lbndFx,lbndFy,lbndDFx,lbndDFy,mx,px,my,py
    real,dimension(lbndFx :,lbndFy :,1:),intent(in) :: F
    real,dimension(lbndDFx:,lbndDFy:,1:),intent(inout) :: DF
    real,intent(in) :: Cdiss
    real,dimension(lbndFx:,lbndFy:),intent(in) :: xq,xr,yq,yr,Jac,JacInv,rho,G

    integer :: i,j,s
    real :: c(0:3)

    c=0d0
    do i = 1,ubound(DI,1)
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
                  )*DI(s)/rho(i,j)*JacInv(i,j)+c(s)*JacInv(i,j)*( &
                  Jac(i+s,j)*F(i+s,j,1)+Jac(i-s,j)*F(i-s,j,1)+ &
                  Jac(i,j+s)*F(i,j+s,1)+Jac(i,j-s)*F(i,j-s,1))

             DF(i,j,2) = DF(i,j,2)+( &
                  0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,1)-0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,1)- &
                  0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,1)+0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,1) &
                  )*DI(s)*G(i,j)*JacInv(i,j)+c(s)*JacInv(i,j)*( &
                  Jac(i+s,j)*F(i+s,j,2)+Jac(i-s,j)*F(i-s,j,2)+ &
                  Jac(i,j+s)*F(i,j+s,2)+Jac(i,j-s)*F(i,j-s,2))

             DF(i,j,3) = DF(i,j,3)+( &
                 -0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,1)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,1)+ &
                  0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,1)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,1) &
                  )*DI(s)*G(i,j)*JacInv(i,j)+c(s)*JacInv(i,j)*( &
                  Jac(i+s,j)*F(i+s,j,3)+Jac(i-s,j)*F(i-s,j,3)+ &
                  Jac(i,j+s)*F(i,j+s,3)+Jac(i,j-s)*F(i,j-s,3))

          end do
       end do
    end do

  end subroutine rates_h_dissipation_interior_mode3


  subroutine rates_h_dissipation_boundary_mode3(lbndFx,lbndFy,lbndDFx,lbndDFy, &
       F,DF,rho,G,Cdiss,xq,xr,yq,yr,Jac,JacInv,mx,px,my,py,bnd)

    use fd_coeff, only : DL,DR,DI,AL,AR,AI,pI,nbnd

    implicit none

    integer,intent(in) :: lbndFx,lbndFy,lbndDFx,lbndDFy,mx,px,my,py
    real,dimension(lbndFx :,lbndFy :,1:),intent(in) :: F
    real,dimension(lbndDFx:,lbndDFy:,1:),intent(inout) :: DF
    real,intent(in) :: Cdiss
    real,dimension(lbndFx:,lbndFy:),intent(in) :: xq,xr,yq,yr,Jac,JacInv,rho,G
    character(*),intent(in) :: bnd

    integer :: i,j,s

    select case(bnd)

    case('L')

       do j = my,py
          do i = mx,px

             do s = 0,nbnd-1 ! one-sided for q
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,2)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,3) &
                     )*DL(s,i-mx)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,1) &
                     )*DL(s,i-mx)*G(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     -0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,1) &
                     )*DL(s,i-mx)*G(i,j)*JacInv(i,j)
                
                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(mx+s,j)*F(mx+s,j,:) &
                     )*AL(s,i-mx)*JacInv(i,j)*Cdiss
                
             end do

             DF(i,j,:) = DF(i,j,:)+F(i,j,:)*AI(0)*Cdiss

             do s = 1,pI ! central for r
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,2)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,3)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,2)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,3) &
                     )*DI(s)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,1)+0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,1) &
                     )*DI(s)*G(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,1)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,1) &
                     )*DI(s)*G(i,j)*JacInv(i,j)

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
                     )*DR(-s,i-px)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,1) &
                     )*DR(-s,i-px)*G(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     -0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,1) &
                     )*DR(-s,i-px)*G(i,j)*JacInv(i,j)

                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(px-s,j)*F(px-s,j,:) &
                     )*AR(-s,i-px)*JacInv(i,j)*Cdiss

             end do

             DF(i,j,:) = DF(i,j,:)+F(i,j,:)*AI(0)*Cdiss

             do s = 1,pI ! central for r
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,2)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,3)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,2)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,3) &
                     )*DI(s)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,1)+0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,1) &
                     )*DI(s)*G(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,1)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,1) &
                     )*DI(s)*G(i,j)*JacInv(i,j)
                
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
                     )*DL(s,j-my)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,1) &
                     )*DL(s,j-my)*G(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,1) &
                     )*DL(s,j-my)*G(i,j)*JacInv(i,j)
                
                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(i,my+s)*F(i,my+s,:) &
                     )*AL(s,j-my)*JacInv(i,j)*Cdiss
                
             end do

             DF(i,j,:) = DF(i,j,:)+F(i,j,:)*AI(0)*Cdiss

             do s = 1,pI ! central for q
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,2)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,3)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,2)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,3) &
                     )*DI(s)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,1)-0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,1) &
                     )*DI(s)*G(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     -0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,1)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,1) &
                     )*DI(s)*G(i,j)*JacInv(i,j)
                
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
                     )*DR(-s,j-py)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,1) &
                     )*DR(-s,j-py)*G(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,1) &
                     )*DR(-s,j-py)*G(i,j)*JacInv(i,j)
                
                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(i,py-s)*F(i,py-s,:) &
                     )*AR(-s,j-py)*JacInv(i,j)*Cdiss
                
             end do

             DF(i,j,:) = DF(i,j,:)+F(i,j,:)*AI(0)*Cdiss

             do s = 1,pI ! central for q
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,2)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,3)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,2)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,3) &
                     )*DI(s)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,1)-0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,1) &
                     )*DI(s)*G(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     -0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,1)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,1) &
                     )*DI(s)*G(i,j)*JacInv(i,j)
                
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
                     )*DL(s,i-mx)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,1) &
                     )*DL(s,i-mx)*G(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     -0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,1) &
                     )*DL(s,i-mx)*G(i,j)*JacInv(i,j)
                
                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(mx+s,j)*F(mx+s,j,:) &
                     )*AL(s,i-mx)*JacInv(i,j)*Cdiss
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,2)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,3) &
                     )*DL(s,j-my)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,1) &
                     )*DL(s,j-my)*G(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,1) &
                     )*DL(s,j-my)*G(i,j)*JacInv(i,j)
                
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
                     )*DR(-s,i-px)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,1) &
                     )*DR(-s,i-px)*G(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     -0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,1) &
                     )*DR(-s,i-px)*G(i,j)*JacInv(i,j)
                
                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(px-s,j)*F(px-s,j,:) &
                     )*AR(-s,i-px)*JacInv(i,j)*Cdiss
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,2)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,3) &
                     )*DL(s,j-my)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,1) &
                     )*DL(s,j-my)*G(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,1) &
                     )*DL(s,j-my)*G(i,j)*JacInv(i,j)

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
                     )*DL(s,i-mx)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,1) &
                     )*DL(s,i-mx)*G(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     -0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,1) &
                     )*DL(s,i-mx)*G(i,j)*JacInv(i,j)

                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(mx+s,j)*F(mx+s,j,:) &
                     )*AL(s,i-mx)*JacInv(i,j)*Cdiss

                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,2)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,3) &
                     )*DR(-s,j-py)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,1) &
                     )*DR(-s,j-py)*G(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,1) &
                     )*DR(-s,j-py)*G(i,j)*JacInv(i,j)
                
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
                     )*DR(-s,i-px)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,1) &
                     )*DR(-s,i-px)*G(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     -0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,1) &
                     )*DR(-s,i-px)*G(i,j)*JacInv(i,j)
                
                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(px-s,j)*F(px-s,j,:) &
                     )*AR(-s,i-px)*JacInv(i,j)*Cdiss
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,2)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,3) &
                     )*DR(-s,j-py)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,1) &
                     )*DR(-s,j-py)*G(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,1) &
                     )*DR(-s,j-py)*G(i,j)*JacInv(i,j)
                
                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(i,py-s)*F(i,py-s,:) &
                     )*AR(-s,j-py)*JacInv(i,j)*Cdiss
                
             end do
          end do
       end do

    end select

  end subroutine rates_h_dissipation_boundary_mode3


  subroutine rates_h_dissipation_interior_mode2(lbndFx,lbndFy,lbndDFx,lbndDFy, &
       F,DF,rho,G,M,gamma,Cdiss,xq,xr,yq,yr,Jac,JacInv,mx,px,my,py)

    use fd_coeff, only : DI,AI

    implicit none

    integer,intent(in) :: lbndFx,lbndFy,lbndDFx,lbndDFy,mx,px,my,py
    real,dimension(lbndFx :,lbndFy :,1:),intent(in) :: F
    real,dimension(lbndDFx:,lbndDFy:,1:),intent(inout) :: DF
    real,intent(in) :: Cdiss
    real,dimension(lbndFx:,lbndFy:),intent(in) :: xq,xr,yq,yr,Jac,JacInv,rho,G,M,gamma

    integer :: i,j,s
    real :: d(0:3)

    d=0d0
    do i = 1,ubound(DI,1)
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

             DF(i,j,1) = DF(i,j,1)+(DI(s)/rho(i,j)*( &
                  0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,3)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,4)- &
                  0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,3)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,4)- &
                  0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,3)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,4)+ &
                  0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,3)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,4) &
                  )+d(s)*( &
                  Jac(i+s,j)*F(i+s,j,1)+Jac(i-s,j)*F(i-s,j,1)+ &
                  Jac(i,j+s)*F(i,j+s,1)+Jac(i,j-s)*F(i,j-s,1)))*JacInv(i,j)

             DF(i,j,2) = DF(i,j,2)+(DI(s)/rho(i,j)*( &
                  0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,4)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,5)- &
                  0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,4)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,5)- &
                  0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,4)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,5)+ &
                  0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,4)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,5) &
                  )+d(s)*( &
                  Jac(i+s,j)*F(i+s,j,2)+Jac(i-s,j)*F(i-s,j,2)+ &
                  Jac(i,j+s)*F(i,j+s,2)+Jac(i,j-s)*F(i,j-s,2)))*JacInv(i,j)

             DF(i,j,3) = DF(i,j,3)+(DI(s)*M(i,j)*( &
                  0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,1)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,2)*gamma(i,j)- &
                  0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,1)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,2)*gamma(i,j)- &
                  0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,1)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,2)*gamma(i,j)+ &
                  0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,1)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,2)*gamma(i,j) &
                  )+d(s)*( &
                  Jac(i+s,j)*F(i+s,j,3)+Jac(i-s,j)*F(i-s,j,3)+ &
                  Jac(i,j+s)*F(i,j+s,3)+Jac(i,j-s)*F(i,j-s,3)))*JacInv(i,j)

             DF(i,j,4) = DF(i,j,4)+(DI(s)*G(i,j)*( &
                  0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,2)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,1)- &
                  0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,2)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,1)- &
                  0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,2)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,1)+ &
                  0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,2)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,1) &
                  )+d(s)*( &
                  Jac(i+s,j)*F(i+s,j,4)+Jac(i-s,j)*F(i-s,j,4)+ &
                  Jac(i,j+s)*F(i,j+s,4)+Jac(i,j-s)*F(i,j-s,4)))*JacInv(i,j)

             DF(i,j,5) = DF(i,j,5)+(DI(s)*M(i,j)*( &
                  0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,1)*gamma(i,j)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,2)- &
                  0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,1)*gamma(i,j)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,2)- &
                  0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,1)*gamma(i,j)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,2)+ &
                  0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,1)*gamma(i,j)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,2) &
                  )+d(s)*( &
                  Jac(i+s,j)*F(i+s,j,5)+Jac(i-s,j)*F(i-s,j,5)+ &
                  Jac(i,j+s)*F(i,j+s,5)+Jac(i,j-s)*F(i,j-s,5)))*JacInv(i,j)
             
          end do
       end do
    end do

  end subroutine rates_h_dissipation_interior_mode2


  subroutine rates_h_dissipation_boundary_mode2(lbndFx,lbndFy,lbndDFx,lbndDFy, &
       F,DF,rho,G,M,gamma,Cdiss,xq,xr,yq,yr,Jac,JacInv,mx,px,my,py,bnd)

    use fd_coeff, only : DL,DR,DI,AL,AR,AI,pI,nbnd

    implicit none

    integer,intent(in) :: lbndFx,lbndFy,lbndDFx,lbndDFy,mx,px,my,py
    real,dimension(lbndFx :,lbndFy :,1:),intent(in) :: F
    real,dimension(lbndDFx:,lbndDFy:,1:),intent(inout) :: DF
    real,intent(in) :: Cdiss
    real,dimension(lbndFx:,lbndFy:),intent(in) :: xq,xr,yq,yr,Jac,JacInv,rho,G,M,gamma
    character(*),intent(in) :: bnd

    integer :: i,j,s

    select case(bnd)

    case('L')

       do j = my,py
          do i = mx,px

             do s = 0,nbnd-1 ! one-sided for q
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,3)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,4) &
                     )*DL(s,i-mx)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,4)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,5) &
                     )*DL(s,i-mx)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,1)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,2)*gamma(i,j) &
                     )*DL(s,i-mx)*M(i,j)*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,2)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,1) &
                     )*DL(s,i-mx)*G(i,j)*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,1)*gamma(i,j)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,2) &
                     )*DL(s,i-mx)*M(i,j)*JacInv(i,j)

                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(mx+s,j)*F(mx+s,j,:) &
                     )*AL(s,i-mx)*JacInv(i,j)*Cdiss
                
             end do

             DF(i,j,:) = DF(i,j,:)+F(i,j,:)*AI(0)*Cdiss

             do s = 1,pI ! central for r
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,3)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,4)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,3)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,4) &
                     )*DI(s)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,4)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,5)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,4)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,5) &
                     )*DI(s)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,1)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,2)*gamma(i,j)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,1)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,2)*gamma(i,j) &
                     )*DI(s)*M(i,j)*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,2)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,1)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,2)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,1) &
                     )*DI(s)*G(i,j)*JacInv(i,j)

                DF(i,j,5) = DF(i,j,5)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,1)*gamma(i,j)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,2)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,1)*gamma(i,j)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,2) &
                     )*DI(s)*M(i,j)*JacInv(i,j)

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
                     )*DR(-s,i-px)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,4)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,5) &
                     )*DR(-s,i-px)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,1)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,2)*gamma(i,j) &
                     )*DR(-s,i-px)*M(i,j)*JacInv(i,j)

                DF(i,j,4) = DF(i,j,4)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,2)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,1) &
                     )*DR(-s,i-px)*G(i,j)*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,1)*gamma(i,j)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,2) &
                     )*DR(-s,i-px)*M(i,j)*JacInv(i,j)

                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(px-s,j)*F(px-s,j,:) &
                     )*AR(-s,i-px)*JacInv(i,j)*Cdiss

             end do

             DF(i,j,:) = DF(i,j,:)+F(i,j,:)*AI(0)*Cdiss

             do s = 1,pI ! central for r
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,3)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,4)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,3)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,4) &
                     )*DI(s)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,4)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,5)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,4)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,5) &
                     )*DI(s)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,1)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,2)*gamma(i,j)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,1)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,2)*gamma(i,j) &
                     )*DI(s)*M(i,j)*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,2)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,1)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,2)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,1) &
                     )*DI(s)*G(i,j)*JacInv(i,j)

                DF(i,j,5) = DF(i,j,5)+( &
                    -0.5d0*(yq(i,j+s)+yq(i,j))*F(i,j+s,1)*gamma(i,j)+0.5d0*(xq(i,j+s)+xq(i,j))*F(i,j+s,2)+ &
                     0.5d0*(yq(i,j-s)+yq(i,j))*F(i,j-s,1)*gamma(i,j)-0.5d0*(xq(i,j-s)+xq(i,j))*F(i,j-s,2) &
                     )*DI(s)*M(i,j)*JacInv(i,j)

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
                     )*DL(s,j-my)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,4)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,5) &
                     )*DL(s,j-my)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,1)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,2)*gamma(i,j) &
                     )*DL(s,j-my)*M(i,j)*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,2)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,1) &
                     )*DL(s,j-my)*G(i,j)*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,1)*gamma(i,j)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,2) &
                     )*DL(s,j-my)*M(i,j)*JacInv(i,j)
                
                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(i,my+s)*F(i,my+s,:) &
                     )*AL(s,j-my)*JacInv(i,j)*Cdiss
                
             end do

             DF(i,j,:) = DF(i,j,:)+F(i,j,:)*AI(0)*Cdiss

             do s = 1,pI ! central for q
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,3)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,4)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,3)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,4) &
                     )*DI(s)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,4)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,5)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,4)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,5) &
                     )*DI(s)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,1)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,2)*gamma(i,j)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,1)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,2)*gamma(i,j) &
                     )*DI(s)*M(i,j)*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,2)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,1)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,2)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,1) &
                     )*DI(s)*G(i,j)*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,1)*gamma(i,j)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,2)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,1)*gamma(i,j)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,2) &
                     )*DI(s)*M(i,j)*JacInv(i,j)
                
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
                     )*DR(-s,j-py)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,4)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,5) &
                     )*DR(-s,j-py)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,1)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,2)*gamma(i,j) &
                     )*DR(-s,j-py)*M(i,j)*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,2)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,1) &
                     )*DR(-s,j-py)*G(i,j)*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,1)*gamma(i,j)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,2) &
                     )*DR(-s,j-py)*M(i,j)*JacInv(i,j)
                
                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(i,py-s)*F(i,py-s,:) &
                     )*AR(-s,j-py)*JacInv(i,j)*Cdiss
                
             end do

             DF(i,j,:) = DF(i,j,:)+F(i,j,:)*AI(0)*Cdiss

             do s = 1,pI ! central for q
                
                DF(i,j,1) = DF(i,j,1)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,3)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,4)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,3)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,4) &
                     )*DI(s)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,4)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,5)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,4)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,5) &
                     )*DI(s)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,1)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,2)*gamma(i,j)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,1)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,2)*gamma(i,j) &
                     )*DI(s)*M(i,j)*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,2)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,1)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,2)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,1) &
                     )*DI(s)*G(i,j)*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                     0.5d0*(yr(i+s,j)+yr(i,j))*F(i+s,j,1)*gamma(i,j)-0.5d0*(xr(i+s,j)+xr(i,j))*F(i+s,j,2)- &
                     0.5d0*(yr(i-s,j)+yr(i,j))*F(i-s,j,1)*gamma(i,j)+0.5d0*(xr(i-s,j)+xr(i,j))*F(i-s,j,2) &
                     )*DI(s)*M(i,j)*JacInv(i,j)
                
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
                     )*DL(s,i-mx)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,4)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,5) &
                     )*DL(s,i-mx)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,1)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,2)*gamma(i,j) &
                     )*DL(s,i-mx)*M(i,j)*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,2)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,1) &
                     )*DL(s,i-mx)*G(i,j)*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,1)*gamma(i,j)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,2) &
                     )*DL(s,i-mx)*M(i,j)*JacInv(i,j)
                
                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(mx+s,j)*F(mx+s,j,:) &
                     )*AL(s,i-mx)*JacInv(i,j)*Cdiss
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,3)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,4) &
                     )*DL(s,j-my)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,4)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,5) &
                     )*DL(s,j-my)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,1)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,2)*gamma(i,j) &
                     )*DL(s,j-my)*M(i,j)*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,2)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,1) &
                     )*DL(s,j-my)*G(i,j)*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,1)*gamma(i,j)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,2) &
                     )*DL(s,j-my)*M(i,j)*JacInv(i,j)
                
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
                     )*DR(-s,i-px)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,4)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,5) &
                     )*DR(-s,i-px)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,1)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,2)*gamma(i,j) &
                     )*DR(-s,i-px)*M(i,j)*JacInv(i,j)

                DF(i,j,4) = DF(i,j,4)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,2)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,1) &
                     )*DR(-s,i-px)*G(i,j)*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,1)*gamma(i,j)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,2) &
                     )*DR(-s,i-px)*M(i,j)*JacInv(i,j)

                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(px-s,j)*F(px-s,j,:) &
                     )*AR(-s,i-px)*JacInv(i,j)*Cdiss
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,3)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,4) &
                     )*DL(s,j-my)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,4)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,5) &
                     )*DL(s,j-my)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,1)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,2)*gamma(i,j) &
                     )*DL(s,j-my)*M(i,j)*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,2)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,1) &
                     )*DL(s,j-my)*G(i,j)*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                    -0.5d0*(yq(i,my+s)+yq(i,j))*F(i,my+s,1)*gamma(i,j)+0.5d0*(xq(i,my+s)+xq(i,j))*F(i,my+s,2) &
                     )*DL(s,j-my)*M(i,j)*JacInv(i,j)
                
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
                     )*DL(s,i-mx)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,4)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,5) &
                     )*DL(s,i-mx)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,1)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,2)*gamma(i,j) &
                     )*DL(s,i-mx)*M(i,j)*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,2)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,1) &
                     )*DL(s,i-mx)*G(i,j)*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                     0.5d0*(yr(mx+s,j)+yr(i,j))*F(mx+s,j,1)*gamma(i,j)-0.5d0*(xr(mx+s,j)+xr(i,j))*F(mx+s,j,2) &
                     )*DL(s,i-mx)*M(i,j)*JacInv(i,j)
                
                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(mx+s,j)*F(mx+s,j,:) &
                     )*AL(s,i-mx)*JacInv(i,j)*Cdiss

                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,3)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,4) &
                     )*DR(-s,j-py)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,4)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,5) &
                     )*DR(-s,j-py)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,1)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,2)*gamma(i,j) &
                     )*DR(-s,j-py)*M(i,j)*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,2)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,1) &
                     )*DR(-s,j-py)*G(i,j)*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,1)*gamma(i,j)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,2) &
                     )*DR(-s,j-py)*M(i,j)*JacInv(i,j)
                
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
                     )*DR(-s,i-px)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,4)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,5) &
                     )*DR(-s,i-px)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,1)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,2)*gamma(i,j) &
                     )*DR(-s,i-px)*M(i,j)*JacInv(i,j)

                DF(i,j,4) = DF(i,j,4)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,2)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,1) &
                     )*DR(-s,i-px)*G(i,j)*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                     0.5d0*(yr(px-s,j)+yr(i,j))*F(px-s,j,1)*gamma(i,j)-0.5d0*(xr(px-s,j)+xr(i,j))*F(px-s,j,2) &
                     )*DR(-s,i-px)*M(i,j)*JacInv(i,j)

                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(px-s,j)*F(px-s,j,:) &
                     )*AR(-s,i-px)*JacInv(i,j)*Cdiss
                
                DF(i,j,1) = DF(i,j,1)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,3)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,4) &
                     )*DR(-s,j-py)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,2) = DF(i,j,2)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,4)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,5) &
                     )*DR(-s,j-py)/rho(i,j)*JacInv(i,j)
                
                DF(i,j,3) = DF(i,j,3)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,1)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,2)*gamma(i,j) &
                     )*DR(-s,j-py)*M(i,j)*JacInv(i,j)
                
                DF(i,j,4) = DF(i,j,4)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,2)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,1) &
                     )*DR(-s,j-py)*G(i,j)*JacInv(i,j)
                
                DF(i,j,5) = DF(i,j,5)+( &
                    -0.5d0*(yq(i,py-s)+yq(i,j))*F(i,py-s,1)*gamma(i,j)+0.5d0*(xq(i,py-s)+xq(i,j))*F(i,py-s,2) &
                     )*DR(-s,j-py)*M(i,j)*JacInv(i,j)
                
                DF(i,j,:) = DF(i,j,:)+( &
                     Jac(i,py-s)*F(i,py-s,:) &
                     )*AR(-s,j-py)*JacInv(i,j)*Cdiss
                
             end do
          end do
       end do

    end select

  end subroutine rates_h_dissipation_boundary_mode2


end module rates_heterogeneous
