subroutine risdual ( ni , nj , U, U_cur , F ,G, n1 , n2 , n3 , n4 , dt , sigma , alpha , Risd , U_next , x_cell , y_cell , tau )
implicit none
integer :: i , j , ip , im , ni , nj , ip2 , im2
real, dimension (4) :: Fdotn_bc
real, dimension (ni , nj) :: sigma , x_cell , y_cell
real, dimension (4 ,ni ,nj) :: U, U_cur, F, G, Risd, U_next
real, dimension (2 ,ni ,nj) :: n1, n2, n3, n4 ! C e l l normals
real :: alpha , dt , pres_func , vx_func , vy_func , c_func , &
        u_plus_c_ds_func
real, dimension (ni, nj) :: tau
real, dimension (4) :: u_plus_c_ds, eps2_ds , eps4_ds
real, dimension (4, 4) :: D_ds
real :: alpha_diff2 , alpha_diff4
real, dimension (ni, nj) :: nui
real, dimension (ni , nj+1) :: nuj
Risd = 0.
alpha_diff2 = 1. / 4.0
alpha_diff4 = 1. / 50.
U_next = U_cur
D_ds = 0.
call nui_diff(ni ,nj ,U , nui)
call nuj_diff(ni ,nj ,U ,nuj)
do i = 1, ni
do j = 1, nj-1 ! Don ’ t i n c l u d e gho st c e l l
if (isnan (U(1 ,i ,j))) then
write (*, *)  "U1 is nan " , i , j
stop
end if
if ( isnan (U(2 ,i ,j ))) then
write (*, *) "U2 is nan " , i , j
stop
end if
if ( isnan (U( 3 , i , j ))) then
write (*, *) "U3 is nan " , i , j
stop
endif
if ( isnan (U(4 ,i ,j ))) then
write (* , *) "U4 is nan" , i , j
stop
endif
ip = i + 1
im = i - 1
ip2 = i + 2
im2 = i - 2
if ( i .EQ. 1) then
im = ni
im2 = ni - 1
else if ( i .EQ. ni ) then
ip = 1
ip2 = 2
else if ( i .EQ. 2) then
im2 = ni
else if ( i .EQ. ni - 1) then
ip2 = 1
endif
! For each cell sum up the 4 flux contributions
!( Note the appropriate minus sign is already included in the normal n#y)
Risd(: , i , j ) = 0.5 * (F(: , ip , j ) + F(: , i , j ) ) * n1(1 , i , j ) & ! i +1/2, j Side
                 + 0.5 * (G(: , ip , j ) + G(: , i , j ) ) * n1(2 , i , j ) &
                 + 0.5 * (F(: , i , j+1) + F(: , i , j ) ) * n2(1 , i , j ) & ! i , j+1/2 Side
                 + 0.5 * (G(: , i , j+1) + G(: , i , j ) ) * n2(2 , i , j ) &
                 + 0.5 * (F(: ,im , j ) + F(: , i , j ) ) * n3(1 , i , j ) & ! i −1/2, j Side
                 + 0.5 * (G(: ,im , j ) + G(: , i , j ) ) * n3(2 , i , j )
! i+1/2 and i −1/2 sides , a r t i f i c i a l diffusivity
u_plus_c_ds(1) = u_plus_c_ds_func(ni ,nj ,U ,n1 ,i ,ip ,j ,j)
u_plus_c_ds(3) = u_plus_c_ds_func(ni ,nj ,U ,n3 ,i ,im ,j ,j)
eps2_ds(1) = 0.5 * alpha_diff2 * u_plus_c_ds(1) &
           * max( nui(im , j ) , nui( i , j ) , nui(ip , j ) , nui( ip2 , j ))
eps2_ds(3) = 0.5 * alpha_diff2 * u_plus_c_ds(3) &
           * max( nui(im2 ,j) , nui(im ,j) , nui(i ,j) , nui(ip ,j))
eps4_ds(1) = max(0. , 0.5 * alpha_diff4 * u_plus_c_ds(1) - eps2_ds(1))
eps4_ds(3) = max(0. , 0.5 * alpha_diff4 * u_plus_c_ds(3) - eps2_ds(3))
D_ds(: ,1) = eps2_ds(1) * (U(: ,ip ,j) - U(: ,i ,j)) &
           - eps4_ds(1) * (U(: ,ip2 ,j) - 3. * U(: ,ip ,j) + 3. * U(: ,i ,j) - U(: ,im ,j))
D_ds(: ,3) = eps2_ds(3) * (U(: ,i ,j) - U(: ,im ,j)) &
           - eps4_ds(3) * (U(: ,ip ,j) - 3. * U(: ,i ,j) + 3. * U(: ,im ,j) - U(: ,im2,j))
Risd (: ,i ,j) = Risd(: ,i ,j) - (D_ds(: ,1) - D_ds(: ,3))
! Wall Boundary Condition
! Set F. n = [0 p∗nx p∗ny 0] at j = 1 at edge ( i , j −1/2) , i = 1 ,2 , . . . ,ni −1 (?)
if ( j .EQ. 1 ) then
!May need to invert the components of normal
Fdotn_bc(1) = 0.
Fdotn_bc(2) = pres_func(U_cur(: ,i ,j)) * n4 (1 ,i ,j)
Fdotn_bc(3) = pres_func(U_cur(: ,i ,j)) * n4 (2 ,i ,j)
Fdotn_bc (4) = 0.
Risd(: ,i ,j) = Risd(: ,i ,j) + Fdotn_bc (:) ! i , j −1/2 Side
else
Risd(: ,i ,j) = Risd (: ,i ,j) &
              + 0.5 * (F(: ,i ,j-1) + F(: ,i ,j)) * n4(1 ,i ,j) & ! i , j −1/2 Side
              + 0.5 * (G(: ,i ,j-1) + G(: ,i ,j)) * n4(2 ,i ,j)
u_plus_c_ds(2) = u_plus_c_ds_func(ni ,nj ,U ,n2 ,i ,i ,j ,j+1)
u_plus_c_ds(4) = u_plus_c_ds_func(ni ,nj ,U ,n4 ,i ,i ,j ,j-1)
eps2_ds(2) = 0.5 * alpha_diff2 * u_plus_c_ds(2) &
           * max(nuj(i ,j-1) ,nuj(i ,j) ,nuj(i ,j+1) , nuj(i ,j+2))
eps4_ds(2) = max( 0. ,0.5 * alpha_diff4 * u_plus_c_ds(2) - eps2_ds(2))
if ( j .EQ. 2) then
eps2_ds(4) = 0.5 * alpha_diff2 * u_plus_c_ds(4) &
           * max(nuj(i ,j-1) ,nuj(i ,j) ,nuj( i ,j+1))
eps4_ds(4) = max(0. ,0.5 * alpha_diff4 * u_plus_c_ds(4) - eps2_ds(4))
D_ds(: ,4) = eps2_ds(4) * (U(: ,i ,j) - U(: ,i ,j-1)) &
           - eps4_ds(4) * (U(: ,i ,j+1) - 3. * U(: ,i ,j) + 3. * U(: ,i ,j-1))
else
eps2_ds(4) = 0.5 * alpha_diff2 * u_plus_c_ds(4) &
           * max(nuj(i ,j-2), nuj(i ,j-1), nuj(i ,j) ,nuj(i ,j+1))
eps4_ds(4) = max(0. ,0.5 * alpha_diff4 * u_plus_c_ds(4) - eps2_ds(4))
D_ds(: ,4) = eps2_ds(4) * (U(: ,i ,j) - U(: ,i ,j-1)) &
           - eps4_ds(4) * (U(: ,i ,j+1) - 3. * U(: ,i ,j ) + 3. * U(: ,i ,j-1) - U(: ,i ,j-2))
endif
if (j .EQ. nj-1) then
D_ds(: ,2) = eps2_ds(2) * (U(: ,i ,j+1) - U(: ,i ,j))
else
D_ds(: ,2) = eps2_ds(2) * (U(: ,i ,j+1) - U(: ,i ,j)) &
           - eps4_ds(2) * (U(: ,i ,j+2) - 3. * U(: ,i ,j+1) + 3. * U(: ,i ,j ) - U(: ,i ,j-1))
endif
Risd(: ,i ,j) = Risd(: ,i ,j) - (D_ds(: ,2) - D_ds(: ,4 ))
endif
! Calculate the U value using given residual
U_next(: ,i ,j ) = U(: ,i ,j ) - tau(i ,j) * alpha * Risd (: ,i ,j)
if (isnan(U(1 ,i ,j))) then
write (*, *) "U1 is nan " , i , j
stop
endif
if (isnan(U(2 ,i ,j))) then
write (*, *) "U2 is nan " , i , j
stop
endif
if (isnan(U(3 ,i ,j ))) then
write (*, *) "U3 is nan " , i , j
stop
endif
if (isnan(U(4 ,i ,j ))) then
write ( * , *) "U4 is nan " , i , j
stop
endif
enddo
enddo
end subroutine risdual

