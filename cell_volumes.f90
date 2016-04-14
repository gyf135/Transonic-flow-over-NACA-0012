subroutine cell_volumes ( ni , nj , x , y , x_cell , y_cell , sigma , n1 , n2 , n3 , n4)
implicit none
integer :: i , j , ip , im
integer :: ni , nj
real , dimension ( ni , nj ) :: x , y , sigma , x_cell , y_cell
real , dimension (2 , ni , nj ) :: n1 , n2 , n3 , n4 ! Cell normals
x_cell =0.
y_cell =0.
sigma =0.
n1 = 0.
n2 = 0.
n3 = 0.
n4 = 0.
do i = 1 , ni
do j = 1 , nj
ip = i + 1
im = i - 1
if ( i .EQ. 1) then
im = ni
else if ( i .EQ. ni ) then
ip = 1
end if
if ( j .LT. nj ) then
! Assume x ( i , j ) i s point D ( lower l e f t c e l l point )
! The i , j of the c e l l i s a c t u a l l y at the c e l l center ( c oord inat es (x_cell , y_cell ) )
! Number of c e l l s : ( ni −1, nj ) , i n c l u d e s ’ ghost c e l l s ’ ( i , nj )
x_cell( i , j ) = 1. / 4. * (x( i , j ) + x( ip , j ) + x( i , j + 1) + x( ip , j + 1))
y_cell( i , j ) = 1. / 4. * ( y( i , j ) + y( ip , j ) + y( i , j + 1) + y( ip , j + 1) )
sigma( i , j ) = 0.5 * abs ( ( x( ip , j ) - x( i , j +1) ) * ( y( i , j ) - y( ip , j + 1) ) &
- ( x( i , j ) - x( ip , j + 1) ) * ( y( ip , j ) - y( i , j + 1)))
! Calculate the normals
n1(1 , i , j ) = y( ip , j + 1) - y( ip , j ) ! ip /2 , j Side
n1(2 , i , j ) = -(x( ip , j + 1) - x( ip , j ) )
n2(1 , i , j ) = y( i , j + 1) - y( ip , j +1) ! i , j +1/2 Side
n2(2 , i , j ) = -(x( i , j + 1) - x( ip , j + 1) )
n3(1 , i , j ) = y( i , j ) - y( i , j + 1) ! im/2 , j Side
n3(2 , i , j ) = -(x( i , j ) - x( i , j + 1) )
end if
n4(1 , i , j ) = y( ip , j ) - y( i , j ) ! i , j −1/2 Side
n4(2 , i , j ) = -(x( ip , j ) - x ( i , j ) )
enddo
enddo
end subroutine cell_volumes
