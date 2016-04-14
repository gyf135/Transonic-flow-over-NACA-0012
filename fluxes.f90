
subroutine fluxes ( ni , nj , U, F, G)
implicit none
integer :: i , j , ni , nj
real , parameter :: gam = 1.4
real , dimension (4 , ni , nj ) :: U, F, G
do i = 1 , ni
do j = 1 , nj
! Assign F Flux Vector
F(1 , i , j ) = U(2 , i , j )
F(2 , i , j ) = (U(2 , i , j )**2 / U(1 , i , j )) + (gam - 1.) * (U(4 , i , j ) &
              - 0.5 / (U(1 , i , j ) ) * (U(2 , i , j )**2 + U(3 , i , j ) ** 2))
F(3 , i , j ) = U(2 , i , j ) * U(3 , i , j ) / U(1 , i , j)
F(4 , i , j ) = U(2 , i , j ) * (gam * U(4 , i , j ) / U(1 , i , j) &
              - (gam - 1.) / 2. / (U(1 , i , j )**2) * (U(2 , i , j )**2 + U(3 , i , j )**2))
! Assign G Flux Vector
G(1 , i , j ) = U(3 , i , j )
G(2 , i , j ) = U(2 , i , j ) * U(3 , i , j ) / U(1 , i , j )
G(3 , i , j ) = (U(3 , i , j )**2 / U(1 , i , j )) + (gam - 1.) * (U(4 , i , j ) &
              - 1. / ( 2. * U(1 , i , j )) * (U(2 , i , j )**2 + U(3 , i , j )**2) )
G(4 , i , j ) = U(3 , i , j ) * (gam * U(4 , i , j ) / U(1 , i , j ) &
              - (gam - 1.) / 2. / (U(1 , i , j )**2) * (U(2 , i , j )**2 + U(3 , i , j )**2))
enddo
enddo
end subroutine
