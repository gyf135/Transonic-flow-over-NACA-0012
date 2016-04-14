subroutine nuj_diff(ni ,nj ,U ,nuj)
implicit none
integer :: i ,j ,ni ,nj
real, dimension(4 , ni , nj ) :: U
real, dimension (ni ,nj + 1) :: nuj ! Add a d d i t i o n a l point f o r j+2
real :: pres_func
nuj = 0.
do i = 1 ,ni
do j = 2 ,nj-1
nuj(i ,j) = abs(pres_func(U(: ,i ,j+1)) - 2. * pres_func(U(: ,i ,j)) + pres_func(U(: ,i ,j-1))) &
          / ( pres_func(U(: ,i ,j+1)) + 2. * pres_func(U(: ,i ,j)) + pres_func(U(: ,i ,j-1)))
enddo
enddo
end subroutine nuj_diff

