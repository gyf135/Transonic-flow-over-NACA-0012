! A r t i f i c i a l D i f f u s i v i t y Subroutines
subroutine nui_diff(ni ,nj ,U ,nui)
implicit none
integer :: i ,j ,ni ,nj ,ip ,im
real, dimension(4 ,ni ,nj ) :: U
real, dimension (ni ,nj ) :: nui
real :: pres_func
nui = 0.
do i = 1, ni
do j = 1, nj
ip = i + 1
im = i - 1
if (i .EQ. 1) then
im = ni
else if ( i .EQ. ni ) then
ip = 1
end if
nui(i ,j) = abs(pres_func(U(: ,ip ,j)) - 2. * pres_func(U(: ,i ,j)) + pres_func(U(: ,im ,j))) &
          / ( pres_func(U(: ,ip ,j)) + 2. * pres_func(U(: ,i ,j)) + pres_func(U(: ,im ,j)))
enddo
enddo
end subroutine nui_diff
