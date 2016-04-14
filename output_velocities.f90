subroutine output_velocities (ni ,nj ,x_cell ,y_cell ,U ,F ,G ,n4 ,p_inf ,rho_inf ,vx_inf ,vy_inf)
implicit none
integer :: i ,j ,ni ,nj
real :: vx ,vy ,speed
real, dimension( ni ,nj ) :: x_cell ,y_cell
real, dimension(4 ,ni ,nj) :: U ,F ,G
real, dimension(2 ,ni ,nj) :: n4 ! Cell normals
real :: e_nx ,e_ny
real :: pres_func ,M_func ,H_func ,S_func ,rho_func ,cp_func
real :: p_inf ,rho_inf ,vx_inf ,vy_inf ! Used in cp calculation
open( unit = 3 , file = "v_2alpha.dat") ;
! Output the velocity magnitude and direction
do i = 1 , ni
do j = 1 , nj !Don¡¯ t want to include the outer ghost cell
vx = U(2 ,i ,j) / U(1 ,i ,j)
vy = U(3 ,i ,j) / U(1 ,i ,j)
speed = (vx**2 + vy**2)**0.5
e_nx = n4(1 ,i ,j) / (n4(1 ,i ,j)**2 + n4(2 ,i ,j)**2)**0.5
e_ny = n4(2 ,i ,j) / (n4(1 ,i ,j)**2 + n4(2 ,i ,j)**2)**0.5
if (j .EQ. nj ) then
endif

if (j .LE. nj-1) then ! Don ¡¯ t include outer ghost c e l l
write (3 ,*) i ,j ,x_cell(i ,j) ,y_cell(i ,j) ,vx ,vy ,speed , &
             pres_func(U( : ,i ,j)) ,cp_func(U( : ,i ,j) ,p_inf ,rho_inf ,vx_inf ,vy_inf ) , &
             M_func(U(: ,i ,j)) ,H_func(U(: ,i ,j)) ,S_func(U(: ,i ,j)), &
             rho_func (U(: ,i ,j))
endif
enddo
enddo
end subroutine output_velocities

