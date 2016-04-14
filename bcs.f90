subroutine inlet_outlet_bcs (ni ,nj ,n4 ,vx_inf ,vy_inf ,p_inf ,rho_inf ,c_inf ,U)
implicit none
integer :: i ,j ,ni ,nj
real, dimension(4 ,ni ,nj) :: U ,F ,G
real, dimension(2 ,ni ,nj) :: n4 ! Outer Ghost Cell normal
real, parameter :: gam = 1.4
real :: vx_inf ,vy_inf ,p_inf ,rho_inf ,c_inf
real :: pres_func ,rho_func ,c_func
! Variables f o r BC’ s
real :: v_inf_dot_n ,e_nx ,e_ny ,e_tx ,e_ty ,p_i ,rho_i ,c_i ,u_n_i ,rni ,rnb
real :: u_n_b,c_b ,u_t_b ,rho_b ,p_b ,E_b
j = nj
do i = 1 , ni
v_inf_dot_n = vx_inf * n4 (1 , i , j ) + vy_inf * n4 (2 , i , j )
! Unit normal : only n4 (1 , i , nj ) and n4 (2 , i , nj ) are required
e_nx = n4(1 ,i ,j) / sqrt(n4(1 ,i ,j)**2 + n4(2 ,i ,j)**2)
e_ny = n4(2 ,i ,j) / sqrt(n4(1 ,i ,j)**2 + n4(2 ,i ,j)**2)
e_tx = e_ny ! tangential , 90 deg clockwise r o t a t i o n
e_ty = -e_nx ! tangential , 90 deg clockwise r o t a t i o n
! Calculate v a r i a b l e s at i n t e r i o r c e l l : (k , i , nj −1) , k = 1 ,2 ,3 ,4
p_i = pres_func(U(: ,i ,j-1))
rho_i = rho_func(U(: ,i ,j-1))
c_i = c_func(U(: ,i ,j-1))
u_n_i = U(2 ,i ,j-1) / U(1 ,i ,j-1) * e_nx + U(3 ,i ,j-1) / U(1 ,i ,j-1) * e_ny
! Inlet BC’ s
if (v_inf_dot_n .GT. 0.) then
rni = u_n_i - 2. * c_i / (gam - 1.) !Eq 4 , lec 14
rnb = vx_inf * e_nx + vy_inf * e_ny + 2. * c_inf / (gam - 1.) !Eq 5 , lec 14
! Calculate variables at boundary ( ’ ghost ’) cell : (k , i, nj ) , k = 1 ,2 ,3 ,4
u_n_b = 0.5 * (rnb + rni ) !Eq A, lec 14
c_b = (rnb - rni ) * (gam - 1.) * 0.25 !Eq B, lec 14
u_t_b = vx_inf * e_tx + vy_inf * e_ty !Eq C, lec 14
rho_b = (c_b**2 * rho_inf ** gam / p_inf / gam)**(1./(gam - 1.) ) ! Derived using S_B = S_infinity
p_b = c_b**2 * rho_b/gam
! Outlet BC’ s
else
rni = -abs(u_n_i) - 2. * c_i/(gam - 1.) !Eq B, lec . 14
rnb = -abs(vx_inf * e_nx + vy_inf * e_ny) + 2. * c_inf / (gam - 1.) !Eq A,lec . 14
! Calculate variables at boundary ( ’ ghost ’) cell : (k ,i ,nj ) , k = 1 ,2 ,3 ,4
u_n_b = 0.5 * (rnb + rni) !Eq A, lec 14
c_b = (rnb - rni) * (gam - 1.) * 0.25 !Eq B, lec 14
u_t_b = U(2 ,i ,j-1) / U(1 ,i ,j-1) * e_tx + U(3 ,i ,j-1) / U(1 ,i ,j-1) * e_ty !Eq C, lec 14
rho_b = (c_b**2 * rho_i**gam / p_i / gam)**(1. / (gam - 1.)) ! Derived using S_B = S_I
p_b = c_b**2 * rho_b / gam
endif
! Assign BC’ s
U(1 ,i ,j) = rho_b
U(2 ,i ,j) = rho_b * (u_n_b * e_ty - u_t_b * e_ny) / (e_nx * e_ty - e_ny * e_tx)
!Solve linear equations : vn = vx*nx +vy*ny ; vt = vx*tx +vy*ty
U(3 ,i ,j) = rho_b *(u_n_b * e_tx - u_t_b * e_nx) / (e_ny * e_tx - e_nx * e_ty)
!Solve linear equations : vn = vx*nx +vy*ny ; vt = vx*tx +vy*ty
E_b = (1. / (gam - 1.)) * (p_b / rho_b) + 0.5 * ((U(2 ,i ,j) / rho_b)**2 + (U(3 ,i ,j) / rho_b)**2)
U(4 ,i ,j) = rho_b * E_b
enddo
end subroutine inlet_outlet_bcs

