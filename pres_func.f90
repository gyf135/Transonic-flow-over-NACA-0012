! Computes pressure at single grid location using components of U
real function pres_func(U_vec)
implicit none
real, parameter :: gam = 1.4
real, dimension(4) :: U_vec
pres_func = (gam - 1.) * (U_vec(4) &
          - 1. / ( 2. * U_vec(1) ) * (U_vec(2)**2 + U_vec(3)**2))
return
end function pres_func

! Computes density at single grid location using components of U
real function rho_func(U_vec)
implicit none
real, dimension(4) :: U_vec
rho_func = U_vec (1)
return
end function rho_func

! Computes speed of sound at s i n g l e grid l o c a t i o n using components of U
real function c_func (U_vec)
implicit none
real, parameter :: gam = 1.4
real, dimension(4) :: U_vec
real :: pres_func, rho_func
c_func = (gam * pres_func (U_vec) / rho_func(U_vec))**0.5
return
end function c_func

real function vx_func (U_vec)
implicit none
real, dimension(4) :: U_vec
vx_func = U_vec(2) / U_vec(1)
return
end function vx_func

real function vy_func(U_vec)
implicit none
real, dimension (4) :: U_vec
vy_func = U_vec(3) / U_vec(1)
return
end function vy_func

real function E_func(U_vec)
implicit none
real, parameter :: gam = 1.4
real, dimension(4) :: U_vec
real :: pres_func , rho_func , vx_func , vy_func
real :: e
e = pres_func(U_vec) / (rho_func(U_vec) * (gam - 1.))
E_func = e + 0.5 * (vx_func(U_vec)**2 + vy_func(U_vec)**2)
return
end function E_func

real function H_func(U_vec)
implicit none
real, dimension(4) :: U_vec
real :: pres_func ,rho_func ,E_func
H_func = E_func(U_vec) + pres_func(U_vec) / rho_func(U_vec)
return
end function H_func

real function M_func(U_vec)
implicit none
real, dimension(4) :: U_vec
real :: vx_func ,vy_func ,c_func
M_func = (vx_func(U_vec)**2+vy_func(U_vec)**2)**0.5 / c_func(U_vec)
return
end function M_func

real function S_func(U_vec)
implicit none
real, parameter :: c_v = 0.718
real, parameter :: gam = 1.4
real, dimension(4) :: U_vec
real :: pres_func ,rho_func
S_func = c_v * log(pres_func(U_vec) / (rho_func(U_vec))**gam)
return
end function S_func

real function cp_func(U_vec ,p_inf ,rho_inf ,vx_inf ,vy_inf )
implicit none
real, dimension(4) :: U_vec
real :: pres_func
real :: p_inf ,rho_inf ,vx_inf ,vy_inf ,vmag_inf
vmag_inf = (vx_inf**2 + vy_inf**2)**0.5
cp_func = ( pres_func(U_vec) - p_inf) / (0.5 * rho_inf * vmag_inf**2)
return
end function cp_func

