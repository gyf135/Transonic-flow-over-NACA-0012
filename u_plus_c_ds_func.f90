real function u_plus_c_ds_func(ni ,nj ,U ,n_out ,i ,i2 ,j ,j2)
implicit none
integer :: ni ,nj ,i ,i2 ,j ,j2
real, dimension (4 ,ni ,nj ) :: U
real, dimension (2 ,ni ,nj ) :: n_out ! Outward normal
real :: vx_func ,vy_func ,c_func
real :: vx_face ,vy_face ,c_face ,ds_mag
vx_face = 0.5 * (vx_func(U(: ,i ,j)) + vx_func(U(: ,i2 ,j2)))
vy_face = 0.5 * (vy_func(U(: ,i ,j)) + vy_func(U(: ,i2 ,j2)))
c_face = 0.5 * ( c_func(U(: ,i ,j)) + c_func(U(: ,i2 ,j2)))
ds_mag = (n_out(1 ,i ,j)**2 + n_out(2 ,i ,j)**2)**0.5
u_plus_c_ds_func = vx_face * n_out(1 ,i ,j ) + vy_face * n_out(2 ,i ,j) &
                 + c_face * ds_mag
return
end function u_plus_c_ds_func

