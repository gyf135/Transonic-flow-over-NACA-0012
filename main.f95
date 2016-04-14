program euler2D
implicit none
integer :: n, i , j , tmpi , tmpj
integer , parameter :: ni = 129-1 ! ni = 16 !128 !Number of grid points in x
integer , parameter :: nj = 65 ! nj = ni !Number of grid points in y
real , parameter :: gam = 1.4
real , parameter :: pi = 4. * atan (1.)
real , dimension ( ni , nj+1) :: x , y , x_cell , y_cell , sigma
real , dimension (2 , ni , nj ) :: n1 , n2 , n3 , n4 ! Cell normals
real , dimension ( ni , nj ) :: tau ! local timestep
real , dimension (4 , ni , nj ) :: U, F, G, Risd , U1, U2, U3, U_next , &
                                    Risd_1 , Risd_2 , Risd_3
real :: alt_inf , vx_inf , vy_inf , p_inf , T_inf , rho_inf , M_inf , c_inf , &
E_inf , H_inf , alpha ! I n i t i a l Condition Variables , S_inf
! Runge Kutta alpha values
real , parameter :: alpha1 = 1./8. , alpha2 = 0.306 , alpha3 = 0.587 , alpha4 = 1.
real , parameter :: dt = 0.0001
! Set Initial Conditions
alpha = 2. *( pi /180.)
! alpha = 0.
!write (* , *) pi , alpha
M_inf = 0.85
! Set US Standard Atmospheric Conditions
! http ://www. engineeringtoolbox . com/ standard−atmosphere−d_604 . html
alt_inf = 1000. ! 1000m above sea l e v e l
p_inf = 101325. * (1. - (2.25577e-5) * alt_inf)**5.25588 ! Pressure in Pa
T_inf = 8.50 + 273.15 ! Temperature in K
rho_inf = 11.12e-1 ! Density in kg/m^3
! V e l o c i t i e s are in m/ s
c_inf = (gam * p_inf / rho_inf ) ** 0.5
vx_inf = M_inf * c_inf * cos ( alpha )
vy_inf = M_inf * c_inf * sin ( alpha )
E_inf = (1. / ( gam - 1.)) * (p_inf / rho_inf) + 0.5 * (vx_inf**2 + vy_inf**2)
H_inf = E_inf + p_inf / rho_inf
U( 1 , : , : ) = rho_inf
U( 2 , : , : ) = rho_inf * vx_inf
U( 3 , : , : ) = rho_inf * vy_inf
U( 4 , : , : ) = rho_inf * E_inf
! Read in the Grid Data :
open ( unit = 2 , file = "coords.dat") ;
x=0.
y=0.
do i =1, ni
do j =1, nj
read ( unit =2, fmt=*) tmpi, tmpj, x(i ,j) ,y(i ,j)
enddo
enddo
! Calculate the c e l l volumes
call cell_volumes ( ni , nj , x , y , x_cell , y_cell , sigma , n1 , n2 , n3 , n4 )
call timestep ( ni , nj , n1 , n2 , n3 , n4 ,U, tau )
! Assign BC’ s at i n l e t and o u t l e t
open ( unit = 4 , file = " res2alpha.dat " ) ;
do n = 1 ,20000
call inlet_outlet_bcs ( ni , nj , n4 , vx_inf , vy_inf , p_inf , rho_inf , c_inf ,U) !Enforce BC’ s
call fluxes (ni ,nj ,U ,F ,G)
call timestep ( ni , nj , n1 , n2 , n3 , n4 ,U, tau )
call risdual ( ni , nj ,U, U, F,G, n1 , n2 , n3 , n4 , dt , sigma , alpha1 , Risd , U1, x_cell, y_cell, tau)
! Step 2 − Runge−Kutta
call inlet_outlet_bcs ( ni , nj , n4 , vx_inf , vy_inf , p_inf , rho_inf , c_inf , U1) !Enforce BC’ s
call fluxes ( ni , nj , U1 , F,G)
call risdual ( ni , nj ,U, U1 , F,G, n1 , n2 , n3 , n4 , dt , sigma , alpha2 , Risd_1 , U2, x_cell , y_cell , tau)
! Step 3 − Runge−Kutta
call inlet_outlet_bcs ( ni , nj , n4 , vx_inf , vy_inf , p_inf , rho_inf , c_inf , U2) !Enforce BC’ s
call fluxes ( ni , nj , U2 , F,G)
call risdual ( ni , nj ,U, U2 , F,G, n1 , n2 , n3 , n4 , dt , sigma , alpha3 , Risd_2 , U3 , x_cell , y_cell , tau)
! Step 4 − Runge−Kutta
call inlet_outlet_bcs ( ni , nj , n4 , vx_inf , vy_inf , p_inf , rho_inf , c_inf , U3) !Enforce BC’ s
call fluxes ( ni , nj , U3 , F, G)
call risdual ( ni , nj ,U, U3 , F,G, n1 , n2 , n3 , n4 , dt , sigma , alpha4 , Risd_3 , U_next ,x_cell , y_cell , tau)
! Output the v e l o c i t y magnitude and d i r e c t i o n
if (n .EQ. 20000) then
call output_velocities ( ni , nj , x_cell , y_cell ,U, F,G, n4 , p_inf , rho_inf , vx_inf , vy_inf)
end if
if ( maxval ( abs ( Risd_3 ( 1 , : , : ) ) ) .LE. 10.**(-5) .OR. &
maxval ( abs ( Risd_3 ( 2 , : , : ) ) ) .LE. 10.**(-5) .OR. &
maxval ( abs ( Risd_3 ( 3 , : , : ) ) ) .LE. 10.**(-5) .OR. &
maxval ( abs ( Risd_3 ( 4 , : , : ) ) ) .LE. 10.**(-5)) then
write ( *, *) n , maxval ( Risd_3 ( 1 , : , : ) ) , maxval ( Risd_3 ( 2 , : , : ) ) , maxval (Risd_3 ( 3 , : , : ) ) &
            , maxval ( Risd_3 ( 4 , : , : ) )
endif
U = U_next
write ( 4 , *) n , maxval ( Risd_3 ( 1 , : , : ) ) , maxval ( Risd_3 ( 2 , : , : ) ) , maxval (Risd_3 ( 3 , : , : ) ) &
             , maxval ( Risd_3 ( 4 , : , : ) )
enddo
end program euler2D
