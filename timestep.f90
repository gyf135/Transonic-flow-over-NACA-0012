subroutine timestep(ni ,nj ,n1 ,n2 ,n3 ,n4 ,U ,tau )
integer :: ni ,nj ,i ,j
real :: cfl
real, dimension (ni ,nj) :: tau
real, dimension (2 ,ni ,nj ) :: n1 , n2 , n3 , n4 , dsi , dsj ! Cell normals
real, dimension (4 ,ni ,nj) :: U
real :: clocal , ulocal , vlocal , plocal , rholocal
real, parameter :: gam = 1.4
cfl = 1.0 !Max value
dsi = 0.5 * (n1 - n3)
dsj = 0.5 * (n2 - n4)
do i =1 ,ni
do j =1 ,nj - 1
plocal = pres_func (U(: ,i ,j))
rholocal = U(1 ,i ,j)
clocal = sqrt(abs(gam * plocal / rholocal))
ulocal = U(2 ,i ,j) / U(1 ,i ,j)
vlocal = U(3 ,i ,j) / U(1 ,i ,j)
tau(i ,j) = cfl / (abs((ulocal + clocal) * dsi(1 ,i ,j) + (vlocal + clocal) * dsi(2 ,i ,j)) &
          + abs((ulocal + clocal) * dsj(1 ,i ,j) + (vlocal + clocal) * dsj(2 ,i ,j)))
enddo
enddo
end subroutine timestep

