module legendre
use regime
implicit none
	
contains
	
! recursive function Associated_Legendre_polynomials(n, m, x) result(P)
! 	integer, intent(in)		:: n, m
! 	real(knd), intent(in)	:: x

! 	integer					:: i

! 	real(knd)				:: P


! 	if (n == m) then
! 		P = (-1)**m * product([(i,i=2*m-1,1,-2)]) * (1-x**2)**(m/2._knd)
! 	elseif (n == m+1) then
! 		P = x * (2*m+1) * Associated_Legendre_polynomials(m, m, x)
! 	else
! 		P = ((2*n-1) * x * Associated_Legendre_polynomials(n-1, m, x) - (n-1+m) * Associated_Legendre_polynomials(n-2, m, x)) / (n-m)
! 	end if
	
! end function Associated_Legendre_polynomials



end module legendre