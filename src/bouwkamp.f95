module bouwkamp
use regime
use matrix_fun
implicit none

abstract interface
	function fun_B(m, r, c) result(f)
		import knd
		integer, intent(in) :: m, r
		complex(knd), intent(in) :: c
		complex(knd) :: f
	end function fun_B
end interface

! procedure(fun_B), pointer :: alpha_fun => B_fun
! procedure(fun_B), pointer :: alpha_fun => B_fun2
! procedure(fun_B), pointer :: alpha_fun => alpha_fun_bouwkamp
! procedure(fun_B), pointer :: alpha_fun => alpha_fun_bouwkamp2
procedure(fun_B), pointer :: alpha_fun => gamma_fun
procedure(fun_B), pointer :: beta_fun => beta_fun1
! procedure(fun_B), pointer :: beta_fun => beta_fun2
	
contains

function beta_fun1(m, r, c) result(f)
	integer, intent(in) :: m, r
	complex(knd), intent(in) :: c
	complex(knd) :: f

	f = C_fun_b(m, r, c) * A_fun_b(m, r-2, c)

end function beta_fun1

function beta_fun2(m, r, c) result(f)
	integer, intent(in) :: m, r
	complex(knd), intent(in) :: c
	complex(knd) :: f

	f = C_fun(m, r, c) * A_fun(m, r-2, c)


	! f = r * (r-1._knd) * (2*m+r) * (2*m+r-1) / (2*m+2*r-1)**2 / (2*m+2*r-3) / (2*m+2*r+1) * c**4

end function beta_fun2

function B_fun2(m, r, c) result(f)
	integer, intent(in) :: m, r
	complex(knd), intent(in) :: c
	complex(knd) :: f

	f = -(2._knd*(m+r)*(m+r+1) - 2*m**2 - 1) / (2*m+2*r-1) / (2*m+2*r+3) * c**2 + (m+r) * (m+r+1)
	
end function B_fun2

function alpha_fun_bouwkamp2(m, r, c) result(f)
	integer, intent(in) :: m, r
	complex(knd), intent(in) :: c
	complex(knd) :: f

	! f = (m+r) * (m+r+1) - (1 + 1._knd/(2*m+2*r-1)/(2*m+2*r+3)) * c**2 / 2
	f = -B_fun_b(m, r, c)

end function alpha_fun_bouwkamp2

function alpha_fun_bouwkamp(m, r, c) result(f)
	integer, intent(in) :: m, r
	complex(knd), intent(in) :: c
	complex(knd) :: f

	f = r*(r+1) - (1 + 1._knd/(2*r-1)/(2*r+3)) * c**2 / 2

end function alpha_fun_bouwkamp


function gamma_fun(m, r, c) result(f)
	integer, intent(in) :: m, r
	complex(knd), intent(in) :: c
	complex(knd) :: f

	f = (m+r)*(m+r+1._knd) + (1._knd - (4*m**2-1._knd)/(2*m+2*r-1)/(2*m+2*r+3)) * c**2 / 2
	! f = B_fun(m, r, c)

end function gamma_fun

! N_n = f(N_n+2)
! recursive function N_fun_down(m, r, c, lambda, i) result(f)
! 	integer, intent(in) :: m, r, i
! 	complex(knd), intent(in) :: c, lambda
! 	complex(knd) :: f

! 	if (i >= 20) then
! 		f = beta_fun(m, r, c) / (alpha_fun(m, r, c) - lambda)
! 	else
! 		f = beta_fun(m, r, c) / (alpha_fun(m, r, c) - lambda - N_fun_down(m, r+2, c, lambda, i+1))
! 	end if
	
! end function N_fun_down
recursive subroutine N_fun_down(m, r, c, lambda, N_arr, beta_arr, i, iter, f) 
	integer, intent(in) :: m, r, i, iter
	complex(knd), intent(in) :: c, lambda
	complex(knd), intent(inout) :: N_arr(:), beta_arr(:)
	complex(knd), intent(out) :: f
	! integer, intent(out) :: i

	complex(knd) :: beta_n

	if (i >= iter) then
		beta_n = beta_fun(m, r, c)
		f = 0

		N_arr(i+1) = f
		beta_arr(i+1) = beta_n
	else
		call N_fun_down(m, r+2, c, lambda, N_arr, beta_arr, i+1, iter, f)
		beta_n = beta_fun(m, r, c)

		f = beta_n / (alpha_fun(m, r, c) - lambda - f)

		N_arr(i+1) = f
		beta_arr(i+1) = beta_n
	end if

	! print*, i, f
	
end subroutine N_fun_down

! N_n+2 = f(N_n)
! recursive function N_fun_up(m, r, c, lambda) result(f)
! 	integer, intent(in) :: m, r
! 	complex(knd), intent(in) :: c, lambda
! 	complex(knd):: f

! 	if (r == 2) then
! 		f = alpha_fun(m, 0, c) - lambda
! 	elseif (r == 3) then
! 		f = alpha_fun(m, 1, c) - lambda
! 	else
! 		f = alpha_fun(m, r-2, c) - lambda - (beta_fun(m, r-2, c) / N_fun_up(m, r-2, c, lambda))
! 	end if
	
! end function N_fun_up
recursive subroutine N_fun_up(m, r, c, lambda, N_arr, beta_arr, i, f) 
	integer, intent(in) :: m, r
	complex(knd), intent(in) :: c, lambda
	complex(knd), intent(inout) :: N_arr(:), beta_arr(:)
	complex(knd), intent(out) :: f
	integer, intent(out) :: i

	complex(knd) :: beta_n

	if (r == 2) then
		i = 1

		f = alpha_fun(m, 0, c) - lambda
		N_arr(i) = f
		! print*, 2, i, f, alpha_fun(m, 0, c)
		! print*, f, i, r
	elseif (r == 3) then
		i = 1

		f = alpha_fun(m, 1, c) - lambda
		N_arr(i) = f
		! print*, 3, i, f, alpha_fun(m, 1, c)
		! print*, f, i, r


	else
		call N_fun_up(m, r-2, c, lambda, N_arr, beta_arr, i, f)
		i = i + 1

		beta_n = beta_fun(m, r-2, c)

		f = alpha_fun(m, r-2, c) - lambda - (beta_n / f)

		N_arr(i) = f
		beta_arr(i-1) = beta_n
		! print*, 0, i, f
		! print*, f, i, r


	end if
	
end subroutine N_fun_up

! function U_fun(m, c, lambda, j) result(f)
subroutine U_fun(m, c, lambda, r, N_up_arr, N_down_arr, beta_up_arr, beta_down_arr, iter, f)
	integer, intent(in) :: m, r, iter
	complex(knd), intent(in) :: c, lambda
	complex(knd), intent(inout) :: N_up_arr(:), N_down_arr(:), beta_up_arr(:), beta_down_arr(:)
	complex(knd), intent(out) :: f

	integer :: i
	complex(knd) :: U_1, U_2
	! complex(knd), allocatable :: N_up_arr(:), N_down_arr(:), beta_up_arr(:), beta_down_arr(:)
	

	! r = j + 2

	! iter = 7

	! allocate(N_up_arr(r/2), N_down_arr(iter+1), beta_up_arr(r/2), beta_down_arr(iter+1))
! print*, "lol"

	! U_1(lambda_j) = N_fun_up(m, r=j+2, c, lambda_j)
	call N_fun_up(m, r, c, lambda, N_up_arr, beta_up_arr, i, U_1) 
! print*, "N_up_arr", N_up_arr
! print*, "beta_up_arr", beta_up_arr

	! U_2(lambda_j) = -N_fun_down(m, r=j+2, c, lambda_j)
	call N_fun_down(m, r, c, lambda, N_down_arr, beta_down_arr, 0, iter, U_2) 
	U_2 = -U_2
! print*, "N_down_arr", N_down_arr
! print*, "beta_down_arr", beta_down_arr

! print*, "lol3"

	! beta_up_arr(r/2) = beta_down_arr(1)
! print*, "lol4"

	! ! U_2(lambda_j) = -N_fun_down(m, r=j+2, c, lambda, i)
	! U_2 = -N_fun_down(m, r, c, lambda, 0)
	! ! U_1(lambda_j) = N_fun_up(m, r=j+2, c, lambda)
	! U_1 = N_fun_up(m, r, c, lambda)
	! print*, N_up_arr
	! print*, "---------------------------"
	! print*, N_down_arr
	! print*, "---------------------------"
	! print*, N_up_arr(size(N_up_arr)) - N_down_arr(1)
	! print*, "---------------------------"
	! print*, beta_up_arr
	! print*, "---------------------------"
	! print*,  beta_down_arr

	! print*, "u1", U_1
	! print*, "u2", U_2
	f = U_1 + U_2
	
end subroutine U_fun
! end function U_fun
	

function dU_fun(m, c, lambda, j, N_up_arr, N_down_arr, beta_up_arr, beta_down_arr) result(f)
	integer, intent(in) :: m, j
	complex(knd), intent(in) :: c, lambda, N_up_arr(:), N_down_arr(:), beta_up_arr(:), beta_down_arr(:)

	integer :: i, len_N
	complex(knd) :: dU_1, dU_2, f

	len_N = size(N_up_arr)-1

	! dU_1
	dU_1 = 1
	do i = 1, size(beta_up_arr)
		! dU_1 = dU_1 + product(beta_up_arr(i:)) / product(N_up_arr(i:len_N)**2)
		dU_1 = dU_1 + product(beta_up_arr(i:) / N_up_arr(i:len_N)**2)
		! print*,product(beta_up_arr(i:)),  product(N_up_arr(i:len_N)**2), product(beta_up_arr(i:) / N_up_arr(i:len_N)**2), dU_1
	end do

	! dU_2
	dU_2 = 0
	do i = 1, size(N_down_arr)
		! dU_2 = dU_2 + product(N_down_arr(:i)**2) / product(beta_down_arr(:i))
		dU_2 = dU_2 + product(N_down_arr(:i)**2 / beta_down_arr(:i))
		! print*, "####", i, product(N_down_arr(:i)**2), product(beta_down_arr(:i)), &
		!  product(N_down_arr(:i)**2) / product(beta_down_arr(:i))
	end do

	f = dU_1 + dU_2

	! print*, "dU_1", dU_1
	! print*, "dU_2", dU_2

end function dU_fun


function delta_lambda(m, c, lambda, j) result(f)
	integer, intent(in) :: m, j
	complex(knd), intent(in) :: c, lambda

	integer :: r, iter, i
	! complex(knd) :: U_1, U_2
	complex(knd), allocatable :: N_up_arr(:), N_down_arr(:), beta_up_arr(:), beta_down_arr(:)
	complex(knd) :: f, U, dU
	

	r = j + 2

	iter = 50 !33 - странное,ь 34-не работает

	allocate(N_up_arr(r/2), N_down_arr(iter+1), beta_up_arr(r/2-1), beta_down_arr(iter+1))

	call U_fun(m, c, lambda, r, N_up_arr, N_down_arr, beta_up_arr, beta_down_arr, iter, U)	

	! print*, "-----------"

	dU = dU_fun(m, c, lambda, j, N_up_arr, N_down_arr, beta_up_arr, beta_down_arr) 

	deallocate(N_up_arr, N_down_arr, beta_up_arr, beta_down_arr)

	f = U / dU

	! print*, " U", U
	! print*, "dU", dU

	
end function delta_lambda


end module bouwkamp