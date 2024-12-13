module matrix_fun
	! Модуль с функциями, которые конструируют необходимые для задачи матрицы,
	! и производными от них вещами
use regime
implicit none
	

contains

! коэффициенты 3х-диагональной матрицы по Ходжу [BEGIN]
! https://doi.org/10.1063/1.1665398 (f.15-17)
function A_fun(m, r, c) result(f)
	! Вход:
	! 	m [integer]			: азимутальное число (порядок), m >= 0
	! 	r [integer]			: индекс коэффициента разложения угловой сфероидальной функции / элемента матрицы, r >= 0 
	! 	c [complex(knd)]	: c = kd/2, k - волновое число в среде, d - интерфокальное расстояние

	! Выход:
	! 	f [complex(knd)]	: значение коэффициента A^m_r(c)

	integer, intent(in)			:: m, r
	complex(knd), intent(in)	:: c

	complex(knd)				:: f

	f = (2._knd*m+r+2) * (2*m+r+1) / (2*m+2*r+3) / (2*m+2*r+5) * c**2
	
end function A_fun

function B_fun(m, r, c) result(f)
	! Вход:
	! 	m [integer]			: азимутальное число (порядок), m >= 0
	! 	r [integer]			: индекс коэффициента разложения угловой сфероидальной функции / элемента матрицы, r >= 0 
	! 	c [complex(knd)]	: c = kd/2, k - волновое число в среде, d - интерфокальное расстояние

	! Выход:
	! 	f [complex(knd)]	: значение коэффициента B^m_r(c)

	integer, intent(in)			:: m, r
	complex(knd), intent(in)	:: c
	complex(knd)				:: f

	f = (2._knd*(m+r)*(m+r+1) - 2*m**2 - 1) / (2*m+2*r-1) / (2*m+2*r+3) * c**2 + (m+r) * (m+r+1)
	
end function B_fun

function C_fun(m, r, c) result(f)
	! Вход:
	! 	m [integer]			: азимутальное число (порядок), m >= 0
	! 	r [integer]			: индекс коэффициента разложения угловой сфероидальной функции / элемента матрицы, r >= 0 
	! 	c [complex(knd)]	: c = kd/2, k - волновое число в среде, d - интерфокальное расстояние

	! Выход:
	! 	f [complex(knd)]	: значение коэффициента C^m_r(c)

	integer, intent(in)			:: m, r
	complex(knd), intent(in)	:: c
	complex(knd)				:: f

	f = r * (r-1._knd) / (2*m+2*r-3) / (2*m+2*r-1) * c**2
	
end function C_fun
! коэффициенты 3х-диагональной матрицы по Ходжу [END]


! составляющие коэффициентов преобразованный 3х-диагональной матрицы по Ходжу [BEGIN]
! https://doi.org/10.1063/1.1665398 (f.18-20)
function D_q(q, m, n, c) result(f)
	! Вход:
	! 	q [integer]			: индекс коэффициента разложения угловой сфероидальной функции / элемента матрицы, q >= 0 
	! 	m [integer]			: азимутальное число (порядок), m >= 0
	! 	n [integer]			: ??? число (степень), n >= m
	! 	c [complex(knd)]	: c = kd/2, k - волновое число в среде, d - интерфокальное расстояние

	! Выход:
	! 	f [complex(knd)]	: значение коэффициента D^m_q(c) = C^m_2q+s(c)

	integer, intent(in)			:: m, n, q
	complex(knd), intent(in)	:: c

	complex(knd)				:: f

	f = C_fun(m, 2*q+mod(n-m,2), c)

end function D_q

function E_q(q, m, n, c) result(f)
	! Вход:
	! 	q [integer]			: индекс коэффициента разложения угловой сфероидальной функции / элемента матрицы, q >= 0 
	! 	m [integer]			: азимутальное число (порядок), m >= 0
	! 	n [integer]			: ??? число (степень), n >= m
	! 	c [complex(knd)]	: c = kd/2, k - волновое число в среде, d - интерфокальное расстояние

	! Выход:
	! 	f [complex(knd)]	: значение коэффициента E^m_q(c) = B^m_2q+s(c)

	integer, intent(in)			:: m, n, q
	complex(knd), intent(in)	:: c

	complex(knd)				:: f

	f = B_fun(m, 2*q+mod(n-m,2), c)

end function E_q

function F_q(q, m, n, c) result(f)
	! Вход:
	! 	q [integer]			: индекс коэффициента разложения угловой сфероидальной функции / элемента матрицы, q >= 0 
	! 	m [integer]			: азимутальное число (порядок), m >= 0
	! 	n [integer]			: ??? число (степень), n >= m
	! 	c [complex(knd)]	: c = kd/2, k - волновое число в среде, d - интерфокальное расстояние

	! Выход:
	! 	f [complex(knd)]	: значение коэффициента F^m_q(c) = A^m_2q+s(c)

	integer, intent(in)			:: m, n, q
	complex(knd), intent(in)	:: c

	complex(knd)				:: f

	f = A_fun(m, 2*q+mod(n-m,2), c)

end function F_q
! составляющие коэффициентов преобразованный 3х-диагональной матрицы по Ходжу [END]


! изначальная 3х-диагональная матрица в полном виде по Ходжу
function ABC_matr_old(NN, m, n, c) result(f)
	integer, intent(in) :: m, n, NN
	integer :: r, s, k
	complex(knd), intent(in) :: c
	complex(knd) :: f(0:NN-1,0:NN-1)

	f = 0

	s = mod(n-m, 2)

	k = 0

	do r=s, 2*(NN-1)-1, 2
		k = k + 1

		f(k,k-1) = C_fun(m, r, c)
		f(k,k) = B_fun(m, r, c)
		f(k,k+1) = A_fun(m, r, c)
	end do

	

end function ABC_matr_old


! изначальная 3х-диагональная матрица в полном виде по Ходжу
function ABC_matr(NN, m, n, c) result(f)
	integer, intent(in) :: m, n, NN
	integer :: r, s, k
	complex(knd), intent(in) :: c
	complex(knd) :: f(0:NN-1,0:NN-1)

	f = 0

	do r=0, (NN-1)-2
		f(r,r+2) = A_fun(m, r, c)
		f(r,r) = B_fun(m, r, c)
		f(r+2,r) = C_fun(m, r+2, c)
	end do

	f((NN-1)-1,(NN-1)-1) = B_fun(m, (NN-1)-1, c)
	f((NN-1),(NN-1)) = B_fun(m, (NN-1), c)

end function ABC_matr

! симметричная (преобразованная) 3х-диагональная матрица в полном виде по Ходжу
function DFE_matr(NN, m, n, c) result(f)
	integer, intent(in) :: m, n, NN
	integer :: r, s
	complex(knd), intent(in) :: c
	complex(knd) :: f(0:NN-1,0:NN-1)

	f = 0

	do r=0, (NN-1)-1
		f(r,r+1) = sqrt(D_q(r+1, m, n, c) * F_q(r, m, n, c))
		f(r,r) = E_q(r, m, n, c)
		f(r+1,r) = f(r,r+1)
	end do
	r = NN-1
	f(r,r) = E_q(r, m, n, c)


end function DFE_matr



! симметричная (преобразованная) 3х-диагональная матрица в полном виде по Ходжу
function DFE_matr_old(NN, m, n, c) result(f)
	integer, intent(in) :: m, n, NN
	integer :: r, s
	complex(knd), intent(in) :: c
	complex(knd) :: f(NN,NN)

	f = 0
	s = mod(n-m, 2)

	r = 1
	! f(r,r-1) = sqrt(D_q(r, m, n, c) * F_q(r-1, m, n, c))
	f(r,r) = E_q(r, m, n, c)
	f(r,r+1) = sqrt(D_q(r+1, m, n, c) * F_q(r, m, n, c))

	r = NN
	f(r,r-1) = sqrt(D_q(r, m, n, c) * F_q(r-1, m, n, c))
	f(r,r) = E_q(r, m, n, c)
	! f(r,r+1) = sqrt(D_q(r+1, m, n, c) * F_q(r, m, n, c))

	do r=2, NN-1
		f(r,r-1) = sqrt(D_q(r, m, n, c) * F_q(r-1, m, n, c))
		f(r,r) = E_q(r, m, n, c)
		f(r,r+1) = sqrt(D_q(r+1, m, n, c) * F_q(r, m, n, c))
	end do

end function DFE_matr_old

! симметричная (преобразованная) 3х-диагональная матрица в компактном виде по Ходжу
function DFE_matr_mini(NN, m, n, c) result(f)
	integer, intent(in)			:: m, n, NN
	complex(knd), intent(in)	:: c

	integer						:: r, s

	complex(knd)				:: f(NN,2)

	f = 0
	s = mod(n-m, 2)

	do r=0, NN-1
		f(r+1,1) = E_q(r, m, n, c)
		f(r+1,2) = sqrt(D_q(r+1, m, n, c) * F_q(r, m, n, c))
	end do

	! IF (ieee_support_inf(r)) print*, .true.
	! f(N,2) = ieee_value(f(N,2),  ieee_negative_inf)
	print*, size(f), NN
	f(NN,2) = complex(-9999, -9999)

end function DFE_matr_mini


function DFE_matr_mini_full(NN, m, n, c) result(f)
	integer, intent(in)			:: m, n, NN
	integer						:: r, s
	complex(knd), intent(in)	:: c
	complex(knd)				:: f(NN/2,4)

	f = 0

	f(:,:2) = DFE_matr_mini(NN/2, m, n, c)
	f(:,3:) = DFE_matr_mini(NN/2, n+1, m, c)

	! print*, "#", f(:,:2)
	! print*, "#", f(:,3:)

	
end function DFE_matr_mini_full



! ----------------------------------------------------------------------------------------------------------------------------------
function from_sim3diagMini_to_sim3diagBase(mtr, N) result(retval)
	integer :: i, N
	complex(knd), intent(in) :: mtr(N,2)
	complex(knd) :: retval(N,N)

	retval = 0

	do i=1, N-1
		retval(i,i) = mtr(i,1)
		retval(i,i+1) = mtr(i,2)
		retval(i+1,i) = mtr(i,2)
	end do

	i = N
	retval(i,i) = mtr(i,1)
	! retval(i,i+1) = mtr(i+1,2)
	! retval(i+1,i) = mtr(i+1,2)
end function from_sim3diagMini_to_sim3diagBase

function from_sim3diagMini_to_sim3diagBase_8(mtr, N) result(retval)
	integer :: i, N
	complex(8), intent(in) :: mtr(N,2)
	complex(8) :: retval(N,N)

	retval = 0

	do i=1, N-1
		retval(i,i) = mtr(i,1)
		retval(i,i+1) = mtr(i,2)
		retval(i+1,i) = mtr(i,2)
	end do

	i = N
	retval(i,i) = mtr(i,1)
	! retval(i,i+1) = mtr(i+1,2)
	! retval(i+1,i) = mtr(i+1,2)
end function from_sim3diagMini_to_sim3diagBase_8



! function coef_d(m, n, c) result(retval)
! 	integer, intent(in) :: m, n
! 	complex(knd), intent(in) :: c

! 	integer :: NN
! 	complex(knd), allocatable :: mtr(:,:), lambdas(:)

! 	integer :: retval


! 	NN = 5

! 	allocate(mtr(NN,2))
! 	mtr = DFE_matr_mini(NN, m, n, c)

! 	allocate(lambdas(NN))
! 	lambdas = eigenvalues_of_matrix(mtr%Re, NN)

! 	print*, lambdas 


	
! end function coef_d


function from_sim3diagBase_to_sim3diagMini(mtr, N) result(mtr_mini)
	integer, intent(in)			:: N
	complex(knd), intent(in)	:: mtr(N, N)

	integer						:: i, j
	complex(knd)				:: mtr_mini(N, 2)
	! integer :: retval

	mtr_mini = 0


	do i=1, N-1
		mtr_mini(i,1) = mtr(i,i)
		mtr_mini(i,2) = mtr(i,i+1)
	end do

	mtr_mini(N,1) = mtr(N,N)
	mtr_mini(N,2) = complex(-9999, -9999)
	
end function from_sim3diagBase_to_sim3diagMini

subroutine print_matr_mini(mtr_mini, N) 
	integer, intent(in)			:: N
	complex(knd), intent(in)	:: mtr_mini(N, 2)

	integer						:: i, j
	complex(knd)				:: mtr(N, N)
	! integer :: retval

	mtr = from_sim3diagMini_to_sim3diagBase(mtr_mini, N)

	do i=1, N
		write(*,*) (mtr(i,j), j=1, N)
	end do

end subroutine print_matr_mini


subroutine print_matr(mtr, N) 
	integer, intent(in)			:: N
	complex(knd), intent(in)	:: mtr(N, N)

	integer						:: i, j

	do i=1, N
		write(*,*) (mtr(i,j), j=1, N)
	end do

end subroutine print_matr




! коэффициенты 3х-диагональной матрицы по Баукампу (begin)
function A_fun_b(m, r, c) result(f)
	integer, intent(in)			:: m, r
	complex(knd), intent(in)	:: c
	complex(knd)				:: f

	! f = (r + 1._knd) * (r + 2) / (2*r + 1) / (2*r + 3) * c**2
	f = C_fun(m, r+2, c)
end function A_fun_b

function B_fun_b(m, r, c) result(f)
	integer, intent(in)			:: m, r
	complex(knd), intent(in) 	:: c
	complex(knd)				:: f

	f = (2._knd*(m+r)*(m+r+1) - 2*m**2 - 1) / (2*m+2*r-1) / (2*m+2*r+3) * c**2 - (m+r) * (m+r+1)
	! f  = (2*r*(r + 1._knd) - 1) / (2*r - 1) / (2*r + 3) * c**2 - r*(r+1)
end function B_fun_b

function C_fun_b(m, r, c) result(f)
	integer, intent(in)			:: m, r
	complex(knd), intent(in)	:: c
	complex(knd)				:: f

	! f = r * (r - 1._knd) / (2*r - 1) / (2*r + 1) * c**2
	f = A_fun(m, r-2, c)
end function C_fun_b
! коэффициенты 3х-диагональной матрицы по Баукампу (end)

function ABC_matr_b(NN, m, n, c) result(f)
	integer, intent(in) :: m, n, NN
	integer :: r, s, k
	complex(knd), intent(in) :: c
	complex(knd) :: f(0:NN-1,0:NN-1)

	f = 0

	do r=0, (NN-1)-2
		f(r,r+2) = C_fun_b(m, r+2, c)
		f(r,r) = B_fun_b(m, r, c)
		f(r+2,r) = A_fun_b(m, r, c)
	end do

	f((NN-1)-1,(NN-1)-1) = B_fun_b(m, (NN-1)-1, c)
	f((NN-1),(NN-1)) = B_fun_b(m, (NN-1), c)

end function ABC_matr_b

function DFE_matr_b(NN, m, n, c) result(f)
	integer, intent(in)			:: m, n, NN
	integer						:: r, s, k
	complex(knd), intent(in)	:: c
	complex(knd)				:: f(0:NN-1,0:NN-1)

	f = 0

	do r=0, (NN-1)-1
		f(r,r+1) = sqrt(D_q_b(r+1, m, n, c) * F_q_b(r, m, n, c))
		f(r,r) = E_q_b(r, m, n, c)
		f(r+1,r) = f(r,r+1)
	end do

	r = NN-1
	f(r,r) = E_q_b(r, m, n, c)

end function DFE_matr_b

function DFE_matr_mini_b(NN, m, n, c) result(f)
	integer, intent(in)			:: m, n, NN
	integer						:: r, s
	complex(knd), intent(in)	:: c
	complex(knd)				:: f(NN,2)

	f = 0
	s = mod(n-m, 2)

	do r=0, NN-1
		f(r+1,1) = E_q_b(r, m, n, c)
		f(r+1,2) = sqrt(D_q_b(r+1, m, n, c) * F_q_b(r, m, n, c))
	end do

	! IF (ieee_support_inf(r)) print*, .true.
	! f(N,2) = ieee_value(f(N,2),  ieee_negative_inf)
	f(NN,2) = complex(-9999, -9999)

end function DFE_matr_mini_b

function DFE_matr_mini_b_full(NN, m, n, c) result(f)
	integer, intent(in)			:: m, n, NN
	integer						:: r, s
	complex(knd), intent(in)	:: c
	complex(knd)				:: f(NN/2,4)

	f = 0

	f(:,:2) = DFE_matr_mini_b(NN/2, m, n, c)
	f(:,3:) = DFE_matr_mini_b(NN/2, n+1, m, c)

	! print*, "#", f(:,:2)
	! print*, "#", f(:,3:)

	
end function DFE_matr_mini_b_full

function D_q_b(q, m, n, c) result(f)
	integer, intent(in) :: m, n, q
	complex(knd), intent(in) :: c
	complex(knd) :: f

	f = C_fun_b(m, 2*q+mod(n-m,2), c)

end function D_q_b

function E_q_b(q, m, n, c) result(f)
	integer, intent(in) :: m, n, q
	complex(knd), intent(in) :: c
	complex(knd) :: f

	f = B_fun_b(m, 2*q+mod(n-m,2), c)

end function E_q_b

function F_q_b(q, m, n, c) result(f)
	integer, intent(in) :: m, n, q
	complex(knd), intent(in) :: c
	complex(knd) :: f

	f = A_fun_b(m, 2*q+mod(n-m,2), c)

end function F_q_b

end module matrix_fun