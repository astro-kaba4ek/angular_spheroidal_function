module ang_func
use regime
use mod_eigenvalues
use matrix_fun
implicit none

contains


! ФУНКЦИИ 1-ГО ПОРЯДКА -------------------------------------------- [BEGIN]
! ВЫЧИСЛЕНИЕ НЕПОСРЕДСТВЕННЫХ СОСТАВЛЯЮЩИХ ФОРМУЛЫ 
! S_mn^norm(с,eta) = SUM_r[d_r(c) * P_m+r^m(eta)] / N_mn(c)
! (ф.10 -- Hodge,1970; ф.A.9 -- Il'in+,2023)
function N_mn(m, n, c, p, NN, d_coef) result(Nmn)
	! Фактор нормализации. 
	! Для различных нормировок S_mn(c,eta) 

	! Вход:
	! 	m [integer]							: азимутальное число (порядок), m >= 0
	! 	n [integer]							: индекс (степень), n >= m
	! 	c [complex(knd)]					: c = kd/2, k -- волновое число в среде, d -- интерфокальное расстояние
	!	p [integer]							: некоторое большое число, такое что d_p(c) почти 0, p > n-m
	! 	NN [integer]						: размер матрицы рекурентного расчета d_r(c), количество собственных значений
	! 	d_coef((p+2)/2+1) [complex(knd)]	: значения коэффициентов d_r(c) 
	
	! Выход:
	! 	Nmn [complex(knd)]					: фактор нормализации
	
	integer, intent(in) 		:: n, m, p, NN
	complex(knd), intent(in)	:: c
	complex(knd), intent(in)	:: d_coef(:)

	integer						:: i, j , r, s
	complex(knd)				:: sum1

	complex(knd)				:: Nmn


	s = mod(n-m, 2)

	if (c == 0) then
		
		Nmn = 1q0
		do j=(n-m)+1, n+m
			Nmn = Nmn * j
		end do
		Nmn = sqrt(2._knd * Nmn / (2*n+1))

	else
		Nmn = 0
		do i=1, size(d_coef)
			r = 2*(i-1)+s

			sum1 = 1q0
			do j=r+1, r+2*m
				sum1 = sum1 * j
			end do
			sum1 = sum1 * 2 / (2*m + 2*r + 1) * d_coef(i)**2

			Nmn = Nmn + sum1
		
		end do
		Nmn = sqrt(Nmn)

	end if 

end function N_mn


function S_mn_eta_n(c, eta, m, n_max, norm) result(Smn)
	! Функция для расчета угловых сфероидальных функций для набора 
	! значений n = [m, n_max] x eta = [eta_1, eta_k]
	! при заданных m, c и norm

	! Вход:
	! 	c [complex(knd)]						: c = kd/2, k -- волновое число в среде, d -- интерфокальное расстояние
	!	eta(eta_1:eta_k) [real(knd)]			: сфероидальная координата (массив) \eta (поверхность -- гиперболоид вращения), [-1, 1]
	! 	m [integer]								: азимутальное число (порядок), m >= 0
	! 	n_max [integer]							: (последний из сетки) индекс (степень), n_max >= m
	! 	norm [character(4)]						: нормировка значения функции -- Фламмера, Майкснера, либо нормированная

	! Выход:
	! 	Smn(m:n_max, size(eta)) [complex(knd)]	: значения функций S_mn(c, eta) для сетки [m, n_max] x [eta_1, eta_k]
	
	integer, intent(in) 				:: n_max, m
	real(knd), intent(in) 				:: eta(:)
	complex(knd), intent(in)			:: c
	character(4), intent(in), optional	:: norm

	integer								:: i, r, s, p, NN, eta_ind, n
	complex(knd), allocatable			:: d_coef(:,:)
	real(knd), allocatable				:: P_Leg(:)

	complex(knd)						:: Smn(m:n_max, size(eta))


	s = mod(n_max-m, 2)
	! p = 3*(n_max-m+1) + (n_max-m+1)/2 
	! if (s /= mod(p, 2)) p = p + 1

	if (n_max-m <= 10) then
		NN = max(4*(n_max-m+1), int(abs(c)))
	else if (n_max-m > 10 .and. n_max-m <= 20) then
		NN = max(3*(n_max-m+1), int(abs(c)))
	else 
		NN = max(2*(n_max-m+1), int(abs(c)))
	end if
	! NN = 1000
	p = max(100, 2*(n_max-m), int(abs(c)))
	! NN = min(1000, 10*((n_max-m+1)+int(abs(c))))
	! NN = n_max-m
	NN = p
	
	if (s /= mod(p, 2)) p = p + 1


	allocate(d_coef((p+2)/2+1, m:n_max), P_Leg(p+2-m+1))
	print*, "d1"
	d_coef = d_coefficients_n(n_max, m, c, p, NN)
	print*, "d2"

	P_Leg = 0

	Smn = 0

	do eta_ind = 1, size(eta)

		P_Leg = Associated_Legendre_polynomials(p+2, m, eta(eta_ind))

		do n=m, n_max
			s = mod(n-m, 2)

			if (eta(eta_ind) == 0) then
				if (s == 0)	Smn(n,eta_ind) = F_func(n-m, n, m, "eta0")
			else		

				do i = 1, (p+2)/2+1

					r = 2*i + (s-1)
					if (r > size(P_Leg)) exit

					Smn(n,eta_ind) = Smn(n,eta_ind) + d_coef(i,n) * P_Leg(r)*(-1)**m
				end do

			end if

			if (present(norm)) then
				if (norm == "Meix")	Smn(n,eta_ind) = Smn(n,eta_ind) / N_mn(m, n, c, p, NN, d_coef(:,n)) &
													* N_mn(m, n, 0*c, p, NN, d_coef(:,n)) 
				if (norm == "Norm")	Smn(n,eta_ind) = Smn(n,eta_ind) / N_mn(m, n, c, p, NN, d_coef(:,n))
			end if

		end do

	end do

	deallocate(d_coef, P_Leg)

end function S_mn_eta_n


function S_mn_eta(c, eta, m, n, norm) result(Smn)
	! Функция для расчета угловых сфероидальных функций для набора 
	! значений eta = [eta_1, eta_k]
	! при заданных m, n, c и norm

	! Вход:
	! 	c [complex(knd)]				: c = kd/2, k -- волновое число в среде, d -- интерфокальное расстояние
	!	eta(eta_1:eta_k) [real(knd)]	: сфероидальная координата (массив) \eta (поверхность -- гиперболоид вращения), [-1, 1]
	! 	m [integer]						: азимутальное число (порядок), m >= 0
	! 	n [integer]						: индекс (степень), n >= m
	! 	norm [character(4)]				: нормировка значения функции -- Фламмера, Майкснера, либо нормированная

	! Выход:
	! 	Smn(size(eta)) [complex(knd)]	: значения функций S_mn(c, eta) для массива [eta_1, eta_k]
	
	integer, intent(in) 				:: n, m
	real(knd), intent(in) 				:: eta(:)
	complex(knd), intent(in)			:: c
	character(4), intent(in), optional	:: norm

	integer								:: i, r, s, p, NN, eta_ind
	complex(knd), allocatable			:: d_coef(:)
	real(knd), allocatable				:: P_Leg(:)

	complex(knd)						:: Smn(size(eta))


	s = mod(n-m, 2)
	! p = 3*(n-m+1) + (n-m+1)/2 
	! if (s /= mod(p, 2)) p = p + 1

	if (n-m <= 10) then
		NN = max(4*(n-m+1), int(abs(c)))
	else if (n-m > 10 .and. n-m <= 20) then
		NN = max(3*(n-m+1), int(abs(c)))
	else 
		NN = max(2*(n-m+1), int(abs(c)))
	end if
	NN = max(100, 2*(n-m), int(abs(c)))

	! NN = 500
	p = max(100, 2*(n-m), int(abs(c)))

	! p = NN
	NN = p
	if (s /= mod(p, 2)) p = p + 1

	! print

	allocate(d_coef((p+2)/2+1), P_Leg(p+2-m+1))
	d_coef = d_coefficients(n, m, c, p, NN)

	Smn = 0

	do eta_ind=1, size(eta)

		if (eta(eta_ind) == 0) then
			if (s == 0)	Smn(eta_ind) = F_func(n-m, n, m, "eta0")
		else		
			P_Leg = Associated_Legendre_polynomials(p+2, m, eta(eta_ind))

			do i=1, (p+2)/2+1
				r = 2*i + (s-1)
				Smn(eta_ind) = Smn(eta_ind) + d_coef(i) * P_Leg(r)*(-1)**m
			end do

		end if

	end do

	if (present(norm)) then
		if (norm == "Meix")	Smn = Smn / N_mn(m, n, c, p, NN, d_coef) * N_mn(m, n, 0*c, p, NN, d_coef) 
		if (norm == "Norm")	Smn = Smn / N_mn(m, n, c, p, NN, d_coef)
	end if

	deallocate(d_coef, P_Leg)

end function S_mn_eta


function S_mn(c, eta, m, n, norm) result(Smn)
	! Функция для расчета угловых сфероидальных функций
	! при заданных m, n, eta, c и norm

	! Вход:
	! 	c [complex(knd)]	: c = kd/2, k -- волновое число в среде, d -- интерфокальное расстояние
	!	eta [real(knd)]		: сфероидальная координата \eta (поверхность -- гиперболоид вращения), [-1, 1]
	! 	m [integer]			: азимутальное число (порядок), m >= 0
	! 	n [integer]			: индекс (степень), n >= m
	! 	norm [character(4)]	: нормировка значения функции -- Фламмера, Майкснера, либо нормированная

	! Выход:
	! 	Smn [complex(knd)]	: значения функций S_mn(c, eta) 
	
	integer, intent(in) 				:: n, m
	real(knd), intent(in) 				:: eta
	complex(knd), intent(in)			:: c
	character(4), intent(in), optional	:: norm

	integer								:: i, r, s, p, NN
	complex(knd), allocatable			:: d_coef(:)
	real(knd), allocatable				:: P_Leg(:)

	complex(knd)						:: Smn


	s = mod(n-m, 2)
	! p = 3*(n-m+1) + (n-m+1)/2 
	! if (s /= mod(p, 2)) p = p + 1

	if (n-m <= 10) then
		NN = max(4*(n-m+1), int(abs(c)))
	else if (m-n > 10 .and. m-n <= 20) then
		NN = max(3*(n-m+1), int(abs(c)))
	else 
		NN = max(2*(n-m+1), int(abs(c)))
	end if

	NN = max(100, 2*(n-m), int(abs(c)))

	! NN = 500
	! p = NN
	NN = p
	if (s /= mod(p, 2)) p = p + 1

	Smn = 0

	if (eta == 0) then
		if (s == 0)	Smn = F_func(n-m, n, m, "eta0")

		if (present(norm)) then
			allocate(d_coef((p+2)/2+1))
			d_coef = d_coefficients(n, m, c, p, NN)

			if (norm == "Meix")	Smn = Smn / N_mn(m, n, c, p, NN, d_coef) * N_mn(m, n, 0*c, p, NN, d_coef) 
			if (norm == "Norm")	Smn = Smn / N_mn(m, n, c, p, NN, d_coef)
		end if

	else
		allocate(d_coef((p+2)/2+1))
		d_coef = d_coefficients(n, m, c, p, NN)

		allocate(P_Leg(p+2-m+1))

		P_Leg = Associated_Legendre_polynomials(p+2, m, eta)
		
		do i=1, (p+2)/2+1
			r = 2*i + (s-1)
			Smn = Smn + d_coef(i) * P_Leg(r)*(-1)**m
		end do

		deallocate(P_Leg)

		if (present(norm)) then
			if (norm == "Meix")	Smn = Smn / N_mn(m, n, c, p, NN, d_coef) * N_mn(m, n, 0*c, p, NN, d_coef) 
			if (norm == "Norm")	Smn = Smn / N_mn(m, n, c, p, NN, d_coef)
		end if

	end if

	if (allocated(d_coef))	deallocate(d_coef)

end function S_mn


function d_coefficients(n, m, c, p, NN) result(d_coef)
	! Функция для расчета коэффициентов d_r(c) разложения 
	! по присоединенным полиномам Лежандра 
	! при заданных m, n, и c

	! Вход:
	! 	m [integer]							: азимутальное число (порядок), m >= 0
	! 	n [integer]							: индекс (степень), n >= m
	! 	c [complex(knd)]					: c = kd/2, k -- волновое число в среде, d -- интерфокальное расстояние
	!	p [integer]							: некоторое большое число, такое что d_p(c) почти 0, p > n-m
	! 	NN [integer]						: размер матрицы рекурентного расчета d_r(c), количество собственных значений

	! Выход:
	! 	d_coef((p+2)/2+1) [complex(knd)]	: значения коэффициентов d_r(c) 

	integer, intent(in) 		:: n, m, p, NN
	complex(knd), intent(in)	:: c

	integer						:: i, s
	complex(knd), allocatable	:: h_coef(:), eig_arr(:)
	complex(knd)				:: g_1

	complex(knd)				:: d_coef((p+2)/2+1)


	s = mod(n-m, 2)
	d_coef = 0

	if (c == 0) then
		i = ((n-m) - s) / 2 + 1
		d_coef(i) = 1
	else
		allocate(eig_arr(m:m+NN-1), h_coef((p+2)/2+1))

		eig_arr = find_eig_corrected(NN, n, m, c)

		h_coef =  h_coefficients_r(n, m, c, eig_arr(n), p)

		g_1 =  g_1_func(n, m, h_coef)

		d_coef = g_1 * h_coef

		deallocate(eig_arr, h_coef)
	end if
	
end function d_coefficients


function d_coefficients_n(n_max, m, c, p_max, NN) result(d_coef)
	! Функция для расчета коэффициентов d_r(c) разложения 
	! по присоединенным полиномам Лежандра 
	! для набора значений n = [m, n_max] при заданных m и c

	! Вход:
	! 	m [integer]										: азимутальное число (порядок), m >= 0
	! 	n_max [integer]									: (последний из сетки) индекс (степень), n_max >= m
	! 	c [complex(knd)]								: c = kd/2, k -- волновое число в среде, d -- интерфокальное расстояние
	!	p_max [integer]									: некоторое большое число, такое что d_p_max(c) почти 0, p_max = max{p_n}, p_n > n-m
	! 	NN [integer]									: размер матрицы рекурентного расчета d_r(c), количество собственных значений

	! Выход:
	! 	d_coef((p_max+2)/2+1, m:n_max) [complex(knd)]	: значения коэффициентов d_r(c)

	integer, intent(in) 		:: n_max, m, p_max, NN
	complex(knd), intent(in)	:: c

	integer						:: i, s, s0, n, p
	complex(knd)				:: h_coef((p_max+2)/2+1), eig_arr(m:m+NN-1)
	complex(knd)				:: g_1

	complex(knd)				:: d_coef((p_max+2)/2+1, m:n_max)


	s0 = mod(n_max-m, 2)

	d_coef = 0

	if (c == 0) then
		do n=m, n_max
			s = mod(n-m, 2)

			i = ((n-m) - s) / 2 + 1
			d_coef(i,n) = 1
		end do
	else

		eig_arr = find_eig_corrected(NN, n_max, m, c)

		do n = m, n_max
			! МУСОР ЛИ ЭТО?? НЕ ПОМНЮ.. НО ВРОДЕ БЫ И ТАК РАБОТАЕТ
			! s = mod(n-m, 2)
			! p = max(100, min(10*(n-m+1), p_max))

			! ! print*, s, p, "pppp"
			! if (s /= mod(p, 2)) p = p + 1
			! ! p = p_max
			! print*, m, n, p, p_max

			! print*, "d. n=", n
			! if (s0 /= mod(n-m, 2)) then
				! h_coef =  h_coefficients_r(n, m, c, eig_arr(n), p-1)
				! g_1 =  g_1_func(n, m, h_coef)
			! else 
				h_coef = h_coefficients_r(n, m, c, eig_arr(n), p_max)
				! print*, "d. lol1"
				g_1 = g_1_func(n, m, h_coef)
				! print*, "d. lol2"

			! end if

			d_coef(:,n) = g_1 * h_coef
			! d_coef(:size(h_coef),n) = g_1 * h_coef

			! ! print*, "d. n=", n
			! ! if (s0 /= mod(n-m, 2)) then
			! 	! h_coef =  h_coefficients_r(n, m, c, eig_arr(n), p-1)
			! 	! g_1 =  g_1_func(n, m, h_coef)
			! ! else 
			! 	h_coef =  h_coefficients_r(n, m, c, eig_arr(n), p_max)
			! 	! print*, "d. lol1"
			! 	g_1 =  g_1_func(n, m, h_coef)
			! 	! print*, "d. lol2"

			! ! end if

			! d_coef(:,n) = g_1 * h_coef

		end do

		! deallocate(eig_arr, h_coef)
	end if
	
end function d_coefficients_n


function Associated_Legendre_polynomials(l, m, x) result(P)
	! Функция для расчета массива присоединенных 
	! полиномов Лежандра P_i^m(x), i = [m:m+l]

	! Вход:
	! 	l [integer]				: m+l -- индекс (степень) многочлена Лежандра, m+l >= m
	! 	m [integer]				: азимутальное число (порядок) многочлена Лежандра, m >= 0
	! 	x [real(knd)]			: аргумент функции, [-1,1]

	! Выход:
	! 	P(m:m+l) [real(knd)]	: значения присоединенных полиномов Лежандра P_i^m(x), i = [m:m+l]

	integer, intent(in)		:: l, m
	real(knd), intent(in)	:: x

	integer					:: i

	real(knd)				:: P(m:m+l)

	P = 0

	! P(m) = (-1)**m * product([(i,i=2*m-1,1,-2)]) * (1-x**2)**(m/2._knd)
	if (.not. ((x == -1 .or. x == 1) .and. m /= 0)) then
		P(m) = 1q0
		do i = 2*m-1, 1, -2
			P(m) = P(m) * i		
		end do
		P(m) = P(m) * (-1)**m * (1q0-x**2)**(m/2._knd)

		P(m+1) = x * (2*m+1) * P(m)

		do i = m+2, m+l
			P(i) = ((2*i-1)*x*P(i-1) - (i-1+m)*P(i-2)) / (i-m)
		end do
	end if
	
end function Associated_Legendre_polynomials
! ФУНКЦИИ 1-ГО ПОРЯДКА ---------------------------------------------- [END]




! ФУНКЦИИ 2-ГО ПОРЯДКА -------------------------------------------- [BEGIN]
! ВЫЧИСЛЕНИЕ СОСТАВЛЯЮЩИХ КОЭФФИЦИЕНТА d_r(c) (ф.38 -- Hodge,1970)
function h_coefficients_r(n, m, c, lambda, p) result(h_coef)
	! Функция для расчета всего набора
	! коэффициентов h_r(c) (~ d_r(c)) при заданных m, n и c

	! Вход:
	! 	n [integer]							: индекс (степень), n_max >= m
	! 	m [integer]							: азимутальное число (порядок), m >= 0
	! 	c [complex(knd)]					: c = kd/2, k -- волновое число в среде, d -- интерфокальное расстояние
	! 	lambda [complex(knd)]				: собственное значение (одно из) матрицы для расчета коэффициентов d_r(c)
	!	p [integer]							: некоторое большое число, такое что d_p(c) почти 0, p > n-m

	! Выход:
	! 	h_coef((p+2)/2+1) [complex(knd)]	: значения коэффициентов h_r(c)

	integer, intent(in) 		:: n, m, p
	complex(knd), intent(in)	:: c, lambda

	complex(knd)				:: hArray(0:(n-m)/2), h_Array(-1:(p-(n-m))/2)

	complex(knd)				:: h_coef((p+2)/2+1)


	hArray = h_bfr_nm(n, m, c, lambda)
	h_Array = h_aft_nm(n, m, c, lambda, p)

	h_Array = hArray((n-m)/2) / h_Array(-1) * h_Array


	h_coef = 0
	h_coef(:(n-m)/2+1) = hArray

	if (size(h_coef((n-m)/2+2:)) == size(h_Array(0:))) then
		h_coef((n-m)/2+2:) = h_Array(0:)
	else
		h_coef((n-m)/2+2:(p+2)/2) = h_Array(0:)
	end if

end function h_coefficients_r


function h_bfr_nm(n, m, c, lambda) result(h_coef)
	! Функция для расчета коэффициентов h_r(c) 
	! при r < n-m и заданных m, n и c

	! Вход:
	! 	n [integer]							: индекс (степень), n_max >= m
	! 	m [integer]							: азимутальное число (порядок), m >= 0
	! 	c [complex(knd)]					: c = kd/2, k -- волновое число в среде, d -- интерфокальное расстояние
	! 	lambda [complex(knd)]				: собственное значение (одно из) матрицы для расчета коэффициентов d_r(c)

	! Выход:
	! 	h_coef(0:(n-m)/2) [complex(knd)]	: значения коэффициентов h_r(c)

	integer, intent(in) 		:: n, m
	complex(knd), intent(in)	:: c, lambda

	integer						:: s, i, r

	complex(knd)				:: h_coef(0:(n-m)/2)


	s = mod(n-m, 2)

	h_coef = 0
	h_coef(0) = 1

	if (size(h_coef) > 1) h_coef(1) = - (B_fun(m, s, c) - lambda) / A_fun(m, s, c)

	if (size(h_coef) > 2) then
		do i = 2, (n-m)/2
			r = 2*i + s
			h_coef(i) = - ((B_fun(m,r-2,c)-lambda)*h_coef(i-1) + C_fun(m,r-2,c)*h_coef(i-2)) / A_fun(m, r-2, c)
		end do
	end if

end function h_bfr_nm


function h_aft_nm(n, m, c, lambda, p) result(h_coef)
	! Функция для расчета вспомогательных коэффициентов h'_r(c) 
	! при r > n-m и заданных m, n и c

	! Вход:
	! 	n [integer]								: индекс (степень), n_max >= m
	! 	m [integer]								: азимутальное число (порядок), m >= 0
	! 	c [complex(knd)]						: c = kd/2, k -- волновое число в среде, d -- интерфокальное расстояние
	! 	lambda [complex(knd)]					: собственное значение (одно из) матрицы для расчета коэффициентов d_r(c)
	!	p [integer]								: некоторое большое число, такое что d_p(c) почти 0, p > n-m

	! Выход:
	! 	h_coef(-1:(p-(n-m))/2) [complex(knd)]	: значения коэффициентов h'_r(c)

	integer, intent(in) 		:: n, m, p
	complex(knd), intent(in)	:: c, lambda

	integer						:: s, r, i

	complex(knd)				:: h_coef(-1:(p-(n-m))/2)


	s = mod(n-m, 2)

	h_coef((p-(n-m))/2) = 0
	h_coef((p-(n-m))/2-1) = 1

	do i = (p-(n-m))/2-2, -1, -1
		r = (n-m+2) + 2*i
		h_coef(i) = - (A_fun(m,r+2,c)*h_coef(i+2) + (B_fun(m,r+2,c)-lambda)*h_coef(i+1)) / C_fun(m, r+2, c)
	end do

end function h_aft_nm


complex(knd) function g_1_func(n, m, h_coef)
	! Функция для расчета нормировочной константы g^-1
	! с помощью нормализации Фламмера

	! Вход:
	! 	n [integer]					: индекс (степень), n_max >= m
	! 	m [integer]					: азимутальное число (порядок), m >= 0
	! 	h_coef(:) [complex(knd)]	: значения коэффициентов h_r(c)

	! Выход:
	! 	g_1_func [complex(knd)]		: значение нормировочной константы

	integer, intent(in) 		:: n, m
	complex(knd), intent(in)	:: h_coef(:)

	integer						:: i, s, r
	complex(knd)				:: summ

	s = mod(n-m, 2)

	summ = 0
	do i = 1, size(h_coef)
		r = 2*i - 2 + s

		summ = summ + F_func(r, n, m) * h_coef(i)
	end do

	g_1_func = F_func(n-m, n, m) / summ

end function g_1_func


real(knd) function F_func(r, n, m, eta0)
	! Функция для расчета нормализации Фламмера F_r^m (для g^-1),
	! либо для расчета полинома Лежандра P_n^m(0) (= S_mn^F(c,0))

	! Вход:
	! 	r [integer]			: (для нормализации Фламмера) индекс, r >= 0
	! 	n [integer]			: (для полинома Лежандра) индекс (степень), n >= m
	! 	m [integer]			: азимутальное число (порядок), m >= 0
	! 	eta0 [character(4)]	: если == "eta0", рассчитывается полином Лежандра

	! Выход:
	! 	F_func [real(knd)]	: значение нормализация Фламмера, либо полинома Лежандра

	integer, intent(in) 				:: r, n, m
	character(4), intent(in), optional 	:: eta0

	integer								:: s, i
	real(16)							:: d1, d2

	s = mod(n-m, 2)

	d1 = 1q0
	do i = (r+2*m+s)/2+1, r+2*m+s
		d1 = d1 * i
	end do
	d1 = d1 * (-1)**((r-s)/2)

	d2 = 1q0
	do i = 1, (r-s)/2
		d2 = d2 * i
	end do

	if (present(eta0)) then
		if (eta0 == "eta0") d2 = d2 * 2q0**n
	else
		d2 = d2 * 2q0**r
	end if

	F_func = d1 / d2

end function F_func
! ФУНКЦИИ 2-ГО ПОРЯДКА ---------------------------------------------- [END]

end module