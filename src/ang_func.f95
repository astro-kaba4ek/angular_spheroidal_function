module ang_func
use regime
use mod_eigenvalues
use matrix_fun
implicit none

contains

function N_mn(m, n, c, p, NN, d_coef) result(summ)
	integer, intent(in) 		:: n, m, p, NN
	complex(knd), intent(in)	:: c
	complex(knd), intent(in)	:: d_coef(:)

	integer						:: i, j , r, s
	complex(knd)				:: sum1

	complex(knd)				:: summ

	s = mod(n-m, 2)

	if (c == 0) then
		
		sum1 = 1q0
		do j=(n-m)+1, (n-m)+2*m
			sum1 = sum1 * j
		end do
		summ = sqrt(sum1 * 2 / (2*m + 2*(n-m) + 1))

	else
		summ = 0
		do i=1, size(d_coef)
			r = 2*(i-1)+s

			sum1 = 1q0
			do j=r+1, r+2*m
				sum1 = sum1 * j
			end do
			sum1 = sum1 * 2 / (2*m + 2*r + 1) * d_coef(i)**2

			summ = summ + sum1
		
		end do
		summ = sqrt(summ)

	end if 

end function N_mn

function S_mn_eta_n(c, eta, m, n, norm) result(summ)
	integer, intent(in) 				:: n, m
	real(knd), intent(in) 				:: eta(:)
	complex(knd), intent(in)			:: c
	character(4), intent(in), optional	:: norm

	integer						:: i, r, s, p, NN, eta_ind, n_ind
	complex(knd), allocatable	:: d_coef(:,:)
	real(knd), allocatable		:: P_Leg(:)

	complex(knd)				:: summ(m:n, size(eta))


	s = mod(n-m, 2)
	! p = 3*(n-m+1) + (n-m+1)/2 
	! if (s /= mod(p, 2)) p = p + 1

	if (m-n <= 10) then
		NN = max(4*(n-m+1), int(abs(c)))
	else if (m-n > 10 .and. m-n <= 20) then
		NN = max(3*(n-m+1), int(abs(c)))
	else 
		NN = max(2*(n-m+1), int(abs(c)))
	end if
	NN = 100
	p = NN
	print*, s, p, "pppp"
	if (s /= mod(p, 2)) p = p + 1
	print*, "p", p

	allocate(d_coef((p+2)/2+1, m:n), P_Leg(p+2-m+1))
	print*, "d1"
	d_coef = d_coef_Arr_eigs(n, m, c, p, NN)
	print*, "d2"

	P_Leg = 0

	summ = 0

	! s = mod(abs(n_ind-m), 2)

	do eta_ind=1, size(eta)

		P_Leg = Associated_Legendre_polynomials_Arr(p+2, m, eta(eta_ind))

		! do i = 1, size(P_Leg)
		! 	print*, i, eta(eta_ind), P_Leg(i)
			
		! end do

		do n_ind=m, n
			! print*,  "n = ", n_ind

			s = mod(n_ind-m, 2)


			if (eta(eta_ind) == 0) then
				if (s == 0)	summ(n_ind,eta_ind) = F_func(n_ind-m, n_ind, m, "eta0")
				print*, "!!"
			else		

				do i=1, (p+2)/2+1

					r = 2*i + (s-1)
					! print*, s, i, r
					if (r > size(P_Leg)) exit

					! summ(eta_ind,n_ind) = summ(eta_ind,n_ind) + d_coef(i,n_ind) * P_Leg(r)
					summ(n_ind,eta_ind) = summ(n_ind,eta_ind) + d_coef(i,n_ind) * P_Leg(r)*(-1)**m
					! write(3,*) m, n_ind, eta(eta_ind), i, r, d_coef(i,n_ind), P_Leg(r), summ(n_ind,eta_ind)
					! write(3,*) m, n_ind, eta(eta_ind), i, r-1, &
					! "---------------------------------------------------------------------------------------------", &
					! P_Leg(r-1), summ(n_ind,eta_ind)
				end do

			end if
			! print*, "#", eta_ind, eta(eta_ind),  n_ind, summ(n_ind,eta_ind) 

			if (present(norm)) then
				if (norm == "Meix")	summ(n_ind,eta_ind) = summ(n_ind,eta_ind) / N_mn(m, n_ind, c, p, NN, d_coef(:,n_ind)) &
														* N_mn(m, n_ind, c*0, p, NN, d_coef(:,n_ind)) 
				if (norm == "Norm")	summ(n_ind,eta_ind) = summ(n_ind,eta_ind) / N_mn(m, n_ind, c, p, NN, d_coef(:,n_ind))
			end if

		end do

	end do

	deallocate(d_coef, P_Leg)

end function S_mn_eta_n

function S_mn_eta(c, eta, m, n, norm) result(summ)
	integer, intent(in) 				:: n, m
	real(knd), intent(in) 				:: eta(:)
	complex(knd), intent(in)			:: c
	character(4), intent(in), optional	:: norm

	integer						:: i, r, s, p, NN, eta_ind
	complex(knd), allocatable	:: d_coef(:)
	real(knd), allocatable		:: P_Leg(:)

	complex(knd)				:: summ(size(eta))


	s = mod(n-m, 2)
	! p = 3*(n-m+1) + (n-m+1)/2 
	! if (s /= mod(p, 2)) p = p + 1

	if (m-n <= 10) then
		NN = max(4*(n-m+1), int(abs(c)))
	else if (m-n > 10 .and. m-n <= 20) then
		NN = max(3*(n-m+1), int(abs(c)))
	else 
		NN = max(2*(n-m+1), int(abs(c)))
	end if
	NN = 50
	p = NN
	if (s /= mod(p, 2)) p = p + 1
	print*, "p", p

	allocate(d_coef((p+2)/2+1), P_Leg(p+2-m+1))
	d_coef = d_coef_Arr(n, m, c, p, NN)

	summ = 0

	do eta_ind=1, size(eta)

		if (eta(eta_ind) == 0) then
			if (s == 0)	summ(eta_ind) = F_func(n-m, n, m, "eta0")
		else		
			P_Leg = Associated_Legendre_polynomials_Arr(p+2, m, eta(eta_ind))

			do i=1, (p+2)/2+1
				r = 2*i + (s-1)
				summ(eta_ind) = summ(eta_ind) + d_coef(i) * P_Leg(r)*(-1)**m
			end do

		end if

	end do

	if (present(norm)) then
		if (norm == "Meix")	summ = summ / N_mn(m, n, c, p, NN, d_coef) * N_mn(m, n, c*0, p, NN, d_coef) 
		if (norm == "Norm")	summ = summ / N_mn(m, n, c, p, NN, d_coef)
	end if

	deallocate(d_coef, P_Leg)

end function S_mn_eta


function S_mn(c, eta, m, n, norm) result(summ)
	integer, intent(in) 				:: n, m
	real(knd), intent(in) 				:: eta
	complex(knd), intent(in)			:: c
	character(4), intent(in), optional	:: norm

	integer						:: i, r, s, p, NN
	complex(knd), allocatable	:: d_coef(:)
	real(knd), allocatable		:: P_Leg(:)

	complex(knd)				:: summ


	s = mod(n-m, 2)
	! p = 3*(n-m+1) + (n-m+1)/2 
	! if (s /= mod(p, 2)) p = p + 1

	if (m-n <= 10) then
		NN = max(4*(n-m+1), int(abs(c)))
	else if (m-n > 10 .and. m-n <= 20) then
		NN = max(3*(n-m+1), int(abs(c)))
	else 
		NN = max(2*(n-m+1), int(abs(c)))
	end if
	NN = 50
	p = NN
	if (s /= mod(p, 2)) p = p + 1
	print*, "p", p

	summ = 0

	if (eta == 0) then
		if (s == 0)	summ = F_func(n-m, n, m, "eta0")

		if (present(norm)) then
			allocate(d_coef((p+2)/2+1))
			d_coef = d_coef_Arr(n, m, c, p, NN)

			if (norm == "Meix")	summ = summ / N_mn(m, n, c, p, NN, d_coef) * N_mn(m, n, c*0, p, NN, d_coef) 
			if (norm == "Norm")	summ = summ / N_mn(m, n, c, p, NN, d_coef)
		end if

	else
		allocate(d_coef((p+2)/2+1))
		d_coef = d_coef_Arr(n, m, c, p, NN)

		! if (m-n <= 10) then
		! 	NN = max(4*(n-m+1), int(abs(c)))
		! else if (m-n > 10 .and. m-n <= 20) then
		! 	NN = max(3*(n-m+1), int(abs(c)))
		! else 
		! 	NN = max(2*(n-m+1), int(abs(c)))
		! end if
		! NN = 50
		! p = NN
		! if (s /= mod(p, 2)) p = p + 1
	
		
		allocate(P_Leg(p+2-m+1))

		P_Leg = Associated_Legendre_polynomials_Arr(p+2, m, eta)
		! print*, P_Leg
		
		print*, "----------------------------------------"
		! print*, "          i         r                 h_coef                                        d_coef         P_leg"
		
		print*, "m", m
		print*, "n", n
		print*, "c", c
		print*, "eta", eta
		print*, "-------"
		print*, "  i   r    d_coef",  d_coef(1),     "                                                            P_leg"
		print*, "-------"
		
		do i=1, (p+2)/2+1
			r = 2*i + (s-1)
			summ = summ + d_coef(i) * P_Leg(r)*(-1)**m
			! print*, i, r, h_coef(i), d_coef(i), P_Leg(r), summ!, Associated_Legendre_polynomials(m+r-1,m,nu)
			write(*,11) i, r, d_coef(i+1)/d_coef(i), d_coef(i), P_Leg(r), P_Leg(r+2)/P_Leg(r)!, Associated_Legendre_polynomials(m+r-1,m,nu)
			11 format(1x,i3, 1x,i3, 1x,e39.30,e39.30, 1x,f39.30,f39.30, 1x,f39.30, 1x,f39.30)
		end do

		deallocate(P_Leg)

		if (present(norm)) then
			if (norm == "Meix")	summ = summ / N_mn(m, n, c, p, NN, d_coef) * N_mn(m, n, c*0, p, NN, d_coef) 
			if (norm == "Norm")	summ = summ / N_mn(m, n, c, p, NN, d_coef)
		end if

	end if

	if (allocated(d_coef))	deallocate(d_coef)

end function S_mn


function d_coef_Arr(n, m, c, p, NN) result(d_coef)
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

		h_coef =  h_coef_Arr(n, m, c, eig_arr(n), p)

		! print*, "----------------------------------------"
		g_1 =  g_1_Arr(n, m, p, h_coef)
		! print*, g_1

		d_coef = g_1 * h_coef
		! do i = 1, size(eig_arr)
		! 	print*, i, eig_arr(i)
		! end do
		! print*, eig_arr

		deallocate(eig_arr, h_coef)
	end if
	
end function d_coef_Arr


function d_coef_Arr_eigs(n_max, m, c, p, NN) result(d_coef)
	integer, intent(in) 		:: n_max, m, p, NN
	complex(knd), intent(in)	:: c

	integer						:: i, s, s0, n
	complex(knd)				:: h_coef((p+2)/2+1), eig_arr(m:m+NN-1)
	complex(knd)				:: g_1

	complex(knd)				:: d_coef((p+2)/2+1, m:n_max)


	s0 = mod(n_max-m, 2)

	! allocate(d_coef((p+2)/2+1, m:m+NN-1)
	d_coef = 0

	if (c == 0) then
		do n=m, n_max
			s = mod(n-m, 2)

			i = ((n-m) - s) / 2 + 1
			d_coef(i,n) = 1
		end do
	else
		! allocate(eig_arr(m:m+NN-1), h_coef((p+2)/2+1))

		eig_arr = find_eig_corrected(NN, n_max, m, c)

		do n=m, n_max
			! print*, "d. n=", n
			! if (s0 /= mod(n-m, 2)) then
				! h_coef =  h_coef_Arr(n, m, c, eig_arr(n), p-1)
				! g_1 =  g_1_Arr(n, m, p-1, h_coef)
			! else 
				h_coef =  h_coef_Arr(n, m, c, eig_arr(n), p)
				! print*, "d. lol1"
				g_1 =  g_1_Arr(n, m, p, h_coef)
				! print*, "d. lol2"

			! end if

			d_coef(:,n) = g_1 * h_coef

		end do

		! deallocate(eig_arr, h_coef)
	end if
	
end function d_coef_Arr_eigs


! function Associated_Legendre_polynomials_Arr2(b, a, x) result(P)
	! 	integer, intent(in)		:: a, b
	! 	real(knd), intent(in)	:: x

	! 	integer					:: i, j

	! 	real(16)				:: P(0:b)

	! 	P = 0

	! 	P(0) = 1q0
	! 	do i=2*a-1, 1, -2
	! 		P(0) = P(0) * i
	! 	end do
	! 	P(0) = P(0) * (-1)**a * (1 - x**2)**(a/2._knd)

	! 	P(1) = x * (2*a+1) * P(0)

	! 	do i=2, b
	! 		j = i + a
	! 		P(i) = (2q0*j-1q0)/(j-a) * x * P(i-1) - (j-1q0+a)/(j-a) * P(i-2)
	! 		! print*, i, (2q0*b-1q0)/(b-a), x , P(i-1), (b-1q0+a)/(b-a), P(i-2)
	! 		print*, i, (2q0*j-1q0)/(j-a), x , P(i-1), (j-1q0+a)/(j-a), P(i-2)
	! 	end do

	! 	print*, "P", P, x

! end function Associated_Legendre_polynomials_Arr2



recursive function Associated_Legendre_polynomials(n, m, x) result(P)
	integer, intent(in)		:: n, m
	real(knd), intent(in)	:: x

	integer					:: i

	real(knd)				:: P


	if (n == m) then
		P = (-1)**m * product([(i,i=2*m-1,1,-2)]) * (1-x**2)**(m/2._knd)
	elseif (n == m+1) then
		P = x * (2*m+1) * Associated_Legendre_polynomials(m, m, x)
	else
		P = ((2*n-1) * x * Associated_Legendre_polynomials(n-1, m, x) - (n-1+m) * Associated_Legendre_polynomials(n-2, m, x)) / (n-m)
	end if
	
end function Associated_Legendre_polynomials

function Associated_Legendre_polynomials_Arr(l, m, x) result(P)
	integer, intent(in)		:: l, m
	real(knd), intent(in)	:: x

	integer					:: i

	real(knd)				:: P(m:m+l)

	P = 0

	! P(m) = (-1)**m * product([(i,i=2*m-1,1,-2)]) * (1-x**2)**(m/2._knd)
	P(m) = 1q0
	do i=2*m-1, 1, -2
		P(m) = P(m) * i		
	end do
	P(m) = P(m) * (-1)**m * (1q0-x**2)**(m/2._knd)

	P(m+1) = x * (2*m+1) * P(m)

	do i=m+2, m+l
		P(i) = ((2*i-1)*x*P(i-1) - (i-1+m)*P(i-2)) / (i-m)
	end do
	
end function Associated_Legendre_polynomials_Arr

! --------------------------------------------------------------------------------------------


function g_1(n, m, c, lambda) result(f)
	integer, intent(in) 		:: n, m
	complex(knd), intent(in)	:: c, lambda

	real(knd)					:: F_r, F_nm
	integer						:: i, p, s

	complex(knd)				:: f, summ

	s = mod(n-m, 2)

	p = 2*(n-m+1) + (n-m+1)/2 
	if (s /= mod(p, 2)) p = p + 1
	print*, "p", p

	summ = 0
	do i=s, p+2, 2
		summ = summ + F_func(i, n, m) * h_coef(i, n, m, c, lambda)
	end do

	F_nm = F_func(n-m, n, m) 

	f = F_nm / summ
	! print*, f

end function g_1



real(knd)	function F_func(r, n, m, eta0) !result(f)
	integer, intent(in) 				:: r, n, m
	character(4), intent(in), optional 	:: eta0

	integer						:: s
	integer					:: i!, d1, d2
	real(16)					::  d1, d2

	s = mod(n-m, 2)

	! product((/(i,i=1,n)/))
	! d1 = (-1)**(r-s) * product([(i,i=1,r+2*m+s)])
	! d2 = 2**r * product([(i,i=1,(r-s)/2)]) * product([(i,i=1,(r+2*m+s)/2)])

	! d1 = (-1)**(r-s) * product([(i, i=(r+2*m+s)/2+1, r+2*m+s)])
	! d2 = 2q0**r * product([(i, i=1, (r-s)/2)]) 

	! print*, "F. lol0"

	d1 = 1q0
	do i=(r+2*m+s)/2+1, r+2*m+s
		d1 = d1 * i
	end do
	d1 = d1 * (-1)**((r-s)/2)

	! print*, "F. lol1"

	d2 = 1q0
	do i=1, (r-s)/2
		d2 = d2 * i
	end do
	! print*, "F. lol2"
	! print*, "F. lol2", present(eta0)
	! print*, "F. lol2", eta0 == "eta0"
	! print*, "F. lol2", present(eta0) .and. eta0 == "eta0"

	if (present(eta0)) then
		if (eta0 == "eta0") d2 = d2 * 2q0**n
	else
		d2 = d2 * 2q0**r
	end if
	! print*, "F. lol3"

	! F_func = d1 / d2

	! d1 = 1q0
	! do i=1, r+2*m+s
	! 	d1 = d1 * i
	! end do
	! d1 = d1 * (-1)**((r-s)/2)


	! d2 = 1q0
	! do i=1, (r-s)/2
	! 	d2 = d2 * i
	! end do
	! do i=1, (r+2*m+s)/2
	! 	d2 = d2 * i
	! end do
	! d2 = d2 * 2q0**r

	F_func = d1 / d2

	! print*, "--", r, n, m, d1, d2

end function F_func

function h_coef(r, n, m, c, lambda) result(f)
	integer, intent(in) 		:: r, n, m
	complex(knd), intent(in)	:: c, lambda

	complex(knd)				:: f

	! print*, "kkkk"

	if (r <= n-m) then
		! write(1, *) "kek0", r
		
		f = h(r, n, m, c, lambda)
		! write(1, *) "kek20", r
		! print*, "kkkk", f

	else
		! print*, "lol"

		f = h(n-m, n, m, c, lambda) / h_(n-m, n, m, c, lambda) * h_(r, n, m, c, lambda)
		! print*, "lollol", f
		
	end if	
	! print*, r, f, h_(r, n, m, c, lambda)
end function h_coef

recursive function h(r, n, m, c, lambda) result(f)
	integer, intent(in) 		:: r, n, m
	complex(knd), intent(in)	:: c, lambda

	integer						:: s

	complex(knd)				:: f

	! write(1, *) "kek", r
	s = mod(n-m, 2)

	if (r == s) then
		f = 1
	elseif (r == s+2) then
		! f = - (B_fun_b(m, r-2, c) - lambda) / A_fun_b(m, r-2, c)
		f = - (B_fun(m, r-2, c) - lambda) / A_fun(m, r-2, c)
		! print*, "a", f, B_fun(m, r-2, c) - lambda,  A_fun(m, r-2, c)
	else
		! f = - ((B_fun_b(m, r-2, c) - lambda)*h(r-2, n, m, c, lambda) + C_fun_b(m, r-2, c)*h(r-4, n, m, c, lambda)) / A_fun_b(m, r-2, c)
		f = - ((B_fun(m, r-2, c) - lambda)*h(r-2, n, m, c, lambda) + C_fun(m, r-2, c)*h(r-4, n, m, c, lambda)) / A_fun(m, r-2, c)
		! print*, "b", f

	end if
	! write(1, *) "kek2", r

end function h

function g_1_Arr(n, m, p, h_coef) result(f)
	integer, intent(in) 		:: n, m, p
	complex(knd), intent(in)	:: h_coef((p+2)/2+1)

	real(knd)					:: F_r, F_nm, ff
	integer						:: i, s, r

	complex(knd)				:: f, summ

	s = mod(n-m, 2)

	! p = 2*(n-m+1) + (n-m+1)/2 
	! if (s /= mod(p, 2)) p = p + 1
	! print*, "p", p

	summ = 0
	! print*, "          i           r                F_mr                               h_coef"
	do i=1, (p+2)/2+1
		! print*, "g_1. i=", i
		r = 2*i - 2 + s
		! print*, "g_1. lol0"

		ff = F_func(r, n, m)
		! print*, "g_1. lol1"
		summ = summ + ff * h_coef(i)
		! print*, "g_1. lol2"

		! print*, i, r,  ff, h_coef(i), "summ", summ
	end do
	! print*, "g_1"

	F_nm = F_func(n-m, n, m) 
	! summ = summ - 45
	f = F_nm / summ
	! f = 1/2.35475616775390750620396589815889165
	! print*, F_nm, summ
	! f = -6.5319660505q-04

end function g_1_Arr

function h_coef_Arr(n, m, c, lambda, p) result(f)
	integer, intent(in) 		:: n, m, p
	complex(knd), intent(in)	:: c, lambda

	complex(knd)				:: hArray(0:(n-m)/2), h_Array(-1:(p-(n-m))/2)

	complex(knd)				:: f((p+2)/2+1)

	! print*, "kkkk"
	hArray = hArr(n, m, c, lambda)
	h_Array = h_Arr(n, m, c, lambda, p)

	h_Array = hArray((n-m)/2) / h_Array(-1) * h_Array

	! print*, n, m, n-m

	! print*, size(hArray), size(h_Array), size(f)
	! print*, size(hArray), size(h_Array(0:))
	! print*, size(f(:(n-m)/2+1)), size(f((n-m)/2+2:))
	! print*, hArray
	! print*, h_Array
	! print*, h_Array(0:)
	f = 0


	f(:(n-m)/2+1) = hArray

	if (size(f((n-m)/2+2:)) == size(h_Array(0:))) then
		f((n-m)/2+2:) = h_Array(0:)
	else
		f((n-m)/2+2:(p+2)/2) = h_Array(0:)
	end if
	! print*, f
	! print*, size(f)

	! print*, hArray
	! print*, h_Array

	! if (r <= n-m) then
	! 	! write(1, *) "kek0", r
		
	! 	f = h(r, n, m, c, lambda)
	! 	! write(1, *) "kek20", r
	! 	! print*, "kkkk", f

	! else
	! 	! print*, "lol"

	! 	f = h(n-m, n, m, c, lambda) / h_(n-m, n, m, c, lambda) * h_(r, n, m, c, lambda)
	! 	! print*, "lollol", f
		
	! end if	
	! print*, r, f, h_(r, n, m, c, lambda)
end function h_coef_Arr

function hArr(n, m, c, lambda) result(f_arr)
	integer, intent(in) 		:: n, m
	complex(knd), intent(in)	:: c, lambda

	integer						:: s, i, r

	complex(knd)				:: f
	complex(knd)				:: f_arr(0:(n-m)/2)
	! logical						:: f_log_arr(0:(n-m)/2)


	! f_log_arr = .False.

	s = mod(n-m, 2)

	f_arr = 0


	f_arr(0) = 1
	if (size(f_arr) > 1) f_arr(1) = - (B_fun(m, s, c) - lambda) / A_fun(m, s, c)

	if (size(f_arr) > 2) then
		do i=2, (n-m)/2
			r = 2*i + s
			f_arr(i) = - ((B_fun(m, r-2, c) - lambda)*f_arr(i-1) + C_fun(m, r-2, c)*f_arr(i-2)) / A_fun(m, r-2, c)
		end do
	end if

end function hArr


function h_Arr(n, m, c, lambda, p) result(f_arr)
	integer, intent(in) 		:: n, m, p
	complex(knd), intent(in)	:: c, lambda

	integer						:: s, r, i

	complex(knd)				:: f

	complex(knd)				:: f_arr(-1:(p-(n-m))/2)
	! logical						:: f_log_arr(0:(n-m)/2)


	s = mod(n-m, 2)
	! p = 2*(n-m+1) + (n-m+1)/2 
	! if (s /= mod(p, 2)) p = p + 1
	f_arr((p-(n-m))/2) = 0
	f_arr((p-(n-m))/2-1) = 1
	! f_arr((p-(n-m))/2) = 1
	! f_arr((p-(n-m))/2-1) = - ((B_fun(m, p, c)-lambda)*f_arr((p-(n-m))/2)) / C_fun(m, p, c)

	! print*, (p-(n-m))/2, (n-m+2) + 2*(p-(n-m))/2, f_arr((p-(n-m))/2)
	! print*, (p-(n-m))/2-1, (n-m+2) + 2*((p-(n-m))/2-1), f_arr((p-(n-m))/2-1)


	do i=(p-(n-m))/2-2, -1, -1
		r = (n-m+2) + 2*i
		
		f_arr(i) = - (A_fun(m, r+2, c)*f_arr(i+2) + (B_fun(m, r+2, c)-lambda)*f_arr(i+1)) / C_fun(m, r+2, c)
		! print*, i, r, f_arr(i+2), f_arr(i+1), f_arr(i)

	end do

	! f_arr = h_mn / f_arr(-1) * f_arr
	! print*, f_arr



end function h_Arr


recursive function h_(r, n, m, c, lambda) result(f)
	integer, intent(in) 		:: r, n, m
	complex(knd), intent(in)	:: c, lambda

	integer						:: s, p

	complex(knd)				:: f

	s = mod(n-m, 2)
	p = 2*(n-m+1) + (n-m+1)/2 
	if (s /= mod(p, 2)) p = p + 1



	if (r >= p+2) then
		f = 0
		! print*, "lol1", r, f
		write(1, *) "lol1", r, f


	else if (r == p) then
		f = 1
		! f = 1
		! print*, "lol2", r, f
		write(1, *) "lol2", r, f


	else
		! print*, "lol31", r, f
		write(1, *) "lol31", r, f

		! f = - (A_fun_b(m, r+2, c)*h_(r+4, n, m, c, lambda) + (B_fun_b(m, r+2, c)-lambda)*h_(r+2, n, m, c, lambda)) / C_fun_b(m, r+2, c)
		f = - (A_fun(m, r+2, c)*h_(r+4, n, m, c, lambda) + (B_fun(m, r+2, c)-lambda)*h_(r+2, n, m, c, lambda)) / C_fun(m, r+2, c)
		! print*, "lol3", r, f
		write(1, *) "lol3", r, f

	
	end if
	! print*, "lol", f

end function h_

end module