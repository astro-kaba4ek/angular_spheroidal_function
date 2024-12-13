module mod_eigenvalues
use regime
use matrix_fun
use bouwkamp
use utils
implicit none

contains


recursive function characteristic_polynomial(lambda, k, mtr) result(res)
	real(knd), intent(in) :: lambda, mtr(:,:)
	integer, intent(in) :: k
	! real(knd), intent(out) :: res
	real(knd) res

	if (k == 0) then
		res = 1
	else if (k == 1) then
		res = mtr(1,1) - lambda
	else 
		res = (mtr(k,1)-lambda) * characteristic_polynomial(lambda, k-1, mtr) &
		 - mtr(k,2)**2 * characteristic_polynomial(lambda, k-2, mtr)
	end if

end function characteristic_polynomial


function bisecn2(aa,bb,epsx,epsy, k, mtr) result(res)
	! function bisecn2(aa,bb,epsx,epsy,kk,ier, k, mtr) result(res)
	real(knd) aa, bb, epsx, epsy, res, f, a, b, c, fa, fb, fc, mtr(:,:)
	integer kk, ier, k, kmax / 10000 / 
	a=aa; b=bb
	fa = characteristic_polynomial(a, k, mtr)
	fb = characteristic_polynomial(b, k, mtr)

	! if ( (fa*fb).le.0_knd) then    ! Если f(x) перескает ось абсцисс, то:
		kk=0                           ! Обнуляем число сужений промежутка.
		do
			c=(a+b)*0.5_knd
			fc=characteristic_polynomial(c, k, mtr)
				! Находим середину промежутка (с) и f(c).
			if ( fa*fc.gt.0_knd) then       ! Если корень в правой половине отрезка,
				a=c; fa=fc
			else              ! смещаем левую границу
				b=c; fb=fc                    ! если в левой - смещаем правую.
			endif !

			kk=kk+1                         ! Увеличиваем счетчик сужений
			if ((kk.gt.kmax).or.&          ! Если число сужений больше допустимого
			& (abs(b-a).lt.epsx.and.&    ! ИЛИ одновременно достигнуты точности
			& abs(fc ).lt.epsy)) exit    ! и по аргументу и по функции, то
		enddo   							! выходим из цикла.

		res=(a+b)*0.5_knd               ! Фиксируем найденное значение корня.

		if (kk.le.kmax) then
			ier=0    ! Если число сужений в норме, то код = 0
		else
			ier=1                   ! иначе код завершения = 1
		endif
	! else
	! 	ier=2                   ! Если f(x) НЕ ПЕРЕСЕКАЕТ ось абсцисс, код = 2
	! endif
end function bisecn2


! собственные значения симметричной 3х-диагональной (в компактном виде) матрицы (методом бисекции)
function eigenvalues_of_matrix(mtr, n) result(res_arr)
	integer, intent(in) :: n
	real(knd), intent(in) :: mtr(n,2)
	real(knd) :: a1, b1, res, intervals(2,n), mtr_norm(n,3)
	integer :: i, k, kk, ier
	! real(knd), intent(out) :: res_arr(n)
	real(knd) :: res_arr(n)


	mtr_norm(:,2:) = mtr
	mtr_norm(:,1) = mtr(:,2)

	a1 = -norm2(mtr_norm)
	b1 = -a1

	intervals(:,1) = [a1, b1]
	! print*, "lol", intervals(:,1)


	do k=1, n

		! print*, k, intervals(1,:)
		! print*, k, intervals(2,:)
		do i=1, k
			! call bisecn1(intervals(1,i),intervals(2,i),1e-5_knd,1e-5_knd,res,kk,ier, k, mtr)
			res = bisecn2(intervals(1,i),intervals(2,i),1e-15_knd,1e-15_knd, k, mtr)
			res_arr(i) = res
		end do

		if (k == n) then
			exit
		else
			do i=1, k
				res = res_arr(i)

				if (i == 1 .and. i == k) then
					intervals(:,i) = [a1, res]
					intervals(:,i+1) = [res, b1]

				else if (i == 1) then 
					intervals(:,i) = [a1, res]
				else if (i == k) then
					intervals(:,i) = [res_arr(i-1), res]
					intervals(:,i+1) = [res, b1]
				else 
					intervals(:,i) = [res_arr(i-1), res]
				end if
			end do
		end if

		! print*, k, res_arr(:k)
		! print*, k, intervals(1,:)
		! print*, k, intervals(2,:)
		! print*, ""

	end do

end function


! function eig_by_cgeev(mtr_mini, n, info) result(eig_arr)
! 	! INFO is INTEGER
! 	! = 0:  successful exit
! 	! < 0:  if INFO = -i, the i-th argument had an illegal value.
! 	! > 0:  if INFO = i, the QR algorithm failed to compute all the
! 	! 		eigenvalues, and no eigenvectors have been computed;
! 	! 		elements i+1:N of W contain eigenvalues which have
! 	! 		converged.
! 	integer, intent(in)			:: n
! 	complex(knd), intent(in)	:: mtr_mini(n,2)

! 	integer						:: lwork, ldvl, ldvr
! 	real(knd)					:: rwork(2*n)
! 	complex(knd) 				:: mtr(n,n) 
! 	complex(knd), allocatable	:: work(:), vl(:,:), vr(:,:)

! 	integer, intent(out)		:: info
! 	complex(knd) 				:: eig_arr(n)


! 	mtr = from_sim3diagMini_to_sim3diagBase(mtr_mini, n)

! 	ldvl = 1
! 	ldvr = 1

! 	allocate(work(1), vl(ldvl,n), vr(ldvr,n))

! 	if (knd == 4) then
! 		call cgeev("N", "N", n, mtr, n, eig_arr, vl, ldvl, vr, ldvr, work, -1, rwork, info)
! 	else 
! 		print*, "lol2"
! 		call zgeev("N", "N", n, mtr, n, eig_arr, vl, ldvl, vr, ldvr, work, -1, rwork, info)
! 	end if

! 	lwork = int(work(1))
! 	deallocate(work)


! 	if (lwork == 0) lwork = max(1, 2*n) + 100
	
! 	print*, "lol", lwork, size(work), info



! 	mtr = from_sim3diagMini_to_sim3diagBase(mtr_mini, n)

! 	allocate(work(max(1,lwork)))
! 	! lwork = size(work)
! 	print*, "lol", lwork,  size(work), info

! 	if (knd == 4) then
! 		call cgeev("N", "N", n, mtr, n, eig_arr, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
! 	else 
! 		print*, "lol2"
! 		call zgeev("N", "N", n, mtr, n, eig_arr, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
! 	end if

! 	deallocate(work, vl, vr)

	
! end function eig_by_cgeev




function eig_by_cgeev2(mtr_mini, n, info) result(eig_arr)
	! INFO is INTEGER
	! = 0:  successful exit
	! < 0:  if INFO = -i, the i-th argument had an illegal value.
	! > 0:  if INFO = i, the QR algorithm failed to compute all the
	! 		eigenvalues, and no eigenvectors have been computed;
	! 		elements i+1:N of W contain eigenvalues which have
	! 		converged.
	integer, intent(in)			:: n
	complex(8), intent(in)	:: mtr_mini(n,2)

	integer						:: lwork, ldvl, ldvr
	real(8)					:: rwork(2*n)
	complex(8) 				:: mtr(n,n) 
	complex(8), allocatable	:: work(:), vl(:,:), vr(:,:)

	integer, intent(out)		:: info
	complex(8) 				:: eig_arr(n)


	mtr = from_sim3diagMini_to_sim3diagBase_8(mtr_mini, n)

	ldvl = 1
	ldvr = 1

	allocate(work(1), vl(ldvl,n), vr(ldvr,n))

	! if (knd == 4) then
	! 	call cgeev("N", "N", n, mtr, n, eig_arr, vl, ldvl, vr, ldvr, work, -1, rwork, info)
	! else 
		! print*, "lol2"
		call zgeev("N", "N", n, mtr, n, eig_arr, vl, ldvl, vr, ldvr, work, -1, rwork, info)
	! end if

	lwork = int(work(1))
	deallocate(work)


	if (lwork == 0) lwork = max(1, 2*n) + 100
	
	print*, "lol", lwork, size(work), info



	mtr = from_sim3diagMini_to_sim3diagBase_8(mtr_mini, n)

	allocate(work(max(1,lwork)))
	! lwork = size(work)
	print*, "lol", lwork,  size(work), info

	! if (knd == 4) then
	! 	call cgeev("N", "N", n, mtr, n, eig_arr, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
	! else 
	! 	print*, "lol2"
		call zgeev("N", "N", n, mtr, n, eig_arr, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
	! end if

	deallocate(work, vl, vr)

	
end function eig_by_cgeev2



function eig_by_cgeev_base(mtr_in, n, info) result(eig_arr)
	! INFO is INTEGER
	! = 0:  successful exit
	! < 0:  if INFO = -i, the i-th argument had an illegal value.
	! > 0:  if INFO = i, the QR algorithm failed to compute all the
	! 		eigenvalues, and no eigenvectors have been computed;
	! 		elements i+1:N of W contain eigenvalues which have
	! 		converged.
	integer, intent(in)		:: n
	complex(8), intent(in)	:: mtr_in(n,n)

	integer					:: lwork, ldvl, ldvr, i,j
	real(8)					:: rwork(2*n)
	complex(8) 				:: mtr(n,n) 
	complex(8), allocatable	:: work(:), vl(:,:), vr(:,:)

	integer, intent(out)	:: info
	complex(8) 				:: eig_arr(n)


	! mtr = from_sim3diagMini_to_sim3diagBase_8(mtr_mini, n)
	mtr = mtr_in

	ldvl = 1
	ldvr = 1

	allocate(work(1), vl(ldvl,n), vr(ldvr,n))

	! if (knd == 4) then
	! 	call cgeev("N", "N", n, mtr, n, eig_arr, vl, ldvl, vr, ldvr, work, -1, rwork, info)
	! else 
		! print*, "lol2"
		call zgeev("N", "N", n, mtr, n, eig_arr, vl, ldvl, vr, ldvr, work, -1, rwork, info)
	! end if

	lwork = int(work(1))
	deallocate(work)


	if (lwork == 0) lwork = max(1, 2*n) + 100
	
	print*, "lol", lwork, size(work), info



	! mtr = from_sim3diagMini_to_sim3diagBase_8(mtr_mini, n)
	! do i=1, N
	! 	write(*,*) (mtr_in(i,j), j=1, N)
	! end do

	mtr = mtr_in


	allocate(work(max(1,lwork)))
	! lwork = size(work)
	print*, "lol", lwork,  size(work), info

	! if (knd == 4) then
	! 	call cgeev("N", "N", n, mtr, n, eig_arr, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
	! else 
	! 	print*, "lol2"
		call zgeev("N", "N", n, mtr, n, eig_arr, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
	! end if

	deallocate(work, vl, vr)

	
end function eig_by_cgeev_base



function find_eig(mtr_mini_all, NN, n, m, info) result(eig_arr)
		! INFO is INTEGER
	! = 0:  successful exit
	! < 0:  if INFO = -i, the i-th argument had an illegal value.
	! > 0:  if INFO = i, the QR algorithm failed to compute all the
	! 		eigenvalues, and no eigenvectors have been computed;
	! 		elements i+1:N of W contain eigenvalues which have
	! 		converged.
	integer, intent(in)		:: NN, n, m 
	complex(8), intent(in)	:: mtr_mini_all(NN/2,4) ! == DFE_matr_mini_b_full(NN, n, m, c)

	integer					:: lwork, ldvl, ldvr, i, j, A_ind(NN), T_ind(NN)
	real(8)					:: rwork(2*NN)
	complex(8) 				:: mtr(NN,NN) 
	complex(8), allocatable	:: work(:), vl(:,:), vr(:,:)

	integer, intent(inout)	:: info
	complex(8) 				:: eig_arr(NN), eig_arr1(NN/2), eig_arr2(NN/2), eig, T((NN+1)/2), T12((NN/2+1)/2)

	! print*, "hhh", mtr_mini_all(:,:2)
	! print*, "hhh", mtr_mini_all(:,3:)

	! eig_arr1 = -eig_by_cgeev2(mtr_mini_all(:,:2), NN/2, info)
	! eig_arr2 = -eig_by_cgeev2(mtr_mini_all(:,3:), NN/2, info)
	eig_arr1 = eig_by_cgeev2(mtr_mini_all(:,:2), NN/2, info)
	eig_arr2 = eig_by_cgeev2(mtr_mini_all(:,3:), NN/2, info)

	! print*, "t0", eig_arr1
	! print*, "t0", eig_arr2

	if (real(eig_arr1(1)) > real(eig_arr1(NN/2))) eig_arr1 = eig_arr1(NN/2:1:-1)
	if (real(eig_arr2(1)) > real(eig_arr2(NN/2))) eig_arr2 = eig_arr2(NN/2:1:-1)

	! do i=1, NN/2-1
	! 	if (real(eig_arr1(i)) > real(eig_arr1(i+1))) then
	! 		eig = eig_arr1(i)
	! 		eig_arr1(i) = eig_arr1(i+1)
	! 		eig_arr1(i+1) = eig
	! 	end if
	! 	if (real(eig_arr2(i)) > real(eig_arr2(i+1))) then
	! 		eig = eig_arr2(i)
	! 		eig_arr2(i) = eig_arr2(i+1)
	! 		eig_arr2(i+1) = eig
	! 	end if
	! end do



	call MergeSort8(eig_arr1, NN/2, T12)
	call MergeSort8(eig_arr2, NN/2, T12)
	! print*, "t0", eig_arr1
	! print*, "t0", eig_arr2




	if (mod(n-m, 2) == 0) then
		! do i=1, NN/2
		! 	eig_arr(2*i-1) = eig_arr1(i)
		! 	eig_arr(2*i) = eig_arr2(i)
		! end do
		eig_arr(1:NN:2) = eig_arr1
		eig_arr(2:NN:2) = eig_arr2
	else
		! do i=1, NN/2
		! 	eig_arr(2*i-1) = eig_arr2(i)
		! 	eig_arr(2*i) = eig_arr1(i)
		! end do
		eig_arr(1:NN:2) = eig_arr2
		eig_arr(2:NN:2) = eig_arr1
	end if



	! ! do j=1, NN+10
	! 	do i=1, NN-1
	! 		if (real(eig_arr(i)) > real(eig_arr(i+1))) then
	! 			eig = eig_arr(i)
	! 			eig_arr(i) = eig_arr(i+1)
	! 			eig_arr(i+1) = eig
	! 		end if
	! 	end do
	! ! end do
	call MergeSort8(eig_arr, NN, T)



	! print*, "tt", eig_arr1
	! print*, "tt", eig_arr2
	! print*, "tt", eig_arr1(ubound(eig_arr1,dim=1)::-1)
	! print*, "tt", eig_arr1(NN/2:1:-1)


	
end function find_eig


function find_eig_corrected(NN, n, m, c) result(eig_arr_corr)
	integer, intent(in)			:: NN, n, m 
	complex(knd), intent(in)	:: c

	integer						:: i, info, j, count
	complex(knd)				:: eig_arr(m:m+NN-1), mtr_mini_all(NN/2,4), dl, eig, dl0, T((NN+1)/2)
	complex(8)					:: eig_arr_8(m:m+NN-1), mtr_mini_all_8(NN/2,4)

	complex(knd)				:: eig_arr_corr(m:m+NN-1), eig_arr_corr0(m:m+NN-1)


	if (c == complex(0_knd, 0_knd)) then
		eig_arr_corr = [complex :: (n*(n+1), n=m, m+NN-1)]
	else
		! mtr_mini_all = DFE_matr_mini_b_full(NN, n, m, c) 
		mtr_mini_all = DFE_matr_mini_full(NN, n, m, c) 
		mtr_mini_all_8 = mtr_mini_all


		eig_arr_8 = find_eig(mtr_mini_all_8, NN, n, m, info)
		eig_arr = eig_arr_8

		! do i = 1, NN
		! 	dl = delta_lambda(m, c, eig_arr(i), i-1)
		! 	print*, i-1, eig_arr(i), eig_arr(i)+dl, dl
		! 	eig_arr_corr(i) = eig_arr(i) + dl
		! end do

		! print*, "          i", "           n", "             eigenvalue", &
		! "                                             true_eigenvalue", &
		! "                                 delta"

		do i = m, m+NN-1
			! if (i == 11) then
			! 	dl = delta_lambda(m, c, eig_arr(13), 13-m)
			! 	eig_arr_corr(i) = eig_arr(13) + dl
			! else if (i == 13) then 
			! 	dl = delta_lambda(m, c, eig_arr(11), 11-m)
			! 	eig_arr_corr(i) = eig_arr(11) + dl
			! else
				dl0 = delta_lambda(m, c, eig_arr(i), i-m)
				eig_arr_corr(i) = eig_arr(i) + dl0
				count = 1
				do while (.true.)
					count = count + 1
					dl = delta_lambda(m, c, eig_arr_corr(i), i-m)
					eig_arr_corr(i) = eig_arr_corr(i) + dl
					! if (abs(dl)/abs(eig_arr_corr(i)) < 0.001_knd .or. count > 50) exit
					if (abs(dl) < 1q-12 .or. count > 100) exit
				end do
				if (abs(dl) > 0.1) eig_arr_corr(i) = eig_arr(i)

				! dl = delta_lambda(m, c, eig_arr_corr(i), i-m)
				! eig_arr_corr(i) = eig_arr_corr(i) + dl
				! dl = delta_lambda(m, c, eig_arr_corr(i), i-m)
				! eig_arr_corr(i) = eig_arr_corr(i) + dl

			! end if
			! write(25,5) i-m, i, eig_arr(i), eig_arr_corr(i), dl, dl0
			! 5	format(1x,i3, 1x,i3," |", 1x,f40.30,f40.30," |", 1x,f40.30,f40.30," |", 1x,f40.30,f40.30," |", 1x,f40.30,f40.30)
			! dl = delta_lambda(m, c, eig_arr(i), i-m)
			! print*, i-m, i, eig_arr(i)+dl, dl
			! ! print*, i, eig_arr(i)+dl
			
			! eig_arr_corr(i) = eig_arr(i) + dl
		end do
		! call MergeSort(eig_arr_corr, NN, T)


		! eig_arr_corr0 = eig_arr_corr
		! eig_arr_corr(10) = eig_arr_corr0(12)
		! eig_arr_corr(11) = eig_arr_corr0(13)
		! eig_arr_corr(12) = eig_arr_corr0(10)
		! eig_arr_corr(13) = eig_arr_corr0(11)
		do i = m, m+NN-1
			
			write(25,5) m, i, eig_arr(i), eig_arr_corr(i), dl
			5	format(1x,i3, 1x,i3," |", 1x,f40.30,f40.30," |", 1x,f40.30,f40.30," |", 1x,f40.30,f40.30)
			
		end do

		! eig_arr_corr = eig_arr


		! do j=1, NN+10
		! 	do i=1, NN-1
		! 		if (real(eig_arr(i)) > real(eig_arr(i+1))) then
		! 			eig = eig_arr(i)
		! 			eig_arr(i) = eig_arr(i+1)
		! 			eig_arr(i+1) = eig
		! 		end if
		! 	end do
		! end do
		
	end if
	
end function find_eig_corrected


end module