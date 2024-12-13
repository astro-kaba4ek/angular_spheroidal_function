module utils
	! Модуль со всевозможными вспомогательными утилитами
use regime
implicit none
	
contains
	

! subroutine Merge(A, A_ind, NA, B, B_ind, NB, C, C_ind, NC)

! 	integer, intent(in)			:: NA, NB, NC       ! Normal usage: NA+NB = NC
! 	real(8), intent(inout)	:: A(NA)        	! B overlays C(NA+1:NC)
! 	real(8), intent(in)		:: B(NB)
! 	real(8), intent(inout)	:: C(NC)

! 	integer, intent(inout)	:: A_ind(NA)        	
! 	integer, intent(in)		:: B_ind(NB)
! 	integer, intent(inout)	:: C_ind(NC)

! 	integer :: I, J, K

! 	I = 1
! 	J = 1
! 	K = 1

! 	do while(I <= NA .and. J <= NB)
! 		if (A(I) <= B(J)) then
! 			C(K) = A(I)
! 			C_ind(K) = A_ind(I)
! 			I = I+1
! 		else
! 			C(K) = B(J)
! 			C_ind(K) = B_ind(J)
! 			J = J+1
! 		endif
! 		K = K + 1
! 	enddo

! 	do while (I <= NA)
! 		C(K) = A(I)
! 		C_ind(K) = A_ind(I)
! 		I = I + 1
! 		K = K + 1
! 	enddo

! end subroutine Merge

! recursive subroutine MergeSort(A, A_ind, N, T, T_ind)

! 	integer, intent(in)							:: N
! 	real(8), dimension(N), intent(inout)		:: A
! 	integer, dimension(N), intent(inout)		:: A_ind
! 	real(8), dimension((N+1)/2), intent (out)	:: T
! 	integer, dimension((N+1)/2), intent(out)	:: T_ind


! 	integer :: NA, NB, V, V_ind

! 	if (N < 2) return

! 	if (N == 2) then
! 		if (A(1) > A(2)) then
! 			V = A(1)
! 			A(1) = A(2)
! 			A(2) = V

! 			V_ind = A_ind(1)
! 			A_ind(1) = A_ind(2)
! 			A_ind(2) = V_ind
! 		endif
! 		return
! 	endif      
! 	NA = (N+1) / 2
! 	NB = N - NA

! 	call MergeSort(A, A_ind, NA, T, T_ind)
! 	call MergeSort(A(NA+1), A_ind(NA+1), NB, T, T_ind)

! 	if (A(NA) > A(NA+1)) then
! 		T(1:NA) = A(1:NA)
! 		T_ind(1:NA) = A_ind(1:NA)
! 		call Merge(T, T_ind, NA, A(NA+1), A_ind(NA+1), NB, A, A_ind, N)
! 	endif

! end subroutine MergeSort




recursive subroutine MergeSort(A, N, T)
	! Рекурсивная подпрограмма для сортировки комплексного массива
	! по возрастанию вещественной части (можно переделать 
	! на мнимую / модуль / др.) элементов
	! Идея была украдена непомню-откуда..

	! Вход:
	! 	N [integer]				: размер массива A 
	! 	A(N) [complex(knd)]		: массив A, который нужно отсортировать

	! Выход:
	! 	T(N+1)/2) [complex(knd)]: размер массива A 
	! 	A(N) [complex(knd)]		: отсортированный массив A

	integer, intent(in)			:: N
	! complex(knd), dimension(N), intent(inout)		:: A
	! complex(knd), dimension((N+1)/2), intent(out)	:: T

	complex(knd)				:: V
	integer 					:: NA, NB

	complex(knd), intent(inout)	:: A(N)
	complex(knd), intent(out)	:: T((N+1)/2)

	if (N < 2) return

	if (N == 2) then
		if (real(A(1)) > real(A(2))) then
			V = A(1)
			A(1) = A(2)
			A(2) = V
		endif
		return
	endif      
	NA = (N+1) / 2
	NB = N - NA

	call MergeSort(A, NA, T)
	call MergeSort(A(NA+1), NB, T)

	if (real(A(NA)) > real(A(NA+1))) then
		T(1:NA) = A(1:NA)
		call Merge(T, NA, A(NA+1), NB, A, N)
	endif

end subroutine MergeSort


subroutine Merge(A, NA, B, NB, C, NC)
	! Вспомогательная подпрограмма для MergeSort,
	! которая что-то там переставляет...

	integer, intent(in)			:: NA, NB, NC       ! Normal usage: NA+NB = NC
	complex(knd), intent(in)	:: B(NB)

	integer 					:: I, J, K

	complex(knd), intent(inout)	:: C(NC)
	complex(knd), intent(inout)	:: A(NA)        	! B overlays C(NA+1:NC)

	I = 1
	J = 1
	K = 1

	do while(I <= NA .and. J <= NB)
		if (real(A(I)) <= real(B(J))) then
			C(K) = A(I)
			I = I + 1
		else
			C(K) = B(J)
			J = J + 1
		endif
		K = K + 1
	enddo

	do while (I <= NA)
		C(K) = A(I)
		I = I + 1
		K = K + 1
	enddo

end subroutine Merge

! ---------------------------------------------------------------------------------------

recursive subroutine MergeSort8(A, N, T)
	! Рекурсивная подпрограмма для сортировки комплексного массива
	! по возрастанию вещественной части (можно переделать 
	! на мнимую / модуль / др.) элементов
	! Идея была украдена непомню-откуда..

	! Вход:
	! 	N [integer]				: размер массива A 
	! 	A(N) [complex(8)]		: массив A, который нужно отсортировать

	! Выход:
	! 	T(N+1)/2) [complex(8)]: размер массива A 
	! 	A(N) [complex(8)]		: отсортированный массив A

	integer, intent(in)			:: N
	! complex(8), dimension(N), intent(inout)		:: A
	! complex(8), dimension((N+1)/2), intent(out)	:: T

	complex(8)					:: V
	integer 					:: NA, NB

	complex(8), intent(inout)	:: A(N)
	complex(8), intent(out)		:: T((N+1)/2)

	if (N < 2) return

	if (N == 2) then
		if (real(A(1)) > real(A(2))) then
			V = A(1)
			A(1) = A(2)
			A(2) = V

		endif
		return
	endif      
	NA = (N+1) / 2
	NB = N - NA

	call MergeSort8(A, NA, T)
	call MergeSort8(A(NA+1), NB, T)

	if (real(A(NA)) > real(A(NA+1))) then
		T(1:NA) = A(1:NA)
		call Merge8(T, NA, A(NA+1), NB, A, N)
	endif

end subroutine MergeSort8


subroutine Merge8(A, NA, B, NB, C, NC)
	! Вспомогательная подпрограмма для MergeSort8,
	! которая что-то там переставляет...

	integer, intent(in)			:: NA, NB, NC       ! Normal usage: NA+NB = NC
	complex(8), intent(in)		:: B(NB)

	integer						:: I, J, K

	complex(8), intent(inout)	:: A(NA)        	! B overlays C(NA+1:NC)
	complex(8), intent(inout)	:: C(NC)

	I = 1
	J = 1
	K = 1

	do while(I <= NA .and. J <= NB)
		if (real(A(I)) <= real(B(J))) then
			C(K) = A(I)
			I = I + 1
		else
			C(K) = B(J)
			J = J + 1
		endif
		K = K + 1
	enddo

	do while (I <= NA)
		C(K) = A(I)
		I = I + 1
		K = K + 1
	enddo

end subroutine Merge8


end module utils