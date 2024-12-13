program prog
use regime
use mod_eigenvalues
use matrix_fun
use bouwkamp
use ang_func
use utils
implicit none

! real(knd)					:: res_arr(5), mtr(5,2)
integer						:: NN, i, j, m, n, n0, lwork, info, ldvl, ldvr, p, s, eta_num, m0, m_step, m_num, i_m
integer			:: ii
real(knd)					:: rr, nu, eta0, eta_step
! complex(knd), allocatable :: mtr_c(:,:), mtr_c2(:,:), mtr_c3(:,:)
! real(knd), allocatable :: mtr_c4(:,:)
complex(knd)				:: c, f, dl, c1, c2(5), c3(4,5)
complex(knd), allocatable	:: mtr_mini(:,:), eig_arr(:), mtr(:,:), mtr2(:,:), work(:), vl(:,:), vr(:,:),&
 mtr_mini_all(:,:), f_arr(:), f_arr0, fff(:), Smn_eta_n(:,:)
! complex(8), allocatable		:: mtr_mini_8(:,:), eig_arr_8(:), mtr_8(:,:), mtr_mini_all_8(:,:)
real(knd), allocatable		:: rwork(:), eta_arr(:)
character(4)				:: norm
! complex(knd)				:: f_arr(-1:(p-(n-m))/2)


print*, "##############################################################"
! c = cmplx(5.45561817985860741941905871499329805_knd, 0_knd) 
! c = cmplx(8.18342726978791112912858807248994708_knd,5.45561817985860741941905871499329805_knd)
! c = cmplx(5_knd, 0_knd)
! c = cmplx(7.49999999999998667732370449812151492_knd,4.99999999999999111821580299874767661_knd)
! c = cmplx(0_knd,-4_knd)
! c = cmplx(0_knd,-5_knd)


! m = 0
! n = 6
! m = 1
! n = 4
! ! c = cmplx(1_knd,0_knd)
! c = cmplx(5_knd, 0_knd)
! ! nu = 0.99619469809175
! nu = 0.5
! nu = 0

open(1, file='../input.txt')
read(1,*) m0, m_step, m_num
read(1,*) n0, c
read(1,*) eta0, eta_step, eta_num
read(1,*) norm
close(1)


allocate(eta_arr(eta_num))
eta_arr = [real(knd) :: (eta0+i*eta_step, i=0, eta_num-1)]
! m = 0
! n = 6
! f = S_mn(c, 0._knd, m, 1, "Meix")
! c2 = S_mn_eta(c, [real(knd) :: 0, 0.25, 0.5, 0.75, 1], m, 1, "Meix")
! open(2, file='/home/sergey/kypca4/vb/build/output.txt')
open(2, file='../output.txt')
! open(3, file='../log.txt')
open(25, file='../eig.txt')

do i_m=0, m_num-1
	m = m0 + i_m*m_step
	n = n0 + m

	allocate(Smn_eta_n(m:n,eta_num), fff(eta_num))

	if (norm == "Flam") then
		Smn_eta_n = S_mn_eta_n(c, eta_arr, m, n)
	else
		Smn_eta_n = S_mn_eta_n(c, eta_arr, m, n, norm)
	end if

	fff =  S_mn_eta(c, eta_arr, m, n, norm)
	print*, "tttt", fff

	ii = 3
	f_arr0 =  S_mn(c, eta_arr(ii), m, n, norm)
	print*, "tttt1", f_arr0, eta_arr(ii)

	! f = S_mn(c, 0.0_knd, 0, 6)
	! print*, f
	! print*, "-------"
	! print*, c2
	! print*, "-------"
	! do i=1, size(Smn_eta_n(:,1))
	! 	print*,m+i-1, Smn_eta_n(i, :)

	! end do



	write(2, 50) c, m

	do i=m, n
		write(2, 150) i
		do j=1, eta_num
			write(2, 740) eta_arr(j)
			write(2, 780) Smn_eta_n(i, j)
			! write(2, *) Smn_eta_n(i, j)
			
		end do

	end do
	50	format(1x/ /'c = ',e39.30, e39.30,'; m = ',i5)
	150	format(1x,'l = ',i6)
	740	format(1x,'eta = ',f17.14$)

	780 format(6x,'s1 = ',f40.30,1x,f40.30)
	! c = cmplx(0_knd, -5_knd)
	deallocate(Smn_eta_n)
end do
close(2)
! close(3)
deallocate(eta_arr)


! s = mod(n-m, 2)

! NN = 100
! p = NN
! if (s /= mod(p, 2)) p = p + 1
! print*, "p", p
! c3 = N_mn(m, n, c, p, NN)
! print*, "-------"
! print*, "-------"

! print*, c3
print*, "-------"
print*, "m", m
print*, "n", n
print*, "c", c
print*, "nu", nu
close(25)



! do i=0, 12
! c1 = S_mn(c, 0.2_knd, 0, i, norm)
! print*, i, c1
! end do
	! c1 = S_mn(c, 0.2_knd, 0, 13, norm)
	! print*, c1
! m = 5
! rr = 1q0
! 	do i=2*m-1, 1, -2
! 		rr = rr * i		
! 	end do
! print*, rr

end program prog




! 1                                               f
! 1                                               number_of_layers
! 5.0q0                                           xv
! 1.5247025799298517701583439572615466q0                                           ab
! (1q0, 0q0) (1.5q0, 1q0)                      ri

! m = 1
! n = 1
! c = cmplx(7.49999999999998667732370449812151492_knd,4.99999999999999111821580299874767661_knd)

! 7.780629822129808:|9,                   4.9732550381294:|898 me8  15+13i
! 7.780629822129808|24468046264379661:645,4.9732550381294|9092451884650122293:889 vb8
! 7.780629822129808:2446804626437|7309238,4.9732550381294:909245188465012|0624493 me16  28+28i знаков после запятой
! 7.7806298221298082446804626437|9661:414,4.9732550381294909245188465012|2293:581 vb16

! ---------------------------------------------------------------------------------

! -1                                               f
! 1                                               number_of_layers
! 5.0q0                                           xv
! 1.7742319565673444813391511263099460q0                                           ab
! (1q0, 0q0) (1.5q0, 1q0)                      ri

! m = 2
! n = 12
! c = cmplx(0_knd,+-5_knd) // cmplx(5_knd,0_knd) vb


! 143.91073470297805526097739953411|3956,0.00000000000000000000000000000000000 me16
! 143.91073470297805526097739953411|4868,0.00000000000000000000000000000000000 vb16   29

! ---------------------------------------------------------------------------------

! -1                                               f
! 1                                               number_of_layers
! 5.0q0                                           xv
! 1.7742319565673444813391511263099460q0                                           ab
! (1q0, 0q0) (1.5q0, 1q0)                      ri

! m = 2
! n = 22
! c = cmplx(5_knd,-7.5_knd) // cmplx(7.5_knd,5_knd) vb

! 490.2151384069147273213052833105|56063,-36.945247445599196970088265220|9442628
! 490.2151384069147273213052833105|42208,-36.945247445599196970088265220|9678301 me16  28+27i
! 490.2151384069147273213052833105|22733,-36.945247445599196970088265220|8812464 vb16

! n=5
! 17.33765872373103296|40136031970200367,-21.5178881169255735731783479158114821
! 17.3376587236741962818413144265248862,-21.5178881169897663585170858244847161
! 17.3376589753963345688885856372827035,-21.5178886130079931955678238330329066
! 17.3385677818205751052889427212161349,-21.5179180378645792332232116192687257
! 17.33765872373103296|39957361490323956,-21.5178881169255735732497875764378877



! 8.33628222 4433673:62802182039884344362,-7.29892754044188595134828884506143600
! 8.33628222|344754566557599859797540958,-7.29892754688861824427202959783408781	
! 8.33628222|4433673:72082783361002577930,-7.29892754044188599031865847822985097