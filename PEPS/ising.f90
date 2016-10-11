
include "tensor.f90"

! H=J\sum_{\langle i,j\rangle}\sigma_{i}^{z}\sigma_{j}^{z}-g\sum_{i}\sigma_{i}^{x}
! J=1, g=0.5: E_exact = -1.063544409973372
!
! ld(1)       A(0)      ld(0)       A(1)      ld(1)
!l1---r1    j0---k0    l0---r0    j1---k1    l1---r1
!              |                     |
!              d0                    d1
!			             H
!					  u0  u1
!					   \  /
!					   /  \
!					  d0  d1

program main
	use omp_lib
	use mkl_service
	use M_tensor
	use M_matrix
	implicit none
	integer, parameter :: chi=6,d=2,N=70000
	integer :: i0,i1,i
	real(8) :: DJ=1d0,g=0.5d0,delta=0.005d0
	type(t_tensor) :: A(0:1),tmp1,tmp2
	real(8), allocatable :: s(:),U(:,:),V(:,:),H(:,:)
	real(8) :: ld(0:1,chi),HU(d**2,d**2),He(d**2),Hexp(d**2,d**2),E

	call omp_set_nested(.false.)
	!call omp_set_max_active_levels(1)
	call omp_set_dynamic(.false.)
	call mkl_set_dynamic(0)
	call omp_set_schedule(omp_sched_static,0)

	call mkl_set_num_threads(1)
	call omp_set_num_threads(mkl_get_max_threads())

	do i=0,1
		call A(i)%new([character(5)::"d"//to_char(i),"l","r"],[d,chi,chi],"r")
	enddo
	call random_number(ld)
	

	! ising
	call tmp1%new(["u0","d0","x"],[d,d,3],"0")
	call tmp2%new(["u1","d1","x"],[d,d,3],"0")
	tmp1%T(tmp1%get_idx([1,1,1, 2,2,1, 1,1,2, 2,2,2, 1,2,3, 2,1,3]))=[1d0,1d0,DJ,-DJ,-g*0.5d0,-g*0.5d0]
	tmp2%T(tmp2%get_idx([1,2,1, 2,1,1, 1,1,2, 2,2,2, 1,1,3, 2,2,3]))=[-g*0.5d0,-g*0.5d0,1d0,-1d0,1d0,1d0]

	tmp1=dot(tmp1,tmp2,["x","x"])
	call tmp1%get_mat(["u0","u1"],["d0","d1"],H)
	HU=H
	call syev(HU,He,"V")

	do i=0,N-1
		i0=mod(i,2)
		i1=mod(i+1,2)

		delta=0.005d0*exp(-19d0*i/N)
		Hexp=matmul(matmul(HU,diag(exp(-He*delta))),transpose(HU))

		tmp2=dot(ld(i1,:),A(i0),["l"])
		tmp2=dot(tmp2,ld(i0,:),["r"])
		tmp2=dot(tmp2,A(i1),["r","l"])
		tmp2=dot(tmp2,ld(i1,:),["r"])

		tmp1=dot(tmp2,Hexp,["d0","d1"])

		call tmp1%svd([character(5)::"l","d"//to_char(i0)],[character::],U,s,V)

		ld(i0,:)=s(:chi)/sqrt(sum(s(:chi)**2))

		call A(i0)%get_tensor([character(5)::tmp1%label(:2),"r"],U(:,:chi))
		call A(i1)%get_tensor([character(5)::"l",tmp1%label(3:)],V(:chi,:))

		ld(i1,:)=1d0/ld(i1,:)
		A(i0)=dot(ld(i1,:),A(i0),["l"])
		A(i1)=dot(A(i1),ld(i1,:),["r"])
		ld(i1,:)=1d0/ld(i1,:)



		if(mod(i,N/1000)==0.or.i==N-1) then
			E=0d0
			do i0=0,1
				i1=mod(i0+1,2)
				tmp2=dot(ld(i1,:),A(i0),["l"])
				tmp2=dot(tmp2,ld(i0,:),["r"])
				tmp2=dot(tmp2,A(i1),["r","l"])
				tmp2=dot(tmp2,ld(i1,:),["r"])
				tmp2=dot(dot(tmp2,H,["d0","d1"]),tmp2,["r","r","l","l","d0","d0","d1","d1"])
				E=E+tmp2%T(1)
			enddo
			write(*,*)E/2d0
		endif
	enddo
	!write(*,*)-log(sum(tmp1%T**2))/delta/2

end program
