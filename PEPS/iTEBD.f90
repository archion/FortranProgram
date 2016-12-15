include "tensor.f90"
!include "tensor.f90"

! H=J\sum_{\langle i,j\rangle}\sigma_{i}^{z}\sigma_{j}^{z}-g\sum_{i}\sigma_{i}^{x}
! J=1, g=0.5: E_exact = -1.063544409973372
!
! ld(1)       A(0)      ld(0)       A(1)      ld(1)
!  -.-       l-.-r       -.-       l-.-r       -.-  
!              |                     |
!              d0                    d1
!			             H
!					  u0  u1
!					   \  /
!					   /  \
!					  d0  d1
! E_h=-1.772588722
! ising: N=70000, alpha=19, delta=0.005, chi=6
! heisenberg: N=20000, alpha=10, delta=0.05, chi=24

program main
	use omp_lib
	use mkl_service
	use M_tensor
	use M_matrix
	implicit none
	integer, parameter :: d=2,chi=24
	integer :: i,j,i0,i1,N=20000
	real(8) :: DJ=1d0,g=0.5d0,delta=0.05d0,alpha=10d0
	type(t_tensor) :: A(0:1),App,tmp1,tmp2,check
	real(8), allocatable :: s(:)
	complex(8), allocatable :: U(:,:),V(:,:),X(:,:),iX(:,:),Y(:,:),iY(:,:),H(:,:)
	complex(8) :: vec(chi*chi),ld(0:1,chi),Hexp(d**2,d**2),HU(size(Hexp,1),size(Hexp,2)),E
	real(8) :: He(size(Hexp,1))
	complex(8) :: w

	call omp_set_nested(.false.)
	!call omp_set_max_active_levels(1)
	call omp_set_dynamic(.false.)
	call mkl_set_dynamic(0)
	call omp_set_schedule(omp_sched_static,0)

	call mkl_set_num_threads(1)
	call omp_set_num_threads(mkl_get_max_threads())

	call A(0)%new(["d0","l","r"],[d,chi,chi],flag="r")
	call A(1)%new(["d1","l","r"],[d,chi,chi],flag="r")

	call random_number(ld(0,:))
	call random_number(ld(1,:))
	ld=real(ld)

	!! ising
	!call tmp1%new(["u0","d0","x"],[d,d,3],"0")
	!call tmp2%new(["u1","d1","x"],[d,d,3],"0")
	!tmp1%rc%T(tmp1%get_idx([1,1,1, 2,2,1, 1,1,2, 2,2,2, 1,2,3, 2,1,3]))=[1d0,1d0,DJ,-DJ,-g*0.5d0,-g*0.5d0]
	!tmp2%rc%T(tmp2%get_idx([1,2,1, 2,1,1, 1,1,2, 2,2,2, 1,1,3, 2,2,3]))=[-g*0.5d0,-g*0.5d0,1d0,-1d0,1d0,1d0]

	! Heisenberg
	call tmp1%new(["u0","d0","x"],[d,d,3],flag="0")
	call tmp2%new(["u1","d1","x"],[d,d,3],flag="0")
	tmp1%rc%T(tmp1%get_idx([1,1,1, 2,2,1, 1,2,2, 2,1,3]))=[1d0,-1d0,1d0,1d0]
	tmp2%rc%T(tmp2%get_idx([1,1,1, 2,2,1, 2,1,2, 1,2,3]))=[1d0,-1d0,2d0,2d0]

	tmp1=dot(tmp1,tmp2,["x","x"])
	call tmp1%get_mat(["u0","u1"],["d0","d1"],H)
	HU=H
	!call syev(HU,He,"V")
	call heev(HU,He,"V")

	j=0
	do i=0,N-1
		i0=mod(i,2)
		i1=mod(i+1,2)
		

		delta=0.05d0*exp(-alpha*i/N)
		Hexp=matmul(matmul(HU,diag(exp(-He*delta))),transpose(conjg(HU)))

		App=dot(dot(dot(A(i0),ld(i0,:),["r"]),A(i1),["r","l"]),Hexp,["d0","d1"])

		tmp1=dot(App,ld(i1,:),["r"])
		tmp1=dot(set_conjg(tmp1%change_label(["l","r"],["L","R"])),tmp1,["d0","d0","d1","d1"])
		call tmp1%get_mat(["l","L"],["r","R"],U)
		call ar_naupd(U,w,vec,lr="r")
		vec=vec*conjg(vec(1))/abs(conjg(vec(1)))
		call allocate(X,reshape(vec,[chi,chi]))
		call allocate(s,chi)
		!call syev(X,s,"V")
		call heev(X,s,"V")
		s=abs(s)
		call allocate(iX,matmul(diag(1d0/sqrt(s)),transpose(conjg(X))))
		X=matmul(X,diag(sqrt(s)))

		tmp1=dot(ld(i1,:),App,["l"])
		tmp1=dot(set_conjg(tmp1%change_label(["l","r"],["L","R"])),tmp1,["d0","d0","d1","d1"])
		call tmp1%get_mat(["l","L"],["r","R"],U)
		call ar_naupd(U,w,vec,lr="l")
		vec=vec*conjg(vec(1))/abs(conjg(vec(1)))
		call allocate(Y,reshape(vec,[chi,chi]))
		!call syev(Y,s,"V")
		call heev(Y,s,"V")
		s=abs(s)
		call allocate(iY,transpose(matmul(diag(1d0/sqrt(s)),transpose(conjg(Y)))))
		Y=transpose(matmul(Y,diag(sqrt(s))))

		call svd(matmul(matmul(Y,diag(ld(i1,:))),X),U,s,V)
		ld(i1,:)=s/sqrt(sum(s**2))

		App=dot(dot(ld(i1,:),dot(dot(matmul(V,iX),App,["l"]),matmul(iY,U),["r"]),["l"]),ld(i1,:),["r"])
		call App%svd([character(5)::"l","d"//to_char(i0)],[character::],U,s,V)

		ld(i0,:)=s(:chi)/sqrt(sum(s(:chi)**2))

		call A(i0)%get_tensor([character(5)::App%label(:2),"r"],U(:,:chi))
		call A(i1)%get_tensor([character(5)::"l",App%label(3:)],V(:chi,:))

		A(i0)=dot(1d0/ld(i1,:),A(i0),["l"])
		A(i1)=dot(A(i1),1d0/ld(i1,:),["r"])

		if(mod(i,N/1000)==0.or.i==N-1) then
			j=j+1
			E=0d0
			do i0=0,1
				i1=mod(i0+1,2)
				App=dot(dot(dot(dot(ld(i1,:),A(i0),["l"]),ld(i0,:),["r"]),A(i1),["r","l"]),ld(i1,:),["r"])
				App=dot(dot(App,H,["d0","d1"]),set_conjg(App),["d0","d0","d1","d1","l","l","r","r"])
				E=E+App%rc%T(1)
			enddo
			write(*,*)E/2d0
			!if(j==10) exit
		endif

	enddo
	write(*,*)N,alpha

end program
