module M_arpack
	use blas95
	implicit none
contains
	subroutine ar_naupd(H,E,V,lr)
		integer, parameter :: ncv = 20, nev = 1
		real(8) :: H(:,:),V(:)
		complex(8) :: E
		real(8) :: sigmar, sigmai, tol
		integer :: iparam(11), ipntr(14), info, ierr, ido ,n
		real(8) :: d(ncv,3), resid(size(H,1)), vec(size(H,1),ncv), workd(3*size(H,1)), workev(3*ncv), workl(3*ncv*ncv+6*ncv)
		logical :: select(ncv)
		integer :: i
		character(1), optional :: lr
		n=size(H,1)
		tol=0d0
		ido=0
		info=0
		iparam(1)=1
		iparam(3)=size(V)
		iparam(7)=1
		do
			call dnaupd ( ido, 'I', n, 'LM', nev, tol, resid, ncv, vec, n, iparam, ipntr, workd, workl, size(workl), info )
			if (ido==-1 .or. ido==1) then
				if(present(lr)) then
					if(lr=="l") then
						call gemv(H,workd(ipntr(1):ipntr(1)+n-1),workd(ipntr(2):ipntr(2)+n-1),trans="t")
						!workd(ipntr(2):ipntr(2)+n-1)=matmul(workd(ipntr(1):ipntr(1)+n-1),H)
						!!$omp parallel do
						!do i=0,n-1
							!workd(ipntr(2)+i)=sum(workd(ipntr(1):ipntr(1)+n-1)*H(:,i+1))
						!enddo
						!!$omp end parallel do
						cycle
					endif
				endif
				call gemv(H,workd(ipntr(1):ipntr(1)+n-1),workd(ipntr(2):ipntr(2)+n-1))
				!workd(ipntr(2):ipntr(2)+n-1)=matmul(H,workd(ipntr(1):ipntr(1)+n-1))
				!!$omp parallel do
				!do i=0,n-1
					!workd(ipntr(2)+i)=sum(H(i+1,:)*workd(ipntr(1):ipntr(1)+n-1))
				!enddo
				!!$omp end parallel do
			else 
				if(info==0) then
					call dneupd ( .true., 'A', select, d, d(1,2), vec, n, sigmar, sigmai, workev, 'I', n, 'LM', nev, tol, resid, ncv, vec, n, iparam, ipntr, workd, workl, size(workl), ierr )
					exit
				else
					stop "arpack info error"
				endif
			endif
		enddo
		E=cmplx(d(1,1),d(1,2))
		V=vec(:,1)
	end subroutine
end module
include "tensor.f90"

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
	use M_arpack
	implicit none
	integer, parameter :: chi=24,d=2
	integer :: i,j,i0,i1,N=20000
	real(8) :: DJ=1d0,g=0.5d0,delta=0.005d0,alpha=11d0
	type(t_tensor) :: A(0:1),App,tmp1,tmp2
	real(8), allocatable :: s(:),U(:,:),V(:,:),X(:,:),iX(:,:),Y(:,:),iY(:,:),H(:,:)
	real(8) :: vec(chi*chi),ld(0:1,chi),Hexp(d**2,d**2),HU(size(Hexp,1),size(Hexp,2)),He(size(Hexp,1)),E
	complex(8) :: w

	call omp_set_nested(.false.)
	!call omp_set_max_active_levels(1)
	call omp_set_dynamic(.false.)
	call mkl_set_dynamic(0)
	call omp_set_schedule(omp_sched_static,0)

	call mkl_set_num_threads(1)
	call omp_set_num_threads(mkl_get_max_threads())

	call A(0)%new(["d0","l","r"],[d,chi,chi],"r")
	call A(1)%new(["d1","l","r"],[d,chi,chi],"r")

	call random_number(ld(0,:))
	call random_number(ld(1,:))

	!allocate(Y(40,30))
	!call random_number(Y)
	!H=Y
	!call qr(Y,U,V)
	!write(*,*)sum(abs(matmul(U,V)-H))
	!write(*,*)sum(abs(matmul(U,transpose(U))-diag(1d0,size(U,1))))
	!stop
	
	!! ising
	!call tmp1%new(["u0","d0","x"],[d,d,3],"0")
	!call tmp2%new(["u1","d1","x"],[d,d,3],"0")
	!tmp1%T(tmp1%get_idx([1,1,1, 2,2,1, 1,1,2, 2,2,2, 1,2,3, 2,1,3]))=[1d0,1d0,DJ,-DJ,-g*0.5d0,-g*0.5d0]
	!tmp2%T(tmp2%get_idx([1,2,1, 2,1,1, 1,1,2, 2,2,2, 1,1,3, 2,2,3]))=[-g*0.5d0,-g*0.5d0,1d0,-1d0,1d0,1d0]

	! Heisenberg
	call tmp1%new(["u0","d0","x"],[d,d,3],"0")
	call tmp2%new(["u1","d1","x"],[d,d,3],"0")
	tmp1%T(tmp1%get_idx([1,1,1, 2,2,1, 1,2,2, 2,1,3]))=[1d0,-1d0,1d0,1d0]
	tmp2%T(tmp2%get_idx([1,1,1, 2,2,1, 2,1,2, 1,2,3]))=[1d0,-1d0,2d0,2d0]

	tmp1=dot(tmp1,tmp2,["x","x"])
	call tmp1%get_mat(["u0","u1"],["d0","d1"],H)
	HU=H
	call syev(HU,He,"V")

	call tmp1%clear()
	call tmp2%clear()


	j=0
	do i=0,N-1
		i0=mod(i,2)
		i1=mod(i+1,2)
		

		delta=0.05d0*exp(-alpha*i/N)
		Hexp=matmul(matmul(HU,diag(exp(-He*delta))),transpose(HU))

		App=dot(A(i0),ld(i0,:),["r"])
		App=dot(App,A(i1),["r","l"])
		App=dot(App,Hexp,["d0","d1"])

		tmp1=dot(App,ld(i1,:),["r"])
		tmp2=tmp1
		call tmp2%change_label(["l","L","r","R"])
		tmp2=dot(tmp1,tmp2,["d0","d0","d1","d1"])
		call tmp2%get_mat(["l","L"],["r","R"],U)
		call ar_naupd(U,w,vec,lr="r")
		call allocate(X,reshape(vec,[chi,chi]))
		if(allocated(s)) deallocate(s)
		allocate(s(chi))
		call syev(X,s,"V")
		s=abs(s)
		call allocate(iX,matmul(diag(1d0/sqrt(s)),transpose(X)))
		X=matmul(X,diag(sqrt(s)))

		tmp1=dot(ld(i1,:),App,["l"])
		tmp2=tmp1
		call tmp2%change_label(["l","L","r","R"])
		tmp2=dot(tmp1,tmp2,["d0","d0","d1","d1"])
		call tmp2%get_mat(["l","L"],["r","R"],U)
		call ar_naupd(U,w,vec,lr="l")
		call allocate(Y,reshape(vec,[chi,chi]))
		call syev(Y,s,"V")
		s=abs(s)
		call allocate(iY,transpose(matmul(diag(1d0/sqrt(s)),transpose(Y))))
		Y=transpose(matmul(Y,diag(sqrt(s))))

		call svd(matmul(matmul(Y,diag(ld(i1,:))),X),U,s,V)
		ld(i1,:)=s/sqrt(sum(s**2))


		App=dot(dot(matmul(V,iX),App,["l"]),matmul(iY,U),["r"])
		App=dot(dot(ld(i1,:),App,["l"]),ld(i1,:),["r"])
		call App%svd([character(5)::"l","d"//to_char(i0)],[character::],U,s,V)

		ld(i0,:)=s(:chi)/sqrt(sum(s(:chi)**2))

		call A(i0)%get_tensor([character(5)::App%label(:2),"r"],U(:,:chi))
		call A(i1)%get_tensor([character(5)::"l",App%label(3:)],V(:chi,:))

		ld(i1,:)=1d0/ld(i1,:)
		A(i0)=dot(ld(i1,:),A(i0),["l"])
		A(i1)=dot(A(i1),ld(i1,:),["r"])
		ld(i1,:)=1d0/ld(i1,:)

		!if(mod(i,100)==0) then
			!App=dot(ld(i1,:),A(i0),["l"])
			!App=dot(App,ld(i0,:),["r"])
			!App=dot(App,A(i1),["r","l"])
			!App=dot(App,ld(i1,:),["r"])
			!App=dot(dot(App,H,["d0","d1"]),App,["d0","d0","d1","d1","l","l","r","r"])
			!write(*,"(es18.10$)")App%T

			!App=dot(ld(i1,:),A(i0),["l"])
			!App=dot(App,ld(i0,:),["r"])
			!App=dot(App,A(i1),["r","l"])
			!App=dot(App,App,["d0","d0","d1","d1","l","l"])
			!write(*,*)sum(abs(App%T-reshape(diag(App%T(1),chi),shape(App%T))))

			!read(*,*)
		!endif
		if(mod(i,N/1000)==0.or.i==N-1) then
			j=j+1
			E=0d0
			do i0=0,1
				i1=mod(i0+1,2)
				App=dot(ld(i1,:),A(i0),["l"])
				App=dot(App,ld(i0,:),["r"])
				App=dot(App,A(i1),["r","l"])
				App=dot(App,ld(i1,:),["r"])
				App=dot(dot(App,H,["d0","d1"]),App,["d0","d0","d1","d1","l","l","r","r"])
				E=E+App%T(1)
			enddo
			write(*,*)E/2d0
			if(j==10) exit
		endif

		call App%clear()
		call tmp1%clear()
		call tmp2%clear()

	enddo
	write(*,*)N,alpha

end program
