include "../lib/utility.f90"
module M_DMFT
	use M_const
	use M_matrix
	use M_utility
	use mkl_service
	implicit none
	real(8), parameter :: onepi=1d0/pi,halfpi=onepi*0.5d0
	integer, parameter :: nlog=1,nl=500,nadd(1)=[0],n=(nl*nlog+sum(nadd*2))*2,nk=128,D=2
	real(8) :: U=0d0,t=1d0,bcut=10d0,base=10d0/1d-3,iomega(n),omega(n),Tk,et=1d-2
	integer, parameter :: use_k=1,use_uk=2
	type gfs(n,D,nk)
		integer, len :: n,D,nk
		real(8) :: k(D,(2*nk)**D)
		real(8), allocatable :: uk(:,:)
		real(8) :: rhoc(n)
		complex(8) :: Gc(n),Delta(n),SEc(n)
		complex(8) :: Susq(n,(2*nk)**D),Susi(n,(2*nk)**D),Gck(n,(2*nk)**D),Gci(n,(2*nk)**D)
		logical :: conv
		integer :: flag
	contains
		procedure :: self_consistent
		procedure :: impurity_solver
		procedure :: export
		procedure, nopass :: ek
	end type
contains
	include "../largeN/mbroyden.f90"
	real(8) function ek(k) result(rt)
		real(8) :: k(D)
		integer :: i
		rt=0d0
		do i=1,size(k)
			rt=rt+2d0*t*cos(k(i))
		enddo
	end function
	subroutine export(self,ut)
		class(gfs(*,*,*)) :: self
		integer :: ut
		integer :: i,k
		write(ut,"(A)")"Tk omega rGc iGc rDelta iDelta rSEc iSEc"
		do i=1,n
			write(ut,"(*(e28.20))")Tk,omega(i),self%Gc(i),self%Delta(i),self%SEc(i)
		enddo
		write(ut,"(x/)")
	end subroutine
	subroutine impurity_solver(self)
		class(gfs(*,*,*)) :: self
		self%Gc=0.5d0/(omega-self%Delta+U/2d0)+0.5d0/(omega-self%Delta-U/2d0)
	end subroutine
	subroutine self_consistent(self,rate,niter,tol)
		class(gfs(*,*,*)) :: self
		real(8) :: rate,tol
		integer :: niter
		integer :: i,j,ik
		real(8) :: x(self%n),x_(size(x)),ek_((2*nk)**D)
		complex(8) :: tmp
		call mbroyden(0,x,x_,rate,40)
		!$OMP PARALLEL DO
		!do ik=1,(2*nk)**D
		!ek_(ik)=ek(self%k(1:D,ik))
		!enddo
		do ik=1,size(self%uk,2)
			ek_(ik)=ek(self%uk(1:D,ik))
		enddo
		!$OMP END PARALLEL DO
		do i=1,niter
			call set_realpart(omega,self%Delta)
			call self%impurity_solver()
			x=self%Delta%im
			self%SEc=omega-self%Delta-1d0/self%Gc
			if(D==0) then
			else
				!$OMP PARALLEL DO
				do j=1,n
					self%Gc(j)=0d0
					!do ik=1,(2*nk)**D
						!self%Gc(j)=self%Gc(j)+(1d0/(omega(j)-ek_(ik)-self%SEc(j)))
					!enddo
					do ik=1,size(self%uk,2)
						self%Gc(j)=self%Gc(j)+self%uk(0,ik)/(omega(j)-ek_(ik)-self%SEc(j)+img*et)
					enddo
				enddo
				!$OMP END PARALLEL DO
			endif
			self%Gc=self%Gc/real((2*nk)**D,kind=8)
			x_=imag(omega-1d0/self%Gc-self%SEc)
			x_=merge(x_,0d0,x_<0d0)
			write(*,*)i,maxval(abs(x-x_))
			if(all(abs(x-x_)<tol).or.i==niter+1) then
				self%conv=(i/=niter+1)
				exit
			else
				call mbroyden(i,x,x_)
			endif
			self%Delta%im=x
		enddo
		call mbroyden(huge(1),x,x_)
	end subroutine
	subroutine get_convolution(omega,A,B,tau,eta,rt)
		real(8) :: omega(:)
		complex(8) :: rt(:)
		complex(8) :: A(:),B(:)
		integer :: tau(2),eta(2)
		integer :: i,j,k,n,iwb(size(omega)),ibw(size(omega)),m
		real(8) :: e1,e2,omg
		complex(8) :: A1,A2,B1,B2,rt_
		procedure(real(8)), pointer :: f1,f2
		n=size(omega)
		select case(eta(1))
		case(1)
			f1 => fb
		case(-1)
			f1 => ff
		end select
		select case(eta(2))
		case(1)
			f2 => fb
		case(-1)
			f2 => ff
		end select
		!$OMP PARALLEL DO PRIVATE(omg,ibw,iwb,A1,B1,A2,B2,e1,e2,m,rt_)
		!do k=1,n/2
			do k=1,n
			omg=omega(k)
			do i=1,n
				ibw(i)=find_sidx(omega,tau(1)*(omg-tau(2)*omega(i)))
			end do
			A1=A(1)
			do i=1,n
				iwb(i)=find_sidx(omega,tau(2)*(omg-tau(1)*omega(i)))
			end do
			i=iwb(1)
			B1=merge(B(i)+(B(i+1)-B(i))/(omega(i+1)-omega(i))*(tau(2)*(omg-tau(1)*omega(1))-omega(i)),cmplx(0d0,kind=8),i>0.and.i<n)
			e1=omega(1)
			rt_=0d0
			do i=1,n-1
				do j=iwb(i),iwb(i+1),sign(1,iwb(i+1)-iwb(i))
					if(j==iwb(i+1)) then
						e2=omega(i+1)
						A2=A(i+1)
						B2=merge(B(j)+(B(j+1)-B(j))/(omega(j+1)-omega(j))*(tau(2)*(omg-tau(1)*e2)-omega(j)),cmplx(0d0,kind=8),j>0.and.j<n)
					else
						m=j+(1-tau(1)*tau(2))/2
						e2=tau(1)*(omg-tau(2)*omega(m))
						B2=B(m)
						m=ibw(m)
						A2=A(m)+(A(m+1)-A(m))/(omega(m+1)-omega(m))*(e2-omega(m))
					endif
					rt_=rt_+0.5d0*(e2-e1)*cmplx(&
						(-tau(1)*f1(beta*tau(1)*e1)*imag(A1)*real(B1)+tau(2)*eta(1)*eta(2)*f2(beta*(tau(1)*e1-omg))*real(A1)*imag(B1))+&
						(-tau(1)*f1(beta*tau(1)*e2)*imag(A2)*real(B2)+tau(2)*eta(1)*eta(2)*f2(beta*(tau(1)*e2-omg))*real(A2)*imag(B2))&
						,&
						(-tau(1)*tau(2)*f1(beta*tau(1)*e1)*imag(A1)*imag(B1)+tau(1)*tau(2)*eta(1)*eta(2)*f2(beta*(tau(1)*e1-omg))*imag(A1)*imag(B1))+&
						(-tau(1)*tau(2)*f1(beta*tau(1)*e2)*imag(A2)*imag(B2)+tau(1)*tau(2)*eta(1)*eta(2)*f2(beta*(tau(1)*e2-omg))*imag(A2)*imag(B2)),kind=8)
					A1=A2
					B1=B2
					e1=e2
				enddo
			enddo
			rt(k)=rt(k)-eta(1)*rt_/pi
			!rt(n+1-k)=rt(n+1-k)+cmplx(-eta(2)*rt_%re,eta(2)*rt_%im,kind=8)/pi
		enddo
		!$OMP END PARALLEL DO
	end subroutine
	subroutine set_realpart(omega,G)
		real(8) :: omega(:)
		complex(8) :: G(:)
		real(8) ::dx
		integer :: i,j,dn(2)
		G(1)%re=0d0
		G(n)%re=0d0
		do i=2,n-1
			G(i)%re=2d0*G(i)%im*log((omega(n)-omega(i))/(omega(i)-omega(1)))
		enddo
		do j=1,n
			dn=[min(n,j+1),max(1,j-1)]
			dx=omega(dn(1))-omega(dn(2))
			G(j)%re=G(j)%re+G(dn(1))%im-G(dn(2))%im
			do i=1,n
				if(j/=i) then
					G(i)%re=G(i)%re+(G(j)%im-G(i)%im)/(omega(j)-omega(i))*dx
				endif
			enddo
		enddo
		do i=1,n
			G(i)%re=G(i)%re*halfpi
		enddo
	end subroutine
	subroutine set_omega(omega,base,nlin,width)
		real(8) :: omega(:),base,width
		integer :: nlin
		integer :: nlog,n,i,j
		real(8) :: u,l,t,tmp
		n=size(omega)/2
		nlog=n/nlin
		if(mod(n,nlin)/=0) stop "check n and nlin"
		u=width
		l=u/base
		do i=nlog,1,-1
			t=(u-l)/real(nlin,8)
			tmp=l
			do j=1,nlin
				tmp=tmp+t
				omega(n+(i-1)*nlin+j)=tmp
				omega(n-((i-1)*nlin+j)+1)=-tmp
			enddo
			u=l
			l=u/base
		enddo
	end subroutine
	subroutine add_omega(omega,center,width,nadd)
		real(8) :: omega(:),center,width
		integer :: nadd,n,m,l,u,i
		real(8) :: s,at1,at2
		n=size(omega)-nadd*4
		l=find_sidx(omega(1:n),center-width)
		u=find_sidx(omega(1:n),center+2d0*width)
		omega(u+1+nadd*4:)=omega(u+1:n)
		omega(n/2+1+nadd*2:l+nadd*2)=omega(n/2+1:l)
		omega(n-l+1+nadd*2:n/2+nadd*2)=omega(n-l+1:n/2)
		n=size(omega)
		l=l+nadd*2
		u=u+nadd*4
		m=u-l
		m=l+m-m/2
		s=width/100d0
		at1=atan(width)
		at2=atan(width)
		do i=l+1,u
			if(i<=m) then
				omega(i)=center-tan(s+(at1-s)/(m-l)*(m-i))
			else
				omega(i)=center+tan(s+(at2-s)/(u-m)*(i-m-1))
			endif
			omega(n-i+1)=-omega(i)
		enddo
	end subroutine
	integer function find_sidx(omega,x) result(l)
		real(8) :: omega(:),x
		integer :: m,u
		l=0
		u=size(omega)+1
		do
			if((u-l)>1) then
				m=(l+u)/2
				if(x>omega(m)) then
					l=m
				else
					u=m
				endif
			else
				exit
			endif
		enddo
	end function
	subroutine get_uk(k,uk)
		real(8) :: k(:,:)
		real(8), allocatable :: uk(:,:)
		real(8) :: tmp(0:size(k,1),size(k,2))
		integer :: n,i,j
		n=1
		tmp(1:D,n)=k(:,1)
		call into_sym(tmp(1:D,n))
		tmp(0,n)=1d0
		do i=2,size(k,2)
			tmp(1:D,n+1)=k(:,i)
			call into_sym(tmp(1:D,n+1))
			do j=1,n
				if(sum(abs(tmp(1:D,n+1)-tmp(1:D,j)))<1d-6) then
					tmp(0,j)=tmp(0,j)+1d0
					exit
				endif
			enddo
			if(j==n+1) then
				n=n+1
				tmp(0,n)=1d0
			endif
		enddo
		allocate(uk(0:D,n))
		uk=tmp(0:D,1:n)
		write(*,*)nint(sum(uk(0,:)))==(2*nk)**D
	end subroutine
	subroutine into_sym(k)
		real(8) :: k(:)
		integer :: i,j
		real(8) :: tmp
		k=abs(k)
		do i=1,D
			do j=i+1,D
				if(k(i)<k(j)) then
					tmp=k(i)
					k(i)=k(j)
					k(j)=tmp
				endif
			enddo
		enddo
	end subroutine
	subroutine get_uk_sort(k,uk)
		real(8) :: k(:,:)
		real(8), allocatable :: uk(:,:)
		real(8) :: tmp(1:size(k,1),size(k,2))
		integer :: n,i,j
		integer :: ord(size(k,2))
		integer, allocatable :: c(:)
		do i=1,size(k,2)
			tmp(1:,i)=k(:,i)
			call into_sym(tmp(1:,i))
		enddo
		ord=[1:size(k,2)]
		call qsort(tmp(:,:),ord)
		!do i=1,size(tmp,2)
			!write(*,"(*(es12.4))")tmp(:,i),tmp(:,ord(i))
		!enddo
		tmp(:,:)=tmp(:,ord)
		call collect(tmp,c)
		!write(*,*)c
		!stop
		allocate(uk(0:D,size(c)-1))
		uk(1:D,:)=tmp(:,c(1:size(c)-1))
		uk(0,:)=real(c(2:)-c(1:size(c)-1),kind=8)
		write(*,*)nint(sum(uk(0,:)))==(2*nk)**D
	end subroutine
end module
include "../lib/serde.f90"
program main
	use M_DMFT
	use M_serde
	implicit none
	type(gfs(n,D,nk)) :: phy
	integer :: j,i,ik
	type(t_serde(D)) :: map

	call omp_set_nested(.false.)
	call omp_set_dynamic(.false.)
	call mkl_set_dynamic(0)

	call mkl_set_num_threads(16)
	call omp_set_num_threads(mkl_get_max_threads())

	open(unit=10,File="../data/fk.dat")
	open(unit=20,file="../data/fk_init.dat")
	call set_omega(omega(1:nl*nlog*2),base,nl,bcut)
	!call add_omega(omega(1:nl*nlog*2+sum(nadd(1:1))*4),2.4d0,0.3d0,nadd(1))
	!call add_omega(omega(1:nl*nlog*2+sum(nadd(1:2))*4),4d0,0.5d0,nadd(2))
	map%shap=[(2*nk,i=1,D)]
	!$OMP PARALLEL DO
	do ik=1,size(phy%k,2)
		phy%k(:,ik)=pi*map%get_idx(ik)/nk-pi
	enddo
	!$OMP END PARALLEL DO
	write(*,*)"stating get unique k point...."
	!call get_uk(phy%k,phy%uk)
	call get_uk_sort(phy%k,phy%uk)
	write(*,*)"finished get unique k point",size(phy%k,2)/size(phy%uk,2)
	!do i=1,size(phy%uk,2)
		!write(*,*)phy%uk(:,i)
		!read(*,*)
	!enddo
	!do i=1,n
		!read(20,"(2(e28.20))")phy%Delta(i)
	!enddo
	phy%Delta%im=-exp(abs(omega))
	phy%flag=use_uk
	call phy%self_consistent(0.3d0,3000,1d-5)
	if(.not.phy%conv) then
		write(*,*)"not converge!!!!"
	else
		rewind(20)
		do i=1,n
			write(20,"(*(e28.20))")phy%Delta(i)
		enddo
	endif

	!$OMP PARALLEL DO
	do j=1,n
		do ik=1,size(phy%k,2)
			phy%Gck(j,ik)=1d0/(omega(j)-ek(phy%k(1:D,ik))-phy%SEc(j)+img*et)
		enddo
		call mfft(phy%Gck(j,:),1,map%shap,phy%Gci(j,:))
	enddo
	!$OMP END PARALLEL DO
	do i=1,size(phy%k,2)
		call get_convolution(omega,phy%Gci(:,i),phy%Gci(:,i),[1,-1],[-1,-1],phy%Susi(:,i))
	enddo
	do j=1,n
		call mfft(phy%Susi(j,:),-1,map%shap,phy%Susq(j,:))
	enddo
	call phy%export(10)
end program
