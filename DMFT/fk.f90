include "../lib/utility.f90"
module M_DMFT
	use M_const
	use M_matrix
	use M_utility
	use mkl_service
	implicit none
	real(8), parameter :: onepi=1d0/pi,halfpi=onepi*0.5d0
	integer, parameter :: nadd(2)=[4000,1000],n=(sum(nadd))*2,nk=256/8,D=2
	real(8) :: U=8d0,t=1d0,bcut=200d0,iomega(n),omega(n),Tk,beta,et=1d-2,mu
	integer, parameter :: use_k=1,use_uk=2,export_w=1,export_k=2
	integer :: iq
	type gfs(n,D,nk)
		integer, len :: n,D,nk
		real(8) :: k(D,nk**D)
		real(8), allocatable :: uk(:,:)
		real(8) :: rhoc(n),B(nk**2)
		complex(8) :: Gc(n),Delta(n),SEc(n),Vertex(n,n)
		complex(8) :: Sus(nk**D),Sus0(n,nk**D),Gck(n,nk**D),Gci(n,nk**D),GM(n)
		logical :: conv
		integer :: flag
	contains
		procedure :: self_consistent
		procedure :: impurity_solver,get_vertex
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
	subroutine export(self,ut,flag)
		class(gfs(*,*,*)) :: self
		integer :: ut(:),flag(:)
		integer :: i,k,ik
		do k=1,size(flag)
			select case(flag(k))
			case(export_w)
				write(ut(k),"(A)")"Tk omega rGc iGc rDelta iDelta rSEc iSEc"
				do i=1,n
					write(ut(k),"(*(e28.20))")Tk,omega(i),self%Gc(i),self%Delta(i),self%SEc(i)
				enddo
				write(ut(k),"(x/)")
			case(export_k)
				write(ut(k),"(A)")"Tk k rSus iSus B"
				do ik=1,nk**D
					write(ut(k),"(*(e28.20))")Tk,self%k(:,ik),self%Sus(ik),self%B(ik)
				enddo
				write(ut(k),"(x/)")
			end select
		end do
	end subroutine
	subroutine impurity_solver(self)
		class(gfs(*,*,*)) :: self
		self%Gc=0.5d0/(omega-self%Delta+mu)+0.5d0/(omega-self%Delta+mu-U)
	end subroutine
	subroutine get_vertex(self)
		class(gfs(*,*,*)) :: self
		integer :: i,j
		self%GM=1d0/(self%Gc**2*(omega-self%Delta+mu)*(omega-self%Delta+mu-U))
		do i=1,n
			do j=1,n
				if(i==j) then
					self%Vertex(i,j)=0d0
				else
					!self%Vertex(i,j)=(self%Gc(i)*self%Gc(j))**-2*(0.5d0/((omega(i)-self%Delta(i)+mu)*(omega(j)-self%Delta(j)+mu))+0.5d0/((omega(i)-self%Delta(i)+mu-U)*(omega(j)-self%Delta(j)+mu-U))-self%Gc(i)*self%Gc(j))
					self%Vertex(i,j)=(0.5d0/((omega(i)-self%Delta(i)+mu)*(omega(j)-self%Delta(j)+mu))+0.5d0/((omega(i)-self%Delta(i)+mu-U)*(omega(j)-self%Delta(j)+mu-U))-self%Gc(i)*self%Gc(j))/((self%Gc(i)*self%Gc(j))**2)
					!write(*,*)self%Vertex(i,j),0.25d0*U**2*GM(i)*GM(j)
					!read(*,*)
				endif
			enddo
		enddo
	end subroutine
	subroutine self_consistent(self,rate,niter,tol)
		class(gfs(*,*,*)) :: self
		real(8) :: rate,tol
		integer :: niter
		integer :: i,j,ik
		real(8) :: x(self%n),x_(size(x)),ek_(nk**D)
		complex(8) :: tmp
		mu=0.5d0*U
		call mbroyden(0,x,x_,rate,40)
		!$OMP PARALLEL DO
		!do ik=1,nk**D
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
					!do ik=1,nk**D
						!self%Gc(j)=self%Gc(j)+(1d0/(omega(j)-ek_(ik)-self%SEc(j)))
					!enddo
					do ik=1,size(self%uk,2)
						self%Gc(j)=self%Gc(j)+self%uk(0,ik)/(omega(j)-ek_(ik)-self%SEc(j)+img*et)
					enddo
				enddo
				!$OMP END PARALLEL DO
			endif
			self%Gc=self%Gc/real(nk**D,kind=8)
			x_=imag(omega-1d0/self%Gc-self%SEc)
			x_=merge(x_,0d0,x_<0d0)
			write(*,*)i,maxval(abs(x-x_))
			if(all(abs(x-x_)<tol).or.i==niter+1) then
				self%conv=(i/=niter+1)
				exit
			else
				call mbroyden(i,x,x_)
			endif
			!do j=1,n
				!self%Delta%im(j)=x(j)
			!enddo
			self%Delta=cmplx(0d0,x,kind=8)
		enddo
		call mbroyden(huge(1),x,x_)
		call self%impurity_solver()
	end subroutine
	subroutine get_convolution(omega,A,B,tau,eta,rt,omega_rt)
		real(8) :: omega(:)
		complex(8) :: A(:),B(:)
		integer :: tau(2),eta(2)
		complex(8) :: rt(:)
		real(8), optional :: omega_rt(:)
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
			do k=1,size(rt)
			if(present(omega_rt)) then
				omg=omega_rt(k)
			else
				omg=omega(k)
			endif
			do i=1,n
				ibw(i)=find(omega,tau(1)*(omg-tau(2)*omega(i)))
			end do
			A1=A(1)
			do i=1,n
				iwb(i)=find(omega,tau(2)*(omg-tau(1)*omega(i)))
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
	subroutine into_sym(k)
		real(8) :: k(:)
		integer :: i,j
		real(8) :: tmp
		k=abs(k-pi)
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
	subroutine get_uk(k,uk)
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
		tmp(:,:)=tmp(:,ord)
		call collect(tmp,c)
		allocate(uk(0:D,size(c)-1))
		uk(1:D,:)=tmp(:,c(1:size(c)-1))
		uk(0,:)=real(c(2:)-c(1:size(c)-1),kind=8)
		write(*,*)nint(sum(uk(0,:)))==nk**D
	end subroutine
end module
include "../lib/serde.f90"
include "../lib/fft.f90"
program main
	use M_DMFT
	use M_serde
	use M_fft
	implicit none
	type(gfs(n,D,nk)) :: phy
	integer :: j,i,ik,jk,iT,iU
	real(8) :: ff_(n),dTk
	type(t_serde(D)) :: map

	call omp_set_nested(.false.)
	call omp_set_dynamic(.false.)
	call mkl_set_dynamic(0)

	call mkl_set_num_threads(1)
	!call omp_set_num_threads(mkl_get_max_threads())
	call omp_set_num_threads(4)

	open(unit=10,File="../data/fk.dat")
	open(unit=20,file="../data/fk_init.dat")
	open(unit=30,File="../data/fk_q.dat")
	open(unit=40,File="../data/fk_phase.dat")
	!call add_omega(omega(1:nl*nlog*2+sum(nadd(1:2))*4),4d0,0.5d0,nadd(2))
	map%shap=[(nk,i=1,D)]
	iq=0
	!!$OMP PARALLEL DO
	do ik=1,size(phy%k,2)
		phy%k(:,ik)=2d0*pi*(map%get_idx(ik,[1:D])-1)/nk
		if(sum(abs(phy%k(:,ik)-pi))<1d-8) then
			if(iq==0) then
				iq=ik
			else
				stop "err"
			endif
		endif
	enddo
	!!$OMP END PARALLEL DO
	write(*,*)"stating get unique k point...."
	call get_uk(phy%k,phy%uk)
	write(*,*)"finished get unique k point",size(phy%k,2)/size(phy%uk,2)
	!do i=1,size(phy%uk,2)
		!write(*,*)phy%uk(:,i)
		!read(*,*)
	!enddo
	!do i=1,n
		!read(20,"(2(e28.20))")phy%Delta(i)
	!enddo
	phy%flag=use_uk
	write(40,"(A)")"U T"
	do iU=1,1,1
		U=iU
		!call set_omega(omega(1:nl*nlog*2),base,nl,bcut)
		!call add_omega(omega(1:nl*nlog*2+sum(nadd(1:1))*4),max(U/2d0-5d0,0d0),min(U/2d0+5d0,bcut),nadd(1))
		call set_grid(omega(n/2+1:),[reset,add_linear,add_linear],from=[0d0,max(U/2d0-5d0,0d0)],to=[bcut,min(U/2d0+5d0,bcut)],n=[nadd])
		!write(*,*)omega(n/2+1:)
		!stop
		omega(1:n/2)=-omega(n:n/2+1:-1)
		!phy%Delta%im=-exp(abs(omega(i)))
		phy%Delta=cmplx(0d0,-exp(abs(omega)),kind=8)
		call phy%self_consistent(0.3d0,3000,1d-6)
		call phy%get_vertex()
		if(.not.phy%conv) then
			write(*,*)"not converge!!!!"
		else
			rewind(20)
			do i=1,n
				write(20,"(*(e28.20))")phy%Delta(i)
			enddo
		endif
		call phy%export([10,30],[export_w,export_k])
		!stop
		Tk=1d-4
		if(next(Tk,0d0,1d-1,1d-4)) then
		endif
		do
			if(U==0d0) then
				write(40,"(*(e28.20))")U,0d0
				exit
			endif
			beta=1d0/Tk

			do i=1,n
				ff_(i)=ff(beta*omega(i))
			enddo

			!$OMP PARALLEL DO
			do j=1,n
				do ik=1,size(phy%k,2)
					phy%Gck(j,ik)=1d0/(omega(j)-ek(phy%k(1:D,ik))-phy%SEc(j)+img*et)
				enddo
				call mfft(phy%Gck(j,:),1,map%shap,phy%Gci(j,:),1d0/product(map%shap))
			enddo
			!$OMP END PARALLEL DO
			!$OMP PARALLEL DO
			do i=1,size(phy%k,2)
				!call get_convolution(omega,phy%Gci(:,i),phy%Gci(:,i),[1,-1],[-1,-1],phy%Sus(i:i),[0d0])
				phy%Sus(i)=-integrate(omega,A=-imag(phy%Gci(:,i)**2)*ff_/pi)
			enddo
			!$OMP END PARALLEL DO
			call mfft(phy%Sus,-1,map%shap)
			!$OMP PARALLEL DO
			do i=1,n
				phy%Sus0(i,:)=-phy%Gci(i,:)**2
				call mfft(phy%Sus0(i,:),-1,map%shap)
			enddo
			!$OMP END PARALLEL DO
			i=iq
			!do i=1,size(phy%k,2)
				phy%Sus0(:,i)=phy%Sus0(:,i)+phy%Gc**2
				phy%B(i)=integrate(omega,A=-imag(0.25d0*U**2*phy%GM**2*phy%Sus0(:,i)/(1d0+0.25d0*U**2*phy%GM**2*phy%Sus0(:,i)))*ff_/pi)/Tk
			!enddo
			write(*,*)Tk,phy%B(iq)
			if(next(Tk,phy%B(iq)-1d0)) then
				write(40,"(*(e28.20))")U,Tk
				exit
			endif
		enddo
	enddo
	!do ik=1,size(phy%k,2)
		!write(*,*)phy%Sus(ik:ik)
		!!phy%Sus(ik:ik)=-phy%Sus(ik:ik)
		!phy%Sus(ik:ik)=0d0
		!do i=1,size(phy%k,2)
			!do j=1,size(phy%k,2)
				!if(sum(abs(mod(phy%k(:,ik)-(phy%k(:,i)-phy%k(:,j))+4d0*pi+1d-9,2d0*pi)))<1d-7) then
					!exit
				!endif
			!enddo
			!if(j==size(phy%k,2)+1) stop "err"
			!!call get_convolution(omega,phy%Gck(:,i)/size(phy%k,2),phy%Gck(:,j),[1,-1],[-1,-1],phy%Sus(ik:ik),[0d0])
			!phy%Sus(ik)=phy%Sus(ik)+integrate(omega,A=-imag(phy%Gck(:,i)*phy%Gck(:,j))*ff_,x=omega)/(pi*size(phy%k,2))
		!enddo
		!!if(abs(phy%Sus(ik))>1d-8) then
			!write(*,*)phy%Sus(ik:ik)
			!read(*,*)
		!!endif
	!enddo
end program
