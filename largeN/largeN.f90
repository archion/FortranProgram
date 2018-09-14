include "fn.f90"
include "../test/pade.f90"
include "../lib/utility.f90"
include "mbroyden.f90"
module global
	use M_fn
	use M_solution
	use M_const
	use M_utility
	use M_matrix
	use mkl_service
	use M_pade
	implicit none
	real(wp) :: alpha=0.3_wp,r=0.3_wp,kp=0.3_wp,fcut=1._wp,fbeta=1._wp,phicut=5e-2_wp,phibeta=200._wp,Jk=0.983_wp,g=0.2_wp
	real(wp) :: Tk,beta,norm_Ac=1._wp,norm_Aphi=1._wp,bcut=10._wp,min_omega=1e-9_wp
	!integer, parameter :: nlog=31,nl=7,nadd(3)=[5,10,25]
	!integer, parameter :: nlog=31,nl=7,nadd(3)=[5,10,1]
	integer, parameter :: nlog=35,nl=14,nadd(3)=[5,10,25]
	!integer, parameter :: nlog=35,nl=40,nadd(3)=[25,25,75]
	!integer, parameter :: nlog=35,nl=34,nadd(2)=[5,10]
	integer, parameter :: n=(nl*nlog+sum(nadd))*2
	integer, parameter :: set_rGB=1,set_Gt=2,set_rG=3,set_rSE=4,set_SE=5,export_Gw=1,export_Gt=2,export_fe=3,export_entropy=4,export_SE=5
	real(wp) :: omega(n),tau(n)
	logical :: is_expinit=.true.
	type gfs(n)
		integer, len :: n
		real(wp) :: rhoc(n),rhophi(n),lambda
		complex(wp) :: Gf(n),GB(n),Gc(n),Gphi(n),Gc0(n),Gphi0(n),Tmat(n),Sus(n)
		complex(wp) :: SEf(n),SEB(n),SEc(n),SEphi(n)
		real(wp) :: Gft(n),GBt(n),Gc0t(n),Gphi0t(n),Gct(n),Gphit(n),SEft(n),SEBt(n),Tmatt(n),Sust(n),SEct(n),SEphit(n)
		logical :: conv=.false.
	contains
		procedure :: self_consistent,analytic,set_Gfunction
		procedure :: entropy,fenergy
		procedure :: export,init
	end type
contains
	subroutine init(self,ut)
		class(gfs(*)) :: self
		integer :: ut
		integer :: i,e
		real(wp) :: omg(2)
		complex(wp) :: dat(2,2)

		is_expinit=.true.
		norm_Ac=1._wp
		do i=1,n
			self%rhoc(i)=A_c(omega(i))
			self%rhophi(i)=A_phi(omega(i))
		enddo
		norm_Ac=1._wp/integrate(omega,A=self%rhoc)
		norm_Aphi=1._wp/integrate(omega(n/2+1:),A=self%rhoc(n/2+1:))
		self%rhoc=self%rhoc*norm_Ac
		self%Gc0%im=-self%rhoc*pi
		self%Gphi0%im=-self%rhophi*pi

		i=1
		rewind(ut)
		read(ut,"(*(e28.20))")omg(1),dat(1,:)
		o:	do 
			read(ut,"(*(e28.20))",iostat=e)omg(2),dat(2,:)
			if(e==-1) then
				backspace(ut)
				backspace(ut)
				backspace(ut)
				read(ut,"(*(e28.20))")omg(1),dat(1,:)
				dat(2,:)=dat(2,:)+(dat(2,:)-dat(1,:))/(omg(2)-omg(1))*(-omega(1)-omg(2))
				omg(2)=-omega(1)
				read(ut,"(*(e28.20))")omg(1),dat(1,:)
			endif
			do 
				if(omg(2)-omega(i)>=0._wp) then
					dat(1,:)=dat(1,:)+(dat(2,:)-dat(1,:))/(omg(2)-omg(1))*(omega(i)-omg(1))
					self%SEf(i)=dat(1,1)
					self%SEB(i)=dat(1,2)
					omg(1)=omega(i)
					i=i+1
					if(i==n+1) exit o
				else
					omg(1)=omg(2)
					dat(1,:)=dat(2,:)
					exit
				endif
			enddo
		enddo o
	end subroutine
	subroutine export(self,ut,flag)
		class(gfs(*)) :: self
		integer :: ut(:)
		integer :: flag(:)
		integer :: i,k
		do k=1,size(flag)
			select case(flag(k))
			case(export_Gw)
				write(ut(k),"(A)")"Tk omega rGf iGf rGB iGB rSEf iSEf rSEB iSEB rGc iGc rGphi iGphi rhoc rhophi rSEc iSEc rSEphi iSEphi rTmat iTmat rSus iSus rGc0 iGc0 rGphi0 iGphi0"
				do i=1,n
					write(ut(k),"(*(e28.20))")Tk,omega(i),self%Gf(i),self%GB(i),self%SEf(i),self%SEB(i),self%Gc(i),self%Gphi(i),self%rhoc(i),self%rhophi(i),self%SEc(i),self%SEphi(i),self%Tmat(i),self%Sus(i),self%Gc0(i),self%Gphi0(i)
				enddo
				write(ut(k),"(x/)")
			case(export_Gt)
				write(ut(k),"(A)")"Tk tau Gft GBt SEft SEBt Gct Gphit SEct SEphit Tmatt Sust Gc0t Gphi0t"
				do i=1,n
					write(ut(k),"(*(e28.20))")Tk,tau(i),self%Gft(i),self%GBt(i),self%SEft(i),self%SEBt(i),self%Gct(i),self%Gphit(i),self%SEct(i),self%SEphit(i),self%Tmatt(i),self%Sust(i),self%Gc0t(i),self%Gphi0t(i)
				enddo
				write(ut(k),"(x/)")
			case(export_SE)
				write(*,*)"export initial...."
				rewind(ut(k))
				is_expinit=.false.
				do i=1,n
					write(ut(k),"(*(e28.20))")omega(i),self%SEf(i),self%SEB(i)
				enddo
			end select
		end do
	end subroutine
	subroutine set_Gfunction(self,flag)
		class(gfs(*)) :: self
		integer :: flag(:),i
		do i=1,size(flag)
			select case(flag(i))
			case(set_rG)
				! make GB a usual Green's function that satisfy kk relation
				call set_realpart(omega,self%GB)
				call set_realpart(omega,self%Gf)
			case(set_Gt)
				call get_Gtau(omega,G=self%Gf,eta=-1,time=tau,rt=self%Gft)
				call get_Gtau(omega,G=self%GB,eta=1,time=tau,rt=self%GBt)
				call get_Gtau(omega,G=self%Gc0,eta=-1,time=tau,rt=self%Gc0t)
				call get_Gtau(omega,G=self%Gphi0,eta=1,time=tau,rt=self%Gphi0t)
				call get_Gtau(omega,G=self%SEf,eta=-1,time=tau,rt=self%SEft)
				call get_Gtau(omega,G=self%SEB,eta=1,time=tau,rt=self%SEBt)
				call get_Gtau(omega,G=self%SEc,eta=-1,time=tau,rt=self%SEct)
				call get_Gtau(omega,G=self%SEphi,eta=1,time=tau,rt=self%SEphit)
				call get_Gtau(omega,G=self%Gc,eta=-1,time=tau,rt=self%Gct)
				call get_Gtau(omega,G=self%Gphi,eta=1,time=tau,rt=self%Gphit)
				call get_Gtau(omega,G=self%Tmat,eta=-1,time=tau,rt=self%Tmatt)
				call get_Gtau(omega,G=self%Sus,eta=1,time=tau,rt=self%Sust)
			case(set_rSE)
				call set_realpart(omega,self%SEf)
				call set_realpart(omega,self%SEB)
			case(set_SE)
				call get_convolution(omega,A=-self%Gf,B=self%GB+Jk,tau=[1,-1],eta=[-1,1],rt=self%SEc)
				call set_realpart(omega,self%SEc)
				call get_convolution(omega,A=g**2*self%Gf,B=self%Gf,tau=[1,-1],eta=[-1,-1],rt=self%SEphi)
				call set_realpart(omega,self%SEphi)
			end select
		enddo
	end subroutine
	real(wp) elemental function A_phi(x)
		real(wp), intent(in) :: x
		if(abs(x)<eps) then
			A_phi=0._wp
		elseif(abs(x)<phicut) then
			A_phi=norm_Aphi*abs(x)**(1._wp-alpha)*sign(1._wp,x)
		else
			A_phi=norm_Aphi*phicut**(1._wp-alpha)*sign(1._wp,x)*1._wp/((exp(phibeta*(x-phicut-4._wp/phibeta))+1._wp)*(exp(phibeta*(-x-phicut-4._wp/phibeta))+1._wp))
		endif
	end function
	real(wp) elemental function A_c(x)
		real(wp), intent(in) :: x
		if(abs(x)<fcut) then
			!A_c=norm_Ac*(abs(x)/fcut)**r*exp(-fcut**2/pi*fbeta)
			A_c=norm_Ac*(abs(x)/fcut)**r
			!A_c=norm_Ac*abs(x)**r
			!A_c=norm_Ac*abs(x)**r
		elseif(abs(x)<=(-omega(1))) then
			!A_c=norm_Ac*exp(-x**2/pi*fbeta)
			A_c=norm_Ac*exp(-(x**2-fcut**2)/pi*fbeta)
			!A_c=norm_Ac*abs(fcut)**r*exp(-fbeta*(abs(x)-fcut)**2)
			!A_c=norm_Ac*fcut**r*ff(fbeta*(x-fcut-15.0/fbeta))*ff(fbeta*(-x-fcut-15.0/fbeta))
			!A_c=norm_Ac*fcut**r*1d0/(exp(fbeta*(x-fcut-15.0/fbeta))+1d0)*1d0/(exp(fbeta*(-x-fcut-15.0/fbeta))+1d0)
		else
			A_c=0._wp
		endif
	end function
	subroutine get_Giomg(omega,G,iomg,rt)
		real(wp) :: omega(:)
		complex(wp) :: G(:),iomg(:),rt(:)
		integer :: i
		do i=1,size(iomg)
			rt(i)=-integrate(omega,G%im/(iomg(i)-omega))/pi
		enddo
	end subroutine
	subroutine get_Gtau(omega,G,Gf,eta,time,rt)
		real(wp) :: omega(:),time(:),rt(:)
		complex(wp), optional :: G(:),Gf
		integer :: eta
		integer :: tau,i,k,n
		real(wp) :: A1,A2,B1,B2,e1,e2,t
		procedure(real(wp)), pointer :: f
		n=size(omega)
		tau=1
		if((time(size(time))+time(1))>0._wp) then
			tau=-1
		endif
		do k=1,size(time)
			t=time(k)
			e1=omega(1)
			if(present(G)) then
				A1=G(1)%im
			else

				A1=imag(Gf(e1))
			endif
			B1=tau*1._wp/(exp((tau*beta+t)*e1)-eta*exp(e1*t))
			rt(k)=0._wp
			do i=1,n-1
				e2=omega(i+1)
				if(present(G)) then
					A2=G(i+1)%im
				else
					A2=imag(Gf(e2))
				endif
				B2=tau*1._wp/(exp((tau*beta+t)*e2)-eta*exp(e2*t))
				rt(k)=rt(k)+0.5_wp*(e2-e1)*(A1*B1+A2*B2)
				A1=A2
				B1=B2
				e1=e2
			enddo
			rt(k)=1._wp/pi*eta*rt(k)
		enddo
	end subroutine
	subroutine get_convolution_new(omega,A,B,tau,eta,rt,noreset)
		real(wp) :: omega(:)
		complex(wp) :: rt(:)
		complex(wp) :: A(:),B(:)
		integer :: tau(2),eta(2)
		logical, optional :: noreset
		integer :: i,j,k,n,m,n1,n2,sg
		real(wp) :: e1,e2,e0,omg,omegaB(size(omega))
		complex(wp) :: A1,A2,B1,B2,rt_
		logical :: flag
		if(.not.present(noreset)) rt=0._wp
		n=size(omega)
		call pade(x=cmplx(omega(:nadd(3)),kind=wp),y=A(:nadd(3)),id=1)
		call pade(x=cmplx(omega(n-nadd(3)+1:n),kind=wp),y=A(n-nadd(3)+1:n),id=2)
		call pade(x=cmplx(omega(:nadd(3)),kind=wp),y=B(:nadd(3)),id=3)
		call pade(x=cmplx(omega(n-nadd(3)+1:n),kind=wp),y=B(n-nadd(3)+1:n),id=4)
		sg=-tau(1)*tau(2)
		!$OMP PARALLEL DO PRIVATE(omg,A1,B1,A2,B2,e1,e2,n1,n2,flag,rt_)
		do k=1,n
			e1=nan
			omg=omega(k)
			n1=0
			n2=merge(0,n+1,sg>0)
			rt_=0d0
			A1=0d0
			B1=0d0
			do
				if(n2==merge(n,1,sg>0).and.n1==n) then
					exit
				elseif(n2==merge(n,1,sg>0)) then
					flag=.true.
				elseif(n1==n) then
					flag=.false.
				elseif(omega(n1+1)<tau(1)*(omg-tau(2)*omega(n2+sg))) then
					flag=.true.
				else
					flag=.false.
				endif
				if(flag) then
					n1=n1+1
					e2=omega(n1)
					A2=A(n1)
					if(any(n2==merge([n,0],[1,n+1],sg>0))) then
						B2=Bf(tau(2)*(omg-tau(1)*e2))
						!B2=0d0
					else
						B2=B(n2)+(B(n2+sg)-B(n2))/(omega(n2+sg)-omega(n2))*(tau(2)*(omg-tau(1)*e2)-omega(n2))
					endif
				else
					n2=n2+sg
					e2=tau(1)*(omg-tau(2)*omega(n2))
					B2=B(n2)
					if(n1==n.or.n1==0) then
						A2=Af(e2)
						!A2=0d0
					else
						A2=A(n1)+(A(n1+1)-A(n1))/(omega(n1+1)-omega(n1))*(e2-omega(n1))
					endif
				endif
				if(isnan(e1)) then
					e0=e2
				!elseif(abs(A1)*abs(B1)>eps.and.abs(A2)*abs(B2)>eps) then
				!elseif((abs(A1)*abs(B1)>eps.and.omg<0d0).or.(abs(A2)*abs(B2)>eps.and.omg>0._wp)) then
				!elseif(abs(A2)*abs(B2)>eps) then
				else
					rt_=rt_+0.5_wp*(e2-e1)*cmplx(&
						(-tau(1)*f1(tau(1)*e1)*imag(A1)*real(B1)+tau(2)*eta(1)*eta(2)*f2((tau(1)*e1-omg))*real(A1)*imag(B1))+&
						(-tau(1)*f1(tau(1)*e2)*imag(A2)*real(B2)+tau(2)*eta(1)*eta(2)*f2((tau(1)*e2-omg))*real(A2)*imag(B2))&
						,&
						(-tau(1)*tau(2)*f1(tau(1)*e1)*imag(A1)*imag(B1)+tau(1)*tau(2)*eta(1)*eta(2)*f2((tau(1)*e1-omg))*imag(A1)*imag(B1))+&
						(-tau(1)*tau(2)*f1(tau(1)*e2)*imag(A2)*imag(B2)+tau(1)*tau(2)*eta(1)*eta(2)*f2((tau(1)*e2-omg))*imag(A2)*imag(B2)),kind=wp)
				endif
				!write(*,*)e2,abs(A2)*abs(B2),rt_%re
				!read(*,*)
				A1=A2
				B1=B2
				e1=e2
				!write(*,"(3es12.4,3i4)")e2,omega(n1),omegaB(n2),n1,n2,n
				!read(*,*)
			enddo
			rt(k)=rt(k)-eta(1)*rt_/pi
			!rt(n+1-k)=rt(n+1-k)+cmplx(-eta(2)*rt_%re,eta(2)*rt_%im,kind=8)/pi
			!rt(k)=rt(k)-eta(1)*(integrate(f,[-inf,e0],1e-5_wp)+integrate(f,[e2,inf],1e-5_wp))/pi
			!rt(k)=rt(k)-eta(1)*(integrate(f,[-50._wp,e0],1e-5_wp)+integrate(f,[e2,50._wp],1e-5_wp))/pi
		enddo
		!$OMP END PARALLEL DO
		call pade(y=[complex(wp)::],id=0)
	contains
		real(wp) function f1(x) result(y)
			real(wp) :: x
			select case(eta(1))
			case(1)
				y=fb(beta*x)
			case(-1)
				y=ff(beta*x)
			end select
		end function
		real(wp) function f2(x) result(y)
			real(wp) :: x
			select case(eta(2))
			case(1)
				y=fb(beta*x)
			case(-1)
				y=ff(beta*x)
			end select
		end function
		complex(wp) function Af(x) result(y)
			real(wp) :: x
			complex(wp) :: y_(1)
			if(x<0._wp) then
				call pade(xp=[cmplx(x,kind=wp)],y=y_,id=1)
			else
				call pade(xp=[cmplx(x,kind=wp)],y=y_,id=2)
			endif
			y=y_(1)
		end function
		complex(wp) function Bf(x) result(y)
			real(wp) :: x
			complex(wp) :: y_(1)
			if(x<0._wp) then
				call pade(xp=[cmplx(x,kind=wp)],y=y_,id=3)
			else
				call pade(xp=[cmplx(x,kind=wp)],y=y_,id=4)
			endif
			y=y_(1)
		end function
		complex(wp) function f(x) result(y)
			real(wp) :: x
			complex(wp) :: A,B
			A=Af(x)
			B=Bf(tau(2)*(omg-tau(1)*x))
			y=cmplx(&
				(-tau(1)*f1(beta*tau(1)*x)*imag(A)*real(B)+tau(2)*eta(1)*eta(2)*f2(beta*(tau(1)*x-omg))*real(A)*imag(B))&
				,&
				(-tau(1)*tau(2)*f1(beta*tau(1)*x)*imag(A)*imag(B)+tau(1)*tau(2)*eta(1)*eta(2)*f2(beta*(tau(1)*x-omg))*imag(A)*imag(B)),kind=wp)
		end function
	end subroutine
	subroutine get_convolution(omega,A,B,tau,eta,rt,noreset)
		real(wp) :: omega(:)
		complex(wp) :: rt(:)
		complex(wp) :: A(:),B(:)
		integer :: tau(2),eta(2)
		logical, optional :: noreset
		integer :: i,j,k,n,iwb(size(omega)),ibw(size(omega)),m
		real(wp) :: e1,e2,omg
		complex(wp) :: A1,A2,B1,B2,rt_
		if(.not.present(noreset)) rt=0._wp
		n=size(omega)
		!$OMP PARALLEL DO PRIVATE(omg,ibw,iwb,A1,B1,A2,B2,e1,e2,m,rt_)
		!do k=1,n/2
			do k=1,n
			omg=omega(k)
			do i=1,n
				ibw(i)=find(omega,tau(1)*(omg-tau(2)*omega(i)))
			end do
			A1=A(1)
			do i=1,n
				iwb(i)=find(omega,tau(2)*(omg-tau(1)*omega(i)))
			end do
			i=iwb(1)
			B1=merge(B(i)+(B(i+1)-B(i))/(omega(i+1)-omega(i))*(tau(2)*(omg-tau(1)*omega(1))-omega(i)),cmplx(0._wp,kind=wp),i>0.and.i<n)
			e1=omega(1)
			rt_=0._wp
			!write(*,*)e1,abs(A1)*abs(B1),rt_%re
			!read(*,*)
			do i=1,n-1
				do j=iwb(i),iwb(i+1),sign(1,iwb(i+1)-iwb(i))
					if(j==iwb(i+1)) then
						e2=omega(i+1)
						A2=A(i+1)
						B2=merge(B(j)+(B(j+1)-B(j))/(omega(j+1)-omega(j))*(tau(2)*(omg-tau(1)*e2)-omega(j)),cmplx(0._wp,kind=wp),j>0.and.j<n)
					else
						m=j+(1-tau(1)*tau(2))/2
						e2=tau(1)*(omg-tau(2)*omega(m))
						B2=B(m)
						m=ibw(m)
						A2=A(m)+(A(m+1)-A(m))/(omega(m+1)-omega(m))*(e2-omega(m))
					endif
					rt_=rt_+0.5_wp*(e2-e1)*cmplx(&
						(-tau(1)*f1(tau(1)*e1)*imag(A1)*real(B1)+tau(2)*eta(1)*eta(2)*f2((tau(1)*e1-omg))*real(A1)*imag(B1))+&
						(-tau(1)*f1(tau(1)*e2)*imag(A2)*real(B2)+tau(2)*eta(1)*eta(2)*f2((tau(1)*e2-omg))*real(A2)*imag(B2))&
						,&
						(-tau(1)*tau(2)*f1(tau(1)*e1)*imag(A1)*imag(B1)+tau(1)*tau(2)*eta(1)*eta(2)*f2((tau(1)*e1-omg))*imag(A1)*imag(B1))+&
						(-tau(1)*tau(2)*f1(tau(1)*e2)*imag(A2)*imag(B2)+tau(1)*tau(2)*eta(1)*eta(2)*f2((tau(1)*e2-omg))*imag(A2)*imag(B2)),kind=wp)
					A1=A2
					B1=B2
					e1=e2
				enddo
			enddo
			rt(k)=rt(k)-eta(1)*rt_/pi
			!rt(n+1-k)=rt(n+1-k)+cmplx(-eta(2)*rt_%re,eta(2)*rt_%im,kind=8)/pi
		enddo
		!$OMP END PARALLEL DO
	contains
		real(wp) function f1(x) result(y)
			real(wp) :: x
			select case(eta(1))
			case(1)
				y=fb(beta*x)
			case(-1)
				y=ff(beta*x)
			end select
		end function
		real(wp) function f2(x) result(y)
			real(wp) :: x
			select case(eta(2))
			case(1)
				y=fb(beta*x)
			case(-1)
				y=ff(beta*x)
			end select
		end function
	end subroutine
	subroutine set_realpart(omega,G)
		real(wp) :: omega(:)
		complex(wp) :: G(:)
		real(wp) ::dx
		integer :: i,j,dn(2)
		G(1)%re=0._wp
		G(n)%re=0._wp
		do i=2,n-1
			G(i)%re=2._wp*G(i)%im*log((omega(n)-omega(i))/(omega(i)-omega(1)))
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
			G(i)%re=0.5_wp/pi*G(i)%re
		enddo
	end subroutine
	subroutine  self_consistent(self,tol,niter,rate)
		class(gfs(*)) :: self
		real(wp) :: tol,rate
		integer :: niter(:)
		real(wp) :: x(n*2*2),x_(n*2*2),nf(1),tmp(1)
		logical :: flag
		integer :: i1,i2,k,i,j,ei
		self%conv=.false.
		self%Gf=0._wp
		self%GB=0._wp
		self%Gc=0._wp
		self%Gphi=0._wp
		call mbroyden(0,x(:n*2),x_(:n*2),rate,nsave=40,id=2)
		self%lambda=0._wp
		do i2=1,niter(2)
			flag=.false.
			ei=0
			call set_realpart(omega,self%Gc0)
			call set_realpart(omega,self%Gphi0)
			!call get_Gtau(omega,G=self%Gc0,eta=-1,time=[eps],rt=tmp)
			do i1=1,niter(1)
				call set_realpart(omega,self%SEf)
				call set_realpart(omega,self%SEB)
				!self%lambda=-(self%SEf(1)%re+self%SEf(n)%re)/2d0
				!underscore=next(0._wp,0._wp,0.1_wp,ytol=1e-5_wp)
				!do 
					self%Gf=1._wp/(omega-self%lambda-self%SEf)
					call get_Gtau(omega,G=self%Gf,eta=-1,time=[-eps],rt=nf)
					!if(next(self%lambda,nf(1)-0.5_wp)) then
						!exit
					!endif
				!enddo
				self%GB=1._wp/(-1._wp/Jk-self%SEB)
				!self%GB=1d0/(-1d0/Jk-self%SEB)+Jk
				!call self%export([233],[export_Gw])
				!call set_realpart(omega,self%GB)
				!call self%export([233],[export_Gw])
				!self%Gc=1._wp/(omega+img*0.0001_wp)
				!call get_convolution_new(omega,A=self%Gc,B=self%Gf,tau=[1,1],eta=[-1,-1],rt=self%SEB)
				!call self%export([233],[export_Gw])
				!call get_convolution_new(omega,A=self%SEB,B=1._wp/self%Gc,tau=[-1,-1],eta=[1,-1],rt=self%Gf)
				!call self%export([233],[export_Gw])
				!stop "debug"



				x=[self%SEf%im,self%SEB%im,self%SEf%re,self%SEB%re]

				self%SEB=0._wp
				self%SEf=0._wp
				call get_convolution_new(omega,A=-kp*self%Gc0,B=self%GB,tau=[1,1],eta=[-1,1],rt=self%SEf,noreset=underscore)
				!call self%export([233],[export_Gw])
				!call set_realpart(omega,self%SEf)
				!call self%export([233],[export_Gw])
				!stop
				call get_convolution_new(omega,A=self%Gc0,B=self%Gf,tau=[1,-1],eta=[-1,-1],rt=self%SEB,noreset=underscore)
				if(abs(g)>0._wp) then
					call get_convolution_new(omega,A=-2d0*g**2*self%Gphi0,B=self%Gf,tau=[1,1],eta=[1,-1],rt=self%SEf,noreset=underscore)
				endif

				x_=[self%SEf%im,self%SEB%im,self%SEf%re,self%SEB%re]
				if(all(abs(x(:n*2)-x_(:n*2))/(abs(x_(:n*2))+eps)<tol)) then
					!if(ei==0) ei=1
					!if(ei>20) then
						write(*,*)"converged"


						if(abs(0.5_wp-nf(1))>1e-2_wp) then
							is_expinit=.false.
							write(*,*)"nf is not correct****************************",nf(1)
							return
						else
							self%conv=.true.
						endif
						exit
					!endif
				endif
				if(i1==niter(1)) then
					write(*,*)"not converged",maxval(abs(x-x_)/(abs(x_)+eps)),"***************"
					is_expinit=.false.
					return
				endif

				!if(flag) then
				!if(ei>0) then
					!ei=ei+1
					!write(*,*)ei,maxval(abs(x-x_)/(abs(x_)+eps))
					!if(abs(Tk-1e-9_wp)<1e-4_wp) then
					!x=(1._wp-rate)*x+rate*x_
					!x=(x+rate*x_)/(1._wp+rate)
					!endif

					!ei=30
					!x=x_
				!else
					call mbroyden(i1,x(:n*2),x_(:n*2),id=1)
				!endif
				!x=x_


				self%SEf%im=x(1:n)
				self%SEB%im=x(n+1:2*n)
				!self%SEf%re=x(2*n+1:3*n)
				!self%SEB%re=x(3*n+1:4*n)
				write(*,*)i1,maxval(abs(x(:n*2)-x_(:n*2))/(abs(x_(:n*2))+eps))
				!write(*,*)"free energy:",free_energy(rhoc,rhophi,Af,AB,iSEf,iSEB,1)
				!read(*,*)
			enddo
			write(*,*)i1,maxval(abs(x-x_)/(abs(x_)+eps)),"nf:",nf*2._wp

			x=[self%Gc0%im,self%Gphi0%im]

			if(niter(2)>1) then
				!call get_convolution(omega,A=-self%Gf,B=self%GB,tau=[1,-1],eta=[-1,1],rt=self%SEc)

				self%Gc=self%Gc0+self%Gc0*self%SEc*self%Gc0

				do i=1,n
					self%Gc0=1._wp/(omega(i)-omega-self%SEc(i))
					self%Gc(i)=integrate(omega,A=self%rhoc*self%Gc0)
				enddo

				self%Gc0=1._wp/(1._wp/self%Gc+self%SEc)
			endif
			x_=[self%Gc0%im,self%Gphi0%im]

			if(all(abs(x-x_)/(abs(x_)+eps)<tol)) then
				write(*,*)"converged"
				exit
			else
				write(*,*)"|___",i2,maxval(abs(x-x_)/(abs(x_)+eps))
				if(i2==niter(2)) then
					write(*,*)"not converged"
				endif
			endif

			call mbroyden(i2,x,x_,id=2)
			!x=rate*x_+(1d0-rate)*x
			!x=x_
			self%Gc0%im=x(1:n)
			self%Gphi0%im=x(n+1:)
		enddo
		write(*,*)i2,maxval(abs(x-x_)/(abs(x_)+eps))
		call mbroyden(huge(1),x,x_)
	end subroutine
	real(wp) function entropy(self) result(rt)
		class(gfs(*)) :: self
		real(wp) :: df(n),db(n),tmp1(n),tmp2(n)
		integer :: i
		df=omega*(0.5_wp*beta/cosh(0.5_wp*omega*beta))**2
		db=omega*(0.5_wp*beta/sinh(0.5_wp*omega*beta))**2

		!rt=-((integrate(omega,A=&
			!ff_*(imag(log(-1d0/Gf)+Gf*(omega-1d0/Gf))+0.5d0*(1d0+sign(1d0,omega))*pi)&
			!+kp*fb_*imag(log(-1d0/GB)+GB*(-1d0/Jk-1d0/GB))&
			!!+kp*fb_*imag(log(-1d0/GB)+GB*SEB)&
			!!ff_*(imag(log(-1d0/Gf))+0.5d0*(1d0+sign(1d0,omega))*pi)&
			!!+kp*fb_*imag(log(-1d0/GB))&
			!)/pi&
				!-1d0/beta*log(2d0)&
				!+integrate(tau,A=-kp*(-Gft)*Gct*GBt)&
			!))
		!call get_convolution(omega,A=-self%Gf,B=self%GB+Jk,tau=[1,-1],eta=[-1,1],rt=self%SEc)
		!call set_realpart(omega,self%SEc)
		rt=-(integrate(omega,A=&
			df*(imag(log(-1d0/self%Gf))+self%Gf%re*imag(omega-1d0/self%Gf)+0.5d0*(1d0+sign(1d0,omega))*pi)&
			!df*(imag(log(-1d0/self%Gf))+self%Gf%re*self%SEf%im+0.5d0*(1d0+sign(1d0,omega))*pi)&
			!df*(imag(log(-1d0/self%Gf))+self%Gf%re*self%SEf%im)&
			+kp*db*(imag(log(-1d0/self%GB))+(self%GB%re+Jk)*imag(-1d0/Jk-1d0/self%GB))&
			!+kp*db*(imag(log(-1d0/self%GB))+(self%GB%re)*self%SEB%im)&
			-kp*df*(self%Gc0%im*self%SEc%re)&
			-db*(self%Gphi0%im*self%SEphi%re)&
			)/pi&
			-log(2d0)&
			)
		!write(*,"(es12.4,A)")-(integrate(omega,A=df*(-kp*self%SEc%re*self%Gc0%im)))/pi," SEc*Gc"
		do i=1,n
			tmp1(i)=integrate(omega,A=fb(beta*(omega-omega(i)))*self%Gf%im)
		enddo
		do i=1,n
			tmp2(i)=integrate(omega,A=fb(beta*(omega(i)-omega))*self%Gc0%im)
		enddo
		rt=rt+integrate(omega,A=df*(tmp1*self%Gc0%im-tmp2*self%Gf%im))/pi**2*Jk*kp
			
	end function
	real(wp) function fenergy(self) result(rt)
		class(gfs(*)) :: self
		integer :: i
		integer :: cs

		write(*,*)integrate(tau,A=(-self%Gft)*self%SEft),"Gft*SEft"
		write(*,*)integrate(tau,A=-kp*self%GBt*self%SEBt),"GBt*SEBt"
		write(*,*)integrate(tau,A=-kp*(-self%Gft)*self%Gc0t*self%GBt-g**2*(-self%Gft)*self%Gphi0t*self%Gft),"GftGc0tGBt"
		write(*,*)integrate(omega,A=omega*ff(beta*omega)*(-self%Gf%im/pi)),"rhoGf0^-1"
		write(*,*)integrate(omega,A=-ff(beta*omega)*imag(self%Gf*self%SEf)/pi),"Gf*SEf"
		write(*,*)integrate(omega,A=fb(beta*omega)*imag((self%GB+Jk)*self%SEB)/pi*(-kp)),"(GB-Jk)*SEB"
		write(*,*)integrate(omega,A=-ff(beta*omega)*imag(self%Gc0*self%SEc)/pi*kp-fb(beta*omega)*imag(self%Gphi0*self%SEphi)/pi),"Gc0*SEc"

		rt=(integrate(omega,A=&
			ff(beta*omega)*(imag(log(-1._wp/self%Gf)+self%Gf*(omega-1._wp/self%Gf))+0.5_wp*(1._wp+sign(1._wp,omega))*pi)&
			!!+kp*fb_*imag(log(-1._wp/self%GB)+self%GB*(-1d0/Jk-1d0/self%GB))&
			!ff(beta*omega)*(imag(log(-1._wp/self%Gf)+self%Gf*self%SEf)+0.5_wp*(1._wp+sign(1._wp,omega))*pi)&
			!+kp*fb(beta*omega)*imag(log(-1._wp/self%GB)+(self%GB+Jk)*self%SEB)&
			+kp*fb(beta*omega)*imag(log(-1._wp/self%GB)+(self%GB+Jk)*(-1._wp/Jk-1._wp/self%GB))&
			!ff(beta*omega)*(imag(log(-1._wp/self%Gf))+0.5_wp*(1._wp+sign(1._wp,omega))*pi)&
			!+kp*fb(beta*omega)*imag(log(-1._wp/self%GB))&
			)/pi&
			-1._wp/beta*log(2._wp)&
			!!+integrate(tau,A=-kp*(-self%Gft)*self%Gct*self%GBt)&
			!-integrate(omega,A=ff(beta*omega)*imag(self%Gf*self%SEf)/pi)&
			-integrate(omega,A=ff(beta*omega)*imag(self%Gc0*self%SEc)/pi*kp)&
			-integrate(omega,A=fb(beta*omega)*imag(self%Gphi0*self%SEphi)/pi)&
			)

		!rt=integrate(omega,A=omega*ff(beta*omega)*(-self%Gf%im/pi))
		!rt=integrate(omega,A=-ff(beta*omega)*imag(self%Gf*self%SEf)/pi)
		!call get_convolution(omega,A=-self%Gf,B=self%GB+Jk,tau=[1,-1],eta=[-1,1],rt=self%SEc)
		!call set_realpart(omega,self%SEc)
		!call set_realpart(omega,self%Gc0)
		write(26,"(e28.20$)")Tk,&
			!integrate(tau,A=-kp*(-self%Gft)*self%Gc0t*self%GBt),&
			!integrate(tau,A=-kp*self%GBt*self%SEBt),&
			!integrate(tau,A=(-self%Gft)*self%SEft),&
			!integrate(omega,A=ff(beta*omega)*(imag(log(-1._wp/self%Gf)+self%Gf*(omega-1._wp/self%Gf))+0.5_wp*(1._wp+sign(1._wp,omega))*pi))/pi-1._wp/beta*log(2._wp),&
			!integrate(omega,A=kp*fb(beta*omega)*imag(log(-1._wp/self%GB)+(self%GB+Jk)*(-1._wp/Jk-1._wp/self%GB)))/pi,&
			integrate(omega,A=ff(beta*omega)*(imag(log(-1._wp/self%Gf)+self%Gf*self%SEf)+0.5_wp*(1._wp+sign(1._wp,omega))*pi))/pi-1._wp/beta*log(2._wp),&
			integrate(omega,A=kp*fb(beta*omega)*imag(log(-1._wp/self%GB)+(self%GB+Jk)*self%SEB))/pi,&
			rt,&
			!integrate(omega,A=ff(beta*omega)*(imag(log(-1._wp/self%Gf))+0.5_wp*(1._wp+sign(1._wp,omega))*pi))/pi-1._wp/beta*log(2._wp),&
			!integrate(omega,A=kp*fb(beta*omega)*imag(log(-1._wp/self%GB)))/pi,&
			integrate(omega,A=-ff(beta*omega)*imag(self%Gf*self%SEf)/pi),&
			!integrate(omega,A=omega*ff(beta*omega)*(-self%Gf%im/pi)),&
			integrate(omega,A=fb(beta*omega)*imag((self%GB+Jk)*self%SEB)/pi*(-kp)),&
			integrate(omega,A=-ff(beta*omega)*imag(self%Gc0*self%SEc)/pi*kp)
		!cs=1
		!write(26,"(e28.20$)")integrate(ft,[0._wp+eps,beta-eps],epsilon(1._dp)*10._wp)
		!cs=2
		!write(26,"(e28.20$)")integrate(ft,[0._wp+eps,beta-eps],epsilon(1._dp)*10._wp)
		!cs=3
		!write(26,"(e28.20$)")integrate(ft,[0._wp+eps,beta-eps],epsilon(1._dp)*10._wp)
		!rt=integrate(omega,A=-ff(beta*omega)*imag(self%Gf*self%SEf)/pi)
		!rt=integrate(tau,A=(-self%Gft)*self%SEft)
		!rt=integrate(tau,A=-kp*self%GBt*self%SEBt)
		!rt=integrate(tau,A=-kp*(-self%Gft)*self%Gc0t*self%GBt)
		!rt=integrate(omega,A=fb(beta*omega)*imag((self%GB+Jk)*self%SEB)/pi*(-kp))
		write(26,"(x)")
	contains
		real(wp) function ft(x) result(rt)
			real(wp) :: x
			real(wp) :: rt_(1)
			rt=1._wp
			select case(cs)
			case(1)
				call get_Gtau(omega,G=self%GB,eta=1,time=[x],rt=rt_)
				rt=rt*rt_(1)
				call get_Gtau(omega,G=self%Gf,eta=-1,time=[x],rt=rt_)
				rt=rt*rt_(1)
				call get_Gtau(omega,G=self%Gc0,eta=-1,time=[x],rt=rt_)
				rt=kp*rt*rt_(1)
			case(2)
				call get_Gtau(omega,G=self%GB,eta=1,time=[x],rt=rt_)
				rt=rt*rt_(1)
				call get_Gtau(omega,G=self%SEB,eta=1,time=[x],rt=rt_)
				rt=-kp*rt*rt_(1)
			case(3)
				call get_Gtau(omega,G=self%Gf,eta=-1,time=[x],rt=rt_)
				rt=rt*rt_(1)
				call get_Gtau(omega,G=self%SEf,eta=-1,time=[x],rt=rt_)
				rt=-rt*rt_(1)
			end select
		end function
	end function
	subroutine analytic(self,a_f,a_b)
		class(gfs(*)) :: self
		real(wp) :: a_f,a_b
		real(wp) :: A1,B1,A0,t,sh,ch,lsh,lch,lgam_a,lgam_b
		integer :: i,j
		A1=1._wp
		B1=1._wp
		A0=1._wp
		t=(gamma(a_b)*gamma(1._wp+r+a_f)/(pi**2*(1._wp+1._wp/tan(pi*(r+a_f)/2._wp)**2)*gamma(1._wp+r))/(A1*B1*A0))**(1._wp/(a_f+a_b))
		do i=1,n
			if(exp(beta*abs(omega(i)))>1e100_wp) then
				lsh=beta*abs(omega(i))/2._wp-log(2._wp)
				lch=beta*abs(omega(i))/2._wp-log(2._wp)
				lgam_a=log(2._wp*pi)-beta*abs(omega(i))/2._wp+(a_f-1._wp)*(log(beta*abs(omega(i)))-log(2._wp*pi))
				lgam_b=log(2._wp*pi)-beta*abs(omega(i))/2._wp+(a_b-1._wp)*(log(beta*abs(omega(i)))-log(2._wp*pi))
			else
				lsh=log(sinh(beta*abs(omega(i))/2._wp))
				lch=log(cosh(beta*abs(omega(i))/2._wp))
				lgam_a=real(lgamma(a_f/2._wp+img*beta*omega(i)/(2._wp*pi))+lgamma(a_f/2._wp-img*beta*omega(i)/(2._wp*pi)))
				lgam_b=real(lgamma(a_b/2._wp+img*beta*omega(i)/(2._wp*pi))+lgamma(a_b/2._wp-img*beta*omega(i)/(2._wp*pi)))
			endif
			sh=real(exp((a_f-1._wp)*(log(2._wp*pi)-log(beta))+lgam_a+lsh))*sign(1._wp,omega(i))
			ch=real(exp((a_f-1._wp)*(log(2._wp*pi)-log(beta))+lgam_a+lch))
			self%Gf(i)=A1*t**a_f/gamma(a_f)*cmplx(1d0/tan(pi*a_f/2._wp)*sh,-ch,kind=wp)
			sh=real(exp((a_b-1._wp)*(log(2._wp*pi)-log(beta))+lgam_b+lsh))*sign(1._wp,omega(i))
			ch=real(exp((a_b-1._wp)*(log(2._wp*pi)-log(beta))+lgam_b+lch))
			self%GB(i)=B1*t**a_b/gamma(a_b)*cmplx(-tan(pi*a_b/2._wp)*ch,-sh,kind=wp)-Jk
			!self%Gc(i)%im=-pi*A0*(abs(omega(i))/fcut)**r
			self%Gc0(i)%im=-pi*A0*(abs(omega(i))/fcut)**r
			if(abs(omega(i))<fcut) then
				self%Gc0(i)%im=-pi*A0*(abs(omega(i))/fcut)**r
			elseif(abs(omega(i))<=(-omega(1))) then
				self%Gc0(i)%im=-pi*A0*exp(-(omega(i)**2-fcut**2)/pi*fbeta)
			else
				self%Gc0(i)%im=0._wp
			endif
		enddo

		self%SEf=omega-1._wp/self%Gf
		self%SEB=-1._wp/Jk-1._wp/self%GB
		call get_convolution(omega,A=-self%Gf,B=self%GB+Jk,tau=[1,-1],eta=[-1,1],rt=self%SEc)
		call get_convolution(omega,A=-kp*self%Gc0,B=self%GB+Jk,tau=[1,1],eta=[-1,1],rt=self%SEf)
		call get_convolution(omega,A=self%Gc0,B=self%Gf,tau=[1,-1],eta=[-1,-1],rt=self%SEB)
		call set_realpart(omega,self%SEc)
		call set_realpart(omega,self%SEf)
		call set_realpart(omega,self%SEB)
		call set_realpart(omega,self%Gc0)

		!call self%export([22],[export_Gw])
		!call set_realpart(omega,self%Gf)
		!call set_realpart(omega,self%GB)
		!call self%export([22],[export_Gw])
		!stop

		self%Gft=-A1*(pi*t/((sin(pi*tau/beta))*beta))**a_f
		self%GBt=-B1*(pi*t/((sin(pi*tau/beta))*beta))**a_b
		self%Gc0t=-A0*(t/tau)**(r+1._wp)
		self%Gft(n/2+1:)=self%Gft(n/2:1:-1)
		self%GBt(n/2+1:)=self%GBt(n/2:1:-1)
		self%Gc0t(n/2+1:)=self%Gc0t(n/2:1:-1)
		self%SEft=-kp*self%Gc0t*self%GBt
		self%SEBt=-self%Gc0t*self%Gft
		self%SEct=-self%GBt*self%Gft
	end subroutine
end module
program main
	use global
	implicit none
	type(gfs(n)) :: phy,aphy
	integer :: i,j,e
	complex(wp) :: p(10)
	logical :: flag

	call omp_set_nested(.false.)
	!call omp_set_max_active_levels(1)
	call omp_set_dynamic(.false.)
	call mkl_set_dynamic(0)
	!call omp_set_schedule(omp_sched_static,0)

	call mkl_set_num_threads(16)
	call omp_set_num_threads(mkl_get_max_threads())
	
	open(10,file="../data/init.dat")
	open(11,file="../data/t=0.dat")
	open(13,file="TEMPLIST.dat")
	open(22,File="../data/largeN_omega_my.dat")
	open(23,File="../data/largeN_tau_my.dat")
	open(24,File="../data/largeN_T.dat")
	open(25,File="../data/largeN_entropy.dat")
	open(26,File="../data/largeN_fenergy.dat")
	open(27,File="../data/largeN_A0.dat")
	open(233,File="../data/check.dat")

	call set_grid(omega(n/2+1:),[reset,add_log_linear,(add_linear,i=1,size(nadd))],from=[min_omega,fcut-0.1_wp,phicut-0.01_wp,bcut/2],to=[bcut,fcut+0.7_wp/fbeta,phicut+0.7_wp/phibeta,bcut],n=[nlog,nl, nadd])
	omega(1:n/2)=-omega(n:n/2+1:-1)
	write(*,"('n: ',i4,', max freq: ',es10.2,', min freq: ',es10.2)")n,maxval(abs(omega)),minval(abs(omega))

	do 
		call phy%init(10)

		rewind(13)
		write(25,"(A)")"T Sn Sa"
		!write(26,"(A)")"T Fn Fa pT pFn pFa"
		write(26,"(A)")"T F1 F2 F3 F4 F5 F6 pT pF1 pF2 pF3 pF4 pF5 pF6"
		flag=.true.
		do
			if(flag) then
				read(13,*,iostat=e)Tk
				if(Tk-1e-4_wp>eps) then
					flag=.false.
				endif
			else
				Tk=Tk-1e-4_wp
				flag=.true.
			endif
			beta=1d0/Tk
			if(e==-1) then
				write(*,*)"End of Tk!!!"
				exit
			endif
			write(*,"(*(A4,es9.2))")"Tk:  ",Tk," Jk:",Jk," kp:",kp," r: ",r," g: ",g

			tau(1:n/2)=omega(n/2+1:n)/(3._wp*omega(n)-omega(n-1))/Tk
			tau(n/2+1:)=1._wp/Tk-tau(n/2:1:-1)

			call phy%self_consistent(1e-5_wp,[6000,1],2e-2_wp)

			call phy%set_Gfunction([set_SE,set_rSE,set_Gt])
			call aphy%analytic(a_f=&
				0.41_wp&
				!0.03_wp&
				,a_b=&
				0.29_wp&
				!0.67_wp&
			)
			call aphy%set_Gfunction([set_Gt])
			if(flag) then
				write(25,"(*(e28.20))")Tk,phy%entropy(),0._wp!aphy%entropy()
			endif
			!underscore=aphy%fenergy()
			underscore=phy%fenergy()
			!write(26,"(*(e28.20))")Tk,phy%fenergy(),aphy%fenergy()

			!write(*,*)(integrate(omega,A=ff(beta*omega)*imag(phy%Gf)/pi)*integrate(omega,A=ff(beta*omega)*imag(phy%Gc0)/pi))
			!call get_convolution(omega,phy%Gf,phy%Gc0,tau=[1,-1],eta=[-1,-1],rt=phy%GB)
			!write(*,*)integrate(omega,A=fb(beta*omega)*imag(phy%GB)/pi)
			!!write(*,*)integrate(omega,A=fb(beta*omega)*phy%SEB%im)
			!stop

			!call phy%set_Gfunction([set_Gt,set_rSE,set_rG])
			!call phy%export([22],[export_Gw])

			if(.not.phy%conv) then
				exit
			endif
			if(flag) then
				call phy%export([22],[export_Gw])
			endif
			!!write(27,"(*(e28.20))")Tk,-phy%Gf(n/2)%im

			if(is_expinit.and.phy%conv) then
				call phy%export([10],[export_SE])
			endif
		enddo
		stop
		Jk=Jk+0.1_wp
		if(Jk>2.2_wp-1e-4_wp) then
			exit
		endif
	enddo
end program
