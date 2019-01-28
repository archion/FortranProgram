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
	real(wp) :: s=nan,r=nan,kp=nan,fcut=1._wp,fbeta=1._wp,phicut=5e-2_wp,phibeta=200._wp,Jk=nan,g=nan
	real(wp) :: Tk,beta,norm_Ac=1._wp,norm_Aphi=1._wp,bcut=10._wp,min_omega=1e-11_wp!7e-11_wp!1e-11_wp!1e-12_wp
	integer, parameter :: nlog=31,nl=7,nadd(3)=[5,10,25]
	!integer, parameter :: nlog=31,nl=14,nadd(3)=[5,10,20]
	!integer, parameter :: nlog=35,nl=14,nadd(3)=[5,10,25]
	!integer, parameter :: nlog=35,nl=40,nadd(3)=[25,25,75]
	!integer, parameter :: nlog=35,nl=34,nadd(2)=[5,10]
	integer, parameter :: n=(nl*nlog+sum(nadd))*2
	integer, parameter :: set_rGB=1,set_Gt=2,set_rG=3,set_rSE=4,set_SE=5,export_Gw=1,export_Gt=2,export_fe=3,export_entropy=4,export_SE=5
	real(wp) :: omega(n),tau(n)
	logical :: is_expinit=.true.
	type gfs(n)
		integer, len :: n
		real(wp) :: rhoc(n),rhophi(n),lambda
		complex(wp) :: Gf(n),GB(n),Gc(n),Gphi(n),Gc0(n),Gphi0(n)
		complex(wp) :: SEf(n),SEB(n),SEc(n),SEphi(n)
		real(wp) :: Gft(n),GBt(n),Gc0t(n),Gphi0t(n),Gct(n),Gphit(n),SEft(n),SEBt(n),SEct(n),SEphit(n)
		logical :: conv=.false.
		character :: label=" "
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
		!norm_Ac=1._wp/integrate(omega,A=merge(self%rhoc,0._wp,abs(omega)<10._wp))
		norm_Ac=1._wp/integrate(omega,A=self%rhoc)
		norm_Aphi=1._wp/integrate(omega(n/2+1:),A=self%rhoc(n/2+1:))
		self%rhoc=self%rhoc*norm_Ac
		self%rhophi=self%rhophi*norm_Aphi
		self%Gc0%im=-self%rhoc*pi
		self%Gphi0%im=-self%rhophi*pi

		i=1
		rewind(ut)
		read(ut,"(*(e29.20e3))")omg(1),dat(1,:)
		o:	do 
			read(ut,"(*(e29.20e3))",iostat=e)omg(2),dat(2,:)
			if(e==-1) then
				backspace(ut)
				backspace(ut)
				backspace(ut)
				read(ut,"(*(e29.20e3))")omg(1),dat(1,:)
				dat(2,:)=dat(2,:)+(dat(2,:)-dat(1,:))/(omg(2)-omg(1))*(-omega(1)-omg(2))
				omg(2)=-omega(1)
				read(ut,"(*(e29.20e3))")omg(1),dat(1,:)
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
				write(ut(k),"(A)")"T omega rGf iGf rGB iGB rSEf iSEf rSEB iSEB rGc iGc rGphi iGphi rhoc rhophi rSEc iSEc rSEphi iSEphi rGc0 iGc0 rGphi0 iGphi0 lgGf lgGB lgSEf iSEfGf"
				do i=1,n
					write(ut(k),"(*(e29.20e3))")Tk,omega(i),self%Gf(i),self%GB(i),self%SEf(i),self%SEB(i),self%Gc(i),self%Gphi(i),self%rhoc(i),self%rhophi(i),self%SEc(i),self%SEphi(i),self%Gc0(i),self%Gphi0(i),imag(log(-1d0/self%Gf(i))),imag(log(-1d0/self%GB(i))),imag(log(self%SEf(i))),imag(self%Gf(i)*self%SEf(i))
				enddo
				write(ut(k),"(x)")
			case(export_Gt)
				write(ut(k),"(A)")"T tau Gft GBt SEft SEBt Gct Gphit SEct SEphit Gc0t Gphi0t"
				do i=1,n
					write(ut(k),"(*(e29.20e3))")Tk,tau(i),self%Gft(i),self%GBt(i),self%SEft(i),self%SEBt(i),self%Gct(i),self%Gphit(i),self%SEct(i),self%SEphit(i),self%Gc0t(i),self%Gphi0t(i)
				enddo
				write(ut(k),"(x)")
			case(export_SE)
				write(*,*)"export initial...."
				rewind(ut(k))
				is_expinit=.false.
				do i=1,n
					write(ut(k),"(*(e29.20e3))")omega(i),self%SEf(i),self%SEB(i)
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
				self%GB=self%GB-Jk
				call set_realpart(omega,self%Gf)
				!write(*,"(4es12.4)")self%GB(n/2:n/2+1)+Jk,1._wp/(-1._wp/Jk-self%SEB(n/2:n/2+1))+Jk,self%SEB(n/2:n/2+1)
				!call set_realpart(omega,self%Gf)
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
			case(set_rSE)
				call set_realpart(omega,self%SEf)
				call set_realpart(omega,self%SEB)
			case(set_SE)
				call get_convolution_new(omega,A=-self%Gf,B=self%GB+Jk,tau=[1,-1],eta=[-1,1],rt=self%SEc)
				call set_realpart(omega,self%SEc)
				call get_convolution_new(omega,A=self%Gf,B=self%Gf,tau=[1,-1],eta=[-1,-1],rt=self%SEphi)
				call set_realpart(omega,self%SEphi)
			end select
		enddo
	end subroutine
	real(wp) elemental function A_phi(x)
		real(wp), intent(in) :: x
		if(abs(x)<eps) then
			A_phi=0._wp
		elseif(abs(x)<phicut) then
			A_phi=norm_Aphi*abs(x)**(1._wp-s)*sign(1._wp,x)
		else
			A_phi=norm_Aphi*phicut**(1._wp-s)*sign(1._wp,x)*1._wp/((exp(phibeta*(x-phicut-4._wp/phibeta))+1._wp)*(exp(phibeta*(-x-phicut-4._wp/phibeta))+1._wp))
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
	!subroutine get_convolution(x,f1,f2,rt_x,rt)
		!real(wp) :: x(:),x_rt(:)
		!complex(wp) :: rt(:),f1,f2
		!external f1,f2
		!integer :: i,j,k,n
		!do k=1,size(rt_x)
			!do 
			!enddo
		!enddo
	!end subroutine
	subroutine get_convolution_new(omega,A,B,tau,eta,rt,noreset,d)
		real(wp) :: omega(:)
		complex(wp) :: rt(:)
		complex(wp) :: A(:),B(:)
		integer :: tau(2),eta(2)
		logical, optional :: noreset,d
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
				if(present(d)) then
					if(abs(x)<eps) then
						y=1._wp
					else
						y=x*(0.5_wp*beta/sinh(0.5_wp*x*beta))**2
					endif
				else
					y=fb(beta*x)
				endif
			case(-1)
				if(present(d)) then
					y=x*(0.5_wp*beta/cosh(0.5_wp*x*beta))**2
				else
					y=ff(beta*x)
				endif
			end select
		end function
		real(wp) function f2(x) result(y)
			real(wp) :: x
			select case(eta(2))
			case(1)
				if(present(d)) then
					if(abs(x)<eps) then
						y=1._wp
					else
						y=x*(0.5_wp*beta/sinh(0.5_wp*x*beta))**2
					endif
				else
					y=fb(beta*x)
				endif
			case(-1)
				if(present(d)) then
					y=x*(0.5_wp*beta/cosh(0.5_wp*x*beta))**2
				else
					y=ff(beta*x)
				endif
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
		real(wp) :: x(n*2*2),x_(n*2*2),nf(1),tmp(1),fe
		logical :: flag
		integer :: i1,i2,k,i,j,ei
		fe=0d0
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
				self%GB=Jk/(-1._wp-Jk*self%SEB)
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

				self%SEf=0._wp
				call get_convolution_new(omega,A=-kp*self%Gc0,B=self%GB,tau=[1,1],eta=[-1,1],rt=self%SEf,noreset=underscore)
				call get_convolution_new(omega,A=self%Gc0,B=self%Gf,tau=[1,-1],eta=[-1,-1],rt=self%SEB)
				if(abs(g)>0._wp) then
					call get_convolution_new(omega,A=-2d0*g**2*self%Gphi0,B=self%Gf,tau=[1,1],eta=[1,-1],rt=self%SEf,noreset=underscore)
				endif


				x_=[self%SEf%im,self%SEB%im,self%SEf%re,self%SEB%re]

				!call self%set_Gfunction([set_SE,set_rSE,set_rG])
				!write(*,*)i1,maxval(abs(x(:n*2)-x_(:n*2))/(abs(x_(:n*2))+eps)),self%fenergy()
				!if(self%fenergy()>fe) then
					!write(*,*)"check"
					!read(*,*)
				!else
					!fe=self%fenergy()
				!endif

				if(all(abs(x(:n*2)-x_(:n*2))/(abs(x_(:n*2))+eps)<tol)) then
				!if(all(abs(x(:n*2)-x_(:n*2))<tol)) then
				!if(all(abs(x(:n*2)-x_(:n*2))/(min(abs(x_(:n*2))+eps,1._wp))<tol)) then
					!if(ei==0) ei=1
					!if(ei>20) then
						write(*,*)"converged"


						if(abs(0.5_wp-nf(1))>1e-2_wp) then
							is_expinit=.false.
							write(*,*)"nf is not correct****************************",nf(1)
							exit
						else
							call set_realpart(omega,self%SEB)
							if(self%SEB(n/2)%re<-1._wp/Jk) then
								write(*,*)"warning from SEB", -1._wp/Jk,self%SEB(n/2)%re
							endif
							write(*,"(i6,es12.4,',nf: ',es12.4,',SEB0: ',es12.4)")i1,maxval(abs(x-x_)/(abs(x_)+eps)),nf*2._wp,-1._wp/Jk-self%SEB(n/2)%re
							self%conv=.true.
						endif
						exit
					!endif
				endif
				if(i1==niter(1)) then
					write(*,*)"not converged",maxval(abs(x-x_)/(abs(x_)+eps)),"***************"
					is_expinit=.false.
					exit
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
			enddo

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
		real(wp) :: df(n),db(n),kp_
		complex(wp) :: ctmp1(n),ctmp2(n),ctmp3(n),ctmp4(n)
		integer :: i
		df=omega*(0.5_wp*beta/cosh(0.5_wp*omega*beta))**2
		db=omega*(0.5_wp*beta/sinh(0.5_wp*omega*beta))**2

		if(self%label=="a") then
			!rt=-(integrate(omega,A=&
				!ff(beta*omega)*(imag(log(-1d0/self%Gf)))&
				!+kp*fb(beta*omega)*(imag(log(-1d0/self%GB)))&
				!)/pi*beta&
				!)
			rt=-(integrate(omega,A=&
				df*(imag(log(-1d0/self%Gf))+self%Gf%re*imag(omega-1d0/self%Gf)+0.5d0*(1d0+sign(1d0,omega))*pi)&
				+kp*db*(imag(log(-1d0/self%GB))+(self%GB%re)*imag(-1d0/Jk-1d0/self%GB))&
				)/pi&
				-log(2d0)&
				)
		else
			rt=-(integrate(omega,A=&
				df*(imag(log(-1d0/self%Gf))+self%Gf%re*imag(omega-1d0/self%Gf)+0.5d0*(1d0+sign(1d0,omega))*pi)&
				+kp*db*(imag(log(-1d0/self%GB))+(self%GB%re)*imag(-1d0/Jk-1d0/self%GB))&
				-kp*df*(self%Gc0%im*(self%SEc%re-Jk/2d0))&
				-db*(g**2*self%Gphi0%im*self%SEphi%re)&
				)/pi&
				-log(2d0)&
				)
		endif

		!call get_convolution_new(omega,A=-self%Gf,B=self%GB+Jk,tau=[1,-1],eta=[-1,1],rt=ctmp1,d=underscore)
		!call set_realpart(omega,ctmp1)
		!call get_convolution_new(omega,A=self%Gf,B=self%Gf,tau=[1,-1],eta=[-1,-1],rt=ctmp2,d=underscore)
		!call set_realpart(omega,ctmp2)
		!rt=-(integrate(omega,A=&
			!df*(imag(log(-1d0/self%Gf)+self%Gf*(omega-1d0/self%Gf))+0.5d0*(1d0+sign(1d0,omega))*pi)&
			!+kp*db*(imag(log(-1d0/self%GB)+self%GB*(-1d0/Jk-1d0/self%GB)))&
			!-kp*df*imag(self%Gc0*self%SEc)&
			!-kp*ff(beta*omega)*imag(self%Gc0*ctmp1)&
			!-db*imag(g**2*self%Gphi0*self%SEphi)&
			!-fb(beta*omega)*imag(g**2*self%Gphi0*ctmp2)&
			!)/pi&
			!-log(2d0)&
			!)

		!call get_convolution_new(omega,A=self%Gc0,B=self%Gf,tau=[1,-1],eta=[-1,-1],rt=ctmp1,d=underscore)
		!call set_realpart(omega,ctmp1)
		!call get_convolution_new(omega,A=-g**2*self%Gphi0,B=self%Gf,tau=[1,1],eta=[1,-1],rt=ctmp2,d=underscore)
		!call set_realpart(omega,ctmp2)
		!call get_convolution_new(omega,A=-g**2*self%Gphi0,B=self%Gf,tau=[1,1],eta=[1,-1],rt=ctmp3)
		!call set_realpart(omega,ctmp3)
		!rt=-(integrate(omega,A=&
			!df*(imag(log(-1d0/self%Gf)+self%Gf*(omega-1d0/self%Gf))+0.5d0*(1d0+sign(1d0,omega))*pi)&
			!+kp*db*(imag(log(-1d0/self%GB)))&
			!-kp*fb(beta*omega)*imag(self%GB*ctmp1)&
			!-ff(beta*omega)*imag(self%Gf*ctmp2)&
			!-df*imag(self%Gf*ctmp3)&
			!)/pi&
			!-log(2d0)&
			!)

		!call get_convolution_new(omega,A=self%Gc0,B=self%Gf,tau=[1,-1],eta=[-1,-1],rt=ctmp1,d=underscore)
		!call set_realpart(omega,ctmp1)
		!call get_convolution_new(omega,A=self%Gf,B=self%Gphi0,tau=[1,1],eta=[-1,1],rt=ctmp2,d=underscore)
		!call set_realpart(omega,ctmp2)
		!ctmp3=0._wp
		!call get_convolution_new(omega,A=-kp*self%Gc0,B=self%GB,tau=[1,1],eta=[-1,1],rt=ctmp3)
		!call set_realpart(omega,ctmp3)
		!ctmp4=0._wp
		!call get_convolution_new(omega,A=-g**2*self%Gphi0,B=self%Gf,tau=[1,1],eta=[1,-1],rt=ctmp4)
		!call set_realpart(omega,ctmp4)
		!write(25,"(e28.20$)")Tk,Jk,g,r,s,kp,rt
		!rt=-(integrate(omega,A=&
			!df*(imag(log(-1d0/self%Gf))+imag(self%Gf*(ctmp3+ctmp4))+0.5d0*(1d0+sign(1d0,omega))*pi)&
			!+kp*db*(imag(log(-1d0/self%GB)))&
			!-kp*fb(beta*omega)*(imag(self%GB*ctmp1))&
			!+g**2*ff(beta*omega)*(imag(self%Gf*ctmp2))&
			!)/pi&
			!-log(2d0)&
			!)
		!write(25,"(e28.20)")rt

		!rt=-(integrate(omega,A=&
			!df*(self%Gf%im*real(ctmp3+ctmp4,kind=wp))&
			!-kp*fb(beta*omega)*(imag(self%GB*ctmp1))&
			!+g**2*ff(beta*omega)*(imag(self%Gf*ctmp2))&
			!)/pi&
			!)&
		!+(integrate(omega,A=&
			!df*(self%Gf%re*imag(ctmp4))&
			!+kp*db*((self%GB%re)*imag(-1d0/Jk-1d0/self%GB))&
			!-kp*df*(self%Gc0%im*(self%SEc%re-Jk/2d0))&
			!-db*(g**2*self%Gphi0%im*self%SEphi%re)&
			!)/pi&
			!)
		!write(25,"(*(e28.20))")Tk,Jk,g,r,s,kp,&
		!integrate(omega,A=kp*fb(beta*omega)*imag(self%GB*ctmp1)),&
		!-integrate(omega,A=kp*db*self%GB%re*(self%SEB%im)-kp*df*(self%Gc0%im*(self%SEc%re-Jk/2d0))-df*self%SEf%re*self%Gf%im)
		!integrate(omega,A=ff(beta*omega)*imag(1._wp/omega*ctmp3)),&
		!-integrate(omega,A=df*imag(1._wp/omega*self%SEf))
		!integrate(omega,A=fb(beta*omega)*imag(ctmp1)),&
		!-integrate(omega,A=db*imag(self%SEB))
		!integrate(omega,A=ff(beta*omega)*(imag(ctmp3))),&
		!integrate(omega,A=df*(imag(self%SEf)))
		!integrate(omega,A=df*(self%Gf%im*real(ctmp3,kind=wp))),&
		!integrate(omega,A=df*(self%Gf%im*real(ctmp4,kind=wp))),&
		!integrate(omega,A=-kp*fb(beta*omega)*(imag(self%GB*ctmp1))),&
		!integrate(omega,A=+g**2*ff(beta*omega)*(imag(self%Gf*ctmp2))),&
		!integrate(omega,A=df*(self%Gf%re*imag(ctmp4))),&
		!integrate(omega,A=+kp*db*((self%GB%re)*imag(-1d0/Jk-1d0/self%GB))),&
		!integrate(omega,A=-kp*df*(self%Gc0%im*(self%SEc%re-Jk/2d0))),&
		!integrate(omega,A=-db*(g**2*self%Gphi0%im*self%SEphi%re))
			
	end function
	real(wp) function fenergy(self) result(rt)
		class(gfs(*)) :: self
		integer :: i
		integer :: cs

		if(self%label=="a") then
			! lnG
			rt=(integrate(omega,A=&
				ff(beta*omega)*(imag(log(-1._wp/self%Gf))+0.5_wp*(1._wp+sign(1._wp,omega))*pi)&
				+kp*fb(beta*omega)*imag(log(-1._wp/self%GB))&
				)/pi&
				-1._wp/beta*log(2._wp)&
				)
		else
			!! large-N
			!rt=(integrate(omega,A=&
				!ff(beta*omega)*(imag(log(-1._wp/self%Gf))+0.5_wp*(1._wp+sign(1._wp,omega))*pi)&
				!+kp*fb(beta*omega)*imag(log(-1._wp/self%GB)+(self%GB+Jk)*(-1._wp/Jk-1._wp/self%GB))&
				!)/pi&
				!-1._wp/beta*log(2._wp)&
				!+integrate(omega,A=fb(beta*omega)*imag(g**2*self%Gphi0*self%SEphi)/pi)&
				!)

			! B-S functional
			rt=(&
				!integrate(omega,A=&
				!ff(beta*omega)*(imag(log(-1._wp/self%Gf)+self%Gf*(omega-1._wp/self%Gf))+0.5_wp*(1._wp+sign(1._wp,omega))*pi)&
				!+kp*fb(beta*omega)*imag(log(-1._wp/self%GB)+(self%GB+Jk)*(-1._wp/Jk-1._wp/self%GB))&
				!)/pi&
				!-1._wp/beta*log(2._wp)&
				-integrate(omega,A=ff(beta*omega)*imag(self%Gc0*self%SEc)/pi*kp)&
				-integrate(omega,A=fb(beta*omega)*imag(g**2*self%Gphi0*self%SEphi)/pi)&
				)
		endif


		!! Sigma G
		!rt=(integrate(omega,A=&
			!ff(beta*omega)*(imag(self%Gf*(omega-1._wp/self%Gf)))&
			!+kp*fb(beta*omega)*imag((self%GB+Jk)*(-1._wp/Jk-1._wp/self%GB))&
			!)/pi&
			!)

		!! LM
		!rt=(integrate(omega,A=&
			!ff(beta*omega)*(imag(log(-1._wp/self%Gf))+0.5_wp*(1._wp+sign(1._wp,omega))*pi)&
			!)/pi&
			!-1._wp/beta*log(2._wp)&
			!!+integrate(omega,A=fb(beta*omega)*imag(g**2*self%Gphi0*self%SEphi)/pi)&
			!!+integrate(omega,A=0.5_wp*ff(beta*omega)*imag(self%Gf*self%SEf)/pi)&
			!)
		!write(26,"(*(e29.20e3))")Tk&
			!!,integrate(omega,A=ff(beta*omega)*(imag(log(-1._wp/self%Gf))))/pi&
			!!,integrate(omega,A=kp*fb(beta*omega)*(imag(log(-1._wp/self%GB))))/pi&
			!,integrate(omega,A=ff(beta*omega)*imag(self%Gf*self%SEf))&
			!,integrate(omega,A=ff(beta*omega)*imag(self%Gf*omega))
			!!,integrate(omega,A=fb(beta*omega)*imag(g**2*self%Gphi0*self%SEphi))
			!rt=integrate(omega,A=ff(beta*omega)*(imag(log(-1._wp/self%Gf))))/pi+integrate(omega,A=kp*fb(beta*omega)*(imag(log(-1._wp/self%GB))))/pi
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
		real(wp) :: M_f,M_B,M_c,M_phi,t,sh,ch,lsh,lch,lgam_a,lgam_b,lgam_sa(2),lgam_sb,lgam_c,lgam_phi,SEB0,w
		integer :: i,j
		!M_f=0.724_wp
		M_f=0.59_wp
		M_B=1._wp
		M_c=1._wp
		!M_phi=0.85_wp
		M_c=norm_Ac
		M_phi=norm_Aphi
		t=1._wp
		if(abs(1._wp-r-a_f-a_b)<1e-8_wp) then
			!t=((r+a_f)*tan(pi*(r+a_f)/2._wp)/(gamma(r+1)*pi*2*M_f*M_B*M_c))**(1._wp/(1-r))
			M_B=(r+a_f)*tan(pi*(r+a_f)/2._wp)/(gamma(r+1)*pi*2*M_f*M_c*t**(1-r))
			SEB0=-1._wp/Jk
		else
			!t=((1._wp-0.5_wp*s)*tan(pi*s/4._wp)/(4._wp*pi*M_phi*M_f**2*g**2*gamma(2-s)))**(0.5_wp/a_f)
			M_f=((1._wp-0.5_wp*s)*tan(pi*s/4._wp)/(4._wp*pi*M_phi*t**(2._wp*a_f)*g**2*gamma(2-s)))**(0.5_wp)
			M_B=gamma(1._wp+r)*Jk**2*M_c*M_f/(t**(1._wp+r))
			SEB0=0._wp
		endif
		do i=1,n
			w=omega(i)!/2.1_wp
			if(exp(beta*abs(w))>1e100_wp) then
				lsh=beta*abs(w)/2._wp-log(2._wp)
				lch=beta*abs(w)/2._wp-log(2._wp)
				lgam_a=log(2._wp*pi)-beta*abs(w)/2._wp+(a_f-1._wp)*(log(beta*abs(w))-log(2._wp*pi))
				lgam_b=log(2._wp*pi)-beta*abs(w)/2._wp+(a_b-1._wp)*(log(beta*abs(w))-log(2._wp*pi))
				lgam_sa(1)=log(2._wp*pi)-beta*abs(w)/2._wp+(r+a_B)*(log(beta*abs(w))-log(2._wp*pi))
				lgam_sa(2)=log(2._wp*pi)-beta*abs(w)/2._wp+(1._wp-s+a_f)*(log(beta*abs(w))-log(2._wp*pi))
				lgam_sb=log(2._wp*pi)-beta*abs(w)/2._wp+(r+a_f)*(log(beta*abs(w))-log(2._wp*pi))
				lgam_c=log(2._wp*pi)-beta*abs(w)/2._wp+r*(log(beta*abs(w))-log(2._wp*pi))
				lgam_phi=log(2._wp*pi)-beta*abs(w)/2._wp+(1._wp-s)*(log(beta*abs(w))-log(2._wp*pi))
			else
				lsh=log(sinh(beta*abs(w)/2._wp))
				lch=log(cosh(beta*abs(w)/2._wp))
				lgam_a=real(lgamma(a_f/2._wp+img*beta*w/(2._wp*pi))+lgamma(a_f/2._wp-img*beta*w/(2._wp*pi)))
				lgam_b=real(lgamma(a_b/2._wp+img*beta*w/(2._wp*pi))+lgamma(a_b/2._wp-img*beta*w/(2._wp*pi)))
				lgam_sa(1)=real(lgamma((1._wp+r+a_B)/2._wp+img*beta*w/(2._wp*pi))+lgamma((1._wp+r+a_B)/2._wp-img*beta*w/(2._wp*pi)))
				lgam_sa(2)=real(lgamma((2._wp-s+a_f)/2._wp+img*beta*w/(2._wp*pi))+lgamma((2._wp-s+a_f)/2._wp-img*beta*w/(2._wp*pi)))
				lgam_sb=real(lgamma((1._wp+r+a_f)/2._wp+img*beta*w/(2._wp*pi))+lgamma((1._wp+r+a_f)/2._wp-img*beta*w/(2._wp*pi)))
				lgam_c=real(lgamma((r+1._wp)/2._wp+img*beta*w/(2._wp*pi))+lgamma((r+1._wp)/2._wp-img*beta*w/(2._wp*pi)))
				lgam_phi=real(lgamma((2._wp-s)/2._wp+img*beta*w/(2._wp*pi))+lgamma((2._wp-s)/2._wp-img*beta*w/(2._wp*pi)))
			endif
			sh=real(exp((a_f-1._wp)*(log(2._wp*pi)-log(beta))+lgam_a+lsh))*sign(1._wp,w)
			ch=real(exp((a_f-1._wp)*(log(2._wp*pi)-log(beta))+lgam_a+lch))
			self%Gf(i)=-M_f*t**a_f/gamma(a_f)*cmplx(-1d0/tan(pi*a_f/2._wp)*sh,ch,kind=wp)
			sh=real(exp(r*(log(2._wp*pi)-log(beta))+lgam_c+lsh))*sign(1._wp,w)
			ch=real(exp(r*(log(2._wp*pi)-log(beta))+lgam_c+lch))
			self%Gc0(i)=-M_c*cmplx(-1d0/tan(pi*(r+1._wp)/2._wp)*sh,ch,kind=wp)
			sh=real(exp((a_b-1._wp)*(log(2._wp*pi)-log(beta))+lgam_b+lsh))*sign(1._wp,w)
			ch=real(exp((a_b-1._wp)*(log(2._wp*pi)-log(beta))+lgam_b+lch))
			self%GB(i)=-M_B*t**a_b/gamma(a_b)*cmplx(tan(pi*a_b/2._wp)*ch,sh,kind=wp)-Jk
			sh=real(exp((1._wp-s)*(log(2._wp*pi)-log(beta))+lgam_phi+lsh))*sign(1._wp,w)
			ch=real(exp((1._wp-s)*(log(2._wp*pi)-log(beta))+lgam_phi+lch))
			self%Gphi0(i)=-M_phi*cmplx(tan(pi*(2._wp-s)/2._wp)*ch,sh,kind=wp)
			sh=real(exp((r+a_B)*(log(2._wp*pi)-log(beta))+lgam_sa(1)+lsh))*sign(1._wp,w)
			ch=real(exp((r+a_B)*(log(2._wp*pi)-log(beta))+lgam_sa(1)+lch))
			self%SEf(i)=0._wp
			self%SEf(i)=self%SEf(i)-kp*M_c*M_B*gamma(1._wp+r)*t**a_B/gamma(1._wp+r+a_B)*cmplx(-1d0/tan(pi*(1._wp+r+a_B)/2._wp)*sh,ch,kind=wp)
			sh=real(exp((1._wp-s+a_f)*(log(2._wp*pi)-log(beta))+lgam_sa(2)+lsh))*sign(1._wp,w)
			ch=real(exp((1._wp-s+a_f)*(log(2._wp*pi)-log(beta))+lgam_sa(2)+lch))
			self%SEf(i)=self%SEf(i)-2._wp*g**2*M_phi*M_f*gamma(2._wp-s)*t**a_f/gamma(2._wp-s+a_f)*cmplx(-1d0/tan(pi*(2._wp-s+a_f)/2._wp)*sh,ch,kind=wp)
			sh=real(exp((r+a_f)*(log(2._wp*pi)-log(beta))+lgam_sb+lsh))*sign(1._wp,w)
			ch=real(exp((r+a_f)*(log(2._wp*pi)-log(beta))+lgam_sb+lch))
			self%SEB(i)=-M_c*M_f*gamma(1._wp+r)*t**a_f/gamma(1._wp+r+a_f)*cmplx(tan(pi*(1._wp+r+a_f)/2._wp)*ch,sh,kind=wp)+SEB0
			if(abs(w)<fcut) then
				!self%Gc0(i)%im=-pi*M_c*(abs(w)/fcut)**r
			elseif(abs(w)<=(-omega(1))) then
				self%Gc0(i)%im=-pi*M_c*exp(-(w**2-fcut**2)/pi*fbeta)
			else
				self%Gc0(i)%im=0._wp
			endif
			if(abs(w)<eps) then
				!self%Gphi0(i)%im=0._wp
			elseif(abs(w)<phicut) then
				!self%Gphi0(i)%im=-pi*M_phi*abs(w)**(1._wp-s)*sign(1._wp,w)
			else
				self%Gphi0(i)%im=-pi*M_phi*phicut**(1._wp-s)*sign(1._wp,w)*1._wp/((exp(phibeta*(w-phicut-4._wp/phibeta))+1._wp)*(exp(phibeta*(-w-phicut-4._wp/phibeta))+1._wp))
			endif
		enddo
		call set_realpart(omega,self%Gc0)
		call set_realpart(omega,self%Gphi0)

		!do i=n/2+1,n
			!if(abs(omega(i))<fcut) then
				!j=i
			!else
				!self%Gf(i)%im=-abs(self%Gf(j)%im)*exp(-(omega(i)**2-omega(j)**2)/pi*fbeta*8._wp)
			!endif
			!self%Gf(n-i+1)%im=self%Gf(i)%im
		!enddo
		!call set_realpart(omega,self%Gf)
		!self%SEf=omega-1._wp/self%Gf
		!call self%export([22,23],[export_Gw])
		!self%SEf=0._wp
		!call get_convolution_new(omega,A=-kp*self%Gc0,B=self%GB,tau=[1,1],eta=[-1,1],rt=self%SEf,noreset=underscore)
		!call get_convolution_new(omega,A=self%Gc0,B=self%Gf,tau=[1,-1],eta=[-1,-1],rt=self%SEB)
		!if(abs(g)>0._wp) then
			!call get_convolution_new(omega,A=-2d0*g**2*self%Gphi0,B=self%Gf,tau=[1,1],eta=[1,-1],rt=self%SEf,noreset=underscore)
		!endif
		!call set_realpart(omega,self%SEf)
		!call set_realpart(omega,self%Gf)
		!self%SEB=-1._wp/Jk-1._wp/self%GB
		!call get_convolution(omega,A=-self%Gf,B=self%GB+Jk,tau=[1,-1],eta=[-1,1],rt=self%SEc)
		!call get_convolution(omega,A=-kp*self%Gc0,B=self%GB+Jk,tau=[1,1],eta=[-1,1],rt=self%SEf)
		!call get_convolution(omega,A=self%Gc0,B=self%Gf,tau=[1,-1],eta=[-1,-1],rt=self%SEB)
		!call set_realpart(omega,self%SEc)
		!call set_realpart(omega,self%SEf)
		!call set_realpart(omega,self%SEB)
		!call set_realpart(omega,self%Gc0)

		!rewind(22)
		!call self%export([22],[export_Gw])
		!!write(*,*)self%Gf(n)%im,j,omega(j)
		!!stop
		!call set_realpart(omega,self%Gf)
		!call set_realpart(omega,self%GB)
		!call self%export([22],[export_Gw])
		!stop

		!self%Gft=-M_f*(pi*t/((sin(pi*tau/beta))*beta))**a_f
		!self%GBt=-M_B*(pi*t/((sin(pi*tau/beta))*beta))**a_b
		!self%Gc0t=-M_c*(t/tau)**(r+1._wp)
		!self%Gft(n/2+1:)=self%Gft(n/2:1:-1)
		!self%GBt(n/2+1:)=self%GBt(n/2:1:-1)
		!self%Gc0t(n/2+1:)=self%Gc0t(n/2:1:-1)
		!self%SEft=-kp*self%Gc0t*self%GBt
		!self%SEBt=-self%Gc0t*self%Gft
		!self%SEct=-self%GBt*self%Gft
	end subroutine
	subroutine xscale(A,x,scal)
		complex(wp) :: A(:)
		real(wp) :: x(:),scal
		complex(wp) :: tmp(size(A))
		integer :: l,i
		do i=1,size(A)
			l=min(max(find(x,x(i)*scal),1),size(A)-2)
			tmp(i)=A(l)+(A(l)-A(l+1))/(x(l)-x(l+1))*(x(i)*scal-x(l))
		enddo
		A=tmp
	end subroutine
	real(wp) function afn() result(rt)
		integer :: i
		integer :: N
		N=1000000
		rt=0._wp
		do i=1,N/2
			rt=rt+log(sin(pi*(N+1-i)/(N+N*kp)))-log(sin(pi*i/(N+N*kp)))
		enddo
		rt=rt/N
	end function
end module
program main
	use global
	implicit none
	type(gfs(n)) :: phy,aphy,mphy
	integer :: i,j,e
	complex(wp) :: p(10)
	logical :: flag,last
	real(wp) :: a_f,a_b,tmp(1)
	phy%label="n"
	aphy%label="a"
	mphy%label="m"

	call omp_set_nested(.false.)
	!call omp_set_max_active_levels(1)
	call omp_set_dynamic(.false.)
	call mkl_set_dynamic(0)
	!call omp_set_schedule(omp_sched_static,0)

	call mkl_set_num_threads(4)
	call omp_set_num_threads(mkl_get_max_threads())
	
	open(10,file="../data/init.dat")
	open(11,file="../data/t=0.dat")
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

	write(25,"(A)")"T J g r s kappa Sn Sa GB g f B c phi"
	Jk=1.8_wp
	g=0.0_wp
	r=0.3_wp
	!r=-0.23_wp
	s=0.3_wp
	kp=0.3_wp
	do 
		!Jk=for_in(x=[0.5_wp,0.978_wp,1.8_wp,1.8_wp,1.153_wp,0.5_wp],id=1)
		!g=for_in(x=[0._wp,0._wp,0._wp,0.2_wp,0.2_wp,0.2_wp],id=2)
		!Jk=for_in(x=[0.978_wp,1.8_wp,1.153_wp,0.5_wp],id=1)
		!g=for_in(x=[0._wp,0._wp,0.2_wp,0.2_wp],id=2)
		!Jk=for_in(x=[1.8_wp,0.5_wp],id=1)
		!g=for_in(x=[0._wp],id=2)
		Jk=for_in(x=[1.8_wp],id=1)
		g=for_in(x=[0._wp],id=2)
		!Jk=for_in(x=[1._wp],id=1)
		!Jk=for_in(x=[&
			![0.5_wp,0.6_wp,0.7_wp,0.8_wp,0.9_wp,0.93_wp,0.96_wp,0.97_wp,0.98_wp,0.99_wp,1._wp,1.05_wp,1.2_wp,1.3_wp,1.4_wp,1.5_wp,1.6_wp,1.7_wp,1.8_wp,1.9_wp,2._wp]&
			!![0.1_wp,0.5_wp,0.7_wp,0.9_wp,1.05_wp,1.1_wp,1.15_wp,1.16_wp,1.19_wp,1.3_wp,1.4_wp,1.5_wp,1.7_wp,2._wp]&
			!![0.1_wp,0.5_wp,0.7_wp,0.8_wp,0.9_wp,0.93_wp,0.96_wp,0.97_wp,0.98_wp,0.99_wp,1._wp,1.05_wp,1.2_wp,1.3_wp,1.4_wp,1.5_wp,1.6_wp,1.7_wp,1.8_wp,1.9_wp,2._wp]&
			!],id=1)
		!g=for_in(x=[&
			!(0._wp,i=1,21)&
			!!(0.02_wp,i=1,21)&
			!],id=2)
		!Jk=for_in(x=[0.1_wp],id=1)
		!g=for_in(x=[0.2_wp],id=2)
		!a_f=for_in([0.001_wp,0.04_wp,0.41_wp,0.41_wp,0.5_wp*s,0.5_wp*s],id=4)
		!a_b=for_in([1._wp-r-0.001_wp,1._wp-r-0.04_wp,1._wp-r-0.41_wp,1._wp-r-0.41_wp,1._wp-r-0.5_wp*s,1._wp+r+0.5_wp*s],id=5)
		!a_f=for_in([0.04_wp,0.41_wp,0.5_wp*s,0.5_wp*s],id=4)
		!a_b=for_in([1._wp-r-0.04_wp,1._wp-r-0.41_wp,1._wp-r-0.5_wp*s,1._wp+r+0.5_wp*s],id=5)
		!a_f=for_in([0.5_wp*s],id=4)
		!a_b=for_in([1._wp+r+0.5_wp*s],id=5)
		!a_f=for_in([1._wp/(1._wp+kp)],id=4)
		!a_b=for_in([1._wp-1._wp/(1._wp+kp)],id=5)
		a_f=for_in([0.413_wp],id=4)
		a_b=for_in([1._wp-r-0.413_wp],id=5)
		!a_f=for_in([0.04_wp],id=4)
		!a_b=for_in([1._wp-r-0.04_wp],id=5)
		!Jk=for_in(x=[0.5_wp,0.978_wp],id=1)
		!g=for_in(x=[0.0_wp,0.0_wp],id=2)
		!kp=for_in(x=[(i*0.1_wp,i=15,1,-1)],id=1)
		if(any(isnan([Jk,g,kp,r,s,a_f,a_b]))) then
			exit
		endif
		call phy%init(10)

		!write(26,"(A)")"T Fn Fa pT pFn pFa"
		write(26,"(A)")"T F1 F2 pT pF1 pF2"! F2 F3 F4 F5 F6 pT pF1 pF2 pF3 pF4 pF5 pF6"
		flag=.true.
		Tk=for_in(id=3)
		do
			if(flag) then
				Tk=for_in(x=[&
					 [90:10:-5]*1e-3_wp&
					,[90:10:-5]*1e-4_wp&
					,[90:10:-5]*1e-5_wp&
					,[90:10:-5]*1e-6_wp&
					,[90:10:-5]*1e-7_wp&
					,[90:10:-5]*1e-8_wp&
					,[90:10:-5]*1e-9_wp&
					,[90:10:-5]*1e-10_wp&
					!,[96:10:-2]*1e-7_wp&
					!,[96:10:-2]*1e-8_wp&
					!,[96:10:-2]*1e-9_wp&
					!,[96:10:-2]*1e-10_wp&
					!,[96:10:-2]*1e-11_wp&
					!,[96:10:-2]*1e-12_wp&
					],id=3)
				if(isnan(Tk)) then
					write(*,*)"End of T!!!"
					exit
				endif
				!if(Tk-5e-5_wp>eps) then
					!flag=.false.
				!else
					!exit
				!endif

				!Tk=1e-7_wp
				!s=for_in([1e-2_wp,(i*0.1_wp,i=1,19),2._wp-1e-2_wp],id=3)
				!r=for_in([0._wp,0.05_wp,0.1_wp,0.15_wp,0.2_wp,0.25_wp,0.3_wp,0.35_wp,0.4_wp],id=3)
				!a_f=for_in([0.00001_wp,0.0008_wp,0.003_wp,0.0076_wp,0.014_wp,0.0247_wp,0.04_wp,0.064_wp,0.116_wp],id=4)
				!a_f=for_in([0.769_wp,0.715_wp,0.66_wp,0.603_wp,0.543_wp,0.481_wp,0.41_wp,0.339_wp,0.237_wp],id=4)
			else
				Tk=Tk-5e-5_wp
				flag=.true.
			endif
			beta=1d0/Tk
			if(any(isnan([Jk,g,kp,r,s,Tk]))) then
				exit
			endif

			tau(1:n/2)=omega(n/2+1:n)/(3._wp*omega(n)-omega(n-1))/Tk
			tau(n/2+1:)=1._wp/Tk-tau(n/2:1:-1)


			write(*,"(*(A4,es9.2))")"T:  ",Tk," Jk:",Jk," kp:",kp," r: ",r," g: ",g
			!call phy%set_Gfunction([set_Gt])
			!call phy%export([22],[export_Gw])

			call aphy%analytic(a_f=a_f,a_b=a_b)
			phy%Gc0=aphy%Gc0
			phy%Gphi0=aphy%Gphi0

			call phy%self_consistent(1e-5_wp,[500,1],3e-3_wp)
			if(.not.phy%conv) then
				write(*,*)"not converge!"
				exit
			else
				!rewind(23)
				!rewind(22)
			endif
			call phy%set_Gfunction([set_SE,set_rSE,set_Gt,set_rG])
			call phy%export([22,23],[export_Gw,export_Gt])
			if(is_expinit) then
				call phy%export([10],[export_SE])
				is_expinit=.false.
			endif

			call aphy%analytic(a_f=a_f,a_b=a_b)
			call aphy%set_Gfunction([set_SE])
			call aphy%export([22],[export_Gw])

			!mphy=aphy
			!call xscale(mphy%Gf(n/2+1:),omega(n/2+1:),0.555_wp)
			!mphy%Gf(1:n/2)=-conjg(mphy%Gf(n:n/2+1:-1))
			!mphy%SEf=omega-1._wp/mphy%Gf
			!call mphy%set_Gfunction([set_SE])
			!call mphy%export([22,23],[export_Gw,export_Gt])
			!!!!call xscale(phy%GB(n/2+1:),omega(n/2+1:),1.02_wp)
			!!!call xscale(phy%GB(n/2+1:),omega(n/2+1:),1._wp)
			!!!phy%GB(1:n/2)=conjg(phy%GB(n:n/2+1:-1))

			!write(*,*)phy%fenergy()
			!write(*,*)aphy%fenergy()
			write(26,"(*(e29.20e3))")Tk,phy%fenergy(),aphy%fenergy()!,mphy%fenergy(),phy%entropy()
			!if((Tk-1e-8_wp)<1e-13_wp) then
			write(25,"(*(e29.20e3))")Tk,Jk,g,r,s,kp,phy%entropy(),aphy%entropy()!,mphy%entropy()
			!endif
			!write(25,"(*(e29.20e3))")Tk,Jk,g,r,s,kp,afn()
			!exit

		enddo
		if(flag) then
			write(26,"(x)")
			write(25,"(x)")
			write(22,"(x)")
			write(23,"(x)")
		endif
	enddo
end program
