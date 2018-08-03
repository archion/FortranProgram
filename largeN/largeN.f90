include "fn.f90"
module global
	use M_fn
	use M_const
	use M_matrix
	use mkl_service
	implicit none
	real(8), parameter :: onepi=1d0/pi,halfpi=onepi*0.5d0
	real(8) :: alpha=0.8d0,r=0.3d0,kp=0.3d0,fcut=1d0,fbeta=1d0,phicut=5d-2,phibeta=200d0,Jk=1.8d0,g=0d0
	real(8) :: Tk,beta,norm_Ac=1d0,norm_Aphi=1d0,base=2.3d0,bcut=10d0
	integer, parameter :: nlog=31,nl=7,nadd(2)=[5,10]
	!integer, parameter :: nlog=35,nl=34,nadd(2)=[5,10]
	integer, parameter :: n=(nl*nlog+sum(nadd*2))*2
	integer, parameter :: set_rGB=1,set_Gt=2,set_rG=3,set_rSE=4,export_Gw=1,export_Gt=2,export_fe=3,export_entropy=4
	real(8) :: omega(n),tau(n)
	logical :: is_expinit=.true.
	type gfs(n)
		integer, len :: n
		real(8) :: rhoc(n),rhophi(n)
		complex(8) :: Gf(n),GB(n),Gc(n),Gphi(n),Gc0(n),Gphi0(n),Tmat(n),Sus(n)
		complex(8) :: SEf(n),SEB(n),SEc(n),SEphi(n)
		real(8) :: Gft(n),GBt(n),Gc0t(n),Gphi0t(n),Gct(n),Gphit(n),SEft(n),SEBt(n),Tmatt(n),Sust(n),SEct(n),SEphit(n)
		logical :: conv
	contains
		procedure :: self_consistent,analytic,set_Gfunction
		procedure :: entropy,fenergy
		procedure :: export
	end type
contains
	include "mbroyden.f90"
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
				write(ut(k),"(A)")"Tk omega Gft GBt SEft SEBt Gct Gphit SEct SEphit Tmatt Sust Gc0t Gphi0t"
				do i=1,n
					write(ut(k),"(*(e28.20))")Tk,tau(i),self%Gft(i),self%GBt(i),self%SEft(i),self%SEBt(i),self%Gct(i),self%Gphit(i),self%SEct(i),self%SEphit(i),self%Tmatt(i),self%Sust(i),self%Gc0t(i),self%Gphi0t(i)
				enddo
			end select
		end do
	end subroutine
	subroutine set_Gfunction(self,flag)
		class(gfs(*)) :: self
		integer :: flag(:),i
		do i=1,size(flag)
			select case(flag(i))
			case(set_rGB)
				! make GB a usual Green's function that satisfy kk relation
				call set_realpart(omega,self%GB)
			case(set_Gt)
				call get_ftau(omega,G=self%Gf,eta=-1,time=tau,rt=self%Gft)
				call get_ftau(omega,G=self%GB,eta=1,time=tau,rt=self%GBt)
				call get_ftau(omega,G=self%Gc0,eta=-1,time=tau,rt=self%Gc0t)
				call get_ftau(omega,G=self%Gphi0,eta=1,time=tau,rt=self%Gphi0t)
				call get_ftau(omega,G=self%SEf,eta=-1,time=tau,rt=self%SEft)
				call get_ftau(omega,G=self%SEB,eta=1,time=tau,rt=self%SEBt)
				call get_ftau(omega,G=self%Gc,eta=-1,time=tau,rt=self%Gct)
				call get_ftau(omega,G=self%Gphi,eta=1,time=tau,rt=self%Gphit)
				call get_ftau(omega,G=self%Tmat,eta=-1,time=tau,rt=self%Tmatt)
				call get_ftau(omega,G=self%Sus,eta=1,time=tau,rt=self%Sust)
			case(set_rSE)
				call set_realpart(omega,self%SEf)
				call set_realpart(omega,self%SEB)
			end select
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
		at2=atan(2d0*width)
		do i=l+1,u
			if(i<=m) then
				omega(i)=center-tan(s+(at1-s)/(m-l)*(m-i))
			else
				omega(i)=center+tan(s+(at2-s)/(u-m)*(i-m-1))
			endif
			omega(n-i+1)=-omega(i)
		enddo
	end subroutine
	real(8) function A_phi(x)
		real(8) :: x
		if(abs(x)<1.0d-19) then
			A_phi=0.0d0
		elseif(abs(x)<phicut) then
			A_phi=norm_Aphi*abs(x)**(1d0-alpha)*sign(1d0,x)
		else
			A_phi=norm_Aphi*phicut**(1d0-alpha)*sign(1d0,x)*1d0/((exp(phibeta*(x-phicut-4d0/phibeta))+1d0)*(exp(phibeta*(-x-phicut-4d0/phibeta))+1d0))
		endif
	end function
	real(8) function A_c(x)
		real(8) :: x
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
			A_c=0d0
		endif
	end function
	real(8) function integrate(Af,A,x) result(f)
		real(8), optional :: Af,A(:)
		real(8), optional :: x(:)
		integer :: i
		f=0d0
		if(present(x)) then
			if(present(A)) then
				do i=1,size(x)-1
					f=f+0.5d0*(x(i+1)-x(i))*(A(i+1)+A(i))
				enddo
			else
				do i=1,size(x)-1
					f=f+0.5d0*(x(i+1)-x(i))*(Af(x(i+1))+Af(x(i)))
				enddo
			endif
		else
			if(present(A)) then
				do i=1,size(omega)-1
					f=f+0.5d0*(omega(i+1)-omega(i))*(A(i+1)+A(i))
				enddo
			else
				do i=1,size(omega)-1
					f=f+0.5d0*(omega(i+1)-omega(i))*(Af(omega(i+1))+Af(omega(i)))
				enddo
			endif
		endif
	end function
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
	pure real(8) function ff(x)
		real(8), intent(in) :: x
		!if(x<omega(1)) then
		!ff=1d0
		!elseif(x>-omega(1)) then
		!ff=0d0
		!else
		ff=1d0/(exp(x)+1d0)
		!endif
	end function
	real(8) function fb(x)
		real(8), intent(in) :: x
		integer :: m
		real(8) :: f1,f2
		if(abs(x)<1d-20) then
			!fb=1d0/(exp(-1d-11)-1d0)+(1d0/(exp(-1d-11)-1d0)-(-exp(-1d-11)/(exp(-1d-11)-1d0)))/(-2d-11)*(x+1d-11)
			fb=-0.5d0
		elseif(x<0d0) then
			!if(x<omega(1)) then
			!fb=-1d0
			!elseif(x>-omega(1)) then
			!fb=0d0
			!else
			fb=1d0/(exp(x)-1d0)
		else
			fb=-exp(-x)/(exp(-x)-1d0)
			!endif
		endif
		!!write(*,*)x/beta,fb
		!m=find_sidx(omega*beta,x)
		!if(m<1) then
			!fb=1d0/(exp(x)-1d0)
		!elseif(m>=n) then
			!fb=1d0/(exp(x)-1d0)
		!else
			!if(omega(m)>0) then
				!f1=-exp(-beta*omega(m))/(exp(-beta*omega(m))-1d0)
				!f2=-exp(-beta*omega(m+1))/(exp(-beta*omega(m+1))-1d0)
			!elseif(omega(m+1)>0) then
				!f1=1d0/(exp(beta*omega(m))-1d0)
				!f2=-exp(-beta*omega(m+1))/(exp(-beta*omega(m+1))-1d0)
			!else
				!f1=1d0/(exp(beta*omega(m))-1d0)
				!f2=1d0/(exp(beta*omega(m+1))-1d0)
			!endif

			!fb=f1+(f2-f1)/(omega(m+1)*beta-omega(m)*beta)*(x-omega(m)*beta)
		!endif
		!!write(*,*)x/beta,fb
		!!read(*,*)
	end function
	subroutine get_ftau(omega,G,Gf,eta,time,rt)
		real(8) :: omega(:),time(:),rt(:)
		complex(8), optional :: G(:),Gf
		integer :: eta
		integer :: tau,i,k,n
		real(8) :: A1,A2,B1,B2,e1,e2,t,onepi
		procedure(real(8)), pointer :: f
		onepi=1d0/pi
		n=size(omega)
		tau=1
		if((time(size(time))+time(1))>0d0) then
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
			B1=tau*1d0/(exp((tau*beta+t)*e1)-eta*exp(e1*t))
			rt(k)=0d0
			do i=1,n-1
				e2=omega(i+1)
				if(present(G)) then
					A2=G(i+1)%im
				else
					A2=imag(Gf(e2))
				endif
				B2=tau*1d0/(exp((tau*beta+t)*e2)-eta*exp(e2*t))
				rt(k)=rt(k)+0.5d0*(e2-e1)*(A1*B1+A2*B2)
				A1=A2
				B1=B2
				e1=e2
			enddo
			rt(k)=eta*rt(k)*onepi
		enddo
	end subroutine
	!subroutine get_convolution(omega,A,B,tau,eta,rt)
		!real(8) :: omega(:)
		!complex(8) :: rt(:)
		!complex(8) :: A(:),B(:)
		!integer :: tau(2),eta(2)
		!integer :: i,j,k,iwb(size(omega)),ibw(size(omega)),m
		!real(8) :: e1,e2,omg
		!complex(8) :: A1,A2,B1,B2,rt_
		!procedure(real(8)), pointer :: f1,f2
		!select case(eta(1))
		!case(1)
			!f1 => fb
		!case(-1)
			!f1 => ff
		!end select
		!select case(eta(2))
		!case(1)
			!f2 => fb
		!case(-1)
			!f2 => ff
		!end select
		!!$OMP PARALLEL DO PRIVATE(omg,ibw,iwb,A1,B1,A2,B2,e1,e2,m,rt_)
		!!do k=1,n/2
			!do k=1,n
			!omg=omega(k)
			!do i=1,n
				!ibw(i)=find_sidx(omega,tau(1)*(omg-tau(2)*omega(i)))
			!end do
			!A1=A(1)
			!do i=1,n
				!iwb(i)=find_sidx(omega,tau(2)*(omg-tau(1)*omega(i)))
			!end do
			!i=iwb(1)
			!!B1=merge(B(i)+(B(i+1)-B(i))/(omega(i+1)-omega(i))*(tau(2)*(omg-tau(1)*omega(1))-omega(i)),cmplx(0d0,kind=8),i>0.and.i<n)
			!B1=merge(B(i)+(B(i+1)-B(i))/(omega(i+1)-omega(i))*(tau(2)*(omg-tau(1)*omega(1))-omega(i)),cmplx(1d0/(tau(2)*(omg-tau(1)*omega(1))),0d0,kind=8),i>0.and.i<n)
			!e1=omega(1)
			!rt_=0d0
			!do i=1,n-1
				!do j=iwb(i),iwb(i+1),sign(1,iwb(i+1)-iwb(i))
					!if(j==iwb(i+1)) then
						!e2=omega(i+1)
						!A2=A(i+1)
						!!B2=merge(B(j)+(B(j+1)-B(j))/(omega(j+1)-omega(j))*(tau(2)*(omg-tau(1)*e2)-omega(j)),cmplx(0d0,kind=8),j>0.and.j<n)
						!B2=merge(B(j)+(B(j+1)-B(j))/(omega(j+1)-omega(j))*(tau(2)*(omg-tau(1)*e2)-omega(j)),cmplx(1d0/(tau(2)*(omg-tau(1)*e2)),0d0,kind=8),j>0.and.j<n)
					!else
						!m=j+(1-tau(1)*tau(2))/2
						!e2=tau(1)*(omg-tau(2)*omega(m))
						!B2=B(m)
						!m=ibw(m)
						!A2=A(m)+(A(m+1)-A(m))/(omega(m+1)-omega(m))*(e2-omega(m))
					!endif
					!rt_=rt_+0.5d0*(e2-e1)*cmplx(&
						!(-tau(1)*f1(beta*tau(1)*e1)*imag(a1)*real(b1)+tau(2)*eta(1)*eta(2)*f2(beta*(tau(1)*e1-omg))*real(a1)*imag(b1))+&
						!(-tau(1)*f1(beta*tau(1)*e2)*imag(a2)*real(b2)+tau(2)*eta(1)*eta(2)*f2(beta*(tau(1)*e2-omg))*real(a2)*imag(b2))&
						!,&
						!(-tau(1)*tau(2)*f1(beta*tau(1)*e1)*imag(A1)*imag(B1)+tau(1)*tau(2)*eta(1)*eta(2)*f2(beta*(tau(1)*e1-omg))*imag(A1)*imag(B1))+&
						!(-tau(1)*tau(2)*f1(beta*tau(1)*e2)*imag(A2)*imag(B2)+tau(1)*tau(2)*eta(1)*eta(2)*f2(beta*(tau(1)*e2-omg))*imag(A2)*imag(B2)),kind=8)
					!!if(k==1.or.k==n) then
						!!write(233,"(*(es12.4))")&
						!!f1(beta*tau(1)*e1),imag(A1)*real(B1),f1(beta*tau(1)*e1)*imag(A1)*real(B1)
					!!endif
					!A1=A2
					!B1=B2
					!e1=e2
				!enddo
			!enddo
			!rt(k)=rt(k)-eta(1)*rt_/pi
			!!rt(n+1-k)=rt(n+1-k)+cmplx(-eta(2)*rt_%re,eta(2)*rt_%im,kind=8)/pi
		!enddo
		!!$OMP END PARALLEL DO
		!!stop
	!end subroutine
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
	subroutine  self_consistent(self,tol,niter,rate)
		class(gfs(*)) :: self
		real(8) :: tol,rate
		integer :: niter(:)
		real(8) :: x(n*2),x_(n*2),nf(1),tmp(1)
		logical :: flag
		integer :: i1,i2,k,i,j,ei
		self%conv=.false.
		self%Gf=0d0
		self%GB=0d0
		self%Gc=0d0
		self%Gphi=0d0
		call mbroyden(0,x,x_,rate,nsave=40,id=2)
		do i2=1,niter(2)
			flag=.false.
			ei=0
			call set_realpart(omega,self%Gc0)
			call set_realpart(omega,self%Gphi0)
			!call get_ftau(omega,G=self%Gc0,eta=-1,time=[1d-70],rt=tmp)
			do i1=1,niter(1)
				call set_realpart(omega,self%SEf)
				call set_realpart(omega,self%SEB)
				self%Gf=1d0/(omega-self%SEf)
				self%GB=1d0/(-1d0/Jk-self%SEB)
				call get_ftau(omega,G=self%Gf,eta=-1,time=[-1d-70],rt=nf)

				x=[self%SEf%im,self%SEB%im]

				self%SEB=0d0
				self%SEf=0d0
				call get_convolution(omega,A=-kp*self%Gc0,B=self%GB+Jk,tau=[1,1],eta=[-1,1],rt=self%SEf)
				call get_convolution(omega,A=self%Gc0,B=self%Gf,tau=[1,-1],eta=[-1,-1],rt=self%SEB)
				call get_convolution(omega,A=-g**2*self%Gphi0,B=self%Gf,tau=[1,1],eta=[1,-1],rt=self%SEf)


				x_=[self%SEf%im,self%SEB%im]
				if(all(abs(x-x_)/(abs(x_)+1d-70)<tol)) then
					!if(ei==0) ei=1
					!if(ei>20) then
						write(*,*)"converged"
						if(abs(0.5d0-nf(1))>1d-2) then
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
					write(*,*)"not converged",maxval(abs(x-x_)/(abs(x_)+1d-70)),"***************"
					is_expinit=.false.
					return
				endif

				!if(flag) then
				if(ei>0) then
					ei=ei+1
					!write(*,*)ei,maxval(abs(x-x_)/(abs(x_)+1d-70))
					!if(abs(Tk-1.d-9)<1d-4) then
					!x=(1d0-rate)*x+rate*x_
					!x=(x+rate*x_)/(1d0+rate)
					!endif

					!ei=30
					!x=x_
				else
					call mbroyden(i1,x,x_,id=1)
				endif
				!x=x_


				self%SEf%im=x(1:n)
				self%SEB%im=x(n+1:)
				!write(*,*)"free energy:",free_energy(rhoc,rhophi,Af,AB,iSEf,iSEB,1)
				!read(*,*)
			enddo
			write(*,*)i1,maxval(abs(x-x_)/(abs(x_)+1d-70)),"nf:",nf*2d0

			x=[self%Gc0%im,self%Gphi0%im]

			if(niter(2)>1) then
				!call get_convolution(omega,A=-self%Gf,B=self%GB,tau=[1,-1],eta=[-1,1],rt=self%SEc)

				self%Gc=self%Gc0+self%Gc0*self%SEc*self%Gc0

				do i=1,n
					self%Gc0=1d0/(omega(i)-omega-self%SEc(i))
					self%Gc(i)=cmplx(integrate(A=self%rhoc*self%Gc0%re),integrate(A=self%rhoc*self%Gc0%im),kind=8)
				enddo

				self%Gc0=1d0/(1d0/self%Gc+self%SEc)
			endif
			x_=[self%Gc0%im,self%Gphi0%im]

			if(all(abs(x-x_)/(abs(x_)+1d-70)<tol)) then
				write(*,*)"converged"
				exit
			else
				write(*,*)"|___",i2,maxval(abs(x-x_)/(abs(x_)+1d-70))
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
		write(*,*)i2,maxval(abs(x-x_)/(abs(x_)+1d-70))
		call mbroyden(huge(1),x,x_)
	end subroutine
	real(8) function entropy(self) result(rt)
		class(gfs(*)) :: self
		real(8) :: df(n),db(n)
		df=omega*(0.5d0*beta/cosh(0.5d0*omega*beta))**2
		db=omega*(0.5d0*beta/sinh(0.5d0*omega*beta))**2

		!rt=-((integrate(A=&
			!ff_*(imag(log(-1d0/Gf)+Gf*(omega-1d0/Gf))+0.5d0*(1d0+sign(1d0,omega))*pi)&
			!+kp*fb_*imag(log(-1d0/GB)+GB*(-1d0/Jk-1d0/GB))&
			!!+kp*fb_*imag(log(-1d0/GB)+GB*SEB)&
			!!ff_*(imag(log(-1d0/Gf))+0.5d0*(1d0+sign(1d0,omega))*pi)&
			!!+kp*fb_*imag(log(-1d0/GB))&
			!)/pi&
				!-1d0/beta*log(2d0)&
				!+integrate(A=-kp*(-Gft)*Gct*GBt,x=tau)&
			!))
		rt=-(integrate(A=&
			df*(imag(log(-1d0/self%Gf))+self%Gf%re*imag(omega-1d0/self%Gf)+0.5d0*(1d0+sign(1d0,omega))*pi)&
			!df*(imag(log(-1d0/self%Gf))+self%Gf%re*self%SEf%im)&
			!+kp*db*(imag(log(-1d0/self%GB))+self%GB%re*self%SEB%im)&
			+kp*db*(imag(log(-1d0/self%GB))+(self%GB%re+Jk)*imag(-1d0/Jk-1d0/self%GB))&
			)/pi&
			-log(2d0)&
			!+integrate(A=-kp*(-self%Gft)*self%Gct*self%GBt,x=tau)&
			)
			
	end function
	real(8) function fenergy(self) result(rt)
		class(gfs(*)) :: self
		real(8) :: ff_(n),fb_(n)
		integer :: i
		do i=1,n
			ff_(i)=ff(beta*omega(i))
			fb_(i)=fb(beta*omega(i))
		enddo

		!write(*,*)integrate(A=(-self%Gft)*self%SEft,x=tau),"Gft*SEft"
		!write(*,*)integrate(A=-kp*self%GBt*self%SEBt,x=tau),"GBt*SEBt"
		!write(*,*)integrate(A=-kp*(-self%Gft)*self%Gc0t*self%GBt,x=tau),"GftGc0tGBt"
		!write(*,*)integrate(A=omega*ff_*(-self%Gf%im/pi)),"rhoGf0^-1"
		!write(*,*)integrate(A=-ff_*imag(self%Gf*self%SEf)/pi),"Gf*SEf"
		!write(*,*)integrate(A=fb_*imag((self%GB+Jk)*self%SEB)/pi*(-kp)),"(GB+J)*SEB"

		rt=(integrate(A=&
			ff_*(imag(log(-1d0/self%Gf)+self%Gf*(omega-1d0/self%Gf))+0.5d0*(1d0+sign(1d0,omega))*pi)&
			!+kp*fb_*imag(log(-1d0/self%GB)+self%GB*(-1d0/Jk-1d0/self%GB))&
			!+kp*fb_*imag(log(-1d0/self%GB)+self%GB*self%SEB)&
			+kp*fb_*imag(log(-1d0/self%GB)+(self%GB+Jk)*(-1d0/Jk-1d0/self%GB))&
			!ff_*(imag(log(-1d0/self%Gf))+0.5d0*(1d0+sign(1d0,omega))*pi)&
			!+kp*fb_*imag(log(-1d0/self%GB))&
			)/pi&
			-1d0/beta*log(2d0)&
			+integrate(A=-kp*(-self%Gft)*self%Gct*self%GBt,x=tau)&
			)
	end function
	subroutine analytic(self,a_f,a_b)
		class(gfs(*)) :: self
		real(8) :: a_f,a_b
		real(8) :: A1,B1,A0,t,sh,ch,lsh,lch,lgam_a,lgam_b
		integer :: i,j
		A1=1d0
		B1=1d0
		A0=1d0
		t=(gamma(a_b)*gamma(1d0+r+a_f)/(pi**2*(1d0+1d0/tan(pi*(r+a_f)/2d0)**2)*gamma(1d0+r))/(A1*B1*A0))**(1d0/(a_f+a_b))
		do i=1,n
			if(exp(beta*abs(omega(i)))>1d100) then
				lsh=beta*abs(omega(i))/2d0-log(2d0)
				lch=beta*abs(omega(i))/2d0-log(2d0)
				lgam_a=log(2d0*pi)-beta*abs(omega(i))/2d0+(a_f-1d0)*(log(beta*abs(omega(i)))-log(2d0*pi))
				lgam_b=log(2d0*pi)-beta*abs(omega(i))/2d0+(a_b-1d0)*(log(beta*abs(omega(i)))-log(2d0*pi))
			else
				lsh=log(sinh(beta*abs(omega(i))/2d0))
				lch=log(cosh(beta*abs(omega(i))/2d0))
				lgam_a=real(lgamma(a_f/2+img*beta*omega(i)/(2d0*pi))+lgamma(a_f/2-img*beta*omega(i)/(2d0*pi)))
				lgam_b=real(lgamma(a_b/2+img*beta*omega(i)/(2d0*pi))+lgamma(a_b/2-img*beta*omega(i)/(2d0*pi)))
			endif
			sh=real(exp((a_f-1d0)*(log(2d0*pi)-log(beta))+lgam_a+lsh))*sign(1d0,omega(i))
			ch=real(exp((a_f-1d0)*(log(2d0*pi)-log(beta))+lgam_a+lch))
			self%Gf(i)=A1*t**a_f/gamma(a_f)*cmplx(1d0/tan(pi*a_f/2d0)*sh,-ch,kind=8)
			sh=real(exp((a_b-1d0)*(log(2d0*pi)-log(beta))+lgam_b+lsh))*sign(1d0,omega(i))
			ch=real(exp((a_b-1d0)*(log(2d0*pi)-log(beta))+lgam_b+lch))
			self%GB(i)=B1*t**a_b/gamma(a_b)*cmplx(-tan(pi*a_b/2d0)*ch,-sh,kind=8)
			self%Gc(i)%im=A0*(abs(omega(i))/fcut)**r
		enddo

		self%SEf=omega-1d0/self%Gf
		self%SEB=-1d0/Jk-1d0/self%GB
	end subroutine
end module
program main
	use global
	implicit none
	type(gfs(n)) :: phy,aphy
	real(8) :: omg(2),dat(2,2)
	integer :: i,j,e

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
	call set_omega(omega(1:nl*nlog*2),base,nl,bcut)
	call add_omega(omega(1:nl*nlog*2+sum(nadd(1:1))*4),1d0,0.1d0,nadd(1))
	call add_omega(omega(1:nl*nlog*2+sum(nadd(1:2))*4),0.8d0,0.16d0,nadd(2))
	write(*,"('n: ',i4,', max freq: ',es10.2,', min freq: ',es10.2)")n,maxval(abs(omega)),minval(abs(omega))

	do 
		is_expinit=.true.
		norm_Ac=1d0
		do i=1,n
			phy%rhoc(i)=A_c(omega(i))
			phy%rhophi(i)=A_phi(omega(i))
		enddo
		norm_Ac=1d0/integrate(A=phy%rhoc)
		phy%rhoc=phy%rhoc*norm_Ac
		phy%Gc0%im=-phy%rhoc*pi
		phy%Gphi0%im=-phy%rhophi*pi

		i=1
		rewind(10)
		read(10,"(*(e28.20))")omg(1),dat(1,:)
		o:	do 
			read(10,"(*(e28.20))",iostat=e)omg(2),dat(2,:)
			if(e==-1) then
				backspace(10)
				backspace(10)
				backspace(10)
				read(10,"(*(e28.20))")omg(1),dat(1,:)
				dat(2,:)=dat(2,:)+(dat(2,:)-dat(1,:))/(omg(2)-omg(1))*(-omega(1)-omg(2))
				omg(2)=-omega(1)
				read(10,"(*(e28.20))")omg(1),dat(1,:)
			endif
			do 
				if(omg(2)-omega(i)>=0d0) then
					dat(1,:)=dat(1,:)+(dat(2,:)-dat(1,:))/(omg(2)-omg(1))*(omega(i)-omg(1))
					phy%SEf(i)%im=dat(1,1)
					phy%SEB(i)%im=dat(1,2)
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
		rewind(13)
		write(25,"(A)")"T Sn Sa"
		write(26,"(A)")"T Fn Fa pT pFn pFa"
		do
			read(13,*,iostat=e)Tk
			beta=1d0/Tk
			if(e==-1) exit
			write(*,"(*(A4,es9.2))")"Tk:  ",Tk," Jk:",Jk," kp:",kp," r: ",r," g: ",g

			tau(1:n/2)=omega(n/2+1:n)/(3d0*omega(n)-omega(n-1))/Tk
			tau(n/2+1:)=1d0/Tk-tau(n/2:1:-1)

			call phy%self_consistent(1d-5,[6000,1],2d-2)
			call phy%set_Gfunction([set_Gt,set_rSE])
			call aphy%analytic(a_f=&
				0.41d0&
				,a_b=&
				0.29d0&
			)
			call aphy%set_Gfunction([set_Gt])
			write(25,"(*(e28.20))")Tk,phy%entropy(),aphy%entropy()
			write(26,"(*(e28.20))")Tk,phy%fenergy(),aphy%fenergy()
			if(.not.phy%conv) then
				exit
			endif

			call phy%export([22],[export_Gw])
			!write(27,"(*(e28.20))")Tk,-phy%Gf(n/2)%im

			if(is_expinit) then
				write(*,*)"export initial...."
				rewind(10)
				do i=1,n
					write(10,"(*(e28.20))")omega(i),phy%SEf(i)%im,phy%SEB(i)%im
				enddo
				is_expinit=.false.
			endif
		enddo
		stop
		Jk=Jk+0.1d0
		if(Jk>2.2d0-1d-4) then
			exit
		endif
	enddo
end program
