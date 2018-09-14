include "fn.f90"
module global
	use M_fn
	use M_const
	use M_matrix
	use mkl_service
	implicit none
	!real(8) :: alpha=0.1d0,r=0.01d0,kp=0.5d0,fcut=1d0,fbeta=1d0,phicut=5d-2,phibeta=200d0,Jk=2.0d0,g=0d0
	!real(8) :: alpha=0.8d0,r=0.3d0,kp=0.5d0,fcut=1d0,fbeta=1d0,phicut=5d-2,phibeta=200d0,Jk=1d0,g=0d0
	real(8) :: alpha=0.8d0,r=0.3d0,kp=0.3d0,fcut=1d0,fbeta=1d0,phicut=5d-2,phibeta=200d0,Jk=1.8d0,g=0d0
	!real(8) :: alpha=0.3d0,r=0.3d0,kp=0.3d0,fcut=1d0,fbeta=20d0,phicut=5d-2,phibeta=200d0,Jk=1.5d0,g=0.2d0
	!real(8) :: alpha=0.3d0,r=0d0,kp=0.5d0,fcut=1d0,fbeta=20d0,phicut=5d-2,phibeta=200d0,Jk=1d0,g=0.5d0
	!real(8) :: alpha=0.3d0,r=0d0,kp=0.3d0,fcut=1d0,fbeta=1d0,phicut=5d-2,phibeta=200d0,Jk=0.7d0,g=0d0
	real(8) :: Tk=4d-3,beta,norm_Ac=1d0,norm_Aphi=1d0,base=2.3d0
	integer, parameter :: nlog=31,nl=7,nadd(2)=[5,10]
	!integer, parameter :: nlog=35,nl=50,nadd(2)=[5,10]
	integer, parameter :: n=(nl*nlog+sum(nadd*2))*2
	!integer, parameter :: n=536
	!real(8) :: omega((n+nadd1*2+nadd2*2)*2)
	real(8) :: omega(n),tau(n)
	logical :: is_expinit=.true.
	real(8) :: rG0(n),iG0(n),pTk=0d0
	real(8) :: &
		a_f=0.41d0,&
		a_b=0.29d0
		!a_f=0.03d0,&
		!a_b=0.67d0
		!a_f=0.7692307692307689d0,&
		!a_b=0.230769d0
		!a_f=0.666666666d0,&
		!a_b=0.333333d0
		!a_f=0.543604996816506d0,&
		!a_b=0.256395d0
		!a_f=0.014569749254921006d0,&
		!a_b=0.78543d0
		!a_f=0.229d0,&
		!a_b=0.47d0
		!a_f=0.106d0,&
		!a_b=0.593d0
		!a_f=0.4d0,&
		!a_b=1.7d0
contains
	include "mbroyden.f90"
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
	SUBROUTINE KRAKRO(x,im,re)
		real(8) :: x(:),im(:),re(:)
		real(8) :: lg(size(x)),dx(size(x)),st(size(x))
		integer :: i,n,j
		n=size(x)
		do i=1,n
			dx(i)=x(min(n,i+1))-x(max(1,i-1))
			st(i)=(im(min(n,i+1))-im(max(1,i-1)))/dx(i)
			if ((i.gt.1).and.(i.lt.n)) then
				lg(i)=dlog(dabs((x(i)-x(1))/(x(i)-x(n))))*im(i)
			endif
		enddo
		do j=1,n
			re(j)=0d0
			if ((j.gt.1).and.(j.lt.n)) re(j)=-2.0d0*lg(j)
			do i=1,n
				if (i.ne.j) then
					re(j)=re(j)+(im(i)-im(j))/(x(i)-x(j))*dx(i)
				else
					re(j)=re(j)+st(i)*dx(i)
				endif
			enddo
		enddo
		re=0.5d0/pi*re
	end
	subroutine get_realpart(omega,im,re)
		real(8) :: omega(:),im(:),re(:)
		real(8) ::dx
		integer :: n,i,j,dn(2)
		n=size(omega)
		re(1)=0d0
		re(n)=0d0
		re(2:n-1)=2d0*im(2:n-1)*log((omega(n)-omega(2:n-1))/(omega(2:n-1)-omega(1)))
		!re(2:n-1)=0d0
		do j=1,n
			dn=[min(n,j+1),max(1,j-1)]
			dx=omega(dn(1))-omega(dn(2))
			re(j)=re(j)+im(dn(1))-im(dn(2))
			do i=1,n
				if(j/=i) then
					re(i)=re(i)+(im(j)-im(i))/(omega(j)-omega(i))*dx
				endif
			enddo
		enddo
		re=0.5d0/pi*re
	end subroutine
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
	subroutine get_ftau(omega,A,Af,eta,time,rt)
		real(8) :: omega(:),time(:),rt(:)
		real(8), optional :: A(:),Af
		integer :: eta
		integer :: tau,i,k,n
		real(8) :: A1,A2,B1,B2,e1,e2,t
		procedure(real(8)), pointer :: f
		n=size(omega)
		tau=1
		if((time(size(time))+time(1))>0d0) then
			tau=-1
		endif
		do k=1,size(time)
			t=time(k)
			e1=omega(1)
			if(present(A)) then
				A1=A(1)
			else
				A1=Af(e1)
			endif
			B1=tau*1d0/(exp((tau*beta+t)*e1)-eta*exp(e1*t))
			rt(k)=0d0
			do i=1,n-1
				e2=omega(i+1)
				if(present(A)) then
					A2=A(i+1)
				else
					A2=Af(e2)
				endif
				B2=tau*1d0/(exp((tau*beta+t)*e2)-eta*exp(e2*t))
				rt(k)=rt(k)+0.5d0*(e2-e1)*(A1*B1+A2*B2)
				A1=A2
				B1=B2
				e1=e2
			enddo
			rt(k)=-eta*rt(k)
		enddo
	end subroutine
	subroutine get_convolution_complex(omega,A,B,tau,eta,rt)
		real(8) :: omega(:)
		complex(8) :: rt(:)
		complex(8) :: A(:),B(:)
		integer :: tau(2),eta(2)
		integer :: i,j,k,n,iwb(size(omega)),ibw(size(omega)),m
		real(8) :: e1,e2,omg
		complex(8) :: A1,A2,B1,B2
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
		!$OMP PARALLEL DO PRIVATE(omg,ibw,iwb,A1,B1,A2,B2,e1,e2,m)
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
			rt(k)=0d0
			do i=1,n-1
				do j=iwb(i),iwb(i+1),sign(1,iwb(i+1)-iwb(i))
					if(j==iwb(i+1)) then
						e2=omega(i+1)
						A2=A(i+1)
						B2=merge(B(j)+(B(j+1)-B(j))/(omega(j+1)-omega(j))*(tau(2)*(omg-tau(1)*e2)-omega(j)),cmplx(0d0,kind=8),j>0.and.j<n)
					else
						m=j+(1-tau(1)*tau(2))/2
						e2=tau(1)*(omg-tau(2)*omega(m))
						B2=merge(B(m),cmplx(0d0,kind=8),m>0.and.m<n)
						m=ibw(m)
						A2=merge(A(m)+(A(m+1)-A(m))/(omega(m+1)-omega(m))*(e2-omega(m)),cmplx(0d0,kind=8),m>0.and.m<n)
					endif
					rt(k)=rt(k)+0.5d0*(e2-e1)*cmplx(&
						(-tau(1)*f1(beta*tau(1)*e1)*imag(A1)*real(B1)+tau(2)*eta(1)*eta(2)*f2(beta*(tau(1)*e1-omg))*real(A1)*imag(B1))+&
						(-tau(1)*f1(beta*tau(1)*e2)*imag(A2)*real(B2)+tau(2)*eta(1)*eta(2)*f2(beta*(tau(1)*e2-omg))*real(A2)*imag(B2))&
						,&
						(-tau(1)*tau(2)*f1(beta*tau(1)*e1)*imag(A1)*imag(B1)+tau(1)*tau(2)*eta(1)*eta(2)*f2(beta*(tau(1)*e1-omg))*imag(A1)*imag(B1))+&
						(-tau(1)*tau(2)*f1(beta*tau(1)*e2)*imag(A2)*imag(B2)+tau(1)*tau(2)*eta(1)*eta(2)*f2(beta*(tau(1)*e2-omg))*imag(A2)*imag(B2)))
					A1=A2
					B1=B2
					e1=e2
				enddo
			enddo
			rt(k)=-eta(1)*rt(k)/pi
			!rt(n+1-k)=-eta(1)*eta(2)*rt(k)
		enddo
		!$OMP END PARALLEL DO
	end subroutine
	subroutine get_convolution(omega,A,B,Af,Bf,tau,eta,rt)
		real(8) :: omega(:),rt(:)
		real(8), optional :: A(:),B(:),Af,Bf
		integer :: tau(2),eta(2)
		integer :: i,j,k,n,iwb(size(omega)),ibw(size(omega)),m
		real(8) :: e1,e2,A1,A2,B1,B2,omg
		procedure(real(8)), pointer :: f1,f2
		if(.not.(present(A).xor.present(Af))) stop "error"
		if(.not.(present(B).xor.present(Bf))) stop "error"
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
		!!$OMP PARALLEL DO PRIVATE(omg,ibw,iwb,A1,B1,A2,B2,e1,e2,m)
		do k=1,n/2
			!do k=1,n
			omg=omega(k)
			if(present(A)) then
				do i=1,n
					ibw(i)=find_sidx(omega,tau(1)*(omg-tau(2)*omega(i)))
				end do
				A1=A(1)
			else
				A1=Af(omega(1))
			endif
			do i=1,n
				iwb(i)=find_sidx(omega,tau(2)*(omg-tau(1)*omega(i)))
			end do
			if(present(B)) then
				i=iwb(1)
				B1=merge(B(i)+(B(i+1)-B(i))/(omega(i+1)-omega(i))*(tau(2)*(omg-tau(1)*omega(1))-omega(i)),0d0,i>0.and.i<n)
			else
				B1=Bf(tau(2)*(omg-tau(1)*omega(1)))
			endif
			e1=omega(1)
			rt(k)=0d0
			do i=1,n-1
				do j=iwb(i),iwb(i+1),sign(1,iwb(i+1)-iwb(i))
					if(j==iwb(i+1)) then
						e2=omega(i+1)
						if(present(A)) then
							A2=A(i+1)
						else
							A2=Af(omega(i+1))
						endif
						if(present(B)) then
							B2=merge(B(j)+(B(j+1)-B(j))/(omega(j+1)-omega(j))*(tau(2)*(omg-tau(1)*e2)-omega(j)),0d0,j>0.and.j<n)
						else
							B2=Bf(tau(2)*(omg-tau(1)*e2))
						endif
					else
						m=j+(1-tau(1)*tau(2))/2
						e2=tau(1)*(omg-tau(2)*omega(m))
						if(present(B)) then
							B2=merge(B(m),0d0,m>0.and.m<n)
						else
							B2=Bf(omega(m))
						endif
						if(present(A)) then
							m=ibw(m)
							A2=merge(A(m)+(A(m+1)-A(m))/(omega(m+1)-omega(m))*(e2-omega(m)),0d0,m>0.and.m<n)
						else
							A2=Af(e2)
						endif
					endif
					rt(k)=rt(k)+0.5d0*(e2-e1)*&
						((-f1(beta*tau(1)*e1)+eta(1)*eta(2)*f2(beta*(tau(1)*e1-omg)))*A1*B1+&
						(-f1(beta*tau(1)*e2)+eta(1)*eta(2)*f2(beta*(tau(1)*e2-omg)))*A2*B2)
					A1=A2
					B1=B2
					e1=e2
				enddo
			enddo
			rt(k)=-pi*eta(1)*tau(1)*tau(2)*rt(k)
			rt(n+1-k)=-eta(1)*eta(2)*rt(k)
		enddo
		!!$OMP END PARALLEL DO
	end subroutine
	logical function  self_consistent(tol,niter,rate,iSEf,iSEB,rhoc,rhophi,Ac0,Aphi0,Af,AB,Ac,Aphi) result(conv)
		real(8) :: tol,rate,Af(:),AB(:),iSEF(:),iSEB(:),rhoc(:),rhophi(:),Ac0(:),Aphi0(:),Ac(:),Aphi(:)
		real(8) :: A0(n),Al(n),tmp(n)
		complex(8) :: Gl(n),G0(n),SE(n),Tmat(n)
		integer :: niter(:)
		real(8) :: iSEf1(n),iSEf2(n),iSEf_(n),rSEf(n),iSEB_(n),rSEB(n),Af_(n),AB_(n)
		real(8) :: x(n*2),x_(n*2),nf(1),tmp1,sx(n*2),sx_(n*2)
		logical :: flag
		integer :: i1,i2,k,i,j,ei
		complex(8) :: Gf(n),GB(n)
		conv=.true.
		beta=1d0/Tk
		Af=0d0
		AB=0d0
		Ac=0d0
		Aphi=0d0
		call mbroyden(0,x,x_,rate,nsave=40,id=2)
		do i2=1,niter(2)
			flag=.false.
			ei=0
			do i1=1,niter(1)
				call get_realpart(omega,iSEf,rSEf)
				call get_realpart(omega,iSEB,rSEB)
				!call krakro(omega,iSEf,rSEf)
				!call krakro(omega,iSEB,rSEB)
				Af_=-1d0/pi*iSEf/((omega-rSEf)**2+iSEf**2)
				AB_=-1d0/pi*iSEB/((-1d0/Jk-rSEB)**2+iSEB**2)
				call get_ftau(omega,A=Af_,eta=-1,time=[-1d-70],rt=nf)

				!x=[Af*pi,AB]
				!x_=[Af_*pi,AB_]
				!write(*,*)maxval(abs(x-x_)/(abs(x_)+1d-70))
				!if(all(abs(x(n/2:find_sidx(omega,1d0))-x_(n/2:find_sidx(omega,1d0)))/(x_(n/2:find_sidx(omega,1d0))+1d-70)<tol)) then
					!if(ei==0) ei=1
					!if(ei>20) then
						!write(*,*)"converged"
						!if(abs(0.5d0-nf(1))>1d-2) then
							!is_expinit=.false.
						!endif
						!exit
					!endif
				!endif
				!if(i1==niter(1)) then
					!write(*,*)"not converged",maxval(abs(x-x_)/(abs(x_)+1d-70)),"***************"
					!is_expinit=.false.
					!stop
				!endif

				Af=Af_
				AB=AB_


				call get_convolution(omega,A=Ac0,B=AB,tau=[1,1],eta=[-1,1],rt=iSEf1)
				call get_convolution(omega,A=Ac0,B=Af,tau=[1,-1],eta=[-1,-1],rt=iSEB_)
				call get_convolution(omega,A=Aphi0,B=Af,tau=[1,1],eta=[1,-1],rt=iSEf2)

				!call get_convolution(omega,Af=A_c,B=AB,tau=[1,1],eta=[-1,1],rt=iSEf1)
				!call get_convolution(omega,Af=A_c,B=Af,tau=[1,1],eta=[-1,-1],rt=iSEB_)
				!call get_convolution(omega,Af=A_phi,B=Af,tau=[1,1],eta=[1,-1],rt=iSEf2)


				iSEf_=-kp*iSEf1-g**2*iSEf2
				iSEB_=iSEB_


				x=[iSEf,iSEB]
				x_=[iSEf_,iSEB_]
				!if(ei>20) then
					!write(*,*)"---->",sum(abs(x-x_))
					!return
				!endif
				!write(*,*)x(n/4::n/2),x_(n/4::n/2)
				!read(*,*)
				if(all(abs(x-x_)/(abs(x_)+1d-70)<tol)) then
				!if(all(abs(x-x_)<tol)) then
					!if(ei==0) ei=1
					!if(ei>20) then
						write(*,*)"converged"
						if(abs(0.5d0-nf(1))>1d-2) then
							is_expinit=.false.
							write(*,*)"nf is not correct****************************",nf(1)
							conv=.false.
							return
						endif
						exit
					!endif
				endif
				if(i1==niter(1)) then
					write(*,*)"not converged",maxval(abs(x-x_)/(abs(x_)+1d-70)),"***************"
					is_expinit=.false.
					conv=.false.
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

				iSEf=x(1:n)
				iSEB=x(n+1:)
				!write(*,*)"free energy:",free_energy(rhoc,rhophi,Af,AB,iSEf,iSEB,1)
				!read(*,*)
			enddo
			write(*,*)i1,maxval(abs(x-x_)/(abs(x_)+1d-70)),"nf:",nf*2d0,rSEB(n/2),1d0/Jk

			x=[Ac0,Aphi0]

			!call get_convolution(omega,A=Af,B=AB,tau=[1,-1],eta=[-1,1],rt=Al)
			!tau(1:n/2)=asin(pi*Tk/(omega(1:n/2)/omega(1)*(1000000d0-pi*Tk)+pi*Tk))/(pi*Tk)
			!tau(n/2+1:)=1d0/Tk-tau(1:n/2)
			!call get_ftau(omega,A=Al,eta=-1,time=tau,rt=tmp)

			!write(23,"(A)")"Tk tau Tmat"
			!do j=1,n/2
			!write(23,"(*(e28.20))")Tk,tau(j),tmp(j)
			!enddo
			!write(23,"(x/)")

			!write(22,"(A)")"omega G1 G2"
			!do j=1,n
				!write(22,"(*(e28.20))")omega(j),tmp(j),Al(j)
			!enddo
			!write(22,"(x/)")

			if(niter(2)>1) then
				!write(*,*)Af(1:n:n-1),AB(1:n:n-1)
				call get_convolution(omega,A=-Af,B=AB,tau=[1,-1],eta=[-1,1],rt=Al)

				!call get_realpart(omega,-pi*Af,tmp)
				!G0=cmplx(tmp,-pi*Af)
				!call get_realpart(omega,-pi*AB,tmp)
				!Gl=cmplx(tmp,-pi*AB)
				!call get_convolution_complex(omega,A=G0,B=Gl,tau=[1,-1],eta=[-1,1],rt=SE)
				!call get_realpart(omega,Al,tmp)
				!write(22,"(A)")"omega kk iG1 conv iG2 this"
				!do j=1,n
					!write(22,"(*(e28.20))")omega(j),tmp(j),Al(j),SE(j)
				!enddo
				!write(22,"(x/)")
				!stop


				!call get_convolution(omega,A=-Af,B=AB,tau=[1,-1],eta=[-1,1],rt=Al)
				call get_realpart(omega,Al,tmp)
				Tmat=cmplx(tmp,Al)


				call get_realpart(omega,-Ac0*pi,tmp)
				G0=cmplx(tmp,-Ac0*pi)
				Gl=G0+G0*Tmat*G0
				Ac=-imag(Gl)/pi



				!Ac=merge(1d-10,Ac,Ac<=0d0)
				!Gl=cmplx(tmp,-Ac*pi)

				!SE=1d0/G0-1d0/Gl
				!write(22,"(A)")"omega rG1 rG2"
				!do j=1,n
					!write(22,"(*(e28.20))")omega(j),SE(j),Gl(j)
				!enddo
				!write(22,"(x/)")
				!stop


				do i=1,n
					G0=1d0/(omega(i)-omega-SE(i))
					Gl(i)=cmplx(integrate(A=rhoc*real(G0)),integrate(A=rhoc*imag(G0)))
				enddo
				write(*,*)sum(abs((Gl+img*sqrt(pi)*wzag(omega-SE,1d-10))))
				Gl=-img*sqrt(pi)*wzag(omega-SE,1d-10)
				!stop

				Gl=cmplx(tmp,imag(Gl))

				G0=1d0/(1d0/Gl+SE)
				Ac0=-imag(G0)/pi
				!Ac0=merge(1d-10,Ac0,Ac0<0d0)
				!Ac0=Abs(Ac0)
				!write(*,*)integrate(A=Ac0)
				Ac0=Ac0/abs(integrate(A=Ac0))
				!write(*,*)integrate(A=Ac0)
				!stop
				!stop
			endif
			x_=[Ac0,Aphi0]
			!write(22,"(A)")"Tk omega Ac0 Ac iTmat iSE Af"
			!do i=1,n
			!write(22,"(*(e28.20))")Tk,omega(i),Ac0(i),Ac(i),-imag(Tmat(i)),-imag(SE(i)),Af(i)
			!enddo
			!write(22,"(x/)")

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
			Ac0=x(1:n)
			Aphi0=x(n+1:)
		enddo
		write(*,*)i2,maxval(abs(x-x_)/(abs(x_)+1d-70))
		call mbroyden(huge(1),x,x_)
	end function
	!real(8) function entropy_coleman(Ac,Af,AB) result(rt)
	real(8) function entropy_coleman(Ac,iSEf,iSEB) result(rt)
		!real(8) :: Ac(:),Af(:),AB(:)
		real(8) :: Ac(:),iSEf(:),iSEB(:)
		complex(8) :: Gf(n),GB(n),SEf(n),SEB(n)
		real(8) :: tmp(n),df(n),db(n),SEc(n),iSEf_(n),iSEB_(n),Gft(n),GBt(n),ff_(n),fb_(n)
		real(8) :: A1,B1,A0,t
		integer :: i,j
		df=omega*(0.5d0*beta/cosh(0.5d0*omega*beta))**2
		db=omega*(0.5d0*beta/sinh(0.5d0*omega*beta))**2
		do i=1,n
			ff_(i)=ff(-beta*omega(i))
			fb_(i)=fb(-beta*omega(i))
		enddo

		!call get_realpart(omega,-pi*Af,tmp)
		!Gf=cmplx(tmp,pi*Af,kind=8)
		!call get_realpart(omega,-pi*AB,tmp)
		!GB=cmplx(tmp,pi*AB,kind=8)
		!call get_convolution(omega,A=-Af,B=AB,tau=[1,-1],eta=[-1,1],rt=tmp)
		!call get_realpart(omega,tmp,SEc)

		call get_realpart(omega,iSEf,tmp)
		SEf=cmplx(tmp,iSEf,kind=8)
		call get_realpart(omega,iSEB,tmp)
		SEB=cmplx(tmp,iSEB,kind=8)
		Gf=conjg(1d0/(omega-SEf))
		GB=conjg(1d0/(-1d0/Jk-SEB))


		call get_convolution(omega,A=-Gf.im/pi,B=GB.im/pi,tau=[1,-1],eta=[-1,1],rt=tmp)
		call get_realpart(omega,tmp,SEc)
		!write(233,"(A)")"Tk omega rGf iGf rGB iGB rSEf iSEf rSEB iSEB df db SEc Ac iSEc ilnGf ilnGB"
		!do i=1,n
			!write(233,"(*(e28.20))")Tk,omega(i),Gf(i),GB(i),SEf(i),SEB(i),df(i),db(i),SEc(i),Ac(i),tmp(i),imag(log(-1d0/Gf(i))),imag(log(-1d0/GB(i)))
		!enddo
		!write(233,"(x/)")
		!write(*,*)integrate(A=&
			!df*imag(log(-1d0/Gf))&
			!-df*Gf.re*SEf.im&
			!!-df*kp*(pi*Ac)*SEc&
			!+kp*db*imag(log(-1d0/GB))&
			!-kp*db*GB.re*SEB.im&
			!)/pi


		!do i=2,n/2
			!if(Gf(i).im<Gf(i-1).im) then
				!Gf(i)=Gf(i-1)
				!Gf(n-i+1)=Gf(n-i+2)
			!endif
		!enddo
		call get_realpart(omega,-Gf.im,tmp)
		Gf=cmplx(tmp,Gf.im,kind=8)
		call get_realpart(omega,-GB.im,tmp)
		GB=cmplx(tmp,GB.im,kind=8)
		call get_convolution(omega,A=Ac,B=GB.im/pi,tau=[1,1],eta=[-1,1],rt=iSEf_)
		call get_convolution(omega,A=Ac,B=Gf.im/pi,tau=[1,1],eta=[-1,-1],rt=iSEB_)
		iSEf_=-kp*iSEf_
		iSEB_=-iSEB_
		call get_realpart(omega,iSEf_,tmp)
		SEf=cmplx(tmp,iSEf_,kind=8)
		call get_realpart(omega,iSEB_,tmp)
		SEB=cmplx(tmp,iSEB_,kind=8)
		call get_convolution(omega,A=-Gf.im/pi,B=GB.im/pi,tau=[1,-1],eta=[-1,1],rt=tmp)
		call get_realpart(omega,tmp,SEc)


		call get_convolution(omega,A=-Gf.im/pi,B=GB.im/pi,tau=[1,-1],eta=[-1,1],rt=tmp)
		call get_realpart(omega,tmp,SEc)
		write(233,"(A)")"Tk omega rGf iGf rGB iGB rSEf iSEf rSEB iSEB df db SEc Ac iSEc ilnGf ilnGB"
		do i=1,n
			write(233,"(*(e28.20))")Tk,omega(i),Gf(i),GB(i),SEf(i),SEB(i),df(i),db(i),SEc(i),Ac(i),tmp(i),imag(log(-1d0/Gf(i))),imag(log(-1d0/GB(i)))
		enddo
		write(233,"(x/)")


		!iSEf=-imag(omega-1d0/Gf)
		!iSEB=-imag(-1d0/Jk-1d0/GB)

		!write(233,"(A)")"Tk omega Af Ab iSEf iSEB rGf rGB df db SEc Ac iSEc ilnGf ilnGB"
		!do i=1,n
			!write(233,"(*(e28.20))")Tk,omega(i),Gf(i).im,GB(i).im,iSEf(i),iSEB(i),Gf(i).re,GB(i).re,df(i),db(i),SEc(i),Ac(i),tmp(i),imag(log(-1d0/Gf(i))),imag(log(-1d0/GB(i)))
		!enddo
		!write(233,"(x/)")

		rt=integrate(A=&
			df*imag(log(-1d0/Gf))&
			-df*Gf.re*SEf.im&
			-df*kp*(pi*Ac)*SEc&
			+kp*db*imag(log(-1d0/GB))&
			-kp*db*GB.re*SEB.im&
			)/pi

	end function
	real(8) function entropy_analytic_T(cs) result(rt)
		real(8) :: Ac(n),Aphi(n),iSEc(n),iSEf(n),iSEB(n),iSEf_(n),iSEB_(n)
		complex(8) :: Gf(n),GB(n),SEf(n),SEB(n)
		real(8) :: df(n),db(n),SEc(n),tmp(n),ff_(n),fb_(n),u(n),Gft(n),GBt(n)
		real(8) :: A1,B1,A0,B0,t,sh,ch,lsh,lch,lgam_a,lgam_b
		integer :: i,j,cs
		df=omega*(0.5d0*beta/cosh(0.5d0*omega*beta))**2
		db=omega*(0.5d0*beta/sinh(0.5d0*omega*beta))**2

		A1=1d0
		B1=1d0
		A0=norm_Ac
		t=(gamma(a_b)*gamma(1d0+r+a_f)/(pi**2*(1d0+1d0/tan(pi*(r+a_f)/2d0)**2)*gamma(1d0+r))/(A1*B1*A0))**(1d0/(a_f+a_b))
		do i=1,n
			if(exp(beta*abs(omega(i)))>1d40) then
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
			!sh=(2d0*pi/beta)**(a_f-1d0)*abs(zgamma(a_f/2+img*beta*omega(i)/(2d0*pi)))**2*sinh(beta*omega(i)/2d0)
			!ch=(2d0*pi/beta)**(a_f-1d0)*abs(zgamma(a_f/2+img*beta*omega(i)/(2d0*pi)))**2*cosh(beta*omega(i)/2d0)
			Gf(i)=A1*t**a_f/gamma(a_f)*cmplx(1d0/tan(pi*a_f/2d0)*sh,ch,kind=8)
			sh=real(exp((a_b-1d0)*(log(2d0*pi)-log(beta))+lgam_b+lsh))*sign(1d0,omega(i))
			ch=real(exp((a_b-1d0)*(log(2d0*pi)-log(beta))+lgam_b+lch))
			!sh=(2d0*pi/beta)**(a_b-1d0)*abs(zgamma(a_b/2+img*beta*omega(i)/(2d0*pi)))**2*sinh(beta*omega(i)/2d0)
			!ch=(2d0*pi/beta)**(a_b-1d0)*abs(zgamma(a_b/2+img*beta*omega(i)/(2d0*pi)))**2*cosh(beta*omega(i)/2d0)
			GB(i)=B1*t**a_b/gamma(a_b)*cmplx(-tan(pi*a_b/2d0)*ch,sh,kind=8)
			!if(abs(omega(i))<fcut) then
				Ac(i)=A0*(abs(omega(i))/fcut)**r
			!else
				!Ac(i)=0d0
			!endif
			if(abs(omega((i)))<1.0d-12) then
				Aphi(i)=0.0d0
			elseif(abs(omega(i))<phicut) then
				Aphi(i)=B0*abs(omega(i))**(1d0-alpha)*sign(1d0,omega(i))
			else
				Aphi(i)=0d0
			endif
		enddo
		call get_convolution(omega,A=-Gf.im/pi,B=GB.im/pi,tau=[1,-1],eta=[-1,1],rt=tmp)
		call get_realpart(omega,tmp,SEc)

		SEf=omega-1d0/conjg(Gf)
		SEB=-1d0/Jk-1d0/conjg(GB)


		!call get_convolution(omega,A=-Gf.im/pi,B=GB.im/pi,tau=[1,-1],eta=[-1,1],rt=tmp)
		!call get_realpart(omega,tmp,SEc)
		!if(cs==1) then
			!write(233,"(A)")"Tk omega rGf iGf rGB iGB rSEf iSEf rSEB iSEB df db SEc Ac iSEc ilnGf ilnGB"
			!do i=1,n
				!write(233,"(*(e28.20))")Tk,omega(i),Gf(i),GB(i),SEf(i),SEB(i),df(i),db(i),SEc(i),Ac(i),tmp(i),imag(log(-1d0/Gf(i))),imag(log(-1d0/GB(i)))
			!enddo
			!write(233,"(x/)")
		!endif

		!call get_convolution(omega,A=Ac,B=GB.im/pi,tau=[1,1],eta=[-1,1],rt=iSEf_)
		!call get_convolution(omega,A=Ac,B=Gf.im/pi,tau=[1,1],eta=[-1,-1],rt=iSEB_)
		!iSEf_=-kp*iSEf_
		!iSEB_=-iSEB_
		!call get_realpart(omega,iSEf_,tmp)
		!SEf=cmplx(tmp,iSEf_,kind=8)
		!call get_realpart(omega,iSEB_,tmp)
		!SEB=cmplx(tmp,iSEB_,kind=8)



		call get_convolution(omega,A=-Gf.im/pi,B=GB.im/pi,tau=[1,-1],eta=[-1,1],rt=tmp)
		call get_realpart(omega,tmp,SEc)
		if(cs==1) then
			write(233,"(A)")"Tk omega rGf iGf rGB iGB rSEf iSEf rSEB iSEB df db SEc Ac iSEc ilnGf ilnGB"
			do i=1,n
				write(233,"(*(e28.20))")Tk,omega(i),Gf(i),GB(i),SEf(i),SEB(i),df(i),db(i),SEc(i),Ac(i),tmp(i),imag(log(-1d0/Gf(i))),imag(log(-1d0/GB(i)))
			enddo
			write(233,"(x/)")
		endif


		!if(cs==1) then
			!call get_ftau(omega,A=-Gf.im/pi,eta=-1,time=tau,rt=Gft)
			!call get_ftau(omega,A=-GB.im/pi,eta=1,time=tau,rt=GBt)
			!write(23,"(A)")"T tau Gf GB"
			!do i=1,n/2
				!write(23,"(*(e28.20))")Tk,tau(i),Gft(i),GBt(i)
			!enddo
			!write(23,"(x)")
			!do i=n/2+1,n
				!write(23,"(*(e28.20))")Tk,tau(i),Gft(i),GBt(i)
			!enddo
			!write(23,"(x/)")
		!endif

		select case(cs)
		case(1)
			rt=integrate(A=&
				df*imag(log(-1d0/Gf))&
				-df*Gf.re*SEf.im&
				!-df*kp*(pi*Ac)*SEc&
				+kp*db*imag(log(-1d0/GB))&
				-kp*db*GB.re*SEB.im&
				)/pi
		case(2)
			do i=1,n
				ff_(i)=ff(-beta*omega(i))
				fb_(i)=fb(-beta*omega(i))
			enddo
			rt=-integrate(A=&
				(2d0*ff_(:n/2)-1d0)*atan(1d0/tan(0.5d0*pi*a_f)*tanh(0.5d0*omega(:n/2)*beta))-atan(1d0/tan(0.5d0*pi*a_f))&
				+kp*((2d0*fb_(:n/2)+1d0)*atan(1d0/tan(0.5d0*pi*a_b)*tanh(0.5d0*omega(:n/2)*beta))+atan(1d0/tan(0.5d0*pi*a_b)))&
				,x=omega(1:n/2))/pi*beta
		end select
	end function
	real(8) function free_energy(Ac,Aphi,Af,AB,iSEf,iSEB,cs) result(rt)
		real(8) :: Ac(:),Aphi(:),Af(:),AB(:),iSEf(:),iSEB(:)
		integer :: cs
		complex(8) :: Gf(n),GB(n),SEf(n),SEB(n)
		real(8) :: tmp(n),ff_(n),fb_(n),SEc(n),iSEf_(n),iSEB_(n),Gft(n),SEft(n),GBt(n),Gct(n),SEBt(n),Gphit(n)
		real(8) :: A1,B1,A0,t
		integer :: i,j
		do i=1,n
			ff_(i)=ff(beta*omega(i))
			fb_(i)=fb(beta*omega(i))
		enddo
		!call get_realpart(omega,-pi*Af,tmp)
		!Gf=cmplx(tmp,pi*Af,kind=8)
		!call get_realpart(omega,-pi*AB,tmp)
		!GB=cmplx(tmp,pi*AB,kind=8)

		call get_convolution(omega,A=-Af,B=AB,tau=[1,-1],eta=[-1,1],rt=tmp)
		call get_realpart(omega,tmp,SEc)
		call get_realpart(omega,iSEf,tmp)
		SEf=cmplx(tmp,iSEf,kind=8)
		call get_realpart(omega,iSEB,tmp)
		SEB=cmplx(tmp,iSEB,kind=8)
		Gf=1d0/(omega-SEf)
		GB=1d0/(-1d0/Jk-SEB)

		! make GB a usual Green's function that satisfy kk relation
		call get_realpart(omega,GB.im,tmp)
		GB=cmplx(tmp,GB.im,kind=8)

		!call get_convolution(omega,A=Ac,B=GB.im/pi,tau=[1,1],eta=[-1,1],rt=iSEf)
		!call get_convolution(omega,A=Ac,B=Gf.im/pi,tau=[1,1],eta=[-1,-1],rt=iSEB)
		!call get_realpart(omega,iSEf,tmp)
		!SEf=cmplx(tmp,iSEf,kind=8)
		!call get_realpart(omega,iSEB,tmp)
		!SEB=cmplx(tmp,iSEB,kind=8)
		call get_ftau(omega,A=-Gf.im/pi,eta=-1,time=tau,rt=Gft)
		call get_ftau(omega,A=-GB.im/pi,eta=1,time=tau,rt=GBt)
		call get_ftau(omega,A=Ac,eta=-1,time=tau,rt=Gct)
		call get_ftau(omega,A=Aphi,eta=1,time=tau,rt=Gphit)

		!call get_convolution(omega,A=Aphi,B=Gf.im/pi,tau=[1,1],eta=[1,-1],rt=iSEf)
		!call get_convolution(omega,A=Gf.im/pi,B=Gf.im/pi,tau=[1,1],eta=[-1,-1],rt=iSEf)
		!call get_realpart(omega,iSEf,tmp)
		!SEf=cmplx(tmp,iSEf,kind=8)

		call get_ftau(omega,A=-SEf.im/pi,eta=-1,time=tau,rt=SEft)
		call get_ftau(omega,A=-SEB.im/pi,eta=1,time=tau,rt=SEBt)
		!write(23,"(A)")"T tau SEf GphiGf"
		!do i=1,n
			!write(23,"(*(e28.20))")Tk,tau(i),SEft(i),Gphit(i)*Gft(i)
		!enddo
		!write(23,"(x/)")
		if(cs==1) then
		write(*,*)integrate(A=(-Gft)*SEft,x=tau)
		write(*,*)integrate(A=-kp*GBt*SEBt,x=tau)
		write(*,*)integrate(A=-kp*(-Gft)*Gct*GBt,x=tau)
		!write(*,*)integrate(A=Gft*Gphit*Gft,x=tau)
		write(*,*)integrate(A=omega*ff_*(-Gf.im/pi))
		write(*,*)integrate(A=-ff_*imag(Gf*SEf)/pi)
		write(*,*)integrate(A=fb_*imag(GB*SEB)/pi*(-kp))
		!write(*,*)integrate(A=Gft*Gft*Gft*Gft,x=tau)
		!write(*,*)integrate(A=1d0/(exp(beta*omega)-1d0)*imag(SEf*SEf)/pi)
		!write(233,"(A)")"Tk omega SEf GphiGf"
		!!do i=1,n
			!!write(233,"(*(e28.20))")Tk,omega(i)
		!!enddo
		!!write(233,"(x/)")
		!!stop
		endif
		!write(*,*)imag(log(-1d0/Gf(1:n:n-1))+Gf(1:n:n-1)*SEf(1:n:n-1)),imag(log(-Jk/GB(1:n:n-1))+GB(1:n:n-1)*SEB(1:n:n-1))
		!write(*,*)imag(log(-1d0/Gf(1:n:n-1))),imag(log(-Jk/GB(1:n:n-1)))
		!write(*,*)integrate(A=ff_*imag(log(-1d0/Gf))/pi),integrate(A=(ff_-0.5d0)*(imag(log(-1d0/Gf))/pi+0.5d0*(1d0+sign(1d0,omega))))-1d0/beta*log(2d0)

		select case(cs)
		case(1)
			!rt=integrate(A=&
				!ff_*(imag(log(-1d0/Gf)+Gf*(omega-1d0/Gf))+0.5d0*(1d0+sign(1d0,omega))*pi)&
				!+fb_*imag(log(-1d0/GB)+GB*(-1d0/Jk-1d0/GB))&
				!!-kp*db*GB.re*SEB.im&
				!)/pi-1d0/beta*log(2d0)
			rt=(integrate(A=&
				ff_*(imag(log(-1d0/Gf)+Gf*(omega-1d0/Gf))+0.5d0*(1d0+sign(1d0,omega))*pi)&
				+kp*fb_*imag(log(-1d0/GB)+GB*(-1d0/Jk-1d0/GB))&
				!+kp*fb_*imag(log(-1d0/GB)+GB*SEB)&
				!ff_*(imag(log(-1d0/Gf))+0.5d0*(1d0+sign(1d0,omega))*pi)&
				!+kp*fb_*imag(log(-1d0/GB))&
				)/pi&
					-1d0/beta*log(2d0)&
					+integrate(A=-kp*(-Gft)*Gct*GBt,x=tau)&
						)
			!write(*,*)"fn1*******",rt
		case(2)
			rt=integrate(A=&
				ff_*(atan(Gf.re/Gf.im)-0.5d0*pi)&
				-fb_*kp*atan(GB.im/GB.re)&
				)/pi
			!write(*,*)"fn2*******",rt
			!stop
		end select
	end function
	real(8) function free_energy_analytic(cs) result(rt)
		integer :: cs
		real(8) :: Ac(n),iSEc(n),iSEf(n),iSEB(n)
		complex(8) :: Gf(n),GB(n),SEf(n),SEB(n)
		real(8) :: ff_(n),fb_(n),SEc(n),tmp(n)
		real(8) :: A1,B1,A0,t,sh,ch,lsh,lch,lgam_a,lgam_b
		integer :: i,j
		do i=1,n
			ff_(i)=ff(beta*omega(i))
			fb_(i)=fb(beta*omega(i))
		enddo

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
			!sh=(2d0*pi/beta)**(a_f-1d0)*abs(zgamma(a_f/2+img*beta*omega(i)/(2d0*pi)))**2*sinh(beta*omega(i)/2d0)
			!ch=(2d0*pi/beta)**(a_f-1d0)*abs(zgamma(a_f/2+img*beta*omega(i)/(2d0*pi)))**2*cosh(beta*omega(i)/2d0)
			Gf(i)=A1*t**a_f/gamma(a_f)*cmplx(1d0/tan(pi*a_f/2d0)*sh,-ch,kind=8)
			sh=real(exp((a_b-1d0)*(log(2d0*pi)-log(beta))+lgam_b+lsh))*sign(1d0,omega(i))
			ch=real(exp((a_b-1d0)*(log(2d0*pi)-log(beta))+lgam_b+lch))
			!sh=(2d0*pi/beta)**(a_b-1d0)*abs(zgamma(a_b/2+img*beta*omega(i)/(2d0*pi)))**2*sinh(beta*omega(i)/2d0)
			!ch=(2d0*pi/beta)**(a_b-1d0)*abs(zgamma(a_b/2+img*beta*omega(i)/(2d0*pi)))**2*cosh(beta*omega(i)/2d0)
			GB(i)=B1*t**a_b/gamma(a_b)*cmplx(-tan(pi*a_b/2d0)*ch,-sh,kind=8)
			!if(abs(omega(i))<fcut) then
				Ac(i)=A0*(abs(omega(i))/fcut)**r
			!else
				!Ac(i)=0d0
			!endif
		enddo
		call get_convolution(omega,A=-Gf.im/pi,B=GB.im/pi,tau=[1,-1],eta=[-1,1],rt=tmp)
		call get_realpart(omega,tmp,SEc)

		SEf=omega-1d0/Gf
		SEB=-1d0/Jk-1d0/GB

		!call get_convolution(omega,A=Ac,B=GB.im/pi,tau=[1,1],eta=[-1,1],rt=iSEf)
		!call get_convolution(omega,A=Ac,B=Gf.im/pi,tau=[1,1],eta=[-1,-1],rt=iSEB)
		!iSEf=-kp*iSEf
		!iSEB=-iSEB
		
		select case(cs)
		case(1)
			!rt=integrate(A=&
				!(0.5d0-ff_)*imag(log(-1d0/Gf))&
				!-(0.5d0-ff_)*Gf.re*SEf.im&
				!!-(0.5d0-ff_)*kp*(pi*Ac)*SEc&
				!-(0.5d0+fb_)*kp*imag(log(-1d0/GB))&
				!+(0.5d0+fb_)*kp*GB.re*SEB.im&
				!)/pi
			rt=integrate(A=&
				!ff_*(imag(log(-1d0/Gf)+Gf*(omega-1d0/Gf))+0.5d0*(1d0+sign(1d0,omega))*pi)&
				!+fb_*imag(log(-1d0/GB)+GB*(-1d0/Jk-1d0/GB))&
				!ff_*(imag(log(-1d0/Gf))+0.5d0*(1d0+sign(1d0,omega))*pi)&
				!+kp*fb_*imag(log(-1d0/GB))&
				ff_*(imag(log(-1d0/Gf))+0.5d0*(1d0+sign(1d0,omega))*pi)&
				+kp*fb_*imag(log(-1d0/GB))&
				)/pi-1d0/beta*log(2d0)
			!write(*,*)"1********",rt
		case(2)
			rt=integrate(A=&
				ff_*(atan(Gf.re/Gf.im)-0.5d0*pi)&
				-fb_*kp*atan(GB.im/GB.re)&
				)/pi
			!write(*,*)"2********",rt
		end select
	end function
end module
program main
	use global
	implicit none
	real(8) :: Af(n),AB(n),Ac(n),Aphi(n),Ac0(n),Aphi0(n),rhoc(n),rhophi(n),rGf(n),rGB(n),rSEF(n),rSEB(n),iSEF(n),iSEB(n),SUS(n),rSUS(n),Tmat(n),SUSt(n),Tmatt(n),Gft(n),GBt(n),Gphit(n),Gct(n),SEft(n),dTk,omg(2),fe(2)=0d0,tmp,mix
	complex(8) :: dat(2,2)
	integer :: i,j,e
	logical :: conv

	call omp_set_nested(.false.)
	!call omp_set_max_active_levels(1)
	call omp_set_dynamic(.false.)
	call mkl_set_dynamic(0)
	!call omp_set_schedule(omp_sched_static,0)

	call mkl_set_num_threads(24)
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
	call set_omega(omega(1:nl*nlog*2),base,nl,10d0)
	!!call add_omega(omega(1:nl*nlog*2+nadd1*4),fcut+sqrt(abs(log(0.2d0))/fbeta),sqrt(abs(log(0.2d0))/fbeta)*2d0,nadd1)
	!!call add_omega(omega,0.8d0,0.16d0,nadd2)
	!call set_omega(omega(1:nl*nlog*2),1.5d0,15,10d0)
	call add_omega(omega(1:nl*nlog*2+sum(nadd(1:1))*4),1d0,0.1d0,nadd(1))
	!call add_omega(omega(1:nll*2+sum(nadd(1:1))*4),2d0,1d0,nadd(1))
	!call add_omega(omega(1:nll*2+sum(nadd(1:1))*4),1d0,0.3d0,nadd(1))
	!call add_omega(omega(1:nll*2+nadd1*4),1.5d0,0.5d0,nadd1)
	call add_omega(omega(1:nl*nlog*2+sum(nadd(1:2))*4),0.8d0,0.16d0,nadd(2))
	!call add_omega(omega(1:nll*2+sum(nadd(1:3))*4),9.8d0,0.09d0,nadd(3))
	write(*,"('n: ',i4,', max freq: ',es10.2,', min freq: ',es10.2)")n,maxval(abs(omega)),minval(abs(omega))

	!iSEf=-1d0
	!iSEB(:size(omega)/2)=-1d0
	!iSEB(size(omega)/2+1:)=1d0

	!do i=1,n
	!write(22,"(*(e28.20))")omega(i),Af(i)*pi,Ab(i)*pi,iSEf(i),iSEB(i),SUS(i)/pi,Tmat(i)*pi,A_c(omega(i))*pi,merge(0d0,A_phi(omega(i))*pi,abs(A_phi(omega(i))*pi)<1d-20)
	!enddo
	!write(22,"(x/)")
	!stop

	Tk=0.09d0
	!Tk=0.01d0

	!fbeta=1d0
	!Jk=Jk*integrate(Af=A_c)/pi
	!norm_Ac=norm_Ac/integrate(Af=A_c)
	!call self_consistent(1d-5,6000,3d-3,Af,AB,iSEf,iSEB)
	!do i=1,n
	!write(22,"(*(e28.20))")omega(i),Af(i)*pi,Ab(i)*pi,iSEf(i),iSEB(i),SUS(i)/pi,Tmat(i)*pi,A_c(omega(i))*pi,merge(0d0,A_phi(omega(i))*pi,abs(A_phi(omega(i))*pi)<1d-20)
	!enddo
	!write(22,"(x/)")
	!fbeta=7d0
	!Jk=1d0
	!norm_Ac=1d0

	!norm_Aphi=norm_Aphi/(integrate(Af=A_phi,l=size(omega)/2)*2)
	!fcut=abs(omega(find_sidx(omega,-fcut)))
	!Jk=Jk*integrate(Af=A_c)/pi
	!norm_Ac=norm_Ac/integrate(Af=A_c)*Jk
	!norm_Ac=norm_Ac/integrate(Af=A_c,l=find_sidx(omega,-fcut),u=find_sidx(omega,fcut))

	!norm_Ac=(r+1d0)/(2d0*fcut)

	!do i=1,n
		!read(10,"(4(e28.20))")omega(i),iSEf(i),iSEB(i),Ac0(i)
	!enddo
	!iSEf=-iSEf
	!iSEB=iSEB
	!Ac0=Ac0/pi
	!rhoc=Ac0
	!rhophi=0d0
	!Aphi0=0d0

	!beta=1d-11
	!write(*,*)-exp(-beta)/(exp(-beta)-1d0),1d0/(exp(-beta)-1d0),-exp(-beta)/(exp(-beta)-1d0)+1d0/(exp(-beta)-1d0)
	!write(233,"(A)")"omega fb"
	!do i=1,n
		!write(233,"(*(e28.20))")omega(i),1d0/(exp(11.111111111111111111d0*omega(i))-1)
	!enddo
	!stop

	write(25,"(A)")"T J Sn Sa"
	do 
		is_expinit=.true.
		norm_Ac=1d0
		do i=1,n
			Ac0(i)=A_c(omega(i))
			!if(abs(omega(i))<=1d0) then
			!Ac0(i)=sqrt(1d0-omega(i)**2)
			!else
			!Ac0(i)=0d0
			!endif
			!Ac0(i)=exp(-omega(i)**2)/sqrt(pi)
			Aphi0(i)=A_phi(omega(i))
		enddo
		norm_Ac=1d0/integrate(A=Ac0)
		Ac0=Ac0*norm_Ac
		rhoc=Ac0
		rhophi=Aphi0

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
					iSEf(i)=imag(dat(1,1))
					iSEB(i)=imag(dat(1,2))
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
		write(26,"(A)")"T Fa1 F1 Fa2 F2 pT pFa1 pF1 pFa2 pF2"
		write(27,"(A)")"T A0 AfAB"
		do
			read(13,*,iostat=e)Tk
			if(e==-1) exit
			write(*,"(*(A4,es9.2))")"Tk:  ",Tk," Jk:",Jk," kp:",kp," r: ",r," g: ",g

			!tau=(omega(n/2+1:)+pi*Tk/10d0)*10d0
			!tau(1:n/2)=asin(pi*Tk/(omega(1:n/2)/omega(1)*(100000d0-pi*Tk)+pi*Tk))/(pi*Tk)
			!tau(n/2+1:)=1d0/Tk-tau(1:n/2)
			tau(1:n/2)=omega(n/2+1:n)/(3d0*omega(n)-omega(n-1))/Tk
			tau(n/2+1:)=1d0/Tk-tau(n/2:1:-1)


			conv=self_consistent(1d-5,[2000,1],3d-3,iSEf,iSEB,rhoc,rhophi,Ac0,Aphi0,Af,AB,Ac,Aphi)
			if(.not.conv) then
				exit
			endif
			write(26,"(es28.20$)")Tk
			write(26,"(es28.20$)")free_energy_analytic(1),free_energy(rhoc,rhophi,Af,AB,iSEf,iSEB,1),free_energy_analytic(2),free_energy(rhoc,rhophi,Af,AB,iSEf,iSEB,2)
			write(26,"(x)")

			write(27,"(*(es28.20))")Tk,Af(n/2),Af(n/2+n/7)*AB(n/2+n/7)
			write(25,"(es28.20$)")Tk,Jk
			!tmp=entropy_coleman(rhoc,Af,AB)
			tmp=entropy_coleman(rhoc,iSEf,iSEB)
			write(*,*)"entropy: ",tmp
			write(25,"(es28.20$)")tmp
			tmp=entropy_analytic_T(1)
			write(*,*)"entropy analytic T: ",tmp
			write(25,"(es28.20)")tmp
			!tmp=entropy_analytic_T(2)
			!write(*,*)"entropy thesis: ",tmp
			!write(25,"(es28.20)")tmp

			!if(Tk<1.1d-9) then
			!stop
			!endif

			call get_convolution(omega,A=Af,B=Af,tau=[1,-1],eta=[-1,-1],rt=SUS)
			call get_convolution(omega,A=Af,B=AB,tau=[1,1],eta=[-1,1],rt=Tmat)

			call get_realpart(omega,-Af*pi,rGf)
			call get_realpart(omega,-AB*pi,rGB)


			write(22,"(A)")"Tk omega Af Ab Ac0 Ac Aphi0 Aphi iSEf iSEB SUS Tmat rGf rGB atan"
			do i=1,n
				write(22,"(*(e28.20))")Tk,omega(i),Af(i)*pi,Ab(i)*pi,Ac0(i),Ac(i),merge(0d0,Aphi0(i),abs(Aphi0(i))<1d-20),Aphi(i),iSEf(i),iSEB(i),SUS(i)/pi,Tmat(i)*pi,rGf(i),rGB(i),atan(tanh((omega(i)/Tk)*0.5d0)/tan((1d0-0.785d0)*pi*0.5d0))
			enddo
			write(22,"(x/)")

			call get_ftau(omega,A=-SUS/pi,eta=1,time=tau,rt=SUSt)
			call get_ftau(omega,A=-Tmat/pi,eta=-1,time=tau,rt=Tmatt)
			call get_ftau(omega,A=Af,eta=-1,time=tau,rt=Gft)
			call get_ftau(omega,A=AB,eta=1,time=tau,rt=GBt)
			call get_ftau(omega,A=rhoc,eta=-1,time=tau,rt=Gct)
			call get_ftau(omega,A=rhophi,eta=1,time=tau,rt=Gphit)
			call get_ftau(omega,A=iSEf/pi,eta=-1,time=tau,rt=SEft)

			write(23,"(A)")"T tau Gf GB Gc Gphi SUS Tmat GfGf GfGB GfSEf"
			do i=1,n/2
				write(23,"(*(e28.20))")Tk,tau(i),Gft(i),GBt(i),Gct(i),Gphit(i),SUSt(i),Tmatt(i),-Gft(i)*Gft(i),Gft(i)*GBt(i),Gft(i)*SEft(i)
			enddo
			do i=n/2+1,n
				write(23,"(*(e28.20))")Tk,tau(i),Gft(i),GBt(i),Gct(i),Gphit(i),SUSt(i),Tmatt(i),-Gft(i)*Gft(i),Gft(i)*GBt(i),Gft(i)*SEft(i)
			enddo
			write(23,"(x/)")
			!stop

			!write(23,"(A)")"T tau SUS Tmat Gf"
			!do i=1,n/2
				!write(23,"(*(e28.20))")Tk,tau(i),SUSt(i),Tmatt(i),Gft(i)
			!enddo
			!write(23,"(x)")
			!do i=n/2+1,n
				!write(23,"(*(e28.20))")Tk,tau(i),SUSt(i),Tmatt(i),Gft(i)
			!enddo
			!write(23,"(x/)")

			!write(24,"(*(e28.20))")Tk,rSUS(n/2)
			if(is_expinit) then
				call get_realpart(omega,iSEf,rSEf)
				call get_realpart(omega,iSEB,rSEB)
				write(*,*)"export initial...."
				rewind(10)
				do i=1,n
					write(10,"(*(e28.20))")omega(i),rSEf(i),iSEf(i),rSEB(i),iSEB(i)
				enddo
				is_expinit=.false.
			endif
			pTk=Tk
			iG0=Af*pi
			rG0=rGf
			fe(1)=fe(2)
		enddo
		stop
		Jk=Jk+0.1d0
		if(Jk>2.2d0-1d-4) then
			exit
		endif
		!kp=kp+0.2d0
		!if(kp>0.8d0) then
		!exit
		!endif
		!fbeta=fbeta+1d0
		!if(fbeta>8d0) then
		!exit
		!endif
		!g=g+0.5d0
		!if(g>5d0) then
		!exit
		!endif
	enddo
end program
