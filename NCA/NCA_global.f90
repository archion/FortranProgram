include "../../DGC/DGC_mod.f90"
module global
	use param
	use M_DGC
	implicit none
	integer, parameter :: in_init=10,out_Gw=20,out_T=30,out_Gt=40,out_init=50,out_info=60,out_param=70,xx_save=1,xx_diff=2,xx_new=-1
	integer, parameter :: omega_grid=1,tau_grid=2,Gc0=3,Gd=4,Gc=5,Q=6,SEc=7,lambda=8,Gf=9,rhoc=10,Gb=11,tSEf=12,tSEb=13,Gft=14,Gbt=15,PMT=16,output=17,SEf=18,SEb=19,tlambda=20,SC1=-1,SC2=-2
	logical :: is_expinit=.true.,is_bathT=.false.
	type gfs
		real(wp) :: rhoc(n),lambda,Q(1)
		complex(wp) :: Gf(n),Gb(n),Gc(n),Gd(n),Gc0(n)
		complex(wp) :: tSEf(n),tSEb(n),SEf(n),SEb(n),SEc(n),SEd(n)
		real(wp) :: Gft(n),Gbt(n)
		logical :: conv=.false.
	contains
		procedure :: io
	end type
	type(gfs) :: self
    ! graph
	integer :: niter(-5:5)
	real(8) :: iter_err(-5:5)=1._wp,iter_rate(-5:-1)
	real(8), allocatable :: x(:),x_(:)
contains
	function compute_(node) result(rt)
		integer :: node
		character(:), allocatable :: rt
		integer, save :: nx,cnd_=0
		integer :: i,j
		logical :: is_
		real(wp) :: rho(n),norm,W
		select case(node)
		case(PMT)
			rt="PMT"; if(debug) return
		case(output)
			rt="output"; if(debug) return
			!if(abs(Tk-0.1_wp)<1e-6_wp) then
				!call self%io([out_Gw,out_Gt],[out_Gw,out_Gt])
			!endif
		case(omega_grid)
			rt="omega_grid"; if(debug) return
			call set_omega_grid()
			do i=1,n-1
				if(omega(i)>=omega(i+1)) then
					write(*,*)"omega error"
					stop
				endif
			enddo
		case(tau_grid)
			rt="tau_grid"; if(debug) return
			tau(1:n/2)=omega(n/2+1:n)/(2.0001_wp*omega(n))/Tk
			tau(n/2+1:)=1._wp/Tk-tau(n/2:1:-1)
		case(rhoc)
			rt="rhoc"; if(debug) return
			do i=1,n
				rho(i)=A_c(omega(i))
			enddo
			!norm=1._wp/integrate(omega,A=rho)
			!self%rhoc=rho*norm
			self%rhoc=rho
			self%Gc0%im=-pi*self%rhoc
			call set_realpart(omega,self%Gc0)
		case(Gf)
			rt="Gf"; if(debug) return
			!$OMP PARALLEL DO
			do i=1,n
				self%Gf(i)%im=self%tSEf(i)%im/((omega(i)+self%lambda-ed-self%SEf(i)%re)**2+(self%SEf(i)%im)**2)
			enddo
			!$OMP END PARALLEL DO
		case(Gb)
			rt="Gb"; if(debug) return
			!$OMP PARALLEL DO
			do i=1,n
				self%Gb(i)%im=self%tSEb(i)%im/((omega(i)+self%lambda-self%SEb(i)%re)**2+(self%SEb(i)%im)**2)
			enddo
			!$OMP END PARALLEL DO
		case(Q)
			rt="Q"; if(debug) return
			call get_Gtau(omega,G=(2._wp*self%Gf+2._wp*kp*self%Gb),eta=-1,time=[-eps*1e3_wp],rt=self%Q)
		case(lambda)
			rt="lambda"; if(debug) return
		case(tlambda)
			rt="tlambda"; if(debug) return
			self%lambda=min(0.5_wp*(self%SEf(n/2)%re+self%SEf(n/2+1)%re)+ed,0.5_wp*(self%SEb(n/2)%re+self%SEb(n/2+1)%re))
		case(tSEf)
			rt="tSEf"; if(debug) return
			call get_convolution(omega,A=2._wp*kp*V2*self%Gb,B=self%Gc0,tau=[1,1],rt=self%tSEf)
		case(tSEb)
			rt="tSEb"; if(debug) return
			call get_convolution(omega,A=2._wp*V2*self%Gf,B=self%Gc0,tau=[1,-1],rt=self%tSEb)
		case(SEf)
			rt="SEf"; if(debug) return
			self%SEf%im=self%tSEf%im*ff(-beta*omega)
			call set_realpart(omega,self%SEf)
		case(SEb)
			rt="SEb"; if(debug) return
			self%SEb%im=self%tSEb%im*ff(-beta*omega)
			call set_realpart(omega,self%SEb)
		case(Gd)
			rt="Gd"; if(debug) return
			call get_convolution(omega,A=self%Gf/self%Q(1),B=self%Gb,tau=[1,-1],rt=self%Gd,two=underscore)
			call set_realpart(omega,self%Gd)
		case(SEc)
			rt="SEc"; if(debug) return
			self%SEc=1._wp/(1._wp/(V2*self%Gd)+self%Gc0)
		case(Gc)
			rt="Gc"; if(debug) return
			self%Gc=1._wp/(1._wp/self%Gc0-self%SEc)
		case(Gft)
			rt="Gft"; if(debug) return
			call get_Gtau(omega,G=self%Gf,eta=-1,time=tau,rt=self%Gft)
		case(Gbt)
			rt="GBt"; if(debug) return
			call get_Gtau(omega,G=self%Gb,eta=1,time=tau,rt=self%GBt)
		case(:-1)
			rt="SC"; if(debug) return
			if(allocated(x)) then
				deallocate(x_,x)
			endif
			do
				nx=0
				call x_x(onode(node,1:onode(node,0)),xx_save)
				if(allocated(x)) then
					exit
				else
					allocate(x_(nx),x(nx))
				endif
			enddo
			cnd_=node
			
			if(iter_exit(node)==1) then
				write(*,"('SC=',i2,' T=',es10.3,', V2: ',es9.2,', ed: ',es9.2,', kp: ',es9.2,', r: ',es9.2,' initial...')")node,Tk,V2,ed,kp,r
				call mbroyden(0,x,x_,iter_rate(node),nsave=20,id=abs(node))
			endif
			write(*,"('SC=',i2,' T=',es10.3,' err=',es10.3,' iter=',i4,' lambda=',es10.3,' Q=',es10.3)")node,Tk,iter_err(abs(node)),iter_exit(node),self%lambda,self%Q(1)
			if(iter_exit(node)==niter(node).or.(iter_err(abs(node))<=iter_err(node))) then
				call mbroyden(huge(1),x,x_,id=abs(node))
				deallocate(x,x_)
				write(*,"('SC=',i2,' T=',es10.3,', V2: ',es9.2,', ed: ',es9.2,', kp: ',es9.2,', r: ',es9.2,' finish...')")node,Tk,V2,ed,kp,r
				cnd_=0
				iter_exit(node)=merge(0,-iter_exit(node),iter_err(abs(node))<=iter_err(node))
				iter_err(abs(node))=1._wp
				if(iter_exit(node)<0) then
					write(*,*)"not convergent!!!!!!"
				endif
			endif
			nx=0
		case default
			write(*,*)"case not implimented",node
			stop
		end select
		if(cnd_<0.and.any(onode(cnd_,1:onode(cnd_,0))==node)) then
			call x_x([node],xx_diff)
			if(nx==size(x)) then
				call x_x(onode(cnd_,1:onode(cnd_,0)),xx_new)
				cnd_=0
			endif
		endif
	contains
		subroutine x_x(nd,flag)
			integer :: nd(:),flag
			integer :: l
			if(flag==xx_new) then
				iter_err(abs(cnd_))=maxval(abs(x_))
				!iter_err(abs(cnd_))=maxval(abs(x_)/merge(abs(x),1._wp,abs(omega(1:n))<2.*fcut.and.abs(x)<1e-4_wp))
				!iter_err(abs(cnd_))=maxval(abs(x_)/merge(abs(x),1._wp,abs(x)<1e-4_wp))
				!iter_err(abs(cnd_))=maxval(abs(x_)/merge(1._wp,abs(x),abs(x)<1e-6_wp))
				!iter_err(abs(cnd_))=maxval(abs(x_)/merge(1e-6_wp,abs(x),abs(x)<1e-6_wp))
				if(niter(cnd_)>0) then
					call mbroyden(iter_exit(cnd_),x,x_,id=abs(cnd_))
					!x=x_+x
				else
					x=x+iter_rate(cnd_)*x_
				endif
				nx=0
			endif
			do l=1,size(nd)
				select case(nd(l))
				case(tSEf)
					! SEf self-consistent
					if(allocated(x)) then
						select case(flag)
						case(xx_save)
							x(nx+1:nx+size(self%tSEf))=self%tSEf%im
						case(xx_diff)
							x_(nx+1:nx+size(self%tSEf))=self%tSEf%im-x(nx+1:nx+size(self%tSEf))
						case(xx_new)
							self%tSEf%im=x(nx+1:nx+size(self%tSEf))
						end select
					endif
					nx=nx+size(self%tSEf)
				case(tSEb)
					! SEb self-consistent
					if(allocated(x)) then
						select case(flag)
						case(xx_save)
							x(nx+1:nx+size(self%tSEb))=self%tSEb%im
						case(xx_diff)
							x_(nx+1:nx+size(self%tSEb))=self%tSEb%im-x(nx+1:nx+size(self%tSEb))
						case(xx_new)
							self%tSEb%im=x(nx+1:nx+size(self%tSEb))
							call set_realpart(omega,self%tSEb)
						end select
					endif
					nx=nx+size(self%tSEb)
				case(lambda)
					if(allocated(x)) then
						select case(flag)
						case(xx_save)
							x(nx+1)=self%lambda
						case(xx_diff)
							!x_(nx+1)=ff(self%Q(1)-1._wp)-0.5_wp
							x_(nx+1)=(self%Q(1)-1._wp)
						case(xx_new)
							if(abs(self%lambda-x(nx+1))>max(0.1_wp*Tk,0._wp)) then
								self%lambda=self%lambda-max(0.1_wp*Tk,0._wp)*sign(1._wp,(self%Q(1)-1._wp))
							else
								self%lambda=x(nx+1)
							endif
						end select
					endif
					nx=nx+1
				case default
					write(*,*)"not complete for self consistent: ",node
					stop
				end select
			enddo
		end subroutine
	end function
	subroutine io(self,ut,flag)
		class(gfs) :: self
		integer :: ut(:)
		integer :: flag(:)
		integer :: i,k,e
		real(wp) :: omg(2)
		real(wp) :: dat(2)
		logical :: isopen
		do k=1,size(flag)
			select case(flag(k))
			case(out_Gw)
				write(ut(k),"(A)")"T omega omegaT ff rGf iGf rGb iGb rSEf iSEf rSEb iSEb rGc0 iGc0 rSEc iSEc rGc iGc rGd iGd"
				do i=1,n
					write(ut(k),"(*(e29.20e3))")Tk,omega(i),omega(i)/Tk,ff(-beta*omega(i)),self%Gf(i),self%Gb(i),self%SEf(i),self%SEb(i),self%Gc0(i),self%SEc(i),self%Gc(i),self%Gd(i)
				enddo
				write(ut(k),"(x)")
			case(out_T)
				write(ut(k),"(*(e29.20e3))")Tk,Tkunit,self%lambda,self%Q(1),(self%Gf(n/2+1)+self%Gf(n/2))/2.,(self%GB(n/2+1)+self%GB(n/2))/2.,(self%SEc(n/2+1)+self%SEc(n/2))/2.
			case(out_Gt)
				write(ut(k),"(A)")"T tau Gft Gbt"
				do i=1,n/2
					write(ut(k),"(*(e29.20e3))")Tk,tau(i),self%Gft(i),self%GBt(i)
				enddo
				write(ut(k),"(x)")
			case(out_init)
				write(*,*)"export initial...."
				rewind(ut(k))
				is_expinit=.false.
				do i=1,n
					write(ut(k),"(*(e29.20e3))")omega(i),self%SEf(i)%im
				enddo
			case(in_init)
				i=1
				rewind(ut(k))
				read(ut(k),"(*(e29.20e3))")omg(1),dat(1)
			o:	do 
					read(ut(k),"(*(e29.20e3))",iostat=e)omg(2),dat(2)
					if(e==-1) then
						backspace(ut(k))
						backspace(ut(k))
						backspace(ut(k))
						read(ut(k),"(*(e29.20e3))")omg(1),dat(1)
						dat(2)=dat(2)+(dat(2)-dat(1))/(omg(2)-omg(1))*(-omega(1)-omg(2))
						omg(2)=-omega(1)
						read(ut(k),"(*(e29.20e3))")omg(1),dat(1)
					endif
					do 
						if(omg(2)-omega(i)>=0._wp) then
							dat(1)=dat(1)+(dat(2)-dat(1))/(omg(2)-omg(1))*(omega(i)-omg(1))
							self%SEf(i)%im=dat(1)
							omg(1)=omega(i)
							i=i+1
							if(i==n+1) exit o
						else
							omg(1)=omg(2)
							dat(1)=dat(2)
							exit
						endif
					enddo
				enddo o
				call set_realpart(omega,self%SEf)
			case(out_info)
				write(*,"('fcut: ',es9.2,', fbeta: ',es9.2,', bcut: ',es9.2)")fcut,fbeta,bcut
				write(*,"('nlog: ',i5,', nl: ',i5,', nadd: ',*(i5))")nlog,nl,nadd
				write(*,"('n: ',i5,', max freq: ',es13.5,', min freq: ',es13.5)")n,maxval(abs(omega)),minval(abs(omega))
				do i=10,1000
					if(i==in_init) cycle
					inquire(unit=i, opened=isopen)
					if(isopen) then
						write(i,"('#fcut: ',es9.2,', fbeta: ',es9.2,', bcut: ',es9.2)")fcut,fbeta,bcut
						write(i,"('#nlog: ',i5,', nl: ',i5,', nadd: ',*(i5))")nlog,nl,nadd
						write(i,"('#n: ',i5,', max freq: ',es13.5,', min freq: ',es13.5)")n,maxval(abs(omega)),minval(abs(omega))
						if(i==out_T) then
							write(i,"(A)")"T unit lambda Q rGf0 iGf0 rGB0 iGB0 rSEc iSEc"
						endif
					endif
				enddo
			case(out_param)
				write(*,"('V2: ',es9.2,', ed: ',es9.2,', kp: ',es9.2,', r: ',es9.2)")V2,ed,kp,r
				do i=10,1000
					if(i==in_init) cycle
					inquire(unit=i, opened=isopen)
					if(isopen) then
						write(i,"('#V2: ',es9.2,', ed: ',es9.2,', kp: ',es9.2,', r: ',es9.2)")V2,ed,kp,r
					endif
				enddo
			end select
		end do
	end subroutine
	real(wp) elemental function A_c(x)
		real(wp), intent(in) :: x
		if(abs(x)<fcut) then
			if(isnan(gap)) then
				!A_c=(abs(x)/fcut)**r
				A_c=0.5_wp
			else
				if(abs(x)<gap) then
					A_c=0._wp
				else
					A_c=abs(x)/sqrt(x**2-gap**2)
				endif
			endif
		!elseif(abs(x)<=(-omega(1))) then
			!A_c=exp(-(x**2-fcut**2)/pi*fbeta)
		else
			A_c=0._wp
		endif
		A_c=exp(-x**2/pi)
	end function
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
	subroutine get_convolution(omega,A,B,tau,rt,two)
		real(wp) :: omega(:)
		complex(wp) :: rt(:)
		complex(wp) :: A(:),B(:)
		integer :: tau(2)
		logical, optional :: two
		integer :: i,j,k,n,iwb(size(omega)),ibw(size(omega)),m
		real(wp) :: e1,e2,omg,M_phi
		complex(wp) :: A1,A2,B1,B2,rt_
		rt=0._wp
		n=size(omega)
		!if(any(isnan(abs(A)))) then
			!write(*,*)"A nan"
		!endif
		!if(any(isnan(abs(B)))) then
			!write(*,*)"B nan"
		!endif
		!$OMP PARALLEL DO PRIVATE(omg,ibw,iwb,A1,B1,A2,B2,e1,e2,m,rt_)
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
			B1=merge(B(i)+(tau(2)*(omg-tau(1)*omega(1))-omega(i))/(omega(i+1)-omega(i))*(B(i+1)-B(i)),cmplx(0._wp,kind=wp),i>0.and.i<n)
			e1=omega(1)
			rt_=0._wp
			do i=1,n-1
				do j=iwb(i),iwb(i+1),sign(1,iwb(i+1)-iwb(i))
					if(j==iwb(i+1)) then
						e2=omega(i+1)
						A2=A(i+1)
						if(j>0.and.j<n) then
							B2=B(j)+(tau(2)*(omg-tau(1)*e2)-omega(j))/(omega(j+1)-omega(j))*(B(j+1)-B(j))
						else
							B2=0._wp
						endif
					else
						m=j+(1-tau(1)*tau(2))/2
						e2=tau(1)*(omg-tau(2)*omega(m))
						B2=B(m)
						m=ibw(m)
						if(m>0.and.m<n) then
							A2=A(m)+(e2-omega(m))/(omega(m+1)-omega(m))*(A(m+1)-A(m))
						else
							A2=0._wp
						endif
					endif
					rt_=rt_+0.5_wp*(e2-e1)*cmplx(0._wp,&
						!(imag(A1)*imag(B1)*((1._wp+exp(-beta*omg))/((1._wp+exp(-beta*e1))*(1._wp+exp(beta*(e1-omg))))+merge(f(e1)*f(omg-e1),0._wp,present(two))))+&
						!(imag(A2)*imag(B2)*((1._wp+exp(-beta*omg))/((1._wp+exp(-beta*e2))*(1._wp+exp(beta*(e2-omg))))+merge(f(e2)*f(omg-e2),0._wp,present(two)))),kind=wp)
						!(imag(A1)*imag(B1)*(f(-e1)*f(e1-omg)/f(-omg)+merge(f(e1)*f(omg-e1),0._wp,present(two))))+&
						!(imag(A2)*imag(B2)*(f(-e2)*f(e2-omg)/f(-omg)+merge(f(e2)*f(omg-e2),0._wp,present(two)))),kind=wp)
						!(imag(A1)*imag(B1)*(f(-e1)*f(e1-omg)+merge(f(e1)*f(omg-e1),0._wp,present(two))))+&
						!(imag(A2)*imag(B2)*(f(-e2)*f(e2-omg)+merge(f(e2)*f(omg-e2),0._wp,present(two)))),kind=wp)
						(imag(A1)*imag(B1)*f(e1,omg)+&
						imag(A2)*imag(B2)*f(e2,omg)),kind=wp)

					
					A1=A2
					B1=B2
					e1=e2
				enddo
			enddo
			rt(k)=rt(k)-rt_/pi
		enddo
		!$OMP END PARALLEL DO
	contains
		real(wp) function f(x,y)
			real(wp) :: x,y
			if(present(two)) then
				if(y>0._wp) then
					f=(1._wp+exp(-beta*y))/(1._wp+exp(-beta*x)+exp(beta*(x-y))+exp(-beta*y))
				else
					f=(1._wp+exp(beta*y))/(1._wp+exp(beta*x)+exp(beta*(y-x))+exp(beta*y))
				endif
			else
				if(y<0._wp) then
					f=1._wp/(1._wp+(exp(-beta*(x-y))+exp(beta*(x)))/(1._wp+exp(beta*y)))
				else
					f=1._wp/(1._wp+(exp(-beta*x)+exp(beta*(x-y)))/(1._wp+exp(-beta*y)))
				endif
			endif
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
					G(i)%re=G(i)%re+(dx/(omega(j)-omega(i)))*(G(j)%im-G(i)%im)
				endif
			enddo
		enddo
		do i=1,n
			G(i)%re=0.5_wp/pi*G(i)%re
		enddo
	end subroutine
	function fit(y,x,nx) result(rt)
		real(wp) :: y(:),rt(size(y))
		real(wp) :: x(:),nx(:)
		integer :: i,n
		n=1
		i=0
		do 
			i=i+1
			if(i>size(nx)) exit
			if(n==size(x)) then
				rt(i)=y(n)
			elseif(nx(i)<=x(n)) then
				rt(i)=y(n)
			elseif(nx(i)>x(n).and.nx(i)<=x(n+1)) then
				rt(i)=y(n)+(y(n+1)-y(n))/(x(n+1)-x(n))*(nx(i)-x(n))
			elseif(nx(i)>x(n+1)) then 
				n=n+1
				i=i-1
			endif
		enddo
	end function
end module
