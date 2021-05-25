include "../../DGC/DGC_mod.f90"
module global
	use param
	use M_DGC
	implicit none
	integer, parameter :: in_init=10,out_Gw=20,out_T=30,out_Gt=40,out_init=50,out_info=60,out_param=70,xx_save=1,xx_diff=2,xx_new=-1
	integer, parameter :: omega_grid=1,tau_grid=2,Gc0=3,Gd=4,Gc=5,Q=6,SEc=7,lambda=8,Gf=9,rhoc=10,Gb=11,tSEf=12,tSEb=13,Gft=14,Gbt=15,PMT=16,output=17,SEf=18,SEb=19,tlambda=20,Fimp=21,Ga=22,SEa=23,tSEa=24,lattice=25,optcond=26,rhophi=27,SC1=-1,SC2=-2
	logical :: is_expinit=.true.,is_bathT=.false.
	type gfs
		real(wp) :: rhoc(n),rhophi(n),rhov(size(omega_ek)),lambda,Q(1)
		complex(wp) :: Gf(n,size(ed),size(ed)),Gb(n),Ga(n)=0._wp,Gc(n,2),Gd(n),Fd(n),Gc0(n),Gphi0(n),Fc0(n),Fc(n)
		complex(wp) :: tSEf(n,size(ed),size(ed)),tSEb(n),tSEa(n)=0._wp,SEf(n,size(ed),size(ed)),SEb(n),SEa(n)=0._wp,SEc(n),SEd(n),sigma(n)
		real(wp) :: Gft(n),Gbt(n),Fimp
		real(wp), allocatable :: uk(:,:)
		logical :: conv=.false.
	contains
		procedure :: io
	end type
	type(gfs) :: self
    ! graph
	integer :: niter(-5:5)
	real(8) :: iter_err(-5:5)=1._wp,iter_rate(-5:-1)
	real(8), allocatable :: x(:),x_(:)
	real(wp) :: k(D,nk**D)
	integer :: ord(nk**D)
contains
	function compute_(node) result(rt)
		integer :: node
		character(:), allocatable :: rt
		integer, save :: nx,cnd_=0,im
		integer :: i,j,ik,i1,i2
		logical :: is_
		real(wp) :: rho(n),norm,W,tmp
		complex(wp) :: A(2,2),B(2,2),ntmp(n),det(2,2),det0
		integer, allocatable :: c(:)
		select case(node)
		case(PMT)
			rt="PMT"; if(debug) return
		case(output)
			rt="output"; if(debug) return
			!if(abs(Tk-0.1_wp)<1e-6_wp) then
				!call self%io([out_Gw,out_Gt],[out_Gw,out_Gt])
			!endif
		case(lattice)
			rt="lattice"; if(debug) return
			!$OMP PARALLEL DO PRIVATE(tmp)
			do ik=1,nk**D
				do i=0,D-1
					k(i+1,ik)=abs((mod((ik-1)/nk**i,nk)-nk/2)*2._wp*pi/nk)
				enddo
				do i=1,D
					do j=i+1,D
						if(k(i,ik)<k(j,ik)) then
							tmp=k(i,ik)
							k(i,ik)=k(j,ik)
							k(j,ik)=tmp
						endif
					enddo
				enddo
			enddo
			!$OMP END PARALLEL DO
			ord=[1:size(k,2)]
			call qsort(k,ord)
			k(:,:)=k(:,ord)
			call collect(k,c)
			allocate(self%uk(0:D,size(c)-1))
			self%uk(1:D,:)=k(:,c(1:size(c)-1))
			self%uk(0,:)=c(2:)-c(1:size(c)-1)

			!do ik=1,nk**D
				!do i=0,D-1
					!k(i+1,ik)=(mod((ik-1)/nk**i,nk)-nk/2)*pi/nk
				!enddo
			!enddo
			!allocate(self%uk(0:D,nk**D))
			!self%uk(1:D,:)=k
			!self%uk(0,:)=1._wp
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
			!$OMP PARALLEL DO
			do i=1,n
				rho(i)=A_c(omega(i))
			enddo
			!$OMP END PARALLEL DO
			norm=1._wp/integrate(omega,A=rho)
			self%rhoc=rho*norm
			self%rhoc=rho

			!do i=1,n
				!if(abs(omega(i))<fcut.and.abs(omega(i))>merge(0._wp,gap,isnan(gap))) then
					!self%rhoc(i)=0.5_wp
				!else
					!self%rhoc(i)=0._wp
				!endif
			!enddo

			self%Gc0%im=-pi*self%rhoc
			call set_realpart(omega,self%Gc0)

			if(.not.isnan(scgap)) then
				do i=1,n
					rho(i)=F_c(omega(i))
				enddo
				self%Fc0%im=-pi*rho*norm
				call set_realpart(omega,self%Fc0)
			endif


			self%rhov=0._wp
			!$OMP PARALLEL DO
			do i=1,size(omega_ek)
				do ik=1,size(self%uk,2)
					self%rhov(i)=self%rhov(i)-self%uk(0,ik)*dek(self%uk(1:D,ik))*imag(1._wp/(omega_ek(i)-ek(self%uk(1:D,ik))+img*eta))/pi
				enddo
			enddo
			!$OMP END PARALLEL DO
			self%rhov=self%rhov/nk**D
		case(rhophi)
			rt="rhophi"; if(debug) return
			do i=1,n
				rho(i)=A_phi(omega(i))
			enddo
			!norm=1._wp/integrate(omega,A=rho)
			!self%rhoc=rho*norm
			self%rhophi=rho

			!do i=1,n
				!if(abs(omega(i))<fcut.and.abs(omega(i))>merge(0._wp,gap,isnan(gap))) then
					!self%rhoc(i)=0.5_wp
				!else
					!self%rhoc(i)=0._wp
				!endif
			!enddo

			self%Gphi0%im=-pi*self%rhophi
			call set_realpart(omega,self%Gphi0)
		case(Gf)
			rt="Gf"; if(debug) return
			!$OMP PARALLEL DO PRIVATE(det,det0,tmp)
			do i=1,n
				select case(size(ed))
				case(1)
					self%Gf(i,:,:)%im=self%tSEf(i,:,:)%im/((diag(omega(i)+self%lambda-ed)-self%SEf(i,:,:)%re)**2+(self%SEf(i,:,:)%im)**2)
					self%Gf(i,:,:)%re=(diag(omega(i)+self%lambda-ed)-self%SEf(i,:,:)%re)/((diag(omega(i)+self%lambda-ed)-self%SEf(i,:,:)%re)**2+(self%SEf(i,:,:)%im)**2)
				case(2)
					det(2,2)=(omega(i)+self%lambda-ed(1)-self%SEf(i,1,1))
					det(1,1)=(omega(i)+self%lambda-ed(2)-self%SEf(i,2,2))
					det(1,2)=self%SEf(i,1,2)
					det(2,1)=self%SEf(i,2,1)
					det0=det(1,1)*det(2,2)-det(1,2)*det(2,1)
					self%Gf(i,1:2,1:2)%re=real(det/det0,kind=8)
					!imaginary part of det0*
					tmp=det(1,1)%re*self%tSEf(i,1,1)%im+self%tSEf(i,2,2)%im*det(2,2)%re+det(1,2)%re*self%tSEf(i,2,1)%im+self%tSEf(i,1,2)%im*det(2,1)%re
					self%Gf(i,1,1)%im=(-self%tSEf(i,2,2)%im*det0%re+det(1,1)%re*tmp)/(det0*conjg(det0))
					self%Gf(i,2,2)%im=(-self%tSEf(i,1,1)%im*det0%re+det(2,2)%re*tmp)/(det0*conjg(det0))
					self%Gf(i,1,2)%im=(self%tSEf(i,1,2)%im*det0%re+det(1,2)%re*tmp)/(det0*conjg(det0))
					self%Gf(i,2,1)%im=(self%tSEf(i,2,1)%im*det0%re+det(2,1)%re*tmp)/(det0*conjg(det0))
				case default
					stop 'ed larger than 2 is not implimented yet'
				end select
			enddo
			!$OMP END PARALLEL DO
		case(Gb)
			rt="Gb"; if(debug) return
			!$OMP PARALLEL DO
			do i=1,n
				self%Gb(i)%im=self%tSEb(i)%im/((omega(i)+self%lambda-self%SEb(i)%re)**2+(self%SEb(i)%im)**2)
				self%Gb(i)%re=(omega(i)+self%lambda-self%SEb(i)%re)/((omega(i)+self%lambda-self%SEb(i)%re)**2+(self%SEb(i)%im)**2)
			enddo
			!$OMP END PARALLEL DO
		case(Ga)
			rt="Ga"; if(debug) return
			if(.not.isnan(U)) then
				if(ch>1) then
					write(*,*)"multichannel for finite U is not considered"
					stop
				endif
				if(size(ed)>1) then
					write(*,*)"ed spliting for finite U is not considered"
					stop
				endif
				!$OMP PARALLEL DO
				do i=1,n
					self%Ga(i)%im=self%tSEa(i)%im/((omega(i)+self%lambda-2._wp*ed(1)-U-self%SEa(i)%re)**2+(self%SEa(i)%im)**2)
					self%Ga(i)%re=(omega(i)+self%lambda-2._wp*ed(1)-U-self%SEa(i)%re)/((omega(i)+self%lambda-2._wp*ed(1)-U-self%SEa(i)%re)**2+(self%SEa(i)%im)**2)
				enddo
				!$OMP END PARALLEL DO
			endif
		case(Q)
			rt="Q"; if(debug) return
			do i=1,n
				ntmp(i)=sum(diag(self%Gf(i,1:2,1:2)))
			enddo
			call get_Gtau(omega,G=merge((2._wp*ntmp+ch*self%Gb),(2._wp*ntmp+ch*self%Gb+self%Ga),isnan(U)),eta=-1,time=[-eps*1e3_wp],rt=self%Q)
		case(lambda)
			rt="lambda"; if(debug) return
			!write(*,*)0.5_wp*(self%SEf(n/2)%re+self%SEf(n/2+1)%re)+ed-0.5_wp*(self%SEb(n/2)%re+self%SEb(n/2+1)%re)
		case(tlambda)
			rt="tlambda"; if(debug) return
			self%lambda=min(minval(0.5_wp*(diag(self%SEf(n/2,:,:)%re+self%SEf(n/2+1,:,:)%re))+ed),(0.5_wp*(self%SEb(n/2)%re+self%SEb(n/2+1)%re)))
		case(tSEf)
			rt="tSEf"; if(debug) return
			self%tSEf=0._wp
			ntmp=0._wp
			call get_convolution(omega,A=self%Gc0([n:1:-1]),B=self%Gb,rt=ntmp)
			do i=1,n
				self%tSEf(i,:,:)=self%tSEf(i,:,:)+ch*V2*ntmp(i)
			enddo
			if(.not.isnan(U)) then
				if(size(ed)>1) then
					write(*,*)"ed spliting for finite U is not considered"
					stop
				endif
				call get_convolution(omega,A=V2(1,1)*self%Gc0,B=self%Ga,rt=self%tSEf(1:n,1,1))
			endif
			if(.not.isnan(g2)) then
				do i1=1,size(ed)
					do i2=1,size(ed)
						call get_convolution(omega,A=-g2*self%Gphi0([n:1:-1]),B=self%Gf(1:n,i1,i2),rt=self%tSEf(1:n,i1,i2),Bbath=underscore)
					enddo
				enddo
			endif
			if(ocd) then
				if(size(ed)>1) then
					write(*,*)"ed spliting for OCD is not considered"
					stop
				endif
				if(.not.isnan(scgap)) then
					call get_convolution2(omega,ga=2._wp*V2(1,1)**2*self%Fc0(n:1:-1),gb=self%Fc0(n:1:-1),A=self%Gb,B=self%Gb,C=self%Gf(1:n,1,1),rt=self%tSEf(1:n,1,1))
				endif
				if(.not.isnan(U)) then
					call get_convolution2(omega,ga=-2._wp*V2(1,1)**2*self%Gc0(n:1:-1),gb=self%Gc0,A=self%Gb,B=self%Ga,C=self%Gf(1:n,1,1),rt=self%tSEf(1:n,1,1))
					if(.not.isnan(scgap)) then
						call get_convolution2(omega,ga=2._wp*V2(1,1)**2*self%Fc0(n:1:-1),gb=self%Fc0(n:1:-1),A=self%Ga,B=self%Ga,C=self%Gf(1:n,1,1),rt=self%tSEf(1:n,1,1))
					endif
				endif
			endif
		case(tSEb)
			rt="tSEb"; if(debug) return
			self%tSEb=0._wp
			do i=1,n
				ntmp(i)=sum(self%Gf(i,:,:)*V2)
			enddo
			call get_convolution(omega,A=2._wp*self%Gc0,B=ntmp,rt=self%tSEb)
			if(ocd) then
				if(size(ed)>1) then
					write(*,*)"ed spliting for OCD is not considered"
					stop
				endif
				if(.not.isnan(scgap)) then
					call get_convolution2(omega,ga=2._wp*V2(1,1)**2*self%Fc0,gb=self%Fc0(n:1:-1),A=self%Gf(1:n,1,1),B=self%Gf(1:n,1,1),C=self%Gb,rt=self%tSEb)
				endif
				if(.not.isnan(U)) then
					call get_convolution2(omega,ga=-2._wp*V2(1,1)**2*self%Gc0,gb=self%Gc0,A=self%Gf(1:n,1,1),B=self%Gf(1:n,1,1),C=self%Ga,rt=self%tSEb)
				endif
			endif
		case(tSEa)
			rt="tSEa"; if(debug) return
			if(.not.isnan(U)) then
				if(ch>1) then
					write(*,*)"multichannel for finite U is not considered"
					stop
				endif
				if(size(ed)>1) then
					write(*,*)"ed spliting for finite U is not considered"
					stop
				endif
				self%tSEa=0._wp
				call get_convolution(omega,A=2._wp*V2(1,1)*self%Gc0([n:1:-1]),B=self%Gf(1:n,1,1),rt=self%tSEa)
				if(.not.isnan(g2)) then
					call get_convolution(omega,A=-2._wp*g2*self%Gphi0([n:1:-1]),B=self%Ga,rt=self%tSEa,Bbath=underscore)
				endif
				if(ocd) then
					if(.not.isnan(scgap)) then
						call get_convolution2(omega,ga=2._wp*V2(1,1)**2*self%Fc0,gb=self%Fc0(n:1:-1),A=self%Gf(1:n,1,1),B=self%Gf(1:n,1,1),C=self%Ga,rt=self%tSEa)
					endif
					!diagram 3
					call get_convolution2(omega,ga=-2._wp*V2(1,1)**2*self%Gc0([n:1:-1]),gb=self%Gc0([n:1:-1]),A=self%Gf(1:n,1,1),B=self%Gf(1:n,1,1),C=self%Gb,rt=self%tSEa)
				endif
			endif
		case(SEf)
			rt="SEf"; if(debug) return
			do i=1,n
				self%SEf(i,:,:)%im=self%tSEf(i,:,:)%im*ff(-beta*omega(i))
			enddo
			do i1=1,size(ed)
				do i2=1,size(ed)
					call set_realpart(omega,self%SEf(1:n,i1,i2))
				enddo
			enddo
		case(SEb)
			rt="SEb"; if(debug) return
			self%SEb%im=self%tSEb%im*ff(-beta*omega)
			call set_realpart(omega,self%SEb)
		case(SEa)
			rt="SEa"; if(debug) return
			if(.not.isnan(U)) then
				if(ch>1) then
					write(*,*)"multichannel for finite U is not considered"
					stop
				endif
				if(size(ed)>1) then
					write(*,*)"ed spliting for finite U is not considered"
					stop
				endif
				self%SEa%im=self%tSEa%im*ff(-beta*omega)
				call set_realpart(omega,self%SEa)
			endif
		case(Gd)
			rt="Gd"; if(debug) return
			self%Gd=0._wp
			self%Fd=0._wp
			do i=1,n
				ntmp(i)=sum(self%Gf(i,:,:)*V2)
				!ntmp(i)=sum(self%Gf(i,:,:))
			enddo
			call get_convolution(omega,A=self%Gb/self%Q(1),B=ntmp,rt=self%Gd,two=underscore)
			if(.not.isnan(U)) then
				if(size(ed)>1) then
					write(*,*)"ed spliting for finite U is not considered"
					stop
				endif
				call get_convolution(omega,A=V2(1,1)/self%Q(1)*self%Gf(1:n,1,1),B=self%Ga,rt=self%Gd,two=underscore)
			endif
			if(ocd) then
				!call get_convolution2p(omega,ga=-2._wp*V2(1,1)**2/self%Q(1)*self%Gc0,gb=self%Gb,A=self%Gf(1:n,1,1),B=self%Gf(1:n,1,1),C=self%Ga,two=underscore,rt=self%Gd)
				if(.not.isnan(U)) then
					call get_convolution2pp(omega,ga=-2._wp*V2(1,1)**2/self%Q(1)*self%Gb,gb=self%Gc0,A=self%Gf(1:n,1,1),B=self%Gf(1:n,1,1),C=self%Ga,two=underscore,rt=self%Gd)
				endif
				if(.not.isnan(scgap)) then
					call get_convolution2pp(omega,ga=2._wp*V2(1,1)**2/self%Q(1)*self%Gf(1:n,1,1),gb=self%Fc0,A=self%Gb,B=self%Gb,C=self%Gf(1:n,1,1),two=underscore,rt=self%Fd)
					if(.not.isnan(U)) then
						call get_convolution2pp(omega,ga=2._wp*V2(1,1)**2/self%Q(1)*self%Gf(1:n,1,1),gb=self%Fc0,A=self%Ga,B=self%Ga,C=self%Gf(1:n,1,1),two=underscore,rt=self%Fd)
					endif
				endif
			endif
			call set_realpart(omega,self%Gd)
			call set_realpart(omega,self%Fd)
		case(Fimp)
			rt="Fimp"; if(debug) return
			self%Fimp=self%lambda-Tk*log(self%Q(1))
		case(optcond)
			rt="optcond"; if(debug) return
			call get_opt_conductivity(omega,self%Gd,self%uk,rt=self%sigma)
		case(SEc)
			rt="SEc"; if(debug) return
			self%SEc=1._wp/(1._wp/(self%Gd)+self%Gc0)
		case(Gc)
			rt="Gc"; if(debug) return
			self%Gc=0._wp
			self%Fc=0._wp
			!!self%Gc=1._wp/(1._wp/self%Gc0-self%SEc)

			!self%Gc=self%Gc0+self%Gc0*(self%Gd)*self%Gc0
			!self%Gc%im=self%Gc0%im-self%Gc0%im*(self%Gd%im)*self%Gc0%im
			!if(.not.isnan(scgap)) then
				!self%Gc=self%Gc+self%Fc0*conjg(-self%Gd(n:1:-1))*self%Fc0
			!endif

			!self%Gc=(-conjg(self%Gc0(n:1:-1)))+self%Fc0*(self%Gd)*self%Fc0
			!if(.not.isnan(scgap)) then
				!self%Gc=self%Gc+(-conjg(self%Gc0(n:1:-1)))*conjg(-self%Gd(n:1:-1))*(-conjg(self%Gc0(n:1:-1)))
			!endif
			!self%Gc=(-conjg(self%Gc(n:1:-1)))

			!do i=1,n
				!A(1,1)=omega(i)-ed-V2*self%Gc0(i)
				!A(2,2)=omega(i)+ed-V2*(-conjg(self%Gc0(n+1-i)))
				!A(1,2)=-V2*self%Fc0(i)
				!A(2,1)=-V2*self%Fc0(i)
				!call mat_inv(A)
				!self%Gd(i)=A(1,1)
				!self%Fd(i)=A(1,2)
			!enddo

			do i=1,n
				B=0._wp
				B(1,1)=self%Gc0(i)
				B(2,2)=(-conjg(self%Gc0(n+1-i)))

				if(.not.isnan(scgap)) then
					B(1,2)=self%Fc0(i)
					B(2,1)=self%Fc0(i)
				endif
				A=0._wp
				A(1,1)=self%Gd(i)
				A(2,2)=(-conjg(self%Gd(n+1-i)))
				A(1,2)=-(self%Fd(i))
				A(2,1)=-(self%Fd(i))

				!A(1,1)=1._wp
				!A(2,2)=1._wp
				!call mat_inv(A)
				!A=A-B
				!call mat_inv(A)

				self%Gc(i,1)=self%Gc(i,1)+B(1,1)+sum(matmul(B(1,:),A)*B(:,1))
				self%Gc(n-i+1,2)=self%Gc(n-i+1,2)-conjg(B(2,2)+sum(matmul(B(2,:),A)*B(:,2)))
				self%Fc(i)=self%Fc(i)+B(1,2)+sum(matmul(B(1,:),A)*B(:,2))
			enddo
			!write(*,*)integrate(omega,A=-self%Gc0%im/pi),integrate(omega,A=-(self%Gc(:,1)%im+self%Gc(:,2)%im)/pi),integrate(omega,A=-self%Gd%im/pi)
		case(Gft)
			rt="Gft"; if(debug) return
			!call get_Gtau(omega,G=self%Gf,eta=-1,time=tau,rt=self%Gft)
		case(Gbt)
			rt="GBt"; if(debug) return
			!call get_Gtau(omega,G=self%Gb,eta=1,time=tau,rt=self%GBt)
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
				write(*,"('SC=',i2,' T=',es10.3,', V2: ',es9.2,', ed: ',es9.2,', ch: ',i3,', r: ',es9.2,' initial...')")node,Tk,V2,ed,ch,r
				if(onode(cnd_,0)==1.and.onode(cnd_,1)==lambda) then
					underscore=next(ytol=iter_err(cnd_),dx=iter_rate(cnd_),id=abs(cnd_))
				else
					call mbroyden(0,x,x_,iter_rate(node),nsave=80,id=abs(node))
				endif
			endif
			write(*,"('SC=',i2,' T=',es10.3,' err=',es10.3,' iter=',i4,' lambda=',es10.3,' Q=',es10.3)")node,Tk,iter_err(abs(node)),iter_exit(node),self%lambda,self%Q(1)
			if(iter_exit(node)==abs(niter(node)).or.(iter_err(abs(node))<=iter_err(node))) then
				call mbroyden(huge(1),x,x_,id=abs(node))
				deallocate(x,x_)
				write(*,"('SC=',i2,' T=',es10.3,', V2: ',es9.2,', ed: ',es9.2,', ch: ',i3,', r: ',es9.2,' finish...')")node,Tk,V2,ed,ch,r
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
				!iter_err(abs(cnd_))=maxval(abs(x_)/merge(1e-4_wp,abs(x),abs(x)<1e-4_wp))
				if(onode(cnd_,0)==1.and.onode(cnd_,1)==lambda) then
					underscore=next(self%lambda,x_(1),id=abs(cnd_))
					return
				else
					if(niter(cnd_)>0) then
						call mbroyden(iter_exit(cnd_),x,x_,id=abs(cnd_))
						!x=x_+x
					else
						x=x+iter_rate(cnd_)*x_
					endif
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
							x(nx+1:nx+size(self%tSEf))=reshape(self%tSEf%im,[n*size(ed)**2])
						case(xx_diff)
							x_(nx+1:nx+size(self%tSEf))=(reshape(self%tSEf%im,[n*size(ed)**2])-x(nx+1:nx+size(self%tSEf)))
						case(xx_new)
							!self%tSEf%im=-abs(x(nx+1:nx+size(self%tSEf)))
							self%tSEf%im=reshape(x(nx+1:nx+size(self%tSEf)),shape(self%tSEf))
							!self%tSEf%im=x(nx+1:nx+size(self%tSEf))+0.5_wp*x_(nx+1:nx+size(self%tSEf))
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
							!self%tSEb%im=-abs(x(nx+1:nx+size(self%tSEb)))
							self%tSEb%im=x(nx+1:nx+size(self%tSEb))
						end select
					endif
					nx=nx+size(self%tSEb)
				case(tSEa)
					! SEb self-consistent
					if(allocated(x)) then
						select case(flag)
						case(xx_save)
							x(nx+1:nx+size(self%tSEa))=self%tSEa%im
						case(xx_diff)
							x_(nx+1:nx+size(self%tSEa))=self%tSEa%im-x(nx+1:nx+size(self%tSEa))
						case(xx_new)
							self%tSEa%im=x(nx+1:nx+size(self%tSEa))
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
							!x_(nx+1)=(self%Q(1)-1._wp)*0.01_wp
							!x_(nx+1)=(minval([0.5_wp*(self%SEf(n/2)%re+self%SEf(n/2+1)%re)+ed,(0.5_wp*(self%SEb(n/2)%re+self%SEb(n/2+1)%re)),merge(1e9_wp,(0.5_wp*(self%SEa(n/2)%re+self%SEa(n/2+1)%re)),isnan(U))])-x(nx+1))
							!x_(nx+1)=0.5_wp*(self%SEb(n/2)%re+self%SEb(n/2+1)%re)-x(nx+1)
							x_(nx+1)=minval(0.5_wp*(diag(self%SEf(n/2,:,:)%re+self%SEf(n/2+1,:,:)%re))+ed)-x(nx+1)
							!x_(nx+1)=0.5_wp*sum(0.5_wp*(diag(self%SEf(n/2,:,:)%re+self%SEf(n/2+1,:,:)%re))+ed)-x(nx+1)
							!do i=n/2,n-1
								!if((omega(i)+0.5_wp*(self%SEf(n/2)%re+self%SEf(n/2+1)%re)+ed-self%SEb(i)%re)*(omega(i+1)+0.5_wp*(self%SEf(n/2)%re+self%SEf(n/2+1)%re)+ed-self%SEb(i+1)%re)<0._wp) then
									!write(*,*)omega(i)
									!exit
								!endif
							!enddo
						case(xx_new)
							!if(abs(self%lambda-x(nx+1))>max(0.01_wp*Tk,0._wp)) then
								!self%lambda=self%lambda-max(0.01_wp*Tk,0._wp)*sign(1._wp,(self%Q(1)-1._wp))
							!else
								self%lambda=x(nx+1)
							!endif
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
				write(ut(k),"(A)")"T omega omegaT ff rGf iGf rGb iGb rGa iGa rSEf iSEf rSEb iSEb rSEa iSEa rGc0 iGc0 rFc0 iFc0 rSEc iSEc rGcu iGcu rGcd iGcd rFc iFc rGd iGd rFd iFd rsigma isigma rGphi0 iGphi0 omega_ek rhov"
				do i=1,max(n,size(omega_ek))
					if(i>n) then
						write(*,*)"io error size of omega_ek is larger than n"
						stop
					endif
					write(ut(k),"(e29.20e3$)")Tk,omega(i),omega(i)/Tk,ff(-beta*omega(i)),sum(self%Gf(i,:,:)),self%Gb(i),self%Ga(i),sum(self%SEf(i,:,:)),self%SEb(i),self%SEa(i),self%Gc0(i),self%Fc0(i),self%SEc(i),self%Gc(i,1),self%Gc(i,2),self%Fc(i),self%Gd(i),self%Fd(i),self%sigma(i),self%Gphi0(i)
					if(i<=size(omega_ek)) then
						write(ut(k),"(e29.20e3$)")omega_ek(i),self%rhov(i)
					endif
					write(ut(k),"(x)")
				enddo
				write(ut(k),"(x)")
			case(out_T)
				write(ut(k),"(*(e29.20e3))")Tk,Tkunit,self%lambda,self%Q(1),self%Fimp,minval(diag((self%Gf(n/2+1,:,:)%re+self%Gf(n/2,:,:)%re)/2._wp)),(self%GB(n/2+1)+self%GB(n/2))/2.,(self%SEc(n/2+1)+self%SEc(n/2))/2.
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
					write(ut(k),"(*(e29.20e3))")omega(i),self%SEf(i,1,1)%im
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
							self%SEf(i,1,1)%im=dat(1)
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
				call set_realpart(omega,self%SEf(1:n,1,1))
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
							write(i,"(A)")"T unit lambda Q Fimp rGf0 iGf0 rGB0 iGB0 rSEc iSEc"
						endif
					endif
				enddo
			case(out_param)
				write(*,"('V2: ',es9.2,', ed: ',es9.2,', ch: ',i3,', r: ',es9.2)")V2,ed,ch,r
				do i=10,1000
					if(i==in_init) cycle
					inquire(unit=i, opened=isopen)
					if(isopen) then
						write(i,"('#V2: ',es9.2,', ed: ',es9.2,', ch: ',i3,', r: ',es9.2)")V2,ed,ch,r
					endif
				enddo
			end select
		end do
	end subroutine
	real(wp) elemental function A_c(x)
		real(wp), intent(in) :: x
		integer :: ik
		!A_c=0._wp
		!do ik=1,size(self%uk,2)
			!A_c=A_c-self%uk(0,ik)*imag(1._wp/(x-ek(self%uk(1:D,ik))+img*eta))/pi
		!enddo
		!A_c=A_c/nk**D
		!return
		A_c=exp(-x**2/pi)
		return
		
		if(abs(x)<fcut) then
			if(isnan(scgap)) then
				A_c=(abs(x)/fcut)**r
			else
				if(abs(x)>scgap) then
					A_c=abs(x)/sqrt(x**2-scgap**2)
				else
					A_c=0._wp
				endif
				!A_c=imag((img*0.00001_wp+x)/sqrt(scgap**2+(0.00001_wp-img*x)**2))
			endif
		elseif(abs(x)<=(-omega(1))) then
			if(isnan(scgap)) then
				A_c=exp(-(x**2-fcut**2)/pi*fbeta)
			else
				A_c=abs(fcut)/sqrt(fcut**2-scgap**2)*exp(-(x**2-fcut**2)/pi*fbeta)
			endif
		else
			A_c=0._wp
		endif
		if(.not.isnan(gap)) then
			if(abs(x)<gap) then
				A_c=0._wp
				!A_c=1e-1_wp
			endif
		endif
	end function
	real(wp) elemental function F_c(x)
		real(wp), intent(in) :: x
		if(abs(x)<fcut) then
			if(abs(x)>scgap) then
				F_c=-sign(1._wp,x)*scgap/sqrt(x**2-scgap**2)
			else
				F_c=0._wp
			endif
		elseif(abs(x)<=(-omega(1))) then
			F_c=-sign(1._wp,x)*scgap/sqrt(fcut**2-scgap**2)*exp(-(x**2-fcut**2)/pi*fbeta)
		else
			F_c=0._wp
		endif
	end function
	real(wp) elemental function A_phi(x)
		real(wp), intent(in) :: x
		A_phi=-imag(1._wp/(x-boson+img*etab)-1._wp/(x+boson+img*etab))/pi
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
	subroutine get_convolution2p(omega,ga,gb,A,B,C,two,rt)
		! calculate ga(x)gb(y)A(omg+x)B(omg+y)C(omg+x+y)
		real(wp) :: omega(:)
		complex(wp) :: ga(:),gb(:),A(:),B(:),C(:),rt(:)
		logical, optional :: two
		integer :: i,j,k,i1,i2,j1,j2,j3,iwb(0:n),iwmb(0:n),iwb1(0:n),iwb2(0:n),d(2),sg(3)
		complex(wp) :: A_,B_,C_,ga_,gb_
		real(wp) :: omg,ei1,ei2,ej1,ej2,rt_,rt1,rt2,rt1_,rt2_
		!$OMP PARALLEL DO PRIVATE(omg,iwb,iwmb,iwb1,iwb2,ga_,gb_,A_,B_,C_,ei1,ei2,ej1,ej2,d,rt1,rt2,rt_,rt1_,rt2_,j1,j2,j3,i1,i2)
		do k=1,n
			omg=omega(k)
			do i=1,n
				iwb(i)=find(omega,(omg+omega(i)))
			enddo
			iwb(0)=iwb(1)
			do i=1,n
				iwmb(i)=find(omega,(omega(i)-omg))
			enddo
			iwmb(0)=iwmb(1)
			ei1=omega(1)
			rt1=0._wp
			rt_=0._wp
			! out integration
			do i=1,n
				i1=iwmb(i-1)+1
				i2=iwb(i-1)+1
				do
					d(1)=0
					if(i1/=iwmb(i)+1) then
						ei2=(omega(i1)+omg)
						d(1)=1
					endif
					if(i2/=iwb(i)+1) then
						if(d(1)==0.or.(d(1)/=0.and.ei2>(omega(i2)-omg))) then
							ei2=(omega(i2)-omg)
							d(1)=2
						endif
					endif
					if(d(1)==0) then
						ei2=omega(i)
						ga_=ga(i)
					else
						ga_=ga(i-1)+(ei2-omega(i-1))/(omega(i)-omega(i-1))*(ga(i)-ga(i-1))
					endif
					do j=1,n
						iwb1(j)=find(omega,(ei2+omega(j)))
					enddo
					iwb1(0)=iwb1(1)
					do j=1,n
						iwb2(j)=find(omega,omg+ei2+omega(j))
					enddo
					iwb2(0)=iwb2(1)
					ej1=omega(1)
					rt1_=0._wp
					rt2=0._wp
					! inner integration
					do j=1,n
						j1=iwb1(j-1)+1
						j2=iwb(j-1)+1
						j3=iwb2(j-1)+1
						do
							d(2)=0
							if(j1/=iwb1(j)+1) then
								ej2=(omega(j1)-ei2)
								!ej2=(omega(j1)-omg)
								d(2)=1
							endif
							if(j2/=iwb(j)+1) then
								if(d(2)==0.or.(d(2)/=0.and.ej2>(omega(j2)-omg))) then
									ej2=(omega(j2)-omg)
									d(2)=2
								endif
							endif
							if(j3/=iwb2(j)+1) then
								if(d(2)==0.or.(d(2)/=0.and.ej2>(omega(j3)-omg-ei2))) then
									ej2=(omega(j3)-omg-ei2)
									d(2)=3
								endif
							endif
							if(d(2)==0) then
								ej2=omega(j)
								gb_=gb(j)
							else
								gb_=gb(j-1)+(ej2-omega(j-1))/(omega(j)-omega(j-1))*(gb(j)-gb(j-1))
							endif
							if(j1>1.and.j1<=n) then
								A_=A(j1-1)+(ei2+ej2-omega(j1-1))/(omega(j1)-omega(j1-1))*(A(j1)-A(j1-1))
							else
								A_=0._wp
							endif
							if(j2>1.and.j2<=n) then
								B_=B(j2-1)+(omg+ej2-omega(j2-1))/(omega(j2)-omega(j2-1))*(B(j2)-B(j2-1))
							else
								B_=0._wp
							endif
							if(j3>1.and.j3<=n) then
								C_=C(j3-1)+(omg+ei2+ej2-omega(j3-1))/(omega(j3)-omega(j3-1))*(C(j3)-C(j3-1))
							else
								C_=0._wp
							endif
							rt2_=sum(f(beta*ei2,beta*ej2,beta*omg)*ga_%im*[gb_%im*A_%re*B_%im*C_%re,gb_%re*A_%im*B_%im*C_%re,gb_%im*A_%re*B_%re*C_%im,gb_%re*A_%im*B_%re*C_%im])
							rt2=rt2+0.5_wp*(ej2-ej1)*(rt1_+rt2_)
							ej1=ej2
							rt1_=rt2_
							select case(d(2))
							case(1)
								j1=j1+1
							case(2)
								j2=j2+1
							case(3)
								j3=j3+1
							case(0)
								exit
							end select
						enddo
					enddo
					rt_=rt_+0.5_wp*(ei2-ei1)*(rt1+rt2)
					ei1=ei2
					rt1=rt2
					select case(d(1))
					case(1)
						i1=i1+1
					case(2)
						i2=i2+1
					case(0)
						exit
					end select
				enddo
			enddo
			rt(k)%im=rt(k)%im-rt_/pi**2
		enddo
		!$OMP END PARALLEL DO
	contains
		function f(x,y,z)
			real(wp) :: f(4)
			real(wp) :: x,y,z
			!IRIR
			f(1)=1._wp/((1._wp+(exp(sign(1._wp,z)*y)+exp((-sign(1._wp,z)*(y+z))))/(1._wp+exp(-abs(z))))*(1._wp+exp(x)))
			!RIIR
			f(2)=1._wp/((1._wp+(exp(sign(1._wp,z)*y)+exp((-sign(1._wp,z)*(y+z))))/(1._wp+exp(-abs(z))))*(1._wp+(exp(sign(1._wp,y)*x)+exp((-sign(1._wp,y)*(x+y))))/(1._wp+exp(-abs(y)))))
			!IRRI
			f(3)=1._wp/((1._wp+exp(sign(1._wp,z)*x))*(1._wp+exp(-sign(1._wp,z)*(z+x+y)))*(1._wp+exp(sign(1._wp,z)*y))/(1._wp+exp(-abs(z))))
			!RIRI
			f(4)=1._wp/((1._wp+exp(sign(1._wp,sign(1._wp,z)+sign(1._wp,y)-1._wp)*x))*(1._wp+exp(sign(1._wp,sign(1._wp,z)-sign(1._wp,y)+1._wp)*y))*(1._wp+exp(-sign(1._wp,sign(1._wp,y)-sign(1._wp,z)+1._wp)*(x+y)))*(1._wp+exp(-sign(1._wp,z)*(x+y+z)))/((1._wp+exp(-abs(z)))*(1._wp+exp(-abs(y)))))
		end function
	end subroutine
	subroutine get_convolution2pp(omega,ga,gb,A,B,C,two,rt)
		! calculate ga(x)gb(y)A(omg+x)B(x+y)C(omg+x+y)
		real(wp) :: omega(:)
		complex(wp) :: ga(:),gb(:),A(:),B(:),C(:),rt(:)
		logical, optional :: two
		integer :: i,j,k,j1,j2,i1,iwb(0:n),iwb1(0:n),iwb2(0:n),d(2)
		complex(wp) :: A_,B_,C_,rt1,rt2,ga_,gb_
		real(wp) :: omg,ei1,ei2,ej1,ej2,rt_,rt1_,rt2_
		!$OMP PARALLEL DO PRIVATE(omg,iwb,iwb1,iwb2,ga_,gb_,A_,B_,C_,ei1,ei2,ej1,ej2,d,rt1,rt2,rt_,rt1_,rt2_,i1,j1,j2)
		do k=1,n
			omg=omega(k)
			do i=1,n
				iwb(i)=find(omega,(omg+omega(i)))
			enddo
			iwb(0)=iwb(1)
			ei1=omega(1)
			rt1=0._wp
			rt_=0._wp
			! out integration
			do i=1,n
				i1=iwb(i-1)+1
				do
					d(1)=0
					if(i1/=iwb(i)+1) then
						ei2=(omega(i1)-omg)
						d(1)=1
					endif
					if(d(1)==0) then
						ei2=omega(i)
						ga_=ga(i)
					else
						ga_=ga(i-1)+(ei2-omega(i-1))/(omega(i)-omega(i-1))*(ga(i)-ga(i-1))
					endif
					if(i1>1.and.i1<=n) then
						A_=A(i1-1)+((omg+ei2)-omega(i1-1))/(omega(i1)-omega(i1-1))*(A(i1)-A(i1-1))
					else
						A_=0._wp
					endif

					do j=1,n
						iwb1(j)=find(omega,(ei2+omega(j)))
						!iwb1(j)=find(omega,(omg+omega(j)))
					enddo
					iwb1(0)=iwb1(1)
					do j=1,n
						iwb2(j)=find(omega,omg+ei2+omega(j))
					enddo
					iwb2(0)=iwb2(1)
					ej1=omega(1)
					rt2=0._wp
					rt1_=0._wp
					! inner integration
					do j=1,n
						j1=iwb1(j-1)+1
						j2=iwb2(j-1)+1
						do
							d(2)=0
							if(j1/=iwb1(j)+1) then
								ej2=(omega(j1)-ei2)
								d(2)=1
							endif
							if(j2/=iwb2(j)+1) then
								if(d(2)==0.or.(d(2)/=0.and.ej2>(omega(j2)-omg-ei2))) then
									ej2=(omega(j2)-omg-ei2)
									d(2)=2
								endif
							endif
							if(d(2)==0) then
								ej2=omega(j)
								gb_=gb(j)
							else
								gb_=gb(j-1)+(ej2-omega(j-1))/(omega(j)-omega(j-1))*(gb(j)-gb(j-1))
							endif
							if(j1>1.and.j1<=n) then
								B_=B(j1-1)+(ei2+ej2-omega(j1-1))/(omega(j1)-omega(j1-1))*(B(j1)-B(j1-1))
							else
								B_=0._wp
							endif
							if(j2>1.and.j2<=n) then
								C_=C(j2-1)+(omg+ei2+ej2-omega(j2-1))/(omega(j2)-omega(j2-1))*(C(j2)-C(j2-1))
							else
								C_=0._wp
							endif
							rt2_=sum(f(beta*ej2,beta*ei2,beta*omg)*gb_%im*[ga_%im*B_%re*A_%im*C_%re,ga_%re*B_%im*A_%im*C_%re,ga_%im*B_%re*A_%re*C_%im,ga_%re*B_%im*A_%re*C_%im])
							rt2=rt2+0.5_wp*(ej2-ej1)*(rt1_+rt2_)
							ej1=ej2
							rt1_=rt2_
							select case(d(2))
							case(1)
								j1=j1+1
							case(2)
								j2=j2+1
							case(0)
								exit
							end select
						enddo
					enddo
					!rt2=rt2
					rt_=rt_+0.5_wp*(ei2-ei1)*(rt1+rt2)
					ei1=ei2
					rt1=rt2
					select case(d(1))
					case(1)
						i1=i1+1
					case(0)
						exit
					end select
				enddo
			enddo
			rt(k)%im=rt(k)%im-rt_/pi**2
		enddo
		!$OMP END PARALLEL DO
	contains
		function f(x,y,z)
			real(wp) :: f(4)
			real(wp) :: x,y,z
			!real(wp) :: sg,esx,esy,oez,esz,ex,ey,ez
			!ex=exp(x)
			!ey=exp(y)
			!ez=exp(-z)
			!if(z>=0._wp) then
				!esx=ex
				!esy=ey
				!esz=ez
			!else
				!esx=1._wp/ex
				!esy=1._wp/ey
				!esz=1._wp/ez
			!endif
			!oez=1._wp/(1._wp+esz)
			!if(y>=0._wp) then
				!esx=ex
				!esy=ey
				!esz=ez
			!else
				!esx=1._wp/ex
				!esy=1._wp/ey
				!esz=1._wp/ez
			!endif
			!!IRIR
			!f(1)=1._wp/((1._wp+(esy+esz/esy)/oez)*(1._wp+ex))
			!!RIIR
			!f(2)=1._wp/((1._wp+(esy+esz/esy)/oez)*(1._wp+(exp(sign(1._wp,y)*x)+exp((-sign(1._wp,y)*(x+y))))/(1._wp+exp(-abs(y)))))
			!!IRRI
			!f(3)=1._wp/((1._wp+exp(sign(1._wp,z)*x))*(1._wp+exp(-sign(1._wp,z)*(z+x+y)))*(1._wp+exp(sign(1._wp,z)*y))/(1._wp+exp(-abs(z))))
			!!RIRI
			!f(4)=1._wp/((1._wp+exp(sign(1._wp,sign(1._wp,z)+sign(1._wp,y)-1._wp)*x))*(1._wp+exp(sign(1._wp,sign(1._wp,z)-sign(1._wp,y)+1._wp)*y))*(1._wp+exp(-sign(1._wp,sign(1._wp,y)-sign(1._wp,z)+1._wp)*(x+y)))*(1._wp+exp(-sign(1._wp,z)*(x+y+z)))/((1._wp+exp(-abs(z)))*(1._wp+exp(-abs(y)))))

			!IRIR
			f(1)=1._wp/((1._wp+(exp(sign(1._wp,z)*y)+exp((-sign(1._wp,z)*(y+z))))/(1._wp+exp(-abs(z))))*(1._wp+exp(x)))
			!RIIR
			f(2)=1._wp/((1._wp+(exp(sign(1._wp,z)*y)+exp((-sign(1._wp,z)*(y+z))))/(1._wp+exp(-abs(z))))*(1._wp+(exp(sign(1._wp,y)*x)+exp((-sign(1._wp,y)*(x+y))))/(1._wp+exp(-abs(y)))))
			!IRRI
			f(3)=1._wp/((1._wp+exp(sign(1._wp,z)*x))*(1._wp+exp(-sign(1._wp,z)*(z+x+y)))*(1._wp+exp(sign(1._wp,z)*y))/(1._wp+exp(-abs(z))))
			!RIRI
			f(4)=1._wp/((1._wp+exp(sign(1._wp,sign(1._wp,z)+sign(1._wp,y)-1._wp)*x))*(1._wp+exp(sign(1._wp,sign(1._wp,z)-sign(1._wp,y)+1._wp)*y))*(1._wp+exp(-sign(1._wp,sign(1._wp,y)-sign(1._wp,z)+1._wp)*(x+y)))*(1._wp+exp(-sign(1._wp,z)*(x+y+z)))/((1._wp+exp(-abs(z)))*(1._wp+exp(-abs(y)))))
		end function
	end subroutine
	subroutine get_convolution2(omega,ga,gb,A,B,C,two,rt)
		! calculate ga(x)gb(y)A(omg+x)B(omg+y)C(omg+x+y)
		real(wp) :: omega(:)
		complex(wp) :: ga(:),gb(:),A(:),B(:),C(:),rt(:)
		logical, optional :: two
		integer :: i,j,k,j1,j2,i1,iwb(0:n),iwb2(0:n),d(2)
		complex(wp) :: A_,B_,C_,ga_,gb_
		real(wp) :: omg,ei1,ei2,ej1,ej2,rt_,rt1_,rt2_,rt1,rt2
		!write(*,*)"start"
		!call start_time(id=1)
		!$OMP PARALLEL DO PRIVATE(omg,iwb,iwb2,ga_,gb_,A_,B_,C_,ei1,ei2,ej1,ej2,d,rt1,rt2,rt_,rt1_,rt2_,i1,j1,j2)
		do k=1,n
			omg=omega(k)
			do i=1,n
				iwb(i)=find(omega,(omg+omega(i)))
			enddo
			iwb(0)=iwb(1)
			ei1=omega(1)
			rt1=0._wp
			rt_=0._wp
			! out integration
			do i=1,n
				i1=iwb(i-1)+1
				do
					! insert point on x based on x+omega
					d(1)=0
					if(i1/=iwb(i)+1) then
						ei2=(omega(i1)-omg)
						d(1)=1
					endif
					if(d(1)==0) then
						ei2=omega(i)
						ga_=ga(i)
					else
						ga_=ga(i-1)+(ei2-omega(i-1))/(omega(i)-omega(i-1))*(ga(i)-ga(i-1))
					endif
					if(i1>1.and.i1<=n) then
						A_=A(i1-1)+((omg+ei2)-omega(i1-1))/(omega(i1)-omega(i1-1))*(A(i1)-A(i1-1))
					else
						A_=0._wp
					endif

					do j=1,n
						iwb2(j)=find(omega,omg+ei2+omega(j))
					enddo
					iwb2(0)=iwb2(1)
					ej1=omega(1)
					rt2=0._wp
					rt1_=0._wp
					! inner integration
					do j=1,n
						j1=iwb(j-1)+1
						j2=iwb2(j-1)+1
						do
							! insert point on y based on y+omega and y+x+omega
							d(2)=0
							if(j1/=iwb(j)+1) then
								ej2=(omega(j1)-omg)
								d(2)=1
							endif
							if(j2/=iwb2(j)+1) then
								if(d(2)==0.or.(d(2)/=0.and.ej2>(omega(j2)-omg-ei2))) then
									ej2=(omega(j2)-omg-ei2)
									d(2)=2
								endif
							endif
							if(d(2)==0) then
								ej2=omega(j)
								gb_=gb(j)
							else
								gb_=gb(j-1)+(ej2-omega(j-1))/(omega(j)-omega(j-1))*(gb(j)-gb(j-1))
							endif
							if(j1>1.and.j1<=n) then
								B_=B(j1-1)+(omg+ej2-omega(j1-1))/(omega(j1)-omega(j1-1))*(B(j1)-B(j1-1))
							else
								B_=0._wp
							endif
							if(j2>1.and.j2<=n) then
								C_=C(j2-1)+(omg+ei2+ej2-omega(j2-1))/(omega(j2)-omega(j2-1))*(C(j2)-C(j2-1))
							else
								C_=0._wp
							endif
							rt2_=ga_%im*gb_%im*sum(f(beta*ei2,beta*ej2,beta*omg)*[A_%im*B_%re*C_%re,A_%re*B_%im*C_%re,A_%re*B_%re*C_%im,-A_%im*B_%im*C_%im])
							rt2=rt2+0.5_wp*(ej2-ej1)*(rt1_+rt2_)
							ej1=ej2
							rt1_=rt2_
							select case(d(2))
							case(1)
								j1=j1+1
							case(2)
								j2=j2+1
							case(0)
								exit
							end select
						enddo
					enddo
					rt_=rt_+0.5_wp*(ei2-ei1)*(rt1+rt2)
					ei1=ei2
					rt1=rt2
					select case(d(1))
					case(1)
						i1=i1+1
					case(0)
						exit
					end select
				enddo
			enddo
			rt(k)%im=rt(k)%im-rt_/pi**2
		enddo
		!$OMP END PARALLEL DO
		!call stop_time(id=1,show=underscore)
	contains
		function f(x,y,z)
			real(wp) :: f(4)
			real(wp) :: x,y,z
			real(wp) :: sg,esx,esy,oez,esz,ex,ey,ez
			sg=sign(1._wp,z)
			esx=exp(sg*x)
			esy=exp(sg*y)
			esz=exp(-sg*z)
			oez=1._wp/(1._wp+esz)
			if(sg>0._wp) then
				ex=esx
				ey=esy
				!ez=esz
			else
				ex=1._wp/esx
				ey=1._wp/esy
				!ez=1._wp/esz
			endif
			!IRR
			f(1)=1._wp/((1._wp+(esx+exp((-sign(1._wp,z)*(x+z))))*oez)*(1._wp+ey))
			!RIR
			f(2)=1._wp/((1._wp+(esy+exp((-sign(1._wp,z)*(y+z))))*oez)*(1._wp+ex))
			!RRI
			f(3)=1._wp/((1._wp+esx)*(1._wp+exp(-sign(1._wp,z)*(z+x+y)))*(1._wp+esy)*oez)
			!III
			!f(4)=f(1)/((1._wp+ez/ey)*(1._wp+ez/(ex*ey)))
			f(4)=f(1)/((1._wp+exp(-z-y))*(1._wp+exp(-z-x-y)))

			!!IRR
			!f(1)=1._wp/((1._wp+(exp(sign(1._wp,z)*x)+exp((-sign(1._wp,z)*(x+z))))/(1._wp+exp(-abs(z))))*(1._wp+exp(y)))
			!!RIR
			!f(2)=1._wp/((1._wp+(exp(sign(1._wp,z)*y)+exp((-sign(1._wp,z)*(y+z))))/(1._wp+exp(-abs(z))))*(1._wp+exp(x)))
			!!RRI
			!f(3)=1._wp/((1._wp+exp(sign(1._wp,z)*x))*(1._wp+exp(-sign(1._wp,z)*(z+x+y)))*(1._wp+exp(sign(1._wp,z)*y))/(1._wp+exp(-abs(z))))
			!!III
			!f(4)=f(1)/((1._wp+exp(-z-y))*(1._wp+exp(-z-x-y)))
		end function
	end subroutine
	subroutine get_opt_conductivity(omega,Tmat,uk,rt)
		real(wp) :: omega(:)
		complex(wp) :: rt(:)
		complex(wp) :: Tmat(:)
		real(wp) :: uk(0:,:)
		integer :: i,j,k,ik,n,iwb(0:size(omega)),d
		real(wp) :: e1,e2,omg,M_phi,rt1_,rt2_
		complex(wp) :: A_,B_,rt_,Gk(2)
		n=size(omega)
		!$OMP PARALLEL DO PRIVATE(omg,iwb,A_,B_,Gk,e1,e2,rt_,rt1_,rt2_,d,j)
		do k=1,n
			omg=omega(k)
			do i=1,n
				iwb(i)=find(omega,(omg+omega(i)))
			end do
			iwb(0)=iwb(1)
			e1=omega(1)
			rt_=0._wp
			rt1_=0._wp
			do i=1,n
				j=iwb(i-1)+1
				do 
					!insert point on x based on omega+x
					d=0
					if(j/=iwb(i)+1) then
						e2=omega(j)-omg
						d=1
					endif
					if(d==0) then
						e2=omega(i)
						A_=Tmat(i)
					else
						A_=Tmat(i-1)+(e2-omega(i-1))/(omega(i)-omega(i-1))*(Tmat(i)-Tmat(i-1))
					endif
					if(j>1.and.j<=n) then
						B_=Tmat(j-1)+((omg+e2)-omega(j-1))/(omega(j)-omega(j-1))*(Tmat(j)-Tmat(j-1))
					else
						B_=0._wp
					endif
					rt2_=0._wp

					rt2_=integrate(omega_ek,self%rhov*imag(1._wp/(e2-omega_ek-A_))*imag(1._wp/(e2+omg-omega_ek-B_)))*(ff(beta*e2)-ff(beta*(e2+omg)))

					!do ik=1,size(uk,2)
						!!Gk(1)=1._wp/(e2-ek(uk(1:D,ik))+img*eta)
						!!Gk(2)=1._wp/(e2+omg-ek(uk(1:D,ik))+img*eta)
						!!Gk=Gk+Gk*[A_,B_]*Gk

						!Gk(1)=1._wp/(e2-ek(uk(1:D,ik))+img*eta-A_)
						!Gk(2)=1._wp/(e2+omg-ek(uk(1:D,ik))+img*eta-B_)
						!rt2_=rt2_+uk(0,ik)*dek(uk(1:D,ik))*imag(Gk(1))*imag(Gk(2))
					!enddo
					!!rt2_=-rt2_/nk**D*ff(-beta*e2)*ff(beta*(e2+omg))
					!rt2_=rt2_/nk**D*(ff(beta*e2)-ff(beta*(e2+omg)))

					rt_%re=rt_%re+0.5_wp*(e2-e1)*(rt1_+rt2_)
					e1=e2
					rt1_=rt2_
					select case(d)
					case(1)
						j=j+1
					case(0)
						exit
					end select
				enddo
			enddo
			!rt(k)%re=rt_%re/(pi*omg*fb(beta*omg))
			rt(k)%re=-rt_%re/(pi*omg)
		enddo
		!$OMP END PARALLEL DO
	end subroutine
	real(wp) pure function ek(k) result(rt)
		real(wp), intent(in) :: k(:)
		rt=-2._wp*sum(cos(k))
	end function
	real(wp) function dek(k) result(rt)
		real(wp), intent(in) :: k(:)
		!real(wp) :: k
		rt=sum((2._wp*sin(k))**2)
		!rt=2._wp*sin(k)
	end function
	subroutine get_convolution(omega,A,B,rt,two,Bbath)
		!A(x)B(omega+x)
		real(wp) :: omega(:)
		complex(wp) :: rt(:)
		complex(wp) :: A(:),B(:)
		logical, optional :: two,Bbath
		integer :: i,j,k,n,iwb(0:size(omega)),d
		real(wp) :: e1,e2,omg,M_phi,rt1_,rt2_
		complex(wp) :: A_,B_,rt_
		n=size(omega)
		!$OMP PARALLEL DO PRIVATE(omg,iwb,A_,B_,e1,e2,rt_,rt1_,rt2_,d,j)
		do k=1,n
			omg=omega(k)
			do i=1,n
				iwb(i)=find(omega,(omg+omega(i)))
			end do
			iwb(0)=iwb(1)
			e1=omega(1)
			rt_=0._wp
			rt1_=0._wp
			do i=1,n
				j=iwb(i-1)+1
				do 
					! insert point on x based on omega+x
					d=0
					if(j/=iwb(i)+1) then
						e2=omega(j)-omg
						d=1
					endif
					if(d==0) then
						e2=omega(i)
						A_=A(i)
					else
						A_=A(i-1)+(e2-omega(i-1))/(omega(i)-omega(i-1))*(A(i)-A(i-1))
					endif
					if(j>1.and.j<=n) then
						B_=B(j-1)+((omg+e2)-omega(j-1))/(omega(j)-omega(j-1))*(B(j)-B(j-1))
					else
						B_=0._wp
					endif
					rt2_=imag(A_)*imag(B_)*f(beta*e2,beta*omg)
					rt_%im=rt_%im+0.5_wp*(e2-e1)*(rt1_+rt2_)
					e1=e2
					rt1_=rt2_
					select case(d)
					case(1)
						j=j+1
					case(0)
						exit
					end select
				enddo
			enddo
			rt(k)%im=rt(k)%im-rt_%im/pi
		enddo
		!$OMP END PARALLEL DO
	contains
		real(wp) function f(x,y)
			real(wp) :: x,y
			if(present(two)) then
				f=(1._wp+exp(-abs(y)))/(1._wp+exp(sign(1._wp,y)*x)+exp(-sign(1._wp,y)*(x+y))+exp(-abs(y)))
			elseif(present(Bbath)) then
				f=fb(x)+ff(x+y)
			else
				f=1._wp/(1._wp+(exp(sign(1._wp,y)*x)+exp((-sign(1._wp,y)*(x+y))))/(1._wp+exp(-abs(y))))
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
