!include "../lib/utility.f90"
include "../largeN/mbroyden.f90"
module M_DMFT
	use lapack95, only : gesv
	use M_const
	use M_solution
	use M_matrix
	use M_utility
	use mkl_service
	use M_serde
	use M_fft
	implicit none
	real(8), parameter :: onepi=1d0/pi,halfpi=onepi*0.5d0,w0w1=0.25_wp
	integer, parameter :: nadd(3)=[200,0,0],n=(sum(nadd))*2,nk=256/32,D=1,nkD=nk**D
	real(8) :: U=nan,t=1d0,bcut=20d0,omega(n),Tk,beta,et=1d-1,mu
	complex(8) :: iomega(n)
	integer, parameter :: export_w=1,export_k=2,export_wk=3,set_omega=0,set_Ver=1,set_dSEk=2,set_dVer=3,set_dSus=4,set_GM=5,set_B=6,set_dGk0=7,set_Gk=8,set_Gi=9,set_Gl=10,set_SEl=11,set_Gkl=12,set_dGkl=13,set_ek=14,set_dVer_D=15,prime=20
	integer :: ipi,i0
	logical, parameter :: is_real=.false.
	type(t_serde(D)) :: map
	type gf(D,nk)
		integer, len :: D,nk
		complex(8) :: r(nk**D)
		complex(8) :: k(nk**D)
	end type
	type gfs(n,D,nk)
		integer, len :: n,D,nk
		real(8) :: k(D,nk**D)
		real(8), allocatable :: uk(:,:)
		integer, allocatable :: uk2k(:,:)
		real(8) :: rhoc(n),B(nk**D),ek(nk**D)
		complex(8) :: Gl(n),Gkl(n),dGl(n),Delta(n),SEl(n),GM(n),Ver(n,n)
		type(gf(D,nk)) :: dVer(n,n),dVer_D(n,n),Gk(n),dGk(n),dG0(n),dSE(n),dSus(n,n)
		!complex(8) :: Sus(n,nk**D),dSus(n,nk**D),Gck(n,nk**D),Gci(n,nk**D),dSEc(n,nk**D),dGck(n,nk**D)
		logical :: conv
	contains
		procedure :: self_consistent
		procedure :: set_Gfunction
		procedure :: export
	end type
contains
	subroutine export(self,ut,flag)
		class(gfs(*,*,*)) :: self
		integer :: ut(:),flag(:)
		integer :: iw,k,ik
		do k=1,size(flag)
			select case(flag(k))
			case(export_w)
				write(ut(k),"(A)")"Tk omega rGl iGl rDelta iDelta rSEl iSEl rG iG"
				do iw=1,n
					write(ut(k),"(*(e29.20e3))")Tk,omega(iw),self%Gl(iw),self%Delta(iw),self%SEl(iw),sum(self%Gk(iw)%k)/nkD
				enddo
				write(ut(k),"(x/)")
			case(export_wk)
				write(ut(k),"(A)")"Tk omega rdG idG rG iG rdVer idVer"
				do iw=1,n
					write(ut(k),"(*(e29.20e3))")Tk,omega(iw),self%dGk(iw)%k(ipi),self%Gk(iw)%k(ipi),self%dVer(iw,iw)%k(ipi)
				enddo
				write(ut(k),"(x/)")
			case(export_k)
				write(ut(k),"(A)")"Tk k rSus iSus B"
				do ik=1,nkD
					write(ut(k),"(*(e29.20e3))")Tk,self%k(:,ik),self%B(ik)
				enddo
				write(ut(k),"(x/)")
			end select
		end do
	end subroutine
	subroutine set_Gfunction(self,flag)
		class(gfs(*,*,*)) :: self
		integer :: flag(:),iw,jw,ik,l,i
		complex(wp), save :: phi(n,nkD),tmp(n,nkD)
		complex(wp) :: x(n**2),x_(size(x)),A(n,n),B(n,n)
		do l=1,size(flag)
			select case(flag(l))
			case(set_omega)
				if(is_real) then
					! real frequency
					call set_grid(omega(n/2+1:),[reset,add_linear,add_linear,add_linear],from=[0d0,max(U/2d0-5d0,0d0),1.5d0],to=[bcut,min(U/2d0+5d0,bcut),3d0],n=[nadd])
					omega(1:n/2)=-omega(n:n/2+1:-1)
					iomega=cmplx(omega,et,kind=8)
				else
					! matsubara frequency
					call set_grid(omega(n/2+1:),[reset,add_linear],from=[0d0],to=[n/2*2d0*pi*Tk],n=[n/2])
					!call set_grid(inu(m/2+1:),[reset,add_linear,add_linear],from=[0d0,m/4*2d0*pi*Tk],to=[m/4*2d0*pi*Tk,m/2*2d0*pi*Tk*16],n=[m/4,m/4])
					omega(1:n/2)=-omega(n:n/2+1:-1)+2d0*pi*Tk
					omega=omega-pi*Tk
					iomega=cmplx(0d0,omega,kind=8)
				endif
			case(set_GM)
				!self%GM=1d0/(self%Gl**2*(iomega%re-self%Delta+mu)*(iomega%re-self%Delta+mu-U))
				self%GM=1d0/(self%Gl**2*(iomega-self%Delta+mu)*(iomega-self%Delta+mu-U))
			case(set_B)
				!ik=ipi
				!$OMP PARALLEL DO
				do ik=1,size(self%uk,2)
					phi(:,self%uk2k(ik,1))=w0w1*U**2*self%GM(:)/(1d0+w0w1*U**2*self%GM(:)**2*[(self%dSus(iw,iw)%k(self%uk2k(ik,1)),iw=1,n)])
					if(is_real) then
						self%B(self%uk2k(ik,1))=integrate(omega,A=-imag(phi(:,self%uk2k(ik,1))*self%GM*[(self%dSus(iw,iw)%k(self%uk2k(ik,1)),iw=1,n)])*ff(beta*iomega%re)/pi)/Tk
					else
						self%B(self%uk2k(ik,1))=sum(phi(:,self%uk2k(ik,1))*self%GM*[(self%dSus(iw,iw)%k(self%uk2k(ik,1)),iw=1,n)])
					endif
				enddo
				!$OMP END PARALLEL DO
			case(set_dSus)
				!$OMP PARALLEL DO
				do ik=1,size(self%uk,2)
					do i=2,self%uk2k(ik,0)
						self%dGk(:)%k(self%uk2k(ik,i))=self%dGk(:)%k(self%uk2k(ik,1))
					enddo
				enddo
				!$OMP END PARALLEL DO

				!$OMP PARALLEL DO
				do iw=1,n
					call mfft(self%dGk(iw)%k(1:nkD),1,map%shap,self%dGk(iw)%r(1:nkD),1d0/product(map%shap))
					!self%dSus(iw)%r=-self%dGk(iw)%r**2
					call mfft(self%dGk(iw)%k(1:nkD),-1,map%shap,tmp(iw,:),1d0/product(map%shap))
					self%dSus(iw,iw)%r=-self%dGk(iw)%r*tmp(iw,:)
					call mfft(self%dSus(iw,iw)%r(1:nkD),-1,map%shap,self%dSus(iw,iw)%k(1:nkD))
					!self%dSus(iw)%k(1:nkD)=self%dSus(iw)%k(1:nkD)*Tk
				enddo
				!$OMP END PARALLEL DO
			case(set_dSus+prime)
				!$OMP PARALLEL DO
				do iw=1,n
					do jw=1,n
						if(iw/=jw) then
							self%dSus(iw,jw)%r=-self%dGk(iw)%r*tmp(jw,:)
							call mfft(self%dSus(iw,jw)%r(1:nkD),-1,map%shap,self%dSus(iw,jw)%k(1:nkD))
						endif
					enddo
				enddo
				!$OMP END PARALLEL DO
			case(set_Ver)
				do iw=1,n
					do jw=1,n
						self%Ver(iw,jw)=w0w1*U**2*self%GM(iw)*self%GM(jw)
					enddo
					self%Ver(iw,iw)=0d0
				enddo
			case(set_dVer)
				!$OMP PARALLEL DO
				do ik=1,size(self%uk,2)
					!if(self%B(self%uk2k(ik,1))>1._wp) then
						!do jw=1,n
							!x(1+(jw-1)*n:jw*n)=self%Ver(:,jw)
						!enddo
						!do i=1,1000
							!do jw=1,n
								!self%dVer(:,jw)%k(self%uk2k(ik,1))=x(1+(jw-1)*n:jw*n)
								!x_(1+(jw-1)*n:jw*n)=self%Ver(:,jw)
								!do iw=1,n
									!x_(1+(jw-1)*n:jw*n)=x_(1+(jw-1)*n:jw*n)+self%Ver(:,iw)*self%dSus(iw,iw)%k(self%uk2k(ik,1))*self%dVer(iw,jw)%k(self%uk2k(ik,1))
								!enddo
							!enddo
							!!call execute_command_line(delend)
							!!call execute_command_line(delstart)
							!!write(*,"(i3,(es14.4)$)")i,maxval(abs(x(1:n**2)-x_(1:n**2)))
							!if(all(abs(x(1:n**2)-x_(1:n**2))<1d-6).or.i==1000) then
								!exit
							!else
								!call mbroyden(i,x(1:n**2),x_(1:n**2),id=3)
							!endif
						!enddo
					!else
						do iw=1,n
							self%dVer(iw,iw)%k(self%uk2k(ik,1))=phi(iw,self%uk2k(ik,1))*(self%B(self%uk2k(ik,1))*self%GM(iw)-phi(iw,self%uk2k(ik,1))*self%GM(iw)**2*self%dSus(iw,iw)%k(self%uk2k(ik,1)))/(1._wp-self%B(self%uk2k(ik,1)))
						enddo
					!endif
				enddo
				!$OMP END PARALLEL DO
			case(set_dVer+3*prime)
				!!$OMP PARALLEL DO PRIVATE(A,B)
				do ik=1,size(self%uk,2)
					!self%dVer(:,:)%k(self%uk2k(ik,1))=self%Ver(:,:)
					B=self%Ver(:,:)
					do iw=1,n
						do jw=1,n
							A(iw,jw)=-self%Ver(iw,jw)*self%dSus(jw,jw)%k(self%uk2k(ik,1))
						enddo
						A(iw,iw)=A(iw,iw)+1._wp
					enddo
					call gesv(A,B)
					self%dVer(:,:)%k(self%uk2k(ik,1))=B
				enddo
				!!$OMP END PARALLEL DO
			case(set_dVer_D)
				!$OMP PARALLEL DO
				do ik=1,size(self%uk,2)
					do iw=1,n
						do jw=1,n
							!if(real(self%dSus(iw,jw)%k(self%uk2k(ik,1))*self%Ver(iw,jw))>-1._wp) then
								self%dVer_D(iw,jw)%k(self%uk2k(ik,1))=self%Ver(iw,jw)*self%dSus(iw,jw)%k(self%uk2k(ik,1))*self%Ver(iw,jw)/(1._wp+self%dSus(iw,jw)%k(self%uk2k(ik,1))*self%Ver(iw,jw))
							!else
								!x(1)=self%Ver(iw,jw)*self%dSus(iw,jw)%k(self%uk2k(ik,1))*self%Ver(iw,jw)
								!do i=1,1000
									!self%dVer_D(iw,jw)%k(self%uk2k(ik,1))=x(1)
									!x_(1)=self%Ver(iw,jw)*self%dSus(iw,jw)%k(self%uk2k(ik,1))*self%Ver(iw,jw)
									!x_(1)=x_(1)+self%Ver(iw,jw)*self%dSus(iw,jw)%k(self%uk2k(ik,1))*self%dVer_D(iw,jw)%k(self%uk2k(ik,1))
									!if(all(abs(x(1:1)-x_(1:1))<1d-6).or.i==1000) then
										!exit
									!else
										!!x(1:1)=x_(1:1)*0.05_wp+x(1:1)*(1._wp-0.05_wp)
										!call mbroyden(i,x(1:1),x_(1:1),id=4)
									!endif
								!enddo
							!endif
						enddo
					enddo
				enddo
				!$OMP END PARALLEL DO
			case(set_dSEk)
				!$OMP PARALLEL DO
				do ik=1,size(self%uk,2)
					do i=2,self%uk2k(ik,0)
						do iw=1,n
							self%dVer(iw,iw)%k(self%uk2k(ik,i))=self%dVer(iw,iw)%k(self%uk2k(ik,1))
						enddo
					enddo
				enddo
				!$OMP END PARALLEL DO
				!$OMP PARALLEL DO
				do iw=1,n
					!call mfft(self%dVer(iw,iw)%k(1:nkD),1,map%shap,self%dVer(iw,iw)%r(1:nkD),1d0/product(map%shap))
					call mfft(self%dVer(iw,iw)%k(1:nkD),-1,map%shap,self%dVer(iw,iw)%r(1:nkD),1d0/product(map%shap))
					!call mfft(self%dVer(iw,iw)%k(1:nkD),-1,map%shap,self%dSE(iw)%r(1:nkD),1d0/product(map%shap))
					call mfft(self%dGk(iw)%k(1:nkD),1,map%shap,self%dGk(iw)%r(1:nkD),1d0/product(map%shap))
					!self%dSE(iw)%r=self%dVer(iw)%r*self%dGk(iw)%r
					!self%dSE(iw)%r=Tk*self%dVer(iw,iw)%r*self%dGk(iw)%r
					!self%dSE(iw)%r=self%dSE(iw)%r*self%dGk(iw)%r
					self%dSE(iw)%r=self%dVer(iw,iw)%r*self%dGk(iw)%r
					call mfft(self%dSE(iw)%r(1:nkD),-1,map%shap,self%dSE(iw)%k(1:nkD))
				enddo
				!$OMP END PARALLEL DO
			case(set_dSEk+prime)
				!$OMP PARALLEL DO
				do ik=1,size(self%uk,2)
					do i=2,self%uk2k(ik,0)
						do iw=1,n
							do jw=1,n
								self%dVer_D(iw,jw)%k(self%uk2k(ik,i))=self%dVer_D(iw,jw)%k(self%uk2k(ik,1))
							enddo
						enddo
					enddo
				enddo
				!$OMP END PARALLEL DO

				!$OMP PARALLEL DO
				do iw=1,n
					do jw=1,n
						!if(iw/=jw) then
							call mfft(self%dVer_D(iw,jw)%k(1:nkD),-1,map%shap,self%dVer_D(iw,jw)%r(1:nkD),1d0/product(map%shap))
							!self%dSE(iw)%r=self%dSE(iw)%r+Tk*self%dVer(iw,jw)%r*self%dGk(jw)%r
							!self%dSE(iw)%r=self%dSE(iw)%r+self%dVer(iw,jw)%r*self%dGk(jw)%r
							self%dSE(iw)%r=self%dSE(iw)%r+self%dVer_D(iw,jw)%r*self%dGk(jw)%r
						!endif
					enddo
					call mfft(self%dSE(iw)%r(1:nkD),-1,map%shap,self%dSE(iw)%k(1:nkD))
				enddo
				!$OMP END PARALLEL DO
			case(set_dGk0)
				!$OMP PARALLEL DO
				do ik=1,size(self%uk,2)
					self%dG0(:)%k(self%uk2k(ik,1))=(self%Gk(:)%k(self%uk2k(ik,1))-self%Gl(:))
				enddo
				!$OMP END PARALLEL DO
			case(set_dGk0+prime)
				do ik=1,size(self%uk,2)
					self%dG0(:)%k(self%uk2k(ik,1))=-1._wp/(1._wp/(self%Gl(:)**2*(self%Delta(:)-self%ek(self%uk2k(ik,1))))+1._wp/self%Gl(:))
				enddo
			case(set_Gk)
				!$OMP PARALLEL DO
				do ik=1,size(self%uk,2)
					self%Gk(:)%k(self%uk2k(ik,1))=1d0/(iomega(:)-self%ek(self%uk2k(ik,1))-self%SEl(:))
				enddo
				!$OMP END PARALLEL DO
			case(set_Gk+prime)
				!$OMP PARALLEL DO
				do ik=1,size(self%uk,2)
					self%Gk(:)%k(self%uk2k(ik,1))=1d0/(self%Delta(:)-self%ek(self%uk2k(ik,1)))+1d0/((self%Delta(:)-self%ek(self%uk2k(ik,1)))*self%Gl(:))**2*self%dGk(:)%k(self%uk2k(ik,1))
				enddo
				!$OMP END PARALLEL DO
			case(set_Gl)
				self%Gl=0.5d0/(iomega-self%Delta+mu)+0.5d0/(iomega-self%Delta+mu-U)
			case(set_SEl)
				self%SEl=iomega-self%Delta-1d0/self%Gl
			case(set_Gkl)
				!$OMP PARALLEL DO
				do iw=1,n
					self%Gkl(iw)=0d0
					do ik=1,size(self%uk,2)
						self%Gkl(iw)=self%Gkl(iw)+self%uk2k(ik,0)*self%Gk(iw)%k(self%uk2k(ik,1))
					enddo
				enddo
				!$OMP END PARALLEL DO
				self%Gkl=self%Gkl/real(nkD,kind=8)
			case(set_dGkl)
				!$OMP PARALLEL DO
				do iw=1,n
					self%dGl(iw)=0d0
					!do ik=1,nkD
					!self%Gl(j)=self%Gl(j)+(1d0/(iomega(j)%re-ek(ik)-self%SEl(j)))
					!enddo
					do ik=1,size(self%uk,2)
						self%dGl(iw)=self%dGl(iw)+self%uk2k(ik,0)*self%dGk(iw)%k(self%uk2k(ik,1))
					enddo
				enddo
				!$OMP END PARALLEL DO
				self%dGl=self%dGl/real(nkD,kind=8)
			case(set_ek)
				!$OMP PARALLEL DO
				do ik=1,size(self%k,2)
					self%ek(ik)=0._wp
					do i=1,D
						self%ek(ik)=self%ek(ik)+2d0*t*cos(self%k(i,ik))
					enddo
				enddo
				!$OMP END PARALLEL DO
			case default
				stop "the case is not implimented"
			end select
		enddo
	end subroutine
	subroutine self_consistent(self,rate,niter,tol)
		class(gfs(*,*,*)) :: self
		real(8) :: rate,tol
		integer :: niter(2)
		integer :: i,i1,i2,ik,iw,nx,nuk
		complex(8) :: x(self%n*max(size(self%uk,2),self%n)),x_(size(x))
		complex(8) :: tmp(n)
		mu=0.5d0*U
		if(.not.is_real) call self%set_Gfunction([set_omega])
		call mbroyden(0,x,x_,rate,20,id=4)
		do i=1,iw
			self%dSE(iw)%k=0._wp
		enddo
		do i1=1,abs(niter(1))
			if(is_real) call set_realpart(iomega%re,self%Delta)
			call self%set_Gfunction([set_Gl,set_GM,set_dGk0+prime])

			!dual fermion
			nx=self%n*size(self%uk,2)
			do ik=1,size(self%uk,2)
				self%dGk(:)%k(self%uk2k(ik,1))=1._wp/(1._wp/self%dG0(:)%k(self%uk2k(ik,1))-self%dSE(:)%k(self%uk2k(ik,1)))
				x((ik-1)*n+1:ik*n)=self%dGk(:)%k(self%uk2k(ik,1))
			enddo
			do i2=1,abs(niter(2))
				call self%export([50],[export_wk])
				!call self%set_Gfunction([set_dSus,set_B,set_Ver,set_dVer,set_dSEk])
				call self%set_Gfunction([set_dSus,set_B,set_Ver,set_dVer,set_dSus+prime,set_dVer_D,set_dSEk,set_dSEk+prime])
				!exit
				!$OMP PARALLEL DO
				do ik=1,size(self%uk,2)
					x_((ik-1)*n+1:ik*n)=1._wp/(1._wp/self%dG0(:)%k(self%uk2k(ik,1))-self%dSE(:)%k(self%uk2k(ik,1)))
				enddo
				!$OMP END PARALLEL DO
				if(i2/=1) call execute_command_line(delend)
				call execute_command_line(delstart)
				write(*,"(i3,3(es14.4)$)")i2,maxval(abs(x(1:nx)-x_(1:nx))),Tk,self%B(ipi)
				if(all(abs(x(1:nx)-x_(1:nx))<tol).or.i2==abs(niter(2))) then
					exit
				else
					if(niter(2)<0) then
						x(1:nx)=x_(1:nx)*rate+x(1:nx)*(1._wp-rate)
					else
						call mbroyden(i2,x(1:nx),x_(1:nx),id=2)
					endif
				endif
				!$OMP PARALLEL DO
				do ik=1,size(self%uk,2)
					if(is_real) call set_realpart(iomega%re,x((ik-1)*n+1:ik*n))
					self%dGk(:)%k(self%uk2k(ik,1))=x((ik-1)*n+1:ik*n)
				enddo
				!$OMP END PARALLEL DO
			enddo

			call self%set_Gfunction([set_Gk+prime,set_Gkl,set_dGkl])
			!call self%set_Gfunction([set_SEl,set_Gk,set_Gkl,set_dGkl])
			nx=self%n
			x(1:nx)=self%Delta
			!x_(1:nx)=iomega-1d0/self%Gkl-self%SEl
			x_(1:nx)=self%Delta+self%dGl/(self%Gkl*self%Gl)
			if(is_real) then
				x(1:nx)%re=0._wp
				x_(1:nx)%re=0._wp
			endif
			write(*,"(i3,*(es14.4))")i1,maxval(abs(x(1:nx)-x_(1:nx))),self%B(ipi)
			if(all(abs(x(1:nx)-x_(1:nx))<tol).or.i1==abs(niter(1))) then
				self%conv=(i1/=abs(niter(1)))
				exit
			else
				if(niter(1)<0) then
					x(1:nx)=x_(1:nx)*rate+x(1:nx)*(1._wp-rate)
				else
					call mbroyden(i1,x(1:nx),x_(1:nx),id=1)
				endif
			endif
			self%Delta=x(1:nx)
			call self%export([10],[export_w])
		enddo
		call self%set_Gfunction([set_dSus,set_B])
		call mbroyden(huge(1),x,x_)
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
	subroutine get_uk(k,uk,uk2k)
		real(8) :: k(:,:)
		real(8), allocatable :: uk(:,:)
		integer, allocatable :: uk2k(:,:)
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
		allocate(uk(1:D,size(c)-1))
		uk(1:D,:)=tmp(:,c(1:size(c)-1))
		!uk(0,:)=real(c(2:)-c(1:size(c)-1),kind=8)
		allocate(uk2k(size(uk,2),0:maxval(c(2:)-c(1:size(c)-1))))
		uk2k(:,0)=c(2:)-c(1:size(c)-1)
		do i=1,size(uk,2)
			uk2k(i,1:c(i+1)-c(i))=ord(c(i):c(i+1)-1)
		enddo
		
		write(*,*)sum(uk2k(:,0))==nkD
	end subroutine
end module
include "../lib/serde.f90"
include "../lib/fft.f90"
program main
	use M_DMFT
	implicit none
	type(gfs(n,D,nk)) :: phy
	integer :: j,i,ik,jk,iT,iU
	real(8) :: ff_(n),dTk

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
	open(unit=50,File="../data/fk_dG.dat")
	!call add_omega(omega(1:nl*nlog*2+sum(nadd(1:2))*4),4d0,0.5d0,nadd(2))
	map%shap=[(nk,i=1,D)]
	ipi=0
	i0=0
	!!$OMP PARALLEL DO
	do ik=1,size(phy%k,2)
		phy%k(:,ik)=2d0*pi*(map%get_idx(ik,[1:D])-1)/nk
		if(sum(abs(phy%k(:,ik)-pi))<1d-8) then
			if(ipi==0) then
				ipi=ik
			else
				stop "err"
			endif
		elseif(sum(abs(phy%k(:,ik)))<1d-8) then
			if(i0==0) then
				i0=ik
			else
				stop "err"
			endif
		endif
	enddo
	!!$OMP END PARALLEL DO
	write(*,*)"stating get unique k point...."
	call get_uk(phy%k,phy%uk,phy%uk2k)
	write(*,*)"finished get unique k point",size(phy%k,2)/size(phy%uk,2)
	!do i=1,size(phy%uk,2)
		!write(*,*)phy%uk(:,i)
		!read(*,*)
	!enddo
	!do i=1,n
		!read(20,"(2(e28.20))")phy%Delta(i)
	!enddo
	write(40,"(A)")"U T"
	call phy%set_Gfunction([set_ek])
	!phy%Delta=1d-1
	phy%Delta=cmplx(1d-1,1d-1)
	Tk=0.207046875000000d0
	!5.00000000000000       0.125046875000000       0.991754771775755
	!5.00000000000000       0.124546875000000        1.15946531278591

	do 
		U=for_in([8]*1d0,id=1)
		if(isnan(U)) exit
		if(next(Tk,0d0,1d-1,1d-4)) then
		endif
		if(is_real) then
			call phy%set_Gfunction([set_omega])
		endif
		do 
			beta=1._wp/Tk

			call phy%self_consistent(0.2d0,[300,0],1d-6)
			!call phy%self_consistent(0.05d0,[300,1000],1d-6)
			!stop
			if(next(Tk,phy%B(ipi)-1d0)) then
				write(*,*)U,Tk,phy%B(ipi)
				write(40,"(*(e29.20e3))")U,Tk
				exit
			endif
		enddo
		!cycle
		if(next(Tk,0d0,1d-3,1d-4)) then
		endif
		do 
			beta=1._wp/Tk

			!call phy%self_consistent(0.05d0,[300,0],1d-6)
			call phy%self_consistent(0.05d0,[100,100],1d-6)
			write(*,*)U,Tk,phy%B(ipi)
			if(next(Tk,phy%B(ipi)-1d0)) then
				write(40,"(*(e29.20e3))")U,Tk
				exit
			endif
		enddo
	enddo
end program
