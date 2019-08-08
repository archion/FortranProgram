!include "../lib/utility.f90"
!include "../test/pade.f90"
include "../largeN/mbroyden.f90"
module M_DMFT
	use lapack95, only : gesv, rsteqr 
	use M_const
	use M_solution
	use M_matrix
	use M_utility
	use mkl_service
	use M_serde
	use M_fft
	!use M_pade
	implicit none
	real(8), parameter :: onepi=1d0/pi,halfpi=onepi*0.5d0,&
		nf=1._wp/3._wp,nc=2._wp/3._wp,&
		!nf=1._wp/4._wp,nc=3._wp/4._wp,&
		!nf=1._wp/7._wp,nc=6._wp/7._wp,&
		!nf=1._wp/5._wp,nc=4._wp/5._wp,&
		!nf=1._wp/2._wp,nc=1._wp/2._wp,&
		!nf=1._wp/2._wp,nc=0.25_wp,&
		!nf=0.82_wp,nc=1._wp-nf,&
		w0w1=(1._wp-nf)*nf
	integer, parameter :: nadd(3)=[300,0,0],n=(sum(nadd))*2,&
		!nk=256/8,D=2,&
		nk=24*2,D=2,&
		nkD=nk**D
	real(8) :: U=nan,t=1d0,bcut=10d0,omega(n),Tk,beta,et=1d-1,mu,a(D,D),Res(n)
	complex(8) :: iomega(n)
	integer, parameter :: export_w=1,export_k=2,export_wk=3,export_check=4,set_omega=0,set_Ver=1,set_dSEk=2,set_dVer=3,set_dSus=4,set_GM=5,set_B=6,set_dGk0=7,set_Gk=8,set_Gi=9,set_Gl=10,set_SEl=11,set_Gkl=12,set_dGkl=13,set_ek=14,set_dVer_D=15,set_nc=16,set_Sus_static=17,set_dGk=18,set_Delta=19,set_phi=20,set_dGr=21,prime=25
	logical :: need_update(-2*prime:4*prime)=.true.
	integer :: ipi,i0,iorder
	logical, parameter :: is_real=.false.,is_tri=.true.,is_normal_mat=.false.
	type(t_serde(D)) :: map
	type gf(D,nk)
		integer, len :: D,nk
		complex(8) :: r(nk**D)
		complex(8) :: k(nk**D)
	end type
	type gf2(D,nk)
		integer, len :: D,nk
		complex(8) :: r(nk**D)
		complex(8) :: rp(nk**D)
		complex(8) :: k(nk**D)
	end type
	type gfk(D,nk)
		integer, len :: D,nk
		complex(8) :: k(nk**D)
	end type
	type gfs(n,D,nk)
		integer, len :: n,D,nk
		real(8) :: k(D,nk**D)
		real(8) :: nc
		real(8), allocatable :: uk(:,:)
		integer, allocatable :: uk2k(:,:)
		real(8) :: rhoc(n),B(nk**D),ek(nk**D),Sus(nk**D)
		complex(8) :: Gl(n),Gkl(n),dGl(n),Delta(n),SEl(n),GM(n),Ver(n,n)
		type(gf(D,nk)) :: Gk(n),dG0(n),dSE(n),phi(n)
		type(gf2(D,nk)) :: dGk(n)
		type(gfk(D,nk)) :: dVer(n,n),dVer_D(n,n),dSus(n,n),Sus0(n)
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
		integer :: iw,k,ik,i
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
					write(ut(k),"(*(e29.20e3))")Tk,omega(iw),self%dGk(iw)%k(self%uk2k(iorder,1)),self%Gk(iw)%k(self%uk2k(iorder,1)),self%dVer(iw,iw)%k(self%uk2k(iorder,1))
				enddo
				write(ut(k),"(x/)")
			case(export_k)
				do ik=1,size(self%uk,2)
					do i=2,self%uk2k(ik,0)
						self%B(self%uk2k(ik,i))=self%B(self%uk2k(ik,1))
						self%Sus(self%uk2k(ik,i))=self%Sus(self%uk2k(ik,1))
					enddo
				enddo
				write(ut(k),"(A)")"U Tk k1 k2 B Sus ek"
				do ik=1,nkD
					write(ut(k),"(*(e29.20e3))")U,Tk,self%k(:,ik)/pi,self%B(ik),self%Sus(ik),self%ek(ik)
				enddo
				write(ut(k),"(x/)")
			case(export_check)
				write(ut(k),"(A)")"i kx1 ky1 kx2 ky2 kx3 ky3 kx4 ky4 kx5 ky5 kx6 ky6 kx7 ky7 kx8 ky8"
				do ik=1,size(self%uk,2)
					write(ut(k),"(i5$)")self%uk2k(ik,0)
					do i=1,self%uk2k(ik,0)
						write(ut(k),"(e29.20e3$)")self%uk(:,ik)/pi,self%k(:,self%uk2k(ik,i))/pi
					enddo
					write(ut(k),"(x)")
				enddo
				write(ut(k),"(x/)")
			end select
		end do
	end subroutine
	recursive subroutine set_Gfunction(self,flag)
		class(gfs(*,*,*)) :: self
		integer :: flag(:),iw,jw,ik,l,i,iord(n)
		!complex(wp) :: x(n**2),x_(size(x)),A(n,n),B(n,n)
		complex(wp) :: r(nkD),rp(nkD),nc,y(1)
		complex(wp), save :: LL(n,nkD)
		real(wp) :: maxv,op,U_(n,n)
		do l=1,size(flag)
			select case(flag(l))
			case(set_omega)
				need_update=.true.
				if(is_real) then
					! real frequency
					call set_grid(omega(n/2+1:),[reset,add_linear,add_linear,add_linear],from=[0d0,max(U/2d0-5d0,0d0),1.5d0],to=[bcut,min(U/2d0+5d0,bcut),3d0],n=[nadd])
					omega(1:n/2)=-omega(n:n/2+1:-1)
					iomega=cmplx(omega,et,kind=8)
				else
					! matsubara frequency
					if(is_normal_mat) then
						call set_grid(omega(n/2+1:),[reset,add_linear],from=[0d0],to=[nadd(1)*2d0*pi],n=[nadd(1)])
						do i=1,nadd(2)
							omega(n/2+nadd(1)+i)=omega(n/2+nadd(1)+i-1)+2d0*pi*2**(min(i-1,9))
						enddo
						!call set_grid(inu(m/2+1:),[reset,add_linear,add_linear],from=[0d0,m/4*2d0*pi],to=[m/4*2d0*pi,m/2*2d0*pi*16],n=[m/4,m/4])
						omega(1:n/2)=-omega(n:n/2+1:-1)+2d0*pi
						omega=omega-pi
						Res=-1._wp
					else
						! continue fraction
						do i=1,n-1
							Res(i)=0.5_wp/(sqrt((2._wp*i-1._wp)*(2._wp*i+1._wp)))
						enddo
						omega=0._wp
						call rsteqr(omega,Res(:n-1),U_)
						omega=1._wp/omega
						iord=[1:n]
						call qsort(omega,iord)
						Res=-0.25_wp*(U_(1,:)*omega)**2
						omega=omega(iord)
						Res=Res(iord)

						!! continue fraction
						!do i=1,n-3
							!Res(i)=0.5_wp/(sqrt((2._wp*i-1._wp)*(2._wp*i+1._wp)))
						!enddo
						!omega=0._wp
						!call rsteqr(omega(2:n-1),Res(:n-3),U_(:n-2,:n-2))
						!omega(2:n-1)=1._wp/omega(2:n-1)
						!Res([1,n])=0._wp
						!omega(n)=max(omega(n-1)*2._wp,1e5_wp)
						!omega(1)=-omega(n)
						!iord=[1:n]
						!call qsort(omega,iord)
						!Res(2:n-1)=-0.25_wp*(U_(1,:n-2)*omega(2:n-1))**2
						!omega=omega(iord)
						!Res=Res(iord)
					endif
				endif
			case(set_ek)
				need_update=.true.
				if(is_tri) then
					! triangular
					!$OMP PARALLEL DO
					do ik=1,size(self%k,2)
						self%ek(ik)=-2d0*t*(cos(self%k(1,ik))+cos(0.5_wp*(self%k(1,ik)+sqrt(3._wp)*self%k(2,ik)))+cos(0.5_wp*(-self%k(1,ik)+sqrt(3._wp)*self%k(2,ik))))
					enddo
					!$OMP END PARALLEL DO
				else
					! cubic
					!$OMP PARALLEL DO
					do ik=1,size(self%k,2)
						self%ek(ik)=0._wp
						do i=1,D
							self%ek(ik)=self%ek(ik)-2d0*t*cos(self%k(i,ik))
						enddo
					enddo
					!$OMP END PARALLEL DO
				endif
			case(set_nc)
				if(is_normal_mat) then
					call self%set_Gfunction([set_Gl])
					if(.not.need_update(flag(l))) cycle
					self%nc=real(msum(Tk*exp(iomega*0.0001_wp)*(self%Gl-&
						((1._wp-nf)/(iomega+mu)+nf/(iomega+mu-U))))+(1._wp-nf)*ff(-mu*beta)+nf*ff((-mu+U)*beta))
				else
					call self%set_Gfunction([set_Gkl])
					if(.not.need_update(flag(l))) cycle
					self%nc=Tk*real(msum(self%Gkl,not_conv=underscore))
				endif
			case(set_dGk0)
				call self%set_Gfunction([set_Gl,set_Delta])
				if(.not.need_update(flag(l))) cycle
				!$OMP PARALLEL DO
				do ik=1,size(self%uk,2)
					self%dG0(:)%k(self%uk2k(ik,1))=-1._wp/(1._wp/(self%Gl(:)**2*(self%Delta(:)-self%ek(self%uk2k(ik,1))))+1._wp/self%Gl(:))
				enddo
				!$OMP END PARALLEL DO
				need_update([set_dGk])=.true.
			case(set_dGk0+prime)
				call self%set_Gfunction([set_Gk,set_Gl])
				if(.not.need_update(flag(l))) cycle
				!$OMP PARALLEL DO
				do ik=1,size(self%uk,2)
					self%dG0(:)%k(self%uk2k(ik,1))=(self%Gk(:)%k(self%uk2k(ik,1))-self%Gl(:))
				enddo
				!$OMP END PARALLEL DO
				need_update([set_dGk])=.true.
			case(set_Gk)
				call self%set_Gfunction([set_Gl,set_Delta,set_dGk])
				if(.not.need_update(flag(l))) cycle
				!$OMP PARALLEL DO
				do ik=1,size(self%uk,2)
					self%Gk(:)%k(self%uk2k(ik,1))=1d0/(self%Delta(:)-self%ek(self%uk2k(ik,1)))+1d0/((self%Delta(:)-self%ek(self%uk2k(ik,1)))*self%Gl(:))**2*self%dGk(:)%k(self%uk2k(ik,1))
				enddo
				!$OMP END PARALLEL DO
				need_update([set_dGk0+prime,set_Gkl,set_Sus_static])=.true.
			case(set_Gk+prime)
				call self%set_Gfunction([set_SEl])
				if(.not.need_update(flag(l))) cycle
				!$OMP PARALLEL DO
				do ik=1,size(self%uk,2)
					self%Gk(:)%k(self%uk2k(ik,1))=1d0/(iomega(:)-self%ek(self%uk2k(ik,1))-self%SEl(:))
				enddo
				!$OMP END PARALLEL DO
				need_update([set_dGk0+prime,set_Gkl,set_Sus_static])=.true.
			case(set_Gl)
				call self%set_Gfunction([set_Delta])
				if(.not.need_update(flag(l))) cycle
				self%Gl=(1._wp-nf)/(iomega-self%Delta+mu)+nf/(iomega-self%Delta+mu-U)
				!self%Gl=1._wp/(iomega+.1_wp)
				need_update([set_dGk0,set_dGk0+prime,set_Gk,set_SEl,set_GM,set_Sus_static])=.true.
				if(is_normal_mat) need_update(set_nc)=.true.
			case(set_SEl)
				call self%set_Gfunction([set_Gl,set_Delta])
				if(.not.need_update(flag(l))) cycle
				self%SEl=iomega-self%Delta-1d0/self%Gl
				need_update([set_Gk+prime])=.true.
			case(set_Gkl)
				call self%set_Gfunction([set_Gk])
				if(.not.need_update(flag(l))) cycle
				!$OMP PARALLEL DO
				do iw=1,n
					self%Gkl(iw)=0d0
					do ik=1,size(self%uk,2)
						self%Gkl(iw)=self%Gkl(iw)+self%uk2k(ik,0)*self%Gk(iw)%k(self%uk2k(ik,1))
					enddo
				enddo
				!$OMP END PARALLEL DO
				self%Gkl=self%Gkl/real(nkD,kind=8)
				if(.not.is_normal_mat) need_update(set_nc)=.true.
			case(set_dGk)
				need_update([set_Gk,set_dGkl,set_dGr,set_Sus_static])=.true.
			case(set_dGkl)
				call self%set_Gfunction([set_dGk])
				if(.not.need_update(flag(l))) cycle
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
			case(set_Delta)
				need_update([set_Gl,set_GM,set_SEl,set_Sus_static,set_Gk,set_dGk0])=.true.
			case(set_GM)
				call self%set_Gfunction([set_Gl,set_Delta])
				if(.not.need_update(flag(l))) cycle
				!self%GM=1d0/(self%Gl**2*(iomega%re-self%Delta+mu)*(iomega%re-self%Delta+mu-U))
				self%GM=1d0/(self%Gl**2*(iomega-self%Delta+mu)*(iomega-self%Delta+mu-U))
				need_update([set_B,set_Ver,set_dVer,set_dVer+prime,set_phi])=.true.
			case(set_phi)
				call self%set_Gfunction([set_GM,set_dSus])
				if(.not.need_update(flag(l))) cycle
				!$OMP PARALLEL DO
				do ik=1,size(self%uk,2)
					self%phi(:)%k(self%uk2k(ik,1))=w0w1*U**2*self%GM(:)/(1d0+w0w1*U**2*self%GM(:)**2*[(self%dSus(iw,iw)%k(self%uk2k(ik,1)),iw=1,n)])
				enddo
				!$OMP END PARALLEL DO
				need_update([set_B,set_dVer,set_dVer+prime])=.true.
			case(set_B)
				call self%set_Gfunction([set_GM,set_dSus,set_phi])
				if(.not.need_update(flag(l))) cycle
				!ik=ipi
				maxv=-(huge(1._wp)-1._wp)
				!$OMP PARALLEL DO
				do ik=1,size(self%uk,2)
					if(is_real) then
						self%B(self%uk2k(ik,1))=integrate(omega,A=-imag(self%phi(:)%k(self%uk2k(ik,1))*self%GM*[(self%dSus(iw,iw)%k(self%uk2k(ik,1)),iw=1,n)])*ff(beta*iomega%re)/pi)/Tk
					else
						!self%B(self%uk2k(ik,1))=sum(self%phi(:)%k(self%uk2k(ik,1))*self%GM*[(self%dSus(iw,iw)%k(self%uk2k(ik,1)),iw=1,n)])
						self%B(self%uk2k(ik,1))=msum(self%phi(:)%k(self%uk2k(ik,1))*self%GM*[(self%dSus(iw,iw)%k(self%uk2k(ik,1)),iw=1,n)])
						!if(abs(0.5_wp*self%phi(n)%k(self%uk2k(ik,1))*self%GM(n)*(self%dSus(n,n)%k(self%uk2k(ik,1))*iomega(n)))>1e-6_wp) then
							!write(*,*)abs(0.5_wp*self%phi(n)%k(self%uk2k(ik,1))*self%GM(n)*(self%dSus(n,n)%k(self%uk2k(ik,1))*iomega(n)))
						!endif
					endif
					!$OMP CRITICAL
					if(maxv<self%B(self%uk2k(ik,1))) then
						maxv=self%B(self%uk2k(ik,1))
						iorder=ik
					endif
					!$OMP END CRITICAL
				enddo
				!$OMP END PARALLEL DO
				need_update([set_dVer,set_dVer+prime])=.true.
			case(set_dGr)
				call self%set_Gfunction([set_dGk])
				if(.not.need_update(flag(l))) cycle
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
					call mfft(self%dGk(iw)%k(1:nkD),-1,map%shap,self%dGk(iw)%rp(1:nkD),1d0/product(map%shap))
				enddo
				need_update([set_dSus,set_dSus+prime,set_dSEk])=.true.
			case(set_dSus)
				call self%set_Gfunction([set_dGr])
				if(.not.need_update(flag(l))) cycle
				!$OMP PARALLEL DO PRIVATE(r)
				do iw=1,n
					r=-self%dGk(iw)%r*self%dGk(iw)%rp
					call mfft(r(1:nkD),-1,map%shap,self%dSus(iw,iw)%k(1:nkD))
					!self%dSus(iw)%k(1:nkD)=self%dSus(iw)%k(1:nkD)*Tk
				enddo
				!$OMP END PARALLEL DO
				need_update([set_B,set_dVer,set_dVer+prime,set_dVer+3*prime,set_dVer_D,set_phi,set_dSus+prime])=.true.
			case(set_dSus+prime)
				call self%set_Gfunction([set_dSus,set_dGr])
				if(.not.need_update(flag(l))) cycle
				!$OMP PARALLEL DO PRIVATE(r)
				do iw=1,n
					do jw=1,n
						if(iw/=jw) then
							r=-self%dGk(iw)%r*self%dGk(jw)%rp
							call mfft(r(1:nkD),-1,map%shap,self%dSus(iw,jw)%k(1:nkD))
						endif
					enddo
				enddo
				!$OMP END PARALLEL DO
				need_update(set_dVer_D)=.true.
			case(set_Ver)
				call self%set_Gfunction([set_GM])
				if(.not.need_update(flag(l))) cycle
				do iw=1,n
					do jw=1,n
						self%Ver(iw,jw)=w0w1*U**2*self%GM(iw)*self%GM(jw)
					enddo
					self%Ver(iw,iw)=0d0
				enddo
				need_update([set_dVer+3*prime,set_dVer_D])=.true.
			case(set_dVer)
				call self%set_Gfunction([set_GM,set_B,set_dSus])
				if(.not.need_update(flag(l))) cycle
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
								!call mbroyden(i,x(1:n**2),x_(1:n**2)-x(1:n**2),id=3)
							!endif
						!enddo
					!else
						do iw=1,n
							self%dVer(iw,iw)%k(self%uk2k(ik,1))=self%phi(iw)%k(self%uk2k(ik,1))*(self%B(self%uk2k(ik,1))*self%GM(iw)-self%phi(iw)%k(self%uk2k(ik,1))*self%GM(iw)**2*self%dSus(iw,iw)%k(self%uk2k(ik,1)))/(1._wp-self%B(self%uk2k(ik,1)))
						enddo
					!endif
				enddo
				!$OMP END PARALLEL DO
				need_update([set_dSEk,set_dSEk+prime,set_Sus_static])=.true.
			case(set_dVer+prime)
				call self%set_Gfunction([set_GM,set_B,set_dSus])
				if(.not.need_update(flag(l))) cycle
				!$OMP PARALLEL DO
				do ik=1,size(self%uk,2)
					do iw=1,n
						do jw=1,n
							if(iw/=jw) then
								self%dVer(iw,jw)%k(self%uk2k(ik,1))=self%phi(iw)%k(self%uk2k(ik,1))*(self%GM(jw)-self%phi(jw)%k(self%uk2k(ik,1))*self%GM(jw)**2*self%dSus(jw,jw)%k(self%uk2k(ik,1)))/(1._wp-self%B(self%uk2k(ik,1)))
							endif
						enddo
					enddo
				enddo
				!$OMP END PARALLEL DO
				need_update(set_Sus_static)=.true.
			!case(set_dVer+3*prime)
				!call self%set_Gfunction([set_dSus,set_Ver])
				!if(.not.need_update(flag(l))) cycle
				!if(nadd(2)/=0) then
					!write(*,*)"warning! logarithm matsubara grid!!"
				!endif
				!!!$OMP PARALLEL DO PRIVATE(A,B)
				!do ik=1,size(self%uk,2)
					!!self%dVer(:,:)%k(self%uk2k(ik,1))=self%Ver(:,:)
					!B=self%Ver(:,:)
					!do iw=1,n
						!do jw=1,n
							!A(iw,jw)=-self%Ver(iw,jw)*self%dSus(jw,jw)%k(self%uk2k(ik,1))
						!enddo
						!A(iw,iw)=A(iw,iw)+1._wp
					!enddo
					!call gesv(A,B)
					!self%dVer(:,:)%k(self%uk2k(ik,1))=B
				!enddo
				!!!$OMP END PARALLEL DO
				!need_update([set_dSEk,set_dSEk+prime,set_Sus_static])=.true.
			case(set_dVer_D)
				call self%set_Gfunction([set_dSus+prime,set_Ver])
				if(.not.need_update(flag(l))) cycle
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
										!call mbroyden(i,x(1:1),x_(1:1)-x(1:1),id=4)
									!endif
								!enddo
							!endif
						enddo
					enddo
				enddo
				!$OMP END PARALLEL DO
				need_update(set_dSEk+prime)=.true.
			case(set_dSEk)
				call self%set_Gfunction([set_dGr,set_dVer])
				if(.not.need_update(flag(l))) cycle
				!$OMP PARALLEL DO
				do ik=1,size(self%uk,2)
					do i=2,self%uk2k(ik,0)
						do iw=1,n
							self%dVer(iw,iw)%k(self%uk2k(ik,i))=self%dVer(iw,iw)%k(self%uk2k(ik,1))
						enddo
					enddo
				enddo
				!$OMP END PARALLEL DO
				!$OMP PARALLEL DO PRIVATE(r)
				do iw=1,n
					!call mfft(self%dVer(iw,iw)%k(1:nkD),-1,map%shap,self%dVer(iw,iw)%r(1:nkD),1d0/product(map%shap))
					call mfft(self%dVer(iw,iw)%k(1:nkD),-1,map%shap,r(1:nkD),1d0/product(map%shap))
					!call mfft(self%dVer(iw,iw)%k(1:nkD),-1,map%shap,self%dSE(iw)%r(1:nkD),1d0/product(map%shap))
					!self%dSE(iw)%r=self%dVer(iw)%r*self%dGk(iw)%r
					!self%dSE(iw)%r=Tk*self%dVer(iw,iw)%r*self%dGk(iw)%r
					!self%dSE(iw)%r=self%dSE(iw)%r*self%dGk(iw)%r
					self%dSE(iw)%r=r*self%dGk(iw)%r
					call mfft(self%dSE(iw)%r(1:nkD),-1,map%shap,self%dSE(iw)%k(1:nkD))
				enddo
				!$OMP END PARALLEL DO
				need_update([set_dGk,set_dSEk+prime])=.true.
			case(set_dSEk+prime)
				call self%set_Gfunction([set_dGr,set_dVer_D,set_dSEk])
				if(.not.need_update(flag(l))) cycle
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

				!$OMP PARALLEL DO PRIVATE(LL)
				do iw=1,n
					do jw=1,n
						!if(iw/=jw) then
							!call mfft(self%dVer_D(iw,jw)%k(1:nkD),-1,map%shap,r(1:nkD),1d0/product(map%shap))
							call mfft(self%dVer_D(iw,jw)%k(1:nkD),-1,map%shap,LL(jw,1:nkD),1d0/product(map%shap))
							!self%dSE(iw)%r=self%dSE(iw)%r+Tk*self%dVer(iw,jw)%r*self%dGk(jw)%r
							!self%dSE(iw)%r=self%dSE(iw)%r+self%dVer(iw,jw)%r*self%dGk(jw)%r
						!endif
					enddo
					do ik=1,nkD
						self%dSE(iw)%r(ik)=self%dSE(iw)%r(ik)+msum(LL(:,ik)*self%dGk(:)%r(ik))
					enddo
					call mfft(self%dSE(iw)%r(1:nkD),-1,map%shap,self%dSE(iw)%k(1:nkD))
				enddo
				!$OMP END PARALLEL DO
				need_update(set_dGk)=.true.
			case(set_Sus_static)
				call self%set_Gfunction([set_Gk,set_Gl,set_dGk,set_dVer,set_dVer+prime,set_Delta])
				if(.not.need_update(flag(l))) cycle
				!$OMP PARALLEL DO
				do ik=1,size(self%uk,2)
					do i=2,self%uk2k(ik,0)
						self%dGk(:)%k(self%uk2k(ik,i))=self%dGk(:)%k(self%uk2k(ik,1))
						self%Gk(:)%k(self%uk2k(ik,i))=self%Gk(:)%k(self%uk2k(ik,1))
						self%dVer(:,:)%k(self%uk2k(ik,i))=self%dVer(:,:)%k(self%uk2k(ik,1))
					enddo
				enddo
				!$OMP END PARALLEL DO
				do iw=1,n
					LL(iw,:)=self%dGk(iw)%k(:)/((self%Delta(iw)-self%ek)*self%Gl(iw))
					call mfft(LL(iw,1:nkD),1,map%shap,r,1d0/product(map%shap))
					call mfft(LL(iw,1:nkD),-1,map%shap,rp,1d0/product(map%shap))
					call mfft(r(1:nkD)*rp(1:nkD),-1,map%shap,LL(iw,1:nkD))
					call mfft(self%Gk(iw)%k(1:nkD),1,map%shap,r(1:nkD),1d0/product(map%shap))
					call mfft(self%Gk(iw)%k(1:nkD),-1,map%shap,rp(1:nkD),1d0/product(map%shap))
					call mfft(-r(1:nkD)*rp(1:nkD),-1,map%shap,self%Sus0(iw)%k(1:nkD))
				enddo
				do ik=1,size(self%uk,2)
					!self%Sus(self%uk2k(ik,1))=Tk*sum(self%Sus0(1:n)%k(self%uk2k(ik,1))+LL(1:n,self%uk2k(ik,1))*matmul(self%dVer(1:n,1:n)%k(self%uk2k(ik,1)),LL(1:n,self%uk2k(ik,1))))
					self%Sus(self%uk2k(ik,1))=Tk*msum(self%Sus0(1:n)%k(self%uk2k(ik,1))+LL(1:n,self%uk2k(ik,1))*[(msum(self%dVer(iw,1:n)%k(self%uk2k(ik,1))*LL(1:n,self%uk2k(ik,1))),iw=1,n)])
					!write(*,*)abs(0.5_wp*self%Sus0(n)%k(self%uk2k(ik,1))*iomega(n))/Tk,abs(0.5_wp*self%dVer(n/2,n)%k(self%uk2k(ik,1))*LL(n,self%uk2k(ik,1))*iomega(n))/Tk
					!write(*,*)sum(self%Sus0(:)%k(self%uk2k(ik,1))),sum(Res*self%Sus0(:)%k(self%uk2k(ik,1)))
				enddo
			case default
				stop "the case is not implimented"
			end select
			!write(*,*)flag(l)
			need_update(flag(l))=.false.
		enddo
	end subroutine
	subroutine self_consistent(self,rate,niter,tol)
		class(gfs(*,*,*)) :: self
		real(8) :: rate,tol
		integer :: niter(2)
		integer :: i,i1,i2,ik,iw,nx,nuk
		complex(8) :: x(self%n*max(size(self%uk,2),self%n)),x_(size(x))
		complex(8) :: tmp(n)
		need_update=.true.
		need_update(set_omega)=.false.
		if(.not.is_real) iomega=cmplx(0d0,omega*Tk,kind=8)
		call mbroyden(0,x,x_,rate,20,id=4)
		do i1=1,abs(niter(1))
			if(is_real) call set_realpart(iomega%re,self%Delta)
			call self%set_Gfunction([set_dGk0])

			!dual fermion
			nx=self%n*size(self%uk,2)
			do ik=1,size(self%uk,2)
				self%dGk(:)%k(self%uk2k(ik,1))=1._wp/(1._wp/self%dG0(:)%k(self%uk2k(ik,1))-self%dSE(:)%k(self%uk2k(ik,1)))
				x((ik-1)*n+1:ik*n)=self%dGk(:)%k(self%uk2k(ik,1))
			enddo
			call self%set_Gfunction([set_dGk])
			do i2=1,abs(niter(2))
				!call self%export([50],[export_wk])
				call self%set_Gfunction([set_dSEk])
				!call self%set_Gfunction([set_dSEk+prime])
				!exit
				!$OMP PARALLEL DO
				do ik=1,size(self%uk,2)
					x_((ik-1)*n+1:ik*n)=1._wp/(1._wp/self%dG0(:)%k(self%uk2k(ik,1))-self%dSE(:)%k(self%uk2k(ik,1)))-x((ik-1)*n+1:ik*n)
				enddo
				!$OMP END PARALLEL DO
				if(is_real) then
					x(1:nx)%re=0._wp
					x_(1:nx)%re=0._wp
				endif
				!if(i2/=1) call execute_command_line(delend)
				!call execute_command_line(delstart)
				!write(*,"(i3,3(es14.4)$)")i2,maxval(abs(x(1:nx)-x_(1:nx))),Tk,self%B(phy%uk2k(iorder,1))
				if(all(abs(x_(1:nx))<tol).or.i2==abs(niter(2))) then
					if(i2==abs(niter(2))) write(*,*)"warning2",maxval(abs(x_(1:nx)))
					exit
				else
					if(niter(2)<0) then
						x(1:nx)=x_(1:nx)*rate+x(1:nx)
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
				call self%set_Gfunction([set_dGk])
			enddo

			call self%set_Gfunction([set_Gk,set_Gkl,set_dGkl])
			!call self%set_Gfunction([set_SEl,set_Gk+prime,set_Gkl])
			nx=self%n
			x(1:nx)=self%Delta
			x_(1:nx)=self%Delta+self%dGl/(self%Gkl*self%Gl)-x(1:nx)
			!x_(1:nx)=iomega-1d0/self%Gkl-self%SEl
			if(is_real) then
				x(1:nx)%re=0._wp
				x_(1:nx)%re=0._wp
			endif
			if(all(abs(x_(1:nx))<tol).or.i1==abs(niter(1))) then
				self%conv=(i1/=abs(niter(1)))
				if(i1==abs(niter(1))) write(*,*)"warning1 ",maxval(abs(x_(1:nx)))
				exit
			else
				if(niter(1)<0) then
					x(1:nx)=x_(1:nx)*rate+x(1:nx)
				else
					call mbroyden(i1,x(1:nx),x_(1:nx),id=1)
				endif
			endif
			self%Delta=x(1:nx)
			call self%set_Gfunction([set_Delta])
		enddo
		!call self%export([10],[export_w])
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
	subroutine gen_k(k,uk,uk2k)
		real(8) :: k(:,:)
		real(8), allocatable :: uk(:,:)
		integer, allocatable :: uk2k(:,:)
		integer :: ik,i
		map%shap=[(nk,i=1,D)]
		!ipi=0
		!i0=0
		if(D==1) then
			if(is_tri) stop "err: D not equal 2 will is_tri is true"
			!a(1,1)=2._wp*pi
		else
			if(is_tri) then
				a(:,1)=[1._wp,1._wp/sqrt(3._wp)]*2._wp*pi
				a(:,2)=[1._wp,-1._wp/sqrt(3._wp)]*2._wp*pi
			else
				a(:,1)=[1._wp,0._wp]*2._wp*pi
				a(:,2)=[0._wp,1._wp]*2._wp*pi
			endif
		endif
		!!$OMP PARALLEL DO
		do ik=1,size(k,2)
			do i=1,D
				k(i,ik)=sum(a(i,:)*(map%get_idx(ik,[1:D])-1)/nk)
			enddo
			!if(sum(abs(k(:,ik)-pi))<1d-8) then
				!if(ipi==0) then
					!ipi=ik
				!else
					!stop "err"
				!endif
			!elseif(sum(abs(k(:,ik)))<1d-8) then
				!if(i0==0) then
					!i0=ik
				!else
					!stop "err"
				!endif
			!endif
		enddo
		!!$OMP END PARALLEL DO
		write(*,*)"stating get unique k point...."
		call get_uk(k,uk,uk2k)
		write(*,*)"finished get unique k point",size(k,2)/size(uk,2)
		contains
			subroutine into_sym(k)
				real(8) :: k(:)
				integer :: i,j
				real(8) :: tmp
				if(is_tri) then
					! 2-D triangular lattice
					k=rotate(abs(k(1:2)-[2._wp*pi,0._wp])-[2._wp*pi/3._wp,0._wp],pi/6._wp)
					k(1)=abs(k(1))
					k=abs(rotate(k,-pi/6._wp))*[-1._wp,1._wp]+[4._wp*pi/3._wp,0._wp]
				else
					! cubic lattice
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
					k=k+pi
				endif
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
					tmp(:,i)=k(:,i)
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
	end subroutine
	function rotate(r,theta) result(rt)
		! anticlock-wise rotation
		real(wp) :: r(2),theta
		real(wp) :: rt(2)
		rt=r*cos(theta)+[-r(2),r(1)]*sin(theta)
	end function
	complex(wp) function msum(G,not_conv) result(rt)
		complex(wp) :: G(:)
		logical, optional :: not_conv
		complex(wp) :: G1,G2,d
		real(wp) :: ipi_
		integer :: i,rg_(2),i1,i2
		!rt=sum(G)
		if(present(not_conv).and.(.not.is_normal_mat)) then
			rt=0.25_wp*(G(n)*img*omega(n)+G(1)*img*omega(1))
			!if(abs(G(n)*img*omega(n)-G(1)*img*omega(1))>1e-6_wp) then
				!write(*,*)"warn!!! msum may not be accurate",abs(G(n)*img*omega(n)-G(1)*img*omega(1))
			!endif
		else
			rt=0._wp
			!if(abs(G(n)*img*omega(n))+abs(G(1)*img*omega(1))>1e-6_wp) then
				!write(*,*)"warn!!! msum may not converge",abs(G(n)*img*omega(n))+abs(G(1)*img*omega(1))
			!endif
		endif
		rt=rt+sum(-Res*G)
		!ipi_=1._wp/pi
		!if(present(rg)) then
			!rg_=rg
		!else
			!rg_=[1,n]
		!endif
		!rt=0._wp
		!do i=rg_(1),rg_(2)-1
			!i1=floor((omega(i)-10000._wp*eps)*ipi_)+1
			!i2=floor((omega(i+1)+10000._wp*eps*merge(-1._wp,1._wp,i+1/=rg_(2)))*ipi_)
			!i1=i1+mod(abs(i1)+1,2)
			!i2=i2-mod(abs(i2)+1,2)
			!d=(G(i)-G(i+1))/(omega(i)-omega(i+1))
			!G1=G(i)+d*(i1*pi-omega(i))
			!G2=G(i)+d*(i2*pi-omega(i))
			!rt=rt+0.5_wp*(G1+G2)*((i2-i1)/2+1)
			!!write(*,*)i1,((i2-i1)/2+1)
		!enddo
	end function
end module
include "../lib/serde.f90"
include "../lib/fft.f90"
program main
	use M_DMFT
	implicit none
	type(gfs(n,D,nk)) :: phy
	complex(8) :: pDelta(n)
	type(gf(D,nk)) :: pdSE(n)
	integer :: j,i,jk,iT,iU,ik,iw,id_
	real(8) :: ff_(n),dTk,k1(D),k2(D)
	complex(8) :: r(nkD),rp(nkD)
	logical :: flag

	call omp_set_nested(.false.)
	call omp_set_dynamic(.false.)
	call mkl_set_dynamic(0)

	call mkl_set_num_threads(1)
	!call omp_set_num_threads(mkl_get_max_threads())
	call omp_set_num_threads(16)

	open(unit=10,File="../data/fk.dat")
	open(unit=20,file="../data/fk_init.dat")
	open(unit=30,File="../data/fk_q.dat")
	open(unit=40,File="../data/fk_phase.dat")
	open(unit=50,File="../data/fk_dG.dat")
	open(unit=60,File="../data/fk_check.dat")
	open(unit=70,File="../data/fk_exp.dat")
	call gen_k(phy%k,phy%uk,phy%uk2k)
	!call phy%export([60],[export_check])
	!stop
	write(40,"(A)")"U T kx ky Td kx1 ky1"
	call phy%set_Gfunction([set_ek])


	if(.not.is_real) then
		call phy%set_Gfunction([set_omega])
	endif
	!iomega=cmplx(0._wp,omega,kind=8)
	!write(*,*)sum(-abs(iomega)),msum(cmplx(-abs(iomega),kind=wp))
	!stop
	phy%Delta=cmplx(1d-1,1d-1)
	Tk=0.135_wp
	!Tk=0.25_wp
	do 
		!U=for_in([[20:4:-1]*1._wp,[10:1:-1]*0.3_wp],id=1)
		!U=for_in([1._wp],id=1)
		U=for_in([[20:6:-2]*1._wp,[5*2:1:-1]*0.5_wp],id=1)
		mu=0.5d0*U
		if(isnan(U)) exit
		if(next(xtol=1d-4,dx=-min(1d-2,Tk*0.3_wp),id=2)) then
		endif
		if(is_real) then
			call phy%set_Gfunction([set_omega])
		endif
		do i=1,n
			phy%dSE(i)%k=0._wp
		enddo
		write(70,"(A)")"U T B Sus Sus0 kx1 ky1"
		flag=.true.
		do 
			if(phy%B(phy%uk2k(iorder,1))<1._wp) then
				pDelta=phy%Delta
			endif
			Tk=max(Tk,1e-4_wp)
			beta=1._wp/Tk

			if(next(ytol=1d-3,dx=max(U/2._wp,1e-1_wp),id=1)) then
			endif
			do
				call phy%self_consistent(0.05d0,[300,0],1d-6)
				call phy%set_Gfunction([set_nc])
				!write(*,"(*(es12.4))")mu,phy%nc
				!exit
				if(next(mu,phy%nc-nc,id=1)) then
					exit
				endif
			enddo
			call phy%set_Gfunction([set_B,set_Sus_static])
			!write(*,*)iorder
			write(*,"('DMFT: ',*(es12.4))")U,Tk,phy%B(phy%uk2k(iorder,1)),phy%Sus(phy%uk2k(iorder,1)),phy%uk(:,iorder)/pi,phy%nc,mu
			if(next(Tk,1._wp-phy%B(phy%uk2k(iorder,1)),id=2)) then
				if(phy%B(phy%uk2k(iorder,1))<1d0) then
					write(70,"(*(e29.20e3))")U,Tk,phy%B(phy%uk2k(iorder,1)),phy%Sus(phy%uk2k(iorder,1)),Tk*msum(phy%Sus0(1:n)%k(phy%uk2k(iorder,1))),phy%uk(:,iorder)/pi
					call phy%export([10,30],[export_w,export_k])
					write(*,"(*(es12.4))")U,Tk,phy%B(phy%uk2k(iorder,1)),phy%uk(:,iorder)/pi,phy%nc
					write(40,"(e29.20e3$)")U,Tk,phy%uk(:,iorder)/pi
					exit
				else
					Tk=Tk+1d-4
				endif
			endif
		enddo
		write(70,"(x)")
		write(40,"(x)"); cycle
		!if(next(xtol=1d-4,dx=-min(1d-2,Tk*0.3_wp),id=1)) then
		!Tk=Tk+0.03_wp
		!Tk=1.7e-1
		flag=.true.
		if(next(xtol=1d-5,dx=-1d-4,id=4)) then
		endif
		if(next(xtol=1d-5,dx=-1d-3,id=3)) then
		endif
		if(next(xtol=1d-5,dx=-5d-3,id=2)) then
		endif
		id_=2
		do 
			Tk=max(Tk,1e-4_wp)
			beta=1._wp/Tk
			pdSE=phy%dSE
			pDelta=phy%Delta
			if(next(ytol=1d-3,dx=max(U/2._wp,1e-1_wp),id=1)) then
			endif
			do
				call phy%self_consistent(0.04d0,[300,300],1d-8)
				call phy%set_Gfunction([set_nc])
				!write(*,"(*(es12.4))")mu,phy%nc
				!exit
				if(next(mu,phy%nc-nc,id=1)) then
					exit
				endif
			enddo
			call phy%set_Gfunction([set_B,set_Sus_static])
			if(phy%B(phy%uk2k(iorder,1))>1d0) then
				phy%dSE=pdSE
				phy%Delta=pDelta
				if(flag) then
					Tk=Tk+1d-3
					cycle
				endif
			else
				write(70,"(*(e29.20e3))")U,Tk,phy%B(phy%uk2k(iorder,1)),phy%Sus(phy%uk2k(iorder,1)),phy%k(:,iorder)/pi
			endif
			flag=.false.
			if(phy%B(phy%uk2k(iorder,1))>0.95_wp) then
				id_=4
			elseif(phy%B(phy%uk2k(iorder,1))>0.9_wp) then
				id_=3
			endif
			if(next(Tk,1d0-phy%B(phy%uk2k(iorder,1)),id=id_)) then
				write(40,"(*(e29.20e3))")Tk,phy%uk(:,iorder)/pi
				exit
			endif
			write(*,"('DF: ',*(es12.4))")U,Tk,phy%B(phy%uk2k(iorder,1)),phy%Sus(phy%uk2k(iorder,1)),Tk*msum(phy%Sus0(1:n)%k(phy%uk2k(iorder,1))),phy%uk(:,iorder)/pi,phy%nc,mu
		enddo
		write(70,"(x)")
	enddo
end program
