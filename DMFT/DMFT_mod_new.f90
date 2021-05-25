#include "../../DGC/DGC_mod.f90"
module M_DMFT
	use global
	use M_DGC
	use lapack95, only : gesv, rsteqr, geev
	use M_const
	use M_solution
	use M_matrix
	use M_utility
	use mkl_service
	use M_serde
	use M_fft
	!use M_pade
	implicit none
	integer, parameter :: DMFT=1,dDMFT=2,dDMFT_simp=3,dDMFT_simp_nc=4,dDMFT_all=5,xx_save=1,xx_diff=2,xx_new=-1
	integer :: scheme=dDMFT
	integer :: which_gen_cfg=1,ipp,ip0,Cnv,Cnvc,pm(Nc,Nc),ica,info
	real(8) :: ba(D,D),bc(D,D),bs(D,D)
	real(8), allocatable :: Ta(:,:),Tc(:,:),Ts(:,:)
	complex(8) :: iomega(n)=0._wp
	integer, parameter :: out_w=20,out_k=2,out_wk=3,out_check=4,out_lattice=50,out_T=60,out_pattern=70,out_kc=80,out_r=90,inout_save=8,&
		omega_grid=1,dSEc=2,dVer=3,dSus0=4,SEc=6,dGc0=7,Gk=8,Gc=10,SEck=11,Gsk=12,Gs=13,lattice=14,Deltas=15,inc=16,Sus=17,dGc=18,dGcp=19,Delta=20,Gck_lp=21,db=22,Gs2c=23,Deltas2c=24,Ef=24,fe=25,V4c=26,dF4=27,Gck_l=28,mu=29,PMT=30,mat_omega=31,dF4d=32,incp=33,dGk=34,fsd=35,fs=36,prime=50,SC_DMFT=-1,SC_dDMFT=-2,SC_cp=-3
	integer :: niter(-5:-1)=0
	real(8) :: iter_err(-5:5)[*]=1._wp,iter_rate(-5:-1)=1._wp
	integer :: ipi,i0,iorder,iw0
	integer, allocatable :: idx(:),crt(:,:)
	type(t_serde(D)) :: map
	type gfs
		integer :: n=n,D=D,nk=nk,Ns=Ns,Nc=Nc
		real(wp) :: k(D,nk**D*Nc),ks(D,nk**D*Nc/Ns),kc(D,nk**D),Qs(D,Ns),Qc(D,Nc/Ns),Ef=0._wp,mu=0._wp
		real(wp) :: rs(D,Ns),rc(D,Nc),r(D,nk**D*Nc)
		real(wp) :: n_c(Nc),n_cimp(Nc),n_f(Nc),db(Nc),ek(nk**D*Nc),fsd,fs(nk**D*Nc)
		integer, allocatable :: uk2k(:,:),ukc2kc(:,:)
		complex(wp) :: Gk(n,nk**D*Nc),dGk(n,nk**D*Nc),dSEk(n,nk**D*Nc),ife,fe,Susr(nkD*Nc),Susc(nkD,Nc,Nc),Sus(nkD*Nc),Sus0(nkD*Nc)
#ifdef cluster
		complex(wp) :: Delta(n,Nc,Nc),Gc(n,Nc,Nc),SEc(n,Nc,Nc)
#else
		complex(wp) :: Delta(n,Nc),Gc(n,Nc),SEc(n,Nc)
#endif
		complex(wp) :: eks(nk**D*Nc/Ns,Ns,Ns)
		complex(wp) :: ekc(nk**D,Nc,Nc),dGc0(nk**D,n,Nc,Nc),dGc(nk**D,n,Nc,Nc),dSEc(nk**D,n,Nc,Nc),SEck(nk**D,n,Nc,Nc),Gck(nk**D,n,Nc,Nc),Gcl(n,Nc,Nc)!,Sus0(nk**D,n,Nc*Nc,Nc*dNs),dSus0(nk**D,n*Nc*Nc,Nc*dNs)
		complex(wp), allocatable :: dF4(:,:),dF4d(:,:,:,:),G4c(:,:),Sus0c(:,:,:,:),dSus0(:,:,:)
		logical :: conv,is_symbk=.true.
	contains
		procedure :: to_k
		procedure :: io
	end type
	complex(wp), allocatable :: x(:),x_(:)
	type(gfs) :: self[*]
#ifdef cluster
	complex(wp) :: LL(n,Nc,Nc**2,nkD),RR(n,Nc**2,Nc,nkD)
#else
	complex(wp) :: LL(n,Nc,Nc,nkD),RR(n,Nc,Nc,nkD)
#endif
contains
	subroutine coarray(flag)
		integer :: flag(:)
		integer :: ik,iw,i,j,iQ,l,i1,j1
		sync all
		do l=1,size(flag)
			select case(flag(l))
			case(lattice)
				if(this_image()/=1) then
					self%uk2k=self[1]%uk2k
					self%ukc2kc=self[1]%ukc2kc
				endif
			case(dF4d)
				do ica=1,num_images()
					if(ica/=this_image()) then
						do ik=this_image(),merge(nkD,size(self%ukc2kc,1),self%is_symbk),num_images()
							i=merge(ik,self%ukc2kc(ik,1),self%is_symbk)
#ifdef cluster
							self[ica]%dF4d(1:n,1:Nc**2,1:Nc**2,i:i)=self%dF4d(1:n,1:Nc**2,1:Nc**2,i:i)
#else
							self[ica]%dF4d(1:n,1:Nc,1:Nc,i:i)=self%dF4d(1:n,1:Nc,1:Nc,i:i)
#endif
						enddo
					endif
				enddo
			case(dF4d+prime)
				if(this_image()/=1) then
					self%dF4d=self[1]%dF4d
				endif
			case(Sus)
				do ica=1,num_images()
					if(ica/=this_image()) then
						do i=this_image(),merge(nkD,size(self%ukc2kc,1),self%is_symbk),num_images()
							do j=1,merge(1,self%ukc2kc(i,0),self%is_symbk)
								ik=merge(i,self%ukc2kc(i,j),self%is_symbk)
								do i1=1,Nc
									do j1=1,Nc
										self[ica]%Susc(ik,i1,j1)=self%Susc(ik,i1,j1)
									enddo
								enddo
							enddo
						enddo
					endif
				enddo
			case(dSEc)
				if(this_image()/=1) then
					self%dSEc=self[1]%dSEc
				endif
			case(dSEc+prime)
				do ica=1,num_images()
					if(ica/=this_image()) then
						do ik=1,nkD
							do i=1,Nc
								do j=1,Nc
									do iw=this_image(),n,num_images()
										!!self[1]%dSEc(1:nkD,iw,1:Nc,1:Nc)=self%dSEc(1:nkD,iw,1:Nc,1:Nc)
										self[ica]%dSEc(ik,iw,i,j)=self%dSEc(ik,iw,i,j)
									enddo
								enddo
							enddo
						enddo
					endif
				enddo
			case(Ef)
				if(this_image()/=1) then
					self%Ef=self[1]%Ef
				endif
			case(Delta)
				if(this_image()/=1) then
					self%Delta=self[1]%Delta
					self%is_symbk=self[1]%is_symbk
				endif
			case(mu)
				if(this_image()/=1) then
					self%mu=self[1]%mu
				endif
			case(:-1)
				if(this_image()/=1) then
					iter_err(abs(flag(l)))=iter_err(abs(flag(l)))[1]
				endif
			case default
				write(*,*)"coarray case not implimented"
				stop
			end select
		enddo
		sync all
	end subroutine
	function compute_(node) result(rt)
		!type(gfs) :: self[*]
		integer :: node
		character(:), allocatable :: rt
		integer, save :: nx,cnd_=0
		real(wp) :: U_(n,n),tf(D,D),dr(D),tmpD(D),rc_(D,Nc/Ns),rc(D),om!,b(D,D),k0(D),norm
		integer :: iord(n),i,j,m,ik,Nc1(Nc),iw,jw,iQ,l,i1,j1
#ifdef cluster
		integer :: Nc2(Nc*Nc)
		complex(wp) :: dSus0d(nkD,Nc**2,Nc**2)
#else
		integer :: Nc2(Nc)
		complex(wp) :: dSus0d(nkD,Nc,Nc)
#endif
		complex(wp) :: iA(n,Nc,Nc),iAs(Ns,Ns),iAc(Nc,Nc),iBc(Nc,Nc),iAc2(Nc**2,Nc**2),iBc2(Nc**2,Nc**2),cNc4(Nc**2,Nc**2),r1(nkD,Nc,Nc),r2(nkD,Nc,Nc),r(nkD),cNc(Nc),cn(n),cnp(n)
		!complex(wp), allocatable :: LL(:,:,:,:),RR(:,:,:,:)

		!write(*,*)"<",node,this_image()
		!sync all
		!if(all([dSus0,V4c,Ef,mu,Delta,dSEc,Sus,lattice,omega_grid,mat_omega,[-5:-1]]/=node).and.this_image()/=1) then
			!return
		!endif
		!write(*,*)">",node,this_image()
		select case(node)
		case(PMT)
			rt="PMT"; if(debug) return
		case(omega_grid)
			rt="omega_grid"; if(debug) return
			if(is_real) then
				!! real frequency
				!call set_grid(omega(n/2+1:),[reset,add_linear,add_linear,add_linear],from=[0d0,max(U/2d0-5d0,0d0),1.5d0],to=[bcut,min(U/2d0+5d0,bcut),3d0],n=[nadd])
				call set_grid(omega(n/2+1:),[reset,add_linear],from=[0d0],to=[sum(bcut*[-1d0,1d0])/2._wp],n=[n/2])
				omega(1:n/2)=-omega(n:n/2+1:-1)
				omega=omega+sum(bcut)/2._wp
				iomega=cmplx(omega,et,kind=8)
				iw0=0
				do i=1,n-1
					if(omega(i)*omega(i+1)<=0._wp) then
						iw0=merge(i,i+1,abs(omega(i))<abs(omega(i+1)))
						exit
					endif
				enddo
				if(iw0==0) then
					stop "error iw0"
				endif
			else
				iw0=n/2
				! matsubara frequency
				if(is_normal_mat) then
					call set_grid(omega(n/2+1:),[reset,add_linear],from=[0d0],to=[nadd(1)*2d0*pi],n=[nadd(1)])
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
					call rsteqr(omega,Res(1:n-1),U_,info=info)
					omega(1:n/2)=(omega(1:n/2)-omega(n:n/2+1:-1))*0.5_wp
					do i=1,n/2
						U_(1:n,i)=sign(1._wp,U_(1:n,i))*0.5_wp*(abs(U_(1:n,i))+abs(U_(1:n,n+1-i)))
						U_(1:n,n+1-i)=U_(1:n,i)*(-1._wp)**([1:n])
					enddo
					if(info/=0) then
						write(*,*)"mat_omega error",info
						stop
					endif

					Res(n/2:1:-1)=-0.25_wp*(abs((U_(1,1:n/2)))/omega(1:n/2))**2
					omega(1:n/2)=1._wp/omega(n/2:1:-1)
					Res(n/2+1:n)=Res(n/2:1:-1)
					omega(n/2+1:n)=-omega(n/2:1:-1)

				endif
			endif
			do i=1,n-1
				if(omega(i)>=omega(i+1)) then
					write(*,*)"omega error"
					stop
				endif
			enddo
		case(mat_omega)
			rt="mat_omega"; if(debug) return
			beta=1._wp/Tk
			if(.not.is_real) then
				if(sum(abs(iomega))>1e-8_wp) then
					self%is_symbk=.false.
					do iw=1,n
#ifdef cluster
						if(any(abs(diag(self%Delta(iw,1:Nc,1:Nc))-self%Delta(iw,1,1))>1e-6_wp)) then
#else
						if(any(abs(self%Delta(iw,1:Nc)-self%Delta(iw,1))>1e-6_wp)) then
#endif
							self%is_symbk=.true.
							exit
						endif
					enddo
					do i=1,Nc
						do j=1,Nc
#ifdef cluster
							self%Delta(1:n,i,j)=fit(self%Delta(1:n,i,j),imag(iomega),omega*Tk)
#endif
							do ik=1,nkD
								self%dSEc(ik,1:n,i,j)=fit(self%dSEc(ik,1:n,i,j),imag(iomega),omega*Tk)
							enddo
							!write(*,*)"fitting"
						enddo
#ifndef cluster
						self%Delta(1:n,i)=fit(self%Delta(1:n,i),imag(iomega),omega*Tk)
#endif
					enddo
				endif
				iomega=cmplx(0d0,omega*Tk,kind=8)
			endif

		case(lattice)
			rt="lattice"; if(debug) return
			! basic checking begin
			!if(is_CDMFT.and.is_dual_DCAmix) then
				!write(*,"(A)")"is_CDMFT and is_dual_DCAmix can not be true at the same time"
				!stop
			!endif
			if(is_CDMFT.and.dNs/=Nc) then
				write(*,"(A)")"in CDMFT case, dNs must be equal to Nc"
				stop
			endif
			! basic checking end
			map%shap=[(nk,i=1,D)]
			ba=recip(a)
			bs=recip(scell)
			bc=recip(clust)
			call get_wigner_seitz(ba,Ta)
			call get_wigner_seitz(bs,Ts)
			call get_wigner_seitz(bc,Tc)
			call gen_grid(bs(:,1),bs(:,2),[0._wp,0._wp],reshape([[0._wp,0._wp],ba(:2,1),ba(:2,1)+ba(:2,2),ba(:2,2)],[2,4]),self%Qs)
			call gen_grid(bc(:,1),bc(:,2),[0._wp,0._wp],reshape([[0._wp,0._wp],bs(:2,1),bs(:2,1)+bs(:2,2),bs(:2,2)],[2,4]),self%Qc)

			m=0
			tf=reshape([scell(1,1),scell(1,2),scell(2,1),scell(2,2)],[2,2])/(2._wp*pi)
			ipp=0
			ip0=0
			do ik=1,nkD
				do i=1,D
					self%kc(i,ik)=sum(bc(i,:)*(map%get_idx(ik,[1:D])-1)/nk)
					self%r(i,ik)=sum(clust(i,:)*(map%get_idx(ik,[1:D])-1))
				enddo
				if(.not.is_in(self%kc(:,ik),Tc,dr)) then
					self%kc(:,ik)=self%kc(:,ik)+dr
				endif
				do i=1,Nc/Ns
					self%ks(:,ik+nkD*(i-1))=self%kc(:,ik)+self%Qc(:,i)
					if(.not.is_in(self%ks(:,ik+nkD*(i-1)),Ts,dr)) then
						self%ks(:,ik+nkD*(i-1))=self%ks(:,ik+nkD*(i-1))+dr
					endif
					do j=1,Ns
						self%k(:,ik+nkD*(i-1)+nkD*Nc/Ns*(j-1))=self%ks(:,ik+nkD*(i-1))+self%Qs(:,j)
						if(.not.is_in(self%k(:,ik+nkD*(i-1)+nkD*Nc/Ns*(j-1)),Ta,dr)) then
							self%k(:,ik+nkD*(i-1)+nkD*Nc/Ns*(j-1))=self%k(:,ik+nkD*(i-1)+nkD*Nc/Ns*(j-1))+dr
						endif
						tmpD=self%k(:,ik+nkD*(i-1)+nkD*Nc/Ns*(j-1))-Ta(:,1)
						if(.not.is_in(tmpD,Ta,dr)) then
							tmpD=tmpD+dr
						endif
						if(sum(abs(tmpD))<1e-6) then
							if(ipp==0) then
								ipp=ik+nkD*(i-1)+nkD*Nc/Ns*(j-1)
							else
								write(*,*)"ipp error"
								stop
							endif
						endif

						tmpD=self%k(:,ik+nkD*(i-1)+nkD*Nc/Ns*(j-1))-(Ta(:,1)+Ta(:,2))*0.5_wp
						if(.not.is_in(tmpD,Ta,dr)) then
							tmpD=tmpD+dr
						endif
						if(sum(abs(tmpD))<1e-6) then
							if(ip0==0) then
								ip0=ik+nkD*(i-1)+nkD*Nc/Ns*(j-1)
							else
								write(*,*)"ip0 error"
								stop
							endif
						endif
					enddo
				enddo
			enddo

			call gen_grid(a(:,1),a(:,2),rc0,reshape([[0._wp,0._wp],scell(:2,1),scell(:2,1)+scell(:2,2),scell(:2,2)],[2,4]),self%rs)
			call gen_grid(scell(:,1),scell(:,2),[0._wp,0._wp],reshape([[0._wp,0._wp],clust(:2,1),clust(:2,1)+clust(:2,2),clust(:2,2)],[2,4]),rc_)
			do i=1,Nc
				self%rc(:,i)=self%rs(:,mod(i-1,Ns)+1)+rc_(:,(i-1)/Ns+1)
				do j=1,D
					self%r(j,(i-1)*nkD+1:(i)*nkD)=self%r(j,(i-1)*nkD+1:(i)*nkD)+self%rc(j,i)
				enddo
				!if(.not.is_in(self%rc(:,i),reshape([[0._wp,0._wp],clust(:2,1),clust(:2,1)+clust(:2,2),clust(:2,2)],[2,4]),dr)) then
					!self%rc(:,i)=self%rc(:,i)+dr
				!endif
			enddo
			!tf=reshape([bs(1,1),bs(1,2),bs(2,1),bs(2,2)],[2,2])/(2._wp*pi*nk)
			do i=1,Nc
				do j=1,Nc
					self%ekc(:,i,j)=0._wp
					dr=self%rc(:,i)-self%rc(:,j)
					do ik=1,nkD
						rc=matmul(clust(:,:),(map%get_idx(ik,[1:D])-1))
						if(any(abs(sum((spread(dr+rc,2,9)+reshape([matmul(clust(:,:),[0,0]*nk),matmul(clust(:,:),[1,0]*nk),matmul(clust(:,:),[0,1]*nk),matmul(clust(:,:),[1,1]*nk),matmul(clust(:,:),[-1,0]*nk),matmul(clust(:,:),[0,-1]*nk),matmul(clust(:,:),[-1,-1]*nk),matmul(clust(:,:),[-1,1]*nk),matmul(clust(:,:),[1,-1]*nk)],[D,9]))**2,dim=1)-1._wp)<1e-8_wp)) then
							self%ekc(:,i,j)=self%ekc(:,i,j)-t*exp(-img*matmul(rc,self%kc(:,:)))
						endif
					enddo
				enddo
			enddo
			do ik=1,nkD*Nc/Ns
				do i=1,Ns
					do j=1,Ns
						self%eks(ik,i,j)=sum(self%ekc(mod(ik-1,nkD)+1,i::Ns,j)*exp(-img*matmul(self%ks(1:D,ik),rc_(1:D,:))))
					enddo
				enddo
			enddo
			do ik=1,nkD*Nc
				self%ek(ik)=sum(self%eks(mod(ik-1,nkD*Nc/Ns)+1,1:Ns,1)*exp(-img*matmul(self%k(1:D,ik),self%rs(1:D,1:Ns)-spread(self%rs(1:D,1),D,Ns))))
			enddo
			!do ik=1,nkD*Nc
				!write(*,*)abs(self%ek(ik)-sum(self%ekc(mod(ik-1,nkD)+1,1:Nc,1)*exp(-img*matmul(self%k(1:D,ik),self%rc(1:D,1:Nc)-spread(self%rc(1:D,1),D,Nc)))))
			!enddo
			!stop
			call get_uk(self%k(1:D,1:nkD*Nc),self%uk2k)
			if(Nc>1) then
				call get_uk(self%kc(1:D,1:nkD),self%ukc2kc)
			endif
			do ik=1,size(self%ukc2kc,1)
				do i=2,self%ukc2kc(ik,0)
					Nc1=crt(1:Nc,self%ukc2kc(ik,i+self%ukc2kc(ik,0)))
					if(sum(abs(self%ekc(self%ukc2kc(ik,1),1:Nc,1:Nc)-self%ekc(self%ukc2kc(ik,i),Nc1,Nc1)))>1e-8_wp) then
						write(*,*)"ukc2kc error!"
						stop
					endif
				enddo
			enddo

			if(&
				sum([(sum(abs(&
				self%to_k(cmplx(self%eks(ik,1:Ns,1:Ns),kind=wp),self%ks(1:D,ik))&
			+diag(&
				merge(2._wp*t*(cos(self%k(1,ik:nkD*Nc:nkD*Nc/Ns))+cos(self%k(2,ik:nkD*Nc:nkD*Nc/Ns))),&
				2d0*t*(cos(self%k(1,ik:nkD*Nc:nkD*Nc/Ns))+cos(0.5_wp*(self%k(1,ik:nkD*Nc:nkD*Nc/Ns)+sqrt(3._wp)*self%k(2,ik:nkD*Nc:nkD*Nc/Ns)))+cos(0.5_wp*(-self%k(1,ik:nkD*Nc:nkD*Nc/Ns)+sqrt(3._wp)*self%k(2,ik:nkD*Nc:nkD*Nc/Ns)))),sum(abs(a(:,1)*a(:,2)))<1e-9)&
			))),ik=1,nkD*Nc/Ns)])>1e-7&
			.or.&
			sum(abs(&
				self%ek(1:nkD*Nc)&
				+merge(2._wp*t*(cos(self%k(1,1:nkD*Nc))+cos(self%k(2,1:nkD*Nc))),&
				2d0*t*(cos(self%k(1,1:nkD*Nc))+cos(0.5_wp*(self%k(1,1:nkD*Nc)+sqrt(3._wp)*self%k(2,1:nkD*Nc)))+cos(0.5_wp*(-self%k(1,1:nkD*Nc)+sqrt(3._wp)*self%k(2,1:nkD*Nc)))),sum(abs(a(:,1)*a(:,2)))<1e-9)))>1e-7&
				) then
				write(*,*)"ek error!"
				stop
			endif
			call get_kline(self%k(1:D,1:nkD*Ns),[0._wp,0._wp,(ba(:,1)+ba(:,2))/3._wp,ba(:,1)/2._wp],[(ba(:,1)+ba(:,2))/3._wp,ba(:,1)/2._wp,0._wp,0._wp],idx)
			!call get_kline(self%k,[(ba(:,1)+ba(:,2))/3._wp],[ba(:,1)/2._wp],idx)


			pm=0
			do i=1,Nc
				do j=1,Nc
					if(i<=Ns) then
						pm(i,j)=j
					else
						do ik=1,Nc
							tmpD=((self%rc(:,j)-self%rc(:,i))-(self%rc(:,ik)-self%rc(:,mod(i-1,Ns)+1)))
							if(.not.is_in(tmpD,reshape([[0._wp,0._wp],clust(:2,1),clust(:2,1)+clust(:2,2),clust(:2,2)],[2,4]),dr)) then
								tmpD=tmpD+dr
							endif
							if(sum(abs(tmpD))<1d-6) then
								if(pm(i,j)/=0) then
									write(*,*)"pm error"
									stop
								else
									pm(i,j)=ik
								endif
							endif
						enddo
					endif
				enddo
			enddo
			call coarray([node])
		case(db)
			rt="db"; if(debug) return
			call solver(self%Delta,self%mu,self%Ef,self%Gc,self%ife,db=self%db)
		case(V4c)
			rt="V4c"; if(debug) return
			if(.not.allocated(self%G4c)) then
#ifdef cluster
				allocate(self%G4c(n*Nc*Nc,n*Nc*Ns))
#else
				allocate(self%G4c(n*Nc,n*Nc))
#endif
			endif
			!call coarray([Delta])
			call solver(self%Delta,self%mu,self%Ef,self%Gc,self%ife,G4c=self%G4c)
			!write(*,*)sum(abs(self%G4c)),sum(abs(self%Gc))
#ifdef cluster
			!$OMP PARALLEL DO
			do iw=1,n
				iA(iw,1:Nc,1:Nc)=self%Gc(iw,1:Nc,1:Nc)
				call mat_inv(iA(iw,1:Nc,1:Nc))
			enddo
			!$OMP END PARALLEL DO
#endif
			!$OMP PARALLEL DO PRIVATE(cNc4)
			do iw=1,n
				do jw=1,n
#ifdef cluster
					if(iw==jw) then
						self%G4c(iw:n*Nc**2:n,iw:n*Nc*Ns:n)=0._wp
						cycle
					endif
					cNc4(1:Nc**2,1:Nc*Ns)=self%G4c(iw:n*Nc**2:n,jw:n*Nc*Ns:n)
					do i=1,Nc**2
						do j=1,Nc*Ns
							self%G4c(iw+(i-1)*n,jw+(j-1)*n)=sum(matmul([(iA(iw,mod(i-1,Nc)+1,:)*iA(iw,ik,(i-1)/Nc+1),ik=1,Nc)],cNc4(1:Nc**2,1:Nc*Ns))*[(iA(jw,mod(j-1,Nc)+1,:)*iA(jw,ik,(j-1)/Nc+1),ik=1,Ns)])&
								-iA(iw,mod(i-1,Nc)+1,(i-1)/Nc+1)*iA(jw,mod(j-1,Nc)+1,(j-1)/Nc+1)
						enddo
					enddo
#else
					if(iw==jw) then
						self%G4c(iw:n*Nc:n,iw:n*Nc:n)=0._wp
						cycle
					endif
					do i=1,Nc
						do j=1,Nc
							self%G4c(iw+(i-1)*n,jw+(j-1)*n)=self%G4c(iw+(i-1)*n,jw+(j-1)*n)/(self%Gc(iw,i)*self%Gc(jw,j))**2-1._wp/(self%Gc(iw,i)*self%Gc(jw,j))
							if(i/=j) then
								if(abs(self%G4c(iw+(i-1)*n,jw+(j-1)*n))>5e-6_wp) then
									write(*,*)"warning having inter site in non cluster case",abs(self%G4c(iw+(i-1)*n,jw+(j-1)*n))
								endif
							endif
						enddo
					enddo
					!do i=2,Nc
						!self%G4c(iw+(1-1)*n,jw+(1-1)*n)=self%G4c(iw+(1-1)*n,jw+(1-1)*n)+self%G4c(iw+(i-1)*n,jw+(i-1)*n)
					!enddo
					!self%G4c(iw+(1-1)*n,jw+(1-1)*n)=self%G4c(iw+(1-1)*n,jw+(1-1)*n)/Nc
					!do i=2,Nc
						!self%G4c(iw+(i-1)*n,jw+(i-1)*n)=self%G4c(iw+(1-1)*n,jw+(1-1)*n)
					!enddo
#endif
				enddo
			enddo
			!$OMP END PARALLEL DO
			!write(*,*)sum(abs(self%G4c))/Nc,sum(abs(self%Gc))/Nc
			!stop
		case(dGc0)
			rt="dGk0"; if(debug) return
			!$OMP PARALLEL DO PRIVATE(iAc,iBc)
			do iw=1,n
				do ik=1,nkD
					if(is_CDMFT.or.(.not.is_dual_pureDCA)) then
						iBc(1:Nc,1:Nc)=self%ekc(ik,1:Nc,1:Nc)
					elseif(is_dual_pureDCA) then
						do i=1,Nc
							do j=1,Nc
								iBc(i,j)=self%ekc(ik,i,j)*exp(-img*sum(self%kc(:,ik)*(self%rc(:,i)-self%rc(:,j))))
							enddo
						enddo
					endif
#ifdef cluster
					iBc=matmul(self%Delta(iw,1:Nc,1:Nc)-iBc,self%Gc(iw,1:Nc,1:Nc))
#else
					do i=1,Nc 
						do j=1,Nc
							iBc(i,j)=(merge(self%Delta(iw,i),cmplx(0._wp,kind=wp),i==j)-iBc(i,j))*self%Gc(iw,j)
						enddo
					enddo
#endif
					iAc=diag(1._wp,Nc)+iBc
					call mat_inv(iAc)
					iAc=matmul(iAc,iBc)
#ifdef cluster
					self%dGc0(ik,iw,1:Nc,1:Nc)=-matmul(self%Gc(iw,1:Nc,1:Nc),iAc)
#else
					do j=1,Nc
						self%dGc0(ik,iw,1:Nc,j)=-self%Gc(iw,1:Nc)*iAc(1:Nc,j)
					enddo
#endif
				enddo
			enddo
			!$OMP END PARALLEL DO
		case(dGcp)
			rt="dGcp"; if(debug) return
			!$OMP PARALLEL DO
			do iw=1,n
				self%dGc(1:nkD,iw,1:Nc,1:Nc)=self%dGc0(1:nkD,iw,1:Nc,1:Nc)
			enddo
			!$OMP END PARALLEL DO
		case(dGc)
			rt="dGc"; if(debug) return
			!$OMP PARALLEL DO PRIVATE(iAc)
			do iw=1,n
				do ik=1,nkD
					iAc(1:Nc,1:Nc)=diag(1._wp,Nc)-matmul(self%dGc0(ik,iw,1:Nc,1:Nc),self%dSEc(ik,iw,1:Nc,1:Nc))
					call mat_inv(iAc)
					self%dGc(ik,iw,1:Nc,1:Nc)=matmul(iAc,self%dGc0(ik,iw,1:Nc,1:Nc))
				enddo
			enddo
			!$OMP END PARALLEL DO
			!call self%io([out_k],[out_k])
			!if(sum(abs(self%dSEc))>1e-2_wp) then
				!stop
			!endif
		case(dSus0)
			rt="dSus0"; if(debug) return
			if(.not.allocated(self%dSus0)) then
#ifdef cluster
				allocate(self%dSus0(nk**D,n*Nc*Nc,Nc*dNs))
#else
				allocate(self%dSus0(nk**D,n*Nc,Nc))
#endif
			endif
			!call coarray([dGc])
			!$OMP PARALLEL DO PRIVATE(r1,r2,r)
			do iw=1,n
				do i=1,Nc
					do j=1,Nc
						call mfft(self%dGc(1:nkD,iw,i,j),1,map%shap,r1(1:nkD,i,j),1._wp/product(map%shap))
						call mfft(self%dGc(1:nkD,iw,i,j),-1,map%shap,r2(1:nkD,i,j),1._wp/product(map%shap))
					enddo
				enddo
#ifdef cluster
				do j=1,Nc*dNs
					do i=1,Nc*Nc
						r=-r2(1:nkD,(i-1)/Nc+1,mod(j-1,Nc)+1)*r1(1:nkD,(j-1)/Nc+1,mod(i-1,Nc)+1)
						call mfft(r(1:nkD),-1,map%shap,self%dSus0(1:nkD,iw+(i-1)*n,j))
					enddo
				enddo
#else
				do j=1,Nc
					do i=1,Nc
						r=-r2(1:nkD,i,j)*r1(1:nkD,j,i)
						call mfft(r(1:nkD),-1,map%shap,self%dSus0(1:nkD,iw+(i-1)*n,j))
					enddo
				enddo
#endif
			enddo
			!$OMP END PARALLEL DO
		case(Gck_l,Gck_lp)
			rt="Gck_l"; if(debug) return
			!$OMP PARALLEL DO PRIVATE(iAc,iBc)
			do iw=1,n
				self%Gcl(iw,1:Nc,1:Nc)=0._wp
				do ik=1,nkD
					if(is_CDMFT.or.(.not.is_dual_pureDCA)) then
						iAc(1:Nc,1:Nc)=-self%ekc(ik,1:Nc,1:Nc)
					elseif(is_dual_pureDCA) then
						do i=1,Nc
							do j=1,Nc
								iAc(i,j)=-self%ekc(ik,i,j)*exp(-img*sum(self%kc(:,ik)*(self%rc(:,i)-self%rc(:,j))))
							enddo
						enddo
					endif
					select case(lse)
					case(0)
#ifdef cluster
						iBc=self%Gc(iw,1:Nc,1:Nc)
						call mat_inv(iBc)
						iBc=matmul(matmul(iBc,self%dGc(ik,iw,1:Nc,1:Nc)),iBc)
						iAc=self%Delta(iw,1:Nc,1:Nc)+iAc
#else
						do i=1,Nc
							do j=1,Nc
								iBc(i,j)=self%dGc(ik,iw,i,j)/(self%Gc(iw,i)*self%Gc(iw,j))
							enddo
							iAc(i,i)=self%Delta(iw,i)+iAc(i,i)
						enddo
#endif
						call mat_inv(iAc)
						self%Gck(ik,iw,1:Nc,1:Nc)=iAc+matmul(matmul(iAc,iBc),iAc)
					case(1:2)
						iAc=diag(iomega(iw)+self%mu,Nc)-self%SEck(ik,iw,1:Nc,1:Nc)+iAc
						call mat_inv(iAc)
						self%Gck(ik,iw,1:Nc,1:Nc)=iAc
					end select
					self%Gcl(iw,1:Nc,1:Nc)=self%Gcl(iw,1:Nc,1:Nc)+self%Gck(ik,iw,1:Nc,1:Nc)
				enddo
				self%Gcl(iw,1:Nc,1:Nc)=self%Gcl(iw,1:Nc,1:Nc)/nkD
			enddo
			!$OMP END PARALLEL DO
		case(SEck)
			rt="SEck"; if(debug) return
			!$OMP PARALLEL DO PRIVATE(iAc)
			do iw=1,n
				do ik=1,nkD
					select case(lse)
					case(0)
						iAc=self%Gck(ik,iw,1:Nc,1:Nc)
						call mat_inv(iAc)
						self%SEck(ik,iw,1:Nc,1:Nc)=diag(iomega(iw)+self%mu,Nc)-iAc
						if(is_CDMFT.or.(.not.is_dual_pureDCA)) then
							self%SEck(ik,iw,1:Nc,1:Nc)=self%SEck(ik,iw,1:Nc,1:Nc)-self%ekc(ik,1:Nc,1:Nc)
						elseif(is_dual_pureDCA) then
							do i=1,Nc
								do j=1,Nc
									self%SEck(ik,iw,1:Nc,1:Nc)=self%SEck(ik,iw,1:Nc,1:Nc)-self%ekc(ik,i,j)*exp(-img*sum(self%kc(:,ik)*(self%rc(:,i)-self%rc(:,j))))
								enddo
							enddo
						endif
					case(1)
#ifdef cluster
						iAc=diag(1._wp,Nc)+matmul(self%Gc(iw,1:Nc,1:Nc),self%dSEc(ik,iw,1:Nc,1:Nc))
						call mat_inv(iAc)
						self%SEck(ik,iw,1:Nc,1:Nc)=self%SEc(iw,1:Nc,1:Nc)+matmul(self%dSEc(ik,iw,1:Nc,1:Nc),iAc)
#else
						if(Nc>1) then
							write(*,*)"SEck not implimented for noncluster"
							stop
						endif
						iAc=diag(1._wp,Nc)+matmul(diag(self%Gc(iw,1:Nc)),self%dSEc(ik,iw,1:Nc,1:Nc))
						call mat_inv(iAc)
						self%SEck(ik,iw,1:Nc,1:Nc)=diag(self%SEc(iw,1:Nc))+matmul(self%dSEc(ik,iw,1:Nc,1:Nc),iAc)
#endif
					case(2)
#ifdef cluster
						self%SEck(ik,iw,1:Nc,1:Nc)=self%SEc(iw,1:Nc,1:Nc)+self%dSEc(ik,iw,1:Nc,1:Nc)
#else
						self%SEck(ik,iw,1:Nc,1:Nc)=diag(self%SEc(iw,1:Nc))+self%dSEc(ik,iw,1:Nc,1:Nc)
#endif
					end select
				enddo
			enddo
			!$OMP END PARALLEL DO
			!write(*,*)U*n_f,self%SEck(1,1,1,1)
		case(SEc)
			rt="SEc"; if(debug) return
			!$OMP PARALLEL DO PRIVATE(iAc)
			do iw=1,n
#ifdef cluster
				iAc=self%Gc(iw,1:Nc,1:Nc)
				call mat_inv(iAc)
				self%SEc(iw,1:Nc,1:Nc)=diag(iomega(iw)+self%mu,Nc)-self%Delta(iw,1:Nc,1:Nc)-iAc
#else
				self%SEc(iw,1:Nc)=iomega(iw)+self%mu-self%Delta(iw,1:Nc)-1._wp/(self%Gc(iw,1:Nc))
#endif
			enddo
			!$OMP END PARALLEL DO
		case(Gk)
			rt="Gk"; if(debug) return
			do iw=1,n
				do ik=1,nkD*Nc
					self%Gk(iw,ik)=0._wp
					do i=1,Nc
						do j=1,Nc
							if(is_CDMFT.or.(.not.is_dual_pureDCA)) then
								self%Gk(iw,ik)=self%Gk(iw,ik)+self%Gck(mod(ik-1,nkD)+1,iw,i,j)*exp(-img*sum(self%k(:,ik)*(self%rc(:,i)-self%rc(:,j))))
							elseif(is_dual_pureDCA) then
								self%Gk(iw,ik)=self%Gk(iw,ik)+self%Gck(mod(ik-1,nkD)+1,iw,i,j)*exp(-img*sum((self%Qc(:,mod((ik-1)/nkD,Nc/Ns)+1)+self%Qs(:,(ik-1)/nkD/(Nc/Ns)+1))*(self%rs(:,i)-self%rs(:,j))))
							endif
						enddo
					enddo
					self%Gk(iw,ik)=self%Gk(iw,ik)/Nc
				enddo
			enddo
		case(fs)
			rt="fs"; if(debug) return
			self%fs=0._wp
			do ik=1,nkD*Nc
				self%fs(ik)=self%fs(ik)-real(sum(-Res(1:n)*self%Gk(1:n,ik)*exp(-0.5_wp*iomega(1:n)/Tk),abs(Res(1:n)+1._wp)<1e-5_wp))/pi
			enddo
		case(dGk)
			rt="dGk"; if(debug) return
			do iw=1,n
				do ik=1,nkD*Nc
					self%dGk(iw,ik)=0._wp
					self%dSEk(iw,ik)=0._wp
					do i=1,Nc
						do j=1,Nc
							if(is_CDMFT.or.(.not.is_dual_pureDCA)) then
								self%dGk(iw,ik)=self%dGk(iw,ik)+self%dGc(mod(ik-1,nkD)+1,iw,i,j)*exp(-img*sum(self%k(:,ik)*(self%rc(:,i)-self%rc(:,j))))
								self%dSEk(iw,ik)=self%dSEk(iw,ik)+self%dSEc(mod(ik-1,nkD)+1,iw,i,j)*exp(-img*sum(self%k(:,ik)*(self%rc(:,i)-self%rc(:,j))))
							elseif(is_dual_pureDCA) then
								self%dGk(iw,ik)=self%dGk(iw,ik)+self%dGc(mod(ik-1,nkD)+1,iw,i,j)*exp(-img*sum((self%Qc(:,mod((ik-1)/nkD,Nc/Ns)+1)+self%Qs(:,(ik-1)/nkD/(Nc/Ns)+1))*(self%rs(:,i)-self%rs(:,j))))
								self%dSEk(iw,ik)=self%dSEk(iw,ik)+self%dSEc(mod(ik-1,nkD)+1,iw,i,j)*exp(-img*sum((self%Qc(:,mod((ik-1)/nkD,Nc/Ns)+1)+self%Qs(:,(ik-1)/nkD/(Nc/Ns)+1))*(self%rs(:,i)-self%rs(:,j))))
							endif
						enddo
					enddo
					self%dGk(iw,ik)=self%dGk(iw,ik)/Nc
					self%dSEk(iw,ik)=self%dSEk(iw,ik)/Nc
				enddo
			enddo
		case(inc,incp)
			rt="inc"; if(debug) return
			if(debug) write(*,*)"inc, image: ",this_image()
			do i=1,Nc
				self%n_c(i)=Tk*real(msum(self%Gcl(1:n,i,i),not_conv=underscore))
#ifdef cluster
				self%n_cimp(i)=Tk*real(msum(self%Gc(1:n,i,i),not_conv=underscore))
#else
				self%n_cimp(i)=Tk*real(msum(self%Gc(1:n,i),not_conv=underscore))
#endif
			enddo
			!stop
		case(fsd)
			rt="fsd"; if(debug) return
			self%fsd=0._wp
			do i=1,Nc
				self%fsd=self%fsd-real(sum(-Res(1:n)*self%Gcl(1:n,i,i)*exp(-0.5_wp*iomega(1:n)/Tk),abs(Res(1:n)+1._wp)<1e-5_wp))/pi
			enddo
			self%fsd=self%fsd/Nc
		case(fe)
			rt="fe"; if(debug) return
			self%fe=nan
			!!$OMP PARALLEL DO PRIVATE(iAc,cNc)
			!do iw=1,n
				!cn(iw)=0._wp
				!do iQ=1,Nc/Ns
					!do ik=1,nkD
						!do i=1,Ns
							!do j=1,Ns
								!!iAs(i,j)=sum(self%Gc(iw,i:Nc:Ns,j)*exp(-img*matmul(self%Qc(:,iQ),rc_(:,:))))
								!if(is_CDMFT) then
									!iAs(i,j)=-self%eks((iQ-1)*nkD+ik,i,j)-self%SEs(iQ,iw,i,j)
								!else
									!iAs(i,j)=-self%eks((iQ-1)*nkD+ik,i,j)-self%SEs(iQ,iw,i,j)*exp(img*sum(self%kc(:,ik)*(self%rs(:,i)-self%rs(:,j))))
								!endif
							!enddo
							!iAs(i,i)=iAs(i,i)+iomega(iw)+self%mu
						!enddo
						!!iAs(1:Ns,1:Ns)=matmul(diag(iomega(iw)+self%mu,Ns)-self%eks(self%ikc(iQ,ik),1:Ns,1:Ns)-self%SEs(iQ,iw,1:Ns,1:Ns),iAs)
						!!iAs(1:Ns,1:Ns)=matmul(iBs,iAs)
						!iAs(1:Ns,1:Ns)=matmul(iAs,self%Gs(iQ,iw,1:Ns,1:Ns))
						!call geev(iAs(1:Ns,1:Ns),cNs)
						!cn(iw)=cn(iw)+sum(log(cNs))
					!enddo
				!enddo
			!enddo
			!!$OMP END PARALLEL DO
			!self%fe=self%ife-Tk*msum(cn,not_conv=underscore)/(Nc*nkD)
		case(Sus)
			rt="Sus"; if(debug) return
			if(.not.allocated(self%Sus0c)) then
				allocate(self%Sus0c(nk**D,n,Nc,Nc))
			endif
#ifdef cluster
			if(.not.allocated(self%dF4)) then
				allocate(self%dF4(n*Nc*Nc,n*Nc*dNs))
			endif
			if(.not.allocated(self%G4c)) then
				allocate(self%G4c(n*Nc*Nc,n*Nc*Ns))
			endif
			if(.not.allocated(self%dSus0)) then
				allocate(self%dSus0(nk**D,n*Nc*Nc,Nc*dNs))
			endif
			!allocate(LL(n,Nc,Nc**2,nkD),RR(n,Nc**2,Nc,nkD))
#else
			if(.not.allocated(self%dF4)) then
				allocate(self%dF4(n*Nc,n*Nc))
			endif
			if(.not.allocated(self%G4c)) then
				allocate(self%G4c(n*Nc,n*Nc))
			endif
			if(.not.allocated(self%dSus0)) then
				allocate(self%dSus0(nk**D,n*Nc,Nc))
			endif
			!allocate(LL(n,Nc,Nc,nkD),RR(n,Nc,Nc,nkD))
	
#endif
			!call coarray([Gck_l,Delta,Gc,dGc])
			!$OMP PARALLEL DO PRIVATE(iAc,iBc,r1,r2,r)
			do iw=1,n
				do ik=1,nkD
					if(is_CDMFT.or.(.not.is_dual_pureDCA)) then
						iAc=-self%ekc(ik,1:Nc,1:Nc)
					elseif(is_dual_pureDCA) then
						do i=1,Nc
							do j=1,Nc
								iAc(i,j)=-self%ekc(ik,i,j)*exp(-img*sum(self%kc(:,ik)*(self%rc(:,i)-self%rc(:,j))))
							enddo
						enddo
					endif
#ifdef cluster
					iBc=self%Gc(iw,1:Nc,1:Nc)
					call mat_inv(iBc)
					iAc=self%Delta(iw,1:Nc,1:Nc)+iAc
#else
					iBc=diag(1._wp/self%Gc(iw,1:Nc))
					iAc=diag(self%Delta(iw,1:Nc))+iAc
#endif
					call mat_inv(iAc)
					r1(ik,1:Nc,1:Nc)=matmul(matmul(iAc,iBc),self%dGc(ik,iw,1:Nc,1:Nc))
					r2(ik,1:Nc,1:Nc)=matmul(self%dGc(ik,iw,1:Nc,1:Nc),matmul(iBc,iAc))
				enddo

				do i=1,Nc
					do j=1,Nc
						r=r1(1:nkD,i,j)
						call mfft(r,1,map%shap,r1(1:nkD,i,j),1d0/product(map%shap))
						r=r2(1:nkD,i,j)
						call mfft(r,-1,map%shap,r2(1:nkD,i,j),1d0/product(map%shap))
					enddo
				enddo
				do i=1,Nc
#ifdef cluster
					do j=1,Nc**2
						r=r1(1:nkD,i,mod(j-1,Nc)+1)*r2(1:nkD,(j-1)/Nc+1,i)
#else
					do j=1,Nc
						r=r1(1:nkD,i,j)*r2(1:nkD,j,i)
#endif
						call mfft(r,1,map%shap,LL(iw,i,j,1:nkD))
#ifdef cluster
						!r=r2(1:nkD,i,(j-1)/Nc+1)*r1(1:nkD,mod(j-1,Nc)+1,i)
						r=r2(1:nkD,(j-1)/Nc+1,i)*r1(1:nkD,i,mod(j-1,Nc)+1)
#else
						!r=r2(1:nkD,i,j)*r1(1:nkD,j,i)
						r=r2(1:nkD,j,i)*r1(1:nkD,i,j)
#endif
						call mfft(r,-1,map%shap,RR(iw,j,i,1:nkD))
					enddo
				enddo
				do i=1,Nc
					do j=1,Nc
						call mfft(self%Gck(1:nkD,iw,i,j),1,map%shap,r1(1:nkD,i,j),1d0/product(map%shap))
						call mfft(self%Gck(1:nkD,iw,j,i),-1,map%shap,r2(1:nkD,j,i),1d0/product(map%shap))
						r=-r1(1:nkD,i,j)*r2(1:nkD,j,i)
						call mfft(r,1,map%shap,self%Sus0c(1:nkD,iw,i,j))
					enddo
				enddo
			enddo
			!$OMP END PARALLEL DO
			self%Sus=0._wp
			self%Susc=0._wp
			if(this_image()==1) call start_time(id=3)
			do i=this_image(),merge(nkD,size(self%ukc2kc,1),self%is_symbk),num_images()
				if(.not.omp_in_parallel()) then
					write(*,"(i5$)")i
				endif
				if(.not.omp_in_parallel()) then
					call mkl_set_num_threads(omp_get_max_threads())
				endif
				!if(.not.omp_in_parallel()) then
					!call start_time(id=4)
				!endif
				ik=merge(i,self%ukc2kc(i,1),self%is_symbk)
				!call get_dF4(self%G4c(1:n*Nc*Nc,1:n*Nc*Ns),self%dSus0(self%ukc2kc(i,1),1:n*Nc*Nc,1:Nc*dNs),self%dF4(1:n*Nc**2,1:n*Nc*dNs))
				call get_dF4(self%G4c,self%dSus0(ik,:,:),self%dF4,[1,-1])
				!if(.not.omp_in_parallel()) then
					!write(*,"(A$)")"inverse of "
					!call stop_time(id=4,show=underscore)
				!endif
				if(.not.omp_in_parallel()) then
					call mkl_set_num_threads(1)
				endif

				!$OMP PARALLEL DO PRIVATE(iAc,cnp,cn,tmpD,ik,Nc2,cNc4)
				do j=1,merge(1,self%ukc2kc(i,0),self%is_symbk)
#ifdef cluster
						Nc2=merge([1:Nc**2],[(crt(1:Nc,self%ukc2kc(i,j+self%ukc2kc(i,0)))+(crt(iQ,self%ukc2kc(i,j+self%ukc2kc(i,0)))-1)*Nc,iQ=1,Nc)],self%is_symbk)
#else
						Nc2=merge([1:Nc],crt(1:Nc,self%ukc2kc(i,j+self%ukc2kc(i,0))),self%is_symbk)
#endif

					ik=merge(i,self%ukc2kc(i,j),self%is_symbk)
					do i1=1,Nc
						do j1=1,Nc
							do iw=1,n
								cn(iw)=0._wp
								do jw=1,n
									cNc4(Nc2,Nc2)=self%dF4(iw::n,jw::n)
#ifdef cluster
									cnp(jw)=sum(LL(iw,i1,:,ik)*matmul(cNc4,RR(jw,:,j1,ik)))
#else
									cnp(jw)=sum(LL(iw,i1,:,ik)*matmul(cNc4(1:Nc,1:Nc),RR(jw,:,j1,ik)))
#endif
								enddo
								cn(iw)=cn(iw)+msum(cnp(1:n))
							enddo
							self%Susc(ik,i1,j1)=Tk*msum(cn)
						enddo
					enddo

				enddo
				!$OMP END PARALLEL DO
			enddo
			call coarray([node])
			do ik=1,nkD
				do iQ=1,Nc
					if(is_CDMFT.or.(.not.is_dual_pureDCA)) then
						tmpD=self%k(:,ik+(iQ-1)*nkD)
					elseif(is_dual_pureDCA) then
						tmpD=self%Qc(:,mod(iQ-1,Nc/Ns)+1)+self%Qs(:,(iQ-1)/(Nc/Ns)+1)
					endif
					self%Sus(ik+(iQ-1)*nkD)=sum(exp(img*matmul(tmpD,self%rc(1:D,1:Nc)))*matmul(self%Susc(ik,1:Nc,1:Nc),exp(-img*matmul(tmpD,self%rc(1:D,1:Nc)))))/Nc
				enddo
			enddo
			do j=1,Nc
				call mfft(self%Susc(1:nkD,1,j),1,map%shap,r(1:nkD),1._wp/product(map%shap))
				self%Susr((j-1)*nkD+1:(j)*nkD)=r(1:nkD)
			enddo
			if(this_image()==1) then
				write(*,"(x)")
				write(*,"(A$)")"Sus of "
				call stop_time(id=3,show=underscore)
				do ik=1,size(self%uk2k,1)
					if(any(abs(self%Sus(self%uk2k(ik,1))-self%Sus(self%uk2k(ik,1:self%uk2k(ik,0))))>1e-7_wp)) then
						write(*,*)"warning!!! rotation symmetry broken",maxval(abs(self%Sus(self%uk2k(ik,1))-self%Sus(self%uk2k(ik,1:self%uk2k(ik,0)))))
					endif
				enddo
			endif
		case(Gc)
			rt="Gc"; if(debug) return
			call solver(self%Delta,self%mu,self%Ef,self%Gc,self%ife,nf=self%n_f)
		case(dSEc)
			rt="dSEc"; if(debug) return
#ifdef cluster
			if(.not.allocated(self%dF4d)) then
				allocate(self%dF4d(n,Nc**2,Nc**2,nkD))
			endif
			if(.not.allocated(self%dF4)) then
				allocate(self%dF4(n*Nc*Nc,n*Nc*dNs))
			endif
			if(.not.allocated(self%dSus0)) then
				allocate(self%dSus0(nk**D,n*Nc*Nc,Nc*dNs))
			endif
			if(.not.allocated(self%G4c)) then
				allocate(self%G4c(n*Nc*Nc,n*Nc*Ns))
			endif
#else
			if(.not.allocated(self%dF4d)) then
				allocate(self%dF4d(n,Nc,Nc,nkD))
			endif
			if(.not.allocated(self%dF4)) then
				allocate(self%dF4(n*Nc,n*Nc))
			endif
			if(.not.allocated(self%dSus0)) then
				allocate(self%dSus0(nk**D,n*Nc,Nc))
			endif
			if(.not.allocated(self%G4c)) then
				allocate(self%G4c(n*Nc,n*Nc))
			endif
#endif
			if(sdiagram(2)/=0) then
			do ik=this_image(),merge(nkD,size(self%ukc2kc,1),self%is_symbk),num_images()
			!do ik=1,size(self%ukc2kc,1)
				if(.not.omp_in_parallel()) then
					call mkl_set_num_threads(omp_get_max_threads())
				endif
				i=merge(ik,self%ukc2kc(ik,1),self%is_symbk)
				!call get_dF4(self%G4c(1:n*Nc*Nc,1:n*Nc*Ns),self%dSus0(self%ukc2kc(ik,1),1:n*Nc*Nc,1:Nc*dNs),self%dF4(1:n*Nc**2,1:n*Nc*dNs))
				call get_dF4(self%G4c,self%dSus0(i,:,:),self%dF4,sdiagram)
				if(.not.omp_in_parallel()) then
					call mkl_set_num_threads(1)
				endif
				!$OMP PARALLEL DO
				do iw=1,n
#ifdef cluster
					self%dF4d(iw,1:Nc*Nc,1:Nc*Nc,i)=self%dF4(iw:n*Nc*Nc:n,iw:n*Nc*Nc:n)
#else
					self%dF4d(iw,1:Nc,1:Nc,i)=self%dF4(iw:n*Nc:n,iw:n*Nc:n)
#endif
				enddo
				!$OMP END PARALLEL DO
			enddo
			call coarray([dF4d])
			!if(this_image()==1) then
			if(.not.self%is_symbk) then
				do ik=1,size(self%ukc2kc,1)
					!$OMP PARALLEL DO PRIVATE(Nc2)
					do i=2,self%ukc2kc(ik,0)
#ifdef cluster
						Nc2=[(crt(1:Nc,self%ukc2kc(ik,i+self%ukc2kc(ik,0)))+(crt(j,self%ukc2kc(ik,i+self%ukc2kc(ik,0)))-1)*Nc,j=1,Nc)]
#else
						Nc2=crt(1:Nc,self%ukc2kc(ik,i+self%ukc2kc(ik,0)))
#endif
						self%dF4d(1:n,Nc2,Nc2,self%ukc2kc(ik,i))=self%dF4d(1:n,:,:,self%ukc2kc(ik,1))
					enddo
					!$OMP END PARALLEL DO
				enddo
			endif
			endif
			self%dSEc=0._wp
			!$OMP PARALLEL DO PRIVATE(r,r1,r2,dSus0d,iAc2,iBc2,iAc,iBc,cn)
			!do iw=1,n
			do iw=this_image(),n,num_images()
				do i=1,Nc**2
					call mfft(self%dGc(1:nkD,iw,mod(i-1,Nc)+1,(i-1)/Nc+1),1,map%shap,r2(1:nkD,mod(i-1,Nc)+1,(i-1)/Nc+1),1d0/product(map%shap))
				enddo
#ifdef cluster
				if(sdiagram(2)/=0) then
				do i=1,Nc**2
					do j=1,Nc*dNs
						call mfft(self%dF4d(iw,mod(i-1,Nc)+1+mod(j-1,Nc)*Nc,(j-1)/Nc+1+((i-1)/Nc)*Nc,1:nkD),-1,map%shap,r1(1:nkD,mod(j-1,Nc)+1,(j-1)/Nc+1),1d0/product(map%shap))
						! need to check the prefactor
						self%dSEc(1:nkD,iw,mod(i-1,Nc)+1,(i-1)/Nc+1)=self%dSEc(1:nkD,iw,mod(i-1,Nc)+1,(i-1)/Nc+1)+r1(1:nkD,mod(j-1,Nc)+1,(j-1)/Nc+1)*r2(1:nkD,mod(j-1,Nc)+1,(j-1)/Nc+1)!*Nc

					enddo
				enddo
				endif
				if(ddiagram(2)/=0) then
					! dynamic contribution
					do jw=1,n
						do i=1,Nc**2
							call mfft(self%dGc(1:nkD,jw,mod(i-1,Nc)+1,(i-1)/Nc+1),-1,map%shap,r1(1:nkD,mod(i-1,Nc)+1,(i-1)/Nc+1),1d0/product(map%shap))
						enddo
						do i=1,Nc**2
							do j=1,Nc**2 
								r=r1(1:nkD,(i-1)/Nc+1,mod(j-1,Nc)+1)*r2(1:nkD,(j-1)/Nc+1,mod(i-1,Nc)+1)
								call mfft(r(1:nkD),-1,map%shap,dSus0d(1:nkD,i,j))
							enddo
						enddo
						do ik=1,nkD
							do i=1,Nc**2
								do j=1,Nc**2
									iAc2((j-1)/Nc*Nc+mod(i-1,Nc)+1,(i-1)/Nc*Nc+mod(j-1,Nc)+1)=self%G4c((i-1)*n+iw,(j-1)*n+jw)
								enddo
							enddo
							iBc2=matmul(iAc2,dSus0d(ik,:,:))
							dSus0d(ik,:,:)=iAc2
							iAc2=diag(1._wp,Nc**2)-iBc2
							call gesv(iAc2,dSus0d(ik,:,:))
							!dSus0d(ik,:,:)=matmul(dSus0d(ik,:,:),iBc2)
						enddo
						do i=1,Nc**2
							do j=1,Nc**2
								call mfft(dSus0d(1:nkD,mod(j-1,Nc)*Nc+mod(i,Nc)+1,(i-1)/Nc*Nc+(j-1)/Nc+1),-1,map%shap,r1(1:nkD,mod(j-1,Nc)+1,(j-1)/Nc+1),1d0/product(map%shap))
								self%dSEc(1:nkD,iw,mod(i-1,Nc)+1,(i-1)/Nc+1)=self%dSEc(1:nkD,iw,mod(i-1,Nc)+1,(i-1)/Nc+1)-(-Res(jw))*r1(1:nkD,mod(j-1,Nc)+1,(j-1)/Nc+1)*r2(1:nkD,mod(j-1,Nc)+1,(j-1)/Nc+1)
							enddo
						enddo
					enddo
				endif

#else
				if(sdiagram(2)/=0) then
				do i=1,Nc
					do j=1,Nc
						call mfft(self%dF4d(iw,i,j,1:nkD),-1,map%shap,r1(1:nkD,i,j),1d0/product(map%shap))
						self%dSEc(1:nkD,iw,i,j)=self%dSEc(1:nkD,iw,i,j)+r1(1:nkD,i,j)*r2(1:nkD,i,j)
					enddo
				enddo
				endif

				if(ddiagram(2)/=0) then
					!! dynamic contribution
					!do jw=1,n
					do jw=1,n
						!if(abs(Res(jw)+1._wp)>1e-6_wp) then
							!cycle
						!endif
						do i=1,Nc**2
							call mfft(self%dGc(1:nkD,jw,mod(i-1,Nc)+1,(i-1)/Nc+1),-1,map%shap,r1(1:nkD,mod(i-1,Nc)+1,(i-1)/Nc+1),1d0/product(map%shap))
						enddo
						do i=1,Nc
							do j=1,Nc
								r=r1(1:nkD,i,j)*r2(1:nkD,j,i)
								!call mfft(r(1:nkD),-1,map%shap,dSus0d(1:nkD,i,j))
								call mfft(r(1:nkD),1,map%shap,dSus0d(1:nkD,i,j))
							enddo
						enddo
						do i=1,Nc
							iAc(i,1)=self%G4c((i-1)*n+iw,(i-1)*n+jw)
						enddo
						do ik=1,nkD
							do i=1,Nc
								iBc(:,i)=iAc(:,1)*dSus0d(ik,:,i)
							enddo
							if(ddiagram(2)<0) then
								dSus0d(ik,:,:)=diag(1._wp,Nc)-iBc
								call mat_inv(dSus0d(ik,:,:))
								do i=1,Nc
									dSus0d(ik,:,i)=iAc(:,1)*dSus0d(ik,:,i)
								enddo
							else
								dSus0d(ik,:,:)=diag(iAc(:,1))
							endif
							do i=2,ddiagram(2)-ddiagram(1)+1
								dSus0d(ik,:,:)=diag(iAc(:,1))+matmul(iBc,dSus0d(ik,:,:))
							enddo
							do i=2,ddiagram(1)
								dSus0d(ik,:,:)=matmul(iBc,dSus0d(ik,:,:))
							enddo
						enddo
						do i=1,Nc
							do j=1,Nc
								call mfft(dSus0d(1:nkD,i,j),-1,map%shap,r(1:nkD),1d0/product(map%shap))
								self%dSEc(1:nkD,iw,i,j)=self%dSEc(1:nkD,iw,i,j)-(-Res(jw))*r(1:nkD)*r1(1:nkD,i,j)
							enddo
						enddo
					enddo
				endif
#endif
				do i=1,Nc
					do j=1,Nc
						r=self%dSEc(1:nkD,iw,i,j)
						call mfft(r,-1,map%shap)
						self%dSEc(1:nkD,iw,i,j)=r
					enddo
				enddo
			enddo
			!$OMP END PARALLEL DO
			call coarray([dSEc+prime])
		case(Ef)
			rt="Ef"; if(debug) return
		case(mu)
			rt="mu"; if(debug) return
		case(Delta)
			rt="Delta"; if(debug) return
			!call self%io([out_w],[out_w])
			select case(sc_scheme)
			case(1)
				!$OMP PARALLEL DO PRIVATE(iAc,iBc)
				do iw=1,n
					iAc=0._wp
					do ik=1,nkD
						iAc=iAc+self%dGc(ik,iw,1:Nc,1:Nc)
					enddo
					iAc=iAc/nkD
					iBc=self%Gcl(iw,1:Nc,1:Nc)
					!do i=1,Nc
						!do j=1,Nc
							!if(i/=j) then
								!iBc(i,j)=0._wp
							!endif
						!enddo
					!enddo
#ifdef cluster
					call mat_inv(iBc)
					iBc=matmul(iAc,iBc)
					iAc=self%Gc(iw,1:Nc,1:Nc)
					call mat_inv(iAc)
					self%Delta(iw,1:Nc,1:Nc)=self%Delta(iw,1:Nc,1:Nc)+matmul(iAc,iBc)
					!do i=1,Nc
						!do j=1,Nc
							!if(i/=j) then
								!self%Delta(iw,i,j)=0._wp
							!endif
						!enddo
					!enddo
#else
					do i=1,Nc
						self%Delta(iw,i)=self%Delta(iw,i)+iAc(i,i)/(iBc(i,i)*self%Gc(iw,i))
					enddo
#endif
				enddo
				!$OMP END PARALLEL DO
			case(2)
				!$OMP PARALLEL DO PRIVATE(iAc,iBc)
				do iw=1,n
					iBc=self%Gcl(iw,1:Nc,1:Nc)
#ifdef cluster
					iAc=self%Gc(iw,1:Nc,1:Nc)
					call mat_inv(iBc)
					call mat_inv(iAc)
					self%Delta(iw,1:Nc,1:Nc)=self%Delta(iw,1:Nc,1:Nc)+iAc-iBc
#else
					do i=1,Nc
						self%Delta(iw,i)=self%Delta(iw,i)+1._wp/self%Gc(iw,i)-1._wp/iBc(i,i)
					enddo
#endif
				enddo
				!$OMP END PARALLEL DO
			end select
			self%is_symbk=.false.
			do iw=1,n
#ifdef cluster
				if(any(abs(diag(self%Delta(iw,1:Nc,1:Nc))-self%Delta(iw,1,1))>1e-6_wp)) then
#else
				if(any(abs(self%Delta(iw,1:Nc)-self%Delta(iw,1))>1e-6_wp)) then
#endif
					self%is_symbk=.true.
					exit
				endif
			enddo
		case(:-1)
			select case(node)
			case(SC_DMFT)
				rt="DMFT"
			case(SC_dDMFT)
				rt="dDMFT"
			case(SC_cp)
				rt="cp"
			end select
			if(debug) return
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
				if(this_image()==1) write(*,"('SC=',A5,' T=',es10.3,' U=',es10.3,' initial...')")pack(["DMFT","dDMFT","cp"],[SC_DMFT,SC_dDMFT,SC_cp]==node),Tk,U
				if(onode(node,1)==Ef.and.size(x)==1) then
					if(next(ytol=iter_err(node),dx=iter_rate(node),id=abs(node))) then
					endif
				else
					call mbroyden(0,x,x_,iter_rate(node),nsave=30,id=abs(node))
				endif
			endif
			if(this_image()==1) write(*,"('SC=',A5,' err=',es10.3,' iter=',i4,' nf=',es10.3,' Ef=',es10.3,' nc=',es10.3,' mu=',es10.3)")pack(["DMFT","dDMFT","cp"],[SC_DMFT,SC_dDMFT,SC_cp]==node),iter_err(abs(node)),iter_exit(node),sum(self%n_f)/Nc,self%Ef,sum(self%n_c)/Nc,self%mu
			if(iter_exit(node)==niter(node).or.iter_err(abs(node))<=iter_err(node)) then
				call mbroyden(huge(1),x,x_,id=abs(node))
				deallocate(x,x_)
				if(this_image()==1) write(*,"('SC=',A5,' T=',es10.3,' U=',es10.3,' finish...')")pack(["DMFT","dDMFT","cp"],[SC_DMFT,SC_dDMFT,SC_cp]==node),Tk,U
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
				!call coarray([cnd_,node])
				call coarray([cnd_,onode(cnd_,1:onode(cnd_,0))])
				!call self%io([inout_save],[inout_save])
				cnd_=0
			endif
		endif

		!do l=1,inode(node,0)
			!if(inode(node,l)<0) then
				!call x_x([node],xx_diff)
				!if(nx==size(x)) then
					!call x_x(onode(inode(node,l),1:onode(inode(node,l),0)),xx_new)
					!call coarray([inode(node,l),node])
				!endif
				!exit
			!endif
		!enddo
		!if(debug) read(*,*)
	contains
		subroutine x_x(nd,flag)
			integer :: nd(:),flag
			integer :: l
			if(flag==xx_new) then
				iter_err(abs(cnd_))=maxval(abs(x_))
				if(this_image()==1) then
					if(size(nd)==1.and.nx==1.and.nd(1)==Ef) then
						underscore=next(self%Ef,-x_(nx)%re,id=abs(cnd_))
					else
						if(niter(cnd_)>0) then
							call mbroyden(iter_exit(cnd_),x,x_,id=abs(cnd_))
							!x=x_+x
						else
							x=x+iter_rate(cnd_)*x_
						endif
					endif
				endif
				nx=0
			endif
			do l=1,size(nd)
				select case(nd(l))
				case(dSEc)
					if(allocated(x)) then
						do iw=1,n
							do i=1,Nc**2
								select case(flag)
								case(xx_save)
									x(nx+iw+(i-1)*n::Nc**2*n)=self%dSEc(1:nkD,iw,mod(i-1,Nc)+1,(i-1)/Nc+1)
								case(xx_diff)
									x_(nx+iw+(i-1)*n::Nc**2*n)=(self%dSEc(1:nkD,iw,mod(i-1,Nc)+1,(i-1)/Nc+1)-x(nx+iw+(i-1)*n::Nc**2*n))
								case(xx_new)
									self%dSEc(1:nkD,iw,mod(i-1,Nc)+1,(i-1)/Nc+1)=x(nx+iw+(i-1)*n::Nc**2*n)
								end select
							enddo
						enddo
						!if(flag==xx_diff) write(*,*)"dSEc",maxval(abs(x_(nx+1:nx+size(self%dSEc))))
					endif
					nx=nx+size(self%dSEc)
				case(Delta)
					if(allocated(x)) then
						do iw=1,n
							select case(flag)
							case(xx_save)
#ifdef cluster
								x(nx+Ns*Nc*(iw-1)+1:nx+Ns*Nc*iw)=reshape(self%Delta(iw,1:Ns,1:Nc),[Ns*Nc])
#else
								x(nx+Ns*(iw-1)+1:nx+Ns*iw)=self%Delta(iw,1:Ns)
#endif
							case(xx_diff)
#ifdef cluster
								x_(nx+Ns*Nc*(iw-1)+1:nx+Ns*Nc*iw)=reshape(self%Delta(iw,1:Ns,1:Nc),[Ns*Nc])-x(nx+Ns*Nc*(iw-1)+1:nx+Ns*Nc*iw)
#else
								x_(nx+Ns*(iw-1)+1:nx+Ns*iw)=self%Delta(iw,1:Ns)-x(nx+Ns*(iw-1)+1:nx+Ns*iw)
#endif
							case(xx_new)
#ifdef cluster
								self%Delta(iw,1:Ns,1:Nc)=reshape(x(nx+Ns*Nc*(iw-1)+1:nx+Ns*Nc*iw),[Ns,Nc])
								do i=Ns+1,Nc
									self%Delta(iw,i,1:Nc)=self%Delta(iw,mod(i-1,Ns)+1,pm(i,1:Nc))
								enddo
#else
								self%Delta(iw,1:Ns)=x(nx+Ns*(iw-1)+1:nx+Ns*iw)
								do i=Ns+1,Nc
									self%Delta(iw,i)=self%Delta(iw,mod(i-1,Ns)+1)
								enddo
#endif
							end select
						enddo
						if(flag==xx_new) then
							self%is_symbk=.false.
							do iw=1,n
#ifdef cluster
								if(any(abs(diag(self%Delta(iw,1:Nc,1:Nc))-self%Delta(iw,1,1))>1e-6_wp)) then
#else
								if(any(abs(self%Delta(iw,1:Nc)-self%Delta(iw,1))>1e-6_wp)) then
#endif
									self%is_symbk=.true.
									exit
								endif
							enddo
						endif
						!if(flag==xx_diff) write(*,*)"Delta",maxval(abs(x_(nx+1:nx+size(self%Delta)/Nc*Ns)))
					endif
					nx=nx+size(self%Delta)/Nc*Ns
				case(Ef)
					! Ef self-consistent
					if(allocated(x)) then
						select case(flag)
						case(xx_save)
							x(nx+1)=self%Ef
							!x(nx+1)%re=self%Ef+self%mu
						case(xx_diff)
							x_(nx+1)=((sum(self%n_f)/Nc-n_f))*1._wp
							!x_(nx+1)=1._wp/(exp(6._wp*(sum(self%n_f)/Nc-n_f))+1._wp)-0.5_wp
							!x_(nx+1)%re=-(sum(self%n_f-self%n_c)/Nc-(n_f-n_c))*10._wp
						case(xx_new)
							if(size(x)/=1) then
								!if(abs(real(x(nx+1))-(self%Ef+self%mu))<0.05_wp) then
									!self%Ef=real(x(nx+1))-self%mu
								!else
									!self%Ef=self%Ef+0.05_wp*sign(1._wp,real(x(nx+1))-(self%Ef+self%mu))
								!endif
								if(abs(real(x(nx+1))-self%Ef)<0.05_wp) then
									self%Ef=real(x(nx+1))
								else
									self%Ef=self%Ef+0.05_wp*sign(1._wp,real(x(nx+1))-self%Ef)
								endif
								!else
									!self%Ef=-U*0.5_wp
								!endif
								!self%Ef=(x(nx+1)%re+x(nx+1)%im)*0.5_wp
							endif
						end select
					endif
					nx=nx+1
				case(mu)
					if(allocated(x)) then
						select case(flag)
						case(xx_save)
							x(nx+1)=self%mu
							!x(nx+1)%im=self%Ef-self%mu
						case(xx_diff)
							x_(nx+1)=-((sum(self%n_c)/Nc-n_c))*1._wp
							!x_(nx+1)=1._wp/(exp(6._wp*(sum(self%n_c)/Nc-n_c))+1._wp)-0.5_wp
							!x_(nx+1)%im=-(sum(self%n_f+self%n_c)/Nc-(n_f+n_c))*10._wp
						case(xx_new)
							if(abs(real(x(nx+1))-self%mu)<0.05_wp) then
								self%mu=real(x(nx+1))
								!self%mu=(x(nx+1)%re-x(nx+1)%im)*0.5_wp
							else
								self%mu=self%mu+0.05_wp*sign(1._wp,real(x(nx+1))-self%mu)
							endif
						end select
					endif
					nx=nx+1
				case default
					write(*,*)"not complete for self consistent: ",nd
					stop
				end select
			enddo
		end subroutine
		subroutine gen_grid(a1,a2,r0,T,grid)
			real(wp) :: grid(:,:)
			real(wp) :: r0(:),a1(:),a2(:),T(:,:)
			real(wp) :: r(2),x(size(T,2)),y(size(T,2)),tf(2,2)
			integer :: i,j,n
			n=0
			tf=reshape((/a2(2),-a1(2),-a2(1),a1(1)/)/(a1(1)*a2(2)-a1(2)*a2(1)),(/2,2/))
			do i=1,size(T,2)
				x(i)=sum(tf(1,:)*(T(:2,i)-r0(:2)))
				y(i)=sum(tf(2,:)*(T(:2,i)-r0(:2)))
			enddo
			do j=nint(minval(y)),nint(maxval(y))
				do i=nint(minval(x)),nint(maxval(x))
					r=a1*i+a2*j+r0
					if(is_in(r,T)) then
						n=n+1
						if(size(grid,2)<n) then
							write(*,*)"error number of site is small"
						else
							grid(:,n)=r
						endif
					endif
				enddo
			enddo
			if(n/=size(grid,2)) then
				write(*,*)"error eq: ",n,size(grid,2)
				stop
			endif
		end subroutine
		subroutine get_wigner_seitz(a,T)
			real(wp) :: a(:,:)
			real(wp), allocatable :: T(:,:)
			real(wp) :: tf(2,2),T_(2,0:24),th,dr1(2),dr2(2)
			real(wp) :: st_dr(2,24),st_val(24)
			integer :: n,i,j,ord(24)
			integer, allocatable :: ist(:)
			n=0
			do i=-2,2
				do j=-2,2
					if(i==0.and.j==0) then
						cycle
					endif
					n=n+1
					st_dr(:,n)=(a(:,1)*i+a(:,2)*j)/2._wp
					st_val(n)=theta(st_dr(:,n))
				enddo
			enddo
			ord(1:n)=[1:n]
			call qsort(st_val(1:n),ord)
			st_val(1:n)=st_val(ord(1:n))
			st_dr(:,1:n)=st_dr(:,ord(1:n))
			call collect(st_val(1:n),ist)
			do i=1,size(ist)-1
				st_dr(:,i)=st_dr(:,ist(i))
				do j=ist(i)+1,ist(i+1)-1
					if(sqrt(sum(st_dr(:,j)**2))<sqrt(sum(st_dr(:,i)**2))) then
						st_dr(:,i)=st_dr(:,j)
					endif
				enddo
			enddo
			n=1
			T_(:2,1)=st_dr(:2,1)
			do i=2,(size(ist)-1)/2
				dr1=T_(:2,n)
				dr2=st_dr(:2,i)
				th=theta(matmul(reshape((/dr2(2),-dr2(1),-dr1(2),dr1(1)/),(/2,2/)),(/sum(dr1**2),sum(dr2**2)/))/(dr1(1)*dr2(2)-dr1(2)*dr2(1)))
				if(sin(th-theta(dr1)-100._wp*eps)<0._wp) then
					n=n-1
				endif
				if(sin(th-theta(dr2)+100._wp*eps)<0._wp) then
					n=n+1
					T_(:2,n)=dr2
				endif
			enddo
			T_(:,n+1:)=-T_(:,1:n)
			n=n*2
			if(allocated(T)) deallocate(T)
			allocate(T(2,n))
			do i=1,n
				j=mod(i,n)+1
				dr1=T_(:2,i)
				dr2=T_(:2,j)
				T(:2,i)=matmul(reshape((/dr2(2),-dr2(1),-dr1(2),dr1(1)/),(/2,2/)),(/sum(dr1**2),sum(dr2**2)/))/(dr1(1)*dr2(2)-dr1(2)*dr2(1))
			enddo
		end subroutine
		subroutine into_sym(k,n)
			real(8) :: k(:)
			integer :: n
			integer :: i
			!return
			i=size(k)-1
			do 
				if(theta(k(1:i))<2._wp*pi/abs(n)-1e-8_wp.and.theta(k(1:i))>=-1e-8_wp) then
					if(n<0) then
						k(1:i)=rotate(k(1:i),-pi/abs(n))
						if(k(2)>1e-8_wp) then
							k(i+1)=k(i+1)+50._wp
						endif
						k(2)=abs(k(2))
						k(1:i)=rotate(k(1:i),pi/abs(n))
					endif
					exit
				else
					k(1:i)=rotate(k(1:i),-2._wp*pi/abs(n))
					k(i+1)=k(i+1)+1._wp
				endif
			enddo
		end subroutine
		subroutine get_uk(k,uk2k,iCnv)
			real(8) :: k(:,:)
			integer, allocatable :: uk2k(:,:)
			integer, optional :: iCnv
			real(8) :: tmp(1:D+1,size(k,2)),T(D,12),rt_(D),rt0(D)
			integer :: n,i,j,m,Cnv_
			integer :: ord(size(k,2))
			integer, allocatable :: c(:)
			if(.not.(Nc>1.and.size(k,2)==nkD)) then
				n=size(Ta,2)
				T(:,1:n)=Ta(:,1:n)
			else
				n=size(Tc,2)
				T(:,1:n)=Tc(:,1:n)
			endif
			if(present(iCnv)) then
				Cnv_=iCnv
			else
				if(any(abs(sum(T(:,1)**2,1)-sum(T(:,1:n)**2,1))>1e-8_wp)) then
					write(*,*)"T is not equal"
					stop
				endif
				do Cnv_=1,n-1
					do i=Cnv_+1,n
						if(abs(sum(T(:,i-Cnv_)*T(:,i))-sum(T(:,1)**2)*cos(2._wp*pi/n*Cnv_))>1e-8_wp) then
							exit
						endif
					enddo
					if(i==n+1) then
						exit
					endif
				enddo
				if(Cnv_>2) then
					write(*,*)"symmetry for large n not considered, please check lattice!"
					stop
				endif
				if(Cnv_==1) then
					Cnv_=-1
				else
					if(abs(sum(abs(T(:,1)+T(:,n)))-sum(abs(T(:,1+Cnv_)+T(:,Cnv_))))<1e-8_wp&
						.and.&
						(abs(sum(((T(:,1)+T(:,n))-(T(:,1+Cnv_)+T(:,Cnv_)))*(T(:,1)+T(:,2))))<1e-8_wp)&

						) then
						Cnv_=-abs(Cnv_)
					endif
				endif
				Cnv_=n/Cnv_
			endif
			if(Nc>1.and.size(k,2)==nkD) then
				rt0=sum(self%rc(1:D,1:Nc),2)/Nc
				i=1
				do 
					if(mod(abs(Cnv_),i)==0) then
						do j=1,Nc
							rt_=self%rc(1:D,j)-rt0
							rt_=rotate(rt_,2._wp*pi/abs(Cnv_)*i)
							do m=1,Nc
								if(sum(abs(self%rc(1:D,m)-rt0-rt_))<1e-8_wp) then
									exit
								endif
							enddo
							if(m==Nc+1) then
								write(*,*)"The symmetry of cluster is smaller than orignal lattice"
								exit
							endif
						enddo
						if(j==Nc+1) exit
					endif
					i=i+1
					if(i>abs(Cnv)) then
						stop
					endif
				enddo
				Cnv_=Cnv_/i
				if(Cnv_<0) then
					do j=1,Nc
						rt_=self%rc(1:D,j)-rt0
						rt_=rotate(rt_,-theta(T(:,1)))
						rt_=rotate(rt_,-pi/abs(Cnv_))
						rt_(2)=-rt_(2)
						rt_=rotate(rt_,pi/abs(Cnv_))
						do m=1,Nc
							if(sum(abs(self%rc(1:D,m)-rt0-rt_))<1e-8_wp) then
								exit
							endif
						enddo
						if(m==Nc+1) then
							Cnv_=abs(Cnv_)
							exit
						endif
					enddo
				endif
			endif

			do i=1,size(k,2)
				tmp(1:D,i)=k(:,i)
				tmp(D+1,i)=0._wp
				if(sum(abs(tmp(:,i)))>1e-8_wp) then
					tmp(1:D,i)=rotate(tmp(1:D,i),-theta(T(:,1)))
					call into_sym(tmp(1:,i),Cnv_)
				endif
				ord(i)=i
			enddo
			!write(*,"(*(i4))")ord
			call qsort(tmp(:,:),ord)
			!write(*,"(*(i4))")ord
			tmp(:,:)=tmp(:,ord)
			call collect(tmp(1:D,:),c)
			allocate(uk2k(size(c)-1,0:maxval(c(2:)-c(1:size(c)-1))*2))
			uk2k(:,0)=c(2:)-c(1:size(c)-1)
			do i=1,size(uk2k,1)
				uk2k(i,1:c(i+1)-c(i))=ord(c(i):c(i+1)-1)
				!Following is done in this way intentionally due to bug in intel compiler. See https://software.intel.com/en-us/forums/intel-fortran-compiler/topic/842733
				do j=1,uk2k(i,0)
					uk2k(i,uk2k(i,0)+j)=(mod(nint(tmp(D+1,c(i)+j-1)),50)+1)*(1-nint(tmp(D+1,c(i)+j-1))/50*2)
				enddo
				!uk2k(i,uk2k(i,0)+1:uk2k(i,0)*2)=(mod(nint(tmp(D+1,c(i):c(i+1)-1)),50)+1)*(1-nint(tmp(D+1,c(i):c(i+1)-1))/50*2)
			enddo
			if(sum(uk2k(:,0))/=size(k,2)) then
				write(*,*)"get_uk error2"
				stop
			endif
			if(Nc>1.and.size(k,2)==nkD) then
				Cnvc=Cnv_
				allocate(crt(Nc,merge(0,Cnv_,Cnv_>0):abs(Cnv_)))
				crt=0
				!write(*,*)theta(T(:,1))/pi*180,180/abs(Cnv_)
				do i=lbound(crt,2),ubound(crt,2)
					if(i/=0) then
						do j=1,Nc
							rt_=self%rc(1:D,j)-rt0
							rt_=rotate(rt_,-theta(T(:,1)))
							if(i<0) then
								rt_=rotate(rt_,-pi/abs(Cnv_))
								rt_(2)=-rt_(2)
								rt_=rotate(rt_,pi/abs(Cnv_))
							endif
							if(abs(i)/=1) then
								rt_=rotate(rt_,2._wp*pi/abs(Cnv_)*(abs(i)-1))
							endif
							rt_=rotate(rt_,theta(T(:,1)))
							do m=1,Nc
								if(sum(abs(self%rc(1:D,m)-rt0-rt_))<1e-8_wp) then
									if(crt(j,i)==0) then
										crt(j,i)=m
									else
										write(*,*)"crt error 1!",i
									endif
								endif
							enddo
							if(crt(j,i)==0) then
								write(*,*)"crt error 2!"
							endif
						enddo
					endif
					!write(*,"(*(i3))")i,crt(1:Nc,i)
				enddo
				!stop
			elseif(Nc==1) then
				Cnv=Cnv_
				Cnvc=Cnv_
				allocate(self%ukc2kc(size(uk2k,1),0:(size(uk2k,2)-1)))
				self%ukc2kc=self%uk2k
				allocate(crt(Nc,-abs(Cnv_):abs(Cnv_)))
				crt=1
			else
				Cnv=Cnv_
			endif
		end subroutine
		function rotate(r,theta) result(rt)
			! anticlock-wise rotation
			real(wp) :: r(2),theta
			real(wp) :: rt(2)
			rt=r*cos(theta)+[-r(2),r(1)]*sin(theta)
		end function
		function theta(r)
			real(wp) :: r(:),theta,d
			d=sqrt(sum(r**2))
			if(d<1e-10_wp) then
				theta=0._wp
				return
			endif
			theta=acos(r(1)/d)
			if(r(2)<-1e-8_wp) then
				theta=2._wp*pi-theta
			endif
		end function
		subroutine get_kline(k,from,to,idx)
			real(wp) :: k(:,:),from(:),to(:)
			integer, allocatable :: idx(:)
			integer :: i,l,idx_(size(k,2)+size(from)/size(k,1)+1),m(2),ord(size(k,2)),sl
			real(wp) :: k_(size(k,2))
			m(1)=0
			m(2)=0
			sl=size(from)/size(k,1)
			idx_(size(idx_))=0
			do l=0,sl-1
				do i=1,size(k,2)
					if(abs((to(2+l*D)-from(2+l*D))*k(1,i)-(to(1+l*D)-from(1+l*D))*k(2,i)+(to(1+l*D)*from(2+l*D)-to(2+l*D)*from(1+l*D)))/sqrt(sum((from(l*D+1:l*D+2)-to(l*D+1:l*D+2))**2))<1e-6_wp) then
						if(sum((from(l*D+1:l*D+2)-k(:,i))*(to(l*D+1:l*D+2)-k(:,i)))>0._wp) cycle
						m(2)=m(2)+1
						idx_(m(2))=i
						if(abs(to(l*D+1)-from(l*D+1))>abs(to(l*D+2)-from(l*D+2))) then
							k_(m(2))=k(1,i)*sign(1._wp,to(l*D+1)-from(l*D+1))
						else
							k_(m(2))=k(2,i)*sign(1._wp,to(l*D+2)-from(l*D+2))
						endif
					endif
				enddo
				if(m(2)-m(1)<2) then
					write(*,*)"line",l+1," has only ",m(2)-m(1),"point, try to increase distance"
				endif
				!call qsort(k_(m(1)+1:m(2)),idx_(m(1)+1:m(2)))
				ord(1:m(2)-m(1))=[1:m(2)-m(1)]
				call qsort(k_(m(1)+1:m(2)),ord(1:m(2)-m(1)))
				idx_(m(1)+1:m(2))=idx_(m(1)+ord(1:m(2)-m(1)))
				m(1)=m(2)
				idx_(size(idx_)-l-1)=m(1)
			enddo
			allocate(idx(m(2)+sl+1))
			idx(1:m(2))=idx_(1:m(2))
			idx(size(idx)-sl:size(idx))=-idx_(size(idx_)-sl:size(idx_))
		end subroutine
	end function
	subroutine solver(Delta,mu,Ef,Gc,ife,nf,db,G4c)
#ifdef cluster
		complex(wp) :: Delta(:,:,:),Gc(:,:,:)
#else
		complex(wp) :: Delta(:,:),Gc(:,:)
#endif
		real(wp) :: Ef,mu
		complex(wp) :: ife
		complex(wp), optional :: G4c(:,:)
		real(wp), optional :: db(:),nf(:)
		complex(wp) :: nf_(Nc),db_(Nc),iAc(Nc,Nc),cNc(Nc),iA(n,Nc,Nc),cn(n),Z0,Z,lnd,lndj
		logical :: cfg(Nc)
		integer :: iw,jw,i,j,dnf
#ifndef cluster
		Z0=0._wp
		Z=0._wp
		do i=1,Nc
			Z0=Z0+(msum(log((iomega(1:n)+mu-Delta(1:n,i))/iomega(1:n)),not_conv=underscore))
			lnd=(-beta*(Ef-mu)+msum(log((iomega(1:n)+mu-Delta(1:n,i)-U)/(iomega(1:n)+mu-Delta(1:n,i))),not_conv=underscore))
			Z=Z+log(1._wp+exp(lnd))
			if(present(nf)) then
				nf(i)=exp(lnd)/(1._wp+exp(lnd))
			endif
			if(present(db)) then
				db(i)=Tk*msum(1._wp/(iomega(1:n)+mu-Delta(1:n,i)-U),not_conv=underscore)*exp(lnd)/(1._wp+exp(lnd))
			endif

			!$OMP PARALLEL DO
			do iw=1,n
				Gc(iw,i)=(1._wp/(iomega(iw)+mu-Delta(iw,i))+1._wp/(iomega(iw)+mu-Delta(iw,i)-U)*exp(lnd))/(1._wp+exp(lnd))
			enddo
			!$OMP END PARALLEL DO
			if(present(G4c)) then
				do j=1,Nc
					if(i/=j) then
						lndj=(-beta*(Ef-mu)+msum(log((iomega(1:n)+mu-Delta(1:n,j)-U)/(iomega(1:n)+mu-Delta(1:n,j))),not_conv=underscore))
					endif
					!$OMP PARALLEL DO
					do iw=1,n
						do jw=1,n
							if(i==j) then
								G4c(iw+(i-1)*n,jw+(j-1)*n)=(1._wp/(iomega(iw)+mu-Delta(iw,i))*1._wp/(iomega(jw)+mu-Delta(jw,j))+1._wp/(iomega(iw)+mu-Delta(iw,i)-U)*1._wp/(iomega(jw)+mu-Delta(jw,j)-U)*exp(lnd))/(1._wp+exp(lnd))
							else
								G4c(iw+(i-1)*n,jw+(j-1)*n)=(1._wp/(iomega(iw)+mu-Delta(iw,i))+1._wp/(iomega(iw)+mu-Delta(iw,i)-U)*exp(lnd))*(1._wp/(iomega(jw)+mu-Delta(jw,j))+1._wp/(iomega(jw)+mu-Delta(jw,j)-U)*exp(lndj))/((1._wp+exp(lnd))*(1._wp+exp(lndj)))
							endif
						enddo
					enddo
					!$OMP END PARALLEL DO
				enddo
			endif
		enddo
		ife=-Tk*Z0/Nc-Tk*log(2._wp)-Tk*Z/Nc+mu*n_c-(Ef-mu)*n_f
#else

		cfg=.false.
		!$OMP PARALLEL DO PRIVATE(iAc,cNc)
		do iw=1,n
			iA(iw,1:Nc,1:Nc)=diag(iomega(iw)+mu,Nc)-Delta(iw,1:Nc,1:Nc)-diag(merge(U,0._wp,cfg))
			iAc(1:Nc,1:Nc)=matmul(iA(iw,1:Nc,1:Nc),diag(1._wp/iomega(iw),Nc))
			!iAc(1:Nc,1:Nc)=matmul(iA(iw,1:Nc,1:Nc),diag(1._wp/(iomega(iw)+mu),Nc))

			call geev(iAc(1:Nc,1:Nc),cNc)
			cn(iw)=sum(log(cNc))

			call mat_inv(iA(iw,1:Nc,1:Nc))
		enddo
		!$OMP END PARALLEL DO
		Z0=-Tk*msum(cn,not_conv=underscore)/Nc-Tk*log(2._wp)
		!Z0=-Tk*msum(cn,not_conv=underscore)/Nc-Tk*log(1._wp+exp(-mu/Tk))-mu
		Gc(1:n,1:Nc,1:Nc)=iA(1:n,1:Nc,1:Nc)
		!do iw=1,n
			!do iQ=1,Nc/Ns
				!do j=1,Ns
					!do i=1,Ns
						!!self%Gs(iQ,iw,i,j)=sum(lexp(iQ,:)*matmul(iA(iw,i:Nc:Ns,j:Nc:Ns),rexp(iQ,:)))
						!self%Gs(iQ,iw,i,j)=sum(lexp(iQ,:)*iA(iw,i:Nc:Ns,j))
					!enddo
				!enddo
			!enddo
		!enddo
		if(present(G4c)) then
			! segment fault if parallel out loop, seems to be a bug
			!$OMP PARALLEL DO
			do iw=1,n
				do jw=1,n
					do i=1,Nc**2
						do j=1,Nc*Ns
							G4c(iw+(i-1)*n,jw+(j-1)*n)=iA(iw,mod(i-1,Nc)+1,(i-1)/Nc+1)*iA(jw,mod(j-1,Nc)+1,(j-1)/Nc+1)
						enddo
					enddo
				enddo
			enddo
			!$OMP END PARALLEL DO
		endif
		!write(*,*)iA(n/2,1,1),self%Gc(n/2,1,1)
		nf_=merge(1._wp,0._wp,cfg)
		Z=0._wp
		dnf=0
		lnd=0._wp

		do i=1,size(cfg)
			call gen_cfg(i)
		enddo
		ife=Z0-Tk*Z/Nc+mu*n_c-(Ef-mu)*n_f

		if(present(nf)) then
			nf=real(nf_)
		endif
		if(present(db)) then
			db=real(db_)
		endif
	contains
		recursive subroutine gen_cfg(m)
			integer :: m
			integer :: k,i,j,iw,jw,NcNs
			real(8) :: sg
			complex(8), save :: tmp(Nc),exp_,norm_
			NcNs=Nc/Ns
			cfg(m)=.not.cfg(m)
			sg=merge(1._wp,-1._wp,cfg(m))
			dnf=dnf+nint(sg)
			lnd=lnd+msum(log((1._wp+iA(1:n,m,m)*(-sg)*U)),not_conv=underscore)
			if(real(-beta*(Ef-mu)*dnf+lnd)>real(Z)) then
				exp_=exp(-(-beta*(Ef-mu)*dnf+lnd-Z))
				norm_=exp_/(1._wp+exp_)
				Z=(-beta*(Ef-mu)*dnf+lnd)+log(1._wp+exp_)
				exp_=1._wp/(1._wp+exp_)
			else
				exp_=exp(-beta*(Ef-mu)*dnf+lnd-Z)
				norm_=1._wp/(1._wp+exp_)
				Z=Z+log(1._wp+exp_)
				exp_=exp_*norm_
			endif
			!write(*,*)cfg,Z,exp_
			!$OMP PARALLEL DO PRIVATE(tmp)
			do iw=1,n
				tmp(:)=iA(iw,:,m)
				do j=1,Nc
					iA(iw,:,j)=iA(iw,:,j)-iA(iw,m,j)*tmp(:)*U/(tmp(m)*U-sg)
					Gc(iw,:,j)=Gc(iw,:,j)*norm_+iA(iw,:,j)*exp_
				enddo
				!do iQ=1,NcNs
					!do i=1,Ns
						!self%Gs(iQ,iw,i,1:Ns)=self%Gs(iQ,iw,i,1:Ns)*norm_+matmul(lexp(iQ,:),iA(iw,i:Nc:Ns,1:Ns))*exp_
					!enddo
				!enddo
			enddo
			!$OMP END PARALLEL DO
			!write(*,*)sum([(sum(abs(matmul(iA(i,1:Nc,1:Nc),diag(iomega(i)+mu,Nc)-self%Delta(i,1:Nc,1:Nc)-diag(merge(U,0._wp,cfg)))-diag(1._wp,Nc))),i=1,n)])
			do j=1,Nc
				nf_(j)=nf_(j)*norm_+merge(exp_,cmplx(0._wp,kind=wp),cfg(j))
			enddo
			if(present(db)) then
				do j=1,Nc
					if(cfg(j)) then
						db_(j)=db_(j)*norm_+Tk*msum(iA(1:n,j,j),not_conv=underscore)*exp_
					else
						db_(j)=db_(j)*norm_
					endif
				enddo
			endif
			if(present(G4c)) then
				!$OMP PARALLEL DO
				do iw=1,n
					do jw=1,n
						do i=1,Nc**2
							do j=1,Nc*Ns
								G4c(iw+(i-1)*n,jw+(j-1)*n)=G4c(iw+(i-1)*n,jw+(j-1)*n)*norm_+&
											iA(iw,mod(i-1,Nc)+1,(i-1)/Nc+1)*iA(jw,mod(j-1,Nc)+1,(j-1)/Nc+1)&
											*exp_
							enddo
						enddo
					enddo
				enddo
				!$OMP END PARALLEL DO
			endif
			do k=1,m-1
				call gen_cfg(k)
			enddo
		end subroutine
#endif
	end subroutine
	function recip(r) result(rt)
		real(wp) :: r(:,:)
		real(wp) :: rt(size(r,1),size(r,2))
		integer :: i,j,D
		D=size(r,1)
		do i=1,D
			do j=1,D
				rt(j,i)=(-1._wp)**(j+1)*det(r(pack([1:D],[1:D]/=j),pack([1:D],[1:D]/=i)))
			enddo
			rt(:,i)=2d0*pi*rt(:,i)/sum(r(:,i)*rt(:,i))
		enddo
	end function
	function fit(y,x,nx) result(rt)
		complex(wp) :: y(:),rt(size(y))
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
	subroutine io(self,ut,flag)
		class(gfs) :: self
		integer :: ut(:),flag(:)
		integer :: iw,jw,k,ik,i,j,l,iQ
		complex(wp) :: Gkl(n),Gc(n),tmp(Ns,3),fe_,iA_(Ns,Ns)
		real(wp) :: val,dr(D)
		logical :: isopen
		if(this_image()==1) then
		write(*,*)"***********************"
		do k=1,size(flag)
			select case(flag(k))
			case(out_T)
				val=abs(real(sum((self%n_f(1:Ns)-sum(self%n_f)/Nc)*exp(img*matmul(self%Qs(:,1),self%rs(:,1:Ns))))))
				iQ=1
				do i=1,Ns
					write(*,"(*(es12.4))")sum((self%n_f(1:Ns)-sum(self%n_f)/Nc)*exp(img*matmul(self%Qs(:,i),self%rs(:,1:Ns))))
					if(abs(real(sum((self%n_f(1:Ns)-sum(self%n_f)/Nc)*exp(img*matmul(self%Qs(:,i),self%rs(:,1:Ns))))))>val) then
						val=abs(real(sum((self%n_f(1:Ns)-sum(self%n_f)/Nc)*exp(img*matmul(self%Qs(:,i),self%rs(:,1:Ns))))))
						iQ=i
					endif
				enddo
				do i=1,Ns
					if(abs(real(sum((self%n_c(1:Ns)-sum(self%n_c)/Nc)*exp(img*matmul(self%Qs(:,i),self%rs(:,1:Ns))))))>abs(real(sum((self%n_c(1:Ns)-sum(self%n_c)/Nc)*exp(img*matmul(self%Qs(:,iQ),self%rs(:,1:Ns))))))) then
						write(*,*)"warning, iQ"
					endif
				enddo
				if(.not.is_in(self%Qs(:,iQ),Ta,dr)) then
				endif
				write(*,"(es12.4$)")Tk,U,self%Ef,self%mu,sum(abs(self%n_c(1:Ns)-sum(self%n_c)/Nc))/Ns,sum(abs(self%n_f(1:Ns)-sum(self%n_f)/Nc))/Ns,(self%Qs(:,iQ)+dr)/pi,self%fe,self%ife
				write(*,"(l3)")self%is_symbk


				if(.not.self%conv) then
					write(ut(k),"(A$)")"#"
				endif
				write(ut(k),"(e29.20e3$)")Tk,U,self%Ef,self%mu,sum(abs(self%n_c(1:Ns)-sum(self%n_c)/Nc))/Ns,sum(abs(self%n_f(1:Ns)-sum(self%n_f)/Nc))/Ns,sum(self%n_c)/Nc,sum(self%n_cimp)/Nc,sum(self%n_f)/Nc,(self%Qs(:,iQ)+dr)/pi,self%fe,self%ife,sum(self%db)/Nc,self%fsd,merge(maxval(real(self%Sus)),minval(real(self%Sus)),all(real(self%Sus)>=0._wp))
				write(ut(k),"(l3)")self%is_symbk
			case(out_pattern)
				if(self%conv) then
					do i=1,Nc
						write(ut(k),"(*(e29.20e3))")Tk,U,self%rc(:,i),self%n_c(i),self%n_f(i)
					enddo
					write(ut(k),"(x)")
				endif
			case(out_w)
				write(ut(k),"(A)")"Tk omega rGcl iGcl rSEck iSEck rdSEc idSEc rGc iGc rDelta iDelta rSEc iSEc"
				!do iw=1,n
					!do ik=1,nkD
						!tmp(iw,1:Ns,1:Ns)=self%to_k(self%Gkl(iw,1:Ns,1:Ns),self%k(1:D,ik))
						!Gkl(iw)=sum([(tmp(iw,i,i),i=1,Ns)])/Ns
					!enddo
				!enddo
				do i=1,Nc
					!write(*,*)integrate(omega,A=-imag(self%Gkl(:,i,i))/pi)
					do j=1,Nc
						do iw=1,n
							write(ut(k),"(*(e29.20e3))")Tk,omega(iw),self%Gcl(iw,i,j),sum(self%SEck(:,iw,i,j))/nkD,sum(self%dSEc(:,iw,i,j))/nkD,&
#ifdef cluster
								self%Gc(iw,i,j),self%Delta(iw,i,j),self%SEc(iw,i,j)
#else
								merge(self%Gc(iw,i),cmplx(0._wp,kind=wp),i==j),merge(self%Delta(iw,i),cmplx(0._wp,kind=wp),i==j),merge(self%SEc(iw,i),cmplx(0._wp,kind=wp),i==j)
#endif
						enddo
						write(ut(k),"(x)")
					enddo
				enddo
				write(ut(k),"(x)")
			!case(out_wk)
				!write(ut(k),"(A)")"k Tk omega kx ky rGk iGk rGk0 iGk0 ek rdVer idVer"
				!l=0
				!do ik=1,size(idx)
					!i=mod(idx(ik)-1,nkD)+1
					!j=(idx(ik)-1)/nkD+1
					!tmp(:,1)=diag(self%to_k(cmplx(self%eks(i,1:Ns,1:Ns),kind=wp),self%ks(1:D,i)))-self%mu
					!!write(*,"((i4,e29.20e3))")ik,l+(real(ik-abs(idx(size(idx)-l))))/(abs(idx(size(idx)-l-1))-abs(idx(size(idx)-l)))
					!do iw=1,n
						!tmp(:,2)=diag(self%to_k(cmplx(self%Gsk(i,iw,1:Ns,1:Ns),kind=wp),self%ks(1:D,i)))
						!!tmp(:,3)=diag(self%to_k(cmplx(self%Gk0(i,iw,1:Ns,1:Ns),kind=wp),self%k(1:D,i)))
						!write(ut(k),"(*(e29.20e3))")l+(real(ik-abs(idx(size(idx)-l))))/(abs(idx(size(idx)-l-1))-abs(idx(size(idx)-l))),&
							!Tk,omega(iw),self%k(:,idx(ik))/pi,tmp(j,2),tmp(j,3),tmp(j,1)!,self%dGk(idx(ik),iw),self%Gk(idx(ik),iw),self%dVer(self%uk2k(iorder,1),iw,iw)
					!enddo
					!write(ut(k),"(x)")
					!if(ik==abs(idx(size(idx)-l-1))) then
						!l=l+1
						!if(idx(size(idx)-l-1)>0) then
							!exit
						!endif
					!endif
				!enddo
				!write(ut(k),"(x)")
			case(out_kc)
				write(ut(k),"(A)")"Tk omega kx ky rdGc idGc"
				do i=1,Nc
					do j=1,Nc
						do ik=1,nkD
							write(ut(k),"(*(e29.20e3))")Tk,omega(n/2+1),self%kc(:,ik)/pi,self%dGc(ik,n/2+1,i,j)
						enddo
						write(ut(k),"(x)")
					enddo
				enddo
				write(ut(k),"(x)")
			case(out_k)
				if(self%conv) then
					write(*,*)compute_(dGk)
					do ik=1,nkD*Nc
						write(ut(k),"(*(e29.20e3))")Tk,omega(n/2+1),self%k(:,ik)/pi,self%Gk(n/2+1,ik),self%dGk(n/2+1,ik),self%dSEk(n/2+1,ik),self%fs(ik),self%Sus(ik),self%Sus0(ik)
					enddo
					write(ut(k),"(x)")
				endif
			case(out_r)
				if(self%conv) then
					do i=1,nkD*Nc
						write(ut(k),"(*(e29.20e3))")Tk,U,self%r(:,i),self%Susr(i)
					enddo
					write(ut(k),"(x)")
				endif
			case(inout_save)
				rewind(inout_save)
				write(inout_save)Tk,U,self%Delta,self%dSEc,self%mu,self%Ef
			case(out_lattice)
				write(*,"('nk:',i5,', Ns:',i5,', Nc:',i5,', n:',i5,', n_c:',es12.4,', n_f:',es12.4', et:',es12.4)")nk,Ns,Nc,n,n_c,n_f,et
				write(*,"('Cnv:',i5,', Cnvc:',i5,', sdiagram:',2i3,', ddiagram:',2i3)")Cnv,Cnvc,sdiagram,ddiagram
				write(*,"('a:    ',4es12.4)")a
				write(*,"('scell:',4es12.4)")scell
				write(*,"('clust:',4es12.4)")clust
				write(*,"('bcut:',4es12.4)")bcut
				write(*,"('is_real:',l2,', which_gen_cfg:',i2,', is_SE:',l2,', sc_scheme:',A,i3,', lse: ',i3,', is_CDMFT:',l2,', is_dF4_iter:',l2,', is_dual_pureDCA:',l2,', is_cluster:',l2)")&
					is_real,which_gen_cfg,is_SE,pack(["DMFT","dDMFT","dDMFT_all","dDMFT_simp","dDMFT_simp_nc"],[DMFT,dDMFT,dDMFT_all,dDMFT_simp,dDMFT_simp_nc]==scheme),sc_scheme,lse,is_CDMFT,is_dF4_iter,is_dual_pureDCA,&
#ifdef cluster
.true.
#else
.false.
#endif
				do i=10,1000
					inquire(unit=i, opened=isopen)
					if(isopen) then
						write(i,"('#nk:',i5,', Ns:',i5,', Nc:',i5,', n:',i5,', n_c:',es12.4,', n_f:',es12.4', et:',es12.4)")nk,Ns,Nc,n,n_c,n_f,et
						write(i,"('#Cnv:',i5,', Cnvc:',i5,', sdiagram:',2i3,', ddiagram:',2i3)")Cnv,Cnvc,sdiagram,ddiagram
						write(i,"('#a:    ',4es12.4)")a
						write(i,"('#scell:',4es12.4)")scell
						write(i,"('#clust:',4es12.4)")clust
						write(i,"('#bcut:',4es12.4)")bcut
						write(i,"('#is_real:',l2,', which_gen_cfg:',i2,', is_SE:',l2,', sc_scheme:',A,i3,', lse: ',i3,', is_CDMFT:',l2,', is_dF4_iter:',l2,', is_dual_pureDCA:',l2,', is_cluster:',l2)")&
							is_real,which_gen_cfg,is_SE,pack(["DMFT","dDMFT","dDMFT_all","dDMFT_simp","dDMFT_simp_nc"],[DMFT,dDMFT,dDMFT_all,dDMFT_simp,dDMFT_simp_nc]==scheme),sc_scheme,lse,is_CDMFT,is_dF4_iter,is_dual_pureDCA,&
#ifdef cluster
.true.
#else
.false.
#endif
					endif
				enddo

				write(out_T,"(A)")"T1 T2 "
				do i=1,size(Ta,2)
					write(out_T,"(*(e29.20e3))")Ta(:,i)/pi
				enddo
				write(out_T,"(*(e29.20e3))")Ta(:,1)/pi
				write(out_T,"(x)")
				do i=1,size(Ts,2)
					write(out_T,"(*(e29.20e3))")Ts(:,i)/pi
				enddo
				write(out_T,"(*(e29.20e3))")Ts(:,1)/pi
				write(out_T,"(x)")
				do i=1,size(Tc,2)
					write(out_T,"(*(e29.20e3))")Tc(:,i)/pi
				enddo
				write(out_T,"(*(e29.20e3))")Tc(:,1)/pi
				write(out_T,"(x/)")
				write(out_T,"(A)")"T U ef mu tnc tnf nc ncimp nf Qx Qy fe _ ife _ db fsd Sus symbrk"

				write(out_pattern,"(A)")"x y"
				write(out_pattern,"(*(e29.20e3))")0._wp,0._wp
				write(out_pattern,"(*(e29.20e3))")scell(:,1)
				write(out_pattern,"(*(e29.20e3))")scell(:,1)+scell(:,2)
				write(out_pattern,"(*(e29.20e3))")scell(:,2)
				write(out_pattern,"(*(e29.20e3))")0._wp,0._wp
				write(out_pattern,"(x)")
				write(out_pattern,"(*(e29.20e3))")0._wp,0._wp
				write(out_pattern,"(*(e29.20e3))")clust(:,1)
				write(out_pattern,"(*(e29.20e3))")clust(:,1)+clust(:,2)
				write(out_pattern,"(*(e29.20e3))")clust(:,2)
				write(out_pattern,"(*(e29.20e3))")0._wp,0._wp
				write(out_pattern,"(x/)")
				write(out_pattern,"(A)")"T U x y nc nf"

				write(out_k,"(A)")"x y"
				do i=1,size(Ta,2)
					write(out_k,"(*(e29.20e3))")Ta(:,i)/pi
				enddo
				write(out_k,"(*(e29.20e3))")Ta(:,1)/pi
				write(out_k,"(x)")
				do i=1,size(Ts,2)
					write(out_k,"(*(e29.20e3))")Ts(:,i)/pi
				enddo
				write(out_k,"(*(e29.20e3))")Ts(:,1)/pi
				write(out_k,"(x)")
				do i=1,size(Tc,2)
					write(out_k,"(*(e29.20e3))")Tc(:,i)/pi
				enddo
				write(out_k,"(*(e29.20e3))")Tc(:,1)/pi
				write(out_k,"(x/)")
				write(out_k,"(A)")"Tk omega kx ky rGk iGk rdGk idGk rdSEk idSEk fs Sus _ Sus0 _"

				write(out_r,"(A)")"Tk U rx ry Sus _"

				! momentum space
				write(ut(k),"(A)")"x0 y0 b1 b2 "
				write(ut(k),"(*(e29.20e3))")0._wp,0._wp,ba(:,1)/pi
				write(ut(k),"(*(e29.20e3))")0._wp,0._wp,ba(:,2)/pi
				write(ut(k),"(x)")
				write(ut(k),"(*(e29.20e3))")0._wp,0._wp,bs(:,1)/pi
				write(ut(k),"(*(e29.20e3))")0._wp,0._wp,bs(:,2)/pi
				write(ut(k),"(x)")
				write(ut(k),"(*(e29.20e3))")0._wp,0._wp,bc(:,1)/pi
				write(ut(k),"(*(e29.20e3))")0._wp,0._wp,bc(:,2)/pi
				write(ut(k),"(x/)")

				write(ut(k),"(A)")"T1 T2 "
				do i=1,size(Ta,2)
					write(ut(k),"(*(e29.20e3))")Ta(:,i)/pi
				enddo
				write(ut(k),"(*(e29.20e3))")Ta(:,1)/pi
				write(ut(k),"(x)")
				do i=1,size(Ts,2)
					write(ut(k),"(*(e29.20e3))")Ts(:,i)/pi
				enddo
				write(ut(k),"(*(e29.20e3))")Ts(:,1)/pi
				write(ut(k),"(x)")
				do i=1,size(Tc,2)
					write(ut(k),"(*(e29.20e3))")Tc(:,i)/pi
				enddo
				write(ut(k),"(*(e29.20e3))")Tc(:,1)/pi
				write(ut(k),"(x/)")
				write(ut(k),"(A)")"k1 k2"
				do ik=1,size(self%k,2)
					write(ut(k),"(*(e29.20e3))")self%k(:,ik)/pi
				enddo
				write(ut(k),"(x)")
				do ik=1,size(self%ks,2)
					write(ut(k),"(*(e29.20e3))")self%ks(:,ik)/pi
				enddo
				write(ut(k),"(x)")
				do ik=1,size(self%kc,2)
					write(ut(k),"(*(e29.20e3))")self%kc(:,ik)/pi
				enddo
				write(ut(k),"(x)")
				do ik=1,size(self%Qs,2)
					write(ut(k),"(*(e29.20e3))")self%Qs(:,ik)/pi
				enddo
				write(ut(k),"(x)")
				do ik=1,size(self%Qc,2)
					write(ut(k),"(*(e29.20e3))")self%Qc(:,ik)/pi
				enddo
				write(ut(k),"(x)")
				do ik=1,size(self%uk2k,1)
					write(ut(k),"(e29.20e3$)")self%k(:,self%uk2k(ik,1))/pi
					write(ut(k),"(i3)")self%uk2k(ik,0)
				enddo
				write(ut(k),"(x/)")
				! real space
				write(ut(k),"(A)")"x0 y0 b1 b2 "
				write(ut(k),"(*(e29.20e3))")0._wp,0._wp,a(:,1)
				write(ut(k),"(*(e29.20e3))")0._wp,0._wp,a(:,2)
				write(ut(k),"(x/)")
				write(ut(k),"(A)")"T1 T2 "
				write(ut(k),"(*(e29.20e3))")[1:D]*0._wp
				write(ut(k),"(*(e29.20e3))")scell(:,1)
				write(ut(k),"(*(e29.20e3))")scell(:,1)+scell(:,2)
				write(ut(k),"(*(e29.20e3))")scell(:,2)
				write(ut(k),"(*(e29.20e3))")[1:D]*0._wp
				write(ut(k),"(x)")
				write(ut(k),"(*(e29.20e3))")[1:D]*0._wp
				write(ut(k),"(*(e29.20e3))")clust(:,1)
				write(ut(k),"(*(e29.20e3))")clust(:,1)+clust(:,2)
				write(ut(k),"(*(e29.20e3))")clust(:,2)
				write(ut(k),"(*(e29.20e3))")[1:D]*0._wp
				write(ut(k),"(x/)")
				write(ut(k),"(A)")"r1 r2"
				do i=1,Ns
					write(ut(k),"(*(e29.20e3))")self%rs(:,i)
				enddo
				write(ut(k),"(x)")
				do i=1,Nc
					write(ut(k),"(*(e29.20e3))")self%rc(:,i)
				enddo
				write(ut(k),"(x)")
			end select
		end do
		write(*,*)"***********************"
		endif
		sync all
	end subroutine
	recursive function is_in(r,T,dr)
		logical :: is_in
		real(wp) :: r(2),T(:,:)
		real(wp), optional :: dr(:)
		real(wp) :: tr(2),EG(2),DT(2),tf(2,2),rerr(2),r_(2),err=1e-6_wp
		integer :: i,j,n
		r_=r
		if(present(dr)) then
			dr=0._wp
		endif
		is_in=.true.
		rerr=0._wp
		n=size(T,2)
		do i=1,n/2
			rerr=rerr+(T(:,i+1)-T(:,i))*err*i/1.34_wp
		enddo
		do 
			do i=1,size(T,2)
				j=mod(i,size(T,2))+1
				EG=T(:,j)-T(:,i)
				DT=(T(:,mod(j-1+n/2,n)+1)-T(:,j))+EG
				tf=reshape((/EG(2),-DT(2),-EG(1),DT(1)/)/(DT(1)*EG(2)-DT(2)*EG(1)),(/2,2/))
				tr=matmul(tf,r_(:2)-(T(:,j)+T(:,i))/2._wp+rerr(:2))
				if(tr(1)<0._wp) then
					is_in=.false.
					if(present(dr)) then
						r_=r_-floor(tr(1))*DT
					else
						return
					endif
				endif
			enddo
			if(present(dr)) then
				if(is_in(r_,T)) then
					dr=r_-r
					exit
				endif
			else
				exit
			endif
		enddo
	end function
	subroutine get_dF4(G4c,dSus0,dF4,order)
		complex(wp) :: G4c(:,:),dSus0(:,:),dF4(:,:)
		integer :: order(2)
		complex(wp), allocatable :: x(:),x_(:)!,V4Sus(:,:)
#ifdef cluster
		complex(wp) :: V4Sus(n*Nc*Nc,n*Nc*dNs)
#else
		complex(wp) :: V4Sus(n*Nc,n*Nc)
#endif
		complex(wp) :: ctmp!,V4Sus(n*Nc**2,n*Nc*dNs)
		real(wp) :: maxv
		integer :: iw,jw,i
#ifndef cluster
		!allocate(V4Sus(n*Nc,n*Nc))
		dF4(1:n*Nc,1:n*Nc)=G4c(1:n*Nc,1:n*Nc)
		if(order(2)<0) then
			do jw=1,n*Nc
				V4Sus(1:n*Nc,jw)=matmul(G4c(1:n*Nc,mod(jw-1,n)+1:n*Nc:n),dSus0(mod(jw-1,n)+1:n*Nc:n,(jw-1)/n+1))*Res(mod(jw-1,n)+1)
				V4Sus(jw,jw)=1._wp+V4Sus(jw,jw)
			enddo
			call gesv(V4Sus,dF4(1:n*Nc,1:n*Nc))
		endif
		if(order(1)>1.or.order(2)>1) then
			do jw=1,n*Nc
				V4Sus(1:n*Nc,jw)=-matmul(G4c(1:n*Nc,mod(jw-1,n)+1:n*Nc:n),dSus0(mod(jw-1,n)+1:n*Nc:n,(jw-1)/n+1))*Res(mod(jw-1,n)+1)
			enddo
		endif
		do i=2,order(2)-order(1)+1
			dF4(1:n*Nc,1:n*Nc)=G4c(1:n*Nc,1:n*Nc)+matmul(V4Sus(1:n*Nc,1:n*Nc),dF4(1:n*Nc,1:n*Nc))
		enddo
		do i=2,order(1)
			dF4(1:n*Nc,1:n*Nc)=matmul(V4Sus(1:n*Nc,1:n*Nc),dF4(1:n*Nc,1:n*Nc))
		enddo
#else
		!allocate(V4Sus(n*Nc**2,n*Nc*dNs))
		if(.not.is_dF4_iter) then
			dF4(1:n*Nc**2,1:n*Nc*dNs)=G4c(1:n*Nc**2,1:n*Nc*dNs)
			do jw=1,n*Nc*dNs
				V4Sus(1:n*Nc**2,jw)=matmul(G4c(1:n*Nc**2,mod(jw-1,n)+1:n*Nc*dNs:n),dSus0(mod(jw-1,n)+1:n*Nc*dNs:n,(jw-1)/n+1))*Res(mod(jw-1,n)+1)
				V4Sus(jw,jw)=1._wp+V4Sus(jw,jw)
			enddo
			call gesv(V4Sus,dF4(1:n*Nc**2,1:n*Nc*dNs))
		else
			allocate(x(n**2*Nc**3*dNs),x_(n**2*Nc**3*dNs))
			!$OMP PARALLEL DO
			do jw=1,n*Nc*dNs
					x((jw-1)*n*Nc**2+1:jw*n*Nc**2)=dF4(1:n*Nc**2,jw)
					V4Sus(1:n*Nc**2,jw)=-matmul(G4c(1:n*Nc**2,mod(jw-1,n)+1:n*Nc*dNs:n),dSus0(mod(jw-1,n)+1:n*Nc*dNs:n,(jw-1)/n+1))*Res(mod(jw-1,n)+1)
			enddo
			!$OMP END PARALLEL DO
			call mbroyden(0,x,x_,0.2_wp,10,id=3)
			do i=1,1000
				maxv=0._wp
				if(this_image()==1) call start_time(id=1)
				!$OMP PARALLEL DO
				do jw=1,n*Nc*dNs
					x_((jw-1)*n*Nc**2+1:jw*n*Nc**2)=G4c(1:n*Nc**2,jw)&
						+matmul(V4Sus(:,:),x((jw-1)*n*Nc**2+1:jw*n*Nc**2))-x((jw-1)*n*Nc**2+1:jw*n*Nc**2)
				enddo
				!$OMP END PARALLEL DO
				if(this_image()==1) then
					write(*,"(A$)")"dF4 of "
					call stop_time(id=1,show=underscore)
				endif
				maxv=maxval(abs(x_))
				write(*,*)maxv
				if(maxv<1e-6_wp) then
					exit
				endif
				!x=x_+x
				if(this_image()==1) call start_time(id=2)
				call mbroyden(i,x,x_,id=3)
				if(this_image()==1) then
					write(*,"(A$)")"broyden of "
					call stop_time(id=2,show=underscore)
				endif
			enddo
			dF4=reshape(x(1:n**2*Nc**4),[n*Nc**2,n*Nc*dNs])
			call mbroyden(huge(1),x,x_,id=3)
		endif

#endif
	end subroutine
	function to_k(self,A,k) result(rt)
		class(gfs) :: self
		complex(wp) :: A(:,:)
		real(wp) :: k(:)
		complex(wp) :: rt(size(A,1),size(A,2))
		integer :: i,j,l
		do i=1,size(A,1)
			do j=1,size(A,2)
				rt(i,j)=0._wp
				do l=1,size(A,1)
					rt(i,j)=rt(i,j)+sum(A(l,:)*exp(-img*(sum((k+self%Qs(:,i))*self%rs(:,l))-matmul(k+self%Qs(:,j),self%rs))))
				enddo
				rt(i,j)=rt(i,j)/Ns
			enddo
		enddo
	end function
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
	complex(wp) function msum(G,not_conv) result(rt)
		complex(wp) :: G(:)
		logical, optional :: not_conv
		complex(wp) :: G1,G2,d
		real(wp) :: ipi_
		integer :: i,rg_(2),i1,i2
		if(is_real) then
			rt=integrate(omega,A=-imag(G)*ff(beta*omega))/(Tk*pi)
		else
			if(present(not_conv).and.(.not.is_normal_mat)) then
				rt=0.25_wp*(G(n)*img*omega(n)+G(1)*img*omega(1))
				!rt=0.25_wp*(G(n-1)*img*omega(n-1)+G(2)*img*omega(2))
			else
				rt=0._wp
			endif
			rt=rt+sum(-Res*G)
			!rt=rt+sum(-Res(2:n-1)*G(2:n-1))
		endif
	end function
	subroutine gen_graph(scheme)
		integer :: scheme
		graph=0
		inode=0
		onode=0
		update=.false.
		inode(mat_omega,1:2)=[omega_grid,PMT]
		inode(Gc,1:4)=[Delta,mu,Ef,mat_omega]
		inode(dGc0,1:2)=[Delta,Gc]
		inode(dSus0,1:2)=[dGc,dGcp]
		inode(V4c,1:3)=[Delta,mu,Ef]
		inode(Sus,1:5)=[Gc,dGc,dGcp,V4c,dSus0]
		if(lse==0) then
			inode(SEck,1:2)=[Gck_l,mu]
			inode(Gck_l,1:3)=[Delta,dGc,dGcp]
		else
			inode(SEck,1:3)=[SEc,Gc,dSEc]
			inode(Gck_l,1:2)=[mu,SEck]
		endif
		inode(inc,1:1)=[Gck_l]
		inode(fsd,1:1)=[Gck_l]
		inode(Gk,1:1)=[Gck_l]
		inode(fs,1:1)=[Gk]
		inode(SEc,1:3)=[Delta,Gc,mu]
		if((.not.fix_cp).and.split_cp) then
			inode(Ef,1:1)=[SC_cp]
			inode(mu,1:1)=[SC_cp]
			inode(SC_cp,1:4)=[inc,Gc,Ef,mu]
		endif
		select case(scheme)
		case(DMFT)
			inode(dGcp,1:1)=[dGc0]
			inode(Delta,1:1)=[SC_DMFT]
			if(fix_cp) then
				!inode(SC_DMFT,1:4)=[Gck_l,dGcp,Gc,Delta]
				!underscore=evaluate(graph=[Sus,V4c,dSus0,fs,fsd,Gk,[SC_DMFT,inc,Gck_l,dGcp,dGc0,Gc,Delta],dSEc,mu,Ef,mat_omega,lattice,omega_grid,PMT])

				!inode(Ef,1:1)=[SC_DMFT]
				!inode(SC_DMFT,1:5)=[Gck_l,dGcp,Gc,Delta,Ef]

				inode(SC_DMFT,1:4)=[Gck_l,dGcp,Gc,Delta]
				if(lse==0) then
					!underscore=evaluate(graph=[SEck,SEc,Sus,V4c,dSus0,fs,fsd,Gk,[SC_DMFT,inc,Gck_l,dGcp,dGc0,Gc,Delta,Ef],dSEc,mu,mat_omega,lattice,omega_grid,PMT])
					underscore=evaluate(graph=[SEck,SEc,Sus,V4c,dSus0,fs,fsd,Gk,[SC_DMFT,inc,Gck_l,dGcp,dGc0,Gc,Delta],Ef,dSEc,mu,mat_omega,lattice,omega_grid,PMT])
				else
					underscore=evaluate(graph=[Sus,V4c,dSus0,fs,fsd,Gk,[SC_DMFT,inc,Gck_l,SEck,dGcp,dGc0,SEc,Gc,Delta,Ef],dSEc,mu,mat_omega,lattice,omega_grid,PMT])
				endif
			else
				if(split_cp) then
					inode(SC_DMFT,1:4)=[Gck_l,dGcp,Gc,Delta]
					underscore=evaluate(graph=[Sus,V4c,dSus0,fs,fsd,Gk,[SC_DMFT,[SC_cp,inc,Gck_l,dGcp,dGc0,Gc,Ef,mu],Delta],dSEc,mat_omega,lattice,omega_grid,PMT])
				else
					inode(mu,1:1)=[SC_DMFT]
					inode(Ef,1:1)=[SC_DMFT]
					inode(SC_DMFT,1:7)=[inc,Gck_l,dGcp,Gc,Delta,mu,Ef]
					underscore=evaluate(graph=[Sus,V4c,dSus0,fs,fsd,Gk,[SC_DMFT,inc,Gck_l,dGcp,dGc0,Gc,Delta,mu,Ef],dSEc,mat_omega,lattice,omega_grid,PMT])
				endif
			endif
		case(dDMFT)
			inode(dGc,1:2)=[dGc0,dSEc]
			inode(Delta,1:1)=[SC_dDMFT]
			inode(dSEc,1:1)=[SC_dDMFT]
			if(fix_cp) then
				!inode(SC_dDMFT,1:7)=[Gck_l,Gc,V4c,dSus0,dGc,dSEc,Delta]
				!underscore=evaluate(graph=[Sus,fs,fsd,Gk,[SC_dDMFT,inc,Gck_l,V4c,dSus0,dGc,dGc0,Gc,dSEc,Delta],mu,Ef,mat_omega,lattice,omega_grid,PMT])

				!inode(Ef,1:1)=[SC_dDMFT]
				!inode(SC_dDMFT,1:9)=[inc,Gck_l,Gc,V4c,dSus0,dGc,dSEc,Delta,Ef]

				inode(SC_dDMFT,1:8)=[inc,Gck_l,Gc,V4c,dSus0,dGc,dSEc,Delta]
				if(lse==0) then
					!underscore=evaluate(graph=[SEc,SEck,Sus,fs,fsd,Gk,[SC_dDMFT,inc,Gck_l,V4c,dSus0,dGc,dGc0,Gc,dSEc,Delta,Ef],mu,lattice,mat_omega,omega_grid,PMT])
					underscore=evaluate(graph=[SEc,SEck,Sus,fs,fsd,Gk,[SC_dDMFT,inc,Gck_l,V4c,dSus0,dGc,dGc0,Gc,dSEc,Delta],Ef,mu,lattice,mat_omega,omega_grid,PMT])
				else
					underscore=evaluate(graph=[Sus,fs,fsd,Gk,[SC_dDMFT,inc,Gck_l,SEck,V4c,dSus0,dGc,dGc0,SEc,Gc,dSEc,Delta,Ef],mu,lattice,mat_omega,omega_grid,PMT])
				endif
			else
				if(split_cp) then
					inode(SC_dDMFT,1:7)=[Gck_l,Gc,V4c,dSus0,dGc,dSEc,Delta]
					if(lse==0) then
						underscore=evaluate(graph=[SEc,SEck,Sus,fs,fsd,Gk,[SC_dDMFT,V4c,dSus0,[SC_cp,inc,Gck_l,dGc,dGc0,Gc,Ef,mu],dSEc,Delta],mat_omega,lattice,omega_grid,PMT])
					else
						underscore=evaluate(graph=[Sus,fs,fsd,Gk,[SC_dDMFT,V4c,dSus0,[SC_cp,inc,Gck_l,SEck,dGc,dGc0,SEc,Gc,Ef,mu],dSEc,Delta],mat_omega,lattice,omega_grid,PMT])
					endif
				else
					inode(mu,1:1)=[SC_dDMFT]
					inode(Ef,1:1)=[SC_dDMFT]
					inode(SC_dDMFT,1:10)=[inc,Gck_l,Gc,V4c,dSus0,dGc,dSEc,Delta,mu,Ef]
					if(lse==0) then
						underscore=evaluate(graph=[SEc,SEck,Sus,fs,fsd,Gk,[SC_dDMFT,inc,Gck_l,V4c,dSus0,dGc,dGc0,Gc,dSEc,Delta,mu,Ef],lattice,mat_omega,omega_grid,PMT])
					else
						underscore=evaluate(graph=[Sus,fs,fsd,Gk,[SC_dDMFT,inc,Gck_l,SEck,V4c,dSus0,dGc,dGc0,SEc,Gc,dSEc,Delta,mu,Ef],lattice,mat_omega,omega_grid,PMT])
					endif
				endif
			endif
		case(dDMFT_simp)
			inode(dGc,1:2)=[dGc0,dSEc]
			inode(dSEc,1:1)=[SC_dDMFT]
			inode(SC_dDMFT,1:4)=[dSus0,V4c,dGc,dSEc]
			if(fix_cp) then
				!underscore=evaluate(graph=[Sus,fs,fsd,Gk,inc,Gck_l,[SC_dDMFT,dSus0,dGc,dSEc],V4c,dGc0,Gc,Delta,mu,Ef,lattice,mat_omega,omega_grid,PMT])

				inode(Ef,1:1)=[SC_dDMFT]
				inode(SC_dDMFT,1:8)=[inc,Gck_l,Gc,V4c,dSus0,dGc,dSEc,Ef]
				if(lse==0) then
					underscore=evaluate(graph=[SEc,SEck,Sus,fs,fsd,Gk,[SC_dDMFT,inc,Gck_l,dSus0,dGc,V4c,dGc0,Gc,dSEc,Ef],Delta,mu,lattice,mat_omega,omega_grid,PMT])
				else
					underscore=evaluate(graph=[Sus,fs,fsd,Gk,[SC_dDMFT,inc,Gck_l,SEck,dSus0,dGc,V4c,dGc0,SEc,Gc,dSEc,Ef],Delta,mu,lattice,mat_omega,omega_grid,PMT])
				endif
			else
				inode(Ef,:)=0
				inode(mu,:)=0
				inode(SC_cp,:)=0
				if(lse==0) then
					underscore=evaluate(graph=[SEc,SEck,Sus,fs,fsd,Gk,inc,Gck_l,[SC_dDMFT,dSus0,dGc,dSEc],V4c,dGc0,Gc,Ef,mu,Delta,mat_omega,lattice,omega_grid,PMT])
				else
					underscore=evaluate(graph=[Sus,fs,fsd,Gk,inc,Gck_l,SEck,[SC_dDMFT,dSus0,dGc,dSEc],V4c,dGc0,SEc,Gc,Ef,mu,Delta,mat_omega,lattice,omega_grid,PMT])
				endif
			endif
		case(dDMFT_simp_nc)
			inode(dGc,1:2)=[dGc0,dSEc]
			inode(dSEc,1:1)=[SC_dDMFT]
			if(split_cp) then
				inode(SC_dDMFT,1:4)=[dSus0,V4c,dGc,dSEc]
				if(lse==0) then
					underscore=evaluate(graph=[SEc,SEck,Sus,fs,fsd,Gk,inc,Gck_l,[SC_dDMFT,dSus0,V4c,[SC_cp,inc,Gck_l,dGc,dGc0,Gc,Ef,mu],dSEc],Delta,mat_omega,lattice,omega_grid,PMT])
				else
					underscore=evaluate(graph=[Sus,fs,fsd,Gk,inc,Gck_l,[SC_dDMFT,dSus0,V4c,[SC_cp,inc,Gck_l,SEck,dGc,dGc0,SEc,Gc,Ef,mu],dSEc],Delta,mat_omega,lattice,omega_grid,PMT])
				endif
			else
				inode(mu,1:1)=[SC_dDMFT]
				inode(Ef,1:1)=[SC_dDMFT]
				inode(SC_dDMFT,1:9)=[inc,Gck_l,dSus0,V4c,dGc,Gc,dSEc,mu,Ef]
				if(lse==0) then
					underscore=evaluate(graph=[SEc,SEck,Sus,fs,fsd,Gk,inc,Gck_l,[SC_dDMFT,inc,Gck_l,dSus0,V4c,dGc,dGc0,Gc,dSEc,Ef,mu],Delta,mat_omega,lattice,omega_grid,PMT])
				else
					underscore=evaluate(graph=[Sus,fs,fsd,Gk,inc,Gck_l,[SC_dDMFT,inc,Gck_l,SEck,dSus0,V4c,dGc,dGc0,SEc,Gc,dSEc,Ef,mu],Delta,mat_omega,lattice,omega_grid,PMT])
				endif
			endif
		end select
	end subroutine
end module
