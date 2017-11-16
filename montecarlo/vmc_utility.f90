include "model_utility.f90"
module vmc_utility
	use coarray
	use model_utility
	use lapack95, only : hegv, geev, hegv
	implicit none
	type t_mc
		integer :: sg,ne(2),hot,samp,step
		integer, allocatable :: cfg(:),Ecfg(:),nn(:,:)
		real(8), allocatable :: Emf(:)
		real(8) :: E=0d0,err=0d0
		type(t_ham) :: sphy,dphy
		complex(8), allocatable :: S(:,:),Ok2(:,:),Ek2(:,:),psi0(:),wf(:,:),dwf(:,:,:)
		real(8), allocatable :: g(:)
		integer :: num
		integer :: delay
		integer :: opt=2 ! 1: Tr(AdA)
						 ! 2: O=cicj
	contains
		procedure :: init
		procedure :: spect
		procedure :: do_mc
		procedure :: do_vmc
		procedure :: do_var
	end type
	real(8) :: dx=0d0
contains
	subroutine init(self,init_cfg)
		class(t_mc) :: self[ica1,*]
		logical :: init_cfg
		complex(8) :: H(Hmf%Hs,Hmf%Hs),cH(size(H,1),size(H,2)),D(size(H,1),size(H,2),sum(Hmf%var(1:)%n)),psi0_(Ns*Ns),Hq(Hmf%Hi)
		real(8) :: E(size(H,1)),dr(3),tmp
		integer :: l,i,j,n,n1,n2,nn_(0:Ns*Ns,2),kq,ord(Hmf%Hs),Nu
		logical :: check=.true.,is_mix
		if(is_ph) then
			change => change_ph
			get_row => get_row_ph
		else
			change => change_noph
			get_row => get_row_noph
		endif
		
		if(init_cfg) then
			if(allocated(self%cfg)) deallocate(self%cfg)
			allocate(self%cfg(Hmf%Hs))
			if(is_ph) then
				do i=1,size(self%ne)
					self%cfg(Hmf%s2H(1:self%ne(i),i))=[sum(self%ne(1:i))-self%ne(i)+1:sum(self%ne(1:i))]
					self%cfg(Hmf%s2H(self%ne(i)+1:Ns,i))=0
				enddo
			else
				self%cfg=0
				do i=1,size(self%ne)
					self%cfg(Hmf%s2H(sum(self%ne(1:i))-self%ne(i)+1:sum(self%ne(1:i)),i))=[sum(self%ne(1:i))-self%ne(i)+1:sum(self%ne(1:i))]
				enddo
			endif
		endif
		if(allocated(self%Ecfg)) deallocate(self%Ecfg)
		allocate(self%Ecfg(Hmf%Hs))
		if(allocated(self%Emf)) deallocate(self%Emf)
		allocate(self%Emf(Hmf%Hs))
		self%Ecfg=0

		H=0d0
		E=0d0
		ord=[1:size(ord)]
		is_mix=any([(Hmf%var(l)%c(1::2,:)%tp+Hmf%var(l)%c(2::2,:)%tp,l=Hmf%rg(1),Hmf%rg(2))]/=0)
		if(is_mix) then
			do i=1,brizon%nk
				call Hmf%Hamilton([Hmf%rg(1):Hmf%rg(2)],H((i-1)*Hmf%Hi+1:i*Hmf%Hi,(i-1)*Hmf%Hi+1:i*Hmf%Hi),brizon%k(i,:))
				call heev(H((i-1)*Hmf%Hi+1:i*Hmf%Hi,(i-1)*Hmf%Hi+1:i*Hmf%Hi),E((i-1)*Hmf%Hi+1:i*Hmf%Hi),"V")
			enddo
			call qsort(E,ord)
			self%Ecfg(ord(1:sum(self%ne)))=[1:sum(self%ne)]
			if(self%sg/=3) then
				if(abs(E(ord(sum(self%ne)))-E(ord(sum(self%ne)+1)))<1d-7) then
					if(this_image()==1) write(*,*)"warning, close-shell condition is not satisfied"
				endif
			endif
		else
			Nu=sum(latt%i2isb(latt%Ni,:),Hmf%mask(Hmf%tp2i(1),:))
			do i=1,brizon%nk
				call Hmf%Hamilton([Hmf%rg(1):Hmf%rg(2)],H((i-1)*Hmf%Hi+1:i*Hmf%Hi,(i-1)*Hmf%Hi+1:i*Hmf%Hi),brizon%k(i,:))
				call heev(H((i-1)*Hmf%Hi+1:(i-1)*Hmf%Hi+Nu,(i-1)*Hmf%Hi+1:(i-1)*Hmf%Hi+Nu),E((i-1)*Hmf%Hi+1:(i-1)*Hmf%Hi+Nu),"V")
				call heev(H((i-1)*Hmf%Hi+Nu+1:i*Hmf%Hi,(i-1)*Hmf%Hi+Nu+1:i*Hmf%Hi),E((i-1)*Hmf%Hi+Nu+1:i*Hmf%Hi),"V")

			enddo
			call qsort(E+(maxval(E)-minval(E)+1d0)*reshape([integer::],[Hmf%Hs],[reshape([integer::],[Nu],[0]),reshape([integer::],[Hmf%Hi-Nu],[1])]),ord)
			Nu=Nu*brizon%nk
			self%Ecfg(ord(1:self%ne(1)))=[1:self%ne(1)]
			self%Ecfg(ord(Nu+1:Nu+self%ne(2)))=[self%ne(1)+1:sum(self%ne)]
			if(self%sg/=3) then
				if(abs(E(ord(self%ne(1)))-E(ord(self%ne(1)+1)))<1d-7.or.abs(E(ord(Nu+self%ne(2)))-E(ord(Nu+self%ne(2)+1)))<1d-7) then
					if(this_image()==1) write(*,*)"warning, close-shell condition is not satisfied, maybe set to:",(self%ne(1)-count(abs((E(ord(self%ne(1)))-E(ord(:self%ne(1)))))<1d-7)),"or",(self%ne(1)+count(abs((-E(ord(self%ne(1)))+E(ord(self%ne(1)+1:Ns))))<1d-7))
				endif
				if(this_image()==1) then
					do i=Hmf%rg(1),Hmf%rg(2)
						if(Hmf%var(i)%label=="cp") then
							write(*,*)"cp is: ",Hmf%var(i)%val(1)-Hmf%var(i)%bd(1)*(E(ord(self%ne(1)))+E(ord(self%ne(1)+1)))/2d0
							exit
						endif
					enddo
				endif
			endif
		endif
		self%Emf=E

		if(self%sg==3) then
			if(is_mix) then
				nn_(0,:)=self%Ecfg(ord(sum(self%ne)-1:sum(self%ne)))
				self%Ecfg(ord(sum(self%ne)-1:sum(self%ne)))=0
			else
				nn_(0,:)=self%Ecfg(ord([self%ne(1),self%ne(2)+Nu]))
				self%Ecfg(ord([self%ne(1),self%ne(2)+Nu]))=0
			endif
			psi0_=0d0
			l=0
			do i=1,brizon%nk
				kq=0
				if(.not.is_in(-brizon%k(i,:)+self%dphy%var(1)%extdat,brizon%Tc,dr)) then
				endif
				do j=1,brizon%nk
					if(sum(abs(-brizon%k(i,:)+self%dphy%var(1)%extdat+dr-brizon%k(j,:)))<1d-5) then
						if(kq/=0) then
							if(this_image(self,2)==1) write(*,*)"warn! kq"
						else
							kq=j
						endif
					endif
				enddo
				
				if(kq==0) then
					if(this_image(self,2)==1) write(*,*)self%dphy%var(1)%extdat(:2)/pi,"warn!! kq2"
					error stop
				endif
				Hq(Hmf%i2H(:,1))=exp(-img*(dr(1)*latt%nb(0)%bd(1:latt%Ni)%r(1)+dr(2)*latt%nb(0)%bd(1:latt%Ni)%r(2)))
				do n1=(i-1)*Hmf%Hi+1,i*Hmf%Hi
					if((is_mix.or.Hmf%H2i(mod(n1-1,Hmf%Hi)+1,2)==2).and.self%Ecfg(n1)==0) then
						do n2=(kq-1)*Hmf%Hi+1,kq*Hmf%Hi
							if((is_mix.or.Hmf%H2i(mod(n2-1,Hmf%Hi)+1,2)==1).and.self%Ecfg(n2)==0.and.n2/=n1) then
								do j=l,1,-1
									if(nn_(j,1)==n2.and.nn_(j,2)==n1) then
										psi0_(j)=psi0_(j)+sum(conjg(H((i-1)*Hmf%Hi+Hmf%i2H(:,2),n1))*conjg(H((kq-1)*Hmf%Hi+Hmf%i2H(:,1),n2))*Hq(Hmf%i2H(:,1)))
										exit
									elseif(nn_(j,1)==n1.and.nn_(j,2)==n2) then
										psi0_(j)=psi0_(j)-sum(conjg(H((i-1)*Hmf%Hi+Hmf%i2H(:,2),n1))*conjg(H((kq-1)*Hmf%Hi+Hmf%i2H(:,1),n2))*Hq(Hmf%i2H(:,1)))
										exit
									endif
								enddo
								if(j==0) then
									l=l+1
									nn_(l,:)=[n2,n1]
									psi0_(l)=psi0_(l)+sum(conjg(H((i-1)*Hmf%Hi+Hmf%i2H(:,2),n1))*conjg(H((kq-1)*Hmf%Hi+Hmf%i2H(:,1),n2))*Hq(Hmf%i2H(:,1)))
								endif
							endif
						enddo
					endif
				enddo
			enddo

			if(allocated(self%nn)) deallocate(self%nn,self%Ok2,self%Ek2,self%psi0)
			allocate(self%nn(0:l,2),self%Ok2(l,l),self%Ek2(l,l),self%psi0(l))
			self%nn=nn_(:l,:)
			self%nn(0,:)=self%nn(1,[2,1])
			self%Ecfg(self%nn(0,:))=nn_(0,:)
			self%psi0=psi0_(:l)
			write(*,*)"************************",l
		endif
		H=matmul(Hmf%Uik,H)
		cH=0d0
		call Hmf%Hamilton([Hmf%rg(1):Hmf%rg(2)],cH)
		if(check) then
			tmp=sum(abs(matmul(transpose(conjg(H)),matmul(cH,H))-diag(E)))
			if(tmp>1d-6) then
				write(*,*)"hamilton err, plz check the supercell!",tmp
				error stop 
			endif
			check=.false.
		endif
		if(allocated(self%wf)) deallocate(self%wf)
		allocate(self%wf(size(H,1),size(H,2)))
		self%wf=conjg(H)


		if(self%sg==2) then
			if(allocated(self%dwf)) deallocate(self%dwf,self%g,self%S)
			allocate(self%dwf(Hmf%Hs,Hmf%Hs,sum(Hmf%var(1:)%n)),self%g(sum(Hmf%var(1:)%n)+sum(Hja%var(1:)%n)))
			allocate(self%S(size(self%g),size(self%g)))
			cH=transpose(conjg(H))
			self%dwf=0d0
			D=0d0
			call Hmf%dHamilton([1:Hmf%rg(2)],H,cH,D)

			do l=1,size(self%dwf,3)
				select case(self%opt)
				case(1)
					cH=0d0
					!$omp parallel do
					do i=1,size(E)
						do j=1,size(E)
							if(abs(E(i)-E(j))<1d-8) then
								cycle
							endif
							cH(:,i)=cH(:,i)+D(j,i,l)*H(:,j)/(E(i)-E(j))
						enddo
					enddo
					!$omp end parallel do
					self%dwf(:,:,l)=cH
				case(2)
					!$omp parallel do
					do i=1,size(E)
						do j=1,size(E)
							if(abs(E(i)-E(j))<1d-8.or.self%Ecfg(i)/=0.or.self%Ecfg(j)==0) then
								D(i,j,l)=0d0
							else
								D(i,j,l)=D(i,j,l)/(E(j)-E(i))
							endif
						enddo
					enddo
					!$omp end parallel do
					self%dwf(:,:,l)=matmul(matmul(H,D(:,:,l)),cH)
				end select
			enddo
		endif
	end subroutine
	subroutine do_vmc(self)
		class(t_mc) :: self[ica1,*]
		real(8) :: eg(size(self%g)),g_(size(self%g))
		integer :: k,l,l1,l2,seed
		self%samp=self%samp/ica2
		call random_number(seed,4285)
		call self%do_mc(seed+7*this_image())
		sync all
		if(this_image(self,2)==1) then
			select case(self%sg)
			case(1:2)
				do l=2,ica2
					do k=1,self%sphy%rg(2)
						self%sphy%var(k)%val=self%sphy%var(k)%val+self[this_image(self,1),l]%sphy%var(k)%val;
					enddo
					self%E=self%E+self[this_image(self,1),l]%E; self%err=self%err+self[this_image(self,1),l]%err
					if(self%sg==2) then
						self%S=self%S+self[this_image(self,1),l]%S; self%g=self%g+self[this_image(self,1),l]%g
					endif
				enddo
				do k=1,self%sphy%rg(2)
					self%sphy%var(k)%val=self%sphy%var(k)%val*1d0/ica2
				enddo
				self%E=self%E*1d0/ica2; self%err=sqrt(abs(self%err*1d0/ica2-self%E**2))
				if(self%sg==2) then
					self%S=self%S*1d0/ica2; self%g=self%g*1d0/ica2
				endif
			case(3)
				do l=2,ica2
					self%Ok2=self%Ok2+self[this_image(self,1),l]%Ok2; self%Ek2=self%Ek2+self[this_image(self,1),l]%Ek2
				enddo
				self%Ok2=self%Ok2*1d0/ica2; self%Ek2=self%Ek2*1d0/ica2
			end select
		endif
		sync all
		if(this_image(self,2)/=1) then
			select case(self%sg)
			case(1:2)
				do k=1,self%sphy%rg(2)
					self%sphy%var(k)%val=self[this_image(self,1),1]%sphy%var(k)%val
				enddo
				self%E=self[this_image(self,1),1]%E; self%err=self[this_image(self,1),1]%err
				if(self%sg==2) then
					self%S=self[this_image(self,1),1]%S; self%g=self[this_image(self,1),1]%g
				endif
			case(3)
				self%Ok2=self[this_image(self,1),1]%Ok2; self%Ek2=self[this_image(self,1),1]%Ek2
			end select
		endif
		sync all
		select case(self%sg)
		case(2)
			call heev(self%S,eg,'V')
			eg=eg+abs(min(eg(1),0d0))+0.2d0
			g_=self%g
			do l=1,size(self%g)
				self%g(l)=0d0
				do l1=1,size(self%g)
					do l2=1,size(self%g)
						if(abs(eg(l2)/eg(size(self%g)))<1d-3) then
							cycle
						endif
						self%g(l)=self%g(l)+real(self%S(l,l2)*conjg(self%S(l1,l2))*g_(l1)/eg(l2))
					enddo
				enddo
			enddo
			!grad=g
			self%g=self%g*Ns
		end select
		sync all
		self%samp=self%samp*ica2
	end subroutine
	subroutine do_mc(self,seed)
		class(t_mc) :: self[ica1,*]
		integer :: seed
		complex(8) :: iA(sum(self%ne),sum(self%ne)),iAl(size(iA,1),self%delay),A(size(iA,1),size(self%wf,2)),WA(size(self%wf,1),size(iA,2)),WAl(size(WA,1),size(iAl,2)),WAr(size(iAl,2),size(WA,2)),AW(size(iA,1),size(self%wf,2)),WAW(size(self%wf,1),size(self%wf,2)),AWr(size(iAl,2),size(AW,2)),wf(size(self%wf,1),size(self%wf,2)),dwf(size(self%dwf,1),size(self%dwf,2),size(self%dwf,3))
		integer :: cfgl(size(self%cfg)),Ecfgl(size(self%Ecfg)),icfg(size(A,1)),iEcfg(size(A,1))
		real(8) :: gp(size(self%g)),gl(size(gp))
		complex(8) :: Sp(size(gp),size(gp)),Op(size(gp)),Ol(size(gp)),Sl(size(gp),size(gp)),El,Ep
		complex(8) :: Ok(size(self%Ok2,1)),Ek(size(Ok)),Ok2p(size(Ok),size(Ok)),Ek2p(size(Ok),size(Ok))
		integer :: nn(lbound(self%nn,1):ubound(self%nn,1),size(self%nn,2))
		type(t_ham) :: sphyl,dphyl
		real(8) :: sphyp(sum(self%sphy%var(1:self%sphy%rg(2))%n)),dphyp(sum(self%dphy%var(1:self%dphy%rg(2))%n))
		real(8) :: rd,isamp,Oq
		complex(8) :: pb
		integer :: i,j,l,l1,l2,info,n,apt,samp,sg,dly,np
		integer :: k(0:4),m(0:4),dcfg(0:4)
		logical :: is_update,is_full,is_accept
		type(randomNumberSequence) :: rnd
		call mt_init_random_seed(rnd,seed)
		sphyl=self%sphy;dphyl=self%dphy
		sphyp=0d0;dphyp=0d0;Ep=0d0
		Sp=0d0; gp=0d0; Op=0d0
		Ek2p=0d0; Ok2p=0d0
		cfgl=self%cfg
		Ecfgl=self%Ecfg
		wf=self%wf
		if(allocated(self%dwf)) dwf=self%dwf
		if(allocated(self%nn)) nn=self%nn
		isamp=1d0/self%samp
		do i=1,size(cfgl)
			if(cfgl(i)/=0) icfg(cfgl(i))=i
			if(Ecfgl(i)/=0) iEcfg(Ecfgl(i))=i
		enddo
		A=wf(icfg,:)
		iA=A(:,iEcfg)
		call mat_inv(iA,info)
		n=0
		do
			n=n+1
			call change(cfgl,icfg,dcfg,k,rnd,.true.)
			cfgl(-dcfg(2:dcfg(0):2))=cfgl(-dcfg(2:dcfg(0):2))+cfgl(dcfg(1:dcfg(0):2))
			cfgl(dcfg(1:dcfg(0):2))=cfgl(-dcfg(2:dcfg(0):2))-cfgl(dcfg(1:dcfg(0):2))
			cfgl(-dcfg(2:dcfg(0):2))=cfgl(-dcfg(2:dcfg(0):2))-cfgl(dcfg(1:dcfg(0):2))
			icfg(k(1:k(0):2))=k(2:k(0):2)
			if(n>latt%Ns*10) then
				A=wf(icfg,:)
				iA=A(:,iEcfg)
				call mat_inv(iA,info)
				if(sum(abs(matmul(A(:,iEcfg),iA)-diag(1d0,size(A,1))))<1d-6) then
					exit
				else
					if(n>latt%Ns*10+200) then
						stop "ini cfg err"
					endif
				endif
			endif
		enddo
		WA=matmul(wf(:,iEcfg),iA)
		AW=matmul(iA,wf(icfg,:))
		WAW=matmul(wf(:,iEcfg),AW)
		WAW=WAW-wf
		n=0; apt=0; samp=0
		dly=0
		np=0
		is_update=.true.
		is_full=.false.
	lm:	do
			n=n+1
			call change(cfgl,icfg,dcfg,k,rnd,is_project)

			if(self%sg==3) then
			!if(.false.) then
				call get_Spb(cfgl,dcfg(1:dcfg(0)),Ecfgl,nn,pb,WA,iA,AW,WAW)
				if(isnan(real(pb))) then
					write(*,*)pb
					stop
				endif
			else
				call get_pb(k(1:k(0)),shape(0),pb,WA=WA,WAl=WAl(:,:dly),WAr=WAr(:dly,:))
				pb=pb*conjg(pb)
			endif
			!call random_number(rd)
			pb=pb*exp(2d0*jast(cfgl,dcfg(1:dcfg(0))))
			!do i=1,size(cfgl)/2
				!if(cfgl(i)>0.and.cfgl(i+Ns)==0) then
					!write(*,*)"double"
				!endif
			!enddo
			!stop
			rd=getrandomreal(rnd)
			is_accept=.false.
			if(rd<real(pb)) then

				if(self%sg==3) then
					i=-1
					call get_pb(k(1:k(0)),shape(0),pb,WA=WA,WAl=WAl(:,:dly),WAr=WAr(:dly,:))
					if(abs(pb)<1d-10) then
						do i=1,ubound(nn,1)
							if(get_row([nn(i,:),-nn(0,:)],Ecfgl,m,sg,[integer::])) then
								call get_pb(k(1:k(0)),m(1:m(0)),pb,WA=WA,WAl=WAl(:,:dly),WAr=WAr(:dly,:),iA=iA,iAl=iAl(:,:dly),AW=AW,WAW=WAW,AWr=AWr(:dly,:))
								if(abs(pb)>1d0) then
									Ecfgl(iEcfg(m(1:m(0):2)))=0
									Ecfgl(m(2:m(0):2))=m(1:m(0):2)
									iEcfg(m(1:m(0):2))=m(2:m(0):2)
									nn(0,:)=nn(i,[2,1])
									exit
								endif
							endif
						enddo
					endif

					if(i==ubound(nn,1)+1) then
						n=n-1
						write(*,*)"warn update, matrix is singularity, skip this configure"
						cycle lm
					endif
				endif
				apt=apt+1
				dly=dly+k(0)/2
				is_update=.true.
				is_accept=.true.
				np=np+1

				if(self%sg==2) rd=jast(cfgl,dcfg(1:dcfg(0)),Ol(sum(Hmf%var(1:)%n)+1:))

				cfgl(-dcfg(2:dcfg(0):2))=cfgl(-dcfg(2:dcfg(0):2))+cfgl(dcfg(1:dcfg(0):2))
				cfgl(dcfg(1:dcfg(0):2))=cfgl(-dcfg(2:dcfg(0):2))-cfgl(dcfg(1:dcfg(0):2))
				cfgl(-dcfg(2:dcfg(0):2))=cfgl(-dcfg(2:dcfg(0):2))-cfgl(dcfg(1:dcfg(0):2))
				icfg(k(1:k(0):2))=k(2:k(0):2)
				if(any(k(2:k(0):2)<=0)) then
					write(*,*)"k error"
					write(*,*)k
					error stop
				endif

				if(self%sg==3) then
					if(i>0.and.i/=ubound(nn,1)+1) then
						A=wf(icfg,:)
						iA=A(:,iEcfg)
						call mat_inv(iA,info)
						WA=matmul(wf(:,iEcfg),iA)
						AW=matmul(iA,wf(icfg,:))
						WAW=matmul(wf(:,iEcfg),AW)
						WAW=WAW-wf
						dly=0
						k(0)=0
					endif
				endif
			else
				k(0)=0
			endif
			if(((n>self%hot*self%step.and.mod(n-self%hot*self%step,self%step)==0).or.(dly>=self%delay-2))) is_full=.true.
			if((is_accept.or.is_full)) then
				if(self%sg==2.and.self%opt==1) then
					call update(k(1:k(0)),WA=WA,WAl=WAl(:,:dly),WAr=WAr(:dly,:),iA=iA,iAl=iAl(:,:dly),full=is_full)
				elseif(self%sg/=3) then
					call update(k(1:k(0)),WA=WA,WAl=WAl(:,:dly),WAr=WAr(:dly,:),full=is_full)
				else
					call update(k(1:k(0)),WA=WA,WAl=WAl(:,:dly),WAr=WAr(:dly,:),iA=iA,iAl=iAl(:,:dly),AW=AW,WAW=WAW,AWr=AWr(:dly,:),cwf=A,wf=wf,full=is_full)
				endif
				if(is_full) then
					is_full=.false.
					dly=0
					if(np>1024) then
						np=0
						pb=sum(abs(wf(:,iEcfg)-matmul(WA,wf(icfg,iEcfg))))
						if(abs(pb)>1d-5) then
							write(*,*)"warn!! ",apt,abs(pb)
							A=wf(icfg,:)
							iA=A(:,iEcfg)
							call mat_inv(iA,info)
							WA=matmul(wf(:,iEcfg),iA)
							if(self%sg==3) then
								AW=matmul(iA,wf(icfg,:))
								WAW=matmul(wf(:,iEcfg),AW)-wf
							endif
							pb=sum(abs(wf(:,iEcfg)-matmul(WA,wf(icfg,iEcfg))))
							if(abs(pb)>1d-5) then
								write(*,*)"update err",abs(pb)
								write(*,*)iEcfg
								write(*,*)icfg
							endif
						endif
					endif
				endif
			endif
			if(n>self%hot*self%step.and.mod(n-self%hot*self%step,self%step)==0) then
				samp=samp+1
				if(is_update) then
					is_update=.false.
					select case(self%sg)
					case(1)
						El=get_energy(cfgl,WA)*1d0/Ns
						do i=1,sphyl%rg(2)
							if(sphyl%var(i)%label(1:1)=="1") then
								call get_phy1(sphyl%var(i),cfgl,WA)
							elseif(sphyl%var(i)%label(1:1)=="2") then
								call get_phy2(sphyl%var(i),cfgl,WA)
							else
								error stop "label must be set to 1 or 2 in sphy"
							endif
						enddo
						if(dphyl%rg(2)==1) then
							call get_phy2(dphyl%var(1),cfgl,WA)
						endif
					case(2)
						El=get_energy(cfgl,WA)*1d0/Ns
						call get_O(icfg,iEcfg,WA,iA,dwf,Ol(:size(dwf,3)),self%opt)
					case(3)
						call get_overlap(cfgl,Ecfgl,nn,WA,iA,AW,WAW,Ok,Ek)
						Oq=1d0/sum(real(Ok*conjg(Ok)))
					end select
				endif
				sphyp=sphyp+sphyl%put([1:sphyl%rg(2)])
				dphyp=dphyp+dphyl%put([1:dphyl%rg(2)])
				Ep=Ep+El

				select case(self%sg)
				case(2)
					do l1=1,size(Ol)
						do l2=1,size(Ol)
							Sp(l2,l1)=Sp(l2,l1)+Ol(l2)*conjg(Ol(l1))
						enddo
					enddo
					gp=gp+real(El*conjg(Ol))
					Op=Op+Ol
				case(3)
					do l1=1,size(Ok)
						do l2=1,size(Ok)
							Ek2p(l2,l1)=Ek2p(l2,l1)+Ek(l2)*conjg(Ok(l2))*Oq
						enddo
					enddo
					do l1=1,size(Ok)
						do l2=1,size(Ok)
							Ok2p(l2,l1)=Ok2p(l2,l1)+Ok(l2)*conjg(Ok(l1))*Oq
						enddo
					enddo
				end select
				!!if(mod(samp,1024)==0.and.self%sg==3) then
				if((mod(samp,512)==0.or.samp==self%samp).and.self%sg==3) then
					!pb=sum(abs(Ek2p-conjg(transpose(Ek2p)))*1d0/samp)
					if(this_image(self,2)==1) then
						write(*,"(i11$)")this_image(),samp,self%samp
						write(*,"(*(es16.5))")real(apt)*1d0/n!,real(pb)
					endif
				endif
				!!read(*,*)
				if(samp==self%samp) then
					exit
				endif
			endif
		enddo lm
		! check
		pb=sum(abs(wf(:,iEcfg)-matmul(WA,wf(icfg,iEcfg))))
		!write(*,*)sum(abs(wf(:,iEcfg)-matmul(WA,wf(icfg,iEcfg)))),sum(abs(diag(1d0,size(iA,1))-matmul(iA,wf(icfg,iEcfg)))),sum(abs(matmul(wf(icfg,iEcfg),AW)-wf(icfg,:))),sum(abs(matmul(wf(:,iEcfg),AW)-WAW))
		if(abs(pb)>1d-5) then
			write(*,*)"warn!!!!!",abs(pb)
		endif

		sphyp=sphyp*isamp;dphyp=dphyp*isamp;Ep=Ep*isamp
		Sp=Sp*isamp; gp=gp*isamp; Op=Op*isamp
		Ek2p=Ek2p*isamp
		Ok2p=Ok2p*isamp
		do l1=1,size(Op)
			do l2=1,size(Op)
				!S(l1,l2)=2d0*SE(l1,l2)-O(l1)*g(l2)-O(l2)*g(l1)-2d0*E*S(l1,l2) ! maybe works
				!S(l1,l2)=2d0*(SE(l1,l2)-S(l1,l2)*E-O(l1)*g(l2)-O(l2)*g(l1))
				Sp(l1,l2)=Sp(l1,l2)-Op(l1)*conjg(Op(l2))
			enddo
		enddo
		gp=2d0*(gp-real(Ep*conjg(Op)))

		self%cfg=cfgl
		select case(self%sg)
		case(1:2)
			call self%sphy%get(sphyp,[1:sphyl%rg(2)])
			call self%dphy%get(dphyp,[1:dphyl%rg(2)])
			self%E=real(Ep)
			self%err=real(Ep)**2

			if(self%sg==2) then
				self%S=Sp; self%g=gp
			endif
		case(3)
			self%Ek2=Ek2p
			self%Ok2=Ok2p
		end select
		call finalize_RandomNumberSequence(rnd)
	end subroutine
	subroutine spect(self,ut,gm,omg,m,E)
		class(t_mc) :: self
		integer :: ut,m
		real(8), optional :: E
		real(8) :: gm,omg(:)
		complex(8) :: Sq(m),H_(size(self%Ek2,1),size(self%Ek2,2)),O_(size(self%Ek2,1),size(self%Ek2,2)),Om(size(self%Ek2,1),size(self%Ek2,2)),psi0_(size(self%psi0)),pnn(size(self%psi0))
		real(8) :: domg,EO(size(O_,1)),Eg(size(O_,1)),norm,Smax(2)
		integer :: l,i,j,n0
		domg=(omg(2)-omg(1))/m
		Sq=0d0
		Smax=0d0
		O_=self%Ok2
		H_=self%Ek2
		!H_=0.5d0*(H_+transpose(conjg(H_)))
		n0=1
		call heev(O_,EO,'V')
		do i=1,size(EO)
			if(EO(i)>1d-8) then
				n0=i
				write(*,*)n0,size(EO)
				exit
			endif
		enddo
		!check truncate
		H_=matmul(transpose(conjg(O_)),matmul(H_,O_))
		!write(*,*)maxval(abs(H_(:n0-1,:)))
		EO(1:n0-1)=0d0
		Om=diag(EO)
		call hegv(H_(n0:,n0:),Om(n0:,n0:),Eg(n0:),jobz="V")

		H_=matmul(O_,H_)
		psi0_=matmul(transpose(conjg(H_)),matmul(self%Ok2,self%psi0))

		!psi0_=EO*matmul(transpose(conjg(O_)),self%psi0)
		!psi0_=matmul(transpose(conjg(H_)),psi0_)

		psi0_=conjg(psi0_)*psi0_
		write(*,*)abs(sum(psi0_(1:n0-1)))

		do l=1,m
			Sq(l)=Sq(l)+sum(psi0_(n0:)/(omg(1)+domg*l-Eg(n0:)+self%E*Ns+img*gm))
		enddo
		norm=self%dphy%var(1)%val(1)/real(sum(psi0_(n0:)))
		write(*,"(es12.4$)")-sum(imag(Sq)*domg)/pi,real(sum(psi0_(n0:))),self%dphy%var(1)%val(1)
		write(ut,"('#q=',3es12.4,2i5)")self%dphy%var(1)%extdat/pi,i,size(EO)
		do l=1,m
			if(Smax(2)<abs(imag(Sq(l)))) Smax=[omg(1)+domg*l,abs(imag(Sq(l)))]
			write(ut,"(*(es17.9))")omg(1)+domg*l,Sq(l)*norm
		enddo
		write(ut,"(x)")
		write(*,*)"max intensity: ",Smax(1)

		if(present(E)) then
			write(*,*)sum(psi0_(n0:)/(E-Eg(n0:)+self%E*Ns+img*gm))
			!psi0_=matmul(self%Ok2,self%psi0)
			pnn=0d0
			psi0_=matmul(transpose(conjg(H_)),matmul(self%Ok2,self%psi0))
			do l=n0,size(Eg)
				
				!pnn=pnn+conjg(H_(:,l))*psi0_*sum(H_(:,l)*conjg(psi0_))/(E-Eg(l)+self%E*Ns+img*gm)
				!pnn=pnn+H_(:,l)*conjg(H_(:,l))*psi0_(l)/(E-Eg(l)+self%E*Ns+img*gm)
				!pnn=pnn+H_(:,l)*conjg(H_(:,l))*psi0_(l)/(E-Eg(l)+self%E*Ns+img*gm)
				pnn=pnn+conjg(matmul(self%Ok2,H_(:,l))*psi0_(l))*matmul(self%Ok2,H_(:,l))*psi0_(l)/(E-Eg(l)+self%E*Ns+img*gm)
			enddo
			write(*,*)sum(pnn)
			do l=1,size(pnn)
				write(ut+1,"(es17.9$)")brizon%k(mod(abs(self%nn(l,1))-1,brizon%nk)+1,:2),brizon%k(mod(abs(self%nn(l,2))-1,brizon%nk)+1,:2),self%Emf(abs(self%nn(l,:))),pnn(l)
				write(ut+1,"(2i4)")abs(self%nn(l,:))
			enddo
		endif
		write(ut+1,"(x)")
	end subroutine
	subroutine do_var(self,n)
		class(t_mc) :: self[ica1,*]
		integer :: omp,n
		real(8) :: x(sum(Hmf%var(1:)%n)+sum(Hja%var(1:)%n),n),er,pgrad(size(x)),El(n)
		integer :: i,hot,ord(n)
		logical :: init_cfg
		self%sg=2
		hot=self%hot
		i=0
		if(size(x)==0) then
			call self%init(.true.)
			call self%do_vmc()
			if(this_image()==1) then
				write(*,"(es12.4$)")Hmf%var(1:)%val(1),Hja%var(1:)%val(1),self%E,self%err
				write(*,"(i3$)")int(sign(1d0,self%g))
			endif
			return
		endif
		x(:sum(Hmf%var(1:)%n),1)=Hmf%put([1:Hmf%rg(2)])
		x(sum(Hmf%var(1:)%n)+1:,1)=Hja%put([1:Hja%rg(2)])
		init_cfg=.true.
		do 
			i=i+1
			if(allocated(self%g)) pgrad=self%g
			call Hmf%get(x(:sum(Hmf%var(1:)%n),i),[1:Hmf%rg(2)])
			call Hja%get(x(sum(Hmf%var(1:)%n)+1:,i),[1:Hja%rg(2)])
			if(this_image()==1) write(*,"(i4$)")i
			call self%init(init_cfg)
			call self%do_vmc()
			if(this_image()==1) then
				write(*,"(es12.4$)")Hmf%var(1:)%val(1),Hja%var(1:)%val(1),self%E,self%err
				write(*,"(i3$)")int(sign(1d0,self%g))
			endif
			El(i)=self%E
			ord(i)=i
			er=er*1.5d0
			!if((El(i)-er)>El(max(i-1,1))) then
				!grad=pgrad
				!x=x+grad*dx
				!dx=dx*0.8d0
				!write(*,"(' err',es10.2$)")dx,(El(i)-er),El(max(i-1,1))
				!i=i-1
			!endif
			!if(all(abs(grad*dx/x)<1d-5).or.dx<1d-3.or.i==size(El)) then
				!if(minval(El(:i))<(El(i)-er)) then
					!write(*,*)"var err, min E is ", minval(El(:i))
				!endif
			if(i==size(El)) then
				exit
			endif
			x(:,i+1)=x(:,i)-self%g*dx
			!self%hot=128
			!init_cfg=.false.
			if(this_image()==1) write(*,"(x)")
			if(this_image()==1) write(20,"(x)")
		enddo
		self%hot=hot
		if(this_image()==1) write(*,*)"finished"
		call qsort(El,ord)
		do i=2,min(n/5,30)
			x(:,ord(1))=x(:,ord(1))+x(:,ord(i))
		enddo
		call Hmf%get(x(:sum(Hmf%var(1:)%n),ord(1))/min(n/5,30),[1:Hmf%rg(2)])
		call Hja%get(x(sum(Hmf%var(1:)%n)+1:,ord(1))/min(n/5,30),[1:Hja%rg(2)])
		!E=minval(El(:i))
	end subroutine
end module
