include "../lib/hamilton_final_1.f90"
module global
	use M_hamilton_final_m
	use ifport, ifort_qsort => qsort
	use M_omp_rand
	use omp_lib
	use mkl_service
	implicit none
	real(8) :: t(2)=[1d0,-0.3d0]
	real(8), parameter :: DJ=0.3d0,V=0d0,U=0d0
	integer, parameter :: opt=2 ! 1: Tr(AdA)
								! 2: O=cicj
	integer, parameter :: iE=1,ier=2,idsc=3,iaf=4,iddw=5,iSq_pm=6,iSq_zz=7,iCq=8
	integer, parameter :: ica1=1
	logical :: is_project=.true.
	integer :: ica2
	integer :: Ns,icp=-111
	type(t_ham) :: Hmf,Ham,Hja
	real(8) :: q(3,ica1)=reshape([&
		pi,0d0,0d0&
		],[3,ica1])
contains
	subroutine initial()
		integer :: i,l,idx
		real(8) :: q(3)
		q=[1d0/8d0,0d0,0d0]*2d0*pi
		!call init_random_seed()
		! lattice
		latt%is_all=.true.
		latt%a1=[1d0,0d0,0d0]
		latt%a2=[0d0,1d0,0d0]
		latt%c1=latt%a1
		latt%c2=latt%a2
		!latt%c1=[1d0,1d0,0d0]
		!latt%c2=[-1d0,1d0,0d0]
		!latt%c1=[4d0,1d0,0d0]
		!latt%c2=[0d0,2d0,0d0]
		latt%T1=[1d0,0d0,0d0]*8d0
		latt%T2=[0d0,1d0,0d0]*8d0
		latt%bdc=[-1d0,1d0,0d0]
		allocate(latt%rsb(1,3))
		latt%rsb(1,:)=[0d0,0d0,0d0]
		call latt%gen_latt(size(t))
		call latt%gen_brizon(brizon)
		if(this_image()==1) then
			call check_lattice(101)
			write(*,*)"Total site number is: ",latt%Ns
		endif

		allocate(Ham%var(-10:10),Hmf%var(-10:10),Hja%var(-10:10))
		!! cp
		!idx=Hmf%add(nb=0,ca=[c("i",1,+1),c("i",1,-1),c("i",1,+2),c("i",1,-2)],n=2,sg=[1d0,-1d0],label="cp")
		!Hmf%var(idx)%bd=-1d0
		!Hmf%var(idx)%val=0d0
		!icp=idx

		! d-wave sc
		idx=Hmf%add(nb=1,ca=[c("i",1,+1),c("j",1,-2),c("j",1,+1),c("i",1,-2),c("j",1,+2),c("i",1,-1),c("i",1,+2),c("j",1,-1)],n=2)
		do i=1,size(Hmf%var(idx)%bd)
			if(abs(latt%nb(Hmf%var(idx)%nb)%bd(i)%dr(1))>1d-6) then
				Hmf%var(idx)%bd(i)=1d0
			else
				Hmf%var(idx)%bd(i)=-1d0
			endif
		enddo
		!Hmf%var(idx)%val=2.92d-1
		Hmf%var(idx)%val=3.24d-1

		!! ddw
		!idx=Hmf%add(nb=1,ca=[c("i",1,+1),c("j",1,-1),c("j",1,+1),c("i",1,-1),c("j",1,+2),c("i",1,-2),c("i",1,+2),c("j",1,-2)],n=2,sg=[1d0,-1d0,-1d0,1d0])
		!do i=1,size(Hmf%var(idx)%bd)
			!if(abs(latt%nb(Hmf%var(idx)%nb)%bd(i)%dr(mod(nint(sum(latt%nb(0)%bd(latt%nb(Hmf%var(idx)%nb)%bd(i)%i(1))%r)),2)+1))>1d-6) then
				!Hmf%var(idx)%bd(i)=img
			!else
				!Hmf%var(idx)%bd(i)=-img
			!endif
		!enddo
		!Hmf%var(idx)%val=3.24d-1

		!! sdw
		!idx=Hmf%add(nb=0,ca=[c("i",1,+1),c("i",1,-1),c("i",1,+2),c("i",1,-2)],n=2)
		!do i=1,size(Hmf%var(idx)%bd)
			!if(mod(nint(sum(latt%nb(0)%bd(latt%nb(Hmf%var(idx)%nb)%bd(i)%i(1))%r)),2)==0) then
				!Hmf%var(idx)%bd(i)=1d0
			!else
				!Hmf%var(idx)%bd(i)=-1d0
			!endif
		!enddo
		!Hmf%var(idx)%val=0d0

		! bond order
		do l=2,size(t)
			idx=Hmf%add(nb=l,ca=[c("i",1,+1),c("j",1,-1),c("j",1,+1),c("i",1,-1),c("i",1,+2),c("j",1,-2),c("j",1,+2),c("i",1,-2)],n=2,sg=[+1d0,+1d0,-1d0,-1d0])
			Hmf%var(idx)%bd=-1d0
			Hmf%var(idx)%val=-t(l)
		enddo

		! hp
		do l=1,size(t)
			idx=Hmf%add(nb=l,ca=[c("i",1,+1),c("j",1,-1),c("j",1,+1),c("i",1,-1),c("i",1,+2),c("j",1,-2),c("j",1,+2),c("i",1,-2)],n=2,sg=[+1d0,+1d0,-1d0,-1d0],is_var=.false.)
			Hmf%var(idx)%bd=-1d0
			Hmf%var(idx)%val=t(l)
		enddo

		call Hmf%init()

		!idx=Hja%add(nb=0,ca=[c("i",1,+1),c("i",1,-1),c("i",1,-2),c("i",1,+2)],n=4)
		!Hja%var(idx)%bd=-1d0
		!Hja%var(idx)%val=0d0

		call Hja%init()

		do l=1,size(t)
			idx=Ham%add(nb=l,ca=[c("i",1,+1),c("j",1,-1),c("j",1,+1),c("i",1,-1),c("i",1,-2),c("j",1,+2),c("j",1,-2),c("i",1,+2)],n=2,label="tp")
			Ham%var(idx)%bd=-1d0
			Ham%var(idx)%val=t(l)
		enddo

		idx=Ham%add(nb=1,ca=[c("i",1,+1),c("i",1,+2),c("j",1,-2),c("j",1,-1),c("j",1,+1),c("j",1,+2),c("i",1,-2),c("i",1,-1)],n=4,V=DJ/2d0,label="J")
		Ham%var(idx)%bd=1d0
		Ham%var(idx)%val=1d0

		idx=Ham%add(nb=1,ca=[c("i",1,+1),c("i",1,-1),c("j",1,+1),c("j",1,-1),c("i",1,-2),c("i",1,+2),c("j",1,-2),c("j",1,+2),c("i",1,+1),c("i",1,-1),c("j",1,-2),c("j",1,+2),c("i",1,-2),c("i",1,+2),c("j",1,+1),c("j",1,-1)],n=4,sg=[V+DJ/4d0,V+DJ/4d0,V-DJ/4d0,V-DJ/4d0],label="4term")
		Ham%var(idx)%bd=1d0
		Ham%var(idx)%val=1d0

		call Ham%init()

		Ns=latt%Ns
	end subroutine
end module

module mc_utility
	use global
	use blas95, only : gemm, gemv, her, gerc
	implicit none
contains
	subroutine get_pb(k,m,pb,WA,WAl,WAr,iA,iAl,AW,WAW,AWr)
		integer :: k(:),m(:)
		complex(8) :: WA(:,:),pb
		complex(8), optional :: WAl(:,:),WAr(:,:),iA(:,:),iAl(:,:),AW(:,:),WAW(:,:),AWr(:,:)
		integer :: i,j,n
		integer :: k_(size(k)/2+size(m)/2),m_(size(k_))
		complex(8) :: A(size(k_,1),size(k_,1))
		n=size(k)/2
		k_=[k(1::2),m(1::2)]
		m_=[k(2::2),m(2::2)]
		do i=1,size(m_)
			do j=1,size(k_)
				if(i<=n.and.j<=n) then
					if(present(WAl)) then
						A(i,j)=WA(m_(i),k_(j))+sum(WAl(m_(i),:)*WAr(:,k_(j)))
					else
						A(i,j)=WA(m_(i),k_(j))
					endif
				elseif(i<=n.and.j>n) then
					if(present(AWr)) then
						A(i,j)=WAW(m_(i),m_(j))+sum(WAl(m_(i),:)*AWr(:,m_(j)))
					else
						A(i,j)=WAW(m_(i),m_(j))
					endif
				elseif(i>n.and.j<=n) then
					if(present(iAl)) then
						A(i,j)=iA(k_(i),k_(j))+sum(iAl(k_(i),:)*WAr(:,k_(j)))
					else
						A(i,j)=iA(k_(i),k_(j))
					endif
				else
					if(present(AWr)) then
						A(i,j)=AW(k_(i),m_(j))+sum(iAl(k_(i),:)*AWr(:,m_(j)))
					else
						A(i,j)=AW(k_(i),m_(j))
					endif
				endif
			enddo
		enddo
		!A(1:n,1:n)=WA(m_(1:n),k_(1:n))
		!A(1:n,n+1:)=WAW(m_(1:n),m_(n+1:))
		!A(n+1:,1:n)=iA(k_(n+1:),k_(1:n))
		!A(n+1:,n+1:)=AW(k_(n+1:),m_(n+1:))
		select case(size(k_))
		case(0)
			pb=1d0
		case(1)
			pb=A(1,1)
		case(2)
			pb=(A(1,1)*A(2,2)-A(1,2)*A(2,1))
		case(3)
			pb=A(1,1)*A(2,2)*A(3,3)-A(1,1)*A(2,3)*A(3,2)-A(1,2)*A(2,1)*A(3,3)+A(1,2)*A(2,3)*A(3,1)+A(1,3)*A(2,1)*A(3,2)-A(1,3)*A(2,2)*A(3,1)
		case(4)
			pb=A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))-&
			   A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))+&
			   A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))-&
			   A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
	   case(5:)
		   write(*,*)"above 5 is not considered"
		   stop
	   end select
	end subroutine
	subroutine update(k,WA,WAl,WAr,iA,iAl,AW,WAW,AWr,cwf,wf,full)
		integer :: k(:)
		complex(8) :: WA(:,:),WAl(:,:),WAr(:,:)
		complex(8), optional :: AW(:,:),WAW(:,:),iA(:,:),cwf(:,:),wf(:,:),iAl(:,:),AWr(:,:)
		logical :: full
		integer :: k_(size(k)/2),m_(size(k_))
		integer :: n,i,j,l,ct=0
		complex(8) :: A(size(k_),size(k_)),beta=1d0
		if(size(k)>0) then
			n=size(WAl,2)-size(k_)
			k_=[k(1::2)]
			m_=[k(2::2)]
			do i=1,size(m_)
				do j=1,size(k_)
					A(i,j)=WA(m_(i),k_(j))+sum(WAl(m_(i),:n)*WAr(:n,k_(j)))
				enddo
			enddo
			select case(size(k_))
			case(0)
				return
			case(1)
				A=-1d0/A
			case(2)
				A=-1d0/(A(1,1)*A(2,2)-A(2,1)*A(1,2))*reshape([A(2,2),-A(2,1),-A(1,2),A(1,1)],[2,2])
			case default
				write(*,*)"err"
				stop
			end select
			WAl(:,n+1:)=WA(:,k_)
			WAr(n+1:,:)=WA(m_,:)
			do i=1,size(k_)
				call gemv(WAl(:,:n),WAr(:n,k_(i)),WAl(:,n+i),beta=beta)
				call gemv(WAr(:n,:),WAl(m_(i),:n),WAr(n+i,:),beta=beta,trans='T')
				WAr(n+i,k_(i))=WAr(n+i,k_(i))-1d0
			enddo
			WAr(n+1:,:)=matmul(A,WAr(n+1:,:))
			if(present(iA)) then
				do i=1,size(k_)
					iAl(:,n+i)=iA(:,k_(i))
					call gemv(iAl(:,:n),WAr(:n,k_(i)),iAl(:,n+i),beta=beta)
				enddo
				if(present(WAW)) then
					cwf(k_,:)=wf(m_,:)-cwf(k_,:)
					AWr(n+1:,:)=WAW(m_,:)+matmul(WAl(m_,n+1:),cwf(k_,:))
					do i=1,size(k_)
						call gemv(AWr(:n,:),WAl(m_(i),:n),AWr(n+i,:),beta=beta,trans='T')
					enddo
					AWr(n+1:,:)=cwf(k_,:)+matmul(A,AWr(n+1:,:))
					cwf(k_,:)=wf(m_,:)
				endif
			endif
		endif
		if(full.and.size(WAl,2)>0) then
			if(size(WAl,2)>2) then
			!if(.true.) then
				call gemm(WAl,WAr,WA,beta=beta)
				if(present(iA)) then
					call gemm(iAl,WAr,iA,beta=beta)
					if(present(WAW)) then
						call gemm(iAl,AWr,AW,beta=beta)
						call gemm(WAl,AWr,WAW,beta=beta)
					endif
				endif
			else
				!$omp parallel
				!$omp do
				do i=1,size(WA,2)
					do j=1,size(WAl,2)
						do l=1,size(WA,1)
							WA(l,i)=WA(l,i)+WAl(l,j)*WAr(j,i)
						enddo
					enddo
				enddo
				!$omp end do
				if(present(iA)) then
					!$omp do
					do i=1,size(iA,2)
						do j=1,size(WAl,2)
							do l=1,size(iA,1)
								iA(l,i)=iA(l,i)+iAl(l,j)*WAr(j,i)
							enddo
						enddo
					enddo
					!$omp end do
					if(present(WAW)) then
						!$omp do
						do i=1,size(AW,2)
							do j=1,size(WAl,2)
								do l=1,size(AW,1)
									AW(l,i)=AW(l,i)+iAl(l,j)*AWr(j,i)
								enddo
								do l=1,size(WAW,1)
									WAW(l,i)=WAW(l,i)+WAl(l,j)*AWr(j,i)
								enddo
							enddo
						enddo
						!$omp end do
					endif
				endif
				!$omp end parallel
			endif
			!ct=ct+1
			!if(mod(ct,1000)==0) write(*,*)ct
		endif
	end subroutine
end module

module model
	use mc_utility
	implicit none
contains
	subroutine change(cfg,icfg,dcfg,k,rnd,is_P)
		integer :: cfg(:),icfg(:)
		integer :: dcfg(0:),k(0:)
		type(randomNumberSequence) :: rnd
		integer :: i1,j1,i2,j2,n
		logical :: is_P
		!call random_number(i,size(icfg))
		i1=mt_random_number(rnd,size(icfg))
		i1=icfg(i1)
		n=Hmf%H2s(i1,2)
		i2=Hmf%s2H(Hmf%H2s(i1,1),mod(n,2)+1)
		do 
			!call random_number(j,Ns)
			j1=mt_random_number(rnd,Ns)
			j2=Hmf%s2H(j1,mod(n,2)+1)
			j1=Hmf%s2H(j1,n)
			if(cfg(j1)==0) then
				if((.not.is_P).or.sign(1,-cfg(j2))==sign(1,-cfg(i2))) then
					dcfg(0:2)=[2,j1,-i1]
					k(0:2)=[2,cfg(abs(dcfg(2))),dcfg(1)]
				else
					dcfg(0:4)=[4,j2,-i2,j1,-i1]
					k(0:4)=[4,cfg(abs(dcfg(4))),dcfg(3),cfg(abs(dcfg(2))),dcfg(1)]
				endif
				exit
			endif
		enddo
	end subroutine
	logical function get_row(ca,cfg,k,sg,P)
		integer, intent(in) :: cfg(:),ca(:)
		integer, intent(out) :: k(0:),sg
		integer, intent(in) :: P(:)
		integer :: i,j,n1,n2,l,tp1,tp2,n,tmp
		logical :: flag(size(ca)),fP(2)
		get_row=.false.
		flag=.true.
		n1=-1; n2=0
		l=0
		sg=1
		n=size(P)
		do i=size(ca),1,-1
			tp1=merge(1,2,Hmf%H2s(abs(ca(i)),2)==1)
			fP(tp1)=(cfg(abs(ca(i)))/=0)
			if(flag(i)) then
				do j=i-1,1,-1
					if(flag(j)) then
						if(ca(i)==-ca(j)) then
							if((ca(i)>0).xor.fP(tp1)) then
								flag(j)=.false.
								sg=sg*sign(1,-mod(count(flag(j+1:i-1)),2))
								exit
							else
								!write(*,*)"already return"
								!read(*,*)
								return
							endif
						elseif(ca(i)==ca(j)) then
							!write(*,*)"double return"
							!read(*,*)
							return 
						endif
					endif
				enddo
				if(j==0) then
					if(ca(i)<0.and.fp(tp1)) then
						n1=n1+2
						k(n1)=cfg(-ca(i))
						if(n2>n1.and.mod((n2-n1)/2,2)==0) sg=-sg
					elseif(ca(i)>0.and.(.not.fp(tp1))) then
						n2=n2+2
						k(n2)=ca(i)
						if(n1>n2.and.mod((n1-n2)/2,2)==0) sg=-sg
					else
						!write(*,*)"already return"
						!read(*,*)
						return
					endif
				endif
			endif
			!perform project, should be edit if changes
			if(n>0) then
				fP(tp1)=.not.(fP(tp1).xor.flag(i))
				tp2=mod(tp1,2)+1
				tmp=Hmf%s2H(Hmf%H2s(abs(ca(i)),1),tp2)
				fP(tp2)=cfg(tmp)/=0
				do j=i+1,size(ca)
					if(abs(ca(j))==tmp) then
						fp(tp2)=.not.fp(tp2)
					endif
				enddo
				if(fP(1).and.(.not.fP(2))) then
					l=l-1
				else
					fP(tp1)=(ca(i)>0)
					if(fP(1).and.(.not.fP(2))) then
						l=l+1
					endif
				endif
				if(P(n)==i) then
					n=n-1
					if(l/=0) then
						!write(*,*)"Project return"
						!read(*,*)
						return
					endif
				endif
			endif
		enddo
		!if(n2==n1+1) then
			get_row=.true.
			k(0)=n2
			!write(*,*)"ca: ",ca
			!write(*,*)"cg: ",cfg(abs(ca))
			!write(*,*)"cg':",(cfg(Hmf%s2H(Hmf%H2s(abs(ca(i)),1),mod(Hmf%H2s(abs(ca(i)),2),2)+1)),i=1,size(ca))
			!write(*,*)"k: ",k(1:k(0))
			!write(*,*)"sg: ",sg
			!read(*,*)
		!else
			!write(*,*)"error",n1,n2,k(1:k(0))
			!stop
		!endif
	end function
	complex(8) function get_spin_pm(cfg,D,q)
		complex(8) :: D(:,:)
		integer :: cfg(:)
		real(8) :: q(3)
		complex(8) :: pb
		integer :: i,j,sg,l
		integer :: k(0:4)
		get_spin_pm=0d0
		!$omp parallel do collapse(2) reduction(+:get_spin_pm) private(pb,k,sg)
		do i=1,Ns
			do j=1,Ns
				if(get_row([i,i+Ns,-j-Ns,-j],cfg,k,sg,pack([1],[is_project]))) then
					call get_pb(k(1:k(0)),shape(0),pb,D)
					get_spin_pm=get_spin_pm+pb*exp(img*sum(q*(latt%nb(0)%bd(i)%r-latt%nb(0)%bd(j)%r)))*sg*exp(jast(cfg,[i,i+Ns,-j-Ns,-j]))
				endif
			enddo
		enddo
		!$omp end parallel do
		get_spin_pm=get_spin_pm/Ns
	end function
	complex(8) function get_spin_zz(cfg,q)
		integer :: cfg(:)
		real(8) :: q(3)
		integer :: i,j
		get_spin_zz=0d0
		!$omp parallel do collapse(2) reduction(+:get_spin_zz)
		do i=1,Ns
			do j=1,Ns
				get_spin_zz=get_spin_zz+0.25d0*sum(((1-sign(1,-cfg(i::Ns)))/2-[0,1]))*sum(((1-sign(1,-cfg(j::Ns)))/2-[0,1]))*exp(img*sum(q*(latt%nb(0)%bd(i)%r-latt%nb(0)%bd(j)%r)))
			enddo
		enddo
		!$omp end parallel do
		get_spin_zz=get_spin_zz/Ns
	end function
	complex(8) function get_af(cfg,q)
		integer :: cfg(:)
		real(8) :: q(3)
		integer :: i
		get_af=0d0
		!$omp parallel do reduction(+:get_af)
		do i=1,Ns
			get_af=get_af+sum(((1-sign(1,-cfg(i::Ns)))/2-[0,1]))*exp(img*sum(q*latt%nb(0)%bd(i)%r))
		enddo
		!$omp end parallel do
		get_af=0.5d0*get_af/Ns
	end function
	subroutine get_phy_n(cfg,phy_n)
		integer :: cfg(:)
		real(8) :: phy_n(:,:)
		integer :: i
		do i=1,Ns
			phy_n(:,i)=((1-sign(1,-cfg(i::Ns)))/2-[0,1])*[1,-1]
		enddo
	end subroutine
	subroutine get_phy_sc(cfg,WA,phy_sc)
		complex(8) :: WA(:,:)
		integer :: cfg(:)
		real(8) :: phy_sc(:)
		complex(8) :: pb,bdc1,bdc2
		integer :: i1,j1,i2,j2,sg,l,n1,n2,p1,p2
		integer :: k(0:6)
		logical :: flag
		phy_sc=0d0
		do n1=1,size(latt%nb(1)%bd)
			if(abs(latt%nb(1)%bd(n1)%r(1))<1d-4) then
				do n2=1,size(latt%nb(1)%bd)
					if(any(latt%nb(1)%bd(n1)%i==latt%nb(1)%bd(n2)%i(1)).or.any(latt%nb(1)%bd(n1)%i==latt%nb(1)%bd(n2)%i(2))) cycle
					i1=latt%nb(1)%bd(n1)%i(1)
					j1=latt%nb(1)%bd(n1)%i(2)
					i2=latt%nb(1)%bd(n2)%i(1)
					j2=latt%nb(1)%bd(n2)%i(2)
					bdc1=conjg(latt%nb(1)%bd(n1)%bdc)
					bdc2=latt%nb(1)%bd(n2)%bdc
					do p1=1,2
						do p2=1,2
							if(get_row([-j1-Ns,i1,-i2,j2+Ns],cfg,k,sg,[5])) then
								call get_pb(k(1:k(0)),shape(0),pb,WA)
								phy_sc(n2)=phy_sc(n2)+real(pb*sg*bdc1*bdc2*exp(jast(cfg,[-j1-Ns,i1,-i2,j2+Ns])))
							endif
							call swap(i2,j2)
							bdc2=conjg(bdc2)
						enddo
						call swap(i1,j1)
						bdc1=conjg(bdc1)
					enddo
				enddo
			endif
		enddo
	end subroutine
	complex(8) function get_dsc(cfg,WA,dcfg)
		complex(8) :: WA(:,:)
		integer :: cfg(:)
		integer, optional :: dcfg(:)
		complex(8) :: pb,bdc1,bdc2
		integer :: i1,j1,i2,j2,sg,l,n1,n2,p1,p2
		integer :: k(0:6),dcfg_(0:2)
		if(present(dcfg)) then
			dcfg_(0:)=[size(dcfg),dcfg]
		else
			dcfg_(0)=0
		endif
		get_dsc=0d0
		!$omp parallel do collapse(2) reduction(+:get_dsc) private(pb,k,sg,i1,j1,i2,j2,bdc1,bdc2)
		do n1=1,size(latt%nb(1)%bd)
			do n2=1,size(latt%nb(1)%bd)
				if(n1==n2) cycle
				i1=latt%nb(1)%bd(n1)%i(1)
				j1=latt%nb(1)%bd(n1)%i(2)
				i2=latt%nb(1)%bd(n2)%i(1)
				j2=latt%nb(1)%bd(n2)%i(2)
				bdc1=conjg(latt%nb(1)%bd(n1)%bdc)
				bdc2=latt%nb(1)%bd(n2)%bdc
				do p1=1,2
					do p2=1,2
						if(get_row([-j1-Ns,i1,-i2,j2+Ns,dcfg_(1:dcfg_(0))],cfg,k,sg,[5])) then
							call get_pb(k(1:k(0)),shape(0),pb,WA)
							if(sum(abs(latt%nb(1)%bd(n1)%dr*latt%nb(1)%bd(n2)%dr))>1d-6) then
								get_dsc=get_dsc+pb*sg*bdc1*bdc2*exp(jast(cfg,[-j1-Ns,i1,-i2,j2+Ns,dcfg_(1:dcfg_(0))]))
							else
								get_dsc=get_dsc-pb*sg*bdc1*bdc2*exp(jast(cfg,[-j1-Ns,i1,-i2,j2+Ns,dcfg_(1:dcfg_(0))]))
							endif
						endif
						call swap(i2,j2)
						bdc2=conjg(bdc2)
					enddo
					call swap(i1,j1)
					bdc1=conjg(bdc1)
				enddo
			enddo
		enddo
		!$omp end parallel do
		get_dsc=get_dsc/(Ns**2)
	end function
	complex(8) function get_energy(cfg,WA,dcfg,m,iA,AW,WAW)
		complex(8) :: WA(:,:)
		integer :: cfg(:)
		integer, optional :: dcfg(:),m(:)
		complex(8), optional :: AW(:,:),WAW(:,:),iA(:,:)
		complex(8) :: pb,bdc,val
		integer :: i,j,n,l,p,sg,ic,jc,sb(2)
		integer :: k(0:4),dcfg_(0:2),dcfg_all(0:4)
		integer :: Pr(0:2)
		logical :: flag
		if(present(dcfg)) then
			dcfg_(0:)=[size(dcfg),dcfg]
			Pr=[2,1,dcfg_(0)+1]
		else
			dcfg_(0)=0
			Pr(0:1)=[1,1]
		endif
		if(.not.is_project) then
			Pr(0)=0
		endif
		get_energy=0d0
		do l=minval(Ham%var%nb),maxval(Ham%var%nb)
			do n=1,size(latt%nb(l)%bd)
				i=latt%nb(l)%bd(n)%i(1)
				j=latt%nb(l)%bd(n)%i(2)
				sb=latt%nb(l)%bd(n)%sb
				do p=Ham%rg(1),Ham%rg(2)
					if(Ham%var(p)%nb==l) then
						if(.not.isnan(Ham%var(p)%V)) then
							val=Ham%var(p)%val(Ham%var(p)%bd2v(n))*Ham%var(p)%bd(n)*Ham%var(p)%V
						else
							val=Ham%var(p)%val(Ham%var(p)%bd2v(n))*Ham%var(p)%bd(n)
						endif
						if(size(Ham%var(p)%c,1)==0) then
							if(present(m)) then
								call get_pb(shape(0),m,pb,WA=WA,iA=iA,AW=AW,WAW=WAW)
								get_energy=get_energy+val*pb
							else
								get_energy=get_energy+val
							endif
							cycle
						endif
						do jc=1,size(Ham%var(p)%c,2)
							dcfg_all(0)=dcfg_(0)+size(Ham%var(p)%c,1)
							dcfg_all(1:dcfg_all(0))=[dcfg_(1:dcfg_(0)),(Hmf%s2H(merge(i,j,Ham%var(p)%c(ic,jc)%l=="i"),abs(Ham%var(p)%c(ic,jc)%tp))*sign(1,Ham%var(p)%c(ic,jc)%tp),ic=1,size(Ham%var(p)%c,1))]
							if(get_row(dcfg_all(1:dcfg_all(0)),cfg,k,sg,Pr(1:Pr(0)))) then
								if(present(m)) then
									call get_pb(k(1:k(0)),m,pb,WA=WA,iA=iA,AW=AW,WAW=WAW)
								else
									call get_pb(k(1:k(0)),shape(0),pb,WA=WA)
								endif
								if(size(Ham%var(p)%c,1)==2.and.Ham%var(p)%c(1,jc)%l/=Ham%var(p)%c(2,jc)%l) then
									if(Ham%var(p)%c(1,jc)%l=="i") then
										bdc=latt%nb(l)%bd(n)%bdc
									else
										bdc=conjg(latt%nb(l)%bd(n)%bdc)
									endif
								else
									bdc=abs(latt%nb(l)%bd(n)%bdc)
								endif
								get_energy=get_energy+sg*val*Ham%var(p)%sg(jc)*pb*bdc*exp(jast(cfg,dcfg_all(1:dcfg_all(0))))
							endif
						enddo
					endif
				enddo
			enddo
		enddo
	end function
	subroutine get_O(icfg,iEcfg,WA,iA,dwf,O)
		complex(8) :: WA(:,:),iA(:,:),dwf(:,:,:),O(:)
		complex(8) :: beta,tmp
		integer :: icfg(:),iEcfg(:)
		integer :: l,i,j,lp
		O=0d0
		beta=1d0
		select case(opt)
		case(1)
			!$omp parallel do reduction(+:O)
			do i=1,size(iA,1)
				call gemv(dwf(icfg,iEcfg(i),:),iA(i,:),O,trans="T",beta=beta)
			enddo
			!$omp end parallel do
		case(2)
			!$omp parallel do reduction(+:O)
			do i=1,size(WA,1)
				call gemv(dwf(icfg,i,:),WA(i,:),O,trans="T",beta=beta)
			enddo
			!$omp end parallel do
		end select
	end subroutine
	subroutine get_Spb(cfg,dcfg,Ecfg,nn,pb,WA,iA,AW,WAW)
		integer :: cfg(:),dcfg(:),Ecfg(:),nn(0:,:)
		complex(8) :: WA(:,:),AW(:,:),WAW(:,:),iA(:,:)
		complex(8) :: pb
		complex(8) :: pb2(2),pbs(2)
		integer :: i,l,sg,ct=0
		integer :: k(0:4),m(0:4)
		if(.not.get_row(dcfg,cfg,k,sg,pack([1],[is_project]))) then
			write(*,"(*(i5))")[1:Ns]
			write(*,"(*(i5))")cfg(1:Ns)
			write(*,"(*(i5))")cfg(Ns+1:)
			write(*,"(*(i5))")dcfg
			write(*,"(*(i5))")k(1:k(0))
			error stop "err get_Spb"
		endif
		pbs=0d0
		!$omp parallel do private(pb2,sg,m) reduction(+:pbs)
		do i=1,ubound(nn,1)
			if(get_row([nn(i,:),-nn(0,:)],Ecfg,m,sg,[integer::])) then
				call get_pb(shape(0),m(1:m(0)),pb2(1),WA=WA,iA=iA,AW=AW,WAW=WAW)
				pbs(1)=pbs(1)+pb2(1)*conjg(pb2(1))
				call get_pb(k(1:k(0)),m(1:m(0)),pb2(2),WA=WA,iA=iA,AW=AW,WAW=WAW)
				pbs(2)=pbs(2)+pb2(2)*conjg(pb2(2))
			endif
		enddo
		!$omp end parallel do
		!write(*,*)real(pbs)
		pb=pbs(2)/pbs(1)
	end subroutine
	subroutine get_overlap(cfg,Ecfg,nn,WA,iA,AW,WAW,Ok,Ek)
		integer :: cfg(:),Ecfg(:),nn(0:,:)
		complex(8) :: WA(:,:),AW(:,:),WAW(:,:),iA(:,:),Ek(:),Ok(:)
		complex(8) :: pb
		integer :: i,sg,j
		integer :: m(0:4)
		!$omp parallel do private(pb,sg,m)
		do i=1,ubound(nn,1)
			if(get_row([nn(i,:),-nn(0,:)],Ecfg,m,sg,[integer::])) then
				call get_pb(shape(0),m(1:m(0)),pb,WA=WA,iA=iA,AW=AW,WAW=WAW)
				Ok(i)=conjg(pb*sg)
				Ek(i)=conjg(get_energy(cfg,WA,m=m(1:m(0)),iA=iA,AW=AW,WAW=WAW)*sg)
			endif
		enddo
		!$omp end parallel do
	end subroutine
	real(8) function jast(cfg,dcfg,Oja)
		integer :: cfg(:),dcfg(:)
		complex(8), optional :: Oja(:)
		integer :: i(0:size(dcfg)),j,plv,lv,ibd,l,l1,l2,ic,jc,nb,k(0:size(dcfg)),sg
		real(8) :: n
		jast=0d0
		if(sum(Hja%var(:)%n)==0) then
			return
		endif
		if(present(Oja)) then
			if(sum(Hja%var(1:)%n)==0) return
		endif
		! the index of jastrew facter
		plv=0
		i(0)=0
		i(1:)=Hmf%H2s(abs(dcfg),1)
		do l=1,size(dcfg)
			if(all(i(1:i(0))/=i(l))) then
				i(0)=i(0)+1
				i(i(0))=i(l)
			endif
		enddo
		do l=1,Hja%rg(2)
			nb=Hja%var(l)%nb
			do l1=1,i(0)
				do l2=1,size(latt%nb(nb)%st(i(l1))%j)
					n=0d0
					j=latt%nb(nb)%st(i(l1))%j(l2)
					if(any(i(1:l1-1)==j)) cycle
					do jc=1,size(Hja%var(l)%c,2)
						if(get_row([(Hmf%s2H(merge(i(l1),j,Hja%var(l)%c(ic,jc)%l=="i"),abs(Hja%var(l)%c(ic,jc)%tp))*sign(1,Hja%var(l)%c(ic,jc)%tp),ic=1,size(Hja%var(l)%c,1)),dcfg],cfg,k,sg,[integer::])) then
							n=n+Hja%var(l)%sg(jc)
						endif
						if(get_row([(Hmf%s2H(merge(i(l1),j,Hja%var(l)%c(ic,jc)%l=="i"),abs(Hja%var(l)%c(ic,jc)%tp))*sign(1,Hja%var(l)%c(ic,jc)%tp),ic=1,size(Hja%var(l)%c,1))],cfg,k,sg,[integer::])) then
							n=n-Hja%var(l)%sg(jc)
						endif
					enddo
					ibd=latt%nb(nb)%st(i(l1))%bd(l2)
					lv=Hja%var(l)%bd2v(ibd)
					n=n*Hja%var(l)%bd(ibd)
					if(present(Oja)) then
						Oja(plv+lv)=Oja(plv+lv)+n
					else
						jast=jast+Hja%var(l)%val(lv)*n
					endif
				enddo
			enddo
			plv=plv+Hja%var(l)%n
		enddo
	end function
end module

module vmc_main
	use model
	use lapack95, only : hegv, geev, hegv
	implicit none
	type t_mc
		integer :: sg,ne(2),hot,samp,step
		integer, allocatable :: cfg(:),Ecfg(:),nn(:,:)
		real(8), allocatable :: Emf(:)
		real(8) :: phy(7)
		real(8), allocatable :: phy_n(:,:),phy_sc(:)
		complex(8), allocatable :: S(:,:),Ok2(:,:),Ek2(:,:),psi0(:),wf(:,:),dwf(:,:,:)
		real(8), allocatable :: g(:)
		real(8) :: q(3)
		integer :: num
		integer :: delay
	contains
		procedure :: init
		procedure :: spect
		procedure :: do_mc
		procedure :: do_vmc
		procedure :: do_var
	end type
contains
	subroutine init(self,init_cfg)
		class(t_mc) :: self[ica1,*]
		logical :: init_cfg
		complex(8) :: H(Hmf%Hs,Hmf%Hs),cH(size(H,1),size(H,2)),D(size(H,1),size(H,2),sum(Hmf%var(1:)%n)),psi0_(Ns*Ns),Hq(Hmf%Hi)
		real(8) :: E(size(H,1)),dr(3),tmp
		integer :: l,i,j,n,n1,n2,nn_(0:Ns*Ns,2),kq,ord(Hmf%Hs),Nu
		logical :: check=.true.,is_mix
		
		if(allocated(self%phy_n)) deallocate(self%phy_n,self%phy_sc)
		allocate(self%phy_n(2,Ns))
		allocate(self%phy_sc(size(latt%nb(1)%bd)))

		if(init_cfg) then
			if(allocated(self%cfg)) deallocate(self%cfg)
			allocate(self%cfg(Hmf%Hs))
			do i=1,size(self%ne)
				self%cfg(Hmf%s2H(1:self%ne(i),i))=[sum(self%ne(1:i))-self%ne(i)+1:sum(self%ne(1:i))]
				self%cfg(Hmf%s2H(self%ne(i)+1:Ns,i))=0
			enddo
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
					if(this_image()==1) write(*,*)"warning, close-shell condition is not satisfied, maybe set to:",Ns/2-(self%ne(1)-count((E(ord(self%ne(1)))-E(ord(:self%ne(1))))<1d-7)),"or",Ns/2-(self%ne(1)+count((-E(ord(self%ne(1)))+E(ord(self%ne(1):Ns)))<1d-7))
				endif
				if(icp/=-111) then
					if(this_image()==1) write(*,*)"cp is: ",Hmf%var(icp)%val(1)-Hmf%var(icp)%bd(1)*(E(ord(self%ne(1)))+E(ord(self%ne(1)+1)))/2d0
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
				if(.not.is_in(-brizon%k(i,:)+self%q,brizon%Tc,dr)) then
				endif
				do j=1,brizon%nk
					if(sum(abs(-brizon%k(i,:)+self%q+dr-brizon%k(j,:)))<1d-5) then
						if(kq/=0) then
							if(this_image(self,2)==1) write(*,*)"warn! kq"
						else
							kq=j
						endif
					endif
				enddo
				
				if(kq==0) then
					if(this_image(self,2)==1) write(*,*)self%q(:2)/pi,"warn!! kq2"
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
		self%wf=H


		if(self%sg==2) then
			if(allocated(self%dwf)) deallocate(self%dwf,self%g,self%S)
			allocate(self%dwf(Hmf%Hs,Hmf%Hs,sum(Hmf%var(1:)%n)),self%g(sum(Hmf%var(1:)%n)+sum(Hja%var(1:)%n)))
			allocate(self%S(size(self%g),size(self%g)))
			cH=transpose(conjg(H))
			self%dwf=0d0
			D=0d0
			call Hmf%dHamilton([1:Hmf%rg(2)],H,cH,D)

			do l=1,size(self%dwf,3)
				select case(opt)
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
		select case(self%sg)
		case(2)
			if(ica1/=1) write(*,*)"ica1 is not 1!!!"
			self%phy=0d0; self%S=0d0; self%g=0d0
			call random_number(seed,4285)
			call self%do_mc(seed+7*this_image())
			sync all
			if(this_image(self,2)==1) then
				do l=2,ica2
					!self[this_image(self,1),1]%phy=self[this_image(self,1),1]%phy+self[this_image(self,1),l]%phy; self[this_image(self,1),1]%S=self[this_image(self,1),1]%S+self[this_image(self,1),l]%S; self[this_image(self,1),1]%g=self[this_image(self,1),1]%g+self[this_image(self,1),l]%g
					self%phy=self%phy+self[this_image(self,1),l]%phy; self%S=self%S+self[this_image(self,1),l]%S; self%g=self%g+self[this_image(self,1),l]%g
				enddo
				self%phy=self%phy*1d0/ica2; self%S=self%S*1d0/ica2; self%g=self%g*1d0/ica2
				self%phy(ier)=sqrt(self%phy(ier)-self%phy(iE)**2)
			endif
			sync all
			if(this_image(self,2)/=1) then
				self%phy=self[this_image(self,1),1]%phy; self%S=self[this_image(self,1),1]%S; self%g=self[this_image(self,1),1]%g
			endif
			sync all
			call heev(self%S,eg,'V')
			if(this_image()==1) write(*,"(es12.4$)")Hmf%var(1:)%val(1),self%phy(iE),self%phy(ier)
			if(this_image()==1) write(*,"(i3$)")int(sign(1d0,self%g))
			if(this_image()==1) write(*,"(A$)")' | '
			if(this_image()==1) write(20,"(es12.4$)")Hmf%var(1:)%val(1),self%phy(iE),self%phy(ier),self%g
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
			if(this_image()==1) write(*,"(i3$)")int(sign(1d0,self%g))
			!grad=g
			self%g=self%g*Ns
			!stop
		case(1)
			self%phy=0d0;self%phy_sc=0d0;self%phy_n=0d0
			call random_number(seed,4285)
			call self%do_mc(seed+7*this_image())
			sync all
			if(this_image(self,2)==1) then
				do l=2,ica2
					!self[this_image(self,1),1]%phy=self[this_image(self,1),1]%phy+self[this_image(self,1),l]%phy;
					self%phy=self%phy+self[this_image(self,1),l]%phy;
					self%phy_sc=self%phy_sc+self[this_image(self,1),l]%phy_sc;
					self%phy_n=self%phy_n+self[this_image(self,1),l]%phy_n;
				enddo
				self%phy=self%phy*1d0/ica2
				self%phy_sc=self%phy_sc*1d0/ica2
				self%phy_n=self%phy_n*1d0/ica2
				self%phy(ier)=sqrt(self%phy(ier)-self%phy(iE)**2)
			endif
			sync all
			if(this_image(self,2)/=1) then
				self%phy=self[this_image(self,1),1]%phy
				self%phy_sc=self[this_image(self,1),1]%phy_sc
				self%phy_n=self[this_image(self,1),1]%phy_n
			endif
			sync all
			!critical
			if(this_image(self,2)==1) write(*,"(*(es12.4))")self%q(:2)/pi,Hmf%var(1:)%val(1),self%phy(iE),self%phy(ier),self%phy(iSq_pm)
			if(this_image(self,2)==1) then
				write(20,"(*(es12.4))")Hmf%var(1:)%val(1)
				do k=1,Ns
					write(30,"(*(es12.4))")latt%nb(0)%bd(k)%r,self%phy_n(:,k)
				enddo
				write(30,"(x/)")
				do k=1,size(latt%nb(1)%bd)
					write(30,"(*(es12.4))")latt%nb(1)%bd(k)%r,latt%nb(1)%bd(k)%dr,sign(1d0,self%phy_sc(k))*sqrt(abs(self%phy_sc(k)))
				enddo
			endif
			!endcritical
		case(3)
			self%Ok2=0d0; self%Ek2=0d0
			call random_number(seed,4285)
			call self%do_mc(seed+7*this_image())
			sync all
			if(this_image(self,2)==1) then
				do l=2,ica2
					!self[this_image(self,1),1]%Ok2=self[this_image(self,1),1]%Ok2+self[this_image(self,1),l]%Ok2; self[this_image(self,1),1]%Ek2=self[this_image(self,1),1]%Ek2+self[this_image(self,1),l]%Ek2
					self%Ok2=self%Ok2+self[this_image(self,1),l]%Ok2; self%Ek2=self%Ek2+self[this_image(self,1),l]%Ek2
				enddo
				self%Ok2=self%Ok2*1d0/ica2; self%Ek2=self%Ek2*1d0/ica2
			endif
			sync all
			if(this_image(self,2)/=1) then
				self%Ok2=self[this_image(self,1),1]%Ok2; self%Ek2=self[this_image(self,1),1]%Ek2
			endif
			sync all
		end select
		self%samp=self%samp*ica2
	end subroutine
	subroutine do_mc(self,seed)
		class(t_mc) :: self[ica1,*]
		integer :: seed
		complex(8) :: iA(sum(self%ne),sum(self%ne)),iAl(size(iA,1),self%delay),A(size(iA,1),size(self%wf,2)),WA(size(self%wf,1),size(iA,2)),WAl(size(WA,1),size(iAl,2)),WAr(size(iAl,2),size(WA,2)),AW(size(iA,1),size(self%wf,2)),WAW(size(self%wf,1),size(self%wf,2)),AWr(size(iAl,2),size(AW,2)),wf(size(self%wf,1),size(self%wf,2)),dwf(size(self%dwf,1),size(self%dwf,2),size(self%dwf,3))
		integer :: cfgl(size(self%cfg)),Ecfgl(size(self%Ecfg)),icfg(size(A,1)),iEcfg(size(A,1))
		real(8) :: gp(size(self%g)),gl(size(gp))
		complex(8) :: Sp(size(gp),size(gp)),Op(size(gp)),Ol(size(gp)),Sl(size(gp),size(gp))
		complex(8) :: Ok(size(self%Ok2,1)),Ek(size(Ok)),Ok2p(size(Ok),size(Ok)),Ek2p(size(Ok),size(Ok))
		integer :: nn(lbound(self%nn,1):ubound(self%nn,1),size(self%nn,2))
		complex(8) :: phyl(size(self%phy)),phyp(size(self%phy))
		real(8) :: phyl_sc(size(self%phy_sc)),phyp_sc(size(self%phy_sc)),phyl_n(2,size(self%phy_n)),phyp_n(2,size(self%phy_n))
		real(8) :: rd,isamp,Oq
		complex(8) :: pb
		integer :: i,j,l,l1,l2,info,n,apt,samp,sg,dly,np
		integer :: k(0:4),m(0:4),dcfg(0:4)
		logical :: is_update,is_full,is_accept
		type(randomNumberSequence) :: rnd
		call mt_init_random_seed(rnd,seed)
		phyp=0d0; phyl=0d0;phyl_sc=0d0;phyp_sc=0d0;phyl_n=0d0;phyp_n=0d0
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

				rd=jast(cfgl,dcfg(1:dcfg(0)),Ol(sum(Hmf%var(1:)%n)+1:))

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
				if(self%sg==2.and.opt==1) then
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
						phyl(iE)=get_energy(cfgl,WA)*1d0/Ns
						!phyl%dsc=real(get_dsc(cfgl,WA))
						phyl(iSq_pm)=get_spin_pm(cfgl,WA,self%q)
						!phyl(iSq_zz)=get_spin_zz(cfgl,q)
						!call get_phy_n(cfgl,phyl_n)
						!call get_phy_sc(cfgl,WA,phyl_sc)
					case(2)
						phyl(iE)=get_energy(cfgl,WA)*1d0/Ns
						call get_O(icfg,iEcfg,WA,iA,dwf,Ol(:size(dwf,3)))
					case(3)
						call get_overlap(cfgl,Ecfgl,nn,WA,iA,AW,WAW,Ok,Ek)
						Oq=1d0/sum(real(Ok*conjg(Ok)))
					end select
				endif
				phyp=phyp+phyl
				phyp_sc=phyp_sc+phyl_sc
				phyp_n=phyp_n+phyl_n
				select case(self%sg)
				case(2)
					do l1=1,size(Ol)
						do l2=1,size(Ol)
							Sp(l2,l1)=Sp(l2,l1)+conjg(Ol(l2))*Ol(l1)
						enddo
					enddo
					gp=gp+real(phyl(iE)*conjg(Ol))
					Op=Op+Ol
				case(3)
					do l1=1,size(Ok)
						do l2=1,size(Ok)
							Ek2p(l2,l1)=Ek2p(l2,l1)+Ok(l2)*conjg(Ek(l1))*Oq
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

		phyp=phyp*isamp;phyp_sc=phyp_sc*isamp;phyp_n=phyp_n*isamp
		Sp=Sp*isamp; gp=gp*isamp; Op=Op*isamp
		Ek2p=Ek2p*isamp
		Ok2p=Ok2p*isamp
		do l1=1,size(Op)
			do l2=1,size(Op)
				!S(l1,l2)=2d0*SE(l1,l2)-O(l1)*g(l2)-O(l2)*g(l1)-2d0*E*S(l1,l2) ! maybe works
				!S(l1,l2)=2d0*(SE(l1,l2)-S(l1,l2)*E-O(l1)*g(l2)-O(l2)*g(l1))
				Sp(l1,l2)=Sp(l1,l2)-conjg(Op(l1))*Op(l2)
			enddo
		enddo
		gp=2d0*(gp-real(phyp(iE)*conjg(Op)))

		self%cfg=cfgl
		select case(self%sg)
		case(1:2)
			self%phy=real(phyp)
			self%phy(ier)=real(phyp(iE))**2
			self%phy_sc=phyp_sc;self%phy_n=phyp_n
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
			Sq(l)=Sq(l)+sum(psi0_(n0:)/(omg(1)+domg*l-Eg(n0:)+self%phy(iE)*Ns+img*gm))
		enddo
		norm=self%phy(iSq_pm)/real(sum(psi0_(n0:)))
		write(*,"(es12.4$)")-sum(imag(Sq)*domg)/pi,real(sum(psi0_(n0:))),self%phy(iSq_pm)
		write(ut,"('#q=',3es12.4,2i5)")self%q/pi,i,size(EO)
		do l=1,m
			if(Smax(2)<abs(imag(Sq(l)))) Smax=[omg(1)+domg*l,abs(imag(Sq(l)))]
			write(ut,"(*(es17.9))")omg(1)+domg*l,Sq(l)*norm
		enddo
		write(ut,"(x)")
		write(*,*)"max intensity: ",Smax(1)

		if(present(E)) then
			write(*,*)sum(psi0_(n0:)/(E-Eg(n0:)+self%phy(iE)*Ns+img*gm))
			!psi0_=matmul(self%Ok2,self%psi0)
			pnn=0d0
			psi0_=matmul(transpose(conjg(H_)),matmul(self%Ok2,self%psi0))
			do l=n0,size(Eg)
				
				!pnn=pnn+conjg(H_(:,l))*psi0_*sum(H_(:,l)*conjg(psi0_))/(E-Eg(l)+self%phy(iE)*Ns+img*gm)
				!pnn=pnn+H_(:,l)*conjg(H_(:,l))*psi0_(l)/(E-Eg(l)+self%phy(iE)*Ns+img*gm)
				!pnn=pnn+H_(:,l)*conjg(H_(:,l))*psi0_(l)/(E-Eg(l)+self%phy(iE)*Ns+img*gm)
				pnn=pnn+conjg(matmul(self%Ok2,H_(:,l))*psi0_(l))*matmul(self%Ok2,H_(:,l))*psi0_(l)/(E-Eg(l)+self%phy(iE)*Ns+img*gm)
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
		real(8) :: x(sum(Hmf%var(1:)%n),n),dx,er,pgrad(size(x)),El(n)
		integer :: i,hot,ord(n)
		logical :: init_cfg
		self%sg=2
		!dx=0.03d0
		dx=0.3d0
		hot=self%hot
		i=0
		if(size(x)==0) then
			call self%init(.true.)
			call self%do_vmc()
			return
		endif
		x(:,1)=Hmf%put([1:Hmf%rg(2)])
		init_cfg=.true.
		do 
			i=i+1
			if(allocated(self%g)) pgrad=self%g
			call Hmf%get(x(:,i),[1:Hmf%rg(2)])
			if(this_image()==1) write(*,"(i4$)")i
			call self%init(init_cfg)
			call self%do_vmc()
			El(i)=self%phy(iE)
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
			self%hot=128
			init_cfg=.false.
			if(this_image()==1) write(*,"(x)")
			if(this_image()==1) write(20,"(x)")
		enddo
		self%hot=hot
		if(this_image()==1) write(*,*)"finished"
		call qsort(El,ord)
		do i=2,min(n/5,30)
			x(:,ord(1))=x(:,ord(1))+x(:,ord(i))
		enddo
		call Hmf%get(x(:,ord(1))/min(n/5,30),[1:Hmf%rg(2)])
		!E=minval(El(:i))
	end subroutine
end module

program main
	use vmc_main
	implicit none
	logical :: f
	integer :: i,j,l
	type(t_mc) :: mc[ica1,*]
	real(8) :: et,E
	character(10) :: hostname
	i=hostnm(hostname)
	ica2=num_images()/ica1
	if(num_images()/=ica2*ica1) stop "plz check the coarray number"
	if(this_image()==1) then
		open(101,file="../data/lattice.dat")
		open(111,file="../data/tmp.dat")
		open(20,file="../data/var.dat")
		open(30,file="../data/phyri.dat")
		open(40,file="../data/phyvar.dat")
		open(50,file="../data/matrix.dat",access="stream")
		open(70,file="../data/spect.dat")
		open(71,file="../data/spect_kmap.dat")
		open(80,file="../data/band.dat")
		write(*,*)"coarray info: ",ica1,ica2,"out of",num_images()

	endif
	sync all

	write(*,*)"runing ",this_image()," in ",hostname

	call initial()
	if(this_image()==1) then
		write(50)size(brizon%k,1),brizon%k,brizon%nk
	endif
	otime=0d0


	call omp_set_nested(.false.)
	!call omp_set_max_active_levels(1)
	call omp_set_dynamic(.false.)
	call mkl_set_dynamic(0)
	call omp_set_schedule(omp_sched_static,0)

	call mkl_set_num_threads(1)
	call omp_set_num_threads(mkl_get_max_threads())


	mc%hot=1
	mc%step=Ns
	mc%delay=2
	mc%ne(1)=Ns/2
	mc%ne(2)=Ns-mc%ne(1)
	!Hmf%var(1:vn)%val(1)=[-6.2487E-01,2.5764E-01,1.0392E-01,-3.6208E-02] ! dsc+mu+t'+SDW E=-3.2721E-01
	mc%samp=1024*8*16*8
	mc%hot=1024*8*2
	mc%step=nint(sqrt(real(Ns)))
	!call mc%do_var(100)
	!stop
	mc%q=q(:,this_image(mc,1))
	mc%num=this_image(mc,1)

	mc%sg=1
	mc%samp=1024*8*8
	mc%hot=1024*8*8
	mc%step=nint(sqrt(real(Ns)))*2

	call mc%init(.true.)
	call mc%do_vmc()
	stop

	if(this_image()==1) then
		do l=2,ica1
			mc[1,1]%phy(iE)=mc[1,1]%phy(iE)+mc[l,1]%phy(iE)
		enddo
		mc[1,1]%phy(iE)=mc[1,1]%phy(iE)*1d0/ica1
		sync images(*)
	else
		sync images(1)
		mc%phy(iE)=mc[1,1]%phy(iE)
	endif
	sync all

	mc%sg=3
	mc%ne=mc%ne+1
	mc%hot=1024*8*8
	mc%samp=1024*8*8*8 !dsc 16x16
	mc%step=nint(sqrt(real(Ns)))*2
	call mc%init(.true.)
	call mc%do_vmc()

	if(this_image()==1) then
		write(*,*)"finished, exporting data....",mc%samp
		do l=1,ica1
			write(50)-1,size(mc[l,1]%psi0),Ns,mc[l,1]%q,mc[l,1]%phy(iE),mc[l,1]%phy(iSq_pm),mc[l,1]%Emf,mc[l,1]%psi0,mc[l,1]%Ok2,mc[l,1]%Ek2,mc[l,1]%nn(1:,:)
			!call mc%spect(70,0.02d0,[-10d0,10d0],5000)
			!mc[l,1]%Ek2=transpose(conjg(mc[l,1]%Ek2))
			!call mc%spect(70,0.02d0,[-10d0,10d0],5000)
		enddo
		rewind(50)
		read(50)i
		deallocate(brizon%k)
		allocate(brizon%k(i,3))
		read(50)brizon%k
		read(50)brizon%nk
		do l=1,ica1
			read(50)i,j,Ns
			if(allocated(mc%psi0)) deallocate(mc%psi0,mc%Ok2,mc%Ek2,mc%Emf,mc%nn)
			allocate(mc%psi0(j),mc%Ok2(j,j),mc%Ek2(j,j),mc%nn(j,2),mc%Emf(Ns*2))
			read(50)mc%q,mc%phy(iE),mc%phy(iSq_pm),mc%Emf,mc%psi0,mc%Ok2,mc%Ek2,mc%nn
			write(*,*)mc%q,mc%phy(iE),Ns,i,j
			call mc%spect(70,0.02d0,[-10d0,10d0],5000)
			!call mc%spect(70,0.1d0,[-10d0,10d0],5000)
			mc%Ek2=transpose(conjg(mc%Ek2))
			call mc%spect(70,0.02d0,[-10d0,10d0],5000)
			!call mc%spect(70,0.1d0,[-10d0,10d0],5000)
			mc%Ek2=0.5d0*(mc%Ek2+transpose(conjg(mc%Ek2)))
			call mc%spect(70,0.02d0,[-10d0,10d0],5000)
			!call mc%spect(70,0.1d0,[-10d0,10d0],5000)
			deallocate(mc%psi0,mc%Ok2,mc%Ek2,mc%Emf,mc%nn)
		enddo
	endif
end program
