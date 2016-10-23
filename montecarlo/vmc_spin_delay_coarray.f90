module global
	use M_hamilton_test1
	use ifport, ifort_qsort => qsort
	use M_omp_rand
	use omp_lib
	use mkl_service
	implicit none
	real(8) :: t(2)=(/1d0,-0.3d0/)
	real(8), parameter :: DJ=0.3d0,V=0d0,U=0d0
	integer, parameter :: opt=2 ! 1: Tr(AdA)
								! 2: O=cicj
	integer, parameter :: iE=1,ier=2,idsc=3,iaf=4,iddw=5,iSq_pm=6,iSq_zz=7,iCq=8
	integer, parameter :: ica1=4
	integer :: ica2
	integer :: Ns,vn
	real(8) :: q(3,ica1)=reshape([&
		!pi*5d0/8d0,pi*5d0/8d0,0d0,&
		!pi*6d0/8d0,pi*6d0/8d0,0d0,&
		!pi*7d0/8d0,pi*7d0/8d0,0d0,&
		!pi,pi,0d0&
		pi,pi*7d0/8d0,0d0,&
		pi,pi*6d0/8d0,0d0,&
		pi,pi*5d0/8d0,0d0,&
		pi,pi*4d0/8d0,0d0&
		],[3,ica1])
contains
	subroutine initial()
		integer :: i,l
		real(8) :: q(3)
		q=(/1d0/8d0,0d0,0d0/)*2d0*pi
		allocate(var(-50:50))
		call init_random_seed()
		! lattice
		latt%a1=(/1d0,0d0,0d0/)
		latt%a2=(/0d0,1d0,0d0/)
		latt%c1=(/8d0,0d0,0d0/)
		latt%c2=(/0d0,2d0,0d0/)
		!latt%c1=(/1d0,1d0,0d0/)
		!latt%c2=(/-1d0,1d0,0d0/)
		!latt%c1=latt%a1
		!latt%c2=latt%a2
		latt%T1=(/1d0,0d0,0d0/)*16
		latt%T2=(/0d0,1d0,0d0/)*16
		latt%bdc=(/1d0,1d0,0d0/)
		allocate(latt%rsb(1,3))
		latt%rsb(1,:)=(/0d0,0d0,0d0/)
		latt%n1=1
		latt%n2=1
		!call latt%gen_latt(size(t))
		call latt%gen_latt(3)
		call latt%gen_brizon(brizon)
		!call check_lattice(101)
		if(this_image()==1) write(*,*)"Total site number is: ",latt%Ns

		! cp
		call gen_var(sg=1,nb=0)
		var(iv(0))%bd=-1d0
		var(iv(0))%val=-1.39752d0

		! dsc
		call gen_var(sg=-2,nb=1)
		do i=1,size(var(iv(0))%bd)
			var(iv(0))%bd(i)=dwave(i)
		enddo
		var(iv(0))%val=1d-4

		!! ssc
		!call gen_var(sg=-2,nb=0)
		!do i=1,size(var(iv(0))%bd)
			!var(iv(0))%bd(i)=1d0
		!enddo
		!var(iv(0))%val=1d-3

		!! ddw
		!call gen_var(sg=3,nb=1)
		!do i=1,size(var(iv(0))%bd)
			!var(iv(0))%bd(i)=ab(latt%nb(var(iv(0))%nb)%bd(i)%i(1))*dwave(i)*img
		!enddo
		!var(iv(0))%val=1d-1

		! sdw
		call gen_var(sg=4,nb=0)
		do i=1,size(var(iv(0))%bd)
			var(iv(0))%bd(i)=ab(i)*sin(sum(q*(latt%nb(var(iv(0))%nb)%bd(i)%r)))
		enddo
		var(iv(0))%val=1d-1

		! on site cdw
		call gen_var(sg=3,nb=0)
		do i=1,size(var(iv(0))%bd)
			var(iv(0))%bd(i)=cos(sum(2d0*q*(latt%nb(var(iv(0))%nb)%bd(i)%r)))
		enddo
		var(iv(0))%val=1.32678d-2

		do l=2,size(t)
			call gen_var(sg=3,nb=l)
			var(iv(0))%bd=-1d0
			var(iv(0))%val=-2.00021d-1
		enddo

		do l=1,size(t)
			call gen_var(sg=-3,nb=l)
			var(iv(0))%bd=-1d0
			var(iv(0))%val=t(l)
		enddo

		vn=iv(1)

		!! n-n jast
		!do i=1,3
			!call gen_var(sg=9,nb=i)
			!var(iv(0))%bd=-1d0
			!var(iv(0))%val=1d0
		!enddo

		!! s-s jast
		!do i=1,1
			!call gen_var(sg=10,nb=i)
			!var(iv(0))%bd=1d0
			!var(iv(0))%val=1d0
		!enddo

		call var_shrink()

		Ns=latt%Ns
	end subroutine
	real(8) function ab(i)
		integer :: i
		ab=(-1d0)**(mod(nint(sum(latt%nb(0)%bd(i)%r(:2))),2))
	end function
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
		k_=(/k(1::2),m(1::2)/)
		m_=(/k(2::2),m(2::2)/)
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
			k_=(/k(1::2)/)
			m_=(/k(2::2)/)
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
				A=-1d0/(A(1,1)*A(2,2)-A(2,1)*A(1,2))*reshape((/A(2,2),-A(2,1),-A(1,2),A(1,1)/),(/2,2/))
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
	subroutine change(cfg,icfg,dcfg,k,rnd)
		integer :: cfg(:),icfg(:)
		integer :: dcfg(0:),k(0:)
		type(randomNumberSequence) :: rnd
		integer :: i,j,ii,jj
		!call random_number(i,size(icfg))
		i=mt_random_number(rnd,size(icfg))
		i=icfg(i)
		ii=i+(1-(i-1)/Ns*2)*Ns
		do 
			!call random_number(j,Ns)
			j=mt_random_number(rnd,Ns)
			j=j+(i-1)/Ns*Ns
			jj=j+(1-(j-1)/Ns*2)*Ns
			if(cfg(j)==0) then
				if(sign(1,-cfg(jj))==sign(1,-cfg(ii))) then
					dcfg(0:2)=(/2,i,-j/)
					k(0:2)=(/2,cfg(dcfg(1)),abs(dcfg(2))/)
				else
					dcfg(0:4)=(/4,i,-j,ii,-jj/)
					k(0:4)=(/4,cfg(dcfg(1)),abs(dcfg(2)),cfg(dcfg(3)),abs(dcfg(4))/)
				endif
				exit
			endif
		enddo
	end subroutine
	logical function get_row(cfg,dcfg,k,sg,P)
		integer, intent(in) :: cfg(:),dcfg(:)
		integer, intent(out) :: k(0:),sg
		integer, intent(in), optional :: P(:)
		integer :: i,j,n,l
		logical :: db(0:1),dcfg_(size(dcfg)),iP(size(dcfg))
		sg=1
		get_row=.false.
		do i=1,size(dcfg)+1
			if(present(P)) then
				if(any(P==i).or.i==size(dcfg)+1) then
					do j=1,i-1
						if(dcfg_(j)) then
							l=sign(1,abs(dcfg(j))-1-Ns)
							db((l+1)/2)=(cfg(abs(dcfg(j)))==0)
							db(mod((l+1)/2+1,2))=((cfg(abs(dcfg(j))-l*Ns)/=0)/=iP(j))
							if(db(0).and.(.not.db(1))) return
						endif
					enddo
				endif
			endif
			if(i==size(dcfg)+1) exit
			l=-sign(1,abs(dcfg(i))-1-Ns)*Ns
			iP(i)=.false.
			dcfg_(i)=.true.
			do j=i-1,1,-1
				if(dcfg_(j)) then
					if(dcfg(i)==-dcfg(j)) then
						if(mod(count(dcfg_(j+1:i-1)),2)==1) sg=-sg
						dcfg_(i)=.false.
						dcfg_(j)=.false.
					elseif(dcfg(i)==dcfg(j)) then
						get_row=.false.
						return
					elseif(abs(dcfg(i))+l==(abs(dcfg(j)))) then
						iP(i)=.not.iP(i)
						iP(j)=.not.iP(j)
					endif
				endif
			enddo
			if(dcfg_(i)) then
				if(sign(1,-cfg(abs(dcfg(i))))*dcfg(i)>0) return
			endif
		enddo
		n=0
		get_row=.true.
		do i=1,size(dcfg)
			if(dcfg_(i)) then
				do j=i,size(dcfg)
					if(dcfg_(j)) then
						if(mod(n,2)==0.and.dcfg(j)>0) then
							n=n+1
							k(n)=cfg(dcfg(j))
						elseif(mod(n,2)/=0.and.dcfg(j)<0) then
							n=n+1
							k(n)=abs(dcfg(j))
						else
							cycle
						endif
						if(i/=j) then
							if(mod(n,2)==0.and.dcfg(i)>0) then
								n=n+1
								k(n)=cfg(dcfg(i))
							elseif(mod(n,2)/=0.and.dcfg(i)<0) then
								n=n+1
								k(n)=abs(dcfg(i))
							endif
							dcfg_(j)=.false.
							if(mod(count(dcfg_(i:j-1)),2)==1) sg=-sg
						endif
						exit
					endif
				enddo
			endif
		enddo
		k(0)=n
	end function
	complex(8) function get_spin_pm(cfg,D,q,ja)
		complex(8) :: D(:,:)
		integer :: cfg(:)
		real(8) :: q(3),ja(:)
		complex(8) :: pb
		integer :: i,j,sg,l
		integer :: k(0:4)
		get_spin_pm=0d0
		!$omp parallel do collapse(2) reduction(+:get_spin_pm) private(pb,k,sg)
		do i=1,Ns
			do j=1,Ns
				if(get_row(cfg,(/i,i+Ns,-j-Ns,-j/),k,sg,shape(0))) then
					call get_pb(k(1:k(0)),shape(0),pb,D)
					get_spin_pm=get_spin_pm+pb*exp(img*sum(q*(latt%nb(0)%bd(i)%r-latt%nb(0)%bd(j)%r)))*sg*exp(jast(cfg,(/i,i+Ns,-j-Ns,-j/),ja))
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
				get_spin_zz=get_spin_zz+0.25d0*sum(((1-sign(1,-cfg(i::Ns)))/2-(/0,1/)))*sum(((1-sign(1,-cfg(j::Ns)))/2-(/0,1/)))*exp(img*sum(q*(latt%nb(0)%bd(i)%r-latt%nb(0)%bd(j)%r)))
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
			get_af=get_af+sum(((1-sign(1,-cfg(i::Ns)))/2-(/0,1/)))*exp(img*sum(q*latt%nb(0)%bd(i)%r))
		enddo
		!$omp end parallel do
		get_af=0.5d0*get_af/Ns
	end function
	complex(8) function get_dsc(cfg,WA,ja,dcfg)
		complex(8) :: WA(:,:)
		integer :: cfg(:)
		real(8) :: ja(:)
		integer, optional :: dcfg(:)
		complex(8) :: pb
		integer :: i1,j1,i2,j2,sg,l,n1,n2,p1,p2
		integer :: k(0:6),dcfg_(0:2)
		if(present(dcfg)) then
			dcfg_(0:)=(/size(dcfg),dcfg/)
		else
			dcfg_(0)=0
		endif
		get_dsc=0d0
		!$omp parallel do collapse(2) reduction(+:get_dsc) private(pb,k,sg,i1,j1,i2,j2)
		do n1=1,size(latt%nb(1)%bd)
			do n2=1,size(latt%nb(1)%bd)
				if(n1==n2) cycle
				i1=latt%nb(1)%bd(n1)%i(1)
				j1=latt%nb(1)%bd(n1)%i(2)
				i2=latt%nb(1)%bd(n2)%i(1)
				j2=latt%nb(1)%bd(n2)%i(2)
				do p1=1,2
					do p2=1,2
						if(get_row(cfg,(/-j1-Ns,i1,-i2,j2+Ns,dcfg_(1:dcfg_(0))/),k,sg,(/5/))) then
							call get_pb(k(1:k(0)),shape(0),pb,WA)
							get_dsc=get_dsc+pb*sg*dwave(n1)*dwave(n2)*latt%nb(1)%bd(n1)%bdc*latt%nb(1)%bd(n2)%bdc*exp(jast(cfg,(/-j1-Ns,i1,-i2,j2+Ns,dcfg_(1:dcfg_(0))/),ja))
						endif
						call swap(i2,j2)
					enddo
					call swap(i1,j1)
				enddo
			enddo
		enddo
		!$omp end parallel do
		get_dsc=get_dsc/(Ns**2)
	end function
	complex(8) function get_energy(cfg,WA,ja,dcfg,m,iA,AW,WAW)
		complex(8) :: WA(:,:)
		integer :: cfg(:)
		real(8) :: ja(:)
		integer, optional :: dcfg(:),m(:)
		complex(8), optional :: AW(:,:),WAW(:,:),iA(:,:)
		complex(8) :: pb
		integer :: i,j,n,l,c(2,2),p,sg
		integer :: k(0:4),dcfg_(0:2)
		logical :: flag
		if(present(dcfg)) then
			dcfg_(0:)=(/size(dcfg),dcfg/)
		else
			dcfg_(0)=0
		endif
		get_energy=0d0
		!$omp parallel do collapse(2) reduction(+:get_energy) private(i,j,c,k,pb,sg)
		do l=1,size(t)
			do n=1,size(latt%nb(l)%bd)
				i=latt%nb(l)%bd(n)%i(1)
				j=latt%nb(l)%bd(n)%i(2)
				do p=1,2
					!if(get_row(cfg,(/i,-j,dcfg_(1:dcfg_(0))/),k,sg,(/3/))) then
					if(cfg(i)/=0.and.cfg(j)==0.and.cfg(j+Ns)/=0.and.cfg(i+Ns)/=0) then
						sg=1
						k(0:2)=(/2,cfg(i),j/)
						if(present(m)) then
							call get_pb(k(1:k(0)),m,pb,WA=WA,iA=iA,AW=AW,WAW=WAW)
						else
							call get_pb(k(1:k(0)),shape(0),pb,WA=WA)
						endif
						get_energy=get_energy-sg*t(l)*pb*latt%nb(l)%bd(n)%bdc*exp(jast(cfg,(/i,-j,dcfg_(1:dcfg_(0))/),ja))
					!elseif(get_row(cfg,(/-i-Ns,j+Ns,dcfg_(1:dcfg_(0))/),k,sg,(/3/))) then
					elseif(cfg(i+Ns)==0.and.cfg(j+Ns)/=0.and.cfg(i)==0.and.cfg(j)==0) then
						sg=-1
						k(0:2)=(/2,cfg(j+Ns),i+Ns/)
						if(present(m)) then
							call get_pb(k(1:k(0)),m,pb,WA=WA,iA=iA,AW=AW,WAW=WAW)
						else
							call get_pb(k(1:k(0)),shape(0),pb,WA=WA)
						endif
						get_energy=get_energy-sg*t(l)*pb*conjg(latt%nb(l)%bd(n)%bdc)*exp(jast(cfg,(/-i-Ns,j+Ns,dcfg_(1:dcfg_(0))/),ja))
					elseif(l==1) then
						!if(get_row(cfg,(/i,i+Ns,-j-Ns,-j,dcfg_(1:dcfg_(0))/),k,sg,(/5/))) then
						if(cfg(i)/=0.and.cfg(i+Ns)/=0.and.cfg(j+Ns)==0.and.cfg(j)==0) then
							sg=1
							k(0:4)=(/4,cfg(i),j,cfg(i+Ns),j+Ns/)
							if(present(m)) then
								call get_pb(k(1:k(0)),m,pb,WA=WA,iA=iA,AW=AW,WAW=WAW)
							else
								call get_pb(k(1:k(0)),shape(0),pb,WA=WA)
							endif
							get_energy=get_energy+0.5d0*DJ*pb*sg*exp(jast(cfg,(/i,i+Ns,-j-Ns,-j,dcfg_(1:dcfg_(0))/),ja))
						endif
					endif
					j=latt%nb(l)%bd(n)%i(1)
					i=latt%nb(l)%bd(n)%i(2)
				enddo
				if(l==1) then
					!if(get_row(cfg,dcfg_(1:dcfg_(0)),k,sg,shape(0))) then
						k(0)=0
						sg=1
						if(present(m)) then
							call get_pb(k(1:k(0)),m,pb,WA=WA,iA=iA,AW=AW,WAW=WAW)
						else
							call get_pb(k(1:k(0)),shape(0),pb,WA=WA)
						endif
						c(:,1)=(1-sign(1,-cfg(i::Ns)))/2-(/0,1/)
						c(:,2)=(1-sign(1,-cfg(j::Ns)))/2-(/0,1/)
						get_energy=get_energy+pb*sg*(V*sum(c(:,1)*(/1,-1/))*sum(c(:,2)*(/1,-1/))+0.25d0*DJ*sum(c(:,1))*sum(c(:,2)))*exp(jast(cfg,dcfg_(1:dcfg_(0)),ja))
					!endif
				endif
			enddo
		enddo
		!$omp end parallel do
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
		if(.not.get_row(cfg,dcfg,k,sg,shape(0))) then
			write(*,*)cfg
			write(*,*)dcfg
			write(*,*)k(1:k(0))
			stop "err get_Spb"
		endif
		pbs=0d0
		!$omp parallel do private(pb2,sg,m) reduction(+:pbs)
		do i=1,ubound(nn,1)
			if(get_row(Ecfg,(/-nn(0,:),nn(i,:)/),m,sg)) then
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
	subroutine get_overlap(cfg,Ecfg,nn,WA,iA,AW,WAW,Ok,Ek,ja)
		integer :: cfg(:),Ecfg(:),nn(0:,:)
		real(8) :: ja(:)
		complex(8) :: WA(:,:),AW(:,:),WAW(:,:),iA(:,:),Ek(:),Ok(:)
		complex(8) :: pb
		integer :: i,sg,j
		integer :: m(0:4)
		!$omp parallel do private(pb,sg,m)
		do i=1,ubound(nn,1)
			if(get_row(Ecfg,(/-nn(0,:),nn(i,:)/),m,sg)) then
				call get_pb(shape(0),m(1:m(0)),pb,WA=WA,iA=iA,AW=AW,WAW=WAW)
				Ok(i)=conjg(pb*sg)
				Ek(i)=conjg(get_energy(cfg,WA,ja,m=m(1:m(0)),iA=iA,AW=AW,WAW=WAW)*sg)
			endif
		enddo
		!$omp end parallel do
	end subroutine
	real(8) function jast(cfg,dcfg,ja,jan)
		integer :: cfg(:),dcfg(:)
		real(8), optional :: ja(:)
		complex(8), optional :: jan(:)
		integer :: ncfg(size(dcfg),3),ni,nja,nb,i,j,l,ii,jj,dn,n,njap,sg(2)
		real(8) :: jabd
		jast=0d0
		if(present(ja)) then
			if(size(ja)<1) return
		endif
		if(present(jan)) then
			if(size(jan)<1) return
		endif
		sg=1
		if((var(vn+1)%sg==9.and.sg(2)==-1).or.(var(vn+1)%sg==10.and.sg(2)==1)) sg(2)=-sg(2)
		nja=0
		do l=vn+1,size(var(1:))
			nb=var(l)%nb
			if((var(l)%sg==9.and.sg(2)/=-1).or.(var(l)%sg==10.and.sg(2)/=1)) then
				sg(2)=-sg(2)
				ni=0
			j1:	do i=1,size(dcfg)
					j=(mod(abs(dcfg(i))-1,Ns)+1)
					dn=-sign(1,dcfg(i))*sg(((abs(dcfg(i))-1)/Ns+1))
					do n=1,ni
						if(ncfg(n,1)==j) then
							ncfg(n,3)=ncfg(n,3)+dn
							cycle j1
						endif
					enddo
					ni=ni+1
					ncfg(ni,1)=j
					ncfg(ni,2)=sum(((1-sign(1,-cfg(j::Ns)))/2-(/0,1/))*sg)
					ncfg(ni,3)=dn
				enddo j1
			endif
			if(all(ncfg(1:ni,3)==0)) then
				nja=nja+var(l)%n
				cycle
			endif
			do i=1,ni
				ii=ncfg(i,1)
			oj:	do j=1,size(latt%nb(nb)%st(ii)%j)
					jj=latt%nb(nb)%st(ii)%j(j)
					njap=nja+var(l)%bd2v(latt%nb(nb)%st(ii)%bd(j))
					jabd=real(var(l)%bd(latt%nb(nb)%st(ii)%bd(j)))
					do n=1,ni
						if(ncfg(n,1)==jj) then
							if(present(ja)) then
								jast=jast+0.5d0*jabd*ja(njap)*(ncfg(i,3)*ncfg(n,3)+ncfg(i,2)*ncfg(n,3)+ncfg(i,3)*ncfg(n,2))
							endif
							if(present(jan)) then
								jan(njap)=jan(njap)+0.5d0*jabd*(ncfg(i,3)*ncfg(n,3)+ncfg(i,2)*ncfg(n,3)+ncfg(i,3)*ncfg(n,2))
							endif
							cycle oj
						endif
					enddo
					if(present(ja)) then
						jast=jast+jabd*ja(njap)*ncfg(i,3)*sum(((1-sign(1,-cfg(jj::Ns)))/2-(/0,1/))*sg)
					endif
					if(present(jan)) then
						jan(njap)=jan(njap)+jabd*ncfg(i,3)*sum(((1-sign(1,-cfg(jj::Ns)))/2-(/0,1/))*sg)
					endif
				enddo oj
			enddo
			nja=nja+var(l)%n
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
	type, extends(t_sort) :: t_mysort
		integer :: idx
	contains
		procedure :: swap_sort => myswap
	end type
contains
	subroutine myswap(self,a)
		class(t_mysort) :: self
		class(t_sort) :: a
		type(t_mysort), allocatable :: tmp
		select type(a)
		type is (t_mysort)
			call self%t_sort%swap_sort(a)
			allocate(tmp)

			tmp%idx=a%idx
			a%idx=self%idx
			self%idx=tmp%idx

			deallocate(tmp)
		end select
	end subroutine
	subroutine init(self,init_cfg)
		class(t_mc) :: self[ica1,*]
		logical :: init_cfg
		complex(8) :: H(Ns*spin,Ns*spin),cH(size(H,1),size(H,2)),D(size(H,1),size(H,2),sum(var(1:vn)%n)),Uk(size(H,1),size(H,2)),psi0_(Ns*Ns)
		real(8) :: E(size(H,1))
		integer :: l,i,j,n,n1,n2,nk(size(H,1),size(H,2)),nn_(0:Ns*Ns,2),kq(Ns),Ecfg_(Ns*spin,2)
		type(t_mysort) :: mysort(size(E))
		logical :: check=.true.
		
		if(init_cfg) then
			if(allocated(self%cfg)) deallocate(self%cfg)
			allocate(self%cfg(Ns*spin))
			do i=1,Ns
				if(i<=self%ne(1)) then
					self%cfg(i)=3
				elseif(i<=self%ne(2)) then
					self%cfg(i)=2
				else
					self%cfg(i)=0
				endif
			enddo
			self%cfg(Ns+1:)=0
			call fisher_yates_shuffle(self%cfg(:Ns))
			n=0
			do i=1,Ns
				if(self%cfg(i)==3) then
					n=n+2
					self%cfg(i)=n-1
					self%cfg(i+Ns)=n
				elseif(self%cfg(i)==2) then
					n=n+1
					self%cfg(i+Ns)=n
					self%cfg(i)=0
				endif
			enddo
		endif
		if(allocated(self%Ecfg)) deallocate(self%Ecfg)
		allocate(self%Ecfg(Ns*spin))
		if(allocated(self%Emf)) deallocate(self%Emf)
		allocate(self%Emf(Ns*spin))
		self%Ecfg=0

		call Hamilton(var,H)

		Uk=0d0
		do i=1,Ns
			do j=1,Ns
				Uk(j,i)=exp(-img*sum(latt%nb(0)%bd(j)%r*brizon%ok(i,:)))
			enddo
		enddo
		Uk(Ns+1:,Ns+1:)=Uk(:Ns,:Ns)
		H=matmul(transpose(conjg(Uk)),matmul(H,Uk))*(1d0/Ns)

		nk=0
		if(any(abs(var(:)%sg)==2)) then
			do i=1,size(brizon%k,1)
				call mat_diag(H(i::size(brizon%k,1),i::size(brizon%k,1)),E(i::size(brizon%k,1)))
				nk(i::size(brizon%k,1),i::size(brizon%k,1))=1
				!write(30,"(*(es12.4))")brizon%k(i,:2),E(i::size(brizon%k,1)),abs(H(i,i::size(brizon%k,1)))
			enddo
			self%Emf=E
			mysort%val=E
			mysort%idx=(/(i,i=1,size(E))/)
			call qsort(mysort)
			do i=1,sum(self%ne)
				self%Ecfg(mysort(i)%idx)=i
			enddo
			if(abs(mysort(sum(self%ne))%val-mysort(sum(self%ne)+1)%val)<1d-7) then
				if(this_image()==1) write(*,*)"warning, close-shell condition is not satisfied"
			endif
			if(self%sg==3) then
				nn_(0,:)=self%Ecfg(mysort(sum(self%ne)-1:sum(self%ne))%idx)
				self%Ecfg(mysort(sum(self%ne)-1:sum(self%ne))%idx)=0
			endif
		else
			do i=1,size(brizon%k,1)
				call mat_diag(H(i:Ns:size(brizon%k,1),i:Ns:size(brizon%k,1)),E(i:Ns:size(brizon%k,1)))
				call mat_diag(H(Ns+i::size(brizon%k,1),Ns+i::size(brizon%k,1)),E(Ns+i::size(brizon%k,1)))
				nk(i:Ns:size(brizon%k,1),i:Ns:size(brizon%k,1))=1
				nk(Ns+i::size(brizon%k,1),Ns+i::size(brizon%k,1))=1
			enddo
			self%Emf=E
			mysort%val=E
			mysort%idx=(/(i,i=1,size(E))/)
			call qsort(mysort(:Ns))
			call qsort(mysort(Ns+1:))
			do i=1,self%ne(1)
				self%Ecfg(mysort(i)%idx)=i
			enddo
			!stop
			do i=1,self%ne(2)
				self%Ecfg(mysort(i+Ns)%idx)=i+self%ne(1)
			enddo
			!write(*,*)mysort(self%ne(1))%val
			if(abs(mysort(self%ne(1))%val-mysort(self%ne(1)+1)%val)<1d-7.or.abs(mysort(Ns+self%ne(2))%val-mysort(Ns+self%ne(2)+1)%val)<1d-7) then
				if(self%sg/=3) then
					if(this_image()==1) write(*,"(A$)")"warning, close-shell condition is not satisfied, maybe set ne(1) to:"
					do i=self%ne(1),1,-1
						if(abs(mysort(self%ne(1))%val-mysort(i)%val)>1d-5) then
							if(this_image()==1) write(*,"(i5$)")i
							exit
						endif
					enddo
					do i=self%ne(1),Ns
						if(abs(mysort(self%ne(1))%val-mysort(i)%val)>1d-5) then
							if(this_image()==1) write(*,"(i5)")i-1
							exit
						endif
					enddo
				endif
			endif
			if(self%sg==3) then
				nn_(0,:)=self%Ecfg(mysort((/self%ne(1),self%ne(2)+Ns/))%idx)
				self%Ecfg(mysort((/self%ne(1),self%ne(2)+Ns/))%idx)=0
				self%Emf(:Ns)=self%Emf(:Ns)-self%Emf(mysort(self%ne(1)-1)%idx)
				self%Emf(Ns+1:)=self%Emf(Ns+1:)-self%Emf(mysort(self%ne(2)+Ns-1)%idx)
			endif
			do i=1,size(brizon%k,1)
				if(this_image()==1) write(71,"(es12.4$)")brizon%k(i,:)
				if(this_image()==1) write(71,"(2i3)")count(self%Ecfg(i:Ns:size(brizon%k,1))/=0),count(self%Ecfg(i+Ns::size(brizon%k,1))==0)
			enddo
		endif

		if(self%sg==3) then
			kq=0
			do i=1,Ns
				do j=1,Ns
					if(sum(abs(mod(brizon%ok(i,:2)/pi+brizon%ok(j,:2)/pi+(/2d0,2d0/)+1d-7,(/2d0,2d0/))-self%q(:2)/pi))<1d-5) then
						if(kq(i)/=0) then
							if(this_image(self,2)==1) write(*,*)"warn! kq"
						else
							kq(i)=j
						endif
					endif
				enddo
				if(kq(i)==0) then
					if(this_image(self,2)==1) write(*,*)self%q(:2)/pi,"warn! kq2"
					error stop
				endif
			enddo
			Ecfg_(:,1)=1-sign(1,-self%Ecfg)
			Ecfg_(:,2)=Ecfg_(:,1)
			psi0_=0d0
			l=0
			nn_(1:,:)=0
			do i=1,Ns
				do n2=1,size(E)
					if(Ecfg_(n2,1)==0.and.nk(kq(i)+Ns,n2)>0) then
						Ecfg_(n2,1)=1+sign(1,-Ecfg_(n2,1))
					o:	do n1=1,size(E)
							if(Ecfg_(n1,1)==0.and.nk(i,n1)>0) then
								Ecfg_(n1,1)=1+sign(1,-Ecfg_(n1,1))
								do j=1,l
									Ecfg_(abs(nn_(j,:)),2)=1-sign(1,nn_(j,:))
									if(all(Ecfg_(:,1)==Ecfg_(:,2))) then
										if(nn_(j,1)==-n2.and.nn_(j,2)==-n1) then
											psi0_(j)=psi0_(j)+conjg(H(i,n1))*conjg(H(kq(i)+Ns,n2))
										elseif(nn_(j,1)==-n1.and.nn_(j,2)==-n2) then
											psi0_(j)=psi0_(j)-conjg(H(i,n1))*conjg(H(kq(i)+Ns,n2))
										endif
										Ecfg_(abs(nn_(j,:)),2)=1-sign(1,-nn_(j,:))
										Ecfg_(n1,1)=1+sign(1,-Ecfg_(n1,1))
										cycle o
									endif
									Ecfg_(abs(nn_(j,:)),2)=1-sign(1,-nn_(j,:))
								enddo
								Ecfg_(n1,1)=1+sign(1,-Ecfg_(n1,1))
								!if(abs(conjg(H(i,n1))*conjg(H(kq(i)+Ns,n2)))>1d-6) then
									l=l+1
									nn_(l,:)=(/-n2,-n1/)
									psi0_(l)=psi0_(l)+conjg(H(i,n1))*conjg(H(kq(i)+Ns,n2))
								!endif
							endif
						enddo o
						Ecfg_(n2,1)=1+sign(1,-Ecfg_(n2,1))
					endif
				enddo
			enddo
			if(allocated(self%nn)) deallocate(self%nn,self%Ok2,self%Ek2,self%psi0)
			allocate(self%nn(0:l,2),self%Ok2(l,l),self%Ek2(l,l),self%psi0(l))
			self%nn=nn_(:l,:)
			self%nn(0,:)=self%nn(1,(/2,1/))
			self%Ecfg(abs(self%nn(0,:)))=nn_(0,:)
			self%psi0=psi0_(:l)
		endif
		H=matmul(Uk,H)*sqrt(1d0/Ns)
		call Hamilton(var,cH)
		if(check) then
			if(sum(abs(matmul(transpose(conjg(H)),matmul(cH,H))-diag(E)))>1d-6) error stop "hamilton err"
			check=.false.
		endif
		if(allocated(self%wf)) deallocate(self%wf)
		allocate(self%wf(Ns*spin,Ns*spin))
		self%wf=H


		if(self%sg==2) then
			if(allocated(self%dwf)) deallocate(self%dwf,self%g,self%S)
			allocate(self%dwf(Ns*spin,Ns*spin,sum(var(1:vn)%n)),self%g(sum(var(1:)%n)))
			allocate(self%S(size(self%g),size(self%g)))
			cH=transpose(conjg(H))
			self%dwf=0d0
			call dHamilton(var(1:vn),H,cH,D)

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
					Ecfg_(:,1)=self%Ecfg
					!$omp parallel do
					do i=1,size(E)
						do j=1,size(E)
							if(abs(E(i)-E(j))<1d-8.or.Ecfg_(i,1)/=0.or.Ecfg_(j,1)==0) then
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
			self%phy=0d0; self%S=0d0; self%g=0d0
			call random_number(seed,4285)
			call self%do_mc(seed+7*this_image())
			sync all
			if(this_image(self,2)==1) then
				do l=2,ica2
					self[this_image(self,1),1]%phy=self[this_image(self,1),1]%phy+self[this_image(self,1),l]%phy; self[this_image(self,1),1]%S=self[this_image(self,1),1]%S+self[this_image(self,1),l]%S; self[this_image(self,1),1]%g=self[this_image(self,1),1]%g+self[this_image(self,1),l]%g
				enddo
				self%phy=self%phy*1d0/ica2; self%S=self%S*1d0/ica2; self%g=self%g*1d0/ica2
				self%phy(ier)=sqrt(self%phy(ier)-self%phy(iE)**2)
			endif
			sync all
			if(this_image(self,2)/=1) then
				self%phy=self[this_image(self,1),1]%phy; self%S=self[this_image(self,1),1]%S; self%g=self[this_image(self,1),1]%g
			endif
			call heev(self%S,eg,'V')
			if(this_image()==1) write(*,"(es12.4$)")var(1:)%val(1),self%phy(iE),self%phy(ier)
			if(this_image()==1) write(*,"(i3$)")int(sign(1d0,self%g))
			if(this_image()==1) write(*,"(A$)")' | '
			if(this_image()==1) write(20,"(es12.4$)")var(1:)%val(1),self%phy(iE),self%phy(ier),self%g
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
			self%phy=0d0
			call random_number(seed,4285)
			call self%do_mc(seed+7*this_image())
			sync all
			if(this_image(self,2)==1) then
				do l=2,ica2
					self[this_image(self,1),1]%phy=self[this_image(self,1),1]%phy+self[this_image(self,1),l]%phy;
				enddo
				self%phy=self%phy*1d0/ica2
				self%phy(ier)=sqrt(self%phy(ier)-self%phy(iE)**2)
			endif
			sync all
			if(this_image(self,2)/=1) then
				self%phy=self[this_image(self,1),1]%phy
			endif
			sync all
			critical
			if(this_image(self,2)==1) write(*,"(*(es12.4))")self%q(:2)/pi,var(1:)%val(1),self%phy(iE),self%phy(ier),self%phy(iSq_pm)
			if(this_image(self,2)==1) write(10,"(*(es12.4))")var(1:)%val(1)
			endcritical
		case(3)
			self%Ok2=0d0; self%Ek2=0d0
			call random_number(seed,4285)
			call self%do_mc(seed+7*this_image())
			sync all
			if(this_image(self,2)==1) then
				do l=2,ica2
					self[this_image(self,1),1]%Ok2=self[this_image(self,1),1]%Ok2+self[this_image(self,1),l]%Ok2; self[this_image(self,1),1]%Ek2=self[this_image(self,1),1]%Ek2+self[this_image(self,1),l]%Ek2
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
		real(8) :: rd,isamp,Oq
		complex(8) :: pb
		integer :: i,j,l,l1,l2,info,n,apt,samp,sg,dly,np
		integer :: k(0:4),m(0:4),dcfg(0:4)
		logical :: is_update,is_full,is_accept
		type(randomNumberSequence) :: rnd
		real(8) :: ja(sum(var(vn+1:)%n))
		ja=put(var(vn+1:))
		call mt_init_random_seed(rnd,seed)
		phyp=0d0; phyl=0d0
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
			if(sum(abs(matmul(A(:,iEcfg),iA)-diag(1d0,size(A,1))))>1d-6) then
				call change(cfgl,icfg,dcfg,k,rnd)
				cfgl(abs(dcfg(2:dcfg(0):2)))=cfgl(abs(dcfg(2:dcfg(0):2)))+cfgl(dcfg(1:dcfg(0):2))
				cfgl(dcfg(1:dcfg(0):2))=cfgl(abs(dcfg(2:dcfg(0):2)))-cfgl(dcfg(1:dcfg(0):2))
				cfgl(abs(dcfg(2:dcfg(0):2)))=cfgl(abs(dcfg(2:dcfg(0):2)))-cfgl(dcfg(1:dcfg(0):2))
				icfg(k(1:k(0):2))=k(2:k(0):2)
				A=wf(icfg,:)
				iA=A(:,iEcfg)
				call mat_inv(iA,info)
				if(n>200) then
					stop "ini cfg err"
				endif
			else
				exit
			endif
		enddo
		WA=matmul(wf(:,iEcfg),iA)
		AW=matmul(iA,wf(icfg,:))
		WAW=matmul(wf(:,iEcfg),AW)
		WAW=WAW-wf
		n=0; apt=0; samp=0
		dly=0
		is_update=.true.
		is_full=.false.
	lm:	do
			n=n+1
			call change(cfgl,icfg,dcfg,k,rnd)

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
			pb=pb*exp(2d0*jast(cfgl,dcfg(1:dcfg(0)),ja))
			rd=getrandomreal(rnd)
			is_accept=.false.
			if(rd<real(pb)) then

				if(self%sg==3) then
					i=-1
					call get_pb(k(1:k(0)),shape(0),pb,WA=WA,WAl=WAl(:,:dly),WAr=WAr(:dly,:))
					if(abs(pb)<1d-10) then
						do i=1,ubound(nn,1)
							if(get_row(Ecfgl,(/-nn(0,:),nn(i,:)/),m,sg)) then
								call get_pb(k(1:k(0)),m(1:m(0)),pb,WA=WA,WAl=WAl(:,:dly),WAr=WAr(:dly,:),iA=iA,iAl=iAl(:,:dly),AW=AW,WAW=WAW,AWr=AWr(:dly,:))
								if(abs(pb)>1d0) then
									Ecfgl(iEcfg(m(1:m(0):2)))=0
									Ecfgl(m(2:m(0):2))=m(1:m(0):2)
									iEcfg(m(1:m(0):2))=m(2:m(0):2)
									nn(0,:)=nn(i,(/2,1/))
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

				rd=jast(cfgl,dcfg(1:dcfg(0)),jan=Ol(sum(var(1:vn)%n)+1:))

				cfgl(abs(dcfg(2:dcfg(0):2)))=cfgl(abs(dcfg(2:dcfg(0):2)))+cfgl(dcfg(1:dcfg(0):2))
				cfgl(dcfg(1:dcfg(0):2))=cfgl(abs(dcfg(2:dcfg(0):2)))-cfgl(dcfg(1:dcfg(0):2))
				cfgl(abs(dcfg(2:dcfg(0):2)))=cfgl(abs(dcfg(2:dcfg(0):2)))-cfgl(dcfg(1:dcfg(0):2))
				icfg(k(1:k(0):2))=k(2:k(0):2)

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
						phyl(iE)=get_energy(cfgl,WA,ja)*1d0/Ns
						!phyl%dsc=real(get_dsc(cfgl,WA,ja))
						phyl(iSq_pm)=get_spin_pm(cfgl,WA,self%q,ja)
						!phyl%Sq_zz=get_spin_zz(cfgl,q)
					case(2)
						phyl(iE)=get_energy(cfgl,WA,ja)*1d0/Ns
						call get_O(icfg,iEcfg,WA,iA,dwf,Ol(:size(dwf,3)))
					case(3)
						call get_overlap(cfgl,Ecfgl,nn,WA,iA,AW,WAW,Ok,Ek,ja)
						Oq=1d0/sum(real(Ok*conjg(Ok)))
					end select
				endif
				phyp=phyp+phyl
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
						!if(mod(samp,1024*8*8)==0) then
							!rewind(50)
							!write(50,*)size(brizon%k,1),brizon%k
							!write(50,*)samp,size(self%psi0),Ns
							!write(50,*)self%phy(iE),self%phy(iSq_pm),self%Emf,self%psi0,Ok2p*1d0/samp,Ek2p*1d0/samp,self%nn(1:,:)
						!endif
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

		phyp=phyp*isamp
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
			self%S=Sp; self%g=gp
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
			if(Smax(2)<abs(imag(Sq(l)))) Smax=(/omg(1)+domg*l,abs(imag(Sq(l)))/)
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
				write(ut+1,"(*(es17.9))")brizon%k(mod(abs(self%nn(l,1))-1,size(brizon%k,1))+1,:2),brizon%k(mod(abs(self%nn(l,2))-1,size(brizon%k,1))+1,:2),self%Emf(abs(self%nn(l,:))),pnn(l)
			enddo
		endif
		write(ut+1,"(x)")
	end subroutine
	subroutine do_var(self,n)
		class(t_mc) :: self[ica1,*]
		integer :: omp,n
		real(8) :: x(sum(var(1:)%n),n),dx,er,pgrad(size(x))
		type(t_mysort) :: El(n)
		integer :: i,hot
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
		x(:,1)=put(var(1:))
		init_cfg=.true.
		do 
			i=i+1
			if(allocated(self%g)) pgrad=self%g
			call get(var(1:),x(:,i))
			if(this_image()==1) write(*,"(i4$)")i
			call self%init(init_cfg)
			call self%do_vmc()
			El(i)%val=self%phy(iE)
			El(i)%idx=i
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
		call qsort(El)
		do i=2,min(n/5,30)
			x(:,El(1)%idx)=x(:,El(1)%idx)+x(:,El(i)%idx)
		enddo
		call get(var(1:),x(:,El(1)%idx)/min(n/5,30))
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
		!open(101,file="../data/lattice.dat")
		!open(111,file="../data/tmp.dat")
		!open(20,file=fn("../data/var.dat"))
		!open(40,file=fn("../data/phyvar.dat"))
		!open(50,file=fn("../data/matrix.dat"),access="stream")
		!open(70,file=fn("../data/spect.dat"))
		!open(71,file=fn("../data/spect_kmap.dat"))
		open(101,file="../data/lattice.dat")
		open(111,file="../data/tmp.dat")
		open(20,file="../data/var.dat")
		open(40,file="../data/phyvar.dat")
		open(50,file="../data/matrix_2.dat",access="stream")
		open(70,file="../data/spect_2.dat")
		open(71,file="../data/spect_kmap.dat")
		write(*,*)"coarray info: ",ica1,ica2,"out of",num_images()



		!write(*,*)"input the energy resolution and energy:"
		!read(*,*)et,E
		!read(50,*)i
		!rewind(50)
		!allocate(brizon%k(i,3))
		!read(50,*)i,brizon%k
		!do 
			!read(50,*)i,j,Ns
			!allocate(mc%psi0(j),mc%Ok2(j,j),mc%Ek2(j,j),mc%nn(j,2),mc%Emf(Ns*2))
			!read(50,*)mc%phy(iE),mc%phy(iSq_pm),mc%Emf,mc%psi0,mc%Ok2,mc%Ek2,mc%nn
			!write(*,*)mc%phy(iE),Ns,i,j
			!call mc%spect(70,et,(/-10d0,10d0/),5000,E)
			!mc%Ek2=transpose(conjg(mc%Ek2))
			!call mc%spect(70,et,(/-10d0,10d0/),5000)
			!deallocate(mc%psi0,mc%Ok2,mc%Ek2)
		!enddo
	endif
	sync all

	write(*,*)"runing ",this_image()," in ",hostname

	!stop

	call initial()
	if(this_image()==1) then
		write(50)size(brizon%k,1),brizon%k
	endif
	otime=0d0


	!var(1:)%val(1)=(/0d0,0.36d0,0.44d0/)
	!var(1:)%val(1)=(/0d0,0.36d0,0d0/)
	!var(1:)%val(1)=(/-0.39d0,0.32d0,-0.205d0/)
	!var(1:)%val(1)=(/-0.46d0,0.300d0,-0.268d0/)
	!var(1:)%val(1)=(/-1.39752E+00,8.88530E-01*2E0,1.32678E-02,-2.00021E-01/)
	!var(1:)%val(1)=(/-1.4051E+00,2.4353E+00,-1.9875E-01,-2.6827E-01/)
	!var(1:)%val(1)=(/-1.2070E+00,1.8269E+00,-1.3504E-01,-1.5381E-01/)
	!var(1:)%val(1)=(/-1.53531E+00,8.07142E-01*2E0,-6.06570E-02,-2.62224E-01,-2.15494E-01/)
	!var(1:)%val(1)=(/-8.0718E-01,2.0919E-01,-9.9164E-02/) ! 16x16
	!var(1:)%val(1)=(/-6.6706E-01,2.1394E-01,-1.3299E-01/) ! 8x8
	var(1:)%val(1)=(/0d0,3.9616d-1,0d0/) ! 8x8
	!var(1:)%val(1)=(/1.9616d-1,0d0/) ! 8x8
	!var(1:)%val(1)=(/-8.6811E-01,2.0951E-01,-9.7832E-02/) ! 12x12
	!call variation()
	!stop
	!call export_data(101)
	call omp_set_nested(.false.)
	!call omp_set_max_active_levels(1)
	call omp_set_dynamic(.false.)
	call mkl_set_dynamic(0)
	call omp_set_schedule(omp_sched_static,0)

	call mkl_set_num_threads(1)
	call omp_set_num_threads(mkl_get_max_threads())


	mc%hot=1024*8
	mc%step=Ns
	mc%delay=2
	mc%ne(1)=Ns/2-16
	mc%ne(2)=Ns-mc%ne(1)
	!var(1:)%val(1)=(/-3.2296E-01,2.3218E-01,7.3575E-02,3.4608E-02/)
	!var(1:)%val(1)=(/-8.1623E-01, 2.5240E-01,-1.2011E-01, 4.2456E-01, 2.0609E-01, 4.0513E-02/)
	var(1:)%val(1)=(/-1.4051E+00,2.4353E+00,-1.9875E-01,-2.6827E-01/)
	!var(1:)%val(1)=(/2.0628E-01,0.0000E-00,1.9582E-01,1.1811E-01,2.6489E-02/)
	!var(1:)%val(1)=(/-8.6811E-01,2.0951E-01,-9.7832E-02/) ! 1/8
	!var(1:)%val(1)=(/-7.5957E-01,2.2205E-01,-8.2460E-02/) ! 16x16 -14
	!var(1:)%val(1)=(/-9.5867E-01,1.6445E-01,-1.0910E-01/) ! 16x16 -22
	!var(1:)%val(1)=(/-7.8290E-01,2.0821E-01,-7.6381E-02/)
	!var(1:)%val(1)=(/2.1445d-01/)!0.125,ddw,compare
	!var(1:)%val(1)=(/2.0268d-01/)!0.125,ddw,compare
	!var(1:)%val(1)=(/1.9268d-01/)!0.125,ddw,compare
	!var(1:)%val(1)=(/0.E+00,2.8386E-01,3.0E-01/)
	!var(1:)%val(1)=(/2.8386E-01/)
	!var(vn+1:)%val(1)=0d0
	!var(1:)%val(1)=(/6.1000d-01/)
	!do i=-10,10
		!mc%sg=2
		!var(4)%val(1)=i*0.03d0
		!call mc%init(.true.)
		!call mc%do_vmc(1)
		!write(*,"(x)")
	!enddo
	!var(1:)%val(1)=(/1.6047E-01/)
	!do i=1,100
		!var(1)%val=1d-1+i/200d0
		!write(*,*)var(1)%val
		!call mc%init(.true.)
	!enddo
	!otime(0)=omp_get_wtime()
	!call mc%do_var(100)
	!otime(1)=otime(1)+omp_get_wtime()-otime(0)
	!call show_time()
	!stop
	mc%q=q(:,this_image(mc,1))
	mc%num=this_image(mc,1)

	mc%sg=1
	mc%samp=1024*8*8
	mc%hot=1024*8*8
	mc%step=Ns

	call mc%init(.true.)
	call mc%do_vmc()

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
	mc%samp=1024*8*8*32*8*4
	mc%step=nint(sqrt(real(Ns)))
	call mc%init(.true.)
	call mc%do_vmc()

	!if(this_image()/=1.and.this_image()<=ica1) sync images(this_image()-1)
	!if(this_image()<=ica1) then
		!if(this_image()==1) then
			!open(50,file="../data/matrix.dat")
			!open(70,file="../data/spect.dat")
		!else
			!open(50,file="../data/matrix.dat",position="append")
			!open(70,file="../data/spect.dat",position="append")
		!endif
		!write(*,*)this_image(mc)
		!write(50,*)-1,size(mc%psi0),Ns
		!write(50,*)mc%phy(iE),mc%phy(iSq_pm),mc%Emf,mc%psi0,mc%Ok2,mc%Ek2,mc%nn(1:,:)
		!call mc%spect(70,0.02d0,(/-10d0,10d0/),5000)
		!mc%Ek2=transpose(conjg(mc%Ek2))
		!call mc%spect(70,0.02d0,(/-10d0,10d0/),5000)
		!close(50)
		!close(70)
	!endif
	!if(this_image()<ica1) sync images(this_image()+1)

	if(this_image()==1) then
		write(*,*)"finished, exporting data....",mc%samp
		do l=1,ica1
			write(50)-1,size(mc[l,1]%psi0),Ns,mc[l,1]%q,mc[l,1]%phy(iE),mc[l,1]%phy(iSq_pm),mc[l,1]%Emf,mc[l,1]%psi0,mc[l,1]%Ok2,mc[l,1]%Ek2,mc[l,1]%nn(1:,:)
			!call mc%spect(70,0.02d0,(/-10d0,10d0/),5000)
			!mc[l,1]%Ek2=transpose(conjg(mc[l,1]%Ek2))
			!call mc%spect(70,0.02d0,(/-10d0,10d0/),5000)
		enddo
		rewind(50)
		read(50)i
		deallocate(brizon%k)
		allocate(brizon%k(i,3))
		read(50)brizon%k
		do l=1,ica1
			read(50)i,j,Ns
			if(allocated(mc%psi0)) deallocate(mc%psi0,mc%Ok2,mc%Ek2,mc%Emf,mc%nn)
			allocate(mc%psi0(j),mc%Ok2(j,j),mc%Ek2(j,j),mc%nn(j,2),mc%Emf(Ns*2))
			read(50)mc%q,mc%phy(iE),mc%phy(iSq_pm),mc%Emf,mc%psi0,mc%Ok2,mc%Ek2,mc%nn
			write(*,*)mc%q,mc%phy(iE),Ns,i,j
			call mc%spect(70,0.02d0,(/-10d0,10d0/),5000)
			mc%Ek2=transpose(conjg(mc%Ek2))
			call mc%spect(70,0.02d0,(/-10d0,10d0/),5000)
			deallocate(mc%psi0,mc%Ok2,mc%Ek2,mc%Emf,mc%nn)
		enddo
	endif
end program
