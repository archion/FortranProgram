module global
	use M_hamilton_test
	use ifport, ifort_qsort => qsort
	use M_omp_rand
	implicit none
	real(8) :: t(2)=(/1d0,-0.3d0/)
	real(8), parameter :: DJ=0.3d0,V=0d0,U=0d0
	integer :: n_omp=4
	integer, parameter :: opt=2 ! 1: Tr(AdA)
								! 2: O=cicj
	integer, parameter :: iE=1,ier=2,idsc=3,iaf=4,iddw=5,iSq_pm=6,iSq_zz=7,iCq=8
	integer :: Ns,vn,Cnb=1
contains
	subroutine initial()
		integer :: i,l
		real(8) :: q(3)
		q=(/1d0/8d0,0d0,0d0/)*2d0*pi
		allocate(var(-1000:1000))
		call init_random_seed()
		! lattice
		latt%a1=(/1d0,0d0,0d0/)
		latt%a2=(/0d0,1d0,0d0/)
		!latt%c1=(/8d0,0d0,0d0/)
		!latt%c2=(/0d0,2d0,0d0/)
		!latt%c1=(/1d0,1d0,0d0/)
		!latt%c2=(/-1d0,1d0,0d0/)
		latt%c1=latt%a1
		latt%c2=latt%a2
		latt%T1=(/1d0,0d0,0d0/)*8
		latt%T2=(/0d0,1d0,0d0/)*8
		latt%bdc=(/1d0,1d0,0d0/)
		allocate(latt%rsb(1,3))
		latt%rsb(1,:)=(/0d0,0d0,0d0/)
		latt%n1=1
		latt%n2=1
		call latt%gen_latt(size(t))
		call latt%gen_brizon(brizon)
		!call check_lattice(101)
		write(*,*)"Total site number is: ",latt%Ns

		! cp
		call gen_var(sg=1,nb=0)
		var(iv(0))%bd=-1d0
		var(iv(0))%val=-1.39752d0

		! dsc
		call gen_var(sg=2,nb=1)
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
			!var(iv(0))%bd(i)=ab(latt%sb(1,1)%nb(var(iv(0))%nb)%bd(i)%i(1))*dwave(i)*img
		!enddo
		!var(iv(0))%val=1d-3

		!! sdw
		!call gen_var(sg=4,nb=0)
		!do i=1,size(var(iv(0))%bd)
			!var(iv(0))%bd(i)=ab(i)!*sin(sum(q*(latt%sb(1,1)%nb(var(iv(0))%nb)%bd(i)%r)))
		!enddo
		!var(iv(0))%val=8.88530d-1*2d0

		!! on site cdw
		!call gen_var(sg=3,nb=0)
		!do i=1,size(var(iv(0))%bd)
			!var(iv(0))%bd(i)=cos(sum(2d0*q*(latt%sb(1,1)%nb(var(iv(0))%nb)%bd(i)%r)))
		!enddo
		!var(iv(0))%val=1.32678d-2

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

		call var_shrink()

		Ns=latt%Ns
	end subroutine
	real(8) function ab(i)
		integer :: i
		ab=(-1d0)**(mod(nint(sum(latt%sb(1,1)%nb(0)%bd(i)%r(:2))),2))
	end function
	real(8) function charge_symmetry(i)
		integer :: i
		!charge_symmetry=1d0
		if(abs(latt%sb(1,1)%nb(Cnb)%bd(i)%dr(2))<1d-6) then
			charge_symmetry=1d0
		else
			charge_symmetry=-1d0
		endif
	end function
end module

module mc_utility
	use global
	use blas95, only : gemm, gemv, her, gerc
	implicit none
contains
	subroutine get_pb(k,m,pb,WA,iA,AW,WAW,wf)
		integer :: k(:)
		integer :: m(:)
		complex(8) :: WA(:,:),pb
		complex(8), optional :: AW(:,:),WAW(:,:),iA(:,:),wf(:,:)
		integer :: i,j,n
		integer :: k_(size(k)/2+size(m)/2),m_(size(k_))
		complex(8) :: A(size(k_,1),size(k_,1))
		n=size(k)/2
		k_=(/k(1::2),m(1::2)/)
		m_=(/k(2::2),m(2::2)/)
		do i=1,size(k_)
			do j=1,size(k_)
				if(i<=n.and.j<=n) then
					A(i,j)=WA(m_(i),k_(j))
				elseif(i<=n.and.j>n) then
					A(i,j)=WAW(m_(i),m_(j))-wf(m_(i),m_(j))
				elseif(i>n.and.j<=n) then
					A(i,j)=iA(k_(i),k_(j))
				else
					A(i,j)=AW(k_(i),m_(j))
				endif
			enddo
		enddo
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
	subroutine update(k,WA,iA,AW,WAW,W2,wf)
		complex(8) :: WA(:,:)
		complex(8), optional :: AW(:,:),WAW(:,:),iA(:,:),W2(:,:),wf(:,:)
		integer :: k(:)
		integer :: k_(size(k)/2),m_(size(k_))
		complex(8) :: WAr(size(k_),size(WA,2)),WAl(size(WA,1),size(k_)),iAl(size(WA,2),size(k_)),WAWr(size(k_),size(WA,1)),ipb(size(k_),size(m_)),dW2(size(k_),size(WA,1)),p1,m1
		integer :: i,j
		p1=1d0
		m1=-1d0
		k_=k(1::2)
		m_=k(2::2)
		do i=1,size(k_)
			do j=1,size(k_)
				ipb(i,j)=WA(m_(i),k_(j))
			enddo
		enddo
		select case(size(k_))
		case(0)
			return
		case(1)
			ipb=1d0/ipb
		case(2)
			ipb=reshape((/ipb(2,2),-ipb(2,1),-ipb(1,2),ipb(1,1)/),(/2,2/))/(ipb(1,1)*ipb(2,2)-ipb(2,1)*ipb(1,2))
		case default
			write(*,*)"err"
			stop
		end select
		WAr=WA(m_,:)
		WAr(1:size(k_),k_)=WAr(1:size(k_),k_)-diag(1d0,size(k_))
		WAr=matmul(ipb,WAr)
		WAl=WA(:,k_)
		!WA=WA-matmul(WAl,WAr)
		call gemm(WAl,WAr,WA,alpha=m1,beta=p1)
		if(present(iA)) then
			iAl=iA(:,k_)
			!iA=iA-matmul(iAl,WAr)
			call gemm(iAl,WAr,iA,alpha=m1,beta=p1)
		endif
		if(present(WAW)) then
			dW2=wf(m_,:)-W2(k_,:)
			WAWr=WAW(m_,:)-W2(k_,:)+matmul(WA(m_,k_)-diag(1d0,size(k_)),dW2)
			WAWr=matmul(ipb,WAWr)
			!AW=AW+matmul(iA(:,k_),dW2)-matmul(iAl,WAWr)
			call gemm(iA(:,k_),dW2,AW,beta=p1)
			call gemm(iAl,WAWr,AW,alpha=m1,beta=p1)
			!WAW=WAW+matmul(WA(:,k_),dW2)-matmul(WAl,WAWr)
			call gemm(WA(:,k_),dW2,WAW,beta=p1)
			call gemm(WAl,WAWr,WAW,alpha=m1,beta=p1)
			W2(k_,:)=wf(m_,:)
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
		integer :: cfg(:),dcfg(:),sg
		integer :: k(0:)
		integer, optional :: P(:)
		integer :: i,j,n,cfg_(size(cfg)),dcfg_(size(dcfg))
		integer :: db(0:10)
		k=0
		get_row=.true.
		db(0)=0
		do i=1,size(dcfg)
			if(all(db(1:db(0))/=(mod(abs(dcfg(i))-1,Ns)+1))) then
				db(0)=db(0)+1
				db(db(0))=mod(abs(dcfg(i))-1,Ns)+1
				cfg_(db(db(0)))=cfg(db(db(0)))
				cfg_(db(db(0))+Ns)=cfg(db(db(0))+Ns)
			endif
		enddo
		dcfg_=0
		sg=1
	o:	do i=1,size(dcfg)+1
			if(present(P)) then
				if(any(P==i).or.i==size(dcfg)+1) then
					do j=1,db(0)
						if(cfg_(db(j))/=0.and.cfg_(db(j)+Ns)==0) then
							get_row=.false.
							return
						endif
					enddo
				endif
			endif
			if(i==size(dcfg)+1) exit
			if(dcfg(i)>0.and.cfg_(abs(dcfg(i)))/=0) then
				cfg_(abs(dcfg(i)))=0
				do j=size(dcfg_),1,-2
					if(dcfg_(j)==-dcfg(i)) then
						dcfg_(j)=0
						sg=sg*(-1)**count(dcfg_(j+1:)/=0)
						do n=j+2,size(dcfg_),2
							if(dcfg_(n)==0) exit
							if(dcfg_(n-1)/=0) sg=-sg
							dcfg_(n-2)=dcfg_(n)
							dcfg_(n)=0
						enddo
						cycle o
					endif
					if(dcfg_(j-1)==0) then
						n=j-1
					endif
				enddo
				sg=sg*(-1)**count(dcfg_(n+1:)/=0)
				dcfg_(n)=dcfg(i)
			elseif(dcfg(i)<0.and.cfg_(abs(dcfg(i)))==0) then
				cfg_(abs(dcfg(i)))=1
				do j=size(dcfg_)-1,1,-2
					if(dcfg_(j)==-dcfg(i)) then
						dcfg_(j)=0
						sg=sg*(-1)**count(dcfg_(j+1:)/=0)
						do n=j+2,size(dcfg_),2
							if(dcfg_(n)==0) exit
							if(dcfg_(n-1)/=0) sg=-sg
							dcfg_(n-2)=dcfg_(n)
							dcfg_(n)=0
						enddo
						cycle o
					endif
					if(dcfg_(j+1)==0) then
						n=j+1
					endif
				enddo
				sg=sg*(-1)**count(dcfg_(n+1:)/=0)
				dcfg_(n)=dcfg(i)
			else
				get_row=.false.
				return
			endif
		enddo o
		n=0
		do i=1,size(dcfg_)
			if(dcfg_(i)>0) then
				n=n+1
				k(n)=abs(cfg(dcfg_(i)))
			elseif(dcfg_(i)<0) then
				n=n+1
				k(n)=abs(dcfg_(i))
			endif
		enddo
		k(0)=n
		!k(1)=raiseqq(sig$int)
	end function
	complex(8) function get_charge(cfg,D,q)
		integer :: cfg(:)
		complex(8) :: D(:,:)
		real(8) :: q(3)
		complex(8) :: pb,sm
		integer :: n1,n2,i1,j1,i2,j2,p1,p2,k(0:4),sg
		get_charge=0d0
		!$omp parallel do collapse(2) reduction(+:get_charge) private(i1,i2,j1,j2,pb,sm,sg,k) if(.not.omp_in_parallel())
		do n1=1,size(latt%sb(1,1)%nb(Cnb)%bd)
			do n2=1,size(latt%sb(1,1)%nb(Cnb)%bd)
				i1=latt%sb(1,1)%nb(Cnb)%bd(n1)%i(1)
				j1=latt%sb(1,1)%nb(Cnb)%bd(n1)%i(2)
				i2=latt%sb(1,1)%nb(Cnb)%bd(n2)%i(1)
				j2=latt%sb(1,1)%nb(Cnb)%bd(n2)%i(2)
				sm=exp(img*sum(q*(latt%sb(1,1)%nb(Cnb)%bd(n1)%r-latt%sb(1,1)%nb(Cnb)%bd(n2)%r)))*charge_symmetry(n1)*charge_symmetry(n2)
				do p1=1,2
					do p2=1,2
						if(get_row(cfg,(/j2,-i2,i1,-j1/),k,sg,(/3/))) then
							call get_pb(k(1:k(0)),shape(0),pb,D)
							get_charge=get_charge+pb*sg*sm
						endif
						call swap(i2,j2)
					enddo
					call swap(i1,j1)
				enddo
			enddo
		enddo
		!$omp end parallel do
		get_charge=get_charge
	end function
	complex(8) function get_spin_pm(cfg,D,q)
		complex(8) :: D(:,:)
		integer :: cfg(:)
		real(8) :: q(3)
		complex(8) :: pb
		integer :: i,j,sg
		integer :: k(0:4)
		get_spin_pm=0d0
		!$omp parallel do collapse(2) reduction(+:get_spin_pm) private(pb,k,sg)
		do i=1,Ns
			do j=1,Ns
				if(get_row(cfg,(/i,i+Ns,-j-Ns,-j/),k,sg,shape(0))) then
					call get_pb(k(1:k(0)),shape(0),pb,D)
					get_spin_pm=get_spin_pm+pb*exp(img*sum(q*(latt%sb(1,1)%nb(0)%bd(i)%r-latt%sb(1,1)%nb(0)%bd(j)%r)))*sg
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
				get_spin_zz=get_spin_zz+0.25d0*((1-sign(1,-cfg(i)))/2-(1-sign(1,-cfg(i+Ns)))/2)*((1-sign(1,-cfg(j)))/2-(1-sign(1,-cfg(j+Ns)))/2)
			enddo
		enddo
		!$omp end parallel do
		get_spin_zz=get_spin_zz/Ns
	end function
	complex(8) function get_dsc(cfg,D)
		complex(8) :: D(:,:)
		integer :: cfg(:)
		complex(8) :: pb
		integer :: i1,j1,i2,j2,sg,l,n1,n2,p1,p2
		integer :: k(0:4)
		get_dsc=0d0
		!$omp parallel do collapse(2) reduction(+:get_dsc) private(pb,k,sg,i1,j1,i2,j2)
		do n1=1,size(latt%sb(1,1)%nb(1)%bd)
			do n2=1,size(latt%sb(1,1)%nb(1)%bd)
				i1=latt%sb(1,1)%nb(1)%bd(n1)%i(1)
				j1=latt%sb(1,1)%nb(1)%bd(n1)%i(2)
				i2=latt%sb(1,1)%nb(1)%bd(n2)%i(1)
				j2=latt%sb(1,1)%nb(1)%bd(n2)%i(2)
				if(any(latt%sb(1,1)%nb(1)%bd(n1)%i==i2).or.any(latt%sb(1,1)%nb(1)%bd(n1)%i==j2)) cycle
				do p1=1,2
					do p2=1,2
						if(get_row(cfg,(/-j1-Ns,i1,-i2,j2+Ns/),k,sg,shape(0))) then
							call get_pb(k(1:k(0)),shape(0),pb,D)
							get_dsc=get_dsc+pb*sg*dwave(n1)*dwave(n2)*latt%sb(1,1)%nb(1)%bd(n1)%bdc*latt%sb(1,1)%nb(1)%bd(n2)%bdc
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
	complex(8) function get_energy(cfg,D)
		complex(8) :: D(:,:)
		integer :: cfg(:)
		complex(8) :: pb
		integer :: i,j,n,l,c(2,2),p,sg
		integer :: k(0:4)
		get_energy=0d0
		!$omp parallel do collapse(2) reduction(+:get_energy) private(i,j,c,k,pb,sg)
		do l=1,size(t)
			do n=1,size(latt%sb(1,1)%nb(l)%bd)
				i=latt%sb(1,1)%nb(l)%bd(n)%i(1)
				j=latt%sb(1,1)%nb(l)%bd(n)%i(2)
				do p=1,2
					if(get_row(cfg,(/i,-j/),k,sg,shape(0))) then
						call get_pb(k(1:k(0)),shape(0),pb,D)
						get_energy=get_energy-sg*t(l)*pb*latt%sb(1,1)%nb(l)%bd(n)%bdc!*-1
					endif
					if(get_row(cfg,(/-i-Ns,j+Ns/),k,sg,shape(0))) then
						call get_pb(k(1:k(0)),shape(0),pb,D)
						get_energy=get_energy-sg*t(l)*pb*conjg(latt%sb(1,1)%nb(l)%bd(n)%bdc)!*-1
					endif
					if(l==1) then
						if(get_row(cfg,(/i,i+Ns,-j-Ns,-j/),k,sg,shape(0))) then
							call get_pb(k(1:k(0)),shape(0),pb,D)
							get_energy=get_energy+0.5d0*DJ*pb*sg!*-2
						endif
					endif
					j=latt%sb(1,1)%nb(l)%bd(n)%i(1)
					i=latt%sb(1,1)%nb(l)%bd(n)%i(2)
				enddo
				if(l==1) then
					c(:,1)=abs((1-sign(1,-cfg(i::Ns)))/2-(/0,1/))
					c(:,2)=abs((1-sign(1,-cfg(j::Ns)))/2-(/0,1/))
					get_energy=get_energy+V*sum(c(:,1))*sum(c(:,2))+0.25d0*DJ*(c(1,1)-c(2,1))*(c(1,2)-c(2,2))
				endif
			enddo
		enddo
		!$omp end parallel do
		get_energy=get_energy/Ns
	end function
	subroutine get_O(icfg,iEcfg,WA,iA,dwf,O)
		complex(8) :: WA(:,:),iA(:,:),dwf(:,:,:),O(:)
		complex(8) :: beta
		integer :: icfg(:),iEcfg(:)
		integer :: l,i
		O=0d0
		beta=1d0
		select case(opt)
		case(1)
			!$omp parallel do reduction(+:O)
			do i=1,size(iA,1)
				!O=O+matmul(iA(i,:),dwf(icfg,iEcfg(i),:))
				call gemv(dwf(icfg,iEcfg(i),:),iA(i,:),O,trans="T",beta=beta)
			enddo
			!$omp end parallel do
		case(2)
			!$omp parallel do reduction(+:O)
			do i=1,size(WA,1)
				!O=O+matmul(WA(i,:),dwf(icfg,i,:))
				call gemv(dwf(icfg,i,:),WA(i,:),O,trans="T",beta=beta)
			enddo
			!$omp end parallel do
		end select
	end subroutine
	subroutine get_Spb(cfg,dcfg,Ecfg,nn,q,pb,WA,iA,AW,WAW,wf)
		integer :: cfg(:),dcfg(0:),Ecfg(:),nn(:,:)
		real(8) :: q(:)
		complex(8) :: WA(:,:),AW(:,:),WAW(:,:),iA(:,:),wf(:,:)
		complex(8) :: pb,sm
		complex(8) :: pb2(2),pbs(2)
		integer :: i,l,n,j,p,sg
		integer :: k(0:4),m(0:4)
		if(.not.get_row(cfg,dcfg(1:dcfg(0)),k,sg,shape(0))) stop "err get_Spb"
		pbs=0d0

		!$omp parallel do private(pb2,sg,m) reduction(+:pbs)
		do i=1,ubound(nn,1)
			if(get_row(Ecfg,nn(i,:),m,sg)) then
				call get_pb(shape(0),m(1:m(0)),pb2(1),WA,iA,AW,WAW,wf)
				pbs(1)=pbs(1)+pb2(1)*conjg(pb2(1))
				call get_pb(k(1:k(0)),m(1:m(0)),pb2(2),WA,iA,AW,WAW,wf)
				pbs(2)=pbs(2)+pb2(2)*conjg(pb2(2))
			endif
		enddo
		!$omp end parallel do
		pb=pbs(2)/pbs(1)

		!!$omp parallel do private(pb2,sg,k,m,i,j,sm) reduction(+:pbs)
		!do n=1,size(latt%sb(1,1)%nb(Cnb)%bd)
			!i=latt%sb(1,1)%nb(Cnb)%bd(n)%i(1)
			!j=latt%sb(1,1)%nb(Cnb)%bd(n)%i(2)
			!sm=exp(img*sum(q*latt%sb(1,1)%nb(Cnb)%bd(n)%r))*charge_symmetry(n)
			!do p=1,2
				!if(get_row(cfg,(/i,-j/),k,sg,shape(0))) then
					!call get_pb(k(1:k(0)),shape(0),pb2(1),WA,iA,AW,WAW,wf)
					!pbs(1)=pbs(1)+pb2(1)*sm*sg
				!endif
				!if(get_row(cfg,(/dcfg(1:dcfg(0)),i,-j/),k,sg,shape(0))) then
					!call get_pb(k(1:k(0)),shape(0),pb2(2),WA,iA,AW,WAW,wf)
					!pbs(2)=pbs(2)+pb2(2)*sm*sg
				!endif
				!j=latt%sb(1,1)%nb(Cnb)%bd(n)%i(1)
				!i=latt%sb(1,1)%nb(Cnb)%bd(n)%i(2)
			!enddo
		!enddo
		!!$omp end parallel do
		!pb=pbs(2)/pbs(1)
		!pb=conjg(pb)*pb
	end subroutine
	subroutine get_overlap(cfg,Ecfg,nn,q,WA,iA,AW,WAW,wf,O,E)
		integer :: cfg(:),Ecfg(:),nn(:,:)
		real(8) :: q(:)
		complex(8) :: WA(:,:),AW(:,:),WAW(:,:),iA(:,:),wf(:,:),E(:,:),O(:,:)
		complex(8) :: Ek(size(E,1)),Ok(size(O,1)),pb,Oq,tmp,O_,E_,sm
		real(8) :: iNs
		integer :: i1,j1,i,j,l,n,n1,p1,c(2,2),p,sg,sgn
		integer :: k(0:6),m(0:4)
		iNs=1d0/Ns
		if(size(nn,1)+1==size(Ok)) then
			Ok=0d0
			Ek=0d0
			E_=0d0
			O_=0d0
			!$omp parallel do private(i,j,i1,j1,c,pb,sg,sm,k,m,tmp) reduction(+:O_,E_)
			do n1=1,size(latt%sb(1,1)%nb(Cnb)%bd)
				i1=latt%sb(1,1)%nb(Cnb)%bd(n1)%i(1)
				j1=latt%sb(1,1)%nb(Cnb)%bd(n1)%i(2)
				sm=exp(img*sum(q*latt%sb(1,1)%nb(Cnb)%bd(n1)%r))*charge_symmetry(n1)
				do p1=1,2
					tmp=0d0
					if(get_row(cfg,(/i1,-j1/),k,sg,shape(0))) then
						call get_pb(k(1:k(0)),shape(0),tmp,WA,iA,AW,WAW,wf)
						tmp=tmp*sg*sm
						O_=O_+tmp
					endif
					do l=1,size(t)
						do n=1,size(latt%sb(1,1)%nb(l)%bd)
							i=latt%sb(1,1)%nb(l)%bd(n)%i(1)
							j=latt%sb(1,1)%nb(l)%bd(n)%i(2)
							do p=1,2
								if(get_row(cfg,(/i,-j,i1,-j1/),k,sg,(/3/))) then
									call get_pb(k(1:k(0)),shape(0),pb,WA,iA,AW,WAW,wf)
									E_=E_-pb*t(l)*latt%sb(1,1)%nb(l)%bd(n)%bdc*sg*sm
								endif
								if(get_row(cfg,(/-i-Ns,j+Ns,i1,-j1/),k,sg,(/3/))) then
									call get_pb(k(1:k(0)),shape(0),pb,WA,iA,AW,WAW,wf)
									E_=E_-pb*t(l)*conjg(latt%sb(1,1)%nb(l)%bd(n)%bdc)*sg*sm
								endif
								if(l==1) then
									if(get_row(cfg,(/i,i+Ns,-j-Ns,-j,i1,-j1/),k,sg,shape(0))) then
										call get_pb(k(1:k(0)),shape(0),pb,WA,iA,AW,WAW,wf)
										E_=E_+pb*0.5d0*DJ*sg*sm
									endif
								endif
								j=latt%sb(1,1)%nb(l)%bd(n)%i(1)
								i=latt%sb(1,1)%nb(l)%bd(n)%i(2)
							enddo
							if(l==1) then
								c(:,1)=abs((1-sign(1,-cfg(i::Ns)))/2-(/0,1/))
								c(:,2)=abs((1-sign(1,-cfg(j::Ns)))/2-(/0,1/))
								E_=E_+tmp*(V*sum(c(:,1))*sum(c(:,2))+0.25d0*DJ*(c(1,1)-c(2,1))*(c(1,2)-c(2,2)))
							endif
						enddo
					enddo
					j1=latt%sb(1,1)%nb(Cnb)%bd(n1)%i(1)
					i1=latt%sb(1,1)%nb(Cnb)%bd(n1)%i(2)
				enddo
			enddo
			!$omp end parallel do
			Ek(ubound(nn,1)+1)=E_
			Ok(ubound(nn,1)+1)=O_
			!write(*,*)E_,O_
		elseif(size(nn,1)+Ns==size(Ok)) then
			Ok=0d0
			Ek=0d0
			!$omp parallel do private(i,j,i1,j1,c,pb,sg,k,O_,E_,sm) reduction(+:Ek,Ok)
			do i1=1,Ns
				do j1=1,Ns
					E_=0d0
					O_=0d0
					if(get_row(cfg,(/i1,-j1/),k,sg,shape(0))) then
						call get_pb(k(1:k(0)),shape(0),pb,WA,iA,AW,WAW,wf)
						O_=O_+pb*sg
					endif
					do l=1,size(t)
						do n=1,size(latt%sb(1,1)%nb(l)%bd)
							i=latt%sb(1,1)%nb(l)%bd(n)%i(1)
							j=latt%sb(1,1)%nb(l)%bd(n)%i(2)
							do p=1,2
								if(get_row(cfg,(/i,-j,i1,-j1/),k,sg,(/3/))) then
									call get_pb(k(1:k(0)),shape(0),pb,WA,iA,AW,WAW,wf)
									E_=E_-pb*t(l)*latt%sb(1,1)%nb(l)%bd(n)%bdc*sg
								endif
								if(get_row(cfg,(/-i-Ns,j+Ns,i1,-j1/),k,sg,(/3/))) then
									call get_pb(k(1:k(0)),shape(0),pb,WA,iA,AW,WAW,wf)
									E_=E_-pb*t(l)*conjg(latt%sb(1,1)%nb(l)%bd(n)%bdc)*sg
								endif
								if(l==1) then
									if(get_row(cfg,(/i,i+Ns,-j-Ns,-j,i1,-j1/),k,sg,shape(0))) then
										call get_pb(k(1:k(0)),shape(0),pb,WA,iA,AW,WAW,wf)
										E_=E_+pb*0.5d0*DJ*sg
									endif
								endif
								j=latt%sb(1,1)%nb(l)%bd(n)%i(1)
								i=latt%sb(1,1)%nb(l)%bd(n)%i(2)
							enddo
							if(l==1) then
								c(:,1)=abs((1-sign(1,-cfg(i::Ns)))/2-(/0,1/))
								c(:,2)=abs((1-sign(1,-cfg(j::Ns)))/2-(/0,1/))
								E_=E_+O_*(V*sum(c(:,1))*sum(c(:,2))+0.25d0*DJ*(c(1,1)-c(2,1))*(c(1,2)-c(2,2)))
							endif
						enddo
					enddo
					do n=1,Ns
						sm=exp(img*sum(brizon%ok(n,:)*(-latt%sb(1,1)%nb(0)%bd(i1)%r+latt%sb(1,1)%nb(0)%bd(j1)%r)+q*latt%sb(1,1)%nb(0)%bd(j1)%r))*iNs
						Ek(size(nn,1)+n)=Ek(size(nn,1)+n)+E_*sm
						Ok(size(nn,1)+n)=Ok(size(nn,1)+n)+O_*sm
					enddo
				enddo
			enddo
			!$omp end parallel do
			!E_=0d0
			!O_=0d0
			!do n=1,Ns
				!sm=2d0*(cos(brizon%ok(n,1))-cos(brizon%ok(n,2)))
				!E_=E_+Ek(size(nn,1)+n)*sm
				!O_=O_+Ok(size(nn,1)+n)*sm
			!enddo
			!write(*,*)E_,O_
		endif
		!$omp parallel do private(i,j,c,pb,O_,E_,sg,sgn,k,m)
		do i1=1,ubound(nn,1)
			Ek(i1)=0d0
			Ok(i1)=0d0
			if(get_row(Ecfg,nn(i1,:),m,sgn)) then
				call get_pb(shape(0),m(1:m(0)),O_,WA,iA,AW,WAW,wf)
				O_=O_*sgn
				E_=0d0
				do l=1,size(t)
					do n=1,size(latt%sb(1,1)%nb(l)%bd)
						i=latt%sb(1,1)%nb(l)%bd(n)%i(1)
						j=latt%sb(1,1)%nb(l)%bd(n)%i(2)
						do p=1,2
							if(get_row(cfg,(/i,-j/),k,sg,shape(0))) then
								call get_pb(k(1:k(0)),m(1:m(0)),pb,WA,iA,AW,WAW,wf)
								E_=E_-pb*t(l)*latt%sb(1,1)%nb(l)%bd(n)%bdc*sg*sgn
							endif
							if(get_row(cfg,(/-i-Ns,j+Ns/),k,sg,shape(0))) then
								call get_pb(k(1:k(0)),m(1:m(0)),pb,WA,iA,AW,WAW,wf)
								E_=E_-pb*t(l)*conjg(latt%sb(1,1)%nb(l)%bd(n)%bdc)*sg*sgn
							endif
							if(l==1) then
								if(get_row(cfg,(/i,i+Ns,-j-Ns,-j/),k,sg,shape(0))) then
									call get_pb(k(1:k(0)),m(1:m(0)),pb,WA,iA,AW,WAW,wf)
									E_=E_+pb*0.5d0*DJ*sg*sgn
								endif
							endif
							j=latt%sb(1,1)%nb(l)%bd(n)%i(1)
							i=latt%sb(1,1)%nb(l)%bd(n)%i(2)
						enddo
						if(l==1) then
							c(:,1)=abs((1-sign(1,-cfg(i::Ns)))/2-(/0,1/))
							c(:,2)=abs((1-sign(1,-cfg(j::Ns)))/2-(/0,1/))
							E_=E_+O_*(V*sum(c(:,1))*sum(c(:,2))+0.25d0*DJ*(c(1,1)-c(2,1))*(c(1,2)-c(2,2)))
						endif
					enddo
				enddo
				Ek(i1)=conjg(E_)
				Ok(i1)=conjg(O_)
			endif
		enddo
		!$omp end parallel do
		Oq=1d0/sum(Ok(:ubound(nn,1))*conjg(Ok(:ubound(nn,1))))
		O=0d0; E=0d0
		call her(O,Ok,alpha=real(Oq))
		call gerc(E,Ok,Ek,alpha=Oq)
	end subroutine
end module

module vmc_main
	use model
	use lapack95, only : hegv, geev, hegv
	implicit none
	type t_mc
		integer :: sg,ne(2),hot,samp,step
		integer, allocatable :: cfg(:),Ecfg(:),nn(:,:)
		real(8) :: phy(8)
		complex(8), allocatable :: S(:,:),Ok2(:,:),Ek2(:,:),psi0(:),wf(:,:),dwf(:,:,:)
		real(8), allocatable :: g(:)
		real(8) :: q(3)
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
		class(t_mc) :: self
		logical :: init_cfg
		complex(8) :: H(Ns*spin,Ns*spin),cH(size(H,1),size(H,2)),D(size(H,1),size(H,2),sum(var(1:vn)%n)),Uk(size(H,1),size(H,2)),psi0_(Ns*Ns)
		real(8) :: E(size(H,1))
		integer :: l,i,j,n,n1,n2,nk(size(H,1),size(H,2)),nn_(0:Ns*Ns,2),kq(Ns),Ecfg_(Ns*spin,2)
		type(t_mysort) :: mysort(size(E))

		if(init_cfg) then
			if(allocated(self%cfg)) deallocate(self%cfg)
			allocate(self%cfg(Ns*spin))
			self%cfg=0
			self%cfg(Ns+1:Ns+self%ne(2))=(/(n,n=1,self%ne(2))/)
			call fisher_yates_shuffle(self%cfg(Ns+1:))
			n=self%ne(2)
			do i=1,Ns
				if(self%cfg(i+Ns)/=0) then
					n=n+1
					self%cfg(i)=n
				endif
				if(n==sum(self%ne)) exit
			enddo
		endif
		if(allocated(self%Ecfg)) deallocate(self%Ecfg)
		allocate(self%Ecfg(Ns*spin))
		self%Ecfg=0

		call Hamilton(var,H)

		Uk=0d0
		do i=1,Ns
			do j=1,Ns
				Uk(j,i)=exp(-img*sum(latt%sb(1,1)%nb(0)%bd(j)%r*brizon%ok(i,:)))
			enddo
		enddo
		Uk(Ns+1:,Ns+1:)=Uk(:Ns,:Ns)
		H=matmul(transpose(conjg(Uk)),matmul(H,Uk))*(1d0/Ns)

		kq=0
		do i=1,Ns
			do j=1,Ns
				if(sum(abs(mod(brizon%ok(i,:2)/pi-brizon%ok(j,:2)/pi+(/2d0,2d0/)+1d-7,(/2d0,2d0/))-self%q(:2)/pi))<1d-5) then
					!if(sum(abs(mod(brizon%ok(i,:2)+brizon%ok(j,:2)+pi*(/2d0,2d0/),pi*(/2d0,2d0/))-q(:2)))<1d-7) then
					if(kq(i)/=0) then
						write(*,*)"warn! kq"
					else
						kq(i)=j
					endif
				endif
			enddo
			if(kq(i)==0) then
				write(*,*)self%q(:2)/pi,"warn! kq2"
				stop
			endif
		enddo
		nk=0
		if(any(abs(var(:)%sg)==2)) then
			do i=1,size(brizon%k,1)
				!call heev(H(i::size(brizon%k,1),i::size(brizon%k,1)),E(i::size(brizon%k,1)),"V")
				call mat_diag(H(i::size(brizon%k,1),i::size(brizon%k,1)),E(i::size(brizon%k,1)))
				nk(i::size(brizon%k,1),i::size(brizon%k,1))=1
				!write(30,"(es12.4$)")brizon%k(i,:2),E(i::size(brizon%k,1)),abs(H(i,i::size(brizon%k,1)))
				!write(30,"(x)")
			enddo
			mysort%val=E
			mysort%idx=(/(i,i=1,size(E))/)
			call qsort(mysort)
			do i=1,sum(self%ne)
				self%Ecfg(mysort(i)%idx)=i
			enddo
			if(abs(mysort(sum(self%ne))%val-mysort(sum(self%ne)+1)%val)<1d-7) then
				write(*,*)"warning, close-shell condition is not satisfied"
			endif
		else
			do i=1,size(brizon%k,1)
				call mat_diag(H(i:Ns:size(brizon%k,1),i:Ns:size(brizon%k,1)),E(i:Ns:size(brizon%k,1)))
				call mat_diag(H(Ns+i::size(brizon%k,1),Ns+i::size(brizon%k,1)),E(Ns+i::size(brizon%k,1)))
				nk(i:Ns:size(brizon%k,1),i:Ns:size(brizon%k,1))=1
				nk(Ns+i::size(brizon%k,1),Ns+i::size(brizon%k,1))=1
			enddo
			mysort%val=E
			mysort%idx=(/(i,i=1,size(E))/)
			call qsort(mysort(:Ns))
			call qsort(mysort(Ns+1:))
			do i=1,self%ne(1)
				self%Ecfg(mysort(i)%idx)=i
			enddo
			do i=1,self%ne(2)
				self%Ecfg(mysort(i+Ns)%idx)=i+self%ne(1)
			enddo
			if(abs(mysort(self%ne(1))%val-mysort(self%ne(1)+1)%val)<1d-7.or.abs(mysort(Ns+self%ne(2))%val-mysort(Ns+self%ne(2)+1)%val)<1d-7) then
				write(*,*)"warning, close-shell condition is not satisfied"
			endif
		endif

		if(self%sg==3) then
			Ecfg_(:,1)=(1-sign(1,-self%Ecfg))/2
			Ecfg_(:,2)=Ecfg_(:,1)
			psi0_=0d0
			l=0
			nn_(1:,:)=0
			do i=1,Ns
				do n2=1,size(E)
					if(Ecfg_(n2,1)==1.and.nk(kq(i),n2)>0) then
						Ecfg_(n2,1)=0
					o:	do n1=1,size(E)
							if(Ecfg_(n1,1)==0.and.nk(i,n1)>0) then
								Ecfg_(n1,1)=1
								do j=1,l
									Ecfg_(abs(nn_(j,:)),2)=(/0,1/)
									if(all(Ecfg_(:,1)==Ecfg_(:,2))) then
										psi0_(j)=psi0_(j)+2d0*(cos(brizon%ok(i,1))-cos(brizon%ok(i,2)))*conjg(H(i,n1))*H(kq(i),n2)
										Ecfg_(nn_(j,:),2)=(/1,0/)
										Ecfg_(n1,1)=0
										cycle o
									endif
									Ecfg_(nn_(j,:),2)=(/1,0/)
								enddo
								Ecfg_(n1,1)=0
								l=l+1
								nn_(l,:)=(/n2,n1/)
								psi0_(l)=psi0_(l)+2d0*(cos(brizon%ok(i,1))-cos(brizon%ok(i,2)))*conjg(H(i,n1))*H(kq(i),n2)
							endif
						enddo o
						Ecfg_(n2,1)=1
					endif
				enddo
			enddo
			!write(*,*)l
			!write(*,"(es12.4)")abs(psi0_(1:l))
			if(allocated(self%nn)) deallocate(self%nn,self%Ok2,self%Ek2,self%psi0)

			!allocate(self%nn(l,2),self%Ok2(l+Ns,l+Ns),self%Ek2(l+Ns,l+Ns),self%psi0(l+Ns))
			!self%psi0=0d0
			!do i=1,Ns
				!self%psi0(i+l)=2d0*(cos(brizon%ok(i,1))-cos(brizon%ok(i,2)))
			!enddo

			allocate(self%nn(l,2),self%Ok2(l+1,l+1),self%Ek2(l+1,l+1),self%psi0(l+1))
			self%psi0=0d0
			self%psi0(l+1)=1d0

			!self%psi0(1:l)=psi0_(1:l)

			self%nn(:,1)=nn_(1:l,1)
			self%nn(:,2)=-nn_(1:l,2)
		endif
		H=matmul(Uk,H)*sqrt(1d0/Ns)
		!write(*,"(2es12.4)")psi0
		!write(*,"(es12.4)")E(:sum(ne))
		!call Hamilton(var,cH)
		!!$omp critical
		!write(*,*)sum(abs(matmul(transpose(conjg(H)),matmul(cH,H))-diag(E)))
		!!$omp end critical
		!stop

		if(allocated(self%wf)) deallocate(self%wf)
		allocate(self%wf(Ns*spin,Ns*spin))
		self%wf=H

		if(self%sg==2) then
			if(allocated(self%dwf)) deallocate(self%dwf,self%g,self%S)
			allocate(self%dwf(Ns*spin,Ns*spin,vn),self%g(sum(var(1:)%n)))
			allocate(self%S(size(self%g),size(self%g)))
			cH=transpose(conjg(H))
			self%dwf=0d0
			call dHamilton(var(1:vn),H,cH,D)

			do l=1,size(self%dwf,3)
				select case(opt)
				case(1)
					!$omp parallel do
					do i=1,size(E)
						do j=1,size(E)
							if(abs(E(i)-E(j))<1d-8) then
								cycle
							endif
							self%dwf(:,i,l)=self%dwf(:,i,l)+D(j,i,l)*H(:,j)/(E(i)-E(j))
						enddo
					enddo
					!$omp end parallel do
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
		!write(*,*)"ini_wf finished"
	end subroutine
	subroutine do_vmc(self)
		class(t_mc) :: self
		real(8) :: eg(size(self%g)),g_(size(self%g))
		integer :: k,l,l1,l2,seed(n_omp)
		type(randomNumberSequence) :: rnd(n_omp)
		do k=1,n_omp
			call random_number(seed(k),4285441)
			call mt_init_random_seed(rnd(k),seed(k))
		enddo
		select case(self%sg)
		case(2)
			self%phy=0d0
			self%S=0d0; self%g=0d0
			!$omp parallel do if(n_omp>1)
			do k=1,n_omp
				call self%do_mc(rnd(k))
			enddo
			!$omp end parallel do
			self%phy(ier)=sqrt(self%phy(ier)-self%phy(iE)**2)
			call heev(self%S,eg,'V')
			write(*,"(es12.4$)")var(1:)%val(1),self%phy(iE),self%phy(ier)
			write(*,"(i3$)")int(sign(1d0,self%g))
			write(*,"(A$)")' | '
			write(20,"(es12.4$)")var(1:)%val(1),self%phy(iE),self%phy(ier),self%g
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
			write(*,"(i3$)")int(sign(1d0,self%g))
			!grad=g
			self%g=self%g*Ns
			!stop
		case(1)
			self%phy=0d0
			!$omp parallel do if(n_omp>1)
			do k=1,n_omp
				call self%do_mc(rnd(k))
			enddo
			!$omp end parallel do
			self%phy(ier)=self%phy(ier)-self%phy(iE)**2
			!$omp critical
			write(*,"(es12.4$)")self%q(:2)/pi,var(1:)%val(1),self%phy(iE),self%phy(ier),self%phy(iCq)
			write(10,"(es12.4$)")var(1:)%val(1)
			write(*,"(x)")
			write(10,"(x)")
			!$omp end critical
		case(3)
			self%Ok2=0d0; self%Ek2=0d0
			!$omp parallel do if(n_omp>1)
			do k=1,n_omp
				call self%do_mc(rnd(k))
			enddo
			!$omp end parallel do
			!write(50,*)-1,size(psi0),Ns
			!write(50,*)E,psi0,Ok2,Ek2
			!!$omp end critical
			!$omp ordered
			call self%spect(70,0.05d0,(/-10d0,10d0/),5000)
			self%Ek2=transpose(conjg(self%Ek2))
			call self%spect(70,0.05d0,(/-10d0,10d0/),5000)
			!$omp end ordered
		end select
		do k=1,n_omp
			call finalize_RandomNumberSequence(rnd(k))
		enddo
	end subroutine
	subroutine do_mc(self,rnd)
		class(t_mc) :: self
		type(randomNumberSequence) :: rnd
		complex(8) :: iA(sum(self%ne),sum(self%ne)),A(size(iA,1),size(self%wf,2)),WA(size(self%wf,1),size(iA,2)),AW(size(iA,1),size(self%wf,2)),WAW(size(self%wf,1),size(self%wf,2)),wf(size(self%wf,1),size(self%wf,2)),dwf(size(self%dwf,1),size(self%dwf,2),size(self%dwf,3))
		integer :: cfgl(size(self%cfg)),Ecfgl(size(self%Ecfg)),icfg(size(A,1)),iEcfg(size(A,1))
		real(8) :: gp(size(self%g)),gl(size(gp))
		complex(8) :: Sp(size(gp),size(gp)),Op(size(gp)),Ol(size(gp)),Sl(size(gp),size(gp))
		complex(8) :: Ok2l(size(self%Ok2,1),size(self%Ok2,2)),Ek2l(size(Ok2l,1),size(Ok2l,2)),Ok2p(size(Ok2l,1),size(Ok2l,2)),Ek2p(size(Ok2l,1),size(Ok2l,2))
		integer :: nn(lbound(self%nn,1):ubound(self%nn,1),size(self%nn,2))
		complex(8) :: phyl(size(self%phy)),phyp(size(self%phy))
		real(8) :: rd,isamp
		complex(8) :: pb
		integer :: i,j,l,l1,l2,info,n,apt,samp,sg
		integer :: k(0:4),m(0:4),dcfg(0:4)
		logical :: is_update
		phyp=0d0; phyl=0d0
		Sp=0d0; gp=0d0; Op=0d0
		Ek2p=0d0; Ok2p=0d0
		cfgl=self%cfg
		Ecfgl=self%Ecfg
		wf=self%wf
		dwf=self%dwf
		nn=self%nn
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
			else
				exit
			endif
		enddo
		WA=matmul(wf(:,iEcfg),iA)
		AW=matmul(iA,wf(icfg,:))
		WAW=matmul(wf(:,iEcfg),AW)
		!!$omp critical
		!write(*,"(es16.5$)")sum(abs(matmul(A(:,iEcfg),iA)-diag(1d0,size(A,1)))),sum(abs(wf(:,iEcfg)-matmul(WA,A(:,iEcfg)))),sum(abs(matmul(A(:,iEcfg),AW)-wf(icfg,:)))
		!!$omp end critical
		n=0; apt=0; samp=0
		is_update=.true.
	lm:	do
			n=n+1
			call change(cfgl,icfg,dcfg,k,rnd)

			if(self%sg==3) then
				!if(.false.) then
				call get_Spb(cfgl,dcfg,Ecfgl,nn,self%q,pb,WA,iA,AW,WAW,wf)
				if(isnan(real(pb))) then
					write(*,*)pb
					stop
				endif
			else
				call get_pb(k(1:k(0)),shape(0),pb,WA)
				pb=pb*conjg(pb)
			endif
			!call random_number(rd)
			rd=getrandomreal(rnd)
			if(rd<real(pb)) then
				apt=apt+1

				cfgl(abs(dcfg(2:dcfg(0):2)))=cfgl(abs(dcfg(2:dcfg(0):2)))+cfgl(dcfg(1:dcfg(0):2))
				cfgl(dcfg(1:dcfg(0):2))=cfgl(abs(dcfg(2:dcfg(0):2)))-cfgl(dcfg(1:dcfg(0):2))
				cfgl(abs(dcfg(2:dcfg(0):2)))=cfgl(abs(dcfg(2:dcfg(0):2)))-cfgl(dcfg(1:dcfg(0):2))
				icfg(k(1:k(0):2))=k(2:k(0):2)

				if(self%sg==2.and.opt==1) then
					call update(k(1:k(0)),WA,iA)
				elseif(self%sg/=3) then
					call update(k(1:k(0)),WA)
				else
					call get_pb(k(1:k(0)),shape(0),pb,WA)
					if(abs(pb)<1d-6) then
						cfgl(abs(dcfg(2:dcfg(0):2)))=cfgl(abs(dcfg(2:dcfg(0):2)))+cfgl(dcfg(1:dcfg(0):2))
						cfgl(dcfg(1:dcfg(0):2))=cfgl(abs(dcfg(2:dcfg(0):2)))-cfgl(dcfg(1:dcfg(0):2))
						cfgl(abs(dcfg(2:dcfg(0):2)))=cfgl(abs(dcfg(2:dcfg(0):2)))-cfgl(dcfg(1:dcfg(0):2))
						do i=1,dcfg(0)
							if(cfgl(abs(dcfg(i)))/=0) icfg(cfgl(abs(dcfg(i))))=abs(dcfg(i))
						enddo
						write(*,*)"warn update, matrix is singularity, skip this configure"
						cycle lm
					else
						call update(k(1:k(0)),WA,iA,AW,WAW,A,wf)
					endif
					!if(mod(n,10000)==0) then
					!write(*,"(i7$)")n
					!write(*,"(es12.4$)")abs(pb),sum(abs(matmul(A(:,:sum(ne)),iA)-diag(1d0,size(A,1)))),sum(abs(wf(:,:sum(ne))-matmul(WA,A(:,:sum(ne))))),sum(abs(matmul(A(:,:sum(ne)),AW)-wf(icfg,:))),sum(abs(WAW-matmul(wf(:,:sum(ne)),AW)))
					!write(*,"(x)")
					!endif
				endif
				is_update=.true.
				!if(mod(n,1024)==0) then
				!!$omp critical
				!write(*,*)n,abs(sum(abs(wf(:,iEcfg)-matmul(WA,wf(icfg,iEcfg)))))
				!!$omp end critical
				!endif
				if(mod(apt,1024)==0) then
					pb=sum(abs(wf(:,iEcfg)-matmul(WA,wf(icfg,iEcfg))))
					if(abs(pb)>1d-5) then
						!$omp critical
						write(*,*)"warn!! ",apt,abs(pb)
						!$omp end critical
						A=wf(icfg,:)
						iA=A(:,iEcfg)
						call mat_inv(iA,info)
						WA=matmul(wf(:,iEcfg),iA)
						if(self%sg==3) then
							AW=matmul(iA,wf(icfg,:))
							WAW=matmul(wf(:,iEcfg),AW)
						endif
						pb=sum(abs(wf(:,iEcfg)-matmul(WA,wf(icfg,iEcfg))))
						if(abs(pb)>1d-5) then
							!$omp critical
							write(*,*)"update err",abs(pb)
							!$omp end critical
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
						phyl(iE)=get_energy(cfgl,WA)
						!phyl%dsc=real(get_dsc(cfgl,WA))
						phyl(iCq)=get_charge(cfgl,WA,self%q)
						!phyl%Sq_zz=get_spin_zz(cfgl,q)
					case(2)
						phyl(iE)=get_energy(cfgl,WA)
						call get_O(icfg,iEcfg,WA,iA,dwf,Ol(:size(dwf,3)))
						do l1=1,size(Ol)
							do l2=1,size(Ol)
								Sl(l1,l2)=conjg(Ol(l1))*Ol(l2)
							enddo
						enddo
						gl=real(phyl(iE)*Ol)
					case(3)
						call get_overlap(cfgl,Ecfgl,nn,self%q,WA,iA,AW,WAW,wf,Ok2l,Ek2l)
					end select
				endif
				phyp=phyp+phyl
				select case(self%sg)
				case(2)
					gp=gp+gl; Op=Op+Ol; Sp=Sp+Sl
				case(3)
					Ek2p=Ek2p+Ek2l
					Ok2p=Ok2p+Ok2l
				end select
				if(mod(samp,1024)==0.and.self%sg==3) then
					!$omp critical
					write(*,"(i7$)")samp
					write(*,"(es16.5$)")real(apt)/n,sum(abs(Ek2p-conjg(transpose(Ek2p)))/samp)!,sum(abs(matmul(A(:,iEcfg),iA)-diag(1d0,size(A,1)))),sum(abs(wf(:,:sum(ne))-matmul(WA,A(:,:sum(ne))))),sum(abs(matmul(A(:,:sum(ne)),AW)-wf(icfg,:))),sum(abs(WAW-matmul(wf(:,:sum(ne)),AW)))
					write(*,"(x)")
					rewind(50)
					write(50,*)samp,size(self%psi0),Ns
					write(50,*)self%phy(iE),self%phy(iCq),self%psi0,Ok2p/samp,Ek2p/samp
					!$omp end critical
				endif
				!read(*,*)
				if(samp==self%samp) then
					exit
				endif
			endif
		enddo lm
		! check
		pb=sum(abs(wf(:,iEcfg)-matmul(WA,wf(icfg,iEcfg))))
		if(abs(pb)>1d-5) then
			!$omp critical
			write(*,*)"warn!!!!!",abs(pb)
			!$omp end critical
		endif

		phyp=phyp*isamp
		Sp=Sp*isamp; gp=gp*isamp; Op=Op*isamp
		do l1=1,size(Op)
			do l2=1,size(Op)
				!S(l1,l2)=2d0*SE(l1,l2)-O(l1)*g(l2)-O(l2)*g(l1)-2d0*E*S(l1,l2) ! maybe works
				!S(l1,l2)=2d0*(SE(l1,l2)-S(l1,l2)*E-O(l1)*g(l2)-O(l2)*g(l1))
				Sp(l1,l2)=Sp(l1,l2)-conjg(Op(l1))*Op(l2)
			enddo
		enddo
		gp=2d0*(gp-real(phyp(iE)*Op))
		Ek2p=Ek2p*isamp
		Ok2p=Ok2p*isamp

		!$omp critical
		self%cfg=cfgl
		select case(self%sg)
		case(1:2)
			self%phy=self%phy+real(phyp)/n_omp
			self%phy(ier)=self%phy(ier)+real(phyp(iE)*conjg(phyp(iE)))/n_omp
			self%S=self%S+Sp/n_omp; self%g=self%g+gp/n_omp
		case(3)
			self%Ek2=self%Ek2+Ek2p/n_omp
			self%Ok2=self%Ok2+Ok2p/n_omp
		end select
		!$omp end critical
	end subroutine
	subroutine spect(self,ut,gm,omg,m)
		class(t_mc) :: self
		integer :: ut,m
		real(8) :: gm,omg(:)
		complex(8) :: Sq(m),H_(size(self%Ek2,1),size(self%Ek2,2)),O_(size(self%Ek2,1),size(self%Ek2,2)),Om(size(self%Ek2,1),size(self%Ek2,2)),psi0_(size(self%psi0))
		real(8) :: domg,EO(size(O_,1)),Eg(size(O_,1)),norm
		integer :: l,i,j,n0
		domg=(omg(2)-omg(1))/m
		Sq=0d0
		O_=self%Ok2
		H_=self%Ek2
		psi0_=self%psi0
		!H_=0.5d0*(H_+transpose(conjg(H_)))
		n0=1
		call heev(O_,EO,'V')
		!$omp critical
		do i=1,size(EO)
			if(EO(i)>1d-8) then
				n0=i
				write(*,*)n0,size(EO)
				exit
			endif
		enddo
		!check truncate
		H_=matmul(transpose(conjg(O_)),matmul(H_,O_))
		write(*,*)maxval(abs(H_(:n0-1,:)))
		!$omp end critical
		Om(n0:,n0:)=diag(EO(n0:))
		call hegv(H_(n0:,n0:),Om(n0:,n0:),Eg(n0:),jobz="V")
		psi0_=matmul(transpose(conjg(O_)),psi0_)
		psi0_(n0:)=matmul(transpose(conjg(H_(n0:,n0:))),EO(n0:)*psi0_(n0:))
		psi0_(n0:)=conjg(psi0_(n0:))*psi0_(n0:)
		do l=1,m
			Sq(l)=Sq(l)+sum(psi0_(n0:)/(omg(1)+domg*l-Eg(n0:)+self%phy(iE)*Ns+img*gm))
		enddo
		norm=self%phy(iCq)/real(sum(psi0_(n0:)))
		write(*,"(es12.4$)")-sum(imag(Sq)*domg)/pi,real(sum(psi0_(n0:))),self%phy(iCq)
		write(ut,"('#q=',3es12.4,2i5)")self%q/pi,i,size(EO)
		do l=1,m
			write(ut,"(es17.9$)")omg(1)+domg*l,Sq(l)*norm
			write(ut,"(x)")
		enddo
		write(ut,"(x)")
	end subroutine
	subroutine do_var(self)
		class(t_mc) :: self
		real(8) :: x(sum(var(1:)%n)),dx,El(150),er,pgrad(size(x))
		integer :: i,hot
		logical :: init_cfg
		!dx=0.03d0
		dx=0.1d0
		hot=self%hot
		i=0
		if(size(x)==0) then
			return
		endif
		x=put(var(1:))
		init_cfg=.true.
		do 
			i=i+1
			if(allocated(self%g)) pgrad=self%g
			call get(var(1:),x)
			write(*,"(i4$)")i
			call self%init(init_cfg)
			call self%do_vmc()
			El(i)=self%phy(iE)
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
			x=x-self%g*dx
			self%hot=128
			init_cfg=.false.
			write(*,"(x)")
			write(20,"(x)")
		enddo
		self%hot=hot
		write(*,*)"finished"
		!E=minval(El(:i))
	end subroutine
end module

program main
	use vmc_main
	implicit none
	logical :: f
	integer :: i,j,k
	type(t_mc) :: mc
	f=openfile(101,"../data/lattice.dat")
	f=openfile(10,"../data/grad.dat")
	f=openfile(20,"../data/var.dat")
	!f=openfile(30,"../data/alg.dat")
	!f=openfile(40,"../data/phyvar.dat")
	f=openfile(50,"../data/matrix_save.dat")
	f=openfile(70,"../data/spect_qusi_stripe_debug.dat")


	rewind(50)
	read(50,*)i,k,Ns
	allocate(mc%psi0(k),mc%Ok2(k,k),mc%Ek2(k,k))
	read(50,*)mc%phy(iE),mc%phy(iCq),mc%psi0,mc%Ok2,mc%Ek2
	write(*,*)mc%phy(iE),mc%phy(iCq),Ns,i,k
	do i=1,size(mc%Ok2,1)
		do j=1,i-1
			mc%Ok2(i,j)=conjg(mc%Ok2(j,i))
		enddo
	enddo
	write(*,*)sum(matmul(mc%psi0(k-Ns+1:k),mc%Ek2(k-Ns+1:k,k-Ns+1:k))*mc%psi0(k-Ns+1:k))/sum(matmul(mc%psi0(k-Ns+1:k),mc%Ok2(k-Ns+1:k,k-Ns+1:k))*mc%psi0(k-Ns+1:k))
	call mc%spect(70,0.05d0,(/-10d0,10d0/),5000)
	mc%Ek2=transpose(conjg(mc%Ek2))
	call mc%spect(70,0.05d0,(/-10d0,10d0/),5000)
	stop

	call initial()

	mc%hot=1024
	mc%step=Ns
	mc%samp=1024

	!var(1:)%val(1)=(/0d0,0.36d0,0.44d0/)
	!var(1:)%val(1)=(/0d0,0.36d0,0d0/)
	!var(1:)%val(1)=(/-0.39d0,0.32d0,-0.205d0/)
	!var(1:)%val(1)=(/-0.46d0,0.300d0,-0.268d0/)
	var(1:)%val(1)=(/-1.39752E+00,8.88530E-01*2E0,1.32678E-02,-2.00021E-01/)
	!var(1:)%val(1)=(/-1.4051E+00,2.4353E+00,-1.9875E-01,-2.6827E-01/)
	!var(1:)%val(1)=(/-1.2070E+00,1.8269E+00,-1.3504E-01,-1.5381E-01/)
	!var(1:)%val(1)=(/-1.53531E+00,8.07142E-01*2E0,-6.06570E-02,-2.62224E-01,-2.15494E-01/)
	!var(1:)%val(1)=(/-8.0718E-01,2.0919E-01,-9.9164E-02/) ! 16x16
	!var(1:)%val(1)=(/-6.4806E-01,1.7347E-01,-1.7155E-01/) ! 8x8
	var(1:)%val(1)=(/-7.7283E-01,2.0321E-01,-1.7129E-01/) ! 8x8
	!var(1:)%val(1)=(/-8.6811E-01,2.0951E-01,-9.7832E-02/) ! 12x12
	!call variation()
	!stop
	!call export_data(101)
	!call omp_set_nested(.true.)
	!call mkl_set_dynamic(.true.)
	!call mkl_domain_num_threads("MKL_DOMAIN_BLAS=2")

	n_omp=12
	mc%ne(1)=Ns/2-4
	mc%ne(2)=Ns-mc%ne(1)
	mc%samp=1024*4
	mc%sg=2
	!call mc%do_var()
	!stop

	j=nint(sqrt(real(Ns)))/2
	n_omp=1
	!!$omp parallel do ordered firstprivate(mc) schedule(static,1)
	!!!do i=1,j*3
	!do i=j/2+1,j+j/2
		!if(i<=j) then
			!mc%q=(/0d0,0d0,0d0/)+((/pi,pi,0d0/)-(/0d0,0d0,0d0/))/j*mod(i-1,j)
		!elseif(i<=2*j) then
			!mc%q=(/pi,pi,0d0/)+((/0d0,pi,0d0/)-(/pi,pi,0d0/))/j*mod(i-1,j)
		!else
			!mc%q=(/0d0,pi,0d0/)+((/0d0,0d0,0d0/)-(/0d0,pi,0d0/))/j*mod(i-1,j)
		!endif
		mc%q=(/0d0,0d0,0d0/)

		mc%sg=1
		mc%ne(1)=Ns/2-4
		mc%ne(2)=Ns-mc%ne(1)
		mc%samp=1024*8*8
		call mc%init(.true.)
		call mc%do_vmc()

		call omp_set_nested(.true.)
		mc%sg=3
		mc%samp=1024*8*8*8
		call mc%init(.true.)
		call mc%do_vmc()

		!stop
	!enddo
	!!$omp end parallel do
end program
