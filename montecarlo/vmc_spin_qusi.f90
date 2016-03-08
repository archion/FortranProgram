module global
	use M_hamilton_m
	use ifport
	implicit none
	real(8) :: t(2)=(/1d0,-0.25d0/)
	real(8), parameter :: DJ=1d0/3d0,V=0d0,U=0d0
	integer :: Nmc(4)
	integer :: n_omp=1
	integer, parameter :: opt=1 ! 1: Tr(AdA)
								! 2: O=cicj
	integer :: ne(2),vn
	integer, allocatable :: cfg(:,:)
	real(8) :: E,dsc,af,ddw,Sq_pm,Sq_zz,er
	real(8), allocatable :: g(:)
	complex(8), allocatable :: S(:,:)
	real(8), allocatable :: grad(:)
	complex(8), allocatable :: Ok2(:,:),Ek2(:,:)
	complex(8), allocatable :: psi0(:)
	integer :: map(0:3),Ns
	real(8) :: q(3)=(/pi,pi,0d0/)
	integer, allocatable :: kq(:)
	integer :: mc_sg=1 ! 1 static physical 
					   ! 2 energy and grad
					   ! 3 dynamic physical
contains
	subroutine initial()
		integer :: i,l
		allocate(var(-1000:1000))
		!call init_random_seed()
		! lattice
		latt%a1=(/1d0,0d0,0d0/)
		latt%a2=(/0d0,1d0,0d0/)
		!latt%T1=(/9d0,1d0,0d0/)
		!latt%T2=(/-1d0,9d0,0d0/)
		latt%T1=(/1d0,0d0,0d0/)*10
		latt%T2=(/0d0,1d0,0d0/)*10
		latt%bdc=(/1d0,-1d0,0d0/)
		allocate(latt%sub(1,3))
		latt%sub(1,:)=(/0d0,0d0,0d0/)
		call latt%gen_latt()
		call latt%gen_bond(size(t))
		brizon%n1=1
		brizon%n2=1
		call latt%gen_brizon(brizon)
		call latt%gen_origin_brizon(brizon_o)
		if(abs(latt%bdc(2)+1d0)<1d-7) then
			brizon_o%k(:,2)=brizon_o%k(:,2)+pi/latt%T2(2)
		endif
		if(abs(latt%bdc(1)+1d0)<1d-7) then
			brizon_o%k(:,1)=brizon_o%k(:,1)+pi/latt%T1(1)
		endif
		!call check_lattice(101)
		write(*,*)"Total site number is: ",latt%Ns

		! cp
		call gen_var(sg=1,nb=0)
		var(iv(0))%bd=-1d0
		var(iv(0))%val=-8.1583d-01

		! dsc
		call gen_var(sg=2,nb=1)
		do i=1,size(var(iv(0))%bd)
			var(iv(0))%bd(i)=dwave(i)
		enddo
		var(iv(0))%val=1.3409d-01

		!! ssc
		!call gen_var(sg=-2,nb=0)
		!do i=1,size(var(iv(0))%bd)
			!var(iv(0))%bd(i)=1d0
		!enddo
		!var(iv(0))%val=4d-4

		!! ddw
		!call gen_var(sg=3,nb=1)
		!do i=1,size(var(iv(0))%bd)
			!var(iv(0))%bd(i)=ab(latt%sb(1)%nb(var(iv(0))%nb)%bd(i)%i(1))*dwave(i)*img
		!enddo
		!var(iv(0))%val=1d-3

		!! sdw
		!call gen_var(sg=4,nb=0)
		!do i=1,size(var(iv(0))%bd)
			!var(iv(0))%bd(i)=ab(i)
		!enddo
		!var(iv(0))%val=1d-1

		do l=2,size(t)
			call gen_var(sg=3,nb=l)
			var(iv(0))%bd=-1d0
			var(iv(0))%val=-3.4320d-01
		enddo

		do l=1,1!size(t)
			call gen_var(sg=-3,nb=l)
			var(iv(0))%bd=-1d0
			var(iv(0))%val=t(l)
		enddo

		vn=iv(1)

		call var_shrink()

		Ns=latt%Ns
		allocate(cfg(Ns,3))
		allocate(grad(sum(var(1:)%n)))
		allocate(g(size(grad)))
		allocate(S(size(grad),size(grad)))
		allocate(Ok2(Ns,Ns),Ek2(Ns,Ns),psi0(Ns))
		allocate(kq(Ns))
		map=(/b"10",b"11",b"00",b"01"/)
	end subroutine
	subroutine ini_wf(wf,dwf)
		complex(8), optional :: wf(:,:)
		complex(8), optional :: dwf(:,:,:)
		complex(8) :: H(Ns*spin,Ns*spin),cH(size(H,1),size(H,2)),D(size(H,1),size(H,2),sum(var(1:vn)%n)),DQ(size(H,1)-sum(ne),sum(ne)),Uk(size(H,1),size(H,2)),tmp
		real(8) :: E(size(H,1))
		integer :: l,i,j
		call Hamilton(var,H)

		Uk=0d0
		do i=1,Ns
			Uk(:Ns,i)=exp(-img*(latt%i2r(:,1)*brizon_o%k(i,1)+latt%i2r(:,2)*brizon_o%k(i,2)))
		enddo
		do i=1,Ns
			do j=1,Ns
				if(sum(abs(mod(brizon_o%k(i,:2)+brizon_o%k(j,:2)+2d0*(/pi,pi/),2d0*(/pi,pi/))-q(:2)))<1d-7) then
					kq(i)=j
				endif
			enddo
		enddo
		Uk(Ns+1:,Ns+1:)=Uk(:Ns,:Ns)
		H=matmul(transpose(conjg(Uk)),matmul(H,Uk))/Ns
		do i=1,Ns
			call heev(H(i::Ns,i::Ns),E(i::Ns),"V")
		enddo
		do i=1,Ns
			psi0(i)=H(i,Ns+i)*H(Ns+kq(i),Ns+kq(i))
		enddo
		H=matmul(Uk,H)/sqrt(real(Ns))

		!if(any(abs(var(:)%sg)==1.and.any(abs(var(:)%sg)==2))) then
			!call heev(H,E,"V")
		!elseif(all(abs(var(:)%sg)/=1).and.any(abs(var(:)%sg)==2)) then
			!wf=H
			!call heev(wf(:Ns,:Ns),E(:Ns))
			!H(:Ns,:Ns)=H(:Ns,:Ns)-diag((E(ne(1))+E(ne(1)+1))/2d0,Ns)
			!H(Ns+1:,Ns+1:)=H(Ns+1:,Ns+1:)+diag((E(ne(1))+E(ne(1)+1))/2d0,Ns)
			!call heev(H,E,"V")
		!else
			!call heev(H(:Ns,:Ns),E(:Ns),"V")
			!call heev(H(Ns+1:,Ns+1:),E(Ns+1:),"V")
			!do i=1,ne(2)
				!call swap(H(:,ne(1)+i),H(:,Ns+i))
				!call swap(E(ne(1)+i),E(Ns+i))
			!enddo
		!endif

		if(present(wf)) wf=H

		if(present(dwf)) then
			cH=transpose(conjg(H))
			dwf=0d0
			call dHamilton(var(1:vn),H,cH,D)

			do l=1,size(dwf,3)
				select case(opt)
				case(1)
					!$omp parallel do
					do i=1,size(E)
						do j=1,size(E)
							if(abs(E(i)-E(j))<1d-8) then
								cycle
							endif
							dwf(:,i,l)=dwf(:,i,l)+D(j,i,l)*H(:,j)/(E(i)-E(j))
						enddo
					enddo
					!$omp end parallel do
				case(2)
					!$omp parallel do collapse(2)
					do i=sum(ne)+1,size(E)
						do j=1,sum(ne)
							if(abs(E(i)-E(j))<1d-8) then
								DQ(i-sum(ne),j)=0d0
								cycle
							endif
							DQ(i-sum(ne),j)=D(i,j,l)/(E(j)-E(i))
						enddo
					enddo
					!$omp end parallel do
					!$omp parallel do
					do i=1,size(dwf,1)
						dwf(i,:,l)=matmul(matmul(H(i,sum(ne)+1:),DQ),cH(:sum(ne),:))
					enddo
					!$omp end parallel do
				end select
			enddo
		endif
		write(*,*)"ini_wf finished"
	end subroutine
	function ab(i)
		integer :: i
		real(8) :: ab
		ab=(-1d0)**(mod(nint(sum(latt%sb(1)%nb(0)%bd(i)%r(:2))),2))
	end function
end module

module mc_utility
	use global
	use blas95, only : gemm, gerc
	implicit none
contains
	subroutine get_pb(k,m,pb,WA,iA,AW,WAW,wf)
		integer, allocatable :: k(:)
		integer :: m(:)
		complex(8) :: WA(:,:),pb
		complex(8), optional :: AW(:,:),WAW(:,:),iA(:,:),wf(:,:)
		integer :: i,j,n
		integer :: k_(size(k)/2+size(m)/2),m_(size(k_))
		complex(8) :: A(size(k_,1),size(k_,1))
		if(.not.allocated(k)) then
			pb=0d0
			return
		endif
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
	subroutine two_all(cfg,dcfg)
		integer :: cfg(:,:)
		integer, allocatable :: dcfg(:)
		integer :: i,j,l
		if(allocated(dcfg)) deallocate(dcfg)
		call random_number(i,size(cfg,1))
		do 
			call random_number(j,size(cfg,1))
			call random_number(l,2)
			l=l-1
			if(cfg(i,1)<cfg(j,1)) call swap(i,j)
			if(ibits(cfg(i,1),l,1)/=ibits(cfg(j,1),l,1)) then
				allocate(dcfg(2))
				dcfg=(/i+Ns*l,-j-Ns*l/)
				exit
			endif
		enddo
	end subroutine
	subroutine two(cfg,dcfg)
		integer :: cfg(:,:)
		integer, allocatable :: dcfg(:)
		integer :: i,j
		if(allocated(dcfg)) deallocate(dcfg)
		call random_number(i,size(cfg,1))
		do 
			call random_number(j,size(cfg,1))
			if(cfg(i,1)<cfg(j,1)) call swap(i,j)
			select case(ieor(cfg(i,1),cfg(j,1)))
			case(1)
				allocate(dcfg(2))
				dcfg=(/i,-j/)
			case(2)
				allocate(dcfg(2))
				dcfg=(/i+Ns,-j-Ns/)
			case(3)
				allocate(dcfg(4))
				dcfg=(/i,-j,i+Ns,-j-Ns/)
			case(0)
				cycle
			end select
			exit
		enddo
	end subroutine
	subroutine get_row(cfg,dcfg,k,sg)
		integer :: cfg(:,:),dcfg(:),sg
		integer, allocatable :: k(:),m(:)
		integer :: i,l,j1,j2,j,n,k_(10),c(size(cfg,1),3),dcfg_(size(dcfg))
		k_=0
		sg=1
		if(allocated(k)) deallocate(k)
		dcfg_=dcfg
		do i=1,size(dcfg_)-1
			if(mod(i,2)==1.and.dcfg_(i)<0) then
				do j=i+1,size(dcfg_)
					if(dcfg_(j)>0) then
						n=dcfg_(j)
						if(any(dcfg_(i:j-1)==-n)) stop "jump errr"
						dcfg_(i+1:j)=dcfg_(i:j-1)
						dcfg_(i)=n
						if(mod(j-i,2)==1) sg=-sg
						exit
					endif
				enddo
			elseif(mod(i,2)==0) then
				l=sign((abs(dcfg_(i-1))-1)/Ns+1,dcfg_(i-1))
				if(sign((abs(dcfg_(i))-1)/Ns+1,dcfg_(i))/=-l) then
					do j=i+1,size(dcfg_)
						if(sign((abs(dcfg_(j))-1)/Ns+1,dcfg_(j))==-l) then
							n=dcfg_(j)
							if(any(dcfg_(i:j-1)==-n)) stop "jump errr"
							dcfg_(i+1:j)=dcfg_(i:j-1)
							dcfg_(i)=n
							if(mod(j-i,2)==1) sg=-sg
							exit
						endif
					enddo
				endif
			endif
		enddo
		c(mod(abs(dcfg)-1,Ns)+1,:)=cfg(mod(abs(dcfg)-1,Ns)+1,:)
		n=0
		do i=1,size(dcfg_),2
			l=(abs(dcfg_(i))-1)/Ns
			j1=abs(dcfg_(i))-Ns*l
			j2=abs(dcfg_(i+1))-Ns*l
			if(j1==j2.and.btest(c(j1,1),l)) cycle
			if(btest(c(j1,1),l).and.(.not.btest(c(j2,1),l))) then
				c(j1,1)=ibclr(c(j1,1),l)
				c(j2,1)=ibset(c(j2,1),l)
				!if(any(map(c((/j1,j2/),1))==3)) then
					!return
				!endif
				k_(n+1)=c(j1,2+l)
				k_(n+2)=j2+Ns*l
				n=n+2
				call swap(c(j2,2+l),c(j1,2+l))
			else
				return
			endif
			do j=2,n-2,2
				if((j1+Ns*l)==k_(j)) then
					k_(j)=k_(n)
					n=n-2
					exit
				endif
			enddo
		enddo
		if(any(map(c(mod(abs(dcfg)-1,Ns)+1,1))==3)) then
			return
		endif
		allocate(k(n))
		k=k_(:n)
		!k(1)=raiseqq(sig$int)
	end subroutine
	function get_spin_pm(cfg,D,q)
		complex(8) :: D(:,:)
		integer :: cfg(:,:)
		real(8) :: q(3)
		complex(8) :: get_spin_pm
		complex(8) :: pb
		integer :: i,j,sg
		integer, allocatable :: k(:)
		get_spin_pm=0d0
		!$omp parallel if(.not.omp_in_parallel())
		!$omp do collapse(2) reduction(+:get_spin_pm) private(pb,k,sg)
		do i=1,Ns
			do j=1,Ns
				call get_row(cfg,(/j,j+Ns,-i-Ns,-i/),k,sg)
				call get_pb(k,shape(0),pb,D)
				get_spin_pm=get_spin_pm+pb*exp(img*sum(q*(latt%i2r(i,:)-latt%i2r(j,:))))*sg
			enddo
		enddo
		!$omp end do
		if(allocated(k)) deallocate(k)
		!$omp end parallel
		get_spin_pm=get_spin_pm/Ns
	end function
	function get_spin_zz(cfg,q)
		integer :: cfg(:,:)
		real(8) :: q(3)
		complex(8) :: get_spin_zz
		integer :: i,j
		get_spin_zz=0d0
		!$omp parallel do collapse(2) reduction(+:get_spin_zz) if(.not.omp_in_parallel())
		do i=1,Ns
			do j=1,Ns
				get_spin_zz=get_spin_zz+0.25d0*(ibits(map(cfg(i,1)),0,1)-ibits(map(cfg(i,1)),1,1))*(ibits(map(cfg(j,1)),0,1)-ibits(map(cfg(j,1)),1,1))*exp(img*sum(q*(latt%i2r(i,:)-latt%i2r(j,:))))
			enddo
		enddo
		!$omp end parallel do
		get_spin_zz=get_spin_zz/Ns
	end function
	function get_dsc(cfg,D)
		complex(8) :: D(:,:)
		integer :: cfg(:,:)
		complex(8) :: get_dsc
		complex(8) :: pb
		integer :: i1,j1,i2,j2,sg,l,n1,n2,p1,p2
		integer, allocatable :: k(:)
		get_dsc=0d0
		!$omp parallel if(.not.omp_in_parallel())
		!$omp do collapse(2) reduction(+:get_dsc) private(pb,k,sg,i1,j1,i2,j2)
		do n1=1,size(latt%sb(1)%nb(1)%bd)
			do n2=1,size(latt%sb(1)%nb(1)%bd)
				i1=latt%sb(1)%nb(1)%bd(n1)%i(1)
				j1=latt%sb(1)%nb(1)%bd(n1)%i(2)
				i2=latt%sb(1)%nb(1)%bd(n2)%i(1)
				j2=latt%sb(1)%nb(1)%bd(n2)%i(2)
				if(any(latt%sb(1)%nb(1)%bd(n1)%i==i2).or.any(latt%sb(1)%nb(1)%bd(n1)%i==j2)) cycle
				do p1=1,2
					do p2=1,2
						call get_row(cfg,(/-j1-Ns,i1,-i2,j2+Ns/),k,sg)
						call get_pb(k,shape(0),pb,D)
						get_dsc=get_dsc+pb*sg*dwave(n1)*dwave(n2)*latt%sb(1)%nb(1)%bd(n1)%bdc*latt%sb(1)%nb(1)%bd(n2)%bdc
						call swap(i2,j2)
					enddo
					call swap(i1,j1)
				enddo
			enddo
		enddo
		!$omp end do
		if(allocated(k)) deallocate(k)
		!$omp end parallel
		get_dsc=get_dsc/(Ns**2)
	end function
	function get_energy(cfg,D)
		complex(8) :: D(:,:)
		integer :: cfg(:,:)
		complex(8) :: get_energy
		complex(8) :: pb
		integer :: i,j,n,l,ci,cj,p,sg
		integer, allocatable :: k(:)
		get_energy=0d0
		!$omp parallel if(.not.omp_in_parallel())
		!$omp do collapse(2) reduction(+:get_energy) private(i,j,ci,cj,k,pb,sg)
		do l=1,size(t)
			do n=1,size(latt%sb(1)%nb(l)%bd)
				i=latt%sb(1)%nb(l)%bd(n)%i(1)
				j=latt%sb(1)%nb(l)%bd(n)%i(2)
				do p=1,2
					call get_row(cfg,(/i,-j/),k,sg)
					call get_pb(k,shape(0),pb,D)
					get_energy=get_energy-t(l)*pb*sg*latt%sb(1)%nb(l)%bd(n)%bdc!*-1
					call get_row(cfg,(/-i-Ns,j+Ns/),k,sg)
					call get_pb(k,shape(0),pb,D)
					get_energy=get_energy-t(l)*pb*sg*latt%sb(1)%nb(l)%bd(n)%bdc!*-1
					if(l==1) then
						call get_row(cfg,(/i,i+Ns,-j-Ns,-j/),k,sg)
						call get_pb(k,shape(0),pb,D)
						get_energy=get_energy+0.5d0*DJ*pb*sg!*-2
					endif
					j=latt%sb(1)%nb(l)%bd(n)%i(1)
					i=latt%sb(1)%nb(l)%bd(n)%i(2)
				enddo
				ci=map(cfg(i,1))
				cj=map(cfg(j,1))
				if(l==1) get_energy=get_energy+V*(ibits(ci,0,1)+ibits(ci,1,1))*(ibits(cj,0,1)+ibits(cj,1,1))+0.25d0*DJ*(ibits(ci,0,1)-ibits(ci,1,1))*(ibits(cj,0,1)-ibits(cj,1,1))!*-4
			enddo
		enddo
		!$omp end do
		if(allocated(k)) deallocate(k)
		!$omp end parallel
		get_energy=get_energy/Ns
	end function
	function get_A(cfg,wf,wcfg)
		integer :: cfg(:,:)
		integer, optional :: wcfg(:)
		complex(8) :: wf(:,:)
		complex(8) :: get_A(sum(ne),size(wf,2))
		integer :: n,l
		do l=0,1
			do n=1,Ns
				if(cfg(n,2+l)>0) then
					get_A(cfg(n,2+l),:)=wf(n+Ns*l,:)
					if(present(wcfg)) then
						wcfg(cfg(n,2+l))=n+Ns*l
					endif
				endif
			enddo
		enddo
	end function
	subroutine get_O(cfg,WA,iA,dwf,O)
		complex(8) :: WA(:,:),iA(:,:),dwf(:,:,:),O(:)
		integer :: cfg(:,:)
		integer :: i,Ns
		Ns=latt%Ns
		O=0d0
		select case(opt)
		case(1)
			!$omp parallel do reduction(+:O) if(.not.omp_in_parallel())
			do i=1,size(cfg,1)
				if(btest(cfg(i,1),0)) then
					O=O+matmul(iA(:,cfg(i,2)),dwf(i,:sum(ne),:))
				endif
				if(btest(cfg(i,1),1)) then
					O=O+matmul(iA(:,cfg(i,3)),dwf(i+Ns,:sum(ne),:))
				endif
			enddo
			!$omp end parallel do
		case(2)
			!$omp parallel do reduction(+:O) if(.not.omp_in_parallel())
			do i=1,size(cfg,1)
				if(btest(cfg(i,1),0)) then
					O=O+matmul(WA(:,cfg(i,2)),dwf(i,:,:))
				endif
				if(btest(cfg(i,1),1)) then
					O=O+matmul(WA(:,cfg(i,3)),dwf(i+Ns,:,:))
				endif
			enddo
			!$omp end parallel do
		end select
	end subroutine
	subroutine get_Spb(cfg,dcfg,pb,WA,iA,AW,WAW,wf)
		integer :: cfg(:,:),dcfg(:)
		complex(8) :: WA(:,:),AW(:,:),WAW(:,:),iA(:,:),wf(:,:)
		complex(8) :: pb
		complex(8) :: pb2(2),pbs(2)
		integer :: i,l,sg
		integer, allocatable :: k(:)
		!write(*,*)"get_Oq"
		pbs=0d0
		!$omp parallel if(.not.omp_in_parallel())
		!$omp do private(pb2,k,sg) reduction(+:pbs)
		do i=1,Ns
			call get_row(cfg,shape(0),k,sg)
			call get_pb(k,(/Ns+1,Ns+i,Ns+2,Ns+kq(i)/),pb2(1),WA,iA,AW,WAW,wf)
			pb2(1)=pb2(1)*sg
			call get_row(cfg,dcfg,k,sg)
			call get_pb(k,(/Ns+1,Ns+i,Ns+2,Ns+kq(i)/),pb2(2),WA,iA,AW,WAW,wf)
			pb2(2)=pb2(2)*sg
			!pbs=pbs+exp(img*sum(q*latt%i2r(i,:)))*pb2
			pbs=pbs+pb2*conjg(pb2)
		enddo
		!$omp end do
		if(allocated(k)) deallocate(k)
		!$omp end parallel
		!write(*,*)real(pbs)
		pb=pbs(2)/pbs(1)
	end subroutine
	subroutine get_overlap(cfg,WA,iA,AW,WAW,wf,O,E)
		integer :: cfg(:,:)
		complex(8) :: WA(:,:),AW(:,:),WAW(:,:),iA(:,:),wf(:,:),E(:,:),O(:,:)
		complex(8) :: Ek(size(E,1)),Ok(size(O,1)),pb,Oq,O_,E_
		real(8) :: iNs
		integer :: cfg_s(size(cfg,1)),i1,j1,i,j,l,n,Si(0:1),ci,cj,p,sg
		integer, allocatable :: k(:)
		iNs=1d0/Ns
		!Ok=0d0
		!Ek=0d0
		!$omp parallel if(.not.omp_in_parallel())
		!$omp do private(i,j,ci,cj,pb,O_,E_,k,sg)
		do i1=1,Ns
			call get_row(cfg,shape(0),k,sg)
			call get_pb(k,(/Ns+1,Ns+i1,Ns+2,Ns+kq(i1)/),O_,WA,iA,AW,WAW,wf)
			O_=O_*sg
			E_=0d0
			do l=1,size(t)
				do n=1,size(latt%sb(1)%nb(l)%bd)
					i=latt%sb(1)%nb(l)%bd(n)%i(1)
					j=latt%sb(1)%nb(l)%bd(n)%i(2)
					do p=1,2
						call get_row(cfg,(/i,-j/),k,sg)
						call get_pb(k,(/Ns+1,Ns+i1,Ns+2,Ns+kq(i1)/),pb,WA,iA,AW,WAW,wf)
						E_=E_-pb*t(l)*latt%sb(1)%nb(l)%bd(n)%bdc*sg
						call get_row(cfg,(/-i-Ns,j+Ns/),k,sg)
						call get_pb(k,(/Ns+1,Ns+i1,Ns+2,Ns+kq(i1)/),pb,WA,iA,AW,WAW,wf)
						E_=E_-pb*t(l)*latt%sb(1)%nb(l)%bd(n)%bdc*sg
						if(l==1) then
							call get_row(cfg,(/i,i+Ns,-j-Ns,-j/),k,sg)
							call get_pb(k,(/Ns+1,Ns+i1,Ns+2,Ns+kq(i1)/),pb,WA,iA,AW,WAW,wf)
							E_=E_+pb*0.5d0*DJ*sg
						endif
						j=latt%sb(1)%nb(l)%bd(n)%i(1)
						i=latt%sb(1)%nb(l)%bd(n)%i(2)
					enddo
					ci=map(cfg(i,1))
					cj=map(cfg(j,1))
					if(l==1) E_=E_+(V*(ibits(ci,0,1)+ibits(ci,1,1))*(ibits(cj,0,1)+ibits(cj,1,1))+0.25d0*DJ*(ibits(ci,0,1)-ibits(ci,1,1))*(ibits(cj,0,1)-ibits(cj,1,1)))*O_
				enddo
			enddo
			Ek(i1)=E_
			Ok(i1)=O_
		enddo
		!$omp end do
		if(allocated(k)) deallocate(k)
		!$omp end parallel
		Oq=1d0/sum(Ok*conjg(Ok))
		!!$omp parallel do if(.not.omp_in_parallel())
		do i=1,Ns
			do j=1,Ns
				O(i,j)=conjg(Ok(i))*Ok(j)*Oq
				E(i,j)=conjg(Ok(i))*Ek(j)*Oq
			enddo
		enddo
		!!$omp end parallel do
	end subroutine
end module

module vmc_main
	use model
	use lapack95, only : hegv, geev, hegv
	implicit none
contains
	subroutine set_cfg()
		integer :: c,n
		cfg=0
		cfg(:ne(1),1)=ibset(cfg(:ne(1),1),0)
		cfg(:ne(2),1)=ibset(cfg(:ne(2),1),1)
		call fisher_yates_shuffle(cfg(:,1))
		c=0
		do n=1,Ns
			select case(cfg(n,1))
			case(1)
				c=c+1
				cfg(n,2)=c
			case(2)
				c=c+1
				cfg(n,3)=c
			case(3)
				c=c+1
				cfg(n,2)=c
				c=c+1
				cfg(n,3)=c
			end select
		enddo
	end subroutine
	subroutine mc(wf,dwf)
		complex(8) :: wf(:,:)
		complex(8), optional :: dwf(:,:,:)
		complex(8) :: Ep
		real(8) :: Eb2,Eb
		complex(8) :: Sp(size(g),size(g))
		real(8) :: gp(size(g)),dscp,ddwp,afp,dscl,ddwl,afl
		complex(8) :: iA(sum(ne),sum(ne)),A(sum(ne),size(wf,2)),WA(size(wf,1),size(iA,2)),AW(size(iA,1),size(wf,2)),WAW(size(wf,1),size(wf,2)),El,pb,pbn,Sq_pml,Sq_zzl,Sq_pmp,Sq_zzp
		complex(8) :: Op(size(g)),Ol(size(g)),Sl(size(g),size(g)),SEp(size(g),size(g)),SEl(size(g),size(g)),Ok2l(size(Ok2,1),size(Ok2,2)),Ek2l(size(Ek2,1),size(Ek2,2)),Ok2p(size(Ok2,1),size(Ok2,2)),Ek2p(size(Ek2,1),size(Ek2,2))
		real(8) :: rd,gl(size(g))
		integer :: i,j,n,l,l1,l2,info,c,cfgl(size(cfg,1),size(cfg,2)),wcfg(size(A,1)),ac,sg
		integer, allocatable :: k(:),dcfg(:)
		logical :: is_update
		Ep=0d0; Eb=0d0; Eb2=0d0
		dscp=0d0; ddwp=0d0; afp=0d0; Sq_pmp=0d0; Sq_zzp=0d0
		Sp=0d0; gp=0d0; Op=0d0; SEp=0d0
		Ek2p=0d0; Ok2p=0d0
		cfgl=cfg
		c=0
		A=get_A(cfgl,wf,wcfg)
		iA=A(:,:sum(ne))
		call mat_inv(iA,info); if(info/=0) stop "err info"
		WA=matmul(wf(:,:sum(ne)),iA)
		AW=matmul(iA,wf(wcfg,:))
		WAW=matmul(wf(:,:sum(ne)),AW)
		c=0
		n=0
		ac=0
		is_update=.true.
		do
			c=c+1
			call two(cfgl,dcfg)
			!call two_all(cfgl,dcfg)
			call get_row(cfgl,dcfg,k,sg)

			if(mc_sg==3) then
				call get_Spb(cfgl,dcfg,pb,WA,iA,AW,WAW,wf)
				if(isnan(real(pb))) then
					write(*,*)pb
					stop
				endif
			else
				call get_pb(k,shape(0),pb,WA)
				pb=pb*conjg(pb)
			endif
			call random_number(rd)
			if(rd<real(pb)) then
				n=n+1
				do i=1,size(dcfg),2
					cfgl(mod(abs(dcfg(i:i+1))-1,Ns)+1,1)=ieor(cfgl(mod(abs(dcfg(i:i+1))-1,Ns)+1,1),(abs(dcfg(i))-1)/Ns+1)
					call swap(cfgl(mod(abs(dcfg(i+1))-1,Ns)+1,2+(dcfg(i)-1)/Ns),cfgl(mod(abs(dcfg(i))-1,Ns)+1,2+(abs(dcfg(i))-1)/Ns))
					wcfg(k(i))=k(i+1)
				enddo
				!call get_pb(k,shape(0),pb,WA)
				if(mc_sg==2.and.opt==1) then
					call update(k,WA,iA)
				elseif(mc_sg/=3) then
					call update(k,WA)
				else
					call update(k,WA,iA,AW,WAW,A,wf)
				endif
				!if(abs(pb)<1d-10.or.abs(pb)>5d2) then
				!if(mod(n,10000)==0) then
					!write(*,"(es12.4$)")sum(abs(matmul(A(:,:sum(ne)),iA)-diag(1d0,size(A,1)))),sum(abs(wf(:,:sum(ne))-matmul(WA,A(:,:sum(ne))))),sum(abs(matmul(A(:,:sum(ne)),AW)-wf(wcfg,:))),sum(abs(WAW-matmul(wf(:,:sum(ne)),AW)))
					!write(*,"(x)")
					!!if(any(abs(matmul(A(:,:sum(ne)),iA)-diag(1d0,size(A,1)))>2d-6)) then
						!iA=A(:,:sum(ne))
						!call mat_inv(iA,info)
						!WA=matmul(wf(:,:sum(ne)),iA)
						!AW=matmul(iA,wf(wcfg,:))
						!WAW=matmul(wf(:,:sum(ne)),AW)
					!!endif
				!endif
				is_update=.true.
			endif
			if(c>Nmc(1).and.mod(c-Nmc(1),Nmc(2))==0) then
				ac=ac+1
				if(is_update) then
					is_update=.false.
					select case(mc_sg)
					case(1)
						El=get_energy(cfgl,WA)
						!dscl=real(get_dsc(cfgl,WA))
						!Sq_pml=get_spin_pm(cfgl,WA,q)
						!Sq_zzl=get_spin_zz(cfgl,q)
					case(2)
						El=get_energy(cfgl,WA)
						call get_O(cfgl,WA,iA,dwf,Ol(:size(dwf,3)))
						do l1=1,size(Ol)
							do l2=1,size(Ol)
								Sl(l1,l2)=conjg(Ol(l1))*Ol(l2)
							enddo
						enddo
						gl=real(El*Ol)
					case(3)
						call get_overlap(cfgl,WA,iA,AW,WAW,wf,Ok2l,Ek2l)
					end select
				endif
				Ep=Ep+El
				Sq_pmp=Sq_pmp+Sq_pml
				Sq_zzp=Sq_zzp+Sq_zzl
				select case(mc_sg)
				case(1)
					dscp=dscp+dscl; ddwp=ddwp+ddwl; afp=afp+afl
				case(2)
					gp=gp+gl; Op=Op+Ol; Sp=Sp+Sl; SEp=SEp+SEl
				case(3)
					Ek2p=Ek2p+Ek2l
					Ok2p=Ok2p+Ok2l
				end select
				if(mod(ac,32)==0.and.mc_sg==3) then
					write(*,"(i7$)")ac
					write(*,"(es16.5$)")real(n)/c,sum(abs(Ek2p-conjg(transpose(Ek2p)))/ac)
					!write(*,"(i6,4es12.4)")ac,real(n)/c,Sq_pmp/ac,real(Ep/ac)
					write(*,"(x)")
					if(mod(ac,256)==0) then
						rewind(50)
						write(50,*)ac
						write(50,*)Ok2p/ac,Ek2p/ac
					endif
				endif
				!read(*,*)
				if(ac==Nmc(3)) then
					! check
					iA=get_A(cfgl,wf(:,:sum(ne)))
					pb=sum(abs(wf(:,:sum(ne))-matmul(WA,iA)))
					if(abs(pb)>1d-5) then
						!$omp critical
						write(*,*)"warn!!!!!",abs(pb)
						!$omp end critical
					endif

					Ep=Ep/Nmc(3)
					dscp=dscp/Nmc(3); ddwp=ddwp/Nmc(3); afp=afp/Nmc(3)
					Sq_pmp=Sq_pmp/Nmc(3); Sq_zzp=Sq_zzp/Nmc(3)
					Sp=Sp/Nmc(3); gp=gp/Nmc(3); Op=Op/Nmc(3); SEp=SEp/Nmc(3)
					do l1=1,size(Op)
						do l2=1,size(Op)
							!S(l1,l2)=2d0*SE(l1,l2)-O(l1)*g(l2)-O(l2)*g(l1)-2d0*E*S(l1,l2) ! maybe works
							!S(l1,l2)=2d0*(SE(l1,l2)-S(l1,l2)*E-O(l1)*g(l2)-O(l2)*g(l1))
							Sp(l1,l2)=Sp(l1,l2)-conjg(Op(l1))*Op(l2)
						enddo
					enddo
					gp=2d0*(gp-real(Ep*Op))
					Ek2p=Ek2p/Nmc(3)
					Ok2p=Ok2p/Nmc(3)

					!$omp critical
					cfg=cfgl
					select case(mc_sg)
					case(1:2)
						E=E+real(Ep)/n_omp; dsc=dsc+dscp/n_omp; ddw=ddw+ddwp/n_omp; af=af+afp/n_omp; S=S+Sp/n_omp; g=g+gp/n_omp
						!er=er+sqrt(abs((real(Eb2)-real(Ep)**2/(Nmc(3)/Nmc(4)-1))))/n_omp
						Sq_pm=Sq_pm+real(Sq_pmp)/n_omp
						Sq_zz=Sq_zz+real(Sq_zzp)/n_omp
					case(3)
						Ek2=Ek2+Ek2p/n_omp
						Ok2=Ok2+Ok2p/n_omp
					end select
					!$omp end critical
					exit
				endif
			endif
		enddo
	end subroutine
	subroutine vmc(sg,grad)
		integer :: sg
		real(8), optional :: grad(:)
		real(8) :: eg(size(g))
		complex(8) :: wf(Ns*spin,Ns*spin)
		complex(8) :: dwf(Ns*spin,Ns*spin,size(g))
		integer :: k,l,l1,l2,info
		mc_sg=sg
		select case(mc_sg)
		case(2)
			call ini_wf(wf,dwf)
			E=0d0; S=0d0; g=0d0; er=0d0
			!$omp parallel do if(n_omp>1)
			do k=1,n_omp
				call mc(wf,dwf)
			enddo
			!$omp end parallel do
			call heev(S,eg,'V')
			write(*,"(es12.4$)")var(1:)%val(1),E,er
			write(*,"(i3$)")int(sign(1d0,g))
			write(*,"(A$)")' | '
			write(20,"(es12.4$)")var(1:)%val(1),E,er,g
			eg=eg+abs(min(eg(1),0d0))+0.2d0
			do l=1,size(g)
				grad(l)=0d0
				do l1=1,size(g)
					do l2=1,size(g)
						if(abs(eg(l2)/eg(size(g)))<1d-3) then
							cycle
						endif
						grad(l)=grad(l)+real(S(l,l2)*conjg(S(l1,l2))*g(l1)/eg(l2))
					enddo
				enddo
			enddo
			write(*,"(i3$)")int(sign(1d0,grad))
			!grad=g
			grad=grad*Ns
			!stop
		case(1)
			call ini_wf(wf)
			E=0d0; er=0d0; dsc=0d0; ddw=0d0; af=0d0; Sq_pm=0d0; Sq_zz=0d0
			!$omp parallel do if(n_omp>1)
			do k=1,n_omp
				call mc(wf)
			enddo
			!$omp end parallel do
			write(*,"(es12.4$)")var(1:)%val(1),E,Sq_pm,Sq_zz,Sq_pm+Sq_zz,E/Sq_pm,dsc,er
			write(10,"(es12.4$)")var(1:)%val(1),E,er
		case(3)
			call ini_wf(wf)
			Ok2=0d0; Ek2=0d0
			!$omp parallel do if(n_omp>1)
			do k=1,n_omp
				call mc(wf)
			enddo
			!$omp end parallel do
			!write(*,"(es12.4$)")var(1:)%val(1),E,sum(Ek2),sum(Ok2),er
			rewind(50)
			write(50,*)-1
			write(50,*)Ok2,Ek2
		end select
	end subroutine
	subroutine variation()
		real(8) :: x(sum(var(1:)%n)),dx,El(150),er,grad(size(x)),pgrad(size(x))
		integer :: i
		!dx=0.03d0
		dx=0.1d0
		i=0
		if(size(x)==0) then
			return
		endif
		x=put(var(1:))
		do 
			i=i+1
			pgrad=grad
			call get(var(1:),x)
			write(*,"(i4$)")i
			call vmc(2,grad)
			El(i)=E
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
			x=x-grad*dx
			write(*,"(x)")
			write(20,"(x)")
		enddo
		write(*,*)"finished"
		!E=minval(El(:i))
	end subroutine
	subroutine spect(ut,gm,omg,m)
		integer :: ut,m
		real(8) :: gm,omg(:)
		complex(8) :: Sq(m),H_(size(Ek2,1),size(Ek2,2)),O_(size(Ek2,1),size(Ek2,2)),Om(size(Ek2,1),size(Ek2,2)),tmp(Ns)
		real(8) :: domg,EO(Ns),Eg(Ns)
		integer :: l,i,j,n0
		call ini_wf()
		domg=(omg(2)-omg(1))/m
		Sq=0d0
		!rewind(80)
		!read(80,*)Ok2,Ek2
		!write(*,*)sum(abs(Ek2-transpose(conjg(Ek2))))
		!call random_number(Ok2)
		!call random_number(Ek2)
		!Ok2=matmul(Ok2,transpose(conjg(Ok2)))
		O_=Ok2
		H_=Ek2
		!H_=0.5d0*(H_+transpose(conjg(H_)))
		call heev(O_,EO,'V')
		do i=1,Ns
			if(EO(i)>1d-9) then
				n0=i
				write(*,*)n0
				exit
			endif
		enddo
		!check truncate
		H_=matmul(transpose(conjg(O_)),matmul(H_,O_))
		write(*,*)maxval(abs(H_(:n0-1,:)))
		Om(n0:,n0:)=diag(EO(n0:))
		call hegv(H_(n0:,n0:),Om(n0:,n0:),Eg(n0:),jobz="V")
		psi0=matmul(transpose(conjg(O_)),psi0)
		tmp(n0:)=matmul(transpose(conjg(H_(n0:,n0:))),EO(n0:)*psi0(n0:))
		tmp(n0:)=conjg(tmp(n0:))*tmp(n0:)
		do l=1,m
			Sq(l)=Sq(l)+sum(tmp(n0:)/(omg(1)+domg*l-Eg(n0:)+img*gm))
		enddo
		write(*,*)-sum(imag(Sq)*domg)/pi
		write(*,*)sum(tmp(n0:))
		do l=1,m
			write(ut,"(es17.9$)")omg(1)+domg*l-E*Ns,Sq(l)
			write(ut,"(x)")
		enddo
		write(ut,"(x)")
	end subroutine
end module

program main
	use vmc_main
	implicit none
	logical :: f
	integer :: i,j,k
	!f=openfile(101,"../data/lattice.dat")
	f=openfile(10,"../data/grad.dat")
	f=openfile(20,"../data/var.dat")
	!f=openfile(30,"../data/alg.dat")
	!f=openfile(40,"../data/phyvar.dat")
	f=openfile(50,"../data/matrix_qusi.dat")
	f=openfile(70,"../data/spect.dat")

	call initial()

	Nmc=(/50000,Ns,1024*4,32/)
	!!20%
	!var(1:)%val(1)=(/-8.1583d-01,1.3409d-01,-3.4320d-01/)
	!E=-5.5474d-1
	!!4%
	!var(1:)%val(1)=(/-4.6190d-01,3.0090d-01,-2.6875d-01/)
	!E=-2.8843E-01
	!var(1:)%val(1)=(/-0.62d0,0.235d0,-0.315d0/)
	!E=-3.9796d-01
	var(1:)%val(1)=(/-0.93d0,0.145d0,-0.34d0/)
	ne(1)=Ns/2-10
	ne(2)=Ns-ne(1)
	ne=ne+1
	E=-5.5710d-01

	!!read(60,*)ne(1),var(1:)%val(1),E
	!rewind(50)
	!read(50,*)i
	!read(50,*)Ok2,Ek2
	!write(*,*)i
	!call spect(70,0.1d0,E*Ns+(/-10d0,10d0/),500)
	!Ek2=transpose(conjg(Ek2))
	!call spect(70,0.1d0,E*Ns+(/-10d0,10d0/),500)
	!stop
	do i=4,-4,-1
		!read(60,*)ne(1),var(1:)%val(1)
		!ne(1)=Ns/2-ne(1)
		!ne(1)=Ns/2-10
		!Nmc(1)=0
		!call set_cfg()
		!call vmc(1)
		q=(/pi,pi,0d0/)-real(i)/10d0*2d0*(/pi,pi,0d0/)
		!write(*,*)q

		call set_cfg()
		Nmc(3)=1024*128
		call vmc(3)
		call spect(70,0.1d0,E*Ns+(/-10d0,10d0/),500)
		Ek2=transpose(conjg(Ek2))
		call spect(70,0.1d0,E*Ns+(/-10d0,10d0/),500)
		write(*,"(x)")
		write(10,"(x)")
		!stop
	enddo
	stop
	ne(1)=Ns/2-20
	ne(2)=Ns-ne(1)
	ne=ne+1
	call set_cfg()
	call vmc(3)
	call spect(70,0.1d0,E*Ns+(/-10d0,10d0/),500)
	Ek2=transpose(conjg(Ek2))
	call spect(70,0.1d0,E*Ns+(/-10d0,10d0/),500)
	stop
	do i=2,2
		ne(1)=Ns/2-i
		ne(2)=Ns-ne(1)
		call set_cfg()
		call variation()
		write(*,"(es12.4$)")var(1:)%val(1),E,Sq_pm,Sq_zz,Sq_pm+Sq_zz,E/Sq_pm,er
		write(10,"(es12.4$)")var(1:)%val(1),E,er
		write(*,"(x)")
		write(10,"(x)")
	enddo
	stop
	do j=0,20
		!var(2)%val(1)=10d0**(-2d0+0.2*j)
		var(2)%val(1)=1d-2+j*0.04
		call vmc(2,grad)
		write(*,"(x)")
		write(10,"(x)")
	enddo
			
end program
