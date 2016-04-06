module global
	use M_hamilton_test
	use ifport, ifort_qsort => qsort
	implicit none
	real(8) :: t(2)=(/1d0,-0.3d0/)
	real(8), parameter :: DJ=0.3d0,V=0d0,U=0d0
	integer :: Nmc(4)
	integer :: n_omp=1
	integer, parameter :: opt=2 ! 1: Tr(AdA)
								! 2: O=cicj
	integer :: ne(2),vn
	integer, allocatable :: cfg(:),Ecfg(:)
	real(8) :: E,dsc,af,ddw,Sq_pm,Sq_zz,er
	real(8), allocatable :: g(:)
	complex(8), allocatable :: S(:,:)
	real(8), allocatable :: grad(:)
	complex(8), allocatable :: Ok2(:,:),Ek2(:,:)
	complex(8), allocatable :: psi0(:)
	integer :: Ns
	real(8) :: q(3)=(/pi,pi,0d0/)
	integer, allocatable :: nn(:,:)
	integer :: mc_sg=1 ! 1 static physical 
					   ! 2 energy and grad
					   ! 3 dynamic physical
	!$omp threadprivate(Ok2,Ek2,q,psi0,nn,cfg,Ecfg,ne,E,Sq_pm,Nmc,mc_sg)
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
		latt%c1=(/8d0,0d0,0d0/)
		latt%c2=(/0d0,2d0,0d0/)
		!latt%c1=latt%a1
		!latt%c2=latt%a2
		latt%T1=(/1d0,0d0,0d0/)*16
		latt%T2=(/0d0,1d0,0d0/)*16
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
		!var(iv(0))%val=1d-4

		!! ddw
		!call gen_var(sg=3,nb=1)
		!do i=1,size(var(iv(0))%bd)
			!var(iv(0))%bd(i)=ab(latt%sb(1,1)%nb(var(iv(0))%nb)%bd(i)%i(1))*dwave(i)*img
		!enddo
		!var(iv(0))%val=1d-3

		! sdw
		call gen_var(sg=4,nb=0)
		do i=1,size(var(iv(0))%bd)
			var(iv(0))%bd(i)=ab(i)*sin(sum(q*(latt%sb(1,1)%nb(var(iv(0))%nb)%bd(i)%r)))
		enddo
		var(iv(0))%val=8.88530d-1*2d0

		! on site cdw
		call gen_var(sg=3,nb=0)
		do i=1,size(var(iv(0))%bd)
			var(iv(0))%bd(i)=cos(sum(2d0*q*(latt%sb(1,1)%nb(var(iv(0))%nb)%bd(i)%r)))
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

		call var_shrink()

		Ns=latt%Ns
	end subroutine
	subroutine export_data(ut)
		integer :: ut,l,i
		do l=2,size(var(1:))
			do i=1,size(var(l)%bd)
				write(ut,"(es17.9$)")latt%sb(1,1)%nb(var(l)%nb)%bd(i)%r,&
					latt%sb(1,1)%nb(var(l)%nb)%bd(i)%dr,&
					var(l)%val(var(l)%bd2v(i)),&
					var(l)%bd(i)
				write(ut,"(x)")
			enddo
			write(ut,"(x/)")
		enddo
		stop
	end subroutine
	subroutine ini_wf(wf,dwf)
		complex(8), optional :: wf(:,:)
		complex(8), optional :: dwf(:,:,:)
		complex(8) :: H(Ns*spin,Ns*spin),cH(size(H,1),size(H,2)),D(size(H,1),size(H,2),sum(var(1:vn)%n)),DQ(size(H,1)-sum(ne),sum(ne)),Uk(size(H,1),size(H,2)),psi0_(Ns*Ns)
		real(8) :: E(size(H,1))
		integer :: l,i,j,n1,n2,nk(size(H,1),size(H,2)),nn_(0:Ns*Ns,2),kq(Ns)

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
				if(sum(abs(mod(brizon%ok(i,:2)/pi+brizon%ok(j,:2)/pi+(/2d0,2d0/)+1d-7,(/2d0,2d0/))-q(:2)/pi))<1d-5) then
				!if(sum(abs(mod(brizon%ok(i,:2)+brizon%ok(j,:2)+pi*(/2d0,2d0/),pi*(/2d0,2d0/))-q(:2)))<1d-7) then
					if(kq(i)/=0) then
						write(*,*)"warn! kq"
					else
						kq(i)=j
					endif
				endif
			enddo
			if(kq(i)==0) then
				write(*,*)q(:2)/pi,"warn! kq2"
				stop
			endif
		enddo
		nk=0
		do i=1,size(brizon%k,1)
			if(any(abs(var(:)%sg)==2).or.any(abs(var(:)%sg)==1)) then
				!call heev(H(i::size(brizon%k,1),i::size(brizon%k,1)),E(i::size(brizon%k,1)),"V")
				call mat_diag(H(i::size(brizon%k,1),i::size(brizon%k,1)),E(i::size(brizon%k,1)))
				nk(i::size(brizon%k,1),i::size(brizon%k,1))=1
				!write(30,"(es12.4$)")brizon%k(i,:2),E(i::size(brizon%k,1)),abs(H(i,i::size(brizon%k,1)))
				!write(30,"(x)")
			else
				call mat_diag(H(i:Ns:size(brizon%k,1),i:Ns:size(brizon%k,1)),E(i:Ns:size(brizon%k,1)))
				call mat_diag(H(Ns+i::size(brizon%k,1),Ns+i::size(brizon%k,1)),E(Ns+i::size(brizon%k,1)))
				nk(i:Ns:size(brizon%k,1),i:Ns:size(brizon%k,1))=1
				nk(Ns+i::size(brizon%k,1),Ns+i::size(brizon%k,1))=1
			endif
		enddo
		l=0
		Ecfg=0
		do i=1,size(E)
			if(E(i)<0d0) then
				l=l+1
				Ecfg(i)=l
			endif
		enddo
		if(count(Ecfg/=0)/=Ns) stop "Ecfg err"

		if(mc_sg==3) then
			psi0_=0d0
			l=0
			nn_=0
			do i=1,Ns
				do n1=1,size(E)
					o:	do n2=1,size(E)
						if(nk(i,n1)*nk(kq(i)+Ns,n2)>0.and.n1/=n2.and.Ecfg(n1)==0.and.Ecfg(n2)==0) then
							do j=1,l
								if(nn_(j,1)==n2.and.nn_(j,2)==n1) then
									psi0_(j)=psi0_(j)-conjg(H(i,n1))*conjg(H(kq(i)+Ns,n2))
									cycle o
								elseif(nn_(j,1)==n1.and.nn_(j,2)==n2) then
									psi0_(j)=psi0_(j)+conjg(H(i,n1))*conjg(H(kq(i)+Ns,n2))
									cycle o
								endif
							enddo
							l=l+1
							nn_(l,:)=(/n1,n2/)
							psi0_(l)=psi0_(l)+conjg(H(i,n1))*conjg(H(kq(i)+Ns,n2))
						endif
					enddo o
				enddo
			enddo
			!write(*,*)l
			!write(*,"(es12.4)")abs(psi0_(1:l))
			if(allocated(nn)) deallocate(nn,Ok2,Ek2,psi0)
			if(allocated(Ok2)) deallocate(Ok2)
			if(allocated(Ek2)) deallocate(Ek2)
			if(allocated(psi0)) deallocate(psi0)
			allocate(nn(0:l,2),Ok2(l,l),Ek2(l,l),psi0(l))
			nn=nn_(:l,:)
			Ecfg(nn(1,1))=Ns+1
			Ecfg(nn(1,2))=Ns+2
			nn(0,:)=nn(1,:)
			psi0=psi0_(:l)
		endif
		H=matmul(Uk,H)*sqrt(1d0/Ns)
		!write(*,"(2es12.4)")psi0
		!write(*,"(es12.4)")E(:sum(ne))
		!call Hamilton(var,cH)
		!write(*,*)sum(abs(matmul(transpose(conjg(H)),matmul(cH,H))-diag(E)))
		!stop

		!write(*,*)sum(abs(matmul(transpose(conjg(H)),matmul(wf,H))-diag(E)))

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
					!$omp parallel do if(.not.omp_in_parallel())
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
					!!$omp parallel do collapse(2) if(.not.omp_in_parallel()) !Ecfg is threadprivate
					do i=1,size(E)
						do j=1,size(E)
							if(abs(E(i)-E(j))<1d-8.or.Ecfg(i)/=0.or.Ecfg(j)==0) then
								D(i,j,l)=0d0
							else
								D(i,j,l)=D(i,j,l)/(E(j)-E(i))
							endif
						enddo
					enddo
					!!$omp end parallel do
					dwf(:,:,l)=matmul(matmul(H,D(:,:,l)),cH)
				end select
			enddo
		endif
		!write(*,*)"ini_wf finished"
	end subroutine
	real(8) function ab(i)
		integer :: i
		ab=(-1d0)**(mod(nint(sum(latt%sb(1,1)%nb(0)%bd(i)%r(:2))),2))
	end function
end module

module mc_utility
	use global
	use blas95, only : gemm, gerc
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
	subroutine change(cfg,icfg,dcfg,k)
		integer :: cfg(:),icfg(:)
		integer, allocatable :: dcfg(:),k(:)
		integer :: i,j,ii,jj
		if(allocated(dcfg)) deallocate(dcfg)
		if(allocated(k)) deallocate(k)
		call random_number(i,size(icfg))
		i=icfg(i)
		ii=i+(1-(i-1)/Ns*2)*Ns
		do 
			call random_number(j,Ns)
			j=j+(i-1)/Ns*Ns
			jj=j+(1-(j-1)/Ns*2)*Ns
			if(cfg(j)==0) then
				if(sign(1,-cfg(jj))==sign(1,-cfg(ii))) then
					allocate(dcfg(2))
					dcfg=(/i,-j/)
					allocate(k(2))
					k=(/cfg(dcfg(1)),abs(dcfg(2))/)
				else
					allocate(dcfg(4))
					dcfg=(/i,-j,ii,-jj/)
					allocate(k(4))
					k=(/cfg(dcfg(1)),abs(dcfg(2)),cfg(dcfg(3)),abs(dcfg(4))/)
				endif
				exit
			endif
		enddo
	end subroutine
	logical function get_row(cfg,dcfg,k,sg,P)
		integer :: cfg(:),dcfg(:),sg
		integer, allocatable :: k(:)
		integer, optional :: P(:)
		integer :: i,j,n,k_(10),cfg_(size(cfg)),dcfg_(size(dcfg))
		k_=0
		get_row=.true.
		cfg_(mod(abs(dcfg)-1,Ns)+1)=cfg(mod(abs(dcfg)-1,Ns)+1)
		cfg_(mod(abs(dcfg)-1,Ns)+1+Ns)=cfg(mod(abs(dcfg)-1,Ns)+1+Ns)
		dcfg_=0
		sg=1
	o:	do i=1,size(dcfg)+1
			if(present(P)) then
				if(any(P==i).or.i==size(dcfg)+1) then
					do j=1,size(dcfg)
						if(cfg_(mod(abs(dcfg(j))-1,Ns)+1)/=0.and.cfg_(mod(abs(dcfg(j))-1,Ns)+1+Ns)==0) then
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
				k_(n)=abs(cfg(dcfg_(i)))
			elseif(dcfg_(i)<0) then
				n=n+1
				k_(n)=abs(dcfg_(i))
			endif
		enddo
		if(allocated(k)) deallocate(k)
		allocate(k(n))
		k=k_(:n)
		!k(1)=raiseqq(sig$int)
	end function
	complex(8) function get_spin_pm(cfg,D,q)
		complex(8) :: D(:,:)
		integer :: cfg(:)
		real(8) :: q(3)
		complex(8) :: pb
		integer :: i,j,sg
		integer, allocatable :: k(:)
		get_spin_pm=0d0
		!$omp parallel do collapse(2) reduction(+:get_spin_pm) private(pb,k,sg) if(.not.omp_in_parallel())
		do i=1,Ns
			do j=1,Ns
				if(get_row(cfg,(/i,i+Ns,-j-Ns,-j/),k,sg,shape(0))) then
					call get_pb(k,shape(0),pb,D)
					get_spin_pm=get_spin_pm+pb*exp(img*sum(q*(latt%sb(1,1)%nb(0)%bd(i)%r-latt%sb(1,1)%nb(0)%bd(j)%r)))*sg
					if(allocated(k)) deallocate(k)
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
		!$omp parallel do collapse(2) reduction(+:get_spin_zz) if(.not.omp_in_parallel())
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
		integer, allocatable :: k(:)
		get_dsc=0d0
		!$omp parallel do collapse(2) reduction(+:get_dsc) private(pb,k,sg,i1,j1,i2,j2) if(.not.omp_in_parallel())
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
							call get_pb(k,shape(0),pb,D)
							get_dsc=get_dsc+pb*sg*dwave(n1)*dwave(n2)*latt%sb(1,1)%nb(1)%bd(n1)%bdc*latt%sb(1,1)%nb(1)%bd(n2)%bdc
						endif
						call swap(i2,j2)
					enddo
					call swap(i1,j1)
				enddo
				if(allocated(k)) deallocate(k)
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
		integer, allocatable :: k(:)
		get_energy=0d0
		!$omp parallel do collapse(2) reduction(+:get_energy) private(i,j,c,k,pb,sg) if(.not.omp_in_parallel())
		do l=1,size(t)
			do n=1,size(latt%sb(1,1)%nb(l)%bd)
				i=latt%sb(1,1)%nb(l)%bd(n)%i(1)
				j=latt%sb(1,1)%nb(l)%bd(n)%i(2)
				do p=1,2
					if(get_row(cfg,(/i,-j/),k,sg,shape(0))) then
						call get_pb(k,shape(0),pb,D)
						get_energy=get_energy-sg*t(l)*pb*latt%sb(1,1)%nb(l)%bd(n)%bdc!*-1
					endif
					if(get_row(cfg,(/-i-Ns,j+Ns/),k,sg,shape(0))) then
						call get_pb(k,shape(0),pb,D)
						get_energy=get_energy-sg*t(l)*pb*conjg(latt%sb(1,1)%nb(l)%bd(n)%bdc)!*-1
					endif
					if(l==1) then
						if(get_row(cfg,(/i,i+Ns,-j-Ns,-j/),k,sg,shape(0))) then
							call get_pb(k,shape(0),pb,D)
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
				if(allocated(k)) deallocate(k)
			enddo
		enddo
		!$omp end parallel do
		get_energy=get_energy/Ns
	end function
	subroutine get_O(icfg,iEcfg,WA,iA,dwf,O)
		complex(8) :: WA(:,:),iA(:,:),dwf(:,:,:),O(:)
		integer :: icfg(:),iEcfg(:)
		integer :: l,i
		O=0d0
		select case(opt)
		case(1)
			!$omp parallel do reduction(+:O) if(.not.omp_in_parallel())
			do i=1,size(iA,1)
				O=O+matmul(iA(i,:),dwf(icfg,iEcfg(i),:))
			enddo
			!$omp end parallel do
		case(2)
			!$omp parallel do reduction(+:O) if(.not.omp_in_parallel())
			do i=1,size(WA,1)
				O=O+matmul(WA(i,:),dwf(icfg,i,:))
			enddo
			!$omp end parallel do
		end select
	end subroutine
	subroutine get_Spb(cfg,dcfg,Ecfg,nn,pb,WA,iA,AW,WAW,wf)
		integer :: cfg(:),dcfg(:),Ecfg(:),nn(0:,:)
		complex(8) :: WA(:,:),AW(:,:),WAW(:,:),iA(:,:),wf(:,:)
		complex(8) :: pb
		complex(8) :: pb2(2),pbs(2)
		integer :: i,l,sg
		integer, allocatable :: k(:),m(:)
		!write(*,*)"get_Oq"
		if(.not.get_row(cfg,dcfg,k,sg,shape(0))) stop "err get_Spb"
		pbs=0d0
		!$omp parallel do private(pb2,sg) firstprivate(m) reduction(+:pbs) if(.not.omp_in_parallel())
		do i=1,ubound(nn,1)
			if(get_row(Ecfg,(/nn(0,1),nn(0,2),-nn(i,2),-nn(i,1)/),m,sg)) then
				call get_pb(shape(0),m,pb2(1),WA,iA,AW,WAW,wf)
				pbs(1)=pbs(1)+pb2(1)*conjg(pb2(1))
				call get_pb(k,m,pb2(2),WA,iA,AW,WAW,wf)
				pbs(2)=pbs(2)+pb2(2)*conjg(pb2(2))
			endif
			if(allocated(m)) deallocate(m)
		enddo
		!$omp end parallel do
		!write(*,*)real(pbs)
		pb=pbs(2)/pbs(1)
	end subroutine
	subroutine get_overlap(cfg,Ecfg,nn,WA,iA,AW,WAW,wf,O,E)
		integer :: cfg(:),Ecfg(:),nn(0:,:)
		complex(8) :: WA(:,:),AW(:,:),WAW(:,:),iA(:,:),wf(:,:),E(:,:),O(:,:)
		complex(8) :: Ek(size(E,1)),Ok(size(O,1)),pb,Oq,O_,E_
		real(8) :: iNs
		integer :: i1,j1,i,j,l,n,Si(0:1),c(2,2),p,sg,sgn
		integer, allocatable :: k(:),m(:)
		iNs=1d0/Ns
		!Ok=0d0
		!Ek=0d0
		Oq=1d0
		!$omp parallel do private(i,j,c,pb,O_,E_,sg,sgn) firstprivate(k,m) if(.not.omp_in_parallel())
		do i1=1,ubound(nn,1)
			if(get_row(Ecfg,(/nn(0,1),nn(0,2),-nn(i1,2),-nn(i1,1)/),m,sgn)) then
				call get_pb(shape(0),m,O_,WA,iA,AW,WAW,wf)
				O_=O_*sgn
				E_=0d0
				do l=1,size(t)
					do n=1,size(latt%sb(1,1)%nb(l)%bd)
						i=latt%sb(1,1)%nb(l)%bd(n)%i(1)
						j=latt%sb(1,1)%nb(l)%bd(n)%i(2)
						do p=1,2
							if(get_row(cfg,(/i,-j/),k,sg,shape(0))) then
								call get_pb(k,m,pb,WA,iA,AW,WAW,wf)
								E_=E_-pb*t(l)*latt%sb(1,1)%nb(l)%bd(n)%bdc*sg*sgn
							endif
							if(get_row(cfg,(/-i-Ns,j+Ns/),k,sg,shape(0))) then
								call get_pb(k,m,pb,WA,iA,AW,WAW,wf)
								E_=E_-pb*t(l)*conjg(latt%sb(1,1)%nb(l)%bd(n)%bdc)*sg*sgn
							endif
							if(l==1) then
								if(get_row(cfg,(/i,i+Ns,-j-Ns,-j/),k,sg,shape(0))) then
									call get_pb(k,m,pb,WA,iA,AW,WAW,wf)
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
				Ek(i1)=E_
				Ok(i1)=O_
			else
				Ek(i1)=0d0
				Ok(i1)=0d0
			endif
			if(allocated(k)) deallocate(k)
			if(allocated(m)) deallocate(m)
		enddo
		!$omp end parallel do
		Oq=1d0/sum(Ok*conjg(Ok))
		!$omp parallel do if(.not.omp_in_parallel())
		do i=1,size(O,1)
			do j=1,size(O,2)
				O(i,j)=conjg(Ok(i))*Ok(j)*Oq
				E(i,j)=conjg(Ok(i))*Ek(j)*Oq
			enddo
		enddo
		!$omp end parallel do
	end subroutine
end module

module vmc_main
	use model
	use lapack95, only : hegv, geev, hegv
	implicit none
contains
	subroutine set_cfg()
		integer :: i,n
		cfg=0
		cfg(Ns+1:Ns+ne(2))=(/(n,n=1,ne(2))/)
		call fisher_yates_shuffle(cfg(Ns+1:))
		n=ne(2)
		do i=1,Ns
			if(cfg(i+Ns)/=0) then
				n=n+1
				cfg(i)=n
			endif
			if(n==sum(ne)) exit
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
		integer :: i,j,n,l,l1,l2,info,c,cfgl(size(cfg)),Ecfgl(size(Ecfg)),icfg(size(A,1)),iEcfg(size(A,1)),ac,sg,icfgt(size(icfg)),iEcfgt(size(iEcfg))
		integer, allocatable :: k(:),m(:),dcfg(:)
		logical :: is_update
		Ep=0d0; Eb=0d0; Eb2=0d0
		dscp=0d0; ddwp=0d0; afp=0d0; Sq_pmp=0d0; Sq_zzp=0d0
		Sp=0d0; gp=0d0; Op=0d0; SEp=0d0
		Ek2p=0d0; Ok2p=0d0
		cfgl=cfg
		Ecfgl=Ecfg
		do i=1,size(cfgl)
			if(cfgl(i)/=0) icfg(cfgl(i))=i
			if(Ecfgl(i)/=0) iEcfg(Ecfgl(i))=i
		enddo
		A=wf(icfg,:)
		iA=A(:,iEcfg)
		call mat_inv(iA,info); if(info/=0) stop "err info"
		WA=matmul(wf(:,iEcfg),iA)
		AW=matmul(iA,wf(icfg,:))
		WAW=matmul(wf(:,iEcfg),AW)
		c=0
		n=0
		ac=0
		is_update=.true.
	lm:	do
			c=c+1
			call change(cfgl,icfg,dcfg,k)

			if(mc_sg==3) then
			!if(.false.) then
				call get_Spb(cfgl,dcfg,Ecfgl,nn,pb,WA,iA,AW,WAW,wf)
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
				
				cfgl(abs(dcfg(2::2)))=cfgl(abs(dcfg(2::2)))+cfgl(dcfg(1::2))
				cfgl(dcfg(1::2))=cfgl(abs(dcfg(2::2)))-cfgl(dcfg(1::2))
				cfgl(abs(dcfg(2::2)))=cfgl(abs(dcfg(2::2)))-cfgl(dcfg(1::2))
				icfg(k(1::2))=k(2::2)

				if(mc_sg==2.and.opt==1) then
					call update(k,WA,iA)
				elseif(mc_sg/=3) then
					call update(k,WA)
				else
					call get_pb(k,shape(0),pb,WA)
					if(abs(pb)<1d-10) then
						do i=1,ubound(nn,1)
							if(get_row(Ecfgl,(/nn(0,1),nn(0,2),-nn(i,2),-nn(i,1)/),m,sg)) then
								call get_pb(k,m,pb,WA,iA,AW,WAW,wf)
								if(abs(pb)>1d0) then
									Ecfgl(iEcfg(m(1::2)))=0
									Ecfgl(m(2::2))=m(1::2)
									iEcfg(m(1::2))=m(2::2)
									nn(0,:)=nn(i,:)
									exit
								endif
							endif
						enddo
						!do i=1,Ns
							!if(cfgl(i)/=0.and.cfgl(i+Ns)/=0) then
								!write(*,"(A$)")" ●"
							!elseif(cfgl(i)==0.and.cfgl(i+Ns)==0) then
								!write(*,"(A$)")" ○"
							!elseif(cfgl(i)==0.and.cfgl(i+Ns)/=0) then
								!write(*,"(A$)")"  "
							!endif
							!if(abs(latt%sb(1,1)%nb(0)%bd(i)%r(2)-latt%sb(1,1)%nb(0)%bd(i+1)%r(2))>1d-5) then
								!write(*,"(x)")
							!endif
						!enddo
						!write(*,"(x)")
						if(i==ubound(nn,1)+1) then
							cfgl(abs(dcfg(2::2)))=cfgl(abs(dcfg(2::2)))+cfgl(dcfg(1::2))
							cfgl(dcfg(1::2))=cfgl(abs(dcfg(2::2)))-cfgl(dcfg(1::2))
							cfgl(abs(dcfg(2::2)))=cfgl(abs(dcfg(2::2)))-cfgl(dcfg(1::2))
							do i=1,size(dcfg)
								if(cfgl(abs(dcfg(i)))/=0) icfg(cfgl(abs(dcfg(i))))=abs(dcfg(i))
							enddo
							write(*,*)"warn update, matrix is singularity, skip this configure"
							cycle lm
						else
							A=wf(icfg,:)
							iA=A(:,iEcfg)
							call mat_inv(iA,info)
							WA=matmul(wf(:,iEcfg),iA)
							AW=matmul(iA,wf(icfg,:))
							WAW=matmul(wf(:,iEcfg),AW)
						endif
					else
						call update(k,WA,iA,AW,WAW,A,wf)
					endif
					!if(mod(n,10000)==0) then
						!write(*,"(i7$)")n
						!write(*,"(es12.4$)")abs(pb),sum(abs(matmul(A(:,:sum(ne)),iA)-diag(1d0,size(A,1)))),sum(abs(wf(:,:sum(ne))-matmul(WA,A(:,:sum(ne))))),sum(abs(matmul(A(:,:sum(ne)),AW)-wf(icfg,:))),sum(abs(WAW-matmul(wf(:,:sum(ne)),AW)))
						!write(*,"(x)")
					!endif
				endif
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
						Sq_pml=get_spin_pm(cfgl,WA,q)
						!Sq_zzl=get_spin_zz(cfgl,q)
					case(2)
						El=get_energy(cfgl,WA)
						call get_O(icfg,iEcfg,WA,iA,dwf,Ol(:size(dwf,3)))
						do l1=1,size(Ol)
							do l2=1,size(Ol)
								Sl(l1,l2)=conjg(Ol(l1))*Ol(l2)
							enddo
						enddo
						gl=real(El*Ol)
					case(3)
						call get_overlap(cfgl,Ecfgl,nn,WA,iA,AW,WAW,wf,Ok2l,Ek2l)
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
				!$omp critical
				if(mod(ac,1024)==0.and.mc_sg==3) then
					write(*,"(i7$)")ac
					write(*,"(es16.5$)")real(n)/c,sum(abs(Ek2p-conjg(transpose(Ek2p)))/ac)!,sum(abs(matmul(A(:,iEcfg),iA)-diag(1d0,size(A,1)))),sum(abs(wf(:,:sum(ne))-matmul(WA,A(:,:sum(ne))))),sum(abs(matmul(A(:,:sum(ne)),AW)-wf(icfg,:))),sum(abs(WAW-matmul(wf(:,:sum(ne)),AW)))
					write(*,"(x)")
					if(mod(ac,10000)==0) then
						if(sum(abs(matmul(A(:,iEcfg),iA)-diag(1d0,size(A,1))))>1d-4) stop "update err"
					endif
					!rewind(50)
					!write(50,*)ac,size(psi0),Ns
					!write(50,*)E,psi0,Ok2p/ac,Ek2p/ac
				endif
				!$omp end critical
				!read(*,*)
				if(ac==Nmc(3)) then
					! check
					iA=wf(icfg,iEcfg)
					pb=sum(abs(wf(:,iEcfg)-matmul(WA,iA)))
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
		enddo lm
	end subroutine
	subroutine vmc(sg,grad)
		integer :: sg
		real(8), optional :: grad(:)
		real(8) :: eg(sum(var(1:)%n))
		complex(8) :: wf(Ns*spin,Ns*spin)
		complex(8) :: dwf(Ns*spin,Ns*spin,size(eg))
		integer :: k,l,l1,l2,info
		if(allocated(cfg)) deallocate(cfg)
		allocate(cfg(Ns*spin))
		if(allocated(Ecfg)) deallocate(Ecfg)
		allocate(Ecfg(Ns*spin))
		call set_cfg()
		mc_sg=sg
		select case(mc_sg)
		case(2)
			if(allocated(g)) deallocate(g)
			allocate(g(sum(var(1:)%n)))
			if(allocated(S)) deallocate(S)
			allocate(S(size(g),size(g)))
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
			!$omp critical
			write(*,"(es12.4$)")q(:2)/pi,var(1:)%val(1),E,Sq_pm,Sq_zz,Sq_pm+Sq_zz,E/Sq_pm,dsc,er
			write(10,"(es12.4$)")var(1:)%val(1),E,er
			write(*,"(x)")
			write(10,"(x)")
			!$omp end critical
		case(3)
			call ini_wf(wf)
			Ok2=0d0; Ek2=0d0
			!$omp parallel do if(n_omp>1)
			do k=1,n_omp
				call mc(wf)
			enddo
			!$omp end parallel do
			!write(*,"(es12.4$)")var(1:)%val(1),E,sum(Ek2),sum(Ok2),er
			!$omp critical
			rewind(50)
			write(50,*)-1,size(psi0),Ns
			write(50,*)E,psi0,Ok2,Ek2
			!$omp end critical
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
		complex(8) :: Sq(m),H_(size(Ek2,1),size(Ek2,2)),O_(size(Ek2,1),size(Ek2,2)),Om(size(Ek2,1),size(Ek2,2)),tmp(size(psi0)),psi0_(size(psi0))
		real(8) :: domg,EO(size(O_,1)),Eg(size(O_,1)),norm
		integer :: l,i,j,n0
		domg=(omg(2)-omg(1))/m
		Sq=0d0
		O_=Ok2
		H_=Ek2
		psi0_=psi0
		!H_=0.5d0*(H_+transpose(conjg(H_)))
		n0=1
		call heev(O_,EO,'V')
		!$omp critical
		do i=1,size(EO)
			if(EO(i)>1d-8) then
				n0=i
				write(*,*)n0,size(Ek2,1)
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
			Sq(l)=Sq(l)+sum(psi0_(n0:)/(omg(1)+domg*l-Eg(n0:)+img*gm))
		enddo
		norm=Sq_pm/real(sum(psi0_(n0:)))
		write(*,"(es12.4$)")-sum(imag(Sq)*domg)/pi,real(sum(psi0_(n0:))),Sq_pm
		write(ut,"('#q=',3es12.4,2i5)")q/pi,i,size(EO)
		do l=1,m
			write(ut,"(es17.9$)")omg(1)+domg*l-E*Ns,Sq(l)*norm
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
	real(8), allocatable :: ql(:,:)
	f=openfile(101,"../data/lattice.dat")
	f=openfile(10,"../data/grad.dat")
	f=openfile(20,"../data/var.dat")
	!f=openfile(30,"../data/alg.dat")
	!f=openfile(40,"../data/phyvar.dat")
	f=openfile(50,"../data/matrix_qusi.dat")
	f=openfile(70,"../data/spect.dat")


	!rewind(50)
	!read(50,*)i,j,Ns
	!allocate(psi0(j),Ok2(j,j),Ek2(j,j))
	!read(50,*)E,psi0,Ok2,Ek2
	!write(*,*)E,Ns,i,j
	!call spect(70,0.02d0,E*Ns+(/-10d0,10d0/),5000)
	!Ek2=transpose(conjg(Ek2))
	!call spect(70,0.02d0,E*Ns+(/-10d0,10d0/),5000)
	!stop

	call initial()

	Nmc=(/50000,Ns,1024,32/)
	!!20%
	!var(1:)%val(1)=(/-8.1583d-01,1.3409d-01,-3.4320d-01/)
	!E=-5.5474d-1
	!!4%
	!E=-2.8843E-01
	!var(1:)%val(1)=(/-0.62d0,0.235d0,-0.315d0/)
	!E=-3.9796d-01
	!ne(1)=Ns/2-10
	!ne(2)=Ns-ne(1)
	!!ne=ne+1
	!E=-5.5710d-01

	!var(1:)%val(1)=(/-1.15d0,0.074d0,-0.29d0/)
	!E=-6.5253E-01
	!ne(1)=Ns/2-14
	!ne(2)=Ns-ne(1)

	!var(1:)%val(1)=(/0d0,0.36d0,0.44d0/)
	!var(1:)%val(1)=(/0d0,0.36d0,0d0/)
	!var(1:)%val(1)=(/-0.39d0,0.32d0,-0.205d0/)
	!var(1:)%val(1)=(/-0.46d0,0.300d0,-0.268d0/)
	var(1:)%val(1)=(/-1.39752E+00,8.88530E-01*2E0,1.32678E-02,-2.00021E-01/)
	var(1:)%val(1)=(/-1.4051E+00,2.4353E+00,-1.9875E-01,-2.6827E-01/)
	!var(1:)%val(1)=(/-1.53531E+00,8.07142E-01*2E0,-6.06570E-02,-2.62224E-01,-2.15494E-01/)
	!var(1:)%val(1)=(/-9.4052E-01,2.1723E-01,-1.1267E-01/)
	!call variation()
	!stop
	!call export_data(101)

	ne(1)=Ns/2-16
	ne(2)=Ns-ne(1)
	Nmc(3)=1024*8*8
	call variation()
	stop

	!j=nint(sqrt(real(Ns)))/2
	!!$omp parallel do ordered copyin(Nmc) schedule(static,1)
	!do i=1,j*3
		!if(i<=j) then
			!q=(/0d0,0d0,0d0/)+((/pi,pi,0d0/)-(/0d0,0d0,0d0/))/j*mod(i-1,j)
		!elseif(i<=2*j) then
			!q=(/pi,pi,0d0/)+((/pi,0d0,0d0/)-(/pi,pi,0d0/))/j*mod(i-1,j)
		!else
			!q=(/pi,0d0,0d0/)+((/0d0,0d0,0d0/)-(/pi,0d0,0d0/))/j*mod(i-1,j)
		!endif

		ne(1)=Ns/2-10
		ne(2)=Ns-ne(1)
		Nmc(3)=1024*8
		call vmc(1)

		ne=ne+1
		Nmc(3)=1024*8*8
		call vmc(3)

		!!$omp ordered
		call spect(70,0.1d0,E*Ns+(/-10d0,10d0/),5000)
		Ek2=transpose(conjg(Ek2))
		call spect(70,0.1d0,E*Ns+(/-10d0,10d0/),5000)
		!!$omp end ordered
		!stop
	!enddo
	!!$omp end parallel do
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
