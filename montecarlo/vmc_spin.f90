module global
	use M_hamilton_m
	use ifport
	implicit none
	real(8) :: t(2)=(/1d0,-0.25d0/)
	real(8), parameter :: DJ=1d0/3d0,V=0d0,U=0d0
	integer :: Nmc(4)
	integer :: n_omp=1
	integer, parameter :: opt=2 ! 1: Tr(AdA)
								! 2: O=cicj
	integer :: ne(2),vn
	integer, allocatable :: cfg(:,:)
	real(8) :: E,dsc,af,ddw,Sq_pm,Sq_zz,er
	real(8), allocatable :: g(:)
	complex(8), allocatable :: S(:,:)
	real(8), allocatable :: grad(:)
	complex(8), allocatable :: Ok2(:,:),Ek2(:,:)
	real(8), allocatable :: E2(:)
	integer :: map(0:3),Ns
	real(8) :: q(3)=(/pi,pi,0d0/)
	integer :: mc_sg=1 ! 1 static physical 
					   ! 2 energy and grad
					   ! 3 dynamic physical
contains
	subroutine initial()
		integer :: i,l
		allocate(var(-1000:1000))
		call init_random_seed()
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
		!brizon_o%k(:,2)=brizon_o%k(:,2)+pi/10d0
		!call check_lattice(101)
		write(*,*)"Total site number is: ",latt%Ns

		! cp
		call gen_var(sg=1,nb=0)
		var(iv(0))%bd=-1d0
		var(iv(0))%val=0d0

		! dsc
		call gen_var(sg=2,nb=1)
		do i=1,size(var(iv(0))%bd)
			var(iv(0))%bd(i)=dwave(i)
		enddo
		var(iv(0))%val=1d-3

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
			var(iv(0))%val=0d0
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
		allocate(g(sum(var(1:)%n)))
		allocate(S(sum(var(1:)%n),sum(var(1:)%n)))
		allocate(Ok2(Ns,Ns),Ek2(Ns,Ns),E2(Ns))
		map=(/b"10",b"11",b"00",b"01"/)
	end subroutine
	subroutine ini_wf(wf,dwf)
		complex(8) :: wf(:,:)
		complex(8), optional :: dwf(:,:,:)
		complex(8) :: H(latt%Ns*spin,latt%Ns*spin),cH(size(H,1),size(H,2)),D(size(H,1),size(H,2),sum(var(1:vn)%n)),Q(size(H,1)-sum(ne),sum(ne))
		real(8) :: E(size(H,1))
		integer :: l,i,j
		call Hamilton(var,H)
		if(any(abs(var(:)%sg)==1.or.abs(var(:)%sg)==2)) then
			call mat_diag(H,E)
		else
			call mat_diag(H(:Ns,:Ns),E(:Ns))
			call mat_diag(H(Ns+1:,Ns+1:),E(Ns+1:))
			!write(*,*)0.5d0*(E(ne(1))+E(ne(1)+1))
			do i=1,ne(2)
				call swap(H(:,ne(1)+i),H(:,Ns+i))
				call swap(E(ne(1)+i),E(Ns+i))
			enddo
		endif
		wf=H(:,:sum(ne))
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
								Q(i-sum(ne),j)=0d0
								cycle
							endif
							Q(i-sum(ne),j)=D(i,j,l)/(E(j)-E(i))
						enddo
					enddo
					!$omp end parallel do
					!$omp parallel do
					do i=1,size(dwf,1)
						dwf(i,:,l)=matmul(matmul(H(i,sum(ne)+1:),Q),cH(:sum(ne),:))
					enddo
					!$omp end parallel do
				end select
			enddo
		endif
	end subroutine
	function ab(i)
		integer :: i
		real(8) :: ab
		ab=(-1d0)**(mod(nint(sum(latt%sb(1)%nb(0)%bd(i)%r(:2))),2))
	end function
end module

module mc_utility
	use global
	implicit none
contains
	subroutine get_pb(k,m,D,pb)
		complex(8) :: D(:,:),pb
		integer, allocatable :: k(:),m(:)
		complex(8) :: A(size(k,1),size(k,1))
		integer :: i,j,n
		if(.not.allocated(k)) then
			pb=0d0
			return
		endif
		n=size(k)
		do i=1,n
			do j=1,n
				A(i,j)=D(m(i),k(j))
			enddo
		enddo
		select case(n)
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
	subroutine update(k,m,D)
		complex(8) :: D(:,:)
		integer :: k(:),m(:)
		complex(8) :: tmp(size(D,1),size(D,2)),ipb(2,2),pb,tmp2(2)
		integer :: i,j,n
		select case(size(k))
		case(0)
			return
		case(1)
			tmp=D
			pb=1d0/D(m(1),k(1))
			!$omp parallel do private(tmp2) if(.not.omp_in_parallel())
			do j=1,size(D,2)
				if(j==k(1)) then
					tmp2(1)=pb*(tmp(m(1),j)-1d0)
				else
					tmp2(1)=pb*tmp(m(1),j)
				endif
				do i=1,size(D,1)
					D(i,j)=tmp(i,j)-tmp(i,k(1))*tmp2(1)
				enddo
			enddo
			!$omp end parallel do
		case(2)
			tmp=D
			pb=1d0/(D(m(1),k(1))*D(m(2),k(2))-D(m(1),k(2))*D(m(2),k(1)))
			ipb(1,1)=pb*D(m(2),k(2))
			ipb(1,2)=-pb*D(m(1),k(2))
			ipb(2,1)=-pb*D(m(2),k(1))
			ipb(2,2)=pb*D(m(1),k(1))
			!$omp parallel do private(tmp2) if(.not.omp_in_parallel())
			do j=1,size(D,2)
				if(k(1)==j) then
					tmp2=matmul(ipb,(/tmp(m(1),j)-1d0,tmp(m(2),j)/))
				elseif(k(2)==j) then
					tmp2=matmul(ipb,(/tmp(m(1),j),tmp(m(2),j)-1d0/))
				else
					tmp2=matmul(ipb,(/tmp(m(1),j),tmp(m(2),j)/))
				endif
				do i=1,size(D,1)
					D(i,j)=tmp(i,j)-sum((/tmp(i,k(1)),tmp(i,k(2))/)*tmp2)
				enddo
			enddo
			!$omp end parallel do
		case default
			write(*,*)"err"
			stop
		end select
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
	subroutine get_km(cfg,dcfg,k,m,sg)
		integer :: cfg(:,:),dcfg(:),sg
		integer, allocatable :: k(:),m(:)
		integer :: i,l,j1,j2,j,n,km(5,2),c(size(cfg,1),3),dcfg_(size(dcfg))
		km=0
		sg=1
		if(allocated(k)) deallocate(k)
		if(allocated(m)) deallocate(m)
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
				n=n+1
				km(n,1)=c(j1,2+l)
				km(n,2)=j2+Ns*l
				call swap(c(j2,2+l),c(j1,2+l))
			else
				return
			endif
			do j=1,n-1
				if((j1+Ns*l)==km(j,2)) then
					km(j,2)=km(n,2)
					n=n-1
					exit
				endif
			enddo
		enddo
		if(any(map(c(mod(abs(dcfg)-1,Ns)+1,1))==3)) then
			return
		endif
		allocate(k(n),m(n))
		k=km(:n,1)
		m=km(:n,2)
		!k(1)=raiseqq(sig$int)
	end subroutine
	function get_spin_kk(cfg,D,q)
		complex(8) :: D(:,:)
		integer :: cfg(:,:)
		real(8) :: q(3)
		complex(8) :: get_spin_kk
		complex(8) :: pb
		integer :: i,j,sg,n
		integer, allocatable :: k(:),m(:)
		real(8) :: k1(3),k2(3)
		k1=brizon_o%k(1,:)
		k2=brizon_o%k(9,:)
		get_spin_kk=0d0
		!do n=1,size(latt%sb(1)%nb(1)%bd)
			!i=latt%sb(1)%nb(1)%bd(n)%i(1)
			!j=latt%sb(1)%nb(1)%bd(n)%i(2)
		do i=1,Ns
			do j=1,Ns
			call get_km(cfg,(/i,-j/),k,m,sg)
			call get_pb(k,m,D,pb)
			get_spin_kk=get_spin_kk+pb*sg*exp(img*sum(k1*latt%i2r(i,:)-k2*latt%i2r(j,:)))!*latt%sb(1)%nb(1)%bd(n)%bdc
		enddo
		enddo
	end function
	function get_spin_pm(cfg,D,q)
		complex(8) :: D(:,:)
		integer :: cfg(:,:)
		real(8) :: q(3)
		complex(8) :: get_spin_pm
		complex(8) :: pb
		integer :: i,j,sg
		integer, allocatable :: k(:),m(:)
		get_spin_pm=0d0
		!$omp parallel do collapse(2) reduction(+:get_spin_pm) private(pb,k,m,sg) if(.not.omp_in_parallel())
		do i=1,Ns
			do j=1,Ns
				call get_km(cfg,(/j,j+Ns,-i-Ns,-i/),k,m,sg)
				call get_pb(k,m,D,pb)
				get_spin_pm=get_spin_pm+pb*exp(img*sum(q*(latt%i2r(i,:)-latt%i2r(j,:))))*sg
			enddo
		enddo
		!$omp end parallel do
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
	function get_energy(cfg,D)
		complex(8) :: D(:,:)
		integer :: cfg(:,:)
		complex(8) :: get_energy
		complex(8) :: pb
		integer :: i,j,n,l,ci,cj,p,sg
		integer, allocatable :: k(:),m(:)
		get_energy=0d0
		!$omp parallel do collapse(2) reduction(+:get_energy) private(i,j,ci,cj,m,k,pb,sg) if(.not.omp_in_parallel())
		do l=1,size(t)
			do n=1,size(latt%sb(1)%nb(l)%bd)
				i=latt%sb(1)%nb(l)%bd(n)%i(1)
				j=latt%sb(1)%nb(l)%bd(n)%i(2)
				do p=1,2
					call get_km(cfg,(/i,-j/),k,m,sg)
					call get_pb(k,m,D,pb)
					get_energy=get_energy-t(l)*pb*sg*latt%sb(1)%nb(l)%bd(n)%bdc!*-1
					call get_km(cfg,(/-i-Ns,j+Ns/),k,m,sg)
					call get_pb(k,m,D,pb)
					get_energy=get_energy-t(l)*pb*sg*latt%sb(1)%nb(l)%bd(n)%bdc!*-1
					if(l==1) then
						call get_km(cfg,(/i,i+Ns,-j-Ns,-j/),k,m,sg)
						call get_pb(k,m,D,pb)
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
		!$omp end parallel do
		get_energy=get_energy/Ns
	end function
	function get_A(cfg,wf)
		integer :: cfg(:,:)
		complex(8) :: wf(:,:)
		complex(8) :: get_A(sum(ne),sum(ne))
		integer :: n,l
		do l=0,1
			do n=1,Ns
				if(cfg(n,2+l)>0) then
					get_A(cfg(n,2+l),:)=wf(n+Ns*l,:)
				endif
			enddo
		enddo
	end function
	subroutine get_Spb(cfg,dcfg,D,pb)
		integer :: cfg(:,:),dcfg(:)
		complex(8) :: D(:,:)
		complex(8) :: pb
		complex(8) :: pb2(2),pbs(2)
		integer :: i,l,Si(0:1),sg
		integer, allocatable :: k(:),m(:)
		!write(*,*)"get_Oq"
		do l=0,1
			do i=1,Ns
				if((.not.btest(cfg(i,1),l)).and.cfg(i,2+l)>0) then
					Si(l)=i
					exit
				endif
			enddo
		enddo
		cfg(Si,1)=ibset(cfg(Si,1),(/0,1/))
		pbs=0d0
		!$omp parallel do private(pb2,k,m,sg) reduction(+:pbs) if(.not.omp_in_parallel())
		do i=1,size(cfg,1)
			call get_km(cfg,(/Si(0),Si(1)+Ns,-i-Ns,-i/),k,m,sg)
			call get_pb(k,m,D,pb2(1))
			pb2(1)=pb2(1)*sg
			call get_km(cfg,(/Si(0),Si(1)+Ns,dcfg,-i-Ns,-i/),k,m,sg)
			call get_pb(k,m,D,pb2(2))
			pb2(2)=pb2(2)*sg
			pbs=pbs+exp(img*sum(q*latt%i2r(i,:)))*pb2
		enddo
		!$omp end parallel do
		cfg(Si,1)=ibclr(cfg(Si,1),(/0,1/))
		pb=pbs(2)/pbs(1)
	end subroutine
	subroutine get_overlap(cfg,D,O,E)
		integer :: cfg(:,:)
		complex(8) :: E(:,:),D(:,:),O(:,:)
		complex(8) :: Ek(size(E,1)),Ok(size(O,1)),pb,Oq,O_,E_
		real(8) :: iNs
		integer :: cfg_s(size(cfg,1)),i1,j1,i,j,l,n,Si(0:1),ci,cj,p,sg
		integer, allocatable :: k(:),m(:)
		iNs=1d0/Ns
		Ok=0d0
		Ek=0d0
		do l=0,1
			do i=1,Ns
				if((.not.btest(cfg(i,1),l)).and.cfg(i,2+l)>0) then
					Si(l)=i
					exit
				endif
			enddo
		enddo
		cfg_s=cfg(:,1)
		cfg(Si,1)=ibset(cfg(Si,1),(/0,1/))
		!$omp parallel do collapse(2) private(i,j,ci,cj,pb,O_,E_,k,m,sg) reduction(+:Ok,Ek) if(.not.omp_in_parallel())
		do j1=1,Ns
			do i1=1,Ns
				!j1=i1
				call get_km(cfg,(/Si(0),Si(1)+Ns,-i1-Ns,-j1/),k,m,sg)
				call get_pb(k,m,D,O_)
				O_=O_*sg
				E_=0d0
				do l=1,size(t)
					do n=1,size(latt%sb(1)%nb(l)%bd)
						i=latt%sb(1)%nb(l)%bd(n)%i(1)
						j=latt%sb(1)%nb(l)%bd(n)%i(2)
						do p=1,2
							call get_km(cfg,(/Si(0),Si(1)+Ns,i,-j,-i1-Ns,-j1/),k,m,sg)
							call get_pb(k,m,D,pb)
							E_=E_-pb*t(l)*latt%sb(1)%nb(l)%bd(n)%bdc*sg
							call get_km(cfg,(/Si(0),Si(1)+Ns,-i-Ns,j+Ns,-i1-Ns,-j1/),k,m,sg)
							call get_pb(k,m,D,pb)
							E_=E_-pb*t(l)*latt%sb(1)%nb(l)%bd(n)%bdc*sg
							if(l==1) then
								call get_km(cfg,(/Si(0),Si(1)+Ns,i,i+Ns,-j-Ns,-j,-i1-Ns,-j1/),k,m,sg)
								call get_pb(k,m,D,pb)
								E_=E_+pb*0.5d0*DJ*sg
							endif
							j=latt%sb(1)%nb(l)%bd(n)%i(1)
							i=latt%sb(1)%nb(l)%bd(n)%i(2)
						enddo
						ci=map(cfg_s(i))
						cj=map(cfg_s(j))
						if(l==1) E_=E_+(V*(ibits(ci,0,1)+ibits(ci,1,1))*(ibits(cj,0,1)+ibits(cj,1,1))+0.25d0*DJ*(ibits(ci,0,1)-ibits(ci,1,1))*(ibits(cj,0,1)-ibits(cj,1,1)))*O_
					enddo
				enddo
				E_=E_*iNs
				O_=O_*iNs
				do n=1,Ns
					Ek(n)=Ek(n)+exp(img*sum(brizon_o%k(n,:)*(latt%i2r(i1,:)-latt%i2r(j1,:))+q*latt%i2r(i1,:)))*E_
					Ok(n)=Ok(n)+exp(img*sum(brizon_o%k(n,:)*(latt%i2r(i1,:)-latt%i2r(j1,:))+q*latt%i2r(i1,:)))*O_
				enddo
			enddo
		enddo
		!$omp end parallel do
		Oq=sum(Ok)
		Oq=1d0/(Oq*conjg(Oq))
		!!$omp parallel do if(.not.omp_in_parallel())
		do i=1,Ns
			do j=1,Ns
				O(i,j)=conjg(Ok(i))*Ok(j)*Oq
				E(i,j)=conjg(Ok(i))*Ek(j)*Oq
			enddo
		enddo
		!!$omp end parallel do
		cfg(Si,1)=ibclr(cfg(Si,1),(/0,1/))
	end subroutine
end module

module vmc_main
	use model
	use lapack95, only : hegv
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
		complex(8) :: iA(sum(ne),sum(ne)),D(size(wf,1),size(wf,2)),El,pb,A(sum(ne),sum(ne)),pbn,Sq_pml,Sq_zzl,Sq_pmp,Sq_zzp
		complex(8) :: Op(size(g)),Ol(size(g)),Sl(size(g),size(g)),SEp(size(g),size(g)),SEl(size(g),size(g)),Ok2l(size(Ok2,1),size(Ok2,2)),Ek2l(size(Ek2,1),size(Ek2,2)),Ok2p(size(Ok2,1),size(Ok2,2)),Ek2p(size(Ek2,1),size(Ek2,2))
		real(8) :: rd,gl(size(g))
		integer :: i,j,n,l,l1,l2,info,c,cfgl(size(cfg,1),size(cfg,2)),ac,sg
		integer, allocatable :: k(:),m(:),dcfg(:)
		logical :: is_update
		Ep=0d0; Eb=0d0; Eb2=0d0
		dscp=0d0; ddwp=0d0; afp=0d0; Sq_pmp=0d0; Sq_zzp=0d0
		Sp=0d0; gp=0d0; Op=0d0; SEp=0d0
		Ek2p=0d0; Ok2p=0d0
		cfgl=cfg
		c=0
		iA=get_A(cfgl,wf)
		call mat_inv(iA,info); if(info/=0) stop "err info"
		D=matmul(wf,iA)
		c=0
		n=0
		ac=0
		is_update=.true.
		do
			c=c+1
			call two(cfgl,dcfg)
			!call two_all(cfgl,dcfg)
			call get_km(cfgl,dcfg,k,m,sg)

			if(mc_sg==3) then
				call get_Spb(cfgl,dcfg,D,pb)
			else
				call get_pb(k,m,D,pb)
			endif
			call random_number(rd)
			if(rd<real(pb*conjg(pb))) then
				n=n+1
				do i=1,size(dcfg),2
					cfgl(mod(abs(dcfg(i:i+1))-1,Ns)+1,1)=ieor(cfgl(mod(abs(dcfg(i:i+1))-1,Ns)+1,1),(abs(dcfg(i))-1)/Ns+1)
					call swap(cfgl(mod(abs(dcfg(i+1))-1,Ns)+1,2+(dcfg(i)-1)/Ns),cfgl(mod(abs(dcfg(i))-1,Ns)+1,2+(abs(dcfg(i))-1)/Ns))
				enddo
				call get_pb(k,m,D,pb)
				call update(k,m,D)
				if(abs(pb)<1d-10.or.abs(pb)>5d2) then
					A=get_A(cfgl,wf)
					if(any(abs(matmul(D,A)-wf)>2d-6)) then
						!write(*,*)"warn!",n
						call mat_inv(A,info)
						D=matmul(wf,A)
					endif
				endif
				is_update=.true.
			endif
			if(c>Nmc(1).and.mod(c-Nmc(1),Nmc(2))==0) then
				ac=ac+1
				if(is_update) then
					is_update=.false.
					select case(mc_sg)
					case(1)
						El=get_energy(cfgl,D)
						!Sq_pml=get_spin_pm(cfgl,D,q)
						!Sq_zzl=get_spin_zz(cfgl,q)
						!Sq_pml=get_spin_kk(cfgl,D,q)
					case(3)
						call get_overlap(cfgl,D,Ok2l,Ek2l)
					end select
				endif
				Ep=Ep+El
				Sq_pmp=Sq_pmp+Sq_pml
				Sq_zzp=Sq_zzp+Sq_zzl
				if(mc_sg==3) then
					Ek2p=Ek2p+Ek2l
					Ok2p=Ok2p+Ok2l
				endif
				if(mod(ac,8)==0) then
					write(*,"(i9,es12.4$)")ac,real(n)/c!,real(sum(Ek2p))/ac,sum(abs(Ek2p-conjg(transpose(Ek2p))))/ac
					!write(*,"(i6,4es12.4)")ac,real(n)/c,Sq_pmp/ac,real(Ep/ac)
					write(*,"(x)")
					if(mod(ac,256)==0) then
						rewind(50)
						write(50,*)ac
						write(50,*)Ok2p/ac,Ek2p/ac
					endif
				endif
			endif
			if(c>=(Nmc(1)+Nmc(2)*Nmc(3))) then
				! check
				A=get_A(cfgl,wf)
				call mat_inv(A,info)
				if(sum(abs(D-matmul(wf,A)))>1d-5) then
					!$omp critical
					write(*,*)"warn!!!!!",sum(abs(D-matmul(wf,A)))
					!$omp end critical
				endif

				Ek2p=Ek2p/Nmc(3)
				Ok2p=Ok2p/Nmc(3)
				Ep=Ep/Nmc(3); Sq_pmp=Sq_pmp/Nmc(3); Sq_zzp=Sq_zzp/Nmc(3)

				!$omp critical
				cfg=cfgl
				E=E+real(Ep)/n_omp
				Sq_pm=Sq_pm+real(Sq_pmp)/n_omp
				Sq_zz=Sq_zz+real(Sq_zzp)/n_omp
				Ek2=Ek2+Ek2p/n_omp
				Ok2=Ok2+Ok2p/n_omp
				!$omp end critical
				exit
			endif
		enddo
	end subroutine
	subroutine vmc(sg,grad)
		integer :: sg
		real(8), optional :: grad(:)
		real(8) :: eg(size(g))
		complex(8) :: wf(Ns*spin,sum(ne))
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
			call mat_diag(S,eg)
			write(*,"(es12.4$)")var(1:)%val(1),E,er,g
			write(10,"(es12.4$)")var(1:)%val(1),E,er,g
			eg=eg+abs(min(eg(1),0d0))+1d-2
			do l=1,size(g)
				grad(l)=0d0
				do l1=1,size(g)
					do l2=1,size(g)
						grad(l)=grad(l)+real(S(l,l2)*conjg(S(l1,l2))*g(l1)/eg(l2))
					enddo
				enddo
			enddo
			!grad=g
		case(1)
			call ini_wf(wf)
			E=0d0; er=0d0; dsc=0d0; ddw=0d0; af=0d0; Sq_pm=0d0; Sq_zz=0d0
			!$omp parallel do if(n_omp>1)
			do k=1,n_omp
				call mc(wf)
			enddo
			!$omp end parallel do
			write(*,"(es12.4$)")var(1:)%val(1),E,Sq_pm,Sq_zz,Sq_pm+Sq_zz,E/Sq_pm,er
			write(10,"(es12.4$)")var(1:)%val(1),E,er
		case(3)
			do k=1,size(cfg)
				if(cfg(k,1)==3) then
					cfg(k,1)=0
					exit
				endif
			enddo
			call ini_wf(wf)
			Ok2=0d0; Ek2=0d0
			!$omp parallel do if(n_omp>1)
			do k=1,n_omp
				call mc(wf)
			enddo
			!$omp end parallel do
			!write(*,"(es12.4$)")var(1:)%val(1),E,sum(Ek2),sum(Ok2),er
			rewind(50)
			write(50,*)"finial"
			write(50,*)Ok2,Ek2
		end select
	end subroutine
	subroutine spect(ut,gm,omg,m)
		integer :: ut,m
		real(8) :: gm,omg(:)
		complex(8) :: Sq(m),Ok2_(Ns,Ns),Ek2_(Ns,Ns),tmp(Ns)
		real(8) :: domg
		integer :: l,i,n0
		domg=(omg(2)-omg(1))/m
		Sq=0d0
		!rewind(80)
		!read(80,*)Ok2,Ek2
		write(*,*)sum(Ek2),sum(Ok2)
		write(*,*)sum(abs(Ek2-transpose(conjg(Ek2))))
		!stop
		!call random_number(Ok2)
		!call random_number(Ek2)
		!Ek2=0.5d0*(Ek2+transpose(conjg(Ek2)))
		Ek2_=Ek2
		!Ok2=matmul(Ok2,transpose(conjg(Ok2)))
		Ok2_=Ok2
		call heev(Ok2_,E2,'V')
		do i=1,Ns
			if(E2(i)>1d-9) then
				n0=i
				exit
			endif
		enddo
		Ok2_(:,n0:)=matmul(Ok2_(:,n0:),diag(1d0/sqrt(E2(n0:))))
		Ek2=matmul(transpose(conjg(Ok2_)),matmul(Ek2,Ok2_))
		!write(*,*)maxval(abs(Ek2(:n0-1,:)))
		!stop
		call heev(Ek2(n0:,n0:),E2(n0:),'V')
		Ek2(:n0-1,:)=0d0
		Ek2(:,:n0-1)=0d0
		Ek2(:n0-1,:n0-1)=diag(1d0,n0-1)
		Ek2=matmul(Ok2_,Ek2)
		!write(*,*)sum(abs(matmul(Ek2_,Ek2(:,n0:))-matmul(matmul(Ok2,Ek2(:,n0:)),diag(E2(n0:)))))
		!write(*,"(es12.4)")E2(n0:)
		!stop
		!write(*,*)sum(abs(matmul(transpose(conjg(Ek2(:,n0:))),matmul(Ok2,Ek2(:,n0:)))-diag(1d0,Ns-n0+1)))
		write(*,*)sum(abs(matmul(matmul(conjg(Ek2(:,n0:)),transpose(Ek2(:,n0:))),transpose(Ok2(:,:)))-diag(1d0,Ns)))
		tmp=0d0
		do i=1,size(Ek2,1)
			tmp(n0:)=tmp(n0:)+Ek2(i,n0:)*sum(Ok2(i,:))
		enddo
		do l=1,m
			Sq(l)=Sq(l)+sum(tmp(n0:)*conjg(tmp(n0:))/(omg(1)+domg*l-E2(n0:)+img*gm))
		enddo
		write(*,*)-sum(imag(Sq)*domg)/pi
		do l=1,m
			write(ut,"(es17.9$)")omg(1)+domg*l,Sq(l)
			write(ut,"(x)")
		enddo
	end subroutine
end module

program main
	use vmc_main
	implicit none
	logical :: f
	integer :: i,j,k
	!f=openfile(101,"../data/lattice.dat")
	f=openfile(10,"../data/grad.dat")
	!f=openfile(20,"../data/var.dat")
	!f=openfile(30,"../data/alg.dat")
	!f=openfile(40,"../data/phyvar.dat")
	f=openfile(50,"../data/matrix.dat")
	open(unit=60,file="../data/input_vmc.dat")
	open(unit=70,file="../data/spect.dat")
	open(unit=80,file="../data/matrix_read.dat")

	call initial()

	!Nmc=(/1,1,100000/)
	!Nmc=(/10000,Ns,100000/)
	!Nmc=(/10000,Ns,16384,2**i/)
	!Nmc=(/10000,Ns,16384,32/)
	!Nmc=(/10000,Ns*50,1024*16,32/)
	!Nmc=(/100000,Ns,1024,32/)
	!Nmc=(/5000000,Ns,1024*1024,32/)
	!Nmc=(/500000,Ns,1024*1024,32/)
	Nmc=(/500000,Ns,1024*1024,32/)
	!Nmc=(/10000,Ns,1000,32/)
	!var(1:)%val(1)=(/-9.3d-1,1.45d-1,-3.4d-1/)

	read(50,*)i
	read(50,*)Ok2,Ek2
	write(*,*)i
	write(*,*)sum(abs(Ek2-transpose(conjg(Ek2))))
	!stop
	call spect(70,0.1d0,(/-50d0,-30d0/),500)
	stop
	do i=0,20
		read(60,*)ne(1),var(1:)%val(1)
		ne(1)=Ns/2-ne(1)
		!var(2)%val=10d0**(-2d0+i/4d0)
		ne(1)=Ns/2-5
		ne(2)=Ns-ne(1)
		call set_cfg()
		call vmc(3)
		!call spect(70,0.1d0,(/E2(1),E2(size(E2))/),100)
		write(*,"(x)")
		write(10,"(x)")
		stop
	enddo
	stop
	!do 
		!read(*,*)var(1:)%val(1)
		!call vmc(1)
	!enddo
	do i=-10,10
		var(1)%val(1)=0.1d0*i
		do j=1,20
			var(2)%val(1)=0.01d0*j
			do k=-8,8
				var(3)%val(1)=-0.02d0*k
				call vmc(1)
				write(*,"(x)")
				write(10,"(x)")
			enddo
		enddo
	enddo
			
end program
