module global
	use M_hamilton_m
	implicit none
	real(8) :: t(2)=(/1d0,-0.25d0/)
	real(8) :: DJ=1d0/3d0,V=-1d0/3d0/4d0,U=0d0
	integer :: Nmc(4)
	integer :: n_omp=1
	integer, parameter :: opt=2 ! 1: Tr(AdA)
								! 2: O=cicj
	integer :: ne(2),vn
	integer, allocatable :: cfg(:,:)
	real(8) :: E,dsc,af,ddw,er
	real(8), allocatable :: g(:)
	complex(8), allocatable :: S(:,:)
	real(8), allocatable :: grad(:)
	include 'nlopt.f'
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
		latt%T1=(/1d0,0d0,0d0/)*16
		latt%T2=(/0d0,1d0,0d0/)*16
		latt%bdc=(/1d0,-1d0,0d0/)
		allocate(latt%sub(1,3))
		latt%sub(1,:)=(/0d0,0d0,0d0/)
		call latt%gen_latt()
		call latt%gen_bond(size(t))
		brizon%n1=1
		brizon%n2=1
		call latt%gen_brizon(brizon)
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
		var(iv(0))%val=1d-01

		! ssc
		call gen_var(sg=-2,nb=0)
		do i=1,size(var(iv(0))%bd)
			var(iv(0))%bd(i)=1d0
		enddo
		var(iv(0))%val=1d-3

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

		do l=1,size(t)
			call gen_var(sg=-3,nb=l)
			var(iv(0))%bd=-1d0
			var(iv(0))%val=t(l)
		enddo

		vn=iv(1)

		call var_shrink()

		allocate(cfg(latt%Ns,3))
		allocate(grad(sum(var(1:)%n)))
		allocate(g(sum(var(1:)%n)))
		allocate(S(sum(var(1:)%n),sum(var(1:)%n)))
	end subroutine
	subroutine ini_wf(wf,dwf)
		complex(8) :: wf(:,:)
		complex(8), optional :: dwf(:,:,:)
		complex(8) :: H(latt%Ns*spin,latt%Ns*spin),cH(size(H,1),size(H,2)),D(size(H,1),size(H,2),sum(var(1:vn)%n)),Q(size(H,1)-sum(ne),sum(ne))
		real(8) :: E(size(H,1))
		integer :: l,i,j
		call var%Hamilton(H)
		if(any(abs(var(:)%sg)==1.or.abs(var(:)%sg)==2)) then
			call mat_diag(H,E)
		else
			call mat_diag(H(:latt%Ns,:latt%Ns),E(:latt%Ns))
			call mat_diag(H(latt%Ns+1:,latt%Ns+1:),E(latt%Ns+1:))
			!write(*,*)0.5d0*(E(ne(1))+E(ne(1)+1))
			do i=1,ne(2)
				call swap(H(:,ne(1)+i),H(:,latt%Ns+i))
				call swap(E(ne(1)+i),E(latt%Ns+i))
			enddo
		endif
		wf=H(:,:sum(ne))
		if(present(dwf)) then
			cH=transpose(conjg(H))
			dwf=0d0
			call var(1:vn)%dHamilton(H,cH,D)

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
		integer :: k(:),m(:)
		if(k(2)==0) then
			pb=D(m(1),k(1))
		else
			pb=(D(m(1),k(1))*D(m(2),k(2))-D(m(1),k(2))*D(m(2),k(1)))
		endif
	end subroutine
	subroutine update(k,m,D,iA)
		complex(8) :: D(:,:)
		complex(8), optional :: iA(:,:)
		integer :: k(:),m(:)
		complex(8) :: tmp(size(D,1),size(D,2)),ipb(2,2),pb,tmp2(2),tmpA(latt%Ns,latt%Ns)
		integer :: i,j
		tmp=D
		if(present(iA)) then
			tmpA=iA
		endif
		if(k(2)==0) then
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
					if(present(iA)) then
						if(i<=size(iA,1)) then
							iA(i,j)=tmpA(i,j)-tmpA(i,k(1))*tmp2(1)
						endif
					endif
				enddo
			enddo
			!$omp end parallel do
		else
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
					if(present(iA)) then
						if(i<=size(iA,1)) then
							iA(i,j)=tmpA(i,j)-sum((/tmpA(i,k(1)),tmpA(i,k(2))/)*tmp2)
						endif
					endif
				enddo
			enddo
			!$omp end parallel do
		endif
	end subroutine
end module

module model
	use mc_utility
	implicit none
contains
	function get_af(cfg)
		integer :: cfg(:,:)
		real(8) :: get_af
		integer :: n
		!$omp parallel do reduction(+:get_af) if(.not.omp_in_parallel())
		do n=1,size(latt%sb(1)%nb(0)%bd)
			select case(cfg(n,1))
			case(1:2)
				get_af=get_af+ab(n)
			case(3)
				get_af=get_af+2d0*ab(n)
			end select
			get_af=get_af-ab(n)
		enddo
		!$omp end parallel do
		get_af=get_af*0.5d0/latt%Ns
	end function
	function get_ddw(cfg,D)
		complex(8) :: D(:,:)
		integer :: cfg(:,:)
		complex(8) :: get_ddw
		complex(8) :: pb
		integer :: n,i,j,k(2),m(2)
		get_ddw=0d0
		!$omp parallel do private(i,j,k,m,pb) reduction(+:get_ddw) if(.not.omp_in_parallel())
		do n=1,size(latt%sb(1)%nb(1)%bd)
			i=latt%sb(1)%nb(1)%bd(n)%i(1)
			j=latt%sb(1)%nb(1)%bd(n)%i(2)
			if(cfg(i,1)<cfg(j,1)) call swap(i,j)
			k=0; m=0
			select case(ieor(cfg(i,1),cfg(j,1)))
			case(1)
				k(1)=cfg(i,2)
				m(1)=j
			case(2)
				k(1)=cfg(i,3)
				m(1)=j+latt%Ns
			end select
			call get_pb(k,m,D,pb)
			get_ddw=get_ddw+pb*dwave(n)*ab(i)*(1.5d0-ieor(cfg(i,1),cfg(j,1)))
		enddo
		!$omp end parallel do
		get_ddw=get_ddw/latt%Ns
	end function
	function get_dsc(cfg,D)
		complex(8) :: D(:,:)
		integer :: cfg(:,:)
		complex(8) :: get_dsc
		complex(8) :: pb
		integer :: i1,i2,j1,j2,n1,n2,n,k(2),m(2)
		get_dsc=0d0
		!$omp parallel do private(i1,i2,j1,j2,k,m,pb) reduction(+:get_dsc) if(.not.omp_in_parallel())
		do n1=1,size(latt%sb(1)%nb(1)%bd)
			i1=latt%sb(1)%nb(1)%bd(n1)%i(1)
			i2=latt%sb(1)%nb(1)%bd(n1)%i(2)
			do n2=1,size(latt%sb(1)%nb(1)%bd)
				j1=latt%sb(1)%nb(1)%bd(n2)%i(1)
				j2=latt%sb(1)%nb(1)%bd(n2)%i(2)
				if(cfg(j1,1)==2.and.cfg(j2,1)==2.and.ieor(cfg(i1,1),cfg(i2,1))==3) then
					do n=1,2
						k=0; m=0
						select case(cfg(i1,1))
						case(3)
							k(1)=cfg(i1,2)
							m(1)=j1
							k(2)=cfg(j2,3)
							m(2)=i2+latt%Ns
						case(0)
							k(1)=cfg(j1,3)
							m(1)=i1+latt%Ns
							k(2)=cfg(i2,2)
							m(2)=j2
						end select
						call get_pb(k,m,D,pb)
						get_dsc=get_dsc+dwave(n1)*dwave(n2)*pb
						call swap(j1,j2)
					enddo
				endif
			enddo
		enddo
		!$omp end parallel do
		get_dsc=get_dsc*0.25d0/latt%Ns**2
	end function
	function get_energy(cfg,D)
		complex(8) :: D(:,:)
		integer :: cfg(:,:)
		complex(8) :: get_energy
		complex(8) :: pb
		integer :: i,j,n,l,k(2),m(2)
		get_energy=0d0
		do l=1,size(t)
			!$omp parallel do private(i,j,k,m,pb) reduction(+:get_energy) if(.not.omp_in_parallel())
			do n=1,size(latt%sb(1)%nb(l)%bd)
				i=latt%sb(1)%nb(l)%bd(n)%i(1)
				j=latt%sb(1)%nb(l)%bd(n)%i(2)
				if(cfg(i,1)<cfg(j,1)) call swap(i,j)
				k=0; m=0
				select case(ieor(cfg(i,1),cfg(j,1)))
				case(1)
					k(1)=cfg(i,2)
					m(1)=j
					call get_pb(k,m,D,pb)
					get_energy=get_energy-t(l)*conjg(pb)*latt%sb(1)%nb(1)%bd(n)%bdc
				case(2)
					k(1)=cfg(i,3)
					m(1)=j+latt%Ns
					call get_pb(k,m,D,pb)
					get_energy=get_energy+t(l)*conjg(pb)*latt%sb(1)%nb(1)%bd(n)%bdc
				case(3)
					if(l==1) then
						k(1)=cfg(i,2)
						m(1)=j
						k(2)=cfg(i,3)
						m(2)=j+latt%Ns
						call get_pb(k,m,D,pb)
						get_energy=get_energy+0.5d0*DJ*conjg(pb)-0.25d0*DJ+V
					endif
				case(0)
					if(l==1.and.cfg(i,1)/=2) then
						get_energy=get_energy+0.25d0*DJ+V
					endif
				end select
			enddo
			!$omp end parallel do
		enddo
		get_energy=get_energy/latt%Ns
	end function
	subroutine get_O(cfg,D,iA,dwf,O)
		complex(8) :: D(:,:),iA(:,:),dwf(:,:,:),O(:)
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
					O=O+matmul(D(:,cfg(i,2)),dwf(i,:,:))
				endif
				if(btest(cfg(i,1),1)) then
					O=O+matmul(D(:,cfg(i,3)),dwf(i+Ns,:,:))
				endif
			enddo
			!$omp end parallel do
		end select
	end subroutine
	subroutine two(cfg,i,j)
		integer :: cfg(:,:),i,j
		call random_number(i,size(cfg,1))
		do 
			call random_number(j,size(cfg,1))
			if(ieor(cfg(i,1),cfg(j,1))/=0) then
				exit
			endif
		enddo
	end subroutine
	function get_A(cfg,wf)
		integer :: cfg(:,:)
		complex(8) :: wf(:,:)
		complex(8) :: get_A(sum(ne),sum(ne))
		integer :: n,Ns
		Ns=latt%Ns
		do n=1,Ns
			select case(cfg(n,1))
			case(1)
				get_A(cfg(n,2),:)=wf(n,:)
			case(2)
				get_A(cfg(n,3),:)=wf(n+Ns,:)
			case(3)
				get_A(cfg(n,2),:)=wf(n,:)
				get_A(cfg(n,3),:)=wf(n+Ns,:)
			end select
		enddo
	end function
end module

module vmc_main
	use model
	implicit none
contains
	subroutine set_cfg()
		integer :: c,n
		cfg=0
		cfg(:ne(1),1)=ibset(cfg(:ne(1),1),0)
		cfg(:ne(2),1)=ibset(cfg(:ne(2),1),1)
		call fisher_yates_shuffle(cfg(:,1))
		c=0
		do n=1,latt%Ns
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
	subroutine mc(sg,wf,dwf)
		integer :: sg
		complex(8) :: wf(:,:)
		complex(8), optional :: dwf(:,:,:)
		complex(8) :: Ep
		real(8) :: Eb2,Eb
		complex(8) :: Sp(size(g),size(g))
		real(8) :: gp(size(g)),dscp,ddwp,afp,dscl,ddwl,afl
		complex(8) :: iA(sum(ne),sum(ne)),D(size(wf,1),size(wf,2)),El,pb,A(sum(ne),sum(ne))
		complex(8) :: Op(size(g)),Ol(size(g)),Sl(size(g),size(g)),SEp(size(g),size(g)),SEl(size(g),size(g))
		real(8) :: rd,gl(size(g))
		integer :: i,j,n,l,l1,l2,Ns,info,k(2),m(2),ncfg(2,size(cfg,2)),c,cfgl(size(cfg,1),size(cfg,2))
		logical :: is_update
		Ns=latt%Ns
		Ep=0d0; Eb=0d0; Eb2=0d0
		dscp=0d0; ddwp=0d0; afp=0d0
		Sp=0d0; gp=0d0; Op=0d0; SEp=0d0
		cfgl=cfg
		!call fisher_yates_shuffle(cfgl)
		iA=get_A(cfgl,wf)
		call mat_inv(iA,info); if(info/=0) stop "err info"
		D=matmul(wf,iA)
		c=0
		is_update=.true.
		do
			c=c+1
			call two(cfgl,i,j)
			!call random_number(n,size(latt%sb(1)%nb(1)%bd))
			!i=latt%sb(1)%nb(1)%bd(n)%i(1)
			!j=latt%sb(1)%nb(1)%bd(n)%i(2)
			if(cfgl(i,1)<cfgl(j,1)) call swap(i,j)
			ncfg(1,:)=cfgl(i,:)
			ncfg(2,:)=cfgl(j,:)
			k=0; m=0
			select case(ieor(cfgl(i,1),cfgl(j,1)))
			case(1)
				ncfg(1,1)=ibclr(ncfg(1,1),0)
				ncfg(2,1)=ibset(ncfg(2,1),0)
				call swap(ncfg(1,2),ncfg(2,2))
				k(1)=cfgl(i,2)
				m(1)=j
				call get_pb(k,m,D,pb)
			case(2)
				ncfg(1,1)=ibclr(ncfg(1,1),1)
				ncfg(2,1)=ibset(ncfg(2,1),1)
				call swap(ncfg(1,3),ncfg(2,3))
				k(1)=cfgl(i,3)
				m(1)=j+Ns
				call get_pb(k,m,D,pb)
			case(3)
				call swap(ncfg(1,1:3),ncfg(2,1:3))
				k(1)=cfgl(i,2)
				m(1)=j
				k(2)=cfgl(i,3)
				m(2)=j+Ns
				call get_pb(k,m,D,pb)
			case(0)
				pb=0d0
			end select

			call random_number(rd)
			if(rd<real(pb*conjg(pb))) then
				if(sg==2.and.opt==1) then
					call update(k,m,D,iA)
				else
					call update(k,m,D)
				endif
				cfgl(i,:)=ncfg(1,:)
				cfgl(j,:)=ncfg(2,:)
				is_update=.true.
			endif
			if(c>Nmc(1).and.mod(c-Nmc(1),Nmc(2))==0) then
				if(is_update) then
					is_update=.false.
					El=get_energy(cfgl,D)
					if(sg==3) then
						dscl=real(get_dsc(cfgl,D))
						ddwl=real(get_ddw(cfgl,D))
						afl=real(get_af(cfgl))
					elseif(sg==2) then
						call get_O(cfgl,D,iA,dwf,Ol(:size(dwf,3)))
						do l1=1,size(Ol)
							do l2=1,size(Ol)
								Sl(l1,l2)=conjg(Ol(l1))*Ol(l2)
								SEl(l1,l2)=conjg(Ol(l1))*Ol(l2)*El
							enddo
						enddo
						gl=real(El*Ol)
					endif
				endif
				Ep=Ep+El
				Eb=Eb+real(El)/Nmc(4)
				if(mod(c-Nmc(1),Nmc(2)*Nmc(4))==0) then
					Eb2=Eb2+Eb**2/(Nmc(3)/Nmc(4))/(Nmc(3)/Nmc(4)-1)
					Eb=0d0
				endif
				if(sg==3) then
					dscp=dscp+dscl; ddwp=ddwp+ddwl; afp=afp+afl
				elseif(sg==2) then
					gp=gp+gl; Op=Op+Ol; Sp=Sp+Sl; SEp=SEp+SEl
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

				Ep=Ep/Nmc(3)
				dscp=dscp/Nmc(3); ddwp=ddwp/Nmc(3); afp=afp/Nmc(3)
				Sp=Sp/Nmc(3); gp=gp/Nmc(3); Op=Op/Nmc(3); SEp=SEp/Nmc(3)
				do l1=1,size(Op)
					do l2=1,size(Op)
						!S(l1,l2)=2d0*SE(l1,l2)-O(l1)*g(l2)-O(l2)*g(l1)-2d0*E*S(l1,l2) ! maybe works
						!S(l1,l2)=2d0*(SE(l1,l2)-S(l1,l2)*E-O(l1)*g(l2)-O(l2)*g(l1))
						Sp(l1,l2)=Sp(l1,l2)-conjg(Op(l1))*Op(l2)
					enddo
				enddo
				gp=2d0*(gp-real(Ep*Op))

				!$omp critical
				cfg=cfgl
				E=E+real(Ep)/n_omp; dsc=dsc+dscp/n_omp; ddw=ddw+ddwp/n_omp; af=af+afp/n_omp; S=S+Sp/n_omp; g=g+gp/n_omp
				er=er+sqrt(abs((real(Eb2)-real(Ep)**2/(Nmc(3)/Nmc(4)-1))))/n_omp
				!$omp end critical
				exit
			endif
		enddo
	end subroutine
	subroutine variation()
		real(8) :: x(sum(var(1:)%n)),dx,El(200),er,grad(size(x)),pgrad(size(x))
		integer :: i
		dx=5d0
		i=0
		if(size(x)==0) then
			return
		endif
		x=var(1:)%put()
		do 
			i=i+1
			pgrad=grad
			call var(1:)%get(x)
			write(*,"(i4$)")i
			call vmc(2,grad)
			El(i)=E
			er=er*1.5d0
			if((El(i)-er)>El(max(i-1,1))) then
				grad=pgrad
				x=x+grad*dx
				dx=dx*0.8d0
				write(*,"(' err',es10.2$)")dx,(El(i)-er),El(max(i-1,1))
				i=i-1
			endif
			if(all(abs(grad*dx/x)<1d-3).or.dx<1d-3.or.i==size(El)) then
				if(minval(El(:i))<(El(i)-er)) then
					write(*,*)"var err, min E is ", minval(El(:i))
				endif
				exit
			endif
			x=x-grad*dx
			write(*,"(x)")
		enddo
		write(*,*)"finished"
		E=minval(El(:i))
	end subroutine
	subroutine vmc(sg,grad)
		integer :: sg
		real(8), optional :: grad(:)
		real(8) :: eg(size(g))
		complex(8) :: wf(latt%Ns*spin,sum(ne))
		complex(8) :: dwf(latt%Ns*spin,latt%Ns*spin,size(g))
		integer :: k,l,l1,l2
		if(sg==2) then
			call ini_wf(wf,dwf)
			E=0d0; S=0d0; g=0d0; er=0d0
			!$omp parallel do if(n_omp>1)
			do k=1,n_omp
				call mc(sg,wf,dwf)
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
		else
			call ini_wf(wf)
			E=0d0; er=0d0; dsc=0d0; ddw=0d0; af=0d0
			!$omp parallel do if(n_omp>1)
			do k=1,n_omp
				call mc(sg,wf)
			enddo
			!$omp end parallel do
			write(*,"(es12.4$)")var(1:)%val(1),E,er
			write(10,"(es12.4$)")var(1:)%val(1),E,er
		endif
	end subroutine
	subroutine variation_nlopt()
		real(8) :: x(sum(var(1:)%n)),minf
		real(8) :: lb(size(x)),ub(size(x))
		real(8) :: cst(1)
		integer(8) opt,opt1
		integer ires
		lb=-1d0
		ub=1d0
		x=0d0

		!if(size(x)>1) then
			!lb(1)=-2d0
			!ub(1)=0d0
			!x(1)=-1d0
		!endif

		!call nlo_create(opt, NLOPT_GN_DIRECT, size(x))
		!call nlo_create(opt, NLOPT_G_MLSL_LDS, size(x))
		!call nlo_create(opt, NLOPT_GD_STOGO, size(x)) ! dosen't works
		call nlo_create(opt, NLOPT_GN_CRS2_LM, size(x)) ! best
		!call nlo_create(opt, NLOPT_GN_ESCH, size(x))
		!call nlo_create(opt, NLOPT_GN_ISRES, size(x))

		!!call nlo_create(opt, NLOPT_LN_BOBYQA, size(x))
		!call nlo_create(opt, NLOPT_LN_NELDERMEAD, size(x))
		!call nlo_create(opt, NLOPT_LN_SBPLX, size(x))

		!call nlo_create(opt, NLOPT_LD_LBFGS, size(x))
		!call nlo_create(opt, NLOPT_LD_TNEWTON_PRECOND_RESTART, size(x))
		!call nlo_create(opt, NLOPT_LD_MMA, size(x))
		!call nlo_create(opt, NLOPT_LD_SLSQP, size(x))
		!call nlo_create(opt, NLOPT_LD_VAR1, size(x))
		!call nlo_create(opt, NLOPT_LD_CCSAQ, size(x))

		call nlo_set_lower_bounds(ires, opt, lb)
		call nlo_set_upper_bounds(ires, opt, ub)
		call nlo_set_min_objective(ires, opt, nlopt_fn, cst)
		!call nlo_set_precond_min_objective(ires, opt, nlopt_fn, nlopt_precond , cst)
		!call nlo_set_xtol_rel(ires, opt, 1d-2)
		call nlo_set_maxeval(ires, opt, min(size(x)*50,200))
		!call nlo_set_local_optimizer(ires, opt, opt1);

		call nlo_optimize(ires, opt, x, minf)
		if(ires<0) then
			write(*,"(A$)")"ires is negetive"
			write(*,"(i3$)")ires
			write(*,"(es12.4)")cst
		endif
		call nlo_destroy(opt)
		E=minf
		write(*,"(i3$)")ires
		write(*,"(es12.4$)")minf,x,cst
		write(*,"(x)")
		call var(1:)%get(x)
	end subroutine
	subroutine nlopt_precond(n, x, v, vpre, cst)
		integer :: n
		real(8) :: x(n),v(n),vpre(n)
		complex(8) :: cst(:),hess(n,n)
		hess=reshape(cst(3:),(/n,n/))
		vpre=real(matmul(hess,v))
	end subroutine
	subroutine nlopt_fn(val, n, x, grad, need_gradient,cst)
		integer :: n
		real(8) :: val
		real(8), optional :: x(n),grad(n)
		real(8), optional :: cst(:)
		integer, optional :: need_gradient
		real(8) :: er
		call var(1:)%get(x)
		if(need_gradient/=0) then
			call vmc(2,grad)
		else
			call vmc(1)
		endif
		write(*,"(es12.4$)")x
		write(*,"(x)")
	end subroutine
end module

program main
	use vmc_main
	implicit none
	logical :: f
	integer :: i,j
	!f=openfile(101,"../data/lattice.dat")
	f=openfile(10,"../data/grad.dat")
	!f=openfile(20,"../data/var.dat")
	!f=openfile(30,"../data/alg.dat")
	!f=openfile(40,"../data/phyvar.dat")

	call initial()

	ne(1)=latt%Ns/2-20
	ne(2)=latt%Ns-ne(1)
	call set_cfg()
	!Nmc=(/1,1,100000/)
	!Nmc=(/10000,latt%Ns,100000/)
	!Nmc=(/10000,latt%Ns,16384,2**i/)
	!Nmc=(/10000,latt%Ns,16384,32/)
	!Nmc=(/10000,latt%Ns*50,1024*16,32/)
	Nmc=(/10000,latt%Ns,1024*4,32/)

	!do i=-10,10
		!!var(2)%val(1)=10d0**(-2d0+i/10d0)
		!var(3)%val(1)=1.8d-01+i/50d0
		call vmc(1,grad)
		write(*,"(x)")
		write(10,"(x)")
	!enddo
	stop
	!var(1:)%val(1)=(/-9.0500d-01, 1.2000d-01, 1.8000d-01/)

	do i=20,20,1
		ne(1)=latt%Ns/2-i
		ne(2)=latt%Ns-ne(1)
		call set_cfg()
		Nmc=(/10000,latt%Ns,1024*4,32/)

		call variation()
		write(30,*)E
		Nmc=(/10000,sum(ne),1024,1/)
		call vmc(3)
		write(*,"(es12.4$)")2d0*ne(1)/latt%Ns,E,er,dsc,af
		write(*,"(x)")
		write(40,"(es12.4$)")2d0*ne(1)/latt%Ns,E,er,var(1:)%put(),dsc,af
		write(40,"(x)")
	enddo
	stop
end program
