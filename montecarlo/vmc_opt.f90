module global
	use M_const
	use M_utility
	use M_hamilton
	implicit none
	type t_phy
		complex(8) :: E
		!complex(8), allocatable :: cbond(:)
		!complex(8), allocatable :: sbond(:)
		real(8), allocatable :: csite(:)
		real(8), allocatable :: ssite(:)
		real(8) :: dsc
		real(8) :: af
		real(8) :: ddw
	end type
	integer :: ne,ne2,vn,n_omp=20
	real(8) :: t(2)=(/1d0,-0.3d0/)
	real(8) :: DJ=0.3d0,V=0.3d0/4d0
	!real(8) :: DJ=0.5d0,V=0.5d0/4d0
	!real(8) :: DJ=1d0/3d0,V=1d0/12d0
	real(8) :: U=10d0
	integer :: opt=1 ! 1 Tr(AdA)
					 ! 2 O=cicj
contains
	subroutine initial()
		integer ::i,l2,l3,l
		real(8), allocatable :: bd0(:),bd1(:)
		real(8) :: q(3)
		allocate(var(-1000:1000))
		call init_random_seed()
		!q=(/1d0/8d0,0d0,0d0/)
		!q=(/1d0/50d0,0d0,0d0/)
		q=0d0
		! lattice 
		latt%a1=(/1d0,1d0,0d0/)
		latt%a2=(/-1d0,1d0,0d0/)
		latt%T1=(/1d0,0d0,0d0/)*16
		latt%T2=(/0d0,1d0,0d0/)*10
		latt%bdc(1)=1d0
		latt%bdc(2)=1d0
		allocate(latt%sub(2,3))
		latt%sub(1,:)=(/0d0,0d0,0d0/)
		latt%sub(2,:)=(/1d0,0d0,0d0/)
		latt%layer=1
		call latt%gen_latt()
		call latt%gen_neb(3)
		call latt%gen_bond(3)
		call check_lattice(101)
		write(*,*)"Total site number is: ",latt%Ns

		allocate(bd0(latt%Ns))
		do i=1,size(bd0)
			!sort_site(k)%val=sin(2d0*pi*sum((latt%bond(0)%bd(k)%r+(/0.5d0,0d0,0d0/))*q))
			bd0(i)=cos(2d0*pi*sum((latt%bond(0)%bd(i)%r)*q))
		enddo

		allocate(bd1(size(latt%bond(1)%bd)))
		do i=1,size(bd1)
			bd1(i)=cos(2d0*pi*sum((latt%bond(1)%bd(i)%r-(/0.5d0,0d0,0d0/))*q))
		enddo

		! cp
		call gen_var(sg=1,nb=0)
		var(iv(0))%val=-1.39752d+00
		var(iv(0))%bd_sg=-1d0

		! ddw
		!call gen_var(sg=3,nb=1,val=bd1)
		call gen_var(sg=3,nb=1)
		do i=1,size(var(iv(0))%bd_sg)
			var(iv(0))%bd_sg(i)=img*ab(latt%bond(var(iv(0))%nb)%bd(i)%i(1))*dwave(i)*&
					!1d0
					!!cos(2d0*pi*sum(q*(latt%bond(var(iv(0))%nb)%bd(i)%r-(/0.5d0,0d0,0d0/))))
		enddo
		var(iv(0))%val=2d-1


		! d-wave sc
		!call gen_var(sg=2,nb=1,val=bd1)
		call gen_var(sg=2,nb=1)
		do i=1,size(var(iv(0))%bd_sg)
			var(iv(0))%bd_sg(i)=dwave(i)
		enddo
		var(iv(0))%val=1d-1

		!call gen_var(sg=2,nb=1)
		!do i=1,size(var(iv(0))%bd_sg)
			!var(iv(0))%bd_sg(i)=dwave(i)*&
					!cos(2d0*pi*sum(q*(latt%bond(var(iv(0))%nb)%bd(i)%r-(/0.5d0,0d0,0d0/))))
		!enddo
		!var(iv(0))%val=1d-5


		!! s-wave sc
		!call gen_var(sg=2,nb=1,val=bd1)
		!var(iv(0))%val=1d-1
		!var(iv(0))%bd_sg=1d0

		! sdw
		!call gen_var(sg=4,nb=0,val=bd0)
		call gen_var(sg=4,nb=0)
		do i=1,size(var(iv(0))%bd_sg)
			var(iv(0))%bd_sg(i)=ab(i)*&
					1d0
					!sin(2d0*pi*sum(q*(latt%bond(var(iv(0))%nb)%bd(i)%r+(/0.5d0,0d0,0d0/))))
					!sin(2d0*pi*sum(q*(latt%bond(var(iv(0))%nb)%bd(i)%r)))
		enddo

		!! d-wave cdw
		!!call gen_var(sg=3,nb=1,bd1)
		!call gen_var(sg=3,nb=1)
		!var(iv(0))%val=1d-1
		!do i=1,size(var(iv(0))%bd_sg)
			!var(iv(0))%bd_sg(i)=dwave(i)
		!enddo

		!! on site cdw
		!!call gen_var(sg=3,nb=0,bd0)
		!call gen_var(sg=3,nb=0)
		!var(iv(0))%val=1.32678d-02
		!do i=1,size(var(iv(0))%bd_sg)
				!var(iv(0))%bd_sg(k)=&
					!!cos(2d0*pi*sum(2d0*q*(latt%bond(var(iv(0))%nb)%bd(i)%r+(/0.5d0,0d0,0d0/))))
					!cos(2d0*pi*sum(2d0*q*(latt%bond(var(iv(0))%nb)%bd(i)%r)))
			!enddo
		!enddo

		! bond order
		call gen_var(sg=3,nb=1,bd1(2:))
		var(iv(0))%val=-2.00021d-01
		var(iv(0))%bd_sg=-1d0

		call gen_var(sg=3,nb=2)
		var(iv(0))%val=-2.00021d-01
		var(iv(0))%bd_sg=-1d0

		! hp
		do l=1,size(t)
			call gen_var(sg=-3,nb=l)
			var(iv(0))%bd_sg=-1d0
			var(iv(0))%val=t(l)
		enddo

		vn=iv(1)
		!! n-n jast
		!do i=0,7
			!call gen_var(sg=11,nb=i)
			!var(iv(0))%val=0d0
		!enddo

		call var_shrink()
	end subroutine
	subroutine export_data(ut)
		integer :: ut,l1,i
		do l1=2,size(var(1:))
			do i=1,size(var(l1)%bd_sg)
				write(ut,"(es17.9$)")latt%bond(var(l1)%nb)%bd(i)%r,&
					latt%bond(var(l1)%nb)%bd(i)%dr,&
					var(l1)%val(var(l1)%i2v(i)),&
					var(l1)%bd_sg(i)
				write(ut,"(x)")
			enddo
			write(ut,"(x/)")
		enddo
	end subroutine
end module

module M_wf
	use global
	use M_matrix
	implicit none
contains
	subroutine iniwf(wf,dwf,cp)
		complex(8) :: H(latt%Ns*2,latt%Ns*2),cH(size(H,1),size(H,2)),D(size(H,1),size(H,2)),wf(:,:),tmp(latt%Ns),Q(latt%Ns,latt%Ns)
		complex(8), optional :: dwf(:,:,:)
		real(8) :: E(size(H,1)),r,cp1
		type(t_var), optional :: cp
		complex(8) :: bd
		integer :: i,j,k,l,n,info,sg
		if(present(cp).and.cp%sg/=0) then
			return
		endif

		call var%Hamilton(H)

		associate(Ns=>latt%Ns)
			if(present(cp)) then
				write(*,"(A$)")"*"
				call heevd(H(:Ns,:Ns),E(:Ns),"N","U",info)
				cp%val=cp%val+(E(ne2)+E(ne2+1))/2d0
				return
			endif
			if(var(1)%sg/=1.and.(.not.is_sc)) then
				call heevd(H(:Ns,:Ns),E(:Ns),"V","U",info)
				call heevd(H(Ns+1:,Ns+1:),E(Ns+1:),"V","U",info)
				E(Ns+1:)=E(Ns+1:)+(E(ne2)+E(ne2+1))/2d0
				E(:Ns)=E(:Ns)-(E(ne2)+E(ne2+1))/2d0
				do i=ne2+1,Ns
					call swap(H(:,i),H(:,i+Ns-ne2))
					call swap(E(i),E(i+Ns-ne2))
				enddo
			else
				call heevd(H,E,"V","U",info)
			endif
			wf=H(:,1:Ns)

			if(present(dwf)) then
				cH=transpose(conjg(H))
				dwf=0d0

				do l=1,vn
					call dHamilton(var(l),H,cH,D)
					select case(opt)
					case(1)
						!method 1
						!$OMP PARALLEL DO
						do i=1,Ns*2
							do j=1,Ns*2
								if(abs(E(i)-E(j))<1d-8) then
									cycle
								endif
								dwf(:,i,l)=dwf(:,i,l)+D(j,i)*H(:,j)/(E(i)-E(j))
							enddo
						enddo
						!$OMP END PARALLEL DO
					case(2)
						!method 2
						!$OMP PARALLEL DO COLLAPSE(2)
						do i=1,Ns
							do j=1,Ns
								if(abs(E(j)-E(i+Ns))<1d-8) then
									Q(i,j)=0d0
									cycle
								endif
								Q(i,j)=D(i+Ns,j)/(E(j)-E(i+Ns))
							enddo
						enddo
						!$OMP END PARALLEL DO
						!$OMP PARALLEL DO PRIVATE(tmp)
						do i=1,Ns*2
							tmp=matmul(H(i,Ns+1:),Q)
							do j=1,Ns*2
								dwf(i,j,l)=sum(tmp*cH(:Ns,j))
							enddo
						enddo
						!$OMP END PARALLEL DO
					end select
				enddo
			endif
		end associate
		!if(any(isnan(real(wf)))) then
			!write(*,*)"NAN in wf"
			!stop
		!endif
		write(*,"(A$)")"*"
	end subroutine
end module

module M_mc_matrix
	use M_matrix
	use global
	implicit none
contains
	subroutine getindex(i,j,cfg,sg,k,m)
		integer :: i,j,sg,cfg(:),k(:),m(:)
		select case(sg) 
		case(1:4) 
			k(1)=i
			m(1)=cfg(j)
		case(5,7) 
			k(1)=j
			m(1)=cfg(i)+latt%Ns
		case(6,8) 
			k(1)=j+ne2
			m(1)=cfg(i)+latt%Ns
		case(0)
			k=(/i,i+ne2/)
			m=(/cfg(j),cfg(j)+latt%Ns/)
		end select
	end subroutine
	subroutine getpb(k,m,sg,D,pb)
		complex(8) :: D(:,:),pb
		integer :: k(:),m(:),sg
		select case(sg) 
		case(0)
			pb=(D(m(1),k(1))*D(m(2),k(2))-D(m(1),k(2))*D(m(2),k(1)))
		case(1:8)
			pb=D(m(1),k(1))
		end select
	end subroutine
	subroutine update(D,k,m,sg,iA)
		complex(8) :: D(:,:),tmp(size(D,1),size(D,2)),ipb(2,2),pb,tmp2(2),tmpA(latt%Ns,latt%Ns)
		complex(8), optional :: iA(:,:)
		integer :: i,j,sg,m(:),k(:)
		tmp=D
		if(present(iA)) then
			tmpA=iA
		endif
		select case(sg)
		case(0)
			pb=1d0/(D(m(1),k(1))*D(m(2),k(2))-D(m(1),k(2))*D(m(2),k(1)))
			ipb(1,1)=pb*D(m(2),k(2))
			ipb(1,2)=-pb*D(m(1),k(2))
			ipb(2,1)=-pb*D(m(2),k(1))
			ipb(2,2)=pb*D(m(1),k(1))
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
					if(present(iA).and.i<=latt%Ns) then
						iA(i,j)=tmpA(i,j)-sum((/tmpA(i,k(1)),tmpA(i,k(2))/)*tmp2)
					endif
				enddo
			enddo
		case(1:8)
			pb=1d0/D(m(1),k(1))
			do j=1,size(D,2)
				if(j==k(1)) then
					tmp2(1)=pb*(tmp(m(1),j)-1d0)
				else
					tmp2(1)=pb*tmp(m(1),j)
				endif
				do i=1,size(D,1)
					D(i,j)=tmp(i,j)-tmp(i,k(1))*tmp2(1)
					if(present(iA).and.i<=latt%Ns) then
						iA(i,j)=tmpA(i,j)-tmpA(i,k(1))*tmp2(1)
					endif
				enddo
			enddo
		end select
	end subroutine
	subroutine update_swap(i,j,nd,cfg,icfg,D,sg,iA)
		complex(8) :: D(:,:)
		complex(8), optional :: iA(:,:)
		integer :: i,j,cfg(:),icfg(:),sg,nd
		select case(sg)
		case(0,2,5)
			call swap(cfg(i),cfg(j))
			call swap(icfg(cfg(i)),icfg(cfg(j)))
		case(1)
			call swap(cfg(i),cfg(j))
			call swap(icfg(cfg(i)),icfg(cfg(j)))
			call swap(D(:,i+ne2),D(:,j))
			if(present(iA)) then
				call swap(iA(:,i+ne2),iA(:,j))
			endif
		case(3)
			call swap(cfg(i),cfg(ne2-nd+1))
			call swap(icfg(cfg(i)),icfg(cfg(ne2-nd+1)))
			call swap(cfg(j),cfg(ne-nd+1))
			call swap(icfg(cfg(j)),icfg(cfg(ne-nd+1)))
			call swap(cfg(ne2-nd+1),cfg(ne-nd+1))
			call swap(icfg(cfg(ne2-nd+1)),icfg(cfg(ne-nd+1)))
			call swap(D(:,i),D(:,ne2-nd+1))
			call swap(D(:,j),D(:,ne-nd+1))
			if(present(iA)) then
				call swap(iA(:,i),iA(:,ne2-nd+1))
				call swap(iA(:,j),iA(:,ne-nd+1))
			endif
			nd=nd-1
		case(4)
			call swap(cfg(i),cfg(ne2-nd))
			call swap(icfg(cfg(i)),icfg(cfg(ne2-nd)))
			call swap(cfg(j),cfg(ne-nd))
			call swap(icfg(cfg(j)),icfg(cfg(ne-nd)))
			call swap(cfg(ne2-nd),cfg(ne-nd))
			call swap(icfg(cfg(ne2-nd)),icfg(cfg(ne-nd)))
			call swap(D(:,i),D(:,ne2-nd))
			call swap(D(:,i+ne2),D(:,ne-nd))
			if(present(iA)) then
				call swap(iA(:,i),iA(:,ne2-nd))
				call swap(iA(:,i+ne2),iA(:,ne-nd))
			endif
			nd=nd+1
		case(6)
			call swap(cfg(i),cfg(j))
			call swap(icfg(cfg(i)),icfg(cfg(j)))
			call swap(D(:,i),D(:,j))
			if(present(iA)) then
				call swap(iA(:,i),iA(:,j))
			endif
		case(7)
			call swap(cfg(i),cfg(ne2-nd+1))
			call swap(icfg(cfg(i)),icfg(cfg(ne2-nd+1)))
			call swap(cfg(j),cfg(ne-nd+1))
			call swap(icfg(cfg(j)),icfg(cfg(ne-nd+1)))
			call swap(D(:,i),D(:,ne2-nd+1))
			call swap(D(:,j),D(:,ne-nd+1))
			if(present(iA)) then
				call swap(iA(:,i),iA(:,ne2-nd+1))
				call swap(iA(:,j),iA(:,ne-nd+1))
			endif
			nd=nd-1
		case(8)
			call swap(cfg(j),cfg(ne2-nd))
			call swap(icfg(cfg(j)),icfg(cfg(ne2-nd)))
			call swap(cfg(i),cfg(ne-nd))
			call swap(icfg(cfg(i)),icfg(cfg(ne-nd)))
			call swap(D(:,j),D(:,ne2-nd))
			call swap(D(:,j+ne2),D(:,ne-nd))
			if(present(iA)) then
				call swap(iA(:,j),iA(:,ne2-nd))
				call swap(iA(:,j+ne2),iA(:,ne-nd))
			endif
			nd=nd+1
		end select
	end subroutine
	function detO(cfg,A,D)
		complex(8) :: A(:,:),D(:,:),detO,pb
		integer :: i,j,k(2),m(2),sg,cfg(:)
		detO=0d0
		do i=1,latt%Ns
			do j=1,latt%Ns*2
				detO=detO+A(cfg(i),j)*D(j,i)
			enddo
		enddo
	end function
	function cs(i,nd)
		integer :: i,cs,nd
		if(i<=(ne2-nd)) then
			cs=1
		elseif(i<=ne2) then
			cs=2
		elseif(i<=ne-nd) then
			cs=-1
		else
			cs=0
		endif
	end function
	subroutine get_case(i,j,ud,nd,sg)
		integer :: i,j,ud,sg,ci,cj,nd
		ci=0
		cj=0
		sg=-1
		ci=cs(i,nd)
		cj=cs(j,nd)
		if(ci>cj) then
			call swap(ci,cj)
		elseif(ci==cj) then
			return
		endif
		if(ci==0.and.cj==1) then
			sg=1
		elseif(ci==-1.and.cj==2.and.ud==1) then
			sg=2
		elseif(ci==0.and.cj==2) then
			select case(ud)
			case(1)
				sg=3
			case(-1)
				sg=7
			end select
		elseif(ci==-1.and.cj==1) then
			select case(ud)
			case(1)
				sg=4
			case(-1)
				sg=8
			case(0)
				sg=0
			end select
		elseif(ci==-1.and.cj==0) then
			sg=5
		elseif(ci==1.and.cj==2.and.ud==-1) then
			sg=6
		else
			sg=-1
		endif
	end subroutine
	subroutine update_jast(i,j,sg,dja,djan,ja,jan,lO)
		integer :: djan(:,:),dja(:,:,:),i,j,sg,l
		real(8) :: ja(:,:),jan(:)
		complex(8) :: lO(:)
		if(sg==0) then
			return
		endif
		do l=1,size(djan,2)
			lO(vn+l)=lO(vn+l)+(-djan(i,l)+djan(j,l)-dja(i,j,l))*2d0+dja(i,i,l)+dja(j,j,l)
			djan(:,l)=djan(:,l)-dja(:,i,l)+dja(:,j,l)
		enddo
		jan=jan-ja(:,i)+ja(:,j)
	end subroutine
	function jast(i,j,sg,ja,jan)
		real(8) :: jast,ja(:,:),jan(:),df
		integer :: i,j,sg,k,l
		jast=0d0
		if(sg==0) then
			return
		endif
		jast=(-jan(i)+jan(j)-ja(i,j))*2d0+ja(i,i)+ja(j,j)
	end function
end module

module M_tJ
	use M_omp_rand
	use M_mc_matrix
	implicit none
contains
	subroutine two(rnd,i,j,nd,sg)
		integer :: i,j,sg,n,n1,n0,nd
		type(randomNumberSequence) :: rnd
		n1=ne2-nd
		n0=latt%Ns-ne+nd
		n=mt_random_number(rnd,n1*n0*2+n1*n1)
		if(n<=n1*n1) then
			sg=0
			i=(n-1)/n1+1
			j=ne2+mod(n-1,n1)+1
		elseif(n<=n1*n1+n1*n0) then
			n=n-n1*n1
			sg=1
			i=(n-1)/n0+1
			j=ne-nd+mod(n-1,n0)+1
		else
			n=n-n1*n1-n1*n0
			sg=5
			i=ne2+(n-1)/n0+1
			j=ne-nd+mod(n-1,n0)+1
		endif
	end subroutine
	subroutine get_energy(cfg,icfg,nd,D,ja,jan,El)
		complex(8) :: pb,D(:,:)
		complex(8) :: El
		real(8) :: ja(:,:),jan(:)
		integer :: cfg(:),icfg(:),i,j,n,l,k(2),m(2),sg,nd
		El=0d0
		do i=1,ne
			do l=1,size(t)
				do n=1,size(latt%neb(cfg(i))%nb(l)%neb)
					j=icfg(latt%neb(cfg(i))%nb(l)%neb(n))
					call get_case(i,j,0,nd,sg)
					select case(sg) 
					case(1) 
						call getindex(i,j,cfg,sg,k,m)
						call getpb(k,m,sg,D,pb)
						El=El-(t(l)*conjg(pb))&
							*exp(jast(cfg(i),cfg(j),sg,ja,jan))&
							*latt%neb(cfg(i))%nb(l)%bdc(n)
					case(5) 
						call getindex(i,j,cfg,sg,k,m)
						call getpb(k,m,sg,D,pb)
						El=El+(t(l)*conjg(pb))&
							*exp(jast(cfg(i),cfg(j),sg,ja,jan))&
							*latt%neb(cfg(i))%nb(l)%bdc(n)
					case(0)
						if(l==1.and.i<=ne2) then
							call getindex(i,j,cfg,sg,k,m)
							call getpb(k,m,sg,D,pb)
							El=El+0.5d0*DJ*conjg(pb)-0.25d0*DJ&
								-0.25d0*DJ&
								+V
						endif
					case(-1)
						if(l==1.and.(j<i)) then
							El=El+0.25d0*DJ&
								-0.25d0*DJ&
								+V
						endif
					end select
				enddo
			enddo
		enddo
		!write(*,*)El
		!El=(El-(4d0*V-DJ)*(ne-latt%Ns/2))/latt%Ns
		!El=(El-(4d0*V)*(ne-latt%Ns/2))/latt%Ns
		El=El/latt%Ns
	end subroutine
	subroutine spin_order(icfg,nd,lphy,af)
		real(8) :: lphy(:),af
		integer :: icfg(:),i,nd
		af=0d0
		do i=1,latt%Ns
			lphy(i)=0.5d0*mod(cs(icfg(i),nd),2)
			af=af+0.5*ab(i)*mod(cs(icfg(i),nd),2)
		enddo
		af=af/latt%Ns
	end subroutine
	subroutine charge_order(icfg,nd,lphy)
		real(8) :: lphy(:)
		integer :: icfg(:),i,nd
		do i=1,latt%Ns
			lphy(i)=abs(cs(icfg(i),nd))
		enddo
	end subroutine
	subroutine ddw_order(cfg,icfg,nd,D,ja,jan,lphy)
		complex(8) :: D(:,:),pb
		real(8) :: lphy,ja(:,:),jan(:)
		integer :: cfg(:),icfg(:),b1,j1,i1,sg,nd,k(2),m(2)
		lphy=0d0
		do b1=1,size(latt%bond(1)%bd)
			i1=latt%bond(1)%bd(b1)%i(1)
			j1=latt%bond(1)%bd(b1)%i(2)
			if(icfg(i1)<=ne2.and.icfg(j1)>ne) then
				k(1)=icfg(i1)
				m(1)=j1
				sg=1
			elseif(icfg(j1)<=ne2.and.icfg(i1)>ne) then
				k(1)=icfg(j1)
				m(1)=i1
				sg=-1
				call swap(i1,j1)
			elseif(icfg(i1)>ne2.and.icfg(i1)<=ne.and.icfg(j1)>ne) then
				k(1)=icfg(j1)
				m(1)=i1+latt%Ns
				sg=-1
			elseif(icfg(j1)>ne2.and.icfg(j1)<=ne.and.icfg(i1)>ne) then
				k(1)=icfg(i1)
				m(1)=j1+latt%Ns
				sg=1
				call swap(i1,j1)
			else
				cycle
			endif
			call getpb(k,m,1,D,pb)
			lphy=lphy+imag(pb)*sg*dwave(b1)*exp(jast(i1,j1,1,ja,jan))*ab(i1)
		enddo
		lphy=lphy/latt%Ns
	end subroutine
	subroutine dsc_corelation(cfg,icfg,nd,D,ja,jan,lphy)
		complex(8) :: D(:,:),pb
		real(8) :: lphy,ja(:,:),jan(:)
		integer :: cfg(:),icfg(:),i1,j1,i2,j2,nd,k(2),m(2),b1,b2
		lphy=0d0
		do b1=1,size(latt%bond(1)%bd)
			i1=latt%bond(1)%bd(b1)%i(1)
			j1=latt%bond(1)%bd(b1)%i(2)
			if(icfg(i1)<=ne2.and.icfg(j1)>ne2.and.icfg(j1)<=ne) then
				k(1)=icfg(i1)
				m(2)=j1+latt%Ns
			elseif(icfg(j1)<=ne2.and.icfg(i1)>ne2.and.icfg(i1)<=ne) then
				k(1)=icfg(j1)
				m(2)=i1+latt%Ns
			else
				cycle
			endif
			do b2=1,size(latt%bond(1)%bd)
				i2=latt%bond(1)%bd(b2)%i(1)
				j2=latt%bond(1)%bd(b2)%i(2)
				if(icfg(i2)>ne.and.icfg(j2)>ne) then
					m(1)=i2
					k(2)=icfg(j2)
				else
					cycle
				endif
				call getpb(k,m,0,D,pb)
				lphy=lphy+real(pb)*dwave(b1)*dwave(b2)&
					*exp(jast(i1,i2,1,ja,jan))&
					*exp(jast(j1,j2,1,ja,jan-ja(:,i1)+ja(:,i2)))
				m(1)=j2
				k(2)=icfg(i2)
				call getpb(k,m,0,D,pb)
				lphy=lphy+real(pb)*dwave(b1)*dwave(b2)&
					*exp(jast(i1,i2,1,ja,jan))&
					*exp(jast(j1,j2,1,ja,jan-ja(:,i1)+ja(:,i2)))
			enddo
		enddo
		lphy=lphy/(latt%Ns**2)
	end subroutine
end module

module M_hubbard
	use M_omp_rand
	use M_mc_matrix
	implicit none
contains
	subroutine two(rnd,i,j,nd,sg)
		integer :: i,j,sg,n,n1,n0,nd
		type(randomNumberSequence) :: rnd
		n1=ne2-nd
		n0=latt%Ns-ne+nd
		n=mt_random_number(rnd,ne2*(latt%Ns-ne2)*2)
		if(n<=n1*n0) then
			sg=1
			i=(n-1)/n0+1
			j=ne-nd+mod(n-1,n0)+1
		elseif(n<=n1*n0+nd*n1) then
			n=n-n1*n0
			sg=2
			i=n1+(n-1)/n1+1
			j=ne2+mod(n-1,n1)+1
		elseif(n<=n1*n0+nd*n1+nd*n0) then
			n=n-n1*n0-nd*n1
			sg=3
			i=n1+(n-1)/n0+1
			j=ne-nd+mod(n-1,n0)+1
		elseif(n<=ne2*(latt%Ns-ne2)) then
			n=n-n1*n0-nd*n1-nd*n0
			sg=4
			i=(n-1)/n1+1
			j=ne2+mod(n-1,n1)+1
		elseif(n<=n1*n0+ne2*(latt%Ns-ne2)) then
			n=n-ne2*(latt%Ns-ne2)
			sg=5
			i=ne2+(n-1)/n0+1
			j=ne-nd+mod(n-1,n0)+1
		elseif(n<=n1*n0+nd*n1+ne2*(latt%Ns-ne2)) then
			n=n-ne2*(latt%Ns-ne2)-n1*n0
			sg=6
			j=(n-1)/nd+1
			i=n1+mod(n-1,nd)+1
		elseif(n<=n1*n0+nd*n1+nd*n0+ne2*(latt%Ns-ne2)) then
			n=n-ne2*(latt%Ns-ne2)-n1*n0-nd*n1
			sg=7
			i=n1+(n-1)/n0+1
			j=ne-nd+mod(n-1,n0)+1
		elseif(n<=ne2*(latt%Ns-ne2)*2) then
			n=n-ne2*(latt%Ns-ne2)-n1*n0-nd*n1-nd*n0
			sg=8
			j=mod(n-1,n1)+1
			i=ne2+(n-1)/n1+1
		endif
	end subroutine
	!subroutine get_energy(cfg,icfg,nd,nn,wf,ja,A,iA,El)
		!complex(8) :: A(:,:),iA(:,:),vu(size(A,1),2),Y(2,2),pb,wf(:,:)
		!complex(8) :: El
		!real(8) :: ja(:)
		!integer :: cfg(:),icfg(:),i,j,k,l,sg,cr(2),nd,n,ud,nn(:),nndif(size(nn))
		!El=0d0
		!do n=1,ne
			!if(n>ne2) then
				!ud=-1
				!i=n-nd
			!else
				!ud=1
				!i=n
			!endif
			!do l=1,size(t)
				!do k=1,size(latt%neb((cfg(i))%nb(l+1)%neb)
					!j=icfg(latt%neb((cfg(i))%nb(l+1)%neb(k))
					!call get_case(i,j,ud,nd,sg)
					!select case(sg)
					!case(1:4)
						!call nn_dif(cfg(i),cfg(j),cfg,icfg,sg,nd,nndif)
						!call getdiff(wf,i,j,cfg(i),cfg(j),sg,vu,cr)
						!call det(vu,cr,A,iA,Y,pb)
						!El=El-(t(l)*conjg(pb)*jast(nndif,ja))
					!case(5:8)
						!call nn_dif(cfg(i),cfg(j),cfg,icfg,sg,nd,nndif)
						!call getdiff(wf,i,j,cfg(i),cfg(j),sg,vu,cr)
						!call det(vu,cr,A,iA,Y,pb)
						!El=El+(t(l)*conjg(pb)*jast(nndif,ja))
					!end select
				!enddo
			!enddo
		!enddo
		!El=El+U*nd+V*nn(2)
		!El=El/latt%Ns
		!!stop
	!end subroutine
	!subroutine spin_order(cfg,nd,lphy)
		!real(8) :: lphy
		!integer :: cfg(:),i,x1(2),nd
		!lphy=0d0
		!do i=1,ne-nd*2
			!call latt_one2two(cfg(i),latt%Ns,x1)
			!!lphy=lphy+sign(-1d0,i-ne2-0.1d0)*(0.5d0-mod(x1(1)+x1(2),2))
			!lphy=lphy+(1-(i-1)/(ne2-nd)*2)*(0.5d0-mod(sum(x1),2))
		!enddo
		!lphy=lphy/latt%Ns
	!end subroutine
	!!subroutine charge_order(cfg,icfg,nd,wf,A,iA,lphy)
		!!complex(8) :: A(:,:),iA(:,:),vu(latt%Ns,2),Y(2,2),pb,wf(:,:),a11
		!!complex(8) :: lphy
		!!integer :: cfg(:),icfg(:),i,n,m,j,k,l,inb,sg,cr(2),x1(2),x2(2),nr(2),dr(2),r(2),nd
		!!lphy=0d0
		!!lphy=lphy/latt%Ns
	!!end subroutine
	!subroutine ddw_order(cfg,icfg,nd,wf,A,iA,ja,lphy)
		!complex(8) :: A(:,:),iA(:,:),vu(latt%Ns,2),Y(2,2),pb,wf(:,:),a11
		!real(8) :: lphy,ja(:)
		!integer :: cfg(:),icfg(:),i,n,m,j,k,l,inb,sg,cr(2),x1(2),nd,ud
		!lphy=0d0
		!do n=1,ne
			!ud=1-(n-1)/ne2*2
			!i=n-nd*(1-ud)/2
			!call latt_one2two(cfg(i),latt%Ns,x1)
			!do inb=1,4
				!j=icfg(latt%neb((cfg(i),inb,1))
				!call get_case(i,j,ud,nd,sg)
				!select case(sg) 
				!case(1:4)
					!call getdiff(wf,i,j,cfg(i),cfg(j),sg,vu,cr)
					!call det(vu,cr,A,iA,Y,pb)
					!lphy=lphy+imag(conjg(pb)*jast(cfg(i),cfg(j),cfg,icfg,sg,nd,ja))*&
						!(0.5d0-mod(inb,2))*(0.5d0-mod(sum(x1),2))
				!case(5:8)
					!call getdiff(wf,i,j,cfg(i),cfg(j),sg,vu,cr)
					!call det(vu,cr,A,iA,Y,pb)
					!lphy=lphy-imag(conjg(pb)*jast(cfg(i),cfg(j),cfg,icfg,sg,nd,ja))*&
						!(0.5d0-mod(inb,2))*(0.5d0-mod(sum(x1),2))
				!end select
			!enddo
		!enddo
		!lphy=lphy/latt%Ns
	!end subroutine
	!!subroutine dsc_order(cfg,icfg,nd,wf,A,iA,ja,lphy)
		!!complex(8) :: A(:,:),iA(:,:),vu(latt%Ns,2),Y(2,2),pb,wf(:,:),a11
		!!real(8) :: lphy,ja(:)
		!!integer :: cfg(:),icfg(:),i,n,m,j,k,l,inb,sg,cr(2),x1(2),x2(2),nr(2),dr(2),r(2),nd
		!!do i=ne2-nd+1,ne-nd
			!!do inb=1,4,1
				!!j=icfg(latt%neb((cfg(i),inb,1))
				!!if(j>ne2) then
					!!vu(:,1)=wf(cfg(j),:)
					!!cr=(/i,-1/)
					!!call det(vu,cr,A,iA,Y,pb)
					!!lphy=lphy+real(conjg(pb)*jast(cfg(i),cfg(j),cfg,icfg,nd,ja))*&
						!!(0.5d0-mod(inb,2))
				!!endif
			!!enddo
		!!enddo
	!!end subroutine
	!!subroutine spin_corelation(cfg,icfg,nd,wf,A,iA,lphy)
		!!complex(8) :: A(:,:),iA(:,:),vu(latt%Ns,2),Y(2,2),pb,wf(:,:),a11
		!!real(8) :: lphy
		!!integer :: cfg(:),icfg(:),i,n,m,j,k,l,inb,sg,cr(2),x1(2),x2(2),nr(2),dr(2),r(2),nd
	!!end subroutine
	!subroutine dsc_corelation(cfg,icfg,nd,wf,A,iA,ja,lphy)
		!complex(8) :: A(:,:),iA(:,:),vu(latt%Ns,2),Y(2,2),pb,wf(:,:),a11
		!real(8) :: lphy,ja(:)
		!integer :: cfg(:),icfg(:),i,j,ii,jj,inb1,inb2,cr(2),x1(2),x2(2),nr(2),dr(2),r(2),nd
		!lphy=0d0
		!do i=ne2-nd+1,ne-nd
			!do inb1=1,4,1
				!j=icfg(latt%neb((cfg(i),inb1,1))
				!if(j<=ne2) then
					!do jj=ne2-nd+1,latt%Ns
						!do inb2=1,4,1
							!ii=icfg(latt%neb((cfg(jj),inb2,1))
							!if(ii<=(ne2-nd)) then
								!cr=(/ii+ne2,j/)
							!elseif(ii>(ne-nd)) then
								!cr=(/ii,j/)
							!else
								!cycle
							!endif
							!vu(:,1)=wf(cfg(i)+latt%Ns,:)
							!vu(:,2)=wf(cfg(jj),:)
							!call det(vu,cr,A,iA,Y,pb)
							!lphy=lphy+real(conjg(pb)*jast(cfg(i),cfg(ii),cfg,icfg,1,nd,ja)*jast(cfg(j),cfg(jj),cfg,icfg,1,nd,ja))*&
								!(0.5d0-mod(inb1+inb2,2))
						!enddo
					!enddo
				!endif
			!enddo
		!enddo
		!lphy=lphy/(latt%Ns*latt%Ns)
	!end subroutine
end module

module M_vmc
	use M_tJ
	!use M_hubbard
	use M_wf
	implicit none
	type, extends(t_sort), private :: t_mysort
		real(8), allocatable :: var(:)
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
			allocate(tmp%var(size(self%var)))
			tmp%var=a%var
			a%var=self%var
			self%var=tmp%var
			deallocate(tmp%var)
			deallocate(tmp)
		end select
	end subroutine
	subroutine mc(rnd,wf,Nmc,cfg,nd,phy,dwf,S,g)
		complex(8) :: pb,A(latt%Ns,latt%Ns),iA(latt%Ns,latt%Ns),dA(latt%Ns,latt%Ns),wf(:,:),D(size(wf,1),size(wf,2))
		complex(8), optional :: S(:,:),dwf(:,:,:)
		real(8), optional :: g(:)
		real(8) :: rpb
		real(8) :: ja(latt%Ns,latt%Ns),jan(latt%Ns)
		integer :: dja(latt%Ns,latt%Ns,size(var(1:))-vn),djan(latt%Ns,size(var(1:))-vn)
		complex(8) :: lO(size(var(1:))),lS(size(lO),size(lO)),lg(size(lO)),O(size(lO))
		complex(8) :: El,E
		type(t_phy) :: phy,lphy
		integer :: icfg(latt%Ns),i,j,ti,tj,l,p,sg,n,Nmc(:),k(2),m(2),acp,info
		integer :: cfg(latt%Ns),nd,nn(max(size(var(1:))-vn,2)),nndif(size(nn))
		logical :: flag
		type(randomNumberSequence) :: rnd

		if(mod(latt%Ns,2)/=0) stop "Number of site is not even, exit!!"
		if(ne/=ne2*2) stop "Number of electrons is not even, exit!!"

		flag=present(S)

		!do i=1,latt%Ns
			!cfg(i)=i
		!enddo
		!nd=0
		!call fisher_yates_shuffle(cfg,latt%Ns)

		do i=1,latt%Ns
			icfg(cfg(i))=i
		enddo

		n=1
		!get initial matrix
		do 
			do i=1,latt%Ns
				if(i<=ne2) then
					A(i,:)=wf(cfg(i),:)
				elseif(i>=ne-nd+1)  then
					A(i,:)=wf(cfg(i)+latt%Ns,:)
				else
					A(i,:)=wf(cfg(i-ne2)+latt%Ns,:)
				endif
			enddo
			call mat_inv(A,info)
			if(info==0) then
				exit
			endif
			if(n>500) then
				!OMP CRITICAL 
				write(*,*)"get initial matrix wronge",info
				!OMP END CRITICAL 
				exit
			endif
			call two(rnd,i,j,nd,sg)
			call swap(cfg(i),cfg(j))
			n=n+1
		enddo
		D=matmul(wf,A)

		

		n=0
		acp=0

		phy%E=0d0
		phy%dsc=0d0
		phy%af=0d0
		phy%ddw=0d0
		if(allocated(phy%csite)) then
			allocate(lphy%csite(size(phy%csite)),lphy%ssite(size(phy%ssite)))
			phy%csite=0d0
			phy%ssite=0d0
		endif
		lphy=phy

		if(flag) then
			O=0d0
			S=0d0
			g=0d0
			lS=0d0
			g=0d0
			lO=0d0
		endif

		ja=0d0
		jan=0d0
		dja=0
		djan=0
		do l=1,size(djan,2)
			do p=1,size(var(vn+l)%n)
				i=latt%bond(var(vn+l)%nb)%bd(var(vn+l)%n(p))%i(1)
				j=latt%bond(var(vn+l)%nb)%bd(var(vn+l)%n(p))%i(2)
				ja(i,j)=-var(vn+l)%val
				dja(i,j,l)=-1
				jan(i)=jan(i)+ja(i,j)*abs(cs(icfg(j),nd))
				djan(i,l)=djan(i,l)+dja(j,i,l)*abs(cs(icfg(j),nd))
				if(i/=j) then
					ja(j,i)=-var(vn+l)%val
					dja(j,i,l)=-1
					jan(j)=jan(j)+ja(j,i)*abs(cs(icfg(i),nd))
					djan(j,l)=djan(j,l)+dja(i,j,l)*abs(cs(icfg(i),nd))
				endif
			enddo
			do p=1,latt%Ns
				lO(vn+l)=lO(vn+l)+abs(cs(icfg(p),nd))*djan(p,l)
			enddo
		enddo

		do
			n=n+1
			call two(rnd,i,j,nd,sg)
			!write(*,*)i,j,nd,sg
			call getindex(i,j,cfg,sg,k,m)
			call getpb(k,m,sg,D,pb)
			!call random_number(rpb)
			rpb=getrandomreal(rnd)
			if(rpb<real(pb*conjg(pb))*exp(jast(cfg(i),cfg(j),sg,ja,jan)*2d0)) then
			!if(rpb<1d0/(1d0+1d0/(real(pb*conjg(pb))*exp(jast(cfg(i),cfg(j),sg,ja,jan)*2d0)))) then
				acp=acp+1
				call update_jast(cfg(i),cfg(j),sg,dja,djan,ja,jan,lO)
				select case(opt)
				case(1)
					if(flag) then
						call update(D,k,m,sg,A)
						call update_swap(i,j,nd,cfg,icfg,D,sg,A)
					else
						call update(D,k,m,sg)
						call update_swap(i,j,nd,cfg,icfg,D,sg)
					endif
				case(2)
					call update(D,k,m,sg)
					call update_swap(i,j,nd,cfg,icfg,D,sg)
				end select
				!if(mod(acp,500)) then
					!do i=1,latt%Ns
						!if(i<=ne2) then
							!iA(i,:)=wf(cfg(i),:)
						!elseif(i>=ne-nd+1)  then
							!iA(i,:)=wf(cfg(i)+latt%Ns,:)
						!else
							!iA(i,:)=wf(cfg(i-ne2)+latt%Ns,:)
						!endif
					!enddo
					!if(sum(abs(matmul(A,iA)-diag(1d0,size(A,1))))>1d-9) then
						!!$OMP CRITICAL
						!write(*,*)sum(abs(matmul(A,iA)-diag(1d0,size(A,1)))),n
						!stop
						!!$OMP END CRITICAL
					!endif
				!endif
			endif
			!!$OMP MASTER
			!if(n==Nmc(1)) write(*,"(A$)")"#"
			!!$OMP END MASTER
			if(n>Nmc(1).and.mod(n-Nmc(1),Nmc(2))==0) then
				call get_energy(cfg,icfg,nd,D,ja,jan,lphy%E)
				if(flag) then
					do ti=1,vn
						select case(opt)
						case(1)
							lO(ti)=0d0
							do i=1,latt%Ns
								if(i<=ne2) then
									lO(ti)=lO(ti)+sum(dwf(cfg(i),:latt%Ns,ti)*A(:,i))
								elseif(i>=ne-nd+1) then
									lO(ti)=lO(ti)+sum(dwf(cfg(i)+latt%Ns,:latt%Ns,ti)*A(:,i))
								else
									lO(ti)=lO(ti)+sum(dwf(cfg(i-ne2)+latt%Ns,:latt%Ns,ti)*A(:,i))
								endif
							enddo
						case(2)
							lO(ti)=detO(cfg,dwf(:,:,ti),D)
						end select
					enddo
					do ti=1,size(lO)
						do tj=1,size(lO)
							lS(ti,tj)=conjg(lO(ti))*lO(tj)
						enddo
					enddo
					lg=real(lphy%E*lO)
					S=S+lS
					g=g+lg
					O=O+lO
				else
					if(allocated(phy%csite)) then
						call spin_order(icfg,nd,lphy%ssite,lphy%af)
						call charge_order(icfg,nd,lphy%csite)
						call ddw_order(cfg,icfg,nd,D,ja,jan,lphy%ddw)
						call dsc_corelation(cfg,icfg,nd,D,ja,jan,lphy%dsc)
					endif
				endif
				phy%E=phy%E+lphy%E
				if(allocated(phy%csite)) then
					phy%csite=phy%csite+lphy%csite
					phy%ssite=phy%ssite+lphy%ssite
					phy%dsc=phy%dsc+lphy%dsc
					phy%af=phy%af+lphy%af
					phy%ddw=phy%ddw+lphy%ddw
				endif
			endif
			if(n>=(Nmc(1)+Nmc(2)*Nmc(3))) then
				phy%E=phy%E/Nmc(3)
				phy%dsc=abs(phy%dsc/Nmc(3))
				phy%af=abs(phy%af/Nmc(3))
				phy%ddw=abs(phy%ddw/Nmc(3))
				if(allocated(phy%csite)) then
					phy%csite=phy%csite/Nmc(3)
					phy%ssite=phy%ssite/Nmc(3)
				endif
				if(flag) then
					S=S/Nmc(3)
					g=g/Nmc(3)
					O=O/Nmc(3)
					do ti=1,size(O)
						do tj=1,size(O)
							S(ti,tj)=S(ti,tj)-conjg(O(ti))*O(tj)
						enddo
					enddo
					g=2d0*(g-real(phy%E*O))*latt%Ns
				endif
				if(allocated(lphy%csite)) then
					deallocate(lphy%csite,lphy%ssite)
				endif
				exit
			endif
		enddo
		!!$OMP MASTER
		!write(*,"('accept/total number is ',I6,'/',I6)")acp,n
		!!$OMP END MASTER
	end subroutine
	subroutine variational(Nmc,cfg,nd)
		use lapack95, only: heevd,heevr
		real(8) :: dvar(sum(var(1:)%n)),eg(size(dvar)),dt,scv
		complex(8) :: wf(latt%Ns*spin,latt%Ns),dwf(latt%Ns*spin,latt%Ns*spin,vn)
		complex(8) :: S(size(dvar),size(dvar),n_omp)
		real(8) :: g(size(dvar),n_omp),cg(size(dvar)),h(size(dvar))
		integer :: m,i,j,k,info,Nmc(:),l,ml,ll
		integer :: cfg(:,:),nd(:)
		real(8) :: er
		type(randomNumberSequence) :: rnd(n_omp)
		type(t_phy) :: phy(n_omp)
		type(t_mysort), allocatable :: allE(:)
		dt=0.2d0
		!dt=1d0
		m=50
		l=1
		ll=0
		scv=0d0
		allocate(allE(100))
		cg=1d0
		h=0d0
		ml=1
		do
			write(*,"(I3$)")ne
			write(*,"(es9.2$)")var(1:)%val(1)
			call iniwf(wf,dwf)
			!$OMP PARALLEL DO
			do k=1,n_omp
				call mt_init_random_seed(rnd(k),seed=k*10)
				call mc(rnd(k),wf,Nmc,cfg(:,k),nd(k),phy(k),dwf,S(:,:,k),g(:,k))
				!call finalize_RandomNumberSequence(rnd(k))
			enddo
			!$OMP END PARALLEL DO
			er=real(phy(1)%E)**2
			do k=2,n_omp
				phy(1)%E=phy(1)%E+phy(k)%E
				S(:,:,1)=S(:,:,1)+S(:,:,k)
				g(:,1)=g(:,1)+g(:,k)
				er=er+real(phy(k)%E)**2
			enddo
			phy(1)%E=phy(1)%E/n_omp
			S(:,:,1)=S(:,:,1)/n_omp
			g(:,1)=g(:,1)/n_omp
			er=sqrt(abs(er/n_omp-real(phy(1)%E)**2))
			!S(:,:,1)=S(:,:,1)+diag(1d-3,size(S(:,:,1),1))
			!write(*,*)g(:,1)
			!write(*,*)S(:,:,1)
			call heevd(S(:,:,1),eg,"V")
			!write(*,*)eg
			eg=eg+abs(min(eg(1),0d0))+1d-2
			!write(*,*)eg
			do i=1,size(g(:,1))
				dvar(i)=0d0
				do j=1,size(g(:,1))
					do k=1,size(g(:,1))
						!if(abs(eg(k)/eg(size(g)))<1d-3) then
							!cycle
						!endif
						dvar(i)=dvar(i)+real(S(i,k,1)*conjg(S(j,k,1))*g(j,1)/eg(k))
					enddo
				enddo
			enddo
			!dvar=g(:,1)*dt
			write(*,"(es11.4$)")real(phy(1)%E)
			write(*,"(es9.2$)")er
			!write(30,"(es13.5$)")real(latt%Ns-ne)/latt%Ns,var(1:)%val,real(phy(1)%E),er,g(:,1),dvar
			!write(30,"(x)")
			!return
			if(.not.allocated(allE(l)%var)) then
				allocate(allE(l)%var(size(var(1:))))
			endif
			allE(l)%val=real(phy(1)%E)
			allE(l)%var=var(1:)%val
			do i=max(l-m+1,1),l-1
				scv=scv+sign(1d0,real(phy(1)%E)-allE(i)%val)
				if(max(l-m+1,1)>1) then
					scv=scv+sign(1d0,allE(max(l-m+1,1)-1)%val-allE(i)%val)
				endif
			enddo
			write(*,"(es9.2$)")scv/(min(l*max(l-1,1),m*(m-1))),dt
			if(allE(ml)%val>allE(l)%val) then
				ml=l
			endif
			!if(l>(size(var(1:))+5)*size(var(1:))*2.or.l==size(allE)) then
			!if(l>size(var(1:))*2.or.l==size(allE)) then
			if(allE(max(l-1,1))%val<(real(phy(1)%E)-er*0.8d0)) then
				write(*,"(A)")"var err"
				if(dt<7d-3) then
					var(1:)%val=allE(l)%var*0.4d0+allE(max(l-1,1))%var*0.6d0
					dt=dt*0.4d0
				elseif(dt<1d-1) then
					!var(1:)%val=allE(max(l-1,1))%var+dvar*dt
					var(1:)%val=allE(l)%var*0.8d0+allE(max(l-1,1))%var*0.2d0
					dt=dt*0.8d0
				else
					var(1:)%val=allE(l)%var*0.6d0+allE(max(l-1,1))%var*0.4d0
					dt=dt*0.6d0
				endif
				cycle
			endif
			!if(abs(var(2)%val)<5d-3) then
				!exit
			!endif
			if(all(abs((allE(max(l-1,1))%var-allE(l)%var)/allE(l)%var)<0.01d0).or.(dt<1d-3).or.(l==size(allE))) then
				ll=ll+1
				if((ll>10).or.(l==size(allE)).or.(dt<1d-3)) then
					write(*,*)"variational finish",l
					call allE(:l)%qsort()
					var(1:)%val=0d0
					i=1
					do
						if(abs(allE(i)%val-allE(1)%val)<er/3d0) then
							var(1:)%val=var(1:)%val+allE(i)%var
							i=i+1
						else
							exit
						endif
					enddo
					var(1:)%val=var(1:)%val/real(i-1)
					write(20,"(x)")
					do l=1,size(allE)
						if(allocated(allE(l)%var)) then
							deallocate(allE(l)%var)
						endif
					enddo
					deallocate(allE)
					exit
				endif
			endif
			!if(l>50.and.abs(scv)/(min(l*max(l-1,1),m*(m-1)))<9d-3.and.scv<0.or.l==size(allE)) then
			!if(l>(ml+10)) then
				!write(*,*)"variational finish",l
				!call qsort(allE(:l))
				!var(1:)%val=0d0
				!i=1
				!do
					!if(abs(allE(i)%val-allE(1)%val)<er/3d0) then
						!var(1:)%val=var(1:)%val+allE(i)%var
						!i=i+1
					!else
						!exit
					!endif
				!enddo
				!var(1:)%val=var(1:)%val/real(i-1)
				!write(20,"(x)")
				!exit
			!endif
			!h=-g(:,1)+sum(-(-g(:,1)-cg)*g(:,1))/sum(cg*cg)*h
			!call oned_search(allE(l)%val,er,h,dt,Nmc,cfg,nd)
			!cg=-g(:,1)
			!if(mod(l,size(var(1:))+5)==0) then
				!ml=l+1
				!cg=1d0
				!h=0d0
				!call random_number(var(1:vn)%val)
				!var(1:vn)%val=sign(1d0,allE(l)%var)*var(1:vn)%val
			!endif
			var(1:)%val=var(1:)%val-dvar*dt
			write(*,"(i4)")l
			write(20,"(I4$)")ne
			write(20,"(es13.5$)")var(1:)%val,real(g(:,1)),real(phy(1)%E),er
			write(20,"(x)")
			!Nmc(1)=500
			!Nmc(1)=1000
			l=l+1
		enddo
	end subroutine
	subroutine vmc_main()
		integer :: n,i,j,l,k,Nmc(3),Nvar(3)
		integer :: cfg(latt%Ns,n_omp),nd(n_omp)
		complex(8) :: wf(latt%Ns*2,latt%Ns)
		real(8) :: er,rd
		type(t_phy) :: phy(n_omp)
		type(randomNumberSequence) :: rnd(n_omp)
		do i=1,n_omp
			allocate(phy(i)%csite(latt%Ns),phy(i)%ssite(latt%Ns))
		enddo
		Nmc(3)=1024*4
		Nvar(3)=128*4
		ne=latt%Ns
		ne2=ne/2
		do i=1,latt%Ns
			cfg(i,:)=i
		enddo
		do i=1,n_omp
			call fisher_yates_shuffle(cfg(:,i),latt%Ns)
		enddo
		nd=0

		open(40,file="../data/phyvar_input.dat",status="old",action="read",form="formatted")
		! calculation started
		!var(:)%val=(/ -1.39752E+00       , 8.88530E-01 , 1.32678E-02       , -2.00021E-01/)
		do i=0,100,2
		!do j=0,-100,-2
			!ne=latt%Ns*8/10
			ne=latt%Ns*7/8
			ne=ne+mod(ne,2)
			ne2=ne/2
			!var(1)%val=0.1d0*i-2d0
			!var(2)%val=10d0**(i/10d0-2d0)
			!var(3)%val=1d-1
			!var(2)%val=real(i)/50d0
			!var(1)%val=real(j)/50d0
			!var(2)%val=real(j)/50d0
			!var(1:)%val=(/-2d0,1d-1,1d0,1d-1,0.1d0/)
			!!do l=1,vn
				!call random_number(rd)
				!!var(l)%val=rd*var(l)%val
				!var(2)%val=rd
			!!enddo
			!call iniwf(wf,cp=var(1))
			!var(1:)%val=(/-1.12991E+00,1.04891E-01,5.62154E-01,-3.75703E-02,-3.84161E-01/)
			!var(1:)%val=(/&
				 !-1.12991d+00, 1.04891d-01, 5.62154d-01,-3.75703d-02,-3.84161d-01-t(2)&
				!/)
			!read(40,fmt="(i4,5e16.8)")ne,var(1:)%val
			!var(2)%val=var(2)%val/2d0
			Nmc(1:2)=(/50000,2**6/)
			Nvar(1:2)=(/10000,ne/)
			!var(1:)%val=(/-1.73E+00,1.87E+00,-4.51E-01,-6.08E-01,-1.84E-03/)
			!var(1:)%val=(/-2.69E-01,2.13E-01,1.35E-02,2.97E-02/)
			call variational(Nvar,cfg,nd)
			!call variational(Nmc,cfg,nd)
			write(*,"(i4$)")ne
			write(*,"(es9.2$)")var(1:)%val
			call iniwf(wf)
			!$OMP PARALLEL DO
			do k=1,n_omp
				call mt_init_random_seed(rnd(k),seed=k*50)
				call mc(rnd(k),wf,Nmc,cfg(:,k),nd(k),phy(k))
				!call mc(rnd(k),wf,Nvar,cfg(:,k),nd(k),phy(k))
				!call finalize_RandomNumberSequence(rnd(k))
			enddo
			!$OMP END PARALLEL DO
			er=real(phy(1)%E)**2
			do k=2,n_omp
				phy(1)%E=phy(1)%E+phy(k)%E
				phy(1)%af=phy(1)%af+phy(k)%af
				phy(1)%dsc=phy(1)%dsc+phy(k)%dsc
				phy(1)%ddw=phy(1)%ddw+phy(k)%ddw
				phy(1)%csite=phy(1)%csite+phy(k)%csite
				phy(1)%ssite=phy(1)%ssite+phy(k)%ssite
				er=er+real(phy(k)%E)**2
			enddo
			phy(1)%E=phy(1)%E/n_omp
			phy(1)%af=phy(1)%af/n_omp
			phy(1)%dsc=phy(1)%dsc/n_omp
			phy(1)%ddw=phy(1)%ddw/n_omp
			phy(1)%csite=phy(1)%csite/n_omp
			phy(1)%ssite=phy(1)%ssite/n_omp
			er=sqrt(abs(er/n_omp-real(phy(1)%E)**2))
			write(*,"(es11.4$)")real(phy(1)%E),er,phy(1)%af,phy(1)%dsc,phy(1)%ddw
			write(*,"(1x)")
			write(30,"(es13.5$)")real(latt%Ns-ne)/latt%Ns,var(1:)%val,real(phy(1)%E),er,phy(1)%af,phy(1)%dsc,phy(1)%ddw
			write(30,"(1x)")
			do l=1,vn
				if(var(l)%sg/=var(max(l-1,1))%sg.or.var(l)%nb/=var(max(l-1,1))%nb) then
					write(101,"(x/)")
				endif
				do k=1,size(var(l)%n)
					write(101,"(es13.2$)")latt%bond(var(l)%nb)%bd(var(l)%n(k))%r,latt%bond(var(l)%nb)%bd(var(l)%n(k))%dr,var(l)%bd_sg(k)*var(l)%val
					write(101,"(x)")
				enddo
			enddo
			write(101,"(x/)")
			do k=1,latt%Ns
				write(101,"(es13.2$)")latt%bond(0)%bd(k)%r,latt%bond(0)%bd(k)%dr,phy(1)%csite(k)
				write(101,"(x)")
			enddo
			write(101,"(x/)")
			do k=1,latt%Ns
				write(101,"(es13.2$)")latt%bond(0)%bd(k)%r,latt%bond(0)%bd(k)%dr,phy(1)%ssite(k)
				write(101,"(i5$)")ab(j)
				write(101,"(x)")
			enddo
			write(101,"(x/)")
			!V=V+2
			stop
		enddo
		!write(30,"(x)")
		!enddo
	end subroutine
end module

program main
	use M_vmc
	implicit none
	logical :: f
	! file
	f=openfile(101,"../data/lattice.dat")
	f=openfile(20,"../data/var.dat")
	!f=openfile(40,"../data/2d.dat")
	!f=openfile(10,"../data/err.dat")
	!f=openfile(30,"../data/phyvar.dat",access="direct",recl=226)
	f=openfile(30,"../data/phyvar.dat")

	call initial()
	call vmc_main()
end program
