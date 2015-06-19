module global
	use M_const
	use M_utility
	use M_lattice
	implicit none
	type t_sort
		real(8) :: val
		real(8), allocatable :: var(:)
		!real(8) :: var(4)
	end type
	type t_var
		real(8) :: val
		integer :: nb
		complex(8), allocatable :: bd_sg(:)
		integer, allocatable :: n(:)
		integer :: sg   ! 1:  pair
						! 2:  on-site charge
						! 5:  chemical potainial
						! 3:  on-site spin
						! 4:  bond charge
						! -1: n-n jast
	end type
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
	integer :: ne,ne2,vn,n_omp=24
	real(8) :: t(2)=(/1d0,-0.3d0/)
	real(8) :: DJ=0.3d0,V=0d0
	!real(8) :: DJ=1d0/3d0,V=1d0/12d0
	real(8) :: U=10d0
	type(t_var), allocatable :: var(:)
	logical :: no_sc=.true.
contains
	subroutine initial()
		integer :: i,l,k,Nb
		type(t_var) :: tmp(20)
		call init_random_seed()
		! lattice 
		a1=(/1d0,0d0/)
		a2=(/0d0,1d0/)
		T1=(/9d0,1d0/)
		T2=(/-1d0,9d0/)
		!T1=a1*10
		!T2=a2*10
		!T1=a1*10
		!T2=a2*10
		bdc(1)=1d0
		bdc(2)=1d0
		allocate(sub(1,2))
		sub(1,:)=(/0d0,0d0/)
		call gen_latt()
		call gen_neb()
		call gen_bond(5)
		! finish i2r(i,2),neb(i)%nb(j)%bond(k)/bdc(k)/r(k,2)
		write(*,*)"Total site number is: ",Ns
		do k=1,Ns
			if(is_a(k)) then
				!write(101,"(2es13.2,2I5)")i2r(k,:),k,1
			else
				!write(101,"(2es13.2,2I5)")i2r(k,:),k,-1
			endif
		enddo
		!write(101,"(1x/)")

		Nb=size(bond(1)%bd)
		i=0
		! cp
		i=i+1
		tmp(i)%sg=5
		tmp(i)%val=0d0
		allocate(tmp(i)%bd_sg(Ns))
		tmp(i)%bd_sg(:)=-1d0

		! d-wave sc
		i=i+1
		tmp(i)%sg=1
		tmp(i)%nb=1
		tmp(i)%val=1d-1
		allocate(tmp(i)%bd_sg(Nb))
		do k=1,Nb
			tmp(i)%bd_sg(k)=dwave(k)*&
				1d0
			!write(101,"(2es13.2,2es13.2)")bond(1)%bd(k)%r,tmp(i)%bd_sg(k)
		enddo
		!write(101,"(1x/)")

		!! stripe d-wave sc
		!i=i+1
		!tmp(i)%sg=1
		!tmp(i)%nb=1
		!tmp(i)%val=1d-1
		!allocate(tmp(i)%bd_sg(Nb))
		!do k=1,Nb
			!tmp(i)%bd_sg(k)=dwave(k)*&
				!cos(4d0*pi*0.125d0*(bond(1)%bd(k)%r(1)+0.5d0))
			!write(101,"(2es13.2,2es13.2)")bond(1)%bd(k)%r,tmp(i)%bd_sg(k)
		!enddo
		!write(101,"(1x/)")

		!! s-wave sc
		!i=i+1
		!tmp(i)%sg=1
		!tmp(i)%nb=1
		!tmp(i)%val=1d-01
		!allocate(tmp(i)%bd_sg(Nb))
		!do k=1,Nb
			!tmp(i)%bd_sg(k)=1d0
			!write(101,"(2es13.2,I5)")bond(1)%bd(k)%r,nint(dwave(k))
		!enddo
		!write(101,"(1x/)")


		! sdw
		i=i+1
		tmp(i)%sg=3
		tmp(i)%val=1d-1
		allocate(tmp(i)%bd_sg(Ns))
		do k=1,Ns
			if(is_a(k)) then
				tmp(i)%bd_sg(k)=&
					1d0
					!sin(2d0*pi*0.125d0*(i2r(k,1)+0.5d0))
			else
				tmp(i)%bd_sg(k)=&
					-1d0
					!-sin(2d0*pi*0.125d0*(i2r(k,1)+0.5d0))
			endif
			!write(101,"(2es13.2,2es13.2)")i2r(k,:),tmp(i)%bd_sg(k)
		enddo
		!write(101,"(1x/)")

		! ddw
		i=i+1
		tmp(i)%sg=4
		tmp(i)%nb=1
		tmp(i)%val=1d-1
		allocate(tmp(i)%bd_sg(Nb))
		do k=1,Nb
			if(is_a(bond(1)%bd(k)%i(1))) then
				tmp(i)%bd_sg(k)=img*dwave(k)
				!write(101,"(4es13.2)")bond(1)%bd(k)%r-0.5d0*bond(1)%bd(k)%dir(:)*dwave(k),bond(1)%bd(k)%dir(:)*dwave(k)
			else
				tmp(i)%bd_sg(k)=-img*dwave(k)
				!write(101,"(4es13.2)")bond(1)%bd(k)%r+0.5d0*bond(1)%bd(k)%dir(:)*dwave(k),-bond(1)%bd(k)%dir(:)*dwave(k)
			endif
		enddo
		!write(101,"(1x/)")

		!! d-wave cdw
		!i=i+1
		!tmp(i)%sg=4
		!tmp(i)%nb=1
		!tmp(i)%val=1d-1
		!allocate(tmp(i)%bd_sg(Nb))
		!do k=1,size(bond(1)%bd)
			!tmp(i)%bd_sg(k)=dwave(k)
		!enddo

		!! on site cdw
		!i=i+1
		!tmp(i)%sg=2
		!tmp(i)%val=1d-1
		!allocate(tmp(i)%bd_sg(Ns))
		!do k=1,Ns
			!tmp(i)%bd_sg(k)=cos(4d0*pi*0.125d0*(i2r(k,1)+0.5d0))
			!write(101,"(2es13.2,2es13.2)")i2r(k,:),tmp(i)%bd_sg(k)
		!enddo
		!write(101,"(1x/)")

		! hp
		i=i+1
		tmp(i)%sg=4
		tmp(i)%nb=2
		tmp(i)%val=0d0
		allocate(tmp(i)%bd_sg(size(bond(2)%bd)))
		tmp(i)%bd_sg(:)=-1d0

		i=i+1
		tmp(i)%sg=4
		tmp(i)%nb=3
		tmp(i)%val=0d0
		allocate(tmp(i)%bd_sg(size(bond(3)%bd)))
		tmp(i)%bd_sg(:)=-1d0

		vn=i
		!! n-n jast
		!i=i+1
		!tmp(i)%sg=-1
		!tmp(i)%val=0d0
		!!allocate(tmp(i)%bd_sg(Ns))
		!i=i+1
		!tmp(i)%sg=-1
		!tmp(i)%val=0d0
		!allocate(tmp(i)%bd_sg(Nb))
		!i=i+1
		!tmp(i)%sg=-1
		!tmp(i)%val=0d0
		!!allocate(tmp(i)%bd_sg(Nb))
		!i=i+1
		!tmp(i)%sg=-1
		!tmp(i)%val=0d0
		!!allocate(tmp(i)%bd_sg(Nb))
		!i=i+1
		!tmp(i)%sg=-1
		!tmp(i)%val=0d0
		!!allocate(tmp(i)%bd_sg(Nb))
		!i=i+1
		!tmp(i)%sg=-1
		!tmp(i)%val=0d0
		!!allocate(tmp(i)%bd_sg(Nb))

		allocate(var(i))
		do l=1,i
			allocate(var(l)%bd_sg(size(tmp(l)%bd_sg)))
			if(allocated(tmp(l)%n)) then
				allocate(var(l)%n(size(tmp(l)%n)))
			endif
		enddo
		var=tmp(:i)
	end subroutine
	function is_a(i)
		logical :: is_a
		integer :: i
		if(mod(sum(nint(i2r(i,:)-i2r(1,:))),2)==0) then
			is_a=.true.
		else
			is_a=.false.
		endif
	end function
	function dwave(i)
		real(8) :: dwave
		integer :: i
		if(nint(bond(1)%bd(i)%dir(2))==0) then
			dwave=1d0
		else
			dwave=-1d0
		endif
	end function
end module

module M_wf
	use global
	use M_matrix
	use lapack95, only: heevd,heevr
	implicit none
contains
	subroutine iniwf(wf,dwf,cp)
		complex(8) :: H(Ns*2,Ns*2),D(Ns*2,Ns*2),wf(:,:),tmp(Ns),Q(Ns,Ns)
		complex(8), optional :: dwf(:,:,:)
		real(8) :: E(Ns*2),r,cp1
		type(t_var), optional :: cp
		complex(8) :: bd
		integer :: i,j,k,l,n,info,sg
		if(present(cp).and.cp%sg/=5) then
			return
		endif
		H=0d0
		do k=1,size(bond(1)%bd)
			i=bond(1)%bd(k)%i(1)
			j=bond(1)%bd(k)%i(2)
			bd=bond(1)%bd(k)%bdc
			H(i,j)=H(i,j)-1d0*conjg(bd)
			H(j,i)=H(j,i)-1d0*bd
			H(i+Ns,j+Ns)=H(i+Ns,j+Ns)+1d0*bd
			H(j+Ns,i+Ns)=H(j+Ns,i+Ns)+1d0*conjg(bd)
		enddo
		do l=1,size(var)
			select case(var(l)%sg)
			case(1)
				! pair channel
				no_sc=.false.
				do k=1,size(var(l)%bd_sg)
					i=bond(var(l)%nb)%bd(k)%i(1)
					j=bond(var(l)%nb)%bd(k)%i(2)
					bd=bond(var(l)%nb)%bd(k)%bdc
					H(i,j+Ns)=H(i,j+Ns)+var(l)%val*var(l)%bd_sg(k)*bd
					H(j,i+Ns)=H(j,i+Ns)+var(l)%val*var(l)%bd_sg(k)*bd
					H(i+Ns,j)=H(i+Ns,j)+var(l)%val*var(l)%bd_sg(k)*conjg(bd)
					H(j+Ns,i)=H(j+Ns,i)+var(l)%val*var(l)%bd_sg(k)*conjg(bd)
				enddo
			case(2,5)
				!on site
				do k=1,size(var(l)%bd_sg)
					i=k
					H(i,i)=H(i,i)+var(l)%val*var(l)%bd_sg(k)
					H(i+Ns,i+Ns)=H(i+Ns,i+Ns)-var(l)%val*var(l)%bd_sg(k)
				enddo
			case(3)
				!on site spin channel
				do k=1,size(var(l)%bd_sg)
					i=k
					H(i,i)=H(i,i)+var(l)%val*var(l)%bd_sg(k)
					H(i+Ns,i+Ns)=H(i+Ns,i+Ns)+var(l)%val*var(l)%bd_sg(k)
				enddo
			case(4)
				!bond charge channel
				do k=1,size(var(l)%bd_sg)
					i=bond(var(l)%nb)%bd(k)%i(1)
					j=bond(var(l)%nb)%bd(k)%i(2)
					bd=bond(var(l)%nb)%bd(k)%bdc
					H(i,j)=H(i,j)+var(l)%val*conjg(var(l)%bd_sg(k))*bd
					H(j,i)=H(j,i)+conjg(var(l)%val*conjg(var(l)%bd_sg(k))*bd)
					H(i+Ns,j+Ns)=H(i+Ns,j+Ns)-var(l)%val*var(l)%bd_sg(k)*bd
					H(j+Ns,i+Ns)=H(j+Ns,i+Ns)-conjg(var(l)%val*var(l)%bd_sg(k)*bd)
				enddo
			end select
		enddo

		if(present(cp)) then
			write(*,"(A$)")"*"
			call heevd(H(:Ns,:Ns),E(:Ns),"N","U",info)
			cp%val=cp%val+(E(ne2)+E(ne2+1))/2d0
			return
		endif
		if(no_sc) then
			do k=1,size(bond(1)%bd)
				i=bond(1)%bd(k)%i(1)
				j=bond(1)%bd(k)%i(2)
				bd=bond(1)%bd(k)%bdc
				H(i,j+Ns)=H(i,j+Ns)+dwave(k)*1d-4*bd
				H(j,i+Ns)=H(j,i+Ns)+dwave(k)*1d-4*bd
				H(i+Ns,j)=H(i+Ns,j)+dwave(k)*1d-4*conjg(bd)
				H(j+Ns,i)=H(j+Ns,i)+dwave(k)*1d-4*conjg(bd)
			enddo
		endif
		if(var(1)%sg/=5) then
			!Q=H(:Ns,:Ns)
			!call heevd(Q,E(:Ns),"N","U",info)
			!H(:Ns,:Ns)=H(:Ns,:Ns)-diag((E(ne2)+E(ne2+1))/2d0,Ns)
			!H(Ns+1:,Ns+1:)=H(Ns+1:,Ns+1:)+diag((E(ne2)+E(ne2+1))/2d0,Ns)

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
		!call mwrite(101,H(:,:))
		!stop
		if(present(dwf)) then
			dwf=0d0
			do l=1,vn
				D=0d0
				Q=0d0
				select case(var(l)%sg)
				case(1)
					! pair channel
					do k=1,size(var(l)%bd_sg)
						i=bond(var(l)%nb)%bd(k)%i(1)
						j=bond(var(l)%nb)%bd(k)%i(2)
						bd=bond(var(l)%nb)%bd(k)%bdc
						D(i,:)=D(i,:)+H(j+Ns,:)*var(l)%bd_sg(k)*bd
						D(j,:)=D(j,:)+H(i+Ns,:)*var(l)%bd_sg(k)*bd
						D(i+Ns,:)=D(i+Ns,:)+H(j,:)*var(l)%bd_sg(k)*conjg(bd)
						D(j+Ns,:)=D(j+Ns,:)+H(i,:)*var(l)%bd_sg(k)*conjg(bd)
					enddo
				case(2,5)
					!on site
					do k=1,size(var(l)%bd_sg)
						i=k
						D(i,:)=D(i,:)+H(i,:)*var(l)%bd_sg(k)
						D(i+Ns,:)=D(i+Ns,:)-H(i+Ns,:)*var(l)%bd_sg(k)
					enddo
				case(3)
					!on site spin channel
					do k=1,size(var(l)%bd_sg)
						i=k
						D(i,:)=D(i,:)+H(i,:)*var(l)%bd_sg(k)
						D(i+Ns,:)=D(i+Ns,:)+H(i+Ns,:)*var(l)%bd_sg(k)
					enddo
				case(4)
					!bond charge channel
					do k=1,size(var(l)%bd_sg)
						i=bond(var(l)%nb)%bd(k)%i(1)
						j=bond(var(l)%nb)%bd(k)%i(2)
						bd=bond(var(l)%nb)%bd(k)%bdc
						D(i,:)=D(i,:)+H(j,:)*conjg(var(l)%bd_sg(k))*bd
						D(j,:)=D(j,:)+H(i,:)*conjg(conjg(var(l)%bd_sg(k))*bd)
						D(i+Ns,:)=D(i+Ns,:)-H(j+Ns,:)*var(l)%bd_sg(k)*bd
						D(j+Ns,:)=D(j+Ns,:)-H(i+Ns,:)*conjg(var(l)%bd_sg(k)*bd)
					enddo
				end select
				!$OMP PARALLEL DO COLLAPSE(2)
				do i=1,Ns
					do j=1,Ns
						if(abs(E(j)-E(i+Ns))<1d-8) then
							cycle
						endif
						Q(i,j)=Q(i,j)+dot_product(H(:,i+Ns),D(:,j))/(E(j)-E(i+Ns))
					enddo
				enddo
				!$OMP END PARALLEL DO
				!$OMP PARALLEL DO PRIVATE(tmp)
				do i=1,Ns*2
					tmp=matmul(H(i,Ns+1:),Q)
					do j=1,Ns*2
						dwf(i,j,l)=dwf(i,j,l)+sum(tmp*conjg(H(j,:Ns)))
					enddo
				enddo
				!$OMP END PARALLEL DO
			enddo
		endif
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
			m(1)=cfg(i)+Ns
		case(6,8) 
			k(1)=j+ne2
			m(1)=cfg(i)+Ns
		case(0)
			k=(/i,i+ne2/)
			m=(/cfg(j),cfg(j)+Ns/)
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
	subroutine update(D,k,m,sg)
		complex(8) :: D(:,:),tmp(size(D,1),size(D,2)),ipb(2,2),pb,tmp2(2)
		integer :: i,j,sg,m(:),k(:)
		tmp=D
		select case(sg) 
		case(0)
			pb=1d0/(D(m(1),k(1))*D(m(2),k(2))-D(m(1),k(2))*D(m(2),k(1)))
			ipb(1,1)=pb*D(m(2),k(2))
			ipb(1,2)=-pb*D(m(1),k(2))
			ipb(2,1)=-pb*D(m(2),k(1))
			ipb(2,2)=pb*D(m(1),k(1))
			do i=1,size(D,1)
				tmp2=matmul((/tmp(i,k(1)),tmp(i,k(2))/),ipb)
				do j=1,size(D,2)
					if(k(1)==j) then
						D(i,j)=tmp(i,j)-sum(tmp2*(/tmp(m(1),j)-1d0,tmp(m(2),j)/))
					elseif(k(2)==j) then
						D(i,j)=tmp(i,j)-sum(tmp2*(/tmp(m(1),j),tmp(m(2),j)-1d0/))
					else
						D(i,j)=tmp(i,j)-sum(tmp2*(/tmp(m(1),j),tmp(m(2),j)/))
					endif
				enddo
			enddo
		case(1:8)
			pb=1d0/D(m(1),k(1))
			do i=1,size(D,1)
				do j=1,size(D,2)
					if(j==k(1)) then
						D(i,j)=tmp(i,j)-tmp(i,k(1))*pb*(tmp(m(1),j)-1d0)
					else
						D(i,j)=tmp(i,j)-tmp(i,k(1))*pb*(tmp(m(1),j))
					endif
				enddo
			enddo
		end select
	end subroutine
	subroutine update_swap(i,j,nd,cfg,icfg,D,sg)
		complex(8) :: D(:,:)
		integer :: i,j,cfg(:),icfg(:),sg,nd
		select case(sg)
		case(0,2,5)
			call swap(cfg(i),cfg(j))
			call swap(icfg(cfg(i)),icfg(cfg(j)))
		case(1)
			call swap(cfg(i),cfg(j))
			call swap(icfg(cfg(i)),icfg(cfg(j)))
			call swap(D(:,i+ne2),D(:,j))
		case(3)
			call swap(cfg(i),cfg(ne2-nd+1))
			call swap(icfg(cfg(i)),icfg(cfg(ne2-nd+1)))
			call swap(cfg(j),cfg(ne-nd+1))
			call swap(icfg(cfg(j)),icfg(cfg(ne-nd+1)))
			call swap(cfg(ne2-nd+1),cfg(ne-nd+1))
			call swap(icfg(cfg(ne2-nd+1)),icfg(cfg(ne-nd+1)))
			call swap(D(:,i),D(:,ne2-nd+1))
			call swap(D(:,j),D(:,ne-nd+1))
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
			nd=nd+1
		case(6)
			call swap(cfg(i),cfg(j))
			call swap(icfg(cfg(i)),icfg(cfg(j)))
			call swap(D(:,i),D(:,j))
		case(7)
			call swap(cfg(i),cfg(ne2-nd+1))
			call swap(icfg(cfg(i)),icfg(cfg(ne2-nd+1)))
			call swap(cfg(j),cfg(ne-nd+1))
			call swap(icfg(cfg(j)),icfg(cfg(ne-nd+1)))
			call swap(D(:,i),D(:,ne2-nd+1))
			call swap(D(:,j),D(:,ne-nd+1))
			nd=nd-1
		case(8)
			call swap(cfg(j),cfg(ne2-nd))
			call swap(icfg(cfg(j)),icfg(cfg(ne2-nd)))
			call swap(cfg(i),cfg(ne-nd))
			call swap(icfg(cfg(i)),icfg(cfg(ne-nd)))
			call swap(D(:,j),D(:,ne2-nd))
			call swap(D(:,j+ne2),D(:,ne-nd))
			nd=nd+1
		end select
	end subroutine
	function detO(cfg,A,D)
		complex(8) :: A(:,:),D(:,:),detO,pb
		integer :: i,j,k(2),m(2),sg,cfg(:)
		detO=0d0
		do i=1,Ns
			do j=1,Ns*2
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
			cs=3
		else
			cs=4
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
		if(ci==1.and.cj==4) then
			sg=1
		elseif(ci==2.and.cj==3.and.ud==1) then
			sg=2
		elseif(ci==2.and.cj==4) then
			select case(ud)
			case(1)
				sg=3
			case(-1)
				sg=7
			end select
		elseif(ci==1.and.cj==3) then
			select case(ud)
			case(1)
				sg=4
			case(-1)
				sg=8
			case(0)
				sg=0
			end select
		elseif(ci==3.and.cj==4) then
			sg=5
		elseif(ci==1.and.cj==2.and.ud==-1) then
			sg=6
		else
			sg=-1
		endif
	end subroutine
	subroutine ini_nn(cfg,icfg,nd,nn)
		integer :: ni,nj,icfg(:),cfg(:),i,j,k,l,nn(:),nd
		if(size(nn)==0) then
			return
		endif
		nn=0
		do k=1,Ns
			select case(cs(icfg(k),nd))
			case(1,3)
				nn(1)=nn(1)+1
			case(2)
				nn(1)=nn(1)+4
			end select
		enddo
		do l=2,size(nn)
			do k=1,size(bond(l-1)%bd)
				i=bond(l-1)%bd(k)%i(1)
				j=bond(l-1)%bd(k)%i(2)
				select case(cs(icfg(i),nd))
				case(1,3)
					ni=1
				case(2)
					ni=2
				case default
					ni=0
				end select
				select case(cs(icfg(j),nd))
				case(1,3)
					nj=1
				case(2)
					nj=2
				case default
					nj=0
				end select
				nn(l)=nn(l)+ni*nj
			enddo
		enddo
	end subroutine
	subroutine nn_dif(ri,rj,cfg,icfg,sg,nd,nndif)
		integer :: ri,rj,icfg(:),cfg(:),sg,i,j,k,l,nndif(:),nd
		if(size(nndif)==0) then
			return
		endif
		nndif=0
		select case(sg)
		case(1:8)
			select case(cs(icfg(ri),nd))
			case(1,3)
				nndif(1)=nndif(1)-1
			case(2)
				nndif(1)=nndif(1)-3
			end select
			select case(cs(icfg(rj),nd))
			case(1,3)
				nndif(1)=nndif(1)+3
			case(4)
				nndif(1)=nndif(1)+1
			end select
			do l=2,size(nndif)
				do k=1,size(neb(ri)%nb(l)%neb)
					i=neb(ri)%nb(l)%neb(k)
					select case(cs(icfg(i),nd))
					case(1,3)
						nndif(l)=nndif(l)-1
					case(2)
						nndif(l)=nndif(l)-2
					end select
					if(i==rj) then
						nndif(l)=nndif(l)-1
					endif
				enddo
				do k=1,size(neb(rj)%nb(l)%neb)
					j=neb(rj)%nb(l)%neb(k)
					select case(cs(icfg(j),nd))
					case(1,3)
						nndif(l)=nndif(l)+1
					case(2)
						nndif(l)=nndif(l)+2
					end select
				enddo
			enddo
		end select
	end subroutine
	function jast(nn_dif,ja)
		integer :: nn_dif(:)
		real(8) :: jast,ja(:)
		if(size(ja)==0) then
			jast=1d0
			return
		endif
		jast=exp(-sum(ja*nn_dif(:size(ja))))
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
		n0=Ns-ne+nd
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
	subroutine energy(cfg,icfg,nd,nn,D,ja,El)
		complex(8) :: pb,D(:,:)
		complex(8) :: El
		real(8) :: ja(:)
		integer :: cfg(:),icfg(:),i,j,n,l,k(2),m(2),sg,nd,nn(:),nndif(size(nn))
		El=0d0
		do i=1,ne
			do l=1,size(t)
				do n=1,size(neb(cfg(i))%nb(l+1)%neb)
					j=icfg(neb(cfg(i))%nb(l+1)%neb(n))
					call get_case(i,j,0,nd,sg)
					select case(sg) 
					case(1) 
						call nn_dif(cfg(i),cfg(j),cfg,icfg,sg,nd,nndif)
						call getindex(i,j,cfg,sg,k,m)
						call getpb(k,m,sg,D,pb)
						El=El-(t(l)*conjg(pb))&
							*jast(nndif,ja)&
							*neb(cfg(i))%nb(l+1)%bdc(n)
					case(5) 
						call nn_dif(cfg(i),cfg(j),cfg,icfg,sg,nd,nndif)
						call getindex(i,j,cfg,sg,k,m)
						call getpb(k,m,sg,D,pb)
						El=El+(t(l)*conjg(pb))&
							*jast(nndif,ja)&
							*neb(cfg(i))%nb(l+1)%bdc(n)
					case(0)
						if(l==1.and.i<=ne2) then
							call getindex(i,j,cfg,sg,k,m)
							call getpb(k,m,sg,D,pb)
							El=El+0.5d0*DJ*conjg(pb)-0.25d0*DJ
						endif
					case(-1)
						if(l==1.and.n<3) then
							El=El+0.25d0*DJ
						endif
					end select
				enddo
			enddo
		enddo
		El=El-0.25d0*DJ*nn(2)
		El=El+V*nn(2)
		!write(*,*)El
		!El=(El-(4d0*V-DJ)*(ne-Ns/2))/Ns
		!El=(El-(4d0*V)*(ne-Ns/2))/Ns
		El=El/Ns
	end subroutine
	subroutine spin_order(cfg,nd,lphy,af)
		real(8) :: lphy(:),af
		integer :: cfg(:),i,nd
		lphy=0d0
		af=0d0
		do i=1,ne2-nd
			lphy(cfg(i))=lphy(cfg(i))+0.5d0
			if(is_a(cfg(i))) then
				af=af+0.5d0
			else
				af=af-0.5d0
			endif
		enddo
		do i=ne2+1,ne-nd
			lphy(cfg(i))=lphy(cfg(i))-0.5d0
			if(is_a(cfg(i))) then
				af=af-0.5d0
			else
				af=af+0.5d0
			endif
		enddo
		af=af/Ns
	end subroutine
	subroutine charge_order(cfg,nd,lphy)
		real(8) :: lphy(:)
		integer :: cfg(:),i,nd
		lphy=0d0
		do i=1,ne2
			lphy(cfg(i))=lphy(cfg(i))+1d0
		enddo
		do i=ne2-nd+1,ne-nd
			lphy(cfg(i))=lphy(cfg(i))+1d0
		enddo
	end subroutine
	subroutine ddw_order(cfg,icfg,nd,D,ja,lphy)
		complex(8) :: D(:,:),pb
		real(8) :: lphy,ja(:)
		integer :: cfg(:),icfg(:),b1,j1,i1,sg,nd,k(2),m(2)
		lphy=0d0
		do b1=1,size(bond(1)%bd)
			i1=bond(1)%bd(b1)%i(1)
			j1=bond(1)%bd(b1)%i(2)
			if(icfg(i1)<=ne2.and.icfg(j1)>ne) then
				k(1)=icfg(i1)
				m(1)=j1
				sg=1
			elseif(icfg(j1)<=ne2.and.icfg(i1)>ne) then
				k(1)=icfg(j1)
				m(1)=i1
				sg=-1
			elseif(icfg(i1)>ne2.and.icfg(i1)<=ne.and.icfg(j1)>ne) then
				k(1)=icfg(j1)
				m(1)=i1+Ns
				sg=-1
			elseif(icfg(j1)>ne2.and.icfg(j1)<=ne.and.icfg(i1)>ne) then
				k(1)=icfg(i1)
				m(1)=j1+Ns
				sg=1
			else
				cycle
			endif
			call getpb(k,m,1,D,pb)
			if(is_a(i1)) then
				lphy=lphy+imag(pb)*sg*dwave(b1)
			else                               
				lphy=lphy-imag(pb)*sg*dwave(b1)
			endif
		enddo
		lphy=lphy/Ns
	end subroutine
	subroutine dsc_corelation(cfg,icfg,nd,D,ja,lphy)
		complex(8) :: D(:,:),pb
		real(8) :: lphy,ja(:)
		integer :: cfg(:),icfg(:),i1,j1,i2,j2,nd,k(2),m(2),b1,b2
		lphy=0d0
		do b1=1,size(bond(1)%bd)
			i1=bond(1)%bd(b1)%i(1)
			j1=bond(1)%bd(b1)%i(2)
			if(icfg(i1)<=ne2.and.icfg(j1)>ne2.and.icfg(j1)<=ne) then
				k(1)=icfg(i1)
				m(2)=j1+Ns
			elseif(icfg(j1)<=ne2.and.icfg(i1)>ne2.and.icfg(i1)<=ne) then
				k(1)=icfg(j1)
				m(2)=i1+Ns
			else
				cycle
			endif
			do b2=1,size(bond(1)%bd)
				i2=bond(1)%bd(b2)%i(1)
				j2=bond(1)%bd(b2)%i(2)
				if(icfg(i2)>ne.and.icfg(j2)>ne) then
					m(1)=i2
					k(2)=icfg(j2)
				else
					cycle
				endif
				call getpb(k,m,0,D,pb)
				lphy=lphy+real(pb)*dwave(b1)*dwave(b2)
				m(1)=j2
				k(2)=icfg(i2)
				call getpb(k,m,0,D,pb)
				lphy=lphy+real(pb)*dwave(b1)*dwave(b2)
			enddo
		enddo
		lphy=lphy/(Ns*Ns)
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
		n0=Ns-ne+nd
		n=mt_random_number(rnd,ne2*(Ns-ne2)*2)
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
		elseif(n<=ne2*(Ns-ne2)) then
			n=n-n1*n0-nd*n1-nd*n0
			sg=4
			i=(n-1)/n1+1
			j=ne2+mod(n-1,n1)+1
		elseif(n<=n1*n0+ne2*(Ns-ne2)) then
			n=n-ne2*(Ns-ne2)
			sg=5
			i=ne2+(n-1)/n0+1
			j=ne-nd+mod(n-1,n0)+1
		elseif(n<=n1*n0+nd*n1+ne2*(Ns-ne2)) then
			n=n-ne2*(Ns-ne2)-n1*n0
			sg=6
			j=(n-1)/nd+1
			i=n1+mod(n-1,nd)+1
		elseif(n<=n1*n0+nd*n1+nd*n0+ne2*(Ns-ne2)) then
			n=n-ne2*(Ns-ne2)-n1*n0-nd*n1
			sg=7
			i=n1+(n-1)/n0+1
			j=ne-nd+mod(n-1,n0)+1
		elseif(n<=ne2*(Ns-ne2)*2) then
			n=n-ne2*(Ns-ne2)-n1*n0-nd*n1-nd*n0
			sg=8
			j=mod(n-1,n1)+1
			i=ne2+(n-1)/n1+1
		endif
	end subroutine
	!subroutine energy(cfg,icfg,nd,nn,wf,ja,A,iA,El)
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
				!do k=1,size(neb(cfg(i))%nb(l+1)%neb)
					!j=icfg(neb(cfg(i))%nb(l+1)%neb(k))
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
		!El=El/Ns
		!!stop
	!end subroutine
	!subroutine spin_order(cfg,nd,lphy)
		!real(8) :: lphy
		!integer :: cfg(:),i,x1(2),nd
		!lphy=0d0
		!do i=1,ne-nd*2
			!call latt_one2two(cfg(i),Ns,x1)
			!!lphy=lphy+sign(-1d0,i-ne2-0.1d0)*(0.5d0-mod(x1(1)+x1(2),2))
			!lphy=lphy+(1-(i-1)/(ne2-nd)*2)*(0.5d0-mod(sum(x1),2))
		!enddo
		!lphy=lphy/Ns
	!end subroutine
	!!subroutine charge_order(cfg,icfg,nd,wf,A,iA,lphy)
		!!complex(8) :: A(:,:),iA(:,:),vu(Ns,2),Y(2,2),pb,wf(:,:),a11
		!!complex(8) :: lphy
		!!integer :: cfg(:),icfg(:),i,n,m,j,k,l,inb,sg,cr(2),x1(2),x2(2),nr(2),dr(2),r(2),nd
		!!lphy=0d0
		!!lphy=lphy/Ns
	!!end subroutine
	!subroutine ddw_order(cfg,icfg,nd,wf,A,iA,ja,lphy)
		!complex(8) :: A(:,:),iA(:,:),vu(Ns,2),Y(2,2),pb,wf(:,:),a11
		!real(8) :: lphy,ja(:)
		!integer :: cfg(:),icfg(:),i,n,m,j,k,l,inb,sg,cr(2),x1(2),nd,ud
		!lphy=0d0
		!do n=1,ne
			!ud=1-(n-1)/ne2*2
			!i=n-nd*(1-ud)/2
			!call latt_one2two(cfg(i),Ns,x1)
			!do inb=1,4
				!j=icfg(neb(cfg(i),inb,1))
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
		!lphy=lphy/Ns
	!end subroutine
	!!subroutine dsc_order(cfg,icfg,nd,wf,A,iA,ja,lphy)
		!!complex(8) :: A(:,:),iA(:,:),vu(Ns,2),Y(2,2),pb,wf(:,:),a11
		!!real(8) :: lphy,ja(:)
		!!integer :: cfg(:),icfg(:),i,n,m,j,k,l,inb,sg,cr(2),x1(2),x2(2),nr(2),dr(2),r(2),nd
		!!do i=ne2-nd+1,ne-nd
			!!do inb=1,4,1
				!!j=icfg(neb(cfg(i),inb,1))
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
		!!complex(8) :: A(:,:),iA(:,:),vu(Ns,2),Y(2,2),pb,wf(:,:),a11
		!!real(8) :: lphy
		!!integer :: cfg(:),icfg(:),i,n,m,j,k,l,inb,sg,cr(2),x1(2),x2(2),nr(2),dr(2),r(2),nd
	!!end subroutine
	!subroutine dsc_corelation(cfg,icfg,nd,wf,A,iA,ja,lphy)
		!complex(8) :: A(:,:),iA(:,:),vu(Ns,2),Y(2,2),pb,wf(:,:),a11
		!real(8) :: lphy,ja(:)
		!integer :: cfg(:),icfg(:),i,j,ii,jj,inb1,inb2,cr(2),x1(2),x2(2),nr(2),dr(2),r(2),nd
		!lphy=0d0
		!do i=ne2-nd+1,ne-nd
			!do inb1=1,4,1
				!j=icfg(neb(cfg(i),inb1,1))
				!if(j<=ne2) then
					!do jj=ne2-nd+1,Ns
						!do inb2=1,4,1
							!ii=icfg(neb(cfg(jj),inb2,1))
							!if(ii<=(ne2-nd)) then
								!cr=(/ii+ne2,j/)
							!elseif(ii>(ne-nd)) then
								!cr=(/ii,j/)
							!else
								!cycle
							!endif
							!vu(:,1)=wf(cfg(i)+Ns,:)
							!vu(:,2)=wf(cfg(jj),:)
							!call det(vu,cr,A,iA,Y,pb)
							!lphy=lphy+real(conjg(pb)*jast(cfg(i),cfg(ii),cfg,icfg,1,nd,ja)*jast(cfg(j),cfg(jj),cfg,icfg,1,nd,ja))*&
								!(0.5d0-mod(inb1+inb2,2))
						!enddo
					!enddo
				!endif
			!enddo
		!enddo
		!lphy=lphy/(Ns*Ns)
	!end subroutine
end module

module M_vmc
	use M_tJ
	!use M_hubbard
	use M_wf
	implicit none
contains
	recursive subroutine qsort(A)
		!From: http://rosettacode.org/wiki/Sorting_algorithms/Quicksort#Fortran
		type(t_sort) :: A(:),tmp
		integer :: left,right,n
		real(8) :: random
		real(8) :: pivot
		integer :: marker
		n=size(A)
		if(n>1) then
			call random_number(random)
			pivot=A(int(random*real(n-1))+1)%val   !random pivor (not best performance, but avoids worst-case)
			left=0
			right=n+1
			do while (left < right)
				right=right-1
				do while(A(right)%val>pivot)
					right=right-1
				enddo
				left=left+1
				do while(A(left)%val<pivot)
					left= left+1
				enddo
				if(left<right) then
					tmp=A(left)
					A(left)=A(right)
					A(right)=tmp
				endif
			end do
			if(left==right) then
				marker=left+1
			else
				marker=left
			endif
			call qsort(A(1:marker-1))
			call qsort(A(marker:n))
		endif
	end subroutine
	subroutine mc(rnd,wf,Nmc,cfg,nd,phy,dwf,S,g)
		complex(8) :: pb,A(Ns,Ns),iA(Ns,Ns),wf(:,:),D(size(wf,1),size(wf,2))
		complex(8), optional :: S(:,:),dwf(:,:,:)
		real(8), optional :: g(:)
		real(8) :: rpb
		real(8) :: ja(size(var)-vn)
		complex(8) :: lO(size(var)),lS(size(lO),size(lO)),lg(size(lO)),O(size(lO))
		complex(8) :: El,E
		type(t_phy) :: phy,lphy
		integer :: icfg(Ns),i,j,ti,tj,l,sg,n,Nmc(:),k(2),m(2),acp,info
		integer :: cfg(Ns),nd,nn(max(size(var)-vn,2)),nndif(size(nn))
		logical :: flag
		type(randomNumberSequence) :: rnd

		if(mod(Ns,2)/=0) stop "Number of site is not even, exit!!"
		if(ne/=ne2*2) stop "Number of electrons is not even, exit!!"

		flag=present(S)

		n=1
		!get initial matrix
		do 
			do i=1,Ns
				if(i<=ne2) then
					A(i,:)=wf(cfg(i),:)
				elseif(i>=ne-nd+1)  then
					A(i,:)=wf(cfg(i)+Ns,:)
				else
					A(i,:)=wf(cfg(i-ne2)+Ns,:)
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

		do i=1,Ns
			icfg(cfg(i))=i
		enddo

		ja=var(vn+1:)%val



		!do i=1,Ns
			!cfg(i)=i
		!enddo
		!nd=0
		!call fisher_yates_shuffle(cfg,Ns)

		n=0
		acp=0
		call ini_nn(cfg,icfg,nd,nn)

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
		endif

		do
			n=n+1
			call two(rnd,i,j,nd,sg)
			!write(*,*)i,j,nd,sg
			call getindex(i,j,cfg,sg,k,m)
			call getpb(k,m,sg,D,pb)
			call nn_dif(cfg(i),cfg(j),cfg,icfg,sg,nd,nndif)
			!call random_number(rpb)
			rpb=getrandomreal(rnd)
			if(rpb<real(pb*conjg(pb))*(jast(nndif,ja))**2) then
				acp=acp+1
				call update(D,k,m,sg)
				!write(*,*)sum(abs(matmul(A,iA)-diag(1d0,size(A,1))))
				call update_swap(i,j,nd,cfg,icfg,D,sg)
				nn=nn+nndif
			endif
			!$OMP MASTER
			if(n==Nmc(1)) write(*,"(A$)")"#"
			!$OMP END MASTER
			if(n>Nmc(1).and.mod(n-Nmc(1),Nmc(2))==0) then
				call energy(cfg,icfg,nd,nn,D,ja,lphy%E)
				if(flag) then
					do ti=1,vn
						lO(ti)=detO(cfg,dwf(:,:,ti),D)
					enddo
					if(size(ja)/=0) then
						lO(vn+1:)=-nn
					endif
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
						call spin_order(cfg,nd,lphy%ssite,lphy%af)
						call charge_order(cfg,nd,lphy%csite)
						call ddw_order(cfg,icfg,nd,D,ja,lphy%ddw)
						call dsc_corelation(cfg,icfg,nd,D,ja,lphy%dsc)
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
					g=2d0*(g-real(phy%E*O))*Ns
				endif
				exit
			endif
		enddo
		!write(*,"('accept/total number is ',I6,'/',I6)")acp,n
	end subroutine
	subroutine variational(Nmc,cfg,nd)
		use lapack95, only: heevd,heevr
		real(8) :: dvar(size(var)),eg(size(var)),dt,scv
		complex(8) :: wf(Ns*2,Ns),dwf(Ns*2,Ns*2,vn)
		complex(8) :: S(size(var),size(var),n_omp)
		real(8) :: g(size(var),n_omp),cg(size(var)),h(size(var))
		integer :: m,i,j,k,info,Nmc(:),l,ml
		integer :: cfg(:,:),nd(:)
		real(8) :: er
		type(randomNumberSequence) :: rnd(n_omp)
		type(t_phy) :: phy(n_omp)
		type(t_sort), allocatable :: allE(:)
		dt=0.1d0
		m=50
		l=1
		scv=0d0
		allocate(allE(600))
		cg=1d0
		h=0d0
		ml=1
		do
			write(*,"(I3$)")ne
			write(*,"(es9.2$)")var(:)%val
			call iniwf(wf,dwf)
			!$OMP PARALLEL DO
			do k=1,n_omp
				call mt_init_random_seed(rnd(k),seed=k*10)
				call mc(rnd(k),wf,Nmc,cfg(:,k),nd(k),phy(k),dwf,S(:,:,k),g(:,k))
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
			S(:,:,1)=S(:,:,1)+diag(1d-3,size(S(:,:,1),1))
			call heevd(S(:,:,1),eg,"V")
			do i=1,size(g(:,1))
				dvar(i)=0d0
				do j=1,size(g(:,1))
					do k=1,size(g(:,1))
						!if(eg(k)/eg(size(g))<1d-3) then
							!cycle
						!endif
						dvar(i)=dvar(i)+real(S(i,k,1)*conjg(S(j,k,1))*g(j,1)/eg(k))
					enddo
				enddo
			enddo
			write(*,"(es11.4$)")real(phy(1)%E)
			write(*,"(es9.2$)")er
			write(20,"(I4$)")ne
			write(20,"(es13.5$)")var(:)%val,real(g(:,1)),real(phy(1)%E),er
			write(20,"(x)")
			allocate(allE(l)%var(size(var)))
			allE(l)%val=real(phy(1)%E)
			allE(l)%var=var(:)%val
			do i=max(l-m+1,1),l-1
				scv=scv+sign(1d0,real(phy(1)%E)-allE(i)%val)
				if(max(l-m+1,1)>1) then
					scv=scv+sign(1d0,allE(max(l-m+1,1)-1)%val-allE(i)%val)
				endif
			enddo
			write(*,"(es9.2$)")scv/(min(l*max(l-1,1),m*(m-1)))
			if(allE(ml)%val>allE(l)%val) then
				ml=l
			endif
			!if(l>(size(var)+5)*size(var)*2.or.l==size(allE)) then
			!if(l>size(var)*2.or.l==size(allE)) then
			if(allE(max(l-1,1))%val<(real(phy(1)%E)-2d0*er)) then
				write(*,"(A$)")"var err"
			endif
			h=-g(:,1)+sum(-(-g(:,1)-cg)*g(:,1))/sum(cg*cg)*h
			if(sum(abs(allE(max(l-1,1))%var-allE(l)%var))<1d-4) then
				dt=dt/2d0
				if(dt<1d-4) then
					write(*,*)"variational finish",l
					call qsort(allE(:l))
					var(:)%val=0d0
					i=1
					do
						if(abs(allE(i)%val-allE(1)%val)<er/3d0) then
							var(:)%val=var(:)%val+allE(i)%var
							i=i+1
						else
							exit
						endif
					enddo
					var(:)%val=var(:)%val/real(i-1)
					write(20,"(x)")
					exit
				endif
			endif
			!if(l>50.and.abs(scv)/(min(l*max(l-1,1),m*(m-1)))<9d-3.and.scv<0.or.l==size(allE)) then
			if(l>(ml+10)) then
				write(*,*)"variational finish",l
				call qsort(allE(:l))
				var(:)%val=0d0
				i=1
				do
					if(abs(allE(i)%val-allE(1)%val)<er/3d0) then
						var(:)%val=var(:)%val+allE(i)%var
						i=i+1
					else
						exit
					endif
				enddo
				var(:)%val=var(:)%val/real(i-1)
				write(20,"(x)")
				exit
			endif
			call oned_search(allE(l)%val,er,h,dt,Nmc,cfg,nd)
			cg=-g(:,1)
			!if(mod(l,size(var)+5)==0) then
				!ml=l+1
				!cg=1d0
				!h=0d0
				!call random_number(var(:vn)%val)
				!var(:vn)%val=sign(1d0,allE(l)%var)*var(:vn)%val
			!endif
			!var(:)%val=var(:)%val-dvar*dt
			write(*,"(i4)")l
			!Nmc(1)=500
			!Nmc(1)=2000
			l=l+1
		enddo
	end subroutine
	subroutine oned_search(E,er,h,dt,Nmc,cfg,nd)
		real(8) :: h(:),E,dt
		complex(8) :: wf(Ns*2,Ns)
		integer :: k,Nmc(:),l
		integer :: cfg(:,:),nd(:)
		real(8) :: er
		type(randomNumberSequence) :: rnd(n_omp)
		type(t_phy) :: phy(n_omp)
		type(t_sort) :: mE
		allocate(mE%var(size(var)))
		mE%var=var(:)%val
		mE%val=E
		write(*,"(x)")
		write(*,"(A$)")"-->"
		write(*,"(es11.4$)")var(:)%val,E,dt
		write(*,"(x)")
		l=0
		var(:)%val=mE%var-h*dt*3.5d0
		do
			l=l+1
			var(:)%val=var(:)%val+h*dt
			call iniwf(wf)
			!$OMP PARALLEL DO
			do k=1,n_omp
				call mt_init_random_seed(rnd(k),seed=k*10)
				call mc(rnd(k),wf,Nmc,cfg(:,k),nd(k),phy(k))
			enddo
			!$OMP END PARALLEL DO
			!er=real(phy(1)%E)**2
			do k=2,n_omp
				phy(1)%E=phy(1)%E+phy(k)%E
				!er=er+real(phy(k)%E)**2
			enddo
			phy(1)%E=phy(1)%E/n_omp
			!er=sqrt(abs(er/n_omp-real(phy(1)%E)**2))
			if(real(phy(1)%E)<mE%val) then
				mE%val=real(phy(1)%E)
				mE%var=var(:)%val
			endif
			if((real(phy(1)%E)>(mE%val+er*4).and.l>5).or.l>20) then
				exit
			endif
			write(*,"(A$)")">"
			write(*,"(es11.4$)")var(maxloc(abs(h)))%val,real(phy(1)%E)
			write(*,"(x)")
		enddo
		var(:)%val=mE%var
		E=mE%val
	end subroutine
	subroutine vmc_main()
		integer :: n,i,j,Nmc(3),Nvar(3)
		integer :: cfg(Ns,n_omp),nd(n_omp)
		complex(8) :: wf(Ns*2,Ns)
		real(8) :: er
		type(t_phy) :: phy(n_omp)
		type(randomNumberSequence) :: rnd(n_omp)
		do i=1,n_omp
			allocate(phy(i)%csite(Ns),phy(i)%ssite(Ns))
		enddo
		Nmc(3)=1024*4
		Nvar(3)=128*4
		ne=Ns
		ne2=ne/2
		do i=1,Ns
			cfg(i,:)=i
		enddo
		do i=1,n_omp
			call fisher_yates_shuffle(cfg(:,i),Ns)
		enddo
		nd=0

		open(40,file="../data/phyvar_input.dat",status="old",action="read",form="formatted")
		! calculation started
		do i=0,40,2
			ne=Ns-i
			!ne=Ns-24
			ne=ne+mod(ne,2)
			ne2=ne/2
			call iniwf(wf,cp=var(1))
			!var(1)%val=-0.1d0*i-0.7d0
			!var(2)%val=1d-4
			!var(3)%val=1d-1
			!var(1:)%val=(/-1.37E+00,8.01E-01,-4.67E-01/)+(/0d0,0d0,1d0/)*i/300d0
			!read(40,fmt="(i4,5e16.8)")ne,var(:)%val
			!var(2)%val=var(2)%val/2d0
			Nmc(1:2)=(/50000,2**6/)
			Nvar(1:2)=(/10000,ne/)
			!var(:)%val=(/-1.73E+00,1.45E-03,-1.22E-02,1.87E+00,-4.51E-01,-6.08E-01,-1.84E-03/)
			call variational(Nvar,cfg,nd)
			write(*,"(i4$)")ne
			write(*,"(es9.2$)")var%val
			call iniwf(wf)
			!$OMP PARALLEL DO
			do j=1,n_omp
				call mt_init_random_seed(rnd(j),seed=j*50)
				call mc(rnd(j),wf,Nmc,cfg(:,j),nd(j),phy(j))
				!call mc(rnd(j),wf,Nvar,cfg(:,j),nd(j),phy(j))
			enddo
			!$OMP END PARALLEL DO
			er=real(phy(1)%E)**2
			do j=2,n_omp
				phy(1)%E=phy(1)%E+phy(j)%E
				phy(1)%af=phy(1)%af+phy(j)%af
				phy(1)%dsc=phy(1)%dsc+phy(j)%dsc
				phy(1)%ddw=phy(1)%ddw+phy(j)%ddw
				phy(1)%csite=phy(1)%csite+phy(j)%csite
				phy(1)%ssite=phy(1)%ssite+phy(j)%ssite
				er=er+real(phy(j)%E)**2
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
			write(30,"(es13.5$)")real(Ns-ne)/Ns,var%val,real(phy(1)%E),er,phy(1)%af,phy(1)%dsc,phy(1)%ddw
			write(30,"(1x)")
			do j=1,Ns
				write(101,"(2es13.2,es13.2)")i2r(j,:),phy(1)%csite(j)
			enddo
			write(101,"(x/)")
			do j=1,Ns
				if(is_a(j)) then
					write(101,"(2es13.2,2es13.2)")i2r(j,:),phy(1)%ssite(j),1d0
				else
					write(101,"(2es13.2,2es13.2)")i2r(j,:),phy(1)%ssite(j),-1d0
				endif
			enddo
			write(101,"(x/)")
			!V=V+2
			!stop
		enddo
	end subroutine
end module

program main
	use M_vmc
	implicit none
	logical :: f
	! file
	f=openfile(101,"../data/check.dat")
	f=openfile(20,"../data/var.dat")
	!f=openfile(40,"../data/2d.dat")
	!f=openfile(10,"../data/err.dat")
	!f=openfile(30,"../data/phyvar.dat",access="direct",recl=226)
	f=openfile(30,"../data/phyvar.dat")

	call initial()
	call vmc_main()
end program
