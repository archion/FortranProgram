module global
	use M_const
	use M_utility
	use M_lattice
	implicit none
	type t_sort
		real(8) :: val
		real(8), allocatable :: var(:)
	end type
	type t_var
		real(8) :: val
		integer :: nb
		complex(8), allocatable :: bd_sg(:)
		integer :: sg   ! 1:  pair
						! 2:  on-site charge
						! 3:  on-site spin
						! 4:  bond charge
						! -1: n-n jast
	end type
	integer :: ne,ne2,vn
	real(8), parameter :: t(1)=(/1d0/)
	type(t_var), allocatable :: var(:)
	logical :: no_sc
contains
	subroutine initial()
		integer :: i,l,k,Nb
		type(t_var) :: tmp(10)
		call init_random_seed()
		! lattice 
		a1=(/1d0,0d0/)
		a2=(/0d0,1d0/)
		!T1=(/9d0,1d0/)
		!T2=(/-1d0,9d0/)
		T1=a1*12
		T2=a2*12
		!T1=a1*4
		!T2=a2*4
		bdc(1)=1d0
		bdc(2)=1d0
		allocate(sub(1,2))
		sub(1,:)=(/0d0,0d0/)
		call gen_latt()
		call gen_neb()
		call gen_bond(3)
		! finish i2r(i,2),neb(i)%nb(j)%bond(k)/bdc(k)/r(k,2)
		write(*,*)"Total site number is: ",Ns
		do k=1,Ns
			if(is_a(k)) then
				write(101,"(2es13.2,2I5)")i2r(k,:),k,1
			else
				write(101,"(2es13.2,2I5)")i2r(k,:),k,-1
			endif
		enddo
		write(101,"(1x/)")

		Nb=size(bond(1)%bd)
		i=0
		! d-wave sc
		i=i+1
		tmp(i)%sg=1
		tmp(i)%nb=1
		tmp(i)%val=2.47d-01
		allocate(tmp(i)%bd_sg(Nb))
		do k=1,Nb
			tmp(i)%bd_sg(k)=dwave(k)
			write(101,"(2es13.2,I5)")bond(1)%bd(k)%r,nint(dwave(k))
		enddo
		write(101,"(1x/)")

		! cp
		i=i+1
		tmp(i)%sg=2
		tmp(i)%val=-5d-1
		allocate(tmp(i)%bd_sg(Ns))
		tmp(i)%bd_sg(:)=-1d0

		! sdw
		i=i+1
		tmp(i)%sg=3
		tmp(i)%val=1.60d-01
		allocate(tmp(i)%bd_sg(Ns))
		do k=1,Ns
			if(is_a(k)) then
				tmp(i)%bd_sg(k)=1d0
			else
				tmp(i)%bd_sg(k)=-1d0
			endif
		enddo

		! ddw
		i=i+1
		tmp(i)%sg=4
		tmp(i)%nb=1
		tmp(i)%val=1.06d-01
		allocate(tmp(i)%bd_sg(Nb))
		do k=1,Nb
			if(is_a(bond(1)%bd(k)%i(1))) then
				tmp(i)%bd_sg(k)=img*dwave(k)
				write(101,"(4es13.2)")bond(1)%bd(k)%r-0.5d0*bond(1)%bd(k)%dir(:)*dwave(k),bond(1)%bd(k)%dir(:)*dwave(k)
			else
				tmp(i)%bd_sg(k)=-img*dwave(k)
				write(101,"(4es13.2)")bond(1)%bd(k)%r+0.5d0*bond(1)%bd(k)%dir(:)*dwave(k),-bond(1)%bd(k)%dir(:)*dwave(k)
			endif
		enddo
		write(101,"(1x/)")

		!! d-wave cdw
		!i=i+1
		!tmp(i)%sg=4
		!tmp(i)%nb=1
		!tmp(i)%val=0d0
		!allocate(tmp(i)%bd_sg(Nb))
		!do k=1,size(bond(1)%bd)
			!tmp(i)%bd_sg(k)=dwave(k)
		!enddo

		!! hp
		!i=i+1
		!tmp(i)%sg=4
		!tmp(i)%nb=2
		!tmp(i)%val=0d0
		!allocate(tmp(i)%bd_sg(size(bond(2)%bd)))
		!tmp(i)%bd_sg(:)=-1d0
		!i=i+1
		!tmp(i)%sg=4
		!tmp(i)%nb=3
		!tmp(i)%val=0d0
		!allocate(tmp(i)%bd_sg(size(bond(3)%bd)))
		!tmp(i)%bd_sg(:)=-1d0

		vn=i
		! n-n jast
		!i=i+1
		!tmp(i)%sg=-1
		!tmp(i)%val=0d0
		!allocate(tmp(i)%bd_sg(Ns))
		!tmp(i)%sg=-1
		!tmp(i)%val=0d0
		!allocate(tmp(i)%bd_sg(Nb))

		allocate(var(i))
		do l=1,i
			allocate(var(l)%bd_sg(size(tmp(l)%bd_sg)))
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
		complex(8) :: H(Ns*2,Ns*2),D(Ns*2,Ns*2),wf(:,:),dH(Ns*2,Ns*2,size(var)),tmp(Ns,Ns)
		complex(8), optional :: dwf(:,:,:)
		real(8) :: E(Ns*2),r
		real(8), optional :: cp
		complex(8) :: bd
		integer :: i,j,k,l,n,info,sg
		logical :: no_sc
		no_sc=.true.
		H=0d0
		do k=1,size(bond(1)%bd)
			i=bond(1)%bd(k)%i(1)
			j=bond(1)%bd(k)%i(2)
			bd=bond(1)%bd(k)%bdc
			H(i,j)=H(i,j)-1d0*conjg(bd)
			H(j,i)=H(j,i)-1d0*bd
			H(i+Ns,j+Ns)=H(i+Ns,j)+1d0*bd
			H(j+Ns,i+Ns)=H(j+Ns,i)+1d0*conjg(bd)
			!!sc
			!no_sc=.false.
			!H(i,j+Ns)=H(i,j+Ns)+1d-3*dwave(k)*bd
			!H(j,i+Ns)=H(j,i+Ns)+1d-3*dwave(k)*bd
			!H(i+Ns,j)=H(i+Ns,j)+1d-3*dwave(k)*conjg(bd)
			!H(j+Ns,i)=H(j+Ns,i)+1d-3*dwave(k)*conjg(bd)
			!if(nint(imag(bd))/=0) then
				!write(*,*)i,j,bd
				!read(*,*)
			!endif
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
			case(2)
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
		if(no_sc) then
			call heevd(H(:Ns,:Ns),E(:Ns),"V","U",info)
			call heevd(H(Ns+1:,Ns+1:),E(Ns+1:),"V","U",info)
			!write(101,"(es13.2$)")E(:Ns)
			!write(101,"(x)")
			!write(101,"(es13.2$)")E(Ns+1:)
		else
			if(present(cp)) then
				call heevd(H(:Ns,:Ns),E(:Ns),"N","U",info)
				cp=cp+(E(ne2)+E(ne2+1))/2
			endif
			do k=1,size(bond(1)%bd)
				i=bond(1)%bd(k)%i(1)
				j=bond(1)%bd(k)%i(2)
				bd=bond(1)%bd(k)%bdc
				call random_number(r)
				H(i,j+Ns)=H(i,j+Ns)+r*1d-5*bd
				H(j,i+Ns)=H(j,i+Ns)+r*1d-5*bd
				H(i+Ns,j)=H(i+Ns,j)+r*1d-5*conjg(bd)
				H(j+Ns,i)=H(j+Ns,i)+r*1d-5*conjg(bd)
			enddo
			call heevd(H,E,"V","U",info)
			!write(101,"(es13.2$)")E(:Ns)
			!write(101,"(x)")
			!write(101,"(es13.2$)")E(Ns+1:)
		endif

		if(present(dwf)) then
			dH=0d0
			do l=1,size(var)
				D=0d0
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
				case(2)
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
				do i=1,Ns*2
					do j=1,Ns*2
						if(abs(E(i)-E(j))<1d-10) then
							cycle
						endif
						dH(:,i,l)=dH(:,i,l)+dot_product(H(:,j),D(:,i))*H(:,j)/((E(i)-E(j)))
					enddo
				enddo
				!$OMP END PARALLEL DO
			enddo
		endif
		if(no_sc) then
			do i=ne2+1,Ns
				call swap(H(:,i),H(:,i+Ns-ne2))
				do l=1,size(var)
					call swap(dH(:,i,l),dH(:,i+Ns-ne2,l))
				enddo
			enddo
		endif
		wf=conjg(H(:,1:Ns))
		if(present(dwf)) then
			dwf=conjg(dH(:,1:Ns,:))
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
	subroutine getdiff(wf,i,j,cfgi,cfgj,sg,vu,cr,flag,dwf,dvu)
		integer :: i,j,sg,cfgi,cfgj,cr(:)
		complex(8) :: vu(:,:),wf(:,:)
		complex(8), optional :: dvu(:,:,:),dwf(:,:,:)
		logical, optional :: flag
		select case(sg) 
		case(1:4) 
			cr=(/i,-1/)
			if(present(flag).and.flag) then
				dvu(:,1,:)=dwf(cfgj,:,:)
			endif
			vu(:,1)=wf(cfgj,:)
		case(5,7) 
			cr=(/j,-1/)
			if(present(flag).and.flag) then
				dvu(:,1,:)=dwf(cfgi+Ns,:,:)
			endif
			vu(:,1)=wf(cfgi+Ns,:)
		case(6,8) 
			cr=(/j+ne2,-1/)
			if(present(flag).and.flag) then
				dvu(:,1,:)=dwf(cfgi+Ns,:,:)
			endif
			vu(:,1)=wf(cfgi+Ns,:)
		case(0)
			cr=(/i,i+ne2/)
			if(present(flag).and.flag) then
				dvu(:,1,:)=dwf(cfgj,:,:)
				dvu(:,2,:)=dwf(cfgj+Ns,:,:)
			endif
			vu(:,1)=wf(cfgj,:)
			vu(:,2)=wf(cfgj+Ns,:)
		end select
	end subroutine
	subroutine update_swap(i,j,nd,cfg,icfg,A,iA,dA,sg)
		complex(8) :: A(:,:),iA(:,:),dA(:,:,:)
		integer :: i,j,cfg(:),icfg(:),sg,ti,nd
		select case(sg) 
		case(0,2,5)
			call swap(cfg(i),cfg(j))
			call swap(icfg(cfg(i)),icfg(cfg(j)))
		case(1)
			call swap(cfg(i),cfg(j))
			call swap(icfg(cfg(i)),icfg(cfg(j)))
			call swap(A(i+ne2,:),A(j,:))
			call swap(iA(:,i+ne2),iA(:,j))
			do ti=1,vn
				call swap(dA(i+ne2,:,ti),dA(j,:,ti))
			enddo
		case(3)
			call swap(cfg(i),cfg(ne2-nd+1))
			call swap(icfg(cfg(i)),icfg(cfg(ne2-nd+1)))
			call swap(cfg(j),cfg(ne-nd+1))
			call swap(icfg(cfg(j)),icfg(cfg(ne-nd+1)))
			call swap(cfg(ne2-nd+1),cfg(ne-nd+1))
			call swap(icfg(cfg(ne2-nd+1)),icfg(cfg(ne-nd+1)))
			call swap(A(i,:),A(ne2-nd+1,:))
			call swap(A(j,:),A(ne-nd+1,:))
			call swap(iA(:,i),iA(:,ne2-nd+1))
			call swap(iA(:,j),iA(:,ne-nd+1))
			do ti=1,vn
				call swap(dA(i,:,ti),dA(ne2-nd+1,:,ti))
				call swap(dA(j,:,ti),dA(ne-nd+1,:,ti))
			enddo
			nd=nd-1
		case(4)
			call swap(cfg(i),cfg(ne2-nd))
			call swap(icfg(cfg(i)),icfg(cfg(ne2-nd)))
			call swap(cfg(j),cfg(ne-nd))
			call swap(icfg(cfg(j)),icfg(cfg(ne-nd)))
			call swap(cfg(ne2-nd),cfg(ne-nd))
			call swap(icfg(cfg(ne2-nd)),icfg(cfg(ne-nd)))
			call swap(A(i,:),A(ne2-nd,:))
			call swap(A(i+ne2,:),A(ne-nd,:))
			call swap(iA(:,i),iA(:,ne2-nd))
			call swap(iA(:,i+ne2),iA(:,ne-nd))
			do ti=1,vn
				call swap(dA(i,:,ti),dA(ne2-nd,:,ti))
				call swap(dA(i+ne2,:,ti),dA(ne-nd,:,ti))
			enddo
			nd=nd+1
		case(6)
			call swap(cfg(i),cfg(j))
			call swap(icfg(cfg(i)),icfg(cfg(j)))
			call swap(A(i,:),A(j,:))
			call swap(iA(:,i),iA(:,j))
			do ti=1,vn
				call swap(dA(i,:,ti),dA(j,:,ti))
			enddo
		case(7)
			call swap(cfg(i),cfg(ne2-nd+1))
			call swap(icfg(cfg(i)),icfg(cfg(ne2-nd+1)))
			call swap(cfg(j),cfg(ne-nd+1))
			call swap(icfg(cfg(j)),icfg(cfg(ne-nd+1)))
			call swap(A(i,:),A(ne2-nd+1,:))
			call swap(A(j,:),A(ne-nd+1,:))
			call swap(iA(:,i),iA(:,ne2-nd+1))
			call swap(iA(:,j),iA(:,ne-nd+1))
			do ti=1,vn
				call swap(dA(i,:,ti),dA(ne2-nd+1,:,ti))
				call swap(dA(j,:,ti),dA(ne-nd+1,:,ti))
			enddo
			nd=nd-1
		case(8)
			call swap(cfg(j),cfg(ne2-nd))
			call swap(icfg(cfg(j)),icfg(cfg(ne2-nd)))
			call swap(cfg(i),cfg(ne-nd))
			call swap(icfg(cfg(i)),icfg(cfg(ne-nd)))
			call swap(A(j,:),A(ne2-nd,:))
			call swap(A(j+ne2,:),A(ne-nd,:))
			call swap(iA(:,j),iA(:,ne2-nd))
			call swap(iA(:,j+ne2),iA(:,ne-nd))
			do ti=1,vn
				call swap(dA(j,:,ti),dA(ne2-nd,:,ti))
				call swap(dA(j+ne2,:,ti),dA(ne-nd,:,ti))
			enddo
			nd=nd+1
		end select
	end subroutine
	subroutine det(vu,cr,A,iA,Y,pb)
		complex(8) :: vu(:,:),A(:,:),iA(:,:),Y(:,:),pb
		integer :: cr(:)
		if(cr(2)>0) then
			call det_ratio_tworow(vu,cr,A,iA,Y,pb)
		else
			call det_ratio_row(vu(:,1),cr(1),A,iA,pb)
		endif
	end subroutine
	subroutine update_matrix(vu,cr,Y,pb,A,iA,flag,dvu,dA)
		complex(8) :: vu(:,:),A(:,:),iA(:,:),Y(:,:),pb
		complex(8) :: dvu(:,:,:),dA(:,:,:)
		integer :: cr(:),sg
		logical :: flag
		if(cr(2)>0) then
			call inv_update_tworow(vu,cr,Y,A,iA)
			if(flag) then
				dA(cr(1),:,:)=dvu(:,1,:)
				dA(cr(2),:,:)=dvu(:,2,:)
			endif
		else
			call inv_update_row(vu(:,1),cr(1),pb,A,iA)
			if(flag) then
				dA(cr(1),:,:)=dvu(:,1,:)
			endif
		endif
	end subroutine
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
	!real(8) :: DJ=1d0/3d0,V=2d0
	real(8) :: DJ=1d0/3d0,V=0d0
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
	subroutine energy(cfg,icfg,nd,nn,wf,ja,A,iA,El)
		complex(8) :: A(:,:),iA(:,:),vu(size(A,1),2),Y(2,2),pb,wf(:,:)
		complex(8) :: El
		real(8) :: ja(:)
		integer :: cfg(:),icfg(:),i,j,k,l,sg,cr(2),nd,nn(:),nndif(size(nn))
		El=0d0
		do i=1,ne
			do l=1,size(t)
				do k=1,size(neb(cfg(i))%nb(l+1)%neb)
					j=icfg(neb(cfg(i))%nb(l+1)%neb(k))
					call get_case(i,j,0,nd,sg)
					select case(sg) 
					case(1) 
						call nn_dif(cfg(i),cfg(j),cfg,icfg,sg,nd,nndif)
						call getdiff(wf,i,j,cfg(i),cfg(j),sg,vu,cr)
						call det(vu,cr,A,iA,Y,pb)
						El=El-(t(l)*conjg(pb))&
							*jast(nndif,ja)&
							*neb(cfg(i))%nb(l+1)%bdc(k)
					case(5) 
						call nn_dif(cfg(i),cfg(j),cfg,icfg,sg,nd,nndif)
						call getdiff(wf,i,j,cfg(i),cfg(j),sg,vu,cr)
						call det(vu,cr,A,iA,Y,pb)
						El=El+(t(l)*conjg(pb))&
							*jast(nndif,ja)&
							*neb(cfg(i))%nb(l+1)%bdc(k)
					case(0)
						if(l==1.and.i<=ne2) then
							call getdiff(wf,i,j,cfg(i),cfg(j),sg,vu,cr)
							call det(vu,cr,A,iA,Y,pb)
							El=El+0.5d0*DJ*conjg(pb)-0.25d0*DJ
						endif
					case(-1)
						if(l==1.and.k<3) then
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
	!subroutine charge_order(cfg,icfg,nd,wf,A,iA,lphy)
		!complex(8) :: A(:,:),iA(:,:),vu(Ns,2),Y(2,2),pb,wf(:,:),a11
		!complex(8) :: lphy
		!integer :: cfg(:),icfg(:),i,n,m,j,k,l,inb,sg,cr(2),x1(2),x2(2),nr(2),dr(2),r(2),nd
		!lphy=0d0
		!lphy=lphy/Ns
	!end subroutine
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
				!case(1)
					!call getdiff(wf,i,j,cfg(i),cfg(j),sg,vu,cr)
					!call det(vu,cr,A,iA,Y,pb)
					!lphy=lphy+imag(conjg(pb))*&
						!(0.5d0-mod(inb,2))*(0.5d0-mod(sum(x1),2))
				!case(5)
					!call getdiff(wf,i,j,cfg(i),cfg(j),sg,vu,cr)
					!call det(vu,cr,A,iA,Y,pb)
					!lphy=lphy-imag(conjg(pb))*&
						!(0.5d0-mod(inb,2))*(0.5d0-mod(sum(x1),2))
				!end select
			!enddo
		!enddo
		!lphy=lphy/Ns
	!end subroutine
	!subroutine dsc_order(cfg,icfg,nd,wf,A,iA,ja,lphy)
		!complex(8) :: A(:,:),iA(:,:),vu(Ns,2),Y(2,2),pb,wf(:,:),a11
		!real(8) :: lphy,ja(:)
		!integer :: cfg(:),icfg(:),i,n,m,j,k,l,inb,sg,cr(2),x1(2),x2(2),nr(2),dr(2),r(2),nd
		!do i=ne2-nd+1,ne-nd
			!do inb=1,4,1
				!j=icfg(neb(cfg(i),inb,1))
				!if(j>ne2) then
					!vu(:,1)=wf(cfg(j),:)
					!cr=(/i,-1/)
					!call det(vu,cr,A,iA,Y,pb)
					!lphy=lphy+real(conjg(pb)*jast(cfg(i),cfg(j),cfg,icfg,nd,ja))*&
						!(0.5d0-mod(inb,2))
				!endif
			!enddo
		!enddo
	!end subroutine
	!subroutine spin_corelation(cfg,icfg,nd,wf,A,iA,lphy)
		!complex(8) :: A(:,:),iA(:,:),vu(Ns,2),Y(2,2),pb,wf(:,:),a11
		!real(8) :: lphy
		!integer :: cfg(:),icfg(:),i,n,m,j,k,l,inb,sg,cr(2),x1(2),x2(2),nr(2),dr(2),r(2),nd
	!end subroutine
	!subroutine dsc_corelation(cfg,icfg,nd,wf,A,iA,ja,lphy)
		!complex(8) :: A(:,:),iA(:,:),vu(Ns,2),Y(2,2),pb,wf(:,:),a11
		!real(8) :: lphy,ja(:)
		!integer :: cfg(:),icfg(:),i,j,ii,jj,inb1,inb2,cr(2),x1(2),x2(2),nr(2),dr(2),r(2),nd
		!lphy=0d0
		!do i=ne2-nd+1,ne-nd
			!do inb1=1,4,1
				!j=icfg(neb(cfg(i),inb1,1))
				!if(j<=ne2) then
					!do jj=ne-nd+1,Ns
						!do inb2=1,4,1
							!ii=icfg(neb(cfg(jj),inb2,1))
							!if(ii>(ne-nd)) then
								!cr=(/ii,j/)
							!else
								!cycle
							!endif
							!vu(:,1)=wf(cfg(i)+Ns,:)
							!vu(:,2)=wf(cfg(jj),:)
							!call det(vu,cr,A,iA,Y,pb)
							!lphy=lphy+real(dconjg(pb))*&
								!(0.5d0-mod(inb1+inb2,2))
						!enddo
					!enddo
				!endif
			!enddo
		!enddo
		!lphy=lphy/(Ns*Ns)
	!end subroutine
end module

module M_hubbard
	use M_omp_rand
	use M_mc_matrix
	implicit none
	real(8), parameter :: U=10d0,V=0d0
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
	subroutine energy(cfg,icfg,nd,nn,wf,ja,A,iA,El)
		complex(8) :: A(:,:),iA(:,:),vu(size(A,1),2),Y(2,2),pb,wf(:,:)
		complex(8) :: El
		real(8) :: ja(:)
		integer :: cfg(:),icfg(:),i,j,k,l,sg,cr(2),nd,n,ud,nn(:),nndif(size(nn))
		El=0d0
		do n=1,ne
			if(n>ne2) then
				ud=-1
				i=n-nd
			else
				ud=1
				i=n
			endif
			do l=1,size(t)
				do k=1,size(neb(cfg(i))%nb(l+1)%neb)
					j=icfg(neb(cfg(i))%nb(l+1)%neb(k))
					call get_case(i,j,ud,nd,sg)
					select case(sg)
					case(1:4)
						call nn_dif(cfg(i),cfg(j),cfg,icfg,sg,nd,nndif)
						call getdiff(wf,i,j,cfg(i),cfg(j),sg,vu,cr)
						call det(vu,cr,A,iA,Y,pb)
						El=El-(t(l)*conjg(pb)*jast(nndif,ja))
					case(5:8)
						call nn_dif(cfg(i),cfg(j),cfg,icfg,sg,nd,nndif)
						call getdiff(wf,i,j,cfg(i),cfg(j),sg,vu,cr)
						call det(vu,cr,A,iA,Y,pb)
						El=El+(t(l)*conjg(pb)*jast(nndif,ja))
					end select
				enddo
			enddo
		enddo
		El=El+U*nd+V*nn(2)
		El=El/Ns
		!stop
	end subroutine
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
							!lphy=lphy+real(dconjg(pb)*jast(cfg(i),cfg(ii),cfg,icfg,1,nd,ja)*jast(cfg(j),cfg(jj),cfg,icfg,1,nd,ja))*&
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
	subroutine mc(rnd,wf,Nmc,cfg,nd,phyval,dwf,S,g)
		complex(8) :: pb,A(Ns,Ns),iA(Ns,Ns),dA(Ns,Ns,vn),vu(Ns,2),dvu(Ns,2,vn),wf(:,:),Y(2,2)
		complex(8), optional :: S(:,:),g(:),dwf(:,:,:)
		real(8) :: rpb
		real(8) :: ja(size(var)-vn)
		complex(8) :: lO(size(var)),lS(size(lO),size(lO)),lg(size(lO)),O(size(lO))
		complex(8) :: El,E
		real(8) :: phyval(:),lphy(size(phyval))
		integer :: icfg(Ns),i,j,ti,tj,l,sg,n,Nmc(:),cr(2),acp,bi,bj,ii,jj
		integer :: cfg(Ns),nd,nn(max(size(var)-vn,2)),nndif(size(nn))
		logical :: flag
		type(randomNumberSequence) :: rnd

		if(mod(Ns,2)/=0) then
			write(*,*)"Number of site is not even, exit!!"
			stop
		endif
		if(ne/=ne2*2) then
			write(*,*)"Number of electrons is not even, exit!!"
			stop
		endif

		ja=var(vn+1:)%val

		do i=1,Ns
			icfg(cfg(i))=i
		enddo

		flag=present(S)

		!do i=1,Ns
			!cfg(i)=i
		!enddo
		!nd=0
		!call fisher_yates_shuffle(cfg,Ns)

		n=0
		acp=0
		call ini_nn(cfg,icfg,nd,nn)
		!!$OMP CRITICAL
			!write(*,*)nn
		!!$OMP END CRITICAL
		!write(*,"(A)")"start monte carlo"

		E=0
		phyval=0d0
		lphy=0d0
		if(flag) then
			O=0d0
			S=0d0
			g=0d0
		endif

		!get initial matrix
		do i=1,Ns
			if(i<=ne2) then
				A(i,:)=wf(cfg(i),:)
				if(flag) then
					dA(i,:,:)=dwf(cfg(i),:,:)
				endif
			elseif(i>=ne-nd+1)  then
				A(i,:)=wf(cfg(i)+Ns,:)
				if(flag) then
					dA(i,:,:)=dwf(cfg(i)+Ns,:,:)
				endif
			else
				A(i,:)=wf(cfg(i-ne2)+Ns,:)
				if(flag) then
					dA(i,:,:)=dwf(cfg(i-ne2)+Ns,:,:)
				endif
			endif
		enddo
		iA=A
		call mat_inv(iA)
		!write(*,*)sum(abs(matmul(A,iA)-diag(1d0,size(A,1))))
		do
			n=n+1
			call two(rnd,i,j,nd,sg)
			!write(*,*)i,j,nd,sg
			call getdiff(wf,i,j,cfg(i),cfg(j),sg,vu,cr,flag,dwf,dvu)
			call det(vu,cr,A,iA,Y,pb)
			call nn_dif(cfg(i),cfg(j),cfg,icfg,sg,nd,nndif)
			!call random_number(rpb)
			rpb=getRandomReal(rnd)
			if(rpb<real(pb*conjg(pb))*(jast(nndif,ja))**2) then
				acp=acp+1
				call update_matrix(vu,cr,Y,pb,A,iA,flag,dvu,dA)
				!write(*,*)sum(abs(matmul(A,iA)-diag(1d0,size(A,1))))
				call update_swap(i,j,nd,cfg,icfg,A,iA,dA,sg)
				nn=nn+nndif
				!if(sum(abs(matmul(A,iA)-diag(1d0,size(A,1))))>1d-7) then
					!write(*,*)i,j,sg,nd,sum(abs(matmul(A,iA)-diag(1d0,size(A,1))))
				!endif
				if(mod(n,500)==0) then
					iA=A
					call mat_inv(iA)
				endif
			endif
			!if(n==Nmc(1)) then
				!write(*,*)"hot finish"
			!endif
			if(n>Nmc(1).and.mod(n-Nmc(1),Nmc(2))==0) then
				call energy(cfg,icfg,nd,nn,wf,ja,A,iA,El)
				if(flag) then
					do ti=1,vn
						lO(ti)=Tr(iA,dA(:,:,ti))
					enddo
					if(size(ja)/=0) then
						lO(vn+1:)=-nn
					endif
					do ti=1,size(lO)
						do tj=1,size(lO)
							lS(ti,tj)=dconjg(lO(ti))*lO(tj)
						enddo
					enddo
					lg=real(El*lO)
					S=S+lS
					g=g+lg
					O=O+lO
				!else
					!call spin_order(cfg,nd,lphy(2))
					!call ddw_order(cfg,icfg,nd,wf,A,iA,ja,lphy(3))
					!call dsc_corelation(cfg,icfg,nd,wf,A,iA,ja,lphy(4))
				endif
				E=E+El
				phyval=phyval+lphy
			endif
			if(n>=(Nmc(1)+Nmc(2)*Nmc(3))) then
				E=E/Nmc(3)
				phyval=phyval/Nmc(3)
				phyval(1)=real(E)
				if(flag) then
					S=S/Nmc(3)
					g=g/Nmc(3)
					O=O/Nmc(3)
					do ti=1,size(O)
						do tj=1,size(O)
							S(ti,tj)=S(ti,tj)-dconjg(O(ti))*O(tj)
						enddo
					enddo
					g=2d0*(g-real(E*O))*Ns
				endif
				exit
			endif
		enddo
		!write(*,"('accept/total number is ',I6,'/',I6)")acp,n
	end subroutine
	subroutine variational(Nmc,cfg,nd)
		use lapack95, only: heevd,heevr
		type(t_sort) :: allE(300)
		real(8) :: dvar(size(var)),pvar(size(var)),eg(size(var)),dt=0.04d0,scv,var_av(size(var))
		complex(8) :: wf(Ns*2,Ns),dwf(Ns*2,Ns,vn),g(size(var)),S(size(g),size(g)),S_omp(size(g),size(g)),g_omp(size(g))
		integer :: n,m,i,j,k,info,Nmc(:),nm,l,nav
		integer :: cfg(Ns),nd,cfg_omp(Ns),nd_omp
		real(8) :: phyval(5),phyval_omp(5),er(5)
		logical :: flag
		type(randomNumberSequence) :: rnd
		m=50
		n=24
		l=1
		scv=0d0
		flag=.true.
		var_av=0d0
		nav=0
		call iniwf(wf,dwf,var(2)%val)
		do
			write(*,"(I3$)")ne
			write(*,"(es9.2$)")var(:)%val
			call iniwf(wf,dwf)
			S=0d0
			g=0d0
			er=0d0
			phyval=0d0
			!$OMP PARALLEL DO REDUCTION(+:phyval,S,g,er) PRIVATE(phyval_omp,S_omp,g_omp,rnd) LASTPRIVATE(cfg_omp,nd_omp)
			do k=1,n
				call mt_init_random_seed(rnd,seed=k*10)
				!call mt_init_random_seed(rnd,seed=k*12301)
				cfg_omp=cfg
				nd_omp=nd
				call mc(rnd,wf,Nmc,cfg_omp,nd_omp,phyval_omp,dwf,S_omp,g_omp)
				phyval(1)=phyval(1)+phyval_omp(1)
				phyval(2:)=phyval(2:)+abs(phyval_omp(2:))
				S=S+S_omp
				g=g+g_omp
				er=er+phyval_omp**2
				call finalize_RandomNumberSequence(rnd)
			enddo
			!$OMP END PARALLEL DO
			cfg=cfg_omp
			nd=nd_omp
			phyval=phyval/n
			S=S/n
			g=g/n
			!if(n<5) then
				!er=sqrt((sga(nm-3,1)-sga(size(sga,1)/2,1))/(2**(nm-(nm-3)+1)-1))
			!else
				er=sqrt(abs(er/n-phyval**2))
			!endif
			S=S+diag(1d-2,size(S,1))
			call heevd(S,eg,"V")
			!if(any(eg<0d0)) then
				!eg=eg+minval(abs(eg))*1.1d0
			!endif
			!eg=eg+1d-1
			!write(*,*)eg
			do i=1,size(g)
				dvar(i)=0d0
				do j=1,size(g)
					do k=1,size(g)
						!if(eg(k)/eg(size(g))<1d-3) then
							!cycle
						!endif
						dvar(i)=dvar(i)+real(S(i,k)*dconjg(S(j,k))*g(j)/eg(k))
					enddo
				enddo
			enddo
			!write(*,"(es9.2$)")real(g)
			write(*,"(es11.4$)")phyval(1)
			write(*,"(es9.2$)")er(1)
			write(20,"(I4$)")ne
			write(20,"(es13.5$)")var(:)%val,real(g),phyval(1),er(1)
			write(20,"(x)")
			allocate(allE(l)%var(size(var)))
			allE(l)%val=phyval(1)
			allE(l)%var=var(:)%val
			if(allE(max(l-1,1))%val<(phyval(1)-2d0*er(1))) then
				write(*,"(A)")"var err"
				!exit
			endif
			do i=max(l-m+1,1),l-1
				scv=scv+sign(1d0,phyval(1)-allE(i)%val)
				if(max(l-m+1,1)>1) then
					scv=scv+sign(1d0,allE(max(l-m+1,1)-1)%val-allE(i)%val)
				endif
			enddo
			write(*,"(es9.2$)")scv/(min(l*max(l-1,1),m*(m-1)))
			if(l>m.and.abs(scv)/(min(l*max(l-1,1),m*(m-1)))<5d-3.and.scv<0.or.l==size(allE)) then
				write(*,*)"variational finish"
				call qsort(allE(:l))
				!write(*,*)allE(:l)%val
				var(:)%val=0d0
				do i=1,10
					var(:)%val=var(:)%val+allE(i)%var
				enddo
				!write(10,"(I4$)")ne
				!write(10,"(es13.5$)")var,Ev,er
				!write(10,"(x)")
				var(:)%val=var(:)%val/10d0
				write(20,"(x)")
				exit
			endif
			var(:)%val=var(:)%val-dvar*dt
			!if(abs(abs(var(2)%val)-abs(var(3)%val))<1d-5) then
				!write(*,*)"cp and sdw may be the same"
				!var(2)%val=var(2)%val-dvar(2)*dt
			!endif
			write(*,"(x)")
			Nmc(1)=500
			l=l+1
		enddo
	end subroutine
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
end module

program main
	use M_vmc
	implicit none
	integer :: n,i,j,Nmc(3),Nvar(3),nm,nd,nd_omp
	integer, allocatable :: cfg(:),cfg_omp(:)
	complex(8), allocatable :: wf(:,:)
	real(8) :: er(5),phyval(5)
	real(8), allocatable :: phyval_omp(:)
	logical :: f
	type(randomNumberSequence) :: rnd
	! file
	f=openfile(101,"../data/check.dat")
	f=openfile(20,"../data/var.dat")
	!f=openfile(40,"../data/2d.dat")
	!f=openfile(10,"../data/err.dat")
	!f=openfile(30,"../data/phyvar.dat",access="direct",recl=226)
	f=openfile(30,"../data/phyvar.dat")

	call initial()
	n=32
	Nmc(3)=1024*4
	Nvar(3)=128*4
	allocate(cfg(Ns),cfg_omp(Ns),wf(Ns*2,Ns),phyval_omp(5))
	ne=Ns
	ne2=ne/2
	do i=1,Ns
		cfg(i)=i
	enddo
	nd=0
	call fisher_yates_shuffle(cfg,Ns)


	! calculation started
	do i=1,28,2
	!do i=5,16,1
		var(1)%val=1d-1
		!var(2)%val=0d0
		var(3)%val=1d-1
		var(4)%val=1d-1
		!var(:)%val=(/1.22835d-01/)
		!var(1)%val=var(1)%val+2d-2
		!ne=Ns*7/8
		!ne=Ns-8
		ne=Ns-i*2
		ne=ne+mod(ne,2)
		!read(30,rec=i+1,fmt="(i4,7es13.5)")ne,var
		ne2=ne/2
		nd=0
		!call self_n(var)
		Nmc(1:2)=(/50000,2**6/)
		Nvar(1:2)=(/50000,2**6/)
		!Nvar(1:2)=(/50000,ne/)
		call variational(Nvar,cfg,nd)
		write(*,"(i4$)")ne
		write(*,"(es9.2$)")var%val
		call iniwf(wf)
		er=0d0
		phyval=0d0
		!$OMP PARALLEL DO REDUCTION(+:phyval,er) PRIVATE(phyval_omp,rnd) LASTPRIVATE(cfg_omp,nd_omp)
		do j=1,n
			call mt_init_random_seed(rnd,seed=j*50)
			cfg_omp=cfg
			nd_omp=nd
			!call mc(rnd,wf,Nvar,cfg_omp,nd_omp,phyval_omp)
			call mc(rnd,wf,Nmc,cfg_omp,nd_omp,phyval_omp)
			phyval(1)=phyval(1)+phyval_omp(1)
			phyval(2:)=phyval(2:)+abs(phyval_omp(2:))
			er=er+phyval_omp**2
		enddo
		!$OMP END PARALLEL DO
		cfg=cfg_omp
		nd=nd_omp
		phyval=phyval/n
		er=sqrt(abs(er/n-phyval**2))
		!do j=1,size(sga,1)/2-1
			!write(10,"(es13.5$)")(sqrt((sga(j,k)-sga(size(sga,1)/2,k))/(2**(nm-j+1)-1)),sqrt((sga(j,k)-sga(size(sga,1)/2,k))/(2**(nm-j+1)-1)*sqrt(2d0/(2**(nm-j+1)-1))),k=1,size(sga,2))
			!write(10,"(x)")
		!enddo
		!write(*,"(es9.2$)")phyval(1),sqrt((sga(nm-5,1)-sga(size(sga,1)/2,1))/(2**(nm-(nm-5)+1)-1))
		write(*,"(es11.4$)")phyval(1)
		write(*,"(es9.2$)")phyval(2:)
		write(*,"(1x)")
		!write(30,rec=i+1,fmt="(i4,17es13.5,A)")ne,var,(phyval(j),sqrt((sga(nm-5,j)-sga(size(sga,1)/2,j))/(2**(nm-(nm-5)+1)-1)),j=1,5),char(10)
		!write(30,"(es13.5$)")real(Ns-ne)/Ns,var,(phyval(j),sqrt((sga(nm-5,j)-sga(size(sga,1)/2,j))/(2**(nm-(nm-5)+1)-1)),j=1,5)
		write(30,"(es13.5$)")real(Ns-ne)/Ns,var%val,(phyval(j),er(j),j=1,size(phyval))
		write(30,"(1x)")
		!V=V+2
	enddo
end program
