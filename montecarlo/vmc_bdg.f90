module M_pmt
	use M_const
	use M_utility
	use M_latt, only : &
		 latt         =>  square_tilda         ,& 
		 latt_bc      =>  square_tilda_bc      ,& 
		 latt_one2two =>  square_tilda_one2two ,& 
		 latt_two2one =>  square_tilda_two2one 
		 !latt         =>  square         ,& 
		 !latt_bc      =>  square_bc      ,& 
		 !latt_one2two =>  square_one2two ,& 
		 !latt_two2one =>  square_two2one 
	implicit none
	integer, parameter :: Ns(2)=(/9,10/),Ns2=Ns(1)**2+1,Tx(2)=(/Ns(1),-1/),Ty(2)=(/1,Ns(1)/),vn=6
	!integer, parameter :: Ns(2)=(/10,10/),Ns2=Ns(1)**2,Tx(2)=(/Ns(1),0/),Ty(2)=(/0,Ns(1)/),vn=6
	integer :: neb(Ns2,4,3),ne=Ns2,ne2=Ns2/2
	real(8), parameter :: t(1)=(/1d0/)
end module

module M_wf
	use M_pmt
	use M_matrix
	use lapack95, only: heevd,heevr
	implicit none
contains
	subroutine self_n(var)
		complex(8) :: H(Ns2*2,Ns2*2)
		real(8) :: n,var(:),eg(Ns2*2),step,pn
		integer :: sg,i
		pn=0d0
		step=0.01d0
		do
			call bdg((/0d0,var(2),0d0,0d0,0d0,0d0/),H)
			call heevd(H,eg,"V")
			n=0d0
			do i=1,Ns2
				n=n+dot_product(H(i,1:Ns2),H(i,1:Ns2))
				n=n+1d0-dot_product(H(i+Ns2,1:Ns2),H(i+Ns2,1:Ns2))
			enddo
			if(abs(n-ne)<0.1d0) then
				exit
			endif
			call find_cross(pn,n-ne,sg)
			if(sg/=0) then
				step=step*0.5d0
			endif
			var(2)=var(2)+step*sign(1d0,ne-n)
		enddo
	end subroutine
	subroutine wf_r(var,wf,dwf)
		complex(8) :: H(Ns2*2,Ns2*2),dH(Ns2*2,Ns2*2,vn),wf(:,:),dtmp(size(wf,1),size(wf,2)),D(Ns2*2,Ns2*2,vn)
		complex(8), optional :: dwf(:,:,:)
		real(8) :: eg(Ns2*2),var(vn),nvar(vn),r
		integer :: i,j,n,ix(2),bc
		call bdg(var,H,D)
		call heevd(H,eg,"V")
		if(present(dwf)) then
			call getDer(D,H,eg,dH)
			dwf=conjg(dH(:,1:Ns2,:))
		endif
		wf=conjg(H(:,1:Ns2))
		if(any(isnan(real(wf)))) then
			write(*,*)"wf is NAN"
			stop
		endif
	end subroutine
	subroutine getDer(D,H,eg,dH)
		complex(8) :: H(:,:),D(:,:,:),dH(:,:,:),X(size(H,2)),C(size(H,1),size(H,2)),tmp(size(D,1)),&
			va(size(D,1)*10)
		integer :: ja(size(va)),ia(size(D,1)+1)
		real(8) :: eg(:)
		integer :: i,j,n,m(size(H,2))
		do i=1,size(H,2)
			X(i)=H(1,i)
			m(i)=1
			do j=2,size(H,1)
				if(abs(X(i))<abs(H(j,i))) then
					X(i)=H(j,i)
					m(i)=j
				endif
			enddo
			H(:,i)=H(:,i)/X(i)
		enddo
		dH=0d0
		do n=1,size(D,3)
			C=0d0
			call crs(D(:,:,n),va,ja,ia)
			!$OMP PARALLEL DO PRIVATE(tmp)
			do i=1,size(H,2)
				!tmp=matmul(D(:,:,n),H(:,i))
				call crsmv(va,ja,ia,H(:,i),tmp)
				do j=1,size(H,1)
					if(abs(eg(i)-eg(j))<1d-10) then
						cycle
					endif
					C(j,i)=dot_product(H(:,j),tmp)*X(j)*conjg(X(j))/((eg(i)-eg(j)))
				enddo
				C(i,i)=-sum(H(m(i),:)*C(:,i))
				dH(:,i,n)=matmul(H(:,:),C(:,i))
			enddo
			!$OMP END PARALLEL DO
			!dH(:,:,n)=matmul(H,C)
		enddo
	end subroutine
	subroutine bdg(var,H,D)
		complex(8) :: H(Ns2*2,Ns2*2)
		complex(8), optional :: D(:,:,:)
		real(8) :: eg(Ns2*2),var(vn)
		integer :: i,j,n,ix(2),bc
		H=0d0
		if(present(D)) then
			D=0d0
		endif
		do i=1,Ns2
			call latt_one2two(i,Ns,ix)
			do j=1,4
				!call random_number(r)
				!r=(r-0.5d0)*1d-8
				bc=1
				!if(abs(neb(i,j,1)-i)>Ns(1)+3) then
					!bc=-1
				!endif
				H(i,neb(i,j,1))=(-t(1)+var(4)*img*(-1)**mod(sum(ix),2)*(-1)**mod(j,2))*bc
				H(i+Ns2,neb(i,j,1)+Ns2)=(t(1)+var(4)*img*(-1)**mod(sum(ix),2)*(-1)**mod(j,2))*bc
				H(i,neb(i,j,2))=-var(5)*bc
				H(i+Ns2,neb(i,j,2)+Ns2)=-H(i,neb(i,j,2))
				H(i,neb(i,j,3))=-var(6)*bc
				H(i+Ns2,neb(i,j,3)+Ns2)=-H(i,neb(i,j,3))
				H(i,neb(i,j,1)+Ns2)=var(1)*(-1)**mod(j,2)*bc
				H(i+Ns2,neb(i,j,1))=conjg(H(i,neb(i,j,1)+Ns2))
				if(present(D)) then
					D(i,neb(i,j,1),4)=(img*(-1)**mod(sum(ix),2)*(-1)**mod(j,2))*bc
					D(i+Ns2,neb(i,j,1)+Ns2,4)=(img*(-1)**mod(sum(ix),2)*(-1)**mod(j,2))*bc
					D(i,neb(i,j,2),5)=-bc
					D(i+Ns2,neb(i,j,2)+Ns2,5)=-D(i,neb(i,j,2),5)
					D(i,neb(i,j,3),6)=-bc
					D(i+Ns2,neb(i,j,3)+Ns2,6)=-D(i,neb(i,j,3),6)
					D(i,neb(i,j,1)+Ns2,1)=(-1)**mod(j,2)*bc
					D(i+Ns2,neb(i,j,1),1)=conjg(D(i,neb(i,j,1)+Ns2,1))
				endif
			enddo
			! on site
			H(i,i)=var(3)*(-1)**mod(sum(ix),2)-var(2)
			H(i+Ns2,i+Ns2)=var(3)*(-1)**mod(sum(ix),2)+var(2)
			if(present(D)) then
				D(i,i,3)=(-1)**mod(sum(ix),2)
				D(i+Ns2,i+Ns2,3)=(-1)**mod(sum(ix),2)
				D(i,i,2)=-1d0
				D(i+Ns2,i+Ns2,2)=1d0
			endif
		enddo
	end subroutine
end module

module M_mc_matrix
	use M_matrix
	use M_pmt
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
				dvu(:,1,:)=dwf(cfgi+Ns2,:,:)
			endif
			vu(:,1)=wf(cfgi+Ns2,:)
		case(6,8) 
			cr=(/j+ne2,-1/)
			if(present(flag).and.flag) then
				dvu(:,1,:)=dwf(cfgi+Ns2,:,:)
			endif
			vu(:,1)=wf(cfgi+Ns2,:)
		case(0)
			cr=(/i,i+ne2/)
			if(present(flag).and.flag) then
				dvu(:,1,:)=dwf(cfgj,:,:)
				dvu(:,2,:)=dwf(cfgj+Ns2,:,:)
			endif
			vu(:,1)=wf(cfgj,:)
			vu(:,2)=wf(cfgj+Ns2,:)
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
	function jast(ri,rj,cfg,icfg,nd,ja)
		integer :: ri,rj,icfg(:),cfg(:),sg,i,nd
		real(8) :: jast,ja(:)
		jast=1d0
		if(cs(icfg(ri),nd)==2) then
			jast=jast*exp(ja(1))
		endif
		if(cs(icfg(rj),nd)/=4) then
			jast=jast*exp(-ja(1))
		endif
		do i=2,size(ja)
		enddo
	end function
	subroutine checkcfg(icfg,nd,s)
		integer :: i,s,icfg(:),r0(2),nd
		call latt_one2two(1,Ns,r0)
		write(*,"(1X)")
		do i=1,r0(1)
			write(*,"('  '$)")
		enddo
		do i=1,Ns2
			if(icfg(i)<=ne2-nd) then
				write(*,"(' ●'$)")
			elseif(icfg(i)<=ne2) then
				write(*,"(' ◐'$)")
			elseif(icfg(i)<=ne-nd) then
				write(*,"(' ○'$)")
			else
				write(*,"('  '$)")
			endif
			if(mod(i+r0(1),Ns(1))==0) then
				write(*,"(1X)")
			endif
		enddo
		call sleepqq(s)
		write(*,"(1X)")
	end subroutine
end module

module M_tJ
	use M_rd
	use M_mc_matrix
	implicit none
	real(8), parameter :: DJ=1d0/3d0,V=0d0!DJ*4d0
	!real(8), parameter :: DJ=0.3d0,V=DJ*4d0
contains
	subroutine two(i,j,nd,sg)
		integer :: i,j,sg,n,n1,n0,nd
		n1=ne2-nd
		n0=Ns2-ne+nd
		call random_number(n,n1*n0*2+n1*n1)
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
	subroutine energy(cfg,icfg,nd,wf,ja,A,iA,El)
		complex(8) :: A(:,:),iA(:,:),vu(size(A,1),2),Y(2,2),pb,wf(:,:)
		real(8) :: El
		real(8) :: ja(:)
		integer :: cfg(:),icfg(:),i,ii,j,l,inb,sg,cr(2),nd
		El=0d0
		do i=1,ne
			do ii=1,4
				do inb=1,size(t)
					j=icfg(neb(cfg(i),ii,inb))
					call get_case(i,j,0,nd,sg)
					select case(sg) 
					case(1) 
						call getdiff(wf,i,j,cfg(i),cfg(j),sg,vu,cr)
						call det(vu,cr,A,iA,Y,pb)
						El=El-real(t(inb)*conjg(pb))
					case(5) 
						call getdiff(wf,i,j,cfg(i),cfg(j),sg,vu,cr)
						call det(vu,cr,A,iA,Y,pb)
						El=El+real(t(inb)*conjg(pb))
					case(0)
						if(inb==1.and.i<=ne2) then
							call getdiff(wf,i,j,cfg(i),cfg(j),sg,vu,cr)
							call det(vu,cr,A,iA,Y,pb)
							El=El+real(0.5d0*DJ*conjg(pb)-0.25d0*DJ&
								-0.25d0*DJ&
								+V)
						endif
					case default
						if(inb==1.and.ii<3) then
							El=El+real(0.25d0*DJ&
								-0.25d0*DJ&
								+V)
						endif
					end select
				enddo
			enddo
		enddo
		El=(El-(4d0*V-DJ)*(ne-Ns2/2))/Ns2
	end subroutine
	subroutine spin_order(cfg,nd,lphy)
		real(8) :: lphy
		integer :: cfg(:),i,x1(2),nd
		lphy=0d0
		do i=1,ne-nd*2
			call latt_one2two(cfg(i),Ns,x1)
			!lphy=lphy+sign(-1d0,i-ne2-0.1d0)*(0.5d0-mod(x1(1)+x1(2),2))
			lphy=lphy+(1-(i-1)/(ne2-nd)*2)*(0.5d0-mod(sum(x1),2))
		enddo
		lphy=lphy/Ns2
	end subroutine
	!subroutine charge_order(cfg,icfg,nd,wf,A,iA,lphy)
		!complex(8) :: A(:,:),iA(:,:),vu(Ns2,2),Y(2,2),pb,wf(:,:),a11
		!complex(8) :: lphy
		!integer :: cfg(:),icfg(:),i,n,m,j,k,l,inb,sg,cr(2),x1(2),x2(2),nr(2),dr(2),r(2),nd
		!lphy=0d0
		!lphy=lphy/Ns2
	!end subroutine
	subroutine ddw_order(cfg,icfg,nd,wf,A,iA,ja,lphy)
		complex(8) :: A(:,:),iA(:,:),vu(Ns2,2),Y(2,2),pb,wf(:,:),a11
		real(8) :: lphy,ja(:)
		integer :: cfg(:),icfg(:),i,n,m,j,k,l,inb,sg,cr(2),x1(2),nd,ud
		lphy=0d0
		do n=1,ne
			ud=1-(n-1)/ne2*2
			i=n-nd*(1-ud)/2
			call latt_one2two(cfg(i),Ns,x1)
			do inb=1,4
				j=icfg(neb(cfg(i),inb,1))
				call get_case(i,j,ud,nd,sg)
				select case(sg) 
				case(1)
					call getdiff(wf,i,j,cfg(i),cfg(j),sg,vu,cr)
					call det(vu,cr,A,iA,Y,pb)
					lphy=lphy+imag(conjg(pb))*&
						(0.5d0-mod(inb,2))*(0.5d0-mod(sum(x1),2))
				case(5)
					call getdiff(wf,i,j,cfg(i),cfg(j),sg,vu,cr)
					call det(vu,cr,A,iA,Y,pb)
					lphy=lphy-imag(conjg(pb))*&
						(0.5d0-mod(inb,2))*(0.5d0-mod(sum(x1),2))
				end select
			enddo
		enddo
		lphy=lphy/Ns2
	end subroutine
	!subroutine dsc_order(cfg,icfg,nd,wf,A,iA,ja,lphy)
		!complex(8) :: A(:,:),iA(:,:),vu(Ns2,2),Y(2,2),pb,wf(:,:),a11
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
		!complex(8) :: A(:,:),iA(:,:),vu(Ns2,2),Y(2,2),pb,wf(:,:),a11
		!real(8) :: lphy
		!integer :: cfg(:),icfg(:),i,n,m,j,k,l,inb,sg,cr(2),x1(2),x2(2),nr(2),dr(2),r(2),nd
	!end subroutine
	subroutine dsc_corelation(cfg,icfg,nd,wf,A,iA,ja,lphy)
		complex(8) :: A(:,:),iA(:,:),vu(Ns2,2),Y(2,2),pb,wf(:,:),a11
		real(8) :: lphy,ja(:)
		integer :: cfg(:),icfg(:),i,j,ii,jj,inb1,inb2,cr(2),x1(2),x2(2),nr(2),dr(2),r(2),nd
		lphy=0d0
		do i=ne2-nd+1,ne-nd
			do inb1=1,4,1
				j=icfg(neb(cfg(i),inb1,1))
				if(j<=ne2) then
					do jj=ne-nd+1,Ns2
						do inb2=1,4,1
							ii=icfg(neb(cfg(jj),inb2,1))
							if(ii>(ne-nd)) then
								cr=(/ii,j/)
							else
								cycle
							endif
							vu(:,1)=wf(cfg(i)+Ns2,:)
							vu(:,2)=wf(cfg(jj),:)
							call det(vu,cr,A,iA,Y,pb)
							lphy=lphy+real(dconjg(pb))*&
								(0.5d0-mod(inb1+inb2,2))
						enddo
					enddo
				endif
			enddo
		enddo
		lphy=lphy/(Ns2*Ns2)
	end subroutine
end module

module M_hubbard
	use M_rd
	use M_mc_matrix
	implicit none
	real(8), parameter :: U=10d0,V=0d0
contains
	subroutine two(i,j,nd,sg)
		integer :: i,j,sg,n,n1,n0,nd
		n1=ne2-nd
		n0=Ns2-ne+nd
		call random_number(n,ne2*(Ns2-ne2)*2)
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
		elseif(n<=ne2*(Ns2-ne2)) then
			n=n-n1*n0-nd*n1-nd*n0
			sg=4
			i=(n-1)/n1+1
			j=ne2+mod(n-1,n1)+1
		elseif(n<=n1*n0+ne2*(Ns2-ne2)) then
			n=n-ne2*(Ns2-ne2)
			sg=5
			i=ne2+(n-1)/n0+1
			j=ne-nd+mod(n-1,n0)+1
		elseif(n<=n1*n0+nd*n1+ne2*(Ns2-ne2)) then
			n=n-ne2*(Ns2-ne2)-n1*n0
			sg=6
			j=(n-1)/nd+1
			i=n1+mod(n-1,nd)+1
		elseif(n<=n1*n0+nd*n1+nd*n0+ne2*(Ns2-ne2)) then
			n=n-ne2*(Ns2-ne2)-n1*n0-nd*n1
			sg=7
			i=n1+(n-1)/n0+1
			j=ne-nd+mod(n-1,n0)+1
		elseif(n<=ne2*(Ns2-ne2)*2) then
			n=n-ne2*(Ns2-ne2)-n1*n0-nd*n1-nd*n0
			sg=8
			j=mod(n-1,n1)+1
			i=ne2+(n-1)/n1+1
		endif
	end subroutine
	subroutine energy(cfg,icfg,nd,wf,ja,A,iA,El)
		complex(8) :: A(:,:),iA(:,:),vu(size(A,1),2),Y(2,2),pb,wf(:,:)
		real(8) :: El
		real(8) :: ja(:)
		integer :: cfg(:),icfg(:),i,ii,j,l,inb,sg,cr(2),nd,n,ud
		El=0d0
		do n=1,ne
			ud=1-(n-1)/ne2*2
			i=n-nd*(1-ud)/2
			do ii=1,4
				do inb=1,size(t)
					j=icfg(neb(cfg(i),ii,inb))
					call get_case(i,j,ud,nd,sg)
					!if(cfg(i)==1.or.cfg(j)==1) then
					!write(*,*)i,j,nd,sg
					!endif
					select case(sg) 
					case(1:4)
						call getdiff(wf,i,j,cfg(i),cfg(j),sg,vu,cr)
						call det(vu,cr,A,iA,Y,pb)
						El=El-real(t(inb)*conjg(pb)*jast(cfg(i),cfg(j),cfg,icfg,nd,ja))
					case(5:8)
						call getdiff(wf,i,j,cfg(i),cfg(j),sg,vu,cr)
						call det(vu,cr,A,iA,Y,pb)
						El=El+real(t(inb)*conjg(pb)*jast(cfg(i),cfg(j),cfg,icfg,nd,ja))
					end select
				enddo
			enddo
		enddo
		El=El+U*nd
		El=El/Ns2
		!stop
	end subroutine
	subroutine spin_order(cfg,nd,lphy)
		real(8) :: lphy
		integer :: cfg(:),i,x1(2),nd
		lphy=0d0
		do i=1,ne-nd*2
			call latt_one2two(cfg(i),Ns,x1)
			!lphy=lphy+sign(-1d0,i-ne2-0.1d0)*(0.5d0-mod(x1(1)+x1(2),2))
			lphy=lphy+(1-(i-1)/(ne2-nd)*2)*(0.5d0-mod(sum(x1),2))
		enddo
		lphy=lphy/Ns2
	end subroutine
	!subroutine charge_order(cfg,icfg,nd,wf,A,iA,lphy)
		!complex(8) :: A(:,:),iA(:,:),vu(Ns2,2),Y(2,2),pb,wf(:,:),a11
		!complex(8) :: lphy
		!integer :: cfg(:),icfg(:),i,n,m,j,k,l,inb,sg,cr(2),x1(2),x2(2),nr(2),dr(2),r(2),nd
		!lphy=0d0
		!lphy=lphy/Ns2
	!end subroutine
	subroutine ddw_order(cfg,icfg,nd,wf,A,iA,ja,lphy)
		complex(8) :: A(:,:),iA(:,:),vu(Ns2,2),Y(2,2),pb,wf(:,:),a11
		real(8) :: lphy,ja(:)
		integer :: cfg(:),icfg(:),i,n,m,j,k,l,inb,sg,cr(2),x1(2),nd,ud
		lphy=0d0
		do n=1,ne
			ud=1-(n-1)/ne2*2
			i=n-nd*(1-ud)/2
			call latt_one2two(cfg(i),Ns,x1)
			do inb=1,4
				j=icfg(neb(cfg(i),inb,1))
				call get_case(i,j,ud,nd,sg)
				select case(sg) 
				case(1:4)
					call getdiff(wf,i,j,cfg(i),cfg(j),sg,vu,cr)
					call det(vu,cr,A,iA,Y,pb)
					lphy=lphy+imag(conjg(pb)*jast(cfg(i),cfg(j),cfg,icfg,nd,ja))*&
						(0.5d0-mod(inb,2))*(0.5d0-mod(sum(x1),2))
				case(5:8)
					call getdiff(wf,i,j,cfg(i),cfg(j),sg,vu,cr)
					call det(vu,cr,A,iA,Y,pb)
					lphy=lphy-imag(conjg(pb)*jast(cfg(i),cfg(j),cfg,icfg,nd,ja))*&
						(0.5d0-mod(inb,2))*(0.5d0-mod(sum(x1),2))
				end select
			enddo
		enddo
		lphy=lphy/Ns2
	end subroutine
	!subroutine dsc_order(cfg,icfg,nd,wf,A,iA,ja,lphy)
		!complex(8) :: A(:,:),iA(:,:),vu(Ns2,2),Y(2,2),pb,wf(:,:),a11
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
		!complex(8) :: A(:,:),iA(:,:),vu(Ns2,2),Y(2,2),pb,wf(:,:),a11
		!real(8) :: lphy
		!integer :: cfg(:),icfg(:),i,n,m,j,k,l,inb,sg,cr(2),x1(2),x2(2),nr(2),dr(2),r(2),nd
	!end subroutine
	subroutine dsc_corelation(cfg,icfg,nd,wf,A,iA,ja,lphy)
		complex(8) :: A(:,:),iA(:,:),vu(Ns2,2),Y(2,2),pb,wf(:,:),a11
		real(8) :: lphy,ja(:)
		integer :: cfg(:),icfg(:),i,j,ii,jj,inb1,inb2,cr(2),x1(2),x2(2),nr(2),dr(2),r(2),nd
		lphy=0d0
		do i=ne2-nd+1,ne-nd
			do inb1=1,4,1
				j=icfg(neb(cfg(i),inb1,1))
				if(j<=ne2) then
					do jj=ne2-nd+1,Ns2
						do inb2=1,4,1
							ii=icfg(neb(cfg(jj),inb2,1))
							if(ii<=(ne2-nd)) then
								cr=(/ii+ne2,j/)
							elseif(ii>(ne-nd)) then
								cr=(/ii,j/)
							else
								cycle
							endif
							vu(:,1)=wf(cfg(i)+Ns2,:)
							vu(:,2)=wf(cfg(jj),:)
							call det(vu,cr,A,iA,Y,pb)
							lphy=lphy+real(dconjg(pb)*jast(cfg(i),cfg(ii),cfg,icfg,nd,ja)*jast(cfg(j),cfg(jj),cfg,icfg,nd,ja))*&
								(0.5d0-mod(inb1+inb2,2))
						enddo
					enddo
				endif
			enddo
		enddo
		lphy=lphy/(Ns2*Ns2)
	end subroutine
end module

module M_vmc
	use M_tJ
	!use M_hubbard
	implicit none
contains
	subroutine mc(wf,ja,Nmc,cfg,nd,phyval,sga,dwf,O,S,g)
		complex(8) :: pb,A(Ns2,Ns2),iA(Ns2,Ns2),dA(Ns2,Ns2,vn),vu(Ns2,2),dvu(Ns2,2,vn),wf(:,:),Y(2,2)
		complex(8), optional :: O(:),S(:,:),g(:),dwf(:,:,:)
		real(8) :: rpb,sga(:,:),ja(:)
		complex(8) :: lO(vn+size(ja)),lS(size(lO),size(lO)),lg(size(lO))
		real(8) :: phyval(:),lphy(size(phyval))
		integer :: icfg(Ns2),i,j,ti,tj,l,sg,n,Nmc(:),cr(2),acp,bi,bj,ii,jj
		integer :: cfg(Ns2),nd
		logical :: flag
		flag=present(O)
		!write(*,"(A)")"start monte carlo"
		if(mod(Ns2,2)/=0) then
			write(*,*)"Number of site is not even, exit!!"
		endif
		bi=1
		bj=size(sga,1)/2+1
		sga=0d0
		n=0
		phyval=0d0
		acp=0
		nd=0
		do i=1,Ns2
			icfg(cfg(i))=i
		enddo
		if(flag) then
			O=0d0
			S=0d0
			g=0d0
		endif
		!get initial matrix
		do i=1,Ns2
			if(i<=ne2) then
				A(i,:)=wf(cfg(i),:)
				if(flag) then
					dA(i,:,:)=dwf(cfg(i),:,:)
				endif
			elseif(i>=ne-nd+1)  then
				A(i,:)=wf(cfg(i)+Ns2,:)
				if(flag) then
					dA(i,:,:)=dwf(cfg(i)+Ns2,:,:)
				endif
			else
				A(i,:)=wf(cfg(i-ne2)+Ns2,:)
				if(flag) then
					dA(i,:,:)=dwf(cfg(i-ne2)+Ns2,:,:)
				endif
			endif
		enddo
		iA=A
		call matrix_inv(iA)
		!write(*,*)sum(abs(matmul(A,iA)-diag(1d0,size(A,1))))
		do
			n=n+1
			call two(i,j,nd,sg)
			!write(*,*)i,j,nd,sg
			call getdiff(wf,i,j,cfg(i),cfg(j),sg,vu,cr,flag,dwf,dvu)
			call det(vu,cr,A,iA,Y,pb)
			call random_number(rpb)
			if(rpb<real(pb*conjg(pb))*(jast(cfg(i),cfg(j),cfg,icfg,nd,ja))**2) then
				acp=acp+1
				call update_matrix(vu,cr,Y,pb,A,iA,flag,dvu,dA)
				!write(*,*)sum(abs(matmul(A,iA)-diag(1d0,size(A,1))))
				call update_swap(i,j,nd,cfg,icfg,A,iA,dA,sg)
				!if(sum(abs(matmul(A,iA)-diag(1d0,size(A,1))))>1d-5) then
					!write(*,*)i,j,sg,nd,sum(abs(matmul(A,iA)-diag(1d0,size(A,1))))
					!stop
				!endif
				if(mod(n,500)==0) then
					iA=A
					call matrix_inv(iA)
				endif
				!if(n>Nmc(1)) then
					!call checkcfg(icfg,nd,100)
				!endif
			endif
			if(n>Nmc(1).and.mod(n-Nmc(1),Nmc(2))==0) then
				call energy(cfg,icfg,nd,wf,ja,A,iA,lphy(1))
				if(flag) then
					do ti=1,vn
						lO(ti)=Tr(iA,dA(:,:,ti))
					enddo
					do ti=vn+1,size(lO)
						lO(ti)=-nd
					enddo
					do ti=1,size(lO)
						do tj=1,size(lO)
							lS(ti,tj)=dconjg(lO(ti))*lO(tj)
						enddo
					enddo
					lg=real(lphy(1)*lO)
					S=S+lS
					g=g+lg
					O=O+lO
				else
					call spin_order(cfg,nd,lphy(2))
					call ddw_order(cfg,icfg,nd,wf,A,iA,ja,lphy(3))
					call dsc_corelation(cfg,icfg,nd,wf,A,iA,ja,lphy(4))
				endif
				call binning(bi,bj,sga,lphy)
				phyval=phyval+lphy
			endif
			if(n>=(Nmc(1)+Nmc(2)*Nmc(3))) then
				phyval=phyval/Nmc(3)
				if(flag) then
					S=S/Nmc(3)
					g=g/Nmc(3)
					O=O/Nmc(3)
				endif
				exit
			endif
		enddo
		!write(*,"('accept/total number is ',I6,'/',I6)")acp,n
	end subroutine
	subroutine variational(var,Nmc,cfg,nd)
		use lapack95, only: heevd,heevr
		use M_wf, iniwf => wf_r
		real(8) :: var(:),dvar(size(var)),pvar(size(var)),eg(size(var)),dt=0.1d0,er(1),Ev,allE(200),scv,var_av(size(var))
		complex(8) :: wf(Ns2*2,Ns2),dwf(Ns2*2,Ns2,vn),O(size(var)),S(size(O),size(O)),g(size(O)),Ov(size(O)),Sv(size(O),size(O)),gv(size(O))
		integer :: n,i,j,k,info,Nmc(:),nm,l,nav
		integer :: cfg(Ns2),nd,cfg_omp(Ns2),nd_omp
		real(8) :: phyval(1)
		real(8), allocatable :: sga(:,:)
		logical :: flag
		nm=log(real(Nmc(3)))/log(2d0)
		allocate(sga((nm+1)*2,1))
		n=24
		l=0
		scv=0d0
		flag=.true.
		pvar=var
		nav=0
		var_av=0d0
		do
			l=l+1
			write(*,"(I3$)")ne
			write(*,"(es9.2$)")var
			call iniwf(var(:vn),wf,dwf)
			do i=1,Ns2
				cfg(i)=i
			enddo
			call fisher_yates_shuffle(cfg,Ns2)
			Ev=0d0
			Sv=0d0
			gv=0d0
			er=0d0 
			!$OMP PARALLEL DO REDUCTION(+:Ev,Sv,gv,er) PRIVATE(phyval,O,S,g) LASTPRIVATE(sga,cfg_omp,nd_omp) SCHEDULE(GUIDED)
			do k=1,n
				cfg_omp=cfg
				nd_omp=nd
				call mc(wf,var(vn+1:),Nmc,cfg_omp,nd_omp,phyval,sga,dwf,O,S,g)
				Ev=Ev+phyval(1)
				do i=1,size(O)
					do j=1,size(O)
						Sv(i,j)=Sv(i,j)+S(i,j)-dconjg(O(i))*O(j)
					enddo
				enddo
				gv=gv+g-real(phyval(1)*O)
				er=er+abs(phyval)**2
			enddo
			!$OMP END PARALLEL DO
			cfg=cfg_omp
			nd=nd_omp
			Ev=Ev/n
			Sv=Sv/n
			gv=2d0*gv/n*Ns2
			!call mwrite(10,real(Sv))
			!write(*,*)real(gv)
			!stop
			er=sqrt(abs(er/n-abs(Ev)**2))
			!Sv=Sv+diag(2d-1,size(Sv,1))
			Sv=Sv+diag(1d-2,size(Sv,1))
			!call zero(Sv,gv,(/1,2,5,6/))
			call heevd(Sv,eg,"V")
			do i=1,size(O)
				dvar(i)=0d0
				do j=1,size(O)
					do k=1,size(dvar)
						if(eg(k)/eg(size(O))<1d-3) then
							cycle
						endif
						dvar(i)=dvar(i)+real(Sv(i,k)*dconjg(Sv(j,k))*gv(j)/eg(k))
					enddo
				enddo
			enddo
			if(n<5) then
				er=sqrt((sga(nm-3,1)-sga(size(sga,1)/2,1))/(2**(nm-(nm-3)+1)-1))
			endif
			write(*,"(es9.2$)")real(gv),Ev,er
			write(20,"(I4$)")ne
			write(20,"(es13.5$)")var,real(gv),Ev,er
			write(20,"(x)")
			allE(l)=Ev
			if(allE(max(l-1,1))<(Ev-2d0*er(1))) then
				write(*,"(x)")
				if(flag) then
					l=0
					var(5:6)=0d0
					pvar=var
					flag=.false.
					cycle
				else
					write(20,"(x)")
					var=pvar
					exit
				endif
			endif
			do i=max(l-39,1),l-1
				scv=scv+sign(1d0,Ev-allE(i))
				if(max(l-39,1)>1) then
					scv=scv+sign(1d0,allE(max(l-39,1)-1)-allE(i))
				endif
			enddo
			write(*,"(es9.2$)")scv/(min(l*max(l-1,1),40*39))
			if(l>40.and.abs(scv)/(min(l*max(l-1,1),40*39))<1d-2.and.scv<0.or.l==size(allE).or.nav>0) then
				nav=nav+1
				var_av=var_av+var
				write(*,*)"variational finish",l,nav
				write(20,"(x)")
				!write(10,"(I4$)")ne
				!write(10,"(es13.5$)")var,Ev,er
				!write(10,"(x)")
				if(nav>100) then
					var=var_av/nav
					exit
				endif
			endif
			pvar=var
			var=var-dvar*dt
			write(*,"(x)")
			!Nmc(1)=500
		enddo
	end subroutine
	subroutine zero(Sv,gv,o)
		integer :: o(:),i
		complex(8) :: Sv(:,:),gv(:)
		do i=1,size(o)
			Sv(:,o(i))=0d0
			Sv(o(i),:)=0d0
			gv(o(i))=0d0
		enddo
	end subroutine
end module

program main
	use M_vmc
	use M_wf, iniwf => wf_r
	implicit none
	integer :: n,i,j,Nmc(3),Nvar(3),nm,cfg(Ns2),nd
	complex(8) :: wf(Ns2*2,Ns2)
	real(8) :: var(vn+1)
	real(8) :: er(5),phyval(5)
	real(8), allocatable :: sga(:,:)
	logical :: f
	!open(10,file="../data/2d.dat")
	f=openfile(20,"../data/var.dat")
	!f=openfile(40,"../data/2d.dat")
	!f=openfile(10,"../data/err.dat")
	!f=openfile(30,"../data/phyvar.dat",access="direct",recl=226)
	f=openfile(30,"../data/phyvar.dat")
	Nmc(3)=1024*4
	Nvar(3)=128
	nm=log(real(Nmc(3)))/log(2d0)
	allocate(sga((nm+1)*2,size(phyval)))
	call latt(Ns,Tx,Ty,neb)
	call init_random_seed()
	nd=0
	n=32
	!var=(/1d-1,1d-1,1d-1,1d-1,1d-1,1d-1/)
	!var=(/1d-1,0d0,1d-1,1d-1,0d0,0d0,1d-1/)
	var=(/1.09E-01,-8.23E-02,-1.84E-03,-1.74E-02,1.84E-01,1.28E-01,1.00E-01/)
	do i=10,16,1
		!var(7)=var(7)+0.1
		!var=(/1d-2,0d0,1d-2,1d-2,0d0,0d0/)
		!var=(/1d-1,0d0,1d-1,1d-1,0d0,0d0/)
		ne=int(Ns2*(40-i)/40d0)
		ne=ne+mod(ne,2)
		!read(30,rec=i+1,fmt="(i4,7es13.5)")ne,var
		ne2=ne/2
		!call self_n(var)
		!Nmc(1:2)=(/5000*Ns2,5*Ns2/)*(i+1)
		!Nvar(1:2)=(/500*Ns2,5*Ns2/)*(i+1)
		!Nvar(1:2)=(/500,5/)*Ns2
		Nvar(1:2)=(/1000,5/)*100
		Nmc(1:2)=(/1000,5/)*100
		call variational(var,Nvar,cfg,nd)
		write(*,"(i4$)")ne
		write(*,"(es9.2$)")var
		call iniwf(var(:vn),wf)
		call mc(wf,var(vn+1:),Nmc,cfg,nd,phyval,sga)
		!do j=1,size(sga,1)/2-1
			!write(10,"(es13.5$)")(sqrt((sga(j,k)-sga(size(sga,1)/2,k))/(2**(nm-j+1)-1)),sqrt((sga(j,k)-sga(size(sga,1)/2,k))/(2**(nm-j+1)-1)*sqrt(2d0/(2**(nm-j+1)-1))),k=1,size(sga,2))
			!write(10,"(x)")
		!enddo
		!write(*,"(es9.2$)")phyval(1),sqrt((sga(nm-5,1)-sga(size(sga,1)/2,1))/(2**(nm-(nm-5)+1)-1))
		write(*,"(es9.2$)")phyval
		write(*,"(1x)")
		!write(30,rec=i+1,fmt="(i4,17es13.5,A)")ne,var,(phyval(j),sqrt((sga(nm-5,j)-sga(size(sga,1)/2,j))/(2**(nm-(nm-5)+1)-1)),j=1,5),char(10)
		write(30,"(es13.5$)")real(Ns2-ne)/Ns2,var,(phyval(j),sqrt((sga(nm-5,j)-sga(size(sga,1)/2,j))/(2**(nm-(nm-5)+1)-1)),j=1,5)
		write(30,"(1x)")
	enddo
end program
