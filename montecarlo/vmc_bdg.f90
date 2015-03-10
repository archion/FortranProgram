module M_pmt
	use M_const
	implicit none
	integer, parameter :: Ns(2)=(/9,10/),Ns2=Ns(1)**2+1,Tx(2)=(/Ns(1),-1/),Ty(2)=(/1,Ns(1)/),vn=6
	!integer, parameter :: Ns(2)=(/10,10/),Ns2=Ns(1)**2,Tx(2)=(/Ns(1),0/),Ty(2)=(/0,Ns(1)/),vn=6
	integer :: neb(Ns2,4,3),ne=Ns2,ne2=Ns2/2
	real(8), parameter :: t(1)=(/1d0/)
end module

module M_wf
	use M_utility
	use M_pmt
	use M_matrix
	use M_latt, only : &
		 latt         =>  square_tilda         ,& 
		 latt_bc      =>  square_tilda_bc      ,& 
		 latt_one2two =>  square_tilda_one2two ,& 
		 latt_two2one =>  square_tilda_two2one 
		 !latt         =>  square         ,& 
		 !latt_bc      =>  square_bc      ,& 
		 !latt_one2two =>  square_one2two ,& 
		 !latt_two2one =>  square_two2one 
contains
	subroutine wf_r(var,wf,dwf)
		use lapack95, only: heevd,heevr
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
		complex(8) :: H(:,:),D(:,:,:),dH(:,:,:),X(size(H,2)),C(size(H,1),size(H,2))
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
			do i=1,size(H,2)
				do j=1,size(H,1)
					if(abs(eg(i)-eg(j))<1d-10) then
						cycle
					endif
					C(j,i)=dot_product(H(:,j),matmul(D(:,:,n),H(:,i)))*X(j)*conjg(X(j))/((eg(i)-eg(j)))
				enddo
				C(i,i)=-sum(H(m(i),:)*C(:,i))
			enddo
			dH(:,:,n)=matmul(H,C)
		enddo
	end subroutine
	subroutine bdg(var,H,D)
		complex(8) :: H(Ns2*2,Ns2*2),D(:,:,:)
		real(8) :: eg(Ns2*2),var(vn)
		integer :: i,j,n,ix(2),bc
		H=0d0
		wf=0d0
		D=0d0
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
				D(i,neb(i,j,1),4)=(img*(-1)**mod(sum(ix),2)*(-1)**mod(j,2))*bc
				D(i+Ns2,neb(i,j,1)+Ns2,4)=(img*(-1)**mod(sum(ix),2)*(-1)**mod(j,2))*bc
				D(i,neb(i,j,2),5)=-bc
				D(i+Ns2,neb(i,j,2)+Ns2,5)=-D(i,neb(i,j,2),5)
				D(i,neb(i,j,3),6)=-bc
				D(i+Ns2,neb(i,j,3)+Ns2,6)=-D(i,neb(i,j,3),6)
				D(i,neb(i,j,1)+Ns2,1)=(-1)**mod(j,2)*bc
				D(i+Ns2,neb(i,j,1),1)=conjg(D(i,neb(i,j,1)+Ns2,1))
			enddo
			! on site
			H(i,i)=var(3)*(-1)**mod(sum(ix),2)-var(2)
			H(i+Ns2,i+Ns2)=var(3)*(-1)**mod(sum(ix),2)+var(2)
			D(i,i,3)=(-1)**mod(sum(ix),2)
			D(i+Ns2,i+Ns2,3)=(-1)**mod(sum(ix),2)
			D(i,i,2)=-1d0
			D(i+Ns2,i+Ns2,2)=1d0
		enddo
	end subroutine
end module

module M_mc_matrix
	use M_matrix
	use M_pmt
	use M_utility
	implicit none
contains
	subroutine getdiff(wf,ii,jj,cfgii,cfgjj,sg,vu,cr,flag,dwf,dvu)
		integer :: i,j,sg,cfgi,cfgj,cr(:),ii,jj,cfgii,cfgjj
		complex(8) :: vu(:,:),wf(:,:)
		complex(8), optional :: dvu(:,:,:),dwf(:,:,:)
		logical, optional :: flag
		if(ii>jj) then
			i=jj
			j=ii
			cfgi=cfgjj
			cfgj=cfgii
		else
			i=ii
			j=jj
			cfgi=cfgii
			cfgj=cfgjj
		endif
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
			cr=(/i+ne2,-1/)
			if(present(flag).and.flag) then
				dvu(:,1,:)=dwf(cfgj+Ns2,:,:)
			endif
			vu(:,1)=wf(cfgj+Ns2,:)
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
			call swap(cfg(i),cfg(ne2-nd))
			call swap(icfg(cfg(i)),icfg(cfg(ne2-nd)))
			call swap(cfg(j),cfg(ne-nd))
			call swap(icfg(cfg(j)),icfg(cfg(ne-nd)))
			call swap(A(i,:),A(ne2-nd,:))
			call swap(A(i+ne2,:),A(ne-nd,:))
			call swap(iA(:,i),iA(:,ne2-nd))
			call swap(iA(:,i+ne2),iA(:,ne-nd))
			do ti=1,vn
				call swap(dA(i,:,ti),dA(ne2-nd,:,ti))
				call swap(dA(i+ne2,:,ti),dA(ne-nd,:,ti))
			enddo
			nd=nd+1
		end select
	end subroutine
	subroutine det(vu,cr,A,iA,Y,pb,sg)
		complex(8) :: vu(:,:),A(:,:),iA(:,:),Y(:,:),pb
		integer :: cr(:),sg
		select case(sg) 
		case(1:8) 
			call det_ratio_row(vu(:,1),cr(1),A,iA,pb)
		case(0)
			call det_ratio_tworow(vu,cr,A,iA,Y,pb)
		end select
	end subroutine
	subroutine update_matrix(vu,cr,Y,pb,A,iA,sg,flag,dvu,dA)
		complex(8) :: vu(:,:),A(:,:),iA(:,:),Y(:,:),pb
		complex(8) :: dvu(:,:,:),dA(:,:,:)
		integer :: cr(:),sg
		logical :: flag
		select case(sg) 
		case(1:8)
			call inv_update_row(vu(:,1),cr(1),pb,A,iA)
			if(flag) then
				dA(cr(1),:,:)=dvu(:,1,:)
			endif
		case(0)
			call inv_update_tworow(vu,cr,Y,A,iA)
			if(flag) then
				dA(cr(1),:,:)=dvu(:,1,:)
				dA(cr(2),:,:)=dvu(:,2,:)
			endif
		end select
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
	function jast(ri,rj,icfg,nd,sg,ja)
		integer :: ri,rj,icfg(:),sg,i,nd
		real(8) :: jast,ja(:)
		jast=1d0
		select case(sg)
		case(3,7)
			jast=jast*exp(ja(1))
			!jast=(Ns2-(1d0-jvar(1))*(nd-1))/(Ns2-(1d0-jvar(1))*nd)
		case(4,8)
			jast=jast*exp(-ja(1))
			!jast=(Ns2-(1d0-jvar(1))*(nd+1))/(Ns2-(1d0-jvar(1))*nd)
		end select
		do i=2,size(ja)
		enddo
	end function
end module

module M_tJ
	use M_pmt
	use M_rd
	use M_mc_matrix
	implicit none
	real(8), parameter :: DJ=0.4d0,V=0d0
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
		complex(8) :: El
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
						call det(vu,cr,A,iA,Y,pb,sg)
						El=El-t(inb)*conjg(pb)
					case(5) 
						call getdiff(wf,i,j,cfg(i),cfg(j),sg,vu,cr)
						call det(vu,cr,A,iA,Y,pb,sg)
						El=El+t(inb)*conjg(pb)
					case(0)
						if(inb==1.and.i<=ne2) then
							call getdiff(wf,i,j,cfg(i),cfg(j),sg,vu,cr)
							call det(vu,cr,A,iA,Y,pb,sg)
							El=El+0.5d0*DJ*conjg(pb)-0.25d0*DJ&
								-0.25d0*DJ&
								+V
						endif
					case default
						if(inb==1.and.ii<3) then
							El=El+0.25d0*DJ&
								-0.25d0*DJ&
								+V
						endif
					end select
				enddo
			enddo
		enddo
		El=El/Ns2
	end subroutine
end module

module M_hubbard
	use M_pmt
	use M_rd
	use M_mc_matrix
	implicit none
	real(8), parameter :: U=10d0
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
			i=(n-1)/nd+1
			j=n1+mod(n-1,nd)+1
		elseif(n<=n1*n0+nd*n1+nd*n0+ne2*(Ns2-ne2)) then
			n=n-ne2*(Ns2-ne2)-n1*n0-nd*n1
			sg=7
			i=n1+(n-1)/n0+1
			j=ne-nd+mod(n-1,n0)+1
		elseif(n<=ne2*(Ns2-ne2)*2) then
			n=n-ne2*(Ns2-ne2)-n1*n0-nd*n1-nd*n0
			sg=8
			i=mod(n-1,n1)+1
			j=ne2+(n-1)/n1+1
		endif
	end subroutine
	subroutine energy(cfg,icfg,nd,wf,ja,A,iA,El)
		complex(8) :: A(:,:),iA(:,:),vu(size(A,1),2),Y(2,2),pb,wf(:,:)
		complex(8) :: El
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
						call det(vu,cr,A,iA,Y,pb,sg)
						El=El-t(inb)*conjg(pb)*jast(cfg(i),cfg(j),icfg,nd,sg,ja)
					case(5:8)
						call getdiff(wf,i,j,cfg(i),cfg(j),sg,vu,cr)
						call det(vu,cr,A,iA,Y,pb,sg)
						El=El+t(inb)*conjg(pb)*jast(cfg(i),cfg(j),icfg,nd,sg,ja)
					end select
				enddo
			enddo
		enddo
		El=El+U*nd
		El=El/Ns2
		!stop
	end subroutine
end module

module M_vmc
	use M_wf, iniwf => wf_r
	!use M_tJ
	use M_hubbard
	implicit none
contains
	subroutine mc(wf,ja,Nmc,phyval,sga,dwf,O,S,g)
		complex(8) :: pb,A(Ns2,Ns2),iA(Ns2,Ns2),dA(Ns2,Ns2,vn),vu(Ns2,2),dvu(Ns2,2,vn),wf(:,:),Y(2,2)
		complex(8), optional :: O(:),S(:,:),g(:),dwf(:,:,:)
		real(8) :: rpb,sga(:,:),ja(:)
		complex(8) :: lO(vn+size(ja)),lS(size(lO),size(lO)),lg(size(lO))
		complex(8) :: phyval(:),lphy(size(phyval))
		integer :: cfg(Ns2),icfg(Ns2),i,j,ti,tj,l,sg,n,Nmc(:),cr(2),acp,bi,bj,ii,jj,nd
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
			cfg(i)=i
		enddo
		call fisher_yates_shuffle(cfg,Ns2)
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
			call det(vu,cr,A,iA,Y,pb,sg)
			call random_number(rpb)
			if(rpb<real(pb*conjg(pb))*(jast(cfg(i),cfg(j),icfg,nd,sg,ja))**2) then
				acp=acp+1
				call update_matrix(vu,cr,Y,pb,A,iA,sg,flag,dvu,dA)
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
					!call phymeasure(cfg,icfg,wf,A,iA,lphy(2:5))
				endif
				call binning(bi,bj,sga,(/real(lphy(1)),real(lphy(2)),imag(lphy(3)),real(lphy(1)),real(lphy(1))/))
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
	!subroutine phymeasure(cfg,icfg,wf,A,iA,lod)
		!complex(8) :: A(:,:),iA(:,:),vu(Ns2,2),Y(2,2),pb,wf(:,:),a11
		!complex(8) :: lod(:)
		!integer :: cfg(:),icfg(:),i,n,m,j,k,l,inb,sg,cr(2),x1(2),x2(2),nr(2),dr(2),r(2)
		!lod=0d0
		!!SDW Order
		!do i=1,ne
			!call latt_one2two(cfg(i),Ns,x1)
			!lod(1)=lod(1)+sign(-1d0,i-ne2-0.1d0)*(0.5d0-mod(x1(1)+x1(2),2))
		!enddo
		!!DDW Order
		!do j=ne+1,Ns2
			!call latt_one2two(cfg(j),Ns,x1)
			!do inb=1,4
				!i=icfg(neb(cfg(j),inb,1))
				!call getdiff(wf,i,j,cfg(i),cfg(j),sg,vu,cr)
				!call det(vu,cr,A,iA,Y,pb,sg)
				!select case(sg) 
				!case(1) 
					!lod(2)=lod(2)+conjg(pb)*(0.5d0-mod(inb,2))*(0.5d0-mod(sum(x1),2))
				!case(2) 
					!lod(2)=lod(2)-conjg(pb)*(0.5d0-mod(inb,2))*(0.5d0-mod(sum(x1),2))
				!end select
			!enddo
		!enddo
		!!d-SC Order
		!!do i=ne+1,Ns2
			!!do inb=1,4,1
				!!j=icfg(neb(cfg(i),inb,1))
				!!if(j>ne) then
					!!vu(:,1)=wf(cfg(i),:)
					!!cr=(/j,-1/)
					!!call det(vu,cr,A,iA,Y,pb,1)
					!!lod(3)=lod(3)+conjg(pb)*(0.5d0-mod(inb,2))
				!!endif
			!!enddo
		!!enddo
		!!d-SC corelation
		!do i=ne+1,Ns2
			!do inb=1,4,1
				!j=icfg(neb(cfg(i),inb,1))
				!if(j>ne) then
					!do k=1,ne2
						!do n=1,4,1
							!l=icfg(neb(cfg(k),n,1))
							!if(l>ne2.and.l<=ne) then
								!vu(:,1)=wf(cfg(i),:)
								!vu(:,2)=wf(Ns2+cfg(l),:)
								!cr=(/k,j/)
								!call det(vu,cr,A,iA,Y,pb,3)
								!lod(4)=lod(4)+dconjg(pb)*(0.5d0-mod(inb+n,2))
							!endif
						!enddo
					!enddo
				!endif
			!enddo
		!enddo
		!lod(4)=lod(4)/Ns2
		!lod=lod/Ns2
	!end subroutine
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
	subroutine variational(var,Nmc)
		use lapack95, only: heevd,heevr
		real(8) :: var(:),dvar(size(var)),pvar(size(var)),eg(size(var)),dt=1d0,er(1),Ev,allE(400),scv
		complex(8) :: wf(Ns2*2,Ns2),dwf(Ns2*2,Ns2,vn),O(size(var)),S(size(O),size(O)),g(size(O)),Ov(size(O)),Sv(size(O),size(O)),gv(size(O))
		integer :: n,i,j,k,info,Nmc(:),nm,l
		complex(8) :: phyval(1)
		real(8), allocatable :: sga(:,:)
		logical :: flag
		nm=log(real(Nmc(3)))/log(2d0)
		allocate(sga((nm+1)*2,1))
		n=24
		l=0
		scv=0d0
		flag=.true.
		pvar=var
		do
			l=l+1
			write(*,"(I3$)")ne
			write(*,"(es9.2$)")var
			call iniwf(var(:vn),wf,dwf)
			Ev=0d0
			Sv=0d0
			gv=0d0
			er=0d0
			!$OMP PARALLEL DO REDUCTION(+:Ev,Sv,gv,er) PRIVATE(phyval,O,S,g) LASTPRIVATE(sga)
			do k=1,n
				call mc(wf,var(vn+1:),Nmc,phyval,sga,dwf,O,S,g)
				Ev=Ev+real(phyval(1))
				do i=1,size(O)
					do j=1,size(O)
						Sv(i,j)=Sv(i,j)+S(i,j)-dconjg(O(i))*O(j)
					enddo
				enddo
				gv=gv+g-real(phyval(1)*O)
				er=er+abs(phyval)**2
			enddo
			!$OMP END PARALLEL DO
			Ev=Ev/n
			Sv=Sv/n
			gv=2d0*gv/n
			!call mwrite(10,real(Sv))
			!write(*,*)real(gv)
			!stop
			er=sqrt(abs(er/n-abs(Ev)**2))
			Sv=Sv+diag(1d-2,size(Sv,1))
			!Sv(:,5:6)=0d0
			!Sv(5:6,:)=0d0
			!Sv(:,1)=0d0
			!Sv(1,:)=0d0
			!Sv(:,3:4)=0d0
			!Sv(3:4,:)=0d0
			!Sv(:,2:5)=0d0
			!Sv(2:5,:)=0d0
			call heevd(Sv,eg,"V")
			!gv(1)=0d0
			!gv(3:4)=0d0
			!gv(5:6)=0d0
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
			do i=max(l-49,1),l-1
				scv=scv+sign(1d0,Ev-allE(i))
				if(max(l-49,1)>1) then
					scv=scv+sign(1d0,allE(max(l-49,1)-1)-allE(i))
				endif
			enddo
			write(*,"(es9.2$)")scv/(min(l*max(l-1,1),50*49))
			if(l>50.and.abs(scv)/(min(l*max(l-1,1),50*49))<1d-2.or.l==size(allE)) then
				write(*,*)"variational finish",sum(abs(dvar))
				write(20,"(x)")
				write(10,"(I4$)")ne
				write(10,"(es13.5$)")var,Ev,er
				write(10,"(x)")
				exit
			endif
			pvar=var
			var=var-dvar*dt
			write(*,"(x)")
		enddo
	end subroutine
end module

program main
	use M_vmc
	use M_utility
	integer :: n,i,j,Nmc(3),Nvar(3),nm
	complex(8) :: wf(Ns2*2,Ns2)
	real(8) :: var(vn+1)
	complex(8) :: er(5),phyval(5)
	real(8), allocatable :: sga(:,:)
	logical :: f
	!open(10,file="../data/2d.dat")
	f=openfile(20,"../data/var.dat")
	f=openfile(40,"../data/2d.dat")
	f=openfile(10,"../data/err.dat")
	!f=openfile(30,"../data/phyvar.dat",access="direct",recl=278)
	f=openfile(30,"../data/phyvar.dat")
	Nmc(3)=1024*4
	Nvar(3)=128
	nm=log(real(Nmc(3)))/log(2d0)
	allocate(sga((nm+1)*2,size(phyval)))
	call init_random_seed()
	call latt(Ns,Tx,Ty,neb)
	n=32
	!var=(/1d-1,1d-1,1d-1,1d-1,1d-1,1d-1/)
	!var=(/9.46E-02,-4.62E-01,5.37E-03,5.19E-12,8.10E-01,4.40E-01/)
	!var=(/9.96E-02,-5.21E-01,4.42E-03,6.81E-12,7.56E-01,3.43E-01/)
	!var=(/9.96E-02,-4.62E-01,5.37E-03,5.19E-12,0.00E-01,3.43E-01/)
	!var=(/1.88764E-01,-4.24821E-02,2.34273E-02,4.23654E-04,5.68984E-02,8.57877E-02/)
	!var=(/1.55689E-01,-5.94975E-02,3.20398E-03,2.28620E-05,8.86424E-02,9.65006E-02/)
	!var=1d-10
	!var(2)=1d0
	!var=(/0.18764E-01,-4.24821E-02,2.34273E-02,4.23654E-04,5.68984E-02,8.57877E-02/)
	var=(/1d-5,0d0,1.5d0,0d0,0d0,0d0,1d-1/)
	do i=0,40,1
		!var(7)=var(7)+0.1
		!var=(/1d-2,0d0,1d-2,1d-2,0d0,0d0/)
		!var=(/1d-1,0d0,1d-1,1d-1,0d0,0d0/)
		ne=Ns2-2*i
		!read(21,rec=i+1,fmt="(i4,6es13.5$)")ne,var
		!read(30,rec=i+1,fmt="(i4,6es13.5)")ne,var
		!write(*,"(i4,6es13.5)")ne,var
		!read(30,"(6es13.5)",advance="no")var
		!write(30,"(+,6es13.5)")1d0,1d0,1d0,1d0,1d0,1d0
		ne2=ne/2
		!Nmc(1:2)=(/5000*Ns2,5*Ns2/)*(i+1)
		!Nvar(1:2)=(/500*Ns2,5*Ns2/)*(i+1)
		Nmc(1:2)=(/5000*Ns2,5*Ns2/)
		Nvar(1:2)=(/500*Ns2,5*Ns2/)
		!write(30,"(i4$)")ne
		call variational(var,Nvar)
		write(*,"(i4$)")ne
		write(*,"(es9.2$)")var
		!write(30,"(es13.5$)")var
		!write(30,"(1x)")
		!write(*,"(1x)")
		!cycle
		call iniwf(var(:vn),wf)
		call mc(wf,var(vn+1:),Nmc,phyval,sga)
		do j=1,size(sga,1)/2-1
			write(10,"(es13.5$)")(sqrt((sga(j,k)-sga(size(sga,1)/2,k))/(2**(nm-j+1)-1)),sqrt((sga(j,k)-sga(size(sga,1)/2,k))/(2**(nm-j+1)-1)*sqrt(2d0/(2**(nm-j+1)-1))),k=1,size(sga,2))
			write(10,"(x)")
		enddo
		!write(*,"(es9.2$)")phyval(1),sqrt((sga(nm-5,1)-sga(size(sga,1)/2,1))/(2**(nm-(nm-5)+1)-1))
		write(*,"(es9.2$)")phyval
		write(*,"(1x)")
		!write(30,rec=i+1,fmt="(i4,21es13.5,A)")ne,var,(phyval(j),sqrt((sga(nm-5,j)-sga(size(sga,1)/2,j))/(2**(nm-(nm-5)+1)-1)),j=1,5),char(10)
		write(30,"(i4$)")ne
		write(30,"(es13.5$)")var,(phyval(j),sqrt((sga(nm-5,j)-sga(size(sga,1)/2,j))/(2**(nm-(nm-5)+1)-1)),j=1,5)
		write(30,"(1x)")
	enddo
end program
