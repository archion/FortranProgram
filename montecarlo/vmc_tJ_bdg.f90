module M_pmt
	use M_const
	implicit none
	integer, parameter :: Ns(2)=(/9,10/),Ns2=Ns(1)**2+1,Tx(2)=(/Ns(1),-1/),Ty(2)=(/1,Ns(1)/),vn=6
	!integer, parameter :: Ns(2)=(/10,10/),Ns2=Ns(1)**2,Tx(2)=(/Ns(1),0/),Ty(2)=(/0,Ns(1)/),vn=6
	integer :: neb(Ns2,4,3),ne=Ns2,ne2=Ns2/2
	real(8), parameter :: t(1)=(/1d0/),DJ=1/3d0,V=4d0/3d0
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
			write(*,*)"NAN"
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
module M_vmc
	use M_pmt
	use M_matrix
	use M_rd
	use M_utility
	use M_wf, iniwf => wf_r
	implicit none
contains
	subroutine mc(wf,Nmc,phyval,sga,dwf,O,S,g)
		complex(8) :: pb,A(Ns2,Ns2),iA(Ns2,Ns2),dA(Ns2,Ns2,vn),vu(Ns2,2),dvu(Ns2,2,vn),wf(:,:),Y(2,2)
		complex(8), optional :: O(vn),S(vn,vn),g(vn),dwf(:,:,:)
		complex(8) :: lO(vn),lS(vn,vn),lg(vn)
		real(8) :: rpb,sga(:,:)
		complex(8) :: phyval(:),lphy(size(phyval,1))
		integer :: cfg(Ns2),ncfg(Ns2),icfg(Ns2),i,j,ti,tj,l,sg,n,Nmc(:),cr(2),acp,bi,bj
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
		do i=1,Ns2
			call getmat(wf,i,cfg,A(i,:),flag,dwf,dA(i,:,:))
		enddo
		iA=A
		call mat_inv(iA)
		do
			n=n+1
			call two(i,j)
			call getdiff(wf,i,j,cfg,sg,vu,cr,flag,dwf,dvu)
			call det(vu,cr,A,iA,Y,pb,sg)
			call random_number(rpb)
			if(rpb<real(pb*conjg(pb))) then
				acp=acp+1
				call swap(cfg(i),cfg(j))
				call swap(icfg(cfg(i)),icfg(cfg(j)))
				call update(vu,cr,Y,pb,A,iA,sg,flag,dvu,dA)
				!write(*,"(es9.2)")sum(abs(matmul(A,iA)-diag(1d0,size(A,1))))
				if(mod(n,500)==0) then
					iA=A
					call mat_inv(iA)
				endif
				!if(n>Nmc(1)) then
					!call checkcfg(icfg,100)
				!endif
			endif
			if(n>Nmc(1).and.mod(n-Nmc(1),Nmc(2))==0) then
				call energy(cfg,icfg,wf,A,iA,lphy(1))
				if(flag) then
					do ti=1,vn
						lO(ti)=Tr(iA,dA(:,:,ti))
					enddo
					do ti=1,vn
						do tj=1,vn
							lS(ti,tj)=dconjg(lO(ti))*lO(tj)
						enddo
					enddo
					lg=real(lphy(1)*lO)
					S=S+lS
					g=g+lg
					O=O+lO
				else
					call phymeasure(cfg,icfg,wf,A,iA,lphy(2:5))
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
					do ti=1,size(O)
						do tj=1,size(O)
							S(ti,tj)=S(ti,tj)-dconjg(O(ti))*O(tj)
						enddo
					enddo
					g=2d0*(g-real(phyval(1)*O))*Ns2
				endif
				exit
			endif
		enddo
		!write(*,"('accept/total number is ',I6,'/',I6)")acp,n
	end subroutine
	subroutine two(i,j)
		integer :: i,j,sg,cg(2,2),n
		call random_number(n,ne2**2+2*ne2*(Ns2-ne))
		if(n<=ne2*(Ns2-ne2)) then
			i=(n-1)/(Ns2-ne2)+1
			j=ne2+mod(n-1,(Ns2-ne2))+1
		else
			n=n-ne2*(Ns2-ne2)
			i=ne2+(n-1)/(Ns2-ne)+1
			j=ne+mod(n-1,(Ns2-ne))+1
		endif
	end subroutine
	subroutine getmat(wf,i,cfg,u,flag,dwf,du)
		complex(8) :: wf(:,:),u(:)
		complex(8), optional :: dwf(:,:,:),du(:,:)
		integer :: i,cfg(:)
		logical, optional :: flag
		if(i<=ne2) then
			u=wf(cfg(i),:)
			if(present(flag).and.flag) then
				du(:,:)=dwf(cfg(i),:,:)
			endif
		elseif(i>ne) then
			u(:)=wf(cfg(i)+Ns2,:)
			if(present(flag).and.flag) then
				du(:,:)=dwf(cfg(i)+Ns2,:,:)
			endif
		else
			u(:)=wf(cfg(i-ne2)+Ns2,:)
			if(present(flag).and.flag) then
				du(:,:)=dwf(cfg(i-ne2)+Ns2,:,:)
			endif
		endif
	end subroutine
	subroutine getdiff(wf,i,j,cfg,sg,vu,cr,flag,dwf,dvu)
		integer :: i,j,sg,cfg(:),ncfg(size(cfg,1)),cr(:)
		complex(8) :: vu(:,:),wf(:,:)
		complex(8), optional :: dvu(:,:,:),dwf(:,:,:)
		logical, optional :: flag
		ncfg=cfg
		if(j>ne.and.i<=ne2) then
			sg=1
			call swap(ncfg(i),ncfg(j))
			cr=(/i,j/)
			if(present(flag).and.flag) then
				call getmat(wf,cr(1),ncfg,vu(:,1),flag,dwf,dvu(:,1,:))
			else
				call getmat(wf,cr(1),ncfg,vu(:,1))
			endif
		elseif(j>ne.and.i<=ne) then
			sg=2
			call swap(ncfg(i),ncfg(j))
			cr=(/-1,j/)
			if(present(flag).and.flag) then
				call getmat(wf,cr(2),ncfg,vu(:,2),flag,dwf,dvu(:,2,:))
			else
				call getmat(wf,cr(2),ncfg,vu(:,2))
			endif
		elseif(j>ne2.and.i<=ne2) then
			sg=3
			call swap(ncfg(i),ncfg(j))
			cr=(/i,i+ne2/)
			if(present(flag).and.flag) then
				call getmat(wf,cr(1),ncfg,vu(:,1),flag,dwf,dvu(:,1,:))
				call getmat(wf,cr(2),ncfg,vu(:,2),flag,dwf,dvu(:,2,:))
			else
				call getmat(wf,cr(1),ncfg,vu(:,1))
				call getmat(wf,cr(2),ncfg,vu(:,2))
			endif
		else
			sg=-1
		endif
	end subroutine
	subroutine det(vu,cr,A,iA,Y,pb,sg)
		complex(8) :: vu(:,:),A(:,:),iA(:,:),Y(:,:),pb
		integer :: cr(:),sg
		select case(sg) 
		case(1) 
			call det_ratio_row(vu(:,1),cr(1),A,iA,pb)
		case(2) 
			call det_ratio_row(vu(:,2),cr(2),A,iA,pb)
		case(3)
			call det_ratio_tworow(vu,cr,A,iA,Y,pb)
		end select
	end subroutine
	subroutine update(vu,cr,Y,pb,A,iA,sg,flag,dvu,dA)
		complex(8) :: vu(:,:),A(:,:),iA(:,:),Y(:,:),pb
		complex(8) :: dvu(:,:,:),dA(:,:,:)
		integer :: cr(:),sg,ti
		logical :: flag
		select case(sg) 
		case(1) 
			call inv_update_row(vu(:,1),cr(1),pb,A,iA)
			call swap(A(cr(1)+ne2,:),A(cr(2),:))
			call swap(iA(:,cr(1)+ne2),iA(:,cr(2)))
			if(flag) then
				dA(cr(1),:,:)=dvu(:,1,:)
				do ti=1,vn
					call swap(dA(cr(1)+ne2,:,ti),dA(cr(2),:,ti))
				enddo
			endif
		case(2) 
			call inv_update_row(vu(:,2),cr(2),pb,A,iA)
			if(flag) then
				dA(cr(2),:,:)=dvu(:,2,:)
			endif
		case(3)
			call inv_update_tworow(vu,cr,Y,A,iA)
			if(flag) then
				dA(cr(1),:,:)=dvu(:,1,:)
				dA(cr(2),:,:)=dvu(:,2,:)
			endif
		end select
	end subroutine
	subroutine energy(cfg,icfg,wf,A,iA,El)
		complex(8) :: A(:,:),iA(:,:),vu(size(A,1),2),Y(2,2),pb,wf(:,:)
		complex(8) :: El
		integer :: cfg(:),icfg(:),i,ii,j,l,inb,sg,cr(2)
		El=0d0
		do i=1,ne
			do ii=1,4
				do inb=1,size(t)
					j=icfg(neb(cfg(i),ii,inb))
					call getdiff(wf,i,j,cfg,sg,vu,cr)
					select case(sg) 
					case(1) 
						call det(vu,cr,A,iA,Y,pb,sg)
						El=El-t(inb)*conjg(pb)
					case(2) 
						call det(vu,cr,A,iA,Y,pb,sg)
						El=El+t(inb)*conjg(pb)
					case default
						if(inb==1) then
							if(i<=ne2.and.j>ne2) then
								call det(vu,cr,A,iA,Y,pb,sg)
								El=El+0.5d0*DJ*conjg(pb)-0.25d0*DJ&
									-0.25d0*DJ&
									+V
							elseif((i<=ne2.or.j>ne2).and.ii<3) then
								El=El+0.25d0*DJ&
									-0.25d0*DJ&
									+V
							endif
						endif
					end select
				enddo
			enddo
		enddo
		El=El/Ns2
	end subroutine
	subroutine phymeasure(cfg,icfg,wf,A,iA,lod)
		complex(8) :: A(:,:),iA(:,:),vu(Ns2,2),Y(2,2),pb,wf(:,:),a11
		complex(8) :: lod(:)
		integer :: cfg(:),icfg(:),i,n,m,j,k,l,inb,sg,cr(2),x1(2),x2(2),nr(2),dr(2),r(2)
		lod=0d0
		!SDW Order
		do i=1,ne
			call latt_one2two(cfg(i),Ns,x1)
			lod(1)=lod(1)+sign(-1d0,i-ne2-0.1d0)*(0.5d0-mod(x1(1)+x1(2),2))
		enddo
		!DDW Order
		do j=ne+1,Ns2
			call latt_one2two(cfg(j),Ns,x1)
			do inb=1,4
				i=icfg(neb(cfg(j),inb,1))
				call getdiff(wf,i,j,cfg,sg,vu,cr)
				call det(vu,cr,A,iA,Y,pb,sg)
				select case(sg) 
				case(1) 
					lod(2)=lod(2)+conjg(pb)*(0.5d0-mod(inb,2))*(0.5d0-mod(sum(x1),2))
				case(2) 
					lod(2)=lod(2)-conjg(pb)*(0.5d0-mod(inb,2))*(0.5d0-mod(sum(x1),2))
				end select
			enddo
		enddo
		!d-SC Order
		!do i=ne+1,Ns2
			!do inb=1,4,1
				!j=icfg(neb(cfg(i),inb,1))
				!if(j>ne) then
					!vu(:,1)=wf(cfg(i),:)
					!cr=(/j,-1/)
					!call det(vu,cr,A,iA,Y,pb,1)
					!lod(3)=lod(3)+conjg(pb)*(0.5d0-mod(inb,2))
				!endif
			!enddo
		!enddo
		!d-SC corelation
		do i=ne+1,Ns2
			do inb=1,4,1
				j=icfg(neb(cfg(i),inb,1))
				if(j>ne) then
					do k=1,ne2
						do n=1,4,1
							l=icfg(neb(cfg(k),n,1))
							if(l>ne2.and.l<=ne) then
								vu(:,1)=wf(cfg(i),:)
								vu(:,2)=wf(Ns2+cfg(l),:)
								cr=(/k,j/)
								call det(vu,cr,A,iA,Y,pb,3)
								lod(4)=lod(4)+dconjg(pb)*(0.5d0-mod(inb+n,2))
							endif
						enddo
					enddo
				endif
			enddo
		enddo
		lod(4)=lod(4)/Ns2
		lod=lod/Ns2
	end subroutine
	subroutine checkcfg(icfg,s)
		integer :: i,s,icfg(:),r0(2)
		call latt_one2two(1,Ns,r0)
		write(*,"(1X)")
		do i=1,r0(1)
			write(*,"('  '$)")
		enddo
		do i=1,Ns2
			if(icfg(i)<=ne2) then
				write(*,"(' ●'$)")
			elseif(icfg(i)<=ne) then
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
		complex(8) :: wf(Ns2*2,Ns2),dwf(Ns2*2,Ns2,vn),O(vn),S(vn,vn),g(vn),Ov(vn),Sv(vn,vn),gv(vn)
		complex(8) :: tmp(vn,vn)
		integer :: n,i,j,k,info,sg(vn),cs,Nmc(:),nm,l
		real(8) :: var(vn),dvar(vn),pvar(vn),eg(vn),dt=0.03d0,er(1),Ev,allE(400),scv
		complex(8) :: phyval(1)
		real(8), allocatable :: sga(:,:)
		logical :: flag
		nm=log(real(Nmc(3)))/log(2d0)
		allocate(sga((nm+1)*2,1))
		n=24
		cs=0
		sg=0
		l=0
		scv=0d0
		flag=.true.
		pvar=var
		do
			l=l+1
			write(*,"(I3$)")ne
			write(*,"(es9.2$)")var
			call iniwf(var,wf,dwf)
			Ev=0d0
			Sv=0d0
			gv=0d0
			er=0d0
			!$OMP PARALLEL DO REDUCTION(+:Ev,Sv,gv,er) PRIVATE(phyval,O,S,g) LASTPRIVATE(sga)
			do k=1,n
				call mc(wf,Nmc,phyval,sga,dwf,O,S,g)
				Ev=Ev+real(phyval(1))
				Sv=Sv+S
				gv=gv+g
				er=er+abs(phyval)**2
			enddo
			!$OMP END PARALLEL DO
			Ev=Ev/n
			Sv=Sv/n
			gv=gv/n
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
			do i=1,vn
				dvar(i)=0d0
				do j=1,vn
					do k=1,vn
						if(eg(k)/eg(vn)<1d-3) then
							cycle
						endif
						dvar(i)=dvar(i)+real(Sv(i,k)*1d0/eg(k)*dconjg(Sv(j,k))*gv(j))
					enddo
				enddo
			enddo
			if(n<5) then
				er=sqrt((sga(nm-3,1)-sga(size(sga,1)/2,1))/(2**(nm-(nm-3)+1)-1))
			endif
			write(*,"(es9.2$)")real(gv)
			write(*,"(es11.4$)")Ev
			write(*,"(es9.2$)")er
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
	real(8) :: var(vn)
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
	do i=4,20,1
		var=(/1d-1,0d0,0d0,0d0,0d0,0d0/)
		var=(/6.72658E-02,-4.57858E-02,7.89528E-02,1.99349E-01,2.38169E-02,1.50484E-01/)
		!var(1)=1d-1
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
		call iniwf(var,wf)
		call mc(wf,Nmc,phyval,sga)
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
