module M_pmt
	use M_const
	implicit none
	integer, parameter :: Ns(2)=(/9,10/),Ns2=Ns(1)**2+1,Tx(2)=(/Ns(1),-1/),Ty(2)=(/1,Ns(1)/),r0(2)=(/Ns(1)/2,0/),vn=6
	integer :: neb(Ns2,4,3),ne=Ns2,ne2=Ns2/2
	real(8), parameter :: t(1)=(/1d0/),DJ=0.3d0,V=1d0
	real(8) :: ncfg(Ns2),scfg(Ns2),e2=0d0
end module
module M_wf
	use M_pmt
	use M_matrix
	use M_latt, only : &
         latt         =>  square_tilda         ,& 
         latt_bc      =>  square_tilda_bc      ,& 
         latt_one2two =>  square_tilda_one2two ,& 
         latt_two2one =>  square_tilda_two2one 
contains
	subroutine wf_k(var,wf,dwf)
		complex(8) :: wf(:,:),dwf(:,:,:),uv,duv(2),tmp(2),dtmp(2,2)
		real(8) :: var(vn),a(2),da(2,vn),gk,eks,eka,k(2),ek(2),E(2),pgk(2),sck,u2,v2,du2(2),dv2(2),etmp
		integer :: i,j,n,m,x,y,en,ik(2),ix(2),iy(2)
		wf=0d0
		dwf=0d0
		en=size(wf)/2
		!var(1)=abs(var(1))
		!var(3)=abs(var(3))
		!var(4)=abs(var(4))
		do n=-(Ns(1)-1)/2,(Ns(1)-1)/2
			do m=-(Ns(1)-1)/2,(Ns(1)+1)/2
				!call square_one2two(j,Ns,ik)
				!k=(/2d0*pi/Ns(1)*ik(1),2*pi/Ns(2)*ik(2)+pi/Ns(2)/)-pi
				k=2d0*pi*(/real(n*Ns(1)+m)/Ns2,real(-n+m*Ns(1))/Ns2/)
				if(abs(k(1))+abs(k(2))>pi) then
					cycle
				endif
				!if(any(abs(abs(k))>pi)) then
				!cycle
				!endif
				eks=-4d0*var(5)*cos(k(1))*cos(k(2))-2d0*var(6)*(cos(2d0*k(1))+cos(2d0*k(2)))-var(2)
				eka=-2d0*t(1)*(cos(k(1))+cos(k(2)))
				gk=(cos(k(1))-cos(k(2)))
				sck=gk*var(1)
				pgk(1)=var(3)
				pgk(2)=gk*var(4)
				etmp=sqrt(pgk(1)**2+pgk(2)**2+eka**2)
				ek=eks+(/1d0,-1d0/)*etmp
				do
					if(abs(abs(eks)-abs(etmp))<1d-10) then
						etmp=etmp+1d-10
					else
						exit
					endif
				enddo
				E=sqrt(ek**2+sck**2)
				do
					if(any(1d0/(E+ek)>1d10)) then
						if(abs(gk)<1d-10) then
							sck=sck+sign(1d-10,sck)
						else
							var(1)=var(1)+1d-10
							sck=gk*var(1)
						endif
						E=sqrt(ek**2+sck**2)
					else
						exit
					endif
				enddo
				a=sck/(ek+E)*(/-1d0,1d0/)
				da(:,1)=ek/E*gk/(E+ek)*(/-1d0,1d0/)
				da(:,2)=a/E
				da(:,3)=da(:,2)*pgk(1)/etmp*(/-1d0,1d0/)
				da(:,4)=da(:,2)*pgk(2)*gk/etmp*(/-1d0,1d0/)
				da(:,5)=da(:,2)*4d0*cos(k(1))*cos(k(2))
				da(:,6)=da(:,2)*2d0*(cos(2d0*k(1))+cos(2d0*k(2)))
				u2=0.5d0*(1d0+eka/etmp)
				v2=0.5d0*(1d0-eka/etmp)
				uv=0.5d0*(pgk(1)/etmp+pgk(2)/etmp*img)
				du2=-0.5d0*(/eka*pgk(1)/etmp**3,eka*pgk(2)*gk/etmp**3/)
				dv2=0.5d0*(/eka*pgk(1)/etmp**3,eka*pgk(2)*gk/etmp**3/)
				duv=0.5d0*(/(eka**2+pgk(2)**2-img*pgk(2)*pgk(1))/etmp**3,(eka**2+pgk(1)**2+img*pgk(2)*pgk(1))*img*gk/etmp**3/)
				if(any(isnan(a)).or.any(isnan(da))) then
					write(*,"(A)")"NAN in wf, quit!!"
					write(*,"(es10.2$)")1d0/(E+ek)
					write(*,"(x)")
					stop
				endif
				do i=1,Ns2
					call latt_one2two(i,Ns,ix)
					x=(-1)**mod(sum(ix),2)
					do j=1,Ns2
						call latt_one2two(j,Ns,iy)
						y=(-1)**mod(sum(iy),2)
						tmp=(/u2-v2*x*y,-v2+u2*x*y/)-uv*y+conjg(uv)*x
						dtmp(:,1)=(/du2(1)-dv2(1)*x*y,-dv2(1)+du2(1)*x*y/)-duv(1)*y+conjg(duv(1))*x
						dtmp(:,2)=(/du2(2)-dv2(2)*x*y,-dv2(2)+du2(2)*x*y/)-duv(2)*y+conjg(duv(2))*x
						wf(i,j)=wf(i,j)+exp(img*sum(k*(ix-iy)))*sum(tmp*a)
						dwf(i,j,:)=dwf(i,j,:)+exp(img*sum(k*(ix-iy)))*&
							(/&
							sum(tmp(:)*da(:,1)),&
							sum(tmp(:)*da(:,2)),&
							sum(tmp(:)*da(:,3)+dtmp(:,1)*a(:)),&
							sum(tmp(:)*da(:,4)+dtmp(:,2)*a(:)),&
							sum(tmp(:)*da(:,5)),&
							sum(tmp(:)*da(:,6))&
							/)
					enddo
				enddo
				!write(10,"(es10.2$)")k,a
				!write(10,"(1x)")
			enddo
		enddo
		dwf=dwf/Ns2
		wf=wf/Ns2
		!do i=1,Ns2
			!call latt_one2two(i,Ns,ix)
			!write(40,"(2i3,2es12.4)")ix,wf(1,i)
		!enddo
		!stop
	end subroutine
	subroutine wf_r(var,wf,dwf)
		use lapack95, only: heevd,heevr
		complex(8) :: H(Ns2*2,Ns2*2),wf(:,:),dwf(:,:,:),D(Ns2*2,Ns2*2)
		real(8) :: eg(Ns2*2),var(vn)
		integer :: i,j,n,ix(2)
		H=0d0
		wf=0d0
		dwf=0d0
		do i=1,Ns2
			call latt_one2two(i,Ns,ix)
			do j=1,4
				do n=2,3
					H(i,neb(i,j,n))=-var(n+3)
					H(i+Ns2,neb(i,j,n)+Ns2)=-H(i,neb(i,j,n))
				enddo
				H(i,neb(i,j,1))=-t(1)+var(4)*img*(-1)**mod(sum(ix),2)*(-1)**mod(j,2)
				H(i+Ns2,neb(i,j,1)+Ns2)=-H(i,neb(i,j,1))
				H(i,neb(i,j,1)+Ns2)=var(1)
				H(i+Ns2,neb(i,j,1))=var(1)
			enddo
			! on site
			H(i,i)=var(3)*(-1)**mod(sum(ix),2)-var(2)
			H(i+Ns2,i+Ns2)=var(3)*(-1)**mod(sum(ix),2)+var(2)
		enddo
		call heevd(H,eg,"V")
		D(1:Ns2,1:Ns2)=transpose(conjg(H(1:Ns2,Ns2+1:Ns2*2)))
		D(Ns2+1:2*Ns2,1:Ns2)=-transpose((H(Ns2+1:2*Ns2,1:Ns2)))
		D(1:Ns2,Ns2+1:Ns2*2)=transpose(conjg(H(Ns2+1:2*Ns2,Ns2+1:2*Ns2)))
		D(Ns2+1:2*Ns2,Ns2+1:Ns2*2)=transpose(H(1:Ns2,1:Ns2))
		call matrix_inv(D(1:Ns2,1:Ns2))
		call matrix_inv(D(Ns2+1:2*Ns2,1:Ns2))
		do i=1,Ns2
			do j=1,Ns2
				!do n=1,Ns2*2
					!wf(i,j)=wf(i,j)+D(n,Ns2-i+1)*D(n,Ns2+j)
				!enddo
				do n=1,Ns2
					wf(i,j)=wf(i,j)+D(i,n)*D(n,Ns2+j)+D(i+Ns2,n)*D(n+Ns2,Ns2+j)
				enddo
			enddo
		enddo
		if(any(isnan(real(wf)))) then
			write(*,*)"NAN"
			stop
		endif
		dwf=0d0
		!do i=1,Ns2
			!call latt_one2two(i,Ns,ix)
			!write(40,"(2i3,2es12.4)")ix,wf(2,i)
			!!write(40,"(2i3,2es12.4)")ix,real(dot_product(H(i,1:Ns2),H(i,1:Ns2))),1d0-real(dot_product(H(i+Ns2,1:Ns2),H(i+Ns2,1:Ns2)))
		!enddo
		!stop
	end subroutine
end module
module M_vmc
	use M_pmt
	use M_matrix
	use M_rd
	use M_utility
	use M_wf, iniwf => wf_k
	implicit none
contains
	subroutine mc(wf,dwf,Nmc,phyval,sga,O,S,g,flag)
		complex(8) :: pb,A(ne2,ne2),iA(ne2,ne2),dA(ne2,ne2,vn),vu(ne2,2),dvu(ne2,2,vn),wf(:,:),dwf(:,:,:),Y(2,2),&
			O(vn),S(vn,vn),g(vn),lO(vn),lS(vn,vn),lg(vn)
		real(8) :: rpb,sga(:,:)
		complex(8) :: phyval(:),lphy(size(phyval,1))
		integer :: cfg(Ns2),icfg(Ns2),i,j,ti,tj,l,sg,n,Nmc(:),cr(2),acp,bi,bj
		logical :: flag
		!write(*,"(A)")"start monte carlo"
		if(mod(Ns2,2)/=0) then
			write(*,*)"Number of site is not even, exit!!"
		endif
		bi=1
		bj=size(sga,1)/2+1
		sga=0d0
		n=0
		phyval=0d0
		O=0d0
		S=0d0
		g=0d0
		acp=0
		do i=1,Ns2
			cfg(i)=i
		enddo
		call fisher_yates_shuffle(cfg,Ns2)
		do i=1,Ns2
			icfg(cfg(i))=i
		enddo
		do i=1,ne2
			do j=1,ne2
				A(i,j)=wf(cfg(i),cfg(j+ne2))
				if(flag) then
					do ti=1,vn
						dA(i,j,ti)=dwf(cfg(i),cfg(j+ne2),ti)
					enddo
				endif
			enddo
		enddo
		iA=A
		call matrix_inv(iA)
		do
			n=n+1
			call rtwosite(i,j,sg)
			select case(sg)
			case(1)
				cr=(/i,-1/)
				do l=1,ne2
					vu(l,1)=wf(cfg(j),cfg(l+ne2))
					if(flag) then
						dvu(l,1,:)=dwf(cfg(j),cfg(l+ne2),:)
					endif
				enddo
			case(2)
				cr=(/-1,i-ne2/)
				do l=1,ne2
					vu(l,2)=wf(cfg(l),cfg(j))
					if(flag) then
						dvu(l,2,:)=dwf(cfg(l),cfg(j),:)
					endif
				enddo
			case(3)
				cr=(/i,j-ne2/)
				do l=1,ne2
					vu(l,1)=wf(cfg(j),cfg(l+ne2))
					vu(l,2)=wf(cfg(l),cfg(i))
					if(flag) then
						dvu(l,1,:)=dwf(cfg(j),cfg(l+ne2),:)
						dvu(l,2,:)=dwf(cfg(l),cfg(i),:)
					endif
				enddo
				vu(cr(2),1)=wf(cfg(j),cfg(i))
				vu(cr(1),2)=wf(cfg(j),cfg(i))
				if(flag) then
					dvu(cr(2),1,:)=dwf(cfg(j),cfg(i),:)
					dvu(cr(1),2,:)=dwf(cfg(j),cfg(i),:)
				endif
			end select
			call det(vu,cr,A,iA,Y,pb,sg)
			call random_number(rpb)
			if(rpb<abs(pb)**2) then
				acp=acp+1
				call swap(icfg(cfg(i)),icfg(cfg(j)))
				call swap(cfg(i),cfg(j))
				call update(vu,cr,Y,pb,A,iA,dvu,dA,sg)
				if(mod(n,500)==0) then
					!write(*,"(I9,E14.3$)")n,rpb
					iA=A
					call matrix_inv(iA)
					!write(*,"(E14.3)")rpb
				endif
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
				else
					call phymeasure(cfg,icfg,wf,A,iA,lphy(2:5))
				endif
				!call realcfg(cfg)
				call binning(bi,bj,sga,(/real(lphy(1)),real(lphy(2)),imag(lphy(3)),real(lphy(1)),real(lphy(1))/))
				phyval=phyval+lphy
				!e2=e2+lphy(1)**2
				S=S+lS
				g=g+lg
				O=O+lO
			endif
			if(n>=(Nmc(1)+Nmc(2)*Nmc(3))) then
				phyval=phyval/Nmc(3)
				!e2=e2/Nmc
				S=S/Nmc(3)
				g=g/Nmc(3)
				O=O/Nmc(3)
				exit
			endif
		enddo
		!write(*,"('accept/total number is ',I6,'/',I6)")acp,n
	end subroutine
	subroutine det(vu,cr,A,iA,Y,pb,sg)
		complex(8) :: vu(:,:),A(:,:),iA(:,:),Y(:,:),pb
		integer :: cr(:),sg
		if(sg==3) then
			call det_ratio_rowcol(vu,cr,A,iA,Y,pb)
		elseif(sg==1) then
			call det_ratio_row(vu(:,1),cr(1),A,iA,pb)
		else
			call det_ratio_col(vu(:,2),cr(2),A,iA,pb)
		endif
	end subroutine
	subroutine update(vu,cr,Y,pb,A,iA,dvu,dA,sg)
		complex(8) :: vu(:,:),A(:,:),iA(:,:),dvu(:,:,:),dA(:,:,:),Y(:,:),pb
		integer :: cr(:),sg
		if(sg==3) then
			call inv_update_rowcol(vu,cr,Y,A,iA)
			dA(cr(1),:,:)=dvu(:,1,:)
			dA(:,cr(2),:)=dvu(:,2,:)
		elseif(sg==1) then
			call inv_update_row(vu(:,1),cr(1),pb,A,iA)
			dA(cr(1),:,:)=dvu(:,1,:)
		else
			call inv_update_col(vu(:,2),cr(2),pb,A,iA)
			dA(:,cr(2),:)=dvu(:,2,:)
		endif
	end subroutine
	subroutine rtwosite(i,j,sg)
		integer :: i,j,sg
		call irandom(i,ne)
		call irandom(j,Ns2-ne2)
		if(j>ne2.and.i<=ne2) then
			sg=1
			j=j+ne2
		elseif(j>ne2.and.i>ne2) then
			sg=2
			j=j+ne2
		else
			sg=3
			if(i>ne2) then
				call swap(i,j)
			else
				j=j+ne2
			endif
		endif
	end subroutine
	subroutine realcfg(cfg)
		integer :: cfg(:),ix(2),i
		do i=1,ne
			call latt_one2two(cfg(i),Ns,ix)
			ncfg(cfg(i))=ncfg(cfg(i))+1d0
			!scfg(cfg(i))=scfg(cfg(i))+(-1d0)**((i-1)/ne2))
			scfg(cfg(i))=scfg(cfg(i))+(-1d0)**((i-1)/ne2)*((-1d0)**(mod(ix(1)+ix(2),2)))
		enddo
	end subroutine
	subroutine phymeasure(cfg,icfg,wf,A,iA,lod)
		complex(8) :: A(:,:),iA(:,:),vu(ne2,2),Y(2,2),pb,wf(:,:),a11
		complex(8) :: lod(:)
		integer :: cfg(:),icfg(:),i,n,m,j,l,inb,sg,cr(2),x1(2),x2(2),nr(2),dr(2),r(2)
		lod=0d0
		return
		do i=1,ne
			call latt_one2two(cfg(i),Ns,x1)
			lod(1)=lod(1)+sign(-1d0,i-ne2-0.1d0)*(0.5d0-mod(x1(1)+x1(2),2))
		enddo
		do i=ne+1,Ns2
			nr(1)=cfg(i)
			call latt_one2two(nr(1),Ns,x1)
			do inb=1,4
				j=icfg(neb(nr(1),inb,1))
				if(j<=ne) then
					if(j>ne2) then
						sg=2
						cr=(/-1,j-ne2/)
						do l=1,ne2
							vu(l,2)=wf(cfg(l),nr(1))
						enddo
					else
						sg=1
						cr=(/j,-1/)
						do l=1,ne2
							vu(l,1)=wf(nr(1),cfg(l+ne2))
						enddo
					endif
					call det(vu,cr,A,iA,Y,pb,sg)
					!if(imag(pb)>1d-10) then
					!write(*,*)imag(pb)
					!endif
					!lod(2)=lod(2)+imag(dconjg(pb))*(0.5d0-mod(inb,2))*(0.5d0-mod(x1(1)+x1(2),2))
					lod(2)=lod(2)+dconjg(pb)*(0.5d0-mod(inb,2))*(0.5d0-mod(x1(1)+x1(2),2))
				endif
			enddo
		enddo
		do i=ne+1,Ns2
			nr(1)=cfg(i)
			do inb=1,4,1
				if(icfg(neb(nr(1),inb,1))>ne) then
					nr(2)=neb(nr(1),inb,1)
					do l=1,ne2
						vu(l,1)=wf(nr(1),cfg(l+ne2))
						vu(l,2)=wf(cfg(l),nr(2))
					enddo
					a11=wf(nr(1),nr(2))
					pb=0d0
					do l=1,ne2
						do m=1,ne2
							pb=pb+vu(l,1)*iA(l,m)*vu(m,2)
						enddo
					enddo
					pb=a11-pb
					!lod(3)=lod(3)+real(dconjg(pb))*(0.5d0-mod(inb,2))
					lod(3)=lod(3)+dconjg(pb)*(0.5d0-mod(inb,2))
				endif
			enddo
		enddo
		do i=ne+1,Ns2
			nr(1)=cfg(i)
			do inb=1,4,1
				nr(2)=neb(nr(1),inb,1)
				if(icfg(nr(2))>ne) then
					do j=1,ne2
						r(1)=cfg(j)
						do n=1,4,1
							r(2)=neb(r(1),n,1)
							if(icfg(r(2))>ne2.and.icfg(r(2))<=ne) then
								!if(any(neb(nr(1),:,:)==r(1)).or.any(neb(nr(1),:,:)==r(2)).or.any(neb(nr(2),:,:)==r(1)).or.any(neb(nr(2),:,:)==r(2))) then
								!cycle
								!endif
								cr=(/icfg(r(1)),icfg(r(2))-ne2/)
								do l=1,ne2
									vu(l,1)=wf(nr(1),cfg(l+ne2))
									vu(l,2)=wf(cfg(l),nr(2))
								enddo
								vu(cr(2),1)=wf(nr(1),nr(2))
								vu(cr(1),2)=wf(nr(1),nr(2))
								call det(vu,cr,A,iA,Y,pb,3)
								!lod(4)=lod(4)+real(dconjg(pb))*(0.5d0-mod(inb+n,2))
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
	subroutine energy(cfg,icfg,wf,A,iA,El)
		complex(8) :: A(:,:),iA(:,:),vu(ne2,2),Y(2,2),pb,wf(:,:)
		complex(8) :: El
		integer :: cfg(:),icfg(:),i,ii,j,k,l,inb,sg,cr(2)
		El=0d0
		do i=1,ne
			do ii=1,4
				do inb=1,size(t)
					k=neb(cfg(i),ii,inb)
					j=icfg(k)
					if(j>ne) then
						!hopping
						if(i<=ne2) then
							sg=1
							cr=(/i,-1/)
							do l=1,ne2
								vu(l,1)=wf(cfg(j),cfg(l+ne2))
							enddo
						else
							sg=2
							cr=(/-1,i-ne2/)
							do l=1,ne2
								vu(l,2)=wf(cfg(l),cfg(j))
							enddo
						endif
						call det(vu,cr,A,iA,Y,pb,sg)
						!El=El-t(inb)*real(dconjg(pb))
						El=El-t(inb)*dconjg(pb)
						cycle
					endif
					if(inb>1) then
						cycle
					endif
					El=El+0.5d0*V
					!El=El+0.125d0*DJ
					if(j>ne2.and.i<=ne2) then
						!diagnal
						El=El-0.5d0*DJ
						!spin flip
						cr=(/i,j-ne2/)
						do l=1,ne2
							vu(l,1)=wf(cfg(j),cfg(l+ne2))
							vu(l,2)=wf(cfg(l),cfg(i))
						enddo
						vu(cr(2),1)=wf(cfg(j),cfg(i))
						vu(cr(1),2)=wf(cfg(j),cfg(i))
						call det(vu,cr,A,iA,Y,pb,3)
						!El=El-0.5d0*DJ*real(dconjg(pb))
						El=El-0.5d0*DJ*dconjg(pb)
					endif
				enddo
			enddo
		enddo
		El=El/Ns2
	end subroutine
	subroutine checkcfg(icfg,s)
		integer :: i,s,icfg(:)
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
		complex(8) :: wf(Ns2,Ns2),dwf(Ns2,Ns2,vn),O(vn),S(vn,vn),g(vn),Ov(vn),Sv(vn,vn),gv(vn)
		complex(8) :: tmp(vn,vn)
		integer :: n,i,j,k,info,sg(vn),cs,Nmc(:),nm,l
		real(8) :: var(vn),dvar(vn),pdvar(vn),eg(vn),dt=3d0,er(1),Ev,mEv
		complex(8) :: phyval(1)
		real(8), allocatable :: sga(:,:)
		logical :: flag
		nm=log(real(Nmc(3)))/log(2d0)
		allocate(sga((nm+1)*2,1))
		n=32
		pdvar=0d0
		cs=0
		sg=0
		mEv=100d0
		flag=.false.
		l=0
		do
			l=l+1
			write(*,"(I3$)")ne
			!if(var(1)<1d-3) then
			!var(1)=abs(var(1))
			!var(3:4)=abs(var(3:4))
			!endif
			write(*,"(es9.2$)")var
			call iniwf(var,wf,dwf)
			Ev=0d0
			Sv=0d0
			gv=0d0
			er=0d0
			!$OMP PARALLEL DO REDUCTION(+:Ev,Sv,gv,er) PRIVATE(phyval,O,S,g) LASTPRIVATE(sga) SCHEDULE(STATIC)
			do k=1,n
				call mc(wf,dwf,Nmc,phyval,sga,O,S,g,.true.)
				Ev=Ev+real(phyval(1))
				do i=1,vn
					do j=1,vn
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
			er=sqrt(abs(er/n-abs(Ev)**2))
			!Sv(:,4)=0d0
			!Sv(4,:)=0d0
			!Sv(:,2:5)=0d0
			!Sv(2:5,:)=0d0
			Sv=Sv+diag((/1d0,1d0,1d0,1d0,1d0,1d0/)*1d-2)
			call heevd(Sv,eg,"V")
			!write(*,*)sum(abs(matmul(matmul(Sv,diag(eg)),dconjg(transpose(Sv)))-tmp)),eg(vn)/eg(1)
			!if(eg(1)<=0d0) then
			!write(*,*)"S matrix is not positive define or singular"
			!endif
			!gv(4)=0d0
			!gv(2:5)=0d0
			!dvar=real(matmul(matmul(matmul(Sv,diag(1d0/eg)),dconjg(transpose(Sv))),gv))
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
			write(*,"(es9.2$)")real(gv),Ev,er
			write(20,"(I4$)")ne
			write(20,"(es13.5$)")var,real(gv),Ev,er
			write(20,"(x)")
			if(sum(abs(dvar))<vn*3d-4) then
				write(*,"(i3$)")cs
				cs=cs+1
				if(cs>=12.or.l>500) then
					!if(sum(abs(gv))<1d-3) then
					write(*,*)"variational finish",sum(abs(dvar))
					write(20,"(x)")
					write(10,"(I4$)")ne
					write(10,"(es13.5$)")var,Ev,er
					write(10,"(x)")
					exit
				endif
			endif
			do while((var(1)-dvar(1)*dt)<0.or.(var(3)-dvar(3)*dt)<0.or.(var(4)-dvar(4)*dt)<0)
				dt=dt*0.9
			enddo
			var=var-dvar*dt
			pdvar=dvar
			write(*,"(x)")
		enddo
	end subroutine
end module
program main
	use M_vmc
	use M_utility
	integer :: n,i,j,info,Nmc(3),Nvar(3),nm
	complex(8) :: wf(Ns2,Ns2),dwf(Ns2,Ns2,vn),O(vn),S(vn,vn),g(vn),tmp(vn,vn),z(vn,vn),&
		Ov(vn),Sv(vn,vn),gv(vn),Ev
	real(8) :: var(vn),dvar(vn),Eb,dp(25),sphy(5)
	complex(8) :: er(5),phyval(5)
	real(8), allocatable :: sga(:,:)
	logical :: f
	!open(10,file="../data/2d.dat")
	!f=openfile(10,"../data/err.dat")
	f=openfile(20,"../data/var.dat")
	f=openfile(30,"../data/phyvar.dat")
	f=openfile(40,"../data/2d.dat")
	!open(21,file="../data/phyvar_save.dat",status="old",action="read",access="direct",form="formatted",recl=245)
	Nmc(3)=1024
	Nvar(3)=128
	nm=log(real(Nmc(3)))/log(2d0)
	allocate(sga((nm+1)*2,size(phyval)))
	call init_random_seed()
	call latt(Ns,Tx,Ty,neb)
	n=32
	!var=(/1d-1,0d0,1d-1,1d-1,0d0,0d0/)
	do i=0,20,1
		!var=(/1d-2,0d0,1d-2,1d-2,0d0,0d0/)
		var=(/1d-1,0d0,1d-1,1d-1,0d0,0d0/)
		!var(3)=1d-2*10**(i/10d0)
		ne=Ns2-i*2
		!read(21,rec=i+1,fmt="(I4,5e16.8)")ne,var
		ne2=ne/2
		!Nmc(1:2)=(/5000*Ns2,5*Ns2/)*(i+1)
		!Nvar(1:2)=(/500*Ns2,5*Ns2/)*(i+1)
		Nmc(1:2)=(/5000*Ns2,5*Ns2/)
		Nvar(1:2)=(/500*Ns2,5*Ns2/)
		write(30,"(i4$)")ne
		call variational(var,Nvar)
		write(*,"(i4$)")ne
		write(*,"(es9.2$)")var
		write(30,"(es13.5$)")var
		!write(30,"(1x)")
		!write(*,"(1x)")
		!cycle
		call iniwf(var,wf,dwf)
		!sphy=0d0
		!er=0d0
		!ncfg=0d0
		!scfg=0d0
		!!$OMP PARALLEL DO REDUCTION(+:sphy,er) PRIVATE(phyval) LASTPRIVATE(sga) SCHEDULE(STATIC)
		!do j=1,n
		!call mc(wf,dwf,Nmc,phyval,sga,O,S,g,.false.)
		!sphy(1)=sphy(1)+phyval(1)
		!sphy(2:5)=sphy(2:5)+abs(phyval(2:5))
		!er=er+phyval**2
		!enddo
		!!$OMP END PARALLEL DO
		!sphy=sphy/n
		!!write(*,*)sqrt((e2-er(1))/(Nmc-1))
		!er=sqrt(abs(er/n-abs(sphy)**2))
		!write(*,"(es9.2$)")sphy
		!write(*,"(1x)")
		!write(30,"(e16.8$)")(sphy(j),er(j),j=1,5)
		!write(30,"(1x)")
		call mc(wf,dwf,Nmc,phyval,sga,O,S,g,.false.)
		do j=1,size(sga,1)/2-1
			write(10,"(es13.5$)")(sqrt((sga(j,k)-sga(size(sga,1)/2,k))/(2**(nm-j+1)-1)),sqrt((sga(j,k)-sga(size(sga,1)/2,k))/(2**(nm-j+1)-1)*sqrt(2d0/(2**(nm-j+1)-1))),k=1,size(sga,2))
			write(10,"(x)")
		enddo
		!write(*,"(es9.2$)")phyval(1),sqrt((sga(nm-5,1)-sga(size(sga,1)/2,1))/(2**(nm-(nm-5)+1)-1))
		write(*,"(es9.2$)")phyval
		write(*,"(1x)")
		write(30,"(es13.5$)")phyval(1),sqrt((sga(nm-5,1)-sga(size(sga,1)/2,1))/(2**(nm-(nm-5)+1)-1))
		write(30,"(es13.5$)")(phyval(j),sqrt((sga(nm-5,j)-sga(size(sga,1)/2,j))/(2**(nm-(nm-5)+1)-1)),j=2,5)
		write(30,"(1x)")
		do j=1,Ns(2)
			write(40,"(es13.5$)")ncfg((j-1)*Ns(1)+1:j*Ns(1))
			write(40,"(1x)")
		enddo
		write(40,"(1x)")
		do j=1,Ns(2)
			write(40,"(es13.5$)")scfg((j-1)*Ns(1)+1:j*Ns(1))
			write(40,"(1x)")
		enddo
		write(40,"(1x)")
	enddo
end program
