module M_pmt
	use M_const
	implicit none
	integer, parameter :: Ns(2)=(/9,10/),Ns2=Ns(1)**2+1,Tx(2)=(/Ns(1),-1/),Ty(2)=(/1,Ns(1)/),r0(2)=(/Ns(1)/2,0/),vn=5
	integer :: neb(Ns2,4,3),ne=Ns2,ne2=Ns2/2
	real(8), parameter :: t(1)=(/1d0/),DJ=0.3d0,V=-1d0
	!character(3) :: pgflag="sdw"
	character(3) :: pgflag="ddw"
	real(8) :: ncfg(Ns2),scfg(Ns2),e2=0d0
end module
module M_vmc
	use M_pmt
	use M_matrix
	use M_rd
	use M_utility
	use M_latt
	implicit none
contains
	subroutine iniwf(var,wf,dwf)
		complex(8) :: wf(:),dwf(:,:),uv,duv,tmpe(2),dtmpe(2),tmpo(2),dtmpo(2)
		real(8) :: var(vn),a(2),da(2,vn),gk,eks,eka,k(2),ek(2),E(2),pgk,sck,u2,v2,du2,dv2,tmp
		integer :: i,j,n,m,en,ik(2),ix(2),eo,r02(2)
		wf=0d0
		dwf=0d0
		en=size(wf)/2
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
				eks=-4d0*var(4)*cos(k(1))*cos(k(2))-2d0*var(5)*(cos(2d0*k(1))+cos(2d0*k(2)))-var(2)
				eka=-2d0*t(1)*(cos(k(1))+cos(k(2)))
				gk=(cos(k(1))-cos(k(2)))
				sck=gk*var(1)
				select case(pgflag)
				case("ddw")
					pgk=gk*var(3)
				case("sdw")
					pgk=var(3)
				end select
				tmp=sqrt(pgk**2+eka**2)
				ek=eks+(/1d0,-1d0/)*tmp
				E=sqrt(ek**2+sck**2)
				do
					if(any(1d0/(E+ek)>1d10)) then
						sck=sck+sign(1d-10,gk)
						E=sqrt(ek**2+sck**2)
					else
						exit
					endif
				enddo
				a=sck/(ek+E)*(/-1d0,1d0/)
				da(:,1)=ek*gk/(E*(E+ek))*(/-1d0,1d0/)
				da(:,2)=a/E
				da(:,4)=da(:,2)*4d0*cos(k(1))*cos(k(2))
				da(:,5)=da(:,2)*2d0*(cos(2d0*k(1))+cos(2d0*k(2)))
				select case(pgflag)
				case("ddw")
					!write(*,*)"ddw"
					da(:,3)=da(:,2)*pgk*gk/tmp*(/-1d0,1d0/)
					u2=0.5d0*(1d0+eka/tmp)
					v2=0.5d0*(1d0-eka/tmp)
					uv=0.5d0*pgk/tmp*img
					du2=-0.5d0*eka*pgk*gk/tmp**3
					dv2=0.5d0*eka*pgk*gk/tmp**3
					duv=0.5d0*(eka**2*gk/tmp**3)*img
				case("sdw")
					!write(*,*)"sdw"
					da(:,3)=da(:,2)*pgk/tmp*(/-1d0,1d0/)
					u2=0.5d0*(1d0+eka/tmp)
					v2=0.5d0*(1d0-eka/tmp)
					uv=0.5d0*pgk/tmp
					du2=-0.5d0*eka*pgk/tmp**3
					dv2=0.5d0*eka*pgk/tmp**3
					duv=0.5d0*(eka**2/tmp**3)
				end select
				if(any(isnan(a)).or.any(isnan(da))) then
					write(*,"(A)")"NAN in wf, quit!!"
					write(*,"(es10.2$)")1d0/(E+ek)
					write(*,"(x)")
					stop
				endif
				do i=1,Ns2
					call square_tilda_one2two(i,Ns,ix)
					ix=ix-r0
					eo=(-1)**mod(sum(ix),2)
					tmpe=(/u2-eo*v2,-v2+eo*u2/)
					dtmpe=(/du2-eo*dv2,-dv2+eo*du2/)
					!tmpo=(-eo*dconjg(uv)+uv)
					tmpo=(/-eo*uv+uv,-eo*conjg(uv)+conjg(uv)/)
					dtmpo=(/-eo*duv+duv,-eo*conjg(duv)+conjg(duv)/)
					wf(i)=wf(i)+exp(img*sum(k*ix))*sum((tmpe+tmpo)*a)
					wf(i+en)=wf(i+en)+exp(img*sum(k*ix))*sum((tmpe-tmpo)*a)
					dwf(i,:)=dwf(i,:)+exp(img*sum(k*ix))*&
						(/&
						(tmpe(1)+tmpo(1))*da(1,1)+(tmpe(2)+tmpo(2))*da(2,1)&
						,(tmpe(1)+tmpo(1))*da(1,2)+(tmpe(2)+tmpo(2))*da(2,2)&
						,(tmpe(1)+tmpo(1))*da(1,3)+(tmpe(2)+tmpo(2))*da(2,3)+(dtmpe(1)+dtmpo(1))*a(1)+(dtmpe(2)+dtmpo(2))*a(2)&
						,(tmpe(1)+tmpo(1))*da(1,4)+(tmpe(2)+tmpo(2))*da(2,4)&
						,(tmpe(1)+tmpo(1))*da(1,5)+(tmpe(2)+tmpo(2))*da(2,5)&
						/)
					dwf(i+en,:)=dwf(i+en,:)+exp(img*sum(k*ix))*&
						(/&
						(tmpe(1)-tmpo(1))*da(1,1)+(tmpe(2)-tmpo(2))*da(2,1)&
						,(tmpe(1)-tmpo(1))*da(1,2)+(tmpe(2)-tmpo(2))*da(2,2)&
						,(tmpe(1)-tmpo(1))*da(1,3)+(tmpe(2)-tmpo(2))*da(2,3)+(dtmpe(1)-dtmpo(1))*a(1)+(dtmpe(2)-dtmpo(2))*a(2)&
						,(tmpe(1)-tmpo(1))*da(1,4)+(tmpe(2)-tmpo(2))*da(2,4)&
						,(tmpe(1)-tmpo(1))*da(1,5)+(tmpe(2)-tmpo(2))*da(2,5)&
						/)
				enddo
				!write(10,"(es10.2$)")k,a
				!write(10,"(1x)")
			enddo
		enddo
		!stop
		dwf=dwf/Ns2
		wf=wf/Ns2
		!do i=1,Ns2
			!call square_tilda_one2two(i,Ns,ix)
			!write(10,"(es18.6$)")real(ix),real(wf(i))
			!write(10,"(x)")
		!enddo
		!stop
		!dwf=dwf/maxval(abs(wf))
		!wf=wf/maxval(abs(wf))
	end subroutine
	subroutine pair(ri,rj,wf,a)
		integer :: ri,rj,ri2(2),rj2(2),dr(2),n,r02(2)
		complex(8) :: wf(:),a
		call square_tilda_one2two(ri,Ns,ri2)
		call square_tilda_one2two(rj,Ns,rj2)
		dr=ri2-rj2+r0
		call square_tilda_bc(dr,Ns,Tx,Ty,Ns2,n)
		a=wf(n+size(wf)/2*mod(sum(ri2),2))
	end subroutine
	subroutine mc(wf,dwf,Nmc,phyval,sga,O,S,g,flag)
		complex(8) :: pb,A(ne2,ne2),iA(ne2,ne2),dA(ne2,ne2,vn),vu(ne2,2),dvu(ne2,2,vn),wf(:),dwf(:,:),Y(2,2),&
			O(vn),S(vn,vn),g(vn),lO(vn),lS(vn,vn),lg(vn)
		real(8) :: rpb,sga(:,:)
		complex(8) :: phyval(:),lphy(size(phyval,1))
		integer :: cfg(Ns2),icfg(Ns2),i,j,ti,tj,l,sg,n,Nmc,Nhot,cr(2),acp,bi,bj
		logical :: flag
		!write(*,"(A)")"start monte carlo"
		if(mod(Ns2,2)/=0) then
			write(*,*)"Number of site is not even, exit!!"
		endif
		n=0
		!Nhot=Nmc
		Nhot=min(50000,Nmc)
		bi=1
		bj=size(sga,1)/2+1
		sga=0d0
		do
			if(n==0) then
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
						call pair(cfg(i),cfg(j+ne2),wf,A(i,j))
						do ti=1,vn
							call pair(cfg(i),cfg(j+ne2),dwf(:,ti),dA(i,j,ti))
						enddo
					enddo
				enddo
				iA=A
				call matrix_inv(iA)
			endif
			n=n+1
			if(n==Nhot) then
				!call checkcfg(icfg,1000)
				iA=A
				call matrix_inv(iA)
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
			endif
			call irandom(i,ne)
			call irandom(j,Ns2-ne2)
			vu=0d0
			cr=0
			if(j>ne2.and.i<=ne2) then
				sg=1
				cr=(/i,-1/)
				j=j+ne2
				do l=1,ne2
					call pair(cfg(j),cfg(l+ne2),wf,vu(l,1))
					do ti=1,vn
						call pair(cfg(j),cfg(l+ne2),dwf(:,ti),dvu(l,1,ti))
					enddo
				enddo
			elseif(j>ne2.and.i>ne2) then
				sg=2
				cr=(/-1,i-ne2/)
				j=j+ne2
				do l=1,ne2
					call pair(cfg(l),cfg(j),wf,vu(l,2))
					do ti=1,vn
						call pair(cfg(l),cfg(j),dwf(:,ti),dvu(l,2,ti))
					enddo
				enddo
			else
				sg=3
				if(i>ne2) then
					call swap(i,j)
				else
					j=j+ne2
				endif
				cr=(/i,j-ne2/)
				do l=1,ne2
					call pair(cfg(j),cfg(l+ne2),wf,vu(l,1))
					call pair(cfg(l),cfg(i),wf,vu(l,2))
					do ti=1,vn
						call pair(cfg(j),cfg(l+ne2),dwf(:,ti),dvu(l,1,ti))
						call pair(cfg(l),cfg(i),dwf(:,ti),dvu(l,2,ti))
					enddo
				enddo
				call pair(cfg(j),cfg(i),wf,vu(cr(2),1))
				call pair(cfg(j),cfg(i),wf,vu(cr(1),2))
				do ti=1,vn
					call pair(cfg(j),cfg(i),dwf(:,ti),dvu(cr(2),1,ti))
					call pair(cfg(j),cfg(i),dwf(:,ti),dvu(cr(1),2,ti))
				enddo
			endif
			call det(vu,cr,A,iA,Y,pb,sg)
			!write(*,*)A
			!if(abs(pb)>1d10.and.n>1) then
			!call mwrite(10,A)
			!write(*,*)"initial configure may be not good",n
			!n=0
			!cycle
			!endif
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
				if(n>Nhot) then
					!call checkcfg(icfg,1000)
				!if(.false.) then
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
				endif
			endif
			if(n>Nhot) then
				call realcfg(cfg)
				call binning(bi,bj,sga,real(lphy))
				phyval=phyval+lphy
				!e2=e2+lphy(1)**2
				S=S+lS
				g=g+lg
				O=O+lO
			endif
			if(n>=(Nmc+Nhot)) then
				phyval=phyval/Nmc
				!e2=e2/Nmc
				S=S/Nmc
				g=g/Nmc
				O=O/Nmc
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
	subroutine realcfg(cfg)
		integer :: cfg(:),ix(2),i
		do i=1,ne
			call square_tilda_one2two(cfg(i),Ns,ix)
			ncfg(cfg(i))=ncfg(cfg(i))+1d0
			!scfg(cfg(i))=scfg(cfg(i))+(-1d0)**((i-1)/ne2))
			scfg(cfg(i))=scfg(cfg(i))+(-1d0)**((i-1)/ne2)*((-1d0)**(mod(ix(1)+ix(2),2)))
		enddo
	end subroutine
	subroutine phymeasure(cfg,icfg,wf,A,iA,lod)
		complex(8) :: A(:,:),iA(:,:),vu(ne2,2),Y(2,2),pb,wf(:),a11
		complex(8) :: lod(:)
		integer :: cfg(:),icfg(:),i,n,m,j,l,inb,sg,cr(2),x1(2),x2(2),nr(2),dr(2),r(2)
		lod=0d0
		do i=1,ne
			call square_tilda_one2two(cfg(i),Ns,x1)
			lod(1)=lod(1)+sign(-1d0,i-ne2-0.1d0)*(0.5d0-mod(x1(1)+x1(2),2))
		enddo
		do i=ne+1,Ns2
			nr(1)=cfg(i)
			call square_tilda_one2two(nr(1),Ns,x1)
			do inb=1,4
				j=icfg(neb(nr(1),inb,1))
				if(j<=ne) then
					if(j>ne2) then
						sg=2d0
						cr=(/-1,j-ne2/)
						do l=1,ne2
							call pair(cfg(l),nr(1),wf,vu(l,2))
						enddo
					else
						sg=1d0
						cr=(/j,-1/)
						do l=1,ne2
							call pair(nr(1),cfg(l+ne2),wf,vu(l,1))
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
						call pair(nr(1),cfg(l+ne2),wf,vu(l,1))
						call pair(cfg(l),nr(2),wf,vu(l,2))
					enddo
					call pair(nr(1),nr(2),wf,a11)
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
									call pair(nr(1),cfg(l+ne2),wf,vu(l,1))
									call pair(cfg(l),nr(2),wf,vu(l,2))
								enddo
								call pair(nr(1),nr(2),wf,vu(cr(2),1))
								call pair(nr(1),nr(2),wf,vu(cr(1),2))
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
		complex(8) :: A(:,:),iA(:,:),vu(ne2,2),Y(2,2),pb,wf(:)
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
								call pair(cfg(j),cfg(l+ne2),wf,vu(l,1))
							enddo
						else
							sg=2
							cr=(/-1,i-ne2/)
							do l=1,ne2
								call pair(cfg(l),cfg(j),wf,vu(l,2))
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
					if(j>ne2.and.i<=ne2) then
						!diagnal
						El=El-0.5d0*DJ
						!spin flip
						cr=(/i,j-ne2/)
						do l=1,ne2
							call pair(cfg(j),cfg(l+ne2),wf,vu(l,1))
							call pair(cfg(l),cfg(i),wf,vu(l,2))
						enddo
						call pair(cfg(j),cfg(i),wf,vu(cr(2),1))
						call pair(cfg(j),cfg(i),wf,vu(cr(1),2))
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
		integer :: i,s,icfg(:),r02(2)
		call square_tilda_one2two(1,Ns,r02)
		write(*,"(1X)")
		do i=1,r02(1)
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
			if(mod(i+r02(1),Ns(1))==0) then
				write(*,"(1X)")
			endif
		enddo
		call sleepqq(s)
		write(*,"(1X)")
	end subroutine
	subroutine variational(var,nm)
		use lapack95, only: heevd,heevr
		complex(8) :: wf(Ns2*2),dwf(Ns2*2,vn),O(vn),S(vn,vn),g(vn),Ov(vn),Sv(vn,vn),gv(vn)
		complex(8) :: tmp(vn,vn)
		integer :: n,i,j,k,info,sg(vn),cs,nm,Nmc
		real(8) :: var(vn),dvar(vn),pdvar(vn),eg(vn),dt=1d0,er(1),sga((nm+1)*2,1),Ev,mEv
		complex(8) :: phyval(1)
		open(20,file="../data/var.dat")
		Nmc=2**nm
		n=32
		pdvar=0d0
		cs=0
		sg=0
		mEv=100d0
		do
			write(*,"(I3$)")ne
			!if(var(1)<1d-3) then
				!var(1)=1d-3
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
			!Sv(:,2:5)=0d0
			!Sv(2:5,:)=0d0
			Sv=Sv+diag((/1d0,1d0,1d0,1d0,1d0/)*2d-1)
			call heevd(Sv,eg,"V")
			!write(*,*)sum(abs(matmul(matmul(Sv,diag(eg)),dconjg(transpose(Sv)))-tmp)),eg(vn)/eg(1)
			!if(eg(1)<=0d0) then
				!write(*,*)"S matrix is not positive define or singular"
			!endif
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
			write(*,"(x)")
			write(20,"(I4$)")ne
			write(20,"(e14.5$)")var,real(gv),Ev,er
			write(20,"(x)")
			if(real(Ev)<real(mEv)) then
				mEv=Ev
			endif
			!do i=1,vn
				!call find_cross(pdvar(i),dvar(i),sg(i))
			!enddo
			!if(all(sg/=0)) then
				cs=cs+1
			!endif
			!if(cs>5.or.sum(abs(dvar))<7d-4) then
			if(sum(abs(dvar))*dt<vn*3d-4.and.cs>40.and.abs(Ev-mEv)<abs(er(1)/2d0)) then
				!if(sum(abs(gv))<1d-3) then
				write(*,*)"variational finish",sum(abs(dvar))
				write(20,"(x)")
				write(10,"(I4$)")ne
				write(10,"(e14.5$)")var,Ev,er
				write(10,"(x)")
				exit
			endif
			var=var-dvar*dt
			pdvar=dvar
		enddo
	end subroutine
end module
program main
	use M_vmc
	integer :: n,i,j,info,nm=16,Nmc
	complex(8) :: wf(Ns2*2),dwf(Ns2*2,vn),O(vn),S(vn,vn),g(vn),tmp(vn,vn),z(vn,vn),&
		Ov(vn),Sv(vn,vn),gv(vn),Ev
	real(8) :: var(vn),dvar(vn),Eb,dp(25),sphy(5)
	complex(8) :: er(5),phyval(5)
	real(8), allocatable:: sga(:,:)
	!open(10,file="../data/2d.dat")
	open(10,file="../data/err.dat")
	open(30,file="../data/phyvar.dat")
	open(40,file="../data/2d.dat")
	open(21,file="../data/phyvar_save.dat",status="old",action="read",access="direct",form="formatted",recl=245)
	allocate(sga((nm+1)*2,size(phyval)))
	Nmc=2**nm
	call init_random_seed()
	call square_tilda(Ns,Tx,Ty,neb)
	n=32
	var=(/1d-1,0d0,1d-1,1d-1,1d-1/)
	do i=20,20,2
		!var(3)=1d-2*10**(i/10d0)
		ne=Ns2-i*2
		!read(21,rec=i+1,fmt="(I4,5e16.8)")ne,var
		ne2=ne/2
		write(30,"(i4$)")ne
		call variational(var,13)
		write(*,"(i4$)")ne
		write(*,"(es9.2$)")var
		write(30,"(e16.8$)")var
		!write(30,"(1x)")
		!write(*,"(1x)")
		!stop
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
			write(10,"(e16.8$)")(sqrt((sga(j,k)-sga(size(sga,1)/2,k))/(2**(nm-j+1)-1)),sqrt((sga(j,k)-sga(size(sga,1)/2,k))/(2**(nm-j+1)-1)*sqrt(2d0/(2**(nm-j+1)-1))),k=1,size(sga,2))
			write(10,"(x)")
		enddo
		!write(*,"(es9.2$)")phyval(1),sqrt((sga(nm-5,1)-sga(size(sga,1)/2,1))/(2**(nm-(nm-5)+1)-1))
		write(*,"(es9.2$)")phyval
		write(*,"(1x)")
		write(30,"(e16.8$)")phyval(1),sqrt((sga(nm-5,1)-sga(size(sga,1)/2,1))/(2**(nm-(nm-5)+1)-1))
		write(30,"(e16.8$)")(phyval(j),sqrt((sga(nm-5,j)-sga(size(sga,1)/2,j))/(2**(nm-(nm-5)+1)-1)),j=2,5)
		write(30,"(1x)")
		do j=1,Ns(2)
			write(40,"(e16.8$)")ncfg((j-1)*Ns(1)+1:j*Ns(1))
			write(40,"(1x)")
		enddo
		write(40,"(1x)")
		do j=1,Ns(2)
			write(40,"(e16.8$)")scfg((j-1)*Ns(1)+1:j*Ns(1))
			write(40,"(1x)")
		enddo
		write(40,"(1x)")
		!stop
	enddo
end program
