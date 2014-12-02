module M_pmt
	use M_const
	implicit none
	integer, parameter :: Ns(2)=(/10,10/),Ns2=Ns(1)*Ns(2),vn=5
	integer :: neb(Ns2,4,3),ne=100,ne2
	real(8), parameter :: t(1)=(/1d0/),DJ=0.5d0,V=-2d0
	character(3) :: pgflag="sdw"
	!character(3) :: pgflag="ddw"
	real(8) :: ncfg(Ns2),scfg(Ns2)
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
		complex(8) :: wf(0:,0:),dwf(0:,0:,:),uv,duv,tmpe(2),dtmpe(2),tmpo(2),dtmpo(2)
		real(8) :: var(vn),a(2),da(2,vn),al,np,pnp,gk,eks,eka,k(2),ek(2),E(2),pgk,sck,u2,v2,du2,dv2,tmp
		integer :: i,j,n,ik(2),ix(2),eo
		logical :: flag
		al=0.01d0
		pnp=0
		wf=0d0
		dwf=0d0
		n=size(wf,2)/2
		do j=1,Ns2
			call square_one2two(j,Ns,ik)
			k=(/2d0*pi/Ns(1)*ik(1),2*pi/Ns(2)*ik(2)+pi/Ns(2)/)-pi
			eks=-4d0*var(4)*cos(k(1))*cos(k(2))-2d0*var(5)*(cos(2d0*k(1))+cos(2d0*k(2)))-var(2)
			eka=-2d0*t(1)*(cos(k(1))+cos(k(2)))
			gk=(cos(k(1))-cos(k(2)))
			if(abs(k(1))+abs(k(2))>pi) then
				cycle
			endif
			if(abs(var(1))<1d-3) then
				var(1)=sign(1d-3,var(1))
			endif
			sck=gk*var(1)
			select case(pgflag)
			case("ddw")
				pgk=gk*var(3)
			case("sdw")
				pgk=var(3)
			end select
			tmp=sqrt(pgk**2+eka**2)
			ek=eks+(/1d0,-1d0/)*tmp
			!if(abs(sck)<1d-6) then
			!sck=sign(1d-6,sck)
			!endif
			E=sqrt(ek**2+sck**2)
			a=sck/(ek+E)*(/-1d0,1d0/)
			a=min(abs(a),1d6)*sign(1d0,a)
			!if(sck<1d-3.and.ek<=0d0) then
				!a=min(2d0*abs(ek)/sck,1d8)*(/-1d0,1d0/)
			!else
			!endif
			!if(any(isnan(a))) then
				!write(*,"(A)")"NAN in wf, quit!!"
				!stop
			!endif
			!if(abs(ek)<1d-5) then
				!write(*,"(A)")"ek is too small"
				!stop
			!endif
			da(:,1)=ek*gk/(E**2+ek*E)*(/-1d0,1d0/)
			!da(:,2)=sck/(E**2+ek*E)*(/-1d0,1d0/)
			da(:,2)=a/E
			da(:,4)=da(:,2)*4d0*cos(k(1))*cos(k(2))
			da(:,5)=da(:,2)*2d0*(cos(2d0*k(1))+cos(2d0*k(2)))
			!if(sck<1d-6) then
			!a=(1d0-sign(1d0,ek))*1d10*(/-1d0,1d0/)
			!da(:,1)=((1d0-sign(1d0,ek))*1d10+1d0/abs(2d0*ek))*sign(1d0,ek)*(/-1d0,1d0/)
			!da(:,2)=(1d0-sign(1d0,ek))*1d10*(/1d0,-1d0/)
			!endif
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
			do i=1,Ns2
				call square_one2two(i,Ns,ix)
				eo=(-1)**mod(ix(1)+ix(2),2)
				tmpe=(/u2-eo*v2,-v2+eo*u2/)
				dtmpe=(/du2-eo*dv2,-dv2+eo*du2/)
				tmpo=(-eo*dconjg(uv)+uv)
				dtmpo=(-eo*dconjg(duv)+duv)
				wf(ix(1),ix(2))=wf(ix(1),ix(2))+exp(img*sum(k*ix))*sum((tmpe+tmpo)*a)
				!((tmpe(1)+tmpo(1))*a(1)+(tmpe(2)+tmpo(2))*a(2))
				wf(ix(1),ix(2)+n)=wf(ix(1),ix(2)+n)+exp(img*sum(k*ix))*sum((tmpe-tmpo)*a)
				!((tmpe(1)-tmpo(1))*a(1)+(tmpe(2)-tmpo(2))*a(2))
				dwf(ix(1),ix(2),:)=dwf(ix(1),ix(2),:)+exp(img*sum(k*ix))*&
					(/&
					(tmpe(1)+tmpo(1))*da(1,1)+(tmpe(2)+tmpo(2))*da(2,1)&
					,(tmpe(1)+tmpo(1))*da(1,2)+(tmpe(2)+tmpo(2))*da(2,2)&
					,(tmpe(1)+tmpo(1))*da(1,3)+(tmpe(2)+tmpo(2))*da(2,3)+(dtmpe(1)+dtmpo(1))*a(1)+(dtmpe(2)+dtmpo(2))*a(2)&
					,(tmpe(1)+tmpo(1))*da(1,4)+(tmpe(2)+tmpo(2))*da(2,4)&
					,(tmpe(1)+tmpo(1))*da(1,5)+(tmpe(2)+tmpo(2))*da(2,5)&
					/)
				dwf(ix(1),ix(2)+n,:)=dwf(ix(1),ix(2)+n,:)+exp(img*sum(k*ix))*&
					(/&
					(tmpe(1)-tmpo(1))*da(1,1)+(tmpe(2)-tmpo(2))*da(2,1)&
					,(tmpe(1)-tmpo(1))*da(1,2)+(tmpe(2)-tmpo(2))*da(2,2)&
					,(tmpe(1)-tmpo(1))*da(1,3)+(tmpe(2)-tmpo(2))*da(2,3)+(dtmpe(1)-dtmpo(1))*a(1)+(dtmpe(2)-dtmpo(2))*a(2)&
					,(tmpe(1)-tmpo(1))*da(1,4)+(tmpe(2)-tmpo(2))*da(2,4)&
					,(tmpe(1)-tmpo(1))*da(1,5)+(tmpe(2)-tmpo(2))*da(2,5)&
					/)
			enddo
		enddo
		dwf=dwf/Ns2
		wf=wf/Ns2
		!dwf=dwf/maxval(abs(wf))
		!wf=wf/maxval(abs(wf))
	end subroutine
	subroutine pair(ri,rj,wf,a)
		integer :: ri,rj,ri2(2),rj2(2),dr(2)
		complex(8) :: wf(0:,0:),a
		real(8) :: s
		call square_one2two(ri,Ns,ri2)
		call square_one2two(rj,Ns,rj2)
		dr=ri2-rj2
		s=sign(1d0,dr(2)+0.1d0)
		dr=mod(dr+Ns,Ns)
		a=wf(dr(1),dr(2)+size(wf,2)/2*mod(ri2(1)+ri2(2),2))*s
	end subroutine
	subroutine mc(wf,dwf,Nmc,phyval,O,S,g,flag)
		complex(8) :: pb,A(ne2,ne2),iA(ne2,ne2),dA(ne2,ne2,vn),vu(ne2,2),dvu(ne2,2,vn),wf(:,:),dwf(:,:,:),Y(2,2),&
			phyval(:),lphy(size(phyval,1)),O(vn),S(vn,vn),g(vn),lO(vn),lS(vn,vn),lg(vn)
		real(8) :: rpb
		integer :: cfg(Ns2),icfg(Ns2),i,j,ti,tj,l,sg,n,Nmc,Nhot,cr(2),acp
		logical :: flag
		!write(*,"(A)")"start monte carlo"
		if(mod(Ns2,2)/=0) then
			write(*,*)"Number of site is not even, exit!!"
		endif
		n=0
		!Nhot=Nmc
		Nhot=min(10000,Nmc)
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
							call pair(cfg(i),cfg(j+ne2),dwf(:,:,ti),dA(i,j,ti))
						enddo
					enddo
				enddo
				iA=A
				call matrix_inv(iA)
			endif
			n=n+1
			if(n==Nhot) then
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
					call phymeasure(cfg,icfg,wf,A,iA,lphy(2),lphy(3))
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
						call pair(cfg(j),cfg(l+ne2),dwf(:,:,ti),dvu(l,1,ti))
					enddo
				enddo
			elseif(j>ne2.and.i>ne2) then
				sg=2
				cr=(/-1,i-ne2/)
				j=j+ne2
				do l=1,ne2
					call pair(cfg(l),cfg(j),wf,vu(l,2))
					do ti=1,vn
						call pair(cfg(l),cfg(j),dwf(:,:,ti),dvu(l,2,ti))
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
						call pair(cfg(j),cfg(l+ne2),dwf(:,:,ti),dvu(l,1,ti))
						call pair(cfg(l),cfg(i),dwf(:,:,ti),dvu(l,2,ti))
					enddo
				enddo
				call pair(cfg(j),cfg(i),wf,vu(cr(2),1))
				call pair(cfg(j),cfg(i),wf,vu(cr(1),2))
				do ti=1,vn
					call pair(cfg(j),cfg(i),dwf(:,:,ti),dvu(cr(2),1,ti))
					call pair(cfg(j),cfg(i),dwf(:,:,ti),dvu(cr(1),2,ti))
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
				!call checkcfg(icfg,80)
				if(mod(n,500)==0) then
					!write(*,"(I9,E14.3$)")n,rpb
					iA=A
					call matrix_inv(iA)
					!write(*,"(E14.3)")rpb
				endif
				if(n>Nhot) then
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
						call phymeasure(cfg,icfg,wf,A,iA,lphy(2),lphy(3))
					endif
				endif
			endif
			if(n>Nhot) then
				call realcfg(cfg)
				phyval=phyval+lphy
				S=S+lS
				g=g+lg
				O=O+lO
			endif
			if(n>=(Nmc+Nhot)) then
				phyval=phyval/(Nmc*Ns2)
				S=S/(Nmc)
				g=g/(Nmc*Ns2)
				O=O/(Nmc)
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
			call square_one2two(cfg(i),Ns,ix)
			ncfg(cfg(i))=ncfg(cfg(i))+1d0
			!scfg(cfg(i))=scfg(cfg(i))+(-1d0)**((i-1)/ne2))
			scfg(cfg(i))=scfg(cfg(i))+(-1d0)**((i-1)/ne2)*((-1d0)**(mod(ix(1)+ix(2),2)))
		enddo
	end subroutine
	subroutine phymeasure(cfg,icfg,wf,A,iA,lpg,lsc)
		complex(8) :: A(:,:),iA(:,:),vu(ne2,2),Y(2,2),pb,wf(:,:),lpg,lsc,a11
		integer :: cfg(:),icfg(:),i,n,m,j,l,inb,sg,cr(2),x1(2),x2(2),nr(2),dr(2),r(2)
		lpg=0d0
		select case(pgflag)
		case("sdw")
			do i=1,ne
				call square_one2two(cfg(i),Ns,x1)
				lpg=lpg+(-1d0)**((i-1)/ne2)*((-1d0)**(mod(x1(1)+x1(2),2)))
			enddo
		case("ddw")
			do i=ne+1,Ns2
				nr(1)=cfg(i)
				call square_one2two(nr(1),Ns,x1)
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
						lpg=lpg+abs(dimag(dconjg(pb)))
						!lpg=lpg+dimag(dconjg(pb))*(-1d0)**inb*((-1d0)**(mod(x1(1)+x1(2),2)))
					endif
				enddo
			enddo
		end select
		lsc=0d0
		return
		!do i=ne+1,Ns2
			!nr(1)=cfg(i)
			!do inb=1,4,1
				!if(icfg(neb(nr(1),inb,1))>ne) then
					!nr(2)=neb(nr(1),inb,1)
					!do l=1,ne2
						!call pair(nr(1),cfg(l+ne2),wf,vu(l,1))
						!call pair(cfg(l),nr(2),wf,vu(l,2))
					!enddo
					!call pair(nr(1),nr(2),wf,a11)
					!pb=0d0
					!do l=1,ne2
						!do m=1,ne2
							!pb=pb+vu(l,1)*iA(l,m)*vu(m,2)
						!enddo
					!enddo
					!pb=a11-pb
					!lsc=lsc+dconjg(pb)*(-1d0)**inb
				!endif
			!enddo
		!enddo
		!return
		do i=ne+1,Ns2
			nr(1)=cfg(i)
			do inb=1,4,1
				if(icfg(neb(nr(1),inb,1))>ne) then
					nr(2)=neb(nr(1),inb,1)
					do j=1,ne2
						r(1)=cfg(j)
						do n=1,4,1
							if(icfg(neb(r(1),n,1))>ne2.and.icfg(neb(r(1),n,1))<=ne) then
								r(2)=neb(r(1),n,1)
								!call square_one2two(abs(nr(1)-r(1)),Ns,x1)
								!if(sum(x1**2)<4d0.or.sum(x1**2)>25d0) then
		if(any(neb(nr(1),:,:)==r(1)).or.any(neb(nr(1),:,:)==r(2)).or.any(neb(nr(2),:,:)==r(1)).or.any(neb(nr(2),:,:)==r(2))) then
									cycle
								endif
								cr=(/icfg(r(1)),icfg(r(2))-ne2/)
								do l=1,ne2
									call pair(nr(1),cfg(l+ne2),wf,vu(l,1))
									call pair(cfg(l),nr(2),wf,vu(l,2))
								enddo
								call pair(nr(1),nr(2),wf,vu(cr(2),1))
								call pair(nr(1),nr(2),wf,vu(cr(1),2))
								call det(vu,cr,A,iA,Y,pb,3)
								lsc=lsc+dconjg(pb)*(-1)**(inb+n)
							endif
						enddo
					enddo
				endif
			enddo
		enddo
		lsc=lsc/Ns2
	end subroutine
	subroutine energy(cfg,icfg,wf,A,iA,El)
		complex(8) :: A(:,:),iA(:,:),vu(ne2,2),Y(2,2),pb,wf(:,:),El
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
						if(abs(k-cfg(i))>3*Ns(1)) then
							El=El+t(inb)*dconjg(pb)
						else
							El=El-t(inb)*dconjg(pb)
						endif
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
						El=El-0.5d0*DJ*dconjg(pb)
					endif
				enddo
			enddo
		enddo
	end subroutine
	subroutine checkcfg(icfg,s)
		integer :: i,s,icfg(:)
		do i=1,Ns2
			if(icfg(i)<=ne2) then
				write(*,"(' ●'$)")
			elseif(icfg(i)<=ne) then
				write(*,"(' ○'$)")
			else
				write(*,"('  '$)")
			endif
			if(mod(i,Ns(1))==0) then
				write(*,"(1X)")
			endif
		enddo
		call sleepqq(s)
		write(*,"(1X)")
	end subroutine
	subroutine variational(var)
		use lapack95, only: heevd,heevr
		complex(8) :: wf(Ns(2),Ns(1)*2),dwf(Ns(2),Ns(1)*2,vn),O(vn),S(vn,vn),g(vn),phyval(1),er(1),Ov(vn),Sv(vn,vn),gv(vn),Ev
		complex(8) :: tmp(vn,vn)
		real(8) :: var(vn),dvar(vn),pdvar(vn),eg(vn),dt=5d0
		integer :: n,i,j,k,info,sg(vn),cs
		open(20,file="../data/var_-sdw.dat")
		n=32
		pdvar=0d0
		cs=0
		sg=0
		do
			write(*,"(I3$)")ne
			if(var(1)<1d-3) then
				var(1)=1d-3
			endif
			write(*,"(es9.2$)")var
			call iniwf(var,wf,dwf)
			Ev=0d0
			!Ov=0d0
			Sv=0d0
			gv=0d0
			er=0d0
			!$OMP PARALLEL DO REDUCTION(+:Ev,Sv,gv,er) PRIVATE(phyval,O,S,g) SCHEDULE(STATIC)
			do k=1,n
				call mc(wf,dwf,5000,phyval,O,S,g,.true.)
				Ev=Ev+phyval(1)
				do i=1,vn
					do j=1,vn
						Sv(i,j)=Sv(i,j)+S(i,j)-dconjg(O(i))*O(j)
					enddo
				enddo
				gv=gv+g-real(phyval(1)*O)
				er=er+abs(phyval(1))**2
			enddo
			!$OMP END PARALLEL DO
			Ev=Ev/n
			Sv=Sv/n
			gv=2d0*gv/n
			er=sqrt(abs(er/n-abs(Ev)**2))
			!Sv(:,4:5)=0d0
			!Sv(4:5,:)=0d0
			Sv=Sv+diag((/1d0,1d0,1d0,1d0,1d0/)*2d-1)
			call heevd(Sv,eg,"V")
			!write(*,*)sum(abs(matmul(matmul(Sv,diag(eg)),dconjg(transpose(Sv)))-tmp)),eg(vn)/eg(1)
			!if(eg(1)<=0d0) then
				!write(*,*)"S matrix is not positive define or singular"
			!endif
			!gv(4)=0d0
			!gv(5)=0d0
			!dvar=real(matmul(matmul(matmul(Sv,diag(1d0/eg)),dconjg(transpose(Sv))),gv))
			do i=1,vn
				dvar(i)=0d0
				do j=1,vn
					do k=1,vn
						!if(eg(k)/eg(vn)<1d-3) then
							!cycle
						!endif
						dvar(i)=dvar(i)+real(Sv(i,k)*1d0/eg(k)*dconjg(Sv(j,k))*gv(j))
					enddo
				enddo
			enddo
			write(*,"(es9.2$)")real(gv),real(Ev),real(er)
			write(*,"(x)")
			write(20,"(I4$)")ne
			write(20,"(e14.5$)")var,real(gv),real(Ev),real(er)
			write(20,"(x)")
			!do i=1,vn
				!call find_cross(pdvar(i),dvar(i),sg(i))
			!enddo
			!if(all(sg/=0)) then
				cs=cs+1
			!endif
			!if(cs>5.or.sum(abs(dvar))<7d-4) then
			if(sum(abs(dvar))<1d-3.and.cs>40) then
				!if(sum(abs(gv))<1d-3) then
				write(*,*)"variational finish",sum(abs(dvar))
				write(20,"(x)")
				write(10,"(I4$)")ne
				write(10,"(e14.5$)")var,real(Ev),real(er)
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
	complex(8) :: wf(Ns(2),Ns(1)*2),dwf(Ns(2),Ns(1)*2,vn),O(vn),S(vn,vn),g(vn),phyval(3),sphy(3),er(3),tmp(vn,vn),z(vn,vn),&
		Ov(vn),Sv(vn,vn),gv(vn),Ev
	real(8) :: var(vn),dvar(vn),Eb,E2,dp(25)
	integer :: n,i,j,info,Nmc=50000
	!open(10,file="../data/2d.dat")
	open(10,file="../data/test.dat")
	open(30,file="../data/phyvar_-sdw.dat")
	open(40,file="../data/2d.dat")
	open(21,file="../data/phyvar_save.dat",status="old",action="read",access="direct",form="formatted",recl=181)
	call init_random_seed()
	call square(Ns,neb)
	n=1
	do i=25,0,-1
		var=(/2d-1,0d0,2d-1,0d0,0d0/)
		!read(21,rec=i+1,fmt="(I4,5e16.8)")ne,var
		ne=Ns2-i*2
		ne2=ne/2
		write(30,"(i4$)")ne
		call variational(var)
		write(*,"(i4$)")ne
		write(*,"(es9.2$)")var
		write(30,"(e16.8$)")var
		write(30,"(1x)")
		cycle
		call iniwf(var,wf,dwf)
		sphy=0d0
		er=0d0
		ncfg=0d0
		scfg=0d0
		!!$OMP PARALLEL DO REDUCTION(+:sphy,er) PRIVATE(phyval) SCHEDULE(STATIC)
		do j=1,n
			call mc(wf,dwf,Nmc,phyval,O,S,g,.false.)
			sphy=sphy+(/real(phyval(1)),abs(phyval(2)),real(phyval(3))/)
			er=er+abs(phyval)**2
		enddo
		!!$OMP END PARALLEL DO
		sphy=sphy/n
		er=sqrt(abs(er/n-abs(sphy)**2))
		ncfg=ncfg/Nmc
		scfg=scfg/Nmc
		write(*,"(es9.2$)")real(sphy)
		write(*,"(1x)")
		write(30,"(e16.8$)")real(sphy(1)),real(er(1)),real(sphy(2)),real(er(2)),abs(sphy(3)),real(er(3))
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
	enddo
end program
