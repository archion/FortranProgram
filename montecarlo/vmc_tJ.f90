module M_pmt
	use M_const
	implicit none
	integer, parameter :: Ns(2)=(/10,10/),Ns2=Ns(1)*Ns(2),ne=88,ne2=ne/2
	integer :: neb(Ns2,4,3)
	real(8), parameter :: t(3)=(/1d0,-0.25d0,0.1d0/),DJ=0.25d0,V=0.d0
	character(3) :: pgflag="sdw"
	!character(3) :: pgflag="ddw"
end module
module M_vmc
	use M_pmt
	use M_matrix
	use M_rd
	use M_utility
	use M_latt
	implicit none
contains
	!subroutine ini(var,wf)
		!complex(8) :: wf(0:,0:)
		!real(8) :: k(2),var(:),Ek,sck,tmp
		!integer :: i,j,ix(2),ik(2)
		!write(*,"(A)")"initial the wavefunction"
		!wf=0d0
		!do j=1,Ns2
			!call square_one2two(j,Ns,ik)
			!k=(/2d0*pi/Ns(1)*ik(2),2*pi/Ns(2)*ik(1)+pi/Ns(2)/)
			!Ek=-2d0*t(1)*(cos(k(1))+cos(k(2)))-var(2)
			!sck=(cos(k(1))-cos(k(2)))*var(1)
			!tmp=sck/(Ek+sqrt(Ek**2+sck**2))
			!do i=1,Ns2
				!call square_one2two(i,Ns,ix)
				!wf(ix(1),ix(2))=wf(ix(1),ix(2))+tmp*exp(img*(k(1)*ix(2)+k(2)*ix(1)))
				!if(isnan(real(wf(ix(1),ix(2))))) then
					!write(*,"(A)")"NAN in wf, quit!!"
					!stop
				!endif
			!enddo
		!enddo
		!!wf=wf/maxval(abs(wf))
	!end subroutine
	!subroutine pair(ri,rj,wf,a)
		!integer :: ri,rj,ri2(2),rj2(2),dr(2)
		!complex(8) :: wf(0:,0:),a
		!real(8) :: s
		!call square_one2two(ri,Ns,ri2)
		!call square_one2two(rj,Ns,rj2)
		!dr=ri2-rj2
		!s=sign(1d0,dr(1)+0.1d0)
		!dr(1)=mod(dr(1)+Ns(2),Ns(2))
		!dr(2)=mod(dr(2)+Ns(1),Ns(1))
		!a=wf(dr(1),dr(2))*s
	!end subroutine
	subroutine ini(var,wf)
		complex(8) :: wf(0:,0:),p(2)
		real(8) :: k(2),var(:),sck,a(2),ek(2),pgk,eks,eka,al
		integer :: i,j,ix(2),ik(2),n,np,pnp,eo
		!write(*,"(A$)")"initial the wavefunction ... "
		wf=0d0
		n=size(wf,2)/2
		al=0.01d0
		do
			np=0
			pnp=0
			do j=1,Ns2
				call square_one2two(j,Ns,ik)
				k=(/2d0*pi/Ns(1)*ik(2),2*pi/Ns(2)*ik(1)+pi/Ns(2)/)-(/pi,pi/)
				eks=-4d0*t(2)*cos(k(1))*cos(k(2))-2d0*t(3)*(cos(2d0*k(1))+cos(2d0*k(2)))-var(2)
				eka=-2d0*t(1)*(cos(k(1))+cos(k(2)))
				if((eka+eks)<0d0) then
					np=np+1
				endif
			enddo
			!write(*,*)np,var(2)
			!call sleepqq(1000)
			if((abs(ne2-np)/abs(ne2-pnp))>0.8d0) then
				al=max(al-0.05d0,0.001d0)
			endif
			if(abs(ne2-np)<1d-3) then
				exit
			endif
			var(2)=var(2)+al*(ne2-np)
			pnp=np
		enddo
		do j=1,Ns2
			call square_one2two(j,Ns,ik)
			k=(/2d0*pi/Ns(1)*ik(2),2*pi/Ns(2)*ik(1)+pi/Ns(2)/)-(/pi,pi/)
			eks=-4d0*t(2)*cos(k(1))*cos(k(2))-2d0*t(3)*(cos(2d0*k(1))+cos(2d0*k(2)))-var(2)
			eka=-2d0*t(1)*(cos(k(1))+cos(k(2)))
			gk=(cos(k(1))-cos(k(2)))
			if(abs(k(1))+abs(k(2))>pi) then
				cycle
			endif
			if(abs(var(1))<1d-5) then
				var(1)=sign(1d-5,var(1))
			endif
			sck=gk*var(1)
			select case(pgflag)
			case("ddw")
				pgk=gk*var(3)
			case("sdw")
				pgk=var(3)
			end select
			tmp=sqrt(pgk**2+eka**2)
			ek=eks+(/1d0,-1d0/)*sqrt(pgk**2+eka**2)
			E=sqrt(ek**2+sck**2)
			a=sck/(ek+E)*(/-1d0,1d0/)
			da(:,1)=gk/(ek+E)+sck*gk/((ek+E)**2*E)*(/1d0,-1d0/)
			da(:,2)=sck*gk*(1d0+ek/E)/(ek+E)**2*(/-1d0,1d0/)
			select case(pgflag)
			case("ddw")
				da(:,3)=da(:,2)*pgk*gk/tmp*(/-1d0,1d0/)
				u2=0.5d0*(1d0+eka/tmp)
				v2=0.5d0*(1d0-eka/tmp)
				uv=0.5d0*pgk/tmp*img
				du2=0.5d0*(1d0-eka*pgk*gk/tmp**3)
				dv2=0.5d0*(1d0+eka*pgk*gk/tmp**3)
				duv=0.5d0*(gk/tmp-pgk*pgk*gk/tmp**3)*img
			case("sdw")
				da(:,3)=da(:,2)*pgk/tmp*(/-1d0,1d0/)
				u2=0.5d0*(1d0+eka/tmp)
				v2=0.5d0*(1d0-eka/tmp)
				uv=0.5d0*pgk/tmp
				du2=0.5d0*(1d0-eka*pgk/tmp**3)
				dv2=0.5d0*(1d0+eka*pgk/tmp**3)
				duv=0.5d0*(1d0/tmp-pgk*pgk/tmp**3)
			end select
			do i=1,Ns2
				call square_one2two(i,Ns,ix)
				eo=(-1)**mod(ix(1)+ix(2),2)
				tmpe=(/u2-eo*v2,-v2+eo*u2/)
				dtmpe=(/du2-eo*dv2,-dv2+eo*du2/)
				tmpo=(-eo*dconjg(uv)+uv)
				dtmpo=(-eo*dconjg(duv)+duv)
				dwf(ix(1),ix(2),:)=dwf(ix(1),ix(2),:)+&
				(/&
				(tmpe(1)+tmpo(1))*da(1,1)+(tmpe(2)+tmpo(2))*da(2,1),&
				(tmpe(1)+tmpo(1))*da(1,2)+(tmpe(2)+tmpo(2))*da(2,2),&
				(tmpe(1)+tmpo(1))*da(1,3)+(tmpe(2)+tmpo(2))*da(2,3)+(dtmpe(1)+dtmpo(1))*a(1)+(dtmpe(2)+dtmpo(2))*a(2)&
				/)*exp(img*(k(1)*ix(2)+k(2)*ix(1)))
				dwf(ix(1),ix(2)+n,:)=dwf(ix(1),ix(2)+n,:)+&
				(/&
				(tmpe(1)-tmpo(1))*da(1,1)+(tmpe(2)-tmpo(2))*da(2,1),&
				(tmpe(1)-tmpo(1))*da(1,2)+(tmpe(2)-tmpo(2))*da(2,2),&
				(tmpe(1)-tmpo(1))*da(1,3)+(tmpe(2)-tmpo(2))*da(2,3)+(dtmpe(1)-dtmpo(1))*a(1)+(dtmpe(2)-dtmpo(2))*a(2)&
				/)*exp(img*(k(1)*ix(2)+k(2)*ix(1)))
				if(isnan(real(wf(ix(1),ix(2))))) then
					write(*,"(A)")"NAN in wf, quit!!"
					stop
				endif
			enddo
		enddo
			call square_one2two(j,Ns,ik)
			k=(/2d0*pi/Ns(1)*ik(2),2*pi/Ns(2)*ik(1)+pi/Ns(2)/)-(/pi,pi/)
			eks=-4d0*t(2)*cos(k(1))*cos(k(2))-2d0*t(3)*(cos(2d0*k(1))+cos(2d0*k(2)))-var(2)
			eka=-2d0*t(1)*(cos(k(1))+cos(k(2)))
			if(abs(k(1))+abs(k(2))>pi) then
				cycle
			endif
			if(abs(var(1))<1d-5) then
				var(1)=sign(1d-5,var(1))
			endif
			sck=(cos(k(1))-cos(k(2)))*var(1)
			select case(pgflag)
			case("ddw")
				pgk=(cos(k(1))-cos(k(2)))*var(3)
			case("sdw")
				pgk=var(3)
			end select
			ek=eks+(/1d0,-1d0/)*sqrt(pgk**2+eka**2)
			a=sck/(ek+sqrt(ek**2+sck**2))*(/-1d0,1d0/)
			if(abs(eka)<1d-8) then
				p=(/1d0,-1d0/)*sign(1d-8,eka)
			else
				p=eka/sqrt(pgk**2+eka**2)*(/1d0,-1d0/)
			endif
			select case(pgflag)
			case("ddw")
				p=(/sqrt(0.5d0*(1d0+p(1))),sqrt(0.5d0*(1d0+p(2)))*sign(1d0,pgk)*img/)
			case("sdw")
				p=(/sqrt(0.5d0*(1d0+p(1))),sqrt(0.5d0*(1d0+p(2)))*sign(1d0,pgk)/)
			end select
			do i=1,Ns2
				call square_one2two(i,Ns,ix)
				eo=(-1)**mod(ix(1)+ix(2),2)
				wf(ix(1),ix(2))=wf(ix(1),ix(2))+((p(1)**2+eo*(-p(2)*dconjg(p(2))-p(1)*dconjg(p(2)))+p(1)*p(2))*a(1)+&
					(-p(2)*dconjg(p(2))+eo*(p(1)**2-p(1)*dconjg(p(2)))+p(1)*p(2))*a(2))*exp(img*(k(1)*ix(2)+k(2)*ix(1)))
				wf(ix(1),ix(2)+n)=wf(ix(1),ix(2)+n)+((p(1)**2+eo*(-p(2)*dconjg(p(2))+p(1)*dconjg(p(2)))-p(1)*p(2))*a(1)+&
					(-p(2)*dconjg(p(2))+eo*(p(1)**2+p(1)*dconjg(p(2)))-p(1)*p(2))*a(2))*exp(img*(k(1)*ix(2)+k(2)*ix(1)))
				if(isnan(real(dwf(ix(1),ix(2))))) then
					write(*,"(A)")"NAN in wf, quit!!"
					stop
				endif
			enddo
		enddo
		!write(*,"(A)")"initial finish"
	end subroutine
	subroutine pair(ri,rj,wf,a)
		integer :: ri,rj,ri2(2),rj2(2),dr(2)
		complex(8) :: wf(0:,0:),a
		real(8) :: s
		call square_one2two(ri,Ns,ri2)
		call square_one2two(rj,Ns,rj2)
		dr=ri2-rj2
		s=sign(1d0,dr(1)+0.1d0)
		dr(1)=mod(dr(1)+Ns(2),Ns(2))
		dr(2)=mod(dr(2)+Ns(1),Ns(1))
		a=wf(dr(1),dr(2)+size(wf,2)/2*mod(ri2(1)+ri2(2),2))*s
	end subroutine
	subroutine iniwf(var,wf)
		complex(8) :: wf(0:,0:),p(2),dp(2,3)
		real(8) :: var(:),a(2),da(2,3)
		wf=0d0
		n=size(wf,2)/2
			al=0.01d0
			do
				np=0
				pnp=0
				do j=1,Ns2
					call square_one2two(j,Ns,ik)
					k=(/2d0*pi/Ns(1)*ik(2),2*pi/Ns(2)*ik(1)+pi/Ns(2)/)-(/pi,pi/)
					eks=-4d0*t(2)*cos(k(1))*cos(k(2))-2d0*t(3)*(cos(2d0*k(1))+cos(2d0*k(2)))-var(2)
					eka=-2d0*t(1)*(cos(k(1))+cos(k(2)))
					if((eka+eks)<0d0) then
						np=np+1
					endif
				enddo
				!write(*,*)np,var(2)
				!call sleepqq(1000)
				if((abs(ne2-np)/abs(ne2-pnp))>0.8d0) then
					al=max(al-0.05d0,0.001d0)
				endif
				if(abs(ne2-np)<1d-3) then
					exit
				endif
				var(2)=var(2)+al*(ne2-np)
				pnp=np
			enddo
		do j=1,Ns2
			call square_one2two(j,Ns,ik)
			k=(/2d0*pi/Ns(1)*ik(2),2*pi/Ns(2)*ik(1)+pi/Ns(2)/)-(/pi,pi/)
			eks=-4d0*t(2)*cos(k(1))*cos(k(2))-2d0*t(3)*(cos(2d0*k(1))+cos(2d0*k(2)))-var(2)
			eka=-2d0*t(1)*(cos(k(1))+cos(k(2)))
			gk=(cos(k(1))-cos(k(2)))
			if(abs(k(1))+abs(k(2))>pi) then
				cycle
			endif
			if(abs(var(1))<1d-5) then
				var(1)=sign(1d-5,var(1))
			endif
			sck=gk*var(1)
			select case(pgflag)
			case("ddw")
				pgk=gk*var(3)
			case("sdw")
				pgk=var(3)
			end select
			tmp=sqrt(pgk**2+eka**2)
			ek=eks+(/1d0,-1d0/)*sqrt(pgk**2+eka**2)
			E=sqrt(ek**2+sck**2)
			a=sck/(ek+E)*(/-1d0,1d0/)
			select case(pgflag)
			case("ddw")
				u2=0.5d0*(1d0+eka/tmp)
				v2=0.5d0*(1d0-eka/tmp)
				uv=0.5d0*pgk/tmp*img
			case("sdw")
				u2=0.5d0*(1d0+eka/tmp)
				v2=0.5d0*(1d0-eka/tmp)
				uv=0.5d0*pgk/tmp
			end select
			do i=1,Ns2
				call square_one2two(i,Ns,ix)
				eo=(-1)**mod(ix(1)+ix(2),2)
				tmpe=(/u2-eo*v2,-v2+eo*u2/)
				tmpo=(-eo*dconjg(uv)+uv)
				wf(ix(1),ix(2),1)=wf(ix(1),ix(2),1)+&
					((tmpe(1)+tmpo(1))*a(1)+(tmpe(2)+tmpo(2))*a(2))*exp(img*(k(1)*ix(2)+k(2)*ix(1)))
				wf(ix(1),ix(2)+n,1)=wf(ix(1),ix(2)+n,1)+&
					((tmpe(1)-tmpo(1))*a(1)+(tmpe(2)-tmpo(2))*a(2))*exp(img*(k(1)*ix(2)+k(2)*ix(1)))
				if(isnan(real(wf(ix(1),ix(2))))) then
					write(*,"(A)")"NAN in wf, quit!!"
					stop
				endif
			enddo
		enddo
	end subroutine
	subroutine inidwf(var,wf)
		complex(8) :: wf(0:,0:,:),p(2),dp(2,3)
		real(8) :: var(:),a(2),da(2,3)
		logical :: flag
		wf=0d0
		n=size(wf,2)/2
		do j=1,Ns2
			call square_one2two(j,Ns,ik)
			k=(/2d0*pi/Ns(1)*ik(2),2*pi/Ns(2)*ik(1)+pi/Ns(2)/)-(/pi,pi/)
			eks=-4d0*t(2)*cos(k(1))*cos(k(2))-2d0*t(3)*(cos(2d0*k(1))+cos(2d0*k(2)))-var(2)
			eka=-2d0*t(1)*(cos(k(1))+cos(k(2)))
			gk=(cos(k(1))-cos(k(2)))
			if(abs(k(1))+abs(k(2))>pi) then
				cycle
			endif
			if(abs(var(1))<1d-5) then
				var(1)=sign(1d-5,var(1))
			endif
			sck=gk*var(1)
			select case(pgflag)
			case("ddw")
				pgk=gk*var(3)
			case("sdw")
				pgk=var(3)
			end select
			tmp=sqrt(pgk**2+eka**2)
			ek=eks+(/1d0,-1d0/)*sqrt(pgk**2+eka**2)
			E=sqrt(ek**2+sck**2)
			a=sck/(ek+E)*(/-1d0,1d0/)
			da(:,1)=gk/(ek+E)+sck*gk/((ek+E)**2*E)*(/1d0,-1d0/)
			da(:,2)=sck*gk*(1d0+ek/E)/(ek+E)**2*(/-1d0,1d0/)
			select case(pgflag)
			case("ddw")
				da(:,3)=da(:,2)*pgk*gk/tmp*(/-1d0,1d0/)
				u2=0.5d0*(1d0+eka/tmp)
				v2=0.5d0*(1d0-eka/tmp)
				uv=0.5d0*pgk/tmp*img
				du2=0.5d0*(1d0-eka*pgk*gk/tmp**3)
				dv2=0.5d0*(1d0+eka*pgk*gk/tmp**3)
				duv=0.5d0*(gk/tmp-pgk*pgk*gk/tmp**3)*img
			case("sdw")
				da(:,3)=da(:,2)*pgk/tmp*(/-1d0,1d0/)
				u2=0.5d0*(1d0+eka/tmp)
				v2=0.5d0*(1d0-eka/tmp)
				uv=0.5d0*pgk/tmp
				du2=0.5d0*(1d0-eka*pgk/tmp**3)
				dv2=0.5d0*(1d0+eka*pgk/tmp**3)
				duv=0.5d0*(1d0/tmp-pgk*pgk/tmp**3)
			end select
			do i=1,Ns2
				call square_one2two(i,Ns,ix)
				eo=(-1)**mod(ix(1)+ix(2),2)
				tmpe=(/u2-eo*v2,-v2+eo*u2/)
				dtmpe=(/du2-eo*dv2,-dv2+eo*du2/)
				tmpo=(-eo*dconjg(uv)+uv)
				dtmpo=(-eo*dconjg(duv)+duv)
				wf(ix(1),ix(2),:)=wf(ix(1),ix(2),:)+&
					(/&
					(tmpe(1)+tmpo(1))*da(1,1)+(tmpe(2)+tmpo(2))*da(2,1),&
					(tmpe(1)+tmpo(1))*da(1,2)+(tmpe(2)+tmpo(2))*da(2,2),&
					(tmpe(1)+tmpo(1))*da(1,3)+(tmpe(2)+tmpo(2))*da(2,3)+(dtmpe(1)+dtmpo(1))*a(1)+(dtmpe(2)+dtmpo(2))*a(2)&
					/)*exp(img*(k(1)*ix(2)+k(2)*ix(1)))
				wf(ix(1),ix(2)+n,:)=wf(ix(1),ix(2)+n,:)+&
					(/&
					(tmpe(1)-tmpo(1))*da(1,1)+(tmpe(2)-tmpo(2))*da(2,1),&
					(tmpe(1)-tmpo(1))*da(1,2)+(tmpe(2)-tmpo(2))*da(2,2),&
					(tmpe(1)-tmpo(1))*da(1,3)+(tmpe(2)-tmpo(2))*da(2,3)+(dtmpe(1)-dtmpo(1))*a(1)+(dtmpe(2)-dtmpo(2))*a(2)&
					/)*exp(img*(k(1)*ix(2)+k(2)*ix(1)))
				if(isnan(real(wf(ix(1),ix(2))))) then
					write(*,"(A)")"NAN in wf, quit!!"
					stop
				endif
			enddo
		enddo
	end subroutine
	subroutine mfO(iA,dA,O)
		complex(8) :: iA(:,:),dA(:,:,:),O(:)
		integer :: n,i,j,k
		n=size(iA,1)
		O=0d0
		do i=1,n
			do j=1,n
				O=O+iA(i,j)*dA(j,i,:)
			enddo
		enddo
	end subroutine
	subroutine mc(wf,Nmc,phyval,S,g,flag)
		complex(8) :: pb,A(ne2,ne2),iA(ne2,ne2),vu(ne2,2),wf(:,:),Y(2,2),phyval(:),lphy(size(phyval,1)),S(:,:),g(:)
		real(8) :: rpb
		integer :: cfg(Ns2),icfg(Ns2),i,j,l,sg,n,Nmc,Nhot,cr(2),acp
		logical :: flag
		!write(*,"(A)")"start monte carlo"
		n=0
		Nhot=Nmc
		do
			if(n==0) then
				E=0d0
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
					enddo
				enddo
				iA=A
				call matrix_inv(iA)
			endif
			n=n+1
			if(n==Nhot) then
				if(flag) then
					call energy(cfg,icfg,wf,A,iA,lphy(1))
				else
					call calSg(lS,lg)
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
				enddo
			elseif(j>ne2.and.i>ne2) then
				sg=2
				cr=(/-1,i-ne2/)
				j=j+ne2
				do l=1,ne2
					call pair(cfg(l),cfg(j),wf,vu(l,2))
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
				enddo
				call pair(cfg(j),cfg(i),wf,vu(cr(2),1))
				call pair(cfg(j),cfg(i),wf,vu(cr(1),2))
			endif
			call det(vu,cr,A,iA,Y,pb,sg)
			!call checkcfg(icfg,80)
			!write(*,*)A
			if(abs(pb)>1d10) then
				write(*,*)"initial configure may be not good"
				n=0
				cycle
			endif
			call random_number(rpb)
			if(rpb<abs(pb)**2) then
				acp=acp+1
				call swap(icfg(cfg(i)),icfg(cfg(j)))
				call swap(cfg(i),cfg(j))
				call update(vu,cr,Y,pb,A,iA,sg)
				!call checkcfg(icfg,80)
				!if(mod(n,500)==0) then
					!call checkinv(A,iA,rpb)
					!write(*,"(I9,E14.3$)")n,rpb
					!iA=A
					!call matrix_inv(iA)
					!call checkinv(A,iA,rpb)
					!write(*,"(E14.3)")rpb
				!endif
				if(n>Nhot) then
					if(flag) then
						call energy(cfg,icfg,wf,A,iA,lphy(1))
					else
						call calcuSg(iA,dA,lS,lg)
					endif
				endif
			endif
			if(n>Nhot) then
				phyval=phyval+lphy
				S=S+lS
				g=g+lg
			endif
			if(n>=(Nmc+Nhot)) then
				phyval=phyval/(Nmc*Ns2)
				S=S/(Nmc*Ns2)
				g=g/(Nmc*Ns2)
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
	subroutine update(vu,cr,Y,pb,A,iA,sg)
		complex(8) :: vu(:,:),A(:,:),iA(:,:),Y(:,:),pb
		integer :: cr(:),sg
		if(sg==3) then
			call inv_update_rowcol(vu,cr,Y,A,iA)
		elseif(sg==1) then
			call inv_update_row(vu(:,1),cr(1),pb,A,iA)
		else
			call inv_update_col(vu(:,2),cr(2),pb,A,iA)
		endif
	end subroutine
	subroutine energy(cfg,icfg,wf,A,iA,El)
		complex(8) :: A(:,:),iA(:,:),vu(ne2,2),Y(2,2),pb,wf(:,:),El
		!real(8) :: 
		integer :: cfg(:),icfg(:),i,ii,j,k,l,inb,sg,cr(2)
		El=0d0
		do i=1,ne
			do ii=1,4
				do inb=1,1
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
						if(abs(k-cfg(i))>2*Ns(1)) then
							El=El+t(1)*dconjg(pb)
						else
							El=El-t(1)*dconjg(pb)
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
	subroutine calcuSg(iA,dA,lS,lg)
		call mfO(iA,dA,O)
		lg=
		lS=
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
end module
program main
	use M_vmc
	complex(8) :: wf(Ns(2),Ns(1)*2),E,S(3,3),g(3),phyval(1)
	real(8) :: var(3)=(/0d0,0d0,0d0/),dvar(3),Eb,E2
	integer :: n,i,j,d
	open(10,file="../data/2d.dat")
	call init_random_seed()
	call square(Ns,neb)
	n=12
	do
		call ini(var,wf)
		call derive(var,dwf)
		call mc(wf,dwf,10000,phyval,S,g,.false.)
		call conjgrad(S,g,dvar)
		if(sum(abs(dvar))<1d-5) then
			exit
		endif
		var=var+dvar
	enddo
	!$OMP PARALLEL DO REDUCTION(+:sphy,err) PRIVATE(phyval) SCHEDULE(STATIC)
	do j=1,n
		call mc(wf,dwf,50000,phyval,S,g,.true.)
		sphy=sphy+phyval
		err=err+abs(phyval)**2
	enddo
	!$OMP END PARALLEL DO
	sphy=sphy/n
	err=sqrt(abs(err/n-sphy))
	write(*,"(e16.8$)")(sphy(i),err(i),i=1,size(sphy))
	write(*,"(1x)")
	write(10,"(e16.8$)")(sphy(i),err(i),i=1,size(sphy))
	write(10,"(1x)")
end program
