module pmt
	use M_const
	implicit none
	real(8), parameter :: t(5)=(/1d0,-0.25d0,0.1d0,0d0,0d0/),&
	!real(8), parameter :: t(5)=(/1d0,0d0,0d0,0d0,0d0/),&
		cvg=1e-6,&
		V=0.12d0,DJ=0.35d0,&
		!V=0.d0,DJ=0.25d0,&
	Vs=DJ-V,Vd=0.5d0*DJ+V
	character(3) :: pgflag="ddw"
	!character(3) :: pgflag="sdw"
end module
module selfcons
	use pmt
	use M_utility
	use M_matrix
	use lapack95, only: heevd
	implicit none
	integer, parameter :: mk=512
contains
	subroutine selfconsist_tg(nf,Tk,sc,pg,ap,cp)
		real(8) :: Tk,nf
		complex(8) :: Uk(4,4)
		real(8) :: k(2),np,pnp,pg,pgp,ek(4),cp,sc,scp,cvg1,al,ap,app
		integer :: i,j,c,sg
		cvg1=0.0001
		pg=max(pg,0.1d0)
		sc=max(sc,0.1d0)
		ap=max(ap,0.1d0)
		do 
			pnp=0d0
			al=0.2d0
			do 			
				pgp=0d0
				scp=0d0
				app=0d0
				np=0d0
				c=0
				!$OMP PARALLEL DO REDUCTION(+:np,pgp,scp,app,c) PRIVATE(k,ek,Uk)
				!do i=0,mk
					!do j=0,min(i,mk-i)
						!k=(/pi/mk*i,pi/mk*j/)
				do i=-mk/2,mk/2-1
					do j=-mk/2+abs(i),mk/2-1-abs(i)
						c=c+1
						k=2d0*pi/mk*(/i,j/)
						call EU2(k,nf,pg,sc,ap,cp,ek,Uk)
						call order(k,ek,Uk,Tk,np,pgp,scp,app)
					enddo
				enddo
				!$OMP END PARALLEL DO
				c=c*2
				!np=np/(mk**2)*2d0
				!pgp=pgp/(mk**2)*2d0
				!scp=scp/(mk**2)*2d0
				!app=app/(mk**2)*2d0
				np=np/c
				pgp=pgp/c
				scp=scp/c
				app=app/c
				if(abs(nf-np)<1d-5) then
					exit
				endif
				call find_cross(pnp,np-nf,sg)
				if(sg/=0) then
					al=al*0.3d0
				endif
				cp=cp-al*sign(1d0,np-nf)
			enddo
			if((abs(scp-sc)+abs(pgp-pg)+abs(app-ap))<cvg) then
				exit
			endif
			sc=scp
			pg=pgp
			ap=app
		enddo
		!write(*,"(e12.4$)")np,Tk,cp,sc,pg,ap
		!write(*,"(1X)")
	end subroutine
	subroutine EU1(k,nf,pg,sc,ap,cp,ek,Uk)
		complex(8) :: Uk(4,4),cth,sth,cfy(2),sfy(2)
		real(8) :: eka,eks,e1(2),e2(2),ek(4),k(2),pg,pgk,cp,sck,sc,gk,th,fy(2),cos2th,sin2th,cos2fy(2),sin2fy(2),ap,dp,nf
		dp=abs(1d0-nf)
		!dp=abs(1d0)
		eks=-4d0*dp*t(2)*cos(k(1))*cos(k(2))-cp-2d0*dp*t(3)*(cos(2d0*k(1))+cos(2d0*k(2)))
		eka=-2d0*(dp*t(1)+ap*Vd)*(cos(k(1))+cos(k(2)))
		gk=0.5d0*(cos(k(1))-cos(k(2)))
		select case(pgflag)
		case("ddw")
			pgk=pg*gk*4d0*Vd
		case("sdw")
			pgk=-pg*DJ*2d0
		end select
		sck=sc*gk*Vs*4d0
		if(abs(sck)<1d-8) then
			sck=1d-8
		endif
		e1=eks+(/1d0,-1d0/)*sqrt(pgk**2+eka**2)*sign(1d0,eka)
		e2=sqrt(e1**2+sck**2)
		ek=(/e2*sign(1d0,e1),-e2*sign(1d0,e1)/)
		cos2th=eka/sqrt(pgk**2+eka**2)*sign(1d0,eka)
		sin2th=pgk/sqrt(pgk**2+eka**2)*sign(1d0,eka)
		cos2fy=e1/e2*sign(1d0,e1)
		sin2fy=sck/e2*sign(1d0,e1)
		cth=dcmplx(sqrt(0.5d0*(1d0+cos2th)))
		select case(pgflag)
		case("ddw")
			sth=img*dcmplx(sqrt(0.5d0*(1d0-cos2th))*sign(1d0,sin2th))
		case("sdw")
			sth=dcmplx(sqrt(0.5d0*(1d0-cos2th))*sign(1d0,sin2th))
		end select
		cfy=dcmplx(sqrt(0.5d0*(1d0+cos2fy)))
		sfy=dcmplx(sqrt(0.5d0*(1d0-cos2fy))*sign(1d0,sin2fy))
		Uk=reshape((/   cth*cfy(1) ,  dconjg(sth)*cfy(1) ,  cth*sfy(1) , -dconjg(sth)*sfy(1) ,   & 
			-sth*cfy(2) ,          cth*cfy(2) , -sth*sfy(2) ,         -cth*sfy(2) ,   & 
			-cth*sfy(1) , -dconjg(sth)*sfy(1) ,  cth*cfy(1) , -dconjg(sth)*cfy(1) ,   & 
			-sth*sfy(2) ,          cth*sfy(2) ,  sth*cfy(2) ,   cth*cfy(2)     /) , (/4 , 4/))
	end subroutine
	subroutine EU_lapack(k,nf,pg,sc,ap,cp,ek,Uk)
		complex(8) :: Uk(4,4),pgk,sck,eka,eks,Hk(4,4)
		real(8) :: ek(4),k(2),pg,cp,sc,gk,ap,dp,nf
		integer :: info
		dp=abs(1d0-nf)
		eks=-4d0*dp*t(2)*cos(k(1))*cos(k(2))-cp-2d0*dp*t(3)*(cos(2d0*k(1))+cos(2d0*k(2)))
		eka=-2d0*(dp*t(1)+ap*Vd)*(cos(k(1))+cos(k(2)))
		gk=0.5d0*(cos(k(1))-cos(k(2)))
		select case(pgflag)
		case("ddw")
			pgk=pg*gk*4d0*Vd*img
		case("sdw")
			pgk=-pg*DJ*2d0
		end select
		sck=sc*gk*Vs*4d0
		Uk=reshape((/ eks+eka,conjg(pgk) , sck , dcmplx(0d0) ,   & 
					 pgk , eks-eka , dcmplx(0d0)  , -sck, & 
					 sck, dcmplx(0d0) , -eks-eka , conjg(pgk) , & 
					 dcmplx(0d0), -sck , pgk , -eks+eka /) , (/4 , 4/))
		!Hk=Uk
		!call EU2(k,nf,pg,sc,ap,cp,ek,Uk)
		call heevd(Uk,ek,"V","U",info)
		!write(40,*)sum(abs(matmul(conjg(transpose(Uk)),matmul(Hk,Uk))-diag(ek)))
	end subroutine
	subroutine EU2(k,nf,pg,sc,ap,cp,ek,Uk)
		complex(8) :: Uk(4,4),cth,sth,cfy(2),sfy(2)
		real(8) :: eka,eks,e1(2),e2(2),ek(4),k(2),pg,pgk,cp,sck,sc,gk,th,fy(2),cos2th,sin2th,cos2fy(2),sin2fy(2),ap,dp,nf
		dp=abs(1d0-nf)
		!dp=abs(1d0)
		eks=-4d0*dp*t(2)*cos(k(1))*cos(k(2))-cp-2d0*dp*t(3)*(cos(2d0*k(1))+cos(2d0*k(2)))
		eka=-2d0*(dp*t(1)+ap*Vd)*(cos(k(1))+cos(k(2)))
		gk=0.5d0*(cos(k(1))-cos(k(2)))
		select case(pgflag)
		case("ddw")
			pgk=pg*gk*4d0*Vd
		case("sdw")
			pgk=-pg*DJ*2d0
		end select
		sck=sc*gk*Vs*4d0
		if(abs(sck)<1d-8) then
			sck=1d-8
		endif
		e1=eks+(/1d0,-1d0/)*sqrt(pgk**2+eka**2)
		e2=sqrt(e1**2+sck**2)
		!ek=(/e2*sign(1d0,e1),-e2*sign(1d0,e1)/)
		ek=(/e2,-e2/)
		cos2th=eka/sqrt(pgk**2+eka**2)
		sin2th=pgk/sqrt(pgk**2+eka**2)
		if(abs(eka)<1d-8) then
			cos2th=0d0
			sin2th=sign(1d0,pgk)
		endif
		!cos2fy=e1/e2*sign(1d0,e1)
		!sin2fy=sck/e2*sign(1d0,e1)
		cos2fy=e1/e2
		sin2fy=sck/e2
		cth=dcmplx(sqrt(0.5d0*(1d0+cos2th)))
		select case(pgflag)
		case("ddw")
			sth=img*dcmplx(sqrt(0.5d0*(1d0-cos2th))*sign(1d0,sin2th))
		case("sdw")
			sth=dcmplx(sqrt(0.5d0*(1d0-cos2th))*sign(1d0,sin2th))
		end select
		cfy=dcmplx(sqrt(0.5d0*(1d0+cos2fy)))
		sfy=dcmplx(sqrt(0.5d0*(1d0-cos2fy))*sign(1d0,sin2fy))
		Uk=reshape((/   cth*cfy(1) ,  dconjg(sth)*cfy(1) ,  cth*sfy(1) , -dconjg(sth)*sfy(1) ,   & 
			-sth*cfy(2) ,          cth*cfy(2) , -sth*sfy(2) ,         -cth*sfy(2) ,   & 
			-cth*sfy(1) , -dconjg(sth)*sfy(1) ,  cth*cfy(1) , -dconjg(sth)*cfy(1) ,   & 
			-sth*sfy(2) ,          cth*sfy(2) ,  sth*cfy(2) ,   cth*cfy(2)     /) , (/4 , 4/))
	end subroutine
	subroutine order(k,ek,Uk,Tk,n,pg,sc,ap)
		complex(8) :: Uk(4,4)
		real(8) :: k(2),ek(4),pg,sc,gk,n,Tk,fk(4),bt,ap
		bt=1d0/Tk
		fk=1d0/(1d0+exp(bt*ek))
		gk=0.5d0*(cos(k(1))-cos(k(2)))
		select case(pgflag)
		case("ddw")
			pg=pg+dimag(dot_product(Uk(1,:)*dconjg(Uk(2,:))-Uk(2,:)*dconjg(Uk(1,:)),fk)*gk)
		case("sdw")
			pg=pg+0.5d0*dot_product(Uk(2,:)*dconjg(Uk(1,:))+Uk(1,:)*dconjg(Uk(2,:))+&
				Uk(4,:)*dconjg(Uk(3,:))+Uk(3,:)*dconjg(Uk(4,:)),fk)
		end select
		n=n+2d0+dot_product(Uk(1,:)*dconjg(Uk(1,:))+Uk(2,:)*dconjg(Uk(2,:))-&
			Uk(3,:)*dconjg(Uk(3,:))-Uk(4,:)*dconjg(Uk(4,:)),fk)
		sc=sc+dot_product(-Uk(1,:)*dconjg(Uk(3,:))+Uk(2,:)*dconjg(Uk(4,:)),fk)*gk
		ap=ap+dot_product(Uk(1,:)*dconjg(Uk(1,:))-Uk(2,:)*dconjg(Uk(2,:)),fk)*0.5d0*(cos(k(1))+cos(k(2)))
	end subroutine
	subroutine freeenergy(pg,sc,ap,cp,nf,Tk,feg,ddf)
		real(8) :: k(2),eka,eks,e1,e2(2),E(2),Tk,bt,dp,feg,pg,sc,ap,cp,gk,pgk,sck,ddf(2),tmp(2),s,nf
		integer :: i,j,l,Ns
		dp=abs(1d0-nf)
		bt=1d0/Tk
		feg=0d0
		ddf=0d0
		Ns=0
		!$OMP PARALLEL DO REDUCTION(+:feg,ddf,Ns) PRIVATE(k,eks,eka,gk,sck,pgk,e1,e2,E,tmp,s)
		do i=1,mk
			do j=0,min(i,mk-i)
				k=(/pi/mk*i,pi/mk*j/)
				eks=-4d0*dp*t(2)*cos(k(1))*cos(k(2))-cp-2d0*dp*t(3)*(cos(2d0*k(1))+cos(2d0*k(2)))
				eka=-2d0*(dp*t(1)+ap*Vd)*(cos(k(1))+cos(k(2)))
				gk=(cos(k(1))-cos(k(2)))
				select case(pgflag)
				case("ddw")
					pgk=pg*gk*2d0*Vd
				case("sdw")
					pgk=-pg*DJ*2d0
				end select
				sck=sc*gk*Vs*2d0
				e1=sqrt(pgk**2+eka**2)
				e2=eks+(/1d0,-1d0/)*e1
				E=sqrt(e2**2+sck**2)
				do l=1,2
					s=-(-1)**l
					tmp(1)=tanh(bt*E(l)/2d0)/E(l)
					tmp(2)=cosh(bt*E(l)/2d0)
					if(bt*E(l)>40d0) then
						feg=feg-E(l)
					else
						feg=feg-2d0*Tk*log(2d0*tmp(2))
					endif
					ddf(1)=ddf(1)&
						+(tmp(1)-bt/(2d0*tmp(2)**2))&
						*((2d0*Vs*gk)**2*sc/E(l))**2&
						-tmp(1)*(2d0*Vs*gk)**2+8d0*Vs
					ddf(2)=ddf(2)&
						+(tmp(1)*((e2(l)/E(l))**2-1+s*e2(l)/e1)-bt/(2d0*tmp(2)**2)*(e2(l)/E(l))**2)&
						*((2d0*Vd*gk)**2*pg/e1)**2&
						-s*tmp(1)*e2(l)/e1*(2d0*Vd*gk)**2+8d0*Vd
					Ns=Ns+1
				enddo
			enddo
		enddo
		!$OMP END PARALLEL DO
		select case(pgflag)
		case("ddw")
			feg=feg/Ns+(4d0*(Vd*(pg**2+ap**2)+Vs*sc**2)-cp)+cp*nf
		case("sdw")
			feg=feg/Ns+(4d0*(DJ*pg**2+Vd*ap**2+Vs*sc**2)-cp)+cp*nf
		end select
		ddf=sign(1d0,ddf)
	end subroutine
	subroutine freeorder(nf,Tk,pg,sc,ap,cp,np,pgp,scp,app)
		real(8) :: k(2),eka,eks,e1,e2(2),E(2),Tk,bt,dp,tmp,pg,sc,ap,cp,np,pgp,scp,app,pgk,sck,gk,s,nf
		integer :: i,j,l,Ns
		dp=abs(1d0-nf)
		bt=1d0/Tk
		Ns=0
		!$OMP PARALLEL DO REDUCTION(+:np,pgp,scp,app,Ns) PRIVATE(k,eks,eka,gk,sck,pgk,e1,e2,E,tmp,s)
		do i=1,mk
			do j=0,min(i,mk-i)
				k=(/pi/mk*i,pi/mk*j/)
				eks=-4d0*dp*t(2)*cos(k(1))*cos(k(2))-cp-2d0*dp*t(3)*(cos(2d0*k(1))+cos(2d0*k(2)))
				eka=-2d0*(dp*t(1)+ap*Vd)*(cos(k(1))+cos(k(2)))
				gk=(cos(k(1))-cos(k(2)))
				select case(pgflag)
				case("ddw")
					pgk=pg*gk*2d0*Vd
				case("sdw")
					pgk=-pg*DJ*2d0
				end select
				sck=sc*gk*Vs*2d0
				e1=sqrt(pgk**2+eka**2)
				e2=eks+(/1d0,-1d0/)*e1
				E=sqrt(e2**2+sck**2)
				do l=1,2
					s=-(-1)**l
					if(abs(e1)<1d-8) then
						e1=1d-8
					endif
					if(abs(E(l))<1d-8) then
						E(l)=1d-8
					endif
					tmp=tanh(bt*E(l)/2d0)/E(l)
					scp=scp+tmp*sck*gk*0.25d0
					tmp=tmp*e2(l)
					np=np+tmp
					tmp=tmp/e1*s
					select case(pgflag)
					case("ddw")
						pgp=pgp+tmp*gk*pgk*0.25d0
					case("sdw")
						pgp=pgp-tmp*pgk*0.5d0
					end select
					app=app-tmp*eka*(cos(k(1))+cos(k(2)))*0.25d0
					Ns=Ns+1
				enddo
			enddo
		enddo
		!$OMP END PARALLEL DO
		np=1d0-np/Ns
		pgp=pgp/Ns
		scp=scp/Ns
		app=app/Ns
	end subroutine
	subroutine selforder(pg,sc,ap,cp,nf,Tk,sig)
		real(8) :: k(2),pg,sc,ap,cp,pgp,scp,app,np,pnp,al,Tk,nf
		integer :: c,sig,sg
		if(Tk<0d0) then
			write(*,*)"Tk<0"
			stop
		endif
		pg=max(pg,0.1d0)
		sc=max(sc,0.1d0)
		!sc=0d0
		!pg=0d0
		ap=max(ap,0.1d0)
		do 
			al=0.2d0
			pnp=0d0
			do 			
				pgp=0d0
				scp=0d0
				app=0d0
				np=0d0
				c=c+1
				call freeorder(nf,Tk,pg,sc,ap,cp,np,pgp,scp,app)
				!call freeorder(0d0,Tk,pg,sc,0d0,cp,np,pgp,scp,app)
				!app=0d0
				if(isnan(np)) then
					write(*,*)"NAN"
					exit
				endif
				!write(*,"(3(a6,e12.4))")"cp=",cp,"n=",np,"al=",al
				if(abs(nf-np)<1d-5) then
					exit
				endif
				call find_cross(pnp,np-nf,sg)
				if(sg/=0) then
					al=al*0.3d0
				endif
				cp=cp-al*sign(1d0,np-nf)
			enddo
			!if(sig==0) then
				!scp=sc
				!pgp=pg
			!endif
			if((abs(scp-sc)+abs(pgp-pg)+abs(app-ap))<cvg*3d0) then
				exit
			endif
			sc=(scp*0.4d0+sc*0.6d0)
			pg=(pgp*0.4d0+pg*0.6d0)
			ap=(app*0.4d0+ap*0.6d0)
		enddo
		!write(*,"(e12.4$)")np,Tk,cp,sc,pg,ap
		!write(*,"(1X)")
		!call sleepqq(1000)
	end subroutine
end module
module phys
	use pmt
	use M_utility
	use M_matrix
	use selfcons
	implicit none
contains
	function superfluid(pg,sc,ap,cp,nf,Tk)
		integer :: i,j,l,m1,m2,psg(2),sg,mr
		complex(8) :: Uk(4,4),Upk(4,4),tmp,id,gm(4,4),mm(4,4)
		real(8) :: k(2),cp,bt,sc,pg,ek(4),dfk(4),fk(4),omg,domg,Tk,DJp,ap,dp,nf,&
			Rs,superfluid,eka,eks
		superfluid=0d0
		mr=512
		bt=1d0/Tk
		dp=abs(1d0-nf)
		!dp=1d0
		!$OMP PARALLEL DO REDUCTION(+:superfluid) PRIVATE(k,gm,ek,Uk,fk,Upk,dfk,mm)
		!do i=-mr/2,mr/2-1
			!do j=-mr/2+abs(i),mr/2-1-abs(i)
				!k=2d0*pi/mr*(/i,j/)
		do i=0,mr
			k(1)=pi/mr*i
			do j=0,mr-i
				k(2)=pi/mr*j
				gm=0d0
				mm=0d0
				gm(1,1)=(4d0*t(2)*dp*cos(k(1))*sin(k(2))+4d0*t(3)*dp*sin(2d0*k(2))+2d0*(t(1)*dp+ap*Vd)*sin(k(2)))
				gm(2,2)=(4d0*t(2)*dp*cos(k(1))*sin(k(2))+4d0*t(3)*dp*sin(2d0*k(2))-2d0*(t(1)*dp+ap*Vd)*sin(k(2)))
				mm(1,1)=(4d0*t(2)*dp*cos(k(1))*cos(k(2))+8d0*t(3)*dp*cos(2d0*k(2))+2d0*(t(1)*dp+ap*Vd)*cos(k(2)))
				mm(2,2)=(4d0*t(2)*dp*cos(k(1))*cos(k(2))+8d0*t(3)*dp*cos(2d0*k(2))-2d0*(t(1)*dp+ap*Vd)*cos(k(2)))
				!superfluid=superfluid+real(mm(1,1)+mm(2,2))
				if(pgflag=="ddw") then
					gm(1,2)=img*pg*2d0*Vd*sin(k(2))
					mm(1,2)=img*pg*2d0*Vd*cos(k(2))
					gm(2,1)=-gm(1,2)
					gm(3,4)=-gm(1,2)
					gm(4,3)=-gm(3,4)
					mm(2,1)=-mm(1,2)
					mm(3,4)=mm(1,2)
					mm(4,3)=-mm(3,4)
				endif
				gm(3,3)=gm(1,1)
				gm(4,4)=gm(2,2)
				mm(3,3)=-mm(1,1)
				mm(4,4)=-mm(2,2)
				call EU2(k,1-dp,pg,sc,ap,cp,ek,Uk)
				!call EU_lapack(k,1-dp,pg,sc,ap,cp,ek,Uk)
				fk=1d0/(1d0+exp(bt*ek))
				dfk=-bt/((1d0+cosh(bt*ek))*2)
				gm=matmul(transpose(conjg(Uk)),matmul(gm,Uk))
				mm=matmul(transpose(conjg(Uk)),matmul(mm,Uk))
				do m1=1,4
					do m2=1,4
						if(m1==m2.or.abs(ek(m1)-ek(m2))<1d-20) then
							cycle
						endif
						superfluid=superfluid+real(gm(m1,m2)*conjg(gm(m1,m2))*(fk(m1)-fk(m2))/(ek(m1)-ek(m2)))
					enddo
					superfluid=superfluid+real(gm(m1,m1)*conjg(gm(m1,m1))*dfk(m1)+mm(m1,m1)*fk(m1))
				enddo
			enddo
		enddo
		!$END OMP PARALLEL DO
		superfluid=superfluid/(mr**2)
	end function
	subroutine raman(pg,sc,ap,cp,nf,Tk,omgr,pk)
		integer :: i,j,l,m1,m2,psg(2),sg,mr
		complex(8) :: Uk(4,4),tmp,id
		real(8) :: k(2),cp,bt,sc,pg,ek(4),gm(4,2),Tr(2),fk(4),omg,omgr(:),domg,Tk,peak(3,2),R_rpa,pk(2),DJp,ap,dp,R(2),pR(2,2),nf,&
			Rs(4,2),pRs(2,4,2)
		!id=(0d0,0.003d0)
		id=(0d0,0.04d0)
		mr=512
		pR=0d0
		pRs=0d0
		dp=abs(1d0-nf)
		bt=1d0/Tk
		omg=omgr(1)
		domg=omgr(3)
		DJp=DJ/(dp*t(1)+DJ*ap)**2
		do while(omg<omgr(2))
			Rs=0d0
			!$OMP PARALLEL DO REDUCTION(+:Rs) PRIVATE(k,gm,ek,Uk,fk,Tr,tmp)
			do i=0,mr
				k(1)=pi/mr*i
				do j=0,mr-i
					k(2)=pi/mr*j
					gm(1:2,1)=(cos(k(1))-cos(k(2)))*(dp*t(1)+ap*Vd)*(/1d0,-1d0/)+&
						4d0*t(3)*dp*(cos(2d0*k(1))-cos(2d0*k(2)))
					gm(1:2,2)=sin(k(1))*sin(k(2))*dp*t(2)*(/-4d0,-4d0/)
					gm(3:4,:)=-gm(1:2,:)
					call EU2(k,nf,pg,sc,ap,cp,ek,Uk)
					fk=1d0/(1d0+exp(bt*ek))
					!do m1=1,4
						!do m2=1,4
							!Tr=(/abs(sum(gm(:,1)*dconjg(Uk(:,m2))*Uk(:,m1)))**2,abs(sum(gm(:,2)*dconjg(Uk(:,m2))*Uk(:,m1)))**2/)
							!!R=R+Tr*(fk(m1)-fk(m2))*1d0/(omg+ek(m1)-ek(m2)+img*0.04d0)
							!tmp=omg+ek(m1)-ek(m2)+id
							!!tmp=omg+ek(m1)-ek(m2)+img*0.001d0
							!!tmp=omg+ek(m1)-ek(m2)+img*0.04d0
							!R=R+Tr*(fk(m1)-fk(m2))*dimag(tmp)/(dreal(tmp)**2+dimag(tmp)**2)
						!enddo
					!enddo
					do l=1,2
						Rs(1,l)=Rs(1,l)+dimag(&
							abs(sum(gm(:,l)*dconjg(Uk(:,1))*Uk(:,3)))**2*(fk(1)-fk(3))*1d0/(omg+ek(3)-ek(1)+id)+&
							abs(sum(gm(:,l)*dconjg(Uk(:,3))*Uk(:,1)))**2*(fk(3)-fk(1))*1d0/(omg+ek(1)-ek(3)+id)&
							)
						Rs(2,l)=Rs(2,l)++dimag(&
							abs(sum(gm(:,l)*dconjg(Uk(:,1))*Uk(:,4)))**2*(fk(1)-fk(4))*1d0/(omg+ek(4)-ek(1)+id)+&
							abs(sum(gm(:,l)*dconjg(Uk(:,4))*Uk(:,1)))**2*(fk(4)-fk(1))*1d0/(omg+ek(1)-ek(4)+id)+&
							abs(sum(gm(:,l)*dconjg(Uk(:,2))*Uk(:,3)))**2*(fk(2)-fk(3))*1d0/(omg+ek(3)-ek(2)+id)+&
							abs(sum(gm(:,l)*dconjg(Uk(:,3))*Uk(:,2)))**2*(fk(3)-fk(2))*1d0/(omg+ek(2)-ek(3)+id)&
							)
						Rs(3,l)=Rs(3,l)++dimag(&
							abs(sum(gm(:,l)*dconjg(Uk(:,2))*Uk(:,4)))**2*(fk(2)-fk(4))*1d0/(omg+ek(4)-ek(2)+id)+&
							abs(sum(gm(:,l)*dconjg(Uk(:,4))*Uk(:,2)))**2*(fk(4)-fk(2))*1d0/(omg+ek(2)-ek(4)+id)&
							)
						Rs(4,l)=Rs(4,l)+dimag(&
							abs(sum(gm(:,l)*dconjg(Uk(:,2))*Uk(:,1)))**2*(fk(2)-fk(1))*1d0/(omg+ek(1)-ek(2)+id)+&
							abs(sum(gm(:,l)*dconjg(Uk(:,1))*Uk(:,2)))**2*(fk(1)-fk(2))*1d0/(omg+ek(2)-ek(1)+id)+&
							abs(sum(gm(:,l)*dconjg(Uk(:,3))*Uk(:,4)))**2*(fk(3)-fk(4))*1d0/(omg+ek(4)-ek(3)+id)+&
							abs(sum(gm(:,l)*dconjg(Uk(:,4))*Uk(:,3)))**2*(fk(4)-fk(3))*1d0/(omg+ek(3)-ek(4)+id)&
							)
					enddo
				enddo
			enddo
			!$OMP END PARALLEL DO
			Rs=Rs/mr**2
			Rs(4,1)=sum(Rs(:,1))
			Rs(4,2)=sum(Rs(:,2))
			Rs(1,1)=Rs(1,1)+Rs(3,1)
			!R_rpa=-dimag(R(1)/(1d0+DJp*R(1)))
			!R_rpa=-dimag(R(1))/((1d0+DJp*dreal(R(1)))**2+(DJp*dimag(R(1)))**2)
			!write(50,"(e16.8$)")omg,-dimag(R),R_rpa,Djp,dreal(R(1)),(1d0+DJp*dreal(R(1)))**2+(DJp*dimag(R(1)))**2
			!write(50,"(e16.5$)")omg,R
			write(50,"(e16.5$)")nf,Tk,pg*Vd*4d0,sc*Vs*4d0,omg,Rs(:,1),Rs(:,2)
			!do i=1,2
				!call find_peak(pR(:,i),R(i),sg)
				!write(50,"(i3$)")sg
			!enddo
			do j=1,2
				do i=1,4
					call find_peak(pRs(:,i,j),Rs(i,j),sg)
					write(50,"(i3$)")sg
					if(i==4.and.sg==1) then
						pk(j)=omg
					endif
				enddo
			enddo
			!call find_peak(pR(:,1),R(1),sg)
			write(50,"(1X)")
			!if(peak(3,2)<R_rpa) then
			!peak(3,:)=(/omg,R_rpa/)
			!endif
			omg=omg+domg
		enddo
		write(50,"(1X)")
	end subroutine
	subroutine fermisurface(nf,pg,sc,ap,cp,omg)
		integer :: i,j,l,m
		complex(8) :: Uk(4,4)
		real(8) :: k(2),ek(4),cp,pg,sc,ap,omg,A,nf
		open(unit=20,file="../data/fermi.dat")
		m=256
		do i=0,m
			do j=0,m
				k=(/pi/m*i,pi/m*j/)
				call EU1(k,nf,pg,sc,ap,cp,ek,Uk)
				A=0d0
				do l=1,4
					A=A-dimag(abs(Uk(1,l))**2/(omg+img*0.01d0-ek(l)))
				enddo
				write(20,"(e16.8$)")A
			enddo
			write(20,"(1X)")
		enddo
		write(20,"(1X/)")
	end subroutine
	subroutine DOS(nf,pg,sc,ap,cp,Tk,omgr)
		integer :: i,j,sg
		complex(8) :: Uk(4,4),id
		real(8) :: k(2),cp,bt,sc,pg,ek(4),fk(4),A(4),pA(2,4),omg,omgr(:),domg,Tk,ap,nf
		open(unit=90,file="../data/DOS.dat")
		bt=1d0/Tk
		id=(0d0,0.003d0)
		omg=omgr(1)
		domg=omgr(3)
		pA=0d0
		do while(omg<omgr(2))
			A=0d0
			omg=omg+domg
			!$OMP PARALLEL DO REDUCTION(+:A) PRIVATE(k,ek,Uk,fk)
			do i=1,mk
				do j=1,mk
					k=(/pi/mk*i,pi/mk*j/)
					call EU2(k,nf,pg,sc,ap,cp,ek,Uk)
					!fk=1d0/(1d0+exp(bt*ek))
					A=A-DIMAG(Uk(1,:)*dconjg(Uk(1,:))/(omg-ek+id))
				enddo
			enddo
			!$OMP END PARALLEL DO
			A=A/(mk**2)
			write(90,"(e16.8$)")omg,A,sum(A)
			do i=1,4
				call find_peak(pA(:,i),A(i),sg)
				write(90,"(i3$)")sg
			enddo
			write(90,"(1X)")
			!write(*,"(e16.8$)")omg,A
			!write(*,"(1X)")
		enddo
		write(90,"(1X)")
	end subroutine
	subroutine EDC(k,pg,sc,ap,cp,nf,Tk,omgr)
		integer :: sg
		complex(8) :: Uk(4,4)
		real(8) :: k(2),cp,bt,sc,pg,ap,Tk,ek(4),fk(4),A,omg,omgr(:),domg,pA(2),nf
		pA=0d0
		bt=1d0/Tk
		omg=omgr(1)
		domg=omgr(3)
		!call EU1(k,nf,pg,sc,ap,cp,ek,Uk)
		call EU2(k,nf,pg,sc,ap,cp,ek,Uk)
		fk=1d0/(1d0+exp(bt*ek))
		do while(omg<omgr(2))
			omg=omg+domg
			!A=-sum(DIMAG(fk*Uk(1,:)*dconjg(Uk(1,:))/(omg-ek+img*domg*400d0)))
			A=-sum(dimag(Uk(1,:)*dconjg(Uk(1,:))/(omg-ek+img*domg*40d0)))
			call find_peak(pA,A,sg)
			write(30,"(e16.8$)")omg,A
			write(30,"(i2$)")sg
			write(30,"(1X)")
		enddo
	end subroutine
	subroutine band(pg,sc,ap,cp,nf,Tk,ki,kf)
		integer :: i,j,l,bd(1),sg,m
		complex(8) :: Uk(4,4)
		real(8) :: ki(:),kf(:),km(2),kn(2),ek(4),pek(2,4),cp,pg,sc,gap,Tk,peak(2),pk(4),eks,eka,ap,nf
		open(unit=10,file="../data/band.dat")
		pek=0d0
		m=256
		!j=0
		do j=1,size(ki),2
			kn=ki(j:j+1)
			do while(all((kn-kf(j:j+1))*(kf(j:j+1)-ki(j:j+1))<=0d0))
				kn=kn+(kf(j:j+1)-ki(j:j+1))/m
				call EU2(kn,nf,pg,sc,ap,cp,ek,Uk)
				!write(10,"(i4$)")j
				write(10,"(e16.8$)")nf,Tk,kn,(ek(l),abs(Uk(1,l)),l=1,4)
				do i=1,4
					call find_peak(pek(:,i),abs(ek(i)),sg)
					write(10,"(i3$)")sg
				enddo
				write(10,"(e16.8$)")&
					!-0.5d0*(cos(kn(1))-cos(kn(2)))*sc*Vs*4d0,-0.5d0*(cos(kn(1))-cos(kn(2)))*pg*Vd*4d0,kn
				-0.5d0*(cos(kn(1))-cos(kn(2)))*sc*Vs*4d0,pg*Vd*4d0
				write(10,"(1X)")
				!j=j+1
			enddo
		enddo
		write(10,"(1X)")
	end subroutine
	subroutine phasediagram(nr)
		real(8) :: nr(:),n,Tc
		open(unit=80,file="../data/phase_tJ.dat")
		write(*,*)"*******phase diagram start********"
		n=nr(1)
		do
			write(80,"(e16.8$)")n
			Tc=5d-2
			call findTc(n,Tc,0)
			write(80,"(e16.8$)")Tc
			call findTc(n,Tc,1)
			write(80,"(e16.8$)")Tc
			if(Tc/=0d0) then
				call findTc(n,Tc,2)
				write(80,"(e16.8$)")Tc
			endif
			write(80,"(1X)")
			if(n>nr(2)) then
				exit
			endif
			n=n+nr(3)
		enddo
		write(*,*)"*******phase diagram finish********"
	end subroutine
	subroutine findTc(nf,Tc,sg)
		real(8) :: pg,sc,ap,cp,Tk,dTk,psc,ppg,Tc,nf,zo,mTk
		integer :: sg,isg
		write(*,"(x)")
		write(*,*)"*******findTc********"
		mTk=2d-4
		zo=1d-3
		Tc=max(Tc,0.01d0)
		dTk=Tc*0.2d0
		Tk=Tc*0.5d0
		psc=0d0
		ppg=0d0
		do
			call selforder(pg,sc,ap,cp,nf,Tk,1)
			write(*,"(e12.4$)")sc,pg,Tk
			if(sg==0) then
				call find_cross(psc,sc-zo,isg)
				if(isg/=0) then
					if(dTk<mTk) then
						Tc=Tk
						exit
					endif
					dTk=dTk*0.5d0
				endif
				Tk=Tk+dTk*sign(1d0,sc-zo)
			elseif(sg==1) then
				call find_cross(ppg,pg-zo,isg)
				if(isg/=0) then
					if(dTk<mTk) then
						Tc=Tk
						exit
					endif
					dTk=dTk*0.5d0
				endif
				Tk=Tk+dTk*sign(1d0,pg-zo)
			else
				call find_cross(ppg,pg-zo,isg)
				if(isg/=0) then
					if(dTk<mTk) then
						Tc=Tk
						exit
					endif
					dTk=dTk*0.5d0
				endif
				Tk=Tk+dTk*sign(1d0,pg-zo)
			endif
			if(Tk<0d0) then
				if(dTk>mTk) then
					dTk=mTk
					Tk=4d0*mTk
					cycle
				else
					Tc=0d0
					exit
				endif
			endif
			write(*,"(e12.4)")dTk
		enddo
	end subroutine
	subroutine rpainstable(pg,sc,ap,cp,nf,Tk,den)
		integer :: i,j,l,m,n,mr
		complex(8) :: Uk(4,4),Ukq(4,4)
		real(8) :: k(2),q(2),cp,bt,pg,sc,ek(4),ekq(4),fk(4),fkq(4),Tk,DJq,Xq(2),den(2),Vrpa(2),tran(2),ap,nf
		mr=512
		bt=1d0/Tk
		Xq=0d0
		select case(pgflag)
		case("ddw")
			Vrpa=(/Vd,Vs/)
		case("sdw")
			Vrpa=(/DJ,Vs/)
		end select
		!$OMP PARALLEL DO REDUCTION(+:Xq) PRIVATE(k,ek,Uk,fk,ekq,Ukq,fkq,tran)
		do i=0,mr
			k(1)=pi/mr*i
			do j=0,mr
				k(2)=pi/mr*j
				call EU1(k,nf,pg,sc,ap,cp,ek,Uk)
				call EU1(k+(/pi,pi/),nf,pg,sc,ap,cp,ekq,Ukq)
				fk=1d0/(1d0+exp(bt*ek))
				fkq=1d0/(1d0+exp(bt*ekq))
				do n=1,4
					do m=1,4
						!if(ek(n)*ekq(m)<0d0) then
						!cycle
						!endif
						select case(pgflag)
						case("ddw")
							tran(1)=(Uk(1,n)*dconjg(Uk(1,n))*Ukq(1,m)*dconjg(Ukq(1,m))+&
								Uk(3,n)*dconjg(Uk(1,n))*Ukq(1,m)*dconjg(Ukq(3,m)))*(cos(k(1))-cos(k(2)))**2
						case("sdw")
							tran(1)=(Uk(1,n)*dconjg(Uk(1,n))*Ukq(1,m)*dconjg(Ukq(1,m))+&
								Uk(3,n)*dconjg(Uk(1,n))*Ukq(1,m)*dconjg(Ukq(3,m)))*2d0
						end select
						tran(2)=(Uk(1,n)*dconjg(Uk(1,n))*Uk(3,m)*dconjg(Uk(3,m))-&
							Uk(2,n)*dconjg(Uk(1,n))*Uk(3,m)*dconjg(Uk(4,m)))*(cos(k(1))-cos(k(2)))**2
						Xq(1)=Xq(1)+tran(1)*(fk(n)-fkq(m))/(ek(n)-ekq(m)+1d-10)
						Xq(2)=Xq(2)+tran(2)*(fk(n)-fk(m))/(ek(n)-ek(m)+1d-10)
					enddo
				enddo
			enddo
		enddo
		!$OMP END PARALLEL DO
		Xq=Xq/mr**2
		den=1d0+Vrpa*Xq
	end subroutine
end module
program main
	use pmt
	use M_utility
	use selfcons
	use phys
	implicit none
	complex(8) :: Uk(4,4)
	real(8) :: kf(2),cp,sc,pg,ap,gap,Tk,pk0(2),pk(2),&
		!nd(3)=(/0.10d0,0.13d0,0.15d0/),&
		nd(1)=(/0.135d0/),&
		!nd(12)=(/0.05d0,0.07d0,0.09d0,0.1d0,0.12d0,0.125d0,0.13d0,0.135d0,0.145d0,0.15d0,0.18d0,0.25d0/),&
		!nd(8)=(/0.1d0,0.125d0,0.13d0,0.135d0,0.145d0,0.15d0,0.18d0,0.2d0/),&
		!nd(9)=(/0.08d0,0.1d0,0.12d0,0.125d0,0.13d0,0.135d0,0.145d0,0.15d0,0.18d0/),&
		!Td(14)=(/0.0d0,0.2d0,0.4d0,0.45d0,0.5d0,0.55d0,0.6d0,0.65d0,0.7d0,0.75d0,0.8d0,0.85d0,0.9d0,1d0/),&
		Td(16),&
			nf,Tc,Tp
	integer :: i,j,l
	logical :: f
	if(pgflag=="ddw".or.pgflag=="sdw") then
		write(*,*)"calculate for: ",pgflag
	else
		write(*,*)"the pg must be ddw or sdw"
		stop
	endif
	f=openfile(30,"../data/EDC.dat")
	f=openfile(40,"../data/gaptemp.dat")
	f=openfile(50,"../data/sraman.dat")
	f=openfile(60,"../data/gaptheta.dat")
	f=openfile(70,"../data/ssum.dat")
	f=openfile(80,"../data/superfluid.dat")
	cp=0d0
	Tc=5d-2
	Tp=5d-2
	!write(*,"(e12.4$)")Tk,superfluid(1d0,1d0,0d0,-0.7d0,0.8d0,0.002d0)
	!stop
	!call phasediagram((/0.9d0,0.83d0,-0.01d0/))
	do i=1,size(nd)
		nf=1d0-nd(i)
		write(*,"(e12.4$)")nf
		!pg=0.049d0*(1d0-(1d0-nf)/0.2d0)
		!sc=sqrt(0.037d0**2-pg**2)
		!cp=-(1-nf)**2
		!write(30,"(e12.4$)")nd(i),superfluid(pg,sc,0d0,cp,0d0,1d-4)
		!write(30,"(x)")
		!cycle
		!call findTc(nf,Tc,0)
		write(*,*)"*******findTc finish********"
		write(*,"(1x,e12.4$)")Tc
		!call findTc(nf,Tp,1)
		!write(*,"(e12.4$)")Tp
		!call selforder(pg,sc,ap,cp,nf,Tp,1)
		!call band(pg,sc,ap,cp,nf,Tp,(/pi,pi/2d0/),(/pi,0d0/))
		!stop
		do j=0,size(Td)-1
			if(Tc<1d-10) then
				cycle
			endif
			!if(j/=0.and.j/=6.and.j/=9) then
				!cycle
			!endif
			!Tk=max(j*max(Tc,Tp)/(size(Td)-1),1d-5)
			Tk=max(j*max(0d0,Tc)/(size(Td)-1),1d-5)
			!Tk=max(j*Tc/(size(Td)-1),1d-5)
			!Tk=max(Tc*Td(j+1),1d-5)
			!Tk=max(j*Tp/(size(Td)-1),1d-5)
			!write(*,"(e12.4$)")Tk/Tc
			!write(70,"(e12.4$)")nd(i),Tk/Tc
			write(*,"(e12.4$)")Tk
			write(70,"(e12.4$)")nd(i),Tk
			!sc=0d0
			call selforder(pg,sc,ap,cp,nf,Tk,1)
			call EDC((/pi,0.52768939d0/),pg,sc,ap,cp,nf,Tk,(/-0.5d0,0.7d0,0.0001d0/))
			!call selfconsist_tg(nf,Tk,sc,pg,ap,cp)
			if(pg<1d-3) then
				pg=1d-10
			endif
			write(*,"(e12.4$)")pg,sc,cp,ap
			!write(40,"(e12.4$)")nf,Tk/Tc,pg,sc,cp
			write(40,"(e12.4$)")Tk,pg*Vd*4d0,sc*Vs*4d0,pg,sc,cp,ap
			!call raman(pg,sc,ap,cp,nf,Tk,(/0d0,0.4d0,0.0002d0/),pk)
			!call raman(pg,sc,ap,cp,nf,Tk,(/0d0,0.4d0,0.0002d0/),pk)
			!write(30,"(e12.4$)")Tk,superfluid(pg,sc,ap,cp,nf,Tk)
			!write(30,"(e12.4$)")Tk,superfluid(pg,sc,ap,cp,nf,Tk),superfluid(0d0,sc,ap,cp,nf,Tk),superfluid(pg,sc,ap,0d0,nf,Tk),superfluid(0d0,sc,ap,0d0,nf,Tk)
			!write(30,"(e12.4$)")Tk,superfluid(pg,sc,0d0,cp,0d0,Tk)
			!if(j==0) then
				!pk0=pk
			!endif
			!write(70,"(e12.4$)")pk/pk0
			!sc=0d0
			!call band(pg,sc,ap,cp,nf,Tk,(/pi,pi/2d0,pi,0d0,pi,0d0/),(/pi,0d0,pi/2d0,pi/2d0,pi/2d0,0d0/))
			!call band(pg,sc,ap,cp,nf,Tk,(/pi,pi/2d0/),(/pi,0d0/))
			!call DOS(nf,pg,sc,ap,cp,Tk,(/-0.4d0,0.4d0,0.0002d0/))
			!write(*,"(e12.4$)")nf,Tk,Tk/Tc,gap
			!do l=0,10
				!call EDC((/pi,pi/20d0*l/),pg,sc,ap,cp,Tk,(/-0.5d0,0.7d0,0.0001d0/))
			!enddo
			write(40,"(1X)")
			write(70,"(1X)")
			write(*,"(1X)")
			write(30,"(1X)")
			write(80,"(1X)")
		enddo
		write(10,"(1X)")
		write(30,"(1X)")
		write(50,"(1X/)")
		write(40,"(1X)")
		write(70,"(1X/)")
		write(80,"(1X)")
	enddo
end program
