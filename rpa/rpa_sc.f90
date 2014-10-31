module pmt
	use M_const
	implicit none
	real(8), parameter :: t(5)=(/2d0,-0.9d0,0.d0,0.d0,0d0/),&
		!V=0.12d0,DJ=0.35d0,&
		V=0.d0,DJ=1.5d0,&
	Vs=0.5d0*(DJ-V),Vd=0.5d0*DJ+V,alp=0.34
	real(8) :: it=0.02d0
	character(3) :: pgflag="ddw"
	!character(3) :: pgflag="sdw"
	!real(8), parameter :: t(5)=(/0.1483d0,-0.035d0,0.0355d0,0d0,0d0/),&
		!!V=0.12d0,DJ=0.35d0,&
		!V=0.d0,DJ=1.2d0,&!DJ=1.5d0,&
	!Vs=0.25,Vd=0.5d0*DJ+V,alp=0.2
end module
module selfcons
	use pmt
	use M_utility
	implicit none
contains
	subroutine selfconsist_tg(nf,Tk,pg,sc,ap,sp)
		real(8) :: Tk,nf
		complex(8) :: Uk(4,4)
		real(8),target :: k(2),np,pnp,pg,pgp,ek(4),sp,sc,scp,cvg1,al,ap,app,cvg
		integer :: i,j,c,sg,mk
		real(8),pointer :: p
		mk=512
		cvg=1e-6
		if(Tk<0d0) then
			write(*,*)"Tk<0"
			stop
		endif
		al=1d0
		al=0.2d0
		pg=max(pg,0.1d0)
		sc=max(sc,0.1d0)
		ap=max(ap,0.1d0)
		c=0
		do 
			pnp=0d0
			do 			
				pgp=0d0
				scp=0d0
				app=0d0
				np=0d0
				c=c+1
				pg=0d0
				!$OMP PARALLEL DO REDUCTION(+:np,pgp,scp,app) PRIVATE(k,ek,Uk) SCHEDULE(GUIDED)
				do i=0,mk
					do j=0,min(i,mk-i)
						k=(/pi/mk*i,pi/mk*j/)
						call EU(k,nf,pg,sc,ap,sp,ek,Uk)
						call order(k,ek,Uk,Tk,np,pgp,scp,app)
					enddo
				enddo
				!$OMP END PARALLEL DO
				np=np/(mk**2)*2d0
				pgp=pgp/(mk**2)*2d0
				scp=scp/(mk**2)*2d0
				app=app/(mk**2)*2d0
				if(abs(nf-np)<1d-5) then
					exit
				endif
				call find_cross(pnp,np-nf,sg)
				if(sg/=0) then
					al=al*0.3d0
				endif
				sp=sp+al*sign(1d0,nf-np)
			enddo
			if((abs(scp-sc)+abs(pgp-pg)+abs(app-ap))<cvg) then
				exit
			endif
			sc=(scp*0.4d0+sc*0.6d0)
			pg=(pgp*0.4d0+pg*0.6d0)
			ap=(app*0.4d0+ap*0.6d0)
		enddo
	end subroutine
	subroutine EU(k,nf,pg,sc,ap,sp,ek,Uk)
		complex(8) :: Uk(4,4),cth,sth,cfy(2),sfy(2)
		real(8) :: eka,eks,e1(2),e2(2),ek(4),k(2),pg,pgk,sp,sck,sc,gk,th,fy(2),cos2th,sin2th,cos2fy(2),sin2fy(2),ap,dp,nf
		dp=abs(1d0-nf)
		eks=-4d0*dp*t(2)*cos(k(1))*cos(k(2))-sp-2d0*dp*t(3)*(cos(2d0*k(1))+cos(2d0*k(2)))
		eka=-2d0*(dp*t(1)+ap*Vd)*(cos(k(1))+cos(k(2)))
		gk=0.5d0*(cos(k(1))-cos(k(2)))
		select case(pgflag)
		case("ddw")
			pgk=pg*gk*4d0*Vd
		case("sdw")
			pgk=-pg*DJ*2d0
		end select
		sck=sc*gk*Vs*4d0
		e1=eks+(/1d0,-1d0/)*sqrt(pgk**2+eka**2)
		e2=sqrt(e1**2+sck**2)
		!ek=(/e2*sign(1d0,e1),-e2*sign(1d0,e1)/)
		ek=(/e2,-e2/)
		cos2th=eka/sqrt(pgk**2+eka**2)
		sin2th=pgk/sqrt(pgk**2+eka**2)
		if(abs(eka)<1d-8) then
			cos2th=0d0
			sin2th=1d0
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
end module
module phys
	use pmt
	use selfcons
	use M_utility
	implicit none
	integer, parameter :: mk=64
contains
	subroutine spinsus(nf,Tk,pg,sc,ap,sp,q,omg,Xq)
		implicit none
		integer :: i,j,l,m,n,mq(2),ki,kj,qi,qj
		complex(8) :: Uk(4,4),Ukq(4,4),omg(:),Xq(:,:,:,:,:),Xtmp1(size(Xq,1),2,2),Xtmp2(size(Xq,1),2,2),Det(size(Xq,1))
		real(8) :: k(2),q(:,:,:),sp,bt,pg,sc,ap,ek(4),ekq(4),fk(4),fkq(4),DJq,DJqq,nf,Tk
		bt=1d0/Tk
		Xq=0d0
		mq(1)=size(Xq,2)
		mq(2)=size(Xq,3)
		write(*,"(A$)")"    "
		!$OMP PARALLEL DO REDUCTION(+:Xq) PRIVATE(k,ek,Uk,fk,ekq,Ukq,fkq) SCHEDULE(GUIDED)
		do ki=0,mk-1
			write(*,"(4A,I3,'%'$)")char(8),char(8),char(8),char(8),int(real(ki)/mk*100)
			do kj=0,mk-1
				k=2d0*pi*(/ki,kj/)/mk
				call EU(k,nf,pg,sc,ap,sp,ek,Uk)
				fk=1d0/(1d0+exp(bt*ek))
				do qi=1,mq(1)
					do qj=1,mq(2)
						call EU(k+q(qi,qj,:),nf,pg,sc,ap,sp,ekq,Ukq)
						fkq=1d0/(1d0+exp(bt*ekq))
						do n=1,4
							do m=1,4
								do i=1,2
									do j=1,2
										Xq(:,qi,qj,i,j)=Xq(:,qi,qj,i,j)+(Ukq(1,n)*Uk(3+i-1,m)+Ukq(2,n)*Uk(3+i-1+(-1)**(i-1),m)-&
											Ukq(3+i-1,n)*Uk(1,m)-Ukq(3+i-1-(-1)**(i-1),n)*Uk(2,m))*&
											dconjg(Ukq(1,n)*Uk(3+j-1,m))*(1-fk(n)-fkq(m))/(omg+ek(n)+ekq(m)+1d-10)
									enddo
								enddo
							enddo
						enddo
					enddo
				enddo
			enddo
		enddo
		!$OMP END PARALLEL DO
		Xq=Xq/mk**2
		!!$OMP PARALLEL DO PRIVATE(DJq,DJqq,Xtmp2,Xtmp1,Det) SCHEDULE(GUIDED)
		!do qi=1,mq(1)
			!do qj=1,mq(2)
				!DJq=alp*(cos(q(qi,qj,1))+cos(q(qi,qj,2)))
				!DJqq=-alp*(cos(q(qi,qj,1))+cos(q(qi,qj,2)))
				!Xtmp1=Xq(:,qi,qj,:,:)
				!Xtmp2(:,1,1)=1d0+DJqq*Xtmp1(:,2,2)
				!Xtmp2(:,2,2)=1d0+DJq*Xtmp1(:,1,1)
				!Xtmp2(:,1,2)=-DJq*Xtmp1(:,1,2)
				!Xtmp2(:,2,1)=-DJqq*Xtmp1(:,2,1)
				!Det=1d0/(Xtmp2(:,1,1)*Xtmp2(:,2,2)-Xtmp2(:,1,2)*Xtmp2(:,2,1))
				!Xq(:,qi,qj,1,1)=(Xtmp1(:,1,1)*Xtmp2(:,1,1)+Xtmp1(:,1,2)*Xtmp2(:,2,1))*Det
				!Xq(:,qi,qj,1,2)=(Xtmp1(:,1,1)*Xtmp2(:,1,2)+Xtmp1(:,1,2)*Xtmp2(:,2,2))*Det
				!Xq(:,qi,qj,2,1)=(Xtmp1(:,2,1)*Xtmp2(:,1,1)+Xtmp1(:,2,2)*Xtmp2(:,2,1))*Det
				!Xq(:,qi,qj,2,2)=(Xtmp1(:,2,1)*Xtmp2(:,1,2)+Xtmp1(:,2,2)*Xtmp2(:,2,2))*Det
			!enddo
		!enddo
		!!$OMP END PARALLEL DO
	end subroutine
	subroutine greenfunc(nf,pg,sc,ap,sp,k,omg,G0,G0_inv)
		implicit none
		integer :: i,j,l,m,n,o,p,ii,jj,sig
		complex(8) :: Uk(4,4),omg(:),G0(:,:,:,:,:),G0_inv(:,:,:,:,:)
		real(8) :: pg,sc,ap,sp,ek(4),k(:,:,:),nf
		write(*,*)"*************greenfunc() start***********"
		G0=0d0
		G0_inv=0d0
		l=size(G0,4)
		o=size(G0,2)
		!$OMP PARALLEL DO PRIVATE(ek,Uk) SCHEDULE(GUIDED)
		do i=1,size(k,1)
			do j=1,size(k,2)
				call EU(k(i,j,:),nf,pg,sc,ap,sp,ek,Uk)
				do m=1,l
					do n=1,l
						do p=1,l
							G0(:,i,j,m,n)=G0(:,i,j,m,n)+Uk(m,p)*dconjg(Uk(n,p))/(omg-ek(p)+1d-10)
							G0_inv(:,i,j,m,n)=G0_inv(:,i,j,m,n)+Uk(m,p)*dconjg(Uk(n,p))*(omg-ek(p))
						enddo
					enddo
				enddo
			enddo
		enddo
		!$OMP END PARALLEL DO
	end subroutine
	subroutine effint(nf,Tk,pg,sc,ap,sp,omg,Veff)
		use M_fft, only : fft3d
		implicit none
		integer :: i,j,n,ii,jj,mq
		complex(8) :: Veff(:,:,:,:,:),omg(:)
		real(8) :: pg,sc,ap,sp,DJq,q(size(Veff,2),size(Veff,3),2),nf,Tk
		write(*,*)"*************effint() start***********"
		!write(*,"(A$)")"    "
		!write(*,"(4A,I3,'%'$)")char(8),char(8),char(8),char(8),int(2d0*i/n*100)
		mq=size(Veff,2)
		!$OMP PARALLEL DO SCHEDULE(GUIDED)
		do i=1,mq
			do j=1,mq
				q(i,j,:)=2d0*pi*(/i-1,j-1/)/mq
			enddo
		enddo
		!$OMP END PARALLEL DO
		call spinsus(nf,Tk,pg,sc,ap,sp,q,omg,Veff)
		!!$OMP PARALLEL DO PRIVATE(DJq) SCHEDULE(GUIDED)
		!do i=1,mq
			!do j=1,mq
				!DJq=(cos(q(i,j,1))+cos(q(i,j,2)))
				!Veff(:,i,j,1,1)=Veff(:,i,j,1,1)*DJq**2
				!Veff(:,i,j,1,2)=-Veff(:,i,j,1,2)*DJq**2
				!Veff(:,i,j,2,1)=-Veff(:,i,j,2,1)*DJq**2
				!Veff(:,i,j,2,2)=Veff(:,i,j,2,2)*DJq**2
			!enddo
		!enddo
		!!$OMP END PARALLEL DO
		!do i=1,2
			!do j=1,2
				!call symfill(Veff(:,:,:,i,j),1)
			!enddo
		!enddo
	end subroutine
	subroutine symfill(A,sig)
		implicit none
		integer :: i,j,n,ii,jj,sig
		complex(8) :: A(:,:,:)
		n=size(A,2)
		!$OMP PARALLEL DO PRIVATE(ii,jj) SCHEDULE(GUIDED)
		do i=1,n/2
			do j=0,i
				ii=mod(n-i,n)
				jj=mod(n-j,n)
				A(: , j  , i ) = A(: , i , j) * sig
				A(: , jj , i ) = A(: , i , j) * sig
				A(: , j  , ii) = A(: , i , j) * sig
				A(: , jj , ii) = A(: , i , j) * sig
				A(: , ii , j ) = A(: , i , j)
				A(: , i  , jj) = A(: , i , j)
				A(: , ii , jj) = A(: , i , j)
			enddo
		enddo
		!$OMP END PARALLEL DO
	end subroutine
	subroutine selfefft(Tk,G0,Veff,sfeg)
		use M_fft
		implicit none
		integer :: i,j,n,m
		complex(8) :: sfeg(:,:,:,:,:),Veff(:,:,:,:,:),G0(:,:,:,:,:),&
			Gtmp(size(G0,1),size(G0,2),size(G0,3),size(G0,4),size(G0,5))
		real(8) :: k(2),q(2),sp,bt,sc,ap,ek(2),ekq(2),fk(2),fkq(2),DJq,Tk
		write(*,*)"*************selfefft() start***********"
		Gtmp=G0
		do i=1,2
			do j=1,2
				call fft3d(Veff(:,:,:,i,j),-1)
			enddo
		enddo
		do i=1,4
			do j=1,4
				call fft3d(Gtmp(:,:,:,i,j),-1)
			enddo
		enddo
		do i=1,4
			do j=1,4
				do n=1,1
					do m=1,2
						sfeg(:,:,:,i,j)=sfeg(:,:,:,i,j)+Veff(:,:,:,n,m)*Gtmp(:,:,:,i+(1-n)*(-1)**i,j+(1-m)*(-1)**j)
					enddo
				enddo
				call fft3d(sfeg(:,:,:,i,j),1)
			enddo
		enddo
		sfeg=sfeg*Tk*size(sfeg,1)
	end subroutine
	!subroutine rpainstable(sc,ap,sp,q)
	!implicit none
	!integer :: i,j,l,m,n
	!complex(8) :: Uk(2,2),Ukq(2,2)
	!real(8) :: k(2),q(2),sp,bt,pg,sc,ap,ek(2),ekq(2),fk(2),fkq(2),omg,omgr(2),domg,DJq,Xq,Xq_rpa
	!bt=1d0/Tk
	!Xq=0d0
	!!$OMP PARALLEL DO REDUCTION(+:Xq) PRIVATE(k,ek,Uk,fk,ekq,Ukq,fkq) SCHEDULE(GUIDED)
	!do i=0,mq
	!k(1)=pi/mq*i
	!do j=0,mq
	!k(2)=pi/mq*j
	!call EU(k,pg,sc,ap,sp,ek,Uk)
	!call EU(k+q,pg,sc,ap,sp,ekq,Ukq)
	!fk=1d0/(1d0+exp(bt*ek))
	!fkq=1d0/(1d0+exp(bt*ekq))
	!do n=1,2
	!do m=1,2
	!Xq=Xq+(Uk(1,n)*Ukq(2,m)*dconjg(Uk(1,n)*Ukq(2,m))-Uk(2,n)*Ukq(1,m)*dconjg(Uk(1,n)*Ukq(2,m)))*&
	!(1-fk(n)-fkq(m))/(ek(n)+ekq(m))
	!enddo
	!enddo
	!enddo
	!enddo
	!!$OMP END PARALLEL DO
	!Xq=Xq/mq**2
	!write(*,"(2e16.3)")1d0-nf,1d0+DJq*Xq
	!end subroutine
	subroutine band(pg,sc,ap,sp,nf,Tk,ki,kf)
		!call band(pg,sc,ap,sp,nf,Tk,(/0d0,0d0,pi,pi,pi,0d0/),(/pi,pi,pi,0d0,0d0,0d0/))
		integer :: i,j,l,bd(1),sg,m
		complex(8) :: Uk(4,4)
		real(8) :: ki(:),kf(:),km(2),kn(2),ek(4),pek(2,4),sp,pg,sc,gap,Tk,peak(2),pk(4),eks,eka,ap,nf
		open(unit=10,file="../data/band.dat")
		pek=0d0
		m=256
		!j=0
		do j=1,size(ki),2
			kn=ki(j:j+1)
			do while(all((kn-kf(j:j+1))*(kf(j:j+1)-ki(j:j+1))<=0d0))
				kn=kn+(kf(j:j+1)-ki(j:j+1))/m
				call EU(kn,nf,pg,sc,ap,sp,ek,Uk)
				!write(10,"(i4$)")j
				write(10,"(e16.8$)")(ek(l),abs(Uk(1,l)),l=1,4)
				do i=1,4
					call find_peak(pek(:,i),abs(ek(i)),sg)
					write(10,"(i3$)")sg
				enddo
				write(10,"(e16.8$)")&
					-0.5d0*(cos(kn(1))-cos(kn(2)))*sc*Vs*4d0,-0.5d0*(cos(kn(1))-cos(kn(2)))*pg*Vd*4d0
				!-0.5d0*(cos(kn(1))-cos(kn(2)))*sc*Vs*4d0,pg*Vd*4d0
				write(10,"(1X)")
				!j=j+1
			enddo
		enddo
		write(10,"(1X)")
	end subroutine
end module
program main
	use M_green
	use M_utility
	use selfcons
	use phys
	use M_matrix
	use M_fft
	implicit none
	integer, parameter :: mq=mk,mo=128,mro=256
	complex(8) :: iomgf(mo),iomgb(mo),romg(mro),&
		G0(mo,mq,mq,4,4),G0_inv(mo,mq,mq,4,4),Veff(mo,mq,mq,2,2),sfeg(mo,mq,mq,4,4),Gf(mo,mq,mq,4,4),&
		rG0(mro,mq,mq,4,4),rVeff(mro,mq,mq,2,2),rsfeg(mro,mq,mq,4,4),rGf(mro,mq,mq,4,4)
	real(8) :: k(mk,mk,2),q(mq,mq,2),sp,sc,pg,ap,&
		Tk=5d-3,nf,&
		nd(1)=(/0.12d0/)
	integer :: i,j,l,m,n,ik(2)
	open(unit=10,file="../data/1d.dat")
	open(unit=20,file="../data/2d.dat")
	open(unit=30,file="../data/EDC.dat")
	open(unit=40,file="../data/pattern.dat")
	!iomgf=arth(-mo/2,mo)*cmplx(0d0,2d0*pi*Tk)+cmplx(0d0,pi*Tk)
	!iomgb=arth(-mo/2,mo)*cmplx(0d0,2d0*pi*Tk)
	iomgf(1:mo/2)=arth(0,mo/2)*cmplx(0d0,2d0*pi*Tk)+cmplx(0d0,pi*Tk)
	iomgf(mo/2+1:mo)=arth(-mo/2,mo/2)*cmplx(0d0,2d0*pi*Tk)+cmplx(0d0,pi*Tk)
	iomgb(1:mo/2)=arth(0,mo/2)*cmplx(0d0,2d0*pi*Tk)
	iomgb(mo/2+1:mo)=arth(-mo/2,mo/2)*cmplx(0d0,2d0*pi*Tk)
	romg=arth(-1.2d0,4d0/mro,mro)+it*img
	!romg=(/0.1d0,0.3d0,0.515d0,0.8d0/)+it*img
	do i=1,mq
		do j=1,mq
			q(i,j,:)=2d0*pi*(/i-1,j-1/)/mq
		enddo
	enddo
	!q(1,1,:)=(/pi,pi/)
	do n=1,size(nd)
		nf=1d0-nd(n)
		call selfconsist_tg(nf,Tk,pg,sc,ap,sp)
		pg=1d-1
		sc=1d-7
		write(*,"(e14.4$)")nf,Tk,pg,sc,ap,sp
		write(*,"(1X)")
		!check spinsus()
		!call spinsus(nf,Tk,pg,sc,ap,sp,q(1:1,1:1,:),iomgb,Veff(:,1:1,1:1,:,:))
		!call pade(iomgb(mo/2:mo),Veff(mo/2:mo,1,1,1,1),romg,rVeff(:,1,1,1,1))
		!write(10,"(2e14.4)")(real(romg(i)),dimag(rVeff(i,1,1,1,1)),i=1,size(romg,1))
		!write(10,"(1x/)")
		!call spinsus(nf,Tk,pg,sc,ap,sp,q(1:1,1:1,:),romg,rVeff(:,1:1,1:1,:,:))
		!write(10,"(2e14.4)")(real(romg(i)),dimag(rVeff(i,1,1,1,1)),i=1,size(romg,1))
		!stop
		!write(10,"(2e16.4)")(real(romg(i)),dimag(rVeff(i,1,1,1,1)),i=1,size(romg,1))
		!call spinsus(nf,Tk,pg,sc,ap,sp,q,romg,rVeff)
		!call mwrite(20,dimag(rVeff(1,:,:,1,1)))
		!write(20,"(1x)")
		!call mwrite(20,dimag(rVeff(2,:,:,1,1)))
		!write(20,"(1x)")
		!call mwrite(20,dimag(rVeff(3,:,:,1,1)))
		!write(20,"(1x)")
		!stop
		!end check
		call greenfunc(nf,pg,sc,ap,sp,q,iomgf,G0,G0_inv)
		!call greenfunc(nf,pg,sc,ap,sp,q,romg,rGf,rG0)
		!do i=1,mq
			!do l=1,mq
				!call pade(iomgf(1:mo/2),G0(1:mo/2,i,l,1,1),romg,rG0(:,i,l,1,1))
			!enddo
		!enddo
		!do l=1,mq/2
			!write(20,"(e14.4$)")(dimag(rG0(i,l,l,1,1)),i=1,mro)
			!write(20,"(1x)")
		!enddo
		!write(20,"(1x)")
		!do l=mq/2,1,-1
			!write(20,"(e14.4$)")(dimag(rG0(i,mq/2,l,1,1)),i=1,mro)
			!write(20,"(1x)")
		!enddo
		!write(20,"(1x)")
		!do l=mq/2,1,-1
			!write(20,"(e14.4$)")(dimag(rG0(i,l,1,1,1)),i=1,mro)
			!write(20,"(1x)")
		!enddo
		!write(20,"(1x)")
		!stop
		!call pade(iomgf,G0(:,mq/2,int(0.15d0/2d0*mq),1,1),romg,rG0(:,1,1,1,1))
		!write(10,"(2e14.4)")(real(romg(i)),dimag(rG0(i,1,1,1,1)),i=1,size(romg,1))
		!write(10,"(1x/)")
		!call greenfunc(nf,pg,sc,ap,sp,q,romg,rG0,rGf)
		!write(10,"(2e14.4)")(real(romg(i)),dimag(rG0(i,mq/2,int(0.15d0/2d0*mq),1,1)),i=1,size(romg,1))
		!stop
		call effint(nf,Tk,pg,sc,ap,sp,iomgb,Veff)
		!call effint(nf,Tk,pg,sc,ap,sp,romg,rVeff)
		!do i=1,mq
			!do j=1,mq
				!call pade(iomgb(1:mo/2),Veff(1:mo/2,i,j,1,1),romg,rVeff(:,i,j,1,1))
			!enddo
		!enddo
		do l=1,mo,8
			call mwrite(20,dimag(Veff(l,:,:,1,1)))
			write(20,"(1x)")
			call mwrite(20,dimag(Veff(l,:,:,2,2)))
			write(20,"(1x)")
			call mwrite(20,dimag(Veff(l,:,:,1,2)))
			write(20,"(1x)")
			call mwrite(20,dimag(Veff(l,:,:,2,1)))
			write(20,"(1x)")
		enddo
		stop
		!do l=1,mq
			!write(20,"(e14.4$)")(dimag(rVeff(i,l,l,1,1)),i=1,mro)
			!write(20,"(1x)")
		!enddo
		!write(20,"(1x)")
		!do l=1,mq
			!write(20,"(e14.4$)")(dimag(rVeff(i,mq/2,l,1,1)),i=1,mro)
			!write(20,"(1x)")
		!enddo
		!write(20,"(1x)")
		!call mwrite(20,dimag(rVeff(1,:,:,1,1)))
		!write(20,"(1x)")
		!call mwrite(20,dimag(rVeff(2,:,:,1,1)))
		!write(20,"(1x)")
		!call mwrite(20,dimag(rVeff(3,:,:,1,1)))
		!write(20,"(1x)")
		!stop
		sfeg=0d0
		call selfefft(Tk,G0,Veff,sfeg)
		!ik=(/mq/2,int(0.15d0/2d0*mq)/)
		!call pade(iomgf(1:mo/2),sfeg(1:mo/2,mq/2,int(0.15d0/2d0*mq),1,1),romg,rsfeg(:,1,1,1,1))
		!write(10,"(3e14.4)")(real(romg(i)),rsfeg(i,1,1,1,1),i=1,mro)
		!write(10,"(1x/)")
		!$OMP PARALLEL DO SCHEDULE(GUIDED)
		do i=1,mq
			do j=1,mq
				Gf(:,i,j,:,:)=G0_inv(:,i,j,:,:)-sfeg(:,i,j,:,:)
				do l=1,mo
					call matrix_inv(Gf(l,i,j,:,:))
				enddo
				call pade(iomgf(1:mo/2),Gf(1:mo/2,i,j,1,1),romg,rGf(:,i,j,1,1))
			enddo
		enddo
		!$OMP END PARALLEL DO
		call mwrite(20,dimag(rGf(int(1.2d0/4d0*mro),:,:,1,1)))
		write(20,"(1x)")
		do l=1,mq/2
			write(20,"(e14.4$)")(dimag(rGf(i,l,l,1,1)),i=1,mro)
			write(20,"(1x)")
		enddo
		do l=mq/2,1,-1
			write(20,"(e14.4$)")(dimag(rGf(i,mq/2,l,1,1)),i=1,mro)
			write(20,"(1x)")
		enddo
		do l=mq/2,1,-1
			write(20,"(e14.4$)")(dimag(rGf(i,l,1,1,1)),i=1,mro)
			write(20,"(1x)")
		enddo
		do l=1,mq/2+1
			write(20,"(e14.4$)")(dimag(rGf(i,l,mq/2+2-l,1,1)),i=1,mro)
			write(20,"(1x)")
		enddo
		write(20,"(1x)")
	enddo
end program
