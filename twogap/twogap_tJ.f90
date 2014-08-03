module global
	implicit none
	save
	real(8), parameter ::t(5)=(/1d0,-0.25d0,0.1d0,0d0,0d0/),escal=0.125d0,pi=3.1415926d0,cvg=1e-6,&
		!V=0.12d0,DJ=0.35d0,&
		V=0.d0,DJ=0.25d0,&
		Vs=DJ-V,Vd=0.5d0*DJ+V
	real(8) :: nf=0.8d0
	!character(3) :: pgflag="ddw"
	character(3) :: pgflag="sdw"
	integer, parameter :: mk=512,ms=100,mf=128,mr=512
	complex(8),parameter :: img=(0d0,1d0)
end module
program main
	use global
	implicit none
	complex(8) :: Uk(4,4)
	real(8) :: kf(2,2),sp=0d0,bt,sc,psc,pg,ppg,ap,th,gap,Tk,Tp(100),pk(3,100),peak(2),&
		nd(3)=(/0.1d0,0.135d0,0.2d0/),Td(3)=(/0d0,40d0,70d0/),dltn=0.0025d0,dltT=-1d0,den(2),pden(2),dTk,feg,ddf(2)
	integer :: i,j,l
	if(pgflag=="ddw".or.pgflag=="sdw") then
		write(*,*)"calculate for: ",pgflag
	else
		write(*,*)"the pg must be ddw or sdw"
		stop
	endif
	call fileopen
	call gnuplot
	!do i=24,96,1
	!do i=24,72,8
	nf=1d0
	sp=0d0
	!do i=1,3
	do
		!nf=1d0-nd(i)
		write(80,"(e12.4$)")nf
		psc=1d0
		ppg=-1d0
		pden=1d0
		pg=0d0
		sc=0d0
		ap=0d0
		!call selfconsist_tg(0.01d0,psc,ppg,sp)
		!call band(ppg,psc,sp)
		!call mingap(ppg,psc,sp,Tk,(/pi,0/),(/pi,pi/),gap,kf)
		!call fermisurface(ppg,psc,sp,0d0)
		!write(*,*)mx,my,gap
		!!write(70,"(F5.2)")1d0-nf
		!Tk=0.001
		!do
			!Tk=Tk+1d0
			!sc=sc+0.1
			!pg=pg+0.1
			!call selfconsist_tg(Tk,sc,pg,sp)
			!if(sc<cvg*100d0) then
				!Tc=Tk
				!exit
			!endif
		!enddo
		Tk=0.0001d0
		dTk=0.002d0
		!Tk=0.12d0
		!dTk=-0.002d0
		i=1
		do
		!do j=1,3
			!Tk=Td(j)+0.001
			sc=sc+0.1d0
			pg=pg+0.1d0
			ap=ap+0.1d0
			!call selfconsist_tg(Tk,sc,pg,ap,sp)
			call selforder(pg,sc,ap,sp,Tk,1)
			!call freeenergy(pg,sc,ap,sp,Tk,feg,ddf)
			write(*,"(e12.4$)")nf,Tk,sp,sc,pg,ap,feg,ddf
			write(*,"(1X)")
			write(40,"(e12.4$)")nf,Tk,sc*Vs,pg*DJ*0.5d0,ap*Vd
			!write(40,"(e12.4$)")Tk,sc*Vs,pg*Vd
			!write(90,"(e12.4$)")Tk,feg
			!sc=0.01d0
			!pg=0.d0
			!sc=sc+0.1d0
			!pg=0d0
			!call selforder(pg,sc,ap,sp,Tk,1)
			!call freeenergy(pg,sc,ap,sp,Tk,feg,ddf)
			!write(*,"(e12.4$)")nf,Tk,sp,sc,pg,ap,feg,ddf
			!write(*,"(1X)")
			!write(40,"(e12.4$)")sc,pg
			!write(90,"(e12.4$)")feg
			!sc=0d0
			!pg=pg+0.1d0
			!call selforder(pg,sc,ap,sp,Tk,1)
			!call freeenergy(pg,sc,ap,sp,Tk,feg,ddf)
			!write(*,"(e12.4$)")nf,Tk,sp,sc,pg,ap,feg,ddf
			!write(*,"(1X)")
			!write(40,"(e12.4$)")sc,pg
			!write(90,"(e12.4$)")feg
			!sc=0d0
			!pg=0d0
			!!ap=0d0
			!call selforder(pg,sc,ap,sp,Tk,1)
			!call freeenergy(pg,sc,ap,sp,Tk,feg,ddf)
			!write(*,"(e12.4$)")nf,Tk,sp,sc,pg,ap,feg,ddf
			!write(*,"(1X)")
			!write(90,"(e12.4$)")feg
			!write(90,"(1X)")
			!write(40,"(e12.4$)")sc,pg
			write(40,"(1X)")
			!!!!!!!!!!!!!!!!!!!!!!meanfield begin!!!!!!!!!!!!!!!!!!!!!!!!
			if((sc-cvg*100)*psc<0d0) then
				psc=-psc
				write(80,"(e12.4$)")Tk
			endif
			if((pg-cvg*100)*ppg<0d0) then
				ppg=-ppg
				write(80,"(e12.4$)")Tk
			endif
			write(50,"(F5.0)")Tk
			!!!!!!!!!!!!!!!!!!!!!!!meanfield end!!!!!!!!!!!!!!!!!!!!!!!!
			!!!!!!!!!!!!!!!!!!!!!!!!!!instable begin!!!!!!!!!!!!!!!!!!!!!!!!
			!call rpainstable(pg,sc,ap,sp,Tk,den)
			!write(*,"(e10.3$)")nf,Tk,sp,den,pg,sc,ap
			!write(*,"(1X)")
			!!write(40,"(e16.4$)")Tk,den
			!if(den(1)*pden(1)<0d0.and.pg<cvg) then
				!pden(1)=-pden(1)
				!write(80,"(e12.4$)")Tk
				!pg=pg+0.1d0
				!if(dTk>0d0) then
					!exit
				!endif
			!endif
			!if(den(2)*pden(2)<0d0.and.sc<cvg) then
				!pden(2)=-pden(2)
				!write(80,"(e12.4$)")Tk
				!sc=sc+0.1d0
			!endif
			!if(sc>cvg.and.dTk<0) then
				!dTk=-dTk
				!Tk=dTk
				!sc=sc+0.1d0
				!pg=0d0
				!pden(1)=1d0
			!endif
			!if(dTk>0d0.and.sc<cvg*100) then
				!exit
			!endif
			!!!!!!!!!!!!!!!!!!!!!!!!!!instable end!!!!!!!!!!!!!!!!!!!!!!!!!
			!Tp(i)=Tk
			!call raman(pg,sc,ap,sp,Tp(i),(/0d0,0.4d0/),0.0004d0,pk(:,i))
			!if(Tk<0.01d0) then
				!pk0=pk
			!endif
			!call mingap(pg,sc,sp,Tk,(/pi,0d0/),(/pi,pi/2d0/),gap,kf)
			!do l=0,10
				!call EDC((/pi,pi/30d0*l/),pg,sc,sp,Tk,(/-0.1d0,0d0/),0.0001d0,peak,.true.)
			!enddo
			!write(*,"(3e12.4)")gap,mx,my
			!write(30,"(F5.0)")Tk
			!call EDC(pi,0d0,pg,sc,sp,Tk,(/-0.1d0,0d0/),0.0001d0,peak)
			!call EDC(mx,my,pg,sc,sp,Tk,(/-0.1d0,0d0/),0.0001d0,peak)
			!write(40,"(6e12.4)")nf,Tk,sp,sc,pg,peak(1)
			!write(40,"(7e12.4)")nf,Tk,sp,sc,pg,gap,kf(2)
			!write(40,"(5e12.4)")nf,Tk,sp,sc,pg
			!call fermisurface(pg,sc,sp,0d0)
			!if(sc<cvg*100d0.and.pg<cvg*100d0) then
			if(sc<cvg*100d0) then
			!if(Tk<abs(dTk)*1.5d0) then
			!if(Tk>0.1) then
				!write(70,"(4e16.3)")(Tp(j)/Tp(i),pk(:,j)/pk(:,1),j=1,i)
				!call DOS(0d0,sc,ap,sp,Tk,(/-0.5d0,0.7d0/),0.0001d0)
				!write(30,"(1X)")
				!call DOS(pg,sc,ap,sp,Tk,(/-0.5d0,0.7d0/),0.0001d0)
				!call fermisurface(pg,sc,ap,sp,0d0)
				!call fermisurface(0.1d0,0d0,ap,sp,0d0)
				!stop
				exit
			endif
			exit
			i=i+1
			Tk=Tk+dTk
			write(30,"(1X)")
		enddo
		!stop
		nf=nf-0.05d0
		!if(pg<cvg*100d0) then
		if(nf<0.2d0) then
			exit
		endif
		write(10,"(1X)")
		write(30,"(1X/)")
		write(50,"(1X)")
		!write(40,"(1X/)")
		write(70,"(1X/)")
		write(80,"(1X)")
	enddo
	call fileclose
end
subroutine rpainstable(pg,sc,ap,sp,Tk,den)
	use global
	implicit none
	integer :: i,j,l,m,n
	complex(8) :: Uk(4,4),Ukq(4,4)
	real(8) :: k(2),q(2),sp,bt,pg,sc,ek(4),ekq(4),fk(4),fkq(4),Tk,DJq,Xq(2),den(2),Vrpa(2),tran(2),ap
	bt=1d0/Tk
	Xq=0d0
	select case(pgflag)
	case("ddw")
		Vrpa=(/Vd,Vs/)
	case("sdw")
		Vrpa=(/DJ,Vs/)
	end select
	!$OMP PARALLEL DO REDUCTION(+:Xq) PRIVATE(k,ek,Uk,fk,ekq,Ukq,fkq,tran) SCHEDULE(GUIDED)
	do i=0,mr
		k(1)=pi/mr*i
		do j=0,mr
			k(2)=pi/mr*j
			call EU1(k,pg,sc,ap,sp,ek,Uk)
			call EU1(k+(/pi,pi/),pg,sc,ap,sp,ekq,Ukq)
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
end
subroutine raman(pg,sc,ap,sp,Tk,omgr,domg,pk)
	use global
	implicit none
	integer :: i,j,l,m1,m2
	complex(8) :: Uk(4,4),R(2)
	real(8) :: k(2),sp,bt,sc,pg,ek(4),gm(4,2),Tr(2),fk(4),omg,omgr(2),domg,Tk,peak(3,2),R_rpa,pk(3),DJp,ap,dp
	dp=abs(1d0-nf)
	!bt=escal/Tk*1.16e4
	bt=1d0/Tk
	omg=omgr(1)
	peak=0d0
	DJp=DJ/(dp*t(1)+DJ*ap)**2
	do while(omg<omgr(2))
		omg=omg+domg
		R=0d0
		!$OMP PARALLEL DO REDUCTION(+:R) PRIVATE(k,gm,ek,Uk,fk,Tr) SCHEDULE(GUIDED)
		do i=0,mr
			k(1)=pi/mr*i
			do j=0,mr-i
				k(2)=pi/mr*j
				gm(1:2,1)=(cos(k(1))-cos(k(2)))*(dp*t(1)+ap*Vd)*(/1d0,-1d0/)+&
					4d0*t(3)*dp*(cos(2d0*k(1))-cos(2d0*k(2)))
				gm(1:2,2)=sin(k(1))*sin(k(2))*dp*t(2)*(/-4d0,-4d0/)
				gm(3:4,:)=-gm(1:2,:)
				call EU1(k,pg,sc,ap,sp,ek,Uk)
				fk=1d0/(1d0+exp(bt*ek))
				do m1=1,4
					do m2=1,4
						Tr=(/abs(sum(gm(:,1)*dconjg(Uk(:,m2))*Uk(:,m1)))**2,abs(sum(gm(:,2)*dconjg(Uk(:,m2))*Uk(:,m1)))**2/)
						R=R+Tr*(fk(m1)-fk(m2))*1d0/(omg+ek(m1)-ek(m2)+img*0.04d0)
					enddo
				enddo
			enddo
		enddo
		!$OMP END PARALLEL DO
		R=1d0*R/mr**2
		!R_rpa=-dimag(R(1)/(1d0+DJp*R(1)))
		!R_rpa=-dimag(R(1))/((1d0+DJp*dreal(R(1)))**2+(DJp*dimag(R(1)))**2)
		!write(50,"(e16.3$)")omg,-dimag(R),R_rpa,Djp,dreal(R(1)),(1d0+DJp*dreal(R(1)))**2+(DJp*dimag(R(1)))**2
		write(50,"(e16.3$)")omg,-dimag(R)
		write(50,"(1X)")
		if(peak(1,2)<-dimag(R(1))) then
			peak(1,:)=(/omg,-dimag(R(1))/)
		endif
		if(peak(2,2)<-dimag(R(2))) then
			peak(2,:)=(/omg,-dimag(R(2))/)
		endif
		!if(peak(3,2)<R_rpa) then
			!peak(3,:)=(/omg,R_rpa/)
		!endif
	enddo
	pk=peak(:,1)
	write(50,"(1X)")
end
!subroutine susp(qx,qy,omg,X)
	!implicit none
	!call sfcsap(Tk,sp,r)
	!call bubb(qx,qy,omg,B)
	!X=(1d0-nf)**2/4d0*(8d0/DJ*r**2-B)
!end
!subroutine bubb(qx,qy,omg,B)
	!implicit none
	!call sfcsap(Tk,sp,r)
	!do i
		!k(1)=
		!kqx=
		!do j
			!k(2)=
			!kqy=
			!ek=
			!ekq=
			!fk=
			!fkq=
			!g=2d0*r*(sin(k(1)-qx/2d0)-sin(k(2)-qy/2d0))
			!B=B+g*(fk-fkq)/(omg+ek-ekq+img*0.01d0)
		!enddo
	!enddo
!end
subroutine EDC(k,pg,sc,ap,sp,Tk,omgr,domg,peak,fw)
	use global
	implicit none
	integer :: i
	logical :: fw
	complex(8) :: Uk(4,4)
	real(8) :: k(2),sp,bt,sc,pg,ek(4),fk(4),A,omg,omgr(2),domg,Tk,peak(2),ap
	!bt=escal/Tk*1.16e4
	bt=1d0/Tk
	omg=omgr(1)
	call EU1(k,pg,sc,ap,sp,ek,Uk)
	fk=1d0/(1d0+exp(bt*ek))
	peak=0d0
	!fk=1d0
	do while(omg<omgr(2))
		omg=omg+domg
		A=-sum(DIMAG(fk*Uk(1,:)*dconjg(Uk(1,:))/(omg-ek+img*domg*200d0)))
		if(peak(2)<A) then
			peak=(/-omg,A/)
		endif
		if(fw) then
			write(30,"(5e16.3)")-omg,A,k(2),nf,Tk
		endif
	enddo
	if(fw) then
		write(30,"(1X)")
	endif
end
subroutine DOS(pg,sc,ap,sp,Tk,omgr,domg)
	use global
	integer :: i,j
	complex(8) :: Uk(4,4)
	real(8) :: k(2),sp,bt,sc,pg,ek(4),fk(4),A,omg,omgr(2),domg,Tk,ap
	bt=1d0/Tk
	omg=omgr(1)
	do while(omg<omgr(2))
		A=0d0
		omg=omg+domg
		!$OMP PARALLEL DO REDUCTION(+:A) PRIVATE(k,ek,Uk,fk) SCHEDULE(STATIC)
		do i=1,mk
			do j=1,mk
				k=(/pi/mk*i,pi/mk*j/)
				call EU1(k,pg,sc,ap,sp,ek,Uk)
				!fk=1d0/(1d0+exp(bt*ek))
				A=A-sum(DIMAG(Uk(1,:)*dconjg(Uk(1,:))/(omg-ek+img*domg*20d0)))
			enddo
		enddo
		!$OMP END PARALLEL DO
		A=A/(mk**2)
		write(30,"(e16.3$)")omg,A
		write(30,"(1X)")
		!write(*,"(e16.3$)")omg,A
		!write(*,"(1X)")
	enddo
end
!subroutine band(pg,sc,sp,ki,kf)
	!use global
	!implicit none
	!integer :: i,l
	!real(8) :: pg,sc,sp,k(2),ek(4)
	!complex(8) :: Uk(4,4)
	!kn=ki
	!!do i=1,3*ms-1
	!do
		!kn=kn+(kf-ki)/sum((kf-ki)**2)*pi/mk
		!!k(1)=pi*(min((i/ms)*ms,ms)+(1-i/ms)*mod(i,ms))/ms
		!!k(2)=pi*min(max((i-ms),0),3*ms-i)/ms
		!call EU1(kn,pg,sc,ap,sp,ek,Uk)
		!eks=-2d0*(1d0-nf)*t(2)*cos(kn(1))*cos(kn(2))-sp
		!eka=-((1d0-nf)*t(1)+ap)*(cos(kn(1))+cos(kn(2)))
		!write(10,"(10e16.3)")(ek(l),abs(Uk(1,l)),l=1,4),&
			!eks+(/1d0,-1d0/)*sqrt((0.5d0*(cos(kn(1))-cos(kn(2)))*pg)**2+eka**2)*sign(1d0,eka)
		!!write(10,"(2e16.3)")(ek(l),abs(Uk(1,l)),l=4,4)
		!if(all((kn-kf)*(kf-ki)>=0d0)) then
			!exit
		!endif
	!enddo
	!write(10,"(1X)")
!end
subroutine mingap(pg,sc,ap,sp,Tk,ki,kf,gap,km)
	use global
	implicit none
	integer :: i,l
	complex(8) :: Uk(4,4)
	real(8) :: ki(2),kf(2),km(2),kn(2),ek(4),sp,pg,sc,gap,Tk,peak(2),eks,eka,ap
	gap=1000d0
	kn=ki
	do
		kn=kn+(kf-ki)/sum((kf-ki)**2)*pi/mk
		!!do l=1,4
		call EU1(kn,pg,sc,ap,sp,ek,Uk)
		!if(abs(Uk(1,2))>5e-2) then
			!if(abs(ek(2))<gap.and.ek(2)<0d0) then
				!gap=abs(ek(2))
				!km=kn
			!endif
		!elseif(abs(ek(1))<gap.and.ek(1)<0d0) then
			!gap=abs(ek(1))
			!km=kn
		!endif
		call EDC(kn,pg,sc,sp,Tk,(/-0.1d0,0d0/),0.0001d0,peak,.false.)
		if(peak(1)<gap) then
			gap=peak(1)
			km=kn
		endif
		!call EDC(kn,pg,sc,sp,Tk,(/-0.2d0,0d0/),0.0001d0,peak,.false.)
		!if(peak(1)<gap) then
			!gap=peak(1)
			!km=kn
		!endif
		!export band
		!eks=-2d0*(1d0-nf)*t(2)*cos(kn(1))*cos(kn(2))-sp
		!eka=-((1d0-nf)*t(1)+ap)*(cos(kn(1))+cos(kn(2)))
		!write(10,"(10e16.3)")(ek(l),abs(Uk(1,l)),l=1,4),&
			!eks+(/1d0,-1d0/)*sqrt((0.5d0*(cos(kn(1))-cos(kn(2)))*pg)**2+eka**2)*sign(1d0,eka)
		write(10,"(11e16.3)")kn(2),(ek(l),abs(Uk(1,l)),l=1,4)
		if(all((kn-kf)*(kf-ki)>=0d0)) then
			exit
		endif
	enddo
	call EDC(km,pg,sc,sp,Tk,(/-0.1d0,0d0/),0.0001d0,peak,.true.)
	write(10,"(1X)")
	write(40,"(9e12.4)")nf,Tk,sp,sc,pg,gap,-0.5d0*(cos(km(1))-cos(km(2)))*sc,-0.5d0*(cos(km(1))-cos(km(2)))*pg,km(2)
end
subroutine fermisurface(pg,sc,ap,sp,omg)
	use global
	implicit none
	integer :: i,j,l
	complex(8) :: Uk(4,4)
	real(8) :: k(2),ek(4),sp,pg,sc,ap,omg,A
	do i=0,mf
		do j=0,mf
			k=(/pi/mf*i,pi/mf*j/)
			call EU1(k,pg,sc,ap,sp,ek,Uk)
			A=0d0
			do l=1,4
				A=A-dimag(abs(Uk(1,l))**2/(omg+img*0.01d0-ek(l)))
			enddo
			write(20,"(e16.3$)")A
		enddo
		write(20,"(1X)")
	enddo
	write(20,"(1X/)")
end
subroutine selfconsist_tg(Tk,sc,pg,ap,sp)
	use global
	implicit none
	real(8) :: Tk
	complex(8) :: Uk(4,4)
	real(8) :: k(2),np,pnp,pg,pgp,ek(4),sp,sa,sb,sp0,sc,scp,wide,cvg1,al,ap,app
	integer :: i,j,c,info
	logical :: flaga,flagb
	!scp(2)=-2d0
	wide=0.5d0
	al=1d0
	sp0=sp+wide
	c=0
	cvg1=0.0001
	do 
		!sa=sp
		!sb=sp
		!flaga=.true.
		!flagb=.true.
		pnp=nf+0.1d0
		do 			
			!sp=0.5d0*(sa+sb)
			pgp=0d0
			scp=0d0
			app=0d0
			np=0d0
			c=c+1
			!call sfcsap(Tk,sp,ap)
			!$OMP PARALLEL DO REDUCTION(+:np,pgp,scp,app) PRIVATE(k,ek,Uk) SCHEDULE(GUIDED)
			do i=0,mk
				do j=0,min(i,mk-i)
					k=(/pi/mk*i,pi/mk*j/)
					call EU1(k,pg,sc,ap,sp,ek,Uk)
					call order(k,ek,Uk,Tk,np,pgp,scp,app)
				enddo
			enddo
			!$OMP END PARALLEL DO
			np=np/(mk**2)*2d0
			pgp=pgp/(mk**2)*2d0
			scp=scp/(mk**2)*2d0
			app=app/(mk**2)*2d0
			!write(*,"(2(a6,e12.4))")"sp=",sp,"n=",np
			!write(*,*)abs(np-nf),c
			if((abs(nf-np)/abs(nf-pnp))>0.8d0) then
				al=max(al-0.05d0,0.05)
			endif
			if(abs(nf-np)<=cvg) then
				exit
			else
				sp=sp+al*(nf-np)
			endif
			pnp=np
			!if(np<nf) then
				!flaga=.false.
				!sa=sp
				!if(flagb) then
					!sb=sp+wide
				!endif
			!else
				!flagb=.false.
				!sb=sp
				!if(flaga) then
					!sa=sp-wide
				!endif
			!endif
		enddo
		!wide=max(abs(sp0-sp),100*cvg)
		!sp0=sp
		!cvg1=max(cvg,min(cvg1,0.1d0*(abs(scp-sc)+abs(pgp-pg))))
		!write(*,"(5(a6,e12.4)a6i4)")"n=",np,",Tk=",Tk,",sp=",sp,",DSC=",sc,",DDW=",pg,"num:",c
		!app=ap
		if((abs(scp-sc)+abs(pgp-pg)+abs(app-ap))<cvg) then
			exit
		endif
		sc=scp
		pg=pgp
		ap=app
	enddo
	!write(*,"(5e12.4,i4,e10.3)")np,Tk,sp,sc,pg,c,al
end
subroutine freeenergy(pg,sc,ap,sp,Tk,feg,ddf)
	use global
	implicit none
	real(8) :: k(2),eka,eks,e1,e2(2),E(2),Tk,bt,dp,feg,pg,sc,ap,sp,gk,pgk,sck,ddf(2),tmp(2),s
	integer :: i,j,l,Ns
	dp=abs(1d0-nf)
	bt=1d0/Tk
	feg=0d0
	ddf=0d0
	Ns=0
	!$OMP PARALLEL DO REDUCTION(+:feg,ddf,Ns) PRIVATE(k,eks,eka,gk,sck,pgk,e1,e2,E,tmp,s) SCHEDULE(GUIDED)
	do i=1,mk
		do j=0,min(i,mk-i)
			k=(/pi/mk*i,pi/mk*j/)
			eks=-4d0*dp*t(2)*cos(k(1))*cos(k(2))-sp-2d0*dp*t(3)*(cos(2d0*k(1))+cos(2d0*k(2)))
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
		feg=feg/Ns+(4d0*(Vd*(pg**2+ap**2)+Vs*sc**2)-sp)+sp*nf
	case("sdw")
		feg=feg/Ns+(4d0*(DJ*pg**2+Vd*ap**2+Vs*sc**2)-sp)+sp*nf
	end select
	ddf=sign(1d0,ddf)
end
subroutine freeorder(Tk,pg,sc,ap,sp,np,pgp,scp,app)
	use global
	implicit none
	real(8) :: k(2),eka,eks,e1,e2(2),E(2),Tk,bt,dp,tmp,pg,sc,ap,sp,np,pgp,scp,app,pgk,sck,gk,s
	integer :: i,j,l,Ns
	dp=abs(1d0-nf)
	bt=1d0/Tk
	Ns=0
	!$OMP PARALLEL DO REDUCTION(+:np,pgp,scp,app,Ns) PRIVATE(k,eks,eka,gk,sck,pgk,e1,e2,E,tmp,s) SCHEDULE(GUIDED)
	do i=1,mk
		do j=0,min(i,mk-i)
			k=(/pi/mk*i,pi/mk*j/)
			eks=-4d0*dp*t(2)*cos(k(1))*cos(k(2))-sp-2d0*dp*t(3)*(cos(2d0*k(1))+cos(2d0*k(2)))
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
end
subroutine selforder(pg,sc,ap,sp,Tk,sig)
	use global
	implicit none
	real(8) :: k(2),pg,sc,ap,sp,pgp,scp,app,np,pnp,al,Tk
	integer :: c,sig
	al=1d0
	do 
		pnp=nf+0.1d0
		do 			
			pgp=0d0
			scp=0d0
			app=0d0
			np=0d0
			c=c+1
			call freeorder(Tk,pg,sc,ap,sp,np,pgp,scp,app)
			!write(*,"(2(a6,e12.4))")"sp=",sp,"n=",np
			if((abs(nf-np)/abs(nf-pnp))>0.8d0) then
				al=max(al-0.05d0,0.05)
			endif
			if(abs(nf-np)<=cvg) then
				exit
			else
				sp=sp+al*(nf-np)
			endif
			pnp=np
		enddo
		if(sig==0) then
			app=ap
			scp=sc
			pgp=pg
		endif
		if((abs(scp-sc)+abs(pgp-pg)+abs(app-ap))<cvg) then
			exit
		endif
		sc=scp
		pg=pgp
		ap=app
	enddo
	!write(*,"(e12.4$)")np,Tk,sp,sc,pg,ap
	!write(*,"(1X)")
end
subroutine EU1(k,pg,sc,ap,sp,ek,Uk)
	use global
	implicit none
	complex(8) :: Uk(4,4),cth,sth,cfy(2),sfy(2)
	real(8) :: eka,eks,e1(2),e2(2),ek(4),k(2),pg,pgk,sp,sck,sc,gk,th,fy(2),cos2th,sin2th,cos2fy(2),sin2fy(2),ap,dp
	dp=abs(1d0-nf)
	!dp=abs(1d0)
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
end
subroutine EU2(k,pg,sc,ap,sp,ek,Uk)
	use global
	implicit none
	complex(8) :: Uk(4,4),cth,sth,cfy(2),sfy(2)
	real(8) :: eka,eks,e1(2),e2(2),ek(4),k(2),pg,pgk,sp,sck,sc,gk,th,fy(2),cos2th,sin2th,cos2fy(2),sin2fy(2),ap,dp
	dp=abs(1d0-nf)
	!dp=abs(1d0)
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
	ek=(/e2,-e2/)
	cos2th=eka/sqrt(pgk**2+eka**2)
	sin2th=pgk/sqrt(pgk**2+eka**2)
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
end
subroutine sfcsap(Tk,sp,ap)
	use global, only : mk,cvg,escal,pi,t,DJ,nf
	implicit none
	real(8) :: k(2),Tk,bt,sp,ap,app,eks,eka
	integer :: i,j
	!bt=escal/Tk*1.16e4
	bt=1d0/Tk
	do
		app=0d0
		!$OMP PARALLEL DO REDUCTION(+:app) PRIVATE(k,eks,eka) SCHEDULE(GUIDED)
		do i=0,mk
			do j=0,mk
				k=(/pi/mk*i,pi/mk*j/)
				eks=-4d0*(1d0-nf)*t(2)*cos(k(1))*cos(k(2))-sp-2d0*(1d0-nf)*t(3)*(cos(2d0*k(1))+cos(2d0*k(2)))
				eka=-2d0*((1d0-nf)*t(1)+ap)*(cos(k(1))+cos(k(2)))
				app=app+cos(k(1))*1d0/(1d0+exp(bt*(eks+eka)))
			enddo
		enddo
		!$OMP END PARALLEL DO
		app=DJ*app/mk**2
		if(abs(app-ap)<cvg) then
			exit
		endif
		ap=app
	enddo
end
subroutine order(k,ek,Uk,Tk,n,pg,sc,ap)
	use global
	implicit none
	complex(8) :: Uk(4,4)
	real(8) :: k(2),ek(4),pg,sc,gk,n,Tk,fk(4),bt,ap
	!bt=escal/Tk*1.16e4
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
end
subroutine fileclose()
	implicit none
	close(10)
	close(20)
	close(30)
	close(40)
	close(50)
	close(60)
	close(70)
	close(80)
end
subroutine fileopen()
	implicit none
	open(unit=10,file="../data/energy.dat")
	open(unit=20,file="../data/fermi.dat")
	open(unit=30,file="../data/EDC.dat")
	open(unit=40,file="../data/gaptemp.dat")
	open(unit=50,file="../data/raman.dat")
	open(unit=60,file="../data/gaptheta.dat")
	open(unit=70,file="../data/ramantemp.dat")
	open(unit=80,file="../data/phase_tJ.dat")
	open(unit=90,file="../data/freeenergy.dat")
end
subroutine gnuplot()
	implicit none
	integer :: i
	do i=10,90,10
		write(i,"(A)")'reset'
		write(i,"(A)")'#custom'
		write(i,"(A)")'#multiplot'
		write(i,"(A)")'ix=15'
		write(i,"(A)")'iy=1'
		write(i,"(A)")'sx=2*ix'
		write(i,"(A)")'sy=4*iy'
		write(i,"(A)")'xlb="T"'
		write(i,"(A)")'ylb="gap size"'
		write(i,"(A)")'ml=0.8/sx'
		write(i,"(A)")'mb=0.6/sy'
		write(i,"(A)")'mr=0.1/sx'
		write(i,"(A)")'mt=0.1/sy'
		write(i,"(A)")'gx=(1.-ml-mr)/ix'
		write(i,"(A)")'gy=(1.-mb-mt)/iy'
		write(i,"(A)")'#term'
		write(i,"(A)")'set term eps font ",18" size sx,sy'
		write(i,"(A)")'set output "-."."eps"'
		write(i,"(A)")'set label ylb font ",20" rotate by 90 center at character 0.8,screen 0.5*(1+mb-mt)'
		write(i,"(A)")'set label xlb font ",20" center at screen 0.5*(1+ml-mr),character 0.5'
		write(i,"(A)")'set xtics 30 out nomirror'
		write(i,"(A)")'set ytics 0.02 out nomirror'
		write(i,"(A)")'#'
		write(i,"(A)")'set size gx,gy'
		write(i,"(A)")'set tmargin 0'
		write(i,"(A)")'set bmargin 0'
		write(i,"(A)")'set lmargin 0'
		write(i,"(A)")'set rmargin 0'
		write(i,"(A)")'set multiplot'
		write(i,"(A)")'set mxtics 5'
		write(i,"(A)")'set mytics 5'
		write(i,"(A)")'unset xlabel'
		write(i,"(A)")'unset ylabel'
		write(i,"(A)")'unset key'
		write(i,"(A)")'do for[i=0:(ix*iy-1)]{'
		write(i,"(A)")'	oxi=i%ix'
		write(i,"(A)")'	oyi=i/ix+1'
		write(i,"(A)")'	set origin ml+oxi*gx,mb+(iy-oyi)*gy'
		write(i,"(A)")'	if(oxi!=0){'
		write(i,"(A)")'		set format y ""'
		write(i,"(A)")'	}'
		write(i,"(A)")'	if(oyi!=iy){'
		write(i,"(A)")'		set format x ""'
		write(i,"(A)")'	}'
		write(i,"(A)")'	if(i==(ix*iy-1)){'
		write(i,"(A)")'		set key font ",16" at screen 1-mr-0.01/sx,1-mt-0.01/sy horizontal maxcols 1 opaque autotitle columnhead'
		write(i,"(A)")'	}'
		write(i,"(A)")'	set label sprintf("(%1.3f)",0.11+0.005*i) at graph 0.1/(gx*sx),1.-0.15/(gy*sy)'
	enddo
	!plot energy
	write(10,"(A)")'	set xzeroaxis'
	write(10,"(A)")'	plot [0:20][-0.1:0.1] for[j=0:20] for[k=0:3] "-" index i every 5::0:j:250:j using ($1+j*1.2):2*k+2:2*k+3 with points pt 7 ps variable lc k+1,for [j=0:20] "gaptemp.dat" index i every ::j:0:j:0 using ($9+j*1.2):(-1*$6) with points pt 28 ps 1 lc 0 lw 2'
	!plot 
	write(20,"(A)")'	plot [:55][0:0.1] for[j=0:2] "-" index i using 2:j+4 with line title word("SC DDW Total",j+1)'
	!plot phase diagram
	write(30,"(A)")'	plot [:55][0:0.1] for[j=0:2] "-" index i using 2:j+4 with line title word("SC DDW Total",j+1)'
	!plot phase diagram
	write(40,"(A)")'	plot [:55][0:0.1] for[j=0:2] "-" index i using 2:j+4 with line title word("SC DDW Total",j+1)'
	!plot phase diagram
	write(50,"(A)")'	plot [:][0:] for[j=0:10] "-" index i/2 every :::j::j using 1:i%2+2:(1e10) smooth acsplines with lp pi 50 title "".(j*5)'
	!plot phase diagram
	write(60,"(A)")'	plot [:55][0:0.1] for[j=0:2] "-" index i using 2:j+4 with line title word("SC DDW Total",j+1)'
	!plot phase diagram
	write(70,"(A)")'	plot [:1.1][0.2:1.3] for[j=0:12:2] "-" using 1:i:(100.) lc (i/2+1) smooth acsplines notitle'
	!plot phase diagram
	write(80,"(A)")'set style fill transparent solid 0.3 border'
	write(80,"(A)")'set style rect fc lt -1 fs solid 0.15 noborder'
	write(80,"(A)")'set obj rect from 0.14, graph 0 to 0.175, graph 1 behind'
	write(80,"(A)")'plot [0.06:0.24][0:100] "-" using (1-$1):3:(1e7) smooth acsplines with filledcu y1=0 lc 3 lw 2, "-" u (1-$1):4:(1e9) smooth acsplines with filledcu x1=0 lc 1 lw 2, "-" u (1-$1):5:(1e9) smooth acsplines with filledcu x1=0 lc 1 lw 2 '
	write(90,"(A)")'	plot [:1.1][0.2:1.3] for[j=0:12:2] "-" using 1:i:(100.) lc (i/2+1) smooth acsplines notitle'
	do i=10,90,10
		write(i,"(A)")'	unset label'
		write(i,"(A)")'	unset format'
		write(i,"(A)")'}'
		write(i,"(A)")'#data'
	enddo
end
