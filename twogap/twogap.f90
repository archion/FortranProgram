module global
	implicit none
	save
	real(8), parameter ::t(5)=(/1d0,-0.2d0,0.0d0,0d0,0d0/),escal=0.25d0,pi=3.1415926d0,cvg=1e-5,U=1.24d0,V=-2.4d0
	integer, parameter :: mk=512,mth=1000,ms=100,mr=128
	complex(8),parameter :: img=(0d0,1d0)
	character :: selt="s"
end module
program main
	use global, only : pi,cvg,t,U,V
	implicit none
	complex(8) :: Uk(4,4)
	real(8) :: kx,ky,sp,bt,nf=0.8d0,dt(2),pdt(2),dd,pdd,ek(4),gap,Tk
	character(100) :: pmt
	integer :: i,j
	write(pmt,"('_DDW=',f6.2,'_SC=',f6.2,'_t1=',f6.2,'_t2=',f6.2,'_t3=',f6.3)")U,V,t(1),t(2),t(3)
	open(unit=10,file="../data/energy.dat")
	open(unit=20,file="../data/fermi.dat")
	open(unit=40,file="../data/gaptemp"//trim(adjustl(pmt))//".dat")
	open(unit=50,file="../data/raman"//trim(adjustl(pmt))//".dat")
	open(unit=60,file="../data/fhase"//trim(adjustl(pmt))//".dat")
	call gnuplot
	!call raman(0.5d0,(/0.00d0,0d0/),0.05d0,170d0,(/0d0,3d0/),0.002d0)
	!call band(0.5d0,(/0.00d0,0d0/),0.05d0)
	!!!call mingap(0d0,(/0.0d0,0d0/),-0.85d0,0d0,gap,selt)
	!stop
	do i=60,90
		nf=1d0-0.0025d0*i
		pdt=0d0
		pdd=0d0
		do j=0,300,5
			Tk=0.01+j
			dt=dt+0.001
			dd=dd+0.001
			!call selfconsist(Tk,nf,dt,dd,sp)
			call selfconsist_tg(Tk,nf,dt,dd,sp)
			if(((dt(1)-cvg*100)*(pdt(1)-cvg*100)<0d0).or.((dd-cvg*100)*(pdd-cvg*100)<0d0).or.dt(1)<cvg*100d0.and.dd<cvg*100d0) then
				write(60,"(5(a6,e12.4))")"n=",nf,",Tk=",Tk,"sp=",sp,",DSC=",dt(1),",DDW=",dd
			endif
			if((nf-0.8075)<cvg) then
				call raman(dd,dt,sp,Tk,(/0d0,3d0/),0.002d0)
			endif
			if(dt(1)<cvg*100d0.and.dd<cvg*100d0) then
				exit
			endif
			pdt=dt
			pdd=dd
			!call band(dd,dt,sp)
		enddo
	enddo
	write(10,*)"e"
	write(40,*)"e"
	write(50,*)"e"
	write(60,*)"e"
	close(10)
	close(20)
	close(30)
	close(40)
	close(50)
	close(60)
end
subroutine raman(dd,dt,sp,Tk,omgr,domg)
	use global
	implicit none
	integer :: i,j,k,l,m1,m2
	complex(8) :: Uk(4,4)
	real(8) :: kx,ky,sp,bt,dt(2),dd,ek(4),gm(4,2),Tr(2),fk(4),R(2),omg,omgr(2),domg,Tk
	bt=escal/Tk*1.16e4
	omg=omgr(1)
	do while(omg<omgr(2))
		omg=omg+domg
		R=0d0
		!$OMP PARALLEL DO REDUCTION(+:R) PRIVATE(kx,ky,gm,ek,Uk,fk,Tr) SCHEDULE(GUIDED)
		do i=0,mr
			kx=pi/mr*i
			do j=0,mr-i
				ky=pi/mr*j
				gm(1:2,1)=(cos(kx)-cos(ky))*(/(t(1)+6d0*t(4)+8d0*t(3)*cos(ky)+cos(kx)*(8d0*t(3)+12d0*t(4)*cos(ky))),&
					-(t(1)+6d0*t(4)-8d0*t(3)*cos(ky)+cos(kx)*(-8d0*t(3)+12d0*t(4)*cos(ky)))/)
				gm(1:2,2)=sin(kx)*sin(ky)*(/(-4d0*t(2)-16d0*t(4)*cos(kx)-16d0*t(4)*cos(ky)),&
					(-4d0*t(2)+16d0*t(4)*cos(kx)+16d0*t(4)*cos(ky))/)
				gm(3:4,:)=-gm(1:2,:)
				call EU(kx,ky,dd,dt,sp,ek,Uk)
				fk=1d0/(1d0+exp(bt*ek))
				do m1=1,4
					do m2=1,4
						Tr=(/abs(sum(gm(:,1)*dconjg(Uk(:,m2))*Uk(:,m1)))**2,abs(sum(gm(:,2)*dconjg(Uk(:,m2))*Uk(:,m1)))**2/)
						R=R+Tr*(-fk(m1)+fk(m2))*DIMAG(1d0/(omg+ek(m1)-ek(m2)+img*domg*5d0))
					enddo
				enddo
			enddo
		enddo
		!$OMP END PARALLEL DO
		R=2d0*R/mr**2
		write(50,"(5e16.3)")omg,R
		write(50,"(1X/)")
	enddo
end
subroutine mingap(dd,dt,sp,th,gap)
	use global
	implicit none
	integer :: i,l
	complex(8) :: Uk(4,4)
	real(8) :: kx,ky,th,ek(4),sp,dd,dt(2),gap
	gap=1000d0
	! !$OMP PARALLEL DO PRIVATE(kx,ky,gk,dk,u2,v2,cos2,sin2,hn,tmp,tmpa,e1,e2,e) SCHEDULE(GUIDED)
	do i=0,mth
		ky=pi-i*pi/mth
		kx=pi-i*pi/mth*tan(th/180d0*pi)
		call EU(kx,ky,dd,dt,sp,ek,Uk)
		do l=1,4
			if(abs(ek(l))<gap.and.ek(l)<0.and.abs(Uk(1,l))>cvg) then
				gap=abs(ek(l))
			endif
		enddo
	enddo
	! !$OMP END PARALLEL DO
end
!subroutine selfconsist(Tk,nf,dt,dd,rsp)
	!use global
	!implicit none
	!real(8) :: bt,Tk,nf
	!complex(8) :: Uk(4,4),cth,sth,cfy(2),sfy(2)
	!real(8) :: kx,ky,n1,dd,ddp,ddk,ek(4),rsp,sp=0d0,sa,sb,sp0,dk,dt(2),dtp(2),dta,dtb,gk(2),wide,th,fy(2),fk(4)
	!integer :: i,j,c,info
	!logical :: flaga,flagb,dtflaga,dtflagb,ddflaga,ddflagb
	!! dtp(2)=-2d0
	!ddp=dd+100*cvg
	!wide=0.5d0
	!sp0=sp+wide
	!bt=escal/Tk*1.16e4
	!c=0
	!do while(abs(ddp-dd)>cvg)
		!dd=ddp
		!dtflaga=.true.
		!dtflagb=.true.
		!dta=dt(1)
		!dtb=dt(1)
		!dtp=0d0
		!do while(abs(dtp(1)-dt(1))>cvg)
			!write(*,*)dta,dt(1),dtb
			!write(*,*)dtflaga,dtflagb,dtp(1)
			!dt(1)=0.5d0*(dta+dtb)
			!flaga=.true.
			!flagb=.true.
			!sa=sp
			!sb=sp
			!n1=0d0
			!do while(abs(n1-nf)>cvg)
				!sp=0.5d0*(sa+sb)
				!ddp=0d0
				!dtp=0d0
				!c=c+1
				!!$OMP PARALLEL DO REDUCTION(+:n1,ddp,dtp) PRIVATE(kx,ky,gk,ek,Uk)&
				!!$OMP SCHEDULE(GUIDED)
				!do i=0,mk
					!do j=0,min(i,mk-i)
						!kx=pi/mk*i
						!ky=pi/mk*j
						!call EU(kx,ky,dd,dt,sp,ek,Uk)
						!fk=1d0/(1d0+dexp(bt*ek))
						!gk(1)=0.5d0*(cos(kx)-cos(ky))
						!n1=n1+2d0+dot_product(Uk(1,:)*dconjg(Uk(1,:))+Uk(2,:)*dconjg(Uk(2,:))-&
							!Uk(3,:)*dconjg(Uk(3,:))-Uk(4,:)*dconjg(Uk(4,:)),fk)
						!ddp=ddp+dimag(dot_product(Uk(1,:)*dconjg(Uk(2,:))-Uk(2,:)*dconjg(Uk(1,:))-&
							!Uk(4,:)*dconjg(Uk(3,:))+Uk(3,:)*dconjg(Uk(4,:)),fk)*gk(1))
						!dtp=dtp+dot_product(Uk(1,:)*dconjg(Uk(3,:))-Uk(2,:)*dconjg(Uk(4,:)),fk)*gk(1)
					!enddo
				!enddo
				!!$OMP END PARALLEL DO
				!n1=n1/(mk**2)*2d0
				!ddp=U*ddp/(mk**2)*2d0
				!dtp=V*dtp/(mk**2)*2d0
				!!write(*,"(4(a6,e12.4))")"sa=",sa,",sp=",sp,",sb=",sb,"n=",n1
				!!write(*,*)flaga,flagb,nf
				!if(abs(n1-nf)<=cvg) then
					!exit
				!endif
				!if(n1<nf) then
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
			!enddo
			!wide=max(abs(sp0-sp),100*cvg)
			!sp0=sp
			!if(abs(dtp(1)-dt(1))<=cvg) then
				!exit
			!endif
			!if(dtp(1)>dt(1)) then
				!dtflaga=.false.
				!dta=dt(1)
				!if(dtflaga.or.dtflagb) then
					!dtb=dt(1)+0.05
				!endif
			!else
				!dtflagb=.false.
				!dtb=dt(1)
				!if(dtflaga.or.dtflagb) then
					!dta=dt(1)-0.05
				!endif
			!endif
		!enddo
	!enddo
	!!write(*,*)"!!!!!!selfconsist return!!!!!!!!"
	!write(*,"(4(a6,e12.4)a6i4)")"n=",n1,",Tk=",Tk,",DSC=",dt(1),",DDW=",dd,"num:",c
	!write(40,"(4(a6,e12.4)a6i4)")"n=",n1,",Tk=",Tk,",DSC=",dt(1),",DDW=",dd,"num:",c
	!!write(*,*)"!!!!!!!!!!!!end!!!!!!!!!!!!!!!!!"
	!rsp=sp
!end
subroutine selfconsist_tg(Tk,nf,dt,dd,rsp)
	use global
	implicit none
	real(8) :: Tk,nf
	complex(8) :: Uk(4,4),cth,sth,cfy(2),sfy(2)
	real(8) :: kx,ky,n1,dd,ddp,ddk,ek(4),rsp,sp=0d0,sa,sb,sp0,dk,dt(2),dtp(2),gk(2),wide,th,fy(2),cvg1
	integer :: i,j,c,info
	logical :: flaga,flagb
	! dtp(2)=-2d0
	wide=0.5d0
	sp0=sp+wide
	c=0
	cvg1=0.001
	do 
		sa=sp
		sb=sp
		flaga=.true.
		flagb=.true.
		do 			
			sp=0.5d0*(sa+sb)
			ddp=0d0
			dtp=0d0
			n1=0d0
			c=c+1
			!$OMP PARALLEL DO REDUCTION(+:n1,ddp,dtp) PRIVATE(kx,ky,gk,ek,Uk)&
			!$OMP SCHEDULE(GUIDED)
			do i=0,mk
				do j=0,min(i,mk-i)
					kx=pi/mk*i
					ky=pi/mk*j
					call EU(kx,ky,dd,dt,sp,ek,Uk)
					call order(kx,ky,ek,Uk,Tk,n1,ddp,dtp)
				enddo
			enddo
			!$OMP END PARALLEL DO
			n1=n1/(mk**2)*2d0
			ddp=U*ddp/(mk**2)*2d0
			dtp=V*dtp/(mk**2)*2d0
			!write(*,"(4(a6,e12.4))")"sa=",sa,",sp=",sp,",sb=",sb,"n=",n1
			!write(*,*)flaga,flagb,nf
			if(abs(n1-nf)<=cvg1) then
				exit
			endif
			if(n1<nf) then
				flaga=.false.
				sa=sp
				if(flaga.or.flagb) then
					sb=sp+wide
					sp=sb
				endif
			else
				flagb=.false.
				sb=sp
				if(flaga.or.flagb) then
					sa=sp-wide
					sp=sa
				endif
			endif
		enddo
		wide=max(abs(sp0-sp),100*cvg)
		sp0=sp
		cvg1=max(cvg,min(cvg1,0.1d0*(abs(dtp(1)-dt(1))+abs(ddp-dd))))
		if(abs(dtp(1)-dt(1))<cvg.and.abs(ddp-dd)<cvg) then
			exit
		endif
		dt=dtp
		dd=ddp
		!write(*,"(5(a6,e12.4)a6i4)")"n=",n1,",Tk=",Tk,",sp=",sp,",DSC=",dt(1),",DDW=",dd,"num:",c
	enddo
	!write(*,*)"!!!!!!selfconsist return!!!!!!!!"
	write(*,"(5(a6,e12.4)a6i4)")"n=",n1,",Tk=",Tk,",sp=",sp,",DSC=",dt(1),",DDW=",dd,"num:",c
	write(40,"(5(a6,e12.4)a6i4)")"n=",n1,",Tk=",Tk,",sp=",sp,",DSC=",dt(1),",DDW=",dd,"num:",c
	!write(*,*)"!!!!!!!!!!!!end!!!!!!!!!!!!!!!!!"
	rsp=sp
end
subroutine band(dd,dt,sp)
	use global
	implicit none
	integer :: i,l
	real(8) :: dd,dt(2),sp,kx,ky,ek(4)
	complex(8) :: Uk(4,4)
	do i=1,3*ms-1
		kx=pi*(min((i/ms)*ms,ms)+(1-i/ms)*mod(i,ms))/ms
		ky=pi*min(max((i-ms),0),3*ms-i)/ms
		call EU(kx,ky,dd,dt,sp,ek,Uk)
		write(10,"(8e16.3)")(ek(l),abs(Uk(1,l)),l=1,4)
	enddo
end
subroutine EU(kx,ky,dd,dt,sp,ek,Uk)
	use global
	implicit none
	complex(8) :: Uk(4,4),cth,sth,cfy(2),sfy(2)
	real(8) :: eka,eks,e1(2),e2(2),ek(4),kx,ky,dd,ddk,sp,dk,dt(2),gk(2),th,fy(2)
	eks=-4d0*t(2)*cos(kx)*cos(ky)-2d0*t(3)*(cos(2d0*kx)+cos(2d0*ky))-4d0*t(5)*cos(2d0*kx)*cos(2d0*ky)-sp
	eka=-2d0*t(1)*(cos(kx)+cos(ky))-4d0*t(4)*(cos(2d0*kx)*cos(ky)+cos(kx)*cos(2d0*ky))
	gk(1)=0.5d0*(cos(kx)-cos(ky))
	! gk(2)=0.5d0*(cos(3d0*kx)-cos(3d0*ky))
	select case(selt)
	case("d")
		ddk=dd*gk(1)
	case("s")
		ddk=-dd
	case default
		write(*,*)"unknown order"
		stop
	end select
	dk=dt(1)*gk(1)
	e1=eks+(/1d0,-1d0/)*sqrt(ddk**2+eka**2)
	e2=sqrt(e1**2+dk**2)
	ek=(/e2,-e2/)
	if(ddk>0) then
		th=acos(eka/sqrt(ddk**2+eka**2))/2d0
	else
		th=(2d0*pi-acos(eka/sqrt(ddk**2+eka**2)))/2d0
	endif
	if(dk>0) then
		fy=acos(e1/e2)/2d0
	else
		fy=(2d0*pi-acos(e1/e2))/2d0
	endif
	cth=dcmplx(cos(th))
	sth=dcmplx(sin(th))
	cfy=dcmplx(cos(fy))
	sfy=dcmplx(sin(fy))
	select case(selt)
	case("d")
	Uk=reshape((/   cth*cfy(1) , -img*sth*cfy(1) ,      cth*sfy(1) ,    img*sth*sfy(1) ,   & 
	           -img*sth*cfy(2) ,      cth*cfy(2) , -img*sth*sfy(2) ,       -cth*sfy(2) ,   & 
	               -cth*sfy(1) ,  img*sth*sfy(1) ,      cth*cfy(1) ,    img*sth*cfy(1) ,   & 
	           -img*sth*sfy(2) ,      cth*sfy(2) ,  img*sth*cfy(2) , cth*cfy(2)     /) , (/4 , 4/))
	case("s")
	Uk=reshape((/   cth*cfy(1) ,  sth*cfy(1) ,  cth*sfy(1) ,       -sth*sfy(1) ,   & 
	               -sth*cfy(2) ,  cth*cfy(2) , -sth*sfy(2) ,       -cth*sfy(2) ,   & 
	               -cth*sfy(1) , -sth*sfy(1) ,  cth*cfy(1) ,       -sth*cfy(1) ,   & 
	               -sth*sfy(2) ,  cth*sfy(2) ,  sth*cfy(2) ,        cth*cfy(2)     /) , (/4 , 4/))
	case default
		write(*,*)"unknown order"
		stop
	end select
end
subroutine order(kx,ky,ek,Uk,Tk,n,dd,dt)
	use global
	implicit none
	complex(8) :: Uk(4,4)
	real(8) :: kx,ky,ek(4),dd,dt(2),gk(2),n,Tk,fk(4),bt
	bt=escal/Tk*1.16e4
	fk=1d0/(1d0+exp(bt*ek))
	gk(1)=0.5d0*(cos(kx)-cos(ky))
	select case(selt)
	case("d")
		dd=dd+dimag(dot_product(Uk(1,:)*dconjg(Uk(2,:))-Uk(2,:)*dconjg(Uk(1,:))-&
			Uk(4,:)*dconjg(Uk(3,:))+Uk(3,:)*dconjg(Uk(4,:)),fk)*gk(1))
	case("s")
		dd=dd+dot_product(Uk(1,:)*dconjg(Uk(2,:))+Uk(2,:)*dconjg(Uk(1,:))+&
			Uk(4,:)*dconjg(Uk(3,:))+Uk(3,:)*dconjg(Uk(4,:)),fk)
	case default
		write(*,*)"unknown order"
		stop
	end select
	n=n+2d0+dot_product(Uk(1,:)*dconjg(Uk(1,:))+Uk(2,:)*dconjg(Uk(2,:))-&
		Uk(3,:)*dconjg(Uk(3,:))-Uk(4,:)*dconjg(Uk(4,:)),fk)
	dt=dt+dot_product(Uk(1,:)*dconjg(Uk(3,:))-Uk(2,:)*dconjg(Uk(4,:)),fk)*gk(1)
end
subroutine gnuplot()
	use global
	implicit none
	character(100) :: tl
	write(tl,"('t_1=',f6.2,',t_2=',f6.2,',t_3=',f6.3,'\nV_{DDW}=',f6.2,'\nV_{SC}=',f6.2)")U,V,t(1),t(2),t(3)
	!plot energy band
	write(10,"(A)")"set term pngcairo"
	write(10,"(A)")"set output 'energy.png'"
	write(10,"(A)")"unset key"
	write(10,"(A)")"set size square"
	write(10,"(A)")"set xlabel 'k'"
	write(10,"(A)")"set ylabel '能量'"
	write(10,"(A)")'set label "n=0.8075\nT=0K" at 5,0.8'
	write(10,"(A)")"plot [:][-1:1] for[i=1:4] '-' using 0:(column(2*i-1)):(column(i*2)) with points lt 1 pt 7 ps variable"
	write(10,"(A)")"#data"
	!plot Raman
	write(50,"(A)")"set term pngcairo"
	write(50,"(A)")"set output 'raman.png'"
	write(50,"(A)")"unset key"
	write(50,"(A)")"set ytics"
	write(50,"(A)")"set y2tics"
	write(50,"(A)")"plot [:][:] '-' using 1:2 with line axis x1y1, '' using 1:3 with line axis x1y2"
	write(50,"(A)")"#data"
	!plot gap_tmp
	write(40,"(A)")"reset"
	write(40,"(A)")"set term pngcairo"
	write(40,"(A)")"set output 'gap.png'"
	write(40,"(A)")"set termoption enhanced"
	write(40,"(A)")"unset key"
	write(40,"(A)")"set xtic 50"
	write(40,"(A)")"set xlabel 'T(K)'"
	write(40,"(A)")"set ylabel 'gap/t'"
	write(40,"(A)")"set title '温度依赖'"
	write(40,"(A)")'set label "'//trim(adjustl(tl))//'" at 150,0.11'
	write(40,"(A)")"plot '-' using 4:8 with linespoints pt 7 lw 5 , '' using 4:10 with linespoints pt 7 lw 5"
	write(40,"(A)")"#data"
	!plot phase diagram
	write(60,"(A)")"reset"
	write(60,"(A)")"set term pngcairo"
	write(60,"(A)")"set output 'phase.png'"
	write(60,"(A)")"set termoption enhanced"
	write(60,"(A)")"unset key"
	write(60,"(A)")"set xlabel 'n'"
	write(60,"(A)")"set ylabel 'T(K)'"
	write(60,"(A)")"set title '相图'"
	write(60,"(A)")'set label "'//trim(adjustl(tl))//'" at 0.61,470'
	write(60,"(A)")"set label ""DDW"" at 0.85,300"
	write(60,"(A)")"set label ""SC"" at 0.75,50"
	write(60,"(A)")"plot '-' using 2:4 with points pt 4 ps 0.6"
	write(60,"(A)")"#data"
end
