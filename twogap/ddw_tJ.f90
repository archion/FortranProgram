module global
	implicit none
	save
	real(8), parameter ::t(5)=(/1d0,-0.35d0,0.000d0,0.d0,0d0/),escal=0.25d0,pi=3.1415926d0,cvg=1e-5,DJ=0.3d0,V=-0.48d0
	real(8) :: nf=0.8d0,ap=0.1d0
	integer, parameter :: mk=512,mth=1000,ms=100,mf=128,mr=512
	complex(8),parameter :: img=(0d0,1d0)
end module
program main
	use global, only : nf,pi,cvg,t,DJ,V
	implicit none
	complex(8) :: Uk(4,4)
	real(8) :: mx,my,sp=0d0,bt,dt(2),pdt(2),dd,pdd,th,gap,Tk,Tc,pk(3),pk0(3),peak(2),&
		nd(3)=(/0.11d0,0.14d0,0.17d0/),Td(3)=(/0d0,40d0,70d0/)
	integer :: i,j,k
	open(unit=10,file="../data/energy.dat")
	open(unit=20,file="../data/fermi.dat")
	open(unit=30,file="../data/EDC.dat")
	open(unit=40,file="../data/gaptemp.dat")
	open(unit=50,file="../data/raman.dat")
	open(unit=60,file="../data/phase_tJ.dat")
	open(unit=70,file="../data/gaptheta.dat")
	open(unit=80,file="../data/ramantemp.dat")
	call gnuplot
	!write(40,"('n T sp DSC DDW num')")
	!write(50,"('omg B1g B2g')")
	!write(60,"('n T sp DSC DDW')")
	!call raman(0.5d0,(/0.00d0,0d0/),0.05d0,170d0,(/0d0,3d0/),0.002d0)
	!call band(0.5d0,(/0.00d0,0d0/),0.05d0)
	!call mingap(0.2d0,(/1d0,0d0/),-0.847d0,0d0,gap,mx,my)
	!call EDC(mx,my,0.2d0,(/1d0,0d0/),-0.847d0,0.0001d0,(/-5d0,5d0/),0.002d0)
	!call fermisurface(0.5d0,(/0.d0,0d0/),-0.0d0,0d0)
	!write(*,*)gap,mx,my
	!stop
	do i=1,3
	!do i=54,58,1
		nf=1d0-0.0025d0*i
		nf=1d0-nd(i)
		write(*,*)nf
		pdt=0d0
		pdd=0d0
		!call selfconsist_tg(0.01d0,pdt,pdd,sp)
		!call band(pdd,pdt,sp)
		!call mingap(pdd,pdt,sp,0d0,gap,mx,my)
		!call fermisurface(pdd,pdt,sp,0d0)
		!write(*,*)mx,my,gap
		!!write(80,"(F5.2)")1d0-nf
		!do j=0,600,1
			!Tk=0.001+j
			!dt=dt+0.1
			!dd=dd+0.1
			!call selfconsist_tg(Tk,dt,dd,sp)
			!if(dt(1)<cvg*100d0) then
				!Tc=Tk
				!exit
			!endif
		!enddo
		do j=0,100,2
		!do j=1,3
			Tk=0.001+j
			!Tk=Td(j)+0.001
			dt=dt+0.1
			dd=dd+0.1
			!call selfconsist(Tk,nf,dt,dd,sp)
			call selfconsist_tg(Tk,dt,dd,sp)
			!if((dt(1)-cvg*100)*(pdt(1)-cvg*100)<0d0) then
				!write(60,"('0',2e12.4)")nf,Tk
			!endif
			!if((dd-cvg*100)*(pdd-cvg*100)<0d0) then
				!write(60,"('1',2e12.4)")nf,Tk
			!endif
			!write(50,"(F5.0)")Tk
			!call raman(dd,dt,sp,Tk,(/0d0,0.2d0/),0.0005d0,pk)
			!if(j==0) then
				!pk0=pk
			!endif
			!write(80,"(4e16.3)")Tk/Tc,pk/pk0
			!call band(dd,dt,sp)
			call mingap(dd,dt,sp,0d0,gap,mx,my)
			!write(*,"(3e12.4)")gap,mx,my
			!write(30,"(F5.0)")Tk
			!call EDC(pi,0d0,dd,dt,sp,Tk,(/-0.1d0,0d0/),0.0001d0,peak)
			!call EDC(mx,my,dd,dt,sp,Tk,(/-0.1d0,0d0/),0.0001d0,peak)
			!write(40,"(6e12.4)")nf,Tk,sp,dt(1),dd,peak(1)
			write(40,"(6e12.4)")nf,Tk,sp,dt(1),dd,gap
			!call fermisurface(dd,dt,sp,0d0)
			!if(dt(1)<cvg*100d0.and.dd<cvg*100d0) then
			if(dt(1)<cvg*100d0) then
				exit
			endif
			pdt=dt
			pdd=dd
		enddo
		write(10,"(1X)")
		write(30,"(1X)")
		write(40,"(1X/)")
		write(80,"(1X/)")
	enddo
	close(10)
	close(20)
	close(30)
	close(40)
	close(50)
	close(60)
end
subroutine raman(dd,dt,sp,Tk,omgr,domg,pk)
	use global
	implicit none
	integer :: i,j,k,l,m1,m2
	complex(8) :: Uk(4,4),R(2)
	real(8) :: kx,ky,sp,bt,dt(2),dd,ek(4),gm(4,2),Tr(2),fk(4),omg,omgr(2),domg,Tk,peak(3,2),R_rpa,pk(3)
	bt=escal/Tk*1.16e4
	omg=omgr(1)
	peak=0d0
	do while(omg<omgr(2))
		omg=omg+domg
		R=0d0
		!$OMP PARALLEL DO REDUCTION(+:R) PRIVATE(kx,ky,gm,ek,Uk,fk,Tr) SCHEDULE(GUIDED)
		do i=0,mr
			kx=pi/mr*i
			do j=0,mr-i
				ky=pi/mr*j
				gm(1:2,1)=(cos(kx)-cos(ky))*((1d0-nf)*t(1)+ap*DJ)*(/0.5d0,-0.5d0/)
				gm(1:2,2)=sin(kx)*sin(ky)*(1d0-nf)*t(2)*(/-2d0,-2d0/)
				gm(3:4,:)=-gm(1:2,:)
				call EU(kx,ky,dd,dt,sp,ek,Uk)
				fk=1d0/(1d0+exp(bt*ek))
				do m1=1,4
					do m2=1,4
						Tr=(/abs(sum(gm(:,1)*dconjg(Uk(:,m2))*Uk(:,m1)))**2,abs(sum(gm(:,2)*dconjg(Uk(:,m2))*Uk(:,m1)))**2/)
						R=R+Tr*(fk(m1)-fk(m2))*1d0/(omg+ek(m1)-ek(m2)+img*domg*40d0)
					enddo
				enddo
			enddo
		enddo
		!$OMP END PARALLEL DO
		R=1d0*R/mr**2
		R_rpa=-dimag(R(1)/(1d0+DJ/((1d0-nf)*t(1)+DJ*ap)**2*R(1)))
		write(50,"(4e16.3)")omg,-dimag(R),R_rpa
		if(peak(1,2)<-dimag(R(1))) then
			peak(1,:)=(/omg,-dimag(R(1))/)
		endif
		if(peak(2,2)<-dimag(R(2))) then
			peak(2,:)=(/omg,-dimag(R(2))/)
		endif
		if(peak(3,2)<R_rpa) then
			peak(3,:)=(/omg,R_rpa/)
		endif
	enddo
	pk=peak(:,1)
	write(50,"(1X/)")
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
		!kx=
		!kqx=
		!do j
			!ky=
			!kqy=
			!ek=
			!ekq=
			!fk=
			!fkq=
			!g=2d0*r*(sin(kx-qx/2d0)-sin(ky-qy/2d0))
			!B=B+g*(fk-fkq)/(omg+ek-ekq+img*0.01d0)
		!enddo
	!enddo
!end
subroutine EDC(kx,ky,dd,dt,sp,Tk,omgr,domg,peak)
	use global
	implicit none
	integer :: i
	complex(8) :: Uk(4,4)
	real(8) :: kx,ky,sp,bt,dt(2),dd,ek(4),fk(4),A,omg,omgr(2),domg,Tk,peak(2)
	bt=escal/Tk*1.16e4
	omg=omgr(1)
	call EU(kx,ky,dd,dt,sp,ek,Uk)
	fk=1d0/(1d0+exp(bt*ek))
	peak=0d0
	!fk=1d0
	do while(omg<omgr(2))
		omg=omg+domg
		A=-sum(DIMAG(fk*Uk(1,:)*dconjg(Uk(1,:))/(omg-ek+img*domg*10d0)))
		if(peak(2)<A) then
			peak=(/-omg,A/)
		endif
		write(30,"(5e16.3)")-omg,A,kx,ky
	enddo
	write(30,"(1X)")
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
		!write(10,"(2e16.3)")(ek(l),abs(Uk(1,l)),l=4,4)
	enddo
end
subroutine mingap(dd,dt,sp,th,gap,mx,my)
	use global
	implicit none
	integer :: i,l
	complex(8) :: Uk(4,4)
	real(8) :: kx,ky,mx,my,th,ek(4),sp,dd,dt(2),gap
	gap=1000d0
	do i=0,mth
		ky=pi-i*pi/mth
		kx=pi-i*pi/mth*tan(th/180d0*pi)
		call EU(kx,ky,dd,dt,sp,ek,Uk)
		!do l=1,4
		do l=1,1
			if(abs(ek(l))<gap.and.ek(l)<0d0.and.abs(Uk(1,l))**2>cvg) then
				gap=abs(ek(l))
				mx=kx
				my=ky
			endif
		enddo
	enddo
end
subroutine fermisurface(dd,dt,sp,omg)
	use global
	implicit none
	integer :: i,j,l
	complex(8) :: Uk(4,4)
	real(8) :: kx,ky,ek(4),sp,dd,dt(2),omg,A
	do i=0,mf
		do j=0,mf
			kx=pi/mf*i
			ky=pi/mf*j
			call EU(kx,ky,dd,dt,sp,ek,Uk)
			A=0d0
			do l=1,4
				A=A-dimag(abs(Uk(1,l))**2/(omg+img*0.01d0-ek(l)))
			enddo
			write(20,"(3e16.3)")kx,ky,A
		enddo
		write(20,"(1X)")
	enddo
	write(20,"(1X/)")
end
subroutine selfconsist_tg(Tk,dt,dd,sp)
	use global
	implicit none
	real(8) :: Tk
	complex(8) :: Uk(4,4)
	real(8) :: kx,ky,n1,dd,ddp,ddk,ek(4),sp,sa,sb,sp0,dk,dt(2),dtp(2),wide,cvg1
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
			call sfcsap(Tk,sp,ap)
			!$OMP PARALLEL DO REDUCTION(+:n1,ddp,dtp) PRIVATE(kx,ky,ek,Uk) SCHEDULE(GUIDED)
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
			ddp=DJ*ddp/(mk**2)*2d0
			dtp=V*dtp/(mk**2)*2d0
			!write(*,"(4(a6,e12.4))")"sa=",sa,",sp=",sp,",sb=",sb,"n=",n1
			!write(*,*)flaga,flagb,nf
			if(abs(n1-nf)<=cvg1) then
				exit
			endif
			if(n1<nf) then
				flaga=.false.
				sa=sp
				if(flagb) then
					sb=sp+wide
				endif
			else
				flagb=.false.
				sb=sp
				if(flaga) then
					sa=sp-wide
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
		!write(*,"(5(a6,e12.4)a6i4)")"n="-,n1,",Tk=",Tk,",sp=",sp,",DSC=",dt(1),",DDW=",dd,"num:",c
	enddo
	!write(*,*)"!!!!!!selfconsist return!!!!!!!!"
	write(*,"(5e12.4,i4)")n1,Tk,sp,dt(1),dd,c
	!write(*,*)"!!!!!!!!!!!!end!!!!!!!!!!!!!!!!!"
end
subroutine EU(kx,ky,dd,dt,sp,ek,Uk)
	use global
	implicit none
	complex(8) :: Uk(4,4),cth,sth,cfy(2),sfy(2)
	real(8) :: eka,eks,e1(2),e2(2),ek(4),kx,ky,dd,ddk,sp,dk,dt(2),gk(2),th,fy(2),cos2th,sin2th,cos2fy(2),sin2fy(2)
	eks=-2d0*(1d0-nf)*t(2)*cos(kx)*cos(ky)-sp
	eka=-((1d0-nf)*t(1)+ap*DJ)*(cos(kx)+cos(ky))
	gk(1)=0.5d0*(cos(kx)-cos(ky))
	! gk(2)=0.5d0*(cos(3d0*kx)-cos(3d0*ky))
	ddk=dd*gk(1)
	dk=dt(1)*gk(1)
	e1=eks+(/1d0,-1d0/)*sqrt(ddk**2+eka**2)*sign(1d0,eka)
	e2=sqrt(e1**2+dk**2)*sign(1d0,e1)
	ek=(/e2,-e2/)
	cos2th=eka/sqrt(ddk**2+eka**2)*sign(1d0,eka)
	sin2th=ddk/sqrt(ddk**2+eka**2)*sign(1d0,eka)
	cos2fy=e1/e2*sign(1d0,e1)*sign(1d0,e1)
	sin2fy=dk/e2*sign(1d0,e1)*sign(1d0,e1)
	cth=dcmplx(sqrt(0.5d0*(1d0+cos2th)))
	sth=dcmplx(sqrt(0.5d0*(1d0-cos2th))*sign(1d0,sin2th))
	cfy=dcmplx(sqrt(0.5d0*(1d0+cos2fy)))
	sfy=dcmplx(sqrt(0.5d0*(1d0-cos2fy))*sign(1d0,sin2fy))
	Uk=reshape((/   cth*cfy(1) , -img*sth*cfy(1) ,      cth*sfy(1) ,    img*sth*sfy(1) ,   & 
	           -img*sth*cfy(2) ,      cth*cfy(2) , -img*sth*sfy(2) ,       -cth*sfy(2) ,   & 
	               -cth*sfy(1) ,  img*sth*sfy(1) ,      cth*cfy(1) ,    img*sth*cfy(1) ,   & 
	           -img*sth*sfy(2) ,      cth*sfy(2) ,  img*sth*cfy(2) , cth*cfy(2)     /) , (/4 , 4/))
end
subroutine sfcsap(Tk,sp,ap)
	use global, only : mk,cvg,escal,pi,t,DJ,nf
	implicit none
	real(8) :: kx,ky,Tk,bt,sp,ap,app
	integer :: i,j
	bt=escal/Tk*1.16e4
	do
		app=0d0
		!$OMP PARALLEL DO REDUCTION(+:app) PRIVATE(kx,ky) SCHEDULE(GUIDED)
		do i=0,mk
			do j=0,mk
				kx=pi/mk*i
				ky=pi/mk*j
				app=app+cos(kx)*1d0/(1d0+exp(bt*(-((1d0-nf)*t(1)+ap*DJ)*(cos(kx)+cos(ky))-&
					2d0*(1d0-nf)*t(2)*cos(kx)*cos(ky)-sp)))
			enddo
		enddo
		!$OMP END PARALLEL DO
		app=app/mk**2
		if(abs(app-ap)<cvg) then
			exit
		endif
		ap=app
	enddo
end
subroutine order(kx,ky,ek,Uk,Tk,n,dd,dt)
	use global
	implicit none
	complex(8) :: Uk(4,4)
	real(8) :: kx,ky,ek(4),dd,dt(2),gk(2),n,Tk,fk(4),bt
	bt=escal/Tk*1.16e4
	fk=1d0/(1d0+exp(bt*ek))
	gk(1)=0.5d0*(cos(kx)-cos(ky))
	dd=dd+dimag(dot_product(Uk(1,:)*dconjg(Uk(2,:))-Uk(2,:)*dconjg(Uk(1,:))-&
		Uk(4,:)*dconjg(Uk(3,:))+Uk(3,:)*dconjg(Uk(4,:)),fk)*gk(1))
	n=n+2d0+dot_product(Uk(1,:)*dconjg(Uk(1,:))+Uk(2,:)*dconjg(Uk(2,:))-&
		Uk(3,:)*dconjg(Uk(3,:))-Uk(4,:)*dconjg(Uk(4,:)),fk)
	dt=dt+dot_product(Uk(1,:)*dconjg(Uk(3,:))-Uk(2,:)*dconjg(Uk(4,:)),fk)*gk(1)
end
subroutine gnuplot()
	implicit none
	!plot energy band
	write(10,"(A)")"set term pngcairo"
	write(10,"(A)")"set output 'energy.png'"
	write(10,"(A)")"unset key"
	write(10,"(A)")"set size square"
	write(10,"(A)")"set xlabel 'k'"
	write(10,"(A)")"set ylabel '能量'"
	write(10,"(A)")'set label "n=0.8075\nT=0K" at 5,0.8'
	write(10,"(A)")"plot [:][-1:1] for[i=1:4] '-' using 0:2*i-1:2*i with points lt i pt 7 ps variable"
	write(10,"(A)")"#data"
	!plot fermi surface
	write(20,"(A)")"set term pdf fontscale 0.7 transparent enhanced size 15,5"
	write(20,"(A)")"set output 'fermi.pdf'"
	write(20,"(A)")"set tmargin 0"
	write(20,"(A)")"set bmargin 0"
	write(20,"(A)")"set lmargin 0"
	write(20,"(A)")"set rmargin 0"
	write(20,"(A)")"set multiplot"
	write(20,"(A)")"unset key"
	write(20,"(A)")"set palette rgbformulae 22,13,-31"
	write(20,"(A)")"set size 0.28,1"
	write(20,"(A)")"#set cbrange [0:1]"
	write(20,"(A)")"set pm3d map"
	write(20,"(A)")"set pm3d interpolate 0,0"
	write(20,"(A)")"do for[i=0:2]{"
	write(20,"(A)")"	set origin 0.1+0.28*i,0"
	write(20,"(A)")"	splot [0:3.14][0:3.14] '-' index i"
	write(20,"(A)")"}"
	write(20,"(A)")"#data"
	!plot EDC
	write(30,"(A)")"set term pdf fontscale 0.7 transparent enhanced size 5,5"
	write(30,"(A)")"set output 'EDC.pdf'"
	write(30,"(A)")"set tmargin 0"
	write(30,"(A)")"set bmargin 0"
	write(30,"(A)")"set multiplot"
	write(30,"(A)")"set key autotitle columnhead"
	write(30,"(A)")"set xtics 0,0.05,1 out nomirror"
	write(30,"(A)")"set mxtics 5"
	write(30,"(A)")"set ytics 0,400,1000 out nomirror"
	write(30,"(A)")"set mytics 5"
	write(30,"(A)")"set size 1,0.28"
	write(30,"(A)")"do for[j=0:2] {"
	write(30,"(A)")"	set origin 0,0.14+0.28*j"
	write(30,"(A)")"	if (j>0) {"
	write(30,"(A)")"		set format x ''"
	write(30,"(A)")"	}"
	write(30,"(A)")"	plot [:0.09][:] for[i=0:100] '-' index j every ::1:i::i u 1:($2+i*100) with lines"
	write(30,"(A)")"}"
	write(30,"(A)")"#data"
	!plot Raman
	write(50,"(A)")"set term pngcairo size 1200,800"
	write(50,"(A)")"set termoption enhanced"
	write(50,"(A)")"set output 'raman.png'"
	write(50,"(A)")"set title 'n=0.855 T_c='"
	write(50,"(A)")"set multiplot layout 3,2"
	write(50,"(A)")"set key autotitle columnhead"
	write(50,"(A)")"plot [:0.2][:] for[i=0:130:5] '-' index i using 1:2 with lp pi 5 ps 0.7"
	write(50,"(A)")"plot [:0.2][:] for[i=0:130:5] '' index i using 1:($2-$5) with lp pi 5 ps 0.7"
	write(50,"(A)")"plot [:0.2][:] for[i=0:130:5] '' index i using 1:4 with lp pi 5 ps 0.7"
	write(50,"(A)")"plot [:0.2][:] for[i=0:130:5] '' index i using 1:($4-$7) with lp pi 5 ps 0.7"
	write(50,"(A)")"plot [:0.2][:] for[i=0:130:5] '' index i using 1:3 with lp pi 5 ps 0.7"
	write(50,"(A)")"plot [:0.2][:] for[i=0:130:5] '' index i using 1:($3-$6) with lp pi 5 ps 0.7"
	write(50,"(A)")"#data"
	!plot Raman temp
	write(80,"(A)")"set term pngcairo size 600,800"
	write(80,"(A)")"set termoption enhanced"
	write(80,"(A)")"set output 'ramantemp.png'"
	write(80,"(A)")"set multiplot layout 3,1"
	write(80,"(A)")"set key autotitle columnhead"
	write(80,"(A)")"plot [:][:] for[i=0:130:5] '-' index i using 1:2 with lp pi 5 ps 0.7"
	write(80,"(A)")"plot [:][:] for[i=0:130:5] '-' index i using 1:3 with lp pi 5 ps 0.7"
	write(80,"(A)")"plot [:][:] for[i=0:130:5] '-' index i using 1:4 with lp pi 5 ps 0.7"
	write(80,"(A)")"#data"
	!plot gap_tmp
	write(40,"(A)")"reset"
	write(40,"(A)")"set term pdf fontscale 0.7 transparent enhanced size 7,5"
	write(40,"(A)")"set output 'gaptmp.pdf'"
	write(40,"(A)")"set termoption enhanced"
	write(40,"(A)")"set lmargin 0"
	write(40,"(A)")"set rmargin 0"
	write(40,"(A)")"set bmargin 6"
	write(40,"(A)")"set tmargin 1"
	write(40,"(A)")"set multiplot"
	write(40,"(A)")"set xtics 0,0.1,2 out nomirror"
	write(40,"(A)")"set mxtics 5"
	write(40,"(A)")"unset xlabel"
	write(40,"(A)")"set ytics 0,0.03,1 out nomirror"
	write(40,"(A)")"set mytics 5"
	write(40,"(A)")"unset ylabel"
	write(40,"(A)")"set key at screen 0.95,0.95"
	write(40,"(A)")"set size 0.28,1"
	write(40,"(A)")"set label 'gap size' rotate by 90 at screen 0.03,0.5"
	write(40,"(A)")"set label '$\delta$' at screen 0.5,0.05"
	write(40,"(A)")"i=0"
	write(40,"(A)")"do for[name in '(c) (b) (a)'] {"
	write(40,"(A)")"	if(i>0){"
	write(40,"(A)")"		set format y ''"
	write(40,"(A)")"		unset key"
	write(40,"(A)")"	}"
	write(40,"(A)")"	set label name at graph 0.05,0.9"
	write(40,"(A)")"	set origin 0.14+i*0.28,0"
	write(40,"(A)")"	plot [:][0:0.1] '-' index i using 2:4 with lines lw 5 title 'dSC', '' index i using 2:5 with lines lw 5 title 'dCDW', '' index i using 2:6 with lines lw 5 title 'gap'"
	write(40,"(A)")"	unset label"
	write(40,"(A)")"	i=i+1"
	write(40,"(A)")"}"
	write(40,"(A)")"#data"
	!plot phase diagram
	write(60,"(A)")"reset"
	write(60,"(A)")"set term pngcairo"
	write(60,"(A)")"set output 'phase.png'"
	write(60,"(A)")"set termoption enhanced"
	write(60,"(A)")"set xlabel '{/Symbol d}'"
	write(60,"(A)")"set ylabel 'T(K)'"
	write(60,"(A)")'#set label "label" at graph 0.61,0.61'
	write(60,"(A)")"#set style rect fc lt -1 fs solid 0.15 noborder"
	write(60,"(A)")"#set obj rect from 0.13, graph 0 to 0.15, graph 1"
	write(60,"(A)")"plot [][0:70] '-' using (1-$2):($1==1?$3:1/0) with points pt 7 title 'DDW', '' using (1-$2):($1==0?$3:1/0) &
		with points pt 7 title 'SC'"
	write(60,"(A)")"#data"
	!plot gap_theta
	write(70,"(A)")"reset"
	write(70,"(A)")"set term pngcairo"
	write(70,"(A)")"set output 'gaptheta.png'"
	write(70,"(A)")"set termoption enhanced"
	write(70,"(A)")"unset key"
	write(70,"(A)")"set xtic 5"
	write(70,"(A)")"set xlabel '\theta'"
	write(70,"(A)")"set ylabel 'gap/t'"
	write(70,"(A)")"set title '角度依赖'"
	write(70,"(A)")'set label "label" at 150,0.11'
	write(70,"(A)")"plot '-' using 1:2 with linespoints pt 7 lw 2"
	write(70,"(A)")"#data"
end
