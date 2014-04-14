module global
	implicit none
	save
	real(8), parameter ::t(5)=(/2d0,-0.9d0,0.000d0,0.d0,0d0/),escal=0.25d0,pi=3.1415926d0,cvg=1e-5,DJ=0.75d0,V=-3d0
	real(8) :: nf=0.8d0
	integer, parameter :: mk=1024,ms=200,mf=128,mr=1024
	complex(8),parameter :: img=(0d0,1d0)
end module
program main
	use global, only : nf,pi,cvg,t,DJ,V
	implicit none
	complex(8) :: Uk(2,2)
	real(8) :: kf(2,2),sp=0d0,bt,dt(2),pdt(2),ap=0.1d0,th,gap,Tk,Tc,pk(3),pk0(3),peak(2),&
		nd(3)=(/0.0025d0,0.114d0,0.178d0/),Td(3)=(/0d0,40d0,70d0/)
	integer :: i,j,l
	call fileopen
	call gnuplot
	!do i=24,96,1
	!do i=24,72,8
	!nf=0.94d0
	nf=0.98d0
	!do i=1,3
	do
		!nf=1d0-nd(i)
		pdt=0d0
		Tk=0.0001d0
		do
		!do j=1,3
			!Tk=Td(j)+0.001
			dt=dt+0.1
			call selfconsist_tg(Tk,dt,ap,sp)
			call spinrpa(dt,ap,sp,Tk,(/pi,pi/),(/-0.05d0,1d0/),0.001d0)
			write(40,"(4e12.4)")nf,Tk,sp,dt(1)
			!call band(dt,ap,sp,Tk,(/0d0,0d0/),(/pi,pi/))
			!call band(dt,ap,sp,Tk,(/pi,pi/),(/pi,0d0/))
			!call band(dt,ap,sp,Tk,(/pi,0d0/),(/0d0,0d0/))
			!if(dt(1)<cvg*100d0.and.dd<cvg*100d0) then
			!!if(dt(1)<cvg*100d0) then
				!exit
			!endif
			if(Tk>0d0) then
				exit
			endif
			Tk=Tk+1d0
			pdt=dt
			write(30,"(1X)")
		enddo
		nf=nf-0.02d0
		if(nf<0.84d0) then
			exit
		endif
		write(10,"(1X)")
		write(50,"(1X)")
		write(40,"(1X/)")
		write(70,"(1X)")
		write(80,"(1X)")
	enddo
	call fileclose
end
subroutine spinrpa(dt,ap,sp,Tk,q,omgr,domg)
	use global
	implicit none
	integer :: i,j,l,m,n
	complex(8) :: Uk(2,2),Ukq(2,2),Xq,Xq_rpa
	real(8) :: k(2),q(2),sp,bt,dt(2),ap,ek(2),ekq(2),fk(2),fkq(2),omg,omgr(2),domg,Tk,DJq
	bt=escal/Tk*1.16e4
	omg=omgr(1)
	DJq=(cos(q(1))+cos(q(2)))*0.34d0
	do while(omg<omgr(2))
		omg=omg+domg
		Xq=0d0
		!$OMP PARALLEL DO REDUCTION(+:Xq) PRIVATE(k,ek,Uk,fk,ekq,Ukq,fkq) SCHEDULE(GUIDED)
		do i=0,mr
			k(1)=pi/mr*i
			do j=0,mr
				k(2)=pi/mr*j
				call EU(k,dt,ap,sp,ek,Uk)
				call EU(k+q,dt,ap,sp,ekq,Ukq)
				fk=1d0/(1d0+exp(bt*ek))
				fkq=1d0/(1d0+exp(bt*ekq))
				do n=1,2
					do m=1,2
						Xq=Xq+(Uk(1,n)*Ukq(2,m)*dconjg(Uk(1,n)*Ukq(2,m))-Uk(2,n)*Ukq(1,m)*dconjg(Uk(1,n)*Ukq(2,m)))*&
							(1-fk(n)-fkq(m))/(omg+ek(n)+ekq(m)+img*0.001d0)
					enddo
				enddo
			enddo
		enddo
		!$OMP END PARALLEL DO
		Xq=Xq/mr**2
		Xq_rpa=Xq/(1d0+DJq*Xq)
		write(50,"(5e16.3)")omg,dimag(Xq),dreal(Xq),(1d0+DJq*dreal(Xq)),dimag(Xq_rpa)
	enddo
	write(50,"(1X)")
end
subroutine band(dt,ap,sp,Tk,ki,kf)
	use global
	implicit none
	integer :: i,l
	complex(8) :: Uk(2,2),Ukq(2,2)
	real(8) :: ki(2),kf(2),kn(2),ek(2),ekq(2),sp,dt(2),Tk,ap
	kn=ki
	i=0
	do
		kn=kn+(kf-ki)/sqrt(sum((kf-ki)**2))*pi/ms
		call EU(kn,dt,ap,sp,ek,Uk)
		call EU(kn+(/pi,pi/),dt,ap,sp,ekq,Ukq)
		!export band
		!write(10,"(10e16.3)")(ek(l),abs(Uk(1,l)),l=1,4)
		write(10,"(12e16.3)")(ek(l),abs(Uk(1,l)),l=1,2),(ekq(l),abs(Ukq(1,l)),l=1,2),(ek(l)+ekq(l),1.0,l=1,2)
		if(all((kn-kf)*(kf-ki)>=0d0)) then
			exit
		endif
		i=i+1
	enddo
	write(*,*)i
end
subroutine selfconsist_tg(Tk,dt,ap,sp)
	use global
	implicit none
	real(8) :: Tk
	complex(8) :: Uk(2,2)
	real(8) :: k(2),n1,ap,app,ek(2),sp,sa,sb,sp0,dk,dt(2),dtp(2),wide,cvg1
	integer :: i,j,c,info
	logical :: flaga,flagb
	! dtp(2)=-2d0
	wide=0.5d0
	sp0=sp+wide
	c=0
	cvg1=0.01
	do 
		sa=sp
		sb=sp
		flaga=.true.
		flagb=.true.
		do 			
			sp=0.5d0*(sa+sb)
			dtp=0d0
			app=0d0
			n1=0d0
			c=c+1
			!call sfcsap(Tk,sp,ap)
			!$OMP PARALLEL DO REDUCTION(+:n1,dtp,app) PRIVATE(k,ek,Uk) SCHEDULE(GUIDED)
			do i=0,mk
				do j=0,mk
					k=(/pi/mk*i,pi/mk*j/)
					call EU(k,dt,ap,sp,ek,Uk)
					call order(k,ek,Uk,Tk,n1,dtp,app)
				enddo
			enddo
			!$OMP END PARALLEL DO
			n1=n1/(mk**2)
			dtp=V*dtp/(mk**2)
			app=app/(mk**2)*DJ
			!write(*,"(4(a6,e12.4))")"sa=",sa,",sp=",sp,",sb=",sb,"n=",n1
			!write(*,*)abs(n1-nf),c
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
		cvg1=max(cvg,min(cvg1,0.1d0*(abs(dtp(1)-dt(1)))))
		!write(*,"(5(a6,e12.4)a6i4)")"n=",n1,",Tk=",Tk,",sp=",sp,",DSC=",dt(1),",DDW=",dd,"num:",c
		if((abs(dtp(1)-dt(1))+abs(ap-app))<cvg) then
			exit
		endif
		dt=dtp
		ap=app
	enddo
	!write(*,*)"!!!!!!selfconsist return!!!!!!!!"
	write(*,"(5e12.4,i4)")n1,Tk,sp,dt(1),ap,c
	!write(*,*)"!!!!!!!!!!!!end!!!!!!!!!!!!!!!!!"
end
subroutine EU(k,dt,ap,sp,ek,Uk)
	use global
	implicit none
	complex(8) :: Uk(2,2),cfy,sfy
	real(8) :: eka,eks,ek(2),k(2),sp,dk,dt(2),ap,gk(2),th,fy(2),cos2fy,sin2fy
	eks=-4d0*(1d0-nf)*t(2)*cos(k(1))*cos(k(2))-sp-2d0*(1d0-nf)*t(3)*(cos(2d0*k(1))+cos(2d0*k(2)))
	eka=-2d0*((1d0-nf)*t(1)+ap)*(cos(k(1))+cos(k(2)))
	gk(1)=0.5d0*(cos(k(1))-cos(k(2)))
	dk=dt(1)*gk(1)
	ek=sqrt((eka+eks)**2+dk**2)*(/1d0,-1d0/)
	cos2fy=(eka+eks)/sqrt((eka+eks)**2+dk**2)
	sin2fy=dk/sqrt((eka+eks)**2+dk**2)
	cfy=dcmplx(sqrt(0.5d0*(1d0+cos2fy)))
	sfy=dcmplx(sqrt(0.5d0*(1d0-cos2fy))*sign(1d0,sin2fy))
	Uk=reshape((/   cfy , sfy , & 
				   -sfy , cfy /) , (/2 , 2/))
end
subroutine order(k,ek,Uk,Tk,n,dt,ap)
	use global
	implicit none
	complex(8) :: Uk(2,2)
	real(8) :: k(2),ek(2),dt(2),gk(2),n,Tk,fk(2),bt,ap
	bt=escal/Tk*1.16e4
	fk=1d0/(1d0+exp(bt*ek))
	gk(1)=0.5d0*(cos(k(1))-cos(k(2)))
	n=n+1d0+dot_product(Uk(1,:)*dconjg(Uk(1,:))-Uk(2,:)*dconjg(Uk(2,:)),fk)
	ap=ap+dot_product(Uk(1,:)*dconjg(Uk(1,:)),fk)*cos(k(1))
	dt=dt+dot_product(Uk(1,:)*dconjg(Uk(2,:)),fk)*gk(1)
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
	open(unit=50,file="../data/rpa.dat")
	open(unit=60,file="../data/gaptheta.dat")
	open(unit=70,file="../data/ramantemp.dat")
	open(unit=80,file="../data/phase_tJ.dat")
end
subroutine gnuplot()
	implicit none
	integer :: i
	do i=10,80,10
		write(i,"(A)")'reset'
		write(i,"(A)")'#custom'
		write(i,"(A)")'#multiplot'
		write(i,"(A)")'ix=1'
		write(i,"(A)")'iy=1'
		write(i,"(A)")'sx=5*ix'
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
		write(i,"(A)")'#set term eps font ",18" size sx,sy'
		write(i,"(A)")'#set output "-."."eps"'
		write(i,"(A)")'set label ylb font ",20" rotate by 90 center at character 0.8,screen 0.5*(1+mb-mt)'
		write(i,"(A)")'set label xlb font ",20" center at screen 0.5*(1+ml-mr),character 0.5'
		write(i,"(A)")'set xtics out nomirror'
		write(i,"(A)")'set ytics out nomirror'
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
	write(10,"(A)")'	plot [0:][:] for[j=0:11] "-" index i every 5 using 0:2*j+1:(column(2*j+2)*0.5) with points pt 7 ps variable lc j+1 notitle'
	!plot 
	write(20,"(A)")'	plot [:55][0:0.1] for[j=0:2] "-" index i using 2:j+4 with line title word("SC DDW Total",j+1)'
	!plot phase diagram
	write(30,"(A)")'	plot [:55][0:0.1] for[j=0:2] "-" index i using 2:j+4 with line title word("SC DDW Total",j+1)'
	!plot phase diagram
	write(40,"(A)")'	plot [:55][0:0.1] for[j=0:2] "-" index i using 2:j+4 with line title word("SC DDW Total",j+1)'
	!plot phase diagram
	write(50,"(A)")'	plot [:1.2][0:] for[j=0:8] "-" index j using 1:2 with line title "".sprintf("%1.3f",0.02+0.02*j)'
	!plot phase diagram
	write(60,"(A)")'	plot [:55][0:0.1] for[j=0:2] "-" index i using 2:j+4 with line title word("SC DDW Total",j+1)'
	!plot phase diagram
	write(70,"(A)")'	plot [:1.1][0.2:1.3] for[j=0:12:2] "-" using 1:i:(100.) lc (i/2+1) smooth acsplines notitle'
	!plot phase diagram
	write(80,"(A)")'set style fill transparent solid 0.3 border'
	write(80,"(A)")'set style rect fc lt -1 fs solid 0.15 noborder'
	write(80,"(A)")'set obj rect from 0.14, graph 0 to 0.175, graph 1 behind'
	write(80,"(A)")'plot [0.06:0.24][0:100] "-" using (1-$1):3:(1e7) smooth acsplines with filledcu y1=0 lc 3 lw 2, "-" u (1-$1):4:(1e9) smooth acsplines with filledcu x1=0 lc 1 lw 2, "-" u (1-$1):5:(1e9) smooth acsplines with filledcu x1=0 lc 1 lw 2 '
	do i=10,80,10
		write(i,"(A)")'	unset label'
		write(i,"(A)")'	unset format'
		write(i,"(A)")'}'
		write(i,"(A)")'	unset multiplot'
		write(i,"(A)")'if(GPVAL_TERM eq "qt"){'
		write(i,"(A)")'	pause -1'
		write(i,"(A)")'}'
		write(i,"(A)")'#data'
	enddo
end
