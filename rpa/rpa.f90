module global
	implicit none
	save
	real(8), parameter :: t(5)=(/0.07415d0,-0.01750d0,0.0116d0,0.d0,0d0/),escal=1d0,pi=3.1415926d0,cvg=1e-4,&
		ph=0.1464d0,V=-0.5d0,DJ=0.4d0,al=0.5d0,nf=0.9d0
	integer, parameter :: mk=512,mq=128,ms=100,momg=8
	complex(8),parameter :: img=(0d0,1d0)
end module
program main
	use global
	implicit none
	complex(8) :: H(4,4),Uk(4,4),Ukq(4,4),Xq(2,2)
	integer :: i,j
	real(8) :: qx,qy,Jq,Tk,sp,Xq_RPA,dt(2),omg
	character(30) :: filename
	open(unit=10,file="../data/rpa.dat")
	open(unit=20,file="../data/check.dat")
	call gnuplot()
	dt=0.01d0
	Tk=0.05d0
	call selfconsist_tg(Tk,dt,ph,sp)
	call band(ph,dt,sp)
	omg=0.3d0
	do i=0,mq
		qx=pi*(1d0+0.5d0*i/mq)
		do j=0,i
			qy=pi*(1d0+0.5d0*j/mq)
			call susp(qx,qy,dt,sp,Tk,omg,Xq)
			Jq=al*DJ*(cos(qx)+cos(qy))
			Xq_RPA=dimag((Xq(1,1)-Jq*(Xq(1,1)*Xq(2,2)+Xq(1,2)*Xq(2,1)))/&
				(1d0+Jq*(Xq(1,1)-Xq(2,2))-Jq**2*(Xq(1,1)*Xq(2,2)-Xq(1,2)*Xq(2,1))))
			write(10,"(4e16.3)")qx-pi,qy-pi,dimag(Xq(1,1)),Xq_RPA
		enddo
		write(10,"(1X)")
	enddo
end
subroutine selfconsist_tg(Tk,dt,dd,rsp)
	use global
	implicit none
	complex(8) :: Uk(4,4),cth,sth,cfy(2),sfy(2)
	real(8) :: kx,ky,n1,dd,ddp,ddk,ek(4),rsp,sp=0d0,sa,sb,sp0,dk,dt(2),dtp(2),gk(2),wide,th,fy(2),Tk,cvg1
	integer :: i,j,c
	logical :: flaga,flagb
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
					call order(kx,ky,ek,Uk,Tk,n1,dtp)
				enddo
			enddo
			!$OMP END PARALLEL DO
			n1=n1/(mk**2)*2d0
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
		if(abs(dtp(1)-dt(1))<cvg) then
			exit
		endif
		dt=dtp
		!write(*,"(5(a6,e12.4)a6i4)")"n="-,n1,",Tk=",Tk,",sp=",sp,",DSC=",dt(1),",DDW=",dd,"num:",c
	enddo
	!write(*,*)"!!!!!!selfconsist return!!!!!!!!"
	write(*,"(4e12.4,i4)")n1,Tk,sp,dt(1),c
	!write(*,*)"!!!!!!!!!!!!end!!!!!!!!!!!!!!!!!"
	rsp=sp
end
subroutine EU(kx,ky,dd,dt,sp,ek,Uk)
	use global
	implicit none
	complex(8) :: Uk(4,4),cth,sth,cfy(2),sfy(2),eka
	real(8) :: eks,e1(2),e2(2),ek(4),kx,ky,dd,ddk,sp,dk,dt(2),gk(2),th,fy(2)
	eka=-2d0*t(1)*(1d0+DJ*30/8)*(cos(kx)*exp(img*ph)+cos(ky)*exp(-img*ph))
	eks=-4d0*t(2)*cos(kx)*cos(ky)-2d0*t(3)*(cos(2*kx)+cos(2*ky))-sp
	dk=dt(1)*0.5d0*(cos(kx)-cos(ky))
	e1=eks+(/1d0,-1d0/)*sqrt(dimag(eka)**2+real(eka)**2)
	e2=sqrt(e1**2+dk**2)
	ek=(/e2,-e2/)
	if(dimag(eka)>0) then
		th=acos(real(eka)/sqrt(dimag(eka)**2+real(eka)**2))/2d0
	else
		th=(2d0*pi-acos(real(eka)/sqrt(dimag(eka)**2+real(eka)**2)))/2d0
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
	Uk=reshape((/   cth*cfy(1) , -img*sth*cfy(1) ,      cth*sfy(1) ,    img*sth*sfy(1) ,   & 
			   -img*sth*cfy(2) ,      cth*cfy(2) , -img*sth*sfy(2) ,       -cth*sfy(2) ,   & 
				   -cth*sfy(1) ,  img*sth*sfy(1) ,      cth*cfy(1) ,    img*sth*cfy(1) ,   & 
			   -img*sth*sfy(2) ,      cth*sfy(2) ,  img*sth*cfy(2) , cth*cfy(2)     /) , (/4 , 4/))
end
subroutine order(kx,ky,ek,Uk,Tk,n,dt)
	use global
	implicit none
	complex(8) :: Uk(4,4)
	real(8) :: kx,ky,ek(4),dd,dt(2),gk(2),n,Tk,fk(4),bt
	bt=escal/Tk*1.16e4
	fk=1d0/(1d0+exp(bt*ek))
	gk(1)=0.5d0*(cos(kx)-cos(ky))
	n=n+2d0+dot_product(Uk(1,:)*dconjg(Uk(1,:))+Uk(2,:)*dconjg(Uk(2,:))-&
		Uk(3,:)*dconjg(Uk(3,:))-Uk(4,:)*dconjg(Uk(4,:)),fk)
	dt=dt+dot_product(Uk(1,:)*dconjg(Uk(3,:))-Uk(2,:)*dconjg(Uk(4,:)),fk)*gk(1)
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
		write(20,"(8e16.3)")(ek(l),abs(Uk(1,l)),l=1,4)
	enddo
end
subroutine susp(qx,qy,dt,sp,Tk,omg,Xq)
	use global
	implicit none
	complex(8) :: Uk(4,4),Ukq(4,4),Xq(2,2),Xq_tmp(2,2)
	real(8) :: kqx,kqy,kx,ky,qx,qy,ek(4),ekq(4),fk(4),fkq(4),Tk,bt,sp,omg,dt(2)
	integer :: i,j,n,m
	Xq=0d0
	bt=escal/Tk*1.16e4
	!$OMP PARALLEL DO REDUCTION(+:Xq) PRIVATE(kx,ky,UK,ek,fk,Ukq,ekq,fkq,kqx,kqy,Xq_tmp) SCHEDULE(GUIDED)
	do i=-mk,mk
		kx=pi/mk*i
		do j=-mk,mk
			ky=pi/mk*j
			kqx=kx+qx
			kqy=ky+qy
			call EU(kx,ky,ph,dt,sp,ek,Uk)
			call EU(kqx,kqy,ph,dt,sp,ekq,Ukq)
			fk=1d0/(1d0+exp(bt*ek))
			fkq=1d0/(1d0+exp(bt*ekq))
			do n=1,4
				do m=1,4
					Xq_tmp(1,1)=&
						(Uk(1,n)*Ukq(3,m)*dconjg(Uk(1,n)*Ukq(3,m))+Uk(1,n)*Ukq(3,m)*dconjg(Uk(2,n)*Ukq(4,m))-&
						Uk(1,n)*Ukq(3,m)*dconjg(Uk(3,n)*Ukq(1,m))-Uk(1,n)*Ukq(3,m)*dconjg(Uk(4,n)*Ukq(2,m)))
					Xq_tmp(1,2)=&
						(Ukq(3,m)*Uk(1,n)*dconjg(Ukq(3,m)*Uk(2,n))+Uk(1,n)*Ukq(3,m)*dconjg(Uk(1,n)*Ukq(4,m))-&
						Uk(1,n)*Ukq(3,m)*dconjg(Uk(3,n)*Ukq(2,m))-Uk(1,n)*Ukq(3,m)*dconjg(Uk(4,n)*Ukq(1,m)))
					Xq_tmp(2,1)=dconjg(Xq_tmp(1,2))
					Xq_tmp(2,2)=&
						(Ukq(4,m)*Uk(1,n)*dconjg(Ukq(4,m)*Uk(1,n))+Uk(1,n)*Ukq(4,m)*dconjg(Uk(2,n)*Ukq(3,m))-&
						Uk(1,n)*Ukq(4,m)*dconjg(Uk(3,n)*Ukq(2,m))-Uk(1,n)*Ukq(4,m)*dconjg(Uk(4,n)*Ukq(1,m)))
					Xq(:,:)=Xq(:,:)+(1-fk(n)-fkq(m))/(omg-ek(n)-ekq(m)+img*0.001d0)*Xq_tmp
				enddo
			enddo
		enddo
	enddo
	!$OMP END PARALLEL DO
	Xq=-Xq/((2*mk)**2)
end
subroutine gnuplot()
	implicit none
	!plot fermi surface
	write(10,"(A)")"set term pngcairo"
	write(10,"(A)")"set output 'rpa.png'"
	write(10,"(A)")"set multiplot layout 1,2"
	write(10,"(A)")"unset key"
	write(10,"(A)")"set palette rgbformulae 22,13,-31"
	write(10,"(A)")"set size square"
	write(10,"(A)")"#set cbrange [0:1]"
	write(10,"(A)")"set pm3d map"
	write(10,"(A)")"set pm3d interpolate 0,0"
	write(10,"(A)")"splot [:][:] '-' u 1:2:3, '' u 2:1:3, '' u (-$2):1:3, '' u (-$1):2:3, '' u 1:(-$2):3, '' u 2:(-$1):3, &
		'' u (-$2):(-$1):3, '' u (-$1):(-$2):3"
	write(10,"(A)")"splot [:][:] '-' u 1:2:4, '' u 2:1:4, '' u (-$2):1:4, '' u (-$1):2:4, '' u 1:(-$2):4, '' u 2:(-$1):4, &
		'' u (-$2):(-$1):4, '' u (-$1):(-$2):4"
	write(10,"(A)")"unset multiplot"
	write(10,"(A)")"#data"
	!plot energy band
	write(20,"(A)")"set term pngcairo"
	write(20,"(A)")"set output 'energy.png'"
	write(20,"(A)")"unset key"
	write(20,"(A)")"set size square"
	write(20,"(A)")"set xlabel 'k'"
	write(20,"(A)")"set ylabel '能量'"
	write(20,"(A)")'set label "n=0.8075\nT=0K" at 5,0.8'
	write(20,"(A)")"plot [:][-1:1] for[i=1:4] '-' using 0:2*i-1:2*i with points lt i pt 7 ps variable"
	write(20,"(A)")"#data"
end
