module global
	implicit none
	save
	real(8), parameter ::t(5)=(/0.25d0,-0.025d0,0.012d0,0d0,0d0/),escal=1d0,pi=3.1415926d0,cvg=1e-5
	integer, parameter :: mk=512,mth=1000,ms=100,mr=128,mo=720
	complex(8),parameter :: img=(0d0,1d0)
end module
program main
	use global, only : pi,cvg,mo
	implicit none
	complex(8) :: Uk(4,4)
	real(8) :: kx,ky,sp,bt,U,V,nf=0.8d0,dt(2),pdt(2),dd,pdd,ek(4),tmp,tmp1,tmpa,z,gap,gm(4,2),rtmp(4,4,4,2),fk(4),R(2),&
		omg,DOS,d,omgp,Rp(2),DOSp,M_omg,DOS1,A(mo),Tk
	character(30) :: filename
	integer :: i,j,k,l,p,b,info,m1,m2,Mn
	logical :: flaga,flagb,flag,Tr(4,4)
	open(unit=10,file="../data/energy.dat")
	open(unit=20,file="../data/fermi.dat")
	open(unit=30,file="../data/temp.dat")
	open(unit=40,file="../data/gap.dat")
	open(unit=50,file="../data/raman.dat")
	open(unit=60,file="../data/fhase.dat")
	!gnuplot command start
	write(40,*)"reset"
	write(40,*)"set term wxt 6"
	write(40,*)"# plot phase diagram"
	write(40,*)"#set term pngcairo font 'AR PL UKai CN,14'"
	write(40,*)"#set output 'precipitation.png'"
	write(40,*)"unset key"
	write(40,*)"plot '-' using 2:4:($8<1e-4?1/0:$8) with labels, phase using 2:($8<1e-4?$4:1/0) with points pt 4"
	write(60,*)"reset"
	write(60,*)"set term wxt 6"
	write(60,*)"# plot phase diagram"
	write(60,*)"#set term pngcairo font 'AR PL UKai CN,14'"
	write(60,*)"#set output 'precipitation.png'"
	write(60,*)"unset key"
	write(60,*)"plot '-' using 2:4:($8<1e-4?1/0:$8) with labels, phase using 2:($8<1e-4?$4:1/0) with points pt 4"
	!gnuplot command end
	flag=.false.
	Tr=reshape((/ .true.,.false.,.false.,.false.,&
				.false.,.true.,.false.,.false.,&
				.false.,.false.,.true.,.false.,&
				.false.,.false.,.false.,.true. /),(/4,4/))

	do p=50,100
		U=0.7d0
		V=-0.0d0
		nf=1d0-0.002d0*p
		pdt=0d0
		pdd=0d0
		do j=0,1000,50
			Tk=0.01+j
			dt=1
			dd=1
			!call selfconsist(U,V,Tk,nf,dt,dd,sp)
			call selfconsist_tg(U,V,Tk,nf,dt,dd,sp)
			if(((dt(1)-cvg*10)*(pdt(1)-cvg*10)<0d0).or.((dd-cvg*10)*(pdd-cvg*10)<0d0)) then
				write(60,"(4(a6,e12.4))")"n=",nf,",Tk=",Tk,",DSC=",dt(1),",DDW=",dd
			endif
			if(dt(1)<cvg*10d0.and.dd<cvg*10d0) then
				exit
			endif
			pdt=dt
			pdd=dd
			!call band(dd,dt,sp)
		enddo
		!write(30,"(3(a5e12.4))")"T=",Tk(p),"DSC=",dt(1),"DDW=",dd

		!do i=0,6
			!write(filename,"(A14,I1,A4)")"../data/spect_",i,".dat"
			!open(unit=100+i,file=filename)
			!A=0d0
			!gap=1000d0
			!th=i*pi/4d0/6d0
			!! !$OMP PARALLEL DO PRIVATE(kx,ky,gk,dk,u2,v2,cos2,sin2,hn,tmp,tmpa,e1,e2,e) SCHEDULE(GUIDED)
			!do j=0,mth
				!ky=pi-j*pi/mth
				!kx=pi-j*pi/mth*tan(th)
				!gk(1)=0.5d0*(cos(kx)-cos(ky))
				!! gk(2)=0.5d0*(cos(3d0*kx)-cos(3d0*ky))
				!dk=dt(1)*gk(1) ! +dt(2)*gk(2)
				!tmpa=eka(kx,ky)
				!tmp=sqrt((tmpa)**2+US**2)
				!e1=eks(kx,ky)+(/tmp,-tmp/)
				!e2=sqrt(e1**2+dk**2)
				!e=(/e2,-e2/)
				!u2=(1d0+e1/e2)*0.5d0
				!v2=(1d0-e1/e2)*0.5d0
				!cos2=(1d0+tmpa/tmp)*0.5d0
				!sin2=(1d0-tmpa/tmp)*0.5d0
				!hn=(/u2(1)*cos2,u2(2)*sin2,v2(1)*cos2,v2(2)*sin2/)
				!! write(*,"(4e16.3)")e
				!do l=1,4
					!do k=1,mo
						!A(k)=A(k)-DIMAG(hn(l)/(-0.1d0*k/mo-e(l)+img*0.01d0))
					!enddo
					!! !$OMP CRITICAL
					!if(abs(e(l))<gap.and.e(l)<0.and.hn(l)>cvg) then
						!gap=abs(e(l))
					!endif
					!! !$OMP END CRITICAL
				!enddo
			!enddo
			!! !$OMP END PARALLEL DO
			!write(40,"(2e16.3)")th/pi*180,gap
			!write(100+i,"(2e16.3)")(-0.1d0*k/mo,A(k)/mth,k=1,mo)
			!write(100+i,"(1x)")
		!enddo
		!write(40,"(1x)")
!  Raman
		! d=0.001d0
		! omgp=0d0
		! Rp=0d0
		! DOSp=0d0
		! M_omg=0.1d0
		! Mn=1000
		! do while(omgp<M_omg)
			! omg=omgp+min(M_omg/Mn,d)
			! R=0d0
			! DOS=0d0
			! DOS1=0d0
			! !$OMP PARALLEL DO REDUCTION(+:R,DOS,DOS1) PRIVATE(kx,ky,gk,dk,gm,h,work,e,info,f,rtmp) SCHEDULE(GUIDED)
			! do i=0,mr
				! kx=pi/mr*i
				! do j=0,mr-i
					! ky=pi/mr*j
					! gk(1)=0.5d0*(cos(kx)-cos(ky))
					! ! gk(2)=0.5d0*(cos(3d0*kx)-cos(3d0*ky))
					! dk=dt(1)*gk(1) ! +dt(2)*gk(2)
					! gm(1:2,1)=(cos(kx)-cos(ky))*(/(t+6d0*tppp+8d0*tpp*cos(ky)+cos(kx)*(8d0*tpp+12d0*tppp*cos(ky))),&
											     ! -(t+6d0*tppp-8d0*tpp*cos(ky)+cos(kx)*(-8d0*tpp+12d0*tppp*cos(ky)))/)
					! gm(1:2,2)=sin(kx)*sin(ky)*(/(-4d0*tp-16d0*tppp*cos(kx)-16d0*tppp*cos(ky)),&
											    ! (-4d0*tp+16d0*tppp*cos(kx)+16d0*tppp*cos(ky))/)
					! gm(3:4,:)=-gm(1:2,:)
					! h=reshape((/ eks(kx,ky)+eka(kx,ky),      -US,         dk,          0d0,           &
									! -US,      eks(kx,ky)-eka(kx,ky),     0d0,          -dk,        &
									! dk,             0d0,           -eks(kx,ky)-eka(kx,ky),    -US,       &
									! 0d0,            -dk,         -US,  -eks(kx,ky)+eka(kx,ky) /),(/4,4/))
					! call dsyev( "v", "u", 4, h, 4, e, work, 12, info )
					! f=1d0/(1d0+exp(bt*e))
					! DOS1=DOS1-DIMAG(1d0/(-omg-eks(kx,ky)-eka(kx,ky)+img*0.0004d0))&
						! -DIMAG(1d0/(-omg-eks(kx,ky)+eka(kx,ky)+img*0.0004d0))
					! do m1=1,4
						! do m2=1,4
							! rtmp(m1,m2,:,1)=h(m1,:)*h(m2,:)
							! rtmp(m1,m2,:,2)=gm(m1,2)*rtmp(m1,m2,:,1)
							! rtmp(m1,m2,:,1)=gm(m1,1)*rtmp(m1,m2,:,1)
						! enddo
					! enddo
					! do m1=1,4
						! DOS=DOS-DIMAG(1d0/(omg-e(m1)+img*0.0004d0))
						! do m2=1,4
							! R=R+(/sum(matmul(rtmp(:,:,m1,1),rtmp(:,:,m2,1)),Tr),sum(matmul(rtmp(:,:,m1,2),rtmp(:,:,m2,2)),Tr)/)*&
								! (-f(m1)+f(m2))*DIMAG(1d0/(omg+e(m1)-e(m2)+img*0.0004d0))
						! enddo
					! enddo
				! enddo
			! enddo
			! !$OMP END PARALLEL DO
			! DOS=2d0*DOS/mr**2
			! R=2d0*R/mr**2
			! d=min(abs(omg-omgp)/abs(1-DOSp/DOS)/10d0,1.5d0*d)
			! DOSp=DOS
			! omgp=omg
			! write(50,"(5e16.3)")omg,R,DOS,DOS1
			! write(*,"(5e16.3)")omg,R,DOS,DOS1
		! enddo
	enddo
	close(10)
	close(20)
	close(30)
	close(40)
	close(50)
	close(60)
end
subroutine selfconsist(U,V,Tk,nf,dt,dd,rsp)
	use global
	implicit none
	real(8) :: bt,Tk,U,V,nf
	complex(8) :: Uk(4,4),cth,sth,cfy(2),sfy(2)
	real(8) :: kx,ky,n1,dd,ddp,ddk,ek(4),rsp,sp=0d0,sa,sb,sp0,dk,dt(2),dtp(2),dta,dtb,gk(2),wide,th,fy(2),fk(4)
	integer :: i,j,c,info
	logical :: flaga,flagb,dtflaga,dtflagb,ddflaga,ddflagb
	! dtp(2)=-2d0
	ddp=dd+100*cvg
	wide=0.5d0
	sp0=sp+wide
	bt=escal/Tk*1.16e4
	c=0
	do while(abs(ddp-dd)>cvg)
		dd=ddp
		dtflaga=.true.
		dtflagb=.true.
		dta=dt(1)
		dtb=dt(1)
		dtp=0d0
		do while(abs(dtp(1)-dt(1))>cvg)
			write(*,*)dta,dt(1),dtb
			write(*,*)dtflaga,dtflagb,dtp(1)
			dt(1)=0.5d0*(dta+dtb)
			flaga=.true.
			flagb=.true.
			sa=sp
			sb=sp
			n1=0d0
			do while(abs(n1-nf)>cvg)
				sp=0.5d0*(sa+sb)
				ddp=0d0
				dtp=0d0
				c=c+1
				!$OMP PARALLEL DO REDUCTION(+:n1,ddp,dtp) PRIVATE(kx,ky,gk,ek,Uk)&
				!$OMP SCHEDULE(GUIDED)
				do i=0,mk
					do j=0,min(i,mk-i)
						kx=pi/mk*i
						ky=pi/mk*j
						call EU(kx,ky,dd,dt,sp,ek,Uk)
						fk=1d0/(1d0+dexp(bt*ek))
						gk(1)=0.5d0*(cos(kx)-cos(ky))
						n1=n1+2d0+dot_product(Uk(1,:)*dconjg(Uk(1,:))+Uk(2,:)*dconjg(Uk(2,:))-&
							Uk(3,:)*dconjg(Uk(3,:))-Uk(4,:)*dconjg(Uk(4,:)),fk)
						ddp=ddp+dimag(dot_product(Uk(1,:)*dconjg(Uk(2,:))-Uk(2,:)*dconjg(Uk(1,:))-&
							Uk(4,:)*dconjg(Uk(3,:))+Uk(3,:)*dconjg(Uk(4,:)),fk)*gk(1))
						dtp=dtp+dot_product(Uk(1,:)*dconjg(Uk(3,:))-Uk(2,:)*dconjg(Uk(4,:)),fk)*gk(1)
					enddo
				enddo
				!$OMP END PARALLEL DO
				n1=n1/(mk**2)*2d0
				ddp=U*ddp/(mk**2)*2d0
				dtp=V*dtp/(mk**2)*2d0
				!write(*,"(4(a6,e12.4))")"sa=",sa,",sp=",sp,",sb=",sb,"n=",n1
				!write(*,*)flaga,flagb,nf
				if(abs(n1-nf)<=cvg) then
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
			if(abs(dtp(1)-dt(1))<=cvg) then
				exit
			endif
			if(dtp(1)>dt(1)) then
				dtflaga=.false.
				dta=dt(1)
				if(dtflaga.or.dtflagb) then
					dtb=dt(1)+0.05
				endif
			else
				dtflagb=.false.
				dtb=dt(1)
				if(dtflaga.or.dtflagb) then
					dta=dt(1)-0.05
				endif
			endif
		enddo
	enddo
	!write(*,*)"!!!!!!selfconsist return!!!!!!!!"
	write(*,"(4(a6,e12.4)a6i4)")"n=",n1,",Tk=",Tk,",DSC=",dt(1),",DDW=",dd,"num:",c
	write(40,"(4(a6,e12.4)a6i4)")"n=",n1,",Tk=",Tk,",DSC=",dt(1),",DDW=",dd,"num:",c
	!write(*,*)"!!!!!!!!!!!!end!!!!!!!!!!!!!!!!!"
	rsp=sp
end
subroutine selfconsist_tg(U,V,Tk,nf,dt,dd,rsp)
	use global
	implicit none
	real(8) :: bt,Tk,U,V,nf
	complex(8) :: Uk(4,4),cth,sth,cfy(2),sfy(2)
	real(8) :: kx,ky,n1,dd,ddp,ddk,ek(4),rsp,sp=0d0,sa,sb,sp0,dk,dt(2),dtp(2),gk(2),wide,th,fy(2),fk(4)
	integer :: i,j,c,info
	logical :: flaga,flagb
	! dtp(2)=-2d0
	wide=0.5d0
	sp0=sp+wide
	bt=escal/Tk*1.16e4
	c=0
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
					fk=1d0/(1d0+exp(bt*ek))
					gk(1)=0.5d0*(cos(kx)-cos(ky))
					n1=n1+2d0+dot_product(Uk(1,:)*dconjg(Uk(1,:))+Uk(2,:)*dconjg(Uk(2,:))-&
						Uk(3,:)*dconjg(Uk(3,:))-Uk(4,:)*dconjg(Uk(4,:)),fk)
					ddp=ddp+dimag(dot_product(Uk(1,:)*dconjg(Uk(2,:))-Uk(2,:)*dconjg(Uk(1,:))-&
						Uk(4,:)*dconjg(Uk(3,:))+Uk(3,:)*dconjg(Uk(4,:)),fk)*gk(1))
					dtp=dtp+dot_product(Uk(1,:)*dconjg(Uk(3,:))-Uk(2,:)*dconjg(Uk(4,:)),fk)*gk(1)
				enddo
			enddo
			!$OMP END PARALLEL DO
			n1=n1/(mk**2)*2d0
			ddp=U*ddp/(mk**2)*2d0
			dtp=V*dtp/(mk**2)*2d0
			!write(*,"(4(a6,e12.4))")"sa=",sa,",sp=",sp,",sb=",sb,"n=",n1
			!write(*,*)flaga,flagb,nf
			!if(.not.flaga.and..not.flagb.and.abs(sa-sb)<=1e-4) then
			if(abs(n1-nf)<=cvg) then
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
		if(abs(dtp(1)-dt(1))<cvg.and.abs(ddp-dd)<cvg) then
			exit
		endif
		dt=dtp
		dd=ddp
	enddo
	!write(*,*)"!!!!!!selfconsist return!!!!!!!!"
	write(*,"(4(a6,e12.4)a6i4)")"n=",n1,",Tk=",Tk,",DSC=",dt(1),",DDW=",dd,"num:",c
	write(40,"(4(a6,e12.4)a6i4)")"n=",n1,",Tk=",Tk,",DSC=",dt(1),",DDW=",dd,"num:",c
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
	ddk=dd*gk(1)
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
	Uk=reshape((/   cth*cfy(1) , -img*sth*cfy(1) ,      cth*sfy(1) ,    img*sth*sfy(1) ,   & 
	           -img*sth*cfy(2) ,      cth*cfy(2) , -img*sth*sfy(2) ,       -cth*sfy(2) ,   & 
	               -cth*sfy(1) ,  img*sth*sfy(1) ,      cth*cfy(1) ,    img*sth*cfy(1) ,   & 
	           -img*sth*sfy(2) ,      cth*sfy(2) ,  img*sth*cfy(2) , cth*cfy(2)     /) , (/4 , 4/))
end
