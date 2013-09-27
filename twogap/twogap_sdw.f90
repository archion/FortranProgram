program main
	implicit none
	real(8),parameter :: t=0.25d0,tp=-0.025d0,tpp=0.012d0,tppp=0.035d0,tiv=0d0,nf=0.8d0,&
		pi=3.1415926d0,cvg=1e-4
	real(8) :: te,bt,u=0.318,v=-0.3
	integer,parameter :: mk=256,mth=1000,ms=100,mr=128,mo=720
	complex(8),parameter :: img=(0d0,1d0)
	character(30) :: filename
	real(8) :: eka,eks,e1(2),e2(2),e(4),kx,ky,n,US,USp,sp=-0.2,sa,sb,sp0,dk,dt(2),dtp(2),gk(2),hn(4),wide,tmp,tmp1,tmpa,z,th,gap,&
		u2(2),v2(2),uv(2),cos2,sin2,sc,h(4,4),work(15),gm(4,2),rtmp(4,4,4,2),f(4),R(2),omg,DOS,d,omgp,Rp(2),DOSp,M_omg,DOS1,A(mo)
	integer :: i,j,k,l,p,b,info,m1,m2,Mn
	logical :: flaga,flagb,flag,Tr(4,4)
	eks(kx,ky)=-4d0*tp*cos(kx)*cos(ky)-2d0*tpp*(cos(2d0*kx)+cos(2d0*ky))-4d0*tiv*cos(2d0*kx)*cos(2d0*ky)-sp
	eka(kx,ky)=-2d0*t*(cos(kx)+cos(ky))-4d0*tppp*(cos(2d0*kx)*cos(ky)+cos(kx)*cos(2d0*ky))
	open(unit=20,file="../data/fermi.dat")
	open(unit=30,file="../data/temp.dat")
	open(unit=40,file="../data/gap.dat")
	open(unit=50,file="../data/raman.dat")
	flag=.false.
	Tr=reshape((/ .true.,.false.,.false.,.false.,&
				.false.,.true.,.false.,.false.,&
				.false.,.false.,.true.,.false.,&
				.false.,.false.,.false.,.true. /),(/4,4/))
	sp=-0.128d0
	wide=1d0
	US=0.1d0
	dt=0.07d0
	do p=90,90
		te=p
		bt=1d0/(te*8.6e-5)
		! u=0.375d0+p/20d0
		u=0.375d0
		v=-0.0d0
		sp0=sp
		dtp(1)=dt(1)+100*cvg
		! dtp(2)=-2d0
		USp=US+100*cvg
		do while(abs(dtp(1)-dt(1))>cvg.or.abs(USp-US)>cvg)
			dt=dtp
			US=USp
			sa=sp
			sb=sp
			n=0d0
			flaga=.true.
			flagb=.true.
			! write(*,"(4e16.3)")US,dt(1),sp
			do while(abs(n-nf)>cvg)
				sp=0.5d0*(sa+sb)
				n=0d0
				USp=0d0
				dtp=0d0
				!$OMP PARALLEL DO REDUCTION(+:n,USp,dtp) PRIVATE(kx,ky,gk,dk,u2,v2,uv,cos2,sin2,sc,hn,tmp,tmpa,e1,e2,e)&
				!$OMP SCHEDULE(GUIDED)
				do i=0,mk
					do j=0,min(i,mk-i)
						kx=pi/mk*i
						ky=pi/mk*j
						gk(1)=0.5d0*(cos(kx)-cos(ky))
						! gk(2)=0.5d0*(cos(3d0*kx)-cos(3d0*ky))
						tmpa=eka(kx,ky)
						dk=dt(1)*gk(1) ! +dt(2)*gk(2)
						tmp=sqrt((tmpa)**2+US**2)
						e1=eks(kx,ky)+(/tmp,-tmp/)
						e2=sqrt(e1**2+dk**2)
						e=(/e2,-e2/)
						e=1d0/(1d0+exp(bt*e))
						u2=(1d0+e1/e2)*0.5d0
						v2=(1d0-e1/e2)*0.5d0
						uv=0.5d0*dk/e2
						cos2=(1d0+tmpa/tmp)*0.5d0
						sin2=(1d0-tmpa/tmp)*0.5d0
						sc=-0.5d0*US/tmp
						hn=(/u2(1)*cos2,u2(2)*sin2,v2(1)*cos2,v2(2)*sin2/)+(/u2(1)*sin2,u2(2)*cos2,v2(1)*sin2,v2(2)*cos2/)-&
						   (/v2(1)*cos2,v2(2)*sin2,u2(1)*cos2,u2(2)*sin2/)-(/v2(1)*sin2,v2(2)*cos2,u2(1)*sin2,u2(2)*cos2/)
						n=n+dot_product(hn,e)
						hn=(/u2(1),-u2(2),v2(1),-v2(2)/)+(/-v2(1),v2(2),-u2(1),u2(2)/)
						USp=USp+2d0*u*dot_product(hn*sc,e)
						hn=(/uv(1)*cos2,uv(2)*sin2,-uv(1)*cos2,-uv(2)*sin2/)-(/-uv(1)*sin2,-uv(2)*cos2,uv(1)*sin2,uv(2)*cos2/)
						tmp=dot_product(hn,e)
						dtp(1)=dtp(1)+v*gk(1)*tmp
						! dtp(2)=dtp(2)+v*gk(2)*tmp
						n=2d0+n
					enddo
				enddo
				!$OMP END PARALLEL DO
				n=n/(mk**2)*2
				USp=USp/(mk**2)*2
				dtp=dtp/(mk**2)*2
				if(abs(n-nf)<=cvg) then
					exit
				endif
				if(n<nf) then
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
			! write(*,*)n,USp,dtp
			wide=max(abs(sp0-sp),100*cvg)
			sp0=sp
		enddo
		write(*,"(4e16.3)")te,US,dt(1),sp
		write(30,"(5e16.3)")te,US,dt(1),sp
		! if(flag) then
			! cycle
		! endif
		! flag=.true.
		! if(abs(USp)<=cvg*2) then
			! US=0d0
		! endif
		do i=1,3*ms-1
			write(filename,"(A15,I3.3,A4)")"../data/energy_",p,".dat"
			open(unit=200+p/15,file=filename)
			kx=pi*(min((i/ms)*ms,ms)+(1-i/ms)*mod(i,ms))/ms
			ky=pi*min(max((i-ms),0),3*ms-i)/ms
			gk(1)=0.5d0*(cos(kx)-cos(ky))
			! gk(2)=0.5d0*(cos(3d0*kx)-cos(3d0*ky))
			dk=dt(1)*gk(1) ! +dt(2)*gk(2)
			tmpa=eka(kx,ky)
			tmp=sqrt((tmpa)**2+US**2)
			e1=eks(kx,ky)+(/tmp,-tmp/)
			e2=sqrt(e1**2+dk**2)
			e=(/e2,-e2/)
			u2=(1d0+e1/e2)*0.5d0
			v2=(1d0-e1/e2)*0.5d0
			cos2=(1d0+tmpa/tmp)*0.5d0
			sin2=(1d0-tmpa/tmp)*0.5d0
			hn=(/u2(1)*cos2,u2(2)*sin2,v2(1)*cos2,v2(2)*sin2/)
			write(200+p/15,"(8e16.3)")(e(l),hn(l),l=1,4)
		enddo
		do i=0,6
			write(filename,"(A14,I1,A4)")"../data/spect_",i,".dat"
			open(unit=100+i,file=filename)
			A=0d0
			gap=1000d0
			th=i*pi/4d0/6d0
			! !$OMP PARALLEL DO PRIVATE(kx,ky,gk,dk,u2,v2,cos2,sin2,hn,tmp,tmpa,e1,e2,e) SCHEDULE(GUIDED)
			do j=0,mth
				ky=pi-j*pi/mth
				kx=pi-j*pi/mth*tan(th)
				gk(1)=0.5d0*(cos(kx)-cos(ky))
				! gk(2)=0.5d0*(cos(3d0*kx)-cos(3d0*ky))
				dk=dt(1)*gk(1) ! +dt(2)*gk(2)
				tmpa=eka(kx,ky)
				tmp=sqrt((tmpa)**2+US**2)
				e1=eks(kx,ky)+(/tmp,-tmp/)
				e2=sqrt(e1**2+dk**2)
				e=(/e2,-e2/)
				u2=(1d0+e1/e2)*0.5d0
				v2=(1d0-e1/e2)*0.5d0
				cos2=(1d0+tmpa/tmp)*0.5d0
				sin2=(1d0-tmpa/tmp)*0.5d0
				hn=(/u2(1)*cos2,u2(2)*sin2,v2(1)*cos2,v2(2)*sin2/)
				! write(*,"(4e16.3)")e
				do l=1,4
					do k=1,mo
						A(k)=A(k)-DIMAG(hn(l)/(-0.1d0*k/mo-e(l)+img*0.01d0))
					enddo
					! !$OMP CRITICAL
					if(abs(e(l))<gap.and.e(l)<0.and.hn(l)>cvg) then
						gap=abs(e(l))
					endif
					! !$OMP END CRITICAL
				enddo
			enddo
			! !$OMP END PARALLEL DO
			write(40,"(2e16.3)")th/pi*180,gap
			write(100+i,"(2e16.3)")(-0.1d0*k/mo,A(k)/mth,k=1,mo)
			write(100+i,"(1x)")
		enddo
		write(40,"(1x)")
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
	close(100)
	close(101)
	close(102)
	close(103)
	close(104)
	close(105)
	close(106)
end
	
			
