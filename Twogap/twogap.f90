program main
	implicit none
	real(8),parameter :: t=0.25d0,tp=-0.025d0,tpp=0.012d0,tppp=0.035d0,tiv=0d0,nf=0.8d0,&
		pi=3.1415926d0,cvg=1e-7
	real(8) :: te,bt,u=0.318,v=-0.3
	integer,parameter :: mp1=1024,mp2=10000
	complex(8) :: img=(0d0,1d0)
	real(8) :: eka,eks,e1(2),e2(2),e(4),kx,ky,n,US,USp,sp=-0.2,sa,sb,sp0,dk,dt(2),dtp(2),gk(2),hn(4),wide,tmp,tmp1,tmpa,z,th,gap,&
		u2(2),v2(2),uv(2),cos2,sin2,sc,h(4,4),work(15)
	integer :: i,j,k,l,p,a,b,info
	logical :: flaga,flagb,flag
	eks(kx,ky)=-4d0*tp*cos(kx)*cos(ky)-2d0*tpp*(cos(2d0*kx)+cos(2d0*ky))-4d0*tiv*cos(2d0*kx)*cos(2d0*ky)-sp
	eka(kx,ky)=-2d0*t*(cos(kx)+cos(ky))-4d0*tppp*(cos(2d0*kx)*cos(ky)+cos(kx)*cos(2d0*ky))
	open(unit=10,file="../data/energy.dat")
	open(unit=20,file="../data/fermi.dat")
	open(unit=30,file="../data/temp.dat")
	open(unit=40,file="../data/gap.dat")
	flag=.false.
	do p=245,245
		te=p
		bt=1d0/(te*8.6e-5)
		! u=0.375d0+p/20d0
		u=0.375d0
		v=-0.39d0
		wide=1d0
		sp=-0.135d0
		sp0=sp+wide
		dt=1d0
		US=1d0
		dtp(1)=2d0
		! dtp(2)=-2d0
		USp=2d0
		do while(any(abs(dtp-dt)>cvg).or.abs(US-USp)>cvg)
			! write(*,*)"*"
			dt=dtp
			US=USp
			sa=sp
			sb=sp
			n=0d0
			flaga=.true.
			flagb=.true.
			do while(abs(n-nf)>cvg)
				sp=0.5d0*(sa+sb)
				n=0d0
				USp=0d0
				dtp=0d0
				!$OMP PARALLEL DO REDUCTION(+:n,USp,dtp) PRIVATE(kx,ky,gk,dk,u2,v2,uv,cos2,sin2,sc,hn,tmp,tmpa,e1,e2,e)&
				!$OMP SCHEDULE(GUIDED)
				do i=0,mp1
					do j=0,min(i,mp1-i)
						kx=pi/mp1*i
						ky=pi/mp1*j
						gk(1)=0.5d0*(cos(kx)-cos(ky))
						! gk(2)=0.5d0*(cos(3d0*kx)-cos(3d0*ky))
						tmpa=eka(kx,ky)
						dk=dt(1)*gk(1) ! +dt(2)*gk(2)
						tmp=sqrt((tmpa)**2+US**2)
						e1=eks(kx,ky)+(/-tmp,tmp/)*sign(1d0,tmpa)
						e2=sqrt(e1**2+dk**2)
						e=(/e2*sign(1d0,-e1),-e2*sign(1d0,-e1)/)
						e=1d0/(1d0+exp(bt*e))
						u2=(1d0-abs(e1)/e2)*0.5d0
						v2=(1d0+abs(e1)/e2)*0.5d0
						uv=0.5d0*dk/e2*sign(1d0,-e1)
						cos2=(1d0-abs(tmpa)/tmp)*0.5d0
						sin2=(1d0+abs(tmpa)/tmp)*0.5d0
						sc=0.5d0*US*sign(1d0,tmpa)/tmp
						hn=(/u2(1)*cos2,u2(2)*sin2,v2(1)*cos2,v2(2)*sin2/)+(/u2(1)*sin2,u2(2)*cos2,v2(1)*sin2,v2(2)*cos2/)-&
						   (/v2(1)*cos2,v2(2)*sin2,u2(1)*cos2,u2(2)*sin2/)-(/v2(1)*sin2,v2(2)*cos2,u2(1)*sin2,u2(2)*cos2/)
						n=n+dot_product(hn,e)
						hn=(/u2(1),-u2(2),v2(1),-v2(2)/)+(/-v2(1),v2(2),-u2(1),u2(2)/)
						USp=USp+u*dot_product(hn*sc,e)
						hn=(/uv(1)*cos2,uv(2)*sin2,-uv(1)*cos2,-uv(2)*sin2/)-(/-uv(1)*sin2,-uv(2)*cos2,uv(1)*sin2,uv(2)*cos2/)
						tmp=dot_product(hn,e)
						dtp(1)=dtp(1)+v*gk(1)*tmp
						! dtp(2)=dtp(2)+v*gk(2)*tmp
						n=2d0+n
					enddo
				enddo
				!$OMP END PARALLEL DO
				n=n/(mp1**2)*2d0
				USp=USp/(mp1**2)*2d0
				dtp=dtp/(mp1**2)*2d0
				IF(ABS(n-nf)<=CVG) THEN
					EXIT
				ENDIF
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
		write(*,"(4e16.3)")te,USp,dtp(1),sp
		write(30,"(5e16.3)")te,USp,dtp(1),sp
		! if(flag) then
			! cycle
		! endif
		! flag=.true.
		do i=0,100
			gap=1000d0
			th=i*pi/4d0/100d0
			!$OMP PARALLEL DO PRIVATE(kx,ky,gk,dk,u2,v2,cos2,sin2,hn,tmp,tmpa,e1,e2,e) SCHEDULE(GUIDED)
			do j=0,mp2
				ky=pi-j*pi/mp2
				kx=pi-j*pi/mp2*tan(th)
				gk(1)=0.5d0*(cos(kx)-cos(ky))
				! gk(2)=0.5d0*(cos(3d0*kx)-cos(3d0*ky))
				dk=dt(1)*gk(1) ! +dt(2)*gk(2)
				tmpa=eka(kx,ky)
				tmp=sqrt((tmpa)**2+US**2)
				e1=eks(kx,ky)+(/-tmp,tmp/)*sign(1d0,tmpa)
				e2=sqrt(e1**2+dk**2)
				e=(/e2*sign(1d0,-e1),-e2*sign(1d0,-e1)/)
				u2=(1d0-abs(e1)/e2)*0.5d0
				v2=(1d0+abs(e1)/e2)*0.5d0
				cos2=(1d0-abs(tmpa)/tmp)*0.5d0
				sin2=(1d0+abs(tmpa)/tmp)*0.5d0
				hn=(/u2(1)*cos2,u2(2)*sin2,v2(1)*cos2,v2(2)*sin2/)
				! write(*,"(4e16.3)")e
				do l=1,4
					!$OMP CRITICAL
					if(abs(e(l))<gap.and.e(l)<0.and.hn(l)>1e-4) then
						gap=abs(e(l))
					endif
					!$OMP END CRITICAL
				enddo
			enddo
			!$OMP END PARALLEL DO
			write(40,"(2e16.3)")th/pi*180,gap
		enddo
		write(40,"(1x)")
		i=0
		j=0
		a=1
		b=0
		do while(i>=0)
			if(i==mp1.and.j==0) then
				a=0
				b=1
			endif
			if(i==mp1.and.j==mp1) then
				a=-1
				b=-1
			endif
			i=i+a
			j=j+b
			kx=pi/mp1*i
			ky=pi/mp1*j
			gk(1)=0.5d0*(cos(kx)-cos(ky))
			! gk(2)=0.5d0*(cos(3d0*kx)-cos(3d0*ky))
			dk=dt(1)*gk(1) ! +dt(2)*gk(2)
			tmpa=eka(kx,ky)
			tmp=sqrt((tmpa)**2+US**2)
			e1=eks(kx,ky)+(/-tmp,tmp/)*sign(1d0,tmpa)
			e2=sqrt(e1**2+dk**2)
			e=(/e2*sign(1d0,-e1),-e2*sign(1d0,-e1)/)
			u2=(1d0-abs(e1)/e2)*0.5d0
			v2=(1d0+abs(e1)/e2)*0.5d0
			cos2=(1d0-abs(tmpa)/tmp)*0.5d0
			sin2=(1d0+abs(tmpa)/tmp)*0.5d0
			hn=(/u2(1)*cos2,u2(2)*sin2,v2(1)*cos2,v2(2)*sin2/)
			write(10,"(8e16.3)")(e(l),hn(l),l=1,4)
		enddo
		! do i=0,mp2
			! kx=pi/mp2*i
			! do j=0,mp2
				! tmp1=0
				! ky=pi/mp2*j
				! gk(1)=0.5d0*(cos(kx)-cos(ky))
				! ! gk(2)=0.5d0*(cos(3d0*kx)-cos(3d0*ky))
				! dk=dt(1)*gk(1) ! +dt(2)*gk(2)
				! tmpa=eka(kx,ky)
				! tmp=sqrt((tmpa)**2+US**2)
				! e1=eks(kx,ky)+(/-tmp,tmp/)*sign(1d0,tmpa)
				! e2=sqrt(e1**2+dk**2)
				! e=(/e2*sign(1d0,-e1),-e2*sign(1d0,-e1)/)
				! u2=(1d0-abs(e1)/e2)*0.5d0
				! v2=(1d0+abs(e1)/e2)*0.5d0
				! cos2=(1d0-abs(tmpa)/tmp)*0.5d0
				! sin2=(1d0+abs(tmpa)/tmp)*0.5d0
				! hn=(/u2(1)*cos2,u2(2)*sin2,v2(1)*cos2,v2(2)*sin2/)
				! do l=1,4
					! if(abs(e(l))<0.01) then
						! tmp1=tmp1+hn(l)
					! endif
				! enddo
				! write(20,"(3e16.3)")kx,ky,tmp1
			! enddo
			! write(20,"(1x)")
		! enddo
	enddo
	close(10)
	close(20)
	close(30)
	close(40)
end
			
