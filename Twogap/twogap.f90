program main
	implicit none
	real(8),parameter :: t=0.25d0,tp=-0.025d0,tpp=0.012d0,tppp=0.035d0,tiv=0d0,nf=1.12d0,&
		pi=3.1415926d0,cvg=1e-4
	real(8) :: te,bt,u=0.318,v=-0.3
	integer,parameter :: mp=128
	complex(8) :: img=(0d0,1d0)
	real(8) :: ek,f,e(4),kx,ky,work(12),n,sq,sqp,sp=-0.2,sa,sb,sp0,dk,dt(2),dtp(2),gk(2),h(4,4),wide,tmp,z
	integer :: i,j,k,l,p,info,a,b
	logical :: flaga,flagb
	ek(kx,ky)=-2d0*t*(cos(kx)+cos(ky))-4d0*tp*cos(kx)*cos(ky)-2d0*tpp*(cos(2d0*kx)+cos(2d0*ky))&
			-4d0*tppp*(cos(2d0*kx)*cos(ky)+cos(kx)*cos(2d0*ky))-4d0*tiv*cos(2d0*kx)*cos(2d0*ky)
	f(n)=1d0/(1d0+exp(bt*n))
	open(unit=10,file="../data/energy.dat")
	open(unit=20,file="../data/fermi.dat")
	open(unit=30,file="../data/search.dat")
	do p=60,60
		te=p
		bt=1d0/(te*8.6e-5)
		! u=0.375d0+p/20d0
		u=0.318d0
		v=-0.3d0
		wide=1d0
		sp=-0.2d0
		sp0=sp+wide
		dt=1d0
		sq=1d0
		dtp(1)=2d0
		! dtp(2)=-2d0
		sqp=2d0
		do while(any(abs(dtp-dt)>cvg).or.abs(sq-sqp)>cvg)
			! write(*,*)"*"
			dt=dtp
			sq=sqp
			sa=sp
			sb=sp
			n=0d0
			flaga=.true.
			flagb=.true.
			do while(abs(n-nf)>cvg)
				sp=0.5d0*(sa+sb)
				n=0d0
				sqp=0d0
				dtp=0d0
				!$OMP PARALLEL DO REDUCTION(+:n,sqp,dtp) PRIVATE(kx,ky,gk,dk,h,work,e,tmp) SCHEDULE(GUIDED)
				do i=-mp+1,mp
					do j=-mp+1,mp
						kx=pi/mp*i*sqrt(2d0)/2d0+pi/mp*j*sqrt(2d0)/2d0
						ky=-pi/mp*i*sqrt(2d0)/2d0+pi/mp*j*sqrt(2d0)/2d0
						gk(1)=0.5d0*(cos(kx)-cos(ky))
						! gk(2)=0.5d0*(cos(3d0*kx)-cos(3d0*ky))
						dk=dt(1)*gk(1) ! +dt(2)*gk(2)
						h=reshape((/ ek(kx,ky)-sp,      -0.5d0*u*sq,         dk,          0d0,           &
									-0.5d0*u*sq,      ek(kx+pi,ky+pi)-sp,     0d0,          -dk,        &
									 dk,             0d0,           -ek(kx,ky)+sp,    -0.5d0*u*sq,       &
									   0d0,            -dk,         -0.5d0*u*sq,  -ek(kx+pi,ky+pi)+sp /),(/4,4/))
						call dsyev( "v", "u", 4, h, 4, e, work, 12, info )
						! write(10,"(4e16.3)")(h(l,:),l=1,4)
						! stop
						! do k=-10*mp,10*mp
							! z=k/(9d0*mp)*e(4)
							! ! write(*,*)z
							! n=n-dot_product((h(1,:)*h(1,:)+h(2,:)*h(2,:)-h(3,:)*h(3,:)-h(4,:)*h(4,:)),dimag(1d0/(z-e(:)+img)))&
								! *f(z)/(2d0*pi**2)/((9d0*mp)*e(4))
						! enddo
						e=1d0/(1d0+exp(bt*e))
						n=n+dot_product((h(1,:)*h(1,:)+h(2,:)*h(2,:)-h(3,:)*h(3,:)-h(4,:)*h(4,:)),e)
						sqp=sqp+2d0*dot_product(h(1,:)*h(2,:)+h(3,:)*h(4,:),e)
						tmp=dot_product((h(1,:)*h(3,:)-h(2,:)*h(4,:)),e)
						dtp(1)=dtp(1)+v*gk(1)*tmp
						! dtp(2)=dtp(2)+v*gk(2)*tmp
						n=2d0+n
					enddo
				enddo
				!$OMP END PARALLEL DO
				n=n/(8d0*mp**2)
				sqp=sqp/(8d0*mp**2)
				dtp=dtp/(8d0*mp**2)
				! write(*,"(2e16.5)")sp,n
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
			wide=max(abs(sp0-sp),100*cvg)
			sp0=sp
		enddo
		write(*,"(2e16.3)")sqp,dtp(1)
		write(30,"(5e16.3)")te,u,v,sqp,dtp(1)
	enddo
	i=0
	j=0
	a=1
	b=0
	! sq=0.5d0
	! dt(1)=0.1d0
	! dt(2)=-0.05d0	
	do while(i>=0)
		if(i==mp.and.j==0) then
			a=0
			b=1
		endif
		if(i==mp.and.j==mp) then
			a=-1
			b=-1
		endif
		i=i+a
		j=j+b
		kx=pi/mp*i
		ky=pi/mp*j
		gk(1)=0.5d0*(cos(kx)-cos(ky))
		! gk(2)=0.5d0*(cos(3d0*kx)-cos(3d0*ky))
		dk=dt(1)*gk(1) ! +dt(2)*gk(2)
		h=reshape((/ ek(kx,ky)-sp,      -0.5d0*u*sq,         dk,          0d0,           &
					-0.5d0*u*sq,      ek(kx+pi,ky+pi)-sp,     0d0,          -dk,        &
					 dk,             0d0,           -ek(kx,ky)+sp,    -0.5d0*u*sq,       &
					   0d0,            -dk,         -0.5d0*u*sq,  -ek(kx+pi,ky+pi)+sp /),(/4,4/))
		call dsyev( "v", "u", 4, h, 4, e, work, 12, info )
		write(10,"(8e16.3)")(e(l),h(1,l)**2,l=1,4)
	enddo
	do i=0,mp
		do j=0,mp
			tmp=0
			kx=pi/mp*i
			ky=pi/mp*j
			gk(1)=0.5d0*(cos(kx)-cos(ky))
			! gk(2)=0.5d0*(cos(3d0*kx)-cos(3d0*ky))
			dk=dt(1)*gk(1) ! +dt(2)*gk(2)
			h=reshape((/ ek(kx,ky)-sp,      -0.5d0*u*sq,         dk,          0d0,           &
						-0.5d0*u*sq,      ek(kx+pi,ky+pi)-sp,     0d0,          -dk,        &
						dk,             0d0,           -ek(kx,ky)+sp,    -0.5d0*u*sq,       &
						0d0,            -dk,         -0.5d0*u*sq,  -ek(kx+pi,ky+pi)+sp /),(/4,4/))
			call dsyev( "v", "u", 4, h, 4, e, work, 12, info )
			do l=1,4
				if(abs(e(l))<0.05) then
					tmp=tmp+h(1,l)**2
				endif
			enddo
			write(20,"(3e16.3)")kx,ky,tmp
		enddo
		write(20,"(1x)")
	enddo
	close(10)
	close(20)
	close(30)
end
			