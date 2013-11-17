module global
	implicit none
	save
	integer,parameter :: dn=24,cof=1000
	real(8),parameter :: dv=1d0,nf0=0.85d0,du=2.44d0,t=1d0,tp=-0.25d0,eff=1d0,pi=3.141592653589793d0
	integer :: pbc(0:dn+1)
end module
program main
	use global
	implicit none
	character(25) :: it,ft,fmat(3)
	real(8) :: time(2)
	logical :: flaga,flagb,flag1,flag2,flag3
	 !general nd(:,1):spin down,nd(:,0):spin up
	real(8) :: nd(dn*dn,0:1),ndp(dn*dn,0:1),dt(dn*dn,4),dt1(dn*dn,4),dtf(dn*dn),rdom(dn*dn),er(3)=1,sa,sb,sp,sp0,nf,cvg=1e-4,wide
	integer :: i,ii,j,k,l,n,rc,im(5)=(/300,-1,-1,-1,-1/),np(4)
	! cpe
	real(8) :: va(2*dn*dn*13),hhn(0:2,2*dn*dn),a,b,emax,emin,ac_ba,tn(cof),ld=0.01,tmp_sh,gm=0.1d0
	integer :: ja(2*dn*dn*13),ia(2*dn*dn+1),n_mod(3),info
	open(unit=10,file='../data/order_crs.dat')
	open(unit=20,file='../data/error.dat')
	open(unit=30,file='../data/rundata.dat')
	open(unit=40,file='../data/initial.dat')
	open(unit=50,file='../data/gap.dat')
	open(unit=60,file='../data/check.dat')
	write(fmat(1),*)"(",dn,"(e15.6))"
	write(fmat(2),*)"(",dn*dn,"(e15.6))"
	write(fmat(3),*)"(",2*dn*dn,"(e15.6))"
200	format(24(e15.6))
300	format(625(e15.6))
500	format(100000e10.3)
	call fdate(it)
	call cpu_time(time(1))
	! call cpu_time(time(1))
	rc=0
	pbc=(/dn,(i,i=1,dn),1/)
	dt=0d0
	dt1=0d0
	call init_random_seed()
	call random_number(rdom)
	! rdom=0.5d0	!uniform initial value
	do i=1,dn*dn
		dt1(i,:)=(/0.07d0,-0.07d0,0.07d0,-0.07d0/)+((/rdom(i)-0.5d0,-rdom(i)+0.5d0,rdom(i)-0.5d0,-rdom(i)+0.5d0/))/10d0
	enddo
	! ndp(:,0)=rdom*nf0
	! ndp(:,1)=nf0-ndp(:,0)
	ndp(:,0)=nf0/2+(rdom-0.5d0)/20d0
	ndp(:,1)=nf0-ndp(:,0)
	do i=1,dn*dn
		dtf(i)=(dt1(i,1)+dt1(i,3)-dt1(i,2)-dt1(i,4))/4.0
	enddo
	! write(40,fmat(1))((dtf(dn*i+j),j=1,dn),i=0,dn-1)
	! write(40,fmat(1))(((-1)**(i+j)*0.5*(ndp(dn*i+j,0)-ndp(dn*i+j,1)),j=1,dn),i=0,dn-1)
	flag1=.true.
	flag2=.true.
	flag3=.true.
	tmp_sh=sinh(ld)
	nd=0
	sp=0.14
	wide=0.2
	sp0=sp+wide
	! im(1)=(dn/2-1)*dn+dn/2
	! im(2)=im(1)+1
	! im(3)=im(2)+dn
	! im(4)=im(1)+dn
	! im(5)=im(1)-1
	! do i=1,dn*dn
		! if(any(im==i)) then
			! nd(i,0)=1
		! endif
	! enddo
	! write(20,fmat(1))((nd(dn*i+j,0),j=1,dn),i=0,dn-1)
ex:	do while(any(er(1:2)>cvg).and.(rc<200))
		rc=rc+1
		dt=dt1
		nd=ndp
		nf=0d0
		sa=sp
		sb=sp
		flaga=.true.
		flagb=.true.
		do while(abs(nf-nf0)>cvg)
			sp=0.5d0*(sa+sb)
			! write(*,*)sp,"***"
			n=1
			ia(1)=1
			l=2
			! write(*,*)"crs va ja ia start"
			do i=1,dn*dn
				do k=1,0,-1
					if(any(im==i)) then
						va(n:n)=(-1)**k*(-1*du*nd(i,k)+sp)+(-1)**minloc(abs(im-i))*eff
						ja(n)=2*i-k
						n=n+1
					else
						va(n)=(-1)**k*(-1*du*nd(i,k)+sp)
						ja(n)=2*i-k
						n=n+1
					endif
					call nnp(i,.true.,np)
					do j=1,4
						va(n)=(-1)**k*t
						ja(n)=2*np(j)-k
						n=n+1
						va(n)=dt(i,j)
						ja(n)=2*np(j)+k-1
						n=n+1
					enddo
					call nnp(i,.false.,np)
					do j=1,4
						va(n)=(-1)**k*tp
						ja(n)=2*np(j)-k
						n=n+1
					enddo
					ia(l)=n
					l=l+1
				enddo
			enddo
			! write(*,*)"crs va ja ia finish,lanmax start"
			call lanmax(va,ja,ia,emax,emin,2*dn*dn,info)
			a=(emax-emin)/2d0+gm
			b=(emax+emin)/2d0
			! write(*,*)"lanmax finish,a=",a,"b=",b,"rescale va start"
			do i=1,2*dn*dn
				do j=ia(i),ia(i+1)-1
					if(ja(j)==i) then
						va(j)=va(j)-b
						exit
					endif
				enddo
			enddo
			va=va/a
			! write(*,*)"rescale va finish, produce tn start"
			ndp=0
			ac_ba=acos(-b/a)
			do n=1,cof
				tn(n)=-2.0d0*sin(n*ac_ba)/(n*pi)*sinh(ld*(1.0d0-real(n,8)/cof))/tmp_sh
			enddo
			!$omp parallel do reduction(+:ndp) private(hhn,n_mod,l,ii) schedule(guided)
			do i=1,2*dn*dn
				hhn=0
				hhn(0,i)=1d0
				l=(mod(i,2)-1)*(-1)
				ii=(i+1)/2
				ndp(ii,l)=ndp(ii,l)+hhn(0,i)*(1.0d0-ac_ba/pi)			
				call crsmv(va,ja,ia,hhn(0,:),hhn(1,:),2*dn*dn)
				ndp(ii,l)=ndp(ii,l)+hhn(1,i)*tn(1)
				do n=2,cof
					n_mod(1)=mod(n-2,3)
					n_mod(2)=mod(n-1,3)
					n_mod(3)=mod(n,3)
					call crsmv(va,ja,ia,hhn(n_mod(2),:),hhn(n_mod(3),:),2*dn*dn)
					hhn(n_mod(3),:)=2d0*hhn(n_mod(3),:)-hhn(n_mod(1),:)
					ndp(ii,l)=ndp(ii,l)+hhn(n_mod(3),i)*tn(n)
				enddo
			enddo
			!$omp end parallel do
			ndp(:,1)=1d0-ndp(:,1)
			nf=0d0
			do i=1,dn*dn
				nf=nf+ndp(i,0)+ndp(i,1)
			enddo
			nf=1.0d0/(dn*dn)*nf
			write(*,"(4e15.8)")nf,sp,a,b
			if(abs(nf-nf0)<=cvg) then
				exit
			endif			
			if(nf>2.or.nf<0) then
				write(*,*)"the a maybe small"
				gm=gm+0.1
				flaga=.true.
				flagb=.true.
				cycle
			endif
			if(nf<nf0) then
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
		wide=max(abs(sp0-sp),10*cvg)
		sp0=sp
		dt1=0d0
		!$omp parallel do reduction(+:dt1) private(hhn,n_mod,np) schedule(guided)
		do i=1,dn*dn
			call nnp(i,.true.,np)
			hhn=0
			hhn(0,2*i)=1d0
			do j=1,4
				dt1(i,j)=dt1(i,j)+hhn(1,2*np(j)-1)*(1-ac_ba/pi)
			enddo
			call crsmv(va,ja,ia,hhn(0,:),hhn(1,:),2*dn*dn)
			do j=1,4
				dt1(i,j)=dt1(i,j)+hhn(1,2*np(j)-1)*tn(1)
			enddo
			do n=2,cof
				n_mod(1)=mod(n-2,3)
				n_mod(2)=mod(n-1,3)
				n_mod(3)=mod(n,3)
				call crsmv(va,ja,ia,hhn(n_mod(2),:),hhn(n_mod(3),:),2*dn*dn)
				hhn(n_mod(3),:)=2*hhn(n_mod(3),:)-hhn(n_mod(1),:)
				do j=1,4
					dt1(i,j)=dt1(i,j)+hhn(n_mod(3),2*np(j)-1)*tn(n)
				enddo
			enddo
		enddo
		! do i=1,dn*dn
			! call nnp(i,.true.,np)
			! hhn=0
			! hhn(0,2*i-1)=1d0
			! do j=1,4
				! dt1(i,j)=dt1(i,j)+hhn(1,2*np(j))*(1-ac_ba/pi)
			! enddo
			! call crsmv(va,ja,ia,hhn(0,:),hhn(1,:),2*dn*dn)
			! do j=1,4
				! dt1(i,j)=dt1(i,j)+hhn(1,2*np(j))*tn(1)
			! enddo
			! do n=2,cof
				! n_mod(1)=mod(n-2,3)
				! n_mod(2)=mod(n-1,3)
				! n_mod(3)=mod(n,3)
				! call crsmv(va,ja,ia,hhn(n_mod(2),:),hhn(n_mod(3),:),2*dn*dn)
				! hhn(n_mod(3),:)=2*hhn(n_mod(3),:)-hhn(n_mod(1),:)
				! do j=1,4
					! dt1(i,j)=dt1(i,j)+hhn(n_mod(3),2*np(j))*tn(n)
				! enddo
			! enddo
		! enddo
		!$omp end parallel do
		dt1=-1d0*dt1*dv
		! rewind(50)
		! write(50,fmat(1))((dtf(dn*i+j),j=1,dn),i=0,dn-1)
		write(*,*)"dt=",dt1(1,1)
		er(1)=sum(abs(dt1-dt))/real(dn*dn,8)
		! er(2)=sum(abs(ndp(:,0)+ndp(:,1)-nd(:,0)-nd(:,1)))/real(2*dn*dn,8)
		er(2)=sum(abs(ndp(:,1)-nd(:,1)))/real(dn*dn,8)
		write(*,"(2e15.8)")er(1:2)
		! if(flag1.and.all(er<cvg*100)) then
			! do i=1,dn*dn
				! dtf(i)=(dt1(i,1)+dt1(i,3)-dt1(i,2)-dt1(i,4))/4.0
			! enddo
			! write(10,fmat(1))((dtf(dn*i+j),j=1,dn),i=0,dn-1)
			! write(10,fmat(1))((nd(dn*i+j,0)+nd(dn*i+j,1),j=1,dn),i=0,dn-1)
			! write(10,fmat(1))(((-1)**(i+j)*0.5*(nd(dn*i+j,0)-nd(dn*i+j,1)),j=1,dn),i=0,dn-1)
			! write(10,*)"error=",cvg*100
			! flag1=.false.
		! endif		
		! if(flag2.and.all(er<cvg*50)) then
			! do i=1,dn*dn
				! dtf(i)=(dt1(i,1)+dt1(i,3)-dt1(i,2)-dt1(i,4))/4.0
			! enddo
			! write(10,fmat(1))((dtf(dn*i+j),j=1,dn),i=0,dn-1)
			! write(10,fmat(1))((nd(dn*i+j,0)+nd(dn*i+j,1),j=1,dn),i=0,dn-1)
			! write(10,fmat(1))(((-1)**(i+j)*0.5*(nd(dn*i+j,0)-nd(dn*i+j,1)),j=1,dn),i=0,dn-1)
			! write(10,*)"error=",cvg*50
			! flag2=.false.
		! endif
		! if(flag3.and.all(er<cvg*10)) then
			! do i=1,dn*dn
				! dtf(i)=(dt1(i,1)+dt1(i,3)-dt1(i,2)-dt1(i,4))/4.0
			! enddo
			! write(10,fmat(1))((dtf(dn*i+j),j=1,dn),i=0,dn-1)
			! write(10,fmat(1))((nd(dn*i+j,0)+nd(dn*i+j,1),j=1,dn),i=0,dn-1)
			! write(10,fmat(1))(((-1)**(i+j)*0.5*(nd(dn*i+j,0)-nd(dn*i+j,1)),j=1,dn),i=0,dn-1)
			! write(10,*)"error=",cvg*10
			! flag3=.false.
		! endif
	enddo ex
	dtf=(dt1(:,1)+dt1(:,3)-dt1(:,2)-dt1(:,4))/4.0
	call fdate(ft)
	call cpu_time(time(2))
	! call cpu_time(time(2))
	write(10,fmat(1))((dtf(dn*i+j),j=1,dn),i=0,dn-1)
	write(10,fmat(1))((nd(dn*i+j,0)+nd(dn*i+j,1),j=1,dn),i=0,dn-1)
	write(10,fmat(1))(((-1)**(i+j)*0.5*(nd(dn*i+j,0)-nd(dn*i+j,1)),j=1,dn),i=0,dn-1)
	write(10,*)"v=",dv,"u=",du,"t'=",tp,"heff=",eff,"nf=",nf0,"n=",dn
	write(10,*)"info=",info
	write(10,*)"the tatal runtime is form ",it," to ",ft,"cost cup time:",time(2)-time(1),"s"
	close(10)
	close(20)
	close(30)
	close(40)
	close(50)
end program main

subroutine nnp(ii,near,jj)
	use global
	implicit none
	integer :: ii,i,j,jj(4)
	logical :: near
	j=pbc(mod(ii,dn))
	i=(ii-j)/dn+1
	if(near) then
		jj(1)=(i-1)*dn+pbc(j+1)
		jj(2)=(pbc(i+1)-1)*dn+j
		jj(3)=(i-1)*dn+pbc(j-1)
		jj(4)=(pbc(i-1)-1)*dn+j
	else
		jj(1)=dn*(pbc(i+1)-1)+pbc(j+1)
		jj(2)=dn*(pbc(i+1)-1)+pbc(j-1)
		jj(3)=dn*(pbc(i-1)-1)+pbc(j+1)
		jj(4)=dn*(pbc(i-1)-1)+pbc(j-1)
	endif
end subroutine nnp
subroutine crsmv(va,ja,ia,x,y,m)
	implicit none
	integer :: i,j,m
	integer :: ja(*),ia(*)
	real(8) :: va(*),y(*),x(*)
	do i=1,m
		y(i)=0
		do j=ia(i),ia(i+1)-1
			y(i)=y(i)+va(j)*x(ja(j))
		enddo
	enddo
end
subroutine init_random_seed()
	integer :: i, n, clock
	integer, dimension(:), allocatable :: seed
	call random_seed(size = n)
	allocate(seed(n))
	call system_clock(count=clock)
	seed = clock + 37 * (/ (i - 1, i = 1, n) /)
	call random_seed(put = seed)
	deallocate(seed)
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                                                                                                                  !!!!!
!!!!!                       use lanczos to detect maximum-minimum eigenvalue store in crs formal                       !!!!!
!!!!!                                                                                                                  !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine lanmax(va,ja,ia,emax,emin,m,info)
	implicit none
	integer :: ja(*),ia(*),m,modi(0:2),i,j,info
	real(8) :: va(*),emax,emin,sun,work
	real(8),allocatable :: v(:,:),ap(:),bt(:),tmp(:),e(:),fe(:)
	allocate(v(m,0:2),ap(m),bt(0:m),tmp(0:m),e(m),fe(m))
	fe(1)=1000
	v(:,0)=0d0
	bt(0)=0d0
	call random_number(v(:,1))
	sun=dot_product(v(:,1),v(:,1))
	v(:,1)=v(:,1)/sqrt(abs(sun))
	do i=1,m
		modi(2)=mod(i+1,3)
		modi(1)=mod(i,3)
		modi(0)=mod(i-1,3)
		call crsmv(va,ja,ia,v(:,modi(1)),v(:,modi(2)),m)
		ap(i)=dot_product(v(:,modi(1)),v(:,modi(2)))
		e(1:i)=ap(1:i)
		tmp(1:i-1)=bt(1:i-1)
		call dsteqr('n',i,e(1:i),tmp(1:i-1),work,i,work,info)
		! write(10,*)e(1:i)
		if(abs(e(1)-fe(1))<1e-3.and.abs(e(i)-fe(i-1))<1e-3.or.info/=0) then
			write(*,*)"info=",info
			exit
		endif
		fe=e
		v(:,modi(2))=v(:,modi(2))-ap(i)*v(:,modi(1))-bt(i-1)*v(:,modi(0))
		sun=dot_product(v(:,modi(2)),v(:,modi(2)))
		bt(i)=sqrt(abs(sun))
		v(:,modi(2))=v(:,modi(2))/bt(i)
	enddo
	emax=e(min(i,m))
	emin=e(1)
	deallocate(v,ap,bt,tmp,e,fe)
end
