module global
	implicit none
	save
	integer,parameter :: dn=24
	real(8),parameter :: bt=1e5,dv=1,nf0=0.85,du=2.44,t=1,tp=-0.25,eff=1,pi=3.14159265d0
	integer :: pbc(0:dn+1)
end module
program main
	use global
	implicit none
	character(25) :: it,ft,fmat(3)
	logical :: flaga,flagb,flag1,flag2,flag3
	 !general nd(:,1):spin down,nd(:,0):spin up
	real(8) :: nd(dn*dn,0:1),ndp(dn*dn,0:1),dt(dn*dn,4),dt1(dn*dn,4),dtf(dn*dn),rdom(dn*dn),er(3)=1,sa,sb,sp,sp0,nf,cvg=1e-3,wide
	integer :: i,j,k,l,n,rc,im(5)=-1,np(4)
	! cpe
	real(8) :: va(2*dn*dn*13),hhn(0:2,2*dn*dn),a,b,emax,emin,ac_ba,tn,an(1000)
	integer :: ja(2*dn*dn*13),ia(2*dn*dn+1),cutoff=1000,n_mod(3),info
	! lapack
	! real(8) :: e(2*dn*dn),h(2*dn*dn,2*dn*dn),work(3*2*dn*dn)
	! real(8) :: f,en
	! integer :: info=0,error=0
	! f(en)=1/(1+exp(bt*en))
	open(unit=10,file='../data/output.dat')
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
	rc=0
	pbc=(/dn,(i,i=1,dn),1/)
	dt=0d0
	dt1=0d0
	call init_random_seed()
	call random_number(rdom)
	! rdom=0.5d0	!uniform initial value
	do i=1,dn*dn
		dt1(i,:)=(/-0.2d0,0.2d0,-0.2d0,0.2d0/) !  +((/rdom(i)-0.5d0,-rdom(i)+0.5d0,rdom(i)-0.5d0,-rdom(i)+0.5d0/))/10d0
	enddo
	! ndp(:,0)=rdom*nf0
	! ndp(:,1)=nf0-ndp(:,0)
	ndp(:,0)=nf0/2+(rdom-0.5d0)/20d0
	ndp(:,1)=nf0-ndp(:,0)
	do i=1,dn*dn
		dtf(i)=(dt1(i,1)+dt1(i,3)-dt1(i,2)-dt1(i,4))/4.0
	enddo
	write(40,fmat(1))((dtf(dn*i+j),j=1,dn),i=0,dn-1)
	write(40,fmat(1))(((-1)**(i+j)*0.5*(ndp(dn*i+j,0)-ndp(dn*i+j,1)),j=1,dn),i=0,dn-1)
	flag1=.true.
	flag2=.true.
	flag3=.true.
	nd=0
	sp=0
	wide=0.1
	sp0=sp+wide
	! im(1)=121
	! im(2)=123
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
ex:	do while(any(er>cvg))
		rc=rc+1
		dt=dt1
		nd=ndp
		nf=0
		sa=sp
		sb=sp
		flaga=.true.
		flagb=.true.
		write(*,*)"next"
		do while(abs(nf-nf0)>cvg)
			sp=0.5*(sa+sb)
			write(*,*)sp,"***"
			n=1
			ia(1)=1
			l=2
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
			call lanmax(va,ja,ia,emax,emin,2*dn*dn,info)
			a=(emax-emin)/2+0.1
			b=(emax+emin)/2
			! write(*,*)emax,emin,info
			do i=1,2*dn*dn
				do j=ia(i),ia(i+1)-1
					if(ja(j)==i) then
						va(j)=va(j)-b
						exit
					endif
				enddo
			enddo
			va=va/a
			nf=0
			ndp=0
			ac_ba=acos(-b/a)
			!$omp parallel do reduction(+:nf) private(hhn) num_threads(42)
			do i=1,2*dn*dn
				hhn=0
				hhn(0,i)=1d0
				l=(mod(i,2)-1)*(-1)
				ndp(i,l)=ndp(i,l)+a*hhn(0,i)*(1-ac_ba/pi)
				call crsmv(va,ja,ia,hhn(0,:),hhn(1,:),2*dn*dn)
				ndp(i,l)=ndp(i,l)-a*hhn(1,i)*2.0*sin(ac_ba)/pi
				do n=2,cutoff
					tn=-2.0*sin(n*ac_ba)/(n*pi)
					n_mod(1)=mod(n-2,3)
					n_mod(2)=mod(n-1,3)
					n_mod(3)=mod(n,3)
					call crsmv(va,ja,ia,hhn(n_mod(2),:),hhn(n_mod(3),:),2*dn*dn)
					hhn(n_mod(3),:)=2*hhn(n_mod(3),:)-hhn(n_mod(1),:)
					ndp(i,l)=ndp(i,l)+a*hhn(n_mod(3),i)*tn
					! nf=nf+(-1)**mod(i,2)*hhn(n_mod(3),i)*tn
					! write(20,*)hhn(n_mod(3),i),i
				enddo
				ndp(i,1)=1-ndp(i,1)
				nf=nf+ndp(i,0)+ndp(i,1)
				! if(ndp(i,0)<0.or.ndp(i,1)<0) then
					! write(*,*)ndp(i,0),ndp(i,1),"below the zero exit!!"
				! endif
			enddo
			!$omp end parallel do
			nf=1.0/(dn*dn)*nf
			write(*,*)nf
			if(abs(nf-nf0)<=cvg) then
				exit
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
		wide=max(abs(sp0-sp),0.01)
		sp0=sp
		ndp=0
		dt1=0
		!$omp parallel do private(hhn) num_threads(42)
		do i=1,dn*dn
			call nnp(i,.true.,np)
			hhn=0
			hhn(0,2*i)=1d0
			do j=1,4
				dt1(i,j)=dt1(i,j)+a*hhn(1,2*np(j)-1)*(1-ac_ba/pi)
			enddo
			call crsmv(va,ja,ia,hhn(0,:),hhn(1,:),2*dn*dn)
			do j=1,4
				dt1(i,j)=dt1(i,j)-a*hhn(1,2*np(j)-1)*2.0*sin(ac_ba)/pi
			enddo
			do n=2,cutoff
				tn=-2.0*sin(n*ac_ba)/(n*pi)
				n_mod(1)=mod(n-2,3)
				n_mod(2)=mod(n-1,3)
				n_mod(3)=mod(n,3)
				call crsmv(va,ja,ia,hhn(n_mod(2),:),hhn(n_mod(3),:),2*dn*dn)
				hhn(n_mod(3),:)=2*hhn(n_mod(3),:)-hhn(n_mod(1),:)
				do j=1,4
					dt1(i,j)=dt1(i,j)+a*hhn(n_mod(3),2*np(j)-1)*tn
				enddo
			enddo
		enddo
		!$omp end parallel do
		! do i=1,dn*dn
			! do n=1,2*dn*dn
				! nd(i,0)=nd(i,0)+h(2*i-1,n)**2*f(e(n))
				! nd(i,1)=nd(i,1)+h(2*i,n)**2*(1-f(e(n)))
			! enddo
		! enddo
		rewind(30)
		write(30,fmat(1))((ndp(dn*i+j,0)+ndp(dn*i+j,1),j=1,dn),i=0,dn-1)
		! dt1=0
		! do i=1,dn*dn
			! call nnp(i,.true.,np)
			! do j=1,4
				! do n=1,2*dn*dn
					! dt1(i,j)=dt1(i,j)+0.25*dv*(h(2*i-1,n)*h(2*np(j),n)+h(2*np(j)-1,n)*h(2*i,n))*tanh(0.5*bt*e(n))
				! enddo
			! enddo
			! write(30,"(4e10.3,5i4)")dt1(i,:),i,np
		! enddo
		er(1)=sum(abs(dt1-dt))/(dn*dn)
		er(2)=sum(abs(ndp(:,0)-nd(:,0)))/(dn*dn)
		er(3)=sum(abs(ndp(:,1)-nd(:,1)))/(dn*dn)
		write(*,"(3e15.6,f15.12)")er,sa
		if(mod(rc,20)==0.and.rc>=100) then
			do i=1,dn*dn
				dtf(i)=(dt1(i,1)+dt1(i,3)-dt1(i,2)-dt1(i,4))/4.0
			enddo
			write(50,*)"n=",rc
			write(50,fmat(1))((dtf(dn*i+j),j=1,dn),i=0,dn-1)
		endif
		if(flag1.and.all(er<cvg*100)) then
			do i=1,dn*dn
				dtf(i)=(dt1(i,1)+dt1(i,3)-dt1(i,2)-dt1(i,4))/4.0
			enddo
			write(10,fmat(1))((dtf(dn*i+j),j=1,dn),i=0,dn-1)
			write(10,fmat(1))((ndp(dn*i+j,0)+ndp(dn*i+j,1),j=1,dn),i=0,dn-1)
			write(10,fmat(1))(((-1)**(i+j)*0.5*(ndp(dn*i+j,0)-ndp(dn*i+j,1)),j=1,dn),i=0,dn-1)
			write(10,*)"error=",cvg*100
			flag1=.false.
		endif		
		if(flag2.and.all(er<cvg*50)) then
			do i=1,dn*dn
				dtf(i)=(dt1(i,1)+dt1(i,3)-dt1(i,2)-dt1(i,4))/4.0
			enddo
			write(10,fmat(1))((dtf(dn*i+j),j=1,dn),i=0,dn-1)
			write(10,fmat(1))((ndp(dn*i+j,0)+ndp(dn*i+j,1),j=1,dn),i=0,dn-1)
			write(10,fmat(1))(((-1)**(i+j)*0.5*(ndp(dn*i+j,0)-ndp(dn*i+j,1)),j=1,dn),i=0,dn-1)
			write(10,*)"error=",cvg*50
			flag2=.false.
		endif
		if(flag3.and.all(er<cvg*10)) then
			do i=1,dn*dn
				dtf(i)=(dt1(i,1)+dt1(i,3)-dt1(i,2)-dt1(i,4))/4.0
			enddo
			write(10,fmat(1))((dtf(dn*i+j),j=1,dn),i=0,dn-1)
			write(10,fmat(1))((ndp(dn*i+j,0)+ndp(dn*i+j,1),j=1,dn),i=0,dn-1)
			write(10,fmat(1))(((-1)**(i+j)*0.5*(ndp(dn*i+j,0)-ndp(dn*i+j,1)),j=1,dn),i=0,dn-1)
			write(10,*)"error=",cvg*10
			flag3=.false.
		endif
	enddo ex
	do i=1,dn*dn
		dtf(i)=(dt(i,1)+dt(i,3)-dt(i,2)-dt(i,4))/4.0
	enddo
	call fdate(ft)
	write(10,fmat(1))((dtf(dn*i+j),j=1,dn),i=0,dn-1)
	write(10,fmat(1))((nd(dn*i+j,0)+nd(dn*i+j,1),j=1,dn),i=0,dn-1)
	write(10,fmat(1))(((-1)**(i+j)*0.5*(nd(dn*i+j,0)-nd(dn*i+j,1)),j=1,dn),i=0,dn-1)
	write(10,*)"the beta is",bt,"v=",dv,"u=",du,"t'=",tp,"heff=",eff,"nf=",nf0,"n=",dn
	write(10,*)"info=",info
	write(10,*)"the tatal runtime is form ",it," to ",ft
	! if(info/=0) then
		! error=mod(info,2*dn*dn+1)
		! do j=1,error+1
			! write(20,500)(h(i,j),i=1,error+1)
		! enddo
	! endif
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
		if(abs(e(1)-fe(1))<1e-3.or.info/=0) then
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
