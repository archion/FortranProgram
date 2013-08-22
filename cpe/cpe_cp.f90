module global
	implicit none
	save
	integer,parameter :: dn=24
	real(8),parameter :: bt=1e5,dv=1,nf0=0.85,du=2.44d0,t=1,tp=-0.25d0,eff=0.3,pi=3.14159265d0
	integer :: pbc(0:dn+1)
end module
program main
	use global
	implicit none
	character(25) :: it,ft,fmat(3)
	logical :: flaga,flagb,flag1,flag2,flag3
	 !general nd(:,1):spin down,nd(:,0):spin up
	real(8) :: nd(dn*dn,0:1),ndp(dn*dn,0:1),dt(dn*dn,4),dt1(dn*dn,4),dtf(dn*dn),rdom(dn*dn),er(3)=1,sa,sb,sp,sp0,nf,cvg=1e-3,wide
	integer :: i,j,k,l,n,rc,im(1)=-1,np(4),shf(dn*dn)
	! cpe
	real(8) :: va(2*dn*dn*13),hhn(0:2,2*dn*dn),a,b,emax,emin,ac_ba,tn
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
	do i=1,dn*dn
		shf(i)=i
	enddo
	call fisher_yates_shuffle(shf,size(shf))
	! im=shf(1)
	do i=1,dn*dn
		dt1(i,:)=(/0.02d0,-0.02d0,0.02d0,-0.02d0/)
	enddo
	nd=nf0/2
	sp=0.142
ex:	do while(abs(er(1))>cvg)
		rc=rc+1
		dt=dt1
		n=1
		ia(1)=1
		l=2
		do i=1,dn*dn
			do k=1,0,-1
				if(any(im==i)) then
					va(n:n)=(-1)**k*(-1*du*nd(i,k)+sp+eff)
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
		a=(emax-emin)/2+0.3
		b=(emax+emin)/2
		write(*,*)emax,emin,info
		do i=1,2*dn*dn
			do j=ia(i),ia(i+1)-1
				if(ja(j)==i) then
					va(j)=va(j)-b
					exit
				endif
			enddo
		enddo
		va=va/a
		dt1=0
		ac_ba=acos(-b/a)
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
				tn=-2.0*a*sin(n*ac_ba)/(n*pi)
				n_mod(1)=mod(n-2,3)
				n_mod(2)=mod(n-1,3)
				n_mod(3)=mod(n,3)
				call crsmv(va,ja,ia,hhn(n_mod(2),:),hhn(n_mod(3),:),2*dn*dn)
				hhn(n_mod(3),:)=2*hhn(n_mod(3),:)-hhn(n_mod(1),:)
				do j=1,4
					dt1(i,j)=dt1(i,j)+hhn(n_mod(3),2*np(j)-1)*tn
				enddo
			enddo
		enddo
		dt1(i,j)=dv*dt1(i,j)
		er(1)=sum(abs(dt1-dt))/(dn*dn)
		rewind(50)
		do i=1,dn*dn
			dtf(i)=(dt1(i,1)+dt1(i,3)-dt1(i,2)-dt1(i,4))/4.0
		enddo
		write(50,fmat(1))((dtf(dn*i+j),j=1,dn),i=0,dn-1)
	enddo ex
	do i=1,dn*dn
		dtf(i)=(dt(i,1)+dt(i,3)-dt(i,2)-dt(i,4))/4.0
	enddo
	call fdate(ft)
	write(10,fmat(1))((dtf(dn*i+j),j=1,dn),i=0,dn-1)
	write(10,fmat(1))((nd(dn*i+j,0)+nd(dn*i+j,1),j=1,dn),i=0,dn-1)
	write(10,fmat(1))(((-1)**(i+j)*0.5*(nd(dn*i+j,0)-nd(dn*i+j,1)),j=1,dn),i=0,dn-1)
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

subroutine fisher_yates_shuffle(a,m)
	implicit none
	integer :: a(*),i,j,tem,m
	real(8) :: ram
	do i=m,1,-1
		call random_number(ram)
		j=ceiling(ram*i)
		tem=a(i)
		a(i)=a(j)
		a(j)=tem
	enddo
end
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
