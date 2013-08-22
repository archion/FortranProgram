module global
	implicit none
	save
	real(8) :: t=1
	integer, parameter :: dn=4
	real(8), parameter :: pi=3.14159
	integer :: pbc(0:dn+1)
end module
program cpe
	use global
	implicit none
	integer,parameter :: sz=1600,cof=1024
	integer :: i,j,ja(dn*dn*5),ia(dn*dn+1),nj(4),n,n_mod(3),shf(sz),info
	real(8) :: va(dn*dn*5),x,w,hhn(0:2,dn*dn),p(0:2),emax,emin,a,b,an(0:2*cof-1),dos(0:2*cof-1),ld=0.01,tmp_sh
	open(unit=10,file="../data/output.dat")
	pbc=(/dn,(i,i=1,dn),1/)
	call init_random_seed()
	do i=1,sz
		shf(i)=-i+100
	enddo
	call fisher_yates_shuffle(shf,size(shf))
	n=1
	ia(1)=1
	do i=1,dn*dn
		! va(n)=-4
		! ja(n)=i
		! n=n+1
		call nnp(i,.true.,nj)
		do j=1,4
			va(n)=t
			ja(n)=nj(j)
			n=n+1
		enddo
		ia(i+1)=n
	enddo
	call lanmax(va,ja,ia,emax,emin,size(ia)-1,info)
	write(*,*)emax,emin,info
	a=(emax-emin)/2+0.1
	b=(emax+emin)/2
	write(*,*)a,b
	do i=1,dn*dn
		do j=ia(i),ia(i+1)-1
			if(ja(j)==i) then
				va(j)=va(j)-b
			endif
		enddo
	enddo
	va=va/a
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                                                                                                                  !!!!!
!!!!!                                     expanded by tchebycheff polynomial                                           !!!!!
!!!!!                                                                                                                  !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! do j=1,sz
		! x=shf(j)/(sz-99d0)
		! if(x*a+b<emax.and.x*a+b>emin) then
			! dos=0
			! w=1.0/sqrt(1-x**2)
			! !$omp parallel do reduction(+:dos) private(hhn) num_threads(42)
			! do i=1,dn*dn
				! hhn=0
				! hhn(0,i)=1d0
				! dos=dos+hhn(0,i)/pi*w
				! call crsmv(va,ja,ia,hhn(0,:),hhn(1,:),dn*dn)
				! dos=dos+hhn(1,i)/(2.0*pi)*w*cos(acos(x))
				! do n=2,cof
					! n_mod(1)=mod(n-2,3)
					! n_mod(2)=mod(n-1,3)
					! n_mod(3)=mod(n,3)
					! call crsmv(va,ja,ia,hhn(n_mod(2),:),hhn(n_mod(3),:),dn*dn)
					! hhn(n_mod(3),:)=2*hhn(n_mod(3),:)-hhn(n_mod(1),:)
					! dos=dos+hhn(n_mod(3),i)/(2.0*pi)*w*cos(n*acos(x))
				! enddo
			! enddo
			! !$omp end parallel do
			! write(10,*)x*a+b,dos
		! endif
	! enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                                                                                                                  !!!!!
!!!!!                                 expanded by tchebycheff polynomial using fft                                     !!!!!
!!!!!                                                                                                                  !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	dos=0
	tmp_sh=sinh(ld)
	!$omp parallel do reduction(+:dos) private(hhn,an,n_mod) num_threads(6) schedule(guided)
	do i=1,dn*dn
		hhn=0
		hhn(0,i)=1d0
		an(0)=hhn(0,i)
		call crsmv(va,ja,ia,hhn(0,:),hhn(1,:),dn*dn)
		an(1)=hhn(1,i) !*sinh(ld*(1.0-1.0/cof))/tmp_sh
		do n=2,cof-1
			n_mod(1)=mod(n-2,3)
			n_mod(2)=mod(n-1,3)
			n_mod(3)=mod(n,3)
			call crsmv(va,ja,ia,hhn(n_mod(2),:),hhn(n_mod(3),:),dn*dn)
			hhn(n_mod(3),:)=2*hhn(n_mod(3),:)-hhn(n_mod(1),:)
			an(n)=hhn(n_mod(3),i) !*sinh(ld*(1.0-real(n)/cof))/tmp_sh
		enddo
		an(cof:2*cof-1)=0
		call fftcos1_2(an,2*cof,-1)
		dos=dos+an
	enddo
	!$omp end parallel do
	write(10,"(2e16.3)")(cos(pi*(i+0.5)/(2*cof-1)),dos(i),i=0,2*cof-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                                                                                                                  !!!!!
!!!!!                                  expanded by legendre ploynomial                                                 !!!!!
!!!!!                                                                                                                  !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! do j=1,sz
		! x=shf(j)/(sz-99d0)
		! if(x*(bw+0.5)<bw.and.x*(bw+0.5)>-1*bw) then
			! dos=0
			! !$omp parallel do reduction(+:dos) private(hhn,p) num_threads(42)
			! do i=1,dn*dn
				! hhn=0
				! hhn(0,i)=1d0
				! p(0)=1
				! p(1)=x
				! dos=dos+hhn(0,i)*0.5*p(0)
				! call crsmv(va,ja,ia,hhn(0,:),hhn(1,:),dn*dn)
				! dos=dos+1.5*hhn(1,i)*p(1)
				! do n=2,cof
					! n_mod(1)=mod(n-2,3)
					! n_mod(2)=mod(n-1,3)
					! n_mod(3)=mod(n,3)
					! call crsmv(va,ja,ia,hhn(n_mod(2),:),hhn(n_mod(3),:),dn*dn)
					! hhn(n_mod(3),:)=(2*n+1)*1.0/(n+1)*hhn(n_mod(3),:)-n*1.0/(n+1)*hhn(n_mod(1),:)
					! p(n_mod(3))=(2*n+1)*1.0/(n+1)*x*p(n_mod(2))-n*1.0/(n+1)*p(n_mod(1))
					! dos=dos+(2*n+1)/2.0*hhn(n_mod(3),i)*p(n_mod(3))
				! enddo
			! enddo
			! !$omp end parallel do
			! write(10,*)x*(bw+0.5),dos
		! endif
	! enddo
	close(10)
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                                                                                                                  !!!!!
!!!!!                                  find nearest and next-nearest site                                              !!!!!
!!!!!                                                                                                                  !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                                                                                                                  !!!!!
!!!!!                                  vector-matrix product use crs storage                                           !!!!!
!!!!!                                                                                                                  !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

subroutine init_random_seed()
	integer :: i, n, clock
	integer, dimension(:), allocatable :: seed
	call random_seed(size = n)
	allocate(seed(n))
	call system_clock(count=clock)
	seed = clock + 37 * (/ (i - 1, i = 1, n) /)
	call random_seed(put = seed)
	deallocate(seed)
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!       fft for complex array in one-demension              !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft1(a,n,sig)
	implicit none
	real(8),parameter :: pi=3.1415926
	integer :: i,j,m,sig,n
	complex(8) :: a(0:n-1),tmp,wp,w
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!          bit-reversal using gold-rader algorithm          !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	j=n/2
	do i=1,n-2
		if(j>i) then
			tmp=a(i)
			a(i)=a(j)
			a(j)=tmp
		endif
		m=n/2
		do while(j>=m)
			j=j-m
			m=m/2
		enddo
		j=j+m
	enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                     fast fourier transform using danielson-lanczos lemma                          !!!!
!!!!     trigonometric recurrence relations using 5.5.6 in book "numerical recipes in fortran 77"      !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	m=2
	do while(m<=n)
		wp=cmplx(-2d0*sin(sig*pi/m)**2,sin(sig*2d0*pi/m))
		w=cmplx(1d0,0d0)
		do j=0,m/2-1
			do i=j,n-1,m
				tmp=w*a(i+m/2)
				a(i+m/2)=a(i)-tmp
				a(i)=a(i)+tmp
			enddo
			w=w*wp+w
		enddo
		m=m*2
	enddo
	if(sig==-1) then
		a=a/n
	endif
end




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!          fft for real array in one-demension              !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fftreal1(a,b,n)
	implicit none
	real(8),parameter :: pi=3.1415926
	integer :: i,m,n
	real(8) :: a(0:n-1)
	complex(8) :: b(0:n-1),h1(0:n/2-1),h2(0:n/2-1),wp,w
	b=cmplx(a(0:n-2:2),a(1:n-1:2))
	m=n/2
	call fft1(b(0:m-1),m,1)
	h1(0)=cmplx(real(b(0)),0)
	h2(0)=cmplx(aimag(b(0)),0)
	h1(1:m/2)=0.5*(b(1:m/2)+conjg(b(m-1:m/2:-1)))
	h2(1:m/2)=cmplx(0,-0.5)*(b(1:m/2)-conjg(b(m-1:m/2:-1)))
	h1(m-1:m/2:-1)=conjg(h1(1:m/2))
	h2(m-1:m/2:-1)=conjg(h2(1:m/2))
	wp=cmplx(-2d0*sin(pi/n)**2,sin(2d0*pi/n))
	w=cmplx(1d0,0d0)
	do i=0,m-1
		b(i)=h1(i)+w*h2(i)
		b(i+m)=h1(i)-w*h2(i)
		w=w*wp+w
	enddo
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                      fft for real array in one-demension                            !!!!
!!!!               using input real array to storage the output complex array            !!!!
!!!!     with even index refer to real parts and odd index refer to imaginary parts      !!!!
!!!!           and a(0) is 0th while a(1) is n/2th of the fft of input data              !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fftreal1_rp(a,n,sig)
	implicit none
	real(8),parameter :: pi=3.1415926d0
	integer :: i,m,n,sig
	real(8) :: a(0:n-1)
	complex(8) :: b(0:n/2-1),wp,w,he,ho,tmp
	m=n/2
	if(sig==1) then
		b=cmplx(a(0:n-2:2),a(1:n-1:2))
		call fft1(b,m,1)
		wp=cmplx(-2d0*sin(pi/n)**2,sin(2d0*pi/n))
		w=cmplx(1d0,0d0)
		w=w*wp+w
		a(0)=real(b(0))+aimag(b(0))
		a(1)=real(b(0))-aimag(b(0))
		do i=1,m-1
			he=0.5*(b(i)+conjg(b(m-i)))
			ho=cmplx(0,-0.5)*(b(i)-conjg(b(m-i)))
			tmp=he+w*ho
			a(2*i)=real(tmp)
			a(2*i+1)=aimag(tmp)
			w=w*wp+w
		enddo
	else
		b(1:m-1)=cmplx(a(2:n-2:2),a(3:n-1:2))
		wp=cmplx(-2d0*sin(-1d0*pi/n)**2,sin(-2d0*pi/n))
		w=cmplx(1d0,0d0)
		w=w*wp+w
		do i=1,m/2
			he=0.5*(b(i)+conjg(b(m-i)))
			ho=cmplx(0d0,0.5d0)*(b(i)-conjg(b(m-i)))
			b(i)=he+w*ho
			b(m-i)=conjg(he)-conjg(w*ho)
			w=w*wp+w
		enddo
		b(0)=0.5*cmplx(a(0)+a(1),a(0)-a(1))
		call fft1(b,m,-1)
		do i=0,m-1
			a(2*i)=real(b(i))
			a(2*i+1)=aimag(b(i))
		enddo
	endif
end

subroutine fftsin1(a,n)
	implicit none
	real(8),parameter :: pi=3.1415926d0
	integer :: i,n
	real(8) :: a(0:n-1),tmp
	complex(8) :: wp,w
	wp=cmplx(-2d0*sin(0.5d0*pi/n)**2,sin(pi/n))
	w=cmplx(1,0)
	a(0)=0
	do i=1,n/2
		w=w*wp+w
		tmp=0.5*(a(i)-a(n-i))
		a(i)=aimag(w)*(a(i)+a(n-i))+tmp
		a(n-i)=a(i)-2*tmp
	enddo
	call fftreal1_rp(a,n,1)
	a(1)=0.5*a(0)
	a(0)=0
	do i=3,n-1,2
		tmp=a(i)
		a(i)=a(i-2)+a(i-1)
		a(i-1)=tmp
	enddo
end
subroutine fftcos1(a,n)
	implicit none
	real(8),parameter :: pi=3.1415926d0
	integer :: i,n
	real(8) :: a(0:n),sun,tmp
	complex(8) :: wp,w
	wp=cmplx(-2d0*sin(0.5d0*pi/n)**2,sin(pi/n))
	w=cmplx(1,0)
	sun=0.5*(a(0)-a(n))
	a(0)=0.5*(a(0)+a(n))
	do i=1,n/2
		w=w*wp+w
		tmp=a(n-i)-a(i)
		a(i)=0.5*(a(i)+a(n-i))+aimag(w)*tmp
		a(n-i)=a(i)-2*aimag(w)*tmp
		sun=sun-tmp*real(w)
	enddo
	call fftreal1_rp(a(0:n-1),n,1)
	a(n)=a(1)
	a(1)=sun
	do i=3,n-1,2
		a(i)=a(i-2)+a(i)
	enddo
end
subroutine fftcos1_2(a,n,sig)
	implicit none
	real(8),parameter :: pi=3.1415926d0
	integer :: i,n,m,sig
	real(8) :: a(0:n-1),tmp,tmp0,tmp1
	complex(8) :: wp,w,wi,wpi
	m=n/2
	w=cmplx(cos(0.5d0*pi/n),sin(0.5d0*pi/n))
	wi=cmplx(cos(0.5d0*pi),sin(0.5d0*pi))
	wp=cmplx(-2d0*sin(0.5d0*pi/n)**2,sin(pi/n))
	wpi=cmplx(-2d0*sin(-0.5d0*pi/n)**2,sin(-1d0*pi/n))
	if(sig==1) then
		do i=0,m-1
			tmp=aimag(w)*(a(n-i-1)-a(i))
			a(i)=0.5*(a(i)+a(n-i-1))+tmp
			a(n-i-1)=a(i)-2*tmp
			w=w*wp+w
		enddo
		call fftreal1_rp(a,n,sig)
		tmp=a(n-1)
		a(n-1)=0.5*a(1)
		do i=n-1,3,-2
			wi=wi*wpi+wi
			tmp0=a(i-1)
			tmp1=a(i-2)
			a(i-1)=real(wi)*a(i-1)-aimag(wi)*tmp
			a(i-2)=aimag(wi)*tmp0+real(wi)*tmp+a(i)
			tmp=tmp1
		enddo
	else
		tmp0=a(n-1)
		do i=n-1,3,-2
			wi=wi*wpi+wi
			a(i)=a(i-2)-a(i)
			tmp=a(i)*real(wi)-a(i-1)*aimag(wi)
			a(i-1)=a(i-1)*real(wi)+a(i)*aimag(wi)
			a(i)=tmp
		enddo
		a(1)=2*tmp0
		call fftreal1_rp(a,n,sig)
		do i=0,m-1
			tmp=a(i)+a(n-i-1)
			a(i)=(a(n-i-1)-a(i))/(2*aimag(w))
			a(n-i-1)=0.5*(tmp-a(i))
			a(i)=0.5*(tmp+a(i))
			w=w*wp+w
		enddo
	endif
end
