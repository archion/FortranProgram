module global
	implicit none
	integer, parameter :: dn(2)=(/24,24/),dn2=dn(1)*dn(2),imp=dn2/2+dn(1)/2
	integer :: latt(dn2,4,3)
	real(8), parameter :: t(3)=(/1d0,-0.25d0,0d0/),nf=0.85d0,U=2.44d0,V=-1d0,Vimp=1d0,pi=3.14159265359d0,cvg=1e-5,bt=1e5
	complex(8), parameter :: img=(0d0,1d0),ipi=2d0*pi*img
	contains
		subroutine gen_latt_square()
!					  j
!			 1 --  2 --  3 --  4
!			 |     |     |     |           4      3   4
!			 5 --  6 --  7 --  8           |       \ /
!		i	 |     |     |     |        3--0--1     0
!			 9 -- 10 -- 11 -- 12           |       / \
!			 |     |     |     |           2      2   1
!			13 -- 14 -- 15 -- 16
			implicit none
			integer :: ii,i,j,pbcx(-dn(1)+1:2*dn(1)),pbcy(-dn(2)+1:2*dn(2))
			pbcx=(/(i,i=1,dn(1)),(i,i=1,dn(1)),(i,i=1,dn(1))/)
			pbcy=(/(i,i=1,dn(2)),(i,i=1,dn(2)),(i,i=1,dn(2))/)
			do ii=1,dn2
				j=mod(ii-1,dn(1))+1
				i=(ii-1)/dn(1)+1
				latt(ii,1,1)=dn(1)*(i-1)+pbcx(j+1)
				latt(ii,2,1)=dn(1)*(pbcy(i+1)-1)+j
				latt(ii,3,1)=dn(1)*(i-1)+pbcx(j-1)
				latt(ii,4,1)=dn(1)*(pbcy(i-1)-1)+j
				latt(ii,1,2)=dn(1)*(pbcy(i+1)-1)+pbcx(j+1)
				latt(ii,2,2)=dn(1)*(pbcy(i+1)-1)+pbcx(j-1)
				latt(ii,3,2)=dn(1)*(pbcy(i-1)-1)+pbcx(j-1)
				latt(ii,4,2)=dn(1)*(pbcy(i-1)-1)+pbcx(j+1)
				latt(ii,1,3)=dn(1)*(i-1)+pbcx(j+2)
				latt(ii,2,3)=dn(1)*(pbcy(i+2)-1)+j
				latt(ii,3,3)=dn(1)*(i-1)+pbcx(j-2)
				latt(ii,4,3)=dn(1)*(pbcy(i-2)-1)+j
			enddo
		end subroutine gen_latt_square
end module
program main
	use global
	implicit none
	character(25) :: fmat(2)
	integer :: i,j
	real(8) :: n(dn2,2)
	complex(8) :: dt(dn2,4)
	open(unit=10,file='../data/order.dat')
	call gnuplot
	write(fmat(1),*)dn(1)
	call inital(dt,n)
	call bdg(dt,n)
	write(10,"("//TRIM(ADJUSTL(fmat(1)))//"e15.6)")abs(dt(:,1)+dt(:,3)-dt(:,2)-dt(:,4))/4.0d0
	write(10,"(1X/)")
	write(10,"("//TRIM(ADJUSTL(fmat(1)))//"e15.6)")n(:,1)+n(:,2)
	write(10,"(1X/)")
	write(10,"("//TRIM(ADJUSTL(fmat(1)))//"e15.6)")(((-1)**(i+j)*0.5*(n(dn(1)*i+j,1)-n(dn(1)*i+j,2)),j=1,dn(1)),i=0,dn(2)-1)
	close(10)
end
subroutine inital(dt,n)
	use global
	implicit none
	complex(8) :: dt(dn2,4)
	real(8) :: n(dn2,2)
	integer :: i
	call init_random_seed()
	!call random_number(rdom)
	do i=1,dn2
		dt(i,:)=(/0.07d0,-0.07d0,0.07d0,-0.07d0/)
	enddo
	n(:,1)=nf/2d0
	n(:,2)=nf-n(:,1)
	call gen_latt_square()
end
subroutine bdg(dt,n)
	use global
	implicit none
	complex(8) :: dt(dn2,4),dtp(dn2,4),va(2*dn2*13)
	real(8) :: gm,n(dn2,2),np(dn2,2),n1,wide,sp0,sp,sa,sb,a,b,emax,emin
	integer :: ja(2*dn2*13),ia(2*dn2+1),info,r=0
	logical :: flaga,flagb
	gm=0.1d0
	sp=0d0
	wide=0.5d0
	sp0=sp+wide
	do 
		sa=sp
		sb=sp
		flaga=.true.
		flagb=.true.
		r=r+1
		do
			sp=0.5d0*(sa+sb)
			call Ham(dt,n,sp,va,ja,ia)
			call lanmax(va,ja,ia,emax,emin,2*dn2,info)
			a=(emax-emin)/2+gm
			b=(emax+emin)/2
			call rescale(va,ja,ia,a,b)
			call cbv(va,ja,ia,a,b,np,dtp)
			n1=1.0d0/dn2*sum(np)
			if(n1>2.or.n1<0) then
				write(*,*)"the a maybe small"
				gm=gm+0.1d0
				flaga=.true.
				flagb=.true.
				cycle
			endif
			if(abs(n1-nf)<cvg) then
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
			write(*,*)n1,sp
		enddo
		dtp=dtp*V
		wide=max(abs(sp0-sp),100*cvg)
		sp0=sp
		if(sum(abs(dt-dtp))/dn2<cvg.or.r>200) then
			exit
		endif
		dt=dtp
		n=np
		write(*,*)"*",dt(1,1)
	enddo 
end
subroutine Ham(dt,n,sp,va,ja,ia)
	use global
	implicit none
	integer :: i,j,k,m,l,sg,ja(2*dn2*13),ia(2*dn2+1)
	complex(8) :: dt(dn2,4),va(2*dn2*13)
	real(8) :: sp,n(dn2,2)
	m=1
	ia(1)=1
	l=2
	sg=-1
	do k=1,2
		do i=1,dn2
			va(m)=sg*(-U*n(i,k)+sp)
			ja(m)=i+dn2*(k-1)
			if(i==imp) then
				va(m)=va(m)+Vimp
			endif
			m=m+1
			do j=1,4
				va(m)=sg*t(1)
				ja(m)=latt(i,j,1)+dn2*(k-1)
				m=m+1
				va(m)=dt(i,j)
				ja(m)=latt(i,j,1)+dn2*(mod(k,2))
				m=m+1
				va(m)=sg*t(2)
				ja(m)=latt(i,j,2)+dn2*(k-1)
				m=m+1
			enddo
			ia(l)=m
			l=l+1
		enddo
		dt=dconjg(dt)
		sg=-sg
	enddo
	dt=dconjg(dt)
end
subroutine cbv(va,ja,ia,a,b,n,dt)
	use global
	implicit none
	integer, parameter :: cf=1000
	complex(8) :: hn(2*dn2,0:2),va(2*dn2*13),dt(dn2,4)
	integer :: i,j,k,m,md(3),ja(2*dn2*13),ia(2*dn2+1)
	real(8) :: a,b,ac,ld=0.01d0,Tn(cf),n(dn2,2)
	n=0d0
	dt=0d0
	ac=acos(-b/a)
	do m=1,cf
		Tn(m)=-2.0d0*sin(m*ac)/(m*pi)*sinh(ld*(1.0d0-real(m,8)/cf))
	enddo
	Tn=Tn/sinh(ld)
	!$OMP PARALLEL DO REDUCTION(+:n,dt) PRIVATE(hn,md) SCHEDULE(GUIDED)
	do i=1,dn2
		do k=1,2
			hn=0d0
			hn(i+dn2*(k-1),0)=1d0
			n(i,k)=n(i,k)+(1d0-ac/pi)
			do m=1,cf
				md(1)=abs(mod(m-2,3))
				md(2)=mod(m-1,3)
				md(3)=mod(m,3)
				call crsmv(va,ja,ia,hn(:,md(2)),hn(:,md(3)),2*dn2)
				hn(:,md(3))=2d0*hn(:,md(3))-hn(:,md(1))
				n(i,k)=n(i,k)+hn(i+dn2*(k-1),md(3))*Tn(m)
				if(k==2) then
					do j=1,4
						dt(i,j)=dt(i,j)+hn(latt(i,j,1),md(3))*Tn(m)
					enddo
				endif
			enddo
		enddo
	enddo
	!$OMP END PARALLEL DO
	n(:,2)=1d0-n(:,2)
end
subroutine rescale(va,ja,ia,a,b)
	use global
	implicit none
	complex(8) :: va(2*dn2*13)
	integer :: i,j,ja(2*dn2*13),ia(2*dn2+1)
	real(8) :: a,b
	do i=1,2*dn2
		do j=ia(i),ia(i+1)-1
			if(ja(j)==i) then
				va(j)=va(j)-b
				exit
			endif
		enddo
	enddo
	va=va/a
end
subroutine lanmax(va,ja,ia,emax,emin,m,info)
	implicit none
	integer :: m,ja(*),ia(*),md(0:2),i,j,info
	real(8) :: emax,emin,s,work,ap(m),bt(0:m),tmp(0:m),e(m),ep(m)
	complex(8) :: va(*),v(m,0:2)
	!real(8),allocatable :: v(:,:),ap(:),bt(:),tmp(:),e(:),ep(:)
	!allocate(v(m,0:2),ap(m),bt(0:m),tmp(0:m),e(m),ep(m))
	ep(1)=1000
	v(:,0)=0d0
	bt(0)=0d0
	call random_number(tmp(1:m))
	v(:,1)=tmp(1:m)
	s=dot_product(v(:,1),v(:,1))
	v(:,1)=v(:,1)/sqrt(abs(s))
	do i=1,m
		md(2)=mod(i+1,3)
		md(1)=mod(i,3)
		md(0)=mod(i-1,3)
		call crsmv(va,ja,ia,v(:,md(1)),v(:,md(2)),m)
		ap(i)=dot_product(v(:,md(1)),v(:,md(2)))
		e(1:i)=ap(1:i)
		tmp(1:i-1)=bt(1:i-1)
		call dsteqr('n',i,e(1:i),tmp(1:i-1),work,i,work,info)
		! write(10,*)e(1:i)
		if(abs(e(1)-ep(1))<1e-3.or.info/=0) then
			exit
		endif
		ep=e
		v(:,md(2))=v(:,md(2))-ap(i)*v(:,md(1))-bt(i-1)*v(:,md(0))
		s=dot_product(v(:,md(2)),v(:,md(2)))
		bt(i)=sqrt(abs(s))
		v(:,md(2))=v(:,md(2))/bt(i)
	enddo
	emax=e(min(i,m))
	emin=e(1)
	!deallocate(v,ap,bt,tmp,e,ep)
end
subroutine crsmv(va,ja,ia,x,y,m)
	implicit none
	integer :: i,j,m
	integer :: ja(*),ia(*)
	complex(8) :: va(*),y(*),x(*)
	do i=1,m
		y(i)=0d0
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
subroutine gnuplot()
	implicit none
	write(10,"(A)")"set term pngcairo"
	write(10,"(A)")"set output 'order.png'"
	write(10,"(A)")"set pm3d map"
	write(10,"(A)")"set pm3d corners2color c1"
	write(10,"(A)")"set cbrange [:]"
	write(10,"(A)")"set pm3d interpolate 0,0"
	write(10,"(A)")"set size square"
	write(10,"(A)")"set palette rgbformulae 22,13,-31"
	write(10,"(A)")"splot '-' matrix index 0"
	write(10,"(A)")"#data"
end


