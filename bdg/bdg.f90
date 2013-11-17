module global
	implicit none
	integer, parameter :: dn(2)=(/24,24/),dn2=dn(1)*dn(2),imp=dn2/2+dn(1)/2
	integer :: latt(dn2,4,3)
	real(8), parameter :: t(3)=(/1d0,-0.25d0,0d0/),nf=0.85d0,U=2.44d0,V=1d0,Vimp=1d0,pi=3.14159265359d0,cvg=1e-5,bt=1e5
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
module lapack
	implicit none
	contains
	subroutine czheev(jobz,uplo,a,w,info)
		implicit none
		character :: jobz,uplo
		integer :: n,lda,lwork,info
		real(8) :: w(:)
		complex(8) :: a(:,:)
		real(8), allocatable :: rwork(:)
		complex(8), allocatable :: work(:)
		n=size(w)
		lda=n
		lwork=3*n
		allocate(work(lwork),rwork(3*n-2))
		call zheev (jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
	end subroutine
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
200	format(24(e15.6))
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
	complex(8) :: dt(dn2,4),dtp(dn2,4),H(2*dn2,2*dn2)
	real(8) :: n(dn2,2),np(dn2,2),n1,wide,sp0,sp,sa,sb,E(2*dn2)
	logical :: flaga,flagb
	wide=0.5d0
	sp0=sp+wide
	do 
		sa=sp
		sb=sp
		flaga=.true.
		flagb=.true.
		do
			sp=0.5d0*(sa+sb)
			call EU(dt,n,sp,H,E)
			call order(H,E,np,dtp)
			n1=1.0/dn2*sum(np)
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
		wide=max(abs(sp0-sp),100*cvg)
		sp0=sp
		if(sum(abs(dt-dtp))/dn2<cvg) then
			exit
		endif
		dt=dtp
		n=np
		write(*,*)"*",dt(1,1)
	enddo 
end
subroutine EU(dt,n,sp,H,E)
	use lapack
	use global
	implicit none
	complex(8) :: H(2*dn2,2*dn2),dt(dn2,4)
	real(8) :: E(2*dn2),n(dn2,2),sp
	integer :: i,j,info
	H=0d0
	do i=1,dn2
		do j=1,4
			! nearest neighbor pairs
			H(i,latt(i,j,1))=-t(1)
			H(i,latt(i,j,1)+dn2)=dt(i,j)
			H(i+dn2,latt(i,j,1)+dn2)=t(1)
			H(i+dn2,latt(i,j,1))=dconjg(dt(i,j))
			! next-nearest neighbor pairs
			H(i,latt(i,j,2))=-t(2)
			H(i+dn2,latt(i,j,2)+dn2)=t(2)
		enddo
		! on site
		H(i,i)=U*n(i,2)-sp
		H(i+dn2,i+dn2)=-U*n(i,1)+sp
	enddo
	! impure
	H(imp,imp)=H(imp,imp)+Vimp
	H(imp+dn2,imp+dn2)=H(imp+dn2,imp+dn2)+Vimp
	call czheev("V","U",H,E,info)
end
subroutine order(H,E,np,dtp)
	use global
	implicit none
	real(8) :: np(dn2,2),f(2*dn2),th(2*dn2),E(2*dn2)
	complex(8) :: dtp(dn2,4),H(2*dn2,2*dn2)
	integer :: i,j,n
	np=0d0
	dtp=0d0
	f=1d0/(exp(bt*E)+1d0)
	th=tanh(0.5d0*bt*E)
	do i=1,dn2
		do n=1,2*dn2
			np(i,1)=np(i,1)+H(i,n)*dconjg(H(i,n))*f(n)
			np(i,2)=np(i,2)+H(i+dn2,n)*dconjg(H(i+dn2,n))*(1d0-f(n))
			do j=1,4
				dtp(i,j)=dtp(i,j)+0.25d0*V*(H(i,n)*dconjg(H(latt(i,j,1)+dn2,n))+H(latt(i,j,1),n)*dconjg(H(i+dn2,n)))*th(n)
			enddo
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

	
