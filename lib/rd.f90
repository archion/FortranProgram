module M_rd
	implicit none
	interface random_number
		module procedure ncrandom,mcrandom,irandom,acrandom
	end interface
	interface init_random_seed
		module procedure init_random_seed_time
	end interface
	interface binning
		module procedure nbinning,vbinning
	end interface
	interface fisher_yates_shuffle
		module procedure sfisher_yates_shuffle,vfisher_yates_shuffle
	end interface
contains
	subroutine init_random_seed_time()
		integer :: i, n, clock
		integer, dimension(:), allocatable :: seed
		call random_seed(size = n)
		allocate(seed(n))
		call system_clock(count=clock)
		seed = clock + 37 * (/ (i - 1, i = 1, n) /)
		call random_seed(put = seed)
		deallocate(seed)
	end subroutine
	subroutine irandom(i,n)
		real(8) :: rn
		integer :: i,n
		call random_number(rn)
		i=1+int(rn*(n))
	end subroutine
	subroutine vfisher_yates_shuffle(a)
		integer :: a(:,:),i,j,tmp(size(a,2))
		integer :: n
		real(8) :: rn
		n=size(a,1)
		do i=n,1,-1
			call random_number(rn)
			j=ceiling(rn*i)
			tmp=a(i,:)
			a(i,:)=a(j,:)
			a(j,:)=tmp
		enddo
	end subroutine
	subroutine sfisher_yates_shuffle(a)
		integer :: a(:),i,j,tmp
		integer :: n
		real(8) :: rn
		n=size(a)
		do i=n,1,-1
			call random_number(rn)
			j=ceiling(rn*i)
			tmp=a(i)
			a(i)=a(j)
			a(j)=tmp
		enddo
	end subroutine
	subroutine ncrandom(r)
		complex(8) :: r
		real(8) :: a,b
		call random_number(a) 
		call random_number(b) 
		r=cmplx(a,b) 
	end subroutine
	subroutine mcrandom(r)
		complex(8) :: r(:,:)
		real(8) :: a(size(r,1),size(r,2)),b(size(r,1),size(r,2))
		call random_number(a) 
		call random_number(b) 
		r=cmplx(a,b) 
	end subroutine
	subroutine acrandom(r)
		complex(8) :: r(:)
		real(8) :: a(size(r)),b(size(r))
		call random_number(a) 
		call random_number(b) 
		r=cmplx(a,b) 
	end subroutine
	subroutine vbinning(i,j,a,x)
		real(8) :: a(:,:),x(:)
		integer :: i,j,n,k
		!i=1
		!j=size(a)/2+1
		!a=0d0
		a(j,:)=x
		a(1,:)=(a(1,:)*(i-1)+a(j,:)**2)/i
		n=i
		k=1
		!write(*,"(es9.1$)")a
		!write(*,"(x)")
		do 
			if(mod(n,2)==0) then
				a(j-1,:)=0.5d0*(a(j-1,:)+a(j,:))
				k=k+1
				a(k,:)=(a(k,:)*(i/2**(k-1)-1)+a(j-1,:)**2)/(i/2**(k-1))
				a(j,:)=0d0
				j=j-1
				n=n/2
				!write(*,"(es9.1$)")a
				!write(*,"(x)")
			else
				j=j+1
				exit
			endif
		enddo
		i=i+1
		!(sqrt((a(i)-a(size(a)/2))/(2**(n-i+1)-1)),sqrt((a(i)-a(size(a)/2))/(2**(n-i+1)-1)*sqrt(2d0/(2**(n-i+1)-1))),i=1,size(a)/2)
	end subroutine
	subroutine nbinning(i,j,a,x)
		real(8) :: a(:),x
		integer :: i,j,n,k
		!i=1
		!j=size(a)/2+1
		!a=0d0
		a(j)=x
		a(1)=(a(1)*(i-1)+a(j)**2)/i
		n=i
		k=1
		!write(*,"(es9.1$)")a
		!write(*,"(x)")
		do 
			if(mod(n,2)==0) then
				a(j-1)=0.5d0*(a(j-1)+a(j))
				k=k+1
				a(k)=(a(k)*(i/2**(k-1)-1)+a(j-1)**2)/(i/2**(k-1))
				a(j)=0d0
				j=j-1
				n=n/2
				!write(*,"(es9.1$)")a
				!write(*,"(x)")
			else
				j=j+1
				exit
			endif
		enddo
		i=i+1
		!(sqrt((a(i)-a(size(a)/2))/(2**(n-i+1)-1)),sqrt((a(i)-a(size(a)/2))/(2**(n-i+1)-1)*sqrt(2d0/(2**(n-i+1)-1))),i=1,size(a)/2)
	end subroutine
end module
