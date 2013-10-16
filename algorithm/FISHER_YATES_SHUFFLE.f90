program main
	implicit none
	integer :: a(100),i
	open(unit=10,file="../data/output1.dat")
	do i=1,size(a)
		a(i)=i
	enddo
	call init_random_seed()
	call fisher_yates_shuffle(a,size(a))
	write(10,"(i10)")a
	close(10)
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
eND SUBROUTINE
