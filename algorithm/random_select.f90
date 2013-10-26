program main
	implicit none
	integer :: cfg(100),i,cont(100)
	open(unit=10,file="random.dat")
	cont=0
	call random_seed()
	do i=1,1000000
		call random_init(cfg,50,100)
		cont=cont+cfg
	enddo
	write(10,"(I7)")cont
end
subroutine random_select(cfg,m,n)
!select a single, random combination of values
!invented by Robert Floyd
!Jon Bentley's Programming Pearls column "A sample of Brilliance"
	implicit none
	integer :: m,n,i,ir,cfg(n)
	real(8) :: r
	cfg=0
	do i=n-m+1,n
		call random_number(r)
		ir=1+int(r*(i))
		if(cfg(ir)==0) then
			cfg(ir)=1
		else
			cfg(i)=1
		endif
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
