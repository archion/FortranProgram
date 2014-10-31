module M_rd
	implicit none
	contains
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
		subroutine irandom(i,n)
			real(8) :: rn
			integer :: i,n
			call random_number(rn)
			i=1+int(rn*(n))
		end subroutine
		subroutine fisher_yates_shuffle(a,n)
			integer :: a(:),i,j,tmp,n
			real(8) :: rn
			do i=n,1,-1
				call random_number(rn)
				j=ceiling(rn*i)
				tmp=a(i)
				a(i)=a(j)
				a(j)=tmp
			enddo
		end subroutine
end module
