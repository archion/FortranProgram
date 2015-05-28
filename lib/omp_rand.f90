module M_omp_rand
	use MersenneTwister_mod
	implicit none
	interface mt_random_number
		module procedure mt_irandom,mt_rrandom
	end interface
contains
	subroutine mt_init_random_seed(rnd,seed)
		integer :: seed
		type(randomNumberSequence) :: rnd
		rnd=new_RandomNumberSequence(seed=seed)
	end subroutine
	function mt_irandom(rnd,n)
		type(randomNumberSequence) :: rnd
		integer :: mt_irandom,n
		mt_irandom=1+int(getRandomReal(rnd)*n)
	end function
	function mt_rrandom(rnd)
		type(randomNumberSequence) :: rnd
		real(8) :: mt_rrandom
		mt_rrandom=getRandomReal(rnd)
	end function
end module
