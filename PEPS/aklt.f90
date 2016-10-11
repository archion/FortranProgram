include "tensor.f90"
program main
	use M_tensor
	implicit none
	integer :: i,N=6
	type(t_tensor) :: M,A
	real(8), allocatable :: s(:)
	real(8), allocatable :: U(:,:),V(:,:)

	call mkl_set_num_threads(1)

	call M%new(["i1","j","k"],[3,2,2],"0")
	M%T(M%get_idx([1,2,1, 2,1,1, 2,2,2, 3,1,2],["i1","j","k"]))=[-sqrt(2d0/3d0),-sqrt(1d0/3d0),sqrt(1d0/3d0),sqrt(2d0/3d0)]
	call M%reorder([2,1,3])

	A=M
	do i=2,N
		call M%change_label(["i"//to_char(i-1),"i"//to_char(i)])
		A=dot(A,M,["k","j"])
	enddo

	call A%svd(["j  ",("i"//to_char(i),i=1,N/2)],[character::],U,s,V)
	s=s**2/sum(s**2)
	write(*,*)s
	write(*,*)size(s)
end program
