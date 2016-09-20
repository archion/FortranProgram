include "tensor.f90"
program main
	use M_tensor
	implicit none
	integer :: i,N
	type(t_tensor) :: M,A
	real(8), allocatable :: s(:)
	real(8), allocatable :: U(:,:),V(:,:)

	call mkl_set_num_threads(1)

	call M%new([3,2,2],"i1 j k","0")
	M%T(M%get_idx([1,2,1]))=-sqrt(2d0/3d0)
	M%T(M%get_idx([2,1,1]))=-sqrt(1d0/3d0)
	M%T(M%get_idx([2,2,2]))=sqrt(1d0/3d0)
	M%T(M%get_idx([3,1,2]))=sqrt(2d0/3d0)

	N=6

	A=M
	do i=2,N
		call M%change_label(["i"//to_char(i-1),"i"//to_char(i)])
		call A%dot(A,M,["k","j"],1)
	enddo

	call A%svd("j",N/2+1,U,s,V)
	s=s**2/sum(s**2)
	write(*,*)s
	write(*,*)size(s)
end program
