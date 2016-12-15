include "tensor.f90"
program main
	use M_tensor
	implicit none
	integer :: i,N=6
	type(t_tensor) :: M,A
	real(8), allocatable :: s(:)
	complex(8), allocatable :: U(:,:),V(:,:)

	call mkl_set_num_threads(1)

	call M%new(["d","l","r"],[3,2,2],flag="0")
	M%rc%T(M%get_idx([1,2,1, 2,1,1, 2,2,2, 3,1,2]))=[-sqrt(2d0/3d0),-sqrt(1d0/3d0),sqrt(1d0/3d0),sqrt(2d0/3d0)]

	A=M%clone()
	A=A%change_label(["d"],["d1"])
	do i=2,N
		A=dot(A,M%change_label(["d"],[character(4)::"d"//to_char(i)]),["r","l"])
	enddo

	call A%svd([character(5)::"l",("d"//to_char(i),i=1,N/2)],[character::],U,s,V)
	s=s**2/sum(s**2)
	write(*,*)s
	write(*,*)size(s)
end program
