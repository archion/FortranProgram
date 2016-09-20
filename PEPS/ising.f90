
include "tensor.f90"

program main
	use M_tensor
	use M_matrix
	implicit none
	integer, parameter :: chi=6,d=2,N=3000
	integer :: i,j,i0,i1
	real(8) :: DJ=1d0,g=0.5d0,delta=0.005
	type(t_tensor) :: H,expH,A(0:1),ld(0:1),tmp,tmp1
	real(8), allocatable :: s(:),U(:,:),V(:,:)

	call mkl_set_num_threads(1) 
	call A(0)%new([d,chi,chi],"d0 j0 k0","r")
	call A(1)%new([d,chi,chi],"d1 j1 k1","r")
	call ld(0)%new([chi,chi],"l0 r0","0")
	call ld(1)%new([chi,chi],"l1 r1","0")

	call random_number(ld(0)%T(ld(0)%get_idx([([i,i],i=1,chi)])))
	call random_number(ld(1)%T(ld(1)%get_idx([([i,i],i=1,chi)])))
	
	call H%new([d,d,3],"d0 u0 x0","0")
	call tmp%new([d,d,3],"u1 d1 x1","0")
	H%T(H%get_idx([1,1,1, 2,2,1, 1,1,2, 2,2,2, 1,2,3, 2,1,3]))=[1d0,1d0,DJ,-DJ,-g*0.5d0,-g*0.5d0]
	tmp%T(tmp%get_idx([1,2,1, 2,1,1, 1,1,2, 2,2,2, 1,1,3, 2,2,3]))=[-g*0.5d0,-g*0.5d0,1d0,-1d0,1d0,1d0]
	call H%dot(H,tmp,["x0","x1"],1)
	call H%get_mat("d0",2,U)
	allocate(s(size(U,1)))
	call syev(U,s,"V")
	s=exp(-s*delta)
	U=matmul(matmul(U,diag(s)),transpose(U))
	expH=H
	call expH%get_tensor("d0",2,U)

	do i=0,N-1
		i0=mod(i,2)
		i1=mod(i+1,2)

		call tmp%dot(ld(i0),A(i0),["r"//to_char(i0),"j"//to_char(i0)],1)
		call tmp%dot(tmp,ld(i1),["k"//to_char(i0),"l"//to_char(i1)],1)
		call tmp%dot(tmp,A(i1),["r"//to_char(i1),"j"//to_char(i1)],1)
		call tmp%dot(tmp,ld(i0),["k"//to_char(i1),"l"//to_char(i0)],1)


		expH%clock=(-1)**i0
		call tmp%dot(tmp,expH,["d"//to_char(i1),"u"//to_char(i1)],2)

		call tmp%svd("l"//to_char(i0),2,U,s,V)

		s(:chi)=s(:chi)/sqrt(sum(s(:chi)**2))

		call A(i0)%get_tensor("j"//to_char(i0),2,U(:,:chi))
		call A(i1)%get_tensor("j"//to_char(i1),1,V(:chi,:))
		call ld(i1)%get_tensor("l"//to_char(i1),1,diag(s(:chi)))

		ld(i0)%T(ld(i0)%get_idx([([i,i],i=1,chi)]))=1d0/ld(i0)%T(ld(i0)%get_idx([([i,i],i=1,chi)]))
		call A(i0)%dot(ld(i0),A(i0),["r"//to_char(i0),"j"//to_char(i0)],1)
		call A(i0)%change_label(["l"//to_char(i0),"j"//to_char(i0)])
		call A(i1)%dot(A(i1),ld(i0),["k"//to_char(i1),"l"//to_char(i0)],1)
		call A(i1)%change_label(["r"//to_char(i0),"k"//to_char(i1)])
		ld(i0)%T(ld(i0)%get_idx([([i,i],i=1,chi)]))=1d0/ld(i0)%T(ld(i0)%get_idx([([i,i],i=1,chi)]))
	enddo

	write(*,*)-log(sum(tmp%T**2))/delta/2


	!call tmp%dot(ld(i0),A(i0),["r"//to_char(i0),"j"//to_char(i0)],1)
	!call tmp%dot(tmp,ld(i1),["k"//to_char(i0),"l"//to_char(i1)],1)
	!call tmp%dot(tmp,A(i1),["r"//to_char(i1),"j"//to_char(i1)],1)
	!call tmp%dot(tmp,ld(i0),["k"//to_char(i1),"l"//to_char(i0)],1)

	!if(i0==1) H%clock=-1 

	!call tmp1%dot(tmp,H,["d"//to_char(i1),"u"//to_char(i1)],2)
	!tmp%clock=-1
	!call tmp1%dot(tmp1,tmp,["d0","d0"],4)
	!write(*,*)tmp1%T

end program
