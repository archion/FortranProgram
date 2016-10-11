include "tensor.f90"
module M_fig
	use M_tensor
	implicit none
	type t_ptensor
		type(t_tensor), pointer :: p
	end type
	type t_fig
		type(t_ptensor) :: T(16)
		integer :: bond(2,24)
		character(5) :: label(2,24)
	end type
contains
	subroutine fig_contract(fig,ibond,T)
		type(t_fig) :: fig
		integer :: ibond(:)
		type(t_tensor) :: T
		integer :: i,j,idx(2)
		logical :: flag(size(fig%T))
		flag=.false.
		do i=1,size(ibond)
			idx=fig%bond(:,ibond(i))
			do j=1,size(fig%T(idx(1))%p%label)
				fig%T(idx(1))%p%label(j)=fig%T(idx(1))%p%label(j)(1:1)//to_char(idx(1))
			enddo
			do j=1,size(fig%T(idx(2))%p%label)
				fig%T(idx(2))%p%label(j)=fig%T(idx(2))%p%label(j)(1:1)//to_char(idx(2))
			enddo

			if(all(flag(idx)==.false.)) then
				T=dot(fig%T(idx(1))%p,fig%T(idx(2))%p,[character(5)::trim(fig%label(1,ibond(i)))//to_char(idx(1)),trim(fig%label(2,ibond(i)))//to_char(idx(2))])
			elseif(all(flag(idx)==.true.)) then
				T=contract(T,[character(5)::trim(fig%label(1,ibond(i)))//to_char(idx(1)),trim(fig%label(2,ibond(i)))//to_char(idx(2))])
			else
				if(flag(idx(1)))  then
					T=dot(T,fig%T(idx(2))%p,[character(5)::trim(fig%label(1,ibond(i)))//to_char(idx(1)),trim(fig%label(2,ibond(i)))//to_char(idx(2))])
				else
					T=dot(T,fig%T(idx(1))%p,[character(5)::trim(fig%label(2,ibond(i)))//to_char(idx(2)),trim(fig%label(1,ibond(i)))//to_char(idx(1))])
				endif
			endif

			do j=1,size(fig%T(idx(1))%p%label)
				fig%T(idx(1))%p%label(j)=fig%T(idx(1))%p%label(j)(1:1)
			enddo
			do j=1,size(fig%T(idx(2))%p%label)
				fig%T(idx(2))%p%label(j)=fig%T(idx(2))%p%label(j)(1:1)
			enddo
			flag(idx)=.true.
		enddo

	end subroutine
end module
program main
	use M_fig
	use M_matrix
	implicit none
	integer, parameter :: chi=60,d=6,DD=d**2,p=2
	real(8) :: delta=0.005d0
	real(8), allocatable :: H(:,:)
	real(8) :: He(p**2),HU(size(He),size(He))
	type(t_tensor), target :: C1,C2,C3,C4,T1(2),T2(2),T3(2),T4(2),A(2),aa(2),Ap(2),aap(2)
	type(t_tensor) :: E6
	type(t_fig)  :: fig1


	! Heisenberg
	call tmp1%new(["u0","d0","x"],[d,d,3],"0")
	call tmp2%new(["u1","d1","x"],[d,d,3],"0")
	tmp1%T(tmp1%get_idx([1,1,1, 2,2,1, 1,2,2, 2,1,3]))=[1d0,-1d0,1d0,1d0]
	tmp2%T(tmp2%get_idx([1,1,1, 2,2,1, 2,1,2, 1,2,3]))=[1d0,-1d0,2d0,2d0]

	tmp1=dot(tmp1,tmp2,["x","x"])
	call tmp1%get_mat(["u0","u1"],["d0","d1"],H)
	HU=H
	call syev(HU,He,"V")

	Hexp=matmul(matmul(HU,diag(exp(-He*delta))),transpose(HU))


	call C1%new(["r","d"],[chi,chi],"r")
	call C2%new(["l","d"],[chi,chi],"r")
	call C3%new(["l","u"],[chi,chi],"r")
	call C4%new(["r","u"],[chi,chi],"r")

	call T1(1)%new(["l","r","d"],[chi,chi,DD],"r")
	call T1(2)%new(["l","r","d"],[chi,chi,DD],"r")
	call T2(1)%new(["u","d","l"],[chi,chi,DD],"r")
	call T2(2)%new(["u","d","l"],[chi,chi,DD],"r")
	call T3(1)%new(["l","r","u"],[chi,chi,DD],"r")
	call T3(2)%new(["l","r","u"],[chi,chi,DD],"r")
	call T4(1)%new(["u","d","r"],[chi,chi,DD],"r")
	call T4(2)%new(["u","d","r"],[chi,chi,DD],"r")



	do i=1,2
		call A(i)%new([character(5)::"l","r","u","d","p"//to_char(i)],[d,d,d,d,p],"r")
		At(i)=A(i)
		call At(i)%change_label(["l","L","r","R","u","U","d","D"])
		aa(i)=dot(A(i),At(i),[character(5)::"p"//to_char(i),"p"//to_char(i)])
		call aa(i)%merge_label(["l","L","r","R","u","U","d","D"],["l","r","u","d"])
	enddo


	fig1%T(1)%p => C1
	fig1%T(2)%p => T1(2)
	fig1%T(3)%p => T1(1)
	fig1%T(4)%p => C2
	fig1%T(5)%p => T4(2)
	fig1%T(6)%p => aa(1)
	fig1%T(7)%p => aa(2)
	fig1%T(8)%p => T2(1)
	fig1%T(9)%p => T4(1)
	fig1%T(10)%p => aa(2)
	fig1%T(11)%p => aa(1)
	fig1%T(12)%p => T2(2)
	fig1%T(13)%p => C4
	fig1%T(14)%p => T3(1)
	fig1%T(15)%p => T3(2)
	fig1%T(16)%p => C3

	fig1%bond=reshape([1,2, 2,3, 3,4, 1,5, 2,6, 3,7, 4,8, 5,6, 6,7, 7,8, 5,9, 6,10, 7,11, 8,12, 9,10, 10,11,&
		11,12, 9,13, 10,14, 11,15, 12,16, 13,14, 14,15, 15,16  ],shape(fig1%bond))
	fig1%label=reshape(["r","l", "r","l", "r","l", "d","u", "d","u", "d","u", "d","u", "r","l", "r","l", "r",&
		"l", "d","u", "d","u", "d","u", "d","u", "r","l", "r","l", "r","l", "d","u", "d","u", "d","u", "d","u",&
		"r","l", "r","l", "r","l"], shape(fig1%label))

	do i=0,N-1

		call fig_contract(fig1,[1,2,3,4,7,11,14,15,16,17,18,19,20,21,22,23,24],E6)
		call E6%split_label([""],["e1","E1","e2","E2","e3","E3","e4","E4","e5","E5","e6","E6"])
		
		A(1)%merge_label([],["l"])
		A(2)%merge_label([],["r"])
		App=dot(A(1),A(2),["r","l"])
		App=dot(App,Hexp,[character(5)::"p1","p2"])
		

		call A(i)%QR([],Q,R)
		call X(i)%get_tensor([],Q)
		E6=dot(E6,X(1),[])
		E6=dot(E6,X(1),[])
		E6=dot(E6,X(2),[])
		E6=dot(E6,X(2),[])

		do
			call al(i)%get_tensor([],R)

		enddo
	enddo
end program
