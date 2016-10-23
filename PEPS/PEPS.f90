include "tensor.f90"
include "list.f90"
module M_fig
	use M_tensor
	use M_list
	implicit none
	type t_dat
		type(t_tensor) :: T
		type(t_dat), pointer :: link1
		type(t_dat), pointer :: link2
		integer(8) :: cost(2)
		integer :: cfg=0
		logical :: flag=.false.
	contains
		final :: delete_dat
	end type
contains
	subroutine delete_dat(self)
		type(t_dat) :: self
		if(.not.self%flag) then
			call self%T%clear()
			nullify(self%link1,self%link2)
		endif
	end subroutine
	function get_dat(self) result(rt)
		type(t_node) :: self
		type(t_dat), pointer :: rt
		associate(p=>self%dat)
			select type(p)
			type is(t_dat)
				rt => p
			end select
		end associate
	end function
	function put_dat(self) result(rt)
		type(t_dat), target :: self
		class(*), pointer :: rt
		rt => self
	end function
	recursive function fig_contract(que) result(rt)
		type(t_dat) :: que
		type(t_tensor) :: rt
		if(.not.associated(que%T%rc)) then
			if(.not.associated(que%link1%T%rc)) then
				que%link1%T=fig_contract(que%link1)
				que%link1%T%is_return=.true.
			endif
			if(.not.associated(que%link2%T%rc)) then
				que%link2%T=fig_contract(que%link2)
				que%link2%T%is_return=.true.
			endif
			rt=dot(que%link1%T,que%link2%T)
		else
			rt=que%T
		endif
		rt%is_return=.true.
	end function
	function merge_node(dat1,dat2,dat) result(flag)
		type(t_dat), pointer :: dat1,dat2
		type(t_dat), pointer :: dat
		logical :: flag
		logical, allocatable :: flag1(:),flag2(:)
		integer :: cost,i,j
		integer(8), allocatable :: shap(:)
		flag=.true.
		dat%cfg=ior(dat1%cfg,dat2%cfg)
		if(popcnt(dat%cfg)==(popcnt(dat1%cfg)+popcnt(dat2%cfg))) then
			allocate(flag1(size(dat1%T%shap)),flag2(size(dat2%T%shap)))
			cost=1
			flag1=.true.
			flag2=.true.
			do i=1,size(dat1%T%shap)
				do j=1,size(dat2%T%shap)
					if(dat1%T%label(i)==dat2%T%label(j)) then
						flag1(i)=.false.
						flag2(j)=.false.
						cost=cost*dat1%T%shap(i)
					endif
				enddo
			enddo
			call allocate(dat%T%shap,[pack(dat1%T%shap,flag1),pack(dat2%T%shap,flag2)])
			allocate(shap(size(dat%T%shap)))
			shap=dat%T%shap
			dat%T%label=[pack(dat1%T%label,flag1),pack(dat2%T%label,flag2)]
			dat%cost(1)=product(shap)
			dat%cost(2)=max(dat%cost(2),dat%cost(1)+dat1%cost(2)+dat1%cost(2))
			dat%cost(1)=dat%cost(1)*cost+dat1%cost(1)+dat2%cost(1)
			dat%link1 => dat1
			dat%link2 => dat2
			dat%flag=.false.
		else
			flag=.false.
		endif
	end function
	subroutine contract_order(que)
		type(t_list) :: que
		type(t_list) :: S(que%num)
		integer :: n,i,j,c,d,k,minx
		integer(8) :: u(2),uc,u0
		type(t_dat), pointer :: dat,dat1,dat2,dat3
		type(t_node), pointer :: pos1,pos2,pos3
		integer :: cm=1
		n=que%num
		minx=10000
		pos1 => que%head
		do i=1,n
			if(i/=1) pos1 => pos1%next
			dat => get_dat(pos1)
			dat%cfg=ibset(0,i-1)
			if(minx>minval(dat%T%shap)) minx=minval(dat%T%shap)
			dat%cost=[0,product(dat%T%shap)]
			call S(1)%insert(put_dat(dat))
		enddo
		do i=1,n
			pos1 => S(i)%head
			do j=1,S(i)%num
				if(j/=1) pos1 => pos1%next
				dat => get_dat(pos1)
				dat%flag=.false.
			enddo
		enddo
		nullify(dat)
		uc=1; u(1)=0
		do while(S(n)%num==0)
			u(2)=100000000
			do c=2,n
				do d=1,c/2
					pos1 => S(d)%head
					do i=1,S(d)%num
						if(i/=1) pos1 => pos1%next
						dat1 => get_dat(pos1)
						pos2 => S(c-d)%head
					o:	do j=1,S(c-d)%num
							if(j/=1) pos2 => pos2%next
							if((d==(c-d)).and.(i>=j)) cycle
							dat2 => get_dat(pos2)
							allocate(dat)
							if(merge_node(dat1,dat2,dat)) then
								if(dat1%flag.or.dat2%flag) then
									u0=0
								else
									u0=u(1)
								endif
								if(dat%cost(cm)>uc.and.dat%cost(cm)<u(2)) then
									u(2)=dat%cost(cm)
								elseif(dat%cost(cm)<=uc.and.dat%cost(cm)>u0) then
									pos3 => S(c)%head
									do k=1,S(c)%num
										if(k/=1) pos3 => pos3%next
										dat3 => get_dat(pos3)
										if(dat%cfg==dat3%cfg) then
											if(dat%cost(cm)<dat3%cost(cm)) then
												dat3%flag=.true.
												dat3%cost=dat%cost
												dat3%link1 => dat%link1
												dat3%link2 => dat%link2
												deallocate(dat)
											endif
											cycle o
										endif
									enddo
									call S(c)%insert(put_dat(dat))
								else
									deallocate(dat)
								endif
							endif
						enddo o
					enddo
				enddo
			enddo
			u(1)=uc
			uc=max(u(2),minx*uc)
			do i=1,n
				pos1 => S(i)%head
				do j=1,S(i)%num
					if(j/=1) pos1 => pos1%next
					dat1 => get_dat(pos1)
					dat1%flag=.false.
				enddo
			enddo
		enddo
		pos1 => S(1)%head
		do i=1,S(1)%num
			if(i/=1) pos1 => pos1%next
			dat1 => get_dat(pos1)
			dat1%flag=.true.
		enddo
		pos1 => S(n)%head
		dat => get_dat(S(n)%head)
		do i=1,S(n)%num
			if(i/=1) pos1 => pos1%next
			dat1 => get_dat(pos1)
			write(*,*)dat1%cost
			if(dat%cost(cm)>dat1%cost(cm)) then
				dat => dat1
			endif
		enddo
		write(*,*)"the cost is: ",dat%cost
		call get_que(que,dat)
	end subroutine
	recursive subroutine get_que(que,dat)
		type(t_list) :: que
		type(t_dat), pointer :: dat
		if(dat%flag) return
		if(.not.dat%link1%flag) then
			call get_que(que,dat%link1)
		endif
		if(.not.dat%link2%flag) then
			call get_que(que,dat%link2)
		endif
		dat%flag=.true.
		call que%insert(put_dat(dat))
	end subroutine
end module
program main
	use M_fig
	use M_matrix
	implicit none
	integer, parameter :: d=10,chi=d**2,DD=d**2,p=2
	real(8) :: delta=0.005d0
	real(8), allocatable :: H(:,:)
	real(8) :: He(p**2),HU(size(He),size(He)),Hexp(size(He),size(He))
	type(t_tensor), target :: C1,C2,C3,C4,T1(2),T2(2),T3(2),T4(2),A(2),aa(2),Ap(2),aap(2),tmp1,tmp2
	type(t_tensor) :: E
	type(t_dat), target :: fig(16)
	type(t_list) :: que1
	integer :: i

	call omp_set_nested(.false.)
	!call omp_set_max_active_levels(1)
	call omp_set_dynamic(.false.)
	call mkl_set_dynamic(0)
	call omp_set_schedule(omp_sched_static,0)

	call mkl_set_num_threads(1)
	call omp_set_num_threads(mkl_get_max_threads())

	! Heisenberg
	call tmp1%new(["u0","d0","x"],[d,d,3],"0")
	call tmp2%new(["u1","d1","x"],[d,d,3],"0")
	tmp1%rc%T(tmp1%get_idx([1,1,1, 2,2,1, 1,2,2, 2,1,3]))=[1d0,-1d0,1d0,1d0]
	tmp2%rc%T(tmp2%get_idx([1,1,1, 2,2,1, 2,1,2, 1,2,3]))=[1d0,-1d0,2d0,2d0]

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
		call A(i)%new(["l","r","u","d","p"],[d,d,d,d,p],"r")
		aa(i)=dot(A(i),A(i)%change_label(["l","r","u","d"],["L","R","U","D"]),["p","p"])
		aa(i)=aa(i)%merge_label(["l","L","r","R","u","U","d","D"],["l","r","u","d"])
	enddo
	call A(1)%qr(["u","l","d"],["r","p"],X(1),y(1))
	call A(2)%lq(["l","p"],["u","r","d"],X(1),y(1))


	fig(1)%T=C1%change_label(["r","d"],["1","4"])
	fig(2)%T=T1(2)%change_label(["l","r","d"],["1","2","5"])
	fig(3)%T=T1(1)%change_label(["l","r","d"],["2","3","6"])
	fig(4)%T=C2%change_label(["l","d"],["3","7"])
	fig(5)%T=T4(2)%change_label(["u","r","d"],["4","8","11"])
	fig(6)%T=aa(1)%change_label(["l","r","u","d"],["8","9","5","12"])
	fig(7)%T=aa(2)%change_label(["l","r","u","d"],["9","10","6","13"])
	fig(8)%T=T2(1)%change_label(["l","u","d"],["10","7","14"])
	fig(9)%T=T4(1)%change_label(["u","r","d"],["11","15","18"])
	fig(10)%T=aa(2)%change_label(["l","r","u","d"],["15","16","12","19"])
	fig(11)%T=aa(1)%change_label(["l","r","u","d"],["16","17","13","20"])
	fig(12)%T=T2(2)%change_label(["l","u","d"],["17","14","21"])
	fig(13)%T=C4%change_label(["u","r"],["18","22"])
	fig(14)%T=T3(1)%change_label(["l","r","u"],["22","23","19"])
	fig(15)%T=T3(2)%change_label(["l","r","u"],["23","24","20"])
	fig(16)%T=C3%change_label(["u","l"],["21","24"])

	do i=1,size(fig)
		if(i/=6.and.i/=7) then
			call que1%insert(put_dat(fig(i)))
		else
			call que1%insert(put_dat(fig(i)))
		endif
	enddo
	write(*,*)"start"
	call contract_order(que1)
	write(*,*)que1%num


	!do i=0,N-1
		

		E=fig_contract(get_dat(que1%tail))
		write(*,*)E%label
		!call E6%split_label([""],["e1","E1","e2","E2","e3","E3","e4","E4","e5","E5","e6","E6"])
		
		!A(1)%merge_label([],["l"])
		!A(2)%merge_label([],["r"])
		!App=dot(A(1),A(2),["r","l"])
		!App=dot(App,Hexp,[character(5)::"p1","p2"])
		

		!call A(i)%QR([],Q,R)
		!call X(i)%get_tensor([],Q)
		!E6=dot(E6,X(1),[])
		!E6=dot(E6,X(1),[])
		!E6=dot(E6,X(2),[])
		!E6=dot(E6,X(2),[])

		!do
			!call al(i)%get_tensor([],R)

		!enddo
	!enddo
end program

