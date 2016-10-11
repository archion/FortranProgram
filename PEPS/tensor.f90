module M_tensor
	use lapack95
	use blas95
	use M_utility
	implicit none
	type t_tensor
		character(:), allocatable :: label(:)
		integer, allocatable :: shap(:)
		real(8), allocatable :: T(:) 
	contains
		procedure :: new
		procedure :: get_idx
		procedure :: change_label
		procedure :: reorder
		procedure :: get_mat
		procedure :: get_tensor
		!procedure :: dot
		procedure :: svd => svd_t
		procedure :: qr => qr_t
		procedure :: clear
	end type
	interface dot
		module procedure dot_tt,dot_tm,dot_mt,dot_td,dot_dt
	end interface
	interface allocate
		module procedure allocate_2,allocate_2copy
	end interface
contains
	subroutine new(self,label,shap,flag)
		class(t_tensor) :: self
		character(*) :: label(:)
		integer :: shap(:)
		character, optional :: flag
		integer :: i,j
		if(product(self%shap)/=product(shap)) then
			if(allocated(self%T)) deallocate(self%T)
			allocate(self%T(product(shap)))
		endif
		if(.not.allocated(self%T)) allocate(self%T(product(shap)))
		allocate(self%shap(size(shap)))
		self%shap=shap
		self%label=label
		if(size(self%label)==0) return
		if(present(flag)) then
			select case(flag)
			case("0")
				self%T=0d0
			case("r")
				call random_number(self%T)
			end select
		endif
	end subroutine
	function get_idx(self,idx,label)
		class(t_tensor) :: self
		character(*), optional :: label(:)
		integer :: idx(:)
		integer :: get_idx(size(idx)/size(self%label))
		integer :: prod(size(self%label)),ord(size(self%label))
		integer :: i,j,n
		prod(1)=1
		do i=2,size(self%shap)
			prod(i)=prod(i-1)*self%shap(i-1)
		enddo

		ord=[1:size(self%label)]
		if(present(label)) then
			do i=1,size(label)
				do j=1,size(self%label)
					if(label(i)==self%label(ord(j))) then
						n=ord(j)
						ord(j)=ord(i)
						ord(i)=n
						exit
					endif
				enddo
			enddo
		endif

		n=size(self%label)
		do j=1,size(get_idx)
			get_idx(j)=1
			do i=1,size(self%label)
				get_idx(j)=get_idx(j)+(idx((j-1)*n+i)-1)*prod(ord(i))
			enddo
		enddo
	end function
	subroutine reorder(self,ord)
		class(t_tensor) :: self
		integer :: ord(:)
		integer :: i,map(size(self%T)),prod(size(ord)),n_prod(size(ord))
		n_prod(1)=1
		prod(1)=1
		do i=2,size(self%shap)
			prod(i)=prod(i-1)*self%shap(i-1)
			n_prod(i)=n_prod(i-1)*self%shap(ord(i-1))
		enddo

		prod=prod(ord)
		self%label=self%label(ord)
		self%shap=self%shap(ord)

		do i=1,size(self%T)
			map(i)=1+sum((mod((i-1)/n_prod,self%shap))*prod)
		enddo

		self%T=self%T(map)
	end subroutine
	!subroutine trace(self,label)
	!end subroutine
	subroutine change_label(self,label)
		class(t_tensor) :: self
		character(*) :: label(:)
		integer :: i,j
		do i=1,size(self%shap)
			do j=1,size(label),2
				if(label(j)==self%label(i)) then
					self%label(i)=label(j+1)
					exit
				endif
			enddo
		enddo
	end subroutine
	subroutine get_mat(self,l,r,H)
		class(t_tensor) :: self
		character(*) :: l(:),r(:)
		real(8), allocatable :: H(:,:)
		integer :: ord(size(self%label)),shap(2),i,j,tmp
		logical :: ll(size(l)),lr(size(r))
		ll=.true.
		lr=.true.
		ord=[1:size(self%label)]
		do i=1,size(self%label)
			do j=1,size(self%label)
				if(size(l)/=0.and.i<=size(l)) then
					if(l(i)==self%label(ord(j))) then
						tmp=ord(j)
						ord(j)=ord(i)
						ord(i)=tmp
						ll(i)=.false.
						exit
					endif
				endif
				if(size(r)/=0.and.i>(size(self%label)-size(r))) then
					if(r(i-(size(self%label)-size(r)))==self%label(ord(j))) then
						tmp=ord(j)
						ord(j)=ord(i)
						ord(i)=tmp
						lr(i-(size(self%label)-size(r)))=.false.
						exit
					endif
				endif
			enddo
		enddo
		if(any(ll).or.any(lr)) then
			write(*,*)"error! check label: ",l," | ",r," in ",self%label
			stop 
		endif
		call self%reorder(ord)
		if(size(l)/=0) then
			shap=[product(self%shap(:size(l))),product(self%shap(size(l)+1:))]
		else
			shap=[product(self%shap(:size(self%label)-size(r))),product(self%shap(size(self%label)-size(r)+1:))]
		endif
		if(any((shape(H)-shap)/=0)) then
			if(allocated(H)) deallocate(H)
			allocate(H(shap(1),shap(2)))
		endif
		H=reshape(self%T,shape(H))
	end subroutine
	subroutine get_tensor(self,label,H)
		class(t_tensor) :: self
		character(*) :: label(:)
		real(8) :: H(size(self%T))
		integer :: ord(size(self%shap)),i,j,tmp
		logical :: l(size(self%label))
		l=.true.
		ord=[1:size(self%label)]
		do i=1,size(label)
			do j=1,size(self%label)
				if(label(i)==self%label(ord(j))) then
					tmp=ord(j)
					ord(j)=ord(i)
					ord(i)=tmp
					l(i)=.false.
					exit
				endif
			enddo
		enddo
		if(any(l)) then
			write(*,*)"error! check label: ",label," in ",self%label
			stop 
		endif
		call self%reorder(ord)
		self%T=H
	end subroutine
	!subroutine split_label(A,label)
		!class(t_tensor) :: A
		!character(*) :: label(:)
	!end subroutine
	!subroutine combin_label(A,label)
		!class(t_tensor) :: A
		!character(*) :: label(:)
	!end subroutine
	function contract(A,label)
		class(t_tensor) :: A
		type(t_tensor) :: contract
		character(*) :: label(:)
		integer :: i,j,n
		real(8), allocatable :: T(:),MA(:,:)
		call A%get_mat([label(1::2),label(2::2)],[character::],MA)

		allocate(T(size(MA,2)))

		n=product(A%shap(1:size(label)/2))
		do i=1,size(T)
			T(i)=sum(MA(j::n,i))
		enddo

		call contract%new(A%label(size(label)+1:),A%shap(size(label)+1:))
		contract%T=reshape(T,shape(contract%T))
	end function
	function dot_tt(A,B,label)
		class(t_tensor) :: A,B
		type(t_tensor) :: dot_tt
		character(*) :: label(:)
		integer :: i,j
		real(8), allocatable :: T(:,:),MA(:,:),MB(:,:)
		character(5) :: r_label(size(A%label)+size(B%label)-size(label))
		integer :: r_shap(size(r_label))
		call A%get_mat([character::],label(1::2),MA)
		r_label(:size(A%shap)-size(label)/2)=A%label(:size(A%shap)-size(label)/2)
		r_shap(:size(A%shap)-size(label)/2)=A%shap(:size(A%shap)-size(label)/2)

		call B%get_mat(label(2::2),[character::],MB)
		r_label(size(r_label)-(size(B%label)-size(label)/2)+1:)=B%label(size(label)/2+1:)
		r_shap(size(r_label)-(size(B%label)-size(label)/2)+1:)=B%shap(size(label)/2+1:)

		allocate(T(size(MA,1),size(MB,2)))

		if(size(T)==1) then
			T=sum(MA(1,:)*MB(:,1))
		else
			call gemm(MA,MB,T)
			!T=matmul(MA,MB)
		endif

		call dot_tt%new(r_label,r_shap)
		dot_tt%T=reshape(T,shape(dot_tt%T))
	end function
	function dot_tm(A,MB,label)
		class(t_tensor) :: A
		real(8) :: MB(:,:)
		type(t_tensor) :: dot_tm
		character(*) :: label(:)
		integer :: i,j
		real(8), allocatable :: T(:,:),MA(:,:)
		call A%get_mat([character::],label,MA)

		allocate(T(size(MA,1),size(MB,2)))

		if(size(T)==1) then
			T=sum(MA(1,:)*MB(:,1))
		else
			call gemm(MA,MB,T)
			!T=matmul(MA,MB)
		endif

		call dot_tm%new(A%label,A%shap)
		dot_tm%T=reshape(T,shape(dot_tm%T))
	end function
	function dot_mt(MA,B,label)
		real(8) :: MA(:,:)
		class(t_tensor) :: B
		type(t_tensor) :: dot_mt
		character(*) :: label(:)
		integer :: i,j
		real(8), allocatable :: T(:,:),MB(:,:)
		call B%get_mat(label,[character::],MB)

		allocate(T(size(MA,1),size(MB,2)))

		if(size(T)==1) then
			T=sum(MA(1,:)*MB(:,1))
		else
			call gemm(MA,MB,T)
			!T=matmul(MA,MB)
		endif

		call dot_mt%new(B%label,B%shap)
		dot_mt%T=reshape(T,shape(dot_mt%T))
	end function
	function dot_td(A,dB,label)
		class(t_tensor) :: A
		real(8) :: dB(:)
		type(t_tensor) :: dot_td
		character(*) :: label(:)
		integer :: i,j
		real(8), allocatable :: T(:,:),MA(:,:)
		call A%get_mat([character::],label,MA)

		allocate(T(size(MA,1),size(dB)))

		!$omp parallel do
		do i=1,size(dB)
			T(:,i)=MA(:,i)*dB(i)
		enddo
		!$omp end parallel do

		call dot_td%new(A%label,A%shap)
		dot_td%T=reshape(T,shape(dot_td%T))
	end function
	function dot_dt(dA,B,label)
		real(8) :: dA(:)
		class(t_tensor) :: B
		type(t_tensor) :: dot_dt
		character(*) :: label(:)
		integer :: i,j
		real(8), allocatable :: T(:,:),MB(:,:)
		call B%get_mat(label,[character::],MB)

		allocate(T(size(dA),size(MB,2)))

		!$omp parallel do
		do i=1,size(dA)
			T(i,:)=dA(i)*MB(i,:)
		enddo
		!$omp end parallel do

		call dot_dt%new(B%label,B%shap)
		dot_dt%T=reshape(T,shape(dot_dt%T))
	end function
	subroutine svd_t(self,l,r,U,s,V)
		class(t_tensor) :: self
		character(*) :: l(:),r(:)
		real(8), allocatable :: H(:,:)
		real(8), allocatable, optional :: s(:)
		real(8), allocatable, optional :: U(:,:),V(:,:)
		call self%get_mat(l,r,H)
		if(present(U)) then
			call svd(H,U,s,V)
		else
			call svd(H,s=s)
		endif
	end subroutine
	subroutine svd(H,U,s,V)
		real(8) :: H(:,:)
		real(8), allocatable, optional :: s(:)
		real(8), allocatable, optional :: U(:,:),V(:,:)
		if(present(U)) then
			if(size(s)/=minval(shape(H))) then
				if(allocated(s)) deallocate(s)
				allocate(s(minval(shape(H))))
			endif
			if(any(shape(U)-[size(H,1),size(H,1)]/=0)) then
				if(allocated(U)) deallocate(U)
				allocate(U(size(H,1),size(H,1)))
			endif
			if(any(shape(V)-[size(H,2),size(H,2)]/=0)) then
				if(allocated(V)) deallocate(V)
				allocate(V(size(H,2),size(H,2)))
			endif
			call gesdd(H,s,U,V,"A")
		else
			if(size(s)/=minval(shape(H))) then
				if(allocated(s)) deallocate(s)
				allocate(s(minval(shape(H))))
			endif
			call gesdd(H,s)
		endif
	end subroutine
	subroutine qr_t(self,l,r,MQ,MR)
		class(t_tensor) :: self
		character(*) :: l(:),r(:)
		real(8), allocatable :: H(:,:)
		real(8), allocatable :: MQ(:,:),MR(:,:)
		call self%get_mat(l,r,H)
		call qr(H,MQ,MR)
	end subroutine
	subroutine qr(H,Q,R)
		real(8) :: H(:,:)
		real(8), allocatable :: Q(:,:),R(:,:),tau(:)
		integer :: i,j
		allocate(Q(size(H,1),size(H,1)),R(size(H,1),size(H,2)),tau(minval(shape(H))))
		call geqrf(H,tau)
		do i=1,size(R,1)
			do j=1,size(R,2)
				if(i>j) then
					R(i,j)=0d0
					Q(i,j)=H(i,j)
				else
					R(i,j)=H(i,j)
				endif
			enddo
		enddo
		call orgqr(Q,tau)
	end subroutine
	!subroutine lq(H,Q,L)
		!real(8) :: H(:,:)
		!real(8), allocatable :: Q(:,:),R(:,:),tau(:)
		!integer :: i,j
		!allocate(Q(size(H,1),size(H,1)),R(size(H,1),size(H,2)),tau(minval(shape(H))))
		!call geqrf(H,tau)
		!do i=1,size(R,1)
			!do j=1,size(R,2)
				!if(i>j) then
					!R(i,j)=0d0
					!Q(i,j)=H(i,j)
				!else
					!R(i,j)=H(i,j)
				!endif
			!enddo
		!enddo
		!call orgqr(Q,tau)
	!end subroutine
	subroutine clear(self)
		class(t_tensor) :: self
		if(allocated(self%shap)) deallocate(self%shap,self%T,self%label)
	end subroutine
	subroutine allocate_2(A,shap)
		real(8), allocatable :: A(:,:)
		integer :: shap(2)
		if(any(shape(A)-shap/=0)) then
			if(allocated(A)) deallocate(A)
		else
			return
		endif
		allocate(A(shap(1),shap(2)))
	end subroutine
	subroutine allocate_2copy(A,B)
		real(8), allocatable :: A(:,:)
		real(8) :: B(:,:)
		if(any(shape(A)-shape(B)/=0)) then
			if(allocated(A)) deallocate(A)
		else
			A=B
			return
		endif
		allocate(A(size(B,1),size(B,2)))
		A=B
	end subroutine
end module
