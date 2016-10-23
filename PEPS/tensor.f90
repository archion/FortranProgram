module M_tensor
	use lapack95
	use blas95
	use M_utility
	implicit none
	type t_rc
		real(8), pointer, contiguous :: T(:)
		integer :: c=0
	end type
	type t_ptensor
		type(t_tensor), pointer :: tg
	end type
	type t_tensor
		character(:), allocatable :: label(:)
		integer, allocatable :: shap(:)
		type(t_rc), pointer :: rc
		type(t_ptensor) :: link(2)
		logical :: is_return=.false.
	contains
		procedure :: new
		procedure :: get_idx
		procedure :: change_label
		procedure :: split_label
		procedure :: merge_label
		procedure :: get_order
		procedure :: reorder
		procedure :: get_mat
		procedure :: get_tensor
		procedure :: contract
		procedure :: svd => svd_t
		procedure :: qr => qr_t
		procedure :: lq => lq_t
		procedure :: clear
		procedure :: insert
		procedure :: remove
		procedure :: equal
		generic :: assignment(=) => equal
	end type
	interface dot
		module procedure dot_tt,dot_tm,dot_mt,dot_td,dot_dt
	end interface
	interface allocate
		module procedure allocate_1i,allocate_1icopy,allocate_1r,allocate_1rcopy,allocate_2r,allocate_2rcopy,allocate_1ip,allocate_1ipcopy,allocate_1rp,allocate_1rpcopy,allocate_2rp,allocate_2rpcopy
	end interface
contains
	subroutine equal(to,from)
		class(t_tensor), intent(inout) :: to
		type(t_tensor) :: from
		if(.not.associated(to%rc,from%rc)) then
			if(associated(to%rc)) then
				if(to%rc%c>1) then
					write(*,*)"equal error"
				else
					deallocate(to%rc%T)
				endif
			endif
		endif
		to%rc => from%rc
		call to%new(from%label,from%shap)
		call from%insert(to)
		if(from%is_return) then
			call from%clear()
		endif
	end subroutine
	subroutine new(self,label,shap,flag)
		class(t_tensor) :: self
		character(*) :: label(:)
		integer :: shap(:)
		character, optional :: flag
		real(8), allocatable :: T(:)
		integer :: i,j
		if(.not.associated(self%rc)) then
			allocate(self%rc)
			self%rc%c=1
		endif
		call allocate(self%rc%T,product(shap))
		call allocate(self%shap,shap)
		self%label=label
		if(size(self%label)==0) return
		if(present(flag)) then
			select case(flag)
			case("0")
				self%rc%T=0d0
			case("r")
				call random_number(self%rc%T)
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
		class(t_tensor), target :: self
		integer :: ord(:)
		integer :: i,map(size(self%rc%T)),prod(size(ord)),n_prod(size(ord))
		type(t_tensor), pointer :: node
		if(all(ord-[1:size(ord)]==0)) return
		n_prod(1)=1
		prod(1)=1
		do i=2,size(self%shap)
			prod(i)=prod(i-1)*self%shap(i-1)
			n_prod(i)=n_prod(i-1)*self%shap(ord(i-1))
		enddo

		prod=prod(ord)
		self%label=self%label(ord)
		self%shap=self%shap(ord)

		do i=1,size(self%rc%T)
			map(i)=1+sum((mod((i-1)/n_prod,self%shap))*prod)
		enddo

		self%rc%T=self%rc%T(map)

		do i=1,2
			node => self
			do
				if(associated(node%link(i)%tg)) then
					node => node%link(i)%tg
					if(size(node%shap)==size(ord)) then
						node%label=node%label(ord)
						node%shap=node%shap(ord)
					else
						stop "alias has different label size"
					endif
				else
					exit
				endif
			enddo
		enddo
	end subroutine
	function split_label(self,l_from,l_to) result(rt)
		class(t_tensor), target :: self
		type(t_tensor) :: rt
		character(*) :: l_from(:),l_to(:)
		integer :: shap(size(self%shap)-size(l_from)+size(l_to)),i,j,n,m
		character(max(len(self%label),len(l_to))) :: label(size(shap))
		integer, allocatable :: ord(:)
		m=size(l_to)/size(l_from)
		call self%get_order(l_from,[character::],ord)
		ord(ord)=[1:size(ord)]
		n=0
		do i=1,size(self%shap)
			if(ord(i)<=size(l_from)) then
				label(n+1:n+m)=l_to((ord(i)-1)*m+1:ord(i)*m)
				shap(n+1:n+m)=nint(self%shap(i)**(1d0/m))
				n=n+m
			else
				n=n+1
				label(n)=self%label(i)
				shap(n)=self%shap(i)
			endif
		enddo

		allocate(rt%rc)
		rt%rc%T => self%rc%T
		call rt%new(label,shap)
		call self%insert(rt)
		rt%is_return=.true.
	end function
	function merge_label(self,l_from,l_to) result(rt)
		class(t_tensor), target :: self
		type(t_tensor), target :: rt
		character(*) :: l_from(:),l_to(:)
		integer :: shap(size(self%shap)-size(l_from)+size(l_to)),i,j,n,m
		character(max(len(self%label),len(l_to))) :: label(size(shap))
		integer, allocatable :: ord(:)
		m=size(l_from)/size(l_to)
		call self%get_order(l_from,[character::],ord)
		do i=1,size(l_to)
			if(any((ord((i-1)*m+1:i*m)-ord((i-1)*m+1))-[0:m-1]/=0)) then
				exit
			endif
		enddo
		if(i==size(l_to)+1) then
			ord(ord)=[1:size(ord)]
			n=0
			do i=1,size(self%shap)
				if(ord(i)<=size(l_from)) then
					if(mod(ord(i),m)==1) then
						n=n+1
						label(n)=l_to((ord(i)-1)/m+1)
						shap(n)=product(self%shap(i:i+m-1))
					endif
				else
					n=n+1
					label(n)=self%label(i)
					shap(n)=self%shap(i)
				endif
			enddo
		else
			call self%reorder(ord)
			do i=1,size(shap)
				if(i<=size(l_to)) then
					label(i)=l_to(i)
					shap(i)=product(self%shap((i-1)*m+1:i*m))
				else
					label(i)=self%label(size(l_from)-size(l_to)+i)
					shap(i)=self%shap(size(l_from)-size(l_to)+i)
				endif
			enddo
		endif
		
		allocate(rt%rc)
		rt%rc%T => self%rc%T
		call rt%new(label,shap)
		call self%insert(rt)
		rt%is_return=.true.
	end function
	function change_label(self,l_from,l_to) result(rt)
		class(t_tensor), target :: self
		type(t_tensor) :: rt
		character(*) :: l_from(:),l_to(:)
		integer :: i
		character(max(len(self%label),len(l_to))) :: label(size(self%shap))
		integer, allocatable :: ord(:)
		call self%get_order(l_from,[character::],ord)
		label=self%label
		do i=1,size(l_from)
			label(ord(i))=l_to(i)
		enddo

		allocate(rt%rc)
		rt%rc%T => self%rc%T
		call rt%new(label,self%shap)
		call self%insert(rt)
		rt%is_return=.true.
	end function
	subroutine get_order(self,l,r,ord)
		class(t_tensor) :: self
		character(*) :: l(:),r(:)
		logical :: ll(size(l)),lr(size(r))
		integer, allocatable :: ord(:)
		integer :: i,j,tmp
		allocate(ord(size(self%label)))
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
			i=RAISEQQ(SIG$ABORT)
		endif
	end subroutine
	subroutine get_mat(self,l,r,H)
		class(t_tensor) :: self
		character(*) :: l(:),r(:)
		real(8), allocatable :: H(:,:)
		integer, allocatable :: ord(:)
		integer :: shap(2)
		logical :: ll(size(l)),lr(size(r))
		call self%get_order(l,r,ord)
		call self%reorder(ord)
		if(size(l)/=0) then
			shap=[product(self%shap(:size(l))),product(self%shap(size(l)+1:))]
		else
			shap=[product(self%shap(:size(self%label)-size(r))),product(self%shap(size(self%label)-size(r)+1:))]
		endif
		call allocate(H,shap)
		H=reshape(self%rc%T,shape(H))
	end subroutine
	subroutine get_tensor(self,label,H)
		class(t_tensor) :: self
		character(*) :: label(:)
		real(8) :: H(size(self%rc%T))
		integer, allocatable :: ord(:)
		call self%get_order(label,[character::],ord)
		call self%reorder(ord)
		self%rc%T=H
	end subroutine
	function contract(self,label) result(rt)
		class(t_tensor) :: self
		type(t_tensor), target :: rt
		character(*) :: label(:)
		integer :: i,j,n
		real(8), pointer, contiguous :: T(:),MA(:,:)
		integer, allocatable :: ord(:)

		call self%get_order([label(1::2),label(2::2)],[character::],ord)
		call self%reorder(ord)

		n=product(self%shap(1:size(label)/2))

		MA(1:n**2,1:product(self%shap(size(label)+1:))) => self%rc%T

		allocate(T(size(MA,2)))

		do i=1,size(T)
			T(i)=0d0
			do j=1,n
				T(i)=T(i)+MA(j*n-n+j,i)
			enddo
		enddo

		allocate(rt%rc)
		rt%rc%c=1
		rt%rc%T => T
		call rt%new(self%label(size(label)+1:),self%shap(size(label)+1:))
		rt%is_return=.true.
		if(self%is_return) then
			call self%clear()
		endif
	end function
	function dot_tt(A,B,label) result(rt)
		class(t_tensor) :: A,B
		type(t_tensor), target :: rt
		character(*), optional :: label(:)
		character(max(len(A%label),len(B%label))) :: label_(size(A%label)+size(A%label))
		integer :: i,j
		real(8), pointer, contiguous :: T(:,:),MA(:,:),MB(:,:)
		real(8), allocatable, target :: tmp(:)
		integer :: n
		integer, allocatable :: orda(:),ordb(:)
		if(.not.present(label)) then
			n=0
			do i=1,size(A%label)
				if(any(B%label==A%label(i))) then
					label_(n+1)=A%label(i)
					label_(n+2)=A%label(i)
					n=n+2
				endif
			enddo
			n=n/2
		else
			n=size(label)/2
			label_(1:n*2)=label
		endif
		call A%get_order(label_(1:n*2:2),[character::],orda)
		call B%get_order(label_(2:n*2:2),[character::],ordb)
		if(associated(A%rc%T,B%rc%T)) then
			if(all(orda-ordb==0)) then
				call A%reorder(orda)
			else
				allocate(tmp(size(A%rc%T)))
				tmp=A%rc%T
				if(A%is_return) then
					call B%remove(A)
					allocate(A%rc)
					A%rc%T => tmp
					A%rc%c=1
				else
					call A%remove(B)
					allocate(B%rc)
					B%rc%T => tmp
					B%rc%c=1
				endif
				call A%reorder(orda)
				call B%reorder(ordb)
			endif
		else
			call A%reorder(orda)
			call B%reorder(ordb)
		endif
		MA(1:max(product(A%shap(1:n)),1),1:max(product(A%shap(n+1:)),1)) => A%rc%T
		MB(1:max(product(B%shap(1:n)),1),1:max(product(B%shap(n+1:)),1)) => B%rc%T

		allocate(T(size(MA,2),size(MB,2)))

		if(size(T)==1) then
			T=sum(MA(:,1)*MB(:,1))
		else
			call gemm(MA,MB,T,transa="t")
		endif

		allocate(rt%rc)
		rt%rc%c=1
		rt%rc%T(1:size(T)) => T
		call rt%new([character(5)::A%label(n+1:),B%label(n+1:)],[A%shap(n+1:),B%shap(n+1:)])
		rt%is_return=.true.

		if(A%is_return) then
			call A%clear()
		endif
		if(B%is_return) then
			call B%clear()
		endif
	end function
	function dot_tm(A,MB,label) result(rt)
		class(t_tensor) :: A
		real(8) :: MB(:,:)
		type(t_tensor), target :: rt
		character(*) :: label(:)
		integer :: n
		real(8), pointer, contiguous :: T(:,:),MA(:,:)
		integer, allocatable :: ord(:)
		n=size(A%shap)-size(label)
		call A%get_order([character::],label,ord)
		call A%reorder(ord)

		MA(1:product(A%shap(1:n)),1:product(A%shap(n+1:))) => A%rc%T

		allocate(T(size(MA,1),size(MB,2)))

		if(size(T)==1) then
			T=sum(MA(1,:)*MB(:,1))
		else
			call gemm(MA,MB,T)
			!T=matmul(MA,MB)
		endif

		allocate(rt%rc)
		rt%rc%c=1
		rt%rc%T(1:size(T)) => T
		call rt%new(A%label,A%shap)
		rt%is_return=.true.
		if(A%is_return) then
			call A%clear()
		endif
	end function
	function dot_mt(MA,B,label) result(rt)
		real(8) :: MA(:,:)
		class(t_tensor) :: B
		type(t_tensor), target :: rt
		character(*) :: label(:)
		integer :: n
		real(8), pointer, contiguous :: T(:,:),MB(:,:)
		integer, allocatable :: ord(:)
		n=size(label)
		call B%get_order(label,[character::],ord)
		call B%reorder(ord)

		MB(1:product(B%shap(1:n)),1:product(B%shap(n+1:))) => B%rc%T

		allocate(T(size(MA,1),size(MB,2)))

		if(size(T)==1) then
			T=sum(MA(1,:)*MB(:,1))
		else
			call gemm(MA,MB,T)
			!T=matmul(MA,MB)
		endif


		allocate(rt%rc)
		rt%rc%c=1
		rt%rc%T(1:size(T)) => T
		call rt%new(B%label,B%shap)
		rt%is_return=.true.
		if(B%is_return) then
			call B%clear()
		endif
	end function
	function dot_td(A,dB,label) result(rt)
		class(t_tensor) :: A
		real(8) :: dB(:)
		type(t_tensor), target :: rt
		character(*) :: label(:)
		integer :: i,n
		real(8), pointer, contiguous :: T(:,:),MA(:,:)
		integer, allocatable :: ord(:)
		n=size(label)
		call A%get_order(label,[character::],ord)
		call A%reorder(ord)

		MA(1:product(A%shap(1:n)),1:product(A%shap(n+1:))) => A%rc%T

		allocate(T(size(dB),size(MA,2)))

		!$omp parallel do
		do i=1,size(T,2)
			T(:,i)=MA(:,i)*dB(:)
		enddo
		!$omp end parallel do

		allocate(rt%rc)
		rt%rc%c=1
		rt%rc%T(1:size(T)) => T
		call rt%new(A%label,A%shap)
		rt%is_return=.true.
		if(A%is_return) then
			call A%clear()
		endif
	end function
	function dot_dt(dA,B,label) result(rt)
		real(8) :: dA(:)
		class(t_tensor) :: B
		type(t_tensor), target :: rt
		character(*) :: label(:)
		integer :: i,n
		real(8), pointer, contiguous :: T(:,:),MB(:,:)
		integer, allocatable :: ord(:)
		n=size(label)
		call B%get_order(label,[character::],ord)
		call B%reorder(ord)

		MB(1:product(B%shap(1:n)),1:product(B%shap(n+1:))) => B%rc%T

		allocate(T(size(dA),size(MB,2)))

		!$omp parallel do
		do i=1,size(T,2)
			T(:,i)=dA(:)*MB(:,i)
		enddo
		!$omp end parallel do

		allocate(rt%rc)
		rt%rc%c=1
		rt%rc%T(1:size(T)) => T
		call rt%new(B%label,B%shap)
		rt%is_return=.true.
		if(B%is_return) then
			call B%clear()
		endif

	end function
	subroutine svd_t(self,l,r,U,s,V)
		class(t_tensor) :: self
		character(*) :: l(:),r(:)
		real(8), allocatable :: H(:,:)
		real(8), allocatable, optional :: s(:)
		real(8), allocatable, optional :: U(:,:),V(:,:)
		call self%get_mat(l,r,H)
		call svd(H,U,s,V)
	end subroutine
	subroutine svd(H,U,s,V)
		real(8) :: H(:,:)
		real(8), allocatable, optional :: s(:)
		real(8), allocatable, optional :: U(:,:),V(:,:)
		if(.not.present(U)) then
			call allocate(s,minval(shape(H)))
			call gesdd(H,s)
		elseif(.not.present(s)) then
		else
			call allocate(s,minval(shape(H)))
			call allocate(U,shape(H))
			call allocate(V,[size(H,2),size(H,2)])
			call gesdd(H,s,U,V,"A")
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
		real(8), allocatable :: Q(:,:),R(:,:)
		real(8) :: tau(minval(shape(H)))
		integer :: i,j
		call allocate(Q,[size(H,1),minval(shape(H))])
		call allocate(R,[minval(shape(H)),size(H,2)])
		call geqrf(H,tau)
		do i=1,size(H,1)
			do j=1,size(H,2)
				if(i>j) then
					if(i<=size(R,1)) R(i,j)=0d0
					Q(i,j)=H(i,j)
				else
					R(i,j)=H(i,j)
				endif
			enddo
		enddo
		call orgqr(Q,tau)
	end subroutine
	subroutine lq(H,Q,L)
		real(8) :: H(:,:)
		real(8), allocatable :: Q(:,:),L(:,:)
		real(8) :: tau(minval(shape(H)))
		integer :: i,j
		call allocate(Q,[minval(shape(H)),size(H,2)])
		call allocate(L,[size(H,1),minval(shape(H))])
		call gelqf(H,tau)
		do i=1,size(H,1)
			do j=1,size(H,2)
				if(i<j) then
					if(j<=size(L,2)) L(i,j)=0d0
					Q(i,j)=H(i,j)
				else
					L(i,j)=H(i,j)
				endif
			enddo
		enddo
		call orglq(Q,tau)
	end subroutine
	subroutine lq_t(self,l,r,MQ,ML)
		class(t_tensor) :: self
		character(*) :: l(:),r(:)
		real(8), allocatable :: H(:,:)
		real(8), allocatable :: MQ(:,:),ML(:,:)
		call self%get_mat(l,r,H)
		call lq(H,MQ,ML)
	end subroutine
	subroutine remove(self,t)
		class(t_tensor), target :: self
		type(t_tensor), target :: t
		type(t_tensor), pointer :: node
		integer :: i
		self%rc%c=self%rc%c-1
		do i=1,2
			node => self
			do
				if(associated(node,t)) then
					if(associated(node%link(1)%tg).and.associated(node%link(2)%tg)) then
						node%link(1)%tg%link(2)%tg => node%link(2)%tg
						node%link(2)%tg%link(1)%tg => node%link(1)%tg
					elseif(associated(node%link(1)%tg)) then
						nullify(node%link(1)%tg%link(2)%tg)
					else
						nullify(node%link(2)%tg%link(1)%tg)
					endif
					nullify(node%link(1)%tg)
					nullify(node%link(2)%tg)
					nullify(node%rc)
					return
				endif
				if(associated(node%link(i)%tg)) then
					node => node%link(i)%tg
				else
					exit
				endif
			enddo
		enddo
		write(*,*)"can not find node"
		i=RAISEQQ(SIG$ABORT)
	end subroutine
	subroutine insert(self,t)
		class(t_tensor), target :: self
		type(t_tensor), target :: t
		type(t_tensor), pointer :: node
		integer :: i,n
		do i=1,2
			node => self
			do
				if(associated(node,t)) then
					return
				endif
				if(associated(node%link(i)%tg)) then
					node => node%link(i)%tg
				else
					exit
				endif
			enddo
		enddo
		t%rc => self%rc
		self%rc%c=self%rc%c+1
		t%link(2)%tg => self%link(2)%tg
		self%link(2)%tg => t
		t%link(1)%tg => self
	end subroutine
	subroutine clear(self)
		class(t_tensor) :: self
		if(allocated(self%shap)) deallocate(self%label,self%shap)
		if(associated(self%rc)) then
			if(self%rc%c==1) then
				deallocate(self%rc%T)
			else
				call self%remove(self)
				nullify(self%rc)
			endif
		endif
	end subroutine
	subroutine ar_naupd(H,E,V,lr)
		integer, parameter :: ncv = 20, nev = 1
		real(8) :: H(:,:),V(:)
		complex(8) :: E
		real(8) :: sigmar, sigmai, tol
		integer :: iparam(11), ipntr(14), info, ierr, ido ,n
		real(8) :: d(ncv,3), resid(size(H,1)), vec(size(H,1),ncv), workd(3*size(H,1)), workev(3*ncv), workl(3*ncv*ncv+6*ncv)
		logical :: select(ncv)
		integer :: i
		character(1), optional :: lr
		n=size(H,1)
		tol=0d0
		ido=0
		info=0
		iparam(1)=1
		iparam(3)=size(V)
		iparam(7)=1
		do
			call dnaupd ( ido, 'I', n, 'LM', nev, tol, resid, ncv, vec, n, iparam, ipntr, workd, workl, size(workl), info )
			if (ido==-1 .or. ido==1) then
				if(present(lr)) then
					if(lr=="l") then
						call gemv(H,workd(ipntr(1):ipntr(1)+n-1),workd(ipntr(2):ipntr(2)+n-1),trans="t")
						!workd(ipntr(2):ipntr(2)+n-1)=matmul(workd(ipntr(1):ipntr(1)+n-1),H)
						!!$omp parallel do
						!do i=0,n-1
							!workd(ipntr(2)+i)=sum(workd(ipntr(1):ipntr(1)+n-1)*H(:,i+1))
						!enddo
						!!$omp end parallel do
						cycle
					endif
				endif
				call gemv(H,workd(ipntr(1):ipntr(1)+n-1),workd(ipntr(2):ipntr(2)+n-1))
				!workd(ipntr(2):ipntr(2)+n-1)=matmul(H,workd(ipntr(1):ipntr(1)+n-1))
				!!$omp parallel do
				!do i=0,n-1
					!workd(ipntr(2)+i)=sum(H(i+1,:)*workd(ipntr(1):ipntr(1)+n-1))
				!enddo
				!!$omp end parallel do
			else 
				if(info==0) then
					call dneupd ( .true., 'A', select, d, d(1,2), vec, n, sigmar, sigmai, workev, 'I', n, 'LM', nev, tol, resid, ncv, vec, n, iparam, ipntr, workd, workl, size(workl), ierr )
					exit
				else
					stop "arpack info error"
				endif
			endif
		enddo
		E=cmplx(d(1,1),d(1,2))
		V=vec(:,1)
	end subroutine
	subroutine allocate_1r(A,shap)
		real(8), allocatable :: A(:)
		integer :: shap
		if(allocated(A)) then
			if(size(A)-shap/=0) then
				deallocate(A)
			else
				return
			endif
		endif
		allocate(A(shap))
	end subroutine
	subroutine allocate_1rcopy(A,B)
		real(8), allocatable :: A(:)
		real(8) :: B(:)
		if(allocated(A)) then
			if(size(A)-size(B)/=0) then
				deallocate(A)
			else
				A=B
				return
			endif
		endif
		allocate(A(size(B)))
		A=B
	end subroutine
	subroutine allocate_2r(A,shap)
		real(8), allocatable :: A(:,:)
		integer :: shap(2)
		if(allocated(A)) then
			if(any(shape(A)-shap/=0)) then
				deallocate(A)
			else
				return
			endif
		endif
		allocate(A(shap(1),shap(2)))
	end subroutine
	subroutine allocate_2rcopy(A,B)
		real(8), allocatable :: A(:,:)
		real(8) :: B(:,:)
		if(allocated(A)) then
			if(any(shape(A)-shape(B)/=0)) then
				deallocate(A)
			else
				A=B
				return
			endif
		endif
		allocate(A(size(B,1),size(B,2)))
		A=B
	end subroutine
	subroutine allocate_1i(A,shap)
		integer, allocatable :: A(:)
		integer :: shap
		if(allocated(A)) then
			if(size(A)-shap/=0) then
				deallocate(A)
			else
				return
			endif
		endif
		allocate(A(shap))
	end subroutine
	subroutine allocate_1icopy(A,B)
		integer, allocatable :: A(:)
		integer :: B(:)
		if(allocated(A)) then
			if(size(A)-size(B)/=0) then
				deallocate(A)
			else
				A=B
				return
			endif
		endif
		allocate(A(size(B)))
		A=B
	end subroutine
	subroutine allocate_1rp(A,shap)
		real(8), pointer :: A(:)
		integer :: shap
		if(associated(A)) then
			if(size(A)-shap/=0) then
				deallocate(A)
			else
				return
			endif
		endif
		allocate(A(shap))
	end subroutine
	subroutine allocate_1rpcopy(A,B)
		real(8), pointer :: A(:)
		real(8) :: B(:)
		if(associated(A)) then
			if(size(A)-size(B)/=0) then
				deallocate(A)
			else
				A=B
				return
			endif
		endif
		allocate(A(size(B)))
		A=B
	end subroutine
	subroutine allocate_2rp(A,shap)
		real(8), pointer :: A(:,:)
		integer :: shap(2)
		if(associated(A)) then
			if(any(shape(A)-shap/=0)) then
				deallocate(A)
			else
				return
			endif
		endif
		allocate(A(shap(1),shap(2)))
	end subroutine
	subroutine allocate_2rpcopy(A,B)
		real(8), pointer :: A(:,:)
		real(8) :: B(:,:)
		if(associated(A)) then
			if(any(shape(A)-shape(B)/=0)) then
				deallocate(A)
			else
				A=B
				return
			endif
		endif
		allocate(A(size(B,1),size(B,2)))
		A=B
	end subroutine
	subroutine allocate_1ip(A,shap)
		integer, pointer :: A(:)
		integer :: shap
		if(associated(A)) then
			if(size(A)-shap/=0) then
				deallocate(A)
			else
				return
			endif
		endif
		allocate(A(shap))
	end subroutine
	subroutine allocate_1ipcopy(A,B)
		integer, pointer :: A(:)
		integer :: B(:)
		if(associated(A)) then
			if(size(A)-size(B)/=0) then
				deallocate(A)
			else
				A=B
				return
			endif
		endif
		allocate(A(size(B)))
		A=B
	end subroutine
end module
