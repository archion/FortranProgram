module M_tensor
	use lapack95
	use blas95
	use M_utility
	use M_matrix
	use M_rd
	implicit none
	integer, allocatable :: grp(:,:,:)
	type t_rc
		!real(8), pointer, contiguous :: T(:)
		complex(8), pointer, contiguous :: T(:)
		integer :: c=1
	end type
	type t_ptensor
		type(t_tensor), pointer :: tg
	end type
	type t_tensor
		character(:), allocatable :: label(:)
		integer, allocatable :: shap(:)
		integer, allocatable :: stp(:)
		type(t_rc), pointer :: rc
		type(t_ptensor) :: link(2)
		logical :: is_return=.false.
		logical :: is_conjg=.false.
	contains
		procedure :: new
		procedure :: clone
		procedure :: set_conjg
		procedure :: get_idx
		procedure :: get_value
		procedure :: get_ldx
		procedure :: change_label
		procedure :: split_label
		procedure :: merge_label
		procedure :: get_order
		generic :: reorder => reorder_l,reorder_n
		procedure :: reorder_l
		procedure :: reorder_n
		procedure :: get_mat
		procedure :: get_tensor
		procedure :: contract
		procedure :: svd => svd_t
		generic :: qr => qr_t,qr_2t
		procedure :: qr_t
		procedure :: qr_2t
		generic :: lq => lq_t,lq_2t
		procedure :: lq_t
		procedure :: lq_2t
		generic :: fact => fact_t,fact_a
		procedure :: fact_a
		procedure :: fact_p
		procedure :: fact_t
		procedure :: clear
		procedure :: insert
		procedure :: remove
		procedure :: equal
		generic :: assignment(=) => equal
		final :: delete,delete_all,delete_all2
	end type
	interface dot
		module procedure dot_tt,dot_tm,dot_mt,dot_td,dot_dt
	end interface
	interface ar_naupd
		module procedure ar_naupd_r,ar_naupd_c
	end interface
	interface allocate
		module procedure allocate_1i,allocate_1icopy,allocate_1r,allocate_1rcopy,allocate_2r,allocate_2rcopy,allocate_1ip,allocate_1ipcopy,allocate_1rp,allocate_1rpcopy,allocate_2rp,allocate_2rpcopy,allocate_1c,allocate_1ccopy,allocate_2c,allocate_2ccopy,allocate_1cp,allocate_1cpcopy,allocate_2cp,allocate_2cpcopy
	end interface
contains
	subroutine equal(to,from)
		class(t_tensor), intent(inout) :: to
		type(t_tensor) :: from
		if(.not.associated(to%rc,from%rc)) then
			call to%clear()
			!if(associated(to%rc)) then
				!if(to%rc%c>1) then
					!!write(*,*)"equal error",RAISEQQ(SIG$ABORT)
					!call to%remove(to)
				!else
					!deallocate(to%rc%T)
					!deallocate(to%rc)
				!endif
			!endif
		endif
		if(allocated(from%stp)) then
			call to%new(from%label,from%shap,from%stp,flag="n")
		else
			call to%new(from%label,from%shap,flag="n")
		endif
		to%rc => from%rc
		to%is_conjg=from%is_conjg
		call from%insert(to)
		if(from%is_return) then
			call from%clear()
		endif
	end subroutine
	subroutine new(self,label,shap,stp,flag)
		class(t_tensor) :: self
		character(*) :: label(:)
		integer :: shap(:)
		integer, optional :: stp(:)
		character, optional :: flag
		integer :: i,j
		integer, allocatable :: map(:),grp(:,:)
		if(product(int(shap,kind=8))>2147483647) then
			write(*,*)"err: the size of tensor is large than integer(4), try integer(8) as default integer"
		endif
		if(associated(self%rc)) call self%remove(self)
		allocate(self%rc)
		if(present(flag)) then
			select case(flag)
			case("n")
			case("0")
				call allocate(self%rc%T,product(shap))
				self%rc%T=0d0
			case("r")
				call allocate(self%rc%T,product(shap))
				if(present(stp)) then
					call smerge(gstp(stp),map,grp)
					!call random_number(self%rc%T)
					do i=1,grp(1,2)
						call random_number(self%rc%T(map(i)))
					enddo
					self%rc%T(map(grp(1,2)+1:))=0d0
				else
					call random_number(self%rc%T)
				endif
			end select
		else
			call allocate(self%rc%T,product(shap))
		endif
		call allocate(self%shap,shap)
		self%label=label
		if(present(stp)) then
			call allocate(self%stp,stp)
		else
			if(allocated(self%stp)) deallocate(self%stp)
		endif
		self%is_conjg=.false.
		self%is_return=.false.
	end subroutine
	function get_ldx(self,label) result(rt)
		class(t_tensor) :: self
		character(*) :: label(:)
		integer :: rt(size(label))
		integer, allocatable :: ord(:)
		call self%get_order(label,[character::],ord)
		rt=ord(1:size(rt))
	end function
	!function get_stp(stp) result(rt)
		!integer :: stp(:)
		!integer, allocatable :: rt(:)
		!integer :: i
		!do i=1,size(stp)
			!do
				!stp
			!enddo
		!enddo
	!end function
	function set_conjg(self) result(rt)
		class(t_tensor) :: self
		type(t_tensor) :: rt
		if(allocated(self%stp)) then
			call rt%new(self%label,self%shap,self%stp,flag="n")
		else
			call rt%new(self%label,self%shap,flag="n")
		endif
		rt%rc => self%rc
		rt%is_conjg=(.not.self%is_conjg)
		rt%is_return=.true.
		call self%insert(rt)
		if(self%is_return) then
			call self%clear()
		endif
	end function
	function clone(self) result(rt)
		class(t_tensor) :: self
		type(t_tensor) :: rt
		if(allocated(self%stp)) then
			call rt%new(self%label,self%shap,self%stp)
		else
			call rt%new(self%label,self%shap)
		endif
		if(self%is_conjg) then
			rt%rc%T=conjg(self%rc%T)
		else
			rt%rc%T=self%rc%T
		endif
		rt%is_return=.true.
		if(self%is_return) then
			call self%clear()
		endif
	end function
	function get_value(self,idx) result(rt)
		class(t_tensor) :: self
		integer :: idx
		complex(8) :: rt
		if(self%is_conjg) then
			rt=conjg(self%rc%T(idx))
		else
			rt=self%rc%T(idx)
		endif
		if(self%is_return) then
			call self%clear()
		endif
	end function
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
	subroutine reorder_l(self,label)
		class(t_tensor), target :: self
		character(*) :: label(:)
		integer, allocatable :: ord(:)
		call self%get_order(label,[character::],ord)
		call self%reorder_n(ord)
	end subroutine
	subroutine reorder_n(self,ord)
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
		if(allocated(self%stp)) self%stp=self%stp(ord)

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
						if(allocated(node%stp)) node%stp=node%stp(ord)
					else
						write(*,*)"alias has different label size",self%label,"|",node%label
						write(*,*)RAISEQQ(SIG$ABORT)
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
		integer :: shap(size(self%shap)-size(l_from)+size(l_to)),stp(size(shap)),i,j,n,m
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
				if(allocated(self%stp)) then
					stp(n+1:n+m)=gstp(self%stp(i:i))
				endif
				n=n+m
			else
				n=n+1
				label(n)=self%label(i)
				shap(n)=self%shap(i)
				if(allocated(self%stp)) then
					stp(n)=self%stp(i)
				endif
			endif
		enddo

		if(allocated(self%stp)) then
			call rt%new(label,shap,stp,flag="n")
		else
			call rt%new(label,shap,flag="n")
		endif
		rt%rc => self%rc
		call self%insert(rt)
		rt%is_return=.true.
		rt%is_conjg=self%is_conjg
		if(self%is_return) then
			call self%clear()
		endif
	end function
	function merge_label(self,l_from,l_to) result(rt)
		class(t_tensor), target :: self
		type(t_tensor), target :: rt
		character(*) :: l_from(:),l_to(:)
		integer :: shap(size(self%shap)-size(l_from)+size(l_to)),stp(size(shap)),i,j,n,m
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
						if(allocated(self%stp)) then
							stp(n)=sum(self%stp(i:i+m-1)*[(10**i,i=0,m-1)])
						endif
					endif
				else
					n=n+1
					label(n)=self%label(i)
					shap(n)=self%shap(i)
					if(allocated(self%stp)) then
						stp(n)=self%stp(i)
					endif
				endif
			enddo
		else
			call self%reorder(ord)
			do i=1,size(shap)
				if(i<=size(l_to)) then
					label(i)=l_to(i)
					shap(i)=product(self%shap((i-1)*m+1:i*m))
					if(allocated(self%stp)) then
						stp(i)=sum(self%stp((i-1)*m+1:i*m)*[(10**i,i=0,m-1)])
					endif
				else
					label(i)=self%label(size(l_from)-size(l_to)+i)
					shap(i)=self%shap(size(l_from)-size(l_to)+i)
					if(allocated(self%stp)) then
						stp(i)=self%stp(size(l_from)-size(l_to)+i)
					endif
				endif
			enddo
		endif
		
		if(allocated(self%stp)) then
			call rt%new(label,shap,stp,flag="n")
		else
			call rt%new(label,shap,flag="n")
		endif
		rt%rc => self%rc
		call self%insert(rt)
		rt%is_return=.true.
		rt%is_conjg=self%is_conjg
		if(self%is_return) then
			call self%clear()
		endif
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

		if(allocated(self%stp)) then
			call rt%new(label,self%shap,self%stp,flag="n")
		else
			call rt%new(label,self%shap,flag="n")
		endif
		rt%rc => self%rc
		call self%insert(rt)
		rt%is_return=.true.
		rt%is_conjg=self%is_conjg
		if(self%is_return) then
			call self%clear()
		endif
	end function
	subroutine get_order(self,l,r,ord)
		class(t_tensor) :: self
		character(*) :: l(:),r(:)
		logical :: ll(size(l)),lr(size(r))
		integer, allocatable :: ord(:)
		integer :: i,j,tmp
		call allocate(ord,size(self%label))
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
		!real(8), allocatable :: H(:,:)
		complex(8), allocatable :: H(:,:)
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
		if(self%is_conjg) then
			H=reshape(conjg(self%rc%T),shape(H))
		else
			H=reshape(self%rc%T,shape(H))
		endif
		if(self%is_return) then
			call self%clear()
		endif
	end subroutine
	subroutine get_tensor(self,label,H)
		class(t_tensor) :: self
		character(*) :: label(:)
		!real(8) :: H(size(self%rc%T))
		complex(8) :: H(size(self%rc%T))
		complex(8), pointer, contiguous :: T(:)
		complex(8), target :: zero(0)
		integer, allocatable :: ord(:)
		call self%get_order(label,[character::],ord)
		T => self%rc%T
		self%rc%T => zero
		call self%reorder(ord)
		self%rc%T => T
		if(self%is_conjg.and.self%rc%c/=1) then
			self%rc%T=conjg(H)
		else
			self%is_conjg=.false.
			self%rc%T=H
		endif
	end subroutine
	function gstp(stp) result(rt)
		integer :: stp(:)
		integer, allocatable :: rt(:)
		integer :: stp_(size(stp))
		integer :: n,i
		stp_=stp
		n=0
		do i=1,size(stp_)
			do
				n=n+1
				if(stp_(i)<10) then
					exit
				else
					stp_(i)=stp_(i)/10
				endif
			enddo
		enddo
		allocate(rt(n))
		n=0
		stp_=stp
		do i=1,size(stp_)
			do 
				n=n+1
				rt(n)=mod(stp_(i),10)
				if(stp_(i)<10) then
					exit
				else
					stp_(i)=stp_(i)/10
				endif
			enddo
		enddo
	end function
	subroutine smerge(stp,map,mgrp)
		integer :: stp(:)
		integer, allocatable, optional :: map(:)
		integer, allocatable, optional :: mgrp(:,:)
		integer :: i,k,n(size(stp)),m(size(stp)),pos(2),prod(size(stp)+1),idx(size(stp))

		if(size(stp)==0) then
			if(present(map)) then
				allocate(map(1))
				map(1)=1
			endif
			if(present(mgrp)) then
				allocate(mgrp(2,2))
				mgrp(:,1)=[1,-1]
				mgrp(1,2)=1
				mgrp(2,2)=1
			endif
			return
		endif
		
		if(any(stp>size(grp,3))) then
			write(*,*)RAISEQQ(SIG$ABORT)
		endif
		
		k=size(stp)
		m=grp(1,2,stp)
		n=grp(2,2,stp)

		prod(1)=1
		do i=1,k
			prod(i+1)=prod(i)*(m(i)+n(i))
		enddo

		pos=1
		do i=1,2**k
			idx=mod(i-1,2**[1:k])/(2**[0:k-1])
			if(product(1-2*idx)==1) then
				pos(2)=pos(2)+product(merge(m,n,idx==0))
			endif
		enddo

		if(present(mgrp)) then
			if(allocated(mgrp)) deallocate(mgrp)
			allocate(mgrp(2,2))
			mgrp(:,1)=[1,-1]
			mgrp(1,2)=pos(2)-1
			mgrp(2,2)=prod(k+1)-mgrp(1,2)
		endif

		if(present(map)) then
			if(allocated(map)) then
				if(size(map)/=(prod(k+1)+1)) then
					deallocate(map)
				endif
			endif
			if(.not.allocated(map)) allocate(map(prod(k+1)))
			do i=1,prod(k+1)
				select case(product(sign(1,-mod(i-1,prod(2:))/prod(:k)/m)))
				case(1)
					map(pos(1))=i
					pos(1)=pos(1)+1
				case(-1)
					map(pos(2))=i
					pos(2)=pos(2)+1
				end select
			enddo
		endif

	end subroutine
	function contract(self,label) result(rt)
		class(t_tensor) :: self
		type(t_tensor), target :: rt
		character(*) :: label(:)
		integer :: i,j,n
		!real(8), pointer, contiguous :: T(:),MA(:,:)
		complex(8), pointer, contiguous :: T(:),MA(:,:)
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

		if(allocated(self%stp)) then
			call rt%new(self%label(size(label)+1:),self%shap(size(label)+1:),self%stp(size(label)+1:),flag="n")
		else
			call rt%new(self%label(size(label)+1:),self%shap(size(label)+1:),flag="n")
		endif
		rt%rc%T => T
		rt%is_return=.true.
		rt%is_conjg=self%is_conjg
		if(self%is_return) then
			call self%clear()
		endif
	end function
	function dot_tt(A,B,label) result(rt)
		class(t_tensor) :: A,B
		type(t_tensor), target :: rt
		character(*), optional :: label(:)
		character(max(len(A%label),len(B%label))) :: label_(size(A%label)+size(A%label))
		!real(8), pointer, contiguous :: T(:,:),MA(:,:),MB(:,:)
		!real(8), allocatable, target :: tmp(:)
		complex(8), pointer, contiguous :: T(:,:),MA(:,:),MB(:,:)
		complex(8), allocatable, target :: tmp(:)
		integer :: n,i,j,mrg(2),lrg(2),rrg(2)
		integer, allocatable :: orda(:),ordb(:)
		integer, allocatable :: lmap(:),rmap(:),lgrp(:,:),rgrp(:,:),mmap(:),mgrp(:,:)
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

		if(allocated(A%stp).and.allocated(B%stp)) then
		!if(.false.) then
			call smerge(gstp(A%stp(1:n)),mmap,mgrp)
			call smerge(gstp(A%stp(n+1:)),lmap,lgrp)
			call smerge(gstp(B%stp(n+1:)),rmap,rgrp)

			lrg=0
			rrg=0
			mrg=0
			T=0d0
			if(size(T)==1) then
				do i=1,size(mgrp,1)
					mrg=[min(mrg(2)+1,size(mmap)),min(mrg(2)+mgrp(i,2),size(mmap))]
					if(A%is_conjg.and.(.not.B%is_conjg)) then
						T=T+sum(conjg(MA(mmap(mrg(1):mrg(2)),1))*MB(mmap(mrg(1):mrg(2)),1))
					elseif(B%is_conjg.and.(.not.A%is_conjg)) then
						T=T+sum(MA(mmap(mrg(1):mrg(2)),1)*conjg(MB(mmap(mrg(1):mrg(2)),1)))
					else
						T=T+sum(MA(mmap(mrg(1):mrg(2)),1)*MB(mmap(mrg(1):mrg(2)),1))
					endif
				enddo
				if(A%is_conjg.and.B%is_conjg) then
					T=conjg(T)
				endif
			else
				do i=1,size(lgrp,1)
					lrg=[min(lrg(2)+1,size(lmap)),min(lrg(2)+lgrp(i,2),size(lmap))]
					rrg=[min(rrg(2)+1,size(rmap)),min(rrg(2)+rgrp(i,2),size(rmap))]
					mrg=[min(mrg(2)+1,size(mmap)),min(mrg(2)+mgrp(i,2),size(mmap))]
					if(A%is_conjg.and.(.not.B%is_conjg)) then
						call gemm(MA(mmap(mrg(1):mrg(2)),lmap(lrg(1):lrg(2))),MB(mmap(mrg(1):mrg(2)),rmap(rrg(1):rrg(2))),T(lrg(1):lrg(2),rrg(1):rrg(2)),transa="c")
					elseif(B%is_conjg.and.(.not.A%is_conjg)) then
						call gemm(MA(mmap(mrg(1):mrg(2)),lmap(lrg(1):lrg(2))),conjg(MB(mmap(mrg(1):mrg(2)),rmap(rrg(1):rrg(2)))),T(lrg(1):lrg(2),rrg(1):rrg(2)),transa="t")
					else
						call gemm(MA(mmap(mrg(1):mrg(2)),lmap(lrg(1):lrg(2))),MB(mmap(mrg(1):mrg(2)),rmap(rrg(1):rrg(2))),T(lrg(1):lrg(2),rrg(1):rrg(2)),transa="t")
					endif
				enddo
				lmap(lmap)=[1:size(lmap)]
				rmap(rmap)=[1:size(rmap)]
				if(A%is_conjg.and.B%is_conjg) then
					T=conjg(T(lmap,rmap))
				else
					T=T(lmap,rmap)
				endif
			endif
		else
			if(size(T)==1) then
				if(A%is_conjg.and.B%is_conjg) then
					T=conjg(sum(MA(:,1)*MB(:,1)))
				elseif(A%is_conjg) then
					T=sum(conjg(MA(:,1))*MB(:,1))
				elseif(B%is_conjg) then
					T=sum(MA(:,1)*conjg(MB(:,1)))
				else
					T=sum(MA(:,1)*MB(:,1))
				endif
			else
				if(A%is_conjg.and.B%is_conjg) then
					call gemm(MA,MB,T,transa="t")
					T=conjg(T)
				elseif(A%is_conjg) then
					call gemm(MA,MB,T,transa="c")
				elseif(B%is_conjg) then
					call gemm(MA,conjg(MB),T,transa="t")
				else
					call gemm(MA,MB,T,transa="t")
				endif
			endif
		endif


		if(allocated(A%stp).and.allocated(B%stp)) then
			call rt%new([character(5)::A%label(n+1:),B%label(n+1:)],[A%shap(n+1:),B%shap(n+1:)],[A%stp(n+1:),B%stp(n+1:)],flag="n")
		else
			call rt%new([character(5)::A%label(n+1:),B%label(n+1:)],[A%shap(n+1:),B%shap(n+1:)],flag="n")
		endif
		rt%rc%T(1:size(T)) => T
		rt%is_return=.true.

		if(A%is_return) then
			call A%clear()
		endif
		if(B%is_return) then
			call B%clear()
		endif
	end function
	function dot_tm(A,MB,label,stp) result(rt)
		class(t_tensor) :: A
		!real(8) :: MB(:,:)
		complex(8) :: MB(:,:)
		type(t_tensor), target :: rt
		character(*) :: label(:)
		integer, optional :: stp
		integer :: n
		!real(8), pointer, contiguous :: T(:,:),MA(:,:)
		complex(8), pointer, contiguous :: T(:,:),MA(:,:)
		integer, allocatable :: ord(:)
		n=size(A%shap)-size(label)
		call A%get_order([character::],label,ord)
		call A%reorder(ord)

		MA(1:product(A%shap(1:n)),1:product(A%shap(n+1:))) => A%rc%T

		allocate(T(size(MA,1),size(MB,2)))

		if(size(T)==1) then
			if(A%is_conjg) then
				T=sum(conjg(MA(:,1))*MB(:,1))
			else
				T=sum(MA(:,1)*MB(:,1))
			endif
		else
			if(A%is_conjg) then
				call gemm(conjg(MA),MB,T)
			else
				call gemm(MA,MB,T)
			endif
			!T=matmul(MA,MB)
		endif

		if(size(MB,1)/=size(MB,2).and.size(label)/=1) write(*,*)RAISEQQ(SIG$ABORT)

		if(allocated(A%stp)) then
			if(present(stp)) then
				call rt%new(A%label,[A%shap(:n),size(MB,2)],[A%stp(:n),stp],flag="n")
			else
				if(size(MB,1)==size(MB,2)) then
					call rt%new(A%label,A%shap,A%stp,flag="n")
				else
					write(*,*)RAISEQQ(SIG$ABORT)
				endif
			endif
		else
			if(size(MB,1)==size(MB,2)) then
				call rt%new(A%label,A%shap,flag="n")
			else
				call rt%new(A%label,[A%shap(:n),size(MB,2)],flag="n")
			endif
		endif
		rt%rc%T(1:size(T)) => T
		rt%is_return=.true.
		if(A%is_return) then
			call A%clear()
		endif
	end function
	function dot_mt(MA,B,label,stp) result(rt)
		!real(8) :: MA(:,:)
		complex(8) :: MA(:,:)
		class(t_tensor) :: B
		type(t_tensor), target :: rt
		character(*) :: label(:)
		integer, optional :: stp
		integer :: n
		!real(8), pointer, contiguous :: T(:,:),MB(:,:)
		complex(8), pointer, contiguous :: T(:,:),MB(:,:)
		integer, allocatable :: ord(:)
		n=size(label)
		call B%get_order(label,[character::],ord)
		call B%reorder(ord)

		MB(1:product(B%shap(1:n)),1:product(B%shap(n+1:))) => B%rc%T

		allocate(T(size(MA,1),size(MB,2)))

		if(size(T)==1) then
			if(B%is_conjg) then
				T=sum(MA(1,:)*conjg(MB(:,1)))
			else
				T=sum(MA(1,:)*MB(:,1))
			endif
		else
			if(B%is_conjg) then
				call gemm(MA,conjg(MB),T)
			else
				call gemm(MA,MB,T)
			endif
			!T=matmul(MA,MB)
		endif

		if(size(MA,1)/=size(MA,2).and.size(label)/=1) write(*,*)RAISEQQ(SIG$ABORT)

		if(allocated(B%stp)) then
			if(present(stp)) then
				call rt%new(B%label,[size(MA,1),B%shap(n+1:)],[stp,B%stp(n+1:)],flag="n")
			else
				if(size(MA,1)==size(MA,2)) then
					call rt%new(B%label,B%shap,B%stp,flag="n")
				else
					write(*,*)RAISEQQ(SIG$ABORT)
				endif
			endif
		else
			if(size(MA,1)==size(MA,2)) then
				call rt%new(B%label,B%shap,flag="n")
			else
				call rt%new(B%label,[size(MA,1),B%shap(n+1:)],flag="n")
			endif
		endif

		rt%rc%T(1:size(T)) => T
		rt%is_return=.true.
		if(B%is_return) then
			call B%clear()
		endif
	end function
	function dot_td(A,dB,label) result(rt)
		class(t_tensor) :: A
		real(8) :: dB(:)
		!complex(8) :: dB(:)
		type(t_tensor), target :: rt
		character(*) :: label(:)
		integer :: i,n
		!real(8), pointer, contiguous :: T(:,:),MA(:,:)
		complex(8), pointer, contiguous :: T(:,:),MA(:,:)
		integer, allocatable :: ord(:)
		n=size(label)
		call A%get_order(label,[character::],ord)
		call A%reorder(ord)

		MA(1:product(A%shap(1:n)),1:product(A%shap(n+1:))) => A%rc%T

		allocate(T(size(dB),size(MA,2)))

		if(A%is_conjg) then
			!$omp parallel do
			do i=1,size(T,2)
				T(:,i)=conjg(MA(:,i))*dB(:)
			enddo
			!$omp end parallel do
		else
			!$omp parallel do
			do i=1,size(T,2)
				T(:,i)=MA(:,i)*dB(:)
			enddo
			!$omp end parallel do
		endif

		if(allocated(A%stp)) then
			call rt%new(A%label,A%shap,A%stp,flag="n")
		else
			call rt%new(A%label,A%shap,flag="n")
		endif
		rt%rc%T(1:size(T)) => T
		rt%is_return=.true.
		if(A%is_return) then
			call A%clear()
		endif
	end function
	function dot_dt(dA,B,label) result(rt)
		real(8) :: dA(:)
		!complex(8) :: dA(:)
		class(t_tensor) :: B
		type(t_tensor), target :: rt
		character(*) :: label(:)
		integer :: i,n
		!real(8), pointer, contiguous :: T(:,:),MB(:,:)
		complex(8), pointer, contiguous :: T(:,:),MB(:,:)
		integer, allocatable :: ord(:)
		n=size(label)
		call B%get_order(label,[character::],ord)
		call B%reorder(ord)

		MB(1:product(B%shap(1:n)),1:product(B%shap(n+1:))) => B%rc%T

		allocate(T(size(dA),size(MB,2)))

		if(B%is_conjg) then
			!$omp parallel do
			do i=1,size(T,2)
				T(:,i)=dA(:)*conjg(MB(:,i))
			enddo
			!$omp end parallel do
		else
			!$omp parallel do
			do i=1,size(T,2)
				T(:,i)=dA(:)*MB(:,i)
			enddo
			!$omp end parallel do
		endif

		if(allocated(B%stp)) then
			call rt%new(B%label,B%shap,B%stp,flag="n")
		else
			call rt%new(B%label,B%shap,flag="n")
		endif
		rt%rc%T(1:size(T)) => T
		rt%is_return=.true.
		if(B%is_return) then
			call B%clear()
		endif

	end function
	!subroutine svd_t(self,l,r,U,s,V)
		!class(t_tensor) :: self
		!character(*) :: l(:),r(:)
		!!real(8), allocatable :: H(:,:)
		!real(8), allocatable, optional :: s(:)
		!!real(8), allocatable, optional :: U(:,:),V(:,:)
		!complex(8), allocatable :: H(:,:)
		!!complex(8), allocatable, optional :: s(:)
		!complex(8), allocatable, optional :: U(:,:),V(:,:)
		!call self%get_mat(l,r,H)
		!call svd(H,U,s,V)
	!end subroutine
	subroutine svd_t(self,l,r,U,s,V,cgrp,ctp)
		class(t_tensor) :: self
		character(*) :: l(:),r(:)
		!real(8), allocatable :: H(:,:)
		real(8), allocatable, optional :: s(:)
		!real(8), allocatable, optional :: U(:,:),V(:,:)
		complex(8), allocatable :: H(:,:)
		!complex(8), allocatable, optional :: s(:)
		complex(8), allocatable, optional :: U(:,:),V(:,:)
		integer, optional :: cgrp(:,:)
		character(*), optional :: ctp
		integer, allocatable :: lmap(:),rmap(:),lgrp(:,:),rgrp(:,:)
		integer :: i,mrg(2),lrg(2),rrg(2)
		logical :: is_return
		!if(any(isnan(abs(self%rc%T)))) then
			!write(*,*)RAISEQQ(SIG$ABORT)
		!endif
		is_return=self%is_return
		self%is_return=.false.
		call self%get_mat(l,r,H)
		if(size(l)/=0) then
			i=size(l)
		else
			i=size(self%shap)-size(r)
		endif

		if(allocated(self%stp)) then
			call smerge(gstp(self%stp(1:i)),lmap,lgrp)
			call smerge(gstp(self%stp(i+1:)),rmap,rgrp)
			
			H=H(lmap,rmap)

			mrg=0
			do i=1,size(lgrp,1)
				mrg(1)=mrg(1)+min(lgrp(i,2),rgrp(i,2))
			enddo
			if(present(cgrp)) then
				mrg(1)=max(sum(cgrp(:,2)),mrg(1))
			endif
			call allocate(s,mrg(1))
			call allocate(U,[size(H,1),size(s)])
			call allocate(V,[size(s),size(H,2)])
			s=0d0

			lrg=0
			rrg=0
			mrg=0
			do i=1,size(lgrp,1)
				lrg=[lrg(2)+1,lrg(2)+lgrp(i,2)]
				rrg=[rrg(2)+1,rrg(2)+rgrp(i,2)]
				mrg=[mrg(2)+1,mrg(2)+min(lgrp(i,2),rgrp(i,2))]
				U(:,mrg(1):mrg(2))=0d0
				V(mrg(1):mrg(2),:)=0d0
				if(present(ctp)) then
					call svd(H(lrg(1):lrg(2),rrg(1):rrg(2)),U(lrg(1):lrg(2),mrg(1):mrg(2)),s(mrg(1):mrg(2)),V(mrg(1):mrg(2),rrg(1):rrg(2)),ctp)
				else
					call svd(H(lrg(1):lrg(2),rrg(1):rrg(2)),U(lrg(1):lrg(2),mrg(1):mrg(2)),s(mrg(1):mrg(2)),V(mrg(1):mrg(2),rrg(1):rrg(2)))
				endif
				if(present(cgrp)) then
					mrg(1)=mrg(2)+1
					mrg(2)=sum(cgrp(1:i,2))
					U(:,mrg(1):mrg(2))=0d0
					V(mrg(1):mrg(2),:)=0d0
				endif
			enddo
			lmap(lmap)=[1:size(lmap)]
			rmap(rmap)=[1:size(rmap)]
			U=U(lmap,:)
			V=V(:,rmap)
		else
			if(present(ctp)) then
				call svd_a(H,U,s,V,ctp)
			else
				call svd_a(H,U,s,V)
			endif
		endif

		if(is_return) then
			call self%clear()
		endif
	end subroutine
	subroutine svd(H,U,s,V,ctp)
		!real(8) :: H(:,:)
		real(8), optional :: s(:)
		!real(8), allocatable, optional :: U(:,:),V(:,:)
		complex(8) :: H(:,:)
		complex(8), optional :: U(:,:),V(:,:)
		character(*), optional :: ctp
		if(.not.present(U)) then
			call gesdd(H,s)
		elseif(.not.present(s)) then
		else
			call gesdd(H,s,U,V,"S")
		endif
		if(present(ctp)) then
			select case(ctp)
			case("s")
				s=s(size(s):1:-1)
				if(present(U)) then
					U=U(:,size(s):1:-1)
					V=V(size(s):1:-1,:)
				endif
			case default
				write(*,*)RAISEQQ(SIG$ABORT)
			end select
		endif
	end subroutine
	subroutine svd_a(H,U,s,V,ctp)
		!real(8) :: H(:,:)
		real(8), allocatable, optional :: s(:)
		!real(8), allocatable, optional :: U(:,:),V(:,:)
		complex(8) :: H(:,:)
		complex(8), allocatable, optional :: U(:,:),V(:,:)
		character(*), optional :: ctp
		if(.not.present(U)) then
			call allocate(s,minval(shape(H)))
			if(present(ctp)) then
				call svd(H,s=s,ctp=ctp)
			else
				call svd(H,s=s)
			endif
		elseif(.not.present(s)) then
		else
			call allocate(s,minval(shape(H)))
			call allocate(U,[size(H,1),size(s)])
			call allocate(V,[size(s),size(H,2)])
			if(present(ctp)) then
				call svd(H,U,s,V,ctp)
			else
				call svd(H,U,s,V)
			endif
		endif
	end subroutine
	subroutine qr_2t(self,l,r,tQ,tR,label,stp)
		class(t_tensor) :: self
		character(*) :: l(:),r(:)
		type(t_tensor) :: tQ,tR
		character(*) :: label(:)
		integer, optional :: stp
		complex(8), allocatable :: MQ(:,:),MR(:,:)
		integer :: i
		logical :: is_return
		is_return=self%is_return
		self%is_return=.false.
		if(size(l)/=0) then
			i=size(l)
		else
			i=size(self%shap)-size(r)
		endif
		call self%qr_t(l,r,MQ,MR)
		if(present(stp)) then
			call tQ%new([self%label(1:i),label(1)],[self%shap(1:i),size(MQ,2)],[self%stp(1:i),stp])
			call tR%new([label(2),self%label(i+1:)],[size(MR,1),self%shap(i+1:)],[stp,self%stp(i+1:)])
		else
			call tQ%new([self%label(1:i),label(1)],[self%shap(1:i),size(MQ,2)])
			call tR%new([label(2),self%label(i+1:)],[size(MR,1),self%shap(i+1:)])
		endif
		call tQ%get_tensor(tQ%label,MQ)
		call tR%get_tensor(tR%label,MR)
		if(is_return) then
			call self%clear()
		endif
	end subroutine
	subroutine lq_2t(self,l,r,tL,tQ,label,stp)
		class(t_tensor) :: self
		character(*) :: l(:),r(:)
		type(t_tensor) :: tQ,tL
		character(*) :: label(:)
		integer, optional :: stp
		complex(8), allocatable :: MQ(:,:),ML(:,:)
		integer :: i
		logical :: is_return
		is_return=self%is_return
		self%is_return=.false.
		if(size(l)/=0) then
			i=size(l)
		else
			i=size(self%shap)-size(r)
		endif
		call self%lq_t(l,r,ML,MQ)
		if(present(stp)) then
			call tL%new([self%label(1:i),label(1)],[self%shap(1:i),size(ML,2)],[self%stp(1:i),stp])
			call tQ%new([label(2),self%label(i+1:)],[size(MQ,1),self%shap(i+1:)],[stp,self%stp(i+1:)])
		else
			call tL%new([self%label(1:i),label(1)],[self%shap(1:i),size(ML,2)])
			call tQ%new([label(2),self%label(i+1:)],[size(MQ,1),self%shap(i+1:)])
		endif
		call tQ%get_tensor(tQ%label,MQ)
		call tL%get_tensor(tL%label,ML)
		if(is_return) then
			call self%clear()
		endif
	end subroutine
	subroutine fact_t(self,l,r,tL,tR,tp,label,stp)
		class(t_tensor) :: self
		character(*) :: l(:),r(:)
		type(t_tensor) :: tL,tR
		character(*) :: label(:),tp
		integer, optional :: stp
		integer, allocatable :: mgrp(:,:)
		complex(8), pointer, contiguous :: ML(:,:),MR(:,:)
		integer :: i
		logical :: is_return
		is_return=self%is_return
		self%is_return=.false.
		nullify(ML,MR)
		if(size(l)/=0) then
			i=size(l)
		else
			i=size(self%shap)-size(r)
		endif
		if(present(stp)) then
			call self%fact_p(l,r,ML,MR,tp,mgrp)
			call tL%new([self%label(1:i),label(1)],[self%shap(1:i),size(ML,2)],[self%stp(1:i),stp],flag="n")
			call tR%new([label(2),self%label(i+1:)],[size(MR,1),self%shap(i+1:)],[stp,self%stp(i+1:)],flag="n")
			if(sum(abs(grp(:,:,stp)-mgrp))/=0) then
				if(sum(abs(grp(:,:,stp)))/=0) then
					write(*,*)"warning, overwrite grp ",stp," from ",grp(:,:,stp)," to ", mgrp
					write(*,*)RAISEQQ(SIG$ABORT)
				endif
				grp(:,:,stp)=mgrp
			endif
		else
			call self%fact_p(l,r,ML,MR,tp)
			call tL%new([self%label(1:i),label(1)],[self%shap(1:i),size(ML,2)],flag="n")
			call tR%new([label(2),self%label(i+1:)],[size(MR,1),self%shap(i+1:)],flag="n")
		endif
		tL%rc%T(1:size(ML)) => ML
		tR%rc%T(1:size(MR)) => MR
		if(is_return) then
			call self%clear()
		endif
	end subroutine
	subroutine fact_p(self,l,r,ML,MR,tp,mgrp)
		class(t_tensor) :: self
		character(*) :: l(:),r(:),tp
		integer, allocatable, optional :: mgrp(:,:)
		!real(8), allocatable :: H(:,:)
		!real(8), allocatable :: MQ(:,:),MR(:,:)
		complex(8), allocatable :: H(:,:)
		complex(8), contiguous, pointer :: ML(:,:),MR(:,:)
		integer, allocatable :: lmap(:),rmap(:),lgrp(:,:),rgrp(:,:),mgrp_(:,:)
		integer :: i,mrg(2),lrg(2),rrg(2)
		logical :: is_return
		is_return=self%is_return
		self%is_return=.false.
		call self%get_mat(l,r,H)
		if(size(l)/=0) then
			i=size(l)
		else
			i=size(self%shap)-size(r)
		endif
		if(allocated(self%stp)) then
			call smerge(gstp(self%stp(1:i)),lmap,lgrp)
			call smerge(gstp(self%stp(i+1:)),rmap,rgrp)

			allocate(mgrp_(size(lgrp,1),size(lgrp,2)))
			mgrp_(:,1)=lgrp(:,1)
			do i=1,size(lgrp,1)
				mgrp_(i,2)=min(lgrp(i,2),rgrp(i,2))
			enddo
			call allocate(ML,[size(H,1),sum(mgrp_(:,2))])
			call allocate(MR,[sum(mgrp_(:,2)),size(H,2)])

			!mrg=0
			!do i=1,size(lgrp,1)
				!mrg(1)=mrg(1)+min(lgrp(i,2),rgrp(i,2))
			!enddo
			!call allocate(ML,[size(H,1),mrg(1)])
			!call allocate(MR,[mrg(1),size(H,2)])

			H=H(lmap,rmap)
			ML=0d0
			MR=0d0

			lrg=0
			rrg=0
			mrg=0
			do i=1,size(lgrp,1)
				lrg=[lrg(2)+1,lrg(2)+lgrp(i,2)]
				rrg=[rrg(2)+1,rrg(2)+rgrp(i,2)]
				mrg=[mrg(2)+1,mrg(2)+min(lgrp(i,2),rgrp(i,2))]
				select case(tp)
				case("qr")
					call qr(H(lrg(1):lrg(2),rrg(1):rrg(2)),ML(lrg(1):lrg(2),mrg(1):mrg(2)),MR(mrg(1):mrg(2),rrg(1):rrg(2)))
				case("lq")
					call lq(H(lrg(1):lrg(2),rrg(1):rrg(2)),ML(lrg(1):lrg(2),mrg(1):mrg(2)),MR(mrg(1):mrg(2),rrg(1):rrg(2)))
				case("rq")
					call rq(H(lrg(1):lrg(2),rrg(1):rrg(2)),ML(lrg(1):lrg(2),mrg(1):mrg(2)),MR(mrg(1):mrg(2),rrg(1):rrg(2)))
				case("ql")
					call ql(H(lrg(1):lrg(2),rrg(1):rrg(2)),ML(lrg(1):lrg(2),mrg(1):mrg(2)),MR(mrg(1):mrg(2),rrg(1):rrg(2)))
				end select
			enddo
			lmap(lmap)=[1:size(lmap)]
			rmap(rmap)=[1:size(rmap)]

			ML=ML(lmap,:)
			MR=MR(:,rmap)
			if(present(mgrp)) then
				allocate(mgrp(size(mgrp_,1),size(mgrp_,2)))
				mgrp=mgrp_
			endif
		else
			call allocate(ML,[size(H,1),minval(shape(H))])
			call allocate(MR,[minval(shape(H)),size(H,2)])
			select case(tp)
			case("qr")
				call qr(H,ML,MR)
			case("lq")
				call lq(H,ML,MR)
			case("rq")
				call rq(H,ML,MR)
			case("ql")
				call ql(H,ML,MR)
			end select
		endif
		if(is_return) then
			call self%clear()
		endif
	end subroutine
	subroutine fact_a(self,l,r,ML,MR,tp,mgrp)
		class(t_tensor) :: self
		character(*) :: l(:),r(:),tp
		integer, allocatable, optional :: mgrp(:,:)
		!real(8), allocatable :: H(:,:)
		!real(8), allocatable :: MQ(:,:),MR(:,:)
		complex(8), allocatable :: H(:,:)
		complex(8), allocatable :: ML(:,:),MR(:,:)
		integer, allocatable :: lmap(:),rmap(:),lgrp(:,:),rgrp(:,:),mgrp_(:,:)
		integer :: i,mrg(2),lrg(2),rrg(2)
		logical :: is_return
		is_return=self%is_return
		self%is_return=.false.
		call self%get_mat(l,r,H)
		if(size(l)/=0) then
			i=size(l)
		else
			i=size(self%shap)-size(r)
		endif
		if(allocated(self%stp)) then
			call smerge(gstp(self%stp(1:i)),lmap,lgrp)
			call smerge(gstp(self%stp(i+1:)),rmap,rgrp)

			allocate(mgrp_(size(lgrp,1),size(lgrp,2)))
			mgrp_(:,1)=lgrp(:,1)
			do i=1,size(lgrp,1)
				mgrp_(i,2)=min(lgrp(i,2),rgrp(i,2))
			enddo
			call allocate(ML,[size(H,1),sum(mgrp_(:,2))])
			call allocate(MR,[sum(mgrp_(:,2)),size(H,2)])

			!mrg=0
			!do i=1,size(lgrp,1)
				!mrg(1)=mrg(1)+min(lgrp(i,2),rgrp(i,2))
			!enddo
			!call allocate(ML,[size(H,1),mrg(1)])
			!call allocate(MR,[mrg(1),size(H,2)])

			H=H(lmap,rmap)
			ML=0d0
			MR=0d0

			lrg=0
			rrg=0
			mrg=0
			do i=1,size(lgrp,1)
				lrg=[lrg(2)+1,lrg(2)+lgrp(i,2)]
				rrg=[rrg(2)+1,rrg(2)+rgrp(i,2)]
				!mrg=[mrg(2)+1,mrg(2)+min(lgrp(i,2),rgrp(i,2))]
				mrg=[mrg(2)+1,mrg(2)+mgrp_(i,2)]
				select case(tp)
				case("qr")
					call qr(H(lrg(1):lrg(2),rrg(1):rrg(2)),ML(lrg(1):lrg(2),mrg(1):mrg(2)),MR(mrg(1):mrg(2),rrg(1):rrg(2)))
				case("lq")
					call lq(H(lrg(1):lrg(2),rrg(1):rrg(2)),ML(lrg(1):lrg(2),mrg(1):mrg(2)),MR(mrg(1):mrg(2),rrg(1):rrg(2)))
				case("rq")
					call rq(H(lrg(1):lrg(2),rrg(1):rrg(2)),ML(lrg(1):lrg(2),mrg(1):mrg(2)),MR(mrg(1):mrg(2),rrg(1):rrg(2)))
				case("ql")
					call ql(H(lrg(1):lrg(2),rrg(1):rrg(2)),ML(lrg(1):lrg(2),mrg(1):mrg(2)),MR(mrg(1):mrg(2),rrg(1):rrg(2)))
				end select
			enddo
			lmap(lmap)=[1:size(lmap)]
			rmap(rmap)=[1:size(rmap)]

			ML=ML(lmap,:)
			MR=MR(:,rmap)
			if(present(mgrp)) then
				allocate(mgrp(size(mgrp_,1),size(mgrp_,2)))
				mgrp=mgrp_
			endif
		else
			select case(tp)
			case("qr")
				call qr_a(H,ML,MR)
			case("lq")
				call lq_a(H,ML,MR)
			case("rq")
				call rq_a(H,ML,MR)
			case("ql")
				call ql_a(H,ML,MR)
			end select
		endif
		if(is_return) then
			call self%clear()
		endif
	end subroutine
	subroutine qr_t(self,l,r,MQ,MR)
		class(t_tensor) :: self
		character(*) :: l(:),r(:)
		!real(8), allocatable :: H(:,:)
		!real(8), allocatable :: MQ(:,:),MR(:,:)
		complex(8), allocatable :: H(:,:)
		complex(8), allocatable :: MQ(:,:),MR(:,:)
		integer, allocatable :: lmap(:),rmap(:),lgrp(:,:),rgrp(:,:)
		integer :: i,mrg(2),lrg(2),rrg(2)
		logical :: is_return
		is_return=self%is_return
		self%is_return=.false.
		call self%get_mat(l,r,H)
		if(size(l)/=0) then
			i=size(l)
		else
			i=size(self%shap)-size(r)
		endif
		if(allocated(self%stp)) then
			call smerge(gstp(self%stp(1:i)),lmap,lgrp)
			call smerge(gstp(self%stp(i+1:)),rmap,rgrp)

			mrg=0
			do i=1,size(lgrp,1)
				mrg(1)=mrg(1)+min(lgrp(i,2),rgrp(i,2))
			enddo
			call allocate(MQ,[size(H,1),mrg(1)])
			call allocate(MR,[mrg(1),size(H,2)])
			H=H(lmap,rmap)
			MQ=0d0
			MR=0d0

			lrg=0
			rrg=0
			mrg=0
			do i=1,size(lgrp,1)
				lrg=[lrg(2)+1,lrg(2)+lgrp(i,2)]
				rrg=[rrg(2)+1,rrg(2)+rgrp(i,2)]
				mrg=[mrg(2)+1,mrg(2)+min(lgrp(i,2),rgrp(i,2))]
				call qr(H(lrg(1):lrg(2),rrg(1):rrg(2)),MQ(lrg(1):lrg(2),mrg(1):mrg(2)),MR(mrg(1):mrg(2),rrg(1):rrg(2)))
			enddo
			lmap(lmap)=[1:size(lmap)]
			rmap(rmap)=[1:size(rmap)]

			MQ=MQ(lmap,:)
			MR=MR(:,rmap)
		else
			call qr_a(H,MQ,MR)
		endif
		if(is_return) then
			call self%clear()
		endif
	end subroutine
	subroutine lq_t(self,l,r,ML,MQ)
		class(t_tensor) :: self
		character(*) :: l(:),r(:)
		!real(8), allocatable :: H(:,:)
		!real(8), allocatable :: MQ(:,:),ML(:,:)
		complex(8), allocatable :: H(:,:)
		complex(8), allocatable :: MQ(:,:),ML(:,:)
		integer, allocatable :: lmap(:),rmap(:),lgrp(:,:),rgrp(:,:)
		integer :: i,mrg(2),lrg(2),rrg(2)
		logical :: is_return
		is_return=self%is_return
		self%is_return=.false.
		call self%get_mat(l,r,H)
		if(size(l)/=0) then
			i=size(l)
		else
			i=size(self%shap)-size(r)
		endif
		if(allocated(self%stp)) then
			call smerge(gstp(self%stp(1:i)),lmap,lgrp)
			call smerge(gstp(self%stp(i+1:)),rmap,rgrp)

			mrg=0
			do i=1,size(lgrp,1)
				mrg(1)=mrg(1)+min(lgrp(i,2),rgrp(i,2))
			enddo
			call allocate(ML,[size(H,1),mrg(1)])
			call allocate(MQ,[mrg(1),size(H,2)])
			H=H(lmap,rmap)
			ML=0d0
			MQ=0d0

			lrg=0
			rrg=0
			mrg=0
			do i=1,size(lgrp,1)
				lrg=[lrg(2)+1,lrg(2)+lgrp(i,2)]
				rrg=[rrg(2)+1,rrg(2)+rgrp(i,2)]
				mrg=[mrg(2)+1,mrg(2)+min(lgrp(i,2),rgrp(i,2))]
				call lq(H(lrg(1):lrg(2),rrg(1):rrg(2)),ML(lrg(1):lrg(2),mrg(1):mrg(2)),MQ(mrg(1):mrg(2),rrg(1):rrg(2)))
			enddo
			lmap(lmap)=[1:size(lmap)]
			rmap(rmap)=[1:size(rmap)]

			ML=ML(lmap,:)
			MQ=MQ(:,rmap)
		else
			call lq_a(H,ML,MQ)
		endif
		if(is_return) then
			call self%clear()
		endif
	end subroutine
	subroutine qr(H,Q,R)
		!real(8) :: H(:,:)
		!real(8), allocatable :: Q(:,:),R(:,:)
		!real(8) :: tau(minval(shape(H)))
		complex(8) :: H(:,:)
		complex(8) :: Q(:,:),R(:,:)
		complex(8) :: tau(minval(shape(H)))
		integer :: i,j,m,n
		m=size(H,1)
		n=size(H,2)
		call geqrf(H,tau)
		do i=1,m
			do j=1,n
				if(i>j) then
					if(i<=size(R,1)) R(i,j)=0d0
					Q(i,j)=H(i,j)
				else
					R(i,j)=H(i,j)
				endif
			enddo
		enddo
		!call orgqr(Q,tau)
		call ungqr(Q,tau)
	end subroutine
	subroutine lq(H,L,Q)
		!real(8) :: H(:,:)
		!real(8), allocatable :: Q(:,:),L(:,:)
		!real(8) :: tau(minval(shape(H)))
		complex(8) :: H(:,:)
		complex(8) :: Q(:,:),L(:,:)
		complex(8) :: tau(minval(shape(H)))
		integer :: i,j,m,n
		m=size(H,1)
		n=size(H,2)
		call gelqf(H,tau)
		do i=1,m
			do j=1,n
				if(i<j) then
					if(j<=size(L,2)) L(i,j)=0d0
					Q(i,j)=H(i,j)
				else
					L(i,j)=H(i,j)
				endif
			enddo
		enddo
		!call orglq(Q,tau)
		call unglq(Q,tau)
	end subroutine
	subroutine rq(H,R,Q)
		!real(8) :: H(:,:)
		!real(8), allocatable :: Q(:,:),R(:,:)
		!real(8) :: tau(minval(shape(H)))
		complex(8) :: H(:,:)
		complex(8) :: Q(:,:),R(:,:)
		complex(8) :: tau(minval(shape(H)))
		integer :: i,j,m,n,d
		m=size(H,1)
		n=size(H,2)
		call gerqf(H,tau)
		if(m<=n) then
			d=n-m
			do i=1,m
				do j=1,n
					if(i<=(j-d)) then
						R(i,j-d)=H(i,j)
					else
						if(j>d) R(i,j-d)=0d0
						Q(i,j)=H(i,j)
					endif
				enddo
			enddo
		else
			d=m-n
			do i=1,m
				do j=1,n
					if((i-d)<=j) then
						R(i,j)=H(i,j)
					else
						R(i,j)=0d0
						if(i>d) Q(i-d,j)=H(i,j)
					endif
				enddo
			enddo
		endif
		!call orgqr(Q,tau)
		call ungrq(Q,tau)
	end subroutine
	subroutine ql(H,Q,L)
		!real(8) :: H(:,:)
		!real(8), allocatable :: Q(:,:),R(:,:)
		!real(8) :: tau(minval(shape(H)))
		complex(8) :: H(:,:)
		complex(8) :: Q(:,:),L(:,:)
		complex(8) :: tau(minval(shape(H)))
		integer :: i,j,m,n,d
		m=size(H,1)
		n=size(H,2)
		call geqlf(H,tau)
		if(m<=n) then
			d=n-m
			do i=1,m
				do j=1,n
					if(i>=(j-d)) then
						L(i,j)=H(i,j)
					else
						L(i,j)=0d0
						if(j>d) Q(i,j-d)=H(i,j)
					endif
				enddo
			enddo
		else
			d=m-n
			do i=1,m
				do j=1,n
					if((i-d)>=j) then
						L(i-d,j)=H(i,j)
					else
						if(i>d) L(i-d,j)=0d0
						Q(i,j)=H(i,j)
					endif
				enddo
			enddo
		endif
		!call orgqr(Q,tau)
		call ungql(Q,tau)
	end subroutine
	subroutine qr_a(H,Q,R)
		!real(8) :: H(:,:)
		!real(8), allocatable :: Q(:,:),R(:,:)
		!real(8) :: tau(minval(shape(H)))
		complex(8) :: H(:,:)
		complex(8), allocatable :: Q(:,:),R(:,:)
		call allocate(Q,[size(H,1),minval(shape(H))])
		call allocate(R,[minval(shape(H)),size(H,2)])
		call qr(H,Q,R)
	end subroutine
	subroutine lq_a(H,L,Q)
		!real(8) :: H(:,:)
		!real(8), allocatable :: Q(:,:),L(:,:)
		complex(8) :: H(:,:)
		complex(8), allocatable :: Q(:,:),L(:,:)
		call allocate(Q,[minval(shape(H)),size(H,2)])
		call allocate(L,[size(H,1),minval(shape(H))])
		call lq(H,L,Q)
	end subroutine
	subroutine rq_a(H,R,Q)
		!real(8) :: H(:,:)
		!real(8), allocatable :: Q(:,:),R(:,:)
		complex(8) :: H(:,:)
		complex(8), allocatable :: Q(:,:),R(:,:)
		call allocate(R,[size(H,1),minval(shape(H))])
		call allocate(Q,[minval(shape(H)),size(H,2)])
		call rq(H,R,Q)
	end subroutine
	subroutine ql_a(H,Q,L)
		!real(8) :: H(:,:)
		!real(8), allocatable :: Q(:,:),R(:,:)
		complex(8) :: H(:,:)
		complex(8), allocatable :: Q(:,:),L(:,:)
		call allocate(Q,[size(H,1),minval(shape(H))])
		call allocate(L,[minval(shape(H)),size(H,2)])
		call ql(H,Q,L)
	end subroutine
	subroutine smat_inv(lstp,rstp,M)
		integer :: lstp(:),rstp(:)
		complex(8) :: M(:,:)
		complex(8), allocatable :: Mp(:,:)
		integer, allocatable :: lmap(:),rmap(:),lgrp(:,:),rgrp(:,:)
		integer :: i,lrg(2),rrg(2)
		call smerge(gstp(lstp),lmap,lgrp)
		call smerge(gstp(rstp),rmap,rgrp)
		call allocate(Mp,M)
		M=M(lmap,rmap)
		lrg=0
		rrg=0
		do i=1,size(lgrp,1)
			lrg=[lrg(2)+1,lrg(2)+lgrp(i,2)]
			rrg=[rrg(2)+1,rrg(2)+rgrp(i,2)]
			call mat_inv(M(lrg(1):lrg(2),rrg(1):rrg(2)))
		enddo
		lmap(lmap)=[1:size(lmap)]
		rmap(rmap)=[1:size(rmap)]
		M=M(rmap,lmap)
	end subroutine
	subroutine smat_eg(stp,M,Eg,cM,cgrp)
		integer :: stp(:)
		complex(8) :: M(:,:)
		real(8), allocatable :: Eg(:)
		complex(8), allocatable, optional :: cM(:,:)
		integer, optional :: cgrp(:,:)
		integer, allocatable :: map(:),grp(:,:)
		integer :: i,rg(2),crg(2)
		call smerge(gstp(stp),map,grp)
		if(present(cgrp)) then
			call allocate(Eg,max(size(M,2),sum(cgrp(:,2))))
			call allocate(cM,[size(M,1),max(size(M,2),sum(cgrp(:,2)))])
			cM=0d0
			Eg=0d0
			M=-M(map,map)
		else
			call allocate(Eg,size(M,1))
			M=M(map,map)
		endif
		rg=0
		crg=0
		do i=1,size(grp,1)
			rg=[rg(2)+1,rg(2)+grp(i,2)]
			crg=[crg(2)+1,crg(2)+grp(i,2)]
			call heev(M(rg(1):rg(2),rg(1):rg(2)),Eg(crg(1):crg(2)),"V")
			if(present(cgrp)) then
				rg(2)=rg(1)-1+min(cgrp(i,2),grp(i,2))
				crg(2)=crg(1)-1+min(cgrp(i,2),grp(i,2))
				cM(rg(1):rg(1)-1+grp(i,2),crg(1):crg(2))=M(rg(1):rg(1)-1+grp(i,2),rg(1):rg(2))
				rg(2)=rg(1)-1+grp(i,2)
				crg(2)=crg(1)-1+cgrp(i,2)
			endif
		enddo
		map(map)=[1:size(map)]
		if(present(cgrp)) then
			cM=cM(map,:)
			Eg=-Eg
		else
			M=M(map,:)
		endif
	end subroutine
	subroutine remove(self,t)
		class(t_tensor), target :: self
		type(t_tensor), target :: t
		type(t_tensor), pointer :: node
		integer :: i
		self%rc%c=self%rc%c-1
		if(self%rc%c==0) then
			deallocate(self%rc%T)
			deallocate(self%rc)
		endif
		do i=1,2
			node => self
			do
				if(associated(node,t)) then
					if(associated(node%link(1)%tg).and.associated(node%link(2)%tg)) then
						node%link(1)%tg%link(2)%tg => node%link(2)%tg
						node%link(2)%tg%link(1)%tg => node%link(1)%tg
					elseif(associated(node%link(1)%tg)) then
						nullify(node%link(1)%tg%link(2)%tg)
					elseif(associated(node%link(2)%tg)) then
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
		if(associated(self%link(2)%tg)) then
			t%link(2)%tg => self%link(2)%tg
			t%link(2)%tg%link(1)%tg => t
		endif
		self%link(2)%tg => t
		t%link(1)%tg => self
	end subroutine
	subroutine clear(self)
		class(t_tensor) :: self
		if(allocated(self%shap)) deallocate(self%label,self%shap)
		if(allocated(self%stp)) deallocate(self%stp)
		if(associated(self%rc)) then
			call self%remove(self)
		endif
		self%is_return=.false.
		self%is_conjg=.false.
	end subroutine
	subroutine delete(self)
		type(t_tensor) :: self
		call self%clear()
	end subroutine
	subroutine delete_all(self)
		type(t_tensor) :: self(:)
		integer :: i
		do i=1,size(self)
			call self(i)%clear()
		enddo
	end subroutine
	subroutine delete_all2(self)
		type(t_tensor) :: self(:,:)
		integer :: i,j
		do i=1,size(self,1)
			do j=1,size(self,2)
				call self(i,j)%clear()
			enddo
		enddo
	end subroutine
	subroutine ar_naupd_r(H,E,V,lr)
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
	subroutine ar_naupd_c(H,E,V,lr)
		integer, parameter :: ncv = 20, nev = 1
		complex(8) :: H(:,:),V(:)
		complex(8) :: E
		real(8) :: tol
		integer :: iparam(11), ipntr(14), info, ierr, ido ,n
		complex(8) :: d(ncv,3), resid(size(H,1)), vec(size(H,1),ncv), workd(3*size(H,1)), workev(3*ncv), workl(3*ncv*ncv+6*ncv),sigma
		real(8) :: rwork(ncv)
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
			call znaupd ( ido, 'I', n, 'LM', nev, tol, resid, ncv, vec, n, iparam, ipntr, workd, workl, size(workl), rwork, info )
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
					call zneupd ( .true., 'A', select, d, vec, n, sigma, workev, 'I', n, 'LM', nev, tol, resid, ncv, vec, n, iparam, ipntr, workd, workl, size(workl), rwork, ierr )
					exit
				else
					stop "arpack info error"
				endif
			endif
		enddo
		E=d(1,1)
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
	subroutine allocate_1c(A,shap)
		complex(8), allocatable :: A(:)
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
	subroutine allocate_1ccopy(A,B)
		complex(8), allocatable :: A(:)
		complex(8) :: B(:)
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
	subroutine allocate_2c(A,shap)
		complex(8), allocatable :: A(:,:)
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
	subroutine allocate_2ccopy(A,B)
		complex(8), allocatable :: A(:,:)
		complex(8) :: B(:,:)
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
	subroutine allocate_1cp(A,shap)
		complex(8), pointer :: A(:)
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
	subroutine allocate_1cpcopy(A,B)
		complex(8), pointer :: A(:)
		complex(8) :: B(:)
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
	subroutine allocate_2cp(A,shap)
		complex(8), pointer :: A(:,:)
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
	subroutine allocate_2cpcopy(A,B)
		complex(8), pointer :: A(:,:)
		complex(8) :: B(:,:)
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
end module
