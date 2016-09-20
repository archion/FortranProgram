module M_tensor
	use lapack95
	use M_utility
	implicit none
	type t_tensor
		character(:), allocatable :: label
		integer, allocatable :: l(:)
		integer, allocatable :: ord(:)
		integer, allocatable :: shap(:) 
		real(8), allocatable :: T(:) 
		integer, allocatable :: prod(:)
		integer :: clock=1
	contains
		procedure :: new
		procedure :: get_idx
		procedure :: reorder
		procedure :: get_mat
		procedure :: get_tensor
		procedure :: dot
		procedure :: clear
		procedure :: find_label
		procedure :: get_label
		procedure :: change_label
		procedure :: svd
		procedure :: show
	end type
contains
	subroutine show(self)
		class(t_tensor) :: self
		integer :: ord(size(self%shap))
		ord(self%ord)=[1:size(self%shap)]
		write(*,*)self%label
		write(*,*)self%l
		write(*,*)self%shap(ord)
	end subroutine
	subroutine new(self,shap,label,flag)
		class(t_tensor) :: self
		integer :: shap(:)
		character(*) :: label
		character, optional :: flag
		integer :: i,j
		self%clock=1
		if(product(self%shap)/=product(shap)) then
			if(allocated(self%T)) deallocate(self%T)
			allocate(self%T(product(shap)))
		endif
		self%label=label
		self%shap=shap
		allocate(self%prod,mold=shap)
		allocate(self%l,mold=shap)
		if(size(self%shap)==0) return
		self%prod(1)=1
		self%l=0
		j=1
		self%ord=[1:size(shap)]
		do i=1,len(label)
			if(label(i:i)/=" ") then
				self%l(j)=self%l(j)+1
			else
				j=j+1
			endif
		enddo
		do i=2,size(shap)
			self%prod(i)=self%prod(i-1)*self%shap(i-1)
		enddo
		if(present(flag)) then
			select case(flag)
			case("0")
				self%T=0d0
			case("r")
				call random_number(self%T)
			end select
		endif
	end subroutine
	function get_idx(self,idx)
		class(t_tensor) :: self
		integer :: idx(:)
		integer :: get_idx(size(idx)/size(self%shap))
		integer :: i,j
		do j=1,size(get_idx)
			get_idx(j)=1
			do i=1,size(self%shap)
				get_idx(j)=get_idx(j)+(idx((j-1)*size(self%shap)+i)-1)*self%prod(self%ord(i))
			enddo
		enddo
	end function
	!subroutine reorder(self,ord)
		!class(t_tensor) :: self
		!integer :: ord(:)
		!integer :: i,map(size(self%T)),prod(size(self%prod))
		!self%ord=self%ord(ord)

		!prod(1)=1
		!do i=2,size(self%shap)
			!prod(i)=prod(i-1)*self%shap(self%ord(i-1))
		!enddo

		!self%shap=self%shap(self%ord)
		!self%prod=self%prod(self%ord)

		!do i=1,size(self%T)
			!map(i)=1+sum((mod((i-1)/prod,self%shap))*self%prod)
		!enddo

		!self%T=self%T(map)
		!self%prod=prod

		!self%ord(ord)=[1:size(ord)]
	!end subroutine
	subroutine reorder(self,ord)
		class(t_tensor) :: self
		integer :: ord(:)
		integer :: i,map(size(self%T)),prod(size(self%prod))
		self%ord(self%ord)=[1:size(ord)]
		!self%ord=ord(self%ord)
		self%ord=self%ord(ord)

		prod(1)=1
		do i=2,size(self%shap)
			prod(i)=prod(i-1)*self%shap(self%ord(i-1))
		enddo

		self%shap=self%shap(self%ord)
		self%prod=self%prod(self%ord)

		do i=1,size(self%T)
			map(i)=1+sum((mod((i-1)/prod,self%shap))*self%prod)
		enddo

		self%T=self%T(map)
		self%prod=prod

		self%ord=ord
	end subroutine
	!subroutine leg_swap(self,label)
	!end subroutine
	!subroutine trace(self,label)
	!end subroutine
	subroutine change_label(self,label)
		class(t_tensor) :: self
		character(*) :: label(:)
		integer :: i,j
		do i=1,size(self%l)
			do j=1,size(label),2
				if(label(j)==self%label((sum(self%l(:i))+i-1)-self%l(i)+1:sum(self%l(:i))+i-1)) then
					self%label=self%label(:(sum(self%l(:i))+i-1)-self%l(i))//label(j+1)//self%label(sum(self%l(:i))+i:)
					self%l(i)=len(label(j+1))
					exit
				endif
			enddo
		enddo
	end subroutine
	subroutine get_mat(self,label,n,H)
		class(t_tensor) :: self
		character(*) :: label
		integer :: n
		real(8), allocatable :: H(:,:)
		integer :: ia(size(self%shap)),shap(2),i,l(2),m(2)
		if(n<0) then
			m=[size(self%shap)-abs(n),abs(n)]
			l(2)=self%find_label(label)
			l(1)=mod(l(2)-self%clock-1+sum(m),sum(m))+1
		else
			m=[abs(n),size(self%shap)-abs(n)]
			l(1)=self%find_label(label)
			l(2)=mod(l(1)+self%clock-1+sum(m),sum(m))+1
		endif
		ia(:m(1))=mod(l(1)-[0:m(1)-1]*self%clock-1+sum(m),sum(m))+1
		ia(m(1)+1:)=mod(l(2)+[0:m(2)-1]*self%clock-1+sum(m),sum(m))+1
		call self%reorder(ia)
		shap=[product(self%shap(:m(1))),product(self%shap(m(1)+1:))]
		if(any((shape(H)-shap)/=0)) then
			if(allocated(H)) deallocate(H)
			allocate(H(shap(1),shap(2)))
		endif
		H=reshape(self%T,shape(H))
	end subroutine
	!subroutine get_mat(self,label,n,H)
		!class(t_tensor) :: self
		!character(*) :: label
		!integer :: n
		!real(8), allocatable :: H(:,:)
		!integer :: ia(size(self%shap)),shap(2),i
		!i=self%find_label(label)
		!if(n>0) then
			!ia=mod(i+[0:size(ia)-1]-1,size(self%shap))+1
		!else
			!ia=mod(i+[[0:n+1:-1]+size(self%shap),[1:size(ia)-n]]-1,size(self%shap))+1
		!endif
		!call self%reorder(ia)
		!shap=[product(self%shap(:abs(n))),product(self%shap(abs(n)+1:))]
		!if(any((shape(H)-shap)/=0)) then
			!if(allocated(H)) deallocate(H)
			!allocate(H(shap(1),shap(2)))
		!endif
		!H=reshape(self%T,shape(H))
	!end subroutine
	!real(8) function trace(A,B,label)
		!class(t_tensor) :: A,B
		!character(*) :: label(:)
		!character(:), allocatable :: nl
		!integer :: n
		!integer :: i,j
		!integer :: ia(size(B%shap))
		!i=self%find_label(label(1))
		!call self%reorder(ia)
		!ib=mod(i+[0:size(B%shap)-1:-1]+size(B%shap)-1,size(B%shap))+1
	!end function
	subroutine get_tensor(self,label,n,H)
		class(t_tensor) :: self
		character(*) :: label
		integer :: n
		real(8) :: H(:,:)
		integer :: ia(size(self%shap)),i,l(2),m(2)
		if(n<0) then
			m=[size(self%shap)-abs(n),abs(n)]
			l(2)=self%find_label(label)
			l(1)=mod(l(2)-self%clock-1+sum(m),sum(m))+1
		else
			m=[abs(n),size(self%shap)-abs(n)]
			l(1)=self%find_label(label)
			l(2)=mod(l(1)+self%clock-1+sum(m),sum(m))+1
		endif
		ia(:m(1))=mod(l(1)-[0:m(1)-1]*self%clock-1+sum(m),sum(m))+1
		ia(m(1)+1:)=mod(l(2)+[0:m(2)-1]*self%clock-1+sum(m),sum(m))+1
		call self%reorder(ia)
		self%T=reshape(H,shape(self%T))
	end subroutine
	subroutine dot(self,A,B,label,n)
		class(t_tensor) :: self,A,B
		character(*) :: label(:)
		character(:), allocatable :: nl
		integer :: n
		integer :: i,j
		real(8), allocatable :: T(:,:),MA(:,:),MB(:,:)
		!i=mod(A%find_label(label(1))-1+n,size(A%shap))+1
		!j=sum(A%l(:i))+i-1
		!call A%get_mat(A%label(j-A%l(A%ord(i))+1:j),size(A%shap)-n,MA)
		call A%get_mat(label(1),-n,MA)
		call B%get_mat(label(2),n,MB)

		allocate(T(size(MA,1),size(MB,2)))

		if(size(T)==1) then
			T=sum(MA(1,:)*MB(:,1))
		else
			T=matmul(MA,MB)
		endif

		nl=''

		do i=size(A%shap)-n,1,-1
			nl=nl//A%get_label(A%ord(i))//" "
		enddo
		do i=n+1,size(B%shap)
			nl=nl//B%get_label(B%ord(i))
			if(i/=size(B%shap)) then
				nl=nl//" "
			endif
		enddo
		i=size(A%shap)-n
		call self%new([A%shap(:i),B%shap(n+1:)],nl)
		self%ord(:i)=[i:1:-1]
		self%T=reshape(T,shape(self%T))
	end subroutine
	function get_label(self,i)
		class(t_tensor) :: self
		integer :: i
		character(:), allocatable :: get_label
		integer :: j
		j=sum(self%l(:i))+i-1
		get_label=self%label(j-self%l(i)+1:j)
	end function
	integer function find_label(self,label)
		class(t_tensor) :: self
		character(*) :: label
		integer :: i
		do i=1,size(self%l)
			if(label==self%label((sum(self%l(:i))+i-1)-self%l(i)+1:sum(self%l(:i))+i-1)) then
				find_label=i
				return
			elseif(i==size(self%l)) then
				write(*,*)"error, can't find label ",label," in ",self%label
				stop 
			endif
		enddo
		find_label=0
	end function
	subroutine svd(self,label,n,U,s,V)
		class(t_tensor) :: self
		character(*) :: label
		integer :: n
		real(8), allocatable :: H(:,:)
		real(8), allocatable, optional :: s(:)
		real(8), allocatable, optional :: U(:,:),V(:,:)
		call self%get_mat(label,n,H)
		if(.not.present(U)) then
			if(size(s)/=minval(shape(H))) then
				if(allocated(s)) deallocate(s)
				allocate(s(minval(shape(H))))
			endif
			call gesdd(H,s)
		elseif(.not.present(s)) then
		else
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
		endif
	end subroutine
	subroutine clear(self)
		class(t_tensor) :: self
		if(allocated(self%shap)) deallocate(self%shap,self%T,self%prod,self%l,self%ord)
	end subroutine
end module
