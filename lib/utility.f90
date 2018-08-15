module M_utility
	USE IFPORT, ifort_qsort => qsort
	use omp_lib
	use mkl_service
	use M_rd
	implicit none
	type t_sort
		real(8) :: val
	contains
		procedure :: swap_sort
	end type
	interface integrate
		module procedure rintegrate,cintegrate
	end interface
	interface get
		module procedure cget,rget,iget
	end interface
	interface find
		module procedure rfind
	end interface
	interface mwrite
		module procedure cmwrite,rmwrite,imwrite
	end interface
	interface arth
		module procedure carth,rarth,iarth
	end interface
	interface swap
		module procedure siswap,srswap,vcswap,viswap,vrswap
	end interface
	interface to_char
		module procedure to_char_single,to_char_array
	end interface
	interface qsort
		module procedure rqsort,sqsort,rmqsort
	end interface
	interface collect
		module procedure rcollect,scollect,rmcollect
	end interface
	real(8) :: otime(0:6)=0d0
	integer, parameter :: reset=0,add_log_linear=1,add_linear=2
contains
	logical function next(x,y,dx,tol) result(rt)
		real(8) :: x,y
		real(8), optional :: dx,tol
		real(8), save :: px,py,dx_,tol_
		logical, save :: flag
		rt=.false.
		if(present(dx)) then 
			py=transfer([Z'00000000',Z'7FF80000'],1d0)
			flag=.false.
			dx_=dx
			tol_=tol
			return
		endif
		if(.not.isnan(py)) then
			if(py*y<0d0) flag=.true.
			if(flag) then
				dx_=dx_/2d0
			endif
			if(py*y<0d0) then
				dx_=-dx_
			elseif(sign(1d0,py)*sign(1d0,(py-y)/(px-x)*dx_)>0d0) then 
				x=px
				dx_=-dx_
			endif
		endif
		px=x
		py=y
		if(abs(dx_)<tol_) then
			rt=.true.
		else
			x=x+dx_
		endif
	end function
	real(8) function rintegrate(x,A) result(f)
		real(8) :: x(:)
		real(8) :: A(:)
		integer :: i
		f=0d0
		do i=1,size(x)-1
			f=f+0.5d0*(x(i+1)-x(i))*(A(i+1)+A(i))
		enddo
	end function
	complex(8) function cintegrate(x,A) result(f)
		real(8) :: x(:)
		complex(8) :: A(:)
		integer :: i
		f=0d0
		do i=1,size(x)-1
			f=f+0.5d0*(x(i+1)-x(i))*(A(i+1)+A(i))
		enddo
	end function
	pure real(8) function ff(x)
		real(8), intent(in) :: x
		ff=1d0/(exp(x)+1d0)
	end function
	pure real(8) function fb(x)
		real(8), intent(in) :: x
		integer :: m
		real(8) :: f1,f2
		if(abs(x)<1d-20) then
			fb=-0.5d0
		elseif(x<0d0) then
			fb=1d0/(exp(x)-1d0)
		else
			fb=-exp(-x)/(exp(-x)-1d0)
		endif
	end function
	subroutine set_grid(grid,flag,from,to,n)
		real(8) :: grid(:),from(:),to(:)
		integer :: flag(:),n(:)
		integer,save :: n_
		integer :: i,j,k,l,u,m,o
		real(8) :: base,from_,dx
		m=1
		o=1
		do k=1,size(flag)
			select case(flag(k))
			case(reset)
				n_=0
			case(add_log_linear)
				l=find(grid(1:n_),from(o))
				u=find(grid(1:n_),to(o))
				grid(u+n(m)*n(m+1)+1:n_+n(m)*n(m+1))=grid(u+1:n_)
				u=u+n(m)*n(m+1)
				n_=n_+n(m)*n(m+1)
				base=(to(o)/from(o))**(1d0/n(m))
				from_=from(o)/base
				do i=l+1,u
					if(mod(i-l-1,n(m+1))==0) then
						from_=from_*base
						dx=from_*(base-1d0)/n(m+1)
					endif
					grid(i)=from_+dx*(mod(i-l-1,n(m+1))+1)
				enddo
				m=m+2
				o=o+1
			case(add_linear)
				l=find(grid(1:n_),from(o))
				u=find(grid(1:n_),to(o))
				grid(u+n(m)+1:n_+n(m))=grid(u+1:n_)
				u=u+n(m)
				n_=n_+n(m)
				dx=(to(o)-from(o))/(u-l)
				do i=l+1,u
					grid(i)=from(o)+dx*(i-l)
				enddo
				m=m+1
				o=o+1
			end select
		enddo
	end subroutine
	pure integer function rfind(A,x) result(rt)
		real(8), intent(in) :: A(:)
		real(8), intent(in) :: x
		integer :: m,u
		rt=0
		u=size(A)+1
		do
			if((u-rt)>1) then
				m=(rt+u)/2
				if(x>=A(m)) then
					rt=m
				else
					u=m
				endif
			else
				exit
			endif
		enddo
	end function
	pure complex(8) function cget(x,i) result(rt)
		complex(8), intent(in) :: x(:)
		integer, intent(in) :: i
		rt=x(i)
	end function
	pure real(8) function rget(x,i) result(rt)
		real(8), intent(in) :: x(:)
		integer, intent(in) :: i
		rt=x(i)
	end function
	pure integer function iget(x,i) result(rt)
		integer, intent(in) :: x(:)
		integer, intent(in) :: i
		rt=x(i)
	end function
	subroutine show_time()
		integer :: i
		write(*,"(A)")"The timing results is:"
		do i=1,ubound(otime,1)
			if(otime(i)>1d-6) then
				write(*,"(i4,'m',f6.3,'s')")int(otime(i)/60d0),otime(i)-int(otime(i)/60d0)*60d0
			endif
		enddo
	end subroutine
	subroutine cmwrite(f,A)
		complex(8) :: A(:,:)
		integer :: f,i
		do i=1,size(A,1)
			!write(f,"(sp,e16.4,e11.4'I'$)")A(i,:)
			!write(f,"(sp,e10.2','e9.2$)")A(i,:)
			write(f,"(es12.4$)")real(A(i,:))
			write(f,"(1X)")
		enddo
		do i=1,size(A,1)
			!write(f,"(sp,e16.4,e11.4'I'$)")A(i,:)
			!write(f,"(sp,e10.2','e9.2$)")A(i,:)
			write(f,"(es12.4$)")imag(A(i,:))
			write(f,"(1X)")
		enddo
		write(f,"(1X)")
	end subroutine
	subroutine rmwrite(f,A)
		real(8) :: A(:,:)
		integer :: f,i
		do i=1,size(A,1)
			write(f,"(es12.4$)")A(i,:)
			write(f,"(1X)")
		enddo
		write(f,"(1X)")
	end subroutine
	subroutine imwrite(f,A)
		integer :: f,i,A(:,:)
		do i=1,size(A,1)
			write(f,"(i4$)")A(i,:)
			write(f,"(1X)")
		enddo
		write(f,"(1X)")
	end subroutine
	subroutine siswap(a,b)
		integer :: a,b,tmp
		tmp=a
		a=b
		b=tmp
	end subroutine
	subroutine srswap(a,b)
		real(8) :: a,b,tmp
		tmp=a
		a=b
		b=tmp
	end subroutine
	subroutine vrswap(a,b)
		real(8) :: a(:),b(:),tmp(size(a))
		tmp=a
		a=b
		b=tmp
	end subroutine
	subroutine vcswap(a,b)
		complex(8) :: a(:),b(:),tmp(size(a))
		tmp=a
		a=b
		b=tmp
	end subroutine
	subroutine viswap(a,b)
		integer :: a(:),b(:),tmp(size(a))
		tmp=a
		a=b
		b=tmp
	end subroutine
	subroutine find_peak(a,ipeak)
		real(8) :: a(:)
		integer, allocatable :: ipeak(:)
		integer :: ipeak_(10),n,i
		n=0
		do i=2,size(a)-1
			if((a(i)-a(i-1))>0d0.and.(a(i)-a(i+1))>0d0) then
				n=n+1
				if(n>10) then
					n=n-1
					exit
				endif
				ipeak_(n)=i
			endif
		enddo
		allocate(ipeak(n))
		ipeak=ipeak_(:n)
	end subroutine
	subroutine is_cross(pa,a,sg)
		real(8) :: pa,a
		integer :: sg
		sg=0
		if(pa/=0d0) then
			if(pa>0d0.and.a<0d0) then
				sg=1
			endif
			if(pa<0d0.and.a>0d0) then
				sg=-1
			endif
		endif
		pa=a
	end subroutine
	function carth(f,d,n)
		integer :: n,i
		complex(8) :: carth(1:n),f,d
		carth(1)=f
		if(n==1) then
			return
		endif
		do i=2,n
			carth(i)=carth(i-1)+d
		enddo
	end function
	function rarth(f,d,n)
		integer :: n,i
		real(8) :: rarth(1:n),f,d
		rarth(1)=f
		if(n==1) then
			return
		endif
		do i=2,n
			rarth(i)=rarth(i-1)+d
		enddo
	end function
	function iarth(f,n)
		integer :: n,i
		integer :: iarth(1:n),f
		iarth(1)=f
		if(n==1) then
			return
		endif
		do i=2,n
			iarth(i)=iarth(i-1)+1
		enddo
	end function
	function openfile(unit,file,access,recl)
		integer :: unit
		character(*) :: file
		character(100) :: nf
		logical :: openfile,flag
		character(*), optional :: access
		integer, optional :: recl
		!inquire(file=file, exist=flag) 
		write(*,"(A$)")"Using default '",file,"' or enter a new name: " 
		read(*,"(A)")nf
		if(len(trim(adjustl(nf)))>0) then
			if(present(access)) then
				open(unit,file=file(1:scan(file,".",.true.)-1)//"_"//trim(adjustl(nf))//".dat",access=access,recl=recl,form="formatted")
			else
				open(unit,file=file(1:scan(file,".",.true.)-1)//"_"//trim(adjustl(nf))//".dat")
			endif
		else
			if(present(access)) then
				open(unit,file=file,access=access,recl=recl,form="formatted")
			else
				open(unit,file=file)
			endif
		endif
		openfile=.true.
	end function
	function to_char_single(i,l) result(to_char)
		integer :: i
		integer, optional :: l
		integer :: j
		character(:), allocatable :: to_char
		j=i
		to_char=''
		do 
			to_char=char(48+mod(j,10))//to_char
			j=j/10
			if(j==0) then
				if(present(l)) then
					do
						if(l>len(to_char)) then
							to_char="0"//to_char
						else
							exit
						endif
					enddo
				endif
				exit
			endif
		enddo
	end function
	function to_char_array(i,l) result(to_char)
		integer :: i(:)
		integer, optional :: l
		integer :: j,n
		character(:), allocatable :: to_char(:)
		character(:), allocatable :: c
		to_char=[character::]
		do n=1,size(i)
			j=i(n)
			c=''
			do 
				c=char(48+mod(j,10))//c
				j=j/10
				if(j==0) then
					if(present(l)) then
						do
							if(l>len(c)) then
								c="0"//c
							else
								exit
							endif
						enddo
					endif
					exit
				endif
			enddo
			to_char=[to_char,c]
		enddo
	end function
	function fn(f)
		character(*) :: f
		character(:), allocatable :: fn
		character(50) :: nf
		character(50), save :: nf_save
		write(*,"(A$)")"Use default: '",f,"' or enter a new: "
		read(*,"(A$)")nf
		if(len(trim(adjustl(nf)))>0) then
			if(scan(nf,"/")/=0) then
				fn=f(1:scan(f,"/",.true.)-1)//trim(adjustl(nf))//".dat"
			elseif(scan(nf,"*")/=0) then
				fn=f(1:scan(f,".",.true.)-1)//"_"//trim(adjustl(nf_save))//".dat"
			else
				fn=f(1:scan(f,".",.true.)-1)//"_"//trim(adjustl(nf))//".dat"
				nf_save=nf
			endif
		else
			fn=f
		endif
	end function
	recursive subroutine rqsort(self,ord)
		!From: http://rosettacode.org/wiki/Sorting_algorithms/Quicksort#Fortran
		real(8) :: self(:)
		integer :: ord(:)
		integer :: left,right,n
		real(8) :: random
		real(8) :: pivot
		integer :: marker
		n=size(ord)
		if(n>1) then
			!call random_number(random)
			!pivot=self(int(random*real(n-1))+1)%val   !random pivor (not best performance, but avoids worst-case)
			pivot=self(ord((n+1)/2))   !random pivor (not best performance, but avoids worst-case)
			left=0
			right=n+1
			do while (left < right)
				right=right-1
				do while(self(ord(right))>pivot)
					right=right-1
				enddo
				left=left+1
				do while(self(ord(left))<pivot)
					left= left+1
				enddo
				if(left<right) then
					marker=ord(left)
					ord(left)=ord(right)
					ord(right)=marker
				endif
			end do
			if(left==right) then
				marker=left+1
			else
				marker=left
			endif
			call rqsort(self,ord(1:marker-1))
			call rqsort(self,ord(marker:n))
		endif
	end subroutine
	subroutine rcollect(self,a)
		real(8) :: self(:)
		integer :: i,j,tmp(10000)
		integer, allocatable :: a(:)
		tmp(1)=1
		j=1
		do i=2,size(self)
			if((abs(self(i)-self(i-1))>1d-6)) then
				j=j+1
				tmp(j)=i
			endif
		enddo
		j=j+1
		tmp(j)=i
		allocate(a(j))
		a=tmp(:j)
	end subroutine
	recursive subroutine rmqsort(self,ord)
		!From: http://rosettacode.org/wiki/Sorting_algorithms/Quicksort#Fortran
		real(8) :: self(:,:)
		integer :: ord(:)
		integer :: left,right,n
		real(8) :: random
		real(8) :: pivot(size(self,1))
		integer :: marker,i
		logical :: is_large,is_small
		n=size(ord)
		if(n>1) then
			!call random_number(random)
			!pivot=self(int(random*real(n-1))+1)%val   !random pivor (not best performance, but avoids worst-case)
			pivot=self(:,ord((n+1)/2))   !random pivor (not best performance, but avoids worst-case)
			left=0
			right=n+1
			do while (left < right)
				right=right-1
				do
					is_large=.false.
					do i=1,size(self,1)
						if((self(i,ord(right))-pivot(i))>1d-8) then
							is_large=.true.
							exit
						elseif((self(i,ord(right))-pivot(i))<-1d-8) then
							exit
						endif
					enddo
					if(.not.is_large) exit
					right=right-1
				enddo
				left=left+1
				do
					is_large=.true.
					do i=1,size(self,1)
						if((self(i,ord(left))-pivot(i))<-1d-8) then
							is_large=.false.
							exit
						elseif((self(i,ord(left))-pivot(i))>1d-8) then
							exit
						endif
					enddo
					if(is_large) exit
					left= left+1
				enddo
				if(left<right) then
					marker=ord(left)
					ord(left)=ord(right)
					ord(right)=marker
				endif
			end do
			if(left==right) then
				marker=left+1
			else
				marker=left
			endif
			call rmqsort(self,ord(1:marker-1))
			call rmqsort(self,ord(marker:n))
		endif
	end subroutine
	subroutine rmcollect(self,a)
		real(8) :: self(:,:)
		integer :: i,j,tmp(size(self,2))
		integer, allocatable :: a(:)
		tmp(1)=1
		j=1
		do i=2,size(self,2)
			if(any(abs(self(:,i)-self(:,i-1))>1d-8)) then
				j=j+1
				tmp(j)=i
			endif
		enddo
		j=j+1
		tmp(j)=i
		allocate(a(j))
		a=tmp(:j)
	end subroutine
	recursive subroutine sqsort(self)
		!From: http://rosettacode.org/wiki/Sorting_algorithms/Quicksort#Fortran
		class(t_sort) :: self(:)
		integer :: left,right,n
		real(8) :: random
		real(8) :: pivot
		integer :: marker
		n=size(self)
		if(n>1) then
			!call random_number(random)
			!pivot=self(int(random*real(n-1))+1)%val   !random pivor (not best performance, but avoids worst-case)
			pivot=self((size(self)+1)/2)%val   !random pivor (not best performance, but avoids worst-case)
			left=0
			right=n+1
			do while (left < right)
				right=right-1
				do while(self(right)%val>pivot)
					right=right-1
				enddo
				left=left+1
				do while(self(left)%val<pivot)
					left= left+1
				enddo
				if(left<right) then
					call self(left)%swap_sort(self(right))
				endif
			end do
			if(left==right) then
				marker=left+1
			else
				marker=left
			endif
			call sqsort(self(1:marker-1))
			call sqsort(self(marker:n))
		endif
	end subroutine
	subroutine scollect(self,a)
		class(t_sort) :: self(:)
		integer :: i,j,tmp(10000)
		integer, allocatable :: a(:)
		tmp(1)=1
		j=1
		do i=2,size(self)
			if((abs(self(i)%val-self(i-1)%val)>1d-6)) then
				j=j+1
				tmp(j)=i
			endif
		enddo
		j=j+1
		tmp(j)=i
		allocate(a(j))
		a=tmp(:j)
	end subroutine
	subroutine swap_sort(self,a)
		class(t_sort) :: self
		class(t_sort) :: a
		type(t_sort), allocatable :: tmp
		allocate(tmp)
		tmp%val=a%val
		a%val=self%val
		self%val=tmp%val
		deallocate(tmp)
	end subroutine
	elemental function abs2(a)
		complex(8), intent(in) :: a
		real(8) :: abs2
		abs2=real(a*conjg(a))
	end function
end module
