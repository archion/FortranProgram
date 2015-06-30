module M_utility
	use M_rd
	implicit none
	type t_sort
		real(8) :: val
	contains
		procedure :: swap_sort
		procedure :: qsort
		procedure :: collect
	end type
	interface mwrite
		module procedure cmwrite,rmwrite,imwrite
	end interface
	interface arth
		module procedure carth,rarth,iarth
	end interface
	interface swap
		module procedure siswap,srswap,vcswap,viswap,vrswap
	end interface
contains
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
	subroutine find_peak(pa,a,sg)
		real(8) :: pa(:),a
		integer :: sg
		sg=0
		if(pa(1)/=0d0) then
			if((pa(1)-pa(2))<0d0.and.(pa(2)-a)>0d0) then
				sg=1
			endif
			if((pa(1)-pa(2))>0d0.and.(pa(2)-a)<0d0) then
				sg=-1
			endif
		endif
		pa(1)=pa(2)
		pa(2)=a
	end subroutine
	subroutine find_cross(pa,a,sg)
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
		inquire(file=file, exist=flag) 
		if(flag) then
			write(*,"(A$)")"File '",file,"' already exit, enter a new name: " 
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
		else
			if(present(access)) then
				open(unit,file=file,access=access,recl=recl,form="formatted")
			else
				open(unit,file=file)
			endif
		endif
		openfile=.true.
	end function
	recursive subroutine qsort(self)
		!From: http://rosettacode.org/wiki/Sorting_algorithms/Quicksort#Fortran
		class(t_sort) :: self(:)
		integer :: left,right,n
		real(8) :: random
		real(8) :: pivot
		integer :: marker
		n=size(self)
		if(n>1) then
			call random_number(random)
			pivot=self(int(random*real(n-1))+1)%val   !random pivor (not best performance, but avoids worst-case)
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
			call self(1:marker-1)%qsort()
			call self(marker:n)%qsort()
		endif
	end subroutine
	subroutine collect(self,a)
		class(t_sort) :: self(:)
		integer :: i,j,tmp(1000)
		integer, allocatable :: a(:)
		tmp(1)=1
		j=1
		do i=2,size(self)
			if(((self(i)%val-self(i-1)%val)>1d-6)) then
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
	end subroutine
end module
