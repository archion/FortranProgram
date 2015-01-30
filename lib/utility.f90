module M_utility
	implicit none
	interface mwrite
		module procedure cmwrite,rmwrite,imwrite
	end interface
	interface arth
		module procedure carth,rarth,iarth
	end interface
	interface swap
		module procedure siswap,vcswap
	end interface
contains
	subroutine cmwrite(f,A)
		complex(8) :: A(:,:)
		integer :: f,i
		do i=1,size(A,1)
			!write(f,"(sp,e16.4,e11.4'I'$)")A(i,:)
			!write(f,"(sp,e10.2','e9.2$)")A(i,:)
			write(f,"(e16.4$)")real(A(i,:))
			write(f,"(1X)")
		enddo
		write(f,"(1X)")
	end subroutine
	subroutine rmwrite(f,A)
		real(8) :: A(:,:)
		integer :: f,i
		do i=1,size(A,1)
			write(f,"(e16.4$)")A(i,:)
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
	subroutine vcswap(a,b)
		complex(8) :: a(:),b(:),tmp(size(a))
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
	function openfile(i,f)
		integer :: i
		character(*) :: f
		character(100) :: nf
		logical :: openfile,flag
		inquire(file=f, exist=flag) 
		if(flag) then
			write(*,"(A$)")"File '",f,"' already exit, enter a new name: " 
			read(*,"(A)")nf
			if(len(trim(adjustl(nf)))>0) then
				open(i,file=f(1:scan(f,".",.true.)-1)//"_"//trim(adjustl(nf))//".dat")
			else
				open(i,file=f)
			endif
		else
			open(i,file=f)
		endif
		openfile=.true.
	end function
end module
