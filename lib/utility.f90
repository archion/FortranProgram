module M_utility
	USE IFPORT, ifort_qsort => qsort
	use omp_lib
	use mkl_service
	use M_rd
	use M_const
	implicit none
	type t_sort
		real(wp) :: val
	contains
		procedure :: swap_sort
	end type
	interface integrate
		module procedure rintegrate,cintegrate,rfintegrate,cfintegrate
	end interface
	interface set_grid
		module procedure rset_grid,cset_grid
	end interface
	interface get
		module procedure cget,rget
	end interface
	interface find
		module procedure rfind
	end interface
	interface mwrite
		module procedure cmwrite,rmwrite,imwrite
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
	real(wp) :: otime(0:6)=0._wp
contains
	logical function next(x,y,dx,xtol,ytol) result(rt)
		real(wp) :: x,y
		real(wp), optional :: dx,xtol,ytol
		real(wp), save :: px,py,dx_,xtol_,ytol_
		logical, save :: flag
		rt=.false.
		if(present(dx)) then 
			py=nan
			flag=.false.
			dx_=dx
			if(present(xtol)) then
				xtol_=xtol
			else
				xtol_=huge(1._wp)
			endif
			if(present(ytol)) then
				ytol_=ytol
			else
				ytol_=huge(1._wp)
			endif
			return
		endif
		if(.not.isnan(py)) then
			if(py*y<0._wp) flag=.true.
			if(flag) then
				dx_=dx_/2._wp
			endif
			if(py*y<0._wp) then
				dx_=-dx_
			elseif(sign(1._wp,py)*sign(1._wp,(py-y)/(px-x)*dx_)>0._wp) then 
				x=px
				dx_=-dx_
			endif
		endif
		px=x
		py=y
		!if(abs(dx_)<tol_) then
		if(abs(y)<ytol_.and.abs(dx_)<xtol_) then
			rt=.true.
		else
			x=x+dx_
		endif
	end function
	real(wp) function rintegrate(x,A) result(f)
		real(wp) :: x(:)
		real(wp) :: A(:)
		integer :: i
		f=0._wp
		do i=1,size(x)-1
			f=f+0.5_wp*(x(i+1)-x(i))*(A(i+1)+A(i))
		enddo
	end function
	complex(wp) function cintegrate(x,A) result(f)
		real(wp) :: x(:)
		complex(wp) :: A(:)
		integer :: i
		f=0._wp
		do i=1,size(x)-1
			f=f+0.5_wp*(x(i+1)-x(i))*(A(i+1)+A(i))
		enddo
	end function
	complex(wp) function cfintegrate(f,bound,tol,singular) result(rt)
		complex(wp) :: f
		real(wp) :: bound(2),tol
		logical, optional :: singular
		integer, parameter :: lim=300,lenw=lim*4
		real(dp) :: abserr,work(lenw),rrt,crt
		integer :: neval,ier,last,iwork(lim)
		external f
		if(bound(1)/=-inf.and.bound(2)/=inf) then
			if(present(singular)) then
				if(singular) then
					call dqags(rf,real(bound(1),kind=dp),real(bound(2),kind=dp),real(tol,kind=dp),real(tol,kind=dp),rrt,abserr,neval,ier,lim,lenw,last,iwork,work)
					call dqags(cf,real(bound(1),kind=dp),real(bound(2),kind=dp),real(tol,kind=dp),real(tol,kind=dp),crt,abserr,neval,ier,lim,lenw,last,iwork,work)
				else
					call dqag(rf,real(bound(1),kind=dp),real(bound(2),kind=dp),real(tol,kind=dp),real(tol,kind=dp),1,rrt,abserr,neval,ier,lim,lenw,last,iwork,work)
					call dqag(cf,real(bound(1),kind=dp),real(bound(2),kind=dp),real(tol,kind=dp),real(tol,kind=dp),1,crt,abserr,neval,ier,lim,lenw,last,iwork,work)
				endif
			else
				call dqag(rf,real(bound(1),kind=dp),real(bound(2),kind=dp),real(tol,kind=dp),real(tol,kind=dp),1,rrt,abserr,neval,ier,lim,lenw,last,iwork,work)
				call dqag(cf,real(bound(1),kind=dp),real(bound(2),kind=dp),real(tol,kind=dp),real(tol,kind=dp),1,crt,abserr,neval,ier,lim,lenw,last,iwork,work)
			endif
		elseif(bound(1)/=-inf) then
			call dqagi(rf,real(bound(1),kind=dp),1,real(tol,kind=dp),real(tol,kind=dp),rrt,abserr,neval,ier,lim,lenw,last,iwork,work)
			call dqagi(cf,real(bound(1),kind=dp),1,real(tol,kind=dp),real(tol,kind=dp),crt,abserr,neval,ier,lim,lenw,last,iwork,work)
		elseif(bound(2)/=inf) then
			call dqagi(rf,real(bound(2),kind=dp),-1,real(tol,kind=dp),real(tol,kind=dp),rrt,abserr,neval,ier,lim,lenw,last,iwork,work)
			call dqagi(cf,real(bound(2),kind=dp),-1,real(tol,kind=dp),real(tol,kind=dp),crt,abserr,neval,ier,lim,lenw,last,iwork,work)
		else
			call dqagi(rf,real(bound(2),kind=dp),2,real(tol,kind=dp),real(tol,kind=dp),rrt,abserr,neval,ier,lim,lenw,last,iwork,work)
			call dqagi(cf,real(bound(2),kind=dp),2,real(tol,kind=dp),real(tol,kind=dp),crt,abserr,neval,ier,lim,lenw,last,iwork,work)
		endif
		rt=cmplx(rrt,crt,kind=wp)
	contains
		real(dp) function rf(x) result(y)
			real(dp) :: x
			y=real(f(real(x,kind=wp)),kind=dp)
		end function
		real(dp) function cf(x) result(y)
			real(dp) :: x
			y=real(imag(f(real(x,kind=wp))),kind=dp)
		end function
	end function
	real(wp) function rfintegrate(f,bound,tol,singular) result(rt)
		real(wp) :: f,bound(2),tol
		logical, optional :: singular
		integer, parameter :: lim=300,lenw=lim*4
		real(dp) :: abserr,work(lenw),rt_
		integer :: neval,ier,last,iwork(lim)
		external f
		if(bound(1)/=-inf.and.bound(2)/=inf) then
			if(present(singular)) then
				if(singular) then
					call dqags(f_,real(bound(1),kind=dp),real(bound(2),kind=dp),real(tol,kind=dp),real(tol,kind=dp),rt_,abserr,neval,ier,lim,lenw,last,iwork,work)
				else
					call dqag(f_,real(bound(1),kind=dp),real(bound(2),kind=dp),real(tol,kind=dp),real(tol,kind=dp),1,rt_,abserr,neval,ier,lim,lenw,last,iwork,work)
				endif
			else
				call dqag(f_,real(bound(1),kind=dp),real(bound(2),kind=dp),real(tol,kind=dp),real(tol,kind=dp),1,rt_,abserr,neval,ier,lim,lenw,last,iwork,work)
			endif
		elseif(bound(1)/=-inf) then
			call dqagi(f_,real(bound(1),kind=dp),1,real(tol,kind=dp),real(tol,kind=dp),rt_,abserr,neval,ier,lim,lenw,last,iwork,work)
		elseif(bound(2)/=inf) then
			call dqagi(f_,real(bound(2),kind=dp),-1,real(tol,kind=dp),real(tol,kind=dp),rt_,abserr,neval,ier,lim,lenw,last,iwork,work)
		else
			call dqagi(f_,real(bound(2),kind=dp),2,real(tol,kind=dp),real(tol,kind=dp),rt_,abserr,neval,ier,lim,lenw,last,iwork,work)
		endif
		rt=real(rt_,kind=wp)
	contains
		real(dp) function f_(x) result(y)
			real(dp) :: x
			y=real(f(real(x,kind=wp)),kind=dp)
		end function
	end function
	pure elemental real(wp) function ff(x)
		real(wp), intent(in) :: x
		ff=1._wp/(exp(x)+1._wp)
	end function
	pure elemental real(wp) function fb(x)
		real(wp), intent(in) :: x
		integer :: m
		real(wp) :: f1,f2
		if(abs(x)<1e-20_wp) then
			fb=-0.5_wp
		elseif(x<0._wp) then
			fb=1._wp/(exp(x)-1._wp)
		else
			fb=-exp(-x)/(exp(-x)-1._wp)
		endif
	end function
	subroutine cset_grid(grid,flag,from,to,n)
		complex(wp) :: grid(:)
		real(wp) :: from(:),to(:)
		integer :: flag(:),n(:)
		integer,save :: n_
		integer :: i,j,k,l,u,m,o
		real(wp) :: base,from_,dx
		m=1
		o=1
		do k=1,size(flag)
			select case(flag(k))
			case(reset)
				n_=0
				grid=0._wp
			case(add_log_linear)
				l=find(grid(1:n_)%re,from(o))
				u=find(grid(1:n_)%re,to(o))
				grid(u+n(m)*n(m+1)+1:n_+n(m)*n(m+1))=grid(u+1:n_)
				u=u+n(m)*n(m+1)
				n_=n_+n(m)*n(m+1)
				base=(to(o)/from(o))**(1._wp/n(m))
				from_=from(o)/base
				do i=l+1,u
					if(mod(i-l-1,n(m+1))==0) then
						from_=from_*base
						dx=from_*(base-1._wp)/n(m+1)
					endif
					grid(i)=from_+dx*(mod(i-l-1,n(m+1))+1)
				enddo
				m=m+2
				o=o+1
			case(add_linear)
				l=find(grid(1:n_)%re,from(o))
				u=find(grid(1:n_)%re,to(o))
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
	subroutine rset_grid(grid,flag,from,to,n)
		real(wp) :: grid(:),from(:),to(:)
		integer :: flag(:),n(:)
		integer,save :: n_
		integer :: i,j,k,l,u,m,o
		real(wp) :: base,from_,dx
		m=1
		o=1
		do k=1,size(flag)
			select case(flag(k))
			case(reset)
				n_=0
				grid=0._wp
			case(add_log_linear)
				l=find(grid(1:n_),from(o))
				u=find(grid(1:n_),to(o))
				grid(u+n(m)*n(m+1)+1:n_+n(m)*n(m+1))=grid(u+1:n_)
				u=u+n(m)*n(m+1)
				n_=n_+n(m)*n(m+1)
				base=(to(o)/from(o))**(1._wp/n(m))
				from_=from(o)/base
				do i=l+1,u
					if(mod(i-l-1,n(m+1))==0) then
						from_=from_*base
						dx=from_*(base-1._wp)/n(m+1)
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
		real(wp), intent(in) :: A(:)
		real(wp), intent(in) :: x
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
	pure complex(wp) function cget(x,i) result(rt)
		complex(wp), intent(in) :: x(:)
		integer, intent(in) :: i
		if(i<=size(x,1)) then
			rt=x(i)
		else
			rt=nan
		endif
	end function
	pure real(wp) function rget(x,i) result(rt)
		real(wp), intent(in) :: x(:)
		integer, intent(in) :: i
		if(i<=size(x,1)) then
			rt=x(i)
		else
			rt=nan
		endif
	end function
	real(wp) function for_in(x,i,id) result(rt)
		real(wp), optional, intent(in) :: x(:)
		integer, optional, intent(in) :: i
		integer, intent(in) :: id
		integer, save :: i_(20)=0
		integer, save :: n_(20)=0
		real(wp), save :: x_(200,20)
		if(present(i)) then
			if(i>0) then
				if(i<=n_(id)) then
					rt=x_(i,id)
				else
					rt=nan
				endif
			else
				if(i+i_(id)>0) then
					rt=x_(i+i_(id),id)
				else
					rt=nan
				endif
			endif
		else
			if(present(x)) then
				n_(id)=size(x,1)
				x_(1:n_(id),id)=x
			else
				i_(id)=0
			endif
			i_(id)=i_(id)+1
			if(i_(id)<=n_(id)) then
				rt=x_(i_(id),id)
			else
				rt=nan
			endif
		endif
	end function
	subroutine show_time()
		integer :: i
		write(*,"(A)")"The timing results is:"
		do i=1,ubound(otime,1)
			if(otime(i)>1e-6_wp) then
				write(*,"(i4,'m',f6.3,'s')")int(otime(i)/60._wp),otime(i)-int(otime(i)/60._wp)*60._wp
			endif
		enddo
	end subroutine
	subroutine cmwrite(f,A)
		complex(wp) :: A(:,:)
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
		real(wp) :: A(:,:)
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
		real(wp) :: a,b,tmp
		tmp=a
		a=b
		b=tmp
	end subroutine
	subroutine vrswap(a,b)
		real(wp) :: a(:),b(:),tmp(size(a))
		tmp=a
		a=b
		b=tmp
	end subroutine
	subroutine vcswap(a,b)
		complex(wp) :: a(:),b(:),tmp(size(a))
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
		real(wp) :: a(:)
		integer, allocatable :: ipeak(:)
		integer :: ipeak_(10),n,i
		n=0
		do i=2,size(a)-1
			if((a(i)-a(i-1))>0._wp.and.(a(i)-a(i+1))>0._wp) then
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
		real(wp) :: pa,a
		integer :: sg
		sg=0
		if(pa/=0._wp) then
			if(pa>0._wp.and.a<0._wp) then
				sg=1
			endif
			if(pa<0._wp.and.a>0._wp) then
				sg=-1
			endif
		endif
		pa=a
	end subroutine
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
		!ifort error
		!to_char=[character::]
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
		real(wp) :: self(:)
		integer :: ord(:)
		integer :: left,right,n
		real(wp) :: random
		real(wp) :: pivot
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
		real(wp) :: self(:)
		integer :: i,j,tmp(size(self))
		integer, allocatable :: a(:)
		tmp(1)=1
		j=1
		do i=2,size(self)
			if((abs(self(i)-self(i-1))>eps*100._wp)) then
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
		real(wp) :: self(:,:)
		integer :: ord(:)
		integer :: left,right,n
		real(wp) :: random
		real(wp) :: pivot(size(self,1))
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
						if((self(i,ord(right))-pivot(i))>eps*100._wp) then
							is_large=.true.
							exit
						elseif((self(i,ord(right))-pivot(i))<-eps*100._wp) then
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
						if((self(i,ord(left))-pivot(i))<-eps*100._wp) then
							is_large=.false.
							exit
						elseif((self(i,ord(left))-pivot(i))>eps*100._wp) then
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
		real(wp) :: self(:,:)
		integer :: i,j,tmp(size(self,2))
		integer, allocatable :: a(:)
		tmp(1)=1
		j=1
		do i=2,size(self,2)
			if(any(abs(self(:,i)-self(:,i-1))>eps*100._wp)) then
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
		real(wp) :: random
		real(wp) :: pivot
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
		complex(wp), intent(in) :: a
		real(wp) :: abs2
		abs2=real(a*conjg(a))
	end function
end module
