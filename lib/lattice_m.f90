module M_lattice_m
	use M_const
	use M_utility
	implicit none
	private
	type, extends(t_sort) :: t_mysort
		integer :: idx
		real(8) :: dr(3)
	contains
		procedure :: swap_sort => myswap
	end type
	type ::  t_bd
		integer :: i(2)
		real(8) :: r(3)
		real(8) :: dr(3)
		complex(8) :: bdc
	end type
	type t_nb
		type(t_bd), allocatable :: bd(:)
	end type
	type t_sb
		type(t_nb), allocatable :: nb(:)
	end type
	real(8) :: err=1d-6
	type t_latt
		real(8) :: T1(3),T2(3),a1(3),a2(3),a3(3)
		integer :: Ns,layer=1
		complex(8) :: bdc(3)=(/1d0,1d0,0d0/)
		type(t_sb), allocatable :: sb(:)
		real(8), allocatable :: sub(:,:)
		real(8), allocatable :: i2r(:,:)
	contains
		procedure gen_latt
		procedure gen_bond
		procedure gen_brizon
		procedure gen_origin_brizon
	end type
	type t_brizon
		integer :: n1,n2
		real(8) :: b1(3),b2(3)
		real(8), allocatable :: k(:,:)
		real(8), allocatable :: T(:,:)
	end type
	integer, allocatable :: nsub(:)
	type(t_latt), public :: latt
	type(t_brizon), public, target :: brizon,brizon_o
	public t_latt,t_brizon,check_lattice,sb,theta
contains
	subroutine myswap(self,a)
		class(t_mysort) :: self
		class(t_sort) :: a
		type(t_mysort), allocatable :: tmp
		select type(a)
		type is (t_mysort)
			call self%t_sort%swap_sort(a)
			allocate(tmp)

			tmp%idx=a%idx
			a%idx=self%idx
			self%idx=tmp%idx

			tmp%dr=a%dr
			a%dr=self%dr
			self%dr=tmp%dr

			deallocate(tmp)
		end select
	end subroutine
	function is_in(r,T,dr)
		logical :: is_in
		real(8) :: r(3),tr(2),T(:,:),T1(3),T2(3),tf(2,2),rerr(3)
		real(8), optional :: dr(:)
		integer :: i
		if(present(dr)) then
			dr=0d0
		endif
		is_in=.true.
		rerr=0d0
		do i=1,size(T,1)/2
			rerr=rerr+(T(i+1,:)-T(i,:))*err*i/1.34d0
		enddo
		do i=1,size(T,1),2
			T1=T(mod(i-2+size(T,1),size(T,1))+1,:)-T(i,:)
			T2=T(mod(i,size(T,1))+1,:)-T(i,:)
			tf=reshape((/T2(2),-T1(2),-T2(1),T1(1)/)/(T1(1)*T2(2)-T1(2)*T2(1)),(/2,2/))
			tr=matmul(tf,r(:2)-T(i,:2)+rerr(:2))
			if(any(tr<0d0)) then
				is_in=.false.
				if(present(dr)) then
					if(tr(1)<0d0) then
						do 
							dr=dr-(T2+2d0*T(i,:2))
							tr=matmul(tf,r(:2)+dr-T(i,:2)+rerr(:2))
							if(tr(1)>=0d0) then
								exit
							endif
						enddo
					endif
					if(tr(2)<0d0) then
						do 
							dr=dr-(T1+2d0*T(i,:2))
							tr=matmul(tf,r(:2)+dr-T(i,:2)+rerr(:2))
							if(tr(2)>=0d0) then
								exit
							endif
						enddo
					endif
				else
					exit
				endif
			endif
		enddo
	end function
	subroutine gen_grid(a1,a2,r0,T,tmp,n)
		integer :: n
		real(8) :: tmp(:,:),r0(:),a1(:),a2(:),T(:,:)
		real(8) :: r(3),x(size(T,1)),y(size(T,1)),tf(2,2)
		integer :: i,j
		tf=reshape((/a2(2),-a1(2),-a2(1),a1(1)/)/(a1(1)*a2(2)-a1(2)*a2(1)),(/2,2/))
		do i=1,size(T,1)
			x(i)=sum(tf(1,:)*(T(i,:2)-r0(:2)))
			y(i)=sum(tf(2,:)*(T(i,:2)-r0(:2)))
		enddo
		do j=nint(minval(y)),nint(maxval(y))
			do i=nint(minval(x)),nint(maxval(x))
				r=a1*i+a2*j+r0
				if(is_in(r,T)) then
					n=n+1
					tmp(n,:)=r
				endif
			enddo
		enddo
	end subroutine
	subroutine gen_latt(self)
		class(t_latt) :: self
		integer :: i,j,l,n
		real(8) :: x(4),y(4),r(3),tmp(10000000,3),tf(2,2),T(4,3)
		allocate(nsub(size(self%sub,1)+1))
		nsub=1
		associate(T1=>self%T1,T2=>self%T2,a1=>self%a1,a2=>self%a2,a3=>self%a3,Ns=>self%Ns,layer=>self%layer,bdc=>self%bdc)
			tf=reshape((/a2(2),-a1(2),-a2(1),a1(1)/)/(a1(1)*a2(2)-a1(2)*a2(1)),(/2,2/))
			T(1,:)=0d0
			T(2,:)=T1
			T(3,:)=T1+T2
			T(4,:)=T2
			n=0
			do l=1,size(self%sub,1)
				if(l>1) then
					nsub(l)=n+1
				endif
				call gen_grid(a1,a2,self%sub(l,:),T,tmp,n)
			enddo
			nsub(size(nsub))=n+1
			Ns=n
			allocate(self%i2r(Ns*layer,3))
			self%i2r(1:Ns,:)=tmp
			do i=2,layer
				do j=1,size(a3)
					self%i2r(1+Ns*(i-1):Ns*i,j)=self%i2r(:Ns,j)+a3(j)*(i-1)
				enddo
			enddo
			Ns=Ns*layer
		end associate
	end subroutine
	subroutine gen_bond(self,l)
		class(t_latt) :: self
		type(t_mysort) :: tmp(10000000)
		real(8) :: r(3),dr(3),T(4,3)
		complex(8) :: bdc_
		integer :: i,j,k,n,m,p,n1,n2,l,l1,l2,nr(2)
		integer, allocatable :: c(:)
		type(t_bd) ::  bd(10000000)
		associate(T1=>self%T1,T2=>self%T2,a1=>self%a1,a2=>self%a2,a3=>self%a3,Ns=>self%Ns,layer=>self%layer,bdc=>self%bdc)
			l1=max(abs(nint(sum(l*a1*T1)/sum(T1**2))),1)
			l2=max(abs(nint(sum(l*a2*T2)/sum(T2**2))),1)
			T(1,:)=0d0-(T1+T2)/2d0
			T(2,:)=T1-(T1+T2)/2d0
			T(3,:)=T1+T2-(T1+T2)/2d0
			T(4,:)=T2-(T1+T2)/2d0
			allocate(self%sb(size(latt%sub,1)))
			do n=1,size(nsub)-1
				r=self%i2r(nsub(n),:)
				k=0
				do j=1,Ns
					do n1=-l1,l1
						do n2=-l2,l2
							dr=self%i2r(j,:)-r+T1*n1+T2*n2
							k=k+1
							tmp(k)%val=sum(dr**2)
							tmp(k)%idx=j
							tmp(k)%dr=dr
						enddo
					enddo
				enddo
				call qsort(tmp(:k))
				call collect(tmp(:k),c)
				allocate(self%sb(n)%nb(0:l))
				do l1=1,l+1
					k=0
					do i=nsub(n),nsub(n+1)-1
					o:	do m=c(l1),c(l1+1)-1
							do p=c(l1),m-1
								if(sum(abs(tmp(m)%dr(:)+tmp(p)%dr(:)))<err.and.sb(tmp(m)%idx)==n) then
									cycle o
								endif
							enddo
							bdc_=1d0
							if(.not.is_in(latt%i2r(i,:)+tmp(m)%dr-(T1+T2)/2d0,T,dr)) then
								nr=-nint(matmul(reshape((/T2(2),-T1(2),-T2(1),T1(1)/)/(T1(1)*T2(2)-T1(2)*T2(1)),(/2,2/)),dr(:2)))
								if((abs(bdc(1))<err.and.nr(1)/=0).or.(abs(bdc(2))<err.and.nr(2)/=0)) then
									cycle
								else
									bdc_=bdc(1)**nr(1)*bdc(2)**nr(2)
								endif
							endif
							r=tmp(m)%dr+latt%i2r(i,:)+dr
							do j=1,Ns
								if(sum(abs(latt%i2r(j,:)-r))<err) then
									k=k+1
									bd(k)%i=(/i,j/)
									bd(k)%r=latt%i2r(i,:)+tmp(m)%dr/2d0
									bd(k)%dr=tmp(m)%dr
									bd(k)%bdc=bdc_
									exit
								endif
							enddo
						enddo o
					enddo
					allocate(self%sb(n)%nb(l1-1)%bd(k))
					self%sb(n)%nb(l1-1)%bd=bd(:k)
				enddo
				deallocate(c)
			enddo
		end associate
	end subroutine
	function sb(i)
		integer, intent(in) :: i
		integer :: sb
		integer :: n
		do n=1,size(nsub)-1
			if(i>=nsub(n).and.i<nsub(n+1)) then
				sb=n
				exit
			endif
		enddo
	end function
	subroutine gen_brizon(self,brizon)
		class(t_latt) :: self
		type(t_brizon) :: brizon
		real(8) :: tf(2,2),T(0:24,3),th,dr1(2),dr2(2)
		type(t_mysort) :: st(24)
		integer :: n,i,j
		integer, allocatable :: ist(:)
		associate(b1=>brizon%b1,b2=>brizon%b2,n1=>brizon%n1,n2=>brizon%n2,T1=>self%T1,T2=>self%T2,bdc=>self%bdc)
			if(allocated(brizon%k)) then
				deallocate(brizon%k)
			endif
			if(all(abs(bdc(:2))>err)) then
				allocate(brizon%k(n1*n2,3))
			elseif(all(abs(bdc(:2))<err)) then
				allocate(brizon%k(1,3))
			elseif(abs(bdc(1))<err) then
				allocate(brizon%k(n2,3))
				n1=1
			else
				allocate(brizon%k(n1,3))
				n2=1
			endif
			b1=(/-T2(2),T2(1),0d0/)
			b2=(/T1(2),-T1(1),0d0/)
			if(abs(bdc(1))<err) b2=T2
			if(abs(bdc(2))<err) b1=T1
			b1=2d0*pi*b1/sum(T1*b1)
			b2=2d0*pi*b2/sum(T2*b2)

			n=0
			do i=-2,2
				do j=-2,2
					if(i==0.and.j==0) then
						cycle
					endif
					n=n+1
					st(n)%dr=(b1*i+b2*j)/2d0
					st(n)%val=theta(st(n)%dr)
				enddo
			enddo
			call qsort(st(1:n))
			call collect(st(1:n),ist)
			do i=1,size(ist)-1
				st(i)=st(ist(i))
				do j=ist(i)+1,ist(i+1)-1
					if(sqrt(sum(st(j)%dr**2))<sqrt(sum(st(i)%dr**2))) then
						st(i)=st(j)
					endif
				enddo
			enddo
			n=1
			T(1,:2)=st(1)%dr(:2)
			do i=2,(size(ist)-1)/2
				dr1=T(n,:2)
				dr2=st(i)%dr(:2)
				th=theta(matmul(reshape((/dr2(2),-dr2(1),-dr1(2),dr1(1)/),(/2,2/)),(/sum(dr1**2),sum(dr2**2)/))/(dr1(1)*dr2(2)-dr1(2)*dr2(1)))
				if(sin(th-theta(dr1)-err)<0d0) then
					n=n-1
				endif
				if(sin(th-theta(dr2)+err)<0d0) then
					n=n+1
					T(n,:2)=dr2
				endif
			enddo
			T(n+1:,:)=-T(1:n,:)
			n=n*2
			allocate(brizon%T(n,3))
			do i=1,n
				j=mod(i,n)+1
				dr1=T(i,:2)
				dr2=T(j,:2)
				brizon%T(i,:2)=matmul(reshape((/dr2(2),-dr2(1),-dr1(2),dr1(1)/),(/2,2/)),(/sum(dr1**2),sum(dr2**2)/))/(dr1(1)*dr2(2)-dr1(2)*dr2(1))
				brizon%T(i,3)=0d0
			enddo

			n=0
			call gen_grid(b1/n1,b2/n2,(/0d0,0d0,0d0/),brizon%T,brizon%k,n)
			deallocate(ist)

			if(abs(bdc(1))<err) then
				b1=0d0
				deallocate(brizon%T)
				allocate(brizon%T(2,3))
				brizon%T(1,:)=-b2/2d0
				brizon%T(2,:)=b2/2d0
			endif
			if(abs(bdc(2))<err) then
				b2=0d0
				deallocate(brizon%T)
				allocate(brizon%T(2,3))
				brizon%T(1,:)=-b1/2d0
				brizon%T(2,:)=b1/2d0
			endif
		end associate
	end subroutine
	subroutine gen_origin_brizon(self,brizon_o)
		class(t_latt) :: self
		type(t_brizon) :: brizon_o
		type(t_latt) :: latt_
		real(8) :: dr(2),a(2,3),tmp(1000,3)
		integer :: nk,i,j,n
		nk=size(brizon%k,1)
		latt_%T1=self%a1
		latt_%T2=self%a2
		brizon_o%n1=1
		brizon_o%n2=1
		call latt_%gen_brizon(brizon_o)
		n=0
		call gen_grid(brizon%b1,brizon%b2,(/0d0,0d0,0d0/),brizon_o%T,tmp,n)
		deallocate(brizon_o%k)
		allocate(brizon_o%k(nk*n,3))
		brizon_o%k(1:nk,:)=brizon%k
		j=1
		do i=1,n
			if(sum(abs(tmp(i,:)))>err) then
				j=j+1
				brizon_o%k(1+nk*(j-1):nk*j,1)=brizon%k(:,1)+tmp(i,1)
				brizon_o%k(1+nk*(j-1):nk*j,2)=brizon%k(:,2)+tmp(i,2)
				brizon_o%k(1+nk*(j-1):nk*j,3)=brizon%k(:,3)+tmp(i,3)
			endif
		enddo
		do i=nk+1,size(brizon_o%k,1)
			if(.not.is_in(brizon_o%k(i,:),brizon_o%T,dr)) then
				brizon_o%k(i,:2)=brizon_o%k(i,:2)+dr
			endif
		enddo
	end subroutine
	function theta(r)
		real(8) :: r(:),theta,d
		d=sqrt(sum(r**2))
		if(d<1d-10) then
			theta=0d0
			return
		endif
		theta=acos(r(1)/d)
		if(r(2)<0d0) then
			theta=2d0*pi-theta
		endif
	end function
	subroutine check_lattice(ut)
		integer :: ut,k,l,m

		write(ut,"(3es13.2)")0d0,0d0,0d0
		write(ut,"(3es13.2)")brizon%b1
		write(ut,"(3es13.2)")brizon%b1+brizon%b2
		write(ut,"(3es13.2)")brizon%b2
		write(ut,"(3es13.2)")0d0,0d0,0d0
		write(ut,"(x/)")

		do k=1,size(brizon%T,1)+1
			write(ut,"(3es13.2)")brizon%T(mod(k-1,size(brizon%T,1))+1,:)
		enddo
		write(ut,"(x/)")

		do k=1,size(brizon%k,1)
			write(ut,"(es13.2$)")brizon%k(k,:),0d0,0d0,0d0
			write(ut,"(i7$)")k
			write(ut,"(x)")
		enddo
		write(ut,"(x/)")


		write(ut,"(3es13.2)")0d0,0d0,0d0
		write(ut,"(3es13.2)")latt%a1
		write(ut,"(3es13.2)")latt%a1+latt%a2
		write(ut,"(3es13.2)")latt%a2
		write(ut,"(3es13.2)")0d0,0d0,0d0
		write(ut,"(x/)")

		write(ut,"(3es13.2)")0d0,0d0,0d0
		write(ut,"(3es13.2)")latt%T1
		write(ut,"(3es13.2)")latt%T1+latt%T2
		write(ut,"(3es13.2)")latt%T2
		write(ut,"(3es13.2)")0d0,0d0,0d0
		write(ut,"(x/)")

		do m=1,size(latt%sb)
			do k=1,size(latt%sb(m)%nb(0)%bd)
				write(ut,"(es13.2$)")latt%sb(m)%nb(0)%bd(k)%r-latt%sb(m)%nb(0)%bd(k)%dr/2d0,latt%sb(m)%nb(0)%bd(k)%dr
				write(ut,"(i7$)")latt%sb(m)%nb(0)%bd(k)%i(1)
				write(ut,"(x)")
			enddo
		enddo
		write(ut,"(x/)")

		do m=1,size(latt%sb)
			do l=1,ubound(latt%sb(m)%nb,1)
				do k=1,size(latt%sb(m)%nb(l)%bd)
					write(ut,"(es13.2$)")latt%sb(m)%nb(l)%bd(k)%r-latt%sb(m)%nb(l)%bd(k)%dr/2d0,latt%sb(m)%nb(l)%bd(k)%dr
					write(ut,"(i7$)")k
					write(ut,"(x)")
				enddo
				write(ut,"(x/)")
			enddo
		enddo
	end subroutine
end module
