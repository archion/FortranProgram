module M_lattice
	use M_const
	use M_utility
	implicit none
	private
	type, extends(t_sort) :: t_mysort
		integer :: idx
		complex(8) :: bdc
		integer :: bond
		real(8) :: dr(3)
	contains
		procedure :: swap_sort => myswap
	end type
	type t_nb
		integer, allocatable :: neb(:)
		integer, allocatable :: bond(:)
		complex(8), allocatable :: bdc(:)
		real(8), allocatable :: dr(:,:)
	end type
	type t_neb
		type(t_nb), allocatable :: nb(:)
	end type
	type ::  t_bond
		integer :: i(2)
		real(8) :: r(3)
		real(8) :: dr(3)
		complex(8) :: bdc
	end type
	type t_bd
		type(t_bond), allocatable :: bd(:)
	end type
	real(8) :: err=1d-6
	type t_latt
		real(8) :: T1(3),T2(3),a1(3),a2(3),a3(3)
		integer :: Ns,layer
		complex(8) :: bdc(3)
		type(t_neb), allocatable :: neb(:)
		type(t_bd), allocatable :: bond(:)
		real(8), allocatable :: sub(:,:)
		real(8), allocatable :: i2r(:,:)
	contains
		procedure gen_latt
		procedure gen_neb
		procedure gen_bond
		procedure gen_brizon
	end type
	type t_brizon
		integer :: n1,n2
		real(8) :: b1(3),b2(3)
		real(8), allocatable :: k(:,:)
		real(8), allocatable :: T(:,:)
	end type
	integer :: nab
	type(t_latt), public :: latt
	type(t_brizon), public :: brizon
	public theta,t_latt,check_lattice,ab
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

			tmp%bdc=a%bdc
			a%bdc=self%bdc
			self%bdc=tmp%bdc

			tmp%dr=a%dr
			a%dr=self%dr
			self%dr=tmp%dr

			tmp%bond=a%bond
			a%bond=self%bond
			self%bond=tmp%bond
			deallocate(tmp)
		end select
	end subroutine
	!function is_in(r,T1,T2)
		!logical :: is_in
		!real(8) :: r(3),tr(2),T1(3),T2(3),tf(2,2)
		!tf=reshape((/T2(2),-T1(2),-T2(1),T1(1)/)/(T1(1)*T2(2)-T1(2)*T2(1)),(/2,2/))
		!tr=matmul(tf,r(:2))
		!if(all(tr>=-err.and.tr<(1d0-err))) then
			!is_in=.true.
		!else
			!is_in=.false.
		!endif
	!end function
	function is_in(r,T)
		logical :: is_in
		real(8) :: r(3),tr(2),T(:,:),T1(3),T2(3),tf(2,2),rerr(3)
		integer :: i
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
				exit
			endif
		enddo
	end function
	subroutine gen_latt(self)
		class(t_latt) :: self
		integer :: i,j,l,n
		real(8) :: x(4),y(4),r(3),tmp(10000000,3),tf(2,2),T(4,3)
		associate(T1=>self%T1,T2=>self%T2,a1=>self%a1,a2=>self%a2,a3=>self%a3,Ns=>self%Ns,layer=>self%layer,bdc=>self%bdc)
			tf=reshape((/a2(2),-a1(2),-a2(1),a1(1)/)/(a1(1)*a2(2)-a1(2)*a2(1)),(/2,2/))
			T(1,:)=0d0
			T(2,:)=T1
			T(3,:)=T1+T2
			T(4,:)=T2
			n=0
			do l=1,size(self%sub,1)
				do i=1,4
					x(i)=sum(tf(1,:)*T(i,:2))-sum(tf(1,:)*self%sub(l,:2))
					y(i)=sum(tf(2,:)*T(i,:2))-sum(tf(2,:)*self%sub(l,:2))
				enddo
				do j=nint(minval(y)),nint(maxval(y))
					do i=nint(minval(x)),nint(maxval(x))
						r=a1*i+a2*j+self%sub(l,:)
						!if(is_in(r,T1,T2)) then
						if(is_in(r,T)) then
							n=n+1
							tmp(n,:)=r
						endif
					enddo
				enddo
				if(l==1) then
					nab=n
				endif
			enddo
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
	subroutine gen_neb(self,l)
		class(t_latt) :: self
		type(t_mysort) :: tmp(10000000)
		real(8) :: p,dr(3)
		integer :: i,j,k,n,n1,n2,l,sg(2),l1,l2
		integer, allocatable :: t(:)
		associate(T1=>self%T1,T2=>self%T2,a1=>self%a1,a2=>self%a2,a3=>self%a3,Ns=>self%Ns,layer=>self%layer,bdc=>self%bdc)
			l1=max(abs(nint(sum(l*a1*T1)/sum(T1**2))),1)
			l2=max(abs(nint(sum(l*a2*T2)/sum(T2**2))),1)
			allocate(self%neb(Ns))
			do i=1,Ns/layer
				k=0
				do j=1,Ns/layer
					do n1=-l1,l1
						do n2=-l2,l2
							if((abs(bdc(1))<err.and.n1/=0).or.(abs(bdc(2))<err.and.n2/=0)) then
								cycle
							endif
							dr=self%i2r(j,:)-self%i2r(i,:)+T1*n1+T2*n2

							k=k+1
							!if(.not.is_in(dr+T1/2d0+T2/2d0,T1,T2,sg,tmp(k)%bdc)) then
								!!dr=dr+T1*sg(1)+T2*sg(2)
							!endif
							tmp(k)%val=sum(dr**2)
							tmp(k)%idx=j
							tmp(k)%dr=dr
							tmp(k)%bdc=1d0
							select case(n1)
							case(1:)
								tmp(k)%bdc=tmp(k)%bdc*bdc(1)**n1
							case(:-1)
								tmp(k)%bdc=tmp(k)%bdc*conjg(bdc(1))**n1
							end select
							select case(n2)
							case(1:)
								tmp(k)%bdc=tmp(k)%bdc*bdc(2)**n2
							case(:-1)
								tmp(k)%bdc=tmp(k)%bdc*conjg(bdc(2))**n2
							end select
						enddo
					enddo
				enddo
				call tmp(:k)%qsort()
				call tmp(:k)%collect(t)
				do k=1,layer
					allocate(self%neb(i+Ns/layer*(k-1))%nb(0:l))
				enddo
				do j=1,l+1
					allocate(self%neb(i)%nb(j-1)%neb(t(j+1)-t(j)))
					allocate(self%neb(i)%nb(j-1)%bond(t(j+1)-t(j)))
					allocate(self%neb(i)%nb(j-1)%bdc(t(j+1)-t(j)))
					allocate(self%neb(i)%nb(j-1)%dr(t(j+1)-t(j),3))
					do k=t(j),t(j+1)-1
						tmp(k)%val=theta(tmp(k)%dr)
					enddo
					call tmp(t(j):t(j+1)-1)%qsort()
					self%neb(i)%nb(j-1)%neb(:)=tmp(t(j):t(j+1)-1)%idx
					self%neb(i)%nb(j-1)%bdc(:)=tmp(t(j):t(j+1)-1)%bdc
					self%neb(i)%nb(j-1)%dr(:,1)=tmp(t(j):t(j+1)-1)%dr(1)
					self%neb(i)%nb(j-1)%dr(:,2)=tmp(t(j):t(j+1)-1)%dr(2)
					self%neb(i)%nb(j-1)%dr(:,3)=tmp(t(j):t(j+1)-1)%dr(3)
					do k=2,layer
						allocate(self%neb(i+Ns/layer*(k-1))%nb(j-1)%neb(t(j+1)-t(j)))
						allocate(self%neb(i+Ns/layer*(k-1))%nb(j-1)%bond(t(j+1)-t(j)))
						allocate(self%neb(i+Ns/layer*(k-1))%nb(j-1)%bdc(t(j+1)-t(j)))
						allocate(self%neb(i+Ns/layer*(k-1))%nb(j-1)%dr(t(j+1)-t(j),3))
						self%neb(i+Ns/layer*(k-1))%nb(j-1)%neb(:)=self%neb(i)%nb(j-1)%neb(:)+Ns/layer*(k-1)
						self%neb(i+Ns/layer*(k-1))%nb(j-1)%bdc(:)=self%neb(i)%nb(j-1)%bdc(:)
						self%neb(i+Ns/layer*(k-1))%nb(j-1)%dr(:,1)=self%neb(i)%nb(j-1)%dr(:,1)
						self%neb(i+Ns/layer*(k-1))%nb(j-1)%dr(:,2)=self%neb(i)%nb(j-1)%dr(:,2)
						self%neb(i+Ns/layer*(k-1))%nb(j-1)%dr(:,3)=self%neb(i)%nb(j-1)%dr(:,3)
					enddo
				enddo
				deallocate(t)
			enddo
		end associate
	end subroutine
	subroutine gen_bond(self,l)
		class(t_latt) :: self
		integer :: l,i,j,k,n,m,p
		type(t_bond) :: tmp(10000)
		associate(T1=>self%T1,T2=>self%T2,a1=>self%a1,a2=>self%a2,a3=>self%a3,Ns=>self%Ns,layer=>self%layer,bdc=>self%bdc)
			allocate(self%bond(0:l))
			do k=0,l
				n=0
				do i=1,Ns
					if(size(self%neb(i)%nb(1:))<k) then
						cycle
					endif
				o:	do m=1,size(self%neb(i)%nb(k)%neb)
						j=self%neb(i)%nb(k)%neb(m)
						do p=1,m-1
							if(abs(sum(self%neb(i)%nb(k)%dr(p,:)+self%neb(i)%nb(k)%dr(m,:)))<err) then
								cycle o
							endif
						enddo
						n=n+1
						tmp(n)%i=(/i,j/)
						tmp(n)%r=self%i2r(i,:)+self%neb(i)%nb(k)%dr(m,:)/2d0
						tmp(n)%dr=self%neb(i)%nb(k)%dr(m,:)
						tmp(n)%bdc=self%neb(i)%nb(k)%bdc(m)
						self%neb(i)%nb(k)%bond(m)=n
						do p=1,size(self%neb(j)%nb(k)%neb)
							if(self%neb(j)%nb(k)%neb(p)==i) then
								self%neb(j)%nb(k)%bond(p)=n
								exit
							endif
						enddo
					enddo o
				enddo
				allocate(self%bond(k)%bd(n))
				self%bond(k)%bd=tmp(:n)
			enddo
		end associate
	end subroutine
	subroutine gen_brizon(self)
		class(t_latt) :: self
		real(8) :: r(3),tf(2,2),T(0:24,3),th,dr1(2),dr2(2)
		real(8), allocatable :: x(:),y(:)
		type(t_mysort) :: st(24)
		integer :: n,i,j
		integer, allocatable :: ist(:)
		associate(b1=>brizon%b1,b2=>brizon%b2,n1=>brizon%n1,n2=>brizon%n2,T1=>self%T1,T2=>self%T2,bdc=>self%bdc)
			if(all(abs(bdc(:2))>err)) then
				allocate(brizon%k(n1*n2,3))
			elseif(all(abs(bdc(:2))<err)) then
				allocate(brizon%k(1,3))
			elseif(abs(bdc(1))<err) then
				allocate(brizon%k(n2,3))
			else
				allocate(brizon%k(n1,3))
			endif
			b1=(/-T2(2),T2(1),0d0/)
			b2=(/T1(2),-T1(1),0d0/)
			if(abs(bdc(1))<err) b2=T2
			if(abs(bdc(2))<err) b1=T1
			b1=2d0*pi*b1/sum(T1*b1)
			b2=2d0*pi*b2/sum(T2*b2)
			tf=reshape((/b2(2)*n1,-b1(2)*n2,-b2(1)*n1,b1(1)*n2/)/(b1(1)*b2(2)-b1(2)*b2(1)),(/2,2/))
			if(abs(bdc(1))<err) tf(1,:)=0d0
			if(abs(bdc(2))<err) tf(2,:)=0d0
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
			call st(1:n)%qsort()
			call st(1:n)%collect(ist)
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
			allocate(brizon%T(n,3),x(n),y(n))
			do i=1,n
				j=mod(i,n)+1
				dr1=T(i,:2)
				dr2=T(j,:2)
				brizon%T(i,:2)=matmul(reshape((/dr2(2),-dr2(1),-dr1(2),dr1(1)/),(/2,2/)),(/sum(dr1**2),sum(dr2**2)/))/(dr1(1)*dr2(2)-dr1(2)*dr2(1))
				brizon%T(i,3)=0d0
			enddo
			do i=1,size(x)
				x(i)=sum(tf(1,:)*brizon%T(i,:2))
				y(i)=sum(tf(2,:)*brizon%T(i,:2))
			enddo
			n=0
			do j=nint(minval(y)),nint(maxval(y))
				do i=nint(minval(x)),nint(maxval(x))
					r=(b1/n1*i+b2/n2*j)
					if(is_in(r,brizon%T)) then
						n=n+1
						brizon%k(n,:)=r
					endif
				enddo
			enddo
			deallocate(x,y,ist)
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
	function ab(i)
		integer, intent(in) :: i
		integer :: ab
		if((i<=nab).or.(latt%Ns/latt%layer<i.and.i<=latt%Ns/latt%layer+nab)) then
			ab=1
		else
			ab=-1
		endif
	end function
	subroutine check_lattice(ut)
		integer :: ut,k,l

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

		do l=0,ubound(latt%bond,1)
			do k=1,size(latt%bond(l)%bd)
				write(ut,"(es13.2$)")latt%bond(l)%bd(k)%r-latt%bond(l)%bd(k)%dr/2d0,latt%bond(l)%bd(k)%dr
				write(ut,"(i7$)")k
				write(ut,"(x)")
			enddo
			write(ut,"(x/)")
		enddo
	end subroutine
end module
