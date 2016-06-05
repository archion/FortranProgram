module M_lattice_test1
	use M_const
	use M_utility
	implicit none
	private
	type, extends(t_sort) :: t_mysort
		integer :: i(2)
		real(8) :: dr(3)
		complex(8) :: bdc
	contains
		procedure :: swap_sort => myswap
	end type
	type ::  t_bd
		integer :: sb(2)
		integer :: i(2)
		real(8) :: r(3)
		real(8) :: dr(3)
		complex(8) :: bdc
	end type
	type ::  t_st
		integer, allocatable :: j(:)
		integer, allocatable :: bd(:)
	end type
	type t_nb
		type(t_bd), allocatable :: bd(:)
		type(t_st), allocatable :: st(:)
	end type
	real(8) :: err=1d-6
	type t_latt
		real(8) :: T1(3),T2(3),a1(3),a2(3),a3(3),c1(3),c2(3),c3(3)
		integer :: Ns,layer=1,n1=1,n2=1
		complex(8) :: bdc(3)=(/1d0,1d0,0d0/)
		type(t_nb), allocatable :: nb(:)
		real(8), allocatable :: rsb(:,:)
	contains
		procedure gen_latt
		procedure gen_brizon
	end type
	type t_brizon
		real(8) :: b1(3),b2(3)
		real(8), allocatable :: T(:,:)
		real(8), allocatable :: k(:,:)
		real(8), allocatable :: ok(:,:)
	end type
	type(t_latt), public :: latt
	type(t_brizon), public, target :: brizon
	public t_latt,t_brizon,check_lattice,theta
contains
	subroutine myswap(self,a)
		class(t_mysort) :: self
		class(t_sort) :: a
		type(t_mysort), allocatable :: tmp
		select type(a)
		type is (t_mysort)
			call self%t_sort%swap_sort(a)
			allocate(tmp)

			tmp%i=a%i
			a%i=self%i
			self%i=tmp%i

			tmp%dr=a%dr
			a%dr=self%dr
			self%dr=tmp%dr

			tmp%bdc=a%bdc
			a%bdc=self%bdc
			self%bdc=tmp%bdc
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
							dr=dr-(T2-T(mod(i+size(T,1)/2-1,size(T,1))+1,:2)+T(i,:2))
							tr=matmul(tf,r(:2)+dr-T(i,:2)+rerr(:2))
							if(tr(1)>=0d0) then
								exit
							endif
						enddo
					endif
					if(tr(2)<0d0) then
						do 
							dr=dr-(T1-T(mod(i+size(T,1)/2-1,size(T,1))+1,:2)+T(i,:2))
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
	subroutine gen_grid(a1,a2,r0,T,grid)
		real(8), allocatable :: grid(:,:)
		real(8) :: r0(:),a1(:),a2(:),T(:,:)
		real(8) :: r(3),x(size(T,1)),y(size(T,1)),tf(2,2),tmp(10000000,3)
		integer :: i,j,n
		n=0
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
		if(allocated(grid)) deallocate(grid)
		allocate(grid(n,3))
		grid=tmp(:n,:)
	end subroutine
	subroutine gen_latt(self,nb)
		class(t_latt) :: self
		integer :: nb
		integer :: i,j,l,l1,l2,n1,n2,m1,m2,n,m,p,k,sb
		real(8) :: r(3),tf(2,2),T(4,3),dr(3)
		real(8), allocatable :: tmp(:,:)
		type(t_mysort) :: sort(10000000)
		type(t_bd) ::  bd(10000000)
		integer ::  st(5000,0:20,2)
		integer, allocatable :: c(:)
		sb=size(latt%rsb,1)
		associate(T1=>self%T1,T2=>self%T2,a1=>self%a1,a2=>self%a2,a3=>self%a3,Ns=>self%Ns,layer=>self%layer,bdc=>self%bdc)
			tf=reshape((/a2(2),-a1(2),-a2(1),a1(1)/)/(a1(1)*a2(2)-a1(2)*a2(1)),(/2,2/))
			T(1,:)=0d0
			T(2,:)=T1
			T(3,:)=T1+T2
			T(4,:)=T2
			Ns=0
			allocate(self%nb(0:nb))
			do n=1,sb
				call gen_grid(a1,a2,self%rsb(n,:),T,tmp)
				do i=1,size(tmp,1)
					bd(Ns+i)%r=tmp(i,:)
					bd(Ns+i)%i=(/Ns+i,Ns+i/)
					bd(Ns+i)%sb=(/n,n/)
					bd(Ns+i)%dr=0d0
					bd(Ns+i)%bdc=1d0
				enddo
				Ns=Ns+size(tmp,1)
			enddo
			allocate(self%nb(0)%bd(Ns),self%nb(0)%st(Ns))
			self%nb(0)%bd(:)=bd(:Ns)
			do i=1,Ns
				allocate(self%nb(0)%st(i)%j(1),self%nb(0)%st(i)%bd(1))
				self%nb(0)%st(i)%j=i
				self%nb(0)%st(i)%bd=i
			enddo

			!do i=2,layer
				!do j=1,size(a3)
					!self%i2r(1+Ns*(i-1):Ns*i,j)=self%i2r(:Ns,j)+a3(j)*(i-1)
				!enddo
			!enddo
			!Ns=Ns*layer

			l1=max(abs(nint(sum(l*a1*T1)/sum(T1**2))),1)
			l2=max(abs(nint(sum(l*a2*T2)/sum(T2**2))),1)
			if(abs(bdc(1))<err) l1=0
			if(abs(bdc(2))<err) l2=0
			do i=1,Ns
				do j=1,Ns
					do m1=-l1,l1
						do m2=-l2,l2
							dr=self%nb(0)%bd(j)%r-self%nb(0)%bd(i)%r+T1*m1+T2*m2
							k=k+1
							sort(k)%val=sum(dr**2)
							sort(k)%i=(/i,j/)
							sort(k)%dr=dr
							sort(k)%bdc=1d0
							if(m1>0) then
								sort(k)%bdc=sort(k)%bdc*bdc(1)
							elseif(m1<0) then
								sort(k)%bdc=sort(k)%bdc*conjg(bdc(1))
							endif
							if(m2>0) then
								sort(k)%bdc=sort(k)%bdc*bdc(2)
							elseif(m2<0) then
								sort(k)%bdc=sort(k)%bdc*conjg(bdc(2))
							endif
						enddo
					enddo
				enddo
			enddo
			call qsort(sort(:k))
			call collect(sort(:k),c)
			do l=1,nb
				st(:,0,:)=0
				k=0
			o:	do m=c(l+1),c(l+2)-1
					do p=c(l+1),m-1
						if(sum(abs(sort(m)%dr(:)+sort(p)%dr(:)))<err.and.sort(m)%i(2)==sort(p)%i(1).and.sort(m)%i(1)==sort(p)%i(2)) then
							cycle o
						endif
					enddo
					k=k+1
					bd(k)%i=sort(m)%i
					bd(k)%sb=(/self%nb(0)%bd(sort(m)%i(1))%sb(1),self%nb(0)%bd(sort(m)%i(2))%sb(1)/)
					bd(k)%r=self%nb(0)%bd(sort(m)%i(1))%r+sort(m)%dr/2d0
					bd(k)%dr=sort(m)%dr
					bd(k)%bdc=sort(m)%bdc
					st(sort(m)%i(1),0,1)=st(sort(m)%i(1),0,1)+1
					st(sort(m)%i(1),st(sort(m)%i(1),0,1),:)=(/sort(m)%i(2),k/)
					st(sort(m)%i(2),0,1)=st(sort(m)%i(2),0,1)+1
					st(sort(m)%i(2),st(sort(m)%i(2),0,1),:)=(/sort(m)%i(1),k/)
				enddo o
				allocate(self%nb(l)%bd(k),self%nb(l)%st(Ns))
				self%nb(l)%bd=bd(:k)
				do i=1,Ns
					allocate(self%nb(l)%st(i)%bd(st(i,0,1)),self%nb(l)%st(i)%j(st(i,0,1)))
					self%nb(l)%st(i)%j=st(i,1:st(i,0,1),1)
					self%nb(l)%st(i)%bd=st(i,1:st(i,0,1),2)
				enddo
			enddo
			deallocate(c)
		end associate
	end subroutine
	subroutine get_brizon(b1,b2,bT)
		real(8) :: b1(:),b2(:)
		real(8), allocatable :: bT(:,:)
		real(8) :: tf(2,2),T(0:24,3),th,dr1(2),dr2(2)
		type(t_mysort) :: st(24)
		integer :: n,i,j
		integer, allocatable :: ist(:)
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
		if(allocated(bT)) deallocate(bT)
		allocate(bT(n,3))
		do i=1,n
			j=mod(i,n)+1
			dr1=T(i,:2)
			dr2=T(j,:2)
			bT(i,:2)=matmul(reshape((/dr2(2),-dr2(1),-dr1(2),dr1(1)/),(/2,2/)),(/sum(dr1**2),sum(dr2**2)/))/(dr1(1)*dr2(2)-dr1(2)*dr2(1))
			bT(i,3)=0d0
		enddo
	end subroutine
	subroutine gen_brizon(self,brizon)
		class(t_latt) :: self
		type(t_brizon) :: brizon
		real(8) :: p1(3),p2(3),k0(3),dr(3)
		real(8), allocatable :: tmp(:,:),T(:,:)
		integer :: n,i,j
		associate(b1=>brizon%b1,b2=>brizon%b2,n1=>self%n1,n2=>self%n2,c1=>self%c1,c2=>self%c2,T1=>self%T1,T2=>self%T2,a1=>self%a1,a2=>self%a2,bdc=>self%bdc)
			k0=0d0
			if(all(abs(bdc(:2))>err)) then
				k0(:2)=matmul(reshape((/T2(2),-T1(2),-T2(1),T1(1)/)/(T1(1)*T2(2)-T1(2)*T2(1)),(/2,2/)),(/theta((/real(bdc(1)),imag(bdc(1))/)),theta((/real(bdc(2)),imag(bdc(2))/))/))
			elseif(all(abs(bdc(:2))<err)) then
				allocate(brizon%k(1,3))
				allocate(brizon%ok(1,3))
				brizon%k=0d0
				brizon%ok=0d0
				return
			elseif(abs(bdc(2))<err) then
				k0(:2)=T1(1:2)
				k0=theta((/real(bdc(1)),imag(bdc(1))/))*k0/sum(T1*k0)
			else
				k0(:2)=T2(1:2)
				k0=theta((/real(bdc(2)),imag(bdc(2))/))*k0/sum(T2*k0)
			endif

			p1=(/-T2(2),T2(1),0d0/)
			p2=(/T1(2),-T1(1),0d0/)
			p1=2d0*pi*p1/sum(T1*n2*p1)
			p2=2d0*pi*p2/sum(T2*n1*p2)

			b1=(/-c2(2),c2(1),0d0/)
			b2=(/c1(2),-c1(1),0d0/)
			b1=2d0*pi*b1/sum(b1*c1)
			b2=2d0*pi*b2/sum(b2*c2)

			call get_brizon(b1,b2,brizon%T)
			call gen_grid(p1,p2,k0,brizon%T,brizon%k)

			p1=(/-a2(2),a2(1),0d0/)
			p2=(/a1(2),-a1(1),0d0/)
			p1=2d0*pi*p1/sum(p1*a1)
			p2=2d0*pi*p2/sum(p2*a2)
			call get_brizon(p1,p2,T)

			call gen_grid(b1,b2,(/0d0,0d0,0d0/),T,tmp)
			allocate(brizon%ok(size(brizon%k,1)*size(tmp,1),3))

			do i=1,size(tmp,1)
				do j=1,size(brizon%k,1)
					brizon%ok(j+(i-1)*size(brizon%k,1),:)=brizon%k(j,:)+tmp(i,:)
					if(.not.is_in(brizon%ok(j+(i-1)*size(brizon%k,1),:),T,dr)) then
						brizon%ok(j+(i-1)*size(brizon%k,1),:2)=brizon%ok(j+(i-1)*size(brizon%k,1),:2)+dr
					endif
				enddo
			enddo

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
	subroutine check_lattice(ut)
		integer :: ut,k,l,m,n

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

		do k=1,size(brizon%ok,1)
			write(ut,"(es13.2$)")brizon%ok(k,:),0d0,0d0,0d0
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

		do k=1,size(latt%nb(0)%bd)
			write(ut,"(es13.2$)")latt%nb(0)%bd(k)%r-latt%nb(0)%bd(k)%dr/2d0,latt%nb(0)%bd(k)%dr
			write(ut,"(i7$)")latt%nb(0)%bd(k)%i(1)
			write(ut,"(x)")
		enddo
		write(ut,"(x/)")

		do l=1,ubound(latt%nb,1)
			do k=1,size(latt%nb(l)%bd)
				write(ut,"(es13.2$)")latt%nb(l)%bd(k)%r-latt%nb(l)%bd(k)%dr/2d0,latt%nb(l)%bd(k)%dr
				write(ut,"(i7$)")k
				write(ut,"(x)")
			enddo
			write(ut,"(x/)")
		enddo
	end subroutine
end module
