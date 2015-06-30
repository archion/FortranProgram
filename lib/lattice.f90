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
		real(8) :: dir(3)
		complex(8) :: bdc
	end type
	type t_bd
		type(t_bond), allocatable :: bd(:)
	end type

	real(8), public :: T1(3),T2(3),a1(3),a2(3),a3(3)
	integer, public :: Ns,layer
	complex(8), public :: bdc(3)
	type(t_neb), allocatable, public :: neb(:)
	type(t_bd), allocatable, public :: bond(:)
	real(8), allocatable, public :: sub(:,:)
	real(8), allocatable, public :: i2r(:,:)
	real(8) :: err=1d-8
	public gen_latt,gen_neb,gen_bond
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
		end select
	end subroutine
	function is_in(r,T1,T2,sg,bc)
		logical :: is_in
		real(8) :: r(3),k1,k2,bx,by,rx,ry,T1(3),T2(3)
		integer, optional :: sg(2)
		complex(8), optional :: bc
		k1=T1(2)/T1(1)
		k2=T2(2)/T2(1)
		if(isnan(k1)) then
			k1=1d100
		endif
		if(isnan(k2)) then
			k2=1d100
		endif
		if(present(sg)) then
			bc=1d0
			sg=0
		endif
		bx=-(T1(2)/k2-T1(1))-err
		by=T2(2)-k1*T2(1)-err
		ry=r(2)-k1*r(1)
		rx=-(r(2)/k2-r(1))
		if(rx>=-err.and.rx<bx.and.ry>=-err.and.ry<by) then
			is_in=.true.
		else
			is_in=.false.
			if(present(sg)) then
				if(rx<-err.and.abs(bdc(2))>err) then
					sg(1)=1
					bc=bc*conjg(bdc(2))
				elseif(rx>=bx.and.abs(bdc(2))>err) then
					sg(1)=-1
					bc=bc*bdc(2)
				else
					sg(1)=0
				endif
				if(ry<-err.and.abs(bdc(1))>err) then
					sg(2)=1
					bc=bc*conjg(bdc(1))
				elseif(ry>=by.and.abs(bdc(1))>err) then
					sg(2)=-1
					bc=bc*bdc(1)
				else
					sg(2)=0
				endif
			endif
		endif
	end function
	subroutine gen_latt()
		integer :: i,j,n,nsub
		real(8) :: x(2),y(2),r(3),tmp(1000,3)
		x(1)=minval((/T1(1),T2(1),T1(1)+T2(1),0d0/))
		x(2)=maxval((/T1(1),T2(1),T1(1)+T2(1),0d0/))
		y(1)=minval((/T1(2),T2(2),T1(2)+T2(2),0d0/))
		y(2)=maxval((/T1(2),T2(2),T1(2)+T2(2),0d0/))
		n=0
		do j=int(y(1)/a2(2)),int(y(2)/a2(2))
			do i=int((x(1)-a2(1)*j)/a1(1)),int((x(2)-a2(1)*j)/a1(1))
				r=a1*i+a2*j
				if(is_in(r,T1,T2)) then
					n=n+1
					tmp(n,:)=r
				endif
			enddo
		enddo
		Ns=n
		nsub=size(sub,1)
		allocate(i2r(Ns*size(sub,1)*layer,3))
		i2r(1:Ns,:)=tmp
		do i=2,nsub
			do j=1,size(sub(i,:))
				i2r(1+Ns*(i-1):Ns*i,j)=i2r(:Ns,j)+sub(i,j)
			enddo
		enddo
		Ns=Ns*nsub
		do i=2,layer
			do j=1,size(a3)
				i2r(1+Ns*(i-1):Ns*i,j)=i2r(:Ns,j)+a3(j)*(i-1)
			enddo
		enddo
		Ns=Ns*layer
	end subroutine
	subroutine gen_neb()
		type(t_mysort) :: tmp(Ns)
		real(8) :: p,dr(3)
		integer :: i,j,k,n,sg(2),nbond(100)
		integer, allocatable :: t(:)
		nbond=0
		allocate(neb(Ns))
		do i=1,Ns/layer
			do j=1,Ns/layer
				dr=i2r(j,:)-i2r(i,:)

				if(.not.is_in(dr+T1/2d0+T2/2d0,T1,T2,sg,tmp(j)%bdc)) then
					dr=dr+T1*sg(1)+T2*sg(2)
				endif
				tmp(j)%val=sum(dr**2)
				tmp(j)%idx=j
				tmp(j)%dr=dr
			enddo
			call tmp(:Ns/layer)%qsort()
			call tmp(:Ns/layer)%collect(t)
			do k=1,layer
				allocate(neb(i+Ns/layer*(k-1))%nb(0:size(t)-2))
			enddo
			do j=1,size(t)-1
				allocate(neb(i)%nb(j-1)%neb(t(j+1)-t(j)))
				allocate(neb(i)%nb(j-1)%bond(t(j+1)-t(j)))
				allocate(neb(i)%nb(j-1)%bdc(t(j+1)-t(j)))
				allocate(neb(i)%nb(j-1)%dr(t(j+1)-t(j),3))
				neb(i)%nb(j-1)%neb(:)=tmp(t(j):t(j+1)-1)%idx
				neb(i)%nb(j-1)%bdc(:)=tmp(t(j):t(j+1)-1)%bdc
				neb(i)%nb(j-1)%dr(:,1)=tmp(t(j):t(j+1)-1)%dr(1)
				neb(i)%nb(j-1)%dr(:,2)=tmp(t(j):t(j+1)-1)%dr(2)
				neb(i)%nb(j-1)%dr(:,3)=tmp(t(j):t(j+1)-1)%dr(3)
				do k=2,layer
					allocate(neb(i+Ns/layer*(k-1))%nb(j-1)%neb(t(j+1)-t(j)))
					allocate(neb(i+Ns/layer*(k-1))%nb(j-1)%bond(t(j+1)-t(j)))
					allocate(neb(i+Ns/layer*(k-1))%nb(j-1)%bdc(t(j+1)-t(j)))
					allocate(neb(i+Ns/layer*(k-1))%nb(j-1)%dr(t(j+1)-t(j),3))
					neb(i+Ns/layer*(k-1))%nb(j-1)%neb(:)=neb(i)%nb(j-1)%neb(:)+Ns/layer*(k-1)
					neb(i+Ns/layer*(k-1))%nb(j-1)%bdc(:)=neb(i)%nb(j-1)%bdc(:)
					neb(i+Ns/layer*(k-1))%nb(j-1)%dr(:,1)=neb(i)%nb(j-1)%dr(:,1)
					neb(i+Ns/layer*(k-1))%nb(j-1)%dr(:,2)=neb(i)%nb(j-1)%dr(:,2)
					neb(i+Ns/layer*(k-1))%nb(j-1)%dr(:,3)=neb(i)%nb(j-1)%dr(:,3)
				enddo
			enddo
			deallocate(t)
		enddo
	end subroutine
	subroutine gen_bond(l)
		integer :: l,i,j,k,n,m,p
		type(t_bond) :: tmp(Ns*4)
		allocate(bond(0:l))
		do k=0,l
			n=0
			do i=1,Ns
				if(size(neb(i)%nb)<k+1) then
					cycle
				endif
				do m=1,size(neb(i)%nb(k)%neb)
					j=neb(i)%nb(k)%neb(m)
					if(j<i) then
						cycle
					endif
					n=n+1
					tmp(n)%i=(/i,j/)
					tmp(n)%r=i2r(i,:)+neb(i)%nb(k)%dr(m,:)/2d0
					tmp(n)%dir=neb(i)%nb(k)%dr(m,:)
					tmp(n)%bdc=neb(i)%nb(k)%bdc(m)
					neb(i)%nb(k)%bond(m)=n
					do p=1,size(neb(j)%nb(k)%neb)
						if(neb(j)%nb(k)%neb(p)==i) then
							neb(j)%nb(k)%bond(p)=n
							exit
						endif
					enddo
				enddo
			enddo
			allocate(bond(k)%bd(n))
			bond(k)%bd=tmp(:n)
		enddo
	end subroutine
end module
