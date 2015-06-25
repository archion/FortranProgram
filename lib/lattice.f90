module M_lattice
	use M_const
	use M_utility
	implicit none
	private
	type, extends(t_sort) :: t_mysort
		integer :: idx
		real(8) :: r(2)
		complex(8) :: bdc
	contains
		procedure :: swap_sort => myswap
	end type
	type t_nb
		integer, allocatable :: neb(:)
		complex(8), allocatable :: bdc(:)
		real(8), allocatable :: dir(:,:)
	end type
	type t_neb
		type(t_nb), allocatable :: nb(:)
	end type
	type ::  t_bond
		integer :: i(2)
		real(8) :: r(2)
		real(8) :: dir(2)
		complex(8) :: bdc
	end type
	type t_bd
		type(t_bond), allocatable :: bd(:)
	end type

	real(8) :: T1(2),T2(2),a1(2),a2(2)
	integer :: Ns
	complex(8) :: bdc(2)
	type(t_neb), allocatable :: neb(:)
	type(t_bd), allocatable :: bond(:)
	real(8), allocatable :: sub(:,:)
	real(8), allocatable :: i2r(:,:)
	real(8) :: err=1d-8
	public sub,bond,neb,i2r,T1,T2,a1,a2,Ns,bdc,gen_latt,gen_neb,gen_bond
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
			tmp%r=a%r
			a%r=self%r
			self%r=tmp%r
			tmp%bdc=a%bdc
			a%bdc=self%bdc
			self%bdc=tmp%bdc
		end select
	end subroutine
	function is_in(r,T1,T2,sg)
		logical :: is_in
		integer, optional :: sg(2)
		real(8) :: r(2),k1,k2,bx,by,rx,ry,T1(2),T2(2)
		k1=T1(2)/T1(1)
		k2=T2(2)/T2(1)
		if(isnan(k1)) then
			k1=1d100
		endif
		if(isnan(k2)) then
			k2=1d100
		endif
		bx=-(T1(2)/k2-T1(1))-err
		by=T2(2)-k1*T2(1)-err
		ry=r(2)-k1*r(1)
		rx=-(r(2)/k2-r(1))
		if(rx>=-err.and.rx<bx.and.ry>=-err.and.ry<by) then
			is_in=.true.
			if(present(sg)) then
				sg=0
			endif
		else
			is_in=.false.
			if(present(sg)) then
				if(rx<-err.and.abs(bdc(2))>err) then
					sg(1)=1
				elseif(rx>=bx.and.abs(bdc(2))>err) then
					sg(1)=-1
				else
					sg(1)=0
				endif
				if(ry<-err.and.abs(bdc(1))>err) then
					sg(2)=1
				elseif(ry>=by.and.abs(bdc(1))>err) then
					sg(2)=-1
				else
					sg(2)=0
				endif
			endif
		endif
	end function
	subroutine gen_latt()
		integer :: i,n,nsub
		real(8) :: ka1,ka2,xrang(2),yrang(2),r(2),tmp(3*(T1(1)**2+T1(2)**2),2)
		i=1
		n=1
		ka1=a1(2)/a1(1)
		ka2=a2(2)/a2(1)
		if(isnan(ka1)) then
			ka1=1d100
		endif
		if(isnan(ka2)) then
			ka2=1d100
		endif
		yrang=(/min(T1(2),0d0),T2(2)+max(T1(2),0d0)/)
		xrang=(/min(T2(1),0d0),T1(1)+max(T2(1),0d0)/)
		xrang=xrang+(/min(0d0,-(yrang(2)-yrang(1))/ka2),-min(0d0,(yrang(2)-yrang(1))/ka2)/)
		!write(*,*)xrang,bx
		!write(*,*)yrang,by
		r=(/xrang(1),yrang(1)/)
		do 
			do
				if((r(1)-xrang(2))>-err) then
					r(1)=xrang(1)
					exit
				endif
				if(is_in(r,T1,T2)) then
					tmp(i,:)=r
					i=i+1
				endif
				r=r+a1
			enddo
			r=r+a2
			xrang=xrang+a2(1)
			if((r(2)-yrang(2))>-err) then
				exit
			endif
		enddo
		Ns=i-1
		nsub=size(sub,1)
		allocate(i2r(Ns*nsub,2))
		do i=1,nsub
			i2r(1+Ns*(i-1):Ns*i,1)=tmp(:Ns,1)+sub(i,1)
			i2r(1+Ns*(i-1):Ns*i,2)=tmp(:Ns,2)+sub(i,2)
		enddo
		Ns=Ns*nsub
	end subroutine
	subroutine gen_neb()
		type(t_mysort) :: tmp(Ns)
		real(8) :: p,dr(2)
		integer :: i,j,k,n,sg(2),t(Ns)
		allocate(neb(Ns))
		do i=1,Ns
			do j=1,Ns
				dr=i2r(j,:)-i2r(i,:)

				dr=dr+T1/2d0+T2/2d0

				if(.not.is_in(dr,T1,T2,sg)) then
					dr=dr+T1*sg(1)+T2*sg(2)
				endif
				dr=dr-T1/2d0-T2/2d0

				tmp(j)%val=sum(dr**2)
				tmp(j)%idx=j
				tmp(j)%r=dr
				tmp(j)%bdc=1d0
				select case(sg(1))
				case(-1)
					tmp(j)%bdc=tmp(j)%bdc*bdc(2)
				case(1)
					tmp(j)%bdc=tmp(j)%bdc*conjg(bdc(2))
				end select
				select case(sg(2))
				case(-1)
					tmp(j)%bdc=tmp(j)%bdc*bdc(1)
				case(1)
					tmp(j)%bdc=tmp(j)%bdc*conjg(bdc(1))
				end select
			enddo
			call tmp%qsort()
			n=1
			t(1)=1
			do j=2,Ns
				if((tmp(j)%val-sum(tmp(j-1)%r**2))>err) then
					n=n+1
					t(n)=j
					call tmp(t(n-1):t(n)-1)%qsort()
				endif
				tmp(j)%val=acos(tmp(j)%r(1)/sqrt(tmp(j)%val))
				if(tmp(j)%r(2)<0d0) then
					tmp(j)%val=2d0*pi-tmp(j)%val
				endif
			enddo
			allocate(neb(i)%nb(n-1))
			do j=1,n-1
				allocate(neb(i)%nb(j)%neb(t(j+1)-t(j)))
				allocate(neb(i)%nb(j)%bdc(t(j+1)-t(j)))
				allocate(neb(i)%nb(j)%dir(t(j+1)-t(j),2))
				neb(i)%nb(j)%neb(:)=tmp(t(j):t(j+1)-1)%idx
				neb(i)%nb(j)%bdc(:)=tmp(t(j):t(j+1)-1)%bdc
				neb(i)%nb(j)%dir(:,1)=tmp(t(j):t(j+1)-1)%r(1)
				neb(i)%nb(j)%dir(:,2)=tmp(t(j):t(j+1)-1)%r(2)
			enddo
		enddo
	end subroutine
	subroutine gen_bond(l)
		integer :: l,i,j,k,n
		type(t_bond) :: tmp(Ns*4)
		allocate(bond(l))
		do j=1,l
			n=0
			do i=1,Ns
				if(size(neb(i)%nb)<j+1) then
					cycle
				endif
				do k=1,size(neb(i)%nb(j+1)%neb)
					if(neb(i)%nb(j+1)%neb(k)<=i) then
						cycle
					endif
					n=n+1
					tmp(n)%i=(/i,neb(i)%nb(j+1)%neb(k)/)
					tmp(n)%r=i2r(i,:)+neb(i)%nb(j+1)%dir(k,:)/2d0
					tmp(n)%dir=neb(i)%nb(j+1)%dir(k,:)
					tmp(n)%bdc=neb(i)%nb(j+1)%bdc(k)
				enddo
			enddo
			allocate(bond(j)%bd(n))
			bond(j)%bd=tmp(:n)
		enddo
	end subroutine
	subroutine test()
		integer :: i,j,n
		open(unit=101,file="../data/test.dat")
		do i=1,Ns
			write(101,"(2es12.3,i10)")i2r(i,:),i
		enddo
		write(101,"(x/)")
		write(101,"(2es12.3)")0d0,0d0
		write(101,"(2es12.3)")T1
		write(101,"(2es12.3)")T1+T2
		write(101,"(2es12.3)")T2
		write(101,"(2es12.3)")0d0,0d0
		write(101,"(x/)")
		do i=1,size(bond(1)%bd)
			write(101,"(2es12.3,i6)")bond(1)%bd(i)%r,i
		enddo
		do
			write(*,*)"enter:"
			read(*,*)i,n
			write(*,"(2es8.1$)")neb(i)%nb(n)%bdc(:)
			write(*,"(x)")
			write(*,"(i16$)")neb(i)%nb(n)%neb(:)
			write(*,"(x)")
		enddo
	end subroutine
end module
