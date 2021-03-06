module M_lattice_final
	use M_const
	use M_utility
	implicit none
	private
	type, extends(t_sort) :: t_mysort
		integer :: i(2)
		real(wp) :: dr(3)
		complex(wp) :: bdc
	contains
		procedure :: swap_sort => myswap
	end type
	type ::  t_bd
		integer :: sb(2)
		integer :: sbc(2)
		integer :: i(2)
		real(wp) :: r(3)
		real(wp) :: dr(3)
		complex(wp) :: bdc
	end type
	type ::  t_st
		integer, allocatable :: j(:)
		integer, allocatable :: bd(:)
		integer, allocatable :: dir(:)
	end type
	type t_nb
		type(t_bd), allocatable :: bd(:)
		type(t_st), allocatable :: st(:)
	end type
	real(wp) :: err=1e-6_wp
	type t_latt
		logical :: is_all=.false.
		real(wp) :: T1(3),T2(3),a1(3),a2(3),a3(3),c1(3),c2(3),c3(3)
		integer :: Ns,Nc,Ni,sb,layer=1
		complex(wp) :: bdc(3)=(/1._wp,1._wp,0._wp/)
		type(t_nb), allocatable :: nb(:)
		real(wp), allocatable :: rsb(:,:)
		integer, allocatable :: i2isb(:,:),isb2i(:,:)
	contains
		procedure gen_latt
		procedure gen_brizon
	end type
	type t_brizon
		integer :: nk,nq
		real(wp) :: a1(3),a2(3)
		real(wp) :: c1(3),c2(3)
		real(wp) :: k0(3)
		real(wp), allocatable :: Ta(:,:)
		real(wp), allocatable :: Tc(:,:)
		real(wp), allocatable :: q(:,:)
		real(wp), allocatable :: k(:,:)
	end type
	type(t_latt), public :: latt
	type(t_brizon), public, target :: brizon
	public t_latt,t_brizon,theta,check_lattice,is_in
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
	recursive function is_in(r,T,dr)
		logical :: is_in
		real(wp) :: r(3),T(:,:)
		real(wp), optional :: dr(:)
		real(wp) :: tr(2),EG(3),DT(3),tf(2,2),rerr(3),r_(3)
		integer :: i,j,n
		r_=r
		if(present(dr)) then
			dr=0._wp
		endif
		is_in=.true.
		rerr=0._wp
		n=size(T,1)
		do i=1,n/2
			rerr=rerr+(T(i+1,:)-T(i,:))*err*i/1.34_wp
		enddo
		do 
			do i=1,size(T,1)
				j=mod(i,size(T,1))+1
				EG=T(j,:)-T(i,:)
				DT=(T(mod(j-1+n/2,n)+1,:)-T(j,:))+EG
				tf=reshape((/EG(2),-DT(2),-EG(1),DT(1)/)/(DT(1)*EG(2)-DT(2)*EG(1)),(/2,2/))
				tr=matmul(tf,r_(:2)-(T(j,:)+T(i,:))/2._wp+rerr(:2))
				if(tr(1)<0._wp) then
					is_in=.false.
					if(present(dr)) then
						r_=r_-floor(tr(1))*DT
					else
						return
					endif
				endif
			enddo
			if(present(dr)) then
				if(is_in(r_,T)) then
					dr=r_-r
					exit
				endif
			else
				exit
			endif
		enddo
	end function
	subroutine gen_grid(a1,a2,r0,T,grid)
		real(wp), allocatable :: grid(:,:)
		real(wp) :: r0(:),a1(:),a2(:),T(:,:)
		real(wp) :: r(3),x(size(T,1)),y(size(T,1)),tf(2,2),tmp(10000000,3)
		integer :: i,j,n
		n=1
		tf=reshape((/a2(2),-a1(2),-a2(1),a1(1)/)/(a1(1)*a2(2)-a1(2)*a2(1)),(/2,2/))
		do i=1,size(T,1)
			x(i)=sum(tf(1,:)*(T(i,:2)-r0(:2)))
			y(i)=sum(tf(2,:)*(T(i,:2)-r0(:2)))
		enddo
		do j=nint(minval(y)),nint(maxval(y))
			do i=nint(minval(x)),nint(maxval(x))
				r=a1*i+a2*j+r0
				if(is_in(r,T)) then
					if(i==0.and.j==0) then
						tmp(1,:)=r
					else
						n=n+1
						tmp(n,:)=r
					endif
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
		integer :: i,j,l,l1,l2,n1,n2,m1,m2,n,m,p,k,Ns_
		real(wp) :: r(3),tf(2,2),T(4,3,3),dr(3),dr2,rnb(0:nb)
		real(wp), allocatable :: tmp(:,:),DT(:,:)
		type(t_mysort) :: sort(100)
		type(t_bd), allocatable ::  bd(:,:,:)
		integer, allocatable ::  st(:,:,:)
		integer, allocatable :: c(:)
		associate(T1=>self%T1,T2=>self%T2,a1=>self%a1,a2=>self%a2,a3=>self%a3,c1=>self%c1,c2=>self%c2,Ns=>self%Ns,Nc=>self%Nc,Ni=>self%Ni,sb=>self%sb,layer=>self%layer,bdc=>self%bdc)
			sb=size(self%rsb,1)

			allocate(self%nb(0:nb))

			T(1,:,1)=0._wp
			T(2,:,1)=a1
			T(3,:,1)=a1+a2
			T(4,:,1)=a2
			T(1,:,2)=0._wp
			T(2,:,2)=c1
			T(3,:,2)=c1+c2
			T(4,:,2)=c2
			T(1,:,3)=0._wp
			T(2,:,3)=T1
			T(3,:,3)=T1+T2
			T(4,:,3)=T2

			if(self%is_all) then
				allocate(bd(0:nb,2*nint(abs((T1(1)*T2(2)-T1(2)*T2(1))/(a1(1)*a2(2)-a1(2)*a2(1))))*sb,0:10))
			else
				allocate(bd(0:nb,2*nint(abs((c1(1)*c2(2)-c1(2)*c2(1))/(a1(1)*a2(2)-a1(2)*a2(1))))*sb,0:10))
			endif
			allocate(st(size(bd,2),0:20,2))


			bd(:,1:sb,0)%sb(1)=0
			bd(:,1:sb,0)%sb(2)=0
			do i=1,sb
				k=0
				do j=1,sb
					do m1=-nb,nb
						do m2=-nb,nb
							dr=self%rsb(j,:)-self%rsb(i,:)+a1*m1+a2*m2
							k=k+1
							sort(k)%val=sum(dr**2)
							sort(k)%i=(/i,j/)
							sort(k)%dr=dr
							sort(k)%bdc=1._wp
						enddo
					enddo
				enddo
				call qsort(sort(:k))
				call collect(sort(:k),c)
				do m1=0,nb
				o:	do m2=c(m1+1),c(m1+2)-1
						do p=1,i
							do k=1,bd(m1,p,0)%sb(1)
								if(sum(abs(sort(m2)%dr(:)+bd(m1,p,k)%dr(:)))<err.and.sort(m2)%i(2)==bd(m1,p,k)%i(1).and.sort(m2)%i(1)==bd(m1,p,k)%i(2)) then
									cycle o
								endif
							enddo
						enddo
						bd(m1,i,0)%sb(1)=bd(m1,i,0)%sb(1)+1
						bd(m1,i,bd(m1,i,0)%sb(1))%dr=sort(m2)%dr
						bd(m1,i,bd(m1,i,0)%sb(1))%sb=sort(m2)%i
						bd(m1,i,bd(m1,i,0)%sb(1))%sbc=sort(m2)%i
						bd(m1,i,bd(m1,i,0)%sb(1))%i=sort(m2)%i
						bd(m1,i,bd(m1,i,0)%sb(1))%bdc=sort(m2)%bdc
						bd(m1,i,bd(m1,i,0)%sb(1))%r=self%rsb(i,:)+sort(m2)%dr/2._wp
					enddo o
				enddo
				deallocate(c)
			enddo
			Ns=sb


			call gen_grid(a1,a2,bd(0,1,1)%r,T(:,:,2),DT)
			Nc=size(DT,1)

			m=0
			do i=1,Nc
				do j=1,Ns
					l=j+(i-1)*Ns
					m=m+1
					do n=0,nb
						bd(n,l,0)%sb=bd(n,j,0)%sb
						do k=1,bd(n,j,0)%sb(1)
							!if(bd(0,bd(n,l,k)%sb(2),1)) then
								!exit
							!endif
								
							bd(n,l,k)%dr=bd(n,j,k)%dr
							bd(n,l,k)%sb=bd(n,j,k)%sb
							bd(n,l,k)%sbc(2)=bd(n,j,k)%sb(2)
							bd(n,l,k)%bdc=1._wp
							if(n==0) then
								if(k==1) then
									bd(n,l,k)%r=bd(n,j,k)%r+DT(i,:)
									if(.not.is_in(bd(n,l,k)%r,T(:,:,2),dr)) then
										bd(n,l,k)%r=bd(n,l,k)%r+dr
									endif
								else
									bd(n,l,k)%r=bd(n,l,1)%r
								endif
								bd(n,l,k)%sbc(2)=bd(n,j,k)%sbc(2)+(i-1)*Ns
								r=bd(n,l,k)%r+bd(0,1,1)%r-bd(0,bd(n,l,k)%sb(2),1)%r
								if(.not.is_in(r,T(:,:,1),dr)) then
									r=r+dr
								endif
							else
								bd(n,l,k)%r=bd(0,l,1)%r+bd(n,j,k)%dr/2._wp
								r=bd(0,l,1)%r+bd(n,j,k)%dr+bd(0,1,1)%r-bd(0,bd(n,l,k)%sb(2),1)%r
								if(.not.is_in(r,T(:,:,2),dr)) then
									r=r+dr
								endif
								if(.not.is_in(r,T(:,:,1),dr)) then
									r=r+dr
									do p=1,Nc
										if(sum(abs(dr+DT(p,:)))<err) then
											bd(n,l,k)%sbc(2)=bd(n,l,k)%sbc(2)+(p-1)*Ns
											exit
										endif
									enddo
								endif
							endif
							!if(sum(abs(r-bd(0,bd(n,l,k)%sb(2),1)%r))>err) then
							if(sum(abs(r-self%rsb(1,:)))>err) then
								bd(n,l,k)%sbc(2)=-1
								if(n==0) then
									m=m-1
								endif
							endif
							bd(n,l,k)%sbc(1)=m
							bd(n,l,k)%i=bd(n,l,k)%sbc
						enddo
					enddo
				enddo
			enddo
			Ni=l
			if(any(abs(bdc(:2))<err)) then
				do n=0,nb
					do i=1,Ni
						m2=0
						do k=1,bd(n,i,0)%sb(1)
							if(bd(n,i,k)%sbc(2)>0) then 
								m2=m2+1
								bd(n,i,k)%sbc(2)=bd(0,bd(n,i,k)%sbc(2),1)%sbc(1)
								bd(n,i,k)%i(2)=bd(n,i,k)%sbc(2)
								bd(n,i,m2)=bd(n,i,k)
							endif
						enddo
						bd(n,i,0)%sb(1)=m2
					enddo
				enddo
				m1=0
				do i=1,Ni
					if(bd(0,i,1)%sbc(2)>0) then 
						m1=m1+1
						do n=0,nb
							do k=0,bd(n,i,0)%sb(1)
								bd(n,m1,k)=bd(n,i,k)
							enddo
						enddo
					endif
				enddo
				Ni=m1
			endif


			if(self%is_all) then
				call gen_grid(c1,c2,bd(0,1,1)%r,T(:,:,3),DT)
				allocate(self%i2isb(Ni*size(DT,1),sb),self%isb2i(Nc*size(DT,1),sb))
				Nc=size(DT,1)
				Ns=Ni*Nc
			else
				allocate(self%i2isb(Ni,sb),self%isb2i(Nc,sb))
				Nc=1
				if(allocated(DT)) deallocate(DT)
				allocate(DT(1,3))
				DT=0._wp
				Ns=Ni*nint(abs((T1(1)*T2(2)-T1(2)*T2(1))/(c1(1)*c2(2)-c1(2)*c2(1))))
			endif

			tf=reshape((/T2(2),-T1(2),-T2(1),T1(1)/)/(T1(1)*T2(2)-T1(2)*T2(1)),(/2,2/))

			do i=1,Nc
				do j=1,Ni
					l=j+(i-1)*Ni
					do n=0,nb
						bd(n,l,0)%sb=bd(n,j,0)%sb
						do k=1,bd(n,j,0)%sb(1)
							bd(n,l,k)%dr=bd(n,j,k)%dr
							bd(n,l,k)%sb=bd(n,j,k)%sb
							bd(n,l,k)%sbc=bd(n,j,k)%sbc
							bd(n,l,k)%i(1)=l
							bd(n,l,k)%i(2)=bd(n,j,k)%sbc(2)
							bd(n,l,k)%bdc=1._wp
							if(n==0) then
								if(k==1) then
									bd(n,l,k)%r=bd(n,j,k)%r+DT(i,:)
									if(.not.is_in(bd(n,l,k)%r,T(:,:,3),dr)) then
										bd(n,l,k)%r=bd(n,l,k)%r+dr
									endif
								else
									bd(n,l,k)%r=bd(n,l,1)%r
								endif
								bd(n,l,k)%i(2)=bd(n,j,k)%i(2)+(i-1)*Ni
							else
								bd(n,l,k)%r=bd(0,l,1)%r+bd(n,j,k)%dr/2._wp
								r=bd(0,l,1)%r+bd(n,j,k)%dr
								if(.not.is_in(r,T(:,:,3),dr)) then
									r=r+dr
									dr(:2)=matmul(tf,dr(:2))
									do p=1,2
										if(nint(dr(p))<0) then
											bd(n,l,k)%bdc=bd(n,l,k)%bdc*bdc(p)**(-nint(dr(p)))
										elseif(nint(dr(p))>0) then
											bd(n,l,k)%bdc=bd(n,l,k)%bdc*conjg(bdc(p))**nint(dr(p))
										endif
									enddo
								endif
								r=r+bd(0,1,1)%r-bd(0,bd(n,l,k)%sbc(2),1)%r
								if(.not.is_in(r,T(:,:,3),dr)) then
									r=r+dr
								endif
								if(.not.is_in(r,T(:,:,2),dr)) then
									do p=1,Nc
										if(sum(abs(dr+DT(p,:)))<err) then
											bd(n,l,k)%i(2)=bd(n,l,k)%i(2)+(p-1)*Ni
											exit
										endif
									enddo
								endif
							endif
						enddo
					enddo
				enddo
			enddo


			do n=0,nb
				allocate(self%nb(n)%bd(Nc*sum(bd(n,1:Ni,0)%sb(1))))
				allocate(self%nb(n)%st(Ni*Nc))
				k=0
				do i=1,Ni*Nc
					do j=1,bd(n,i,0)%sb(1)
						k=k+1
						self%nb(n)%bd(k)=bd(n,i,j)
						st(bd(n,i,j)%i(1),j,:)=[bd(n,i,j)%i(2),k]
						if(n/=0.or.(n==0.and.bd(n,i,j)%i(1)/=bd(n,i,j)%i(2))) then
							bd(n,bd(n,i,j)%i(2),0)%sb(2)=bd(n,bd(n,i,j)%i(2),0)%sb(2)+1
							st(bd(n,i,j)%i(2),sum(bd(n,bd(n,i,j)%i(2),0)%sb),:)=[i,k]
						endif
					enddo
				enddo
				do i=1,Ni*Nc
					k=sum(bd(n,i,0)%sb)
					allocate(self%nb(n)%st(i)%bd(k),self%nb(n)%st(i)%j(k),self%nb(n)%st(i)%dir(k))
					self%nb(n)%st(i)%j=st(i,1:k,1)
					self%nb(n)%st(i)%bd=st(i,1:k,2)
					self%nb(n)%st(i)%dir(1:bd(n,i,0)%sb(1))=1
					self%nb(n)%st(i)%dir(bd(n,i,0)%sb(1)+1:k)=-1
				enddo
			enddo
			Nc=Ns/Ni

			self%i2isb(1,:)=0
			self%i2isb(1,self%nb(0)%bd(1)%sb(1))=1
			self%isb2i(1,self%nb(0)%bd(1)%sb(1))=1
			do i=2,size(self%i2isb,1)
				self%i2isb(i,:)=self%i2isb(i-1,:)
				self%i2isb(i,self%nb(0)%bd(i)%sb(1))=self%i2isb(i,self%nb(0)%bd(i)%sb(1))+1
				self%isb2i(self%i2isb(i,self%nb(0)%bd(i)%sb(1)),self%nb(0)%bd(i)%sb(1))=i
			enddo

		end associate
	end subroutine
	subroutine get_brizon(b1,b2,bT)
		real(wp) :: b1(:),b2(:)
		real(wp), allocatable :: bT(:,:)
		real(wp) :: tf(2,2),T(0:24,3),th,dr1(2),dr2(2)
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
				st(n)%dr=(b1*i+b2*j)/2._wp
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
			if(sin(th-theta(dr1)-err)<0._wp) then
				n=n-1
			endif
			if(sin(th-theta(dr2)+err)<0._wp) then
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
			bT(i,3)=0._wp
		enddo
	end subroutine
	subroutine gen_brizon(self,brizon)
		class(t_latt) :: self
		type(t_brizon) :: brizon
		real(wp) :: p1(3),p2(3),dr(3)
		real(wp), allocatable :: tmp(:,:)
		integer :: n,i,j
		associate(ba1=>brizon%a1,ba2=>brizon%a2,bc1=>brizon%c1,bc2=>brizon%c2,c1=>self%c1,c2=>self%c2,T1=>self%T1,T2=>self%T2,a1=>self%a1,a2=>self%a2,bdc=>self%bdc,k0=>brizon%k0)
			p1=(/-T2(2),T2(1),0._wp/)
			p2=(/T1(2),-T1(1),0._wp/)
			p1=2._wp*pi*p1/sum(T1*p1)
			p2=2._wp*pi*p2/sum(T2*p2)

			ba1=(/-a2(2),a2(1),0._wp/)
			ba2=(/a1(2),-a1(1),0._wp/)
			ba1=2._wp*pi*ba1/sum(ba1*a1)
			ba2=2._wp*pi*ba2/sum(ba2*a2)
			call get_brizon(ba1,ba2,brizon%Ta)

			bc1=(/-c2(2),c2(1),0._wp/)
			bc2=(/c1(2),-c1(1),0._wp/)
			bc1=2._wp*pi*bc1/sum(bc1*c1)
			bc2=2._wp*pi*bc2/sum(bc2*c2)
			call get_brizon(bc1,bc2,brizon%Tc)


			k0=0._wp
			if(all(abs(bdc(:2))>err)) then
				k0(:2)=matmul(reshape((/T2(2),-T2(1),-T1(2),T1(1)/)/(T1(1)*T2(2)-T1(2)*T2(1)),(/2,2/)),(/theta((/real(bdc(1)),imag(bdc(1))/)),theta((/real(bdc(2)),imag(bdc(2))/))/))
			elseif(all(abs(bdc(:2))<err)) then
				allocate(brizon%k(1,3))
				allocate(brizon%q(1,3))
				brizon%k=0._wp
				brizon%q=0._wp
				brizon%nk=1
				brizon%nq=1
				return
			elseif(abs(bdc(2))<err) then
				k0(:2)=T1(1:2)
				k0=theta((/real(bdc(1)),imag(bdc(1))/))*k0/sum(T1*k0)
			else
				k0(:2)=T2(1:2)
				k0=theta((/real(bdc(2)),imag(bdc(2))/))*k0/sum(T2*k0)
			endif

			call gen_grid(p1,p2,k0,brizon%Tc,tmp)
			brizon%nk=size(tmp,1)

			if(any(abs(bdc(:2))<err)) then
				allocate(brizon%q(1,3))
				brizon%q=0._wp
				brizon%nq=1
			else
				call gen_grid(bc1,bc2,(/0._wp,0._wp,0._wp/),brizon%Ta,brizon%q)
				brizon%nq=size(brizon%q,1)
			endif
			allocate(brizon%k(brizon%nk*brizon%nq,3))

			do i=1,brizon%nq
				do j=1,brizon%nk
					brizon%k(j+(i-1)*brizon%nk,:)=tmp(j,:)+brizon%q(i,:)
					if(.not.is_in(brizon%k(j+(i-1)*brizon%nk,:),brizon%Ta,dr)) then
						brizon%k(j+(i-1)*brizon%nk,:2)=brizon%k(j+(i-1)*brizon%nk,:2)+dr
					endif
				enddo
			enddo

		end associate
	end subroutine
	function theta(r)
		real(wp) :: r(:),theta,d
		d=sqrt(sum(r**2))
		if(d<1e-10_wp) then
			theta=0._wp
			return
		endif
		theta=acos(r(1)/d)
		if(r(2)<0._wp) then
			theta=2._wp*pi-theta
		endif
	end function
	subroutine check_lattice(ut)
		integer :: ut,k,l,m,n

		write(ut,"(*(es13.2))")0._wp,0._wp,0._wp,brizon%a1
		write(ut,"(*(es13.2))")0._wp,0._wp,0._wp,brizon%a2
		write(ut,"(x)")

		do k=1,size(brizon%Ta,1)+1
			write(ut,"(*(es13.2))")brizon%Ta(mod(k-1,size(brizon%Ta,1))+1,:)
		enddo
		write(ut,"(x)")

		write(ut,"(*(es13.2))")0._wp,0._wp,0._wp,brizon%c1
		write(ut,"(*(es13.2))")0._wp,0._wp,0._wp,brizon%c2
		write(ut,"(x)")

		do k=1,size(brizon%Tc,1)+1
			write(ut,"(*(es13.2))")brizon%Tc(mod(k-1,size(brizon%Tc,1))+1,:)
		enddo
		write(ut,"(x)")


		do k=1,size(brizon%k,1)
			write(ut,"(es13.2$)")brizon%k(k,:)
			write(ut,"(i7)")k
			if(mod(k,brizon%nk)==0) then
				write(ut,"(x)")
			endif
		enddo
		write(ut,"(x)")

		write(ut,"(6es13.2)")0._wp,0._wp,0._wp,latt%a1
		write(ut,"(6es13.2)")0._wp,0._wp,0._wp,latt%a2
		write(ut,"(x)")

		write(ut,"(6es13.2)")0._wp,0._wp,0._wp,latt%c1
		write(ut,"(6es13.2)")0._wp,0._wp,0._wp,latt%c2
		write(ut,"(x)")

		if(latt%is_all) then
			write(ut,"(6es13.2)")0._wp,0._wp,0._wp,latt%T1
			write(ut,"(6es13.2)")0._wp,0._wp,0._wp,latt%T2
			write(ut,"(x)")
		else
			write(ut,"(6es13.2)")0._wp,0._wp,0._wp,0._wp,0._wp,0._wp
			write(ut,"(6es13.2)")0._wp,0._wp,0._wp,0._wp,0._wp,0._wp
			write(ut,"(x)")
		endif

		do k=1,size(latt%nb(0)%bd)
			write(ut,"(es13.2$)")latt%nb(0)%bd(k)%r-latt%nb(0)%bd(k)%dr/2._wp,latt%nb(0)%bd(k)%dr
			write(ut,"(i7$)")latt%nb(0)%bd(k)%sbc(1)
			write(ut,"(x)")
			if(latt%is_all) then
				if(mod(k,latt%Ns/latt%Nc)==0) then
					write(ut,"(x)")
				endif
			elseif(k==size(latt%nb(0)%bd)) then
				write(ut,"(x)")
			endif
		enddo
		write(ut,"(x)")

		do l=1,ubound(latt%nb,1)
			do k=1,size(latt%nb(l)%bd)
				write(ut,"(es13.2$)")latt%nb(l)%bd(k)%r-latt%nb(l)%bd(k)%dr/2._wp,latt%nb(l)%bd(k)%dr,abs(latt%nb(l)%bd(k)%bdc-1._wp)
				write(ut,"(i7$)")k
				write(ut,"(x)")
			enddo
			write(ut,"(x)")
		enddo
		!do l=1,ubound(latt%nb,1)
		!do k=1,size(latt%nb(l)%st)
		!do i=1,size(latt%nb(l)%st(k)%bd)
		!if(abs(latt%nb(l)%bd(latt%nb(l)%st(k)%bd(i))%bdc)>err) then
		!write(ut,"(es13.2$)")latt%nb(l)%bd(latt%nb(l)%st(k)%bd(i))%r-latt%nb(l)%st(k)%dir(i)*latt%nb(l)%bd(latt%nb(l)%st(k)%bd(i))%dr/2._wp,latt%nb(l)%st(k)%dir(i)*latt%nb(l)%bd(latt%nb(l)%st(k)%bd(i))%dr
		!write(ut,"(i7$)")k
		!write(ut,"(x)")
		!endif
		!enddo
		!enddo
		!write(ut,"(x)")
		!enddo
	end subroutine
end module
