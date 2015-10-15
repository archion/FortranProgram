module M_Hamilton
	use M_lattice
	use M_const
	use M_utility
	use M_matrix
	implicit none
	type, extends(t_sort), private :: t_mysort
		integer :: idx
	contains
		procedure :: swap_sort => myswap1
	end type
	type t_v2i
		integer, allocatable :: i(:)
	end type
	type t_var
		real(8), allocatable :: val(:)
		type(t_v2i), allocatable :: v2i(:)
		integer, allocatable :: i2v(:)
		complex(8), allocatable :: bd(:)
		real(8), allocatable :: Vbd(:)
		integer :: n
		integer :: Vn
		integer :: nb
		integer :: sg   
						! 1:  chemical potainial
						! 2:  singlet pair
						! 3:  charge
						! 4:  spin
						! 7: n-n jast
						! -: don't do variation
		procedure(update), pointer, nopass :: update
	contains
		procedure :: Hamilton
		procedure :: dHamilton
		procedure :: put
		procedure :: get
	end type
	type(t_var), allocatable :: var(:)
	integer :: spin=2
	logical :: is_sc=.true.
	integer :: iv(-1:1)=(/1,0,0/)
	real(8) :: Tk,nf
	private myswap1
	interface 
		subroutine update()
		end subroutine
	end interface
contains
	function dwave(i)
		integer, intent(in) :: i
		real(8) :: dwave
		if(nint(latt%bond(1)%bd(i)%dr(2))==0) then
			dwave=1d0
		else
			dwave=-1d0
		endif
	end function
	subroutine gen_var(sg,nb,val,V,Vn)
		integer, intent(in) :: sg,nb
		real(8), intent(in), optional :: val(:)
		real(8), intent(in), optional :: V
		integer, intent(in), optional :: Vn
		type(t_mysort), allocatable :: sort(:)
		integer, allocatable :: collect(:)
		integer :: i,j
		iv(sign(1,sg))=iv(sign(1,sg))+sign(1,sg)
		iv(0)=iv(sign(1,sg))
		var(iv(0))%sg=sg
		var(iv(0))%nb=nb
		if(present(Vn)) then
			var(iv(0))%Vn=Vn
		else
			var(iv(0))%Vn=nb
		endif
		if(present(val)) then

			allocate(sort(size(val)))
			sort(:)%val=val
			sort(:)%idx=(/(i,i=1,size(sort))/)
			call sort%qsort()
			call sort%collect(collect)
			allocate(var(iv(0))%i2v(size(val)))
			allocate(var(iv(0))%bd(size(val)))
			if(present(Vn)) then
				allocate(var(iv(0))%Vbd(size(latt%bond(Vn)%bd)))
			else
				allocate(var(iv(0))%Vbd(size(val)))
			endif
			allocate(var(iv(0))%val(size(collect)-1))
			allocate(var(iv(0))%v2i(size(collect)-1))
			var(iv(0))%n=size(collect)-1
			if(present(V)) then
				var(iv(0))%Vbd=V
			else
				var(iv(0))%Vbd=1d0
			endif

			do i=1,size(collect)-1
				allocate(var(iv(0))%v2i(i)%i(collect(i+1)-collect(i)))
				var(iv(0))%v2i(i)%i=sort(collect(i):collect(i+1)-1)%idx
				do j=collect(i),collect(i+1)-1
					var(iv(0))%i2v(sort(j)%idx)=i
				enddo
			enddo
		else
			allocate(var(iv(0))%i2v(size(latt%bond(nb)%bd)))
			allocate(var(iv(0))%bd(size(latt%bond(nb)%bd)))
			allocate(var(iv(0))%Vbd(size(latt%bond(var(iv(0))%Vn)%bd)))
			allocate(var(iv(0))%val(1))
			allocate(var(iv(0))%v2i(1))
			allocate(var(iv(0))%v2i(1)%i(size(latt%bond(nb)%bd)))
			var(iv(0))%v2i(1)%i=(/(i,i=1,size(latt%bond(nb)%bd))/)
			var(iv(0))%i2v=1
			var(iv(0))%n=1
			if(present(V)) then
				var(iv(0))%Vbd=V
			else
				var(iv(0))%Vbd=1d0
			endif
		endif
	end subroutine
	subroutine get(self,val)
		class(t_var), intent(inout) :: self(:)
		real(8), intent(in) :: val(:)
		integer :: n,l1,l2
		n=0
		do l1=1,size(self)
			do l2=1,size(self(l1)%val)
				n=n+1
				if(n>size(val)) then
					return
				endif
				self(l1)%val(l2)=val(n)
			enddo
		enddo
	end subroutine
	function put(self)
		class(t_var), intent(in) :: self(:)
		real(8) :: put(sum(self%n))
		integer :: n,l1,l2
		n=0
		do l1=1,size(self)
			do l2=1,size(self(l1)%val)
				n=n+1
				put(n)=self(l1)%val(l2)
			enddo
		enddo
	end function
	subroutine default_update()
	end subroutine
	subroutine var_shrink()
		integer :: lb
		type(t_var) :: tmp(iv(-1):iv(1))
		tmp(iv(-1):iv(1))=var(iv(-1):iv(1))
		deallocate(var)
		allocate(var(iv(-1):iv(1)))
		var(iv(-1):iv(1))=tmp(iv(-1):iv(1))
		var(iv(-1))%update => default_update
	end subroutine
	subroutine Hamilton(self,H,k,fn,info)
		class(t_var), intent(in) :: self(:)
		complex(8), intent(inout) :: H(:,:)
		real(8), optional, intent(in) :: k(:)
		complex(8), external, optional :: fn
		integer, optional, intent(in) :: info
		complex(8) :: bdc,expk,bd
		integer :: i,j,l1,l2,l3,n,isc
		! TO-DO: bd is incorrect when nambu space is not used
		if(.not.present(info)) then
			H=0d0
		endif
		expk=1d0
		if(is_sc) then
			isc=-1
		else
			isc=1
		endif
		associate(Ns=>latt%Ns)
			do l1=1,size(self)
				do l2=1,size(self(l1)%val)
					do l3=1,size(self(l1)%v2i(l2)%i)
						n=self(l1)%v2i(l2)%i(l3)
						i=latt%bond(self(l1)%nb)%bd(n)%i(1)
						j=latt%bond(self(l1)%nb)%bd(n)%i(2)
						bdc=latt%bond(self(l1)%nb)%bd(n)%bdc
						bd=self(l1)%val(l2)*self(l1)%bd(n)*self(l1)%Vbd(n)
						if(present(k)) then
							expk=exp(-img*sum(k*latt%bond(self(l1)%nb)%bd(n)%dr))
						endif
						if(present(fn)) then
							bd=bd*fn(latt%bond(self(l1)%nb)%bd(n)%dr)
						endif
						select case(mod(abs(self(l1)%sg),10))
						case(2)
							! pair channel
							if(.not.is_sc) then
								write(*,*)"is_sc is not setted as true"
								stop
							endif
							if(present(fn)) then
								cycle
							endif
							if(self(l1)%nb/=0) then
								H(i,j+Ns)=H(i,j+Ns)+bd*bdc*expk
								H(i+Ns,j)=H(i+Ns,j)+bd*conjg(bdc)*expk
								H(j+Ns,i)=H(j+Ns,i)+conjg(bd*bdc*expk)
								H(j,i+Ns)=H(j,i+Ns)+conjg(bd*conjg(bdc)*expk)
							else
								do n=1,size(latt%neb(i)%nb(self(l1)%Vn)%neb)
									j=latt%neb(i)%nb(self(l1)%Vn)%neb(n)
									bd=self(l1)%val(l2)*self(l1)%bd(i)*self(l1)%Vbd(latt%neb(i)%nb(self(l1)%Vn)%bond(n))
									H(j,j+Ns)=H(j,j+Ns)+bd*bdc*expk
									H(j+Ns,j)=H(j+Ns,j)+bd*conjg(bdc)*expk
								enddo
							endif
						case(3,1,4)
							! charge and spin channel
							if(mod(abs(self(l1)%sg),10)==4) then
								isc=-isc
							endif
							if(self(l1)%nb/=0) then
								H(i,j)=H(i,j)+conjg(bd)*bdc*expk
								H(j,i)=H(j,i)+conjg(conjg(bd)*bdc*expk)
								if(spin==2) then
									H(i+Ns,j+Ns)=H(i+Ns,j+Ns)+isc*bd*bdc*expk
									H(j+Ns,i+Ns)=H(j+Ns,i+Ns)+isc*conjg(bd*bdc*expk)
								endif
							else
								do n=1,size(latt%neb(i)%nb(self(l1)%Vn)%neb)
									j=latt%neb(i)%nb(self(l1)%Vn)%neb(n)
									bd=self(l1)%val(l2)*self(l1)%bd(i)*self(l1)%Vbd(latt%neb(i)%nb(self(l1)%Vn)%bond(n))
									H(j,j)=H(j,j)+conjg(bd)*bdc*expk
									if(spin==2) then
										H(j+Ns,j+Ns)=H(j+Ns,j+Ns)+isc*bd*bdc*expk
									endif
								enddo
							endif
							if(mod(abs(self(l1)%sg),10)==4) then
								isc=-isc
							endif
						end select
					enddo
				enddo
			enddo
		end associate
	end subroutine
	subroutine dHamilton(self,H,cH,D,k)
		class(t_var), intent(in) :: self(:)
		complex(8), intent(in) :: H(:,:),cH(:,:)
		real(8), optional, intent(in) :: k(:)
		complex(8), intent(inout) :: D(:,:,:)
		complex(8) :: bdc,expk,bd,tmp
		integer :: i,j,l1,l2,l3,l,m1,m2,m2p,n,isc
		D=0d0
		if(is_sc) then
			isc=-1
		else
			isc=1
		endif
		expk=1d0
		l=0
		associate(Ns=>latt%Ns)
			do l1=1,size(self)
				do l2=1,size(self(l1)%val)
					l=l+1
					if(l>size(D,3)) then
						return
					endif
					do l3=1,size(self(l1)%v2i(l2)%i)
						n=self(l1)%v2i(l2)%i(l3)
						i=latt%bond(self(l1)%nb)%bd(n)%i(1)
						j=latt%bond(self(l1)%nb)%bd(n)%i(2)
						bdc=latt%bond(self(l1)%nb)%bd(n)%bdc
						bd=self(l1)%bd(n)*self(l1)%Vbd(n)
						if(present(k)) then
							expk=exp(-img*sum(k*latt%bond(self(l1)%nb)%bd(n)%dr))
						endif
						do m2p=1,size(D,2)
							do m1=1,size(D,1)
								if(size(D,2)==1) then
									m2=m1 
								else
									m2=m2p
								endif
								tmp=0d0
								select case(mod(self(l1)%sg,10))
								case(2)
									! pair channel
									if(self(l1)%nb/=0) then
										tmp=tmp+cH(m1,i)*H(j+Ns,m2)*bd*bdc*expk
										tmp=tmp+cH(m1,i+Ns)*H(j,m2)*bd*conjg(bdc)*expk
										tmp=tmp+cH(m1,j+Ns)*H(i,m2)*conjg(bd*bdc*expk)
										tmp=tmp+cH(m1,j)*H(i+Ns,m2)*conjg(bd*conjg(bdc)*expk)
									else
										do n=1,size(latt%neb(i)%nb(self(l1)%Vn)%neb)
											!j=latt%neb(i)%nb(self(l1)%Vn)%neb(n)
											bd=self(l1)%bd(latt%neb(i)%nb(self(l1)%Vn)%neb(n))*self(l1)%Vbd(latt%neb(i)%nb(self(l1)%Vn)%bond(n))
											tmp=tmp+cH(m1,j)*H(j+Ns,m2)*bd*bdc*expk
											tmp=tmp+cH(m1,j+Ns)*H(j,m2)*bd*conjg(bdc)*expk
										enddo
									endif
								case(3,1,4)
									! charge and spin channel
									if(mod(self(l1)%sg,10)==4) then
										isc=-isc
									endif
									if(self(l1)%nb/=0) then
										tmp=tmp+cH(m1,i)*H(j,m2)*conjg(bd)*bdc*expk
										tmp=tmp+cH(m1,j)*H(i,m2)*conjg(conjg(bd)*bdc*expk)
										if(spin==2) then
											tmp=tmp+cH(m1,i+Ns)*H(j+Ns,m2)*isc*bd*bdc*expk
											tmp=tmp+cH(m1,j+Ns)*H(i+Ns,m2)*isc*conjg(bd*bdc*expk)
										endif
									else
										do n=1,size(latt%neb(i)%nb(self(l1)%Vn)%neb)
											!j=latt%neb(i)%nb(self(l1)%Vn)%neb(n)
											bd=self(l1)%bd(latt%neb(i)%nb(self(l1)%Vn)%neb(n))*self(l1)%Vbd(latt%neb(i)%nb(self(l1)%Vn)%bond(n))
											tmp=tmp+cH(m1,j)*H(j,m2)*conjg(bd)*bdc*expk
											if(spin==2) then
												tmp=tmp+cH(m1,j+Ns)*H(j+Ns,m2)*isc*bd*bdc*expk
											endif
										enddo
									endif
									if(mod(self(l1)%sg,10)==4) then
										isc=-isc
									endif
								end select
								D(m1,m2p,l)=D(m1,m2p,l)+tmp
							enddo
						enddo
					enddo
				enddo
			enddo
		end associate
	end subroutine
	subroutine MF_var(n,x,v,info)
		integer, intent(in) :: n
		integer, intent(in) :: info
		real(8), intent(inout) :: x(n),v(n)
		complex(8) :: H(latt%Ns*spin,latt%Ns*spin)
		real(8) :: E(size(H,1))
		complex(8) :: D(size(H,1),1,n),cH(size(H,1),size(H,2))
		real(8) :: f(size(E)),err,dE(size(H,1),n),fE,tmp,mt
		integer :: l,lp,l1,l2,l3,m,i,j,k
		common fE,err
		v=0d0
		fE=0d0
		call var(1:)%get(x)
		if(size(x)>1) then
			call var%update()
		endif
		!$OMP PARALLEL DO REDUCTION(+:fE,v) PRIVATE(H,E,cH,f,dE,D)
		do k=1,size(brizon%k,1)
			call var%Hamilton(H,brizon%k(k,:))
			call mat_diag(H,E)
			cH=transpose(conjg(H))
			f=1d0/(exp(E/Tk)+1d0)
			call var(1:)%dHamilton(H,cH,D(:,1:1,:),brizon%k(k,:))
			dE=real(D(:,1,:))
			do l=1,size(E)
				if((-E(l))>1d2*Tk) then
					fE=fE+E(l)
				else
					fE=fE-Tk*log(1d0+exp(-E(l)/Tk))
				endif
			enddo
			do l=1,n
				v(l)=v(l)+sum(dE(:,l)*f(:))
			enddo
		enddo
		!$OMP END PARALLEL DO
		fE=fE/size(brizon%k,1)
		v=v/size(brizon%k,1)

		if(var(1)%sg==1) then
			fE=fE+nf*latt%Ns*var(1)%val(1)
			v(1)=v(1)+nf*latt%Ns
		endif

		!Nambu space
		l=0
		if(is_sc) then
			do l1=1,size(var(1:))
				do l2=1,size(var(l1)%val)
					l=l+1
					if(var(l1)%nb==0) then
						do l3=1,size(var(l1)%v2i(l2)%i)
							i=var(l1)%v2i(l2)%i(l3)
							do m=1,size(latt%neb(i)%nb(var(l1)%Vn)%neb)
								j=latt%neb(i)%nb(var(l1)%Vn)%neb(m)
								tmp=real(var(l1)%bd(j))*var(l1)%Vbd(latt%neb(i)%nb(var(l1)%Vn)%bond(m))
								if(mod(abs(var(l1)%sg),10)==4) then
									tmp=-tmp
								endif
								fE=fE+tmp*var(l1)%val(var(l1)%i2v(j))
								if(l<=n) v(l)=v(l)+tmp
							enddo
						enddo
					endif
				enddo
			enddo
		endif

		l=0
		do l1=1,size(var(1:))
			do l2=1,size(var(l1)%val)
				l=l+1
				mt=0d0
				if(var(l1)%nb==var(l1)%Vn) then
					if(var(l1)%sg/=1) then
						do l3=1,size(var(l1)%v2i(l2)%i)
							i=var(l1)%v2i(l2)%i(l3)
							tmp=abs2(var(l1)%bd(i))*var(l1)%Vbd(i)
							if(var(l1)%nb/=0) then
								tmp=tmp*spin
							endif
							fE=fE-tmp*var(l1)%val(l2)**2
							if(l<=n) then
								mt=mt+tmp*2d0
							endif
						enddo
					endif
				else
					do l3=1,size(var(l1)%v2i(l2)%i)
						i=var(l1)%v2i(l2)%i(l3)
						do m=1,size(latt%neb(i)%nb(var(l1)%Vn)%neb)
							j=latt%neb(i)%nb(var(l1)%Vn)%neb(m)
							tmp=real(var(l1)%bd(i)*conjg(var(l1)%bd(j)))*var(l1)%Vbd(latt%neb(i)%nb(var(l1)%Vn)%bond(m))
							fE=fE-tmp*var(l1)%val(l2)*var(l1)%val(var(l1)%i2v(j))/2d0
							if(l<=n) then
								mt=mt+tmp
							endif
						enddo
					enddo
				endif
				if(l<=n) then
					tmp=(v(l)-mt*x(l))/size(var(l1)%v2i(l2)%i)
					if(info<0.and.var(l1)%sg/=1) then
						x(l)=1d0/mt*v(l)
					endif
					v(l)=tmp
				endif
			enddo
		enddo
		fE=fE/latt%Ns
		err=sum(abs(v))/size(v)

		!if(size(x)==1) then
			!write(*,"(es12.4$)")fE
			!write(*,"(es12.4$)")err
			!write(*,"(es12.4$)")x
			!write(*,"(es12.4$)")v
			!write(*,"(x)")
		!endif
	end subroutine
	!function free_energy()
		!integer :: l,l1,l2,l3,m,i,j
		!real(8) :: free_energy,tmp
		!free_energy=free_energy+nf*latt%Ns*var(1)%val(1)
		!!Nambu space
		!l=0
		!if(is_sc) then
			!do l1=1,size(var(1:))
				!do l2=1,size(var(l1)%val)
					!l=l+1
					!if(var(l1)%nb==0) then
						!do l3=1,size(var(l1)%v2i(l2)%i)
							!i=var(l1)%v2i(l2)%i(l3)
							!do m=1,size(latt%neb(i)%nb(var(l1)%Vn)%neb)
								!j=latt%neb(i)%nb(var(l1)%Vn)%neb(m)
								!tmp=real(var(l1)%bd(j))*var(l1)%Vbd(latt%neb(i)%nb(var(l1)%Vn)%bond(m))
								!if(mod(abs(var(l1)%sg),10)==4) then
									!tmp=-tmp
								!endif
								!free_energy=free_energy+tmp*var(l1)%val(var(l1)%i2v(j))
							!enddo
						!enddo
					!endif
				!enddo
			!enddo
		!endif

		!l=0
		!do l1=1,size(var(1:))
			!do l2=1,size(var(l1)%val)
				!l=l+1
				!if(var(l1)%nb==var(l1)%Vn) then
					!if(var(l1)%sg/=1) then
						!do l3=1,size(var(l1)%v2i(l2)%i)
							!i=var(l1)%v2i(l2)%i(l3)
							!tmp=abs2(var(l1)%bd(i))*var(l1)%Vbd(i)
							!if(var(l1)%nb/=0) then
								!tmp=tmp*spin
							!endif
							!free_energy=free_energy-tmp*var(l1)%val(l2)**2
						!enddo
					!endif
				!else
					!do l3=1,size(var(l1)%v2i(l2)%i)
						!i=var(l1)%v2i(l2)%i(l3)
						!do m=1,size(latt%neb(i)%nb(var(l1)%Vn)%neb)
							!j=latt%neb(i)%nb(var(l1)%Vn)%neb(m)
							!tmp=real(var(l1)%bd(i)*conjg(var(l1)%bd(j)))*var(l1)%Vbd(latt%neb(i)%nb(var(l1)%Vn)%bond(m))
							!free_energy=free_energy-tmp*var(l1)%val(l2)*var(l1)%val(var(l1)%i2v(j))/2d0
						!enddo
					!enddo
				!endif
			!enddo
		!enddo
	!end function
	function Green(gm,k,omg)
		real(8) :: k(:),gm,omg
		complex(8) :: H(latt%Ns*spin,latt%Ns*spin),Green
		real(8) :: E(size(H,1))
		integer :: n,i,j
		call var%Hamilton(H,k)
		call mat_diag(H,E)
		Green=0d0
		do i=1,latt%Ns
			do j=1,latt%Ns
				Green=Green+sum(H(i,:)*conjg(H(j,:))/(omg-E+img*gm))/latt%Ns
			enddo
		enddo
	end function
	subroutine LDOS(ut,gm,i,omg,m)
		integer :: ut,i,m
		real(8) :: gm,omg(2)
		complex(8) :: H(latt%Ns*spin,latt%Ns*spin)
		real(8) :: E(size(H,1)),domg,D(m)
		integer :: k,l
		domg=(omg(2)-omg(1))/m
		D=0d0
		!$OMP PARALLEL DO PRIVATE(H,E) REDUCTION(+:D)
		do k=1,size(brizon%k,1)
			call Hamilton(var,H,brizon%k(k,:))
			call mat_diag(H,E)
			do l=1,m
				D(l)=D(l)+imag(sum(H(i,:)*conjg(H(i,:))/(omg(1)+domg*l-E+img*gm)))
			enddo
		enddo
		!$OMP END PARALLEL DO
		do l=1,m
			write(ut,"(es17.9$)")omg(1)+domg*l,D(l)/(size(brizon%k,1))
			write(ut,"(x)")
		enddo
	end subroutine
	subroutine EDC(ut,gm,k,omg,m)
		integer :: ut,m
		real(8) :: gm,omg(2),k(:)
		complex(8) :: H(latt%Ns*spin,latt%Ns*spin)
		real(8) :: E(size(H,1)),domg,A(m)
		integer :: i,j,l
		call var%Hamilton(H,k)
		call mat_diag(H,E)
		domg=(omg(2)-omg(1))/m
		A=0d0
		!$OMP PARALLEL DO
		do l=1,m
			do i=1,latt%Ns
				do j=1,latt%Ns
					A(l)=A(l)+imag(sum(H(i,:)*conjg(H(j,:))/(omg(1)+domg*l-E+img*gm)))
				enddo
			enddo
		enddo
		!$OMP END PARALLEL DO
		do l=1,m
			write(ut,"(2es17.9)")omg(1)+domg*l,A(l)
			write(ut,"(x)")
		enddo
	end subroutine
	subroutine DOS(ut,gm,omg,m,peak)
		integer :: ut,m
		real(8) :: gm,omg(:)
		real(8), allocatable, optional :: peak(:)
		complex(8) :: H(latt%Ns*spin,latt%Ns*spin)
		real(8) :: E(size(H,1)),D(m),x,domg
		integer :: k,i,l
		integer, allocatable :: ipeak(:)
		domg=(omg(2)-omg(1))/m
		D=0d0
		!$OMP PARALLEL DO PRIVATE(H,E) REDUCTION(+:D)
		do k=1,size(brizon%k,1)
			call var%Hamilton(H,brizon%k(k,:))
			call mat_diag(H,E)
			do l=1,m
				do i=1,latt%Ns
					D(l)=D(l)-imag(sum(H(i,:)*conjg(H(i,:))/(omg(1)+domg*l-E+img*gm)))
				enddo
			enddo
		enddo
		!$OMP END PARALLEL DO
		do l=1,m
			write(ut,"(es17.9$)")omg(1)+domg*l,D(l)/(latt%Ns*size(brizon%k,1))
			write(ut,"(x)")
		enddo
		if(present(peak)) then
			if(allocated(peak)) then
				deallocate(peak)
			endif
			call find_peak(D,ipeak)
			allocate(peak(size(ipeak)))
			do i=1,size(ipeak)
				peak(i)=omg(1)+domg*ipeak(i)
			enddo
		endif
	end subroutine
	subroutine fermis(ut,gm,k,omg)
		integer :: ut
		real(8) :: gm,omg,k(:,:)
		complex(8) :: A(size(k,1)),G(latt%Ns,latt%Ns),H(latt%Ns*spin,latt%Ns*spin)
		real(8) :: q(3),E(size(H,1))
		integer :: i,j,l,m,nk
		A=0d0
		nk=size(brizon%k,1)
		!$OMP PARALLEL DO PRIVATE(H,E,G,q)
		do m=1,nk
			call var%Hamilton(H,brizon%k(m,:))
			call mat_diag(H,E)
			do i=1,latt%Ns
				do j=1,latt%Ns
					G(i,j)=sum(H(i,:)*conjg(H(j,:))/(omg-E+img*gm))
				enddo
			enddo
			do l=0,size(k,1)/nk-1
				q=k(l*nk+m,:)-brizon%k(m,:)
				do i=1,latt%Ns
					do j=1,latt%Ns
						A(l*nk+m)=A(l*nk+m)+G(i,j)*exp(img*sum(q*(latt%i2r(i,:)-latt%i2r(j,:))))/latt%Ns
					enddo
				enddo
			enddo
		enddo
		!$OMP END PARALLEL DO
		do i=1,size(k,1)
			write(ut,"(es17.9$)")k(i,:),A(i)/(latt%Ns*size(brizon%k,1))
			write(ut,"(x)")
		enddo
		write(ut,"(x/)")
	end subroutine
	subroutine energy(ut,k)
		integer :: ut
		real(8) :: k(:,:)
		complex(8) :: H(latt%Ns*spin,latt%Ns*spin),G(latt%Ns,latt%Ns,size(H,1)*2)
		real(8) :: E(size(H,1)),A(size(k,1),size(H,1)*2),q(3)
		integer :: nk,m,l,i,j
		A=0d0
		nk=size(brizon%k,1)
		!$OMP PARALLEL DO PRIVATE(H,E,G,q)
		do m=1,nk
			call var%Hamilton(H,brizon%k(m,:))
			call mat_diag(H,E)
			do i=1,latt%Ns
				do j=1,latt%Ns
					G(i,j,:)=H(i,:)*conjg(H(j,:))
				enddo
			enddo
			do l=0,size(k,1)/nk-1
				q=k(l*nk+m,:)-brizon%k(m,:)
				do i=1,latt%Ns
					do j=1,latt%Ns
						A(l*nk+m,size(E)+1:)=A(l*nk+m,size(E)+1:)+real(G(i,j,:)*exp(img*sum(q*(latt%i2r(i,:)-latt%i2r(j,:)))))
					enddo
				enddo
				A(l*nk+m,:size(E))=E
			enddo
		enddo
		!$OMP END PARALLEL DO
		do i=1,size(k,1)
			do j=1,size(E)
				write(ut,"(es17.9$)")k(i,:),A(i,j),A(i,j+size(E))
				write(ut,"(x)")
			enddo
			write(ut,"(x)")
		enddo
	end subroutine
	subroutine band_e(ut,gm,ki,kf,n,omg,m)
		integer :: ut,n,m
		real(8) :: ki(3),kf(3),gm,omg(2)
		complex(8) :: H(latt%Ns*spin,latt%Ns*spin)
		real(8) :: A(n,m),dk(3),domg,E(size(H,1))
		integer :: l1,l2,i,j
		domg=(omg(2)-omg(1))/m
		dk=(kf-ki)/n
		A=0d0
		!$OMP PARALLEL DO PRIVATE(H,E) REDUCTION(+:A)
		do l1=1,n
			call var%Hamilton(H,ki+dk*l1)
			call mat_diag(H,E)
			do l2=1,m
				do i=1,latt%Ns
					do j=1,latt%Ns
						A(l1,l2)=A(l1,l2)+imag(sum(H(i,:)*conjg(H(j,:))/(omg(1)+domg*l2-E+img*gm)))
					enddo
				enddo
			enddo
		enddo
		!$OMP END PARALLEL DO
		do l1=1,n
			do l2=1,m
				write(ut,"(es17.9$)")ki+dk*l1,omg(1)+domg*l2,A(l1,l2)
				write(ut,"(x)")
			enddo
			write(ut,"(x)")
		enddo
	end subroutine
	subroutine band(ut,ki,kf,n)
		integer :: ut,n,i,l,j,m
		real(8) :: ki(3),kf(3),k(3),dk(3)
		complex(8) :: H(latt%Ns*spin,latt%Ns*spin)
		real(8) :: E(size(H,1)),G
		dk=(kf-ki)/n
		do m=0,n-1
			k=ki+dk*m
			call var%Hamilton(H,k)
			call mat_diag(H,E)
			do l=1,size(E)
				G=0d0
				do i=1,latt%Ns
					do j=1,latt%Ns
						G=G+real(H(i,l)*conjg(H(j,l)))
					enddo
				enddo
				write(ut,"(es17.9$)")k,E(l),G
				write(ut,"(x)")
			enddo
			write(ut,"(x)")
		enddo
	end subroutine
	subroutine myswap1(self,a)
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
			deallocate(tmp)
		end select
	end subroutine
end module
