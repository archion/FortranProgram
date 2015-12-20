module M_Hamilton_test
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
		integer :: Vnb
		integer :: nb
		integer :: sg   
						! 1:  chemical potainial
						! 2:  singlet pair
						! 3:  spin symmetry charge
						! 4:  spin z
						! 5:  spin asymmetry charge
						! 6:  triplet pair
						! 7:  spin xy
						! 9: n-n jast
						! -: don't do variation
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
	procedure(real(8)), pointer :: free_energy
	procedure(mat_diag_interface), pointer :: mat_diag
	procedure(self_consist_interface), pointer :: self_consist
	interface
		subroutine self_consist_interface(fE)
			real(8), optional :: fE
		end subroutine
		subroutine mat_diag_interface(H,E,info)
			complex(8) :: H(:,:)
			real(8) :: E(:)
			integer, optional :: info
		end subroutine
	end interface
contains
	subroutine default_mat_diag(H,E,info)
		complex(8) :: H(:,:)
		real(8) :: E(:)
		integer, optional :: info
		select case(size(H,1))
		case(2)
			call diag2(H,E)
		case(100:)
			call heevd(H,E,"V")
		case default
			call heev(H,E,"V")
		end select
	end subroutine
	function dwave(i,n,th)
		integer, intent(in) :: i
		integer, optional, intent(in) :: n
		real(8), optional, intent(in) :: th
		real(8) :: dwave
		integer :: o_n
		real(8) :: o_th
		if(present(n)) then
			o_n=n
		else
			o_n=1
		endif
		if(present(th)) then
			o_th=th
		else
			o_th=0d0
		endif
		dwave=cos((theta(latt%bond(o_n)%bd(i)%dr)-o_th)*2d0)
		!write(*,"(es12.4$)")dwave,latt%bond(o_n)%bd(i)%dr
		!read(*,*)
	end function
	subroutine gen_var(sg,nb,val,V,Vnb)
		integer, intent(in) :: sg,nb
		real(8), intent(in), optional :: val(:)
		real(8), intent(in), optional :: V
		integer, intent(in), optional :: Vnb
		type(t_mysort), allocatable :: sort(:)
		integer, allocatable :: collect(:)
		integer :: i,j
		iv(sign(1,sg))=iv(sign(1,sg))+sign(1,sg)
		iv(0)=iv(sign(1,sg))
		var(iv(0))%sg=sg
		var(iv(0))%nb=nb
		if(present(Vnb)) then
			var(iv(0))%Vnb=Vnb
		else
			var(iv(0))%Vnb=nb
		endif
		if(present(val)) then

			allocate(sort(size(val)))
			sort(:)%val=val
			sort(:)%idx=(/(i,i=1,size(sort))/)
			call sort%qsort()
			call sort%collect(collect)
			allocate(var(iv(0))%i2v(size(val)))
			allocate(var(iv(0))%bd(size(val)))
			if(present(V)) then
				if(var(iv(0))%Vnb==var(iv(0))%nb) then
					allocate(var(iv(0))%Vbd(size(val)))
				else
					allocate(var(iv(0))%Vbd(size(latt%bond(var(iv(0))%Vnb)%bd)))
				endif
				var(iv(0))%Vbd=V
			endif
			allocate(var(iv(0))%val(size(collect)-1))
			allocate(var(iv(0))%v2i(size(collect)-1))
			var(iv(0))%n=size(collect)-1

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
			if(present(V)) then
				allocate(var(iv(0))%Vbd(size(latt%bond(var(iv(0))%Vnb)%bd)))
				var(iv(0))%Vbd=V
			endif
			allocate(var(iv(0))%val(1))
			allocate(var(iv(0))%v2i(1))
			allocate(var(iv(0))%v2i(1)%i(size(latt%bond(nb)%bd)))
			var(iv(0))%v2i(1)%i=(/(i,i=1,size(latt%bond(nb)%bd))/)
			var(iv(0))%i2v=1
			var(iv(0))%n=1
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
	subroutine var_shrink()
		integer :: lb
		type(t_var) :: tmp(iv(-1):iv(1))
		tmp(iv(-1):iv(1))=var(iv(-1):iv(1))
		deallocate(var)
		allocate(var(iv(-1):iv(1)))
		var(iv(-1):iv(1))=tmp(iv(-1):iv(1))
		free_energy => default_free_energy
		mat_diag => default_mat_diag
	end subroutine
	subroutine Hamilton(self,H,k,fn,info)
		class(t_var), intent(in) :: self(:)
		complex(8), intent(inout) :: H(:,:)
		real(8), optional, intent(in) :: k(:)
		complex(8), external, optional :: fn
		integer, optional, intent(in) :: info
		complex(8) :: bdc,expk,bd
		integer :: i,j,l1,l2,l3,n
		real(8) :: isc
		! TO-DO: bd is incorrect when nambu space is not used
		if(.not.present(info)) then
			H=0d0
		endif
		expk=1d0
		if(is_sc) then
			isc=-1d0
		else
			isc=1d0
		endif
		associate(Ns=>latt%Ns)
			do l1=1,size(self)
				do l2=1,size(self(l1)%val)
					do l3=1,size(self(l1)%v2i(l2)%i)
						n=self(l1)%v2i(l2)%i(l3)
						i=latt%bond(self(l1)%nb)%bd(n)%i(1)
						j=latt%bond(self(l1)%nb)%bd(n)%i(2)
						bdc=latt%bond(self(l1)%nb)%bd(n)%bdc
						bd=self(l1)%val(l2)*self(l1)%bd(n)
						if(allocated(self(l1)%Vbd)) then
							bd=bd*self(l1)%Vbd(n)
						elseif(self(l1)%nb==0) then
							bd=bd/2d0
						endif
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
							if(self(l1)%nb==self(l1)%Vnb) then
								H(i,j+Ns)=H(i,j+Ns)+conjg(bd)*bdc*expk
								H(i+Ns,j)=H(i+Ns,j)+conjg(bd)*conjg(bdc)*expk
								H(j+Ns,i)=H(j+Ns,i)+bd*conjg(bdc*expk)
								H(j,i+Ns)=H(j,i+Ns)+bd*conjg(conjg(bdc)*expk)
							else
								do n=1,size(latt%neb(i)%nb(self(l1)%Vnb)%neb)
									j=latt%neb(i)%nb(self(l1)%Vnb)%neb(n)
									bd=self(l1)%val(l2)*self(l1)%bd(i)*self(l1)%Vbd(latt%neb(i)%nb(self(l1)%Vnb)%bond(n))
									H(j,j+Ns)=H(j,j+Ns)+conjg(bd)*bdc*expk
									H(j+Ns,j)=H(j+Ns,j)+bd*conjg(bdc*expk)
								enddo
							endif
						case(3,1,4,5)
							! charge and spin z channel
							if(mod(abs(self(l1)%sg),10)==4) then
								isc=-isc*0.5d0
							endif
							if(mod(abs(self(l1)%sg),10)==5) then
								isc=-isc
							endif
							if(self(l1)%nb==self(l1)%Vnb) then
								H(i,j)=H(i,j)+abs(isc)*conjg(bd)*bdc*expk
								H(j,i)=H(j,i)+abs(isc)*bd*conjg(bdc*expk)
								if(spin==2) then
									if(is_sc) then
										H(i+Ns,j+Ns)=H(i+Ns,j+Ns)+isc*bd*bdc*expk
										H(j+Ns,i+Ns)=H(j+Ns,i+Ns)+isc*conjg(bd*bdc*expk)
									else
										H(i+Ns,j+Ns)=H(i+Ns,j+Ns)+isc*conjg(bd)*bdc*expk
										H(j+Ns,i+Ns)=H(j+Ns,i+Ns)+isc*bd*conjg(bdc*expk)
									endif
								endif
							else
								do n=1,size(latt%neb(i)%nb(self(l1)%Vnb)%neb)
									j=latt%neb(i)%nb(self(l1)%Vnb)%neb(n)
									bd=self(l1)%val(l2)*self(l1)%bd(i)*self(l1)%Vbd(latt%neb(i)%nb(self(l1)%Vnb)%bond(n))
									H(j,j)=H(j,j)+abs(isc)*conjg(bd)*bdc*expk
									if(spin==2) then
										if(is_sc) then
											H(j+Ns,j+Ns)=H(j+Ns,j+Ns)+isc*bd*bdc*expk
										else
											H(j+Ns,j+Ns)=H(j+Ns,j+Ns)+isc*conjg(bd)*bdc*expk
										endif
									endif
								enddo
							endif
							if(mod(abs(self(l1)%sg),10)==4) then
								isc=-isc*2d0
							endif
							if(mod(abs(self(l1)%sg),10)==5) then
								isc=-isc
							endif
						case(7)
							! spin xy channel
							if(is_sc==.true.) then
								write(*,*)"is_sc is true while spin xy channel is used"
							endif
							if(self(l1)%nb==self(l1)%Vnb) then
								H(i,i+Ns)=H(i,i+Ns)+0.5d0*bd*bdc*expk
								H(i+Ns,i)=H(i+Ns,i)+0.5d0*conjg(bd*bdc*expk)
								H(j,j+Ns)=H(j,j+Ns)+0.5d0*bd*bdc*expk
								H(j+Ns,j)=H(j+Ns,j)+0.5d0*conjg(bd*bdc*expk)
							else
								do n=1,size(latt%neb(i)%nb(self(l1)%Vnb)%neb)
									j=latt%neb(i)%nb(self(l1)%Vnb)%neb(n)
									bd=self(l1)%val(l2)*self(l1)%bd(i)*self(l1)%Vbd(latt%neb(i)%nb(self(l1)%Vnb)%bond(n))
									H(j,j+Ns)=H(j,j+Ns)+0.5d0*bd*bdc*expk
									H(j+Ns,j)=H(j+Ns,j)+0.5d0*conjg(bd*bdc*expk)
								enddo
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
		integer :: i,j,l1,l2,l3,l,m1,m2,m2p,n
		real(8) :: isc
		D=0d0
		if(is_sc) then
			isc=-1d0
		else
			isc=1d0
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
						bd=self(l1)%bd(n)
						if(allocated(self(l1)%Vbd)) then
							bd=bd*self(l1)%Vbd(n)
						elseif(self(l1)%nb==0) then
							bd=bd/2d0
						endif
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
								select case(mod(abs(self(l1)%sg),10))
								case(2)
									! pair channel
									if(self(l1)%nb==self(l1)%Vnb) then
										tmp=tmp+cH(m1,i)*H(j+Ns,m2)*conjg(bd)*bdc*expk
										tmp=tmp+cH(m1,i+Ns)*H(j,m2)*conjg(bd)*conjg(bdc)*expk
										tmp=tmp+cH(m1,j+Ns)*H(i,m2)*bd*conjg(bdc*expk)
										tmp=tmp+cH(m1,j)*H(i+Ns,m2)*bd*conjg(conjg(bdc)*expk)
									else
										do n=1,size(latt%neb(i)%nb(self(l1)%Vnb)%neb)
											bd=self(l1)%bd(latt%neb(i)%nb(self(l1)%Vnb)%neb(n))*self(l1)%Vbd(latt%neb(i)%nb(self(l1)%Vnb)%bond(n))
											tmp=tmp+cH(m1,i)*H(i+Ns,m2)*conjg(bd)*bdc*expk
											tmp=tmp+cH(m1,i+Ns)*H(i,m2)*bd*conjg(bdc*expk)
										enddo
									endif
								case(3,1,4,5)
									! charge and spin z channel
									if(mod(abs(self(l1)%sg),10)==4) then
										isc=-isc*0.5d0
									endif
									if(mod(abs(self(l1)%sg),10)==5) then
										isc=-isc
									endif
									if(self(l1)%nb==self(l1)%Vnb) then
										tmp=tmp+cH(m1,i)*H(j,m2)*abs(isc)*conjg(bd)*bdc*expk
										tmp=tmp+cH(m1,j)*H(i,m2)*abs(isc)*bd*conjg(bdc*expk)
										if(spin==2) then
											if(is_sc) then
												tmp=tmp+cH(m1,i+Ns)*H(j+Ns,m2)*isc*bd*bdc*expk
												tmp=tmp+cH(m1,j+Ns)*H(i+Ns,m2)*isc*conjg(bd*bdc*expk)
											else
												tmp=tmp+cH(m1,i+Ns)*H(j+Ns,m2)*isc*conjg(bd)*bdc*expk
												tmp=tmp+cH(m1,j+Ns)*H(i+Ns,m2)*isc*bd*conjg(bdc*expk)
											endif
										endif
									else
										do n=1,size(latt%neb(i)%nb(self(l1)%Vnb)%neb)
											bd=self(l1)%bd(latt%neb(i)%nb(self(l1)%Vnb)%neb(n))*self(l1)%Vbd(latt%neb(i)%nb(self(l1)%Vnb)%bond(n))
											tmp=tmp+cH(m1,i)*H(i,m2)*abs(isc)*conjg(bd)*bdc*expk
											if(spin==2) then
												if(is_sc) then
													tmp=tmp+cH(m1,i+Ns)*H(i+Ns,m2)*isc*bd*bdc*expk
												else
													tmp=tmp+cH(m1,i+Ns)*H(i+Ns,m2)*isc*conjg(bd)*bdc*expk
												endif
											endif
										enddo
									endif
									if(mod(abs(self(l1)%sg),10)==4) then
										isc=-isc*2d0
									endif
									if(mod(abs(self(l1)%sg),10)==5) then
										isc=-isc
									endif
								case(7)
									! spin xy channel
									if(self(l1)%nb==self(l1)%Vnb) then
										tmp=tmp+cH(m1,i)*H(i+Ns,m2)*0.5d0*bd*bdc*expk
										tmp=tmp+cH(m1,i+Ns)*H(i,m2)*0.5d0*conjg(bd*bdc*expk)
										tmp=tmp+cH(m1,j)*H(j+Ns,m2)*0.5d0*bd*bdc*expk
										tmp=tmp+cH(m1,j+Ns)*H(j,m2)*0.5d0*conjg(bd*bdc*expk)
									else
										do n=1,size(latt%neb(i)%nb(self(l1)%Vnb)%neb)
											bd=self(l1)%bd(latt%neb(i)%nb(self(l1)%Vnb)%neb(n))*self(l1)%Vbd(latt%neb(i)%nb(self(l1)%Vnb)%bond(n))
											tmp=tmp+cH(m1,i)*H(i+Ns,m2)*0.5d0*bd*bdc*expk
											tmp=tmp+cH(m1,i+Ns)*H(i,m2)*0.5d0*conjg(bd*bdc*expk)
										enddo
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
	subroutine  MF_val(ll,v,fE,info)
		integer :: ll
		real(8) :: v(:)
		real(8) :: fE
		integer :: info ! 1: F=<H>-TS; else F=-TlnZ

		complex(8) :: H(latt%Ns*spin,latt%Ns*spin)
		real(8) :: E(size(H,1)),av(sum(var%n))
		complex(8) :: D(size(H,1),1,size(av)),cH(size(H,1),size(H,2))
		real(8) :: f(size(E)),dE(size(H,1),size(D,3)),tmp,mt,fE1,x(size(av)),px(sum(var(:0)%n))
		integer :: l,lp,l1,l2,l3,m,i,j,k
		av=0d0
		fE=0d0
		!$OMP PARALLEL DO REDUCTION(+:fE,av) PRIVATE(H,E,cH,f,dE,D)
		do k=1,size(brizon%k,1)
			call var%Hamilton(H,brizon%k(k,:))
			call mat_diag(H,E)
			cH=transpose(conjg(H))
			f=1d0/(exp(min(E/Tk,1d2))+1d0)
			!call var(ll:)%dHamilton(H,cH,D(:,1:1,:),brizon%k(k,:))
			call var(:)%dHamilton(H,cH,D(:,1:1,:),brizon%k(k,:))
			dE=real(D(:,1,:))
			if(info==1) then
				fE=fE+entropy(E,f)
			else
				fE=fE+fenergy(E,f)
			endif
			do l=1,size(av)
				av(l)=av(l)+sum(dE(:,l)*f(:))
			enddo
		enddo
		!$OMP END PARALLEL DO
		fE=fE/size(brizon%k,1)
		av=av/size(brizon%k,1)

		fE1=fE/latt%Ns

		!Nambu space
		l=0
		if(is_sc) then
			do l1=lbound(var,1),ubound(var,1)
				do l2=1,size(var(l1)%val)
					l=l+1
					if(var(l1)%nb==0) then
						do l3=1,size(var(l1)%v2i(l2)%i)
							i=var(l1)%v2i(l2)%i(l3)
							do m=1,size(latt%neb(i)%nb(var(l1)%Vnb)%neb)
								j=latt%neb(i)%nb(var(l1)%Vnb)%neb(m)
								tmp=real(var(l1)%bd(j))
								if(allocated(var(l1)%Vbd)) tmp=tmp*var(l1)%Vbd(latt%neb(i)%nb(var(l1)%Vnb)%bond(m))
								if(mod(abs(var(l1)%sg),10)==4) then
									tmp=-tmp*0.5d0
								endif
								if(var(l1)%Vnb==0.and.allocated(var(l1)%Vbd)) then
									tmp=tmp*2d0
								endif
								fE=fE+tmp*var(l1)%val(var(l1)%i2v(j))
								av(l)=av(l)+tmp
							enddo
						enddo
					endif
				enddo
			enddo
		endif

		l=0
		x=0d0
		do l1=lbound(var,1),ubound(var,1)
			do l2=1,size(var(l1)%val)
				l=l+1
				if(var(l1)%sg==1) then
					fE=fE-nf*var(l1)%val(1)*real(sum(var(1)%bd))
					av(l)=av(l)-nf*real(sum(var(l1)%bd))
					x(l)=var(l1)%val(1)
					cycle
				endif
				mt=0d0
				if(var(l1)%nb==var(l1)%Vnb) then
					do l3=1,size(var(l1)%v2i(l2)%i)
						i=var(l1)%v2i(l2)%i(l3)
						if(allocated(var(l1)%Vbd)) then
							tmp=real(conjg(var(l1)%bd(i))*var(l1)%bd(i))*var(l1)%Vbd(i)
							if(var(l1)%nb/=0) then
								tmp=tmp*spin
							endif
							if(l1>0) fE=fE-tmp*var(l1)%val(l2)**2
							mt=mt+tmp*2d0
						else
							i=var(l1)%v2i(l2)%i(1)
							if(var(l1)%nb/=0) then
								if(abs(real(var(l1)%bd(i)))>1d-10) then
									tmp=real(var(l1)%bd(i))*spin*2d0
								else
									tmp=imag(var(l1)%bd(i))*spin*2d0
								endif
							else
								if(abs(real(var(l1)%bd(i)))>1d-10) then
									tmp=real(var(l1)%bd(i))*spin
								else
									tmp=imag(var(l1)%bd(i))*spin
								endif
							endif
							mt=mt+tmp
						endif
					enddo
				else
					do l3=1,size(var(l1)%v2i(l2)%i)
						i=var(l1)%v2i(l2)%i(l3)
						do m=1,size(latt%neb(i)%nb(var(l1)%Vnb)%neb)
							j=latt%neb(i)%nb(var(l1)%Vnb)%neb(m)
							tmp=real(var(l1)%bd(i)*var(l1)%bd(j))*var(l1)%Vbd(latt%neb(i)%nb(var(l1)%Vnb)%bond(m))
							if(l1>0) fE=fE-tmp*var(l1)%val(l2)*var(l1)%val(var(l1)%i2v(j))/2d0
							mt=mt+tmp
						enddo
					enddo
				endif
				if(abs(mt)>1d-10) then
					x(l)=1d0/mt*av(l)
					av(l)=av(l)-mt*var(l1)%val(l2)
				elseif(abs(av(l))>1d-10) then
					write(*,*)"err: mt is zero, may be you could try to distinct more bound"
				endif
			enddo
		enddo
		av=av/latt%Ns
		px=var(:0)%put()
		do l1=lbound(var,1),0
			if(.not.allocated(var(l1)%Vbd)) then
				do l2=1,size(var(l1)%val)
					do l3=1,size(var(l1)%v2i(l2)%i)
						if(abs(var(l1)%val(l2))>1d-17) var(l1)%bd(var(l1)%v2i(l2)%i(l3))=var(l1)%bd(var(l1)%v2i(l2)%i(l3))*var(l1)%val(l2)
					enddo
				enddo
			endif
		enddo
		call var%get(x)
		!fE=fE+MFterm()
		fE=fE/latt%Ns
		if(info==1) then
			fE=fE1+free_energy()
		endif
		call var(:0)%get(px)
		do l1=lbound(var,1),0
			if(.not.allocated(var(l1)%Vbd)) then
				do l2=1,size(var(l1)%val)
					do l3=1,size(var(l1)%v2i(l2)%i)
						if(abs(var(l1)%val(l2))>1d-17) var(l1)%bd(var(l1)%v2i(l2)%i(l3))=var(l1)%bd(var(l1)%v2i(l2)%i(l3))/var(l1)%val(l2)
					enddo
				enddo
			endif
		enddo
		v=av(sum(var(:ll-1)%n)+1:sum(var(:ll-1)%n)+size(v))
	end subroutine
	function fenergy(E,f)
		real(8) :: f(:),E(:)
		integer :: l
		real(8) :: fenergy
		real(8) :: tmp
		tmp=0d0
		do l=1,size(E)
			if((-E(l))>1d2*Tk) then
				tmp=tmp+E(l)
			else
				tmp=tmp-Tk*log(1d0+exp(-E(l)/Tk))
			endif
		enddo
		fenergy=tmp
	end function
	function entropy(E,f)
		real(8) :: f(:),E(:)
		integer :: l
		real(8) :: entropy
		real(8) :: tmp
		tmp=0d0
		do l=1,size(E)
			tmp=tmp-(f(l)*log(max(f(l),1d-17))+(1d0-f(l))*log(max(1d0-f(l),1d-17)))
		enddo
		entropy=-tmp*Tk
	end function
	function MFterm()
		integer :: l1,l2,l3,i,j,m
		real(8) :: MFterm,tmp,fE
		do l1=1,ubound(var,1)
			do l2=1,size(var(l1)%val)
				if(var(l1)%sg==1) then
					fE=fE-nf*var(l1)%val(1)*real(sum(var(1)%bd))
					cycle
				endif
				if(var(l1)%nb==var(l1)%Vnb) then
					do l3=1,size(var(l1)%v2i(l2)%i)
						i=var(l1)%v2i(l2)%i(l3)
						if(allocated(var(l1)%Vbd)) then
							tmp=real(conjg(var(l1)%bd(i))*var(l1)%bd(i))*var(l1)%Vbd(i)
							if(var(l1)%nb/=0) then
								tmp=tmp*spin
							endif
							if(l1>0) fE=fE-tmp*var(l1)%val(l2)**2
						endif
					enddo
				else
					do l3=1,size(var(l1)%v2i(l2)%i)
						i=var(l1)%v2i(l2)%i(l3)
						do m=1,size(latt%neb(i)%nb(var(l1)%Vnb)%neb)
							j=latt%neb(i)%nb(var(l1)%Vnb)%neb(m)
							tmp=real(var(l1)%bd(i)*var(l1)%bd(j))*var(l1)%Vbd(latt%neb(i)%nb(var(l1)%Vnb)%bond(m))
							if(l1>0) fE=fE-tmp*var(l1)%val(l2)*var(l1)%val(var(l1)%i2v(j))/2d0
						enddo
					enddo
				endif
			enddo
		enddo
		MFterm=fE
	end function
	function default_free_energy()
		integer :: l1,n
		real(8) :: default_free_energy
		default_free_energy=0d0
		do l1=lbound(var,1),ubound(var,1)
			if(var(l1)%sg==1) then
				cycle
			endif
			do n=1,size(latt%bond(var(l1)%Vnb)%bd)
				default_free_energy=default_free_energy+Eavg(l1,n)
			enddo
		enddo
		default_free_energy=default_free_energy/latt%Ns
	end function
	function Eavg(l1,n)
		integer :: l1,n
		real(8) :: Eavg
		integer :: i,j,vi,vj,m
		i=latt%bond(var(l1)%Vnb)%bd(n)%i(1)
		j=latt%bond(var(l1)%Vnb)%bd(n)%i(2)
		if(var(l1)%nb==var(l1)%Vnb) then
			vi=var(l1)%i2v(n)
			if(allocated(var(l1)%Vbd)) then
				if(var(l1)%nb/=0) then
					Eavg=real(var(l1)%val(vi)**2*var(l1)%bd(n)*conjg(var(l1)%bd(n))*var(l1)%Vbd(n))*spin
				else
					Eavg=real(var(l1)%val(vi)**2*var(l1)%bd(n)*conjg(var(l1)%bd(n))*var(l1)%Vbd(n))
				endif
			else
				if(mod(abs(var(l1)%sg),10)==4.or.mod(abs(var(l1)%sg),10)==7) then
					if(var(l1)%nb/=0) then
						write(*,*)"spin bond channel case is not considered"
					else
						Eavg=var(l1)%val(vi)*var(l1)%bd(n)
					endif
				else
					m=var(l1)%v2i(vi)%i(1)
					if(var(l1)%nb/=0) then
						if(abs(real(var(l1)%bd(m)))>1d-10) then
							Eavg=var(l1)%val(vi)*real(var(l1)%bd(m))*spin*2d0
						else
							Eavg=var(l1)%val(vi)*imag(var(l1)%bd(m))*spin*2d0
						endif
					else
						if(abs(real(var(l1)%bd(m)))>1d-10) then
							Eavg=var(l1)%val(vi)*real(var(l1)%bd(m))*spin
						else
							Eavg=var(l1)%val(vi)*imag(var(l1)%bd(m))*spin
						endif
					endif
				endif
			endif
		else
			vi=var(l1)%i2v(i)
			vj=var(l1)%i2v(j)
			Eavg=real(var(l1)%bd(i)*var(l1)%bd(j))*var(l1)%Vbd(n)*var(l1)%val(vi)*var(l1)%val(vj)
		endif
	end function
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
		write(ut,"(x/)")
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
			write(ut,"(es17.9$)")k(i,:)
			write(ut,"(es17.9$)")(A(i,j),j=1,size(E)),(A(i,j+size(E)),j=1,size(E))
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
		write(ut,"(x/)")
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
		!write(ut,"(x/)")
	end subroutine
	function find_order(l,is,x,rang,z)
		integer :: l,is
		real(8) :: find_order,x,rang(3)
		real(8) :: z
		real(8) :: pod,order,dx,mdx
		integer :: isg
		dx=rang(3)
		mdx=dx/20d0
		call self_consist()
		pod=sum(abs(var(l)%val(:)))/(size(var(l)%val(:)))-z
		do 
			x=x-dx*sign(1d0,pod)*is
			if(x<rang(1).or.x>rang(2)) then
				if(dx>mdx) then
					x=x+dx*sign(1d0,pod)*is
					dx=dx*0.3333d0
					cycle
				else
					find_order=min(max(x,rang(1)),rang(2))
					exit
				endif
			endif
			call self_consist()
			order=sum(abs(var(l)%val(:)))/(size(var(l)%val))
			call is_cross(pod,order-z,isg)
			if(isg/=0) then
				if(abs(dx)<mdx) then
					find_order=x
					exit
				endif
				dx=dx*0.3333d0
			endif
		enddo
	end function
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
