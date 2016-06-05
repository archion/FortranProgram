module M_Hamilton_test1
	use M_lattice_test1
	use M_const
	use M_utility
	use M_matrix
	use omp_lib
	implicit none
	type, extends(t_sort), private :: t_mysort
		integer :: idx
	contains
		procedure :: swap_sort => myswap1
	end type
	type t_var
		real(8), allocatable :: val(:)
		complex(8), allocatable :: bd(:)
		integer, allocatable :: bd2v(:)
		integer :: nb
		logical :: is_V
		real(8) :: V
		integer :: Vnb
		integer :: sg   
						! 1:  chemical potainial
						! 2:  singlet pair
						! 3:  spin symmetry charge
						! 4:  spin z
						! 5:  spin asymmetry charge
						! 6:  triplet pair
						! 7:  spin xy
						! 9: n-n jast
						!10: s-s jast
						! -: don't do variation
		integer :: n
	!contains
		!procedure :: Hamilton
		!procedure :: dHamilton
		!procedure :: put
		!procedure :: get
	end type
	type(t_var), allocatable :: var(:)
	integer :: spin=2
	logical :: is_sc=.true.
	integer :: iv(-1:1)=(/1,0,0/)
	real(8) :: Tk,nf
	real(8), pointer :: brizon_k(:,:)
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
	function dwave(i,nb,th)
		integer, intent(in) :: i
		integer, optional, intent(in) :: nb
		real(8), optional, intent(in) :: th
		real(8) :: dwave
		integer :: o_nb
		real(8) :: o_th
		if(present(nb)) then
			o_nb=nb
		else
			o_nb=1
		endif
		if(present(th)) then
			o_th=th
		else
			o_th=0d0
		endif
		dwave=cos((theta(latt%nb(o_nb)%bd(i)%dr)-o_th)*2d0)
	end function
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
	subroutine gen_var(sg,nb,sb,val,V,Vnb)
		integer, intent(in) :: sg,nb
		real(8), intent(in), optional :: val(:)
		real(8), intent(in), optional :: V
		integer, intent(in), optional :: Vnb,sb(:)
		type(t_mysort), allocatable :: sort(:)
		integer, allocatable :: coll(:)
		integer :: i,j
		iv(sign(1,sg))=iv(sign(1,sg))+sign(1,sg)
		iv(0)=iv(sign(1,sg))
		var(iv(0))%sg=sg
		var(iv(0))%nb=nb
		if(present(V)) then
			var(iv(0))%V=V
			var(iv(0))%is_V=.true.
		else
			var(iv(0))%V=1d0
			var(iv(0))%is_V=.false.
		endif
		if(present(Vnb)) then
			var(iv(0))%Vnb=Vnb
		else
			var(iv(0))%Vnb=nb
		endif
		if(present(val)) then
			allocate(sort(size(val)))
			sort(:)%val=val
			sort(:)%idx=(/(i,i=1,size(sort))/)
			call qsort(sort)
			call collect(sort,coll)

			allocate(var(iv(0))%bd2v(size(val)))
			allocate(var(iv(0))%bd(size(val)))
			allocate(var(iv(0))%val(size(coll)-1))
			var(iv(0))%n=size(coll)-1
			do i=1,size(coll)-1
				do j=coll(i),coll(i+1)-1
					var(iv(0))%bd2v(sort(j)%idx)=i
				enddo
			enddo
		else
			allocate(var(iv(0))%bd2v(size(latt%nb(nb)%bd)))
			allocate(var(iv(0))%bd(size(var(iv(0))%bd2v)))
			allocate(var(iv(0))%val(1))
			var(iv(0))%bd2v=1
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
		mat_diag => default_mat_diag
		free_energy => default_free_energy
		brizon_k => brizon%k
	end subroutine
	subroutine Hamilton(self,H,k,fn,info)
		class(t_var), intent(in) :: self(:)
		complex(8), intent(inout) :: H(:,:)
		real(8), optional, intent(in) :: k(:)
		complex(8), external, optional :: fn
		integer, optional, intent(in) :: info
		complex(8) :: bdc,expk,bd
		integer :: i,j,l,n,nb,ii,Ns
		real(8) :: isc,dr(3)
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
		Ns=latt%Ns
		do l=1,size(self)
			do n=1,size(latt%nb(self(l)%Vnb)%bd)
				i=latt%nb(self(l)%Vnb)%bd(n)%i(1)
				j=latt%nb(self(l)%Vnb)%bd(n)%i(2)
				dr=latt%nb(self(l)%Vnb)%bd(n)%dr
				if(self(l)%bd2v(n)<0) cycle
				if(self(l)%is_V) then
					bdc=1d0
				else
					bdc=latt%nb(self(l)%Vnb)%bd(n)%bdc
				endif
				if(.not.self(l)%is_V.and.self(l)%nb==0) then
					bd=0.5d0
				else
					bd=1d0
				endif
				if(present(k)) then
					expk=exp(-img*sum(k*dr))
				endif
				if(present(fn)) then
					bd=bd*fn(dr)
				endif
				select case(mod(abs(self(l)%sg),10))
				case(2)
					! pair channel
					if(.not.is_sc) then
						write(*,*)"is_sc is not setted as true"
						stop
					endif
					if(present(fn)) then
						cycle
					endif
					if(self(l)%nb==self(l)%Vnb) then
						bd=bd*self(l)%val(self(l)%bd2v(n))*self(l)%bd(n)*self(l)%V
						H(i,j+Ns)=H(i,j+Ns)+conjg(bd)*bdc*expk
						H(i+Ns,j)=H(i+Ns,j)+conjg(bd)*conjg(bdc)*expk
						H(j+Ns,i)=H(j+Ns,i)+bd*conjg(bdc*expk)
						H(j,i+Ns)=H(j,i+Ns)+bd*conjg(conjg(bdc)*expk)
					else
						write(*,*)"error: not considered"
					endif
				case(3,1,4,5)
					! charge and spin z channel
					if(mod(abs(self(l)%sg),10)==4) then
						isc=-0.5d0*isc
						if(self(l)%nb/=0) then
							write(*,*)"error: not considered"
						endif
					elseif(mod(abs(self(l)%sg),10)==5) then
						isc=-isc
					endif
					if(self(l)%nb==self(l)%Vnb) then
						bd=bd*self(l)%val(self(l)%bd2v(n))*self(l)%bd(n)*self(l)%V
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
						if(present(fn)) then
							cycle
						endif
						bdc=1d0
						do ii=1,2
							bd=self(l)%val(self(l)%bd2v(i))*self(l)%bd(i)*self(l)%V
							H(j,j)=H(j,j)+abs(isc)*conjg(bd)*bdc
							if(spin==2) then
								if(is_sc) then
									H(j+Ns,j+Ns)=H(j+Ns,j+Ns)+isc*bd*bdc
								else
									H(j+Ns,j+Ns)=H(j+Ns,j+Ns)+isc*conjg(bd)*bdc
								endif
							endif
							call swap(i,j)
						enddo
					endif
					if(mod(abs(self(l)%sg),10)==4) then
						isc=-2d0*isc
					elseif(mod(abs(self(l)%sg),10)==5) then
						isc=-isc
					endif
				case(7)
					! spin xy channel
					if(is_sc==.true.) then
						write(*,*)"is_sc is true while spin xy channel is used"
					endif
					if(self(l)%nb==self(l)%Vnb) then
						if(self(l)%nb/=0) then
							write(*,*)"error: not considered"
						endif
						bd=bd*self(l)%val(self(l)%bd2v(n))*self(l)%bd(n)*self(l)%V
						H(i,i+Ns)=H(i,i+Ns)+0.5d0*bd*bdc*expk
						H(i+Ns,i)=H(i+Ns,i)+0.5d0*conjg(bd*bdc*expk)
						H(j,j+Ns)=H(j,j+Ns)+0.5d0*bd*bdc*expk
						H(j+Ns,j)=H(j+Ns,j)+0.5d0*conjg(bd*bdc*expk)
					else
						if(present(fn)) then
							cycle
						endif
						bdc=1d0
						do ii=1,2
							bd=self(l)%val(self(l)%bd2v(i))*self(l)%bd(i)*self(l)%V
							H(j,j+Ns)=H(j,j+Ns)+0.5d0*bd*bdc
							H(j+Ns,j)=H(j+Ns,j)+0.5d0*conjg(bd*bdc)
							call swap(i,j)
						enddo
					endif
				end select
			enddo
		enddo
	end subroutine
	subroutine dHamilton(self,H,cH,D,k)
		class(t_var), intent(in) :: self(:)
		complex(8), intent(in) :: H(:,:),cH(:,:)
		real(8), optional, intent(in) :: k(:)
		complex(8), intent(inout) :: D(:,:,:)
		complex(8) :: bdc,expk,bd,tmp
		integer :: i,j,l,lv,m1,m2,m2p,n,ii,Ns
		real(8) :: isc,dr(3)
		D=0d0
		Ns=latt%Ns
		do l=1,size(self)
			!$omp parallel do reduction(+:d) private(i,j,dr,bdc,expk,m2,tmp,bd,lv,isc) schedule(runtime) if (.not.omp_in_parallel())
			do n=1,size(latt%nb(self(l)%Vnb)%bd)
				expk=1d0
				if(is_sc) then
					isc=-1d0
				else
					isc=1d0
				endif
				i=latt%nb(self(l)%Vnb)%bd(n)%i(1)
				j=latt%nb(self(l)%Vnb)%bd(n)%i(2)
				dr=latt%nb(self(l)%Vnb)%bd(n)%dr
				if(self(l)%bd2v(n)<0) cycle
				if(self(l)%is_V) then
					bdc=1d0
				else
					bdc=latt%nb(self(l)%Vnb)%bd(n)%bdc
				endif
				if(present(k)) then
					expk=exp(-img*sum(k*dr))
				endif
				do m2p=1,size(D,2)
					do m1=1,size(D,1)
						if(size(D,2)==1) then
							m2=m1 
						else
							m2=m2p
						endif
						tmp=0d0
						select case(mod(abs(self(l)%sg),10))
						case(2)
							! pair channel
							if(self(l)%nb==self(l)%Vnb) then
								bd=self(l)%bd(n)*self(l)%V
								tmp=tmp+cH(m1,i)*H(j+Ns,m2)*conjg(bd)*bdc*expk
								tmp=tmp+cH(m1,i+Ns)*H(j,m2)*conjg(bd)*conjg(bdc)*expk
								tmp=tmp+cH(m1,j+Ns)*H(i,m2)*bd*conjg(bdc*expk)
								tmp=tmp+cH(m1,j)*H(i+Ns,m2)*bd*conjg(conjg(bdc)*expk)
								lv=sum(self(:l)%n)-self(l)%n+self(l)%bd2v(n)
								D(m1,m2p,lv)=D(m1,m2p,lv)+tmp
							else
								write(*,*)"error: not considered"
							endif
						case(3,1,4,5)
							! charge and spin z channel
							if(mod(abs(self(l)%sg),10)==4) then
								isc=-isc*0.5d0
							elseif(mod(abs(self(l)%sg),10)==5) then
								isc=-isc
							endif
							if(self(l)%nb==self(l)%Vnb) then
								bd=self(l)%bd(n)*self(l)%V
								if(.not.self(l)%is_V.and.self(l)%nb==0) then
									bd=bd*0.5d0
								endif
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
								lv=sum(self(:l)%n)-self(l)%n+self(l)%bd2v(n)
								D(m1,m2p,lv)=D(m1,m2p,lv)+tmp
							else
								do ii=1,2
									tmp=0d0
									bdc=1d0
									bd=self(l)%bd(j)*self(l)%V
									tmp=tmp+cH(m1,i)*H(i,m2)*abs(isc)*conjg(bd)*bdc
									if(spin==2) then
										if(is_sc) then
											tmp=tmp+cH(m1,i+Ns)*H(i+Ns,m2)*isc*bd*bdc
										else
											tmp=tmp+cH(m1,i+Ns)*H(i+Ns,m2)*isc*conjg(bd)*bdc
										endif
									endif
									lv=sum(self(:l)%n)-self(l)%n+self(l)%bd2v(i)
									D(m1,m2p,lv)=D(m1,m2p,lv)+tmp
									call swap(i,j)
								enddo
							endif
							if(mod(abs(self(l)%sg),10)==4) then
								isc=-isc*2d0
							elseif(mod(abs(self(l)%sg),10)==5) then
								isc=-isc
							endif
						case(7)
							! spin xy channel
							if(self(l)%nb==self(l)%Vnb) then
								bd=self(l)%bd(n)*self(l)%V
								if(.not.self(l)%is_V.and.self(l)%nb==0) then
									bd=bd*0.5d0
								endif
								tmp=tmp+cH(m1,i)*H(i+Ns,m2)*0.5d0*bd*bdc*expk
								tmp=tmp+cH(m1,i+Ns)*H(i,m2)*0.5d0*conjg(bd*bdc*expk)
								tmp=tmp+cH(m1,j)*H(j+Ns,m2)*0.5d0*bd*bdc*expk
								tmp=tmp+cH(m1,j+Ns)*H(j,m2)*0.5d0*conjg(bd*bdc*expk)
								lv=sum(self(:l)%n)-self(l)%n+self(l)%bd2v(n)
								D(m1,m2p,lv)=D(m1,m2p,lv)+tmp
							else
								do ii=1,2
									tmp=0d0
									bdc=1d0
									bd=self(l)%bd(j)*self(l)%V
									tmp=tmp+cH(m1,i)*H(i+Ns,m2)*0.5d0*bd*bdc
									tmp=tmp+cH(m1,i+Ns)*H(i,m2)*0.5d0*conjg(bd*bdc)
									lv=sum(self(:l)%n)-self(l)%n+self(l)%bd2v(i)
									D(m1,m2p,lv)=D(m1,m2p,lv)+tmp
									call swap(i,j)
								enddo
							endif
						end select
					enddo
				enddo
			enddo
			!$omp end parallel do
		enddo
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
	subroutine  MF_val(ll,v,fE,info)
		integer :: ll
		real(8) :: v(:)
		real(8) :: fE
		integer :: info ! 1: F=<H>-TS; else F=-TlnZ

		complex(8) :: H(latt%Ns*spin,latt%Ns*spin)
		real(8) :: E(size(H,1)),av(sum(var%n))
		complex(8) :: D(size(H,1),1,size(av)),cH(size(H,1),size(H,2))
		real(8) :: f(size(E)),tmp,mt(size(av)),fE1,x(size(av)),px(sum(var(:0)%n))
		integer :: l,lv,n,m,i,j,ii,k
		av=0d0
		fE=0d0
		!$omp parallel do reduction(+:fe,av) private(h,e,ch,f,d) schedule(runtime) if (size(brizon_k,1)>1)
		do k=1,size(brizon_k,1)
			call Hamilton(var,H,brizon_k(k,:))
			call mat_diag(H,E)
			cH=transpose(conjg(H))
			f=1d0/(exp(min(E/Tk,1d2))+1d0)
			call dHamilton(var,H,cH,D(:,1:1,:),brizon_k(k,:))
			if(info==1) then
				fE=fE+entropy(E,f)
			else
				fE=fE+fenergy(E,f)
			endif
			do l=1,size(av)
				av(l)=av(l)+sum(real(D(:,1,l))*f(:))
			enddo
		enddo
		!$omp end parallel do
		fE=fE/size(brizon_k,1)
		av=av/size(brizon_k,1)

		fE1=fE/latt%Ns

		!Nambu space
		if(is_sc) then
			do l=lbound(var,1),ubound(var,1)
				if(var(l)%nb==0) then
					do n=1,size(latt%nb(var(l)%Vnb)%bd)
						i=latt%nb(var(l)%Vnb)%bd(n)%i(1)
						j=latt%nb(var(l)%Vnb)%bd(n)%i(2)
						do ii=1,2
							tmp=var(l)%bd(j)*var(l)%V
							if(mod(abs(var(l)%sg),10)==4) then
								tmp=-tmp*0.5d0
							elseif(mod(abs(var(l)%sg),10)==5) then
								tmp=-tmp
							endif
							fE=fE+tmp*var(l)%val(var(l)%bd2v(j))
							lv=sum(var(:l)%n)-var(l)%n+var(l)%bd2v(i)
							av(lv)=av(lv)+tmp
							if(.not.var(l)%is_V) then
								exit
							else
								call swap(i,j)
							endif
						enddo
					enddo
				endif
			enddo
		endif

		x=0d0
		mt=0d0
		do l=lbound(var,1),ubound(var,1)
			if(var(l)%sg==1) then
				lv=sum(var(:l)%n)-var(l)%n+1
				fE=fE-nf*var(l)%val(1)*real(sum(var(1)%bd))
				av(lv)=av(lv)-nf*real(sum(var(l)%bd))
				x(lv)=var(l)%val(1)
				cycle
			endif
			do n=1,size(latt%nb(var(l)%Vnb)%bd)
				i=latt%nb(var(l)%Vnb)%bd(n)%i(1)
				j=latt%nb(var(l)%Vnb)%bd(n)%i(2)
				if(var(l)%nb==var(l)%Vnb) then
					lv=sum(var(:l)%n)-var(l)%n+var(l)%bd2v(n)
					if(var(l)%is_V) then
						tmp=real(conjg(var(l)%bd(n))*var(l)%bd(n))*var(l)%V
						if(var(l)%nb/=0) then
							tmp=tmp*spin
						endif
						fE=fE-tmp*var(l)%val(var(l)%bd2v(n))**2
						mt(lv)=mt(lv)+tmp*2d0
					else
						do m=1,size(var(l)%bd2v)
							if(var(l)%bd2v(m)==var(l)%bd2v(n)) then
								exit
							endif
						enddo
						if(abs(real(var(l)%bd(m)))>1d-10) then
							tmp=real(var(l)%bd(m))*spin
						else
							tmp=imag(var(l)%bd(m))*spin
						endif
						if(var(l)%nb/=0) then
							tmp=tmp*2d0
						endif
						mt(lv)=mt(lv)+tmp
					endif
				else
					tmp=real(var(l)%bd(i)*var(l)%bd(j))*var(l)%V
					fE=fE-tmp*var(l)%val(var(l)%bd2v(i))*var(l)%val(var(l)%bd2v(j))
					lv=sum(var(:l)%n)-var(l)%n+var(l)%bd2v(i)
					mt(lv)=mt(lv)+tmp
					lv=sum(var(:l)%n)-var(l)%n+var(l)%bd2v(j)
					mt(lv)=mt(lv)+tmp
				endif
			enddo
			lv=sum(var(:l)%n)-var(l)%n
			x(lv+1:lv+var(l)%n)=1d0/mt(lv+1:lv+var(l)%n)*av(lv+1:lv+var(l)%n)
			av(lv+1:lv+var(l)%n)=av(lv+1:lv+var(l)%n)-mt(lv+1:lv+var(l)%n)*var(l)%val(:)
		enddo
		av=av/latt%Ns
		v=av(sum(var(:ll)%n)-var(ll)%n+1:sum(var(:ll)%n)-var(ll)%n+size(v))

		px=put(var(:0))
		do l=lbound(var,1),0
			do n=1,size(var(l)%bd)
				if(abs(var(l)%val(var(l)%bd2v(n)))>1d-17) then
					var(l)%bd(n)=var(l)%bd(n)*var(l)%val(var(l)%bd2v(n))
				endif
			enddo
		enddo
		call get(var,x)
		if(info==1) then
			fE=fE1+free_energy()
		else
			fE=fE/latt%Ns
		endif
		call get(var(:0),px)
		do l=lbound(var,1),0
			do n=1,size(var(l)%bd)
				if(abs(var(l)%val(var(l)%bd2v(n)))>1d-17) then
					var(l)%bd(n)=var(l)%bd(n)/var(l)%val(var(l)%bd2v(n))
				endif
			enddo
		enddo
	end subroutine
	function default_free_energy()
		integer :: l,n
		real(8) :: default_free_energy
		default_free_energy=0d0
		do l=lbound(var,1),ubound(var,1)
			if(var(l)%sg==1) then
				cycle
			endif
			do n=1,size(latt%nb(var(l)%Vnb)%bd)
				default_free_energy=default_free_energy+Eavg(l,n)
			enddo
		enddo
		default_free_energy=default_free_energy/latt%Ns
	end function
	function Eavg(l,n)
		integer :: l,n
		real(8) :: Eavg
		integer :: i,j,m
		i=latt%nb(var(l)%Vnb)%bd(n)%i(1)
		j=latt%nb(var(l)%Vnb)%bd(n)%i(2)
		if(var(l)%nb==var(l)%Vnb) then
			if(var(l)%is_V) then
				if(var(l)%nb/=0) then
					Eavg=real(conjg(var(l)%bd(n))*var(l)%bd(n))*var(l)%V*var(l)%val(var(l)%bd2v(n))**2*spin
				else
					Eavg=real(conjg(var(l)%bd(n))*var(l)%bd(n))*var(l)%V*var(l)%val(var(l)%bd2v(n))**2
				endif
			else
				do m=1,size(var(l)%bd2v)
					if(var(l)%bd2v(m)==var(l)%bd2v(n)) then
						exit
					endif
				enddo
				if(abs(real(var(l)%bd(m)))>1d-10) then
					Eavg=real(var(l)%bd(m))*var(l)%val(var(l)%bd2v(n))*spin
				else
					Eavg=imag(var(l)%bd(m))*var(l)%val(var(l)%bd2v(n))*spin
				endif
				if(var(l)%nb/=0) then
					Eavg=Eavg*2d0
				endif
			endif
		else
			Eavg=real(var(l)%bd(i)*var(l)%bd(j))*var(l)%V*var(l)%val(var(l)%bd2v(i))*var(l)%val(var(l)%bd2v(j))
		endif
	end function
	function Green(gm,k,omg)
		real(8) :: k(:),gm,omg
		complex(8) :: H(latt%Ns*spin,latt%Ns*spin),Green
		real(8) :: E(size(H,1))
		integer :: n,i,j
		call Hamilton(var,H,k)
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
		!$omp parallel do private(h,e) reduction(+:d)
		do k=1,size(brizon_k,1)
			call Hamilton(var,H,brizon_k(k,:))
			call mat_diag(H,E)
			do l=1,m
				D(l)=D(l)+imag(sum(H(i,:)*conjg(H(i,:))/(omg(1)+domg*l-E+img*gm)))
			enddo
		enddo
		!$omp end parallel do
		do l=1,m
			write(ut,"(es17.9$)")omg(1)+domg*l,D(l)/(size(brizon_k,1))
			write(ut,"(x)")
		enddo
	end subroutine
	subroutine EDC(ut,gm,k,omg,m)
		integer :: ut,m
		real(8) :: gm,omg(2),k(:)
		complex(8) :: H(latt%Ns*spin,latt%Ns*spin)
		real(8) :: E(size(H,1)),domg,A(m)
		integer :: i,j,l
		call Hamilton(var,H,k)
		call mat_diag(H,E)
		domg=(omg(2)-omg(1))/m
		A=0d0
		!$omp parallel do
		do l=1,m
			do i=1,latt%Ns
				do j=1,latt%Ns
					A(l)=A(l)+imag(sum(H(i,:)*conjg(H(j,:))/(omg(1)+domg*l-E+img*gm)))
				enddo
			enddo
		enddo
		!$omp end parallel do
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
		!$omp parallel do private(h,e) reduction(+:d)
		do k=1,size(brizon_k,1)
			call Hamilton(var,H,brizon_k(k,:))
			call mat_diag(H,E)
			do l=1,m
				do i=1,latt%Ns
					D(l)=D(l)-imag(sum(H(i,:)*conjg(H(i,:))/(omg(1)+domg*l-E+img*gm)))
				enddo
			enddo
		enddo
		!$omp end parallel do
		do l=1,m
			write(ut,"(es17.9$)")omg(1)+domg*l,D(l)/(latt%Ns*size(brizon_k,1))
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
		nk=size(brizon_k,1)
		!$omp parallel do private(h,e,g,q)
		do m=1,nk
			call Hamilton(var,H,brizon_k(m,:))
			call mat_diag(H,E)
			do i=1,latt%Ns
				do j=1,latt%Ns
					G(i,j)=sum(H(i,:)*conjg(H(j,:))/(omg-E+img*gm))
				enddo
			enddo
			do l=0,size(k,1)/nk-1
				q=k(l*nk+m,:)-brizon_k(m,:)
				do i=1,latt%Ns
					do j=1,latt%Ns
						!A(l*nk+m)=A(l*nk+m)+G(i,j)*exp(img*sum(q*(latt%i2r(i,:)-latt%i2r(j,:))))/latt%Ns
						write(*,*)"not implied"
					enddo
				enddo
			enddo
		enddo
		!$omp end parallel do
		do i=1,size(k,1)
			write(ut,"(es17.9$)")k(i,:),A(i)/(latt%Ns*size(brizon_k,1))
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
		nk=size(brizon_k,1)
		!$omp parallel do private(h,e,g,q)
		do m=1,nk
			call Hamilton(var,H,brizon_k(m,:))
			call mat_diag(H,E)
			do i=1,latt%Ns
				do j=1,latt%Ns
					G(i,j,:)=H(i,:)*conjg(H(j,:))
				enddo
			enddo
			do l=0,size(k,1)/nk-1
				q=k(l*nk+m,:)-brizon_k(m,:)
				do i=1,latt%Ns
					do j=1,latt%Ns
						!A(l*nk+m,size(E)+1:)=A(l*nk+m,size(E)+1:)+real(G(i,j,:)*exp(img*sum(q*(latt%i2r(i,:)-latt%i2r(j,:)))))
						write(*,*)"not implied"
					enddo
				enddo
				A(l*nk+m,:size(E))=E
			enddo
		enddo
		!$omp end parallel do
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
		!$omp parallel do private(h,e) reduction(+:a)
		do l1=1,n
			call Hamilton(var,H,ki+dk*l1)
			call mat_diag(H,E)
			do l2=1,m
				do i=1,latt%Ns
					do j=1,latt%Ns
						A(l1,l2)=A(l1,l2)+imag(sum(H(i,:)*conjg(H(j,:))/(omg(1)+domg*l2-E+img*gm)))
					enddo
				enddo
			enddo
		enddo
		!$omp end parallel do
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
			call Hamilton(var,H,k)
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
