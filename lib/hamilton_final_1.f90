! test
module M_Hamilton_final_M
	use M_lattice_final
	use M_const
	use M_utility
	use M_matrix
	use omp_lib
	implicit none
	type c
		character :: l
		integer :: sb=0
		integer :: tp=0
	end type
	type t_var
		integer :: nb
		type(c), allocatable :: c(:,:)
		real(8), allocatable :: sg(:)
		logical, allocatable :: cg(:)
		real(8) :: V
		real(8), allocatable :: val(:)
		complex(8), allocatable :: bd(:)
		integer, allocatable :: bd2v(:)
		integer :: n
		character(:), allocatable :: label
		real(8), allocatable :: extdat(:)
	end type
	type t_ham
		type(t_var), allocatable :: var(:)
		complex(8), allocatable :: Uqi(:,:),Uik(:,:)
		logical, allocatable :: mask(:,:)
		integer, allocatable :: i2tp(:)
		integer, allocatable :: tp2i(:)
		integer, allocatable :: ptp(:)
		integer, allocatable :: H2i(:,:)
		integer, allocatable :: H2s(:,:)
		integer, allocatable :: i2H(:,:)
		integer, allocatable :: s2H(:,:)
		logical :: is_nambu=.false.
		integer :: fer_bos=-1
		integer :: rg(2)=[1,0]
		integer :: Hi=0,Hs=0
	contains
		procedure :: add
		procedure :: init
		procedure :: put
		procedure :: get_val
		procedure :: find_label
		!procedure :: i2H
		procedure :: Hamilton
		procedure :: dHamilton
		procedure :: band
	end type
	real(8) :: Tk,nf
	type(t_brizon), pointer :: pbrizon
	!procedure(real(8)), pointer :: free_energy
	procedure(mat_diag_interface), pointer :: mat_diag
	!procedure(self_consist_interface), pointer :: self_consist
	private init, add, put, get_val, Hamilton, dHamilton, band
	interface
		!subroutine self_consist_interface(fE)
			!real(8), optional :: fE
		!end subroutine
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
	function add(self,nb,ca,n,cg,sg,val,V,label,is_var,extdat) result(idx)
		class(t_ham) :: self
		integer, intent(in) :: nb
		type(c) :: ca(:)
		integer :: n
		logical, optional :: cg(:)
		real(8), optional :: sg(:)
		real(8), intent(in), optional :: val(:)
		real(8), intent(in), optional :: V
		character(*), optional :: label
		logical, optional :: is_var
		real(8), intent(in), optional :: extdat(:)
		integer :: idx
		integer, allocatable :: c(:),ord(:)
		integer :: i,j
		if(present(is_var)) then
			if(is_var) then
				self%rg(2)=self%rg(2)+1
				idx=self%rg(2)
			else
				self%rg(1)=self%rg(1)-1
				idx=self%rg(1)
			endif
		else
			self%rg(2)=self%rg(2)+1
			idx=self%rg(2)
		endif
		if(present(V)) then
			self%var(idx)%V=V
		else
			self%var(idx)%V=nan
		endif
		self%var(idx)%nb=abs(nb)
		if(present(label)) then
			self%var(idx)%label=label
		endif
		allocate(self%var(idx)%c(n,merge(1,size(ca)/n,n==0)))
		self%var(idx)%c=reshape(ca,shape(self%var(idx)%c))
		allocate(self%var(idx)%cg(size(self%var(idx)%c,2)))
		if(present(cg)) then
			self%var(idx)%cg=cg
		else
			self%var(idx)%cg=.false.
		endif
		allocate(self%var(idx)%sg(size(self%var(idx)%c,2)))
		if(present(sg)) then
			self%var(idx)%sg(:)=sg
		else
			self%var(idx)%sg=1d0
		endif
		if(present(extdat)) then
			allocate(self%var(idx)%extdat(size(extdat)))
			self%var(idx)%extdat=extdat
		endif
		if(present(val)) then
			allocate(ord(size(val)))
			ord=[1:size(ord)]
			call qsort(val,ord)
			call collect(val(ord),c)

			allocate(self%var(idx)%bd2v(size(val)))
			allocate(self%var(idx)%bd(size(val)))
			allocate(self%var(idx)%val(size(c)-1))
			self%var(idx)%n=size(c)-1
			do i=1,size(c)-1
				do j=c(i),c(i+1)-1
					self%var(idx)%bd2v(ord(j))=i
				enddo
			enddo
		else
			allocate(self%var(idx)%bd2v(size(latt%nb(self%var(idx)%nb)%bd)))
			allocate(self%var(idx)%bd(size(self%var(idx)%bd2v)))
			allocate(self%var(idx)%val(1))
			self%var(idx)%bd2v=1
			self%var(idx)%n=1
		endif
		self%var(idx)%bd=1d0
		self%var(idx)%val=0d0
	end function
	integer function find_label(self,label) result(i)
		class(t_ham), intent(inout) :: self
		character(*) :: label
		do i=self%rg(1),self%rg(2)
			if(self%var(i)%label==label) then
				exit
			endif
		enddo
		if(i==self%rg(2)+1) then
			i=1/0
		endif
	end function
	subroutine get_val(self,val,idx)
		class(t_ham), intent(inout) :: self
		integer :: idx(:)
		real(8), intent(in) :: val(:)
		integer :: n,l1,l2
		n=0
		do l1=1,size(idx)
			do l2=1,size(self%var(idx(l1))%val)
				n=n+1
				if(n>size(val)) then
					return
				endif
				self%var(idx(l1))%val(l2)=val(n)
			enddo
		enddo
	end subroutine
	function put(self,idx)
		class(t_ham), intent(in) :: self
		integer :: idx(:)
		real(8) :: put(sum(self%var(idx)%n))
		integer :: n,l1,l2
		n=0
		do l1=1,size(idx)
			do l2=1,size(self%var(idx(l1))%val)
				n=n+1
				put(n)=self%var(idx(l1))%val(l2)
			enddo
		enddo
	end function
	!function i2H(self,i,tp,Nc)
		!class(t_ham), intent(in) :: self
		!integer :: i,tp
		!integer :: n
		!integer, optional :: Nc
		!integer :: i2H
		!n=self%tp2i(tp)
		!i2H=self%ptp(n-1)*merge(Nc,1,present(Nc))+sum(latt%i2isb(i,:),self%mask(n,:))
	!end function
	subroutine init(self)
		class(t_ham), intent(inout) :: self
		integer :: lb,i,j,m,n
		type(t_var) :: tmp(self%rg(1):self%rg(2))
		integer :: i2tp(20),tp2i(20)
		logical :: mask(20,8)
		tmp(self%rg(1):self%rg(2))=self%var(self%rg(1):self%rg(2))
		deallocate(self%var)
		allocate(self%var(self%rg(1):self%rg(2)))
		self%var(self%rg(1):self%rg(2))=tmp
		m=0
		i2tp=0
		tp2i=0
		mask=.false.
		do i=self%rg(1),self%rg(2)
			if(size(self%var(i)%c,1)==2) then
				do j=1,size(self%var(i)%c,2)
					do n=1,m
						if(self%var(i)%c(1,j)%tp==i2tp(n)) then
							mask(n,self%var(i)%c(1,j)%sb)=.true.
							exit
						endif
					enddo
					if(n==m+1) then
						m=m+1
						mask(m,self%var(i)%c(1,j)%sb)=.true.
						i2tp(m)=self%var(i)%c(1,j)%tp
					endif
				enddo
			else
				return
			endif
		enddo
		allocate(self%mask(m,latt%sb),self%i2tp(m),self%tp2i(minval(i2tp+merge(100,0,i2tp==0)):maxval(i2tp+merge(-100,0,i2tp==0))))
		self%is_nambu=lbound(self%tp2i,1)<0
		self%mask=mask(:m,:latt%sb)
		self%i2tp=i2tp(:m)
		self%tp2i(self%i2tp)=[1:size(self%i2tp)]
		allocate(self%ptp(0:size(self%mask,1)))
		self%ptp(0)=0
		do i=1,size(self%mask,1)
			self%ptp(i)=self%ptp(i-1)+sum(latt%i2isb(latt%Ni,:),self%mask(i,:))
		enddo
		self%Hi=self%ptp(size(self%mask,1))
		self%Hs=self%Hi*latt%Nc
		allocate(self%H2i(self%Hi,2),self%i2H(latt%Ni,minval(i2tp+merge(100,0,i2tp==0)):maxval(i2tp+merge(-100,0,i2tp==0))))
		if(latt%is_all) allocate(self%H2s(self%Hs,2),self%s2H(latt%Ns,minval(i2tp+merge(100,0,i2tp==0)):maxval(i2tp+merge(-100,0,i2tp==0))))
		do i=1,latt%Ns
			if(.not.latt%is_all.and.i>latt%Ni) exit
			do n=1,size(self%mask,1)
				if(self%mask(n,latt%nb(0)%bd(i)%sb(1))) then
					if(i<=latt%Ni) then
						self%i2H(i,self%i2tp(n))=self%ptp(n-1)+sum(latt%i2isb(i,:),self%mask(n,:))
						self%H2i(self%i2H(i,self%i2tp(n)),1)=i
						self%H2i(self%i2H(i,self%i2tp(n)),2)=self%i2tp(n)
					endif
					if(latt%is_all) then
						self%s2H(i,self%i2tp(n))=self%ptp(n-1)*latt%Nc+sum(latt%i2isb(i,:),self%mask(n,:))
						self%H2s(self%s2H(i,self%i2tp(n)),1)=i
						self%H2s(self%s2H(i,self%i2tp(n)),2)=self%i2tp(n)
					endif
				endif
			enddo
		enddo
		!write(*,*)self%mask
		!write(*,*)self%ptp

		mat_diag => default_mat_diag
		!free_energy => default_free_energy
		pbrizon => brizon
		
		associate(Ns=>latt%Ns,Nc=>latt%Nc,Ni=>latt%Ni,sb=>latt%sb,nq=>pbrizon%nq,nk=>pbrizon%nk,ptp=>self%ptp)
			if(all(abs(latt%bdc(:2))>1d-6)) then
				allocate(self%Uqi(self%Hi,self%Hi))
				self%Uqi=0d0
				!do i=1,nq
					!do n=1,size(self%mask,1)
						!do j=1,self%Hi
							!self%Uqi(ptp(n-1)+count(self%mask(n,:))*(i-1)+count(self%mask(n,:latt%nb(0)%bd(self%H2i(j,1))%sb(1))),j)=exp(-img*sum(pbrizon%q(i,:)*latt%nb(0)%bd(self%H2i(j,1))%r))
						!enddo
					!enddo
				!enddo
				do i=1,self%Hi
					do n=1,nq
						j=latt%isb2i(n,latt%nb(0)%bd(self%H2i(i,1))%sb(1))
						if(latt%nb(0)%bd(j)%sb(1)==latt%nb(0)%bd(self%H2i(i,1))%sb(1)) then
							self%Uqi(i,self%i2H(j,self%H2i(i,2)))=exp(-img*sum(pbrizon%q(latt%i2isb(self%H2i(i,1),latt%nb(0)%bd(self%H2i(i,1))%sb(1)),:)*latt%nb(0)%bd(j)%r))
						else
							stop "sb is not arranged by order"
						endif
					enddo
				enddo
				self%Uqi=self%Uqi*(1d0/sqrt(real(Ni/sb)))
			endif
			if(latt%is_all) then
				allocate(self%Uik(self%Hs,self%Hs))
				self%Uik=0d0
				do i=1,self%Hs
					!do j=1,nk
						!do n=1,size(self%mask,1)
							!self%Uik(i,ptp(n-1)*Nc+count(self%mask(n,:))*(j-1)+count(self%mask(n,:latt%nb(0)%bd(self%H2s(i,1))%sb(1))))=exp(img*sum(latt%nb(0)%bd(self%H2s(i,1))%r*pbrizon%k(j,:)))
						!enddo
					!enddo
					j=latt%nb(0)%bd(self%H2s(i,1))%sbc(1)
					do n=1,nk
						self%Uik(i,self%i2H(j,self%H2s(i,2))+(n-1)*self%Hi)=exp(img*sum(latt%nb(0)%bd(self%H2s(i,1))%r*pbrizon%k(n,:)))
					enddo
				enddo
				self%Uik=self%Uik*(1d0/sqrt(real(Nc)))
			endif
		end associate
	end subroutine
	subroutine Hamilton(self,idx,H,k)
		class(t_ham), intent(in), target :: self
		integer :: idx(:)
		complex(8), intent(inout) :: H(:,:)
		real(8), optional, intent(in) :: k(:)
		complex(8) :: bdc,expk,bd,tmp
		integer :: i,j,l,n,m,nb,hi,hj,sb(2),il,Nc,ir(2)
		real(8) :: dr(3)
		real(8), allocatable :: sg(:)
		integer, pointer :: p2H(:,:)
		character :: flag
		expk=1d0
		if((.not.present(k)).and.latt%is_all) then
			p2H(1:size(self%s2H,1),lbound(self%s2H,2):ubound(self%s2H,2)) => self%s2H
		else
			p2H(1:size(self%i2H,1),lbound(self%i2H,2):ubound(self%i2H,2)) => self%i2H
		endif
		do l=1,size(idx)
			il=idx(l)
			if(size(self%var(il)%c,1)/=2) then
				write(*,*)"warning!!, only quadratic term is considered" 
				cycle
			endif
			nb=abs(self%var(il)%nb)
			do n=1,size(latt%nb(nb)%bd)/merge(latt%Nc,1,(present(k).and.latt%is_all))
				bd=conjg(self%var(il)%bd(n))
				dr=latt%nb(nb)%bd(n)%dr
				! the bdc for spin down in Nambu is conjg(bdc)
				bdc=latt%nb(nb)%bd(n)%bdc
				sb=latt%nb(nb)%bd(n)%sb

				if(abs(bd)<1d-6) cycle
				if((.not.isnan(self%var(il)%V)).or.present(k)) bdc=abs(bdc)
				if(present(k)) then
					i=latt%nb(nb)%bd(n)%sbc(1)
					j=latt%nb(nb)%bd(n)%sbc(2)
					expk=exp(img*sum(k*dr))
				else
					i=latt%nb(nb)%bd(n)%i(1)
					j=latt%nb(nb)%bd(n)%i(2)
				endif
				bd=bd*self%var(il)%val(self%var(il)%bd2v(n))*merge(1d0,self%var(il)%V,isnan(self%var(il)%V))

				flag=" "
				do m=1,size(self%var(il)%c,2)
					if(all(merge(sb(1),sb(2),self%var(il)%c(:,m)%l=="i")-self%var(il)%c(:,m)%sb==0)) then
						if(flag==" ") flag="i"
					elseif(all(merge(sb(2),sb(1),self%var(il)%c(:,m)%l=="i")-self%var(il)%c(:,m)%sb==0)) then
						if(flag==" ") flag="j"
					else
						cycle
					endif
					ir=merge(i,j,self%var(il)%c(1:2,m)%l==flag)
					if(self%var(il)%c(1,m)%l==flag.and.self%var(il)%c(2,m)%l/=flag) then
						tmp=bdc*expk
					elseif(self%var(il)%c(2,m)%l==flag.and.self%var(il)%c(1,m)%l/=flag) then
						tmp=conjg(bdc*expk)
					else
						tmp=1d0
					endif
					ir(1)=p2H(ir(1),self%var(il)%c(1,m)%tp)
					ir(2)=p2H(ir(2),-self%var(il)%c(2,m)%tp)
					H(ir(1),ir(2))=H(ir(1),ir(2))+tmp*merge(conjg(bd),bd,self%var(il)%cg(m))*self%var(il)%sg(m)
				enddo
			enddo
		enddo
	end subroutine
	subroutine dHamilton(self,idx,H,cH,D,k)
		class(t_ham), intent(in), target :: self
		integer, intent(in) :: idx(:)
		complex(8), intent(in) :: H(:,:),cH(:,:)
		complex(8), intent(inout) :: D(:,:,:)
		real(8), optional, intent(in) :: k(:)
		complex(8) :: bdc,expk,bd,tmp
		integer :: i,j,l,n,m,nb,hi,hj,sb(2),il,Nc,m1,m2,lv,ir(2)
		real(8) :: dr(3)
		real(8), allocatable :: sg(:)
		integer, pointer :: p2H(:,:)
		character :: flag
		expk=1d0
		if((.not.present(k)).and.latt%is_all) then
			p2H(1:size(self%s2H,1),lbound(self%s2H,2):ubound(self%s2H,2)) => self%s2H
		else
			p2H(1:size(self%i2H,1),lbound(self%i2H,2):ubound(self%i2H,2)) => self%i2H
		endif
		do l=1,size(idx)
			il=idx(l)
			if(size(self%var(il)%c,1)/=2) then
				write(*,*)"warning!!, only quadratic term is considered" 
				cycle
			endif
			nb=abs(self%var(il)%nb)
			do n=1,size(latt%nb(nb)%bd)/merge(latt%Nc,1,(present(k).and.latt%is_all))
				bd=conjg(self%var(il)%bd(n))
				dr=latt%nb(nb)%bd(n)%dr
				! the bdc for spin down in Nambu is conjg(bdc)
				bdc=latt%nb(nb)%bd(n)%bdc
				sb=latt%nb(nb)%bd(n)%sb

				if(abs(bd)<1d-6) cycle
				if((.not.isnan(self%var(il)%V)).or.present(k)) bdc=abs(bdc)
				if(present(k)) then
					i=latt%nb(nb)%bd(n)%sbc(1)
					j=latt%nb(nb)%bd(n)%sbc(2)
					expk=exp(img*sum(k*dr))
				else
					i=latt%nb(nb)%bd(n)%i(1)
					j=latt%nb(nb)%bd(n)%i(2)
				endif
				bd=bd*merge(1d0,self%var(il)%V,isnan(self%var(il)%V))

				flag=" "
				do m=1,size(self%var(il)%c,2)
					if(all(merge(sb(1),sb(2),self%var(il)%c(:,m)%l=="i")-self%var(il)%c(:,m)%sb==0)) then
						if(flag==" ") flag="i"
					elseif(all(merge(sb(2),sb(1),self%var(il)%c(:,m)%l=="i")-self%var(il)%c(:,m)%sb==0)) then
						if(flag==" ") flag="j"
					else
						cycle
					endif
					ir=merge(i,j,self%var(il)%c(1:2,m)%l==flag)
					if(self%var(il)%c(1,m)%l==flag.and.self%var(il)%c(2,m)%l/=flag) then
						tmp=bdc*expk
					elseif(self%var(il)%c(2,m)%l==flag.and.self%var(il)%c(1,m)%l/=flag) then
						tmp=conjg(bdc*expk)
					else
						tmp=1d0
					endif
					ir(1)=p2H(ir(1),self%var(il)%c(1,m)%tp)
					ir(2)=p2H(ir(2),-self%var(il)%c(2,m)%tp)
					tmp=tmp*merge(conjg(bd),bd,self%var(il)%cg(m))*self%var(il)%sg(m)
					lv=sum(self%var(idx(1:l))%n)-self%var(il)%n+self%var(il)%bd2v(n)
					do m1=1,size(D,1)
						do m2=1,size(D,2)
							D(m1,m2,lv)=D(m1,m2,lv)+cH(m1,ir(1))*H(ir(2),m2)*tmp
						enddo
					enddo
				enddo
			enddo
		enddo
	end subroutine
	!function fenergy(E,f)
		!real(8) :: f(:),E(:)
		!integer :: l
		!real(8) :: fenergy
		!real(8) :: tmp
		!tmp=0d0
		!do l=1,size(E)
			!if((-E(l))>1d2*Tk) then
				!tmp=tmp+E(l)
			!else
				!tmp=tmp-Tk*log(1d0+exp(-E(l)/Tk))
			!endif
		!enddo
		!fenergy=tmp
	!end function
	!function entropy(E,f)
		!real(8) :: f(:),E(:)
		!integer :: l
		!real(8) :: entropy
		!real(8) :: tmp
		!tmp=0d0
		!do l=1,size(E)
			!tmp=tmp-(f(l)*log(max(f(l),1d-17))+(1d0-f(l))*log(max(1d0-f(l),1d-17)))
		!enddo
		!entropy=-tmp*Tk
	!end function
	!subroutine  MF_val(ll,v,fE,info)
		!integer :: ll
		!real(8) :: v(:)
		!real(8) :: fE
		!integer :: info ! 1: F=<H>-TS; else F=-TlnZ

		!complex(8) :: H(latt%Ns*spin,latt%Ns*spin)
		!real(8) :: E(size(H,1)),av(sum(var%n))
		!complex(8) :: D(size(H,1),1,size(av)),cH(size(H,1),size(H,2))
		!real(8) :: f(size(E)),tmp,mt(size(av)),fE1,x(size(av)),px(sum(var(:0)%n))
		!integer :: l,lv,n,m,i,j,ii,k
		!av=0d0
		!fE=0d0
		!!$omp parallel do reduction(+:fe,av) private(h,e,ch,f,d) schedule(runtime) if (size(pbrizon%k,1)>1)
		!do k=1,size(pbrizon%k,1)
			!call self%Hamilton(var,H,pbrizon%k(k,:))
			!call mat_diag(H,E)
			!cH=transpose(conjg(H))
			!f=1d0/(exp(min(E/Tk,1d2))+1d0)
			!call dHamilton(var,H,cH,D(:,1:1,:),pbrizon%k(k,:))
			!if(info==1) then
				!fE=fE+entropy(E,f)
			!else
				!fE=fE+fenergy(E,f)
			!endif
			!do l=1,size(av)
				!av(l)=av(l)+sum(real(D(:,1,l))*f(:))
			!enddo
		!enddo
		!!$omp end parallel do
		!fE=fE/size(pbrizon%k,1)
		!av=av/size(pbrizon%k,1)

		!fE1=fE/latt%Ns

		!!Nambu space
		!if(is_nambu) then
			!do l=lbound(var,1),ubound(var,1)
				!if(var(l)%nb==0) then
					!do n=1,size(latt%nb(var(l)%Vnb)%bd)
						!i=latt%nb(var(l)%Vnb)%bd(n)%i(1)
						!j=latt%nb(var(l)%Vnb)%bd(n)%i(2)
						!do ii=1,2
							!tmp=var(l)%bd(j)*var(l)%V
							!if(mod(abs(var(l)%tp),10)==4) then
								!tmp=-tmp*0.5d0
							!elseif(mod(abs(var(l)%tp),10)==5) then
								!tmp=-tmp
							!endif
							!fE=fE+tmp*var(l)%val(var(l)%bd2v(j))
							!lv=sum(var(:l)%n)-var(l)%n+var(l)%bd2v(i)
							!av(lv)=av(lv)+tmp
							!if(.not.var(l)%is_V) then
								!exit
							!else
								!call swap(i,j)
							!endif
						!enddo
					!enddo
				!endif
			!enddo
		!endif

		!x=0d0
		!mt=0d0
		!do l=lbound(var,1),ubound(var,1)
			!if(var(l)%tp==1) then
				!lv=sum(var(:l)%n)-var(l)%n+1
				!fE=fE-nf*var(l)%val(1)*real(sum(var(1)%bd))
				!av(lv)=av(lv)-nf*real(sum(var(l)%bd))
				!x(lv)=var(l)%val(1)
				!cycle
			!endif
			!do n=1,size(latt%nb(var(l)%Vnb)%bd)
				!i=latt%nb(var(l)%Vnb)%bd(n)%i(1)
				!j=latt%nb(var(l)%Vnb)%bd(n)%i(2)
				!if(var(l)%nb==var(l)%Vnb) then
					!lv=sum(var(:l)%n)-var(l)%n+var(l)%bd2v(n)
					!if(var(l)%is_V) then
						!tmp=real(conjg(var(l)%bd(n))*var(l)%bd(n))*var(l)%V
						!if(var(l)%nb/=0) then
							!tmp=tmp*spin
						!endif
						!fE=fE-tmp*var(l)%val(var(l)%bd2v(n))**2
						!mt(lv)=mt(lv)+tmp*2d0
					!else
						!do m=1,size(var(l)%bd2v)
							!if(var(l)%bd2v(m)==var(l)%bd2v(n)) then
								!exit
							!endif
						!enddo
						!if(abs(real(var(l)%bd(m)))>1d-10) then
							!tmp=real(var(l)%bd(m))*spin
						!else
							!tmp=imag(var(l)%bd(m))*spin
						!endif
						!if(var(l)%nb/=0) then
							!tmp=tmp*2d0
						!endif
						!mt(lv)=mt(lv)+tmp
					!endif
				!else
					!tmp=real(var(l)%bd(i)*var(l)%bd(j))*var(l)%V
					!fE=fE-tmp*var(l)%val(var(l)%bd2v(i))*var(l)%val(var(l)%bd2v(j))
					!lv=sum(var(:l)%n)-var(l)%n+var(l)%bd2v(i)
					!mt(lv)=mt(lv)+tmp
					!lv=sum(var(:l)%n)-var(l)%n+var(l)%bd2v(j)
					!mt(lv)=mt(lv)+tmp
				!endif
			!enddo
			!lv=sum(var(:l)%n)-var(l)%n
			!x(lv+1:lv+var(l)%n)=1d0/mt(lv+1:lv+var(l)%n)*av(lv+1:lv+var(l)%n)
			!av(lv+1:lv+var(l)%n)=av(lv+1:lv+var(l)%n)-mt(lv+1:lv+var(l)%n)*var(l)%val(:)
		!enddo
		!av=av/latt%Ns
		!v=av(sum(var(:ll)%n)-var(ll)%n+1:sum(var(:ll)%n)-var(ll)%n+size(v))

		!px=put(var(:0))
		!do l=lbound(var,1),0
			!do n=1,size(var(l)%bd)
				!if(abs(var(l)%val(var(l)%bd2v(n)))>1d-17) then
					!var(l)%bd(n)=var(l)%bd(n)*var(l)%val(var(l)%bd2v(n))
				!endif
			!enddo
		!enddo
		!call get(var,x)
		!if(info==1) then
			!fE=fE1+free_energy()
		!else
			!fE=fE/latt%Ns
		!endif
		!call get(var(:0),px)
		!do l=lbound(var,1),0
			!do n=1,size(var(l)%bd)
				!if(abs(var(l)%val(var(l)%bd2v(n)))>1d-17) then
					!var(l)%bd(n)=var(l)%bd(n)/var(l)%val(var(l)%bd2v(n))
				!endif
			!enddo
		!enddo
	!end subroutine
	!function default_free_energy()
		!integer :: l,n
		!real(8) :: default_free_energy
		!default_free_energy=0d0
		!do l=lbound(var,1),ubound(var,1)
			!if(var(l)%tp==1) then
				!cycle
			!endif
			!do n=1,size(latt%nb(var(l)%Vnb)%bd)
				!default_free_energy=default_free_energy+Eavg(l,n)
			!enddo
		!enddo
		!default_free_energy=default_free_energy/latt%Ns
	!end function
	!function Eavg(l,n)
		!integer :: l,n
		!real(8) :: Eavg
		!integer :: i,j,m
		!i=latt%nb(var(l)%Vnb)%bd(n)%i(1)
		!j=latt%nb(var(l)%Vnb)%bd(n)%i(2)
		!if(var(l)%nb==var(l)%Vnb) then
			!if(var(l)%is_V) then
				!if(var(l)%nb/=0) then
					!Eavg=real(conjg(var(l)%bd(n))*var(l)%bd(n))*var(l)%V*var(l)%val(var(l)%bd2v(n))**2*spin
				!else
					!Eavg=real(conjg(var(l)%bd(n))*var(l)%bd(n))*var(l)%V*var(l)%val(var(l)%bd2v(n))**2
				!endif
			!else
				!do m=1,size(var(l)%bd2v)
					!if(var(l)%bd2v(m)==var(l)%bd2v(n)) then
						!exit
					!endif
				!enddo
				!if(abs(real(var(l)%bd(m)))>1d-10) then
					!Eavg=real(var(l)%bd(m))*var(l)%val(var(l)%bd2v(n))*spin
				!else
					!Eavg=imag(var(l)%bd(m))*var(l)%val(var(l)%bd2v(n))*spin
				!endif
				!if(var(l)%nb/=0) then
					!Eavg=Eavg*2d0
				!endif
			!endif
		!else
			!Eavg=real(var(l)%bd(i)*var(l)%bd(j))*var(l)%V*var(l)%val(var(l)%bd2v(i))*var(l)%val(var(l)%bd2v(j))
		!endif
	!end function
	!function Green(gm,k,omg)
		!real(8) :: k(:),gm,omg
		!complex(8) :: H(latt%Ns*spin,latt%Ns*spin),Green
		!real(8) :: E(size(H,1))
		!integer :: n,i,j
		!call Hamilton(var,H,k)
		!call mat_diag(H,E)
		!Green=0d0
		!do i=1,latt%Ns
			!do j=1,latt%Ns
				!Green=Green+sum(H(i,:)*conjg(H(j,:))/(omg-E+img*gm))/latt%Ns
			!enddo
		!enddo
	!end function
	!subroutine LDOS(ut,gm,i,omg,m)
		!integer :: ut,i,m
		!real(8) :: gm,omg(2)
		!complex(8) :: H(latt%Ns*spin,latt%Ns*spin)
		!real(8) :: E(size(H,1)),domg,D(m)
		!integer :: k,l
		!domg=(omg(2)-omg(1))/m
		!D=0d0
		!!$omp parallel do private(h,e) reduction(+:d)
		!do k=1,size(pbrizon%k,1)
			!call Hamilton(var,H,pbrizon%k(k,:))
			!call mat_diag(H,E)
			!do l=1,m
				!D(l)=D(l)+imag(sum(H(i,:)*conjg(H(i,:))/(omg(1)+domg*l-E+img*gm)))
			!enddo
		!enddo
		!!$omp end parallel do
		!do l=1,m
			!write(ut,"(es17.9$)")omg(1)+domg*l,D(l)/(size(pbrizon%k,1))
			!write(ut,"(x)")
		!enddo
	!end subroutine
	!subroutine EDC(ut,gm,k,omg,m)
		!integer :: ut,m
		!real(8) :: gm,omg(2),k(:)
		!complex(8) :: H(latt%Ns*spin,latt%Ns*spin)
		!real(8) :: E(size(H,1)),domg,A(m)
		!integer :: i,j,l
		!call Hamilton(var,H,k)
		!call mat_diag(H,E)
		!domg=(omg(2)-omg(1))/m
		!A=0d0
		!!$omp parallel do
		!do l=1,m
			!do i=1,latt%Ns
				!do j=1,latt%Ns
					!A(l)=A(l)+imag(sum(H(i,:)*conjg(H(j,:))/(omg(1)+domg*l-E+img*gm)))
				!enddo
			!enddo
		!enddo
		!!$omp end parallel do
		!do l=1,m
			!write(ut,"(2es17.9)")omg(1)+domg*l,A(l)
			!write(ut,"(x)")
		!enddo
	!end subroutine
	!subroutine DOS(ut,gm,omg,m,peak)
		!integer :: ut,m
		!real(8) :: gm,omg(:)
		!real(8), allocatable, optional :: peak(:)
		!complex(8) :: H(latt%Ns*spin,latt%Ns*spin)
		!real(8) :: E(size(H,1)),D(m),x,domg
		!integer :: k,i,l
		!integer, allocatable :: ipeak(:)
		!domg=(omg(2)-omg(1))/m
		!D=0d0
		!!$omp parallel do private(h,e) reduction(+:d)
		!do k=1,size(pbrizon%k,1)
			!call Hamilton(var,H,pbrizon%k(k,:))
			!call mat_diag(H,E)
			!do l=1,m
				!do i=1,latt%Ns
					!D(l)=D(l)-imag(sum(H(i,:)*conjg(H(i,:))/(omg(1)+domg*l-E+img*gm)))
				!enddo
			!enddo
		!enddo
		!!$omp end parallel do
		!do l=1,m
			!write(ut,"(es17.9$)")omg(1)+domg*l,D(l)/(latt%Ns*size(pbrizon%k,1))
			!write(ut,"(x)")
		!enddo
		!write(ut,"(x/)")
		!if(present(peak)) then
			!if(allocated(peak)) then
				!deallocate(peak)
			!endif
			!call find_peak(D,ipeak)
			!allocate(peak(size(ipeak)))
			!do i=1,size(ipeak)
				!peak(i)=omg(1)+domg*ipeak(i)
			!enddo
		!endif
	!end subroutine
	!subroutine fermis(ut,gm,k,omg)
		!integer :: ut
		!real(8) :: gm,omg,k(:,:)
		!complex(8) :: A(size(k,1)),G(latt%Ns,latt%Ns),H(latt%Ns*spin,latt%Ns*spin)
		!real(8) :: q(3),E(size(H,1))
		!integer :: i,j,l,m,nk
		!A=0d0
		!nk=size(pbrizon%k,1)
		!!$omp parallel do private(h,e,g,q)
		!do m=1,nk
			!call Hamilton(var,H,pbrizon%k(m,:))
			!call mat_diag(H,E)
			!do i=1,latt%Ns
				!do j=1,latt%Ns
					!G(i,j)=sum(H(i,:)*conjg(H(j,:))/(omg-E+img*gm))
				!enddo
			!enddo
			!do l=0,size(k,1)/nk-1
				!q=k(l*nk+m,:)-pbrizon%k(m,:)
				!do i=1,latt%Ns
					!do j=1,latt%Ns
						!!A(l*nk+m)=A(l*nk+m)+G(i,j)*exp(img*sum(q*(latt%i2r(i,:)-latt%i2r(j,:))))/latt%Ns
						!write(*,*)"not implied"
					!enddo
				!enddo
			!enddo
		!enddo
		!!$omp end parallel do
		!do i=1,size(k,1)
			!write(ut,"(es17.9$)")k(i,:),A(i)/(latt%Ns*size(pbrizon%k,1))
			!write(ut,"(x)")
		!enddo
		!write(ut,"(x/)")
	!end subroutine
	!subroutine energy(ut,k)
		!integer :: ut
		!real(8) :: k(:,:)
		!complex(8) :: H(latt%Ns*spin,latt%Ns*spin),G(latt%Ns,latt%Ns,size(H,1)*2)
		!real(8) :: E(size(H,1)),A(size(k,1),size(H,1)*2),q(3)
		!integer :: nk,m,l,i,j
		!A=0d0
		!nk=size(pbrizon%k,1)
		!!$omp parallel do private(h,e,g,q)
		!do m=1,nk
			!call Hamilton(var,H,pbrizon%k(m,:))
			!call mat_diag(H,E)
			!do i=1,latt%Ns
				!do j=1,latt%Ns
					!G(i,j,:)=H(i,:)*conjg(H(j,:))
				!enddo
			!enddo
			!do l=0,size(k,1)/nk-1
				!q=k(l*nk+m,:)-pbrizon%k(m,:)
				!do i=1,latt%Ns
					!do j=1,latt%Ns
						!!A(l*nk+m,size(E)+1:)=A(l*nk+m,size(E)+1:)+real(G(i,j,:)*exp(img*sum(q*(latt%i2r(i,:)-latt%i2r(j,:)))))
						!write(*,*)"not implied"
					!enddo
				!enddo
				!A(l*nk+m,:size(E))=E
			!enddo
		!enddo
		!!$omp end parallel do
		!do i=1,size(k,1)
			!write(ut,"(es17.9$)")k(i,:)
			!write(ut,"(es17.9$)")(A(i,j),j=1,size(E)),(A(i,j+size(E)),j=1,size(E))
			!write(ut,"(x)")
		!enddo
	!end subroutine
	!subroutine band_e(ut,gm,ki,kf,n,omg,m)
		!integer :: ut,n,m
		!real(8) :: ki(3),kf(3),gm,omg(2)
		!complex(8) :: H(latt%Ns*spin,latt%Ns*spin)
		!real(8) :: A(n,m),dk(3),domg,E(size(H,1))
		!integer :: l1,l2,i,j
		!domg=(omg(2)-omg(1))/m
		!dk=(kf-ki)/n
		!A=0d0
		!!$omp parallel do private(h,e) reduction(+:a)
		!do l1=1,n
			!call Hamilton(var,H,ki+dk*l1)
			!call mat_diag(H,E)
			!do l2=1,m
				!do i=1,latt%Ns
					!do j=1,latt%Ns
						!A(l1,l2)=A(l1,l2)+imag(sum(H(i,:)*conjg(H(j,:))/(omg(1)+domg*l2-E+img*gm)))
					!enddo
				!enddo
			!enddo
		!enddo
		!!$omp end parallel do
		!do l1=1,n
			!do l2=1,m
				!write(ut,"(es17.9$)")ki+dk*l1,omg(1)+domg*l2,A(l1,l2)
				!write(ut,"(x)")
			!enddo
			!write(ut,"(x)")
		!enddo
		!write(ut,"(x/)")
	!end subroutine
	subroutine band(self,ut,ki,kf,n)
		class(t_ham) :: self
		integer :: ut
		integer :: n
		real(8) :: ki(3),kf(3)
		real(8) :: k(3),dk(3)
		integer :: i,l,j,m
		complex(8) :: H(self%Hi,self%Hi)
		real(8) :: E(size(H,1)),G
		dk=(kf-ki)/n
		do m=0,n-1
			k=ki+dk*m
			H=0d0
			call self%Hamilton([self%rg(1):self%rg(2)],H,k)
			call mat_diag(H,E)
			write(ut,"((es17.9)$)")k
			do l=1,size(E)
				G=0d0
				do i=1,latt%Ni
					do j=1,latt%Ni
						if(latt%nb(0)%bd(i)%sb(1)==latt%nb(0)%bd(j)%sb(1)) then
							G=G+real(H(i,l)*conjg(H(j,l)))
						endif
					enddo
				enddo
				write(ut,"((es17.9)$)")E(l),G*1d0/pbrizon%nq
			enddo
			write(ut,"(x)")
		enddo
	end subroutine
	!function find_order(l,is,x,rang,z)
		!integer :: l,is
		!real(8) :: find_order,x,rang(3)
		!real(8) :: z
		!real(8) :: pod,order,dx,mdx
		!integer :: isg
		!dx=rang(3)
		!mdx=dx/20d0
		!call self_consist()
		!pod=sum(abs(var(l)%val(:)))/(size(var(l)%val(:)))-z
		!do 
			!x=x-dx*sign(1d0,pod)*is
			!if(x<rang(1).or.x>rang(2)) then
				!if(dx>mdx) then
					!x=x+dx*sign(1d0,pod)*is
					!dx=dx*0.3333d0
					!cycle
				!else
					!find_order=min(max(x,rang(1)),rang(2))
					!exit
				!endif
			!endif
			!call self_consist()
			!order=sum(abs(var(l)%val(:)))/(size(var(l)%val))
			!call is_cross(pod,order-z,isg)
			!if(isg/=0) then
				!if(abs(dx)<mdx) then
					!find_order=x
					!exit
				!endif
				!dx=dx*0.3333d0
			!endif
		!enddo
	!end function
end module
