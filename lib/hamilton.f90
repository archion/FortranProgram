module M_Hamilton
	use M_lattice
	use M_const
	use M_utility
	use lapack95, only: heevd
	implicit none
	type, extends(t_sort), private :: t_mysort
		integer :: idx
	contains
		procedure :: swap_sort => myswap1
	end type
	type t_var
		real(8) :: val
		integer :: nb
		complex(8), allocatable :: bd_sg(:)
		integer, allocatable :: n(:)
		real(8) :: V
		integer :: sg   
						! 1:  chemical potainial
						! 2:  pair
						! 3:  charge
						! 4:  spin
						! 11: n-n jast
						! -: don't do variation
	end type
	type(t_var), allocatable :: var(:)
	integer :: spin=2
	logical :: is_sc=.true.
	integer :: iv(-1:1)=(/1,0,0/)
	integer, private :: getn(-10:10)=0,ni(-1:1)=(/1,1,1/)
	private myswap1 
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
	subroutine gen_var(n,val,sg,nb,V,var)
		integer, intent(in) :: sg,nb
		real(8), intent(in), optional :: val(:)
		real(8), intent(in) :: V
		type(t_var), intent(inout) :: var(:)
		integer, intent(inout) :: n
		integer :: i,bs
		type(t_mysort), allocatable :: sort(:)
		integer, allocatable :: collect(:)
		n=0
		bs=(size(var)-1)/2+1
		if(present(val)) then

			allocate(sort(size(val)))
			sort(:)%val=val
			sort(:)%idx=(/(i,i=1,size(sort))/)
			call sort%qsort()
			call sort%collect(collect)

			do i=1,size(collect)-1
				iv(sign(1,sg))=iv(sign(1,sg))+sign(1,sg)
				n=n+sign(1,sg)
				var(iv(sign(1,sg))+bs)%sg=sg
				var(iv(sign(1,sg))+bs)%nb=nb
				var(iv(sign(1,sg))+bs)%V=V
				allocate(var(iv(sign(1,sg))+bs)%bd_sg(collect(i+1)-collect(i)),var(iv(sign(1,sg))+bs)%n(collect(i+1)-collect(i)))
				var(iv(sign(1,sg))+bs)%n=sort(collect(i):collect(i+1)-1)%idx
			enddo
		else
			iv(sign(1,sg))=iv(sign(1,sg))+sign(1,sg)
			n=n+sign(1,sg)
			var(iv(sign(1,sg))+bs)%sg=sg
			var(iv(sign(1,sg))+bs)%nb=nb
			var(iv(sign(1,sg))+bs)%V=V
			allocate(var(iv(sign(1,sg))+bs)%bd_sg(size(latt%bond(nb)%bd)),var(iv(sign(1,sg))+bs)%n(size(latt%bond(nb)%bd)))
			var(iv(sign(1,sg))+bs)%n=(/(i,i=1,size(latt%bond(nb)%bd))/)
		endif
		iv(0)=iv(sign(1,sg))
		ni(sign(1,sg))=ni(sign(1,sg))+sign(1,sg)
		getn(ni(sign(1,sg)))=iv(0)
	end subroutine
	function get_var(sg,nb,l)
		integer, intent(in) :: sg,nb
		integer, intent(in) :: l
		integer :: get_var(2)
		integer :: n,i,j
		n=0
		do i=1,size(getn(1:))
			if(var(getn(i)+1)%sg==sg.and.var(getn(i)+1)%nb==nb) then
				n=n+1
				if(n==l) then
					get_var(1)=getn(i)+1
					get_var(2)=getn(i+1)
					exit
				endif
				if(i>1.and.getn(i)==0) then
					write(*,*)"get_var fail"
					exit
				endif
			endif
		enddo
	end function
	subroutine Hamilton(H,rk)
		complex(8), intent(inout) :: H(:,:)
		real(8), optional, intent(in) :: rk(:)
		complex(8) :: bd,expk
		integer :: i,j,k,l,isc
		H=0d0
		if(is_sc) then
			isc=-1
		else
			isc=1
		endif
		expk=1d0
		associate(Ns=>latt%Ns)
			do l=lbound(var,1),ubound(var,1)
				do k=1,size(var(l)%bd_sg)
					i=latt%bond(var(l)%nb)%bd(var(l)%n(k))%i(1)
					j=latt%bond(var(l)%nb)%bd(var(l)%n(k))%i(2)
					bd=latt%bond(var(l)%nb)%bd(var(l)%n(k))%bdc
					if(present(rk)) then
						expk=exp(img*sum(rk*latt%bond(var(l)%nb)%bd(var(l)%n(k))%dr))
					endif
					select case(abs(var(l)%sg))
					case(2)
						! pair channel
						if(.not.is_sc) then
							write(*,*)"is_sc is not setted as true"
							stop
						endif
						H(i,j+Ns)=H(i,j+Ns)+var(l)%val*var(l)%bd_sg(k)*bd*var(l)%V*expk
						H(i+Ns,j)=H(i+Ns,j)+var(l)%val*var(l)%bd_sg(k)*conjg(bd)*var(l)%V*expk
						if(var(l)%nb/=0) then
							H(j+Ns,i)=H(j+Ns,i)+conjg(var(l)%val*var(l)%bd_sg(k)*bd*var(l)%V*expk)
							H(j,i+Ns)=H(j,i+Ns)+conjg(var(l)%val*var(l)%bd_sg(k)*conjg(bd)*var(l)%V*expk)
						endif
					case(3,1)
						! charge channel
						H(i,j)=H(i,j)+var(l)%val*conjg(var(l)%bd_sg(k))*bd*var(l)%V*expk
						if(var(l)%nb/=0) then
							H(j,i)=H(j,i)+conjg(var(l)%val*conjg(var(l)%bd_sg(k))*bd*var(l)%V*expk)
						endif
						if(spin==2) then
							H(i+Ns,j+Ns)=H(i+Ns,j+Ns)+isc*var(l)%val*var(l)%bd_sg(k)*bd*var(l)%V*expk
							if(var(l)%nb/=0) then
								H(j+Ns,i+Ns)=H(j+Ns,i+Ns)+isc*conjg(var(l)%val*var(l)%bd_sg(k)*bd*var(l)%V*expk)
							endif
						endif
					case(4)
						! spin channel
						H(i,j)=H(i,j)+var(l)%val*conjg(var(l)%bd_sg(k))*bd*var(l)%V*expk
						H(i+Ns,j+Ns)=H(i+Ns,j+Ns)-isc*var(l)%val*var(l)%bd_sg(k)*bd*var(l)%V*expk
						if(var(l)%nb/=0) then
							H(j,i)=H(j,i)+conjg(var(l)%val*conjg(var(l)%bd_sg(k))*bd*var(l)%V*expk)
							H(j+Ns,i+Ns)=H(j+Ns,i+Ns)-isc*conjg(var(l)%val*var(l)%bd_sg(k)*bd*var(l)%V*expk)
						endif
					end select
				enddo
			enddo
		end associate
	end subroutine
	pure subroutine dHamilton(var,H,cH,D,rk)
		type(t_var), intent(in) :: var
		complex(8), intent(in) :: H(:,:),cH(:,:)
		real(8), optional, intent(in) :: rk(:)
		complex(8), intent(inout) :: D(:,:)
		complex(8) :: bd,expk
		integer :: i,j,k,l,m,n,u,isc
		D=0d0
		if(is_sc) then
			isc=-1
		else
			isc=1
		endif
		expk=1d0
		associate(Ns=>latt%Ns)
			do k=1,size(var%bd_sg)
				i=latt%bond(var%nb)%bd(var%n(k))%i(1)
				j=latt%bond(var%nb)%bd(var%n(k))%i(2)
				bd=latt%bond(var%nb)%bd(var%n(k))%bdc
				if(present(rk)) then
					expk=exp(img*sum(rk*latt%bond(var%nb)%bd(var%n(k))%dr))
				endif
				do u=1,size(D,2)
					do m=1,size(D,1)
						if(size(D,2)==1) then
							n=m 
						else
							n=u
						endif
						select case(var%sg)
						case(2)
							! pair channel
							D(m,u)=D(m,u)+cH(m,i)*H(j+Ns,n)*var%bd_sg(k)*bd*var%V*expk
							D(m,u)=D(m,u)+cH(m,i+Ns)*H(j,n)*var%bd_sg(k)*conjg(bd)*var%V*expk
							if(var%nb/=0) then
								D(m,u)=D(m,u)+cH(m,j+Ns)*H(i,n)*conjg(var%bd_sg(k)*bd*var%V*expk)
								D(m,u)=D(m,u)+cH(m,j)*H(i+Ns,n)*conjg(var%bd_sg(k)*conjg(bd)*var%V*expk)
							endif
						case(3,1)
							! charge channel
							D(m,u)=D(m,u)+cH(m,i)*H(j,n)*conjg(var%bd_sg(k))*bd*var%V*expk
							if(var%nb/=0) then
								D(m,u)=D(m,u)+cH(m,j)*H(i,n)*conjg(conjg(var%bd_sg(k))*bd*var%V*expk)
							endif
							if(spin==2) then
								D(m,u)=D(m,u)+isc*cH(m,i+Ns)*H(j+Ns,n)*var%bd_sg(k)*bd*var%V*expk
								if(var%nb/=0) then
									D(m,u)=D(m,u)+isc*cH(m,j+Ns)*H(i+Ns,n)*conjg(var%bd_sg(k)*bd*var%V*expk)
								endif
							endif
						case(4)
							! spin channel
							D(m,u)=D(m,u)+cH(m,i)*H(j,n)*conjg(var%bd_sg(k))*bd*var%V*expk
							D(m,u)=D(m,u)-isc*cH(m,i+Ns)*H(j+Ns,n)*var%bd_sg(k)*bd*var%V*expk
							if(var%nb/=0) then
								D(m,u)=D(m,u)+cH(m,j)*H(i,n)*conjg(conjg(var%bd_sg(k))*bd*var%V*expk)
								D(m,u)=D(m,u)-isc*cH(m,j+Ns)*H(i+Ns,n)*conjg(var%bd_sg(k)*bd*var%V*expk)
							endif
						end select
					enddo
				enddo
			enddo
		end associate
	end subroutine
	!subroutine ddEnergy(E,H,cH,D,ddE)
		!complex(8), intent(in) :: H(:,:),cH(:,:),D(:,:,:)
		!real(8), intent(in) :: E(:)
		!real(8), intent(inout) :: ddE(:,:,:)
		!real(8) :: tmp,iE(size(E),size(E))
		!integer :: l,l1,l2,n,m
		!ddE(:,:,:)=0d0
		!do n=1,size(E)
			!do m=1,size(E)
				!if(abs(E(n)-E(m))<1d-7) then
					!iE(m,n)=0d0
					!cycle
				!endif
				!iE(m,n)=1d0/(E(n)-E(m))
			!enddo
		!enddo
		!!$OMP PARALLEL DO PRIVATE(tmp)
		!do l1=1,size(ddE,1)
			!!do l2=1,size(ddE,1)
			!do l2=l1,l1
				!do n=1,size(E)
					!tmp=0d0
					!do m=1,size(E)
						!tmp=tmp+real(conjg(D(m,n,l1))*D(m,n,l2))*iE(m,n)
					!enddo
					!ddE(n,l1,l2)=tmp
				!enddo
			!enddo
			!ddE(:,l1,l1)=ddE(:,l1,l1)*2d0
		!enddo
		!!$OMP END PARALLEL DO
	!end subroutine
	subroutine energy(ut)
		integer :: ut,n,i
		real(8) :: ki(3),kf(3),k(3),dk(3)
		complex(8) :: H(latt%Ns*spin,latt%Ns*spin)
		real(8) :: E(size(H,1))
		do i=1,size(brizon%k,1)
			call Hamilton(H,brizon%k(i,:))
			call heevd(H,E,"V")
			write(ut,"(es17.9$)")brizon%k(i,:),E,abs(H(1,:))
			write(ut,"(x)")
		enddo
	end subroutine
	subroutine band(ki,kf,n,ut)
		integer :: ut,n,i
		real(8) :: ki(3),kf(3),k(3),dk(3)
		complex(8) :: H(latt%Ns*spin,latt%Ns*spin)
		real(8) :: E(size(H,1))
		dk=(kf-ki)/n
		do i=0,n-1
			k=ki+dk*i
			call Hamilton(H,k)
			call heevd(H,E,"V")
			write(ut,"(es17.9$)")k,E,abs(H(1,:))
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
