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
	contains
		procedure :: Hamilton
	end type
	type(t_var), allocatable :: var(:)
	integer :: spin=2
	logical :: is_sc=.true.
	integer :: iv(-1:1)=(/1,0,0/)
	integer, private :: getn(-20:20,2)=0,ni(-1:1)=(/1,0,0/)
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
		getn(ni(sign(1,sg)),:)=(/iv(0)-n+sign(1,n),iv(0)/)
	end subroutine
	function get_var(sg,nb,l)
		integer, intent(in) :: sg,nb
		integer, intent(in) :: l
		integer :: get_var(2)
		integer :: n,i,j
		n=0
		do i=ni(-1),ni(1)
			if(var(getn(i,1))%sg==sg.and.var(getn(i,1))%nb==nb) then
				n=n+1
				if(n==l) then
					get_var=getn(i,:)
					exit
				endif
			endif
		enddo
		if(n/=l) then
			write(*,*)"get_var fail"
			stop
		endif
	end function
	subroutine Hamilton(self,H,rk,fn,info)
		class(t_var), intent(in) :: self(:)
		complex(8), intent(inout) :: H(:,:)
		real(8), optional, intent(in) :: rk(:)
		complex(8), external, optional :: fn
		integer, optional, intent(in) :: info
		complex(8) :: bdc,expk,bd
		integer :: i,j,k,l,isc
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
			do l=lbound(self,1),ubound(self,1)
				do k=1,size(self(l)%bd_sg)
					i=latt%bond(self(l)%nb)%bd(self(l)%n(k))%i(1)
					j=latt%bond(self(l)%nb)%bd(self(l)%n(k))%i(2)
					bdc=latt%bond(self(l)%nb)%bd(self(l)%n(k))%bdc
					bd=self(l)%bd_sg(k)
					if(present(rk)) then
						expk=exp(-img*sum(rk*latt%bond(self(l)%nb)%bd(self(l)%n(k))%dr))
					endif
					if(present(fn)) then
						bd=bd*fn(latt%bond(self(l)%nb)%bd(self(l)%n(k))%dr)
					else
					endif
					select case(abs(self(l)%sg))
					case(2)
						! pair channel
						if(.not.is_sc) then
							write(*,*)"is_sc is not setted as true"
							stop
						endif
						if(present(fn)) then
							cycle
						endif
						H(i,j+Ns)=H(i,j+Ns)+self(l)%val*bd*bdc*self(l)%V*expk
						H(i+Ns,j)=H(i+Ns,j)+self(l)%val*bd*conjg(bdc)*self(l)%V*expk
						if(self(l)%nb/=0) then
							H(j+Ns,i)=H(j+Ns,i)+conjg(self(l)%val*bd*bdc*self(l)%V*expk)
							H(j,i+Ns)=H(j,i+Ns)+conjg(self(l)%val*bd*conjg(bdc)*self(l)%V*expk)
						endif
					case(3,1)
						! charge channel
						H(i,j)=H(i,j)+self(l)%val*conjg(bd)*bdc*self(l)%V*expk
						if(self(l)%nb/=0) then
							H(j,i)=H(j,i)+conjg(self(l)%val*conjg(bd)*bdc*self(l)%V*expk)
						endif
						if(spin==2) then
							H(i+Ns,j+Ns)=H(i+Ns,j+Ns)+isc*self(l)%val*bd*bdc*self(l)%V*expk
							if(self(l)%nb/=0) then
								H(j+Ns,i+Ns)=H(j+Ns,i+Ns)+isc*conjg(self(l)%val*bd*bdc*self(l)%V*expk)
							endif
						endif
					case(4)
						! spin channel
						H(i,j)=H(i,j)+self(l)%val*conjg(bd)*bdc*self(l)%V*expk
						H(i+Ns,j+Ns)=H(i+Ns,j+Ns)-isc*self(l)%val*bd*bdc*self(l)%V*expk
						if(self(l)%nb/=0) then
							H(j,i)=H(j,i)+conjg(self(l)%val*conjg(bd)*bdc*self(l)%V*expk)
							H(j+Ns,i+Ns)=H(j+Ns,i+Ns)-isc*conjg(self(l)%val*bd*bdc*self(l)%V*expk)
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
		complex(8) :: bdc,expk,bd
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
				bdc=latt%bond(var%nb)%bd(var%n(k))%bdc
				bd=var%bd_sg(k)
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
							D(m,u)=D(m,u)+cH(m,i)*H(j+Ns,n)*bd*bdc*var%V*expk
							D(m,u)=D(m,u)+cH(m,i+Ns)*H(j,n)*bd*conjg(bdc)*var%V*expk
							if(var%nb/=0) then
								D(m,u)=D(m,u)+cH(m,j+Ns)*H(i,n)*conjg(bd*bdc*var%V*expk)
								D(m,u)=D(m,u)+cH(m,j)*H(i+Ns,n)*conjg(bd*conjg(bdc)*var%V*expk)
							endif
						case(3,1)
							! charge channel
							D(m,u)=D(m,u)+cH(m,i)*H(j,n)*conjg(bd)*bdc*var%V*expk
							if(var%nb/=0) then
								D(m,u)=D(m,u)+cH(m,j)*H(i,n)*conjg(conjg(bd)*bdc*var%V*expk)
							endif
							if(spin==2) then
								D(m,u)=D(m,u)+isc*cH(m,i+Ns)*H(j+Ns,n)*bd*bdc*var%V*expk
								if(var%nb/=0) then
									D(m,u)=D(m,u)+isc*cH(m,j+Ns)*H(i+Ns,n)*conjg(bd*bdc*var%V*expk)
								endif
							endif
						case(4)
							! spin channel
							D(m,u)=D(m,u)+cH(m,i)*H(j,n)*conjg(bd)*bdc*var%V*expk
							D(m,u)=D(m,u)-isc*cH(m,i+Ns)*H(j+Ns,n)*bd*bdc*var%V*expk
							if(var%nb/=0) then
								D(m,u)=D(m,u)+cH(m,j)*H(i,n)*conjg(conjg(bd)*bdc*var%V*expk)
								D(m,u)=D(m,u)-isc*cH(m,j+Ns)*H(i+Ns,n)*conjg(bd*bdc*var%V*expk)
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
