module M_pmt
	use M_const
	implicit none
	real(8), parameter :: t(3)=(/1d0,-0.25d0,0d0/),nf=0.85d0,U=2.44d0,DJ=0d0,V=-1d0,cvg=1d-6,Tk=1d-5,Vimp=1d0
	integer :: imp
	complex(8), allocatable :: H(:,:)
	real(8), allocatable :: E(:)
end module
module M_bdg
	use M_Hamilton
	use M_rd
	use M_utility
	use M_pmt
	use lapack95, only: heevd
	implicit none
	type, extends(t_sort), private :: t_mysort
		integer :: idx
	contains
		procedure :: swap_sort => myswap1
	end type
contains
	subroutine initial()
		integer :: iv(-1:1),i,j,l,n,k
		type(t_var), allocatable :: tmp(:)
		type(t_mysort), allocatable :: sort_site(:),sort_bd1(:)
		integer, allocatable :: collect_site(:),collect_bd1(:)
		real(8) :: q(3)
		allocate(tmp(-10000:10000))
		call init_random_seed()
		! lattice 
		latt%a1=(/1d0,0d0,0d0/)
		latt%a2=(/0d0,1d0,0d0/)
		latt%T1=latt%a1*24
		latt%T2=latt%a2*24
		latt%bdc(1)=1d0
		latt%bdc(2)=1d0
		allocate(latt%sub(1,2))
		latt%sub(1,:)=(/0d0,0d0/)
		latt%layer=1
		call latt%gen_latt()
		call latt%gen_neb(3)
		call latt%gen_bond(3)
		! finish i2r(i,2),neb(i)%nb(j)%bond(k)/bdc(k)/r(k,2)
		write(*,*)"Total site number is: ",latt%Ns
		imp=nint(sqrt(real(latt%Ns)))/2*nint(sqrt(real(latt%Ns)))+nint(sqrt(real(latt%Ns)))/2
		!do k=1,latt%Ns
			!write(101,"(3es13.2,2I5)")i2r(k,:),k,ab(k)
		!enddo
		!write(101,"(1x/)")
		allocate(sort_site(latt%Ns))
		do k=1,size(sort_site)
			!sort_site(k)%val=cos(2d0*pi*sum((bond(0)%bd(k)%r)*q))
			sort_site(k)%val=k
			!sort_site(k)%val=sqrt(sum((latt%bond(0)%bd(k)%r-latt%i2r(imp,:))**2))+abs(theta(abs(latt%bond(0)%bd(k)%r-latt%i2r(imp,:)))-pi/4d0)
			sort_site(k)%idx=k
			write(101,"(es13.2$)")latt%bond(0)%bd(k)%r,latt%bond(0)%bd(k)%dir,sort_site(k)%val
			write(101,"(i5$)")k
			write(101,"(x)")
		enddo
		write(101,"(1x/)")
		call sort_site%qsort()
		call sort_site%collect(collect_site)

		allocate(sort_bd1(size(latt%bond(1)%bd)))
		do k=1,size(sort_bd1)
			!sort_bd1(k)%val=cos(2d0*pi*sum((latt%bond(1)%bd(k)%r-(/0.5d0,0d0,0d0/))*q))
			sort_bd1(k)%val=k
			!sort_bd1(k)%val=sqrt(sum((latt%bond(1)%bd(k)%r-latt%i2r(imp,:))**2))+abs(theta(abs(latt%bond(1)%bd(k)%r-latt%i2r(imp,:)))-pi/4d0)
			sort_bd1(k)%idx=k
			write(101,"(es13.2$)")latt%bond(1)%bd(k)%r,latt%bond(1)%bd(k)%dir,sort_bd1(k)%val
			write(101,"(i5$)")k
			write(101,"(x)")
		enddo
		write(101,"(1x/)")
		call sort_bd1%qsort()
		call sort_bd1%collect(collect_bd1)

		iv=(/1,0,0/)
		! cp
		call gen_var(iv,n,sg=1,nb=0,V=1d0,var=tmp)
		tmp(iv(0))%val=1.52d-01
		tmp(iv(0))%bd_sg=-1d0

		! impure
		call gen_var(iv,n,(/1,2/),(/imp/),sg=-4,nb=0,V=1d0,var=tmp)
		tmp(iv(0))%val=Vimp
		tmp(iv(0))%bd_sg=1d0

		! d-wave sc
		call gen_var(iv,n,collect_bd1,sort_bd1%idx,sg=2,nb=1,V=V,var=tmp)
		!call gen_var(iv,n,sg=2,nb=1,V=V,var=tmp)
		do i=iv(0)-n+sign(1,n),iv(0),sign(1,n)
			call random_number(tmp(i)%val)
			!tmp(i)%val=2d-1
			tmp(i)%val=(tmp(i)%val-0.5d0)*0.1d0+5d-1
			do k=1,size(tmp(i)%n)
				tmp(i)%bd_sg(k)=dwave(tmp(i)%n(k))
			enddo
		enddo

		! sdw
		call gen_var(iv,n,collect_site,sort_site%idx,sg=4,nb=0,V=-U,var=tmp)
		!call gen_var(iv,n,sg=4,nb=0,V=-U,var=tmp)
		do i=iv(0)-n+sign(1,n),iv(0),sign(1,n)
			call random_number(tmp(i)%val)
			!tmp(i)%val=1d-1
			tmp(i)%val=(tmp(i)%val-0.5d0)*0.1d0+1d-1
			do k=1,size(tmp(i)%n)
				tmp(i)%bd_sg(k)=ab(tmp(i)%n(k))*&
					1d0
			enddo
		enddo


		! on site cdw
		call gen_var(iv,n,collect_site,sort_site%idx,sg=3,nb=0,V=U,var=tmp)
		!call gen_var(iv,n,sg=3,nb=0,V=U,var=tmp)
		!deallocate(tmp(iv)%bd_sg,tmp(iv)%n)
		!iv=iv-1
		do i=iv(0)-n+sign(1,n),iv(0),sign(1,n)
			call random_number(tmp(i)%val)
			!tmp(i)%val=nf/2d0
			tmp(i)%val=(tmp(i)%val-0.5d0)*0.1d0+nf/2d0
			!tmp(i)%val=(tmp(i)%val-0.5d0)*0.1d0
			do k=1,size(tmp(i)%n)
				tmp(i)%bd_sg(k)=&
					1d0
			enddo
		enddo

		! hp
		do l=1,size(t)
			call gen_var(iv,n,sg=-3,nb=l,V=1d0,var=tmp)
			do i=iv(0)-n+sign(1,n),iv(0),sign(1,n)
				tmp(i)%bd_sg=-1d0
				tmp(i)%val=t(l)
			enddo
		enddo


		allocate(var(iv(-1):iv(1)))
		do l=iv(-1),iv(1)
			call move_alloc(tmp(l)%bd_sg,var(l)%bd_sg)
			call move_alloc(tmp(l)%n,var(l)%n)
			var(l)%val=tmp(l)%val
			var(l)%nb=tmp(l)%nb
			var(l)%V=tmp(l)%V
			var(l)%sg=tmp(l)%sg      
		enddo
		deallocate(tmp,sort_site,sort_bd1,collect_site,collect_bd1)
	end subroutine
	subroutine self_consist()
		integer :: info
		real(8) :: x(size(var(1:))),v(size(var(1:))),wa(nint((size(var(1:))*(3*size(var(1:))+13))/2.)+10)
		allocate(H(latt%Ns*2,latt%Ns*2),E(latt%Ns*2))
		!call hybrd1(do_var,size(var(1:)),var(1:)%val,v,1d-7,info,wa,size(wa))
		do 
			if(var(1)%sg==1) then
				x(1)=var(1)%val
				!call hybrd1(do_var,size(var(1:1)),var(1:1)%val,v(1:1),1d-7,info,wa(:16),16)
				call hybrd1(do_var,size(x(1:1)),x(1:1),v(1:1),1d-7,info,wa(:16),16)
				!var(1)%val=x(1)
			endif
			info=-1
			write(*,"(A$)")"****"
			call do_var(size(var(1:)),var(1:)%val,v(1:),info=info)
			if(info<0) then
				exit
			endif
		enddo
	end subroutine
	subroutine do_var(n,x,v,info)
		integer, intent(in) :: n
		integer, intent(inout) :: info
		real(8), intent(inout) :: x(n),v(n)
		complex(8) :: D(size(H,1),size(H,2),n),cH(size(H,1),size(H,2)),bd
		real(8) :: f(size(E)),dvar,fE,err,dE(size(H,1),n)
		integer :: l,i=0
		i=i+1
		var(1:n)%val=x
		if(info>0) then
			call Hamilton(H)
			call heevd(H,E,"V")
		endif
		cH=transpose(conjg(H))
		f=1d0/(exp(E/Tk)+1d0)
		call fenergy(E,fE)
		write(*,"(i5$)"),i
		write(*,"(es15.7$)")fE
		if(n<8) then
			write(*,"(es12.4$)")var(1:n)%val
		else
			write(*,"(es12.4$)")var(1)%val,var(5)%val,var(latt%Ns*2+5)%val,var(latt%Ns*2+latt%Ns+5)%val
		endif
		do l=1,n
			call dHamilton(var(l),H,cH,D(:,1:1,l))
		enddo
		dE=real(D(:,1,:))
		!!$OMP PARALLEL DO PRIVATE(dvar)
		do l=1,n
			v(l)=sum(dE(:,l)*f(:))
			dvar=0d0
			select case(var(l)%sg)
			case(1)
				v(l)=v(l)+real(sum(var(l)%bd_sg)*(1d0-nf))
			case default
				dvar=-2d0*sum(abs(var(l)%bd_sg)**2)*var(l)%V
				if(var(l)%nb==0) then
					if(var(l)%sg==3) then
						v(l)=v(l)+real(sum(var(l)%bd_sg)*var(l)%V)
					else
						v(l)=v(l)-real(sum(var(l)%bd_sg)*var(l)%V)
					endif
				else
					dvar=dvar*2d0
				endif
			end select
			if(info<0.and.abs(dvar)>1d-8) then
				x(l)=-v(l)/dvar
			endif
			v(l)=(v(l)+dvar*var(l)%val)/(size(var(l)%n))
		enddo
		!!$OMP END PARALLEL DO
		err=sum(abs(v))/size(v)
		if(err<1d-7) then
			info=-1
		else
			info=1
		endif
		if(n<8) then
			write(*,"(es12.4$)")v(1:n)
		else
			write(*,"(es12.4$)")err
		endif
		write(*,"(x)")
		if(n>1) then
			write(20,"(es17.9$)")fE
			write(20,"(es17.9)")err
			call export_data(10)
		endif
	end subroutine
	subroutine fenergy(E,fE)
		real(8), intent(in) :: E(:)
		real(8), intent(inout) :: fE
		real(8) :: v2
		integer :: k,l
		fE=-Vimp
		do k=1,size(E)
			if((-E(k)/Tk)>1d2) then
				fE=fE+E(k)
			else
				fE=fE-Tk*log(1d0+exp(-E(k)/Tk))
			endif
		enddo
		do l=1,ubound(var,1)
			select case(var(l)%sg)
			case(1)
				v2=0d0
				fE=fE+real(var(l)%val*sum(var(l)%bd_sg)*(1d0-nf))
			case default
				v2=-sum(abs(var(l)%val*var(l)%bd_sg)**2)*var(l)%V
				if(var(l)%nb==0) then
					if(var(l)%sg==3) then
						fE=fE+real(var(l)%val*sum(var(l)%bd_sg)*var(l)%V)
					else
						fE=fE-real(var(l)%val*sum(var(l)%bd_sg)*var(l)%V)
					endif
				else
					v2=v2*2d0
				endif
				fE=fE+v2
			end select
		enddo
		fE=fE/latt%Ns
	end subroutine
	subroutine export_data(ut)
		integer :: ut,l,k
		!rewind(ut)
		!open(unit=ut,file='../data/order.dat')
		do l=1,ubound(var,1)
			if(var(l)%sg/=var(max(l-1,1))%sg.or.var(l)%nb/=var(max(l-1,1))%nb) then
				write(ut,"(x/)")
			endif
			do k=1,size(var(l)%n)
				write(ut,"(es17.9$)")latt%bond(var(l)%nb)%bd(var(l)%n(k))%r,latt%bond(var(l)%nb)%bd(var(l)%n(k))%dir,var(l)%val,var(l)%bd_sg(k)
				write(ut,"(x)")
			enddo
		enddo
		write(ut,"(x/)")
		!close(ut)
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
program main
	use M_bdg
	implicit none
	logical :: f
	real(8) :: tmp(1152)
	integer :: i,l,err
	f=openfile(unit=10,file='../data/order.dat')
	f=openfile(unit=20,file='../data/fenergy.dat')
	f=openfile(unit=101,file='../data/lattice.dat')
	call initial()
	!var(1)%val=1.519266194D-01
	!open(30,file="../data/order_save.dat",status="old",action="read")
	!read(30,"(1152es17.9)",iostat=err)tmp(1:1152)
	!l=2
	!do
		!var(l)%val=0d0
		!do i=1,size(var(l)%n)
			!var(l)%val=var(l)%val+tmp(var(l)%n(i))
		!enddo
		!var(l)%val=var(l)%val/size(var(l)%n)
		!l=l+1
		!if(var(l)%sg/=var(max(l-1,1))%sg.or.var(l)%nb/=var(max(l-1,1))%nb) then
			!exit
		!endif
	!enddo
	!read(30,"(576es17.9)",iostat=err)tmp(1:576)
	!do 
		!var(l)%val=0d0
		!do i=1,size(var(l)%n)
			!var(l)%val=var(l)%val+tmp(var(l)%n(i))
		!enddo
		!var(l)%val=var(l)%val/size(var(l)%n)
		!l=l+1
		!if(var(l)%sg/=var(max(l-1,1))%sg.or.var(l)%nb/=var(max(l-1,1))%nb) then
			!exit
		!endif
	!enddo
	!read(30,"(576es17.9)",iostat=err)tmp(1:576)
	!do 
		!var(l)%val=0d0
		!do i=1,size(var(l)%n)
			!var(l)%val=var(l)%val+tmp(var(l)%n(i))
		!enddo
		!var(l)%val=var(l)%val/size(var(l)%n)
		!l=l+1
		!if(var(l)%sg/=var(max(l-1,1))%sg.or.var(l)%nb/=var(max(l-1,1))%nb.or.l>ubound(var,1)) then
			!exit
		!endif
	!enddo
	call self_consist()
	call export_data(10)
end program
