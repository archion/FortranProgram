module M_pmt
	use M_const
	implicit none
	real(8), parameter :: t(3)=(/1d0,-0.25d0,0d0/),nf=0.85d0,U=2.44d0,DJ=0d0,V=-1d0,cvg=1d-6,Tk=1d-5,Vimp=1d0
	integer :: imp
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
		integer :: iv,i,j,l,ip,k
		type(t_var), allocatable :: tmp(:)
		type(t_mysort), allocatable :: sort_site(:),sort_bd1(:)
		integer, allocatable :: collect_site(:),collect_bd1(:)
		real(8) :: q(3)
		allocate(tmp(10000))
		call init_random_seed()
		q=(/1d0/8d0,0d0,0d0/)
		!q=(/1d0/50d0,0d0,0d0/)
		!q=0d0
		! lattice 
		a1=(/1d0,0d0,0d0/)
		a2=(/0d0,1d0,0d0/)
		!T1=(/8d0,8d0,0d0/)
		!T2=(/-6d0,6d0,0d0/)
		T1=a1*24
		T2=a2*24
		bdc(1)=1d0
		bdc(2)=1d0
		allocate(sub(1,2))
		sub(1,:)=(/0d0,0d0/)
		layer=1
		call gen_latt()
		call gen_neb()
		call gen_bond(3)
		! finish i2r(i,2),neb(i)%nb(j)%bond(k)/bdc(k)/r(k,2)
		write(*,*)"Total site number is: ",Ns
		imp=nint(sqrt(real(Ns)))/2*nint(sqrt(real(Ns)))+nint(sqrt(real(Ns)))/2
		!do k=1,Ns
			!write(101,"(3es13.2,2I5)")i2r(k,:),k,ab(k)
		!enddo
		!write(101,"(1x/)")
		allocate(sort_site(Ns))
		do k=1,size(sort_site)
			!sort_site(k)%val=cos(2d0*pi*sum((bond(0)%bd(k)%r)*q))
			!sort_site(k)%val=k
			sort_site(k)%val=sqrt(sum((bond(0)%bd(k)%r-i2r(imp,:))**2))+abs(theta(abs(bond(0)%bd(k)%r-i2r(imp,:)))-pi/4d0)
			sort_site(k)%idx=k
			write(101,"(es13.2$)")bond(0)%bd(k)%r,bond(0)%bd(k)%dir,sort_site(k)%val
			write(101,"(x)")
		enddo
		write(101,"(1x/)")
		call sort_site%qsort()
		call sort_site%collect(collect_site)

		allocate(sort_bd1(size(bond(1)%bd)))
		do k=1,size(sort_bd1)
			!sort_bd1(k)%val=cos(2d0*pi*sum((bond(1)%bd(k)%r-(/0.5d0,0d0,0d0/))*q))
			!sort_bd1(k)%val=k
			sort_bd1(k)%val=sqrt(sum((bond(1)%bd(k)%r-i2r(imp,:))**2))+abs(theta(abs(bond(1)%bd(k)%r-i2r(imp,:)))-pi/4d0)
			sort_bd1(k)%idx=k
			!write(101,"(es13.2$)")bond(1)%bd(k)%r,bond(1)%bd(k)%dir,sort_bd1(k)%val
			!write(101,"(x)")
		enddo
		!write(101,"(1x/)")
		call sort_bd1%qsort()
		call sort_bd1%collect(collect_bd1)

		iv=0
		! cp
		call gen_var(iv,ip,sg=1,nb=0,V=1d0,var=tmp)
		tmp(iv)%val=1.52d-01
		tmp(iv)%bd_sg=-1d0

		! impure
		call gen_var(iv,ip,(/1,2/),(/imp/),sg=-4,nb=0,V=1d0,var=tmp)
		tmp(iv)%val=Vimp
		tmp(iv)%bd_sg=1d0

		!! ddw
		!!call gen_var(iv,ip,collect_bd1,sort_bd1%idx,sg=3,nb=1,V=1d0,var=tmp)
		!call gen_var(iv,ip,sg=3,nb=1,V=1d0,var=tmp)
		!do i=ip,iv
			!tmp(i)%val=2d-1
			!do k=1,size(tmp(i)%n)
				!tmp(i)%bd_sg(k)=ab(bond(tmp(i)%nb)%bd(tmp(i)%n(k))%i(1))*img*dwave(tmp(i)%n(k))*&
					!1d0
					!!cos(2d0*pi*sum(q*(bond(tmp(i)%nb)%bd(tmp(i)%n(k))%r-(/0.5d0,0d0,0d0/))))
			!enddo
		!enddo

		! d-wave sc
		call gen_var(iv,ip,collect_bd1,sort_bd1%idx,sg=2,nb=1,V=V,var=tmp)
		!call gen_var(iv,ip,sg=2,nb=1,V=V,var=tmp)
		do i=ip,iv
			call random_number(tmp(i)%val)
			!tmp(i)%val=2d-1
			tmp(i)%val=(tmp(i)%val-0.5d0)*0.1d0+5d-1
			do k=1,size(tmp(i)%n)
				tmp(i)%bd_sg(k)=dwave(tmp(i)%n(k))
			enddo
		enddo

		!call gen_var(iv,ip,sg=2,nb=1,V=1d0,var=tmp)
		!do i=ip,iv
			!tmp(i)%val=1d-5
			!do k=1,size(tmp(i)%n)
				!tmp(i)%bd_sg(k)=dwave(tmp(i)%n(k))*&
					!cos(2d0*pi*sum(q*(bond(tmp(i)%nb)%bd(tmp(i)%n(k))%r-(/0.5d0,0d0,0d0/))))
			!enddo
		!enddo


		!! s-wave sc
		!call gen_var(iv,ip,collect_bd1,sort_bd1%idx,sg=2,nb=1,V=1d0,var=tmp)
		!do i=ip,iv
			!tmp(i)%val=1d-1
			!do k=1,size(tmp(i)%n)
				!tmp(i)%bd_sg(k)=1d0
			!enddo
		!enddo

		! sdw
		call gen_var(iv,ip,collect_site,sort_site%idx,sg=4,nb=0,V=-U,var=tmp)
		!call gen_var(iv,ip,sg=4,nb=0,V=-U,var=tmp)
		do i=ip,iv
			call random_number(tmp(i)%val)
			!tmp(i)%val=1d-1
			tmp(i)%val=(tmp(i)%val-0.5d0)*0.1d0+1d-1
			do k=1,size(tmp(i)%n)
				tmp(i)%bd_sg(k)=ab(tmp(i)%n(k))*&
					1d0
					!sin(2d0*pi*sum(q*(bond(tmp(i)%nb)%bd(tmp(i)%n(k))%r+(/0.5d0,0d0,0d0/))))
					!sin(2d0*pi*sum(q*(bond(tmp(i)%nb)%bd(tmp(i)%n(k))%r)))
			enddo
		enddo

		!! d-wave cdw
		!!call gen_var(iv,ip,collect_bd1,sort_bd1%idx,sg=3,nb=1,V=1d0,var=tmp)
		!call gen_var(iv,ip,sg=3,nb=1,V=1d0,var=tmp)
		!do i=ip,iv
			!tmp(i)%val=1d-1
			!do k=1,size(tmp(i)%n)
				!tmp(i)%bd_sg(k)=dwave(tmp(i)%n(k))
			!enddo
		!enddo

		! on site cdw
		call gen_var(iv,ip,collect_site,sort_site%idx,sg=3,nb=0,V=U,var=tmp)
		!call gen_var(iv,ip,sg=3,nb=0,V=U,var=tmp)
		!deallocate(tmp(iv)%bd_sg,tmp(iv)%n)
		!iv=iv-1
		do i=ip,iv
			call random_number(tmp(i)%val)
			!tmp(i)%val=nf/2d0
			tmp(i)%val=(tmp(i)%val-0.5d0)*0.1d0+nf/2d0
			!tmp(i)%val=(tmp(i)%val-0.5d0)*0.1d0
			do k=1,size(tmp(i)%n)
				tmp(i)%bd_sg(k)=&
					1d0
					!cos(2d0*pi*sum(2d0*q*(bond(tmp(i)%nb)%bd(tmp(i)%n(k))%r+(/0.5d0,0d0,0d0/))))
					!cos(2d0*pi*sum(2d0*q*(bond(tmp(i)%nb)%bd(tmp(i)%n(k))%r)))
			enddo
		enddo

		! hp
		!call gen_var(iv,ip,collect_bd1,sort_bd1%idx,sg=3,nb=1,V=1d0,var=tmp)
		!deallocate(tmp(iv)%bd_sg,tmp(iv)%n)
		!iv=iv-1
		call gen_var(iv,ip,sg=-3,nb=1,V=1d0,var=tmp)
		do i=ip,iv
			tmp(i)%val=t(1)
			do k=1,size(tmp(i)%n)
				tmp(i)%bd_sg(k)=&
					-1d0
					!cos(2d0*pi*sum(q*(bond(tmp(i)%nb)%bd(tmp(i)%n(k))%r-(/0.5d0,0d0,0d0/))))
			enddo
		enddo

		k=tmp(iv)%nb+1
		do l=k,size(t)
			call gen_var(iv,ip,sg=-3,nb=l,V=1d0,var=tmp)
			do i=ip,iv
				tmp(i)%bd_sg=-1d0
				tmp(i)%val=t(l)
			enddo
		enddo


		ip=-count(tmp%sg<0)+1
		k=1
		allocate(var(ip:ip+iv-1))
		do l=1,iv
			if(tmp(l)%sg<0) then
				call move_alloc(tmp(l)%bd_sg,var(ip)%bd_sg)
				call move_alloc(tmp(l)%n,var(ip)%n)
				var(ip)%val=tmp(l)%val
				var(ip)%nb=tmp(l)%nb
				var(ip)%V=tmp(l)%V
				var(ip)%sg=tmp(l)%sg      
				ip=ip+1
			else
				call move_alloc(tmp(l)%bd_sg,var(k)%bd_sg)
				call move_alloc(tmp(l)%n,var(k)%n)
				var(k)%val=tmp(l)%val
				var(k)%nb=tmp(l)%nb
				var(k)%V=tmp(l)%V
				var(k)%sg=tmp(l)%sg      
				k=k+1
			endif
		enddo
		deallocate(tmp,sort_site,sort_bd1,collect_site,collect_bd1)
	end subroutine
	function theta(r)
		real(8) :: r(:),theta,d
		d=sqrt(sum(r**2))
		if(d<1d-10) then
			theta=0d0
			return
		endif
		theta=acos(r(1)/d)
		if(r(2)<0d0) then
			theta=2d0*pi-theta
		endif
	end function
	subroutine bdg_var()
		integer :: info
		real(8) :: x(size(var(1:))),v(size(var(1:))),wa(nint((size(var(1:))*(3*size(var(1:))+13))/2.)+10)
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
			if(all(abs(v)<1d-6)) then
				exit
			endif
		enddo
	end subroutine
	subroutine do_var(n,x,v,info)
		integer, intent(in) :: n
		integer, intent(inout) :: info
		real(8), intent(inout) :: x(n),v(n)
		complex(8) :: H(Ns*2,Ns*2),D(size(H,1),size(H,2)),cH(size(H,1),size(H,2)),bd
		real(8) :: E(Ns*2),f(size(E)),dvar,fE,err
		integer :: l,i=0
		i=i+1
		var(1:n)%val=x
		call Hamilton(H)
		call heevd(H,E,"V")
		cH=transpose(conjg(H))
		f=1d0/(exp(E/Tk)+1d0)
		call fenergy(E,fE)
		write(*,"(i5$)"),i
		write(*,"(es15.7$)")fE
		if(n<8) then
			write(*,"(es12.4$)")var(1:n)%val
		else
			write(*,"(es12.4$)")var(1)%val,var(5)%val,var(Ns*2+5)%val,var(Ns*2+Ns+5)%val
		endif
		!!$OMP PARALLEL DO PRIVATE(D)
		do l=1,n
			call dHamilton(var(l),H,cH,D(:,1:1))
			v(l)=sum(real(D(:,1))*f(:))
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
				fE=fE+real(var(l)%val*sum(var(l)%bd_sg)*(1d0-nf))
			case default
				if(var(l)%nb==0) then
					if(var(l)%sg==3) then
						fE=fE-sum(abs(var(l)%val*var(l)%bd_sg)**2)*var(l)%V+real(var(l)%val*sum(var(l)%bd_sg)*var(l)%V)
					else
						fE=fE-sum(abs(var(l)%val*var(l)%bd_sg)**2)*var(l)%V-real(var(l)%val*sum(var(l)%bd_sg)*var(l)%V)
					endif
				else
					fE=fE-2d0*sum(abs(var(l)%val*var(l)%bd_sg)**2)*var(l)%V
				endif
			end select
		enddo
		fE=fE/Ns
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
				write(ut,"(es17.9$)")bond(var(l)%nb)%bd(var(l)%n(k))%r,bond(var(l)%nb)%bd(var(l)%n(k))%dir,var(l)%val
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
	i=1
	call initial()
	!var(1)%val=1.519265227d-01
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
	!call export_data(10)
	!!stop
	call bdg_var()
	!call bdg()
	!export data
	call export_data(10)
end program
