module M_pmt
	use M_const
	implicit none
	real(8), parameter :: t(2)=(/1d0,-0.25d0/),nf=0.85d0,U=2.54d0,DJ=0d0,V=-1d0,cvg=1d-6,Tk=1d-5,Vimp=1d0
	complex(8), allocatable :: H(:,:)
	real(8), allocatable :: E(:)
end module
module M_bdg
	use M_Hamilton
	use M_rd
	use M_utility
	use M_pmt
	implicit none
contains
	subroutine initial()
		integer :: i,j,l,n,k
		type(t_var), allocatable :: tmp(:)
		real(8), allocatable :: bd0(:),bd1(:)
		real(8) :: q(3)
		allocate(tmp(-30000:30000))
		call init_random_seed()

		! lattice 
		latt%a1=(/1d0,1d0,0d0/)
		latt%a2=(/-1d0,1d0,0d0/)
		latt%T1=(/1d0,0d0,0d0/)*24
		latt%T2=(/0d0,1d0,0d0/)*24
		latt%bdc(1)=1d0
		latt%bdc(2)=1d0
		allocate(latt%sub(2,3))
		latt%sub(1,:)=(/0d0,0d0,0d0/)
		latt%sub(2,:)=(/1d0,0d0,0d0/)
		latt%layer=1
		call latt%gen_latt()
		call latt%gen_neb(size(t))
		call latt%gen_bond(size(t))
		brizon%n1=8
		brizon%n2=8
		call latt%gen_brizon()
		write(*,*)"Total site number is: ",latt%Ns
		call check_lattice(101)

		allocate(bd0(latt%Ns))
		do k=1,size(bd0)
			bd0(k)=k
			!sort_site(k)%val=sqrt(sum((latt%bond(0)%bd(k)%r-latt%i2r(imp,:))**2))+abs(theta(abs(latt%bond(0)%bd(k)%r-latt%i2r(imp,:)))-pi/4d0)
		enddo

		allocate(bd1(size(latt%bond(1)%bd)))
		do k=1,size(bd1)
			bd1(k)=k
			!sort_bd1(k)%val=sqrt(sum((latt%bond(1)%bd(k)%r-latt%i2r(imp,:))**2))+abs(theta(abs(latt%bond(1)%bd(k)%r-latt%i2r(imp,:)))-pi/4d0)
		enddo

		! cp
		call gen_var(n,sg=1,nb=0,V=1d0,var=tmp)
		tmp(iv(0))%val=0d0
		tmp(iv(0))%bd_sg=-1d0

		! impure
		call gen_var(n,(/1d0/),sg=-4,nb=0,V=1d0,var=tmp)
		tmp(iv(0))%n=145
		tmp(iv(0))%val=Vimp
		tmp(iv(0))%bd_sg=1d0

		! d-wave sc
		call gen_var(n,bd1,sg=2,nb=1,V=V,var=tmp)
		!call gen_var(n,sg=2,nb=1,V=V,var=tmp)
		do i=iv(0)-n+sign(1,n),iv(0),sign(1,n)
			call random_number(tmp(i)%val)
			tmp(i)%val=(tmp(i)%val-0.5d0)*0.1d0+5d-1
			tmp(i)%val=2d-1
			do k=1,size(tmp(i)%n)
				tmp(i)%bd_sg(k)=dwave(tmp(i)%n(k))
			enddo
		enddo

		! sdw
		call gen_var(n,bd0,sg=4,nb=0,V=-U,var=tmp)
		!call gen_var(n,sg=4,nb=0,V=-U,var=tmp)
		do i=iv(0)-n+sign(1,n),iv(0),sign(1,n)
			call random_number(tmp(i)%val)
			tmp(i)%val=(tmp(i)%val-0.5d0)*0.1d0+1d-1
			tmp(i)%val=1d-2
			do k=1,size(tmp(i)%n)
				tmp(i)%bd_sg(k)=ab(tmp(i)%n(k))*&
					1d0
			enddo
		enddo


		! on site cdw
		call gen_var(n,bd0,sg=3,nb=0,V=U,var=tmp)
		!call gen_var(n,sg=3,nb=0,V=U,var=tmp)
		do i=iv(0)-n+sign(1,n),iv(0),sign(1,n)
			call random_number(tmp(i)%val)
			tmp(i)%val=(tmp(i)%val-0.5d0)*0.1d0+nf/2d0
			tmp(i)%val=nf/2d0
			do k=1,size(tmp(i)%n)
				tmp(i)%bd_sg(k)=&
					1d0
			enddo
		enddo

		! hp
		do l=1,size(t)
			call gen_var(n,sg=-3,nb=l,V=1d0,var=tmp)
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
		deallocate(tmp,bd0,bd1)
	end subroutine
	subroutine self_consist()
		integer :: info
		real(8) :: x(size(var(1:))),v(size(var(1:))),wa(nint((size(var(1:))*(3*size(var(1:))+13))/2.)+10)
		allocate(H(latt%Ns*spin,latt%Ns*spin),E(latt%Ns*spin))
		x=var(1:)%val
		info=1
		call export_data(10)
		call hybrd1(do_var,size(x(:)),x(:),v(:),1d-7,info,wa,size(wa))
		var(1:)%val=x
		call export_data(10)
		stop
		do 
			if(var(1)%sg==1) then
				call hybrd1(do_var,size(x(1:1)),x(1:1),v(1:1),1d-7,info,wa(:16),16)
				x=var(1:)%val
			endif
			write(*,"(A$)")"****"
			info=-1
			call do_var(size(x),x,v,info=info)
			var(1:)%val=x
			if(info<0) then
				exit
			endif
		enddo
	end subroutine
	subroutine do_var(n,x,v,info)
		integer, intent(in) :: n
		integer, intent(inout) :: info
		real(8), intent(inout) :: x(n),v(n)
		complex(8) :: D(size(H,1),1,n),cH(size(H,1),size(H,2)),bd
		real(8) :: f(size(E)),dvar(size(var(1:))),fE,err,dE(size(H,1),n)
		integer :: l,i=0,k
		v=0d0
		dvar=0d0
		fE=0d0
		var(1:n)%val=x
		!$OMP PARALLEL DO REDUCTION(+:fE,v,dvar) PRIVATE(H,E,cH,f,dE,D) if(size(brizon%k,1)/=1)
		do k=1,size(brizon%k,1)
			if(info>0.or.size(brizon%k,1)/=1) then
				call Hamilton(H,brizon%k(k,:))
				call heevd(H,E,"V")
				!write(20,"(es17.9)")E
			endif
			cH=transpose(conjg(H))
			f=1d0/(exp(E/Tk)+1d0)
			do l=1,n
				call dHamilton(var(l),H,cH,D(:,1:1,l),brizon%k(k,:))
			enddo
			dE=real(D(:,1,:))
			do l=1,size(E)
				if((-E(l)/Tk)>1d2) then
					fE=fE+E(l)
				else
					fE=fE-Tk*log(1d0+exp(-E(l)/Tk))
				endif
			enddo
			do l=1,size(var(1:))
				if(l<=n) v(l)=v(l)+sum(dE(:,l)*f(:))
				select case(var(l)%sg)
				case(1)
					if(l<=n) v(l)=v(l)+real(sum(var(l)%bd_sg)*(1d0-nf))
					fE=fE+var(l)%val*real(sum(var(l)%bd_sg)*(1d0-nf))
				case default
					if(var(l)%nb==0) then
						dvar(l)=dvar(l)-2d0*sum(abs(var(l)%bd_sg)**2)*var(l)%V
						if(var(l)%sg==3) then
							if(l<=n) v(l)=v(l)+real(sum(var(l)%bd_sg)*var(l)%V)
							fE=fE+var(l)%val*real(sum(var(l)%bd_sg)*var(l)%V)
						else
							if(l<=n) v(l)=v(l)-real(sum(var(l)%bd_sg)*var(l)%V)
							fE=fE-var(l)%val*real(sum(var(l)%bd_sg)*var(l)%V)
						endif
					else
						dvar(l)=dvar(l)-4d0*sum(abs(var(l)%bd_sg)**2)*var(l)%V
					endif
				end select
			enddo
		enddo
		!$OMP END PARALLEL DO
		do l=1,n
			if(info<0.and.abs(dvar(l))>1d-8) then
				x(l)=-v(l)/dvar(l)
				!x(l)=-v(l)/dvar(l)*0.6d0+var(l)%val*0.4d0
			endif
			v(l)=(v(l)+dvar(l)*var(l)%val)/(size(var(l)%n)*size(brizon%k,1))
		enddo
		fE=(fE+sum(dvar*var(1:)%val**2)/2d0)/(latt%Ns*size(brizon%k,1))
		err=sum(abs(v))/size(v)
		info=1
		if(err<1d-7) then
			info=-1
		endif

		if(n>1) then
			i=i+1
			write(20,"(es17.9$)")fE
			write(20,"(es17.9)")err
			!call export_data(10)
			write(*,"(i5$)"),i
		endif
		write(*,"(es15.7$)")fE
		if(n<8) then
			write(*,"(es12.4$)")var(1:n)%val,v(1:n)
		else
			write(*,"(es12.4$)")var(1)%val,var(2)%val,var(latt%Ns*2+2)%val,var(latt%Ns*2+latt%Ns+2)%val,err
		endif
		write(*,"(x)")

	end subroutine
	subroutine export_data(ut)
		integer :: ut,l,k,i(2),s(3,3)=(/2,1,1,  4,0,1,  3,0,1/),j
		do j=1,size(s,2)
			i=get_var(s(1,j),s(2,j),s(3,j))
			do l=i(1),i(2)
				do k=1,size(var(l)%n)
					write(ut,"(es17.9$)")latt%bond(var(l)%nb)%bd(var(l)%n(k))%r,latt%bond(var(l)%nb)%bd(var(l)%n(k))%dr,var(l)%val,var(l)%bd_sg(k)
					write(ut,"(x)")
				enddo
			enddo
			write(ut,"(x/)")
		enddo
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
