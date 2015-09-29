module M_pmt
	use M_const
	implicit none
	real(8), parameter :: t(2)=(/1d0,-0.25d0/),U=2.54d0,DJ=0d0,V=-1d0,cvg=1d-6,Vimp=1d0
end module
module M_bdg
	use M_Hamilton
	use M_pmt
	implicit none
contains
	subroutine initial()
		integer :: i,l,l2,imp
		real(8), allocatable :: bd0(:),bd1(:)
		real(8) :: q(3),r(3)
		allocate(var(-10000:10000))
		call init_random_seed()
		q=(/0d0,1d0/8d0,0d0/)*2d0*pi

		! lattice 
		latt%a1=(/1d0,1d0,0d0/)
		latt%a2=(/-1d0,1d0,0d0/)
		latt%T1=(/1d0,0d0,0d0/)*2
		latt%T2=(/0d0,1d0,0d0/)*8
		!latt%T1=latt%a1
		!latt%T2=latt%a2
		latt%bdc(1)=1d0
		latt%bdc(2)=1d0
		allocate(latt%sub(2,3))
		latt%sub(1,:)=(/0d0,0d0,0d0/)
		latt%sub(2,:)=(/1d0,0d0,0d0/)
		latt%layer=1
		call latt%gen_latt()
		call latt%gen_neb(size(t))
		call latt%gen_bond(size(t))
		brizon%n1=96
		brizon%n2=24
		call latt%gen_brizon(brizon)
		call check_lattice(101)
		write(*,*)"Total site number is: ",latt%Ns

		imp=145
		r=latt%i2r(imp,:)+(/0.5d0,0.5d0,0d0/)
		allocate(bd0(latt%Ns))
		do i=1,size(bd0)
			bd0(i)=i
			bd0(i)=latt%bond(0)%bd(i)%r(2)
			!bd0(k)=(sum((latt%bond(0)%bd(k)%r-r)**2))+0.1d0*abs(theta(abs(latt%bond(0)%bd(k)%r-r))-pi/4d0)
		enddo

		allocate(bd1(size(latt%bond(1)%bd)))
		do i=1,size(bd1)
			bd1(i)=i
			bd1(i)=latt%bond(1)%bd(i)%r(2)
			!bd1(k)=(sum((latt%bond(1)%bd(k)%r-r)**2))+0.1d0*abs(theta(abs(latt%bond(1)%bd(k)%r-r))-pi/4d0)
		enddo

		! cp
		call gen_var(sg=1,nb=0)
		var(iv(0))%val(:)=0d0
		var(iv(0))%bd_sg=-1d0

		!! impure
		!call gen_var(n,(/1d0/),sg=-4,nb=0,V=1d0,var=tmp)
		!tmp(iv(0))%n=imp
		!tmp(iv(0))%val=Vimp
		!tmp(iv(0))%bd_sg=1d0

		! d-wave sc
		call gen_var(sg=2,nb=1,val=bd1,V=V)
		!call gen_var(sg=2,nb=1,V=V)
		do i=1,size(var(iv(0))%bd_sg)
			var(iv(0))%bd_sg(i)=dwave(i)
		enddo
		do l2=1,size(var(iv(0))%val)
			call random_number(var(iv(0))%val(l2))
			var(iv(0))%val(l2)=(var(iv(0))%val(l2)-0.5d0)*0.1d0+5d-1
			var(iv(0))%val(l2)=sin(sum(2d0*q*(latt%bond(var(iv(0))%nb)%bd(var(iv(0))%v2i(l2)%i(1))%r+(/0d0,0.5d0,0d0/))))+6d-2
			!var(iv(0))%val(l2)=0.7d-1
		enddo

		! sdw
		call gen_var(sg=4,nb=0,val=bd0,V=-U)
		!call gen_var(sg=4,nb=0,V=-U)
		do i=1,size(var(iv(0))%bd_sg)
			var(iv(0))%bd_sg(i)=ab(i)
		enddo
		do l2=1,size(var(iv(0))%val)
			call random_number(var(iv(0))%val(l2))
			var(iv(0))%val(l2)=(var(iv(0))%val(l2)-0.5d0)*0.1d0
			var(iv(0))%val(l2)=cos(sum(q*(latt%bond(var(iv(0))%nb)%bd(var(iv(0))%v2i(l2)%i(1))%r+(/0d0,0.5d0,0d0/))))
			!var(iv(0))%val(l2)=0.5d-1
		enddo


		! on site cdw
		call gen_var(sg=3,nb=0,val=bd0,V=U)
		!call gen_var(sg=3,nb=0,V=U)
		do i=1,size(var(iv(0))%bd_sg)
			var(iv(0))%bd_sg(i)=0.5d0
		enddo
		do l2=1,size(var(iv(0))%val)
			call random_number(var(iv(0))%val(l2))
			var(iv(0))%val(l2)=(var(iv(0))%val(l2)-0.5d0)*0.1d0+nf
			var(iv(0))%val(l2)=sin(sum(2d0*q*(latt%bond(var(iv(0))%nb)%bd(var(iv(0))%v2i(l2)%i(1))%r+(/0d0,0.5d0,0d0/))))+nf
			!var(iv(0))%val(l2)=nf+0.0001d0
		enddo

		! hp
		do l=1,size(t)
			call gen_var(sg=-3,nb=l)
			var(iv(0))%bd_sg=-1d0
			var(iv(0))%val(:)=t(l)
		enddo

		call var_shrink()
		call export_data(10)
	end subroutine
	subroutine self_consist()
		integer :: info,i
		real(8) :: x(sum(var(1:)%n)),v(size(x)),wa(nint((size(x)*(3*size(x)+13))/2.)+10),fE,err
		!complex(8) :: H(latt%Ns*spin,latt%Ns*spin)
		!real(8) :: E(size(H,1))
		common fE,err
		x=var(1:)%put()
		info=1
		!call hybrd1(do_var,size(x(:)),x(:),v(:),1d-7,info,wa,size(wa))
		!var(1:)%val=x
		!return

		do 
			if(var(1)%sg==1) then
				call hybrd1(MF_var,size(x(1:1)),x(1:1),v(1:1),1d-7,info,wa(:16),16)
				x=var(1:)%put()
			endif
			write(*,"(A$)")"****"
			call MF_var(size(x),x,v,info=-1)

			call export_data(10)
			write(20,"(es17.9$)")fE,err,var(1)%val(1)
			write(20,"(x)")
			write(*,"(es12.4$)")fE,err,var(1:)%val(1)
			write(*,"(x)")

			if(err<1d-7) then
				!call hybrd1(MF_var,size(x(:)),x(:),v(:),1d-7,info,wa,size(wa))
				!call var(1:)%get(x)

				write(*,"(es12.4$)")fE,err,var(1:)%val(1)
				write(*,"(x)")
				write(20,"(es17.9$)")fE,err,var(1)%val(1)
				write(20,"(x)")
				call export_data(10)

				exit
			endif
			call var(1:)%get(x)

		enddo
	end subroutine
	subroutine export_data(ut)
		integer :: ut,l1,i
		do l1=2,size(var(1:))
			do i=1,size(var(l1)%bd_sg)
				write(ut,"(es17.9$)")latt%bond(var(l1)%nb)%bd(i)%r,&
					latt%bond(var(l1)%nb)%bd(i)%dr,&
					var(l1)%val(var(l1)%i2v(i)),&
					var(l1)%bd_sg(i)
				write(ut,"(x)")
			enddo
			write(ut,"(x/)")
		enddo
	end subroutine
end module
program main
	use M_bdg
	implicit none
	logical :: f
	integer :: i(2)
	type(t_brizon) :: brizon_
	nf=0.85d0
	Tk=1d-5
	open(unit=101,file='../data/lattice.dat')
	f=openfile(unit=10,file='../data/order.dat')
	f=openfile(unit=20,file='../data/fenergy.dat')
	f=openfile(unit=30,file='../data/band.dat')
	f=openfile(unit=40,file='../data/fermis.dat')
	f=openfile(unit=50,file='../data/DOS.dat')
	f=openfile(unit=60,file='../data/EDC.dat')
	f=openfile(unit=70,file='../data/map_band.dat')
	call initial()
	call self_consist()
	!call band_e(70,0.01d0,(/0d0,0d0,0d0/),(/pi,0d0,0d0/),128,(/-2d0,2d0/),512)
	!call band_e(70,0.01d0,(/pi,0d0,0d0/),(/pi,pi,0d0/),128,(/-2d0,2d0/),512)
	!call band_e(70,0.01d0,(/pi,pi,0d0/),(/0d0,0d0,0d0/),128,(/-2d0,2d0/),512)
	call band(30,(/0d0,0d0,0d0/),(/pi,0d0,0d0/),128)
	call band(30,(/pi,0d0,0d0/),(/pi,pi,0d0/),128)
	call band(30,(/pi,pi,0d0/),(/0d0,0d0,0d0/),128)
	call DOS(50,0.01d0,(/-1d0,1d0/),512)
	stop

	call latt%gen_origin_brizon((/1d0,0d0,0d0/),(/0d0,1d0,0d0/),o_brizon)
	call fermis(40,0.01d0,o_brizon%k,1d0)
	call fermis(40,0.01d0,o_brizon%k,0.5d0)
	call fermis(40,0.01d0,o_brizon%k,0d0)
	call fermis(40,0.01d0,o_brizon%k,-0.5d0)
	call fermis(40,0.01d0,o_brizon%k,-1d0)
end program
