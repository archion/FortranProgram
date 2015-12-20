module pmt
	use M_utility
	use M_const
	use M_hamilton_test
	implicit none
	real(8), parameter :: t(2)=(/1d0,1d0/3d0/),phi=pi/6d0,dm=1d0
contains
	subroutine initial()
		integer :: l2,l3,l,i
		real(8), allocatable :: bd0(:),bd1(:)
		allocate(var(-10000:10000))
		call init_random_seed()

		! lattice 
		latt%a1=(/sqrt(3d0),0d0,0d0/)
		latt%a2=(/sqrt(3d0)/2d0,1.5d0,0d0/)
		!latt%T1=(/0d0,3d0,0d0/)
		!latt%T2=(/30.5d0*sqrt(3d0),0d0,0d0/)
		latt%T1=(/sqrt(3d0),0d0,0d0/)
		latt%T2=(/0d0,31d0,0d0/)
		latt%bdc(1)=1d0
		latt%bdc(2)=0d0
		allocate(latt%sub(2,3))
		latt%sub(1,:)=(/0d0,0d0,0d0/)
		latt%sub(2,:)=(/0.5d0*sqrt(3d0),0.5d0,0d0/)
		latt%layer=1
		call latt%gen_latt()
		call latt%gen_neb(size(t))
		call latt%gen_bond(size(t))
		brizon%n1=32
		brizon%n2=32
		call latt%gen_brizon(brizon)
		spin=1
		call check_lattice(30)
		write(*,*)"Total site number is: ",latt%Ns


		! cp
		call gen_var(sg=1,nb=0)
		var(iv(0))%bd=-1d0
		var(iv(0))%val=0d0

		! hp
		call gen_var(sg=-3,nb=1)
		var(iv(0))%bd=-1d0
		var(iv(0))%val=t(1)

		call gen_var(sg=-3,nb=2)
		do i=1,size(var(iv(0))%bd)
			var(iv(0))%bd(i)=exp(-1d0*ab(latt%bond(var(iv(0))%nb)%bd(i)%i(1))&
				*phase_rule(latt%bond(var(iv(0))%nb)%bd(i)%dr)&
				*img*phi)
			!write(30,"(es17.9$)")latt%bond(tmp(i)%nb)%bd(tmp(i)%n(k))%r-sign(0.5d0,imag(tmp(i)%bd_sg(k)))*latt%bond(tmp(i)%nb)%bd(tmp(i)%n(k))%dr,sign(1d0,imag(tmp(i)%bd_sg(k)))*latt%bond(tmp(i)%nb)%bd(tmp(i)%n(k))%dr
			!write(30,"(x)")
		enddo
		var(iv(0))%val=t(2)
		!write(30,"(x/)")

		call gen_var(sg=-3,nb=0)
		do i=1,size(var(iv(0))%bd)
			var(iv(0))%bd(i)=ab(latt%bond(var(iv(0))%nb)%bd(i)%i(1))
		enddo
		var(iv(0))%val=dm

		call var_shrink()

	end subroutine
	function phase_rule(dr)
		real(8) :: dr(:),er,phase_rule
		integer :: i
		phase_rule=1d0
		do i=0,2
			if(0d0-1d-6+pi*2d0/3d0*i<theta(dr).and.theta(dr)<=pi/3d0-1d-6+pi*2d0/3d0*i) then
				phase_rule=-1d0
			endif
		enddo
	end function
end module

program main
	use pmt
	implicit none
	logical :: f
	f=openfile(unit=10,file='../data/energy.dat')
	f=openfile(unit=20,file='../data/band.dat')
	f=openfile(unit=30,file='../data/lattice.dat')
	call initial()
	!call energy(10)
	call band(20,brizon%T(1,:),brizon%T(2,:),512)
	!call band((/0d0,0d0,0d0/),brizon%T(1,:),128,20)
	!call band(brizon%T(1,:),brizon%T(2,:),128,20)
end program
