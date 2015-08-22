module pmt
	use M_utility
	use M_const
	use M_hamilton
	implicit none
	real(8), parameter :: t(2)=(/1d0,1d0/3d0/),phi=pi/3d0,dm=1d0
contains
	subroutine initial()
		integer :: iv(-1:1),n,i,j,l,k
		type(t_var), allocatable :: tmp(:)
		allocate(tmp(-10000:10000))
		call init_random_seed()

		! lattice 
		latt%a1=(/sqrt(3d0),0d0,0d0/)
		latt%a2=(/sqrt(3d0)/2d0,1.5d0,0d0/)
		latt%T1=(/0d0,3d0,0d0/)
		latt%T2=(/30d0*sqrt(3d0),0d0,0d0/)
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
		call latt%gen_brizon()
		spin=1
		call check_lattice(30)
		write(*,*)"Total site number is: ",latt%Ns

		iv=(/1,0,0/)

		! cp
		call gen_var(iv,n,sg=-1,nb=0,V=1d0,var=tmp)
		do i=iv(0)-n+sign(1,n),iv(0),sign(1,n)
			tmp(i)%bd_sg=-1d0
			tmp(i)%val=1d0
		enddo

		! hp
		call gen_var(iv,n,sg=-3,nb=1,V=1d0,var=tmp)
		do i=iv(0)-n+sign(1,n),iv(0),sign(1,n)
			tmp(i)%bd_sg=-1d0
			tmp(i)%val=t(1)
		enddo

		call gen_var(iv,n,sg=-3,nb=2,V=1d0,var=tmp)
		do i=iv(0)-n+sign(1,n),iv(0),sign(1,n)
			do k=1,size(tmp(i)%n)
				tmp(i)%bd_sg(k)=exp(-1d0*ab(latt%bond(tmp(i)%nb)%bd(tmp(i)%n(k))%i(1))&
					*phase_rule(latt%bond(tmp(i)%nb)%bd(tmp(i)%n(k))%dr)&
					*img*phi)
				!write(30,"(es17.9$)")latt%bond(tmp(i)%nb)%bd(tmp(i)%n(k))%r-sign(0.5d0,imag(tmp(i)%bd_sg(k)))*latt%bond(tmp(i)%nb)%bd(tmp(i)%n(k))%dr,sign(1d0,imag(tmp(i)%bd_sg(k)))*latt%bond(tmp(i)%nb)%bd(tmp(i)%n(k))%dr
				!write(30,"(x)")
			enddo
			tmp(i)%val=t(2)
		enddo
		!write(30,"(x/)")

		call gen_var(iv,n,sg=-3,nb=0,V=1d0,var=tmp)
		do i=iv(0)-n+sign(1,n),iv(0),sign(1,n)
			do k=1,size(tmp(i)%n)
				tmp(i)%bd_sg(k)=1d0&
					*ab(latt%bond(tmp(i)%nb)%bd(tmp(i)%n(k))%i(1))
			enddo
			tmp(i)%val=dm
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
		deallocate(tmp)
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
	call band(brizon%T(1,:),brizon%T(2,:),1024,20)
	!call band((/0d0,0d0,0d0/),brizon%T(1,:),128,20)
	!call band(brizon%T(1,:),brizon%T(2,:),128,20)
end program
