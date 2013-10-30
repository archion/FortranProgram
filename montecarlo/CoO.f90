module global
	implicit none
	integer, parameter :: dn(2)=(/36,36/),dn2=dn(1)*dn(2),ne=dn2/3
	integer :: latt(dn2,4,3),pbcx(-dn(1)+1:2*dn(1)),pbcy(-dn(2)+1:2*dn(2))
	real(8), parameter :: V(3)=(/4d0,-2d0,2d0/),pi=3.14159265359d0
	complex(8), parameter :: img=(0d0,1d0)
	contains
		subroutine gen_latt_square()
!					  j
!			 1 --  2 --  3 --  4
!			 |     |     |     |           4      3   4
!			 5 --  6 --  7 --  8           |       \ /
!		i	 |     |     |     |        3--0--1     0
!			 9 -- 10 -- 11 -- 12           |       / \
!			 |     |     |     |           2      2   1
!			13 -- 14 -- 15 -- 16
			implicit none
			integer :: ii,i,j
			pbcx=(/(i,i=1,dn(1)),(i,i=1,dn(1)),(i,i=1,dn(1))/)
			pbcy=(/(i,i=1,dn(2)),(i,i=1,dn(2)),(i,i=1,dn(2))/)
			do ii=1,dn2
				j=mod(ii-1,dn(1))+1
				i=(ii-1)/dn(1)+1
				latt(ii,1,1)=dn(1)*(i-1)+pbcx(j+1)
				latt(ii,2,1)=dn(1)*(pbcy(i+1)-1)+j
				latt(ii,3,1)=dn(1)*(i-1)+pbcx(j-1)
				latt(ii,4,1)=dn(1)*(pbcy(i-1)-1)+j
				latt(ii,1,2)=dn(1)*(pbcy(i+1)-1)+pbcx(j+1)
				latt(ii,2,2)=dn(1)*(pbcy(i+1)-1)+pbcx(j-1)
				latt(ii,3,2)=dn(1)*(pbcy(i-1)-1)+pbcx(j-1)
				latt(ii,4,2)=dn(1)*(pbcy(i-1)-1)+pbcx(j+1)
				latt(ii,1,3)=dn(1)*(i-1)+pbcx(j+2)
				latt(ii,2,3)=dn(1)*(pbcy(i+2)-1)+j
				latt(ii,3,3)=dn(1)*(i-1)+pbcx(j-2)
				latt(ii,4,3)=dn(1)*(pbcy(i-2)-1)+j
			enddo
		end subroutine gen_latt_square
end module
module phyval
	use global
	implicit none
	contains
		subroutine strfact(cfg,S)
			implicit none
			integer :: qi,qj,i,j,ip,jp,cfg(dn(1),dn(2)),n(dn(1),dn(2))
			complex(8) :: S(dn(1),dn(2)),e
			S=0d0
			n=0
			do i=1,dn(2)
				do j=1,dn(1)
					do ip=1,dn(2)
						do jp=1,dn(1)
							n(i,j)=n(i,j)+(1-cfg(ip,jp))*(1-cfg(pbcy(i+ip),pbcx(j+jp)))
						enddo
					enddo
				enddo
			enddo
			do qi=1,dn(2)
				do qj=1,dn(1)
					do i=1,dn(2)
						do j=1,dn(1)
							e=exp(2d0*pi*img*(1d0/dn(1)*qj*j+1d0/dn(2)*qi*i))
							S(qi,qj)=S(qi,qj)+n(i,j)*e
						enddo
					enddo
				enddo
			enddo
			S=S/(dn2-ne)
		end subroutine
end module
module MonteCarlo
	use global
	implicit none
	contains
		subroutine MC_init(A,m,cdt)
			implicit none
			integer :: m,A(dn2)
			character :: cdt
			call gen_latt_square()
			select case(cdt)
			case("r")
				call random_init(A,m)
			case("s")
				call stripe_init(A,m)
			case default
				write(*,*)"Failed to initialize: unknown pattern"
				stop
			end select
		end subroutine
		subroutine MC_step(bt,cfg)
			implicit none
			integer :: i,j,tmp,cfg(dn2)
			real(8) :: de,bt,r,rn
			call random_number(rn)
			i=1+int((dn2-1)*rn)
			call random_number(rn)
			j=1+int((dn2)*rn)
			if(j==i) then
				j=dn2
			endif
			if(cfg(i)/=cfg(j)) then
				call energydiff(cfg,de,i,j)
				r=exp(bt*de)
				call random_number(rn)
				if(r>1d0.or.rn<r) then
					tmp=cfg(i)
					cfg(i)=cfg(j)
					cfg(j)=tmp
				endif
			endif
		end subroutine
		subroutine energydiff(cfg,de,i,j)
			implicit none
			integer :: i,j,k,l,m,tmp,cfg(dn2)
			real(8) :: de
			de=0d0
			m=(-1)**cfg(i)
			do l=1,size(V)
				tmp=0
				do k=1,4
					tmp=tmp+(cfg(latt(j,k,l))-cfg(latt(i,k,l)))
				enddo
				de=de+V(l)*tmp
			enddo
			de=de*m
		end subroutine
		subroutine random_init(A,m)
			implicit none
			integer :: m,i,ir,A(dn2)
			real(8) :: rn
			A=0
			do i=dn2-m+1,dn2
				call random_number(rn)
				ir=1+int(rn*(i))
				if(A(ir)==0) then
					A(ir)=1
				else
					A(i)=1
				endif
			enddo
		end subroutine
		subroutine stripe_init(A,m)
			implicit none
			integer :: m,i,j,A(dn2),S(m*(dn(1)-1)/dn2+1)
			S=(/(i,i=1,dn(1),3)/)
			A=0
			do i=1,dn(2)
				do j=1,size(S)
					S(j)=pbcx(S(j)-1)
					A((i-1)*dn(2)+S(j))=1
				enddo
			enddo
		end subroutine
		subroutine init_random_seed()
			integer :: i, n, clock
			integer, dimension(:), allocatable :: seed
			call random_seed(size = n)
			allocate(seed(n))
			call system_clock(count=clock)
			seed = clock + 37 * (/ (i - 1, i = 1, n) /)
			call random_seed(put = seed)
			deallocate(seed)
		end subroutine
end module
program main
	use global
	use MonteCarlo
	use phyval
	implicit none
	integer :: i,j,Nmax=2000000,Nhot=10000,cfg(dn2),ret(5)=(/25000,100000,200000,500000,2000000/)
	real(8) :: bt=1e8
	complex(8) :: S(dn2),Savg(dn2)
	open(unit=10,file="../data/montecarlo_i.dat")
	open(unit=20,file="../data/pattern_i.dat")
	write(10,"(A)")"set term pngcairo size 300,300"
	write(10,"(A)")"set output 'montecarlo.png'"
	write(10,"(A)")'set title "step=0"'
	write(10,"(A)")"set pm3d map"
	write(10,"(A)")"set pm3d interpolate 0,0"
	write(10,"(A)")"set size square"
	write(10,"(A)")"set palette rgbformulae 22,13,-31"
	write(10,"(A)")"splot '-' u 1"
	write(10,"('#beta=',e12.4,'hot step=',I7,'particl number=',I5)")bt,Nhot,ne
	write(20,"(A)")"set term pngcairo size 300,300"
	write(20,"(A)")"set output 'pattern.png'"
	write(20,"(A)")"set pm3d map"
	write(20,"(A)")"set size square"
	write(20,"(A)")"set palette rgbformulae 22,13,-31"
	write(20,"(A)")"splot '-' u 1"
	!Initialization
	Savg=0d0
	call init_random_seed()
	call MC_init(cfg,ne,"r")
	!write(*,"(36I2)")cfg
	!call strfact(cfg,S)
	!do i=1,dn(2)
		!write(10,"(2e13.4)")S((i-1)*dn(1)+1:i*dn(1))
		!write(10,"(1X)")
		!write(20,"(I2)")cfg((i-1)*dn(1)+1:i*dn(1))
		!write(20,"(1X)")
	!enddo
	!stop
	!Initialization finish
	do i=1,Nmax
		call MC_step(bt,cfg)
		if(mod(i,1000)==0) then
			write(*,*)"Monte Carlo Step: ",i
		endif
		if(i>Nhot) then
			call strfact(cfg,S)
			Savg=Savg+S
		endif
		if(any(ret==i)) then
			write(10,"('#Monte Carlo Step: ',I8)")i
			write(20,"('#Monte Carlo Step: ',I8)")i
			do j=1,dn(2)
				write(10,"(2e13.4)")Savg((j-1)*dn(1)+1:j*dn(1))/(i-Nhot)
				write(10,"(1X)")
				write(20,"(I2)")cfg((j-1)*dn(1)+1:j*dn(1))
				write(20,"(1X)")
			enddo
			write(10,"(1X)")
			write(20,"(1X)")
		endif
	enddo
end

	


