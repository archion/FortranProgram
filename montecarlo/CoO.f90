module global
	use M_const
	use M_utility
	use M_latt, only : &
		 latt         =>  square         ,& 
		 latt_bc      =>  square_bc      ,& 
		 latt_one2two =>  square_one2two ,& 
		 latt_two2one =>  square_two2one 
	implicit none
	integer, parameter :: Ns(2)=(/36,36/),Ns2=Ns(1)**2,Tx(2)=(/Ns(1),0/),Ty(2)=(/0,Ns(1)/),vn=6,ne=Ns2/3
	integer :: neb(Ns2,4,3)
	real(8), parameter :: V(3)=(/4d0,-2d0,2d0/)
	complex(8), parameter :: ipiN=2d0*pi*img/Ns(1)
end module
module phyval
	use M_fft
	use global
	implicit none
	contains
		subroutine strfact(cfg,S)
			integer :: qi,qj,i,j,cfg(Ns(1),Ns(2))
			complex(8) :: S(0:Ns(1)-1,0:Ns(2)-1),ipi=2d0*pi*img,q(2)
			real(8) :: iNs(2)
			S=0d0
			!S=1d0-cfg
			!call fft2d(S,1)
			!$OMP PARALLEL DO REDUCTION(+:S) PRIVATE(q)
			do qi=0,Ns(2)-1
				do qj=0,Ns(1)-1
					q=ipiN*(/qi,qj/)
					do i=1,Ns(2)
						do j=1,Ns(1)
							S(qi,qj)=S(qi,qj)+(1-cfg(i,j))*exp(sum(q*(/j,i/)))
						enddo
					enddo
				enddo
			enddo
			!$OMP END PARALLEL DO
			S=S*conjg(S)/(Ns2-ne)
			S(0,0)=0d0
		end subroutine
		!subroutine strfact(cfg,S)
			!implicit none
			!integer :: qi,qj,i,j,ip,jp,cfg(dn(1),dn(2)),n(dn(1),dn(2))
			!complex(8) :: S(dn(1),dn(2)),e
			!S=0d0
			!n=0
			!do i=1,Nn(2)
				!do j=1,dn(1)
					!do ip=1,dn(2)
						!do jp=1,dn(1)
							!n(i,j)=n(i,j)+(1-cfg(ip,jp))*(1-cfg(pbcy(i+ip),pbcx(j+jp)))
						!enddo
					!enddo
				!enddo
			!enddo
			!do qi=1,dn(2)
				!do qj=1,dn(1)
					!do i=1,dn(2)
						!do j=1,dn(1)
							!e=exp(2d0*pi*img*(1d0/dn(1)*qj*j+1d0/dn(2)*qi*i))
							!S(qi,qj)=S(qi,qj)+n(i,j)*e
						!enddo
					!enddo
				!enddo
			!enddo
			!S=S/(Ns2-ne)
		!end subroutine
end module
module MonteCarlo
	use global
	implicit none
	contains
		subroutine MC_init(A,m,cdt)
			implicit none
			integer :: m,A(Ns2)
			character :: cdt
			call latt(Ns,Tx,Ty,neb)
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
			integer :: i,j,tmp,cfg(Ns2)
			real(8) :: de,bt,r,rn
			call random_number(rn)
			i=1+int((Ns2-1)*rn)
			call random_number(rn)
			j=1+int((Ns2)*rn)
			if(j==i) then
				j=Ns2
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
			integer :: i,j,k,l,m,tmp,cfg(Ns2)
			real(8) :: de
			de=0d0
			m=(-1)**cfg(i)
			do l=1,size(V)
				tmp=0
				do k=1,4
					tmp=tmp+(cfg(neb(j,k,l))-cfg(neb(i,k,l)))
				enddo
				de=de+V(l)*tmp
			enddo
			de=de*m
		end subroutine
		subroutine random_init(A,m)
			implicit none
			integer :: m,i,ir,A(Ns2)
			real(8) :: rn
			A=0
			do i=Ns2-m+1,Ns2
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
			integer :: m,i,j,A(Ns2),S(m*(Ns(1)-1)/Ns2+1)
			S=(/(i,i=1,Ns(1),3)/)
			A=0
			do i=1,Ns(2)
				do j=1,size(S)
					!S(j)=pbcx(S(j)-1)
					A((i-1)*Ns(2)+S(j))=1
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
	integer :: i,j,Nmax=200000000,Nhot=10000,cfg(Ns2),ret(5)=(/25000,100000,500000,2000000,200000000/)
	real(8) :: bt=1
	complex(8) :: S(Ns2),Savg(Ns2)
	open(unit=10,file="../data/montecarlo.dat")
	open(unit=20,file="../data/pattern.dat")
	!Initialization
	Savg=0d0
	call init_random_seed()
	call MC_init(cfg,ne,"r")
	call strfact(cfg,S)
	write(10,"(36e13.4)")real(S)
	write(20,"(36I2)")cfg
	write(10,"(1X/)")
	write(20,"(1X/)")
	!Initialization finish
	do i=1,Nmax
		call MC_step(bt,cfg)
		if(mod(i,1000)==0) then
			write(*,*)"Monte Carlo Step: ",i
		endif
		if(any((i+1000)>ret.and.i<ret)) then
			call strfact(cfg,S)
			Savg=Savg+S
		endif
		if(any(ret==i)) then
			write(10,"('#Monte Carlo Step: ',I8)")i
			write(20,"('#Monte Carlo Step: ',I8)")i
			write(10,"(36e13.4)")real(Savg/1000)
			write(20,"(36I2)")cfg
			!write(10,"(32e13.4)")real(Savg/(i-Nhot))
			!write(20,"(32I2)")cfg
			write(10,"(1X/)")
			write(20,"(1X/)")
			Savg=0d0
		endif
	enddo
end

	


