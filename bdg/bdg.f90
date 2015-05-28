module M_pmt
	use M_const
	use M_lattice
	implicit none
	real(8), parameter :: t(3)=(/1d0,-0.25d0,0d0/),nf=0.85d0,U=2.44d0,DJ=0d0,V=-1d0,cvg=1d-5,Tk=1d-5,Vimp=1d0
end module
module M_bdg
	use M_rd
	use M_utility
	use M_pmt
	use lapack95, only: heevd
	implicit none
contains
	subroutine bdg(dt,dd,ub,m,n)
		complex(8) :: dt(:),dtp(size(dt)),H(2*Ns,2*Ns)
		real(8) :: n(Ns,2),np(Ns,2),n1,wide,cp,E(2*Ns),step,pn1,dd(:),ddp(size(dd)),ub(:),ubp(size(ub)),m(Ns),mp(size(m))
		integer :: sg
		cp=0d0
		do 
			step=0d0
			pn1=0d0
			do
				call EU(dt,dd,ub,m,n,cp,H,E)
				call phy(H,E,np,dtp,ddp,ubp,mp)
				n1=sum(np)/Ns
				if(abs(n1-nf)<cvg*5d0) then
					exit
				endif
				if(step<1d-10) then
					step=abs(n1-nf)
				endif
				call find_cross(pn1,n1-nf,sg)
				if(sg/=0) then
					step=max(step*0.3d0,1d-8)
				endif
				write(*,*)n1,cp,step*sign(1d0,nf-n1)
				cp=cp-step*sign(1d0,n1-nf)
			enddo
			call export_data(10,dt,m,n)
			if(sum(abs(dt-dtp))/Ns<cvg.and.sum(abs(m-mp))/Ns<cvg.and.sum(abs(dd-ddp))/Ns<cvg) then
				write(20,"(e12.4)")E
				exit
			endif
			dt=dtp
			dd=ddp
			n=np
			ub=ubp
			m=mp
			!write(*,"(e12.4$)")real(dt(1))/(DJ-V),m(1)/(4d0*DJ),dd(1)/(0.5d0*DJ+V),ub(1)/(0.5d0*DJ+V)
			write(*,"(e12.4$)")real(dt(1)),m(1),dd(1),ub(1)
			write(*,"(1X)")
		enddo 
	end subroutine
	subroutine EU(dt,dd,ub,m,ne,cp,H,E)
		complex(8) :: H(:,:),dt(:)
		real(8) :: ne(:,:),cp,E(:),dd(:),ub(:),m(:),dp
		integer :: i,j,k,n,info
		H=0d0
		dp=abs(1d0-nf)
		dp=1d0
		! on site
		do i=1,Ns
			!H(i,i)=H(i,i)+m(i)*0.5d0-cp
			!H(i+Ns,i+Ns)=H(i+Ns,i+Ns)+m(i)*0.5d0+cp
			H(i,i)=H(i,i)+U*ne(i,2)-cp
			H(i+Ns,i+Ns)=H(i+Ns,i+Ns)-U*ne(i,1)+cp
			! impurity
			if(i==nint(sqrt(real(Ns))/2d0)+Ns/2) then
				H(i,i)=H(i,i)+Vimp
				H(i+Ns,i+Ns)=H(i+Ns,i+Ns)-Vimp
			endif
		enddo
		! 1st near bond
		do k=1,size(bond(1)%bd)
			i=bond(1)%bd(k)%i(1)
			j=bond(1)%bd(k)%i(2)
			! uniform bond order
			H(i,j)=H(i,j)-t(1)-ub(k)
			H(j,i)=H(j,i)-t(1)-ub(k)
			H(i+Ns,j+Ns)=H(i+Ns,j+Ns)+t(1)+ub(k)
			H(j+Ns,i+Ns)=H(j+Ns,i+Ns)+t(1)+ub(k)
			! d-wave pairs
			H(i,j+Ns)=H(i,j+Ns)+dt(k)
			H(j,i+Ns)=H(j,i+Ns)+dt(k)
			H(i+Ns,j)=H(i+Ns,j)+conjg(dt(k))
			H(j+Ns,i)=H(j+Ns,i)+conjg(dt(k))
			! DDW
			!H(i,j)=H(i,j)+dd(k)*img
			!H(j,i)=H(j,i)-dd(k)*img
			!H(i+Ns,j+Ns)=H(i+Ns,j+Ns)+dd(k)*img
			!H(j+Ns,i+Ns)=H(j+Ns,i+Ns)-dd(k)*img
		enddo
		! 2nd near bond
		do k=1,size(bond(2)%bd)
			i=bond(2)%bd(k)%i(1)
			j=bond(2)%bd(k)%i(2)
			H(i,j)=H(i,j)-t(2)
			H(j,i)=H(j,i)-t(2)
			H(i+Ns,j+Ns)=H(i+Ns,j+Ns)+t(2)
			H(j+Ns,i+Ns)=H(j+Ns,i+Ns)+t(2)
		enddo
		! 3nd near bond
		do k=1,size(bond(3)%bd)
			i=bond(3)%bd(k)%i(1)
			j=bond(3)%bd(k)%i(2)
			H(i,j)=H(i,j)-t(3)
			H(j,i)=H(j,i)-t(3)
			H(i+Ns,j+Ns)=H(i+Ns,j+Ns)+t(3)
			H(j+Ns,i+Ns)=H(j+Ns,i+Ns)+t(3)
		enddo
		call heevd(H,E,"V","U",info)
	end subroutine
	subroutine phy(H,E,ne,dt,dd,ub,m)
		real(8) :: ne(:,:),E(:),f(size(E)),th(size(E)),dd(:),ub(:),m(:)
		complex(8) :: H(:,:),dt(:)
		integer :: i,k,j,n
		ne=0d0
		dt=0d0
		dd=0d0
		ub=0d0
		m=0d0
		f=1d0/(exp(E/Tk)+1d0)
		!!$OMP PARALLEL DO REDUCTION(+:ne) COLLAPSE(2)
		do i=1,Ns
			do n=1,Ns*2
				ne(i,1)=ne(i,1)+conjg(H(i,n))*H(i,n)*f(n)
				ne(i,2)=ne(i,2)+H(i+Ns,n)*conjg(H(i+Ns,n))*(1d0-f(n))
			enddo
		enddo
		!!$OMP END PARALLEL DO
		!!$OMP PARALLEL DO REDUCTION(+:dt) PRIVATE(i,j) COLLAPSE(2)
		do k=1,size(bond(1)%bd)
			do n=1,Ns*2
				i=bond(1)%bd(k)%i(1)
				j=bond(1)%bd(k)%i(2)
				dt(k)=dt(k)+0.5d0*(H(i,n)*conjg(H(j+Ns,n))+H(j,n)*conjg(H(i+Ns,n)))*(1d0-f(n))
				!dd(k)=dd(k)+0.25d0*imag(conjg(H(i,n))*H(j,n)-conjg(H(j,n))*H(i,n)&
					!+conjg(H(i+Ns,n))*H(j+Ns,n)-conjg(H(j+Ns,n))*H(i+Ns,n))*f(n)
				!ub(k)=ub(k)+0.25d0*real(conjg(H(i,n))*H(j,n)+conjg(H(j,n))*H(i,n)&
					!-conjg(H(i+Ns,n))*H(j+Ns,n)-conjg(H(j+Ns,n))*H(i+Ns,n))*f(n)
			enddo
		enddo
		!!$OMP END PARALLEL DO
		dt=(DJ-V)*dt
		!dd=(0.5d0*DJ+V)*dd
		!ub=(0.5d0*DJ+V)*ub
		do i=1,Ns
			do k=1,size(neb(i)%nb(2)%neb)
				m(i)=m(i)+DJ*0.5d0*sum(ne(neb(i)%nb(1)%neb(k),:)*(/1d0,-1d0/))
			enddo
			m(i)=m(i)+(ne(i,1)-ne(i,2))
		enddo
	end subroutine
	function is_a(i)
		logical :: is_a
		integer :: i
		if(mod(sum(nint(i2r(i,:)-i2r(1,:))),2)==0) then
			is_a=.true.
		else
			is_a=.false.
		endif
	end function
	subroutine export_data(ut,dt,m,n)
		complex(8) :: dt(:)
		real(8) :: m(:),n(:,:)
		integer :: ut,i
		rewind(ut)
		do i=1,size(bond(1)%bd)
			write(ut,"(e15.6$)")bond(1)%bd(i)%r
			if(nint(bond(1)%bd(i)%dir(2))==0) then
				write(ut,"(e15.6$)")dt(i)
			else
				write(ut,"(e15.6$)")-dt(i)
			endif
			write(ut,"(X)")
		enddo
		write(ut,"(X/)")
		do i=1,Ns
			write(ut,"(e15.6$)")i2r(i,:),n(i,1)+n(i,2)
			write(ut,"(X)")
		enddo
		write(ut,"(X/)")
		do i=1,Ns
			if(is_a(i)) then
				write(ut,"(e15.6$)")i2r(i,:),m(i)
			else
				write(ut,"(e15.6$)")i2r(i,:),-m(i)
			endif
			write(ut,"(X)")
		enddo
		write(ut,"(X/)")
	end subroutine
end module
program main
	use M_bdg
	implicit none
	integer :: i,j,sg
	real(8), allocatable :: n(:,:),m(:)
	real(8), allocatable :: dd(:),ub(:)
	complex(8), allocatable :: dt(:),tmp
	logical :: f
	f=openfile(unit=10,file='../data/order.dat')
	f=openfile(unit=20,file='../data/bdg_energy.dat')
	call init_random_seed()
	! lattice 
	a1=(/1d0,0d0/)
	a2=(/0d0,1d0/)
	T1=a1*24
	T2=a2*24
	bdc=(/1d0,1d0/)
	allocate(sub(1,2))
	sub(1,:)=(/0d0,0d0/)
	call gen_latt()
	call gen_neb()
	call gen_bond(3)
	! finish i2r(i,2),neb(i)%nb(j)%bond(k)/bdc(k)/r(k,2)
	allocate(&
		dt(size(bond(1)%bd)),&
		dd(size(bond(1)%bd)),&
		ub(size(bond(1)%bd)),&
		n(Ns,2),&
		m(Ns)&
	)
	do i=1,size(bond(1)%bd)
		if(nint(bond(1)%bd(i)%dir(2))==0) then
			dt(i)=0.1d0
			!if(is_a(bond(1)%bd(i)%i(1))) then
				!dd(i)=0.1d0
			!else
				!dd(i)=-0.1d0
			!endif
		else
			dt(i)=-0.1d0
			!if(is_a(bond(1)%bd(i)%i(1))) then
				!dd(i)=-0.1d0
			!else
				!dd(i)=0.1d0
			!endif
		endif
		ub(i)=0d0
	enddo
	do i=1,Ns
		if(is_a(i)) then
			m(i)=0.1d0
		else
			m(i)=-0.1d0
		endif
		n(i,1)=(nf-m(i))/2d0
		n(i,2)=(nf+m(i))/2d0
	enddo
	n=nf/2d0
	call bdg(dt,dd,ub,m,n)
	!export data
	call export_data(10,dt,m,n)
end program
