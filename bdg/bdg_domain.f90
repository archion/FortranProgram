module global
	implicit none
	integer, parameter :: dn(2)=(/16,16/),dn2=dn(1)*dn(2),imp=dn2/2+dn(1)/2
	integer :: latt(dn2,4,3),fg=1
	real(8), parameter :: t(3)=(/1d0,-0.3d0,0d0/),nf=0.875d0,V=0d0,Vimp=0d0,pi=3.14159265359d0,cvg=1e-4,bt=1e5
	real(8) :: U=2.44d0
	complex(8), parameter :: img=(0d0,1d0),ipi=2d0*pi*img
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
		integer :: ii,i,j,pbcx(-dn(1)+1:2*dn(1)),pbcy(-dn(2)+1:2*dn(2))
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
program main
	use global
	implicit none
	integer :: i,j,k,sg
	real(8) :: n(dn2,2),te
	complex(8) :: dt(dn2,4)
	open(unit=10,file='../data/order.dat')
	open(unit=20,file='../data/energy.dat')
	call gnuplot
	!!$OMP PARALLEL DO PRIVATE(U,dt,n,te,fg) SCHEDULE(GUIDED)
	do i=6,16,2
		U=i*0.5d0
		write(*,"(e15.6,$)")U
		write(20,"(e15.6,$)")U
		fg=1
		call inital(dt,n)
		call bdg(dt,n,te)
		write(20,"(e15.6,$)")te
		fg=0
		call inital(dt,n)
		call bdg(dt,n,te)
		write(20,"(e15.6)")te
		call exportdata(dt,n)
	enddo
	!!$OMP END PARALLEL DO
	!export data
	close(10)
end
subroutine exportdata(dt,n)
	use global
	implicit none
	integer :: i,sg
	complex(8) :: dt(dn2,4)
	real(8) :: n(dn2,2)
	do i=1,dn2
		write(10,"(e15.6,$)")abs(dt(i,1)+dt(i,3)-dt(i,2)-dt(i,4))/4.0d0
		if(mod(i,dn(1))==0) then
			write(10,"(X)")
		endif
	enddo
	write(10,"(X)")
	do i=1,dn2
		write(10,"(e15.6,$)")n(i,1)+n(i,2)
		if(mod(i,dn(1))==0) then
			write(10,"(X)")
		endif
	enddo
	write(10,"(X)")
	sg=-1
	do i=1,dn2
		sg=-sg
		write(10,"(e15.6,$)")sg*0.5d0*(n(i,1)-n(i,2))
		if(mod(i,dn(1))==0) then
			sg=-sg
			write(10,"(X)")
		endif
	enddo
end
subroutine inital(dt,n)
	use global
	implicit none
	complex(8) :: dt(dn2,4)
	real(8) :: n(dn2,2)
	integer :: i,sg
	logical :: lc
	call init_random_seed()
	!call random_number(rdom)
	do i=1,dn2
		dt(i,:)=(/0.0d0,-0.0d0,0.0d0,-0.0d0/)
	enddo
	sg=-1
	do i=1,dn2
		sg=-sg
		n(i,1)=sg*0.1+0.5d0
		n(i,2)=-sg*0.1+0.5d0
		if(fg==0) then
			lc=mod(mod(i-1,dn(2)),8)==0
		else
			lc=mod((i-1)/dn(2)-mod(i-1,dn(2)),8)==0
		endif
		if(mod(i,dn(1))==0.or.lc) then
			sg=-sg
			n(i,:)=0d0
		endif
	enddo
	call gen_latt_square()
end
subroutine bdg(dt,n,te)
	use global
	implicit none
	complex(8) :: dt(dn2,4),dtp(dn2,4),H(2*dn2,2*dn2)
	real(8) :: n(dn2,2),np(dn2,2),n1,al,sp,pn1,E(2*dn2),te
	logical :: flaga,flagb
	do 
		pn1=nf+0.1d0
		al=1d0
		do
			call EU(dt,n,sp,H,E)
			call order(H,E,np,dtp)
			n1=1.0/dn2*sum(np)
			!write(*,*)n1,sp
			if((nf-n1)/(nf-pn1)<-0.8d0) then
				al=al*0.7d0
			endif
			if(abs(n1-nf)<cvg) then
				exit
			else
				sp=sp+al*(nf-n1)
			endif
			pn1=n1
		enddo
		if((sum(abs(dt-dtp))/dn2<cvg.and.sum(abs(np(:,1)-n(:,1)))/dn2<cvg)) then
			exit
		endif
		dt=dtp
		n=np*0.6d0+n*0.4d0
		!write(*,*)fg,n(2,1)
	enddo 
	call totaleng(dt,n,sp,E,te)
end
subroutine EU(dt,n,sp,H,E)
	use f95_lapack, only: la_heevd
	!use lapack
	use global
	implicit none
	complex(8) :: H(2*dn2,2*dn2),dt(dn2,4)
	real(8) :: E(2*dn2),n(dn2,2),sp
	integer :: i,j,info
	logical :: lc
	H=0d0
	do i=1,dn2
		do j=1,4
			! nearest neighbor pairs
			H(i,latt(i,j,1))=-t(1)
			H(i,latt(i,j,1)+dn2)=dt(i,j)
			H(i+dn2,latt(i,j,1)+dn2)=t(1)
			H(i+dn2,latt(i,j,1))=dconjg(dt(i,j))
			! next-nearest neighbor pairs
			H(i,latt(i,j,2))=-t(2)
			H(i+dn2,latt(i,j,2)+dn2)=t(2)
		enddo
		! on site
		!if(fg==0) then
			!lc=mod(mod(i-1,dn(2)),8)==0
		!else
			!lc=mod((i-1)/dn(2)-mod(i-1,dn(2)),8)==0
		!endif
		!if(lc) then
			!H(i,i)=H(i,i)+100d0
			!H(i+dn2,i+dn2)=H(i+dn2,i+dn2)-100d0
		!else
			H(i,i)=U*n(i,2)-sp
			H(i+dn2,i+dn2)=-U*n(i,1)+sp
		!endif
	enddo
	call la_heevd(H,E,"V","U",info)
	!call czheev("V","U",H,E,info)
end
subroutine order(H,E,np,dtp)
	use global
	implicit none
	real(8) :: np(dn2,2),f(2*dn2),th(2*dn2),E(2*dn2)
	complex(8) :: dtp(dn2,4),H(2*dn2,2*dn2)
	integer :: i,j,n
	np=0d0
	dtp=0d0
	f=1d0/(exp(bt*E)+1d0)
	th=tanh(0.5d0*bt*E)
	do i=1,dn2
		do n=1,2*dn2
			np(i,1)=np(i,1)+H(i,n)*dconjg(H(i,n))*f(n)
			np(i,2)=np(i,2)+H(i+dn2,n)*dconjg(H(i+dn2,n))*(1d0-f(n))
			do j=1,4
				dtp(i,j)=dtp(i,j)+0.25d0*V*(H(i,n)*dconjg(H(latt(i,j,1)+dn2,n))+H(latt(i,j,1),n)*dconjg(H(i+dn2,n)))*th(n)
			enddo
		enddo
	enddo
end
subroutine totaleng(dt,n,sp,E,te)
	use global
	implicit none
	real(8) :: n(dn2,2),sp,E(2*dn2),te,f(2*dn2)
	complex(8) :: dt(dn2,4)
	f=1d0/(exp(bt*E)+1d0)
	te=sum(E*f)-U*sum(n(:,1)*n(:,2))+U*sum(n(:,1))-dn2*sp
end
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
subroutine gnuplot()
	implicit none
	write(10,"(A)")'set term pngcairo'
	write(10,"(A)")'set output "order.png"'
	write(10,"(A)")'set pm3d map'
	write(10,"(A)")'set pm3d corners2color c1'
	write(10,"(A)")'set cbrange [:]'
	write(10,"(A)")'set pm3d interpolate 0,0'
	write(10,"(A)")'set size square'
	write(10,"(A)")'set palette rgbformulae 22,13,-31'
	write(10,"(A)")'splot "-" matrix index 0'
	write(10,"(A)")'#data'
	write(20,"(A)")'#data'
end

	
