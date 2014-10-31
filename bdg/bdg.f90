module M_pmt
	use M_const
	implicit none
	integer, parameter :: Ns(2)=(/24,24/),Ns2=Ns(1)*Ns(2),imp=Ns2/2+Ns(1)/2
	integer :: neb(Ns2,4,3)
	real(8), parameter :: t(3)=(/1d0,-0.25d0,0d0/),nf=0.85d0,U=2.44d0,V=1d0,Vimp=1d0,cvg=1e-4,bt=1e5
	complex(8), parameter :: ipi=2d0*pi*img
end module
module M_bdg
	use M_latt
	use M_rd
	use M_pmt
	use f95_precision, only : wp => dp
	use lapack95, only: heevd
	implicit none
contains
	subroutine inital(dt,n)
		complex(8) :: dt(Ns2,4)
		real(8) :: n(Ns2,2)
		integer :: i
		call init_random_seed()
		!call random_number(rdom)
		do i=1,Ns2
			dt(i,:)=(/0.05d0,-0.05d0,0.05d0,-0.05d0/)
		enddo
		n(:,1)=nf/2d0
		n(:,2)=nf-n(:,1)
		call square(Ns,neb)
	end subroutine
	subroutine bdg(dt,n)
		complex(8) :: dt(Ns2,4),dtp(Ns2,4),H(2*Ns2,2*Ns2)
		real(8) :: n(Ns2,2),np(Ns2,2),n1,wide,sp0,sp,sa,sb,E(2*Ns2)
		logical :: flaga,flagb
		wide=0.5d0
		sp0=sp+wide
		do 
			sa=sp
			sb=sp
			flaga=.true.
			flagb=.true.
			do
				sp=0.5d0*(sa+sb)
				call EU(dt,n,sp,H,E)
				call order(H,E,np,dtp)
				n1=1.0/Ns2*sum(np)
				if(abs(n1-nf)<cvg) then
					exit
				endif
				if(n1<nf) then
					flaga=.false.
					sa=sp
					if(flagb) then
						sb=sp+wide
					endif
				else
					flagb=.false.
					sb=sp
					if(flaga) then
						sa=sp-wide
					endif
				endif
				write(*,*)n1,sp
			enddo
			wide=max(abs(sp0-sp),100*cvg)
			sp0=sp
			if(sum(abs(dt-dtp))/Ns2<cvg) then
				exit
			endif
			dt=dtp
			n=np
			write(*,*)"*",dt(1,1)
		enddo 
	end subroutine
	subroutine EU(dt,n,sp,H,E)
		complex(8) :: H(2*Ns2,2*Ns2),dt(Ns2,4)
		real(8) :: E(2*Ns2),n(Ns2,2),sp
		integer :: i,j,info
		H=0d0
		do i=1,Ns2
			do j=1,4
				! nearest neighbor pairs
				H(i,neb(i,j,1))=-t(1)
				H(i,neb(i,j,1)+Ns2)=dt(i,j)
				H(i+Ns2,neb(i,j,1)+Ns2)=t(1)
				H(i+Ns2,neb(i,j,1))=dconjg(dt(i,j))
				! next-nearest neighbor pairs
				H(i,neb(i,j,2))=-t(2)
				H(i+Ns2,neb(i,j,2)+Ns2)=t(2)
			enddo
			! on site
			H(i,i)=U*n(i,2)-sp
			H(i+Ns2,i+Ns2)=-U*n(i,1)+sp
		enddo
		! impure
		H(imp,imp)=H(imp,imp)+Vimp
		H(imp+Ns2,imp+Ns2)=H(imp+Ns2,imp+Ns2)+Vimp
		call heevd(H,E,"V","U",info)
	end subroutine
	subroutine order(H,E,np,dtp)
		real(8) :: np(Ns2,2),f(2*Ns2),th(2*Ns2),E(2*Ns2)
		complex(8) :: dtp(Ns2,4),H(2*Ns2,2*Ns2)
		integer :: i,j,n
		np=0d0
		dtp=0d0
		f=1d0/(exp(bt*E)+1d0)
		th=tanh(0.5d0*bt*E)
		do i=1,Ns2
			do n=1,2*Ns2
				np(i,1)=np(i,1)+H(i,n)*dconjg(H(i,n))*f(n)
				np(i,2)=np(i,2)+H(i+Ns2,n)*dconjg(H(i+Ns2,n))*(1d0-f(n))
				do j=1,4
					dtp(i,j)=dtp(i,j)+0.25d0*V*(H(i,n)*dconjg(H(neb(i,j,1)+Ns2,n))+H(neb(i,j,1),n)*dconjg(H(i+Ns2,n)))*th(n)
				enddo
			enddo
		enddo
	end subroutine
end module
program main
	use M_bdg
	implicit none
	character(25) :: fmat(2)
	integer :: i,j,sg
	real(8) :: n(Ns2,2)
	complex(8) :: dt(Ns2,4)
	open(unit=10,file='../data/order.dat')
	call gnuplot
	write(fmat(1),*)Ns(1)
	200	format(24(e15.6))
	call inital(dt,n)
	call bdg(dt,n)
	!export data
	do i=1,Ns2
		write(10,"(e15.6,$)")abs(dt(i,1)+dt(i,3)-dt(i,2)-dt(i,4))/4.0d0
		if(mod(i,Ns(1))==0) then
			write(10,"(X)")
		endif
	enddo
	write(10,"(X)")
	do i=1,Ns2
		write(10,"(e15.6,$)")n(i,1)+n(i,2)
		if(mod(i,Ns(1))==0) then
			write(10,"(X)")
		endif
	enddo
	write(10,"(X)")
	sg=-1
	do i=1,Ns2
		sg=-sg
		write(10,"(e15.6,$)")sg*0.5d0*(n(i,1)-n(i,2))
		if(mod(i,Ns(1))==0) then
			sg=-sg
			write(10,"(X)")
		endif
	enddo
	close(10)
end program
subroutine gnuplot()
	implicit none
	write(10,"(A)")"set term pngcairo"
	write(10,"(A)")"set output 'order.png'"
	write(10,"(A)")"set pm3d map"
	write(10,"(A)")"set pm3d corners2color c1"
	write(10,"(A)")"set cbrange [:]"
	write(10,"(A)")"set pm3d interpolate 0,0"
	write(10,"(A)")"set size square"
	write(10,"(A)")"set palette rgbformulae 22,13,-31"
	write(10,"(A)")"splot '-' matrix index 0"
	write(10,"(A)")"#data"
end subroutine
