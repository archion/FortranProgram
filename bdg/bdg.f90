module M_pmt
	use M_const
	implicit none
	integer, parameter :: Ns(2)=(/24,24/),Tx(2)=(/Ns(1),0/),Ty(2)=(/0,Ns(2)/),Ns2=Ns(1)*Ns(2),imp=Ns(2)/2*Ns(1)+Ns(1)/2
	integer :: neb(Ns2,4,3)
	real(8), parameter :: t(3)=(/1d0,-0.25d0,0d0/),nf=0.85d0,U=2.44d0,V=1d0,Vimp=1d0,cvg=1d-5,bt=3d5
end module
module M_bdg
	use M_latt, only : latt => square
	use M_rd
	use M_utility
	use M_pmt
	use lapack95, only: heevd
	implicit none
contains
	subroutine bdg(dt,n)
		complex(8) :: dt(Ns2,2),dtp(Ns2,2),H(2*Ns2,2*Ns2)
		real(8) :: n(Ns2,2),np(Ns2,2),n1,wide,sp0,sp,sa,sb,E(2*Ns2),step,pn1
		integer :: sg
		step=0.2d0
		do 
			pn1=0d0
			step=max(step,1d-4)
			do
				call EU(dt,n,sp,H,E)
				call phy(H,E,np,dtp)
				n1=sum(np)/Ns2
				if(abs(n1-nf)<cvg*5d0) then
					exit
				endif
				call find_cross(pn1,n1-nf,sg)
				if(sg/=0) then
					step=step*0.3d0
				endif
				sp=sp+step*sign(1d0,nf-n1)
				write(*,*)n1,sp
			enddo
			if(sum(abs(dt-dtp))/Ns2<cvg) then
				exit
			endif
			dt=dtp
			n=np
			write(*,*)"*",dt(1,1)
		enddo 
	end subroutine
	subroutine EU(dt,ne,sp,H,E)
		complex(8) :: H(:,:),dt(:,:)
		real(8) :: ne(:,:),sp,E(:)
		integer :: i,j,k,n,info
		H=0d0
		do i=1,Ns2
			do n=1,2
				! hopping
				do k=1,size(t)
					H(i,neb(i,n,k))=-t(k)
					H(neb(i,n,k),i)=-t(k)
					H(i+Ns2,neb(i,n,k)+Ns2)=t(k)
					H(neb(i,n,k)+Ns2,i+Ns2)=t(k)
				enddo
				! d-wave pairs
				H(i,neb(i,n,1)+Ns2)=dt(i,n)
				H(neb(i,n,1),i+Ns2)=dt(i,n)
				H(i+Ns2,neb(i,n,1))=conjg(dt(i,n))
				H(neb(i,n,1)+Ns2,i)=conjg(dt(i,n))
			enddo
			! on site Hubbard U and chemical potential
			H(i,i)=U*ne(i,2)-sp
			H(i+Ns2,i+Ns2)=-U*ne(i,1)+sp
		enddo
		! impurity
		H(imp,imp)=H(imp,imp)+Vimp
		H(imp+Ns2,imp+Ns2)=H(imp+Ns2,imp+Ns2)+Vimp

		call heevd(H,E,"V","U",info)
	end subroutine
	subroutine phy(H,E,ne,dt)
		real(8) :: ne(:,:),E(:),f(size(E)),th(size(E))
		complex(8) :: dt(:,:),H(:,:)
		integer :: i,j,n
		ne=0d0
		dt=0d0
		f=1d0/(exp(bt*E)+1d0)
		do i=1,Ns2
			do n=1,2*Ns2
				ne(i,1)=ne(i,1)+conjg(H(i,n))*H(i,n)*f(n)
				ne(i,2)=ne(i,2)+conjg(H(i+Ns2,n))*H(i+Ns2,n)*f(n)
				do j=1,2
					dt(i,j)=dt(i,j)+(H(i,n)*conjg(H(neb(i,j,1)+Ns2,n))+H(neb(i,j,1),n)*conjg(H(i+Ns2,n)))*f(n)
				enddo
			enddo
			ne(i,2)=1d0-ne(i,2)
		enddo
		dt=-0.5d0*V*dt
	end subroutine
end module
program main
	use M_bdg
	implicit none
	integer :: i,j,sg
	real(8) :: n(Ns2,2)
	complex(8) :: dt(Ns2,2)
	logical :: f
	f=openfile(unit=10,file='../data/order.dat')
	call init_random_seed()
	do i=1,Ns2
		dt(i,:)=(/0.05d0,-0.05d0/)
	enddo
	n(:,1)=nf/2d0
	n(:,2)=nf-n(:,1)
	call latt(Ns,Tx,Ty,neb)
	call bdg(dt,n)
	!export data
	do i=1,Ns2
		write(10,"(e15.6,$)")abs(dt(neb(i,1,1),1)-dt(neb(i,2,1),2)+dt(neb(i,3,1),1)-dt(neb(i,4,1),2))/4.0d0
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
end program
