module M_pmt
	use M_const
	implicit none
	integer, parameter :: Ns(2)=(/24,24/),Tx(2)=(/Ns(1),0/),Ty(2)=(/0,Ns(2)/),Ns2=Ns(1)*Ns(2),imp=Ns(2)/2*Ns(1)+Ns(1)/2
	integer :: neb(Ns2,4,3)
	real(8), parameter :: t(3)=(/1d0,-0.25d0,0d0/),nf=0.85d0,U=2.44d0,DJ=0d0,V=-1d0,Vimp=1d0,cvg=1d-5,bt=3d5
end module
module M_bdg
	use M_latt, only : &
		latt         =>  square         ,& 
		latt_bc      =>  square_bc      ,& 
		latt_one2two =>  square_one2two ,& 
		latt_two2one =>  square_two2one 
	use M_rd
	use M_utility
	use M_pmt
	use lapack95, only: heevd
	implicit none
contains
	subroutine bdg(dt,dd,ub,S,n)
		complex(8) :: dt(Ns2,2),dtp(Ns2,2),H(2*Ns2,2*Ns2)
		real(8) :: n(Ns2,2),np(Ns2,2),n1,wide,cp,E(2*Ns2),step,pn1,dd(Ns2,2),ddp(Ns2,2),ub(Ns2,2),ubp(Ns2,2),S(Ns2),Sp(Ns2)
		integer :: sg
		step=0.2d0
		do 
			pn1=0d0
			step=max(step,1d-2)
			do
				call EU(dt,dd,ub,S,n,cp,H,E)
				call phy(H,E,np,dtp,ddp,ubp,Sp)
				n1=sum(np)/Ns2
				if(abs(n1-nf)<cvg*5d0) then
					exit
				endif
				call find_cross(pn1,n1-nf,sg)
				if(sg/=0) then
					step=step*0.3d0
				endif
				cp=cp+step*sign(1d0,nf-n1)
				write(*,*)n1,cp
			enddo
			if(sum(abs(dt-dtp))/Ns2<cvg.and.sum(abs(S-Sp))/Ns2<cvg.and.sum(abs(dd-ddp))/Ns2<cvg) then
				exit
			endif
			dt=dtp
			dd=ddp
			n=np
			ub=ubp
			S=Sp
			write(*,*)"*",real(dt(1,1)),S(1),dd(1,1),ub(1,1),n(1,1),n(1,2)
		enddo 
	end subroutine
	subroutine EU(dt,dd,ub,S,ne,cp,H,E)
		complex(8) :: H(:,:),dt(:,:)
		real(8) :: ne(:,:),cp,E(:),dd(:,:),ub(:,:),S(:)
		integer :: i,j,k,n,info
		H=0d0
		do i=1,Ns2
			do n=1,2
				do k=1,size(t)
					j=neb(i,n,k)
					! hopping
					H(i,j)=H(i,j)-t(k)
					H(j,i)=H(j,i)-t(k)
					H(i+Ns2,j+Ns2)=H(i+Ns2,j+Ns2)+t(k)
					H(j+Ns2,i+Ns2)=H(j+Ns2,i+Ns2)+t(k)
					if(k>1) then
						cycle
					endif
					! uniform bond order
					H(i,j)=H(i,j)+ub(i,n)
					H(j,i)=H(j,i)+ub(i,n)
					H(i+Ns2,j+Ns2)=H(i+Ns2,j+Ns2)-ub(i,n)
					H(j+Ns2,i+Ns2)=H(j+Ns2,i+Ns2)-ub(i,n)
					! d-wave pairs
					H(i,j+Ns2)=H(i,j+Ns2)+dt(i,n)
					H(j,i+Ns2)=H(j,i+Ns2)+dt(i,n)
					H(i+Ns2,j)=H(i+Ns2,j)+conjg(dt(i,n))
					H(j+Ns2,i)=H(j+Ns2,i)+conjg(dt(i,n))
					! DDW
					H(i,j)=H(i,j)+dd(i,n)*img
					H(j,i)=H(j,i)-dd(i,n)*img
					H(i+Ns2,j+Ns2)=H(i+Ns2,j+Ns2)+dd(i,n)*img
					H(j+Ns2,i+Ns2)=H(j+Ns2,i+Ns2)-dd(i,n)*img
				enddo
			enddo
			! on site Hubbard U and chemical potential
			!H(i,i)=H(i,i)+S(i)-cp
			!H(i+Ns2,i+Ns2)=H(i+Ns2,i+Ns2)+S(i)+cp
			H(i,i)=H(i,i)+U*ne(i,2)-cp
			H(i+Ns2,i+Ns2)=H(i+Ns2,i+Ns2)-U*ne(i,1)+cp
			! impurity
			if(i==imp) then
				H(i,i)=H(i,i)+Vimp
				H(i+Ns2,i+Ns2)=H(i+Ns2,i+Ns2)+Vimp
			endif
		enddo
		call heevd(H,E,"V","U",info)
	end subroutine
	subroutine phy(H,E,ne,dt,dd,ub,S)
		real(8) :: ne(:,:),E(:),f(size(E)),th(size(E)),dd(:,:),ub(:,:),S(:)
		complex(8) :: H(:,:),dt(:,:)
		integer :: i,k,j,n
		ne=0d0
		dt=0d0
		dd=0d0
		ub=0d0
		S=0d0
		f=1d0/(exp(bt*E)+1d0)
		do i=1,Ns2
			do n=1,2*Ns2
				ne(i,1)=ne(i,1)+conjg(H(i,n))*H(i,n)*f(n)
				ne(i,2)=ne(i,2)+H(i+Ns2,n)*conjg(H(i+Ns2,n))*(1d0-f(n))
				do k=1,2
					j=neb(i,k,1)
					dt(i,k)=dt(i,k)+(H(i,n)*conjg(H(j+Ns2,n))+H(j,n)*conjg(H(i+Ns2,n)))*(1d0-f(n))
					!dd(i,k)=dd(i,k)+imag(conjg(H(i,n))*H(j,n)-conjg(H(j,n))*H(i,n)&
						!+conjg(H(i+Ns2,n))*H(j+Ns2,n)-conjg(H(j+Ns2,n))*H(i+Ns2,n))*f(n)
					!ub(i,k)=ub(i,k)+real(conjg(H(i,n))*H(j,n)+conjg(H(j,n))*H(i,n)&
						!-conjg(H(i+Ns2,n))*H(j+Ns2,n)-conjg(H(j+Ns2,n))*H(i+Ns2,n))*f(n)
				enddo
			enddo
		enddo
		dt=0.5d0*(DJ-V)*dt
		!dd=0.25d0*(0.5d0*DJ+V)*dd
		!ub=0.25d0*(0.5d0*DJ+V)*ub
		do i=1,Ns2
			do k=1,4
				S(i)=S(i)+0.25d0*DJ*(ne(neb(i,k,1),1)-ne(neb(i,k,1),2))
			enddo
			S(i)=S(i)-0.5d0*U*(ne(i,1)-ne(i,2))
		enddo
	end subroutine
end module
program main
	use M_bdg
	implicit none
	integer :: i,j,sg,ix(2)
	real(8) :: n(Ns2,2),dd(Ns2,2),ub(Ns2,2),S(Ns2)
	complex(8) :: dt(Ns2,2)
	logical :: f
	f=openfile(unit=10,file='../data/order.dat')
	call init_random_seed()
	call latt(Ns,Tx,Ty,neb)
	do i=1,Ns2
		call latt_one2two(i,Ns,ix)
		dt(i,:)=(/1d0,-1d0/)*0.05d0
		dd(i,:)=(/1d0,-1d0/)*(-1)**mod(sum(ix),2)*0d0
		ub(i,:)=(/1d0,1d0/)*0d0
		S(i)=(-1)**mod(sum(ix),2)*0.1d0
		n(i,1)=(nf-S(i))/2d0
		n(i,2)=(nf+S(i))/2d0
	enddo
	call bdg(dt,dd,ub,S,n)
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
	write(10,"(X)")
	do i=1,Ns2
		write(10,"(e15.6,$)")S(i)
	enddo
end program
