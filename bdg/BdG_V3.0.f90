module global
	implicit none
	real(8),save :: bt=1e5,dv=1,nf0=0.85,du=2.44,t=1,tp=-0.25,eff=0.3
	integer, parameter :: dn=24
	integer,save :: nn(dn,dn),pbc(0:dn+1)
end module
program main
	use global
	implicit none
	character(24) :: it,ft
	integer :: i,j,k,n,im(5)=-1,info=0,np(4),error=0
	real(8) :: n1(dn*dn),n2(dn*dn),n1p(dn*dn),n2p(dn*dn),e(2*dn*dn),h(2*dn*dn,2*dn*dn),dt(dn*dn,dn*dn),&
		dt1(dn*dn,dn*dn),dtf(dn*dn),work(3*2*dn*dn),rdom(dn*dn)
	real(8) :: f,en,sa,sb,sp,nf,er=1
	f(en)=1/(1+exp(bt*en))
	open(unit=10,file='../data/output.dat')
	open(unit=20,file='../data/error.dat')
200	format(24(e15.6))
300	format(625(e15.6))
	call fdate(it)
	k=1
	pbc=(/dn,(i,i=1,dn),1/)
	do i=1,dn
		do j=1,dn
			nn(i,j)=k
			k=k+1
		enddo
	enddo
	dt=0d0
	dt1=0d0
	call init_random_seed()
	call random_number(rdom)
	rdom=0.5d0	!uniform initial value
	do i=1,dn
		call nnp(i,0,np)
		dt1(i,np)=(/-0.0795d0,0.0795d0,-0.0795d0,0.0795d0/)+(rdom(4*i-3:4*i)-0.5d0)/10d0
	enddo
	n1p=nf0/2+(rdom-0.5d0)/20d0
	n2p=nf0-n1p
	n1=0
	n2=0
	! im(1)=121
	! im(2)=123
	! im(1)=(dn/2-1)*dn+dn/2
	! im(2)=im(1)+1
	! im(3)=im(2)+dn
	! im(4)=im(1)+dn
	! im(5)=im(1)-1
	! do i=1,dn*dn
		! if(any(im==i)) then
			! n1(i)=1
		! endif
	! enddo
	! write(20,200)((n1(dn*i+j),j=1,dn),i=0,dn-1)
ex:	do while(er>1e-4)
		dt=dt1
		n1=n1p
		n2=n2p
		nf=0
		sa=0
		sb=0.5
		do while(abs(nf-nf0)>1e-5)
			sp=0.5*(sa+sb)
			h=0d0
			do i=1,dn*dn
				! nearest neighbor pairs
				call nnp(i,0,np)
				do j=1,4
					h(2*i-1:2*i,2*np(j)-1:2*np(j))=h(2*i-1:2*i,2*np(j)-1:2*np(j))+reshape((/-t,dt(i,np(j)),dt(i,np(j)),&
						t/),(/2,2/))
				enddo
				! next-nearest neighbor pairs
				call nnp(i,1,np)
				do j=1,4
					h(2*i-1:2*i,2*np(j)-1:2*np(j))=h(2*i-1:2*i,2*np(j)-1:2*np(j))+reshape((/-tp,0.0d0,0.0d0,tp/),(/2,2/))
				enddo
				! on site
				h(2*i-1:2*i,2*i-1:2*i)=h(2*i-1:2*i,2*i-1:2*i)+reshape((/du*n2(i)-sp,0.0d0,0.0d0,-du*n1(i)+sp/),&
					(/2,2/))
				if(any(im==i)) then
					h(2*i-1:2*i,2*i-1:2*i)=h(2*i-1:2*i,2*i-1:2*i)+reshape((/(-1)**minloc(abs(im-i))*eff,0.0d0,&
						0.0d0,(-1)**minloc(abs(im-i))*eff/),(/2,2/))
				endif
			enddo
			call dsyev( "v", "u", 2*dn*dn, h, 2*dn*dn, e, work, 3*2*dn*dn, info )
			if(info/=0) then
				exit ex
			endif
			nf=0
			do i=1,dn*dn
				do n=1,2*dn*dn
					nf=nf+(h(2*i-1,n)**2*f(e(n))+h(2*i,n)**2*(1-f(e(n))))
				enddo
			enddo
			nf=1.0/(dn*dn)*nf
			if (nf<nf0) then
				sa=sp
			else
				sb=sp
			endif
		enddo
		n1p=0
		n2p=0
		do i=1,dn*dn
			do n=1,2*dn*dn
				n1p(i)=n1p(i)+h(2*i-1,n)**2*f(e(n))
				n2p(i)=n2p(i)+h(2*i,n)**2*(1-f(e(n)))
			enddo
		enddo
		er=sum(abs(n1p-n1)+abs(n2p-n2))/(dn*dn)
		dt1=0
		do i=1,dn*dn
			call nnp(i,0,np)
			do j=1,4
				do n=1,2*dn*dn
					dt1(i,np(j))=dt1(i,np(j))+(h(2*i-1,n)*h(2*np(j),n)+h(2*np(j)-1,n)*h(2*i,n))*tanh(0.5*bt*e(n))
				enddo
			enddo
		enddo
		dt1=0.25*dv*dt1
		write(*,"(e15.6,f15.12)")er,sa
		rewind(10)
		write(10,200)((n1(dn*i+j)+n2(dn*i+j),j=1,dn),i=0,dn-1)
		rewind(20)
		do j=1,dn*dn
			write(20,500)(dt(i,j),i=1,dn*dn)
		enddo
	enddo ex
	rewind(10)
	do i=1,dn*dn
		call nnp(i,0,np)
		dtf(i)=(dt(i,np(1))+dt(i,np(3))-dt(i,np(2))-dt(i,np(4)))/4.0
	enddo
	call fdate(ft)
	write(10,200)((dtf(dn*i+j),j=1,dn),i=0,dn-1)
	write(10,200)((n1(dn*i+j)+n2(dn*i+j),j=1,dn),i=0,dn-1)
	write(10,200)(((-1)**(i+j)*0.5*(n1(dn*i+j)-n2(dn*i+j)),j=1,dn),i=0,dn-1)
	write(10,*)"the beta is",bt,"v=",dv,"u=",du,"t'=",tp,"heff=",eff,"nf=",nf0,"n=",dn
	write(10,*)"error=",error,"info=",info
	write(10,*)"the tatal runtime is form ",it," to ",ft
	if(info/=0) then
		error=mod(info,2*dn*dn+1)
		do j=1,error+1
			write(20,500)(h(i,j),i=1,error+1)
		enddo
	endif
	write(20,500)(dt(i,:),i=1,dn*dn)
500	format(1000000(e10.3))
	close(10)
	close(20)
end program main

subroutine nnp(ii,p,jj)
	use global
	implicit none
	integer :: ii,i,j,jj(4),p
lp:	do i=1,dn
		do j=1,dn
			if(ii==nn(i,j)) then
				exit lp
			endif
		enddo
	enddo lp
	if(p==0) then
		jj(1)=nn(i,pbc(j+1))
		jj(2)=nn(pbc(i+1),j)
		jj(3)=nn(i,pbc(j-1))
		jj(4)=nn(pbc(i-1),j)
	endif
	if(p==1) then
		jj(1)=nn(pbc(i+1),pbc(j+1))
		jj(2)=nn(pbc(i+1),pbc(j-1))
		jj(3)=nn(pbc(i-1),pbc(j+1))
		jj(4)=nn(pbc(i-1),pbc(j-1))
	endif
end subroutine nnp

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

	
