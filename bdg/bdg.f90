module global
	implicit none
	save
	real(8) :: bt=1e5,dv=1,nf0=0.85,du=2.44,t=1,tp=-0.25,eff=1
	integer, parameter :: dn=26
	integer :: nn(dn,dn),pbc(0:dn+1)
end module
program main
	use global
	implicit none
	character(25) :: it,ft,fmat(2)
	real(8) :: time(2)
	integer :: i,j,k,n,rc,im(5)=-1,info=0,np(4),error=0
	real(8) :: n1(dn*dn),n2(dn*dn),n1p(dn*dn),n2p(dn*dn),e(2*dn*dn),h(2*dn*dn,2*dn*dn),dt(dn*dn,4),&
		dt1(dn*dn,4),dtf(dn*dn),work(3*2*dn*dn),rdom(dn*dn),er(3)=1
	real(8) :: f,en,sa,sb,sp,sp0,nf,cvg=1e-4,wide
	logical :: flaga,flagb,flag1,flag2,flag3
	f(en)=1/(1+exp(bt*en))
	open(unit=10,file='../data/output.dat')
	open(unit=20,file='../data/error.dat')
	open(unit=30,file='../data/rundata.dat')
	open(unit=40,file='../data/initial.dat')
	open(unit=50,file='../data/gap.dat')
	write(fmat(1),*)"(",dn,"(e15.6))"
	write(fmat(2),*)"(",dn*dn,"(e15.6))"
200	format(24(e15.6))
300	format(625(e15.6))
500	format(100000e10.3)
	call fdate(it)
	rc=0
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
	! rdom=0.5d0	!uniform initial value
	do i=1,dn*dn
		dt1(i,:)=(/0.07d0,-0.07d0,0.07d0,-0.07d0/)+((/rdom(i)-0.5d0,-rdom(i)+0.5d0,rdom(i)-0.5d0,-rdom(i)+0.5d0/))/10d0
	enddo
	! n1p=rdom*nf0
	! n2p=nf0-n1p
	n1p=nf0/2+(rdom-0.5d0)/20d0
	n2p=nf0-n1p
	! call random_number(n1p)
	! call random_number(n2p)
	do i=1,dn*dn
		dtf(i)=(dt1(i,1)+dt1(i,3)-dt1(i,2)-dt1(i,4))/4.0
	enddo
	! write(40,fmat(1))((dtf(dn*i+j),j=1,dn),i=0,dn-1)
	! write(40,fmat(1))(((-1)**(i+j)*0.5*(n1p(dn*i+j)-n2p(dn*i+j)),j=1,dn),i=0,dn-1)
	flag1=.true.
	flag2=.true.
	flag3=.true.
	n1=0
	n2=0
	sp=0.14
	wide=0.2
	sp0=sp+wide
	! im(1)=121
	! im(2)=123
	im(1)=(dn/2-1)*dn+dn/2
	! im(2)=im(1)+1
	! im(3)=im(2)+dn
	! im(4)=im(1)+dn
	! im(5)=im(1)-1
	! do i=1,dn*dn
		! if(any(im==i)) then
			! n1(i)=1
		! endif
	! enddo
	! write(20,fmat(1))((n1(dn*i+j),j=1,dn),i=0,dn-1)
ex:	do while(any(er>cvg))
		rc=rc+1
		! if(er(1)>cvg) then
			dt=dt1
		! endif
		! if(any(er(2:3)>cvg)) then
			n1=n1p
			n2=n2p
		! endif
		nf=0
		sa=sp
		sb=sp
		flaga=.true.
		flagb=.true.
		! do while(abs(nf-nf0)>1e-3)
		do while(abs(nf-nf0)>cvg)
			call cpu_time(time(1))
			sp=0.5*(sa+sb)
			h=0d0
			do i=1,dn*dn
				! nearest neighbor pairs
				call nnp(i,.true.,np)
				do j=1,4
					h(2*i-1:2*i,2*np(j)-1:2*np(j))=h(2*i-1:2*i,2*np(j)-1:2*np(j))+reshape((/-t,dt(i,j),dt(i,j),&
						t/),(/2,2/))
				enddo
				! next-nearest neighbor pairs
				call nnp(i,.false.,np)
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
			write(*,*)nf
			if(abs(nf-nf0)<=cvg) then
				exit
			endif
			if(nf<nf0) then
				flaga=.false.
				sa=sp
				if(flaga.or.flagb) then
					sb=sp+wide
					sp=sb
				endif
			else
				flagb=.false.
				sb=sp
				if(flaga.or.flagb) then
					sa=sp-wide
					sp=sa
				endif
			endif
			call cpu_time(time(2))
			write(10,*)"cost cup time:",time(2)-time(1),"s"
		enddo
		wide=max(abs(sp0-sp),10*cvg)
		sp0=sp
		n1p=0
		n2p=0
		do i=1,dn*dn
			do n=1,2*dn*dn
				n1p(i)=n1p(i)+h(2*i-1,n)**2*f(e(n))
				n2p(i)=n2p(i)+h(2*i,n)**2*(1-f(e(n)))
			enddo
		enddo
		! rewind(30)
		! write(30,fmat(1))((n1p(dn*i+j)+n2p(dn*i+j),j=1,dn),i=0,dn-1)
		dt1=0
		do i=1,dn*dn
			call nnp(i,.true.,np)
			do j=1,4
				do n=1,2*dn*dn
					dt1(i,j)=dt1(i,j)+0.25*dv*(h(2*i-1,n)*h(2*np(j),n)+h(2*np(j)-1,n)*h(2*i,n))*tanh(0.5*bt*e(n))
				enddo
			enddo
			write(30,"(4e10.3,5i4)")dt1(i,:),i,np
		enddo
		write(*,*)dt1(1,1)
		er(1)=sum(abs(dt1-dt))/(dn*dn)
		er(2)=sum(abs(n1p-n1))/(dn*dn)
		er(3)=sum(abs(n2p-n2))/(dn*dn)
		write(*,"(3e15.6,f15.12)")er,sa
		! if(mod(rc,20)==0.and.rc>=100) then
			! do i=1,dn*dn
				! dtf(i)=(dt1(i,1)+dt1(i,3)-dt1(i,2)-dt1(i,4))/4.0
			! enddo
			! write(50,*)"n=",rc
			! write(50,fmat(1))((dtf(dn*i+j),j=1,dn),i=0,dn-1)
		! endif
		! if(flag1.and.all(er<cvg*100)) then
			! do i=1,dn*dn
				! dtf(i)=(dt1(i,1)+dt1(i,3)-dt1(i,2)-dt1(i,4))/4.0
			! enddo
			! write(10,fmat(1))((dtf(dn*i+j),j=1,dn),i=0,dn-1)
			! write(10,fmat(1))((n1p(dn*i+j)+n2p(dn*i+j),j=1,dn),i=0,dn-1)
			! write(10,fmat(1))(((-1)**(i+j)*0.5*(n1p(dn*i+j)-n2p(dn*i+j)),j=1,dn),i=0,dn-1)
			! write(10,*)"error=",cvg*100
			! flag1=.false.
		! endif		
		! if(flag2.and.all(er<cvg*50)) then
			! do i=1,dn*dn
				! dtf(i)=(dt1(i,1)+dt1(i,3)-dt1(i,2)-dt1(i,4))/4.0
			! enddo
			! write(10,fmat(1))((dtf(dn*i+j),j=1,dn),i=0,dn-1)
			! write(10,fmat(1))((n1p(dn*i+j)+n2p(dn*i+j),j=1,dn),i=0,dn-1)
			! write(10,fmat(1))(((-1)**(i+j)*0.5*(n1p(dn*i+j)-n2p(dn*i+j)),j=1,dn),i=0,dn-1)
			! write(10,*)"error=",cvg*50
			! flag2=.false.
		! endif
		! if(flag3.and.all(er<cvg*10)) then
			! do i=1,dn*dn
				! dtf(i)=(dt1(i,1)+dt1(i,3)-dt1(i,2)-dt1(i,4))/4.0
			! enddo
			! write(10,fmat(1))((dtf(dn*i+j),j=1,dn),i=0,dn-1)
			! write(10,fmat(1))((n1p(dn*i+j)+n2p(dn*i+j),j=1,dn),i=0,dn-1)
			! write(10,fmat(1))(((-1)**(i+j)*0.5*(n1p(dn*i+j)-n2p(dn*i+j)),j=1,dn),i=0,dn-1)
			! write(10,*)"error=",cvg*10
			! flag3=.false.
		! endif
	enddo ex
	do i=1,dn*dn
		dtf(i)=(dt(i,1)+dt(i,3)-dt(i,2)-dt(i,4))/4.0
	enddo
	call fdate(ft)
	write(10,fmat(1))((dtf(dn*i+j),j=1,dn),i=0,dn-1)
	write(10,fmat(1))((n1(dn*i+j)+n2(dn*i+j),j=1,dn),i=0,dn-1)
	write(10,fmat(1))(((-1)**(i+j)*0.5*(n1(dn*i+j)-n2(dn*i+j)),j=1,dn),i=0,dn-1)
	write(10,*)"the beta is",bt,"v=",dv,"u=",du,"t'=",tp,"heff=",eff,"nf=",nf0,"n=",dn
	write(10,*)"error=",error,"info=",info,"rc=",rc
	write(10,*)"the tatal runtime is form ",it," to ",ft,"cost cup time:",time(2)-time(1),"s"
	if(info/=0) then
		error=mod(info,2*dn*dn+1)
		do j=1,error+1
			write(20,500)(h(i,j),i=1,error+1)
		enddo
	endif
	close(10)
	close(20)
	close(30)
	close(40)
	close(50)
end program main

subroutine nnp(ii,near,jj)
	use global
	implicit none
	integer :: ii,i,j,jj(4)
	logical :: near
	j=pbc(mod(ii,dn))
	i=(ii-j)/dn+1
	if(near) then
		jj(1)=(i-1)*dn+pbc(j+1)
		jj(2)=(pbc(i+1)-1)*dn+j
		jj(3)=(i-1)*dn+pbc(j-1)
		jj(4)=(pbc(i-1)-1)*dn+j
	else
		jj(1)=dn*(pbc(i+1)-1)+pbc(j+1)
		jj(2)=dn*(pbc(i+1)-1)+pbc(j-1)
		jj(3)=dn*(pbc(i-1)-1)+pbc(j+1)
		jj(4)=dn*(pbc(i-1)-1)+pbc(j-1)
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

	
