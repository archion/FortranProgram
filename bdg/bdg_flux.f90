module global
	implicit none
	save
	integer, parameter :: dnx=17,dny=2*dnx
	real(8), parameter :: pi=3.1415926535898
	real(8) :: bt=1e5,dv=1.0d0,nf0=0.84,du=2.0d0,t=1,tp=-0.25,eff=1,fy0=2d0*pi/(dnx*dny)
	complex(8), parameter :: img=(0d0,1d0)
	integer :: nn(dnx,dny),pbcx(0:dnx+1),pbcy(0:dny+1)
end module
module lapack
	implicit none
	contains
	subroutine czheev(jobz,uplo,a,w,info)
		implicit none
		character :: jobz,uplo
		integer :: n,lda,lwork,info
		real(8) :: w(:)
		complex(8) :: a(:,:)
		real(8), allocatable :: rwork(:)
		complex(8), allocatable :: work(:)
		n=size(w)
		lda=n
		lwork=3*n
		allocate(work(lwork),rwork(3*n-2))
		call zheev (jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
	end subroutine
end module
program main
	use global
	use lapack
	implicit none
	character(25) :: it,ft,fmat(2)
	real(8) :: time(2)
	integer :: i,j,k,n,rc,im(5)=-1,info=0,np(4),error=0
	real(8) :: n1(dnx*dny),n2(dnx*dny),n1p(dnx*dny),n2p(dnx*dny),e(2*dnx*dny),dtf(dnx*dny),&
		rdom(dnx*dny),er(3)=1
	complex(8) :: h(2*dnx*dny,2*dnx*dny),dt(dnx*dny,4),dt1(dnx*dny,4),fy
	real(8) :: f,en,sa,sb,sp,sp0,nf,cvg=1e-4,wide
	logical :: flaga,flagb,flag1,flag2,flag3
	f(en)=1/(1+exp(bt*en))
	open(unit=10,file='../data/dsc.dat')
	open(unit=20,file='../data/sdw.dat')
	open(unit=30,file='../data/cdw.dat')
	open(unit=40,file='../data/initial.dat')
	open(unit=50,file='../data/gap.dat')
	write(fmat(1),*)"(2i3e15.6)"
	write(fmat(2),*)"(",dnx*dny,"(e15.6))"
200	format(24(e15.6))
300	format(625(e15.6))
500	format(100000e10.3)
	call fdate(it)
	rc=0
	k=1
	pbcx=(/dnx,(i,i=1,dnx),1/)
	pbcy=(/dny,(i,i=1,dny),1/)
	do i=1,dny
		do j=1,dnx
			nn(i,j)=k
			k=k+1
		enddo
	enddo
	dt=0d0
	dt1=0d0
	call init_random_seed()
	call random_number(rdom)
	! rdom=0.5d0	!uniform initial value
	do i=1,dnx*dny
		dt1(i,:)=(/0.07d0,-0.07d0,0.07d0,-0.07d0/)+((/rdom(i)-0.5d0,-rdom(i)+0.5d0,rdom(i)-0.5d0,-rdom(i)+0.5d0/))/10d0
	enddo
	! n1p=rdom*nf0
	! n2p=nf0-n1p
	n1p=nf0/2+(rdom-0.5d0)/20d0
	n2p=nf0-n1p
	! call random_number(n1p)
	! call random_number(n2p)
	do i=1,dnx*dny
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
	im(1)=(dny/2-1)*dnx+dnx/2
	! im(2)=im(1)+1
	! im(3)=im(2)+dn
	! im(4)=im(1)+dn
	! im(5)=im(1)-1
	! do i=1,dnx*dny
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
			h=(0d0,0d0)
			do i=1,dnx*dny
				! nearest neighbor pairs
				call nnp(i,.true.,np)
				do j=1,4
					fy=t*exp(img*mod(j,2)*(mod(i,dnx)+1)*fy0*(-1)**(j/3))
					h(2*i-1:2*i,2*np(j)-1:2*np(j))=h(2*i-1:2*i,2*np(j)-1:2*np(j))+reshape((/-fy,dconjg(dt(i,j)),dt(i,j),&
						fy/),(/2,2/))
				enddo
				! next-nearest neighbor pairs
				!call nnp(i,.false.,np)
				!do j=1,4
					!h(2*i-1:2*i,2*np(j)-1:2*np(j))=h(2*i-1:2*i,2*np(j)-1:2*np(j))+reshape((/-tp,0.0d0,0.0d0,tp/),(/2,2/))
				!enddo
				! on site
				h(2*i-1:2*i,2*i-1:2*i)=h(2*i-1:2*i,2*i-1:2*i)+reshape((/du*n2(i)-sp,0.0d0,0.0d0,-du*n1(i)+sp/),&
					(/2,2/))
			   ! if(any(im==i)) then
					!h(2*i-1:2*i,2*i-1:2*i)=h(2*i-1:2*i,2*i-1:2*i)+reshape((/(-1)**minloc(abs(im-i))*eff,0.0d0,&
						!0.0d0,(-1)**minloc(abs(im-i))*eff/),(/2,2/))
				!endif
			enddo
			call czheev( "v", "u", h, e, info )
			if(info/=0) then
				exit ex
			endif
			nf=0
			do i=1,dnx*dny
				do n=1,2*dnx*dny
					nf=nf+(h(2*i-1,n)*dconjg(h(2*i-1,n))*f(e(n))+h(2*i,n)*dconjg(h(2*i,n))*(1-f(e(n))))
				enddo
			enddo
			nf=1.0/(dnx*dny)*nf
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
			!write(10,*)"cost cup time:",time(2)-time(1),"s"
		enddo
		wide=max(abs(sp0-sp),10*cvg)
		sp0=sp
		n1p=0
		n2p=0
		do i=1,dnx*dny
			do n=1,2*dnx*dny
				n1p(i)=n1p(i)+h(2*i-1,n)*dconjg(h(2*i-1,n))*f(e(n))
				n2p(i)=n2p(i)+h(2*i,n)*dconjg(h(2*i,n))*(1-f(e(n)))
			enddo
		enddo
		! rewind(30)
		! write(30,fmat(1))((n1p(dn*i+j)+n2p(dn*i+j),j=1,dn),i=0,dn-1)
		dt1=0
		do i=1,dnx*dny
			call nnp(i,.true.,np)
			do j=1,4
				do n=1,2*dnx*dny
					dt1(i,j)=dt1(i,j)+0.25*dv*(h(2*i-1,n)*dconjg(h(2*np(j),n))+h(2*np(j)-1,n)*dconjg(h(2*i,n)))*tanh(0.5*bt*e(n))
				enddo
			enddo
			!write(30,"(4e10.3,5i4)")dt1(i,:),i,np
		enddo
		write(*,*)dt1(1,1)
		er(1)=sum(abs(dt1-dt))/(dnx*dny)
		er(2)=sum(abs(n1p-n1))/(dnx*dny)
		er(3)=sum(abs(n2p-n2))/(dnx*dny)
		write(*,"(3e15.6,f15.12)")er,sa
		enddo ex
	do i=1,dnx*dny
		dtf(i)=real(dt(i,1)+dt(i,3)-dt(i,2)-dt(i,4))/4.0
	enddo
	call fdate(ft)
	do i=0,dny-1
		do j=1,dnx
			write(10,"(2i3,e16.3)")i+1,j,dtf(dnx*i+j)
			write(20,"(2i3,e16.3)")i+1,j,(-1)**(i+j)*0.5*(n1(dnx*i+j)-n2(dnx*i+j))
			write(30,"(2i3,e16.3)")i+1,j,n1(dnx*i+j)+n2(dnx*i+j)
		enddo
		write(10,"(1x)")
		write(20,"(1x)")
		write(30,"(1x)")
	enddo
	!write(10,fmat(1))((dtf(dnx*i+j),j=1,dnx),i=0,dny-1)
	!write(30,fmat(1))((n1(dnx*i+j)+n2(dnx*i+j),j=1,dnx),i=0,dny-1)
	!write(20,fmat(1))(((-1)**(i+j)*0.5*(n1(dnx*i+j)-n2(dnx*i+j)),j=1,dnx),i=0,dny-1)
	!write(10,*)"the beta is",bt,"v=",dv,"u=",du,"t'=",tp,"heff=",eff,"nf=",nf0,"n=",dn
	!write(10,*)"error=",error,"info=",info,"rc=",rc
	!write(10,*)"the tatal runtime is form ",it," to ",ft,"cost cup time:",time(2)-time(1),"s"
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
	j=pbcx(mod(ii,dnx))
	i=(ii-j)/dnx+1
	if(near) then
		jj(1)=(i-1)*dnx+pbcx(j+1)
		jj(2)=(pbcy(i+1)-1)*dnx+j
		jj(3)=(i-1)*dnx+pbcx(j-1)
		jj(4)=(pbcy(i-1)-1)*dnx+j
	else
		jj(1)=dnx*(pbcy(i+1)-1)+pbcx(j+1)
		jj(2)=dnx*(pbcy(i+1)-1)+pbcx(j-1)
		jj(3)=dnx*(pbcy(i-1)-1)+pbcx(j+1)
		jj(4)=dnx*(pbcy(i-1)-1)+pbcx(j-1)
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

	
