module global
	implicit none
	real(8),save :: bt=1e5,dv=1.0,nf0=0.85,du=2.44,t=1.0,tp=-0.25,eff=1.0
	integer, parameter :: dn=24
	integer,save,allocatable :: nn(:,:),pbc(:)
end module
program main
	use global
	implicit none
	integer :: i,j,k=1,cc,n,im(2),info=0,np(4),error=0
	integer, allocatable :: iwork(:)
	real(8), allocatable :: n1(:),n2(:),n1p(:),n2p(:),e(:),f1(:,:),f2(:,:),rwork(:),dtr(:,:)
	complex(8), allocatable :: h(:,:),dt(:,:),dt1(:,:),dtf(:),u(:,:),v(:,:),ap(:,:),work(:)
	real(8) :: f,en,sa,sb,sp,s0,nf,cache,time1,time2
	f(en)=1/(1+exp(bt*en))
	allocate(h(2*dn*dn,2*dn*dn),n1(dn*dn),n2(dn*dn),n1p(dn*dn),n2p(dn*dn),dt(dn*dn,dn*dn),dt1(dn*dn,dn*dn),dtr(dn*dn,dn*dn),&
		dtf(dn*dn),u(dn*dn,2*dn*dn),v(dn*dn,2*dn*dn),e(2*dn*dn),work(4*dn*dn+(2*dn*dn)**2),rwork(1+10*dn*dn+2*(2*dn*dn)**2),&
		iwork(3+10*dn*dn),ap(2*dn*dn,2*dn*dn),nn(dn,dn),pbc(0:dn+1),f1(dn*dn,dn*dn),f2(dn*dn,dn*dn))
	open(unit=10,file='../data/output.dat')
	open(unit=20,file='../data/error.dat')
200	format(24("(",e15.6,",",e15.6,")  "))
300	format(576("(",e15.6,",",e15.6,")  "))
400	format(24e15.6)
	pbc=(/dn,(i,i=1,dn),1/)
	do i=1,dn
		do j=1,dn
			nn(i,j)=k
			k=k+1
		enddo
	enddo
	dt=0d0
	dt1=0.001d0
	call init_random_seed()
	call random_number(dtr)
	dt1=0.09
	! call random_number(dtr)
	! dt1=dt1+cmplx(0,dtr)
	! dt1=-0.075
	! do i=1,dn*dn
		! dt1(i,i)=0.08
	! enddo
	n1p=nf0/2
	n2p=nf0/2
	n1=0
	n2=0
	im(1)=dn/2*dn-dn/2
	im(2)=im(1)
	cc=0
	call cpu_time(time1)
ex:	do while(any(abs(n1p-n1)>5e-5).or.any(abs(n2p-n2)>5e-5))
		! f1=dt-dt1
		dt=dt1
		n1=n1p
		n2=n2p
		nf=0
		! s0=sa
		! if(cc<=1) then
			sa=-1
			sb=1
		! endif
		do while(abs(nf-nf0)>5e-5)
			sp=0.5*(sa+sb)
			h=0
			do i=1,dn*dn
				! nearest neighbor pairs
				call nnp(i,0,np)
				do j=1,4
					h(2*i-1:2*i,2*np(j)-1:2*np(j))=h(2*i-1:2*i,2*np(j)-1:2*np(j))+reshape((/cmplx(-t,0d0,8),dt(i,np(j)),dt(i,np(j))&
						,cmplx(t,0d0,8)/),(/2,2/))
				enddo
				! next-nearest neighbor pairs
				call nnp(i,1,np)
				do j=1,4
					h(2*i-1:2*i,2*np(j)-1:2*np(j))=h(2*i-1:2*i,2*np(j)-1:2*np(j))+reshape((/-tp,0.0d0,0.0d0,tp/),(/2,2/))
				enddo
				! on site
				h(2*i-1:2*i,2*i-1:2*i)=h(2*i-1:2*i,2*i-1:2*i)+reshape((/du*n2(i)-sp,0.0d0,0.0d0,-du*n1(i)+sp/),(/2,2/))
				if(any(im==i)) then
					h(2*i-1:2*i,2*i-1:2*i)=h(2*i-1:2*i,2*i-1:2*i)+reshape((/eff,0.0d0,0.0d0,eff/),(/2,2/))
				endif
			enddo
			! do i=1,dn*dn
				! do j=1,dn*dn
					! ! nearest neighbor pairs
					! call nnp(i,0,np)
					! if(any(np==j)) then
						! h(2*i-1:2*i,2*j-1:2*j)=h(2*i-1:2*i,2*j-1:2*j)+reshape((/cmplx(-t,0),dt(i,j),dt(i,j),cmplx(t,0)/),(/2,2/))
						! cycle
					! endif
					! ! next-nearest neighbor pairs
					! call nnp(i,1,np)
					! if(any(np==j)) then
						! h(2*i-1:2*i,2*j-1:2*j)=h(2*i-1:2*i,2*j-1:2*j)+reshape((/-tp,0.0,0.0,tp/),(/2,2/))
						! cycle
					! endif
					! ! on site
					! if(i==j) then
						! h(2*i-1:2*i,2*i-1:2*i)=h(2*i-1:2*i,2*i-1:2*i)+reshape((/du*n2(i)-sp,0.0,0.0,-du*n1(i)+sp/),&
							! (/2,2/))
						! if(any(im==i)) then
							! h(2*i-1:2*i,2*i-1:2*i)=h(2*i-1:2*i,2*i-1:2*i)+reshape((/eff,0.0,0.0,eff/),(/2,2/))
						! endif
					! endif
				! enddo
			! enddo
			ap=0
			do i=1,2*dn*dn
				do j=i,2*dn*dn
					ap(i,j)=h(i,j)
				enddo
			enddo
			call zheevd("v","u",2*dn*dn,ap,2*dn*dn,e,work,4*dn*dn+(2*dn*dn)**2,rwork,&
				1+10*dn*dn+2*(2*dn*dn)**2,iwork,3+10*dn*dn,info)
			! call dspev('v','u',2*dn*dn,ap,e,z,2*dn*dn,work,info)
			if(info/=0) then
				exit ex
			endif
			do i=1,2*dn*dn
				if(mod(i,2)==1) then
					u((i+1)/2,:)=ap(i,:)
				else
					v(i/2,:)=ap(i,:)
				endif
			enddo
			! write(*,"(e15.6)")e
			nf=0
			do i=1,dn*dn
				do n=1,2*dn*dn
					nf=nf+(u(i,n)*conjg(u(i,n))*f(e(n))+v(i,n)*conjg(v(i,n))*(1-f(e(n))))
				enddo
			enddo
			nf=1.0/(dn*dn)*nf
			if (nf<nf0) then
				sa=sp
			else
				sb=sp
			endif
			! write(*,"(3e11.3)")nf,sa,sb
		enddo
		n1p=0
		n2p=0
		do i=1,dn*dn
			do n=1,2*dn*dn
				n1p(i)=n1p(i)+u(i,n)*conjg(u(i,n))*f(e(n))
				n2p(i)=n2p(i)+v(i,n)*conjg(v(i,n))*(1-f(e(n)))
			enddo
		enddo
		dt1=0
		do i=1,dn*dn
			call nnp(i,0,np)
			do j=1,4
				do n=1,2*dn*dn
					dt1(i,np(j))=dt1(i,np(j))+(u(i,n)*conjg(v(np(j),n))+conjg(v(i,n))*u(np(j),n))*tanh(0.5*bt*e(n))
				enddo
			enddo
		enddo
		dt1=0.25*dv*dt1
		! write(*,"(2e15.6,f15.12)")n1(1),n1p(1),sa
		! if(cc>=1) then
			! sb=s0
			! sa=sa+sa-s0
			! if(sa>sb) then
				! sb=sa
				! sa=s0
			! endif
		! endif
		! cc=cc+1
		rewind(10)
		write(10,400)((n1(dn*i+j)+n2(dn*i+j),j=1,dn),i=0,dn-1)
		! do i=1,dn*dn
			! call nnp(i,0,np)
			! dtf(i)=(dt(i,np(1))+dt(i,np(3))-dt(i,np(2))-dt(i,np(4)))/4.0
		! enddo
		! write(10,200)((dtf(24*i+j),j=1,24),i=0,23)
		! f2=dt-dt1
		! if(all(f1+f2<1e-5)) then
			! error=1
			! exit
		! endif
	enddo ex
	call cpu_time(time2)
	rewind(10)
	do i=1,dn*dn
		call nnp(i,0,np)
		dtf(i)=(dt(i,np(1))+dt(i,np(3))-dt(i,np(2))-dt(i,np(4)))/4.0
	enddo
	! write(10,300)(dt(i,:),i=1,dn*dn)
	write(10,200)((dtf(dn*i+j),j=1,dn),i=0,dn-1)
	write(10,400)((n1(dn*i+j)+n2(dn*i+j),j=1,dn),i=0,dn-1)
	write(10,400)(((-1)**(i+j)*0.5*(n1(dn*i+j)-n2(dn*i+j)),j=1,dn),i=0,dn-1)
	write(10,*)"the beta is",bt,"v=",dv,"u=",du,"t'=",tp,"heff=",eff,"nf=",nf0,"impure",mod(1,im(2)-im(1)),"n=",dn
	write(10,*)"error=",error,"info=",info
	write(10,*)"the tatal runtime is:",time2-time1
	if(info/=0) then
		error=mod(info,2*dn*dn+1)
		do j=1,error+1
			write(20,500)(h(i,j),i=1,error+1)
		enddo
500		format(1000000("(",e10.3,",",e10.3,")  "))
	endif
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

	
	

	
