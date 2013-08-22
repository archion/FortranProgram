module global
	implicit none
	real(8),save :: dv=1.0,nf0=0.85,du=0d0,t=1,tp=0d0,eff=1.0
	integer, parameter :: dn=24
	real(8), parameter :: pi=3.14159
	integer,save,allocatable :: nn(:,:),pbc(:)
end module
program main
	use global
	implicit none
	character(24) :: now
	integer :: i,j,k=1,cutoff,n,n_mod(3),im(2)=-1,info=0,np(4),error=0,ini=0
	real(8), allocatable :: n1(:),n2(:),n1p(:),n2p(:),e(:),dtr(:,:),hhn(:,:),een(:,:)
	real(8), allocatable :: h(:,:),h_1(:,:),dt(:,:),dt1(:,:),dtf(:)
	real(8) :: en,sa,sb,sp,s0,nf,cache,a,b,ac_ba,l=4
	allocate(h(2*dn*dn,2*dn*dn),h_1(2*dn*dn,2*dn*dn),n1(dn*dn),n2(dn*dn),n1p(dn*dn),n2p(dn*dn),dt(dn*dn,dn*dn),&
		dt1(dn*dn,dn*dn),dtr(dn*dn,dn*dn),dtf(dn*dn),nn(dn,dn),pbc(0:dn+1),hhn(0:2,2*dn*dn),&
		een(0:2,2*dn*dn))
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
	dt=0
	dt1=0.5
	! call init_random_seed()
	! call random_number(dtr)
	! dt1=dtr
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
	! im(2)=im(1)-4
	cutoff=1000
	! call cpu_time(time1)
! ex:	do while(any(abs(n1p-n1)>5e-5).or.any(abs(n2p-n2)>5e-5))
ex:	do while(any(abs(dt1-dt)>5e-5).and.(ini<25))
		ini=ini+1
		call fdate(now)
		write(*,*)dt(1,2),dt1(1,2),now
		dt=dt1
		! n1=n1p
		! n2=n2p
		nf=0
		sa=0
		sb=1
		! do while(abs(nf-nf0)>5e-5)
			! write(*,*)sa,sb
			! sp=0.5*(sa+sb)
			sp=-1.5*t
			dv=-2.2*t
			eff=0.3*t
			a=10
			b=0
			ac_ba=acos(-b/a)
			h=0
			do i=1,dn*dn
				! nearest neighbor pairs
				call nnp(i,0,np)
				do j=1,4
					h(2*i-1:2*i,2*np(j)-1:2*np(j))=h(2*i-1:2*i,2*np(j)-1:2*np(j))+reshape((/-t,dt(i,np(j)),dt(i,np(j)),t/),(/2,2/))
				enddo
				! next-nearest neighbor pairs
				! call nnp(i,1,np)
				! do j=1,4
					! h(2*i-1:2*i,2*np(j)-1:2*np(j))=h(2*i-1:2*i,2*np(j)-1:2*np(j))+reshape((/-tp,0.0d0,0.0d0,tp/),(/2,2/))
				! enddo
				! on site
				h(2*i-1:2*i,2*i-1:2*i)=h(2*i-1:2*i,2*i-1:2*i)+reshape((/-sp,0.0d0,0.0d0,sp/),(/2,2/))
				if(any(im==i)) then
					h(2*i-1:2*i,2*i-1:2*i)=h(2*i-1:2*i,2*i-1:2*i)+reshape((/eff,0.0d0,0.0d0,eff/),(/2,2/))
				endif
			enddo
			! do i=1,2*dn*dn
				! do j=1,2*dn*dn
					! if((h(i,j)-h(j,i))>5e-5) then
						! write(*,*)i,j,"error"
					! endif
				! enddo
			! enddo
			h_1=(h-b)/a
			n1p=0
			n2p=0
			dt1=0
			!$omp parallel do firstprivate(een,hhn) schedule(guided,1) num_threads(32)
			do i=1,dn*dn
				een=0
				hhn=0
				een(0,2*i-1)=1d0
				hhn(0,2*i)=1d0
				! n1p(i)=n1p(i)+een(0,2*i-1)*(1-ac_ba/pi)
				! n2p(i)=n2p(i)+hhn(0,2*i)*(1-ac_ba/pi)
				! call crsmv(h_1,een(0,:),een(1,:),2*dn*dn,2*dn*dn)
				call crsmv(h_1,hhn(0,:),hhn(1,:),2*dn*dn)
				! n1p(i)=n1p(i)+een(1,2*i-1)*2*sin(ac_ba)/pi
				! n2p(i)=n2p(i)+hhn(1,2*i)*2*sin(ac_ba)/pi
				call nnp(i,0,np)
				do j=1,4
					dt1(i,np(j))=dt1(i,np(j))+hhn(1,2*np(j)-1)*sin(ac_ba)*sinh(l*(1-1.0/real(cutoff)))/sinh(l)
				enddo
				do n=2,cutoff
					n_mod(1)=mod(n-2,3)
					n_mod(2)=mod(n-1,3)
					n_mod(3)=mod(n,3)
					! call crsmv(h_1,een(n_mod(2),:),een(n_mod(3),:),2*dn*dn,2*dn*dn)
					! een(n_mod(3),:)=2*een(n_mod(3),:)-een(n_mod(1),:)
					call crsmv(h_1,hhn(n_mod(2),:),hhn(n_mod(3),:),2*dn*dn)
					hhn(n_mod(3),:)=2*hhn(n_mod(3),:)-hhn(n_mod(1),:)
					! n1p(i)=n1p(i)+een(n_mod(3),2*i-1)*2*sin(nc*ac_ba)/nc/pi
					! n2p(i)=n2p(i)+hhn(n_mod(3),2*i)*2*sin(nc*ac_ba)/nc/pi
					call nnp(i,0,np)
					do j=1,4
						dt1(i,np(j))=dt1(i,np(j))+hhn(n_mod(3),2*np(j)-1)*sin(n*ac_ba)/n*sinh(l*(1-real(n)/real(cutoff)))/sinh(l)
					enddo
				enddo
			enddo
			!$omp end parallel do
			! n2p=1-n2p
			! if(any(n2p<0)) then
				! exit ex
			! endif
			! nf=0
			! do i=1,dn*dn
				! nf=nf+n1p(i)+n2p(i)
			! enddo
			! nf=nf/(dn*dn)
			! write(*,"(e15.6)")nf
			! if (nf<nf0) then
				! sa=sp
			! else
				! sb=sp
			! endif
			! ! write(*,"(3e11.3)")nf,sa,sb
		! enddo
		dt1=-2*dv*dt1/pi
		! write(*,"(2e15.6,f15.12)")n1(1),n1p(1),sa
		! rewind(10)
		! write(10,400)((n1(dn*i+j)+n2(dn*i+j),j=1,dn),i=0,dn-1)
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
		rewind(10)
		do i=1,dn*dn
			call nnp(i,0,np)
			dtf(i)=(dt(i,np(1))+dt(i,np(3))+dt(i,np(2))+dt(i,np(4)))/4.0
		enddo
		write(10,400)((dtf(dn*i+j),j=1,dn),i=0,dn-1)
	enddo ex
	! call cpu_time(time2)
	rewind(10)
	do i=1,dn*dn
		call nnp(i,0,np)
		dtf(i)=(dt(i,np(1))+dt(i,np(3))-dt(i,np(2))-dt(i,np(4)))/4.0
	enddo
	! write(10,300)(dt(i,:),i=1,dn*dn)
	write(10,400)((dtf(dn*i+j),j=1,dn),i=0,dn-1)
	! write(10,400)((n1(dn*i+j)+n2(dn*i+j),j=1,dn),i=0,dn-1)
	! write(10,400)(((-1)**(i+j)*0.5*(n1(dn*i+j)-n2(dn*i+j)),j=1,dn),i=0,dn-1)
	! write(10,*)"the tatal runtime is:",time1
	write(10,"(16e12.3)")(dt(i,:),i=1,dn*dn)
	write(10,"(32e12.3)")(h(i,:),i=1,2*dn*dn)
	! if(info/=0) then
		! error=mod(info,2*dn*dn+1)
		! do j=1,error+1
			! write(20,500)(h(i,j),i=1,error+1)
		! enddo
! 500		format(1000000("(",e10.3,",",e10.3,")  "))
	! endif
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

	
	

	
! subroutine crsmv(a,x,y,m)
	! implicit none
	! integer :: i,job(6),info,m
	! real(8) :: a(m,m),x(*),y(*)
	! integer :: ja(10*m),ia(m+1)
	! real(8) :: acsr(10*m)
	! ! integer, allocatable :: ja(:),ia(:)
	! ! real(8), allocatable :: acsr(:)
	! job=(/0,1,1,2,m*m,1/)
	! ! allocate(acsr(m*m),ja(m*m),ia(m+1))
	! ! allocate(acsr(m*m))
	! call mkl_ddnscsr(job, m, m, a, m, acsr, ja, ia, info)
	! call mkl_dcsrgemv("n", m, acsr, ia, ja, x, y)
! end
subroutine crsmv(a,x,y,m)
	implicit none
	integer :: i=0,j=0,k1,k2,m
	integer :: ja(10*m+1),ia(m+1)
	real(8) :: va(10*m)
	real(8) :: a(m,m),x(*),y(*)
	k1=1
	k2=1
	do i=1,m
		ia(k2)=k1
		k2=k2+1
		do j=1,m
			if(abs(a(i,j))>5e-5) then
				va(k1)=a(i,j)
				ja(k1)=j
				k1=k1+1
			endif
		enddo
	enddo
	ia(k2)=k1
	do i=1,m
		y(i)=0
		do j=ia(i),ia(i+1)-1
			y(i)=y(i)+va(j)*x(ja(j))
		enddo
	enddo
end
! subroutine crsmv(a,x,y,m)
	! implicit none
	! integer :: i,job(6),info,m
	! real(8) :: a(m,m),x(*),y(*)
	! integer, allocatable :: ja(:),ia(:)
	! real(8), allocatable :: acsr(:)
	! job=(/0,1,1,2,m*m,1/)
	! allocate(acsr(m*m),ja(m*m),ia(m+1))
	! call mkl_ddnscsr(job, m, m, a, m, acsr, ja, ia, info)
	! call mkl_dcsrgemv("n", m, acsr, ia, ja, x, y)
! end
