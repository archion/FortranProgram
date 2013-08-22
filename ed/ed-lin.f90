! this is a exactly diagonalize program for 2*2 lattice with t' using hubbard model
! for more information see h. q. lin, prb.35.3359
! writen by z. d. yu in 2012.5
! 1---2
! |   |
! 3---4
! n=2: 3,5,6,9,10,12
! n=1: 1,2,4,8
! n=3: 7,11,13,14
module global
	implicit none
	integer, parameter :: n=36,dn=2
	real(8),save :: v(0:n+1,n)=0d0,t0=1d0,u=5d0,tp=0
	integer,save :: tf(n),pbc(0:dn+1),nn(dn,dn)
end module
program main
	use global
	real(8) :: en(n),e
	integer :: i,j,k=1,a1(6)=(/ 3,5,6,9,10,12 /),a2(6)=(/ 3,5,6,9,10,12 /)
	open(unit=10,file='../data/doutput.dat')
	call init_random_seed()
	pbc=(/1,(i,i=1,dn),dn/)
	do i=1,dn
		do j=1,dn
			nn(i,j)=k
			k=k+1
		enddo
	enddo
	do j=1,size(a2)
		do k=1,size(a1)
			tf((j-1)*6+k)=a2(j)*(2**(2*dn))+a1(k)
		enddo
	enddo
	do i=1,25
		u=u+1
		call lanczos(en)
		e=en(1)
		write(10,*)en
		write(10,"(f4.1,e13.5)")u,en
	enddo
	write(10,"(e13.5)")en
	close(10)
end program main

subroutine lanczos(en)
	use global
	implicit none
	real(8) :: nv(n)=0d0,ap(n)=0d0,bt(0:n)=0d0,cache,z(n,n),work(2*n-2),en(n)
	integer :: i,j,k,l,info
	! v(1,4)=1
	! call rn(v(1,:),cache)
	! call hemiltion(1,tf,nv)
	! write(*,"(f5.1f5.1b10.8)")(v(1,i),nv(i),tf(i),i=1,n)
	call random_number(v(1,:))
	call rn(v(1,:),cache)
	! write(10,"(36e9.2)")(v(k,:),k=0,3)
	! write(*,*)v(1,:)
	do i=1,n
		call hemiltion(i,nv)
		ap(i)=dot_product(v(i,:),nv)
		! write(10,"(2e11.4)")(v(i,j),nv(j),j=1,36)
		! write(10,*)ap(i)
		v(i+1,:)=nv(:)-ap(i)*v(i,:)-bt(i-1)*v(i-1,:)
		if(mod(i,4)==0) then
			call ro(i+1)
		endif
		call rn(v(i+1,:),bt(i))
		! write(10,"(e11.3)")bt(i)
		if(bt(i)<1e-5) then
			bt(i)=0
			! do while(.true.)
				call random_number(v(i+1,:))
				call ro(i+1)
				call rn(v(i+1,:),cache)
				! if(cache>1e-10) then
					! exit
				! endif
			! enddo
		endif
		! write(10,*)dot_product(v(i+1,:),v(i-1,:)),dot_product(v(i+1,:),v(i-5,:)),bt(i)
		! write(10,*)dot_product(v(i+1,:),v(i,:)),bt(i),dot_product(nv(:),v(i+1,:))
		! write(10,"(36e9.2)")v(i+1,:)
	enddo
	! write(10,"(e13.5,e13.5)")(ap(i),bt(i),i=1,n)
	call dsteqr('i',n,ap,bt(1:n-1),z,n,work,info)
	! if(info==0) then
		! write(10,"(e13.5)")ap
	! endif
	en=ap
	! close(10)
end subroutine lanczos


subroutine hemiltion(p,nv)
	use global
	implicit none
	real(8) :: nv(n),cache
	integer :: np(4),i,j,k,ll,l=0,p,ctf,nos,uos
	nv=0
	do i=1,n
		do j=0,1
			do k=0,dn*dn-1
				if(btest(tf(i),j*4+k)) then
					call nnp(k+1,0,np)
					do ll=1,4
						ctf=tf(i)
						if((.not.(btest(tf(i),j*4+np(ll)-1))).and.(np(ll)/=k+1)) then
							ctf=ibset(ctf,j*4+np(ll)-1)
							ctf=ibclr(ctf,j*4+k)
							call fd(ctf,l)
							nv(l)=nv(l)+t0*v(p,i)
						endif
					enddo
					call nnp(k+1,1,np)
					do ll=1,4
						ctf=tf(i)
						if((.not.(btest(tf(i),j*4+np(ll)-1))).and.(np(ll)/=k+1)) then
							ctf=ibset(ctf,j*4+np(ll)-1)
							ctf=ibclr(ctf,j*4+k)
							call fd(ctf,l)
							nv(l)=nv(l)+tp*v(p,i)
						endif
					enddo
				endif
			enddo
		enddo
		uos=iand(ibits(tf(i),0,4),ibits(tf(i),4,4))
		nos=0
		do j=0,dn*dn-1
			if(btest(uos,j)) then
				nos=nos+1
			endif
		enddo
		nv(i)=nv(i)+v(p,i)*u*nos !+v(p,i)*4*tp*3
	enddo
end subroutine hemiltion

subroutine fd(ctf,l)
	use global, only : n,tf
	implicit none
	integer :: i,l,ctf
	do i=1,n
		if(tf(i)==ctf) then
			l=i
			exit
		endif
	enddo
end subroutine fd

subroutine rn(v,x)
	use global, only : n
	implicit none
	real(8) :: v(n),sun=0,x
	integer :: i
	sun=dot_product(v,v)
	x=sqrt(abs(sun))
	v=v/x
end subroutine rn

subroutine ro(i)
	use global
	implicit none
	integer :: i,j
	do j=1,i-1
		v(i,:)=v(i,:)-dot_product(v(j,:),v(i,:))*v(j,:)
	enddo
end

subroutine nnp(ii,p,jj)
	use global, only : pbc,nn,dn
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

