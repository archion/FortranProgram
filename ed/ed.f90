module global
	implicit none
	real(8),save :: v(0:37,36)=0d0,t0=1d0,u=6d0
end module
program main
	use global
	implicit none
	real(8) :: nv(36)=0d0,ap(36)=0d0,bt(0:36)=0d0,cache,z(36,36),work(70)
	integer :: tf(36),a(6)=(/ 3,5,6,9,10,12 /),i,j,k,l,info,dg=0
	open(unit=10,file='../data/output.dat')
	do j=1,6
		do k=1,6
			tf((j-1)*6+k)=a(j)*16+a(k)
		enddo
	enddo
	! v(1,2)=1
	! call rn(v(1,:),cache)
	! call hemiltion(1,tf,nv)
	! write(*,"(f5.1f5.1b10.8)")(v(1,i),nv(i),tf(i),i=1,36)
	call random_seed
	call random_number(v(1,:))
	! call rn(v(1,:),cache)
	! call hemiltion(1,tf,nv)
	! write(*,"(f5.1f5.1b10.8)")(v(1,i),nv(i),tf(i),i=1,36)
	call rn(v(1,:),cache)
	! write(10,"(36e9.2)")(v(k,:),k=0,3)
	! write(*,*)v(1,:)
	do i=1,36
		call hemiltion(i,tf,nv)
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
			dg=dg+1
			! do while(.true.)
				call random_number(v(i+1,:))
				call ro(i+1)
				call rn(v(i+1,:),cache)
				! if(cache>1e-10) then
					! exit
				! endif
			! enddo
		endif
		write(10,*)dot_product(v(i+1,:),v(i-1,:)),dot_product(v(i+1,:),v(i-5,:)),bt(i)
		! write(10,*)dot_product(v(i+1,:),v(i,:)),bt(i),dot_product(nv(:),v(i+1,:))
		! write(10,"(36e9.2)")v(i+1,:)
	enddo
	write(10,"(e13.5,e13.5)")(ap(i),bt(i),i=1,36)
	call dsteqr('i',36,ap,bt(1:35),z,36,work,info)
	if(info==0) then
		write(10,"(e13.5)")ap
	endif

	close(10)
end program main


subroutine hemiltion(p,tf,nv)
	use global
	implicit none
	real(8) :: nv(36),cache
	integer :: tf(36),i,j,k,l=0,p,ctf,nos,uos
	nv=0
	do i=1,36
		do j=0,1
			do k=0,3
				ctf=tf(i)
				if(ibits(tf(i),j*4+k,1)/=ibits(tf(i),mod((k+1),4)+j*4,1)) then
					if(btest(tf(i),j*4+k)) then
						ctf=ibset(ctf,mod((j*4+k+1),4)+j*4)
						ctf=ibclr(ctf,j*4+k)
					else
						ctf=ibclr(ctf,mod((j*4+k+1),4)+j*4)
						ctf=ibset(ctf,j*4+k)
					endif
					call fd(tf,ctf,l)
					nv(l)=nv(l)+t0*v(p,i)
				endif
			enddo
		enddo
		uos=iand(ibits(tf(i),0,4),ibits(tf(i),4,4))
		nos=0
		do j=0,3
			if(btest(uos,j)) then
				nos=nos+1
			endif
		enddo
		nv(i)=nv(i)+v(p,i)*u*nos
	enddo
end subroutine hemiltion

subroutine fd(tf,ctf,l)
	implicit none
	integer :: tf(36),i,l,ctf
	do i=1,36
		if(tf(i)==ctf) then
			l=i
			exit
		endif
	enddo
end subroutine fd

subroutine rn(v,x)
	implicit none
	real(8) :: v(36),sun=0,x
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




