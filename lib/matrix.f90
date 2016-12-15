module M_matrix
	use lapack95, only : getrf, getri, heevd, heev, heevx, heevr, sysv
	use ifport, ifort_qsort => qsort
	use M_const
	implicit none
	interface diag
		module procedure mdiag, ndiag, mcdiag, ncdiag
	end interface
	interface crsmv
		module procedure c_crsmv,r_crsmv
	end interface
	interface mat_inv
		module procedure cmat_inv,rmat_inv
	end interface
	interface outprod
		module procedure outprod_i,outprod_r,outprod_c
	end interface
contains
	subroutine det_ratio_row(c,ci,A,iA,pb)
		complex(8) :: c(:),A(:,:),iA(:,:),pb
		integer :: ci,i,n
		pb=0d0
		n=size(A,1)
		c=c-A(ci,:)
		!!$OMP PARALLEL DO REDUCTION(+:pb)
		do i=1,n
			pb=pb+c(i)*iA(i,ci)
		enddo
		!!$OMP END PARALLEL DO
		pb=pb+1d0
	end subroutine
	subroutine det_ratio_col(c,cj,A,iA,pb)
		complex(8) :: c(:),A(:,:),iA(:,:),pb
		integer :: cj,i,n
		pb=0d0
		n=size(A,1)
		c=c-A(:,cj)
		!!$OMP PARALLEL DO REDUCTION(+:pb)
		do i=1,n
			pb=pb+iA(cj,i)*c(i)
		enddo
		!!$OMP END PARALLEL DO
		pb=pb+1
	end subroutine
	subroutine det_ratio_rowcol(c,ij,A,iA,iY,pb)
		complex(8) :: c(:,:),A(:,:),iA(:,:),pb,iY(:,:)
		integer :: ij(:),i,j,n
		iY=0d0
		n=size(A,1)
		c(:,1)=c(:,1)-A(ij(1),:)
		c(ij(2),1)=0d0
		c(:,2)=c(:,2)-A(:,ij(2))
		iY(2,1)=-iA(ij(2),ij(1))
		!!$OMP PARALLEL DO REDUCTION(+:iY)
		do i=1,n
			iY(1,1)=iY(1,1)+iA(ij(2),i)*c(i,2)
			iY(2,2)=iY(2,2)+c(i,1)*iA(i,ij(1))
			do j=1,n
				iY(1,2)=iY(1,2)-c(i,1)*iA(i,j)*c(j,2)
			enddo
		enddo
		!!$OMP END PARALLEL DO
		iY(1,1)=iY(1,1)+1d0
		iY(2,2)=iY(2,2)+1d0
		pb=iY(1,1)*iY(2,2)-iY(1,2)*iY(2,1)
		iY=iY/pb
	end subroutine
	subroutine det_ratio_tworow(c,ij,A,iA,iY,pb)
		complex(8) :: c(:,:),A(:,:),iA(:,:),pb,iY(:,:)
		integer :: ij(:),i,j,n
		iY=0d0
		n=size(A,1)
		c(:,1)=c(:,1)-A(ij(1),:)
		c(:,2)=c(:,2)-A(ij(2),:)
		!!$OMP PARALLEL DO REDUCTION(+:iY)
		do i=1,n
			iY(1,1)=iY(1,1)+iA(i,ij(2))*c(i,2)
			iY(2,2)=iY(2,2)+c(i,1)*iA(i,ij(1))
			iY(1,2)=iY(1,2)-c(i,1)*iA(i,ij(2))
			iY(2,1)=iY(2,1)-c(i,2)*iA(i,ij(1))
		enddo
		!!$OMP END PARALLEL DO
		iY(1,1)=iY(1,1)+1d0
		iY(2,2)=iY(2,2)+1d0
		pb=iY(1,1)*iY(2,2)-iY(1,2)*iY(2,1)
		iY=iY/pb
	end subroutine
	subroutine inv_update_row(c,ci,pb,A,iA)
		complex(8) :: c(:),A(:,:),iA(:,:),pb,ipb
		complex(8) :: tmp1(size(A,1)),tmp2(size(A,1),size(A,1))
		integer :: ci,i,j,n
		ipb=1d0/pb
		n=size(iA,1)
		!allocate(tmp1(n),tmp2(n,n))
		tmp1=0d0
		A(ci,:)=A(ci,:)+c
		!!$OMP PARALLEL DO REDUCTION(+:tmp1)
		do i=1,n
			do j=1,n
				tmp1(i)=tmp1(i)+c(j)*iA(j,i)
			enddo
		enddo
		!!$OMP END PARALLEL DO
		!!$OMP PARALLEL DO REDUCTION(+:tmp2)
		do i=1,n
			do j=1,n
				tmp2(i,j)=iA(i,j)-ipb*iA(i,ci)*tmp1(j)
			enddo
		enddo
		!!$OMP END PARALLEL DO
		iA=tmp2
	end subroutine
	subroutine inv_update_col(c,cj,pb,A,iA)
		complex(8) :: c(:),A(:,:),iA(:,:),pb,ipb
		complex(8) :: tmp1(size(A,1)),tmp2(size(A,1),size(A,1))
		!complex(8), allocatable :: tmp1(:),tmp2(:,:)
		integer :: cj,i,j,n
		ipb=1d0/pb
		n=size(iA,1)
		!allocate(tmp1(n),tmp2(n,n))
		tmp1=0d0
		A(:,cj)=A(:,cj)+c
		!!$OMP PARALLEL DO REDUCTION(+:tmp1)
		do i=1,n
			do j=1,n
				tmp1(i)=tmp1(i)+iA(i,j)*c(j)
			enddo
		enddo
		!!$OMP END PARALLEL DO
		!!$OMP PARALLEL DO REDUCTION(+:tmp2)
		do i=1,n
			do j=1,n
				tmp2(i,j)=iA(i,j)-ipb*tmp1(i)*iA(cj,j)
			enddo
		enddo
		!!$OMP END PARALLEL DO
		iA=tmp2
	end subroutine
	subroutine inv_update_rowcol(c,ij,iY,A,iA)
		complex(8) :: c(:,:),A(:,:),iA(:,:),iY(:,:)
		complex(8) :: tmp1(size(A,1),2),tmp2(size(A,1),size(A,1))
		integer :: ij(:),i,j,n
		n=size(iA,1)
		!allocate(tmp1(n,2),tmp2(n,n))
		tmp1=0d0
		A(ij(1),:)=A(ij(1),:)+c(:,1)
		A(:,ij(2))=A(:,ij(2))+c(:,2)
		!!$OMP PARALLEL DO REDUCTION(+:tmp1)
		do i=1,n
			do j=1,n
				tmp1(i,1)=tmp1(i,1)+c(j,1)*iA(j,i)
				tmp1(i,2)=tmp1(i,2)+iA(i,j)*c(j,2)
			enddo
		enddo
		!!$OMP END PARALLEL DO
		!!$OMP PARALLEL DO REDUCTION(+:tmp2)
		do i=1,n
			do j=1,n
				tmp2(i,j)=iA(i,j)-((iA(i,ij(1))*iY(1,1)+tmp1(i,2)*iY(2,1))*tmp1(j,1)+&
					(iA(i,ij(1))*iY(1,2)+tmp1(i,2)*iY(2,2))*iA(ij(2),j))
			enddo
		enddo
		!!$OMP END PARALLEL DO
		iA=tmp2
	end subroutine
	subroutine inv_update_tworow(c,ij,iY,A,iA)
		complex(8) :: c(:,:),A(:,:),iA(:,:),iY(:,:)
		complex(8) :: tmp1(size(A,1),2),tmp2(size(A,1),size(A,1))
		integer :: ij(:),i,j,n
		n=size(iA,1)
		!allocate(tmp1(n,2),tmp2(n,n))
		tmp1=0d0
		A(ij(1),:)=A(ij(1),:)+c(:,1)
		A(ij(2),:)=A(ij(2),:)+c(:,2)
		!!$OMP PARALLEL DO REDUCTION(+:tmp1)
		do i=1,n
			do j=1,n
				tmp1(i,1)=tmp1(i,1)+c(j,1)*iA(j,i)
				tmp1(i,2)=tmp1(i,2)+c(j,2)*iA(j,i)
			enddo
		enddo
		!!$OMP END PARALLEL DO
		!!$OMP PARALLEL DO REDUCTION(+:tmp2)
		do i=1,n
			do j=1,n
				tmp2(i,j)=iA(i,j)-((iA(i,ij(1))*iY(1,1)+iA(i,ij(2))*iY(2,1))*tmp1(j,1)+&
					(iA(i,ij(1))*iY(1,2)+iA(i,ij(2))*iY(2,2))*tmp1(j,2))
			enddo
		enddo
		!!$OMP END PARALLEL DO
		iA=tmp2
	end subroutine
	subroutine rmat_inv(A,info)
		real(8) :: A(:,:)
		integer :: ipiv(size(A,1)),info1
		integer, optional :: info
		if(present(info)) then
			info=0
		endif
		call getrf(A,ipiv,info1)
		if(info1/=0) then
			if(present(info)) then
				info=info1
				return
			else
				write(*,*)"inverse matrix err1"
				write(*,*)RAISEQQ(SIG$ABORT)
			endif
		endif
		call getri(A,ipiv,info1)
		if(info1/=0) then
			if(present(info)) then
				info=info1
				return
			else
				write(*,*)"inverse matrix err1"
				write(*,*)RAISEQQ(SIG$ABORT)
			endif
		endif
	end subroutine
	subroutine cmat_inv(A,info)
		complex(8) :: A(:,:)
		integer :: ipiv(size(A,1)),info1
		integer, optional :: info
		if(present(info)) then
			info=0
		endif
		call getrf(A,ipiv,info1)
		if(info1/=0) then
			if(present(info)) then
				info=info1
				return
			else
				write(*,*)"inverse matrix err1"
				write(*,*)RAISEQQ(SIG$ABORT)
			endif
		endif
		call getri(A,ipiv,info1)
		if(info1/=0) then
			if(present(info)) then
				info=info1
				return
			else
				write(*,*)"inverse matrix err2"
				write(*,*)RAISEQQ(SIG$ABORT)
			endif
		endif
	end subroutine
	function det(A,info)
		complex(8) :: A(:,:),A_(size(A,1),size(A,2))
		integer, optional :: info
		complex(8) :: det
		integer :: i,ipiv(size(A,1))
		A_=A
		if(present(info)) then
			call getrf(A_,ipiv,info)
			if(info/=0) then
				stop "determinant err"
			endif
		else
			call getrf(A_,ipiv)
		endif
		det=1d0
		do i=1,size(A,1)
			if(ipiv(i)/=i) then
				det=-det*A_(i,i)
			else
				det=det*A_(i,i)
			endif
		enddo
	end function
	function mdiag(a)
		real(8) :: a(:) 
		real(8) :: mdiag(size(a),size(a))
		integer :: i
		mdiag=0d0
		do i=1,size(a)
			mdiag(i,i)=a(i)
		enddo
	end function
	function ndiag(a,n)
		real(8) :: a
		integer :: n
		real(8) :: ndiag(n,n)
		integer :: i
		ndiag=0d0
		do i=1,n
			ndiag(i,i)=a
		enddo
	end function
	function mcdiag(a)
		complex(8) :: a(:) 
		complex(8) :: mcdiag(size(a),size(a))
		integer :: i
		mcdiag=0d0
		do i=1,size(a)
			mcdiag(i,i)=a(i)
		enddo
	end function
	function ncdiag(a,n)
		complex(8) :: a
		integer :: n
		complex(8) :: ncdiag(n,n)
		integer :: i
		ncdiag=0d0
		do i=1,n
			ncdiag(i,i)=a
		enddo
	end function
	function Tr(A,B)
		complex(8) :: A(:,:),B(:,:),Tr
		integer :: n,i,j
		n=size(A,1)
		Tr=0d0
		do i=1,n
			do j=1,n
				Tr=Tr+A(i,j)*B(j,i)
			enddo
		enddo
	end function
	function check_diag(A,U,E,err)
		complex(8) :: A(:,:),U(:,:)
		real(8) :: E(:)
		real(8), optional :: err
		real(8) :: err_
		logical :: check_diag
		check_diag=.true.
		if(present(err)) then
			err_=err
		else
			err_=1d-6
		endif
		if(any(abs((matmul(transpose(conjg(U)),matmul(A,U))-diag(E)))>err_)) then
			check_diag=.false.
		endif
	end function
	!subroutine mat_inv(a)
		!writen by myself, not quite work
		!implicit none
		!complex(8) :: a(:,:)
		!complex(8), allocatable :: ctmp(:)
		!real(8) :: tmp
		!integer :: i,j,k,n,m,l,mx(size(a,1),2)
		!l=size(a,1)
		!allocate(ctmp(l))
		!do i=1,l
			!tmp=0d0
			!do j=i,l
				!do k=i,l
					!if(tmp<abs(a(j,k))) then
						!tmp=abs(a(j,k))
						!mx(i,1)=j
						!mx(i,2)=k
					!endif
				!enddo
			!enddo
			!ctmp=a(i,:)
			!a(i,:)=a(mx(i,1),:)
			!a(mx(i,1),:)=ctmp
			!ctmp=a(:,i)
			!a(:,i)=a(:,mx(i,2))
			!a(:,mx(i,2))=ctmp
			!a(i,i)=1d0/a(i,i)
			!do j=1,l-1
				!n=mod(i+j-1,l)+1
				!do k=1,l-1
					!m=mod(k+i-1,l)+1
					!a(i,m)=a(i,i)*a(i,m)
					!a(n,m)=a(n,m)-a(n,i)*a(i,m)
				!enddo
				!a(n,i)=-a(n,i)*a(i,i)
			!enddo
		!enddo
		!do i=l,1,-1
			!ctmp=a(:,i)
			!a(:,i)=a(:,mx(i,1))
			!a(:,mx(i,1))=ctmp
			!ctmp=a(i,:)
			!a(i,:)=a(mx(i,2),:)
			!a(mx(i,2),:)=ctmp
		!enddo
	!end subroutine
	subroutine r_conjgrad(A,b,x)
		real(8) :: A(:,:),b(:),x(:),r(size(x)),p(size(x)),Ap(size(x)),al,tmp0,tmp1,cvg=1d-10
		r=b-matmul(A,x)
		p=r
		tmp0=dot_product(r,r)
		do
			Ap=matmul(A,p)
			al=tmp0/dot_product(p,Ap)
			x=x+al*p
			r=r-al*Ap
			tmp1=dot_product(r,r)
			if(tmp1<cvg) then
				exit
			endif
			p=r+tmp1/tmp0*p
			tmp0=tmp1
		enddo
	end subroutine
	subroutine r_crsmv(va,ja,ia,x,y)
		integer :: i,j
		integer :: ja(:),ia(:)
		real(8) :: va(:),y(:),x(:)
		do i=1,size(x)
			y(i)=0
			do j=ia(i),ia(i+1)-1
				y(i)=y(i)+va(j)*x(ja(j))
			enddo
		enddo
	end subroutine
	subroutine c_crsmv(va,ja,ia,x,y)
		integer :: i,j
		integer :: ja(:),ia(:)
		complex(8) :: va(:),y(:),x(:)
		do i=1,size(x)
			y(i)=0
			do j=ia(i),ia(i+1)-1
				y(i)=y(i)+va(j)*x(ja(j))
			enddo
		enddo
	end subroutine
	subroutine crs(v,va,ja,ia)
		integer :: n
		complex(8) :: va(:),v(:,:)
		integer :: ja(:),ia(size(v,1)+1)
		integer ::  i,j
		ia(1)=1
		n=1
		do i=1,size(v,1)
			do j=1,size(v,2)
				if((real(v(i,j))**2+imag(v(i,j))**2)>1d-18) then
					va(n)=v(i,j)
					ja(n)=j
					n=n+1
				endif
			enddo
			ia(i+1)=n
		enddo
	end subroutine
	subroutine diag2(A,E,info)
		complex(8) :: A(:,:)
		real(8) :: E(:)
		integer, optional :: info
		real(8) :: tmp
		if(abs(A(2,1))<1d-7) then
			E=real((/A(1,1),A(2,2)/))
			if(E(1)<=E(2)) then
				A(1,1)=1d0
				A(2,2)=1d0
			else
				E=real((/A(2,2),A(1,1)/))
				A(1,2)=1d0
				A(2,1)=1d0
				A(1,1)=0d0
				A(2,2)=0d0
			endif
		else
			tmp=sqrt((A(1,1)-A(2,2))**2+4d0*A(1,2)*A(2,1))
			E=0.5d0*(A(1,1)+A(2,2)+(/-tmp,tmp/))
			A(1,:)=E-A(2,2)
			A(2,:)=A(2,1)
			A(:,1)=A(:,1)/sqrt(A(1,1)**2+A(2,1)*conjg(A(2,1)))
			A(:,2)=A(:,2)/sqrt(A(1,2)**2+A(2,2)*conjg(A(2,2)))
		endif
	end subroutine
	subroutine diag4(A,E,info)
		complex(8) :: A(:,:)
		real(8) :: E(:)
		integer, optional :: info
		complex(8) :: tmp2(2,2),tmp4(4,4)
		tmp4=0d0
		call diag2(A(:2,:2),E(:2))
		A(4,4)=A(1,1)
		A(3,3)=A(2,2)
		A(3,4)=A(2,1)
		A(4,3)=A(1,2)

		tmp4(1,1)=E(1)
		tmp4(2,2)=-E(1)
		tmp4(1,2)=A(1,4)
		tmp4(2,1)=A(4,1)
		call diag2(tmp4(:2,:2),E(3:4))
		tmp4(4,1)=tmp4(2,1)
		tmp4(1,4)=tmp4(1,2)
		tmp4(4,4)=tmp4(2,2)

		tmp4(2,2)=E(2)
		tmp4(3,3)=-E(2)
		tmp4(2,3)=A(2,3)
		tmp4(3,2)=A(3,2)
		call diag2(tmp4(2:3,2:3),E(1:2))

		E=(/E(3),E(1),E(2),E(4)/)
		A(1:2,3:4)=0d0
		A(3:4,1:2)=0d0
		tmp4(2,1)=0d0
		tmp4(1,2)=0d0
		A=matmul(A,tmp4)
	end subroutine
	function outprod_i(A,B) result(rt)
		integer :: A(:),B(:)
		integer :: rt(size(A)*size(B))
		integer :: i,n
		n=size(A)
		do i=1,size(B)
			rt(n*(i-1)+1:n*i)=A*B(i)
		enddo
	end function
	function outprod_r(A,B) result(rt)
		real(8) :: A(:),B(:)
		real(8) :: rt(size(A)*size(B))
		integer :: i,n
		n=size(A)
		do i=1,size(B)
			rt(n*(i-1)+1:n*i)=A*B(i)
		enddo
	end function
	function outprod_c(A,B) result(rt)
		complex(8) :: A(:),B(:)
		complex(8) :: rt(size(A)*size(B))
		integer :: i,n
		n=size(A)
		do i=1,size(B)
			rt(n*(i-1)+1:n*i)=A*B(i)
		enddo
	end function
	!subroutine mat_diag(H,E,info)
		!complex(8) :: H(:,:),tmp(size(H,1),size(H,2))
		!real(8) :: E(:)
		!integer, optional :: info
		!if(present(info)) then
			!tmp=H
		!endif
		!select case(size(E))
		!case(2)
			!call diag2(H,E)
		!case(4)
			!call diag4(H,E)
		!case(100:)
			!call heevd(H,E,"V")
		!case default
			!call heev(H,E,"V")
		!end select
		!if(present(info)) then
			!if(.not.check_diag(tmp,H,E)) then
				!write(*,*)"diag4 err, A is"
				!write(*,*)"real part"
				!write(*,"(4es12.4)")real(tmp)
				!write(*,*)"imag part"
				!write(*,"(4es12.4)")imag(tmp)
				!write(*,*)"E is"
				!write(*,"(4es12.4)")E
				!write(*,*)"UAU is"
				!write(*,"(4es12.4)")real(matmul(transpose(conjg(H)),matmul(tmp,H)))
				!stop
			!endif
		!endif
	!end subroutine
end module
