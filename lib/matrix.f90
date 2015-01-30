module M_matrix
	use lapack95, only : getrf, getri
	implicit none
	interface diag
		module procedure mdiag, ndiag
	end interface
contains
	subroutine det_ratio_row(c,ci,A,iA,pb)
		complex(8) :: c(:),A(:,:),iA(:,:),pb
		integer :: ci,i,n
		pb=0d0
		n=size(A,1)
		c=c-A(ci,:)
		do i=1,n
			pb=pb+c(i)*iA(i,ci)
		enddo
		pb=pb+1d0
	end subroutine
	subroutine det_ratio_col(c,cj,A,iA,pb)
		complex(8) :: c(:),A(:,:),iA(:,:),pb
		integer :: cj,i,n
		pb=0d0
		n=size(A,1)
		c=c-A(:,cj)
		do i=1,n
			pb=pb+iA(cj,i)*c(i)
		enddo
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
		do i=1,n
			iY(1,1)=iY(1,1)+iA(ij(2),i)*c(i,2)
			iY(2,2)=iY(2,2)+c(i,1)*iA(i,ij(1))
			do j=1,n
				iY(1,2)=iY(1,2)-c(i,1)*iA(i,j)*c(j,2)
			enddo
		enddo
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
		do i=1,n
			iY(1,1)=iY(1,1)+iA(i,ij(2))*c(i,2)
			iY(2,2)=iY(2,2)+c(i,1)*iA(i,ij(1))
			iY(1,2)=iY(1,2)-c(i,1)*iA(i,ij(2))
			iY(2,1)=iY(2,1)-c(i,2)*iA(i,ij(1))
		enddo
		iY(1,1)=iY(1,1)+1d0
		iY(2,2)=iY(2,2)+1d0
		pb=iY(1,1)*iY(2,2)-iY(1,2)*iY(2,1)
		iY=iY/pb
	end subroutine
	subroutine inv_update_row(c,ci,pb,A,iA)
		complex(8) :: c(:),A(:,:),iA(:,:),pb,ipb
		complex(8), allocatable :: tmp1(:),tmp2(:,:)
		integer :: ci,i,j,n
		ipb=1d0/pb
		n=size(iA,1)
		allocate(tmp1(n),tmp2(n,n))
		tmp1=0d0
		A(ci,:)=A(ci,:)+c
		do i=1,n
			do j=1,n
				tmp1(i)=tmp1(i)+c(j)*iA(j,i)
			enddo
		enddo
		do i=1,n
			do j=1,n
				tmp2(i,j)=iA(i,j)-ipb*iA(i,ci)*tmp1(j)
			enddo
		enddo
		iA=tmp2
	end subroutine
	subroutine inv_update_col(c,cj,pb,A,iA)
		complex(8) :: c(:),A(:,:),iA(:,:),pb,ipb
		complex(8), allocatable :: tmp1(:),tmp2(:,:)
		integer :: cj,i,j,n
		ipb=1d0/pb
		n=size(iA,1)
		allocate(tmp1(n),tmp2(n,n))
		tmp1=0d0
		A(:,cj)=A(:,cj)+c
		do i=1,n
			do j=1,n
				tmp1(i)=tmp1(i)+iA(i,j)*c(j)
			enddo
		enddo
		do i=1,n
			do j=1,n
				tmp2(i,j)=iA(i,j)-ipb*tmp1(i)*iA(cj,j)
			enddo
		enddo
		iA=tmp2
	end subroutine
	subroutine inv_update_rowcol(c,ij,iY,A,iA)
		complex(8) :: c(:,:),A(:,:),iA(:,:),iY(:,:)
		complex(8), allocatable :: tmp1(:,:),tmp2(:,:)
		integer :: ij(:),i,j,n
		n=size(iA,1)
		allocate(tmp1(n,2),tmp2(n,n))
		tmp1=0d0
		A(ij(1),:)=A(ij(1),:)+c(:,1)
		A(:,ij(2))=A(:,ij(2))+c(:,2)
		do i=1,n
			do j=1,n
				tmp1(i,1)=tmp1(i,1)+c(j,1)*iA(j,i)
				tmp1(i,2)=tmp1(i,2)+iA(i,j)*c(j,2)
			enddo
		enddo
		do i=1,n
			do j=1,n
				tmp2(i,j)=iA(i,j)-((iA(i,ij(1))*iY(1,1)+tmp1(i,2)*iY(2,1))*tmp1(j,1)+&
					(iA(i,ij(1))*iY(1,2)+tmp1(i,2)*iY(2,2))*iA(ij(2),j))
			enddo
		enddo
		iA=tmp2
	end subroutine
	subroutine inv_update_tworow(c,ij,iY,A,iA)
		complex(8) :: c(:,:),A(:,:),iA(:,:),iY(:,:)
		complex(8), allocatable :: tmp1(:,:),tmp2(:,:)
		integer :: ij(:),i,j,n
		n=size(iA,1)
		allocate(tmp1(n,2),tmp2(n,n))
		tmp1=0d0
		A(ij(1),:)=A(ij(1),:)+c(:,1)
		A(ij(2),:)=A(ij(2),:)+c(:,2)
		do i=1,n
			do j=1,n
				tmp1(i,1)=tmp1(i,1)+c(j,1)*iA(j,i)
				tmp1(i,2)=tmp1(i,2)+c(j,2)*iA(j,i)
			enddo
		enddo
		do i=1,n
			do j=1,n
				tmp2(i,j)=iA(i,j)-((iA(i,ij(1))*iY(1,1)+iA(i,ij(2))*iY(2,1))*tmp1(j,1)+&
					(iA(i,ij(1))*iY(1,2)+iA(i,ij(2))*iY(2,2))*tmp1(j,2))
			enddo
		enddo
		iA=tmp2
	end subroutine
	subroutine matrix_inv(A)
		complex(8) :: A(:,:)
		integer :: ipiv(size(A,1)),info
		call getrf(A,ipiv,info)
		if(info/=0) then
			write(*,*)"error1",info
			stop
		endif
		call getri(A,ipiv,info)
		if(info/=0) then
			write(*,*)"error2",info
			stop
		endif
	end subroutine
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
	!subroutine matrix_inv(a)
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
end module
