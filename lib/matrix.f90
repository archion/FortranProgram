module M_matrix
	use lapack95, only : getrf, getri, heevd, heev, heevx, heevr, sysv
	use ifport, ifort_qsort => qsort
	use M_const
	implicit none
	interface det
		module procedure rdet,cdet
	end interface
	interface diag
		module procedure mdiag, ndiag, mcdiag, ncdiag, mmdiag, mmcdiag
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
	subroutine rmat_inv(A,info)
		real(wp) :: A(:,:)
		real(wp) :: det
		real(min(dp,wp)) :: A_(size(A,1),size(A,2))
		integer :: ipiv(size(A,1)),info1
		integer, optional :: info
		if(size(A,1)==1.and.size(A,2)==1) then
			A=1._wp/A
			return
		endif
		if(size(A,1)==2.and.size(A,2)==2) then
			det=(A(1,1)*A(2,2)-A(1,2)*A(2,1))
			if(abs(det)<eps*100._wp) then
				write(*,*)"inverse matrix err det=0"
				write(*,*)RAISEQQ(SIG$ABORT)
			endif
			det=1._wp/det
			A=det*reshape([A(2,2),-A(2,1),-A(1,2),A(1,1)],[2,2])
			return
		endif
		if(present(info)) then
			info=0
		endif
		A_=real(A,kind=min(dp,wp))
		call getrf(A_,ipiv,info1)
		if(info1/=0) then
			if(present(info)) then
				info=info1
				return
			else
				write(*,*)"inverse matrix err1"
				write(*,*)RAISEQQ(SIG$ABORT)
			endif
		endif
		call getri(A_,ipiv,info1)
		A=real(A_,kind=wp)
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
		complex(wp) :: A(:,:)
		complex(wp) :: det
		complex(min(dp,wp)) :: A_(size(A,1),size(A,2))
		integer :: ipiv(size(A,1)),info1
		integer, optional :: info
		if(size(A,1)==1.and.size(A,2)==1) then
			A=1._wp/A
			return
		endif
		if(size(A,1)==2.and.size(A,2)==2) then
			det=(A(1,1)*A(2,2)-A(1,2)*A(2,1))
			if(abs(det)<eps*100._wp) then
				write(*,*)"inverse matrix err det=0"
				write(*,*)RAISEQQ(SIG$ABORT)
			endif
			det=1._wp/det
			A=det*reshape([A(2,2),-A(2,1),-A(1,2),A(1,1)],[2,2])
			return
		endif
		if(present(info)) then
			info=0
		endif
		A_=cmplx(A,kind=min(dp,wp))
		call getrf(A_,ipiv,info1)
		if(info1/=0) then
			if(present(info)) then
				info=info1
				return
			else
				write(*,*)info1
				write(*,*)"inverse matrix err1"
				write(*,*)RAISEQQ(SIG$ABORT)
			endif
		endif
		call getri(A_,ipiv,info1)
		A=cmplx(A_,kind=min(wp,dp))
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
	function cdet(A,info) result(rt)
		complex(wp) :: A(:,:)
		complex(min(wp,dp)) :: A_(size(A,1),size(A,2))
		integer, optional :: info
		complex(wp) :: rt
		integer :: i,ipiv(size(A,1))
		A_=cmplx(A,kind=min(wp,dp))
		if(present(info)) then
			call getrf(A_,ipiv,info)
			if(info/=0) then
				stop "determinant err"
			endif
		else
			call getrf(A_,ipiv)
		endif
		rt=1._wp
		do i=1,size(A,1)
			if(ipiv(i)/=i) then
				rt=-rt*A_(i,i)
			else
				rt=rt*A_(i,i)
			endif
		enddo
	end function
	function rdet(A,info) result(rt)
		real(wp) :: A(:,:)
		real(min(wp,dp)) :: A_(size(A,1),size(A,2))
		integer, optional :: info
		real(wp) :: rt
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
		rt=1._wp
		do i=1,size(A,1)
			if(ipiv(i)/=i) then
				rt=-rt*A_(i,i)
			else
				rt=rt*A_(i,i)
			endif
		enddo
	end function
	function mmdiag(a)
		real(wp) :: a(:,:) 
		real(wp) :: mmdiag(size(a,1))
		integer :: i
		do i=1,size(a,1)
			mmdiag(i)=a(i,i)
		enddo
	end function
	function mmcdiag(a)
		complex(wp) :: a(:,:) 
		complex(wp) :: mmcdiag(size(a,1))
		integer :: i
		do i=1,size(a,1)
			mmcdiag(i)=a(i,i)
		enddo
	end function
	function mdiag(a)
		real(wp) :: a(:) 
		real(wp) :: mdiag(size(a),size(a))
		integer :: i
		mdiag=0._wp
		do i=1,size(a)
			mdiag(i,i)=a(i)
		enddo
	end function
	function ndiag(a,n)
		real(wp) :: a
		integer :: n
		real(wp) :: ndiag(n,n)
		integer :: i
		ndiag=0._wp
		do i=1,n
			ndiag(i,i)=a
		enddo
	end function
	function mcdiag(a)
		complex(wp) :: a(:) 
		complex(wp) :: mcdiag(size(a),size(a))
		integer :: i
		mcdiag=0._wp
		do i=1,size(a)
			mcdiag(i,i)=a(i)
		enddo
	end function
	function ncdiag(a,n)
		complex(wp) :: a
		integer :: n
		complex(wp) :: ncdiag(n,n)
		integer :: i
		ncdiag=0._wp
		do i=1,n
			ncdiag(i,i)=a
		enddo
	end function
	function Tr(A,B)
		complex(wp) :: A(:,:),B(:,:),Tr
		integer :: n,i,j
		n=size(A,1)
		Tr=0._wp
		do i=1,n
			do j=1,n
				Tr=Tr+A(i,j)*B(j,i)
			enddo
		enddo
	end function
	function check_diag(A,U,E,err)
		complex(wp) :: A(:,:),U(:,:)
		real(wp) :: E(:)
		real(wp), optional :: err
		real(wp) :: err_
		logical :: check_diag
		check_diag=.true.
		if(present(err)) then
			err_=err
		else
			err_=1e-6_wp
		endif
		if(any(abs((matmul(transpose(conjg(U)),matmul(A,U))-diag(E)))>err_)) then
			check_diag=.false.
		endif
	end function
	subroutine r_conjgrad(A,b,x)
		real(wp) :: A(:,:),b(:),x(:),r(size(x)),p(size(x)),Ap(size(x)),al,tmp0,tmp1,cvg=1e-10_wp
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
		real(wp) :: va(:),y(:),x(:)
		do i=1,size(x)
			y(i)=0._wp
			do j=ia(i),ia(i+1)-1
				y(i)=y(i)+va(j)*x(ja(j))
			enddo
		enddo
	end subroutine
	subroutine c_crsmv(va,ja,ia,x,y)
		integer :: i,j
		integer :: ja(:),ia(:)
		complex(wp) :: va(:),y(:),x(:)
		do i=1,size(x)
			y(i)=0._wp
			do j=ia(i),ia(i+1)-1
				y(i)=y(i)+va(j)*x(ja(j))
			enddo
		enddo
	end subroutine
	subroutine crs(v,va,ja,ia)
		integer :: n
		complex(wp) :: va(:),v(:,:)
		integer :: ja(:),ia(size(v,1)+1)
		integer ::  i,j
		ia(1)=1
		n=1
		do i=1,size(v,1)
			do j=1,size(v,2)
				if((real(v(i,j))**2+imag(v(i,j))**2)>1e-18_wp) then
					va(n)=v(i,j)
					ja(n)=j
					n=n+1
				endif
			enddo
			ia(i+1)=n
		enddo
	end subroutine
	subroutine diag2(A,E,info)
		complex(wp) :: A(:,:)
		real(wp) :: E(:)
		integer, optional :: info
		real(wp) :: tmp
		if(abs(A(2,1))<1e-7_wp) then
			E=real((/A(1,1),A(2,2)/))
			if(E(1)<=E(2)) then
				A(1,1)=1._wp
				A(2,2)=1._wp
			else
				E=real((/A(2,2),A(1,1)/))
				A(1,2)=1._wp
				A(2,1)=1._wp
				A(1,1)=0._wp
				A(2,2)=0._wp
			endif
		else
			tmp=sqrt((A(1,1)-A(2,2))**2+4._wp*A(1,2)*A(2,1))
			E=0.5_wp*(A(1,1)+A(2,2)+(/-tmp,tmp/))
			A(1,:)=E-A(2,2)
			A(2,:)=A(2,1)
			A(:,1)=A(:,1)/sqrt(A(1,1)**2+A(2,1)*conjg(A(2,1)))
			A(:,2)=A(:,2)/sqrt(A(1,2)**2+A(2,2)*conjg(A(2,2)))
		endif
	end subroutine
	subroutine diag4(A,E,info)
		complex(wp) :: A(:,:)
		real(wp) :: E(:)
		integer, optional :: info
		complex(wp) :: tmp2(2,2),tmp4(4,4)
		tmp4=0._wp
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
		A(1:2,3:4)=0._wp
		A(3:4,1:2)=0._wp
		tmp4(2,1)=0._wp
		tmp4(1,2)=0._wp
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
		real(wp) :: A(:),B(:)
		real(wp) :: rt(size(A)*size(B))
		integer :: i,n
		n=size(A)
		do i=1,size(B)
			rt(n*(i-1)+1:n*i)=A*B(i)
		enddo
	end function
	function outprod_c(A,B) result(rt)
		complex(wp) :: A(:),B(:)
		complex(wp) :: rt(size(A)*size(B))
		integer :: i,n
		n=size(A)
		do i=1,size(B)
			rt(n*(i-1)+1:n*i)=A*B(i)
		enddo
	end function
	subroutine mat_diag(H,E,check)
		complex(wp) :: H(:,:)
		real(wp) :: E(:)
		complex(min(dp,wp)) :: H_(size(H,1),size(H,2))
		logical, optional :: check
		if(present(check).or.(wp>dp)) then
			H_=cmplx(H,kind=min(dp,wp))
		endif
		select case(size(E))
		case(2)
			call diag2(H,E)
		!case(4)
			!call diag4(H,E)
		case(100:)
			if(wp>dp) then
				call heevd(H_,E,"V")
				H=cmplx(H_,kind=wp)
			else
				call heevd(H,E,"V")
			endif
		case default
			if(wp>dp) then
				call heev(H_,E,"V")
				H=cmplx(H_,kind=wp)
			else
				call heev(H,E,"V")
			endif
		end select
		if(present(check)) then
			if(.not.check_diag(H_,H,E)) then
				write(*,*)"diag err"
			endif
		endif
	end subroutine
end module
