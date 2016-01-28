module det_utility
	use M_rd
	use M_matrix
	implicit none
contains
	!subroutine ratio(u,v,rt,D0,pn,n,k,m)
		!complex(8) :: u(:,:),v(:,:),rt(:),D0(:,:)
		!pn(n)=k
		!u(:,n)=D0(:,k)
		!v(:,n)=wf(m,:)
		!do i=0,n-1
			!u(:,n)=u(:,n)-rt(i)*sum(v(:,i)*u(:,n))*u(:,i)
		!enddo
		!rt(n)=1d0+sum(v(:,n)*u(:,n))
	!end subroutine
	subroutine get_pb(k,m,pb,WA,AW,WAW,iA,wf)
		complex(8) :: WA(:,:),pb
		complex(8), optional :: AW(:,:),WAW(:,:),iA(:,:),wf(:,:)
		!integer, allocatable :: k(:),m(:)
		!integer, allocatable :: k(:)
		integer :: k(:)
		integer :: m(:)
		integer :: i,j,n
		integer :: k_(size(k)/2+size(m)/2),m_(size(k_))
		complex(8) :: A(size(k_,1),size(k_,1))
		!if(.not.allocated(k)) then
			!pb=0d0
			!return
		!endif
		n=size(k)/2
		k_=(/k(1::2),m(1::2)/)
		m_=(/k(2::2),m(2::2)/)
		do i=1,size(k_)
			do j=1,size(k_)
				if(i<=n.and.j<=n) then
					A(i,j)=WA(m_(i),k_(j))
				elseif(i<=n.and.j>n) then
					A(i,j)=WAW(m_(i),m_(j))-wf(m_(i),m_(j))
				elseif(i>n.and.j<=n) then
					A(i,j)=iA(k_(i),k_(j))
				else
					A(i,j)=AW(k_(i),m_(j))
				endif
			enddo
		enddo
		select case(size(k_))
		case(0)
			pb=1d0
		case(1)
			pb=A(1,1)
		case(2)
			pb=(A(1,1)*A(2,2)-A(1,2)*A(2,1))
		case(3)
			pb=A(1,1)*A(2,2)*A(3,3)-A(1,1)*A(2,3)*A(3,2)-A(1,2)*A(2,1)*A(3,3)+A(1,2)*A(2,3)*A(3,1)+A(1,3)*A(2,1)*A(3,2)-A(1,3)*A(2,2)*A(3,1)
		case(4)
			pb=A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))-&
			   A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))+&
			   A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))-&
			   A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
	   case(5:)
		   write(*,*)"above 5 is not considered"
		   stop
	   end select
	end subroutine
	subroutine update(k,WA,AW,WAW,iA,W2,wf)
		complex(8) :: WA(:,:)
		complex(8), optional :: AW(:,:),WAW(:,:),iA(:,:),W2(:,:),wf(:,:)
		integer :: k(:)
		integer :: k_(size(k)/2),m_(size(k_))
		complex(8) :: WAr(size(k_),size(WA,2)),WAl(size(WA,1),size(k_)),iAl(size(WA,2),size(k_)),WAWr(size(k_),size(WA,1)),ipb(size(k_),size(m_)),dW2(size(k_),size(WA,1))
		integer :: i,j
		k_=k(1::2)
		m_=k(2::2)
		do i=1,size(k_)
			do j=1,size(k_)
				ipb(i,j)=WA(m_(i),k_(j))
			enddo
		enddo
		select case(size(k_))
		case(0)
			return
		case(1)
			ipb=1d0/ipb
			WAr=WA(m_,:)
			WAr(1,k_(1))=WAr(1,k_(1))-1d0
		case(2)
			ipb=reshape((/ipb(2,2),-ipb(2,1),-ipb(1,2),ipb(1,1)/),(/2,2/))/(ipb(1,1)*ipb(2,2)-ipb(2,1)*ipb(1,2))
			WAr=WA(m_,:)
			WAr(1,k_(1))=WAr(1,k_(1))-1d0
			WAr(2,k_(2))=WAr(2,k_(2))-1d0
		case default
			write(*,*)"err"
			stop
		end select
		WAr=matmul(ipb,WAr)
		WAl=WA(:,k_)
		WA=WA-matmul(WAl,WAr)
		if(present(iA)) then
			iAl=iA(:,k_)
			iA=iA-matmul(iAl,WAr)
		endif
		if(present(WAW)) then
			dW2=wf(m_,:)-W2(k_,:)
			WAWr=WAW(m_,:)-W2(k_,:)+matmul(WA(m_,k_)-diag(1d0,size(k_)),dW2)
			WAWr=matmul(ipb,WAWr)
			AW=AW+matmul(iA(:,k_),dW2)-matmul(iAl,WAWr)
			WAW=WAW+matmul(WA(:,k_),dW2)-matmul(WAl,WAWr)
			W2(k_,:)=wf(m_,:)
		endif
	end subroutine
	subroutine test()
		complex(8) :: wf(100,100),A(size(wf,1)/2,size(wf,2)),iA(size(A,1),size(A,2)/2),tmp(size(iA,1),size(iA,2))
		complex(8) :: WA(size(wf,1),size(iA,2)),AW(size(iA,1),size(wf,2)),WAW(size(wf,1),size(wf,2))
		integer :: Ns
		complex(8) :: pb,det1,det2
		call random_seed()
		call random_number(wf)
		Ns=size(A,1)
		A=wf(::2,:)
		iA=A(:,:Ns)
		call mat_inv(iA)
		WA=matmul(wf(:,:Ns),iA)
		AW=matmul(iA,A)
		WAW=matmul(wf(:,:Ns),AW)

		tmp=A(:,:Ns)
		det1=det(tmp)
		!call get_pb(shape(0),(/7,Ns+1/),pb,WA,iA=iA,AW=AW,WAW=WAW,wf=wf)
		!call get_pb((/5,Ns+4/),shape(0),pb,WA,iA=iA,AW=AW,WAW=WAW,wf=wf)
		!call get_pb((/5,Ns+4/),(/7,Ns+1,9,Ns+7/),pb,WA,iA=iA,AW=AW,WAW=WAW,wf=wf)
		call update((/5,Ns+4/),WA,iA=iA,AW=AW,WAW=WAW,W2=A,wf=wf)
		!A(5,:)=wf(Ns+4,:)
		!tmp=A(:,:Ns)
		!tmp(:,7)=A(:,Ns+1)
		!tmp(:,9)=A(:,Ns+7)
		!det2=det(tmp)

		tmp=A(:,:Ns)
		det1=det(tmp)
		call get_pb((/7,Ns+10,10,Ns+16/),shape(0),pb,WA,iA=iA,AW=AW,WAW=WAW,wf=wf)
		call update((/7,Ns+10,10,Ns+16/),WA,iA=iA,AW=AW,WAW=WAW,W2=A,wf=wf)
		tmp=A(:,:Ns)
		det2=det(tmp)

		write(*,*)sum(abs(matmul(A(:,:Ns),iA)-diag(1d0,Ns)))
		write(*,*)sum(abs(WA-matmul(wf(:,:Ns),iA)))
		write(*,*)sum(abs(AW-matmul(iA,A)))
		write(*,*)sum(abs(WAW-matmul(wf(:,:Ns),AW)))
		write(*,*)pb
		write(*,*)det2/det1
	end subroutine
end module

program main
	use det_utility
	implicit none
	call test()
end program
