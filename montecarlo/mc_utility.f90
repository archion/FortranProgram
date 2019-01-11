module mc_utility
	use blas95, only : gemm, gemv, her, gerc
	use M_const
	implicit none
contains
	subroutine get_pb(k,m,pb,WA,WAl,WAr,iA,iAl,AW,WAW,AWr)
		integer :: k(:),m(:)
		complex(wp) :: WA(:,:),pb
		complex(wp), optional :: WAl(:,:),WAr(:,:),iA(:,:),iAl(:,:),AW(:,:),WAW(:,:),AWr(:,:)
		integer :: i,j,n
		integer :: k_(size(k)/2+size(m)/2),m_(size(k_))
		complex(wp) :: A(size(k_,1),size(k_,1))
		n=size(k)/2
		k_=[k(1::2),m(1::2)]
		m_=[k(2::2),m(2::2)]
		do i=1,size(m_)
			do j=1,size(k_)
				if(i<=n.and.j<=n) then
					if(present(WAl)) then
						A(i,j)=WA(m_(i),k_(j))+sum(WAl(m_(i),:)*WAr(:,k_(j)))
					else
						A(i,j)=WA(m_(i),k_(j))
					endif
				elseif(i<=n.and.j>n) then
					if(present(AWr)) then
						A(i,j)=WAW(m_(i),m_(j))+sum(WAl(m_(i),:)*AWr(:,m_(j)))
					else
						A(i,j)=WAW(m_(i),m_(j))
					endif
				elseif(i>n.and.j<=n) then
					if(present(iAl)) then
						A(i,j)=iA(k_(i),k_(j))+sum(iAl(k_(i),:)*WAr(:,k_(j)))
					else
						A(i,j)=iA(k_(i),k_(j))
					endif
				else
					if(present(AWr)) then
						A(i,j)=AW(k_(i),m_(j))+sum(iAl(k_(i),:)*AWr(:,m_(j)))
					else
						A(i,j)=AW(k_(i),m_(j))
					endif
				endif
			enddo
		enddo
		!A(1:n,1:n)=WA(m_(1:n),k_(1:n))
		!A(1:n,n+1:)=WAW(m_(1:n),m_(n+1:))
		!A(n+1:,1:n)=iA(k_(n+1:),k_(1:n))
		!A(n+1:,n+1:)=AW(k_(n+1:),m_(n+1:))
		select case(size(k_))
		case(0)
			pb=1._wp
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
	subroutine update(k,WA,WAl,WAr,iA,iAl,AW,WAW,AWr,cwf,wf,full)
		integer :: k(:)
		complex(wp) :: WA(:,:),WAl(:,:),WAr(:,:)
		complex(wp), optional :: AW(:,:),WAW(:,:),iA(:,:),cwf(:,:),wf(:,:),iAl(:,:),AWr(:,:)
		logical :: full
		integer :: k_(size(k)/2),m_(size(k_))
		integer :: n,i,j,l,ct=0
		complex(wp) :: A(size(k_),size(k_)),beta=1._wp
		if(size(k)>0) then
			n=size(WAl,2)-size(k_)
			k_=[k(1::2)]
			m_=[k(2::2)]
			do i=1,size(m_)
				do j=1,size(k_)
					A(i,j)=WA(m_(i),k_(j))+sum(WAl(m_(i),:n)*WAr(:n,k_(j)))
				enddo
			enddo
			select case(size(k_))
			case(0)
				return
			case(1)
				A=-1._wp/A
			case(2)
				A=-1._wp/(A(1,1)*A(2,2)-A(2,1)*A(1,2))*reshape([A(2,2),-A(2,1),-A(1,2),A(1,1)],[2,2])
			case default
				write(*,*)"err"
				stop
			end select
			WAl(:,n+1:)=WA(:,k_)
			WAr(n+1:,:)=WA(m_,:)
			do i=1,size(k_)
				call gemv(WAl(:,:n),WAr(:n,k_(i)),WAl(:,n+i),beta=beta)
				call gemv(WAr(:n,:),WAl(m_(i),:n),WAr(n+i,:),beta=beta,trans='T')
				WAr(n+i,k_(i))=WAr(n+i,k_(i))-1._wp
			enddo
			WAr(n+1:,:)=matmul(A,WAr(n+1:,:))
			if(present(iA)) then
				do i=1,size(k_)
					iAl(:,n+i)=iA(:,k_(i))
					call gemv(iAl(:,:n),WAr(:n,k_(i)),iAl(:,n+i),beta=beta)
				enddo
				if(present(WAW)) then
					cwf(k_,:)=wf(m_,:)-cwf(k_,:)
					AWr(n+1:,:)=WAW(m_,:)+matmul(WAl(m_,n+1:),cwf(k_,:))
					do i=1,size(k_)
						call gemv(AWr(:n,:),WAl(m_(i),:n),AWr(n+i,:),beta=beta,trans='T')
					enddo
					AWr(n+1:,:)=cwf(k_,:)+matmul(A,AWr(n+1:,:))
					cwf(k_,:)=wf(m_,:)
				endif
			endif
		endif
		if(full.and.size(WAl,2)>0) then
			if(size(WAl,2)>2) then
			!if(.true.) then
				call gemm(WAl,WAr,WA,beta=beta)
				if(present(iA)) then
					call gemm(iAl,WAr,iA,beta=beta)
					if(present(WAW)) then
						call gemm(iAl,AWr,AW,beta=beta)
						call gemm(WAl,AWr,WAW,beta=beta)
					endif
				endif
			else
				!$omp parallel
				!$omp do
				do i=1,size(WA,2)
					do j=1,size(WAl,2)
						do l=1,size(WA,1)
							WA(l,i)=WA(l,i)+WAl(l,j)*WAr(j,i)
						enddo
					enddo
				enddo
				!$omp end do
				if(present(iA)) then
					!$omp do
					do i=1,size(iA,2)
						do j=1,size(WAl,2)
							do l=1,size(iA,1)
								iA(l,i)=iA(l,i)+iAl(l,j)*WAr(j,i)
							enddo
						enddo
					enddo
					!$omp end do
					if(present(WAW)) then
						!$omp do
						do i=1,size(AW,2)
							do j=1,size(WAl,2)
								do l=1,size(AW,1)
									AW(l,i)=AW(l,i)+iAl(l,j)*AWr(j,i)
								enddo
								do l=1,size(WAW,1)
									WAW(l,i)=WAW(l,i)+WAl(l,j)*AWr(j,i)
								enddo
							enddo
						enddo
						!$omp end do
					endif
				endif
				!$omp end parallel
			endif
			!ct=ct+1
			!if(mod(ct,1000)==0) write(*,*)ct
		endif
	end subroutine
end module
