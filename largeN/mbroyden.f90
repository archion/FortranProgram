module M_solution
	use M_matrix
	type t_id
		real(wp), allocatable :: dF(:,:),dV(:,:),w(:),beta(:,:)
		real(wp) :: alpha
		integer :: nsave
	end type
	interface mbroyden
		procedure :: rmbroyden,cmbroyden
	end interface
	type(t_id), private :: bid(5)
contains
	subroutine rmbroyden(n,vo,F,alpha,nsave,id)
		real(wp) :: vo(:),F(:)
		real(wp), optional :: alpha
		integer :: n
		integer, optional :: id,nsave
		integer :: i,j,k,m,md
		real(wp) :: v(size(vo,1))
		integer :: id_,s
		real(wp) :: norm,w0=0.01_wp,tmp
		if(present(id)) then
			id_=id
		else
			id_=1
		endif
		if(n==0) then
			if(.not.(present(alpha).and.present(nsave))) then
				stop "must provide alpha and nsave"
			endif
			bid(id_)%nsave=nsave
			bid(id_)%alpha=alpha
			allocate(bid(id_)%dF(size(vo),nsave),bid(id_)%dV(size(vo),nsave),bid(id_)%w(nsave),bid(id_)%beta(nsave,nsave))
			bid(id_)%dF=0._wp
			bid(id_)%dV=0._wp
		elseif(n==huge(1)) then
			if(allocated(bid(id_)%dF)) deallocate(bid(id_)%beta,bid(id_)%dF,bid(id_)%dV,bid(id_)%w)
		else
			s=size(vo)
			m=min(n-1,bid(id_)%nsave)
			md=bid(id_)%nsave

			v=vo
			if(n>1) then
				k=n-1-(n-2)/md*md
				bid(id_)%w(k)=1._wp
				bid(id_)%dV(:s,k)=vo-bid(id_)%dV(:s,k)
				bid(id_)%dF(:s,k)=F-bid(id_)%dF(:s,k)
				norm=sqrt(sum(bid(id_)%dF(:s,k)**2))
				if(norm<eps) then
					write(*,*)"not changing in r"
					norm=1._wp
				endif
				bid(id_)%dV(:s,k)=bid(id_)%dV(:s,k)/norm
				bid(id_)%dF(:s,k)=bid(id_)%dF(:s,k)/norm
			else
				k=0
			endif

			do i=1,m
				do j=1,m
					bid(id_)%beta(i,j)=bid(id_)%w(i)*bid(id_)%w(j)*sum(bid(id_)%dF(:s,i)*bid(id_)%dF(:s,j))
					if(i==j) bid(id_)%beta(i,j)=w0**2+bid(id_)%beta(i,j)
				enddo
			enddo
			if(m>0) call mat_inv(bid(id_)%beta(1:m,1:m))
			vo=v+bid(id_)%alpha*F
			do i=1,m
				tmp=0._wp
				do j=1,m
					tmp=tmp+bid(id_)%w(j)*sum(bid(id_)%dF(:s,j)*F)*bid(id_)%beta(j,i)
				enddo
				vo=vo-bid(id_)%w(i)*tmp*(bid(id_)%alpha*bid(id_)%dF(:s,i)+bid(id_)%dV(:s,i))
			enddo

			if((n-1)==0) then
				k=k+1
			else
				k=n-(n-1)/md*md
			endif
			bid(id_)%dV(:s,k)=v
			bid(id_)%dF(:s,k)=F
		endif
	end subroutine
	subroutine cmbroyden(n,vo,F,alpha,nsave,id)
		complex(wp) :: vo(:),F(:)
		real(wp), optional :: alpha
		integer :: n
		integer, optional :: id,nsave
		integer :: i,j,k,m,md,l
		real(wp) :: Fm(size(vo,1)*2),v(size(vo,1)*2)
		integer :: id_,s
		real(wp) :: norm,w0=0.01_wp,tmp
		if(present(id)) then
			id_=id
		else
			id_=1
		endif
		if(n==0) then
			if(.not.(present(alpha).and.present(nsave))) then
				stop "must provide alpha and nsave"
			endif
			bid(id_)%nsave=nsave
			bid(id_)%alpha=alpha
			allocate(bid(id_)%dF(size(vo)*2,nsave),bid(id_)%dV(size(vo)*2,nsave),bid(id_)%w(nsave),bid(id_)%beta(nsave,nsave))
		elseif(n==huge(1)) then
			if(allocated(bid(id_)%dF)) deallocate(bid(id_)%beta,bid(id_)%dF,bid(id_)%dV,bid(id_)%w)
		else
			s=size(vo)*2
			m=min(n-1,bid(id_)%nsave)
			md=bid(id_)%nsave
			
			Fm=[F%re,F%im]
			v=[vo%re,vo%im]
			if(n>1) then
				k=n-1-(n-2)/md*md
				bid(id_)%w(k)=1._wp
				bid(id_)%dV(:s,k)=[vo%re,vo%im]-bid(id_)%dV(:s,k)
				bid(id_)%dF(:s,k)=Fm-bid(id_)%dF(:s,k)
				norm=sqrt(sum(bid(id_)%dF(:s,k)**2))
				if(norm<eps) then
					write(*,*)"not changing"
					norm=1._wp
				endif
				bid(id_)%dV(:s,k)=bid(id_)%dV(:s,k)/norm
				bid(id_)%dF(:s,k)=bid(id_)%dF(:s,k)/norm
			else
				k=0
			endif

			do i=1,m
				do j=1,m
					bid(id_)%beta(i,j)=bid(id_)%w(i)*bid(id_)%w(j)*sum(bid(id_)%dF(:s,i)*bid(id_)%dF(:s,j))
					if(i==j) bid(id_)%beta(i,j)=w0**2+bid(id_)%beta(i,j)
				enddo
			enddo
			if(m>0) call mat_inv(bid(id_)%beta(1:m,1:m))
			vo%re=v(:s/2)+bid(id_)%alpha*Fm(:s/2)
			vo%im=v(s/2+1:)+bid(id_)%alpha*Fm(s/2+1:)
			do i=1,m
				tmp=0._wp
				do j=1,m
					tmp=tmp+bid(id_)%w(j)*sum(bid(id_)%dF(:s,j)*Fm)*bid(id_)%beta(j,i)
				enddo
				do l=1,s/2
					vo(l)%re=vo(l)%re-bid(id_)%w(i)*tmp*(bid(id_)%alpha*bid(id_)%dF(l,i)+bid(id_)%dV(l,i))
					vo(l)%im=vo(l)%im-bid(id_)%w(i)*tmp*(bid(id_)%alpha*bid(id_)%dF(size(v)/2+l,i)+bid(id_)%dV(size(v)/2+l,i))
				enddo
			enddo

			if((n-1)==0) then
				k=k+1
			else
				k=n-(n-1)/md*md
			endif
			bid(id_)%dV(:s,k)=v
			bid(id_)%dF(:s,k)=Fm
		endif
	end subroutine
end module
