module M_solution
	use M_matrix
	interface mbroyden
		procedure :: rmbroyden,cmbroyden
	end interface
contains
	subroutine rmbroyden(n,vo,vn,alpha,nsave,id)
		real(wp) :: vo(:),vn(:)
		real(wp), optional :: alpha
		integer :: n
		integer, optional :: id,nsave
		integer :: i,j,k,m,md
		real(wp) :: Fm(size(vo,1)),v(size(vo,1))
		real(wp), allocatable, save :: dF(:,:,:),dV(:,:,:),w(:,:),beta(:,:)
		real(wp), allocatable :: dF_(:,:,:),dV_(:,:,:),w_(:,:)
		real(wp), save :: alpha_
		integer, save :: nsave_
		integer, save :: mid=1
		integer :: id_
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
			nsave_=nsave
			alpha_=alpha
			mid=id_
			allocate(dF(size(vo),nsave,mid),dV(size(vo),nsave,mid),w(nsave,mid),beta(nsave,nsave))
		elseif(n==huge(1)) then
			if(allocated(dF)) deallocate(beta,dF,dV,w)
		else
			m=min(n-1,nsave_)
			md=nsave_

			Fm=vn-vo
			v=vo
			if(n>1) then
				k=n-1-(n-2)/md*md
				w(k,id_)=1._wp
				dV(:,k,id_)=vo-dV(:,k,id_)
				dF(:,k,id_)=Fm-dF(:,k,id_)
				norm=sqrt(sum(dF(:,k,id_)**2))
				dV(:,k,id_)=dV(:,k,id_)/norm
				dF(:,k,id_)=dF(:,k,id_)/norm
			else
				k=0
			endif

			do i=1,m
				do j=1,m
					beta(i,j)=w(i,id_)*w(j,id_)*sum(dF(:,i,id_)*dF(:,j,id_))
					if(i==j) beta(i,j)=w0**2+beta(i,j)
				enddo
			enddo
			if(m>0) call mat_inv(beta(1:m,1:m))
			vo=v+alpha_*Fm
			do i=1,m
				tmp=0._wp
				do j=1,m
					tmp=tmp+w(j,id_)*sum(dF(:,j,id_)*Fm)*beta(j,i)
				enddo
				vo=vo-w(i,id_)*tmp*(alpha_*dF(:,i,id_)+dV(:,i,id_))
			enddo

			if((n-1)==0) then
				k=k+1
			else
				k=n-(n-1)/md*md
			endif
			dV(:,k,id_)=v
			dF(:,k,id_)=Fm
		endif
	end subroutine
	subroutine cmbroyden(n,vo,vn,alpha,nsave,id)
		complex(wp) :: vo(:),vn(:)
		real(wp), optional :: alpha
		integer :: n
		integer, optional :: id,nsave
		integer :: i,j,k,m,md,l
		real(wp) :: Fm(size(vo,1)*2),v(size(vo,1)*2)
		real(wp), allocatable, save :: dF(:,:,:),dV(:,:,:),w(:,:),beta(:,:)
		real(wp), allocatable :: dF_(:,:,:),dV_(:,:,:),w_(:,:)
		real(wp), save :: alpha_
		integer, save :: nsave_
		integer, save :: mid=1
		integer :: id_
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
			nsave_=nsave
			alpha_=alpha
			mid=id_
			allocate(dF(size(vo)*2,nsave,mid),dV(size(vo)*2,nsave,mid),w(nsave,mid),beta(nsave,nsave))
		elseif(n==huge(1)) then
			if(allocated(dF)) deallocate(beta,dF,dV,w)
		else
			m=min(n-1,nsave_)
			md=nsave_

			Fm=[vn%re-vo%re,vn%im-vo%im]
			v=[vo%re,vo%im]
			if(n>1) then
				k=n-1-(n-2)/md*md
				w(k,id_)=1._wp
				dV(:,k,id_)=[vo%re,vo%im]-dV(:,k,id_)
				dF(:,k,id_)=Fm-dF(:,k,id_)
				norm=sqrt(sum(dF(:,k,id_)**2))
				dV(:,k,id_)=dV(:,k,id_)/norm
				dF(:,k,id_)=dF(:,k,id_)/norm
			else
				k=0
			endif

			do i=1,m
				do j=1,m
					beta(i,j)=w(i,id_)*w(j,id_)*sum(dF(:,i,id_)*dF(:,j,id_))
					if(i==j) beta(i,j)=w0**2+beta(i,j)
				enddo
			enddo
			if(m>0) call mat_inv(beta(1:m,1:m))
			vo%re=v(:size(v)/2)+alpha_*Fm(:size(v)/2)
			vo%im=v(size(v)/2+1:)+alpha_*Fm(size(v)/2+1:)
			do i=1,m
				tmp=0._wp
				do j=1,m
					tmp=tmp+w(j,id_)*sum(dF(:,j,id_)*Fm)*beta(j,i)
				enddo
				do l=1,size(vo)
					vo(l)%re=vo(l)%re-w(i,id_)*tmp*(alpha_*dF(l,i,id_)+dV(l,i,id_))
					vo(l)%im=vo(l)%im-w(i,id_)*tmp*(alpha_*dF(size(v)/2+l,i,id_)+dV(size(v)/2+l,i,id_))
				enddo
			enddo

			if((n-1)==0) then
				k=k+1
			else
				k=n-(n-1)/md*md
			endif
			dV(:,k,id_)=v
			dF(:,k,id_)=Fm
		endif
	end subroutine
end module
