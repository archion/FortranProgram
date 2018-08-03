subroutine mbroyden(n,vo,vn,alpha,nsave,id)
	real(8) :: vo(:),vn(:)
	real(8), optional :: alpha
	integer :: n
	integer, optional :: id,nsave
	integer :: i,j,k,m,md
	real(8) :: Fm(size(vo,1)),v(size(vo,1))
	real(8), allocatable, save :: dF(:,:,:),dV(:,:,:),w(:,:),beta(:,:)
	real(8), allocatable :: dF_(:,:,:),dV_(:,:,:),w_(:,:)
	real(8), save :: alpha_
	integer, save :: nsave_
	integer, save :: mid=1
	integer :: id_
	real(8) :: norm,w0=0.01d0,tmp
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
			w(k,id_)=1d0
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
			tmp=0d0
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
