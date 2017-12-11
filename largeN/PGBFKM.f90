module global
	use M_const
	use M_matrix
	implicit none
	!real(8) :: alpha=0.3d0,r=0.3d0,kp=0.5d0,fcut=1d0,fbeta=20d0,phicut=5d-2,phibeta=200d0,Jk=1.25d0,g=0d0
	real(8) :: alpha=0.3d0,r=0d0,kp=0.5d0,fcut=1d0,fbeta=20d0,phicut=5d-2,phibeta=200d0,Jk=2d0,g=0.5d0
	real(8) :: Tk=4d-3,beta,norm_Ac=1d0
	!integer, parameter :: nll=31*7,nadd1=5,nadd2=20,n=(n+nadd1*2+nadd2*2)*2
	integer, parameter :: n=536
	!real(8) :: omega((n+nadd1*2+nadd2*2)*2)
	real(8) :: omega(n)
contains
	subroutine set_omega(omega,base,nlin,width)
		real(8) :: omega(:),base,width
		integer :: nlin
		integer :: nlog,n,i,j
		real(8) :: u,l,t,tmp
		n=size(omega)/2
		nlog=n/nlin
		if(mod(n,nlin)/=0) stop "check n and nlin"
		u=width
		l=u/base
		do i=nlog,1,-1
			t=(u-l)/real(nlin,8)
			tmp=l
			do j=1,nlin
				tmp=tmp+t
				omega(n+(i-1)*nlin+j)=tmp
				omega(n-((i-1)*nlin+j)+1)=-tmp
			enddo
			u=l
			l=u/base
		enddo
	end subroutine
	subroutine add_omega(omega,center,width,nadd)
		real(8) :: omega(:),center,width
		integer :: nadd,n,m,l,u,i
		real(8) :: s,at1,at2
		n=size(omega)-nadd*4
		l=find_sidx(omega(1:n),center-width)
		u=find_sidx(omega(1:n),center+2d0*width)
		omega(u+1+nadd*4:)=omega(u+1:n)
		omega(n/2+1+nadd*2:l+nadd*2)=omega(n/2+1:l)
		omega(n-l+1+nadd*2:n/2+nadd*2)=omega(n-l+1:n/2)
		n=size(omega)
		l=l+nadd*2
		u=u+nadd*4
		m=u-l
		m=l+m-m/2
		s=width/100d0
		at1=atan(width)
		at2=atan(2d0*width)
		do i=l+1,u
			if(i<=m) then
				omega(i)=center-tan(s+(at1-s)/(m-l)*(m-i))
			else
				omega(i)=center+tan(s+(at2-s)/(u-m)*(i-m-1))
			endif
			omega(n-i+1)=-omega(i)
		enddo
	end subroutine
	real(8) function A_phi(x)
		real(8) :: x
		if(abs(x)<1.0d-12) then
			A_phi=0.0d0
		elseif(abs(x)<phicut) then
			A_phi=pi*abs(x)**(1d0-alpha)*sign(1d0,x)
		else
			A_phi=pi*phicut**(1d0-alpha)*sign(1d0,x)*1d0/((exp(phibeta*(x-phicut-4d0/phibeta))+1d0)*(exp(phibeta*(-x-phicut-4d0/phibeta))+1d0))
		endif
	end function
    real(8) function A_c(x)
		real(8) :: x
		if(abs(x)<fcut) then
			A_c=norm_Ac*exp(-fcut**2/pi)*(abs(x)/fcut)**r
		elseif(abs(x)<=(-omega(1))) then
			A_c=norm_Ac*exp(-x**2/pi)
		else
			A_c=0d0
		endif
	end function
	real(8) function integrate(A,l,u) result(weight)
		real(8) :: A
		integer :: i,il,iu
		integer, optional :: l,u
		il=1
		iu=size(omega)-1
		if(present(l)) then
			il=l
		endif
		if(present(u)) then
			iu=u
		endif
		weight=0d0
		do i=il,iu
			weight=weight+0.5d0*(omega(i+1)-omega(i))*(A(omega(i+1))+A(omega(i)))
		enddo
	end function
	integer function find_sidx(omega,x) result(l)
		real(8) :: omega(:),x
		integer :: m,u
		l=0
		u=size(omega)+1
		do
			if((u-l)>1) then
				m=(l+u)/2
				if(x>omega(m)) then
					l=m
				else
					u=m
				endif
			else
				exit
			endif
		enddo
	end function
	subroutine get_realpart(omega,im,re)
		real(8) :: omega(:),im(:),re(:)
		real(8) ::dx
		integer :: n,i,j,dn(2)
		n=size(omega)
		re(1)=0d0
		re(n)=0d0
		!re(2:n-1)=im(2:n-1)*log((omega(n)-omega(2:n-1))/(omega(2:n-1)-omega(1)))
		re(2:n-1)=2d0*im(2:n-1)*log((omega(n)-omega(2:n-1))/(omega(2:n-1)-omega(1)))
		do j=1,n
			dn=[min(n,j+1),max(1,j-1)]
			dx=omega(dn(1))-omega(dn(2))
			re(j)=re(j)+im(dn(1))-im(dn(2))
			do i=1,n
				if(j/=i) then
					re(i)=re(i)+(im(j)-im(i))/(omega(j)-omega(i))*dx
				endif
			enddo
		enddo
		!re=1d0/pi*re
		re=-0.5d0/pi*re
	end subroutine
	subroutine get_SEf1(omega,A,SE)
		real(8) :: omega(:),A(:),SE(:)
		integer :: i,j,m,n,iwe(size(omega))
		real(8) :: e1,e2,A1,A2,w,f0,x0,y0,D
		D=-omega(1)+1d-10
		n=size(omega)
		do i=1,n/2
			w=omega(i)
			x0=(A(n/2)+(A(n/2+1)-A(n/2))/(omega(n/2+1)-omega(n/2))*(0.0d0-omega(n/2)))*A_c(w)
			y0=w-1d0/beta*log(1d0-exp(beta*(w-D)))+1d0/beta*log(1d0-exp(-beta*(w+D)))
			f0=x0*y0
			do j=1,n
				iwe(j)=find_sidx(omega,w+omega(j))
			end do
			m=iwe(1)
			A1=merge(A(m)+(A(m+1)-A(m))/(omega(m+1)-omega(m))*(w+omega(1)-omega(m)),0d0,m>0.and.m<n)
			e1=omega(1)
			SE(i)=0d0
			do j=1,n-1
				do m=iwe(j),iwe(j+1)
					if(m==iwe(j+1)) then
						e2=omega(j+1)
						A2=merge(A(m)+(A(m+1)-A(m))/(omega(m+1)-omega(m))*(e2+w-omega(m)),0d0,m>0.and.m<n)
					else
						e2=omega(m+1)-w
						A2=merge(A(m+1),0d0,m+1>0.and.m+1<n)
					endif
					SE(i)=SE(i)+0.5d0*(e2-e1)*&
						!((fb(-e1-w)+ff(-e1))*(A_c(-e1)*A1-x0)+&
						!(fb(-e2-w)+ff(-e2))*(A_c(-e2)*A2-x0))
						((fb(e1+w)+ff(e1))*(A_c(-e1)*A1-x0)+&
						(fb(e2+w)+ff(e2))*(A_c(-e2)*A2-x0))
						!((fb(e1+w)+ff(e1))*(A_c(-e1)*A1)+&
						!(fb(e2+w)+ff(e2))*(A_c(-e2)*A2))
					A1=A2
					e1=e2
				enddo
			enddo
			SE(i)=SE(i)+f0
			SE(n+1-i)=SE(i)
		enddo
	end subroutine
	subroutine get_SEf2(omega,A,SE)
		real(8) :: omega(:),A(:),SE(:)
		integer :: i,j,m,n,iwe(size(omega))
		real(8) :: e1,e2,A1,A2,w,D,f0,x0,y0
		D=-omega(1)+1d-10
		n=size(omega)
		do i=1,n/2
			w=omega(i)
			x0=A(i)*A_phi(0.0d0)
			y0=w-1d0/beta*log(1d0+exp(beta*(w-D)))+1d0/beta*log(1d0+exp(-beta*(w+D)))
			f0=x0*y0
			do j=1,n
				iwe(j)=find_sidx(omega,w+omega(j))
			end do
			m=iwe(1)
			A1=merge(A(m)+(A(m+1)-A(m))/(omega(m+1)-omega(m))*(w+omega(1)-omega(m)),0d0,m>0.and.m<n)
			e1=omega(1)
			SE(i)=0d0
			do j=1,n-1
				do m=iwe(j),iwe(j+1)
					if(m==iwe(j+1)) then
						e2=omega(j+1)
						A2=merge(A(m)+(A(m+1)-A(m))/(omega(m+1)-omega(m))*(e2+w-omega(m)),0d0,m>0.and.m<n)
					else
						e2=omega(m+1)-w
						A2=merge(A(m+1),0d0,m+1>0.and.m+1<n)
					endif
					SE(i)=SE(i)+0.5d0*(e2-e1)*&
						((ff(-e1-w)+fb(-e1))*(A_phi(-e1)*A1-x0)+&
						(ff(-e2-w)+fb(-e2))*(A_phi(-e2)*A2-x0))
						!((ff(-e1-w)+fb(-e1))*(A_phi(-e1)*A1)+&
						!(ff(-e2-w)+fb(-e2))*(A_phi(-e2)*A2))
					A1=A2
					e1=e2
				enddo
			enddo
			SE(i)=SE(i)+f0
			SE(n+1-i)=SE(i)
		enddo
	end subroutine
	pure real(8) function ff(x)
		real(8), intent(in) :: x
		!if(x<omega(1)) then
			!ff=1d0
		!elseif(x>-omega(1)) then
			!ff=0d0
		!else
			ff=1d0/(exp(beta*x)+1d0)
		!endif
	end function
	pure real(8) function fb(x)
		real(8), intent(in) :: x
		if(abs(x)<1d-11) then
			fb=0d0
		else
			!if(x<omega(1)) then
				!fb=-1d0
			!elseif(x>-omega(1)) then
				!fb=0d0
			!else
				fb=1d0/(exp(beta*x)-1d0)
			!endif
		endif
	end function
	subroutine get_SEB(omega,A,SE)
		real(8) :: omega(:),A(:),SE(:)
		integer :: i,j,m,n,iwe(size(omega))
		real(8) :: e1,e2,A1,A2,w
		n=size(omega)
		!do i=1,n/2
		do i=1,n
			w=omega(i)
			do j=1,n
				iwe(j)=find_sidx(omega,w+omega(j))
			end do
			m=iwe(1)
			A1=merge(A(m)+(A(m+1)-A(m))/(omega(m+1)-omega(m))*(w+omega(1)-omega(m)),0d0,m>0.and.m<n)
			e1=omega(1)
			SE(i)=0d0
			do j=1,n-1
				do m=iwe(j),iwe(j+1)
					if(m==iwe(j+1)) then
						e2=omega(j+1)
						A2=merge(A(m)+(A(m+1)-A(m))/(omega(m+1)-omega(m))*(e2+w-omega(m)),0d0,m>0.and.m<n)
					else
						e2=omega(m+1)-w
						A2=merge(A(m+1),0d0,m+1>0.and.m+1<n)
					endif
					SE(i)=SE(i)+0.5d0*(e2-e1)*&
						!((ff(e1)-ff(e1+w))*A_c(e1)*A1+&
						!(ff(e2)-ff(e2+w))*A_c(e2)*A2)
						((ff(e1)-ff(e1+w))*A_c(-e1)*A1+&
						(ff(e2)-ff(e2+w))*A_c(-e2)*A2)
					A1=A2
					e1=e2
				enddo
			enddo
			!SE(n+1-i)=-SE(i)
		enddo
	end subroutine
	subroutine self_consistent(tol,niter,rate,Af,AB,iSEf,iSEB)
		real(8) :: tol,rate,Af(:),AB(:),iSEF(:),iSEB(:)
		integer :: niter
		real(8) :: SEf1(size(Af)),SEf2(size(Af)),iSEf_(size(Af)),rSEf(size(Af)),iSEB_(size(Af)),rSEB(size(Af)),Af_(size(Af)),AB_(size(Af))
		real(8) :: x(size(iSEf)*2),x_(size(iSEf)*2)
		integer :: i,n
		n=size(omega)
		beta=1d0/Tk
		Af=0d0
		AB=0d0
		!call mbroyden([real(8)::],[real(8)::],1d0,0,fin=.false.,id=2)
		do i=1,niter
			call get_realpart(omega,iSEf,rSEf)
			call get_realpart(omega,iSEB,rSEB)
			Af_=iSEf/((omega-rSEf)**2+iSEf**2)
			!if(i==2) then
				!do n=1,size(omega)
					!write(12,"(*(e28.20))")omega(n),iSEf(n),rSEf(n),Af_(n),Af(n)
				!enddo
				!stop
			!endif
			AB_=iSEB/((1d0/Jk-rSEB)**2+iSEB**2)
			if(all(abs(Af(n/2+1:)-Af_(n/2+1:))/(Af_(n/2+1:)+1d-70)<tol)) then
			!if(all(abs(Af(n/2+1:)-Af_(n/2+1:))<tol)) then
				!if(all(abs(AB(n/2+1:)-AB_(n/2+1:))/(AB_(n/2+1:)+1d-70)<tol)) then
					write(*,*)"converged"
					exit
				!endif
			endif
			if(i==niter) then
				write(*,*)"not converged",maxval(abs(Af(n/2+1:)-Af_(n/2+1:))/(Af_(n/2+1:)+1d-70))
			endif
			write(*,*)i,maxval(abs(Af(n/2+1:)-Af_(n/2+1:))/(Af_(n/2+1:)+1d-70))
			!write(*,*)i,maxval(abs(Af(n/2+1:)-Af_(n/2+1:)))
			Af=Af_
			AB=AB_

			call get_SEf1(omega,AB,SEf1)
			call get_SEf2(omega,Af,SEf2)
			call get_SEB(omega,Af,iSEB_)
			!iSEf_=pi*(kp*SEf1-g**2*SEf2)
			!iSEB_=pi*iSEB_
			iSEf_=-1d0/pi*(kp*SEf1-g**2*SEf2)
			iSEB_=-1d0/pi*iSEB_

			x=[iSEf,iSEB]
			x_=[iSEf_,iSEB_]
			!if(i<=8) then
				call mbroyden(x,x_,1d0,i,8)
			!else
				!call mbroyden(x,x_,1d0,i-8)
			!endif
			iSEf=x(1:size(x)/2)
			iSEB=x(size(x)/2+1:)
		enddo
		call mbroyden([real(8)::],[real(8)::],1d0,0,fin=.true.)
		do i=1,size(omega)
			write(12,"(*(e28.20))")omega(i),Af(i),Ab(i)
		enddo
		write(12,"(x/)")
	end subroutine
	subroutine mbroyden(vo,vn,alpha,n,n_max,fin,id)
		real(8) :: vo(:),vn(:),alpha
		integer :: n
		integer, optional :: n_max,id
		logical, optional :: fin
		integer :: i,j,k,n1,n2,m,md
		real(8) :: Fm(size(vo,1)),v(size(vo,1))
		real(8), allocatable, save :: dF(:,:,:),dV(:,:,:),w(:,:),beta(:,:)
		real(8), allocatable :: dF_(:,:,:),dV_(:,:,:),w_(:,:)
		integer, save :: mid=1
		integer :: id_
		real(8) :: norm,w0=0.01d0,tmp
		if(present(id)) then
			id_=id
		else
			id_=1
		endif
		if(present(fin)) then
			if(fin) then
				if(allocated(dF)) deallocate(beta,dF,dV,w)
				return
			else
				mid=id_
				return
			endif
		endif
		n1=size(vo)
		if(present(n_max)) then
			m=min(n-1,n_max)
			md=n_max
			if(n==1.and..not.allocated(dF)) allocate(dF(n1,n_max,mid),dV(n1,n_max,mid),w(n_max,mid))
		else
			m=n-1
			md=n
			if(n==1) then
				if(.not.allocated(dF)) allocate(dF(n1,8,mid),dV(n1,8,mid),w(8,mid))
			else
				n2=size(dF,2)
				if(n2<n) then
					allocate(dF_(n1,n2,mid),dV_(n1,n2,mid),w_(n2,mid))
					dF_=dF
					dV_=dV
					w_=w
					deallocate(dF,dV,w)
					n2=n2+max(8,nint(n2*0.2d0))
					allocate(dF(n1,n2,mid),dV(n1,n2,mid),w(n2,mid))
					dF(:,1:size(dF_,2),:)=dF_
					dV(:,1:size(dV_,2),:)=dV_
					w(1:size(w_,1),:)=w_
					deallocate(dF_,dV_,w_)
				endif
			endif
		endif
		if(allocated(beta)) then
			if(m/=size(beta,1)) then
				deallocate(beta)
				allocate(beta(m,m))
			endif
		else
			allocate(beta(m,m))
		endif

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
		if(m>0) call mat_inv(beta)
		vo=v+alpha*Fm
		do i=1,m
			tmp=0d0
			do j=1,m
				tmp=tmp+w(j,id_)*sum(dF(:,j,id_)*Fm)*beta(j,i)
			enddo
			vo=vo-w(i,id_)*tmp*(alpha*dF(:,i,id_)+dV(:,i,id_))
		enddo

		if((n-1)==0) then
			k=k+1
		else
			k=n-(n-1)/md*md
		endif
		dV(:,k,id_)=v
		dF(:,k,id_)=Fm
	end subroutine
end module
program main
	use global
	implicit none
	real(8) :: Af(n),AB(n),rGf(n),rGB(n),iSEF(n),iSEB(n),dTk
	integer :: i
	open(10,file="init.dat")
	open(12,file="/tmp/a_2.dat")
	open(13,file="TEMPLIST.dat")
	!call set_omega(omega(1:nll*2),2.3d0,7,10d0)
	!call add_omega(omega(1:nll*2+nadd1*4),1d0,0.1d0,nadd1)
	!call add_omega(omega,0.8d0,0.16d0,nadd2)
	do i=1,size(omega)
		read(10,"(5e28.20)")omega(i),iSEf(i),iSEB(i)
	enddo

	Tk=0.09d0
	!Tk=0.01d0

	fcut=abs(omega(find_sidx(omega,-fcut)))
	Jk=Jk*integrate(A_c)/pi
	norm_Ac=norm_Ac/integrate(A_c)*pi

	!Tkn=[1d-6,1d-7,1d-4,1d-5,1d-6,1d-7,1d-8,1d-9]
	!Tk=1d-5
	do i=1,100
		read(13,*)Tk
		write(*,*)Tk
		call self_consistent(1d-6,1000,3d-3,Af,AB,iSEf,iSEB)
		!stop
	enddo
end program
