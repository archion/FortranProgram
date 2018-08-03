module M_DMFT
	use M_const
	use M_matrix
	use mkl_service
	implicit none
	real(8), parameter :: onepi=1d0/pi,halfpi=onepi*0.5d0
	integer, parameter :: nlog=1,nl=800,nadd(1)=[0],n=(nl*nlog+sum(nadd*2))*2,nk=256,D=2
	real(8) :: U=9d0,t=1d0,bcut=10d0,base=10d0/1d-3,iomega(n),omega(n),Tk
	type gfs(n,D,nk)
		integer, len :: n,D,nk
		real(8) :: k(D,(2*nk)**D)
		real(8) :: rhoc(n)
		complex(8) :: Gc(n),Delta(n),SEc(n)
		logical :: conv
	contains
		procedure :: self_consistent
		procedure :: impurity_solver
		procedure :: export
		procedure, nopass :: ek
	end type
contains
	include "../largeN/mbroyden.f90"
	real(8) function ek(k) result(rt)
		real(8) :: k(D)
		integer :: i
		rt=0d0
		do i=1,size(k)
			rt=rt+2d0*t*cos(k(i))
		enddo
	end function
	subroutine export(self,ut)
		class(gfs(*,*,*)) :: self
		integer :: ut
		integer :: i,k
		write(ut,"(A)")"Tk omega rGc iGc rDelta iDelta rSEc iSEc"
		do i=1,n
			write(ut,"(*(e28.20))")Tk,omega(i),self%Gc(i),self%Delta(i),self%SEc(i)
		enddo
		write(ut,"(x/)")
	end subroutine
	subroutine impurity_solver(self)
		class(gfs(*,*,*)) :: self
		self%Gc=0.5d0/(omega-self%Delta+U/2d0)+0.5d0/(omega-self%Delta-U/2d0)
	end subroutine
	subroutine self_consistent(self,rate,niter,tol)
		class(gfs(*,*,*)) :: self
		real(8) :: rate,tol
		integer :: niter
		integer :: i,j,ik
		real(8) :: x(self%n),x_(size(x))
		call mbroyden(0,x,x_,rate,40)
		do i=1,niter
			call set_realpart(omega,self%Delta)
			call self%impurity_solver()
			x=self%Delta%im
			self%SEc=omega-self%Delta-1d0/self%Gc
			if(D==0) then
			else
				!$OMP PARALLEL DO 
				do j=1,n
					self%Gc(j)=0d0
					do ik=1,(2*nk)**D
						self%Gc(j)=self%Gc(j)+(1d0/(omega(j)-ek(self%k(:,ik))-self%SEc(j)))
					enddo
				enddo
				!$OMP END PARALLEL DO
			endif
			self%Gc=self%Gc/real((2*nk)**D,kind=8)
			x_=imag(omega-1d0/self%Gc-self%SEc)
			write(*,*)i,maxval(abs(x-x_))
			if(all(abs(x-x_)<tol).or.i==niter+1) then
				self%conv=(i/=niter+1)
				exit
			else
				call mbroyden(i,x,x_)
			endif
			self%Delta%im=x
		enddo
		call mbroyden(huge(1),x,x_)
	end subroutine
	subroutine set_realpart(omega,G)
		real(8) :: omega(:)
		complex(8) :: G(:)
		real(8) ::dx
		integer :: i,j,dn(2)
		G(1)%re=0d0
		G(n)%re=0d0
		do i=2,n-1
			G(i)%re=2d0*G(i)%im*log((omega(n)-omega(i))/(omega(i)-omega(1)))
		enddo
		do j=1,n
			dn=[min(n,j+1),max(1,j-1)]
			dx=omega(dn(1))-omega(dn(2))
			G(j)%re=G(j)%re+G(dn(1))%im-G(dn(2))%im
			do i=1,n
				if(j/=i) then
					G(i)%re=G(i)%re+(G(j)%im-G(i)%im)/(omega(j)-omega(i))*dx
				endif
			enddo
		enddo
		do i=1,n
			G(i)%re=G(i)%re*halfpi
		enddo
	end subroutine
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
		at2=atan(width)
		do i=l+1,u
			if(i<=m) then
				omega(i)=center-tan(s+(at1-s)/(m-l)*(m-i))
			else
				omega(i)=center+tan(s+(at2-s)/(u-m)*(i-m-1))
			endif
			omega(n-i+1)=-omega(i)
		enddo
	end subroutine
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
end module
include "../lib/serde.f90"
program main
	use M_DMFT
	use M_serde
	implicit none
	type(gfs(n,D,nk)) :: phy
	integer :: j,i,ik
	type(t_serde(D)) :: map

	call omp_set_nested(.false.)
	call omp_set_dynamic(.false.)
	call mkl_set_dynamic(0)

	call mkl_set_num_threads(24)
	call omp_set_num_threads(mkl_get_max_threads())

	open(unit=10,File="../data/fk.dat")
	open(unit=20,file="../data/fk_init.dat")
	call set_omega(omega(1:nl*nlog*2),base,nl,bcut)
	!call add_omega(omega(1:nl*nlog*2+sum(nadd(1:1))*4),4d0,0.5d0,nadd(1))
	!call add_omega(omega(1:nl*nlog*2+sum(nadd(1:2))*4),4d0,0.5d0,nadd(2))
	map%shap=[(2*nk,i=1,D)]
	do ik=1,size(phy%k,2)
		phy%k(:,ik)=pi*map%get_idx(ik)/nk-pi
	enddo
	!do i=1,n
		!read(20,"(2(e28.20))")phy%Delta(i)
	!enddo
	phy%Delta%im=-exp(abs(omega))
	call phy%self_consistent(0.3d0,3000,1d-5)
	if(.not.phy%conv) then
		write(*,*)"not converge!!!!"
	else
		rewind(20)
		do i=1,n
			write(20,"(*(e28.20))")phy%Delta(i)
		enddo
	endif
	call phy%export(10)
end program
