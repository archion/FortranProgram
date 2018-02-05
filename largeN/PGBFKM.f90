module global
	use M_const
	use M_matrix
	use mkl_service
	implicit none
	real(8) :: alpha=0.1d0,r=0.3d0,kp=0.5d0,fcut=1d0,fbeta=1d0,phicut=5d-2,phibeta=200d0,Jk=1.5d0,g=0d0
	!real(8) :: alpha=0.3d0,r=0.3d0,kp=0.3d0,fcut=1d0,fbeta=20d0,phicut=5d-2,phibeta=200d0,Jk=1.5d0,g=0.2d0
	!real(8) :: alpha=0.3d0,r=0d0,kp=0.5d0,fcut=1d0,fbeta=20d0,phicut=5d-2,phibeta=200d0,Jk=2d0,g=0.5d0
	!real(8) :: alpha=0.3d0,r=0d0,kp=0.3d0,fcut=1d0,fbeta=1d0,phicut=5d-2,phibeta=200d0,Jk=1d0,g=0d0
	real(8) :: Tk=4d-3,beta,norm_Ac=1d0,norm_Aphi=1d0
	integer, parameter :: nll=31*7,nadd(2)=[5,10],n=(nll+sum(nadd*2))*2
	!integer, parameter :: nll=33*14,nadd1=5,nadd2=10,n=(nll+nadd1*2+nadd2*2)*2
	!integer, parameter :: nll=31*15*2,nadd1=50,nadd2=10,n=(nll+nadd1*2+nadd2*2)*2
	!integer, parameter :: n=536
	!real(8) :: omega((n+nadd1*2+nadd2*2)*2)
	real(8) :: omega(n),tau(n)
	logical :: is_expinit=.true.
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
			A_phi=norm_Aphi*abs(x)**(1d0-alpha)*sign(1d0,x)
		else
			A_phi=norm_Aphi*phicut**(1d0-alpha)*sign(1d0,x)*1d0/((exp(phibeta*(x-phicut-4d0/phibeta))+1d0)*(exp(phibeta*(-x-phicut-4d0/phibeta))+1d0))
		endif
	end function
    real(8) function A_c(x)
		real(8) :: x
		if(abs(x)<fcut) then
			!A_c=norm_Ac*(abs(x)/fcut)**r*exp(-fcut**2/pi*fbeta)
			A_c=norm_Ac*(abs(x)/fcut)**r
			!A_c=norm_Ac*abs(x)**r
			!A_c=norm_Ac*abs(x)**r
		elseif(abs(x)<=(-omega(1))) then
			!A_c=norm_Ac*exp(-x**2/pi*fbeta)
			A_c=norm_Ac*exp(-(x**2-fcut**2)/pi*fbeta)
			!A_c=norm_Ac*abs(fcut)**r*exp(-fbeta*(abs(x)-fcut)**2)
			!A_c=norm_Ac*fcut**r*ff(fbeta*(x-fcut-15.0/fbeta))*ff(fbeta*(-x-fcut-15.0/fbeta))
			!A_c=norm_Ac*fcut**r*1d0/(exp(fbeta*(x-fcut-15.0/fbeta))+1d0)*1d0/(exp(fbeta*(-x-fcut-15.0/fbeta))+1d0)
		else
			A_c=0d0
		endif
	end function
	real(8) function integrate(Af,A,l,u) result(weight)
		real(8), optional :: Af,A(:)
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
		if(present(Af)) then
			do i=il,iu
				weight=weight+0.5d0*(omega(i+1)-omega(i))*(Af(omega(i+1))+Af(omega(i)))
			enddo
		else
			do i=il,iu
				weight=weight+0.5d0*(omega(i+1)-omega(i))*(A(i+1)+A(i))
			enddo
		endif
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
		re=0.5d0/pi*re
	end subroutine
	pure real(8) function ff(x)
		real(8), intent(in) :: x
		!if(x<omega(1)) then
			!ff=1d0
		!elseif(x>-omega(1)) then
			!ff=0d0
		!else
			ff=1d0/(exp(x)+1d0)
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
				fb=1d0/(exp(x)-1d0)
			!endif
		endif
	end function
	subroutine get_ftau(omega,A,Af,eta,time,rt)
		real(8) :: omega(:),time(:),rt(:)
		real(8), optional :: A(:),Af
		integer :: eta
		integer :: tau,i,k,n
		real(8) :: A1,A2,B1,B2,e1,e2,t
		procedure(real(8)), pointer :: f
		n=size(omega)
		tau=1
		if((time(size(time))+time(1))>0d0) then
			tau=-1
		endif
		do k=1,size(time)
			t=time(k)
			e1=omega(1)
			if(present(A)) then
				A1=A(1)
			else
				A1=Af(e1)
			endif
			B1=tau*1d0/(exp((tau*beta+t)*e1)-eta*exp(e1*t))
			rt(k)=0d0
			do i=1,n-1
				e2=omega(i+1)
				if(present(A)) then
					A2=A(i+1)
				else
					A2=Af(e2)
				endif
				B2=tau*1d0/(exp((tau*beta+t)*e2)-eta*exp(e2*t))
				rt(k)=rt(k)+0.5d0*(e2-e1)*(A1*B1+A2*B2)
				A1=A2
				B1=B2
				e1=e2
			enddo
			rt(k)=-eta*rt(k)
		enddo
	end subroutine
	subroutine get_convolution(omega,A,B,Af,Bf,tau,eta,rt)
		real(8) :: omega(:),rt(:)
		real(8), optional :: A(:),B(:),Af,Bf
		integer :: tau(2),eta(2)
		integer :: i,j,k,n,iwb(size(omega)),ibw(size(omega)),m
		real(8) :: e1,e2,A1,A2,B1,B2,omg
		procedure(real(8)), pointer :: f1,f2
		if(.not.(present(A).xor.present(Af))) stop "error"
		if(.not.(present(B).xor.present(Bf))) stop "error"
		n=size(omega)
		select case(eta(1))
		case(1)
			f1 => fb
		case(-1)
			f1 => ff
		end select
		select case(eta(2))
		case(1)
			f2 => fb
		case(-1)
			f2 => ff
		end select
		!$OMP PARALLEL DO PRIVATE(omg,ibw,iwb,A1,B1,A2,B2,e1,e2,m)
		do k=1,n/2
		!do k=1,n
			omg=omega(k)
			if(present(A)) then
				do i=1,n
					ibw(i)=find_sidx(omega,tau(1)*(omg-tau(2)*omega(i)))
				end do
				A1=A(1)
			else
				A1=Af(omega(1))
			endif
			do i=1,n
				iwb(i)=find_sidx(omega,tau(2)*(omg-tau(1)*omega(i)))
			end do
			if(present(B)) then
				i=iwb(1)
				B1=merge(B(i)+(B(i+1)-B(i))/(omega(i+1)-omega(i))*(tau(2)*(omg-tau(1)*omega(1))-omega(i)),0d0,i>0.and.i<n)
			else
				B1=Bf(tau(2)*(omg-tau(1)*omega(1)))
			endif
			e1=omega(1)
			rt(k)=0d0
			do i=1,n-1
				do j=iwb(i),iwb(i+1),sign(1,iwb(i+1)-iwb(i))
					if(j==iwb(i+1)) then
						e2=omega(i+1)
						if(present(A)) then
							A2=A(i+1)
						else
							A2=Af(omega(i+1))
						endif
						if(present(B)) then
							B2=merge(B(j)+(B(j+1)-B(j))/(omega(j+1)-omega(j))*(tau(2)*(omg-tau(1)*e2)-omega(j)),0d0,j>0.and.j<n)
						else
							B2=Bf(tau(2)*(omg-tau(1)*e2))
						endif
					else
						m=j+(1-tau(1)*tau(2))/2
						e2=tau(1)*(omg-tau(2)*omega(m))
						if(present(B)) then
							B2=merge(B(m),0d0,m>0.and.m<n)
						else
							B2=Bf(omega(m))
						endif
						if(present(A)) then
							m=ibw(m)
							A2=merge(A(m)+(A(m+1)-A(m))/(omega(m+1)-omega(m))*(e2-omega(m)),0d0,m>0.and.m<n)
						else
							A2=Af(e2)
						endif
					endif
					rt(k)=rt(k)+0.5d0*(e2-e1)*&
						((-f1(beta*tau(1)*e1)+eta(1)*eta(2)*f2(beta*(tau(1)*e1-omg)))*A1*B1+&
						(-f1(beta*tau(1)*e2)+eta(1)*eta(2)*f2(beta*(tau(1)*e2-omg)))*A2*B2)
					A1=A2
					B1=B2
					e1=e2
				enddo
			enddo
			rt(k)=-pi*eta(1)*tau(1)*tau(2)*rt(k)
			rt(n+1-k)=-eta(1)*eta(2)*rt(k)
		enddo
		!$OMP END PARALLEL DO
	end subroutine
	real(8) function fenergy(Af,rGf,AB,rGB) result(rt)
		real(8) :: Af(:),AB(:),rGf(:),rGB(:)
		real(8) :: F1(size(Af)),F2(size(AB))
		integer :: i
		do i=1,size(F1)
			F1(i)=(atan(rGf(i)/(pi*Af(i)))+pi/2d0)*ff(beta*omega(i))
			F2(i)=atan(-pi*AB(i)/rGB(i))*fb(beta*omega(i))
		enddo
		!rt=-integrate(A=F1)/pi-kp*integrate(A=F2)/pi
		rt=-integrate(A=F1)/pi
	end function
	subroutine self_consistent(tol,niter,rate,Af,AB,iSEf,iSEB)
		real(8) :: tol,rate,Af(:),AB(:),iSEF(:),iSEB(:)
		integer :: niter
		real(8) :: SEf1(size(Af)),SEf2(size(Af)),iSEf_(size(Af)),rSEf(size(Af)),iSEB_(size(Af)),rSEB(size(Af)),Af_(size(Af)),AB_(size(Af))
		real(8) :: x(size(iSEf)*2),x_(size(iSEf)*2),nf(1)
		logical, save :: flag=.true.
		integer :: i,n,k
		n=size(omega)
		beta=1d0/Tk
		Af=0d0
		AB=0d0
		!call mbroyden([real(8)::],[real(8)::],1d0,0,fin=.false.,id=2)
		do i=1,niter
			call get_realpart(omega,iSEf,rSEf)
			call get_realpart(omega,iSEB,rSEB)
			Af_=-1d0/pi*iSEf/((omega-rSEf)**2+iSEf**2)
			!AB_=-1d0/pi*iSEB/((1d0/Jk-rSEB)**2+iSEB**2)
			AB_=-1d0/pi*iSEB/((-1d0/Jk-rSEB)**2+iSEB**2)
			call get_ftau(omega,A=Af_,eta=-1,time=[-1d-70],rt=nf)
			!if(all(abs(Af(n/2+1:)-Af_(n/2+1:))/(abs(Af_(n/2+1:))+1d-70)<tol).and.all(abs(AB(n/2+1:)-AB_(n/2+1:))/(abs(AB_(n/2+1:))+1d-70)<tol)) then
			!if(all(abs(Af(n/2+1:)-Af_(n/2+1:))/(abs(Af_(n/2+1:))+1d-70)<tol)) then
				!write(*,*)"converged",1d0/Jk+rSEB(n/2)
				!exit
			!endif
			!if(i==niter) then
				!write(*,*)"not converged",max(maxval(abs(Af(n/2+1:)-Af_(n/2+1:))/(Af_(n/2+1:)+1d-70)),maxval(abs(AB(n/2+1:)-AB_(n/2+1:))/(AB_(n/2+1:)+1d-70)))
			!endif
			Af=Af_
			AB=AB_

			!call get_convolution(omega,Af=A_c,B=AB,tau=[1,1],eta=[-1,1],rt=SEf1)
			!call get_convolution(omega,A=Af,Bf=A_phi,tau=[1,1],eta=[-1,1],rt=SEf2)
			!call get_convolution(omega,Af=A_c,B=Af,tau=[1,1],eta=[-1,-1],rt=iSEB_)


			call get_convolution(omega,Af=A_c,B=AB,tau=[1,1],eta=[-1,1],rt=SEf1)
			call get_convolution(omega,Af=A_phi,B=Af,tau=[1,1],eta=[1,-1],rt=SEf2)
			call get_convolution(omega,Af=A_c,B=Af,tau=[1,1],eta=[-1,-1],rt=iSEB_)

			!iSEf_=(kp*SEf1-g**2*SEf2)

			iSEf_=(-kp*SEf1-g**2*SEf2)
			iSEB_=-iSEB_

			!write(*,*)nint(sign(1d0,[SEf1(1+n/4),SEf1(n-n/4),SEf2(1+n/4),SEf2(n-n/4),iSEB_(1+n/4),iSEB_(n-n/4)]))
			!stop

			x=[iSEf,iSEB]
			x_=[iSEf_,iSEB_]
			if(all(abs(x-x_)/(abs(x_)+1d-70)<tol)) then
				write(*,*)"converged"
				if(abs(0.5d0-nf(1))>1d-2) then
					is_expinit=.false.
				endif
				exit
			endif
			if(i==niter) then
				write(*,*)"not converged",maxval(abs(x-x_)/(abs(x_)+1d-70))
				is_expinit=.false.
			endif

			!call mbroyden(x,x_,1d0,i-800,8)
			!if(flag) then
				!call mbroyden(x,x_,1d0,i,8)
				call mbroyden(x,x_,rate,i,40)
			!else
				!x=(x_*rate+x)/(1d0+rate)
			!endif


			iSEf=x(1:size(x)/2)
			iSEB=x(size(x)/2+1:)
		enddo
		!write(*,*)i,maxval(abs(Af(n/2+1:)-Af_(n/2+1:))/(abs(Af_(n/2+1:))+1d-70)),"nf:",nf*2d0
		write(*,*)i,maxval(abs(x-x_)/(abs(x_)+1d-70)),"nf:",nf*2d0
		call mbroyden([real(8)::],[real(8)::],1d0,0,fin=.true.)
		flag=.false.
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
	real(8) :: Af(n),AB(n),rGf(n),rGB(n),iSEF(n),iSEB(n),SUS(n),rSUS(n),Tmat(n),SUSt(n),Tmatt(n),Gft(n),dTk,omg(2),dat(2,2),fe(2)=0d0,pTk=0d0
	integer :: i,j,e

	call omp_set_nested(.false.)
	!call omp_set_max_active_levels(1)
	call omp_set_dynamic(.false.)
	call mkl_set_dynamic(0)
	!call omp_set_schedule(omp_sched_static,0)

	call mkl_set_num_threads(16)
	call omp_set_num_threads(mkl_get_max_threads())

	open(10,file="../data/init.dat")
	open(13,file="TEMPLIST.dat")
	open(22,File="../data/largeN_omega_my.dat")
	open(23,File="../data/largeN_tau_my.dat")
	open(24,File="../data/largeN_T.dat")
	open(25,File="../data/largeN_entropy.dat")
	call set_omega(omega(1:nll*2),2.3d0,7,10d0)
	!!call add_omega(omega(1:nll*2+nadd1*4),fcut+sqrt(abs(log(0.2d0))/fbeta),sqrt(abs(log(0.2d0))/fbeta)*2d0,nadd1)
	!!call add_omega(omega,0.8d0,0.16d0,nadd2)
	!call set_omega(omega(1:nll*2),1.5d0,15,10d0)
	call add_omega(omega(1:nll*2+sum(nadd(1:1))*4),1d0,0.1d0,nadd(1))
	!call add_omega(omega(1:nll*2+nadd1*4),1.5d0,0.5d0,nadd1)
	call add_omega(omega(1:nll*2+sum(nadd(1:2))*4),0.8d0,0.16d0,nadd(2))
	!call add_omega(omega(1:nll*2+sum(nadd(1:3))*4),9.8d0,0.1d0,nadd(3))
	write(*,"('n: ',i4,', max freq: ',es10.2,', min freq: ',es10.2)")n,maxval(abs(omega)),minval(abs(omega))

	!stop
	!iSEf=-1d0
	!iSEB(:size(omega)/2)=-1d0
	!iSEB(size(omega)/2+1:)=1d0

	!do i=1,n
		!write(22,"(*(e28.20))")omega(i),Af(i)*pi,Ab(i)*pi,iSEf(i),iSEB(i),SUS(i)/pi,Tmat(i)*pi,A_c(omega(i))*pi,merge(0d0,A_phi(omega(i))*pi,abs(A_phi(omega(i))*pi)<1d-20)
	!enddo
	!write(22,"(x/)")
	!stop

	Tk=0.09d0
	!Tk=0.01d0

	!fbeta=1d0
	!Jk=Jk*integrate(Af=A_c)/pi
	!norm_Ac=norm_Ac/integrate(Af=A_c)
	!call self_consistent(1d-5,6000,3d-3,Af,AB,iSEf,iSEB)
	!do i=1,n
		!write(22,"(*(e28.20))")omega(i),Af(i)*pi,Ab(i)*pi,iSEf(i),iSEB(i),SUS(i)/pi,Tmat(i)*pi,A_c(omega(i))*pi,merge(0d0,A_phi(omega(i))*pi,abs(A_phi(omega(i))*pi)<1d-20)
	!enddo
	!write(22,"(x/)")
	!fbeta=7d0
	!Jk=1d0
	!norm_Ac=1d0

	!norm_Aphi=norm_Aphi/(integrate(Af=A_phi,l=size(omega)/2)*2)
	!fcut=abs(omega(find_sidx(omega,-fcut)))
	!Jk=Jk*integrate(Af=A_c)/pi
	!norm_Ac=norm_Ac/integrate(Af=A_c)*Jk
	!norm_Ac=norm_Ac/integrate(Af=A_c,l=find_sidx(omega,-fcut),u=find_sidx(omega,fcut))

	!norm_Ac=(r+1d0)/(2d0*fcut)

	do 
		is_expinit=.true.
		norm_Ac=norm_Ac/integrate(Af=A_c)
		i=1
		rewind(10)
		read(10,"(*(e28.20))")omg(1),dat(1,:)
	o:	do 
			read(10,"(*(e28.20))",iostat=e)omg(2),dat(2,:)
			if(e==-1) then
				backspace(10)
				backspace(10)
				backspace(10)
				read(10,"(*(e28.20))")omg(1),dat(1,:)
				dat(2,:)=dat(2,:)+(dat(2,:)-dat(1,:))/(omg(2)-omg(1))*(-omega(1)-omg(2))
				omg(2)=-omega(1)
				read(10,"(*(e28.20))")omg(1),dat(1,:)
			endif
			do 
				if(omg(2)-omega(i)>=-1d-12) then
					dat(1,:)=dat(1,:)+(dat(2,:)-dat(1,:))/(omg(2)-omg(1))*(omega(i)-omg(1))
					iSEf(i)=dat(1,1)
					iSEB(i)=dat(1,2)
					omg(1)=omega(i)
					i=i+1
					if(i==n+1) exit o
				else
					omg(1)=omg(2)
					dat(1,:)=dat(2,:)
					exit
				endif
			enddo
		enddo o
		rewind(13)
		do
			read(13,*,iostat=e)Tk
			if(e==-1) exit
			write(*,*)"Tk: ",Tk
			!call self_consistent(1d-5,6000,3d-3,Af,AB,iSEf,iSEB)
			!if(Tk>1d-4) then
				!call self_consistent(1d-6,6000,1d-1,Af,AB,iSEf,iSEB)
			!else
				call self_consistent(1d-5,1200,5d-3,Af,AB,iSEf,iSEB)
			!endif

			call get_convolution(omega,A=Af,B=Af,tau=[1,-1],eta=[-1,-1],rt=SUS)
			call get_convolution(omega,A=Af,B=AB,tau=[1,1],eta=[-1,1],rt=Tmat)

			call get_realpart(omega,-Af*pi,rGf)
			call get_realpart(omega,-AB*pi,rGB)
			!fe(2)=fenergy(Af,rGf,AB,rGB)
			!write(*,*)"entropy: ",-(fe(2)-fe(1))/(Tk-pTk),fe(2)-fe(1)
			!write(25,"(*(e28.20))")Tk,-(fe(2)-fe(1))/(Tk-pTk),fe(2)


			do i=1,n
				write(22,"(*(e28.20))")Tk,omega(i),Af(i)*pi,Ab(i)*pi,iSEf(i),iSEB(i),SUS(i)/pi,Tmat(i)*pi,rGf(i),rGB(i),A_c(omega(i))*pi,merge(0d0,A_phi(omega(i))*pi,abs(A_phi(omega(i))*pi)<1d-20),atan(tanh((omega(i)/Tk)*0.5d0)/tan((1d0-0.785d0)*pi*0.5d0))
			enddo
			write(22,"(x/)")

			!tau=(omega(n/2+1:)+pi*Tk/10d0)*10d0
			tau(1:n/2)=asin(pi*Tk/(omega(1:n/2)/omega(1)*(100d0-pi*Tk)+pi*Tk))/(pi*Tk)
			tau(n/2+1:)=1d0/Tk-tau(1:n/2)
			call get_ftau(omega,A=SUS,eta=1,time=tau,rt=SUSt)
			call get_ftau(omega,A=Tmat,eta=-1,time=tau,rt=Tmatt)
			call get_ftau(omega,A=Af,eta=-1,time=tau,rt=Gft)

			do i=1,n/2
				write(23,"(*(e28.20))")Tk,tau(i),SUSt(i),Tmatt(i),Gft(i)
			enddo
			write(23,"(x)")
			do i=n/2+1,n
				write(23,"(*(e28.20))")Tk,tau(i),SUSt(i),Tmatt(i),Gft(i)
			enddo
			write(23,"(x/)")

			!call get_realpart(omega,SUS,rSUS)
			!write(24,"(*(e28.20))")Tk,rSUS(n/2)
			if(is_expinit) then
				write(*,*)"export initial...."
				rewind(10)
				do i=1,n
					write(10,"(*(e28.20))")omega(i),iSEf(i),iSEB(i)
				enddo
				is_expinit=.false.
			endif
			pTk=Tk
			fe(1)=fe(2)
		enddo
		write(*,*)"Jk: ",Jk,"kp: ",kp,"r: ",r,"g: ",g
		stop
		Jk=Jk+0.05d0
		if(Jk>3d0) then
			exit
		endif
		!kp=kp+0.2d0
		!if(kp>0.8d0) then
			!exit
		!endif
		!fbeta=fbeta+1d0
		!if(fbeta>8d0) then
			!exit
		!endif
		!g=g+0.5d0
		!if(g>5d0) then
			!exit
		!endif
	enddo
end program
