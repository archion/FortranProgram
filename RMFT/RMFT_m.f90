module pmt
	use M_const
	implicit none
	real(8) :: t(3)=(/1d0,-0.2d0,0.1d0/),&
		V=0.02d0,DJ=0.3d0
	integer :: cs=2
end module
module selfcons
	use pmt
	use M_Hamilton_m
	use M_utility
	implicit none
	include 'nlopt.f'
	type(t_brizon), target :: brizon_sc
contains
	subroutine initial()
		integer :: l2,l3,l,i
		real(8), allocatable :: bd0(:),bd1(:)
		allocate(var(-10000:10000))
		call init_random_seed()

		! lattice 
		latt%a1=(/1d0,0d0,0d0/)
		latt%a2=(/0d0,1d0,0d0/)
		latt%T1=(/1d0,1d0,0d0/)
		latt%T2=(/-1d0,1d0,0d0/)
		latt%bdc(1)=1d0
		latt%bdc(2)=1d0
		allocate(latt%sub(1,3))
		latt%sub(1,:)=(/0d0,0d0,0d0/)
		latt%layer=1
		call latt%gen_latt()
		call latt%gen_bond(size(t))
		!brizon_sc%n1=32
		brizon_sc%n1=64
		brizon_sc%n2=brizon_sc%n1
		call latt%gen_brizon(brizon_sc)
		brizon%n1=128
		brizon%n2=brizon%n1
		call latt%gen_brizon(brizon)
		!call latt%gen_origin_brizon((/1d0,0d0,0d0/),(/0d0,1d0,0d0/),o_brizon)
		!call check_lattice(101)
		write(*,*)"Total site number is: ",latt%Ns

		! cp
		call gen_var(sg=1,nb=0)
		var(iv(0))%val(1)=-4.1878d-02
		var(iv(0))%bd=-1d0

		! d-wave sc
		call gen_var(sg=2,nb=1,V=-1d0)
		do i=1,size(var(iv(0))%bd)
			var(iv(0))%bd(i)=dwave(i)
		enddo
		var(iv(0))%val=1d-4

		!! ddw
		!call gen_var(sg=13,nb=1,V=-1d0)
		!do i=1,size(var(iv(0))%bd)
			!var(iv(0))%bd(i)=img*ab(latt%sb(1)%nb(var(iv(0))%nb)%bd(i)%i(1))*dwave(i)
		!enddo
		!var(iv(0))%val=4d-2

		!! sdw
		!call gen_var(sg=4,nb=0,V=1d0,Vnb=1)
		!do i=1,size(var(iv(0))%bd)
			!var(iv(0))%bd(i)=ab(i)
		!enddo
		!var(iv(0))%val=4d-2


		! bond order
		call gen_var(sg=3,nb=1,V=-1d0)
		var(iv(0))%bd=1d0
		var(iv(0))%val=0d0

		call gen_var(sg=-3,nb=1)
		var(iv(0))%bd=-t(1)
		var(iv(0))%val=0d0

		do l=2,size(t)
			call gen_var(sg=3,nb=l)
			var(iv(0))%bd=-t(l)
			var(iv(0))%val=1d0

			!call gen_var(sg=5,nb=l)
			!do i=1,size(var(iv(0))%bd)
				!var(iv(0))%bd(i)=-t(l)*ab(latt%bond(l)%bd(i)%i(1))
			!enddo
			!var(iv(0))%val=0d0
		enddo


		call var_shrink()
		!call export_data(30)
		free_energy => Havg
		self_consist => self_consist_nlopt
		mat_diag => diag4
	end subroutine
	function ab(i)
		integer :: i
		real(8) :: ab
		ab=(-1d0)**(mod(nint(sum(latt%sb(1)%nb(0)%bd(i)%r(:2))),2))
	end function
	subroutine export_data(ut)
		integer :: ut,l,n
		do l=2,size(var(1:))
			do n=1,size(var(l)%bd)
				write(ut,"(es17.9$)")latt%sb(var(l)%sb)%nb(var(l)%nb)%bd(n)%r,&
					latt%sb(var(l)%sb)%nb(var(l)%nb)%bd(n)%dr,&
					var(l)%val(var(l)%bd2v(n)),&
					var(l)%bd(n)
				write(ut,"(x)")
			enddo
			write(ut,"(x/)")
		enddo
	end subroutine
	function gt(dt,ci,mi,mj,sigma)
		real(8) :: ci,dt,mi,mj,sigma
		real(8) :: gt
		real(8) :: X2,a,X,dp
		dp=abs(1-nf)
		X2=2d0*(dt**2+ci**2)
		X=2d0*dp**2*(dt**2-ci**2)-4d0*mi*mj*X2+X2**2/2d0
		select case(cs)
		case(4)
			a=1d0+4d0*X/((1d0-dp**2+4d0*mi**2)*(1d0-dp**2+4d0*mj**2))
		case(1,2,3)
			X2=0d0
			a=1d0
		end select
		gt=sqrt((2d0*dp*(1d0-dp))**2/((1d0-dp**2+4d0*mi**2)*(1d0-dp**2+4d0*mj**2)))*a*&
			0.5d0*(sqrt(((1d0+dp+2d0*mi)*(1d0+dp+2d0*mj))/((1d0+dp-2d0*mi)*(1d0+dp-2d0*mj)))+sigma*sqrt(((1d0+dp-2d0*mi)*(1d0+dp-2d0*mj))/((1d0+dp+2d0*mi)*(1d0+dp+2d0*mj))))
	end function
	function gxy(dt,ci,mi,mj)
		real(8) :: ci,dt,mi,mj
		real(8) :: gxy
		real(8) :: X2,a,X,dp
		dp=abs(1-nf)
		X2=2d0*(dt**2+ci**2)
		X=2d0*dp**2*(dt**2-ci**2)-4d0*mi*mj*X2+X2**2/2d0
		select case(cs)
		case(2,4)
			a=1d0+4d0*X/((1d0-dp**2+4d0*mi**2)*(1d0-dp**2+4d0*mj**2))
		case(1)
			X2=0d0
			a=1d0
		case(3)
			a=1d0
		end select
		gxy=(2d0*(1d0-dp))**2/((1d0-dp**2+4d0*mi**2)*(1d0-dp**2+4d0*mj**2)*a**7)
	end function
	function gz(dt,ci,mi,mj)
		real(8) :: ci,dt,mi,mj
		real(8) :: gz
		real(8) :: X2,a,X,dp
		dp=abs(1-nf)
		X2=2d0*(dt**2+ci**2)
		X=2d0*dp**2*(dt**2-ci**2)-4d0*mi*mj*X2+X2**2/2d0
		select case(cs)
		case(2,4)
			a=1d0+4d0*X/((1d0-dp**2+4d0*mi**2)*(1d0-dp**2+4d0*mj**2))
		case(1)
			X2=0d0
			a=1d0
		case(3)
			a=1d0
		end select
		if(abs(X2-4d0*mi*mj)<1d-7) then
			gz=(2d0*(1d0-dp))**2/((1d0-dp**2+4d0*mi**2)*(1d0-dp**2+4d0*mj**2)*a**7)
		else
			gz=(2d0*(1d0-dp))**2/((1d0-dp**2+4d0*mi**2)*(1d0-dp**2+4d0*mj**2)*a**7)*(X2-4d0*mi*mj*(1d0+6d0*X2*(1d0-dp)**2/(sqrt((1d0-dp**2+4d0*mi**2)*(1d0-dp**2+4d0*mj**2))*a**3))**2)/(X2-4d0*mi*mj)
		endif
	end function
	function gutzwiller(l,n)
		integer :: l,n
		real(8) :: gutzwiller
		real(8) :: Vbd
		real(8) :: gt_,gxy_,gz_
		real(8) :: dt,ci,mi,mj
		integer :: lp,i,j
		dt=0d0
		ci=0d0
		mi=0d0
		mj=0d0
		Vbd=1d0
		i=latt%sb(var(l)%sb)%nb(var(l)%Vnb)%bd(n)%i(1)
		j=latt%sb(var(l)%sb)%nb(var(l)%Vnb)%bd(n)%i(2)
		do lp=lbound(var,1),ubound(var,1)
			select case(var(lp)%sg)
			case(2,-2)
				if(var(lp)%nb==var(l)%Vnb) dt=var(lp)%val(var(lp)%bd2v(n))
			case(3)
				if(var(lp)%nb==var(l)%Vnb) ci=ci+var(lp)%val(var(lp)%bd2v(n))**2
			case(13)
				if(var(lp)%nb==var(l)%Vnb) ci=ci+var(lp)%val(var(lp)%bd2v(n))**2
			case(4,-4)
				mi=var(lp)%val(var(lp)%bd2v(i))*var(lp)%bd(i)
				mj=var(lp)%val(var(lp)%bd2v(j))*var(lp)%bd(j)
			end select
		enddo
		ci=sqrt(ci)
		gxy_=gxy(dt,ci,mi,mj)
		gz_=gz(dt,ci,mi,mj)
		select case(var(l)%sg)
		case(2,-2)
			if(var(l)%nb==1) Vbd=(DJ*(gxy_/2d0+gz_/4d0)-V)
		case(3)
			if(var(l)%nb==1) then
				Vbd=(DJ*(gxy_/2d0+gz_/4d0)+V)
			else
				Vbd=gt(dt,ci,mi,mj,1d0)
			endif
		case(-3)
			Vbd=gt(dt,ci,mi,mj,1d0)
		case(5)
			Vbd=gt(dt,ci,mi,mj,-1d0)
		case(13)
			Vbd=(DJ*(gxy_/2d0+gz_/4d0)+V)
		case(4,-4)
			Vbd=gz_*DJ
		end select
		gutzwiller=Vbd
	end function
	subroutine checkorder(l,rg)
		integer :: l
		real(8) :: rg(:),fE
		real(8) :: x(sum(var(2:)%n)),grad(size(x)),f_data(1)
		call self_consist()
		read(*,*)
		x=put(var(2:))
		x(l-1)=rg(1)
		do while(x(l-1)<=rg(2))
			call nlopt_fn(fE, size(x), x, grad, 0, f_data)
			x(l-1)=x(l-1)+rg(3)
			write(*,"(es12.4$)")x(l-1)
			write(*,"(es17.9$)")fE
			write(*,"(x)")
			write(110,"(es12.4$)")x(l-1)
			write(110,"(es17.9$)")fE
			write(110,"(x)")
		enddo
	end subroutine
	subroutine self_consist_nlopt(fE)
		real(8), optional :: fE
		real(8) :: x(sum(var(2:)%n)),minf
		real(8) :: lb(size(x)),ub(size(x)),cst(1),n
		integer(8) opt
		integer ires
		if(size(x)==0) then
			call nlopt_fn(minf, 0, (/1d0/), (/1d0/), 1, (/1d0/))
			return
		endif
		brizon_k => brizon_sc%k

		lb(1:)=1d-4
		!lb(1:)=0d0
		ub(1:)=1.5d0
		!lb(1)=1d-6
		!x=abs(var(2:)%put())
		x=put(var(2:))
		!call myfunc(minf,size(x),x)
		!call random_number(x)
		x=1d-1

		call nlo_create(opt, NLOPT_LN_BOBYQA, size(x))
		call nlo_set_lower_bounds(ires, opt, lb)
		call nlo_set_upper_bounds(ires, opt, ub)
		call nlo_set_min_objective(ires, opt, nlopt_fn, cst)
		!call nlo_set_xtol_abs(ires, opt, 1d-5)
		call nlo_set_xtol_rel(ires, opt, 1d-7)
		call nlo_set_maxeval(ires, opt, size(x)*200)

		call nlo_optimize(ires, opt, x, minf)
		if(ires<0) then
			write(*,"(A$)")"ires is negetive"
			write(*,"(i3$)")ires
			write(*,"(es12.4)")cst
			!stop
			!x=var(2:)%put()
			!call nlo_optimize(ires, opt, x, minf)
			!if(cst(1)>1d-4) then
				!stop
			!endif
		endif
		call nlo_destroy(opt)
		if(present(fE)) fE=minf
		!call var(2:)%get(x)
		!write(*,"(es12.4$)")x
		!write(*,"(x)")
		!call var%update()
		!write(*,"(es12.4$)")var(2:)%val(1)
		!write(*,"(x)")
		!write(*,"(A$)")"------------------->"
		write(*,"(es12.4$)")Tk,minf,x
		call update_var(x)
		write(*,"(x)")
		brizon_k => brizon%k
		!stop
	end subroutine
	subroutine nlopt_fn(val, n, x, grad, need_gradient, f_data)
		integer :: n
		real(8) :: val, x(n)
		real(8) :: grad(n),f_data(1)
		integer :: need_gradient
		integer :: info,i,j,t
		real(8) :: mu(1),vmu(1),wa(16),v(n)
		call get(var(2:),x)
		mu(1)=var(1)%val(1)-0.01d0
		t=0
		do
			info=1
			call hybrd1(minpack_fn,1,mu(1:1),vmu(1:1),1d-7,info,wa(:16),16)
			if(info/=1.or.abs(vmu(1))>1d-6) then
				t=t+1
				call random_number(mu)
				mu=(mu-0.5d0)*2d0
				if(t>100) then
					write(*,"(A$)")"chemical portianal is err"
					write(*,"(es12.4$)")var(:)%val(1)
					write(*,"(i3$)")info
					write(*,"(x)")
					stop
				endif
			else
				exit
			endif
		enddo
		call get(var(1:1),mu)
		if(size(x)==0) then
			return
		endif
		call MF_val(2,v,val,1)
		!write(*,"(es12.4$)")val,v
		!call var(2:)%get(x)
		!call MF_val(2,v,val,0)
		!val=sum(abs(v))
		!write(*,"(es12.4$)")val,var(:)%val(1)
		if(need_gradient/=0) then
			grad=v
			!write(*,"(es12.4$)")grad
		endif
		f_data(1)=sum(abs(v))/size(v)
		!write(*,"(es12.4$)")f_data(1)
		!write(*,"(es12.4$)")v
		!write(*,"(x)")
	end subroutine
	subroutine minpack_fn(n,x,v,info)
		integer, intent(in) :: n
		integer, intent(in) :: info
		real(8), intent(inout) :: x(n),v(n)
		real(8) :: fE,px(sum(var(2:)%n))
		common fE
		px=put(var(2:))
		call get(var(1:),x)
		call MF_val(1,v,fE,0)
		if(n==1) then
			call get(var(2:),px)
		endif
		!write(*,"(es12.4$)")x,v
		!write(*,"(x)")
	end subroutine
	function Havg()
		real(8) :: Havg
		real(8) :: dp,m,gta
		integer :: l,n,i,j,vi,vj
		!<E>-TS
		Havg=0d0
		do l=lbound(var,1),ubound(var,1)
			if(var(l)%sg==1) then
				cycle
			endif
			if(var(l)%sg==5) then
				gta=sign(1d0,gutzwiller(l,1))
			else
				gta=1d0
			endif
			do n=1,size(latt%sb(var(l)%sb)%nb(var(l)%Vnb)%bd)
				Havg=Havg+Eavg(l,n)*abs(gutzwiller(l,n))*gta
			enddo
		enddo

		!TdS
		dp=abs(1d0-nf)
		do n=1,latt%Ns
			if(minval(abs(abs(var(:)%sg)-4))==0) then
				l=minloc(abs(abs(var(:)%sg)-4),1)
				vi=var(l)%bd2v(n)
				m=var(l)%val(vi)*var(l)%bd(n)
			else
				m=0d0
			endif
			Havg=Havg+Tk*(dp*log(4d0*min(dp,1d-10)/((1d0+dp)**2+4d0*m**2))+nf*log(2d0*nf/(1d0-dp**2+4d0*m**2)))
		enddo

		Havg=Havg/latt%Ns
	end function
	subroutine update_var(x)
		integer :: l1,l
		real(8) :: x(:),scal,px(sum(var%n))
		px=put(var(:))
		do l1=lbound(var,1),ubound(var,1)
			if(var(l1)%sg==-3.and.var(l1)%nb==1) then
				l=l1
			endif
			if(var(l1)%sg==3.and.var(l1)%nb==1) then
				var(l1)%val(1)=var(l1)%val(1)+1d-6
				var(l)%val(1)=var(l1)%val(1)
				scal=Havg()
				var(l1)%val(1)=var(l1)%val(1)-2d-6
				var(l)%val(1)=var(l1)%val(1)
				scal=scal-Havg()
				scal=abs(scal/(8d0*2d-6*x(l1-1)))
				exit
			endif
		enddo
		write(*,"(', scal='es11.4$)")scal
		px=px*scal
		call get(var(:),px)
		call get(var(2:),x*scal)
	end subroutine
	subroutine raman_k(ut,gm,k,omg)
		integer :: ut
		real(8) :: gm,omg,k(:,:)
		complex(8) :: H(latt%Ns*spin,latt%Ns*spin),Bg(size(H,1),size(H,2),2),Bgk(size(H,1),size(H,2),2),BgU(size(H,1),size(H,2),2)
		complex(8) :: cH(size(H,1),size(H,1)),eH(size(H,1),size(H,1)),ceH(size(H,1),size(H,1))
		real(8) :: E(size(H,1)),f(size(E)),bt,q(3)
		complex(8) :: R(size(k,1),2),expm(size(H,1),size(H,2))
		integer :: i,j,l,n,m,nk,info
		bt=1d0/Tk
		R=0d0
		nk=size(brizon%k,1)
		!$OMP PARALLEL DO PRIVATE(E,H,cH,f,Bg,BgU,Bgk,q,eH,ceH)
		do n=1,nk
			call Hamilton(var,H,brizon%k(n,:))
			call mat_diag(H,E)
			cH=transpose(conjg(H))
			f=1d0/(1d0+exp(bt*E))
			call Hamilton(var,Bg(:,:,1),brizon%k(n,:),b1g)
			call Hamilton(var,Bg(:,:,2),brizon%k(n,:),b2g)
			do m=1,2
				BgU(:,:,m)=matmul(cH,matmul(Bg(:,:,m),H))
			enddo
			do l=0,size(k,1)/nk-1
				q=k(l*nk+n,:)-brizon%k(n,:)
				eH=0d0
				do i=1,latt%Ns
					do j=1,latt%Ns
						eH(i,j)=exp(-img*sum(q*(latt%i2r(i,:)-latt%i2r(j,:))))/latt%Ns
					enddo
				enddo
				eH(latt%Ns+1:,latt%Ns+1:)=eH(:latt%Ns,:latt%Ns)
				eH=matmul(eH,H)
				ceH=transpose(conjg(eH))
				do m=1,2
					Bgk(:,:,m)=matmul(ceH,matmul(Bg(:,:,m),eH))
				enddo
				do i=1,size(E)
					do j=1,size(E)
						R(l*nk+n,:)=R(l*nk+n,:)+Bgk(i,j,:)*BgU(j,i,:)*(f(i)-f(j))/(omg+E(j)-E(i)+img*gm)
					enddo
				enddo
			enddo
		enddo
		!$OMP END PARALLEL DO
		!write(*,*)sum(imag(R(:,1)))/(size(brizon%k,1)*latt%Ns)
		do i=1,size(k,1)
			write(ut,"(e16.5$)")k(i,:),R(i,:)/(size(brizon%k,1)*latt%Ns)
			write(ut,"(1X)")
		enddo
		write(ut,"(1X/)")
	end subroutine
	subroutine raman(ut,gm,omg,m,peak)
		integer :: ut,m
		real(8) :: gm,omg(:)
		real(8), optional, allocatable :: peak(:)
		complex(8) :: H(latt%Ns*spin,latt%Ns*spin),Bg(size(H,1),size(H,2),2)
		complex(8) :: cH(size(H,1),size(H,1))
		real(8) :: E(size(H,1)),f(size(E)),domg,bt
		complex(8) :: R(m,2),Rsc(m,2),Rpg(m,2)
		integer :: i,l,n,m1,m2
		integer, allocatable :: ipeak(:)
		bt=1d0/Tk
		domg=(omg(2)-omg(1))/m
		R=0d0
		Rsc=0d0
		Rpg=0d0
		!$OMP PARALLEL DO REDUCTION(+:R,Rsc,Rpg) PRIVATE(E,H,cH,f,Bg)
		do n=1,size(brizon%k,1)
			call Hamilton(var,H,brizon%k(n,:))
			call mat_diag(H,E)
			cH=transpose(conjg(H))
			f=1d0/(1d0+exp(bt*E))

			call Hamilton(var(:),Bg(:,:,1),brizon%k(n,:),b1g)
			call Hamilton(var(:),Bg(:,:,2),brizon%k(n,:),b2g)
			do l=1,2
				Bg(:,:,l)=matmul(cH,matmul(Bg(:,:,l),H))
			enddo
			do m1=1,4
				do m2=1,4
					do l=1,m
						R(l,:)=R(l,:)+Bg(m1,m2,:)*Bg(m2,m1,:)*(f(m1)-f(m2))/(omg(1)+domg*l+E(m2)-E(m1)+img*gm)
						if(abs(abs(E(m1))-abs(E(m2)))>1d-10) then
							Rpg(l,:)=Rpg(l,:)+Bg(m1,m2,:)*Bg(m2,m1,:)*(f(m1)-f(m2))/(omg(1)+domg*l+E(m2)-E(m1)+img*gm)
						else
							Rsc(l,:)=Rsc(l,:)+Bg(m1,m2,:)*Bg(m2,m1,:)*(f(m1)-f(m2))/(omg(1)+domg*l+E(m2)-E(m1)+img*gm)
						endif
					enddo
				enddo
			enddo
		enddo
		!$OMP END PARALLEL DO
		do i=1,m
			write(ut,"(es17.9$)")omg(1)+domg*i,R(i,:)/(size(brizon%k,1)*latt%Ns),Rsc(i,:)/(size(brizon%k,1)*latt%Ns),Rpg(i,:)/(size(brizon%k,1)*latt%Ns)
			write(ut,"(1X)")
		enddo
		write(ut,"(1X)")
		if(present(peak)) then
			if(allocated(peak)) then
				deallocate(peak)
			endif
			call find_peak(imag(R(:,1)),ipeak)
			allocate(peak(size(ipeak)))
			do i=1,size(ipeak)
				peak(i)=omg(1)+domg*ipeak(i)
			enddo
		endif
	end subroutine
	function b1g(dr)
		real(8) :: dr(3)
		complex(8) :: b1g
		b1g=(-dr(1)**2+dr(2)**2)/2d0
	end function
	function b2g(dr)
		real(8) :: dr(3)
		complex(8) :: b2g
		b2g=(-dr(1)*dr(2))
	end function
end module
program main
	use selfcons
	implicit none
	logical :: f
	integer :: l,m,i,od
	real(8) :: n(0:10)=(/(min(1d0-0.005d0*i,0.999d0),i=0,10)/)
	real(8) :: Ts(size(n),2),Td(size(n),2),Tc(2),Tnc,fE,kf(3),Vbd,pnf,dnf,pk
	real(8), allocatable :: peak(:)
	logical :: flag
	open(unit=10,file=fn('../data/phase.dat'))
	!open(unit=20,file=fn('../data/band.dat'))
	!f=openfile(unit=30,file='../data/order.dat')
	!open(unit=40,file=fn('../data/fermis.dat'))
	!f=openfile(unit=50,file='../data/DOS.dat')
	!open(unit=60,file=fn('../data/raman.dat'))
	!f=openfile(unit=70,file='../data/map_raman.dat')
	!f=openfile(unit=80,file='../data/map_band.dat')
	!f=openfile(unit=90,file='../data/energy.dat')
	!open(unit=100,file=fn('../data/gap.dat'))
	!f=openfile(unit=101,file='../data/lattice.dat')
	!open(unit=110,file=fn('../data/peak.dat'))

	flag=.true.
	call initial()

	!Tk=0.095d0
	!nf=1d0-0.15d0
	!call checkorder(3,(/0d0,0.05d0,0.0001d0/))
	!stop

	!Tk=1d-4
	!do l=0,ubound(n,1),2
		!nf=n(l)
		!call self_consist(fE)
		!!call fermis(40,0.01d0,o_brizon%k,0d0)
		!write(*,"(es12.4$)")nf,Tk,fE,var(2:)%val(1)
		!write(*,"(x)")
		!write(100,"(es12.4$)")nf,Tk,fE,var(2:)%val(1)
		!write(100,"(x)")
		!!call band(20,(/0d0,0d0,0d0/),(/pi,pi,0d0/),128)
		!!call band(20,(/pi,pi,0d0/),(/0d0,pi,0d0/),128)
		!!call band(20,(/0d0,pi,0d0/),(/0d0,0d0,0d0/),128)
		!!call raman(60,0.01d0,(/0d0,0.8d0/),256,peak)
		!!call DOS(50,0.01d0,(/-1.8d0,1.8d0/),256)
		!!write(*,*)peak
		!!call raman_k(70,0.01d0,o_brizon%k,peak(1))
		!!stop
	!enddo
	!stop
	!do l=0,10
		!nf=8.6445d-01
		!Tk=max(5.4570d-02/10*l,1d-4)
		!call self_consist()
		!write(*,"(es12.4$)")nf,Tk,var(2:3)%val(1)*8d0
		!write(*,"(x)")
		!write(100,"(es12.4$)")nf,Tk,var(2:3)%val(1)*8d0
		!write(100,"(x)")
		!!call raman(60,0.05d0,(/0d0,0.8d0/),256,peak)
		!call raman(60,0.008d0,(/0d0,0.5d0/),256,peak)
		!if(l==1) then
			!pk=peak(1)
		!endif
		!write(110,"(es12.4$)")nf,Tk,0.1d0*l,peak(1),peak(1)/pk
		!write(110,"(x)")
		!!!call DOS(50,0.01d0,(/-1.8d0,1.8d0/),256)
		!!!call fermis(40,0.01d0,o_brizon%k,0d0)
		!!!stop
		!!!call band(20,(/0d0,0d0,0d0/),(/pi,pi,0d0/),128)
		!!!call band(20,(/pi,pi,0d0/),(/0d0,pi,0d0/),128)
		!!!call band(20,(/0d0,pi,0d0/),(/0d0,0d0,0d0/),128)
		!!kf=(/pi,3.926990817d-01,0d0/)
		!!call band(20,kf,kf,1)
	!enddo
	!!call band(20,(/pi,0d0,0d0/),(/pi,pi,0d0/),128)
	!!call fermis(40,0.01d0,o_brizon%k,0d0)
	!stop

	!od=2
	!Tc=(/1d-1,1d-4/)
	!do l=0,ubound(n,1),2
		!nf=n(l)
		!Tk=Tc(1)
		!Tc(1)=find_order(od,-1,Tk,(/1d-4,0.25d0,5d-3/),1d-3)
		!Tk=Tc(2)
		!Tc(2)=find_order(od,1,Tk,(/1d-4,Tc(1),5d-3/),1d-3)
		!write(10,"(es12.4$)")nf,Tc
		!write(*,"(es12.4$)")nf,Tc
		!write(10,"(x)")
		!write(*,"(x)")
	!enddo
	!write(10,"(x/)")

	od=2
	Tc=(/1d-1,1d-4/)
	do l=0,ubound(n,1),2
		nf=n(l)
		Tk=Tc(1)
		Tc(1)=find_order(od,-1,Tk,(/1d-4,0.3d0,5d-3/),1d-3)
		Tk=Tc(2)
		Tc(2)=find_order(od,1,Tk,(/1d-4,Tc(1),5d-3/),1d-3)
		if(Tc(1)<2d-4) then
			call swap(nf,pnf)
			dnf=abs(nf-pnf)
			Tc(1)=1d-1
			do
				nf=nf-dnf/10d0
				Tk=Tc(1)
				Tc(1)=find_order(od,-1,Tk,(/1d-4,0.25d0,5d-3/),1d-3)
				Tk=Tc(2)
				Tc(2)=find_order(od,1,Tk,(/1d-4,Tc(1),5d-3/),1d-3)
				write(10,"(es12.4$)")nf,Tc
				write(*,"(es12.4$)")nf,Tc
				write(10,"(x)")
				write(*,"(x)")
				if((nf<=pnf).or.(Tc(1)<2d-4).or.(Tc(1)>2d-4.and.(Tc(1)-Tc(2))<5d-4)) then
					exit
				endif
			enddo
			exit
		endif
		pnf=nf
		if(Tc(2)>2d-4.and.flag) then
			flag=.false.
			Tk=1d-4
			nf=find_order(od,1,nf,(/nf-0.1d0,1d0,2d-2/),1d-3)
			Tk=0.5d0*(Tc(1)+Tc(2))
			Tnc=find_order(2,-1,Tk,(/1d-4,0.25d0,5d-3/),1d-3)
			Tk=find_order(od,-1,Tk,(/1d-4,0.25d0,5d-3/),1d-3)
			write(10,"(es12.4$)")nf,Tk,3d-4
			write(*,"(es12.4$)")nf,Tk,3d-4
			write(10,"(x)")
			write(*,"(x)")
			!do i=1,10
				!Tk=Tnc/10d0*i
				!call self_consist()
				!write(*,"(es12.4$)")nf,Tk,var(2:3)%val(1)*8d0
				!write(*,"(x)")
				!write(100,"(es12.4$)")nf,Tk,var(2:3)%val(1)*8d0
				!write(100,"(x)")
				!call raman(60,0.05d0,(/0d0,0.8d0/),256,peak)
			!enddo
		endif
		write(10,"(es12.4$)")n(l),Tc
		write(*,"(es12.4$)")n(l),Tc
		write(10,"(x)")
		write(*,"(x)")
	enddo
	write(10,"(x/)")


	!Tc=(/1.5d-1,1d-4/)
	!!do l=18,size(n)-28
	!do l=1,size(n)-28,2
		!nf=n(l)
		!Tk=Tc(1)
		!Tc(1)=find_order(4,-1,Tk,(/1d-4,0.25d0,5d-3/),1d-3)
		!Tk=Tc(2)
		!Tc(2)=find_order(4,1,Tk,(/1d-4,Tc(1),5d-3/),1d-3)
		!write(10,"(es12.4$)")nf,Tc
		!write(*,"(es12.4$)")nf,Tc
		!write(10,"(x)")
		!write(*,"(x)")
		!if(Tc(1)<2d-4) then
			!exit
		!endif
	!enddo
	!write(10,"(x/)")
end program
