module pmt
	use M_const
	implicit none
	real(8), parameter :: t(2)=(/1d0,-0.3d0/),&
		V=0.12d0,DJ=0.30d0
		!V=0.04d0,DJ=0.30d0
end module
module selfcons
	use pmt
	use M_Hamilton
	use M_utility
	implicit none
	include 'nlopt.f'
contains
	subroutine initial()
		integer :: l2,l3,l,i
		real(8), allocatable :: bd0(:),bd1(:)
		allocate(var(-10000:10000))
		call init_random_seed()

		! lattice 
		latt%a1=(/1d0,1d0,0d0/)
		latt%a2=(/-1d0,1d0,0d0/)
		latt%T1=(/1d0,1d0,0d0/)
		latt%T2=(/-1d0,1d0,0d0/)
		latt%bdc(1)=1d0
		latt%bdc(2)=1d0
		allocate(latt%sub(2,3))
		latt%sub(1,:)=(/0d0,0d0,0d0/)
		latt%sub(2,:)=(/1d0,0d0,0d0/)
		latt%layer=1
		call latt%gen_latt()
		call latt%gen_neb(size(t))
		call latt%gen_bond(size(t))
		!brizon%n1=32
		brizon%n1=64
		brizon%n2=brizon%n1
		call latt%gen_brizon(brizon)
		call latt%gen_origin_brizon((/1d0,0d0,0d0/),(/0d0,1d0,0d0/),o_brizon)
		!call check_lattice(101)
		write(*,*)"Total site number is: ",latt%Ns

		! cp
		call gen_var(sg=1,nb=0)
		var(iv(0))%val(1)=0d0
		var(iv(0))%bd=-1d0

		! d-wave sc
		call gen_var(sg=2,nb=1)
		do i=1,size(var(iv(0))%bd)
			var(iv(0))%bd(i)=dwave(i)
		enddo
		do l2=1,size(var(iv(0))%val)
			var(iv(0))%val(l2)=1d-1
		enddo

		! ddw
		call gen_var(sg=13,nb=1)
		do i=1,size(var(iv(0))%bd)
			var(iv(0))%bd(i)=img*ab(latt%bond(var(iv(0))%nb)%bd(i)%i(1))*dwave(i)
		enddo
		do l2=1,size(var(iv(0))%val)
			var(iv(0))%val(l2)=1d-1
		enddo

		!! sdw
		!call gen_var(sg=4,nb=0)
		!do i=1,size(var(iv(0))%bd)
			!var(iv(0))%bd(i)=ab(i)
		!enddo
		!do l2=1,size(var(iv(0))%val)
			!var(iv(0))%val(l2)=1d-1
		!enddo


		! bond order
		do l=1,size(t)
			call gen_var(sg=3,nb=l)
			var(iv(0))%bd=-sign(1d0,t(l))
			!var(iv(0))%bd=1d0
			do l2=1,size(var(iv(0))%val)
				!var(iv(0))%val(l2)=t(l)
				var(iv(0))%val(l2)=abs(t(l))
			enddo
		enddo


		call var_shrink()
		!var(lbound(var,1))%update => update_var
		call export_data(30)
	end subroutine
	subroutine export_data(ut)
		integer :: ut,l1,i
		do l1=2,size(var(1:))
			do i=1,size(var(l1)%bd)
				write(ut,"(es17.9$)")latt%bond(var(l1)%nb)%bd(i)%r,&
					latt%bond(var(l1)%nb)%bd(i)%dr,&
					var(l1)%val(var(l1)%i2v(i)),&
					var(l1)%bd(i)
				write(ut,"(x)")
			enddo
			write(ut,"(x/)")
		enddo
	end subroutine
	recursive function gt(arg)
		real(8) :: arg(:)
		real(8) :: gt
		real(8) :: cr,ci,dt,m
		real(8) :: X2,a,X,dp,arg_(size(arg))
		cr=arg(1)
		ci=arg(2)
		dt=arg(3)
		m=arg(4)
		cr=0d0
		ci=0d0
		dt=0d0
		dp=abs(1-nf)
		X2=2d0*(dt**2+cr**2+ci**2)
		X=2d0*dp**2*(dt**2-cr**2-ci**2)+m**2*X2+X2**2/2d0
		a=1d0+4d0*X/(1d0-dp**2+m**2)**2
		gt=2d0*dp*(1d0-dp)/(1d0-dp**2+m**2)*((1d0+dp)**2-2d0*X2-m**2)/((1d0+dp)**2-m**2)*a
	end function
	recursive function gxy(arg)
		real(8) :: arg(:)
		real(8) :: gxy
		real(8) :: cr,ci,dt,m
		real(8) :: X2,a,X,dp,arg_(size(arg))
		cr=arg(1)
		ci=arg(2)
		dt=arg(3)
		m=arg(4)
		!cr=0d0
		!ci=0d0
		!dt=0d0
		dp=abs(1-nf)
		X2=2d0*(dt**2+cr**2+ci**2)
		X=2d0*dp**2*(dt**2-cr**2-ci**2)+m**2*X2+X2**2/2d0
		a=1d0+4d0*X/(1d0-dp**2+m**2)**2
		gxy=(2d0*(1d0-dp)/(1d0-dp**2+m**2))**2/a**7
	end function
	recursive function gz(arg)
		real(8) :: arg(:)
		real(8) :: gz
		real(8) :: cr,ci,dt,m
		real(8) :: X2,a,X,dp,arg_(size(arg))
		cr=arg(1)
		ci=arg(2)
		dt=arg(3)
		m=arg(4)
		!cr=0d0
		!ci=0d0
		!dt=0d0
		dp=abs(1-nf)
		X2=2d0*(dt**2+cr**2+ci**2)
		X=2d0*dp**2*(dt**2-cr**2-ci**2)+m**2*X2+X2**2/2d0
		a=1d0+4d0*X/(1d0-dp**2+m**2)**2
		if(abs(m**2+X2)<1d-7) then
			gz=(2d0*(1d0-dp)/(1d0-dp**2+m**2))**2/a**7
		else
			gz=(2d0*(1d0-dp)/(1d0-dp**2+m**2))**2/a**7/(m**2+X2)*(X2+m**2*(1d0+6d0*X2*(1d0-dp)**2/(1d0-dp**2+m**2)/a**3)**2)
		endif
	end function
	subroutine self_consist(fE)
		real(8), optional :: fE
		real(8) :: x(sum(var(2:)%n)),minf
		real(8) :: lb(size(x)),ub(size(x)),cst(1)
		integer(8) opt
		integer ires

		lb(1:)=1d-6
		!lb(1:)=0d0
		ub(1:)=1.5d0
		!lb(1)=1d-6
		!x=abs(var(2:)%put())
		x=var(2:)%put()
		!call myfunc(minf,size(x),x)
		x=1d-1
		!call random_number(x(1:2))
		!x(1)=0.1628d0
		!x(2)=0.1628d0
		!write(*,"(es12.4$)")minf
		!write(*,"(x)")
		!write(*,"(es12.4$)")var(1:)%val(1)
		!write(*,"(x)")
		!write(*,"(es12.4$)")var(1)%val(1),x
		!write(*,"(x)")
		!stop

		call nlo_create(opt, NLOPT_LN_BOBYQA, size(x))
		call nlo_set_lower_bounds(ires, opt, lb)
		call nlo_set_upper_bounds(ires, opt, ub)
		call nlo_set_min_objective(ires, opt, myfunc, cst)
		call nlo_set_xtol_abs(ires, opt, 1d-5)
		!call nlo_set_xtol_rel(ires, opt, 1d-5)
		call nlo_set_maxeval(ires, opt, size(x)*200)

		call nlo_optimize(ires, opt, x, minf)
		call nlo_destroy(opt)
		if(present(fE)) fE=minf
		!call var(2:)%get(x)
		!write(*,"(es12.4$)")x
		!write(*,"(x)")
		!call var%update()
		!write(*,"(es12.4$)")var(2:)%val(1)
		!write(*,"(x)")
		!write(*,"(A$)")"------------------->"
		!write(*,"(es12.4$)")var(1)%val(1),x
		!write(*,"(x)")
		call update_var(x)
	end subroutine
	subroutine myfunc(val, n, x, grad, need_gradient, f_data)
		integer :: n
		real(8) :: val, x(n)
		real(8), optional :: grad(n),f_data(1)
		integer, optional :: need_gradient

		complex(8) :: H(latt%Ns*spin,latt%Ns*spin)
		real(8) :: E(size(H,1)),vr(n),tmp
		complex(8) :: D(size(H,1),1,n),cH(size(H,1),size(H,2))
		real(8) :: f(size(E)),err,dE(size(H,1),n)
		integer :: l,lp,l1,l2,l3,m,i,j,k,info
		real(8) :: wa(16),mu(1)

		real(8) :: arg(4,size(t))
		call var(2:)%get(x)
		!call var%update()
		info=1
		mu(1)=0d0
		call hybrd1(MF_var,1,mu(1),vr(1:1),1d-7,info,wa(:16),16)
		call var(1:)%get(mu)
		vr=0d0
		val=0d0
		!$OMP PARALLEL DO REDUCTION(+:val,vr) PRIVATE(H,E,cH,f,dE,D)
		do k=1,size(brizon%k,1)
			call var%Hamilton(H,brizon%k(k,:))
			call mat_diag(H,E)
			cH=transpose(conjg(H))
			f=1d0/(exp(E/Tk)+1d0)
			call var(2:)%dHamilton(H,cH,D(:,1:1,:),brizon%k(k,:))
			dE=real(D(:,1,:))
			do l=1,size(E)
				val=val-(f(l)*log(max(f(l),1d-17))+(1d0-f(l))*log(max(1d0-f(l),1d-17)))
			enddo
			do l=1,n
				vr(l)=vr(l)+sum(dE(:,l)*f(:))
			enddo
		enddo
		!$OMP END PARALLEL DO
		val=val/(size(brizon%k,1)*latt%Ns)
		vr=vr/size(brizon%k,1)

		l=0
		do l1=2,ubound(var,1)
			do l2=1,size(var(l1)%val)
				l=l+1
				if(var(l1)%sg==3) then
					vr(l)=vr(l)*var(l1)%bd(1)
				endif
				tmp=0d0
				if(is_sc.and.var(l1)%nb==0) then
					do l3=1,size(var(l1)%v2i(l2)%i)
						i=var(l1)%v2i(l2)%i(l3)
						tmp=tmp+real(var(l1)%bd(i))*var(l1)%Vbd(i)
					enddo
					if(mod(abs(var(l1)%sg),10)==4) tmp=-tmp
				endif
				if(var(l1)%nb==0) then
					var(l1)%val(l2)=(vr(l)+tmp)/(size(var(l1)%v2i(l2)%i))
				else
					var(l1)%val(l2)=(vr(l)+tmp)/(2*spin*size(var(l1)%v2i(l2)%i))
				endif
			enddo
		enddo

		val=Eavg()-Tk*val
		!write(*,"(es12.4$)")val,x
		!write(*,"(x)")
	end subroutine
	recursive function Eavg(l)
		integer, optional :: l
		real(8) :: Eavg
		real(8) :: arg(4,size(t))
		integer :: l1
		if(present(l)) then
			var(l)%val(1)=var(l)%val(1)+1d-6
			Eavg=Eavg()
			var(l)%val(1)=var(l)%val(1)-2d-6
			Eavg=Eavg-Eavg()
			Eavg=Eavg/2d-6
			var(l)%val(1)=var(l)%val(1)+1d-6
			Eavg=Eavg/8d0
		else
			Eavg=0d0
			arg=0d0
			do l1=lbound(var,1),ubound(var,1)
				select case(abs(var(l1)%sg))
				case(2)
					arg(3,1)=var(l1)%val(1)
				case(3)
					arg(1,var(l1)%nb)=var(l1)%val(1)
				case(13)
					arg(2,1)=var(l1)%val(1)
				case(4)
					arg(4,1)=var(l1)%val(1)
				end select
			enddo
			do l1=lbound(var,1),ubound(var,1)
				select case(abs(var(l1)%sg))
				case(2)
					Eavg=Eavg-4d0*(DJ*(gxy(arg(:,1))/2d0+gz(arg(:,1))/4d0)-V)*arg(3,1)**2
				case(3)
					Eavg=Eavg-8d0*t(var(l1)%nb)*gt(arg(:,var(l1)%nb))*arg(1,var(l1)%nb)
					if(var(l1)%nb==1) then
						Eavg=Eavg-4d0*(DJ*(gxy(arg(:,1))/2d0+gz(arg(:,1))/4d0)+V)*arg(1,1)**2
					endif
				case(13)
					Eavg=Eavg-4d0*(DJ*(gxy(arg(:,1))/2d0+gz(arg(:,1))/4d0)+V)*arg(2,1)**2
				case(4)
					Eavg=Eavg-2d0*arg(4,1)**2*gz(arg(:,1))/4d0*DJ
				end select
			enddo
		endif
	end function
	subroutine update_var(x)
		integer :: l1,l
		real(8) :: x(:),scal
		do l1=1,size(var)
			if(var(l1)%sg==3.and.var(l1)%nb==1) then
				scal=abs(Eavg(l1)/x(l1-1))
				exit
			endif
		enddo
		var(1)%val=var(1)%val*scal
		do l1=2,size(var)
			var(l1)%val=x(l1-1)*scal
		enddo
	end subroutine
	function findTc(l,is,Tm)
		integer :: l,is
		real(8) :: findTc,Tm
		real(8) :: dTk,er,pod,order
		integer :: isg
		er=1d-3
		if(is<0) then
			Tk=1d-5
		else
			Tk=Tm
		endif
		dTk=1d-3
		call self_consist()
		pod=sum(abs(var(l)%val(:)))/(size(var(l)%val(:)))-er
		do 
			Tk=Tk+dTk*sign(1d0,pod)*is
			!write(*,"(es12.4$)")Tk,dTk
			if(Tk<1d-4.or.Tk>Tm) then
				if(dTk>1d-4) then
					Tk=Tk-dTk*sign(1d0,pod)*is
					dTk=dTk*0.3333d0
					!write(*,"(1x)")
					cycle
				else
					findTc=abs(Tk)
					exit
				endif
			endif
			call self_consist()
			order=sum(abs(var(l)%val(:)))/(size(var(l)%val))
			call is_cross(pod,order-er,isg)
			!write(*,"(es12.4$)")order-er
			!write(*,"(i3)"),isg
			if(isg/=0) then
				if(abs(dTk)<1d-4) then
					findTc=Tk
					exit
				endif
				dTk=dTk*0.3333d0
			endif
		enddo
	end function
	subroutine raman(ut,gm,omg,m,peak)
		integer :: ut,m
		real(8) :: gm,omg(:)
		real(8), optional, allocatable :: peak(:)
		complex(8) :: H(latt%Ns*spin,latt%Ns*spin),Bg(size(H,1),size(H,2),2)
		complex(8) :: cH(size(H,1),size(H,1))
		real(8) :: E(size(H,1)),f(size(E)),domg,bt
		complex(8) :: R(m,2)
		integer :: i,l,n,m1,m2
		integer, allocatable :: ipeak(:)
		bt=1d0/Tk
		domg=(omg(2)-omg(1))/m
		R=0d0
		!$OMP PARALLEL DO REDUCTION(+:R) PRIVATE(E,H,cH,f,Bg)
		do n=1,size(brizon%k,1)
			call var%Hamilton(H,brizon%k(n,:))
			call mat_diag(H,E)
			cH=transpose(conjg(H))
			f=1d0/(1d0+exp(bt*E))

			call var(:)%Hamilton(Bg(:,:,1),brizon%k(n,:),b1g)
			call var(:)%Hamilton(Bg(:,:,2),brizon%k(n,:),b2g)
			do l=1,2
				Bg(:,:,l)=matmul(cH,matmul(Bg(:,:,l),H))
			enddo
			do m1=1,4
				do m2=1,4
					do l=1,m
						R(l,:)=R(l,:)+Bg(m1,m2,:)*Bg(m2,m1,:)*(f(m1)-f(m2))/(omg(1)+domg*l+E(m2)-E(m1)+img*gm)
					enddo
				enddo
			enddo
		enddo
		!$OMP END PARALLEL DO
		do i=1,m
			write(ut,"(es17.9$)")omg(1)+domg*i,R(i,:)/(size(brizon%k,1)*latt%Ns)
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
	integer :: l,m,i
	real(8) :: n(45)=(/(0.999d0-i/100d0,i=0,44)/)
	real(8) :: Ts(size(n),2),Td(size(n),2),Tc(2),fE,kf(3)
	real(8), allocatable :: peak(:)
	f=openfile(unit=10,file='../data/phase.dat')
	f=openfile(unit=20,file='../data/band.dat')
	f=openfile(unit=30,file='../data/order.dat')
	f=openfile(unit=40,file='../data/fermis.dat')
	f=openfile(unit=50,file='../data/DOS.dat')
	f=openfile(unit=60,file='../data/raman.dat')
	f=openfile(unit=70,file='../data/map_raman.dat')
	f=openfile(unit=80,file='../data/map_band.dat')
	f=openfile(unit=90,file='../data/energy.dat')
	f=openfile(unit=100,file='../data/gap.dat')
	f=openfile(unit=101,file='../data/lattice.dat')

	call initial()

	!Tk=1d-4
	!do l=1,size(n)
		!nf=n(l)
		!call self_consist(fE)
		!write(*,"(es12.4$)")nf,Tk,fE,var(2:)%val(1)
		!write(*,"(x)")
		!write(100,"(es12.4$)")nf,Tk,fE,var(2:)%val(1)
		!write(100,"(x)")
		!!call band(20,(/0d0,0d0,0d0/),(/pi,pi,0d0/),128)
		!!call band(20,(/pi,pi,0d0/),(/0d0,pi,0d0/),128)
		!!call band(20,(/0d0,pi,0d0/),(/0d0,0d0,0d0/),128)
	!enddo
	!stop
	!do l=1,10
		!nf=8.6d-1
		!!Tk=6.8d-02/10d0*l
		!Tk=5.5629d-02/10d0*l
		!call self_consist()
		!!call var%update()
		!write(*,"(es12.4$)")nf,Tk,var(1:)%val(1)
		!write(*,"(x)")
		!write(100,"(es12.4$)")nf,Tk,var(1:)%val(1)
		!write(100,"(x)")
		!call raman(60,0.05d0,(/0d0,0.8d0/),256,peak)
		!!call DOS(50,0.01d0,(/-1.8d0,1.8d0/),256)
		!!call fermis(40,0.01d0,o_brizon%k,0d0)
		!!stop
		!!call band(20,(/0d0,0d0,0d0/),(/pi,pi,0d0/),128)
		!!call band(20,(/pi,pi,0d0/),(/0d0,pi,0d0/),128)
		!!call band(20,(/0d0,pi,0d0/),(/0d0,0d0,0d0/),128)
		!!kf=(/pi,4.663301595d-01,0d0/)
		!!call band(20,kf,kf,1)
	!enddo
	!stop

	Tc=2d-1
	do l=9,size(n)-14
		nf=n(l)
		Tc(1)=findTc(2,1,Tc(1)+0.02)
		!Tc(1)=findTc(2,1,0.12d0)
		Tc(2)=findTc(2,-1,Tc(1))
		write(10,"(es12.4$)")nf,Tc
		write(*,"(es12.4$)")nf,Tc
		write(10,"(x)")
		write(*,"(x)")
	enddo
	write(10,"(x/)")

	Tc=2d-1
	do l=9,size(n)-14
		nf=n(l)
		Tc(1)=findTc(3,1,Tc(1)+0.02)
		!Tc(1)=findTc(3,1,0.12d0)
		!if(Tc(1)<1d-4) then
			!Td(l:,:)=0d0
			!exit
		!endif
		Tc(2)=findTc(3,-1,Tc(1))
		write(10,"(es12.4$)")nf,Tc
		write(*,"(es12.4$)")nf,Tc
		write(10,"(x)")
		write(*,"(x)")
	enddo
	write(10,"(x/)")
end program
