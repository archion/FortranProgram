module pmt
	use M_const
	implicit none
	real(8), parameter :: t(3)=(/1d0,-0.3d0,0.1d0/),&
		!V=0.0325d0,DJ=0.35d0
		V=0.0325d0,DJ=0.35d0
		!V=-0.25d0/4d0,DJ=0.25d0
		!V=0d0,DJ=0.25d0
end module
module selfcons
	use pmt
	use M_Hamilton_m
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
		brizon%n1=128
		brizon%n2=brizon%n1
		call latt%gen_brizon(brizon)
		!call check_lattice(101)
		write(*,*)"Total site number is: ",latt%Ns

		allocate(bd0(latt%Ns))
		do i=1,size(bd0)
			!bd0(i)=latt%bond(0)%bd(i)%i(2)
			bd0(i)=1
		enddo

		allocate(bd1(size(latt%sb(1)%nb(1)%bd)))
		do i=1,size(bd1)
			!bd1(i)=latt%bond(1)%bd(i)%i(2)
			bd1(i)=1
		enddo

		! cp
		call gen_var(sg=1,nb=0)
		var(iv(0))%val=0d0
		var(iv(0))%bd=-1d0

		! d-wave sc
		call gen_var(sg=2,nb=1,val=bd1,V=(-DJ*0.75d0+V))
		do i=1,size(var(iv(0))%bd)
			var(iv(0))%bd(i)=dwave(i)
		enddo
		var(iv(0))%val=1d-1

		! ddw
		call gen_var(sg=3,nb=1,V=(-DJ*0.75d0-V))
		do i=1,size(var(iv(0))%bd)
			var(iv(0))%bd(i)=img*ab(latt%sb(1)%nb(var(iv(0))%nb)%bd(i)%i(1))*dwave(i)
		enddo
		var(iv(0))%val=1d-1

		!! sdw
		!call gen_var(sg=4,nb=0,val=bd0,V=DJ,Vnb=1)
		!do i=1,size(var(iv(0))%bd)
			!var(iv(0))%bd(i)=ab(i)
		!enddo
		!var(iv(0))%val=1d-1

		! bond order
		call gen_var(sg=3,nb=1,val=bd1,V=(-DJ*0.75d0-V))
		var(iv(0))%bd=1d0
		var(iv(0))%val=0d0

		! hp
		do l=1,size(t)
			call gen_var(sg=-3,nb=l)
			var(iv(0))%bd=-t(l)
			var(iv(0))%val=1d0
		enddo

		call var_shrink()
		!call export_data(30)
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
	subroutine update()
		integer :: l
		do l=lbound(var,1),0
			if(var(l)%sg==-3) then
				var(l)%bd(:)=-t(var(l)%nb)*abs(1d0-nf)
			endif
		enddo
	end subroutine
	!subroutine self_consist(fE_)
		!real(8), optional :: fE_
		!real(8) :: x(sum(var(1:)%n)),v(size(x)),wa(nint((size(x)*(3*size(x)+13))/2.)+10),minf
		!integer :: info
		!real(8) :: fE
		!common fE
		!call update()
		!info=1
		!var(2)%val=1d-1
		!var(3)%val=1d-1
		!x=put(var(1:))
		!x=x+1d-5
		!call hybrd1(minpack_fn,size(x),x,v,1d-7,info,wa,size(wa))
		!write(*,"(i3$)")info
		!write(*,"(es12.4$)")Tk,fE,var(1:)%val(1),sum(abs(v))
		!write(*,"(x)")
	!end subroutine
	subroutine checkorder(l,rg)
		integer :: l
		real(8) :: rg(:),fE
		real(8) :: x(sum(var(2:)%n)),grad(size(x)),f_data(1)
		call self_consist()
		x=put(var(2:))
		x(l-1)=rg(1)
		do while(x(l-1)<=rg(2))
			call nlopt_fn(fE, size(x), x, grad, 0, f_data)
			x(l-1)=x(l-1)+rg(3)
			write(*,"(es12.4$)")x(l-1)
			write(*,"(es17.9$)")fE
			write(*,"(x)")
		enddo
	end subroutine
	subroutine self_consist_nlopt(fE)
		real(8), optional :: fE
		real(8) :: x(sum(var(2:)%n)),minf
		real(8) :: lb(size(x)),ub(size(x)),cst(1)
		integer(8) opt
		integer ires
		call update()

		lb=1d-5
		ub=1d0
		var(2)%val=1d-1
		var(3)%val=1d-1
		x=put(var(2:))
		x=1d-1

		!call nlo_create(opt, NLOPT_LN_BOBYQA, size(x))
		!call nlo_create(opt, NLOPT_LN_SBPLX, size(x))

		call nlo_create(opt, NLOPT_LD_LBFGS, size(x))
		!call nlo_create(opt, NLOPT_LD_TNEWTON_PRECOND_RESTART, size(x))

		call nlo_set_lower_bounds(ires, opt, lb)
		call nlo_set_upper_bounds(ires, opt, ub)
		call nlo_set_min_objective(ires, opt, nlopt_fn, cst)
		call nlo_set_xtol_rel(ires, opt, 1d-7)
		!call nlo_set_maxeval(ires, opt, size(x)*200)

		call nlo_optimize(ires, opt, x, minf)
		if(ires<0) then
			write(*,"(A$)")"ires is negetive"
			write(*,"(i3$)")ires
			write(*,"(es12.4)")cst
			!if(cst(1)>1d-5) then
				!stop
			!endif
			!x=put(var(2:))
			!call nlo_optimize(ires, opt, x, minf)
			!if(cst(1)>1d-4) then
				!stop
			!endif
		endif
		call nlo_destroy(opt)
		if(present(fE)) fE=minf
		write(*,"(i3$)")ires
		write(*,"(es12.4$)")nf,Tk,minf,x,cst
		write(*,"(x)")
	end subroutine
	subroutine nlopt_fn(val, n, x, grad, need_gradient, f_data)
		integer :: n
		real(8) :: val, x(n)
		real(8) :: grad(n),f_data(1)
		integer :: need_gradient
		integer :: info,i,j,t
		real(8) :: mu(1),vmu(1),wa(16),v(n),val_
		call get(var(2:),x)
		mu(1)=var(1)%val(1)-1d-1
		t=0
		do
			info=1
			call hybrd1(minpack_fn,1,mu(1:1),vmu(1:1),1d-7,info,wa(:16),16)
			if(info/=1.or.abs(vmu(1))>1d-6) then
				t=t+1
				call random_number(mu)
				mu=(mu-0.5d0)
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
		call MF_val(2,v,val,0)
		!write(*,"(es12.4$)")val,v
		!call get(var(2:),x)
		!call MF_val(2,v,val,1)
		!i=2
		!j=1
		!write(*,"(es12.4$)")val,v(sum(var(2:i-1)%n)+j)
		!call get(var(2:),x)
		!var(i)%val(j)=var(i)%val(j)-1d-7
		!call MF_val(2,v,val,1)
		!call get(var(2:),x)
		!var(i)%val(j)=var(i)%val(j)+1d-7
		!call MF_val(2,v,val_,1)
		!write(*,"(es12.4$)")(val_-val)/2d-7
		!write(*,"(x)")
		!stop
		!val=sum(abs(v))
		!write(*,"(es12.4$)")val
		if(need_gradient/=0) then
			grad=v
			!write(*,"(es12.4$)")grad
		endif
		f_data(1)=sum(abs(v))/size(v)
		!write(*,"(es12.4$)")f_data(1)
		!write(*,"(es12.4$)")v
		!write(*,"(x)")
		!stop
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
			!call Hamilton(var,Bg(:,:,1),brizon%k(n,:),b1g)
			!call Hamilton(var,Bg(:,:,2),brizon%k(n,:),b2g)
			call Hamilton(var(:0),Bg(:,:,1),brizon%k(n,:),b1g)
			call Hamilton(var(:0),Bg(:,:,2),brizon%k(n,:),b2g)
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
		complex(8) :: R(m,2)
		integer :: i,l,n,m1,m2
		integer, allocatable :: ipeak(:)
		bt=1d0/Tk
		domg=(omg(2)-omg(1))/m
		R=0d0
		!$OMP PARALLEL DO REDUCTION(+:R) PRIVATE(E,H,cH,f,Bg)
		do n=1,size(brizon%k,1)
			call Hamilton(var,H,brizon%k(n,:))
			call mat_diag(H,E)
			cH=transpose(conjg(H))
			f=1d0/(1d0+exp(bt*E))

			call Hamilton(var(:0),Bg(:,:,1),brizon%k(n,:),b1g)
			call Hamilton(var(:0),Bg(:,:,2),brizon%k(n,:),b2g)
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
	function superfluid()
		complex(8) :: H(latt%Ns*spin,latt%Ns*spin),vv(size(H,1),size(H,2)),mm(size(H,1),size(H,2)),cH(size(H,1),size(H,1))
		real(8) :: E(size(H,1)),f(size(E)),df(size(E)),k(3),bt,superfluid
		integer :: n,m1,m2
		superfluid=0d0
		bt=1d0/Tk
		!$OMP PARALLEL DO REDUCTION(+:superfluid) PRIVATE(k,mm,vv,E,H,cH,f,df)
		do n=1,size(brizon%k,1)
			k=brizon%k(n,:)
			call Hamilton(var,H,k)
			call mat_diag(H,E)
			cH=transpose(conjg(H))
			f=1d0/(1d0+exp(bt*E))
			df=-bt/((1d0+cosh(bt*E))*2)

			call Hamilton(var,vv,k,V)
			call Hamilton(var,mm,k,M)
			vv(:,:)=matmul(cH,matmul(vv,H))
			vv(:,:)=vv*transpose(vv)
			mm(:,:)=matmul(cH,matmul(mm,H))

			do m1=1,4
				do m2=1,4
					if(m1==m2.or.abs(E(m1)-E(m2))<1d-20) then
						cycle
					endif
					superfluid=superfluid+vv(m1,m2)*(f(m1)-f(m2))/(E(m1)-E(m2))
				enddo
				superfluid=superfluid+vv(m1,m1)*df(m1)+mm(m1,m1)*f(m1)
			enddo
		enddo
		!$END OMP PARALLEL DO
	contains
		function V(dr)
			real(8) :: dr(3)
			complex(8) :: V
			V=(-img*dr(1))
		end function
		function M(dr)
			real(8) :: dr(3)
			complex(8) :: M
			M=(-dr(1)**2)
		end function
	end function
end module
program main
	use selfcons
	implicit none
	logical :: f,flag=.true.
	integer :: l,m,i,nopt
	real(8) :: n(0:80)=(/(min(1d0-0.005d0*i,0.999d0),i=0,80)/)
	real(8) :: Ts(size(n),2),Td(size(n),2),Tc(2),Tnc,dnf,pnf,pTk(2)
	real(8), allocatable :: peak(:)
	f=openfile(unit=10,file='../data/phase.dat')
	!f=openfile(unit=20,file='../data/band.dat')
	!f=openfile(unit=30,file='../data/order.dat')
	!f=openfile(unit=40,file='../data/fermis.dat')
	!f=openfile(unit=50,file='../data/DOS.dat')
	!f=openfile(unit=60,file='../data/raman.dat')
	!f=openfile(unit=70,file='../data/map_raman.dat')
	!f=openfile(unit=80,file='../data/map_band.dat')
	!f=openfile(unit=90,file='../data/energy.dat')
	!f=openfile(unit=100,file='../data/gap.dat')
	f=openfile(unit=110,file='../data/bb.dat')
	!f=openfile(unit=101,file='../data/lattice.dat')

	call initial()

	!Tk=1d-3
	!nf=1d0-0.12d0
	!call self_consist()
	!read(*,*)
	!call checkorder(2,(/0d0,0.05d0,0.0005d0/))
	!stop
	!do l=1,size(n)
		!nf=n(l)
		!call self_consist()
		!write(*,"(es12.4$)")nf,Tk,var(1:)%val(1)
		!write(*,"(x)")
		!write(100,"(es12.4$)")nf,Tk,var(1:)%val(1)
		!write(100,"(x)")
		!stop
	!enddo
	!stop

	!do l=1,100
		!nf=8.6648d-1
		!Tk=4.5d-2/100d0*l
		!call self_consist()
		!write(*,"(es12.4$)")nf,Tk,var(2:3)%val(1),var(2:3)%val(1)*(/8d0,8d0/)*abs(var(2:3)%V)
		!write(*,"(x)")
		!write(100,"(es12.4$)")nf,Tk,var(2:3)%val(1),var(2:3)%val(1)*(/8d0,8d0/)*abs(var(2:3)%V)
		!write(100,"(x)")
		!!call raman(60,0.04d0,(/0d0,0.3d0/),256,peak)
	!enddo
	!stop
	!call latt%gen_origin_brizon((/1d0,0d0,0d0/),(/0d0,1d0,0d0/),o_brizon)
	!write(*,"(A$)")"raman peak"
	!write(*,"(es12.4$)")peak
	!write(*,"(x)")
	!do l=1,size(peak)
		!call raman_k(70,0.003d0,o_brizon%k,peak(l))
	!enddo

	!!call band_e(80,0.01d0,(/0d0,0d0,0d0/),(/pi,pi,0d0/),128,(/-1d0,1d0/),512)
	!!call band_e(80,0.01d0,(/pi,pi,0d0/),(/0d0,pi,0d0/),128,(/-1d0,1d0/),512)
	!!call band_e(80,0.01d0,(/0d0,pi,0d0/),(/0d0,0d0,0d0/),128,(/-1d0,1d0/),512)
	!call band(20,(/0d0,0d0,0d0/),(/pi,pi,0d0/),128)
	!call band(20,(/pi,pi,0d0/),(/pi,0d0,0d0/),128)
	!call band(20,(/pi,0d0,0d0/),(/0d0,0d0,0d0/),128)
	!call band(20,(/pi/2d0,pi/2d0,0d0/),(/pi,0d0,0d0/),128)
	!call energy(90,o_brizon%k)
	!call DOS(50,0.005d0,(/-0.4d0,0.4d0/),512,peak)
	!write(*,"(A$)")"DOS peak"
	!write(*,"(es12.4$)")peak
	!write(*,"(x)")
	!do l=1,size(peak)
		!call fermis(40,0.005d0,o_brizon%k,peak(l))
	!enddo
	!var(2)%val(:)=0d0
	!call fermis(40,0.005d0,o_brizon%k,0d0)
	!stop


	Tc=(/1d-1,1d-4/)
	do l=0,ubound(n,1),1
		nf=n(l)
		Tk=Tc(1)
		Tc(1)=find_order(3,-1,Tk,(/1d-4,0.3d0,2d-3/),1d-3)
		Tk=Tc(2)
		Tc(2)=find_order(3,1,Tk,(/1d-4,Tc(1),2d-3/),1d-3)
		if((Tc(2)>2d-4.or.Tc(1)<2d-4).and.flag) then
			flag=.false.
			Tk=1d-4
			nf=find_order(3,1,nf,(/abs(nf),pnf,abs(pnf-nf)/20d0/),1d-3)
			write(110,"(es12.4$)")nf
			Tk=pTk(1)
			!Tnc=find_order(2,-1,Tk,(/1d-4,0.25d0,2d-3/),1d-3)
			Tk=find_order(3,-1,Tk,(/1d-4,0.25d0,1d-4/),1d-3)
			write(10,"(es12.4$)")1d0-nf,Tk,1d0-nf,3d-4
			write(*,"(es12.4$)")1d0-nf,Tk,1d0-nf,3d-4
			write(10,"(x)")
			write(*,"(x)")
			pnf=n(l)
			!do i=1,10
				!Tk=Tnc/10d0*i
				!call self_consist()
				!write(*,"(es12.4$)")nf,Tk,var(2:3)%val(1)*(/8d0,4d0/)*abs(var(2:3)%V)
				!write(*,"(x)")
				!write(100,"(es12.4$)")nf,Tk,var(2:3)%val(1)*(/8d0,4d0/)*abs(var(2:3)%V)
				!write(100,"(x)")
				!call raman(60,0.05d0,(/0d0,0.8d0/),256,peak)
			!enddo
		endif
		if(Tc(1)<2d-4) then
			nopt=l
			if(nf<pnf) call swap(nf,pnf)
			dnf=abs(nf-pnf)
			Tc=pTk
			do
				nf=nf-dnf/10d0
				Tk=Tc(1)
				Tc(1)=find_order(3,-1,Tk,(/Tc(2),0.25d0,1d-4/),1d-3)
				Tk=Tc(2)
				Tc(2)=find_order(3,1,Tk,(/1d-4,Tc(1),1d-4/),1d-3)
				write(10,"(es12.4$)")1d0-nf,Tc(1),1d0-nf,Tc(2)
				write(*,"(es12.4$)")1d0-nf,Tc(1),1d0-nf,Tc(2)
				write(10,"(x)")
				write(*,"(x)")
				if((nf<=pnf).or.(Tc(1)-Tc(2))<1d-4) then
					write(110,"(es12.4$)")nf
					exit
				endif
			enddo
			exit
		endif
		pnf=nf
		pTk=Tc
		write(10,"(es12.4$)")1d0-n(l),Tc(1),1d0-n(l),Tc(2)
		write(*,"(es12.4$)")1d0-n(l),Tc(1),1d0-n(l),Tc(2)
		write(10,"(x)")
		write(*,"(x)")
	enddo
	write(10,"(x/)")

	pnf=nf
	Tc=(/1d-1,1d-4/)
	do l=nopt-1,-1,-1
		Tk=Tc(1)
		Tc(1)=find_order(2,-1,Tk,(/1d-4,0.25d0,2d-3/),1d-3)
		Tk=Tc(2)
		Tc(2)=find_order(2,1,Tk,(/1d-4,Tc(1),2d-3/),1d-3)
		write(10,"(es12.4$)")1d0-nf,Tc(1)
		write(*,"(es12.4$)")1d0-nf,Tc(1)
		write(10,"(x)")
		write(*,"(x)")
		nf=n(l)
	enddo
	write(10,"(x/)")
	write(*,"(x/)")

	nf=pnf
	Tc=(/1d-1,1d-4/)
	do l=nopt,ubound(n,1)+1,1
		Tk=Tc(1)
		Tc(1)=find_order(2,-1,Tk,(/1d-4,0.25d0,2d-3/),1d-3)
		Tk=Tc(2)
		Tc(2)=find_order(2,1,Tk,(/1d-4,Tc(1),2d-3/),1d-3)
		write(10,"(es12.4$)")1d0-nf,Tc(1)
		write(*,"(es12.4$)")1d0-nf,Tc(1)
		write(10,"(x)")
		write(*,"(x)")
		nf=n(l)
	enddo
	write(10,"(x/)")
end program
