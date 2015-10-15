module pmt
	use M_const
	implicit none
	real(8), parameter :: t(1)=(/1d0/),&
		V=0.00d0,DJ=0.30d0
	integer :: cs=2
end module
module selfcons
	use pmt
	use M_Hamilton
	use M_utility
	implicit none
contains
	subroutine initial()
		integer :: l2,l3,l,i
		real(8), allocatable :: bd0(:),bd1(:)
		allocate(var(-10000:10000))
		call init_random_seed()

		! lattice 
		latt%a1=(/1d0,1d0,0d0/)
		latt%a2=(/-1d0,1d0,0d0/)
		!latt%T1=(/1d0,0d0,0d0/)*2
		!latt%T2=(/0d0,1d0,0d0/)*8
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
		brizon%n1=256
		brizon%n2=256
		call latt%gen_brizon(brizon)
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

		!! ddw
		!call gen_var(sg=13,nb=1)
		!do i=1,size(var(iv(0))%bd)
			!var(iv(0))%bd(i)=img*ab(latt%bond(var(iv(0))%nb)%bd(i)%i(1))*dwave(i)
		!enddo
		!do l2=1,size(var(iv(0))%val)
			!var(iv(0))%val(l2)=1d-1
		!enddo

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
			var(iv(0))%bd=1d0
			do l2=1,size(var(iv(0))%val)
				var(iv(0))%val(l2)=t(l)
			enddo
		enddo


		call var_shrink()
		var(lbound(var,1))%update => update_var
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
	subroutine update_var()
		integer :: l1,l
		real(8) :: arg(4,size(t))
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
				l=3
				var(l1)%val=-t(1)*gt(arg(:,1),l)*arg(1,1)&
					-(DJ*(gxy(arg(:,1))/2d0+gz(arg(:,1))/4d0)-V)*arg(l,1)&
					-0.5d0*DJ*(gxy(arg(:,1),l)/2d0+gz(arg(:,1),l)/4d0)*sum(arg(1:3,1)**2)&
					-0.25d0*arg(4,1)**2*gz(arg(:,1),l)/4d0*DJ
			case(3)
				l=1
				var(l1)%val=-t(var(l1)%nb)*gt(arg(:,var(l1)%nb))-t(var(l1)%nb)*gt(arg(:,var(l1)%nb),l)*arg(1,var(l1)%nb)
				if(var(l1)%nb==1) then
					var(l1)%val(:)=var(l1)%val(:)&
						-((gxy(arg(:,1))/2d0+gz(arg(:,1))/4d0)*DJ+V)*arg(l,1)&
						-0.5d0*DJ*(gxy(arg(:,1),l)/2d0+gz(arg(:,1),l)/4d0)*sum(arg(1:3,1)**2)&
						-0.25d0*arg(4,1)**2/4d0*gz(arg(:,1),l)*DJ
				endif
			case(13)
				l=2
				var(l1)%val(:)=-t(var(l1)%nb)*gt(arg(:,var(l1)%nb),l)*arg(1,var(l1)%nb)&
					-((gxy(arg(:,1))/2d0+gz(arg(:,1))/4d0)*DJ+V)*arg(l,1)&
					-0.5d0*DJ*(gxy(arg(:,1),l)/2d0+gz(arg(:,1),l)/4d0)*sum(arg(1:3,1)**2)&
					-0.25d0*arg(4,1)**2/4d0*gz(arg(:,1),l)*DJ
			case(4)
				l=4
				var(l1)%val=-8d0*t(var(l1)%nb)*gt(arg(:,1),l)*arg(1,1)&
					-4d0*DJ*(gxy(arg(:,1),l)/2d0+gz(arg(:,1),l)/4d0)*sum(arg(1:3,1)**2)&
					-4d0*arg(4,1)*gz(arg(:,1))/4d0*DJ&
					-2d0*arg(4,1)**2*gz(arg(:,1),l)/4d0*DJ
			end select
		enddo
		!write(*,"(es14.2$)")arg(:,1),var(:)%val(1)
		!write(*,"(x)")
	end subroutine
	recursive function gt(arg,i)
		real(8) :: arg(:)
		integer, optional :: i
		real(8) :: gt
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
		if(present(i)) then
			arg_=arg
			arg_(i)=arg(i)+1d-6
			gt=gt(arg_)
			arg_(i)=arg(i)-1d-6
			gt=gt-gt(arg_)
			gt=gt/2d-6
		else
			X2=2d0*(dt**2+cr**2+ci**2)
			X=2d0*dp**2*(dt**2-cr**2-ci**2)+m**2*X2+X2**2/2d0
			a=1d0+4d0*X/(1d0-dp**2+m**2)**2
			gt=2d0*dp*(1d0-dp)/(1d0-dp**2+m**2)*((1d0+dp)**2-2d0*X2-m**2)/((1d0+dp)**2-m**2)*a
		endif
	end function
	recursive function gxy(arg,i)
		real(8) :: arg(:)
		integer, optional :: i
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
		if(present(i)) then
			arg_=arg
			arg_(i)=arg(i)+1d-6
			gxy=gxy(arg_)
			arg_(i)=arg(i)-1d-6
			gxy=gxy-gxy(arg_)
			gxy=gxy/2d-6
		else
			X2=2d0*(dt**2+cr**2+ci**2)
			X=2d0*dp**2*(dt**2-cr**2-ci**2)+m**2*X2+X2**2/2d0
			a=1d0+4d0*X/(1d0-dp**2+m**2)**2
			gxy=(2d0*(1d0-dp)/(1d0-dp**2+m**2))**2/a**7
		endif
	end function
	recursive function gz(arg,i)
		real(8) :: arg(:)
		integer, optional :: i
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
		if(present(i)) then
			arg_=arg
			arg_(i)=arg(i)+1d-6
			gz=gz(arg_)
			arg_(i)=arg(i)-1d-6
			gz=gz-gz(arg_)
			gz=gz/2d-6
		else
			X2=2d0*(dt**2+cr**2+ci**2)
			X=2d0*dp**2*(dt**2-cr**2-ci**2)+m**2*X2+X2**2/2d0
			a=1d0+4d0*X/(1d0-dp**2+m**2)**2
			if(abs(m**2+X2)<1d-7) then
				gz=(2d0*(1d0-dp)/(1d0-dp**2+m**2))**2/a**7
			else
				gz=(2d0*(1d0-dp)/(1d0-dp**2+m**2))**2/a**7/(m**2+X2)*(X2+m**2*(1d0+6d0*X2*(1d0-dp)**2/(1d0-dp**2+m**2)/a**3)**2)
			endif
		endif
	end function
	subroutine self_consist()
		integer :: info,i,l
		real(8) :: x(sum(var(1:)%n)),v(size(x)),wa(nint((size(x)*(3*size(x)+13))/2.)+10),fE,mfE,mx(size(x)),px(size(x)),err,merr,s
		common fE,err
		info=1
		
		select case(cs)
		case(1)
			mfE=0d0
			x=0d0
			do i=1,20
				call hybrd1(MF_var,size(x(:)),x(:),v(:),1d-7,info,wa,size(wa))
				call var(1:)%get(x)
				fE=free_energy()
				if(fE<mfE) then
					mfE=fE
					mx=x
					merr=err
				endif
				write(*,"(i5$)")i
				write(*,"(es12.4$)")fE,err,var(1:)%val(1)
				write(*,"(x)")
				call export_data(30)
				do l=1,ubound(var,1)
					call random_number(var(l)%val(:))
					var(l)%val(:)=(var(l)%val(:)-0.5d0)*0.2d0
				enddo
				x=var(1:)%put()
			enddo
			call var(1:)%get(mx)
			write(*,"(es12.4$)")mfE,merr,var(1:)%val(1)
			write(*,"(x)")
		case(2)
			!do l=1,ubound(var,1)
			do l=2,3
				call random_number(var(l)%val(:))
				var(l)%val(:)=(var(l)%val(:)+0.1d0)*0.2d0
			enddo
			x=var(1:)%put()
			i=0
			do 
				call var%update()
				i=i+1
				!write(*,"(i5$)")i
				!write(*,"(es12.4$)")v,x
				!write(*,"(x)")
				if(var(1)%sg==1) then
					call hybrd1(MF_var,size(x(1:1)),x(1:1),v(1:1),1d-7,info,wa(:16),16)
					if(info/=1) then
						write(*,*)"hybrd1 err, info=",info
						stop
					endif
				endif
				call MF_var(size(x),x,v,info=-1)
				if(err<1d-4) then
					call hybrd1(MF_var,size(x(:)),x(:),v(:),1d-7,info,wa,size(wa))
					call var(1:)%get(x)
					!write(*,"(i5$)")i
					!write(*,"(es12.4$)")fE,err,var(1:)%val(1)
					!write(*,"(x)")
					exit
				endif
				call var(1:)%get(x)
			enddo
		end select
	end subroutine
	function free_energy()
		integer :: l1
		real(8) :: arg(4,size(t)),free_energy
		free_energy=0d0
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
				free_energy=free_energy-4d0*(DJ*(gxy(arg(:,1))/2d0+gz(arg(:,1))/4d0)-V)*arg(3,1)**2
			case(3)
				free_energy=free_energy-8d0*t(var(l1)%nb)*gt(arg(:,var(l1)%nb))*arg(1,var(l1)%nb)
				if(var(l1)%nb==1) then
					free_energy=free_energy-4d0*(DJ*(gxy(arg(:,1))/2d0+gz(arg(:,1))/4d0)+V)*arg(1,1)**2
				endif
			case(13)
				free_energy=free_energy-4d0*(DJ*(gxy(arg(:,1))/2d0+gz(arg(:,1))/4d0)+V)*arg(2,1)**2
			case(4)
				free_energy=free_energy-2d0*arg(4,1)**2*gz(arg(:,1))/4d0*DJ
			end select
		enddo
	end function
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
	real(8) :: Ts(size(n),2),Td(size(n),2),Tc(2)
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

	Tk=1d-4
	do l=1,size(n)
		nf=n(l)
		call self_consist()
		write(*,"(es12.4$)")nf,Tk,var(1:)%val(1)
		write(*,"(x)")
		write(100,"(es12.4$)")nf,Tk,var(1:)%val(1)
		write(100,"(x)")
	enddo
	stop
	do l=1,10
		nf=8.00d-01
		Tk=8.0d-02/10d0*l
		call self_consist()
		var(2)%val=0d0
		call var%update()
		write(*,"(es12.4$)")nf,Tk,var(1:)%val(1)
		write(*,"(x)")
		write(100,"(es12.4$)")nf,Tk,var(1:)%val(1)
		write(100,"(x)")
		call raman(60,0.01d0,(/0d0,0.8d0/),256,peak)
		call DOS(50,0.01d0,(/-1.8d0,1.8d0/),256)
		stop
	enddo
	stop

	Tc=1d-1
	do l=10,size(n)-15
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

	Tc=1d-1
	do l=10,size(n)-15
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
