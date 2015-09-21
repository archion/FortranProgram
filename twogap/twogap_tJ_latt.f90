module pmt
	use M_const
	implicit none
	real(8), parameter :: t(3)=(/1d0,-0.25d0,0.1d0/),&
		V=0.12d0,DJ=0.35d0
		!V=0.d0,DJ=0.25d0
	integer :: ipg(2)
	real(8) :: Vs=DJ-V,Vd=0.5d0*DJ+V,Tk=0.05,nf=0.8d0
end module
module selfcons
	use pmt
	use M_Hamilton
	use M_utility
	implicit none
contains
	subroutine initial()
		integer :: n,i,j,l,k
		type(t_var), allocatable :: tmp(:)
		allocate(tmp(-10000:10000))
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
		call latt%gen_neb(3)
		call latt%gen_bond(3)
		brizon%n1=512
		brizon%n2=brizon%n1
		call latt%gen_brizon(brizon)
		!call check_lattice(30)
		write(*,*)"Total site number is: ",latt%Ns

		! cp
		call gen_var(n,sg=1,nb=0,V=1d0,var=tmp)
		tmp(iv(0))%val=0d0
		tmp(iv(0))%bd_sg=-1d0

		! ddw
		call gen_var(n,sg=3,nb=1,V=Vd,var=tmp)
		ipg=(/3,1/)
		do i=iv(0)-n+sign(1,n),iv(0),sign(1,n)
			tmp(i)%val=1d-1
			do k=1,size(tmp(i)%n)
				tmp(i)%bd_sg(k)=ab(latt%bond(tmp(i)%nb)%bd(tmp(i)%n(k))%i(1))*img*dwave(tmp(i)%n(k))
			enddo
		enddo

		! d-wave sc
		call gen_var(n,sg=2,nb=1,V=Vs,var=tmp)
		do i=iv(0)-n+sign(1,n),iv(0),sign(1,n)
			tmp(i)%val=0d0
			do k=1,size(tmp(i)%n)
				tmp(i)%bd_sg(k)=dwave(tmp(i)%n(k))
			enddo
		enddo

		!! sdw
		!call gen_var(n,sg=4,nb=0,V=DJ*2d0,var=tmp)
		!ipg=(/4,0/)
		!do i=iv(0)-n+sign(1,n),iv(0),sign(1,n)
			!tmp(i)%val=1d-1
			!do k=1,size(tmp(i)%n)
				!tmp(i)%bd_sg(k)=ab(tmp(i)%n(k))
			!enddo
		!enddo

		! bond order
		call gen_var(n,sg=3,nb=1,V=Vd,var=tmp)
		do i=iv(0)-n+sign(1,n),iv(0),sign(1,n)
			tmp(i)%val=1d-1
			do k=1,size(tmp(i)%n)
				tmp(i)%bd_sg(k)=1d0
			enddo
		enddo

		! hp
		do l=1,size(t)
			call gen_var(n,sg=-3,nb=l,V=1d0,var=tmp)
			do i=iv(0)-n+sign(1,n),iv(0),sign(1,n)
				tmp(i)%bd_sg=-1d0
				tmp(i)%val=t(l)*(1d0-nf)
			enddo
		enddo


		allocate(var(iv(-1):iv(1)))
		do l=iv(-1),iv(1)
			call move_alloc(tmp(l)%bd_sg,var(l)%bd_sg)
			call move_alloc(tmp(l)%n,var(l)%n)
			var(l)%val=tmp(l)%val
			var(l)%nb=tmp(l)%nb
			var(l)%V=tmp(l)%V
			var(l)%sg=tmp(l)%sg      
		enddo
		deallocate(tmp)
	end subroutine
	subroutine update_var()
		integer :: getn(2),l
		do l=1,size(t)
			getn=get_var(-3,l,1)
			if(nf<0.8d0) then
				var(getn(2):getn(1))%val=t(l)
			else
				var(getn(2):getn(1))%val=t(l)*(1d0-nf)
			endif
		enddo
	end subroutine
	subroutine self_consist()
		integer :: info,i
		real(8) :: x(size(var(1:))),v(size(var(1:))),wa(nint((size(var(1:))*(3*size(var(1:))+13))/2.)+10),fE,mfE,mx(size(var(1:)))
		common fE
		call update_var()
		info=1

		mfE=0d0
		x=0d0
		do i=1,10
			call hybrd1(do_var,size(x(:)),x(:),v(:),1d-7,info,wa,size(wa))
			!write(*,"('->'es12.4$)")Tk,var(2:3)%val,fE
			if(fE<mfE) then
				mfE=fE
				mx=x
				!write(*,"(A$)")"*"
			endif
			call random_number(x(2:3))
			x(2:3)=x(2:3)+0.01d0
			!write(*,"(x)")
		enddo
		var(1:)%val=mx
		!write(*,"(x)")
		return

		i=0
		do 
			i=i+1
			if(var(1)%sg==1) then
				call hybrd1(do_var,size(x(1:1)),x(1:1),v(1:1),1d-7,info,wa(:16),16)
				x=var(1:)%val
			endif
			info=-1
			call do_var(size(x),x,v,info=info)
			write(*,"('->'es12.4$)")Tk,var(1:)%val
			if(sum(abs(var(1:)%val-x))/size(x)<1d-4) then
				info=1
				call hybrd1(do_var,size(x(:)),x(:),v(:),1d-7,info,wa,size(wa))
				var(1:)%val=x
				write(*,"('->'es12.4$)")Tk,var(2:3)%val
				write(*,*)i
				write(*,"(x)")
				exit
			endif
			write(*,*)i
			var(1:)%val=x
		enddo
	end subroutine
	subroutine do_var(n,x,v,info)
		integer, intent(in) :: n
		integer, intent(inout) :: info
		real(8), intent(inout) :: x(n),v(n)
		complex(8) :: H(latt%Ns*spin,latt%Ns*spin)
		real(8) :: E(size(H,1))
		complex(8) :: D(size(H,1),1,n),cH(size(H,1),size(H,2)),bd
		real(8) :: f(size(E)),dvar(size(var(1:))),err,dE(size(H,1),n),fE
		integer :: l,i=0,k
		common fE
		v=0d0
		dvar=0d0
		fE=0d0
		var(1:n)%val=x
		!$OMP PARALLEL DO REDUCTION(+:fE,v,dvar) PRIVATE(H,E,cH,f,dE,D) if(size(brizon%k,1)/=1)
		do k=1,size(brizon%k,1)
			if(info>0.or.size(brizon%k,1)/=1) then
				call var%Hamilton(H,brizon%k(k,:))
				call mat_diag(H,E)
			endif
			cH=transpose(conjg(H))
			f=1d0/(exp(E/Tk)+1d0)
			do l=1,n
				call dHamilton(var(l),H,cH,D(:,1:1,l),brizon%k(k,:))
			enddo
			dE=real(D(:,1,:))
			do l=1,size(E)
				if((-E(l))>1d2*Tk) then
					fE=fE+E(l)
				else
					fE=fE-Tk*log(1d0+exp(-E(l)/Tk))
				endif
			enddo
			do l=1,size(var(1:))
				if(l<=n) v(l)=v(l)+sum(dE(:,l)*f(:))
				select case(var(l)%sg)
				case(1)
					if(l<=n) v(l)=v(l)+real(sum(var(l)%bd_sg)*(1d0-nf))
					fE=fE+var(l)%val*real(sum(var(l)%bd_sg)*(1d0-nf))
				case default
					if(var(l)%nb==0) then
						dvar(l)=dvar(l)+2d0*sum(abs(var(l)%bd_sg)**2)*var(l)%V
						if(var(l)%sg==3) then
							if(l<=n) v(l)=v(l)+real(sum(var(l)%bd_sg)*var(l)%V)
							fE=fE+var(l)%val*real(sum(var(l)%bd_sg)*var(l)%V)
						else
							if(l<=n) v(l)=v(l)-real(sum(var(l)%bd_sg)*var(l)%V)
							fE=fE-var(l)%val*real(sum(var(l)%bd_sg)*var(l)%V)
						endif
					else
						dvar(l)=dvar(l)+4d0*sum(abs(var(l)%bd_sg)**2)*var(l)%V
					endif
				end select
			enddo
		enddo
		!$OMP END PARALLEL DO
		do l=1,n
			if(info<0.and.abs(dvar(l))>1d-8) then
				x(l)=-v(l)/dvar(l)
				!x(l)=-v(l)/dvar(l)*0.6d0+var(l)%val*0.4d0
			endif
			v(l)=(v(l)+dvar(l)*var(l)%val)/(size(var(l)%n)*size(brizon%k,1))
		enddo
		fE=(fE+sum(dvar*var(1:)%val**2)/2d0)/(latt%Ns*size(brizon%k,1))
		err=sum(abs(v))/size(v)

		!if(n>1) then
			!i=i+1
			!write(20,"(es17.9$)")fE
			!write(20,"(es17.9)")err,var(1)%val
			!!if(info<0) then
				!!call export_data(10)
			!!endif
			!write(*,"(i5$)"),i
		!endif
		!write(*,"(es15.7$)")fE
		!if(n<8) then
			!write(*,"(es12.4$)")var(1:n)%val,err
		!else
			!write(*,"(es12.4$)")var(1)%val,var(2)%val,var(latt%Ns*2+2)%val,var(latt%Ns*2+latt%Ns+2)%val,err
		!endif
		!write(*,"(x)")

		info=1
		!if(err<1d-7) then
			!info=-1
		!endif
	end subroutine
	function findTc(sg,nb,l,is,Tm)
		integer :: sg,nb,l,is
		real(8) :: findTc,Tm
		real(8) :: dTk,er,pod,order
		integer :: getn(2),isg
		er=1d-3
		getn=get_var(sg,nb,l)
		if(is<0) then
			Tk=1d-5
		else
			Tk=Tm
		endif
		dTk=1d-3
		call self_consist()
		pod=sum(abs(var(getn(1):getn(2))%val))/(getn(2)-getn(1)+1)-er
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
			order=sum(abs(var(getn(1):getn(2))%val))/(getn(2)-getn(1)+1)
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
			call var%Hamilton(H,brizon%k(n,:))
			call mat_diag(H,E)
			cH=transpose(conjg(H))
			f=1d0/(1d0+exp(bt*E))

			!call var%Hamilton(Bg(:,:,1),brizon%k(n,:),b1g)
			!call var%Hamilton(Bg(:,:,2),brizon%k(n,:),b2g)
			call var(:0)%Hamilton(Bg(:,:,1),brizon%k(n,:),b1g)
			call var(:0)%Hamilton(Bg(:,:,2),brizon%k(n,:),b2g)
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
			call var%Hamilton(H,brizon%k(n,:))
			call mat_diag(H,E)
			cH=transpose(conjg(H))
			f=1d0/(1d0+exp(bt*E))

			call var(:0)%Hamilton(Bg(:,:,1),brizon%k(n,:),b1g)
			call var(:0)%Hamilton(Bg(:,:,2),brizon%k(n,:),b2g)
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
			call var%Hamilton(H,k)
			call mat_diag(H,E)
			cH=transpose(conjg(H))
			f=1d0/(1d0+exp(bt*E))
			df=-bt/((1d0+cosh(bt*E))*2)

			call var%Hamilton(vv,k,V)
			call var%Hamilton(mm,k,M)
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
	logical :: f
	integer :: l,m
	real(8) :: n(63)=&
		(/9.5000d-01, 9.4500d-01, 9.4000d-01, 9.3500d-01, 9.3000d-01, 9.2500d-01, 9.2000d-01, 9.1500d-01, 9.1000d-01, 9.0500d-01, 9.0000d-01, 8.9500d-01, 8.9000d-01, 8.8500d-01, 8.8000d-01, 8.7500d-01, 8.7000d-01, 8.6900d-01, 8.6800d-01, 8.6700d-01, 8.6600d-01, 8.6500d-01, 8.6400d-01, 8.6300d-01, 8.6200d-01, 8.6100d-01, 8.6000d-01, 8.5900d-01, 8.5800d-01, 8.5700d-01, 8.5600d-01, 8.5500d-01, 8.5400d-01, 8.5300d-01, 8.5200d-01, 8.5100d-01, 8.5000d-01, 8.4900d-01, 8.4800d-01, 8.4700d-01, 8.4600d-01, 8.4500d-01, 8.4400d-01, 8.4300d-01, 8.4200d-01, 8.4100d-01, 8.4000d-01, 8.3900d-01, 8.3800d-01, 8.3700d-01, 8.3600d-01, 8.3500d-01, 8.3400d-01, 8.3300d-01, 8.3200d-01, 8.3100d-01, 8.3000d-01, 8.2500d-01, 8.2000d-01, 8.1500d-01, 8.1000d-01, 8.0500d-01, 8.0000d-01/)
	real(8) :: Ts(size(n),2),Td(size(n),2),Tc(2)
	real(8), allocatable :: peak(:)
	f=openfile(unit=10,file='../data/phase.dat')
	f=openfile(unit=20,file='../data/band.dat')
	f=openfile(unit=30,file='../data/lattice.dat')
	f=openfile(unit=40,file='../data/fermis.dat')
	f=openfile(unit=50,file='../data/DOS.dat')
	f=openfile(unit=60,file='../data/raman.dat')
	f=openfile(unit=70,file='../data/map_raman.dat')
	f=openfile(unit=80,file='../data/map_band.dat')
	f=openfile(unit=90,file='../data/energy.dat')

	call initial()

	do l=6,6
		nf=0.865d0
		Tk=0.022d0/10d0*l
		write(*,"(es12.4$)")nf,Tk
		call self_consist()
		write(*,"(es12.4$)")var(1:)%val
		write(*,"(x)")
		call raman(60,0.003d0,(/0d0,0.3d0/),256,peak)
	enddo
	call latt%gen_origin_brizon((/1d0,0d0,0d0/),(/0d0,1d0,0d0/),o_brizon)
	write(*,"(A$)")"raman peak"
	write(*,"(es12.4$)")peak
	write(*,"(x)")
	do l=1,size(peak)
		call raman_k(70,0.003d0,o_brizon%k,peak(l))
	enddo

	!call band_e(80,0.01d0,(/0d0,0d0,0d0/),(/pi,pi,0d0/),128,(/-1d0,1d0/),512)
	!call band_e(80,0.01d0,(/pi,pi,0d0/),(/0d0,pi,0d0/),128,(/-1d0,1d0/),512)
	!call band_e(80,0.01d0,(/0d0,pi,0d0/),(/0d0,0d0,0d0/),128,(/-1d0,1d0/),512)
	call band(20,(/0d0,0d0,0d0/),(/pi,pi,0d0/),128)
	call band(20,(/pi,pi,0d0/),(/0d0,pi,0d0/),128)
	call band(20,(/0d0,pi,0d0/),(/0d0,0d0,0d0/),128)
	call band(20,(/pi/2d0,pi/2d0,0d0/),(/pi,0d0,0d0/),128)
	call energy(90,o_brizon%k)
	call DOS(50,0.005d0,(/-0.4d0,0.4d0/),512,peak)
	write(*,"(A$)")"DOS peak"
	write(*,"(es12.4$)")peak
	write(*,"(x)")
	do l=1,size(peak)
		call fermis(40,0.005d0,o_brizon%k,peak(l))
	enddo
	var(3)%val=0d0
	call fermis(40,0.005d0,o_brizon%k,0d0)
	stop

	Tc=1d-1
	do l=1,size(n)
		nf=n(l)
		Tc(1)=findTc(ipg(1),ipg(2),1,1,Tc(1)+0.02)
		if(Tc(1)<1d-4) then
			Td(l:,:)=0d0
			exit
		endif
		Tc(2)=findTc(ipg(1),ipg(2),1,-1,Tc(1))
		write(10,"(es12.4$)")nf,Tc
		write(*,"(es12.4$)")nf,Tc
		write(10,"(x)")
		write(*,"(x)")
	enddo
	write(10,"(x/)")

	Tc=1d-1
	do l=1,size(n)
		nf=n(l)
		Tc(1)=findTc(2,1,1,1,Tc(1)+0.02)
		if(Tc(1)<1d-4) then
			exit
		endif
		Tc(2)=findTc(2,1,1,-1,Tc(1))
		write(10,"(es12.4$)")nf,Tc
		write(*,"(es12.4$)")nf,Tc
		write(10,"(x)")
		write(*,"(x)")
	enddo
	write(10,"(x/)")

end program
