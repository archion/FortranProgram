module M_pmt
	use M_const
	implicit none
	real(8), parameter :: t(3)=(/1d0,-0.25d0,0d0/),nf=0.85d0,U=2.55d0,DJ=0d0,V=-1d0,cvg=1d-6,Tk=1d-5,Vimp=2d0
end module
module M_bdg
	use M_Hamilton
	use M_rd
	use M_utility
	use M_pmt
	use lapack95, only: heevd
	implicit none
	type, extends(t_sort), private :: t_mysort
		integer :: idx
	contains
		procedure :: swap_sort => myswap1
	end type
contains
	subroutine initial()
		integer :: iv,i,j,l,ip,k
		type(t_var), allocatable :: tmp(:)
		type(t_mysort), allocatable :: sort_site(:),sort_bd1(:)
		integer, allocatable :: collect_site(:),collect_bd1(:)
		real(8) :: q(3)
		allocate(tmp(10000))
		call init_random_seed()
		q=(/1d0/8d0,0d0,0d0/)
		!q=(/1d0/50d0,0d0,0d0/)
		!q=0d0
		! lattice 
		a1=(/1d0,0d0,0d0/)
		a2=(/0d0,1d0,0d0/)
		!T1=(/8d0,8d0,0d0/)
		!T2=(/-6d0,6d0,0d0/)
		T1=a1*30
		T2=a2*30
		bdc(1)=1d0
		bdc(2)=1d0
		allocate(sub(1,2))
		sub(1,:)=(/0d0,0d0/)
		layer=1
		call gen_latt()
		call gen_neb()
		call gen_bond(3)
		! finish i2r(i,2),neb(i)%nb(j)%bond(k)/bdc(k)/r(k,2)
		write(*,*)"Total site number is: ",Ns
		!do k=1,Ns
			!write(101,"(3es13.2,2I5)")i2r(k,:),k,ab(k)
		!enddo
		!write(101,"(1x/)")
		allocate(sort_site(Ns))
		do k=1,size(sort_site)
			!sort_site(k)%val=cos(2d0*pi*sum((bond(0)%bd(k)%r)*q))
			sort_site(k)%val=k
			sort_site(k)%idx=k
			!write(101,"(es13.2$)")bond(0)%bd(k)%r,bond(0)%bd(k)%dir,sort_site(k)%val
			!write(101,"(x)")
		enddo
		!write(101,"(1x/)")
		call sort_site%qsort()
		call sort_site%collect(collect_site)

		allocate(sort_bd1(size(bond(1)%bd)))
		do k=1,size(sort_bd1)
			sort_bd1(k)%val=k
			!sort_bd1(k)%val=cos(2d0*pi*sum((bond(1)%bd(k)%r-(/0.5d0,0d0,0d0/))*q))
			sort_bd1(k)%idx=k
			!write(101,"(es13.2$)")bond(1)%bd(k)%r,bond(1)%bd(k)%dir,sort_bd1(k)%val
			!write(101,"(x)")
		enddo
		!write(101,"(1x/)")
		call sort_bd1%qsort()
		call sort_bd1%collect(collect_bd1)

		iv=0
		! cp
		call gen_var(iv,ip,sg=1,nb=0,V=1d0,var=tmp)
		tmp(iv)%val=-1d-1
		tmp(iv)%bd_sg=-1d0

		call gen_var(iv,ip,(/1,2/),(/nint(sqrt(real(Ns)))/2*nint(sqrt(real(Ns)))+nint(sqrt(real(Ns)))/2/),sg=-4,nb=0,V=1d0,var=tmp)
		tmp(iv)%val=Vimp
		tmp(iv)%bd_sg=1d0

		!! ddw
		!!call gen_var(iv,ip,collect_bd1,sort_bd1%idx,sg=3,nb=1,V=1d0,var=tmp)
		!call gen_var(iv,ip,sg=3,nb=1,V=1d0,var=tmp)
		!do i=ip,iv
			!tmp(i)%val=2d-1
			!do k=1,size(tmp(i)%n)
				!tmp(i)%bd_sg(k)=ab(bond(tmp(i)%nb)%bd(tmp(i)%n(k))%i(1))*img*dwave(tmp(i)%n(k))*&
					!1d0
					!!cos(2d0*pi*sum(q*(bond(tmp(i)%nb)%bd(tmp(i)%n(k))%r-(/0.5d0,0d0,0d0/))))
			!enddo
		!enddo

		! d-wave sc
		call gen_var(iv,ip,collect_bd1,sort_bd1%idx,sg=2,nb=1,V=V,var=tmp)
		!call gen_var(iv,ip,sg=2,nb=1,V=V,var=tmp)
		do i=ip,iv
			call random_number(tmp(i)%val)
			!tmp(i)%val=2d-1
			tmp(i)%val=(tmp(i)%val-0.5d0)*0.1d0+1d-1
			do k=1,size(tmp(i)%n)
				tmp(i)%bd_sg(k)=dwave(tmp(i)%n(k))
			enddo
		enddo

		!call gen_var(iv,ip,sg=2,nb=1,V=1d0,var=tmp)
		!do i=ip,iv
			!tmp(i)%val=1d-5
			!do k=1,size(tmp(i)%n)
				!tmp(i)%bd_sg(k)=dwave(tmp(i)%n(k))*&
					!cos(2d0*pi*sum(q*(bond(tmp(i)%nb)%bd(tmp(i)%n(k))%r-(/0.5d0,0d0,0d0/))))
			!enddo
		!enddo


		!! s-wave sc
		!call gen_var(iv,ip,collect_bd1,sort_bd1%idx,sg=2,nb=1,V=1d0,var=tmp)
		!do i=ip,iv
			!tmp(i)%val=1d-1
			!do k=1,size(tmp(i)%n)
				!tmp(i)%bd_sg(k)=1d0
			!enddo
		!enddo

		! sdw
		call gen_var(iv,ip,collect_site,sort_site%idx,sg=4,nb=0,V=-U,var=tmp)
		!call gen_var(iv,ip,sg=4,nb=0,V=-U,var=tmp)
		do i=ip,iv
			call random_number(tmp(i)%val)
			!tmp(i)%val=1d-1
			tmp(i)%val=(tmp(i)%val-0.5d0)*0.1d0+1d-1
			do k=1,size(tmp(i)%n)
				tmp(i)%bd_sg(k)=ab(tmp(i)%n(k))*&
					1d0
					!sin(2d0*pi*sum(q*(bond(tmp(i)%nb)%bd(tmp(i)%n(k))%r+(/0.5d0,0d0,0d0/))))
					!sin(2d0*pi*sum(q*(bond(tmp(i)%nb)%bd(tmp(i)%n(k))%r)))
			enddo
		enddo

		!! d-wave cdw
		!!call gen_var(iv,ip,collect_bd1,sort_bd1%idx,sg=3,nb=1,V=1d0,var=tmp)
		!call gen_var(iv,ip,sg=3,nb=1,V=1d0,var=tmp)
		!do i=ip,iv
			!tmp(i)%val=1d-1
			!do k=1,size(tmp(i)%n)
				!tmp(i)%bd_sg(k)=dwave(tmp(i)%n(k))
			!enddo
		!enddo

		! on site cdw
		call gen_var(iv,ip,collect_site,sort_site%idx,sg=3,nb=0,V=U,var=tmp)
		!call gen_var(iv,ip,sg=3,nb=0,V=U,var=tmp)
		!deallocate(tmp(iv)%bd_sg,tmp(iv)%n)
		!iv=iv-1
		do i=ip,iv
			call random_number(tmp(i)%val)
			!tmp(i)%val=nf/2d0
			tmp(i)%val=(tmp(i)%val-0.5d0)*0.1d0+nf/2d0
			!tmp(i)%val=(tmp(i)%val-0.5d0)*0.1d0
			do k=1,size(tmp(i)%n)
				tmp(i)%bd_sg(k)=&
					1d0
					!cos(2d0*pi*sum(2d0*q*(bond(tmp(i)%nb)%bd(tmp(i)%n(k))%r+(/0.5d0,0d0,0d0/))))
					!cos(2d0*pi*sum(2d0*q*(bond(tmp(i)%nb)%bd(tmp(i)%n(k))%r)))
			enddo
		enddo

		! hp
		!call gen_var(iv,ip,collect_bd1,sort_bd1%idx,sg=3,nb=1,V=1d0,var=tmp)
		!deallocate(tmp(iv)%bd_sg,tmp(iv)%n)
		!iv=iv-1
		call gen_var(iv,ip,sg=-3,nb=1,V=1d0,var=tmp)
		do i=ip,iv
			tmp(i)%val=t(1)
			do k=1,size(tmp(i)%n)
				tmp(i)%bd_sg(k)=&
					-1d0
					!cos(2d0*pi*sum(q*(bond(tmp(i)%nb)%bd(tmp(i)%n(k))%r-(/0.5d0,0d0,0d0/))))
			enddo
		enddo

		k=tmp(iv)%nb+1
		do l=k,size(t)
			call gen_var(iv,ip,sg=-3,nb=l,V=1d0,var=tmp)
			do i=ip,iv
				tmp(i)%bd_sg=-1d0
				tmp(i)%val=t(l)
			enddo
		enddo


		ip=-count(tmp%sg<0)+1
		k=1
		allocate(var(ip:ip+iv-1))
		do l=1,iv
			if(tmp(l)%sg<0) then
				call move_alloc(tmp(l)%bd_sg,var(ip)%bd_sg)
				call move_alloc(tmp(l)%n,var(ip)%n)
				var(ip)%val=tmp(l)%val
				var(ip)%nb=tmp(l)%nb
				var(ip)%V=tmp(l)%V
				var(ip)%sg=tmp(l)%sg      
				ip=ip+1
			else
				call move_alloc(tmp(l)%bd_sg,var(k)%bd_sg)
				call move_alloc(tmp(l)%n,var(k)%n)
				var(k)%val=tmp(l)%val
				var(k)%nb=tmp(l)%nb
				var(k)%V=tmp(l)%V
				var(k)%sg=tmp(l)%sg      
				k=k+1
			endif
		enddo
		deallocate(tmp,sort_site,sort_bd1,collect_site,collect_bd1)
	end subroutine
	subroutine bdg_var()
		complex(8) :: H(2*Ns,2*Ns)
		real(8) :: wide,E(2*Ns),step,pn,fE,cfE,dt
		integer :: sg,ic,m,i
		logical :: flag
		type(t_var) :: pvar(ubound(var,1))
		dt=0.01d0
		pvar=var(1:)
		m=0
		if(var(1)%sg==1) then
			flag=.true.
			var(1)%val=0d0
			ic=2
		else
			flag=.false.
			ic=1
		endif
		do 
			m=m+1
			if(flag) then
				step=0d0
				pn=0d0
				do
					call Hamilton(H)
					call heevd(H,E,"V")
					call phy(H,E,pvar(1:1))
					call fenergy(E,fE=fE)
					!write(20,"(es15.7)")fE
					write(*,"(es15.7)")fE
					if(abs(pvar(1)%val-nf)<cvg*5d0) then
						exit
					endif
					if(step<1d-10) then
						step=abs(pvar(1)%val-nf)
					endif
					call find_cross(pn,pvar(1)%val-nf,sg)
					if(sg/=0) then
						step=max(step*0.3d0,1d-8)
					endif
					write(*,"(es13.5$)")pvar(1)%val,var(1)%val,step*sign(1d0,nf-pvar(1)%val)
					var(1)%val=var(1)%val-step*sign(1d0,pvar(1)%val-nf)
				enddo
			endif
			do 
				fE=cfE
				pvar(ic:)%val=var(ic:)%val
				call do_var(H,E,var(ic:),dt)
				call fenergy(E,fE=cfE)
				write(*,"(e12.4$)")var(1)%val,var(5)%val,var(Ns*2+5)%val,var(Ns*2+Ns+5)%val
				write(20,"(es15.7)")cfE
				write(*,"(es15.7)")cfE
				if(mod(m,4)==0) then
					call export_data(10)
				endif
				if(cfE>fE) then
					dt=dt/2d0
					var(ic:)%val=pvar(ic:)%val
					cfE=fE
					cycle
				endif
				if(abs(sum(var(ubound(var,1)-Ns+1:ubound(var,1))%val))>1d-1) then
					exit
				endif
				if(dt<1d-5) then
					dt=0.1d0
					exit
				endif
				call Hamilton(H)
				call heevd(H,E,"V")
			enddo
			if(m>100) then
				exit
			endif
		enddo 
	end subroutine
	subroutine bdg()
		complex(8) :: H(2*Ns,2*Ns)
		real(8) :: wide,E(2*Ns),step,pn,fE,cfE,dt
		integer :: sg,ic,m,i
		logical :: flag
		type(t_var) :: pvar(ubound(var,1))
		dt=0.01d0
		pvar=var(1:)
		m=0
		if(var(1)%sg==1) then
			flag=.true.
			var(1)%val=0d0
			ic=2
		else
			flag=.false.
			ic=1
		endif
		do 
			m=m+1
			if(flag) then
				step=0d0
				pn=0d0
				do
					call Hamilton(H)
					call heevd(H,E,"V")
					call phy(H,E,pvar(1:1))
					call fenergy(E,fE=fE)
					!write(20,"(es15.7)")fE
					write(*,"(es15.7)")fE
					if(abs(pvar(1)%val-nf)<cvg*5d0) then
						exit
					endif
					if(step<1d-10) then
						step=abs(pvar(1)%val-nf)
					endif
					call find_cross(pn,pvar(1)%val-nf,sg)
					if(sg/=0) then
						step=max(step*0.3d0,1d-8)
					endif
					write(*,"(es13.5$)")pvar(1)%val,var(1)%val,step*sign(1d0,nf-pvar(1)%val)
					var(1)%val=var(1)%val-step*sign(1d0,pvar(1)%val-nf)
				enddo
			endif
			!do 
				fE=cfE
				pvar(ic:)%val=var(ic:)%val
				call phy(H,E,var(ic:))
				call Hamilton(H)
				call heevd(H,E,"V")
				call fenergy(E,fE=cfE)
				write(*,"(e12.4$)")var(1)%val,var(5)%val,var(Ns*2+5)%val,var(Ns*2+Ns+5)%val
				write(20,"(es15.7)")cfE
				write(*,"(es15.7)")cfE
				if(mod(m,4)==0) then
					call export_data(10)
				endif
				!if(sum(abs(var(ic:)%val-pvar(ic:)%val))/size(var(ic:))<cvg*200) then
				!!if(abs(cfE-fE)<1d-6) then
					!exit
				!endif
			!enddo
			if(m>1000) then
				exit
			endif
		enddo 
	end subroutine
	subroutine do_var(H,E,var,dt)
		complex(8), intent(in) :: H(:,:)
		real(8), intent(in) :: E(:)
		real(8), intent(in) :: dt
		type(t_var), intent(inout) :: var(:)
		complex(8) :: bd,D(size(H,1),size(H,2)),cH(size(H,1),size(H,2))
		real(8) :: f(size(E)),dp,dvar
		integer :: l
		f=1d0/(exp(E/Tk)+1d0)
		cH=transpose(conjg(H))
		!!$OMP PARALLEL DO PRIVATE(D,dvar)
		do l=1,size(var)
			call dHamilton(var(l),H,cH,D(:,1:1))
			dvar=0d0
			select case(var(l)%sg)
			case(1)
				dvar=sum(real(D(:,1))*f(:))
				dvar=dvar+real(sum(var(l)%bd_sg)*(1d0-nf))
			case default
				dvar=sum(real(D(:,1))*f(:))
				if(var(l)%nb==0) then
					if(var(l)%sg==3) then
						dvar=dvar-2d0*(var(l)%val+nf/2d0)*sum(abs(var(l)%bd_sg)**2)*var(l)%V+real(sum(var(l)%bd_sg)*var(l)%V)
						!dvar=dvar/2d0+real(sum(var(l)%bd_sg)*var(l)%V)/2d0
					else
						dvar=dvar-2d0*var(l)%val*sum(abs(var(l)%bd_sg)**2)*var(l)%V-real(sum(var(l)%bd_sg)*var(l)%V)
						!dvar=dvar/2d0-real(sum(var(l)%bd_sg)*var(l)%V)/2d0
					endif
				else
					dvar=dvar-4d0*var(l)%val*sum(abs(var(l)%bd_sg))*var(l)%V
					!dvar=dvar/4d0
				endif
			end select
			var(l)%val=var(l)%val-dvar/size(var(l)%n)*dt
			!var(l)%val=dvar/var(l)%V
		enddo
		!!$OMP END PARALLEL DO
	end subroutine
	subroutine fenergy(E,fE)
		real(8), intent(in) :: E(:)
		real(8), intent(inout) :: fE
		integer :: k,l
		!fE=-Vimp+0.5d0*U*nf**2*Ns
		fE=-Vimp
		do k=1,size(E)
			if((-E(k)/Tk)>1d2) then
				fE=fE+E(k)
			else
				fE=fE-Tk*log(1d0+exp(-E(k)/Tk))
			endif
		enddo
		do l=1,ubound(var,1)
			select case(var(l)%sg)
			case(1)
				fE=fE+real(var(l)%val*sum(var(l)%bd_sg)*(1d0-nf))
			case default
				if(var(l)%nb==0) then
					if(var(l)%sg==3) then
						!fE=fE-sum(abs((var(l)%val+nf/2d0)*var(l)%bd_sg)**2)*var(l)%V+real((var(l)%val)*sum(var(l)%bd_sg)*var(l)%V)
						fE=fE-sum(abs(var(l)%val*var(l)%bd_sg)**2)*var(l)%V+real(var(l)%val*sum(var(l)%bd_sg)*var(l)%V)
					else
						fE=fE-sum(abs(var(l)%val*var(l)%bd_sg)**2)*var(l)%V-real(var(l)%val*sum(var(l)%bd_sg)*var(l)%V)
					endif
				else
					fE=fE-2d0*sum(abs(var(l)%val*var(l)%bd_sg)**2)*var(l)%V
				endif
			end select
		enddo
		fE=fE/Ns
	end subroutine
	subroutine phy(H,E,var)
		complex(8), intent(in) :: H(:,:)
		real(8), intent(in) :: E(:)
		type(t_var), intent(inout) :: var(:)
		real(8) :: f(size(E))
		complex(8) :: bd,val
		integer :: i,k,l,j,n
		f=1d0/(exp(E/Tk)+1d0)
		do l=1,ubound(var,1)
			n=0
			val=0d0
			do k=1,size(var(l)%bd_sg)
				i=bond(var(l)%nb)%bd(var(l)%n(k))%i(1)
				j=bond(var(l)%nb)%bd(var(l)%n(k))%i(2)
				bd=bond(var(l)%nb)%bd(var(l)%n(k))%bdc
				n=n+Ns*2
				select case(var(l)%sg)
				case(2)
					! pair channel
					val=val+var(l)%bd_sg(k)*bd*&
						(0.5d0*sum((H(i,:)*conjg(H(j+Ns,:))+H(j,:)*conjg(H(i+Ns,:)))*f(:)))
				case(3)
					! charge channel
					val=val+var(l)%bd_sg(k)*bd*0.5d0*&
						!(sum(conjg(H(i,:))*H(j,:)*f(:)+H(i+Ns,:)*conjg(H(j+Ns,:))*(1d0-f(:)))-nf)
						(sum(conjg(H(i,:))*H(j,:)*f(:)+H(i+Ns,:)*conjg(H(j+Ns,:))*(1d0-f(:))))
				case(1)
					! chemical potential
					val=val+bd*&
						sum(conjg(H(i,:))*H(j,:)*f(:)+H(i+Ns,:)*conjg(H(j+Ns,:))*(1d0-f(:)))
				case(4)
					! spin channel
					val=val+var(l)%bd_sg(k)*bd*0.5d0*&
						sum(conjg(H(i,:))*H(j,:)*f(:)-H(i+Ns,:)*conjg(H(j+Ns,:))*(1d0-f(:)))
				end select
			enddo
			var(l)%val=real(val)/size(var(l)%n)
		enddo
	end subroutine
	subroutine export_data(ut)
		integer :: ut,l,k
		!rewind(ut)
		!open(unit=ut,file='../data/order.dat')
		do l=1,ubound(var,1)
			if(var(l)%sg/=var(max(l-1,1))%sg.or.var(l)%nb/=var(max(l-1,1))%nb) then
				write(ut,"(x/)")
			endif
			do k=1,size(var(l)%n)
				write(ut,"(es13.2$)")bond(var(l)%nb)%bd(var(l)%n(k))%r,bond(var(l)%nb)%bd(var(l)%n(k))%dir,var(l)%val
				write(ut,"(x)")
			enddo
		enddo
		write(ut,"(x/)")
		!close(ut)
	end subroutine
	subroutine myswap1(self,a)
		class(t_mysort) :: self
		class(t_sort) :: a
		type(t_mysort), allocatable :: tmp
		select type(a)
		type is (t_mysort)
			call self%t_sort%swap_sort(a)
			allocate(tmp)
			tmp%idx=a%idx
			a%idx=self%idx
			self%idx=tmp%idx
			deallocate(tmp)
		end select
	end subroutine
end module
program main
	use M_bdg
	implicit none
	logical :: f
	f=openfile(unit=10,file='../data/order.dat')
	f=openfile(unit=20,file='../data/fenergy.dat')
	call initial()
	call bdg()
	!call bdg_var()
	!export data
	call export_data(10)
end program
