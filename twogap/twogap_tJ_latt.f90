module pmt
	use M_const
	implicit none
	real(8), parameter :: t(3)=(/1d0,-0.25d0,0.1d0/),&
		V=0.12d0,DJ=0.35d0,&
		!V=0.d0,DJ=0.25d0,&
	Vs=DJ-V,Vd=0.5d0*DJ+V,Tk=0.01,nf=0.8d0
end module
module selfcons
	use pmt
	use M_Hamilton
	use M_utility
	implicit none
contains
	subroutine initial()
		integer :: iv(-1:1),n,i,j,l,ip,k
		type(t_var), allocatable :: tmp(:)
		allocate(tmp(-10000:10000))
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
		call latt%gen_neb(3)
		call latt%gen_bond(3)
		brizon%n1=128
		brizon%n2=128
		call latt%gen_brizon()
		call check_lattice(30)
		write(*,*)"Total site number is: ",latt%Ns

		iv=(/1,0,0/)
		! cp
		call gen_var(iv,n,sg=1,nb=0,V=1d0,var=tmp)
		tmp(iv(0))%val=1.52d-01
		tmp(iv(0))%bd_sg=-1d0

		! ddw
		call gen_var(iv,n,sg=3,nb=1,V=Vd,var=tmp)
		do i=iv(0)-n+sign(1,n),iv(0),sign(1,n)
			tmp(i)%val=2d-1
			do k=1,size(tmp(i)%n)
				tmp(i)%bd_sg(k)=ab(latt%bond(tmp(i)%nb)%bd(tmp(i)%n(k))%i(1))*img*dwave(tmp(i)%n(k))
			enddo
		enddo

		! d-wave sc
		call gen_var(iv,ip,sg=2,nb=1,V=Vs,var=tmp)
		do i=iv(0)-n+sign(1,n),iv(0),sign(1,n)
			tmp(i)%val=2d-1
			do k=1,size(tmp(i)%n)
				tmp(i)%bd_sg(k)=dwave(tmp(i)%n(k))
			enddo
		enddo

		!! sdw
		!call gen_var(iv,ip,sg=4,nb=0,V=DJ,var=tmp)
		!do i=iv(0)-n+sign(1,n),iv(0),sign(1,n)
			!tmp(i)%val=2d-1
			!do k=1,size(tmp(i)%n)
				!tmp(i)%bd_sg(k)=ab(tmp(i)%n(k))
			!enddo
		!enddo

		! bond order
		call gen_var(iv,ip,sg=3,nb=1,V=Vd,var=tmp)
		do i=iv(0)-n+sign(1,n),iv(0),sign(1,n)
			tmp(i)%val=0d0
			do k=1,size(tmp(i)%n)
				tmp(i)%bd_sg(k)=1d0
			enddo
		enddo

		! hp
		do l=1,size(t)
			call gen_var(iv,ip,sg=-3,nb=l,V=1d0,var=tmp)
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
	subroutine self_consist()
		integer :: info
		real(8) :: x(size(var(1:))),v(size(var(1:))),wa(nint((size(var(1:))*(3*size(var(1:))+13))/2.)+10)
		info=1
		x=var(1:)%val
		call hybrd1(do_var,size(x),x,v,1d-7,info,wa,size(wa))
		var(1:)%val=x
		write(*,*)info
	end subroutine
	subroutine do_var(n,x,v,info)
		integer, intent(in) :: n
		integer, intent(inout) :: info
		real(8), intent(inout) :: x(n),v(n)
		complex(8) :: H(latt%Ns*2,latt%Ns*2)
		real(8) :: E(latt%Ns*2)
		complex(8) :: D(size(H,1),size(H,2),n),cH(size(H,1),size(H,2)),bd,tmp(2,2)
		real(8) :: f(size(E)),dvar,fE,err,dE(size(H,1),n)
		integer :: k,l,i=0
		i=i+1
		var(1:n)%val=x
		fE=0d0
		v=0d0
		write(*,"(es12.4$)")x
		!$OMP PARALLEL DO PRIVATE(H,E,f,cH,fE,D,dE,dvar) REDUCTION(+:v)
		do k=1,size(brizon%k,1)
			call Hamilton(H,brizon%k(k,:))
			call heevd(H,E,"V")
			f=1d0/(exp(E/Tk)+1d0)
			cH=transpose(conjg(H))
			call fenergy(E,fE)
			do l=1,n
				call dHamilton(var(l),H,cH,D(:,1:1,l),brizon%k(k,:))
			enddo
			dE=real(D(:,1,:))
			do l=1,n
				v(l)=v(l)+sum(dE(:,l)*f(:))
				dvar=0d0
				select case(var(l)%sg)
				case(1)
					v(l)=v(l)+real(sum(var(l)%bd_sg)*(1d0-nf))
				case default
					dvar=2d0*sum(abs(var(l)%bd_sg)**2)*var(l)%V
					if(var(l)%nb==0) then
						if(var(l)%sg==3) then
							v(l)=v(l)+real(sum(var(l)%bd_sg)*var(l)%V)
						else
							v(l)=v(l)-real(sum(var(l)%bd_sg)*var(l)%V)
						endif
					else
						dvar=dvar*2d0
					endif
				end select
				v(l)=v(l)+dvar*var(l)%val
			enddo
		enddo
		!$OMP END PARALLEL DO
		do l=1,n
			v(l)=v(l)/(size(var(l)%n)*size(brizon%k,1))
		enddo
		err=sum(abs(v))/size(v)
		write(*,*)err
		if(err<1d-7) then
			info=-1
		else
			info=1
		endif
	end subroutine
	subroutine fenergy(E,fE)
		real(8), intent(in) :: E(:)
		real(8), intent(inout) :: fE
		real(8) :: v2
		integer :: k,l
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
				v2=0d0
				fE=fE+real(var(l)%val*sum(var(l)%bd_sg)*(1d0-nf))
			case default
				v2=sum(abs(var(l)%bd_sg)**2)*var(l)%V*var(l)%val**2
				if(var(l)%nb==0) then
					if(var(l)%sg==3) then
						fE=fE+real(var(l)%val*sum(var(l)%bd_sg)*var(l)%V)
					else
						fE=fE-real(var(l)%val*sum(var(l)%bd_sg)*var(l)%V)
					endif
				else
					v2=v2*2d0
				endif
				fE=fE+v2
			end select
		enddo
		fE=fE/latt%Ns
	end subroutine
end module
program main
	use selfcons
	implicit none
	logical :: f
	f=openfile(unit=20,file='../data/band.dat')
	f=openfile(unit=30,file='../data/lattice.dat')

	call initial()
	call self_consist()
	call band((/0d0,0d0,0d0/),brizon%T(1,:),128,20)
	call band(brizon%T(1,:),brizon%T(1,:)+brizon%T(2,:),128,20)
	call band(brizon%T(1,:)+brizon%T(2,:),(/0d0,0d0,0d0/),128,20)
end program
