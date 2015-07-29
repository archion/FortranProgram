module pmt
	use M_const
	implicit none
	real(8), parameter :: t(3)=(/1d0,-0.25d0,0.1d0/),&
		V=0.12d0,DJ=0.35d0,&
		!V=0.d0,DJ=0.25d0,&
	Vs=DJ-V,Vd=0.5d0*DJ+V,Tk=1d-5,nf=0.2d0
end module
module selfcons
	use pmt
	use M_Hamilton
	use M_utility
	use lapack95, only: heevd
	implicit none
	type(t_latt) :: brilzone
contains
	subroutine initial()
		integer :: iv,i,j,l,ip,k
		type(t_var), allocatable :: tmp(:)
		allocate(tmp(10000))
		call init_random_seed()

		! lattice 
		latt%a1=(/1d0,0d0,0d0/)
		latt%a2=(/0d0,1d0,0d0/)
		latt%T1=(/1d0,1d0,0d0/)
		latt%T2=(/-1d0,1d0,0d0/)
		latt%bdc(1)=1d0
		latt%bdc(2)=1d0
		allocate(latt%sub(1,2))
		latt%sub(1,:)=(/0d0,0d0/)
		latt%layer=1
		call latt%gen_latt()
		call latt%gen_neb(3)
		call latt%gen_bond(3)
		write(*,*)"Total site number is: ",latt%Ns

		iv=0
		! cp
		call gen_var(iv,ip,sg=1,nb=0,V=1d0,var=tmp)
		tmp(iv)%val=1.52d-01
		tmp(iv)%bd_sg=-1d0

		! ddw
		call gen_var(iv,ip,sg=3,nb=1,V=Vd,var=tmp)
		do i=ip,iv
			tmp(i)%val=2d-1
			do k=1,size(tmp(i)%n)
				tmp(i)%bd_sg(k)=ab(latt%bond(tmp(i)%nb)%bd(tmp(i)%n(k))%i(1))*img*dwave(tmp(i)%n(k))
			enddo
		enddo

		! d-wave sc
		call gen_var(iv,ip,sg=2,nb=1,V=Vs,var=tmp)
		do i=ip,iv
			tmp(i)%val=2d-1
			do k=1,size(tmp(i)%n)
				tmp(i)%bd_sg(k)=dwave(tmp(i)%n(k))
			enddo
		enddo

		! sdw
		call gen_var(iv,ip,sg=4,nb=0,V=DJ,var=tmp)
		do i=ip,iv
			tmp(i)%val=1d-1
			do k=1,size(tmp(i)%n)
				tmp(i)%bd_sg(k)=ab(tmp(i)%n(k))
			enddo
		enddo

		! bond order
		call gen_var(iv,ip,sg=3,nb=1,V=Vd,var=tmp)
		do i=ip,iv
			tmp(i)%val=2d-1
			do k=1,size(tmp(i)%n)
				tmp(i)%bd_sg(k)=1d0
			enddo
		enddo


		! hp
		do l=1,size(t)
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
		deallocate(tmp)
	end subroutine
	subroutine self_consist()
		integer :: info
		real(8) :: x(size(var(1:))),v(size(var(1:))),wa(nint((size(var(1:))*(3*size(var(1:))+13))/2.)+10)
		call hybrd1(do_var,size(var(1:)),var(1:)%val,v,1d-7,info,wa,size(wa))
	end subroutine
	subroutine do_var(n,x,v,info)
		integer, intent(in) :: n
		integer, intent(inout) :: info
		real(8), intent(inout) :: x(n),v(n)
		complex(8) :: H(latt%Ns*2,latt%Ns*2)
		real(8) :: E(latt%Ns*2)
		complex(8) :: D(size(H,1),size(H,2),n),cH(size(H,1),size(H,2)),bd
		real(8) :: f(size(E)),dvar,fE,err,dE(size(H,1),n)
		integer :: k,l,i=0
		i=i+1
		var(1:n)%val=x
		fE=0d0
		do k=1,brilzone%Ns
			call Hamilton(H,brilzone%i2r(k,:))
			call heevd(H,E,"V")
			f=1d0/(exp(E/Tk)+1d0)
			cH=transpose(conjg(H))
			call fenergy(E,fE)
			do l=1,n
				call dHamilton(var(l),H,cH,D(:,1:1,l),brilzone%i2r(k,:))
			enddo
			dE=real(D(:,1,:))
			do l=1,n
				v(l)=sum(dE(:,l)*f(:))
				dvar=0d0
				select case(var(l)%sg)
				case(1)
					v(l)=v(l)+real(sum(var(l)%bd_sg)*(1d0-nf))
				case default
					dvar=-2d0*sum(abs(var(l)%bd_sg)**2)*var(l)%V
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
				if(info<0.and.abs(dvar)>1d-8) then
					x(l)=-v(l)/dvar
				endif
				v(l)=(v(l)+dvar*var(l)%val)/(size(var(l)%n))
			enddo
		enddo
		err=sum(abs(v))/size(v)
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
				v2=-sum(abs(var(l)%bd_sg)**2)*var(l)%V*var(l)%val**2
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
	use pmt
	use selfcons
	implicit none
	! brillouin zone
	brilzone%T1=(/1d0,1d0,0d0/)*pi
	brilzone%T2=(/-1d0,1d0,0d0/)*pi
	brilzone%a1=brilzone%T1/512d0
	brilzone%a2=brilzone%T2/512d0
	brilzone%bdc(1)=1d0
	brilzone%bdc(2)=1d0
	allocate(brilzone%sub(1,2))
	brilzone%sub(1,:)=(/0d0,0d0/)
	brilzone%layer=1
	call brilzone%gen_latt()

	call initial()
	call self_consist()
end program
