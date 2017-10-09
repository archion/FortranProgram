module pmt
	use M_const
	use mkl_service
	implicit none
	real(8), parameter :: t(2)=(/1d0,-0.22d0/),&
		V=0d0,U=2.1d0
	integer :: ic,is,isc
end module
module selfcons
	use pmt
	use M_hamilton_final
	use M_utility
	implicit none
	include 'nlopt.f'
contains
	subroutine initial()
		integer :: l2,l3,l,i
		real(8), allocatable :: bd0(:),bd1(:)
		allocate(var(-10:10))
		call init_random_seed()

		! lattice 
		latt%a1=(/1d0,0d0,0d0/)
		latt%a2=(/0d0,1d0,0d0/)
		latt%c1=(/8d0,0d0,0d0/)
		latt%c2=(/0d0,2d0,0d0/)
		latt%T1=(/1d0,0d0,0d0/)*96
		latt%T2=(/0d0,1d0,0d0/)*96
		latt%bdc(1)=1d0
		latt%bdc(2)=1d0
		allocate(latt%rsb(1,3))
		latt%rsb(1,:)=(/0d0,0d0,0d0/)
		latt%layer=1
		call latt%gen_latt(size(t))
		call latt%gen_brizon(brizon)
		call check_lattice(101)
		write(*,*)"Total site number is: ",latt%Ns

		! cp
		call gen_var(tp=1,nb=0)
		var(iv(0))%bd=-1d0
		var(iv(0))%val=0.02d0

		! d-wave sc
		call gen_var(tp=2,nb=1)
		do i=1,size(var(iv(0))%bd)
			if(abs(latt%nb(var(iv(0))%nb)%bd(i)%dr(1))>1d-6) then
				var(iv(0))%bd(i)=1d0
			else
				var(iv(0))%bd(i)=-1d0
			endif
			!PDW
			if(any(abs(latt%nb(var(iv(0))%nb)%bd(i)%r(1)-[4d0,8d0]+1d0)<1d-5)) then
				var(iv(0))%bd(i)=0d0
			elseif(latt%nb(var(iv(0))%nb)%bd(i)%r(1)>3.1d0) then
				var(iv(0))%bd(i)=-var(iv(0))%bd(i)
			endif
		enddo
		var(iv(0))%val=0.05d0
		isc=iv(0)


		! sdw
		call gen_var(tp=5,nb=0)
		var(iv(0))%bd(:)=&
			![0.4070,0.4372,0.4680,0.3948,0.4680,0.4372,0.4070,0.4808,0.4680,0.4372,0.4070,0.4808,0.4070,0.4372,0.4680,0.3948]-&
			![0.4680,0.4372,0.4070,0.4808,0.4070,0.4372,0.4680,0.3948,0.4070,0.4372,0.4680,0.3948,0.4680,0.4372,0.4070,0.4808
			[0.3650,0.4300,0.5100,0.3250,0.5100,0.4300,0.3650,0.5650,0.5100,0.4300,0.3650,0.5650,0.3650,0.4300,0.5100,0.3250]-&
			[0.5100,0.4300,0.3650,0.5650,0.3650,0.4300,0.5100,0.3250,0.3650,0.4300,0.5100,0.3250,0.5100,0.4300,0.3650,0.5650]
			![0.2200,0.4400,0.6600,0.1800,0.6600,0.4400,0.2200,0.7200,0.6600,0.4400,0.2200,0.7200,0.2200,0.4400,0.6600,0.1800]-&
			![0.6600,0.4400,0.2200,0.7200,0.2200,0.4400,0.6600,0.1800,0.2200,0.4400,0.6600,0.1800,0.6600,0.4400,0.2200,0.7200]
		var(iv(0))%val=-0.5d0*U
		is=iv(0)

		! onsite cdw
		call gen_var(tp=3,nb=0)
		var(iv(0))%bd(:)=&
			![0.4070,0.4372,0.4680,0.3948,0.4680,0.4372,0.4070,0.4808,0.4680,0.4372,0.4070,0.4808,0.4070,0.4372,0.4680,0.3948]+&
			![0.4680,0.4372,0.4070,0.4808,0.4070,0.4372,0.4680,0.3948,0.4070,0.4372,0.4680,0.3948,0.4680,0.4372,0.4070,0.4808]
			[0.3650,0.4300,0.5100,0.3250,0.5100,0.4300,0.3650,0.5650,0.5100,0.4300,0.3650,0.5650,0.3650,0.4300,0.5100,0.3250]+&
			[0.5100,0.4300,0.3650,0.5650,0.3650,0.4300,0.5100,0.3250,0.3650,0.4300,0.5100,0.3250,0.5100,0.4300,0.3650,0.5650]
			![0.2200,0.4400,0.6600,0.1800,0.6600,0.4400,0.2200,0.7200,0.6600,0.4400,0.2200,0.7200,0.2200,0.4400,0.6600,0.1800]+&
			![0.6600,0.4400,0.2200,0.7200,0.2200,0.4400,0.6600,0.1800,0.2200,0.4400,0.6600,0.1800,0.6600,0.4400,0.2200,0.7200]
		var(iv(0))%val=0.5d0*U
		ic=iv(0)

		! hp
		do l=1,size(t)
			call gen_var(tp=-3,nb=l)
			var(iv(0))%bd=-1d0
			var(iv(0))%val=t(l)
		enddo

		call var_init()
	end subroutine
	subroutine rpa_spin(ut,q,gm,omg,l)
		integer :: ut,l
		real(8) :: gm,omg(:),q(:)
		complex(8) :: Uk(latt%Ni*spin,latt%Ni*spin),Ukq(size(Uk,1),size(Uk,2)),Xq(l,brizon%nq*latt%sb,brizon%nq*latt%sb),iomg(l),tmp(size(Xq,2),size(Xq,3)),UXq
		real(8) :: ek(size(Uk,1)),ekq(size(Uk,1)),fk(size(Uk,1)),fkq(size(Uk,1)),bt,Vrpa(size(Xq,2),size(Xq,3)),dr(3)
		integer :: i,j,m,n,nk,l1,l2,isb,jsb,ik,map1(brizon%nq),map2(brizon%nq,brizon%nq)
		type(t_var) :: Vvar(1)
		write(*,*)q
		bt=1d0/Tk
		iomg=(omg(1)+[0:l-1]*(omg(2)-omg(1))/l)+img*gm
		Xq=0d0
		Vrpa=0d0
		spin=1
		call set_var(Vvar(1),tp=3,nb=0)
		Vvar(1)%bd=-1d0
		Vvar(1)%val=U
		call Hamilton(Vvar,tmp,q)
		Vrpa=real(matmul(matmul(Uqi(1:latt%Ni,1:latt%Ni),tmp),conjg(transpose(Uqi(1:latt%Ni,1:latt%Ni)))))
		spin=2

		map1=0; map2=0
		do i=1,brizon%nq
			do j=1,brizon%nq
				if(.not.is_in(brizon%q(i,:)+brizon%q(j,:),brizon%Ta,dr)) then
					dr=dr+brizon%q(i,:)+brizon%q(j,:)
				else
					dr=brizon%q(i,:)+brizon%q(j,:)
				endif
				do m=1,brizon%nq
					if(sum(abs(dr-brizon%q(m,:)))<1d-6) then
						map2(i,j)=m
						exit
					endif
				enddo
				if(sum(abs(dr))<1d-6) then
					if(map1(i)/=0) stop "map1 error"
					map1(i)=j
				endif
			enddo
		enddo
		if(any(map1==0).or.any(map2==0)) stop "map error"
		associate(nq=>brizon%nq,sb=>latt%sb,Ni=>latt%Ni)
		!$OMP PARALLEL DO REDUCTION(+:Xq) PRIVATE(Uk,Ek,Ukq,Ekq,fk,fkq,UXq)
		do nk=1,brizon%nk
			call Hamilton(var,Uk,-brizon%k(nk,:))
			call Hamilton(var,Ukq,brizon%k(nk,:)+q)
			call heev(Uk,Ek,"V")
			call heev(Ukq,Ekq,"V")
			Uk=matmul(Uqi,Uk)
			Ukq=matmul(Uqi,Ukq)
			fk=1d0/(1d0+exp(bt*Ek))
			fkq=1d0/(1d0+exp(bt*Ekq))

			do isb=1,sb
				do jsb=1,sb
					do n=1,size(Uk,1)
						do m=1,size(Uk,2)
							do i=1,nq
								do j=1,nq
									UXq=0d0
									do l1=1,nq
										do l2=1,nq
											UXq=UXq+Ukq(sb*(map2(i,l1)-1)+isb,n)*conjg(Ukq(sb*(map2(j,l2)-1)+jsb,n))*Uk(Ni+sb*(map1(l1)-1)+isb,m)*conjg(Uk(Ni+sb*(map1(l2)-1)+jsb,m))&
												-Ukq(Ni+sb*(map1(map2(i,l1))-1)+isb,n)*conjg(Ukq(sb*(map2(j,l2)-1)+jsb,n))*Uk(sb*(l1-1)+isb,m)*conjg(Uk(Ni+sb*(map1(l2)-1)+jsb,m))
										enddo
									enddo
									Xq(:,sb*(i-1)+isb,sb*(j-1)+jsb)=Xq(:,sb*(i-1)+isb,sb*(j-1)+jsb)+(fk(m)+fkq(n)-1d0)/(iomg-Ek(m)-Ekq(n))*UXq
								enddo
							enddo
						enddo
					enddo
				enddo
			enddo
		enddo
		!$OMP END PARALLEL DO
		Xq=Xq/(brizon%nk*brizon%nq)
		write(ut,"('#q=',3es12.4)")q/pi
		do i=1,l
			tmp=diag(1d0,size(Vrpa,1))+matmul(Vrpa,Xq(i,:,:))
			call mat_inv(tmp)
			write(ut,"(*(es17.9))")real(iomg(i)),-sum(Xq(i,1:sb,1:sb)),-sum(matmul(Xq(i,1:sb,:),tmp(:,1:sb)))
		enddo
		end associate
		write(ut,"(x)")
	end subroutine
	subroutine rpa_spin_i(ut,q,gm,omg,l)
		integer :: ut,l
		real(8) :: gm,omg(:),q(:)
		complex(8) :: Uk(latt%Ni*spin,latt%Ni*spin),Ukq(size(Uk,1),size(Uk,2)),Xq(l,latt%Ni,latt%Ni),iomg(l),tmp(size(Xq,2),size(Xq,3)),UXq
		real(8) :: ek(size(Uk,1)),ekq(size(Uk,1)),fk(size(Uk,1)),fkq(size(Uk,1)),bt,Vrpa(size(Xq,2),size(Xq,3)),dr(3)
		integer :: i,j,m,n,nk,l1,l2,isb,jsb,ik
		type(t_var) :: Vvar(1)
		write(*,*)q
		bt=1d0/Tk
		iomg=(omg(1)+[0:l-1]*(omg(2)-omg(1))/l)+img*gm
		Xq=0d0
		Vrpa=0d0
		spin=1
		call set_var(Vvar(1),tp=3,nb=0)
		Vvar(1)%bd=-1d0
		Vvar(1)%val=U
		call Hamilton(Vvar,tmp,q)
		Vrpa=real(tmp)*brizon%nq
		spin=2

		associate(nq=>brizon%nq,sb=>latt%sb,Ni=>latt%Ni)
		!$OMP PARALLEL DO REDUCTION(+:Xq) PRIVATE(Uk,Ek,Ukq,Ekq,fk,fkq,UXq)
		do nk=1,brizon%nk
			call Hamilton(var,Uk,-brizon%k(nk,:))
			call Hamilton(var,Ukq,brizon%k(nk,:)+q)
			call heev(Uk,Ek,"V")
			call heev(Ukq,Ekq,"V")
			fk=1d0/(1d0+exp(bt*Ek))
			fkq=1d0/(1d0+exp(bt*Ekq))

			do n=1,size(Uk,1)
				do m=1,size(Uk,2)
					do i=1,Ni
						do j=1,Ni
							Xq(:,i,j)=Xq(:,i,j)+(fk(m)+fkq(n)-1d0)/(iomg-Ek(m)-Ekq(n))*(&
								Ukq(i,n)*conjg(Ukq(j,n))*Uk(Ni+i,m)*conjg(Uk(Ni+j,m))&
								-Ukq(Ni+i,n)*conjg(Ukq(j,n))*Uk(i,m)*conjg(Uk(Ni+j,m))&
								)
						enddo
					enddo
				enddo
			enddo
		enddo
		!$OMP END PARALLEL DO
		Xq=Xq/(brizon%nk*brizon%nq)
		write(ut,"('#q=',3es12.4)")q/pi
		do i=1,l
			tmp=diag(1d0,size(Vrpa,1))+matmul(Vrpa,Xq(i,:,:))
			call mat_inv(tmp)
			write(ut,"(*(es17.9))")real(iomg(i)),-sum(Xq(i,:,:)),-sum(matmul(Xq(i,:,:),tmp(:,:)))
		enddo
		end associate
		write(ut,"(x)")
	end subroutine
end module
program main
	use selfcons
	implicit none
	logical :: f,flag=.true.
	integer :: l,m,i,nopt,od,j
	open(20,file="../data/band.dat")
	open(50,file="../data/phyri.dat")
	open(70,file="../data/spect_rpa.dat")
	open(60,file="../data/EDCmap_rpa.dat")
	open(101,file='../data/lattice.dat')
	call omp_set_nested(.false.)
	!!call omp_set_max_active_levels(2)
	call omp_set_dynamic(.false.)
	call mkl_set_dynamic(0)
	call omp_set_schedule(omp_sched_static,0)

	call mkl_set_num_threads(40)
	call omp_set_num_threads(mkl_get_max_threads())


	call initial()
	!do i=1,size(latt%nb(0)%bd)
		!write(50,"(*(es12.4))")latt%nb(0)%bd(i)%r,real(var(ic)%val(1)*var(ic)%bd(i)+var(is)%val(1)*var(is)%bd(i)),real(var(ic)%val(1)*var(ic)%bd(i)-var(is)%val(1)*var(is)%bd(i))
	!enddo
	!write(50,"(x/)")
	!do i=1,size(latt%nb(1)%bd)
		!write(50,"(*(es12.4))")latt%nb(1)%bd(i)%r,latt%nb(1)%bd(i)%dr,real(var(isc)%val(1)*var(isc)%bd(i))
	!enddo
	!stop

	Tk=0.0001d0

	nf=1d0-0.12d0

	!call band(20,[0d0,0d0,0d0],brizon%Ta(1,:),512)
	!call band(20,brizon%Ta(1,:),0.5d0*(brizon%Ta(1,:)+brizon%Ta(2,:)),512)
	!call band(20,0.5d0*(brizon%Ta(1,:)+brizon%Ta(2,:)),[0d0,0d0,0d0],512)
	!stop

	j=nint(sqrt(real(size(brizon%k,1))))
	do i=j/4,3*j/4
		if(i<j/2) then
			!call rpa_spin(70,[2d0*pi*i/j,2d0*pi*i/j,0d0],0.002d0,[0d0,0.5d0],1000)
			!call rpa_spin_i(70,[2d0*pi*i/j,2d0*pi*i/j,0d0],0.02d0,[0d0,1d0],500)
		else
			!call rpa_spin(70,[pi-2d0*pi*(i-j/2)/j,pi,0d0],0.002d0,[0d0,0.5d0],1000)
			call rpa_spin_i(70,[pi-2d0*pi*(i-j/2)/j,pi,0d0],0.02d0,[0d0,1d0],500)
			!call rpa_spin_i(70,[pi,pi-2d0*pi*(i-j/2)/j,0d0],0.02d0,[0d0,1d0],500)
		endif
	enddo

end program
