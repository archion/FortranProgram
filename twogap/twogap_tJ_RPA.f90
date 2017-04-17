module pmt
	use M_const
	use mkl_service
	implicit none
	real(8), parameter :: t(3)=(/1d0,-0.45d0,0d0/),&
		!DJ=0.25d0,V=0.0625d0-DJ/4d0
		V=0.5d0/4d0,DJ=0.5d0
		!V=-0.25d0/4d0,DJ=0.25d0
		!V=0d0,DJ=0.25d0
	integer, parameter :: icp=1,isc=2,iddw=3,iubo=4
end module
module selfcons
	use pmt
	use M_hamilton_super
	use M_utility
	implicit none
	include 'nlopt.f'
contains
	subroutine initial()
		integer :: l2,l3,l,i
		real(8), allocatable :: bd0(:),bd1(:)
		allocate(var(1:4))
		call init_random_seed()

		! lattice 
		latt%a1=(/1d0,0d0,0d0/)
		latt%a2=(/0d0,1d0,0d0/)
		!latt%c1=(/1d0,0d0,0d0/)
		!latt%c2=(/0d0,1d0,0d0/)
		latt%c1=(/-1d0,1d0,0d0/)
		latt%c2=(/1d0,1d0,0d0/)
		latt%T1=(/1d0,0d0,0d0/)*256
		latt%T2=(/0d0,1d0,0d0/)*256
		latt%bdc(1)=1d0
		latt%bdc(2)=1d0
		allocate(latt%rsb(1,3))
		latt%rsb(1,:)=(/0d0,0d0,0d0/)
		latt%layer=1
		call latt%gen_latt(3)
		call latt%gen_brizon(brizon)
		!call check_lattice(101)
		write(*,*)"Total site number is: ",latt%Ns

		allocate(var(icp)%val(1))
		var(isc)%V=DJ-V
		allocate(var(isc)%val(1))
		var(iddw)%V=0.5d0*DJ+V
		allocate(var(iddw)%val(1))
		var(iubo)%V=0.5d0*DJ+V
		allocate(var(iubo)%val(1))
		mat_diag => diag2
	end subroutine
	function ab(i)
		integer :: i
		real(8) :: ab
		ab=(-1d0)**(mod(nint(sum(latt%nb(0)%bd(i)%r(:2))),2))
	end function
	subroutine Hamiltons(H,k)
		complex(8) :: H(:,:)
		real(8) :: eks,eka,k(:),dp
		dp=abs(1d0-nf)
		eks=-4d0*dp*t(2)*cos(k(1))*cos(k(2))-2d0*dp*t(3)*(cos(2d0*k(1))+cos(2d0*k(2)))-var(icp)%val(1)
		eka=-2d0*(dp*t(1)+var(iubo)%val(1)*var(iubo)%V)*(cos(k(1))+cos(k(2)))
		H=0d0
		if(size(H,1)==2) then
			H(1,1)=eks+eka
			H(2,2)=-H(1,1)
			H(1,2)=var(isc)%val(1)*var(isc)%V*(cos(k(1))-cos(k(2)))*2d0
			H(2,1)=H(1,2)
		else
			H(1,1)=eks+eka
			H(2,2)=eks-eka
			H(3,3)=-H(1,1)
			H(4,4)=-H(2,2)

			H(1,2)=img*var(iddw)%val(1)*var(iddw)%V*(cos(k(1))-cos(k(2)))*2d0
			H(2,1)=-H(1,2)
			H(3,4)=H(1,2)
			H(4,3)=-H(3,4)

			H(1,3)=var(isc)%val(1)*var(isc)%V*(cos(k(1))-cos(k(2)))*2d0
			H(2,4)=-H(1,3)
			H(3,1)=H(1,3)
			H(4,2)=H(2,4)
		endif
	end subroutine
	subroutine selfconsist()
		complex(8) :: Uk(4,4)
		real(8) :: x(size(var)),Ek(4),al,fk(4),bt
		!complex(8) :: Uk(2,2)
		!real(8) :: x(size(var)),Ek(2),al,fk(2),bt
		integer :: i,j,c
		logical :: is_cross(2)
		is_cross=.false.
		bt=1d0/Tk
		!var([iubo])%val(1)=1d-1
		!var([isc,iubo])%val(1)=1d-1
		!var([isc,iddw,iubo])%val(1)=1d-1
		var([iddw,iubo])%val(1)=1d-1
		do 
			c=0
			al=0.2d0
			do 			
				x=0d0

				!$OMP PARALLEL DO REDUCTION(+:x) PRIVATE(Ek,Uk,fk)
				do i=1,size(brizon%k,1)
					call Hamiltons(Uk,brizon%k(i,:))
					!call mat_diag(Uk,Ek)
					call heev(Uk,Ek,"V")
					fk=1d0/(1d0+exp(bt*Ek))
					!x(icp)=x(icp)+1d0+dot_product(Uk(1,:)*dconjg(Uk(1,:))-Uk(2,:)*dconjg(Uk(2,:)),fk)
					!x(isc)=x(isc)-dot_product(Uk(1,:)*dconjg(Uk(2,:)),fk)*(cos(brizon%k(i,1))-cos(brizon%k(i,2)))/2d0
					!x(iubo)=x(iubo)+dot_product(Uk(1,:)*dconjg(Uk(1,:)),fk)*(cos(brizon%k(i,1))+cos(brizon%k(i,2)))/2d0

					x(icp)=x(icp)+2d0+dot_product(Uk(1,:)*dconjg(Uk(1,:))+Uk(2,:)*dconjg(Uk(2,:))-Uk(3,:)*dconjg(Uk(3,:))-Uk(4,:)*dconjg(Uk(4,:)),fk)
					!x(isc)=x(isc)-dot_product(Uk(1,:)*dconjg(Uk(3,:))-Uk(2,:)*dconjg(Uk(4,:)),fk)*(cos(brizon%k(i,1))-cos(brizon%k(i,2)))/2d0
					x(iddw)=x(iddw)-real(img*dot_product(Uk(1,:)*dconjg(Uk(2,:))-Uk(2,:)*dconjg(Uk(1,:)),fk)*(cos(brizon%k(i,1))-cos(brizon%k(i,2)))/2d0)
					x(iubo)=x(iubo)+dot_product(Uk(1,:)*dconjg(Uk(1,:))-Uk(2,:)*dconjg(Uk(2,:)),fk)*(cos(brizon%k(i,1))+cos(brizon%k(i,2)))/2d0
				enddo
				!$OMP END PARALLEL DO
				c=c+1
				x=x/size(brizon%ok,1)
				!x=x/(2d0*size(brizon%k,1))
				!write(*,*)x(icp)
				if(abs(nf-x(icp))<1d-6) then
					exit
				endif
				is_cross(1)=is_cross(2)
				is_cross(2)=((x(icp)-nf)>0d0)
				if(is_cross(1).neqv.is_cross(2)) then
					al=al*0.3d0
				endif
				var(icp)%val(1)=var(icp)%val(1)-al*sign(1d0,x(icp)-nf)
			enddo
			!if(sum(abs((var([iubo])%val(1)-x([iubo]))))<1d-6) then
			!if(sum(abs((var([isc,iubo])%val(1)-x([isc,iubo]))))<1d-6) then
			!if(sum(abs((var([isc,iddw,iubo])%val(1)-x([isc,iddw,iubo]))))<1d-6) then
			if(sum(abs((var([iddw,iubo])%val(1)-x([iddw,iubo]))))<1d-6) then
				exit
			endif
			!var([iubo])%val(1)=x([iubo])
			!var([isc,iubo])%val(1)=x([isc,iubo])
			!var([isc,iddw,iubo])%val(1)=x([isc,iddw,iubo])
			var([iddw,iubo])%val(1)=x([iddw,iubo])
		enddo
	end subroutine
	subroutine rpainstable(Tk,q,den)
		real(8) :: Tk,q(:),den
		integer :: i,j,m,n
		!complex(8) :: Uk(4,4),Ukq(4,4),gm(4,4)
		!real(8) :: ek(4),ekq(4),fk(4),fkq(4),Vrpa,Xq,bt
		complex(8) :: Uk(2,2),Ukq(2,2),gm(2,2)
		real(8) :: ek(2),ekq(2),fk(2),fkq(2),Vrpa,Xq,bt
		bt=1d0/Tk
		Xq=0d0
		!$OMP PARALLEL DO REDUCTION(+:Xq) PRIVATE(Ek,Uk,fk,Ekq,Ukq,fkq,gm)
		do i=1,size(brizon%k,1)
			call Hamiltons(Uk,brizon%k(i,:))
			call Hamiltons(Ukq,brizon%k(i,:)+q)
			!call heev(Uk,Ek,"V")
			call mat_diag(Uk,Ek)
			!call heev(Ukq,Ekq,"V")
			call mat_diag(Ukq,Ekq)
			fk=1d0/(1d0+exp(bt*Ek))
			fkq=1d0/(1d0+exp(bt*Ekq))
			gm=0d0
			gm(1,1)=(sin(brizon%k(i,1)-q(1)/2d0)-sin(brizon%k(i,2)-q(2)/2d0))*img
			gm(2,2)=-(sin(-brizon%k(i,1)-3d0*q(1)/2d0)-sin(-brizon%k(i,2)-3d0*q(2)/2d0))*img

			gm=matmul(transpose(conjg(Ukq)),matmul(gm,Uk))
			do n=1,size(gm,1)
				do m=1,size(gm,2)
					if(abs(Ek(n)-Ekq(m))<1d-10) then
						Xq=Xq+gm(n,m)*conjg(gm(n,m))*(fk(n)-1d0/(1d0+exp(bt*(Ek(n)-1d-10))))/1d-10
					else
						Xq=Xq+gm(n,m)*conjg(gm(n,m))*(fk(n)-fkq(m))/(Ek(n)-Ekq(m))
					endif
				enddo
			enddo
		enddo
		!$OMP END PARALLEL DO
		Xq=Xq/size(brizon%k,1)
		den=1d0+var(iddw)%V*Xq/2d0
	end subroutine
	subroutine rpa_EDC(ut,q,gm,E)
		integer :: ut
		real(8) :: gm,E,q(:)
		complex(8) :: Uk(4,4),Ukq(4,4),Xq
		real(8) :: ek(4),ekq(4),fk(4),fkq(4),bt
		integer :: i,j,m,n
		bt=1d0/Tk
		Xq=0d0
		do i=1,size(brizon%k,1)
			call Hamiltons(Uk,-brizon%k(i,:))
			call Hamiltons(Ukq,brizon%k(i,:)+q)
			call heev(Uk,Ek,"V")
			call heev(Ukq,Ekq,"V")
			fk=1d0/(1d0+exp(bt*Ek))
			fkq=1d0/(1d0+exp(bt*Ekq))

			do n=1,size(Uk,1)
				do m=1,size(Uk,2)
					Xq=Xq+(fk(m)+fkq(n)-1d0)/(E-Ek(m)-Ekq(n)+img*gm)*(&
						conjg(Ukq(1,n)*Uk(3,m))*Ukq(1,n)*Uk(3,m)&
						+conjg(Ukq(2,n)*Uk(4,m))*Ukq(2,n)*Uk(4,m)&
						-conjg(Ukq(3,n)*Uk(1,m))*Ukq(1,n)*Uk(3,m)&
						-conjg(Ukq(4,n)*Uk(2,m))*Ukq(2,n)*Uk(4,m))
				enddo
			enddo
		enddo
		Xq=Xq/size(brizon%ok,1)
		write(ut,"(*(es17.9))")q(:2),-Xq,-Xq/(1d0+0.34d0*DJ*(cos(q(1))+cos(q(2)))*Xq)
	end subroutine
	subroutine rpa_spin(ut,q,gm,omg,l)
		integer :: ut,l
		real(8) :: gm,omg(:),q(:)
		complex(8) :: Uk(4,4),Ukq(4,4),Xq(l,2,2),iomg(l),tmp(2,2)
		real(8) :: ek(4),ekq(4),fk(4),fkq(4),bt,Vrpa(2,2)
		integer :: i,j,m,n,nk,ll,nq
		bt=1d0/Tk
		nq=size(Uk,1)/2
		iomg=(omg(1)+[0:l-1]*(omg(2)-omg(1))/l)+img*gm
		Xq=0d0
		Vrpa=0d0
		do i=0,nq-1
			!Vrpa(i+1,i+1)=0.34d0*DJ*(cos(q(1)+pi*i)+cos(q(2)+pi*i))
			Vrpa(i+1,i+1)=0.45d0*DJ*(cos(q(1)+pi*i)+cos(q(2)+pi*i))
			!Vrpa(i+1,i+1)=-0.4d0*(1d0-0.2d0*(cos(q(1)+pi*i)+cos(q(2)+pi*i)))
		enddo
		!$OMP PARALLEL DO REDUCTION(+:Xq) PRIVATE(Uk,Ek,Ukq,Ekq,fk,fkq)
		do nk=1,size(brizon%k,1)
			call Hamiltons(Uk,-brizon%k(nk,:))
			call Hamiltons(Ukq,brizon%k(nk,:)+q)
			call heev(Uk,Ek,"V")
			call heev(Ukq,Ekq,"V")
			fk=1d0/(1d0+exp(bt*Ek))
			fkq=1d0/(1d0+exp(bt*Ekq))

			do n=1,size(Uk,1)
				do m=1,size(Uk,2)
					do i=1,nq
						do j=1,nq
							do ll=0,nq-1
								Xq(:,i,j)=Xq(:,i,j)+(fk(m)+fkq(n)-1d0)/(iomg-Ek(m)-Ekq(n))*(&
									+Ukq(1+mod(i+ll-1,nq),n)*conjg(Ukq(1+mod(j+ll-1,nq),n))*Uk(3+mod(nq-ll,nq),m)*conjg(Uk(3+mod(nq-ll,nq),m))&
									-Ukq(3+mod(i+ll-1,nq),n)*conjg(Ukq(1+mod(j+ll-1,nq),n))*Uk(1+mod(nq-ll,nq),m)*conjg(Uk(3+mod(nq-ll,nq),m))&
									+Ukq(1+mod(i+ll-1,nq),n)*conjg(Ukq(1+mod(j+ll,nq),n))*Uk(3+mod(nq-ll,nq),m)*conjg(Uk(3+mod(nq-ll-1,nq),m))&
									-Ukq(3+mod(i+ll-1,nq),n)*conjg(Ukq(1+mod(j+ll,nq),n))*Uk(1+mod(nq-ll,nq),m)*conjg(Uk(3+mod(nq-ll-1,nq),m)))
							enddo
						enddo
					enddo
				enddo
			enddo
		enddo
		!$OMP END PARALLEL DO
		Xq=Xq/size(brizon%ok,1)
		write(ut,"('#q=',3es12.4)")q/pi
		do i=1,l
			tmp=diag(1d0,nq)+matmul(Vrpa,Xq(i,:,:))
			call mat_inv(tmp)
			write(ut,"(*(es17.9))")real(iomg(i)),-Xq(i,1,1),-sum(Xq(i,1,:)*tmp(:,1))
		enddo
		write(ut,"(x)")
	end subroutine
	subroutine rpainstable_spin(Tk,q,den)
		real(8) :: Tk,q(:),den
		integer :: i,j,m,n,l,nq,nk
		complex(8) :: Uk(4,4),Ukq(4,4),gm(4,4),Xq(2,2)
		real(8) :: ek(4),ekq(4),fk(4),fkq(4),bt
		bt=1d0/Tk
		Xq=0d0
		nq=size(Uk,1)/2
		do nk=1,size(brizon%k,1)
			call Hamiltons(Uk,-brizon%k(nk,:))
			call Hamiltons(Ukq,brizon%k(nk,:)+q)
			call heev(Uk,Ek,"V")
			call heev(Ukq,Ekq,"V")
			fk=1d0/(1d0+exp(bt*Ek))
			fkq=1d0/(1d0+exp(bt*Ekq))

			do n=1,size(Uk,1)
				do m=1,size(Uk,2)
					if(abs(-Ek(m)-Ekq(n))>1d-7) then
						do i=1,nq
							do j=1,nq
								do l=0,nq-1
									Xq(i,j)=Xq(i,j)+(fk(m)+fkq(n)-1d0)/(-Ek(m)-Ekq(n))*(&
										+Ukq(1+mod(i+l-1,nq),n)*conjg(Ukq(1+mod(j+l-1,nq),n))*Uk(3+mod(nq-l,nq),m)*conjg(Uk(3+mod(nq-l,nq),m))&
										-Ukq(3+mod(i+l-1,nq),n)*conjg(Ukq(1+mod(j+l-1,nq),n))*Uk(1+mod(nq-l,nq),m)*conjg(Uk(3+mod(nq-l,nq),m))&
										+Ukq(1+mod(i+l-1,nq),n)*conjg(Ukq(1+mod(j+l,nq),n))*Uk(3+mod(nq-l,nq),m)*conjg(Uk(3+mod(nq-l-1,nq),m))&
										-Ukq(3+mod(i+l-1,nq),n)*conjg(Ukq(1+mod(j+l,nq),n))*Uk(1+mod(nq-l,nq),m)*conjg(Uk(3+mod(nq-l-1,nq),m)))
								enddo
							enddo
						enddo
					endif
				enddo
			enddo
		enddo
		Xq=Xq/size(brizon%ok,1)
		!den=(1d0+DJ*(cos(q(1))+cos(q(2)))*real(Xq))
	end subroutine
end module
program main
	use selfcons
	implicit none
	logical :: f,flag=.true.
	integer :: l,m,i,nopt,od,j
	real(8) :: n(0:40)=(/(min(1d0-0.005d0*i,0.999d0),i=0,80,2)/)
	real(8) :: Ts(size(n),2),Td(size(n),2),Tc(2),Tnc,dnf,pnf,pTk(2),den,E(4)
	real(8), allocatable :: peak(:)
	complex(8) :: H(4,4)
	open(unit=10,file=fn('../data/phase.dat'))
	open(unit=30,file=fn('../data/order.dat'))
	open(unit=50,file=fn('../data/incomm.dat'))
	open(unit=110,file=fn('../data/bb.dat'))
	open(70,file="../data/spect_rpa.dat")
	open(60,file="../data/EDCmap_rpa.dat")
	f=openfile(unit=101,file='../data/lattice.dat')

	call initial()
	call omp_set_nested(.false.)
	!!call omp_set_max_active_levels(2)
	call omp_set_dynamic(.false.)
	call mkl_set_dynamic(0)
	call omp_set_schedule(omp_sched_static,0)

	call mkl_set_num_threads(20)
	call omp_set_num_threads(mkl_get_max_threads())

	!!nf=1d0-0.165d0
	!!Tk=0.02d0
	!var([icp,isc,iubo])%val(1)=[-1.4781d-01,8.8716d-06,1.9090d-01]
	!call Hamiltons(H,[3.141592654d+00,4.908738521d-02,0d0])
	!!call mat_diag(H,E)
	!call heev(H,E,"V")
	!write(*,*)E
	!stop
	var(iddw)%val=0d0
	var(isc)%val=0d0

	Tk=0.0001d0
	n(1:4)=[0.23d0,0.24d0,0.245d0,0.25d0]

	nf=1d0-0.12d0
	!call selfconsist()
	!var([icp,isc,iubo])%val(1)=[-1.5799d-01,1.2158d-01,1.8820d-01]
	!var([icp,iddw,iubo])%val(1)=[-1.5799d-01,1.2158d-01,1.8820d-01]
	var([icp,iddw,iubo])%val(1)=[-1.5799d-01,0d0,1.8820d-01]
	!var([icp,iddw,iubo])%val(1)=[-1.8795d-01,1.2178d-01,1.9021d-01]
	!var([icp,isc,iubo])%val(1)=[-1.083d0,0.168d0,0d0]*0.12d0
	write(*,"(*(es12.4))")Tk,nf,var([icp,isc,iubo])%val(1),var(iddw)%val(1)
	j=nint(sqrt(real(size(brizon%ok,1))))
	!do i=0,j/5
		!!do j=0,0
		!call rpainstable_spin(Tk,[pi-2d0*pi*i/j,pi-2d0*pi*i/j,0d0],den)
		!write(*,"(*(es12.4))")1d0-nf,[pi,pi-2d0*pi*i/j],den
		!write(50,"(*(es12.4))")1d0-nf,[pi,pi-2d0*pi*i/j],den
	!enddo
	!stop
	!do i=1,size(brizon%ok,1)
	!do i=128,256,8
		!do j=128,256,8
			!!if(brizon%ok(i,1)>0.2d0*2d0*pi.and.brizon%ok(i,2)>0.2d0*2d0*pi) then
				!!call rpa_EDC(60,brizon%ok(i,:),0.01d0,0.05d0)
			!!endif
			!call rpa_EDC(60,[2d0*pi/512d0*i,2d0*pi/512d0*j],0.01d0,0.3d0)
		!enddo
	!enddo
	!stop
	do i=j/4,3*j/4
		if(i<=j/2) then
			call rpa_spin(70,[2d0*pi*i/j,2d0*pi*i/j,0d0],0.002d0,[0d0,0.5d0],1000)
		else
			call rpa_spin(70,[pi-2d0*pi*(i-j/2)/j,pi,0d0],0.002d0,[0d0,0.5d0],1000)
		endif
	enddo
	!call rpa_spin(70,[pi,pi,0d0],0.02d0,[-10d0,10d0],5000)
	!call rpa_spin(70,[pi,pi,0d0],0.001d0,[0d0,0.5d0],1000)
	stop

	stop
	do i=1,4
		nf=1d0-n(i)
		!read(30,"(5es12.4)")Tk,nf,var([icp,isc,iubo])%val(1)
		!write(*,"(5es12.4)")Tk,nf,var([icp,isc,iubo])%val(1)
		call selfconsist()
		write(*,"(6es12.4)")Tk,nf,var([icp,isc,iubo])%val(1),var(iddw)%val(1)
		!read(*,*)
		!if(i==54) then
			do j=0,1024/4
			!do j=0,0
				call rpainstable(Tk,[pi,pi-pi*j/1024d0,0d0],den)
				write(*,"(5es12.4)")1d0-nf,[pi,pi-pi*j/1024d0],den
				write(50,"(5es12.4)")1d0-nf,[pi,pi-pi*j/1024d0],den
			enddo
		!endif
		write(50,"(x)")
	enddo
	stop

end program
