module pmt
	use M_const
	use mkl_service
	implicit none
	real(8), parameter :: t(3)=(/1d0,-0.45d0,0d0/),&
		!DJ=0.25d0,V=0.0625d0-DJ/4d0
		V=0.5d0/4d0,DJ=0.5d0
		!V=-0.25d0/4d0,DJ=0.25d0
		!V=0d0,DJ=0.25d0
	integer :: icp=-1,isc=-1,iddw=-1,iubo=-1
	!integer :: icp=1,isc=2,iddw=3,iubo=4
end module
include "../test/hamilton_final.f90"
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
		latt%c1=(/1d0,0d0,0d0/)
		latt%c2=(/0d0,1d0,0d0/)
		!latt%c1=(/-1d0,1d0,0d0/)
		!latt%c2=(/1d0,1d0,0d0/)
		latt%c1=(/4d0,0d0,0d0/)
		latt%c2=(/0d0,2d0,0d0/)
		latt%T1=(/1d0,0d0,0d0/)*72
		latt%T2=(/0d0,1d0,0d0/)*72
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
		var(iv(0))%val=0d0
		icp=iv(0)

		! d-wave sc
		call gen_var(tp=2,nb=1,V=-DJ+V)
		do i=1,size(var(iv(0))%bd)
			if(abs(latt%nb(var(iv(0))%nb)%bd(i)%dr(1))>1d-6) then
				var(iv(0))%bd(i)=1d0
			else
				var(iv(0))%bd(i)=-1d0
			endif
		enddo
		var(iv(0))%val=0d0
		isc=iv(0)

		!! ddw
		!call gen_var(tp=3,nb=1,V=-0.5d0*DJ-V)
		!do i=1,size(var(iv(0))%bd)
			!if(abs(latt%nb(var(iv(0))%nb)%bd(i)%dr(mod(nint(sum(latt%nb(0)%bd(latt%nb(var(iv(0))%nb)%bd(i)%i(1))%r)),2)+1))>1d-6) then
				!var(iv(0))%bd(i)=img
			!else
				!var(iv(0))%bd(i)=-img
			!endif
		!enddo
		!var(iv(0))%val=0d0
		!iddw=iv(0)

		!! sdw
		!call gen_var(tp=4,nb=0,V=DJ,Vnb=1)
		!do i=1,size(var(iv(0))%bd)
			!if(latt%nb(var(iv(0))%nb)%bd(i)%sb(1)==1) then
				!var(iv(0))%bd(i)=1
			!else
				!var(iv(0))%bd(i)=-1
			!endif
			!var(iv(0))%bd(i)=ab(i)
		!enddo
		!var(iv(0))%val=1d-1

		! bond order
		call gen_var(tp=3,nb=1,V=-0.5d0*DJ-V)
		var(iv(0))%bd=1d0
		var(iv(0))%val=0d0
		iubo=iv(0)

		! hp
		do l=1,size(t)
			call gen_var(tp=-3,nb=l)
			var(iv(0))%bd=-1d0
			var(iv(0))%val=t(l)
		enddo

		call var_init()
	end subroutine
	subroutine update()
		integer :: l
		do l=lbound(var,1),0
			if(var(l)%tp==-3.and.(.not.var(l)%is_V)) then
				var(l)%bd=-abs(1d0-nf)
			endif
		enddo
	end subroutine
	subroutine Hamiltons(H,k)
		complex(8) :: H(:,:)
		real(8) :: eks,eka,k(:),dp
		dp=abs(1d0-nf)
		eks=-4d0*dp*t(2)*cos(k(1))*cos(k(2))-2d0*dp*t(3)*(cos(2d0*k(1))+cos(2d0*k(2)))-var(icp)%val(1)
		eka=-2d0*(dp*t(1)-var(iubo)%val(1)*var(iubo)%V)*(cos(k(1))+cos(k(2)))
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

			!H(1,3)=var(isc)%val(1)*var(isc)%V*(cos(k(1))-cos(k(2)))*2d0
			!H(2,4)=-H(1,3)
			!H(3,1)=H(1,3)
			!H(4,2)=H(2,4)
		endif
	end subroutine
	subroutine selfconsist()
		complex(8) :: Uk(4,4)
		real(8) :: x(size(var)),Ek(4),al,fk(4),bt
		!complex(8) :: Uk(2,2)
		!real(8) :: x(size(var)),Ek(2),al,fk(2),bt
		integer :: i,j,c
		logical :: is_cross(2)
		call update()
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
				x=x/size(brizon%k,1)
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
		Xq=Xq/size(brizon%k,1)
		write(ut,"(*(es17.9))")q(:2),-Xq,-Xq/(1d0+0.34d0*DJ*(cos(q(1))+cos(q(2)))*Xq)
	end subroutine
	subroutine rpa_spin(ut,q,gm,omg,l)
		integer :: ut,l
		real(8) :: gm,omg(:),q(:)
		complex(8) :: Uk(latt%Ni*spin,latt%Ni*spin),Ukq(size(Uk,1),size(Uk,2)),Xq(l,brizon%nq*latt%sb,brizon%nq*latt%sb),iomg(l),tmp(size(Xq,2),size(Xq,3)),UXq
		real(8) :: ek(size(Uk,1)),ekq(size(Uk,1)),fk(size(Uk,1)),fkq(size(Uk,1)),bt,Vrpa(size(Xq,2),size(Xq,3)),dr(3)
		integer :: i,j,m,n,nk,l1,l2,isb,jsb,ik,map1(brizon%nq),map2(brizon%nq,brizon%nq)
		type(t_var) :: Vvar(1)
		call update()
		write(*,*)q
		bt=1d0/Tk
		iomg=(omg(1)+[0:l-1]*(omg(2)-omg(1))/l)+img*gm
		Xq=0d0
		Vrpa=0d0
		spin=1
		call set_var(Vvar(1),tp=3,nb=1)
		Vvar(1)%bd=0.5d0
		Vvar(1)%val=0.34d0*DJ
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
		write(ut,"(x)")
		end associate
	end subroutine
	subroutine rpa_spin_i(ut,q,gm,omg,l)
		integer :: ut,l
		real(8) :: gm,omg(:),q(:)
		complex(8) :: Uk(latt%Ni*spin,latt%Ni*spin),Ukq(size(Uk,1),size(Uk,2)),Xq(l,latt%Ni,latt%Ni),iomg(l),tmp(size(Xq,2),size(Xq,3)),UXq
		real(8) :: ek(size(Uk,1)),ekq(size(Uk,1)),fk(size(Uk,1)),fkq(size(Uk,1)),bt,Vrpa(size(Xq,2),size(Xq,3)),dr(3)
		integer :: i,j,m,n,nk,l1,l2,isb,jsb,ik
		type(t_var) :: Vvar(1)
		call update()
		write(*,*)q
		bt=1d0/Tk
		iomg=(omg(1)+[0:l-1]*(omg(2)-omg(1))/l)+img*gm
		Vrpa=0d0
		Xq=0d0
		spin=1
		call set_var(Vvar(1),tp=3,nb=1)
		Vvar(1)%bd=0.5d0
		Vvar(1)%val=0.34d0*DJ
		call Hamilton(Vvar,tmp,q)
		Vrpa=real(tmp)*brizon%nq
		spin=2

		associate(nq=>brizon%nq,sb=>latt%sb,Ni=>latt%Ni)
		!$OMP PARALLEL DO REDUCTION(+:Xq) PRIVATE(Uk,Ek,Ukq,Ekq,fk,fkq)
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
		Xq=Xq/size(brizon%k,1)
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
	open(unit=20,file=fn('../data/band.dat'))
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

	call mkl_set_num_threads(32)
	call omp_set_num_threads(mkl_get_max_threads())

	!!nf=1d0-0.165d0
	!!Tk=0.02d0
	!var([icp,isc,iubo])%val(1)=[-1.4781d-01,8.8716d-06,1.9090d-01]
	!call Hamiltons(H,[3.141592654d+00,4.908738521d-02,0d0])
	!!call mat_diag(H,E)
	!call heev(H,E,"V")
	!write(*,*)E
	!stop

	Tk=0.0001d0
	n(1:4)=[0.23d0,0.24d0,0.245d0,0.25d0]

	nf=1d0-0.12d0
	!call selfconsist()
	var([icp,isc,iubo])%val(1)=[-1.5799d-01,1.2158d-01,1.8820d-01]
	!var([icp,iddw,iubo])%val(1)=[-1.5799d-01,1.2158d-01,1.8820d-01]
	!var([icp,iddw,iubo])%val(1)=[-1.5799d-01,0d0,1.8820d-01]
	!var([icp,iddw,iubo])%val(1)=[-1.8795d-01,1.2178d-01,1.9021d-01]
	!var([icp,isc,iubo])%val(1)=[-1.083d0,0.168d0,0d0]*0.12d0
	write(*,"(*(es12.4))")Tk,nf,var([icp,isc,iubo])%val(1),var(iddw)%val(1)
	!write(*,*)icp,iddw,iubo,isc
	!write(*,"(*(es12.4))")Tk,nf,var([icp,iddw,iubo])%val(1)

	!call update()
	!call band(20,[0d0,0d0,0d0],brizon%Ta(1,:),512)
	!call band(20,brizon%Ta(1,:),0.5d0*(brizon%Ta(1,:)+brizon%Ta(2,:)),512)
	!call band(20,0.5d0*(brizon%Ta(1,:)+brizon%Ta(2,:)),[0d0,0d0,0d0],512)
	!stop

	j=nint(sqrt(real(size(brizon%k,1))))
	!do i=0,j/5
		!!do j=0,0
		!call rpainstable_spin(Tk,[pi-2d0*pi*i/j,pi-2d0*pi*i/j,0d0],den)
		!write(*,"(*(es12.4))")1d0-nf,[pi,pi-2d0*pi*i/j],den
		!write(50,"(*(es12.4))")1d0-nf,[pi,pi-2d0*pi*i/j],den
	!enddo
	!stop
	!do i=1,size(brizon%k,1)
	!do i=128,256,8
		!do j=128,256,8
			!!if(brizon%k(i,1)>0.2d0*2d0*pi.and.brizon%k(i,2)>0.2d0*2d0*pi) then
				!!call rpa_EDC(60,brizon%k(i,:),0.01d0,0.05d0)
			!!endif
			!call rpa_EDC(60,[2d0*pi/512d0*i,2d0*pi/512d0*j],0.01d0,0.3d0)
		!enddo
	!enddo
	!stop
	do i=j/4,3*j/4
		if(i<j/2) then
			!call rpa_spin(70,[2d0*pi*i/j,2d0*pi*i/j,0d0],0.002d0,[0d0,0.5d0],1000)
			call rpa_spin_i(70,[2d0*pi*i/j,2d0*pi*i/j,0d0],0.02d0,[0d0,0.5d0],500)
		else
			!call rpa_spin(70,[pi-2d0*pi*(i-j/2)/j,pi,0d0],0.002d0,[0d0,0.5d0],1000)
			call rpa_spin_i(70,[pi-2d0*pi*(i-j/2)/j,pi,0d0],0.02d0,[0d0,0.5d0],500)
		endif
	enddo
	!call rpa_spin(70,[pi,pi,0d0],0.02d0,[0d0,0.5d0],500)
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
