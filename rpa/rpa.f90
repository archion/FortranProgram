program rpa
	implicit none
	complex(8),parameter :: img=(0d0,1d0)
	integer,parameter :: mk=128,mq=32,ms=100
	real(8),parameter :: t1=0.07415d0,t2=-0.0175d0,t3=0.0116d0,DJ=0.4d0,ph=0.1464d0,al=0.5d0,V=0.25d0,pi=3.1415926d0,n0=0.9,cvg=1e-4
	complex(8) :: H(4,4),Uk(4,4),Ukq(4,4),UEk(4),WORK(20),Xq(2,2,-mq:mq,-mq:mq),Xq_tmp(2,2),eka,cth,sth,cfy(2),sfy(2)
	integer :: i,j,k,l,m,n,o,p,INFO
	real(8) :: kx,ky,qx,qy,kqx,kqy,w,e1(2),e2(2),ek(4),ekq(4),fk(4),fkq(4),RWORK(10),Jq,&
		eks,bt,sp,sa,sb,wide,n1,Xq_RPA,u2(2),v2(2),cs2(2),th,fy(2),dk
	logical :: flag0,flag1,flaga,flagb
	open(unit=10,file="../data/rpa.dat")
	open(unit=20,file="../data/check.dat")
	flag0=.true.
	flag1=.true.
	flaga=.true.
	flagb=.true.
	sp=0d0
	sa=sp
	sb=sp
	wide=0.1d0
	dk=0d0
	bt=1e5
	do while(flag0)
		sp=(sa+sb)*0.5d0
		n1=0d0
		Xq=0d0
		!$OMP PARALLEL DO REDUCTION(+:n1,Xq) PRIVATE(kx,ky,eka,eks,e1,e2,th,fy,cth,sth,cfy,sfy,fk,UK,H,ek,ekq,work,rwork,info,&
		!$OMP qx,qy,kqx,kqy,fkq,Ukq,w,Xq_tmp) SCHEDULE(GUIDED)
		do i=-mk,mk
			kx=pi/mk*i
			do j=-mk,mk
				ky=pi/mk*j
				eka=-2d0*t1*(1d0+DJ*30/8)*(cos(kx)*exp(img*ph)+cos(ky)*exp(-img*ph))
				eks=-4d0*t2*cos(kx)*cos(ky)-2d0*t3*(cos(2*kx)+cos(2*ky))-sp
				e1=eks+(/1d0,-1d0/)*sqrt(dimag(eka)**2+real(eka)**2)
				e2=sqrt(e1**2+dk**2)
				ek=(/e2,-e2/)
				if(dimag(eka)>0) then
					th=acos(real(eka)/sqrt(dimag(eka)**2+real(eka)**2))/2d0
				else
					th=(2d0*pi-acos(real(eka)/sqrt(dimag(eka)**2+real(eka)**2)))/2d0
				endif
				if(dk>0) then
					fy=acos(e1/e2)/2d0
				else
					fy=(2d0*pi-acos(e1/e2))/2d0
				endif
				cth=dcmplx(cos(th))
				sth=dcmplx(sin(th))
				cfy=dcmplx(cos(fy))
				sfy=dcmplx(sin(fy))
				Uk=reshape((/cth*cfy(1),-img*sth*cfy(1),cth*sfy(1),img*sth*sfy(1),&
							-img*sth*cfy(2),cth*cfy(2),-img*sth*sfy(2),cth*sfy(2),&
							-cth*sfy(1),img*sth*sfy(1),cth*cfy(1),img*sth*cfy(1),&
							-img*sth*sfy(2),cth*sfy(2),img*sth*cfy(2),cth*cfy(2) /)&
							,(/4,4/))
				!H=reshape((/ dcmplx(real(eka)+eks,0d0),    dcmplx(0d0,dimag(eka)),         dk,             (0d0,0d0),     &
					!dcmplx(0d0,-dimag(eka)),      dcmplx(eks-real(eka),0d0),     (0d0,0d0),          -dk,        &
					!dconjg(dk),            (0d0,0d0),         dcmplx(-eks-real(eka),0d0),    dcmplx(0d0,dimag(eka)),&
					!(0d0,0d0),       dconjg(-dk),         dcmplx(0d0,-dimag(eka)),  dcmplx(-eks+real(eka),0d0) /),(/4,4/))
				!call ZHEEV( "V", "U", 4, H, 4, ek, WORK, 20, RWORK, INFO )
				! e1=eks+(/-1d0,1d0/)*abs(eka)
				! e2=sqrt(e1**2+dk**2)
				! ek=(/e2,-e2/)
				! uk2=(1d0+e1/e2)*0.5d0
				! vk2=(1d0-e1/e2)*0.5d0
				! csk2=(1d0+(/1d0,-1d0/)*real(eka)/abs(eka))*0.5d0
				! H=(/(u2-v2)*cs2,(v2-u2)*cs2/)
				fk=1d0/(1d0+exp(bt*ek))
				!Uk=H
				if(flag1) then
					! n1=n1+dot_product(H,fk)
					n1=n1+dot_product(Uk(1,:)*dconjg(Uk(1,:))-Uk(3,:)*dconjg(Uk(3,:)),fk)
				else
					! uvk=0.5d0*dk/e2
					! csk=(/-img,img/)0.5d0*dimag(eka)/abs(eka)
					!$OMP CRITICAL
					if(flag0) then
						do o=1,3*ms-1
							kx=pi*(min((o/ms)*ms,ms)+(1-o/ms)*mod(o,ms))/ms
							ky=pi*min(max((o-ms),0),3*ms-o)/ms
							eka=-2d0*t1*(1d0+DJ*30/8)*(cos(kx)*exp(img*ph)+cos(ky)*exp(-img*ph))
							eks=-4d0*t2*cos(kx)*cos(ky)-2d0*t3*(cos(2*kx)+cos(2*ky))-sp
							e1=eks+(/1d0,-1d0/)*sqrt(dimag(eka)**2+real(eka)**2)
							e2=sqrt(e1**2+dk**2)
							ek=(/e2,-e2/)
							if(dimag(eka)>0) then
								th=acos(real(eka)/sqrt(dimag(eka)**2+real(eka)**2))/2d0
							else
								th=(2d0*pi-acos(real(eka)/sqrt(dimag(eka)**2+real(eka)**2)))/2d0
							endif
							if(dk>0) then
								fy=acos(e1/e2)/2d0
							else
								fy=(2d0*pi-acos(e1/e2))/2d0
							endif
							cth=dcmplx(cos(th))
							sth=dcmplx(sin(th))
							cfy=dcmplx(cos(fy))
							sfy=dcmplx(sin(fy))
							UEk=(/cth*cfy(1), -img*sth*cfy(2),-cth*sfy(1),-img*sth*sfy(2)/)
							! u2=(1d0+e1/e2)*0.5d0
							! v2=(1d0-e1/e2)*0.5d0
							! cs2=(1d0+(/1d0,-1d0/)*real(eka)/sqrt(dimag(eka)**2+real(eka)**2))*0.5d0
							! H=(/u2*cs2,v2*cs2/)
							! write(10,"(8e16.3)")(ek(l),real(H(l)),l=1,4)
							!H=reshape((/ dcmplx(real(eka)+eks,0d0), dcmplx(0d0,dimag(eka)),    dk,     (0d0,0d0),        &
								!dcmplx(0d0,-dimag(eka)),     dcmplx(eks-real(eka),0d0),  (0d0,0d0),          -dk,   &
								!dconjg(dk),  (0d0,0d0), dcmplx(-eks-real(eka),0d0), dcmplx(0d0,dimag(eka)),&
								!(0d0,0d0), dconjg(-dk), dcmplx(0d0,-dimag(eka)),  dcmplx(-eks+real(eka),0d0) /),(/4,4/))
							!call ZHEEV( "V", "U", 4, H, 4, ek, WORK, 20, RWORK, INFO )
							write(20,"(8e16.3)")(ek(l),abs(UEk(l)),l=1,4)
						enddo
					endif
					flag0=.false.
					!$OMP END CRITICAL
					do k=-mq,mq
						qx=pi*(1d0+0.5d0*k/mq)
						do l=-mq,mq
							qy=pi*(1d0+0.5d0*l/mq)
							kqx=kx+qx
							kqy=ky+qy
							eka=-2d0*t1*(1d0+DJ*30/8)*(cos(kqx)+cos(kqy))*exp(img*ph)
							eks=-4d0*t2*cos(kqx)*cos(kqy)-2d0*t3*(cos(2*kqx)+cos(2*kqy))-sp
							e1=eks+(/1d0,-1d0/)*sqrt(dimag(eka)**2+real(eka)**2)
							e2=sqrt(e1**2+dk**2)
							ekq=(/e2,-e2/)
							if(dimag(eka)>0) then
								th=acos(real(eka)/sqrt(dimag(eka)**2+real(eka)**2))/2d0
							else
								th=(2d0*pi-acos(real(eka)/sqrt(dimag(eka)**2+real(eka)**2)))/2d0
							endif
							if(dk>0) then
								fy=acos(e1/e2)/2d0
							else
								fy=(2d0*pi-acos(e1/e2))/2d0
							endif
							cth=dcmplx(cos(th))
							sth=dcmplx(sin(th))
							cfy=dcmplx(cos(fy))
							sfy=dcmplx(sin(fy))
							Ukq=reshape((/cth*cfy(1),-img*sth*cfy(1),cth*sfy(1),img*sth*sfy(1),&
										-img*sth*cfy(2),cth*cfy(2),-img*sth*sfy(2),cth*sfy(2),&
										-cth*sfy(1),img*sth*sfy(1),cth*cfy(1),img*sth*cfy(1),&
										-img*sth*sfy(2),cth*sfy(2),img*sth*cfy(2),cth*cfy(2) /)&
										,(/4,4/))
						   ! Ukq=reshape((/dcmplx(cos(th)*cos(fy(1))),-img*sin(th)*cos(fy(1)),dcmplx(cos(th)*sin(fy(1))),&
								!img*sin(th)*sin(fy(1)),&
										!-img*sin(th)*cos(fy(2)),dcmplx(cos(th)*cos(fy(2))),-img*sin(th)*sin(fy(2)),&
								!dcmplx(-cos(th)*sin(fy(2))),&
										!dcmplx(-cos(th)*sin(fy(1))),img*sin(th)*sin(fy(1)),dcmplx(cos(th)*cos(fy(1))),&
								!img*sin(th)*cos(fy(1)),&
										!-img*sin(th)*sin(fy(2)),dcmplx(cos(th)*sin(fy(2))),img*sin(th)*cos(fy(2)),&
								!dcmplx(cos(th)*cos(fy(2))) /),(/4,4/))
							! e1=eks+(/-1d0,1d0/)*abs(eka)
							! e2=sqrt(e1**2+dk**2)
							! ekq=(/e2,-e2/)
							! ukq2=(1d0+e1/e2)*0.5d0
							! vkq2=(1d0-e1/e2)*0.5d0
							! cskq2=(1d0+(/1d0,-1d0/)*real(eka)/abs(eka))*0.5d0
							! uvkq=0.5d0*dk/e2
							! cskq=(/-img,img/)0.5d0*dimag(eka)/abs(eka)
							!H=reshape((/ dcmplx(real(eka)+eks,0d0), dcmplx(0d0,dimag(eka)),    dk,     (0d0,0d0),        &
								!dcmplx(0d0,-dimag(eka)),     dcmplx(eks-real(eka),0d0),  (0d0,0d0),          -dk,   &
								!dconjg(dk),  (0d0,0d0), dcmplx(-eks-real(eka),0d0), dcmplx(0d0,dimag(eka)),&
								!(0d0,0d0), dconjg(-dk), dcmplx(0d0,-dimag(eka)),  dcmplx(-eks+real(eka),0d0) /),(/4,4/))
							!call ZHEEV( "V", "U", 4, H, 4, ekq, WORK, 20, RWORK, INFO )
							fkq=1d0/(1d0+exp(bt*ekq))
							!Ukq=H
							do p=1,1
								w=0.05d0
								do n=1,4
									do m=1,4
										Xq_tmp(1,1)=&
											(Uk(1,n)*Ukq(3,m)*dconjg(Uk(1,n)*Ukq(3,m))+Uk(1,n)*Ukq(3,m)*dconjg(Uk(2,n)*Ukq(4,m))-&
											Uk(1,n)*Ukq(3,m)*dconjg(Uk(3,n)*Ukq(1,m))-Uk(1,n)*Ukq(3,m)*dconjg(Uk(4,n)*Ukq(2,m)))
										Xq_tmp(1,2)=&
											(Ukq(3,m)*Uk(1,n)*dconjg(Ukq(3,m)*Uk(2,n))+Uk(1,n)*Ukq(3,m)*dconjg(Uk(1,n)*Ukq(4,m))-&
											Uk(1,n)*Ukq(3,m)*dconjg(Uk(3,n)*Ukq(2,m))-Uk(1,n)*Ukq(3,m)*dconjg(Uk(4,n)*Ukq(1,m)))
										Xq_tmp(2,1)=dconjg(Xq_tmp(1,2))
										Xq_tmp(2,2)=&
											(Ukq(4,m)*Uk(1,n)*dconjg(Ukq(4,m)*Uk(1,n))+Uk(1,n)*Ukq(4,m)*dconjg(Uk(2,n)*Ukq(3,m))-&
											Uk(1,n)*Ukq(4,m)*dconjg(Uk(3,n)*Ukq(2,m))-Uk(1,n)*Ukq(4,m)*dconjg(Uk(4,n)*Ukq(1,m)))
										Xq(:,:,k,l)=Xq(:,:,k,l)+(1-fk(n)-fkq(m))/(w-ek(n)-ekq(m)+img*0.005d0)*Xq_tmp
									enddo
								enddo
							enddo
						enddo
					enddo
				endif
			enddo
		enddo
		!$OMP END PARALLEL DO
		n1=n1/((2*mk)**2)+1d0
		write(*,*)sp,n1
		if(abs(n1-n0)<=cvg) then
			flag1=.false.
		else
			if(n1<n0) then
				flaga=.false.
				sa=sp
				if(flaga.or.flagb) then
					sb=sp+wide
					sp=sb
				endif
			else
				flagb=.false.
				sb=sp
				if(flaga.or.flagb) then
					sa=sp-wide
					sp=sa
				endif
			endif
		endif
	enddo
	Xq=-Xq/((2*mk)**2)
	! Xq_RPA=dimag((Xq(1,1,:,:)+Jq*(Xq(1,1,:,:)*Xq(2,2,:,:)-Xq(1,2,:,:)*Xq(2,1,:,:)))/&
	! (1d0-Jq*(Xq(1,1,:,:)-Xq(2,2,:,:))-Jq**2*(Xq(1,1,:,:)*Xq(2,2,:,:)-Xq(1,2,:,:)*Xq(2,1,:,:))))
	! Xq_RPA=dimag((Xq(1,1,:,:)-Jq*(Xq(1,1,:,:)*Xq(2,2,:,:)+Xq(1,2,:,:)*Xq(2,1,:,:)))/&
	! (1d0+Jq*(Xq(1,1,:,:)-Xq(2,2,:,:))-Jq**2*(Xq(1,1,:,:)*Xq(2,2,:,:)-Xq(1,2,:,:)*Xq(2,1,:,:))))
	do k=-mq,mq
		qx=pi*(1d0+0.5d0*k/mq)
		do l=-mq,mq
			qy=pi*(1d0+0.5d0*l/mq)
			Jq=al*DJ*(cos(qx)+cos(qy))
			Xq_RPA=dimag((Xq(1,1,k,l)-Jq*(Xq(1,1,k,l)*Xq(2,2,k,l)+Xq(1,2,k,l)*Xq(2,1,k,l)))/&
				(1d0+Jq*(Xq(1,1,k,l)-Xq(2,2,k,l))-Jq**2*(Xq(1,1,k,l)*Xq(2,2,k,l)-Xq(1,2,k,l)*Xq(2,1,k,l))))
			write(10,"(4e16.3)")qx,qy,Xq_RPA,dimag(Xq(1,1,k,l))
		enddo
		write(10,"(1x)")
	enddo
end
