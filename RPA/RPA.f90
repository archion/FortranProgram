program RPA
	implicit none
	complex(8),parameter :: img=(0d0,1d0)
	integer,parameter :: mk=512,mq=128,ms=100
	real(8),parameter :: t1=0.07415d0,t2=-0.0175d0,t3=0.0116d0,DJ=0.4d0,ph=0.1464d0,al=0.5d0,V=0.25d0,pi=3.1415926d0,n0=0.9,cvg=1e-4
	complex(8) :: H(4,4),Uk(4,4),Ukq(4,4),WORK(20),Xq(2,2,-mq:mq,-mq:mq),eka,dk
	integer :: i,j,k,l,m,n,o,p,INFO
	real(8) :: kx,ky,qx,qy,w,ek(4),ekq(4),fk(4),fkq(4),RWORK(10),Jq(-mq:mq,-mq:mq),&
			eks,bt,sp,sa,sb,wide,n1,Xq_RPA(-mq:mq,-mq:mq)
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
		do i=-mk,mk
			kx=pi/mk*i
			do j=-mk,mk
				ky=pi/mk*j
				eka=-2d0*t1*(1d0+DJ*30/8)*(cos(kx)*exp(img*ph)+cos(ky)*exp(-img*ph))
				eks=-4d0*t2*cos(kx)*cos(ky)-2d0*t3*(cos(2*kx)+cos(2*ky))-sp
				H=reshape((/ dcmplx(real(eka)+eks,0d0),    dcmplx(0d0,dimag(eka)),         dk,             (0d0,0d0),     &
							 dcmplx(0d0,-dimag(eka)),      dcmplx(eks-real(eka),0d0),     (0d0,0d0),          -dk,        &
							 dconjg(dk),            (0d0,0d0),         dcmplx(-eks-real(eka),0d0),    dcmplx(0d0,dimag(eka)),&
						     (0d0,0d0),       dconjg(-dk),         dcmplx(0d0,-dimag(eka)),  dcmplx(-eks+real(eka),0d0) /),(/4,4/))
				call ZHEEV( "V", "U", 4, H, 4, ek, WORK, 20, RWORK, INFO )
				fk=1d0/(1d0+exp(bt*ek))
				Uk=H
				if(flag1) then
					n1=n1+dot_product(Uk(1,:)*dconjg(Uk(1,:))-Uk(3,:)*dconjg(Uk(3,:)),fk)
				else
					if(flag0) then
						do o=1,3*ms-1
							kx=pi*(min((o/ms)*ms,ms)+(1-o/ms)*mod(o,ms))/ms
							ky=pi*min(max((o-ms),0),3*ms-o)/ms
							eka=-2d0*t1*(1d0+DJ*30/8)*(cos(kx)*exp(img*ph)+cos(ky)*exp(-img*ph))
							eks=-4d0*t2*cos(kx)*cos(ky)-2d0*t3*(cos(2*kx)+cos(2*ky))-sp
							H=reshape((/ dcmplx(real(eka)+eks,0d0), dcmplx(0d0,dimag(eka)),    dk,     (0d0,0d0),        &
										 dcmplx(0d0,-dimag(eka)),     dcmplx(eks-real(eka),0d0),  (0d0,0d0),          -dk,   &
										 dconjg(dk),  (0d0,0d0), dcmplx(-eks-real(eka),0d0), dcmplx(0d0,dimag(eka)),&
										 (0d0,0d0), dconjg(-dk), dcmplx(0d0,-dimag(eka)),  dcmplx(-eks+real(eka),0d0) /),(/4,4/))
							call ZHEEV( "V", "U", 4, H, 4, ek, WORK, 20, RWORK, INFO )
							write(10,"(8e16.3)")(ek(l),abs(H(1,l)),l=1,4)
						enddo
					endif
					flag0=.false.
					do k=-mq,mq
						qx=pi*(1d0+0.5d0*k/mq)
						do l=-mq,mq
							qy=pi*(1d0+0.5d0*l/mq)
							Jq(k,l)=al*DJ*(cos(qx)+cos(qy))
							kx=kx+qx
							ky=ky+qy
							eka=-2d0*t1*(1d0+DJ*30/8)*(cos(kx)+cos(ky))*exp(img*ph)
							eks=-4d0*t2*cos(kx)*cos(ky)-2d0*t3*(cos(2*kx)+cos(2*ky))-sp
							H=reshape((/ dcmplx(real(eka)+eks,0d0), dcmplx(0d0,dimag(eka)),    dk,     (0d0,0d0),        &
										 dcmplx(0d0,-dimag(eka)),     dcmplx(eks-real(eka),0d0),  (0d0,0d0),          -dk,   &
										 dconjg(dk),  (0d0,0d0), dcmplx(-eks-real(eka),0d0), dcmplx(0d0,dimag(eka)),&
										 (0d0,0d0), dconjg(-dk), dcmplx(0d0,-dimag(eka)),  dcmplx(-eks+real(eka),0d0) /),(/4,4/))
							call ZHEEV( "V", "U", 4, H, 4, ekq, WORK, 20, RWORK, INFO )
							fkq=1d0/(1d0+exp(bt*ekq))
							Ukq=H
							do p=1,1
								w=0.05d0
								do n=1,4
									do m=1,4
										Xq(1,1,k,l)=Xq(1,1,k,l)+&
											(-Uk(1,n)*Ukq(3,m)*dconjg(Uk(1,n)*Ukq(3,m))-Uk(1,n)*Ukq(3,m)*dconjg(Uk(2,n)*Ukq(4,m))+&
											Uk(1,n)*Ukq(3,m)*dconjg(Uk(3,n)*Ukq(1,m))+Uk(1,n)*Ukq(3,m)*dconjg(Uk(4,n)*Ukq(2,m)))
										Xq(1,2,k,l)=Xq(1,2,k,l)+&
											(-Ukq(1,m)*Uk(1,n)*dconjg(Ukq(1,m)*Uk(2,n))-Uk(1,n)*Ukq(3,m)*dconjg(Uk(1,n)*Ukq(4,m))+&
											Uk(1,n)*Ukq(3,m)*dconjg(Uk(3,n)*Ukq(2,m))+Uk(1,n)*Ukq(3,m)*dconjg(Uk(4,n)*Ukq(1,m)))
										Xq(2,1,k,l)=Xq(2,1,k,l)+&
											(-Ukq(2,m)*Uk(1,n)*dconjg(Ukq(2,m)*Uk(2,n))-Uk(1,n)*Ukq(4,m)*dconjg(Uk(1,n)*Ukq(3,m))+&
											Uk(1,n)*Ukq(4,m)*dconjg(Uk(4,n)*Ukq(2,m))+Uk(1,n)*Ukq(4,m)*dconjg(Uk(3,n)*Ukq(1,m)))
										Xq(2,2,k,l)=Xq(2,2,k,l)+&
											(-Ukq(4,m)*Uk(1,n)*dconjg(Ukq(4,m)*Uk(1,n))-Uk(1,n)*Ukq(4,m)*dconjg(Uk(2,n)*Ukq(3,m))+&
											Uk(1,n)*Ukq(4,m)*dconjg(Uk(3,n)*Ukq(2,m))+Uk(1,n)*Ukq(4,m)*dconjg(Uk(4,n)*Ukq(1,m)))
										Xq(:,:,k,l)=(1-fk(n)-fkq(m))/(w-ek(n)-ekq(m)+img*0.1d0)*Xq(:,:,k,l)
									enddo
								enddo
							enddo
						enddo
					enddo
				endif
			enddo
		enddo
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
	Xq_RPA=dimag((Xq(1,1,:,:)-Jq*(Xq(1,1,:,:)*Xq(2,2,:,:)-Xq(1,2,:,:)*Xq(2,1,:,:)))/&
			(1d0+Jq*(Xq(1,1,:,:)-Xq(2,2,:,:))-Jq**2*(Xq(1,1,:,:)*Xq(2,2,:,:)-Xq(1,2,:,:)*Xq(2,1,:,:))))
	do k=-mq,mq
		do l=-mq,mq
			write(10,"(3e16.3)")pi*(1d0+0.05d0*k/mq),pi*(1d0+0.05d0*l/mq),Xq_RPA(k,l)
		enddo
		write(10,"(1x)")
	enddo
end