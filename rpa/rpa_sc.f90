module global
	implicit none
	save
	real(8), parameter ::t(5)=(/2d0,-0.9d0,0.d0,0.d0,0d0/),escal=0.25d0,pi=3.1415926d0,cvg=1d-5,DJ=0.75d0,V=-3d0,alp=0.34d0
	real(8) :: nf=0.8d0,Tk,it=0.005d0,omgrg(2)=(/-2d0,2d0/),qrg(2)=(/0d0,2d0/)*pi
	character(3) :: pgflag="ddw"
	!character(3) :: pgflag="sdw"
	integer, parameter :: mk=64,mq=64,mo=512,mro=512
	complex(8),parameter :: img=(0d0,1d0)
	contains
	subroutine spinsus(pg,sc,ap,sp,q,omg,Xq)
		implicit none
		integer :: i,j,l,m,n
		complex(8) :: Uk(4,4),Ukq(4,4),omg(:),Xq(:,:,:),Xtmp1(size(Xq,1),2,2),Xtmp2(size(Xq,1),2,2),Det(size(Xq,1))
		real(8) :: k(2),q(2),sp,bt,pg,sc,ap,ek(4),ekq(4),fk(4),fkq(4),DJq,DJqq
		bt=1d0/Tk
		Xq=0d0
		!$OMP PARALLEL DO REDUCTION(+:Xq) PRIVATE(k,ek,Uk,fk,ekq,Ukq,fkq) SCHEDULE(GUIDED)
		do i=1,mk
			k(1)=2d0*pi/mk*i
			do j=1,mk
				k(2)=2d0*pi/mk*j
				call EU(k,pg,sc,ap,sp,ek,Uk)
				call EU(k+q,pg,sc,ap,sp,ekq,Ukq)
				fk=1d0/(1d0+exp(bt*ek))
				fkq=1d0/(1d0+exp(bt*ekq))
				do n=1,4
					do m=1,4
						Xq(:,1,1)=Xq(:,1,1)+(Uk(1,n)*Ukq(3,m)*dconjg(Uk(1,n)*Ukq(3,m))+Uk(1,n)*Ukq(3,m)*dconjg(Uk(2,n)*Ukq(4,m))-&
							Uk(1,n)*Ukq(3,m)*dconjg(Uk(3,n)*Ukq(1,m))-Uk(1,n)*Ukq(3,m)*dconjg(Uk(4,n)*Ukq(2,m)))*&
							(1-fk(n)-fkq(m))/(omg-ek(n)-ekq(m)+1d-10)
						Xq(:,1,2)=Xq(:,1,2)+(Uk(1,n)*Ukq(3,m)*dconjg(Uk(1,n)*Ukq(4,m))+Uk(1,n)*Ukq(3,m)*dconjg(Uk(2,n)*Ukq(3,m))-&
							Uk(1,n)*Ukq(3,m)*dconjg(Uk(3,n)*Ukq(2,m))-Uk(1,n)*Ukq(3,m)*dconjg(Uk(4,n)*Ukq(1,m)))*&
							(1-fk(n)-fkq(m))/(omg-ek(n)-ekq(m)+1d-10)
						Xq(:,2,1)=Xq(:,2,1)+(Uk(1,n)*Ukq(4,m)*dconjg(Uk(1,n)*Ukq(3,m))+Uk(1,n)*Ukq(4,m)*dconjg(Uk(2,n)*Ukq(4,m))-&
							Uk(1,n)*Ukq(4,m)*dconjg(Uk(3,n)*Ukq(1,m))-Uk(1,n)*Ukq(4,m)*dconjg(Uk(4,n)*Ukq(2,m)))*&
							(1-fk(n)-fkq(m))/(omg-ek(n)-ekq(m)+1d-10)
						Xq(:,2,2)=Xq(:,2,2)+(Uk(1,n)*Ukq(4,m)*dconjg(Uk(1,n)*Ukq(4,m))+Uk(1,n)*Ukq(4,m)*dconjg(Uk(2,n)*Ukq(3,m))-&
							Uk(1,n)*Ukq(4,m)*dconjg(Uk(3,n)*Ukq(2,m))-Uk(1,n)*Ukq(4,m)*dconjg(Uk(4,n)*Ukq(1,m)))*&
							(1-fk(n)-fkq(m))/(omg-ek(n)-ekq(m)+1d-10)
					enddo
				enddo
			enddo
		enddo
		!$OMP END PARALLEL DO
		Xq=-Xq/mk**2
		DJq=alp*(cos(q(1))+cos(q(2)))
		DJqq=-alp*(cos(q(1))+cos(q(2)))
		!DJq=-DJ*2d0
		!DJqq=-DJ*2d0
		Xtmp1=Xq
		Xtmp2(:,1,1)=1d0+DJqq*Xtmp1(:,2,2)
		Xtmp2(:,2,2)=1d0+DJq*Xtmp1(:,1,1)
		Xtmp2(:,1,2)=-DJq*Xtmp1(:,1,2)
		Xtmp2(:,2,1)=-DJqq*Xtmp1(:,2,1)
		Det=1d0/(Xtmp2(:,1,1)*Xtmp2(:,2,2)-Xtmp2(:,1,2)*Xtmp2(:,2,1))
		Xq(:,1,1)=(Xtmp1(:,1,1)*Xtmp2(:,1,1)+Xtmp1(:,1,2)*Xtmp2(:,2,1))*Det
		Xq(:,1,2)=(Xtmp1(:,1,1)*Xtmp2(:,1,2)+Xtmp1(:,1,2)*Xtmp2(:,2,2))*Det
		Xq(:,2,1)=(Xtmp1(:,2,1)*Xtmp2(:,1,1)+Xtmp1(:,2,2)*Xtmp2(:,2,1))*Det
		Xq(:,2,2)=(Xtmp1(:,2,1)*Xtmp2(:,1,2)+Xtmp1(:,2,2)*Xtmp2(:,2,2))*Det
	end subroutine
	subroutine greenfunc(pg,sc,ap,sp,omg,G0,G0_inv)
		implicit none
		integer :: i,j,l,m,n,o,p,ii,jj,sig
		complex(8) :: Uk(4,4),omg(:),G0(:,0:,0:,:,:),G0_inv(:,0:,0:,:,:)
		real(8) :: pg,sc,ap,sp,ek(4),k(2)
		G0=0d0
		G0_inv=0d0
		l=size(G0,4)
		o=size(G0,2)
		!$OMP PARALLEL DO PRIVATE(k,ek,Uk) SCHEDULE(GUIDED)
		do i=0,o/2
			do j=0,i
				k=(/i,j/)*2d0*pi/o
				call EU(k,pg,sc,ap,sp,ek,Uk)
				do m=1,l
					do n=1,l
						do p=1,l
							!G0(:,i,j,m,n)=G0(:,i,j,m,n)-(-1)**((m-1)/2+(n-1)/2)*Uk(mod(n+1,4)+1,p)*dconjg(Uk(mod(m+1,4)+1,p))/&
								!(-omg-ek(p))
							G0(:,i,j,m,n)=G0(:,i,j,m,n)+Uk(m,p)*dconjg(Uk(n,p))/(omg-ek(p)+1d-10)
							G0_inv(:,i,j,m,n)=G0_inv(:,i,j,m,n)+Uk(m,p)*dconjg(Uk(n,p))*(omg-ek(p))
						enddo
					enddo
				enddo
			enddo
		enddo
		!$OMP END PARALLEL DO
		do i=1,l
			do j=1,l
				call symfill(G0(:,:,:,i,j),(-1)**((i-1)/2+(j-1)/2))
				call symfill(G0_inv(:,:,:,i,j),(-1)**((i-1)/2+(j-1)/2))
			enddo
		enddo
	end subroutine
	subroutine effint(pg,sc,ap,sp,omg,Veff)
		use fft, only : fft3d
		implicit none
		integer :: i,j,n,ii,jj
		real(8) :: pg,sc,ap,sp,DJq,q(2),dq
		complex(8) :: Veff(:,0:,0:,:,:),omg(:)
		n=size(Veff,2)
		dq=(qrg(2)-qrg(1))/n
		!!$OMP PARALLEL DO PRIVATE(q,Veff,DJq) SCHEDULE(GUIDED)
		do i=0,n/2
			do j=0,i
				!q=((/i,j/)*2d0*pi/size(Veff,2)-(/pi,pi/))*0.5+(/pi,pi/)
				q=qrg(1)+(/i,j/)*dq
				call spinsus(pg,sc,ap,sp,q,omg,Veff(:,i,j,:,:))
				DJq=(cos(q(1))+cos(q(2)))
				Veff(:,i,j,1,1)=Veff(:,i,j,1,1)*DJq**2
				Veff(:,i,j,1,2)=-Veff(:,i,j,1,2)*DJq**2
				Veff(:,i,j,2,1)=-Veff(:,i,j,2,1)*DJq**2
				Veff(:,i,j,2,2)=Veff(:,i,j,2,2)*DJq**2
			enddo
		enddo
		!!$OMP END PARALLEL DO
		do i=1,2
			do j=1,2
				call symfill(Veff(:,:,:,i,j),1)
			enddo
		enddo
	end subroutine
	subroutine symfill(A,sig)
		implicit none
		integer :: i,j,n,ii,jj,sig
		complex(8) :: A(:,0:,0:)
		n=size(A,2)
		!$OMP PARALLEL DO PRIVATE(ii,jj) SCHEDULE(GUIDED)
		do i=0,n/2
			do j=0,i
				ii=mod(n-i,n)
				jj=mod(n-j,n)
				A(: , j  , i ) = A(: , i , j) * sig
				A(: , jj , i ) = A(: , i , j) * sig
				A(: , j  , ii) = A(: , i , j) * sig
				A(: , jj , ii) = A(: , i , j) * sig
				A(: , ii , j ) = A(: , i , j)
				A(: , i  , jj) = A(: , i , j)
				A(: , ii , jj) = A(: , i , j)
			enddo
		enddo
		!$OMP END PARALLEL DO
	end subroutine
	subroutine selfe(G0,Veff,sfeg)
		use fft
		implicit none
		integer :: i,j,n
		complex(8) :: sfeg(:,:,:,:,:),Veff(:,:,:,:,:),G0(:,:,:,:,:),&
			Gtmp(size(G0,1),size(G0,2),size(G0,3),size(G0,4),size(G0,5))
		real(8) :: k(2),q(2),sp,bt,sc,ap,ek(2),ekq(2),fk(2),fkq(2),DJq
		n=size(G0,4)
		Gtmp=G0
		do i=1,2
			do j=1,2
				call fft3d(Veff(:,:,:,i,j),-1)
			enddo
		enddo
		do i=1,n
			do j=1,n
				call fft3d(Gtmp(:,:,:,i,j),-1)
				sfeg(:,:,:,i,j)=Gtmp(:,:,:,i,j)*Veff(:,:,:,1,1)+Gtmp(:,:,:,i,j-(-1)**j)*Veff(:,:,:,1,2)+&
					Gtmp(:,:,:,i-(-1)**i,j)*Veff(:,:,:,2,1)+Gtmp(:,:,:,i-(-1)**i,j-(-1)**j)*Veff(:,:,:,2,2)
				call fft3d(sfeg(:,:,:,i,j),1)
			enddo
		enddo
		sfeg=sfeg*Tk*size(Gtmp,1)
	end subroutine
	!subroutine rpainstable(sc,ap,sp,q)
		!implicit none
		!integer :: i,j,l,m,n
		!complex(8) :: Uk(2,2),Ukq(2,2)
		!real(8) :: k(2),q(2),sp,bt,pg,sc,ap,ek(2),ekq(2),fk(2),fkq(2),omg,omgr(2),domg,DJq,Xq,Xq_rpa
		!bt=1d0/Tk
		!Xq=0d0
		!!$OMP PARALLEL DO REDUCTION(+:Xq) PRIVATE(k,ek,Uk,fk,ekq,Ukq,fkq) SCHEDULE(GUIDED)
		!do i=0,mq
			!k(1)=pi/mq*i
			!do j=0,mq
				!k(2)=pi/mq*j
				!call EU(k,pg,sc,ap,sp,ek,Uk)
				!call EU(k+q,pg,sc,ap,sp,ekq,Ukq)
				!fk=1d0/(1d0+exp(bt*ek))
				!fkq=1d0/(1d0+exp(bt*ekq))
				!do n=1,2
					!do m=1,2
						!Xq=Xq+(Uk(1,n)*Ukq(2,m)*dconjg(Uk(1,n)*Ukq(2,m))-Uk(2,n)*Ukq(1,m)*dconjg(Uk(1,n)*Ukq(2,m)))*&
							!(1-fk(n)-fkq(m))/(ek(n)+ekq(m))
					!enddo
				!enddo
			!enddo
		!enddo
		!!$OMP END PARALLEL DO
		!Xq=Xq/mq**2
		!write(*,"(2e16.3)")1d0-nf,1d0+DJq*Xq
	!end subroutine
	subroutine selfconsist_tg(pg,sc,ap,sp,sig)
		implicit none
		complex(8) :: Uk(4,4)
		real(8) :: k(2),np,al,ek(4),sp,pg,pgp,sc,scp,ap,app,pnp
		integer :: i,j,c,info,sig
		logical :: flaga,flagb
		c=0
		al=1d0
		do 
			pnp=nf+0.1d0
			do 			
				scp=0d0
				pgp=0d0
				app=0d0
				np=0d0
				c=c+1
				!call sfcsap(sp,ap)
				!$OMP PARALLEL DO REDUCTION(+:np,pgp,scp,app) PRIVATE(k,ek,Uk) SCHEDULE(GUIDED)
				do i=0,mk
					do j=0,min(i,mk-i)
						k=(/pi/mk*i,pi/mk*j/)
						call EU(k,pg,sc,ap,sp,ek,Uk)
						call order(k,ek,Uk,np,pgp,scp,app)
					enddo
				enddo
				!$OMP END PARALLEL DO
				np=np/(mk**2)*2d0
				pgp=2d0*DJ*pgp/(mk**2)*2d0
				app=2d0*DJ*app/(mk**2)*2d0
				scp=V*scp/(mk**2)*2d0
				if((abs(nf-np)/abs(nf-pnp))>0.8d0) then
					al=max(al-0.05d0,0.05)
				endif
				if(abs(nf-np)<=cvg) then
					exit
				else
					sp=sp+al*(nf-np)
				endif
				pnp=np
			enddo
			if(sig/=0) then
				scp=sc
				app=ap
				pgp=pg
			endif
			if((abs(scp-sc)+abs(pgp-pg)+abs(ap-app))<cvg) then
				exit
			endif
			sc=scp
			ap=app
			pg=pgp
		enddo
		!write(*,*)"!!!!!!selfconsist return!!!!!!!!"
		write(*,"(6e12.4,i4)")np,Tk,sp,pg,sc,ap,c
		!write(*,*)"!!!!!!!!!!!!end!!!!!!!!!!!!!!!!!"
	end subroutine
	subroutine EU(k,pg,sc,ap,sp,ek,Uk)
		implicit none
		complex(8) :: Uk(4,4),cth,sth,cfy(2),sfy(2)
		real(8) :: eka,eks,e1(2),e2(2),ek(4),k(2),pg,pgk,sp,sck,sc,ap,gk,th,fy(2),cos2th,sin2th,cos2fy(2),sin2fy(2)
		eks=-4d0*(1d0-nf)*t(2)*cos(k(1))*cos(k(2))-sp-2d0*(1d0-nf)*t(3)*(cos(2d0*k(1))+cos(2d0*k(2)))
		eka=-2d0*((1d0-nf)*t(1)+ap)*(cos(k(1))+cos(k(2)))
		!eks=-4d0*t(2)*cos(k(1))*cos(k(2))-2d0*t(3)*(cos(2d0*k(1))+cos(2d0*k(2)))-sp
		!eka=-2d0*t(1)*(cos(k(1))+cos(k(2)))
		gk=0.5d0*(cos(k(1))-cos(k(2)))
		select case(pgflag)
		case("ddw")
			pgk=pg*gk
		case("sdw")
			pgk=-pg
		end select
		sck=sc*gk
		e1=eks+(/1d0,-1d0/)*sqrt(pgk**2+eka**2)*sign(1d0,eka)
		e2=sqrt(e1**2+sck**2)
		ek=(/e2*sign(1d0,e1),-e2*sign(1d0,e1)/)
		cos2th=eka/sqrt(pgk**2+eka**2)*sign(1d0,eka)
		sin2th=pgk/sqrt(pgk**2+eka**2)*sign(1d0,eka)
		cos2fy=e1/e2*sign(1d0,e1)
		sin2fy=sck/e2*sign(1d0,e1)
		cth=dcmplx(sqrt(0.5d0*(1d0+cos2th)))
		select case(pgflag)
		case("ddw")
			sth=img*dcmplx(sqrt(0.5d0*(1d0-cos2th))*sign(1d0,sin2th))
		case("sdw")
			sth=dcmplx(sqrt(0.5d0*(1d0-cos2th))*sign(1d0,sin2th))
		end select
		cfy=dcmplx(sqrt(0.5d0*(1d0+cos2fy)))
		sfy=dcmplx(sqrt(0.5d0*(1d0-cos2fy))*sign(1d0,sin2fy))
		Uk=reshape((/   cth*cfy(1) ,  dconjg(sth)*cfy(1) ,  cth*sfy(1) , -dconjg(sth)*sfy(1) ,   & 
					   -sth*cfy(2) ,          cth*cfy(2) , -sth*sfy(2) ,         -cth*sfy(2) ,   & 
					   -cth*sfy(1) , -dconjg(sth)*sfy(1) ,  cth*cfy(1) , -dconjg(sth)*cfy(1) ,   & 
					   -sth*sfy(2) ,          cth*sfy(2) ,  sth*cfy(2) ,   cth*cfy(2)     /) , (/4 , 4/))
	end subroutine
	subroutine order(k,ek,Uk,n,pg,sc,ap)
		implicit none
		complex(8) :: Uk(4,4)
		real(8) :: k(2),ek(4),pg,sc,gk,n,fk(4),ap
		fk=1d0/(1d0+exp(ek/Tk))
		gk=0.5d0*(cos(k(1))-cos(k(2)))
		n=n+2d0+dot_product(Uk(1,:)*dconjg(Uk(1,:))+Uk(2,:)*dconjg(Uk(2,:))-&
			Uk(3,:)*dconjg(Uk(3,:))-Uk(4,:)*dconjg(Uk(4,:)),fk)
		ap=ap+dot_product(Uk(1,:)*dconjg(Uk(1,:))-Uk(2,:)*dconjg(Uk(2,:)),fk)*cos(k(1))
		select case(pgflag)
		case("ddw")
			pg=pg+dimag(dot_product(Uk(1,:)*dconjg(Uk(2,:))-Uk(2,:)*dconjg(Uk(1,:))-&
				Uk(4,:)*dconjg(Uk(3,:))+Uk(3,:)*dconjg(Uk(4,:)),fk)*gk)
		case("sdw")
			pg=pg+0.5d0*dot_product(Uk(1,:)*dconjg(Uk(2,:))+Uk(2,:)*dconjg(Uk(1,:))+&
				Uk(4,:)*dconjg(Uk(3,:))+Uk(3,:)*dconjg(Uk(4,:)),fk)
		end select
		sc=sc+dot_product(Uk(1,:)*dconjg(Uk(3,:))-Uk(2,:)*dconjg(Uk(4,:)),fk)*gk
	end subroutine
	!subroutine matrix_inv(A)
		!implicit none
		!complex(8) :: A(:,:)
		!complex(8), allocatable :: ctmp(:)
		!real(8) :: tmp
		!integer :: i,j,k,n,m,l,mx(size(A,1),2)
		!l=size(A,1)
		!allocate(ctmp(l))
		!do i=1,l
			!tmp=0d0
			!do j=i,l
				!do k=i,l
					!if(tmp<abs(A(j,k))) then
						!tmp=abs(A(j,k))
						!mx(i,1)=j
						!mx(i,2)=k
					!endif
				!enddo
			!enddo
			!ctmp=A(i,:)
			!A(i,:)=A(mx(i,1),:)
			!A(mx(i,1),:)=ctmp
			!ctmp=A(:,i)
			!A(:,i)=A(:,mx(i,2))
			!A(:,mx(i,2))=ctmp
			!A(i,i)=1d0/A(i,i)
			!do j=1,l-1
				!n=mod(i+j-1,l)+1
				!do k=1,l-1
					!m=mod(k+i-1,l)+1
					!A(i,m)=A(i,i)*A(i,m)
					!A(n,m)=A(n,m)-A(n,i)*A(i,m)
				!enddo
				!A(n,i)=-A(n,i)*A(i,i)
			!enddo
		!enddo
		!do i=l,1,-1
			!ctmp=A(:,i)
			!A(:,i)=A(:,mx(i,1))
			!A(:,mx(i,1))=ctmp
			!ctmp=A(i,:)
			!A(i,:)=A(mx(i,2),:)
			!A(mx(i,2),:)=ctmp
		!enddo
	!end subroutine
	subroutine matrix_inv(A)
		use lapack95
		implicit none
		complex(8) :: A(:,:,:,:,:),work(45)
		integer :: ipiv(4),LDA=4,n,n1,n2,n3,i,j,k,info,lwork=45
		n1=size(A,1)
		n2=size(A,2)
		n3=size(A,3)
		n=4
		!$OMP PARALLEL DO PRIVATE(ipiv,info) SCHEDULE(GUIDED)
		do i=1,n1
			do j=1,n2
				do k=1,n3
					!write(*,*)"start"
					!call zgetrf(n,n,A(i,j,k,:,:),LDA,ipiv,info)
					call getrf(A(i,j,k,:,:),ipiv,info)
					if(info/=0) then
						write(*,*)"error1",info
						stop
					endif
					!call zgetri(n,A(i,j,k,:,:),LDA,ipiv,work,lwork,info)
					call getri(A(i,j,k,:,:),ipiv,info)
					!write(*,*)"finish"
					if(info/=0) then
						write(*,*)"error2",info
						stop
					endif
				enddo
			enddo
		enddo
		!$OMP END PARALLEL DO
	end subroutine
end module
program main
	use global
	use green
	use fft, only: fft1d, fft3d
	implicit none
	complex(8) :: iomgf(mo),iomgb(mo),romg(mro),G0(size(iomgb),mq,mq,4,4),G0_inv(size(iomgb),mq,mq,4,4),rG0(size(romg),mq,mq,4,4),&
		Veff(size(G0,1),size(G0,2),size(G0,3),2,2),sfeg(size(G0,1),size(G0,2),size(G0,3),4,4),&
		Gf(size(G0,1),size(G0,2),size(G0,3),4,4), rVeff(size(rG0,1),size(rG0,2),size(rG0,3),2,2),&
		rsfeg(size(rG0,1),size(rG0,2),size(rG0,3),4,4),rGf(size(rG0,1),size(rG0,2),size(rG0,3),4,4)
	real(8) :: k(2),sp=0.1d0,bt,sc,psc,pg,ppg,ap=0.1d0,th,gap,Tc,pk(3),pk0(3),peak(2),&
		nd(3)=(/0.0025d0,0.114d0,0.178d0/),Td(3)=(/0d0,40d0,70d0/),DOS(mro),DJq,ek
	integer :: i,j,l,ntmp
	call fileopen
	call gnuplot
	!do i=24,96,1
	!do i=24,72,8
	!nf=0.94d0
	!Tk=0.02d0
	Tk=0.005d0
	!iomgf=(/(cmplx(0d0,pi*Tk*(2d0*(i-1)+1)),i=-size(iomgf)/2+1,size(iomgf)/2)/)
	ntmp=size(iomgf)
	iomgf=(/(cmplx(0d0,pi*Tk*(2*(i-i/(ntmp/2)*ntmp)+1)),i=0,ntmp-1)/)
	!iomgb=(/(cmplx(0d0,pi*Tk*(2d0*i)),i=-size(iomgb)/2+1,size(iomgb)/2)/)
	ntmp=size(iomgb)
	iomgb=(/(cmplx(0d0,pi*Tk*(2*(i-i/(ntmp/2)*ntmp))),i=0,ntmp-1)/)
	ntmp=size(romg)
	romg=(/(cmplx(omgrg(1)+real(i)/ntmp*(omgrg(2)-omgrg(1)),it),i=1,ntmp)/)
	!romg=(/cmplx(0.01d0,it),cmplx(0.3d0,it),cmplx(0.51d0,it),cmplx(0.7d0,it)/)
	nf=0.88d0
	!nf=1.15d0
	!do i=1,3
	do
		!nf=1d0-nd(i)
		psc=0d0
		!Tk=0.05d0
		do
		!do j=1,3
			!Tk=Td(j)+0.001
			sc=sc+0.1d0
			pg=pg+0.1d0
			sc=0d0
			!pg=0d0
			!ap=0d0
			call selfconsist_tg(pg,sc,ap,sp,0)
			write(*,"('calculate greenfunction')")
			call greenfunc(pg,sc,ap,sp,iomgf,G0,G0_inv)
			write(*,"('finish')")
			write(*,"('calculate effint')")
			call effint(pg,sc,ap,sp,iomgb,Veff)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!call effint(pg,sc,ap,sp,romg,rVeff)
			!call spinsus(pg,sc,ap,sp,(/pi,pi/),romg,rVeff(:,1,1,:,:))
			!write(*,"('export data')")
			!do i=1,mq
				!write(50,"(e16.4,$)")dimag(rVeff(:,i,i,1,1))
				!write(50,"(1X)")
			!enddo
			!write(50,"(1X)")
			!stop
			!do l=1,4
				!do i=1,mq
					!do j=1,mq
						!write(50,"(e16.4,$)")dimag(rVeff(l,i,j,1,1))
					!enddo
					!write(50,"(1X)")
				!enddo
				!write(50,"(1X)")
			!enddo
			!stop
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			write(*,"('calculate self energy')")
			call selfe(G0,Veff,sfeg)
			!ntmp=size(sfeg,1)/2
			!!call pade(iomgf(ntmp+1:ntmp*2),sfeg(ntmp+1:ntmp*2,mq/2,int(0.075*mq),1),romg,rsfeg(:,mq/2,int(0.075*mq),1))
			!call pade(iomgf(1:ntmp),sfeg(1:ntmp,mq/2,int(0.075*mq),1,2),romg,rsfeg(:,mq/2,int(0.075*mq),1,2))
			!call pade(iomgf(1:ntmp),sfeg(1:ntmp,mq/2,int(0.075*mq),1,1),romg,rsfeg(:,mq/2,int(0.075*mq),1,1))
			!write(30,"(5e16.4)")(real(romg(i)),rsfeg(i,mq/2,int(0.075*mq),1,1),rsfeg(i,mq/2,int(0.075*mq),1,2),i=1,size(romg))
			!stop
			!Gf=G0_inv-sfeg
			Gf=G0_inv
			call matrix_inv(Gf)
			!write(30,"(8e16.4)")Gf(mo/2,mq/2,int(0.075*mq),:,:)
			!write(30,"(8e16.4)")Gf(mo/2,mq/2,int(0.075*mq),:,:)
			!write(30,"(8e16.4)")G0(mo/2,mq/2,int(0.075*mq),:,:)
			!stop
			write(*,"('calculate pade')")
			ntmp=size(Gf,1)/2
			!$OMP PARALLEL DO SCHEDULE(GUIDED)
			do i=1,mq
				do j=1,mq
					call pade(iomgf(1:ntmp),Gf(1:ntmp,i,j,1,1),romg,rGf(:,i,j,1,1))
				enddo
			enddo
			!$OMP END PARALLEL DO
			!call pade(iomgf(1:ntmp),Gf(1:ntmp,mq/2,int(0.075*mq),1,1),romg,rGf(:,mq/2,int(0.075*mq),1,1))
			!call pade(iomgf(ntmp+1:ntmp*2),Gf(ntmp+1:ntmp*2,mq/2,int(0.075*mq),1,1),romg,rGf(:,mq/2,int(0.075*mq),1,1))
			write(*,"('export data')")
			!write(30,"(3e16.4)")(dimag(iomgf(i)),sfeg(i,mq/2,int(0.75*mq)),i=1,size(iomgf))
			write(30,"(3e16.4)")(real(romg(i)),rGf(i,mq/2,int(0.075*mq),1,1),i=1,size(romg))
			!write(30,"(3e16.4)")(dimag(iomgf(i)),Gf(i,mq/2,int(0.075*mq),1,1),i=1,size(iomgf))
			do i=1,mq/2
				do j=1,mq/2
					write(50,"(e16.4,$)")sum(dimag(rGf(mro/2-10:mro/2+1,i,j,1,1)))
				enddo
				write(50,"(1X)")
			enddo
			write(50,"(1X)")
			stop
			!call rpainstable(sc,ap,sp,(/pi,pi/))
			write(40,"(4e12.4)")nf,Tk,sp,sc
			!if(sc<cvg*100d0.and.dd<cvg*100d0) then
			!!if(sc<cvg*100d0) then
				!exit
			!endif
			if(Tk>0d0) then
				exit
			endif
			Tk=Tk+1d0
			psc=sc
			write(30,"(1X)")
		enddo
		nf=nf-0.002d0
		if(nf<0.99d0) then
			exit
		endif
		write(10,"(1X)")
		write(50,"(1X)")
		write(40,"(1X/)")
		write(70,"(1X)")
		write(80,"(1X)")
	enddo
	call fileclose
end
subroutine fileclose()
	implicit none
	close(10)
	close(20)
	close(30)
	close(40)
	close(50)
	close(60)
	close(70)
	close(80)
end
subroutine fileopen()
	implicit none
	open(unit=10,file="../data/energy.dat")
	open(unit=20,file="../data/fermi.dat")
	open(unit=30,file="../data/EDC.dat")
	open(unit=40,file="../data/gaptemp.dat")
	open(unit=50,file="../data/rpa.dat")
	open(unit=60,file="../data/gaptheta.dat")
	open(unit=70,file="../data/ramantemp.dat")
	open(unit=80,file="../data/phase_tJ.dat")
end
subroutine gnuplot()
	implicit none
	integer :: i
	do i=10,80,10
		write(i,"(A)")'reset'
		write(i,"(A)")'#custom'
		write(i,"(A)")'#multiplot'
		write(i,"(A)")'ix=1'
		write(i,"(A)")'iy=1'
		write(i,"(A)")'sx=5*ix'
		write(i,"(A)")'sy=4*iy'
		write(i,"(A)")'xlb="T"'
		write(i,"(A)")'ylb="gap size"'
		write(i,"(A)")'ml=0.8/sx'
		write(i,"(A)")'mb=0.6/sy'
		write(i,"(A)")'mr=0.1/sx'
		write(i,"(A)")'mt=0.1/sy'
		write(i,"(A)")'gx=(1.-ml-mr)/ix'
		write(i,"(A)")'gy=(1.-mb-mt)/iy'
		write(i,"(A)")'#term'
		write(i,"(A)")'#set term eps font ",18" size sx,sy'
		write(i,"(A)")'#set output "-."."eps"'
		write(i,"(A)")'set label ylb font ",20" rotate by 90 center at character 0.8,screen 0.5*(1+mb-mt)'
		write(i,"(A)")'set label xlb font ",20" center at screen 0.5*(1+ml-mr),character 0.5'
		write(i,"(A)")'set xtics out nomirror'
		write(i,"(A)")'set ytics out nomirror'
		write(i,"(A)")'#'
		write(i,"(A)")'set size gx,gy'
		write(i,"(A)")'set tmargin 0'
		write(i,"(A)")'set bmargin 0'
		write(i,"(A)")'set lmargin 0'
		write(i,"(A)")'set rmargin 0'
		write(i,"(A)")'set multiplot'
		write(i,"(A)")'set mxtics 5'
		write(i,"(A)")'set mytics 5'
		write(i,"(A)")'unset xlabel'
		write(i,"(A)")'unset ylabel'
		write(i,"(A)")'unset key'
		write(i,"(A)")'do for[i=0:(ix*iy-1)]{'
		write(i,"(A)")'	oxi=i%ix'
		write(i,"(A)")'	oyi=i/ix+1'
		write(i,"(A)")'	set origin ml+oxi*gx,mb+(iy-oyi)*gy'
		write(i,"(A)")'	if(oxi!=0){'
		write(i,"(A)")'		set format y ""'
		write(i,"(A)")'	}'
		write(i,"(A)")'	if(oyi!=iy){'
		write(i,"(A)")'		set format x ""'
		write(i,"(A)")'	}'
		write(i,"(A)")'	if(i==(ix*iy-1)){'
		write(i,"(A)")'		set key font ",16" at screen 1-mr-0.01/sx,1-mt-0.01/sy horizontal maxcols 1 opaque autotitle columnhead'
		write(i,"(A)")'	}'
		write(i,"(A)")'	set label sprintf("(%1.3f)",0.11+0.005*i) at graph 0.1/(gx*sx),1.-0.15/(gy*sy)'
	enddo
	!plot energy
	write(10,"(A)")'	set xzeroaxis'
	write(10,"(A)")'	plot [0:][:] for[j=0:11] "-" index i every 5 using 0:2*j+1:(column(2*j+2)*0.5) with points pt 7 ps variable lc j+1 notitle'
	!plot 
	write(20,"(A)")'	plot [:55][0:0.1] for[j=0:2] "-" index i using 2:j+4 with line title word("SC DDW Total",j+1)'
	!plot phase diagram
	write(30,"(A)")'	plot [:55][0:0.1] for[j=0:2] "-" index i using 2:j+4 with line title word("SC DDW Total",j+1)'
	!plot phase diagram
	write(40,"(A)")'	plot [:55][0:0.1] for[j=0:2] "-" index i using 2:j+4 with line title word("SC DDW Total",j+1)'
	!plot phase diagram
	write(50,"(A)")'	plot [:1.2][0:] for[j=0:8] "-" index j using 1:2 with line title "".sprintf("%1.3f",0.02+0.02*j)'
	!plot phase diagram
	write(60,"(A)")'	plot [:55][0:0.1] for[j=0:2] "-" index i using 2:j+4 with line title word("SC DDW Total",j+1)'
	!plot phase diagram
	write(70,"(A)")'	plot [:1.1][0.2:1.3] for[j=0:12:2] "-" using 1:i:(100.) lc (i/2+1) smooth acsplines notitle'
	!plot phase diagram
	write(80,"(A)")'set style fill transparent solid 0.3 border'
	write(80,"(A)")'set style rect fc lt -1 fs solid 0.15 noborder'
	write(80,"(A)")'set obj rect from 0.14, graph 0 to 0.175, graph 1 behind'
	write(80,"(A)")'plot [0.06:0.24][0:100] "-" using (1-$1):3:(1e7) smooth acsplines with filledcu y1=0 lc 3 lw 2, "-" u (1-$1):4:(1e9) smooth acsplines with filledcu x1=0 lc 1 lw 2, "-" u (1-$1):5:(1e9) smooth acsplines with filledcu x1=0 lc 1 lw 2 '
	do i=10,80,10
		write(i,"(A)")'	unset label'
		write(i,"(A)")'	unset format'
		write(i,"(A)")'}'
		write(i,"(A)")'	unset multiplot'
		write(i,"(A)")'if(GPVAL_TERM eq "qt"){'
		write(i,"(A)")'	pause -1'
		write(i,"(A)")'}'
		write(i,"(A)")'#data'
	enddo
end

