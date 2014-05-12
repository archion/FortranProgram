module global
	implicit none
	save
	real(8), parameter ::t(5)=(/2d0,-0.9d0,0.000d0,0.d0,0d0/),escal=0.25d0,pi=3.1415926d0,cvg=1e-5,DJ=0.75d0,V=-3d0
	real(8) :: nf=0.8d0,Tk
	integer, parameter :: mk=128,mq=128,mo=1024,mro=512
	complex(8),parameter :: img=(0d0,1d0)
	contains
	subroutine spinsus(dt,ap,sp,q,omg,Xq)
		implicit none
		integer :: i,j,l,m,n
		complex(8) :: Uk(2,2),Ukq(2,2),omg(:),Xq(:)
		real(8) :: k(2),q(2),sp,bt,dt,ap,ek(2),ekq(2),fk(2),fkq(2),DJq
		bt=escal/Tk*1.16e4
		Xq=0d0
		!$OMP PARALLEL DO REDUCTION(+:Xq) PRIVATE(k,ek,Uk,fk,ekq,Ukq,fkq) SCHEDULE(GUIDED)
		do i=1,mq
			k(1)=2d0*pi/mq*i
			do j=1,mq
				k(2)=2d0*pi/mq*j
				call EU(k,dt,ap,sp,ek,Uk)
				call EU(k+q,dt,ap,sp,ekq,Ukq)
				fk=1d0/(1d0+exp(bt*ek))
				fkq=1d0/(1d0+exp(bt*ekq))
				do n=1,2
					do m=1,2
						Xq=Xq+(Uk(1,n)*Ukq(2,m)*dconjg(Uk(1,n)*Ukq(2,m))-Uk(2,n)*Ukq(1,m)*dconjg(Uk(1,n)*Ukq(2,m)))*&
							(1-fk(n)-fkq(m))/(omg+ek(n)+ekq(m)+1d-10)
					enddo
				enddo
			enddo
		enddo
		!$OMP END PARALLEL DO
		Xq=Xq/mq**2
	end subroutine
	subroutine greenfunc(dt,ap,sp,omg,sfeg,G0)
		implicit none
		integer :: i,j,l,m,n,o,ii,jj,sig
		complex(8) :: Uk(2,2),omg(:),G0(:,0:,0:,:),sfeg(:,0:,0:,:)
		real(8) :: dt,ap,sp,ek(2),k(2)
		G0=0d0
		l=size(G0,4)
		o=size(G0,2)
		!$OMP PARALLEL DO PRIVATE(ii,jj,sig,k,ek,Uk) SCHEDULE(GUIDED)
		do i=0,o/2
			do j=0,i
				k=(/i,j/)*2d0*pi/o
				call EU(k,dt,ap,sp,ek,Uk)
				do m=1,l
					do n=1,l
						G0(:,i,j,m)=G0(:,i,j,m)+Uk(1,n)*dconjg(Uk(m,n))/(omg-ek(n)-sfeg(:,i,j,m))
					enddo
				enddo
			enddo
		enddo
		!$OMP END PARALLEL DO
		call symfill(G0(:,:,:,1),1)
		call symfill(G0(:,:,:,2),-1)
	end subroutine
	subroutine effint(dt,ap,sp,omg,Veff)
		use fft, only : fft3d
		implicit none
		integer :: i,j,n,ii,jj
		real(8) :: dt,ap,sp,DJq,q(2)
		complex(8) :: Veff(:,0:,0:),omg(:)
		n=size(Veff,2)
			!Gtmp(size(G0,1),size(G0,2),size(G0,3),size(G0,4),size(G0,5))
		!Gtmp=G0
		!call fft3d(Gtmp(:,:,:,1,1),-1)
		!call fft3d(Gtmp(:,:,:,2,2),-1)
		!call fft3d(Gtmp(:,:,:,1,2),-1)
		!call fft3d(Gtmp(:,:,:,2,1),-1)
		!Veff=(Gtmp(:,:,:,1,1)*Gtmp(:,:,:,2,2)-Gtmp(:,:,:,1,2)*Gtmp(:,:,:,2,1))
		!call fft3d(Veff,1)
		!!$OMP PARALLEL DO PRIVATE(q,Veff,DJq) SCHEDULE(GUIDED)
		do i=0,n/2
			do j=0,i
				!q=((/i,j/)*2d0*pi/size(Veff,2)-(/pi,pi/))*0.5+(/pi,pi/)
				q=(/i,j/)*2d0*pi/n
				call spinsus(dt,ap,sp,q,omg,Veff(:,i,j))
				DJq=0.34d0*(cos(q(1))+cos(q(2)))
				Veff(:,i,j)=Veff(:,i,j)/(1d0+DJq*Veff(:,i,j))
				Veff(:,i,j)=DJq**2*Veff(:,i,j)
			enddo
		enddo
		!!$OMP END PARALLEL DO
		call symfill(Veff(:,:,:),1)
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
	subroutine effintfft(omg,G0,Veff)
		use fft, only : fft3d
		implicit none
		integer :: i,j,n,ii,jj
		real(8) :: dt,ap,sp,DJq,q(2)
		complex(8) :: Veff(:,0:,0:),omg(:),G0(:,:,:,:),&
			Gtmp(size(G0,1),size(G0,2),size(G0,3),size(G0,4))
		n=size(Veff,2)
		Gtmp=G0
		call fft3d(Gtmp(:,:,:,1),-1)
		call fft3d(Gtmp(:,:,:,2),-1)
		call fermimats(Gtmp(:,:,:,1),-1)
		call fermimats(Gtmp(:,:,:,2),-1)
		Veff=Gtmp(:,:,:,1)*Gtmp(:,:,:,1)-Gtmp(:,:,:,2)*Gtmp(:,:,:,2)
		call fft3d(Veff,1)
	end subroutine
	subroutine selfe(G0,Veff,sfeg)
		use fft
		implicit none
		complex(8) :: sfeg(:,:,:,:),Veff(:,:,:),G0(:,:,:,:),&
			Gtmp(size(G0,1),size(G0,2),size(G0,3),size(G0,4))
		real(8) :: k(2),q(2),sp,bt,dt,ap,ek(2),ekq(2),fk(2),fkq(2),DJq
		Gtmp=G0
		call fft3d(Gtmp(:,:,:,1),-1)
		call fft3d(Gtmp(:,:,:,2),-1)
		call fft3d(Veff,-1)
		sfeg(:,:,:,1)=Gtmp(:,:,:,1)*Veff
		sfeg(:,:,:,2)=Gtmp(:,:,:,2)*Veff
		call fft3d(sfeg(:,:,:,1),1)
		call fft3d(sfeg(:,:,:,2),1)
		sfeg=sfeg*Tk*size(Gtmp,1)
	end subroutine
	subroutine fermimats(G,sig)
		implicit none
		integer :: i,sig,m
		complex(8) :: w,wp,G(:,:,:)
		m=size(G,1)
		wp=cmplx(-2d0*sin(sig*pi/(2*m))**2,sin(sig*pi/m))
		w=cmplx(1d0,0d0)
		do i=1,m
			G(i,:,:)=G(i,:,:)*w
			w=w*wp+w
		enddo
	end subroutine
	subroutine rpainstable(dt,ap,sp,q)
		implicit none
		integer :: i,j,l,m,n
		complex(8) :: Uk(2,2),Ukq(2,2)
		real(8) :: k(2),q(2),sp,bt,dt,ap,ek(2),ekq(2),fk(2),fkq(2),omg,omgr(2),domg,DJq,Xq,Xq_rpa
		bt=escal/Tk*1.16e4
		Xq=0d0
		!$OMP PARALLEL DO REDUCTION(+:Xq) PRIVATE(k,ek,Uk,fk,ekq,Ukq,fkq) SCHEDULE(GUIDED)
		do i=0,mq
			k(1)=pi/mq*i
			do j=0,mq
				k(2)=pi/mq*j
				call EU(k,dt,ap,sp,ek,Uk)
				call EU(k+q,dt,ap,sp,ekq,Ukq)
				fk=1d0/(1d0+exp(bt*ek))
				fkq=1d0/(1d0+exp(bt*ekq))
				do n=1,2
					do m=1,2
						Xq=Xq+(Uk(1,n)*Ukq(2,m)*dconjg(Uk(1,n)*Ukq(2,m))-Uk(2,n)*Ukq(1,m)*dconjg(Uk(1,n)*Ukq(2,m)))*&
							(1-fk(n)-fkq(m))/(ek(n)+ekq(m))
					enddo
				enddo
			enddo
		enddo
		!$OMP END PARALLEL DO
		Xq=Xq/mq**2
		write(*,"(2e16.3)")1d0-nf,1d0+DJq*Xq
	end subroutine
	subroutine band(dt,ap,sp,ki,kf)
		implicit none
		integer :: i,l
		complex(8) :: Uk(2,2),Ukq(2,2)
		real(8) :: ki(2),kf(2),kn(2),ek(2),ekq(2),sp,dt,ap
		kn=ki
		i=0
		do
			kn=kn+(kf-ki)/sqrt(sum((kf-ki)**2))*pi/mq
			call EU(kn,dt,ap,sp,ek,Uk)
			call EU(kn+(/pi,pi/),dt,ap,sp,ekq,Ukq)
			!export band
			!write(10,"(10e16.3)")(ek(l),abs(Uk(1,l)),l=1,4)
			write(10,"(12e16.3)")(ek(l),abs(Uk(1,l)),l=1,2),(ekq(l),abs(Ukq(1,l)),l=1,2),(ek(l)+ekq(l),1.0,l=1,2)
			if(all((kn-kf)*(kf-ki)>=0d0)) then
				exit
			endif
			i=i+1
		enddo
		write(*,*)i
	end subroutine
	subroutine selfconsist_tg(dt,ap,sp)
		implicit none
		complex(8) :: Uk(2,2)
		real(8) :: k(2),n1,ap,app,ek(2),sp,sa,sb,sp0,dk,dt,dtp,wide,cvg1
		integer :: i,j,c,info
		logical :: flaga,flagb
		! dtp(2)=-2d0
		wide=0.5d0
		sp0=sp+wide
		c=0
		cvg1=0.01
		do 
			sa=sp
			sb=sp
			flaga=.true.
			flagb=.true.
			do 			
				sp=0.5d0*(sa+sb)
				dtp=0d0
				app=0d0
				n1=0d0
				c=c+1
				!call sfcsap(sp,ap)
				!$OMP PARALLEL DO REDUCTION(+:n1,dtp,app) PRIVATE(k,ek,Uk) SCHEDULE(GUIDED)
				do i=0,mk
					do j=0,mk
						k=(/pi/mk*i,pi/mk*j/)
						call EU(k,dt,ap,sp,ek,Uk)
						call order(k,ek,Uk,n1,dtp,app)
					enddo
				enddo
				!$OMP END PARALLEL DO
				n1=n1/(mk**2)
				dtp=V*dtp/(mk**2)
				app=app/(mk**2)*DJ
				!write(*,"(4(a6,e12.4))")"sa=",sa,",sp=",sp,",sb=",sb,"n=",n1
				!write(*,*)abs(n1-nf),c
				if(abs(n1-nf)<=cvg1) then
					exit
				endif
				if(n1<nf) then
					flaga=.false.
					sa=sp
					if(flagb) then
						sb=sp+wide
					endif
				else
					flagb=.false.
					sb=sp
					if(flaga) then
						sa=sp-wide
					endif
				endif
			enddo
			wide=max(abs(sp0-sp),100*cvg)
			sp0=sp
			cvg1=max(cvg,min(cvg1,0.1d0*(abs(dtp-dt))))
			!write(*,"(5(a6,e12.4)a6i4)")"n=",n1,",Tk=",Tk,",sp=",sp,",DSC=",dt,",DDW=",dd,"num:",c
			if((abs(dtp-dt)+abs(ap-app))<cvg) then
				exit
			endif
			dt=dtp
			ap=app
		enddo
		!write(*,*)"!!!!!!selfconsist return!!!!!!!!"
		write(*,"(5e12.4,i4)")n1,Tk,sp,dt,ap,c
		!write(*,*)"!!!!!!!!!!!!end!!!!!!!!!!!!!!!!!"
	end subroutine
	subroutine EU(k,dt,ap,sp,ek,Uk)
		implicit none
		complex(8) :: Uk(2,2),cfy,sfy
		real(8) :: eka,eks,ek(2),k(2),sp,dk,dt,ap,gk,th,fy(2),cos2fy,sin2fy
		eks=-4d0*(1d0-nf)*t(2)*cos(k(1))*cos(k(2))-sp-2d0*(1d0-nf)*t(3)*(cos(2d0*k(1))+cos(2d0*k(2)))
		eka=-2d0*((1d0-nf)*t(1)+ap)*(cos(k(1))+cos(k(2)))
		gk=0.5d0*(cos(k(1))-cos(k(2)))
		dk=dt*gk
		ek=sqrt((eka+eks)**2+dk**2)*(/1d0,-1d0/)
		cos2fy=(eka+eks)/sqrt((eka+eks)**2+dk**2)
		sin2fy=dk/sqrt((eka+eks)**2+dk**2)
		cfy=dcmplx(sqrt(0.5d0*(1d0+cos2fy)))
		sfy=dcmplx(sqrt(0.5d0*(1d0-cos2fy))*sign(1d0,sin2fy))
		Uk=reshape((/   cfy , sfy , & 
					   -sfy , cfy /) , (/2 , 2/))
	end subroutine
	subroutine order(k,ek,Uk,n,dt,ap)
		implicit none
		complex(8) :: Uk(2,2)
		real(8) :: k(2),ek(2),dt,gk,n,fk(2),bt,ap
		bt=escal/Tk*1.16e4
		fk=1d0/(1d0+exp(bt*ek))
		gk=0.5d0*(cos(k(1))-cos(k(2)))
		n=n+1d0+dot_product(Uk(1,:)*dconjg(Uk(1,:))-Uk(2,:)*dconjg(Uk(2,:)),fk)
		ap=ap+dot_product(Uk(1,:)*dconjg(Uk(1,:)),fk)*cos(k(1))
		dt=dt+dot_product(Uk(1,:)*dconjg(Uk(2,:)),fk)*gk
	end subroutine
end module
program main
	use global
	use green
	use fft, only: fft1d, fft3d
	implicit none
	complex(8) :: iomgf(mo),iomgb(mo),romg(mro),G0(size(iomgb),mq,mq,2),rG0(size(romg),mq,mq,2),&
		Veff(size(G0,1),size(G0,2),size(G0,3)),sfeg(size(G0,1),size(G0,2),size(G0,3),2),Gf(size(G0,1),size(G0,2),size(G0,3),2),&
		rVeff(size(rG0,1),size(rG0,2),size(rG0,3)),rsfeg(size(rG0,1),size(rG0,2),size(rG0,3),2),&
		rGf(size(rG0,1),size(rG0,2),size(rG0,3),2)
	real(8) :: k(2),sp=0d0,bt,dt,pdt,ap=0.1d0,th,gap,Tc,pk(3),pk0(3),peak(2),&
		nd(3)=(/0.0025d0,0.114d0,0.178d0/),Td(3)=(/0d0,40d0,70d0/),it=0.02d0,omgrg(2)=(/-2d0,2d0/),DOS(mro),DJq,ek
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
	!romg=(/cmplx(0.1d0,it),cmplx(0.3d0,it),cmplx(0.51d0,it)/)
	nf=0.88d0
	!do i=1,3
	do
		!nf=1d0-nd(i)
		pdt=0d0
		!Tk=0.05d0
		do
		!do j=1,3
			!Tk=Td(j)+0.001
			dt=dt+0.1
			call selfconsist_tg(dt,ap,sp)
			write(*,"('calculate greenfunction')")
			sfeg=0d0
			call greenfunc(dt,ap,sp,iomgf,sfeg,G0)
			write(*,"('finish')")
			!ntmp=size(G0,1)/2
			!do i=1,mq
				!do j=1,mq
					!call pade(iomgf(1:ntmp),G0(1:ntmp,i,j,2),romg,rG0(:,i,j,2))
					!write(50,"(e16.4,$)")sum(dimag(rG0(1:100,i,j,2)))
				!enddo
				!write(50,"(1X)")
			!enddo
			!write(50,"(1X)")
			!stop
			write(*,"('calculate effint')")
			!call effint(dt,ap,sp,romg,Veff)
			call effint(dt,ap,sp,iomgb,Veff)
			!call spinsus(dt,ap,sp,(/pi,pi/),romg,Xq)
			!DJq=-0.68d0
			!write(50,"(5e16.3)")(dimag(iomgb(i)),dimag(Xq(i)),dimag(Xq(i)/(1d0+DJq*Xq(i))Xq_rpa(i)),i=1,size(iomgb))
			!write(50,"(3e16.3)")(dimag(iomgb(i)),dreal(Xq(i)),dreal(Veff(i,mr/2,mr/2)),i=1,size(iomgb))
			!write(50,"(1X)")
			!call pade(iomgb,Veff(:,mr/2,mr/2),romg,Veff(:,mr/2,mr/2))
			!call pade(iomgb,Xq,romg,Xq)
			!call spinsus(dt,ap,sp,(/pi,pi/),romg,Xq)
			!write(30,"(3e16.3)")(dreal(romg(i)),dimag(Xq(i)),dimag(Veff(i,mr/2,mr/2)),i=1,size(romg))
			!write(30,"(1X)")
			!write(*,"('pade')")
			!ntmp=size(Veff,1)/2
			!!$OMP PARALLEL DO SCHEDULE(GUIDED)
			!do i=1,mq/2+1
				!do j=1,i
					!!call pade(iomgb(1:ntmp),Veff(1:ntmp,mq/2,mq/2),romg,rVeff(:,mq/2,mq/2))
					!!write(30,"(2e16.3)")Veff(:,1,1)
					!call pade(iomgb(1:ntmp),Veff(1:ntmp,i,j),romg,rVeff(:,i,j))
					!!write(30,"(1X/)")
					!!write(30,"(2e16.3)")rVeff(:,1,1)
				!enddo
			!enddo
			!!$OMP END PARALLEL DO
			!call symfill(rVeff,1)
			!write(*,"('export data')")
			!do i=1,mq
				!do j=1,mq
					!write(50,"(e16.4,$)")dimag(rVeff(0.1*mro,i,j))
				!enddo
				!write(50,"(1X)")
			!enddo
			!write(50,"(1X)")
			!do i=1,mq
				!do j=1,mq
					!write(50,"(e16.4,$)")dimag(rVeff(0.3*mro,i,j))
				!enddo
				!write(50,"(1X)")
			!enddo
			!write(50,"(1X)")
			!do i=1,mq
				!do j=1,mq
					!write(50,"(e16.4,$)")dimag(rVeff(0.51*mro,i,j))
				!enddo
				!write(50,"(1X)")
			!enddo
			!stop
			!write(50,"(1X/)")
			!call spinsus(dt,ap,sp,(/pi,pi/),romg,Xq)
			!write(30,"(2e16.3)")(real(romg(i)),dimag(Xq(i)),i=1,size(romg))
			!write(30,"(1X)")
			!do i=1,mk
				!do j=1,mk
					!call pade(iomgb,Veff(:,i,j),romg,rVeff(:,i,j))
				!enddo
			!enddo
			!!$OMP END PARALLEL DO
			!do i=1,ms/2
				!tmp=Veff(ms-i+1,:,:)
				!Veff(ms-i+1,:,:)=Veff(i,:,:)
				!Veff(i,:,:)=tmp
			!enddo
			write(*,"('calculate self energy')")
			call selfe(G0,Veff,sfeg)
			call greenfunc(0d0,ap,sp,iomgf,sfeg,Gf)
			do i=1,mq
				do j=1,mq
					k=(/i,j/)*2d0*pi/mq
					sfeg(:,i,j,2)=dt*0.5d0*(cos(k(1))-cos(k(2)))+sfeg(:,i,j,2)
				enddo
			enddo
			do i=1,mo
				Gf(i,:,:,1)=(1d0-nf)*Gf(i,:,:,1)/(1d0+sfeg(i,:,:,2)*dconjg(sfeg(i,:,:,2))*Gf(mo-i+1,:,:,1)*Gf(i,:,:,1))
			enddo
			write(*,"('calculate pade')")
			!ntmp=size(sfeg,1)/2
			!call pade(iomgf(ntmp+1:ntmp*2),sfeg(ntmp+1:ntmp*2,mq/2,int(0.75*mq)),romg,rsfeg(:,mq/2,int(0.75*mq)))
			ntmp=size(Gf,1)/2
			!call pade(iomgf(ntmp+1:ntmp*2),Gf(ntmp+1:ntmp*2,mq/2,int(0.75*mq),1),romg,Gf(:,mq/2,int(0.75*mq),1))
			call pade(iomgf(1:ntmp),Gf(1:ntmp,mq/2,int(0.75*mq),1),romg,rGf(:,mq/2,int(0.75*mq),1))
			write(*,"('export data')")
			!write(30,"(3e16.4)")(real(romg(i)),rsfeg(i,mq/2,int(0.75*mq)),i=1,size(romg))
			write(30,"(3e16.4)")(real(romg(i)),rGf(i,mq/2,int(0.75*mq),1),i=1,size(romg))
			!write(30,"(3e16.4)")(dimag(iomgf(i)),sfeg(i,mq/2,int(0.75*mq)),i=1,size(iomgf))
			stop
			!do i=1,mk
				!do j=1,mk
					!write(50,"(e16.4,$)")dimag(rVeff(38,i,j))
				!enddo
				!write(50,"(1X)")
			!enddo
			!write(50,"(1X/)")
			!do i=1,mk
				!do j=1,mk
					!write(50,"(e16.4,$)")dimag(rVeff(12,i,j))
				!enddo
				!write(50,"(1X)")
			!enddo
			!call rpainstable(dt,ap,sp,(/pi,pi/))
			write(40,"(4e12.4)")nf,Tk,sp,dt
			!call band(dt,ap,sp,(/0d0,0d0/),(/pi,pi/))
			!call band(dt,ap,sp,(/pi,pi/),(/pi,0d0/))
			!call band(dt,ap,sp,(/pi,0d0/),(/0d0,0d0/))
			!if(dt<cvg*100d0.and.dd<cvg*100d0) then
			!!if(dt<cvg*100d0) then
				!exit
			!endif
			if(Tk>0d0) then
				exit
			endif
			Tk=Tk+1d0
			pdt=dt
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

