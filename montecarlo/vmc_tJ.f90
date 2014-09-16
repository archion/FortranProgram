module global
	implicit none
	integer, parameter :: Ns(2)=(/10,10/),Ns2=Ns(1)*Ns(2),ne=Ns2,ne2=ne/2
	integer :: latt(Ns2,4,3),pbcx(-Ns(1)+1:2*Ns(1)),pbcy(-Ns(2)+1:2*Ns(2))
	real(8), parameter :: t(3)=(/1d0,0d0,0d0/),DJ=0.5d0,pi=atan(1d0)*4d0,cvg=1d-5
	complex(8), parameter :: img=(0d0,1d0)
	contains
		subroutine one2two(n,i,j)
			integer :: n,i,j
			i=(n-1)/Ns(1)
			j=mod(n-1,Ns(1))
		end subroutine
		subroutine two2one(i,j,n)
			integer :: n,i,j
			n=i*Ns(1)+j+1
		end subroutine
		subroutine gen_latt_square()
!					j=0...n-1
!			 1 --  2 --  3 --  4
!			 |     |     |     |           4      3   4
!			 5 --  6 --  7 --  8           |       \ /
!i=0...n-1	 |     |     |     |        3--0--1     0
!			 9 -- 10 -- 11 -- 12           |       / \
!			 |     |     |     |           2      2   1
!			13 -- 14 -- 15 -- 16
			implicit none
			integer :: ii,i,j
			pbcx=(/(i,i=1,Ns(1)),(i,i=1,Ns(1)),(i,i=1,Ns(1))/)
			pbcy=(/(i,i=1,Ns(2)),(i,i=1,Ns(2)),(i,i=1,Ns(2))/)
			do ii=1,Ns2
				call one2two(ii,i,j)
				call two2one(i,mod(j+1,Ns(1)),latt(ii,1,1))
				call two2one(mod(i+1,Ns(2)),j,latt(ii,2,1))
				call two2one(i,mod(j-1+Ns(1),Ns(1)),latt(ii,3,1))
				call two2one(mod(i-1+Ns(2),Ns(2)),j,latt(ii,4,1))
				!latt(ii,1,1)=Ns(1)*(i-1)+pbcx(j+1)
				!latt(ii,2,1)=Ns(1)*(pbcy(i+1)-1)+j
				!latt(ii,3,1)=Ns(1)*(i-1)+pbcx(j-1)
				!latt(ii,4,1)=Ns(1)*(pbcy(i-1)-1)+j
				!latt(ii,1,2)=Ns(1)*(pbcy(i+1)-1)+pbcx(j+1)
				!latt(ii,2,2)=Ns(1)*(pbcy(i+1)-1)+pbcx(j-1)
				!latt(ii,3,2)=Ns(1)*(pbcy(i-1)-1)+pbcx(j-1)
				!latt(ii,4,2)=Ns(1)*(pbcy(i-1)-1)+pbcx(j+1)
				!latt(ii,1,3)=Ns(1)*(i-1)+pbcx(j+2)
				!latt(ii,2,3)=Ns(1)*(pbcy(i+2)-1)+j
				!latt(ii,3,3)=Ns(1)*(i-1)+pbcx(j-2)
				!latt(ii,4,3)=Ns(1)*(pbcy(i-2)-1)+j
			enddo
			!write(*,*)latt(13,:,1)
			!stop
		end subroutine gen_latt_square
		subroutine cmwrite(f,A)
			complex(8) :: A(:,:)
			integer :: f,i
			do i=1,size(A,1)
				!write(f,"(sp,e16.4,e11.4'I'$)")A(i,:)
				!write(f,"(sp,e10.2','e9.2$)")A(i,:)
				write(f,"(e16.4$)")real(A(i,:))
				write(f,"(1X)")
			enddo
			write(f,"(1X)")
		end subroutine
		subroutine rmwrite(f,A)
			real(8) :: A(:,:)
			integer :: f,i
			do i=1,size(A,1)
				write(f,"(e16.4$)")A(i,:)
				write(f,"(1X)")
			enddo
			write(f,"(1X)")
		end subroutine
		subroutine imwrite(f,A)
			integer :: f,i,A(:,:)
			do i=1,size(A,1)
				write(f,"(i4$)")A(i,:)
				write(f,"(1X)")
			enddo
			write(f,"(1X)")
		end subroutine
end module
module rd
	implicit none
	contains
		subroutine irandom(i,n)
			real(8) :: rn
			integer :: i,n
			call random_number(rn)
			i=1+int(rn*(n))
		end subroutine
		subroutine fisher_yates_shuffle(a,n)
			integer :: a(:),i,j,tmp,n
			real(8) :: rn
			do i=n,1,-1
				call random_number(rn)
				j=ceiling(rn*i)
				tmp=a(i)
				a(i)=a(j)
				a(j)=tmp
			enddo
		end subroutine
		subroutine init_random_seed()
			integer :: i, n, clock
			integer, dimension(:), allocatable :: seed
			call random_seed(size = n)
			allocate(seed(n))
			call system_clock(count=clock)
			seed = clock + 37 * (/ (i - 1, i = 1, n) /)
			call random_seed(put = seed)
			deallocate(seed)
		end subroutine
end module
module vmc
	use global
	use rd
	implicit none
	contains
		subroutine ini(var,wf)
			complex(8) :: wf(0:,0:)
			real(8) :: k(2),var(:),Ek,sck,tmp
			integer :: i,j,ix,jx,ik,jk
			wf=0d0
			do j=1,Ns2
				call one2two(j,ik,jk)
				k=(/2d0*pi/Ns(1)*jk,2*pi/Ns(2)*ik+pi/Ns(2)/)
				!k=(/2d0*pi/Ns(1)*jk,2d0*pi/Ns(2)*ik/)
				Ek=-2d0*t(1)*(cos(k(1))+cos(k(2)))-var(2)
				!sck=(cos(k(1))-cos(k(2)))*var(1)
				sck=(cos(k(1))+cos(k(2)))*var(1)
				!sck=var(1)
				tmp=sck/(Ek+sqrt(Ek**2+sck**2))
				!if(Ek<0d0) then
					!tmp=1d0
				!else
					!tmp=0d0
				!endif
				!tmp=sck
				!wf(ik,jk)=tmp
				do i=1,Ns2
					!exit
					call one2two(i,ix,jx)
					wf(ix,jx)=wf(ix,jx)+tmp*exp(img*(k(1)*jx+k(2)*ix))
					!wf(ix,jx)=wf(ix,jx)+tmp*cos(k(1)*jx+k(2)*ix)
					if(isnan(real(wf(ix,jx)))) then
						write(*,"(A)")"NAN in wf, quit!!"
						stop
					endif
				enddo
			enddo
			!write(*,*)wf
			!stop
			wf=wf/maxval(abs(wf))
			!call rmwrite(10,imag(wf))
			!call rmwrite(10,real(wf))
			!stop
		end subroutine
		subroutine mc(wf,Nmc,E)
			complex(8) :: pb,A(ne2,ne2),iA(ne2,ne2),vu(ne2,2),wf(:,:),Y(2,2),jw,E,El
			real(8) :: rpb
			integer :: cfg(Ns2),icfg(Ns2),i,j,sg,n,Nmc,Nhot,cr(2),acp
			E=0d0
			do i=1,Ns2
				cfg(i)=i
			enddo
			call fisher_yates_shuffle(cfg,Ns2)
			Nhot=Nmc
			!write(*,*)"***********hot start**********"
			!n=0
			!do 
				!n=n+1
				!call random_number(i,ne)
				!call random_number(j,Ns2-ne2)
				!if(j>ne2) then
					!j=ne2+j
					!call swap(cfg(i),cfg(j))
				!else
					!if(i>ne2) then
						!call swap(i,j)
					!else
						!j=j+ne2
					!endif
					!call swap(cfg(i),cfg(j))
					!cycle
				!endif
				!if(n>Nhot) then
					!exit
				!endif
			!enddo
			!write(*,*)"**********hot end**************"
			do i=1,Ns2
				icfg(cfg(i))=i
			enddo
			do i=1,ne2
				call pair(vu,(/-10,i/),wf,cfg,0)
				A(i,:)=vu(:,1)
			enddo
			iA=A
			!call checkcfg(icfg,80)
			!write(*,*)real(A)
			!stop
			call matrix_inv(iA)
			if(n>Nhot) then
				call energy(cfg,icfg,wf,A,iA,El)
			endif
			!iA=matmul(iA,A)
			!call cmwrite(10,iA)
			!stop
			n=0
			acp=0
			do
				n=n+1
				call irandom(i,ne)
				call irandom(j,Ns2-ne2)
				vu=0d0
				cr=0
				if(j>ne2) then
					j=ne2+j
					sg=(i-1)/ne2+1
					cr(sg)=i
					cr(mod(sg,2)+1)=j
					call pair(vu,cr,wf,cfg,sg-1)
				else
					sg=3
					if(i>ne2) then
						call swap(i,j)
					else
						j=j+ne2
					endif
					cr=(/i,j/)
					call pair(vu,cr,wf,cfg,sg-1)
					!write(*,*)vu(cr(2),1)-vu(cr(1),2)
				endif
				jw=1d0
				call det(vu,cr,A,iA,Y,pb,sg)
				!call jast(cfg,cr,sg,v,jw)
				!call checkcfg(icfg,80)
				!write(*,*)A
				!write(*,*)real(pb*jw*dconjg(pb*jw))
				if(abs(pb)>1d10) then
					write(*,*)"initial configure may be not good"
					call fisher_yates_shuffle(cfg,Ns2)
					do i=1,Ns2
						icfg(cfg(i))=i
					enddo
					do i=1,ne2
						call pair(vu,(/-10,i/),wf,cfg,0)
						A(i,:)=vu(:,1)
					enddo
					iA=A
					call matrix_inv(iA)
					n=0
					cycle
				endif
				call random_number(rpb)
				if(rpb<abs(pb*jw)**2) then
				!if(1<abs(pb*jw)**2) then
					acp=acp+1
					call swap(icfg(cfg(i)),icfg(cfg(j)))
					call swap(cfg(i),cfg(j))
					call update(vu,cr,A,iA,Y,pb,sg)
					call checkinv(A,iA,rpb)
					!if(rpb>1d-8) then
						!write(*,*)n,rpb
						!iA=A
						!call matrix_inv(iA)
						!call checkinv(A,iA,rpb)
						!write(*,*)n,rpb
					!endif
					!call checkcfg(icfg,80)
					!write(*,*)A
					!write(*,*)pb
					!exit
					if(mod(n,500)==0) then
						iA=A
						call matrix_inv(iA)
					endif
					if(n>Nhot) then
						!call checkcfg(icfg,80)
						call energy(cfg,icfg,wf,A,iA,El)
					endif
					!write(*,*)El
					!exit
				endif
				if(n>Nhot) then
					E=E+El
				endif
				if(n>=(Nmc+Nhot)) then
					E=E/(Nmc*Ns2)
					exit
				endif
			enddo
			write(*,"('accept/total number is ',I6,'/',I6)")acp,n
		end subroutine
		subroutine pair(vu,r,wf,cfg,sg)
			complex(8) :: vu(:,:),wf(0:,0:)
			integer :: sg,i,ir,jr,ii,jj,r(2),cfg(:),n
			real(8) :: s
			s=1d0
			vu=0d0
			if(sg==0.or.sg==2) then
				call one2two(cfg(r(2)),ir,jr)
				do i=ne2+1,ne
					if(r(2)==i) then
						call one2two(cfg(r(1)),ii,jj)
					else
						call one2two(cfg(i),ii,jj)
					endif
					ii=(ir-ii)
					jj=(jr-jj)
					s=sign(1d0,ii+0.1d0)
					ii=mod(ii+Ns(2),Ns(2))
					jj=mod(jj+Ns(1),Ns(1))
					if(ii==0.and.jj==0) then
						write(*,*)"some thing go wrong"
					endif
					vu(i-ne2,1)=wf(ii,jj)*s
				enddo
			endif
			if(sg==1.or.sg==2) then
				call one2two(cfg(r(1)),ir,jr)
				do i=1,ne2
					if(r(1)==i) then
						call one2two(cfg(r(2)),ii,jj)
					else
						call one2two(cfg(i),ii,jj)
					endif
					ii=-(ir-ii)
					jj=-(jr-jj)
					s=sign(1d0,ii+0.1d0)
					ii=mod(ii+Ns(2),Ns(2))
					jj=mod(jj+Ns(1),Ns(1))
					if(ii==0.and.jj==0) then
						write(*,*)"some thing go wrong"
					endif
					vu(i,2)=wf(ii,jj)*s
				enddo
			endif
		end subroutine
		subroutine det(vu,cr,A,iA,Y,pb,sg)
			complex(8) :: vu(:,:),A(:,:),iA(:,:),Y(:,:),pb
			integer :: cr(:),sg,i,j,n,ii,jj
			n=size(A,1)
			Y=0d0
			if(sg==3) then
				ii=cr(1)
				jj=cr(2)-ne2
				vu(:,1)=vu(:,1)-A(ii,:)
				vu(jj,1)=0d0
				vu(:,2)=vu(:,2)-A(:,jj)
				Y(2,1)=-iA(jj,ii)
				do i=1,n
					Y(1,1)=Y(1,1)+iA(jj,i)*vu(i,2)
					Y(2,2)=Y(2,2)+vu(i,1)*iA(i,ii)
					do j=1,n
						Y(1,2)=Y(1,2)-vu(i,1)*iA(i,j)*vu(j,2)
					enddo
				enddo
				Y(1,1)=Y(1,1)+1d0
				Y(2,2)=Y(2,2)+1d0
			else
				if(sg==1) then
					ii=cr(1)
					vu(:,1)=vu(:,1)-A(ii,:)
					do i=1,n
						Y(1,1)=Y(1,1)+vu(i,1)*iA(i,ii)
					enddo
					Y(1,1)=Y(1,1)+1d0
					Y(2,2)=1d0
				else
					jj=cr(2)-ne2
					vu(:,2)=vu(:,2)-A(:,jj)
					Y(1,1)=1d0
					do i=1,n
						Y(2,2)=Y(2,2)+iA(jj,i)*vu(i,2)
					enddo
					Y(2,2)=Y(2,2)+1d0
				endif
			endif
			pb=Y(1,1)*Y(2,2)-Y(1,2)*Y(2,1)
		end subroutine
		subroutine update(vu,cr,A,iA,Y,pb,sg)
			complex(8) :: vu(:,:),A(:,:),iA(:,:),tmp(size(iA,1),size(iA,2)),tmp1(size(iA,1),2),Y(:,:),pb
			integer :: cr(:),n,s,i,j,sg,ii,jj
			tmp1=0d0
			n=size(iA,1)
			if(sg==3) then
				ii=cr(1)
				jj=cr(2)-ne2
				A(ii,:)=A(ii,:)+vu(:,1)
				A(:,jj)=A(:,jj)+vu(:,2)
				do i=1,n
					do j=1,n
						tmp1(i,1)=tmp1(i,1)+iA(j,i)*vu(j,1)
						tmp1(i,2)=tmp1(i,2)+iA(i,j)*vu(j,2)
					enddo
				enddo
				do i=1,n
					do j=1,n
						tmp(i,j)=iA(i,j)-1d0/pb*&
							((iA(i,ii)*Y(1,1)+tmp1(i,2)*Y(2,1))*tmp1(j,1)+(iA(i,ii)*Y(1,2)+tmp1(i,2)*Y(2,2))*iA(jj,j))
					enddo
				enddo
			else
				if(sg==1) then
					ii=cr(1)
					A(ii,:)=A(ii,:)+vu(:,1)
					do i=1,n
						do j=1,n
							tmp1(i,1)=tmp1(i,1)+iA(j,i)*vu(j,1)
						enddo
					enddo
					do i=1,n
						do j=1,n
							tmp(i,j)=iA(i,j)-1d0/pb*iA(i,ii)*tmp1(j,1)
						enddo
					enddo
				else
					jj=cr(2)-ne2
					A(:,jj)=A(:,jj)+vu(:,2)
					do i=1,n
						do j=1,n
							tmp1(i,2)=tmp1(i,2)+iA(i,j)*vu(j,2)
						enddo
					enddo
					do i=1,n
						do j=1,n
							tmp(i,j)=iA(i,j)-1d0/pb*tmp1(i,2)*iA(jj,j)
						enddo
					enddo
				endif
			endif
			iA=tmp
		end subroutine
		!subroutine jast(cfg,cr,sg,v,jw)
			!if(sg==3) then
				!do i=1,ne
					!s=(-1)*(mod(i,ne2))
					!jw=jw+2d0*(v(abs(cfg(i)-cfg(cr(1))))-v(abs(cfg(i)-cfg(cr(2)))))*s
				!enddo
				!jw=jw-2d0*v(abs(cfg(cr(1))-cfg(cr(2))))
			!else
				!do i=1,ne
					!s=(-1)*(mod(i,ne2))
					!jw=jw+(v(abs(cfg(i)-cfg(cr(1))))-v(abs(cfg(i)-cfg(cr(2)))))*s
				!enddo
			!endif
		!end subroutine
		subroutine energy(cfg,icfg,wf,A,iA,El)
			complex(8) :: A(:,:),iA(:,:),vu(ne2,2),Y(2,2),pb,wf(:,:),jw,El
			!real(8) :: 
			integer :: cfg(:),icfg(:),i,ii,j,k,sg,cr(2)
			El=0d0
			do i=1,ne
				do ii=1,4
					k=latt(cfg(i),ii,1)
					j=icfg(k)
					if(j>ne) then
						!hopping
						sg=(i-1)/ne2+1
						cr(sg)=mod(i-1,ne2)+1
						call pair(vu,cr,wf,cfg,sg-1)
						call det(vu,cr,A,iA,Y,pb,sg)
						jw=1d0
						!call jast(cfg,cr,sg,v,jw)
						if(abs(k-cfg(i))>2*Ns(1)) then
							El=El+t(1)*dconjg(pb)*jw
						else
							El=El-t(1)*dconjg(pb)*jw
						endif
						write(*,*)"hopping"
						cycle
					endif
					!diagnal
					!El=El-0.125d0*DJ*sign(1d0,j-ne2-0.5d0)
					if(i>ne2) then
						cycle
					endif
					if(j>ne2) then
						!diagnal
						El=El-0.5d0*DJ
						!spin flip
						cr=(/i,j/)
						call pair(vu,cr,wf,cfg,2)
						call det(vu,cr,A,iA,Y,pb,3)
						jw=1d0
						!call jast(cfg,cr,sg,v,jw)
						El=El-0.5d0*DJ*dconjg(pb)*jw
					endif
				enddo
			enddo
		end subroutine
		subroutine checkcfg(icfg,s)
			integer :: i,s,icfg(:)
			do i=1,Ns2
				if(icfg(i)<=ne2) then
					write(*,"(' ●'$)")
				elseif(icfg(i)<=ne) then
					write(*,"(' ○'$)")
				else
					write(*,"('  '$)")
				endif
				if(mod(i,Ns(1))==0) then
					write(*,"(1X)")
				endif
			enddo
			call sleepqq(s)
			write(*,"(1X)")
		end subroutine
		subroutine checkinv(A,B,p)
			complex(8) :: A(:,:),B(:,:),tmp(size(A,1),size(A,2))
			real(8) :: p
			integer :: i
			tmp=matmul(A,B)
			do i=1,size(A,1)
				tmp(i,i)=tmp(i,i)-1d0
			enddo
			tmp=abs(tmp)
			p=sum(tmp)
		end subroutine
		subroutine swap(a,b)
			integer a,b,tmp
			tmp=a
			a=b
			b=tmp
		end subroutine
		subroutine matrix_inv(A)
			use f95_precision, only : wp => dp
			use lapack95
			implicit none
			complex(wp) :: A(:,:)
			integer :: ipiv(size(A,1)),info
			call getrf(A,ipiv,info)
			if(info/=0) then
				write(*,*)"error1",info
				stop
			endif
			call getri(A,ipiv,info)
			if(info/=0) then
				write(*,*)"error2",info
				stop
			endif
		end subroutine
end module
program main
	use vmc
	complex(8) :: wf(Ns(2),Ns(1)),E
	real(8) :: var(2)=(/0.01d0,0d0/),dvar(2),Eb,E2,l
	integer :: n,i,j
	open(10,file="../data/test.dat")
	write(10,"(A)")"#data"
	call init_random_seed()
	call gen_latt_square()
	n=40
	l=-2d0
	do
		Eb=0d0
		E2=0d0
		var(1)=10d0**l
		call ini(var,wf)
		!$OMP PARALLEL DO REDUCTION(+:Eb,E2) PRIVATE(E) SCHEDULE(STATIC)
		do j=1,n
			call mc(wf,300000,E)
			Eb=Eb+real(E)
			E2=E2+real(E)**2
			!if(all(abs(dvar)<cvg)) then
			!exit
			!else
			!var=var+dvar
			!endif
			!exit
		enddo
		!$OMP END PARALLEL DO
		write(*,*)var(1),Eb/n,sqrt(abs((Eb/n)**2-E2/n))
		write(10,*)var(1),Eb/n,sqrt(abs((Eb/n)**2-E2/n))
		l=l+0.1d0
		if(l>2d0) then
			exit
		endif
	enddo
end
