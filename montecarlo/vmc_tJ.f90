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
			real(8) :: A(:,:)
			integer :: f,i
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
				k=(/2d0*pi/Ns(1)*jk,2d0*pi/Ns(2)*ik-pi/Ns(2)/)-pi
				!k=(/2d0*pi/Ns(1)*jk,2d0*pi/Ns(2)*ik/)-pi
				Ek=-2d0*t(1)*(cos(k(1))+cos(k(2)))-var(2)
				sck=(cos(k(1))-cos(k(2)))*var(1)
				!sck=var(1)
				tmp=sck/(Ek+sqrt(Ek**2+sck**2))
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
				call pair(vu(:,1),cfg(i),wf,cfg,0)
				A(i,:)=vu(:,1)
			enddo
			iA=A
			call matrix_inv(iA)
			call energy(cfg,icfg,wf,A,iA,El)
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
					cr(sg)=mod(i-1,ne2)+1
					call pair(vu(:,sg),cfg(j),wf,cfg,sg-1)
				else
					sg=3
					if(i>ne2) then
						call swap(i,j)
					else
						j=j+ne2
					endif
					cr=(/i,j-ne2/)
					call pair(vu(:,1),cfg(j),wf,cfg,0)
					call pair(vu(:,2),cfg(i),wf,cfg,1)
					!write(*,*)vu(cr(2),1)-vu(cr(1),2)
				endif
				!A=0d0
				!do i=1,ne2
					!do j=1,ne2
						!call random_number(rpb)
						!A(i,j)=cmplx(0d0,rpb)
						!call random_number(rpb)
						!A(i,j)=A(i,j)+rpb
						!!A(i,j)=rpb
					!enddo
					!call random_number(rpb)
					!!vu(i,2)=rpb
					!vu(i,2)=cmplx(0d0,rpb)
					!call random_number(rpb)
					!vu(i,2)=vu(i,2)+rpb
					!call random_number(rpb)
					!!vu(i,1)=rpb
					!vu(i,1)=cmplx(0d0,rpb)
					!call random_number(rpb)
					!vu(i,1)=vu(i,1)+rpb
				!enddo
				!iA=A
				!call matrix_inv(iA)
				!sg=3
				!cr(1)=1
				!cr(2)=1
				!call det(vu,cr,A,iA,Y,pb,sg)
				!call cmwrite(10,A)
				!call cmwrite(10,iA)
				!call ifinv(A,iA,rpb)
				!write(*,*)rpb
				!call update(vu,cr,A,iA,Y,pb,sg)
				!call cmwrite(10,A)
				!call cmwrite(10,iA)
				!call cmwrite(10,vu)
				!write(*,*)rpb
				!stop
				jw=1d0
				call det(vu,cr,A,iA,Y,pb,sg)
				!call jast(cfg,cr,sg,v,jw)
				call random_number(rpb)
				!write(*,*)real(pb*jw*dconjg(pb*jw))
				if(rpb<abs(pb*jw)**2) then
					acp=acp+1
					call swap(icfg(cfg(i)),icfg(cfg(j)))
					call swap(cfg(i),cfg(j))
					call update(vu,cr,A,iA,Y,pb,sg)
					!call checkinv(A,iA,rpb)
					!if(rpb>1d-6) then
						!write(*,*)n,rpb
						!iA=A
						!call matrix_inv(iA)
						!call checkinv(A,iA,rpb)
						!write(*,*)n,rpb
					!endif
					!call checkcfg(icfg,80)
					if(mod(Nmc,500)==0) then
						iA=A
						call matrix_inv(iA)
					endif
					call energy(cfg,icfg,wf,A,iA,El)
					!write(*,*)El
					!exit
				endif
				E=E+El
				if(n>=Nmc) then
					E=E/(Nmc*Ns2)
					exit
				endif
			enddo
			!write(*,"('accept/total number is ',I5,'/',I5)")acp,n
		end subroutine
		subroutine pair(vu,r,wf,cfg,sg)
			complex(8) :: vu(:),wf(0:,0:)
			integer :: sg,s1,i,ir,jr,ii,jj,r,cfg(:),n
			real(8) :: s2
			s1=(-1)**sg
			s2=1d0
			n=abs(sg-1)*size(vu,1)
			call one2two(r,ir,jr)
			do i=1,size(vu,1)
				call one2two(cfg(n+i),ii,jj)
				ii=s1*(ir-ii)
				jj=s1*(jr-jj)
				s2=sign(1d0,jj+0.1d0)
				ii=mod(ii+Ns(2),Ns(2))
				jj=mod(jj+Ns(1),Ns(1))
				vu(i)=wf(ii,jj)*s2
			enddo
		end subroutine
		subroutine det(vu,cr,A,iA,Y,pb,sg)
			complex(8) :: vu(:,:),A(:,:),iA(:,:),Y(:,:),pb
			integer :: cr(:),sg,i,j,n
			n=size(A,1)
			Y=0d0
			if(sg==3) then
				vu(:,1)=vu(:,1)-A(cr(1),:)
				vu(cr(2),1)=0d0
				vu(:,2)=vu(:,2)-A(:,cr(2))
				Y(2,1)=-iA(cr(2),cr(1))
				do i=1,n
					Y(1,1)=Y(1,1)+iA(cr(2),i)*vu(i,2)
					Y(2,2)=Y(2,2)+vu(i,1)*iA(i,cr(1))
					do j=1,n
						Y(1,2)=Y(1,2)-vu(i,1)*iA(i,j)*vu(j,2)
					enddo
				enddo
				Y(1,1)=Y(1,1)+1d0
				Y(2,2)=Y(2,2)+1d0
			else
				if(sg==1) then
					vu(:,1)=vu(:,1)-A(cr(1),:)
					do i=1,n
						Y(1,1)=Y(1,1)+vu(i,1)*iA(i,cr(1))
					enddo
					Y(1,1)=Y(1,1)+1d0
					Y(2,2)=1d0
				else
					vu(:,2)=vu(:,2)-A(:,cr(2))
					Y(1,1)=1d0
					do i=1,n
						Y(2,2)=Y(2,2)+iA(cr(2),i)*vu(i,2)
					enddo
					Y(2,2)=Y(2,2)+1d0
				endif
			endif
			pb=Y(1,1)*Y(2,2)-Y(1,2)*Y(2,1)
		end subroutine
		subroutine update(vu,cr,A,iA,Y,pb,sg)
			complex(8) :: vu(:,:),A(:,:),iA(:,:),tmp(size(iA,1),size(iA,2)),tmp1(size(iA,1),2),Y(:,:),pb
			integer :: cr(:),n,s,i,j,sg
			tmp1=0d0
			n=size(iA,1)
			if(sg==3) then
				A(cr(1),:)=A(cr(1),:)+vu(:,1)
				A(:,cr(2))=A(:,cr(2))+vu(:,2)
				do i=1,n
					do j=1,n
						tmp1(i,1)=tmp1(i,1)+iA(j,i)*vu(j,1)
						tmp1(i,2)=tmp1(i,2)+iA(i,j)*vu(j,2)
					enddo
				enddo
				do i=1,n
					do j=1,n
						tmp(i,j)=iA(i,j)-1d0/pb*&
							((iA(i,cr(1))*Y(1,1)+tmp1(i,2)*Y(2,1))*tmp1(j,1)+(iA(i,cr(1))*Y(1,2)+tmp1(i,2)*Y(2,2))*iA(cr(2),j))
					enddo
				enddo
			else
				if(sg==1) then
					A(cr(1),:)=A(cr(1),:)+vu(:,1)
					do i=1,n
						do j=1,n
							tmp1(i,1)=tmp1(i,1)+iA(j,i)*vu(j,1)
						enddo
					enddo
					do i=1,n
						do j=1,n
							tmp(i,j)=iA(i,j)-1d0/pb*iA(i,cr(1))*tmp1(j,1)
						enddo
					enddo
				else
					A(:,cr(2))=A(:,cr(2))+vu(:,2)
					do i=1,n
						do j=1,n
							tmp1(i,2)=tmp1(i,2)+iA(i,j)*vu(j,2)
						enddo
					enddo
					do i=1,n
						do j=1,n
							tmp(i,j)=iA(i,j)-1d0/pb*tmp1(i,2)*iA(cr(2),j)
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
						call pair(vu(:,sg),cfg(j),wf,cfg,sg-1)
						call det(vu,cr,A,iA,Y,pb,sg)
						jw=1d0
						!call jast(cfg,cr,sg,v,jw)
						if(abs(k-cfg(i))>2*Ns(1)) then
							El=El+t(1)*dconjg(pb)*jw
						else
							El=El-t(1)*dconjg(pb)*jw
						endif
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
						cr=(/i,j-ne2/)
						call pair(vu(:,1),cfg(j),wf,cfg,0)
						call pair(vu(:,2),cfg(i),wf,cfg,1)
						call det(vu,cr,A,iA,Y,pb,3)
						jw=1d0
						!call jast(cfg,cr,sg,v,jw)
						El=El+0.5d0*DJ*dconjg(pb)*jw
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
	n=2
	l=-2d0
	do
		Eb=0d0
		E2=0d0
		var(1)=10d0**l
		call ini(var,wf)
		!$OMP PARALLEL DO REDUCTION(+:Eb,E2) PRIVATE(E) SCHEDULE(STATIC)
		do j=1,n
			call mc(wf,50000,E)
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
		exit
		l=l+0.2d0
	enddo
end
