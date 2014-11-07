module M_pmt
	use M_const
	implicit none
	integer, parameter :: Ns(2)=(/10,10/),Ns2=Ns(1)*Ns(2),ne=100,ne2=ne/2
	integer :: neb(Ns2,4,3)
	real(8), parameter :: t(3)=(/1d0,0d0,0d0/),DJ=0.5d0
end module
module M_vmc
	use M_pmt
	use M_matrix
	use M_rd
	use M_latt
	use M_utility
	implicit none
contains
	subroutine ini(var,wf)
		complex(8) :: wf(0:,0:)
		real(8) :: k(2),var(:),Ek,sck,tmp
		integer :: i,j,ix(2),ik(2)
		wf=0d0
		do j=1,Ns2
			call square_one2two(j,Ns,ik)
			k=(/2d0*pi/Ns(1)*ik(2),2*pi/Ns(2)*ik(1)+pi/Ns(2)/)
			Ek=-2d0*t(1)*(cos(k(1))+cos(k(2)))-var(2)
			tmp=(cos(k(1))-cos(k(2)))*var(1)
			tmp=tmp/(Ek+sqrt(tmp**2+sck**2))
			do i=1,Ns2
				call square_one2two(i,Ns,ix)
				wf(ix(1),ix(2))=wf(ix(1),ix(2))+tmp*exp(img*(k(1)*ix(2)+k(2)*ix(1)))
				if(isnan(real(wf(ix(1),ix(2))))) then
					write(*,"(A)")"NAN in wf, quit!!"
					stop
				endif
			enddo
		enddo
		wf=wf/maxval(abs(wf))
	end subroutine
	subroutine mc(wf,Nmc,E)
		complex(8) :: pb,A(ne2,ne2),iA(ne2,ne2),vu(ne2,2),wf(:,:),Y(2,2),jw,E,El
		real(8) :: rpb
		integer :: cfg(Ns2),icfg(Ns2),i,j,sg,n,Nmc,Nhot,cr(2),acp
		character(3) :: flag
		n=0
		Nhot=Nmc
		do
			if(n==0) then
				E=0d0
				acp=0
				do i=1,Ns2
					cfg(i)=i
				enddo
				call fisher_yates_shuffle(cfg,Ns2)
				do i=1,Ns2
					icfg(cfg(i))=i
				enddo
				do i=1,ne2
					call diffsc(vu,(/-10,i/),wf,cfg,1)
					A(i,:)=vu(:,1)
				enddo
				iA=A
				call matrix_inv(iA)
			endif
			n=n+1
			if(n==Nhot) then
				call energy(cfg,icfg,wf,A,iA,El)
			endif
			call irandom(i,ne)
			call irandom(j,Ns2-ne2)
			vu=0d0
			cr=0
			if(j>ne2.and.i<=ne2) then
				sg=1
				cr(1)=i
				cr(2)=j+ne2
				call diffsc(vu,cr,wf,cfg,sg)
			elseif(j>ne2.and.i>ne2) then
				sg=2
				cr(2)=i
				cr(1)=j+ne2
				call diffsc(vu,cr,wf,cfg,sg)
			else
				sg=3
				if(i>ne2) then
					call swap(i,j)
				else
					j=j+ne2
				endif
				cr=(/i,j/)
				call diffsc(vu,cr,wf,cfg,sg)
			endif
			jw=1d0
			call det(vu,cr,A,iA,Y,pb,sg)
			!call jast(cfg,cr,sg,v,jw)
			!call checkcfg(icfg,80)
			!write(*,*)A
			!write(*,*)real(pb*jw*dconjg(pb*jw))
			if(abs(pb)>1d10) then
				write(*,*)"initial configure may be not good"
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
	subroutine diffsc(vu,r,wf,cfg,sg)
		complex(8) :: vu(:,:),wf(0:,0:)
		integer :: sg,i,ir(2),ii(2),r(2),cfg(:),n
		real(8) :: s
		s=1d0
		vu=0d0
		if(sg==1.or.sg==3) then
			call square_one2two(cfg(r(2)),Ns,ir)
			do i=ne2+1,ne
				if(r(2)==i) then
					call square_one2two(cfg(r(1)),Ns,ii)
				else
					call square_one2two(cfg(i),Ns,ii)
				endif
				ii=(ir-ii)
				s=sign(1d0,ii(1)+0.1d0)
				ii(1)=mod(ii(1)+Ns(2),Ns(2))
				ii(2)=mod(ii(2)+Ns(1),Ns(1))
				if(ii(1)==0.and.ii(2)==0) then
					write(*,*)"some thing go wrong"
					stop
				endif
				vu(i-ne2,1)=wf(ii(1),ii(2))*s
			enddo
		endif
		if(sg==2.or.sg==3) then
			call square_one2two(cfg(r(1)),Ns,ir)
			do i=1,ne2
				if(r(1)==i) then
					call square_one2two(cfg(r(2)),Ns,ii)
				else
					call square_one2two(cfg(i),Ns,ii)
				endif
				ii=-(ir-ii)
				s=sign(1d0,ii(1)+0.1d0)
				ii(1)=mod(ii(1)+Ns(2),Ns(2))
				ii(2)=mod(ii(2)+Ns(1),Ns(1))
				if(ii(1)==0.and.ii(2)==0) then
					write(*,*)"some thing go wrong"
					stop
				endif
				vu(i,2)=wf(ii(1),ii(2))*s
			enddo
		endif
	end subroutine
	subroutine det(vu,cr,A,iA,Y,pb,sg)
		complex(8) :: vu(:,:),A(:,:),iA(:,:),Y(:,:),pb
		integer :: cr(:),sg,i,j,n,ii,jj
		if(sg==3) then
			call det_ratio_rowcol(vu,(/cr(1),cr(2)-ne2/),A,iA,Y,pb)
		else
			if(sg==1) then
				call det_ratio_row(vu(:,1),cr(1),A,iA,pb)
			else
				call det_ratio_col(vu(:,2),cr(2)-ne2,A,iA,pb)
			endif
		endif
	end subroutine
	subroutine update(vu,cr,A,iA,Y,pb,sg)
		complex(8) :: vu(:,:),A(:,:),iA(:,:),Y(:,:),pb
		integer :: cr(:),n,s,i,j,sg,ii,jj
		if(sg==3) then
			call inv_update_rowcol(vu,(/cr(1),cr(2)-ne2/),Y,A,iA)
		else
			if(sg==1) then
				call inv_update_row(vu(:,1),cr(1),pb,A,iA)
			else
				call inv_update_col(vu(:,2),cr(2)-ne2,pb,A,iA)
			endif
		endif
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
				k=neb(cfg(i),ii,1)
				j=icfg(k)
				if(j>ne) then
					!hopping
					sg=(i-1)/ne2+1
					cr(sg)=i
					cr(mod(sg,2)+1)=j
					call diffsc(vu,cr,wf,cfg,sg)
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
					cr=(/i,j/)
					call diffsc(vu,cr,wf,cfg,3)
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
end module
program main
	use M_vmc
	complex(8) :: wf(Ns(2),Ns(1)),E
	real(8) :: var(2)=(/0.01d0,0d0/),dvar(2),Eb,E2,l
	integer :: n,i,j
	open(10,file="../data/test.dat")
	write(10,"(A)")"#data"
	call init_random_seed()
	call square(Ns,neb)
	n=20
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
		l=l+0.1d0
		if(l>2d0) then
			exit
		endif
	enddo
end program
