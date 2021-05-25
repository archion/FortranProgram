include "mc_utility.f90"
module model_utility
	use M_hamilton_final_m
	use mc_utility
	use M_omp_rand
	implicit none
	procedure(change_interface), pointer :: change => null()
	procedure(get_row_interface), pointer :: get_row => null()
	interface
		subroutine change_interface(cfg,icfg,dcfg,k,rnd,is_P)
			import :: randomNumberSequence
			integer :: cfg(:),icfg(:)
			integer :: dcfg(0:),k(0:)
			type(randomNumberSequence) :: rnd
			logical :: is_P
		end subroutine
		logical function get_row_interface(ca,cfg,k,sg,P)
			integer, intent(in) :: cfg(:),ca(:)
			integer, intent(out) :: k(0:),sg
			integer, intent(in) :: P(:)
		end function
	end interface
	integer :: Ns=0
	type(t_ham) :: Hmf,Ham,Hja
	logical :: is_project=.false.,is_ph=.true.
contains
	subroutine change_ph(cfg,icfg,dcfg,k,rnd,is_P)
		!cfg: from site index to particle index
		!icfg: from particle index to site index
		!dcfg: index pair: particle index, new site index
		integer :: cfg(:),icfg(:)
		integer :: dcfg(0:),k(0:)
		type(randomNumberSequence) :: rnd
		logical :: is_P
		integer :: i1,j1,i2,j2,n
		!choose a random electron
		!call random_number(i,size(icfg))
		i1=mt_random_number(rnd,size(icfg))
		!get the index in H
		i1=icfg(i1)
		!get the spin index
		n=Hmf%H2s(i1,2)
		!get the index in H for same site but opposite spin
		i2=Hmf%s2H(Hmf%H2s(i1,1),mod(n,2)+1)
		do 
			!choose another random site
			!call random_number(j,Ns)
			j1=mt_random_number(rnd,Ns)

			!!choose another random nn site
			!j1=mt_random_number(rnd,size(latt%nb(1)%st(Hmf%H2s(i1,1))%j))
			!j1=latt%nb(1)%st(Hmf%H2s(i1,1))%j(j1)


			!get the index in H for the opposite spin for this site
			j2=Hmf%s2H(j1,mod(n,2)+1)
			!get the index in H for the same spin for this site
			j1=Hmf%s2H(j1,n)
			if(cfg(j1)==0) then
				if((.not.is_P).or.sign(1,-cfg(j2))==sign(1,-cfg(i2))) then
				!hop
					!set the index for creation and annilation 
					dcfg(0:2)=[2,j1,-i1]
					!set the index for changing row in A to row in H 
					k(0:2)=[2,cfg(abs(dcfg(2))),dcfg(1)]
				else
				!exchange opposite spin case
					dcfg(0:4)=[4,j2,-i2,j1,-i1]
					k(0:4)=[4,cfg(abs(dcfg(4))),dcfg(3),cfg(abs(dcfg(2))),dcfg(1)]
				endif
				exit
			!else
				!k(0)=0
			endif
		enddo
	end subroutine
	logical function get_row_ph(ca,cfg,k,sg,P) result(get_row)
		integer, intent(in) :: cfg(:),ca(:)
		integer, intent(out) :: k(0:),sg
		integer, intent(in) :: P(:)
		integer :: i,j,n1,n2,l,tp1,tp2,n,tmp
		logical :: flag(size(ca)),fP(2)
		get_row=.false.
		flag=.true.
		n1=-1; n2=0
		l=0
		sg=1
		n=size(P)
		do i=size(ca),1,-1
			tp1=merge(1,2,Hmf%H2s(abs(ca(i)),2)==1)
			fP(tp1)=(cfg(abs(ca(i)))/=0)
			if(flag(i)) then
				do j=i-1,1,-1
					if(flag(j)) then
						if(ca(i)==-ca(j)) then
							if((ca(i)>0).xor.fP(tp1)) then
								flag(j)=.false.
								sg=sg*sign(1,-mod(count(flag(j+1:i-1)),2))
								exit
							else
								return
							endif
						elseif(ca(i)==ca(j)) then
							return 
						endif
					endif
				enddo
				if(j==0) then
					if(ca(i)<0.and.fp(tp1)) then
						n1=n1+2
						k(n1)=cfg(-ca(i))
						if(n2>n1.and.mod((n2-n1)/2,2)==0) sg=-sg
					elseif(ca(i)>0.and.(.not.fp(tp1))) then
						n2=n2+2
						k(n2)=ca(i)
						if(n1>n2.and.mod((n1-n2)/2,2)==0) sg=-sg
					else
						return
					endif
				endif
			endif
			if(n>0) then
				fP(tp1)=.not.(fP(tp1).xor.flag(i))
				tp2=mod(tp1,2)+1
				tmp=Hmf%s2H(Hmf%H2s(abs(ca(i)),1),tp2)
				fP(tp2)=cfg(tmp)/=0
				do j=i+1,size(ca)
					if(abs(ca(j))==tmp) then
						fp(tp2)=.not.fp(tp2)
					endif
				enddo
				if(fP(1).and.(.not.fP(2))) then
					l=l-1
				else
					fP(tp1)=(ca(i)>0)
					if(fP(1).and.(.not.fP(2))) then
						l=l+1
					endif
				endif
				if(P(n)==i) then
					n=n-1
					if(l/=0) then
						return
					endif
				endif
			endif
		enddo
		get_row=.true.
		k(0)=n2
	end function
	subroutine change_noph(cfg,icfg,dcfg,k,rnd,is_P)
		integer :: cfg(:),icfg(:)
		integer :: dcfg(0:),k(0:)
		type(randomNumberSequence) :: rnd
		integer :: i1,j1,i2,j2,n
		logical :: is_P
		!call random_number(i,size(icfg))
		i1=mt_random_number(rnd,size(icfg))
		i1=icfg(i1)
		n=Hmf%H2s(i1,2)
		i2=Hmf%s2H(Hmf%H2s(i1,1),mod(n,2)+1)
		do 
			!call random_number(j,Ns)
			j1=mt_random_number(rnd,Ns)
			j2=Hmf%s2H(j1,mod(n,2)+1)
			j1=Hmf%s2H(j1,n)
			if(cfg(j1)==0) then
				if((.not.is_P).or.sign(1,-cfg(j2))==sign(1,-cfg(i2))) then
					dcfg(0:2)=[2,j1,-i1]
					k(0:2)=[2,cfg(abs(dcfg(2))),dcfg(1)]
				else
					dcfg(0:4)=[4,i2,-j2,j1,-i1]
					k(0:4)=[4,cfg(abs(dcfg(4))),dcfg(3),cfg(abs(dcfg(2))),dcfg(1)]
				endif
				exit
			endif
		enddo
	end subroutine
	logical function get_row_noph(ca,cfg,k,sg,P) result(get_row)
		integer, intent(in) :: cfg(:),ca(:)
		integer, intent(out) :: k(0:),sg
		integer, intent(in) :: P(:)
		integer :: i,j,n1,n2,l,tp1,tp2,n,tmp
		logical :: flag(size(ca)),fP(2)
		get_row=.false.
		flag=.true.
		n1=-1; n2=0
		l=0
		sg=1
		n=size(P)
		do i=size(ca),1,-1
			tp1=merge(1,2,Hmf%H2s(abs(ca(i)),2)==1)
			fP(tp1)=(cfg(abs(ca(i)))/=0)
			if(flag(i)) then
				do j=i-1,1,-1
					if(flag(j)) then
						if(ca(i)==-ca(j)) then
							if((ca(i)>0).xor.fP(tp1)) then
								flag(j)=.false.
								sg=sg*sign(1,-mod(count(flag(j+1:i-1)),2))
								exit
							else
								return
							endif
						elseif(ca(i)==ca(j)) then
							return 
						endif
					endif
				enddo
				if(j==0) then
					if(ca(i)<0.and.fp(tp1)) then
						n1=n1+2
						k(n1)=cfg(-ca(i))
						if(n2>n1.and.mod((n2-n1)/2,2)==0) sg=-sg
					elseif(ca(i)>0.and.(.not.fp(tp1))) then
						n2=n2+2
						k(n2)=ca(i)
						if(n1>n2.and.mod((n1-n2)/2,2)==0) sg=-sg
					else
						return
					endif
				endif
			endif
			if(n>0) then
				fP(tp1)=.not.(fP(tp1).xor.flag(i))
				tp2=mod(tp1,2)+1
				tmp=Hmf%s2H(Hmf%H2s(abs(ca(i)),1),tp2)
				fP(tp2)=cfg(tmp)/=0
				do j=i+1,size(ca)
					if(abs(ca(j))==tmp) then
						fp(tp2)=.not.fp(tp2)
					endif
				enddo
				if(fP(1).and.fP(2)) then
					l=l-1
				else
					fP(tp1)=(ca(i)>0)
					if(fP(1).and.fP(2)) then
						l=l+1
					endif
				endif
				if(P(n)==i) then
					n=n-1
					if(l/=0) then
						return
					endif
				endif
			endif
		enddo
		get_row=.true.
		k(0)=n2
	end function
	subroutine get_phy2(var,cfg,WA)
		type(t_var) :: var
		complex(wp) :: WA(:,:)
		integer :: cfg(:)
		complex(wp) :: pb,bd1,bd2,expq
		integer :: n1,n2,sg,ic,jc1,jc2,nb,iextdat
		integer :: k(0:4),dcfg_all(0:4)
		integer :: Pr(0:1)
		character :: flag1,flag2
		real(wp) :: q(3),r1
		Pr(0:1)=[1,1]
		if(.not.is_project) then
			Pr(0)=0
		endif
		var%val=0._wp
		nb=var%nb
		dcfg_all(0)=size(var%c,1)*2
		iextdat=1
		q=0._wp
		r1=0._wp
		if(allocated(var%extdat)) then
			if(var%extdat(iextdat)==real(ichar("q"),8)) then
				q=var%extdat(2:4)
				iextdat=iextdat+4
			endif
			if(var%extdat(iextdat)==real(ichar("r"),8)) then
				r1=var%extdat(iextdat+1)
				iextdat=iextdat+2
			endif
		endif
		do n2=1,size(latt%nb(nb)%bd)
			if(allocated(var%extdat)) then
				if(var%extdat(iextdat)==real(ichar("f"),8)) then
					if(nint(var%extdat(iextdat+n2))==0) then
						cycle
					endif
				endif
			endif
			do n1=1,size(latt%nb(nb)%bd)
				if(sqrt(sum((latt%nb(nb)%bd(n1)%r-latt%nb(nb)%bd(n2)%r)**2))<r1) cycle
				expq=exp(img*sum(q*(latt%nb(nb)%bd(n1)%r-latt%nb(nb)%bd(n2)%r)))
				flag1=" "
				do jc1=1,size(var%c,2)
					if(all(merge(latt%nb(nb)%bd(n1)%sb(1),latt%nb(nb)%bd(n1)%sb(2),var%c(:,jc1)%l=="i")-var%c(:,jc1)%sb==0)) then
						if(flag1==" ") flag1="i"
					elseif(all(merge(latt%nb(nb)%bd(n1)%sb(2),latt%nb(nb)%bd(n1)%sb(1),var%c(:,jc1)%l=="i")-var%c(:,jc1)%sb==0)) then
						if(flag1==" ") flag1="j"
					else
						cycle
					endif
					dcfg_all(1:size(var%c,1))=[(Hmf%s2H(merge(latt%nb(nb)%bd(n1)%i(1),latt%nb(nb)%bd(n1)%i(2),var%c(ic,jc1)%l==flag1),abs(var%c(ic,jc1)%tp))*sign(1,-var%c(ic,jc1)%tp),ic=size(var%c,1),1,-1)]
					if(var%c(1,jc1)%l=="i") then
						bd1=conjg(latt%nb(nb)%bd(n1)%bdc)
					else
						bd1=latt%nb(nb)%bd(n1)%bdc
					endif
					bd1=bd1*merge(var%bd(n1),conjg(var%bd(n1)),var%cg(jc1))*var%sg(jc1)
					flag2=" "
					do jc2=1,size(var%c,2)
						if(all(merge(latt%nb(nb)%bd(n2)%sb(1),latt%nb(nb)%bd(n2)%sb(2),var%c(:2,jc2)%l=="i")-var%c(:,jc2)%sb==0)) then
							if(flag2==" ") flag2="i"
						elseif(all(merge(latt%nb(nb)%bd(n2)%sb(2),latt%nb(nb)%bd(n2)%sb(1),var%c(:2,jc2)%l=="i")-var%c(:,jc2)%sb==0)) then
							if(flag2==" ") flag2="j"
						else
							cycle
						endif
						dcfg_all(size(var%c,1)+1:dcfg_all(0))=[(Hmf%s2H(merge(latt%nb(nb)%bd(n2)%i(1),latt%nb(nb)%bd(n2)%i(2),var%c(ic,jc2)%l==flag2),abs(var%c(ic,jc2)%tp))*sign(1,var%c(ic,jc2)%tp),ic=1,size(var%c,1))]
						if(var%c(1,jc2)%l=="i") then
							bd2=latt%nb(nb)%bd(n2)%bdc
						else
							bd2=conjg(latt%nb(nb)%bd(n2)%bdc)
						endif
						bd2=bd2*merge(conjg(var%bd(n2)),var%bd(n2),var%cg(jc2))*var%sg(jc2)
						if(get_row(dcfg_all(1:dcfg_all(0)),cfg,k,sg,Pr(1:Pr(0)))) then
							call get_pb(k(1:k(0)),shape(0),pb,WA=WA)
							var%val(var%bd2v(n1))=var%val(var%bd2v(n1))+real(sg*pb*bd1*bd2*exp(jast(cfg,dcfg_all(1:dcfg_all(0))))*expq)
						endif
					enddo
				enddo
			enddo
		enddo
		if(.not.isnan(var%V)) then
			var%val=var%val*var%V
		endif
	end subroutine
	subroutine get_phy1(var,cfg,WA)
		type(t_var) :: var
		complex(wp) :: WA(:,:)
		integer :: cfg(:)
		complex(wp) :: pb,bdc
		integer :: i,j,n,l,p,sg,ic,jc,sb(2)
		integer :: k(0:4),dcfg_all(0:4)
		integer :: Pr(0:1)
		character :: flag
		Pr(0:1)=[1,1]
		if(.not.is_project) then
			Pr(0)=0
		endif
		var%val=0._wp
		l=var%nb
		dcfg_all(0)=size(var%c,1)
		do n=1,size(latt%nb(l)%bd)
			i=latt%nb(l)%bd(n)%i(1)
			j=latt%nb(l)%bd(n)%i(2)
			sb=latt%nb(l)%bd(n)%sb
			flag=" "
			do jc=1,size(var%c,2)
				if(all(merge(sb(1),sb(2),var%c(:,jc)%l=="i")-var%c(:,jc)%sb==0)) then
					if(flag==" ") flag="i"
				elseif(all(merge(sb(2),sb(1),var%c(:,jc)%l=="i")-var%c(:,jc)%sb==0)) then
					if(flag==" ") flag="j"
				else
					cycle
				endif
				dcfg_all(1:dcfg_all(0))=[(Hmf%s2H(merge(i,j,var%c(ic,jc)%l==flag),abs(var%c(ic,jc)%tp))*sign(1,var%c(ic,jc)%tp),ic=1,size(var%c,1))]
				if(get_row(dcfg_all(1:dcfg_all(0)),cfg,k,sg,Pr(1:Pr(0)))) then
					call get_pb(k(1:k(0)),shape(0),pb,WA=WA)
					if(size(var%c,1)==2.and.var%c(1,jc)%l/=var%c(2,jc)%l) then
						if(var%c(1,jc)%l=="i") then
							bdc=latt%nb(l)%bd(n)%bdc
						else
							bdc=conjg(latt%nb(l)%bd(n)%bdc)
						endif
					else
						bdc=abs(latt%nb(l)%bd(n)%bdc)
					endif
					var%val(var%bd2v(n))=var%val(var%bd2v(n))+real(sg*var%bd(n)*var%sg(jc)*pb*bdc*exp(jast(cfg,dcfg_all(1:dcfg_all(0)))))
				endif
			enddo
		enddo
		if(.not.isnan(var%V)) then
			var%val=var%val*var%V
		endif
	end subroutine
	complex(wp) function get_energy(cfg,WA,dcfg,m,iA,AW,WAW)
		complex(wp) :: WA(:,:)
		integer :: cfg(:)
		integer, optional :: dcfg(:),m(:)
		complex(wp), optional :: AW(:,:),WAW(:,:),iA(:,:)
		complex(wp) :: pb,bdc,val,tmp
		integer :: i,j,n,l,p,sg,ic,jc,sb(2)
		integer :: k(0:4),dcfg_(0:2),dcfg_all(0:4)
		integer :: Pr(0:2)
		character :: flag
		if(present(dcfg)) then
			dcfg_(0:)=[size(dcfg),dcfg]
			Pr=[2,1,dcfg_(0)+1]
		else
			dcfg_(0)=0
			Pr(0:1)=[1,1]
		endif
		if(.not.is_project) then
			Pr(0)=0
		endif
		get_energy=0._wp
		do p=Ham%rg(1),Ham%rg(2)
			tmp=0._wp
			l=Ham%var(p)%nb
			if(size(Ham%var(p)%c,1)==0) then
				if(present(m)) then
					call get_pb(shape(0),m,pb,WA=WA,iA=iA,AW=AW,WAW=WAW)
					tmp=tmp+val*pb
				else
					tmp=tmp+val
				endif
				if(.not.isnan(Ham%var(p)%V)) then
					tmp=tmp*Ham%var(p)%V
				endif
				get_energy=get_energy+tmp
				cycle
			endif
			do n=1,size(latt%nb(l)%bd)
				i=latt%nb(l)%bd(n)%i(1)
				j=latt%nb(l)%bd(n)%i(2)
				sb=latt%nb(l)%bd(n)%sb
				val=Ham%var(p)%val(Ham%var(p)%bd2v(n))*Ham%var(p)%bd(n)
				flag=" "
				do jc=1,size(Ham%var(p)%c,2)
					if(all(merge(sb(1),sb(2),ham%var(p)%c(:,jc)%l=="i")-ham%var(p)%c(:,jc)%sb==0)) then
						if(flag==" ") flag="i"
					elseif(all(merge(sb(2),sb(1),ham%var(p)%c(:,jc)%l=="i")-ham%var(p)%c(:,jc)%sb==0)) then
						if(flag==" ") flag="j"
					else
						cycle
					endif
					dcfg_all(0)=dcfg_(0)+size(Ham%var(p)%c,1)
					dcfg_all(1:dcfg_all(0))=[dcfg_(1:dcfg_(0)),(Hmf%s2H(merge(i,j,ham%var(p)%c(ic,jc)%l==flag),abs(Ham%var(p)%c(ic,jc)%tp))*sign(1,Ham%var(p)%c(ic,jc)%tp),ic=1,size(Ham%var(p)%c,1))]
					if(get_row(dcfg_all(1:dcfg_all(0)),cfg,k,sg,Pr(1:Pr(0)))) then
						if(present(m)) then
							call get_pb(k(1:k(0)),m,pb,WA=WA,iA=iA,AW=AW,WAW=WAW)
						else
							call get_pb(k(1:k(0)),shape(0),pb,WA=WA)
						endif
						if(size(Ham%var(p)%c,1)==2.and.Ham%var(p)%c(1,jc)%l/=Ham%var(p)%c(2,jc)%l) then
							if(Ham%var(p)%c(1,jc)%l=="i") then
								bdc=latt%nb(l)%bd(n)%bdc
							else
								bdc=conjg(latt%nb(l)%bd(n)%bdc)
							endif
						else
							bdc=abs(latt%nb(l)%bd(n)%bdc)
						endif
						tmp=tmp+sg*val*Ham%var(p)%sg(jc)*pb*bdc*exp(jast(cfg,dcfg_all(1:dcfg_all(0))))
					endif
				enddo
			enddo
			if(.not.isnan(Ham%var(p)%V)) then
				tmp=tmp*Ham%var(p)%V
			endif
			get_energy=get_energy+tmp
		enddo
	end function
	subroutine get_O(icfg,iEcfg,WA,iA,dwf,O,opt)
		complex(wp) :: WA(:,:),iA(:,:),dwf(:,:,:),O(:)
		complex(wp) :: beta,tmp
		integer :: icfg(:),iEcfg(:),opt
		integer :: l,i,j,lp
		O=0._wp
		beta=1._wp
		select case(opt)
		case(1)
			!$omp parallel do reduction(+:O)
			do i=1,size(iA,1)
				call gemv(dwf(icfg,iEcfg(i),:),iA(i,:),O,trans="T",beta=beta)
			enddo
			!$omp end parallel do
		case(2)
			!$omp parallel do reduction(+:O)
			do i=1,size(WA,1)
				call gemv(dwf(icfg,i,:),WA(i,:),O,trans="T",beta=beta)
			enddo
			!$omp end parallel do
		end select
	end subroutine
	subroutine get_Spb(cfg,dcfg,Ecfg,nn,pb,WA,iA,AW,WAW)
		integer :: cfg(:),dcfg(:),Ecfg(:),nn(0:,:)
		complex(wp) :: WA(:,:),AW(:,:),WAW(:,:),iA(:,:)
		complex(wp) :: pb
		complex(wp) :: pb2(2),pbs(2)
		integer :: i,l,sg,ct=0
		integer :: k(0:4),m(0:4)
		if(.not.get_row(dcfg,cfg,k,sg,pack([1],[is_project]))) then
			write(*,"(*(i5))")[1:Ns]
			write(*,"(*(i5))")cfg(1:Ns)
			write(*,"(*(i5))")cfg(Ns+1:)
			write(*,"(*(i5))")dcfg
			write(*,"(*(i5))")k(1:k(0))
			error stop "err get_Spb"
		endif
		pbs=0._wp
		!$omp parallel do private(pb2,sg,m) reduction(+:pbs)
		do i=1,ubound(nn,1)
			if(get_row([nn(i,:),-nn(0,:)],Ecfg,m,sg,[integer::])) then
				call get_pb(shape(0),m(1:m(0)),pb2(1),WA=WA,iA=iA,AW=AW,WAW=WAW)
				pbs(1)=pbs(1)+pb2(1)*conjg(pb2(1))
				call get_pb(k(1:k(0)),m(1:m(0)),pb2(2),WA=WA,iA=iA,AW=AW,WAW=WAW)
				pbs(2)=pbs(2)+pb2(2)*conjg(pb2(2))
			endif
		enddo
		!$omp end parallel do
		!write(*,*)real(pbs)
		pb=pbs(2)/pbs(1)
	end subroutine
	subroutine get_overlap(cfg,Ecfg,nn,WA,iA,AW,WAW,Ok,Ek)
		integer :: cfg(:),Ecfg(:),nn(0:,:)
		complex(wp) :: WA(:,:),AW(:,:),WAW(:,:),iA(:,:),Ek(:),Ok(:)
		complex(wp) :: pb
		integer :: i,sg,j
		integer :: m(0:4)
		!$omp parallel do private(pb,sg,m)
		do i=1,ubound(nn,1)
			if(get_row([nn(i,:),-nn(0,:)],Ecfg,m,sg,[integer::])) then
				call get_pb(shape(0),m(1:m(0)),pb,WA=WA,iA=iA,AW=AW,WAW=WAW)
				Ok(i)=pb*sg
				Ek(i)=get_energy(cfg,WA,m=m(1:m(0)),iA=iA,AW=AW,WAW=WAW)*sg
			endif
		enddo
		!$omp end parallel do
	end subroutine
	real(wp) function jast(cfg,dcfg,Oja)
		integer :: cfg(:),dcfg(:)
		complex(wp), optional :: Oja(:)
		integer :: i(0:size(dcfg)),j,plv,lv,ibd,l,l1,l2,l3,ic,jc,nb,k(0:size(dcfg)*10),sg,sb(2),ij(2),m
		real(wp) :: n(2)
		character :: flag
		logical :: is_dh
		jast=0._wp
		if(sum(Hja%var(:)%n)==0) then
			return
		endif
		if(present(Oja)) then
			if(sum(Hja%var(1:)%n)==0) return
			Oja=0._wp
		endif
		! the index of jastrew facter
		plv=0
		! get the list of site that changes
		i(0)=0
		i(1:)=Hmf%H2s(abs(dcfg),1)
		do l=1,size(dcfg)
			if(all(i(1:i(0))/=i(l))) then
				i(0)=i(0)+1
				i(i(0))=i(l)
			endif
		enddo

		do l=Hja%rg(1),Hja%rg(2)
			nb=Hja%var(l)%nb
			is_dh=.false.
			if(allocated(Hja%var(l)%label)) then
				if(Hja%var(l)%label(:)=="dh") then
					is_dh=.true.
					!m=ichar(Hja%var(l)%label(3:3))-48
				endif
			endif
			if(is_dh) then
				do l1=1,i(0)
					do l2=0,size(latt%nb(nb)%st(i(l1))%j)
						ij(1)=merge(i(l1),latt%nb(nb)%st(i(l1))%j(l2),l2==0)
						n=1._wp
						do l3=1,size(latt%nb(nb)%st(ij(1))%j)
							ij(2)=latt%nb(nb)%st(ij(1))%j(l3)
							ibd=latt%nb(nb)%st(ij(1))%bd(l3)
							if(abs(Hja%var(l)%bd(ibd))>1e-6_wp) then
								if(abs(n(1))>1e-6_wp) then
									if(get_row([Hmf%s2H(ij(1),1),[-Hmf%s2H(ij(1),2)]*merge(-1,1,is_ph),dcfg],cfg,k,sg,[integer::])) then
										n(1)=0._wp
									elseif(get_row([-Hmf%s2H(ij(1),1),[Hmf%s2H(ij(1),2)]*merge(-1,1,is_ph),dcfg],cfg,k,sg,[integer::])) then
										n(1)=0._wp
									elseif(get_row([Hmf%s2H(ij(1),1),-Hmf%s2H(ij(2),1),[Hmf%s2H(ij(1),2),-Hmf%s2H(ij(2),2)]*merge(-1,1,is_ph),dcfg],cfg,k,sg,[integer::])) then
										n(1)=0._wp
									elseif(get_row([Hmf%s2H(ij(2),1),-Hmf%s2H(ij(1),1),[Hmf%s2H(ij(2),2),-Hmf%s2H(ij(1),2)]*merge(-1,1,is_ph),dcfg],cfg,k,sg,[integer::])) then
										n(1)=0._wp
									endif
								endif
								if(abs(n(2))>1e-6_wp) then
									if(get_row([Hmf%s2H(ij(1),1),[-Hmf%s2H(ij(1),2)]*merge(-1,1,is_ph)],cfg,k,sg,[integer::])) then
										n(2)=0._wp
									elseif(get_row([-Hmf%s2H(ij(1),1),[Hmf%s2H(ij(1),2)]*merge(-1,1,is_ph)],cfg,k,sg,[integer::])) then
										n(2)=0._wp
									elseif(get_row([Hmf%s2H(ij(1),1),-Hmf%s2H(ij(2),1),[Hmf%s2H(ij(1),2),-Hmf%s2H(ij(2),2)]*merge(-1,1,is_ph)],cfg,k,sg,[integer::])) then
										n(2)=0._wp
									elseif(get_row([Hmf%s2H(ij(2),1),-Hmf%s2H(ij(1),1),[Hmf%s2H(ij(2),2),-Hmf%s2H(ij(1),2)]*merge(-1,1,is_ph)],cfg,k,sg,[integer::])) then
										n(2)=0._wp
									endif
								endif
							endif
						enddo
						lv=Hja%var(l)%bd2v(ibd)
						if(present(Oja).and.l>=1) then
							Oja(plv+lv)=Oja(plv+lv)+(n(1)-n(2))
						else
							jast=jast+Hja%var(l)%val(lv)*(n(1)-n(2))
						endif
					enddo
				enddo
			else
				do l1=1,i(0)
					do l2=1,size(latt%nb(nb)%st(i(l1))%j)
						n=0._wp
						j=latt%nb(nb)%st(i(l1))%j(l2)
						! check whether the nn site is already considerd
						if(any(i(1:l1-1)==j)) cycle
						sb=latt%nb(0)%bd([i(l1),j])%sb(1)
						flag=" "
						!for all terms
						do jc=1,size(Hja%var(l)%c,2)
							!assign the right i and j according to sb, if none, cycle
							if(all(merge(sb(1),sb(2),Hja%var(l)%c(:,jc)%l=="i")-Hja%var(l)%c(:,jc)%sb==0)) then
								if(flag==" ") flag="i"
							elseif(all(merge(sb(2),sb(1),Hja%var(l)%c(:,jc)%l=="i")-Hja%var(l)%c(:,jc)%sb==0)) then
								if(flag==" ") flag="j"
							else
								cycle
							endif
							!count terms for new cfg
							if(get_row([(Hmf%s2H(merge(i(l1),j,Hja%var(l)%c(ic,jc)%l==flag),abs(Hja%var(l)%c(ic,jc)%tp))*sign(1,Hja%var(l)%c(ic,jc)%tp),ic=1,size(Hja%var(l)%c,1)),dcfg],cfg,k,sg,[integer::])) then
								n(1)=n(1)+Hja%var(l)%sg(jc)
							endif
							!count terms for old cfg
							if(get_row([(Hmf%s2H(merge(i(l1),j,Hja%var(l)%c(ic,jc)%l==flag),abs(Hja%var(l)%c(ic,jc)%tp))*sign(1,Hja%var(l)%c(ic,jc)%tp),ic=1,size(Hja%var(l)%c,1))],cfg,k,sg,[integer::])) then
								n(2)=n(2)+Hja%var(l)%sg(jc)
							endif
						enddo
						! get bond index for this nb
						ibd=latt%nb(nb)%st(i(l1))%bd(l2)
						! get value index for this bond
						lv=Hja%var(l)%bd2v(ibd)
						! multiply the bond value for this bond
						n=n*Hja%var(l)%bd(ibd)%re
						if(present(Oja).and.l>=1) then
							Oja(plv+lv)=Oja(plv+lv)+(n(1)-n(2))
						else
							jast=jast+Hja%var(l)%val(lv)*(n(1)-n(2))
						endif
					enddo
				enddo
			endif
			if(l>=0) then
				plv=plv+Hja%var(l)%n
			endif
		enddo
	end function
end module
