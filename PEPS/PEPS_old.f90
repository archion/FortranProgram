include "tensor.f90"
include "list.f90"
module M_fig
	use M_tensor
	use M_list
	use M_matrix
	use M_lattice_test1
	implicit none
	integer, parameter :: d=2,chi=d*10,DD=d**2,p=2,Nstep=400
	type t_figst
		integer :: tp
		integer :: isub
		character(:), allocatable :: label(:)
		integer, allocatable :: bd(:)
		integer, allocatable :: j(:)
	end type
	type t_fig
		integer, allocatable :: bd(:,:)
		integer, allocatable :: bd2st(:,:)
		type(t_figst), allocatable :: st(:)
	end type
	type t_dat
		type(t_tensor) :: T
		type(t_dat), pointer :: link1 => null()
		type(t_dat), pointer :: link2 => null()
		integer(8) :: cost(2)
		!integer(8) :: icost(200)
		integer :: cfg=0
		logical :: flag=.false.
	contains
		final :: delete_dat
	end type
contains
	subroutine delete_dat(self)
		type(t_dat) :: self
		!write(*,*)"delete"
		!if(.not.self%flag) then
		call self%T%clear()
		nullify(self%link1,self%link2)
		!endif
	end subroutine
	function get_dat(self) result(rt)
		type(t_node) :: self
		type(t_dat), pointer :: rt
		associate(p=>self%dat)
			select type(p)
			type is(t_dat)
				rt => p
			end select
		end associate
	end function
	function put_dat(self) result(rt)
		type(t_dat), pointer :: self
		class(*), pointer :: rt
		!allocate(rt,source=self)
		rt => self
	end function
	recursive function fig_contract(que) result(rt)
		type(t_dat) :: que
		type(t_tensor) :: rt
		if(.not.associated(que%T%rc)) then
			if(.not.associated(que%link1%T%rc)) then
				que%link1%T=fig_contract(que%link1)
				que%link1%T%is_return=.true.
			endif
			if(.not.associated(que%link2%T%rc)) then
				que%link2%T=fig_contract(que%link2)
				que%link2%T%is_return=.true.
			endif
			rt=dot(que%link1%T,que%link2%T)
		else
			rt=que%T
		endif
		rt%is_return=.true.
	end function
	function merge_node(dat1,dat2,dat) result(flag)
		type(t_dat), pointer :: dat1,dat2
		type(t_dat), pointer :: dat
		logical :: flag
		logical, allocatable :: flag1(:),flag2(:)
		integer :: cost,i,j
		integer(8), allocatable :: shap(:)
		flag=.true.
		dat%cfg=ior(dat1%cfg,dat2%cfg)
		if(popcnt(dat%cfg)==(popcnt(dat1%cfg)+popcnt(dat2%cfg))) then
			allocate(flag1(size(dat1%T%shap)),flag2(size(dat2%T%shap)))
			cost=1
			flag1=.true.
			flag2=.true.
			do i=1,size(dat1%T%shap)
				do j=1,size(dat2%T%shap)
					if(dat1%T%label(i)==dat2%T%label(j)) then
						flag1(i)=.false.
						flag2(j)=.false.
						cost=cost*dat1%T%shap(i)
					endif
				enddo
			enddo
			call allocate(dat%T%shap,[pack(dat1%T%shap,flag1),pack(dat2%T%shap,flag2)])
			allocate(shap(size(dat%T%shap)))
			shap=dat%T%shap
			dat%T%label=[pack(dat1%T%label,flag1),pack(dat2%T%label,flag2)]
			dat%cost(1)=product(shap)
			dat%cost(2)=max(dat%cost(2),dat%cost(1)+dat1%cost(2)+dat2%cost(2))
			dat%cost(1)=dat%cost(1)*cost+dat1%cost(1)+dat2%cost(1)
			dat%link1 => dat1
			dat%link2 => dat2
			dat%flag=.false.
		else
			flag=.false.
		endif
	end function
	subroutine contract_order(que)
		type(t_list) :: que
		type(t_list) :: S(que%num)
		integer :: n,i,j,c,d,k,minx
		integer(8) :: u(2),uc,u0
		type(t_dat), pointer :: dat,dat1,dat2,dat3
		type(t_node), pointer :: pos1,pos2,pos3
		integer :: cm=1
		n=que%num
		minx=10000
		pos1 => que%head
		do i=1,n
			if(i/=1) pos1 => pos1%next
			dat => get_dat(pos1)
			dat%cfg=ibset(0,i-1)
			if(minx>minval(dat%T%shap)) minx=minval(dat%T%shap)
			dat%cost=[0,product(dat%T%shap)]
			call S(1)%insert(put_dat(dat))
		enddo
		do i=1,n
			pos1 => S(i)%head
			do j=1,S(i)%num
				if(j/=1) pos1 => pos1%next
				dat => get_dat(pos1)
				dat%flag=.false.
			enddo
		enddo
		uc=1; u(1)=0
		do while(S(n)%num==0)
			u(2)=100000000000
			do c=2,n
				do d=1,c/2
					pos1 => S(d)%head
					do i=1,S(d)%num
						if(i/=1) pos1 => pos1%next
						dat1 => get_dat(pos1)
						pos2 => S(c-d)%head
						do j=1,S(c-d)%num
							if(j/=1) pos2 => pos2%next
							if((d==(c-d)).and.(i>=j)) cycle
							dat2 => get_dat(pos2)
							nullify(dat)
							allocate(dat)
							if(merge_node(dat1,dat2,dat)) then
								if(dat1%flag.or.dat2%flag) then
									u0=0
								else
									u0=u(1)
								endif
								if(dat%cost(cm)>uc.and.dat%cost(cm)<u(2)) then
									u(2)=dat%cost(cm)
								elseif(dat%cost(cm)<=uc.and.dat%cost(cm)>u0) then
									pos3 => S(c)%head
									do k=1,S(c)%num
										if(k/=1) pos3 => pos3%next
										dat3 => get_dat(pos3)
										if(dat%cfg==dat3%cfg) then
											if(dat%cost(cm)<dat3%cost(cm)) then
												dat3%flag=.true.
												dat3%cost=dat%cost
												dat3%link1 => dat%link1
												dat3%link2 => dat%link2
											endif
											exit
										endif
									enddo
									if(k==S(c)%num+1) then
										call S(c)%insert(put_dat(dat))
										cycle
									endif
								endif
							endif
							deallocate(dat)
						enddo
					enddo
				enddo
			enddo
			u(1)=uc
			uc=max(u(2),minx*uc)
			do i=1,n
				pos1 => S(i)%head
				do j=1,S(i)%num
					if(j/=1) pos1 => pos1%next
					dat1 => get_dat(pos1)
					dat1%flag=.false.
				enddo
			enddo
		enddo
		pos1 => S(1)%head
		do i=1,S(1)%num
			if(i/=1) pos1 => pos1%next
			dat1 => get_dat(pos1)
			dat1%flag=.true.
		enddo
		dat => get_dat(S(n)%head)
		!write(*,*)"the cost is: ",dat%cost
		call get_que(que,dat)

		do i=1,n
			pos1 => S(i)%head
			do j=1,S(i)%num
				if(j/=1) pos1 => pos1%next
				dat => get_dat(pos1)
				if(dat%flag) nullify(pos1%dat)
			enddo
			!call delete_list(S(i))
		enddo
	end subroutine
	recursive subroutine get_que(que,dat)
		type(t_list) :: que
		type(t_dat), pointer :: dat
		if(dat%flag) return
		if(.not.dat%link1%flag) then
			call get_que(que,dat%link1)
		endif
		if(.not.dat%link2%flag) then
			call get_que(que,dat%link2)
		endif
		dat%flag=.true.
		call que%insert(put_dat(dat))
	end subroutine
	subroutine gen_fig(fig,bd,st)
		type(t_fig) :: fig
		integer :: bd(:),st(:)
		integer :: Ns,n,i
		Ns=size(st)/2
		allocate(fig%bd(size(bd)/3,2),fig%bd2st(size(bd)/3,2),fig%st(Ns))
		fig%st(:)%isub=st(1::2)
		fig%st(:)%tp=st(2::2)
		do i=1,Ns
			fig%st(i)%label=[character::]
		enddo
		do n=1,size(fig%bd,1)
			fig%bd(n,:)=bd((n-1)*3+2:n*3)
			do i=1,2
				if(.not.allocated(fig%st(fig%bd(n,i))%bd)) then
					select case(fig%st(fig%bd(n,i))%tp)
					case(0)
						allocate(fig%st(fig%bd(n,i))%bd(4),fig%st(fig%bd(n,i))%j(4))
					case(1:4)
						allocate(fig%st(fig%bd(n,i))%bd(3),fig%st(fig%bd(n,i))%j(3))
					case(5:8)
						allocate(fig%st(fig%bd(n,i))%bd(2),fig%st(fig%bd(n,i))%j(2))
					end select
				endif
			enddo
			select case(bd((n-1)*3+1))
			case(-1)
				fig%st(fig%bd(n,1))%label=[fig%st(fig%bd(n,1))%label,"r"]
				fig%st(fig%bd(n,2))%label=[fig%st(fig%bd(n,2))%label,"l"]
			case(-2)
				fig%st(fig%bd(n,1))%label=[fig%st(fig%bd(n,1))%label,"d"]
				fig%st(fig%bd(n,2))%label=[fig%st(fig%bd(n,2))%label,"u"]
			end select
			fig%bd2st(n,:)=[size(fig%st(fig%bd(n,1))%label),size(fig%st(fig%bd(n,2))%label)]
			fig%st(fig%bd(n,1))%bd(size(fig%st(fig%bd(n,1))%label))=n
			fig%st(fig%bd(n,2))%bd(size(fig%st(fig%bd(n,2))%label))=n
			fig%st(fig%bd(n,1))%j(size(fig%st(fig%bd(n,1))%label))=2
			fig%st(fig%bd(n,2))%j(size(fig%st(fig%bd(n,2))%label))=1
		enddo
	end subroutine
	subroutine insert_que(a,que)
		type(t_tensor) :: a
		type(t_list) :: que
		type(t_dat), pointer :: dat
		allocate(dat)
		dat%T=a
		call que%insert(put_dat(dat))
	end subroutine
	function cut_label(st,cut) result(rt)
		type(t_figst) :: st
		integer :: cut(:)
		integer :: i
		character(:), allocatable :: rt(:)
		rt=[character::]
		do i=1,size(st%bd)
			if(any(st%bd(i)==cut)) then
				rt=[rt,trim(st%label(i))//to_char(st%bd(i))]
			else
				rt=[rt,to_char(st%bd(i))]
			endif
		enddo
	end function
	subroutine ini_que(fig,cut,rm,aa,que)
		type(t_fig) :: fig
		integer :: cut(:),rm(:)
		type(t_tensor) :: aa(:,0:)
		type(t_list) :: que
		!type(t_dat), pointer :: dat
		integer :: i,n
		if(associated(que%head)) then
			call delete_list(que)
		endif
		do i=1,size(fig%st)
			if(any(i==rm)) cycle
			!allocate(dat)
			!dat%T=aa(fig%st(i)%isub,fig%st(i)%tp)%change_label(fig%st(i)%label,to_char(fig%st(i)%bd))
			!do n=1,size(fig%st(i)%bd)
				!if(any(fig%st(i)%bd(n)==cut).or.any(fig%bd(fig%st(i)%bd(n),fig%st(i)%j(n))==rm)) then
					!dat%T=dat%T%change_label([to_char(fig%st(i)%bd(n))],&
						![trim(fig%st(i)%label(n))//to_char(fig%st(i)%bd(n))])
				!endif
			!enddo
			!call que%insert(put_dat(dat))
			call insert_que(aa(fig%st(i)%isub,fig%st(i)%tp)%change_label(fig%st(i)%label,cut_label(fig%st(i),cut)),que)
		enddo
	end subroutine
	subroutine CTM(fig,pair,aa)
		type(t_fig) :: fig
		integer :: pair(:,:)
		type(t_tensor) :: aa(:,0:)
		type(t_tensor) :: tmp(size(pair)),a(2),P(size(pair,1)-1),E
		integer :: n,i,j,cut(2,size(pair,1)-1),m,k,i2(2),n2,st(2),c(2)
		integer :: pair_(size(pair))
		character(:), allocatable :: label(:),lbd(:)
		complex(8), allocatable :: U(:,:),V(:,:)
		real(8), allocatable :: s(:)
		type(t_list) :: que
		pair_=reshape(pair,shape(pair_))
		do n=1,size(pair)
			st=fig%bd(pair_(n),:)
			if(mod(n-1,size(pair,1))+1==1) then
				cut=0
				m=0
			endif

			do k=1,2
				a(k)=aa(fig%st(st(k))%isub,fig%st(st(k))%tp)%change_label(&
					[fig%st(st(k))%label(fig%bd2st(pair_(n),k))],[to_char(pair_(n))])
			enddo

			if(size(a(1)%shap)>size(a(2)%shap)) then
				c=[2,1]
			else
				c=[1,2]
			endif
			st=st(c)

			n2=0
			do i=1,size(fig%st(st(1))%bd)
				label=[fig%st(st(1))%label(i)]
				if(fig%st(st(1))%bd(i)/=pair_(n)) then
					do j=1,size(fig%st(st(2))%bd)
						if(label(1)==fig%st(st(2))%label(j)) then
							lbd=[label(1),fig%st(fig%bd(fig%st(st(1))%bd(i),fig%st(st(1))%j(i)))%label(fig%bd2st(fig%st(st(1))%bd(i),fig%st(st(1))%j(i)))]
							label=[to_char(fig%st(st(1))%bd(i)),to_char(fig%st(st(2))%bd(j))]
							a(c(1))=a(c(1))%change_label([lbd(1)],[trim(lbd(1))//label(1)])
							a(c(2))=a(c(2))%change_label([lbd(1)],[trim(lbd(1))//label(2)])
							do k=1,m
								if(all(cut(:,k)-[fig%st(st(1))%bd(i),fig%st(st(2))%bd(j)]==0)) then
									n2=n2+1
									P(k)=set_conjg(P(k)%change_label(trim(lbd(2))//[label(1),label(2),''],trim(lbd(1))//[label(1),label(2),'']))
									i2(n2)=k
									exit
								endif
							enddo
							if(k==m+1) then
								m=m+1
								n2=n2+1
								cut(:,m)=[fig%st(st(1))%bd(i),fig%st(st(2))%bd(j)]
								call ini_que(fig,cut(:,m),[integer::],aa,que)
								call contract_order(que)
								E=fig_contract(get_dat(que%tail))
								call E%svd(trim(lbd(1))//[label(1),label(2)],trim(lbd(2))//[label(1),label(2)],U,s,V)
								call P(m)%new(trim(lbd(1))//[label(1),label(2),''],[E%shap(1:2),min(E%shap(1)*E%shap(2),chi)])
								call P(m)%get_tensor([P(m)%label],U(:,:min(E%shap(1)*E%shap(2),chi)))
								P(m)=set_conjg(P(m))
								i2(n2)=m
								call E%clear()
							endif
						endif
					enddo
				endif
			enddo

			if(n2==1) then
				tmp(n)=dot(dot(a(1),a(2)),P(i2(1)))
				!write(*,*)a(1)%label,"|",a(2)%label,"|",P(i2(1))%label,P(i2(1))%is_conjg
			else
				tmp(n)=dot(dot(dot(a(1),a(2)),P(i2(1))),P(i2(2)))
				!write(*,*)a(1)%label,"|",a(2)%label,"|",P(i2(1))%label,"|",P(i2(2))%label,P(i2(1:2))%is_conjg
			endif
			call a(1)%clear()
			call a(2)%clear()
			if(mod(n,size(pair,1))==0) then
				do i=1,size(P)
					call P(i)%clear()
				enddo
			endif
		enddo

		do n=1,size(pair)
			tmp(n)%is_return=.true.
			st=fig%bd(pair_(n),:)
			if(size(aa(fig%st(st(1))%isub,fig%st(st(1))%tp)%shap)>size(aa(fig%st(st(2))%isub,fig%st(st(2))%tp)%shap)) st([1,2])=st([2,1])
			aa(fig%st(st(2))%isub,fig%st(st(1))%tp)=tmp(n)
		enddo
		!write(*,*)"CTM finished"
	end subroutine
	subroutine simp_update(H,delta0,A,n)
		complex(8) :: H(:,:)
		real(8) :: delta0(:)
		type(t_tensor) :: A(:)
		integer :: n
		type(t_tensor) :: App,b(2)
		real(8) :: ld(d,2*size(A)),alpha
		type(t_figst) :: st(size(A))
		integer :: ist(2),i,j,step
		real(8) :: ld1(d**3),ld2(d**3),ld0(d**4),He(size(H,1)),delta
		real(8), allocatable :: s(:)
		complex(8), allocatable :: U(:,:),V(:,:),Q1(:,:),Q2(:,:),R(:,:),L(:,:)
		complex(8) :: HU(size(H,1),size(H,2)),Hexp(size(H,1),size(H,2))
		integer :: ord1(4),ord2(4)
		character(:), allocatable :: lgate

		call random_number(ld)

		HU=H
		call heev(HU,He,"V")


		st(1)%label=["u","r","d","l"]
		call allocate(st(1)%bd,[1,2,3,4])
		st(2)%label=["u","r","d","l"]
		call allocate(st(2)%bd,[3,4,1,2])

		do i=1,size(A)
			ld0=1d0/sqrt(outprod(outprod(ld(:,st(i)%bd(1)),ld(:,st(i)%bd(2))),outprod(ld(:,st(i)%bd(3)),ld(:,st(i)%bd(4)))))
			A(i)=dot(cmplx(ld0,kind=8),A(i),st(i)%label)
		enddo

		ist=[1,2]
		call b(1)%new(["l1","p1","r1"],[d*p,p,d])
		call b(2)%new(["l2","p2","r2"],[d,p,d*p])
		do step=1,n
			!delta=(delta0(1)-delta0(2))*(1d0/(1d0+exp(8d0*step/n-4d0))-1d0/(1d0+exp(4d0)))/(1d0/(1d0+exp(-4d0))-1d0/(1d0+exp(4d0)))+delta0(2)
			!delta=max(delta0(1)-(delta0(1)-delta0(2))/(n/2d0)*step,delta0(2))
			delta=delta0(1)*(exp(-6d0*step/n)-exp(-6d0))+delta0(2)
			Hexp=matmul(matmul(HU,diag(exp(-He*delta))),transpose(conjg(HU)))
			do i=1,size(ld,2)

				ord1=[1:4]
				ord2=[1:4]
				ist=[1,2]
				do j=1,size(st(1)%bd)
					if(st(1)%bd(j)==i) then
						if(st(1)%label(j)/="r".and.st(1)%label(j)/="d") then
							ist=[2,1]
						endif
						exit
					endif
				enddo
				do j=1,3
					if(st(ist(1))%bd(j)==i) then
						ord1([j,4])=ord1([4,j])
					endif
					if(st(ist(2))%bd(j)==i) then
						ord2([j,4])=ord2([4,j])
					endif
				enddo
				ld1=outprod(outprod(ld(:,st(ist(1))%bd(ord1(1))),ld(:,st(ist(1))%bd(ord1(2)))),ld(:,st(ist(1))%bd(ord1(3))))
				ld2=outprod(outprod(ld(:,st(ist(2))%bd(ord2(1))),ld(:,st(ist(2))%bd(ord2(2)))),ld(:,st(ist(2))%bd(ord2(3))))
				A(ist(1))=dot(cmplx(ld1,kind=8),A(ist(1)),st(ist(1))%label(ord1(1:3)))
				A(ist(2))=dot(A(ist(2)),cmplx(ld2,kind=8),st(ist(2))%label(ord2(1:3)))

				call A(ist(1))%qr([st(ist(1))%label(ord1(1:3))],["p",st(ist(1))%label(ord1(4))],Q1,R)
				call b(1)%get_tensor(["l1","p1","r1"],R)
				call A(ist(2))%lq([st(ist(2))%label(ord2(4)),"p"],[st(ist(2))%label(ord2(1:3))],Q2,L)
				call b(2)%get_tensor(["l2","p2","r2"],L)

				App=dot(dot(dot(b(1),cmplx(ld(:,i),kind=8),["r1"]),b(2),["r1","l2"]),Hexp,["p1","p2"])
				call App%svd(["l1","p1"],["p2","r2"],U,s,V)
				ld(:,i)=s(:d)/sqrt(sum(s(:d)**2))
				call b(1)%get_tensor(["l1","p1","r1"],U(:,:d))
				call b(2)%get_tensor(["l2","p2","r2"],V(:d,:))
				do j=1,size(ld1)
					Q1(j,:)=1d0/ld1(j)*Q1(j,:)
					Q2(:,j)=1d0/ld2(j)*Q2(:,j)
				enddo
				A(ist(1))=split_label(change_label(dot(Q1,b(1),["l1"]),["p1","r1"],["p",st(ist(1))%label(ord1(4))]),["l1"],st(ist(1))%label(ord1(1:3)))
				A(ist(2))=split_label(change_label(dot(b(2),Q2,["r2"]),["p2","l2"],["p",st(ist(2))%label(ord2(4))]),["r2"],st(ist(2))%label(ord2(1:3)))
			enddo
		enddo
		ld=sqrt(ld)
		do i=1,size(A)
			ld0=outprod(outprod(ld(:,st(i)%bd(1)),ld(:,st(i)%bd(2))),outprod(ld(:,st(i)%bd(3)),ld(:,st(i)%bd(4))))
			A(i)=dot(cmplx(ld0,kind=8),A(i),st(i)%label)
		enddo
		call App%clear()
		call b(1)%clear()
		call b(2)%clear()
	end subroutine
	subroutine full_update(E,Hexp,lbd,b)
		type(t_tensor) :: E,b(:)
		complex(8) :: Hexp(:,:)
		character(*) :: lbd(:)
		type(t_tensor) :: Q,Qt,T,tmp,Z
		complex(8), allocatable :: MS(:,:),MR(:,:),Mtmp(:,:),MZ(:,:),MQ(:,:),R(:,:),iR(:,:),L(:,:),iL(:,:),U(:,:),V(:,:)
		real(8), allocatable :: Eg(:),s(:)
		integer :: i,j,n
		complex(8) :: cost,cost0,pcost

		call E%get_mat(lbd,lbd//"'",MZ)
		MZ=0.5d0*(MZ+transpose(conjg(MZ)))
		call allocate(Eg,size(MZ,1))
		call heev(MZ,Eg,"V")
		do i=1,size(MZ,2)
			if(Eg(i)>0d0) then
				MZ(:,i)=MZ(:,i)*sqrt(Eg(i))
			else
				MZ(:,i)=0d0
			endif
		enddo

		!!! Gauge Fixing
		call Z%new([lbd,"p"],[E%shap(1:2),size(MZ,2)])
		call Z%get_tensor([lbd,"p"],MZ)
		call Z%qr([lbd(1),"p"],[lbd(2)],MQ,R)
		call Z%lq([lbd(1)],[lbd(2),"p"],MQ,L)
		call allocate(iL,L)
		call mat_inv(iL)
		call allocate(iR,R)
		call mat_inv(iR)
		b(1)=dot(b(1),L,[lbd(2)])
		b(2)=dot(R,b(2),[lbd(1)])
		Z=dot(dot(iL,Z,[lbd(1)]),iR,[lbd(2)])
		!call E%get_tensor(E%label,matmul(MZ,transpose(conjg(MZ))))
		E=dot(Z,set_conjg(Z%change_label(lbd,lbd//"'")))
		call Z%clear()

		T=dot(dot(b(1)%change_label(["p"],["p1"]),b(2)%change_label(["p"],["p2"]),lbd),Hexp,["p1","p2"])
		Q=dot(E,T,[lbd,lbd(2:1:-1)])
		Qt=dot(E,set_conjg(T%change_label(lbd,lbd//"'")),[lbd//"'",lbd(2:1:-1)//"'"])
		call T%svd([lbd(2),"p1"],[lbd(1),"p2"],U,s,V)
		do i=1,d
			U(:,i)=U(:,i)*sqrt(s(i))
			V(i,:)=V(i,:)*sqrt(s(i))
		enddo
		call b(1)%get_tensor([lbd(2),"p",lbd(1)],U(:,:d))
		call b(2)%get_tensor([lbd(2),lbd(1),"p"],V(:d,:))
		T=dot(E,dot(T,set_conjg(T%change_label(lbd,lbd//"'"))),[lbd,lbd(2:1:-1),lbd//"'",lbd(2:1:-1)//"'"])
		cost0=T%rc%T(1)
		pcost=-1d0
		n=0
		!write(*,*)"cost: "
		do 
			n=n+1
			cost=cost0
			tmp=dot(b(1)%change_label(["p"],["p1"]),b(2)%change_label(["p"],["p2"]),lbd)
			T=dot(Q,set_conjg(tmp%change_label(lbd,lbd//"'")),["p1","p1","p2","p2",lbd//"'",lbd(2:1:-1)//"'"])
			cost=cost-T%rc%T(1)
			T=dot(Qt,tmp,["p1","p1","p2","p2",lbd,lbd(2:1:-1)])
			cost=cost-T%rc%T(1)
			T=dot(E,dot(tmp,set_conjg(tmp%change_label(lbd,lbd//"'"))),[lbd,lbd(2:1:-1),lbd//"'",lbd(2:1:-1)//"'"])
			cost=cost+T%rc%T(1)
			!write(*,"(2es14.6','2es14.6)")T%rc%T(1),cost
			!read(*,*)
			if(abs(cost-pcost)<1d-10) then
				write(*,*)"update finished",cost,n
				exit
			endif
			pcost=cost
			do i=1,2
				j=mod(i,2)+1
				call get_mat(dot(Q,set_conjg(b(j)%change_label([lbd],lbd//"'")),[character(5)::"p"//to_char(j),"p",lbd(j)//"'",lbd(i)//"'"]),lbd//"'",[character::],MS)
				call get_mat(dot(E,dot(b(j),set_conjg(b(j)%change_label([lbd],lbd//"'"))),[character(5)::lbd(j),lbd(i),lbd(j)//"'",lbd(i)//"'"]),lbd//"'",lbd,MR)
				call allocate(MZ,MR)
				call mat_inv(MR)
				cost=sum(abs(matmul(MR,MZ)-diag(1d0,size(MR,1))))
				if(abs(cost)>1d-5) then
					write(*,*)"inverse: ",cost
					stop
				endif
				call b(i)%get_tensor([lbd(2),lbd(1),"p"],matmul(MR,MS))
			enddo
		enddo
		b(1)=dot(b(1),iL,[lbd(2)])
		b(2)=dot(iR,b(2),[lbd(1)])
	end subroutine
	subroutine rescale(fig,aa,A)
		type(t_fig) :: fig
		type(t_tensor) :: aa(:,0:),A(:)
		type(t_list) :: que
		integer :: i,loc(size(fig%st)),ne
		complex(8) :: scale,cp
		real(8) :: norm
		type(t_tensor) :: E
		logical :: flag(size(aa,1),0:size(aa,2)-1)
		call ini_que(fig,[integer::],[integer::],aa,que)
		call contract_order(que)
		E=fig_contract(get_dat(que%tail))
		norm=1d0
		ne=0
		do i=1,size(fig%st)
			loc(i:i)=maxloc(abs(aa(fig%st(i)%isub,fig%st(i)%tp)%rc%T))
			norm=norm*(abs(aa(fig%st(i)%isub,fig%st(i)%tp)%rc%T(loc(i))))
			if(fig%st(i)%tp/=0) then
				ne=ne+1
			endif
		enddo
		cp=(conjg(E%rc%T(1))/abs(E%rc%T(1)))**(1d0/ne)
		norm=(norm/abs(E%rc%T(1)))**(1d0/16d0)
		flag=.false.
		do i=1,size(fig%st)
			if(flag(fig%st(i)%isub,fig%st(i)%tp)) then
				cycle
			endif
			flag(fig%st(i)%isub,fig%st(i)%tp)=.true.

			scale=norm/abs(aa(fig%st(i)%isub,fig%st(i)%tp)%rc%T(loc(i)))
			if(fig%st(i)%tp==0) then
				A(fig%st(i)%isub)%rc%T=A(fig%st(i)%isub)%rc%T*scale**0.5d0
			else
				scale=scale*cp
			endif
			aa(fig%st(i)%isub,fig%st(i)%tp)%rc%T=aa(fig%st(i)%isub,fig%st(i)%tp)%rc%T*scale
			if(fig%st(i)%tp==0) then
				E=merge_label(dot(A(fig%st(i)%isub),set_conjg(A(fig%st(i)%isub)%change_label(["u","d","l","r"],["u","d","l","r"]//"'"))),["u","u'","d","d'","l","l'","r","r'"],["u","d","l","r"])
				call E%reorder(aa(fig%st(i)%isub,fig%st(i)%tp)%label)
				!write(*,*)"check: ",sum(abs(E%rc%T-aa(fig%st(i)%isub,fig%st(i)%tp)%rc%T))
			endif
		enddo
		call E%clear()
	end subroutine
	function get_vsite(fig,A,aa,st,H,sg) result(rt)
		type(t_fig) :: fig
		type(t_tensor) :: A(:),aa(:,0:)
		integer :: st(:)
		complex(8) :: H(:,:)
		real(8), optional :: sg(:)
		type(t_list) :: que
		complex(8) :: rt
		integer :: i,n
		type(t_tensor) :: E
		character(:), allocatable :: label(:)
		type(t_dat), pointer :: dat
		label=["u","u'","r","r'","d","d'","l","l'"]
		rt=0d0
		do n=1,size(st)
			do i=1,size(fig%st)
				if(i==st(n)) then
					call insert_que(change_label(merge_label(&
						dot(dot(A(fig%st(i)%isub),H,["p"]),set_conjg(A(fig%st(i)%isub)%change_label([label(::2)],[label(2::2)])))&
						,label,label(::2)),fig%st(i)%label,cut_label(fig%st(i),[integer::])),que)
				else
					call insert_que(aa(fig%st(i)%isub,fig%st(i)%tp)%change_label(fig%st(i)%label,cut_label(fig%st(i),[integer::])),que)
				endif
			enddo
			call contract_order(que)
			E=fig_contract(get_dat(que%tail))
			if(present(sg)) then
				rt=rt+real(E%rc%T(1))*sg(i)
			else
				rt=rt+real(E%rc%T(1))
			endif
			call delete_list(que)
		enddo
		call E%clear()
	end function
	function get_vbond(fig,A,aa,bd,H,sg) result(rt)
		type(t_fig) :: fig
		type(t_tensor) :: A(:),aa(:,0:)
		integer :: bd(:)
		complex(8) :: H(:,:)
		real(8), optional :: sg(:)
		type(t_list) :: que
		complex(8) :: rt
		integer :: i,n,k
		type(t_tensor) :: E
		character(:), allocatable :: label(:),lgate
		label=["u","u'","r","r'","d","d'","l","l'"]
		rt=0d0
		do n=1,size(bd)
			do i=1,size(fig%st)
				if(any(i==fig%bd(bd(n),:))) then
					call insert_que(change_label(merge_label(&
						dot(A(fig%st(i)%isub)%change_label(["p"],["p"//to_char(i)]),set_conjg(A(fig%st(i)%isub)%change_label([label(::2),"p"],[label(2::2),"P"//to_char(i)])))&
						,label,label(::2)),fig%st(i)%label,cut_label(fig%st(i),[integer::])),que)
				else
					call insert_que(aa(fig%st(i)%isub,fig%st(i)%tp)%change_label(fig%st(i)%label,cut_label(fig%st(i),[integer::])),que)
				endif
			enddo
			call contract_order(que)
			E=fig_contract(get_dat(que%tail))
			E=contract(dot(E,H,"p"//[to_char(fig%bd(bd(n),1)),to_char(fig%bd(bd(n),2))]),[["p","P"]//to_char(fig%bd(bd(n),1)),["p","P"]//to_char(fig%bd(bd(n),2))])
			if(present(sg)) then
				rt=rt+real(E%rc%T(1))*sg(i)
			else
				rt=rt+real(E%rc%T(1))
			endif
			call delete_list(que)
		enddo
		call E%clear()
	end function
	!function get_vsite(fig,A,aa,st,H,sg) result(rt)
		!type(t_fig) :: fig
		!type(t_tensor) :: A(:),aa(:,0:)
		!integer :: st(:)
		!complex(8) :: H(:,:)
		!real(8), optional :: sg(:)
		!type(t_list) :: que
		!complex(8) :: rt
		!integer :: i,n
		!type(t_tensor) :: E,tmp
		!character(:), allocatable :: label(:)
		!integer :: flag(size(st)*4)
		!label=["u","u'","r","r'","d","d'","l","l'"]
		!rt=0d0
		!n=0
		!do i=1,size(fig%st)
			!if(any(i==st)) then
				!call insert_que(change_label(merge_label(dot(A(fig%st(i)%isub)%change_label(["p"],["p"//to_char(i)]),set_conjg(A(fig%st(i)%isub)%change_label([label(::2),"p"],[label(2::2),"P"//to_char(i)]))),label,label(::2)),fig%st(i)%label,cut_label(fig%st(i),[integer::])),que)
				!flag(n+1:n+2)=[i,-i]
				!n=n+2
			!else
				!call insert_que(aa(fig%st(i)%isub,fig%st(i)%tp)%change_label(fig%st(i)%label,cut_label(fig%st(i),[integer::])),que)
			!endif
		!enddo
		!call contract_order(que)
		!E=fig_contract(get_dat(que%tail))
		!do i=1,size(st)
			!tmp=contract(E,pack(merge("p","P",flag(:n)>0)//to_char(abs(flag(:n))),abs(flag(:n))/=st(i)))
			!tmp=contract(dot(tmp,H,"p"//[to_char(st(i))]),[["p","P"]//to_char(st(i))])
			!if(present(sg)) then
				!rt=rt+real(tmp%rc%T(1))*sg(i)
			!else
				!rt=rt+real(tmp%rc%T(1))
			!endif
		!enddo
		!call E%clear()
		!call tmp%clear()
	!end function
	!function get_vbond(fig,A,aa,bd,H,sg) result(rt)
		!type(t_fig) :: fig
		!type(t_tensor) :: A(:),aa(:,0:)
		!integer :: bd(:)
		!complex(8) :: H(:,:)
		!real(8), optional :: sg(:)
		!type(t_list) :: que
		!complex(8) :: rt
		!integer :: i,n
		!type(t_tensor) :: E,tmp
		!character(:), allocatable :: label(:)
		!integer :: flag(size(bd)*4)
		!label=["u","u'","r","r'","d","d'","l","l'"]
		!rt=0d0
		!n=0
		!do i=1,size(fig%st)
			!if(any(i==fig%bd(bd,:))) then
				!call insert_que(change_label(merge_label(dot(A(fig%st(i)%isub)%change_label(["p"],["p"//to_char(i)]),set_conjg(A(fig%st(i)%isub)%change_label([label(::2),"p"],[label(2::2),"P"//to_char(i)]))),label,label(::2)),fig%st(i)%label,cut_label(fig%st(i),[integer::])),que)
				!flag(n+1:n+2)=[i,-i]
				!n=n+2
			!else
				!call insert_que(aa(fig%st(i)%isub,fig%st(i)%tp)%change_label(fig%st(i)%label,cut_label(fig%st(i),[integer::])),que)
			!endif
		!enddo
		!call contract_order(que)
		!E=fig_contract(get_dat(que%tail))
		!!tmp=E%clone()
		!!call E%reorder(["p"//to_char(flag(:n:2)),"P"//to_char(flag(:n:2))])
		!!call tmp%reorder(["P"//to_char(flag(:n:2)),"p"//to_char(flag(:n:2))])
		!!E%rc%T=0.5d0*(E%rc%T+conjg(tmp%rc%T))
		!do i=1,size(bd)
			!tmp=contract(E,pack(merge("p","P",flag(:n)>0)//to_char(abs(flag(:n))),abs(flag(:n))/=fig%bd(bd(i),1).and.abs(flag(:n))/=fig%bd(bd(i),2)))
			!tmp=contract(dot(tmp,H,"p"//[to_char(fig%bd(bd(i),1)),to_char(fig%bd(bd(i),2))]),[["p","P"]//to_char(fig%bd(bd(i),1)),["p","P"]//to_char(fig%bd(bd(i),2))])
			!if(present(sg)) then
				!rt=rt+real(tmp%rc%T(1))*sg(i)
			!else
				!rt=rt+real(tmp%rc%T(1))
			!endif
		!enddo
		!call E%clear()
		!call tmp%clear()
	!end function
end module


program main
	use M_fig
	implicit none
	real(8) :: delta=0.05d0,alpha=10d0,g=3.10d0
	complex(8), allocatable :: HU(:,:),Q(:,:),R(:,:),L(:,:),H(:,:)
	real(8) :: He(p**2)
	real(8), allocatable :: s(:),ps(:)
	complex(8) :: Hexp(size(He),size(He)),energy,Hsx(p,p),Hsy(p,p),Hsz(p,p),sx,sy,sz
	type(t_tensor), target :: A(2),aa(2,0:9),X(2),b(2),tmp1,tmp2,tH
	type(t_tensor) :: E
	type(t_fig) :: fig,figx,figy
	type(t_list) :: que
	type(t_dat), pointer :: dat
	integer :: flag(8)=[1,1,2,2,3,3,4,4]
	character(:), allocatable :: label(:),lbd(:)
	integer :: step,n,i,j,k,bd

	call random_seed()
	call omp_set_nested(.false.)
	!call omp_set_max_active_levels(1)
	call omp_set_dynamic(.false.)
	call mkl_set_dynamic(0)
	call omp_set_schedule(omp_sched_static,0)

	call mkl_set_num_threads(1)
	call omp_set_num_threads(mkl_get_max_threads())

	! Heisenberg
	call tmp1%new(["u0","d0","r"],[p,p,3],flag="0")
	call tmp2%new(["u1","d1","l"],[p,p,3],flag="0")
	tmp1%rc%T(tmp1%get_idx([1,1,1, 2,2,1, 1,2,2, 2,1,3]))=[1d0,-1d0,1d0,1d0]/2d0
	tmp2%rc%T(tmp2%get_idx([1,1,1, 2,2,1, 2,1,2, 1,2,3]))=[1d0,-1d0,2d0,2d0]/2d0
	Hsx=0d0;Hsx(1,2)=1d0;Hsx(2,1)=1d0
	Hsy=0d0;Hsy(1,2)=img;Hsy(2,1)=-img
	Hsz=0d0;Hsz(1,1)=1d0;Hsz(2,2)=-1d0

	!! ising
	!call tmp1%new(["u0","d0","r"],[p,p,3],flag="0")
	!call tmp2%new(["u1","d1","l"],[p,p,3],flag="0")
	!tmp1%rc%T(tmp1%get_idx([1,1,1, 2,2,1, 1,1,2, 2,2,2, 1,2,3, 2,1,3]))=[1d0,1d0,-1d0,1d0,-g*0.25d0,-g*0.25d0]
	!tmp2%rc%T(tmp2%get_idx([1,2,1, 2,1,1, 1,1,2, 2,2,2, 1,1,3, 2,2,3]))=[-g*0.25d0,-g*0.25d0,1d0,-1d0,1d0,1d0]

	tH=dot(tmp1,tmp2,["r","l"])
	call tH%get_mat(["u0","u1"],["d0","d1"],H)
	call allocate(HU,H)
	!call syev(HU,He,"V")
	call heev(HU,He,"V")


	label=["u","u'","r","r'","d","d'","l","l'"]
	do i=1,size(A)
		call A(i)%new([label(::2),"p"],[d,d,d,d,p],flag="r")
	enddo
	call simp_update(H,[0.2d0,0.001d0],A,1000)


    !tp:        
    !
	! 5---1---6
	! |   |   |
	! 4--0/9--2
	! |   |   |
	! 8---3---7
	!

	do i=1,size(A)
		aa(i,0)=dot(A(i),set_conjg(A(i)%change_label(label(::2),label(2::2))))
		do j=1,4
			aa(i,j)=contract(aa(i,0),pack(label,flag==j))
			aa(i,j)=aa(i,j)%merge_label(pack(label,.not.(flag==j)),pack(label(::2),.not.(flag(::2)==j)))
		enddo
		do j=1,4
			aa(i,j+4)=contract(aa(i,0),pack(label,flag==j.or.flag==mod(j-2+4,4)+1))
			aa(i,j+4)=aa(i,j+4)%merge_label(pack(label,.not.(flag==j.or.flag==mod(j-2+4,4)+1)),&
				pack(label(::2),.not.(flag(::2)==j.or.flag(::2)==mod(j-2+4,4)+1)))
		enddo
		aa(i,0)=aa(i,0)%merge_label(label,label(::2))
		aa(i,9)=aa(i,0)%clone()
	enddo

    !fig:        
    !         1         2         3
	!    1---------2---------3---------4
	!    |         |         |         |
	!  13|       14|       15|       16|
	!    |         |         |         |
	!    |    4    |    5    |    6    |
	!    5---------6---------7---------8
	!    |         |         |         |
	!  17|       18|       19|       20|
	!    |         |         |         |
	!    |    7    |    8    |    9    |
	!    9--------10--------11--------12
	!    |         |         |         |
	!  21|       22|       23|       24|
	!    |         |         |         |
	!    |    10   |    11   |    12   |
	!    13-------14--------15--------16
	!

	call gen_fig(fig,[ -1,1,2, -1,2,3, -1,3,4, -1,5,6, -1,6,7, -1,7,8, -1,9,10, -1,10,11, -1,11,12, -1,13,14, -1,14,15, -1,15,16,&
		-2,1,5, -2,2,6, -2,3,7, -2,4,8, -2,5,9, -2,6,10, -2,7,11, -2,8,12, -2,9,13, -2,10,14, -2,11,15, -2,12,16],&
			[1,5, 2,1, 1,1, 2,6, 2,4, 1,0, 2,0, 1,2, 1,4, 2,0, 1,0, 2,2, 2,8, 1,3, 2,3, 1,7])


    !fig:           
    !         1     
	!    1---------2
	!    |         |
	!   5|        6|
	!    |         |
	!    |    2    |
	!    3---------4
	!    |         |
	!   7|        8|
	!    |         |
	!    |    3    |
	!    5---------6
	!    |         |
	!   9|       10|
	!    |         |
	!    |    4    |
	!    7---------8
	!

	call gen_fig(figx,[ -1,1,2, -1,3,4, -1,5,6, -1,7,8,&
		-2,1,3, -2,2,4, -2,3,5, -2,4,6, -2,5,7, -2,6,8],&
			[1,5, 2,6, 2,4, 1,2, 1,4, 2,2, 2,8, 1,7])

    !fig:           
    !         1         2         3     
	!    1---------2---------3---------4
	!    |         |         |         |
	!   7|        8|        9|       10|
	!    |         |         |         |
	!    |    4    |    5    |    6    |
	!    5---------6---------7---------8
	!

	call gen_fig(figy,[ -1,1,2, -1,2,3, -1,3,4, -1,5,6, -1,6,7, -1,7,8,&
		-2,1,5, -2,2,6, -2,3,7, -2,4,8],&
			[1,5, 2,1, 1,1, 2,6, 2,8, 1,3, 2,3, 1,7])

	call rescale(fig,aa(:,0:8),A)

	do step=0,Nstep-1
		!delta=0.05d0*exp(-alpha*step/Nstep)
		!delta=0.1d0
		delta=max(0.1d0-(0.1d0-0.005d0)/(Nstep/4)*step,0.001d0)
		Hexp=matmul(matmul(HU,diag(exp(-He*delta))),transpose(conjg(HU)))
		do n=1,4

			!call ini_que(fig,[integer::],[integer::],aa,que)
			!call contract_order(que)
			!E=fig_contract(get_dat(que%tail))
			!write(*,*)E%rc%T(1)


			!call ini_que(fig,[integer::],[integer::],aa,que)
			!call contract_order(que)
			!E=fig_contract(get_dat(que%tail))
			!write(*,*)E%rc%T(1)

			!stop

			select case(n)
			case(1:2)
				bd=5
			case(3:4)
				bd=18
			end select

			lbd=[fig%st(fig%bd(bd,1))%label(fig%bd2st(bd,1)),fig%st(fig%bd(bd,2))%label(fig%bd2st(bd,2))]
			do i=1,2
				k=fig%st(fig%bd(bd,i))%isub

				if(i==1) then
					call X(i)%new([label(::2)],merge(d*p,d,label(::2)==lbd(1)))
					call b(i)%new([lbd,"p"],[d,d*p,p])
					call A(k)%qr([character::],[lbd(1),"p"],Q,R)
					call X(i)%get_tensor(A(k)%label(1:4),Q)
					call b(i)%get_tensor([lbd(2),A(k)%label(4:)],R)
				else
					call X(i)%new([label(::2)],merge(d*p,d,label(::2)==lbd(2)))
					call b(i)%new([lbd,"p"],[d*p,d,p])
					call A(k)%lq([lbd(2),"p"],[character::],Q,L)
					call X(i)%get_tensor([A(k)%label(1),A(k)%label(3:)],Q)
					call b(i)%get_tensor([A(k)%label(1:2),lbd(1)],L)
				endif

				aa(k,9)=dot(X(i),set_conjg(X(i)%change_label(label(::2),label(2::2))))
				aa(k,9)=aa(k,9)%merge_label(label,label(::2))
			enddo

			fig%st(fig%bd(bd,:))%tp=9


			call ini_que(fig,[bd],[integer::],aa,que)
			call contract_order(que)
			E=fig_contract(get_dat(que%tail))
			E=E%split_label(E%label,[character(2)::E%label(1)(1:1),E%label(1)(1:1)//"'",E%label(2)(1:1),E%label(2)(1:1)//"'"])
			!write(*,*)E%label,E%shap
			!stop


			call full_update(E,Hexp,lbd,b)

			do i=1,size(A)
				k=fig%st(fig%bd(bd,i))%isub
				if(i==1) then
					A(k)=dot(X(i),b(i),lbd)
				else
					A(k)=dot(b(i),X(i),lbd)
				endif
				aa(k,9)=merge_label(dot(A(k),set_conjg(A(k)%change_label(label(::2),label(2::2)))),label,label(::2))
			enddo


			!call allocate(ps,aa(1,5)%shap(1))
			!ps=0d0
			!k=0
			!do 
				!k=k+1
				if(n==1.or.n==2) then
					!call CTM(fig,[1,4,7,10],aa)
					!call CTM(fig,[3,6,9,12],aa)
					call CTM(fig,reshape([1,4,7,10, 3,6,9,12],[4,2]),aa)
					!fig%st%isub=mod(fig%st%isub,2)+1
					!call rescale(fig,aa(:,0:8),A)
					!call ini_que(figx,[integer::],[integer::],aa,que)
					!call contract_order(que)
					!E=fig_contract(get_dat(que%tail))
					!write(*,*)E%label,E%shap,E%rc%T(1),figx%st(1)%isub
					!figx%st%isub=mod(figx%st%isub,2)+1
				else
					!call CTM(fig,[13,14,15,16],aa)
					!call CTM(fig,[21,22,23,24],aa)
					call CTM(fig,reshape([13,14,15,16, 21,22,23,24],[4,2]),aa)
					!fig%st%isub=mod(fig%st%isub,2)+1
					!call rescale(fig,aa(:,0:8),A)
					!call ini_que(figy,[integer::],[integer::],aa,que)
					!call contract_order(que)
					!E=fig_contract(get_dat(que%tail))
					!write(*,*)E%label,E%shap,E%rc%T(1),figy%st(1)%isub
					!figy%st%isub=mod(figy%st%isub,2)+1
				endif
				!call aa(1,5)%svd(["d"],["r"],s=s)
				!if(sum(abs(ps/ps(1)-s/s(1)))<1d-4.or.k>40) then
					!exit
				!endif
				!call allocate(ps,s)
				!call aa(2,5)%svd(["d"],["r"],s=s)
				!write(*,*)s(1),s(size(s))
				!read(*,*)
			!enddo

			if(fig%st(fig%bd(bd,1))%tp==9) then
				fig%st(fig%bd(bd,:))%tp=0
				do i=1,size(A)
					aa(i,0)=aa(i,9)%clone()
				enddo
			endif
			fig%st%isub=mod(fig%st%isub,2)+1

			call rescale(fig,aa(:,0:8),A)
			!read(*,*)
		enddo
		write(*,*)"Energy: ",get_vbond(fig,A,aa,[5,8,18,19],H)/2d0
		!write(*,*)"Mz: ",get_vsite(fig,A,aa,[6,7,10,11],cmplx(reshape([1d0,0d0,0d0,-1d0],[2,2]),kind=8))/4d0
		!sx=get_vsite(fig,A,aa,[6],Hsx)**2
		!sy=get_vsite(fig,A,aa,[6],Hsy)**2
		!sz=get_vsite(fig,A,aa,[6],Hsz)**2
		!write(*,*)"M: ",0.5d0*sqrt(abs(sx+sy+sz))
	enddo
end program

