include "tensor_real.f90"
include "list.f90"
module M_fig
	use M_tensor
	use M_list
	use M_lattice_test1
	implicit none
	integer :: d=3,chi,p,Nstep=300
	integer, parameter :: gp=1,gd=2,gchi=3,gpd=4,gdd=22,gtmp=5
	real(8) :: t(1)=1d0,DJ=1d0,DV=2d0,DU=0d0,mu=0d0,norm=1d0
	type t_figst
		integer, allocatable :: bd(:)
		integer, allocatable :: ltp(:)
	end type
	type t_fig
		character(:), allocatable :: label(:)
		integer, allocatable :: bd(:,:)
		integer, allocatable :: btp(:)
		type(t_figst), allocatable :: st(:)
		integer, allocatable :: stp(:)
		integer, allocatable :: ttp(:)
	end type
	type t_dat
		type(t_tensor) :: T
		type(t_dat), pointer :: link1 => null()
		type(t_dat), pointer :: link2 => null()
		integer(8) :: cost(2)
		!integer(8) :: icost(200)
		integer :: cfg=0
		integer :: n
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
	function fig_contract(que,ord) result(rt)
		type(t_tensor) :: que(:)
		integer :: ord(2,*)
		type(t_tensor) :: tmp(size(que)-2)
		type(t_tensor) :: rt
		integer :: i,j
		integer(8) :: cost(2),m,sz(3)
		cost=0
		m=0
		do i=1,size(tmp)
			if(ord(1,i)<0.and.ord(2,i)<0) then
				sz(1:2)=[size(que(abs(ord(1,i)))%rc%T,kind=8),size(que(abs(ord(2,i)))%rc%T,kind=8)]
				tmp(i)=dot(que(abs(ord(1,i))),que(abs(ord(2,i))))
				tmp(i)%is_return=.true.
				sz(3)=0
			elseif(ord(1,i)<0) then
				sz(1:2)=[size(que(abs(ord(1,i)))%rc%T,kind=8),size(tmp(ord(2,i))%rc%T,kind=8)]
				tmp(i)=dot(que(abs(ord(1,i))),tmp(ord(2,i)))
				tmp(i)%is_return=.true.
				sz(3)=sz(2)
			else
				sz(1:2)=[size(tmp(ord(1,i))%rc%T,kind=8),size(tmp(ord(2,i))%rc%T,kind=8)]
				tmp(i)=dot(tmp(ord(1,i)),tmp(ord(2,i)))
				tmp(i)%is_return=.true.
				sz(3)=sum(sz(1:2))
			endif
			cost(1)=cost(1)+size(tmp(i)%rc%T,kind=8)*nint(sqrt(real(product(sz(1:2))/size(tmp(i)%rc%T,kind=8))))
			cost(2)=max(cost(2),m+size(tmp(i)%rc%T,kind=8))
			m=m-sz(3)+size(tmp(i)%rc%T,kind=8)
		enddo
		if(ord(1,i)<0.and.ord(2,i)<0) then
			sz(1:2)=[size(que(abs(ord(1,i)))%rc%T,kind=8),size(que(abs(ord(2,i)))%rc%T,kind=8)]
			rt=dot(que(abs(ord(1,i))),que(abs(ord(2,i))))
			rt%is_return=.true.
			sz(3)=0
		elseif(ord(1,i)<0) then
			sz(1:2)=[size(que(abs(ord(1,i)))%rc%T,kind=8),size(tmp(ord(2,i))%rc%T,kind=8)]
			rt=dot(que(abs(ord(1,i))),tmp(ord(2,i)))
			rt%is_return=.true.
			sz(3)=sz(2)
		else
			sz(1:2)=[size(tmp(ord(1,i))%rc%T,kind=8),size(tmp(ord(2,i))%rc%T,kind=8)]
			rt=dot(tmp(ord(1,i)),tmp(ord(2,i)))
			rt%is_return=.true.
			sz(3)=sum(sz(1:2))
		endif
		cost(1)=cost(1)+size(rt%rc%T,kind=8)*nint(sqrt(real(product(sz(1:2))/size(rt%rc%T,kind=8))))
		cost(2)=max(cost(2),m+size(rt%rc%T,kind=8))
		!write(*,*)cost
	end function
	function merge_node(dat1,dat2,dat) result(flag)
		type(t_dat), pointer :: dat1,dat2
		type(t_dat), pointer :: dat
		logical :: flag
		logical, allocatable :: flag1(:),flag2(:)
		integer :: i,j
		integer(8), allocatable :: shap(:),shap1(:),shap2(:)
		integer(8) :: cost
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

			!dat%cost(2)=max(dat%cost(2),dat%cost(1)+dat1%cost(2)+dat2%cost(2))

			dat%cost(2)=dat%cost(1)
			allocate(shap1(size(dat1%T%shap)))
			if(dat1%cost(2)/=0) then 
				shap1=dat1%T%shap
			else
				shap1=0
			endif
			dat%cost(2)=dat%cost(2)+product(shap1)
			allocate(shap2(size(dat2%T%shap)))
			if(dat2%cost(2)/=0) then 
				shap2=dat2%T%shap
			else
				shap2=0
			endif
			dat%cost(2)=dat%cost(2)+product(shap2)
			dat%cost(2)=max(max(dat%cost(2),dat1%cost(2)+product(shap2)),dat2%cost(2)+product(shap1))

			dat%cost(1)=dat%cost(1)*cost+dat1%cost(1)+dat2%cost(1)
			dat%link1 => dat1
			dat%link2 => dat2
			dat%flag=.false.
		else
			flag=.false.
		endif
	end function
	subroutine get_corder(que,ord)
		type(t_tensor) :: que(:)
		integer, allocatable :: ord(:)
		type(t_list) :: S(size(que))
		integer :: n,i,j,c,d,k,minx
		integer(8) :: u(2),uc,u0
		type(t_dat), pointer :: dat,dat1,dat2,dat3
		type(t_node), pointer :: pos1,pos2,pos3
		integer :: cm
		!write(*,*)"corder start"
		if(allocated(ord)) deallocate(ord)
		allocate(ord(2*(size(que)-1)))
		cm=1
		n=size(que)
		minx=10000
		do i=1,n
			allocate(dat)
			dat%T=que(i)
			dat%cfg=ibset(0,i-1)
			dat%n=-i
			if(minx>minval(dat%T%shap)) minx=minval(dat%T%shap)
			!dat%cost=[0,product(dat%T%shap)]
			dat%cost=0
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
			u(2)=100000000000000
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
								elseif(dat%cost(cm)<=uc.and.dat%cost(cm)>u0.and.dat%cost(2)<2000000000) then
									!write(*,"(A$)")dat1%T%label,"|",dat2%T%label,"->",dat%T%label
									!write(*,*)"cost:",dat%cost
									pos3 => S(c)%head
									do k=1,S(c)%num
										if(k/=1) pos3 => pos3%next
										dat3 => get_dat(pos3)
										if(dat%cfg==dat3%cfg) then
											!if((dat%cost(cm)<=dat3%cost(cm).and.dat%cost(mod(cm,2)+1)<dat3%cost(mod(cm,2)+1)).or.(c==n.and.dat3%cost(mod(cm,2)+1)/dat%cost(mod(cm,2)+1)>10.and.dat3%cost(mod(cm,2)+1)>40000000)) then
											if((dat%cost(cm)<dat3%cost(cm)).or.(dat%cost(cm)==dat3%cost(cm).and.dat%cost(mod(cm,2)+1)<dat3%cost(mod(cm,2)+1))) then
												dat3%flag=.true.
												dat3%cost=dat%cost
												dat3%link1 => dat%link1
												dat3%link2 => dat%link2
											endif
											exit
										endif
									enddo
									if(k==S(c)%num+1) then
										!write(*,*)c,k,dat%cost
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
			!write(*,*)S%num
			!write(*,*)uc
			!read(*,*)
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

		k=0
		dat => get_dat(S(n)%head)
		!write(*,*)"the cost is: ",dat%cost
		call gen_ord(ord,k,dat)

	end subroutine
	recursive subroutine gen_ord(ord,n,dat)
		integer :: ord(2,*)
		integer :: n
		type(t_dat), pointer :: dat
		if(dat%flag) return
		if(.not.dat%link1%flag) then
			call gen_ord(ord,n,dat%link1)
		endif
		if(.not.dat%link2%flag) then
			call gen_ord(ord,n,dat%link2)
		endif
		n=n+1
		ord(1,n)=dat%link1%n
		ord(2,n)=dat%link2%n
		if(ord(1,n)>0.and.ord(2,n)<0) then
			ord([1,2],n)=ord([2,1],n)
		endif
		dat%n=n
		dat%flag=.true.
	end subroutine
	subroutine gen_fig(fig,bd)
		type(t_fig) :: fig
		integer :: bd(:)
		integer, allocatable :: num(:)
		integer :: i,n,k
		fig%label=["r","l","d","u"]
		allocate(fig%bd(2,size(bd)/3),fig%btp(size(bd)/3))
		fig%bd(1,:)=bd(2::3)
		fig%bd(2,:)=bd(3::3)
		fig%btp=bd(1::3)
		allocate(fig%st(1:maxval(fig%bd)))
		allocate(fig%stp(size(fig%st)),fig%ttp(size(fig%st)),num(size(fig%st)))
		num=0
		do n=1,size(fig%bd,2)
			num(fig%bd(:,n))=num(fig%bd(:,n))+1
		enddo
		do i=1,size(fig%st)
			allocate(fig%st(i)%bd(num(i)),fig%st(i)%ltp(num(i)))
		enddo
		num=0
		do n=1,size(fig%bd,2)
			num(fig%bd(:,n))=num(fig%bd(:,n))+1
			do k=1,2
				fig%st(fig%bd(k,n))%bd(num(fig%bd(k,n)))=n
				fig%st(fig%bd(k,n))%ltp(num(fig%bd(k,n)))=2*(fig%btp(n)-1)+k
				do i=1,num(fig%bd(k,n))-1
					if(fig%st(fig%bd(k,n))%ltp(i)>fig%st(fig%bd(k,n))%ltp(num(fig%bd(k,n)))) then
						fig%st(fig%bd(k,n))%ltp([i,num(fig%bd(k,n))])=fig%st(fig%bd(k,n))%ltp([num(fig%bd(k,n)),i])
						fig%st(fig%bd(k,n))%bd([i,num(fig%bd(k,n))])=fig%st(fig%bd(k,n))%bd([num(fig%bd(k,n)),i])
						exit
					endif
				enddo
			enddo
		enddo
		do i=1,size(fig%st)
			select case(size(fig%st(i)%ltp))
			case(2)
				if(all(fig%st(i)%ltp-[1,3]==0)) then
					fig%ttp(i)=5
				elseif(all(fig%st(i)%ltp-[2,3]==0)) then
					fig%ttp(i)=6
				elseif(all(fig%st(i)%ltp-[2,4]==0)) then
					fig%ttp(i)=7
				else
					fig%ttp(i)=8
				endif
			case(3)
				if(all(fig%st(i)%ltp-[1,2,3]==0)) then
					fig%ttp(i)=1
				elseif(all(fig%st(i)%ltp-[2,3,4]==0)) then
					fig%ttp(i)=2
				elseif(all(fig%st(i)%ltp-[1,2,4]==0)) then
					fig%ttp(i)=3
				else
					fig%ttp(i)=4
				endif
			case(4)
				fig%ttp(i)=0
			end select
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
	function cut_label(fig,i,cut) result(rt)
		type(t_fig) :: fig
		integer :: i
		integer :: cut(:)
		integer :: n
		character(:), allocatable :: rt(:)
		rt=[character::]
		do n=1,size(fig%st(i)%bd)
			if(any(fig%st(i)%bd(n)==cut)) then
				rt=[rt,trim(fig%label(fig%st(i)%ltp(n)))//to_char(fig%st(i)%bd(n))]
			else
				rt=[rt,to_char(fig%st(i)%bd(n))]
			endif
		enddo
	end function
	subroutine ini_que(fig,cut,rm,aa,que)
		type(t_fig) :: fig
		integer :: cut(:),rm(:)
		type(t_tensor) :: aa(:,-1:)
		type(t_tensor), allocatable :: que(:)
		!type(t_dat), pointer :: dat
		integer :: i,n,pos,k
		character(:), allocatable :: to_label(:),label(:)
		do k=1,2
			pos=0
			do i=1,size(fig%st)
				if((any(i==abs(rm)).and.any(rm>0)).or.(all(i/=abs(rm)).and.any(rm<0))) cycle
				label=[character::]
				to_label=[character::]
				select case(fig%ttp(i))
				case(-1)
					pos=pos+2
					if(k==2) then
						label=fig%label(fig%st(i)%ltp)
						to_label=cut_label(fig,i,cut)
						que(pos-1)=aa(fig%stp(i),fig%ttp(i))%change_label([label,"p"],[to_label,"p"//to_char(i)])
						que(pos)=aa(fig%stp(i),fig%ttp(i))%change_label([label,"p"],[(trim(to_label(n))//"'",n=1,size(to_label)),"p"//to_char(i)])


					endif
				!case(1:4)
					!pos=pos+1
					!if(k==2) then
						!do n=1,size(fig%st(i)%ltp)
							!label=[label,fig%label(fig%st(i)%ltp(n))]
							!if(any(fig%st(i)%bd(n)==cut)) then
								!to_label=[to_label,trim(fig%label(fig%st(i)%ltp(n)))//to_char(fig%st(i)%bd(n))]
							!else
								!to_label=[to_label,to_char(fig%st(i)%bd(n))]
							!endif
							!if(size(fig%st(merge(fig%bd(1,fig%st(i)%bd(n)),fig%bd(2,fig%st(i)%bd(n)),fig%bd(2,fig%st(i)%bd(n))==i))%bd)==4) then
								!label=[label,trim(label(size(label)))//"'"]
								!to_label=[to_label,trim(to_label(size(to_label)))//"'"]
							!endif
						!enddo
						!que(pos)=aa(fig%stp(i),fig%ttp(i))%change_label(label,to_label)
					!endif
				case default
					pos=pos+1
					if(k==2) then
						!write(*,*)i,fig%stp(i),fig%ttp(i),aa(fig%stp(i),fig%ttp(i))%label,fig%label(fig%st(i)%ltp),cut_label(fig,i,cut)
						que(pos)=aa(fig%stp(i),fig%ttp(i))%change_label(fig%label(fig%st(i)%ltp),cut_label(fig,i,cut))
					endif
				end select
			enddo
			if(k==2) return
			if(allocated(que)) then
				deallocate(que)
			endif
			allocate(que(pos))
		enddo
	end subroutine
	subroutine CTM(fig,pair,aa)
		type(t_fig) :: fig
		integer :: pair(:,:)
		type(t_tensor) :: aa(:,-1:)
		type(t_tensor) :: P(2),E
		integer :: i,n1,n2,ist(2,size(pair)+size(pair,2)),n
		character(:), allocatable :: lbd(:)
		real(8), allocatable :: U(:,:),V(:,:)
		real(8), allocatable :: s(:)
		type(t_tensor), allocatable :: que(:)
		type(t_tensor) :: aa_(size(ist,2))
		integer, allocatable :: cord(:)
		select case(fig%btp(pair(1,1)))
		case(2)
			lbd=["d","u","r","l"]
		case(1)
			lbd=["r","l","d","u"]
		end select
		n=0
		do n2=1,size(pair,2)
			if(.not.(any(fig%label(fig%st(fig%bd(1,pair(1,n2)))%ltp)==lbd(3)))) then
				lbd([3,4])=lbd([4,3])
			endif
			do n1=1,size(pair,1),2
				call ini_que(fig,pair(n1:n1+1,n2),[integer::],aa,que)
				call get_corder(que,cord)
				E=fig_contract(que,cord)
				deallocate(que)
				call E%svd([lbd(1)//to_char(pair(n1,n2)),lbd(1)//to_char(pair(n1+1,n2))],[lbd(2)//to_char(pair(n1,n2)),lbd(2)//to_char(pair(n1+1,n2))],U,s,V,grp(:,:,gchi))
				call P(1)%new(lbd(1)//["1","2",""],[E%shap(1:2),chi],[E%stp(1:2),gchi])
				call P(1)%get_tensor([P(1)%label],U(:,:chi))
				P(1)=P(1)
				P(2)=P(1)%change_label(lbd(1)//["1","2",""],lbd(2)//["1","2",""])
				P%is_return=.true.
				do i=1,2
					if(fig%bd(i,pair(n1,n2))/=ist(1,n)) then
						n=n+1
						ist(:,n)=fig%bd(i,pair(n1:n1+1,n2))
					endif
					if(allocated(aa_(n)%shap)) then
						aa_(n)=dot(aa_(n),P(i))
					else
						aa_(n)=dot(dot(&
							aa(fig%stp(ist(1,n)),fig%ttp(ist(1,n)))%change_label(pack([lbd(i),lbd(3),lbd(mod(i,2)+1)],[1:3]<=size(fig%st(ist(1,n))%bd)),pack([lbd(i)//"1","m ",lbd(mod(i,2)+1)//"1"],[1:3]<=size(fig%st(ist(1,n))%bd))),P(i)),&
							aa(fig%stp(ist(2,n)),fig%ttp(ist(2,n)))%change_label(pack([lbd(i),lbd(4),lbd(mod(i,2)+1)],[1:3]<=size(fig%st(ist(1,n))%bd)),pack([lbd(i)//"2","m ",lbd(mod(i,2)+1)//"2"],[1:3]<=size(fig%st(ist(1,n))%bd))))
					endif
				enddo
			enddo
		enddo

		do i=1,n
			aa(fig%stp(ist(2,i)),fig%ttp(ist(1,i)))=aa_(i)
		enddo

	end subroutine
	subroutine dCTM(fig,label,aa)
		type(t_fig) :: fig
		character(*) :: label(:)
		type(t_tensor) :: aa(:,-1:)
		type(t_tensor) :: C(2),E(2),a(4)
		integer :: i,j,n,ist(8),l(8),st
		character(:), allocatable :: bd(:),bdp(:),lbd(:)
		real(8), allocatable :: U1(:,:),U2(:,:),Z(:,:),W(:,:),cZ(:,:),cW(:,:)
		real(8), allocatable :: Eg(:)
		type(t_tensor), allocatable :: que(:,:)
		integer, allocatable :: cord(:)
		write(*,*)"dCTM"
		l=[1,2,5,6,5,6,1,2]
		do n=1,size(label)
			allocate(que(7,2))
			select case(label(n))
			case("r")
				ist(:)=[1,2,5,6,9,10,13,14]
				bd=["13","14","17","18","21","22"]
				bdp=["2","5","8","11"]
				lbd=["d","u"]
			case("l")
				ist(:)=[4,3,8,7,12,11,16,15]
				bd=["16","15","20","19","24","23"]
				bdp=["2","5","8","11"]
				lbd=["d","u"]
			case("d")
				ist(:)=[1,5,2,6,3,7,4,8]
				bd=["1","4","2","5","3","6"]
				bdp=["17","18","19","20"]
				lbd=["r","l"]
			case("u")
				ist(:)=[13,9,14,10,15,11,16,12]
				bd=["10","7","11","8","12","9"]
				bdp=["17","18","19","20"]
				lbd=["r","l"]
			end select
			do i=1,2
				do j=1,4
					st=ist(j+(i-1)*4)
					que(l(j+(i-1)*4),i)=aa(fig%stp(st),fig%ttp(st))%change_label(fig%label(fig%st(st)%ltp),to_char(fig%st(st)%bd))
				enddo
				call get_corder(que(1:2,i),cord)
				C(i)=fig_contract(que(1:2,i),cord)
				que(4,i)=C(i)
				call get_corder(que(4:6,i),cord)
				E(i)=fig_contract(que(4:6,i),cord)
			enddo
			
			call get_mat(dot(C(1),C(1)%change_label(bd(1:2),bd(1:2)//"'")),bd(1:2),bd(1:2)//"'",U1)
			call get_mat(dot(C(2),C(2)%change_label(bd(5:6),bd(5:6)//"'")),bd(5:6),bd(5:6)//"'",U2)
			!call get_mat(dot(C(2),C(2)%change_label(bd(5:6),bd(5:6)//"'")),bd(5:6),bd(5:6)//"'",U2)
			!U2=matmul(matmul(diag(fgate(gstp(C(2)%stp(C(2)%get_ldx(bd(5:6)))))),U2),diag(fgate(gstp(C(2)%stp(C(2)%get_ldx(bd(5:6)))))))
			call allocate(Z,U1+U2)
			call que(3,1)%new([bd(1:2),lbd(1)],[C(1)%shap(C(1)%get_ldx(bd(1:2))),chi],[C(1)%stp(C(1)%get_ldx(bd(1:2))),gchi])
			call smat_eg(que(3,1)%stp(1:2),Z,Eg,cZ,grp(:,:,gchi))
			call que(3,1)%get_tensor(que(3,1)%label,cZ(:,:chi))
			que(3,1)=que(3,1)
			que(3,2)=que(3,1)%change_label([bd(1:2),lbd(1)],[bd(5:6),lbd(2)])
			que(4,1)=que(3,1)%change_label([lbd(1)],[lbd(2)])
			que(4,2)=que(3,1)%change_label(bd(1:2),bd(5:6))

			call get_mat(dot(E(1),E(1)%change_label(bd(3:4),bd(3:4)//"'")),bd(3:4),bd(3:4)//"'",U1)
			call get_mat(dot(E(2),E(2)%change_label(bd(3:4),bd(3:4)//"'")),bd(3:4),bd(3:4)//"'",U2)
			!call get_mat(dot(E(2),E(2)%change_label(bd(3:4),bd(3:4)//"'")),bd(3:4),bd(3:4)//"'",U2)
			!U2=matmul(matmul(diag(fgate(gstp(E(2)%stp(E(2)%get_ldx(bd(3:4)))))),U2),diag(fgate(gstp(E(2)%stp(E(2)%get_ldx(bd(3:4)))))))
			call allocate(W,U1+U2)
			call que(7,1)%new([bd(3:4),lbd(1)],[E(1)%shap(E(1)%get_ldx(bd(3:4))),chi],[E(1)%stp(E(1)%get_ldx(bd(3:4))),gchi])
			call smat_eg(que(7,1)%stp(1:2),W,Eg,cW,grp(:,:,gchi))
			call que(7,1)%get_tensor(que(7,1)%label,cW(:,:chi))
			que(7,2)=que(7,1)%change_label([lbd(1)],[lbd(2)])
			que(7,1)=que(7,1)

			call get_corder(que(1:3,1),cord)
			a(1)=change_label(fig_contract(que(1:3,1),cord),[bdp(1)],[label(n)])
			call get_corder(que(1:3,2),cord)
			a(4)=change_label(fig_contract(que(1:3,2),cord),[bdp(4)],[label(n)])
			call get_corder(que(4:7,1),cord)
			a(2)=change_label(fig_contract(que(4:7,1),cord),[bdp(2)],[label(n)])
			call get_corder(que(4:7,2),cord)
			a(3)=change_label(fig_contract(que(4:7,2),cord),[bdp(3)],[label(n)])
			deallocate(que)
			aa(fig%stp(ist(2)),fig%ttp(ist(1)))=a(1)
			aa(fig%stp(ist(4)),fig%ttp(ist(3)))=a(2)
			aa(fig%stp(ist(6)),fig%ttp(ist(5)))=a(3)
			aa(fig%stp(ist(8)),fig%ttp(ist(7)))=a(4)
		enddo

	end subroutine
	subroutine dCTMp(fig,label,aa)
		type(t_fig) :: fig
		character(*) :: label(:)
		type(t_tensor) :: aa(:,-1:)
		type(t_tensor) :: E(2),a(4),P(6)
		integer :: i,j,n,k,l,ist(4,2),istp(4,2),mst(4,4),Pst(2,3),stP(2,4)
		character(:), allocatable :: bd(:,:),lbd(:)
		real(8), allocatable :: U1(:,:),U2(:,:),W(:,:),cW(:,:)
		real(8), allocatable :: Eg(:)
		type(t_tensor), allocatable :: que(:)
		integer, allocatable :: cord(:)
		Pst=reshape([1,2,2,3,3,4],[2,3])
		stP=reshape([1,0,2,3,4,5,6,0],[2,4])
		do n=1,size(label)
			allocate(que(4))
			select case(label(n))
			case("r")
				ist(:,1)=[5,6,1,2]
				ist(:,2)=[9,10,13,14]
				mst=reshape([1,2,4,3,5,6,8,7,9,10,12,11,13,14,16,15],[4,4])
				bd=reshape(["17","18","13","14", "17","18","21","22"],[4,2])
				lbd=["r","l","d","u"]
			case("l")
				ist(:,1)=[8,7,4,3]
				ist(:,2)=[12,11,16,15]
				mst=reshape([4,3,1,2,8,7,5,6,12,11,9,10,16,15,13,14],[4,4])
				bd=reshape(["20","19","16","15", "20","19","24","23"],[4,2])
				lbd=["l","r","d","u"]
			case("d")
				ist(:,1)=[2,6,1,5]
				ist(:,2)=[3,7,4,8]
				mst=reshape([1,5,13,9,2,6,14,10,3,7,15,11,4,8,16,12],[4,4])
				bd=reshape(["2","5","1","4", "2","5","3","6"],[4,2])
				lbd=["d","u","r","l"]
			case("u")
				ist(:,1)=[14,10,13,9]
				ist(:,2)=[15,11,16,12]
				mst=reshape([13,9,1,5,14,10,2,6,15,11,3,7,16,12,4,8],[4,4])
				bd=reshape(["11","8","10","7", "11","8","12","9"],[4,2])
				lbd=["u","d","r","l"]
			end select
			do i=1,3
				if((.not.allocated(P((i-1)*2+1)%shap))) then
					if(i/=2) then
						l=3
					else
						l=1
					endif
					do k=1,2
						istp(l:,k)=mod(fig%stp(ist(l:,k))+1+fig%stp(mst(1,Pst(k,i)))-fig%stp(ist(l,k)),2)+1
						do j=l,4
							que(j)=aa(istp(j,k),fig%ttp(ist(j,k)))%change_label(fig%label(fig%st(ist(j,k))%ltp),to_char(fig%st(ist(j,k))%bd))
						enddo
						call get_corder(que(l:4),cord)
						E(k)=change_label(fig_contract(que(l:4),cord),bd(l:l+1,k),["c1","c2"])
					enddo

					call get_mat(dot(E(1),E(1)%change_label(["c1","c2"],["c1'","c2'"])),["c1","c2"],["c1'","c2'"],U1)
					call get_mat(dot(E(2),E(2)%change_label(["c1","c2"],["c1'","c2'"])),["c1","c2"],["c1'","c2'"],U2)
					!call get_mat(dot(E(2),E(2)%change_label(["c1","c2"],["c1'","c2'"])),["c1","c2"],["c1'","c2'"],U2)
					call allocate(W,U1+U2)
					call smat_eg(E(1)%stp(E(1)%get_ldx(["c1","c2"])),W,Eg,cW,grp(:,:,gchi))

					call P((i-1)*2+1)%new(lbd(3)//["1","2"," "],[E(1)%shap(E(1)%get_ldx(["c1","c2"])),chi],[E(1)%stp(E(1)%get_ldx(["c1","c2"])),gchi])
					call P((i-1)*2+1)%get_tensor(lbd(3)//["1","2"," "],cW(:,:chi))

					P(i*2)=P((i-1)*2+1)%change_label(lbd(3)//["1","2"," "],lbd(4)//["1","2"," "])
					P((i-1)*2+1)=P((i-1)*2+1)
					do j=i+1,3
						if(fig%stp(mst(1,Pst(1,j)))==fig%stp(mst(1,Pst(1,i)))) then
							if((Pst(2,j)==0).eqv.(Pst(2,i)==0)) then
								P((j-1)*2+1)=P((i-1)*2+1)
								P(j*2)=P(i*2)
							endif
						endif
					enddo
				endif
			enddo
			
			P%is_return=.true.

			do i=1,4
				if(stP(2,i)/=0) then
					j=1
					que(4)=P(stP(2,i))
				else
					j=0
				endif
				que(1)=aa(fig%stp(mst(1,i)),fig%ttp(mst(1,i)))%change_label([lbd(3+i/4:3+i/4+j),lbd(1)],[lbd(3+i/4:3+i/4+j)//"1","m "])
				que(2)=aa(fig%stp(mst(2,i)),fig%ttp(mst(2,i)))%change_label([lbd(3+i/4:3+i/4+j),lbd(2)],[lbd(3+i/4:3+i/4+j)//"2","m "])
				que(3)=P(stP(1,i))
				call get_corder(que(1:3+j),cord)
				a(i)=fig_contract(que(1:3+j),cord)
				a(i)%is_return=.true.
			enddo
			do i=1,4
				aa(fig%stp(mst(2,i)),fig%ttp(mst(1,i)))=a(i)
			enddo
			deallocate(que)
		enddo

	end subroutine
	subroutine vCTM(fig,label,aa)
		type(t_fig) :: fig
		character(*) :: label(:)
		type(t_tensor) :: aa(:,-1:)
		type(t_tensor) :: R,E(2,2),Q,a(4),P(6,size(label))
		integer :: i,j,n,k,ist(8,2),istp(8,2),mst(4,4),st,Pst(2,3),stP(2,4)
		character(:), allocatable :: bd(:),lbd(:)
		real(8), allocatable :: U(:,:),V(:,:)
		real(8), allocatable :: s(:)
		type(t_tensor), allocatable :: que(:)
		integer, allocatable :: cord(:),mgrp(:,:)
		Pst=reshape([1,2,2,3,3,4],[2,3])
		stP=reshape([1,0,2,3,4,5,6,0],[2,4])
		do n=1,size(label)
			allocate(que(8))
			select case(label(n))
			case("r")
				ist(:,1)=[5,6,8,7,1,2,4,3]
				ist(:,2)=[9,10,12,11,13,14,16,15]
				mst=reshape([1,2,4,3,5,6,8,7,9,10,12,11,13,14,16,15],[4,4])
				bd=["17","18","20","19"]
				lbd=["r","l","d","u"]
			case("l")
				ist(:,1)=[8,7,5,6,4,3,1,2]
				ist(:,2)=[12,11,9,10,16,15,13,14]
				mst=reshape([4,3,1,2,8,7,5,6,12,11,9,10,16,15,13,14],[4,4])
				bd=["20","19","17","18"]
				lbd=["l","r","d","u"]
			case("d")
				ist(:,1)=[2,6,14,10,1,5,13,9]
				ist(:,2)=[3,7,15,11,4,8,16,12]
				mst=reshape([1,5,13,9,2,6,14,10,3,7,15,11,4,8,16,12],[4,4])
				bd=["2","5","11","8"]
				lbd=["d","u","r","l"]
			case("u")
				ist(:,1)=[14,10,2,6,13,9,1,5]
				ist(:,2)=[15,11,3,7,16,12,4,8]
				mst=reshape([13,9,1,5,14,10,2,6,15,11,3,7,16,12,4,8],[4,4])
				bd=["11","8","2","5"]
				lbd=["u","d","r","l"]
			end select
			do i=1,3
				if((.not.allocated(P((i-1)*2+1,n)%shap))) then
					do k=1,2
						istp(:,k)=mod(fig%stp(ist(:,k))+1+fig%stp(mst(1,Pst(k,i)))-fig%stp(ist(1,k)),2)+1
						do j=1,8
							que(j)=aa(istp(j,k),fig%ttp(ist(j,k)))%change_label(fig%label(fig%st(ist(j,k))%ltp),to_char(fig%st(ist(j,k))%bd))
						enddo
						call get_corder(que(1:8),cord)
						R=change_label(fig_contract(que(1:8),cord),bd,["c1","c2","c3","c4"])
						call R%fact(["c1","c2"],[character::],E(k,1),Q,"rq",[lbd(5-k),lbd(5-k)],gtmp)
						if(size(label)>1) then
							call R%fact([character::],["c3","c4"],Q,E(k,2),"qr",[lbd(5-k),lbd(5-k)],gtmp)
							E(k,2)=E(k,2)%change_label(["c3","c4"],["c1","c2"])
						endif
					enddo
					do k=1,size(label)
						!call svd_t(dot(E(1,k),E(2,k)),[lbd(4)],[lbd(3)],U,s,V,grp(:,:,gchi),"s")
						call svd_t(dot(E(1,k),E(2,k)),[lbd(4)],[lbd(3)],U,s,V,grp(:,:,gchi))
						do j=1,chi
							s(j)=1d0/sqrt(s(j))
							U(:,j)=s(j)*U(:,j)
							V(j,:)=s(j)*V(j,:)
						enddo
						P((i-1)*2+1,k)=dot(V(:chi,:),E(2,k)%change_label(["c1","c2"],lbd(3)//["1","2"]),[lbd(3)],gchi)
						P(i*2,k)=dot(E(1,k)%change_label(["c1","c2"],lbd(4)//["1","2"]),U(:,:chi),[lbd(4)],gchi)
						do j=i+1,3
							if(fig%stp(mst(1+(k-1)*2,Pst(1,j)))==fig%stp(mst(1+(k-1)*2,Pst(1,i)))) then
								P((j-1)*2+1,k)=P((i-1)*2+1,k)
								P(j*2,k)=P(i*2,k)
							endif
						enddo
					enddo
				endif
			enddo

			P(:,n)%is_return=.true.

			do i=1,4
				if(stP(2,i)/=0) then
					j=1
					que(4)=P(stP(2,i),n)
				else
					j=0
				endif
				que(1)=aa(fig%stp(mst(1,i)),fig%ttp(mst(1,i)))%change_label([lbd(3+i/4:3+i/4+j),lbd(1)],[lbd(3+i/4:3+i/4+j)//"1","m "])
				que(2)=aa(fig%stp(mst(2,i)),fig%ttp(mst(2,i)))%change_label([lbd(3+i/4:3+i/4+j),lbd(2)],[lbd(3+i/4:3+i/4+j)//"2","m "])
				que(3)=P(stP(1,i),n)
				call get_corder(que(1:3+j),cord)
				a(i)=fig_contract(que(1:3+j),cord)
				a(i)%is_return=.true.
			enddo
			do i=1,4
				aa(fig%stp(mst(2,i)),fig%ttp(mst(1,i)))=a(i)
			enddo
			deallocate(que)
		enddo
		grp(:,:,gtmp:gtmp+1)=0

	end subroutine
	subroutine vCTMp(fig,label,aa)
		type(t_fig) :: fig
		character(*) :: label(:)
		type(t_tensor) :: aa(:,-1:)
		type(t_tensor) :: R,E(2),Q,a(4),P(6)
		integer :: i,j,n,k,ist(4,2),istp(4,2),mst(4,4),st,Pst(2,3),stP(2,4)
		character(:), allocatable :: bd(:),lbd(:)
		real(8), allocatable :: U(:,:),V(:,:)
		real(8), allocatable :: s(:)
		type(t_tensor), allocatable :: que(:)
		integer, allocatable :: cord(:),mgrp(:,:)
		Pst=reshape([1,2,2,3,3,4],[2,3])
		stP=reshape([1,0,2,3,4,5,6,0],[2,4])
		do n=1,size(label)
			allocate(que(4))
			select case(label(n))
			case("r")
				ist(:,1)=[5,6,1,2]
				ist(:,2)=[9,10,13,14]
				mst=reshape([1,2,4,3,5,6,8,7,9,10,12,11,13,14,16,15],[4,4])
				bd=["17","18","20","19"]
				lbd=["r","l","d","u"]
			case("l")
				ist(:,1)=[8,7,4,3]
				ist(:,2)=[12,11,16,15]
				mst=reshape([4,3,1,2,8,7,5,6,12,11,9,10,16,15,13,14],[4,4])
				bd=["20","19","17","18"]
				lbd=["l","r","d","u"]
			case("d")
				ist(:,1)=[2,6,1,5]
				ist(:,2)=[3,7,4,8]
				mst=reshape([1,5,13,9,2,6,14,10,3,7,15,11,4,8,16,12],[4,4])
				bd=["2","5","11","8"]
				lbd=["d","u","r","l"]
			case("u")
				ist(:,1)=[14,10,13,9]
				ist(:,2)=[15,11,16,12]
				mst=reshape([13,9,1,5,14,10,2,6,15,11,3,7,16,12,4,8],[4,4])
				bd=["11","8","2","5"]
				lbd=["u","d","r","l"]
			end select
			do i=1,3
				if((.not.allocated(P((i-1)*2+1)%shap))) then
					do k=1,2
						istp(:,k)=mod(fig%stp(ist(:,k))+1+fig%stp(mst(1,Pst(k,i)))-fig%stp(ist(1,k)),2)+1
						do j=1,4
							que(j)=aa(istp(j,k),fig%ttp(ist(j,k)))%change_label(fig%label(fig%st(ist(j,k))%ltp),to_char(fig%st(ist(j,k))%bd))
						enddo
						call get_corder(que(1:4),cord)
						R=change_label(fig_contract(que(1:4),cord),bd(1:2),["c1","c2"])
						call R%fact(["c1","c2"],[character::],E(k),Q,"rq",[lbd(5-k),lbd(5-k)],gtmp)
					enddo
					call svd_t(dot(E(1),E(2)),[lbd(4)],[lbd(3)],U,s,V,grp(:,:,gchi))
					do j=1,chi
						s(j)=1d0/sqrt(s(j))
						U(:,j)=s(j)*U(:,j)
						V(j,:)=s(j)*V(j,:)
					enddo
					P((i-1)*2+1)=dot(V(:chi,:),E(2)%change_label(["c1","c2"],lbd(3)//["1","2"]),[lbd(3)],gchi)
					P(i*2)=dot(E(1)%change_label(["c1","c2"],lbd(4)//["1","2"]),U(:,:chi),[lbd(4)],gchi)
					do j=i+1,3
						if(fig%stp(mst(1+(k-1)*2,Pst(1,j)))==fig%stp(mst(1+(k-1)*2,Pst(1,i)))) then
							P((j-1)*2+1)=P((i-1)*2+1)
							P(j*2)=P(i*2)
						endif
					enddo
				endif
			enddo

			P%is_return=.true.

			do i=1,4
				if(stP(2,i)/=0) then
					j=1
					que(4)=P(stP(2,i))
				else
					j=0
				endif
				que(1)=aa(fig%stp(mst(1,i)),fig%ttp(mst(1,i)))%change_label([lbd(3+i/4:3+i/4+j),lbd(1)],[lbd(3+i/4:3+i/4+j)//"1","m "])
				que(2)=aa(fig%stp(mst(2,i)),fig%ttp(mst(2,i)))%change_label([lbd(3+i/4:3+i/4+j),lbd(2)],[lbd(3+i/4:3+i/4+j)//"2","m "])
				que(3)=P(stP(1,i))
				call get_corder(que(1:3+j),cord)
				a(i)=fig_contract(que(1:3+j),cord)
				a(i)%is_return=.true.
			enddo
			do i=1,4
				aa(fig%stp(mst(2,i)),fig%ttp(mst(1,i)))=a(i)
			enddo
			deallocate(que)
		enddo
		grp(:,:,gtmp:gtmp+1)=0

	end subroutine
	subroutine cvCTM(fig,label,aa)
		type(t_fig) :: fig
		character(*) :: label(:)
		type(t_tensor) :: aa(:,-1:)
		type(t_tensor) :: R(2),E(2,2),Q(2),a(4),P(6,size(label)),tmp
		integer :: i,j,n,k,ist(4,2),istp(4,2),mst(4,4),st,Pst(2,3),stP(2,4)
		character(:), allocatable :: bd(:,:),lbd(:)
		real(8), allocatable :: U(:,:),V(:,:),Up(:,:),Vp(:,:)
		real(8), allocatable :: s(:)
		type(t_tensor), allocatable :: que(:)
		integer, allocatable :: cord(:),mgrp(:,:)
		Pst=reshape([1,2,2,3,3,4],[2,3])
		stP=reshape([1,0,2,3,4,5,6,0],[2,4])
		do n=1,size(label)
			allocate(que(4))
			select case(label(n))
			case("r")
				ist(:,1)=[1,2,4,3]
				ist(:,2)=[13,14,16,15]
				mst=reshape([1,2,4,3,5,6,8,7,9,10,12,11,13,14,16,15],[4,4])
				bd=reshape(["13","14","16","15","21","22","24","23"],[4,2])
				lbd=["r","l","d","u"]
			case("l")
				ist(:,1)=[4,3,1,2]
				ist(:,2)=[16,15,13,14]
				mst=reshape([4,3,1,2,8,7,5,6,12,11,9,10,16,15,13,14],[4,4])
				bd=reshape(["16","15","13","14","24","23","21","22"],[4,2])
				lbd=["l","r","d","u"]
			case("d")
				ist(:,1)=[1,5,13,9]
				ist(:,2)=[4,8,16,12]
				mst=reshape([1,5,13,9,2,6,14,10,3,7,15,11,4,8,16,12],[4,4])
				bd=reshape(["1","4","10","7","3","6","12","9"],[4,2])
				lbd=["d","u","r","l"]
			case("u")
				ist(:,1)=[13,9,1,5]
				ist(:,2)=[16,12,4,8]
				mst=reshape([13,9,1,5,14,10,2,6,15,11,3,7,16,12,4,8],[4,4])
				bd=reshape(["10","7","1","4","12","9","3","6"],[4,2])
				lbd=["u","d","r","l"]
			end select
			do i=1,3
				if((.not.allocated(P((i-1)*2+1,n)%shap))) then
					do k=1,2
						istp(:,k)=mod(fig%stp(ist(:,k))+1+fig%stp(mst(1,Pst(k,i)))-fig%stp(ist(1,k)),2)+1
						do j=1,4
							que(j)=aa(istp(j,k),fig%ttp(ist(j,k)))%change_label(fig%label(fig%st(ist(j,k))%ltp),to_char(fig%st(ist(j,k))%bd))
						enddo
						call get_corder(que(1:2),cord)
						E(k,1)=change_label(fig_contract(que(1:2),cord),bd(1:2,k),["c1","c2"])
						call get_corder(que(3:4),cord)
						E(k,2)=change_label(fig_contract(que(3:4),cord),bd(3:4,k),["c1","c2"])
						call E(k,1)%fact(["c1","c2"],[character::],Q(1),R(1),"qr",["m1","m1"],gtmp)
						call E(k,2)%fact([character::],["c1","c2"],R(2),Q(2),"rq",["m2","m2"],gtmp)
						call svd_t(dot(R(1),R(2)),["m1"],["m2"],U,s,V)
						do j=1,size(s)
							s(j)=sqrt(s(j))
							U(:,j)=s(j)*U(:,j)
							V(j,:)=s(j)*V(j,:)
						enddo
						!E(k,1)=change_label(dot(Q(1),U,["m1"]),["m1"],[lbd(5-k)])
						R(1)=dot(Q(1),U,["m1"])
						call R(1)%fact(["c1","c2"],[character::],E(k,1),Q(1),"rq",[lbd(5-k),lbd(5-k)],gtmp+1)
						if(size(label)>1) then
							!E(k,2)=change_label(dot(V,Q(2),["m2"]),["m2"],[lbd(5-k)])
							R(2)=dot(V,Q(2),["m2"])
							call R(2)%fact([character::],["c1","c2"],Q(2),E(k,2),"qr",[lbd(5-k),lbd(5-k)],gtmp+1)
						endif
					enddo
					do k=1,size(label)
						!call svd_t(dot(E(1,k),E(2,k)),[lbd(4)],[lbd(3)],U,s,V,grp(:,:,gchi),"s")
						call svd_t(dot(E(1,k),E(2,k)),[lbd(4)],[lbd(3)],U,s,V,grp(:,:,gchi))
						do j=1,chi
							s(j)=1d0/sqrt(s(j))
							U(:,j)=s(j)*U(:,j)
							V(j,:)=s(j)*V(j,:)
						enddo
						!call E(1,k)%get_mat([lbd(4)],["c1","c2"],Up)
						!call E(2,k)%get_mat(["c1","c2"],[lbd(3)],Vp)
						!!write(*,*)matmul(matmul(Vp,transpose(V)),matmul(transpose(U),transpose(Up)))
						!write(*,*)sum(abs(matmul(Vp,matmul(matmul(transpose(V),transpose(U)),Up))-diag(1d0,size(Up,2))))
						!write(*,*)sum(abs(matmul(matmul(transpose(V),transpose(U)),matmul(Up,Vp))-diag(1d0,size(Vp,2))))
						P((i-1)*2+1,k)=dot(V(:chi,:),E(2,k)%change_label(["c1","c2"],lbd(3)//["1","2"]),[lbd(3)],gchi)
						P(i*2,k)=dot(E(1,k)%change_label(["c1","c2"],lbd(4)//["1","2"]),U(:,:chi),[lbd(4)],gchi)
						!R(1)=dot(P((i-1)*2+1,k),P(i*2,k),[lbd(3),lbd(4)])
						!call R(1)%reorder([lbd(3)//["1","2"],lbd(4)//["1","2"]])
						!read(*,*)
						!write(*,*)sum(abs(R(1)%rc%T-reshape(diag(1d0,product(R(1)%shap(1:2))),shape(R(1)%rc%T))))
						!read(*,*)
						do j=i+1,3
							if(fig%stp(mst(1+(k-1)*2,Pst(1,j)))==fig%stp(mst(1+(k-1)*2,Pst(1,i)))) then
								P((j-1)*2+1,k)=P((i-1)*2+1,k)
								P(j*2,k)=P(i*2,k)
							endif
						enddo
					enddo
				endif
			enddo

			P(:,n)%is_return=.true.

			do i=1,4
				if(stP(2,i)/=0) then
					j=1
					que(4)=P(stP(2,i),n)
				else
					j=0
				endif
				que(1)=aa(fig%stp(mst(1,i)),fig%ttp(mst(1,i)))%change_label([lbd(3+i/4:3+i/4+j),lbd(1)],[lbd(3+i/4:3+i/4+j)//"1","m "])
				que(2)=aa(fig%stp(mst(2,i)),fig%ttp(mst(2,i)))%change_label([lbd(3+i/4:3+i/4+j),lbd(2)],[lbd(3+i/4:3+i/4+j)//"2","m "])
				que(3)=P(stP(1,i),n)
				call get_corder(que(1:3+j),cord)
				a(i)=fig_contract(que(1:3+j),cord)
				a(i)%is_return=.true.
			enddo
			do i=1,4
				aa(fig%stp(mst(2,i)),fig%ttp(mst(1,i)))=a(i)
			enddo
			deallocate(que)
		enddo
		grp(:,:,gtmp:gtmp+1)=0

	end subroutine
	subroutine simp_update(H,delta0,A,n)
		real(8) :: H(:,:)
		real(8) :: delta0(:)
		type(t_tensor) :: A(:)
		integer :: n
		type(t_tensor) :: b(2),Q(2),E
		real(8) :: ld(d,2*size(A)),alpha
		type(t_fig) :: fig
		integer :: i,j,step,ist(2),m,i1,i2
		real(8) :: ld1(d**3),ld2(d**3),ld0(d**4),delta
		real(8), allocatable :: s(:),He(:)
		real(8), allocatable :: U(:,:),V(:,:),Q1(:,:),Q2(:,:),R(:,:),L(:,:)
		real(8) :: HU(size(H,1),size(H,2)),Hexp(size(H,1),size(H,2))
		integer :: bord1(4),bord2(4),lord1(4),lord2(4),g2(9,2)
		character(:), allocatable :: lgate,lbd(:)
		integer, allocatable :: mgrp(:,:)

		call smerge([gp,gd],mgrp=mgrp)
		grp(:,:,gpd)=mgrp

		g2(:,1)=[4,3,4,2,1,2,4,3,4]
		g2(:,2)=[4,2,4,2,1,2,4,2,4]

		call random_number(ld)
		do i=1,size(ld,2)
			ld(:,i)=ld(:,i)/sqrt(sum(ld(:,i)**2))
		enddo

		HU=H
		!call allocate(He,size(HU,1))
		!call heev(HU,He,"V")
		call smat_eg([1,1],HU,He)

		call gen_fig(fig,[1,1,2, 1,2,1, 2,1,2, 2,2,1])
		fig%stp=[1,2]

		do i=1,size(A)
			ld0=1d0/sqrt(outprod(outprod(ld(:,fig%st(i)%bd(1)),ld(:,fig%st(i)%bd(2))),outprod(ld(:,fig%st(i)%bd(3)),ld(:,fig%st(i)%bd(4)))))
			A(i)=dot(ld0,A(i),fig%label(fig%st(i)%ltp))
		enddo

		do step=1,n
			!delta=(delta0(1)-delta0(2))*(1d0/(1d0+exp(8d0*step/n-4d0))-1d0/(1d0+exp(4d0)))/(1d0/(1d0+exp(-4d0))-1d0/(1d0+exp(4d0)))+delta0(2)
			!delta=max(delta0(1)-(delta0(1)-delta0(2))/(n/2d0)*step,delta0(2))
			!delta=delta0(1)-(delta0(1)-delta0(2))/n*step
			delta=delta0(1)*(exp(-6d0*step/n)-exp(-6d0))+delta0(2)
			Hexp=matmul(matmul(HU,diag(exp(-He*delta))),transpose(HU))
			!Hexp=diag(1d0,size(Hexp,1))-delta*H
			!do i=1,size(ld,2)
			do m=1,4
				if((step==1.and.m<=2).or.(step==n.and.m>2)) then
					cycle
				endif
				if(mod(step,2)==0) then
					i=m
				else
					i=5-m
				endif
				select case(fig%btp(i))
				case(1)
					lgate="d"
					lbd=["r","l"]
				case(2)
					lgate="l"
					lbd=["d","u"]
				end select

				bord1=fig%st(fig%bd(1,i))%bd
				lord1=fig%st(fig%bd(1,i))%ltp
				bord2=fig%st(fig%bd(2,i))%bd
				lord2=fig%st(fig%bd(2,i))%ltp
				do j=1,3
					if(fig%st(fig%bd(1,i))%bd(j)==i) then
						bord1([j,4])=bord1([4,j])
						lord1([j,4])=lord1([4,j])
					endif
					if(fig%st(fig%bd(2,i))%bd(j)==i) then
						bord2([j,4])=bord2([4,j])
						lord2([j,4])=lord2([4,j])
					endif
				enddo
				ld1=outprod(outprod(ld(:,bord1(1)),ld(:,bord1(2))),ld(:,bord1(3)))
				ld2=outprod(outprod(ld(:,bord2(1)),ld(:,bord2(2))),ld(:,bord2(3)))
				ist=fig%stp(fig%bd(:,i))

				A(ist(fig%btp(i)))=dot(fgate([gd,gp]),A(ist(fig%btp(i))),[lgate,"p"])

				A(ist(1))=dot(ld1,A(ist(1)),fig%label(lord1(1:3)))
				A(ist(2))=dot(A(ist(2)),ld2,fig%label(lord2(1:3)))
				call A(ist(1))%fact([character::],["p",lbd(1)],Q(1),b(1),"qr",lbd,gpd)
				call A(ist(2))%fact([lbd(2),"p"],[character::],b(2),Q(2),"lq",lbd,gpd)
				E=dot(dot(dot(b(1)%change_label(["p"],["p1"]),ld(:,i),[lbd(1)]),b(2)%change_label(["p"],["p2"]),lbd),Hexp,["p1","p2"])
				call E%svd([lbd(2),"p1"],[lbd(1),"p2"],U,s,V,grp(:,:,gd))
				ld(:,i)=s(:d)/sqrt(sum(s(:d)**2))
				call b(1)%get_tensor([lbd(2),"p",lbd(1)],U(:,:d))
				call b(2)%get_tensor([lbd(2),lbd(1),"p"],V(:d,:))
				A(ist(1))=dot(Q(1),b(1),lbd)
				A(ist(2))=dot(b(2),Q(2),lbd)
				A(ist(1))=dot(1d0/ld1,A(ist(1)),fig%label(lord1(1:3)))
				A(ist(2))=dot(A(ist(2)),1d0/ld2,fig%label(lord2(1:3)))

				A(ist(fig%btp(i)))=dot(fgate([gd,gp]),A(ist(fig%btp(i))),[lgate,"p"])

			enddo
		enddo
		ld=sqrt(ld)
		do i=1,size(A)
			ld0=outprod(outprod(ld(:,fig%st(i)%bd(1)),ld(:,fig%st(i)%bd(2))),outprod(ld(:,fig%st(i)%bd(3)),ld(:,fig%st(i)%bd(4))))
			A(i)=dot(ld0,A(i),fig%label(fig%st(i)%ltp))
		enddo
	end subroutine
	subroutine full_update(E,Hexp,b)
		type(t_tensor) :: E,b(:)
		real(8) :: Hexp(:,:)
		type(t_tensor) :: Q,Qt,T,Z
		real(8), allocatable :: MS(:,:),MR(:,:),Mtmp(:,:),MZ(:,:),MQ(:,:),R(:,:),L(:,:),U(:,:),V(:,:)
		real(8), allocatable :: Eg(:),s(:)
		integer :: i,j,n
		integer, allocatable :: mgrp(:,:)
		real(8) :: cost,cost0,pcost

		!T=E%contract(["r","l","r'","l'"])
		!E%rc%T=E%rc%T/T%rc%T(1)

		call E%get_mat(["r","l"],["r'","l'"],MZ)
		MZ=0.5d0*(MZ+transpose(MZ))
		!call allocate(Eg,size(MZ,1))
		!call heev(MZ,Eg,"V")
		call smat_eg([gpd,gpd],MZ,Eg)
		do i=1,size(MZ,2)
			if(Eg(i)>0d0) then
				MZ(:,i)=MZ(:,i)*sqrt(Eg(i))
			else
				MZ(:,i)=0d0
			endif
		enddo
		call smerge([gpd,gpd],mgrp=mgrp)
		grp(:,:,gtmp)=mgrp

		!! Gauge Fixing
		!call Z%new(["r","l","p"],[E%shap(1:2),size(MZ,2)],[gpd,gpd,gtmp])
		!call Z%get_tensor(["r","l","p"],MZ)
		!call Z%fact(["r","p"],["l"],MQ,R,"qr")
		!call Z%fact(["r"],["l","p"],L,MQ,"lq")
		!b(2)=dot(R,b(2),["r"])
		!b(1)=dot(b(1),L,["l"])
		!!call allocate(U,R)
		!!call allocate(V,L)
		!!call mat_inv(L)
		!call smat_inv([gpd],[gpd],L)
		!!call mat_inv(R)
		!call smat_inv([gpd],[gpd],R)
		!Z=dot(dot(L,Z,["r"]),R,["l"])
		call E%get_tensor(E%label,matmul(MZ,transpose(MZ)))
		!E=dot(Z,Z%change_label(["r","l"],["r'","l'"]))
		!call Z%clear()
		grp(:,:,gtmp)=0

		T=dot(dot(b(1)%change_label(["p"],["p1"]),b(2)%change_label(["p"],["p2"]),["r","l"]),Hexp,["p1","p2"])

		!call T%svd(["l","p1"],["p2","r"],U,s,V,grp(:,:,gd))
		!do i=1,d
			!U(:,i)=U(:,i)*sqrt(s(i))
			!V(i,:)=V(i,:)*sqrt(s(i))
		!enddo
		!call b(1)%get_tensor(["l","p","r"],U(:,:d))
		!call b(2)%get_tensor(["l","p","r"],V(:d,:))

		Q=dot(E,T,["r","l","l","r"])
		Qt=dot(E,T,["r'","l","l'","r"])
		cost0=get_value(dot(Qt,T,["p1","p1","p2","p2","r","l","l","r"]),1)
		pcost=-1d0
		n=0
		!write(*,*)"cost: "
		do 
			n=n+1
			cost=cost0
			T=dot(b(1)%change_label(["p"],["p1"]),b(2)%change_label(["p"],["p2"]),["r","l"])
			cost=cost-get_value(dot(Q,T,["p1","p1","p2","p2","r'","l","l'","r"]),1)
			cost=cost-get_value(dot(Qt,T,["p1","p1","p2","p2","r","l","l","r"]),1)
			cost=cost+get_value(dot(dot(E,T,["r","l","l","r"]),T,["p1","p1","p2","p2","r'","l","l'","r"]),1)
			if(n==1) then
				write(*,"(es12.4,' -> '$)")abs(cost)
			endif
			!write(*,"(2es14.6','2es14.6)")T%rc%T(1),cost
			!read(*,*)
			if(abs((cost-pcost)/pcost)<1d-6.or.n>100) then
				write(*,"(es12.4,i5)")abs(cost),n
				exit
			endif
			!write(*,*)cost
			pcost=cost
			do i=1,2
				j=mod(i,2)+1
				call get_mat(dot(Q,b(j)%change_label(["r","l"],["r'","l'"]),merge(["p2","p","l'","r'"],["p1","p","r'","l'"],i==1)),["r'","l'"],[character::],MS)
				call get_mat(dot(dot(E,b(j),merge(["l","r"],["r","l"],i==1)),b(j)%change_label(["r","l"],["r'","l'"]),merge(["p","p","l'","r'"],["p","p","r'","l'"],i==1)),["r'","l'"],["r","l"],MR)
				!call allocate(MZ,MR)
				!call mat_inv(MR)
				call smat_inv(b(i)%stp(b(i)%get_ldx(["l","r"])),b(i)%stp(b(i)%get_ldx(["l","r"])),MR)
				!call smat_eg(b(i)%stp(b(i)%get_ldx(lbd(2:1:-1))),b(i)%stp(b(i)%get_ldx(lbd(2:1:-1))),MR,Eg)
				!MR=matmul(matmul(MR,diag(merge(0d0,1d0/Eg,abs(Eg)<1d-8))),transpose(MR))
				cost=sum(abs(matmul(MZ,MR)-diag(1d0,size(MR))))
				!if(abs(cost)>1d-5) then
					!write(*,*)"inverse: ",cost
					!!stop
				!endif
				call b(i)%get_tensor(["l","r","p"],matmul(MR,MS))
			enddo
		enddo
		!b(2)=dot(R,b(2),["r"])
		!b(1)=dot(b(1),L,["l"])
	end subroutine
	subroutine rescale(fig,aa)
		type(t_fig) :: fig
		type(t_tensor) :: aa(:,-1:)
		integer :: i,loc(1),j
		real(8) :: scal
		!do i=1,size(aa,1)
			!do j=0,8
				!loc=maxloc(abs(aa(i,j)%rc%T))
				!if(j==0) then
					!scal=norm/abs(aa(i,j)%rc%T(loc(1)))
					!aa(i,-1)%rc%T=aa(i,-1)%rc%T*(scal**0.5)
				!else
					!scal=norm/aa(i,j)%rc%T(loc(1))
				!endif
				!aa(i,j)%rc%T=aa(i,j)%rc%T*scal
			!enddo
		!enddo

		do i=1,size(aa,1)
			do j=0,8
				loc=maxloc(abs(aa(i,j)%rc%T))
				scal=norm/abs(aa(i,j)%rc%T(loc(1)))
				aa(i,j)%rc%T=aa(i,j)%rc%T*scal
				if(j==0) then
					aa(i,-1)%rc%T=aa(i,-1)%rc%T*(scal**0.5)
				endif
			enddo
		enddo

		!do i=1,size(aa,1)
			!do j=0,8
				!if(j==0) then
					!scal=norm/sqrt(real(sum(aa(i,-1)%rc%T*aa(i,-1)%rc%T)))
					!aa(i,-1)%rc%T=aa(i,-1)%rc%T*scal
					!scal=scal**2
				!else
					!scal=norm/sqrt(real(sum(aa(i,j)%rc%T*aa(i,j)%rc%T)))
				!endif
				!aa(i,j)%rc%T=aa(i,j)%rc%T*scal
				!!if(j==0.or.j==1.or.j==5) then
					!!write(*,*)i,j,scal
				!!endif
			!enddo
		!enddo
		!write(*,*)"----------------------------------"
	end subroutine
	function get_vsite(fig,aa,st,H,sg) result(rt)
		type(t_fig) :: fig
		type(t_tensor) :: aa(:,-1:)
		integer :: st(:)
		real(8) :: H(:,:)
		real(8), optional :: sg(:)
		type(t_tensor), allocatable :: que(:)
		real(8) :: rt
		real(8) :: rt_(size(st))
		integer :: i,n,pos
		type(t_tensor) :: E
		character(:), allocatable :: label(:)
		type(t_dat), pointer :: dat
		integer, allocatable :: cord(:)
		label=["u","u'","r","r'","d","d'","l","l'"]
		do n=1,size(st)
			pos=0
			allocate(que(16))
			do i=1,size(fig%st)
				if(i==st(n)) then
					pos=pos+1
					que(pos)=change_label(merge_label(&
						dot(dot(fgate([gd,gd,gd]),&
						dot(aa(fig%stp(i),-1),aa(fig%stp(i),-1)%change_label([label(::2),"p"],[label(2::2),"P "]))&
						,["u'","l'","l"]),fgate([gd,gd,gd]),["r'","d'","d"])&
						,label,label(::2)),fig%label(fig%st(i)%ltp),to_char(fig%st(i)%bd))
				else
					pos=pos+1
					que(pos)=aa(fig%stp(i),fig%ttp(i))%change_label(fig%label(fig%st(i)%ltp),to_char(fig%st(i)%bd))
				endif
			enddo
			call get_corder(que,cord)
			E=fig_contract(que,cord)
			rt_(n)=get_value(contract(dot(E,H,["p"]),["p","P"]),1)
			deallocate(que)
		enddo
		E=contract(E,["p","P"])
		rt_=rt_/E%rc%T(1)
		write(*,*)rt_,E%rc%T(1)
		if(present(sg)) then
			rt=sum(real(rt_)*sg)
		else
			rt=sum(real(rt_))
		endif
		call E%clear()
	end function
	function get_vbond(fig,aa,bd,H,sg) result(rt)
		type(t_fig) :: fig
		type(t_tensor) :: aa(:,-1:)
		integer :: bd(:)
		real(8) :: H(:,:)
		real(8), optional :: sg(:)
		type(t_tensor), allocatable :: que(:)
		real(8) :: rt
		real(8) :: rt_(size(bd))
		integer :: i,n,k,fst(2),pos
		type(t_tensor) :: E
		character(:), allocatable :: label(:),lgate
		integer, allocatable :: cord(:)
		label=["u","u'","r","r'","d","d'","l","l'"]
		rt=0d0
		do n=1,size(bd)
			pos=0
			allocate(que(16))
			do i=1,size(fig%st)
				do k=1,size(fig%st(i)%bd)
					if(fig%st(i)%bd(k)==bd(n)) then
						exit
					endif
				enddo
				if(k/=(size(fig%st(i)%bd)+1)) then
					select case(fig%label(fig%st(i)%ltp(k)))
					case("r")
						lgate="d"
					case("l")
						lgate="u'"
					case("d")
						lgate="r'"
					case("u")
						lgate="l"
					end select
					pos=pos+1
					que(pos)=change_label(merge_label(&
						dot(fgate([gd,gp,gp]),&
						dot(dot(fgate([gd,gd,gd]),&
						dot(aa(fig%stp(i),-1)%change_label(["p"],["p"//to_char(i)]),aa(fig%stp(i),-1)%change_label([label(::2),"p"],[label(2::2),"P"//to_char(i)]))&
						,["u'","l'","l"]),fgate([gd,gd,gd]),["r'","d'","d"])&
						,[lgate,"p"//to_char(i),"P"//to_char(i)])&
						,label,label(::2)),fig%label(fig%st(i)%ltp),to_char(fig%st(i)%bd))
				else
					pos=pos+1
					que(pos)=aa(fig%stp(i),fig%ttp(i))%change_label(fig%label(fig%st(i)%ltp),to_char(fig%st(i)%bd))
				endif
			enddo
			call get_corder(que,cord)
			E=fig_contract(que,cord)
			rt_(n)=get_value(contract(dot(E,H,"p"//[to_char(fig%bd(1,bd(n))),to_char(fig%bd(2,bd(n)))]),[["p","P"]//to_char(fig%bd(1,bd(n))),["p","P"]//to_char(fig%bd(2,bd(n)))]),1)
			deallocate(que)
		enddo
		E=contract(E,[["p","P"]//to_char(fig%bd(1,bd(size(bd)))),["p","P"]//to_char(fig%bd(2,bd(size(bd))))])
		rt_=rt_/E%rc%T(1)
		write(*,*)rt_,E%rc%T(1)
		if(present(sg)) then
			rt=sum(real(rt_)*sg)
		else
			rt=sum(real(rt_))
		endif
		call E%clear()
	end function
	function fgate(legs) result(rt)
		integer :: legs(:)
		real(8), allocatable :: rt(:)
		integer :: n(size(legs)),m(size(legs)),prod(size(legs)+1),idx(size(legs)),i,k,l
		k=size(legs)
		m=grp(1,2,abs(legs(:)))
		n=grp(2,2,abs(legs(:)))
		allocate(rt(product(m+n)))

		prod(1)=1
		l=1
		do i=1,k
			prod(i+1)=prod(i)*(m(i)+n(i))
			if(legs(i)<0) then
				l=i
			endif
		enddo

		do i=1,size(rt)
			rt(i)=1d0
			if(sign(1,-mod(i-1,prod(l+1))/prod(l)/m(l))==-1) then
				if(product(sign(1,-mod(i-1,prod(2:))/prod(:k)/m))==1) then
					rt(i)=-1d0
				endif
			endif
		enddo
	end function
	function is_sym(self) result(rt)
		type(t_tensor) :: self
		integer, allocatable :: map(:),mgrp(:,:)
		logical :: rt
		call smerge(gstp(self%stp),map,mgrp)
		rt=(sum(abs(self%rc%T(map(mgrp(1,2)+1:))))<1d-17)
	end function
end module


program main
	use M_fig
	implicit none
	real(8) :: delta=0.05d0,alpha=10d0,g=3.10d0
	real(8), allocatable :: HU(:,:),Q(:,:),R(:,:),L(:,:),H(:,:),U(:,:),V(:,:),Hexp(:,:),Hn(:,:),Hdb(:,:)
	real(8), allocatable :: s(:),ps(:),He(:)
	real(8) :: val
	type(t_tensor), target :: aa(2,-1:9),X(2),b(2),tmp1,tmp2,tH
	type(t_tensor) :: E
	type(t_fig) :: fig,figx,figy
	type(t_tensor), allocatable :: que(:)
	type(t_dat), pointer :: dat
	integer :: flag(8)=[1,1,2,2,3,3,4,4]
	character(:), allocatable :: label(:),lbd(:)
	integer :: step,n,i,j,k,bd,st(2),nctm
	integer, allocatable :: map(:),mgrp(:,:)
	integer, allocatable :: cord(:)
	!open(10,file=fn("../data/check.dat"))

	chi=min(d*10,64)

	!call random_seed()
	call omp_set_nested(.false.)
	!call omp_set_max_active_levels(1)
	call omp_set_dynamic(.false.)
	call mkl_set_dynamic(0)
	call omp_set_schedule(omp_sched_static,0)

	call mkl_set_num_threads(1)
	call omp_set_num_threads(mkl_get_max_threads())

	allocate(grp(2,2,8))
	grp=0

	!! spinless tight-binding
	!p=2
	!grp(:,:,gp)=reshape([1,-1, 1,1],[2,2])
	!call tH%new(["u0","u1","d0","d1"],[p,p,p,p],[gp,gp,gp,gp],"0")
	!tH%rc%T(tH%get_idx([1,2,2,1, 2,1,1,2],["u0","u1","d0","d1"]))=-t(1)*[1,1]
	!tH%rc%T(tH%get_idx([1,2,1,2, 2,1,2,1, 2,2,2,2],["u0","u1","d0","d1"]))=[-mu/4d0,-mu/4d0,DV-mu/2d0]
	!allocate(Hn(p,p),Hdb(p,p),Hexp(p**2,p**2))
	!Hn=0d0; Hn(1,1)=1d0

	! Hubbard
	p=4
	grp(:,:,gp)=reshape([1,-1, 2,2],[2,2])
	call tH%new(["u0","u1","d0","d1"],[p,p,p,p],[gp,gp,gp,gp],"0")
	tH%rc%T(tH%get_idx([4,3,2,1, 4,3,1,2, 2,3,3,2, 2,4,4,2, 2,1,4,3, 1,2,4,3, 3,2,2,3, 4,2,2,4],["u0","u1","d0","d1"]))=t(1)*[1,1,1,1,1,1,1,1]
	tH%rc%T(tH%get_idx([3,4,2,1, 1,3,3,1, 1,4,4,1, 3,4,1,2, 2,1,3,4, 3,1,1,3, 4,1,1,4, 1,2,3,4],["u0","u1","d0","d1"]))=-t(1)*[1,1,1,1,1,1,1,1]
	tH%rc%T(tH%get_idx([2,1,2,1, 1,2,1,2, 2,2,2,2, 3,2,3,2, 4,2,4,2, 2,3,2,3, 2,4,2,4],["u0","u1","d0","d1"]))=0.25d0*DU*[1,1,2,1,1,1,1]
	tH%rc%T(tH%get_idx([2,1,2,1, 3,1,3,1, 4,1,4,1, 1,2,1,2, 2,2,2,2, 3,2,3,2, 4,2,4,2, 1,3,1,3, 2,3,2,3, 3,3,3,3, 4,3,4,3, 1,4,1,4, 2,4,2,4, 3,4,3,4, 4,4,4,4],["u0","u1","d0","d1"]))=tH%rc%T(tH%get_idx([2,1,2,1, 3,1,3,1, 4,1,4,1, 1,2,1,2, 2,2,2,2, 3,2,3,2, 4,2,4,2, 1,3,1,3, 2,3,2,3, 3,3,3,3, 4,3,4,3, 1,4,1,4, 2,4,2,4, 3,4,3,4, 4,4,4,4],["u0","u1","d0","d1"]))-0.25d0*mu*[2,1,1,2,4,3,3,1,3,2,2,1,3,2,2]
	allocate(Hn(p,p),Hdb(p,p),Hexp(p**2,p**2))
	Hn=0d0; Hn(2,2)=2d0; Hn(3,3)=1d0; Hn(4,4)=1d0
	Hdb=0d0; Hdb(2,2)=1d0

	!! t-J
	!p=3
	!grp(:,:,gp)=reshape([1,-1, 1,2],[2,2])
	!call tH%new(["u0","u1","d0","d1"],[p,p,p,p],[gp,gp,gp,gp],"0")
	!tH%rc%T(tH%get_idx([1,2,2,1, 2,1,1,2, 1,3,3,1, 3,1,1,3],["u0","u1","d0","d1"]))=-t(1)*[1,1,1,1]
	!tH%rc%T(tH%get_idx([2,2,2,2, 3,2,3,2, 2,3,2,3, 3,3,3,3],["u0","u1","d0","d1"]))=[DV+0.25d0*DJ,DV-0.25d0*DJ,DV-0.25d0*DJ,DV+0.25d0*DJ]
	!tH%rc%T(tH%get_idx([2,3,3,2, 3,2,2,3],["u0","u1","d0","d1"]))=[0.5d0*DJ,0.5d0*DJ]
	!tH%rc%T(tH%get_idx([2,1,2,1, 3,1,3,1, 1,2,1,2, 2,2,2,2, 3,2,3,2, 1,3,1,3, 2,3,2,3, 3,3,3,3],["u0","u1","d0","d1"]))=tH%rc%T(tH%get_idx([2,1,2,1, 3,1,3,1, 1,2,1,2, 2,2,2,2, 3,2,3,2, 1,3,1,3, 2,3,2,3, 3,3,3,3],["u0","u1","d0","d1"]))-0.25d0*mu*[1,1,1,2,2,1,2,2]
	!allocate(Hn(p,p),Hdb(p,p),Hexp(p**2,p**2))
	!Hn=0d0; Hn(2,2)=1d0; Hn(3,3)=1d0
	!Hdb=0d0

	grp(:,:,gd)=reshape([1,-1, d/2,d-d/2],[2,2])
	grp(:,:,gchi)=reshape([1,-1, chi/2,chi-chi/2],[2,2])

	call tH%get_mat(["u0","u1"],["d0","d1"],H)
	call allocate(HU,H)
	!call syev(HU,He,"V")
	!call heev(HU,He,"V")
	call smat_eg([1,1],HU,He)

	label=["u","u'","r","r'","d","d'","l","l'"]

	do i=1,size(aa,1)
		call aa(i,-1)%new([label(::2),"p"],[d,d,d,d,p],[gd,gd,gd,gd,gp],"r")
		!aa(i,-1)%rc%T=real(aa(i,-1)%rc%T)
	enddo
	!call aa(1,-1)%new([label(::2),"p"],[d,d,d,d,p],[gd,gd,gd,gd,gp],"r")
	!aa(2,-1)=aa(1,-1)%clone()

	call simp_update(H,[0.05d0,0.05d0],aa(:,-1),1000)

    !tp:        
    !
	! 5---1---6
	! |   |   |
	! 4--0/9--2
	! |   |   |
	! 8---3---7
	!

	do i=1,size(aa,1)
		aa(i,0)=dot(aa(i,-1),aa(i,-1)%change_label(label(::2),label(2::2)))
		aa(i,0)=dot(dot(fgate([gd,gd,gd]),aa(i,0),["u'","l'","l"]),fgate([gd,gd,gd]),["r'","d'","d"])
		!do j=1,4
			!aa(i,j)=contract(aa(i,0),pack(label,flag==j))
			!!aa(i,j)=aa(i,j)%merge_label(pack(label,.not.((flag==j).or.(flag==mod(j+1,4)+1))),pack(label(::2),.not.((flag(::2)==j).or.(flag(::2)==mod(j+1,4)+1))))
			!aa(i,j)=aa(i,j)%merge_label(pack(label,.not.(flag==j)),pack(label(::2),.not.(flag(::2)==j)))
		!enddo
		!do j=1,4
			!aa(i,j+4)=contract(aa(i,0),pack(label,flag==j.or.flag==mod(j-2+4,4)+1))
			!aa(i,j+4)=aa(i,j+4)%merge_label(pack(label,.not.(flag==j.or.flag==mod(j-2+4,4)+1)),&
				!pack(label(::2),.not.(flag(::2)==j.or.flag(::2)==mod(j-2+4,4)+1)))
		!enddo
		call aa(i,1)%new(["l","r","d"],[chi,chi,d*d],[gchi,gchi,gdd],"r")
		call aa(i,2)%new(["u","d","l"],[chi,chi,d*d],[gchi,gchi,gdd],"r")
		call aa(i,3)%new(["l","r","u"],[chi,chi,d*d],[gchi,gchi,gdd],"r")
		call aa(i,4)%new(["u","d","r"],[chi,chi,d*d],[gchi,gchi,gdd],"r")
		call aa(i,5)%new(["r","d"],[chi,chi],[gchi,gchi],"r")
		call aa(i,6)%new(["l","d"],[chi,chi],[gchi,gchi],"r")
		call aa(i,7)%new(["l","u"],[chi,chi],[gchi,gchi],"r")
		call aa(i,8)%new(["u","r"],[chi,chi],[gchi,gchi],"r")
		aa(i,0)=aa(i,0)%merge_label(label,label(::2))
		!do j=1,8
			!aa(i,j)%rc%T=real(aa(i,j)%rc%T)
		!enddo
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

	call gen_fig(fig,[1,1,2, 1,2,3, 1,3,4, 1,5,6, 1,6,7, 1,7,8, 1,9,10, 1,10,11, 1,11,12, 1,13,14, 1,14,15, 1,15,16,&
		2,1,5, 2,2,6, 2,3,7, 2,4,8, 2,5,9, 2,6,10, 2,7,11, 2,8,12, 2,9,13, 2,10,14, 2,11,15, 2,12,16])
	fig%stp=[1,2,1,2,2,1,2,1,1,2,1,2,2,1,2,1]
	!fig%ttp([6,7,10,11])=-1

	call rescale(fig,aa(:,-1:8))

	delta=0.05d0
	do step=0,Nstep
		!delta=0.05d0*exp(-alpha*step/Nstep)
		!delta=0.1d0
		!delta=max(0.1d0-(0.1d0-0.1d0)/(Nstep)*step,0.1d0)
		!if(mod(step+1,30)==0) then
			!delta=delta/2d0
		!endif
		if(mod(step+1,10)==0) then
			delta=max(delta*0.8d0,0.001d0)
		endif
		Hexp=matmul(matmul(HU,diag(exp(-He*delta))),transpose(HU))
		!Hexp=diag(1d0,size(Hexp,1))
		write(*,"('step: ',i4,es12.4)")step,delta
		nctm=0
		do 
			do n=1,4

				select case(mod(step,2)*5+sign(1,-mod(step,2))*n)
				case(1:2)
					bd=5
					lbd=["r","l"]
				case(3:4)
					bd=18
					lbd=["d","u"]
				end select
				st=fig%bd(:,bd)

				!if(.false.) then
				if(nctm==0.and.step/=0.and.(.not.((step==1.and.n<=2).or.(step==Nstep.and.n>2)))) then
					if(fig%btp(bd)==1) then
						aa(fig%stp(st(1)),-1)=dot(fgate([gd,gp]),aa(fig%stp(st(1)),-1),["d","p"])
					else                                                       
						aa(fig%stp(st(2)),-1)=dot(fgate([gd,gp]),aa(fig%stp(st(2)),-1),["l","p"])
					endif

					call aa(fig%stp(st(1)),-1)%fact([character::],[lbd(1),"p"],X(1),b(1),"qr",lbd,gpd)
					call aa(fig%stp(st(2)),-1)%fact([lbd(2),"p"],[character::],b(2),X(2),"lq",lbd,gpd)

					do i=1,2
						aa(fig%stp(st(i)),9)=dot(X(i),X(i)%change_label(label(::2),label(2::2)))
						aa(fig%stp(st(i)),9)=dot(dot(fgate(aa(fig%stp(st(i)),9)%stp(aa(fig%stp(st(i)),9)%get_ldx(["u'","l'","l"]))),aa(fig%stp(st(i)),9),["u'","l'","l"]),fgate(aa(fig%stp(st(i)),9)%stp(aa(fig%stp(st(i)),9)%get_ldx(["r'","d'","d"]))),["r'","d'","d"])
						aa(fig%stp(st(i)),9)=aa(fig%stp(st(i)),9)%merge_label(label,label(::2))
					enddo

					fig%ttp(fig%bd(:,bd))=9


					call ini_que(fig,[bd],[integer::],aa,que)
					call get_corder(que,cord)
					E=fig_contract(que,cord)
					E=E%split_label([lbd(1)//to_char(bd),lbd(2)//to_char(bd)],["r","r'","l","l'"])
					deallocate(que)


					if(fig%btp(bd)==2) then
						do i=1,2
							b(i)=b(i)%change_label(["u","d"],["l","r"])
						enddo
					endif
					call full_update(E,Hexp,b)
					if(fig%btp(bd)==2) then
						do i=1,2
							b(i)=b(i)%change_label(["l","r"],["u","d"])
						enddo
					endif

					aa(fig%stp(st(1)),-1)=dot(X(1),b(1),lbd)
					aa(fig%stp(st(2)),-1)=dot(b(2),X(2),lbd)
					if(fig%btp(bd)==1) then
						aa(fig%stp(st(1)),-1)=dot(fgate([gd,gp]),aa(fig%stp(st(1)),-1),["d","p"])
					else                                                       
						aa(fig%stp(st(2)),-1)=dot(fgate([gd,gp]),aa(fig%stp(st(2)),-1),["l","p"])
					endif
					do i=1,2
						aa(fig%stp(st(i)),9)=dot(aa(fig%stp(st(i)),-1),aa(fig%stp(st(i)),-1)%change_label(label(::2),label(2::2)))
						aa(fig%stp(st(i)),9)=dot(dot(fgate([gd,gd,gd]),aa(fig%stp(st(i)),9),["u'","l'","l"]),fgate([gd,gd,gd]),["r'","d'","d"])
						aa(fig%stp(st(i)),9)=aa(fig%stp(st(i)),9)%merge_label(label,label(::2))
					enddo
				endif
				if(fig%btp(bd)==1) then
					!call CTM(fig,reshape([13,14,17,18,21,22, 16,15,20,19,24,23],[6,2]),aa)
					!call dCTMp(fig,["r","l"],aa)
					!if(step<=20) then
						call vCTM(fig,["r","l"],aa)
						!call vCTMp(fig,["r","l"],aa)
					!else
						!call cvCTM(fig,["r","l"],aa)
					!endif
				else
					!call CTM(fig,reshape([1,4,2,5,3,6, 10,7,11,8,12,9],[6,2]),aa)
					!call dCTMp(fig,["u","d"],aa)
					!if(step<=20) then
						call vCTM(fig,["u","d"],aa)
						!call vCTMp(fig,["u","d"],aa)
					!else
						!call cvCTM(fig,["u","d"],aa)
					!endif
				endif

				if(fig%ttp(st(1))==9) then
					fig%ttp(st)=0
					do i=1,2
						aa(fig%stp(st(i)),0)=aa(fig%stp(st(i)),9)
					enddo
				endif
				!do i=1,size(fig%stp)
					!write(*,*)aa(fig%stp(i),fig%ttp(i))%label,aa(fig%stp(i),fig%ttp(i))%shap,aa(fig%stp(i),fig%ttp(i))%stp
				!enddo

				fig%stp=mod(fig%stp,2)+1

				call rescale(fig,aa(:,-1:8))
			enddo
			nctm=nctm+1
			if((step<Nstep-2.and.nctm>4).or.(step>=Nstep-2.and.nctm>30)) then
				exit
			endif
		enddo
		val=get_vsite(fig,aa,[6,7,10,11],Hn)/4d0
		write(*,*)"n: ",val
		!write(10,"(i4,es12.4$)")step,real(val)
		val=get_vbond(fig,aa,[5,8,18,19],H)/2d0+mu*val
		write(*,*)"Energy: ",val
		!write(10,"(es12.4)")real(val)
		!write(*,*)"Mz: ",get_vsite(fig,aa,[6,7,10,11],cmplx(reshape([1d0,0d0,0d0,-1d0],[2,2]),kind=8))/4d0
	enddo
end program

