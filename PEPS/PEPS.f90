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
											if((dat%cost(cm)<=dat3%cost(cm).and.dat%cost(mod(cm,2)+1)<dat3%cost(mod(cm,2)+1)).or.(c==n.and.dat3%cost(mod(cm,2)+1)/dat%cost(mod(cm,2)+1)>10.and.dat3%cost(mod(cm,2)+1)>40000000)) then
											!if((dat%cost(cm)<=dat3%cost(cm).and.dat%cost(mod(cm,2)+1)<dat3%cost(mod(cm,2)+1))) then
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
			!!read(*,*)
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
				if(any(i==rm)) cycle
				label=[character::]
				to_label=[character::]
				select case(fig%ttp(i))
				case(-1)
					pos=pos+2
					if(k==2) then
						label=fig%label(fig%st(i)%ltp)
						to_label=cut_label(fig,i,cut)
						que(pos-1)=aa(fig%stp(i),fig%ttp(i))%change_label([label,"p"],[to_label,"p"//to_char(i)])
						que(pos)=set_conjg(aa(fig%stp(i),fig%ttp(i))%change_label([label,"p"],[(trim(to_label(n))//"'",n=1,size(to_label)),"p"//to_char(i)]))


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
		complex(8), allocatable :: U(:,:),V(:,:)
		real(8), allocatable :: s(:)
		type(t_tensor), allocatable :: que(:)
		type(t_tensor) :: aa_(size(ist,2))
		integer, allocatable :: cord(:)
		!write(*,*)pair
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
				call E%svd([lbd(1)//to_char(pair(n1,n2)),lbd(1)//to_char(pair(n1+1,n2))],[lbd(2)//to_char(pair(n1,n2)),lbd(2)//to_char(pair(n1+1,n2))],U,s,V)
				call P(1)%new(lbd(1)//["1","2",""],[E%shap(1:2),min(chi,product(E%shap(1:2)))])
				call P(1)%get_tensor([P(1)%label],U(:,:P(1)%shap(3)))
				P(1)=set_conjg(P(1))
				P(2)=set_conjg(P(1)%change_label(lbd(1)//["1","2",""],lbd(2)//["1","2",""]))
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

		!read(*,*)
	end subroutine
	subroutine dCTM(fig,label,aa)
		type(t_fig) :: fig
		character(*) :: label(:)
		type(t_tensor) :: aa(:,-1:)
		type(t_tensor) :: C(2),E(2),a(4)
		integer :: i,j,n,ist(8),l(8),st
		character(:), allocatable :: bd(:),bdp(:),lbd(:)
		complex(8), allocatable :: U1(:,:),U2(:,:),Z(:,:),W(:,:),cZ(:,:),cW(:,:)
		real(8), allocatable :: Eg(:)
		type(t_tensor), allocatable :: que(:,:)
		integer, allocatable :: cord(:)
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

			call get_mat(dot(C(1),set_conjg(C(1)%change_label(bd(1:2),bd(1:2)//"'"))),bd(1:2),bd(1:2)//"'",U1)
			call get_mat(dot(set_conjg(C(2)),C(2)%change_label(bd(5:6),bd(5:6)//"'")),bd(5:6),bd(5:6)//"'",U2)
			call allocate(Z,U1+U2)
			call que(3,1)%new([bd(1:2),lbd(1)],[C(1)%shap(C(1)%get_ldx(bd(1:2))),chi])
			call allocate(Eg,size(Z,1))
			call heev(Z,Eg,"V")
			call que(3,1)%get_tensor(que(3,1)%label,Z(:,size(Z,2)-chi+1:))
			que(3,1)=set_conjg(que(3,1))
			que(3,2)=set_conjg(que(3,1)%change_label([bd(1:2),lbd(1)],[bd(5:6),lbd(2)]))
			que(4,1)=set_conjg(que(3,1)%change_label([lbd(1)],[lbd(2)]))
			que(4,2)=que(3,1)%change_label(bd(1:2),bd(5:6))

			call get_mat(dot(E(1),set_conjg(E(1)%change_label(bd(3:4),bd(3:4)//"'"))),bd(3:4),bd(3:4)//"'",U1)
			call get_mat(dot(set_conjg(E(2)),E(2)%change_label(bd(3:4),bd(3:4)//"'")),bd(3:4),bd(3:4)//"'",U2)
			call allocate(W,U1+U2)
			call que(7,1)%new([bd(3:4),lbd(1)],[E(1)%shap(E(1)%get_ldx(bd(3:4))),chi])
			call allocate(Eg,size(W,1))
			call heev(W,Eg,"V")
			call que(7,1)%get_tensor(que(7,1)%label,W(:,size(W,2)-chi+1:))
			que(7,2)=que(7,1)%change_label([lbd(1)],[lbd(2)])
			que(7,1)=set_conjg(que(7,1))

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
	subroutine simp_update(H,delta0,A,n)
		complex(8) :: H(:,:)
		real(8) :: delta0(:)
		type(t_tensor) :: A(:)
		integer :: n
		type(t_tensor) :: App,b(2),Q(2)
		real(8) :: ld(d,2*size(A)),alpha
		type(t_fig) :: fig
		integer :: i,j,step,ist(2)
		real(8) :: ld1(d**3),ld2(d**3),ld0(d**4),delta
		real(8), allocatable :: s(:),He(:)
		complex(8), allocatable :: U(:,:),V(:,:),Q1(:,:),Q2(:,:),R(:,:),L(:,:)
		complex(8) :: HU(size(H,1),size(H,2)),Hexp(size(H,1),size(H,2))
		integer :: bord1(4),bord2(4),lord1(4),lord2(4)
		character(:), allocatable :: lgate,lbd(:)

		call random_number(ld)

		HU=H
		call allocate(He,size(HU,1))
		call heev(HU,He,"V")

		call gen_fig(fig,[1,1,2, 1,2,1, 2,1,2, 2,2,1])
		fig%stp=[1,2]

		do i=1,size(A)
			ld0=1d0/sqrt(outprod(outprod(ld(:,fig%st(i)%bd(1)),ld(:,fig%st(i)%bd(2))),outprod(ld(:,fig%st(i)%bd(3)),ld(:,fig%st(i)%bd(4)))))
			A(i)=dot(cmplx(ld0,kind=8),A(i),fig%label(fig%st(i)%ltp))
		enddo

		ist=[1,2]
		do step=1,n
			!delta=(delta0(1)-delta0(2))*(1d0/(1d0+exp(8d0*step/n-4d0))-1d0/(1d0+exp(4d0)))/(1d0/(1d0+exp(-4d0))-1d0/(1d0+exp(4d0)))+delta0(2)
			!delta=max(delta0(1)-(delta0(1)-delta0(2))/(n/2d0)*step,delta0(2))
			!delta=delta0(1)-(delta0(1)-delta0(2))/n*step
			delta=delta0(1)*(exp(-6d0*step/n)-exp(-6d0))+delta0(2)
			Hexp=matmul(matmul(HU,diag(exp(-He*delta))),transpose(conjg(HU)))
			do i=1,size(ld,2)

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

				select case(fig%btp(i))
				case(1)
					lbd=["r","l"]
				case(2)
					lbd=["d","u"]
				end select

				A(ist(1))=dot(cmplx(ld1,kind=8),A(ist(1)),fig%label(lord1(1:3)))
				A(ist(2))=dot(A(ist(2)),cmplx(ld2,kind=8),fig%label(lord2(1:3)))
				call A(ist(1))%qr(fig%label(lord1(1:3)),["p",fig%label(lord1(4))],Q(1),b(1),lbd)
				call A(ist(2))%lq([fig%label(lord2(4)),"p"],fig%label(lord2(1:3)),Q(2),b(2),lbd)
				App=dot(dot(dot(b(1)%change_label(["p"],["p1"]),cmplx(ld(:,i),kind=8),[lbd(1)]),b(2)%change_label(["p"],["p2"]),lbd),Hexp,["p1","p2"])
				call App%svd([lbd(2),"p1"],["p2",lbd(1)],U,s,V)
				ld(:,i)=s(:d)/sqrt(sum(s(:d)**2))
				call b(1)%get_tensor([lbd(2),"p",lbd(1)],U(:,:d))
				call b(2)%get_tensor([lbd(2),"p",lbd(1)],V(:d,:))
				A(ist(1))=dot(Q(1),b(1),lbd)
				A(ist(2))=dot(b(2),Q(2),lbd)
				A(ist(1))=dot(cmplx(1d0/ld1,kind=8),A(ist(1)),fig%label(lord1(1:3)))
				A(ist(2))=dot(A(ist(2)),cmplx(1d0/ld2,kind=8),fig%label(lord2(1:3)))
			enddo
		enddo
		ld=sqrt(ld)
		do i=1,size(A)
			ld0=outprod(outprod(ld(:,fig%st(i)%bd(1)),ld(:,fig%st(i)%bd(2))),outprod(ld(:,fig%st(i)%bd(3)),ld(:,fig%st(i)%bd(4))))
			A(i)=dot(cmplx(ld0,kind=8),A(i),fig%label(fig%st(i)%ltp))
		enddo
	end subroutine
	!subroutine simp_update(H,delta0,A,n)
		!complex(8) :: H(:,:)
		!real(8) :: delta0(:)
		!type(t_tensor) :: A(:)
		!integer :: n
		!type(t_tensor) :: App,b(2)
		!real(8) :: ld(d,2*size(A)),alpha
		!type(t_fig) :: fig
		!integer :: i,j,step,ist(2)
		!real(8) :: ld1(d**3),ld2(d**3),ld0(d**4),delta
		!real(8), allocatable :: s(:),He(:)
		!complex(8), allocatable :: U(:,:),V(:,:),Q1(:,:),Q2(:,:),R(:,:),L(:,:)
		!complex(8) :: HU(size(H,1),size(H,2)),Hexp(size(H,1),size(H,2))
		!integer :: bord1(4),bord2(4),lord1(4),lord2(4)
		!character(:), allocatable :: lgate

		!call random_number(ld)

		!HU=H
		!call allocate(He,size(HU,1))
		!call heev(HU,He,"V")

		!call gen_fig(fig,[1,1,2, 1,2,1, 2,1,2, 2,2,1])
		!fig%stp=[1,2]

		!do i=1,size(A)
			!ld0=1d0/sqrt(outprod(outprod(ld(:,fig%st(i)%bd(1)),ld(:,fig%st(i)%bd(2))),outprod(ld(:,fig%st(i)%bd(3)),ld(:,fig%st(i)%bd(4)))))
			!A(i)=dot(cmplx(ld0,kind=8),A(i),fig%label(fig%st(i)%ltp))
		!enddo

		!ist=[1,2]
		!call b(1)%new(["l1","p1","r1"],[d*p,p,d])
		!call b(2)%new(["l2","p2","r2"],[d,p,d*p])
		!do step=1,n
			!!delta=(delta0(1)-delta0(2))*(1d0/(1d0+exp(8d0*step/n-4d0))-1d0/(1d0+exp(4d0)))/(1d0/(1d0+exp(-4d0))-1d0/(1d0+exp(4d0)))+delta0(2)
			!!delta=max(delta0(1)-(delta0(1)-delta0(2))/(n/2d0)*step,delta0(2))
			!!delta=delta0(1)-(delta0(1)-delta0(2))/n*step
			!delta=delta0(1)*(exp(-6d0*step/n)-exp(-6d0))+delta0(2)
			!Hexp=matmul(matmul(HU,diag(exp(-He*delta))),transpose(conjg(HU)))
			!do i=1,size(ld,2)

				!bord1=fig%st(fig%bd(1,i))%bd
				!lord1=fig%st(fig%bd(1,i))%ltp
				!bord2=fig%st(fig%bd(2,i))%bd
				!lord2=fig%st(fig%bd(2,i))%ltp
				!do j=1,3
					!if(fig%st(fig%bd(1,i))%bd(j)==i) then
						!bord1([j,4])=bord1([4,j])
						!lord1([j,4])=lord1([4,j])
					!endif
					!if(fig%st(fig%bd(2,i))%bd(j)==i) then
						!bord2([j,4])=bord2([4,j])
						!lord2([j,4])=lord2([4,j])
					!endif
				!enddo
				!ld1=outprod(outprod(ld(:,bord1(1)),ld(:,bord1(2))),ld(:,bord1(3)))
				!ld2=outprod(outprod(ld(:,bord2(1)),ld(:,bord2(2))),ld(:,bord2(3)))
				!ist=fig%stp(fig%bd(:,i))

				!A(ist(1))=dot(cmplx(ld1,kind=8),A(ist(1)),fig%label(lord1(1:3)))
				!A(ist(2))=dot(A(ist(2)),cmplx(ld2,kind=8),fig%label(lord2(1:3)))

				!call A(ist(1))%qr(fig%label(lord1(1:3)),["p",fig%label(lord1(4))],Q1,R)
				!call b(1)%get_tensor(["l1","p1","r1"],R)
				!call A(ist(2))%lq([fig%label(lord2(4)),"p"],fig%label(lord2(1:3)),Q2,L)
				!call b(2)%get_tensor(["l2","p2","r2"],L)

				!App=dot(dot(dot(b(1),cmplx(ld(:,i),kind=8),["r1"]),b(2),["r1","l2"]),Hexp,["p1","p2"])
				!call App%svd(["l1","p1"],["p2","r2"],U,s,V)
				!ld(:,i)=s(:d)/sqrt(sum(s(:d)**2))
				!call b(1)%get_tensor(["l1","p1","r1"],U(:,:d))
				!call b(2)%get_tensor(["l2","p2","r2"],V(:d,:))
				!do j=1,size(ld1)
					!Q1(j,:)=1d0/ld1(j)*Q1(j,:)
					!Q2(:,j)=1d0/ld2(j)*Q2(:,j)
				!enddo
				!A(ist(1))=split_label(change_label(dot(Q1,b(1),["l1"]),["p1","r1"],["p",fig%label(lord1(4))]),["l1"],fig%label(lord1(1:3)))
				!A(ist(2))=split_label(change_label(dot(b(2),Q2,["r2"]),["p2","l2"],["p",fig%label(lord2(4))]),["r2"],fig%label(lord2(1:3)))

			!enddo
		!enddo
		!ld=sqrt(ld)
		!do i=1,size(A)
			!ld0=outprod(outprod(ld(:,fig%st(i)%bd(1)),ld(:,fig%st(i)%bd(2))),outprod(ld(:,fig%st(i)%bd(3)),ld(:,fig%st(i)%bd(4))))
			!A(i)=dot(cmplx(ld0,kind=8),A(i),fig%label(fig%st(i)%ltp))
		!enddo
	!end subroutine
	subroutine full_update(E,Hexp,b)
		type(t_tensor) :: E,b(:)
		complex(8) :: Hexp(:,:)
		type(t_tensor) :: Q,Qt,T,Z
		complex(8), allocatable :: MS(:,:),MR(:,:),Mtmp(:,:),MZ(:,:),MQ(:,:),R(:,:),L(:,:),U(:,:),V(:,:)
		real(8), allocatable :: Eg(:),s(:)
		integer :: i,j,n
		complex(8) :: cost,cost0,pcost

		call E%get_mat(["r","l"],["r'","l'"],MZ)
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
		call Z%new(["r","l","p"],[E%shap(1:2),size(MZ,2)])
		call Z%get_tensor(["r","l","p"],MZ)
		call Z%qr(["r","p"],["l"],MQ,R)
		call Z%lq(["r"],["l","p"],MQ,L)
		b(2)=dot(R,b(2),["r"])
		b(1)=dot(b(1),L,["l"])
		!call allocate(U,R)
		!call allocate(V,L)
		call mat_inv(L)
		call mat_inv(R)
		Z=dot(dot(L,Z,["r"]),R,["l"])
		!!call E%get_tensor(E%label,matmul(MZ,transpose(conjg(MZ))))
		E=dot(Z,set_conjg(Z%change_label(["r","l"],["r'","l'"])))
		call Z%clear()

		T=dot(dot(b(1)%change_label(["p"],["p1"]),b(2)%change_label(["p"],["p2"]),["r","l"]),Hexp,["p1","p2"])

		call T%svd(["l","p1"],["p2","r"],U,s,V)
		do i=1,d
			U(:,i)=U(:,i)*sqrt(s(i))
			V(i,:)=V(i,:)*sqrt(s(i))
		enddo
		call b(1)%get_tensor(["l","p","r"],U(:,:d))
		call b(2)%get_tensor(["l","p","r"],V(:d,:))

		Q=dot(E,T,["r","l","l","r"])
		Qt=dot(E,set_conjg(T),["r'","l","l'","r"])
		cost0=get_value(dot(Qt,T,["p1","p1","p2","p2","r","l","l","r"]),1)
		pcost=-1d0
		n=0
		!write(*,*)"cost: "
		do 
			n=n+1
			cost=cost0
			T=dot(b(1)%change_label(["p"],["p1"]),b(2)%change_label(["p"],["p2"]),["r","l"])
			cost=cost-get_value(dot(Q,set_conjg(T),["p1","p1","p2","p2","r'","l","l'","r"]),1)
			cost=cost-get_value(dot(Qt,T,["p1","p1","p2","p2","r","l","l","r"]),1)
			cost=cost+get_value(dot(dot(E,T,["r","l","l","r"]),set_conjg(T),["p1","p1","p2","p2","r'","l","l'","r"]),1)
			!write(*,"(2es14.6','2es14.6)")T%rc%T(1),cost
			!read(*,*)
			if(abs(cost-pcost)<1d-10) then
				write(*,*)"update finished",cost,n
				exit
			endif
			!write(*,*)cost
			pcost=cost
			do i=1,2
				j=mod(i,2)+1
				call get_mat(dot(Q,set_conjg(b(j)%change_label(["r","l"],["r'","l'"])),merge(["p2","p","l'","r'"],["p1","p","r'","l'"],i==1)),["r'","l'"],[character::],MS)
				call get_mat(dot(dot(E,b(j),merge(["l","r"],["r","l"],i==1)),set_conjg(b(j)%change_label(["r","l"],["r'","l'"])),merge(["p","p","l'","r'"],["p","p","r'","l'"],i==1)),["r'","l'"],["r","l"],MR)
				!call allocate(MZ,MR)
				call mat_inv(MR)
				!call smat_eg(b(i)%stp(b(i)%get_ldx(lbd(2:1:-1))),b(i)%stp(b(i)%get_ldx(lbd(2:1:-1))),MR,Eg)
				!MR=matmul(matmul(MR,diag(merge(0d0,1d0/Eg,abs(Eg)<1d-8))),transpose(conjg(MR)))
				cost=sum(abs(matmul(MZ,MR)-diag(1d0,size(MR))))
				!if(abs(cost)>1d-5) then
					!write(*,*)"inverse: ",cost
					!!stop
				!endif
				call b(i)%get_tensor(["l","r","p"],matmul(MR,MS))
			enddo
		enddo
		b(2)=dot(R,b(2),["r"])
		b(1)=dot(b(1),L,["l"])
		call E%clear()
	end subroutine
	subroutine rescale(fig,aa)
		type(t_fig) :: fig
		type(t_tensor) :: aa(:,-1:)
		type(t_tensor), allocatable :: que(:)
		integer :: i,loc(size(fig%st)),ne
		complex(8) :: scale,cp
		real(8) :: norm
		type(t_tensor) :: E
		logical :: flag(size(aa,1),0:size(aa,2)-1)
		integer, allocatable :: cord(:)
		call ini_que(fig,[integer::],[integer::],aa,que)
		call get_corder(que,cord)
		E=fig_contract(que,cord)
		norm=1d0
		ne=0
		do i=1,size(fig%st)
			loc(i:i)=maxloc(abs(aa(fig%stp(i),fig%ttp(i))%rc%T))
			norm=norm*(abs(aa(fig%stp(i),fig%ttp(i))%rc%T(loc(i))))
			if(fig%ttp(i)/=0) then
				ne=ne+1
			endif
		enddo
		cp=(conjg(E%rc%T(1))/abs(E%rc%T(1)))**(1d0/ne)
		norm=(norm/abs(E%rc%T(1)))**(1d0/16d0)
		flag=.false.
		do i=1,size(fig%st)
			if(flag(fig%stp(i),fig%ttp(i))) then
				cycle
			endif
			flag(fig%stp(i),fig%ttp(i))=.true.

			scale=norm/abs(aa(fig%stp(i),fig%ttp(i))%rc%T(loc(i)))
			if(fig%ttp(i)==0) then
				aa(fig%stp(i),-1)%rc%T=aa(fig%stp(i),-1)%rc%T*scale**0.5d0
			else
				scale=scale*cp
			endif
			aa(fig%stp(i),fig%ttp(i))%rc%T=aa(fig%stp(i),fig%ttp(i))%rc%T*scale
			if(fig%ttp(i)==0) then
				E=merge_label(dot(aa(fig%stp(i),-1),set_conjg(aa(fig%stp(i),-1)%change_label(["u","d","l","r"],["u","d","l","r"]//"'"))),["u","u'","d","d'","l","l'","r","r'"],["u","d","l","r"])
				call E%reorder(aa(fig%stp(i),fig%ttp(i))%label)
				!write(*,*)"check: ",sum(abs(E%rc%T-aa(fig%st(i)%isub,fig%st(i)%tp)%rc%T))
			endif
		enddo
		call E%clear()
	end subroutine
	function get_vsite(fig,aa,st,H,sg) result(rt)
		type(t_fig) :: fig
		type(t_tensor) :: aa(:,-1:)
		integer :: st(:)
		complex(8) :: H(:,:)
		real(8), optional :: sg(:)
		type(t_tensor), allocatable :: que(:)
		complex(8) :: rt
		complex(8) :: rt_(size(st))
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
						dot(dot(aa(fig%stp(i),-1),H,["p"]),set_conjg(aa(fig%stp(i),-1)%change_label([label(::2)],[label(2::2)])))&
						,label,label(::2)),fig%label(fig%st(i)%ltp),to_char(fig%st(i)%bd))
				else
					pos=pos+1
					que(pos)=aa(fig%stp(i),fig%ttp(i))%change_label(fig%label(fig%st(i)%ltp),to_char(fig%st(i)%bd))
				endif
			enddo
			call get_corder(que,cord)
			rt_(n)=get_value(fig_contract(que,cord),1)
			deallocate(que)
		enddo
		write(*,*)rt_
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
		complex(8) :: H(:,:)
		real(8), optional :: sg(:)
		type(t_tensor), allocatable :: que(:)
		complex(8) :: rt
		complex(8) :: rt_(size(bd))
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
					pos=pos+1
					que(pos)=change_label(merge_label(&
						dot(aa(fig%stp(i),-1)%change_label(["p"],["p"//to_char(i)]),set_conjg(aa(fig%stp(i),-1)%change_label([label(::2),"p"],[label(2::2),"P"//to_char(i)])))&
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
		write(*,*)rt_
		if(present(sg)) then
			rt=sum(real(rt_)*sg)
		else
			rt=sum(real(rt_))
		endif
		call E%clear()
	end function
end module


program main
	use M_fig
	implicit none
	real(8) :: delta=0.05d0,alpha=10d0,g=3.10d0
	complex(8), allocatable :: HU(:,:),Q(:,:),R(:,:),L(:,:),H(:,:)
	real(8) :: He(p**2)
	real(8), allocatable :: s(:),ps(:)
	complex(8) :: Hexp(size(He),size(He)),energy,Hsx(p,p),Hsy(p,p),Hsz(p,p),sx,sy,sz,val
	type(t_tensor), target :: aa(2,-1:9),X(2),b(2),tmp1,tmp2,tH
	type(t_tensor) :: E
	type(t_fig) :: fig,figx,figy
	type(t_tensor), allocatable :: que(:)
	type(t_dat), pointer :: dat
	integer :: flag(8)=[1,1,2,2,3,3,4,4]
	character(:), allocatable :: label(:),lbd(:)
	integer :: step,n,i,j,k,bd,st(2)
	integer, allocatable :: cord(:)

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
	do i=1,size(aa,1)
		call aa(i,-1)%new([label(::2),"p"],[d,d,d,d,p],flag="r")
	enddo

	call simp_update(H,[0.5d0,0.005d0],aa(:,-1),3000)

    !tp:        
    !
	! 5---1---6
	! |   |   |
	! 4--0/9--2
	! |   |   |
	! 8---3---7
	!

	do i=1,size(aa,1)
		aa(i,0)=dot(aa(i,-1),set_conjg(aa(i,-1)%change_label(label(::2),label(2::2))))
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


	call rescale(fig,aa(:,-1:8))

	do step=0,Nstep-1
		delta=max(0.5d0-(0.5d0-0.005d0)/(Nstep/4)*step,0.001d0)
		delta=0d0
		Hexp=matmul(matmul(HU,diag(exp(-He*delta))),transpose(conjg(HU)))
		do n=1,4

			select case(n)
			case(1:2)
				bd=5
				lbd=["r","l"]
			case(3:4)
				bd=18
				lbd=["d","u"]
			end select
			st=fig%bd(:,bd)

			call aa(fig%stp(st(1)),-1)%qr([character::],[lbd(1),"p"],X(1),b(1),lbd)
			call aa(fig%stp(st(2)),-1)%lq([lbd(2),"p"],[character::],X(2),b(2),lbd)

			do i=1,2
				aa(fig%stp(st(i)),9)=dot(X(i),set_conjg(X(i)%change_label(label(::2),label(2::2))))
				aa(fig%stp(st(i)),9)=aa(fig%stp(st(i)),9)%merge_label(label,label(::2))
			enddo

			fig%ttp(fig%bd(:,bd))=9


			call ini_que(fig,[bd],[integer::],aa,que)
			call get_corder(que,cord)
			E=fig_contract(que,cord)
			E=E%split_label([lbd(1)//to_char(bd),lbd(2)//to_char(bd)],["r","r'","l","l'"])


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
			do i=1,2
				aa(fig%stp(st(i)),9)=dot(aa(fig%stp(st(i)),-1),set_conjg(aa(fig%stp(st(i)),-1)%change_label(label(::2),label(2::2))))
				aa(fig%stp(st(i)),9)=aa(fig%stp(st(i)),9)%merge_label(label,label(::2))
			enddo

			if(fig%btp(bd)==1) then
				!call CTM(fig,reshape([13,14,17,18,21,22, 16,15,20,19,24,23],[6,2]),aa)
				call dCTM(fig,["r","l"],aa)
			else
				!call CTM(fig,reshape([1,4,2,5,3,6, 10,7,11,8,12,9],[6,2]),aa)
				call dCTM(fig,["d","u"],aa)
			endif

			if(fig%ttp(st(1))==9) then
				fig%ttp(st)=0
				do i=1,2
					aa(fig%stp(st(i)),0)=aa(fig%stp(st(i)),9)
				enddo
			endif
			fig%stp=mod(fig%stp,2)+1

			call rescale(fig,aa(:,-1:8))
		enddo

		val=get_vbond(fig,aa,[5,8,18,19],H)/2d0
		write(*,*)"Energy: ",val
		!write(*,*)"Mz: ",get_vsite(fig,aa,[6,7,10,11],cmplx(reshape([1d0,0d0,0d0,-1d0],[2,2]),kind=8))/4d0
	enddo
end program

