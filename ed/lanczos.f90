module global
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    �趨����ģ�Ͳ�������ʼ��                                           !
!		����:                                                           !
!     		dn*dn�������С                                             !
!     		nt����������                                                !
!     		ns��������z����                                             !
!     		t��hoppingֵ                                                !
!     		u��hubbard u                                                !
!     		pbc�����ڸ����Ϣ                                           !
!    	�ӳ���:                                                         !
!     		init_lattice�������i��pbc����Ϊ������iλ��Ϊ1������Ϊ��    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	save
	integer(8),parameter :: dn=4,dn2=16,nt=2
	real(8),parameter :: pi=3.14159265358979d0
	integer(8) :: pbc(dn2,0:1)=0,ns=0
	real(8) :: t=-1d0,u=4d0
	contains
	subroutine init_lattice()
		implicit none
		integer :: i
		pbc(dn2,0:1)=0
		do i=1,dn2
			pbc(i,:)=ibset(0,i-1)
			pbc(i,0)=ibset(pbc(i,0),mod(mod(i-1,dn)+1,dn)+(i-1)/dn*dn)
			pbc(i,1)=ibset(pbc(i,1),mod(i+dn-1,dn2))
		enddo
	end subroutine
end module
module lanczos
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       lanczos����                                                                                                  !
!       ������global                                                                                                 !
!       ������                                                                                                       !
!       	ad��au�������£������ϵ���������                                                                         !
!       	ca����ԭ���Ļ����ϲ���һ�����ӵ�������������Ҫ�ֶ��ͷ��ڴ�                                               !
!       	abitd��abitu���ֱ�ָ�������£�����������������ָ��                                                       !
!       	abit��ָ��abitd��abitu��ָ�룬���ڲ���̬                                                                 !
!       �ӳ���                                                                                                     !
!       	initialstate����ʼ����������                                                                             !
!       		������                                                                                               !
!       			a���������������                                                                                !
!       			sz�����������������С                                                                           !
!       			n�����Ӹ���                                                                                      !
!       	findstate������ֵ̬����Ӧ��ָ�꣬���ҵ�����������ָ��abitָ��                                            !
!       		������                                                                                               !
!       			a�����룺���ҵ�ֵ̬                                                                              !
!       			sz�����룺���ҵ�����������С                                                                     !
!       			n�������̬��ָ��                                                                                !
!       	hubbard�����ܶ������õ������ϣ�ʹ��ǰ��Ҫָ��abitd��abitu��hopping������pbcʵ��                          !
!       		������                                                                                               !
!       			va�����룺���õ�����                                                                             !
!       			vb�����룺û�г�ʼ������������õĽ��                                                           !
!       			sz��������ά��                                                                                   !
!       	lanc����׼lanczos���̣���Ҫ����hubbard�ͶԽǻ�����                                                       !
!       		������                                                                                               !
!       			v�����룺v(:,1:2)Ϊ��һ���ĳ�ʼ�����������v(:,2)Ϊ�����������Ѿ���һ����                        !
!       			e�����룺�����������������ֵ                                                                    !
!       			m�����룺����lanczos���̴�������ap��ά���������������                                         !
!       			mo�������ʵ��ִ��lanczos���̵Ĵ�������e����Чά��                                               !
!       			flag�����룺�Ƿ���㱾��������.true.Ϊ���㣬���������                                           !
!       			info��������Խǻ��ɹ���ָʾ��Ϣ                                                                 !
!       	create��������һ��̬�ϲ���������һ�����ӵ�̬����Ҫ�õ�ָ��abitu��abitd���������������һָ��ָ��ca     !
!       		������                                                                                               !
!       			v�����룺��ʼ̬����̬�������������                                                              !
!       			cv�������ָ�룬ָ����̬����Ҫ�ֶ��ͷ��ڴ�                                                     !
!       			i��spin�����룺�������ӵ����������                                                              !
!       			sz�����룺ԭʼ̬��ά������������̬��ά��                                                       !
!       			flag�����룺������������.true.Ϊ����                                                           !
!       	calcu�����㶯��ѧ�������ڴ��ͷ�cv���ڴ棩����Ҫ�õ�ָ��abitu��abitd                                      !
!       		������                                                                                               !
!       			cv���ض��ĳ�ʼlanczos̬v����create�еĽ����                                                     !
!       			sz�����룺̬��ά��                                                                               !
!       			m�����룺lanczos���̵Ĵ���                                                                       !
!       			fq�����룺�����Ƶ������                                                                         !
!       			y�����룺�����㣬���������Ķ���ѧ�����������                                                  !
!       			img�����룺ȡ����չ��                                                                            !
!       	init_random_seed�����������                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use global
	implicit none
	integer(8),pointer :: abit(:),abitd(:),abitu(:)
	integer(8),allocatable,target :: ad(:),au(:),ca(:)
	contains
	subroutine initialstate(p,sz,n)
		implicit none
		integer(8),allocatable,target :: p(:)
		integer(8) :: tmp(10000000),i,n,j,k,sz
		tmp=0
		do i=1,n
			tmp(1)=ibset(tmp(1),i-1)
		enddo
		i=1
ex:		do while(.true.)
			tmp(i+1)=tmp(i)
			i=i+1
			do j=0,dn2-1
				if(j==dn2-1) then
					exit ex
				endif
				if(ibits(tmp(i),j,2)==1) then
					tmp(i)=ibset(tmp(i),j+1)
					tmp(i)=ibclr(tmp(i),j)
					do k=0,j/2-1
						if(btest(tmp(i),k)) then
							exit
						endif
						if(btest(tmp(i),j-k-1)) then
							tmp(i)=ibset(tmp(i),k)
							tmp(i)=ibclr(tmp(i),j-k-1)
						endif
					enddo
					exit
				endif
			enddo
		enddo ex
		sz=i-1
		allocate(p(sz))
		p=tmp(1:sz)
	end subroutine
	subroutine findstate(a,n,sz)
		implicit none
		integer(8) :: a,n,i,ia,ib,sz
		ia=1
		ib=sz+1
		i=(ia+ib)/2
		do while(a/=abit(i))
			if(a<abit(i)) then
				ib=i
			else
				ia=i
			endif
			i=(ia+ib)/2
			! write(*,*)ia,i,ib
		enddo
		n=i
	end subroutine
	subroutine hubbard(va,vb,sz)
		implicit none
		integer(8) :: i,ii,iii,ia(0:1),m,k,l,a,j,ion,aon,n,tmp,sz(2),sig
		logical :: flag
		real(8) :: va(:),vb(:)
		!$omp parallel do reduction(+:vb) private(ia,aon,ion,m,k,l,a,j,tmp,abit,flag,sig) schedule(guided)
		do n=1,sz(1)*sz(2)
			ia(0)=(n-1)/sz(2)+1
			ia(1)=n-(n-1)/sz(2)*sz(2)
			aon=iand(abitd(ia(0)),abitu(ia(1)))
			ion=0
			do i=1,dn2
				do ii=0,3
					if(ii/2==0) then
						abit=>abitd
					else
						abit=>abitu
					endif
					flag=.false.
					sig=1
					m=pbc(i,mod(ii,2))
					k=iand(abit(ia(ii/2)),m)
					if(k==m.or.k==0) then
						cycle
					endif
					l=ieor(k,m)
					a=abit(ia(ii/2))-k+l
					tmp=ia(ii/2)
					call findstate(a,ia(ii/2),sz(ii/2+1))
					do iii=0,dn2-1
						if(btest(m,iii)) then
							flag=.not.flag
						elseif(flag.and.btest(abit(tmp),iii)) then
							sig=-1*sig
						endif
					enddo
					j=(ia(0)-1)*sz(2)+ia(1)
					vb(j)=vb(j)+t*sig*va(n)
					ia(ii/2)=tmp
				enddo
				if(btest(aon,i-1)) then
					ion=ion+1
				endif
			enddo
			vb(n)=vb(n)+u*ion*va(n)
		enddo
		!$omp end parallel do
	end subroutine
	subroutine lanc(v2,e,sz,m,mo,flag,info)
		implicit none
		integer(8) :: sz(2),modi(0:1)
		integer :: info,is,m,mo
		real(8) :: v(sz(1)*sz(2),0:1),ap(m),bt(0:m),e(:),tmp(0:m),fe(m),work(2*m-2),sun,v2(:)
		real(8),allocatable :: z(:,:),ev(:)
		logical :: flag
		v(:,1)=v2
		sun=dot_product(v(:,1),v(:,1))
		v(:,1)=v(:,1)/sqrt(abs(sun))
		v(:,0)=0d0
		fe(1)=1000
		bt(0)=0d0
		do is=1,m
			modi(1)=mod(is,2)
			modi(0)=mod(is-1,2)
			v(:,modi(0))=-1*bt(is-1)*v(:,modi(0))
			call hubbard(v(:,modi(1)),v(:,modi(0)),sz)
			ap(is)=dot_product(v(:,modi(1)),v(:,modi(0)))
			if(is>3) then
				e(1:is)=ap(1:is)
				tmp(1:is-1)=bt(1:is-1)
				call dsteqr('n',is,e(1:is),tmp(1:is-1),work,is,work,info)
				if(abs(e(1)-fe(1))<1e-12.or.info/=0.or.bt(is-1)<1e-8) then
					if(flag) then
						e(1:is)=ap(1:is)
						tmp(1:is-1)=bt(1:is-1)
						allocate(z(is,is),ev(is))
						call dsteqr('i',is,e(1:is),tmp(1:is-1),z,is,work,info)
						ev=z(:,1)
					endif
					mo=is
					write(*,*)"info=",info
					write(*,"(2e17.6,i5)")bt(is-1),abs(e(1)-fe(1)),is
					exit
				endif
				fe=e
			endif
			v(:,modi(0))=v(:,modi(0))-ap(is)*v(:,modi(1))
			sun=dot_product(v(:,modi(0)),v(:,modi(0)))
			bt(is)=sqrt(abs(sun))
			v(:,modi(0))=v(:,modi(0))/bt(is)
		enddo
		if(flag) then
			v(:,0)=0d0
			v(:,1)=v2
			v2=0d0
			do is=1,mo
				modi(1)=mod(is,2)
				modi(0)=mod(is-1,2)
				v2=v2+ev(is)*v(:,modi(1))
				v(:,modi(0))=-1*bt(is-1)*v(:,modi(0))
				call hubbard(v(:,modi(1)),v(:,modi(0)),sz)
				v(:,modi(0))=v(:,modi(0))-ap(is)*v(:,modi(1))
				v(:,modi(0))=v(:,modi(0))/bt(is)
			enddo
			sun=dot_product(v2,v2)
			v2=v2/sqrt(sun)
		endif
	end subroutine
	subroutine create(v,cv,i,spin,sz,flag)
		implicit none
		integer(8) :: i,spin,sz(2),ii,j,n,sig,tmp,x,c,ci
		real(8),pointer :: cv(:)
		real(8) :: v(:)
		logical :: flag
		if(flag) then
			c=1
		else
			c=-1
		endif
		ci=ibset(0,i-1)
		if(spin==0) then
			tmp=sz(1)
			call initialstate(ca,sz(1),(nt-ns)/2+c)
			allocate(cv(sz(2)*sz(1)))
			cv=0d0
			abit=>ca
			do n=1,tmp
				if(flag.neqv.btest(abitd(n),i-1)) then
					sig=-1**((nt+ns)/2)
					x=ieor(abitd(n),ci)
					call findstate(x,j,sz(1))
					do ii=i,dn2-1
						if(btest(abitd(n),ii)) then
							sig=-1*sig
						endif
					enddo
					cv((j-1)*sz(2)+1:j*sz(2))=cv((j-1)*sz(2)+1:j*sz(2))+v((n-1)*sz(2)+1:n*sz(2))*sig
				endif
			enddo
			abitd=>ca
		else
			tmp=sz(2)
			call initialstate(ca,sz(2),(nt+ns)/2+c)
			allocate(cv(sz(2)*sz(1)))
			cv=0d0
			abit=>ca
			do n=1,tmp
				if(flag.neqv.btest(abitu(n),i-1)) then
					sig=1
					x=ieor(abitu(n),ci)
					call findstate(x,j,sz(2))
					do ii=i,dn2-1
						if(btest(abitu(n),ii)) then
							sig=-1*sig
						endif
					enddo
					cv(j:sz(2)*(sz(1)-1)+j:sz(2))=cv(j:sz(2)*(sz(1)-1)+j:sz(2))+v(n:tmp*(sz(1)-1)+n:tmp)*sig
				endif
			enddo
			abitu=>ca
		endif
	end subroutine
	subroutine calcu(cv,sz,m,fq,y,img)
		implicit none
		integer(8) :: sz(2),modi(0:1)
		integer :: is,m,i,j
		real(8) :: sun,fq(:),y(:),img,pii,ap(1:m),bt(0:m),v(sz(1)*sz(2),0:1)
		real(8),pointer :: cv(:)
		complex(8) :: x,z
		v(:,0)=0d0
		v(:,1)=cv
		sun=dot_product(v(:,1),v(:,1))
		v(:,1)=v(:,1)/sqrt(abs(sun))
		bt(0)=0d0
		pii=-1d0/pi
		do is=1,m
			modi(1)=mod(is,2)
			modi(0)=mod(is-1,2)
			v(:,modi(0))=-1*bt(is-1)*v(:,modi(0))
			call hubbard(v(:,modi(1)),v(:,modi(0)),sz)
			ap(is)=dot_product(v(:,modi(1)),v(:,modi(0)))
			v(:,modi(0))=v(:,modi(0))-ap(is)*v(:,modi(1))
			sun=dot_product(v(:,modi(0)),v(:,modi(0)))
			bt(is)=sqrt(abs(sun))
			if(bt(is)<1e-8.and.is>10) then
				write(*,"(e17.6,i5)")bt(is),is
				exit
			endif
			v(:,modi(0))=v(:,modi(0))/bt(is)
		enddo
		bt(0)=dot_product(cv,cv)
		deallocate(cv)
		do i=1,size(fq)
			x=0d0
			z=cmplx(fq(i),img)
			do j=is,1,-1
				x=bt(j-1)**2/(z-ap(j)-x)
			enddo
			y(i)=y(i)+dimag(x)*pii
		enddo
		deallocate(ca)
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


program main
	use global
	use lanczos
	implicit none
	real(8),allocatable :: v2(:),v_low(:),e(:),lowe(:)
	real(8),pointer :: cv(:)
	integer(8) :: i,j,k,sz(2),tmp(2),nslow
	real(8),allocatable :: fq(:),y(:)
	real(8) :: omin,omax
	integer :: m=200,info,mo
	logical :: flag
	open(unit=10,file="../data/output.dat")
	call init_random_seed()
	call init_lattice()
	allocate(v_low(1))
	allocate(e(m),lowe(m))
	lowe(1)=100
	do ns=nt,0,-2
		call initialstate(ad,sz(1),(nt-ns)/2)
		if(ns==0) then
			sz(2)=sz(1)
			allocate(au(sz(2)))
			au=ad
		else
			call initialstate(au,sz(2),(nt+ns)/2)
		endif
		allocate(v2(sz(1)*sz(2)))
		call random_number(v2)
		abitd=>ad
		abitu=>au
		call lanc(v2,e,sz,m,mo,.true.,info)
		if(e(1)<lowe(1)) then
			deallocate(lowe,v_low)
			allocate(lowe(mo),v_low(sz(1)*sz(2)))
			lowe=e(1:mo)
			nslow=ns
			v_low=v2
		endif
		deallocate(v2,ad,au)
	enddo
	ns=nslow
	call initialstate(ad,sz(1),(nt-ns)/2)
	if(ns==0) then
		sz(2)=sz(1)
		allocate(au(sz(2)))
		au=ad
	else
		call initialstate(au,sz(2),(nt+ns)/2)
	endif
	write(10,"(e17.6)")lowe
	write(10,*)"----------------------------------------------"
	if(nt==0) then
		v_low=1d0
	endif
	allocate(fq(1237),y(1237))
	m=100
	deallocate(e)
	allocate(e(m))
	! abitu=>au
	! abitd=>ad
	! allocate(v2(sz(1)*sz(2)))
	! v2=0d0
	! call hubbard(v_low,v2,sz)
	! write(10,"(e17.6)")v_low
	! write(10,*)"----------------------------------------------"
	flag=.true.
	do k=0,3
		flag=.not.flag
		write(*,*)flag,k,k/2
		y=0d0
		tmp=sz
		abitu=>au
		abitd=>ad
		call create(v_low,cv,k/10+1,k/2,sz,flag)
		call lanc(cv,e,sz,m,mo,.false.,info)
		omin=e(1)-0.5
		omax=e(mo)+0.5
		do i=1,size(fq)
			fq(i)=omin+(omax-omin)/size(fq)*i
		enddo
		call calcu(cv,sz,m,fq,y,0.02d0)
		sz=tmp
		write(10,"(2e17.6)")(fq(j)-lowe(1),y(j),j=1,size(fq))
	enddo
end
