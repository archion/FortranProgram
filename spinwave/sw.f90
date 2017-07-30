include "../test/hamilton_final.f90"
module sw
	use lapack95, only : geev,heev
	use M_utility
	use M_matrix
	use M_hamilton_final_M
	type(t_ham) :: ham
	real(8), allocatable :: var(:),fi(:),th(:)
contains
	subroutine diag_B(H,E,info)
		! in Nanbum space
		complex(8) :: H(:,:)
		real(8) :: E(:)
		integer, optional :: info
		complex(8) :: W(size(H,1)),Vr(size(H,1),size(H,2)),L(size(H,1),size(H,2)),Z(size(H,1),size(H,2))
		integer :: Ns,i,j,rg(2),info_,ord(size(H,1))
		integer, allocatable :: c(:)
		logical :: flag(size(H,1),size(H,2))
		Ns=size(H,1)/2
		!if(sum(abs(H(:Ns,:Ns)-conjg(H(Ns+1:,Ns+1:))))>1d-6.or.sum(abs((H(Ns+1:,:Ns)-transpose(H(Ns+1:,:Ns)))))>1d-6.or.sum(abs(conjg(transpose(H))-H))>1d-6) then
			!write(*,*)"H error"
			!stop
		!endif
		Z(:Ns,:)=H(:Ns,:)
		Z(Ns+1:,:)=-H(Ns+1:,:)
		call geev(Z,W,vr=Vr,info=info_)
		if(any(abs(imag(W))>1d-6)) then
			write(*,*)"E is not real"
		endif
		if(info_/=0) then
			write(*,*)info_
			stop
		endif
		ord=[1:size(ord)]
		call qsort(-real(W),ord)
		W=W(ord)
		Vr(:,:)=Vr(:,ord)

		do i=1,size(W)
			if(abs(real(W(i)))<1d-6) then
				Vr(:,i)=0d0
			endif
		enddo
		L(:Ns,:)=Vr(:Ns,:)
		L(Ns+1:,:)=-Vr(Ns+1:,:)
		L=matmul(conjg(transpose(Vr)),L)
		call collect(real(W),c)
		do i=1,size(c)-1
			rg=[c(i),c(i+1)-1]
			if(abs(real(W(rg(1))))<1d-6) then
				E(rg(1):rg(2))=1d0
				L(:,rg(1):rg(2))=0d0
				do j=rg(1),rg(2)
					L(j,j)=1d0
				enddo
				cycle
			endif
			if(rg(1)==rg(2)) then
				E(rg(1))=real(L(rg(1),rg(1)))
				L(rg(1),rg(1))=1d0
			else
				call heev(L(rg(1):rg(2),rg(1):rg(2)),E(rg(1):rg(2)),"v")
			endif
			do j=rg(1),rg(2)
				L(rg(1):rg(2),j)=L(rg(1):rg(2),j)/sqrt(abs(E(j)))
			enddo
		enddo
		if(present(info)) then
			Z=H
			do i=1,size(flag,1)
				do j=1,size(flag,2)
					if(i==j) then
						flag(i,j)=.true.
					else
						flag(i,j)=.false.
					endif
				enddo
			enddo
		endif
		H=matmul(Vr,L)
		!write(*,"(8i2)")nint(sign(1d0,E))
		E=sign(1d0,E)*real(W)
		if(present(info)) then
			if(sum(abs(matmul(transpose(conjg(H)),matmul(Z,H))-diag(E)))/size(H,1)>1d-6) then
				write(*,*)"diag err1",sum(abs(matmul(transpose(conjg(H)),matmul(Z,H))-diag(E)))
				stop
			endif
			if(sum(abs(matmul(transpose(conjg(H)),matmul(diag([(merge(1d0,-1d0,i<=ham%Hi/2),i=1,ham%Hi)]),H))),.not.flag)>1d-6.or.&
				any(abs(pack(abs(matmul(transpose(conjg(H)),matmul(diag([(merge(1d0,-1d0,i<=ham%Hi/2),i=1,ham%Hi)]),H))),flag)-1d0)>1d-6.and.abs(pack(abs(matmul(transpose(conjg(H)),matmul(diag([(merge(1d0,-1d0,i<=ham%Hi/2),i=1,ham%Hi)]),H))),flag))>1d-6)) then
				write(*,*)"diag err2"
				stop
			endif
		endif
	end subroutine
	!subroutine diag_B(H,E,info)
		!! in Nanbum space
		!complex(8) :: H(:,:)
		!real(8) :: E(:)
		!integer, optional :: info
		!complex(8) :: W(size(H,1)),Vr(size(H,1),size(H,2)),L(size(H,1),size(H,2)),Z(size(H,1),size(H,2))
		!integer :: Ns,i,j,rg(2),info_,ord(size(H,1))
		!integer, allocatable :: c(:)
		!logical :: flag(size(H,1),size(H,2))
		!Ns=size(H,1)/2
		!if(sum(abs(H(:Ns,:Ns)-conjg(H(Ns+1:,Ns+1:))))>1d-6.or.sum(abs((H(Ns+1:,:Ns)-transpose(H(Ns+1:,:Ns)))))>1d-6.or.sum(abs(conjg(transpose(H))-H))>1d-6) then
			!write(*,*)"H error"
			!stop
		!endif
		!Z(:Ns,:)=H(:Ns,:)
		!Z(Ns+1:,:)=-H(Ns+1:,:)
		!call geev(Z,W,vr=Vr,info=info_)
		!if(any(abs(imag(W))>1d-6)) then
			!write(*,*)"E is not real"
		!endif
		!if(info_/=0) then
			!write(*,*)info_
			!stop
		!endif
		!ord=[1:size(ord)]
		!call qsort(-real(W),ord)
		!W=W(ord)
		!Vr(:,:)=Vr(:,ord)

		!W(Ns+1:)=-W(:Ns)
		!do i=1,size(W)
			!if(abs(real(W(i)))<1d-6) then
				!Vr(:,i:Ns)=0d0
				!Vr(:,Ns+i:Ns+Ns)=0d0
				!exit
			!endif
		!enddo
		!Vr(:Ns,Ns+1:)=conjg(Vr(Ns+1:,:Ns))
		!Vr(Ns+1:,Ns+1:)=conjg(Vr(:Ns,:Ns))
		!L(:Ns,:)=Vr(:Ns,:)
		!L(Ns+1:,:)=-Vr(Ns+1:,:)
		!L=matmul(conjg(transpose(Vr)),L)
		!call collect(real(W),c)
		!do i=1,size(c)-1
			!if(c(i)>Ns) exit
			!rg=[c(i),c(i+1)-1]
			!if(abs(real(W(rg(1))))<1d-6) then
				!E(rg(1):rg(2))=1d0
				!E(Ns+rg(1):Ns+rg(2))=-1d0
				!L(:,rg(1):rg(2))=0d0
				!L(:,Ns+rg(1):Ns+rg(2))=0d0
				!do j=rg(1),rg(2)
					!L(j,j)=1d0
					!L(Ns+j,Ns+j)=1d0
				!enddo
				!cycle
			!endif
			!if(rg(1)==rg(2)) then
				!E(rg(1))=real(L(rg(1),rg(1)))
				!L(rg(1),rg(1))=1d0
			!else
				!call heev(L(rg(1):rg(2),rg(1):rg(2)),E(rg(1):rg(2)),"v")
			!endif
			!do j=rg(1),rg(2)
				!if(E(j)<0d0) then
					!L(rg(1):rg(2),Ns+j)=L(rg(1):rg(2),j)/sqrt(abs(E(j)))
					!L(Ns+rg(1):Ns+rg(2),j)=conjg(L(rg(1):rg(2),Ns+j))
					!L(Ns+rg(1):Ns+rg(2),Ns+j)=0d0
					!L(rg(1):rg(2),j)=0d0
					!E(Ns+j)=E(j)
					!E(j)=-E(Ns+j)
				!else
					!L(rg(1):rg(2),j)=L(rg(1):rg(2),j)/sqrt(abs(E(j)))
					!L(Ns+rg(1):Ns+rg(2),Ns+j)=conjg(L(rg(1):rg(2),j))
					!E(Ns+j)=-E(j)
				!endif
			!enddo
		!enddo
		!if(present(info)) then
			!Z=H
			!do i=1,size(flag,1)
				!do j=1,size(flag,2)
					!if(i==j) then
						!flag(i,j)=.true.
					!else
						!flag(i,j)=.false.
					!endif
				!enddo
			!enddo
		!endif
		!H=matmul(Vr,L)
		!E=sign(1d0,E)*real(W)
		!if(present(info)) then
			!if(sum(abs(matmul(transpose(conjg(H)),matmul(Z,H))-diag(E)))>1d-6.or.sum(abs(real(matmul(transpose(conjg(H)),matmul(diag([(merge(1d0,-1d0,i<=ham%Hi/2),i=1,ham%Hi)]),H)))),.not.flag)>1d-6) then
				!write(*,*)"diag err"
				!stop
			!endif
		!endif
	!end subroutine
	subroutine init_sw()
		complex(8) :: H(4,4),U(size(H,1),size(H,2))
		real(8) :: E(size(H,1)),Ui(3,3),Uj(3,3),DJ(2)=[1d0,0.1d0]
		real(8), allocatable :: V(:)
		integer, allocatable :: stp(:)
		integer :: Ns,i,j,n,m,Vnb(2),ln=4,tmp(100)
		logical :: flag(size(H,1),size(H,2))
		j=0
		n=0
		do
			!call random_number(i,3)
			!i=i+2
			call random_number(i,2)
			i=1+i*2
			!i=4
			j=j+i
			n=n+1
			tmp(n)=i
			if(n/=1) then
				if((mod(abs(tmp(n-1))-1,2)==1.and.tmp(n-1)>0).or.(mod(abs(tmp(n-1))-1,2)==0.and.tmp(n-1)<0)) then
					tmp(n)=-tmp(n)
				endif
			endif
			if(j>50.and.((mod(abs(tmp(n))-1,2)==0.and.tmp(n)>0).or.(mod(abs(tmp(n))-1,2)==1.and.tmp(n)<0)).and.sum(mod(abs(tmp(:n))-1,2)*sign(1,tmp(:n)))==0) then
				exit
			endif
		enddo
		allocate(stp(n))
		stp=tmp(:n)
		write(*,"(*(i3))")stp,sum(abs(stp))

		allocate(ham%var(-100:100))
		open(101,file="../data/lattice.dat")
		! lattice 
		!latt%is_all=.true.
		latt%a1=[1d0,0d0,0d0]
		latt%a2=[0d0,1d0,0d0]
		latt%c1=latt%a1
		latt%c2=latt%a2
		!latt%c1=[-1d0,1d0,0d0]
		!latt%c2=[1d0,1d0,0d0]
		!latt%c1=[2d0,0d0,0d0]
		!latt%c2=[0d0,1d0,0d0]
		!latt%c1=[4d0,1d0,0d0]
		!latt%c2=[0d0,2d0,0d0]
		!latt%c1=[2d0*ln,0d0,0d0]
		latt%T1=latt%a1*real(sum(abs(stp)))
		latt%T2=latt%a2*32d0
		latt%c1=latt%T1
		latt%c2=[0d0,2d0,0d0]

		latt%bdc(1)=1d0
		latt%bdc(2)=1d0
		allocate(latt%rsb(1,3))
		latt%rsb(1,:)=[0d0,0d0,0d0]
		call latt%gen_latt(3)
		write(*,*)"lattice size is ",latt%Ns,"Nc: ",latt%Nc
		call latt%gen_brizon(brizon)
		write(*,*)"brizon%k: ",size(brizon%k,1),"nk: ",brizon%nk,"nq: ",brizon%nq
		call check_lattice(101)

		allocate(var(latt%Ni),fi(latt%Ni),th(latt%Ni))
		do i=1,latt%Ni
			n=0
			do j=1,size(stp)
				if(n<=nint(latt%nb(0)%bd(i)%r(1)).and.(n+abs(stp(j)))>nint(latt%nb(0)%bd(i)%r(1))) then
					n=n-nint(latt%nb(0)%bd(i)%r(1))
					if(n==0) then
						var(i)=0d0
						fi(i)=0d0
						th(i)=0d0
					else
						var(i)=0.5d0
						fi(i)=0d0
						if(mod(nint(latt%nb(0)%bd(i)%r(2))+n+(1-sign(1,stp(j)))/2,2)==0) then
							th(i)=0d0
						else
							th(i)=pi
						endif
					endif
					exit
				endif
				n=n+abs(stp(j))
			enddo

		enddo
		do i=1,latt%Ni
			write(30,"(*(es12.4))")latt%nb(0)%bd(i)%r,var(i),fi(i),th(i)
		enddo
		write(30,"(x/)")

		Vnb=[1,3]

		! sdw
		do inb=1,size(Vnb)
			if(allocated(V)) deallocate(V)
			allocate(V(size(latt%nb(Vnb(inb))%bd)))
			do n=1,size(latt%nb(Vnb(inb))%bd)
				if(Vnb(inb)==1) then
					V(n)=-DJ(inb)*sign(1d0,cos(th(latt%nb(Vnb(inb))%bd(n)%i(1)))*cos(th(latt%nb(Vnb(inb))%bd(n)%i(2))))
					!V(n)=DJ(inb)
				elseif(Vnb(inb)==3) then
					m=latt%nb(Vnb(inb))%bd(n)%i(1)
					V(n)=0d0
					do i=1,size(latt%nb(1)%st(m)%j)
						if((sum(abs(latt%nb(1)%bd(latt%nb(1)%st(m)%bd(i))%dr*latt%nb(1)%st(m)%dir(i)-latt%nb(Vnb(inb))%bd(n)%dr/2d0))+abs(var(latt%nb(1)%st(m)%j(i))))<1d-6) then
							if(abs(latt%nb(Vnb(inb))%bd(n)%dr(2))<1d-6) then
								V(n)=-DJ(inb)*sign(1d0,cos(th(latt%nb(Vnb(inb))%bd(n)%i(1)))*cos(th(latt%nb(Vnb(inb))%bd(n)%i(2))))
								exit
							endif
						endif
					enddo
				endif
				write(30,"(*(es12.4))")latt%nb(Vnb(inb))%bd(n)%r,latt%nb(Vnb(inb))%bd(n)%dr,V(n)
			enddo
			write(30,"(x/)")
			idx=ham%gen_var(dir=[0,0],left=[+1,-2],right=[-1,+2])
			do i=1,size(ham%var(idx)%bd)
				Ui(1,:)=[cos(th(i))*cos(fi(i)),-sin(fi(i)),sin(th(i))*cos(fi(i))]
				Ui(2,:)=[cos(th(i))*sin(fi(i)),cos(fi(i)),sin(th(i))*sin(fi(i))]
				Ui(3,:)=[-sin(th(i)),0d0,cos(th(i))]
				ham%var(idx)%bd(i)=0d0
				do n=1,size(latt%nb(Vnb(inb))%st(i)%j)
					j=latt%nb(Vnb(inb))%st(i)%j(n)
					Uj(1,:)=[cos(th(j))*cos(fi(j)),-sin(fi(j)),sin(th(j))*cos(fi(j))]
					Uj(2,:)=[cos(th(j))*sin(fi(j)),cos(fi(j)),sin(th(j))*sin(fi(j))]
					Uj(3,:)=[-sin(th(j)),0d0,cos(th(j))]
					ham%var(idx)%bd(i)=ham%var(idx)%bd(i)+var(j)*sum(Uj(:,3)*(-Ui(:,3)))*V(latt%nb(Vnb(inb))%st(i)%bd(n))*abs(latt%nb(Vnb(inb))%bd(latt%nb(Vnb(inb))%st(i)%bd(n))%bdc)
				enddo
			enddo
			ham%var(idx)%val=1d0

			! sdw
			idx=ham%gen_var(dir=[Vnb(inb),Vnb(inb)],left=[+1,-2],right=[-1,+2],cg=[.true.,.true.])
			do n=1,size(ham%var(idx)%bd)
				i=latt%nb(ham%var(idx)%nb)%bd(n)%i(1)
				Ui(1,:)=[cos(th(i))*cos(fi(i)),-sin(fi(i)),sin(th(i))*cos(fi(i))]
				Ui(2,:)=[cos(th(i))*sin(fi(i)),cos(fi(i)),sin(th(i))*sin(fi(i))]
				Ui(3,:)=[-sin(th(i)),0d0,cos(th(i))]
				j=latt%nb(ham%var(idx)%nb)%bd(n)%i(2)
				Uj(1,:)=[cos(th(j))*cos(fi(j)),-sin(fi(j)),sin(th(j))*cos(fi(j))]
				Uj(2,:)=[cos(th(j))*sin(fi(j)),cos(fi(j)),sin(th(j))*sin(fi(j))]
				Uj(3,:)=[-sin(th(j)),0d0,cos(th(j))]
				ham%var(idx)%bd(n)=V(n)*sum(matmul(sqrt(var(j)/2d0)*[complex(8)::1d0,img],transpose(Uj(:,:2)))*matmul(Ui(:,:2),sqrt(var(i)/2d0)*[complex(8)::1d0,-img]))
			enddo
			ham%var(idx)%val=1d0

			! sdw
			idx=ham%gen_var(dir=[Vnb(inb),-Vnb(inb)],left=[+1,+1],right=[+2,+2],cg=[.true.,.true.])
			do n=1,size(ham%var(idx)%bd)
				i=latt%nb(ham%var(idx)%nb)%bd(n)%i(1)
				Ui(1,:)=[cos(th(i))*cos(fi(i)),-sin(fi(i)),sin(th(i))*cos(fi(i))]
				Ui(2,:)=[cos(th(i))*sin(fi(i)),cos(fi(i)),sin(th(i))*sin(fi(i))]
				Ui(3,:)=[-sin(th(i)),0d0,cos(th(i))]
				j=latt%nb(ham%var(idx)%nb)%bd(n)%i(2)
				Uj(1,:)=[cos(th(j))*cos(fi(j)),-sin(fi(j)),sin(th(j))*cos(fi(j))]
				Uj(2,:)=[cos(th(j))*sin(fi(j)),cos(fi(j)),sin(th(j))*sin(fi(j))]
				Uj(3,:)=[-sin(th(j)),0d0,cos(th(j))]
				ham%var(idx)%bd(n)=V(n)*sum(matmul(sqrt(var(j)/2d0)*[complex(8)::1d0,img],transpose(Uj(:,:2)))*matmul(Ui(:,:2),sqrt(var(i)/2d0)*[complex(8)::1d0,img]))
			enddo
			ham%var(idx)%val=1d0
		enddo


		call ham%init([2])
		mat_diag => diag_B
	end subroutine
	subroutine bandhot_B(self,ut,ki,kf,n,omg,m,gm)
		class(t_ham) :: self
		integer :: ut,n,m
		real(8) :: omg(2),gm,ki(3),kf(3)
		complex(8) :: H(self%Hi,self%Hi),U(size(H,1),size(H,2))
		real(8) :: E(size(H,1)),Ui(3,3),W
		logical :: f,flag(size(H,1),size(H,2))
		real(8) :: k(3),dk(3),domg
		integer :: i,l,l1,l2,j
		complex(8) :: G(size(E),3)
		domg=(omg(2)-omg(1))/m
		dk=(kf-ki)/n
		do l1=0,n-1
			k=ki+dk*l1
			call self%Hamilton([self%rg(1):self%rg(2)],H,k)
			call mat_diag(H,E,1)
			G=0d0
			do l=latt%Ni+1,2*latt%Ni
				do i=1,latt%Ni
					Ui(1,:)=[cos(th(i))*cos(fi(i)),-sin(fi(i)),sin(th(i))*cos(fi(i))]
					Ui(2,:)=[cos(th(i))*sin(fi(i)),cos(fi(i)),sin(th(i))*sin(fi(i))]
					Ui(3,:)=[-sin(th(i)),0d0,cos(th(i))]
					G(l-latt%Ni,:)=G(l-latt%Ni,:)+matmul(Ui(:,:2),sqrt(var(i)/2d0)*[complex(8)::1d0,-img])*H(i,l)+matmul(Ui(:,:2),sqrt(var(i)/2d0)*[complex(8)::1d0,img])*H(latt%Ni+i,l)
				enddo
			enddo
			G(:,1)=sum(abs2(G),2)/pbrizon%nq
			G(:,1)=merge(real(G(:,1)),3d0,abs(real(G(:,1)))<3d0)
			do l2=0,m-1
				write(ut,"((es17.9)$)")k,omg(1)+domg*l2,sum(real(G(:,1))/(omg(1)+domg*l2-E(latt%Ni+1:)+img*gm))
				write(ut,"(x)")
			enddo
			write(ut,"(x)")
		enddo
	end subroutine
	subroutine band_B(self,ut,ki,kf,n)
		class(t_ham) :: self
		integer :: ut,n
		complex(8) :: H(self%Hi,self%Hi),U(size(H,1),size(H,2))
		real(8) :: E(size(H,1)),Ui(3,3)
		logical :: f,flag(size(H,1),size(H,2))
		real(8) :: ki(3),kf(3)
		real(8) :: k(3),dk(3)
		integer :: i,l,j,m
		complex(8) :: G(3)
		dk=(kf-ki)/n
		do m=0,n-1
			k=ki+dk*m
			call self%Hamilton([self%rg(1):self%rg(2)],H,k)
			call mat_diag(H,E,1)
			write(ut,"((es17.9)$)")k
			do l=latt%Ni+1,2*latt%Ni
				G=0d0
				do i=1,latt%Ni
					Ui(1,:)=[cos(th(i))*cos(fi(i)),-sin(fi(i)),sin(th(i))*cos(fi(i))]
					Ui(2,:)=[cos(th(i))*sin(fi(i)),cos(fi(i)),sin(th(i))*sin(fi(i))]
					Ui(3,:)=[-sin(th(i)),0d0,cos(th(i))]
					G=G+matmul(Ui(:,:2),sqrt(var(i)/2d0)*[complex(8)::1d0,-img])*H(i,l)+matmul(Ui(:,:2),sqrt(var(i)/2d0)*[complex(8)::1d0,img])*H(latt%Ni+i,l)
				enddo
				write(ut,"((es17.9)$)")E(l),merge(real(sum(G*conjg(G)))/pbrizon%nq,0d0,abs(E(l))>1d-3)
			enddo
			write(ut,"(x)")
		enddo
	end subroutine
end module
program main
	use sw
	implicit none
	integer :: i
	logical :: f


	f=openfile(unit=20,file='../data/band.dat')
	f=openfile(unit=30,file='../data/order.dat')
	f=openfile(unit=40,file='../data/bandhot.dat')
	call init_sw()

	!call band_B(ham,20,[3d0*pi/2d0,pi/2d0,0d0],[pi,pi,0d0],32*8)
	!call band_B(ham,20,[pi,pi,0d0],[pi,0d0,0d0],32*8)
	!call band_B(ham,20,[pi,0d0,0d0],[3d0*pi/2d0,pi/2d0,0d0],32*8)
	!call band_B(ham,20,[3d0*pi/2d0,pi/2d0,0d0],[2d0*pi,0d0,0d0],32*8)
	!call band_B(ham,20,[2d0*pi,0d0,0d0],[pi,0d0,0d0],32*8)
	!stop

	!call band_B(ham,20,[pi,0d0,0d0],[pi,pi,0d0],32*8)
	!call band_B(ham,20,[pi,pi,0d0],[0d0,pi,0d0],32*8)
	!stop

	!call band_B(ham,20,[0d0,0d0,0d0],[2d0*pi,0d0,0d0],32*8*2)
	!stop

	!call band_B(ham,20,[0d0,0d0,0d0],[pi,pi,0d0],32*8)
	!call band_B(ham,20,[pi,pi,0d0],[0d0,pi,0d0],32*8)
	call bandhot_B(ham,40,[0d0,0d0,0d0],[pi,pi,0d0],32*8,[0d0,1.5d0],60,0.1d0)
	call bandhot_B(ham,40,[pi,pi,0d0],[0d0,pi,0d0],32*8,[0d0,1.5d0],60,0.1d0)
	stop

	call band_B(ham,20,[0d0,pi,0d0],[pi,pi,0d0],32*8)
	call band_B(ham,20,[pi,pi,0d0],[pi,0d0,0d0],32*8)
	call band_B(ham,20,[pi,pi,0d0],[0d0,0d0,0d0],32*8)
	!do i=1,size(flag,1)
		!do j=1,size(flag,2)
			!if(i==j) then
				!flag(i,j)=.true.
			!else
				!flag(i,j)=.false.
			!endif
		!enddo
	!enddo

	!H=0d0
	!Ns=size(H,1)/2

	!!H(1:Ns,1:Ns)=reshape([1d0,2d0,0d0,2d0,1d0,3d0,0d0,3d0,1d0],[Ns,Ns])
	!!H(Ns+1:,Ns+1:)=H(1:Ns,1:Ns)
	!!H(Ns+1:,:Ns)=reshape([1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0],[Ns,Ns])
	!!H(:Ns,Ns+1:)=conjg(H(Ns+1:,:Ns))

	!H(1:Ns,1:Ns)=reshape([1d0,2d0,2d0,1d0],[Ns,Ns])
	!H(Ns+1:,Ns+1:)=H(1:Ns,1:Ns)
	!H(Ns+1:,:Ns)=reshape([1d0,0d0,0d0,1d0],[Ns,Ns])
	!H(:Ns,Ns+1:)=conjg(H(Ns+1:,:Ns))

	!write(*,"("//to_char(Ns*2)//"es9.2)")real(H)

	!U=H
	!call diag_B(U,E)
	!read(*,*)
	!write(*,"(es9.2)")sum(abs(matmul(transpose(conjg(U)),matmul(H,U))-diag(E)))
	!write(*,*)"---------UHU------------"
	!write(*,"(es9.2)")sum(abs(real(matmul(transpose(conjg(U)),matmul(diag([(merge(1d0,-1d0,i<=Ns),i=1,Ns*2)]),U)))),.not.flag)
	!write(*,*)"------------------------"
	!write(*,"("//to_char(Ns*2)//"es9.2)")pack(real(matmul(transpose(conjg(U)),matmul(diag([(merge(1d0,-1d0,i<=Ns),i=1,Ns*2)]),U))),flag)
	!write(*,*)"---------UU-------------"
	!write(*,"("//to_char(Ns*2)//"es9.2)")E
	!write(*,*)"----------E-------------"
	
end program
