module global
	implicit none
	integer, parameter :: Ns(2)=(/36,36/),Ns2=Ns(1)*Ns(2),ne=16,ne2=ne/2
	integer :: latt(Ns2,4,3),pbcx(-Ns(1)+1:2*Ns(1)),pbcy(-Ns(2)+1:2*Ns(2))
	real(8), parameter :: t(3)=(/4d0,-2d0,2d0/),J=0.5d0
	complex(8), parameter :: img=(0d0,1d0)
	contains
		subroutine gen_latt_square()
!					  j
!			 1 --  2 --  3 --  4
!			 |     |     |     |           4      3   4
!			 5 --  6 --  7 --  8           |       \ /
!		i	 |     |     |     |        3--0--1     0
!			 9 -- 10 -- 11 -- 12           |       / \
!			 |     |     |     |           2      2   1
!			13 -- 14 -- 15 -- 16
			implicit none
			integer :: ii,i,j
			pbcx=(/(i,i=1,Ns(1)),(i,i=1,Ns(1)),(i,i=1,Ns(1))/)
			pbcy=(/(i,i=1,Ns(2)),(i,i=1,Ns(2)),(i,i=1,Ns(2))/)
			do ii=1,Ns2
				j=mod(ii-1,Ns(1))+1
				i=(ii-1)/Ns(1)+1
				latt(ii,1,1)=Ns(1)*(i-1)+pbcx(j+1)
				latt(ii,2,1)=Ns(1)*(pbcy(i+1)-1)+j
				latt(ii,3,1)=Ns(1)*(i-1)+pbcx(j-1)
				latt(ii,4,1)=Ns(1)*(pbcy(i-1)-1)+j
				latt(ii,1,2)=Ns(1)*(pbcy(i+1)-1)+pbcx(j+1)
				latt(ii,2,2)=Ns(1)*(pbcy(i+1)-1)+pbcx(j-1)
				latt(ii,3,2)=Ns(1)*(pbcy(i-1)-1)+pbcx(j-1)
				latt(ii,4,2)=Ns(1)*(pbcy(i-1)-1)+pbcx(j+1)
				latt(ii,1,3)=Ns(1)*(i-1)+pbcx(j+2)
				latt(ii,2,3)=Ns(1)*(pbcy(i+2)-1)+j
				latt(ii,3,3)=Ns(1)*(i-1)+pbcx(j-2)
				latt(ii,4,3)=Ns(1)*(pbcy(i-2)-1)+j
			enddo
		end subroutine gen_latt_square
end module
module vmc
	use global
	use rn
	implicit none
	contains
		subroutine ini(pg,sc,wf)
			use fft
			do i=1,Ns
				do j=1,Ns
					alp(i,j)=
				enddo
			enddo
			call fft2d(alp)
			ep(:,1)=tmp(1:n/2)
			ep(:,2)=tmp(n/2+1:n)
			hp(:)=tmp(n+1:Ns)
			do i=1,n/2
				do j=1,n/2
					psi(i,j)=alp(ep(i,1)-ep(j,2))
				enddo
			enddo
		end subroutine
		subroutine mc(wf)
			call fisher_yates_shuffle(cfg,Ns)
			icfg=0
			do i=1,ne2
				icfg(cfg(i),1)=1
				icfg(cfg(i),2)=i
				icfg(cfg(i+ne2),1)=-1
				icfg(cfg(i+ne2),2)=i+ne2
			enddo
			do i=1,ne2
				do j=1,ne2
					iA(i,j)=wf(cfg(i)-cfg(ne2+j))
				enddo
			enddo
			call matrix_inv(iA)
			call phy(psinv,E)
			do
				U=0d0
				V=0d0
				call random_number(i,ne)
				call random_number(j,Ns-ne2)
				cfgp=cfg
				if(j>ne2) then
					sg=i/(ne2+1)+1
					cr=(/(2-sg),sg-1/)*i
					j=ne+j
					cfgp(i)=cfg(j)
					cfgp(j)=cfg(i)
					do k=1,ne/2
						uv(k,sg)=wf(cfg(i)-cfg(ne2*(2-sg)+k))-wf(cfgp(i)-cfg(ne2*(2-sg)+k))
					enddo
				else
					sg=3
					cr=(/i,j/)
					cfgp(i)=cfg(j)
					cfgp(j)=cfg(i)
					j=j+i/(ne2+1)*ne2
					do k=1,ne/2
					enddo
				endif
				call det(uv,cr,iA,dv,pb,sg)
				call random_number(tmp)
				if(tmp<pb*dconjg(pb)) then
					icfg(cfg(i),1)=1
					icfg(cfg(i),2)=i
					icfg(cfg(i+ne2),1)=-1
					icfg(cfg(i+ne2),2)=i+ne2
					call invupdate(U,V,cr,iA,dv,sg)
					if(mod(Nmc,100)==0) then
						do i=1,ne2
							do j=1,ne2
								iA(i,j)=wf(cfg(i)-cfg(ne2+j))
							enddo
						enddo
						call matrix_inv(iA)
					endif
					call phy(psi,psinv,E)
				endif
				Eavg=Eavg+E
				if(Nmc>Nmax) then
					exit
				endif
				Nmc=Nmc+1
			enddo
		end subroutine
		subroutine det(uv,cr,iA,dv,pb,sg)
			complex(8) :: uv(:,:),iA(:,:),dv(:,:),pb
			integer :: cr(:),sg,i
			dv=0d0
			if(sg==3) then
				dv(1,1)=dot_product(iA(cr(1),:),uv(:,1))+1d0
				dv(2,2)=dot_product(iA(:,cr(2)),uv(:,2))+1d0
				dv(1,2)=A(cr(1),cr(2))
				do i=1,ne2
					dv(2,1)=dv(2,1)+dot_product(iA(:,i),uv(:,1))*uv(i,2)
				enddo
			else
				if(sg==1) then
					dv(1,1)=dot_product(iA(cr(1),:),uv(:,1))+1d0
					dv(2,2)=1d0
				else
					dv(1,1)=1d0
					dv(2,2)=dot_product(iA(:,cr(2)),uv(:,2))+1d0
				endif
			endif
			pb=dv(1,1)*dv(2,2)-dv(1,2)*dv(2,1)
		end subroutine
		subroutine invupdate(uv,cr,,dv,sg)
			complex(8) :: uv(:,:),iA(:,:),dv(:,:),pb
			integer :: cr(:),n,s,i,j
			n=size(iA,1)
			if(sg==3) then
				do i=1,n
					uv(i,1)=dot_product(iA(i,:),uv(:,1))
					uv(i,2)=dot_product(iA(:,i),uv(:,2))
				enddo
				do i=1,n
					do j=1,n
						iA(i,j)=iA(i,j)-1d0/pb*((U(i,1)*dv(2,2)-iA(i,cr(1))*dv(2,1))*iA(cr(2),j)+(-uv(i,1)*dv(1,2)+iA(i,cr(1))*dv(1,1))*uv(j,2))
					enddo
				enddo
			else
				do i=1,n
					do j=1,n
						if(sg==1) then
							iA(i,j)=iA(i,j)-1d0/pb*dot_product(iA(i,:),uv(:,1))*iA(cr(1),j)
						else
							iA(i,j)=iA(i,j)-1d0/pb*dot_product(iA(:,j),uv(:,2))*iA(i,cr(2))
						endif
					enddo
				enddo
			endif
		end subroutine
		subroutine jastrow()
		end subroutine
		subroutine energy(cfg)
			E=0d0
			do i=1,ne
				do j=1,4
					do k=1,size(t)
						!hopping
						if(icfg(latt(cfg(i),j,k))==0) then
							sg=i/(ne2+1)+1
							cr=(/(2-sg),sg-1/)*i
							do k=1,ne/2
								uv(k,sg)=wf(cfg(i)-cfg(ne2*(2-sg)+k))-wf(latt(cfg(i),j,1)-cfg(ne2*(2-sg)+k))
							enddo
							call det(uv,cr,iA,dv,pb,sg)
							E=E+t(k)*pb
						endif
					enddo
					!spin flip
					if(icfg(latt(cfg(i),j,1))*icfg(cfg(i))==-1) then
						cr=(/i,j/)
						cfgp(i)=cfg(j)
						cfgp(j)=cfg(i)
						j=j+i/(ne2+1)*ne2
						do k=1,ne/2
						enddo
						call det(uv,cr,iA,dv,sg)
						E=E-DJ*pb
					else
						!diagnal
						E=E+DJ*icfg(latt(cfg(i),j,1))*icfg(cfg(i))
					endif
				enddo
			enddo
		end subroutine
		subroutine matrix_inv(A)
			use lapack95
			implicit none
			complex(8) :: A(:,:)
			integer :: ipiv(size(A,1)),info
			call getrf(A,ipiv,info)
			if(info/=0) then
				write(*,*)"error1",info
				stop
			endif
			call getri(A(i,j,k,:,:),ipiv,info)
			if(info/=0) then
				write(*,*)"error2",info
				stop
			endif
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
program main
	call
end
