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
		subroutine v(pg,sc)
			use fft
			do i=1,Ns
				do j=1,Ns
					alp(i,j)=
				enddo
			enddo
			call fft2d(alp)
			call fisher_yates_shuffle(tmp,Ns)
			ep(:,1)=tmp(1:n/2)
			ep(:,2)=tmp(n/2+1:n)
			hp(:)=tmp(n+1:Ns)
			do i=1,n/2
				do j=1,n/2
					psi(i,j)=alp(ep(i,1)-ep(j,2))
				enddo
			enddo
			call mc(psi,E)
		end subroutine
		subroutine mc(pg,sc)
			psinv=psi
			call matrix_inv(psinv)
			call phy(psinv,E)
			do
				call irandom(ri,ne)
				call irandom(ni,Ns-ne)
				npsi=psi
				if(ri/ne2==0) then
					ep(mod(ri,ne2),1)=hp(ni)
					do i=1,ne2
						npsi(mod(ri,ne2),i)=alp(ep(mod(ri,ne2),1)-ep(i,2))
					enddo
				else
					ep(mod(ri,ne2),2)=hp(ni)
					do i=1,ne2
						npsi(i,mod(ri,ne2))=alp(ep(mod(ri,ne2),1)-ep(i,2))
					enddo
				endif
				call det(npsi,psinv,pb)
				call random_number(tmp)
				if(tmp<pb) then
					call update_inv(psinv,pb)
					if(Nck>100) then
						Nck=0
						psinv=psi
						call matrix_inv(psinv)
					endif
					call phy(psi,psinv,E)
				endif
				Eavg=Eavg+E
				if(Nmc>Nmax) then
					exit
				endif
				Nmc=Nmc+1
				Nck=Nck+1
			enddo
		end subroutine
		subroutine det(npsi,psinv,pb)
			pb=0d0
			do i=1,ne2
				pb=pb+npsi(i,l)*psinv(i,l)
			enddo
		end subroutine
		subroutine update_inv(psinv,pb)

			psinv=psinv-1d0/pb*tmp
		end subroutine
		subroutine energy(wf,j)
			call det

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
