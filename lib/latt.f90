module M_latt
	use M_utility
	implicit none
contains
	subroutine square_one2two(n,Ns,i)
		integer :: n,i(:),Ns(:)
		i(1)=mod(n-1,Ns(1))
		i(2)=(n-1)/Ns(1)
	end subroutine
	subroutine square_two2one(i,Ns,n)
		integer :: n,i(:),Ns(:)
		n=i(2)*Ns(1)+i(1)+1
	end subroutine
	subroutine square_bc(i,Ns,Tx,Ty,Ns2,n)
		integer :: n,i(:),ti(2),Ns(:),v(2,4,2),Tx(:),Ty(:),Ns2,p
		v=reshape((/&
			1,0, 0,1, -1,0, 0,-1, &
			1,1, -1,1, -1,-1, 1,-1 &
			/),shape(v))
		ti=i
		p=0
		do
			call square_two2one(ti,Ns,n)
			if((ti(1)<0.or.ti(1)>Ns(1)-1.or.ti(2)<0.or.ti(2)>Ns(2)-1.or.n<1.or.n>Ns2)) then
				ti=i+v(1,mod(p,4)+1,p/4+1)*Tx+v(2,mod(p,4)+1,p/4+1)*Ty
				p=p+1
			else
				i=ti
				exit
			endif
		enddo
	end subroutine
	subroutine square(Ns,Tx,Ty,neb)
		!						i=0...Ns(1)-1
		!					 1 --  2 --  3 --  4
		!					 |     |     |     |           4      3   4
		!					 5 --  6 --  7 --  8           |       \ /
		!	j=0...Ns(2)-1	 |     |     |     |        3--0--1     0
		!					 9 -- 10 -- 11 -- 12           |       / \
		!					 |     |     |     |           2      2   1
		!					13 -- 14 -- 15 -- 16
		implicit none
		integer :: n,m,k,i(2),j(2),Ns(:),neb(:,:,:),v(2,4,size(neb,3)),Tx(:),Ty(:),Ns2
		v=reshape((/&
			1,0, 0,1, -1,0, 0,-1, &
			1,1, -1,1, -1,-1, 1,-1, &
			2,0, 0,2, -2,0, 0,-2 &
			/),shape(v))
		Ns2=size(neb,1)
		do n=1,Ns(1)*Ns(2)
			call square_one2two(n,Ns,i)
			do m=1,size(neb,3)
				do k=1,4
					j=i+v(:,k,m)
					call square_bc(j,Ns,Tx,Ty,Ns2,neb(n,k,m))
					!call square_two2one((/mod(j(1)+Ns(1),Ns(1)),mod(j(2)+Ns(2),Ns(2))/),Ns,neb(n,k,m))
				enddo
			enddo
		enddo
	end subroutine
	subroutine square_tilda_one2two(n,Ns,i)
		integer :: n,i(:),Ns(:)
		i(1)=mod(n+Ns(1)/2-1,Ns(1))
		i(2)=(n+Ns(1)/2-1)/Ns(1)
	end subroutine
	subroutine square_tilda_two2one(i,Ns,n)
		integer :: n,i(:),Ns(:)
		n=i(2)*Ns(1)+i(1)+1-Ns(1)/2
	end subroutine
	subroutine square_tilda_bc(i,Ns,Tx,Ty,Ns2,n)
		integer :: n,i(:),ti(2),Ns(:),v(2,4,2),Tx(:),Ty(:),Ns2,p
		v=reshape((/&
			1,0, 0,1, -1,0, 0,-1, &
			1,1, -1,1, -1,-1, 1,-1 &
			/),shape(v))
		ti=i
		p=0
		do
			call square_tilda_two2one(ti,Ns,n)
			if((ti(1)<0.or.ti(1)>Ns(1)-1.or.ti(2)<0.or.ti(2)>Ns(2)-1.or.n<1.or.n>Ns2)) then
				ti=i+v(1,mod(p,4)+1,p/4+1)*Tx+v(2,mod(p,4)+1,p/4+1)*Ty
				p=p+1
			else
				i=ti
				exit
			endif
		enddo
	end subroutine
	subroutine square_tilda(Ns,Tx,Ty,neb)
		implicit none
		integer :: n,m,k,i(2),j(2),Tx(:),Ty(:),neb(:,:,:),v(2,4,size(neb,3)),p,tj(2),Ns(:),Ns2
		v=reshape((/&
			1,0, 0,1, -1,0, 0,-1, &
			1,1, -1,1, -1,-1, 1,-1, &
			2,0, 0,2, -2,0, 0,-2 &
			/),shape(v))
		Ns2=size(neb,1)
		do n=1,Ns2
			call square_tilda_one2two(n,Ns,i)
			do m=1,size(neb,3)
				do k=1,4
					j=i+v(:,k,m)
					call square_tilda_bc(j,Ns,Tx,Ty,Ns2,neb(n,k,m))
				enddo
			enddo
		enddo
	end subroutine
end module
