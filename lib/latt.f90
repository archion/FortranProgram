module M_latt
	implicit none
contains
	subroutine square_one2two(n,Ns,i)
		integer :: n,i(:),Ns(:)
		i(1)=(n-1)/Ns(1)
		i(2)=mod(n-1,Ns(1))
	end subroutine
	subroutine square_two2one(i,Ns,n)
		integer :: n,i(:),Ns(:)
		n=i(1)*Ns(1)+i(2)+1
	end subroutine
	subroutine square(Ns,neb)
		!					j=0...
		!			 1 --  2 --  3 --  4
		!			 |     |     |     |           4      3   4
		!			 5 --  6 --  7 --  8           |       \ /
		!	i=0...	 |     |     |     |        3--0--1     0
		!			 9 -- 10 -- 11 -- 12           |       / \
		!			 |     |     |     |           2      2   1
		!			13 -- 14 -- 15 -- 16
		implicit none
		integer :: n,m,k,i(2),j(2),Ns(:),neb(:,:,:),v(2,4,size(neb,3))
		v=reshape((/&
			0,1, 1,0, 0,-1, -1,0, &
			1,1, 1,-1, -1,-1, -1,1, &
			0,2, 2,0, 0,-2, -2,0 &
			/),shape(v))
		do n=1,Ns(1)*Ns(2)
			call square_one2two(n,Ns,i)
			do m=1,size(neb,3)
				do k=1,4
					j=i+v(:,k,m)
					call square_two2one((/mod(j(1)+Ns(2),Ns(2)),mod(j(2)+Ns(1),Ns(1))/),Ns,neb(n,k,m))
				enddo
			enddo
		enddo
	end subroutine
end module
