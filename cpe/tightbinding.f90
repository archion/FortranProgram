module global
	implicit none
	real(8),save :: t=1
	integer, parameter :: dn=120
	real(8), parameter :: pi=3.14159
	integer,save :: nn(dn,dn),pbc(0:dn+1)
end module
program tightbinding
	use global
	implicit none
	complex(8), parameter :: ii=(0,1)
	complex(8) :: g(dn*dn)
	real(8) :: e(dn*dn),sg=0.01,h(dn*dn,dn*dn),work(3*dn*dn),fq
	integer :: i,j,k,nb(4),info
	open(unit=10,file='../data/output.dat')
	pbc=(/dn,(i,i=1,dn),1/)
	k=1
	do i=1,dn
		do j=1,dn
			nn(i,j)=k
			k=k+1
		enddo
	enddo
	h=0
	do i=1,dn*dn
		call nnp(i,.true.,nb)
		do j=1,4
			h(i,nb(j))=t
		enddo
	enddo
	call dsyev( "n", "u", dn*dn, h, dn*dn, e, work, 3*dn*dn, info )
	do i=1,1100
		fq=minval(e)+i*abs(minval(e))/1000d0
		g=1.0/(fq-e-ii*sg)
		write(10,"(2e15.6)")fq,dimag(sum(g))
	enddo
end
subroutine nnp(ii,near,jj)
	use global
	implicit none
	integer :: ii,i,j,jj(4)
	logical :: near
	j=pbc(mod(ii,dn))
	i=(ii-j)/dn+1
	if(near) then
		jj(1)=nn(i,pbc(j+1))
		jj(2)=nn(pbc(i+1),j)
		jj(3)=nn(i,pbc(j-1))
		jj(4)=nn(pbc(i-1),j)
	else
		jj(1)=nn(pbc(i+1),pbc(j+1))
		jj(2)=nn(pbc(i+1),pbc(j-1))
		jj(3)=nn(pbc(i-1),pbc(j+1))
		jj(4)=nn(pbc(i-1),pbc(j-1))
	endif
end subroutine nnp
