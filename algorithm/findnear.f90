module findnear
	implicit none
	save
	integer, parameter :: dn=1000
	integer :: pbc(0:dn+1)
	contains
	subroutine nnp(ii,near,jj)
		implicit none
		integer :: ii,i,j,jj(4)
		logical :: near
		pbc=(/dn,(i,i=1,dn),1/)
		j=pbc(mod(ii,dn))
		i=(ii-j)/dn+1
		if(near) then
			jj(1)=(i-1)*dn+pbc(j+1)
			jj(2)=(pbc(i+1)-1)*dn+j
			jj(3)=(i-1)*dn+pbc(j-1)
			jj(4)=(pbc(i-1)-1)*dn+j
		else
			jj(1)=dn*(pbc(i+1)-1)+pbc(j+1)
			jj(2)=dn*(pbc(i+1)-1)+pbc(j-1)
			jj(3)=dn*(pbc(i-1)-1)+pbc(j+1)
			jj(4)=dn*(pbc(i-1)-1)+pbc(j-1)
		endif
	end subroutine
end module
program main
	use findnear
	implicit none
	integer :: i,nb(4)
	open(unit=10,file="../data/output.dat")
	do i=1,dn*dn
		call nnp(i,.false.,nb)
		write(10,*)nb
	enddo
end
