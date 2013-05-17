PROGRAM IO
	IMPLICIT NONE
	REAL :: A=0
	INTEGER :: i,j,mp1
	CHARACTER :: S="AB"
	CHARACTER(50) :: FMTC(2)
	OPEN(UNIT=10,FILE="../DATA/IO.DAT")
	mp1=5
	do i=0,mp1
		do j=0,min(i,mp1-i)
			write(*,*)i,j
		enddo
	enddo
END