MODULE FINDNEAR
	IMPLICIT NONE
	SAVE
	INTEGER, PARAMETER :: DN=1000
	INTEGER :: PBC(0:DN+1)
	CONTAINS
	SUBROUTINE NNP(ii,near,jj)
		IMPLICIT NONE
		INTEGER :: ii,i,j,jj(4)
		LOGICAL :: near
		PBC=(/DN,(i,i=1,DN),1/)
		j=PBC(MOD(ii,DN))
		i=(ii-j)/DN+1
		IF(near) THEN
			jj(1)=(i-1)*DN+PBC(j+1)
			jj(2)=(PBC(i+1)-1)*DN+j
			jj(3)=(i-1)*DN+PBC(j-1)
			jj(4)=(PBC(i-1)-1)*DN+j
		ELSE
			jj(1)=DN*(PBC(i+1)-1)+PBC(j+1)
			jj(2)=DN*(PBC(i+1)-1)+PBC(j-1)
			jj(3)=DN*(PBC(i-1)-1)+PBC(j+1)
			jj(4)=DN*(PBC(i-1)-1)+PBC(j-1)
		ENDIF
	END SUBROUTINE
END MODULE
PROGRAM MAIN
	USE FINDNEAR
	IMPLICIT NONE
	INTEGER :: i,NB(4)
	OPEN(UNIT=10,FILE="../DATA/OUTPUT.DAT")
	DO i=1,DN*DN
		CALL NNP(i,.FALSE.,NB)
		WRITE(10,*)NB
	ENDDO
END