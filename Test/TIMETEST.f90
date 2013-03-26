PROGRAM TIMETEST
	IMPLICIT NONE
	CHARACTER(24) :: CDATE
	INTEGER(8) :: i,n=0,J,A(100000),DN2=16,k
	LOGICAL :: FLAG
	OPEN(UNIT=10,FILE='../DATA/output.dat')
	A(1)=255
	i=1
EX:	DO WHILE(.TRUE.)
		WRITE(10,"(B17.16)")A(i)
		A(i+1)=A(i)
		i=i+1
		DO j=0,DN2-1
			IF(j==DN2-1) THEN
				EXIT EX
			ENDIF
			IF(IBITS(A(i),j,2)==1) THEN
				A(i)=IBSET(A(i),j+1)
				A(i)=IBCLR(A(i),j)
				DO k=0,j/2-1
					IF(BTEST(A(i),k)) THEN
						EXIT
					ENDIF
					IF(BTEST(A(i),j-k-1)) THEN
						A(i)=IBSET(A(i),k)
						A(i)=IBCLR(A(i),j-k-1)
					ENDIF
				ENDDO
				EXIT
			ENDIF
		ENDDO
	ENDDO EX
	WRITE(*,*)i,A(i),A(i-1),A(i-2)
END
	