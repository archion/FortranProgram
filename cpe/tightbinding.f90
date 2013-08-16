MODULE GLOBAL
	IMPLICIT NONE
	REAL(8),SAVE :: T=1
	INTEGER, PARAMETER :: DN=120
	REAL(8), PARAMETER :: PI=3.14159
	INTEGER,SAVE :: nn(DN,DN),PBC(0:DN+1)
END MODULE
PROGRAM TIGHTBINDING
	USE GLOBAL
	IMPLICIT NONE
	COMPLEX(8), PARAMETER :: II=(0,1)
	COMPLEX(8) :: G(DN*DN)
	REAL(8) :: E(DN*DN),sg=0.01,H(DN*DN,DN*DN),WORK(3*DN*DN),FQ
	INTEGER :: i,j,k,NB(4),INFO
	OPEN(UNIT=10,FILE='../DATA/output.dat')
	PBC=(/DN,(i,i=1,DN),1/)
	k=1
	DO i=1,DN
		DO j=1,DN
			nn(i,j)=k
			k=k+1
		ENDDO
	ENDDO
	H=0
	DO i=1,DN*DN
		CALL NNP(i,.TRUE.,NB)
		DO j=1,4
			H(i,NB(j))=T
		ENDDO
	ENDDO
	CALL DSYEV( "N", "U", DN*DN, H, DN*DN, E, WORK, 3*DN*DN, INFO )
	DO i=1,1100
		FQ=MINVAL(E)+i*ABS(MINVAL(E))/1000D0
		G=1.0/(FQ-E-II*sg)
		WRITE(10,"(2E15.6)")FQ,DIMAG(SUM(G))
	ENDDO
END
SUBROUTINE NNP(ii,near,jj)
	USE GLOBAL
	IMPLICIT NONE
	INTEGER :: ii,i,j,jj(4)
	LOGICAL :: near
	j=PBC(MOD(ii,DN))
	i=(ii-j)/DN+1
	IF(near) THEN
		jj(1)=nn(i,PBC(j+1))
		jj(2)=nn(PBC(i+1),j)
		jj(3)=nn(i,PBC(j-1))
		jj(4)=nn(PBC(i-1),j)
	ELSE
		jj(1)=nn(PBC(i+1),PBC(j+1))
		jj(2)=nn(PBC(i+1),PBC(j-1))
		jj(3)=nn(PBC(i-1),PBC(j+1))
		jj(4)=nn(PBC(i-1),PBC(j-1))
	ENDIF
END SUBROUTINE NNP