PROGRAM RPA
	IMPLICIT NONE
	REAL(8), PARAMETER :: PI=3.1415926,t=1,TMP=0.1
	INTEGER, PARAMETER :: N=24
	COMPLEX(8), PARAMETER :: II=(0,1)
	INTEGER :: i,j,k,l,m
	REAL(8) :: U=8,q(2),p(2),fq,z=0.01,E,f,kx,ky,sa,sb,sp
	COMPLEX(8) :: SX0,SX
	f(U)=1/(1/(1+EXP(U/T)))
	E(kx,ky,sp)=-2*t*(COS(kx)+COS(ky))-sp
	OPEN(UNIT=10,FILE="../DATA/OUTPUT.DAT")
	DO i=-N/2,N/2
		q(1)=2*PI/N*i
		DO j=-N/2,N/2
			q(2)=2*PI/N*j
			DO k=1,100
				fq=k*0.003
				SX0=(0,0)
				DO l=-N/2,N/2
					p(1)=2*PI/N*l
					DO m=-N/2,N/2
						p(2)=2*PI/N*m
						SX0=SX0+(f(E(p(1),p(2),sp))-f(E(p(1)+q(1),p(2)+q(2),sp)))/(fq+E(p(1),p(2),sp)-
							&E(p(1)+q(1),p(2)+q(2),sp)+II*z)
					ENDDO
				ENDDO
				SX=SX0/(1-U*SX0)
				WRITE(10,"(5E15.6)")q(1),q(2),fq,SX
			ENDDO
		ENDDO
	ENDDO
	CLOSE(10)
END