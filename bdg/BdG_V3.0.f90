MODULE GLOBAL
	IMPLICIT NONE
	REAL(8),SAVE :: BT=1E5,DV=1,nf0=0.85,DU=2.44,T=1,TP=-0.25,eff=0.3
	INTEGER, PARAMETER :: DN=24
	INTEGER,SAVE :: nn(DN,DN),PBC(0:DN+1)
END MODULE
PROGRAM MAIN
	USE GLOBAL
	IMPLICIT NONE
	CHARACTER(24) :: IT,FT
	INTEGER :: i,j,k,n,Im(5)=-1,INFO=0,NP(4),ERROR=0
	REAL(8) :: n1(DN*DN),n2(DN*DN),n1p(DN*DN),n2p(DN*DN),E(2*DN*DN),H(2*DN*DN,2*DN*DN),DT(DN*DN,DN*DN),&
		DT1(DN*DN,DN*DN),DTF(DN*DN),WORK(3*2*DN*DN),RDOM(DN*DN)
	REAL(8) :: f,En,sa,sb,sp,nf,er=1
	f(En)=1/(1+EXP(BT*En))
	OPEN(UNIT=10,FILE='../DATA/output.dat')
	OPEN(UNIT=20,FILE='../DATA/ERROR.dat')
200	FORMAT(24(E15.6))
300	FORMAT(625(E15.6))
	CALL FDATE(IT)
	k=1
	PBC=(/DN,(i,i=1,DN),1/)
	DO i=1,DN
		DO j=1,DN
			nn(i,j)=k
			k=k+1
		ENDDO
	ENDDO
	DT=0D0
	DT1=0D0
	CALL INIT_RANDOM_SEED()
	CALL RANDOM_NUMBER(RDOM)
	RDOM=0.5D0	!Uniform initial value
	DO i=1,DN
		CALL NNP(i,0,NP)
		DT1(i,NP)=(/-0.0795D0,0.0795D0,-0.0795D0,0.0795D0/)+(RDOM(4*i-3:4*i)-0.5D0)/10D0
	ENDDO
	n1p=nf0/2+(RDOM-0.5D0)/20D0
	n2p=nf0-n1p
	n1=0
	n2=0
	! Im(1)=121
	! Im(2)=123
	! Im(1)=(DN/2-1)*DN+DN/2
	! Im(2)=Im(1)+1
	! Im(3)=Im(2)+DN
	! Im(4)=Im(1)+DN
	! Im(5)=Im(1)-1
	! DO i=1,DN*DN
		! IF(ANY(Im==i)) THEN
			! n1(i)=1
		! ENDIF
	! ENDDO
	! WRITE(20,200)((n1(DN*i+j),j=1,DN),i=0,DN-1)
EX:	DO WHILE(er>1E-4)
		DT=DT1
		n1=n1p
		n2=n2p
		nf=0
		sa=0
		sb=0.5
		DO WHILE(ABS(nf-nf0)>1E-5)
			sp=0.5*(sa+sb)
			H=0D0
			DO i=1,DN*DN
				! nearest neighbor pairs
				CALL NNP(i,0,NP)
				DO j=1,4
					H(2*i-1:2*i,2*NP(j)-1:2*NP(j))=H(2*i-1:2*i,2*NP(j)-1:2*NP(j))+RESHAPE((/-T,DT(i,NP(j)),DT(i,NP(j)),&
						T/),(/2,2/))
				ENDDO
				! next-nearest neighbor pairs
				CALL NNP(i,1,NP)
				DO j=1,4
					H(2*i-1:2*i,2*NP(j)-1:2*NP(j))=H(2*i-1:2*i,2*NP(j)-1:2*NP(j))+RESHAPE((/-TP,0.0D0,0.0D0,TP/),(/2,2/))
				ENDDO
				! on site
				H(2*i-1:2*i,2*i-1:2*i)=H(2*i-1:2*i,2*i-1:2*i)+RESHAPE((/DU*n2(i)-sp,0.0D0,0.0D0,-DU*n1(i)+sp/),&
					(/2,2/))
				IF(ANY(Im==i)) THEN
					H(2*i-1:2*i,2*i-1:2*i)=H(2*i-1:2*i,2*i-1:2*i)+RESHAPE((/(-1)**MINLOC(ABS(Im-i))*eff,0.0D0,&
						0.0D0,(-1)**MINLOC(ABS(Im-i))*eff/),(/2,2/))
				ENDIF
			ENDDO
			CALL DSYEV( "V", "U", 2*DN*DN, H, 2*DN*DN, E, WORK, 3*2*DN*DN, INFO )
			IF(INFO/=0) THEN
				EXIT EX
			ENDIF
			nf=0
			DO i=1,DN*DN
				DO n=1,2*DN*DN
					nf=nf+(H(2*i-1,n)**2*f(E(n))+H(2*i,n)**2*(1-f(E(n))))
				ENDDO
			ENDDO
			nf=1.0/(DN*DN)*nf
			IF (nf<nf0) THEN
				sa=sp
			ELSE
				sb=sp
			ENDIF
		ENDDO
		n1p=0
		n2p=0
		DO i=1,DN*DN
			DO n=1,2*DN*DN
				n1p(i)=n1p(i)+H(2*i-1,n)**2*f(E(n))
				n2p(i)=n2p(i)+H(2*i,n)**2*(1-f(E(n)))
			ENDDO
		ENDDO
		er=SUM(ABS(n1p-n1)+ABS(n2p-n2))/(DN*DN)
		DT1=0
		DO i=1,DN*DN
			CALL NNP(i,0,NP)
			DO j=1,4
				DO n=1,2*DN*DN
					DT1(i,NP(j))=DT1(i,NP(j))+(H(2*i-1,n)*H(2*NP(j),n)+H(2*NP(j)-1,n)*H(2*i,n))*TANH(0.5*BT*E(n))
				ENDDO
			ENDDO
		ENDDO
		DT1=0.25*DV*DT1
		WRITE(*,"(E15.6,F15.12)")er,sa
		REWIND(10)
		WRITE(10,200)((n1(DN*i+j)+n2(DN*i+j),j=1,DN),i=0,DN-1)
		REWIND(20)
		DO j=1,DN*DN
			WRITE(20,500)(DT(i,j),i=1,DN*DN)
		ENDDO
	ENDDO EX
	REWIND(10)
	DO i=1,DN*DN
		CALL NNP(i,0,NP)
		DTF(i)=(DT(i,NP(1))+DT(i,NP(3))-DT(i,NP(2))-DT(i,NP(4)))/4.0
	ENDDO
	CALL FDATE(FT)
	WRITE(10,200)((DTF(DN*i+j),j=1,DN),i=0,DN-1)
	WRITE(10,200)((n1(DN*i+j)+n2(DN*i+j),j=1,DN),i=0,DN-1)
	WRITE(10,200)(((-1)**(i+j)*0.5*(n1(DN*i+j)-n2(DN*i+j)),j=1,DN),i=0,DN-1)
	WRITE(10,*)"The beta is",BT,"V=",DV,"U=",DU,"t'=",TP,"heff=",eff,"nf=",nf0,"N=",DN
	WRITE(10,*)"ERROR=",ERROR,"INFO=",INFO
	WRITE(10,*)"THE TATAL RUNTIME IS FORM ",IT," TO ",FT
	IF(INFO/=0) THEN
		ERROR=MOD(INFO,2*DN*DN+1)
		DO j=1,ERROR+1
			WRITE(20,500)(H(i,j),i=1,ERROR+1)
		ENDDO
	ENDIF
	WRITE(20,500)(DT(i,:),i=1,DN*DN)
500	FORMAT(1000000(E10.3))
	CLOSE(10)
	CLOSE(20)
END PROGRAM MAIN

SUBROUTINE NNP(ii,p,jj)
	USE GLOBAL
	IMPLICIT NONE
	INTEGER :: ii,i,j,jj(4),p
LP:	DO i=1,DN
		DO j=1,DN
			IF(ii==nn(i,j)) THEN
				EXIT LP
			ENDIF
		ENDDO
	ENDDO LP
	IF(p==0) THEN
		jj(1)=nn(i,PBC(j+1))
		jj(2)=nn(PBC(i+1),j)
		jj(3)=nn(i,PBC(j-1))
		jj(4)=nn(PBC(i-1),j)
	ENDIF
	IF(p==1) THEN
		jj(1)=nn(PBC(i+1),PBC(j+1))
		jj(2)=nn(PBC(i+1),PBC(j-1))
		jj(3)=nn(PBC(i-1),PBC(j+1))
		jj(4)=nn(PBC(i-1),PBC(j-1))
	ENDIF
END SUBROUTINE NNP

SUBROUTINE INIT_RANDOM_SEED()
	INTEGER :: i, n, clock
	INTEGER, DIMENSION(:), ALLOCATABLE :: seed
	CALL RANDOM_SEED(size = n)
	ALLOCATE(seed(n))
	CALL SYSTEM_CLOCK(COUNT=clock)
	seed = clock + 37 * (/ (i - 1, i = 1, n) /)
	CALL RANDOM_SEED(PUT = seed)
	DEALLOCATE(seed)
END SUBROUTINE

	