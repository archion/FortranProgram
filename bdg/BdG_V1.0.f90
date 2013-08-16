MODULE GLOBAL
	IMPLICIT NONE
	REAL,SAVE :: BT=1E5,DV=1.0,nf0=0.85,DU=2.44,T=1,TP=-0.25,eff=1.0
	INTEGER, PARAMETER :: DN=24
	INTEGER,SAVE,ALLOCATABLE :: nn(:,:),PBC(:)
END MODULE
PROGRAM MAIN
	USE GLOBAL
	IMPLICIT NONE
	INTEGER :: i,j,k=1,CC,n,Im(2),INFO=0,NP(4),ERROR=0
	INTEGER, ALLOCATABLE :: IWORK(:)
	REAL, ALLOCATABLE :: n1(:),n2(:),n1p(:),n2p(:),E(:),F1(:,:),F2(:,:),RWORK(:),DTR(:,:)
	COMPLEX, ALLOCATABLE :: H(:,:),DT(:,:),DT1(:,:),DTF(:),u(:,:),v(:,:),AP(:,:),WORK(:)
	REAL :: f,En,sa,sb,sp,s0,nf,CACHE,TIME1,TIME2
	f(En)=1/(1+EXP(BT*En))
	ALLOCATE(H(2*DN*DN,2*DN*DN),n1(DN*DN),n2(DN*DN),n1p(DN*DN),n2p(DN*DN),DT(DN*DN,DN*DN),DT1(DN*DN,DN*DN),DTR(DN*DN,DN*DN),&
		DTF(DN*DN),u(DN*DN,2*DN*DN),v(DN*DN,2*DN*DN),E(2*DN*DN),WORK(4*DN*DN+(2*DN*DN)**2),RWORK(1+10*DN*DN+2*(2*DN*DN)**2),&
		IWORK(3+10*DN*DN),AP(2*DN*DN,2*DN*DN),nn(DN,DN),PBC(0:DN+1),F1(DN*DN,DN*DN),F2(DN*DN,DN*DN))
	OPEN(UNIT=10,FILE='../DATA/output.dat')
	OPEN(UNIT=20,FILE='../DATA/ERROR.dat')
200	FORMAT(24("(",E15.6,",",E15.6,")  "))
300	FORMAT(576("(",E15.6,",",E15.6,")  "))
400	FORMAT(24E15.6)
	PBC=(/DN,(i,i=1,DN),1/)
	DO i=1,DN
		DO j=1,DN
			nn(i,j)=k
			k=k+1
		ENDDO
	ENDDO
	DT=0
	DT1=0.001
	CALL INIT_RANDOM_SEED()
	CALL RANDOM_NUMBER(DTR)
	DT1=DTR
	! CALL RANDOM_NUMBER(DTR)
	! DT1=DT1+CMPLX(0,DTR)
	! DT1=-0.075
	! DO i=1,DN*DN
		! DT1(i,i)=0.08
	! ENDDO
	n1p=nf0/2
	n2p=nf0/2
	n1=0
	n2=0
	Im(1)=DN/2*DN-DN/2
	Im(2)=Im(1)
	CC=0
	CALL CPU_TIME(TIME1)
EX:	DO WHILE(ANY(ABS(n1p-n1)>5E-5).OR.ANY(ABS(n2p-n2)>5E-5))
		! F1=DT-DT1
		DT=DT1
		n1=n1p
		n2=n2p
		nf=0
		! s0=sa
		! IF(CC<=1) THEN
			sa=0
			sb=1
		! ENDIF
		DO WHILE(ABS(nf-nf0)>5E-5)
			sp=0.5*(sa+sb)
			H=0
			DO i=1,DN*DN
				DO j=1,DN*DN
					! nearest neighbor pairs
					CALL NNP(i,0,NP)
					IF(ANY(NP==j)) THEN
						H(2*i-1:2*i,2*j-1:2*j)=H(2*i-1:2*i,2*j-1:2*j)+RESHAPE((/CMPLX(-T,0),DT(i,j),DT(i,j),CMPLX(T,0)/),(/2,2/))
						CYCLE
					ENDIF
					! next-nearest neighbor pairs
					CALL NNP(i,1,NP)
					IF(ANY(NP==j)) THEN
						H(2*i-1:2*i,2*j-1:2*j)=H(2*i-1:2*i,2*j-1:2*j)+RESHAPE((/-TP,0.0,0.0,TP/),(/2,2/))
						CYCLE
					ENDIF
					! on site
					IF(i==j) THEN
						H(2*i-1:2*i,2*i-1:2*i)=H(2*i-1:2*i,2*i-1:2*i)+RESHAPE((/DU*n2(i)-sp,0.0,0.0,-DU*n1(i)+sp/),&
							(/2,2/))
						IF(ANY(Im==i)) THEN
							H(2*i-1:2*i,2*i-1:2*i)=H(2*i-1:2*i,2*i-1:2*i)+RESHAPE((/eff,0.0,0.0,eff/),(/2,2/))
						ENDIF
					ENDIF
				ENDDO
			ENDDO
			AP=0
			DO i=1,2*DN*DN
				DO j=i,2*DN*DN
					AP(i,j)=H(i,j)
				ENDDO
			ENDDO
			CALL CHEEVD("V","U",2*DN*DN,AP,2*DN*DN,E,WORK,4*DN*DN+(2*DN*DN)**2,RWORK,&
				1+10*DN*DN+2*(2*DN*DN)**2,IWORK,3+10*DN*DN,INFO)
			! CALL DSPEV('V','U',2*DN*DN,AP,E,Z,2*DN*DN,WORK,INFO)
			IF(INFO/=0) THEN
				EXIT EX
			ENDIF
			DO i=1,2*DN*DN
				IF(MOD(i,2)==1) THEN
					u((i+1)/2,:)=AP(i,:)
				ELSE
					v(i/2,:)=AP(i,:)
				ENDIF
			ENDDO
			! WRITE(*,"(E15.6)")E
			nf=0
			DO i=1,DN*DN
				DO n=1,2*DN*DN
					nf=nf+(u(i,n)*CONJG(u(i,n))*f(E(n))+v(i,n)*CONJG(v(i,n))*(1-f(E(n))))
				ENDDO
			ENDDO
			nf=1.0/(DN*DN)*nf
			IF (nf<nf0) THEN
				sa=sp
			ELSE
				sb=sp
			ENDIF
			! WRITE(*,"(3E11.3)")nf,sa,sb
		ENDDO
		n1p=0
		n2p=0
		DO i=1,DN*DN
			DO n=1,2*DN*DN
				n1p(i)=n1p(i)+u(i,n)*CONJG(u(i,n))*f(E(n))
				n2p(i)=n2p(i)+v(i,n)*CONJG(v(i,n))*(1-f(E(n)))
			ENDDO
		ENDDO
		DT1=0
		DO i=1,DN*DN
			CALL NNP(i,0,NP)
			DO j=1,DN*DN
				IF(ANY(NP==j)) THEN
					DO n=1,2*DN*DN
						DT1(i,j)=DT1(i,j)+(u(i,n)*CONJG(v(j,n))+CONJG(v(i,n))*u(j,n))*TANH(0.5*BT*E(n))
					ENDDO
				ELSE
					DT1(i,j)=0
				ENDIF
			ENDDO
		ENDDO
		DT1=0.25*DV*DT1
		WRITE(*,"(2E15.6,F15.12)")n1(1),n1p(1),sa
		! IF(CC>=1) THEN
			! sb=s0
			! sa=sa+sa-s0
			! IF(sa>sb) THEN
				! sb=sa
				! sa=s0
			! ENDIF
		! ENDIF
		! CC=CC+1
		REWIND(10)
		WRITE(10,400)((n1(DN*i+j)+n2(DN*i+j),j=1,DN),i=0,DN-1)
		! DO i=1,DN*DN
			! CALL NNP(i,0,NP)
			! DTF(i)=(DT(i,NP(1))+DT(i,NP(3))-DT(i,NP(2))-DT(i,NP(4)))/4.0
		! ENDDO
		! WRITE(10,200)((DTF(24*i+j),j=1,24),i=0,23)
		! F2=DT-DT1
		! IF(ALL(F1+F2<1E-5)) THEN
			! ERROR=1
			! EXIT
		! ENDIF
	ENDDO EX
	CALL CPU_TIME(TIME2)
	REWIND(10)
	DO i=1,DN*DN
		CALL NNP(i,0,NP)
		DTF(i)=(DT(i,NP(1))+DT(i,NP(3))-DT(i,NP(2))-DT(i,NP(4)))/4.0
	ENDDO
	! WRITE(10,300)(DT(i,:),i=1,DN*DN)
	WRITE(10,200)((DTF(DN*i+j),j=1,DN),i=0,DN-1)
	WRITE(10,400)((n1(DN*i+j)+n2(DN*i+j),j=1,DN),i=0,DN-1)
	WRITE(10,400)(((-1)**(i+j)*0.5*(n1(DN*i+j)-n2(DN*i+j)),j=1,DN),i=0,DN-1)
	WRITE(10,*)"The beta is",BT,"V=",DV,"U=",DU,"t'=",TP,"heff=",eff,"nf=",nf0,"impure",MOD(1,Im(2)-Im(1)),"N=",DN
	WRITE(10,*)"ERROR=",ERROR,"INFO=",INFO
	WRITE(10,*)"THE TATAL RUNTIME IS:",TIME2-TIME1
	IF(INFO/=0) THEN
		ERROR=MOD(INFO,2*DN*DN+1)
		DO j=1,ERROR+1
			WRITE(20,500)(H(i,j),i=1,ERROR+1)
		ENDDO
500		FORMAT(1000000("(",E10.3,",",E10.3,")  "))
	ENDIF
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

	
	

	