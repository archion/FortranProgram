MODULE GLOBAL
	IMPLICIT NONE
	REAL(8),SAVE :: DV=1.0,nf0=0.85,DU=0D0,T=1,TP=0D0,eff=1.0
	INTEGER, PARAMETER :: DN=24
	REAL(8), PARAMETER :: PI=3.14159
	INTEGER,SAVE,ALLOCATABLE :: nn(:,:),PBC(:)
END MODULE
PROGRAM MAIN
	USE GLOBAL
	IMPLICIT NONE
	CHARACTER(24) :: NOW
	INTEGER :: i,j,k=1,Cutoff,n,n_MOD(3),Im(2)=-1,INFO=0,NP(4),ERROR=0,INI=0
	REAL(8), ALLOCATABLE :: n1(:),n2(:),n1p(:),n2p(:),E(:),DTR(:,:),hhn(:,:),een(:,:)
	REAL(8), ALLOCATABLE :: H(:,:),H_1(:,:),DT(:,:),DT1(:,:),DTF(:)
	REAL(8) :: En,sa,sb,sp,s0,nf,CACHE,a,b,AC_ba,L=4
	ALLOCATE(H(2*DN*DN,2*DN*DN),H_1(2*DN*DN,2*DN*DN),n1(DN*DN),n2(DN*DN),n1p(DN*DN),n2p(DN*DN),DT(DN*DN,DN*DN),&
		DT1(DN*DN,DN*DN),DTR(DN*DN,DN*DN),DTF(DN*DN),nn(DN,DN),PBC(0:DN+1),hhn(0:2,2*DN*DN),&
		een(0:2,2*DN*DN))
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
	DT1=0.5
	! CALL INIT_RANDOM_SEED()
	! CALL RANDOM_NUMBER(DTR)
	! DT1=DTR
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
	! Im(2)=Im(1)-4
	Cutoff=1000
	! CALL CPU_TIME(TIME1)
! EX:	DO WHILE(ANY(ABS(n1p-n1)>5E-5).OR.ANY(ABS(n2p-n2)>5E-5))
EX:	DO WHILE(ANY(ABS(DT1-DT)>5E-5).AND.(INI<25))
		INI=INI+1
		CALL FDATE(NOW)
		WRITE(*,*)DT(1,2),DT1(1,2),NOW
		DT=DT1
		! n1=n1p
		! n2=n2p
		nf=0
		sa=0
		sb=1
		! DO WHILE(ABS(nf-nf0)>5E-5)
			! WRITE(*,*)sa,sb
			! sp=0.5*(sa+sb)
			sp=-1.5*T
			DV=-2.2*T
			eff=0.3*T
			a=10
			b=0
			AC_ba=ACOS(-b/a)
			H=0
			DO i=1,DN*DN
				! nearest neighbor pairs
				CALL NNP(i,0,NP)
				DO j=1,4
					H(2*i-1:2*i,2*NP(j)-1:2*NP(j))=H(2*i-1:2*i,2*NP(j)-1:2*NP(j))+RESHAPE((/-T,DT(i,NP(j)),DT(i,NP(j)),T/),(/2,2/))
				ENDDO
				! next-nearest neighbor pairs
				! CALL NNP(i,1,NP)
				! DO j=1,4
					! H(2*i-1:2*i,2*NP(j)-1:2*NP(j))=H(2*i-1:2*i,2*NP(j)-1:2*NP(j))+RESHAPE((/-TP,0.0D0,0.0D0,TP/),(/2,2/))
				! ENDDO
				! on site
				H(2*i-1:2*i,2*i-1:2*i)=H(2*i-1:2*i,2*i-1:2*i)+RESHAPE((/-sp,0.0D0,0.0D0,sp/),(/2,2/))
				IF(ANY(Im==i)) THEN
					H(2*i-1:2*i,2*i-1:2*i)=H(2*i-1:2*i,2*i-1:2*i)+RESHAPE((/eff,0.0D0,0.0D0,eff/),(/2,2/))
				ENDIF
			ENDDO
			! DO i=1,2*DN*DN
				! DO j=1,2*DN*DN
					! IF((H(i,j)-H(j,i))>5E-5) THEN
						! WRITE(*,*)i,j,"ERROR"
					! ENDIF
				! ENDDO
			! ENDDO
			H_1=(H-b)/a
			n1p=0
			n2p=0
			DT1=0
			!$OMP PARALLEL DO FIRSTPRIVATE(een,hhn) SCHEDULE(GUIDED,1) NUM_THREADS(32)
			DO i=1,DN*DN
				een=0
				hhn=0
				een(0,2*i-1)=1D0
				hhn(0,2*i)=1D0
				! n1p(i)=n1p(i)+een(0,2*i-1)*(1-AC_ba/PI)
				! n2p(i)=n2p(i)+hhn(0,2*i)*(1-AC_ba/PI)
				! CALL CRSMV(H_1,een(0,:),een(1,:),2*DN*DN,2*DN*DN)
				CALL CRSMV(H_1,hhn(0,:),hhn(1,:),2*DN*DN)
				! n1p(i)=n1p(i)+een(1,2*i-1)*2*SIN(AC_ba)/PI
				! n2p(i)=n2p(i)+hhn(1,2*i)*2*SIN(AC_ba)/PI
				CALL NNP(i,0,NP)
				DO j=1,4
					DT1(i,NP(j))=DT1(i,NP(j))+hhn(1,2*NP(j)-1)*SIN(AC_ba)*SINH(L*(1-1.0/REAL(Cutoff)))/SINH(L)
				ENDDO
				DO n=2,Cutoff
					n_MOD(1)=MOD(n-2,3)
					n_MOD(2)=MOD(n-1,3)
					n_MOD(3)=MOD(n,3)
					! CALL CRSMV(H_1,een(n_MOD(2),:),een(n_MOD(3),:),2*DN*DN,2*DN*DN)
					! een(n_MOD(3),:)=2*een(n_MOD(3),:)-een(n_MOD(1),:)
					CALL CRSMV(H_1,hhn(n_MOD(2),:),hhn(n_MOD(3),:),2*DN*DN)
					hhn(n_MOD(3),:)=2*hhn(n_MOD(3),:)-hhn(n_MOD(1),:)
					! n1p(i)=n1p(i)+een(n_MOD(3),2*i-1)*2*SIN(nc*AC_ba)/nc/PI
					! n2p(i)=n2p(i)+hhn(n_MOD(3),2*i)*2*SIN(nc*AC_ba)/nc/PI
					CALL NNP(i,0,NP)
					DO j=1,4
						DT1(i,NP(j))=DT1(i,NP(j))+hhn(n_MOD(3),2*NP(j)-1)*SIN(n*AC_ba)/n*SINH(L*(1-REAL(n)/REAL(Cutoff)))/SINH(L)
					ENDDO
				ENDDO
			ENDDO
			!$OMP END PARALLEL DO
			! n2p=1-n2p
			! IF(ANY(n2p<0)) THEN
				! EXIT EX
			! ENDIF
			! nf=0
			! DO i=1,DN*DN
				! nf=nf+n1p(i)+n2p(i)
			! ENDDO
			! nf=nf/(DN*DN)
			! WRITE(*,"(E15.6)")nf
			! IF (nf<nf0) THEN
				! sa=sp
			! ELSE
				! sb=sp
			! ENDIF
			! ! WRITE(*,"(3E11.3)")nf,sa,sb
		! ENDDO
		DT1=-2*DV*DT1/PI
		! WRITE(*,"(2E15.6,F15.12)")n1(1),n1p(1),sa
		! REWIND(10)
		! WRITE(10,400)((n1(DN*i+j)+n2(DN*i+j),j=1,DN),i=0,DN-1)
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
		REWIND(10)
		DO i=1,DN*DN
			CALL NNP(i,0,NP)
			DTF(i)=(DT(i,NP(1))+DT(i,NP(3))+DT(i,NP(2))+DT(i,NP(4)))/4.0
		ENDDO
		WRITE(10,400)((DTF(DN*i+j),j=1,DN),i=0,DN-1)
	ENDDO EX
	! CALL CPU_TIME(TIME2)
	REWIND(10)
	DO i=1,DN*DN
		CALL NNP(i,0,NP)
		DTF(i)=(DT(i,NP(1))+DT(i,NP(3))-DT(i,NP(2))-DT(i,NP(4)))/4.0
	ENDDO
	! WRITE(10,300)(DT(i,:),i=1,DN*DN)
	WRITE(10,400)((DTF(DN*i+j),j=1,DN),i=0,DN-1)
	! WRITE(10,400)((n1(DN*i+j)+n2(DN*i+j),j=1,DN),i=0,DN-1)
	! WRITE(10,400)(((-1)**(i+j)*0.5*(n1(DN*i+j)-n2(DN*i+j)),j=1,DN),i=0,DN-1)
	! WRITE(10,*)"THE TATAL RUNTIME IS:",TIME1
	WRITE(10,"(16E12.3)")(DT(i,:),i=1,DN*DN)
	WRITE(10,"(32E12.3)")(H(i,:),i=1,2*DN*DN)
	! IF(INFO/=0) THEN
		! ERROR=MOD(INFO,2*DN*DN+1)
		! DO j=1,ERROR+1
			! WRITE(20,500)(H(i,j),i=1,ERROR+1)
		! ENDDO
! 500		FORMAT(1000000("(",E10.3,",",E10.3,")  "))
	! ENDIF
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

	
	

	
! SUBROUTINE CRSMV(A,X,Y,M)
	! IMPLICIT NONE
	! INTEGER :: i,JOB(6),INFO,M
	! REAL(8) :: A(M,M),X(*),Y(*)
	! INTEGER :: JA(10*M),IA(M+1)
	! REAL(8) :: ACSR(10*M)
	! ! INTEGER, ALLOCATABLE :: JA(:),IA(:)
	! ! REAL(8), ALLOCATABLE :: ACSR(:)
	! JOB=(/0,1,1,2,M*M,1/)
	! ! ALLOCATE(ACSR(M*M),JA(M*M),IA(M+1))
	! ! ALLOCATE(ACSR(M*M))
	! CALL MKL_DDNSCSR(JOB, M, M, A, M, ACSR, JA, IA, INFO)
	! CALL MKL_DCSRGEMV("N", M, ACSR, IA, JA, X, Y)
! END
SUBROUTINE CRSMV(A,X,Y,M)
	IMPLICIT NONE
	INTEGER :: i=0,j=0,k1,k2,M
	INTEGER :: JA(10*M+1),IA(M+1)
	REAL(8) :: VA(10*M)
	REAL(8) :: A(M,M),X(*),Y(*)
	k1=1
	k2=1
	DO i=1,M
		IA(k2)=K1
		k2=k2+1
		DO j=1,M
			IF(ABS(A(i,j))>5E-5) THEN
				VA(K1)=A(i,j)
				JA(k1)=j
				k1=k1+1
			ENDIF
		ENDDO
	ENDDO
	IA(k2)=k1
	DO i=1,M
		Y(i)=0
		DO j=IA(i),IA(i+1)-1
			Y(i)=Y(i)+VA(j)*X(JA(j))
		ENDDO
	ENDDO
END
! SUBROUTINE CRSMV(A,X,Y,M)
	! IMPLICIT NONE
	! INTEGER :: i,JOB(6),INFO,M
	! REAL(8) :: A(M,M),X(*),Y(*)
	! INTEGER, ALLOCATABLE :: JA(:),IA(:)
	! REAL(8), ALLOCATABLE :: ACSR(:)
	! JOB=(/0,1,1,2,M*M,1/)
	! ALLOCATE(ACSR(M*M),JA(M*M),IA(M+1))
	! CALL MKL_DDNSCSR(JOB, M, M, A, M, ACSR, JA, IA, INFO)
	! CALL MKL_DCSRGEMV("N", M, ACSR, IA, JA, X, Y)
! END
