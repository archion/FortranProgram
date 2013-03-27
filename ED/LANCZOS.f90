MODULE GLOBAL
	IMPLICIT NONE
	SAVE
	INTEGER(8),PARAMETER :: DN=4,DN2=16,NS=2
	INTEGER(8) :: PBC(DN2,0:1)=0
	REAL(8) :: T=-1D0,U=4D0
END MODULE
PROGRAM LANCZOS
	USE GLOBAL
	IMPLICIT NONE
	REAL(8),ALLOCATABLE :: V(:,:),AP(:),BT(:),FE(:),E(:),TMP(:),WORK(:),Z(:,:),EV(:)
	INTEGER(8) :: i,j,k,im,MODI(0:2)
	INTEGER(8),ALLOCATABLE :: A(:),ABIT(:)
	REAL(8) :: SUN
	INTEGER :: INFO,is,M=400
	OPEN(UNIT=10,FILE="../DATA/OUTPUT.DAT")
	ALLOCATE(A(1000000))
	DO i=1,DN2
		PBC(i,:)=IBSET(0,i-1)
		PBC(i,0)=IBSET(PBC(i,0),MOD(MOD(i-1,DN)+1,DN)+(i-1)/DN*DN)
		PBC(i,1)=IBSET(PBC(i,1),MOD(i+DN-1,DN2))
	ENDDO
	DO i=1,NS
		A(1)=IBSET(A(1),i-1)
	ENDDO
	i=1
EX:	DO WHILE(.TRUE.)
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
	im=i-1
	ALLOCATE(ABIT(im),V(im*im,0:2),AP(M),BT(0:M),TMP(0:M),E(M),FE(M),WORK(2*M-2))
	ABIT=A(1:im)
	DEALLOCATE(A)
	V=0D0
	FE(1)=1000
	V(:,0)=0D0
	BT(0)=0D0
	CALL INIT_RANDOM_SEED()
	CALL RANDOM_NUMBER(V(:,2))
	V(:,1)=V(:,2)
	SUN=DOT_PRODUCT(V(:,1),V(:,1))
	V(:,1)=V(:,1)/SQRT(ABS(SUN))
	DO is=1,M
		MODI(1)=MOD(is,2)
		MODI(0)=MOD(is-1,2)
		V(:,MODI(0))=-1*BT(is-1)*V(:,MODI(0))
		CALL HUBBARD(V(:,MODI(1)),V(:,MODI(0)),ABIT,im)
		AP(is)=DOT_PRODUCT(V(:,MODI(1)),V(:,MODI(0)))
		IF(is>10) THEN
			E(1:is)=AP(1:is)
			TMP(1:is-1)=BT(1:is-1)
			CALL DSTEQR('N',is,E(1:is),TMP(1:is-1),WORK,is,WORK,INFO)
			IF(ABS(E(1)-FE(1))<1E-12.OR.INFO/=0) THEN
				E(1:is)=AP(1:is)
				TMP(1:is-1)=BT(1:is-1)
				ALLOCATE(Z(is,is),EV(is))
				CALL DSTEQR('I',is,E(1:is),TMP(1:is-1),Z,is,WORK,INFO)
				EV=Z(:,1)
				WRITE(*,*)"INFO=",INFO,is
				EXIT
			ENDIF
			FE=E
		ENDIF
		V(:,MODI(0))=V(:,MODI(0))-AP(is)*V(:,MODI(1))
		SUN=DOT_PRODUCT(V(:,MODI(0)),V(:,MODI(0)))
		BT(is)=SQRT(ABS(SUN))
		V(:,MODI(0))=V(:,MODI(0))/BT(is)
	ENDDO
	! WRITE(10,"(E17.6)")E(1:is)
	V(:,1)=V(:,2)
	V(:,2)=0D0
	DO is=1,SIZE(EV)
		MODI(1)=MOD(is,2)
		MODI(0)=MOD(is-1,2)
		V(:,2)=V(:,2)+EV(is)*V(:,MODI(1))
		V(:,MODI(0))=-1*BT(is-1)*V(:,MODI(0))
		CALL HUBBARD(V(:,MODI(1)),V(:,MODI(0)),ABIT,im)
		V(:,MODI(0))=V(:,MODI(0))-AP(is)*V(:,MODI(1))
		V(:,MODI(0))=V(:,MODI(0))/BT(is)
	ENDDO
	WRITE(10,"(E17.6)")E(1:SIZE(EV))
	WRITE(10,*)"----------------------------------------------"
	WRITE(10,"(E17.6)")V(:,2)
	! V(:,1)=0D0
	! CALL HUBBARD(V(:,2),V(:,1),ABIT,im)
	! WRITE(10,"(E17.6)")V(:,1)-E(1)*V(:,2)
	! V(2,1)=1
	! CALL HUBBARD(V(:,1),V(:,2),ABIT,im)
	! DO i=1,im*im
		! IF(ABS(V(i,2))>1E-5) THEN
			! WRITE(10,"(E17.6,B17.16,B17.16)")V(i,2),ABIT(i-(i-1)/im*im),ABIT((i-1)/im+1)
		! ENDIF
	! ENDDO
END
SUBROUTINE FINDSTATE(A,n,ABIT,SZ)
	IMPLICIT NONE
	INTEGER(8) :: A,n,i,ia,ib,ABIT(*),SZ
	ia=1
	ib=SZ+1
	i=(ia+ib)/2
	DO WHILE(A/=ABIT(i))
		IF(A<ABIT(i)) THEN
			ib=i
		ELSE
			ia=i
		ENDIF
		i=(ia+ib)/2
		! WRITE(10,*)ia,i,ib
	ENDDO
	n=i
END
SUBROUTINE HUBBARD(VA,VB,ABIT,SZ)
	USE GLOBAL
	IMPLICIT NONE
	INTEGER(8) :: i,ii,iii,iA(0:1),M,K,L,A,j,iON,AON,n,TMP,SZ,ABIT(SZ),SIG
	LOGICAL :: FLAG
	REAL(8) :: VA(SZ*SZ),VB(SZ*SZ)
	!$OMP PARALLEL DO REDUCTION(+:VB) PRIVATE(iA,AON,iON,M,K,L,A,j,TMP) SCHEDULE(GUIDED)
	DO n=1,SZ*SZ
		! WRITE(*,*)n
		iA(0)=(n-1)/SZ+1
		iA(1)=n-(n-1)/SZ*SZ
		! WRITE(10,"(2B17.16)")ABIT(iA(0)),ABIT(iA(1))
		AON=IAND(ABIT(iA(0)),ABIT(iA(1)))
		iON=0
		DO i=1,DN2
			DO ii=0,3
				FLAG=.FALSE.
				SIG=1
				M=PBC(i,MOD(ii,2))
				K=IAND(ABIT(iA(ii/2)),M)
				IF(K==M.OR.K==0) THEN
					CYCLE
				ENDIF
				L=IEOR(K,M)
				A=ABIT(iA(ii/2))-K+L
				TMP=iA(ii/2)
				CALL FINDSTATE(A,iA(ii/2),ABIT,SZ)
				DO iii=0,DN2-1
					IF(BTEST(M,iii)) THEN
						FLAG=.NOT.FLAG
					ELSEIF(FLAG.AND.BTEST(ABIT(TMP),iii)) THEN
						SIG=-1*SIG
					ENDIF
				ENDDO
				! WRITE(10,"(2B17.16I2)")ABIT(TMP),ABIT(iA(ii/2)),SIG
				j=(iA(0)-1)*SZ+iA(1)
				VB(j)=VB(j)+T*SIG*VA(n)
				iA(ii/2)=TMP
			ENDDO
			IF(BTEST(AON,i-1)) THEN
				iON=iON+1
			ENDIF
		ENDDO
		VB(n)=VB(n)+U*iON*VA(n)
	ENDDO
	!$OMP END PARALLEL DO
END
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