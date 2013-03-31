MODULE GLOBAL
	IMPLICIT NONE
	SAVE
	INTEGER(8),PARAMETER :: DN=4,DN2=16,NT=3,NS=3
	INTEGER(8) :: PBC(DN2,0:1)=0
	REAL(8) :: T=-1D0,U=0D0
END MODULE
MODULE LANCZOS
	USE GLOBAL
	IMPLICIT NONE
	INTEGER(8),POINTER :: ABIT(:),ABITD(:),ABITU(:)
	INTEGER(8),ALLOCATABLE,TARGET :: AD(:),AU(:),CA(:)
	CONTAINS
	SUBROUTINE INITIALSTATE(A,SZ,n)
		IMPLICIT NONE
		INTEGER(8) :: A(:),i,n,j,k,SZ
		A=0
		DO i=1,n
			A(1)=IBSET(A(1),i-1)
		ENDDO
		i=1
EX:		DO WHILE(.TRUE.)
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
		SZ=i-1
	END SUBROUTINE
	SUBROUTINE FINDSTATE(A,n,SZ)
		IMPLICIT NONE
		INTEGER(8) :: A,n,i,ia,ib,SZ
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
			! WRITE(*,*)ia,i,ib
		ENDDO
		n=i
	END SUBROUTINE
	SUBROUTINE HUBBARD(VA,VB,SZ)
		IMPLICIT NONE
		INTEGER(8) :: i,ii,iii,iA(0:1),M,K,L,A,j,iON,AON,n,TMP,SZ(2),SIG
		LOGICAL :: FLAG
		REAL(8) :: VA(SZ(1)*SZ(2)),VB(SZ(1)*SZ(2))
		!$OMP PARALLEL DO REDUCTION(+:VB) PRIVATE(iA,AON,iON,M,K,L,A,j,TMP,ABIT,FLAG,SIG) SCHEDULE(GUIDED)
		DO n=1,SZ(1)*SZ(2)
			! WRITE(*,*)n
			iA(0)=(n-1)/SZ(2)+1
			iA(1)=n-(n-1)/SZ(2)*SZ(2)
			! WRITE(*,"(2B17.16,2I5)")ABITD(iA(0)),ABITU(iA(1)),iA
			AON=IAND(ABITD(iA(0)),ABITU(iA(1)))
			iON=0
			DO i=1,DN2
				DO ii=0,3
					IF(ii/2==0) THEN
						ABIT=>ABITD
					ELSE
						ABIT=>ABITU
					ENDIF
					FLAG=.FALSE.
					SIG=1
					! WRITE(*,"(B17.16)")ABIT(iA(ii/2))
					M=PBC(i,MOD(ii,2))
					! WRITE(*,"(B17.16)")M
					K=IAND(ABIT(iA(ii/2)),M)
					! WRITE(*,"(B17.16)")K
					IF(K==M.OR.K==0) THEN
						CYCLE
					ENDIF
					L=IEOR(K,M)
					! WRITE(*,"(B17.16)")L
					A=ABIT(iA(ii/2))-K+L
					TMP=iA(ii/2)
					! WRITE(*,*)"a=",a
					CALL FINDSTATE(A,iA(ii/2),SZ(ii/2+1))
					! WRITE(*,*)iA(ii/2)
					DO iii=0,DN2-1
						IF(BTEST(M,iii)) THEN
							FLAG=.NOT.FLAG
						ELSEIF(FLAG.AND.BTEST(ABIT(TMP),iii)) THEN
							SIG=-1*SIG
						ENDIF
					ENDDO
					! WRITE(10,"(2B17.16I2)")ABIT(TMP),ABIT(iA(ii/2)),SIG
					j=(iA(0)-1)*SZ(2)+iA(1)
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
	END SUBROUTINE
	SUBROUTINE LANC(V,AP,BT,E,SZ,M,FLAG,INFO)
		IMPLICIT NONE
		INTEGER(8) :: SZ(2),MODI(0:1)
		INTEGER :: INFO,is,M
		REAL(8) :: V(SZ(1)*SZ(2),0:2),AP(M),BT(0:M),E(M),TMP(0:M),FE(M),WORK(2*M-2),SUN
		REAL(8),ALLOCATABLE :: Z(:,:),EV(:)
		LOGICAL :: FLAG
		M=400
		V(:,0)=0D0
		FE(1)=1000
		BT(0)=0D0
		DO is=1,M
			MODI(1)=MOD(is,2)
			MODI(0)=MOD(is-1,2)
			V(:,MODI(0))=-1*BT(is-1)*V(:,MODI(0))
			CALL HUBBARD(V(:,MODI(1)),V(:,MODI(0)),SZ)
			AP(is)=DOT_PRODUCT(V(:,MODI(1)),V(:,MODI(0)))
			IF(is>10) THEN
				E(1:is)=AP(1:is)
				TMP(1:is-1)=BT(1:is-1)
				CALL DSTEQR('N',is,E(1:is),TMP(1:is-1),WORK,is,WORK,INFO)
				IF(ABS(E(1)-FE(1))<1E-12.OR.INFO/=0) THEN
					IF(FLAG) THEN
						E(1:is)=AP(1:is)
						TMP(1:is-1)=BT(1:is-1)
						ALLOCATE(Z(is,is),EV(is))
						CALL DSTEQR('I',is,E(1:is),TMP(1:is-1),Z,is,WORK,INFO)
						EV=Z(:,1)
					ENDIF
					M=is
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
		IF(FLAG) THEN
			V(:,0)=0D0
			V(:,1)=V(:,2)
			V(:,2)=0D0
			DO is=1,M
				MODI(1)=MOD(is,2)
				MODI(0)=MOD(is-1,2)
				V(:,2)=V(:,2)+EV(is)*V(:,MODI(1))
				V(:,MODI(0))=-1*BT(is-1)*V(:,MODI(0))
				CALL HUBBARD(V(:,MODI(1)),V(:,MODI(0)),SZ)
				V(:,MODI(0))=V(:,MODI(0))-AP(is)*V(:,MODI(1))
				V(:,MODI(0))=V(:,MODI(0))/BT(is)
			ENDDO
		ENDIF
	END SUBROUTINE
	SUBROUTINE CALCU(V,AP,BT,SZ,M,INFO)
		IMPLICIT NONE
		INTEGER(8) :: SZ(2),MODI(0:1)
		INTEGER :: INFO,is,M
		REAL(8) :: V(SZ(1)*SZ(2),0:2),AP(M),BT(0:M),SUN
		REAL(8),ALLOCATABLE :: Z(:,:),EV(:)
		V(:,0)=0D0
		BT(0)=0D0
		DO is=1,M
			MODI(1)=MOD(is,2)
			MODI(0)=MOD(is-1,2)
			V(:,MODI(0))=-1*BT(is-1)*V(:,MODI(0))
			CALL HUBBARD(V(:,MODI(1)),V(:,MODI(0)),SZ)
			AP(is)=DOT_PRODUCT(V(:,MODI(1)),V(:,MODI(0)))
			V(:,MODI(0))=V(:,MODI(0))-AP(is)*V(:,MODI(1))
			SUN=DOT_PRODUCT(V(:,MODI(0)),V(:,MODI(0)))
			BT(is)=SQRT(ABS(SUN))
			V(:,MODI(0))=V(:,MODI(0))/BT(is)
		ENDDO
	END SUBROUTINE
	SUBROUTINE CREATE(V,CV,i,SPIN,SZ)
		IMPLICIT NONE
		INTEGER(8) :: i,SPIN,SZ(2),ii,j,n,SIG,TMP,X
		REAL(8),POINTER :: CV(:)
		REAL(8) :: V(:)
		INTEGER(8),ALLOCATABLE :: A(:)
		ALLOCATE(A(1000000))
		IF(SPIN==0) THEN
			TMP=SZ(1)
			CALL INITIALSTATE(A,SZ(1),(NT-NS)/2+1)
			ALLOCATE(CA(SZ(1)))
			CA=A(1:SZ(1))
			DEALLOCATE(A)
			ALLOCATE(CV(SZ(2)*SZ(1)))
			CV=0D0
			DO n=1,TMP
				ABIT=>CA
				IF(.NOT.BTEST(ABITD(n),i-1)) THEN
					SIG=1**((NT+NS)/2)
					X=IBSET(ABITD(n),i-1)
					CALL FINDSTATE(X,j,SZ(1))
					DO ii=i,DN2
						IF(BTEST(ABIT(n),i)) THEN
							SIG=-1*SIG
						ENDIF
					ENDDO
					CV((j-1)*SZ(2)+1:j*SZ(2))=CV((j-1)*SZ(2)+1:j*SZ(2))+V((n-1)*SZ(2)+1:n*SZ(2))*SIG
				ENDIF
			ENDDO
			ABITD=>CA
		ELSE
			TMP=SZ(2)
			CALL INITIALSTATE(A,SZ(2),(NT+NS)/2+1)
			ALLOCATE(CA(SZ(2)))
			CA=A(1:SZ(2))
			DEALLOCATE(A)
			ALLOCATE(CV(SZ(2)*SZ(1)))
			CV=0D0
			ABIT=>CA
			DO n=1,TMP
				IF(.NOT.BTEST(ABITU(n),i-1)) THEN
					SIG=1
					X=IBSET(ABITU(n),i-1)
					CALL FINDSTATE(X,j,SZ(2))
					DO ii=i,DN2
						IF(BTEST(ABIT(n),i)) THEN
							SIG=-1*SIG
						ENDIF
					ENDDO
					CV(j:SZ(2)*(SZ(1)-1)+j:SZ(2))=CV(j:SZ(2)*(SZ(1)-1)+j:SZ(2))+V(n:TMP*(SZ(1)-1)+n:TMP)*SIG
				ENDIF
			ENDDO
			ABITU=>CA
		ENDIF
	END SUBROUTINE
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
END MODULE
PROGRAM MAIN
	USE GLOBAL
	USE LANCZOS
	IMPLICIT NONE
	REAL(8),ALLOCATABLE :: V(:,:),AP(:),BT(:),E(:),V2(:)
	REAL(8),POINTER :: CV(:)
	INTEGER(8) :: i,j,k,im(2),TMP(2),SPIN
	INTEGER(8),ALLOCATABLE :: A(:)
	REAL(8) :: SUN,SITE(DN2,2),Y(-1000:500)
	INTEGER :: M=400,INFO
	COMPLEX(8) :: Z,X
	OPEN(UNIT=10,FILE="../DATA/OUTPUT.DAT")
	ALLOCATE(A(1000000))
	DO i=1,DN2
		PBC(i,:)=IBSET(0,i-1)
		PBC(i,0)=IBSET(PBC(i,0),MOD(MOD(i-1,DN)+1,DN)+(i-1)/DN*DN)
		PBC(i,1)=IBSET(PBC(i,1),MOD(i+DN-1,DN2))
	ENDDO
	CALL INITIALSTATE(A,im(1),(NT-NS)/2)
	ALLOCATE(AD(im(1)))
	AD=A(1:im(1))
	IF(NS==0) THEN
		im(2)=im(1)
		ALLOCATE(AU(im(2)))
		AU=AD
	ELSE
		CALL INITIALSTATE(A,im(2),(NT+NS)/2)
		ALLOCATE(AU(im(2)))
		AU=A(1:im(2))
	ENDIF
	DEALLOCATE(A)
	ALLOCATE(V(im(1)*im(2),0:2),AP(M),BT(0:M),E(M))
	CALL INIT_RANDOM_SEED()
	CALL RANDOM_NUMBER(V(:,2))
	V(:,1)=V(:,2)
	SUN=DOT_PRODUCT(V(:,1),V(:,1))
	V(:,1)=V(:,1)/SQRT(ABS(SUN))
	ABITD=>AD
	ABITU=>AU
	CALL LANC(V,AP,BT,E,im,M,.TRUE.,INFO)
	WRITE(10,"(E17.6)")E(1:M)
	WRITE(10,*)"----------------------------------------------"
	! WRITE(10,"(E17.6)")V(:,2)
	! WRITE(10,*)"----------------------------------------------"
	! V(:,1)=0D0
	! CALL HUBBARD(V(:,2),V(:,1),im)
	! WRITE(10,"(E17.6)")V(:,1)-E(1)*V(:,2)
	SUN=DOT_PRODUCT(V(:,2),V(:,2))
	ALLOCATE(V2(im(1)*im(2)))
	V2=V(:,2)/SQRT(ABS(SUN))
	Y=0D0
	DO k=1,DN2
		TMP=im
		ABITU=>AU
		ABITD=>AD
		CALL CREATE(V2,CV,MOD(k-1,16)+1,(k-1)/16,im)
		DEALLOCATE(V)
		ALLOCATE(V(im(2)*im(1),0:2))
		V(:,2)=CV
		DEALLOCATE(CV)
		V(:,1)=V(:,2)
		M=400
		CALL CALCU(V,AP,BT,im,M,INFO)
		DEALLOCATE(CA)
		WRITE(*,*)k
		DO j=-1000,500
			X=0D0
			BT(0)=DOT_PRODUCT(V(:,2),V(:,2))
			Z=CMPLX(10D0/1000D0*j+E(1),0.0001)
			DO i=M,1,-1
				X=BT(i-1)**2/(Z-AP(i)-X)
			ENDDO
			Y(j)=Y(j)+DIMAG(X)
		ENDDO
		im=TMP
	ENDDO
	WRITE(10,"(2E17.6)")(10D0/1000D0*j+E(1),Y(j),j=-1000,500)
END