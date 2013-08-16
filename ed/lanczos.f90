MODULE GLOBAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    �趨����ģ�Ͳ�������ʼ��                                           !
!		����:                                                           !
!     		DN*DN�������С                                             !
!     		NT����������                                                !
!     		NS��������z����                                             !
!     		T��hoppingֵ                                                !
!     		U��Hubbard U                                                !
!     		PBC�����ڸ����Ϣ                                           !
!    	�ӳ���:                                                         !
!     		INIT_LATTICE�������i��PBC����Ϊ������iλ��Ϊ1������Ϊ��    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	IMPLICIT NONE
	SAVE
	INTEGER(8),PARAMETER :: DN=4,DN2=16,NT=2
	REAL(8),PARAMETER :: PI=3.14159265358979D0
	INTEGER(8) :: PBC(DN2,0:1)=0,NS=0
	REAL(8) :: T=-1D0,U=4D0
	CONTAINS
	SUBROUTINE INIT_LATTICE()
		IMPLICIT NONE
		INTEGER :: i
		PBC(DN2,0:1)=0
		DO i=1,DN2
			PBC(i,:)=IBSET(0,i-1)
			PBC(i,0)=IBSET(PBC(i,0),MOD(MOD(i-1,DN)+1,DN)+(i-1)/DN*DN)
			PBC(i,1)=IBSET(PBC(i,1),MOD(i+DN-1,DN2))
		ENDDO
	END SUBROUTINE
END MODULE
MODULE LANCZOS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       LANCZOS����                                                                                                  !
!       ������GLOBAL                                                                                                 !
!       ������                                                                                                       !
!       	AD��AU�������£������ϵ���������                                                                         !
!       	CA����ԭ���Ļ����ϲ���һ�����ӵ�������������Ҫ�ֶ��ͷ��ڴ�                                               !
!       	ABITD��ABITU���ֱ�ָ�������£�����������������ָ��                                                       !
!       	ABIT��ָ��ABITD��ABITU��ָ�룬���ڲ���̬                                                                 !
!       �ӳ���                                                                                                     !
!       	INITIALSTATE����ʼ����������                                                                             !
!       		������                                                                                               !
!       			A���������������                                                                                !
!       			SZ�����������������С                                                                           !
!       			n�����Ӹ���                                                                                      !
!       	FINDSTATE������ֵ̬����Ӧ��ָ�꣬���ҵ�����������ָ��ABITָ��                                            !
!       		������                                                                                               !
!       			A�����룺���ҵ�ֵ̬                                                                              !
!       			SZ�����룺���ҵ�����������С                                                                     !
!       			n�������̬��ָ��                                                                                !
!       	HUBBARD�����ܶ������õ������ϣ�ʹ��ǰ��Ҫָ��ABITD��ABITU��hopping������PBCʵ��                          !
!       		������                                                                                               !
!       			VA�����룺���õ�����                                                                             !
!       			VB�����룺û�г�ʼ������������õĽ��                                                           !
!       			SZ��������ά��                                                                                   !
!       	LANC����׼LANCZOS���̣���Ҫ����HUBBARD�ͶԽǻ�����                                                       !
!       		������                                                                                               !
!       			V�����룺V(:,1:2)Ϊ��һ���ĳ�ʼ�����������V(:,2)Ϊ�����������Ѿ���һ����                        !
!       			E�����룺�����������������ֵ                                                                    !
!       			M�����룺����LANCZOS���̴�������AP��ά���������������                                         !
!       			Mo�������ʵ��ִ��LANCZOS���̵Ĵ�������E����Чά��                                               !
!       			FLAG�����룺�Ƿ���㱾��������.TRUE.Ϊ���㣬���������                                           !
!       			INFO��������Խǻ��ɹ���ָʾ��Ϣ                                                                 !
!       	CREATE��������һ��̬�ϲ���������һ�����ӵ�̬����Ҫ�õ�ָ��ABITU��ABITD���������������һָ��ָ��CA     !
!       		������                                                                                               !
!       			V�����룺��ʼ̬����̬�������������                                                              !
!       			CV�������ָ�룬ָ����̬����Ҫ�ֶ��ͷ��ڴ�                                                     !
!       			i��SPIN�����룺�������ӵ����������                                                              !
!       			SZ�����룺ԭʼ̬��ά������������̬��ά��                                                       !
!       			FLAG�����룺������������.TRUE.Ϊ����                                                           !
!       	CALCU�����㶯��ѧ�������ڴ��ͷ�CV���ڴ棩����Ҫ�õ�ָ��ABITU��ABITD                                      !
!       		������                                                                                               !
!       			CV���ض��ĳ�ʼLANCZOS̬V����CREATE�еĽ����                                                     !
!       			SZ�����룺̬��ά��                                                                               !
!       			M�����룺LANCZOS���̵Ĵ���                                                                       !
!       			FQ�����룺�����Ƶ������                                                                         !
!       			Y�����룺�����㣬���������Ķ���ѧ�����������                                                  !
!       			img�����룺ȡ����չ��                                                                            !
!       	INIT_RANDOM_SEED�����������                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE GLOBAL
	IMPLICIT NONE
	INTEGER(8),POINTER :: ABIT(:),ABITD(:),ABITU(:)
	INTEGER(8),ALLOCATABLE,TARGET :: AD(:),AU(:),CA(:)
	CONTAINS
	SUBROUTINE INITIALSTATE(P,SZ,n)
		IMPLICIT NONE
		INTEGER(8),ALLOCATABLE,TARGET :: P(:)
		INTEGER(8) :: TMP(10000000),i,n,j,k,SZ
		TMP=0
		DO i=1,n
			TMP(1)=IBSET(TMP(1),i-1)
		ENDDO
		i=1
EX:		DO WHILE(.TRUE.)
			TMP(i+1)=TMP(i)
			i=i+1
			DO j=0,DN2-1
				IF(j==DN2-1) THEN
					EXIT EX
				ENDIF
				IF(IBITS(TMP(i),j,2)==1) THEN
					TMP(i)=IBSET(TMP(i),j+1)
					TMP(i)=IBCLR(TMP(i),j)
					DO k=0,j/2-1
						IF(BTEST(TMP(i),k)) THEN
							EXIT
						ENDIF
						IF(BTEST(TMP(i),j-k-1)) THEN
							TMP(i)=IBSET(TMP(i),k)
							TMP(i)=IBCLR(TMP(i),j-k-1)
						ENDIF
					ENDDO
					EXIT
				ENDIF
			ENDDO
		ENDDO EX
		SZ=i-1
		ALLOCATE(P(SZ))
		P=TMP(1:SZ)
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
		REAL(8) :: VA(:),VB(:)
		!$OMP PARALLEL DO REDUCTION(+:VB) PRIVATE(iA,AON,iON,M,K,L,A,j,TMP,ABIT,FLAG,SIG) SCHEDULE(GUIDED)
		DO n=1,SZ(1)*SZ(2)
			iA(0)=(n-1)/SZ(2)+1
			iA(1)=n-(n-1)/SZ(2)*SZ(2)
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
					M=PBC(i,MOD(ii,2))
					K=IAND(ABIT(iA(ii/2)),M)
					IF(K==M.OR.K==0) THEN
						CYCLE
					ENDIF
					L=IEOR(K,M)
					A=ABIT(iA(ii/2))-K+L
					TMP=iA(ii/2)
					CALL FINDSTATE(A,iA(ii/2),SZ(ii/2+1))
					DO iii=0,DN2-1
						IF(BTEST(M,iii)) THEN
							FLAG=.NOT.FLAG
						ELSEIF(FLAG.AND.BTEST(ABIT(TMP),iii)) THEN
							SIG=-1*SIG
						ENDIF
					ENDDO
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
	SUBROUTINE LANC(V2,E,SZ,M,Mo,FLAG,INFO)
		IMPLICIT NONE
		INTEGER(8) :: SZ(2),MODI(0:1)
		INTEGER :: INFO,is,M,Mo
		REAL(8) :: V(SZ(1)*SZ(2),0:1),AP(M),BT(0:M),E(:),TMP(0:M),FE(M),WORK(2*M-2),SUN,V2(:)
		REAL(8),ALLOCATABLE :: Z(:,:),EV(:)
		LOGICAL :: FLAG
		V(:,1)=V2
		SUN=DOT_PRODUCT(V(:,1),V(:,1))
		V(:,1)=V(:,1)/SQRT(ABS(SUN))
		V(:,0)=0D0
		FE(1)=1000
		BT(0)=0D0
		DO is=1,M
			MODI(1)=MOD(is,2)
			MODI(0)=MOD(is-1,2)
			V(:,MODI(0))=-1*BT(is-1)*V(:,MODI(0))
			CALL HUBBARD(V(:,MODI(1)),V(:,MODI(0)),SZ)
			AP(is)=DOT_PRODUCT(V(:,MODI(1)),V(:,MODI(0)))
			IF(is>3) THEN
				E(1:is)=AP(1:is)
				TMP(1:is-1)=BT(1:is-1)
				CALL DSTEQR('N',is,E(1:is),TMP(1:is-1),WORK,is,WORK,INFO)
				IF(ABS(E(1)-FE(1))<1E-12.OR.INFO/=0.OR.BT(is-1)<1E-8) THEN
					IF(FLAG) THEN
						E(1:is)=AP(1:is)
						TMP(1:is-1)=BT(1:is-1)
						ALLOCATE(Z(is,is),EV(is))
						CALL DSTEQR('I',is,E(1:is),TMP(1:is-1),Z,is,WORK,INFO)
						EV=Z(:,1)
					ENDIF
					Mo=is
					WRITE(*,*)"INFO=",INFO
					WRITE(*,"(2E17.6,I5)")BT(is-1),ABS(E(1)-FE(1)),is
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
			V(:,1)=V2
			V2=0D0
			DO is=1,Mo
				MODI(1)=MOD(is,2)
				MODI(0)=MOD(is-1,2)
				V2=V2+EV(is)*V(:,MODI(1))
				V(:,MODI(0))=-1*BT(is-1)*V(:,MODI(0))
				CALL HUBBARD(V(:,MODI(1)),V(:,MODI(0)),SZ)
				V(:,MODI(0))=V(:,MODI(0))-AP(is)*V(:,MODI(1))
				V(:,MODI(0))=V(:,MODI(0))/BT(is)
			ENDDO
			SUN=DOT_PRODUCT(V2,V2)
			V2=V2/SQRT(SUN)
		ENDIF
	END SUBROUTINE
	SUBROUTINE CREATE(V,CV,i,SPIN,SZ,FLAG)
		IMPLICIT NONE
		INTEGER(8) :: i,SPIN,SZ(2),ii,j,n,SIG,TMP,X,C,Ci
		REAL(8),POINTER :: CV(:)
		REAL(8) :: V(:)
		LOGICAL :: FLAG
		IF(FLAG) THEN
			C=1
		ELSE
			C=-1
		ENDIF
		Ci=IBSET(0,i-1)
		IF(SPIN==0) THEN
			TMP=SZ(1)
			CALL INITIALSTATE(CA,SZ(1),(NT-NS)/2+C)
			ALLOCATE(CV(SZ(2)*SZ(1)))
			CV=0D0
			ABIT=>CA
			DO n=1,TMP
				IF(FLAG.NEQV.BTEST(ABITD(n),i-1)) THEN
					SIG=-1**((NT+NS)/2)
					X=IEOR(ABITD(n),Ci)
					CALL FINDSTATE(X,j,SZ(1))
					DO ii=i,DN2-1
						IF(BTEST(ABITD(n),ii)) THEN
							SIG=-1*SIG
						ENDIF
					ENDDO
					CV((j-1)*SZ(2)+1:j*SZ(2))=CV((j-1)*SZ(2)+1:j*SZ(2))+V((n-1)*SZ(2)+1:n*SZ(2))*SIG
				ENDIF
			ENDDO
			ABITD=>CA
		ELSE
			TMP=SZ(2)
			CALL INITIALSTATE(CA,SZ(2),(NT+NS)/2+C)
			ALLOCATE(CV(SZ(2)*SZ(1)))
			CV=0D0
			ABIT=>CA
			DO n=1,TMP
				IF(FLAG.NEQV.BTEST(ABITU(n),i-1)) THEN
					SIG=1
					X=IEOR(ABITU(n),Ci)
					CALL FINDSTATE(X,j,SZ(2))
					DO ii=i,DN2-1
						IF(BTEST(ABITU(n),ii)) THEN
							SIG=-1*SIG
						ENDIF
					ENDDO
					CV(j:SZ(2)*(SZ(1)-1)+j:SZ(2))=CV(j:SZ(2)*(SZ(1)-1)+j:SZ(2))+V(n:TMP*(SZ(1)-1)+n:TMP)*SIG
				ENDIF
			ENDDO
			ABITU=>CA
		ENDIF
	END SUBROUTINE
	SUBROUTINE CALCU(CV,SZ,M,FQ,Y,img)
		IMPLICIT NONE
		INTEGER(8) :: SZ(2),MODI(0:1)
		INTEGER :: is,M,i,j
		REAL(8) :: SUN,FQ(:),Y(:),img,PII,AP(1:M),BT(0:M),V(SZ(1)*SZ(2),0:1)
		REAL(8),POINTER :: CV(:)
		COMPLEX(8) :: X,Z
		V(:,0)=0D0
		V(:,1)=CV
		SUN=DOT_PRODUCT(V(:,1),V(:,1))
		V(:,1)=V(:,1)/SQRT(ABS(SUN))
		BT(0)=0D0
		PII=-1D0/PI
		DO is=1,M
			MODI(1)=MOD(is,2)
			MODI(0)=MOD(is-1,2)
			V(:,MODI(0))=-1*BT(is-1)*V(:,MODI(0))
			CALL HUBBARD(V(:,MODI(1)),V(:,MODI(0)),SZ)
			AP(is)=DOT_PRODUCT(V(:,MODI(1)),V(:,MODI(0)))
			V(:,MODI(0))=V(:,MODI(0))-AP(is)*V(:,MODI(1))
			SUN=DOT_PRODUCT(V(:,MODI(0)),V(:,MODI(0)))
			BT(is)=SQRT(ABS(SUN))
			IF(BT(is)<1E-8.AND.is>10) THEN
				WRITE(*,"(E17.6,I5)")BT(is),is
				EXIT
			ENDIF
			V(:,MODI(0))=V(:,MODI(0))/BT(is)
		ENDDO
		BT(0)=DOT_PRODUCT(CV,CV)
		DEALLOCATE(CV)
		DO i=1,SIZE(FQ)
			X=0D0
			Z=CMPLX(FQ(i),img)
			DO j=is,1,-1
				X=BT(j-1)**2/(Z-AP(j)-X)
			ENDDO
			Y(i)=Y(i)+DIMAG(X)*PII
		ENDDO
		DEALLOCATE(CA)
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
	REAL(8),ALLOCATABLE :: V2(:),V_LOW(:),E(:),LOWE(:)
	REAL(8),POINTER :: CV(:)
	INTEGER(8) :: i,j,k,SZ(2),TMP(2),NSLOW
	REAL(8),ALLOCATABLE :: FQ(:),Y(:)
	REAL(8) :: OMIN,OMAX
	INTEGER :: M=200,INFO,Mo
	LOGICAL :: FLAG
	OPEN(UNIT=10,FILE="../DATA/OUTPUT.DAT")
	CALL INIT_RANDOM_SEED()
	CALL INIT_LATTICE()
	ALLOCATE(V_LOW(1))
	ALLOCATE(E(M),LOWE(M))
	LOWE(1)=100
	DO NS=NT,0,-2
		CALL INITIALSTATE(AD,SZ(1),(NT-NS)/2)
		IF(NS==0) THEN
			SZ(2)=SZ(1)
			ALLOCATE(AU(SZ(2)))
			AU=AD
		ELSE
			CALL INITIALSTATE(AU,SZ(2),(NT+NS)/2)
		ENDIF
		ALLOCATE(V2(SZ(1)*SZ(2)))
		CALL RANDOM_NUMBER(V2)
		ABITD=>AD
		ABITU=>AU
		CALL LANC(V2,E,SZ,M,Mo,.TRUE.,INFO)
		IF(E(1)<LOWE(1)) THEN
			DEALLOCATE(LOWE,V_LOW)
			ALLOCATE(LOWE(Mo),V_LOW(SZ(1)*SZ(2)))
			LOWE=E(1:Mo)
			NSLOW=NS
			V_LOW=V2
		ENDIF
		DEALLOCATE(V2,AD,AU)
	ENDDO
	NS=NSLOW
	CALL INITIALSTATE(AD,SZ(1),(NT-NS)/2)
	IF(NS==0) THEN
		SZ(2)=SZ(1)
		ALLOCATE(AU(SZ(2)))
		AU=AD
	ELSE
		CALL INITIALSTATE(AU,SZ(2),(NT+NS)/2)
	ENDIF
	WRITE(10,"(E17.6)")LOWE
	WRITE(10,*)"----------------------------------------------"
	IF(NT==0) THEN
		V_LOW=1D0
	ENDIF
	ALLOCATE(FQ(1237),Y(1237))
	M=100
	DEALLOCATE(E)
	ALLOCATE(E(M))
	! ABITU=>AU
	! ABITD=>AD
	! ALLOCATE(V2(SZ(1)*SZ(2)))
	! V2=0D0
	! CALL HUBBARD(V_LOW,V2,SZ)
	! WRITE(10,"(E17.6)")V_LOW
	! WRITE(10,*)"----------------------------------------------"
	FLAG=.TRUE.
	DO k=0,3
		FLAG=.NOT.FLAG
		WRITE(*,*)FLAG,k,k/2
		Y=0D0
		TMP=SZ
		ABITU=>AU
		ABITD=>AD
		CALL CREATE(V_LOW,CV,k/10+1,k/2,SZ,FLAG)
		CALL LANC(CV,E,SZ,M,Mo,.FALSE.,INFO)
		OMIN=E(1)-0.5
		OMAX=E(Mo)+0.5
		DO i=1,SIZE(FQ)
			FQ(i)=OMIN+(OMAX-OMIN)/SIZE(FQ)*i
		ENDDO
		CALL CALCU(CV,SZ,M,FQ,Y,0.02D0)
		SZ=TMP
		WRITE(10,"(2E17.6)")(FQ(j)-LOWE(1),Y(j),j=1,SIZE(FQ))
	ENDDO
END