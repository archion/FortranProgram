! This is a exactly diagonalize program for 2*2 lattice with t' using hubbard model
! for more information see H. Q. LIN, PRB.35.3359
! Writen by Z. D. Yu in 2012.5
! 1---2
! |   |
! 3---4
! n=2: 3,5,6,9,10,12
! n=1: 1,2,4,8
! n=3: 7,11,13,14
MODULE GLOBAL
	IMPLICIT NONE
	INTEGER, PARAMETER :: N=36,DN=2
	REAL(8),SAVE :: V(0:N+1,N)=0D0,T0=1D0,U=5D0,TP=0
	INTEGER,SAVE :: TF(N),PBC(0:DN+1),nn(DN,DN)
END MODULE
PROGRAM MAIN
	USE GLOBAL
	REAL(8) :: EN(N),E
	INTEGER :: i,j,k=1,A1(6)=(/ 3,5,6,9,10,12 /),A2(6)=(/ 3,5,6,9,10,12 /)
	OPEN(UNIT=10,FILE='../DATA/doutput.dat')
	CALL INIT_RANDOM_SEED()
	PBC=(/1,(i,i=1,DN),DN/)
	DO i=1,DN
		DO j=1,DN
			nn(i,j)=k
			k=k+1
		ENDDO
	ENDDO
	DO j=1,SIZE(A2)
		DO k=1,SIZE(A1)
			TF((j-1)*6+k)=A2(j)*(2**(2*DN))+A1(k)
		ENDDO
	ENDDO
	DO i=1,25
		U=U+1
		CALL LANCZOS(EN)
		E=EN(1)
		WRITE(10,*)EN
		WRITE(10,"(F4.1,E13.5)")U,EN
	ENDDO
	WRITE(10,"(E13.5)")EN
	CLOSE(10)
END PROGRAM MAIN

SUBROUTINE LANCZOS(EN)
	USE GLOBAL
	IMPLICIT NONE
	REAL(8) :: NV(N)=0D0,AP(N)=0D0,BT(0:N)=0D0,CACHE,Z(N,N),WORK(2*N-2),EN(N)
	INTEGER :: i,j,k,l,INFO
	! V(1,4)=1
	! CALL RN(V(1,:),CACHE)
	! CALL HEMILTION(1,TF,NV)
	! WRITE(*,"(F5.1F5.1B10.8)")(V(1,i),NV(i),TF(i),i=1,N)
	CALL RANDOM_NUMBER(V(1,:))
	CALL RN(V(1,:),CACHE)
	! WRITE(10,"(36E9.2)")(V(k,:),k=0,3)
	! WRITE(*,*)V(1,:)
	DO i=1,N
		CALL HEMILTION(i,NV)
		AP(i)=DOT_PRODUCT(V(i,:),NV)
		! WRITE(10,"(2E11.4)")(V(i,j),NV(j),j=1,36)
		! WRITE(10,*)AP(i)
		V(i+1,:)=NV(:)-AP(i)*V(i,:)-BT(i-1)*V(i-1,:)
		IF(MOD(i,4)==0) THEN
			CALL RO(i+1)
		ENDIF
		CALL RN(V(i+1,:),BT(i))
		! WRITE(10,"(E11.3)")BT(i)
		IF(BT(i)<1E-5) THEN
			BT(i)=0
			! DO WHILE(.TRUE.)
				CALL RANDOM_NUMBER(V(i+1,:))
				CALL RO(i+1)
				CALL RN(V(i+1,:),CACHE)
				! IF(CACHE>1E-10) THEN
					! EXIT
				! ENDIF
			! ENDDO
		ENDIF
		! WRITE(10,*)DOT_PRODUCT(V(i+1,:),V(i-1,:)),DOT_PRODUCT(V(i+1,:),V(i-5,:)),BT(i)
		! WRITE(10,*)DOT_PRODUCT(V(i+1,:),V(i,:)),BT(i),DOT_PRODUCT(NV(:),V(i+1,:))
		! WRITE(10,"(36E9.2)")V(i+1,:)
	ENDDO
	! WRITE(10,"(E13.5,E13.5)")(AP(i),BT(i),i=1,N)
	CALL DSTEQR('I',N,AP,BT(1:N-1),Z,N,WORK,INFO)
	! IF(INFO==0) THEN
		! WRITE(10,"(E13.5)")AP
	! ENDIF
	EN=AP
	! CLOSE(10)
END SUBROUTINE LANCZOS


SUBROUTINE HEMILTION(p,NV)
	USE GLOBAL
	IMPLICIT NONE
	REAL(8) :: NV(N),CACHE
	INTEGER :: NP(4),i,j,k,ll,l=0,p,CTF,NOS,UOS
	NV=0
	DO i=1,N
		DO j=0,1
			DO k=0,DN*DN-1
				IF(BTEST(TF(i),j*4+k)) THEN
					CALL NNP(k+1,0,NP)
					DO ll=1,4
						CTF=TF(i)
						IF((.NOT.(BTEST(TF(i),j*4+NP(ll)-1))).AND.(NP(ll)/=k+1)) THEN
							CTF=IBSET(CTF,j*4+NP(ll)-1)
							CTF=IBCLR(CTF,j*4+k)
							CALL FD(CTF,l)
							NV(l)=NV(l)+T0*V(p,i)
						ENDIF
					ENDDO
					CALL NNP(k+1,1,NP)
					DO ll=1,4
						CTF=TF(i)
						IF((.NOT.(BTEST(TF(i),j*4+NP(ll)-1))).AND.(NP(ll)/=k+1)) THEN
							CTF=IBSET(CTF,j*4+NP(ll)-1)
							CTF=IBCLR(CTF,j*4+k)
							CALL FD(CTF,l)
							NV(l)=NV(l)+TP*V(p,i)
						ENDIF
					ENDDO
				ENDIF
			ENDDO
		ENDDO
		UOS=IAND(IBITS(TF(i),0,4),IBITS(TF(i),4,4))
		NOS=0
		DO j=0,DN*DN-1
			IF(BTEST(UOS,j)) THEN
				NOS=NOS+1
			ENDIF
		ENDDO
		NV(i)=NV(i)+V(p,i)*U*NOS !+V(p,i)*4*TP*3
	ENDDO
END SUBROUTINE HEMILTION

SUBROUTINE FD(CTF,l)
	USE GLOBAL, ONLY : N,TF
	IMPLICIT NONE
	INTEGER :: i,l,CTF
	DO i=1,N
		IF(TF(i)==CTF) THEN
			l=i
			EXIT
		ENDIF
	ENDDO
END SUBROUTINE FD

SUBROUTINE RN(V,X)
	USE GLOBAL, ONLY : N
	IMPLICIT NONE
	REAL(8) :: V(N),SUN=0,X
	INTEGER :: i
	SUN=DOT_PRODUCT(V,V)
	X=SQRT(ABS(SUN))
	V=V/X
END SUBROUTINE RN

SUBROUTINE RO(i)
	USE GLOBAL
	IMPLICIT NONE
	INTEGER :: i,j
	DO j=1,i-1
		V(i,:)=V(i,:)-DOT_PRODUCT(V(j,:),V(i,:))*V(j,:)
	ENDDO
END

SUBROUTINE NNP(ii,p,jj)
	USE GLOBAL, ONLY : PBC,nn,DN
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

