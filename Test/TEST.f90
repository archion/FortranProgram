! MODULE GLOBAL
	! INTEGER, PARAMETER :: M=10
! END MODULE
PROGRAM TEST
	USE GLOBAL
	IMPLICIT NONE
	INTEGER :: jj,INFO,i,j
	REAL(8) :: x1,A(10,10),A1(10,10),Y(10)=0,X(10)=(/ 1,5,6,9,8,5,-8,-3,5,9 /),WORK(30)
	DATA A /   1,  4,  0,  0,  0,  2,  0,  0,  0,6.5&
			,  4, -1,  0,  0,  0,  0,  0,  0, -8,  4&
			,  0,  0,  5,  6, -9,  0,  0, -4,  0,  0&
			,  0,  0,  6, -5,  0,  0,  0, -2, -5,  0&
			,  0,  0, -9,  0,  7,  8,  9,  0,  0,  0&
			,  2,  0,  0,  0,  8, -7,  0,  0,  0,  0&
			,  0,  0,  0,  0,  9,  0,  9,4.5,  0,  0&
			,  0,  0, -4, -2,  0,  0,4.5, -9,  0,  0&
			,  0, -8,  0, -5,  0,  0,  0,  0,  4,  5&
			,6.5,  4,  0,  0,  0,  0,  0,  0,  5, -4/
	OPEN(UNIT=10,FILE="../DATA/OUTPUT.DAT")
	! A(1,:)=(/ 4,5,0,0 /)
	! A(2,:)=(/ 5,7,0,6 /)
	! A(3,:)=(/ 0,0,0,2 /)
	! A(4,:)=(/ 0,6,2,0 /)
	A1=A
	x1=sum(a)
	write(*,*)sum(a)
	CALL DSYEV( "V", "U", 10, A, 10, Y, WORK, 30, INFO )
	
	! DO jj=1,2
		! CALL CRSMV(A,X,Y,10)
	! ENDDO
	WRITE(10,"(10F8.2)")(A(i,:),i=1,10)
	WRITE(10,"(F8.2)")Y
	WRITE(10,*)"INFO=",INFO
	x=0
	A1=0
	DO i=1,10
		DO j=1,10
			DO jj=1,10
				A1(i,j)=A1(i,j)+A(jj,i)*A(jj,j)
			ENDDO
		ENDDO
	ENDDO
	WRITE(10,"(10F8.2)")(A1(i,:),i=1,10)
END

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
	INTEGER :: i=0,j=0,k1=1,k2=1,M,JA(10*M+1),IA(M+1)
	REAL(8) :: VA(10*M)
	REAL(8) :: A(M,M),X(*),Y(*)
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

! SUBROUTINE CRSMV(A,X,Y)
	! USE GLOBAL
	! IMPLICIT NONE
	! INTEGER :: JA(5*N2+1),IA(N1+1)
	! INTEGER :: i,j,k1=1,k2=1
	! REAL(8) :: VA(5*N2)
	! REAL(8) :: A(N1,N2),X(N1),Y(N1)
	! DO i=1,N1
		! IA(k2)=K1
		! k2=k2+1
		! DO j=1,N2
			! IF(ABS(A(i,j))>=1E-5) THEN
				! VA(K1)=A(i,j)
				! JA(k1)=j
				! k1=k1+1
			! ENDIF
		! ENDDO
	! ENDDO
	! IA(k2)=k1
	! DO i=1,N1
		! Y(i)=0
		! DO j=IA(i),IA(i+1)-1
			! Y(i)=Y(i)+VA(j)*X(JA(j))
		! ENDDO
	! ENDDO
! END

! FOR SYMMETRY MATRIX

! SUBROUTINE CRSMV(A,X,Y)
	! USE GLOBAL
	! IMPLICIT NONE
	! INTEGER :: JA(4*N2+1),IA(N1+1)
	! INTEGER :: i,j,k,k1=1,k2=1
	! REAL(8) :: VA(4*N2)
	! REAL(8) :: A(N1,N2),X(N1),Y(N1)
	! DO i=1,N1
		! IA(k2)=K1
		! k2=k2+1
		! DO j=i,N2
			! IF(ABS(A(i,j))>=1E-8) THEN
				! VA(K1)=A(i,j)
				! JA(k1)=j
				! k1=k1+1
			! ENDIF
		! ENDDO
	! ENDDO
	! IA(k2)=K1
	! DO i=1,N1
		! Y(i)=0
		! DO j=IA(i),IA(i+1)-1
			! Y(i)=Y(i)+VA(j)*X(JA(j))
		! ENDDO
		! DO k=1,IA(i)-1
			! IF(JA(k)==i) THEN
				! DO j=1,i
					! IF(IA(j)-1>=k) THEN
						! Y(i)=Y(i)+VA(k)*X(j-1)
						! EXIT
					! ENDIF
				! ENDDO
			! ENDIF
		! ENDDO
	! ENDDO
! END

! CALL MKL FUNCTION

! SUBROUTINE CRSMV(A,X,Y,M,N)
	! IMPLICIT NONE
	! INTEGER :: i,JOB(6),INFO,M,N
	! REAL(8) :: A(M,N),X(*),Y(*)
	! INTEGER, ALLOCATABLE :: JA(:),IA(:)
	! REAL(8), ALLOCATABLE :: ACSR(:)
	! JOB=(/0,1,1,2,M*N,1/)
	! ALLOCATE(ACSR(M*N),JA(M*N),IA(M+1))
	! CALL MKL_DDNSCSR(JOB, M, N, A, M, ACSR, JA, IA, INFO)
	! CALL MKL_DCSRGEMV("N", M, ACSR, IA, JA, X, Y)
! END
! SUBROUTINE CRSMV(A,X,Y)
	! USE GLOBAL
	! IMPLICIT NONE
	! INTEGER :: i,JOB(6),INFO
	! REAL(8) :: A(M,M),X(*),Y(*)
	! INTEGER, ALLOCATABLE :: JA(:),IA(:)
	! REAL(8), ALLOCATABLE :: ACSR(:)
	! JOB=(/0,1,1,2,M*M,1/)
	! ALLOCATE(ACSR(M*M),JA(M*M),IA(M+1))
	! CALL MKL_DDNSCSR(JOB, M, M, A, M, ACSR, JA, IA, INFO)
	! CALL MKL_DCSRGEMV("N", M, ACSR, IA, JA, X, Y)
! END