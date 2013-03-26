MODULE GLOBAL
	IMPLICIT NONE
	SAVE
	REAL(8) :: T=1
	INTEGER, PARAMETER :: DN=4
	REAL(8), PARAMETER :: PI=3.14159
	INTEGER :: PBC(0:DN+1)
END MODULE
PROGRAM CPE
	USE GLOBAL
	IMPLICIT NONE
	INTEGER,PARAMETER :: SZ=1600,COF=1024
	INTEGER :: i,j,JA(DN*DN*5),IA(DN*DN+1),NJ(4),n,n_MOD(3),SHF(SZ),INFO
	REAL(8) :: VA(DN*DN*5),x,W,hhn(0:2,DN*DN),P(0:2),EMAX,EMIN,a,b,An(0:2*COF-1),DOS(0:2*COF-1),LD=0.01,TMP_SH
	OPEN(UNIT=10,FILE="../DATA/output.dat")
	PBC=(/DN,(i,i=1,DN),1/)
	CALL INIT_RANDOM_SEED()
	DO i=1,SZ
		SHF(i)=-i+100
	ENDDO
	CALL FISHER_YATES_SHUFFLE(SHF,SIZE(SHF))
	n=1
	IA(1)=1
	DO i=1,DN*DN
		! VA(n)=-4
		! JA(n)=i
		! n=n+1
		CALL NNP(i,.TRUE.,NJ)
		DO j=1,4
			VA(n)=T
			JA(n)=NJ(j)
			n=n+1
		ENDDO
		IA(i+1)=n
	ENDDO
	CALL LANMAX(VA,JA,IA,EMAX,EMIN,SIZE(IA)-1,INFO)
	WRITE(*,*)EMAX,EMIN,INFO
	a=(EMAX-EMIN)/2+0.1
	b=(EMAX+EMIN)/2
	WRITE(*,*)a,b
	DO i=1,DN*DN
		DO j=IA(i),IA(i+1)-1
			IF(JA(j)==i) THEN
				VA(j)=VA(j)-b
			ENDIF
		ENDDO
	ENDDO
	VA=VA/a
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                                                                                                                  !!!!!
!!!!!                                     Expanded by Tchebycheff polynomial                                           !!!!!
!!!!!                                                                                                                  !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! DO j=1,SZ
		! x=SHF(j)/(SZ-99D0)
		! IF(x*a+b<EMAX.AND.x*a+b>EMIN) THEN
			! DOS=0
			! W=1.0/SQRT(1-x**2)
			! !$OMP PARALLEL DO REDUCTION(+:DOS) PRIVATE(hhn) NUM_THREADS(42)
			! DO i=1,DN*DN
				! hhn=0
				! hhn(0,i)=1D0
				! DOS=DOS+hhn(0,i)/PI*W
				! CALL CRSMV(VA,JA,IA,hhn(0,:),hhn(1,:),DN*DN)
				! DOS=DOS+hhn(1,i)/(2.0*PI)*W*COS(ACOS(x))
				! DO n=2,COF
					! n_MOD(1)=MOD(n-2,3)
					! n_MOD(2)=MOD(n-1,3)
					! n_MOD(3)=MOD(n,3)
					! CALL CRSMV(VA,JA,IA,hhn(n_MOD(2),:),hhn(n_MOD(3),:),DN*DN)
					! hhn(n_MOD(3),:)=2*hhn(n_MOD(3),:)-hhn(n_MOD(1),:)
					! DOS=DOS+hhn(n_MOD(3),i)/(2.0*PI)*W*COS(n*ACOS(x))
				! ENDDO
			! ENDDO
			! !$OMP END PARALLEL DO
			! WRITE(10,*)x*a+b,DOS
		! ENDIF
	! ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                                                                                                                  !!!!!
!!!!!                                 Expanded by Tchebycheff polynomial using FFT                                     !!!!!
!!!!!                                                                                                                  !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DOS=0
	TMP_SH=SINH(LD)
	!$OMP PARALLEL DO REDUCTION(+:DOS) PRIVATE(hhn,An,n_MOD) NUM_THREADS(6) SCHEDULE(GUIDED)
	DO i=1,DN*DN
		hhn=0
		hhn(0,i)=1D0
		An(0)=hhn(0,i)
		CALL CRSMV(VA,JA,IA,hhn(0,:),hhn(1,:),DN*DN)
		An(1)=hhn(1,i) !*SINH(LD*(1.0-1.0/COF))/TMP_SH
		DO n=2,COF-1
			n_MOD(1)=MOD(n-2,3)
			n_MOD(2)=MOD(n-1,3)
			n_MOD(3)=MOD(n,3)
			CALL CRSMV(VA,JA,IA,hhn(n_MOD(2),:),hhn(n_MOD(3),:),DN*DN)
			hhn(n_MOD(3),:)=2*hhn(n_MOD(3),:)-hhn(n_MOD(1),:)
			An(n)=hhn(n_MOD(3),i) !*SINH(LD*(1.0-REAL(n)/COF))/TMP_SH
		ENDDO
		An(COF:2*COF-1)=0
		CALL FFTCOS1_2(An,2*COF,-1)
		DOS=DOS+An
	ENDDO
	!$OMP END PARALLEL DO
	WRITE(10,"(2E16.3)")(COS(PI*(i+0.5)/(2*COF-1)),DOS(i),i=0,2*COF-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                                                                                                                  !!!!!
!!!!!                                  Expanded by Legendre ploynomial                                                 !!!!!
!!!!!                                                                                                                  !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! DO j=1,SZ
		! x=SHF(j)/(SZ-99D0)
		! IF(x*(BW+0.5)<BW.AND.x*(BW+0.5)>-1*BW) THEN
			! DOS=0
			! !$OMP PARALLEL DO REDUCTION(+:DOS) PRIVATE(hhn,P) NUM_THREADS(42)
			! DO i=1,DN*DN
				! hhn=0
				! hhn(0,i)=1D0
				! P(0)=1
				! p(1)=x
				! DOS=DOS+hhn(0,i)*0.5*P(0)
				! CALL CRSMV(VA,JA,IA,hhn(0,:),hhn(1,:),DN*DN)
				! DOS=DOS+1.5*hhn(1,i)*P(1)
				! DO n=2,COF
					! n_MOD(1)=MOD(n-2,3)
					! n_MOD(2)=MOD(n-1,3)
					! n_MOD(3)=MOD(n,3)
					! CALL CRSMV(VA,JA,IA,hhn(n_MOD(2),:),hhn(n_MOD(3),:),DN*DN)
					! hhn(n_MOD(3),:)=(2*n+1)*1.0/(n+1)*hhn(n_MOD(3),:)-n*1.0/(n+1)*hhn(n_MOD(1),:)
					! P(n_MOD(3))=(2*n+1)*1.0/(n+1)*x*P(n_MOD(2))-n*1.0/(n+1)*P(n_MOD(1))
					! DOS=DOS+(2*n+1)/2.0*hhn(n_MOD(3),i)*P(n_MOD(3))
				! ENDDO
			! ENDDO
			! !$OMP END PARALLEL DO
			! WRITE(10,*)x*(BW+0.5),DOS
		! ENDIF
	! ENDDO
	CLOSE(10)
END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                                                                                                                  !!!!!
!!!!!                                  Find nearest and next-nearest site                                              !!!!!
!!!!!                                                                                                                  !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
SUBROUTINE NNP(ii,near,jj)
	USE GLOBAL
	IMPLICIT NONE
	INTEGER :: ii,i,j,jj(4)
	LOGICAL :: near
	j=PBC(MOD(ii,DN))
	i=(ii-j)/DN+1
	IF(near) THEN
		jj(1)=(i-1)*DN+PBC(j+1)
		jj(2)=(PBC(i+1)-1)*DN+j
		jj(3)=(i-1)*DN+PBC(j-1)
		jj(4)=(PBC(i-1)-1)*DN+j
	ELSE
		jj(1)=DN*(PBC(i+1)-1)+PBC(j+1)
		jj(2)=DN*(PBC(i+1)-1)+PBC(j-1)
		jj(3)=DN*(PBC(i-1)-1)+PBC(j+1)
		jj(4)=DN*(PBC(i-1)-1)+PBC(j-1)
	ENDIF
END SUBROUTINE NNP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                                                                                                                  !!!!!
!!!!!                                  Vector-Matrix product use CRS storage                                           !!!!!
!!!!!                                                                                                                  !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CRSMV(VA,JA,IA,X,Y,M)
	IMPLICIT NONE
	INTEGER :: i,j,M
	INTEGER :: JA(*),IA(*)
	REAL(8) :: VA(*),Y(*),X(*)
	DO i=1,M
		Y(i)=0
		DO j=IA(i),IA(i+1)-1
			Y(i)=Y(i)+VA(j)*X(JA(j))
		ENDDO
	ENDDO
END

SUBROUTINE FISHER_YATES_SHUFFLE(A,M)
	IMPLICIT NONE
	INTEGER :: A(*),i,j,TEM,M
	REAL(8) :: RAM
	DO i=M,1,-1
		CALL RANDOM_NUMBER(RAM)
		j=CEILING(RAM*i)
		TEM=A(i)
		A(i)=A(j)
		A(j)=TEM
	ENDDO
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
END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                                                                                                                  !!!!!
!!!!!                       Use Lanczos to detect maximum-minimum eigenvalue store in CRS formal                       !!!!!
!!!!!                                                                                                                  !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LANMAX(VA,JA,IA,EMAX,EMIN,M,INFO)
	IMPLICIT NONE
	INTEGER :: JA(*),IA(*),M,MODI(0:2),i,j,INFO
	REAL(8) :: VA(*),EMAX,EMIN,SUN,WORK
	REAL(8),ALLOCATABLE :: V(:,:),AP(:),BT(:),TMP(:),E(:),FE(:)
	ALLOCATE(V(M,0:2),AP(M),BT(0:M),TMP(0:M),E(M),FE(M))
	FE(1)=1000
	V(:,0)=0D0
	BT(0)=0D0
	CALL RANDOM_NUMBER(V(:,1))
	SUN=DOT_PRODUCT(V(:,1),V(:,1))
	V(:,1)=V(:,1)/SQRT(ABS(SUN))
	DO i=1,M
		MODI(2)=MOD(i+1,3)
		MODI(1)=MOD(i,3)
		MODI(0)=MOD(i-1,3)
		CALL CRSMV(VA,JA,IA,V(:,MODI(1)),V(:,MODI(2)),M)
		AP(i)=DOT_PRODUCT(V(:,MODI(1)),V(:,MODI(2)))
		E(1:i)=AP(1:i)
		TMP(1:i-1)=BT(1:i-1)
		CALL DSTEQR('N',i,E(1:i),TMP(1:i-1),WORK,i,WORK,INFO)
		! WRITE(10,*)E(1:i)
		IF(ABS(E(1)-FE(1))<1E-3.OR.INFO/=0) THEN
			EXIT
		ENDIF
		FE=E
		V(:,MODI(2))=V(:,MODI(2))-AP(i)*V(:,MODI(1))-BT(i-1)*V(:,MODI(0))
		SUN=DOT_PRODUCT(V(:,MODI(2)),V(:,MODI(2)))
		BT(i)=SQRT(ABS(SUN))
		V(:,MODI(2))=V(:,MODI(2))/BT(i)
	ENDDO
	EMAX=E(MIN(i,M))
	EMIN=E(1)
	DEALLOCATE(V,AP,BT,TMP,E,FE)
END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!       FFT FOR COMPLEX ARRAY IN ONE-DEMENSION              !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FFT1(A,n,SIG)
	IMPLICIT NONE
	REAL(8),PARAMETER :: PI=3.1415926
	INTEGER :: i,j,m,SIG,n
	COMPLEX(8) :: A(0:n-1),TMP,WP,W
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!          bit-reversal using Gold-Rader algorithm          !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	j=n/2
	DO i=1,n-2
		IF(j>i) THEN
			TMP=A(i)
			A(i)=A(j)
			A(j)=TMP
		ENDIF
		m=n/2
		DO WHILE(j>=m)
			j=j-m
			m=m/2
		ENDDO
		j=j+m
	ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                     Fast Fourier Transform using Danielson-Lanczos Lemma                          !!!!
!!!!     trigonometric recurrence relations using 5.5.6 in book "Numerical Recipes in FORTRAN 77"      !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	m=2
	DO WHILE(m<=n)
		WP=CMPLX(-2D0*SIN(SIG*PI/m)**2,SIN(SIG*2D0*PI/m))
		W=CMPLX(1D0,0D0)
		DO j=0,m/2-1
			DO i=j,N-1,m
				TMP=W*A(i+m/2)
				A(i+m/2)=A(i)-TMP
				A(i)=A(i)+TMP
			ENDDO
			W=W*WP+W
		ENDDO
		m=m*2
	ENDDO
	IF(SIG==-1) THEN
		A=A/n
	ENDIF
END




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!          FFT FOR REAL ARRAY IN ONE-DEMENSION              !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FFTREAL1(A,B,n)
	IMPLICIT NONE
	REAL(8),PARAMETER :: PI=3.1415926
	INTEGER :: i,m,n
	REAL(8) :: A(0:n-1)
	COMPLEX(8) :: B(0:n-1),H1(0:n/2-1),H2(0:n/2-1),WP,W
	B=CMPLX(A(0:n-2:2),A(1:n-1:2))
	m=n/2
	CALL FFT1(B(0:m-1),m,1)
	H1(0)=CMPLX(REAL(B(0)),0)
	H2(0)=CMPLX(AIMAG(B(0)),0)
	H1(1:m/2)=0.5*(B(1:m/2)+CONJG(B(m-1:m/2:-1)))
	H2(1:m/2)=CMPLX(0,-0.5)*(B(1:m/2)-CONJG(B(m-1:m/2:-1)))
	H1(m-1:m/2:-1)=CONJG(H1(1:m/2))
	H2(m-1:m/2:-1)=CONJG(H2(1:m/2))
	WP=CMPLX(-2D0*SIN(PI/n)**2,SIN(2D0*PI/n))
	W=CMPLX(1D0,0D0)
	DO i=0,m-1
		B(i)=H1(i)+W*H2(i)
		B(i+m)=H1(i)-W*H2(i)
		W=W*WP+W
	ENDDO
END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                      FFT FOR REAL ARRAY IN ONE-DEMENSION                            !!!!
!!!!               Using input real array to storage the output complex array            !!!!
!!!!     with even index refer to real parts and odd index refer to imaginary parts      !!!!
!!!!           and a(0) is 0th while a(1) is n/2th of the FFT of input data              !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FFTREAL1_RP(A,n,SIG)
	IMPLICIT NONE
	REAL(8),PARAMETER :: PI=3.1415926D0
	INTEGER :: i,m,n,SIG
	REAL(8) :: A(0:n-1)
	COMPLEX(8) :: B(0:n/2-1),WP,W,HE,HO,TMP
	m=n/2
	IF(SIG==1) THEN
		B=CMPLX(A(0:n-2:2),A(1:n-1:2))
		CALL FFT1(B,m,1)
		WP=CMPLX(-2D0*SIN(PI/n)**2,SIN(2D0*PI/n))
		W=CMPLX(1D0,0D0)
		W=W*WP+W
		A(0)=REAL(B(0))+AIMAG(B(0))
		A(1)=REAL(B(0))-AIMAG(B(0))
		DO i=1,m-1
			HE=0.5*(B(i)+CONJG(B(m-i)))
			HO=CMPLX(0,-0.5)*(B(i)-CONJG(B(m-i)))
			TMP=HE+W*HO
			A(2*i)=REAL(TMP)
			A(2*i+1)=AIMAG(TMP)
			W=W*WP+W
		ENDDO
	ELSE
		B(1:m-1)=CMPLX(A(2:n-2:2),A(3:n-1:2))
		WP=CMPLX(-2D0*SIN(-1D0*PI/n)**2,SIN(-2D0*PI/n))
		W=CMPLX(1D0,0D0)
		W=W*WP+W
		DO i=1,m/2
			HE=0.5*(B(i)+CONJG(B(m-i)))
			HO=CMPLX(0D0,0.5D0)*(B(i)-CONJG(B(m-i)))
			B(i)=HE+W*HO
			B(m-i)=CONJG(HE)-CONJG(W*HO)
			W=W*WP+W
		ENDDO
		B(0)=0.5*CMPLX(A(0)+A(1),A(0)-A(1))
		CALL FFT1(B,m,-1)
		DO i=0,m-1
			A(2*i)=REAL(B(i))
			A(2*i+1)=AIMAG(B(i))
		ENDDO
	ENDIF
END

SUBROUTINE FFTSIN1(A,n)
	IMPLICIT NONE
	REAL(8),PARAMETER :: PI=3.1415926D0
	INTEGER :: i,n
	REAL(8) :: A(0:n-1),TMP
	COMPLEX(8) :: WP,W
	WP=CMPLX(-2D0*SIN(0.5D0*PI/n)**2,SIN(PI/n))
	W=CMPLX(1,0)
	A(0)=0
	DO i=1,n/2
		W=W*WP+W
		TMP=0.5*(A(i)-A(n-i))
		A(i)=AIMAG(W)*(A(i)+A(n-i))+TMP
		A(n-i)=A(i)-2*TMP
	ENDDO
	CALL FFTREAL1_RP(A,n,1)
	A(1)=0.5*A(0)
	A(0)=0
	DO i=3,N-1,2
		TMP=A(i)
		A(i)=A(i-2)+A(i-1)
		A(i-1)=TMP
	ENDDO
END
SUBROUTINE FFTCOS1(A,n)
	IMPLICIT NONE
	REAL(8),PARAMETER :: PI=3.1415926D0
	INTEGER :: i,n
	REAL(8) :: A(0:n),SUN,TMP
	COMPLEX(8) :: WP,W
	WP=CMPLX(-2D0*SIN(0.5D0*PI/n)**2,SIN(PI/n))
	W=CMPLX(1,0)
	SUN=0.5*(A(0)-A(n))
	A(0)=0.5*(A(0)+A(n))
	DO i=1,n/2
		W=W*WP+W
		TMP=A(n-i)-A(i)
		A(i)=0.5*(A(i)+A(n-i))+AIMAG(W)*TMP
		A(n-i)=A(i)-2*AIMAG(W)*TMP
		SUN=SUN-TMP*REAL(W)
	ENDDO
	CALL FFTREAL1_RP(A(0:n-1),n,1)
	A(N)=A(1)
	A(1)=SUN
	DO i=3,N-1,2
		A(i)=A(i-2)+A(i)
	ENDDO
END
SUBROUTINE FFTCOS1_2(A,n,SIG)
	IMPLICIT NONE
	REAL(8),PARAMETER :: PI=3.1415926D0
	INTEGER :: i,n,m,SIG
	REAL(8) :: A(0:n-1),TMP,TMP0,TMP1
	COMPLEX(8) :: WP,W,WI,WPI
	m=n/2
	W=CMPLX(COS(0.5D0*PI/n),SIN(0.5D0*PI/n))
	WI=CMPLX(COS(0.5D0*PI),SIN(0.5D0*PI))
	WP=CMPLX(-2D0*SIN(0.5D0*PI/n)**2,SIN(PI/n))
	WPI=CMPLX(-2D0*SIN(-0.5D0*PI/n)**2,SIN(-1D0*PI/n))
	IF(SIG==1) THEN
		DO i=0,m-1
			TMP=AIMAG(W)*(A(n-i-1)-A(i))
			A(i)=0.5*(A(i)+A(n-i-1))+TMP
			A(n-i-1)=A(i)-2*TMP
			W=W*WP+W
		ENDDO
		CALL FFTREAL1_RP(A,n,SIG)
		TMP=A(n-1)
		A(n-1)=0.5*A(1)
		DO i=n-1,3,-2
			WI=WI*WPI+WI
			TMP0=A(i-1)
			TMP1=A(i-2)
			A(i-1)=REAL(WI)*A(i-1)-AIMAG(WI)*TMP
			A(i-2)=AIMAG(WI)*TMP0+REAL(WI)*TMP+A(i)
			TMP=TMP1
		ENDDO
	ELSE
		TMP0=A(N-1)
		DO i=n-1,3,-2
			WI=WI*WPI+WI
			A(i)=A(i-2)-A(i)
			TMP=A(i)*REAL(WI)-A(i-1)*AIMAG(WI)
			A(i-1)=A(i-1)*REAL(WI)+A(i)*AIMAG(WI)
			A(i)=TMP
		ENDDO
		A(1)=2*TMP0
		CALL FFTREAL1_RP(A,n,SIG)
		DO i=0,m-1
			TMP=A(i)+A(n-i-1)
			A(i)=(A(n-i-1)-A(i))/(2*AIMAG(W))
			A(n-i-1)=0.5*(TMP-A(i))
			A(i)=0.5*(TMP+A(i))
			W=W*WP+W
		ENDDO
	ENDIF
END