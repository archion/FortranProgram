PROGRAM MAIN
	IMPLICIT NONE
	REAL(8),PARAMETER :: PI=3.1415926
	INTEGER,PARAMETER :: N=8
	INTEGER :: i,j
	REAL(8) :: A(0:N-1)
	REAL(8) :: B(0:N-1)
	OPEN(UNIT=10,FILE="../DATA/OUTPUT.DAT")
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                     CHACK THE BIT-REVERSAL                !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! DO i=0,N-1
		! A(i)=i
	! ENDDO
	! CALL FFT(A,1,N)
	! WRITE(*,"(I6,B19.16,I6,B19.16)")(INT(A(i)),INT(A(i)),i,i,i=0,N-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                     CHACK THE FFT                         !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! A=(/0,25,2,3,0,3,2,1/)
	! CALL FFTREAL1_RP(A,N,1)
	! CALL FFTREAL1_RP(A,N,-1)
	! ! CALL FFT1(A,N)
	! WRITE(10,"(E11.3)")A
	B=(/1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
	CALL FFTCOS1_2(B,N,1)
	WRITE(10,"(E11.3)")B
	CALL FFTCOS1_2(B,N,-1)
	! CALL FFTCOS1_2(B,N,-1)
	! WRITE(10,"(E11.3)")B
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!      Obtain the whole FFT array in Replaces FFT subroutine        !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! WRITE(10,"(2E11.3)")B(0),0D0
	! WRITE(10,"(2E11.3)")B(2:N-1)
	! WRITE(10,"(2E11.3)")B(1),0D0
	! B(2:N-2:2)=B(N-2:2:-2)
	! B(3:N-1:2)=-1*B(N-1:3:-2)
	! WRITE(10,"(2E11.3)")B(2:N-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                 CHACK THE FFT of sine                     !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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















