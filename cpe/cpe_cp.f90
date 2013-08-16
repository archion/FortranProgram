MODULE GLOBAL
	IMPLICIT NONE
	SAVE
	INTEGER,PARAMETER :: DN=24
	REAL(8),PARAMETER :: BT=1E5,DV=1,nf0=0.85,DU=2.44D0,T=1,TP=-0.25D0,eff=0.3,PI=3.14159265D0
	INTEGER :: PBC(0:DN+1)
END MODULE
PROGRAM MAIN
	USE GLOBAL
	IMPLICIT NONE
	CHARACTER(25) :: IT,FT,fmat(3)
	LOGICAL :: FLAGA,FLAGB,FLAG1,FLAG2,FLAG3
	 !GENERAL nd(:,1):spin down,nd(:,0):spin up
	REAL(8) :: nd(DN*DN,0:1),ndp(DN*DN,0:1),DT(DN*DN,4),DT1(DN*DN,4),DTF(DN*DN),RDOM(DN*DN),er(3)=1,sa,sb,sp,sp0,nf,CVG=1E-3,wide
	INTEGER :: i,j,k,l,n,rc,Im(1)=-1,NP(4),SHF(DN*DN)
	! CPE
	REAL(8) :: VA(2*DN*DN*13),hhn(0:2,2*DN*DN),a,b,EMAX,EMIN,AC_ba,Tn
	INTEGER :: JA(2*DN*DN*13),IA(2*DN*DN+1),Cutoff=1000,n_MOD(3),INFO
	! LAPACK
	! REAL(8) :: E(2*DN*DN),H(2*DN*DN,2*DN*DN),WORK(3*2*DN*DN)
	! REAL(8) :: f,En
	! INTEGER :: INFO=0,ERROR=0
	! f(En)=1/(1+EXP(BT*En))
	OPEN(UNIT=10,FILE='../DATA/output.dat')
	OPEN(UNIT=20,FILE='../DATA/ERROR.dat')
	OPEN(UNIT=30,FILE='../DATA/rundata.dat')
	OPEN(UNIT=40,FILE='../DATA/initial.dat')
	OPEN(UNIT=50,FILE='../DATA/gap.dat')
	OPEN(UNIT=60,FILE='../DATA/check.dat')
	WRITE(fmat(1),*)"(",DN,"(E15.6))"
	WRITE(fmat(2),*)"(",DN*DN,"(E15.6))"
	WRITE(fmat(3),*)"(",2*DN*DN,"(E15.6))"
200	FORMAT(24(E15.6))
300	FORMAT(625(E15.6))
500	FORMAT(100000E10.3)
	CALL FDATE(IT)
	rc=0
	PBC=(/DN,(i,i=1,DN),1/)
	DT=0D0
	DT1=0D0
	CALL INIT_RANDOM_SEED()
	DO i=1,DN*DN
		SHF(i)=i
	ENDDO
	CALL FISHER_YATES_SHUFFLE(SHF,SIZE(SHF))
	! Im=SHF(1)
	DO i=1,DN*DN
		DT1(i,:)=(/0.02D0,-0.02D0,0.02D0,-0.02D0/)
	ENDDO
	nd=nf0/2
	sp=0.142
EX:	DO WHILE(ABS(er(1))>CVG)
		rc=rc+1
		DT=DT1
		n=1
		IA(1)=1
		l=2
		DO i=1,DN*DN
			DO k=1,0,-1
				IF(ANY(Im==i)) THEN
					VA(n:n)=(-1)**k*(-1*DU*nd(i,k)+sp+eff)
					JA(n)=2*i-k
					n=n+1
				ELSE
					VA(n)=(-1)**k*(-1*DU*nd(i,k)+sp)
					JA(n)=2*i-k
					n=n+1
				ENDIF
				CALL NNP(i,.TRUE.,NP)
				DO j=1,4
					VA(n)=(-1)**k*T
					JA(n)=2*NP(j)-k
					n=n+1
					VA(n)=DT(i,j)
					JA(n)=2*NP(j)+k-1
					n=n+1
				ENDDO
				CALL NNP(i,.FALSE.,NP)
				DO j=1,4
					VA(n)=(-1)**k*TP
					JA(n)=2*NP(j)-k
					n=n+1
				ENDDO
				IA(l)=n
				l=l+1
			ENDDO
		ENDDO
		CALL LANMAX(VA,JA,IA,EMAX,EMIN,2*DN*DN,INFO)
		a=(EMAX-EMIN)/2+0.3
		b=(EMAX+EMIN)/2
		WRITE(*,*)EMAX,EMIN,INFO
		DO i=1,2*DN*DN
			DO j=IA(i),IA(i+1)-1
				IF(JA(j)==i) THEN
					VA(j)=VA(j)-b
					EXIT
				ENDIF
			ENDDO
		ENDDO
		VA=VA/a
		DT1=0
		AC_ba=ACOS(-b/a)
		DO i=1,DN*DN
			CALL NNP(i,.TRUE.,NP)
			hhn=0
			hhn(0,2*i)=1D0
			DO j=1,4
				DT1(i,j)=DT1(i,j)+a*hhn(1,2*NP(j)-1)*(1-AC_ba/PI)
			ENDDO
			CALL CRSMV(VA,JA,IA,hhn(0,:),hhn(1,:),2*DN*DN)
			DO j=1,4
				DT1(i,j)=DT1(i,j)-a*hhn(1,2*NP(j)-1)*2.0*SIN(AC_ba)/PI
			ENDDO
			DO n=2,Cutoff
				Tn=-2.0*a*SIN(n*AC_ba)/(n*PI)
				n_MOD(1)=MOD(n-2,3)
				n_MOD(2)=MOD(n-1,3)
				n_MOD(3)=MOD(n,3)
				CALL CRSMV(VA,JA,IA,hhn(n_MOD(2),:),hhn(n_MOD(3),:),2*DN*DN)
				hhn(n_MOD(3),:)=2*hhn(n_MOD(3),:)-hhn(n_MOD(1),:)
				DO j=1,4
					DT1(i,j)=DT1(i,j)+hhn(n_MOD(3),2*NP(j)-1)*Tn
				ENDDO
			ENDDO
		ENDDO
		DT1(i,j)=DV*DT1(i,j)
		er(1)=SUM(ABS(DT1-DT))/(DN*DN)
		REWIND(50)
		DO i=1,DN*DN
			DTF(i)=(DT1(i,1)+DT1(i,3)-DT1(i,2)-DT1(i,4))/4.0
		ENDDO
		WRITE(50,fmat(1))((DTF(DN*i+j),j=1,DN),i=0,DN-1)
	ENDDO EX
	DO i=1,DN*DN
		DTF(i)=(DT(i,1)+DT(i,3)-DT(i,2)-DT(i,4))/4.0
	ENDDO
	CALL FDATE(FT)
	WRITE(10,fmat(1))((DTF(DN*i+j),j=1,DN),i=0,DN-1)
	WRITE(10,fmat(1))((nd(DN*i+j,0)+nd(DN*i+j,1),j=1,DN),i=0,DN-1)
	WRITE(10,fmat(1))(((-1)**(i+j)*0.5*(nd(DN*i+j,0)-nd(DN*i+j,1)),j=1,DN),i=0,DN-1)
	CLOSE(10)
	CLOSE(20)
	CLOSE(30)
	CLOSE(40)
	CLOSE(50)
END PROGRAM MAIN

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