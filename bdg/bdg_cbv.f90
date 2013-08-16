MODULE GLOBAL
	IMPLICIT NONE
	SAVE
	INTEGER,PARAMETER :: DN=24
	REAL(8),PARAMETER :: BT=1E5,DV=1,nf0=0.85,DU=2.44,T=1,TP=-0.25,eff=1,PI=3.14159265D0
	INTEGER :: PBC(0:DN+1)
END MODULE
PROGRAM MAIN
	USE GLOBAL
	IMPLICIT NONE
	CHARACTER(25) :: IT,FT,fmat(3)
	LOGICAL :: FLAGA,FLAGB,FLAG1,FLAG2,FLAG3
	 !GENERAL nd(:,1):spin down,nd(:,0):spin up
	REAL(8) :: nd(DN*DN,0:1),ndp(DN*DN,0:1),DT(DN*DN,4),DT1(DN*DN,4),DTF(DN*DN),RDOM(DN*DN),er(3)=1,sa,sb,sp,sp0,nf,CVG=1E-3,wide
	INTEGER :: i,j,k,l,n,rc,Im(5)=-1,NP(4)
	! CPE
	REAL(8) :: VA(2*DN*DN*13),hhn(0:2,2*DN*DN),a,b,EMAX,EMIN,AC_ba,Tn,An(1000)
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
	CALL RANDOM_NUMBER(RDOM)
	! RDOM=0.5D0	!Uniform initial value
	DO i=1,DN*DN
		DT1(i,:)=(/-0.2D0,0.2D0,-0.2D0,0.2D0/) !  +((/RDOM(i)-0.5D0,-RDOM(i)+0.5D0,RDOM(i)-0.5D0,-RDOM(i)+0.5D0/))/10D0
	ENDDO
	! ndp(:,0)=RDOM*nf0
	! ndp(:,1)=nf0-ndp(:,0)
	ndp(:,0)=nf0/2+(RDOM-0.5D0)/20D0
	ndp(:,1)=nf0-ndp(:,0)
	DO i=1,DN*DN
		DTF(i)=(DT1(i,1)+DT1(i,3)-DT1(i,2)-DT1(i,4))/4.0
	ENDDO
	WRITE(40,fmat(1))((DTF(DN*i+j),j=1,DN),i=0,DN-1)
	WRITE(40,fmat(1))(((-1)**(i+j)*0.5*(ndp(DN*i+j,0)-ndp(DN*i+j,1)),j=1,DN),i=0,DN-1)
	FLAG1=.TRUE.
	FLAG2=.TRUE.
	FLAG3=.TRUE.
	nd=0
	sp=0
	wide=0.1
	sp0=sp+wide
	! Im(1)=121
	! Im(2)=123
	! Im(1)=(DN/2-1)*DN+DN/2
	! Im(2)=Im(1)+1
	! Im(3)=Im(2)+DN
	! Im(4)=Im(1)+DN
	! Im(5)=Im(1)-1
	! DO i=1,DN*DN
		! IF(ANY(Im==i)) THEN
			! nd(i,0)=1
		! ENDIF
	! ENDDO
	! WRITE(20,fmat(1))((nd(DN*i+j,0),j=1,DN),i=0,DN-1)
EX:	DO WHILE(ANY(er>CVG))
		rc=rc+1
		DT=DT1
		nd=ndp
		nf=0
		sa=sp
		sb=sp
		FLAGA=.TRUE.
		FLAGB=.TRUE.
		WRITE(*,*)"NEXT"
		DO WHILE(ABS(nf-nf0)>CVG)
			sp=0.5*(sa+sb)
			WRITE(*,*)sp,"***"
			n=1
			IA(1)=1
			l=2
			DO i=1,DN*DN
				DO k=1,0,-1
					IF(ANY(Im==i)) THEN
						VA(n:n)=(-1)**k*(-1*DU*nd(i,k)+sp)+(-1)**MINLOC(ABS(Im-i))*eff
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
			a=(EMAX-EMIN)/2+0.1
			b=(EMAX+EMIN)/2
			! WRITE(*,*)EMAX,EMIN,INFO
			DO i=1,2*DN*DN
				DO j=IA(i),IA(i+1)-1
					IF(JA(j)==i) THEN
						VA(j)=VA(j)-b
						EXIT
					ENDIF
				ENDDO
			ENDDO
			VA=VA/a
			nf=0
			ndp=0
			AC_ba=ACOS(-b/a)
			!$OMP PARALLEL DO REDUCTION(+:nf) PRIVATE(hhn) NUM_THREADS(42)
			DO i=1,2*DN*DN
				hhn=0
				hhn(0,i)=1D0
				l=(MOD(i,2)-1)*(-1)
				ndp(i,l)=ndp(i,l)+a*hhn(0,i)*(1-AC_ba/PI)
				CALL CRSMV(VA,JA,IA,hhn(0,:),hhn(1,:),2*DN*DN)
				ndp(i,l)=ndp(i,l)-a*hhn(1,i)*2.0*SIN(AC_ba)/PI
				DO n=2,Cutoff
					Tn=-2.0*SIN(n*AC_ba)/(n*PI)
					n_MOD(1)=MOD(n-2,3)
					n_MOD(2)=MOD(n-1,3)
					n_MOD(3)=MOD(n,3)
					CALL CRSMV(VA,JA,IA,hhn(n_MOD(2),:),hhn(n_MOD(3),:),2*DN*DN)
					hhn(n_MOD(3),:)=2*hhn(n_MOD(3),:)-hhn(n_MOD(1),:)
					ndp(i,l)=ndp(i,l)+a*hhn(n_MOD(3),i)*Tn
					! nf=nf+(-1)**MOD(i,2)*hhn(n_MOD(3),i)*Tn
					! WRITE(20,*)hhn(n_MOD(3),i),i
				ENDDO
				ndp(i,1)=1-ndp(i,1)
				nf=nf+ndp(i,0)+ndp(i,1)
				! IF(ndp(i,0)<0.OR.ndp(i,1)<0) THEN
					! WRITE(*,*)ndp(i,0),ndp(i,1),"BELOW THE ZERO EXIT!!"
				! ENDIF
			ENDDO
			!$OMP END PARALLEL DO
			nf=1.0/(DN*DN)*nf
			WRITE(*,*)nf
			IF(ABS(nf-nf0)<=CVG) THEN
				EXIT
			ENDIF
			IF(nf<nf0) THEN
				FLAGA=.FALSE.
				sa=sp
				IF(FLAGA.OR.FLAGB) THEN
					sb=sp+wide
					sp=sb
				ENDIF
			ELSE
				FLAGB=.FALSE.
				sb=sp
				IF(FLAGA.OR.FLAGB) THEN
					sa=sp-wide
					sp=sa
				ENDIF
			ENDIF
		ENDDO
		wide=MAX(ABS(sp0-sp),0.01)
		sp0=sp
		ndp=0
		DT1=0
		!$OMP PARALLEL DO PRIVATE(hhn) NUM_THREADS(42)
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
				Tn=-2.0*SIN(n*AC_ba)/(n*PI)
				n_MOD(1)=MOD(n-2,3)
				n_MOD(2)=MOD(n-1,3)
				n_MOD(3)=MOD(n,3)
				CALL CRSMV(VA,JA,IA,hhn(n_MOD(2),:),hhn(n_MOD(3),:),2*DN*DN)
				hhn(n_MOD(3),:)=2*hhn(n_MOD(3),:)-hhn(n_MOD(1),:)
				DO j=1,4
					DT1(i,j)=DT1(i,j)+a*hhn(n_MOD(3),2*NP(j)-1)*Tn
				ENDDO
			ENDDO
		ENDDO
		!$OMP END PARALLEL DO
		! DO i=1,DN*DN
			! DO n=1,2*DN*DN
				! nd(i,0)=nd(i,0)+H(2*i-1,n)**2*f(E(n))
				! nd(i,1)=nd(i,1)+H(2*i,n)**2*(1-f(E(n)))
			! ENDDO
		! ENDDO
		REWIND(30)
		WRITE(30,fmat(1))((ndp(DN*i+j,0)+ndp(DN*i+j,1),j=1,DN),i=0,DN-1)
		! DT1=0
		! DO i=1,DN*DN
			! CALL NNP(i,.TRUE.,NP)
			! DO j=1,4
				! DO n=1,2*DN*DN
					! DT1(i,j)=DT1(i,j)+0.25*DV*(H(2*i-1,n)*H(2*NP(j),n)+H(2*NP(j)-1,n)*H(2*i,n))*TANH(0.5*BT*E(n))
				! ENDDO
			! ENDDO
			! WRITE(30,"(4E10.3,5I4)")DT1(i,:),i,NP
		! ENDDO
		er(1)=SUM(ABS(DT1-DT))/(DN*DN)
		er(2)=SUM(ABS(ndp(:,0)-nd(:,0)))/(DN*DN)
		er(3)=SUM(ABS(ndp(:,1)-nd(:,1)))/(DN*DN)
		WRITE(*,"(3E15.6,F15.12)")er,sa
		IF(MOD(rc,20)==0.AND.rc>=100) THEN
			DO i=1,DN*DN
				DTF(i)=(DT1(i,1)+DT1(i,3)-DT1(i,2)-DT1(i,4))/4.0
			ENDDO
			WRITE(50,*)"N=",rc
			WRITE(50,fmat(1))((DTF(DN*i+j),j=1,DN),i=0,DN-1)
		ENDIF
		IF(FLAG1.AND.ALL(er<CVG*100)) THEN
			DO i=1,DN*DN
				DTF(i)=(DT1(i,1)+DT1(i,3)-DT1(i,2)-DT1(i,4))/4.0
			ENDDO
			WRITE(10,fmat(1))((DTF(DN*i+j),j=1,DN),i=0,DN-1)
			WRITE(10,fmat(1))((ndp(DN*i+j,0)+ndp(DN*i+j,1),j=1,DN),i=0,DN-1)
			WRITE(10,fmat(1))(((-1)**(i+j)*0.5*(ndp(DN*i+j,0)-ndp(DN*i+j,1)),j=1,DN),i=0,DN-1)
			WRITE(10,*)"ERROR=",CVG*100
			FLAG1=.FALSE.
		ENDIF		
		IF(FLAG2.AND.ALL(er<CVG*50)) THEN
			DO i=1,DN*DN
				DTF(i)=(DT1(i,1)+DT1(i,3)-DT1(i,2)-DT1(i,4))/4.0
			ENDDO
			WRITE(10,fmat(1))((DTF(DN*i+j),j=1,DN),i=0,DN-1)
			WRITE(10,fmat(1))((ndp(DN*i+j,0)+ndp(DN*i+j,1),j=1,DN),i=0,DN-1)
			WRITE(10,fmat(1))(((-1)**(i+j)*0.5*(ndp(DN*i+j,0)-ndp(DN*i+j,1)),j=1,DN),i=0,DN-1)
			WRITE(10,*)"ERROR=",CVG*50
			FLAG2=.FALSE.
		ENDIF
		IF(FLAG3.AND.ALL(er<CVG*10)) THEN
			DO i=1,DN*DN
				DTF(i)=(DT1(i,1)+DT1(i,3)-DT1(i,2)-DT1(i,4))/4.0
			ENDDO
			WRITE(10,fmat(1))((DTF(DN*i+j),j=1,DN),i=0,DN-1)
			WRITE(10,fmat(1))((ndp(DN*i+j,0)+ndp(DN*i+j,1),j=1,DN),i=0,DN-1)
			WRITE(10,fmat(1))(((-1)**(i+j)*0.5*(ndp(DN*i+j,0)-ndp(DN*i+j,1)),j=1,DN),i=0,DN-1)
			WRITE(10,*)"ERROR=",CVG*10
			FLAG3=.FALSE.
		ENDIF
	ENDDO EX
	DO i=1,DN*DN
		DTF(i)=(DT(i,1)+DT(i,3)-DT(i,2)-DT(i,4))/4.0
	ENDDO
	CALL FDATE(FT)
	WRITE(10,fmat(1))((DTF(DN*i+j),j=1,DN),i=0,DN-1)
	WRITE(10,fmat(1))((nd(DN*i+j,0)+nd(DN*i+j,1),j=1,DN),i=0,DN-1)
	WRITE(10,fmat(1))(((-1)**(i+j)*0.5*(nd(DN*i+j,0)-nd(DN*i+j,1)),j=1,DN),i=0,DN-1)
	WRITE(10,*)"The beta is",BT,"V=",DV,"U=",DU,"t'=",TP,"heff=",eff,"nf=",nf0,"N=",DN
	WRITE(10,*)"INFO=",INFO
	WRITE(10,*)"THE TATAL RUNTIME IS FORM ",IT," TO ",FT
	! IF(INFO/=0) THEN
		! ERROR=MOD(INFO,2*DN*DN+1)
		! DO j=1,ERROR+1
			! WRITE(20,500)(H(i,j),i=1,ERROR+1)
		! ENDDO
	! ENDIF
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