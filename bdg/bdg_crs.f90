MODULE GLOBAL
	IMPLICIT NONE
	SAVE
	INTEGER,PARAMETER :: DN=26,COF=1000
	REAL(8),PARAMETER :: DV=1D0,nf0=0.85D0,DU=2.44D0,T=1D0,TP=-0.25D0,eff=1D0,PI=3.141592653589793D0
	INTEGER :: PBC(0:DN+1)
END MODULE
PROGRAM MAIN
	USE GLOBAL
	IMPLICIT NONE
	CHARACTER(25) :: IT,FT,fmat(3)
	REAL(8) :: TIME(2)
	LOGICAL :: FLAGA,FLAGB,FLAG1,FLAG2,FLAG3
	 !GENERAL nd(:,1):spin down,nd(:,0):spin up
	REAL(8) :: nd(DN*DN,0:1),ndp(DN*DN,0:1),DT(DN*DN,4),DT1(DN*DN,4),DTF(DN*DN),RDOM(DN*DN),er(3)=1,sa,sb,sp,sp0,nf,CVG=1E-4,wide
	INTEGER :: i,ii,j,k,l,n,rc,Im(5)=-1,NP(4)
	! CPE
	REAL(8) :: VA(2*DN*DN*13),hhn(0:2,2*DN*DN),a,b,EMAX,EMIN,AC_ba,Tn(COF),LD=0.01,TMP_SH,gm=0.1D0
	INTEGER :: JA(2*DN*DN*13),IA(2*DN*DN+1),n_MOD(3),INFO
	OPEN(UNIT=10,FILE='../DATA/output_CRS.dat')
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
	CALL CPU_TIME(TIME(1))
	! CALL CPU_TIME(TIME(1))
	rc=0
	PBC=(/DN,(i,i=1,DN),1/)
	DT=0D0
	DT1=0D0
	CALL INIT_RANDOM_SEED()
	CALL RANDOM_NUMBER(RDOM)
	! RDOM=0.5D0	!Uniform initial value
	DO i=1,DN*DN
		DT1(i,:)=(/0.07D0,-0.07D0,0.07D0,-0.07D0/)+((/RDOM(i)-0.5D0,-RDOM(i)+0.5D0,RDOM(i)-0.5D0,-RDOM(i)+0.5D0/))/10D0
	ENDDO
	! ndp(:,0)=RDOM*nf0
	! ndp(:,1)=nf0-ndp(:,0)
	ndp(:,0)=nf0/2+(RDOM-0.5D0)/20D0
	ndp(:,1)=nf0-ndp(:,0)
	DO i=1,DN*DN
		DTF(i)=(DT1(i,1)+DT1(i,3)-DT1(i,2)-DT1(i,4))/4.0
	ENDDO
	! WRITE(40,fmat(1))((DTF(DN*i+j),j=1,DN),i=0,DN-1)
	! WRITE(40,fmat(1))(((-1)**(i+j)*0.5*(ndp(DN*i+j,0)-ndp(DN*i+j,1)),j=1,DN),i=0,DN-1)
	FLAG1=.TRUE.
	FLAG2=.TRUE.
	FLAG3=.TRUE.
	TMP_SH=SINH(LD)
	nd=0
	sp=0.14
	wide=0.2
	sp0=sp+wide
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
EX:	DO WHILE(ANY(er(1:2)>CVG).AND.(rc<200))
		rc=rc+1
		DT=DT1
		nd=ndp
		nf=0D0
		sa=sp
		sb=sp
		FLAGA=.TRUE.
		FLAGB=.TRUE.
		DO WHILE(ABS(nf-nf0)>CVG)
			sp=0.5D0*(sa+sb)
			! WRITE(*,*)sp,"***"
			n=1
			IA(1)=1
			l=2
			! WRITE(*,*)"CRS va ja ia start"
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
			! WRITE(*,*)"CRS va ja ia finish,LANMAX start"
			CALL LANMAX(VA,JA,IA,EMAX,EMIN,2*DN*DN,INFO)
			a=(EMAX-EMIN)/2D0+gm
			b=(EMAX+EMIN)/2D0
			! WRITE(*,*)"LANMAX finish,a=",a,"b=",b,"rescale va start"
			DO i=1,2*DN*DN
				DO j=IA(i),IA(i+1)-1
					IF(JA(j)==i) THEN
						VA(j)=VA(j)-b
						EXIT
					ENDIF
				ENDDO
			ENDDO
			VA=VA/a
			! WRITE(*,*)"rescale va finish, produce Tn start"
			ndp=0
			AC_ba=ACOS(-b/a)
			DO n=1,COF
				Tn(n)=-2.0D0*SIN(n*AC_ba)/(n*PI)*SINH(LD*(1.0D0-REAL(n,8)/COF))/TMP_SH
			ENDDO
			!$OMP PARALLEL DO REDUCTION(+:ndp) PRIVATE(hhn,n_MOD,l,ii) SCHEDULE(GUIDED)
			DO i=1,2*DN*DN
				hhn=0
				hhn(0,i)=1D0
				l=(MOD(i,2)-1)*(-1)
				ii=(i+1)/2
				ndp(ii,l)=ndp(ii,l)+hhn(0,i)*(1.0D0-AC_ba/PI)			
				CALL CRSMV(VA,JA,IA,hhn(0,:),hhn(1,:),2*DN*DN)
				ndp(ii,l)=ndp(ii,l)+hhn(1,i)*Tn(1)
				DO n=2,COF
					n_MOD(1)=MOD(n-2,3)
					n_MOD(2)=MOD(n-1,3)
					n_MOD(3)=MOD(n,3)
					CALL CRSMV(VA,JA,IA,hhn(n_MOD(2),:),hhn(n_MOD(3),:),2*DN*DN)
					hhn(n_MOD(3),:)=2D0*hhn(n_MOD(3),:)-hhn(n_MOD(1),:)
					ndp(ii,l)=ndp(ii,l)+hhn(n_MOD(3),i)*Tn(n)
				ENDDO
			ENDDO
			!$OMP END PARALLEL DO
			ndp(:,1)=1D0-ndp(:,1)
			nf=0D0
			DO i=1,DN*DN
				nf=nf+ndp(i,0)+ndp(i,1)
			ENDDO
			nf=1.0D0/(DN*DN)*nf
			WRITE(*,"(4E15.8)")nf,sp,a,b
			IF(ABS(nf-nf0)<=CVG) THEN
				EXIT
			ENDIF			
			IF(nf>2.OR.nf<0) THEN
				WRITE(*,*)"THE A MAYBE SMALL"
				gm=gm+0.1
				FLAGA=.TRUE.
				FLAGB=.TRUE.
				CYCLE
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
		wide=MAX(ABS(sp0-sp),10*CVG)
		sp0=sp
		DT1=0D0
		!$OMP PARALLEL DO REDUCTION(+:DT1) PRIVATE(hhn,n_MOD,NP) SCHEDULE(GUIDED)
		DO i=1,DN*DN
			CALL NNP(i,.TRUE.,NP)
			hhn=0
			hhn(0,2*i)=1D0
			DO j=1,4
				DT1(i,j)=DT1(i,j)+hhn(1,2*NP(j)-1)*(1-AC_ba/PI)
			ENDDO
			CALL CRSMV(VA,JA,IA,hhn(0,:),hhn(1,:),2*DN*DN)
			DO j=1,4
				DT1(i,j)=DT1(i,j)+hhn(1,2*NP(j)-1)*Tn(1)
			ENDDO
			DO n=2,COF
				n_MOD(1)=MOD(n-2,3)
				n_MOD(2)=MOD(n-1,3)
				n_MOD(3)=MOD(n,3)
				CALL CRSMV(VA,JA,IA,hhn(n_MOD(2),:),hhn(n_MOD(3),:),2*DN*DN)
				hhn(n_MOD(3),:)=2*hhn(n_MOD(3),:)-hhn(n_MOD(1),:)
				DO j=1,4
					DT1(i,j)=DT1(i,j)+hhn(n_MOD(3),2*NP(j)-1)*Tn(n)
				ENDDO
			ENDDO
		ENDDO
		! DO i=1,DN*DN
			! CALL NNP(i,.TRUE.,NP)
			! hhn=0
			! hhn(0,2*i-1)=1D0
			! DO j=1,4
				! DT1(i,j)=DT1(i,j)+hhn(1,2*NP(j))*(1-AC_ba/PI)
			! ENDDO
			! CALL CRSMV(VA,JA,IA,hhn(0,:),hhn(1,:),2*DN*DN)
			! DO j=1,4
				! DT1(i,j)=DT1(i,j)+hhn(1,2*NP(j))*Tn(1)
			! ENDDO
			! DO n=2,COF
				! n_MOD(1)=MOD(n-2,3)
				! n_MOD(2)=MOD(n-1,3)
				! n_MOD(3)=MOD(n,3)
				! CALL CRSMV(VA,JA,IA,hhn(n_MOD(2),:),hhn(n_MOD(3),:),2*DN*DN)
				! hhn(n_MOD(3),:)=2*hhn(n_MOD(3),:)-hhn(n_MOD(1),:)
				! DO j=1,4
					! DT1(i,j)=DT1(i,j)+hhn(n_MOD(3),2*NP(j))*Tn(n)
				! ENDDO
			! ENDDO
		! ENDDO
		!$OMP END PARALLEL DO
		DT1=-1D0*DT1*DV
		! REWIND(50)
		! WRITE(50,fmat(1))((DTF(DN*i+j),j=1,DN),i=0,DN-1)
		WRITE(*,*)"DT=",DT1(1,1)
		er(1)=SUM(ABS(DT1-DT))/REAL(DN*DN,8)
		! er(2)=SUM(ABS(ndp(:,0)+ndp(:,1)-nd(:,0)-nd(:,1)))/REAL(2*DN*DN,8)
		er(2)=SUM(ABS(ndp(:,1)-nd(:,1)))/REAL(DN*DN,8)
		WRITE(*,"(2E15.8)")er(1:2)
		! IF(FLAG1.AND.ALL(er<CVG*100)) THEN
			! DO i=1,DN*DN
				! DTF(i)=(DT1(i,1)+DT1(i,3)-DT1(i,2)-DT1(i,4))/4.0
			! ENDDO
			! WRITE(10,fmat(1))((DTF(DN*i+j),j=1,DN),i=0,DN-1)
			! WRITE(10,fmat(1))((nd(DN*i+j,0)+nd(DN*i+j,1),j=1,DN),i=0,DN-1)
			! WRITE(10,fmat(1))(((-1)**(i+j)*0.5*(nd(DN*i+j,0)-nd(DN*i+j,1)),j=1,DN),i=0,DN-1)
			! WRITE(10,*)"ERROR=",CVG*100
			! FLAG1=.FALSE.
		! ENDIF		
		! IF(FLAG2.AND.ALL(er<CVG*50)) THEN
			! DO i=1,DN*DN
				! DTF(i)=(DT1(i,1)+DT1(i,3)-DT1(i,2)-DT1(i,4))/4.0
			! ENDDO
			! WRITE(10,fmat(1))((DTF(DN*i+j),j=1,DN),i=0,DN-1)
			! WRITE(10,fmat(1))((nd(DN*i+j,0)+nd(DN*i+j,1),j=1,DN),i=0,DN-1)
			! WRITE(10,fmat(1))(((-1)**(i+j)*0.5*(nd(DN*i+j,0)-nd(DN*i+j,1)),j=1,DN),i=0,DN-1)
			! WRITE(10,*)"ERROR=",CVG*50
			! FLAG2=.FALSE.
		! ENDIF
		! IF(FLAG3.AND.ALL(er<CVG*10)) THEN
			! DO i=1,DN*DN
				! DTF(i)=(DT1(i,1)+DT1(i,3)-DT1(i,2)-DT1(i,4))/4.0
			! ENDDO
			! WRITE(10,fmat(1))((DTF(DN*i+j),j=1,DN),i=0,DN-1)
			! WRITE(10,fmat(1))((nd(DN*i+j,0)+nd(DN*i+j,1),j=1,DN),i=0,DN-1)
			! WRITE(10,fmat(1))(((-1)**(i+j)*0.5*(nd(DN*i+j,0)-nd(DN*i+j,1)),j=1,DN),i=0,DN-1)
			! WRITE(10,*)"ERROR=",CVG*10
			! FLAG3=.FALSE.
		! ENDIF
	ENDDO EX
	DTF=(DT1(:,1)+DT1(:,3)-DT1(:,2)-DT1(:,4))/4.0
	CALL FDATE(FT)
	CALL CPU_TIME(TIME(2))
	! CALL CPU_TIME(TIME(2))
	WRITE(10,fmat(1))((DTF(DN*i+j),j=1,DN),i=0,DN-1)
	WRITE(10,fmat(1))((nd(DN*i+j,0)+nd(DN*i+j,1),j=1,DN),i=0,DN-1)
	WRITE(10,fmat(1))(((-1)**(i+j)*0.5*(nd(DN*i+j,0)-nd(DN*i+j,1)),j=1,DN),i=0,DN-1)
	WRITE(10,*)"V=",DV,"U=",DU,"t'=",TP,"heff=",eff,"nf=",nf0,"N=",DN
	WRITE(10,*)"INFO=",INFO
	WRITE(10,*)"THE TATAL RUNTIME IS FORM ",IT," TO ",FT,"COST CUP TIME:",TIME(2)-TIME(1),"S"
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
		IF(ABS(E(1)-FE(1))<1E-3.AND.ABS(E(i)-FE(i-1))<1E-3.OR.INFO/=0) THEN
			WRITE(*,*)"INFO=",INFO
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