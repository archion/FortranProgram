MODULE GLOBAL
	IMPLICIT NONE
	SAVE
	REAL(8) :: BT=1E5,DV=1,nf0=0.85,DU=2.44,T=1,TP=-0.25,eff=1
	INTEGER, PARAMETER :: DN=26
	INTEGER :: nn(DN,DN),PBC(0:DN+1)
END MODULE
PROGRAM MAIN
	USE GLOBAL
	IMPLICIT NONE
	CHARACTER(25) :: IT,FT,fmat(2)
	REAL(8) :: TIME(2)
	INTEGER :: i,j,k,n,rc,Im(5)=-1,INFO=0,NP(4),ERROR=0
	REAL(8) :: n1(DN*DN),n2(DN*DN),n1p(DN*DN),n2p(DN*DN),E(2*DN*DN),H(2*DN*DN,2*DN*DN),DT(DN*DN,4),&
		DT1(DN*DN,4),DTF(DN*DN),WORK(3*2*DN*DN),RDOM(DN*DN),er(3)=1
	REAL(8) :: f,En,sa,sb,sp,sp0,nf,CVG=1E-4,wide
	LOGICAL :: FLAGA,FLAGB,FLAG1,FLAG2,FLAG3
	f(En)=1/(1+EXP(BT*En))
	OPEN(UNIT=10,FILE='../DATA/output.dat')
	OPEN(UNIT=20,FILE='../DATA/ERROR.dat')
	OPEN(UNIT=30,FILE='../DATA/rundata.dat')
	OPEN(UNIT=40,FILE='../DATA/initial.dat')
	OPEN(UNIT=50,FILE='../DATA/gap.dat')
	WRITE(fmat(1),*)"(",DN,"(E15.6))"
	WRITE(fmat(2),*)"(",DN*DN,"(E15.6))"
200	FORMAT(24(E15.6))
300	FORMAT(625(E15.6))
500	FORMAT(100000E10.3)
	CALL FDATE(IT)
	rc=0
	k=1
	PBC=(/DN,(i,i=1,DN),1/)
	DO i=1,DN
		DO j=1,DN
			nn(i,j)=k
			k=k+1
		ENDDO
	ENDDO
	DT=0D0
	DT1=0D0
	CALL INIT_RANDOM_SEED()
	CALL RANDOM_NUMBER(RDOM)
	! RDOM=0.5D0	!Uniform initial value
	DO i=1,DN*DN
		DT1(i,:)=(/0.07D0,-0.07D0,0.07D0,-0.07D0/)+((/RDOM(i)-0.5D0,-RDOM(i)+0.5D0,RDOM(i)-0.5D0,-RDOM(i)+0.5D0/))/10D0
	ENDDO
	! n1p=RDOM*nf0
	! n2p=nf0-n1p
	n1p=nf0/2+(RDOM-0.5D0)/20D0
	n2p=nf0-n1p
	! call random_number(n1p)
	! call random_number(n2p)
	DO i=1,DN*DN
		DTF(i)=(DT1(i,1)+DT1(i,3)-DT1(i,2)-DT1(i,4))/4.0
	ENDDO
	! WRITE(40,fmat(1))((DTF(DN*i+j),j=1,DN),i=0,DN-1)
	! WRITE(40,fmat(1))(((-1)**(i+j)*0.5*(n1p(DN*i+j)-n2p(DN*i+j)),j=1,DN),i=0,DN-1)
	FLAG1=.TRUE.
	FLAG2=.TRUE.
	FLAG3=.TRUE.
	n1=0
	n2=0
	sp=0.14
	wide=0.2
	sp0=sp+wide
	! Im(1)=121
	! Im(2)=123
	Im(1)=(DN/2-1)*DN+DN/2
	! Im(2)=Im(1)+1
	! Im(3)=Im(2)+DN
	! Im(4)=Im(1)+DN
	! Im(5)=Im(1)-1
	! DO i=1,DN*DN
		! IF(ANY(Im==i)) THEN
			! n1(i)=1
		! ENDIF
	! ENDDO
	! WRITE(20,fmat(1))((n1(DN*i+j),j=1,DN),i=0,DN-1)
EX:	DO WHILE(ANY(er>CVG))
		rc=rc+1
		! IF(er(1)>CVG) THEN
			DT=DT1
		! ENDIF
		! IF(ANY(er(2:3)>CVG)) THEN
			n1=n1p
			n2=n2p
		! ENDIF
		nf=0
		sa=sp
		sb=sp
		FLAGA=.TRUE.
		FLAGB=.TRUE.
		! DO WHILE(ABS(nf-nf0)>1E-3)
		DO WHILE(ABS(nf-nf0)>CVG)
			CALL CPU_TIME(TIME(1))
			sp=0.5*(sa+sb)
			H=0D0
			DO i=1,DN*DN
				! nearest neighbor pairs
				CALL NNP(i,.TRUE.,NP)
				DO j=1,4
					H(2*i-1:2*i,2*NP(j)-1:2*NP(j))=H(2*i-1:2*i,2*NP(j)-1:2*NP(j))+RESHAPE((/-T,DT(i,j),DT(i,j),&
						T/),(/2,2/))
				ENDDO
				! next-nearest neighbor pairs
				CALL NNP(i,.FALSE.,NP)
				DO j=1,4
					H(2*i-1:2*i,2*NP(j)-1:2*NP(j))=H(2*i-1:2*i,2*NP(j)-1:2*NP(j))+RESHAPE((/-TP,0.0D0,0.0D0,TP/),(/2,2/))
				ENDDO
				! on site
				H(2*i-1:2*i,2*i-1:2*i)=H(2*i-1:2*i,2*i-1:2*i)+RESHAPE((/DU*n2(i)-sp,0.0D0,0.0D0,-DU*n1(i)+sp/),&
					(/2,2/))
				IF(ANY(Im==i)) THEN
					H(2*i-1:2*i,2*i-1:2*i)=H(2*i-1:2*i,2*i-1:2*i)+RESHAPE((/(-1)**MINLOC(ABS(Im-i))*eff,0.0D0,&
						0.0D0,(-1)**MINLOC(ABS(Im-i))*eff/),(/2,2/))
				ENDIF
			ENDDO
			CALL DSYEV( "V", "U", 2*DN*DN, H, 2*DN*DN, E, WORK, 3*2*DN*DN, INFO )
			IF(INFO/=0) THEN
				EXIT EX
			ENDIF
			nf=0
			DO i=1,DN*DN
				DO n=1,2*DN*DN
					nf=nf+(H(2*i-1,n)**2*f(E(n))+H(2*i,n)**2*(1-f(E(n))))
				ENDDO
			ENDDO
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
			CALL CPU_TIME(TIME(2))
			WRITE(10,*)"COST CUP TIME:",TIME(2)-TIME(1),"S"
		ENDDO
		wide=MAX(ABS(sp0-sp),10*CVG)
		sp0=sp
		n1p=0
		n2p=0
		DO i=1,DN*DN
			DO n=1,2*DN*DN
				n1p(i)=n1p(i)+H(2*i-1,n)**2*f(E(n))
				n2p(i)=n2p(i)+H(2*i,n)**2*(1-f(E(n)))
			ENDDO
		ENDDO
		! REWIND(30)
		! WRITE(30,fmat(1))((n1p(DN*i+j)+n2p(DN*i+j),j=1,DN),i=0,DN-1)
		DT1=0
		DO i=1,DN*DN
			CALL NNP(i,.TRUE.,NP)
			DO j=1,4
				DO n=1,2*DN*DN
					DT1(i,j)=DT1(i,j)+0.25*DV*(H(2*i-1,n)*H(2*NP(j),n)+H(2*NP(j)-1,n)*H(2*i,n))*TANH(0.5*BT*E(n))
				ENDDO
			ENDDO
			WRITE(30,"(4E10.3,5I4)")DT1(i,:),i,NP
		ENDDO
		WRITE(*,*)DT1(1,1)
		er(1)=SUM(ABS(DT1-DT))/(DN*DN)
		er(2)=SUM(ABS(n1p-n1))/(DN*DN)
		er(3)=SUM(ABS(n2p-n2))/(DN*DN)
		WRITE(*,"(3E15.6,F15.12)")er,sa
		! IF(MOD(rc,20)==0.AND.rc>=100) THEN
			! DO i=1,DN*DN
				! DTF(i)=(DT1(i,1)+DT1(i,3)-DT1(i,2)-DT1(i,4))/4.0
			! ENDDO
			! WRITE(50,*)"N=",rc
			! WRITE(50,fmat(1))((DTF(DN*i+j),j=1,DN),i=0,DN-1)
		! ENDIF
		! IF(FLAG1.AND.ALL(er<CVG*100)) THEN
			! DO i=1,DN*DN
				! DTF(i)=(DT1(i,1)+DT1(i,3)-DT1(i,2)-DT1(i,4))/4.0
			! ENDDO
			! WRITE(10,fmat(1))((DTF(DN*i+j),j=1,DN),i=0,DN-1)
			! WRITE(10,fmat(1))((n1p(DN*i+j)+n2p(DN*i+j),j=1,DN),i=0,DN-1)
			! WRITE(10,fmat(1))(((-1)**(i+j)*0.5*(n1p(DN*i+j)-n2p(DN*i+j)),j=1,DN),i=0,DN-1)
			! WRITE(10,*)"ERROR=",CVG*100
			! FLAG1=.FALSE.
		! ENDIF		
		! IF(FLAG2.AND.ALL(er<CVG*50)) THEN
			! DO i=1,DN*DN
				! DTF(i)=(DT1(i,1)+DT1(i,3)-DT1(i,2)-DT1(i,4))/4.0
			! ENDDO
			! WRITE(10,fmat(1))((DTF(DN*i+j),j=1,DN),i=0,DN-1)
			! WRITE(10,fmat(1))((n1p(DN*i+j)+n2p(DN*i+j),j=1,DN),i=0,DN-1)
			! WRITE(10,fmat(1))(((-1)**(i+j)*0.5*(n1p(DN*i+j)-n2p(DN*i+j)),j=1,DN),i=0,DN-1)
			! WRITE(10,*)"ERROR=",CVG*50
			! FLAG2=.FALSE.
		! ENDIF
		! IF(FLAG3.AND.ALL(er<CVG*10)) THEN
			! DO i=1,DN*DN
				! DTF(i)=(DT1(i,1)+DT1(i,3)-DT1(i,2)-DT1(i,4))/4.0
			! ENDDO
			! WRITE(10,fmat(1))((DTF(DN*i+j),j=1,DN),i=0,DN-1)
			! WRITE(10,fmat(1))((n1p(DN*i+j)+n2p(DN*i+j),j=1,DN),i=0,DN-1)
			! WRITE(10,fmat(1))(((-1)**(i+j)*0.5*(n1p(DN*i+j)-n2p(DN*i+j)),j=1,DN),i=0,DN-1)
			! WRITE(10,*)"ERROR=",CVG*10
			! FLAG3=.FALSE.
		! ENDIF
	ENDDO EX
	DO i=1,DN*DN
		DTF(i)=(DT(i,1)+DT(i,3)-DT(i,2)-DT(i,4))/4.0
	ENDDO
	CALL FDATE(FT)
	WRITE(10,fmat(1))((DTF(DN*i+j),j=1,DN),i=0,DN-1)
	WRITE(10,fmat(1))((n1(DN*i+j)+n2(DN*i+j),j=1,DN),i=0,DN-1)
	WRITE(10,fmat(1))(((-1)**(i+j)*0.5*(n1(DN*i+j)-n2(DN*i+j)),j=1,DN),i=0,DN-1)
	WRITE(10,*)"The beta is",BT,"V=",DV,"U=",DU,"t'=",TP,"heff=",eff,"nf=",nf0,"N=",DN
	WRITE(10,*)"ERROR=",ERROR,"INFO=",INFO,"rc=",rc
	WRITE(10,*)"THE TATAL RUNTIME IS FORM ",IT," TO ",FT,"COST CUP TIME:",TIME(2)-TIME(1),"S"
	IF(INFO/=0) THEN
		ERROR=MOD(INFO,2*DN*DN+1)
		DO j=1,ERROR+1
			WRITE(20,500)(H(i,j),i=1,ERROR+1)
		ENDDO
	ENDIF
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

	