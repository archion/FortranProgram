! �Ż��Ķ��ַ�
! ���������x���ҳ�x
! by Zuodong Yu at 2012.12.18
! Ҫ�㣺
! 1.�ӳ�ʼֵ������һ������ȷ������ֵ�ķ�Χ
! 2.ȷ������ֵ��Χ�����ж��ַ�
! 3.�õ��Ĵ���ֵ������һ�εĳ�ʼֵ�����ʵ���������
! 4.����������δ���ֵ���ϴδ���ֵ��ͬ��ֱ���˳����ַ�ѭ��
PROGRAM BISECTION_METHOD
	IMPLICIT NONE
	INTEGER :: i,n
	REAL(8) :: sa,sb,sp,sp0,wide,CVG=1E-5,x
	LOGICAL :: FLAGA,FLAGB
	OPEN(UNIT=10,FILE="../DATA/FIND.DAT")
	OPEN(UNIT=20,FILE="../DATA/INPUT.DAT")
	sp=0 ! ��ʼֵ
	wide=0.1 ! ����
	sp0=sp+wide
	DO i=1,13
		READ(20,*)x
		WRITE(10,"(2G12.6)")"x=",x
		sa=sp
		sb=sp
		sp=0
		n=0
		FLAGA=.TRUE. ! ��ǣ���ʾû���ҵ�����ֵ��Χ
		FLAGB=.TRUE. ! ��ǣ���ʾû���ҵ�����ֵ��Χ
		! DO WHILE(ABS(sa-sb)>CVG.OR.ABS(sa-sb)<1E-8)
		DO WHILE(ABS(x-sp)>CVG)
			n=n+1
			! ! ����
			! CALL SLEEP(1)
			sp=0.5*(sa+sb)
			IF(ABS(x-sp)<CVG) THEN ! ����������δ���ֵ���ϴδ���ֵ��ͬ��ֱ���˳����ַ�ѭ��
				EXIT
			ENDIF
			IF(sp<x) THEN !����ֵ����ֵ̽�ұ�
				FLAGA=.FALSE.
				sa=sp ! �������ַ�
				IF(FLAGA.OR.FLAGB) THEN ! �鿴��ǣ����Ƿ��Ѿ��ҵ�����ֵ��Χ������ı���ҷ�Χ
					sb=sp+wide
					sp=sb
				ENDIF
			ELSE !����ֵ����ֵ̽�ұ�
				FLAGB=.FALSE.
				sb=sp ! �������ַ�
				IF(FLAGA.OR.FLAGB) THEN ! �鿴��ǣ����Ƿ��Ѿ��ҵ�����ֵ��Χ������ı���ҷ�Χ
					sa=sp-wide
					sp=sa
				ENDIF
			ENDIF
		ENDDO
		wide=MAX(ABS(sp0-sp),0.02) ! �������Ҳ���
		sp0=sp ! ��Ϊ�´εĳ�ʼֵ
		IF(ABS(sp-x)>CVG) THEN
			WRITE(10,*)"WRONG!!!!!!!!!!!!!!!!!!!!!"
		ENDIF
		WRITE(*,"(4G12.6/G66.2)")"sp=",sp,"n=",n,"******************************************************************"
	ENDDO
END