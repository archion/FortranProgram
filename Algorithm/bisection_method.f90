! 优化的二分法
! 对于任意的x，找出x
! by Zuodong Yu at 2012.12.18
! 要点：
! 1.从初始值出发以一定步长确定待求值的范围
! 2.确定待求值范围并进行二分法
! 3.得到的待求值进行下一次的初始值，并适当调整步长
! 4.如果遇到本次待求值与上次待求值相同则直接退出二分法循环
PROGRAM BISECTION_METHOD
	IMPLICIT NONE
	INTEGER :: i,n
	REAL(8) :: sa,sb,sp,sp0,wide,CVG=1E-5,x
	LOGICAL :: FLAGA,FLAGB
	OPEN(UNIT=10,FILE="../DATA/FIND.DAT")
	OPEN(UNIT=20,FILE="../DATA/INPUT.DAT")
	sp=0 ! 初始值
	wide=0.1 ! 步长
	sp0=sp+wide
	DO i=1,13
		READ(20,*)x
		WRITE(10,"(2G12.6)")"x=",x
		sa=sp
		sb=sp
		sp=0
		n=0
		FLAGA=.TRUE. ! 标记，表示没有找到待求值范围
		FLAGB=.TRUE. ! 标记，表示没有找到待求值范围
		! DO WHILE(ABS(sa-sb)>CVG.OR.ABS(sa-sb)<1E-8)
		DO WHILE(ABS(x-sp)>CVG)
			n=n+1
			! ! 计算
			! CALL SLEEP(1)
			sp=0.5*(sa+sb)
			IF(ABS(x-sp)<CVG) THEN ! 如果遇到本次待求值与上次待求值相同则直接退出二分法循环
				EXIT
			ENDIF
			IF(sp<x) THEN !待求值在试探值右边
				FLAGA=.FALSE.
				sa=sp ! 正常二分法
				IF(FLAGA.OR.FLAGB) THEN ! 查看标记，看是否已经找到待求值范围，否则改变查找范围
					sb=sp+wide
					sp=sb
				ENDIF
			ELSE !待求值在试探值右边
				FLAGB=.FALSE.
				sb=sp ! 正常二分法
				IF(FLAGA.OR.FLAGB) THEN ! 查看标记，看是否已经找到待求值范围，否则改变查找范围
					sa=sp-wide
					sp=sa
				ENDIF
			ENDIF
		ENDDO
		wide=MAX(ABS(sp0-sp),0.02) ! 调整查找步长
		sp0=sp ! 作为下次的初始值
		IF(ABS(sp-x)>CVG) THEN
			WRITE(10,*)"WRONG!!!!!!!!!!!!!!!!!!!!!"
		ENDIF
		WRITE(*,"(4G12.6/G66.2)")"sp=",sp,"n=",n,"******************************************************************"
	ENDDO
END