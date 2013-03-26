PROGRAM RANDOM
	IMPLICIT NONE
	INTEGER :: m=2147483647,q=127773,r=2836,a=16807,z=245698,x,y,i
	OPEN(UNIT=10,FILE='D:\FortranProgram\DATA\RANDOM.DAT')
	DO i=1,30000000
		z=a*MOD(z,q)-r*INT(z/q)
		IF(z>=0) THEN
			x=z
		ELSE
			x=z+m
		ENDIF
		z=a*MOD(z,q)-r*INT(z/q)
		IF(z>=0) THEN
			y=z
		ELSE
			y=z+m
		ENDIF
		WRITE(*,*)i
		WRITE(10,*)x,y
	ENDDO	
END
	
	
	