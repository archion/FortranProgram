! 优化的二分法
! 对于任意的x，找出x
! by zuodong yu at 2012.12.18
! 要点：
! 1.从初始值出发以一定步长确定待求值的范围
! 2.确定待求值范围并进行二分法
! 3.得到的待求值进行下一次的初始值，并适当调整步长
! 4.如果遇到本次待求值与上次待求值相同则直接退出二分法循环
program bisection_method
	implicit none
	integer :: i,n
	real(8) :: sa,sb,sp,sp0,wide,cvg=1e-5,x
	logical :: flaga,flagb
	open(unit=10,file="../data/find.dat")
	open(unit=20,file="../data/input.dat")
	sp=0 ! 初始值
	wide=0.1 ! 步长
	sp0=sp+wide
	do i=1,13
		read(20,*)x
		write(10,"(2g12.6)")"x=",x
		sa=sp
		sb=sp
		sp=0
		n=0
		flaga=.true. ! 标记，表示没有找到待求值范围
		flagb=.true. ! 标记，表示没有找到待求值范围
		! do while(abs(sa-sb)>cvg.or.abs(sa-sb)<1e-8)
		do while(abs(x-sp)>cvg)
			n=n+1
			! ! 计算
			! call sleep(1)
			sp=0.5*(sa+sb)
			if(abs(x-sp)<cvg) then ! 如果遇到本次待求值与上次待求值相同则直接退出二分法循环
				exit
			endif
			if(sp<x) then !待求值在试探值右边
				flaga=.false.
				sa=sp ! 正常二分法
				if(flaga.or.flagb) then ! 查看标记，看是否已经找到待求值范围，否则改变查找范围
					sb=sp+wide
					sp=sb
				endif
			else !待求值在试探值右边
				flagb=.false.
				sb=sp ! 正常二分法
				if(flaga.or.flagb) then ! 查看标记，看是否已经找到待求值范围，否则改变查找范围
					sa=sp-wide
					sp=sa
				endif
			endif
		enddo
		wide=max(abs(sp0-sp),0.02) ! 调整查找步长
		sp0=sp ! 作为下次的初始值
		if(abs(sp-x)>cvg) then
			write(10,*)"wrong!!!!!!!!!!!!!!!!!!!!!"
		endif
		write(*,"(4g12.6/g66.2)")"sp=",sp,"n=",n,"******************************************************************"
	enddo
end
