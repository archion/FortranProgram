! �Ż��Ķ��ַ�
! ���������x���ҳ�x
! by zuodong yu at 2012.12.18
! Ҫ�㣺
! 1.�ӳ�ʼֵ������һ������ȷ������ֵ�ķ�Χ
! 2.ȷ������ֵ��Χ�����ж��ַ�
! 3.�õ��Ĵ���ֵ������һ�εĳ�ʼֵ�����ʵ���������
! 4.����������δ���ֵ���ϴδ���ֵ��ͬ��ֱ���˳����ַ�ѭ��
program bisection_method
	implicit none
	integer :: i,n
	real(8) :: sa,sb,sp,sp0,wide,cvg=1e-5,x
	logical :: flaga,flagb
	open(unit=10,file="../data/find.dat")
	open(unit=20,file="../data/input.dat")
	sp=0 ! ��ʼֵ
	wide=0.1 ! ����
	sp0=sp+wide
	do i=1,13
		read(20,*)x
		write(10,"(2g12.6)")"x=",x
		sa=sp
		sb=sp
		sp=0
		n=0
		flaga=.true. ! ��ǣ���ʾû���ҵ�����ֵ��Χ
		flagb=.true. ! ��ǣ���ʾû���ҵ�����ֵ��Χ
		! do while(abs(sa-sb)>cvg.or.abs(sa-sb)<1e-8)
		do while(abs(x-sp)>cvg)
			n=n+1
			! ! ����
			! call sleep(1)
			sp=0.5*(sa+sb)
			if(abs(x-sp)<cvg) then ! ����������δ���ֵ���ϴδ���ֵ��ͬ��ֱ���˳����ַ�ѭ��
				exit
			endif
			if(sp<x) then !����ֵ����ֵ̽�ұ�
				flaga=.false.
				sa=sp ! �������ַ�
				if(flaga.or.flagb) then ! �鿴��ǣ����Ƿ��Ѿ��ҵ�����ֵ��Χ������ı���ҷ�Χ
					sb=sp+wide
					sp=sb
				endif
			else !����ֵ����ֵ̽�ұ�
				flagb=.false.
				sb=sp ! �������ַ�
				if(flaga.or.flagb) then ! �鿴��ǣ����Ƿ��Ѿ��ҵ�����ֵ��Χ������ı���ҷ�Χ
					sa=sp-wide
					sp=sa
				endif
			endif
		enddo
		wide=max(abs(sp0-sp),0.02) ! �������Ҳ���
		sp0=sp ! ��Ϊ�´εĳ�ʼֵ
		if(abs(sp-x)>cvg) then
			write(10,*)"wrong!!!!!!!!!!!!!!!!!!!!!"
		endif
		write(*,"(4g12.6/g66.2)")"sp=",sp,"n=",n,"******************************************************************"
	enddo
end
