include "pade.f90"
include "pade2.f90"
program main
	! This is a test program
	use fft
	use green
	implicit none
	real :: start,finish
	!integer,parameter :: n=1048576
	!integer :: i,j
	!complex(8) :: a(0:n-1,0:n-1),c(0:n-1)
	!real(8) :: b(0:n-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                     chack the bit-reversal                !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! do i=0,n-1
		!! a(i)=i
	!! enddo
	!c=0d0
	!!write(10,"('#data')")
	!do i=1,20
		!c(1)=1d0
		!!call fft1d_2(c,1)
		!call fft1d_2(c,1)
	!enddo
	!!write(10,"(2(e11.3))")c
	!!write(*,"(i6,b19.16,i6,b19.16)")(int(a(i)),int(a(i)),i,i,i=0,n-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                     chack the fft                         !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! a=(/0,25,2,3,0,3,2,1/)
	!! call fftreal1_rp(a,n,1)
	!! call fftreal1_rp(a,n,-1)
	!! ! call fft1(a,n)
	!! write(10,"(e11.3)")a
	!!b=(/1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
	!!call fftcos1_2(b,1)
	!!write(10,"(e11.3)")b
	!!call fftcos1_2(b,-1)
	!! call fftcos1_2(b,n,-1)
	!! write(10,"(e11.3)")b
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!      obtain the whole fft array in replaces fft subroutine        !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! write(10,"(2e11.3)")b(0),0d0
	!! write(10,"(2e11.3)")b(2:n-1)
	!! write(10,"(2e11.3)")b(1),0d0
	!! b(2:n-2:2)=b(n-2:2:-2)
	!! b(3:n-1:2)=-1*b(n-1:3:-2)
	!! write(10,"(2e11.3)")b(2:n-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                 chack pade approximate                    !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer, parameter :: n=512,m=200
	complex(8) :: lh(n),lhf(m),lhc(m),z(n),omg(m)
	complex(8) :: img=(0d0,1d0)
	real(8) :: Tk,k(2),fk,fkq,ek,ekq,omg1(m)
	integer :: i,j,l,mk=10
	open(unit=10,file="../data/test.dat")
	write(10,"(A)")"#data"
	lh=0d0
	lhc=0d0
	Tk=0.001d0
	omg1=(/(real(i,8)/80d0,i=1,m)/)
	omg=(/(dcmplx(i/80d0,0.001d0),i=1,m)/)
	z=(/(cmplx(0d0,pi*Tk*(2d0*(i-1)+1d0)),i=1,n)/)
	do j=1,mk
		do l=1,mk
			k=(/pi/mk*j,pi/mk*l/)
			ek=(cos(k(1))+cos(k(2)))
			ekq=(cos(k(1)+pi)+cos(k(2)+pi))
			fk=1d0/(1d0+exp(ek/Tk))
			fkq=1d0/(1d0+exp(ekq/Tk))
			do i=1,n
				!lh(i)=lh(i)+(fk-fkq)/(cmplx(0d0,pi*Tk*(2d0*(i-1)+1d0))+ek-ekq)
				lh(i)=lh(i)+1d0/(cmplx(0d0,pi*Tk*(2d0*(i-1)+1d0))-ek)
			enddo
			do i=1,m
				!lhc(i)=lhc(i)+(fk-fkq)/(omg(i)+ek-ekq+img*0.01d0)
				!lhc(i)=lhc(i)+1d0/(omg1(i)-ek+img*0.001d0)
				lhc(i)=lhc(i)+1d0/(omg(i)-ek)
			enddo
		enddo
	enddo
	lh=lh/(mk*mk)
	lhc=lhc/(mk*mk)
	call cpu_time(start)
	do i=1,2000
		call pade(lh,z,omg,0.001d0,lhf)
		!call pade2(n,m,lh,z,omg,0.001d0,lhf)
		!call pade1(size(lh),z,lh,size(omg),omg,lhf)
	enddo
	call cpu_time(finish)
	write(10,"(3e16.4)")(real(omg(i)),dimag(lhf(i)),dimag(lhc(i)),i=1,m)
	write(*,*)finish-start
end
