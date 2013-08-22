program main
	implicit none
	real(8),parameter :: pi=3.1415926
	integer,parameter :: n=8
	integer :: i,j
	real(8) :: a(0:n-1)
	real(8) :: b(0:n-1)
	open(unit=10,file="../data/output.dat")
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                     chack the bit-reversal                !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! do i=0,n-1
		! a(i)=i
	! enddo
	! call fft(a,1,n)
	! write(*,"(i6,b19.16,i6,b19.16)")(int(a(i)),int(a(i)),i,i,i=0,n-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                     chack the fft                         !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! a=(/0,25,2,3,0,3,2,1/)
	! call fftreal1_rp(a,n,1)
	! call fftreal1_rp(a,n,-1)
	! ! call fft1(a,n)
	! write(10,"(e11.3)")a
	b=(/1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
	call fftcos1_2(b,n,1)
	write(10,"(e11.3)")b
	call fftcos1_2(b,n,-1)
	! call fftcos1_2(b,n,-1)
	! write(10,"(e11.3)")b
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!      obtain the whole fft array in replaces fft subroutine        !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! write(10,"(2e11.3)")b(0),0d0
	! write(10,"(2e11.3)")b(2:n-1)
	! write(10,"(2e11.3)")b(1),0d0
	! b(2:n-2:2)=b(n-2:2:-2)
	! b(3:n-1:2)=-1*b(n-1:3:-2)
	! write(10,"(2e11.3)")b(2:n-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                 chack the fft of sine                     !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!       fft for complex array in one-demension              !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft1(a,n,sig)
	implicit none
	real(8),parameter :: pi=3.1415926
	integer :: i,j,m,sig,n
	complex(8) :: a(0:n-1),tmp,wp,w
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!          bit-reversal using gold-rader algorithm          !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	j=n/2
	do i=1,n-2
		if(j>i) then
			tmp=a(i)
			a(i)=a(j)
			a(j)=tmp
		endif
		m=n/2
		do while(j>=m)
			j=j-m
			m=m/2
		enddo
		j=j+m
	enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                     fast fourier transform using danielson-lanczos lemma                          !!!!
!!!!     trigonometric recurrence relations using 5.5.6 in book "numerical recipes in fortran 77"      !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	m=2
	do while(m<=n)
		wp=cmplx(-2d0*sin(sig*pi/m)**2,sin(sig*2d0*pi/m))
		w=cmplx(1d0,0d0)
		do j=0,m/2-1
			do i=j,n-1,m
				tmp=w*a(i+m/2)
				a(i+m/2)=a(i)-tmp
				a(i)=a(i)+tmp
			enddo
			w=w*wp+w
		enddo
		m=m*2
	enddo
	if(sig==-1) then
		a=a/n
	endif
end




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!          fft for real array in one-demension              !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fftreal1(a,b,n)
	implicit none
	real(8),parameter :: pi=3.1415926
	integer :: i,m,n
	real(8) :: a(0:n-1)
	complex(8) :: b(0:n-1),h1(0:n/2-1),h2(0:n/2-1),wp,w
	b=cmplx(a(0:n-2:2),a(1:n-1:2))
	m=n/2
	call fft1(b(0:m-1),m,1)
	h1(0)=cmplx(real(b(0)),0)
	h2(0)=cmplx(aimag(b(0)),0)
	h1(1:m/2)=0.5*(b(1:m/2)+conjg(b(m-1:m/2:-1)))
	h2(1:m/2)=cmplx(0,-0.5)*(b(1:m/2)-conjg(b(m-1:m/2:-1)))
	h1(m-1:m/2:-1)=conjg(h1(1:m/2))
	h2(m-1:m/2:-1)=conjg(h2(1:m/2))
	wp=cmplx(-2d0*sin(pi/n)**2,sin(2d0*pi/n))
	w=cmplx(1d0,0d0)
	do i=0,m-1
		b(i)=h1(i)+w*h2(i)
		b(i+m)=h1(i)-w*h2(i)
		w=w*wp+w
	enddo
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                      fft for real array in one-demension                            !!!!
!!!!               using input real array to storage the output complex array            !!!!
!!!!     with even index refer to real parts and odd index refer to imaginary parts      !!!!
!!!!           and a(0) is 0th while a(1) is n/2th of the fft of input data              !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fftreal1_rp(a,n,sig)
	implicit none
	real(8),parameter :: pi=3.1415926d0
	integer :: i,m,n,sig
	real(8) :: a(0:n-1)
	complex(8) :: b(0:n/2-1),wp,w,he,ho,tmp
	m=n/2
	if(sig==1) then
		b=cmplx(a(0:n-2:2),a(1:n-1:2))
		call fft1(b,m,1)
		wp=cmplx(-2d0*sin(pi/n)**2,sin(2d0*pi/n))
		w=cmplx(1d0,0d0)
		w=w*wp+w
		a(0)=real(b(0))+aimag(b(0))
		a(1)=real(b(0))-aimag(b(0))
		do i=1,m-1
			he=0.5*(b(i)+conjg(b(m-i)))
			ho=cmplx(0,-0.5)*(b(i)-conjg(b(m-i)))
			tmp=he+w*ho
			a(2*i)=real(tmp)
			a(2*i+1)=aimag(tmp)
			w=w*wp+w
		enddo
	else
		b(1:m-1)=cmplx(a(2:n-2:2),a(3:n-1:2))
		wp=cmplx(-2d0*sin(-1d0*pi/n)**2,sin(-2d0*pi/n))
		w=cmplx(1d0,0d0)
		w=w*wp+w
		do i=1,m/2
			he=0.5*(b(i)+conjg(b(m-i)))
			ho=cmplx(0d0,0.5d0)*(b(i)-conjg(b(m-i)))
			b(i)=he+w*ho
			b(m-i)=conjg(he)-conjg(w*ho)
			w=w*wp+w
		enddo
		b(0)=0.5*cmplx(a(0)+a(1),a(0)-a(1))
		call fft1(b,m,-1)
		do i=0,m-1
			a(2*i)=real(b(i))
			a(2*i+1)=aimag(b(i))
		enddo
	endif
end

subroutine fftsin1(a,n)
	implicit none
	real(8),parameter :: pi=3.1415926d0
	integer :: i,n
	real(8) :: a(0:n-1),tmp
	complex(8) :: wp,w
	wp=cmplx(-2d0*sin(0.5d0*pi/n)**2,sin(pi/n))
	w=cmplx(1,0)
	a(0)=0
	do i=1,n/2
		w=w*wp+w
		tmp=0.5*(a(i)-a(n-i))
		a(i)=aimag(w)*(a(i)+a(n-i))+tmp
		a(n-i)=a(i)-2*tmp
	enddo
	call fftreal1_rp(a,n,1)
	a(1)=0.5*a(0)
	a(0)=0
	do i=3,n-1,2
		tmp=a(i)
		a(i)=a(i-2)+a(i-1)
		a(i-1)=tmp
	enddo
end
subroutine fftcos1(a,n)
	implicit none
	real(8),parameter :: pi=3.1415926d0
	integer :: i,n
	real(8) :: a(0:n),sun,tmp
	complex(8) :: wp,w
	wp=cmplx(-2d0*sin(0.5d0*pi/n)**2,sin(pi/n))
	w=cmplx(1,0)
	sun=0.5*(a(0)-a(n))
	a(0)=0.5*(a(0)+a(n))
	do i=1,n/2
		w=w*wp+w
		tmp=a(n-i)-a(i)
		a(i)=0.5*(a(i)+a(n-i))+aimag(w)*tmp
		a(n-i)=a(i)-2*aimag(w)*tmp
		sun=sun-tmp*real(w)
	enddo
	call fftreal1_rp(a(0:n-1),n,1)
	a(n)=a(1)
	a(1)=sun
	do i=3,n-1,2
		a(i)=a(i-2)+a(i)
	enddo
end
subroutine fftcos1_2(a,n,sig)
	implicit none
	real(8),parameter :: pi=3.1415926d0
	integer :: i,n,m,sig
	real(8) :: a(0:n-1),tmp,tmp0,tmp1
	complex(8) :: wp,w,wi,wpi
	m=n/2
	w=cmplx(cos(0.5d0*pi/n),sin(0.5d0*pi/n))
	wi=cmplx(cos(0.5d0*pi),sin(0.5d0*pi))
	wp=cmplx(-2d0*sin(0.5d0*pi/n)**2,sin(pi/n))
	wpi=cmplx(-2d0*sin(-0.5d0*pi/n)**2,sin(-1d0*pi/n))
	if(sig==1) then
		do i=0,m-1
			tmp=aimag(w)*(a(n-i-1)-a(i))
			a(i)=0.5*(a(i)+a(n-i-1))+tmp
			a(n-i-1)=a(i)-2*tmp
			w=w*wp+w
		enddo
		call fftreal1_rp(a,n,sig)
		tmp=a(n-1)
		a(n-1)=0.5*a(1)
		do i=n-1,3,-2
			wi=wi*wpi+wi
			tmp0=a(i-1)
			tmp1=a(i-2)
			a(i-1)=real(wi)*a(i-1)-aimag(wi)*tmp
			a(i-2)=aimag(wi)*tmp0+real(wi)*tmp+a(i)
			tmp=tmp1
		enddo
	else
		tmp0=a(n-1)
		do i=n-1,3,-2
			wi=wi*wpi+wi
			a(i)=a(i-2)-a(i)
			tmp=a(i)*real(wi)-a(i-1)*aimag(wi)
			a(i-1)=a(i-1)*real(wi)+a(i)*aimag(wi)
			a(i)=tmp
		enddo
		a(1)=2*tmp0
		call fftreal1_rp(a,n,sig)
		do i=0,m-1
			tmp=a(i)+a(n-i-1)
			a(i)=(a(n-i-1)-a(i))/(2*aimag(w))
			a(n-i-1)=0.5*(tmp-a(i))
			a(i)=0.5*(tmp+a(i))
			w=w*wp+w
		enddo
	endif
end















