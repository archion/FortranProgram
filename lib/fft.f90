module M_fft
	use M_const
	implicit none
contains
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!       fft for complex array in one-demension              !!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine fft1d(a,sig)
		integer :: i,j,m,sig,n
		complex(8) :: a(0:),tmp,wp,w
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!          bit-reversal using gold-rader algorithm          !!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		n=size(a)
		if(iand(n,n-1)/=0) then
			write(*,*)"the size must be power of 2"
			stop
		endif
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
	end subroutine

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!        one-demension fft for two-demension complex data   !!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine fftrow2(a,sig)
		integer :: i,j,m,sig,n
		complex(8) :: a(0:,0:),tmp(size(a,1)),wp,w
		n=size(a,2)
		if(iand(n,n-1)/=0) then
			write(*,*)"the size must be power of 2"
			stop
		endif
		j=n/2
		do i=1,n-2
			if(j>i) then
				tmp=a(:,i)
				a(:,i)=a(:,j)
				a(:,j)=tmp
			endif
			m=n/2
			do while(j>=m)
				j=j-m
				m=m/2
			enddo
			j=j+m
		enddo
		m=2
		do while(m<=n)
			wp=cmplx(-2d0*sin(sig*pi/m)**2,sin(sig*2d0*pi/m))
			w=cmplx(1d0,0d0)
			do j=0,m/2-1
				do i=j,n-1,m
					tmp=w*a(:,i+m/2)
					a(:,i+m/2)=a(:,i)-tmp
					a(:,i)=a(:,i)+tmp
				enddo
				w=w*wp+w
			enddo
			m=m*2
		enddo
		if(sig==-1) then
			a=a/n
		endif
	end subroutine
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!        one-demension fft for three-demension complex data   !!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine fftrow3(a,sig)
		integer :: i,j,m,sig,n
		complex(8) :: a(0:,0:,0:),tmp(size(a,1),size(a,2)),wp,w
		n=size(a,3)
		if(iand(n,n-1)/=0) then
			write(*,*)"the size must be power of 2"
			stop
		endif
		j=n/2
		do i=1,n-2
			if(j>i) then
				tmp=a(:,:,i)
				a(:,:,i)=a(:,:,j)
				a(:,:,j)=tmp
			endif
			m=n/2
			do while(j>=m)
				j=j-m
				m=m/2
			enddo
			j=j+m
		enddo
		m=2
		do while(m<=n)
			wp=cmplx(-2d0*sin(sig*pi/m)**2,sin(sig*2d0*pi/m))
			w=cmplx(1d0,0d0)
			do j=0,m/2-1
				do i=j,n-1,m
					tmp=w*a(:,:,i+m/2)
					a(:,:,i+m/2)=a(:,:,i)-tmp
					a(:,:,i)=a(:,:,i)+tmp
				enddo
				w=w*wp+w
			enddo
			m=m*2
		enddo
		if(sig==-1) then
			a=a/n
		endif
	end subroutine
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!          fft for real array in one-demension              !!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine fftreal1d(a,b)
		integer :: i,m,n
		real(8) :: a(0:)
		complex(8) :: b(0:),h1(0:size(a)/2-1),h2(0:size(a)/2-1),wp,w
		n=size(a)
		b=cmplx(a(0:n-2:2),a(1:n-1:2))
		m=n/2
		call fft1d(b(0:m-1),1)
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
	end subroutine
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!                      fft for real array in one-demension                            !!!!
	!!!!               using input real array to storage the output complex array            !!!!
	!!!!     with even index refer to real parts and odd index refer to imaginary parts      !!!!
	!!!!           and a(0) is 0th while a(1) is n/2th of the fft of input data              !!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine fftreal1d_rp(a,sig)
		integer :: i,m,n,sig
		real(8) :: a(0:)
		complex(8) :: b(0:size(a)/2-1),wp,w,he,ho,tmp
		n=size(a)
		m=n/2
		if(sig==1) then
			b=cmplx(a(0:n-2:2),a(1:n-1:2))
			call fft1d(b,1)
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
			call fft1d(b,-1)
			do i=0,m-1
				a(2*i)=real(b(i))
				a(2*i+1)=aimag(b(i))
			enddo
		endif
	end subroutine

	subroutine fftsin1d(a)
		integer :: i,n
		real(8) :: a(0:),tmp
		complex(8) :: wp,w
		n=size(a)
		wp=cmplx(-2d0*sin(0.5d0*pi/n)**2,sin(pi/n))
		w=cmplx(1,0)
		a(0)=0
		do i=1,n/2
			w=w*wp+w
			tmp=0.5*(a(i)-a(n-i))
			a(i)=aimag(w)*(a(i)+a(n-i))+tmp
			a(n-i)=a(i)-2*tmp
		enddo
		call fftreal1d_rp(a,1)
		a(1)=0.5*a(0)
		a(0)=0
		do i=3,n-1,2
			tmp=a(i)
			a(i)=a(i-2)+a(i-1)
			a(i-1)=tmp
		enddo
	end subroutine
	subroutine fftcos1d(a)
		integer :: i,n
		real(8) :: a(0:),sun,tmp
		complex(8) :: wp,w
		n=size(a)-1
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
		call fftreal1d_rp(a(0:n-1),1)
		a(n)=a(1)
		a(1)=sun
		do i=3,n-1,2
			a(i)=a(i-2)+a(i)
		enddo
	end subroutine
	subroutine fftcos1d_2(a,sig)
		integer :: i,n,m,sig
		real(8) :: a(0:),tmp,tmp0,tmp1
		complex(8) :: wp,w,wi,wpi
		n=size(a)
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
			call fftreal1d_rp(a,sig)
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
			call fftreal1d_rp(a,sig)
			do i=0,m-1
				tmp=a(i)+a(n-i-1)
				a(i)=(a(n-i-1)-a(i))/(2*aimag(w))
				a(n-i-1)=0.5*(tmp-a(i))
				a(i)=0.5*(tmp+a(i))
				w=w*wp+w
			enddo
		endif
	end subroutine
	subroutine fft2d(a,sig)
		complex(8) :: a(0:,0:),tmp(size(a,2),size(a,1))
		integer :: sig,n1,n2
		n1=size(a,1)
		n2=size(a,2)
		call fftrow2(a,sig)
		tmp=transpose(a)
		call fftrow2(tmp,sig)
		a=transpose(tmp)
	end subroutine
	subroutine fft3d(a,sig)
		complex(8) :: a(0:,0:,0:),tmp1(size(a,2),size(a,3),size(a,1)),tmp2(size(a,3),size(a,1),size(a,2))
		integer :: sig,n1,n2,n3
		n1=size(a,1)
		n2=size(a,2)
		n3=size(a,3)
		call fftrow3(a,sig)
		tmp1=reshape(a,shape=(/n2,n3,n1/),order=(/3,1,2/))
		call fftrow3(tmp1,sig)
		tmp2=reshape(tmp1,shape=(/n3,n1,n2/),order=(/3,1,2/))
		call fftrow3(tmp2,sig)
		a=reshape(tmp2,shape=(/n1,n2,n3/),order=(/3,1,2/))
	end subroutine
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!        one-demension fft that transform the input         !!!!
	!!!!       to two-demension fft, this subroutine should much   !!!!
	!!!!                    faster than fft1d                      !!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine fft1d_2(a,sig)
		complex(8) :: a(0:)
		complex(8), allocatable	:: tmp1(:,:),tmp2(:,:),w(:),wp(:)
		real(8), allocatable	:: m(:)
		integer :: sig,i,n,n1,n2
		n=size(a)
		if(iand(n,n-1)/=0) then
			write(*,*)"the size must be power of 2"
			stop
		endif
		n1=2**ceiling(0.6d0*log(real(n,8))/0.693147d0)
		n2=n/n1
		allocate(tmp1(n1,n2),tmp2(n2,n1),m(n1),w(n1),wp(n1))
		tmp1=reshape(a,(/n1,n2/))
		call fftrow2(tmp1,sig)
		m=(/(real(i,8),i=0,n1-1)/)
		wp=cmplx(-2d0*sin(sig*pi*m/n)**2,sin(sig*2d0*pi*m/n))
		w=cmplx(1d0,0d0)
		do i=2,n2
			w=w*wp+w
			tmp1(:,i)=tmp1(:,i)*w
		enddo
		tmp2=transpose(tmp1)
		call fftrow2(tmp2,sig)
		a=reshape(tmp2,(/n/))
		deallocate(tmp1,tmp2,m,w,wp)
	end subroutine
end module
