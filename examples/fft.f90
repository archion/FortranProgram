program fft
	use mkl_dfti
	implicit none
	real, parameter :: pi=3.1415927
	complex, parameter :: i=(0.0,1.0)
	integer, parameter :: n=2048
	complex :: x(0:n-1)=0
	integer :: k=0,status
	real :: y(0:n-1)=0,xx(0:n-1)=0
	type(dfti_descriptor), pointer :: handle
	open(unit=10,file='../data/output.dat')
	x(0)=1
	do k=2,n/2-1
		x(k)=-x(k-2)
	enddo
	do k=1,n/2-1
		x(k)=2*x(k)*exp(i*pi*k/(2*n))*sinh(4.0*(1.0-real(k)/real(n)*2.0))/sinh(4.0)
	enddo
	status=dfticreatedescriptor(handle,dfti_single,dfti_complex,1,n)
	status=dfticommitdescriptor(handle)
	status=dfticomputebackward(handle,x)
	status=dftifreedescriptor(handle)
	do k=0,n/2-1
		y(2*k)=real(x(k))
		y(2*k+1)=real(x(n-1-k))
	enddo
	do k=0,n-1
		xx(k)=cos(pi*(k+0.5)/n)
	enddo
	y=y/(pi*sqrt(1-xx*xx))
	write(10,"(2e15.7)")(xx(k),y(k),k=0,n-1)
	close(10)
end program fft



	
	
