include "mkl_dfti.f90"
include "../lib/fft.f90"
include "../lib/serde.f90"
program fft
	use M_serde
	use M_fft
	implicit none
	integer, parameter :: n=3
	!complex(8) :: x(0:n*n-1)
	real(8) :: x(0:n*n-1)
	complex(8) :: y(0:n*n-1)
	integer :: k=0,status
	type(dfti_descriptor), pointer :: handle
	type(t_serde(2)) :: map
	map%shap=[n,n]
	!do k=0,n**2-1
		!x(k:k)=merge(1d0,0d0,(map%get_idx(k+1,[1])-1)==1)
	!enddo
	x=[1d0,2d0,0.5d0,0.7d0,3d0,0d0,2d0,4d0,9d0]

	!call mfft(x,-1,[n,n],y,scal=1d0/size(x))
	!call mfft(x,1,[n,n])

	!call mfft(x,-1,[n,n],y,scal=1d0/size(x))
	!call mfft(y,1,[n,n],x)

	!call mfft(x,1,[n,n],y,scal=1d0/size(x))
	!call mfft(y,-1,[n,n],x)

	call mfft(x,1,[n,n],y,scal=1d0/size(x))
	call mfft(y,-1,[n,n],x)
	!write(*,*)merge(cmplx(0d0,kind=8),x,abs(x)<1d-10)
	!write(*,*)merge(cmplx(0d0,kind=8),y,abs(y)<1d-10)
	write(*,*)merge(0d0,x,abs(x)<1d-10)
end program fft



	
	
