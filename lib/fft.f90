include "mkl_dfti.f90"
module M_fft
	use M_const
	use mkl_dfti
	implicit none
	interface mfft
		module procedure fft_r2c,fft_c2r,fft_c2c
	end interface
contains
	subroutine fft_c2c(x,sig,n,y,scal)
		complex(8) :: x(:)
		integer :: sig,n(:)
		complex(8), optional :: y(:)
		real(8), optional :: scal
		integer :: status
		type(dfti_descriptor), pointer :: handle
		status=dfticreatedescriptor(handle,dfti_double,dfti_complex,size(n),n)
		if(present(scal)) then
			if(sig==1) then
				status=dftisetvalue(handle,dfti_backward_scale,scal)
			else
				status=dftisetvalue(handle,dfti_forward_scale,scal)
			endif
		endif
		if(present(y)) then
			status=dftisetvalue(handle,dfti_placement,dfti_not_inplace)
		endif
		status=dfticommitdescriptor(handle)
		if(sig==1) then
			if(present(y)) then
				status=dfticomputebackward(handle,x,y)
			else
				status=dfticomputebackward(handle,x)
			endif
		else
			if(present(y)) then
				status=dfticomputeforward(handle,x,y)
			else
				status=dfticomputeforward(handle,x)
			endif
		endif
		status=dftifreedescriptor(handle)
	end subroutine
	subroutine fft_r2c(x,sig,n,y,scal)
		real(8) :: x(:)
		integer :: sig,n(:)
		real(8), optional :: scal
		complex(8) :: y(:)
		integer :: status,i,NN,prod
		type(dfti_descriptor), pointer :: handle
		status=dfticreatedescriptor(handle,dfti_double,dfti_real,size(n),n)
		if(present(scal)) then
			status=dftisetvalue(handle,dfti_forward_scale,scal)
		endif
		status=dftisetvalue(handle,dfti_placement,dfti_not_inplace)
		status=dftisetvalue(handle,dfti_conjugate_even_storage,dfti_complex_complex)
		status=dfticommitdescriptor(handle)
		status=dfticomputeforward(handle,x,y)
		if(sig==1) then
			y=conjg(y)
		endif
		NN=2
		prod=1
		do i=1,size(n)
			prod=prod*n(i)
			NN=NN+prod
		enddo
		do i=n(1)/2+2,n(1)
			y(i)=conjg(y(n(1)-i+2))
			y(i+n(1)::n(1))=conjg(y(NN-i-n(1)::-n(1)))
		enddo
		status=dftifreedescriptor(handle)
	end subroutine
	subroutine fft_c2r(x,sig,n,y,scal)
		complex(8) :: x(:)
		integer :: sig,n(:)
		real(8), optional :: scal
		real(8) :: y(:)
		integer :: status
		type(dfti_descriptor), pointer :: handle
		status=dfticreatedescriptor(handle,dfti_double,dfti_real,size(n),n)
		if(present(scal)) then
			status=dftisetvalue(handle,dfti_backward_scale,scal)
		endif
		status=dftisetvalue(handle,dfti_placement,dfti_not_inplace)
		status=dftisetvalue(handle,dfti_conjugate_even_storage,dfti_complex_complex)
		status=dfticommitdescriptor(handle)
		if(sig==-1) then
			status=dfticomputebackward(handle,conjg(x),y)
		else
			status=dfticomputebackward(handle,x,y)
		endif
		status=dftifreedescriptor(handle)
	end subroutine
end module
