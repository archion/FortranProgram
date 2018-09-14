module M_const
	use, intrinsic :: iso_fortran_env
	implicit none
	integer, parameter :: sp = REAL32
	integer, parameter :: dp = REAL64
	integer, parameter :: qp = REAL128
	integer, parameter :: wp=dp
	real(wp), parameter :: pi=atan(1._wp)*4._wp,zero=0._wp,one=1._wp,two=2._wp,half=0.5_wp,pi2=pi*pi
	complex(wp), parameter :: img=(0._wp,1._wp)
	real(dp), parameter :: nan=transfer([Z'00000000',Z'7FF80000'],1._dp)
	real(wp), parameter :: inf=huge(1._wp)+huge(1._wp)
	real(wp), parameter :: eps=epsilon(1._wp)
	integer, parameter :: reset=0,add_log_linear=1,add_linear=2
	logical :: underscore
end module
