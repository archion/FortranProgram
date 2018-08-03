module M_fn
	use M_const
	real(8), parameter :: zero = 0.0d0, one = 1.0d0,    &
		two = 2.0d0, half = 0.5d0, rmin = tiny( one ),              &
		eps0 = epsilon( one ), sqrt_log_rmin = sqrt( -log( rmin ) ),        &
		pi2 = pi * pi, one_sqrt_pi = one / sqrt( pi )
	complex(8), parameter ::                                    &
		cmplxj = cmplx( zero, one, kind=8 ),                            &
		cmplx0 = cmplx( zero, zero, kind=8 )
contains
	complex(8) function lgamma(z) result(rt)
		complex(8) :: z
		real(8) :: reim(2),err
		call zgam([z.re,z.im],reim,err,0)
		rt=cmplx(reim(1),reim(2))
	end function
	complex(8) function zgamma(z) result(rt)
		complex(8) :: z
		real(8) :: reim(2),err
		call zgam([z.re,z.im],reim,err,0)
		rt=exp(cmplx(reim(1),reim(2)))
	end function
	!*********************************************************************72
	elemental complex(8) function wzag( z, err )

		! This is the Faddeyeva or the plasma dispersion function,
		! w(z) = exp(-z**2) * erfc(-i*z). erfc(z) is the complex complementary
		! error function of z. z is a complex number.
		!
		! Adapted from the Matlab code implementing TOMS
		! Algorithm 916: http://www.netlib.org/toms/
		!
		! file: 916.zip
		! ref: TOMS 38,2 (Dec 2011) Article: 15
		! for: Computing the Faddeyeva and Voigt Functions
		! by: Mofreh R. Zaghloul and Ahmed N. Ali
		!
		! Most of the code is calculation of equations (13)-(19) of the
		! above paper.
		!
		! Inputs:
		!    z - the argument of the function 
		!  err - The desired accuracy (positive). err must be (err .le. 1.0e-4)
		!        and  (err .ge. 0.06447*epsilon). For efficiency,
		!        no checks are made!
		!        The lowest accuracy and the fastest calculation are obtained
		!        with err .eq. 1.0e-4. For higher accuracy use smaller err.

		complex(8), intent(in) :: z
		real(8), intent(in) :: err 

		real(8) :: a, half_a, a_sqr, two_a, two_a_sqr, &
			four_a_sqr, myerr, a_pi, two_a_pi, x, y, erfcsy, &
			xsign, ysign, x_sqr, y_sqr, two_yx, two_a_x, exp_x_sqr, &
			cos_2yx, sin_2yx, v_old, l_old, sigma1, sigma2, sigma3, &
			sigma4, sigma5, sigma4_5, delta3, delta5, exp1, exp2, exp3, &
			del2_tmp, del3_tmp, del3_3_tmp, den1, exp_del1, &
			exp3_den, exp3_3_den, del5, del3, two_exp_x_sqr_ysqr, aux13
		integer :: n, n3, n3_3

		x = real( z )
		y = aimag( z )

		! For purely imaginary z, use intrinsic scaled complement of
		! the error function, erfc_scaled (F2008 and beyond).
		! Return immediately.
		if ( abs( x ) .eq. zero ) then
			wzag = erfc_scaled( y )
			return
		end if

		myerr = max( err, eps0 )
		a = sqrt( - pi2 / log( err / 2d0 ) )
		half_a = half * a
		a_sqr = a**2
		two_a = 2 * a
		two_a_sqr = 2 * a_sqr
		four_a_sqr = 4 * a_sqr
		a_pi = a / pi
		two_a_pi = 2 * a_pi
		erfcsy = erfc_scaled( abs( y ) )
		xsign = sign( one, x )
		ysign = sign( one, y )
		x = abs( x )
		y = max( rmin, abs( y ) )
		x_sqr = x**2
		y_sqr = y**2
		two_yx = 2 * y * x
		two_a_x = two_a * x
		exp_x_sqr = exp( -x_sqr )
		cos_2yx = cos( two_yx )
		sin_2yx = sin( two_yx )
		v_old = exp_x_sqr * &
			( erfcsy * cos_2yx + two_a_pi * sin( two_yx / 2)**2 / y )
		l_old = - erfcsy + a_pi / y
		sigma3 = rmin
		sigma5 = rmin
		sigma4_5 = zero
		delta3 = one
		delta5 = one
		n = 0
		n3 = ceiling( x / a )
		n3_3 = n3 - 1

		outer: if ( (sqrt_log_rmin - x) .gt. 0 ) then
			sigma1 = rmin
			sigma2 = rmin
			sigma4 = rmin
			exp1 = exp( - two_a_x )
			exp2 = exp( four_a_sqr * n3 - 2 * two_a_x - two_a_sqr )
			exp3 = exp( - ( two_a_sqr * n3 - two_a_x - two_a_sqr ) )
			del2_tmp = one
			del3_tmp = exp( - ( a_sqr * n3**2 - two_a_x * n3                  &
				- two_a_sqr * n3 + x_sqr + two_a_x + a_sqr ) )
			del3_3_tmp = exp( a_sqr - ( two_a_sqr * n3 - two_a_x ) )

			loop1: do
				if ( delta3 .lt. myerr .and. delta5 .lt. myerr .and.              &
					n .gt. 50 ) exit loop1
				n = n + 1
				den1 = a_sqr * n**2 + y_sqr
				exp_del1 = exp( -( a_sqr * n**2) ) / den1
				del2_tmp = del2_tmp * exp1

				minor: if ( n3_3 .ge. 1 ) then
					del3_tmp = del3_tmp * exp3
					exp3_den = del3_tmp * exp_del1 *                              &
						( den1 / ( a_sqr * n3**2 + y_sqr) )
					del3_3_tmp = del3_3_tmp * exp2
					exp3_3_den = exp3_den * del3_3_tmp *                            &
						( (a_sqr * n3**2 + y_sqr ) /                      &
						(a_sqr * n3_3**2 + y_sqr ) )
					del5 = n3_3 * exp3_3_den + n3 * exp3_den
					del3 = exp3_3_den + exp3_den
				else
					del3_tmp = del3_tmp * exp3
					del3 = del3_tmp * exp_del1 *                              &
						( den1 / ( a_sqr * n3**2 + y_sqr ) )
					del5 = n3 * del3
				end if minor

				delta3 = del3 / sigma3
				delta5 = del5 / sigma5
				sigma1 = sigma1 + exp_del1
				sigma2 = sigma2 + del2_tmp * exp_x_sqr * exp_del1
				sigma3 = sigma3 + del3
				sigma4 = sigma4 + n * del2_tmp * exp_x_sqr * exp_del1
				sigma5 = sigma5 + del5

				if ( x .ge. 5d-4 ) then
					sigma4_5 = - sigma4 + sigma5
				else
					sigma4_5 = sigma4_5 + 2 * n**2 * two_a_x * exp_x_sqr            &
						* exp_del1 * ( one + 1.666666666666667d-1     &
						* ( two_a_x * n )**2 + 8.333333333333333d-3   &
						* ( two_a_x * n )**4 )
				end if

				n3 = n3 + 1
				n3_3 = n3_3 - 1

			end do loop1

			! Second line of Eqn (13)
			aux13 = y * two_a_pi *                                              &
				( - cos_2yx*exp_x_sqr*sigma1 + half*(sigma2 + sigma3) )
			mumu: if ( y .le. 5d0 .and. two_yx .gt. rmin ) then
				wzag = v_old + aux13 + cmplxj * xsign             &
					* ( sin_2yx * exp_x_sqr * ( l_old + two_a_pi * y * sigma1 ) &
					+ two_a_pi * half_a * sigma4_5 )
			else if ( y .le. 5d0 .and. two_yx .le. rmin ) then
				wzag = v_old + aux13 + cmplxj * xsign             &
					* ( two * y * exp_x_sqr * ( x * l_old + x * two_a_pi * y    &
					* sigma1 ) + two_a_pi * half_a * sigma4_5 )
			else
				wzag = v_old + aux13 + cmplxj * xsign             &
					* ( sin_2yx * exp_x_sqr * min( zero, abs( l_old             &
					+ ( two_a_pi * y * sigma1) ) ) + two_a_pi *half_a           &
					* sigma4_5 )
			end if mumu

		else if ( x .ge. sqrt_log_rmin .and. x .lt. 1.0d15 ) then

			exp2 = exp( four_a_sqr * n3 - 2 * two_a_x - two_a_sqr )
			del3_3_tmp = exp( a_sqr + two_a_x - two_a_sqr * n3 )

			loop2: do
				if ( delta3 .lt. myerr .and. delta5 .lt. myerr .and.              &
					n .gt. 50 ) exit loop2
				n = n + 1
				if ( n3_3 .ge. 1 ) then
					exp3_den = exp( - ( a * n3 - x ) * ( a * n3 - x ) )           &
						/ ( a_sqr * n3**2 + y_sqr )
					del3_3_tmp = del3_3_tmp * exp2
					exp3_3_den = exp3_den * del3_3_tmp * ( ( a_sqr * n3**2 + y_sqr) &
						/ ( a_sqr * n3_3**2 + y_sqr ) )
					del5 = n3_3 * exp3_3_den + n3 * exp3_den
					del3 = exp3_3_den + exp3_den
				else
					del3 = exp( - (a * n3 - x)**2) / ( a_sqr * n3**2 + y_sqr )
					del5 = n3 * del3
				end if

				delta3 = del3 / sigma3
				delta5 = del5 / sigma5
				sigma3 = sigma3 + del3
				sigma5 = sigma5 + del5
				n3 = n3 + 1
				n3_3 = n3_3 - 1

			end do loop2

			wzag = v_old + y * a_pi * sigma3 + cmplxj * xsign * ( sin_2yx        &
				* exp_x_sqr * l_old + two_a_pi * half_a * sigma5 )

		else
			wzag = one_sqrt_pi * ( ( y + cmplxj * xsign * x ) / ( x_sqr + y_sqr ))
		end if outer

		if ( ysign .lt. zero ) then
			two_exp_x_sqr_ysqr = two * exp( - x_sqr + y_sqr )
			wzag = two_exp_x_sqr_ysqr * cos_2yx - real( wzag ) - cmplxj           &
				* ( - xsign * two_exp_x_sqr_ysqr * sin_2yx - aimag( wzag ) ) 
		end if

	end function
end module
