module stdlib_stats_distribution_uniform
    use stdlib_kinds
    use stdlib_error, only : error_stop
    use stdlib_stats_distribution_PRNG, only : dist_rand

    implicit none
    private

    real(dp), parameter  :: MESENNE_NUMBER = 1.0_dp / (2.0_dp ** 53 - 1.0_dp)
    integer(int64), parameter :: INT_ONE = 1_int64

    public :: rvs_uniform
    public :: pdf_uniform
    public :: cdf_uniform
    public :: shuffle


    interface rvs_uniform
    !! Version experimental
    !!
    !! Get uniformly distributed random variate for integer, real and complex
    !! variables.
    !! ([Specification](../page/specs/stdlib_stats_distribution_uniform.html#
    !! description))

        module procedure unif_dist_rvs_0_rsp                 ! 0 dummy variable

        module procedure unif_dist_rvs_1_iint8     ! 1 dummy variable
        module procedure unif_dist_rvs_1_iint16     ! 1 dummy variable
        module procedure unif_dist_rvs_1_iint32     ! 1 dummy variable
        module procedure unif_dist_rvs_1_iint64     ! 1 dummy variable
        module procedure unif_dist_rvs_1_rsp     ! 1 dummy variable
        module procedure unif_dist_rvs_1_rdp     ! 1 dummy variable
        module procedure unif_dist_rvs_1_rqp     ! 1 dummy variable
        module procedure unif_dist_rvs_1_csp     ! 1 dummy variable
        module procedure unif_dist_rvs_1_cdp     ! 1 dummy variable
        module procedure unif_dist_rvs_1_cqp     ! 1 dummy variable

        module procedure unif_dist_rvs_iint8       ! 2 dummy variables
        module procedure unif_dist_rvs_iint16       ! 2 dummy variables
        module procedure unif_dist_rvs_iint32       ! 2 dummy variables
        module procedure unif_dist_rvs_iint64       ! 2 dummy variables
        module procedure unif_dist_rvs_rsp       ! 2 dummy variables
        module procedure unif_dist_rvs_rdp       ! 2 dummy variables
        module procedure unif_dist_rvs_rqp       ! 2 dummy variables
        module procedure unif_dist_rvs_csp       ! 2 dummy variables
        module procedure unif_dist_rvs_cdp       ! 2 dummy variables
        module procedure unif_dist_rvs_cqp       ! 2 dummy variables

        module procedure unif_dist_rvs_array_iint8 ! 3 dummy variables
        module procedure unif_dist_rvs_array_iint16 ! 3 dummy variables
        module procedure unif_dist_rvs_array_iint32 ! 3 dummy variables
        module procedure unif_dist_rvs_array_iint64 ! 3 dummy variables
        module procedure unif_dist_rvs_array_rsp ! 3 dummy variables
        module procedure unif_dist_rvs_array_rdp ! 3 dummy variables
        module procedure unif_dist_rvs_array_rqp ! 3 dummy variables
        module procedure unif_dist_rvs_array_csp ! 3 dummy variables
        module procedure unif_dist_rvs_array_cdp ! 3 dummy variables
        module procedure unif_dist_rvs_array_cqp ! 3 dummy variables
    end interface rvs_uniform


    interface pdf_uniform
    !! Version experimental
    !!
    !! Get uniform distribution probability density (pdf) for integer, real and
    !! complex variables.
    !! ([Specification](../page/specs/stdlib_stats_distribution_uniform.html#
    !! description))

        module procedure unif_dist_pdf_iint8
        module procedure unif_dist_pdf_iint16
        module procedure unif_dist_pdf_iint32
        module procedure unif_dist_pdf_iint64
        module procedure unif_dist_pdf_rsp
        module procedure unif_dist_pdf_rdp
        module procedure unif_dist_pdf_rqp
        module procedure unif_dist_pdf_csp
        module procedure unif_dist_pdf_cdp
        module procedure unif_dist_pdf_cqp
    end interface pdf_uniform


    interface cdf_uniform
    !! Version experimental
    !!
    !! Get uniform distribution cumulative distribution function (cdf) for integer,
    !! real and complex variables.
    !! ([Specification](../page/specs/stdlib_stats_distribution_uniform.html#
    !! description))
    !!
        module procedure unif_dist_cdf_iint8
        module procedure unif_dist_cdf_iint16
        module procedure unif_dist_cdf_iint32
        module procedure unif_dist_cdf_iint64
        module procedure unif_dist_cdf_rsp
        module procedure unif_dist_cdf_rdp
        module procedure unif_dist_cdf_rqp
        module procedure unif_dist_cdf_csp
        module procedure unif_dist_cdf_cdp
        module procedure unif_dist_cdf_cqp
    end interface cdf_uniform


    interface shuffle
    !! Version experimental
    !!
    !! Fisher-Yates shuffle algorithm for a rank one array of integer, real and
    !! complex variables.
    !! ([Specification](../page/specs/stdlib_stats_distribution_uniform.html#
    !! description))
    !!
        module procedure shuffle_iint8
        module procedure shuffle_iint16
        module procedure shuffle_iint32
        module procedure shuffle_iint64
        module procedure shuffle_rsp
        module procedure shuffle_rdp
        module procedure shuffle_rqp
        module procedure shuffle_csp
        module procedure shuffle_cdp
        module procedure shuffle_cqp
    end interface shuffle






contains


    impure elemental function unif_dist_rvs_1_iint8(scale) result(res)
    !
    ! Uniformly distributed integer in [0, scale]
    ! Bitmask with rejection
    ! https://www.pcg-random.org/posts/bounded-rands.html
    !
    ! Fortran 90 translated from c by Jim-215-fisher
    !
        integer(int8), intent(in) :: scale
        integer(int8) ::  res, u, mask
        integer :: zeros, bits_left, bits

        if(scale <= 0_int8) call error_stop("Error(unif_dist_rvs_1): Uniform"&
            //" distribution scale parameter must be positive")
        zeros = leadz(scale)
        bits = bit_size(scale) - zeros
        mask = shiftr(not(0_int8), zeros)
        L1 : do
            u = dist_rand(scale)
            res = iand(u, mask)
            if(res <= scale) exit L1
            bits_left = zeros
            L2 : do
                if(bits_left < bits) exit L2
                u = shiftr(u, bits)
                res = iand(u, mask)
                if(res <= scale) exit L1
                bits_left = bits_left - bits
            end do L2
        end do L1
    end function unif_dist_rvs_1_iint8

    impure elemental function unif_dist_rvs_1_iint16(scale) result(res)
    !
    ! Uniformly distributed integer in [0, scale]
    ! Bitmask with rejection
    ! https://www.pcg-random.org/posts/bounded-rands.html
    !
    ! Fortran 90 translated from c by Jim-215-fisher
    !
        integer(int16), intent(in) :: scale
        integer(int16) ::  res, u, mask
        integer :: zeros, bits_left, bits

        if(scale <= 0_int16) call error_stop("Error(unif_dist_rvs_1): Uniform"&
            //" distribution scale parameter must be positive")
        zeros = leadz(scale)
        bits = bit_size(scale) - zeros
        mask = shiftr(not(0_int16), zeros)
        L1 : do
            u = dist_rand(scale)
            res = iand(u, mask)
            if(res <= scale) exit L1
            bits_left = zeros
            L2 : do
                if(bits_left < bits) exit L2
                u = shiftr(u, bits)
                res = iand(u, mask)
                if(res <= scale) exit L1
                bits_left = bits_left - bits
            end do L2
        end do L1
    end function unif_dist_rvs_1_iint16

    impure elemental function unif_dist_rvs_1_iint32(scale) result(res)
    !
    ! Uniformly distributed integer in [0, scale]
    ! Bitmask with rejection
    ! https://www.pcg-random.org/posts/bounded-rands.html
    !
    ! Fortran 90 translated from c by Jim-215-fisher
    !
        integer(int32), intent(in) :: scale
        integer(int32) ::  res, u, mask
        integer :: zeros, bits_left, bits

        if(scale <= 0_int32) call error_stop("Error(unif_dist_rvs_1): Uniform"&
            //" distribution scale parameter must be positive")
        zeros = leadz(scale)
        bits = bit_size(scale) - zeros
        mask = shiftr(not(0_int32), zeros)
        L1 : do
            u = dist_rand(scale)
            res = iand(u, mask)
            if(res <= scale) exit L1
            bits_left = zeros
            L2 : do
                if(bits_left < bits) exit L2
                u = shiftr(u, bits)
                res = iand(u, mask)
                if(res <= scale) exit L1
                bits_left = bits_left - bits
            end do L2
        end do L1
    end function unif_dist_rvs_1_iint32

    impure elemental function unif_dist_rvs_1_iint64(scale) result(res)
    !
    ! Uniformly distributed integer in [0, scale]
    ! Bitmask with rejection
    ! https://www.pcg-random.org/posts/bounded-rands.html
    !
    ! Fortran 90 translated from c by Jim-215-fisher
    !
        integer(int64), intent(in) :: scale
        integer(int64) ::  res, u, mask
        integer :: zeros, bits_left, bits

        if(scale <= 0_int64) call error_stop("Error(unif_dist_rvs_1): Uniform"&
            //" distribution scale parameter must be positive")
        zeros = leadz(scale)
        bits = bit_size(scale) - zeros
        mask = shiftr(not(0_int64), zeros)
        L1 : do
            u = dist_rand(scale)
            res = iand(u, mask)
            if(res <= scale) exit L1
            bits_left = zeros
            L2 : do
                if(bits_left < bits) exit L2
                u = shiftr(u, bits)
                res = iand(u, mask)
                if(res <= scale) exit L1
                bits_left = bits_left - bits
            end do L2
        end do L1
    end function unif_dist_rvs_1_iint64




    impure elemental function unif_dist_rvs_iint8(loc, scale)        &
        result( res )
    !
    ! Uniformly distributed integer in [loc, loc + scale]
    !
        integer(int8), intent(in) :: loc, scale
        integer(int8)  ::  res

        if(scale <= 0_int8) call error_stop("Error(unif_dist_rvs): Uniform"  &
            //" distribution scale parameter must be positive")
        res = loc + unif_dist_rvs_1_iint8(scale)
    end function unif_dist_rvs_iint8

    impure elemental function unif_dist_rvs_iint16(loc, scale)        &
        result( res )
    !
    ! Uniformly distributed integer in [loc, loc + scale]
    !
        integer(int16), intent(in) :: loc, scale
        integer(int16)  ::  res

        if(scale <= 0_int16) call error_stop("Error(unif_dist_rvs): Uniform"  &
            //" distribution scale parameter must be positive")
        res = loc + unif_dist_rvs_1_iint16(scale)
    end function unif_dist_rvs_iint16

    impure elemental function unif_dist_rvs_iint32(loc, scale)        &
        result( res )
    !
    ! Uniformly distributed integer in [loc, loc + scale]
    !
        integer(int32), intent(in) :: loc, scale
        integer(int32)  ::  res

        if(scale <= 0_int32) call error_stop("Error(unif_dist_rvs): Uniform"  &
            //" distribution scale parameter must be positive")
        res = loc + unif_dist_rvs_1_iint32(scale)
    end function unif_dist_rvs_iint32

    impure elemental function unif_dist_rvs_iint64(loc, scale)        &
        result( res )
    !
    ! Uniformly distributed integer in [loc, loc + scale]
    !
        integer(int64), intent(in) :: loc, scale
        integer(int64)  ::  res

        if(scale <= 0_int64) call error_stop("Error(unif_dist_rvs): Uniform"  &
            //" distribution scale parameter must be positive")
        res = loc + unif_dist_rvs_1_iint64(scale)
    end function unif_dist_rvs_iint64




    impure elemental function unif_dist_rvs_0_rsp( ) result(res)
    !
    ! Uniformly distributed float in [0,1]
    ! Based on the paper by Frederic Goualard, "Generating Random Floating-
    ! Point Numbers By Dividing Integers: a Case Study", Proceedings of
    ! ICCS 2020, June 2020, Amsterdam, Netherlands
    !
        real(sp)  ::  res
        integer(int64) :: tmp

        tmp = shiftr(dist_rand(INT_ONE), 11)        ! Get random from [0,2^53-1]
        res = real(tmp * MESENNE_NUMBER, kind = sp) ! convert to [0,1]
    end function unif_dist_rvs_0_rsp

    impure elemental function unif_dist_rvs_0_rdp( ) result(res)
    !
    ! Uniformly distributed float in [0,1]
    ! Based on the paper by Frederic Goualard, "Generating Random Floating-
    ! Point Numbers By Dividing Integers: a Case Study", Proceedings of
    ! ICCS 2020, June 2020, Amsterdam, Netherlands
    !
        real(dp)  ::  res
        integer(int64) :: tmp

        tmp = shiftr(dist_rand(INT_ONE), 11)        ! Get random from [0,2^53-1]
        res = real(tmp * MESENNE_NUMBER, kind = dp) ! convert to [0,1]
    end function unif_dist_rvs_0_rdp

    impure elemental function unif_dist_rvs_0_rqp( ) result(res)
    !
    ! Uniformly distributed float in [0,1]
    ! Based on the paper by Frederic Goualard, "Generating Random Floating-
    ! Point Numbers By Dividing Integers: a Case Study", Proceedings of
    ! ICCS 2020, June 2020, Amsterdam, Netherlands
    !
        real(qp)  ::  res
        integer(int64) :: tmp

        tmp = shiftr(dist_rand(INT_ONE), 11)        ! Get random from [0,2^53-1]
        res = real(tmp * MESENNE_NUMBER, kind = qp) ! convert to [0,1]
    end function unif_dist_rvs_0_rqp




    impure elemental function unif_dist_rvs_1_rsp(scale) result(res)
    !
    ! Uniformly distributed float in [0, scale]
    !
        real(sp), intent(in) :: scale
        real(sp)  ::  res

        if(scale == 0._sp) call error_stop("Error(unif_dist_rvs_1): "      &
            //"Uniform distribution scale parameter must be non-zero")
        res = scale * unif_dist_rvs_0_rsp( )
    end function unif_dist_rvs_1_rsp

    impure elemental function unif_dist_rvs_1_rdp(scale) result(res)
    !
    ! Uniformly distributed float in [0, scale]
    !
        real(dp), intent(in) :: scale
        real(dp)  ::  res

        if(scale == 0._dp) call error_stop("Error(unif_dist_rvs_1): "      &
            //"Uniform distribution scale parameter must be non-zero")
        res = scale * unif_dist_rvs_0_rdp( )
    end function unif_dist_rvs_1_rdp

    impure elemental function unif_dist_rvs_1_rqp(scale) result(res)
    !
    ! Uniformly distributed float in [0, scale]
    !
        real(qp), intent(in) :: scale
        real(qp)  ::  res

        if(scale == 0._qp) call error_stop("Error(unif_dist_rvs_1): "      &
            //"Uniform distribution scale parameter must be non-zero")
        res = scale * unif_dist_rvs_0_rqp( )
    end function unif_dist_rvs_1_rqp




    impure elemental function unif_dist_rvs_rsp(loc, scale)        &
        result(res)
    !
    ! Uniformly distributed float in [loc, loc + scale]
    !
        real(sp), intent(in) :: loc, scale
        real(sp)  ::  res

        if(scale == 0._sp) call error_stop("Error(unif_dist_rvs): "        &
           //"Uniform distribution scale parameter must be non-zero")
        res = loc + scale * unif_dist_rvs_0_rsp( )
    end function unif_dist_rvs_rsp

    impure elemental function unif_dist_rvs_rdp(loc, scale)        &
        result(res)
    !
    ! Uniformly distributed float in [loc, loc + scale]
    !
        real(dp), intent(in) :: loc, scale
        real(dp)  ::  res

        if(scale == 0._dp) call error_stop("Error(unif_dist_rvs): "        &
           //"Uniform distribution scale parameter must be non-zero")
        res = loc + scale * unif_dist_rvs_0_rdp( )
    end function unif_dist_rvs_rdp

    impure elemental function unif_dist_rvs_rqp(loc, scale)        &
        result(res)
    !
    ! Uniformly distributed float in [loc, loc + scale]
    !
        real(qp), intent(in) :: loc, scale
        real(qp)  ::  res

        if(scale == 0._qp) call error_stop("Error(unif_dist_rvs): "        &
           //"Uniform distribution scale parameter must be non-zero")
        res = loc + scale * unif_dist_rvs_0_rqp( )
    end function unif_dist_rvs_rqp




    impure elemental function unif_dist_rvs_1_csp(scale)           &
        result(res)
    !
    ! Uniformly distributed complex in [(0,0i), (scale, i(scale))]
    ! The real part and imaginary part are independent of each other, so that
    ! the joint distribution is on an unit square [(0,0i), (scale,i(scale))]
    !
        complex(sp), intent(in) :: scale
        complex(sp) :: res
        real(sp) :: r1, tr, ti

        if(scale == (0.0_sp, 0.0_sp)) call error_stop("Error(uni_dist_"&
           //"rvs_1): Uniform distribution scale parameter must be non-zero")
        r1 = unif_dist_rvs_0_rsp( )
        if(scale % re == 0.0_sp) then
            ti = scale % im * r1
            tr = 0.0_sp
        else if(scale % im == 0.0_sp) then
            tr = scale % re * r1
            ti = 0.0_sp
        else
            tr = scale % re * r1
            r1 = unif_dist_rvs_0_rsp( )
            ti = scale % im * r1
        end if
        res = cmplx(tr, ti, kind=sp)
    end function unif_dist_rvs_1_csp

    impure elemental function unif_dist_rvs_1_cdp(scale)           &
        result(res)
    !
    ! Uniformly distributed complex in [(0,0i), (scale, i(scale))]
    ! The real part and imaginary part are independent of each other, so that
    ! the joint distribution is on an unit square [(0,0i), (scale,i(scale))]
    !
        complex(dp), intent(in) :: scale
        complex(dp) :: res
        real(dp) :: r1, tr, ti

        if(scale == (0.0_dp, 0.0_dp)) call error_stop("Error(uni_dist_"&
           //"rvs_1): Uniform distribution scale parameter must be non-zero")
        r1 = unif_dist_rvs_0_rdp( )
        if(scale % re == 0.0_dp) then
            ti = scale % im * r1
            tr = 0.0_dp
        else if(scale % im == 0.0_dp) then
            tr = scale % re * r1
            ti = 0.0_dp
        else
            tr = scale % re * r1
            r1 = unif_dist_rvs_0_rdp( )
            ti = scale % im * r1
        end if
        res = cmplx(tr, ti, kind=dp)
    end function unif_dist_rvs_1_cdp

    impure elemental function unif_dist_rvs_1_cqp(scale)           &
        result(res)
    !
    ! Uniformly distributed complex in [(0,0i), (scale, i(scale))]
    ! The real part and imaginary part are independent of each other, so that
    ! the joint distribution is on an unit square [(0,0i), (scale,i(scale))]
    !
        complex(qp), intent(in) :: scale
        complex(qp) :: res
        real(qp) :: r1, tr, ti

        if(scale == (0.0_qp, 0.0_qp)) call error_stop("Error(uni_dist_"&
           //"rvs_1): Uniform distribution scale parameter must be non-zero")
        r1 = unif_dist_rvs_0_rqp( )
        if(scale % re == 0.0_qp) then
            ti = scale % im * r1
            tr = 0.0_qp
        else if(scale % im == 0.0_qp) then
            tr = scale % re * r1
            ti = 0.0_qp
        else
            tr = scale % re * r1
            r1 = unif_dist_rvs_0_rqp( )
            ti = scale % im * r1
        end if
        res = cmplx(tr, ti, kind=qp)
    end function unif_dist_rvs_1_cqp




    impure elemental function unif_dist_rvs_csp(loc, scale)        &
        result(res)
    !
    ! Uniformly distributed complex in [(loc,iloc), (loc + scale, i(loc +
    ! scale))].
    ! The real part and imaginary part are independent of each other, so that
    ! the joint distribution is on an unit square [(loc,iloc), (loc + scale,
    ! i(loc + scale))]
    !
        complex(sp), intent(in) :: loc, scale
        complex(sp) :: res
        real(sp) :: r1, tr, ti

        if(scale == (0.0_sp, 0.0_sp)) call error_stop("Error(uni_dist_"&
            //"rvs): Uniform distribution scale parameter must be non-zero")
        r1 = unif_dist_rvs_0_rsp( )
        if(scale % re == 0.0_sp) then
            tr = loc % re
            ti = loc % im + scale % im * r1
        else if(scale % im == 0.0_sp) then
            tr = loc % re + scale % re * r1
            ti = loc % im
        else
            tr = loc % re + scale % re * r1
            r1 = unif_dist_rvs_0_rsp( )
            ti = loc % im + scale % im * r1
        end if
        res = cmplx(tr, ti, kind=sp)
    end function unif_dist_rvs_csp

    impure elemental function unif_dist_rvs_cdp(loc, scale)        &
        result(res)
    !
    ! Uniformly distributed complex in [(loc,iloc), (loc + scale, i(loc +
    ! scale))].
    ! The real part and imaginary part are independent of each other, so that
    ! the joint distribution is on an unit square [(loc,iloc), (loc + scale,
    ! i(loc + scale))]
    !
        complex(dp), intent(in) :: loc, scale
        complex(dp) :: res
        real(dp) :: r1, tr, ti

        if(scale == (0.0_dp, 0.0_dp)) call error_stop("Error(uni_dist_"&
            //"rvs): Uniform distribution scale parameter must be non-zero")
        r1 = unif_dist_rvs_0_rdp( )
        if(scale % re == 0.0_dp) then
            tr = loc % re
            ti = loc % im + scale % im * r1
        else if(scale % im == 0.0_dp) then
            tr = loc % re + scale % re * r1
            ti = loc % im
        else
            tr = loc % re + scale % re * r1
            r1 = unif_dist_rvs_0_rdp( )
            ti = loc % im + scale % im * r1
        end if
        res = cmplx(tr, ti, kind=dp)
    end function unif_dist_rvs_cdp

    impure elemental function unif_dist_rvs_cqp(loc, scale)        &
        result(res)
    !
    ! Uniformly distributed complex in [(loc,iloc), (loc + scale, i(loc +
    ! scale))].
    ! The real part and imaginary part are independent of each other, so that
    ! the joint distribution is on an unit square [(loc,iloc), (loc + scale,
    ! i(loc + scale))]
    !
        complex(qp), intent(in) :: loc, scale
        complex(qp) :: res
        real(qp) :: r1, tr, ti

        if(scale == (0.0_qp, 0.0_qp)) call error_stop("Error(uni_dist_"&
            //"rvs): Uniform distribution scale parameter must be non-zero")
        r1 = unif_dist_rvs_0_rqp( )
        if(scale % re == 0.0_qp) then
            tr = loc % re
            ti = loc % im + scale % im * r1
        else if(scale % im == 0.0_qp) then
            tr = loc % re + scale % re * r1
            ti = loc % im
        else
            tr = loc % re + scale % re * r1
            r1 = unif_dist_rvs_0_rqp( )
            ti = loc % im + scale % im * r1
        end if
        res = cmplx(tr, ti, kind=qp)
    end function unif_dist_rvs_cqp




    function unif_dist_rvs_array_iint8(loc, scale, array_size)       &
        result(res)

        integer, intent(in) :: array_size
        integer(int8), intent(in) :: loc, scale
        integer(int8) :: res(array_size)
        integer(int8) :: u, mask, nn
        integer :: i, zeros, bits_left, bits

        if(scale == 0_int8) call error_stop("Error(unif_dist_rvs_array): "   &
            //"Uniform distribution scale parameter must be non-zero")
        zeros = leadz(scale)
        bits = bit_size(scale) - zeros
        mask = shiftr(not(0_int8), zeros)
        do i = 1, array_size
            L1 : do
                u = dist_rand(scale)
                nn = iand(u, mask)
                if(nn <= scale) exit L1
                bits_left = zeros
                L2 : do
                    if(bits_left < bits) exit L2
                    u = shiftr(u, bits)
                    nn = iand(u, mask)
                    if(nn <= scale) exit L1
                    bits_left = bits_left - bits
                end do L2
            end do L1
            res(i) = loc + nn
        end do
    end function unif_dist_rvs_array_iint8

    function unif_dist_rvs_array_iint16(loc, scale, array_size)       &
        result(res)

        integer, intent(in) :: array_size
        integer(int16), intent(in) :: loc, scale
        integer(int16) :: res(array_size)
        integer(int16) :: u, mask, nn
        integer :: i, zeros, bits_left, bits

        if(scale == 0_int16) call error_stop("Error(unif_dist_rvs_array): "   &
            //"Uniform distribution scale parameter must be non-zero")
        zeros = leadz(scale)
        bits = bit_size(scale) - zeros
        mask = shiftr(not(0_int16), zeros)
        do i = 1, array_size
            L1 : do
                u = dist_rand(scale)
                nn = iand(u, mask)
                if(nn <= scale) exit L1
                bits_left = zeros
                L2 : do
                    if(bits_left < bits) exit L2
                    u = shiftr(u, bits)
                    nn = iand(u, mask)
                    if(nn <= scale) exit L1
                    bits_left = bits_left - bits
                end do L2
            end do L1
            res(i) = loc + nn
        end do
    end function unif_dist_rvs_array_iint16

    function unif_dist_rvs_array_iint32(loc, scale, array_size)       &
        result(res)

        integer, intent(in) :: array_size
        integer(int32), intent(in) :: loc, scale
        integer(int32) :: res(array_size)
        integer(int32) :: u, mask, nn
        integer :: i, zeros, bits_left, bits

        if(scale == 0_int32) call error_stop("Error(unif_dist_rvs_array): "   &
            //"Uniform distribution scale parameter must be non-zero")
        zeros = leadz(scale)
        bits = bit_size(scale) - zeros
        mask = shiftr(not(0_int32), zeros)
        do i = 1, array_size
            L1 : do
                u = dist_rand(scale)
                nn = iand(u, mask)
                if(nn <= scale) exit L1
                bits_left = zeros
                L2 : do
                    if(bits_left < bits) exit L2
                    u = shiftr(u, bits)
                    nn = iand(u, mask)
                    if(nn <= scale) exit L1
                    bits_left = bits_left - bits
                end do L2
            end do L1
            res(i) = loc + nn
        end do
    end function unif_dist_rvs_array_iint32

    function unif_dist_rvs_array_iint64(loc, scale, array_size)       &
        result(res)

        integer, intent(in) :: array_size
        integer(int64), intent(in) :: loc, scale
        integer(int64) :: res(array_size)
        integer(int64) :: u, mask, nn
        integer :: i, zeros, bits_left, bits

        if(scale == 0_int64) call error_stop("Error(unif_dist_rvs_array): "   &
            //"Uniform distribution scale parameter must be non-zero")
        zeros = leadz(scale)
        bits = bit_size(scale) - zeros
        mask = shiftr(not(0_int64), zeros)
        do i = 1, array_size
            L1 : do
                u = dist_rand(scale)
                nn = iand(u, mask)
                if(nn <= scale) exit L1
                bits_left = zeros
                L2 : do
                    if(bits_left < bits) exit L2
                    u = shiftr(u, bits)
                    nn = iand(u, mask)
                    if(nn <= scale) exit L1
                    bits_left = bits_left - bits
                end do L2
            end do L1
            res(i) = loc + nn
        end do
    end function unif_dist_rvs_array_iint64




    function unif_dist_rvs_array_rsp(loc, scale, array_size)       &
        result(res)

        integer, intent(in) :: array_size
        real(sp), intent(in) :: loc, scale
        real(sp) :: res(array_size)
        real(sp) :: t
        integer(int64) :: tmp
        integer :: i


        if(scale == 0._sp) call error_stop("Error(unif_dist_rvs_array):"   &
           //" Uniform distribution scale parameter must be non-zero")
        do i = 1, array_size
            tmp = shiftr(dist_rand(INT_ONE), 11)
            t = real(tmp * MESENNE_NUMBER, kind = sp)
            res(i) = loc + scale * t
        end do
    end function unif_dist_rvs_array_rsp

    function unif_dist_rvs_array_rdp(loc, scale, array_size)       &
        result(res)

        integer, intent(in) :: array_size
        real(dp), intent(in) :: loc, scale
        real(dp) :: res(array_size)
        real(dp) :: t
        integer(int64) :: tmp
        integer :: i


        if(scale == 0._dp) call error_stop("Error(unif_dist_rvs_array):"   &
           //" Uniform distribution scale parameter must be non-zero")
        do i = 1, array_size
            tmp = shiftr(dist_rand(INT_ONE), 11)
            t = real(tmp * MESENNE_NUMBER, kind = dp)
            res(i) = loc + scale * t
        end do
    end function unif_dist_rvs_array_rdp

    function unif_dist_rvs_array_rqp(loc, scale, array_size)       &
        result(res)

        integer, intent(in) :: array_size
        real(qp), intent(in) :: loc, scale
        real(qp) :: res(array_size)
        real(qp) :: t
        integer(int64) :: tmp
        integer :: i


        if(scale == 0._qp) call error_stop("Error(unif_dist_rvs_array):"   &
           //" Uniform distribution scale parameter must be non-zero")
        do i = 1, array_size
            tmp = shiftr(dist_rand(INT_ONE), 11)
            t = real(tmp * MESENNE_NUMBER, kind = qp)
            res(i) = loc + scale * t
        end do
    end function unif_dist_rvs_array_rqp




    function unif_dist_rvs_array_csp(loc, scale, array_size)       &
        result(res)

        integer, intent(in) :: array_size
        complex(sp), intent(in) :: loc, scale
        complex(sp) :: res(array_size)
        real(sp) :: r1, tr, ti
        integer(int64) :: tmp
        integer :: i


        if(scale == (0.0_sp, 0.0_sp)) call error_stop("Error(unif_dist"&
           //"_rvs_array): Uniform distribution scale parameter must be "      &
           //"non-zero")
        do i = 1, array_size
            tmp = shiftr(dist_rand(INT_ONE), 11)
            r1 = real(tmp * MESENNE_NUMBER, kind = sp)
            if(scale % re == 0.0_sp) then
                tr = loc % re
                ti = loc % im + scale % im * r1
            else if(scale % im == 0.0_sp) then
                tr = loc % re + scale % re * r1
                ti = loc % im
            else
                tr = loc % re + scale % re * r1
                tmp = shiftr(dist_rand(INT_ONE), 11)
                r1 = real(tmp * MESENNE_NUMBER, kind = sp)
                ti = loc % im + scale % im * r1
            end if
            res(i) = cmplx(tr, ti, kind=sp)
        end do
    end function unif_dist_rvs_array_csp

    function unif_dist_rvs_array_cdp(loc, scale, array_size)       &
        result(res)

        integer, intent(in) :: array_size
        complex(dp), intent(in) :: loc, scale
        complex(dp) :: res(array_size)
        real(dp) :: r1, tr, ti
        integer(int64) :: tmp
        integer :: i


        if(scale == (0.0_dp, 0.0_dp)) call error_stop("Error(unif_dist"&
           //"_rvs_array): Uniform distribution scale parameter must be "      &
           //"non-zero")
        do i = 1, array_size
            tmp = shiftr(dist_rand(INT_ONE), 11)
            r1 = real(tmp * MESENNE_NUMBER, kind = dp)
            if(scale % re == 0.0_dp) then
                tr = loc % re
                ti = loc % im + scale % im * r1
            else if(scale % im == 0.0_dp) then
                tr = loc % re + scale % re * r1
                ti = loc % im
            else
                tr = loc % re + scale % re * r1
                tmp = shiftr(dist_rand(INT_ONE), 11)
                r1 = real(tmp * MESENNE_NUMBER, kind = dp)
                ti = loc % im + scale % im * r1
            end if
            res(i) = cmplx(tr, ti, kind=dp)
        end do
    end function unif_dist_rvs_array_cdp

    function unif_dist_rvs_array_cqp(loc, scale, array_size)       &
        result(res)

        integer, intent(in) :: array_size
        complex(qp), intent(in) :: loc, scale
        complex(qp) :: res(array_size)
        real(qp) :: r1, tr, ti
        integer(int64) :: tmp
        integer :: i


        if(scale == (0.0_qp, 0.0_qp)) call error_stop("Error(unif_dist"&
           //"_rvs_array): Uniform distribution scale parameter must be "      &
           //"non-zero")
        do i = 1, array_size
            tmp = shiftr(dist_rand(INT_ONE), 11)
            r1 = real(tmp * MESENNE_NUMBER, kind = qp)
            if(scale % re == 0.0_qp) then
                tr = loc % re
                ti = loc % im + scale % im * r1
            else if(scale % im == 0.0_qp) then
                tr = loc % re + scale % re * r1
                ti = loc % im
            else
                tr = loc % re + scale % re * r1
                tmp = shiftr(dist_rand(INT_ONE), 11)
                r1 = real(tmp * MESENNE_NUMBER, kind = qp)
                ti = loc % im + scale % im * r1
            end if
            res(i) = cmplx(tr, ti, kind=qp)
        end do
    end function unif_dist_rvs_array_cqp




    elemental function unif_dist_pdf_iint8(x, loc, scale) result(res)

        integer(int8), intent(in) :: x, loc, scale
        real :: res

        if(scale == 0_int8) then
            res = 0.0
        else if(x < loc .or. x > (loc + scale)) then
            res = 0.0
        else
            res = 1. / (scale + 1_int8)
        end if
    end function unif_dist_pdf_iint8

    elemental function unif_dist_pdf_iint16(x, loc, scale) result(res)

        integer(int16), intent(in) :: x, loc, scale
        real :: res

        if(scale == 0_int16) then
            res = 0.0
        else if(x < loc .or. x > (loc + scale)) then
            res = 0.0
        else
            res = 1. / (scale + 1_int16)
        end if
    end function unif_dist_pdf_iint16

    elemental function unif_dist_pdf_iint32(x, loc, scale) result(res)

        integer(int32), intent(in) :: x, loc, scale
        real :: res

        if(scale == 0_int32) then
            res = 0.0
        else if(x < loc .or. x > (loc + scale)) then
            res = 0.0
        else
            res = 1. / (scale + 1_int32)
        end if
    end function unif_dist_pdf_iint32

    elemental function unif_dist_pdf_iint64(x, loc, scale) result(res)

        integer(int64), intent(in) :: x, loc, scale
        real :: res

        if(scale == 0_int64) then
            res = 0.0
        else if(x < loc .or. x > (loc + scale)) then
            res = 0.0
        else
            res = 1. / (scale + 1_int64)
        end if
    end function unif_dist_pdf_iint64




    elemental function unif_dist_pdf_rsp(x, loc, scale) result(res)

        real(sp), intent(in) :: x, loc, scale
        real :: res

        if(scale == 0.0_sp) then
            res = 0.0
        else if(x < loc .or. x > (loc + scale)) then
            res = 0.0
        else
            res = 1.0 / scale
        end if
    end function unif_dist_pdf_rsp

    elemental function unif_dist_pdf_rdp(x, loc, scale) result(res)

        real(dp), intent(in) :: x, loc, scale
        real :: res

        if(scale == 0.0_dp) then
            res = 0.0
        else if(x < loc .or. x > (loc + scale)) then
            res = 0.0
        else
            res = 1.0 / scale
        end if
    end function unif_dist_pdf_rdp

    elemental function unif_dist_pdf_rqp(x, loc, scale) result(res)

        real(qp), intent(in) :: x, loc, scale
        real :: res

        if(scale == 0.0_qp) then
            res = 0.0
        else if(x < loc .or. x > (loc + scale)) then
            res = 0.0
        else
            res = 1.0 / scale
        end if
    end function unif_dist_pdf_rqp




    elemental function unif_dist_pdf_csp(x, loc, scale) result(res)

        complex(sp), intent(in) :: x, loc, scale
        real :: res
        real(sp) :: tr, ti

        tr = loc % re + scale % re; ti = loc % im + scale % im
        if(scale == (0.0_sp,0.0_sp)) then
            res = 0.0
        else if((x % re >= loc % re .and. x % re <= tr) .and.                  &
            (x % im >= loc % im .and. x % im <= ti)) then
            res = 1.0 / (scale % re * scale % im)
        else
            res = 0.0
        end if
    end function unif_dist_pdf_csp

    elemental function unif_dist_pdf_cdp(x, loc, scale) result(res)

        complex(dp), intent(in) :: x, loc, scale
        real :: res
        real(dp) :: tr, ti

        tr = loc % re + scale % re; ti = loc % im + scale % im
        if(scale == (0.0_dp,0.0_dp)) then
            res = 0.0
        else if((x % re >= loc % re .and. x % re <= tr) .and.                  &
            (x % im >= loc % im .and. x % im <= ti)) then
            res = 1.0 / (scale % re * scale % im)
        else
            res = 0.0
        end if
    end function unif_dist_pdf_cdp

    elemental function unif_dist_pdf_cqp(x, loc, scale) result(res)

        complex(qp), intent(in) :: x, loc, scale
        real :: res
        real(qp) :: tr, ti

        tr = loc % re + scale % re; ti = loc % im + scale % im
        if(scale == (0.0_qp,0.0_qp)) then
            res = 0.0
        else if((x % re >= loc % re .and. x % re <= tr) .and.                  &
            (x % im >= loc % im .and. x % im <= ti)) then
            res = 1.0 / (scale % re * scale % im)
        else
            res = 0.0
        end if
    end function unif_dist_pdf_cqp




    elemental function unif_dist_cdf_iint8(x, loc, scale) result(res)

        integer(int8), intent(in) :: x, loc, scale
        real :: res

        if(scale == 0_int8) then
            res = 0.0
        else if(x < loc) then
            res = 0.0
        else if(x >= loc .and. x <= (loc + scale)) then
            res = real((x - loc + 1_int8)) / real((scale + 1_int8))
        else
            res = 1.0
        end if
    end function unif_dist_cdf_iint8

    elemental function unif_dist_cdf_iint16(x, loc, scale) result(res)

        integer(int16), intent(in) :: x, loc, scale
        real :: res

        if(scale == 0_int16) then
            res = 0.0
        else if(x < loc) then
            res = 0.0
        else if(x >= loc .and. x <= (loc + scale)) then
            res = real((x - loc + 1_int16)) / real((scale + 1_int16))
        else
            res = 1.0
        end if
    end function unif_dist_cdf_iint16

    elemental function unif_dist_cdf_iint32(x, loc, scale) result(res)

        integer(int32), intent(in) :: x, loc, scale
        real :: res

        if(scale == 0_int32) then
            res = 0.0
        else if(x < loc) then
            res = 0.0
        else if(x >= loc .and. x <= (loc + scale)) then
            res = real((x - loc + 1_int32)) / real((scale + 1_int32))
        else
            res = 1.0
        end if
    end function unif_dist_cdf_iint32

    elemental function unif_dist_cdf_iint64(x, loc, scale) result(res)

        integer(int64), intent(in) :: x, loc, scale
        real :: res

        if(scale == 0_int64) then
            res = 0.0
        else if(x < loc) then
            res = 0.0
        else if(x >= loc .and. x <= (loc + scale)) then
            res = real((x - loc + 1_int64)) / real((scale + 1_int64))
        else
            res = 1.0
        end if
    end function unif_dist_cdf_iint64




    elemental function unif_dist_cdf_rsp(x, loc, scale) result(res)

        real(sp), intent(in) :: x, loc, scale
        real :: res

        if(scale == 0.0_sp) then
            res = 0.0
        else if(x < loc) then
            res = 0.0
        else if(x >= loc .and. x <= (loc + scale)) then
            res = (x - loc) / scale
        else
            res = 1.0
        end if
    end function unif_dist_cdf_rsp

    elemental function unif_dist_cdf_rdp(x, loc, scale) result(res)

        real(dp), intent(in) :: x, loc, scale
        real :: res

        if(scale == 0.0_dp) then
            res = 0.0
        else if(x < loc) then
            res = 0.0
        else if(x >= loc .and. x <= (loc + scale)) then
            res = (x - loc) / scale
        else
            res = 1.0
        end if
    end function unif_dist_cdf_rdp

    elemental function unif_dist_cdf_rqp(x, loc, scale) result(res)

        real(qp), intent(in) :: x, loc, scale
        real :: res

        if(scale == 0.0_qp) then
            res = 0.0
        else if(x < loc) then
            res = 0.0
        else if(x >= loc .and. x <= (loc + scale)) then
            res = (x - loc) / scale
        else
            res = 1.0
        end if
    end function unif_dist_cdf_rqp




    elemental function unif_dist_cdf_csp(x, loc, scale) result(res)

        complex(sp), intent(in) :: x, loc, scale
        real :: res
        logical :: r1, r2, i1, i2

        if(scale == (0.0_sp,0.0_sp)) then
            res = 0.0
            return
        end if
        r1 = x % re < loc % re
        r2 = x % re > (loc % re + scale % re)
        i1 = x % im < loc % im
        i2 = x % im > (loc % im + scale % im)
        if(r1 .or. i1) then
            res = 0.0
        else if((.not. r1) .and. (.not. r2) .and. i2) then
            res = (x % re - loc % re) / scale % re
        else if((.not. i1) .and. (.not. i2) .and. r2) then
            res = (x % im - loc % im) / scale % im
        else if((.not. r1) .and. (.not. r2) .and. (.not. i1) .and. (.not. i2)) &
            then
            res = (x % re - loc % re) * (x % im - loc % im) /                  &
                   (scale % re * scale % im)
        else if(r2 .and. i2)then
             res = 1.0
        end if
    end function unif_dist_cdf_csp

    elemental function unif_dist_cdf_cdp(x, loc, scale) result(res)

        complex(dp), intent(in) :: x, loc, scale
        real :: res
        logical :: r1, r2, i1, i2

        if(scale == (0.0_dp,0.0_dp)) then
            res = 0.0
            return
        end if
        r1 = x % re < loc % re
        r2 = x % re > (loc % re + scale % re)
        i1 = x % im < loc % im
        i2 = x % im > (loc % im + scale % im)
        if(r1 .or. i1) then
            res = 0.0
        else if((.not. r1) .and. (.not. r2) .and. i2) then
            res = (x % re - loc % re) / scale % re
        else if((.not. i1) .and. (.not. i2) .and. r2) then
            res = (x % im - loc % im) / scale % im
        else if((.not. r1) .and. (.not. r2) .and. (.not. i1) .and. (.not. i2)) &
            then
            res = (x % re - loc % re) * (x % im - loc % im) /                  &
                   (scale % re * scale % im)
        else if(r2 .and. i2)then
             res = 1.0
        end if
    end function unif_dist_cdf_cdp

    elemental function unif_dist_cdf_cqp(x, loc, scale) result(res)

        complex(qp), intent(in) :: x, loc, scale
        real :: res
        logical :: r1, r2, i1, i2

        if(scale == (0.0_qp,0.0_qp)) then
            res = 0.0
            return
        end if
        r1 = x % re < loc % re
        r2 = x % re > (loc % re + scale % re)
        i1 = x % im < loc % im
        i2 = x % im > (loc % im + scale % im)
        if(r1 .or. i1) then
            res = 0.0
        else if((.not. r1) .and. (.not. r2) .and. i2) then
            res = (x % re - loc % re) / scale % re
        else if((.not. i1) .and. (.not. i2) .and. r2) then
            res = (x % im - loc % im) / scale % im
        else if((.not. r1) .and. (.not. r2) .and. (.not. i1) .and. (.not. i2)) &
            then
            res = (x % re - loc % re) * (x % im - loc % im) /                  &
                   (scale % re * scale % im)
        else if(r2 .and. i2)then
             res = 1.0
        end if
    end function unif_dist_cdf_cqp




    function shuffle_iint8( list ) result(res)

        integer(int8), intent(in) :: list(:)
        integer(int8) :: res(size(list))
        integer(int8) :: tmp
        integer :: n, i, j

        n = size(list)
        res = list
        do i = 1, n - 1
            j = rvs_uniform(n - i) + i
            tmp = res(i)
            res(i) = res(j)
            res(j) = tmp
        end do
    end function shuffle_iint8

    function shuffle_iint16( list ) result(res)

        integer(int16), intent(in) :: list(:)
        integer(int16) :: res(size(list))
        integer(int16) :: tmp
        integer :: n, i, j

        n = size(list)
        res = list
        do i = 1, n - 1
            j = rvs_uniform(n - i) + i
            tmp = res(i)
            res(i) = res(j)
            res(j) = tmp
        end do
    end function shuffle_iint16

    function shuffle_iint32( list ) result(res)

        integer(int32), intent(in) :: list(:)
        integer(int32) :: res(size(list))
        integer(int32) :: tmp
        integer :: n, i, j

        n = size(list)
        res = list
        do i = 1, n - 1
            j = rvs_uniform(n - i) + i
            tmp = res(i)
            res(i) = res(j)
            res(j) = tmp
        end do
    end function shuffle_iint32

    function shuffle_iint64( list ) result(res)

        integer(int64), intent(in) :: list(:)
        integer(int64) :: res(size(list))
        integer(int64) :: tmp
        integer :: n, i, j

        n = size(list)
        res = list
        do i = 1, n - 1
            j = rvs_uniform(n - i) + i
            tmp = res(i)
            res(i) = res(j)
            res(j) = tmp
        end do
    end function shuffle_iint64

    function shuffle_rsp( list ) result(res)

        real(sp), intent(in) :: list(:)
        real(sp) :: res(size(list))
        real(sp) :: tmp
        integer :: n, i, j

        n = size(list)
        res = list
        do i = 1, n - 1
            j = rvs_uniform(n - i) + i
            tmp = res(i)
            res(i) = res(j)
            res(j) = tmp
        end do
    end function shuffle_rsp

    function shuffle_rdp( list ) result(res)

        real(dp), intent(in) :: list(:)
        real(dp) :: res(size(list))
        real(dp) :: tmp
        integer :: n, i, j

        n = size(list)
        res = list
        do i = 1, n - 1
            j = rvs_uniform(n - i) + i
            tmp = res(i)
            res(i) = res(j)
            res(j) = tmp
        end do
    end function shuffle_rdp

    function shuffle_rqp( list ) result(res)

        real(qp), intent(in) :: list(:)
        real(qp) :: res(size(list))
        real(qp) :: tmp
        integer :: n, i, j

        n = size(list)
        res = list
        do i = 1, n - 1
            j = rvs_uniform(n - i) + i
            tmp = res(i)
            res(i) = res(j)
            res(j) = tmp
        end do
    end function shuffle_rqp

    function shuffle_csp( list ) result(res)

        complex(sp), intent(in) :: list(:)
        complex(sp) :: res(size(list))
        complex(sp) :: tmp
        integer :: n, i, j

        n = size(list)
        res = list
        do i = 1, n - 1
            j = rvs_uniform(n - i) + i
            tmp = res(i)
            res(i) = res(j)
            res(j) = tmp
        end do
    end function shuffle_csp

    function shuffle_cdp( list ) result(res)

        complex(dp), intent(in) :: list(:)
        complex(dp) :: res(size(list))
        complex(dp) :: tmp
        integer :: n, i, j

        n = size(list)
        res = list
        do i = 1, n - 1
            j = rvs_uniform(n - i) + i
            tmp = res(i)
            res(i) = res(j)
            res(j) = tmp
        end do
    end function shuffle_cdp

    function shuffle_cqp( list ) result(res)

        complex(qp), intent(in) :: list(:)
        complex(qp) :: res(size(list))
        complex(qp) :: tmp
        integer :: n, i, j

        n = size(list)
        res = list
        do i = 1, n - 1
            j = rvs_uniform(n - i) + i
            tmp = res(i)
            res(i) = res(j)
            res(j) = tmp
        end do
    end function shuffle_cqp

end module stdlib_stats_distribution_uniform
