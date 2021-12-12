module stdlib_stats_distribution_exponential
    use stdlib_kinds, only : sp, dp, xdp, qp, int32
    use stdlib_error, only : error_stop
    use stdlib_random, only : dist_rand
    use stdlib_stats_distribution_uniform, only : uni=>rvs_uniform

    implicit none
    private

    real(dp), parameter  :: ONE = 1.0_dp
    integer :: ke(0:255)
    real(dp) :: we(0:255), fe(0:255)
    logical  :: zig_exp_initialized = .false.

    public :: rvs_expon
    public :: pdf_expon
    public :: cdf_expon



    interface rvs_expon
    !! Version experimental
    !!
    !! Exponential Distribution Random Variates
    !! ([Specification](../page/specs/stdlib_stats_distribution_exponential.html#
    !! rvs_expon-exponential-distribution-random-variates))
    !!
        module procedure rvs_expon_0_rsp                 !0 dummy variable

        module procedure rvs_expon_rsp       !1 dummy variable
        module procedure rvs_expon_rdp       !1 dummy variable
        module procedure rvs_expon_csp       !1 dummy variable
        module procedure rvs_expon_cdp       !1 dummy variable

        module procedure rvs_expon_array_rsp !2 dummy variables
        module procedure rvs_expon_array_rdp !2 dummy variables
        module procedure rvs_expon_array_csp !2 dummy variables
        module procedure rvs_expon_array_cdp !2 dummy variables
    end interface rvs_expon



    interface pdf_expon
    !! Version experimental
    !!
    !! Exponential Distribution Probability Density Function
    !! ([Specification](../page/specs/stdlib_stats_distribution_exponential.html#
    !! pdf_expon-exponential-distribution-probability-density-function))
    !!
        module procedure pdf_expon_rsp
        module procedure pdf_expon_rdp
        module procedure pdf_expon_csp
        module procedure pdf_expon_cdp
    end interface pdf_expon



    interface cdf_expon
    !! Version experimental
    !!
    !! Exponential Distribution Cumulative Distribution Function
    !! ([Specification](../page/specs/stdlib_stats_distribution_exponential.html#
    !! cdf_expon-exponential-distribution-cumulative-distribution-function))
    !!
        module procedure cdf_expon_rsp
        module procedure cdf_expon_rdp
        module procedure cdf_expon_csp
        module procedure cdf_expon_cdp
    end interface cdf_expon





contains

    subroutine zigset
    ! Marsaglia & Tsang generator for random normals & random exponentials.
    ! Translated from C by Alan Miller (amiller@bigpond.net.au)
    !
    ! Marsaglia, G. & Tsang, W.W. (2000) `The ziggurat method for generating
    ! random variables', J. Statist. Software, v5(8).
    !
    ! This is an electronic journal which can be downloaded from:
    ! http://www.jstatsoft.org/v05/i08
    !
    ! Latest version - 1 January 2001
    !
        real(dp), parameter :: M2 = 2147483648.0_dp, ve = 0.003949659822581572_dp
        real(dp)            :: de, te, q
        integer :: i

        de = 7.697117470131487_dp
        te = de
    !tables for random exponetials
        q = ve * exp(de)
        ke(0) = int((de / q) * M2, kind = int32)
        ke(1) = 0
        we(0) = q / M2
        we(255) = de / M2
        fe(0) = ONE
        fe(255) = exp(- de)
        do  i = 254, 1, -1
            de = -log(ve / de + exp(- de))
            ke(i+1) = int(M2 * (de / te), kind = int32)
            te = de
            fe(i) = exp(- de)
            we(i) = de / M2
        end do
        zig_exp_initialized = .true.
    end subroutine zigset




    function rvs_expon_0_rsp( ) result(res)
    !
    ! Standard exponential random variate (lambda=1)
    !
        real(sp) :: res, x
        real(sp), parameter :: r = 7.69711747013104972_sp
        integer :: jz, iz

        if(.not. zig_exp_initialized ) call zigset
        iz = 0
        jz = dist_rand(1_int32)                !32bit random integer
        iz = iand( jz, 255 )                   !random integer in [0, 255]
        if( abs( jz ) < ke(iz) ) then
            res = abs(jz) * we(iz)
        else
            L1: do
                if( iz == 0 ) then
                    res = r - log( uni(1.0_sp) )
                    exit L1
                end if
                x = abs( jz ) * we(iz)
                if(fe(iz) + uni(1.0_sp) * (fe(iz-1) - fe(iz)) < exp(-x)) then
                    res = x
                    exit L1
                end if
                jz = dist_rand(1_int32)
                iz = iand( jz, 255 )
                if( abs( jz ) < ke(iz) ) then
                    res = abs( jz ) * we(iz)
                    exit L1
                end if
           end do L1
       endif
    end function rvs_expon_0_rsp

    function rvs_expon_0_rdp( ) result(res)
    !
    ! Standard exponential random variate (lambda=1)
    !
        real(dp) :: res, x
        real(dp), parameter :: r = 7.69711747013104972_dp
        integer :: jz, iz

        if(.not. zig_exp_initialized ) call zigset
        iz = 0
        jz = dist_rand(1_int32)                !32bit random integer
        iz = iand( jz, 255 )                   !random integer in [0, 255]
        if( abs( jz ) < ke(iz) ) then
            res = abs(jz) * we(iz)
        else
            L1: do
                if( iz == 0 ) then
                    res = r - log( uni(1.0_dp) )
                    exit L1
                end if
                x = abs( jz ) * we(iz)
                if(fe(iz) + uni(1.0_dp) * (fe(iz-1) - fe(iz)) < exp(-x)) then
                    res = x
                    exit L1
                end if
                jz = dist_rand(1_int32)
                iz = iand( jz, 255 )
                if( abs( jz ) < ke(iz) ) then
                    res = abs( jz ) * we(iz)
                    exit L1
                end if
           end do L1
       endif
    end function rvs_expon_0_rdp





    function rvs_expon_rsp(lambda) result(res)
    !
    ! Exponential distributed random variate
    !
        real(sp), intent(in) :: lambda
        real(sp) :: res


        if(lambda <= 0.0_sp) call error_stop("Error(rvs_expon): Exponen"   &
            //"tial distribution lambda parameter must be greater than zero")
        res = rvs_expon_0_rsp(  )
        res = res / lambda
    end function rvs_expon_rsp

    function rvs_expon_rdp(lambda) result(res)
    !
    ! Exponential distributed random variate
    !
        real(dp), intent(in) :: lambda
        real(dp) :: res


        if(lambda <= 0.0_dp) call error_stop("Error(rvs_expon): Exponen"   &
            //"tial distribution lambda parameter must be greater than zero")
        res = rvs_expon_0_rdp(  )
        res = res / lambda
    end function rvs_expon_rdp





    function rvs_expon_csp(lambda) result(res)
        complex(sp), intent(in) :: lambda
        complex(sp) :: res
        real(sp) :: tr, ti

        tr = rvs_expon_rsp(lambda % re)
        ti = rvs_expon_rsp(lambda % im)
        res = cmplx(tr, ti, kind=sp)
    end function rvs_expon_csp

    function rvs_expon_cdp(lambda) result(res)
        complex(dp), intent(in) :: lambda
        complex(dp) :: res
        real(dp) :: tr, ti

        tr = rvs_expon_rdp(lambda % re)
        ti = rvs_expon_rdp(lambda % im)
        res = cmplx(tr, ti, kind=dp)
    end function rvs_expon_cdp





    function rvs_expon_array_rsp(lambda, array_size) result(res)
        real(sp), intent(in) :: lambda
        integer, intent(in) :: array_size
        real(sp) :: res(array_size), x, re
        real(sp), parameter :: r = 7.69711747013104972_sp
        integer :: jz, iz, i

        if(lambda <= 0.0_sp) call error_stop("Error(rvs_expon_array): Exp" &
            //"oonential distribution lambda parameter must be greater than zero")

        if(.not. zig_exp_initialized) call zigset
        do i = 1, array_size
            iz = 0
            jz = dist_rand(1_int32)
            iz = iand( jz, 255 )
            if( abs( jz ) < ke(iz) ) then
                re = abs(jz) * we(iz)
            else
                L1: do
                    if( iz == 0 ) then
                        re = r - log( uni(1.0_sp) )
                        exit L1
                    end if
                    x = abs( jz ) * we(iz)
                    if(fe(iz) + uni(1.0_sp)*(fe(iz-1)-fe(iz)) < exp(-x)) then
                        re = x
                        exit L1
                    end if
                    jz = dist_rand(1_int32)
                    iz = iand( jz, 255 )
                    if( abs( jz ) < ke(iz) ) then
                        re = abs( jz ) * we(iz)
                        exit L1
                    end if
               end do L1
            endif
            res(i) = re / lambda
        end do
    end function rvs_expon_array_rsp

    function rvs_expon_array_rdp(lambda, array_size) result(res)
        real(dp), intent(in) :: lambda
        integer, intent(in) :: array_size
        real(dp) :: res(array_size), x, re
        real(dp), parameter :: r = 7.69711747013104972_dp
        integer :: jz, iz, i

        if(lambda <= 0.0_dp) call error_stop("Error(rvs_expon_array): Exp" &
            //"oonential distribution lambda parameter must be greater than zero")

        if(.not. zig_exp_initialized) call zigset
        do i = 1, array_size
            iz = 0
            jz = dist_rand(1_int32)
            iz = iand( jz, 255 )
            if( abs( jz ) < ke(iz) ) then
                re = abs(jz) * we(iz)
            else
                L1: do
                    if( iz == 0 ) then
                        re = r - log( uni(1.0_dp) )
                        exit L1
                    end if
                    x = abs( jz ) * we(iz)
                    if(fe(iz) + uni(1.0_dp)*(fe(iz-1)-fe(iz)) < exp(-x)) then
                        re = x
                        exit L1
                    end if
                    jz = dist_rand(1_int32)
                    iz = iand( jz, 255 )
                    if( abs( jz ) < ke(iz) ) then
                        re = abs( jz ) * we(iz)
                        exit L1
                    end if
               end do L1
            endif
            res(i) = re / lambda
        end do
    end function rvs_expon_array_rdp





    function rvs_expon_array_csp(lambda, array_size) result(res)
        complex(sp), intent(in) :: lambda
        integer, intent(in) :: array_size
        complex(sp) :: res(array_size)
        integer :: i
        real(sp) :: tr, ti

        do i = 1, array_size
            tr = rvs_expon_rsp(lambda % re)
            ti = rvs_expon_rsp(lambda % im)
            res(i) = cmplx(tr, ti, kind=sp)
        end do
    end function rvs_expon_array_csp

    function rvs_expon_array_cdp(lambda, array_size) result(res)
        complex(dp), intent(in) :: lambda
        integer, intent(in) :: array_size
        complex(dp) :: res(array_size)
        integer :: i
        real(dp) :: tr, ti

        do i = 1, array_size
            tr = rvs_expon_rdp(lambda % re)
            ti = rvs_expon_rdp(lambda % im)
            res(i) = cmplx(tr, ti, kind=dp)
        end do
    end function rvs_expon_array_cdp





    impure elemental function pdf_expon_rsp(x, lambda) result(res)
    !
    ! Exponential Distribution Probability Density Function
    !
        real(sp), intent(in) :: x, lambda
        real(sp) :: res

        if(lambda <= 0.0_sp) call error_stop("Error(pdf_expon): Expon"     &
            //"ential distribution lambda parameter must be greater than zero")
        if(x < 0.0_sp) call error_stop("Error(pdf_expon): Exponential"     &
            //" distribution variate x must be non-negative")
        res = exp(- x * lambda) * lambda
    end function pdf_expon_rsp

    impure elemental function pdf_expon_rdp(x, lambda) result(res)
    !
    ! Exponential Distribution Probability Density Function
    !
        real(dp), intent(in) :: x, lambda
        real(dp) :: res

        if(lambda <= 0.0_dp) call error_stop("Error(pdf_expon): Expon"     &
            //"ential distribution lambda parameter must be greater than zero")
        if(x < 0.0_dp) call error_stop("Error(pdf_expon): Exponential"     &
            //" distribution variate x must be non-negative")
        res = exp(- x * lambda) * lambda
    end function pdf_expon_rdp





    impure elemental function pdf_expon_csp(x, lambda) result(res)
        complex(sp), intent(in) :: x, lambda
        real(sp) :: res

        res = pdf_expon_rsp(x % re, lambda % re)
        res = res * pdf_expon_rsp(x % im, lambda % im)
    end function pdf_expon_csp

    impure elemental function pdf_expon_cdp(x, lambda) result(res)
        complex(dp), intent(in) :: x, lambda
        real(dp) :: res

        res = pdf_expon_rdp(x % re, lambda % re)
        res = res * pdf_expon_rdp(x % im, lambda % im)
    end function pdf_expon_cdp





    impure elemental function cdf_expon_rsp(x, lambda) result(res)
    !
    ! Exponential Distribution Cumulative Distribution Function
    !
        real(sp), intent(in) :: x, lambda
        real(sp) :: res

        if(lambda <= 0.0_sp) call error_stop("Error(cdf_expon): Expon"     &
            //"ential distribution lambda parameter must be greater than zero")
        if(x < 0.0_sp) call error_stop("Error(cdf_expon): Exponential"     &
            //" distribution variate x must be non-negative")
        res = 1.0_sp - exp(- x * lambda)
    end function cdf_expon_rsp

    impure elemental function cdf_expon_rdp(x, lambda) result(res)
    !
    ! Exponential Distribution Cumulative Distribution Function
    !
        real(dp), intent(in) :: x, lambda
        real(dp) :: res

        if(lambda <= 0.0_dp) call error_stop("Error(cdf_expon): Expon"     &
            //"ential distribution lambda parameter must be greater than zero")
        if(x < 0.0_dp) call error_stop("Error(cdf_expon): Exponential"     &
            //" distribution variate x must be non-negative")
        res = 1.0_dp - exp(- x * lambda)
    end function cdf_expon_rdp





    impure elemental function cdf_expon_csp(x, lambda) result(res)
        complex(sp), intent(in) :: x, lambda
        real(sp) :: res

        res = cdf_expon_rsp(x % re, lambda % re)
        res = res * cdf_expon_rsp(x % im, lambda % im)
    end function cdf_expon_csp

    impure elemental function cdf_expon_cdp(x, lambda) result(res)
        complex(dp), intent(in) :: x, lambda
        real(dp) :: res

        res = cdf_expon_rdp(x % re, lambda % re)
        res = res * cdf_expon_rdp(x % im, lambda % im)
    end function cdf_expon_cdp


end module stdlib_stats_distribution_exponential
