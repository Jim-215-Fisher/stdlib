program test_distribution_expon
    use stdlib_kinds
    use stdlib_error, only : check
    use stdlib_stats_distribution_PRNG, only : random_seed
    use stdlib_stats_distribution_expon, only :                                 &
        expon_rvs => exponential_distribution_rvs,                              &
        expon_pdf => exponential_distribution_pdf,                              &
        expon_cdf => exponential_distribution_cdf

    implicit none
    real(sp), parameter :: sptol = 1000 * epsilon(1.0_sp)
    real(dp), parameter :: dptol = 1000 * epsilon(1.0_dp)
    real(qp), parameter :: qptol = 1000 * epsilon(1.0_qp)
    logical ::  warn = .true.
    integer :: put, get

    put = 12345678
    call random_seed(put, get)

    call test_exponential_random_generator

      call test_expon_rvs_rsp
      call test_expon_rvs_rdp
      call test_expon_rvs_rqp
      call test_expon_rvs_csp
      call test_expon_rvs_cdp
      call test_expon_rvs_cqp

    call test_expon_pdf_rsp
    call test_expon_pdf_rdp
    call test_expon_pdf_rqp
    call test_expon_pdf_csp
    call test_expon_pdf_cdp
    call test_expon_pdf_cqp

    call test_expon_cdf_rsp
    call test_expon_cdf_rdp
    call test_expon_cdf_rqp
    call test_expon_cdf_csp
    call test_expon_cdf_cdp
    call test_expon_cdf_cqp


    contains

    subroutine test_exponential_random_generator
        integer :: i, j, freq(0:1000), num=10000000
        real(dp) :: chisq, expct

        print *, ""
        print *, "Test exponential random generator with chi-squared"
        freq = 0
        do i = 1, num
            j = 1000 * (1.0 - exp(- expon_rvs(1.0)))
            freq(j) = freq(j) + 1
        end do
        chisq = 0.0_dp
        expct = num / 1000
        do i = 0, 999
           chisq = chisq + (freq(i) - expct) ** 2 / expct
        end do
        write(*,*) "The critical values for chi-squared with 1000 d. of f. is"  &
                   //" 1143.92"
        write(*,*) "Chi-squared for exponential random generator is : ", chisq
        call check((chisq < 1143.9), &
               msg="exponential randomness failed chi-squared test", warn=warn)
    end subroutine test_exponential_random_generator

    subroutine test_expon_rvs_rsp
        real(sp) :: res(10), scale
        integer :: i, n, k = 5
        integer :: seed, get
        real(sp) :: ans(10) = [1.37178108290154243675829093263018876_sp,      &
                            0.309467303067628429769087006206973456_sp,      &
                            0.303573393522570872393728791394096334_sp,      &
                             3.00973636550766943109636031294940040_sp,      &
                            0.250084097046766984373533659891108982_sp,      &
                             1.20139122141795795517538181229610927_sp,      &
                             4.43019214257137261547825346497120336_sp,      &
                            0.835001950484080046610557701569632627_sp,      &
                             1.82681711031524329769126779865473509_sp,      &
                            0.910435173630070204708886194566730410_sp]

        print *, "Test exponential_distribution_rvs_rsp"
        seed = 593742186
        call random_seed(seed, get)
        scale = 1.5_sp
        do i = 1, 5
            res(i) = expon_rvs(scale)
        end do
        res(6:10) = expon_rvs(scale, k)
        call check(all(abs(res - ans) < sptol), &
            msg="exponential_distribution_rvs_rsp failed", warn=warn)
    end subroutine test_expon_rvs_rsp

    subroutine test_expon_rvs_rdp
        real(dp) :: res(10), scale
        integer :: i, n, k = 5
        integer :: seed, get
        real(dp) :: ans(10) = [1.37178108290154243675829093263018876_dp,      &
                            0.309467303067628429769087006206973456_dp,      &
                            0.303573393522570872393728791394096334_dp,      &
                             3.00973636550766943109636031294940040_dp,      &
                            0.250084097046766984373533659891108982_dp,      &
                             1.20139122141795795517538181229610927_dp,      &
                             4.43019214257137261547825346497120336_dp,      &
                            0.835001950484080046610557701569632627_dp,      &
                             1.82681711031524329769126779865473509_dp,      &
                            0.910435173630070204708886194566730410_dp]

        print *, "Test exponential_distribution_rvs_rdp"
        seed = 593742186
        call random_seed(seed, get)
        scale = 1.5_dp
        do i = 1, 5
            res(i) = expon_rvs(scale)
        end do
        res(6:10) = expon_rvs(scale, k)
        call check(all(abs(res - ans) < dptol), &
            msg="exponential_distribution_rvs_rdp failed", warn=warn)
    end subroutine test_expon_rvs_rdp

    subroutine test_expon_rvs_rqp
        real(qp) :: res(10), scale
        integer :: i, n, k = 5
        integer :: seed, get
        real(qp) :: ans(10) = [1.37178108290154243675829093263018876_qp,      &
                            0.309467303067628429769087006206973456_qp,      &
                            0.303573393522570872393728791394096334_qp,      &
                             3.00973636550766943109636031294940040_qp,      &
                            0.250084097046766984373533659891108982_qp,      &
                             1.20139122141795795517538181229610927_qp,      &
                             4.43019214257137261547825346497120336_qp,      &
                            0.835001950484080046610557701569632627_qp,      &
                             1.82681711031524329769126779865473509_qp,      &
                            0.910435173630070204708886194566730410_qp]

        print *, "Test exponential_distribution_rvs_rqp"
        seed = 593742186
        call random_seed(seed, get)
        scale = 1.5_qp
        do i = 1, 5
            res(i) = expon_rvs(scale)
        end do
        res(6:10) = expon_rvs(scale, k)
        call check(all(abs(res - ans) < qptol), &
            msg="exponential_distribution_rvs_rqp failed", warn=warn)
    end subroutine test_expon_rvs_rqp

    subroutine test_expon_rvs_csp
        complex(sp) :: res(10), scale
        integer :: i, n, k = 5
        integer :: seed, get
        complex(sp) :: ans(10) = [(0.640164494514465332031250000000000000_sp,    &
                             0.268204987049102783203125000000000000_sp),    &
                             (0.141667589545249938964843750000000000_sp,    &
                              2.60843825340270996093750000000000000_sp),    &
                             (0.116705909371376037597656250000000000_sp,    &
                              1.04120576381683349609375000000000000_sp),    &
                              (2.06742310523986816406250000000000000_sp,    &
                             0.723668336868286132812500000000000000_sp),    &
                             (0.852514624595642089843750000000000000_sp,    &
                             0.789043843746185302734375000000000000_sp),    &
                              (1.09098446369171142578125000000000000_sp,    &
                              1.48569476604461669921875000000000000_sp),    &
                              (4.29633092880249023437500000000000000_sp,    &
                             0.338216394186019897460937500000000000_sp),    &
                             (0.340462744235992431640625000000000000_sp,    &
                             0.172319442033767700195312500000000000_sp),    &
                       (6.932352483272552490234375000000000000E-0002_sp,    &
                       6.742518395185470581054687500000000000E-0002_sp),    &
                              (1.03231632709503173828125000000000000_sp,    &
                             0.421413004398345947265625000000000000_sp)]

        print *, "Test exponential_distribution_rvs_csp"
        seed = 593742186
        call random_seed(seed, get)
        scale = (0.7_sp, 1.3_sp)
        do i = 1, 5
            res(i) = expon_rvs(scale)
        end do
        res(6:10) = expon_rvs(scale, k)
        call check(all(abs(res - ans) < sptol), &
            msg="exponential_distribution_rvs_csp failed", warn=warn)
    end subroutine test_expon_rvs_csp

    subroutine test_expon_rvs_cdp
        complex(dp) :: res(10), scale
        integer :: i, n, k = 5
        integer :: seed, get
        complex(dp) :: ans(10) = [(0.640164494514465332031250000000000000_dp,    &
                             0.268204987049102783203125000000000000_dp),    &
                             (0.141667589545249938964843750000000000_dp,    &
                              2.60843825340270996093750000000000000_dp),    &
                             (0.116705909371376037597656250000000000_dp,    &
                              1.04120576381683349609375000000000000_dp),    &
                              (2.06742310523986816406250000000000000_dp,    &
                             0.723668336868286132812500000000000000_dp),    &
                             (0.852514624595642089843750000000000000_dp,    &
                             0.789043843746185302734375000000000000_dp),    &
                              (1.09098446369171142578125000000000000_dp,    &
                              1.48569476604461669921875000000000000_dp),    &
                              (4.29633092880249023437500000000000000_dp,    &
                             0.338216394186019897460937500000000000_dp),    &
                             (0.340462744235992431640625000000000000_dp,    &
                             0.172319442033767700195312500000000000_dp),    &
                       (6.932352483272552490234375000000000000E-0002_dp,    &
                       6.742518395185470581054687500000000000E-0002_dp),    &
                              (1.03231632709503173828125000000000000_dp,    &
                             0.421413004398345947265625000000000000_dp)]

        print *, "Test exponential_distribution_rvs_cdp"
        seed = 593742186
        call random_seed(seed, get)
        scale = (0.7_dp, 1.3_dp)
        do i = 1, 5
            res(i) = expon_rvs(scale)
        end do
        res(6:10) = expon_rvs(scale, k)
        call check(all(abs(res - ans) < dptol), &
            msg="exponential_distribution_rvs_cdp failed", warn=warn)
    end subroutine test_expon_rvs_cdp

    subroutine test_expon_rvs_cqp
        complex(qp) :: res(10), scale
        integer :: i, n, k = 5
        integer :: seed, get
        complex(qp) :: ans(10) = [(0.640164494514465332031250000000000000_qp,    &
                             0.268204987049102783203125000000000000_qp),    &
                             (0.141667589545249938964843750000000000_qp,    &
                              2.60843825340270996093750000000000000_qp),    &
                             (0.116705909371376037597656250000000000_qp,    &
                              1.04120576381683349609375000000000000_qp),    &
                              (2.06742310523986816406250000000000000_qp,    &
                             0.723668336868286132812500000000000000_qp),    &
                             (0.852514624595642089843750000000000000_qp,    &
                             0.789043843746185302734375000000000000_qp),    &
                              (1.09098446369171142578125000000000000_qp,    &
                              1.48569476604461669921875000000000000_qp),    &
                              (4.29633092880249023437500000000000000_qp,    &
                             0.338216394186019897460937500000000000_qp),    &
                             (0.340462744235992431640625000000000000_qp,    &
                             0.172319442033767700195312500000000000_qp),    &
                       (6.932352483272552490234375000000000000E-0002_qp,    &
                       6.742518395185470581054687500000000000E-0002_qp),    &
                              (1.03231632709503173828125000000000000_qp,    &
                             0.421413004398345947265625000000000000_qp)]

        print *, "Test exponential_distribution_rvs_cqp"
        seed = 593742186
        call random_seed(seed, get)
        scale = (0.7_qp, 1.3_qp)
        do i = 1, 5
            res(i) = expon_rvs(scale)
        end do
        res(6:10) = expon_rvs(scale, k)
        call check(all(abs(res - ans) < qptol), &
            msg="exponential_distribution_rvs_cqp failed", warn=warn)
    end subroutine test_expon_rvs_cqp



    subroutine test_expon_pdf_rsp
        real(sp) :: x1, x2(3,4), scale
        integer :: i, n
        integer :: seed, get
        real :: res(3,5)
        real :: ans(15) = [6.14960417E-02, 6.14960417E-02, 6.14960417E-02,      &
                           1.38718796, 0.729365528, 1.67107172E-02, 0.296734482,&
                           2.36971565E-02, 4.06475700E-02, 1.18497872,          &
                           8.39852914E-02, 1.36920142, 1.54058458E-02,          &
                           3.20194475E-02, 0.603879571]

        print *, "Test exponential_distribution_pdf_rsp"
        seed = 123987654
        call random_seed(seed, get)
        scale = 1.5_sp
        x1 = expon_rvs(scale)
        x2 = reshape(expon_rvs(scale, 12), [3,4])
        res(:,1) = expon_pdf(x1, scale)
        res(:, 2:5) = expon_pdf(x2, scale)
        call check(all(abs(res - reshape(ans, [3,5])) < sptol),             &
            msg="exponential_distribution_pdf_rsp failed", warn=warn)
    end subroutine test_expon_pdf_rsp

    subroutine test_expon_pdf_rdp
        real(dp) :: x1, x2(3,4), scale
        integer :: i, n
        integer :: seed, get
        real :: res(3,5)
        real :: ans(15) = [6.14960417E-02, 6.14960417E-02, 6.14960417E-02,      &
                           1.38718796, 0.729365528, 1.67107172E-02, 0.296734482,&
                           2.36971565E-02, 4.06475700E-02, 1.18497872,          &
                           8.39852914E-02, 1.36920142, 1.54058458E-02,          &
                           3.20194475E-02, 0.603879571]

        print *, "Test exponential_distribution_pdf_rdp"
        seed = 123987654
        call random_seed(seed, get)
        scale = 1.5_dp
        x1 = expon_rvs(scale)
        x2 = reshape(expon_rvs(scale, 12), [3,4])
        res(:,1) = expon_pdf(x1, scale)
        res(:, 2:5) = expon_pdf(x2, scale)
        call check(all(abs(res - reshape(ans, [3,5])) < dptol),             &
            msg="exponential_distribution_pdf_rdp failed", warn=warn)
    end subroutine test_expon_pdf_rdp

    subroutine test_expon_pdf_rqp
        real(qp) :: x1, x2(3,4), scale
        integer :: i, n
        integer :: seed, get
        real :: res(3,5)
        real :: ans(15) = [6.14960417E-02, 6.14960417E-02, 6.14960417E-02,      &
                           1.38718796, 0.729365528, 1.67107172E-02, 0.296734482,&
                           2.36971565E-02, 4.06475700E-02, 1.18497872,          &
                           8.39852914E-02, 1.36920142, 1.54058458E-02,          &
                           3.20194475E-02, 0.603879571]

        print *, "Test exponential_distribution_pdf_rqp"
        seed = 123987654
        call random_seed(seed, get)
        scale = 1.5_qp
        x1 = expon_rvs(scale)
        x2 = reshape(expon_rvs(scale, 12), [3,4])
        res(:,1) = expon_pdf(x1, scale)
        res(:, 2:5) = expon_pdf(x2, scale)
        call check(all(abs(res - reshape(ans, [3,5])) < qptol),             &
            msg="exponential_distribution_pdf_rqp failed", warn=warn)
    end subroutine test_expon_pdf_rqp

    subroutine test_expon_pdf_csp
        complex(sp) :: x1, x2(3,4), scale
        integer :: i, n
        integer :: seed, get
        real :: res(3,5)
        real :: ans(15) = [0.386471182, 0.386471182,0.386471182, 2.79592816E-03,&
                           4.01333207E-03, 0.317740440, 0.385551631,            &
                           5.02163824E-03, 0.372386932, 6.09764457E-03,         &
                           0.273956627, 0.407586545, 1.59074634E-03,0.136133015,&
                           0.399842113]

        print *, "Test exponential_distribution_pdf_csp"
        seed = 123987654
        call random_seed(seed, get)
        scale = (0.3_sp, 1.6_sp)
        x1 = expon_rvs(scale)
        x2 = reshape(expon_rvs(scale, 12), [3,4])
        res(:,1) = expon_pdf(x1, scale)
        res(:, 2:5) = expon_pdf(x2, scale)
        call check(all(abs(res - reshape(ans, [3,5])) < sptol),             &
            msg="exponential_distribution_pdf_csp failed", warn=warn)
    end subroutine test_expon_pdf_csp

    subroutine test_expon_pdf_cdp
        complex(dp) :: x1, x2(3,4), scale
        integer :: i, n
        integer :: seed, get
        real :: res(3,5)
        real :: ans(15) = [0.386471182, 0.386471182,0.386471182, 2.79592816E-03,&
                           4.01333207E-03, 0.317740440, 0.385551631,            &
                           5.02163824E-03, 0.372386932, 6.09764457E-03,         &
                           0.273956627, 0.407586545, 1.59074634E-03,0.136133015,&
                           0.399842113]

        print *, "Test exponential_distribution_pdf_cdp"
        seed = 123987654
        call random_seed(seed, get)
        scale = (0.3_dp, 1.6_dp)
        x1 = expon_rvs(scale)
        x2 = reshape(expon_rvs(scale, 12), [3,4])
        res(:,1) = expon_pdf(x1, scale)
        res(:, 2:5) = expon_pdf(x2, scale)
        call check(all(abs(res - reshape(ans, [3,5])) < dptol),             &
            msg="exponential_distribution_pdf_cdp failed", warn=warn)
    end subroutine test_expon_pdf_cdp

    subroutine test_expon_pdf_cqp
        complex(qp) :: x1, x2(3,4), scale
        integer :: i, n
        integer :: seed, get
        real :: res(3,5)
        real :: ans(15) = [0.386471182, 0.386471182,0.386471182, 2.79592816E-03,&
                           4.01333207E-03, 0.317740440, 0.385551631,            &
                           5.02163824E-03, 0.372386932, 6.09764457E-03,         &
                           0.273956627, 0.407586545, 1.59074634E-03,0.136133015,&
                           0.399842113]

        print *, "Test exponential_distribution_pdf_cqp"
        seed = 123987654
        call random_seed(seed, get)
        scale = (0.3_qp, 1.6_qp)
        x1 = expon_rvs(scale)
        x2 = reshape(expon_rvs(scale, 12), [3,4])
        res(:,1) = expon_pdf(x1, scale)
        res(:, 2:5) = expon_pdf(x2, scale)
        call check(all(abs(res - reshape(ans, [3,5])) < qptol),             &
            msg="exponential_distribution_pdf_cqp failed", warn=warn)
    end subroutine test_expon_pdf_cqp


    subroutine test_expon_cdf_rsp
        real(sp) :: x1, x2(3,4), scale
        integer :: i, n
        integer :: seed, get
        real :: res(3,5)
        real :: ans(15) = [0.370481908, 0.370481908, 0.370481908, 0.993631542,  &
                           0.246292457, 0.497423291, 0.999946773, 0.995101511,  &
                           0.305115670, 0.932980001, 0.967154086, 0.777015686,  &
                           0.983209372, 8.37164521E-02, 0.275721848]

        print *, "Test exponential_distribution_cdf_rsp"
        seed = 621957438
        call random_seed(seed, get)

        scale = 2.0_sp
        x1 = expon_rvs(scale)
        x2 = reshape(expon_rvs(scale, 12), [3,4])
        res(:,1) = expon_cdf(x1, scale)
        res(:, 2:5) = expon_cdf(x2, scale)
        call check(all(abs(res - reshape(ans,[3,5])) < sptol),             &
            msg="exponential_distribution_cdf_rsp failed", warn=warn)
    end subroutine test_expon_cdf_rsp

    subroutine test_expon_cdf_rdp
        real(dp) :: x1, x2(3,4), scale
        integer :: i, n
        integer :: seed, get
        real :: res(3,5)
        real :: ans(15) = [0.370481908, 0.370481908, 0.370481908, 0.993631542,  &
                           0.246292457, 0.497423291, 0.999946773, 0.995101511,  &
                           0.305115670, 0.932980001, 0.967154086, 0.777015686,  &
                           0.983209372, 8.37164521E-02, 0.275721848]

        print *, "Test exponential_distribution_cdf_rdp"
        seed = 621957438
        call random_seed(seed, get)

        scale = 2.0_dp
        x1 = expon_rvs(scale)
        x2 = reshape(expon_rvs(scale, 12), [3,4])
        res(:,1) = expon_cdf(x1, scale)
        res(:, 2:5) = expon_cdf(x2, scale)
        call check(all(abs(res - reshape(ans,[3,5])) < dptol),             &
            msg="exponential_distribution_cdf_rdp failed", warn=warn)
    end subroutine test_expon_cdf_rdp

    subroutine test_expon_cdf_rqp
        real(qp) :: x1, x2(3,4), scale
        integer :: i, n
        integer :: seed, get
        real :: res(3,5)
        real :: ans(15) = [0.370481908, 0.370481908, 0.370481908, 0.993631542,  &
                           0.246292457, 0.497423291, 0.999946773, 0.995101511,  &
                           0.305115670, 0.932980001, 0.967154086, 0.777015686,  &
                           0.983209372, 8.37164521E-02, 0.275721848]

        print *, "Test exponential_distribution_cdf_rqp"
        seed = 621957438
        call random_seed(seed, get)

        scale = 2.0_qp
        x1 = expon_rvs(scale)
        x2 = reshape(expon_rvs(scale, 12), [3,4])
        res(:,1) = expon_cdf(x1, scale)
        res(:, 2:5) = expon_cdf(x2, scale)
        call check(all(abs(res - reshape(ans,[3,5])) < qptol),             &
            msg="exponential_distribution_cdf_rqp failed", warn=warn)
    end subroutine test_expon_cdf_rqp

    subroutine test_expon_cdf_csp
        complex(sp) :: x1, x2(3,4), scale
        integer :: i, n
        integer :: seed, get
        real :: res(3,5)
        real :: ans(15) = [0.176930442,0.176930442,0.176930442, 5.98644949E-02, &
                           0.981560826,0.135309443,0.617795825, 7.55468607E-02, &
                           0.118341736, 0.484595388, 0.794088185, 0.912919402,  &
                           0.914170802, 0.370377690, 0.793968201]

        print *, "Test exponential_distribution_cdf_csp"
        seed = 621957438
        call random_seed(seed, get)

        scale = (1.3_sp, 2.1_sp)
        x1 = expon_rvs(scale)
        x2 = reshape(expon_rvs(scale, 12), [3,4])
        res(:,1) = expon_cdf(x1, scale)
        res(:, 2:5) = expon_cdf(x2, scale)
        call check(all(abs(res - reshape(ans,[3,5])) < sptol),             &
            msg="exponential_distribution_cdf_csp failed", warn=warn)
    end subroutine test_expon_cdf_csp

    subroutine test_expon_cdf_cdp
        complex(dp) :: x1, x2(3,4), scale
        integer :: i, n
        integer :: seed, get
        real :: res(3,5)
        real :: ans(15) = [0.176930442,0.176930442,0.176930442, 5.98644949E-02, &
                           0.981560826,0.135309443,0.617795825, 7.55468607E-02, &
                           0.118341736, 0.484595388, 0.794088185, 0.912919402,  &
                           0.914170802, 0.370377690, 0.793968201]

        print *, "Test exponential_distribution_cdf_cdp"
        seed = 621957438
        call random_seed(seed, get)

        scale = (1.3_dp, 2.1_dp)
        x1 = expon_rvs(scale)
        x2 = reshape(expon_rvs(scale, 12), [3,4])
        res(:,1) = expon_cdf(x1, scale)
        res(:, 2:5) = expon_cdf(x2, scale)
        call check(all(abs(res - reshape(ans,[3,5])) < dptol),             &
            msg="exponential_distribution_cdf_cdp failed", warn=warn)
    end subroutine test_expon_cdf_cdp

    subroutine test_expon_cdf_cqp
        complex(qp) :: x1, x2(3,4), scale
        integer :: i, n
        integer :: seed, get
        real :: res(3,5)
        real :: ans(15) = [0.176930442,0.176930442,0.176930442, 5.98644949E-02, &
                           0.981560826,0.135309443,0.617795825, 7.55468607E-02, &
                           0.118341736, 0.484595388, 0.794088185, 0.912919402,  &
                           0.914170802, 0.370377690, 0.793968201]

        print *, "Test exponential_distribution_cdf_cqp"
        seed = 621957438
        call random_seed(seed, get)

        scale = (1.3_qp, 2.1_qp)
        x1 = expon_rvs(scale)
        x2 = reshape(expon_rvs(scale, 12), [3,4])
        res(:,1) = expon_cdf(x1, scale)
        res(:, 2:5) = expon_cdf(x2, scale)
        call check(all(abs(res - reshape(ans,[3,5])) < qptol),             &
            msg="exponential_distribution_cdf_cqp failed", warn=warn)
    end subroutine test_expon_cdf_cqp

end program test_distribution_expon