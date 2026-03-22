module Sparity
    contains
    subroutine Amplitudes_sp(G, P, q, gamval, rho, Amp, Jval, Amp_p, Jval_p, nzval, &
                              rho_min_user, rho_max_user)
        use stdlib_sorting, only: sort_index
        implicit none
    
        integer, parameter :: dp = selected_real_kind(15, 307)
    
        ! Arguments
        complex(dp), dimension(:,:), intent(in) :: P
        complex(dp), dimension(:), intent(in) :: q
        real(dp), dimension(:,:), intent(in) :: G
        real(dp), dimension(:), intent(in) :: gamval
        real(dp), intent(in) :: rho
        complex(dp), dimension(:,:), intent(inout) :: Amp, Amp_p
        real(dp), dimension(:), intent(inout) :: Jval, Jval_p
        integer, dimension(:), intent(inout) :: nzval
        real(dp), intent(in) :: rho_min_user, rho_max_user

        ! Internal variables
        complex(dp), allocatable :: Iden(:,:), Prho(:,:), Plow(:,:), y(:), z(:), u(:), B(:), xnew(:), v(:), &
        znew(:), GT_G(:,:), E(:,:), KKT(:,:), right_hand_side(:), xpol(:), sol(:)
        integer, allocatable :: ipiv(:)
        integer :: n, m, i, j, k, jj, info, count, mm
        real(dp) :: a, res_dual, res_prim, eps_prim, eps_dual, trace_GG
        real(dp), parameter :: eps_abs = 1.0e-6_dp, eps_rel = 1.e-4_dp 
        integer, parameter :: Max_ADMM_Iter = 10000000
        integer, allocatable :: index_s(:)
        real, allocatable :: abz(:)

        ! Adaptive rho variables (Boyd et al. 2011, Section 3.4.1)
        real(dp) :: rho_k
        real(dp), parameter :: mu_adapt  = 10.0_dp
        real(dp), parameter :: tau_adapt = 2.0_dp
        logical :: need_refactor
        real(dp) :: res_prim_prev
        integer :: stagnation_count
        integer, parameter :: stagnation_limit = 50000
        integer :: cooldown_counter
        integer, parameter :: cooldown_iters = 100
        integer, parameter :: print_every = 5000    ! Print residuals every N iterations
    
        ! Size variables
        n = size(G, 1)
        m = size(gamval, 1)
    
        allocate(Iden(n,n), Prho(n,n), Plow(n,n))
        allocate(y(n), z(n), u(n), B(n), xnew(n), v(n), znew(n), xpol(n), sol(n))
        allocate(GT_G(n,n))
        allocate(abz(n), index_s(n))

        ! Identity matrix
        Iden = 0.0_dp
        do i = 1, n
           Iden(i,i) = 1.0_dp
        end do
    
        ! Trace of G'*G
        GT_G = matmul(transpose(G), G)
        trace_GG = 0.0_dp
        do i = 1, n
            trace_GG = trace_GG + GT_G(i,i)
        end do
        deallocate(GT_G)

        ! Initialize ADMM variables once — warm-started across gamma values
        y(:) = 0.0_dp
        z(:) = 0.0_dp
        rho_k = rho
        need_refactor = .true.
        cooldown_counter = cooldown_iters  ! allow immediate rho change at start

        open(unit=2000, file="residuals.dat", status="replace")

        ! =====================================================================
        ! Loop over gamma values — warm-starting from previous solution
        ! =====================================================================
        do k = 1, m

            stagnation_count = 0
            res_prim_prev = huge(1.0_dp)
            need_refactor = .true.

            write(*, '(A,I4,A,ES12.4)') ' --- ADMM for Gamma index ', k, ',  gamma = ', gamval(k)

            ! ADMM iterations for this gamma
            do jj = 1, Max_ADMM_Iter

                if (need_refactor) then
                    rho_k = max(rho_min_user, min(rho_max_user, rho_k))
                    Prho = P + ((rho_k / 2.0_dp) * Iden)
                    Plow = Prho
                    call zpotrf('L', n, Plow, n, info)
                    if (info /= 0) then
                        rho_k = rho_k * 10.0_dp
                        rho_k = max(rho_min_user, min(rho_max_user, rho_k))
                        Prho = P + ((rho_k / 2.0_dp) * Iden)
                        Plow = Prho
                        call zpotrf('L', n, Plow, n, info)
                    end if
                    need_refactor = .false.
                end if

                ! x-update
                u = z - (1.0_dp / rho_k) * y
                B = q + (rho_k / 2.0_dp) * u
                call zpotrs('L', n, 1, Plow, n, B, n, info)
                xnew = B

                ! z-update (soft thresholding)
                a = gamval(k) / rho_k
                v = xnew + (1.0_dp / rho_k) * y
                do i = 1, n
                    if (abs(v(i)) > a) then
                        znew(i) = (1.0_dp - a / abs(v(i))) * v(i)
                    else
                        znew(i) = 0.0_dp
                    end if
                end do

                ! Residuals
                res_prim = sqrt(sum(abs(xnew - znew)**2))
                res_dual = rho_k * sqrt(sum(abs(znew - z)**2))

                ! Dual variable update
                y = y + rho_k * (xnew - znew)

                ! Convergence tolerances
                eps_prim = sqrt(n * 1.0_dp) * eps_abs &
                         + eps_rel * max(sqrt(sum(abs(xnew)**2)), sqrt(sum(abs(znew)**2)))
                eps_dual = sqrt(n * 1.0_dp) * eps_abs &
                         + eps_rel * sqrt(sum(abs(y)**2))

                write(2000, *) jj, res_prim, eps_prim, res_dual, eps_dual, rho_k

                ! Print progress every print_every iterations
                if (mod(jj, print_every) == 0) then
                    count = 0
                    do i = 1, n
                        if (abs(znew(i)) > 1.0e-6_dp) count = count + 1
                    end do
                    write(*, '(A,I8,A,ES9.2,A,ES9.2,A,ES9.2,A,ES9.2,A,ES9.2,A,I4)') &
                        '  it=', jj, &
                        ' rp=', res_prim, ' ep=', eps_prim, &
                        ' rd=', res_dual, ' ed=', eps_dual, &
                        ' rho=', rho_k, ' nnz=', count
                end if

                ! Convergence check
                if (res_prim < eps_prim .and. res_dual < eps_dual) then
                    count = 0
                    do i = 1, n
                        if (abs(znew(i)) > 1.0e-6_dp) count = count + 1
                    end do
                    write(*, '(A,I8,A,ES9.2,A,ES9.2,A,ES9.2,A,I4)') &
                        '  CONVERGED iter=', jj, &
                        ' rp=', res_prim, ' rd=', res_dual, &
                        ' rho=', rho_k, ' nnz=', count
                    z = znew
                    exit
                end if

                ! Adaptive rho with cooldown (Boyd et al. 2011, Section 3.4.1)
                ! Guard: only change rho if the new value stays within bounds
                ! Otherwise y gets scaled but rho is clamped → corrupts u = y/rho
                cooldown_counter = cooldown_counter + 1
                if (cooldown_counter >= cooldown_iters) then
                    if (res_prim > mu_adapt * res_dual .and. rho_k * tau_adapt <= rho_max_user) then
                        rho_k = rho_k * tau_adapt
                        y = y * tau_adapt        ! Keep u = y/rho constant
                        need_refactor = .true.
                        cooldown_counter = 0
                    else if (res_dual > mu_adapt * res_prim .and. rho_k / tau_adapt >= rho_min_user) then
                        rho_k = rho_k / tau_adapt
                        y = y / tau_adapt        ! Keep u = y/rho constant
                        need_refactor = .true.
                        cooldown_counter = 0
                    end if
                end if

                ! Stagnation detection
                if (abs(res_prim_prev - res_prim) < 1.0e-10_dp * res_prim_prev) then
                    stagnation_count = stagnation_count + 1
                else
                    stagnation_count = 0
                end if
                res_prim_prev = res_prim

                if (stagnation_count >= stagnation_limit) then
                    count = 0
                    do i = 1, n
                        if (abs(znew(i)) > 1.0e-6_dp) count = count + 1
                    end do
                    write(*, '(A,I8,A,ES9.2,A,ES9.2,A,I4)') &
                        '  STAGNATED iter=', jj, &
                        ' rp=', res_prim, ' rd=', res_dual, ' nnz=', count
                    z = znew
                    exit
                end if

                z = znew
            end do

            if (jj > Max_ADMM_Iter) then
                write(*, '(A,I10,A)') '  WARNING: Max iterations (', Max_ADMM_Iter, &
                    ') reached without convergence'
                z = znew
            end if

            ! Count nonzeros
            count = 0
            do i = 1, n
                if (abs(z(i)) > 1.0e-6_dp) count = count + 1
            end do
            nzval(k) = count

            ! Store ADMM solution
            Amp(:,k) = z

            ! Jval for ADMM solution
            Jval(k) = real(dot_product(conjg(matmul(conjg(z), P)), z)) &
                     - 2.0_dp * real(dot_product(q, z)) + trace_GG

            ! ---- Polishing step ----
            mm = n - count
            abz = abs(z)
            call sort_index(abz, index_s, reverse=.true.)

            if (mm /= 0) then
                allocate(E(n,mm))
                E = 0.0_dp
                do i = 1, mm
                    E(index_s(n - mm + i), i) = 1
                end do

                allocate(KKT(n+mm, n+mm), right_hand_side(n+mm), ipiv(n+mm))
                KKT = 0.0_dp
                KKT(1:n, 1:n) = P
                KKT(1:n, n+1:n+mm) = E
                KKT(n+1:n+mm, 1:n) = transpose(E)
                right_hand_side(1:n) = q
                right_hand_side(n+1:n+mm) = 0.0_dp

                call zgesv(n+mm, 1, KKT, n+mm, ipiv, right_hand_side, n+mm, info)
                if (info /= 0) then
                    sol = z
                else
                    sol = right_hand_side(1:n)
                end if

                deallocate(E, KKT, right_hand_side, ipiv)
            else
                allocate(KKT(n,n), right_hand_side(n), ipiv(n))
                KKT = P
                right_hand_side = q
                call zgesv(n, 1, KKT, n, ipiv, right_hand_side, n, info)
                if (info /= 0) then
                    sol = z
                else
                    sol = right_hand_side
                end if
                deallocate(KKT, right_hand_side, ipiv)
            end if

            xpol = sol

            Amp_p(:,k) = xpol
            Jval_p(k) = real(dot_product(conjg(matmul(conjg(xpol), P)), xpol)) &
                       - 2.0_dp * real(dot_product(q, xpol)) + trace_GG

            print *, "Gamma Value", k
            flush(6)

        end do

        close(2000)

        deallocate(Iden, Prho, Plow, y, z, u, B, xnew, v, znew, xpol, sol)
        deallocate(abz, index_s)

        print *, "Done (sparse amplitude optimization)"

    end subroutine Amplitudes_sp

end module Sparity
