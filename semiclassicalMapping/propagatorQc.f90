module propagator
    !! lamda as a function of Qc
    use params
    use functions
    implicit none

    type :: system_state
        real(kind=wp1) :: x(jlevel)
        real(kind=wp1) :: px(jlevel)
        real(kind=wp1) :: y(jlevel)
        real(kind=wp1) :: py(jlevel)
        real(kind=wp1) :: xbar(glevel, nel)
        real(kind=wp1) :: pxbar(glevel, nel)
        real(kind=wp1) :: ybar(glevel, nel)
        real(kind=wp1) :: pybar(glevel, nel)
        real(kind=wp1) :: q, p
        real(kind=wp1) :: qt, pt
        real(kind=wp1) :: qc, pc
    end type system_state

    real(kind=wp1), parameter :: q_left_def  = -15.0_wp1
    real(kind=wp1), parameter :: q_right_def =  15.0_wp1
    real(kind=wp1), parameter :: l_min_def   =  0.15_wp1
    real(kind=wp1), parameter :: l_max_def   =  0.2_wp1
    real(kind=wp1), parameter :: period = one

    real(kind=wp1), save :: mu_left = one, mus(nel) = [-one, one]

    type(system_state) :: k1, k2, k3, k4
    type(system_state) :: state_temp

    contains

    subroutine set_mu_from_arg(arg_idx)
        integer, intent(in), optional :: arg_idx
        integer :: idx, nargs
        character(len=100) :: arg

        idx = merge(arg_idx, 1, present(arg_idx))
        nargs = command_argument_count()
        if (nargs >= idx) then
            call get_command_argument(idx, arg)
            read(arg, *) mu_left
            write(*,'(a,2x,f10.4)') "Overriding the default value of mu with", mu_left
        else
            write(*,'(a,2x,f10.4)') "No arg provided for mu, using default:", mu_left
        end if

        mus = [-one, one] * mu_left
    end subroutine set_mu_from_arg

    subroutine deriv(derivatives_out, state_in, tgs, ek)
        implicit none
        type(system_state), intent(out) :: derivatives_out
        type(system_state), intent(in) :: state_in
        real(kind=wp1), intent(in) :: tgs(glevel, jlevel, nel)
        real(kind=wp1), intent(in) :: ek(glevel, nel)
        !! local
        real(kind=wp1) :: tgs_ybar_sum, tgs_pybar_sum, tgs_xbar_sum, tgs_pxbar_sum, &
                          stql, enj, ekl, local_pop(jlevel), h0, hj(jlevel), total_pop, &
                          off_diag, lamda_deriv, lam_qc, stql_qc
        integer :: j, g, l, jbar

        lam_qc = q_lambda(state_in%qc)
        stql = stwo * lam_qc * state_in%q
        lamda_deriv = stwo * derivative_q_lambda(state_in%qc) * state_in%q

        h0 = half * omega_t * (state_in%pt**2 + state_in%qt**2) &
            + half * omega_c * (state_in%pc**2 + state_in%qc**2)
        off_diag = (state_in%x(1) * state_in%py(2) - state_in%y(1) * state_in%px(2) &
                    + state_in%x(2) * state_in%py(1) - state_in%y(2) * state_in%px(1))

        do j = 1, jlevel
            local_pop(j) = state_in%x(j)*state_in%py(j) - state_in%y(j)*state_in%px(j)
            hj(j) = h0 + en(j) + kappa(j) * state_in%qt
        end do

        total_pop = sum(local_pop)
        stql_qc = stql + lamda_bar * state_in%qc

        do j = 1, jlevel
            enj = hj(j)
            jbar = 3 - j

            tgs_ybar_sum = mydot(tgs(:, j, 1), state_in%ybar(:, 1)) + mydot(tgs(:, j, 2), state_in%ybar(:, 2))
            tgs_pybar_sum = mydot(tgs(:, j, 1), state_in%pybar(:, 1)) + mydot(tgs(:, j, 2), state_in%pybar(:, 2))
            tgs_xbar_sum = mydot(tgs(:, j, 1), state_in%xbar(:, 1)) + mydot(tgs(:, j, 2), state_in%xbar(:, 2))
            tgs_pxbar_sum = mydot(tgs(:, j, 1), state_in%pxbar(:, 1)) + mydot(tgs(:, j, 2), state_in%pxbar(:, 2))
            derivatives_out%x(j)  = -enj * state_in%y(j) - tgs_ybar_sum - stql_qc * state_in%y(jbar)
            derivatives_out%px(j) = -enj * state_in%py(j) - tgs_pybar_sum - stql_qc * state_in%py(jbar)
            derivatives_out%y(j)  = enj * state_in%x(j) + tgs_xbar_sum + stql_qc * state_in%x(jbar)
            derivatives_out%py(j) = enj * state_in%px(j) + tgs_pxbar_sum + stql_qc * state_in%px(jbar)
        end do

        derivatives_out%q = omega * state_in%p
        derivatives_out%p = -omega * state_in%q - stwo * lam_qc * off_diag
        derivatives_out%qc = omega_c * state_in%pc * total_pop
        derivatives_out%qt = omega_t * state_in%pt * total_pop
        derivatives_out%pc = -omega_c * state_in%qc * total_pop - (lamda_bar + lamda_deriv) * off_diag
        derivatives_out%pt = -omega_t * state_in%qt * total_pop - kappa(1)*local_pop(1) - kappa(2)*local_pop(2)

        do l = 1, nel
            do g = 1, glevel
                ekl = ek(g, l)
                derivatives_out%xbar(g, l)  = -ekl * state_in%ybar(g, l) - sum(tgs(g, :, l) * state_in%y)
                derivatives_out%pxbar(g, l) = -ekl * state_in%pybar(g, l) - sum(tgs(g, :, l) * state_in%py)
                derivatives_out%ybar(g, l)  = ekl * state_in%xbar(g, l) + sum(tgs(g, :, l) * state_in%x)
                derivatives_out%pybar(g, l) = ekl * state_in%pxbar(g, l) + sum(tgs(g, :, l) * state_in%px)
            end do
        end do
    end subroutine deriv

    subroutine single_euler_step(s_out, s_in, k_in, step_size)
        type(system_state), intent(out) :: s_out
        type(system_state), intent(in)  :: s_in, k_in
        real(kind=wp1), intent(in)      :: step_size

        s_out%x = s_in%x + step_size * k_in%x
        s_out%px = s_in%px + step_size * k_in%px
        s_out%y = s_in%y + step_size * k_in%y
        s_out%py = s_in%py + step_size * k_in%py

        s_out%xbar = s_in%xbar + step_size * k_in%xbar
        s_out%pxbar = s_in%pxbar + step_size * k_in%pxbar
        s_out%ybar = s_in%ybar + step_size * k_in%ybar
        s_out%pybar = s_in%pybar + step_size * k_in%pybar

        s_out%q = s_in%q + step_size * k_in%q
        s_out%p = s_in%p + step_size * k_in%p

        s_out%qc = s_in%qc + step_size * k_in%qc
        s_out%pc = s_in%pc + step_size * k_in%pc

        s_out%qt = s_in%qt + step_size * k_in%qt
        s_out%pt = s_in%pt + step_size * k_in%pt
    end subroutine single_euler_step

    subroutine rk4(state, tgs, ek)
        implicit none
        type(system_state), intent(inout) :: state
        real(kind=wp1), intent(in) :: tgs(glevel, jlevel, nel)
        real(kind=wp1), intent(in) :: ek(glevel, nel)

        ! Step 1
        call deriv(k1, state, tgs, ek)
        call single_euler_step(state_temp, state, k1, dt2)

        ! Step 2
        call deriv(k2, state_temp, tgs, ek)
        call single_euler_step(state_temp, state, k2, dt2)

        ! Step 3
        call deriv(k3, state_temp, tgs, ek)
        call single_euler_step(state_temp, state, k3, dt)

        call deriv(k4, state_temp, tgs, ek)

        state%x  = state%x  + dt6 * (k1%x  + two * (k2%x  + k3%x)  + k4%x)
        state%px = state%px + dt6 * (k1%px + two * (k2%px + k3%px) + k4%px)
        state%y  = state%y  + dt6 * (k1%y  + two * (k2%y  + k3%y)  + k4%y)
        state%py = state%py + dt6 * (k1%py + two * (k2%py + k3%py) + k4%py)
        state%xbar  = state%xbar  + dt6 * (k1%xbar  + two * (k2%xbar  + k3%xbar)  + k4%xbar)
        state%pxbar = state%pxbar + dt6 * (k1%pxbar + two * (k2%pxbar + k3%pxbar) + k4%pxbar)
        state%ybar  = state%ybar  + dt6 * (k1%ybar  + two * (k2%ybar  + k3%ybar)  + k4%ybar)
        state%pybar = state%pybar + dt6 * (k1%pybar + two * (k2%pybar + k3%pybar) + k4%pybar)
        state%q = state%q + dt6 * (k1%q + two * (k2%q + k3%q) + k4%q)
        state%p = state%p + dt6 * (k1%p + two * (k2%p + k3%p) + k4%p)
        state%qt = state%qt + dt6 * (k1%qt + two * (k2%qt + k3%qt) + k4%qt)
        state%pt = state%pt + dt6 * (k1%pt + two * (k2%pt + k3%pt) + k4%pt)
        state%qc = state%qc + dt6 * (k1%qc + two * (k2%qc + k3%qc) + k4%qc)
        state%pc = state%pc + dt6 * (k1%pc + two * (k2%pc + k3%pc) + k4%pc)
    end subroutine rk4


    ! pure function q_lambda(q, q_left_max, q_right_max, lambda_min, lambda_max) result(val)
    !     real(kind=wp1), intent(in) :: q
    !     real(kind=wp1), intent(in), optional :: q_left_max, q_right_max, lambda_min, lambda_max
    !     real(kind=wp1) :: val

    !     real(kind=wp1) :: loc_q_left, loc_q_right, loc_l_min, loc_l_max

    !     loc_q_left  = q_left_def;  if (present(q_left_max))  loc_q_left  = q_left_max
    !     loc_q_right = q_right_def; if (present(q_right_max)) loc_q_right = q_right_max
    !     loc_l_min   = l_min_def;   if (present(lambda_min))  loc_l_min   = lambda_min
    !     loc_l_max   = l_max_def;   if (present(lambda_max))  loc_l_max   = lambda_max

    !     if (q < loc_q_left .or. q > loc_q_right) then
    !         val = loc_l_min
    !         return
    !     end if

    !     if (q < rzero) then
    !         val = ((loc_l_max - loc_l_min) / (- loc_q_left)) * (q - loc_q_left) + loc_l_min
    !     else
    !         val = ((loc_l_min - loc_l_max) / loc_q_right) * q + loc_l_max
    !     end if
    ! end function q_lambda


    ! pure function derivative_q_lambda(q, q_left_max, q_right_max, lambda_min, lambda_max) result(dval)
    !     real(kind=wp1), intent(in) :: q
    !     real(kind=wp1), intent(in), optional :: q_left_max, q_right_max, lambda_min, lambda_max
    !     real(kind=wp1) :: dval

    !     real(kind=wp1) :: loc_q_left, loc_q_right, loc_l_min, loc_l_max

    !     loc_q_left  = q_left_def;  if (present(q_left_max))  loc_q_left  = q_left_max
    !     loc_q_right = q_right_def; if (present(q_right_max)) loc_q_right = q_right_max
    !     loc_l_min   = l_min_def;   if (present(lambda_min))  loc_l_min   = lambda_min
    !     loc_l_max   = l_max_def;   if (present(lambda_max))  loc_l_max   = lambda_max

    !     if (q < loc_q_left .or. q > loc_q_right) then
    !         dval = rzero
    !         return
    !     end if

    !     if (q < rzero) then
    !         dval = (loc_l_max - loc_l_min) / (- loc_q_left)
    !     else
    !         dval = (loc_l_min - loc_l_max) / loc_q_right
    !     end if
    ! end function derivative_q_lambda

    pure function q_lambda(q, q_right_max, lambda_max) result(val)
        real(kind=wp1), intent(in) :: q
        real(kind=wp1), intent(in), optional :: q_right_max
        real(kind=wp1), intent(in), optional :: lambda_max
        real(kind=wp1) :: val

        real(kind=wp1) :: loc_q_right, loc_l_max

        if (present(q_right_max)) then
            loc_q_right = q_right_max
        else
            loc_q_right = 15.0_wp1
        end if

        if (present(lambda_max)) then
            loc_l_max = lambda_max
        else
            loc_l_max = 0.2_wp1
        end if

        val = loc_l_max * half * (cos(q/loc_q_right * pi * period) + one)
    end function q_lambda

    pure function derivative_q_lambda(q, q_right_max, lambda_max) result(dval)
        real(kind=wp1), intent(in) :: q
        real(kind=wp1), intent(in), optional :: q_right_max
        real(kind=wp1), intent(in), optional :: lambda_max
        real(kind=wp1) :: dval

        real(kind=wp1) :: loc_q_right, loc_l_max


        if (present(q_right_max)) then
            loc_q_right = q_right_max
        else
            loc_q_right = 15.0_wp1
        end if

        if (present(lambda_max)) then
            loc_l_max = lambda_max
        else
            loc_l_max = 0.2_wp1
        end if

        dval = -loc_l_max * period * pi/(two * loc_q_right) * sin(q/loc_q_right *pi*period)
    end function derivative_q_lambda

end module propagator
