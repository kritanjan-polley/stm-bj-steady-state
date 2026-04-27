module initial
    use params
    use functions
    use propagator
    implicit none

    contains

    subroutine get_initial_state(state, ek, mus)
        implicit none
        type(system_state), intent(out) :: state
        real(kind=wp1), intent(in) :: ek(glevel, nel), mus(nel)
        ! local
        real(kind=wp1) :: rand_val, occu, randn(2), zpe_amp
        real(kind=wp1) :: rand_phases(glevel, nel)
        real(kind=wp1) :: rand_occ(glevel, nel)
        real(kind=wp1) :: occupation_mask(glevel, nel)
        integer :: j, g, l
        ! system
        do j = 1, jlevel
            call random_number(rand_val)
            rand_val = pi2 * (one - rand_val)

            state%x(j) = cos(rand_val)
            state%y(j) = sin(rand_val)
            state%px(j) = -dirac_delta(j, initial_state) * state%y(j)
            state%py(j) = dirac_delta(j, initial_state) * state%x(j)
        end do

        ! photon
        if (to_lower(photon_sampling) == "aa") then
            zpe_amp = sqrt(two * (rzero + gamma_zpe))
            call random_number(rand_val)
            rand_val = pi2 * (one - rand_val)
            state%q = zpe_amp * cos(rand_val)
            state%p = -zpe_amp * sin(rand_val)
            call random_number(rand_val)
            rand_val = pi2 * (one - rand_val)
            state%qt = zpe_amp * cos(rand_val)
            state%pt = -zpe_amp * sin(rand_val)
            call random_number(rand_val)
            rand_val = pi2 * (one - rand_val)
            state%qc = zpe_amp * cos(rand_val)
            state%pc = -zpe_amp * sin(rand_val)
        elseif (to_lower(photon_sampling) == "wigner") then
            state%q = rand_normal(one/stwo)
            state%p = rand_normal(one/stwo)
            state%qt = rand_normal(one/stwo)
            state%pt = rand_normal(one/stwo)
            state%qc = rand_normal(one/stwo)
            state%pc = rand_normal(one/stwo)
        elseif (to_lower(photon_sampling) == "spin") then
            do j = 1, 2
                randn(j) = rand_normal(one)
            end do
            randn = randn/sqrt(sum(randn*randn))
            state%q = randn(1)
            state%p = randn(2)
            do j = 1, 2
                randn(j) = rand_normal(one)
            end do
            randn = randn/sqrt(sum(randn*randn))
            state%qt = randn(1)
            state%pt = randn(2)
            do j = 1, 2
                randn(j) = rand_normal(one)
            end do
            randn = randn/sqrt(sum(randn*randn))
            state%qc = randn(1)
            state%pc = randn(2)
        else
            call print_banner("Invalid sampling mode for photon mode")
            call print_banner("Choose either action-angle (AA) or Wigner (Wigner) or spin (spin)")
            stop
        end if
        !

        ! bath
        ! do g = 1, glevel
        !     do l = 1, nel
        !         call random_number(rand_val)
        !         rand_val = one - rand_val
        !         rand_val = pi2 * rand_val
        !         state%xbar(g, l) = cos(rand_val)
        !         state%ybar(g, l) = sin(rand_val)

        !         occu = occupation(ek(g, l), mus(l))
        !         state%pxbar(g, l) = -occu * state%ybar(g, l)
        !         state%pybar(g, l) = occu * state%xbar(g, l)
        !     end do
        ! end do
        call random_number(rand_phases)
        call random_number(rand_occ)

        rand_phases = pi2 * (one - rand_phases)
        state%xbar = cos(rand_phases)
        state%ybar = sin(rand_phases)

        do concurrent (l = 1:nel)
            occupation_mask(:, l) = logistic(-beta * (ek(:, l) - mus(l)))
        end do
        occupation_mask = merge(rzero, one, rand_occ > occupation_mask)

        state%pxbar = -occupation_mask * state%ybar
        state%pybar =  occupation_mask * state%xbar

    end subroutine get_initial_state

    subroutine get_eks(ek)
        implicit none
        real(kind=wp1), intent(out) :: ek(glevel, nel)
        integer :: g, l

        do l = 1, nel
            do g = 1, glevel
                ek(g, l) = -delE + (g-1) * eInterval + mus(l)
            end do
        end do
    end subroutine get_eks

    subroutine get_tgs(tgs, ek)
        implicit none
        real(kind=wp1), intent(out) :: tgs(glevel, jlevel, nel)
        real(kind=wp1), intent(in) :: ek(glevel, nel)
        real(kind=wp1) :: common_factor
        integer :: g

        common_factor = eInterval / pi2
        ! left electrode
        do g = 1, glevel
            tgs(g, 1, 1) = sqrt(spec_fermi(ek(g, 1), mus(1)) * common_factor)
            tgs(g, 2, 1) = rzero !! sqrt(spec_fermi(ek(g, 1), mus(1)) * common_factor) !! rzero
        end do
        ! right electrode
        do g = 1, glevel
            tgs(g, 1, 2) = rzero !! sqrt(spec_fermi(ek(g, 2), mus(2)) * common_factor) !! rzero
            tgs(g, 2, 2) = sqrt(spec_fermi(ek(g, 2), mus(2)) * common_factor)
        end do
    end subroutine get_tgs

    subroutine print_params()
        implicit none

        write(*,'(a)') 'The parameters used here are:'
        write(*,'(a,2x,f10.4)') 'Fermi energy (left) =', mus(1)
        write(*,'(a,2x,f10.4)') 'Fermi energy (right) =', mus(2)
        write(*,'(a,2x,f10.4)') 'Width of the band =', two * delE
        write(*,'(a,2x,f10.4)') 'Delta (diabatic energy bias) =', abs(en(2) - en(1))
        write(*,'(a,2x,f10.4)') 'Zero point energy =', gamma_zpe
        write(*,'(a,2x,f10.4)') 'Oscillator frequency =', omega
        write(*,'(a,2x,f10.4)') 'Oscillator-system coupling =', lamda
        write(*,'(a,2x,f8.4)') 'Coupling mode frequency', omega_c
        write(*,'(a,2x,f8.4)') 'Tuning mode frequency', omega_t
        write(*,'(a,2x,f8.4)') 'Diabatic coupling', lamda_bar
        write(*,'(a,2x,f8.4,2x,a,2x,f8.4)') 'Tuning coupling, kappa1', kappa(1), "kappa2", kappa(2)
        write(*,'(a,2x,f10.4)') 'Width of spectral density =', w_val
        write(*,'(a,2x,f10.4)') 'Max time:', tmax
        write(*,'(a,2x,f10.4)') 'Time step:', dt
        write(*,'(a,2x,i8)') 'No of steps in a single trajectory', nsteps
        write(*,'(a,2x,i8)') 'No of trajectories', ntraj
        write(*,'(a,2x,a)') 'Photon mode sampling method is', photon_sampling
        if (debug .eqv. .true.) then
            call print_banner('Running in debug mode, extra files will be generated')
        else
            call print_banner('Running in production mode')
        end if
    end subroutine print_params


end module initial
