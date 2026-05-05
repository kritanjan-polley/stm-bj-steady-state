program main
    use params
    use functions
    use propagator
    use initial
    use iso_c_binding
    implicit none

    interface
        function get_peak_memory_bytes() bind(c, name="get_peak_memory_bytes")
            import c_long_long
            integer(c_long_long) :: get_peak_memory_bytes
        end function get_peak_memory_bytes
    end interface

    integer(c_long_long) :: peak_bytes
    real :: peak_mb
    type(system_state) :: current_state

    real(kind=wp1) :: ek(glevel, nel)
    real(kind=wp1) :: tgs(glevel, jlevel, nel)

    real(kind=wp1) :: current(nsteps + 1, nel), check_energy(nsteps + 1), &
                      population(nsteps + 1, jlevel), photon_number(nsteps + 1), &
                      nuclear_pop(nsteps + 1, 2), check_qc(ntraj), check_qt(ntraj)
    real(kind=wp1) :: pop(jlevel)
    real(kind=wp1) :: curr(nel)
    integer :: i, j, int64_number
    real(kind=wp1) :: start_time, end_time
    logical :: if_check_qc

    call get_walltime(start_time)
    population = rzero; current = rzero; photon_number = rzero

    call set_lambda_from_arg(1)
    call set_mu_from_arg(2)
    call get_eks(ek)
    call get_tgs(tgs, ek)
    call print_params()

    if_check_qc = .true.

    call get_walltime(end_time)
    write(*,'(a,2x,f10.4,2x,a)') 'Starting loop propagation after', end_time - start_time, 'seconds'

    do i = 1, ntraj
        call get_initial_state(current_state, ek, mus)

        call get_population(pop, current_state)
        population(1, :) = population(1, :) + pop

        call get_current(curr, current_state, tgs)
        current(1, :) = current(1, :) + curr

        photon_number(1) = photon_number(1) + get_photon_number(current_state)

        if (debug) then
            check_energy(1) = check_energy(1) + get_energy(current_state, tgs, ek)
            nuclear_pop(1,:) = nuclear_pop(1,:) + get_nuclear_pop(current_state)
        end if

        do j = 1, nsteps
            call rk4(current_state, tgs, ek)

            call get_population(pop, current_state)
            population(j + 1, :) = population(j + 1, :) + pop

            call get_current(curr, current_state, tgs)
            current(j + 1, :) = current(j + 1, :) + curr

            photon_number(j + 1) = photon_number(j + 1) + get_photon_number(current_state)

            if (debug) then
                check_energy(j + 1) = check_energy(j + 1) + get_energy(current_state, tgs, ek)
                nuclear_pop(j + 1, :) = nuclear_pop(j + 1, :) + get_nuclear_pop(current_state)
            end if
        end do

        if (if_check_qc) then
            check_qc(i) = current_state%qc
            check_qt(i) = current_state%qt
        end if

        if (mod(i, ntraj / 20) == 0) then
            call get_walltime(end_time)
            write(*,'(a,2x,i8,2x,a,2x,f14.4,2x,a)') 'Done with', i, 'trajectory after', &
                    end_time - start_time, 'seconds'
            if (any(isnan(current)) .or. any(isnan(population))) then
                write(*,'(a)') 'Population or current is blowing up'
                stop
            end if
        end if
    end do

    int64_number = rand_uniform(1000, 100000000)
    population = population / real(ntraj, kind=wp1)
    current = current / real(ntraj, kind=wp1)
    photon_number = photon_number / real(ntraj, kind=wp1)

    if (debug) then
        check_energy = check_energy / real(ntraj, kind=wp1)
        nuclear_pop = nuclear_pop / real(ntraj, kind=wp1)

        open(unit=8, file='nuclear_pop.'//trim(str(int64_number))//'.txt', status='new')
        ! open(unit=9, file='nuclear_pop.txt', status='replace')
        do i = 1, nsteps + 1
            write(8, '(10f20.8)', advance='no') (i - 1) * dt, nuclear_pop(i,:)
            write(8, *)
        end do
        close(8)

        open(unit=9, file='energy.'//trim(str(int64_number))//'.txt', status='new')
        ! open(unit=9, file='energy.txt', status='replace')
        do i = 1, nsteps + 1
            write(9, '(10f20.8)', advance='no') (i - 1) * dt, check_energy(i)
            write(9, *)
        end do
        close(9)
    end if

    open(unit=10, file='population.'//trim(str(int64_number))//'.txt', status='new')
    ! open(unit=10, file='population.txt', status='replace')
    do i = 1, nsteps + 1
        write(10, '(10f12.8)', advance='no') (i - 1) * dt, population(i, :)
        write(10, *)
    end do
    close(10)

    open(unit=11, file='current.'//trim(str(int64_number))//'.txt', status='new')
    ! open(unit=11, file='current.txt', status='replace')
    do i = 1, nsteps + 1
        write(11, '(10f12.8)', advance='no') (i - 1) * dt, current(i, :)
        write(11, *)
    end do
    close(11)

    open(unit=12, file='photon.'//trim(str(int64_number))//'.txt', status='new')
    ! open(unit=12, file='photon.txt', status='replace')
    do i = 1, nsteps + 1
        write(12, '(10f20.8)', advance='no') (i - 1) * dt, photon_number(i)
        write(12, *)
    end do
    close(12)

    if (if_check_qc) then
        open(unit=13, file='qc.qt.'//trim(str(int64_number))//'.txt', status='new')
        ! open(unit=13, file='qc.qt.txt', status='replace')
        do i = 1, ntraj
            write(13, '(2f20.8)', advance='no') check_qc(i), check_qt(i)
            write(13, *)
        end do
        close(13)
    end if

    peak_bytes = get_peak_memory_bytes()
    if (peak_bytes > 0) then
        peak_mb = real(peak_bytes) / (1024.0 * 1024.0)
        write(*, '(a,2x,f8.2,2x,a)') "Peak memory usage:", peak_mb, "MB"
    else
        call print_banner("Could not determine peak memory usage.")
    end if

    call get_walltime(end_time)
    write(*,'(a,2x,f20.4,2x,a)') 'wall time : ',(end_time-start_time),'seconds'

contains

    subroutine get_population(pop, state)
        implicit none
        real(kind=wp1), intent(out) :: pop(jlevel)
        type(system_state), intent(in) :: state

        pop = state%x * state%py - state%y * state%px
    end subroutine get_population

    subroutine get_current(curr, state, tgs)
        implicit none
        real(kind=wp1), intent(out) :: curr(nel)
        type(system_state), intent(in) :: state
        real(kind=wp1), intent(in) :: tgs(glevel, jlevel, nel)
        ! local
        real(kind=wp1) :: current_left_right(nel)
        integer :: l, j, g

        current_left_right = rzero
        do l = 1, nel
            do j = 1, jlevel
                do g = 1, glevel
                    current_left_right(l) = current_left_right(l) + &
                        tgs(g, j, l) * (state%y(j) * state%pybar(g, l) - state%py(j) * state%ybar(g, l) + &
                                        state%x(j) * state%pxbar(g, l) - state%px(j) * state%xbar(g, l))
                end do
            end do
        end do
        curr = current_left_right
    end subroutine get_current

    real(kind=wp1) function get_photon_number(state) result(num)
        implicit none
        type(system_state), intent(in) :: state
        num = half * ((state%p)**2 + (state%q)**2) - gamma_zpe
    end function get_photon_number

    real(kind=wp1) function get_energy(state, tgs, ek) result(res)
        implicit none
        type(system_state), intent(in) :: state
        real(kind=wp1), intent(in) :: tgs(glevel, jlevel, nel), ek(glevel, nel)
        ! local
        integer :: j
        res = rzero
        ! hs
        res = res + sum(en * (state%x * state%py - state%y * state%px))
        ! hb
        res = res + sum(ek * (state%xbar * state%pybar - state%ybar * state%pxbar))
        ! hsb
        do j = 1, jlevel
            res = res + sum(tgs(:, j, :) * (state%xbar * state%py(j) - state%ybar * state%px(j) &
                                          + state%x(j) * state%pybar - state%y(j) * state%pxbar))
        end do
        ! hp
        res = res + omega * half * (state%q * state%q + state%p * state%p)
        ! hps
        res = res + stwo * lamda * state%q * (state%x(1) * state%py(2) - state%y(1)*state%px(2) &
                                            + state%x(2) * state%py(1) - state%y(2) * state%px(1))
    end function get_energy

    function get_nuclear_pop(state) result(res)
        implicit none
        type(system_state), intent(in) :: state
        real(kind=wp1), dimension(2) :: res

        res(1) = (state%qc**2 + state%pc**2) * half - gamma_zpe
        res(2) = (state%qt**2 + state%pt**2) * half - gamma_zpe
    end function get_nuclear_pop
end program main
