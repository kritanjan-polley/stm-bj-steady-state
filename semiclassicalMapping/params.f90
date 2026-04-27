module params
    implicit none
    use iso_fortran_env, only : real64, real32
    ! integer, parameter :: wp1 = selected_real_kind(15,37)
    integer, parameter :: wp1 = real64 !! selected_real_kind(15,37)


    !!
    real(kind=wp1), parameter :: rzero = 0.0_wp1
    real(kind=wp1), parameter :: quater = 0.25_wp1
    real(kind=wp1), parameter :: half = 0.5_wp1
    real(kind=wp1), parameter :: one = 1.0_wp1
    real(kind=wp1), parameter :: oah = 1.5_wp1
    real(kind=wp1), parameter :: two = 2.0_wp1
    real(kind=wp1), parameter :: three = 3.0_wp1
    real(kind=wp1), parameter :: four = 4.0_wp1
    real(kind=wp1), parameter :: five = 5.0_wp1
    real(kind=wp1), parameter :: six = 6.0_wp1
    real(kind=wp1), parameter :: ten = 10.0_wp1
    real(kind=wp1), parameter :: pi = four*atan(one)
    real(kind=wp1), parameter :: pi2 = pi + pi
    real(kind=wp1), parameter :: stwo = sqrt(two)


    integer, parameter :: glevel = 500
    integer, parameter :: jlevel = 2
    integer, parameter :: nel = 2
    integer, parameter :: ntraj = 100000
    real(kind=wp1), parameter :: tmax = 20.0_wp1
    real(kind=wp1), parameter :: dt = 0.005_wp1
    real(kind=wp1), parameter :: dt2 = dt / two
    real(kind=wp1), parameter :: dt6 = dt / six
    integer, parameter :: nsteps = int(tmax / dt)
    real(kind=wp1), parameter :: beta = five
    logical, parameter :: debug = .false.
    integer :: initial_state = -1


    !!
    real(kind=wp1), parameter :: en(jlevel) = [-one, one], abs_en = abs(en(2) - en(1))
    real(kind=wp1), parameter :: delE = 30.0_wp1
    real(kind=wp1), parameter :: eInterval = two * delE / real(glevel-1, kind=wp1)
    real(kind=wp1), parameter :: gamma_val = one, gamma_zpe = half, &
                                 w_val = ten, w_val2 = w_val * w_val, &
                                 omega = one * abs_en, omega_c = 0.2_wp1 * abs_en, omega_t = 0.1_wp1 * abs_en, &
                                 kappa(jlevel) = [-0.2_wp1, 0.3_wp1] * abs_en, lamda_bar = 0.05_wp1 * abs_en

    character(len=50), parameter :: photon_sampling = "AA"

end module params
