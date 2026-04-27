module functions
    use params
    implicit none

    interface
        pure function ddot(n, dx, incx, dy, incy)
            use params, only: wp1
            implicit none
            real(kind=wp1) :: ddot
            integer, intent(in) :: n
            real(kind=wp1), intent(in) :: dx(*) ! Assumed-size array
            integer, intent(in) :: incx
            real(kind=wp1), intent(in) :: dy(*) ! Assumed-size array
            integer, intent(in) :: incy
        end function ddot
    end interface

    real(kind=wp1), save :: conversion = rzero
    logical, save :: first = .true.
    integer, save ::  count_rate, count_max

    contains

    function str(k)
        implicit none
        integer, intent(in) :: k
        character(len=20) :: str
        write (str, *) k
        str = adjustl(str)
    end function str

    integer function rand_uniform(a,b) result(res)
        implicit none
        integer, intent(in) :: a, b
        real(kind=wp1) :: temp

        call random_number(temp)
        temp = one - temp
        res = min(a,b) + floor(abs(a-b+1)*temp)
    end function rand_uniform

    elemental function my_exp(r) result(res)
        implicit none
        real(kind=wp1), intent(in) :: r
        real(kind=wp1) :: res
        res = merge(rzero, exp(r), r <= -50.0_wp1)
    end function my_exp

    ! function rand_normal(std)
    !     implicit none
    !     real(kind=wp1), intent(in)::std
    !     real(kind=wp1):: rand_normal,r,theta,temp(2)

    !     call random_number(temp)
    !     temp = one - temp
    !     r = sqrt(-two*log(temp(1)))
    !     theta = pi2*temp(2)
    !     rand_normal = std*r*cos(theta)
    ! end function rand_normal
    !
    function rand_normal(std) result(val)
        implicit none
        real(kind=wp1), intent(in) :: std
        real(kind=wp1) :: val
        !! local
        real(kind=wp1), save :: next_val
        logical, save :: have_next = .false.
        real(kind=wp1) :: r, theta, u(2)

        if (have_next) then
            val = next_val * std
            have_next = .false.
        else
            call random_number(u)
            if (u(1) < tiny(u(1))) u(1) = tiny(u(1))

            r = sqrt(-two * log(u(1)))
            theta = pi2 * u(2)

            val = r * cos(theta) * std
            next_val = r * sin(theta)
            have_next = .true.
        end if
    end function rand_normal

    pure real(kind=wp1) function dirac_delta(i, j)
        integer, intent(in) :: i, j
        dirac_delta = merge(one, rzero, i == j)
    end function dirac_delta

    elemental real(kind=wp1) function logistic(x)
        real(kind=wp1), intent(in) :: x
        logistic = one / (one + exp(-x))
    end function logistic

    pure real(kind=wp1) function spec_fermi(energy, mu)
        real(kind=wp1), intent(in) :: energy, mu
        spec_fermi = gamma_val * w_val2 /((energy-mu)**2 + w_val2)
    end function spec_fermi

    pure real(kind=wp1) function fermi(energy, mu)
        real(kind=wp1), intent(in) :: energy, mu
        fermi = logistic(-beta * (energy - mu))
    end function fermi

    real(kind=wp1) function occupation(energy, mu)
        real(kind=wp1), intent(in) :: energy, mu
        real(kind=wp1) :: rand_val
        call random_number(rand_val)
        occupation = merge(rzero, one, rand_val > fermi(energy, mu))
    end function occupation

    subroutine get_walltime(wctime)
        real(kind=wp1), intent(inout) :: wctime
        integer :: count
        if (first) then
            first = .false.
            call system_clock(count, count_rate, count_max)
            conversion = one / real(count_rate, kind=wp1)
        else
            call system_clock(count)
        end if
        wctime = count * conversion
    end subroutine get_walltime

    ! real(kind=wp1) function mydot(dx, dy) result(res)
    !     !! BLAS dot
    !     implicit none
    !     real(kind=wp1), intent(in) :: dx(:), dy(:)
    !     !! local
    !     real(kind=wp1) :: dtemp
    !     integer :: i, m, mp1, n

    !     n = size(dx)
    !     res = rzero
    !     dtemp = rzero

    !     m = mod(n, 5)
    !     if (m .ne. 0) then
    !         do i = 1, m
    !             dtemp = dtemp + dx(i)*dy(i)
    !         end do
    !         if (n .lt. 5) then
    !             res = dtemp
    !             return
    !         end if
    !     end if

    !     mp1 = m + 1
    !     do i = mp1, n, 5
    !         dtemp = dtemp + dx(i)*dy(i) + dx(i+1)*dy(i+1) + &
    !             dx(i+2)*dy(i+2) + dx(i+3)*dy(i+3) + dx(i+4)*dy(i+4)
    !     end do
    !     res = dtemp
    ! end function mydot

    function mydot(x, y)
        implicit none
        real(kind=wp1), intent(in) :: x(:), y(:)
        real(kind=wp1) :: mydot

        mydot = ddot(size(x), x, 1, y, 1)
    end function mydot

    ! function myaxpy(da, dx, dy) result(res)
    !     !! blas axpy. res = da * dx + dy
    !     implicit none
    !     real(kind=wp1), intent(in) :: da
    !     real(kind=wp1), intent(in) :: dx(:), dy(:)
    !     real(kind=wp1), intent(out) :: res(:)
    !     ! local
    !     integer :: i, m, mp1, n

    !     n = size(dx)
    !     m = mod(n, 4)
    !     if (m .ne. 0) then
    !         do i = 1, m
    !             res(i) = dy(i) + da*dx(i)
    !         end do
    !     end if
    !     if (n .lt. 4) return
    !     mp1 = m + 1
    !     do i = mp1, n, 4
    !        res(i) = dy(i) + da*dx(i)
    !        res(i+1) = dy(i+1) + da*dx(i+1)
    !        res(i+2) = dy(i+2) + da*dx(i+2)
    !        res(i+3) = dy(i+3) + da*dx(i+3)
    !     end do
    ! end function myaxpy

    subroutine print_banner(text)
        implicit none
        character(len=*), intent(in) :: text
        !! local
        integer :: text_length
        character(len=:), allocatable :: border

        text_length = len_trim(text)
        border = repeat("#",8) // repeat("#", text_length) // repeat("#",8)

        print '(/a)'
        print '(a)', border
        print '(a, a, a)', "|  ", trim(text), "  |"
        print '(a)', border
        print '(a/)'
    end subroutine print_banner

    function to_upper(str) result(upper_str)
        implicit none
        character(len=*), intent(in) :: str
        character(len=len(str)) :: upper_str
        integer :: i, ich

        upper_str = str
        do i = 1, len(str)
            ich = iachar(str(i:i))
            if (ich >= iachar('a') .and. ich <= iachar('z')) then
                upper_str(i:i) = achar(ich - 32)
            end if
        end do
    end function to_upper

    function to_lower(str) result(lower_str)
        implicit none
        character(len=*), intent(in) :: str
        character(len=len(str)) :: lower_str
        integer :: i, ich

        lower_str = str
        do i = 1, len(str)
            ich = iachar(str(i:i))
            if (ich >= iachar('A') .and. ich <= iachar('Z')) then
                lower_str(i:i) = achar(ich + 32)
            end if
        end do
    end function to_lower

end module functions
