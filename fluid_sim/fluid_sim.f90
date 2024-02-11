module fluid_sim
    use, intrinsic :: iso_c_binding 
    implicit none
    include 'fftw3.f03'

    public :: vort_solve
    private :: meshgrid, step, ETD_func, nonlinear, pad

    contains

    function vort_solve(Nx, Ny, Lx, Ly, nu, dt, steps, w0) result(w_cum)
        double precision, parameter :: pi = 4.d0 * atan(1.d0)
        integer, intent(in) :: Nx, Ny, steps
        double precision, intent(in) :: Lx, Ly, nu, dt
        type(C_PTR) :: forward, backward, forward_padded, backward_padded
        double precision, intent(in) :: w0(Ny, Nx)
        real(C_DOUBLE) :: w(Ny, Nx)
        complex(C_DOUBLE_COMPLEX) :: w_hat(Ny/2+1, Nx)
        real(C_DOUBLE) :: p(3*Ny/2, 3*Nx/2)
        complex(C_DOUBLE_COMPLEX) :: p_hat(3*Ny/4+1, 3*Nx/2)
        double precision :: w_cum(steps, Ny, Nx)
        integer :: ix(Nx), iy(Ny/2+1)
        double precision :: kx(Nx), ky(Ny/2+1)
        double precision :: Kxg(Ny/2+1, Nx), Kyg(Ny/2+1, Nx), K2(Ny/2+1, Nx)
        integer :: i

        ! Set up wave vectors
        ix(:Nx/2+1) = [(i, i = 0, Nx/2)]
        ix(Nx/2+1:) = [(i, i = -Nx/2, -1)]
        kx = 2 * pi * ix / Lx
        iy = [(i, i = 0, Ny/2)] 
        ky = 2 * pi * iy / Ly
        call meshgrid(kx, ky, Kxg, Kyg)
        K2 = Kxg ** 2 + Kyg ** 2

        ! Create Fourier transform plans
        forward  = fftw_plan_dft_r2c_2d(Nx, Ny, w, w_hat, FFTW_ESTIMATE)
        backward = fftw_plan_dft_c2r_2d(Nx, Ny, w_hat, w, FFTW_ESTIMATE)
        forward_padded  = fftw_plan_dft_r2c_2d(3*Nx/2, 3*Ny/2, p, p_hat, FFTW_ESTIMATE)
        backward_padded = fftw_plan_dft_c2r_2d(3*Nx/2, 3*Ny/2, p_hat, p, FFTW_ESTIMATE)
        
        ! Integration loop
        w = w0
        w_cum(1, :, :) = w
        do i = 2, steps
            call fftw_execute_dft_r2c(forward, w, w_hat)
            call step(w_hat, Nx, Ny, Kxg, Kyg, K2, nu, dt, forward_padded, backward_padded)
            call fftw_execute_dft_c2r(backward, w_hat, w)
            w = w / (Nx * Ny) ! Normalisation
            w_cum(i, :, :) = w
        end do
        
        ! Destroy Fourier transform plans
        call fftw_destroy_plan(forward)
        call fftw_destroy_plan(backward)
        call fftw_destroy_plan(forward_padded)
        call fftw_destroy_plan(backward_padded)

    end function vort_solve

    ! 2D meshgrid of 1D arrays
    subroutine meshgrid(x, y, Xg, Yg)
        double precision :: x(:), y(:)
        double precision :: Xg(size(y), size(x)), Yg(size(y), size(x))
        integer :: i, j

        ! Generating meshgrid
        do i = 1, size(x)
            do j = 1, size(y)
                Xg(j, i) = x(i)
                Yg(j, i) = y(j)
            end do
        end do
        
    end subroutine meshgrid

    ! Single time step as per ETD scheme
    subroutine step(w_hat, Nx, Ny, kx, ky, k2, nu, dt, forward_padded, backward_padded)
        complex(C_DOUBLE_COMPLEX) :: w_hat(:, :)
        integer :: Nx, Ny
        double precision :: kx(:, :), ky(:, :), k2(:, :), nu, dt
        type(C_PTR) :: forward_padded, backward_padded

        w_hat = exp(-nu * k2 * dt) * w_hat +&
            ETD_func(-nu * k2, dt) * nonlinear(w_hat, Nx, Ny, kx, ky, k2, forward_padded, backward_padded)

    end subroutine step

    ! Numerically problematic function in ETD scheme
    function ETD_func(x, a) result(phi1)
        double precision :: x(:, :), a
        double precision :: phi1(size(x, 1), size(x, 2))
        double precision, dimension(7), parameter :: &
        p = (/1.d0, 1.d0/26.d0, 5.d0/156.d0, 1.d0/858.d0, 1.d0/5720.d0, 1.d0/205920.d0, 1.d0/8648640.d0/)
        double precision, dimension(7), parameter :: &
        q = (/1.d0, -6.d0/13.d0, 5.d0/52.d0, -5.d0/429.d0, 1.d0/1144.d0, -1.d0/25740.d0, 1.d0/1235520.d0/)
        double precision :: pval_p, pval_q
        double precision :: x_powers(7)
        integer :: i, j, k

        do i = 1, size(x, 1)
            do j = 1, size(x, 2)
                ! Padé approximant
                if (abs(a * x(i, j)) < 1) then 
                    x_powers = [((a * x(i, j))**k, k = 0, 6)]
                    pval_p = sum(p * x_powers)
                    pval_q = sum(q * x_powers)
                    phi1(i, j) = a * pval_p / pval_q
                ! Regular evaluation
                else 
                    phi1(i, j) = (exp(a * x(i, j)) - 1) / x(i, j)
                end if
            end do
        end do

    end function ETD_func

    ! Nonlinear term of the PDE
    function nonlinear(w_hat, Nx, Ny, kx, ky, k2, forward_padded, backward_padded) result(N)
        complex(C_DOUBLE_COMPLEX) :: w_hat(:, :)
        integer :: Nx, Ny
        double precision :: kx(:, :), ky(:, :), k2(:, :)
        type(C_PTR) :: forward_padded, backward_padded
        double complex :: N(Ny/2+1, Nx)
        complex(C_DOUBLE_COMPLEX) :: alpha(3*Ny/4+1, 3*Nx/2), beta(3*Ny/4+1, 3*Nx/2)
        complex(C_DOUBLE_COMPLEX) :: gamm(3*Ny/4+1, 3*Nx/2), delta(3*Ny/4+1, 3*Nx/2)
        real(C_DOUBLE) :: A(3*Ny/2, 3*Nx/2), B(3*Ny/2, 3*Nx/2)
        real(C_DOUBLE) :: C(3*Ny/2, 3*Nx/2), D(3*Ny/2, 3*Nx/2)
        real(C_DOUBLE) :: x(3*Ny/2, 3*Nx/2), y(3*Ny/2, 3*Nx/2)
        complex(C_DOUBLE_COMPLEX) :: conv1(3*Ny/4+1, 3*Nx/2), conv2(3*Ny/4+1, 3*Nx/2)

        ! CONVOLUTIONS
        ! Inverse transforms
        k2(1, 1) = 1 ! Avoid division by zero
        alpha = pad(ky * w_hat / k2, Nx, Ny)
        beta  = pad(kx * w_hat     , Nx, Ny)
        gamm  = pad(kx * w_hat / k2, Nx, Ny)
        delta = pad(ky * w_hat     , Nx, Ny)
        call fftw_execute_dft_c2r(backward_padded, alpha, A)
        call fftw_execute_dft_c2r(backward_padded, beta , B)
        call fftw_execute_dft_c2r(backward_padded, gamm , C)
        call fftw_execute_dft_c2r(backward_padded, delta, D)
        ! Forward transforms
        A = A / size(A) ! Normalisation
        B = B / size(B)
        C = C / size(C)
        D = D / size(D)
        x = A * B 
        y = C * D
        call fftw_execute_dft_r2c(forward_padded, x, conv1)
        call fftw_execute_dft_r2c(forward_padded, y, conv2)

        N(:, 1:Nx/2) = conv1(1:Ny/2+1, 1:Nx/2) - conv2(1:Ny/2+1, 1:Nx/2)
        N(:, Nx/2+1:) = conv1(1:Ny/2+1, Nx+1:) - conv2(1:Ny/2+1, Nx+1:)

    end function nonlinear

    ! Zero padding arrays for dealiasing
    function pad(arr, Nx, Ny) result(arr_padded)
        double complex :: arr(:, :)
        integer :: Nx, Ny
        double complex :: arr_padded(3*Ny/4+1, 3*Nx/2)

        arr_padded(1:Ny/2+1, 1:Nx/2) = arr(:, 1:Nx/2)
        arr_padded(1:Ny/2+1, Nx+1:) = arr(:, Nx/2+1:)
        arr_padded(Ny/2+2:, :) = 0
        arr_padded(:, Nx/2+1:Nx) = 0
        
    end function pad

end module fluid_sim