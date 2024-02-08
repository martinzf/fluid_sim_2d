module fluid_sim
    use, intrinsic :: iso_c_binding 
    implicit none
    include 'fftw3.f03'

    public :: vort_solve
    private :: meshgrid, step, ETD_func, nonlinear

    contains

    function vort_solve(Nx, Ny, Lx, Ly, nu, dt, steps, w0) result(w_cum)
        double precision, parameter :: pi = 4.d0 * atan(1.d0)
        integer, intent(in) :: Nx, Ny, steps
        double precision, intent(in) :: Lx, Ly, nu, dt
        type(C_PTR) :: plan_forward, plan_backward
        double precision, intent(in) :: w0(Ny, Nx)
        real(C_DOUBLE) :: w(Ny, Nx)
        complex(C_DOUBLE_COMPLEX) :: w_hat(Ny/2+1, Nx)
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
        plan_forward = fftw_plan_dft_r2c_2d(Nx, Ny, w, w_hat, FFTW_ESTIMATE)
        plan_backward = fftw_plan_dft_c2r_2d(Nx, Ny, w_hat, w, FFTW_ESTIMATE)
        
        ! Integration loop
        w = w0
        w_cum(1, :, :) = w
        do i = 2, steps
            call fftw_execute_dft_r2c(plan_forward, w, w_hat)
            call step(w_hat, Kxg, Kyg, K2, nu, dt, plan_forward, plan_backward)
            call fftw_execute_dft_c2r(plan_backward, w_hat, w)
            w = w / (Nx * Ny) ! Normalisation
            w_cum(i, :, :) = w
        end do
        
        ! Destroy Fourier transform plans
        call fftw_destroy_plan(plan_forward)
        call fftw_destroy_plan(plan_backward)

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
    subroutine step(w, kx, ky, k2, nu, dt, plan_forward, plan_backward)
        complex(C_DOUBLE_COMPLEX) :: w(:, :)
        double precision :: kx(:, :), ky(:, :), k2(:, :), nu, dt
        type(C_PTR) :: plan_forward, plan_backward

        w = exp(-nu * k2 * dt) * w +&
            ETD_func(-nu * k2, dt) * nonlinear(w, kx, ky, k2, plan_forward, plan_backward)

    end subroutine step

    ! Numerically problematic function in ETD scheme
    function ETD_func(x, a) result(phi)
        double precision :: x(:, :), a
        double precision :: phi(size(x, 1), size(x, 2))

        x(1, 1) = epsilon(1.d0)
        phi = (exp(a * x) - 1) / x

    end function ETD_func

    ! Nonlinear term of the PDE
    function nonlinear(w, kx, ky, k2, plan_forward, plan_backward) result(N)
        complex(C_DOUBLE_COMPLEX) :: w(:, :)
        double precision :: kx(:, :), ky(:, :), k2(:, :)
        type(C_PTR) :: plan_forward, plan_backward
        double complex :: N(size(w, 1), size(w, 2))
        complex(C_DOUBLE_COMPLEX) :: alpha(size(w, 1), size(w, 2)), beta(size(w, 1), size(w, 2))
        complex(C_DOUBLE_COMPLEX) :: gamma(size(w, 1), size(w, 2)), delta(size(w, 1), size(w, 2))
        real(C_DOUBLE) :: A(2*(size(w, 1) - 1), size(w, 2)), B(2*(size(w, 1) - 1), size(w, 2))
        real(C_DOUBLE) :: C(2*(size(w, 1) - 1), size(w, 2)), D(2*(size(w, 1) - 1), size(w, 2))
        real(C_DOUBLE) :: x(2*(size(w, 1) - 1), size(w, 2)), y(2*(size(w, 1) - 1), size(w, 2))
        complex(C_DOUBLE_COMPLEX) :: conv1(size(w, 1), size(w, 2)), conv2(size(w, 1), size(w, 2))

        ! CONVOLUTIONS
        ! Inverse transforms
        k2(1, 1) = 1 ! Avoid divide by zero
        alpha = ky * w / k2
        beta = kx * w
        gamma = kx * w / k2
        delta = ky * w
        call fftw_execute_dft_c2r(plan_backward, alpha, A)
        call fftw_execute_dft_c2r(plan_backward, beta, B)
        call fftw_execute_dft_c2r(plan_backward, gamma, C)
        call fftw_execute_dft_c2r(plan_backward, delta, D)
        ! Forward transforms
        A = A / size(A) ! Normalisation
        B = B / size(B)
        C = C / size(C)
        D = D / size(D)
        x = A * B 
        y = C * D
        call fftw_execute_dft_r2c(plan_forward, x, conv1)
        call fftw_execute_dft_r2c(plan_forward, y, conv2)

        N = conv1 - conv2

    end function nonlinear

end module fluid_sim