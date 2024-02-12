module fluid_sim
    use, intrinsic :: iso_c_binding 
    implicit none
    include 'fftw3.f03'

    public :: vort_solve
    private :: alloc_real, alloc_complex, meshgrid, step, ETD_func, nonlinear, pad

    ! FFT plans
    type(C_PTR) :: forward, backward, forward_padded, backward_padded 
    ! FFT variables
    real(C_DOUBLE), pointer :: w(:, :)
    complex(C_DOUBLE_COMPLEX), pointer :: w_hat(:, :)
    complex(C_DOUBLE_COMPLEX), pointer :: A_hat(:, :), B_hat(:, :)
    complex(C_DOUBLE_COMPLEX), pointer :: C_hat(:, :), D_hat(:, :)
    real(C_DOUBLE), pointer :: A(:, :), B(:, :)
    real(C_DOUBLE), pointer :: C(:, :), D(:, :)
    real(C_DOUBLE), pointer :: x(:, :), y(:, :)
    complex(C_DOUBLE_COMPLEX), pointer :: conv1(:, :), conv2(:, :)

    contains

    function vort_solve(Nx, Ny, Lx, Ly, nu, dt, steps, w0) result(w_cum)
        double precision, parameter :: pi = 4.d0 * atan(1.d0)
        integer, intent(in) :: Nx, Ny, steps
        double precision, intent(in) :: Lx, Ly, nu, dt
        double precision, intent(in) :: w0(Ny, Nx)
        double precision :: w_cum(steps, Ny, Nx)
        ! Wavenumbers
        integer :: ix(Nx), iy(Ny/2+1)
        double precision :: kx(Nx), ky(Ny/2+1)
        double precision :: Kxg(Ny/2+1, Nx), Kyg(Ny/2+1, Nx), K2(Ny/2+1, Nx)
        ! Indexing
        integer :: i

        ! SET UP FOURIER TRANSFORMS
        ! Allocate FFT memory
        type(C_PTR) :: p_w, p_w_hat ! Pointers
        type(C_PTR) :: p_A_hat, p_B_hat, p_C_hat, p_D_hat, p_A, p_B, p_C, p_D, p_x, p_y, p_conv1, p_conv2
        call alloc_real(w, p_w, Ny, Nx)
        call alloc_complex(w_hat, p_w_hat, Ny, Nx)
        call alloc_complex(A_hat, p_A_hat, 3*Ny/2, 3*Nx/2)
        call alloc_complex(B_hat, p_B_hat, 3*Ny/2, 3*Nx/2)
        call alloc_complex(C_hat, p_C_hat, 3*Ny/2, 3*Nx/2)
        call alloc_complex(D_hat, p_D_hat, 3*Ny/2, 3*Nx/2)
        call alloc_real(A, p_A, 3*Ny/2, 3*Nx/2)
        call alloc_real(B, p_B, 3*Ny/2, 3*Nx/2)
        call alloc_real(C, p_C, 3*Ny/2, 3*Nx/2)
        call alloc_real(D, p_D, 3*Ny/2, 3*Nx/2)
        call alloc_real(x, p_x, 3*Ny/2, 3*Nx/2)
        call alloc_real(y, p_y, 3*Ny/2, 3*Nx/2)
        call alloc_complex(conv1, p_conv1, 3*Ny/2, 3*Nx/2)
        call alloc_complex(conv2, p_conv2, 3*Ny/2, 3*Nx/2)
        ! Create Fourier transform plans
        forward         = fftw_plan_dft_r2c_2d(Nx    , Ny    , w, w_hat, FFTW_ESTIMATE)
        backward        = fftw_plan_dft_c2r_2d(Nx    , Ny    , w_hat, w, FFTW_ESTIMATE)
        forward_padded  = fftw_plan_dft_r2c_2d(3*Nx/2, 3*Ny/2, A, A_hat, FFTW_ESTIMATE)
        backward_padded = fftw_plan_dft_c2r_2d(3*Nx/2, 3*Ny/2, A_hat, A, FFTW_ESTIMATE)

        ! SET UP WAVE VECTORS
        ix(:Nx/2+1) = [(i, i = 0, Nx/2)]
        ix(Nx/2+1:) = [(i, i = -Nx/2, -1)]
        kx = 2 * pi * ix / Lx
        iy = [(i, i = 0, Ny/2)] 
        ky = 2 * pi * iy / Ly
        call meshgrid(kx, ky, Kxg, Kyg)
        K2 = Kxg ** 2 + Kyg ** 2
        
        ! INTEGRATION LOOP
        w = w0
        w_cum(1, :, :) = w
        do i = 2, steps
            call fftw_execute_dft_r2c(forward, w, w_hat)
            call step(Nx, Ny, Kxg, Kyg, K2, nu, dt)
            call fftw_execute_dft_c2r(backward, w_hat, w)
            w = w / (Nx * Ny) ! Normalisation
            w_cum(i, :, :) = w
        end do
        
        ! END
        ! Destroy Fourier transform plans
        call fftw_destroy_plan(forward)
        call fftw_destroy_plan(backward)
        call fftw_destroy_plan(forward_padded)
        call fftw_destroy_plan(backward_padded)
        ! Release FFT variables memory
        call fftw_free(p_w)
        call fftw_free(p_w_hat)
        call fftw_free(p_A_hat)
        call fftw_free(p_B_hat)
        call fftw_free(p_C_hat)
        call fftw_free(p_D_hat)
        call fftw_free(p_A)
        call fftw_free(p_B)
        call fftw_free(p_C)
        call fftw_free(p_D)
        call fftw_free(p_x)
        call fftw_free(p_y)
        call fftw_free(p_conv1)
        call fftw_free(p_conv2)

    end function vort_solve

    ! Allocating memory of real array
    subroutine alloc_real(var, p_var, L, M)
        real(C_DOUBLE), pointer, intent(out) :: var(:, :)
        type(C_PTR), intent(out) :: p_var
        integer, intent(in) :: L, M

        p_var =  fftw_alloc_real(int(L * M, C_SIZE_T))
        call c_f_pointer(p_var, var, [L, M])
        
    end subroutine alloc_real

    ! Allocating memory of complex array
    subroutine alloc_complex(var, p_var, L, M)
        complex(C_DOUBLE_COMPLEX), pointer, intent(out) :: var(:, :)
        type(C_PTR), intent(out) :: p_var
        integer, intent(in) :: L, M

        p_var =  fftw_alloc_complex(int((L/2 + 1) * M, C_SIZE_T))
        call c_f_pointer(p_var, var, [L/2 + 1, M])
        
    end subroutine alloc_complex

    ! 2D meshgrid of 1D arrays
    subroutine meshgrid(x, y, Xg, Yg)
        double precision, intent(in) :: x(:), y(:)
        double precision, intent(out) :: Xg(size(y), size(x)), Yg(size(y), size(x))
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
    subroutine step(Nx, Ny, kx, ky, k2, nu, dt)
        integer, intent(in) :: Nx, Ny
        double precision, intent(in) :: kx(:, :), ky(:, :), nu, dt
        double precision, intent(inout) :: k2(:, :)

        w_hat = exp(-nu * k2 * dt) * w_hat + &
                ETD_func(-nu * k2, dt) * nonlinear(Nx, Ny, kx, ky, k2)

    end subroutine step

    ! Numerically problematic function in ETD scheme
    function ETD_func(x, a) result(phi1)
        double precision, intent(in) :: x(:, :), a
        double precision :: phi1(size(x, 1), size(x, 2))
        integer, parameter :: n_terms = 7
        double precision, dimension(n_terms), parameter :: &
        p = (/1.d0, 1.d0/26.d0, 5.d0/156.d0, 1.d0/858.d0, 1.d0/5720.d0, 1.d0/205920.d0, 1.d0/8648640.d0/)
        double precision, dimension(n_terms), parameter :: &
        q = (/1.d0, -6.d0/13.d0, 5.d0/52.d0, -5.d0/429.d0, 1.d0/1144.d0, -1.d0/25740.d0, 1.d0/1235520.d0/)
        double precision :: pval_p, pval_q
        double precision :: x_powers(n_terms)
        integer :: i, j, k

        do i = 1, size(x, 1)
            do j = 1, size(x, 2)
                ! Padé approximant
                if (abs(a * x(i, j)) < 1) then 
                    x_powers = [((a * x(i, j))**k, k = 0, n_terms-1)]
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
    function nonlinear(Nx, Ny, kx, ky, k2) result(N)
        integer, intent(in) :: Nx, Ny
        double precision, intent(in) :: kx(:, :), ky(:, :)
        double precision, intent(inout) :: k2(:, :)
        double complex :: N(Ny/2+1, Nx)

        ! CONVOLUTIONS
        ! Inverse transforms
        k2(1, 1) = 1 ! Avoid division by zero
        A_hat = pad(ky * w_hat / k2, Nx, Ny)
        B_hat = pad(kx * w_hat     , Nx, Ny)
        C_hat = pad(kx * w_hat / k2, Nx, Ny)
        D_hat = pad(ky * w_hat     , Nx, Ny)
        k2(1, 1) = 0.d0 ! Restore value
        call fftw_execute_dft_c2r(backward_padded, A_hat, A)
        call fftw_execute_dft_c2r(backward_padded, B_hat, B)
        call fftw_execute_dft_c2r(backward_padded, C_hat, C)
        call fftw_execute_dft_c2r(backward_padded, D_hat, D)
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
        double complex, intent(in) :: arr(:, :)
        integer, intent(in) :: Nx, Ny
        double complex :: arr_padded(3*Ny/4+1, 3*Nx/2)

        arr_padded(1:Ny/2+1, 1:Nx/2) = arr(:, 1:Nx/2)
        arr_padded(1:Ny/2+1, Nx+1:) = arr(:, Nx/2+1:)
        arr_padded(Ny/2+2:, :) = 0
        arr_padded(:, Nx/2+1:Nx) = 0
        
    end function pad

end module fluid_sim