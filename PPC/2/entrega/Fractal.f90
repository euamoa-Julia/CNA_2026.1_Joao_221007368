program Bairstow_Fractal
    implicit none

    integer :: i_grid, j_grid, unit_csv
    integer, parameter :: res = 1000
    real :: r0_start, s0_start
    real, parameter :: r_min = -10.0, r_max = 10.0
    real, parameter :: s_min = -10.0, s_max = 10.0
    real a(0:100), b(0:100), c(0:100), a_base(0:100)
    real r, s, delta_r, delta_s
    real det, det_r, det_s
    real erro_r, erro_s
    real discriminante
    real tol
    real raiz(1:2)
    integer n, i, n_base
    integer iter, total_iter
    real assinatura_cor

    n_base = 6

    tol = 1.0d-5

    a_base(0) = 53
    a_base(1) = 4
    a_base(2) = 23
    a_base(3) = -1067
    a_base(4) = -12
    a_base(5) = 8
    a_base(6) = -2
    
    open(newunit=unit_csv, file='fractal_bairstow.csv', status='replace')

    write(unit_csv, *) "r0, s0, cor, iteracoes"

    print *, "Iniciando o mapeamento do fractal com limite de 1500 iterações..."

    do i_grid = 0, res

        r0_start = r_min + (r_max - r_min) * (real(i_grid) / real(res))

        do j_grid = 0, res

            s0_start = s_min + (s_max - s_min) * (real(j_grid) / real(res))

            n = n_base
            a = a_base
            r = r0_start
            s = s0_start
            assinatura_cor = 0.0
            total_iter = 0

            do while (n .gt. 2)

                iter = 0
                erro_r = 1.0
                erro_s = 1.0

                do while ((abs(erro_r).gt.tol .and. abs(erro_s).gt.tol) .and. iter .lt. 1500)

                    iter = iter + 1
                    total_iter = total_iter + 1

                    b(n) = a(n)
                    b(n-1) = a(n-1) + r * b(n)
                    c(n) = b(n)
                    c(n-1) = b(n-1) + r * c(n)

                    do i = n - 2, 0, -1

                        b(i) = a(i) + r * b(i+1) + s * b(i+2)

                    end do

                    do i = n - 2, 1, -1

                        c(i) = b(i) + r * c(i+1) + s * c(i+2)

                    end do

                    det = (c(1)*c(3)) - (c(2)*c(2))

                    if (abs(det) .lt. 1.0e-12) then

                        r = r + 0.1
                        s = s + 0.1

                        cycle

                    end if

                    det_r = ((-b(0))*c(3)) - (c(2)*(-b(1)))
                    det_s = (c(1)*(-b(1))) - ((-b(0))*c(2))

                    delta_r = det_r / det
                    delta_s = det_s / det

                    r = r + delta_r
                    s = s + delta_s

                    erro_r = (delta_r / r)
                    erro_s = (delta_s / s)

                end do

                if (iter .ge. 1500) then

                    assinatura_cor = -1.0

                    exit

                end if

                assinatura_cor = assinatura_cor + abs(r) + abs(s)

                n = n - 2

                do i = 0, n

                    a(i) = b(i+2)

                end do

            end do

            if (assinatura_cor .eq. -1.0) then

                write(unit_csv, '(4(F12.6, ","))') r0_start, s0_start, assinatura_cor, real(total_iter)

                cycle

            end if

            if (n.eq.2) then

                discriminante = (r**2) + (4.0 * s)

                if (discriminante .ge. 0.0) then

                    raiz(1) = (r + sqrt(discriminante)) / 2.0
                    raiz(2) = (r - sqrt(discriminante)) / 2.0
                    assinatura_cor = assinatura_cor + abs(raiz(1)) + abs(raiz(2))

                else

                    assinatura_cor = assinatura_cor + abs(r / 2.0) + abs(sqrt(abs(discriminante)) / 2.0)

                end if

            end if

            if (n.eq.1) then

                raiz(1) = -(a(0)/a(1))
                assinatura_cor = assinatura_cor + abs(raiz(1))

            end if

            write(unit_csv, '(4(F12.6, ","))') r0_start, s0_start, assinatura_cor, real(total_iter)

        end do

        if (mod(i_grid, 50) == 0) print *, "Calculando linha X = ", i_grid, " de ", res
    end do

    close(unit_csv)

    print *, "Mapeamento Completo Finalizado."

end program Bairstow_Fractal
