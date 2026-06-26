program bisseccao_falsa_posicao
    implicit none

    real(8) :: xi, xs, xm, fm
    real(8) :: fi, fs, xr, fr
    real(8) :: tol, erro, errov
    real(8) :: xmd, xrd, raizv
    integer :: iter

    ! --- 1. Definições Globais ---
    ! Raiz verdadeira conhecida (para fins de cálculo do erro verdadeiro)
    raizv = 14.780208593679468d0

    write(*,*) "================================================================================"
    write(*,*) "                 PARTE 1: METODO DA BISSECCAO (BUSCA BINARIA)                   "
    write(*,*) "================================================================================"

    ! --- 2. Inicialização: Bissecção ---
    xi   = 12.0d0
    xs   = 16.0d0
    tol  = 1.0d-5
    erro = 1.0d0
    iter = 0

    write(*,"(A, F8.4, F8.4)") " Limites iniciais (xi, xs): ", xi, xs
    write(*,"(A, E8.1)")       " Tolerancia:                ", tol
    write(*,*) "--------------------------------------------------------------------------------"
    write(*,"(A)") " Iter |       x_medio      |        f(xm)       |    Erro Relativo   "
    write(*,*) "--------------------------------------------------------------------------------"

    open(unit=1, file='bisseccao.dat', status='replace')

    ! --- 3. LOOP PRINCIPAL: Bissecção ---
    do while(erro .gt. tol)

        ! Avalia os valores da função nas fronteiras do intervalo
        fi = f(xi)
        fs = f(xs)

        ! Averígua se a função muda de sinal no intervalo (Garantia de raiz)
        if ((fi * fs) .lt. 0.0d0) then

            ! Encontra o ponto médio do intervalo
            xm = (xi + xs) / 2.0d0
            fm = f(xm)

            ! Verifica em qual das duas metades a raiz se encontra
            if ((fi * fm) .lt. 0.0d0) then
                ! A raiz está na metade inferior (entre xi e xm)
                xs = xm
            else
                ! A raiz está na metade superior (entre xm e xs)
                xi = xm
            end if

            ! Calcula o novo ponto médio para estimar o erro iterativo
            xmd = (xi + xs) / 2.0d0

            erro  = abs(xm - xmd)       ! Erro relativo aproximado
            errov = abs(xmd - raizv)    ! Erro verdadeiro

            iter = iter + 1
            write(*,"(I5, ' | ', F18.8, ' | ', E18.8, ' | ', E18.8)") iter, xm, fm, erro
            write(1,*) iter, xi, xs, xm, erro, errov
        else
            write(*,*) "Erro: Não há raiz garantida neste intervalo (f(xi) * f(xs) > 0)."
            exit
        end if
    end do
    close(1)

    write(*,*) "--------------------------------------------------------------------------------"
    write(*, '(A, F15.8)') ' Raiz final pela Bisseccao: ', xm
    write(*,*) ""


    write(*,*) "================================================================================"
    write(*,*) "                     PARTE 2: METODO DA FALSA POSICAO                           "
    write(*,*) "================================================================================"

    ! --- 4. Inicialização: Falsa Posição ---
    ! Resetando as variáveis para a segunda parte do programa
    xi    = 12.0d0
    xs    = 16.0d0
    tol   = 1.0d-5
    erro  = 1.0d0
    errov = 0.0d0
    iter  = 0

    write(*,"(A, F8.4, F8.4)") " Limites iniciais (xi, xs): ", xi, xs
    write(*,"(A, E8.1)")       " Tolerancia:                ", tol
    write(*,*) "--------------------------------------------------------------------------------"
    write(*,"(A)") " Iter |        x_reta      |        f(xr)       |    Erro Relativo   "
    write(*,*) "--------------------------------------------------------------------------------"

    open(unit=2, file='falsa_posicao.dat', status='replace')

    ! --- 5. LOOP PRINCIPAL: Falsa Posição ---
    do while (erro .gt. tol)

        ! Avalia a função nos limites atuais
        fi = f(xi)
        fs = f(xs)

        if ((fi * fs) .lt. 0.0d0) then

            ! Em vez do ponto médio, traça uma reta entre (xi, fi) e (xs, fs)
            ! e encontra o ponto onde essa reta cruza o eixo X (xr)
            xr = xs - (fs * (xi - xs)) / (fi - fs)
            fr = f(xr)

            ! Verifica em qual subintervalo a raiz se encontra
            if ((fi * fr) .lt. 0.0d0) then
                xs = xr
            else
                xi = xr
            end if

            ! Recalcula as imagens e o próximo cruzamento para definir o erro
            fi = f(xi)
            fs = f(xs)
            xrd = xs - (fs * (xi - xs)) / (fi - fs)

            erro  = abs(xr - xrd)       ! Erro relativo aproximado
            errov = abs(xrd - raizv)    ! Erro verdadeiro

            iter = iter + 1
            write(*,"(I5, ' | ', F18.8, ' | ', E18.8, ' | ', E18.8)") iter, xr, fr, erro
            write(2,*) iter, xi, xs, xr, erro, errov
        else
            write(*,*) "Erro: Não há raiz garantida neste intervalo (f(xi) * f(xs) > 0)."
            exit
        end if
    end do
    close(2)

    write(*,*) "--------------------------------------------------------------------------------"
    write(*, '(A, F15.8)') ' Raiz final pela Falsa Pos.:', xr
    write(*,*) "================================================================================"

contains

    ! --- 6. Função Objetivo ---
    ! f(x) = (667.38/x)*(1.0 - exp(-0.146843*x)) - 40.0
    real(8) function f(x_func)
        real(8), intent(in) :: x_func
        f = (667.38d0 / x_func) * (1.0d0 - exp(-0.146843d0 * x_func)) - 40.0d0
    end function f

end program bisseccao_falsa_posicao

