program metodo_muller
    implicit none

    ! --- Variáveis em Dupla Precisão (64 bits) ---
    real(8) :: x0, x1, x2, x3
    real(8) :: h0, h1, a, b, c
    real(8) :: delta_0, delta_1, discriminante, den
    real(8) :: tol, erro, passo
    integer :: iter, max_iter

    ! --- Parâmetros Iniciais ---
    passo    = 0.5d0
    x0       = 6.5d0
    tol      = 1.0d-6
    max_iter = 100
    erro     = 1.0d0
    iter     = 0

    ! O Método de Muller precisa de 3 pontos para começar.
    ! Usamos o 'passo' apenas ANTES do laço para inicializar a fila.
    x1 = x0 - passo
    x2 = x1 - passo

    write(*,*) "--------------------------------------------------------------------------------"
    write(*,*) "                              METODO DE MULLER                                  "
    write(*,*) "--------------------------------------------------------------------------------"
    write(*,"(A, F8.4, F8.4, F8.4)") " Chutes iniciais (x0, x1, x2): ", x0, x1, x2
    write(*,"(A, E8.1)")             " Tolerancia:                   ", tol
    write(*,*) "--------------------------------------------------------------------------------"
    write(*,"(A)") " Iter |        x3         |       f(x3)       |       Erro        "
    write(*,*) "--------------------------------------------------------------------------------"

    open(unit=10, file='raizes.dat', status='replace')

    ! --- Laço Mestre ---
    do while (erro .gt. tol .and. iter .lt. max_iter)
        iter = iter + 1

        h0 = x1 - x0
        h1 = x2 - x1

        delta_0 = (f(x1) - f(x0)) / h0
        delta_1 = (f(x2) - f(x1)) / h1

        a = (delta_1 - delta_0) / (h1 + h0)
        b = (a * h1) + delta_1
        c = f(x2)

        discriminante = (b**2) - (4.0d0 * a * c)

        ! Bloqueio de segurança para evitar crash caso o discriminante fique negativo
        ! (O Muller pode achar raízes complexas, mas como as variáveis aqui são reais, limitamos).
        if (discriminante .lt. 0.0d0) then
            write(*,*) "Aviso: Discriminante negativo. Raízes complexas detectadas."
            exit
        end if

        ! --- O Truque Numérico do Denominador ---
        ! Escolhemos o sinal que MAXIMIZA o denominador para evitar erro de cancelamento
        if (b .ge. 0.0d0) then
            den = b + sqrt(discriminante)
        else
            den = b - sqrt(discriminante)
        end if

        ! Cálculo final da nova estimativa
        x3 = x2 - ((2.0d0 * c) / den)

        ! O erro é a diferença entre a raiz nova e a última raiz estimada (x2)
        erro = abs(x3 - x2)

        write(*,"(I5, ' | ', F15.8, ' | ', E17.8, ' | ', E17.8)") iter, x3, f(x3), erro

        ! --- A Fila Anda ---
        ! Empurramos os pontos para a próxima iteração
        x0 = x1
        x1 = x2
        x2 = x3

    end do

    write(*,*) "--------------------------------------------------------------------------------"
    write(*, '(A, F15.8)') ' Raiz final encontrada: ', x3
    write(*, '(A, I5)')    ' Iteracoes realizadas:  ', iter
    write(*,*) "--------------------------------------------------------------------------------"

    write(10, '(A, F12.8)') 'Raiz encontrada:', x3
    close(10)

contains

    ! ------------------------------------------------------------------
    ! Função Polinomial Principal: f(x) = x^3 - 13x - 12
    ! As raízes reais exatas para esta função são: -3, -1, e 4.
    ! ------------------------------------------------------------------
    real(8) function f(x_func)
        real(8), intent(in) :: x_func
        f = (x_func**3) - (13.0d0 * x_func) - 12.0d0
    end function

end program metodo_muller

! Roteiro Sumário: Este programa implementa o Método de Muller, uma técnica avançada de busca de raízes que substitui as aproximações lineares por uma interpolação quadrática (parábolas) baseada em três estimativas iniciais e suas diferenças divididas. A grande vantagem computacional deste método é a sua convergência superlinear e a capacidade de encontrar raízes complexas de forma nativa, tudo isso sem exigir o cálculo de derivadas analíticas da função principal.
