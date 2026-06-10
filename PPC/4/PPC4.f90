program PPC4
    implicit none

    ! Declaração de variáveis do programa principal
    real(8) :: x_0, y_0                 ! Coordenadas do ponto inicial
    real(8) :: x_grid, y_grid, z_grid   ! Coordenadas espaciais para a malha 3D
    integer :: metodo, i, j             ! Opção do usuário e contadores para a malha

    ! ==============================================================================
    ! MENU INTERATIVO: Seleção do Método de Otimização
    ! ==============================================================================
    print *, "Selecione o metodo de otimizacao para maximizar a funcao f(x,y):"
    print *, "Digite '1' para o metodo do Aclive Maximo (Steepest Ascent)"
    print *, "Digite '2' para o metodo de Fletcher-Reeves (Gradientes Conjugados)"
    read(*, *) metodo

    ! Ponto de partida (chute inicial) fixado para a simulação
    x_0 = -2.0d0
    y_0 = 3.0d0

    ! ==============================================================================
    ! GERENCIAMENTO DE ARQUIVOS DE LOG
    ! ==============================================================================
    ! Para evitar que um método sobrescreva o log do outro quando o programa rodar
    ! duas vezes seguidas, abrimos apenas a unidade correspondente à escolha do usuário.
    if (metodo == 1) then
        open(unit=10, file='output1.dat', status='replace')
    else if (metodo == 2) then
        open(unit=11, file='output2.dat', status='replace')
    end if

    ! Chama o motor principal de otimização, passando o ponto inicial e o método escolhido
    call otimizador_de_passo(x_0, y_0, metodo)

    ! Fecha o arquivo de log que foi utilizado
    if (metodo == 1) then
        close(10)
    else if (metodo == 2) then
        close(11)
    end if

    ! ==============================================================================
    ! GERAÇÃO DA MALHA PARA AS CURVAS DE NÍVEL (function.dat)
    ! ==============================================================================
    open(unit=12, file='function.dat', status='replace')

    ! Utilizamos contadores inteiros (i, j) para evitar erros de arredondamento e avisos
    ! do compilador Fortran ao criar a malha regular.
    do i = 0, 60
        x_grid = -3.0d0 + real(i, 8) * 0.1d0 ! Varredura de X: -3.0 a 3.0

        do j = 0, 50
            y_grid = -1.0d0 + real(j, 8) * 0.1d0 ! Varredura de Y: -1.0 a 4.0

            z_grid = f_objetivo(x_grid, y_grid)  ! Calcula a altura Z da montanha
            write(12, *) x_grid, y_grid, z_grid
        end do
        ! A linha vazia '(A)' é crucial para o software de plotagem entender a quebra de linha da matriz 3D
        write(12, '(A)') ""
    end do

    close(12)

contains

    ! ==============================================================================
    ! SUBROTINA PRINCIPAL: Otimizador Multidimensional
    ! ==============================================================================
    subroutine otimizador_de_passo(x_inicial, y_inicial, metodo)
        implicit none

        ! Argumentos de entrada
        real(8), intent(in) :: x_inicial, y_inicial
        integer, intent(in) :: metodo

        ! Variáveis de controle de iteração e erro
        integer :: iter, iter_interna
        real(8) :: erro, tol
        real(8) :: x_k, y_k

        ! Variáveis para a interpolação quadrática (Busca unidimensional)
        real(8) :: g0, g1, g2, h_star, g_star
        real(8) :: h_lower, h_mid, h_upper
        real(8) :: x_low, y_low, x_mid, y_mid, x_up, y_up
        real(8) :: num_quad, den_quad

        ! Vetores de gradiente e direção de busca
        real(8) :: x_grad, y_grad
        real(8) :: dx, dy
        real(8) :: dx_ant, dy_ant
        real(8) :: grad_ant_sq, grad_atual_sq
        real(8) :: beta

        ! Inicializações
        iter = 0
        erro = 1.0d0
        tol = 1.0d-4  ! Condição de parada (tolerância para o módulo do gradiente)

        x_k = x_inicial
        y_k = y_inicial

        ! O laço principal continua até chegarmos ao topo da montanha (onde o erro se aproxima de 0)
        do while (erro > tol)

            iter = iter + 1

            ! 1. Calcula as derivadas parciais no ponto atual
            x_grad = f_grad_x(x_k, y_k)
            y_grad = f_grad_y(x_k, y_k)

            ! 2. Calcula o módulo do gradiente ao quadrado e define o erro atual
            grad_atual_sq = x_grad**2 + y_grad**2
            erro = sqrt(grad_atual_sq)

            ! Se atingimos o ponto ótimo plano, abortamos a busca para não dar um passo falso
            if (erro <= tol) exit

            ! ==============================================================================
            ! DEFINIÇÃO DA DIREÇÃO DE BUSCA (dx, dy) COM BASE NO MÉTODO ESCOLHIDO
            ! ==============================================================================
            ! Regra de ouro: A primeira iteração de ambos os métodos é sempre o Aclive Máximo puro
            if (metodo == 1 .or. iter == 1) then
                dx = x_grad
                dy = y_grad
            else
                ! Método 2: Fletcher-Reeves (Gradientes Conjugados)
                ! O parâmetro Beta guarda a "memória" da inércia da iteração anterior
                beta = grad_atual_sq / grad_ant_sq
                dx = x_grad + beta * dx_ant
                dy = y_grad + beta * dy_ant
            end if

            ! Armazena o estado atual para ser usado como histórico na próxima rodada do FR
            dx_ant = dx
            dy_ant = dy
            grad_ant_sq = grad_atual_sq

            ! ==============================================================================
            ! BUSCA LINEAR: Interpolação Quadrática (Ajuste do passo ótimo h*)
            ! ==============================================================================
            ! Chutes iniciais para o tamanho do passo
            h_lower = 0.0d0
            h_mid   = 0.05d0
            h_upper = 0.1d0

            iter_interna = 0

            ! Refina o passo h até o intervalo se estreitar ou atingir um limite de iterações
            do while (abs(h_upper - h_lower) > 1.0d-6 .and. iter_interna < 10)
                iter_interna = iter_interna + 1

                ! Transforma o passo h unidimensional nas coordenadas reais (X, Y) do mapa
                x_low = f_xk(h_lower, x_k, dx)
                y_low = f_yk(h_lower, y_k, dy)
                x_mid = f_xk(h_mid, x_k, dx)
                y_mid = f_yk(h_mid, y_k, dy)
                x_up  = f_xk(h_upper, x_k, dx)
                y_up  = f_yk(h_upper, y_k, dy)

                ! Mede a altitude da montanha em cada um dos 3 pontos do chute
                g0 = f_objetivo(x_low, y_low)
                g1 = f_objetivo(x_mid, y_mid)
                g2 = f_objetivo(x_up, y_up)

                ! Calcula o vértice da parábola que se ajusta a esses 3 pontos
                num_quad = g0*(h_mid**2 - h_upper**2) + g1*(h_upper**2 - h_lower**2) + g2*(h_lower**2 - h_mid**2)
                den_quad = g0*(h_mid - h_upper) + g1*(h_upper - h_lower) + g2*(h_lower - h_mid)

                ! Fallback de segurança: Proteção contra divisão por zero se a parábola degenerar
                if (abs(den_quad) < 1.0d-8) then
                    if (g0 >= g1 .and. g0 >= g2) then
                        h_star = h_lower
                    else if (g1 >= g0 .and. g1 >= g2) then
                        h_star = h_mid
                    else
                        h_star = h_upper
                    end if
                    exit ! Se degenerou, assumimos o melhor ponto conhecido e saímos do refinamento
                else
                    ! Divisão segura para encontrar o pico exato da parábola
                    h_star = 0.5d0 * (num_quad / den_quad)
                end if

                ! Avalia a altura exata no novo ponto ótimo encontrado
                g_star = f_objetivo(f_xk(h_star, x_k, dx), f_yk(h_star, y_k, dy))

                ! Lógica de Bracketing: Re-centraliza o intervalo ao redor do novo pico
                if (h_star > h_mid) then
                    h_lower = h_mid
                    h_mid = h_star
                else
                    h_upper = h_mid
                    h_mid = h_star
                end if

            end do

            ! ==============================================================================
            ! ATUALIZAÇÃO DO PONTO REAL E GRAVAÇÃO DO HISTÓRICO
            ! ==============================================================================
            ! Dá o salto no mapa usando a direção escolhida (dx, dy) multiplicada pelo melhor passo (h_star)
            x_k = x_k + h_star * dx
            y_k = y_k + h_star * dy

            ! Imprime o log na tela
            write(*, "(I4, 6F12.6)") iter, erro, h_star, x_k, y_k, x_grad, y_grad

            ! Grava o log da iteração no arquivo correspondente ao método
            if (metodo == 1) then
                write(10, "(I4, 6F12.6)") iter, erro, h_star, x_k, y_k, x_grad, y_grad
            else if (metodo == 2) then
                write(11, "(I4, 6F12.6)") iter, erro, h_star, x_k, y_k, x_grad, y_grad
            end if

        end do

    end subroutine otimizador_de_passo

    ! ==============================================================================
    ! FUNÇÕES MATEMÁTICAS: O Terreno e suas Derivadas
    ! ==============================================================================
    ! Função objetivo (a montanha que queremos escalar)
    real(8) function f_objetivo(x_func, y_func)
        real(8), intent(in) :: x_func, y_func
        f_objetivo = 2.0d0*x_func*y_func + 2.0d0*x_func - (x_func**2) - 2.0d0*(y_func**2)
    end function

    ! Derivada parcial em X (Componente horizontal do gradiente)
    real(8) function f_grad_x(x_func, y_func)
        real(8), intent(in) :: x_func, y_func
        f_grad_x = 2.0d0*y_func + 2.0d0 - 2.0d0*x_func
    end function

    ! Derivada parcial em Y (Componente vertical do gradiente)
    real(8) function f_grad_y(x_func, y_func)
        real(8), intent(in) :: x_func, y_func
        f_grad_y = 2.0d0*x_func - 4.0d0*y_func
    end function

    ! Funções paramétricas: Transformam o passo (h) em coordenadas reais X e Y
    real(8) function f_xk(h, x_atual, dir_x)
        real(8), intent(in) :: h, x_atual, dir_x
        f_xk = x_atual + h * dir_x
    end function

    real(8) function f_yk(h, y_atual, dir_y)
        real(8), intent(in) :: h, y_atual, dir_y
        f_yk = y_atual + h * dir_y
    end function

end program PPC4
