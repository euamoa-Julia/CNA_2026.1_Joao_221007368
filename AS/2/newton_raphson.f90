program Newton_Raphson
    implicit none

    real(8) :: x_atual, x_novo
    real(8) :: fx, erro, tol
    integer :: iter, max_iter
    real(8) :: raizes(10)

    ! --- 1. Inicialização de Variáveis e Parâmetros ---
    x_atual  = 2.5d0    ! Chute inicial (ponto de partida no eixo X)
    tol      = 1.0d-6   ! Tolerância para o critério de parada
    max_iter = 100      ! Limite de segurança para evitar loops infinitos
    iter     = 0        ! Contador de iterações
    erro     = 1.0d0    ! Erro inicializado propositalmente alto para entrar no laço

    write(*,*) "--------------------------------------------------------------------------------"
    write(*,*) "                       METODO DE NEWTON-RAPHSON                                 "
    write(*,*) "--------------------------------------------------------------------------------"
    write(*,"(A, F8.4)") " Chute inicial: ", x_atual
    write(*,"(A, E8.1)") " Tolerancia:    ", tol
    write(*,*) "--------------------------------------------------------------------------------"
    write(*,"(A)") " Iter |      x_atual      |        f(x)       |        Erro        "
    write(*,*) "--------------------------------------------------------------------------------"

    ! --- 2. Preparação de Exportação de Dados ---
    ! Abre (ou cria) o arquivo para registrar a raiz encontrada ao final
    open (unit=10, file='raizes.dat', status='replace')

    ! --- 3. LOOP PRINCIPAL: Iterações de Newton-Raphson ---
    do while (erro .gt. tol .and. iter .lt. max_iter)
        iter = iter + 1

        ! Avalia a imagem da função no ponto atual
        fx = f(x_atual)

        ! --- 4. Cálculo do Próximo Ponto (Projeção da Tangente) ---
        ! A função f_NR traça uma reta tangente à curva no ponto (x_atual, fx).
        ! O 'x_novo' é exatamente a coordenada onde essa reta tangente cruza o eixo X (y=0).
        x_novo = f_NR(x_atual)

        ! --- 5. Atualização e Critério de Parada ---
        ! O erro é medido pela distância (em valor absoluto) entre o novo ponto e o anterior
        erro = abs(x_novo - x_atual)

        write(*,"(I5, ' | ', F15.8, ' | ', E17.8, ' | ', E17.8)") iter, x_atual, fx, erro

        ! O ponto recém-calculado torna-se o ponto de partida para a próxima iteração
        x_atual = x_novo
    end do

    ! --- 6. Resultados Finais e Gravação ---
    write(*,*) "--------------------------------------------------------------------------------"

    ! Armazena a resposta final na primeira posição do vetor de raízes
    raizes(1) = x_atual

    write(*, "(A, F12.8)") " Raiz encontrada:      ", raizes(1)
    write(*, "(A, I5)")    " Iteracoes realizadas: ", iter
    write(*,*) "--------------------------------------------------------------------------------"

    ! Escreve os resultados no arquivo externo e o fecha com segurança
    write(10, '(A, F12.8)') 'Raiz encontrada:', raizes(1)
    write(10, '(A, I4)') 'Iteracoes realizadas:', iter
    close(10)

contains

    ! --- 7. Funções Analíticas e Fórmulas de Newton ---

    ! Função Principal Analisada
    real(8) function f(x_func)
        real(8), intent(in) :: x_func
        f = (x_func**3) - (6.0d0 * x_func**2) + (9.0d0 * x_func) - 4.0d0
    end function

    ! Primeira Derivada (Declive da reta tangente)
    real(8) function f_derivada1(x_func)
        real(8), intent(in) :: x_func
        f_derivada1 = (3.0d0 * x_func**2) - (12.0d0 * x_func) + 9.0d0
    end function

    ! Segunda Derivada (Taxa de variação do declive - concavidade)
    real(8) function f_derivada2(x_func)
        real(8), intent(in) :: x_func
        f_derivada2 = (6.0d0 * x_func) - 12.0d0
    end function

    ! Fórmula Clássica de Newton-Raphson: x_novo = x - f(x)/f'(x)
    ! Utiliza apenas a primeira derivada. Excelente para raízes simples.
    real(8) function f_NR(x_func)
        real(8), intent(in) :: x_func
        f_NR = x_func - (f(x_func) / f_derivada1(x_func))
    end function

    ! Fórmula Modificada de Newton-Raphson
    ! Incorpora a segunda derivada para estabilizar o cálculo em casos de raízes
    ! múltiplas (onde f'(x) se aproxima de zero, o que faria o método clássico explodir).
    real(8) function f_NRmod(x_func)
        real(8), intent(in) :: x_func
        f_NRmod = x_func - ( (f(x_func) * f_derivada1(x_func)) / &
                  (f_derivada1(x_func)**2 - (f(x_func) * f_derivada2(x_func))) )
    end function

end program Newton_Raphson

