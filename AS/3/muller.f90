program metodo_muller
    implicit none

    real(8) :: x0, x1, x2, x3
    real(8) :: h0, h1, a, b, c
    real(8) :: delta_0, delta_1, discriminante, den
    real(8) :: tol, erro, passo
    integer :: iter, max_iter

    ! --- 1. Inicialização e Geração dos Chutes Iniciais ---
    passo    = 0.5d0
    x0       = 6.5d0
    tol      = 1.0d-6
    max_iter = 100
    erro     = 1.0d0
    iter     = 0

    ! Gera x1 e x2 automaticamente recuando a partir do chute base (x0)
    x1 = x0 - passo
    x2 = x1 - passo

    write(*,*) "--------------------------------------------------------------------------------"
    write(*,*) "                               METODO DE MULLER                                 "
    write(*,*) "--------------------------------------------------------------------------------"
    write(*,"(A, F8.4, F8.4, F8.4)") " Chutes iniciais (x0, x1, x2): ", x0, x1, x2
    write(*,"(A, E8.1)")             " Tolerancia:                   ", tol
    write(*,*) "--------------------------------------------------------------------------------"
    write(*,"(A)") " Iter |        x3         |        f(x3)       |        Erro        "
    write(*,*) "--------------------------------------------------------------------------------"

    ! Prepara o arquivo externo para salvar o resultado final
    open(unit=10, file='raizes.dat', status='replace')

    ! --- 2. LOOP PRINCIPAL: Interpolação Quadrática ---
    do while (erro .gt. tol .and. iter .lt. max_iter)
        iter = iter + 1

        ! Cálculo dos espaçamentos entre os pontos no eixo X
        h0 = x1 - x0
        h1 = x2 - x1

        ! Cálculo das primeiras diferenças divididas (taxas de variação)
        delta_0 = (f(x1) - f(x0)) / h0
        delta_1 = (f(x2) - f(x1)) / h1

        ! --- 3. Coeficientes da Parábola ---
        ! Define os coeficientes a, b e c para a equação quadrática P(x)
        ! que passa exatamente pelos três pontos atuais.
        a = (delta_1 - delta_0) / (h1 + h0)
        b = (a * h1) + delta_1
        c = f(x2)

        ! --- 4. Discriminante e Trava de Segurança ---
        discriminante = (b**2) - (4.0d0 * a * c)

        ! Como as variáveis são 'real(8)', raízes negativas no discriminante
        ! causariam erro de compilação. O código aborta se a parábola não cruzar o eixo X.
        if (discriminante .lt. 0.0d0) then
            write(*,*) "Aviso: Discriminante negativo. Raízes complexas detectadas."
            exit
        end if

        ! --- 5. Escolha do Denominador e Cálculo da Raiz ---
        ! O Método de Muller soma ou subtrai a raiz do discriminante no denominador
        ! O objetivo é escolher o sinal que resulta no MAIOR denominador em valor absoluto.
        ! Isso evita erros de cancelamento catastrófico e garante que o algoritmo
        ! encontre a raiz mais próxima da estimativa atual (x2).
        if (b .ge. 0.0d0) then
            den = b + sqrt(discriminante)
        else
            den = b - sqrt(discriminante)
        end if

        ! Calcula o novo ponto estimado (x3) onde a parábola cruza o eixo
        x3 = x2 - ((2.0d0 * c) / den)

        ! O erro é a diferença entre a nova estimativa e a mais recente
        erro = abs(x3 - x2)

        write(*,"(I5, ' | ', F15.8, ' | ', E17.8, ' | ', E17.8)") iter, x3, f(x3), erro

        ! --- 6. Atualização da Janela (Deslizamento) ---
        ! Descarta o ponto mais antigo (x0) e insere a nova estimativa (x3)
        x0 = x1
        x1 = x2
        x2 = x3

    end do

    ! --- 7. Finalização e Exportação de Resultados ---
    write(*,*) "--------------------------------------------------------------------------------"
    write(*, '(A, F15.8)') ' Raiz final encontrada: ', x3
    write(*, '(A, I5)')    ' Iteracoes realizadas:  ', iter
    write(*,*) "--------------------------------------------------------------------------------"

    write(10, '(A, F12.8)') 'Raiz encontrada:', x3
    close(10)

contains

    ! Função Analisada
    real(8) function f(x_func)
        real(8), intent(in) :: x_func
        f = (x_func**3) - (13.0d0 * x_func) - 12.0d0
    end function

end program metodo_muller

! ------------------------------------------------------------------------------------------------
! Roteiro Sumário: Este programa implementa o Método de Muller, uma técnica avançada de busca
! de raízes que substitui as aproximações lineares por uma interpolação quadrática (parábolas)
! baseada em três estimativas iniciais e suas diferenças divididas. A grande vantagem
! computacional deste método é a sua convergência superlinear e a capacidade de encontrar raízes
! complexas de forma nativa, tudo isso sem exigir o cálculo de derivadas analíticas da função principal.
!
! Aplicação: Ideal para funções sem derivadas analíticas fáceis, oferecendo convergência
! superlinear (mais rápida que a Secante) e estabilidade em zonas de alta curvatura. Esta
! implementação em específico restringe-se a raízes reais devido à trava do discriminante.
!
! Inputs: x0 (chute inicial), passo (distância para gerar x1 e x2), tol (tolerância de erro),
! max_iter (limite de segurança) e f(x) (função analisada).
!
! Outputs: Tabela de convergência no terminal (iteração, estimativa, f(x), erro) e a gravação
! do valor da raiz no arquivo 'raizes.dat'.
! ------------------------------------------------------------------------------------------------

! ------------------------------------------------------------------------------------------------
! COMO ADAPTAR ESTE CÓDIGO PARA OUTRAS FUNÇÕES OU PROBLEMAS:
! Para utilizar este método na busca de raízes de outras equações, siga os passos abaixo:
!
! 1. Alterar a Função Analisada: Modifique a função 'f' (localizada no final, dentro
!    do bloco 'contains'). Substitua a equação polinomial atual pela função matemática
!    do seu novo problema.
!
! 2. Ajustar o Ponto de Partida e o Passo: O Método de Muller exige três estimativas
!    iniciais. Este código facilita isso gerando os três pontos automaticamente. Você
!    só precisa alterar o valor de 'x0' (ponto base) e 'passo' (distância entre os
!    pontos) no início do programa para focar na região onde suspeita haver uma raiz.
!
! 3. Gerenciar o Arquivo de Saída: O programa exporta a raiz encontrada para o arquivo
!    'raizes.dat'. Se for executar o código para diferentes problemas sequencialmente,
!    lembre-se de alterar o nome do arquivo na instrução 'open' para não sobrescrever
!    os seus resultados anteriores.
!
! ATENÇÃO (Limitação sobre Raízes Complexas): O Método de Muller original é famoso
! por encontrar raízes complexas de forma nativa. No entanto, ESTA implementação foi
! projetada apenas para raízes reais (utilizando variáveis 'real(8)') e possui um
! bloqueio de segurança ('if discriminante .lt. 0.0d0') que aborta a execução. Para
! buscar raízes complexas, será necessário reescrever as declarações de variáveis
! para 'complex(8)' e remover essa trava de interrupção.
! ------------------------------------------------------------------------------------------------
