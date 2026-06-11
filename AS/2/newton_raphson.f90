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

! ------------------------------------------------------------------------------------------------
! Roteiro Sumário: Este programa implementa o Método de Newton-Raphson para otimização e busca
! de raízes de funções não-lineares, utilizando a derivada analítica para traçar retas tangentes
! que convergem iterativamente para o eixo das abscissas. O código destaca-se pela inclusão do
! Método de Newton Modificado, que emprega a segunda derivada para estabilizar o denominador,
! garantindo uma convergência rápida e segura mesmo em regiões matemáticas críticas que contêm
! raízes múltiplas ou repetidas.
!
! Aplicação: Excelente para funções onde a derivada e a segunda derivada analíticas são
! conhecidas e fáceis de calcular. Oferece convergência muito rápida (quadrática) se o chute
! inicial for bom. Especialmente útil para raízes múltiplas/repetidas utilizando a versão
! modificada f_NRmod presente no código.
!
! Inputs: x_atual (chute inicial), tol (tolerância de erro), max_iter (limite de iterações),
! f(x) (função analisada), f_derivada1(x) (1ª derivada) e f_derivada2(x) (2ª derivada).
!
! Outputs: Tabela detalhada de convergência no terminal (Iteração, estimativa, f(x), erro),
! além da gravação da raiz e do total de iterações no arquivo 'raizes.dat'.
! ------------------------------------------------------------------------------------------------

! ------------------------------------------------------------------------------------------------
! COMO ADAPTAR ESTE CÓDIGO PARA OUTRAS FUNÇÕES OU PROBLEMAS:
! Para utilizar este método na busca de raízes de outras equações, siga os passos abaixo:
!
! 1. Alterar a Função e suas Derivadas: Modifique a função 'f' e as suas derivadas
!    analíticas 'f_derivada1' e 'f_derivada2' (localizadas no final do código, no
!    bloco 'contains'). Certifique-se de que as derivadas matemáticas correspondam
!    exatamente à função principal para que a convergência não seja destruída.
!
! 2. Alternar entre Newton Padrão e Modificado: O loop iterativo atual chama a
!    função 'f_NR' (Newton clássico). Se você estiver lidando com raízes múltiplas
!    (onde a curva tangencia o eixo X em vez de cruzá-lo diretamente), altere a linha
!    'x_novo = f_NR(x_atual)' para 'x_novo = f_NRmod(x_atual)'.
!
! 3. Ajustar o Chute Inicial: A convergência do método de Newton depende fortemente
!    de um bom ponto de partida. Altere o valor de 'x_atual' no início do programa
!    para uma estimativa razoavelmente próxima da raiz desejada.
!
! 4. Gerenciar o Arquivo de Saída: O programa exporta o resultado para 'raizes.dat'.
!    Se for executar o código para diferentes problemas, lembre-se de alterar o nome
!    do arquivo na instrução 'open' para não sobrescrever dados de rodadas anteriores.
!
! ATENÇÃO (Limitações Matemáticas): Este método falhará por divisão por zero se a
! primeira derivada for nula no ponto avaliado (tangente horizontal). Além disso,
! chutes iniciais ruins podem levar o algoritmo a divergir para o infinito ou entrar
! em ciclos oscilatórios. Sempre avalie o comportamento da função antes de rodar.
! ------------------------------------------------------------------------------------------------
