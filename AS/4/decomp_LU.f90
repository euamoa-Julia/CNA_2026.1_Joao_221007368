program eliminacao
    implicit none

    real(8) :: A(10, 10)
    real(8) :: U(10, 10)
    real(8) :: L(10, 10)
    real(8) :: B(10)
    real(8) :: X(10)
    real(8) :: Y(10)
    real(8) :: matriz_LU(10, 10)
    real(8) :: verificando(10)
    real(8) :: fator, soma
    integer :: i, k, j, n

    ! --- 1. Inicialização e Entrada de Dados ---
    n = 10

    ! Matriz de coeficientes do sistema (n x n)
    A = reshape([ &
        45.0d0, 12.0d0, 28.0d0,  5.0d0, 33.0d0, 18.0d0, 41.0d0,  7.0d0, 22.0d0, 39.0d0, &
        14.0d0, 50.0d0,  9.0d0, 26.0d0,  4.0d0, 37.0d0, 11.0d0, 48.0d0, 16.0d0, 25.0d0, &
        31.0d0,  8.0d0, 42.0d0, 19.0d0, 47.0d0,  2.0d0, 35.0d0, 13.0d0, 29.0d0,  6.0d0, &
         3.0d0, 24.0d0, 17.0d0, 49.0d0, 30.0d0, 10.0d0, 46.0d0, 21.0d0, 38.0d0,  1.0d0, &
        27.0d0, 15.0d0, 40.0d0, 23.0d0, 44.0d0, 32.0d0,  5.0d0, 36.0d0,  9.0d0, 20.0d0, &
        11.0d0, 39.0d0,  4.0d0, 28.0d0, 16.0d0, 45.0d0, 22.0d0,  1.0d0, 43.0d0, 18.0d0, &
        41.0d0,  2.0d0, 33.0d0, 10.0d0, 29.0d0,  6.0d0, 48.0d0, 14.0d0, 31.0d0, 25.0d0, &
        19.0d0, 46.0d0,  7.0d0, 34.0d0, 12.0d0, 27.0d0,  4.0d0, 42.0d0, 15.0d0, 38.0d0, &
         8.0d0, 21.0d0, 49.0d0,  3.0d0, 37.0d0, 13.0d0, 26.0d0, 30.0d0, 47.0d0, 11.0d0, &
        35.0d0,  1.0d0, 20.0d0, 43.0d0,  9.0d0, 40.0d0, 17.0d0, 24.0d0,  5.0d0, 50.0d0  &
        ], shape(A), order = [2, 1])

    ! Vetor de termos independentes
    B = [11.0d0, 25.0d0, 36.0d0, 89.0d0, 45.0d0, 12.0d0, 32.0d0, 42.0d0, 21.0d0, 1.0d0]

    ! Prepara as matrizes L e U para receberem a decomposição
    U = A
    L = 0.0d0

    write(*,*) "----------------------------------------------------------------------------------------------------"
    write(*,*) "Matriz A Original:"
    do i = 1, n
        write(*,"(10F8.2)") A(i,:)
    end do
    write(*,*) "----------------------------------------------------------------------------------------------------"

    ! --- 2. PASSO 1: Fatoração LU (Decomposição) ---
    ! Preenche a diagonal principal da matriz inferior (L) com números 1
    do i = 1, n
        L(i,i) = 1.0d0
    end do

    ! Laço de eliminação: transforma a matriz U em uma triangular superior
    do k = 1, n-1
        ! Percorre as linhas abaixo da diagonal (pivô)
        do i = k+1, n
            ! Calcula o multiplicador para zerar o elemento atual
            fator = U(i,k) / U(k,k)

            ! Armazena o multiplicador na matriz inferior (L)
            L(i,k) = fator

            ! Atualiza o restante da linha da matriz U subtraindo a linha do pivô
            do j = k+1, n
                U(i,j) = U(i,j) - fator * U(k,j)
            end do

            ! Zera explicitamente o elemento abaixo do pivô para evitar resíduos flutuantes
            U(i,k) = 0.0d0
        end do
    end do

    ! --- 3. PASSO 2: Substituição Progressiva (L * Y = B) ---
    ! Resolve o sistema intermediário de cima para baixo
    Y(1) = B(1)
    do i = 2, n
        soma = B(i)
        do j = 1, i-1
            soma = soma - L(i,j) * Y(j)
        end do
        Y(i) = soma
    end do

    ! --- 4. PASSO 3: Substituição Regressiva (U * X = Y) ---
    ! Resolve o sistema final de baixo para cima utilizando o vetor Y recém-calculado
    X(n) = Y(n) / U(n,n)
    do i = n-1, 1, -1
        soma = Y(i)
        do j = i+1, n
            soma = soma - U(i,j) * X(j)
        end do
        ! Isola e encontra a incógnita final da respectiva linha
        X(i) = soma / U(i,i)
    end do

    ! --- 5. Impressão dos Resultados e Provas Reais ---
    write(*,*) "Matriz L (Lower):"
    do i = 1, n
        write(*,"(10F8.2)") L(i,:)
    end do
    write(*,*) "----------------------------------------------------------------------------------------------------"

    write(*,*) "Matriz U (Upper):"
    do i = 1, n
        write(*,"(10F8.2)") U(i,:)
    end do
    write(*,*) "----------------------------------------------------------------------------------------------------"

    write(*, "(A, 10F8.4)") "Vetor Solução X: ", X
    write(*,*) "----------------------------------------------------------------------------------------------------"

    ! Validação 1: Reconstrução da matriz original multiplicando as faturadas
    matriz_LU = matmul(L, U)
    write(*,*) "Prova Real 1 (Verificando se L * U = A):"
    write(*,*) "Se os valores abaixo baterem com a Matriz A original, a decomposição está perfeita."
    do i = 1, n
        write(*,"(10F8.2)") matriz_LU(i,:)
    end do
    write(*,*) "----------------------------------------------------------------------------------------------------"

    ! Validação 2: Aplicação do vetor X na matriz original para retornar as forças (B)
    verificando = matmul(A, X)
    write (*,*) "Prova Real 2 (Verificando se A * X = B):"
    write(*,*) "Valores esperados de B: 11.0, 25.0, 36.0, 89.0, 45.0, 12.0, 32.0, 42.0, 21.0, 1.0"
    write(*, "(A, 10F8.4)") "Valores obtidos:      ", verificando

end program eliminacao

! ------------------------------------------------------------------------------------------------
! Roteiro Sumário: Este programa resolve sistemas de equações lineares de ordem n através da
! Fatoração LU, decompondo a matriz de coeficientes original (A) em uma matriz triangular
! inferior (L) e uma superior (U). Do ponto de vista da engenharia, esta topologia é
! computacionalmente superior à Eliminação de Gauss clássica, pois desvincula a matriz
! principal do vetor de forças (B), permitindo que o sistema seja resolvido rapidamente
! para múltiplas condições de contorno através de simples substituições progressivas e
! regressivas intermediadas por um vetor auxiliar (Y).
!
! Aplicação: Ideal para resolver sistemas lineares (AX = B), especialmente quando a mesma
! matriz de coeficientes (A) precisa ser avaliada para múltiplos vetores de resultados (B),
! economizando processamento. Também útil para cálculos eficientes de determinantes.
!
! Inputs: n (ordem do sistema e das matrizes), A (matriz quadrada de coeficientes do sistema)
! e B (vetor de termos independentes).
!
! Outputs: Impressão da matriz A espelhada, as matrizes fatoradas L e U, o vetor solução
! final X, e provas reais de validação (L*U = A e A*X = B) diretamente no terminal.
! ------------------------------------------------------------------------------------------------

! ------------------------------------------------------------------------------------------------
! COMO ADAPTAR ESTE CÓDIGO PARA OUTROS SISTEMAS LINEARES:
! Para utilizar este método em matrizes ou problemas diferentes, siga os passos abaixo:
!
! 1. Alterar as Dimensões do Sistema: O código está configurado para uma ordem
!    n = 10. Para um sistema diferente, atualize o valor da variável 'n' e redimensione
!    todas as matrizes (A, L, U, matriz_LU) para (n, n) e os vetores (B, X, Y,
!    verificando) para (n) no bloco de declarações.
!
! 2. Atualizar os Dados de Entrada: Substitua os valores fixos da matriz 'A' e
!    do vetor 'B' pelos do seu novo sistema. Para matrizes grandes, recomenda-se
!    apagar a atribuição manual (reshape) e ler os dados de um arquivo .txt ou .csv.
!
! 3. Ajustar a Formatação de Saída: Os comandos de impressão estão configurados
!    para 10 colunas (ex: "(10F8.2)"). Se você alterar a dimensão 'n', lembre-se
!    de alterar esse número no 'write' para que a matriz seja impressa corretamente
!    no terminal (ex: para n=5, use "(5F8.2)").
!
! ATENÇÃO (Limitação Matemática): Esta implementação utiliza a Fatoração LU
! "Ingênua" (sem pivotamento parcial de linhas). Isso significa que o algoritmo
! será interrompido por divisão por zero se encontrar um valor nulo na diagonal
! principal (U(k,k)). Este código é 100% seguro apenas para matrizes estritamente
! diagonais dominantes ou definidas positivas.
! ------------------------------------------------------------------------------------------------
