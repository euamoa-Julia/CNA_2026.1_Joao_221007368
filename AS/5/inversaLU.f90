program sala5
    implicit none

    real(8) :: A(10, 10)
    real(8) :: L(10, 10)
    real(8) :: U(10, 10)
    real(8) :: B(10)
    real(8) :: X(10)
    real(8) :: Y(10)
    real(8) :: matriz_LU(10, 10)
    real(8) :: erro_relativo(10,10)
    real(8) :: verificando(10)

    ! --- Novas variáveis para o cálculo da Inversa ---
    real(8) :: inversa(10, 10)
    real(8) :: prova_inversa(10, 10)
    real(8) :: identidade_col(10)
    integer :: col
    ! -------------------------------------------------

    real(8) :: fator, soma
    integer :: i, k, j, n

    n = 10

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

    B = [11.0d0, 25.0d0, 36.0d0, 89.0d0, 45.0d0, 12.0d0, 32.0d0, 42.0d0, 21.0d0, 1.0d0]

    U = A
    L = 0.0d0

    write (*,*) "Matriz A:"
    do i = 1, n
        write(*,"(10F8.2)") A(i,:)
    end do
    write(*,*) "----------------------------------------------------------------------------------------------------"

    ! --- PASSO 1: DECOMPOSIÇÃO LU ---
    do i = 1, n
        L(i,i) = 1.0d0
    end do

    do k = 1, n-1
        do i = k+1, n
            fator = U(i,k) / U(k,k)
            L(i,k) = fator
            do j = k+1, n
                U(i,j) = U(i,j) - fator * U(k,j)
            end do
            U(i,k) = 0.0d0
        end do
    end do

    ! --- PASSO 2: RESOLUÇÃO DO SISTEMA ORIGINAL (AX = B) ---
    Y(1) = B(1)
    do i = 2, n
        soma = B(i)
        do j = 1, i-1
            soma = soma - L(i,j) * Y(j)
        end do
        Y(i) = soma
    end do

    X(n) = Y(n) / U(n,n)
    do i = n-1, 1, -1
        soma = Y(i)
        do j = i+1, n
            soma = soma - U(i,j) * X(j)
        end do
        X(i) = soma / U(i,i)
    end do

    ! Armazena a prova real de AX = B antes de X ser sobrescrito pelo próximo passo
    verificando = matmul(A, X)

    write(*, "(A, 10F8.2)") "Vetor Solucao X: ", X
    write(*,*) "----------------------------------------------------------------------------------------------------"

    ! --- PASSO 3: CÁLCULO DA MATRIZ INVERSA ---
    ! Resolvemos L*U*X = Identidade coluna por coluna
    do col = 1, n
        ! Cria o vetor coluna da matriz identidade para a iteração atual
        identidade_col = 0.0d0
        identidade_col(col) = 1.0d0

        ! Substituição Progressiva (L * Y = Identidade_col)
        Y(1) = identidade_col(1)
        do i = 2, n
            soma = identidade_col(i)
            do j = 1, i-1
                soma = soma - L(i,j) * Y(j)
            end do
            Y(i) = soma
        end do

        ! Substituição Regressiva (U * X = Y)
        X(n) = Y(n) / U(n,n)
        do i = n-1, 1, -1
            soma = Y(i)
            do j = i+1, n
                soma = soma - U(i,j) * X(j)
            end do
            X(i) = soma / U(i,i)
        end do

        ! Armazena a solução (X) como a respectiva coluna da matriz inversa
        inversa(:, col) = X
    end do

    write (*,*) "Matriz Inversa (A^-1):"
    do i = 1, n
        write(*,"(10F8.4)") inversa(i,:)
    end do
    write(*,*) "----------------------------------------------------------------------------------------------------"

    ! --- VALIDAÇÕES E PROVAS REAIS ---
    matriz_LU = matmul(L, U)
    erro_relativo = matriz_LU / A

    write(*,*) "Prova Real 1 (L x U = A) - Valores de razao esperados ~ 1.00:"
    do i = 1, n
        write(*,"(10F8.2)") erro_relativo(i,:)
    end do
    write(*,*) "----------------------------------------------------------------------------------------------------"

    write (*,*) "Prova Real 2 (Verificando se A x X_original = B):"
    write(*, "(10F8.2)") verificando
    write(*,*) "----------------------------------------------------------------------------------------------------"

    prova_inversa = matmul(A, inversa)
    write(*,*) "Prova Real 3 (A x A^-1 = Identidade):"
    write(*,*) "A diagonal principal deve ser ~ 1.00 e os demais elementos ~ 0.00."
    do i = 1, n
        write(*,"(10F8.4)") prova_inversa(i,:)
    end do
    write(*,*) "----------------------------------------------------------------------------------------------------"

end program sala5

