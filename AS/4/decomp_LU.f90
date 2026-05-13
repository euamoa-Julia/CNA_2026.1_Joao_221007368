program eliminacao
    implicit none

    ! --- Variáveis em Dupla Precisão (64 bits) ---
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

    n = 10

    ! --- Definição das Matrizes ---
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

    write(*,*) "----------------------------------------------------------------------------------------------------"
    write(*,*) "Matriz A Original:"
    do i = 1, n
        write(*,"(10F8.2)") A(i,:)
    end do
    write(*,*) "----------------------------------------------------------------------------------------------------"

    ! --- Decomposição LU ---
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

    ! --- 1. Substituição Progressiva (L * Y = B) ---
    Y(1) = B(1)
    do i = 2, n
        soma = B(i)
        do j = 1, i-1
            soma = soma - L(i,j) * Y(j)
        end do
        Y(i) = soma
    end do

    ! --- 2. Substituição Regressiva (U * X = Y) ---
    X(n) = Y(n) / U(n,n)
    do i = n-1, 1, -1
        soma = Y(i)
        do j = i+1, n
            soma = soma - U(i,j) * X(j)
        end do
        X(i) = soma / U(i,i)
    end do

    ! --- Resultados e Verificações ---
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

    matriz_LU = matmul(L, U)
    write(*,*) "Prova Real 1 (Verificando se L * U = A):"
    write(*,*) "Se os valores abaixo baterem com a Matriz A original, a decomposição está perfeita."
    do i = 1, n
        write(*,"(10F8.2)") matriz_LU(i,:)
    end do
    write(*,*) "----------------------------------------------------------------------------------------------------"

    verificando = matmul(A, X)
    write (*,*) "Prova Real 2 (Verificando se A * X = B):"
    write(*,*) "Valores esperados de B: 11.0, 25.0, 36.0, 89.0, 45.0, 12.0, 32.0, 42.0, 21.0, 1.0"
    write(*, "(A, 10F8.4)") "Valores obtidos:      ", verificando

end program eliminacao

! Roteiro Sumário: Este programa resolve sistemas de equações lineares de ordem n através da Fatoração LU, decompondo a matriz de coeficientes original (A) em uma matriz triangular inferior (L) e uma superior (U). Do ponto de vista da engenharia, esta topologia é computacionalmente superior à Eliminação de Gauss clássica, pois desvincula a matriz principal do vetor de forças (B), permitindo que o sistema seja resolvido rapidamente para múltiplas condições de contorno através de simples substituições progressivas e regressivas intermediadas por um vetor auxiliar (Y).
