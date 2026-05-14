program sala5

implicit none

real A(10, 10)
real L(10, 10)
real U(10, 10)
real B(10)
real X(10)
real Y(10)
real matriz_LU(10, 10)
real erro_relativo(10,10)
real verificando(10)
real fator, soma
integer i, k, j, n

n = 10

A = reshape([ &
    45.0, 12.0, 28.0,  5.0, 33.0, 18.0, 41.0,  7.0, 22.0, 39.0, & ! Linha 1
    14.0, 50.0,  9.0, 26.0,  4.0, 37.0, 11.0, 48.0, 16.0, 25.0, & ! Linha 2
    31.0,  8.0, 42.0, 19.0, 47.0,  2.0, 35.0, 13.0, 29.0,  6.0, & ! Linha 3
     3.0, 24.0, 17.0, 49.0, 30.0, 10.0, 46.0, 21.0, 38.0,  1.0, & ! Linha 4
    27.0, 15.0, 40.0, 23.0, 44.0, 32.0,  5.0, 36.0,  9.0, 20.0, & ! Linha 5
    11.0, 39.0,  4.0, 28.0, 16.0, 45.0, 22.0,  1.0, 43.0, 18.0, & ! Linha 6
    41.0,  2.0, 33.0, 10.0, 29.0,  6.0, 48.0, 14.0, 31.0, 25.0, & ! Linha 7
    19.0, 46.0,  7.0, 34.0, 12.0, 27.0,  4.0, 42.0, 15.0, 38.0, & ! Linha 8
     8.0, 21.0, 49.0,  3.0, 37.0, 13.0, 26.0, 30.0, 47.0, 11.0, & ! Linha 9
    35.0,  1.0, 20.0, 43.0,  9.0, 40.0, 17.0, 24.0,  5.0, 50.0  & ! Linha 10
    ], shape(A), order = [2, 1])

B(1) = 11
B(2) = 25
B(3) = 36
B(4) = 89
B(5) = 45
B(6) = 12
B(7) = 32
B(8) = 42
B(9) = 21
B(10) = 1

U = A
L = 0


        write (*,*) "Matriz A:"
        do i = 1, n
                write(*,"(10F8.2)") A(i,:)
        end do        

write(*,*) "----------------------------------------------------------------------------------------------------"

do i = 1, n
    L(i,i) = 1.0 
end do
        
! Eliminação progressiva

do k = 1,n-1					! linha anterior
   do i=k+1,n					! linha atual
      fator=U(i,k)/U(k,k)			! fator
      L(i,k) = fator
      do j=k+1,n				! colunas
        U(i,j) = U(i,j) - fator*U(k,j)	! eliminando [A]
      end do
      U(i,k) = 0.0
   end do
end do

! ALTERAÇÃO DO GEMINI

! --- PASSO INTERMEDIÁRIO: Resolvendo L * Y = B ---
Y(1) = B(1) 
do i = 2, n
   soma = B(i)
   do j = 1, i-1
      soma = soma - L(i,j) * Y(j)
   end do
   Y(i) = soma  ! Não divide por L(i,i) porque a diagonal de L é tudo 1.0
end do

! --- AGORA a sua Substituição Regressiva usa o Y ---
X(n) = Y(n) / U(n,n)
do i = n-1, 1, -1
   soma = Y(i)  ! ANTES ESTAVA B(i), AGORA É Y(i)
   do j = i+1, n
      soma = soma - U(i,j) * X(j)
   end do
   X(i) = soma / U(i,i)
end do

! ALTERAÇÃO DO GEMINI

! Substituição regressiva
    
    write(*, "(A, 10F8.2)") "Vetor X: ", X

write(*,*) "----------------------------------------------------------------------------------------------------"

        write (*,*) "Matriz Lower:"
        do i = 1, n
                write(*,"(10F8.2)") L(i,:)
        end do
                        
write(*,*) "----------------------------------------------------------------------------------------------------"
        
        write(*,*) "Matriz Upper:"                
        do i = 1, n
            write(*,"(10F8.2)") U(i,:)
        end do
            
write(*,*) "----------------------------------------------------------------------------------------------------"

matriz_lU = matmul(L, U)
erro_relativo = matriz_LU/A

    write(*,*) "Verificando se L x U = A:"        
    do i = 1, n
        write(*,"(10F8.2)") erro_relativo(i,:)          
    end do        

write(*,*) "----------------------------------------------------------------------------------------------------"

verificando = matmul(A, X)
    write (*,*) "Verificando se A x X = B:"
    write(*, "(10F8.2)") verificando

end program sala5

