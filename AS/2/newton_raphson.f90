program Newton_Raphson
    implicit none

    real(8) :: x_atual, x_novo
    real(8) :: fx, erro, tol
    integer :: iter, max_iter
    real(8) :: raizes(10)

    x_atual  = 2.5d0
    tol      = 1.0d-6
    max_iter = 100
    iter     = 0
    erro     = 1.0d0

    write(*,*) "--------------------------------------------------------------------------------"
    write(*,*) "                       METODO DE NEWTON-RAPHSON                                 "
    write(*,*) "--------------------------------------------------------------------------------"
    write(*,"(A, F8.4)") " Chute inicial: ", x_atual
    write(*,"(A, E8.1)") " Tolerancia:    ", tol
    write(*,*) "--------------------------------------------------------------------------------"
    write(*,"(A)") " Iter |      x_atual      |        f(x)       |        Erro        "
    write(*,*) "--------------------------------------------------------------------------------"

    open (unit=10, file='raizes.dat', status='replace')

    do while (erro .gt. tol .and. iter .lt. max_iter)
        iter = iter + 1

        fx = f(x_atual)

        x_novo = f_NR(x_atual)

        erro = abs(x_novo - x_atual)

        write(*,"(I5, ' | ', F15.8, ' | ', E17.8, ' | ', E17.8)") iter, x_atual, fx, erro

        x_atual = x_novo
    end do

    write(*,*) "--------------------------------------------------------------------------------"

    raizes(1) = x_atual

    write(*, "(A, F12.8)") " Raiz encontrada:      ", raizes(1)
    write(*, "(A, I5)")    " Iteracoes realizadas: ", iter
    write(*,*) "--------------------------------------------------------------------------------"

    write(10, '(A, F12.8)') 'Raiz encontrada:', raizes(1)
    write(10, '(A, I4)') 'Iteracoes realizadas:', iter
    close(10)

contains

    real(8) function f(x_func)
        real(8), intent(in) :: x_func
        f = (x_func**3) - (6.0d0 * x_func**2) + (9.0d0 * x_func) - 4.0d0
    end function

    real(8) function f_derivada1(x_func)
        real(8), intent(in) :: x_func
        f_derivada1 = (3.0d0 * x_func**2) - (12.0d0 * x_func) + 9.0d0
    end function

    real(8) function f_derivada2(x_func)
        real(8), intent(in) :: x_func
        f_derivada2 = (6.0d0 * x_func) - 12.0d0
    end function

    real(8) function f_NR(x_func)
        real(8), intent(in) :: x_func
        f_NR = x_func - (f(x_func) / f_derivada1(x_func))
    end function

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
