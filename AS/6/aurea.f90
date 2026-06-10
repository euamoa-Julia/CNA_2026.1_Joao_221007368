program aurea
    implicit none

    real(8) :: ra_opt, p_max, ra1, ra2
    real(8) :: r1, r2, r3
    real(8) :: v
    real(8) :: razao
    real(8) :: ra_lower, ra_upper, length_1
    real(8) :: erro, tol
    integer(4) :: iter

    iter = 0

    tol  = 1.0d-4
    erro = 1.0d0

    r1 = 8.0d0
    r2 = 12.0d0
    r3 = 10.0d0
    v  = 80.0d0

    razao = (sqrt(5.0d0) - 1.0d0) / 2.0d0

    ra_upper = 100.0d0
    ra_lower = 0.0d0

    write(*,*) "--------------------------------------------------------------------------------"
    write(*,*) "                      OTIMIZACAO: METODO DA RAZAO AUREA                         "
    write(*,*) "--------------------------------------------------------------------------------"
    write(*,"(A)") " Iter |    ra_lower    |    ra_upper    |   ra_otimizado |      Erro      "
    write(*,*) "--------------------------------------------------------------------------------"

    do while (abs(erro) .gt. tol)

        iter = iter + 1

        length_1 = (ra_upper - ra_lower) * razao
        ra1 = ra_lower + length_1
        ra2 = ra_upper - length_1

        if (p(ra2) .gt. p(ra1)) then
            ra_upper = ra1
        else
            ra_lower = ra2
        end if

        erro = ra1 - ra2
        ra_opt = (ra_upper + ra_lower) / 2.0d0

        write(*,"(I5, ' | ', F14.8, ' | ', F14.8, ' | ', F14.8, ' | ', F14.8)") iter, ra_lower, ra_upper, ra_opt, erro

    end do

    ra_opt = (ra_upper + ra_lower) / 2.0d0
    p_max = p(ra_opt)

    write(*,*) "--------------------------------------------------------------------------------"
    write(*, '(A, F15.8)') ' Ra otimizado:    ', ra_opt
    write(*, '(A, F15.8)') ' Potencia maxima: ', p_max
    write(*,*) "--------------------------------------------------------------------------------"

contains

    real(8) function p(ra)
        real(8), intent(in) :: ra
        real(8) :: num, den
        num = v * r3 * ra
        den = r1 * (ra + r2 + r3) + r3 * ra + r3 * r2
        p = ((num / den)**2.0d0) / ra
    end function

end program aurea

! ------------------------------------------------------------------------------------------------
! Roteiro Sumário: Este programa implementa o Método da Razão Áurea para otimização
! unidimensional. É um algoritmo de busca por eliminação de regiões que utiliza a proporção
! áurea para encolher continuamente, e de forma simétrica, o intervalo de análise, isolando
! o vértice da função até atingir a precisão desejada sem depender de avaliações de derivadas.
!
! Aplicação: Perfeito para localizar máximos ou mínimos em funções unimodais complexas
! num intervalo fechado conhecido. Neste contexto, é aplicado para descobrir a resistência
! elétrica exata que maximiza a potência dissipada pelo circuito, garantindo convergência
! estável graças à razão constante de redução do intervalo.
!
! Inputs: ra_lower e ra_upper (limites iniciais da resistência), tol (tolerância), v, r1,
! r2, r3 (parâmetros do circuito) e a função p(ra).
!
! Outputs: Histórico iterativo impresso no terminal mostrando o estreitamento dos limites,
! seguido pelos valores finais do ponto otimizado no eixo X (resistência ideal) e a sua
! imagem no eixo Y (potência máxima atingida).
! ------------------------------------------------------------------------------------------------
