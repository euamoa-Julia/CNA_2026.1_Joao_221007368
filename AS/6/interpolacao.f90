program quadratica
    implicit none

    real(8) :: v, r1, r2, r3
    real(8) :: p_0, p_1, p_2, p_3
    real(8) :: ra_0, ra_1, ra_2, ra_3
    real(8) :: num_quad, den_quad
    real(8) :: erro, tol
    real(8) :: ra_lower, ra_mid, ra_upper
    real(8) :: ra_anterior

    erro = 1.0d0
    tol = 1.0d-4

    v = 80.0d0
    r1 = 8.0d0
    r2 = 12.0d0
    r3 = 10.0d0

    ra_0 = 10.0d0
    ra_1 = 39.0d0
    ra_2 = 91.0d0

    if (ra_0 .eq. 0.0d0) then
        ra_0 = ra_0 + 1.0d-4
    end if

    num_quad = p(ra_0)*((ra_1**2.0d0)-(ra_2**2.0d0)) + p(ra_1)*((ra_2**2.0d0)-(ra_0**2.0d0)) + p(ra_2)*((ra_0**2.0d0)-(ra_1**2.0d0))
    den_quad = p(ra_0)*(ra_1-ra_2)*2.0d0 + p(ra_1)*(ra_2-ra_0)*2.0d0 + p(ra_2)*(ra_0-ra_1)*2.0d0
    ra_3 = num_quad / den_quad

    p_0 = p(ra_0)
    p_1 = p(ra_1)
    p_2 = p(ra_2)
    p_3 = p(ra_3)

    do while (abs(erro) .gt. tol)

        ra_anterior = ra_3

        if (ra_3 .gt. ra_1) then
            if (p_3 .gt. p_1) then
                ra_lower = ra_1
                ra_mid = ra_3
                ra_upper = ra_2
            else
                ra_lower = ra_0
                ra_mid = ra_1
                ra_upper = ra_3
            end if
        else
            if (p_3 .gt. p_1) then
                ra_lower = ra_0
                ra_mid = ra_3
                ra_upper = ra_1
            else
                ra_lower = ra_3
                ra_mid = ra_1
                ra_upper = ra_2
            end if
        end if

        ra_0 = ra_lower
        ra_1 = ra_mid
        ra_2 = ra_upper

        p_0 = p(ra_0)
        p_1 = p(ra_1)
        p_2 = p(ra_2)

        num_quad = p(ra_0)*((ra_1**2.0d0)-(ra_2**2.0d0)) + p(ra_1)*((ra_2**2.0d0)-(ra_0**2.0d0)) + p(ra_2)*((ra_0**2.0d0)-(ra_1**2.0d0))
        den_quad = p(ra_0)*(ra_1-ra_2)*2.0d0 + p(ra_1)*(ra_2-ra_0)*2.0d0 + p(ra_2)*(ra_0-ra_1)*2.0d0
        ra_3 = num_quad / den_quad

        p_3 = p(ra_3)

        erro = (ra_3 - ra_anterior) / ra_3

    end do

    write(*,*) ra_3, p_3

contains

    function p(ra)
        real(8) :: p, ra, num, den
        num = v * r3 * ra
        den = r1 * (ra + r2 + r3) + r3 * ra + r3 * r2
        p = ((num / den)**2.0d0) / ra
    end function

end program quadratica

! ------------------------------------------------------------------------------------------------
! Roteiro Sumário: Este programa implementa o Método da Interpolação Quadrática Sucessiva para
! otimização unidimensional (busca de pontos de máximo ou mínimo). Em vez de buscar as raízes
! de uma função, o algoritmo constrói parábolas iterativas passando por três estimativas
! para encontrar analiticamente a coordenada do vértice, estreitando o intervalo de busca
! até atingir o cume ou o vale da curva.
!
! Aplicação: Excelente para problemas de otimização onde a derivada da função objetivo não
! está disponível ou é muito complexa. No contexto deste código em específico, é utilizado
! para determinar o valor ótimo de resistência (ra) que resulta na máxima transferência de
! potência elétrica (p) em um circuito parametrizado.
!
! Inputs: ra_0, ra_1, ra_2 (chutes iniciais triplos de resistência), parâmetros do circuito
! (v, r1, r2, r3), tol (tolerância do critério de parada) e a função de potência p(ra).
!
! Outputs: Impressão direta no terminal contendo dois valores numéricos finais: o valor
! ótimo encontrado para a variável independente (resistência ótima) e o valor máximo
! avaliado na função objetivo (potência máxima).
! ------------------------------------------------------------------------------------------------
