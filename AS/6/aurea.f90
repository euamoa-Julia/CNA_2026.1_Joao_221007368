program aurea
    implicit none

    real(8) :: ra_opt, p_max, ra1, ra2
    real(8) :: r1, r2, r3
    real(8) :: v
    real(8) :: razao
    real(8) :: ra_lower, ra_upper, length_1
    real(8) :: erro, tol
    integer(4) :: iter

    ! --- 1. Inicialização e Parâmetros ---
    iter = 0

    tol  = 1.0d-4
    erro = 1.0d0

    ! Constantes do circuito elétrico
    r1 = 8.0d0
    r2 = 12.0d0
    r3 = 10.0d0
    v  = 80.0d0

    ! Definição da constante da Razão Áurea (~0.61803)
    razao = (sqrt(5.0d0) - 1.0d0) / 2.0d0

    ! Limites iniciais da janela de busca
    ra_upper = 100.0d0
    ra_lower = 0.0d0

    write(*,*) "--------------------------------------------------------------------------------"
    write(*,*) "                      OTIMIZACAO: METODO DA RAZAO AUREA                         "
    write(*,*) "--------------------------------------------------------------------------------"
    write(*,"(A)") " Iter |    ra_lower    |    ra_upper    |   ra_otimizado |      Erro      "
    write(*,*) "--------------------------------------------------------------------------------"

    ! --- 2. LOOP PRINCIPAL: Estreitamento do Intervalo ---
    do while (abs(erro) .gt. tol)

        iter = iter + 1

        ! Calcula a distância proporcional à razão áurea a partir das bordas
        length_1 = (ra_upper - ra_lower) * razao

        ! Gera os dois pontos de teste internos de forma simétrica
        ! ra1 fica mais à direita, ra2 fica mais à esquerda dentro do intervalo
        ra1 = ra_lower + length_1
        ra2 = ra_upper - length_1

        ! --- 3. Lógica de Redução da Janela de Busca ---
        ! Avalia em qual subintervalo está o pico da função (Maximização).
        ! O algoritmo descarta a região que definitivamente não contém o máximo.
        if (p(ra2) .gt. p(ra1)) then
            ! O pico está mais próximo de ra2 (à esquerda), então descartamos
            ! a região à direita de ra1, puxando o limite superior para cá.
            ra_upper = ra1
        else
            ! O pico está mais próximo de ra1 (à direita), então descartamos
            ! a região à esquerda de ra2, puxando o limite inferior para cá.
            ra_lower = ra2
        end if

        ! --- 4. Atualização de Erro e Ponto Ótimo ---
        ! O erro é a distância entre os dois pontos internos atuais
        erro = ra1 - ra2

        ! O ponto ótimo estimado é o ponto médio da nova janela encolhida
        ra_opt = (ra_upper + ra_lower) / 2.0d0

        write(*,"(I5, ' | ', F14.8, ' | ', F14.8, ' | ', F14.8, ' | ', F14.8)") iter, ra_lower, ra_upper, ra_opt, erro

    end do

    ! --- 5. Resultados Finais ---
    ! Refina o valor ótimo final e calcula a sua respectiva potência máxima
    ra_opt = (ra_upper + ra_lower) / 2.0d0
    p_max = p(ra_opt)

    write(*,*) "--------------------------------------------------------------------------------"
    write(*, '(A, F15.8)') ' Ra otimizado:    ', ra_opt
    write(*, '(A, F15.8)') ' Potencia maxima: ', p_max
    write(*,*) "--------------------------------------------------------------------------------"

contains

    ! Função Objetivo: Potência em função da resistência
    real(8) function p(ra)
        real(8), intent(in) :: ra
        real(8) :: num, den
        num = v * r3 * ra
        den = r1 * (ra + r2 + r3) + r3 * ra + r3 * r2
        p = ((num / den)**2.0d0) / ra
    end function

end program aurea

