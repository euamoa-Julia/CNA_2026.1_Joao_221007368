program quadratica
    implicit none

    real(8) :: v, r1, r2, r3
    real(8) :: p_0, p_1, p_2, p_3
    real(8) :: ra_0, ra_1, ra_2, ra_3
    real(8) :: num_quad, den_quad
    real(8) :: erro, tol
    real(8) :: ra_lower, ra_mid, ra_upper
    real(8) :: ra_anterior

    ! --- 1. Inicialização de Variáveis e Parâmetros ---
    erro = 1.0d0
    tol = 1.0d-4

    ! Parâmetros do circuito
    v = 80.0d0
    r1 = 8.0d0
    r2 = 12.0d0
    r3 = 10.0d0

    ! Chutes iniciais (resistências)
    ra_0 = 10.0d0
    ra_1 = 39.0d0
    ra_2 = 91.0d0

    ! Prevenção contra divisão por zero no primeiro ponto
    if (ra_0 .eq. 0.0d0) then
        ra_0 = ra_0 + 1.0d-4
    end if

    ! --- 2. Cálculo do Primeiro Vértice da Parábola ---
    ! As equações abaixo encontram analiticamente a coordenada (ra_3) do vértice
    ! da parábola que passa pelas três estimativas iniciais (ra_0, ra_1, ra_2).
    num_quad = p(ra_0)*((ra_1**2.0d0)-(ra_2**2.0d0)) + p(ra_1)*((ra_2**2.0d0)-(ra_0**2.0d0)) + p(ra_2)*((ra_0**2.0d0)-(ra_1**2.0d0))
    den_quad = p(ra_0)*(ra_1-ra_2)*2.0d0 + p(ra_1)*(ra_2-ra_0)*2.0d0 + p(ra_2)*(ra_0-ra_1)*2.0d0
    ra_3 = num_quad / den_quad

    ! Avaliação da função nos pontos iniciais e no vértice encontrado
    p_0 = p(ra_0)
    p_1 = p(ra_1)
    p_2 = p(ra_2)
    p_3 = p(ra_3)

    ! --- 3. LOOP PRINCIPAL: Refinamento Iterativo ---
    do while (abs(erro) .gt. tol)

        ra_anterior = ra_3

        ! --- Lógica de Redução do Intervalo ---
        ! Compara a posição do novo vértice (ra_3) com o ponto do meio (ra_1) e
        ! verifica qual deles resultou na maior potência (p_3 ou p_1).
        ! O objetivo é reter os três pontos que melhor "cercam" o pico da função.
        if (ra_3 .gt. ra_1) then
            if (p_3 .gt. p_1) then
                ! O pico está entre ra_1 e ra_2
                ra_lower = ra_1
                ra_mid = ra_3
                ra_upper = ra_2
            else
                ! O pico está entre ra_0 e ra_3
                ra_lower = ra_0
                ra_mid = ra_1
                ra_upper = ra_3
            end if
        else
            if (p_3 .gt. p_1) then
                ! O pico está entre ra_0 e ra_1
                ra_lower = ra_0
                ra_mid = ra_3
                ra_upper = ra_1
            else
                ! O pico está entre ra_3 e ra_2
                ra_lower = ra_3
                ra_mid = ra_1
                ra_upper = ra_2
            end if
        end if

        ! Atualiza os três pontos principais para a próxima iteração
        ra_0 = ra_lower
        ra_1 = ra_mid
        ra_2 = ra_upper

        ! Recalcula as imagens para os novos pontos
        p_0 = p(ra_0)
        p_1 = p(ra_1)
        p_2 = p(ra_2)

        ! --- 4. Cálculo do Novo Vértice ---
        ! Traça uma nova parábola com o intervalo reduzido
        num_quad = p(ra_0)*((ra_1**2.0d0)-(ra_2**2.0d0)) + p(ra_1)*((ra_2**2.0d0)-(ra_0**2.0d0)) + p(ra_2)*((ra_0**2.0d0)-(ra_1**2.0d0))
        den_quad = p(ra_0)*(ra_1-ra_2)*2.0d0 + p(ra_1)*(ra_2-ra_0)*2.0d0 + p(ra_2)*(ra_0-ra_1)*2.0d0
        ra_3 = num_quad / den_quad

        p_3 = p(ra_3)

        ! Critério de parada: Erro relativo aproximado
        erro = (ra_3 - ra_anterior) / ra_3

    end do

    ! --- 5. Resultados Finais ---
    write(*,*) ra_3, p_3

contains

    ! Função Objetivo: Potência em função da resistência
    function p(ra)
        real(8) :: p, ra, num, den
        num = v * r3 * ra
        den = r1 * (ra + r2 + r3) + r3 * ra + r3 * r2
        p = ((num / den)**2.0d0) / ra
    end function

end program quadratica

