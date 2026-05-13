! Damos início ao código com o comando "program", Em seguida iremos digitar "implicit none" de forma que
! teremos de declarar todas as nossa variáveis, evitando que o código interprete algo de forma errônea
! e também nos dando mais autonomia sobre o código, permitindo total controle sobre as variáveis.

program Bairstow_Completo
    implicit none

    ! ==============================================================================================
    ! DECLARAÇÃO DE VARIÁVEIS (PARTE 1 E PARTE 2)
    ! ==============================================================================================

    ! Aqui iremos definir todas as variáveis que serão utilizadas no código e nos cálculos, como os
    ! determinantes, discriminantes e os valores de (r, s) que serão nossos chutes iniciais. a(0:100)
    ! é a notação de vetores em Fortran, ou seja, estamos permitindo que o vetor a(0:100) armazene
    ! até 100 valores, assim como b e c. Adicionamos a_base para guardar o polinômio original intacto.
    ! Utilizamos real(8) para garantir Dupla Precisão (64 bits) em todos os cálculos, evitando erros de arredondamento.
    real(8) a(0:100), b(0:100), c(0:100), a_base(0:100)
    real(8) r, s, delta_r, delta_s
    real(8) det, det_r, det_s
    real(8) discriminante

    ! As variáveis de erro e tolerâncias serão utilizadas como critérios de parada para alguns de nossos
    ! laços para que não fiquem rodando infinitamente, para mais do que o necessário, ou para valores indesejados.
    real(8) erro_r, erro_s
    real(8) tol

    ! Definiremos o vetor raiz com 2 espaços, assim como o raiz_complexa, em que estaremos preparados
    ! para encontrar tanto as raízes reais quanto as raízes complexas (que aparecerão em pares conjugados).
    real(8) raiz(1:2)
    complex(8) raiz_complexa(1:2)

    integer n, i, n_base
    integer iter

    ! Variáveis extras exclusivas para a geração do Fractal (Malha Cartesiana e Arquivo CSV)
    integer i_grid, j_grid, unit_csv
    integer, parameter :: res = 1000
    real(8) r0_start, s0_start
    real(8), parameter :: r_min = -10.0d0, r_max = 10.0d0
    real(8), parameter :: s_min = -10.0d0, s_max = 10.0d0
    integer total_iter
    real(8) assinatura_cor


    ! ==============================================================================================
    ! SETUP DO POLINÔMIO (Válido para ambas as partes)
    ! ==============================================================================================

    tol = 1.0d-5

    ! Aqui iremos definir o grau e os coeficientes de nossos polinômios. Para cada polinômio que queira
    ! calcular as raízes, basta alterar os valores seguintes de acordo. Os valores abaixo representam
    ! a equação característica de ordem 4 encontrada no sistema do APC2.
    n_base = 4

    a_base(0) = 28.0d0
    a_base(1) = 595.0d0
    a_base(2) = 3193.0d0
    a_base(3) = 440.0d0
    a_base(4) = 6.0d0

    ! ==============================================================================================
    ! PARTE 1: RESOLUÇÃO ÚNICA DO SISTEMA (Impressão no Terminal)
    ! ==============================================================================================

    print *, "========================================================================="
    print *, " PARTE 1: EXTRACAO DE RAIZES DO SISTEMA (CHUTE UNICO)"
    print *, "========================================================================="

    ! Carregamos o polinômio base para as variáveis de trabalho
    n = n_base
    do i = 0, n
        a(i) = a_base(i)
    end do

    do while (n .gt. 2) ! Este é o nosso primeiro laço, e o mais importante, que define que o método deve rodar enquanto o grau do polinômio for maior ou igual a 2.

        iter = 0
        r = 5.0d0
        s = -3.0d0
        erro_r = 1.0d0
        erro_s = 1.0d0

        do while ((abs(erro_r).gt.tol .and. abs(erro_s).gt.tol)) ! Aqui já entramos no método de Bairstow, em que os valores de r e s serão utilizados para o cálculo de seus respectivos deltas. Dentro dessas fórmulas está embutido o método de Newton Raphson, em que retas tangentes serão definidas a partir dos pontos inicias (r e s) e serão feitas até uma raíz ser encontrada. O delta é nada mais do que a distância até a raiz estimada.

            iter = iter + 1

            b(n) = a(n)
            b(n-1) = a(n-1) + r * b(n)
            c(n) = b(n)
            c(n-1) = b(n-1) + r * c(n)

            do i = n - 2, 0, -1
                b(i) = a(i) + r * b(i+1) + s * b(i+2)
            end do

            do i = n - 2, 1, -1
                c(i) = b(i) + r * c(i+1) + s * c(i+2)
            end do

            det = (c(1)*c(3)) - (c(2)*c(2))
            det_r = ((-b(0))*c(3)) - (c(2)*(-b(1)))
            det_s = (c(1)*(-b(1))) - ((-b(0))*c(2))

            delta_r = det_r / det
            delta_s = det_s / det

            r = r + delta_r
            s = s + delta_s

            erro_r = (delta_r / r)
            erro_s = (delta_s / s)

        end do

        discriminante = (r**2) + (4.0d0 * s) ! Após a deflação polinomial ter sido realizada, podemos calcular o nosso delta da fórmula de Bhaskara, ou o discriminante, de forma a definir se nossas raízes após cada deflação serão reais ou complexas, e iremos printar na tela os pares de raízes a cada deflação.

        if (discriminante .ge. 0.0d0) then
            raiz(1) = (r + sqrt(discriminante)) / 2.0d0
            raiz(2) = (r - sqrt(discriminante)) / 2.0d0
            print *, "Raízes: ", raiz(1), "e", raiz(2)
        else
            raiz_complexa(1) = (r + sqrt(dcmplx(discriminante, 0.0d0))) / 2.0d0
            raiz_complexa(2) = (r - sqrt(dcmplx(discriminante, 0.0d0))) / 2.0d0
            print *, "Raízes:", raiz_complexa(1), "e", raiz_complexa(2)
        end if

        n = n - 2

        do i = 0, n ! Nesta região do código estamos apenas rearranjando os coeficientes após cada deflação.
            a(i) = b(i+2)
        end do

    end do

    if (n .eq. 2) then ! Aqui chegamos ao fim do nosso método, em que todas as deflações foram realizadas e chegaremos em um polinômio de grau 1 ou grau 2, em que poderemos finalmente aplicar algumas fórmulas básicas (fórmula de Bhaskara e equação da reta) para encontrar as última raízes. Lembrando que no caso da fórmula de bhaskara, devemos calcular o discriminante novamente de forma a compreender a natureza das nossas raízes e direcionar as mesmas para o método de solução adequado.

        discriminante = (a(1)**2) - (4.0d0 * a(2) * a(0))

        if (discriminante .ge. 0.0d0) then
            raiz(1) = (-a(1) + sqrt(discriminante)) / (2.0d0 * a(2))
            raiz(2) = (-a(1) - sqrt(discriminante)) / (2.0d0 * a(2))
            print *, "Raízes Finais =", raiz(1), " e ", raiz(2)
            print *, "Iterações no último bloco:", iter
        else
            raiz_complexa(1) = (-a(1) + sqrt(dcmplx(discriminante, 0.0d0))) / (2.0d0 * a(2))
            raiz_complexa(2) = (-a(1) - sqrt(dcmplx(discriminante, 0.0d0))) / (2.0d0 * a(2))
            print *, "Raízes Finais =", raiz_complexa(1), " e ", raiz_complexa(2)
            print *, "Iterações no último bloco:", iter
        end if

    end if

    if (n .eq. 1) then ! No caso em que temos um polinômio de grau 1, ou seja, uma reta.
        raiz(1) = -(a(0)/a(1))
        print *, "Raíz Finais = ", raiz(1)
        print *, "Iterações no último bloco:", iter
    end if


    ! ==============================================================================================
    ! PARTE 2: MAPEAMENTO DO FRACTAL (Geração de dados CSV em Malha)
    ! ==============================================================================================

    print *, " "
    print *, "========================================================================="
    print *, " PARTE 2: GERACAO DO FRACTAL DE BAIRSTOW (VARREDURA)"
    print *, "========================================================================="
    print *, "Iniciando o mapeamento cartesianano com limite de 1500 iterações..."

    ! O comando open criará o arquivo fractal_bairstow.csv. Como estamos dentro de dois laços gigantes
    ! (i_grid e j_grid) que simulam os pixels de uma imagem, salvaremos os dados linha por linha.
    open(newunit=unit_csv, file='fractal_bairstow.csv', status='replace')
    write(unit_csv, *) "r0, s0, cor, iteracoes"

    ! Aqui inicia-se a varredura do "terreno" matemático. Para cada coordenada X (i_grid) e Y (j_grid),
    ! o programa estipulará um chute inicial diferente de 'r' e 's'.
    do i_grid = 0, res
        r0_start = r_min + (r_max - r_min) * (real(i_grid, 8) / real(res, 8))

        do j_grid = 0, res
            s0_start = s_min + (s_max - s_min) * (real(j_grid, 8) / real(res, 8))

            ! Como o método de Bairstow destrói o polinômio durante a deflação, precisamos "recarregar"
            ! o polinômio original (a_base) para testar o próximo pixel da imagem com a equação limpa.
            n = n_base
            do i = 0, n
                a(i) = a_base(i)
            end do

            r = r0_start
            s = s0_start
            assinatura_cor = 0.0d0
            total_iter = 0

            ! O Laço Mestre do Bairstow repetido aqui, mas de forma silenciosa (sem os prints)
            ! e com a válvula de segurança de 1500 iterações para lidar com as zonas caóticas do fractal.
            do while (n .gt. 2)
                iter = 0
                erro_r = 1.0d0
                erro_s = 1.0d0

                do while ((abs(erro_r) .gt. tol .and. abs(erro_s) .gt. tol) .and. iter .lt. 1500)
                    iter = iter + 1
                    total_iter = total_iter + 1

                    b(n) = a(n)
                    b(n-1) = a(n-1) + r * b(n)
                    c(n) = b(n)
                    c(n-1) = b(n-1) + r * c(n)

                    do i = n - 2, 0, -1
                        b(i) = a(i) + r * b(i+1) + s * b(i+2)
                    end do

                    do i = n - 2, 1, -1
                        c(i) = b(i) + r * c(i+1) + s * c(i+2)
                    end do

                    det = (c(1)*c(3)) - (c(2)*c(2))

                    ! Válvula de segurança: se o determinante for quase zero (pico da montanha),
                    ! damos um leve empurrão no dardo para evitar erro de divisão por zero.
                    if (abs(det) .lt. 1.0e-12_8) then
                        r = r + 0.1d0
                        s = s + 0.1d0
                        cycle
                    end if

                    det_r = ((-b(0))*c(3)) - (c(2)*(-b(1)))
                    det_s = (c(1)*(-b(1))) - ((-b(0))*c(2))

                    delta_r = det_r / det
                    delta_s = det_s / det

                    r = r + delta_r
                    s = s + delta_s

                    erro_r = (delta_r / r)
                    erro_s = (delta_s / s)
                end do

                if (iter .ge. 1500) then
                    assinatura_cor = -1.0d0
                    exit
                end if

                ! Se o algoritmo achou raízes, começamos a somar os valores absolutos para criar
                ! uma "assinatura" matemática para colorir aquele pixel no Gnuplot.
                assinatura_cor = assinatura_cor + abs(r) + abs(s)

                n = n - 2
                do i = 0, n
                    a(i) = b(i+2)
                end do
            end do

            ! Se divergiu no caminho, aborta a matemática e salva como área caótica (-1)
            if (assinatura_cor .eq. -1.0d0) then
                write(unit_csv, '(4(F12.6, ","))') r0_start, s0_start, assinatura_cor, real(total_iter, 8)
                cycle
            end if

            if (n .eq. 2) then
                discriminante = (r**2) + (4.0d0 * s)
                if (discriminante .ge. 0.0d0) then
                    raiz(1) = (r + sqrt(discriminante)) / 2.0d0
                    raiz(2) = (r - sqrt(discriminante)) / 2.0d0
                    assinatura_cor = assinatura_cor + abs(raiz(1)) + abs(raiz(2))
                else
                    assinatura_cor = assinatura_cor + abs(r / 2.0d0) + abs(sqrt(abs(discriminante)) / 2.0d0)
                end if
            end if

            if (n .eq. 1) then
                raiz(1) = -(a(0)/a(1))
                assinatura_cor = assinatura_cor + abs(raiz(1))
            end if

            ! Salva a linha final no CSV: Coordenada X, Coordenada Y, Cor e Iterações Gastas
            write(unit_csv, '(4(F12.6, ","))') r0_start, s0_start, assinatura_cor, real(total_iter, 8)

        end do

        ! Imprime o progresso no terminal para o usuário acompanhar o processamento da malha
        if (mod(i_grid, 50) == 0) print *, "Calculando linha X = ", i_grid, " de ", res
    end do

    close(unit_csv)
    print *, "Mapeamento Completo Finalizado! Arquivo pronto para o Gnuplot."

end program Bairstow_Completo

! ----------------------------------------------------------------------------------------------------
! ROTEIRO E ANÁLISES FINAIS
! ----------------------------------------------------------------------------------------------------

! Roteiro 1. O código por inteiro é a resposta da nossa primeira pergunta do roteiro, ele está preparado para encontrar as raízes de polinômios de grau qualquer, assim como trata de forma satisfatória raízes complexas e seus pares conjugados. Ao longo do desenvolvimento do código, diversos polinômios foram testados, mas os coeficientes presentes na forma em que o código será entregue se referem ao polinômio encontrado na Atividade Para Casa 2 (APC2).

! Roteiro 2. No caso atual, em que chutamos (r=5) e (s=-3), são necessárias 7 iterações para que todas as raízes sejam encontradas, mas como iremos perceber ao gerarmos os fractais de Bairstow deste polinômio, perceberemos que existem zonas em que serão necessárias quantidades muito maiores de iterações, e zonas "caóticas" em que pequenas alterações nos valores de r e s nos levarão para resultados completamente diferentes, em que encontraremos as mesmas raízes, mas com diferentes quantidades de iterações e em ordens diferentes.

! Roteiro 3. Ao trocarmos os valores de r e s para 516 e -15, respectivamente, saímos de 7 iterações para 25 iterações. Encontramos as raízes, mas foi necessário uma quantidade mais do que quatro vezes maior de iterações para completar o método. Também é importante notar que encontramos as raízes na ordem inversa, de forma que mostramos que alterar as condições inicias afeta significativamente o fluxo do código.

! Roteiro 4. Como mencionado anteriormente, os valores presentes no código atual são o polinômio de quarta ordem encontrado na APC2, em que definimos um sistema massa mola de ordem 2 que resultou num polinômio característico de ordem 4, em que também provamos que para um sistema de ordem n, iremos obter um polinômio característico de ordem 2n. As raízes encontradas foram -65.19, -7.94, -0.100 e -0.0897, que quando analisadas com as raízes reais do sistema, apresentam diferenças muito pequenas, mostrando a eficiência do método e precisão acima de 95%.

! Roteiro 5. Em resposta à etapa de interpretação dos resultados à luz do comportamento físico do sistema (estabilidade, amortecimento, frequência de oscilação), analisamos as quatro raízes encontradas (-65.19, -7.94, -0.100 e -0.0897). Primeiramente, quanto à estabilidade, notamos que todos os autovalores são reais e estritamente negativos. Isso garante que o nosso sistema massa-mola-amortecedor é assintoticamente estável, ou seja, qualquer perturbação inicial introduzida nas massas decairá exponencialmente ao longo do tempo, fazendo o sistema retornar à sua posição de equilíbrio original. Com relação ao amortecimento e à frequência de oscilação, a ausência de raízes complexas conjugadas indica que não há formação de ciclos periódicos. Sendo assim, o sistema possui frequência de oscilação nula e atua em um regime superamortecido, em que a dissipação de energia imposta pelos amortecedores é tão "potente" que freia o sistema de forma extremamente brusca.

! Roteiro 6 e 7. A parte relacionada aos fractais de Bairstow exigiu uma lógica extra, mas conseguimos unificá-la perfeitamente ao código principal na 'Parte 2'. O programa agora replica o "coração" deflacionador em um loop massivo para varrer os valores de r e s, mapeando a estabilidade da convergência em um arquivo .csv que é interpretado pelo gnuplot. As imagens resultantes comprovam visualmente as bacias de atração do polinômio, e as rotinas de compilação estão documentadas no README.md.

! Roteiro Sumário: Este programa implementa o Método de Bairstow para extração iterativa de raízes (reais e complexas) de polinômios de grau n. O algoritmo utiliza deflação polinomial sucessiva e um sistema de Newton-Raphson bidimensional para ajustar os coeficientes de um fator quadrático divisor até que o resto da divisão seja nulo. Esta abordagem puramente real permite a análise direta da estabilidade de sistemas dinâmicos e o mapeamento do esforço computacional (caos e convergência) através da geração de fractais em malhas cartesianas.
