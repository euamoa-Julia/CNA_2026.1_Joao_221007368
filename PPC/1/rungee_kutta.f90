! Damos início ao código com o comando "program", Em seguida iremos digitar "implicit none" de forma que
! teremos de declarar todas as nossa variáveis, evitando que o código interprete algo de forma errônea
! e também nos dando mais autonomia sobre o código, permitindo total controle sobre as variáveis.

program RK4_sedimentacao
    implicit none

    ! Agora iremos declarar todas as variáveis que serão utilizadas nos cálculos:
    real(8) :: t       ! tempo inicial
    real(8) :: v       ! velocidade inicial
    real(8) :: h       ! passo ou avanço no tempo
    real(8) :: t_final ! tempo final
    real(8) :: S       ! número de Stokes
    real(8) :: Re      ! número de Reynolds, que no caso deste código será variado de 0 à 1 em cinco intervalos com avanço de 0,2
    real(8) :: k1      ! K1 Runge-Kutta
    real(8) :: k2      ! K2 Runge-Kutta
    real(8) :: k3      ! K3 Runge-Kutta
    real(8) :: k4      ! K4 Runge-Kutta

    ! Nestas próximas linhas iremos definir os valores inicias do tempo e da velocidade, assim como o passo, o tempo final.
    ! Manteremos S constante e Re será variado conforme comentário anterior. Iremos também variar o h (passo) para entender seu efeito no comportamento do método.
    t       = 0.0d0  ! Tempo inicial.
    v       = 0.0d0  ! Velocidade inicial.
    h       = 0.1d0  ! Passo variando no intervalo de 0.01 à 1 (0.01, 0.1, 0.5 e 1).
    t_final = 10.0d0 ! Optei pelo intervalo de tempo de 0 a 10.
    S       = 0.8d0  ! Stokes variando no intervalo de 0.2 à 1 (0.2, 0.4, 0.6, 0.8 e 1)
    Re      = 0.0d0  ! Reynolds variando no intervalo de 0 à 1 (0, 0.2, 0.4, 0.6, 0.8 e 1).

    ! O comando "open" será utilizado para criar um arquivo com os dados obtidos ao rodarmos o código.
    ! O comando "unit" será utilizado para definir onde o arquivo será guardado.
    ! O comando "file" irá definir o nome e formato do arquivo com os resultados, sendo fundamental para confecção dos gráficos.
    ! O comando "replace" define que a cada vez que rodarmos o código, os dados do arquivo .dat serão sobreescritos.
    open (unit=10, file='dados_velocidade_S08.dat', status='replace')

    write(*,*) "--------------------------------------------------------------------------------"
    write(*,*) "                 SIMULACAO RK4 - SEDIMENTACAO DA PARTICULA                      "
    write(*,*) "--------------------------------------------------------------------------------"
    write(*,"(A, F8.4)") " Passo (h):          ", h
    write(*,"(A, F8.4)") " Stokes (S):         ", S
    write(*,"(A, F8.4)") " Reynolds (Re):      ", Re
    write(*,*) "--------------------------------------------------------------------------------"
    write(*,*) " Calculando e salvando dados no arquivo 'dados_velocidade_S08.dat'..."

    ! Agora iremos iniciar a "captura" dos dados com o laço "while", que define que enquanto t não
    ! chegar até o valor de t_final, que no caso é 10, o código irá continuar rodando e capturando dados em t=t+0.1.
    do while(t <= t_final)

        write (10, '(2F15.8)') t, v

        ! As linhas seguintes serão utilizadas para calcular a parte principal do código, o método númerico.
        ! "f_derivada" é a função adimensionalizada que encontramos na primeira parte do trabalho, ou seja, no APC.
        ! Definimos essa função através do comando "contains", que definirá funções que podem ou não serem utilizadas
        ! ao longo do código. Funciona como uma espécie de biblioteca dentro do código com funções que não necessariamente
        ! serão utilizadas, mas estarão disponíveis.

        ! CORREÇÃO: Substituição de (1/2) por 0.5d0 para evitar divisão inteira que zera o avanço do RK4.
        k1 = f_derivada(t, v, S, Re)
        k2 = f_derivada(t + 0.5d0 * h, v + 0.5d0 * k1 * h, S, Re)
        k3 = f_derivada(t + 0.5d0 * h, v + 0.5d0 * k2 * h, S, Re)
        k4 = f_derivada(t + h, v + k3 * h, S, Re)

        ! CORREÇÃO: (h/6) alterado para (h/6.0d0) pela mesma razão.
        v = v + (h / 6.0d0) * (k1 + 2.0d0 * k2 + 2.0d0 * k3 + k4)
        t = t + h

    end do

    close (10)

    write(*,*) " Simulação concluída com sucesso!"
    write(*,*) "--------------------------------------------------------------------------------"

contains

    real(8) function f_derivada(t_func, v_func, S_func, Re_func)
        ! Respondendo ao comentário: Mesmo que o 't' não apareça na equação da derivada (equação autônoma),
        ! o formato padrão do RK4 exige que a função receba o tempo para manter a generalidade do método!
        real(8), intent(in) :: t_func
        real(8), intent(in) :: v_func
        real(8), intent(in) :: S_func
        real(8), intent(in) :: Re_func

        ! Função adimensionalizada encontrada na primeira parte do trabalho.
        f_derivada = (-v_func/S_func) - (3.0d0/(8.0d0*S_func)) * Re_func * (v_func**2) + (1.0d0/S_func)

    end function f_derivada

! Agora finalizaremos o programa e deixaremos o código pronto para ser executado através do terminal, pelo comando "gfortran RK_sedimentacao.f90 -o simulacao"

end program RK4_sedimentacao

! ----------------------------------------------------------------------------------------------------------------------
! RESPOSTAS DO ROTEIRO
! ----------------------------------------------------------------------------------------------------------------------
! Este código foi utilizado para a confecção de dois gráficos, um em que mostramos o efeito da variação de Reynolds e
! outro que mostramos a diferença do comportamento do método quando variamos o passo, e em ambos os casos utilizamos os
! mesmos parâmetros de ínicio e fim, assim como consideramos (Stokes=1) para os casos da variação de Reynolds e do passo.

! Roteiro 5. No caso do grafico_reynolds.png, percebemos que quanto menor o Reynolds, maior a velocidade com que a
! partícula atinge sua velocidade terminal, assim como sua velocidade terminal também é maior. Isso ocorre por que o
! numero de Reynolds indica a razão entre as forças inérciais e viscosas, de forma que o regime ditado por (Re=0) é
! dado que as forças viscosas são praticamente infinitamente maiores que as forças inérciais, removendo parte da
! resistência que a partícula encontraria ao sendimentar, sendo de fácil visualização que um número de Reynolds maior
! dificultaria a sua movimentação no meio. Quando analisamos as curvas dos números de Reynolds maiores, percebemos
! claramente uma dificuldade maior da partícula ao se mover no meio, com uma menor aceleração e velocidade terminal.

! Roteiro 4. Comparação realizada logo abaixo + gráfico_stokes.png.

! Roteiro 3. Código + Gráficos.

! Roteiro 2. Analisando agora o grafico_passo.png em que variamos o passo, percebemos uma grande variação na precisão
! do "meio" do método, de forma que os resultados possuem variação significativa (exceto entre 0.1 e 0.01). Quando
! analisamos o resultado final, percebemos que existe uma grande variação a partir da quinta casa decimal, e isso
! ocorre por que o passo possui um efeito direto no cálculo dos K's utilizados no método, mas também por que o cálculo
! de K1 afeta o cálculo de K2 e assim por diante, de forma que a alteração no passo pode aumentar ou diminuir a
! precisão do método, mas também aumenta o custo computacional, de forma que passos menores exigem um tempo maior de
! cálculo, assim como processamento adicional, com o passo de (h=0.1) sendo o caso ideal em que percebemos uma precisão
! tão grande quanto (h=0.01) e uma diferença significativa para os outros métodos.

! Roteiro 1. Já com o grafico_stokes.png em que mantemos o regime de (Re=0) e variamos o número de Stokes, e o que
! percebemos é que o número Stokes tem um efeito "negativo" na velocidade da partícula, dificultando sua aceleração
! e aumentando o tempo necessário para atingir a sua velocidade terminal. O número de Stokes pode ser definido como a
! latência do fluido/partícula, de forma que indica o tempo que um leva para perceber o outro, e quanto maior essa
! latência, mais inércia existe na relação. A comparação com a solução analítica permite que o método numérico seja
! validado com uma base real para comparação, ou seja, temos uma forma de validar cada ponto no gráfico. Ao fazer
! isso, percebemos que a solução numérica é extremamente próxima da solução analítica, com os gráficos se sobrepondo
! sobre a grande maioria do percurso da curva (Velocidade x Tempo). O argumento de que as soluções numéricas são
! fundamentais na eficiência de equações diferenciais se mantém.

! ______________________________________________________________________________________________________________________

! Roteiro 5. Comparando os resultados para diferentes Reynolds e com passo h=0.1, obtemos o seguinte:

! Reynolds | Velocidade Terminal Analítica | Velocidade Terminal Numérica
! 0        | 1                             | 0.99996455496668524
! 0.2      | 0.9345                        | 0.93449487948105447
! 0.4      | 0.8830                        | 0.88303481301055242
! 0.6      | 0.8409                        | 0.84089912057000715
! 0.8      | 0.8054                        | 0.80539930587404507
! 1        | 0.7749                        | 0.77485171027814892

! Roteiro 5. Apesar dos resultados apresentarem grande proximidade, mostrando que a solução analítica também é viável,
! ela se mostra consideravelmente mais complicada de se desenvolver, exigindo o cálculo de fatores como Q e P, fatores
! estes presentes na equação de Ricatti que deve ser resolvida para encontrarmos o valor analítico da expressão. O
! método numérico permite uma implementação em código de fácil desenvolvimento, poupando tempo e sendo consideravelmente
! mais simples, além de que o código pode ser alterado sem grandes dificuldades para a resolução de outros problemas
! semelhantes em questão de minutos.

