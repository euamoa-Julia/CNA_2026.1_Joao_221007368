program reactor

implicit none

real(8), dimension(7) :: e, f, g, r, x ! Vetores das diagonais, termos independentes e temperaturas
integer k, n, num_nos                  ! Contadores do laço e quantidade de nós
real(8) fo, bi, alpha                  ! Adimensionais (Fourier, Biot) e Difusividade térmica
real(8) k_therm, rho, cp_mat, h, l, t_inf, t_init, tempo_total, tempo_passado
real(8) dt, dx                         ! Tamanho do passo de tempo e passo espacial
real(8) q_dot, termo_A                 ! Geração de calor interna e termo modificado
character (len=3) :: geracao           ! Entrada do usuário ("com" ou "sem")

logical :: entrada_valida ! Variável lógica para controle do menu interativo

entrada_valida = .false.

! ==============================================================================
! PASSOS 1 e 2: Seleção do limite assintótico (q_dot -> 0) e termo de geração
! ==============================================================================
! Este laço roda até o usuário fornecer uma opção válida, permitindo usar
! a mesma malha (dx, dt) para validar o caso com e sem geração interna.
do while (.not. entrada_valida)
    
    print *, "Digite 'com' ou 'sem' para simular com ou sem geracao:"
    read(*, *) geracao

    if (trim(geracao) == "sem") then
        ! Passo 1: Limite assintótico sem geração para validação com solução exata
        q_dot = 0.0d0 
        print *, "Simulando SEM geracao de calor. Q_dot =", q_dot
        entrada_valida = .true. 

    else if (TRIM(geracao) == "com") then
        ! Passo 2: Adição do termo de geração interna de calor
        q_dot = 50000000.0d0 
        print *, "Simulando COM geracao de calor. Q_dot =", q_dot
        entrada_valida = .true.

    else
        print *, "Opcao invalida! Tente novamente."
        print *, "" 
    end if

end do

print *, "Sucesso! Continuando a simulacao do reator..."
 
! ==============================================================================
! PASSO 4: Parâmetros realistas para o problema do reator nuclear
! ==============================================================================
num_nos = 7
n = num_nos

k_therm = 3.0d0      ! Condutividade térmica [W/m.K]
rho = 10300.0d0      ! Massa específica [kg/m³] (Típico para pastilhas de UO2)
cp_mat = 300.0d0     ! Calor específico [J/kg.K]
h = 20000.0d0        ! Coeficiente de convecção [W/m².K] (Alta transferência típica de reatores)
l = 0.01d0           ! Espessura característica da parede [m]
t_inf = 270.0d0      ! Temperatura do fluido refrigerante externo
t_init = 270.0d0     ! Temperatura inicial uniforme do material

! Malha de discretização espaço-tempo
tempo_passado = 0.0d0
tempo_total = 300.0d0
dx = l / real(num_nos-1) ! Passo espacial entre os 7 nós
dt = 1.0d0               ! Passo de tempo

! Cálculo dos adimensionais e parâmetros para o esquema numérico
termo_A = (q_dot * dt) / (rho * cp_mat) ! Aumento da temperatura por passo via geração

alpha = k_therm/(rho*cp_mat) ! Difusividade térmica
bi = (h*l)/k_therm           ! Número de Biot
fo = (alpha*dt)/(dx**2.0d0)  ! Fator de Fourier por passo de tempo

! Condição Inicial: Atribui a t_init para todos os 7 nós no tempo zero
x = [t_init, t_init, t_init, t_init, t_init, t_init, t_init]

open(unit=10, file='mapa.dat', status='replace')
open(unit=11, file='linhas.dat', status='replace')

! ==============================================================================
! LAÇO TEMPORAL: Evolução no tempo usando ESQUEMA IMPLÍCITO
! ==============================================================================
do while (tempo_passado.lt.tempo_total)

    tempo_passado = tempo_passado + dt

    ! 1. Montagem da matriz tridiagonal: Redefinição a cada passo, 
    ! pois o Algoritmo de Thomas sobrescreve os vetores originais.
    ! f = Diagonal principal
    ! e = Sub-diagonal (abaixo da principal)
    ! g = Super-diagonal (acima da principal)
    
    ! O nó 1 (x=0) tem condição de simetria (fluxo nulo). O nó 7 (x=L) sofre convecção (Biot).
    f = [1.0d0+2.0d0*fo, 1.0d0+2.0d0*fo, 1.0d0+2.0d0*fo, 1.0d0+2.0d0*fo, 1.0d0+2.0d0*fo, 1.0d0+2.0d0*fo, 1.0d0+2.0d0*fo + 2.0d0*fo*bi]
    e = [0.0d0, -fo, -fo, -fo, -fo, -fo, -2.0d0*fo] ! O primeiro termo não é utilizado, então igualamos a zero para manter os vetores com mesma dimensão.
    g = [-2.0d0*fo, -fo, -fo, -fo, -fo, -fo, 0.0d0] ! O último termo não é utilizado, então igualamos a zero para manter os vetores com mesma dimensão.
    
    ! r = Vetor independente contendo o histórico T(p) e condições de contorno fixas
    r = [t_init, t_init, t_init, t_init, t_init, t_init, 2.0d0*fo*bi*(dx/l)*t_inf]
    
    ! 2. Atualiza R com base no vetor das temperaturas calculadas no passo anterior
    r(1:6) = x(1:6) + termo_A
    r(7) = x(7) + (2.0d0 * fo * bi * t_inf) + ((fo * q_dot * dx**2.0d0) / k_therm)

    ! ==============================================================================
    ! ALGORITMO DE THOMAS (Código auxiliar do Prof. Rafael Gabler)
    ! Resolução de sistemas tridiagonais Ax = R.
    ! ==============================================================================

    ! Etapa A: Eliminação progressiva (Zera a sub-diagonal 'e')
    do k = 2, n
        e(k) = e(k)/f(k-1)
        f(k) = f(k)-e(k)*g(k-1)
    end do

    ! Modificação simultânea do vetor independente 'r'
    do k = 2, n
        r(k) = r(k)-e(k)*r(k-1)
    end do

    ! Etapa B: Substituição regressiva
    x(n) = r(n)/f(n) ! Calcula a temperatura do último nó (nó 7)

    ! Encontra a temperatura dos nós restantes voltando do final para o início
    do k = n-1, 1, -1
        x(k) = (r(k)-g(k)*x(k+1))/f(k)
    end do
    ! ==============================================================================
    ! FIM DO ALGORITMO DE THOMAS
    ! ==============================================================================
    
    ! Output do laço no terminal a cada passo de tempo
    write(*, "(A, F8.2, A, 7F8.2)") "Tempo: ", tempo_passado, " | Temp: ", x(1:7)
    
    ! ==============================================================================
    ! Armazenamento para representação visual dos resultados
    ! ==============================================================================
    
    ! Arquivo 1: Formatação para Mapa Térmico/Cores no Gnuplot (comando pm3d)
    ! Requer blocos tridimensionais (Posição, Tempo, Temperatura) separados por linhas em branco
    do k = 1, num_nos
        write(10, "(F10.5, X, F8.2, X, F10.4)") real(k-1)*dx, tempo_passado, x(k)
    end do
    write(10, *) ! Linha vazia obrigatória para separação de blocos no Gnuplot
    
    ! Arquivo 2: Formatação tabular para Gráficos de Perfis e Linhas
    ! Coluna 1 = Tempo; Colunas de 2 a 8 = Temperaturas nos nós 1 a 7
    write(11, "(F8.2, 7(X, F10.4))") tempo_passado, x(1:7)

end do

close(10)
close(11)

end program reactor
