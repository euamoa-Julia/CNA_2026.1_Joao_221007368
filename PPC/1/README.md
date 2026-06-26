# PPC 1: Solução de EDOs com Runge-Kutta 4

## Sobre o Código
O programa `rungee_kutta.f90` (presente na subpasta e nas ramificações) modela matematicamente o problema físico da velocidade de sedimentação de uma partícula decaindo dentro de um fluido. Sendo um problema dependente do tempo, resulta em uma Equação Diferencial Ordinária (EDO) de valor inicial, a qual é resolvida no código via o consagrado método numérico de **Runge-Kutta de 4ª Ordem (RK4)**.

### Matemática Envolvida
O RK4 trabalha projetando passos ao longo da curva baseados não em apenas uma inclinação (derivada), mas em uma média ponderada de quatro inclinações distintas coletadas durante um único passo:
1. $k_1$: A derivada calculada no ponto inicial (Euler explícito).
2. $k_2$: A derivada simulada no ponto médio, com a altura influenciada por $k_1$.
3. $k_3$: Outra derivada simulada no ponto médio, com a altura influenciada por $k_2$.
4. $k_4$: A derivada calculada ao final do passo de tempo, influenciada por $k_3$.

A nova velocidade é calculada pela combinação: $v_{i+1} = v_i + \frac{h}{6}(k_1 + 2k_2 + 2k_3 + k_4)$.

A derivada que alimenta a EDO ($f(t, v)$) foi previamente adimensionalizada baseando-se no Número de Stokes ($S$) e no Número de Reynolds ($Re$):
$$\frac{dv}{dt} = f(t,v) = -\frac{v}{S} - \frac{3}{8S} \cdot Re \cdot v^2 + \frac{1}{S}$$

### Entradas e Saídas
- **Entradas (Ajustáveis dentro do código):**
  - Passo de tempo (`h`), Número de Stokes (`S`), e Número de Reynolds (`Re`).
  - Tempo Inicial e Final da simulação.
  - Velocidade inicial $v(0) = 0$.
- **Saídas (Terminal):** Exibe logs interativos de "Simulação Concluída".
- **Saídas (Arquivos):** Gera um arquivo de banco de dados (`dados_velocidade_S08.dat`) contendo colunas de Tempo e Velocidade. Este arquivo armazena o histórico da aceleração até que a partícula atinja sua Velocidade Terminal (velocidade constante).

### Gráficos
No diretório deste PPC, são criados scripts para gerar gráficos comparativos variando os parâmetros físicos e computacionais:
- **`grafico_reynolds.png`:** Gráfico que sobrepõe a simulação para diferentes números de Reynolds (de 0 a 1). Revela fisicamente como regimes com maiores inércias relativas prejudicam a velocidade terminal da partícula.
- **`grafico_passo.png`:** Gráfico analítico sobre o método numérico. Variando o passo computacional $h$, ilustra o "trade-off" entre custo computacional e estabilidade. Passos maiores causam saltos imprecisos na curva.
- **`grafico_stokes.png`:** Foca no impacto da latência das partículas (Stokes).
