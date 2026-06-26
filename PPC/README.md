# Programas Para Casa (PPC)

Esta pasta contém programas completos e robustos que representam os grandes trabalhos de casa e projetos (Assignments) exigidos ao longo da disciplina. Eles unem múltiplos conceitos numéricos, modelagem física e análise crítica de desempenho computacional.

## Projetos

- **PPC 1 (`PPC/1/rungee_kutta.f90`) - Simulação de Sedimentação**:
  - Resolve uma Equação Diferencial Ordinária (EDO) modelando a velocidade de decaimento de uma partícula no fluido usando o método **Runge-Kutta de 4ª Ordem (RK4)**.
  - O código estuda os impactos dos adimensionais de Reynolds e Stokes na convergência e nas propriedades físicas do escoamento, salvando logs para análise gráfica das curvas de velocidade vs tempo.
- **PPC 2 (`PPC/2/bairstow.f90`) - Raízes Complexas e Mapeamento de Fractais**:
  - Implementa o **Método de Bairstow**, uma técnica avançada que extrai raízes de polinômios (incluindo raízes complexas não conjugadas) por meio da deflação iterativa por fatores quadráticos.
  - A Parte 2 deste código adapta o método para mapear o comportamento caótico da convergência de múltiplos chutes iniciais em uma malha bidimensional, viabilizando a geração do famoso Fractal de Bairstow em alta resolução via CSV.
- **PPC 3 (`PPC/3/nuclear_internal.f90`) - Transferência de Calor em Reator Nuclear**:
  - Modela o perfil de temperaturas de uma parede de um reator nuclear empregando o método de **Diferenças Finitas em Formulação Implícita**.
  - O complexo sistema tridiagonal originado na formulação temporal implícita é resolvido matematicamente usando o eficiente **Algoritmo de Thomas (TDMA)**.
  - O código tem um menu interativo que permite estudar o cenário estacionário e o impacto do termo de Geração Interna de Calor no transiente térmico.
- **PPC 4 (`PPC/4/PPC4.f90`) - Otimização Multidimensional Topográfica**:
  - Um algoritmo massivo que mescla buscas lineares e direcionais para navegar por uma superfície 3D em busca do cume da montanha (ponto de máximo global de uma função de duas variáveis).
  - Inclui e compara o desempenho de dois algoritmos de direção de passo: **Aclive Máximo (Steepest Ascent)** e o algoritmo iterativo de **Fletcher-Reeves (Gradientes Conjugados)**.
  - A fase de Busca Linear (o quão longe saltar na direção) é feita acoplando a rotina de *Interpolação Quadrática* abordada nas atividades de sala.

## Como utilizar

Ao compilar esses projetos, esteja ciente de que muitos deles foram desenhados para gravar grandes massas de dados. O `PPC2`, por exemplo, gera malhas de fractais CSV que podem ter dezenas de megabytes dependendo da resolução (variável `res`).

1. Navegue até o projeto desejado:
   ```bash
   cd 1
   ```
2. Compile com o gfortran:
   ```bash
   gfortran nome_do_arquivo.f90 -o exec_ppc
   ```
3. Execute o programa:
   ```bash
   ./exec_ppc
   ```
4. Se o programa possuir rotinas gráficas atreladas ou exportar arquivos de dados (`.dat`, `.csv`), você poderá analisá-los utilizando scripts auxiliares.
