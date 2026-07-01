# PPC 3: Transferência de Calor em Reator Nuclear

## Sobre o Código
O programa `nuclear_internal.f90` aborda a solução computacional de uma Equação Diferencial Parcial (EDP) parabólica que governa a condução de calor transiente em uma dimensão. Ele simula o perfil de temperaturas através da espessura de uma parede (pastilha) de urânio em um reator nuclear, submetida a resfriamento por fluido em uma face.

### Matemática Envolvida
O domínio espacial (a espessura da parede $L$) foi discretizado em $7$ nós usando o método de **Diferenças Finitas**. O domínio temporal também foi discretizado e resolvido usando a **Formulação Implícita**.

A grande sacada da formulação implícita é que a equação do nó atual no tempo futuro ($T_i^{p+1}$) depende das temperaturas adjacentes *também* no tempo futuro ($T_{i-1}^{p+1}$ e $T_{i+1}^{p+1}$). Isso gera um Sistema Linear para cada instante de tempo, onde todos os nós devem ser resolvidos simultaneamente.
A vantagem é que o esquema Implícito é **Incondicionalmente Estável**, não explodindo para o infinito mesmo usando passos de tempo ($\Delta t$) muito agressivos, escapando das restrições de Fourier e de Biot.

O sistema tridiagonal resultante ($AX=B$) é superotimizado e resolvido matematicamente por meio do **Algoritmo de Thomas (TDMA)**, que é uma versão especializada e enxuta da Fatoração LU para matrizes tridiagonais, reduzindo a complexidade algorítmica de $O(n^3)$ para apenas $O(n)$.

A equação fundamental engloba o termo $\dot{q}$ de Geração Interna de Calor no centro do sistema:
$$\frac{\partial^2 T}{\partial x^2} + \frac{\dot{q}}{k} = \frac{1}{\alpha} \frac{\partial T}{\partial t}$$

### Entradas e Saídas
- **Entradas (Menu Interativo):** Ao rodar o executável, o terminal pergunta se o usuário quer rodar a simulação `"com"` geração de calor interna (a reação nuclear ativada $\dot{q} > 0$) ou `"sem"` (um limite assintótico usado para testes).
- **Entradas (Ajustáveis no Código):** Passos $\Delta t$ e $\Delta x$, condutividade térmica, densidade, $h$ (convecção), e Temperatura ambiente do fluido.
- **Saídas (Terminal):** Um relatório dinâmico a cada segundo da simulação revelando a subida de temperatura nos 7 pontos avaliados.
- **Saídas (Arquivos):** Grava um arquivo formatado no estilo Tabela (`linhas.dat`) e um arquivo para mapas topográficos (`mapa.dat`).

### Gráficos
O script `plot_reactor.gp` compila esses dados de saída em dois excelentes modos de visualização:
- **`curvas.png`:** Um gráfico de linhas comum 2D onde o Eixo X é o tempo e as linhas coloridas representam a evolução térmica isolada de cada um dos 7 nós (do centro quente até a borda resfriada).
- **`campo.png`:** Utiliza o comando `pm3d map` para projetar uma vista superior do reator. O eixo Y é o tempo, o Eixo X é o espaço na espessura da parede, e a terceira dimensão (Temperatura) vira um **Mapa de Cores** contínuo (De Azul escuro a Vermelho Brilhante).
