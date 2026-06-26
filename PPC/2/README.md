# PPC 2: Método de Bairstow e Geração de Fractais

## Sobre o Código
O programa `bairstow.f90` é um algoritmo imenso dividido em duas "Engrenagens" (Partes).
A Parte 1 é um resolvedor focado em implementar o **Método de Bairstow**, capaz de arrancar iterativamente **todas as raízes (reais ou complexas)** de um polinômio genérico usando divisões sintéticas sucessivas (deflação) e um sistema 2D de Newton-Raphson. A Parte 2 usa esse motor deflacionário num loop infinito para "varrer" todas as coordenadas do plano e mapear as bacias de atração do polinômio, gerando um **Fractal de Bairstow**.

### Matemática Envolvida (O Método de Bairstow)
O programa foge da divisão complexa por polinômios de grau 1 e prefere dividir o Polinômio original ($P_n(x)$) por um Fator Quadrático Real ($x^2 - rx - s$).
O objetivo é calibrar os chutes das variáveis $r$ e $s$ até que o resto da divisão polinomial dê zero (sendo exata). Quando a divisão fica exata, significa que encontramos 2 raízes do problema na forma do fator quadrático extraído. As raízes podem ser obtidas rapidamente jogando Bhaskara em $x^2 - rx - s = 0$ (o que lida perfeitamente com discriminantes reais ou negativos geradores de complexos conjugados).
Para calibrar as correções $\Delta r$ e $\Delta s$ até o limite exato, o código constrói um Sistema Linear montado via derivadas parciais (Newton-Raphson 2D) sobre as funções de "resto" da divisão.
Isso tudo ocorre em um super laço: A cada 2 raízes extraídas, o polinômio é deflacionado (seu grau cai em 2) e o ciclo se repete até ser aniquilado em um polinômio de grau 1 ou grau 2 puro.

### Entradas e Saídas
#### Parte 1 (Resolvendo o Sistema)
- **Entradas:** Vetor de coeficientes da equação em ordem decrescente (Polinômio base de grau 4 derivado da equação massa-mola). E as estimativas iniciais de "calibração", $r$ e $s$.
- **Saídas (Terminal):** Exibe a extração das raízes na tela durante as deflações graduais, mostrando inclusive os pares complexos e a quantidade de passos e iterações.

#### Parte 2 (Mapeamento Cartesiano do Fractal)
- **Entradas:** Usa o mesmo polinômio base, mas atua definindo o tamanho da "Grade" / "Resolução" da imagem simulada com as variáveis `r_min`, `r_max`, `s_min`, e `s_max`.
- **Saídas (Arquivos):** Exporta o colossal arquivo bidimensional `fractal_bairstow.csv`. Nele constam:
  - As coordenadas espaciais `r0` e `s0` mapeadas.
  - O cômputo da "Cor" (Uma assinatura geométrica baseada em quais raízes foram descobertas e sua posição).
  - O índice caótico "Iterações" (A quantidade de suor que o algoritmo derramou para convergir começando daquela coordenada específica).

Esses dados CSV são os blocos de montagem renderizados visualmente usando o módulo Gnuplot externo em um mapa de cores (pm3d).
