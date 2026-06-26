# AS 3: Método de Muller

## Sobre o Código
O programa `muller.f90` implementa uma técnica avançada de busca de raízes chamada de **Método de Muller**. Em vez de aproximar a curva original por linhas retas, como fazem as Secantes ou a Falsa Posição, o Método de Muller aproxima a curva por parábolas (interpolação quadrática).

### Matemática Envolvida
O método exige três aproximações iniciais ($x_0, x_1, x_2$) e calcula uma parábola $P(x) = a x^2 + b x + c$ que passa por esses três pontos usando diferenças divididas. Onde essa parábola cruza o eixo $X$ torna-se a nova estimativa da raiz ($x_3$).

O grande trunfo de Muller é que, ao usar a fórmula de Bhaskara para encontrar as raízes da parábola:
$$x_3 = x_2 - \frac{2c}{b \pm \sqrt{b^2 - 4ac}}$$
Ele pode produzir naturalmente **raízes complexas**, algo que o Newton-Raphson não faz a menos que o chute inicial seja complexo. Além disso, a convergência deste método é superlinear.

*Obs: A versão em código `real(8)` deste repositório foi construída com um bloqueio (parada) se o discriminante for negativo. Para atuar na zona complexa, basta converter os tipos para `complex(8)`.*

O código busca a raiz do polinômio: $f(x) = x^3 - 13x - 12$.

### Entradas e Saídas
- **Entradas:**
  - Um chute inicial base ($x_0 = 6.5$) e um passo (`passo = 0.5`) que o programa usa para recuar e gerar $x_1$ e $x_2$ automaticamente.
  - Critérios de parada `tol = 1.0d-6` e `max_iter = 100`.
- **Saídas (Terminal):** Exibe o recálculo dos novos pontos, os valores da função nestes pontos e o erro relativo caindo ao longo das iterações, confirmando a alta eficiência.
- **Saídas (Arquivos):** Exporta o valor final da raiz computada para o arquivo externo `raizes.dat`.
