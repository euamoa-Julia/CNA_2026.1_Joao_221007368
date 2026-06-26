# AS 5: Matriz Inversa por LU

## Sobre o Código
O programa `inversaLU.f90` é uma extensão do algoritmo estudado na Atividade 4. Além de resolver o sistema base $AX=B$, ele demonstra o verdadeiro poder da **Fatoração LU**: calcular a Matriz Inversa ($A^{-1}$) com extremo baixo custo computacional.

### Matemática Envolvida
O cálculo convencional de uma Matriz Inversa (via matriz adjunta e determinantes) é muito lento e custoso para matrizes grandes.
A matemática por trás desse método diz que se nós resolvermos o sistema $AX = B$, onde $B$ é uma das colunas da Matriz Identidade ($I$), o vetor solução resultante será exatamente uma das colunas correspondentes da Matriz Inversa.

A Fatoração faz com que isso seja muito ágil, porque as matrizes $L$ e $U$ são calculadas **apenas uma vez** no Passo 1. Depois, o programa monta um loop e passa cada coluna da matriz Identidade nas rotinas de substituição ($LY=B$ e $UX=Y$) que são praticamente instantâneas, colhendo coluna por coluna da matriz $A^{-1}$.

### Entradas e Saídas
- **Entradas (Hardcoded):** Usa a mesma matriz simétrica de ordem $10$ do `decomp_LU.f90`.
- **Saídas (Terminal):**
  - O **Vetor Solução $X$** da equação original.
  - A matriz inteira contendo os coeficientes da **Matriz Inversa**.
  - **Provas Reais (Validações):** Para validar se o cálculo da inversa funcionou, o programa encerra multiplicando a Matriz original $A$ pela Inversa calculada ($A^{-1}$) e imprime o resultado na tela. Se o código estiver correto, o resultado exibido tem que ser perfeitamente a matriz Identidade ($1.0$ na diagonal e $\approx 0.0$ no restante).
