# AS 4: Resolução de Sistemas Lineares - Decomposição LU

## Sobre o Código
O código `decomp_LU.f90` tem o objetivo de solucionar grandes sistemas de equações lineares simétricos (formato matricial $AX = B$) sem usar a tradicional inversão de matriz ou Eliminação de Gauss bruta. Para isso, ele utiliza o algoritmo de **Fatoração LU**.

### Matemática Envolvida
Do ponto de vista da engenharia de software e análise numérica, resolver sistemas mudando apenas o vetor de forças $B$ usando Gauss é muito ineficiente, pois a matriz $A$ precisaria ser reescalonada inteira do zero.
O método da Decomposição LU divide a matriz global de coeficientes $A$ no produto de duas matrizes triangulares:
$$A = L \times U$$
- **Matriz L (Lower):** Matriz triangular inferior contendo os fatores multiplicadores da eliminação, com $1$s na diagonal principal.
- **Matriz U (Upper):** Matriz triangular superior, que é basicamente o que sobra da matriz $A$ após a Eliminação de Gauss tradicional.

Após decompor $A$ (apenas uma vez), resolvemos o sistema para qualquer vetor $B$ por meio de substituições simples e extremamente rápidas com o auxílio de um vetor temporário ($Y$):
1. $LY = B$ (Substituição progressiva, resolve de cima para baixo).
2. $UX = Y$ (Substituição regressiva, resolve de baixo para cima entregando a resposta final).

### Entradas e Saídas
- **Entradas (Hardcoded):**
  - A ordem da matriz quadrada, neste caso, configurado para `n = 10`.
  - Matriz de Coeficientes (`A` de dimensão $10 \times 10$).
  - Vetor de Forças/Termos Independentes (`B` de dimensão 10).
- **Saídas (Terminal):**
  - Mostra as matrizes `L` e `U` separadas após a fatoração.
  - Imprime o **Vetor Solução $X$**.
  - **Provas Reais (Validações):** O código conta com lógicas internas de prova. Ele reconstrói e multiplica as matrizes $L \times U$ no terminal para provar matematicamente que o produto das duas retorna perfeitamente aos coeficientes da matriz $A$ original, e cruza o vetor $X$ encontrado no sistema para garantir que recuperamos os valores exatos de $B$.
