# AS 2: Método de Newton-Raphson

## Sobre o Código
O programa `newton_raphson.f90` implementa o método aberto de **Newton-Raphson** para encontrar raízes de equações polinomiais e não-lineares. O código está preparado com a fórmula padrão de Newton-Raphson e com a modificação estabilizada para lidar com raízes múltiplas/repetidas.

### Matemática Envolvida
Diferente dos métodos de confinamento que requerem dois pontos iniciais envolvendo a raiz, o Newton-Raphson é um método **aberto** que exige apenas um "chute inicial" ($x_0$).
A partir do ponto, ele usa o cálculo da **Derivada Analítica** ($f'(x)$) para traçar uma reta tangente. Onde a reta tangente intercepta o eixo horizontal ($X$) é a próxima estimativa ($x_{novo}$).

- **Fórmula Padrão:** $x_{i+1} = x_i - \frac{f(x_i)}{f'(x_i)}$
- **Fórmula Modificada (Raízes Múltiplas):** Quando a raiz não "cruza" o eixo $X$ mas apenas "tangencia" (com a inclinação se aproximando de $0$), a fórmula padrão causa divisão por zero. A fórmula modificada insere a Segunda Derivada ($f''(x)$) no denominador para contornar isso e manter a convergência rápida:
  $$x_{i+1} = x_i - \frac{f(x_i)f'(x_i)}{[f'(x_i)]^2 - f(x_i)f''(x_i)}$$

O código resolve especificamente a equação polinomial: $f(x) = x^3 - 6x^2 + 9x - 4$.

### Entradas e Saídas
- **Entradas:**
  - Um chute inicial $x_{atual} = 2.5$.
  - Tolerância de parada `tol = 1.0d-6` e um limite de iterações seguro `max_iter = 100`.
- **Saídas (Terminal):** Uma tabela com as colunas das iterações iterativas ($x_{atual}$, $f(x)$ e $Erro$). No final do loop, o programa exibe a raiz encontrada e a quantidade de iterações gastas.
- **Saídas (Arquivos):** Grava o log com a raiz otimizada final no arquivo de texto `raizes.dat`.
