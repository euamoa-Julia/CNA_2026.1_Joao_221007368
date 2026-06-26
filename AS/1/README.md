# AS 1: Métodos de Confinamento - Bissecção e Falsa Posição

## Sobre o Código
O programa `atividade1.f90` implementa e compara dois métodos clássicos de confinamento para encontrar raízes de equações não-lineares (funções $f(x)=0$). Os métodos são aplicados na solução da temperatura de um corpo baseado em uma equação de transferência de calor.

### Matemática Envolvida
1. **Método da Bissecção (Busca Binária):**
   O algoritmo divide o intervalo $[x_i, x_s]$ exatamente ao meio ($x_m$). Usando o Teorema de Bolzano ($f(x_i) \times f(x_s) < 0$), verifica em qual metade houve a inversão de sinal (a raiz está lá) e repete a divisão. Possui convergência lenta, mas matematicamente garantida.
2. **Método da Falsa Posição:**
   Em vez de dividir ao meio, o método traça uma reta secante entre os pontos $(x_i, f(x_i))$ e $(x_s, f(x_s))$. Onde a reta cruza o eixo $X$ ($x_r$), torna-se a nova estimativa da raiz. Este método utiliza a inclinação (peso) dos valores da função para tentar se aproximar mais rapidamente da raiz do que a Bissecção.

A equação avaliada pelo código é:
$$f(x) = \frac{667.38}{x} \cdot (1.0 - e^{-0.146843 x}) - 40.0$$

### Entradas e Saídas
- **Entradas:** Os limites iniciais $x_i = 12.0$ e $x_s = 16.0$, e a tolerância de erro `tol = 1.0d-5`.
- **Saídas (Terminal):** Tabelas interativas mostrando o passo a passo da convergência, os pontos calculados e o Erro Relativo, com a raiz final exibida.
- **Saídas (Arquivos):** Os arquivos `bisseccao.dat` e `falsa_posicao.dat` (gerados automaticamente ou simulados como `saida.dat` e `saida2.dat`) salvam o registro de todas as iterações e do erro verdadeiro comparado ao erro numérico relativo.

### Gráficos
Nesta pasta, estão incluídos scripts `.gnu` para o Gnuplot utilizarem os arquivos `.dat`:
- **`plot1.gnu` (Comparação de Erros - Bissecção):** Gera a imagem `comparacao_erros.png` plotando o Erro Verdadeiro e o Erro Relativo da Bissecção em escala logarítmica com o passar das iterações. Mostra o decaimento em zigue-zague natural do método.
- **`plot2.gnu` (Comparação de Métodos):** Gera a imagem `comparacao_metodos.png` agrupando as curvas de erro de ambos os métodos no mesmo plano cartesiano. Permite visualizar o quão mais rápido (ou não, dependendo da função) a Falsa Posição chega ao limite de tolerância desejado.
