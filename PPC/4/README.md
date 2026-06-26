# PPC 4: Otimização Multidimensional Topográfica

## Sobre o Código
O programa `PPC4.f90` tem um escopo de exploração e otimização. Ele projeta um algoritmo capaz de "subir" uma montanha matemática descrita pela função bidimensional $f(x,y)$, buscando encontrar, iterativamente, o ponto mais alto do terreno (o Ponto de Máximo Global).

O código é uma demonstração de força comparando duas inteligências de direção de passo: o método do **Aclive Máximo (Steepest Ascent)** e o algoritmo de **Fletcher-Reeves (Gradientes Conjugados)**.

### Matemática Envolvida
Independente do método de direção, toda Otimização Multidimensional segue a fórmula base iterativa:
$$\vec{r}_{k+1} = \vec{r}_k + h^* \cdot \vec{d}_k$$
- $\vec{r}_k$: Coordenadas atuais $(X, Y)$.
- $\vec{d}_k$: O vetor Direção de busca.
- $h^*$: O tamanho ideal do passo (o quão longe pular naquela direção).

**1. O Vetor Direção ($\vec{d}_k$):**
- **Aclive Máximo:** Sempre escolhe o vetor do Gradiente analítico ($\vec{\nabla}f$). Isso significa que o algoritmo aponta estritamente para a rampa mais íngreme naquele instante. Apesar de lógico, isso força o caminho a fazer curvas de 90 graus em vales, criando uma caminhada em "zigue-zague" extremamente devagar perto do cume.
- **Fletcher-Reeves:** Acopla na equação do Gradiente um parâmetro de "Inércia" ($\beta$). Ele "lembra" o caminho da iteração passada e conjuga esse vetor antigo com a ladeira atual, "cortando caminho" nas curvas e atacando diretamente o centro da montanha.

**2. A Busca Linear ($h^*$):**
Uma vez escolhida a direção, para não dar um pulo no vazio ou além da montanha, o código precisa descobrir o passo ótimo ($h^*$). Isso é feito reaproveitando o **Método da Interpolação Quadrática** embutido dentro de uma sub-rotina do código: Ele cria uma parábola na linha da direção traçada, e define que o pulo exato vai terminar na coordenada X do vértice dessa parábola.

### Entradas e Saídas
- **Entradas (Menu Interativo):** O console solicitará que se digite `'1'` para rodar com o Aclive Máximo, ou `'2'` para Fletcher-Reeves.
- **Entradas (Hardcoded):** Ponto inicial (Chute) $X = -2.0, Y = 3.0$. O método encerra quando o tamanho do vetor gradiente cair para $\approx 0$ (o terreno ficou perfeitamente plano - o cume).
- **Saídas (Terminal):** Histórico da caminhada em colunas: Iteração, Erro, Pulo $h^*$, coordenadas atuais ($X, Y$) e os deltas.
- **Saídas (Arquivos):** Grava a "Trilha" do algoritmo no arquivo log correspondente (`output1.dat` ou `output2.dat`). Além disso, varre e exporta também as curvas de nível da montanha inteira no arquivo `function.dat`.

### Gráficos (Python e Matplotlib)
O script `graficos.py` importa todos esses dados para a biblioteca Matplotlib e desenha uma visualização topográfica esplêndida (`comparacao_matplotlib.png`):
- Ele constrói uma visão superior da montanha (Curvas de Nível / Contour Map).
- Ele plota, ponto a ponto, as linhas do caminho trilhado pela Inteligência de "Aclive Máximo" (em preto) contra a "Fletcher-Reeves" (em vermelho) partindo do ponto inicial $(-2, 3)$ em direção ao pico branco em $(2, 1)$, provando de forma esmagadora o traçado muito mais direto e eficiente do algoritmo dos Gradientes Conjugados.
