# Otimização Multidimensional Irrestrita em Fortran: Aclive Máximo vs Fletcher-Reeves

Este projeto implementa e compara dois algoritmos clássicos de otimização irrestrita para encontrar o ponto máximo de uma função matemática de duas variáveis: $f(x,y)$. O código foi construído em Fortran e modularizado para evidenciar a arquitetura do método.

## A Função Objetivo

A "montanha" matemática que desejamos escalar é definida pela função:

$$f(x,y) = 2xy + 2x - x^2 - 2y^2$$

O ótimo analítico desta função encontra-se nas coordenadas **(2, 1)**. Iniciando no ponto **(-2, 3)**, testamos a eficiência de dois métodos numéricos distintos para alcançar o pico. O ponto foi escolhido de acordo com a atividade proposta para a Atividade Para Casa 4 (APC4), mas pode ser alterado facilmente.

---

## O Motor Comum: Busca Linear

Independentemente do método escolhido para decidir *para onde* ir, precisamos saber *o quanto* andar. Para isso, o código utiliza uma busca unidimensional (Busca Linear). 

A cada iteração, avaliamos o terreno em três pontos ao longo da reta de direção escolhida, traçamos uma parábola virtual através deles e calculamos analiticamente o vértice dessa parábola ($h^*$). Isso garante que cada passo seja da magnitude ideal. Realizamos essa busca com o uso da interpolação quadrática, método utilizado em uma das atividades de sala recentes para encontrar o ponto otimizado de uma função de resistência elétrica.

---

## Os Métodos Comparados

### 1. O Método do Aclive Máximo (Steepest Ascent)
**A Lógica:**
É o método mais cru. A cada passo, o algoritmo calcula o gradiente ($\nabla f$) e vira exatamente para a direção da subida mais íngreme. 

**O Desempenho:**
É um método relativamente cego, já que a cada iteração, a nova direção de subida é ortogonal a anterior, realizando uma busca em **ziguezague**. Ele fica batendo nas laterais das curvas de nível, necessitando de uma quantidade significativa de iterações até chegar ao topo de $f(x,y)$. Note a quantidade de iterações é diretamente relacionada com a tolerância definida, com tolerâncias mais baixas resultando em resultados mais imprecisos, mas que convergem mais rápido. Em cinco iterações os valores já estavam próximos, mas para o valor utilizado no código anexado, mais 40 iterações foram realizadas até o ponto final. 

### 2. O Método de Fletcher-Reeves (Gradientes Conjugados)
**A Lógica:**
É uma evolução inteligente do Aclive Máximo. O Fletcher-Reeves não olha apenas para o gradiente atual, ele carrega uma "memória" da inércia do passo anterior (através do parâmetro $\beta_k$). Isso permite que ele ajuste a trajetória, cortando caminho pelas diagonais da função em vez de virar abruptamente 90 graus.

**O Desempenho:**
Em funções puramente quadráticas (como a do nosso projeto), a teoria das Direções Conjugadas garante que o método encontrará a solução em no máximo $n$ passos, onde $n$ é o número de dimensões. Como o nosso problema tem duas variáveis (X e Y), o Fletcher-Reeves resolve em **apenas 2 iterações**, indo quase em linha reta para o topo (2, 1).

---

## Estrutura dos Arquivos Gerados

Após a execução, o programa Fortran exporta arquivos `.dat` estruturados nativamente para leitura por softwares de plotagem (como Gnuplot, Matplotlib/Python ou GNU Octave):

* `function.dat`: Contém a malha tridimensional da topografia (curvas de nível).
* `output1.dat`: Log completo com as 46 iterações do Aclive Máximo.
* `output2.dat`: Log completo com as 2 iterações do Fletcher-Reeves.

Ao plotar os dados em conjunto, a diferença de eficiência entre as trajetórias (Ziguezague vs Direto) fica visualmente clara sobre o mapa de calor. O mapa de calor anexado foi realizado com o Matplotlib, já que enfrentei problemas para fazer o código funcionar em gnuplot. Apesar de gostar mais do gnuplot, a criação de arquivos .dat e o uso de bibliotecas em python para a plotagem dos gráficos é altamente recomendada.
