# AS 6: Otimização Unidimensional

## Sobre os Códigos
Nesta pasta, estão implementados dois algoritmos voltados à Otimização Unidimensional: a busca por pontos críticos (máximos ou mínimos) de uma função sem o uso de derivadas analíticas. Os códigos são `interpolacao.f90` (Método Quadrático) e `aurea.f90` (Método da Razão Áurea).
Eles resolvem um problema elétrico em comum: Qual o valor da resistência $Ra$ do circuito que **maximiza** a dissipação de potência elétrica?

### Matemática Envolvida

#### Interpolação Quadrática Sucessiva (`interpolacao.f90`)
Este método requer três chutes iniciais ($x_0, x_1, x_2$) que cercam a área do pico.
O algoritmo calcula o polinômio quadrático (parábola) que passa por esses três pontos e, com os coeficientes encontrados, ele descobre matematicamente a exata coordenada $x_3$ do vértice da parábola.
Com o novo ponto $x_3$, o algoritmo descarta o pior ponto antigo e redesenha uma nova parábola nos novos limites. Ele repete isso até o vértice estagnar no pico real da curva da função.

#### Razão Áurea (`aurea.f90`)
Este método atua encurtando o intervalo iterativamente na infame "Proporção Áurea" ($\approx 0.618$).
Dado o intervalo de busca inferior e superior ($x_{lower}, x_{upper}$):
1. Ele cria dois pontos de teste simétricos dentro do intervalo baseados na distância da proporção áurea.
2. Avalia a função em ambos os pontos internos.
3. Se o pico estiver visivelmente do lado do ponto esquerdo, o programa "corta" e joga fora a janela que estava após o ponto direito, definindo um novo limite superior no lugar, e vice-versa.
Com isso, ele vai encurralando o cume da montanha em um espaço minúsculo até atingir a precisão.

A função avaliada é a lei de potência de um circuito paramétrico:
$$P(R_a) = \frac{\left[ \frac{V \cdot R_3 \cdot R_a}{R_1(R_a+R_2+R_3) + R_3 R_a + R_3 R_2} \right]^2}{R_a}$$

### Entradas e Saídas
- **Entradas:** As constantes do circuito ($V, R_1, R_2, R_3$) e a tolerância de parada de $1.0\times 10^{-4}$.
  - Para a Quadrática: Os três chutes da resistência $Ra_0, Ra_1, Ra_2$.
  - Para a Áurea: Os dois limites extremos $Ra_{lower}$ e $Ra_{upper}$.
- **Saídas (Terminal):** Ambos os algoritmos imprimem em tempo real o histórico do encurtamento da janela de busca e as estimativas do novo ponto. No fim de ambos, são exibidos o **Valor Ótimo** na coordenada X (A resistência ideal) e o **Ponto Máximo** na coordenada Y (A potência dissipada resultante).
