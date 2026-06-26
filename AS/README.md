# Atividades de Sala (AS)

Esta pasta contém pequenos programas desenvolvidos durante as aulas para praticar a implementação de métodos numéricos fundamentais em Fortran.

## Atividades e Conteúdos

- **AS 1 (`AS/1/atividade1.f90`)**: Implementação dos métodos de confinamento para busca de raízes: **Bissecção** e **Falsa Posição**. Compara a convergência de ambos os métodos na solução de uma equação de transferência de calor.
- **AS 2 (`AS/2/newton_raphson.f90`)**: Implementação do método aberto de **Newton-Raphson** (Padrão e Modificado). Utiliza a derivada analítica para buscar raízes com convergência muito mais rápida. O método modificado estabiliza problemas de raízes múltiplas.
- **AS 3 (`AS/3/muller.f90`)**: Implementação do **Método de Muller**. Técnica que utiliza uma interpolação quadrática em três pontos (parábolas em vez de retas) para encontrar raízes de funções polinomiais com alta eficiência, sem precisar de derivadas analíticas.
- **AS 4 (`AS/4/decomp_LU.f90`)**: Implementação da **Fatoração (Decomposição) LU** para a resolução de sistemas de equações lineares. Decompõe a matriz do sistema $AX=B$ nas matrizes Lower (L) e Upper (U), permitindo a solução ágil de múltiplos vetores de forças por meio de substituições progressivas e regressivas.
- **AS 5 (`AS/5/inversaLU.f90`)**: Expansão do algoritmo de Fatoração LU. Além de resolver o sistema, o programa aproveita a matriz fatorada para determinar computacionalmente a **Matriz Inversa ($A^{-1}$)**.
- **AS 6 (`AS/6/interpolacao.f90` e `AS/6/aurea.f90`)**: Programas focados em **Otimização Unidimensional** (Busca do máximo ou mínimo de uma função):
  - *Interpolação Quadrática* (quadratica): Estima analiticamente o vértice do cume/vale de parábolas sucessivas.
  - *Razão Áurea* (aurea): Reduz gradualmente o intervalo de busca dividindo o domínio na proporção de ouro.

## Como utilizar

A maioria dos algoritmos exporta dados em tabelas interativas no terminal e paralelamente salva logs de convergência em arquivos de dados (`.dat`), os quais podem ser lidos para plotagem de gráficos com o gnuplot ou python.

1. Navegue até a atividade desejada:
   ```bash
   cd 1
   ```
2. Compile com o gfortran:
   ```bash
   gfortran nome_do_arquivo.f90 -o exec_as
   ```
3. Execute o binário gerado:
   ```bash
   ./exec_as
   ```
