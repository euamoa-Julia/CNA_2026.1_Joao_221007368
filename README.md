# Cálculo Numérico Aplicado - Universidade de Brasília (UnB)

Bem-vindo(a) ao repositório da disciplina de Cálculo Numérico Aplicado da Universidade de Brasília (UnB).
Este repositório contém os códigos e programas desenvolvidos ao longo do curso, organizados em duas categorias principais: Atividades de Sala (AS) e Programas Para Casa (PPC).

## Estrutura do Repositório (Branches e Pastas)

O repositório é estruturado em duas pastas principais:

- **AS (Atividades de Sala)**: Contém todos os códigos desenvolvidos em sala de aula. Estes são programas menores, geralmente focados em aprender e aplicar métodos numéricos específicos abordados durante as aulas teóricas.
- **PPC (Programas Para Casa)**: Contém os códigos relacionados aos trabalhos e projetos extraclasse (Assignments). Estes são programas mais completos e robustos, que exigem a aplicação de múltiplos conceitos, análises de desempenho, geração de gráficos e interpretação física dos resultados.

Cada pasta (AS e PPC) possui uma numeração sequencial (1, 2, 3...) correspondente a cada atividade. Para obter mais detalhes sobre o que cada atividade ou programa aborda, visite os arquivos `README.md` específicos dentro de cada pasta.

## Como utilizar os códigos

A grande maioria dos algoritmos deste repositório foi desenvolvida na linguagem **Fortran 90 (`.f90`)**.
Para compilar e executar os códigos em seu próprio ambiente, siga os passos abaixo:

### Pré-requisitos
Certifique-se de ter um compilador Fortran instalado. O mais comum é o `gfortran`.
- No Ubuntu/Debian: `sudo apt install gfortran`
- No macOS (via Homebrew): `brew install gcc`
- No Windows: Pode-se utilizar o [MinGW](https://www.msys2.org/) ou o [Windows Subsystem for Linux (WSL)](https://docs.microsoft.com/pt-br/windows/wsl/install).

### Compilação e Execução

1. Navegue pelo terminal até o diretório do código que deseja executar:
   ```bash
   cd AS/1
   ```

2. Compile o código Fortran utilizando o `gfortran`:
   ```bash
   gfortran nome_do_arquivo.f90 -o executavel
   ```

3. Execute o programa gerado:
   ```bash
   ./executavel
   ```

*(Nota: Alguns programas podem solicitar a entrada de parâmetros pelo usuário via terminal (stdin) e outros geram arquivos de dados (`.dat`) que podem ser usados para plotagem com ferramentas como o `gnuplot` ou `Python/Matplotlib`.)*

## Conteúdo

### AS (Atividades de Sala)
Explore a pasta [`AS/`](AS/README.md) para encontrar implementações de métodos numéricos fundamentais:
- Busca de raízes (Bissecção, Falsa Posição, Newton-Raphson, Muller).
- Resolução de sistemas lineares (Fatoração LU).
- Otimização Unidimensional (Interpolação Quadrática, Busca Áurea).

### PPC (Programas Para Casa)
Explore a pasta [`PPC/`](PPC/README.md) para visualizar os projetos mais avançados:
- Resolução de Equações Diferenciais Ordinárias (EDOs) com Runge-Kutta de 4ª ordem.
- Extração de raízes complexas de polinômios via Método de Bairstow e geração de fractais.
- Simulação de reatores nucleares (diferenças finitas e esquemas implícitos/Algoritmo de Thomas).
- Otimização Multidimensional com Gradientes Conjugados e Aclive Máximo.

---
*Este repositório serve como um portfólio e um registro dos estudos e implementações de algoritmos matemáticos aplicados à engenharia e ciências exatas.*
