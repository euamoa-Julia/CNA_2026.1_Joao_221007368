import os
import numpy as np
import matplotlib.pyplot as plt

pasta_script = os.path.dirname(os.path.abspath(__file__))

caminho_function = os.path.join(pasta_script, 'function.dat')
caminho_output1 = os.path.join(pasta_script, 'output1.dat')
caminho_output2 = os.path.join(pasta_script, 'output2.dat')

malha = np.genfromtxt(caminho_function)
x_grid = malha[:, 0]
y_grid = malha[:, 1]
z_grid = malha[:, 2]

aclive = np.genfromtxt(caminho_output1)
reeves = np.genfromtxt(caminho_output2)

x_aclive = np.insert(aclive[:, 3], 0, -2.0)
y_aclive = np.insert(aclive[:, 4], 0, 3.0)

x_reeves = np.insert(reeves[:, 3], 0, -2.0)
y_reeves = np.insert(reeves[:, 4], 0, 3.0)

plt.figure(figsize=(10, 8))

contorno = plt.tricontourf(x_grid, y_grid, z_grid, levels=20, cmap='magma')
plt.colorbar(contorno, label='Valor da Função Objetivo f(x, y)')

plt.plot(x_aclive, y_aclive, color='black', marker='o', linestyle='-', 
         linewidth=1.2, markersize=3, label='Aclive Máximo')

plt.plot(x_reeves, y_reeves, color='red', marker='s', linestyle='-', 
         linewidth=1.2, markersize=3, label='Fletcher-Reeves')

plt.plot(-2, 3, 'bo', markersize=7, markeredgecolor='black', label='Ponto Inicial (-2, 3)')
plt.plot(2, 1, 'wo', markersize=7, markeredgecolor='black', label='Ótimo Analítico (2, 1)')

plt.title('Otimização: Aclive Máximo vs Fletcher-Reeves', fontsize=16, pad=15)
plt.xlabel('Eixo X', fontsize=12)
plt.ylabel('Eixo Y', fontsize=12)
plt.grid(color='white', linestyle='--', alpha=0.3)
plt.legend(loc='lower right', fontsize=11, framealpha=0.9)

plt.xlim(-3.0, 3.0)
plt.ylim(-1.0, 4.0)

caminho_imagem = os.path.join(pasta_script, 'comparacao_matplotlib.png')
plt.savefig(caminho_imagem, dpi=300, bbox_inches='tight')
