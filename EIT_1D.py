# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 11:06:49 2025

@author: rodgr
"""

import numpy as np
import matplotlib.pyplot as plt
import meshio
from matplotlib.collections import LineCollection
import time

from EIT_functions import calc_Y_local_1D
from EIT_functions import calc_Y_global_1D
from EIT_functions import aplica_cond_contorno
from EIT_functions import solucao_Vmedido_b
from EIT_functions import solucao_Vcalculado_b
from EIT_functions import calc_Y_jacobiano
from EIT_functions import calc_L2_gauss_1D
from EIT_functions import calc_Jacobiano_b
from EIT_functions import calc_delta_sigma_b
from EIT_functions import calc_length_elemento
from EIT_functions import plotar_grafico
from EIT_functions import plotar_iteracoes
from EIT_functions import calc_L2_Adler1
from EIT_functions import calc_L2_Adler2
from EIT_functions import calc_Gaussian_HPF_1D
from EIT_mesh_1D import files_mesh
from EIT_Plot_1D import plot_EIT


inicio = time.time()





###############################################################################





###############################################################################
########################## CARREGA MALHA 1D ###################################
###############################################################################

#opcao = input('Escolha a malha(71, 100, 200, 300, 1000): ')
opcao = '100'

if opcao not in ['71','100', '100', '200', '300', '1000']:
    raise ValueError("Opção inválida.")

caminho, sigma_real_b, noh_med = files_mesh(opcao)

malha_msh = meshio.read(caminho)                            # Lê o arquivo .msh

matriz_coordenadas_b = malha_msh.points[:, :1]       # Monta matriz coordenadas
matriz_topologia_b = malha_msh.cells_dict["line"]      # Monta Matriz topologia
matriz_topologia_b = matriz_topologia_b +1

n_nodes = matriz_coordenadas_b.shape[0]                             # Nr de nós
n_elements = matriz_topologia_b.shape[0]                      # Nr de elementos


###############################################################################
########################## DADOS INICIAIS #####################################
###############################################################################

alpha_b = 0.010                               # inicia variavel de ajuste alpha
lambda_b = 0.010                                # inicia variavel de ajuste alpha
max_iter = 200                                         # Nr máximo de iterações
lista_i = []                                        # Lista armazenar iterações
lista_plotar = []                                      # Lista Valores de sigma
convergencia = False                       # Inicia varriável para convergência
#n =n_elements
std = 0.1
noh_cond_contorno_b = 0           # nó de referência para condições de contorno
noh_med_op =  noh_med                     # Vetor elementos medidos (eletrodos)
noh_medidos = noh_med_op.copy()           # Vetor elementos medidos (eletrodos)
#noh_medidos = [0,  10, 20, 30, 40, 50, 60, 70, 80, 90,  n_elements]
#noh_medidos = [0,  7, 14, 21, 28, 35, 42, 49, 56, 63,  n_elements]      

# Solicita ao usuário os nós medidos
#entrada = input("Digite os nós medidos separados por vírgula (ex.: 0,10,20: ")

# Converte a string em uma lista de inteiros
#noh_medidos = [int(x.strip()) for x in entrada.split(',')]


sigma_inicial_b = np.full(n_elements, 1.0)          # Monta vetor sigma inicial
chute = np.full(n_elements, 1.0)                           # Monta vetor chute
Area_b = np.full(n_elements, 0.001)                             # área da barra
vetor_corrente = np.zeros(n_nodes)                    # Monta vetor de corrente
vetor_corrente[0] = -0.001                            # Nó de saída de corrente
vetor_corrente[n_elements] = 0.001                  # Nó de entrada de corrente



# Solicita ao usuário os nós medidos
#entrada = input("Digite os nós medidos separados por vírgula (ex.: 0,10,20: ")

# Converte a string em uma lista de inteiros
#noh_medidos = [int(x.strip()) for x in entrada.split(',')]
#print("Nó(s) medido(s):", noh_medidos)



jacob_ii = np.zeros((n_nodes, n_elements))   # inicia a variável 1 do jacobiano
J_ii = np.zeros((n_nodes, n_nodes))          # inicia a variável 2 do jacobiano

V_imposto_b = [[0, 0.0]]              # no formato [nó,valor], ou seja, V_2 = 5

x_coords_b = matriz_coordenadas_b.flatten()
topologia_bc = matriz_topologia_b.astype(int) - 1

centroids_1D = np.mean(x_coords_b[topologia_bc], axis=1)

mdl_dim =1.0                                         # comprimento total 1 metro
diam_frac = 0.30*mdl_dim                               # 30% co comprimrnto total
beta = 2.769 / (diam_frac * mdl_dim)**2
n_points = 3
s_k = np.linspace(0, 1, n_points)

#dx = 1.0 / nelements
#centroids_1D = np.linspace(dx/2, 1 - dx/2, nelements)
covariance_vector = np.ones(n_elements) * 0.001

comprimento = calc_length_elemento(matriz_coordenadas_b, 
                                   matriz_topologia_b, n_elements)

#L2 = calc_L2_gauss_1D(std, centroids_1D) #, covariance_vector)
#L2 = calc_L2_Adler1(diam_frac, mdl_dim, beta, comprimento, matriz_coordenadas_b, 
#                 centroids_1D, s_k, n_elements, limite=1e-4, n_points=n_points)

#L2 = calc_L2_Adler2(diam_frac, mdl_dim, beta, comprimento, matriz_coordenadas_b, 
#                 centroids_1D, s_k, n_elements, limite=1e-4, n_points=n_points)

L2 = calc_Gaussian_HPF_1D(x_coords_b, topologia_bc, diam_frac)
###############################################################################





###############################################################################
########################## PROBLEMA DIRETO ####################################
###############################################################################

Y_global_med_b = calc_Y_global_1D(matriz_coordenadas_b,
                               matriz_topologia_b,
                               sigma_real_b,
                               Area_b,
                               n_nodes,
                               n_elements)    # calc M Y ajustada cond contorno

vetor_corrente_cond_contorno_b,Y_cond_contorno_b = (
    aplica_cond_contorno(vetor_corrente,
                         Y_global_med_b,
                         n_nodes,
                         V_imposto_b)
)                              # calcula vetor corrente ajustado cond contorno

Vmedido = (
    solucao_Vmedido_b(vetor_corrente_cond_contorno_b,
                      Y_cond_contorno_b)
)                                                      # calcula valor medido



#print('Vmedido \n',Vmedido)

Vmedido_b = Vmedido[noh_medidos]
#print('matriz_coordenadas_b',matriz_coordenadas_b)
#print('matriz_topologia_b',matriz_topologia_b)
print('Vmedido (eletrotos) \n',Vmedido_b)
###############################################################################
###############################################################################





###############################################################################
########################## PROBLEMA INVERSO ####################################
###############################################################################
# %%time

Y_jacobiano = (
    calc_Y_jacobiano(matriz_coordenadas_b,
                     matriz_topologia_b,
                     Area_b,
                     n_nodes,
                     n_elements)
)                                                       # calc matriz Jacobiana

vetor_corrente_cond_contorno_J,Y_cond_contorno_J = (
    aplica_cond_contorno(vetor_corrente,
                         Y_jacobiano,n_nodes,
                         V_imposto_b)
)                                     # Aplica condições de contorno no Y atual

#np.savetxt('/content/drive/MyDrive/Colab Notebooks/Y_jacobiano_1000e.txt', \
#           Y_jacobiano)





ini_loopFor = time.time()              #########################################
 # LOOP PRINCIPAL

for i in range(max_iter):
  Y_Vcalc_b = (
      calc_Y_global_1D(matriz_coordenadas_b,
                    matriz_topologia_b,
                    sigma_inicial_b,
                    Area_b,
                    n_nodes,
                    n_elements)
  )                                   # cal Y para novo sigma da iteração atual

  vetor_corrente_cond_contorno_c,Y_cond_contorno_c = (
      aplica_cond_contorno(vetor_corrente,
                           Y_Vcalc_b,
                           n_nodes,
                           V_imposto_b)
          )                           # Aplica condições de contorno no Y atual

  V_calc= (
      solucao_Vcalculado_b(vetor_corrente_cond_contorno_c,
                           Y_cond_contorno_c)
  )                                    # calcula vetor tensão da iteração atual
  V_calc_b = V_calc[noh_medidos]


  residuo_b = Vmedido_b - V_calc_b        # calc dif entre Vmedido e VCalculado

  inv_Y_cond_contorno_J = np.linalg.inv(Y_cond_contorno_J)

  jacobiano_b = (
      calc_Jacobiano_b(inv_Y_cond_contorno_J,
                       Area_b, comprimento,
                       noh_cond_contorno_b,
                       vetor_corrente_cond_contorno_J,
                       jacob_ii,
                       J_ii,
                       n_elements,
                       noh_medidos)
  )                                                       # calcula o jacobiano


  delta_sig = (calc_delta_sigma_b(jacobiano_b,
                                  alpha_b,
                                  lambda_b,
                                  residuo_b,
                                  L2,
                                  sigma_inicial_b,
                                  chute)
  )                    # calcula a diferença entre sigma atual e sigma desejado

  plotar = np.linalg.norm(delta_sig)  # calc norma vetor delta_sigama para plot
  lista_i.append(i)                             # Armazena o índice da iteração
  lista_plotar.append(plotar)                  # Armazena o valor a ser plotado

  if np.linalg.norm(delta_sig) < 1e-3:    # Convergência atingida se a norma de
                                          # delta_sigam < que  1e-6
    print(f'Convergência atingida após {i+1} iterações.')
    #print('Vmedido_b \n', Vmedido_b)
    #print('Valor calculado \n',V_calc_b)
    print('Sigma k+1 \n', sigma_inicial_b)
    convergencia = True
    break                                   # interrompe o processo de iteração

  Sig_kMais1_b = sigma_inicial_b + delta_sig           # ajusta vetor sigma k+1
  sigma_inicial_b = Sig_kMais1_b            # armazena vetor sigma k+1 anterior




fim_loopFor = time.time() - ini_loopFor ########################################
fim = time.time()

#Salva vetor sigma_inicial_b
#np.savetxt('/content/drive/MyDrive/Colab Notebooks/curva_std_001_lambda_100_alfa_01.txt',  sigma_inicial_b)

print(f'Tempo de execução TOTAL: {fim - inicio} segundos')
#print('V_calc_b \n', V_calc_b)
#print(f"Vmedido: {Vmedido_b}")
#sigma_inicial_b = np.full(n_elements, 1.0)          # Monta vetor chute inicial
#chute = np.full(n_elements, 1.0)

plotar_iteracoes(lista_i, lista_plotar)

plotar_grafico(matriz_coordenadas_b,
                   matriz_topologia_b,
                   sigma_inicial_b,
                   sigma_real_b,
                   noh_medidos,
                   alpha_b,
                   lambda_b,
                   std)
print('comprimento', comprimento)
plot_EIT(matriz_coordenadas_b,
                   matriz_topologia_b,
                   sigma_inicial_b,
                   sigma_real_b,
                   noh_medidos,
                   alpha_b,
                   lambda_b,
                   std)

print('matriz_coordenadas_b', matriz_coordenadas_b.shape)
print('sigma_inicial_b', sigma_inicial_b.shape)
print('matriz_topologia_b', matriz_topologia_b.shape)




sigma_inicial_b = np.full(n_elements, 1.0)          # Monta vetor chute inicial
chute = np.full(n_elements, 1.0)
#noh_medidos = [0,  10, 20, 30, 40, 50, 60, 70, 80, 90,  n_elements]    # Vetor elementos medidos (eletrodos)
###############################################################################
###############################################################################


































