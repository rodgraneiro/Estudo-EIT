# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 11:04:01 2025

@author: rodgr
"""
import numpy as np
import matplotlib.pyplot as plt
#import meshio
from matplotlib.collections import LineCollection
import time


###############################################################################
# Esta função calcula a matriz local de condutividade da malha de elementos
# finitos
#
#        A_i * sigma_i                                      
# Y_local = ------------- * [[1, -1], [-1, 1]]             
#             L                                             
###############################################################################
def calc_Y_local_1D(xl,
                    xm,
                    Area,
                    sigma): # função p/ calc a matriz local p/ elemento 1D
  return ((Area*sigma)/(xm-xl))*np.array([[1, -1], [-1, 1]])
###############################################################################







###############################################################################
# Esta função calcula a matriz global de condutividade da malha de elementos
# finitos
  
###############################################################################
def calc_Y_global_1D(coord,
                  topologia,
                  sigma,
                  area,
                  nodes,
                  N_ele):               # função para calcular a matriz global
  Y_temp = np.zeros((nodes, nodes))
  for i in range(0, N_ele):                 #  laço para montar a matriz global

        node_l = int(topologia[i][0])         # pegar dados da matriz elementos
        node_m = int(topologia[i][1])

        coord_1 = [coord[node_l-1][0]]                # coordenada x do ponto 1
        coord_2 = [coord[node_m-1][0]]                # coordenada x do ponto 2

        Y_local = (
            calc_Y_local_1D(coord_1[0],
                            coord_2[0],
                            area[i],
                            sigma[i])
        )                                              # calcula a matriz local
        #print('Y_local_b',Y_local_b)

        Y_temp[node_l - 1, node_l - 1] += Y_local[0, 0] # monta a matriz global
        Y_temp[node_m - 1, node_l - 1] += Y_local[0, 1]
        Y_temp[node_l - 1, node_m - 1] += Y_local[1, 0]
        Y_temp[node_m - 1, node_m - 1] += Y_local[1, 1]
  return Y_temp
###############################################################################




###############################################################################
# Esta função aplica condições de contorno conforme exemplo abaixo:
#
# [ 1   0          0          0          0        ] [u1]   [u1]    
# [ 0  k1+k2+k3    -k3        0         -k2       ] [u2] = [F2 +k1*u1] 
# [ 0   -k3      k3+k5+k4     -k5        0        ] [u3]   [F3 +k4*u4] 
# [ 0   0          -k5      k2+k6       -k6       ] [u4]   [F4 +k6*u5] 
# [ 0   0           0          0         1        ] [u5]   [u5]
###############################################################################  
def aplica_cond_contorno(vetor_corrente_b,
                         Y_b,
                         n_nodes,
                         V_imposto_b):           # Aplica condições de contorno

  vetor_corrente_cond_contorno = vetor_corrente_b[:].copy()

  for [noh_cond_contorno_b,valor_cond_contorno_b] in V_imposto_b:
    for i in range(0, n_nodes):                    # corrige matriz de corrente
      vetor_corrente_cond_contorno[i] = (
          vetor_corrente_cond_contorno[i] - \
          Y_b[i][noh_cond_contorno_b]*valor_cond_contorno_b
      )

  for [noh_cond_contorno_b,valor_cond_contorno_b] in V_imposto_b:
    vetor_corrente_cond_contorno[noh_cond_contorno_b] = (
        valor_cond_contorno_b )            # coloca valor de contorno conhecido

  Y_cond_contorno = Y_b[:].copy()                        # Criar matriz solução
  for [noh_cond_contorno_b,valor_cond_contorno_b] in V_imposto_b:
    for k in range(0, n_nodes):                # laço para zerar linha e coluna
        Y_cond_contorno[noh_cond_contorno_b][k] = 0
        Y_cond_contorno[k][noh_cond_contorno_b] = 0

    Y_cond_contorno[noh_cond_contorno_b][noh_cond_contorno_b] = 1
  return vetor_corrente_cond_contorno,Y_cond_contorno
###############################################################################




###############################################################################
# Essa função calcula o valor medido conforme equação abaixo
#
# { V_medido } = [ Y_cond_contorno ]^(-1) * { C_cond_contorno }          
###############################################################################
def solucao_Vmedido_b(vetor_correnteVM,
                      Y_cond_contorno_VM):        # calc valor observado/medido
  Yinversa_b = np.linalg.inv(Y_cond_contorno_VM)
  Vmedido_b = np.dot(Yinversa_b,
                     vetor_correnteVM)
  return Vmedido_b                             # retorna valor observado/medido
###############################################################################



###############################################################################
# Essa função calcula o valor estimado
###############################################################################
def solucao_Vcalculado_b(vetor_corrente_VC,
                         Y_cond_contorno_VC):     # calc valor observado/medido

  Yinversa_VC = np.linalg.inv(Y_cond_contorno_VC)

  Vcalc_b = np.dot(Yinversa_VC, vetor_corrente_VC)
  return Vcalc_b
###############################################################################




###############################################################################
# Essa função calcula as derivadas paraciais para montar a matriz Jacobiana p/ 
# método direto de soma e sobreposição das matrizes locais.
# A função pega as coordenadas do elemento atual e calcula a derivada parcial  
# da matriz local atual  por meio da equação:
#
#        ∂                                                              
# ---------- [ Y_local ] = ( A_i / L_i ) * [ [ 1  -1 ], [ -1  1 ] ]    
#      ∂σ_i                                                            
###############################################################################
def calc_Y_jacobiano(coordenadas_c,
                     topologia_c,
                     A_c,
                     nro_nohs_c,
                     nro_elementos_J):   # função para calcular a matriz global
  Y_jacobiano = np.zeros((nro_nohs_c, nro_nohs_c))
  #global comprimento
  #comprimento = np.zeros(nro_elementos_J)
  for i in range(0, nro_elementos_J):         # laço para montar a matriz global

        node_l_J = int(topologia_c[i][0])     # pegar dados da matriz elementos
        node_m_J = int(topologia_c[i][1])

        coord_1_J = [coordenadas_c[node_l_J-1][0]]    # coordenada x do ponto 1
        coord_2_J = [coordenadas_c[node_m_J-1][0]]    # coordenada x do ponto 2


        Y_local_c =(
            ((A_c[i])/(coord_2_J[0]-coord_1_J[0]))*np.array([[1, -1], [-1, 1]])
        )                                    # derivada parcial da matriz local
        #print('Y_local_c',Y_local_c)

        Y_jacobiano[node_l_J - 1, node_l_J - 1] += Y_local_c[0, 0] # monta a
        Y_jacobiano[node_m_J - 1, node_l_J - 1] += Y_local_c[0, 1] # M global
        Y_jacobiano[node_l_J - 1, node_m_J - 1] += Y_local_c[1, 0]
        Y_jacobiano[node_m_J - 1, node_m_J - 1] += Y_local_c[1, 1]
        #comprimento[i] = coord_2_J[0] - coord_1_J[0]  # calc vetor comprimento
        #print('Y_jacobiano',Y_jacobiano)
  return Y_jacobiano
###############################################################################


###############################################################################
# Essa função calcula o comprimento de cada elemento
###############################################################################
def calc_length_elemento(coordenadas_c,
                     topologia_c,
                     nro_elementos_J):   # função para calcular a matriz global
  comprimento = np.zeros(nro_elementos_J)
  for i in range(0, nro_elementos_J):        # laço para montar a matriz global

        node_l_J = int(topologia_c[i][0])     # pegar dados da matriz elementos
        node_m_J = int(topologia_c[i][1])

        coord_1_J = [coordenadas_c[node_l_J-1][0]]    # coordenada x do ponto 1
        coord_2_J = [coordenadas_c[node_m_J-1][0]]    # coordenada x do ponto 2

        comprimento[i] = coord_2_J[0] - coord_1_J[0]   # calc vetor comprimento
  return comprimento
###############################################################################






###############################################################################
# Essa função monta a matriz Jacobiana
###############################################################################
def calc_Jacobiano_b(Y_inv_b, Area_b,
                     comprimento,
                     noh_cond_contorno_b,
                     corrente_b,
                     jacob_ii,
                     J_ii,
                     n_elements, 
                     #no_med=noh_medidos):  # Calc jacobiano em função de sigma
                     no_med):  # Calc jacobiano em função de sigma


  limite = 1e-12                               # limite inferior para jacobiano
  YI_corrente = Y_inv_b @ corrente_b
  for i in range(n_elements):         # loop para montar a matriz jacobiana 1D
    termo_ii = (Area_b[i]/comprimento[i])       # calc termo da matriz da
                                                # derivada parcial do jacobiano

    J_ii[i, i] = termo_ii
    J_ii[i, (i+1)] = -termo_ii
    J_ii[(i+1), i] = -termo_ii
    J_ii[(i+1), (i+1)] = termo_ii

    jacob_ii[:, i] = Y_inv_b @ (J_ii.T @ YI_corrente) # calc matriz jacobiana x
                                                      # vetor de corrente
    jacob_ii[np.abs(jacob_ii) < limite] = 0      # se valor < 1e-12 forçar zero
                                                 # na matriz jacobiana


    J_ii[:] = 0.0                                                # zerar matriz
  return jacob_ii[no_med, :]     # retorna matriz jacobiana x vetor de corrente

###############################################################################






###############################################################################
# Essa função monta a matriz do Filtro Passa Alta
# FPA = M - I
# Onde I é a matriz Identidade e M a matriz gaussiana.
###############################################################################
def calc_L2_gauss_1D(std, centroids_1D):#, covariance_vector):
    """
    Calcula L2 Gaussiana para malha 1D de EIT em uma barra de 1 metro.

    Parâmetros:
    - std: desvio padrão do kernel Gaussiano
    - centroids_1D: vetor (nelements,) com posições dos centroides
    - covariance_vector: vetor (nelements,) com covariância anatômica

    Retorna:
    - F: matriz Filtro Passa Alta
    """
    nelements = len(centroids_1D)
    tol = 1e-9

    F = np.zeros((nelements, nelements), dtype=np.float32)

    for i in range(nelements):
        soma = 0.0
        ci = centroids_1D[i]

        for j in range(nelements):
            cj = centroids_1D[j]
            dist = np.abs(cj - ci)

            if dist <= 1.5 * std:
                fator = 1.0 if i == j else np.exp(-dist**2 / (2 * std**2))
                soma += fator
                F[i, j] = fator

        for j in range(nelements):
            if i == j:
                aux = (1.0 - F[i, j] / soma)
                F[i, j] = aux if np.abs(aux) > tol else 0.0
            else:
                aux = (-F[i, j] / soma)
                F[i, j] = aux if np.abs(aux) > tol else 0.0

    return F
###############################################################################




###############################################################################
# Essa função calcula a variação de sigma estimado para a próxima iteração
#
# Δ = α_k * [ (J_kᵗ W₁ J_k + λ² L₂ᵗ L₂ )⁻¹ ] *                                             
#         [ J_kᵗ W₁ ( z - h(θ̂_k) )                                                        
#           - λ² L₂ᵗ L₂ ( θ̂_k - θ* ) ]                                                    
###############################################################################
def calc_delta_sigma_b(J_b,
                       alpha_b,
                       lambda_b,
                       residuo_b,
                       L2,
                       sigma_hat, 
                       sig_estrela):  # calc  correção delta próxima iteração
                       #sigma_hat = sigma_inicial_b,
                       #sig_estrela = chute):  # calc  correção delta próxima iteração
    zW1=np.eye(J_b.shape[0])
    JT = J_b.T
    zJTW = JT @ zW1
    zJTWJ = zJTW @ J_b
    zLTL = L2.T @ L2
    termo_L = (lambda_b**2) * zLTL
    primeiroTermo = zJTWJ + termo_L
    inv_primeiroTermo = np.linalg.inv(primeiroTermo)

    JTW_zh = zJTW @ residuo_b
    ztermo_reg = (sigma_hat - sig_estrela)
    zregularizacao = (lambda_b**2)*zLTL @ ztermo_reg
    segundoTermo = JTW_zh - zregularizacao
    return -alpha_b*(inv_primeiroTermo @ segundoTermo)
###############################################################################



###############################################################################
# Essa função plota o gráfico convergência das iterações
###############################################################################

def plotar_iteracoes(lista_indice, lista_valor):
    plt.plot(lista_indice,
            lista_valor,
            marker='.',
            linestyle='-',
            color='b',
            label='Norma $\Delta\sigma$')        # plota gráfico das iterações,
    plt.xlabel("Iteração", fontsize=15)
    plt.ylabel("Norma Delta_Sigma", fontsize=15)
    plt.title("Otimização da Barra de Cobre", fontsize=15)
    plt.legend()
    plt.grid(True)
    plt.show()
###############################################################################




###############################################################################
# Essa função plota a curva de valores reais de sigma versus valores de sigma
# calculados nas iterações
###############################################################################
def plotar_grafico(matriz_coordenadas_b,
                   matriz_topologia_b,
                   sigma_inicial_b,
                   sigma_real_b,
                   noh_medidos,
                   alpha_b,
                   lambda_b,
                   std
                   ):
    # Coordenadas x dos nós
    x_coords = matriz_coordenadas_b.flatten()
    topologia_b = matriz_topologia_b.astype(int) - 1

    centros = np.mean(x_coords[topologia_b], axis=1)
    #centros = (x_coords[:-1] + x_coords[1:]) / 2
    valores = sigma_inicial_b  # ou delta_sig_b
    valores_real = sigma_real_b  # ou delta_sig_b
    noh_medidos[-1] -= 1
    pos_medidas = [centros[i] for i in noh_medidos]
    pos_valores_med = [valores[i] for i in noh_medidos]
    pos_valores_real = [valores_real[i] for i in noh_medidos]
    # 4. Gráfico tipo steam
    plt.figure(figsize=(16, 4))
    #plt.stem(centros, valores)
    plt.plot(centros, valores ,
            marker='None',
            label='$\sigma$ calculado')
    plt.plot(centros, sigma_real_b,
            marker='None',
            linestyle=':',
            label='$\sigma$ real' )
    plt.plot(pos_medidas,
            pos_valores_real,
            marker='x',
            linestyle='None',
            markersize=10,
            color='red',
            label='Ptos medidos')
    plt.xlim(0, 1.01)
    plt.ylim(0.0, 0.60)
    plt.xlabel('Posição [m]')
    plt.ylabel('Condutividade σ')
    plt.title('Distribuição de  σ nos elementos 1D')
    plt.text(0.5, 0.87,
            f'Para α = {alpha_b} , λ = {lambda_b} e std={std}',
            ha='center',
            va='center',
            transform=plt.gcf().transFigure)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()
###############################################################################







































