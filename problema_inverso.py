# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 16:04:22 2025

@author: rodgr
"""
import numpy as np

from EIT_functions import calc_Y_local_1D
from EIT_functions import calc_Y_global
from EIT_functions import aplica_cond_contorno
from EIT_functions import solucao_Vmedido_b
from EIT_functions import solucao_Vcalculado_b
from EIT_functions import calc_Y_jacobiano
from EIT_functions import calc_L2_gauss_1D
from EIT_functions import calc_Jacobiano_b
from EIT_functions import calc_delta_sigma_b
from EIT_functions import calc_length_elemento
from EIT_functions import plotar_grafico

from dados_conhecidos import (n_iter, noh_medidos, noh_cond_contorno_b, jacob_ii, 
                              J_ii, alpha_b, lambda_b, chute, lista_i, lista_plotar, std)
                              #sigma_inicial_b, std)
from EIT_1D import (matriz_coordenadas_b,
              matriz_topologia_b,
              Area_b,
              n_nodes,
              n_elements,
              vetor_corrente,
              V_imposto_b,
              Vmedido_b,
              Y_cond_contorno_J,
              comprimento,
              vetor_corrente_cond_contorno_J,
              residuo_b,
              L2,
              sigma_real_b)

def solve_inverse_problem_1D(sigma_inicial_b):
     # LOOP PRINCIPAL
    for i in range(n_iter):
      Y_Vcalc_b = (
          calc_Y_global(matriz_coordenadas_b,
                        matriz_topologia_b,
                        sigma_inicial_b,
                        Area_b,
                        n_nodes,
                        n_elements,)
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
      
    plotar_grafico(matriz_coordenadas_b,
                       matriz_topologia_b,
                       sigma_inicial_b,
                       sigma_real_b,
                       noh_medidos,
                       alpha_b,
                       lambda_b,
                       std)
      
    return Sig_kMais1_b   