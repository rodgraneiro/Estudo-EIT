# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 22:42:12 2025

@author: rodgr
"""
from EIT_functions import (calc_Y_global_1D, solucao_Vmedido_b)
from EIT_functions import aplica_cond_contorno
from EIT_1D import(matriz_coordenadas_b, matriz_topologia_b, sigma_real_b,
                               Area_b, n_nodes, n_elements, vetor_corrente, V_imposto_b, noh_medidos)

###############################################################################
########################## PROBLEMA DIRETO ####################################
###############################################################################
def solve_forward_problem_1D():
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
    return Vmedido
    #print('Vmedido \n',Vmedido)
    
    #Vmedido_b = Vmedido[noh_medidos]
    #print('matriz_coordenadas_b',matriz_coordenadas_b)
    #print('matriz_topologia_b',matriz_topologia_b)
    #print('Vmedido (eletrotos) \n',Vmedido_b)
###############################################################################
###############################################################################