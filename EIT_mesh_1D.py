# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 15:11:48 2025

@author: rodgr
"""
import numpy as np

def files_mesh(opcao):
    caminho_malhas = {
    '71' : 'unidimensional_71e_py.msh',     # malha elem de tamanhos diferentes    
    # '100' : 'unidimensional_100e_py.msh',
    '100' : 'unidimensional_100T_dif_py.msh',# malha ele de tamanhos diferentes
    '200' : 'unidimensional_200e_py.msh',
    '300' : 'unidimensional_300e_py.msh',
    '1000' : 'unidimensional_1000e_py.msh',
    }

    sigma_reais = {
        '71': np.concatenate([np.full(48, 0.25), np.full(23, 0.5)]),
        '100': np.concatenate([np.full(66, 0.25), np.full(34, 0.5)]),
        '200': np.concatenate([np.full(133, 0.25), np.full(67, 0.5)]),
        '300': np.concatenate([np.full(200, 0.25), np.full(100, 0.5)]),
        '1000': np.concatenate([np.full(667, 0.25), np.full(333, 0.5)])
    }
    
    noh_med = {
    '71':[0,  7, 14, 21, 28, 35, 42, 49, 56, 63,  71],
    '100':[0,  10, 20, 30, 40, 50, 60, 70, 80, 90,  100],
    '200': [0,  20, 40, 60, 90, 100, 120, 140, 160, 180,  200],
    '300': [0,  30, 60, 90, 120, 150, 180, 210, 240, 270,  300],
    '1000': [0,  100, 200, 300, 400, 500, 600, 700, 800, 900,  1000]
}
    return caminho_malhas[opcao], sigma_reais[opcao], noh_med[opcao]
