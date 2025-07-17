# -*- coding: utf-8 -*-
"""
Created on Wed Jul 16 14:21:36 2025

@author: rodgr
"""

import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore

def plot_EIT(matriz_coordenadas_b,
                   matriz_topologia_b,
                   sigma_inicial_b,
                   sigma_real_b,
                   noh_medidos,
                   alpha_b,
                   lambda_b,
                   std):
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
    
    
    app = pg.mkQApp("Plotting Example")
    #mw = QtWidgets.QMainWindow()
    #mw.resize(800,800)
    
    win = pg.GraphicsLayoutWidget(show=True, title="Basic plotting examples")
    win.resize(1000,600)
    win.setWindowTitle('pyqtgraph example: Plotting')
    
    # Enable antialiasing for prettier plots
    pg.setConfigOptions(antialias=True)
    
    
    
    p = win.addPlot(title="Distribuição de σ nos elementos 1D")

    # Linha sigma calculado (sem marcadores)
    p.plot(centros, valores, pen=pg.mkPen(color=(255,0,0), width=2), name='σ calculado')

    # Linha sigma real (linestyle pontilhado)
    p.plot(centros, sigma_real_b, pen=pg.mkPen(color=(0,0,255), style=pg.QtCore.Qt.DotLine, width=2), name='σ real')

    # Pontos medidos (marcador em X)
    p.plot(pos_medidas, pos_valores_real, pen=None, symbol='x', symbolSize=15, symbolBrush='r', name='Ptos medidos')

    # Limites do gráfico
    p.setXRange(0, 1.01)
    p.setYRange(0.0, 0.60)

    # Labels dos eixos
    
    p.setLabel('left', 'Condutividade σ [S/m]', units='')
    p.setLabel('bottom', 'Posição [m]', units='')
    # Título (no topo da janela)
    win.setWindowTitle(f'Distribuição de σ nos elementos 1D')

    # Texto complementar (substituto do plt.text)
    texto = pg.TextItem(html=f'<div style="text-align: center;">α = {alpha_b}, λ = {lambda_b}, std = {std}</div>', anchor=(0.5,0.5))
    texto.setPos(0.5, 0.55)  # define posição do texto (X, Y)
    p.addItem(texto)

    # Legenda
    legend = pg.LegendItem(offset=(70, 20))
    legend.setParentItem(p.graphicsItem())
    legend.addItem(p.listDataItems()[0], 'σ calculado')
    legend.addItem(p.listDataItems()[1], 'σ real')
    legend.addItem(p.listDataItems()[2], 'Ptos medidos')
    pg.exec()
    #if __name__ == '__main__':
    #    pg.exec()
    
