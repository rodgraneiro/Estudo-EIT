# Estudo-EIT

Este código resolve o 'Problema Direto' e  'Problema Inverso' para determinar a condutividade  $\sigma$ de uma barra de 1m, discretizada numa lalha de Elementos Finitos 1D, por meio da minimização pelo método 'Gauss-Newton' da equação: 



$$
\hat{\sigma}_{k+1} = \hat{\sigma}_k + \alpha_k
\left( J_k^T W_1 J_k + \lambda^2 L_2^T L_2 \right)^{-1}
\cdot \left( J_k^T W_1 (z - h(\hat{\sigma}_k)) - \lambda^2 L_2^T L_2 (\hat{\sigma}_k - \sigma^*) \right)
$$


onde,
- $\sigma$ é condutividade $S/m$;
- $\alpha$ é uma constante de regularização;
- $\lambda$ é uma constante de regularização;
- $W_1$ é a matriz Identidade;
- $J$ é a matriz jacobiana;
- $L_2$ é a matriz de um Filtro Passa-Alta.

Dimenssões da barra:
- comprimento = 1 $sigma$
- Área $1\times 10^{-3}\: m^2$.

A barra está submetida a uma corrente de $1\times 10^{-3} A$ nas extremidades e pontos aleatórios de tensão são medidos/calculado ao longo da barra.

Para rodar esse código em python são necessários as seguintes bibliotecas;
- numpy;
- matplotlib;
- meshio;
- time;
- PyQtgraph;
- QtCore;
- PyQt5;
- PySide6.

Para rodar esse código também são necessários os sequites
arquivos:
- EIT_1D.py programa principal;
- EIT_functions.py contém as funções necessárias;
- EIT_mesh_1D.py com as seguintes malhas 1D:
  - unidimensional_71e_py.msh 71 elementos;
  - unidimensional_100e_py.msh 100 elementos;
  - unidimensional_200e_py.msh 200 elementos;
  - unidimensional_300e_py.msh 300 elementos;
  - unidimensional_1000e_py.msh 1000 elementos.
A malha 'unidimensional_71e_py' contém 71 elementos com tamanhos diferentes de 0.01, 0.02 e 0.03 metros. Nos demais arquivos os elementos têm tamanhos iguais.

---

Para escolher uma malha, descomente a linha de comando:

` ` `python
#opcao = input('Escolha a malha(71, 100, 200, 300, 1000): ')
` ` ` 

ou altere a seguinte linha de comando como desejado.

` ` `
opcao = '71'
` ` ` 

---

Para escolher os nós de medição descomente a linha de comando:

` ` `
#entrada = input("Digite os nós medidos separados por vírgula (ex.: 0,10,20: ")
` ` ` 

ou altere a seguinte linha de comando como desejado.

` ` `
noh_medidos = [0,  7, 14, 21, 28, 35, 42, 49, 56, 63,  n_elements] 
` ` ` 

---




