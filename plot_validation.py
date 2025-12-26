import matplotlib.pyplot as plt
import numpy as np
from prody import *
from scipy.stats import pearsonr

# 1. Cargar datos reales de 6LU7
protein = parsePDB('data/6lu7.pdb')
calphas = protein.select('protein and name CA')
b_factors_exp = calphas.getBetas() # Datos de Rayos X

# 2. Tu cálculo ANM (lo mismo que hace tu main.py)
anm = ANM('6LU7 Validation')
anm.buildHessian(calphas)
anm.calcModes(n_modes=20)
sq_flucts_calc = calcSqFlucts(anm) # Tu predicción

# 3. Normalizar para que ambas curvas estén en la misma escala (0 a 1)
def normalize(data):
    return (data - np.min(data)) / (np.max(data) - np.min(data))

calc_norm = normalize(sq_flucts_calc)
exp_norm = normalize(b_factors_exp)
r_val, _ = pearsonr(calc_norm, exp_norm)

# 4. Crear la gráfica profesional
plt.figure(figsize=(12, 5))
residuos = calphas.getResnums()

plt.plot(residuos, calc_norm, label=f'VALLY-Scan (Calculated)', color='firebrick', linewidth=1.5)
plt.plot(residuos, exp_norm, label=f'Experimental B-factors (Rayos-X)', color='black', linestyle='--', alpha=0.4, linewidth=1)

plt.title(f'Figure 2: Dynamic Validation for SARS-CoV-2 Mpro (r = {r_val:.2f})', fontsize=14)
plt.xlabel('Residue Number', fontsize=12)
plt.ylabel('Normalized Fluctuations', fontsize=12)
plt.legend(loc='upper right')
plt.grid(True, linestyle=':', alpha=0.6)

# Guardar en alta resolución para la revista
plt.savefig('Figure_2_Validation_6LU7.png', dpi=300, bbox_inches='tight')
plt.show()