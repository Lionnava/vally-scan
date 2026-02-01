import numpy as np
import pandas as pd
from prody import *
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
from scipy.stats import pearsonr

"""
VALLY-Scan v1.6 (Summit Edition 2026)
Autor: Ing. Lionell E. Nava Ramos
Metodología: Triple Factor de Validación
1. Dinámica Física (ANM)
2. Heurística Geométrica (Exclusión Alostérica)
3. Validación Cristalográfica (B-Factors)
"""

def vally_triple_factor_analysis(pdb_file, active_site_residues):
    print(f"\n=== Iniciando Análisis de Triple Factor: {pdb_file} ===")
    
    # --- FACTOR 1: Dinámica Física (ANM) ---
    structure = parsePDB(pdb_file)
    calpha = structure.select('protein and name CA')
    anm = ANM(pdb_file)
    anm.buildHessian(calpha)
    anm.calcModes(n_modes=20)
    
    msf_predicted = anm.getMSFs()
    residues = calpha.getResnums()
    coords = calpha.getCoords()
    
    # --- FACTOR 3: Validación Cristalográfica (B-Factors) ---
    # Extraemos la vibración real medida en el experimento de Rayos X
    b_factors_experimental = calpha.getBetas()
    
    # Calcular Correlación de Pearson (r) entre Predicción vs Realidad
    correlation, _ = pearsonr(msf_predicted, b_factors_experimental)
    print(f"[VALIDACIÓN] Correlación con Cristalografía (Pearson r): {correlation:.2f}")

    # --- FACTOR 2: Heurística de Exclusión Alostérica ---
    active_coords = []
    for res_num in active_site_residues:
        sel = calpha.select(f'resnum {res_num}')
        if sel:
            active_coords.append(sel.getCoords()[0])
    active_coords = np.array(active_coords)

    msf_threshold = np.percentile(msf_predicted, 85)
    dist_threshold = 15.0
    
    allosteric_results = []
    for i in range(len(residues)):
        min_dist = np.min(cdist(coords[i].reshape(1, -1), active_coords))
        
        # Filtro de Doble Factor: Mucha flexibilidad + Lejos del sitio activo
        if msf_predicted[i] >= msf_threshold and min_dist > dist_threshold:
            allosteric_results.append({
                'Residuo': residues[i],
                'Nombre': calpha[i].getResname(),
                'MSF_Pred': round(msf_predicted[i], 4),
                'B-Factor_Exp': round(b_factors_experimental[i], 2),
                'Distancia_Activo': round(min_dist, 2)
            })

    # --- SALIDA Y GRAFICACIÓN ---
    df_allosteric = pd.DataFrame(allosteric_results).sort_values(by='MSF_Pred', ascending=False)
    
    # Gráfica de Triple Validación para el Poster de Denver
    fig, ax1 = plt.subplots(figsize=(12, 6))

    ax1.set_xlabel('Número de Residuo')
    ax1.set_ylabel('MSF (Predicción VALLY)', color='blue')
    ax1.plot(residues, msf_predicted, label='Predicción Física (ANM)', color='blue', alpha=0.7)
    ax1.tick_params(axis='y', labelcolor='blue')

    ax2 = ax1.twinx()
    ax2.set_ylabel('B-Factors (Cristalografía Real)', color='green')
    ax2.plot(residues, b_factors_experimental, label='Realidad Experimental', color='green', linestyle='--', alpha=0.5)
    ax2.tick_params(axis='y', labelcolor='green')

    plt.title(f'Triple Validación VALLY-Scan: {pdb_file}\nPearson r = {correlation:.2f}')
    fig.tight_layout()
    plt.savefig(f'Triple_Validacion_{pdb_file}.png')
    
    print("\n--- Top Candidatos Alostéricos Detectados ---")
    print(df_results_summary := df_allosteric.head(5))
    
    return correlation, df_allosteric

if __name__ == "__main__":
    # Prueba con SARS-CoV-2 (6LU7)
    r_6lu7, sitios_6lu7 = vally_triple_factor_analysis('6LU7.pdb', [41, 144, 145, 163])
    
    # Prueba con Dengue (2fom)
    r_2fom, sitios_2fom = vally_triple_factor_analysis('2fom.pdb', [75, 135, 157])
