import numpy as np
import pandas as pd
from prody import *
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
from scipy.stats import pearsonr

"""
VALLY-Scan v1.7 (Universal & Summit Edition 2026)
Autor: Ing. Lionell E. Nava Ramos
Metodología: Triple Factor de Validación (Biofísica & Ciencia de Materiales)
"""

def vally_universal_engine(pdb_file, active_site_residues=None, mode='viral'):
    print(f"\n=== Iniciando Motor VALLY v1.7 | Modo: {mode.upper()} | Archivo: {pdb_file} ===")
    
    # Carga de estructura
    structure = parsePDB(pdb_file)
    
    # --- CONFIGURACIÓN DE SELECTOR (Universalidad) ---
    # Para virus usamos Carbonos Alfa; para minerales usamos todos los átomos de la red
    selection_query = 'protein and name CA' if mode == 'viral' else 'all'
    nodes = structure.select(selection_query)
    
    if nodes is None:
        print("Error: No se encontraron nodos con el selector actual.")
        return
    
    # --- FACTOR 1: Dinámica Física (ANM) ---
    anm = ANM(pdb_file)
    anm.buildHessian(nodes)
    anm.calcModes(n_modes=20)
    
    msf_predicted = anm.getMSFs()
    node_indices = nodes.getResnums() if mode == 'viral' else np.arange(len(nodes))
    coords = nodes.getCoords()
    
    # --- FACTOR 3: Validación Cristalográfica (B-Factors) ---
    b_factors_experimental = nodes.getBetas()
    correlation, _ = pearsonr(msf_predicted, b_factors_experimental)
    
    # --- FACTOR 2: Heurística de Exclusión (Solo para modo Viral/Bio) ---
    allosteric_results = []
    if mode == 'viral' and active_site_residues:
        active_coords = []
        for res_num in active_site_residues:
            sel = nodes.select(f'resnum {res_num}')
            if sel: active_coords.append(sel.getCoords()[0])
        
        active_coords = np.array(active_coords)
        msf_threshold = np.percentile(msf_predicted, 85)
        
        for i in range(len(nodes)):
            min_dist = np.min(cdist(coords[i].reshape(1, -1), active_coords))
            if msf_predicted[i] >= msf_threshold and min_dist > 15.0:
                allosteric_results.append({
                    'ID': node_indices[i],
                    'Flex': round(msf_predicted[i], 4),
                    'B-Factor': round(b_factors_experimental[i], 2),
                    'Distancia': round(min_dist, 2)
                })

    # --- VISUALIZACIÓN DE TRIPLE FACTOR ---
    plt.figure(figsize=(12, 5))
    plt.plot(node_indices, msf_predicted, label='Teoría (ANM)', color='blue', lw=2)
    plt.plot(node_indices, b_factors_experimental, label='Realidad (Exp)', color='green', alpha=0.4, linestyle='--')
    plt.title(f'VALLY-Scan v1.7: {pdb_file} | Correlación r = {correlation:.2f}')
    plt.xlabel('Índice del Nodo (Átomo/Residuo)')
    plt.ylabel('Amplitud de Vibración')
    plt.legend()
    plt.grid(True, alpha=0.2)
    plt.savefig(f'VALLY_v1_7_{pdb_file}.png')
    
    print(f"Validación r={correlation:.2f} completada.")
    return pd.DataFrame(allosteric_results)

if __name__ == "__main__":
    # EJEMPLO VIRAL (SARS-CoV-2)
    vally_universal_engine('6LU7.pdb', [41, 145], mode='viral')
    
    # EJEMPLO MINERAL (Estructura cristalina)
    # Nota: Requiere un archivo PDB de un mineral, ej: 'quartz.pdb'
    # vally_universal_engine('quartz.pdb', mode='mineral')
