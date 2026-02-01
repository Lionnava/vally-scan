import numpy as np
import pandas as pd
from prody import *
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
import os

"""
VALLY-Scan v1.6: Framework Computacional de Doble Factor (R2)
Autor: Ing. Lionell E. Nava Ramos
Descripción: Implementación de Heurística Informada para la detección 
automática de sitios alostéricos combinando Dinámica ANM y Geometría.
"""

def vally_framework_execution(pdb_id, active_site_list):
    print(f"\n--- Iniciando VALLY-Scan v1.6 para: {pdb_id} ---")
    
    # 1. Carga de estructura y Modelo de Red Anisotrópica (ANM)
    structure = parsePDB(pdb_id)
    calpha = structure.select('protein and name CA')
    anm = ANM(pdb_id)
    anm.buildHessian(calpha)
    anm.calcModes(n_modes=20)
    
    # 2. Cálculo de Fluctuaciones (Factor 1: Dinámica Física)
    msf = list(abs(anm.getMSFs()))
    residues = calpha.getResnums()
    coords = calpha.getCoords()
    
    # 3. Módulo de Heurística de Doble Factor (R2) - NUEVO
    # Extraemos coordenadas del sitio activo para calcular distancias
    active_coords = []
    for res_num in active_site_list:
        sel = calpha.select(f'resnum {res_num}')
        if sel:
            active_coords.append(sel.getCoords()[0])
    
    active_coords = np.array(active_coords)
    
    # Umbrales para la IA/Heurística
    msf_threshold = np.percentile(msf, 85) # Top 15% de flexibilidad
    dist_threshold = 15.0 # Mínimo 15 Amstrongs del sitio activo
    
    allosteric_candidates = []
    
    print("Analizando picos de flexibilidad y exclusión geométrica...")
    
    for i in range(len(residues)):
        res_flex = msf[i]
        res_coord = coords[i].reshape(1, -1)
        
        # Calcular distancia mínima al sitio catalítico
        min_dist = np.min(cdist(res_coord, active_coords))
        
        # Lógica de Doble Factor
        if res_flex >= msf_threshold and min_dist > dist_threshold:
            allosteric_candidates.append({
                'Residuo': residues[i],
                'Nombre': calpha[i].getResname(),
                'Flexibilidad': round(res_flex, 4),
                'Distancia_Activo': round(min_dist, 2)
            })

    # 4. Generación de Reporte y Validación
    df_results = pd.DataFrame(allosteric_candidates)
    df_results = df_results.sort_values(by='Flexibilidad', ascending=False)
    
    print("\n[RESULTADO] Sitios Alostéricos Potenciales Detectados:")
    print(df_results.head(5))
    
    # Graficar resultados para validación visual (Poster Denver 2026)
    plt.figure(figsize=(10, 6))
    plt.plot(residues, msf, label='Flexibilidad Intrínseca (ANM)', color='blue')
    plt.axhline(y=msf_threshold, color='red', linestyle='--', label='Umbral Alostérico')
    plt.title(f'VALLY-Scan v1.6: Análisis de Doble Factor - {pdb_id}')
    plt.xlabel('Número de Residuo')
    plt.ylabel('MSF (Fluctuación)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(f'Validacion_VALLY_{pdb_id}.png')
    
    return df_results

# --- EJECUCIÓN DE PRUEBA (Basado en tus datos validados) ---
if __name__ == "__main__":
    # Caso 1: SARS-CoV-2 (6LU7) - Sitio activo aprox. 41, 145, 163
    # Tus reportes validados muestran picos en 277, 278, 279
    active_6lu7 = [41, 144, 145, 163, 164]
    resultados_covid = vally_framework_execution('6LU7.pdb', active_6lu7)
    
    # Caso 2: Dengue (2FOM) - Sitio activo aprox. 135, 157, 75
    # Tus reportes muestran picos en 43, 44, 45
    active_2fom = [75, 135, 157]
    resultados_dengue = vally_framework_execution('2fom.pdb', active_2fom)
