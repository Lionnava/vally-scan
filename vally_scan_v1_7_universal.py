import numpy as np
import matplotlib.pyplot as plt
from prody import parsePDB, ANM, calcSqFlucts

def vally_universal_engine(pdb_file, active_site_residues=None, mode='viral'):
    print(f"\n=== Iniciando Motor VALLY v1.7 | Modo: {mode.upper()} | Archivo: {pdb_file} ===")
    
    # 1. Carga de estructura
    structure = parsePDB(pdb_file)
    calpha = structure.select('protein and name CA') if mode == 'viral' else structure.select('all')
    
    # 2. Cálculo de Dinámica Física (Factor 1)
    anm = ANM('Analysis')
    anm.buildHessian(calpha)
    anm.calcModes(20)
    
    # 3. Validación de Triple Factor (Factor 3 - Experimental)
    # calcSqFlucts es la función exacta para obtener MSF de un objeto ANM
    msf_predicted = calcSqFlucts(anm)
    b_factors = calpha.getBetas()
    
    # 4. Cálculo de Correlación Pearson (r)
    # Ajustamos longitudes por si acaso
    correlation = np.corrcoef(msf_predicted, b_factors)[0, 1]
    
    # 5. Visualización
    plt.figure(figsize=(11, 5))
    
    # Normalizamos para comparar visualmente mejor
    msf_norm = msf_predicted / msf_predicted.max()
    bfac_norm = b_factors / b_factors.max()
    
    plt.plot(msf_norm, label=f'Predicción VALLY (Simulación)', color='#1f77b4', linewidth=2)
    plt.plot(bfac_norm, label='Validación Experimental (B-factors)', color='#2ca02c', alpha=0.6, linestyle='--')
    
    plt.title(f'V.A.L.L.Y. v1.7 - Validación de Triple Factor (R2 System)\nEstructura: {pdb_file} | Pearson r = {correlation:.3f}', fontsize=12)
    plt.xlabel('Índice del Residuo / Átomo', fontsize=10)
    plt.ylabel('Fluctuaciones Normalizadas', fontsize=10)
    plt.legend()
    plt.grid(True, linestyle=':', alpha=0.6)
    
    print(f"--> [OK] Validación finalizada.")
    print(f"--> Coeficiente Pearson r: {correlation:.3f}")
    
    # Guardar automáticamente para tu reporte
    plt.savefig('resultado_validacion_vally.png', dpi=300)
    print("--> [INFO] Gráfico guardado como 'resultado_validacion_vally.png'")
    
    plt.show()

    return correlation
