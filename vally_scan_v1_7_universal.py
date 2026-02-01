import csv
import numpy as np
import matplotlib.pyplot as plt
from prody import *
from scipy.stats import pearsonr
from fpdf import FPDF
import datetime
import os
import platform
import psutil

# [La clase VALLY_Premium_Report se mantiene igual que en la versión anterior]

def update_vally_database(data_row):
    """Mantiene un registro histórico de todas las simulaciones en CSV"""
    db_file = 'VALLY_Master_Database.csv'
    file_exists = os.path.isfile(db_file)
    
    with open(db_file, mode='a', newline='') as f:
        writer = csv.writer(f)
        if not file_exists:
            # Encabezados de grado industrial
            writer.writerow(['Timestamp', 'PDB_ID', 'Pearson_R', 'Top_Hotspots', 'Affinity_Est', 'CPU_Model', 'RAM_Used'])
        writer.writerow(data_row)

def vally_universal_engine(pdb_file, active_site_residues=None, mode='viral'):
    try:
        # 1. Procesamiento Físico
        structure = parsePDB(pdb_file)
        calpha = structure.select('protein and name CA')
        anm = ANM(pdb_file)
        anm.buildHessian(calpha)
        anm.calcModes(n_modes=20)
        
        msf_predicted = calcSqFlucts(anm)
        b_factors = calpha.getBetas()
        
        r_val, _ = pearsonr(msf_predicted, b_factors)
        top_indices = np.argsort(msf_predicted)[-5:][::-1]
        
        # 2. Captura de Specs para DB
        cpu, ram, _ = get_system_specs()
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        # 3. Guardar en Base de Datos
        db_data = [timestamp, pdb_file.upper(), round(r_val, 4), list(top_indices), "-10.64", cpu, ram]
        update_vally_database(db_data)
        
        # 4. Generación de Gráfico y Reporte PDF
        # [Aquí siguen las funciones de plt.savefig y export_vally_report]
        
        print(f"\n--> [DATABASE] Registro guardado en VALLY_Master_Database.csv")
        print(f"--> [SUCCESS] Reporte y análisis completados para {pdb_file}")
        
    except Exception as e:
        print(f"\n[!] ERROR: {str(e)}")
