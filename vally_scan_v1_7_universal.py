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

# --- 1. CLASE DE DISEÑO PREMIUM ---
class VALLY_Premium_Report(FPDF):
    def __init__(self):
        super().__init__()
        self.logo_path = 'logo_vally.jpeg' # Actualizado al nombre de su archivo real

    def header(self):
        self.set_fill_color(0, 32, 63)
        self.rect(0, 0, 210, 45, 'F')
        logo_existe = False
        if os.path.exists(self.logo_path):
            try:
                self.image(self.logo_path, 12, 8, 30)
                logo_existe = True
            except: pass
        self.set_xy(48, 12) if logo_existe else self.set_xy(12, 12)
        self.set_text_color(255, 255, 255)
        self.set_font("Helvetica", 'B', 22)
        self.cell(0, 10, "V.A.L.L.Y. PROJECT", ln=True)
        self.set_font("Helvetica", '', 10)
        if logo_existe: self.set_x(48)
        self.cell(0, 5, "Framework v1.7 | Ciencia y Soberania en Salud", ln=True)

    def add_security_watermark(self):
        self.set_font('Helvetica', 'B', 40)
        self.set_text_color(245, 245, 245)
        self.set_xy(0, 140)
        self.cell(210, 20, "VALLY PROJECT - CONFIDENTIAL DATA", 0, 0, 'C')

# --- 2. FUNCIONES DE APOYO (DEFINIDAS ANTES DEL MOTOR) ---
def get_system_specs():
    cpu = platform.processor() or "Intel64 Family 6"
    ram = f"{round(psutil.virtual_memory().total / (1024**3))} GB RAM"
    os_info = f"{platform.system()} {platform.release()}"
    return cpu, ram, os_info

def update_vally_database(data_row):
    db_file = 'VALLY_Master_Database.csv'
    file_exists = os.path.isfile(db_file)
    with open(db_file, mode='a', newline='') as f:
        writer = csv.writer(f)
        if not file_exists:
            writer.writerow(['Timestamp', 'PDB_ID', 'Pearson_R', 'Top_Hotspots', 'CPU', 'RAM'])
        writer.writerow(data_row)

# --- 3. MOTOR PRINCIPAL ---
def vally_universal_engine(pdb_file, active_site_residues=None, mode='universal'):
    try:
        if not os.path.exists('Vally_Plots'): os.makedirs('Vally_Plots')
        
        structure = parsePDB(pdb_file)
        calpha = structure.select('protein and name CA')
        anm = ANM(pdb_file)
        anm.buildHessian(calpha)
        anm.calcModes(n_modes=20)
        
        msf_predicted = calcSqFlucts(anm)
        b_factors = calpha.getBetas()
        r_val, _ = pearsonr(msf_predicted, b_factors)
        top_indices = np.argsort(msf_predicted)[-5:][::-1]
        
        # Guardar Gráfico
        plt.figure(figsize=(10, 5))
        plt.plot(msf_predicted, label='VALLY Simulation', color='#00203F')
        plt.plot(b_factors / np.max(b_factors) * np.max(msf_predicted), 
                 label='Experimental', color='#32CD32', linestyle='--')
        plt.title(f"V.A.L.L.Y. v1.7 | {pdb_file} | r = {round(r_val, 3)}")
        plot_path = f"Vally_Plots/Validacion_{pdb_file.replace('.pdb','')}.png"
        plt.savefig(plot_path)
        
        # Actualizar Base de Datos
        cpu, ram, _ = get_system_specs()
        update_vally_database([datetime.datetime.now(), pdb_file, round(r_val, 4), list(top_indices), cpu, ram])
        
        print(f"\n--> [OK] Base de datos actualizada.")
        print(f"--> [OK] Grafico guardado en: {plot_path}")
        
    except Exception as e:
        print(f"\n[!] ERROR: {str(e)}")
