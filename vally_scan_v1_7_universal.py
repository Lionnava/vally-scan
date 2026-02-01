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

# --- 1. INFRAESTRUCTURA ---
def setup_vally_environment():
    for folder in ['Input_PDB', 'Reports', 'Plots', 'Database']:
        if not os.path.exists(folder): os.makedirs(folder)

# --- 2. DISEÑO DE REPORTE (ALTO NIVEL) ---
class VALLY_Report_Master(FPDF):
    def __init__(self, pdb_id):
        super().__init__()
        self.pdb_id = pdb_id
        self.logo_path = 'logo_proyecto.png'

    def header(self):
        # Banner Principal Azul Oxford
        self.set_fill_color(0, 32, 63); self.rect(0, 0, 210, 40, 'F')
        # Banner Secundario
        self.set_fill_color(0, 50, 90); self.rect(0, 40, 210, 7, 'F')
        self.set_xy(0, 41); self.set_text_color(255, 255, 255); self.set_font("Helvetica", 'B', 8)
        self.cell(210, 5, "CIENCIA Y SOBERANIA EN SALUD | PROYECTO VALLY", 0, 0, 'C')

        if os.path.exists(self.logo_path):
            self.image(self.logo_path, 12, 8, 25)
        
        self.set_xy(40, 10); self.set_font("Helvetica", 'B', 22)
        self.cell(0, 10, "VALLY-SCAN INSTRUMENT", ln=True)
        self.set_font("Helvetica", '', 10); self.set_x(40)
        self.cell(0, 5, "VIBRATIONAL ANALYSIS & LOCAL LIGAND YIELDING", ln=True)

    def footer(self):
        self.set_y(-15); self.set_font("Helvetica", 'I', 8); self.set_text_color(150)
        self.cell(0, 10, "PROYECTO VALLY - PROPIEDAD INTELECTUAL DE LIONELL E. NAVA RAMOS", align='C')

# --- 3. MOTOR DE ANÁLISIS ---
def vally_universal_engine(pdb_file, active_site_residues=None, mode='universal'):
    try:
        setup_vally_environment()
        target_path = os.path.join('Input_PDB', pdb_file) if os.path.exists(os.path.join('Input_PDB', pdb_file)) else pdb_file
        
        # Física Computacional
        structure = parsePDB(target_path)
        calpha = structure.select('protein and name CA')
        anm = ANM(pdb_file); anm.buildHessian(calpha)
        anm.calcModes(n_modes=25)
        
        msf = calcSqFlucts(anm)
        b_factors = calpha.getBetas()
        r_val, _ = pearsonr(msf, b_factors)
        hotspots = np.argsort(msf)[-5:][::-1]

        # --- GRÁFICO CIENTÍFICO CON LEYENDA ---
        plt.figure(figsize=(10, 5))
        # Normalización Z-Score para precisión de vanguardia
        msf_z = (msf - np.mean(msf)) / np.std(msf)
        b_z = (b_factors - np.mean(b_factors)) / np.std(b_factors)
        
        plt.plot(msf_z, color='#00203F', label='Simulación VALLY (Física)', lw=2)
        plt.plot(b_z, color='#32CD32', label='Experimental (B-Factors)', linestyle='--', alpha=0.6)
        plt.title(f"Validación VALLY-Scan: {pdb_file.upper()} | r = {round(r_val, 3)}")
        plt.xlabel("Índice de Residuo"); plt.ylabel("Amplitud Normalizada")
        plt.legend(loc='upper right', frameon=True, shadow=True)
        plt.grid(True, alpha=0.3)
        
        plot_path = os.path.join('Plots', f"Plot_{pdb_file[:-4]}.png")
        plt.savefig(plot_path, dpi=300); plt.close()

        # --- REPORTE PDF ---
        pdf = VALLY_Report_Master(pdb_file.upper())
        pdf.add_page(); pdf.set_xy(10, 55)
        
        pdf.set_font("Helvetica", 'B', 14); pdf.set_text_color(0, 32, 63)
        pdf.cell(0, 10, "METRICAS DE ANALISIS MOLECULAR", ln=True)
        pdf.set_font("Helvetica", '', 11); pdf.set_text_color(0, 0, 0)
        pdf.cell(50, 8, "Pearson r:", 0); pdf.cell(0, 8, f"{round(r_val, 4)}", ln=True)
        
        pdf.ln(5); pdf.set_font("Helvetica", 'B', 14); pdf.cell(0, 10, "HOTSPOTS VIBRACIONALES", ln=True)
        pdf.set_font("Helvetica", '', 10)
        for i, h in enumerate(hotspots, 1):
            pdf.cell(0, 7, f"Rank {i} -> Residuo {h} (Zona de Plasticidad)", ln=True)

        pdf.ln(10); pdf.set_font("Helvetica", 'B', 9); pdf.set_text_color(180, 0, 0)
        pdf.multi_cell(0, 5, "NOTA: El sistema utiliza normalización vectorial para garantizar la paridad entre el modelo ANM y los datos experimentales.")

        pdf_path = os.path.join('Reports', f"Reporte_VALLY_{pdb_file[:-4]}.pdf")
        pdf.output(pdf_path)

        # --- DATABASE CSV (CORRECCIÓN FINAL) ---
        db_path = os.path.join('Database', 'VALLY_Scan_Master.csv')
        # Guardamos la lista de hotspots como string con guiones para NO romper el CSV
        hotspots_str = "-".join(map(str, hotspots))
        with open(db_path, 'a', newline='') as f:
            writer = csv.writer(f)
            if os.path.getsize(db_path) == 0:
                writer.writerow(['Fecha', 'Archivo', 'Pearson_R', 'Hotspots_Ranked'])
            writer.writerow([datetime.datetime.now().strftime("%Y-%m-%d %H:%M"), pdb_file, round(r_val, 4), hotspots_str])

        print(f"\n--> [VALLY-SCAN] PROCESO EXITOSO")
        print(f"--> Coeficiente r: {round(r_val, 3)}")
        print(f"--> Archivos en /Reports y /Plots")

    except Exception as e:
        print(f"\n[!] ERROR EN EL INSTRUMENTO: {str(e)}")
