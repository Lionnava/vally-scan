import csv
import numpy as np
import matplotlib.pyplot as plt
from prody import *
from scipy.stats import pearsonr
from fpdf import FPDF
import datetime
import os

# --- 1. NÚCLEO ESTABLE DE INFRAESTRUCTURA ---
def setup_vally_environment():
    for folder in ['Input_PDB', 'Reports', 'Plots', 'Database']:
        if not os.path.exists(folder): os.makedirs(folder)

# --- 2. MOTOR DE RENDERIZADO PREMIUM (Fiel a la captura) ---
class VALLY_Premium_Report(FPDF):
    def __init__(self, pdb_id):
        super().__init__()
        self.pdb_id = pdb_id
        self.logo_path = 'logo_proyecto.png'

    def header(self):
        # BANNER PRINCIPAL (Azul Oxford)
        self.set_fill_color(0, 32, 63)
        self.rect(0, 0, 210, 45, 'F')
        
        # BANNER DE IDENTIDAD (Azul Soberanía)
        self.set_fill_color(0, 50, 90)
        self.rect(0, 45, 210, 8, 'F')
        self.set_xy(0, 46.5)
        self.set_text_color(255, 255, 255)
        self.set_font("Helvetica", 'B', 9)
        self.cell(210, 5, "CIENCIA Y SOBERANIA EN SALUD | PROYECTO VALLY", 0, 0, 'C')

        # LOGOTIPO E INSTRUMENTO
        logo_x = 15
        if os.path.exists(self.logo_path):
            self.image(self.logo_path, 15, 10, 28)
            logo_x = 48
            
        self.set_xy(logo_x, 12)
        self.set_text_color(255, 255, 255)
        self.set_font("Helvetica", 'B', 22)
        self.cell(0, 10, "VALLY-SCAN INSTRUMENT", ln=True)
        self.set_font("Helvetica", '', 10); self.set_x(logo_x)
        self.cell(0, 5, "VIBRATIONAL ANALYSIS & LOCAL LIGAND YIELDING", ln=True)
        self.set_font("Helvetica", 'I', 9); self.set_x(logo_x)
        self.cell(0, 5, "TECHNICAL DATA SHEET | SUMMIT EDITION v3.1", ln=True)

    def footer(self):
        self.set_y(-25)
        self.set_draw_color(0, 32, 63)
        self.line(10, self.get_y(), 200, self.get_y())
        self.set_font("Helvetica", 'I', 8); self.set_text_color(100)
        self.ln(2)
        self.cell(0, 10, "PROYECTO VALLY - PROPIEDAD INTELECTUAL DE LIONELL E. NAVA RAMOS", align='C')
        self.set_y(-15); self.cell(0, 10, f"Pagina {self.page_no()}", align='R')

# --- 3. PROCESO UNIVERSAL (LÓGICA BLINDADA) ---
def vally_universal_engine(pdb_file, active_site_residues=None, mode='universal'):
    try:
        setup_vally_environment()
        target = os.path.join('Input_PDB', pdb_file) if os.path.exists(os.path.join('Input_PDB', pdb_file)) else pdb_file
        
        # --- CÁLCULO FÍSICO APROBADO ---
        structure = parsePDB(target)
        calpha = structure.select('protein and name CA')
        anm = ANM(pdb_file); anm.buildHessian(calpha)
        anm.calcModes(n_modes=30)
        msf = calcSqFlucts(anm)
        b_factors = calpha.getBetas()
        r_val, _ = pearsonr(msf, b_factors)
        top_indices = np.argsort(msf)[-5:][::-1]

        # --- GENERACIÓN DE GRÁFICO CIENTÍFICO ---
        plt.figure(figsize=(10, 5))
        m_z = (msf - np.mean(msf)) / np.std(msf)
        b_z = (b_factors - np.mean(b_factors)) / np.std(b_factors)
        plt.plot(m_z, color='#00203F', label='VALLY Simulation (ANM)', lw=2)
        plt.plot(b_z, color='#32CD32', label='Experimental (B-factors)', ls='--', alpha=0.6)
        plt.title(f"VALLY-SCAN Validation | {pdb_file.upper()} | r = {round(r_val, 3)}")
        plt.xlabel("Residue Index"); plt.ylabel("Standardized Fluctuation")
        plt.legend(loc='best', frameon=True, shadow=True)
        plt.grid(True, alpha=0.25)
        plt.savefig(os.path.join('Plots', f"Plot_{pdb_file[:-4]}.png"), dpi=300)
        plt.close()

        # --- CONSTRUCCIÓN DE REPORTE TÉCNICO ---
        pdf = VALLY_Premium_Report(pdb_file.upper())
        pdf.add_page()
        pdf.set_xy(10, 60)
        
        # Bloque I: Validación
        pdf.set_font("Helvetica", 'B', 14); pdf.set_text_color(0, 32, 63)
        pdf.cell(0, 10, "I. SCIENTIFIC VALIDATION (R2 SYSTEM)", ln=True)
        pdf.set_font("Helvetica", '', 11); pdf.set_text_color(0, 0, 0)
        pdf.cell(55, 8, "Target PDB Structure:", 0); pdf.cell(0, 8, pdb_file.upper(), ln=True)
        pdf.cell(55, 8, "Pearson Correlation (r):", 0); pdf.cell(0, 8, f"{round(r_val, 4)}", ln=True)
        pdf.cell(55, 8, "Confidence Level:", 0)
        conf = "OPTIMAL" if r_val > 0.5 else "CONSISTENT" if r_val > 0.2 else "DYNAMIC DISORDER"
        pdf.cell(0, 8, conf, ln=True)

        # Bloque II: Hotspots
        pdf.ln(5); pdf.set_font("Helvetica", 'B', 14); pdf.set_text_color(0, 32, 63)
        pdf.cell(0, 10, "II. VIBRATIONAL HOTSPOTS MAPPING", ln=True)
        pdf.set_font("Helvetica", '', 10); pdf.set_text_color(0, 0, 0)
        for i, idx in enumerate(top_indices, 1):
            pdf.cell(0, 7, f"Rank {i} -> Residue Index: {idx} | Allosteric Site Propensity", ln=True)

        # Bloque III: Gráfico Integrado (Enriquecimiento)
        pdf.ln(5); pdf.set_font("Helvetica", 'B', 14); pdf.set_text_color(0, 32, 63)
        pdf.cell(0, 10, "III. DYNAMIC VARIANCE PLOT", ln=True)
        pdf.image(os.path.join('Plots', f"Plot_{pdb_file[:-4]}.png"), x=15, y=pdf.get_y()+2, w=180)

        # SALIDA FINAL
        pdf.output(os.path.join('Reports', f"VALLY_Scan_Report_{pdb_file[:-4]}.pdf"))

        # REGISTRO CSV SEGURO
        db_path = os.path.join('Database', 'VALLY_Scan_Master.csv')
        hot_str = "-".join(map(str, top_indices))
        with open(db_path, 'a', newline='') as f:
            writer = csv.writer(f)
            if os.path.getsize(db_path) == 0:
                writer.writerow(['Timestamp', 'PDB', 'Pearson_R', 'Hotspots'])
            writer.writerow([datetime.datetime.now().strftime("%Y-%m-%d %H:%M"), pdb_file, round(r_val, 4), hot_str])

        print(f"--> [VALLY-SCAN SUCCESS] {pdb_file} | r={round(r_val, 3)}")

    except Exception as e:
        print(f"--> [CRITICAL ERROR] {str(e)}")
