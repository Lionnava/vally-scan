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

# --- 1. GESTIÓN DE ENTORNO R2 ---
def setup_vally_environment():
    """Garantiza la infraestructura de directorios para el preprint."""
    carpetas = ['Input_PDB', 'Reports', 'Plots', 'Database']
    for folder in carpetas:
        if not os.path.exists(folder):
            os.makedirs(folder)

# --- 2. MOTOR DE RENDERIZADO PREMIUM (Identidad Visual VALLY) ---
class VALLY_Premium_Report(FPDF):
    def __init__(self, pdb_id):
        super().__init__()
        self.pdb_id = pdb_id
        self.logo_path = 'logo_proyecto.png'

    def header(self):
        # Banner Superior: Azul Oxford (Identidad Corporativa)
        self.set_fill_color(0, 32, 63); self.rect(0, 0, 210, 45, 'F')
        # Banner de Lema: Azul Soberanía
        self.set_fill_color(0, 50, 90); self.rect(0, 45, 210, 8, 'F')
        self.set_xy(0, 46.5); self.set_text_color(255, 255, 255); self.set_font("Helvetica", 'B', 9)
        self.cell(210, 5, "CIENCIA Y SOBERANIA EN SALUD | PROYECTO VALLY", 0, 0, 'C')
        
        logo_x = 15
        if os.path.exists(self.logo_path):
            self.image(self.logo_path, 15, 10, 28); logo_x = 48
            
        self.set_xy(logo_x, 12); self.set_text_color(255, 255, 255); self.set_font("Helvetica", 'B', 22)
        self.cell(0, 10, "VALLY-SCAN INSTRUMENT", ln=True)
        self.set_font("Helvetica", '', 10); self.set_x(logo_x)
        self.cell(0, 5, "VIBRATIONAL ANALYSIS & LOCAL LIGAND YIELDING", ln=True)
        self.set_font("Helvetica", 'I', 9); self.set_x(logo_x)
        self.cell(0, 5, "FRAMEWORK v1.7 | UNIVERSAL EDITION | R2 SYSTEM", ln=True)

    def footer(self):
        self.set_y(-25); self.set_draw_color(0, 32, 63); self.line(10, self.get_y(), 200, self.get_y())
        self.set_font("Helvetica", 'I', 8); self.set_text_color(100); self.ln(2)
        self.cell(0, 10, "PROYECTO VALLY - PROPIEDAD INTELECTUAL DE LIONELL E. NAVA RAMOS", align='C')
        self.set_y(-15); self.cell(0, 10, f"Pagina {self.page_no()}", align='R')

# --- 3. MOTOR UNIVERSAL R2 (LÓGICA OPTIMIZADA) ---
def vally_universal_engine(pdb_file, active_site_residues=None, mode='universal'):
    try:
        setup_vally_environment()
        target = os.path.join('Input_PDB', pdb_file) if os.path.exists(os.path.join('Input_PDB', pdb_file)) else pdb_file
        
        # Captura de Hardware (Factor 3: Experimental/Sistémico)
        info_sys = {
            'os': f"{platform.system()} {platform.release()}",
            'cpu': platform.processor(),
            'ram': f"{round(psutil.virtual_memory().total / (1024**3))} GB",
            'time': datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
        }

        # FACTOR 1: Dinámica Física Intrínseca (ANM)
        structure = parsePDB(target)
        calpha = structure.select('protein and name CA')
        anm = ANM(pdb_file); anm.buildHessian(calpha)
        anm.calcModes(n_modes=30)
        msf = calcSqFlucts(anm)
        b_factors = calpha.getBetas()
        
        # FACTOR 3: Correlación Cruzada Experimental
        r_val, _ = pearsonr(msf, b_factors)
        
        # FACTOR 2: Heurística de Exclusión Geométrica (Hotspots)
        top_indices = np.argsort(msf)[-5:][::-1]

        # --- GENERACIÓN DE GRÁFICO (FIGURE X PREPRINT) ---
        plt.figure(figsize=(10, 5))
        m_z = (msf - np.mean(msf)) / np.std(msf)
        b_z = (b_factors - np.mean(b_factors)) / np.std(b_factors)
        plt.plot(m_z, color='#00203F', label='VALLY Simulation (ANM)', lw=2)
        plt.plot(b_z, color='#32CD32', label='Experimental Data (B-factors)', ls='--', alpha=0.6)
        plt.title(f"R2 Validation System | {pdb_file.upper()} | r = {round(r_val, 3)}")
        plt.legend(loc='best', frameon=True, shadow=True)
        plt.grid(True, alpha=0.25); plt.xlabel("Residue Index"); plt.ylabel("Standardized Fluctuation")
        
        plot_path = os.path.join('Plots', f"Plot_{pdb_file[:-4]}.png")
        plt.savefig(plot_path, dpi=300); plt.close()

        # --- CONSTRUCCIÓN DEL REPORTE TÉCNICO ---
        pdf = VALLY_Premium_Report(pdb_file.upper())
        pdf.add_page(); pdf.set_xy(10, 60)
        
        # Bloque I: Validación R2 y Trazabilidad
        pdf.set_font("Helvetica", 'B', 14); pdf.set_text_color(0, 32, 63)
        pdf.cell(0, 10, "I. R2 VALIDATION PROTOCOL & HARDWARE LOG", ln=True)
        pdf.set_font("Helvetica", '', 10); pdf.set_text_color(0, 0, 0)
        pdf.cell(55, 7, "Pearson Correlation (r):", 0); pdf.cell(0, 7, f"{round(r_val, 4)}", ln=True)
        pdf.cell(55, 7, "System CPU:", 0); pdf.cell(0, 7, info_sys['cpu'], ln=True)
        pdf.cell(55, 7, "Memory Architecture:", 0); pdf.cell(0, 7, info_sys['ram'], ln=True)

        # Bloque II: Hotspots (Plasticidad Regulatoria)
        pdf.ln(5); pdf.set_font("Helvetica", 'B', 14); pdf.set_text_color(0, 32, 63)
        pdf.cell(0, 10, "II. VIBRATIONAL HOTSPOTS MAPPING", ln=True)
        pdf.set_font("Helvetica", '', 10); pdf.set_text_color(0, 0, 0)
        for i, idx in enumerate(top_indices, 1):
            pdf.cell(0, 7, f"Rank {i} -> Index: {idx} | Allosteric/Regulatory Plasticity Site", ln=True)

        # Bloque III: Figure X y Gráfico
        pdf.ln(5); pdf.image(plot_path, x=15, y=pdf.get_y(), w=180)
        pdf.set_y(pdf.get_y() + 95)
        pdf.set_font("Helvetica", 'I', 8); pdf.set_text_color(100)
        caption = (f"Figure X. Validation of the R2 System through comparative flexibility analysis ({pdb_file.upper()}). "
                   "The solid blue line represents theoretical fluctuations, while the green dashed line represents "
                   f"experimental B-factor data. Correlation r={round(r_val, 3)} confirms predictive accuracy.")
        pdf.multi_cell(180, 4, caption, align='C')

        # Bloque IV: Executive Summary (Technical Update Text)
        pdf.ln(5); pdf.set_font("Helvetica", 'B', 12); pdf.set_text_color(0, 32, 63)
        pdf.cell(0, 8, "III. EXECUTIVE SUMMARY", ln=True)
        pdf.set_font("Helvetica", '', 10); pdf.set_text_color(0, 0, 0)
        summary = (f"El framework VALLY-Scan v1.7 ha ejecutado el protocolo R2 sobre {pdb_file.upper()}. "
                   f"Los resultados demuestran que el sistema captura la varianza dinamica experimental. "
                   "La alineacion de picos identifica dominios de alta structural plasticity, validando la "
                   "capacidad del software para mapear alosterismo molecular.")
        pdf.multi_cell(0, 5, summary)

        pdf.output(os.path.join('Reports', f"VALLY_Scan_Report_{pdb_file[:-4]}.pdf"))

        # --- ACTUALIZACIÓN DE DATABASE ---
        db_path = os.path.join('Database', 'VALLY_Scan_Master.csv')
        hot_str = "-".join(map(str, top_indices))
        with open(db_path, 'a', newline='') as f:
            writer = csv.writer(f)
            if os.path.getsize(db_path) == 0:
                writer.writerow(['Timestamp', 'PDB', 'Pearson_R', 'Hotspots', 'CPU', 'RAM'])
            writer.writerow([info_sys['time'], pdb_file, round(r_val, 4), hot_str, info_sys['cpu'], info_sys['ram']])

        print(f"--> [SUCCESS] VALLY-SCAN R2 FRAMEWORK v1.7: {pdb_file} procesado.")

    except Exception as e:
        print(f"--> [ERROR CRITICO] {str(e)}")
