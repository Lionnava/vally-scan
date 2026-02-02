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

# --- 1. INFRAESTRUCTURA AUTOMÁTICA ---
def setup_vally_environment():
    """Garantiza la existencia de las carpetas Plots, Reports y Database"""
    carpetas = ['Input_PDB', 'Reports', 'Plots', 'Database']
    for folder in carpetas:
        if not os.path.exists(folder):
            os.makedirs(folder)

# --- 2. DISEÑO PREMIUM DE VANGUARDIA (Identidad Proyecto Vally) ---
class VALLY_Premium_Report(FPDF):
    def __init__(self, pdb_id):
        super().__init__()
        self.pdb_id = pdb_id
        self.logo_path = 'logo_proyecto.png'

    def header(self):
        # Banners Sólidos (Captura Aprobada)
        self.set_fill_color(0, 32, 63); self.rect(0, 0, 210, 45, 'F')
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

    def footer(self):
        self.set_y(-25); self.set_draw_color(0, 32, 63); self.line(10, self.get_y(), 200, self.get_y())
        self.set_font("Helvetica", 'I', 8); self.set_text_color(100); self.ln(2)
        self.cell(0, 10, "PROYECTO VALLY - PROPIEDAD INTELECTUAL DE LIONELL E. NAVA RAMOS", align='C')

# --- 3. MOTOR UNIVERSAL (v3.2.1) ---
def vally_universal_engine(pdb_file, active_site_residues=None, mode='universal'):
    try:
        setup_vally_environment()
        target = os.path.join('Input_PDB', pdb_file) if os.path.exists(os.path.join('Input_PDB', pdb_file)) else pdb_file
        
        # Metadatos del Sistema
        info_sys = {
            'os': f"{platform.system()} {platform.release()}",
            'cpu': platform.processor(),
            'ram': f"{round(psutil.virtual_memory().total / (1024**3))} GB",
            'time': datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
        }

        # Núcleo ProDy (Física Blindada)
        struct = parsePDB(target)
        calpha = struct.select('protein and name CA')
        anm = ANM(pdb_file); anm.buildHessian(calpha)
        anm.calcModes(n_modes=30)
        msf = calcSqFlucts(anm)
        b_factors = calpha.getBetas()
        r_val, _ = pearsonr(msf, b_factors)
        top_indices = np.argsort(msf)[-5:][::-1]

        # Gráfico Científico
        plt.figure(figsize=(10, 5))
        m_z = (msf - np.mean(msf)) / np.std(msf); b_z = (b_factors - np.mean(b_factors)) / np.std(b_factors)
        plt.plot(m_z, color='#00203F', label='Simulación VALLY (ANM)', lw=2)
        plt.plot(b_z, color='#32CD32', label='Datos Experimentales (PDB)', ls='--', alpha=0.6)
        plt.title(f"Validación VALLY-Scan | {pdb_file.upper()} | r = {round(r_val, 3)}")
        plt.legend(); plt.grid(True, alpha=0.2); plt.xlabel("Índice de Residuo")
        
        plot_path = os.path.join('Plots', f"Plot_{pdb_file[:-4]}.png")
        plt.savefig(plot_path, dpi=300); plt.close()

        # Reporte PDF (Enriquecido)
        pdf = VALLY_Premium_Report(pdb_file.upper())
        pdf.add_page(); pdf.set_xy(10, 60)
        
        pdf.set_font("Helvetica", 'B', 14); pdf.set_text_color(0, 32, 63)
        pdf.cell(0, 10, "I. SCIENTIFIC VALIDATION & HARDWARE LOG", ln=True)
        pdf.set_font("Helvetica", '', 10); pdf.set_text_color(0, 0, 0)
        pdf.cell(50, 7, "Pearson Correlation:", 0); pdf.cell(0, 7, f"{round(r_val, 4)}", ln=True)
        pdf.cell(50, 7, "CPU Architecture:", 0); pdf.cell(0, 7, info_sys['cpu'], ln=True)
        pdf.cell(50, 7, "Memory RAM:", 0); pdf.cell(0, 7, info_sys['ram'], ln=True)

        pdf.ln(5); pdf.set_font("Helvetica", 'B', 14); pdf.set_text_color(0, 32, 63)
        pdf.cell(0, 10, "II. VIBRATIONAL HOTSPOTS MAPPING", ln=True)
        pdf.set_font("Helvetica", '', 10); pdf.set_text_color(0, 0, 0)
        for i, idx in enumerate(top_indices, 1):
            pdf.cell(0, 7, f"Rank {i} -> Residue: {idx} | Potential Allosteric Candidate", ln=True)

        # Gráfico embebido
        pdf.ln(5); pdf.image(plot_path, x=15, y=pdf.get_y(), w=180)

        pdf.output(os.path.join('Reports', f"VALLY_Scan_Report_{pdb_file[:-4]}.pdf"))

        # Base de Datos (Trazabilidad)
        db_path = os.path.join('Database', 'VALLY_Scan_Master.csv')
        hot_str = "-".join(map(str, top_indices))
        with open(db_path, 'a', newline='') as f:
            writer = csv.writer(f)
            if os.path.getsize(db_path) == 0:
                writer.writerow(['Timestamp', 'PDB', 'Pearson_R', 'Hotspots', 'CPU', 'RAM'])
            writer.writerow([info_sys['time'], pdb_file, round(r_val, 4), hot_str, info_sys['cpu'], info_sys['ram']])

        print(f"--> [VALLY-SCAN OK] {pdb_file} procesado bajo v1.7.4-Universal (Core v3.2.1)")

    except Exception as e:
        print(f"--> [ERROR] {str(e)}")
