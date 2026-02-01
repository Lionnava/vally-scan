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

# --- 1. DICCIONARIO BILINGÜE DE ALTA PRECISIÓN ---
LEXICON = {
    'es': {
        'title': "PROYECTO VALLY",
        'subtitle': "VIBRATIONAL ANALYSIS & LOCAL LIGAND YIELDING",
        'instrument': "INSTRUMENTO VALLY-SCAN | V2.0 PRO",
        'sec1': "I. VALIDACIÓN CIENTÍFICA (SISTEMA R2)",
        'sec2': "II. ARQUITECTURA COMPUTACIONAL",
        'sec3': "III. MAPEO DE HOTSPOTS VIBRACIONALES",
        'sec4': "IV. RESUMEN EJECUTIVO",
        'pearson': "Coeficiente Pearson (r):",
        'status': "Nivel de Confianza:",
        'note': "NOTA DEL PROCESO: Datos normalizados mediante Z-Score. Las discrepancias pueden indicar desorden intrínseco.",
        'plot_label_sim': "Simulación VALLY (ANM)",
        'plot_label_exp': "Datos Experimentales (PDB)",
        'plot_title': "Validación VALLY-Scan",
        'rank': "Rango"
    },
    'en': {
        'title': "VALLY PROJECT",
        'subtitle': "VIBRATIONAL ANALYSIS & LOCAL LIGAND YIELDING",
        'instrument': "VALLY-SCAN INSTRUMENT | V2.0 PRO",
        'sec1': "I. SCIENTIFIC VALIDATION (R2 SYSTEM)",
        'sec2': "II. COMPUTATIONAL ARCHITECTURE",
        'sec3': "III. VIBRATIONAL HOTSPOTS MAPPING",
        'sec4': "IV. EXECUTIVE SUMMARY",
        'pearson': "Pearson Coefficient (r):",
        'status': "Confidence Level:",
        'note': "PROCESS NOTE: Data normalized via Z-Score. Discrepancies may indicate intrinsic disorder.",
        'plot_label_sim': "VALLY Simulation (ANM)",
        'plot_label_exp': "Experimental Data (PDB)",
        'plot_title': "VALLY-Scan Validation",
        'rank': "Rank"
    }
}

# --- 2. MOTOR DE REPORTES REFORMADO ---
class VALLY_Report(FPDF):
    def __init__(self, lang='es'):
        super().__init__()
        self.L = LEXICON[lang]
        self.logo_path = 'logo_proyecto.png'

    def header(self):
        self.set_fill_color(0, 32, 63)
        self.rect(0, 0, 210, 40, 'F')
        
        logo_x = 12
        if os.path.exists(self.logo_path):
            try:
                self.image(self.logo_path, 12, 8, 25)
                logo_x = 42
            except: pass
            
        self.set_xy(logo_x, 10)
        self.set_text_color(255, 255, 255)
        self.set_font("Helvetica", 'B', 20)
        self.cell(0, 10, self.L['title'], ln=True)
        self.set_font("Helvetica", '', 9)
        self.set_x(logo_x)
        self.cell(0, 5, self.L['subtitle'], ln=True)
        self.set_x(logo_x)
        self.cell(0, 5, self.L['instrument'], ln=True)

    def footer(self):
        self.set_y(-20)
        self.set_draw_color(0, 32, 63)
        self.line(10, self.get_y(), 200, self.get_y())
        self.set_font("Helvetica", 'I', 8)
        self.set_text_color(120, 120, 120)
        self.cell(0, 10, f"{self.L['title']} - Proprietary Technology", align='C')
        self.cell(0, 10, f"{self.page_no()}", align='R')

# --- 3. MOTOR ANALÍTICO INTEGRADO ---
def run_vally_scan(pdb_id, lang='es'):
    L = LEXICON[lang]
    folders = ['Reports', 'Plots', 'Database', 'Input_PDB']
    for f in folders: 
        if not os.path.exists(f): os.makedirs(f)

    try:
        # Procesamiento
        path = os.path.join('Input_PDB', pdb_id) if os.path.exists(os.path.join('Input_PDB', pdb_id)) else pdb_id
        struct = parsePDB(path)
        calpha = struct.select('protein and name CA')
        anm = ANM(pdb_id)
        anm.buildHessian(calpha)
        anm.calcModes(n_modes=20)
        
        msf = calcSqFlucts(anm)
        b_factors = calpha.getBetas()
        r_val, _ = pearsonr(msf, b_factors)
        hotspots = np.argsort(msf)[-5:][::-1]

        # Gráfico Bilingüe
        plt.figure(figsize=(10, 5))
        m_n = (msf - np.mean(msf)) / np.std(msf)
        b_n = (b_factors - np.mean(b_factors)) / np.std(b_factors)
        
        plt.plot(m_n, color='#00203F', label=L['plot_label_sim'], lw=2)
        plt.plot(b_n, color='#32CD32', label=L['plot_label_exp'], ls='--', alpha=0.6)
        plt.title(f"{L['plot_title']} | {pdb_id.upper()} (r={round(r_val,3)})")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plot_fn = f"Plot_{lang.upper()}_{pdb_id.replace('.pdb','')}.png"
        plt.savefig(os.path.join('Plots', plot_fn), dpi=300)
        plt.close()

        # Generar PDF
        pdf = VALLY_Report(lang)
        pdf.add_page()
        pdf.set_xy(10, 50)
        
        # Secciones
        pdf.set_font("Helvetica", 'B', 12); pdf.cell(0, 10, L['sec1'], ln=True)
        pdf.set_font("Helvetica", '', 10)
        pdf.cell(50, 8, L['pearson']); pdf.cell(0, 8, f"{round(r_val, 4)}", ln=True)
        
        pdf.ln(5); pdf.set_font("Helvetica", 'B', 12); pdf.cell(0, 10, L['sec3'], ln=True)
        pdf.set_font("Helvetica", '', 10)
        for i, h in enumerate(hotspots, 1):
            pdf.cell(0, 7, f"{L['rank']} {i}: Residue {h}", ln=True)

        pdf.ln(5); pdf.set_font("Helvetica", 'B', 9); pdf.set_text_color(170, 0, 0)
        pdf.multi_cell(0, 5, L['note'])

        pdf_fn = f"VALLY_Scan_{lang.upper()}_{pdb_id.replace('.pdb','')}.pdf"
        pdf.output(os.path.join('Reports', pdf_fn))

        # --- CORRECCIÓN DATABASE CSV (Precisión total) ---
        db_path = os.path.join('Database', 'VALLY_Scan_Master.csv')
        file_exists = os.path.isfile(db_path)
        with open(db_path, 'a', newline='') as f:
            writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            if not file_exists:
                writer.writerow(['Timestamp', 'PDB_ID', 'Pearson_R', 'Hotspots', 'Language'])
            # Convertimos la lista de hotspots a string para evitar errores de delimitador
            str_hotspots = "-".join(map(str, hotspots))
            writer.writerow([datetime.datetime.now().strftime("%Y-%m-%d %H:%M"), pdb_id, round(r_val, 4), str_hotspots, lang])

        print(f"--> [SUCCESS] {pdb_id} processed in {lang.upper()}")

    except Exception as e:
        print(f"--> [CRITICAL ERROR] {str(e)}")

# Ejemplo de ejecución doble para validación total:
# run_vally_scan('6LU7.pdb', lang='es')
# run_vally_scan('6LU7.pdb', lang='en')
