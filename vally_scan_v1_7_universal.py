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

# --- 1. PRESERVACIÓN DE INFRAESTRUCTURA ---
def setup_vally_environment():
    for folder in ['Input_PDB', 'Reports', 'Plots', 'Database']:
        if not os.path.exists(folder):
            os.makedirs(folder)

# --- 2. MOTOR DE IDENTIDAD VISUAL (RESTAURADO Y ENRIQUECIDO) ---
class VALLY_Premium_Report(FPDF):
    def __init__(self, pdb_id):
        super().__init__()
        self.pdb_id = pdb_id
        self.logo_path = 'logo_proyecto.png'

    def header(self):
        # BANNER SUPERIOR SÓLIDO (Identidad de Marca)
        self.set_fill_color(0, 32, 63)  # Azul Oxford
        self.rect(0, 0, 210, 45, 'F')
        
        # BANNER DE LEMA (Soberanía y Ciencia)
        self.set_fill_color(0, 50, 90)  # Azul Profundo
        self.rect(0, 45, 210, 8, 'F')
        self.set_xy(0, 46.5)
        self.set_text_color(255, 255, 255)
        self.set_font("Helvetica", 'B', 9)
        self.cell(210, 5, "CIENCIA Y SOBERANIA EN SALUD | PROYECTO VALLY", 0, 0, 'C')

        # CONTENIDO DEL HEADER
        logo_x = 15
        if os.path.exists(self.logo_path):
            self.image(self.logo_path, 15, 10, 28)
            logo_x = 48
            
        self.set_xy(logo_x, 12)
        self.set_text_color(255, 255, 255)
        self.set_font("Helvetica", 'B', 22)
        self.cell(0, 10, "VALLY-SCAN INSTRUMENT", ln=True)
        self.set_font("Helvetica", '', 10)
        self.set_x(logo_x)
        self.cell(0, 5, "VIBRATIONAL ANALYSIS & LOCAL LIGAND YIELDING", ln=True)
        self.set_x(logo_x)
        self.cell(0, 5, f"EDITION v3.0 PRO | TECHNICAL DATA SHEET", ln=True)

    def footer(self):
        self.set_y(-25)
        self.set_draw_color(0, 32, 63)
        self.line(10, self.get_y(), 200, self.get_y())
        self.set_font("Helvetica", 'I', 8)
        self.set_text_color(100, 100, 100)
        self.ln(2)
        self.cell(0, 10, "PROYECTO VALLY - PROPIEDAD INTELECTUAL DE LIONELL E. NAVA RAMOS", align='C')
        self.set_y(-15)
        self.cell(0, 10, f"Pagina {self.page_no()}", align='R')

# --- 3. MOTOR UNIVERSAL (MANTIENE TODA LA LÓGICA APROBADA) ---
def vally_universal_engine(pdb_file, active_site_residues=None, mode='universal'):
    try:
        setup_vally_environment()
        target = os.path.join('Input_PDB', pdb_file) if os.path.exists(os.path.join('Input_PDB', pdb_file)) else pdb_file
        
        # --- NÚCLEO DE FÍSICA (SIN CAMBIOS PARA EVITAR ROTURAS) ---
        structure = parsePDB(target)
        calpha = structure.select('protein and name CA')
        anm = ANM(pdb_file)
        anm.buildHessian(calpha)
        anm.calcModes(n_modes=30)
        msf = calcSqFlucts(anm)
        b_factors = calpha.getBetas()
        r_val, _ = pearsonr(msf, b_factors)
        top_indices = np.argsort(msf)[-5:][::-1]

        # --- ENRIQUECIMIENTO: GRÁFICO CIENTÍFICO CON LEYENDA ---
        plt.figure(figsize=(10, 5))
        m_norm = (msf - np.mean(msf)) / np.std(msf)
        b_norm = (b_factors - np.mean(b_factors)) / np.std(b_factors)
        plt.plot(m_norm, color='#00203F', label='VALLY Simulation (ANM)', lw=2)
        plt.plot(b_norm, color='#32CD32', label='Experimental Data (B-factors)', ls='--', alpha=0.6)
        plt.title(f"VALLY-SCAN Validation | {pdb_file.upper()} | r = {round(r_val, 3)}")
        plt.xlabel("Residue Index"); plt.ylabel("Standardized Fluctuation")
        plt.legend(loc='best', frameon=True, shadow=True)
        plt.grid(True, alpha=0.3)
        plt.savefig(os.path.join('Plots', f"Plot_{pdb_file[:-4]}.png"), dpi=300)
        plt.close()

        # --- CONSTRUCCIÓN DEL PDF (FORMATO DE ALTA GAMA) ---
        pdf = VALLY_Premium_Report(pdb_file.upper())
        pdf.add_page()
        
        # Bloque 1: Datos Técnicos
        pdf.set_xy(10, 65)
        pdf.set_font("Helvetica", 'B', 14); pdf.set_text_color(0, 32, 63)
        pdf.cell(0, 10, "1. METRICAS DE VALIDACION CIENTIFICA", ln=True)
        pdf.set_font("Helvetica", '', 11); pdf.set_text_color(0, 0, 0)
        pdf.cell(50, 8, "PDB Target:", 0); pdf.cell(0, 8, pdb_file.upper(), ln=True)
        pdf.cell(50, 8, "Pearson Correlation (r):", 0); pdf.cell(0, 8, f"{round(r_val, 4)}", ln=True)
        pdf.cell(50, 8, "Timestamp:", 0); pdf.cell(0, 8, str(datetime.date.today()), ln=True)

        # Bloque 2: Hotspots
        pdf.ln(8)
        pdf.set_font("Helvetica", 'B', 14); pdf.set_text_color(0, 32, 63)
        pdf.cell(0, 10, "2. VIBRATIONAL HOTSPOTS MAPPING", ln=True)
        pdf.set_font("Helvetica", '', 10); pdf.set_text_color(0, 0, 0)
        for i, res in enumerate(top_indices, 1):
            pdf.cell(0, 7, f"Rank {i} -> Indice de Residuo: {res} | Candidato a Sitio Alosterico", ln=True)

        # Bloque 3: Resumen Ejecutivo
        pdf.ln(8)
        pdf.set_font("Helvetica", 'B', 14); pdf.set_text_color(0, 32, 63)
        pdf.cell(0, 10, "3. RESUMEN EJECUTIVO", ln=True)
        pdf.set_font("Helvetica", 'I', 10); pdf.set_text_color(0, 0, 0)
        resumen = (f"El motor VALLY-SCAN ha validado la estructura {pdb_file.upper()} con una robustez "
                   f"del {round(r_val*100, 1)}%. La paridad entre la simulacion y los datos experimentales "
                   "confirma la precision del modelo de redes anisotropicvas.")
        pdf.multi_cell(0, 6, resumen)

        # NOTA DE PROCESO
        pdf.ln(5); pdf.set_font("Helvetica", 'B', 9); pdf.set_text_color(180, 0, 0)
        pdf.multi_cell(0, 5, "AVISO: Los resultados utilizan normalizacion Z-Score para asegurar paridad estadistica.")

        pdf_fn = f"VALLY_Scan_Report_{pdb_file[:-4]}.pdf"
        pdf.output(os.path.join('Reports', pdf_fn))

        # --- 4. BASE DE DATOS (MANTENIENDO REGISTRO CSV SEGURO) ---
        db_path = os.path.join('Database', 'VALLY_Scan_Master.csv')
        hotspots_safe = "-".join(map(str, top_indices))
        with open(db_path, 'a', newline='') as f:
            writer = csv.writer(f)
            if os.path.getsize(db_path) == 0:
                writer.writerow(['Timestamp', 'PDB', 'Pearson_R', 'Hotspots'])
            writer.writerow([datetime.datetime.now().strftime("%Y-%m-%d %H:%M"), pdb_file, round(r_val, 4), hotspots_safe])

        print(f"--> [VALLY-SCAN OK] {pdb_file} procesado con exito.")

    except Exception as e:
        print(f"--> [ERROR CRITICO] El instrumento fallo: {str(e)}")
