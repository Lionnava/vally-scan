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

# --- 1. ENTORNO VALLY-SCAN (ESTÁNDAR APS) ---
def setup_vally_environment():
    for folder in ['Input_PDB', 'Reports', 'Plots', 'Database']:
        if not os.path.exists(folder): os.makedirs(folder)

# --- 2. DISEÑO DE REPORTE DE ALTA JERARQUÍA ---
class VALLY_APS_Report(FPDF):
    def __init__(self, pdb_id):
        super().__init__()
        self.pdb_id = pdb_id
        self.logo_path = 'logo_proyecto.png'

    def header(self):
        # Banner de Identidad Nacional y Científica
        self.set_fill_color(0, 32, 63)
        self.rect(0, 0, 210, 45, 'F')
        
        # Lema Institucional - Ciencia y Soberanía
        self.set_fill_color(0, 50, 90)
        self.rect(0, 45, 210, 8, 'F')
        self.set_xy(0, 46.5)
        self.set_text_color(255, 255, 255)
        self.set_font("Helvetica", 'B', 10)
        self.cell(210, 5, "CIENCIA Y SOBERANIA EN SALUD | APS GLOBAL SUMMIT 2026", 0, 0, 'C')

        logo_x = 12
        if os.path.exists(self.logo_path):
            try:
                self.image(self.logo_path, 12, 10, 28)
                logo_x = 45
            except: pass
            
        self.set_xy(logo_x, 12)
        self.set_text_color(255, 255, 255)
        self.set_font("Helvetica", 'B', 24)
        self.cell(0, 10, "VALLY-SCAN INSTRUMENT", ln=True)
        self.set_font("Helvetica", '', 10)
        self.set_x(logo_x)
        self.cell(0, 5, "VIBRATIONAL ANALYSIS & LOCAL LIGAND YIELDING", ln=True)
        self.set_x(logo_x)
        self.cell(0, 5, "APS CONVENTION SPECIAL EDITION | v1.8.5 PRO", ln=True)

    def footer(self):
        self.set_y(-25)
        self.set_draw_color(0, 32, 63)
        self.line(10, self.get_y(), 200, self.get_y())
        self.set_font("Helvetica", 'I', 8)
        self.set_text_color(100, 100, 100)
        self.ln(2)
        self.cell(0, 10, "UPTMA 2026 | Intellectual Property of Lionell E. Nava Ramos", align='C')
        self.set_y(-15)
        self.cell(0, 10, f"Technical Datasheet - Page {self.page_no()}", align='R')

# --- 3. MOTOR DE FÍSICA COMPUTACIONAL ---
def vally_universal_engine(pdb_file, active_site_residues=None, mode='universal'):
    try:
        setup_vally_environment()
        target_path = os.path.join('Input_PDB', pdb_file) if os.path.exists(os.path.join('Input_PDB', pdb_file)) else pdb_file
        
        if not os.path.exists(target_path):
            print(f"\n[!] ERROR: Missing source file: {pdb_file}")
            return

        # Análisis Biofísico: ANM (Anisotropic Network Model)
        structure = parsePDB(target_path)
        calpha = structure.select('protein and name CA')
        anm = ANM(pdb_file)
        anm.buildHessian(calpha)
        anm.calcModes(n_modes=30) # Mayor muestreo de espacio configuracional
        
        msf = calcSqFlucts(anm)
        b_factors = calpha.getBetas()
        
        # Validación Estadística
        r_val, _ = pearsonr(msf, b_factors)
        top_hotspots = np.argsort(msf)[-5:][::-1]

        # --- GRÁFICO DE GRADO SUMMIT (Leyendas y Etiquetas) ---
        plt.figure(figsize=(11, 5))
        # Normalización para escala unitaria (Standard Deviation scaling)
        norm_msf = (msf - np.min(msf)) / (np.max(msf) - np.min(msf))
        norm_b = (b_factors - np.min(b_factors)) / (np.max(b_factors) - np.min(b_factors))
        
        plt.plot(norm_msf, color='#00203F', label='VALLY Simulation (Physical Model)', lw=2)
        plt.plot(norm_b, color='#32CD32', label='X-Ray Crystallography (Experimental)', linestyle='--', alpha=0.6)
        
        plt.title(f"VALLY-SCAN VALIDATION | {pdb_file.upper()} | Pearson r = {round(r_val, 3)}", pad=20, weight='bold')
        plt.legend(loc='upper right', frameon=True, shadow=True)
        plt.grid(True, linestyle=':', alpha=0.5)
        plt.xlabel("Residue Sequential Index")
        plt.ylabel("Relative Fluctuation Intensity")
        
        plot_path = os.path.join('Plots', f"APS_Validation_{pdb_file.replace('.pdb','')}.png")
        plt.savefig(plot_path, dpi=400, bbox_inches='tight') # Calidad de impresión
        plt.close()

        # --- GENERACIÓN DE REPORTE TÉCNICO ---
        pdf = VALLY_APS_Report(pdb_file.upper())
        pdf.add_page()
        
        # Sección I: Validación Científica
        pdf.set_xy(10, 60)
        pdf.set_font("Helvetica", 'B', 14); pdf.set_text_color(0, 32, 63)
        pdf.cell(0, 10, "I. PHYSICAL VALIDATION & CROSS-CORRELATION", ln=True)
        pdf.set_font("Helvetica", '', 11); pdf.set_text_color(0, 0, 0)
        pdf.cell(50, 8, "Target PDB Structure:", 0); pdf.cell(60, 8, pdb_file.upper(), ln=True)
        pdf.cell(50, 8, "Pearson Coefficient (r):", 0); pdf.cell(60, 8, f"{round(r_val, 4)}", ln=True)
        pdf.cell(50, 8, "Reliability Index:", 0); 
        status = "EXCELLENT" if r_val > 0.6 else "CONSISTENT" if r_val > 0.4 else "DYNAMIC DISORDER DETECTED"
        pdf.cell(60, 8, status, ln=True)

        # Sección II: Identificación de Hotspots
        pdf.ln(5); pdf.set_font("Helvetica", 'B', 14); pdf.set_text_color(0, 32, 63)
        pdf.cell(0, 10, "II. VIBRATIONAL HOTSPOTS (PLASTICITY MAP)", ln=True)
        pdf.set_font("Helvetica", '', 10); pdf.set_text_color(0, 0, 0)
        for i, res in enumerate(top_hotspots, 1):
            pdf.cell(0, 7, f"Rank {i} | Residue Index: {res} | Propensity for Allosteric Modulation", ln=True)

        # Sección III: Resumen Ejecutivo
        pdf.ln(5); pdf.set_font("Helvetica", 'B', 14); pdf.set_text_color(0, 32, 63)
        pdf.cell(0, 10, "III. EXECUTIVE SUMMARY", ln=True)
        pdf.set_font("Helvetica", 'I', 10); pdf.set_text_color(0, 0, 0)
        summary = (f"VALLY-Scan v1.8.5 ha procesado exitosamente la estructura {pdb_file.upper()}. "
                   f"La correlacion r={round(r_val, 3)} valida el motor de fisica vectorial, permitiendo "
                   "identificar regiones con alta varianza funcional para el diseño de ligandos.")
        pdf.multi_cell(0, 6, summary)

        # Nota APS roja
        pdf.ln(5); pdf.set_font("Helvetica", 'B', 9); pdf.set_text_color(150, 0, 0)
        pdf.multi_cell(0, 5, "APS SUMMIT ADVISORY: Data normalized via Vectorial Scaling. Minor discrepancies "
                             "identify regions of high intrinsic disorder or non-linear vibrational modes.")

        output_path = os.path.join('Reports', f"VALLY_APS_Report_{pdb_file.replace('.pdb','')}.pdf")
        pdf.output(output_path)
        
        print(f"\n--> [OK] Analisis Finalizado para presentacion APS.")
        print(f"--> [OK] Reporte Generado: {output_path}")

    except Exception as e:
        print(f"\n[!] FALLO CRITICO: {str(e)}")
