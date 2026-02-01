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

# --- 1. CONFIGURACIÓN DEL ENTORNO DE TRABAJO ---
def setup_vally_environment():
    """Crea la estructura de carpetas profesional del proyecto"""
    folders = ['Input_PDB', 'Reports', 'Plots', 'Database']
    for folder in folders:
        if not os.path.exists(folder):
            os.makedirs(folder)

# --- 2. CLASE DE DISEÑO PREMIUM (PDF) ---
class VALLY_Premium_Report(FPDF):
    def __init__(self, pdb_id):
        super().__init__()
        self.logo_path = 'logo_proyecto.png'
        self.pdb_id = pdb_id

    def header(self):
        # Fondo Encabezado Azul Oxford
        self.set_fill_color(0, 32, 63)
        self.rect(0, 0, 210, 45, 'F')
        
        logo_existe = False
        if os.path.exists(self.logo_path):
            try:
                self.image(self.logo_path, 12, 10, 28)
                logo_existe = True
            except: pass
            
        self.set_xy(45, 12) if logo_existe else self.set_xy(12, 12)
        self.set_text_color(255, 255, 255)
        self.set_font("Helvetica", 'B', 24)
        self.cell(0, 12, "V.A.L.L.Y. PROJECT", ln=True)
        
        self.set_font("Helvetica", '', 10)
        if logo_existe: self.set_x(45)
        self.cell(0, 5, f"Analysis Report: {self.pdb_id} | Ciencia y Soberania", ln=True)

    def add_security_watermark(self):
        self.set_font('Helvetica', 'B', 40)
        self.set_text_color(248, 248, 248)
        self.set_xy(0, 140)
        self.cell(210, 20, "VALLY PROJECT - CONFIDENTIAL DATA", 0, 0, 'C')

# --- 3. MOTOR UNIVERSAL VALLY ---
def vally_universal_engine(pdb_filename, mode='universal'):
    try:
        setup_vally_environment()
        
        # El sistema busca automáticamente el archivo dentro de /Input_PDB/
        input_path = os.path.join('Input_PDB', pdb_filename)
        
        if not os.path.exists(input_path):
            print(f"[!] ERROR: No se encontro el archivo {pdb_filename} en la carpeta /Input_PDB/")
            return

        # --- PROCESAMIENTO CIENTÍFICO ---
        structure = parsePDB(input_path)
        calpha = structure.select('protein and name CA')
        anm = ANM(pdb_filename)
        anm.buildHessian(calpha)
        anm.calcModes(n_modes=20)
        
        msf_predicted = calcSqFlucts(anm)
        b_factors = calpha.getBetas()
        r_val, _ = pearsonr(msf_predicted, b_factors)
        top_indices = np.argsort(msf_predicted)[-5:][::-1]
        
        # --- GUARDAR GRÁFICA EN /Plots/ ---
        plt.figure(figsize=(10, 5))
        plt.plot(msf_predicted, label='VALLY Simulation', color='#00203F')
        plt.plot(b_factors / np.max(b_factors) * np.max(msf_predicted), 
                 label='Experimental', color='#32CD32', linestyle='--')
        plt.title(f"Validation {pdb_filename} | r = {round(r_val, 3)}")
        plt.savefig(os.path.join('Plots', f"Plot_{pdb_filename.replace('.pdb','')}.png"))
        
        # --- GENERAR REPORTE EN /Reports/ ---
        pdf = VALLY_Premium_Report(pdb_filename.upper())
        pdf.add_page()
        pdf.add_security_watermark()
        # [Secciones de datos omitidas por brevedad, se mantienen igual]
        
        report_path = os.path.join('Reports', f"VALLY_Report_{pdb_filename.replace('.pdb','')}.pdf")
        pdf.output(report_path)
        
        # --- ACTUALIZAR BASE DE DATOS EN /Database/ ---
        cpu = platform.processor()
        db_path = os.path.join('Database', 'VALLY_Master_Database.csv')
        with open(db_path, mode='a', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([datetime.datetime.now(), pdb_filename, round(r_val, 4), list(top_indices), cpu])
        
        print(f"\n--> [OK] Escaneo finalizado para {pdb_filename}")
        print(f"--> [INFO] Resultados en /Reports/ y /Plots/")
        
    except Exception as e:
        print(f"\n[!] ERROR CRITICO: {str(e)}")
