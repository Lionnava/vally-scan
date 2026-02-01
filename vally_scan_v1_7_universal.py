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

# --- 1. CONFIGURACIÓN ESTRUCTURAL ---
def setup_vally_environment():
    """Crea la jerarquía de directorios para una gestión profesional de datos"""
    folders = ['Input_PDB', 'Reports', 'Plots', 'Database']
    for folder in folders:
        if not os.path.exists(folder):
            os.makedirs(folder)

# --- 2. DISEÑO DE REPORTE PREMIUM ---
class VALLY_Premium_Report(FPDF):
    def __init__(self, pdb_id):
        super().__init__()
        self.logo_path = 'logo_proyecto.png'
        self.pdb_id = pdb_id

    def header(self):
        # Encabezado Industrial Azul Oxford
        self.set_fill_color(0, 32, 63)
        self.rect(0, 0, 210, 45, 'F')
        
        logo_x = 12
        if os.path.exists(self.logo_path):
            try:
                self.image(self.logo_path, 12, 10, 28)
                logo_x = 45
            except: pass
            
        self.set_xy(logo_x, 12)
        self.set_text_color(255, 255, 255)
        self.set_font("Helvetica", 'B', 24)
        self.cell(0, 12, "V.A.L.L.Y. PROJECT", ln=True)
        
        self.set_font("Helvetica", '', 10)
        self.set_x(logo_x)
        self.cell(0, 5, f"TECHNICAL DATA SHEET: {self.pdb_id} | Framework v1.7.3", ln=True)

    def footer(self):
        self.set_y(-25)
        self.set_draw_color(0, 32, 63)
        self.line(10, self.get_y(), 200, self.get_y())
        self.set_font("Helvetica", 'I', 8)
        self.set_text_color(100, 100, 100)
        self.cell(0, 10, f"Propiedad Intelectual de Lionell E. Nava Ramos | UPTMA - 2026", align='C')

    def add_security_watermark(self):
        # Protección de datos confidenciales
        self.set_font('Helvetica', 'B', 40)
        self.set_text_color(248, 248, 248)
        self.set_xy(0, 140)
        self.cell(210, 20, "VALLY PROJECT - CONFIDENTIAL DATA", 0, 0, 'C')

# --- 3. FUNCIONES TÉCNICAS ---
def get_system_specs():
    """Captura el hardware de ejecución para trazabilidad científica"""
    cpu = platform.processor() or "Intel64 Family 6"
    ram = f"{round(psutil.virtual_memory().total / (1024**3))} GB RAM"
    return cpu, ram

def update_vally_database(data_row):
    """Actualiza el archivo maestro de inteligencia de datos"""
    db_path = os.path.join('Database', 'VALLY_Master_Database.csv')
    file_exists = os.path.isfile(db_path)
    with open(db_path, mode='a', newline='') as f:
        writer = csv.writer(f)
        if not file_exists:
            writer.writerow(['Timestamp', 'PDB_ID', 'Pearson_R', 'Top_Hotspots', 'CPU', 'RAM'])
        writer.writerow(data_row)

# --- 4. MOTOR DE ANÁLISIS UNIVERSAL ---
def vally_universal_engine(pdb_file, active_site_residues=None, mode='universal'):
    """
    Motor principal V.A.L.L.Y. v1.7.3
    Soporta búsqueda inteligente de archivos y registro automático.
    """
    try:
        setup_vally_environment()
        
        # Lógica de búsqueda de archivo (Carpeta Input o Raíz)
        path_in_folder = os.path.join('Input_PDB', pdb_file)
        target_path = path_in_folder if os.path.exists(path_in_folder) else pdb_file
        
        if not os.path.exists(target_path):
            print(f"\n[!] ERROR: No se encuentra '{pdb_file}' en /Input_PDB/ ni en la raiz.")
            return

        # --- FASE 1: PROCESAMIENTO FÍSICO ---
        structure = parsePDB(target_path)
        calpha = structure.select('protein and name CA')
        anm = ANM(pdb_file)
        anm.buildHessian(calpha)
        anm.calcModes(n_modes=20)
        
        msf_predicted = calcSqFlucts(anm)
        b_factors = calpha.getBetas()
        
        # Validación Estadística
        r_val, _ = pearsonr(msf_predicted, b_factors)
        top_indices = np.argsort(msf_predicted)[-5:][::-1]

        # --- FASE 2: GENERACIÓN DE GRÁFICO ---
        plt.figure(figsize=(10, 5))
        plt.plot(msf_predicted, label='VALLY Simulation', color='#00203F', lw=2)
        plt.plot(b_factors / np.max(b_factors) * np.max(msf_predicted), 
                 label='Experimental (B-factors)', color='#32CD32', linestyle='--', alpha=0.6)
        plt.title(f'V.A.L.L.Y. Scan | {pdb_file} | r = {round(r_val, 3)}')
        plt.legend()
        plt.grid(True, alpha=0.2)
        
        plot_name = f"Validacion_{pdb_file.replace('.pdb','')}.png"
        plt.savefig(os.path.join('Plots', plot_name))
        
        # --- FASE 3: GENERACIÓN DE REPORTE PDF ---
        pdf = VALLY_Premium_Report(pdb_file.upper())
        pdf.add_page()
        pdf.add_security_watermark()
        
        # Sección I: Métricas
        pdf.set_xy(10, 55)
        pdf.set_font("Helvetica", 'B', 14)
        pdf.set_text_color(0, 32, 63)
        pdf.cell(0, 10, "I. SCIENTIFIC METRICS & VALIDATION", ln=True)
        pdf.set_font("Helvetica", '', 11)
        pdf.cell(0, 8, f"Target PDB: {pdb_file.upper()} | Pearson r: {round(r_val, 3)}", ln=True)
        
        # Sección II: Hardware
        cpu, ram = get_system_specs()
        pdf.ln(5)
        pdf.set_font("Helvetica", 'B', 14)
        pdf.cell(0, 10, "II. COMPUTATIONAL ENVIRONMENT", ln=True)
        pdf.set_font("Helvetica", '', 10)
        pdf.multi_cell(0, 6, f"Processor: {cpu}\nMemory: {ram}\nSoftware Core: VALLY v1.7.3 Engine")

        # Sección III: Hotspots
        pdf.ln(5)
        pdf.set_font("Helvetica", 'B', 14)
        pdf.cell(0, 10, "III. TOP VIBRATIONAL HOTSPOTS", ln=True)
        for i, res in enumerate(top_indices, 1):
            pdf.cell(0, 7, f"#{i} Residue Index: {res} | Potential Allosteric Site", ln=True)

        report_name = f"VALLY_Report_{pdb_file.replace('.pdb','')}.pdf"
        pdf.output(os.path.join('Reports', report_name))
        
        # --- FASE 4: REGISTRO EN BASE DE DATOS ---
        update_vally_database([datetime.datetime.now(), pdb_file, round(r_val, 4), list(top_indices), cpu, ram])
        
        print(f"\n--> [EXITO] Escaneo de {pdb_file} completado.")
        print(f"--> [DATOS] Reporte: Reports/{report_name}")
        print(f"--> [DATOS] Grafico: Plots/{plot_name}")
        print(f"--> [DATOS] Registro: Database/VALLY_Master_Database.csv")
        
        plt.show()

    except Exception as e:
        print(f"\n[!] ERROR EN EL SISTEMA: {str(e)}")
