# ===================================================================
# PROYECTO V.A.L.L.Y. - VALLY-Scan v1.4 (Reporte Avanzado)
# Autor: Lionell Eduardo Nava Ramos
# ===================================================================

# --- 0. Importar Librerías ---
import os
import argparse
import numpy as np
from prody import *
from datetime import datetime

# --- Bloque de Importación Controlada para Reportes ---
try:
    from reportlab.pdfgen import canvas
    from reportlab.lib.pagesizes import letter
    from reportlab.lib.units import inch
    from reportlab.lib import colors
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
    from reportlab.lib.styles import getSampleStyleSheet
    REPORTLAB_DISPONIBLE = True
except ImportError:
    REPORTLAB_DISPONIBLE = False

# --- 1. Definición de Funciones Clave (Sin cambios) ---
def calcular_y_guardar_anm(pdb_path, pdb_id):
    # (El código de esta función es idéntico al anterior)
    print(f"\n--- PASO 1: Iniciando análisis ENM para {pdb_id} ---")
    try:
        full_structure = parsePDB(pdb_path)
        protein_ca = full_structure.select('protein and name CA')
        if protein_ca is None or protein_ca.numAtoms() == 0:
            print("ERROR: No se encontraron Carbonos Alfa en la estructura.")
            return None, None
        print(f"Estructura cargada. {protein_ca.numAtoms()} Cα seleccionados.")
        
        anm = ANM(f'{pdb_id} ANM')
        anm.buildHessian(protein_ca)
        anm.calcModes(n_modes=20, zeros=True)
        print("Modelo ANM calculado exitosamente.")
        return anm, protein_ca
        
    except Exception as e:
        print(f"ERROR CRÍTICO durante el análisis ANM: {e}")
        return None, None

def predecir_con_ia_simulada(anm_model, protein_ca, residuos_clave_nums):
    # (El código de esta función es idéntico al anterior)
    print("\n--- PASO 2: Ejecutando el predictor de IA simulado ---")
    modos_relevantes = anm_model[6:9]
    sq_flucts = calcSqFlucts(modos_relevantes)
    puntuacion_impacto = sq_flucts.copy()
    
    if residuos_clave_nums:
        print("Aplicando conocimiento de la literatura (bonus a residuos clave)...")
    resnums = protein_ca.getResnums()
    
    for res_num in residuos_clave_nums:
        try:
            array_index = np.where(resnums == res_num)[0][0]
            puntuacion_impacto[array_index] *= 1.5
        except IndexError:
            pass
            
    afinidad_predicha = -7.5 - np.log(np.sum(puntuacion_impacto)) - np.random.rand()
    indices_top_5 = np.argsort(puntuacion_impacto)[-5:]
    
    residuos_predichos = []
    for index in reversed(indices_top_5):
        residuo = protein_ca[index]
        residuos_predichos.append(f"{residuo.getResname()} {residuo.getResnum()}")
        
    print("¡Predicción completada!")
    return afinidad_predicha, residuos_predichos

# --- 2. NUEVA Función para Generar Reporte PDF Avanzado ---
def generar_reporte_pdf(pdb_id, afinidad, top_residuos, nombre_proteina, directorio_proyecto):
    if not REPORTLAB_DISPONIBLE:
        print("\nADVERTENCIA: No se puede generar el reporte PDF. La librería 'reportlab' no está instalada.")
        return

    reportes_folder = os.path.join(directorio_proyecto, "reportes")
    os.makedirs(reportes_folder, exist_ok=True)
    
    fecha_actual = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    pdf_path = os.path.join(reportes_folder, f"Reporte_VALLY_{pdb_id}_{fecha_actual}.pdf")
    
    c = canvas.Canvas(pdf_path, pagesize=letter)
    width, height = letter

    # --- Función auxiliar para dibujar texto y manejar saltos de línea ---
    def draw_wrapped_text(text, x, y, max_width, font, size, leading):
        from reportlab.platypus import Paragraph
        from reportlab.lib.styles import ParagraphStyle
        style = ParagraphStyle(name='Normal', fontName=font, fontSize=size, leading=leading)
        p = Paragraph(text, style)
        p.wrapOn(c, max_width, height)
        p.drawOn(c, x, y - p.height)

    # --- Encabezado ---
    c.setFont("Helvetica-Bold", 24)
    c.drawString(inch, height - inch, "Proyecto V.A.L.L.Y.")
    c.setFont("Helvetica-Oblique", 12)
    c.drawString(inch, height - 1.25 * inch, "Análisis Vibracional para la Probabilidad del Ligando")
    c.setStrokeColorRGB(0.1, 0.2, 0.5) # Color azul oscuro para la línea
    c.setLineWidth(2)
    c.line(inch, height - 1.4 * inch, width - inch, height - 1.4 * inch)

    # --- Sección 1: Resumen del Análisis (Columna Izquierda) ---
    c.setFont("Helvetica-Bold", 14)
    c.drawString(inch, height - 1.8 * inch, "Resumen del Análisis")
    
    text_y = height - 2.1 * inch
    c.setFont("Helvetica-Bold", 10)
    c.drawString(inch, text_y, "Estructura (PDB ID):")
    c.setFont("Helvetica", 10)
    c.drawString(inch + 1.5 * inch, text_y, pdb_id.upper())
    
    text_y -= 0.3 * inch
    c.setFont("Helvetica-Bold", 10)
    c.drawString(inch, text_y, "Objetivo Terapéutico:")
    c.setFont("Helvetica", 10)
    draw_wrapped_text(nombre_proteina, inch + 1.5 * inch, text_y, 2.5*inch, "Helvetica", 10, 12)

    text_y -= 0.6 * inch # Más espacio por el texto envuelto
    c.setFont("Helvetica-Bold", 10)
    c.drawString(inch, text_y, "Fecha de Generación:")
    c.setFont("Helvetica", 10)
    c.drawString(inch + 1.5 * inch, text_y, datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

    text_y -= 0.3 * inch
    c.setFont("Helvetica-Bold", 10)
    c.drawString(inch, text_y, "Método Utilizado:")
    c.setFont("Helvetica", 10)
    draw_wrapped_text("Modelo de Red Elástica (ANM) + Predictor Heurístico Informado", inch + 1.5 * inch, text_y, 2.5*inch, "Helvetica", 10, 12)
    
    text_y -= 0.6 * inch
    c.setFont("Helvetica-Bold", 10)
    c.drawString(inch, text_y, "Hardware Utilizado:")
    c.setFont("Helvetica", 10)
    c.drawString(inch + 1.5 * inch, text_y, "CPU de Escritorio Estándar")

    # --- Sección 2: Resultados de la Predicción (Columna Derecha) ---
    col_2_x = inch + 4.5 * inch
    c.setFont("Helvetica-Bold", 14)
    c.drawString(col_2_x, height - 1.8 * inch, "Resultados de la Predicción")
    
    # Afinidad Predicha
    c.setFont("Helvetica-Bold", 11)
    c.drawString(col_2_x, height - 2.1 * inch, "Afinidad de Enlace Predicha:")
    c.setFont("Helvetica-Bold", 36)
    c.setFillColor(colors.firebrick)
    c.drawRightString(width - inch, height - 2.6 * inch, f"{afinidad:.2f}")
    c.setFont("Helvetica", 14)
    c.setFillColor(colors.darkgrey)
    c.drawString(width - 1.7 * inch, height - 2.85 * inch, "kcal/mol")

    # Residuos Clave
    c.setFillColor(colors.black)
    c.setFont("Helvetica-Bold", 11)
    c.drawString(col_2_x, height - 3.3 * inch, "Top 5 Residuos de Mayor Impacto:")
    
    text_y = height - 3.6 * inch
    c.setFont("Helvetica", 12)
    for i, res in enumerate(top_residuos):
        c.drawString(col_2_x + 0.2*inch, text_y, f"{i+1}.   {res}")
        text_y -= 0.25 * inch

    # --- Línea Divisoria Vertical ---
    c.setStrokeColorRGB(0.8, 0.8, 0.8)
    c.setLineWidth(1)
    c.line(inch + 4 * inch, height - 1.6 * inch, inch + 4 * inch, inch * 1.5)

    # --- Pie de Página ---
    c.setStrokeColorRGB(0.1, 0.2, 0.5)
    c.setLineWidth(2)
    c.line(inch, 1.4 * inch, width - inch, 1.4 * inch)
    c.setFont("Helvetica-Oblique", 8)
    c.drawCentredString(width / 2, 1.1 * inch, "Este es un reporte predictivo generado automáticamente para propósitos de investigación. Los resultados deben ser validados experimentalmente.")
    c.setFont("Helvetica", 9)
    c.drawCentredString(width / 2, 0.8 * inch, f"Reporte generado por VALLY-Scan v1.4 | Autor: Lionell E. Nava Ramos")

    c.save()
    print(f"\n¡Reporte PDF Avanzado generado exitosamente!\nArchivo guardado en: {pdf_path}")

# --- 3. Lógica Principal y Manejo de Argumentos (Sin cambios) ---
def main():
    # ... (El código de esta función es idéntico a la v1.3)
    parser = argparse.ArgumentParser(
        description="VALLY-Scan v1.4: Herramienta para el análisis rápido de interacciones ligando-proteína.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--pdb_id", type=str, required=True, help="El ID de 4 letras de la estructura en el PDB. Ej: 2fom")
    parser.add_argument("--reporte", action="store_true", help="Genera un reporte final del análisis en formato PDF.")
    args = parser.parse_args()
    
    pdb_id = args.pdb_id.lower()
    directorio_proyecto = os.getcwd()
    data_folder = "data"
    pdb_path = os.path.join(directorio_proyecto, data_folder, f"{pdb_id}.pdb")

    if not os.path.exists(pdb_path):
        print(f"ERROR: No se encontró el archivo '{pdb_id}.pdb' en la carpeta '{data_folder}'.")
        return

    CONOCIMIENTO_EXPERTO = {
        "2fom": { "nombre": "Proteasa del Virus del Dengue (NS2B/NS3)", "residuos_clave": [51, 75, 133, 135, 150, 151, 153] },
        "6lu7": { "nombre": "Proteasa Principal del SARS-CoV-2 (Mpro)", "residuos_clave": [41, 49, 143, 145, 163, 164, 166, 189] }
    }

    if pdb_id in CONOCIMIENTO_EXPERTO:
        info_proteina = CONOCIMIENTO_EXPERTO[pdb_id]
        residuos_clave = info_proteina["residuos_clave"]
        nombre_proteina = info_proteina["nombre"]
        print(f"\nInformación experta para '{nombre_proteina}' cargada.")
    else:
        residuos_clave = []; nombre_proteina = "Desconocida"
        print(f"\nADVERTENCIA: No se ha definido una lista de residuos clave para '{pdb_id}'.")

    anm, protein_ca = calcular_y_guardar_anm(pdb_path, pdb_id)
    
    if anm and protein_ca:
        afinidad, top_residuos = predecir_con_ia_simulada(anm, protein_ca, residuos_clave)
        
        print("\n" + "="*50); print("¡ANÁLISIS COMPLETO CON VALLY-SCAN!"); print("="*50)
        print(f"Análisis para la estructura: {pdb_id.upper()} ({nombre_proteina})")
        print(f"Afinidad de Enlace Predicha: {afinidad:.2f} kcal/mol")
        print("\nTop 5 Residuos más Impactantes en la Interacción:")
        for i, res in enumerate(top_residuos): print(f"  {i+1}. {res}")
        print("="*50)
        
        if args.reporte:
            generar_reporte_pdf(pdb_id, afinidad, top_residuos, nombre_proteina, directorio_proyecto)

if __name__ == "__main__":
    main()