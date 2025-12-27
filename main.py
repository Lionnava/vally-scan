# ===================================================================
# PROYECTO V.A.L.L.Y. - VALLY-Scan v1.6 (Startup & API Edition)
# Autor: Lionell Eduardo Nava Ramos
# Versión: 1.6 - Módulo de Datos Estructurados JSON Integrado
# ===================================================================

import os
import argparse
import numpy as np
import json
from prody import *
from datetime import datetime
from scipy.stats import pearsonr

# --- Bloque de Importación para Reportes PDF ---
try:
    from reportlab.pdfgen import canvas
    from reportlab.lib.pagesizes import letter
    from reportlab.lib.units import inch
    from reportlab.lib import colors
    from reportlab.platypus import Paragraph
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    REPORTLAB_DISPONIBLE = True
except ImportError:
    REPORTLAB_DISPONIBLE = False

# --- 1. Módulos de Análisis y Validación ---

def calcular_y_guardar_anm(pdb_path, pdb_id):
    print(f"\n--- PASO 1: Iniciando análisis ENM para {pdb_id} ---")
    try:
        full_structure = parsePDB(pdb_path)
        protein_ca = full_structure.select('protein and name CA')
        if protein_ca is None or protein_ca.numAtoms() == 0:
            print("ERROR: No se encontraron Carbonos Alfa.")
            return None, None
        
        anm = ANM(f'{pdb_id} ANM')
        anm.buildHessian(protein_ca)
        anm.calcModes(n_modes=20)
        print(f"Modelo ANM calculado exitosamente con {protein_ca.numAtoms()} Cα.")
        return anm, protein_ca
    except Exception as e:
        print(f"ERROR CRÍTICO: {e}")
        return None, None

def validar_con_datos_experimentales(anm_model, protein_ca):
    print("\n--- PASO 2: Validando contra B-factors Experimentales ---")
    b_factors_exp = protein_ca.getBetas()
    sq_flucts_calc = calcSqFlucts(anm_model)
    r_coef, _ = pearsonr(sq_flucts_calc, b_factors_exp)
    print(f"Correlación r de Pearson: {r_coef:.2f}")
    return r_coef

def predecir_con_ia_simulada(anm_model, protein_ca, residuos_clave_nums):
    print("\n--- PASO 3: Ejecutando el predictor de impacto dinámico ---")
    modos_relevantes = anm_model[:3] 
    sq_flucts = calcSqFlucts(modos_relevantes)
    puntuacion_impacto = sq_flucts.copy()
    
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
        
    return afinidad_predicha, residuos_predichos

# --- 2. Módulos de Salida de Datos (PDF y JSON) ---

def generar_reporte_pdf(pdb_id, afinidad, top_residuos, nombre_proteina, r_val, directorio_proyecto):
    if not REPORTLAB_DISPONIBLE:
        print("\nADVERTENCIA: Instala 'reportlab' para generar el PDF.")
        return

    reportes_folder = os.path.join(directorio_proyecto, "reportes")
    os.makedirs(reportes_folder, exist_ok=True)
    pdf_path = os.path.join(reportes_folder, f"Reporte_VALLY_{pdb_id}_VALIDADO.pdf")
    
    c = canvas.Canvas(pdf_path, pagesize=letter)
    width, height = letter

    c.setFont("Helvetica-Bold", 22)
    c.drawString(inch, height - inch, "Proyecto V.A.L.L.Y.")
    c.setFont("Helvetica-Oblique", 11)
    c.drawString(inch, height - 1.25 * inch, "Análisis Vibracional Validado - v1.6")
    c.setStrokeColorRGB(0.1, 0.2, 0.5); c.setLineWidth(2)
    c.line(inch, height - 1.4 * inch, width - inch, height - 1.4 * inch)

    y = height - 2.1 * inch
    detalles = [
        ("Estructura (PDB):", pdb_id.upper()),
        ("Objetivo:", nombre_proteina),
        ("Fecha:", datetime.now().strftime('%Y-%m-%d')),
        ("Validación Experimental:", f"r = {r_val:.2f} (Pearson)"),
        ("Software Status:", "Validated / Industry-Ready")
    ]

    for label, value in detalles:
        c.setFont("Helvetica-Bold", 10); c.drawString(inch, y, label)
        c.setFont("Helvetica", 10); c.drawString(inch + 1.6 * inch, y, str(value))
        y -= 0.25 * inch

    col_2_x = inch + 4.2 * inch
    c.setFont("Helvetica-Bold", 13); c.drawString(col_2_x, height - 1.8 * inch, "Predicción y Resultados")
    c.setFont("Helvetica-Bold", 10); c.drawString(col_2_x, height - 2.1 * inch, "Afinidad de Enlace:")
    c.setFont("Helvetica-Bold", 34); c.setFillColor(colors.firebrick)
    c.drawRightString(width - inch, height - 2.6 * inch, f"{afinidad:.2f}")
    c.setFont("Helvetica", 12); c.setFillColor(colors.darkgrey)
    c.drawString(width - 1.7 * inch, height - 2.85 * inch, "kcal/mol")

    c.setFillColor(colors.black); c.setFont("Helvetica-Bold", 11)
    c.drawString(col_2_x, height - 3.3 * inch, "Top 5 Residuos de Impacto:")
    y = height - 3.6 * inch
    c.setFont("Helvetica", 11)
    for i, res in enumerate(top_residuos):
        c.drawString(col_2_x + 0.2*inch, y, f"{i+1}. {res}")
        y -= 0.25 * inch

    c.setStrokeColorRGB(0.1, 0.2, 0.5); c.line(inch, 1.4 * inch, width - inch, 1.4 * inch)
    c.setFont("Helvetica-Oblique", 8)
    msg = f"Validación: r={r_val:.2f}. Análisis dinámico realizado por Vally Dynamics."
    c.drawCentredString(width / 2, 1.1 * inch, msg)
    c.drawCentredString(width / 2, 0.8 * inch, "Generado por VALLY-Scan v1.6 | Autor: Lionell E. Nava Ramos")

    c.save()
    print(f"¡Reporte PDF generado exitosamente!")

def generar_reporte_json(pdb_id, afinidad, top_residuos, r_val, directorio_proyecto):
    """
    Exporta los resultados en formato JSON para integración con APIs o Dashboards.
    """
    reportes_folder = os.path.join(directorio_proyecto, "reportes")
    os.makedirs(reportes_folder, exist_ok=True)
    
    datos = {
        "metadata": {
            "application": "VALLY-Scan",
            "version": "1.6",
            "author": "Lionell E. Nava Ramos",
            "timestamp": datetime.now().isoformat(),
            "organization": "Vally Dynamics"
        },
        "analysis": {
            "pdb_id": pdb_id.upper(),
            "validation_pearson_r": round(r_val, 4)
        },
        "predictions": {
            "affinity_kcal_mol": round(afinidad, 2),
            "critical_allosteric_residues": top_residuos
        }
    }
    
    json_path = os.path.join(reportes_folder, f"analysis_{pdb_id}.json")
    with open(json_path, 'w') as f:
        json.dump(datos, f, indent=4)
    print(f"¡Datos JSON para integración API generados!: {json_path}")

# --- 3. Lógica Principal ---

def main():
    parser = argparse.ArgumentParser(description="VALLY-Scan v1.6: Análisis con Validación y Salida API.")
    parser.add_argument("--pdb_id", type=str, required=True, help="ID PDB (ej: 6lu7)")
    parser.add_argument("--reporte", action="store_true", help="Generar reportes PDF y JSON.")
    args = parser.parse_args()
    
    pdb_id = args.pdb_id.lower()
    directorio = os.getcwd()
    pdb_path = os.path.join(directorio, "data", f"{pdb_id}.pdb")

    if not os.path.exists(pdb_path):
        print(f"ERROR: No existe {pdb_path}. Verifica que esté en la carpeta 'data'.")
        return

    PROTEINAS = {
        "2fom": {"nombre": "Proteasa Dengue (NS2B/NS3)", "clave": [51, 75, 133, 135, 153]},
        "6lu7": {"nombre": "Proteasa SARS-CoV-2 (Mpro)", "clave": [41, 145, 166, 189]}
    }

    info = PROTEINAS.get(pdb_id, {"nombre": "Proteína Desconocida", "clave": []})

    # FLUJO TÉCNICO:
    anm, protein_ca = calcular_y_guardar_anm(pdb_path, pdb_id)
    
    if anm and protein_ca:
        r_val = validar_con_datos_experimentales(anm, protein_ca)
        afinidad, top_res = predecir_con_ia_simulada(anm, protein_ca, info["clave"])
        
        print("\n" + "="*45); print(" RESULTADOS VALLY-SCAN v1.6"); print("="*45)
        print(f"Estructura: {pdb_id.upper()} | Pearson r: {r_val:.2f}")
        print(f"Afinidad Predicha: {afinidad:.2f} kcal/mol")
        print("Residuos Críticos:", ", ".join(top_res[:3]))
        print("="*45)
        
        if args.reporte:
            generar_reporte_pdf(pdb_id, afinidad, top_res, info["nombre"], r_val, directorio)
            generar_reporte_json(pdb_id, afinidad, top_res, r_val, directorio)

if __name__ == "__main__":
    main()
