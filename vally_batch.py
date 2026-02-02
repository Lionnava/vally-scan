import os
# Importamos su motor sin modificarlo
from vally_scan_v1_7_universal import vally_universal_engine

def run_full_inventory():
    pdb_folder = 'Input_PDB'
    # Listamos los archivos sin alterar nada
    pdb_files = [f for f in os.listdir(pdb_folder) if f.endswith('.pdb')]
    
    print(f"--- INICIANDO PROCESAMIENTO R2 (Total: {len(pdb_files)} archivos) ---")
    
    for pdb in pdb_files:
        try:
            # Llamada segura al motor que ya validamos
            vally_universal_engine(pdb)
        except Exception as e:
            print(f"--> [AVISO] Saltando {pdb} por error de formato: {e}")

if __name__ == "__main__":
    run_full_inventory()
