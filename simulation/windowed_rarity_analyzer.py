import os
import hyperloglog
import pandas as pd
import pyfastx
import pickle
import numpy as np
import time
import psutil  

# --- CONFIGURACIÓN ---
INPUT_MUESTRAS_DIR = "muestras_validacion" 
SKETCH_DIR = "hll_sim_catalog"
OUTPUT_RESULTS_DIR = "resultados_performance"

WINDOW_SIZE = 200000 
K_SIZE = 31
HLL_ERROR_RATE = 0.01

# Umbral Dinámico
CALIBRATION_WINDOWS = 5
SIGMA_MULTIPLIER = 2.0

CATALOG_MAP = {
    "DISTINTOS": "distintos.sketch",
    "SIMILARES": "similares.sketch",
    "MEDIOS":    "medios.sketch"
}

def load_hll_sketch(sketch_name):
    path = os.path.join(SKETCH_DIR, sketch_name)
    with open(path, 'rb') as f:
        return pickle.load(f)

def get_canonical_kmer(kmer):
    trans = str.maketrans("ATCGN", "TAGCN")
    rc = kmer.translate(trans)[::-1]
    return min(kmer, rc)

def get_current_memory_mb():
    """Devuelve el uso actual de memoria RAM en MB"""
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / (1024 * 1024)

def analyze_sample_performance(fastq_path, sketch_path, truth_csv_path):
    # Cargar Catálogo
    hll_R = load_hll_sketch(sketch_path)
    if hll_R is None: return None
    card_R = len(hll_R)

    hll_S = hyperloglog.HyperLogLog(HLL_ERROR_RATE)
    results = []
    
    # Variables Umbral
    calib_data = []
    threshold = 0.0
    calibrated = False
    
    # Variables de Ventana
    window_counter = 0
    window_id = 1
    
    print(f"   -> Procesando: {os.path.basename(fastq_path)}")
    
    # Iniciar timer para la primera ventana
    window_start_time = time.time()
    
    fq = pyfastx.Fastq(fastq_path, build_index=False)
    
    for name, seq, qual in fq:
        seq = str(seq).upper()
        
        # Procesar k-mers
        for i in range(len(seq) - K_SIZE + 1):
            kmer = seq[i:i+K_SIZE]
            if 'N' not in kmer:
                hll_S.add(get_canonical_kmer(kmer))
        
        window_counter += 1
        
        # --- FIN DE VENTANA ---
        if window_counter >= WINDOW_SIZE:
            # Medir Tiempo
            window_end_time = time.time()
            elapsed_seconds = window_end_time - window_start_time
            
            # Medir Memoria
            current_ram = get_current_memory_mb()
            
            # Cálculo HLL (Rareza)
            card_S = len(hll_S)
            hll_union = pickle.loads(pickle.dumps(hll_R))
            hll_union.update(hll_S)
            novelty = max(0, len(hll_union) - card_R)
            rho = novelty / card_S if card_S > 0 else 0.0
            
            # Umbral Dinámico
            if window_id <= CALIBRATION_WINDOWS:
                calib_data.append(rho)
            elif not calibrated:
                mu = np.mean(calib_data)
                sigma = np.std(calib_data)
                threshold = mu + (SIGMA_MULTIPLIER * sigma)
                calibrated = True
            
            # Guardar Resultados Completos
            results.append({
                "ventana": window_id,
                "rareza_estimada_hll": rho,
                "umbral_dinamico": threshold if calibrated else None,
                "segundos_procesamiento": elapsed_seconds, 
                "memoria_ram_mb": current_ram             
            })
            
            print(f"      V{window_id}: Rho={rho:.1%} | RAM={current_ram:.1f}MB | T={elapsed_seconds:.2f}s", end='\r')
            
            # Reset para siguiente ventana
            hll_S = hyperloglog.HyperLogLog(HLL_ERROR_RATE)
            window_counter = 0
            window_id += 1
            window_start_time = time.time() 
            
    # --- CRUZAR CON VERDAD ---
    df_est = pd.DataFrame(results)
    if os.path.exists(truth_csv_path):
        df_truth = pd.read_csv(truth_csv_path)
        df_final = pd.merge(df_truth, df_est, on="ventana", how="inner")
        
        # Rellenar umbral visualmente
        if calibrated:
             df_final['umbral_dinamico'] = df_final['umbral_dinamico'].fillna(threshold)
             
        return df_final
    return df_est

def main():
    if not os.path.exists(OUTPUT_RESULTS_DIR): os.makedirs(OUTPUT_RESULTS_DIR)
    
    # Detectar archivos
    files = [f for f in os.listdir(INPUT_MUESTRAS_DIR) if f.endswith(".fastq")]
    files.sort()
    
    if not files:
        print(f"No encontré archivos .fastq en '{INPUT_MUESTRAS_DIR}'. Revisa la ruta.")
        return

    for f in files:
        tipo = next((k for k in CATALOG_MAP if k in f), None)
        if tipo:
            truth_csv = os.path.join(INPUT_MUESTRAS_DIR, f.replace(".fastq", "_VERDAD.csv"))
            sketch = CATALOG_MAP[tipo]
            
            df = analyze_sample_performance(
                os.path.join(INPUT_MUESTRAS_DIR, f), 
                sketch, 
                truth_csv
            )
            
            if df is not None:
                out_name = f.replace(".fastq", "_PERFORMANCE.csv")
                df.to_csv(os.path.join(OUTPUT_RESULTS_DIR, out_name), index=False)
                print(f"\n      -> Guardado: {out_name}")

if __name__ == "__main__":
    main()