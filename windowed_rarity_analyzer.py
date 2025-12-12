import os
import pyfastx
import time
import hyperloglog
import pickle
import csv
import copy
import statistics 

# --- 1. Parámetros de Configuración ---
K_SIZE = 31
SKETCH_FILE_IN = "HLL_R.sketch"
HLL_ERROR_RATE = 0.01
MAX_SAFE_BASELINE = 0.20  
SAMPLE_S_FASTQ = "S_zymo.fq.gz" 

WINDOW_SIZE = 200000 
CSV_OUTPUT_FILE = "rarity_log_dinamico.csv"

# --- CONFIGURACIÓN DEL UMBRAL DINÁMICO ---
CALIBRATION_WINDOWS = 5 
SIGMA_MULTIPLIER = 3.0  

# Banderas de Logging
LOG_PROGRESO_VENTANA = True
INCREMENTO_PROGRESO = 50000 

# --- Funciones ---
def reverse_complement(kmer):
    trans = str.maketrans("ATCGN", "TAGCN")
    return kmer.upper().translate(trans)[::-1]

def get_canonical_kmer(kmer):
    rc_kmer = reverse_complement(kmer)
    return min(kmer, rc_kmer)

def get_kmers(sequence, k):
    seq = str(sequence).upper()
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if 'N' not in kmer:
            yield get_canonical_kmer(kmer)

def read_generator(file_path):
    fq = pyfastx.Fastq(file_path)
    for read in fq:
        yield read

# --- Lógica Principal ---
def main():
    print(f"--- Analizador de Rareza con UMBRAL DINÁMICO ---")
    print(f"Ventana: {WINDOW_SIZE} lecturas")
    print(f"Calibración: Primeras {CALIBRATION_WINDOWS} ventanas")
    
    # Cargar Catálogo
    if not os.path.exists(SKETCH_FILE_IN):
        print(f"Error: No existe el catálogo {SKETCH_FILE_IN}")
        return
    
    print(f"Cargando Catálogo R...")
    with open(SKETCH_FILE_IN, "rb") as f:
        hll_R = pickle.load(f)
    print(f"Catálogo cargado. |R| = {len(hll_R)}")

    # Preparar Muestra
    try:
        lecturas_gen = read_generator(SAMPLE_S_FASTQ)
        print(f"Abriendo flujo de datos: {SAMPLE_S_FASTQ}")
    except Exception as e:
        print(f"Error: {e}")
        return

    # Variables para el Umbral Dinámico
    calibration_data = []
    dynamic_threshold = None
    is_calibrated = False

    with open(CSV_OUTPUT_FILE, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        # Agregamos columna de estado para saber si estaba calibrando o detectando
        csv_writer.writerow(["ventana", "rareza_rho", "umbral_tau", "estado", "alerta"])
        
        window_count = 0
        finished = False
        
        while not finished:
            window_count += 1
            print(f"\n--- Procesando Ventana {window_count} ---")
            
            # --- Lógica de Ventana (Tu código estándar) ---
            hll_S = hyperloglog.HyperLogLog(HLL_ERROR_RATE)
            reads = 0
            t_start = time.time()
            
            while reads < WINDOW_SIZE:
                try:
                    read = next(lecturas_gen)
                    reads += 1
                    for kmer in get_kmers(read.seq, K_SIZE):
                        hll_S.add(kmer)
                except StopIteration:
                    finished = True
                    break
            
            if reads == 0: break

            # --- Cálculos ---
            card_S = len(hll_S)
            hll_U = copy.deepcopy(hll_R)
            hll_U.update(hll_S)
            
            rho = 0.0
            if card_S > 0:
                rho = (len(hll_U) - len(hll_R)) / card_S
            
            print(f"  Rareza actual (ρ): {rho:.6f} ({rho*100:.4f}%)")

            # --- Umbral dinámico  ---
            status = "CALIBRANDO"
            alerta = False
            current_threshold_display = 0.0

            if not is_calibrated:
                # Fase 1: Acumular datos
                calibration_data.append(rho)
                print(f"  [Calibración] Recopilando datos ({len(calibration_data)}/{CALIBRATION_WINDOWS})...")
                
                # Si se termina de recopilar, se calcula el umbral
                if len(calibration_data) >= CALIBRATION_WINDOWS:
                    print("  --- FIN DE CALIBRACIÓN ---")
                    mean_val = statistics.mean(calibration_data)
                    
                    if len(calibration_data) > 1:
                        stdev_val = statistics.stdev(calibration_data)
                    else:
                        stdev_val = 0.0
                    
                    # --- CORRECCIÓN DE SEGURIDAD (SANITY CHECK) ---
                    if mean_val > MAX_SAFE_BASELINE:
                        print(f"  ⚠️ ADVERTENCIA: La rareza inicial ({mean_val*100:.2f}%) es anormalmente alta.")
                        print(f"  ⚠️ Forzando umbral de seguridad fijo.")
                        
                        dynamic_threshold = MAX_SAFE_BASELINE

                    else:
                        dynamic_threshold = mean_val + (SIGMA_MULTIPLIER * stdev_val)
                        if dynamic_threshold < mean_val + 0.01: 
                             dynamic_threshold = mean_val + 0.01

                    print(f"  >> Ruido Base detectado: {mean_val*100:.4f}%")
                    print(f"  >> NUEVO UMBRAL FIJADO (Tau): {dynamic_threshold*100:.4f}%")
                    
                    is_calibrated = True
            else:
                status = "MONITOREO"
                current_threshold_display = dynamic_threshold
                
                if rho > dynamic_threshold:
                    print(f"  🚨 ¡¡¡ALERTA DE ANOMALÍA!!! ({rho*100:.4f}% > {dynamic_threshold*100:.4f}%)")
                    alerta = True
                else:
                    print(f"  [OK] Dentro del rango normal (< {dynamic_threshold*100:.4f}%)")

            csv_writer.writerow([window_count, rho, current_threshold_display, status, "SI" if alerta else "NO"])

    print("\nAnálisis finalizado.")

if __name__ == "__main__":
    main()