import os
import argparse  # <--- NUEVO: Para recibir argumentos
import pyfastx
import time
import hyperloglog
import pickle
import csv
import copy
import statistics 

# --- PARÁMETROS FIJOS ---
K_SIZE = 31
HLL_ERROR_RATE = 0.01
WINDOW_SIZE = 10000 
# Umbral Dinámico
CALIBRATION_WINDOWS = 5 
SIGMA_MULTIPLIER = 3.0 
MAX_SAFE_BASELINE = 0.20  

# --- FUNCIONES AUXILIARES ---
def get_canonical_kmer(kmer):
    # Versión optimizada sin llamar a función externa para velocidad
    trans = str.maketrans("ATCGN", "TAGCN")
    rc_kmer = kmer.upper().translate(trans)[::-1]
    return min(kmer, rc_kmer)

def get_kmers(sequence, k):
    seq = str(sequence).upper()
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if 'N' not in kmer:
            yield get_canonical_kmer(kmer)

def read_generator(file_path):
    # Pyfastx maneja gzip automáticamente
    fq = pyfastx.Fastq(file_path)
    for read in fq:
        yield read

# --- LÓGICA PRINCIPAL ---
def main():
    # 1. CONFIGURACIÓN DE ARGUMENTOS (CLI)
    parser = argparse.ArgumentParser(description="Analizador de Rareza para Benchmark HLL")
    parser.add_argument("-c", "--catalog", required=True, help="Ruta al archivo .sketch (Catálogo)")
    parser.add_argument("-i", "--input", required=True, help="Ruta al archivo .fq o .fq.gz (Muestra)")
    parser.add_argument("-o", "--output", required=True, help="Nombre del archivo CSV de salida")
    
    args = parser.parse_args()

    print(f"--- INICIANDO BENCHMARK ---")
    print(f"Catálogo: {args.catalog}")
    print(f"Muestra:  {args.input}")
    print(f"Salida:   {args.output}")

    # 2. CARGAR CATÁLOGO
    if not os.path.exists(args.catalog):
        print(f"❌ Error: No encuentro el catálogo {args.catalog}")
        return

    print(f"Cargando Sketch de Referencia...")
    with open(args.catalog, "rb") as f:
        hll_R = pickle.load(f)
    
    card_R = len(hll_R)
    print(f"✅ Catálogo cargado. Cardinalidad Base |R| = {card_R}")

    # 3. PREPARAR CSV CON MÁS DATOS
    # Agregamos columnas 'card_S' y 'card_U' para ver el comportamiento del HLL
    with open(args.output, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        headers = ["ventana", "rareza_rho", "umbral_tau", "estado", "alerta", "card_S", "card_U", "card_R"]
        csv_writer.writerow(headers)

        # Variables de estado
        calibration_data = []
        dynamic_threshold = 0.0
        is_calibrated = False
        window_count = 0
        reads_in_window = 0
        
        # HLL Temporal para la ventana
        hll_S = hyperloglog.HyperLogLog(HLL_ERROR_RATE)
        
        try:
            generator = read_generator(args.input)
            
            for read in generator:
                # Agregar k-mers al HLL de la ventana
                for kmer in get_kmers(read.seq, K_SIZE):
                    hll_S.add(kmer)
                
                reads_in_window += 1

                # --- FIN DE VENTANA ---
                if reads_in_window >= WINDOW_SIZE:
                    window_count += 1
                    
                    # A. Cálculos de Cardinalidad
                    card_S = len(hll_S)
                    
                    # Clonamos R para hacer la unión (operación costosa pero necesaria)
                    hll_Union = copy.deepcopy(hll_R)
                    hll_Union.update(hll_S)
                    card_U = len(hll_Union)

                    # B. Cálculo de Rareza (Protegido contra división por cero)
                    rho = 0.0
                    if card_S > 0:
                        rho = (card_U - card_R) / card_S
                    
                    # C. Lógica de Umbral Dinámico
                    status = "CALIBRANDO"
                    alerta_flag = "NO"
                    
                    if not is_calibrated:
                        calibration_data.append(rho)
                        if len(calibration_data) >= CALIBRATION_WINDOWS:
                            mean_val = statistics.mean(calibration_data)
                            stdev_val = statistics.stdev(calibration_data) if len(calibration_data) > 1 else 0.0
                            
                            # Safety Check
                            if mean_val > MAX_SAFE_BASELINE:
                                dynamic_threshold = MAX_SAFE_BASELINE
                            else:
                                dynamic_threshold = mean_val + (SIGMA_MULTIPLIER * stdev_val)
                                # Evitar umbrales demasiado pegados al promedio
                                if dynamic_threshold < mean_val + 0.01: dynamic_threshold = mean_val + 0.01
                            
                            is_calibrated = True
                    else:
                        status = "MONITOREO"
                        if rho > dynamic_threshold:
                            alerta_flag = "SI"
                            print(f"Ventana {window_count}: 🚨 ALERTA (ρ={rho:.4f})")
                        else:
                            if window_count % 5 == 0: # Imprimir cada 5 para no saturar
                                print(f"Ventana {window_count}: OK (ρ={rho:.4f})")

                    # D. Guardar Datos
                    csv_writer.writerow([
                        window_count, 
                        f"{rho:.6f}", 
                        f"{dynamic_threshold:.6f}", 
                        status, 
                        alerta_flag,
                        card_S,   # <--- Dato clave para tu benchmark
                        card_U,
                        card_R
                    ])
                    
                    # E. Reiniciar Ventana
                    reads_in_window = 0
                    hll_S = hyperloglog.HyperLogLog(HLL_ERROR_RATE)

        except Exception as e:
            print(f"❌ Error procesando archivo: {e}")

    print(f"\n✅ Análisis completado. Resultados en: {args.output}")

if __name__ == "__main__":
    main()