import os
import pyfastx
import hyperloglog
import pickle
import csv
import copy
import statistics 
import time

K_SIZE = 31
HLL_ERROR_RATE = 0.01

WINDOW_SIZE = 200000       

CALIBRATION_WINDOWS = 5    
SIGMA_MULTIPLIER = 3.0     
MAX_SAFE_BASELINE = 0.30   

TAREAS = [
    {
        "id": "DISTINTOS (Ventana Grande)",
        "catalog": "HLL_Catalog/GT_distintssos.sketch",   
        "stream":  "stream_distintos.fastq",
        "output":  "resultado_distintos.csv"
    },
    {
        "id": "ESPECIES (Ventana Grande)",
        "catalog": "HLL_Catalog/GT_medios_especies.sketch",
        "stream":  "stream_especies.fastq",
        "output":  "resultado_especies.csv"
    },
    {
        "id": "CEPAS (Ventana Grande)",
        "catalog": "HLL_Catalog/GT_similares_cepas.sketch",
        "stream":  "stream_cepas.fastq",
        "output":  "resultado_cepas.csv"
    }
]

def get_canonical_kmer(kmer):
    trans = bytes.maketrans(b"ATCGN", b"TAGCN")
    kmer_bytes = kmer.encode('ascii')
    rc_kmer = kmer_bytes.translate(trans)[::-1]
    return min(kmer_bytes, rc_kmer).decode('ascii')

def get_kmers(sequence, k):
    seq = str(sequence).upper()
    largo = len(seq)
    if largo < k: return
    for i in range(largo - k + 1):
        kmer = seq[i:i+k]
        if 'N' not in kmer:
            yield get_canonical_kmer(kmer)

def procesar_stream(tarea):
    print(f"\n--- {tarea['id']} ---")
    if not os.path.exists(tarea["catalog"]) or not os.path.exists(tarea["stream"]):
        print("Faltan archivos: Revisa nombres o rutas.")
        return

    # Cargar HLL Referencia
    with open(tarea["catalog"], "rb") as f:
        hll_R = pickle.load(f)
    card_R = len(hll_R)
    print(f"   Ref Cardinalidad: {card_R:,}")
    
    hll_S = hyperloglog.HyperLogLog(HLL_ERROR_RATE)
    
    with open(tarea["output"], 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["ventana", "rareza_rho", "umbral_tau", "estado", "alerta", "card_S", "card_U", "card_R"])

        calibration_data = []
        dynamic_threshold = 0.0
        is_calibrated = False
        reads_in_window = 0
        window_count = 0
        
        # Leemos el FASTQ
        fq = pyfastx.Fastq(tarea["stream"])
        for read in fq:
            for kmer in get_kmers(read.seq, K_SIZE):
                hll_S.add(kmer)
            reads_in_window += 1

            if reads_in_window >= WINDOW_SIZE:
                window_count += 1
                card_S = len(hll_S)
                
                # Unión y Cálculo
                hll_Union = copy.deepcopy(hll_R)
                hll_Union.update(hll_S)
                card_U = len(hll_Union)

                # Cálculo de Rho
                rho = 0.0
                if card_S > 0:
                    rho = (card_U - card_R) / card_S

                # Lógica de Detección
                status = "CALIBRANDO"
                alerta = 0
                
                if not is_calibrated:
                    calibration_data.append(rho)
                    if len(calibration_data) >= CALIBRATION_WINDOWS:
                        mean_val = statistics.mean(calibration_data)
                        stdev_val = statistics.stdev(calibration_data) if len(calibration_data) > 1 else 0.0
                        dynamic_threshold = mean_val + (SIGMA_MULTIPLIER * stdev_val)
                        if dynamic_threshold > MAX_SAFE_BASELINE: dynamic_threshold = MAX_SAFE_BASELINE
                        is_calibrated = True
                        print(f"   [Calibrado] Base esperada: {mean_val:.3f} | Umbral: {dynamic_threshold:.3f}")
                else:
                    status = "MONITOREO"
                    if rho > dynamic_threshold: alerta = 1
                
                writer.writerow([window_count, f"{rho:.6f}", f"{dynamic_threshold:.6f}", status, alerta, card_S, card_U, card_R])
                
                # Reiniciar ventana
                reads_in_window = 0
                hll_S = hyperloglog.HyperLogLog(HLL_ERROR_RATE)
                
                print(f"   V{window_count}: ρ={rho:.4f} {'🚨' if alerta else ''}", end='\r')

    print(f"\nGuardado: {tarea['output']}")

if __name__ == "__main__":
    for t in TAREAS: procesar_stream(t)