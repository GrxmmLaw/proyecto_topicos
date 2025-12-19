import os
import random
import time
import pandas as pd
import math

# Configuracion inicial
OUTPUT_DIR = "muestras_validacion"
WINDOW_SIZE = 200000    
TOTAL_READS = 6000000   
READ_LENGTH = 150
K_SIZE = 31

# Configuracion de Onda de Ruido 
NOISE_CYCLES = 4       
MAX_NOISE_RARITY = 0.05 

# Configuracion de Evento
EVENT_START = 2000000
EVENT_END = 4000000
PEAK_PROB = 0.10        

RUTAS_GENOMAS = {
    "DISTINTOS": {
        "fondo": "sim_catalogs/datasets_distintos/genome_diff_0001.fna",
        "anomalo": "sim_catalogs/datasets_medios/genome_spp_0001.fna"
    },
    "MEDIOS": {
        "fondo": "sim_catalogs/datasets_medios/genome_spp_0001.fna",
        "anomalo": "sim_catalogs/datasets_distintos/genome_diff_0001.fna"
    },
    "SIMILARES": {
        "fondo": "sim_catalogs/datasets_similares/genome_strain_0001.fna",
        "anomalo": "sim_catalogs/datasets_distintos/genome_diff_0001.fna"
    }
}

def cargar_secuencia_limpia(ruta):
    if not os.path.exists(ruta):
        print(f"[ERROR] No se encuentra: {ruta}")
        return "N" * READ_LENGTH 
    with open(ruta, 'r') as f:
        lines = f.readlines()
        seq = "".join([l.strip().upper() for l in lines if not l.startswith(">")])
    return seq

def simular_read_con_error(seq_fuente, error_rate):
    if len(seq_fuente) < READ_LENGTH: return "N" * READ_LENGTH
    
    if error_rate < 1e-6:
        start = random.randint(0, len(seq_fuente) - READ_LENGTH)
        return seq_fuente[start : start + READ_LENGTH]

    start = random.randint(0, len(seq_fuente) - READ_LENGTH)
    read_list = list(seq_fuente[start : start + READ_LENGTH])
    bases = ['A', 'C', 'G', 'T']
    
    n_mutaciones = int(READ_LENGTH * error_rate)
    
    if n_mutaciones == 0 and random.random() < (READ_LENGTH * error_rate):
        n_mutaciones = 1

    if n_mutaciones > 0:
        indices = random.sample(range(READ_LENGTH), n_mutaciones)
        for idx in indices:
            read_list[idx] = random.choice(bases)
            
    return "".join(read_list)

def get_wave_error_rate(read_index, total_reads):

    phase = (read_index / total_reads) * (2 * math.pi * NOISE_CYCLES)

    wave_position = (math.sin(phase) + 1) / 2
 
    target_rarity = wave_position * MAX_NOISE_RARITY

    if target_rarity <= 0: return 0.0, 0.0
    
    base_error = 1 - (1 - target_rarity)**(1/K_SIZE)
    
    return base_error, target_rarity

def get_prob_triangular(idx, start, end, peak):
    if not (start <= idx < end): return 0.0
    mid = (start + end) / 2
    if idx < mid:
        return peak * (idx - start) / (mid - start)
    else:
        return peak * (end - idx) / (end - mid)

def generar_muestra(tipo_catalogo, modo, rutas):
    print(f"\n--- Generando {tipo_catalogo} [{modo}]  ---")
    
    if not os.path.exists(OUTPUT_DIR): os.makedirs(OUTPUT_DIR)
    nombre_base = f"{tipo_catalogo}_{modo}"
    fastq_path = os.path.join(OUTPUT_DIR, f"{nombre_base}.fastq")
    csv_path = os.path.join(OUTPUT_DIR, f"{nombre_base}_VERDAD.csv")

    seq_fondo = cargar_secuencia_limpia(rutas["fondo"])
    seq_anomalo = ""
    if modo == "TEST":
        seq_anomalo = cargar_secuencia_limpia(rutas["anomalo"])

    window_stats = []
    sum_valor_grafico = 0 
    count_window_reads = 0
    
    with open(fastq_path, 'w') as fq:
        for i in range(TOTAL_READS):
            
            if modo == "CONTROL":
                tasa_error_actual, rareza_esperada = get_wave_error_rate(i, TOTAL_READS)
                usar_anomalia = False
                valor_grafico_verdad = rareza_esperada 
                
            else:  
                tasa_error_actual = 0.001 
                prob_anomalia = get_prob_triangular(i, EVENT_START, EVENT_END, PEAK_PROB)
                usar_anomalia = (random.random() < prob_anomalia)
                valor_grafico_verdad = 1.0 if usar_anomalia else 0.0

            if usar_anomalia:
                read = simular_read_con_error(seq_anomalo, tasa_error_actual)
                label = "ANOMALIA"
            else:
                read = simular_read_con_error(seq_fondo, tasa_error_actual)
                label = "FONDO"

            fq.write(f"@READ_{i}_{label}\n{read}\n+\n{'I'*READ_LENGTH}\n")
            
            sum_valor_grafico += valor_grafico_verdad
            count_window_reads += 1

            if (i + 1) % WINDOW_SIZE == 0:
                ventana_idx = (i + 1) // WINDOW_SIZE
                promedio_verdad = sum_valor_grafico / count_window_reads
                
                window_stats.append({
                    "ventana": ventana_idx,
                    "rareza_real": promedio_verdad,
                    "modo": modo
                })
                
                sum_valor_grafico = 0
                count_window_reads = 0
                
                if ventana_idx % 5 == 0:
                    print(f"   Generada V{ventana_idx}...", end='\r')

    pd.DataFrame(window_stats).to_csv(csv_path, index=False)
    print(f"\n   -> Guardado: {fastq_path}")

def main():
    start = time.time()
    for catalogo, rutas in RUTAS_GENOMAS.items():
        generar_muestra(catalogo, "CONTROL", rutas)
        generar_muestra(catalogo, "TEST", rutas)
    print(f"\nTIEMPO TOTAL: {(time.time()-start)/60:.2f} min.")

if __name__ == "__main__":
    main()