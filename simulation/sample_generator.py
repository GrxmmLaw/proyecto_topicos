
import os
import random
import time
import csv

# --- CONFIGURACIÓN PARA QUE SE VEAN LOS PICOS ---
TOTAL_READS = 5000000     
LONGITUD_READ = 150

# Configuración de Rareza Biológica (Genomas del catálogo)
PEAK_RARITY = 0.10        
N_DOMINANTES = 5          

# Configuración de Novedad Sintética (La que detectará el HLL)
PROB_NOVEDAD = 0.15         # Altura del PICO (15%)
DIVERGENCIA_NOVEDAD = 0.10  # Diferencia de la cepa (10%)

TASA_ERROR_TARGET = 0.0      

TAREAS = [
    {
        "input": "simulated_catalog/dataset_2GB_distintos",
        "fastq": "stream_distintos.fastq",
        "csv":   "meta_distintos.csv"
    },
    {
        "input": "simulated_catalog/dataset_2GB_medios_especies",
        "fastq": "stream_especies.fastq",
        "csv":   "meta_especies.csv"
    },
    {
        "input": "simulated_catalog/dataset_2GB_similares_cepas",
        "fastq": "stream_cepas.fastq",
        "csv":   "meta_cepas.csv"
    }
]

def cargar_genomas(carpeta):
    print(f"Cargando genomas de {carpeta}...")
    if not os.path.exists(carpeta):
        print(f"Error: No existe {carpeta}")
        return [], []
        
    archivos = sorted([os.path.join(carpeta, f) for f in os.listdir(carpeta) if f.endswith(".fna")])
    return genomas_memoria(archivos[:N_DOMINANTES]), genomas_memoria(archivos[N_DOMINANTES:])

def genomas_memoria(lista_archivos):
    seqs = []
    for ruta in lista_archivos:
        s = []
        with open(ruta, 'r') as f:
            for line in f:
                if not line.startswith('>'): s.append(line.strip())
        seqs.append("".join(s))
    return seqs

def calcular_probabilidad_rareza(progreso, pico):
    """Calcula una probabilidad en forma de ola triangular /\ """
    if progreso < 0.5:
        return pico * (progreso / 0.5)
    else:
        return pico * ((1.0 - progreso) / 0.5)

def mutar_genoma_sintetico(secuencia_base, tasa_divergencia):
    seq_list = list(secuencia_base)
    largo = len(seq_list)
    n_mutaciones = int(largo * tasa_divergencia)
    
    indices = random.sample(range(largo), n_mutaciones)
    bases = ['A', 'C', 'G', 'T']
    
    for i in indices:
        base_original = seq_list[i]
        posibles = [b for b in bases if b != base_original]
        seq_list[i] = random.choice(posibles)
        
    return "".join(seq_list)

def introducir_error_calibrado(secuencia, tasa):
    if tasa <= 0: return secuencia
    seq_list = list(secuencia)
    largo = len(seq_list)
    expected_errors = largo * tasa
    n_errores = int(expected_errors)
    
    if random.random() < (expected_errors - n_errores):
        n_errores += 1
    
    if n_errores > 0:
        indices = random.sample(range(largo), n_errores)
        bases = ['A', 'C', 'G', 'T', 'N']
        for i in indices:
            base_original = seq_list[i]
            val = [b for b in bases if b != base_original]
            seq_list[i] = random.choice(val)
    return "".join(seq_list)

def generar_stream(tarea):
    print(f"\n--- GENERANDO {tarea['fastq']} (5M reads) con NOVEDAD SINTÉTICA EN OLA ---")
    comunes, raros = cargar_genomas(tarea["input"])
    
    if not comunes: 
        print("⚠️ Genomas no encontrados.")
        return

    if comunes:
        print(f"   -> Creando cepa sintética mutante (Div: {DIVERGENCIA_NOVEDAD*100}%)...")
        cepa_novedosa_base = mutar_genoma_sintetico(comunes[0], DIVERGENCIA_NOVEDAD)
    else:
        cepa_novedosa_base = ""

    calidad = 'I' * LONGITUD_READ
    start = time.time()
    
    with open(tarea['fastq'], 'w') as fq, open(tarea['csv'], 'w', newline='') as log:
        w = csv.writer(log)
        w.writerow(["read_index", "prob_rareza_biologica", "prob_novedad_sintetica", "origen"])
        
        for i in range(TOTAL_READS):
            progreso = i / TOTAL_READS
            
            # --- MODIFICACIÓN CLAVE AQUÍ ---
            # Ahora la probabilidad de novedad TAMBIÉN sigue una ola, igual que la biológica
            prob_actual_novedad = calcular_probabilidad_rareza(progreso, PROB_NOVEDAD)
            prob_bio = calcular_probabilidad_rareza(progreso, PEAK_RARITY)

            rand_val = random.random()
 
            # 1. Intentar insertar Novedad Sintética (Variable)
            if rand_val < prob_actual_novedad and cepa_novedosa_base:
                seq = cepa_novedosa_base
                tipo = "NOVEDAD_SINTETICA"
    
            # 2. Intentar insertar Rareza Biológica (Variable)
            # Se suma prob_actual_novedad para no solapar los rangos de probabilidad
            elif rand_val < (prob_actual_novedad + prob_bio):
                if raros:
                    seq = random.choice(raros)
                    tipo = "RARO"
                else:
                    seq = random.choice(comunes)
                    tipo = "COMUN"
            
            # 3. Relleno Común
            else:
                seq = random.choice(comunes)
                tipo = "COMUN"
            
            # Recortar Read
            ini = random.randint(0, len(seq) - LONGITUD_READ)
            read_raw = seq[ini : ini + LONGITUD_READ]
            
            read_final = introducir_error_calibrado(read_raw, TASA_ERROR_TARGET)

            fq.write(f"@SIM_{i}:{tipo}\n{read_final}\n+\n{calidad}\n")
   
            if i % 5000 == 0:
                w.writerow([i, f"{prob_bio:.5f}", f"{prob_actual_novedad:.5f}", tipo])
            
            if i % 100000 == 0:
                print(f"   Progreso: {progreso*100:.1f}% | Novedad actual: {prob_actual_novedad:.3f}", end='\r')

    print(f"\n   -> Listo en {(time.time() - start)/60:.2f} min")

if __name__ == "__main__":
    for t in TAREAS: generar_stream(t)