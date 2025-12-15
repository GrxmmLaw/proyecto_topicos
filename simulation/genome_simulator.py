import os
import random
import time

TARGET_SIZE_GB = 2.0  

LONGITUD_BASE = 4000000  

GC_CONTENT = 0.50        

BYTES_TOTALES = TARGET_SIZE_GB * (1024**3)
N_GENOMAS = int(BYTES_TOTALES / LONGITUD_BASE)

BASE_CHARS = ['A', 'C', 'G', 'T']


def generar_secuencia_rapida(longitud, gc_target):
    prob_gc = gc_target / 2.0
    prob_at = (1.0 - gc_target) / 2.0
    pesos = [prob_at, prob_gc, prob_gc, prob_at]
    bases = random.choices(BASE_CHARS, weights=pesos, k=longitud)
    return "".join(bases)

def mutar_secuencia(secuencia, tasa_mutacion):
    seq_list = list(secuencia)
    num_mutaciones = int(len(seq_list) * tasa_mutacion)
    indices = random.sample(range(len(seq_list)), num_mutaciones)
    
    for idx in indices:
        base_original = seq_list[idx]
        posibles = [b for b in BASE_CHARS if b != base_original]
        seq_list[idx] = random.choice(posibles)
        
    return "".join(seq_list)

def guardar_fna(carpeta, nombre, cabecera, secuencia):
    ruta = os.path.join(carpeta, nombre + ".fna")
    with open(ruta, 'w') as f:
        f.write(f">{cabecera}\n")
        for i in range(0, len(secuencia), 80):
            f.write(secuencia[i:i+80] + "\n")

def imprimir_progreso(actual, total, inicio):
    porcentaje = (actual / total) * 100
    tiempo_pasado = time.time() - inicio
    if actual > 0:
        tiempo_por_item = tiempo_pasado / actual
        restante = (total - actual) * tiempo_por_item
        print(f"\rProgreso: [{actual}/{total}] {porcentaje:.1f}% - Restante: {restante/60:.1f} min  ", end="")

def main():
    print(f"=== Generador de genomas sintéticos ===")
    print(f"Configuración: {TARGET_SIZE_GB} GB por catálogo.")
    print(f"Tamaño Genoma: {LONGITUD_BASE/1e6} MB.")
    print(f"Archivos a generar por carpeta: {N_GENOMAS}")
    print(f"Total datos a escribir en disco: {TARGET_SIZE_GB * 3} GB")
    print("-" * 50)

    dirs = {
        "distintos": "dataset_distintos",
        "similares": "dataset_similares_cepas",
        "medios":    "dataset_medios_especies"
    }

    for d in dirs.values():
        if not os.path.exists(d): os.makedirs(d)

    t_inicio_total = time.time()

    print(f"\n[1/3] Generando {N_GENOMAS} genomas DISTINTOS...")
    t_fase = time.time()
    for i in range(1, N_GENOMAS + 1):
        seq = generar_secuencia_rapida(LONGITUD_BASE, GC_CONTENT)
        guardar_fna(dirs["distintos"], f"genome_diff_{i:04d}", f"synth_diff_{i}", seq)
        imprimir_progreso(i, N_GENOMAS, t_fase)


    print(f"\n\n[2/3] Generando {N_GENOMAS} genomas SIMILARES (Cepas)...")
    t_fase = time.time()
    
    semilla = generar_secuencia_rapida(LONGITUD_BASE, GC_CONTENT)
    guardar_fna(dirs["similares"], "genome_strain_0000_SEED", "seed_strain", semilla)

    for i in range(1, N_GENOMAS): 
        seq_mut = mutar_secuencia(semilla, 0.001) 
        guardar_fna(dirs["similares"], f"genome_strain_{i:04d}", f"synth_strain_{i}", seq_mut)
        imprimir_progreso(i, N_GENOMAS, t_fase)

    print(f"\n\n[3/3] Generando {N_GENOMAS} genomas MEDIOS (Especies)...")
    t_fase = time.time()

    semilla2 = generar_secuencia_rapida(LONGITUD_BASE, GC_CONTENT)
    guardar_fna(dirs["medios"], "genome_spp_0000_SEED", "seed_spp", semilla2)

    for i in range(1, N_GENOMAS):
        seq_mut = mutar_secuencia(semilla2, 0.03)
        guardar_fna(dirs["medios"], f"genome_spp_{i:04d}", f"synth_spp_{i}", seq_mut)
        imprimir_progreso(i, N_GENOMAS, t_fase)

    duracion = (time.time() - t_inicio_total) / 60
    print(f"Tiempo total: {duracion:.2f} minutos.")

if __name__ == "__main__":
    main()