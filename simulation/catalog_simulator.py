import os
import random
import time

# --- Configuracion inicial ---
TARGET_SIZE_GB = 1.0         # Tamaño objetivo por cada carpeta
LONGITUD_BASE = 4000000      # 4 millones de pares de bases por genoma 
GC_CONTENT = 0.50            # Contenido de Guanina-Citosina

# Calculo automatico de cantidad de archivos necesarios
BYTES_TOTALES = TARGET_SIZE_GB * (1024**3)
N_GENOMAS = int(BYTES_TOTALES / LONGITUD_BASE)

BASE_CHARS = ['A', 'C', 'G', 'T']
BASE_DIR = "sim_catalogs" 

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
    print(f"=== Generador de catlogos sinteticos ===")
    print(f"Carpeta base: {BASE_DIR}/")
    print(f"Configuración: {TARGET_SIZE_GB} GB por set.")
    print(f"Genomas a generar: {N_GENOMAS} por set.")
    print("-" * 50)

    dirs = {
        "distintos": os.path.join(BASE_DIR, "datasets_distintos"),
        "similares": os.path.join(BASE_DIR, "datasets_similares"),
        "medios":    os.path.join(BASE_DIR, "datasets_medios")
    }

    for d in dirs.values():
        if not os.path.exists(d): 
            os.makedirs(d)
            print(f"Creada carpeta: {d}")

    t_inicio_total = time.time()

    # Genomas distintos
    print(f"\n[1/3] Generando {N_GENOMAS} genomas distintos")
    t_fase = time.time()
    for i in range(1, N_GENOMAS + 1):
        seq = generar_secuencia_rapida(LONGITUD_BASE, GC_CONTENT)
        guardar_fna(dirs["distintos"], f"genome_diff_{i:04d}", f"synth_diff_{i}", seq)
        imprimir_progreso(i, N_GENOMAS, t_fase)


    # Genomas similares
    print(f"\n\n[2/3] Generando {N_GENOMAS} genomas SIMILARES (Mutación 0.1%)...")
    t_fase = time.time()
    
    # Nueva semilla
    semilla = generar_secuencia_rapida(LONGITUD_BASE, GC_CONTENT)
    guardar_fna(dirs["similares"], "genome_strain_0000_SEED", "seed_strain", semilla)

    for i in range(1, N_GENOMAS): 
        # Mutacion del 0.1%
        seq_mut = mutar_secuencia(semilla, 0.001) 
        guardar_fna(dirs["similares"], f"genome_strain_{i:04d}", f"synth_strain_{i}", seq_mut)
        imprimir_progreso(i, N_GENOMAS, t_fase)


    # Genomas medianamente similares
    print(f"\n\n[3/3] Generando {N_GENOMAS} genomas MEDIOS (Mutación 3.0%)...")
    t_fase = time.time()

    # Nueva semilla
    semilla2 = generar_secuencia_rapida(LONGITUD_BASE, GC_CONTENT)
    guardar_fna(dirs["medios"], "genome_spp_0000_SEED", "seed_spp", semilla2)

    for i in range(1, N_GENOMAS):
        # Mutacion del 3% 
        seq_mut = mutar_secuencia(semilla2, 0.03)
        guardar_fna(dirs["medios"], f"genome_spp_{i:04d}", f"synth_spp_{i}", seq_mut)
        imprimir_progreso(i, N_GENOMAS, t_fase)

    duracion = (time.time() - t_inicio_total) / 60
    print(f"\n\n¡Proceso finalizado! Tiempo total: {duracion:.2f} minutos.")
    print(f"Los archivos están en '{os.path.abspath(BASE_DIR)}'")

if __name__ == "__main__":
    main()