import os
import pyfastx
import time
import hyperloglog
import pickle
import multiprocessing
from functools import partial

K_SIZE = 31
HLL_ERROR_RATE = 0.01
INPUT_DIR = "UHGG_reps"   
OUTPUT_SKETCH = "HLL_R_UHGG.sketch"

def get_canonical_kmer(kmer):
    kmer_upper = kmer.upper()
    trans = str.maketrans("ATCGN", "TAGCN")
    rc = kmer_upper.translate(trans)[::-1]
    return min(kmer_upper, rc)

def process_single_genome(file_info):
    """
    Esta función es ejecutada por cada núcleo del procesador de forma independiente.
    Recibe: (ruta_archivo, k_size, error_rate)
    Devuelve: Un objeto HLL con los k-mers de ese archivo.
    """
    file_path, k, error = file_info
    local_hll = hyperloglog.HyperLogLog(error)
    kmers_count = 0
    
    try:
        fa = pyfastx.Fasta(file_path, build_index=False)
        for name, seq in fa:
            seq_str = str(seq).upper()
            for i in range(len(seq_str) - k + 1):
                kmer = seq_str[i:i+k]
                if 'N' not in kmer:
                    local_hll.add(get_canonical_kmer(kmer))
                    kmers_count += 1
        return local_hll
    except Exception as e:
        print(f"Error procesando {file_path}: {e}")
        return None

def main():
    print(f"--- Constructor de Catálogo Paralelo (UHGG) ---")
    print(f"Buscando archivos .fna en '{INPUT_DIR}'...")
    
    files_to_process = []
    for root, dirs, files in os.walk(INPUT_DIR):
        for file in files:
            if file.endswith(".fna") or file.endswith(".fasta") or file.endswith(".fa"):
                full_path = os.path.join(root, file)
                files_to_process.append((full_path, K_SIZE, HLL_ERROR_RATE))
    
    total_files = len(files_to_process)
    print(f"Encontrados {total_files} genomas.")
    
    if total_files == 0:
        print("No hay archivos para procesar.")
        return

    num_cores = multiprocessing.cpu_count()
    print(f"Iniciando procesamiento con {num_cores} núcleos paralelos...")
    
    start_time = time.time()

    master_hll = hyperloglog.HyperLogLog(HLL_ERROR_RATE)
    

    with multiprocessing.Pool(processes=num_cores) as pool:
        results = pool.imap_unordered(process_single_genome, files_to_process, chunksize=10)
        
        processed_count = 0
        for hll_result in results:
            processed_count += 1
            if hll_result is not None:
                master_hll.update(hll_result)
            
            if processed_count % 100 == 0:
                print(f"  Progres: {processed_count}/{total_files} genomas integrados...")

    elapsed = time.time() - start_time
    print(f"\n--- Procesamiento completado en {elapsed:.2f} segundos ---")
    print(f"Promedio: {elapsed/total_files:.3f} s/genoma")
    print(f"Cardinalidad Final del Catálogo: {len(master_hll)} k-mers únicos")
    print(f"Guardando {OUTPUT_SKETCH}...")
    with open(OUTPUT_SKETCH, "wb") as f:
        pickle.dump(master_hll, f)
    print("¡Listo!")

if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()