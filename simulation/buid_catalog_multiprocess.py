import os
import pyfastx
import time
import hyperloglog
import pickle
import multiprocessing

K_SIZE = 31
HLL_ERROR_RATE = 0.01

TAREAS = [
    ("dataset_2GB_distintos",       "GT_distintos.sketch"),
    ("dataset_2GB_similares_cepas", "GT_similares_cepas.sketch"),
    ("dataset_2GB_medios_especies", "GT_medios_especies.sketch")
]


def get_canonical_kmer(kmer):
    kmer_upper = kmer.upper()
    trans = str.maketrans("ATCGN", "TAGCN")
    rc = kmer_upper.translate(trans)[::-1]
    return min(kmer_upper, rc)

def process_single_genome(args):

    file_path, k_size, error_rate = args

    local_hll = hyperloglog.HyperLogLog(error_rate)
    
    try:

        fa = pyfastx.Fasta(file_path, build_index=False)
        
        for name, seq in fa:
            seq_str = str(seq).upper()
            for i in range(len(seq_str) - k_size + 1):
                kmer = seq_str[i:i+k_size]
                if 'N' not in kmer:
                    local_hll.add(get_canonical_kmer(kmer))
                    
        return local_hll
        
    except Exception as e:
        print(f"\n[ERROR] Falló el archivo {os.path.basename(file_path)}: {e}")
        return None

def construir_sketch_para_carpeta(input_dir, output_sketch):
    """
    Función orquestadora para una carpeta específica.
    """
    print(f"\n" + "="*60)
    print(f"PROCESANDO CATÁLOGO: {input_dir}")
    print(f"="*60)

    if not os.path.exists(input_dir):
        print(f"Error: La carpeta '{input_dir}' no existe. ¿Ejecutaste el generador?")
        return


    files_to_process = []
    for root, dirs, files in os.walk(input_dir):
        for file in files:
            if file.endswith((".fna", ".fasta", ".fa")):
                full_path = os.path.join(root, file)
                files_to_process.append((full_path, K_SIZE, HLL_ERROR_RATE))
    
    total_files = len(files_to_process)
    print(f"-> Se encontraron {total_files} genomas en la carpeta.")
    
    if total_files == 0:
        return

    num_cores = multiprocessing.cpu_count()

    cores_to_use = max(1, num_cores - 1) 
    
    print(f"-> Iniciando Pool con {cores_to_use} núcleos...")
    start_time = time.time()

    master_hll = hyperloglog.HyperLogLog(HLL_ERROR_RATE)

  
    with multiprocessing.Pool(processes=cores_to_use) as pool:

        results = pool.imap_unordered(process_single_genome, files_to_process, chunksize=5)
        
        processed_count = 0
        for hll_result in results:
            processed_count += 1
            if hll_result is not None:
    
                master_hll.update(hll_result)
            

            if processed_count % 50 == 0 or processed_count == total_files:
                print(f"   Progreso: {processed_count}/{total_files} genomas...", end='\r')

    elapsed = time.time() - start_time
    print(f"\n-> Completado en {elapsed:.2f} segundos ({elapsed/total_files:.3f} s/genoma).")
    

    cardinalidad = len(master_hll)
    print(f"-> Cardinalidad Estimada (K-mers únicos): {cardinalidad:,}")
    
    print(f"-> Guardando sketch en: {output_sketch}")
    with open(output_sketch, "wb") as f:
        pickle.dump(master_hll, f)
    print("-> Guardado exitoso.\n")

def main():
    print(f"--- GENERADOR DE GROUND TRUTH (HLL) ---")
    print(f"K-mer size: {K_SIZE}")
    print(f"Error rate: {HLL_ERROR_RATE}")
    
    total_start = time.time()
    

    for carpeta, archivo_salida in TAREAS:
        construir_sketch_para_carpeta(carpeta, archivo_salida)
        
    print(f"--- TODO FINALIZADO en {(time.time() - total_start)/60:.2f} minutos ---")

if __name__ == "__main__":
    multiprocessing.freeze_support() 
    main()