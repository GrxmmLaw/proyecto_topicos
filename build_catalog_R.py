import os
import pyfastx
import time
import hyperloglog
import pickle

# --- Configuración ---
K_SIZE = 31
GENOMES_DIR = "genomas_test/Experimento1/ZymoBIOMICS.STD.refseq.v2/Genomes"
SKETCH_FILE_OUT = "HLL_R_T2.sketch"
HLL_ERROR_RATE = 0.1
LOG_CADA_ARCHIVO = True 

# --- Funciones Canónicas ---
def reverse_complement(kmer):
    trans = str.maketrans("ATCGN", "TAGCN")
    return kmer.upper().translate(trans)[::-1]

def get_canonical_kmer(kmer):
    rc_kmer = reverse_complement(kmer)
    return min(kmer, rc_kmer)

# --- Logica de cálculo ---
def get_kmers(sequence, k):
    """Extrae k-mers canónicos. Ignora cualquier k-mer con 'N'."""
    seq = str(sequence).upper()
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if 'N' not in kmer: 
            yield get_canonical_kmer(kmer)

# --- Main ---
def main():
    print(f"--- Iniciando construcción del catálogo R (k={K_SIZE}) ---")
    print(f"Directorio: {os.path.abspath(GENOMES_DIR)}")
    print(f"Salida: {SKETCH_FILE_OUT}")
    
    start_time = time.time()
    hll_R = hyperloglog.HyperLogLog(HLL_ERROR_RATE)

    total_files = 0
    total_kmers = 0

    for root, dirs, files in os.walk(GENOMES_DIR):
        if "pan-genome" in dirs: dirs.remove("pan-genome")
        
        for file in files:
            if file.endswith(".fna") or file.endswith(".fna.gz") or file.endswith(".fasta"):
                total_files += 1
                path = os.path.join(root, file)
                kmers_file = 0
                
                try:
                    fa = pyfastx.Fasta(path)
                    for seq in fa:
                        for kmer in get_kmers(seq.seq, K_SIZE):
                            hll_R.add(kmer)
                            kmers_file += 1
                    
                    total_kmers += kmers_file
                    if LOG_CADA_ARCHIVO:
                        print(f"  [Archivo {total_files:04d}] K-mers: {kmers_file: <10} | {file}")

                except Exception as e:
                    print(f"Error en {file}: {e}")

    print("-" * 30)
    print(f"Total Archivos: {total_files}")
    print(f"Cardinalidad Final: {len(hll_R)}")
    
    with open(SKETCH_FILE_OUT, "wb") as f:
        pickle.dump(hll_R, f)
    
    print(f"Tiempo total: {time.time() - start_time:.2f}s")

if __name__ == "__main__":
    main()