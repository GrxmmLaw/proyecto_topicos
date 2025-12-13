import os
import random
import gzip
import shutil

# --- CONFIGURACIÓN ---
OUTPUT_BASE = "entornos_reales_v2"
LONGITUD_LECTURA = 150
COBERTURA = 10  # Profundidad de la muestra

# Definimos los ambientes con los nombres CIENTÍFICOS exactos para NCBI
AMBIENTES = {
    # 1. ENTORNO INFECCIÓN (Variantes de S. aureus)
    # Estrategia: Buscamos cepas famosas por su nombre corto
    "1_INFECCION_HOSPITAL": [
        "Staphylococcus aureus subsp. aureus NCTC 8325",
        "Staphylococcus aureus subsp. aureus USA300",  # Nombre acortado
        "Staphylococcus aureus subsp. aureus COL",
        "Staphylococcus aureus subsp. aureus Mu50",    # Reemplazo de MRSA252 (que falló)
        "Staphylococcus aureus subsp. aureus MW2"
    ],
    
    # 2. ENTORNO INTESTINO (Diversidad Media)
    # Estrategia: Nombres más genéricos para asegurar hits.
    # Reemplazamos Prevotella (que suele ser Draft) por B. fragilis (Complete)
    "2_INTESTINO_HUMANO": [
        "Escherichia coli str. K-12",             # Quitamos "MG1655" para asegurar hit
        "Bacteroides thetaiotaomicron",
        "Faecalibacterium prausnitzii",
        "Bifidobacterium longum",                 # Quitamos "subsp." y cepa
        "Bacteroides fragilis"                    # Reemplazo seguro para Prevotella
    ],
    
    # 3. ENTORNO MARINO (Alta Diversidad)
    "3_OCEANO_MAR": [
        "Candidatus Pelagibacter ubique",
        "Prochlorococcus marinus",                # Quitamos cepa específica
        "Synechococcus elongatus",                # WH 8102 a veces falla, elongatus es seguro
        "Vibrio cholerae",
        "Alteromonas macleodii"
    ]
}

def obtener_url_descarga(nombre_bacteria, datos_ncbi):
    """Busca la bacteria en el índice de NCBI y devuelve la URL"""
    # 1. Intento exacto
    for entry in datos_ncbi:
        if nombre_bacteria.lower() in entry['name'].lower():
            return entry['ftp']
    
    # 2. Intento parcial
    nombre_especie = " ".join(nombre_bacteria.split()[:2])
    for entry in datos_ncbi:
        if nombre_especie.lower() in entry['name'].lower() and entry['category'] == 'representative genome':
            return entry['ftp']
            
    return None

def generar_muestra_fastq(directorio_env, archivos_fna):
    """Genera el archivo .fq simulando secuenciación de esos genomas"""
    ruta_fq = os.path.join(directorio_env, "muestra_simulada.fq")
    print(f"   -> Generando muestra sintética en: {ruta_fq} ...")
    
    with open(ruta_fq, "w") as f_out:
        read_id = 0
        for fna_path in archivos_fna:
            seq_total = ""
            try:
                # Leemos el genoma (descomprimiendo al vuelo)
                with gzip.open(fna_path, "rt", encoding='utf-8', errors='ignore') as f_in:
                    for line in f_in:
                        if not line.startswith(">"):
                            seq_total += line.strip().upper()
            except Exception as e:
                print(f"     Error leyendo {fna_path}: {e}")
                continue

            if not seq_total: continue

            # Generamos reads para esta bacteria
            n_reads = min(int(len(seq_total) * COBERTURA / LONGITUD_LECTURA), 50000)
            nombre_bact = os.path.basename(fna_path).split("_")[0]

            for _ in range(n_reads):
                start = random.randint(0, len(seq_total) - LONGITUD_LECTURA - 1)
                read = seq_total[start : start + LONGITUD_LECTURA]
                
                # Introducir error de secuenciación (1%)
                read_list = list(read)
                for i in range(len(read_list)):
                    if random.random() < 0.01:
                        read_list[i] = random.choice("ACGT")
                read_sucio = "".join(read_list)

                f_out.write(f"@MUESTRA_{nombre_bact}_{read_id}\n")
                f_out.write(f"{read_sucio}\n")
                f_out.write("+\n")
                f_out.write("I" * LONGITUD_LECTURA + "\n")
                read_id += 1
    print("   -> ¡Muestra lista!")

def main():
    if not os.path.exists(OUTPUT_BASE): os.makedirs(OUTPUT_BASE)

    print("1. Descargando índice de NCBI (assembly_summary.txt)...")
    
    # MODIFICACIÓN: Usamos wget para mayor robustez
    # Si el archivo existe pero está corrupto (menos de 100MB), lo borramos
    if os.path.exists("summary.txt"):
        if os.path.getsize("summary.txt") < 100 * 1024 * 1024:
            print("   (Archivo previo corrupto detectado, borrando...)")
            os.remove("summary.txt")
            
    if not os.path.exists("summary.txt"):
        # Usamos -q para quiet (menos ruido) y --show-progress para ver la barra
        ret = os.system("wget -q --show-progress https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt -O summary.txt")
        if ret != 0:
            print("❌ Error crítico: wget no pudo descargar la lista maestra.")
            return
    
    print("2. Procesando índice...")
    datos_ncbi = []
    with open("summary.txt", "r") as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.split("\t")
            if parts[11] in ["Complete Genome", "Chromosome"]:
                datos_ncbi.append({
                    "name": parts[7],
                    "ftp": parts[19],
                    "category": parts[4]
                })

    for nombre_env, bacterias in AMBIENTES.items():
        print(f"\n=== Construyendo Entorno: {nombre_env} ===")
        dir_env = os.path.join(OUTPUT_BASE, nombre_env)
        if not os.path.exists(dir_env): os.makedirs(dir_env)
        
        genomas_descargados = []
        
        for bact in bacterias:
            url_ftp = obtener_url_descarga(bact, datos_ncbi)
            if url_ftp:
                filename = url_ftp.split("/")[-1] + "_genomic.fna.gz"
                url_descarga = f"{url_ftp}/{filename}"
                ruta_local = os.path.join(dir_env, filename)
                
                if not os.path.exists(ruta_local):
                    print(f"   Descargando: {bact}...")
                    # MODIFICACIÓN: Usamos wget para las bacterias también
                    ret = os.system(f"wget -q {url_descarga} -O {ruta_local}")
                    if ret == 0:
                        genomas_descargados.append(ruta_local)
                    else:
                        print(f"   ❌ Error descargando {bact}")
                else:
                    print(f"   Ya existe: {bact}")
                    genomas_descargados.append(ruta_local)
            else:
                print(f"   ⚠️ No encontrada en NCBI: {bact}")

        if genomas_descargados:
            generar_muestra_fastq(dir_env, genomas_descargados)
        else:
            print("   Error: No se descargaron genomas para este entorno.")

    print(f"\n✅ PROCESO TERMINADO. Revisa la carpeta '{OUTPUT_BASE}'")

if __name__ == "__main__":
    main()