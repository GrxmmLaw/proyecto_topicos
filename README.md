# Proyecto Tópicos: Analisis de rareza y novedad usando HyperLogLog

**Integrantes:** Maximiliano Lopez - Iván Zapata

## Archivos de pruebas

Las siguientes carpetas de datos **no están incluidas en este repositorio** debido a su tamaño:

Descargala: [Este enlace](https://drive.google.com/drive/folders/1cUjYGXsyV7hDEPq4DcF2QRbBSNzZ4wtz?usp=sharing)

Hay que colocar estas carpetas en la raíz del proyecto.

---

## Descripción General

Este proyecto implementa un sistema de generación y análisis de datos sintéticos de genomas. El objetivo es crear catálogos de secuencias genómicas simuladas con diferentes grados de similitud, extraer k-mers característicos de cada catálogo y generar muestras (reads) para análisis de anomalías. El proyecto funciona en tres etapas principales:

---

## Etapa 1: Generación de Catálogos Sintéticos

### Descripción

El script `catalog_simulator.py` genera tres catálogos de genomas sintéticos con diferentes características:

#### **Configuración Base**

- **Tamaño objetivo:** 1.0 GB por catálogo
- **Longitud de genoma:** 4 millones de pares de bases (bp)
- **Contenido GC:** 50% (Guanina-Citosina)
- **Número de genomas:** ~250 por catálogo (calculado automáticamente)

#### **Tipos de Catálogos Generados**

| Catálogo      | Descripción              | Cantidad | Características                         |
| ------------- | ------------------------ | -------- | --------------------------------------- |
| **DISTINTOS** | Completamente aleatorios | ~250     | Cada genoma es completamente diferente  |
| **SIMILARES** | Variantes de una semilla | ~250     | Mutación del 0.1% respecto a la semilla |
| **MEDIOS**    | Variantes medianas       | ~250     | Mutación del 3.0% respecto a la semilla |

### Uso

```bash
# Ejecutar desde la carpeta simulation/
python catalog_simulator.py
```

**Salida esperada:**

```
=== Generador de catlogos sinteticos ===
Carpeta base: sim_catalogs/
Configuración: 1.0 GB por set.
Genomas a generar: 250 por set.
--------------------------------------------------

[1/3] Generando 250 genomas distintos
Progreso: [250/250] 100.0% - Restante: 0.0 min

[2/3] Generando 250 genomas SIMILARES (Mutación 0.1%)...
Progreso: [250/250] 100.0% - Restante: 0.0 min

[3/3] Generando 250 genomas MEDIOS (Mutación 3.0%)...
Progreso: [250/250] 100.0% - Restante: 0.0 min

¡Proceso finalizado! Tiempo total: 12.45 minutos.
Los archivos están en '/ruta/al/proyecto/sim_catalogs'
```

**Archivos generados:**

```
sim_catalogs/
├── datasets_distintos/
│   ├── genome_diff_0001.fna
│   ├── genome_diff_0002.fna
│   └── ... (248 más)
├── datasets_similares/
│   ├── genome_strain_0000_SEED.fna
│   ├── genome_strain_0001.fna
│   └── ... (248 más)
└── datasets_medios/
    ├── genome_spp_0000_SEED.fna
    ├── genome_spp_0001.fna
    └── ... (248 más)
```

---

## Etapa 2: Construcción de Sketches HyperLogLog

### Descripción

El script `build_catalog.py` procesa todos los genomas de cada catálogo para extraer k-mers y crear un sketch HyperLogLog (HLL) que representa la "firma" del catálogo.

### Configuración

- **Tamaño de k-mer:** 31 nucleótidos
- **Tasa de error HLL:** 1% (0.01)
- **Salida:** Archivos `.sketch` (binarios, serializados con pickle)

### Uso

```bash
# Ejecutar desde la carpeta simulation/
python build_catalog.py
```

**Salida esperada:**

```
--- Generador de catalogos HLL ---
K-mer size: 31
Error rate: 0.01
Carpeta de salida: hll_sim_catalog

============================================================
PROCESANDO CATALOGO: sim_catalogs/datasets_distintos
============================================================
-> Se encontraron 250 genomas en la carpeta.
-> Iniciando Pool con 7 núcleos...
   Progreso: 250/250 genomas...
-> Completado en 45.32 segundos (0.181 s/genoma).
-> Cardinalidad Estimada (K-mers únicos): 1,000,245,367

============================================================
PROCESANDO CATALOGO: sim_catalogs/datasets_similares
============================================================
-> Se encontraron 250 genomas en la carpeta.
-> Iniciando Pool con 7 núcleos...
   Progreso: 250/250 genomas...
-> Completado en 42.15 segundos (0.169 s/genoma).
-> Cardinalidad Estimada (K-mers únicos): 987,123,456

============================================================
PROCESANDO CATALOGO: sim_catalogs/datasets_medios
============================================================
-> Se encontraron 250 genomas en la carpeta.
-> Iniciando Pool con 7 núcleos...
   Progreso: 250/250 genomas...
-> Completado en 43.98 segundos (0.176 s/genoma).
-> Cardinalidad Estimada (K-mers únicos): 1,102,567,890

--- Tiempo total: 3.16 minutos ---
```

**Archivos generados:**

```
hll_sim_catalog/
├── distintos.sketch     (HLL del catálogo DISTINTOS)
├── similares.sketch     (HLL del catálogo SIMILARES)
└── medios.sketch        (HLL del catálogo MEDIOS)
```

---

## Etapa 3: Generación de Muestras (Reads)

### Descripción

El script `sample_simulator.py` genera muestras de secuenciación (reads) basadas en los catálogos. Crea dos tipos de muestras:

- **CONTROL:** Solo contiene ruido (onda senoidal de 0-5% de rareza)
- **TEST:** Contiene un evento anómalo (anomalía triangular del 0-10% de rareza)

### Configuración

| Parámetro           | Valor         | Descripción                   |
| ------------------- | ------------- | ----------------------------- |
| **Total reads**     | 6,000,000     | Cantidad de reads generados   |
| **Longitud read**   | 150 bp        | Longitud de cada read         |
| **Window size**     | 200,000 reads | Ventanas para análisis        |
| **K-size**          | 31            | Tamaño de k-mer               |
| **Ciclos de ruido** | 4             | Repeticiones de onda senoidal |
| **Max ruido**       | 5%            | Rareza máxima del ruido       |

### Combinaciones Generadas

El script genera **6 muestras** (combinación de 3 catálogos × 2 modos):

1. **DISTINTOS_CONTROL** - Reads del catálogo DISTINTOS con solo ruido
2. **DISTINTOS_TEST** - Reads DISTINTOS con evento anómalo (usa catálogo MEDIOS como anomalía)
3. **MEDIOS_CONTROL** - Reads del catálogo MEDIOS con solo ruido
4. **MEDIOS_TEST** - Reads MEDIOS con evento anómalo (usa catálogo DISTINTOS como anomalía)
5. **SIMILARES_CONTROL** - Reads del catálogo SIMILARES con solo ruido
6. **SIMILARES_TEST** - Reads SIMILARES con evento anómalo (usa catálogo DISTINTOS como anomalía)

### Uso

```bash
# Ejecutar desde la carpeta simulation/
python sample_simulator.py
```

**Salida esperada:**

```
--- Generando DISTINTOS [CONTROL] ---
   Generada V30...
   -> Guardado: muestras_validacion/DISTINTOS_CONTROL.fastq

--- Generando DISTINTOS [TEST] ---
   Generada V30...
   -> Guardado: muestras_validacion/DISTINTOS_TEST.fastq

--- Generando MEDIOS [CONTROL] ---
   Generada V30...
   -> Guardado: muestras_validacion/MEDIOS_CONTROL.fastq

--- Generando MEDIOS [TEST] ---
   Generada V30...
   -> Guardado: muestras_validacion/MEDIOS_TEST.fastq

--- Generando SIMILARES [CONTROL] ---
   Generada V30...
   -> Guardado: muestras_validacion/SIMILARES_CONTROL.fastq

--- Generando SIMILARES [TEST] ---
   Generada V30...
   -> Guardado: muestras_validacion/SIMILARES_TEST.fastq

TIEMPO TOTAL: 15.32 min.
```

---

## Etapa 4: Análisis de Rareza y Novedad

### Descripción

El script `windowed_rarity_analyzer.py` analiza las muestras generadas calculando la rareza de k-mers en ventanas deslizantes. Utiliza los sketches HyperLogLog para detectar anomalías y calcular métricas de rendimiento.

### Configuración

| Parámetro               | Valor         | Descripción                         |
| ----------------------- | ------------- | ----------------------------------- |
| **Window size**         | 200,000 reads | Ventanas para análisis              |
| **K-size**              | 31            | Tamaño de k-mer                     |
| **HLL Error rate**      | 0.01          | Tasa de error del sketch            |
| **Calibration windows** | 5             | Ventanas iniciales para calibración |
| **Sigma multiplier**    | 2.0           | Multiplicador para umbral dinámico  |

### Uso

```bash
# Ejecutar desde la carpeta simulation/
python windowed_rarity_analyzer.py
```
