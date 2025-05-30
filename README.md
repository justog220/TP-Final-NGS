# Análisis de Expresión Génica en *Saccharomyces cerevisiae* bajo Diferentes Condiciones de Crecimiento

**Autor:** Justo Garcia
**Seminario:** Análisis de Secuencias de Nuevas Tecnologías de Secuenciación en Paralelo
**Docente:** Ileana Tossolini
**Carrera:** Licenciatura en Bioinformática - Facultad de Ingeniería, Universidad Nacional de Entre Ríos

## 🎯 Objetivo

Este repositorio contiene el análisis bioinformático del transcriptoma de la levadura *Saccharomyces cerevisiae* (cepa CEN.PK 113-7D) bajo dos condiciones metabólicas:
- Batch (respiro-fermentativo)
- Quimiostato (completamente respiratorio)

Se busca identificar genes diferencialmente expresados y realizar análisis funcionales con un enfoque reproducible usando Snakemake y RMarkdown.

## 🧰 Tecnologías y herramientas

- Snakemake – para orquestar el pipeline.
- Conda – para gestión de entornos y dependencias.
- R / RMarkdown – para análisis de expresión diferencial y generación de reportes.
- Quarto – para generar reportes HTML interactivos. Apartir de RMd

## 🗂️ Estructura del proyecto

```plaintext
[4.0K]  ./
├── [4.0K]  data/ # Directorio de datos
│   ├── [4.0K]  reads/ # Lecturas de secuenciación
│   │   ├── [9.3M]  batch1_1.fastq
│   │   ├── [8.6M]  batch1_2.fastq
│   │   ├── [ 13M]  batch2_1.fastq
│   │   ├── [ 12M]  batch2_2.fastq
│   │   ├── [9.9M]  batch3_1.fastq
│   │   ├── [9.1M]  batch3_2.fastq
│   │   ├── [8.1M]  quimiostato1_1.fastq
│   │   ├── [7.6M]  quimiostato1_2.fastq
│   │   ├── [9.9M]  quimiostato2_1.fastq
│   │   ├── [9.4M]  quimiostato2_2.fastq
│   │   ├── [ 12M]  quimiostato3_1.fastq
│   │   └── [ 11M]  quimiostato3_2.fastq
│   ├── [4.0K]  reference/ # Genoma de referencia y anotaciones
│   │   ├── [229K]  sacCer_ChrI.fa
│   │   ├── [  20]  sacCer_ChrI.fa.fai
│   │   └── [5.1M]  sacCer_genes.gtf
│   └── [4.0K]  templates/ # Plantillas para el reporte de IGV
│       ├── [ 507]  igv_base.html
│       └── [  66]  igv_config.js
├── [4.0K]  notebooks/ # Notebooks de análisis de expresión diferencial
│   ├── [4.0K]  differential_expression_files/
│   ├── [209K]  differential_expression.html # Reporte HTML del análisis
│   ├── [ 24K]  differential_expression.qmd # Archivo Quarto para el análisis
│   └── [1.8M]  graphics_generation.ipynb # Generación de gráficos para el artículo
├── [4.0K]  results/ # Resultados del pipeline
├── [4.0K]  rules/ # Reglas de Snakemake
│   ├── [7.3K]  alignment_with_genome.smk # Reglas relativas al mapeo y su exploración
│   └── [ 592]  quality_control.smk # Reglas para el control de calidad
├── [196K]  Consigna Levadura.pdf # Consigna del trabajo práctico
├── [ 279]  environment.yml # Archivo de entorno de Conda
├── [2.6K]  README.md # Este archivo
└── [1.6K]  Snakefile # Archivo principal de Snakemake
```


## 🚀 Cómo ejecutar el análisis

### 1. Clonar el repositorio

```bash
git clone https://github.com/justog220/TP-Final-NGS.git
cd TP-Final-NGS
```

### 2. Crear el entorno de Conda

```bash
conda env create -f environment.yml
```

Activar el entorno creado:

```bash
conda activate tp-final-ngs
```

### 3. Ejecutar etapas del pipeline

#### Control de calidad

```bash
snakemake --cores all qc_all
```

#### Mapeo al genoma de referencia

```bash
snakemake --cores all mapping_all
```

## 📊 Análisis de expresión diferencial

El análisis de expresión diferencial se realiza en R mediante un archivo `.Rmd` ubicado en `notebooks/`, utilizando DESeq2 o edgeR.

Este paso genera visualizaciones, listas de genes diferencialmente expresados y otros análisis exploratorios.

## 📦 Material incluido

- Scripts y reglas de Snakemake
- Archivos de configuración
- Notebooks de análisis en RMarkdown
- Instrucciones de ejecución reproducible
- Documentación del análisis

## 🧪 Datos

Los datos utilizados corresponden a lecturas *paired-end* alineadas al cromosoma I de *S. cerevisiae*, accesibles en SRA bajo el número de acceso SRS307298.
Por cuestiones de tiempo de cómputo, se utilizará un set reducido de lecturas previamente filtradas.

## :warning: Aclaracion

Algunos análisis no se incluyen en este pipeline, sino que están disponibles a través del informe del trabajo práctico.

##  👥 Autor
<a href="https://github.com/justog220"><img src="https://avatars.githubusercontent.com/u/85772318?v=4" title="justog220" width="50" height="50"></a>