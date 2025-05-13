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

## 🗂️ Estructura del proyecto

⚠️ TODO: Completar una vez definida la estructura del directorio.



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

##  👥 Autor
<a href="https://github.com/justog220"><img src="https://avatars.githubusercontent.com/u/85772318?v=4" title="justog220" width="50" height="50"></a>