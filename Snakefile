import os

# Obtener muestras automáticamente
def get_samples():
    return [f.replace(".fastq", "") for f in os.listdir("data") if f.endswith(".fastq")]

SAMPLES = get_samples()

print(SAMPLES)
# --- Reglas --- #
rule all:
    input:
        "quality_control/multiqc/multiqc_report.html",
        expand("quality_control/fastqc/{sample}_fastqc.html", sample=SAMPLES)  # <- ¡Añadido!

# Crear directorios (opcional, ya que Snakemake los crea automáticamente)
rule create_dirs:
    output:
        directory("quality_control/fastqc"),
        directory("quality_control/multiqc")
    shell:
        "mkdir -p {output}"

# Regla FastQC
rule fastqc:
    input:
        "data/{sample}.fastq"
    output:
        html = "quality_control/fastqc/{sample}_fastqc.html",
        zip = "quality_control/fastqc/{sample}_fastqc.zip"
    conda:
        "environment.yml"
    shell:
        "fastqc {input} --outdir quality_control/fastqc/"

# Regla MultiQC
rule multiqc:
    input:
        expand("quality_control/fastqc/{sample}_fastqc.html", sample=SAMPLES)
    output:
        "quality_control/multiqc/multiqc_report.html"
    conda:
        "environment.yml"
    shell:
        "multiqc quality_control/fastqc/ --outdir quality_control/multiqc"