import os

# Obtener muestras automáticamente
def get_samples():
    return [f.replace(".fastq", "") for f in os.listdir("data") if f.endswith(".fastq")]

SAMPLES = get_samples()


rule all:
    input:
        "quality_control/multiqc/multiqc_report.html",
        expand("quality_control/fastqc/{sample}_fastqc.html", sample=SAMPLES)

rule create_dirs:
    output:
        directory("quality_control/fastqc"),
        directory("quality_control/multiqc"),
        directory("mapping/index_STAR"),
        directory("mapping/alignment_STAR")
    shell:
        "mkdir -p {output}"

include: "rules/quality_control.smk"
include: "rules/alignment_with_genome.smk"