import os

# Obtener muestras automáticamente
def get_samples():
    return [f.replace(".fastq", "") for f in os.listdir("data/reads") if f.endswith(".fastq")]

SAMPLES = get_samples()

SAMPLES_DIR = "data/reads"


rule qc_all:
    input:
        "results/quality_control/multiqc/multiqc_report.html",

rule mapping_all:
    input:
        expand("results/quality_control/fastqc/{sample}_fastqc.html", sample=SAMPLES),
        expand("results/mapping/alignments_STAR/{sample}Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        expand("results/mapping/alignments_STAR/{sample}ReadsPerGene.out.tab", sample=SAMPLES)

rule create_dirs:
    output:
        directory("results/quality_control/fastqc"),
        directory("results/quality_control/multiqc"),
        directory("results/mapping/index_STAR"),
        directory("results/mapping/alignment_STAR")
    shell:
        "mkdir -p {output}"

include: "rules/quality_control.smk"
include: "rules/alignment_with_genome.smk"