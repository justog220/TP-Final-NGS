import os

# Obtener muestras automáticamente
def get_samples():
    return [f.replace(".fastq", "") for f in os.listdir("data/reads") if f.endswith(".fastq")]

def get_replics():
    return set([sample.split("_")[0] for sample in get_samples()])

def get_paired_samples():
    samples = sorted(get_samples())

    paired_samples = [samples[i:i+2] for i in range(0, len(samples), 2)]

    paired_samples_dict = {}
    for pair in paired_samples:
        key = pair[0].split("_")[0]
        paired_samples_dict[key] = [pair[0], pair[1]]

    return paired_samples_dict


SAMPLES = get_samples()

REPLICS = get_replics()

SAMPLES_DIR = "data/reads"


rule qc_all:
    input:
        "results/quality_control/multiqc/multiqc_report.html",

rule mapping_all:
    input:
        expand("results/mapping/alignments_STAR/{replic}/{replic}Aligned.sortedByCoord.out.bam", replic=REPLICS),
        expand("results/mapping/alignments_STAR/{replic}/{replic}ReadsPerGene.out.tab", replic=REPLICS),
        expand("results/mapping/alignments_STAR/{replic}/{replic}Coverage.bw", replic=REPLICS),
        expand("results/qc_qualimap/{replic}/qualimapReport.html", replic=REPLICS),
        "results/igv_visualization/igv_report.html"


rule create_dirs:
    output:
        directory("results/quality_control/fastqc"),
        directory("results/quality_control/multiqc"),
        directory("results/mapping/index_STAR"),
        directory("results/mapping/alignment_STAR"),
        directory("results/mapping_quality_control")
    shell:
        "mkdir -p {output}"

include: "rules/quality_control.smk"
include: "rules/alignment_with_genome.smk"