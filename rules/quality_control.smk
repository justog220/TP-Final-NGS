rule fastqc:
    input:
        os.path.join(SAMPLES_DIR, "{sample}.fastq")
    output:
        html = "results/quality_control/fastqc/{sample}_fastqc.html",
        zip = "results/quality_control/fastqc/{sample}_fastqc.zip"
    shell:
        "fastqc {input} --outdir results/quality_control/fastqc/ -q"

rule multiqc:
    input:
        expand("results/quality_control/fastqc/{sample}_fastqc.html", sample=SAMPLES)
    output:
        "results/quality_control/multiqc/multiqc_report.html"
    shell:
        "multiqc results/quality_control/fastqc/ --outdir results/quality_control/multiqc"