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

rule multiqc:
    input:
        expand("quality_control/fastqc/{sample}_fastqc.html", sample=SAMPLES)
    output:
        "quality_control/multiqc/multiqc_report.html"
    conda:
        "environment.yml"
    shell:
        "multiqc quality_control/fastqc/ --outdir quality_control/multiqc"