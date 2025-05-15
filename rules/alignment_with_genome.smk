import math
import Bio.SeqIO

rule make_star_index:
    input:
        expand("data/reads/{sample}.fastq", sample=SAMPLES),
        genome_reference = "data/reference/sacCer_ChrI.fa",
        gtf_reference = "data/reference/sacCer_genes.gtf"
    output:
        touch("results/mapping/index_STAR/.index_complete")  # Archivo bandera
    params:
        genomeDir="results/mapping/index_STAR",
        genomeFastaFiles=input.genome_reference,
        sjdbGTFfile=input.gtf_reference,
        sjdbOverhang=48, #TODO: revisar, la profe en la guía dice que es igual al menor tamaño de lectura menos uno, esto confirmaria la necesidad de filtrar por longitud
        genomeSAindexNbases=lambda wildcards, input: min(14, math.log2(sum(len(record) for record in Bio.SeqIO.parse(input.genome_reference, "fasta")))/2 - 1),
        outFileNamePrefix="sacCer_chr1",
    threads: 4
    shell:
        """
        STAR --runMode genomeGenerate --genomeDir {params.genomeDir} \
                    --genomeFastaFiles {params.genomeFastaFiles} \
                    --sjdbGTFfile {params.sjdbGTFfile} \
                    --sjdbOverhang {params.sjdbOverhang} \
                    --genomeSAindexNbases {params.genomeSAindexNbases} \
                    --outFileNamePrefix {params.outFileNamePrefix} \
                     --runThreadN {threads}
        """

rule align_STAR:
    input:
        # pair_of_fastq = lambda wildcards: [f"data/reads/{replic}_1.fastq", f"data/reads/{replic}_2.fastq"],
        fastq_1 = "data/reads/{replic}_1.fastq",
        fastq_2 = "data/reads/{replic}_2.fastq",
        index_flag = "results/mapping/index_STAR/.index_complete"
    output:
        bam = "results/mapping/alignments_STAR/{replic}/{replic}Aligned.sortedByCoord.out.bam",
        counts = "results/mapping/alignments_STAR/{replic}/{replic}ReadsPerGene.out.tab",
        log = "results/mapping/alignments_STAR/{replic}/{replic}Log.final.out"
    params:
        prefix = "results/mapping/alignments_STAR/{replic}/{replic}",
        genomeDir = "results/mapping/index_STAR"
    threads: 4
    shell:
        """
        STAR --genomeDir {params.genomeDir} \
             --readFilesIn {input.fastq_1},{input.fastq_2} \
             --outSAMtype BAM SortedByCoordinate \
             --quantMode GeneCounts \
             --outFileNamePrefix {params.prefix} \
             --runThreadN {threads}
        """

rule quality_control_mapping:
    input:
        alignment="results/mapping/alignments_STAR/{sample}/{sample}Aligned.sortedByCoord.out.bam",
        gtf="data/reference/sacCer_genes.gtf",
        gene_counts="results/mapping/alignments_STAR/{sample}/{sample}ReadsPerGene.out.tab"
    output:
        "results/qc_qualimap/{sample}/qualimapReport.html"
    params:
        outdir="results/qc_qualimap/{sample}",

    shell:
            """
            protocol=$(grep -v "N_" {input.gene_counts} | \
                                awk '{{unst+=$2;forw+=$3;rev+=$4}} \
                                    END{{if (forw > 2*rev) print "strand-specific-forward"; \
                                        else if (rev > 2*forw) print "strand-specific-reverse"; \
                                        else print "non-strand-specific"}}')

            qualimap rnaseq -bam {input.alignment} \
                                                -gtf {input.gtf} \
                                                -outdir {params.outdir} \
                                                -p $protocol \
                                                -pe
            """

rule indexing_bam:
    input:
        bam_file="results/mapping/alignments_STAR/{replic}/{replic}Aligned.sortedByCoord.out.bam"
    output:
        "results/mapping/alignments_STAR/{replic}/{replic}Aligned.sortedByCoord.out.bam.bai"
    shell:
        """
        samtools index {input.bam_file}
        """

rule calculate_coverage:
    input:
        bam_file="results/mapping/alignments_STAR/{replic}/{replic}Aligned.sortedByCoord.out.bam",
        bai_file="results/mapping/alignments_STAR/{replic}/{replic}Aligned.sortedByCoord.out.bam.bai"
    output:
        coverage_file="results/mapping/alignments_STAR/{replic}/{replic}Coverage.bw"
    shell:
        """
        bamCoverage -b {input.bam_file} -o {output.coverage_file}
        """

rule generate_igv_tracks:
    input:
        bws = expand("results/mapping/alignments_STAR/{replic}/{replic}Coverage.bw", replic=REPLICS),
        bams = expand("results/mapping/alignments_STAR/{replic}/{replic}Aligned.sortedByCoord.out.bam", replic=REPLICS),
        bais = expand("results/mapping/alignments_STAR/{replic}/{replic}Aligned.sortedByCoord.out.bam.bai", replic=REPLICS),
        fasta_reference = "data/reference/sacCer_ChrI.fa",
        fai_reference = "data/reference/sacCer_ChrI.fa.fai",
        gtf_reference = "data/reference/sacCer_genes.gtf"
    output:
        "results/igv_visualization/igv_tracks.json"
    run:
        import json

        tracks_config = {
            "genome": {
                "fastaURL": input.fasta_reference,
                "indexURL": input.fai_reference,
                "gtfURL": input.gtf_reference
            },
            "tracks": []
        }

        for replic in REPLICS:
            tracks_config["tracks"].append({
                "name": replic,
                "bw": f"results/mapping/alignments_STAR/{replic}/{replic}Coverage.bw",
                "bam": f"results/mapping/alignments_STAR/{replic}/{replic}Aligned.sortedByCoord.out.bam",
                "bai": f"results/mapping/alignments_STAR/{replic}/{replic}Aligned.sortedByCoord.out.bam.bai"
            })

        with open(output[0], 'w') as f:
            json.dump(tracks_config, f, indent=4)


rule generate_igv_report:
    input:
        bws = expand("results/mapping/alignments_STAR/{replic}/{replic}Coverage.bw", replic=REPLICS),
        bams = expand("results/mapping/alignments_STAR/{replic}/{replic}Aligned.sortedByCoord.out.bam", replic=REPLICS),
        fasta = "data/reference/sacCer_ChrI.fa",
        gtf = "data/reference/sacCer_genes.gtf",
        template_base = "data/templates/igv_base.html",
        template_config = "data/templates/igv_config.js"
    output:
        html = "results/igv_visualization/igv_report.html",
        config = "results/igv_visualization/igv_config.js"
    run:
        import json
        import shutil
        import base64

        # 1. Generar configuración JSON
        config = {
            "genome": {
                "id": "yeast_custom",
                "name": "Yeast Chromosome",
                "fastaURL": "data:base64," + base64.b64encode(open(input.fasta, "rb").read()).decode(),
                "tracks": [{
                    "name": "Annotations",
                    "type": "annotation",
                    "format": "gtf",
                    "url": "data:base64," + base64.b64encode(open(input.gtf, "rb").read()).decode(),
                    "displayMode": "EXPANDED"
                }]
            },
            "tracks": [
                {
                    "name": f"{replic} Coverage",
                    "type": "wig",
                    "format": "bigwig",
                    "url": "data:base64," + base64.b64encode(open(f"results/mapping/alignments_STAR/{replic}/{replic}Coverage.bw", "rb").read()).decode()
                }
                for replic in sorted(REPLICS)
            ]
        }

        # 2. Guardar configuración como JS
        with open(output.config, 'w') as f:
            f.write(f"const igvConfig = {json.dumps(config, indent=2)};")

        # 3. Copiar plantilla base y modificar
        shutil.copy2(input.template_base, output.html)

        # 4. Inyectar configuración directamente en el HTML
        with open(output.html, 'a') as f_html:
            with open(output.config, 'r') as f_config:
                f_html.write(f"\n<script>\n{f_config.read()}\n</script>")