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
