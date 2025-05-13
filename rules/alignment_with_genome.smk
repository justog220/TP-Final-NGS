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
        fastq = "data/reads/{sample}.fastq",
        index_flag = "results/mapping/index_STAR/.index_complete"
    output:
        bam = "results/mapping/alignments_STAR/{sample}Aligned.sortedByCoord.out.bam",
        counts = "results/mapping/alignments_STAR/{sample}ReadsPerGene.out.tab",
        log = "results/mapping/alignments_STAR/{sample}Log.final.out"
    params:
        prefix = "results/mapping/alignments_STAR/{sample}",
        genomeDir = "results/mapping/index_STAR"
    threads: 8
    shell:
        """
        STAR --genomeDir {params.genomeDir} \
             --readFilesIn {input.fastq} \
             --outSAMtype BAM SortedByCoordinate \
             --quantMode GeneCounts \
             --outFileNamePrefix {params.prefix} \
             --runThreadN {threads}
        """
