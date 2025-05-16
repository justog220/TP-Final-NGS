# ---- Instalación de paquetes Bioconductor ----

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("biomaRt")

# ---- Importar librerías ----

library(DESeq2)
library(biomaRt)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(ggrepel)

#---- Archivos de recuentos de entrada - STAR ----
# Definimos carpetas de entrada y salida
input_dir  <- "differential_expression/full_data_counts/counts_STAR"
output_dir <- "differential_expression/counts_STAR_selected"
# Aseguramos de que existe la carpeta de salida
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
# Listamos solo los ficheros que nos interesan
files <- list.files(input_dir, pattern = "\\.tab$", full.names = TRUE)
# Iteramos
for (archivo in files) {
  # Leemos saltándonos las 4 líneas iniciales y nombrando las 4 columnas
  data <- read_tsv(archivo,
                   skip = 4,
                   col_names = c("Gene_ID", "Unstranded", "Forward", "Reverse"))
  # Construimos el nombre de muestra a partir del nombre de fichero
  sample <- str_split(basename(archivo), "_")[[1]][1]
  # Exportamos solo Gene_ID + Unstranded
  out_file <- file.path(output_dir, paste0(sample, "_1_counts.txt"))
  write_tsv(data[, c("Gene_ID", "Unstranded")],
            out_file,
            col_names = FALSE)
    message("Escrito: ", out_file)
}

#---- Tabla de muestras ----

sampletable <- data.frame(SampleName=c("5p4_25c", "5p4_27c", "5p4_28c", "5p4_29c", "5p4_30c", "5p4_31cfoxc1", "5p4_32cfoxc1", "5p4_33cfoxc1", "5p4_34cfoxc1", "5p4_35cfoxc1"),
                          FileName=c("SRR3091420_1_counts.txt", "SRR3091421_1_counts.txt", "SRR3091422_1_counts.txt", "SRR3091423_1_counts.txt", "SRR3091424_1_counts.txt", "SRR3091425_1_counts.txt", "SRR3091426_1_counts.txt", "SRR3091427_1_counts.txt", "SRR3091428_1_counts.txt", "SRR3091429_1_counts.txt"),
                          Differentiation=c(rep("undiff", 2), rep("diff5days", 3), rep("undiff", 3), rep("diff5days", 2)),
                          Condition=c(rep("WT", 5), rep("KO", 5)))
rownames(sampletable) <- gsub("_counts.txt", "", sampletable$FileName)

#---- Importar recuentos STAR ----

setwd("./differential_expression/")
se_star <- DESeqDataSetFromHTSeqCount(sampleTable = sampletable,
                                      directory = "counts_STAR_selected",
                                      design = ~ Condition) # El diseño en base al cual va a comparar.

# ----Filtrado de genes de baja expresión----

#Numero de genes antes del filtrado:
nrow(se_star)
#Filtrado (por conteo, que tengan por lo menos más de 10 reads):
se_star <- se_star[rowSums(counts(se_star)) > 10, ]
#Número de genes luego del filtrado:
nrow(se_star)

# ----Preparar anotación----

#Lista de archivos ENSEMBL
listEnsemblArchives()
#En esta práctica se usa la versión 98 de ENSEMBL, correspondiente a la URL https://sep2019.archive.ensembl.org
#Cargar la base de datos correspondiente
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="https://sep2019.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
#Cabecera de la lista de filtros para obtener la anotación
head(listFilters(mart))
#Se buscan los disponibles provenientes de ensembl
grep("ensembl", listFilters(mart)[,1], value=TRUE)
#ensembl_gene_id se encuentra!
#"attributes" corresponde al tipo de anotación que se quiere obtener
head(listAttributes(mart))
#Podes buscar los que te interesen. En esta practica, usaremos 'ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position', 'description', 'external_gene_name'
#Lista de ENSEMBL IDs que queremos anotar
gene_ids <- rownames(se_star)
#18025 IDs
#Anotar!
annot <- getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position', 'description', 'external_gene_name'), filters ='ensembl_gene_id', values = gene_ids, mart = mart)
dim(annot)
#17933 renglones
head(annot)

#---- Ajustar modelo estadístico DESeq2----

se_star2 <- DESeq(se_star)
#Calcular conteo normalizado (transformación log2); + 1
norm_counts <- log2(counts(se_star2, normalized = TRUE)+1)
#Agregar anotación
norm_counts_symbols <- merge(data.frame(ID=rownames(norm_counts), norm_counts, check.names=FALSE), annot, by.x="ID", by.y="ensembl_gene_id", all=F)
#Escribir recuentos normalizados en un archivo
write.table(norm_counts_symbols, "normalized_counts_log2_star.txt", quote=F, col.names=T, row.names=F, sep="\t")



#----Análisis de expresión diferencial----

#Chequeo de nombre de resultados: depende de qué se haya modelado. Aquí se modeló la "Condition"
resultsNames(se_star2)
#Extraer resultados para WT vs KO
de <- results(object = se_star2,
              name="Condition_WT_vs_KO")
#Es equivalente a:
# de <- results(object = se_star2, contrast=c("Condition", "WT", "KO"))

#Si quieren ver los resultados expresados como "KO vs WT", ejecutar:
# de <- results(object = se_star2, contrast=c("Condition", "KO", "WT"))

#Chequeo de filas
head(de)
#Agregar la anotación
de_symbols <- merge(data.frame(ID=rownames(de), de, check.names=FALSE), annot, by.x="ID", by.y="ensembl_gene_id", all=F)
#Escribir los resultados de expresión diferencial en el archivo
write.table(de_symbols, "deseq2_results.txt", quote=F, col.names=T, row.names=F, sep="\t")

 

#---- Visualización----

# Transformar los recuentos sin procesar para poder visualizar los datos
se_rlog <- rlog(se_star2)

#---- Visualización - Correlación de muestras----

#Calcular la matriz de distancia entre muestras
sampleDistMatrix <- as.matrix(dist(t(assay(se_rlog))))
#Preparar metadata
metadata <- sampletable[,c("Differentiation", "Condition")]
rownames(metadata) <- sampletable$SampleName
#Gráfico
pheatmap(sampleDistMatrix, annotation_col=metadata)

#-----Visualización - Análisis de componentes principales (PCA)----
plotPCA(object = se_rlog,
        intgroup = c("Condition", "Differentiation"))

#---- Visualización - FOXC1 es ENSG00000054598----

# Obtener recuentos normalizados para FOXC1 / ENSG00000054598
tmp <- norm_counts[rownames(norm_counts)=="ENSG00000054598",]
#Convertir a formato largo
mygenelong <- melt(tmp)
#Nombre de muestra
mygenelong$name <- rownames(mygenelong)
#Condición de muestra y diferenciación: mezclar
mygenelong <- merge(mygenelong, sampletable, by.x="name", by.y="SampleName", all=F)
# Dot plot
pdot <- ggplot(data=mygenelong, mapping=aes(x=Condition, y=value, col=Differentiation, shape=Condition, label=name)) +
  geom_point() +
  geom_text(nudge_x=0.2) +
  xlab(label="Experimental group") +
  ylab(label="Normalized expression (log2)") +
  theme_bw()
pdot

#---- Visualización - Volcano plot----

#Agregar columna de diferenciación
de_symbols$diffexpressed <- "NO"
#si log2Foldchange > 0.5 y pvalue < 0.05, setear "UP"
de_symbols$diffexpressed[de_symbols$log2FoldChange > 0.5 & de_symbols$pvalue < 0.05] <- "UP"
#si log2Foldchange < -0.5 y pvalue < 0.05, setear "DOWN"
de_symbols$diffexpressed[de_symbols$log2FoldChange < -0.5 & de_symbols$pvalue < 0.05] <- "DOWN"
# Crear columna con los Ids de los genes que tienen expresión diferencial
de_symbols$delabel <- NA
de_symbols$delabel[de_symbols$diffexpressed != "NO"] <- de_symbols$external_gene_name[de_symbols$diffexpressed != "NO"]
#Graficar
ggplot(data = de_symbols,
       aes(
         x = log2FoldChange,
         y = -log10(pvalue),
         color = diffexpressed,
         label = delabel
       )) +
  geom_point() +
  geom_text_repel() +
  scale_color_manual(values = c("blue", "black", "red")) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "red") +
  geom_hline(yintercept = -log10(0.05), col = "red") +
  labs(
    x = "log2FoldChange",
    y = "-log10(pvalue)",
    title = "Volcano plot",
    subtitle = "Se observa el -log10(pvalue) en función del log2FoldChange",
    caption = "Seminario NGS",
    color = "Expresión"
  ) +
  theme(
    legend.background = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 2
    )
  )