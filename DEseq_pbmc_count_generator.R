library("DESeq2")
library("data.table")
dir.create("../05.DESeq2_analysis")
setwd("../05.DESeq2_analysis")
library("progress")

# annot<-read.table("/share/studies/AMASS/RNA-seq/index_data/Mus_musculus.GRCm38.100.annot.txt",header=T,sep=" ")
# annot<-read.table("/share/studies/Dermatology_Data/index_data/gencode.mm10.annotation.txt",header=T,sep=" ")
annot <- read.table("/share/studies/Dermatology_Data/index_data/gencode.v46.primary_assembly.gene_annotation.txt", header = T, sep = " ")

rnaseq_samples <- read.table("poly_a_sample_list.csv", sep = ",", header = T, row.names = 1)
rnaSeqCount.matrix <- matrix(0, nrow = 63140, ncol = nrow(rnaseq_samples)) # nrow is the total number of rows in the HTSeq read counts file (in this case 55487). ncol is the number of samples
pb <- progress_bar$new(total = nrow(rnaseq_samples), format = "  Progress [:bar] :percent")

for (i in 1:nrow(rnaseq_samples)) {
    sample_id <- row.names(rnaseq_samples)[i]
    dir_id <- paste("../04.Alignment/", sample_id, "_1", "/02_featureCounts2/read_counts_s0_PE.txt", sep = "") # Possible input values for read counting include: read_counts_s0 (unstranded), read_counts_s1 (stranded) and read_counts_s2 (reversely stranded).
    rnaSeq <- data.frame(fread(dir_id, header = T))
    rnaSeqCount.matrix[, i] <- as.vector(rnaSeq[, 7])
    pb$tick()
}
# rnaSeq<-read.table(dir_id,header=T)
rnaSeqCount.matrix <- data.frame(rnaSeqCount.matrix)
row.names(rnaSeqCount.matrix) <- rnaSeq$Geneid
names(rnaSeqCount.matrix) <- row.names(rnaseq_samples)

head(rnaSeqCount.matrix)
dim(rnaSeqCount.matrix)

annot.rnaSeqCount.matrix <- rnaSeqCount.matrix[as.character(annot$ENS_ID), ]
write.table(annot.rnaSeqCount.matrix, "pbmc_poly_a_counts.csv", quote = F, row.names = T, sep = ",")

dim(rnaseq.res2)
rnaseq.res.old <- rnaseq.res2 # 9 here is the number of samples
# keep_cpm <- rowSums(cpm(rnaseq.res.old)>2) >= 2 #  keep transcripts that have minimum one count per million for at least 27 samples (the lowest number of members across groups)
row_res <- rowSums(rnaseq.res.old > 0) > 3
# rnaseq.res<-rnaseq.res.old[keep_cpm,]
rnaseq.res <- rnaseq.res.old[row_res, ]
dim(rnaseq.res)
colnames(rnaseq.res) <- NULL
colnames(rnaseq.res) <- rownames(condition)
# keep<-c("PB1WT","PB7WT","PB3BRAF","PB5BRAF","PB9HD3KO","PB10HD3KO")
# rnaseq.res<-rnaseq.res[,keep]
# dim(rnaseq.res)
# condition<-condition[keep,]
########

dds <- DESeqDataSetFromMatrix(countData = rnaseq.res, colData = condition, design = ~Condition)
dds <- dds[rowSums(counts(dds)) > 30, ]
dim(counts(dds))
design(dds)
dds <- DESeq(dds)
res <- results(dds)
head(res)
save(res, dds, file = "res_all.RData")
dds_counts <- counts(dds, normalized = TRUE)

write.table(dds_counts, "dds_counts_all.txt", quote = F, row.names = T)

res <- results(dds, contrast = c("Condition", "HS", "control"))
write.table(dds_counts, "dds_counts_HS_vs_control.txt", quote = F, row.names = T)

ntd <- normTransform(dds)

png(filename = "boxplot_all.png", width = 12, height = 6, units = "in", res = 600)
par(mar = c(8, 2, 2, 2))
boxplot(log10(assays(dds)[["cooks"]]), range = 0, las = 2)
dev.off()

library("pheatmap")

library("vsn")
vsd <- vst(dds, blind = TRUE)
rld <- rlog(dds, blind = TRUE)
head(assay(vsd), 3)
library(ggplot2)
# pcaData <- plotPCA(assay(vsd), intgroup=c("Condition"), returnData=TRUE)
pca_prcomp <- prcomp(t(assay(vsd)))
x.var <- pca_prcomp$sdev^2
x.pvar <- x.var / sum(x.var)

# percentVar <- round(100 * attr(pcaData, "percentVar"))
percentVar <- round(100 * x.pvar)
Condition <- condition$Condition
png(filename = "PCA_plot_all.png", width = 4, height = 5, units = "in", res = 1200)
ggplot(data.frame(t(assay(vsd))), aes(pca_prcomp$x[, 1], pca_prcomp$x[, 2], color = Condition)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() +
    theme_bw()
dev.off()

############### Comparing HS vs control

res <- results(dds, contrast = c("Condition", "HS", "control"))
res <- data.frame(res)
id <- row.names(res)
res <- cbind(id, res)
res_annot <- merge(annot, res, by.x = "ENS_ID", by.y = "id")
head(res_annot)
rownames(res_annot) <- res_annot$ENS_ID
head(res_annot)
write.table(res_annot, "all_results_HS_vs_control.csv", quote = F, row.names = F, sep = ",")
top_sig <- (res_annot$log2FoldChange > 1 | res_annot$log2FoldChange < (-1)) & res_annot$padj < 0.05 & !is.na(res_annot$padj) & res_annot$baseMean > 10
# top_sig<-(res_annot$log2FoldChange>1 | res_annot$log2FoldChange<(-1)) & res_annot$pvalue <0.05 & !is.na(res_annot$pvalue) & res_annot$baseMean>10

res_annot_top <- res_annot[top_sig, ]
dim(res_annot_top)
write.table(res_annot_top, "DE_results_HS_vs_control.csv", quote = F, row.names = F, sep = ",")
res_annot_HS_vs_control <- res_annot
res_annot_top_HS_vs_control <- res_annot_top

############

topVarGenes <- unique(as.character(row.names(res_annot_top_HS_vs_control)))
library(gplots)
library("RColorBrewer")

condition_select <- condition

# Use this if you want to order the samples in the heatmap
# sample_order<-c("PB1WT","PB2WT","PB7WT","PB6HD3KO","PB9HD3KO","PB10HD3KO","PB3BRAF","PB4BRAF","PB5BRAF")
# sample_order<-c("PB1WT","PB7WT","PB9HD3KO","PB10HD3KO")

# condition_select<-condition[sample_order,]
condition_select <- condition[, ]

# select_matrix<-assay(ntd)[topVarGenes,sample_order]
select_matrix <- assay(ntd)[topVarGenes, ]

rownames(select_matrix) <- as.character(res_annot_top$Gene_Name)

library(magrittr)
library(dplyr)
mat_colors <- list(Condition = c("red", "green"))
names(mat_colors$Condition) <- c("HS", "Control")
test1 <- scale(t(select_matrix))
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "black", "yellow"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(
    seq(min(test1), 0, length.out = ceiling(paletteLength / 2) + 1),
    seq(max(test1) / paletteLength, max(test1), length.out = floor(paletteLength / 2))
)


# Plot the heatmap
png(filename = "heatmap_rnaseq_HS_vs_control.png", width = 10, height = 20, units = "in", res = 600, type = "cairo") # adjust the height if the DE genes are plenty
pheatmap(select_matrix,
    cluster_rows = T, show_rownames = T, color = myColor, breaks = myBreaks, annotation_col = condition,
    annotation_colors = mat_colors,
    cluster_cols = F, scale = "row"
)
dev.off()

pdf("heatmap_rnaseq_HS_vs_control.pdf", width = 10, height = 200) # adjust the height if the DE genes are plenty.
pheatmap(select_matrix,
    cluster_rows = T, show_rownames = T, color = myColor, breaks = myBreaks, annotation_col = condition,
    annotation_colors = mat_colors,
    cluster_cols = T, scale = "row"
)
dev.off()

library(EnhancedVolcano)

png("EnhancedVolcano_HS_vs_control.png", width = 11, height = 11, units = "in", res = 1200, type = "cairo")
EnhancedVolcano(res_annot,
    lab = as.character(res_annot$Gene_Name), pCutoff = 0.05,
    x = "log2FoldChange",
    y = "padj", max.overlaps = 30,
    colAlpha = 1,
    col = c("grey80", "forestgreen", "royalblue", "red2"),
    legendPosition = "right",
    drawConnectors = TRUE,
    gridlines.minor = FALSE,
    gridlines.major = FALSE
)
dev.off()

sessionInfo()
