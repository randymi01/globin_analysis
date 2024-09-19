suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(DESeq2)
  library(vsn)
  library(ggplot2)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(hexbin)
  library(iSEE)
  library(dplyr)
  library(patchwork)
  library(Seurat)
})

se <- readRDS('data/GSE96870_se.RDS')

# QC

# remove all genes with less than or equal to 5 counts across all samples
# can also use counts per million, though 5/10 is an accepted cutoff
# need a minimum unique pieces of evidence for each gene
se <- se[rowSums(assay(se, "counts")) > 5,]

# number of gene counts by sample
se$lib_size <- colSums(assay(se))

# show gene counts by sample by group
colData(se) %>% as.data.frame() %>%
  ggplot(aes(x = label, y = lib_size / 1e6, fill = Group)) +
  geom_bar(stat = "identity") +
  theme_bw() + 
  labs(x = "Sample", y = "count (million)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 5))
# similar group counts

# here our variable factor is time and sex
dds <- DESeq2::DESeqDataSet(se, design = ~ sex + time)

# measure representing frequency of gene
dds <- estimateSizeFactors(dds)

ggplot(data.frame(libSize = colSums(assay(dds)),
                  sizeFactor = sizeFactors(dds),
                  Group = dds$Group),
       aes(x = libSize, y = sizeFactor, col = Group)) + 
  geom_point(size = 5) + 
  theme_bw() + 
  labs(x = "Library size", y = "Size factor")

# not homoskedastic
meanSdPlot(assay(dds), ranks = F)

# normalize using vsd
vsd <- DESeq2::vst(dds, blind = TRUE)

# now normalized, so can use glms
meanSdPlot(assay(vsd), ranks = F)

# euclidean distance matrix (can be used since data is normalized)
dst <- dist(t(assay(vsd)))

colors <- colorRampPalette(brewer.pal(9, "Blues"))(255)
ComplexHeatmap::Heatmap(
  as.matrix(dst), 
  col = colors,
  name = "Euclidean\ndistance",
  cluster_rows = hclust(dst),
  cluster_columns = hclust(dst),
  bottom_annotation = columnAnnotation(
    sex = vsd$sex,
    time = vsd$time,
    col = list(sex = c(Female = "red", Male = "blue"),
               time = c(Day0 = "yellow", Day4 = "forestgreen", Day8 = "purple")))
)

# pretty good clusters
# sex is linearlly seperable along pc 1
# pc 2 seperates date well
pcaData <- DESeq2::plotPCA(vsd, intgroup = c("sex", "time"),
                           returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = sex, shape = time), size = 5) +
  theme_minimal() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() + 
  scale_color_manual(values = c(Male = "blue", Female = "red"))

# look at dendogram on heatmap for sources of variance
# first branches have the most variance

sce <- as(dds, "SingleCellExperiment")


# only difference between singlecell and summarized experiment is an
# extra slot for lower dimensional representation
## Add PCA to the 'reducedDim' slot
stopifnot(rownames(pcaData) == colnames(sce))
reducedDim(sce, "PCA") <- as.matrix(pcaData[, c("PC1", "PC2")])

## Add variance-stabilized data as a new assay
stopifnot(colnames(vsd) == colnames(sce))
assay(sce, "vsd") <- assay(vsd)

app <- iSEE(sce)
shiny::runApp(app)


df <- readRDS('data/HS_merged_lite.RDS')
hs <- as(df, "SingleCellExperiment")
app1 <- iSEE(hs)
shiny::runApp(app1)
hs <- as.SingleCellExperiment(df)
