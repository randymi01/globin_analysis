# Load packages
library(org.Mm.eg.db)
library(hgu95av2.db)
library(SummarizedExperiment)
library(DESeq2)
library(tidyverse)
library(cowplot)
library(ExploreModelMatrix)
library(ComplexHeatmap)

# load annotationdbi select last
AnnotationDbi::select()

# load data
se <- readRDS("data/GSE96870_se.RDS")

# make deseqdataset
# last value in design is the variable of interest
dds <- DESeqDataSet(se, design = ~ sex + time)

# Normalization
dds <- estimateSizeFactors(dds)

# Estimate dispersion parameter for negative binomial distribution
dds <- estimateDispersions(dds)

# pops out new graphics window
x11()

# from graph looks like higher frequency genes are less dispersed
# estimatedispersions moves the dispersion of genes closer to the fitted mean
# genes with low dispersion are moved up which

# stein's phenomonon, add bias term when number of covariates >= 4
plotDispEsts(dds)

# do differential expression testing
# fits a negative binomial glm to each gene
# derives p value from each model using wald test

# negative binom test outputs:
# log fold 2 change is statistic measuring change between 2 populations
# baseMean
dds <- nbinomWaldTest(dds)

# equivalent to running
dds <- DESeq(dds)

resultsNames(dds)

# get results for one pairwise comparison
res <- results(dds)

# check what the wald test is testing
# does by default the order of the factors
# log2 fold change (MLE): time Day8 vs Day0 (numerator / denominator)
res

# get explicit comparisons (factors) (use name for continuous)
# contrast (variable, level in numerator, level in denominator)
resTime <- results(dds, contrast = c('time', 'Day8', 'Day0'))

# padj using FDR adjustment
# additional genes are removed if not enough counts for fdr correction
summary(resTime)

# order by raw since FDR has no tie breaker

resTime[order(resTime$pvalue),] %>% head()

# nice p value curve with no spikes mean you have a good chunk of
# differentially expressed genes
hist(resTime$padj, 1000)

dds <- DESeqDataSet(se, design = ~ time + sex)
dds <- DESeq(dds)
res_sex <- results(dds)
res_sex[order(res_sex$pvalue),]


dds <- DESeqDataSet(se, design = ~ sex + time)
dds <- DESeq(dds)
# results is where params for statistical tests are set
resTime <- results(dds)

# Visualize results
# blue dots are significant
# points should be centered on zero
# can get lines on upper and lower portions (fanning out) that are not significant
# probably because increaing from low number
# can get groupings at left side since mean is too low
plotMA(resTime)

# shrink log fold changes
resTimeLfc <- lfcShrink(dds, coef = "sex_Male_vs_Female",
                        res = res_sex)

x11()
plotMA(resTimeLfc)

# results
results(dds)

# Challenge
# RE-TEST resTime, but change independent filtering to false
# not filtering out genes before hypothesis test make the test
# more conservative since the multiple measures tests punishes
# more measures when calculating FDR adjusted p value
resTimeNoFilter <- results(dds, independentFiltering = F)
summary(resTime)
summary(resTimeNoFilter)

# heatmap
# vst is for visualization, not for stats
vsd <- vst(dds, blind = TRUE)
genes <- resTime[order(resTime$pvalue),]|> head(10) |> rownames()

heatmapData <- assay(vsd)[genes,]
View(heatmapData)
# scale data first for heatmap
hist(heatmapData)

# scale works by col data, so must first transpose since we are interested
# in scaling the gene data
heatmapData <- t(heatmapData) |> scale() |> t()

# add annotation

# coldata is a bioc DFrame, must transform to normal dataframe
heatmapColAnnot <- data.frame(colData(vsd)[, c("time", "sex")])
heatmapColAnnot <- HeatmapAnnotation(df = heatmapColAnnot)

x11()

ComplexHeatmap::Heatmap(heatmapData,
                       top_annotation = heatmapColAnnot,
                       cluster_rows = TRUE,
                       cluster_columns = FALSE)
scale(heatmap.T)

# In heatmap, gene expression starts to change during day 4, but lots
# of difference by day 8
# genes on bottom get upregulated as disease continues

# output results
resTime

rowRanges(se)

temp <- cbind(data.frame(rowRanges(se)), data.frame(resTime))
# issues with reading in excel (date genes)
write.csv(temp, 'outputs/Day8vsDay0.csv')

# alternative, write to tab delimited txt file, then read into excel setting dtype by col
temp <- tibble::rownames_to_column(temp, var = "gene")
write.table(temp, file = "outputs/Day8vsDay0.txt",
            row.names = F, sep = "\t")

# glimma interactive MA plot
library(Glimma)

# gene explorer
glimmaMA(dds, html = "outputs/maplotDay8vsDay0.html")

# showing PCs
glimmaMDS(dds)
