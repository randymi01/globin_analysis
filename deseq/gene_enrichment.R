library(SummarizedExperiment)
library(DESeq2)
library(gplots)
library(microbenchmark)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(simplifyEnrichment)

se <- readRDS('data/GSE96870_se.RDS')

# enrichment analysis = over representation analysis
# are more genes in a gene set over represented in a specific group
# ie are there more sex genes over represented in the Differentially Expressed group

# remove non mrna genes
se <- se[rowData(se)$gbkey == 'mRNA',]

# convert se into deseq2 dataset
dds <- DESeq2::DESeqDataSet(se, design = ~ sex + time)
dds <- DESeq(dds)
resSex <- results(dds, contrast = c("sex", "Male", "Female"))
resSex

# get subset of these genes which have sufficient fdr adj p-value
sexDE = data.frame(resSex[resSex$padj <= 0.05 & !isNA(resSex$padj),])

DEgenes <- rownames(sexDE)

# is there association between DE genes and sex

# grange with chromosome and bp range, find only DEgenes
geneGR <- rowRanges(se)
totalGenes <- rownames(se)
XYGeneSet <- totalGenes[as.vector(seqnames(geneGR) %in% c("X","Y"))]

# total number of filtered genes
length(totalGenes)

# total number of sex genes
length(XYGeneSet)

# 13/54 DE genes are in sex chromosomes
library(gplots)
plot(venn(list("sexDEgenes" = DEgenes,
               "XY gene set" = XYGeneSet)))
title("Universe = 21,198")

## Fisher's exact test ----
# measure whether DE gene existence in specified group is statistically significant or not
n <- nrow(se)
n_01 <- length(XYGeneSet)
n_10 <- length(DEgenes)
n_11 <- length(intersect(DEgenes, XYGeneSet))
n_12 <- n_10 - n_11
n_21 <- n_01 - n_11
n_20 <- n    - n_10
n_02 <- n    - n_01
n_22 <- n_02 - n_12
matrix(c(n_11, n_12, n_10, n_21, n_22, n_20, n_01, n_02, n),
       nrow = 3, byrow = TRUE)

# p value represents likelihood of having more DE genes in n_11 (sex gene set)
fisher.test(matrix(c(n_11, n_12, n_21, n_22), nrow = 2, byrow = TRUE),
            alternative = "greater")

## phyper function - hypergeometric distribution
# q - the observed
# m - number of degenes
# n - number of non-degenes
# k - m+n

# same p-value as fisher.test
1 - phyper(n_11 - 1, m = n_10, n = n_20, k = n_01)


lt <- list(gene_set_1 = c("gene_1", "gene_2", "gene_3"),
           gene_set_2 = c("gene_1", "gene_3", "gene_4"),
           gene_set_3 = c("gene_4", "gene_7", "gene_5")
)

# df of gene set and gene
df = data.frame(gene_set = rep(names(lt), times = sapply(lt, length)),
                gene = unname(unlist(lt)))

## gene ontology databases ----
# BP = biological process
# CC = cellular component
# MF - molecular function

library(org.Mm.eg.db)
keytypes(org.Mm.eg.db)

# multivals when multiple keys
# returns the equivalent values from column
# GO: https://geneontology.org/ (bio niacs)
BP_Id <- mapIds(org.Mm.eg.db, keys = c("BP"), keytype = "ONTOLOGY", column = "GO", multiVals = "list")[[1]]

BPGeneSets <- mapIds(org.Mm.eg.db, keys = BP_Id, keytype = "GOALL", column = "ENTREZID", multiVals = "list")

# KEGG - pathway
keggNames = read.table(url("https://rest.kegg.jp/list/pathway/hsa"), sep = "\t")
head(keggNames)

# remove hsa start string
keggNames = gsub("hsa","",keggNames[,1])
keggNames = gsub("","",keggNames[,2])

# clusterprofiler ----
# same thing we did for sex
# define genes that change between day0 and day8
resTime <- DESeq2::results(dds, contrast = c("time", "Day8", "Day0"))
timeDE <- as.data.frame(subset(resTime, 
                               padj < 0.05 & abs(log2FoldChange) > log2(1.5)
))
timeDEgenes <- rownames(timeDE)
head(timeDEgenes)

# automates process from above
library(clusterProfiler)
# expects entrezid but our genes are symbols, use keytype arg
# only outputs signficant genes
enric <- enrichGO(gene = timeDEgenes,
         ont = "BP",
         OrgDb = org.Mm.eg.db,
         keyType = "SYMBOL")


resTimeGOTable = as.data.frame(enric)
head(resTimeGOTable)

## KEGG enrichment analysis
# keys, what you have, keytype, keytype of what your keys are
# column, what you want to output
EntrezIDs <- mapIds(org.Mm.eg.db, keys = timeDEgenes, keytype = "SYMBOL",
                    column = "ENTREZID")

## remove NA
EntrezIDs <- EntrezIDs[!is.na(EntrezIDs)]

resTimeKEGG = enrichKEGG(gene = EntrezIDs, 
                         organism = "mmu",
                         pvalueCutoff = 1,
                         qvalueCutoff = 1)
resTimeKEGGTable = as.data.frame(resTimeKEGG)
head(resTimeKEGGTable)
