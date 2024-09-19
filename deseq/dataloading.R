# download.file("https://github.com/carpentries-incubator/bioc-rnaseq/raw/main/episodes/data/GSE96870_coldata_all.csv",
#     "data/GSE96870_coldata_all.csv")
# download.file("https://github.com/carpentries-incubator/bioc-rnaseq/raw/main/episodes/data/GSE96870_rowranges.tsv",
#     "data/GSE96870_rowranges.tsv")
# download.file("https://github.com/carpentries-incubator/bioc-rnaseq/raw/main/episodes/data/GSE96870_coldata_cerebellum.csv",
#     "data/GSE96870_coldata_cerebellum.csv")
 
library(AnnotationDbi)
library(org.Mm.eg.db)
library(hgu95av2.db)
library(SummarizedExperiment)
library(dplyr)

# counts
cerebellum_counts = read.csv("data/GSE96870_counts_cerebellum.csv",
                      row.names = 1)

# sample info
cerebellum_coldata = read.csv("data/GSE96870_coldata_cerebellum.csv",
                              row.names = 1)

# info about each gene
row_ranges = read.delim("data/GSE96870_rowranges.tsv",
                        colClasses = c(ENTREZID = "character"),
                        row.names = 5)

table(row_ranges$gbkey)

# Summarized Experiment

# check that gene / samples in same order
all.equal(rownames(row_ranges), rownames(cerebellum_counts))
all.equal(colnames(cerebellum_counts), rownames(cerebellum_coldata))

# reorder samples to be same as counts

# returns index of first match of item in first arg in second arg
index <- match(colnames(cerebellum_counts), rownames(cerebellum_coldata))
coldata <- cerebellum_coldata[index, ]

# generate summarized experiment

# Granges, S4 wrapper for row_ranges
se <- SummarizedExperiment(assays = list(counts = as.matrix(cerebellum_counts)),
                           rowRanges = as(row_ranges, "GRanges"),
                           colData = cerebellum_coldata)

# get gene names
rownames(se)

# get sample names
colnames(se)

# get counts
assay(se)

# get sample info
colData(se)

# get gene info
table(rowData(se)$gbkey)

# generate label for each sample
se$label <- paste(se$sex, se$time, se$mouse, sep = "_")

# replace GSM with label
colnames(se) = se$label

# generate factors
se$Group = paste(se$sex, se$time, sep = "_") %>% factor(levels = c("Female_Day0",
                                                                   "Male_Day0",
                                                                   "Female_Day4",
                                                                   "Male_Day4",
                                                                   "Female_Day8",
                                                                   "Male_Day8"))
se

# order samples by group
se <- se[, order(se$Group)]

# Make label a factor
se$label <- factor(se$label, levels = se$label)

#How many samples are there for each level of the Infection variable?
table(se$infection)

# Create 2 objects named se_infected and se_noninfected containing a subset of se with only infected and non-infected samples, respectively. Then, calculate the mean expression levels of the first 500 genes for each object, and use the summary() function to explore the distribution of expression levels for infected and non-infected samples based on these genes.

# subset function only for rows
# for column subsetting use select or bracket
se_infected = se[, se$infection == "InfluenzaA"]
se_noninfected = se[, se$infection == "NonInfected"]

# mean expression of first 500 genes
means_infected = rowMeans(assay(se_infected)[1:500,])
means_noninfected = rowMeans(assay(se_noninfected)[1:500,])

summary(means_infected)
summary(means_noninfected)

## cols ($) reference col data names
# How many samples represent female mice infected with Influenza A on day 8?
sum(colData(se_infected)$time == "Day8" & colData(se_infected)$sex == "Female")
sub_se = se[,(se$infection == "InfluenzaA") & (se$time == "Day8") & (se$sex == "Female")]

# save summarized experiment
saveRDS(se, "data/GSE96870_se.RDS")

# output R session info
sessionInfo()

# gene annotation
str(rowData(se))

# mapIds
# these genes are mapped to the entrez id 497097
# alias, all genes mapped to entrez id
# symbol, single gene
# use keytype = ENSEMBL

mapIds(org.Mm.eg.db, keys = "497097", keytype = "ENTREZID", column = "ALIAS", multiVals = "list")


