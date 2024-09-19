suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(ExploreModelMatrix)
  library(dplyr)
  library(DESeq2)
})

# pull in metadata ----
meta <- read.csv('data/GSE96870_coldata_all.csv', row.names = 1)

# table with labels
with(meta, table(time, infection))

vs = VisualizeDesign(sampleData = meta, designFormula = ~tissue + time + sex)
vs

x11()

meta

vs$plotlist

meta_f <- meta %>% filter(tissue == "Spinalcord" & time == "Day0")

male_spine <- filter(meta, sex == "Male", tissue == "Spinalcord")

meta_non_inf <- filter(meta, infection == "NonInfected")

# plot list shows the formula for mean of each group
# shows observation table
vd = VisualizeDesign(sampleData = meta, designFormula = ~sex + tissue)
vd

# Model Design
# interactive
ExploreModelMatrix(sampleData = meta, design = ~sex + tissue)

# interaction terms
cd <- VisualizeDesign(sampleData = meta_non_inf, designFormula = ~sex + tissue + sex:tissue)
cd$designmatrix

# make combined term with 4 levels
meta_noninf <- meta %>% filter(time == "Day0")
meta_noninf$sex_tissue <- paste0(meta_noninf$sex, "_", meta_noninf$tissue)
meta_noninf
vd <- VisualizeDesign(sampleData = meta_noninf, 
                      designFormula = ~sex_tissue)
vd$plotlist

# paired design, each mouse contributed cerebellum and spinal cord. Treat mouse as factor
# to take into account and get info
cd <- VisualizeDesign(sampleData = meta_non_inf, designFormula = ~sex + factor(mouse) + sex:tissue)
cd


meta_fem_day0 = meta %>% filter(time == "Day0", sex == "Female")
meta_fem_day4 = meta %>% filter(time == "Day4", sex == "Female")
meta_fem_not_day8 = meta %>% filter(time != "Day8", sex == "Female")

vd0 <- VisualizeDesign(sampleData = meta_fem_day0, designFormula = ~ factor(mouse) + tissue)
vd4 <- VisualizeDesign(sampleData = meta_fem_day4, designFormula = ~ factor(mouse) + tissue)
vd8 <- VisualizeDesign(sampleData = meta_fem_not_day8, designFormula = ~ factor(mouse) + tissue)

# show groups and formulas by group VisualizeDesign(sampleData = , designFormula = )
# show design matrix model.matrix(designFormula, data = )
md0 <- model.matrix(~ mouse, data = meta_fem_day0)
md4 <- model.matrix(~mouse, data = meta_fem_day4)
mdn8 <- model.matrix(~mouse, data = meta_fem_not_day8)

# read back in summarized experiment
se <- readRDS("data/GSE96870_se.RDS")

# do qc again
se <- se[rowSums(assay(se, "counts"))>5,]
dds <- DESeq2::DESeqDataSet(se, design = ~ sex + time)
dds <- DESeq(dds)
res_dds <- results(dds, contrast = c("time", "Day8", "Day0"))

# model matrix made by DESEQ2 same as design matrix from above even though DESEQ uses negbinom instead
# of glm
attr(dds, "modelMatrix")

vd <- VisualizeDesign(sampleData = colData(dds)[, c("sex", "time")], 
                      designMatrix = attr(dds, "modelMatrix"),
                      flipCoordFitted = T)
vd



