#####################################################################
# STAT 540: Seminar 7
# RNA-seq Analysis: Differential expression analysis
# Date: Wed Feb 26
# By: Jessica Pilsworth
#####################################################################

# load library
library(edgeR)

# load data
dat <- read.table("bottomly_count_table.tsv", header = TRUE, row.names = 1)
des <- read.table("bottomly_phenodata.tsv", header = TRUE, row.names = 1)
str(dat)
show(des)
all(rownames(des) == colnames(dat))

# to create a 'group' object describing which group each sample belongs to
with(des, table(strain))
group <- factor(c(rep("1", 10), rep("2", 11)))
group

# to produce an object of type DGEList 
dge.glm <- DGEList(counts = dat, group = group)
str(dge.glm)

# DGELIst object has two components, 
# one is a matrix call 'counts' storing the count data and 
# the other is a data.frame called 'samples' storing information for samples.
names(dge.glm)
dge.glm[["samples"]]
nrow(dge.glm[[1]])
ncol(dge.glm[[1]])

# create design model
design <- model.matrix(~group)
design

# to estimate common dispersion for negative binomial GLMs
dge.glm.com.disp <- estimateGLMCommonDisp(dge.glm, design, verbose = TRUE)

# to estimate trended dispersion for negative binomial GLMs
dge.glm.trend.disp <- estimateGLMTrendedDisp(dge.glm.com.disp, design)

# empirical bayes tagwise dispersions ofr negative binomial GLMs
dge.glm.tag.disp <- estimateGLMTagwiseDisp(dge.glm.trend.disp, design)

# to plot the tagwise dispersion against log2-CPM (counts per million)
plotBCV(dge.glm.tag.disp)

# to fit the genewise negative binomial generlaized linear model
fit <- glmFit(dge.glm.tag.disp, design)
colnames(coef(fit))
lrt <- glmLRT(fit, coef = 2)

# table of the top differentially expressed tags
topTags(lrt)
tt.glm <- topTags(lrt, n = Inf)
class(tt.glm)

# number of rows in the table
nrow(tt.glm$table[tt.glm$table$FDR < 0.01, ])

# to find the samples with an FDR less than 1-e50
interestingSamples <- rownames(tt.glm$table[tt.glm$table$FDR < 1e-50, ])
cpm(dge.glm.tag.disp)[interestingSamples, ]

# summary of multiple testing across genes and contrasts
summary(de.glm <- decideTestsDGE(lrt, p = 0.05, adjust = "BH"))

# to plot the tagwise log fold changes against log-cpm
tags.glm <- rownames(dge.glm.tag.disp)[as.logical(de.glm)]
plotSmear(lrt, de.tags = tags.glm)
abline(h = c(-2, 2), col = "blue")

# load library 
library(DESeq)

# read the same count table data and grouping information
deSeqDat <- newCountDataSet(dat, group)
head(counts(deSeqDat))

# to estimate the size factors to account for differences in library coverage and estimate the variance
deSeqDat <- estimateSizeFactors(deSeqDat)
sizeFactors(deSeqDat)

# to estimate dispersions
deSeqDat <- estimateDispersions(deSeqDat)

# to plot the estimated dispersions against the mean normalized counts
plotDispEsts(deSeqDat)

# to fit the model and examine the results
results <- nbinomTest(deSeqDat, levels(group)[1], levels(group)[2])
str(results)

# to plot the results
plotMA(results)

# load library
library(limma)

# to calculate normalization factors to align columns of a count matrix & plot
norm.factor <- calcNormFactors(dat)
dat.voomed <- voom(dat, design, plot = TRUE, lib.size = colSums(dat) * norm.factor)

dat.voomed

# to fit a linear model 
fit <- lmFit(dat.voomed, design)
fit <- eBayes(fit)

# to get top hits
topTable(fit)
