#####################################################################
# STAT 540: Seminar 8
# Methylation analysis
# Date: Wed Mar 5
# By: Jessica Pilsworth
#####################################################################

# load libraries 
library(GEOquery)
library(wateRmelon)
library(IlluminaHumanMethylation450k.db)

# load data
if (file.exists("methyl_ALL.Rdata")) {
  # if previously downloaded
  load("methyl_ALL.Rdata")
} else {
  # if downloading for the first time
  GSE39141 <- getGEO("GSE39141")
  show(GSE39141)  ## 33 samples (29 ALL and 4 healthy B cells)
  GSE42865 <- getGEO("GSE42865") 
  show(GSE42865)  ## 16 samples (9 healthy cells B cells and 7 other cells)
  
  # Extract expression matrices (turn into data frames at once)
  ALL.dat <- as.data.frame(exprs(GSE39141[[1]]))
  CTRL.dat <- as.data.frame(exprs(GSE42865[[1]]))
  
  # Obtain the meta-data for the samples and rename them perhaps?
  ALL.meta <- pData(phenoData(GSE39141[[1]]))
  CTRL.meta <- pData(phenoData(GSE42865[[1]]))
  
  # create some labels
  ALL.meta$Group <- c(rep("ALL", 29), rep("HBC", 4))
  ## ALL: Case; HBC: Healthy B Cells
  
  # Subset both meta-data and data for control (healthy) donors
  CTRL.meta <- droplevels(subset(CTRL.meta, grepl("Healthy donor", characteristics_ch1.1)))
  CTRL.dat <- subset(CTRL.dat, select = as.character(CTRL.meta$geo_accession))
  
  # Rename variables
  names(ALL.dat) <- paste(ALL.meta$Group, gsub("GSM", "", names(ALL.dat)), 
                          sep = "_")
  names(CTRL.dat) <- paste("HBC", gsub("GSM", "", names(CTRL.dat)), sep = "_")
  
  # save the data to avoid future re-downloading
  save(ALL.dat, CTRL.dat, ALL.meta, CTRL.meta, file = "methyl_ALL.Rdata")
}

# load library
library(ggplot2)

# to calculate the average Beta values for the probes in the two datasets 
pmDat <- c(rowMeans(ALL.dat, na.rm = T), rowMeans(CTRL.dat, na.rm = T)) 

# to create data.frame for density plot
probe.meansDat <- data.frame(Beta = pmDat,
                             Dataset = rep(c('ALL', 'CTRL'), each = nrow(ALL.dat)))

# to create a density plot of the probe mean beta values 
(probe.meansPlot <- ggplot(data = probe.meansDat, aes(x = Beta, col = Dataset)) +
    geom_density() + ggtitle("Density Plot of Mean Beta Values") + 
    xlab("Beta") + ylab("Density"))

# Normalization

# combine data from two experiments into one matrix, each column represents beta values of one sample
beta.matrix <- as.matrix(cbind(ALL.dat, CTRL.dat))
str(beta.matrix, max.level = 0)

# to perform quantile normalization
# to calculate normalized betas from Illumina 450K methylation arrays
system.time(beta.norm <- betaqn(beta.matrix))

# to calculate row means
probe.meansDatNorm <- c(rowMeans(beta.norm[, 1:ncol(ALL.dat)], na.rm = TRUE),
                    rowMeans(beta.norm[, ncol(ALL.dat):ncol(beta.norm)], na.rm = TRUE)) 

# to create a normalized plot 
probe.meansNormPlot <- rbind(data.frame(probe.meansDat, Norm = "Before"),
                  data.frame(Beta = pmDat,
                   Dataset = rep(c('ALL', 'CTRL'), each = nrow(ALL.dat)),
                   Norm = "After"))

# to add a level
probe.meansNormPlot$Norm <- factor(probe.meansNormPlot$Norm, levels = c("Before", "After"))

# to create a density plot before and after normalization
(probe.meansNorm <- ggplot(data = probe.meansNormPlot, aes(x = Beta, col = Dataset)) +
   geom_density() + facet_grid(Norm ~ .) + 
   ggtitle("Density Plot of Mean Beta Values Before and After Normalization") + 
   xlab("Beta") + ylab("Density"))

# M values

# to apply a logit transformation on the data to convert it to a continuous variable that spans (−∞,∞), 
# i.e. compute the M value
M.norm <- beta2m(beta.norm)


# CpG Islands

# to extract probe ID to CpG islands association
cginame <- as.data.frame(IlluminaHumanMethylation450kCPGINAME)
names(cginame) <- c("Probe_ID", "cginame")
rownames(cginame) <- cginame$Probe_ID
length(levels(factor(cginame$cginame))) 

# to restrict probes to those within CGIs
beta.inCGI <- beta.norm[cginame$Probe_ID, ]
M.inCGI <- M.norm[cginame$Probe_ID, ]
nrow(M.inCGI) 

# to aggregate probes to CGIs
beta.CGI <- aggregate(beta.inCGI, by = list(cginame$cginame), mean, na.rm = T)
rownames(beta.CGI) <- beta.CGI[, "Group.1"]
beta.CGI <- subset(beta.CGI, select = -Group.1)
str(beta.CGI, max.level = 0)

M.CGI <- aggregate(M.inCGI, by = list(cginame$cginame), mean, na.rm = T)
rownames(M.CGI) <- M.CGI[, "Group.1"]
M.CGI <- subset(M.CGI, select = -Group.1)
str(M.CGI, max.level = 0)

# to reshape the data
library(reshape2)
M.CGI.tall <- melt(t(M.CGI), value.name = "M", varnames = c("Sample", "CGI"))
M.CGI.tall$Group <- substring(M.CGI.tall$Sample, 1, 3)

# to create of bloxplot of CGI M values
(M.boxplot <- ggplot(data = M.CGI.tall, aes(Sample, M, colour= Group)) + 
   geom_boxplot() + ggtitle("Bloxplot of CGI M values") + 
   xlab("Samples") + ylab("M Values") + 
     scale_x_discrete(labels = NULL))

# Differential methylation analysis with limma
library(limma)
design <- data.frame(Group = relevel(factor(gsub("_[0-9]+", "", colnames(M.CGI))), 
                                     ref = "HBC"), row.names = colnames(M.CGI))
str(design)
(DesMat <- model.matrix(~Group, design))
DMRfit <- lmFit(M.CGI, DesMat)
DMRfitEb <- eBayes(DMRfit)
cutoff <- 0.01
DMR <- topTable(DMRfitEb, coef = "GroupALL", number = Inf, p.value = cutoff)
head(DMR) # top hits

# create plots to check these hits
library(gplots)
DMR100 <- topTable(DMRfitEb, coef = "GroupALL", number = 100)
DMR.CGI <- t(as.matrix(subset(beta.CGI, rownames(beta.CGI) %in% rownames(DMR100))))
str(DMR.CGI, max.level = 0)

# create heatmap of beta values of top 100 hits
col <- c(rep("darkgoldenrod1", times = nrow(DMR.CGI)))
col[grepl("HBC", rownames(DMR.CGI))] <- "forestgreen"
op <- par(mai = rep(0.5, 4))
#heatmap.2(DMR.CGI, col = redblue(256), RowSideColors = col, density.info = "none", 
 #         trace = "none", Rowv = TRUE, Colv = TRUE, labCol = FALSE, labRow = FALSE, 
  #        dendrogram = "row", margins = c(1, 5))
#legend("topright", c("ALL", "HBC"), col = c("darkgoldenrod1", "forestgreen"), 
 #      pch = 15)

# to create a stripplot of beta values of probes within top 5 CGI hits
DMR5 <- topTable(DMRfitEb, coef = "GroupALL", number = 5)
beta.DMR5probe <- beta.inCGI[cginame[rownames(beta.inCGI), ]$cginame %in% rownames(DMR5), 
                             ]
beta.DMR5probe.tall <- melt(beta.DMR5probe, value.name = "M", varnames = c("Probe_ID", 
                                                                           "Sample"))
beta.DMR5probe.tall$Group <- factor(gsub("_[0-9]+", "", beta.DMR5probe.tall$Sample))
beta.DMR5probe.tall$CGI <- factor(cginame[as.character(beta.DMR5probe.tall$Probe_ID), 
                                          ]$cginame)
(beta.DMR5.stripplot <- ggplot(data = beta.DMR5probe.tall, aes(x = Group, y = M, 
                                                               color = Group)) + geom_point(position = position_jitter(width = 0.05), na.rm = T) + 
   stat_summary(fun.y = mean, aes(group = 1), geom = "line", color = "black") + 
   facet_grid(. ~ CGI) + ggtitle("Probe beta values within top 5 DM CGIs") + 
   xlab("Group") + ylab("beta") + theme_bw())

# to plot location of differential methylated probes along each chromosome
# get the length of chromosome 1-22 and X
chrlen <- unlist(as.list(IlluminaHumanMethylation450kCHRLENGTHS)[c(as.character(1:22), 
                                                                   "X")])
chrlen <- data.frame(chr = factor(names(chrlen)), length = chrlen)
chr <- IlluminaHumanMethylation450kCHR  # get the chromosome of each probe
# get the probe identifiers that are mapped to chromosome
chr <- unlist(as.list(chr[mappedkeys(chr)]))
# get chromosome coordinate of each probe
coord <- IlluminaHumanMethylation450kCPGCOORDINATE
# get the probe identifiers that are mapped to coordinate
coord <- unlist(as.list(coord[mappedkeys(coord)]))
coord <- data.frame(chr = chr[intersect(names(chr), names(coord))], coord = coord[intersect(names(chr), 
                                                                                            names(coord))])
# coordinates of probes in DM CGIs
coordDMRprobe <- droplevels(na.omit(coord[cginame[cginame$cginame %in% rownames(DMR), 
                                                  ]$Probe_ID, ]))
(coord.plot <- ggplot(data = coordDMRprobe) + geom_linerange(aes(factor(chr, 
                                                                        levels = c("X", as.character(22:1))), ymin = 0, ymax = length), data = chrlen, 
                                                             alpha = 0.5) + geom_point(aes(x = factor(chr, levels = c("X", as.character(22:1))), 
                                                                                           y = coord), position = position_jitter(width = 0.03), na.rm = T) + ggtitle("DMR positions on chromosomes") + 
   ylab("Position of DMRs") + xlab("chr") + coord_flip() + theme_bw())