#####################################################################
# STAT 540: Seminar 7
# RNA-seq Analysis: From BAM File To Count Data
# Date: Wed Feb 26
# By: Jessica Pilsworth
#####################################################################
 
# load libraries
library(ShortRead)
library(Rsamtools)
library(easyRNASeq)

# load data
bamDat <- readAligned("drosophilaMelanogasterSubset.bam", 
                      type = "BAM")

str(bamDat)

# to index a bam file
indexFile <- indexBam("drosophilaMelanogasterSubset.bam")

# to filter reads with at most 2 N's
nFilt <- nFilter(2)

# to filter reads which have been aligned to the reference genome
chrFilt <- chromosomeFilter(regex = "chr")

# to create a new filter which checks that both conditions are satisfied 
# i.e. a read has 2 or fewer N's and is aligned to the reference genome
filt <- compose(nFilt, chrFilt)

# to apply the filter to extract the relevant subset of data
bamDatFiltered <- bamDat[filt(bamDat)]

# to examine the filtered bam file
str(bamDatFiltered)

# to look at which chromosomes are present and their read IDs
levels(chromosome(bamDatFiltered))

id(bamDatFiltered)[1:10]

# to view the DNA sequence for the read
sread(bamDatFiltered)[1:10]

# to look at the base qualities
quality(bamDatFiltered)[1:10]

# to find out the starting (left most) position where the read was aligned on the chromosome
position(bamDatFiltered)[1:10]

# to see which strand (forward/+ or reverse/-) the read was aligned to 
strand(bamDatFiltered)[1:10]

# load genomic data
library(BSgenome.Dmelanogaster.UCSC.dm3)

# to get the lengths of the chromosomes in the genome
(chrSizes <- seqlengths(Dmelanogaster))

# load library and specifiy Drosophila database
library(biomaRt)

ensembl <- useMart("ensembl", dataset = "dmelanogaster_gene_ensembl")

# to define a set of fields we are interested for each gene
annotation.fields <- c("ensembl_gene_id", "strand", "chromosome_name", "start_position", 
                       "end_position")

# to restrict to chr2L
gene.annotation <- getBM(annotation.fields, mart = ensembl, filters = "chromosome_name", 
                         values = c("2L"))

str(gene.annotation)

# to check we only downloaded annotations from chr2L
levels(as.factor(gene.annotation$chromosome))

# to add chr to the annotation data
gene.annotation$chromosome <- paste("chr", gene.annotation$chromosome_name, 
                                    sep = "")

levels(as.factor(gene.annotation$chromosome))

# to store the gene annotation information in an IRanges object
gene.range <- RangedData(IRanges(start = gene.annotation$start_position, end = gene.annotation$end_position), 
                         space = gene.annotation$chromosome, strand = gene.annotation$strand, gene = gene.annotation$ensembl_gene_id, 
                         universe = "Dm3")

show(gene.range)

# to find out how many bases cover each position in every chromosome
(cover <- coverage(bamDatFiltered, width = chrSizes))

# to aggregate the coverage for each gene
gene.coverage <- aggregate(cover[match(names(gene.range), names(cover))], ranges(gene.range), 
                           sum)

# to find the number of reads covering each gene
gene.coverage <- ceiling(gene.coverage/unique(width(bamDat)))
gene.coverage

# to restrict to chr2L
length(gene.coverage[["chr2L"]])

length(ranges(gene.range)$chr2L)

# to build a count table and store it in a data frame
countTable <- data.frame(chromosome = gene.range$space, gene_start = start(gene.range$ranges), 
                         gene_end = end(gene.range$ranges), strand = gene.range$strand, gene = gene.range$gene, 
                         count = as.vector(gene.coverage[["chr2L"]]))

dim(countTable)

head(countTable)

# to add RPKM
countTable <- data.frame(chromosome = gene.range$space, gene_start = start(gene.range$ranges), 
gene_end = end(gene.range$ranges), strand = gene.range$strand, gene = gene.range$gene, 
count = as.vector(gene.coverage[["chr2L"]]), RPKM = (as.vector(gene.coverage[["chr2L"]])/(end(gene.range$ranges) - 
                                                                                            start(gene.range$ranges))) * (1e+09/length(bamDat)))
head(countTable)
