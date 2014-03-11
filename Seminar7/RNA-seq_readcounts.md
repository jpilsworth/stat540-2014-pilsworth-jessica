# STAT 540: Seminar 7
# RNA-seq Analysis: Getting Read Counts
# Date: Wed Feb 26
# By: Jessica Pilsworth
========================================================

Load libraries

```r
library(ShortRead)
```

```
## Loading required package: BiocGenerics
## Loading required package: parallel
## 
## Attaching package: 'BiocGenerics'
## 
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
## 
## The following object is masked from 'package:stats':
## 
##     xtabs
## 
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, as.vector, cbind,
##     colnames, duplicated, eval, evalq, Filter, Find, get,
##     intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rep.int, rownames, sapply, setdiff, sort,
##     table, tapply, union, unique, unlist
## 
## Loading required package: IRanges
## Loading required package: GenomicRanges
## Loading required package: XVector
## Loading required package: Biostrings
## Loading required package: lattice
## Loading required package: Rsamtools
```

```r
library(Rsamtools)
library(easyRNASeq)
```

```
## Loading required package: genomeIntervals
## Loading required package: intervals
## 
## Attaching package: 'intervals'
## 
## The following object is masked from 'package:Biostrings':
## 
##     type
## 
## The following object is masked from 'package:GenomicRanges':
## 
##     reduce
## 
## The following objects are masked from 'package:IRanges':
## 
##     expand, reduce
```

```
## Warning: replacing previous import by 'intervals::coerce' when loading 'genomeIntervals'
## Warning: replacing previous import by 'intervals::initialize' when loading 'genomeIntervals'
```

```
## Loading required package: Biobase
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
## 
## Loading required package: biomaRt
## Loading required package: edgeR
## Loading required package: limma
## 
## Attaching package: 'limma'
## 
## The following object is masked from 'package:BiocGenerics':
## 
##     plotMA
## 
## Loading required package: DESeq
## Loading required package: locfit
## locfit 1.5-9.1 	 2013-03-22
## 
## Attaching package: 'locfit'
## 
## The following objects are masked from 'package:GenomicRanges':
## 
##     left, right
## 
##     Welcome to 'DESeq'. For improved performance, usability and
##     functionality, please consider migrating to 'DESeq2'.
```


Load data

```r
bamDat <- readAligned("drosophilaMelanogasterSubset.bam", type = "BAM")
str(bamDat)
```

```
## Formal class 'AlignedRead' [package "ShortRead"] with 8 slots
##   ..@ chromosome  : Factor w/ 15 levels "chrYHet","chrM",..: 2 2 2 2 2 2 2 2 2 2 ...
##   ..@ position    : int [1:64206] 548 1497 1506 1528 1540 1552 1552 1555 1559 1566 ...
##   ..@ strand      : Factor w/ 3 levels "+","-","*": 2 1 1 1 1 1 1 1 2 2 ...
##   ..@ alignQuality:Formal class 'NumericQuality' [package "ShortRead"] with 1 slots
##   .. .. ..@ quality: int [1:64206] 132 132 127 130 130 122 132 132 132 132 ...
##   ..@ alignData   :Formal class 'AlignedDataFrame' [package "ShortRead"] with 4 slots
##   .. .. ..@ varMetadata      :'data.frame':	1 obs. of  1 variable:
##   .. .. .. ..$ labelDescription: chr "Type of read; see ?scanBam"
##   .. .. ..@ data             :'data.frame':	64206 obs. of  1 variable:
##   .. .. .. ..$ flag: int [1:64206] 16 0 0 0 0 0 0 0 16 16 ...
##   .. .. ..@ dimLabels        : chr [1:2] "readName" "alignColumn"
##   .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
##   .. .. .. .. ..@ .Data:List of 1
##   .. .. .. .. .. ..$ : int [1:3] 1 1 0
##   ..@ quality     :Formal class 'FastqQuality' [package "ShortRead"] with 1 slots
##   .. .. ..@ quality:Formal class 'BStringSet' [package "Biostrings"] with 5 slots
##   .. .. .. .. ..@ pool           :Formal class 'SharedRaw_Pool' [package "XVector"] with 2 slots
##   .. .. .. .. .. .. ..@ xp_list                    :List of 1
##   .. .. .. .. .. .. .. ..$ :<externalptr> 
##   .. .. .. .. .. .. ..@ .link_to_cached_object_list:List of 1
##   .. .. .. .. .. .. .. ..$ :<environment: 0x10a1e7830> 
##   .. .. .. .. ..@ ranges         :Formal class 'GroupedIRanges' [package "XVector"] with 7 slots
##   .. .. .. .. .. .. ..@ group          : int [1:64206] 1 1 1 1 1 1 1 1 1 1 ...
##   .. .. .. .. .. .. ..@ start          : int [1:64206] 1 37 73 109 145 181 217 253 289 325 ...
##   .. .. .. .. .. .. ..@ width          : int [1:64206] 36 36 36 36 36 36 36 36 36 36 ...
##   .. .. .. .. .. .. ..@ NAMES          : NULL
##   .. .. .. .. .. .. ..@ elementType    : chr "integer"
##   .. .. .. .. .. .. ..@ elementMetadata: NULL
##   .. .. .. .. .. .. ..@ metadata       : list()
##   .. .. .. .. ..@ elementType    : chr "BString"
##   .. .. .. .. ..@ elementMetadata: NULL
##   .. .. .. .. ..@ metadata       : list()
##   ..@ sread       :Formal class 'DNAStringSet' [package "Biostrings"] with 5 slots
##   .. .. ..@ pool           :Formal class 'SharedRaw_Pool' [package "XVector"] with 2 slots
##   .. .. .. .. ..@ xp_list                    :List of 1
##   .. .. .. .. .. ..$ :<externalptr> 
##   .. .. .. .. ..@ .link_to_cached_object_list:List of 1
##   .. .. .. .. .. ..$ :<environment: 0x10a1e7830> 
##   .. .. ..@ ranges         :Formal class 'GroupedIRanges' [package "XVector"] with 7 slots
##   .. .. .. .. ..@ group          : int [1:64206] 1 1 1 1 1 1 1 1 1 1 ...
##   .. .. .. .. ..@ start          : int [1:64206] 1 37 73 109 145 181 217 253 289 325 ...
##   .. .. .. .. ..@ width          : int [1:64206] 36 36 36 36 36 36 36 36 36 36 ...
##   .. .. .. .. ..@ NAMES          : NULL
##   .. .. .. .. ..@ elementType    : chr "integer"
##   .. .. .. .. ..@ elementMetadata: NULL
##   .. .. .. .. ..@ metadata       : list()
##   .. .. ..@ elementType    : chr "DNAString"
##   .. .. ..@ elementMetadata: NULL
##   .. .. ..@ metadata       : list()
##   ..@ id          :Formal class 'BStringSet' [package "Biostrings"] with 5 slots
##   .. .. ..@ pool           :Formal class 'SharedRaw_Pool' [package "XVector"] with 2 slots
##   .. .. .. .. ..@ xp_list                    :List of 1
##   .. .. .. .. .. ..$ :<externalptr> 
##   .. .. .. .. ..@ .link_to_cached_object_list:List of 1
##   .. .. .. .. .. ..$ :<environment: 0x10a1e7830> 
##   .. .. ..@ ranges         :Formal class 'GroupedIRanges' [package "XVector"] with 7 slots
##   .. .. .. .. ..@ group          : int [1:64206] 1 1 1 1 1 1 1 1 1 1 ...
##   .. .. .. .. ..@ start          : int [1:64206] 1 29 57 85 114 144 173 202 230 260 ...
##   .. .. .. .. ..@ width          : int [1:64206] 28 28 28 29 30 29 29 28 30 30 ...
##   .. .. .. .. ..@ NAMES          : NULL
##   .. .. .. .. ..@ elementType    : chr "integer"
##   .. .. .. .. ..@ elementMetadata: NULL
##   .. .. .. .. ..@ metadata       : list()
##   .. .. ..@ elementType    : chr "BString"
##   .. .. ..@ elementMetadata: NULL
##   .. .. ..@ metadata       : list()
```


To index a bam file

```r
indexFile <- indexBam("drosophilaMelanogasterSubset.bam")
```


To filter reads with at most 2N's

```r
nFilt <- nFilter(2)
```


To filter reads which have been aligned to the reference genome

```r
chrFilt <- chromosomeFilter(regex = "chr")
```


To create a new filter which checks that both conditions are satisfied i.e. a read has 2 or fewer N's and is aligned to the reference genome

```r
filt <- compose(nFilt, chrFilt)
```


Apply the filter to extract the relevant subset of data

```r
bamDatFiltered <- bamDat[filt(bamDat)]
```


To examine the filter bam file

```r
str(bamDatFiltered)
```

```
## Formal class 'AlignedRead' [package "ShortRead"] with 8 slots
##   ..@ chromosome  : Factor w/ 7 levels "chrM","chr2L",..: 1 1 1 1 1 1 1 1 1 1 ...
##   ..@ position    : int [1:56883] 548 1497 1506 1528 1540 1552 1552 1555 1559 1566 ...
##   ..@ strand      : Factor w/ 3 levels "+","-","*": 2 1 1 1 1 1 1 1 2 2 ...
##   ..@ alignQuality:Formal class 'NumericQuality' [package "ShortRead"] with 1 slots
##   .. .. ..@ quality: int [1:56883] 132 132 127 130 130 122 132 132 132 132 ...
##   ..@ alignData   :Formal class 'AlignedDataFrame' [package "ShortRead"] with 4 slots
##   .. .. ..@ varMetadata      :'data.frame':	1 obs. of  1 variable:
##   .. .. .. ..$ labelDescription: chr "Type of read; see ?scanBam"
##   .. .. ..@ data             :'data.frame':	56883 obs. of  1 variable:
##   .. .. .. ..$ flag: int [1:56883] 16 0 0 0 0 0 0 0 16 16 ...
##   .. .. ..@ dimLabels        : chr [1:2] "readName" "alignColumn"
##   .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
##   .. .. .. .. ..@ .Data:List of 1
##   .. .. .. .. .. ..$ : int [1:3] 1 1 0
##   ..@ quality     :Formal class 'FastqQuality' [package "ShortRead"] with 1 slots
##   .. .. ..@ quality:Formal class 'BStringSet' [package "Biostrings"] with 5 slots
##   .. .. .. .. ..@ pool           :Formal class 'SharedRaw_Pool' [package "XVector"] with 2 slots
##   .. .. .. .. .. .. ..@ xp_list                    :List of 1
##   .. .. .. .. .. .. .. ..$ :<externalptr> 
##   .. .. .. .. .. .. ..@ .link_to_cached_object_list:List of 1
##   .. .. .. .. .. .. .. ..$ :<environment: 0x10a1e7830> 
##   .. .. .. .. ..@ ranges         :Formal class 'GroupedIRanges' [package "XVector"] with 7 slots
##   .. .. .. .. .. .. ..@ group          : int [1:56883] 1 1 1 1 1 1 1 1 1 1 ...
##   .. .. .. .. .. .. ..@ start          : int [1:56883] 1 37 73 109 145 181 217 253 289 325 ...
##   .. .. .. .. .. .. ..@ width          : int [1:56883] 36 36 36 36 36 36 36 36 36 36 ...
##   .. .. .. .. .. .. ..@ NAMES          : NULL
##   .. .. .. .. .. .. ..@ elementType    : chr "integer"
##   .. .. .. .. .. .. ..@ elementMetadata: NULL
##   .. .. .. .. .. .. ..@ metadata       : list()
##   .. .. .. .. ..@ elementType    : chr "BString"
##   .. .. .. .. ..@ elementMetadata: NULL
##   .. .. .. .. ..@ metadata       : list()
##   ..@ sread       :Formal class 'DNAStringSet' [package "Biostrings"] with 5 slots
##   .. .. ..@ pool           :Formal class 'SharedRaw_Pool' [package "XVector"] with 2 slots
##   .. .. .. .. ..@ xp_list                    :List of 1
##   .. .. .. .. .. ..$ :<externalptr> 
##   .. .. .. .. ..@ .link_to_cached_object_list:List of 1
##   .. .. .. .. .. ..$ :<environment: 0x10a1e7830> 
##   .. .. ..@ ranges         :Formal class 'GroupedIRanges' [package "XVector"] with 7 slots
##   .. .. .. .. ..@ group          : int [1:56883] 1 1 1 1 1 1 1 1 1 1 ...
##   .. .. .. .. ..@ start          : int [1:56883] 1 37 73 109 145 181 217 253 289 325 ...
##   .. .. .. .. ..@ width          : int [1:56883] 36 36 36 36 36 36 36 36 36 36 ...
##   .. .. .. .. ..@ NAMES          : NULL
##   .. .. .. .. ..@ elementType    : chr "integer"
##   .. .. .. .. ..@ elementMetadata: NULL
##   .. .. .. .. ..@ metadata       : list()
##   .. .. ..@ elementType    : chr "DNAString"
##   .. .. ..@ elementMetadata: NULL
##   .. .. ..@ metadata       : list()
##   ..@ id          :Formal class 'BStringSet' [package "Biostrings"] with 5 slots
##   .. .. ..@ pool           :Formal class 'SharedRaw_Pool' [package "XVector"] with 2 slots
##   .. .. .. .. ..@ xp_list                    :List of 1
##   .. .. .. .. .. ..$ :<externalptr> 
##   .. .. .. .. ..@ .link_to_cached_object_list:List of 1
##   .. .. .. .. .. ..$ :<environment: 0x10a1e7830> 
##   .. .. ..@ ranges         :Formal class 'GroupedIRanges' [package "XVector"] with 7 slots
##   .. .. .. .. ..@ group          : int [1:56883] 1 1 1 1 1 1 1 1 1 1 ...
##   .. .. .. .. ..@ start          : int [1:56883] 1 29 57 85 114 144 173 202 230 260 ...
##   .. .. .. .. ..@ width          : int [1:56883] 28 28 28 29 30 29 29 28 30 30 ...
##   .. .. .. .. ..@ NAMES          : NULL
##   .. .. .. .. ..@ elementType    : chr "integer"
##   .. .. .. .. ..@ elementMetadata: NULL
##   .. .. .. .. ..@ metadata       : list()
##   .. .. ..@ elementType    : chr "BString"
##   .. .. ..@ elementMetadata: NULL
##   .. .. ..@ metadata       : list()
```


To look at which chromosomes are present and read IDs

```r
levels(chromosome(bamDatFiltered))
```

```
## [1] "chrM"  "chr2L" "chrX"  "chr3L" "chr4"  "chr2R" "chr3R"
```

```r
id(bamDatFiltered)[1:10]
```

```
##   A BStringSet instance of length 10
##      width seq
##  [1]    28 HWI-EAS225_90320:3:1:141:680
##  [2]    28 HWI-EAS225_90320:3:1:660:332
##  [3]    28 HWI-EAS225_90320:3:1:164:226
##  [4]    29 HWI-EAS225_90320:3:1:1088:176
##  [5]    30 HWI-EAS225_90320:3:1:1038:1484
##  [6]    29 HWI-EAS225_90320:3:1:850:1742
##  [7]    29 HWI-EAS225_90320:3:1:1319:586
##  [8]    28 HWI-EAS225_90320:3:1:103:631
##  [9]    30 HWI-EAS225_90320:3:1:1353:1498
## [10]    30 HWI-EAS225_90320:3:1:1092:1016
```


To view the DNA sequence for the read

```r
sread(bamDatFiltered)[1:10]
```

```
##   A DNAStringSet instance of length 10
##      width seq
##  [1]    36 GGAAATCAAAAATGGAAAGGAGCGGCTCCACTTTTT
##  [2]    36 AAATCATAAAGATATTGGAACTTTATATTTTATTTT
##  [3]    36 AGATATTGGAACTTTATATTTTATTTTTGGAGCTTG
##  [4]    36 ATTTTTGGAGCTTGAGCTGGAATAGTTGGAACATCT
##  [5]    36 TGAGCTGGAATAGTTGGAACATCTTTAAGAATTTTA
##  [6]    36 GTTGGAACATCTTTAAGAATTTTAATTAGAGCTGAA
##  [7]    36 GTTGGAACATCTTTAAGAATTTTAATTCGAGCTGAA
##  [8]    36 GGAACATCTTTAAGAATTTTAATTCGAGCTGAATTA
##  [9]    36 GTCCTAATTCAGCTCGAATTAAAATTCTTAAAGATG
## [10]    36 CCAGGATGTCCTAATTCAGCTCGAATTAAAATTCTT
```


To look at the base qualities

```r
quality(bamDatFiltered)[1:10]
```

```
## class: FastqQuality
## quality:
##   A BStringSet instance of length 10
##      width seq
##  [1]    36 BACBCCABBBBBA>8@B@@==>;5-9A<;=7A@@B@
##  [2]    36 BCBBABBA@@B;B>AB@@<>:AAA9?>?A@A<?A@@
##  [3]    36 @?8AB>A?=)A=@*8>6/@3>A)/@4>?BA'(-1B=
##  [4]    36 BBCACCA@-4ABC62?*;A?BBA?B@.8B9?33;+=
##  [5]    36 ?5@A4::@@55;;89<'6?A8@A=4@=>54>76);A
##  [6]    36 A8=B;462>;7BCBAA>1;=</?BA94%<:?(7@9=
##  [7]    36 BBB>ABB@@BBBBCBCC@7ABBBAABB@B?AAA@=@
##  [8]    36 =B=ACCBBC8ACCCBBBCCCBB=CAB9=BBBB@2?:
##  [9]    36 @7ABBBBABBB?BAAB=@CBB;7ABAABBA?@;2=A
## [10]    36 BB>B?:ABABBBBABCB@@@BB@:@BA;>@;@B?AB
```


To find out the starting (left most) position where the read was aligned on the chromosome

```r
position(bamDatFiltered)[1:10]
```

```
##  [1]  548 1497 1506 1528 1540 1552 1552 1555 1559 1566
```


To see which strand (forward/+ or reverse/-) the read was aligned to 

```r
strand(bamDatFiltered)[1:10]
```

```
##  [1] - + + + + + + + - -
## Levels: + - *
```


Load genomic data

```r
library(BSgenome.Dmelanogaster.UCSC.dm3)
```

```
## Loading required package: BSgenome
```


To get the lengths of the chromosomes in the genome

```r
(chrSizes <- seqlengths(Dmelanogaster))
```

```
##     chr2L     chr2R     chr3L     chr3R      chr4      chrX      chrU 
##  23011544  21146708  24543557  27905053   1351857  22422827  10049037 
##      chrM  chr2LHet  chr2RHet  chr3LHet  chr3RHet   chrXHet   chrYHet 
##     19517    368872   3288761   2555491   2517507    204112    347038 
## chrUextra 
##  29004656
```


Load library and specifiy Drosophila database

```r
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "dmelanogaster_gene_ensembl")
```



To define a set of fields we are interested for each gene

```r
annotation.fields <- c("ensembl_gene_id", "strand", "chromosome_name", "start_position", 
    "end_position")
```


To restrict to chr2L

```r
gene.annotation <- getBM(annotation.fields, mart = ensembl, filters = "chromosome_name", 
    values = c("2L"))
str(gene.annotation)
```

```
## 'data.frame':	2986 obs. of  5 variables:
##  $ ensembl_gene_id: chr  "FBgn0031208" "FBgn0002121" "FBgn0031209" "FBgn0263584" ...
##  $ strand         : int  1 -1 -1 1 -1 1 1 1 -1 1 ...
##  $ chromosome_name: chr  "2L" "2L" "2L" "2L" ...
##  $ start_position : int  7529 9839 21823 21952 25402 66584 71757 76348 82421 94739 ...
##  $ end_position   : int  9484 21376 25155 24237 65404 71390 76211 77783 87387 102086 ...
```


To check we only downloaded annotations from chr2L

```r
levels(as.factor(gene.annotation$chromosome))
```

```
## [1] "2L"
```


Add chr to the annotation data

```r
gene.annotation$chromosome <- paste("chr", gene.annotation$chromosome_name, 
    sep = "")
levels(as.factor(gene.annotation$chromosome))
```

```
## [1] "chr2L"
```


To store the gene annotation information in an IRanges object

```r
gene.range <- RangedData(IRanges(start = gene.annotation$start_position, end = gene.annotation$end_position), 
    space = gene.annotation$chromosome, strand = gene.annotation$strand, gene = gene.annotation$ensembl_gene_id, 
    universe = "Dm3")
show(gene.range)
```

```
## RangedData with 2986 rows and 2 value columns across 1 space
##         space               ranges   |    strand        gene
##      <factor>            <IRanges>   | <integer> <character>
## 1       chr2L       [ 7529,  9484]   |         1 FBgn0031208
## 2       chr2L       [ 9839, 21376]   |        -1 FBgn0002121
## 3       chr2L       [21823, 25155]   |        -1 FBgn0031209
## 4       chr2L       [21952, 24237]   |         1 FBgn0263584
## 5       chr2L       [25402, 65404]   |        -1 FBgn0051973
## 6       chr2L       [66584, 71390]   |         1 FBgn0067779
## 7       chr2L       [71757, 76211]   |         1 FBgn0031213
## 8       chr2L       [76348, 77783]   |         1 FBgn0031214
## 9       chr2L       [82421, 87387]   |        -1 FBgn0002931
## ...       ...                  ... ...       ...         ...
## 2978    chr2L [22690251, 22691008]   |        -1 FBgn0058439
## 2979    chr2L [22735486, 22736297]   |        -1 FBgn0262947
## 2980    chr2L [22736952, 22747273]   |         1 FBgn0041004
## 2981    chr2L [22811944, 22834955]   |         1 FBgn0002566
## 2982    chr2L [22841770, 22843208]   |        -1 FBgn0058005
## 2983    chr2L [22874534, 22885080]   |         1 FBgn0000384
## 2984    chr2L [22892306, 22918647]   |        -1 FBgn0250907
## 2985    chr2L [22959606, 22961179]   |        -1 FBgn0086683
## 2986    chr2L [22961737, 22963456]   |         1 FBgn0262887
```


To find out how many bases cover each position in every chromosome

```r
(cover <- coverage(bamDatFiltered, width = chrSizes))
```

```
## RleList of length 7
## $chrM
## integer-Rle of length 19517 with 953 runs
##   Lengths:  547   36  913    9   22    5 ...    3    5   13  258   36 5793
##   Values :    0    1    0    1    2    3 ...    4    3    2    0    1    0
## 
## $chr2L
## integer-Rle of length 23011544 with 17850 runs
##   Lengths:  6777    36  2316    36  1621 ...   107    36   499    36 50474
##   Values :     0     1     0     1     0 ...     0     1     0     1     0
## 
## $chrX
## integer-Rle of length 22422827 with 16522 runs
##   Lengths:  18996     36  12225     36 ...     36    130     36   6180
##   Values :      0      1      0      1 ...      1      0      1      0
## 
## $chr3L
## integer-Rle of length 24543557 with 17396 runs
##   Lengths: 135455     36   6783     23 ...     36  82251     36  12469
##   Values :      0      1      0      1 ...      1      0      1      0
## 
## $chr4
## integer-Rle of length 1351857 with 1255 runs
##   Lengths:  59510     36   2019     36 ...     36    267     36 118808
##   Values :      0      1      0      1 ...      1      0      1      0
## 
## ...
## <2 more elements>
```


To aggregate the coverage for each gene

```r
gene.coverage <- aggregate(cover[match(names(gene.range), names(cover))], ranges(gene.range), 
    sum)
```


To find the number of reads covering each gene

```r
gene.coverage <- ceiling(gene.coverage/unique(width(bamDat)))
gene.coverage
```

```
## NumericList of length 1
## [["chr2L"]] 1 47 0 0 1 6 0 0 8 11 1 1 58 ... 17 16 1 0 0 0 15 4 0 1 0 6 0
```


To restrict to chr2L

```r
length(gene.coverage[["chr2L"]])
```

```
## [1] 2986
```

```r
length(ranges(gene.range)$chr2L)
```

```
## [1] 2986
```


To build a count table and store it in a data frame

```r
countTable <- data.frame(chromosome = gene.range$space, gene_start = start(gene.range$ranges), 
    gene_end = end(gene.range$ranges), strand = gene.range$strand, gene = gene.range$gene, 
    count = as.vector(gene.coverage[["chr2L"]]))
dim(countTable)
```

```
## [1] 2986    6
```

```r
head(countTable)
```

```
##   chromosome gene_start gene_end strand        gene count
## 1      chr2L       7529     9484      1 FBgn0031208     1
## 2      chr2L       9839    21376     -1 FBgn0002121    47
## 3      chr2L      21823    25155     -1 FBgn0031209     0
## 4      chr2L      21952    24237      1 FBgn0263584     0
## 5      chr2L      25402    65404     -1 FBgn0051973     1
## 6      chr2L      66584    71390      1 FBgn0067779     6
```


To use RPKM (reads per kilobase per million mapped reads)

```r
countTable <- data.frame(chromosome = gene.range$space, gene_start = start(gene.range$ranges), 
    gene_end = end(gene.range$ranges), strand = gene.range$strand, gene = gene.range$gene, 
    count = as.vector(gene.coverage[["chr2L"]]), RPKM = (as.vector(gene.coverage[["chr2L"]])/(end(gene.range$ranges) - 
        start(gene.range$ranges))) * (1e+09/length(bamDat)))
head(countTable)
```

```
##   chromosome gene_start gene_end strand        gene count    RPKM
## 1      chr2L       7529     9484      1 FBgn0031208     1  7.9667
## 2      chr2L       9839    21376     -1 FBgn0002121    47 63.4497
## 3      chr2L      21823    25155     -1 FBgn0031209     0  0.0000
## 4      chr2L      21952    24237      1 FBgn0263584     0  0.0000
## 5      chr2L      25402    65404     -1 FBgn0051973     1  0.3894
## 6      chr2L      66584    71390      1 FBgn0067779     6 19.4443
```


