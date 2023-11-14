# 2022 Benjamin J Perry - Attribution-NonCommercial-ShareAlike 4.0 International
# (CC BY-NC-SA 4.0)
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz
# Status: Dev.

library(tidyverse)
library(Rsubread)
library(DESeq2)

list.files()

makeCounts <- function(refIndex,read1,read2,bamOut,refGTF,threads) {
  align.stat <- subjunc(
    index = refIndex,
    readfile1 = read1,
    readfile2 = read2,
    output_file = bamOut,
    input_format = "gzFASTQ",
    output_format = "BAM",
    maxFragLength = 800,
    nthreads = threads,
    sortReadsByCoordinates = T,
    detectSV = T,
    reportAllJunctions = T,
    
  )
  
  counts.stat <- featureCounts(
    files = bamOut,
    annot.ext = refGTF,
    isGTFAnnotationFile = T,
    minMQS = 10,
    strandSpecific = 1,
    juncCounts = T,
    isPairedEnd = T,
    requireBothEndsMapped = T,
    checkFragLength = T,
    maxFragLength = 800,
    allowMultiOverlap = T,
    nthreads = threads
  )
  
  return(counts.stat)
}


### Make the feature counts table for the RNASeq data ###
# Build reference for alignment
GenomeIndex <- "ref/reference.fasta" #Path to reference genome
buildindex(basename = "ref/Reference",
           reference = "ref/Reference.fna")
GenomeGTF <- "ref/Reference.RSEM.gtf" #Rath to ref GTF format


# Process the experimental metadata file
RNASEQ <- read_csv(
  file = "experiment.metadata.csv",
  col_names = T,
  col_types = "cfcc",
  trim_ws = T
)

# Create file paths based on the samples in the metadata file for the output bam files
RNASEQ$BAM <- paste(RNASEQ$SamplID, "subjunc", "sort", "bam", sep = ".")


# Compute feature counts object for  samples
countsObj <- makeCounts(
  refIndex = GenomeIndex,
  read1 = RNASEQ$R1,
  read2 = RNASEQ$R2,
  bamOut = RNASEQ$BAM,
  refGTF = GenomeGTF,
  threads = 24
  )


# Print our featureCounts stats
StatsTable <- as.data.frame(x = countsObj$stat)
write_tsv(StatsTable, path = "featureCounts.summary.txt")


### Prepare DESeqDataSet object ###
# Tidy the count data
countData <- as.data.frame(countsObj$counts)
colnames(countData) <- RNASEQ$SampleID
countData$ID <- rownames(countData)
countData <- countData %>% select(ID, everything())
condition <- factor(RNASEQ$condition)


# Merge the dds Object
ddsObj <- DESeqDataSetFromMatrix(countData, DataFrame(condition), ~ condition, tidy = T)
head(ddsObj@assays@data@listData[["counts"]])


### Start DESeq2 DE Analysis ###
# Filter out features with no transcripts
keep <- rowSums(counts(ddsObj)) > 1
ddsObjFilter <- ddsObj[keep,]
nrow(ddsObj)
nrow(ddsObjFilter)



### Exploratory Analysis to Check Nothing Looks Unusual ###
# Make a regular log normal counts table for cursory comparison
rlogddsObj <- rlog(ddsR7AFilter, blind = F)
head(assay(rlogddsObj), 10)
rlogddsObjsampleDist <- dist(t(assay(rlogddsObj)))
# Show the euclidean distance between samples at all genes after normalization
rlogddsObjsampleDist



### Update DESeq2 commands with the appropriate options ###
ddsObjFilter <- DESeq(ddsObjFilter)
# Compute DE expression results object
resddsObjFilter <- results(
  ddsObjFilter,
  contrast = c("condition", "Treatment1", "Treatment2"), #Update Conditions
  cooksCutoff = F,
  test = "Wald",
  independentFiltering = F,
  pAdjustMethod = 'BH'
)



### Prepare Files to Export ###
resOut <- as.data.frame(resddsObjFilter)
resOut$ID <- row.names(resOut)
resOut <- resOut %>% select(ID, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)
normCounts <- as.data.frame(counts(ddsObjFilter, normalized=T))
normCounts$ID <- row.names(normCounts)

### Update column names in table for your samples ###
normCounts <- normCounts %>% select(ID,
                                    normalize.counts.SAMPLE1 = sample1,
                                    normalize.counts.SAMPLE2 = sample2,
                                    normalize.counts.SAMPLE3 = sample3,
                                    normalize.counts.SAMPLE4 = sample4)


DEout <- left_join(x = countData, y = normCounts, by = "ID")
DEout <- left_join(x = DEout, y = resOut, by = "ID")

