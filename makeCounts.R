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
    reportAllJunctions = T,
    
  )
  
  counts <- featureCounts(
    files = bamOut,
    annot.ext = refGTF,
    isGTFAnnotationFile = T,
    minMQS = 10,
    juncCounts = T,
    isPairedEnd = T,
    requireBothEndsMapped = T,
    checkFragLength = T,
    maxFragLength = 800,
    strandSpecific = 1,
    useMetaFeatures = T,
    nthreads = threads
  )
  
  return(counts)
}


### Make the feature counts table for the RNASeq data ###

# Build reference for alignment
genomeIndex <- "ref/ARS-UI_Ramb_v2_0"

if (file.exists(paste(genomeIndex, ".00.b.tab", sep =''))) {
  cat(paste("Using index ", genomeIndex))
} else {
  print("Building genome index...")
  buildindex(basename = genomeIndex,
             reference = "ref/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna")
}

genomeGTF <- "ref/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.gtf" #path to ref GTF format


# Process the experimental metadata file
RNASEQ <- read_csv(
  file = "NS0016/NS0016.MD.csv",
  col_names = T,
  col_types = "cfccfffdf",
  trim_ws = T
)

RNASEQ$R1 <- paste("fastq", RNASEQ$R1, sep = "/")
RNASEQ$R2 <- paste("fastq", RNASEQ$R2, sep = "/")

# Create file paths based on the samples in the metadata file for the output bam files
dir.create("bam")
RNASEQ$BAM <- paste(RNASEQ$SampleID, "subjunc", "sort", "bam", sep = ".")
RNASEQ$BAM <- paste("bam", RNASEQ$BAM, sep = "/")


# Compute feature counts object for  samples
countsObj <- makeCounts(
  refIndex = genomeIndex,
  read1 = RNASEQ$R1,
  read2 = RNASEQ$R2,
  bamOut = RNASEQ$BAM,
  refGTF = genomeGTF,
  threads = 24
)


# Print our featureCounts stats
statsTable <- as.data.frame(x = countsObj$stat)
write_csv(statsTable, path = "RNA-seq.alignment.stats.csv")


### Prepare DESeqDataSet object ###
# Tidy the count data
countData <- as.data.frame(countsObj$counts)
colnames(countData) <- RNASEQ$SampleID
countData$ID <- rownames(countData)
countData <- countData %>% select(ID, everything())
write_csv(countData, path = "RNA-seq.counts.csv")
