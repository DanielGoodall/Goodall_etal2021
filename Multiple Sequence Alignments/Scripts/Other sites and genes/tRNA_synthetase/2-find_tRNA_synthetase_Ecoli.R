# 1. Read in FASTA file created with 1st program
# 2. Turn into FASTQ format ready for bowtie2
# 3. align those genes against the 17 E.coli genomes (only 5 needed but may as well)
#     this index already exists in the directory as Ecoli_genomes
# 4. Export to SAM file using Rsamtools
# 5. clean the data and put in a spreadsheet as CSV
# 6. Parse through genoplotR to map their positions then compare



#_______________________________________________________________________________

#BiocManager::install("Rbowtie2")
library(Rbowtie2)
#BiocManager::install("Rsamtools")
library(Rsamtools)
#install.packages("genoPlotR")
library(genoPlotR)
library(Biostrings)
library(seqinr)
library(tidyverse)
library(stringr) # for regular expressionn wrapping

# setwd
file_in_dir <- file.choose()
dir <- dirname(file_in_dir)
setwd(dir)
getwd()

# read in the tRNA synthetase genes FASTA
tRNA_fas <- read.fasta('tRNA_synthetase_genes.fasta', seqtype = 'DNA', as.string = TRUE)
tRNA_fas

tRNA_fas2 <- readDNAStringSet('tRNA_synthetase_genes.fasta')
tRNA_fas2

# Write FASTQ format for BT2
fastq <- writeXStringSet(tRNA_fas2, filepath = 'tRNA_synthetase_genes.fastq', format = 'FASTQ')
#fastq


# specify the working dir of index location
# Ecoli_genomes
file_in_dir <- file.choose()
dir <- dirname(file_in_dir)
setwd(dir)
getwd()

# bowtie2 alignment with SAM output
bowtie2(bt2Index ='Ecoli_genomes', 'tRNA_Ecoli.sam',
        seq1 = 'tRNA_synthetase_genes.fastq',
        seq2 = NULL,
        "-D 15 -R 3 -N 1 -L 20 -i S,1,0.50 -M 2", "--all", "--reorder",
        overwrite = TRUE)


##########################
#     SAMTOOLS 
#     sam to csv
##########################
# https://gist.github.com/davetang/6460320
# https://bioconductor.org/packages/release/bioc/vignettes/Rsamtools/inst/doc/Rsamtools-Overview.pdf

library(Rsamtools)

# processing to a dataframe
scan <- scanBam(asBam('tRNA_Ecoli.sam'))
bam <- scanBam('tRNA_Ecoli.bam')
bam_field <- c("rname","qname","pos","qwidth","strand","seq","flag")
bam_field

names(bam[[1]])

#function for collapsing the list of lists into a single list
#as per the Rsamtools vignette
.unlist <- function (x){
  ## do.call(c, ...) coerces factor to integer, which is undesired
  x1 <- x[[1L]]
  if (is.factor(x1)){
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}

#go through each BAM field and unlist
list <- lapply(bam_field, function(y) .unlist(lapply(bam, "[[", y)))

#store as data frame
bam_df <- do.call("DataFrame", list)
names(bam_df) <- bam_field

#Sanity
bam_df

# FINISH AND WRITE TO SPREADSHEET
csv <- write.csv(bam_df, 'tRNA_Ecoli.csv')
csv










