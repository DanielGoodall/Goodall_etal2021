# find position of the pseudo ter sites and see if theres any change in position


# library
library(Biostrings)
library(Rbowtie2)
library(Rsamtools)

# specify the pseudo ter sites
terK <- DNAString('CGATTGAGAGTTGTAATGAAGTC')
terL <- DNAString('GCACTGGGTGTTGTAATGACGCA')
terY <- DNAString('TATGGGTACGTTGTAATTAGGGA')
terZ <- DNAString('TACCCGCAGGTTGTAACGAGAGC')


#turn into DNAStringSet to write to a fastq ready for BT2
pseudo_set <- DNAStringSet(list(terK, terL, terY, terZ))
names(pseudo_set) <- c('terK', 'terL', 'terY', 'terZ')
pseudo_set # sanity

#set wd
Z4_m
# write to fastq
#writeXStringSet(pseudo_set, 'pseudo_ter.fastq', format = 'fastq')

#also write larger fastq with all the ter sites in
ter <- readDNAStringSet(file.choose(), format = 'fastq')
total_ter <- c(ter, pseudo_set)
#writeXStringSet(total_ter, 'total_ter.fastq', format = 'fastq')

#______________Bowtie2__________________
#         Align pseudo ter with Ecoli
# just the 4 pseudo ter sites
bowtie2(bt2Index ='Ecoli_genomes', 'pseudo_ter_Ecoli_V2.sam',
        seq1 = 'pseudo_ter.fastq',
        seq2 = NULL,
        "-D 15 -R 3 -N 1 -L 12 -i S,1,0.50 -M 2", '--mp 3', "--all", "--reorder",
        overwrite = TRUE)

# all 14 ter sites NOT NEEDED as we already have that information and we can just rbind the dna_segs later


#_________________SAMTOOLS ______________
#               sam to csv

# https://gist.github.com/davetang/6460320
# https://bioconductor.org/packages/release/bioc/vignettes/Rsamtools/inst/doc/Rsamtools-Overview.pdf

########################################################################################
#
#                           Pseudo Ter sites Only
#
########################################################################################

# processing to a dataframe
scan <- scanBam(asBam('pseudo_ter_Ecoli_V2.sam'))
bam <- scanBam('pseudo_ter_Ecoli_V2.bam')
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
csv <- write.csv(bam_df, 'pseudo_ter_Ecoli_V2.csv')
csv



########################################################################################
#
#                           Total 14 Ter sites 
#
########################################################################################

# processing to a dataframe
scan <- scanBam(asBam('total_ter_Ecoli.sam'))
bam <- scanBam('total_ter_Ecoli.bam')
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
csv <- write.csv(bam_df, 'total_ter_Ecoli.csv')
csv