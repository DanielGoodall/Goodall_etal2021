# purpose: find the 23 MatS termination sites and overlap that intot the super figure
# date created: 08/05/21
# used with other prgrams: Not needed


library(Biostrings)

mats <- c('GTGACAATGTCAC')
DNAString(mats)
mats_dnaset <- DNAStringSet(mats)
names(mats_dnaset) <- 'MatS'
mats_dnaset

# setwd
file_in_dir <- file.choose()
dir <- dirname(file_in_dir)
setwd(dir)
getwd()

# write to fastq for BT2 to handle
#writeXStringSet(mats_dnaset, 'mats.fastq', format = 'fastq') 

# Part1
#                                     Using Biostrings, IRanges and GRanges
#____________________________________________________________________________________________________

library(IRanges)
library(GenomicRanges)

Ecoli <- readDNAStringSet(file.choose())
Ecoli


# MG1655
MG_dna <- DNAString(Ecoli$MG1655)
mats2 <- c('GTGACNNNGTCAC')


match2_pos <- matchPattern(mats2, MG_dna, max.mismatch = 0, fixed=FALSE)
match2_pos
length(match2_pos)

MG_ir1 <- IRanges(match2_pos)
names(MG_ir1) <- paste('MatS', 1:length(MG_ir1), sep = '')
MG_df <- data.frame(MG_ir1)


# Part 
#________________________________Bowtie2_________________________________
library(Rbowtie2)

# bowtie2 alignment with SAM output
bowtie2(bt2Index ='Ecoli_genomes', 'MatS_Ecoli.sam',
        seq1 = 'MatS.fastq',
        seq2 = NULL,
        "-D 15 -R 3 -N 1 -L 20 -i S,1,0.50 -M 2", "--all", "--reorder", 
        overwrite = TRUE)


#________________________________Samtools_________________________________
library(Rsamtools)

# processing to a dataframe
scan <- scanBam(asBam('MatS_Ecoli.sam'))
bam <- scanBam('MatS_Ecoli.bam')
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
#csv <- write.csv(bam_df, 'MatS_Ecoli.csv')




#_____________________________Clean CSV___________________________
mats_df <- data.frame(read.csv('MatS_Ecoli.csv'))

# drop rows we do not want and rename
mats_df <- mats_df[c(2,3,4,6)]
colnames(mats_df) <- c('strain','name', 'start', 'strand')

#sanity
mats_df



# subset the large df into each strain
# then later test with just the 5 phylogroups
#A
MG1655 <- mats_df[which(mats_df$strain == 'MG1655'),]
BW2952 <- mats_df[which(mats_df$strain == 'BW2952'),]
REL606 <- mats_df[which(mats_df$strain == 'REL606'),]
#B1
APEC078 <- mats_df[which(mats_df$strain == 'APECO78'),]
IAI1 <- mats_df[which(mats_df$strain == 'IAI1'),]
E11368 <- mats_df[which(mats_df$strain == '11368'),]
#B2
S88 <- mats_df[which(mats_df$strain == 'S88'),]
UTI89 <- mats_df[which(mats_df$strain == 'UTI89'),]
E2348 <- mats_df[which(mats_df$strain == 'E2348/69'),]
#D
IAI39 <- mats_df[which(mats_df$strain == 'IAI39'),]
SMS35 <- mats_df[which(mats_df$strain == 'SMS-3-5'),]
UMN026 <- mats_df[which(mats_df$strain == 'UMN026'),]
CE10 <- mats_df[which(mats_df$strain == 'CE10'),]
D042 <- mats_df[which(mats_df$strain == '042'),]
#E
TW14359 <- mats_df[which(mats_df$strain == 'TW14359'),]
Sakai <- mats_df[which(mats_df$strain == 'Sakai'),]
EDL933 <- mats_df[which(mats_df$strain == 'EDL933'),]


# reset indexes
#A
rownames(MG1655) <- NULL
rownames(BW2952) <- NULL
rownames(REL606) <- NULL

#B1
rownames(APEC078) <- NULL
rownames(IAI1) <- NULL
rownames(E11368) <- NULL

#B2
rownames(S88) <- NULL
rownames(UTI89) <- NULL
rownames(E2348) <- NULL

#D
rownames(IAI39) <- NULL
rownames(SMS35) <- NULL
rownames(UMN026) <- NULL
rownames(CE10) <- NULL
rownames(D042) <- NULL

#E
rownames(TW14359) <- NULL
rownames(Sakai) <- NULL
rownames(EDL933) <- NULL


# SANIITY
#A
MG1655
BW2952
REL606
#B1
APEC078
IAI1
E11368
#B2
S88
UTI89
E2348
#D
IAI39
SMS35
UMN026
CE10
D042
#E
TW14359
Sakai
EDL933

# Drop the Strain column and add the end column
# without specifying [-1] for every df
list_mats <- list(
  MG1655,
  REL606,
  APEC078,
  IAI1,
  E11368,
  S88,
  UTI89,
  E2348,
  IAI39,
  SMS35,
  UMN026,
  CE10,
  D042,
  TW14359,
  Sakai,
  EDL933)
# set names of lists
names(list_mats) <- c(
  'MG1655',
  'REL606',
  'APEC078',
  'IAI1',
  'E11368',
  'S88',
  'UTI89',
  'E2348',
  'IAI39',
  'SMS35',
  'UMN026',
  'CE10',
  'D042',
  'TW14359',
  'Sakai',
  'EDL933')

list_mats


# drop the duplicate start positions
list_mats_clean <- lapply(list_mats, function(x) x[!duplicated(x$start),])

# reset the indexes using lapply!
#https://stackoverflow.com/questions/16938064/r-lapply-on-list-of-dataframes-resetting-rownames
list_mats_clean2 <- lapply(list_mats_clean, function(x) {
  rownames(x) <- NULL;x
})

list_mats_clean2[1]

# Part 3 
#_______________________BT2 + Biostrings_____________________

# BT2 MG1655  dataframe
x <- list_mats_clean2$MG1655[17,]$start
list_mats_clean2$MG1655[18,]$start

# Biostrings df
MG_df


