# 1. Parse a DAT file containing the ribosomes genes using RegEx
# 2. Create a FASTA file with these genes




#_______________________________________________________________________________


library(Biostrings)
library(seqinr)
library(tidyverse)
library(stringr) # for regular expressionn wrapping
library(ShortRead) # for writeFASTQ

# read in tRNA synthetase DAT file as text
# use read.delim as it is effectively a tabbed data file
fname <- toString(read.delim(file.choose(), header = TRUE, sep= "\t")[[1]])
fname



# remove the cahracter we do not want:
# \\par
dropped1 <- gsub('[{}*0-9"d"]', '', fname)
dropped1




# split on the empty space using tidyverse and stringr
list0 <- str_split(dropped1, ' ')
list1 <- list0[[1]]
list1

# drop the first 10 elements as they are not the tRNA synthetase seqs
list2 <- list1[-c(1:4)]
list2

# drop all the blank ''
# no RegEx needed
list3 <- list2[list2 != '']
list3

# delete all the recurring ,,,,,
#
list4 <- gsub(',,,,,|,,,,', '', list3)
list4

# next try to extract all the names into a list
## these begin with 2-4 letters and always lowercase then uppercase 
### e.g cycS or glyS or lysU or thS
names <- gsub(',[AGTC]|[ATGC]', '', list4)
names

# then try to extract all seqs into a list
seqs <- gsub('[a-z]|[A-Z],', '', list4)
seqs

# manually replace all the ATGC that got deleted from names
## needs to be cleaned up by not deleting Capitals before the ,
###  but it works
# names df
names_df <- data.frame(names)
names_df[1,] <- 'rplA' # add back in the A
names_df[3,] <- 'rplC' # add back in the C
names_df[18,] <- 'rplT' # add back in the T 
names_df[24,] <- 'rpmA' # add back in the A
names_df[26,] <- 'rpmC' # add back in the C
names_df[30,] <- 'rpmG' # add back in the G
names_df[34,] <- 'rpsA'# add back the A
names_df[36,] <- 'rpsC' # add back the C
names_df[40,] <- 'rpsG' # add back the G

names_df # sanity

colnames(names_df) <- 'ribosome_names'
# sanity
names_df

# seqs df
seqs_df <- data.frame(seqs)
colnames(seqs_df) <- 'ribosome_seqs'
#sanity
seqs_df

# add together into one df
ribo_df <- cbind(names_df, seqs_df)
ribo_df

# Biostrings objects
ribo_dnastringset <- DNAStringSet(ribo_df$ribosome_seqs)
names(ribo_dnastringset) <-  ribo_df$ribosome_names
as.list(ribo_dnastringset)
ribo_dnastringset #sanity

#                           Create FASTA
#____________________________________________________________

# set wd of where to put the fasta file (in bowtie2 then in ribosome_proteins)
file_in_dir <- file.choose()
dir <- dirname(file_in_dir)
setwd(dir)

write_fasta <-  writeXStringSet(ribo_dnastringset, 'ribosome_genes.fasta', format = 'fasta')
#write_fasta

# write fastq in the 2nd script











