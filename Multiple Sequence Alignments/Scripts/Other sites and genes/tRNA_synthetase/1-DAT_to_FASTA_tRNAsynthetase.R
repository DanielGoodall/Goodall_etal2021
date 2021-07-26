# 1. Parse a DAT file containing the tRNA synthetase genes using RegEx
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
dropped1 <- gsub('[\\pard{}*0-9"d"]', '', fname)
dropped1




# split on the empty space using tidyverse and stringr
list0 <- str_split(dropped1, ' ')
list1 <- list0[[1]]
list1

# drop the first 10 elements as they are not the tRNA synthetase seqs
list2 <- list1[-c(1:10)]
list2

# drop all the blank ''
# no RegEx needed
list3 <- list2[list2 != '']
list3

# delete all the recurring ,,,,,
#
list4 <- gsub(',,,,,', '', list3)
list4

# next try to extract all the names into a list
## these begin with 2-4 letters and always lowercase then uppercase 
### e.g cycS or glyS or lysU or thS
names <- gsub(',[AGTC]|[ATGC]', '', list4)
names

# then try to extract all seqs into a list
seqs <- gsub('[a-z]|[A-Z],', '', list4)
seqs

# names df
names_df <- data.frame(names)
names_df <- data.frame(names_df[-24,]) #drop the last value which only has , in
names_df[15,] <- 'metG' # add back in the G after met
names_df[17,] <- 'heT' # add back in the T after he
colnames(names_df) <- 'tRNA_syn_genes'
# sanity
names_df

# seqs df
seqs_df <- data.frame(seqs)
seqs_df <- data.frame(seqs_df[-24,]) #drop the last value which only has , in
colnames(seqs_df) <- 'tRNA_syn_seq'
#sanity
seqs_df

# add together into one df
tRNAsyn_df <- cbind(names_df, seqs_df)
tRNAsyn_df

# Biostrings objects
tRNA_dnastringset <- DNAStringSet(tRNAsyn_df$tRNA_syn_seq)
names(tRNA_dnastringset) <-  tRNAsyn_df$tRNA_syn_genes
as.list(tRNA_dnastringset)


#                           Create FASTA
#____________________________________________________________

write_fasta <-  writeXStringSet(tRNA_dnastringset, 'tRNA_synthetase_genes.fasta', format = 'fasta')
#write_fasta

# prepare a fastq for bowtie2 alignment later
object <-  readFasta('tRNA_synthetase_genes.fasta')
writeFastq(object, 'new', mode="w", full=FALSE, compress=TRUE)
object

tRNA_dnastringset












