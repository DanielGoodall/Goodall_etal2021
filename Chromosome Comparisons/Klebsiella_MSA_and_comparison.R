# Purpose: Klebsiella bowtie ter
# Date created: 11/06/21

##############################################################################################
#                       ALIGNING TEMPLATE TER SEQS AGAINST Salmonella GENOME - BLAST EQUIVALENT
############################################################################################

# Get the seq alignment file of ter first
library(Biostrings)
library(tidyverse)

setwd(dirname(file.choose()))

# create fasta with Biostrings
DX120E <- readDNAStringSet(file.choose())
Kpneu <- readDNAStringSet(file.choose())
Kpneu

genomes_kleb <- c(DX120E, Kpneu)
genomes_kleb

setwd(dirname(file.choose()))

writeXStringSet(genomes, 'Klebsiella_2_genomes.fasta', format = 'fasta')


##########################
#     BOWTIE2 INDEX BUILD
##########################
library(Rbowtie2)

# set wd as the Klebsiella folder
setwd(dirname(file.choose()))
# Build the index (database) of ecoli genomes
bowtie2_build(file.choose(),'Klebsiella_2_genomes')


##########################
#   BOWTIE2 ALIGNMENT
bowtie2_usage()
###########################

#_______________________________________________________________________________________________
#                                             All 10 ter sites

bowtie2(bt2Index = 'Klebsiella_2_genomes', 
        samOutput = 'Klebsiella_ter_11-06.sam', 
        seq1 = 'ter_seq.fastq',
        "-D 15 -R 3 -N 1 -L 15 -i S,1,0.50 -M 2", "--all", "--reorder", "--mp 4,1",
        seq2 = NULL, 
        overwrite = TRUE)





##########################
#     SAMTOOLS VIEW 
##########################
# https://gist.github.com/davetang/6460320
# https://bioconductor.org/packages/release/bioc/vignettes/Rsamtools/inst/doc/Rsamtools-Overview.pdf

library(Rsamtools)

# processing to a dataframe
scan <- scanBam(asBam(file.choose()))
bam <- scanBam(file.choose())
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
write.csv(bam_df, 'Klebsiella_2_genomes.csv')



#                             Clean DF
#__________________________________________________________________________________________

setwd(dirname(file.choose()))

csv <- read.csv('Klebsiella_2_genomes.csv')
csv

# filter each genome
DX120E <- csv %>% 
  filter(rname == 'CP009274.2') %>% 
  mutate(strain = rname, name = qname, start = pos, end = qwidth) %>% 
  select(name, start, end, strand) %>% 
  arrange(name) %>% 
  mutate(end = start+end)
DX120E


Kpneu <- csv %>% 
  filter(rname == 'Klebsiella') %>% 
  mutate(strain = rname, name = qname, start = pos, end = qwidth) %>% 
  select(name, start, end, strand) %>% 
  arrange(name) %>% 
  mutate(end = start+end)
Kpneu

# drop duplicates
DX120E <- DX120E %>% 
  filter(!duplicated(start))
Kpneu <- Kpneu %>% 
  filter(!duplicated(start))

DX120E
Kpneu

# mid point calculation
# (oric-(gen_size-oric))/2
DX120E_mid <- (DX120E %>% filter(name == 'oriC') %>% select(start) - (width(genomes[1]) - DX120E %>% filter(name == 'oriC') %>% select(start)))/2
Kpneu_mid <- (Kpneu %>% filter(name == 'oriC') %>% select(start)-(width(genomes[2]) - Kpneu %>% filter(name == 'oriC') %>% select(start)))/2

# mi point sanity
DX120E_mid
Kpneu_mid


# add mid-point into dfs
DX120E[14,] <- c(name = 'mid', start = DX120E_mid, end = DX120E_mid+23, strand = '+')
DX120E

Kpneu[13,] <- c(name = 'mid', start = Kpneu_mid, end = Kpneu_mid, strand = '+')
Kpneu




#____________________________________________________________________________________________ GenoPlotR
library(genoPlotR)

# sanity
DX120E
Kpneu

# dna_seg, without the outlier 'terI' and OriC
DX120E_dna <- dna_seg(DX120E[-c(2),])
Kpneu_dna <- dna_seg(Kpneu[-c(2),])



# Use plot_gene_map to see which ter sites are true ter sites.
DX120E_dna
Kpneu_dna


#____________________________________________________________________________________ Triangle for ter sites
## Creates a permissible triangle shape grob for ter sites
# plus grobs
plusGrob <- function(gene, ...) {
  plus_x <- c((gene$start+gene$end)/2, (gene$start+gene$end)/2, gene$end+20000)
  plus_y <- c(1.1, -0.1, 0.5)
  polygonGrob(plus_x, plus_y,gp=gpar(fill=gene$fill, col=gene$col, lty=gene$lty,
                                     lwd=gene$lwd), default.units="native")
}

# min grob
minGrob <- function(gene, ...) {
  min_x <- c(gene$start-20000, (gene$start+gene$end)/2, (gene$start+gene$end)/2)
  min_y <- c(0.5, -0.1, 1.1)
  polygonGrob(min_x, min_y,gp=gpar(fill=gene$fill, col=gene$col, lty=gene$lty,
                                   lwd=gene$lwd), default.units="native")
}
#_____________________________________________________________________________________


# specify the grobs
DX120E_dnaR <- reverse(DX120E_dna)
Kpneu_dnaR <- reverse(Kpneu_dna)


# change into trianlges
DX120E_dnaR <- rbind(DX120E_dnaR %>% 
  filter(strand =='1') %>% 
  mutate(gene_type = replace(gene_type, gene_type == 'arrows', 'minGrob')),
  DX120E_dnaR %>% 
    filter(strand =='-1') %>% 
    mutate(gene_type = replace(gene_type, gene_type == 'arrows', 'plusGrob')))
DX120E_dnaR



Kpneu_dnaR <- rbind(Kpneu_dnaR %>% 
                     filter(strand =='1') %>% 
                     mutate(gene_type = replace(gene_type, gene_type == 'arrows', 'minGrob')),
                    Kpneu_dnaR %>% 
                     filter(strand =='-1') %>% 
                     mutate(gene_type = replace(gene_type, gene_type == 'arrows', 'plusGrob')))
Kpneu_dnaR


# colours of the left and right replichores
DX120E_dnaR <- rbind(
  DX120E_dnaR %>% 
    filter(strand =='1') %>%
    mutate(col = replace(col, col == 'blue', 'black'), fill = replace(fill, fill == 'blue', 'royalblue')),
  DX120E_dnaR %>% 
    filter(strand =='-1') %>%
    mutate(col = replace(col, col == 'blue', 'black'), fill = replace(fill, fill == 'blue', 'orangered2'))
) %>% 
  arrange(name)
# sanity
DX120E_dnaR


Kpneu_dnaR <- rbind(
  Kpneu_dnaR %>% 
    filter(strand =='1') %>%
    mutate(col = replace(col, col == 'blue', 'black'), fill = replace(fill, fill == 'blue', 'royalblue')),
  Kpneu_dnaR %>% 
    filter(strand =='-1') %>%
    mutate(col = replace(col, col == 'blue', 'black'), fill = replace(fill, fill == 'blue', 'orangered2'))
) %>% 
  arrange(name)
# sanity
Kpneu_dnaR


# turn dif into its own gene type and colour
DX120E_dnaR[c(1,2),]$gene_type <- 'side_blocks'
Kpneu_dnaR[c(1,2),]$gene_type <- 'side_blocks'


# last step is to invert the mid point for K.pneumoniae as oric is on the other side
Kpneu_dnaR[2,] <- reverse(Kpneu_dnaR[2,])
Kpneu_dnaR

DX120E_dnaR
Kpneu_dnaR

# remove the terI pseudoTerr sites I have identified without GC6
DX120E_dnaR <- DX120E_dnaR[-c(10,11,12),]
Kpneu_dnaR <- Kpneu_dnaR[-c(11,12),]

# change the side block to + for mid in DX120E
DX120E_dnaR[2,]$strand <- 1

# final sanity
DX120E_dnaR
Kpneu_dnaR

# change order of the ter sites by strand
DX120E_dnaR <- DX120E_dnaR %>% arrange(strand)
Kpneu_dnaR <- Kpneu_dnaR %>% arrange(strand)

temp1$name
Kpneu_dnaR$name

# update the name of ter sites to match the MSA names
DX120E_dnaR_names <- c("dif",  "ter1", "ter2", "ter3", "mid",  "ter4", "ter5", "ter6", "ter7", "ter8")
DX120E_dnaR$name <- DX120E_dnaR_names
DX120E_dnaR

Kpneu_dnaR_names <- c("ter1", "ter2", "ter3", "ter4", "dif" , "mid",  "ter5", "ter6", "ter7", "ter8")
Kpneu_dnaR$name <- Kpneu_dnaR_names
Kpneu_dnaR

# middle pos
midpos1 <- middle(DX120E_dnaR)
midpos2 <- middle(Kpneu_dnaR)


# annotations
annot1 <- annotation(midpos1, text = DX120E_dnaR$name, rot = 90, col = 'black')
annot2 <- annotation(midpos2, text = Kpneu_dnaR$name, rot = 90, col = 'black')

# list dna_segs

list_dna <- list(DX120E_dnaR,
                 Kpneu_dnaR)

list_dna

names(list_dna) <- c('K.variicola DX120E','K.pneumoniae')
# list annots
annots <- list(annot1, annot2)

list_dna

plot_gene_map(list_dna, 
              annotations = annots, 
              annotation_height = 5, 
              main = 'Klebsiella', 
              scale = FALSE,
              dna_seg_scale = FALSE)



#________________________________________________________________________________________________ MSA

#_____________________________________________________________________________________________________

#                                   MSA alignment clean
#                   May not need as we can just drop the duplicated start positions (above)
#____________________________________________________________________________________________

library(DECIPHER)
library(msa)

# ______________________________________________ DX120E
DX120E
# drop dif and mid and oric, and the terI pseudoter sites
DX120E_temp <- DX120E[-c(1,2,10,11,12, 14),]

DX120E_plus_temp <- DX120E_temp %>% filter(strand =='+')
DX120E_min_temp <- DX120E_temp %>% filter(strand =='-')

DX120E_plus_temp
DX120E_min_temp

# set empty list
DX120E_seq_plus_list <- list()
DX120E_seq_min_list <- list()


# extract the sequence from the BT2 positions and update the empty list
#plus
for (number in DX120E_plus_temp$start) {
  seq <- list(subseq(genomes_kleb[1], start= number, width = 23))
  DX120E_seq_plus_list <- append(DX120E_seq_plus_list, seq)
}
#minus
for (number in DX120E_min_temp$start) {
  seq <- list(subseq(genomes_kleb[1], start= number, width = 23))
  DX120E_seq_min_list <- append(DX120E_seq_min_list, seq)
}


# collapse into one DNAstringSet
DX120E_seq_plus_list <- do.call(c, DX120E_seq_plus_list)
DX120E_seq_plus_list

DX120E_seq_min_list <- do.call(c, DX120E_seq_min_list)
DX120E_seq_min_list <- reverseComplement(DX120E_seq_min_list)
DX120E_seq_min_list

#__________________________
#       See which ter sites are what
names(DX120E_seq_plus_list) <- DX120E_plus_temp$name
names(DX120E_seq_min_list) <- DX120E_min_temp$name

# combine
combo_DX120E <- c(DX120E_seq_plus_list, DX120E_seq_min_list)
combo_DX120E


ter_numbers <- c("ter1(red)","ter2(red)","ter3(red)","ter4(blue)", "ter5(blue)","ter6(blue)","ter7(blue)", "ter8(blue)")

names(combo_DX120E) <- ter_numbers
combo_DX120E
# MSA on these and exclude based on alignment
# PLUS
# need to not include dif, hence the indexing
DX120E_msa <- msa(combo_DX120E, method = 'Muscle', order = 'input')
DX120E_msa


msaPrettyPrint(x=DX120E_msa, output = 'pdf', file='V2_klebsiella_DX120E_MSA_excludedTer_12-06.pdf',
               shadingMode = 'identical', verbose = TRUE,
               shadingColors = 'grays', shadingModeArg=NA,
               consensusThreshold=80,
               paperWidth=5, paperHeight=4, margins=c(0.3, 0.3))



# ______________________________________________ K.pneumoniae
Kpneu
# drop dif and mid and oric, and the terI pseudoter sites
Kpneu_temp <- Kpneu[-c(1,2,11,12,13),]

Kpneu_plus_temp <- Kpneu_temp %>% filter(strand =='+')
Kpneu_plus_temp
Kpneu_min_temp <- Kpneu_temp %>% filter(strand =='-')
Kpneu_min_temp
# set empty list
Kpneu_seq_plus_list <- list()
Kpneu_seq_min_list <- list()


# extract the sequence from the BT2 positions and update the empty list
#plus
for (number in Kpneu_plus_temp$start) {
  seq <- list(subseq(genomes_kleb[2], start= number, width = 23))
  Kpneu_seq_plus_list <- append(Kpneu_seq_plus_list, seq)
}
#minus
for (number in Kpneu_min_temp$start) {
  seq <- list(subseq(genomes_kleb[2], start= number, width = 23))
  Kpneu_seq_min_list <- append(Kpneu_seq_min_list, seq)
}


# collapse into one DNAstringSet
Kpneu_seq_plus_list <- do.call(c, Kpneu_seq_plus_list)
Kpneu_seq_plus_list

Kpneu_seq_min_list <- do.call(c, Kpneu_seq_min_list)
Kpneu_seq_min_list <- reverseComplement(Kpneu_seq_min_list)
Kpneu_seq_min_list

# combine
combo_Kpneu <- c(Kpneu_seq_plus_list, Kpneu_seq_min_list)
combo_Kpneu


ter_numbers2 <- c("ter1(red)","ter2(red)","ter3(red)","ter4(red)", "ter5(blue)","ter6(blue)","ter7(blue)", "ter8(blue)")

names(combo_Kpneu) <- ter_numbers2
combo_Kpneu
# MSA on these and exclude based on alignment
# PLUS
# need to not include dif, hence the indexing
Kpneu_msa <- msa(combo_Kpneu, method = 'Muscle', order = 'input')
Kpneu_msa


msaPrettyPrint(x=Kpneu_msa, output = 'pdf', file='V2_klebsiella_Pneumoniae_excludedTerMSA_12-06.pdf',
               shadingMode = 'identical', verbose = TRUE,
               shadingColors = 'grays', shadingModeArg=NA,
               consensusThreshold=80,
               paperWidth=5, paperHeight=4, margins=c(0.3, 0.3))
