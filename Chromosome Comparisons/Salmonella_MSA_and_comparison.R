# Purpose: Salmonella bowtie ter
# Date created: 10/06/21

##############################################################################################
#                       ALIGNING TEMPLATE TER SEQS AGAINST Salmonella GENOME
############################################################################################

# Get the seq alignment file of ter first
library(Biostrings)
library(tidyverse)


setwd(dirname(file.choose()))

# create fasta with Biostrings
DT104_fas <- readDNAStringSet(file.choose())
FDAARGOS_878_fas <- readDNAStringSet(file.choose())

genomes <- c(DT104_fas,FDAARGOS_878_fas)
genomes

writeXStringSet(genomes, 'Salmonella_3_genomes.fasta', format = 'fasta')


##########################
#     BOWTIE2 INDEX BUILD
##########################
library(Rbowtie2)

# set wd as the salmonella folder
setwd(dirname(file.choose()))
# Build the index (database) of ecoli genomes
bowtie2_build(file.choose(),'Salmonella_3_genomes')

# Build the index (database) of ter sequences
bowtie2_build(file.choose(),'ter_seq')

##########################
#   BOWTIE2 ALIGNMENT
bowtie2_usage()
###########################

#_______________________________________________________________________________________________
#                                             All 10 ter sites

bowtie2(bt2Index = 'Salmonella_3_genomes', 
        samOutput = 'Salmonella_ter_11-06.sam', 
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
write.csv(bam_df, 'Salmonella_3_genomes_ter_dif.csv')



#                             Clean DF
#__________________________________________________________________________________________

setwd(dirname(file.choose()))

csv <- read.csv('Salmonella_3_genomes_ter_dif.csv')
csv

# filter each genome
DT104 <- csv %>% 
  filter(rname == 'HF937208.1') %>% 
  mutate(strain = rname, name = qname, start = pos, end = qwidth) %>% 
  select(name, start, end, strand) %>% 
  arrange(name) %>% 
  mutate(end = start+end)
DT104



FDAARGOS_878 <- csv %>% 
  filter(rname == 'NZ_CP065718.1') %>% 
  mutate(strain = rname, name = qname, start = pos, end = qwidth) %>% 
  select(name, start, end, strand) %>% 
  arrange(name) %>% 
  mutate(end = start+end)
FDAARGOS_878

# drop duplicates
DT104 <- DT104 %>% 
  filter(!duplicated(start))
FDAARGOS_878 <- FDAARGOS_878 %>% 
  filter(!duplicated(start))
# sanity
DT104
FDAARGOS_878

# issue with FDAARGOS_878 mid point being too far out, oric is on other strand so could explain it.
width(genomes[2]) - FDAARGOS_878 %>% filter(name == 'oriC') %>% select(start)

# mid point calculation
# (oric-(gen_size-oric))/2
DT104_mid <- (DT104 %>% filter(name == 'oriC') %>% select(start) - (width(genomes[1]) - DT104 %>% filter(name == 'oriC') %>% select(start)))/2
FDAARGOS878_mid <- (FDAARGOS_878 %>% filter(name == 'oriC') %>% select(start) - (width(genomes[2]) - FDAARGOS_878 %>% filter(name == 'oriC') %>% select(start)))/2 + width(genomes[2])
width(genomes[2])

# mi point sanity
DT104_mid
FDAARGOS878_mid


# add mid-point into dfs
DT104[13,] <- c(name = 'mid', start = DT104_mid, end = DT104_mid+23, strand = '+')
FDAARGOS_878[13,] <- c(name = 'mid', start = FDAARGOS878_mid, end = FDAARGOS878_mid+23, strand = '+')

DT104
FDAARGOS_878

# drop ter sites not in place (terC for FDAARGOS_878)
# also drop OriC for DT104, keep for FDAARGOS_878
DT104 <- DT104[-c(2),]
FDAARGOS_878 <- FDAARGOS_878[-c(10),]


# sanity
rownames(DT104) <- NULL
rownames(FDAARGOS_878) <- NULL



DT104$start

# change to numeric as the position were characters
FDAARGOS_878$start <- as.numeric(FDAARGOS_878$start)
FDAARGOS_878$end <- as.numeric(FDAARGOS_878$end)

#____________________________________________________________________________________________
library(genoPlotR)

# dna_seg
DT104_dna <- dna_seg(DT104)
FDAARGOS878_dna <- dna_seg(FDAARGOS_878)


# Use plot_gene_map to see which ter sites are true ter sites.
DT104_dna
FDAARGOS878_dna





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
DT104_dnaR <- reverse(DT104_dna)
FDAARGOS878_dnaR <- reverse(FDAARGOS878_dna)


# change into trianlges
DT104_dnaR <- rbind(DT104_dnaR %>% 
                       filter(strand =='1') %>% 
                       mutate(gene_type = replace(gene_type, gene_type == 'arrows', 'minGrob')),
                    DT104_dnaR %>% 
                       filter(strand =='-1') %>% 
                       mutate(gene_type = replace(gene_type, gene_type == 'arrows', 'plusGrob')))
DT104_dnaR



FDAARGOS878_dnaR <- rbind(FDAARGOS878_dnaR %>% 
                      filter(strand =='1') %>% 
                      mutate(gene_type = replace(gene_type, gene_type == 'arrows', 'minGrob')),
                      FDAARGOS878_dnaR %>% 
                      filter(strand =='-1') %>% 
                      mutate(gene_type = replace(gene_type, gene_type == 'arrows', 'plusGrob')))
FDAARGOS878_dnaR


# colours of the left and right replichores
DT104_dnaR <- rbind(
  DT104_dnaR %>% 
    filter(strand =='1') %>%
    mutate(col = replace(col, col == 'blue', 'black'), fill = replace(fill, fill == 'blue', 'royalblue')),
  DT104_dnaR %>% 
    filter(strand =='-1') %>%
    mutate(col = replace(col, col == 'blue', 'black'), fill = replace(fill, fill == 'blue', 'orangered2'))
)
# sanity
DT104_dnaR


FDAARGOS878_dnaR <- rbind(
  FDAARGOS878_dnaR %>% 
    filter(strand =='1') %>%
    mutate(col = replace(col, col == 'blue', 'black'), fill = replace(fill, fill == 'blue', 'royalblue')),
  FDAARGOS878_dnaR %>% 
    filter(strand =='-1') %>%
    mutate(col = replace(col, col == 'blue', 'black'), fill = replace(fill, fill == 'blue', 'orangered2'))
)
# sanity
FDAARGOS878_dnaR <- FDAARGOS878_dnaR[-c(6),]
FDAARGOS878_dnaR

# turn dif into its own gene type and colour
DT104_dnaR[c(5,12),]$gene_type <- 'side_blocks'
FDAARGOS878_dnaR[c(5,11),]$gene_type <- 'side_blocks'

# set the names of each ter site as the same as is shown in the MSA figure.
DT104_ter_num <- c("ter1", "ter2", "ter3", "ter4", "dif", "ter5", "ter6", "ter7", "ter8", "ter9", "ter10", "mid")
DT104_dnaR$name <- DT104_ter_num
DT104_dnaR

FDAARGOS878_ter_num <- c("ter1", "ter2", "ter3", "ter4", "dif", "ter5", "ter6", "ter7", "ter8", "ter9", "mid")
FDAARGOS878_dnaR$name <- FDAARGOS878_ter_num
FDAARGOS878_dnaR

# last step is to invert the mid point for FDAARGOS as oric is on the other side (if its a negative value)
#FDAARGOS878_dnaR[c(12),] <- reverse(FDAARGOS878_dnaR[c(12),])
#FDAARGOS878_dnaR

# middle pos
midpos1 <- middle(DT104_dnaR)
midpos2 <- middle(FDAARGOS878_dnaR)

# annotations
annot1 <- annotation(midpos1, text = DT104_dnaR$name, rot = 90, col = 'black')
annot2 <- annotation(midpos2, text = FDAARGOS878_dnaR$name, rot = 90, col = 'black')

# list dna_segs
list_dna <- list(DT104_dnaR,
                 FDAARGOS878_dnaR)

list_dna

names(list_dna) <- c('DT104','FDAARGOS878')
# list annots
annots <- list(annot1, annot2)


plot_gene_map(list_dna, 
              annotations = annots, 
              annotation_height = 6, 
              main = 'Salmonella', 
              dna_seg_scale = FALSE, 
              scale = FALSE)



#________________________________________________________________________________________________ MSA

#_____________________________________________________________________________________________________

#                                   MSA alignment clean
#                   May not need as we can just drop the duplicated start positions (above)
#____________________________________________________________________________________________

library(DECIPHER)
library(msa)

# ______________________________________________ DT104
DT104
DT104[-c(1,12),]
DT104_temp <- DT104[-c(1,12),]

DT104_plus_temp <- DT104_temp %>% filter(strand =='+')
DT104_min_temp <- DT104_temp %>% filter(strand =='-')

DT104_plus_temp
DT104_min_temp

# set empty list
DT104_seq_plus_list <- list()
DT104_seq_min_list <- list()


# extract the sequence from the BT2 positions and update the empty list
#plus
for (number in DT104_plus_temp$start) {
  seq <- list(subseq(genomes[1], start= number, width = 23))
  DT104_seq_plus_list <- append(DT104_seq_plus_list, seq)
}
#minus
for (number in DT104_min_temp$start) {
  seq <- list(subseq(genomes[1], start= number, width = 23))
  DT104_seq_min_list <- append(DT104_seq_min_list, seq)
}


# collapse into one DNAstringSet
DT104_seq_plus_list <- do.call(c, DT104_seq_plus_list)
DT104_seq_plus_list

DT104_seq_min_list <- do.call(c, DT104_seq_min_list)
DT104_seq_min_list <- reverseComplement(DT104_seq_min_list)
DT104_seq_min_list

# combine
combo_DT104 <- c(DT104_seq_plus_list, DT104_seq_min_list)
combo_DT104


ter_numbers <- c("ter1","ter2","ter3","ter4", "ter5","ter6","ter7", "ter8", "ter9","ter10")

names(combo_DT104) <- ter_numbers
combo_DT104
# MSA on these and exclude based on alignment
# PLUS
# need to not include dif, hence the indexing
DT104_msa <- msa(combo_DT104, method = 'Muscle', order = 'input')
DT104_msa

setwd(dirname(file.choose()))

msaPrettyPrint(x=DT104_msa, output = 'pdf', file='Salmonella_DT104_MSA_12-06.pdf',
               shadingMode = 'identical', verbose = TRUE,
               shadingColors = 'grays', shadingModeArg=NA,
               consensusThreshold=90,
               paperWidth=5, paperHeight=4, margins=c(0.3, 0.3))



# ______________________________________________
FDAARGOS_878

FDAARGOS878_temp <- FDAARGOS_878[-c(1,11),]

FDAARGOS878_plus_temp <- FDAARGOS878_temp %>% filter(strand =='+')
FDAARGOS878_min_temp <- FDAARGOS878_temp %>% filter(strand =='-')

# set empty list
FDAARGOS878_seq_plus_list <- list()
FDAARGOS878_seq_min_list <- list()


# extract the sequence from the BT2 positions and update the empty list
#plus
for (number in FDAARGOS878_plus_temp$start) {
  seq <- list(subseq(genomes[2], start= number, width = 23))
  FDAARGOS878_seq_plus_list <- append(FDAARGOS878_seq_plus_list, seq)
}
#minus
for (number in FDAARGOS878_min_temp$start) {
  seq <- list(subseq(genomes[2], start= number, width = 23))
  FDAARGOS878_seq_min_list <- append(FDAARGOS878_seq_min_list, seq)
}


# collapse into one DNAstringSet
FDAARGOS878_seq_plus_list <- do.call(c, FDAARGOS878_seq_plus_list)
FDAARGOS878_seq_plus_list

FDAARGOS878_seq_min_list <- do.call(c, FDAARGOS878_seq_min_list)
FDAARGOS878_seq_min_list <- reverseComplement(FDAARGOS878_seq_min_list)
FDAARGOS878_seq_min_list

# combine
combo_FDAARGOS878 <- c(FDAARGOS878_seq_plus_list, FDAARGOS878_seq_min_list)
combo_FDAARGOS878


ter_numbers2 <- c("ter1","ter2","ter3","ter4", "ter5","ter6","ter7", "ter8", "ter9")

names(combo_FDAARGOS878) <- ter_numbers2
combo_FDAARGOS878
# MSA on these and exclude based on alignment
# PLUS
# need to not include dif, hence the indexing
FDAARGOS878_msa <- msa(combo_FDAARGOS878, method = 'Muscle', order = 'input')
FDAARGOS878_msa


msaPrettyPrint(x=FDAARGOS878_msa, output = 'pdf', file='Salmonella_FDAARGOS878_MSA_12-06.pdf',
               shadingMode = 'identical', verbose = TRUE,
               shadingColors = 'grays', shadingModeArg=NA,
               consensusThreshold=85,
               paperWidth=5, paperHeight=4, margins=c(0.3, 0.3))
