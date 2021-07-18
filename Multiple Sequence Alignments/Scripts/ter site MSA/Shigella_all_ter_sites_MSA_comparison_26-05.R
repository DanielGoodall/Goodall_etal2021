# Purpose: 4 shigella phylogroups figures with trigrobs - from BT2 large output
# Date updated: 26/05/21

# May have to stop half way through coding as the time pressure is too great

# I need to read in 1000 rows of ter sites, split into terA-J, then add the MG1655 seq to each df,
# then align using msa and select only the site which lines up best to get the true ter site,
# then merge back together to get a complete terA-J df from BT2 origin,
# then do that for the 4 serotypes, when it would be easier to simply read in the BLAST csv files and plot them...


library(tidyverse)
library(Biostrings)
library(genoPlotR)
library(msa)

setwd(dirname(file.choose()))

# read csv of BT2 all ter sites
shig_csv <- read.csv('Shigella_ter_26-05.csv')
shig_csv

# keep the columns we need for genoplotR
shig_df <- shig_csv[,c(1,2,3,4,5,6)]
shig_df

# get rid of the strains we do not need
flex2a <- shig_df %>% 
  filter(rname == 'flex_2a')

sd197 <-  shig_df %>% 
  filter(rname == 'sd197')

ss046<-  shig_df %>% 
  filter(rname == 'ss046')

boy3038 <-  shig_df %>% 
  filter(rname == 'boy_3038')

# sort alphabetically by ter name
flex2a <- arrange(flex2a, qname)
sd197 <- arrange(sd197, qname)
ss046 <- arrange(ss046, qname)
boy3038 <- arrange(boy3038, qname)


# change column names with dplyr rename()
flex2a <- flex2a %>% 
  dplyr::rename(strain = rname, name = qname, start = pos, end = qwidth)

sd197 <- sd197 %>% 
  dplyr::rename(strain = rname, name = qname, start = pos, end = qwidth)

ss046 <- ss046 %>% 
  dplyr::rename(strain = rname, name = qname, start = pos, end = qwidth)

boy3038 <- boy3038 %>% 
  dplyr::rename(strain = rname, name = qname, start = pos, end = qwidth)

# Split into terA-J
# flex2a
terA_flex2a <- 
  flex2a %>% 
  filter(name == 'terA')

terB_flex2a <- 
  flex2a %>% 
  filter(name == 'terB')

terC_flex2a <- 
  flex2a %>% 
  filter(name == 'terC')

terD_flex2a <- 
  flex2a %>% 
  filter(name == 'terD')

terE_flex2a <- 
  flex2a %>% 
  filter(name == 'terE')

terF_flex2a <- 
  flex2a %>% 
  filter(name == 'terF')

terG_flex2a <- 
  flex2a %>% 
  filter(name == 'terG')

terH_flex2a <- 
  flex2a %>% 
  filter(name == 'terH')

terI_flex2a <- 
  flex2a %>% 
  filter(name == 'terI')

terJ_flex2a <- 
  flex2a %>% 
  filter(name == 'terJ')

#sanity
terA_flex2a
terB_flex2a
terC_flex2a
terD_flex2a
terE_flex2a
terF_flex2a
terG_flex2a
terH_flex2a
terI_flex2a
terJ_flex2a

#____________________________________________________________________________________________ Read all ter sequences
ter_set_dna <- readDNAStringSet(file.choose(), format = 'fastq')


#____________________________________________________________________________________________ DNA_SEGS

#_________________________________________________________________________________ flex2a
terA_flex2a_dna <- DNAStringSet(terA_flex2a$seq)
terB_flex2a_dna <- DNAStringSet(terB_flex2a$seq)
terC_flex2a_dna <- DNAStringSet(terC_flex2a$seq)
terD_flex2a_dna <- DNAStringSet(terD_flex2a$seq)
terE_flex2a_dna <- DNAStringSet(terE_flex2a$seq)
terF_flex2a_dna <- DNAStringSet(terF_flex2a$seq)
terG_flex2a_dna <- DNAStringSet(terG_flex2a$seq)
terH_flex2a_dna <- DNAStringSet(terH_flex2a$seq)
terI_flex2a_dna <- DNAStringSet(terI_flex2a$seq)
terJ_flex2a_dna <- DNAStringSet(terJ_flex2a$seq)

# set names 
names(terA_flex2a_dna) <- 1:length(terA_flex2a_dna)
names(terB_flex2a_dna) <- 1:length(terB_flex2a_dna)
names(terC_flex2a_dna) <- 1:length(terC_flex2a_dna)
names(terD_flex2a_dna) <- 1:length(terD_flex2a_dna)
names(terE_flex2a_dna) <- 1:length(terE_flex2a_dna)
names(terF_flex2a_dna) <- 1:length(terF_flex2a_dna)
names(terG_flex2a_dna) <- 1:length(terG_flex2a_dna)
names(terH_flex2a_dna) <- 1:length(terH_flex2a_dna)
names(terI_flex2a_dna) <- 1:length(terI_flex2a_dna)
names(terJ_flex2a_dna) <- 1:length(terJ_flex2a_dna)

# bind MG1655 ter site to the DNA_seg in first position with their name for reference
terA_flex2a_dna <- c(ter_set_dna[1], terA_flex2a_dna)
terB_flex2a_dna <- c(ter_set_dna[2], terB_flex2a_dna)
terC_flex2a_dna <- c(ter_set_dna[3], terC_flex2a_dna)
terD_flex2a_dna <- c(ter_set_dna[4], terD_flex2a_dna)
terE_flex2a_dna <- c(ter_set_dna[5], terE_flex2a_dna)
terF_flex2a_dna <- c(ter_set_dna[6], terF_flex2a_dna)
terG_flex2a_dna <- c(ter_set_dna[7], terG_flex2a_dna)
terH_flex2a_dna <- c(ter_set_dna[8], terH_flex2a_dna)
terI_flex2a_dna <- c(ter_set_dna[9], terI_flex2a_dna)
terJ_flex2a_dna <- c(ter_set_dna[10], terJ_flex2a_dna)


# sanity
terA_flex2a_dna
terB_flex2a_dna
terC_flex2a_dna
terD_flex2a_dna 
terE_flex2a_dna
terF_flex2a_dna 
terG_flex2a_dna
terH_flex2a_dna
terI_flex2a_dna
terJ_flex2a_dna


# MSA and update the real ter site to a df
msa(terA_flex2a_dna, method = 'ClustalW', order='aligned')


















#_________________________________________________________________________________ sd197
terA_sd197
terB_sd197
terC_sd197
terD_sd197
terE_sd197
terF_sd197
terG_sd197
terH_sd197
terI_sd197
terJ_sd197

#_________________________________________________________________________________ ss046
terA_ss046
terB_ss046
terC_ss046
terD_ss046
terE_ss046
terF_ss046
terG_ss046
terH_ss046
terI_ss046
terJ_ss046

#_________________________________________________________________________________ boy3038
terA_boy3038
terB_boy3038
terC_boy3038
terD_boy3038
terE_boy3038
terF_boy3038
terG_boy3038
terH_boy3038
terI_boy3038
terJ_boy3038






# delete rows with duplicates manually using the BLAST csv spreadsheetn for indexes
# Become laborious when we have 100 hits for all ter sites
flex2a %>% 
  dplyr::slice(-c(2,3,7,8,10,12,13))


sd197 %>% 
  dplyr::slice(-c(1,3,4,5,8,10,))

# sanity
flex2a
sd197
ss046
boy3038










# ________________________________test through genoplotR to see if positions are accurate 

# clean
flex2a_clean1 <- flex2a[!duplicated(flex2a$pos),][2:5]

colnames(flex2a_clean1) <- c('name', 'start', 'end', 'strand')

flex2a_clean1$end <- flex2a_clean1$start + 5000
flex2a_clean1
#___________________________________________ dna_seg

flex2a_dna <- dna_seg(flex2a_clean1)
flex2a_dna

##################################################################################
## Creates a permissible triangle shape grob for ter sites
# plus grobs
plusGrob <- function(gene, ...) {
  plus_x <- c((gene$start+gene$end)/2, (gene$start+gene$end)/2, gene$end)
  plus_y <- c(1.1, -0.1, 0.5)
  polygonGrob(plus_x, plus_y,gp=gpar(fill=gene$fill, col=gene$col, lty=gene$lty,
                                     lwd=gene$lwd), default.units="native")
}

# min grob
minGrob <- function(gene, ...) {
  min_x <- c(gene$start, (gene$start+gene$end)/2, (gene$start+gene$end)/2)
  min_y <- c(0.5, -0.1, 1.1)
  polygonGrob(min_x, min_y,gp=gpar(fill=gene$fill, col=gene$col, lty=gene$lty,
                                   lwd=gene$lwd), default.units="native")
}
###################################################################################



# Turn gene_type into plusGrob if + strand
flex2a_dna$gene_type[flex2a_dna$strand==1] <- 'plusGrob'
flex2a_dna$gene_type[flex2a_dna$strand==-1] <- 'minGrob'

dna_segs <- list(reverse(flex2a_dna))
names(dna_segs) <- 'flex'

dna_segs
mid_pos1 <- middle(dna_segs[[1]])
mid_pos1

annot1 <- annotation(x1 = mid_pos1, text = dna_segs[[1]]$name, rot = 90)
annot_list <- list(annot1)
annot1


#plot
plot_gene_map(dna_segs = list(reverse(flex2a_dna)),
              annotations = annot_list)

