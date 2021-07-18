# Purpose: pTer with terFGJ all Ecoli genomes
# Date updated: 25/05/21


library(Biostrings)
library(seqinr)
library(genoPlotR)
library(tibble)
library(tidyverse) # to clean duplicates
library(msa)

# set working directory before proceeding
setwd(dirname(file.choose()))


#______________________________________________________________________________________
#
#   Subseq the ter sites to see if they are correct before moving on
#                                 MSA matching
#_____________________________________________________________________________________

terK <- DNAStringSet('CGATTGAGAGTTGTAATGAAGTC')
terL <- DNAStringSet('GCACTGGGTGTTGTAATGACGCA')
terY <- DNAStringSet('TATGGGTACGTTGTAATTAGGGA')
terZ <- DNAStringSet('GCTCTCGTTACAACCTGCGGGTA')

# Would be easier to join the two csv from BT2 and BioC



# read in the csv from BT2 and samtools
df <- read.csv('pseudo_ter_Ecoli.csv')
df #sanity

#             Clean the tRNA synthetase csv 
#           Cols = rname[2], qname[3], pos[4], strand[6]
#                 (Strain)    (Gene)      
#___________________________________________________________________________

seqs <- df[c(2,3,4,5,6)] # subset the columns we want to work with genoPlotR
colnames(seqs) <- c('strain','name', 'start', 'end', 'strand') # reset column names
seqs$end <- seqs$start + 1000
seqs #sanity

toString(seqs$strain)

# subset the large df into each strain
# then later test with just the 5 phylogroups
#A
MG1655 <- seqs %>% filter(strain=='MG1655')
BW2952 <- seqs %>% filter(strain=='BW2952')
REL606 <- seqs %>% filter(strain=='REL606')
#B1
APEC078 <- seqs %>% filter(strain=='APECO78')
IAI1 <- seqs %>% filter(strain=='IAI1')
E11368 <- seqs %>% filter(strain=='11368')
#B2
S88 <- seqs %>% filter(strain=='S88')
UTI89 <- seqs %>% filter(strain=='UTI89')
E2348 <- seqs %>% filter(strain=='E2348/69')
#D
IAI39 <- seqs %>% filter(strain=='IAI39')
SMS35 <- seqs %>% filter(strain=='SMS-3-5')
UMN026 <- seqs %>% filter(strain=='UMN026')
CE10 <- seqs %>% filter(strain=='CE10')
D042 <- seqs %>% filter(strain=='042')
#E
TW14359 <- seqs %>% filter(strain=='TW14359')
Sakai <- seqs %>% filter(strain=='Sakai')
EDL933 <- seqs %>% filter(strain=='EDL933')

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



#__________________________________________________________________________ Read terY from BioC matching

terY_BioC <- read.csv('terY_Bioc_Ecoli.csv')
terY_BioC

# Clean using dplyr inot dna_seg format for merging two csvs
#A
MG1655_terY <- 
  terY_BioC %>% 
  filter(names=='MG1655') %>% 
  mutate(name ='terY') %>% 
  select(names, name, start, end, strand) %>% 
  dplyr::rename(strain = names)
  

BW2952_terY <- 
  terY_BioC %>% 
  filter(names=='BW2952') %>% 
  mutate(name ='terY') %>% 
  select(names, name, start, end, strand) %>% 
  dplyr::rename(strain = names)


REL606_terY <- 
  terY_BioC %>% 
  filter(names=='REL606') %>% 
  mutate(name ='terY') %>% 
  select(names, name, start, end, strand) %>% 
  dplyr::rename(strain = names)
#B1
APEC078_terY <- 
  terY_BioC %>% 
  filter(names=='APECO78') %>% 
  mutate(name ='terY') %>% 
  select(names, name, start, end, strand) %>% 
  dplyr::rename(strain = names)

IAI1_terY <- 
  terY_BioC %>% 
  filter(names=='IAI1') %>% 
  mutate(name ='terY') %>% 
  select(names, name, start, end, strand) %>% 
  dplyr::rename(strain = names)

E11368_terY <- 
  terY_BioC %>% 
  filter(names=='11368') %>% 
  mutate(name ='terY') %>% 
  select(names, name, start, end, strand) %>% 
  dplyr::rename(strain = names)

#B2
S88_terY <- 
  terY_BioC %>% 
  filter(names=='S88') %>% 
  mutate(name ='terY') %>% 
  select(names, name, start, end, strand) %>% 
  dplyr::rename(strain = names)

UTI89_terY <- 
  terY_BioC %>% 
  filter(names=='UTI89') %>% 
  mutate(name ='terY') %>% 
  select(names, name, start, end, strand) %>% 
  dplyr::rename(strain = names)

E2348_terY <- 
  terY_BioC %>% 
  filter(names=='E2348/69 ') %>% 
  mutate(name ='terY') %>% 
  select(names, name, start, end, strand) %>% 
  dplyr::rename(strain = names)

#D
IAI39_terY <- 
  terY_BioC %>% 
  filter(names=='IAI39') %>% 
  mutate(name ='terY') %>% 
  select(names, name, start, end, strand) %>% 
  dplyr::rename(strain = names)

SMS35_terY <- 
  terY_BioC %>% 
  filter(names=='SMS-3-5') %>% 
  mutate(name ='terY') %>% 
  select(names, name, start, end, strand) %>% 
  dplyr::rename(strain = names)

# UMN026 not present

CE10_terY <- 
  terY_BioC %>% 
  filter(names=='CE10') %>% 
  mutate(name ='terY') %>% 
  select(names, name, start, end, strand) %>% 
  dplyr::rename(strain = names)

D042_terY <- 
  terY_BioC %>% 
  filter(names=='D042') %>% 
  mutate(name ='terY') %>% 
  select(names, name, start, end, strand) %>% 
  dplyr::rename(strain = names)

#E
TW14359_terY <- 
  terY_BioC %>% 
  filter(names=='TW14359') %>% 
  mutate(name ='terY') %>% 
  select(names, name, start, end, strand) %>% 
  dplyr::rename(strain = names)

Sakai_terY <- 
  terY_BioC %>% 
  filter(names=='Sakai') %>% 
  mutate(name ='terY') %>% 
  select(names, name, start, end, strand) %>% 
  dplyr::rename(strain = names)

EDL933_terY <- 
  terY_BioC %>% 
  filter(names=='EDL933') %>% 
  mutate(name ='terY') %>% 
  select(names, name, start, end, strand) %>% 
  dplyr::rename(strain = names)



# change names to strain then drop strain later
MG1655_terY
BW2952_terY
REL606_terY
APEC078_terY
IAI1_terY
E11368_terY
S88_terY
UTI89_terY
E2348_terY
IAI39_terY
SMS35_terY
CE10_terY
D042_terY
TW14359_terY
Sakai_terY
EDL933_terY

# Merge two dfs showing all pter sites and drop duplicates, keeping only dna_seg cols
#A
MG1655_all <- rbind(MG1655, MG1655_terY) %>% 
  distinct(start, .keep_all = T) %>% 
  select(name, start, end, strand)

BW2952_all <- rbind(BW2952, BW2952_terY) %>% 
  distinct(start, .keep_all = T) %>% 
  select(name, start, end, strand)

REL606_all <- rbind(REL606, REL606_terY) %>% 
  distinct(start, .keep_all = T) %>% 
  select(name, start, end, strand)

#B1
APEC078_all <- rbind(APEC078, APEC078_terY) %>%
  distinct(start, .keep_all = T) %>% 
  select(name, start, end, strand)

IAI1_all <- rbind(IAI1, IAI1_terY) %>% 
  distinct(start, .keep_all = T) %>% 
  select(name, start, end, strand)

E11368_all <- rbind(E11368, E11368_terY) %>% 
  distinct(start, .keep_all = T) %>% 
  select(name, start, end, strand)

#B2
S88_all <- rbind(S88, S88_terY) %>% 
  distinct(start, .keep_all = T) %>% 
  select(name, start, end, strand)

UTI89_all <- rbind(UTI89, UTI89_terY) %>% 
  distinct(start, .keep_all = T) %>% 
  select(name, start, end, strand)

E2348_all <- rbind(E2348, E2348_terY) %>% 
  distinct(start, .keep_all = T) %>% 
  select(name, start, end, strand)

#D
IAI39_all <- rbind(IAI39, IAI39_terY) %>% 
  distinct(start, .keep_all = T) %>% 
  select(name, start, end, strand)

SMS35_all <- rbind(SMS35, SMS35_terY) %>% 
  distinct(start, .keep_all = T) %>% 
  select(name, start, end, strand)

UMN026_all <- UMN026 %>% 
  distinct(start, .keep_all = T) %>% 
  select(name, start, end, strand)

CE10_all <- rbind(CE10, CE10_terY) %>% 
  distinct(start, .keep_all = T) %>% 
  select(name, start, end, strand)


D042_all <- rbind(D042, D042_terY) %>% 
  distinct(start, .keep_all = T) %>% 
  select(name, start, end, strand)

#E
TW14359_all <- rbind(TW14359, TW14359_terY) %>% 
  distinct(start, .keep_all = T) %>% 
  select(name, start, end, strand)

Sakai_all <- rbind(Sakai, Sakai_terY) %>% 
  distinct(start, .keep_all = T) %>% 
  select(name, start, end, strand)

EDL933_all <- rbind(EDL933, EDL933_terY) %>% 
  distinct(start, .keep_all = T) %>% 
  select(name, start, end, strand)



# SANIITY
MG1655_all
BW2952_all
REL606_all
APEC078_all
IAI1_all
E11368_all
S88_all
UTI89_all
E2348_all
IAI39_all
SMS35_all
UMN026_all
CE10_all
D042_all
TW14359_all
Sakai_all
EDL933_all


##########################################
#
#   NEW WAY TO CLEAN DATAFRAMES WITH LAPPLY
#       
##########################################

# Drop the Strain column and add the end column
# without specifying [-1] for every df
list_pter <- list(
  dna_seg(MG1655_all),
  dna_seg(BW2952_all),
  dna_seg(REL606_all),
  dna_seg(APEC078_all),
  dna_seg(IAI1_all),
  dna_seg(E11368_all),
  dna_seg(S88_all),
  dna_seg(UTI89_all),
  dna_seg(E2348_all),
  dna_seg(IAI39_all),
  dna_seg(SMS35_all),
  dna_seg(UMN026_all),
  dna_seg(CE10_all),
  dna_seg(D042_all),
  dna_seg(TW14359_all),
  dna_seg(Sakai_all),
  dna_seg(EDL933_all))
# set names of lists
names(list_pter) <- c(
  'MG1655',
  'BW2952',
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

list_pter




# plot pTer tosee what we are working with, no labels etc just yet

plot_gene_map(dna_segs = list_pter,
              gene_type= 'side_blocks')






#______________________________________________________________________________ Add terFGJ onto pTer

setwd(dirname(file.choose()))

# read in all csvs from BT2 into on tibble!
# https://stackoverflow.com/questions/11433432/how-to-import-multiple-csv-files-at-once

tbl <-
  list.files(pattern = "*.csv") %>% 
  map_df(~read_csv(., col_types = cols(.default = "c")))
tbl
# clean tbl in steps to get only terFGJ in dna_seg format

# MG1655
MG1655_terFGJ <- data.frame(tbl) %>% 
  filter(rname == 'MG1655') %>%
  dplyr::rename(strain = rname, name = qname, start = pos) %>% 
  mutate(end = as.numeric(start) + 23) %>% 
  select(name, start, end, strand) %>% 
  filter(name %in% c('terF', 'terJ' , 'terG'))

MG1655_terFGJ # sanity

# BW2952
BW2952_terFGJ <- data.frame(tbl) %>% 
  filter(rname == 'BW2952') %>%
  dplyr::rename(strain = rname, name = qname, start = pos) %>% 
  mutate(end = as.numeric(start) + 23) %>% 
  select(name, start, end, strand) %>% 
  filter(name %in% c('terF', 'terJ' , 'terG'))

BW2952_terFGJ # sanity

# REL606
REL606_terFGJ <- data.frame(tbl) %>% 
  filter(rname == 'REL606') %>%
  dplyr::rename(strain = rname, name = qname, start = pos) %>% 
  mutate(end = as.numeric(start) + 23) %>% 
  select(name, start, end, strand) %>% 
  filter(name %in% c('terF', 'terJ' , 'terG'))

REL606_terFGJ # sanity

# APEC078
APEC078_terFGJ <- data.frame(tbl) %>% 
  filter(rname == 'APECO78') %>%
  dplyr::rename(strain = rname, name = qname, start = pos) %>% 
  mutate(end = as.numeric(start) + 23) %>% 
  select(name, start, end, strand) %>% 
  filter(name %in% c('terF', 'terJ' , 'terG'))

APEC078_terFGJ # sanity

# IAI1
IAI1_terFGJ <- data.frame(tbl) %>% 
  filter(rname == 'IAI1') %>%
  dplyr::rename(strain = rname, name = qname, start = pos) %>% 
  mutate(end = as.numeric(start) + 23) %>% 
  select(name, start, end, strand) %>% 
  filter(name %in% c('terF', 'terJ' , 'terG'))

IAI1_terFGJ # sanity

# E11368
E11368_terFGJ <- data.frame(tbl) %>% 
  filter(rname == '11368') %>%
  dplyr::rename(strain = rname, name = qname, start = pos) %>% 
  mutate(end = as.numeric(start) + 23) %>% 
  select(name, start, end, strand) %>% 
  filter(name %in% c('terF', 'terJ' , 'terG'))

E11368_terFGJ # sanity

# S88
S88_terFGJ <- data.frame(tbl) %>% 
  filter(rname == 'S88') %>%
  dplyr::rename(strain = rname, name = qname, start = pos) %>% 
  mutate(end = as.numeric(start) + 23) %>% 
  select(name, start, end, strand) %>% 
  filter(name %in% c('terF', 'terJ' , 'terG'))

S88_terFGJ # sanity

# UTI89
UTI89_terFGJ <- data.frame(tbl) %>% 
  filter(rname == 'UTI89') %>%
  dplyr::rename(strain = rname, name = qname, start = pos) %>% 
  mutate(end = as.numeric(start) + 23) %>% 
  select(name, start, end, strand) %>% 
  filter(name %in% c('terF', 'terJ' , 'terG'))

UTI89_terFGJ # sanity

# E2348
E2348_terFGJ <- data.frame(tbl) %>% 
  filter(rname == 'E2348/69') %>%
  dplyr::rename(strain = rname, name = qname, start = pos) %>% 
  mutate(end = as.numeric(start) + 23) %>% 
  select(name, start, end, strand) %>% 
  filter(name %in% c('terF', 'terJ' , 'terG'))

E2348_terFGJ # sanity

# IAI39
IAI39_terFGJ <- data.frame(tbl) %>% 
  filter(rname == 'IAI39') %>%
  dplyr::rename(strain = rname, name = qname, start = pos) %>% 
  mutate(end = as.numeric(start) + 23) %>% 
  select(name, start, end, strand) %>% 
  filter(name %in% c('terF', 'terJ' , 'terG'))

IAI39_terFGJ # sanity

# SMS35
SMS35_terFGJ <- data.frame(tbl) %>% 
  filter(rname == 'SMS-3-5') %>%
  dplyr::rename(strain = rname, name = qname, start = pos) %>% 
  mutate(end = as.numeric(start) + 23) %>% 
  select(name, start, end, strand) %>% 
  filter(name %in% c('terF', 'terJ' , 'terG'))

SMS35_terFGJ # sanity

# UMN026
UMN026_terFGJ <- data.frame(tbl) %>% 
  filter(rname == 'UMN026') %>%
  dplyr::rename(strain = rname, name = qname, start = pos) %>% 
  mutate(end = as.numeric(start) + 23) %>% 
  select(name, start, end, strand) %>% 
  filter(name %in% c('terF', 'terJ' , 'terG'))

UMN026_terFGJ # sanity

# CE10
CE10_terFGJ <- data.frame(tbl) %>% 
  filter(rname == 'CE10') %>%
  dplyr::rename(strain = rname, name = qname, start = pos) %>% 
  mutate(end = as.numeric(start) + 23) %>% 
  select(name, start, end, strand) %>% 
  filter(name %in% c('terF', 'terJ' , 'terG'))

CE10_terFGJ # sanity

# D042
D042_terFGJ <- data.frame(tbl) %>% 
  filter(rname == '42') %>%
  dplyr::rename(strain = rname, name = qname, start = pos) %>% 
  mutate(end = as.numeric(start) + 23) %>% 
  select(name, start, end, strand) %>% 
  filter(name %in% c('terF', 'terJ' , 'terG'))

D042_terFGJ # sanity

# TW14359
TW14359_terFGJ <- data.frame(tbl) %>% 
  filter(rname == 'TW14359') %>%
  dplyr::rename(strain = rname, name = qname, start = pos) %>% 
  mutate(end = as.numeric(start) + 23) %>% 
  select(name, start, end, strand) %>% 
  filter(name %in% c('terF', 'terJ' , 'terG'))

TW14359_terFGJ # sanity

# Sakai
Sakai_terFGJ <- data.frame(tbl) %>% 
  filter(rname == 'Sakai') %>%
  dplyr::rename(strain = rname, name = qname, start = pos) %>% 
  mutate(end = as.numeric(start) + 23) %>% 
  select(name, start, end, strand) %>% 
  filter(name %in% c('terF', 'terJ' , 'terG'))

Sakai_terFGJ # sanity

# EDL933
EDL933_terFGJ <- data.frame(tbl) %>% 
  filter(rname == 'EDL933') %>%
  dplyr::rename(strain = rname, name = qname, start = pos) %>% 
  mutate(end = as.numeric(start) + 23) %>% 
  select(name, start, end, strand) %>% 
  filter(name %in% c('terF', 'terJ' , 'terG'))

EDL933_terFGJ # sanity



#SANIITY
MG1655_terFGJ
BW2952_terFGJ
REL606_terFGJ
APEC078_terFGJ
IAI1_terFGJ
E11368_terFGJ
S88_terFGJ
UTI89_terFGJ
E2348_terFGJ
IAI39_terFGJ
SMS35_terFGJ
UMN026_terFGJ
CE10_terFGJ
D042_terFGJ
TW14359_terFGJ
Sakai_terFGJ
EDL933_terFGJ


#_____________________________________________________________________________ merge pTer and terFGJ


# merge
MG1655_ter <- rbind(MG1655_terFGJ, MG1655_all)
BW2952_ter <- rbind(BW2952_terFGJ, BW2952_all)
REL606_ter <- rbind(REL606_terFGJ, REL606_all)
APEC078_ter <- rbind(APEC078_terFGJ, APEC078_all)
IAI1_ter <- rbind(IAI1_terFGJ, IAI1_all)
E11368_ter <- rbind(E11368_terFGJ, E11368_all)
S88_ter <- rbind(S88_terFGJ, S88_all)
UTI89_ter <- rbind(UTI89_terFGJ, UTI89_all)
E2348_ter <- rbind(E2348_terFGJ, E2348_all)
IAI39_ter <- rbind(IAI39_terFGJ, IAI39_all)
SMS35_ter <- rbind(SMS35_terFGJ, SMS35_all)
UMN026_ter <- rbind(UMN026_terFGJ, UMN026_all)
CE10_ter <- rbind(CE10_terFGJ, CE10_all)
D042_ter <- rbind(D042_terFGJ, D042_all)
TW14359_ter <- rbind(TW14359_terFGJ, TW14359_all)
Sakai_ter <- rbind(Sakai_terFGJ, Sakai_all)
EDL933_ter <- rbind(EDL933_terFGJ, EDL933_all)

# change start and end into numeric classes
MG1655_ter <- MG1655_ter %>% 
  mutate(start = as.numeric(start), end = as.numeric(end))

BW2952_ter <- BW2952_ter %>% 
  mutate(start = as.numeric(start), end = as.numeric(end))

REL606_ter <- REL606_ter %>% 
  mutate(start = as.numeric(start), end = as.numeric(end))

APEC078_ter <- APEC078_ter %>% 
  mutate(start = as.numeric(start), end = as.numeric(end))

IAI1_ter <- IAI1_ter %>% 
  mutate(start = as.numeric(start), end = as.numeric(end))

E11368_ter <- E11368_ter %>% 
  mutate(start = as.numeric(start), end = as.numeric(end))

S88_ter <- S88_ter %>% 
  mutate(start = as.numeric(start), end = as.numeric(end))

UTI89_ter <- UTI89_ter %>% 
  mutate(start = as.numeric(start), end = as.numeric(end))

E2348_ter <- E2348_ter %>% 
  mutate(start = as.numeric(start), end = as.numeric(end))

IAI39_ter <- IAI39_ter %>% 
  mutate(start = as.numeric(start), end = as.numeric(end))

SMS35_ter <- SMS35_ter %>% 
  mutate(start = as.numeric(start), end = as.numeric(end))

UMN026_ter <- UMN026_ter %>% 
  mutate(start = as.numeric(start), end = as.numeric(end))

CE10_ter <- CE10_ter %>% 
  mutate(start = as.numeric(start), end = as.numeric(end))

D042_ter <- D042_ter %>% 
  mutate(start = as.numeric(start), end = as.numeric(end))

TW14359_ter <- TW14359_ter %>% 
  mutate(start = as.numeric(start), end = as.numeric(end))

Sakai_ter <- Sakai_ter %>% 
  mutate(start = as.numeric(start), end = as.numeric(end))

EDL933_ter <- EDL933_ter %>% 
  mutate(start = as.numeric(start), end = as.numeric(end))


#________________________________________________________________________________ dna_seg object and plot

list_ter <- list(dna_seg(MG1655_ter),
                 dna_seg(BW2952_ter),
                 dna_seg(REL606_ter),
                 dna_seg(APEC078_ter),
                 dna_seg(IAI1_ter),
                 dna_seg(E11368_ter),
                 dna_seg(S88_ter),
                 dna_seg(UTI89_ter),
                 dna_seg(E2348_ter),
                 dna_seg(IAI39_ter),
                 dna_seg(SMS35_ter),
                 dna_seg(UMN026_ter),
                 dna_seg(CE10_ter),
                 dna_seg(D042_ter),
                 dna_seg(TW14359_ter),
                 dna_seg(Sakai_ter),
                 dna_seg(EDL933_ter))


# assign names to elements in list
names(list_ter) <- c('MG1655',
                      'BW2952',
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

# sanity
list_ter


# outline
# change colours of terFGJ to grey with black outlines
list_ter$MG1655[5][1:3,] <- 'black'
list_ter$BW2952[5][1:3,] <- 'black'
list_ter$REL606[5][1:3,] <- 'black'
list_ter$APEC078[5][1:3,] <- 'black'
list_ter$IAI1[5][1:3,] <- 'black'
list_ter$E11368[5][1:3,] <- 'black'
list_ter$S88[5][1:3,] <- 'black'
list_ter$UTI89[5][1:3,] <- 'black'
list_ter$E2348[5][1:3,] <- 'black'
list_ter$IAI39[5][1:3,] <- 'black'
list_ter$SMS35[5][1:3,] <- 'black'
list_ter$UMN026[5][1:3,] <- 'black'
list_ter$CE10[5][1:3,] <- 'black'
list_ter$D042[5][1:3,] <- 'black'
list_ter$TW14359[5][1:3,] <- 'black'
list_ter$Sakai[5][1:3,] <- 'black'
list_ter$EDL933[5][1:3,] <- 'black'

# change terKLYZ to no outlines and light grey
list_ter$MG1655[5][4:7,] <- 'lightgrey'
list_ter$BW2952[5][4:7,] <- 'lightgrey'
list_ter$REL606[5][4:7,] <- 'lightgrey'
list_ter$APEC078[5][4:7,] <- 'lightgrey'
list_ter$IAI1[5][4:7,] <- 'lightgrey'
list_ter$E11368[5][4:7,] <- 'lightgrey'
list_ter$S88[5][4:7,] <- 'lightgrey'
list_ter$UTI89[5][4:7,] <- 'lightgrey'
list_ter$E2348[5][4:7,] <- 'lightgrey'
list_ter$IAI39[5][4:7,] <- 'lightgrey'
list_ter$SMS35[5][4:7,] <- 'lightgrey'
list_ter$UMN026[5][4:6,] <- 'lightgrey'
list_ter$CE10[5][4:7,] <- 'lightgrey'
list_ter$D042[5][4:7,] <- 'lightgrey'
list_ter$TW14359[5][4:7,] <- 'lightgrey'
list_ter$Sakai[5][4:7,] <- 'lightgrey'
list_ter$EDL933[5][4:7,] <- 'lightgrey'

# fill
# change colours of terFGJ to grey with black outlines
list_ter$MG1655[6][1:3,] <- 'grey'
list_ter$BW2952[6][1:3,] <- 'grey'
list_ter$REL606[6][1:3,] <- 'grey'
list_ter$APEC078[6][1:3,] <- 'grey'
list_ter$IAI1[6][1:3,] <- 'grey'
list_ter$E11368[6][1:3,] <- 'grey'
list_ter$S88[6][1:3,] <- 'grey'
list_ter$UTI89[6][1:3,] <- 'grey'
list_ter$E2348[6][1:3,] <- 'grey'
list_ter$IAI39[6][1:3,] <- 'grey'
list_ter$SMS35[6][1:3,] <- 'grey'
list_ter$UMN026[6][1:3,] <- 'grey'
list_ter$CE10[6][1:3,] <- 'grey'
list_ter$D042[6][1:3,] <- 'grey'
list_ter$TW14359[6][1:3,] <- 'grey'
list_ter$Sakai[6][1:3,] <- 'grey'
list_ter$EDL933[6][1:3,] <- 'grey'

# change terKLYZ to no outlines and light grey
list_ter$MG1655[6][4:7,] <- 'lightgrey'
list_ter$BW2952[6][4:7,] <- 'lightgrey'
list_ter$REL606[6][4:7,] <- 'lightgrey'
list_ter$APEC078[6][4:7,] <- 'lightgrey'
list_ter$IAI1[6][4:7,] <- 'lightgrey'
list_ter$E11368[6][4:7,] <- 'lightgrey'
list_ter$S88[6][4:7,] <- 'lightgrey'
list_ter$UTI89[6][4:7,] <- 'lightgrey'
list_ter$E2348[6][4:7,] <- 'lightgrey'
list_ter$IAI39[6][4:7,] <- 'lightgrey'
list_ter$SMS35[6][4:7,] <- 'lightgrey'
list_ter$UMN026[6][4:6,] <- 'lightgrey'
list_ter$CE10[6][4:7,] <- 'lightgrey'
list_ter$D042[6][4:7,] <- 'lightgrey'
list_ter$TW14359[6][4:7,] <- 'lightgrey'
list_ter$Sakai[6][4:7,] <- 'lightgrey'
list_ter$EDL933[6][4:7,] <- 'lightgrey'

#                                 Trigrobs functions
#________________________________________________________________________________
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


# filter out the strands plus or minus using dplyr
list_ter$MG1655 <- rbind(
  (list_ter$MG1655 %>% 
  filter(strand == 1) %>% 
  mutate(gene_type  = replace(gene_type, gene_type == 'arrows', 'plusGrob'))),
  list_ter$MG1655 %>% 
    filter(strand == -1) %>% 
    mutate(gene_type  = replace(gene_type, gene_type == 'arrows', 'minGrob'))
  )


list_ter$BW2952 <- rbind(
  (list_ter$BW2952 %>% 
     filter(strand == 1) %>% 
     mutate(gene_type  = replace(gene_type, gene_type == 'arrows', 'plusGrob'))),
  list_ter$BW2952 %>% 
    filter(strand == -1) %>% 
    mutate(gene_type  = replace(gene_type, gene_type == 'arrows', 'minGrob'))
)


list_ter$REL606 <- rbind(
  (list_ter$REL606 %>% 
     filter(strand == 1) %>% 
     mutate(gene_type  = replace(gene_type, gene_type == 'arrows', 'plusGrob'))),
  list_ter$REL606 %>% 
    filter(strand == -1) %>% 
    mutate(gene_type  = replace(gene_type, gene_type == 'arrows', 'minGrob'))
)


list_ter$APEC078 <- rbind(
  (list_ter$APEC078 %>% 
     filter(strand == 1) %>% 
     mutate(gene_type  = replace(gene_type, gene_type == 'arrows', 'plusGrob'))),
  list_ter$APEC078 %>% 
    filter(strand == -1) %>% 
    mutate(gene_type  = replace(gene_type, gene_type == 'arrows', 'minGrob'))
)


list_ter$IAI1 <- rbind(
  (list_ter$IAI1 %>% 
     filter(strand == 1) %>% 
     mutate(gene_type  = replace(gene_type, gene_type == 'arrows', 'plusGrob'))),
  list_ter$IAI1 %>% 
    filter(strand == -1) %>% 
    mutate(gene_type  = replace(gene_type, gene_type == 'arrows', 'minGrob'))
)


list_ter$E11368 <- rbind(
  (list_ter$E11368 %>% 
     filter(strand == 1) %>% 
     mutate(gene_type  = replace(gene_type, gene_type == 'arrows', 'plusGrob'))),
  list_ter$E11368 %>% 
    filter(strand == -1) %>% 
    mutate(gene_type  = replace(gene_type, gene_type == 'arrows', 'minGrob'))
)


list_ter$S88 <- rbind(
  (list_ter$S88 %>% 
     filter(strand == 1) %>% 
     mutate(gene_type  = replace(gene_type, gene_type == 'arrows', 'plusGrob'))),
  list_ter$S88 %>% 
    filter(strand == -1) %>% 
    mutate(gene_type  = replace(gene_type, gene_type == 'arrows', 'minGrob'))
)


list_ter$UTI89 <- rbind(
  (list_ter$UTI89 %>% 
     filter(strand == 1) %>% 
     mutate(gene_type  = replace(gene_type, gene_type == 'arrows', 'plusGrob'))),
  list_ter$UTI89 %>% 
    filter(strand == -1) %>% 
    mutate(gene_type  = replace(gene_type, gene_type == 'arrows', 'minGrob'))
)


list_ter$E2348 <- rbind(
  (list_ter$E2348 %>% 
     filter(strand == 1) %>% 
     mutate(gene_type  = replace(gene_type, gene_type == 'arrows', 'plusGrob'))),
  list_ter$E2348 %>% 
    filter(strand == -1) %>% 
    mutate(gene_type  = replace(gene_type, gene_type == 'arrows', 'minGrob'))
)


list_ter$IAI39 <- rbind(
  (list_ter$IAI39 %>% 
     filter(strand == 1) %>% 
     mutate(gene_type  = replace(gene_type, gene_type == 'arrows', 'plusGrob'))),
  list_ter$IAI39 %>% 
    filter(strand == -1) %>% 
    mutate(gene_type  = replace(gene_type, gene_type == 'arrows', 'minGrob'))
)


list_ter$SMS35 <- rbind(
  (list_ter$SMS35 %>% 
     filter(strand == 1) %>% 
     mutate(gene_type  = replace(gene_type, gene_type == 'arrows', 'plusGrob'))),
  list_ter$SMS35 %>% 
    filter(strand == -1) %>% 
    mutate(gene_type  = replace(gene_type, gene_type == 'arrows', 'minGrob'))
)


list_ter$UMN026 <- rbind(
  (list_ter$UMN026 %>% 
     filter(strand == 1) %>% 
     mutate(gene_type  = replace(gene_type, gene_type == 'arrows', 'plusGrob'))),
  list_ter$UMN026 %>% 
    filter(strand == -1) %>% 
    mutate(gene_type  = replace(gene_type, gene_type == 'arrows', 'minGrob'))
)


list_ter$CE10 <- rbind(
  (list_ter$CE10 %>% 
     filter(strand == 1) %>% 
     mutate(gene_type  = replace(gene_type, gene_type == 'arrows', 'plusGrob'))),
  list_ter$CE10 %>% 
    filter(strand == -1) %>% 
    mutate(gene_type  = replace(gene_type, gene_type == 'arrows', 'minGrob'))
)


list_ter$D042 <- rbind(
  (list_ter$D042 %>% 
     filter(strand == 1) %>% 
     mutate(gene_type  = replace(gene_type, gene_type == 'arrows', 'plusGrob'))),
  list_ter$D042 %>% 
    filter(strand == -1) %>% 
    mutate(gene_type  = replace(gene_type, gene_type == 'arrows', 'minGrob'))
)


list_ter$TW14359 <- rbind(
  (list_ter$TW14359 %>% 
     filter(strand == 1) %>% 
     mutate(gene_type  = replace(gene_type, gene_type == 'arrows', 'plusGrob'))),
  list_ter$TW14359 %>% 
    filter(strand == -1) %>% 
    mutate(gene_type  = replace(gene_type, gene_type == 'arrows', 'minGrob'))
)


list_ter$Sakai <- rbind(
  (list_ter$Sakai %>% 
     filter(strand == 1) %>% 
     mutate(gene_type  = replace(gene_type, gene_type == 'arrows', 'plusGrob'))),
  list_ter$Sakai %>% 
    filter(strand == -1) %>% 
    mutate(gene_type  = replace(gene_type, gene_type == 'arrows', 'minGrob'))
)


list_ter$EDL933 <- rbind(
  (list_ter$EDL933 %>% 
     filter(strand == 1) %>% 
     mutate(gene_type  = replace(gene_type, gene_type == 'arrows', 'plusGrob'))),
  list_ter$EDL933 %>% 
    filter(strand == -1) %>% 
    mutate(gene_type  = replace(gene_type, gene_type == 'arrows', 'minGrob'))
)





# reverse the information to get a RFT view
list_ter$MG1655 <- genoPlotR::reverse(list_ter$MG1655)
list_ter$BW2952 <- genoPlotR::reverse(list_ter$BW2952)
list_ter$REL606 <- genoPlotR::reverse(list_ter$REL606)
list_ter$APEC078 <- genoPlotR::reverse(list_ter$APEC078)
list_ter$IAI1 <- genoPlotR::reverse(list_ter$IAI1)
list_ter$E11368 <- genoPlotR::reverse(list_ter$E11368)
list_ter$S88 <- genoPlotR::reverse(list_ter$S88)
list_ter$UTI89 <- genoPlotR::reverse(list_ter$UTI89)
list_ter$E2348 <- genoPlotR::reverse(list_ter$E2348)
list_ter$IAI39 <- genoPlotR::reverse(list_ter$IAI39)
list_ter$SMS35 <- genoPlotR::reverse(list_ter$SMS35)
list_ter$UMN026 <- genoPlotR::reverse(list_ter$UMN026)
list_ter$CE10 <- genoPlotR::reverse(list_ter$CE10)
list_ter$D042 <- genoPlotR::reverse(list_ter$D042)
list_ter$TW14359 <- genoPlotR::reverse(list_ter$TW14359)
list_ter$Sakai <- genoPlotR::reverse(list_ter$Sakai)
list_ter$EDL933 <- genoPlotR::reverse(list_ter$EDL933)

# keep the strand as in the original


# sanity
list_ter


## Calculating middle positions
mid_pos1 <- middle(list_ter[[1]])
mid_pos2 <- middle(list_ter[[2]])
mid_pos3 <- middle(list_ter[[3]])
mid_pos4 <- middle(list_ter[[4]])
mid_pos5 <- middle(list_ter[[5]])
mid_pos6 <- middle(list_ter[[6]])
mid_pos7 <- middle(list_ter[[7]])
mid_pos8 <- middle(list_ter[[8]])
mid_pos9 <- middle(list_ter[[9]])
mid_pos10 <- middle(list_ter[[10]])
mid_pos11 <- middle(list_ter[[11]])
mid_pos12 <- middle(list_ter[[12]])
mid_pos13 <- middle(list_ter[[13]])
mid_pos14 <- middle(list_ter[[14]])
mid_pos15 <- middle(list_ter[[15]])
mid_pos16 <- middle(list_ter[[16]])
mid_pos17 <- middle(list_ter[[17]])



# Annotations
annot1 <- genoPlotR::annotation(x1 = mid_pos1, text = list_ter[[1]]$name, rot = 90, col = 'black')
annot2 <- genoPlotR::annotation(x1 = mid_pos2, text = list_ter[[2]]$name, rot = 90, col = 'black')
annot3 <- genoPlotR::annotation(x1 = mid_pos3, text = list_ter[[3]]$name, rot = 90, col = 'black')
annot4 <- genoPlotR::annotation(x1 = mid_pos4, text = list_ter[[4]]$name, rot = 90, col = 'black')
annot5 <- genoPlotR::annotation(x1 = mid_pos5, text = list_ter[[5]]$name, rot = 90, col = 'black')
annot6 <- genoPlotR::annotation(x1 = mid_pos6, text = list_ter[[6]]$name, rot = 90, col = 'black')
annot7 <- genoPlotR::annotation(x1 = mid_pos7, text = list_ter[[7]]$name, rot = 90, col = 'black')
annot8 <- genoPlotR::annotation(x1 = mid_pos8, text = list_ter[[8]]$name, rot = 90, col = 'black')
annot9 <- genoPlotR::annotation(x1 = mid_pos9, text = list_ter[[9]]$name, rot = 90, col = 'black')
annot10 <- genoPlotR::annotation(x1 = mid_pos10, text = list_ter[[10]]$name, rot = 90, col = 'black')
annot11 <- genoPlotR::annotation(x1 = mid_pos11, text = list_ter[[11]]$name, rot = 90, col = 'black')
annot12 <- genoPlotR::annotation(x1 = mid_pos12, text = list_ter[[12]]$name, rot = 90, col = 'black')
annot13 <- genoPlotR::annotation(x1 = mid_pos13, text = list_ter[[13]]$name, rot = 90, col = 'black')
annot14 <- genoPlotR::annotation(x1 = mid_pos14, text = list_ter[[14]]$name, rot = 90, col = 'black')
annot15 <- genoPlotR::annotation(x1 = mid_pos15, text = list_ter[[15]]$name, rot = 90, col = 'black')
annot16 <- genoPlotR::annotation(x1 = mid_pos16, text = list_ter[[16]]$name, rot = 90, col = 'black')
annot17 <- genoPlotR::annotation(x1 = mid_pos17, text = list_ter[[17]]$name, rot = 90, col = 'black')

#list of annotation objects
annots <- list(annot1, 
                 annot2, 
                 annot3, 
                 annot4, 
                 annot5, 
                 annot6, 
                 annot7, 
                 annot8, 
                 annot9, 
                 annot10, 
                 annot11, 
                 annot12, 
                 annot13, 
                 annot14, 
                 annot15, 
                 annot16, 
                 annot17)



# plot as trigrobs
plot_gene_map(dna_segs = list_ter,
              annotations = annots,
              comparisons = NULL,
              scale = FALSE,
              main = 'Secondary ter sites with terFGJ',
              annotation_cex = 0.5,
              dna_seg_line = TRUE,
              dna_seg_scale = FALSE)




# change pTer into side_blocks
list_ter$MG1655[4:7,11] <- 'side_blocks'
list_ter$BW2952[4:7,11] <- 'side_blocks'
list_ter$REL606[4:7,11] <- 'side_blocks'
list_ter$APEC078[4:7,11] <- 'side_blocks'
list_ter$IAI1[4:7,11] <- 'side_blocks'
list_ter$E11368[4:7,11] <- 'side_blocks'
list_ter$S88[4:7,11] <- 'side_blocks'
list_ter$UTI89[4:7,11] <- 'side_blocks'
list_ter$E2348[4:7,11] <- 'side_blocks'
list_ter$IAI39[4:7,11] <- 'side_blocks'
list_ter$SMS35[4:7,11] <- 'side_blocks'
list_ter$UMN026[4:6,11] <- 'side_blocks'
list_ter$CE10[4:7,11] <- 'side_blocks'
list_ter$D042[4:7,11] <- 'side_blocks'
list_ter$TW14359[4:7,11] <- 'side_blocks'
list_ter$Sakai[4:7,11] <- 'side_blocks'
list_ter$EDL933[4:7,11] <- 'side_blocks'


# change terKLYZ to blue to make it easier to track
list_ter$MG1655[5][4:7,] <- 'red3'

list_ter$BW2952[5][4:7,] <- 'red3'
list_ter$REL606[5][4:7,] <- 'red3'
list_ter$APEC078[5][4:7,] <- 'red3'
list_ter$IAI1[5][4:7,] <- 'red3'
list_ter$E11368[5][4:7,] <- 'red3'
list_ter$S88[5][4:7,] <- 'red3'
list_ter$UTI89[5][4:7,] <- 'red3'
list_ter$E2348[5][4:7,] <- 'red3'
list_ter$IAI39[5][4:7,] <- 'red3'
list_ter$SMS35[5][4:7,] <- 'red3'
list_ter$UMN026[5][4:6,] <- 'red3'
list_ter$CE10[5][4:7,] <- 'red3'
list_ter$D042[5][4:7,] <- 'red3'
list_ter$TW14359[5][4:7,] <- 'red3'
list_ter$Sakai[5][4:7,] <- 'red3'
list_ter$EDL933[5][4:7,] <- 'red3'

list_ter$E2348
# plot as terFGJ as trigrobs and pTer as side_blocks
plot_gene_map(dna_segs = list_ter,
              annotations = annots,
              comparisons = NULL,
              scale = FALSE,
              main = 'Secondary ter sites with terFGJ in Permissible Orientation',
              annotation_cex = 0.5,
              dna_seg_label_cex = 0.8,
              annotation_height = 1,
              dna_seg_line = TRUE,
              dna_seg_scale = FALSE)







