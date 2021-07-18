# 6. Parse through genoplotR to map their positions then compare
# 7. PROFIT

library(Biostrings)
library(seqinr)
library(genoPlotR)
library(tibble)
library(tidyverse)
library(stringr)
library(DataCombine) # for grepl.sub to subset df on ribo proteins

# set working directory before proceeding
file_in_dir <- file.choose() #choose a file in the desired dir
dir <- dirname(file_in_dir) #get the directory of where that file is
setwd(dir) #set the dir
getwd() #sanity

# read in the csv from BT2 and samtools
df <- read.csv('ribosome_Ecoli.csv')
df #sanity

#             Clean the tRNA synthetase csv 
#           Cols = rname[2], qname[3], pos[4], strand[6]
#                 (Strain)    (Gene)      
#___________________________________________________________________________

seqs <- df[c(2,3,4,6)]
colnames(seqs) <- c('strain','name', 'start', 'strand')
seqs #sanity
toString(seqs$strain) # show all the strain names

# subset the large df into each strain
# then later test with just the 5 phylogroups
#A
MG1655 <- seqs[which(seqs$strain == 'MG1655'),]
BW2952 <- seqs[which(seqs$strain == 'BW2952'),]
REL606 <- seqs[which(seqs$strain == 'REL606'),]
#B1
APEC078 <- seqs[which(seqs$strain == 'APECO78'),]
IAI1 <- seqs[which(seqs$strain == 'IAI1'),]
E11368 <- seqs[which(seqs$strain == '11368'),]
#B2
S88 <- seqs[which(seqs$strain == 'S88'),]
UTI89 <- seqs[which(seqs$strain == 'UTI89'),]
E2348 <- seqs[which(seqs$strain == 'E2348/69'),]
#D
IAI39 <- seqs[which(seqs$strain == 'IAI39'),]
SMS35 <- seqs[which(seqs$strain == 'SMS-3-5'),]
UMN026 <- seqs[which(seqs$strain == 'UMN026'),]
CE10 <- seqs[which(seqs$strain == 'CE10'),]
D042 <- seqs[which(seqs$strain == '042'),]
#E
TW14359 <- seqs[which(seqs$strain == 'TW14359'),]
Sakai <- seqs[which(seqs$strain == 'Sakai'),]
EDL933 <- seqs[which(seqs$strain == 'EDL933'),]

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

##########################################
#
#   NEW WAY TO CLEAN DATAFRAMES WITH LAPPLY
#               16/04/21
##########################################

# Drop the Strain column and add the end column
# without specifying [-1] for every df
list_strains <- list(
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
names(list_strains) <- c(
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

# add another column conntaining the end point
# using tibble library
list_strains <- lapply(list_strains, function(x) add_column(x, end = x$start + 500, .after = 3))

# drop the 1st column in all
list_strains <- lapply(list_strains, function(x) x[names(x) != 'strain'])

#sanity
list_strains

#                 GenoplotR and dna_seg objects
#_________________________________________________________________

MG_dna <- dna_seg(list_strains$MG1655)
AP_dna <- dna_seg(list_strains$APEC078)
S8_dna <- dna_seg(list_strains$S88)
UM_dna <- dna_seg(list_strains$UMN026)
TW_dna <- dna_seg(list_strains$TW14359)

#sanity
MG_dna
AP_dna
S8_dna
UM_dna
TW_dna

#dna_seg list
list_dna <- list(
  MG_dna,
  AP_dna,
  S8_dna,
  UM_dna,
  TW_dna
)

# use lapply to sort list_dna based on name
# this will allow the colours to apply to all properly
list_dna <- lapply(list_dna, function(y){
  y[order(y$name),]
})

#sanity
list_dna

# name
names(list_dna) <- c(
  'MG1655',
  'APEC078',
  'S88',
  'UMN026',
  'TW14359')

#sanity
list_dna

# how many genes in each list
length(list_dna$MG1655)

# different colours for each name (tRNA gene)
# random colour generator
col = sample(colours(), 50)
col #sanity

# manually change the fill to the randomly generated colours
# can change colours later if they are not pretty enough.. do not worry about them
list_dna$MG1655[c(5,6)] <- col
list_dna$APEC078[c(5,6)] <- col
list_dna$S88[c(5,6)] <- col
list_dna$UMN026[c(5,6)] <- col
list_dna$TW14359[c(5,6)] <- col

#final sanity before creating the plot
list_dna

#               Specify plot 
#________________________________________________________

# middle positions
mid1 <- middle(list_dna[[1]])
mid2 <- middle(list_dna[[2]])
mid3 <- middle(list_dna[[3]])
mid4 <- middle(list_dna[[4]])
mid5 <- middle(list_dna[[5]])

# annotations
annot1 <- annotation(x1 = mid1, text = list_dna[[1]]$name, rot = 90, col = col)
annot2 <- annotation(x1 = mid2, text = list_dna[[2]]$name, rot = 90, col = col)
annot3 <- annotation(x1 = mid3, text = list_dna[[3]]$name, rot = 90, col = col)
annot4 <- annotation(x1 = mid4, text = list_dna[[4]]$name, rot = 90, col = col)
annot5 <- annotation(x1 = mid5, text = list_dna[[5]]$name, rot = 90, col = col)
annot_list <- list(annot1, annot2, annot3, annot4, annot5)

# plot
plot_gene_map(dna_segs = list_dna, 
              comparisons = NULL,
              annotations = annot_list,
              main = 'Ribosomal Protein Positions in E.coli Phylogroups',
              annotation_height = 3,
              dna_seg_scale=FALSE,
              dna_seg_label_cex=0.8,
              scale = FALSE,
              minimum_gap_size = 0.1, 
              gene_type = 'side_blocks')


#________________________________________________________________________________
#                                 50 proteins  is too dense
#                                 Separate into 3 figures
#                               rpl, rpm, rps ribosome proteins
#________________________________________________________________________________

#                         Example from 
# https://www.rdocumentation.org/packages/DataCombine/versions/0.2.21/topics/grepl.sub

# use DataCombine package to use grepl.sub to extract classes of ribo proteins
rpl <- grepl.sub(data = df, pattern = 'rpl', Var = 'qname')
rpm <- grepl.sub(data = df, pattern = 'rpm', Var = 'qname')
rps <- grepl.sub(data = df, pattern = 'rps', Var = 'qname')

#sanity
# this is super powerful
rpl
rpm
rps

# only keep strain name start and strand information
rpl <- rpl[c(2,3,4,6)]
colnames(rpl) <- c('strain','name', 'start', 'strand')
rpl

rpm <- rpm[c(2,3,4,6)]
colnames(rpm) <- c('strain','name', 'start', 'strand')
rpm

rps <- rps[c(2,3,4,6)]
colnames(rps) <- c('strain','name', 'start', 'strand')
rps

###############################################
#
#                     rpl
#
###############################################

# take the strains out of each ribo df rpl, rpm and rps
# subset the large df into each strain
# then later test with just the 5 phylogroups
#A
MG1655_rpl <- rpl[which(rpl$strain == 'MG1655'),]
BW2952_rpl <- rpl[which(rpl$strain == 'BW2952'),]
REL606_rpl <- rpl[which(rpl$strain == 'REL606'),]
#B1
APEC078_rpl <- rpl[which(rpl$strain == 'APECO78'),]
IAI1_rpl <- rpl[which(rpl$strain == 'IAI1'),]
E11368_rpl <- rpl[which(rpl$strain == '11368'),]
#B2
S88_rpl <- rpl[which(rpl$strain == 'S88'),]
UTI89_rpl <- rpl[which(rpl$strain == 'UTI89'),]
E2348_rpl <- rpl[which(rpl$strain == 'E2348/69'),]
#D
IAI39_rpl <- rpl[which(rpl$strain == 'IAI39'),]
SMS35_rpl <- rpl[which(rpl$strain == 'SMS-3-5'),]
UMN026_rpl <- rpl[which(rpl$strain == 'UMN026'),]
CE10_rpl <- rpl[which(rpl$strain == 'CE10'),]
D042_rpl <- rpl[which(rpl$strain == '042'),]
#E
TW14359_rpl <- rpl[which(rpl$strain == 'TW14359'),]
Sakai_rpl <- rpl[which(rpl$strain == 'Sakai'),]
EDL933_rpl <- rpl[which(rpl$strain == 'EDL933'),]

# SANIITY
#A
MG1655_rpl
BW2952_rpl
REL606_rpl
#B1
APEC078_rpl
IAI1_rpl
E11368_rpl
#B2
S88_rpl
UTI89_rpl
E2348_rpl
#D
IAI39_rpl
SMS35_rpl
UMN026_rpl
CE10_rpl
D042_rpl
#E
TW14359_rpl
Sakai_rpl
EDL933_rpl

# rpl 5 phylogroups
rpl_list <- list(MG1655_rpl,
                APEC078_rpl,
                S88_rpl,
                UMN026_rpl,
                TW14359_rpl)
names(rpl_list) <-  c('MG1655',
'APEC078',
'S88',
'UMN026',
'TW14359')

rpl_list

# add another column conntaining the end point
# using tibble library
rpl_list <- lapply(rpl_list, function(x) add_column(x, end = x$start + 30000, .after = 3))

# drop the 1st column in all
rpl_list <- lapply(rpl_list, function(x) x[names(x) != 'strain'])

#sanity
list_strains

#                 GenoplotR and dna_seg objects
#_________________________________________________________________

MGrpl_dna <- dna_seg(rpl_list$MG1655)
APrpl_dna <- dna_seg(rpl_list$APEC078)
S8rpl_dna <- dna_seg(rpl_list$S88)
UMrpl_dna <- dna_seg(rpl_list$UMN026)
TWrpl_dna <- dna_seg(rpl_list$TW14359)

#sanity
MGrpl_dna
APrpl_dna
S8rpl_dna
UMrpl_dna
TWrpl_dna

#dna_seg list
list_dna_rpl <- list(
  MGrpl_dna,
  APrpl_dna,
  S8rpl_dna,
  UMrpl_dna,
  TWrpl_dna
)

# use lapply to sort list_dna based on name
# this will allow the colours to apply to all properly
list_dna_rpl <- lapply(list_dna_rpl, function(y){
  y[order(y$name),]
})

#sanity
list_dna_rpl

# name
names(list_dna_rpl) <- c(
  'MG1655',
  'APEC078',
  'S88',
  'UMN026',
  'TW14359')

#sanity
list_dna_rpl


# different colours for each name (tRNA gene)
# random colour generator
col_rpl = sample(colours(), 23)
col_rpl #sanity

# manually change the fill to the randomly generated colours
# can change colours later if they are not pretty enough.. do not worry about them
list_dna_rpl$MG1655[c(5,6)] <- col_rpl
list_dna_rpl$APEC078[c(5,6)] <- col_rpl
list_dna_rpl$S88[c(5,6)] <- col_rpl
list_dna_rpl$UMN026[c(5,6)] <- col_rpl
list_dna_rpl$TW14359[c(5,6)] <- col_rpl

#final sanity before creating the plot
list_dna_rpl

#               Specify plot 
#________________________________________________________

# middle positions
mid1rpl <- middle(list_dna_rpl[[1]])
mid2rpl <- middle(list_dna_rpl[[2]])
mid3rpl <- middle(list_dna_rpl[[3]])
mid4rpl <- middle(list_dna_rpl[[4]])
mid5rpl <- middle(list_dna_rpl[[5]])

# annotations
annot1rpl <- annotation(x1 = mid1rpl, text = list_dna_rpl[[1]]$name, rot = 90, col = col_rpl)
annot2rpl <- annotation(x1 = mid2rpl, text = list_dna_rpl[[2]]$name, rot = 90, col = col_rpl)
annot3rpl <- annotation(x1 = mid3rpl, text = list_dna_rpl[[3]]$name, rot = 90, col = col_rpl)
annot4rpl <- annotation(x1 = mid4rpl, text = list_dna_rpl[[4]]$name, rot = 90, col = col_rpl)
annot5rpl <- annotation(x1 = mid5rpl, text = list_dna_rpl[[5]]$name, rot = 90, col = col_rpl)
annot_listrpl <- list(annot1rpl, annot2rpl, annot3rpl, annot4rpl, annot5rpl)

# plot
plot_gene_map(dna_segs = list_dna_rpl, 
              comparisons = NULL,
              annotations = annot_listrpl,
              main = 'Ribosomal Protein (rpl) Positions in E.coli Phylogroups',
              annotation_height = 3,
              dna_seg_scale=FALSE,
              dna_seg_label_cex=0.8,
              scale = FALSE,
              minimum_gap_size = 0.1,
              gene_type = 'arrows',
              arrow_head_len = 100000)


###############################################
#
#                     rpm
#
###############################################

# take the strains out of each ribo df rpl, rpm and rps
# subset the large df into each strain
# then later test with just the 5 phylogroups

#A
MG1655_rpm <- rpm[which(rpm$strain == 'MG1655'),]
BW2952_rpm <- rpm[which(rpm$strain == 'BW2952'),]
REL606_rpm <- rpm[which(rpm$strain == 'REL606'),]
#B1
APEC078_rpm <- rpm[which(rpm$strain == 'APECO78'),]
IAI1_rpm <- rpm[which(rpm$strain == 'IAI1'),]
E11368_rpm <- rpm[which(rpm$strain == '11368'),]
#B2
S88_rpm <- rpm[which(rpm$strain == 'S88'),]
UTI89_rpm <- rpm[which(rpm$strain == 'UTI89'),]
E2348_rpm <- rpm[which(rpm$strain == 'E2348/69'),]
#D
IAI39_rpm <- rpm[which(rpm$strain == 'IAI39'),]
SMS35_rpm <- rpm[which(rpm$strain == 'SMS-3-5'),]
UMN026_rpm <- rpm[which(rpm$strain == 'UMN026'),]
CE10_rpm <- rpm[which(rpm$strain == 'CE10'),]
D042_rpm <- rpm[which(rpm$strain == '042'),]
#E
TW14359_rpm <- rpm[which(rpm$strain == 'TW14359'),]
Sakai_rpm <- rpm[which(rpm$strain == 'Sakai'),]
EDL933_rpm <- rpm[which(rpm$strain == 'EDL933'),]

# SANIITY
#A
MG1655_rpm
BW2952_rpm
REL606_rpm
#B1
APEC078_rpm
IAI1_rpm
E11368_rpm
#B2
S88_rpm
UTI89_rpm
E2348_rpm
#D
IAI39_rpm
SMS35_rpm
UMN026_rpm
CE10_rpm
D042_rpm
#E
TW14359_rpm
Sakai_rpm
EDL933_rpm

# rpl 5 phylogroups
rpm_list <- list(MG1655_rpm,
                 APEC078_rpm,
                 S88_rpm,
                 UMN026_rpm,
                 TW14359_rpm)
names(rpm_list) <-  c('MG1655',
                      'APEC078',
                      'S88',
                      'UMN026',
                      'TW14359')

rpm_list

# add another column conntaining the end point
# using tibble library
rpm_list <- lapply(rpm_list, function(x) add_column(x, end = x$start + 30000, .after = 3))

# drop the 1st column in all
rpm_list <- lapply(rpm_list, function(x) x[names(x) != 'strain'])

#sanity
rpm_list

#                 GenoplotR and dna_seg objects
#_________________________________________________________________

MGrpm_dna <- dna_seg(rpm_list$MG1655)
APrpm_dna <- dna_seg(rpm_list$APEC078)
S8rpm_dna <- dna_seg(rpm_list$S88)
UMrpm_dna <- dna_seg(rpm_list$UMN026)
TWrpm_dna <- dna_seg(rpm_list$TW14359)

#sanity
MGrpm_dna
APrpm_dna
S8rpm_dna
UMrpm_dna
TWrpm_dna

#dna_seg list
list_dna_rpm <- list(
  MGrpm_dna,
  APrpm_dna,
  S8rpm_dna,
  UMrpm_dna,
  TWrpm_dna
)

# use lapply to sort list_dna based on name
# this will allow the colours to apply to all properly
list_dna_rpm <- lapply(list_dna_rpm, function(y){
  y[order(y$name),]
})

#sanity
list_dna_rpm

# name
names(list_dna_rpm) <- c(
  'MG1655',
  'APEC078',
  'S88',
  'UMN026',
  'TW14359')


#sanity
list_dna_rpm


# different colours for each name (tRNA gene)
# random colour generator
col_rpm = sample(colours(), 10)
col_rpm #sanity

# manually change the fill to the randomly generated colours
# can change colours later if they are not pretty enough.. do not worry about them
list_dna_rpm$MG1655[c(5,6)] <- col_rpm
list_dna_rpm$APEC078[c(5,6)] <- col_rpm
list_dna_rpm$S88[c(5,6)] <- col_rpm
list_dna_rpm$UMN026[c(5,6)] <- col_rpm
list_dna_rpm$TW14359[c(5,6)] <- col_rpm

#final sanity before creating the plot
list_dna_rpm

#               Specify plot 
#________________________________________________________

# middle positions
mid1rpm <- middle(list_dna_rpm[[1]])
mid2rpm <- middle(list_dna_rpm[[2]])
mid3rpm <- middle(list_dna_rpm[[3]])
mid4rpm <- middle(list_dna_rpm[[4]])
mid5rpm <- middle(list_dna_rpm[[5]])

# annotations
annot1rpm <- annotation(x1 = mid1rpm, text = list_dna_rpm[[1]]$name, rot = 90, col = col_rpm)
annot2rpm <- annotation(x1 = mid2rpm, text = list_dna_rpm[[2]]$name, rot = 90, col = col_rpm)
annot3rpm <- annotation(x1 = mid3rpm, text = list_dna_rpm[[3]]$name, rot = 90, col = col_rpm)
annot4rpm <- annotation(x1 = mid4rpm, text = list_dna_rpm[[4]]$name, rot = 90, col = col_rpm)
annot5rpm <- annotation(x1 = mid5rpm, text = list_dna_rpm[[5]]$name, rot = 90, col = col_rpm)
annot_listrpm <- list(annot1rpm, annot2rpm, annot3rpm, annot4rpm, annot5rpm)

# plot
plot_gene_map(dna_segs = list_dna_rpm, 
              comparisons = NULL,
              annotations = annot_listrpm,
              main = 'Ribosomal Protein (rpm) Positions in E.coli Phylogroups',
              annotation_height = 3,
              dna_seg_scale=FALSE,
              dna_seg_label_cex=0.8,
              scale = FALSE,
              minimum_gap_size = 0.1, 
              gene_type = 'arrows',
              arrow_head_len = 100000)






###############################################
#
#                     rps
#
###############################################

# take the strains out of each ribo df rpl, rpm and rps
# subset the large df into each strain
# then later test with just the 5 phylogroups

#A
MG1655_rps <- rps[which(rps$strain == 'MG1655'),]
BW2952_rps <- rps[which(rps$strain == 'BW2952'),]
REL606_rps <- rps[which(rps$strain == 'REL606'),]
#B1
APEC078_rps <- rps[which(rps$strain == 'APECO78'),]
IAI1_rps <- rps[which(rps$strain == 'IAI1'),]
E11368_rps <- rps[which(rps$strain == '11368'),]
#B2
S88_rps <- rps[which(rps$strain == 'S88'),]
UTI89_rps <- rps[which(rps$strain == 'UTI89'),]
E2348_rps <- rps[which(rps$strain == 'E2348/69'),]
#D
IAI39_rps <- rps[which(rps$strain == 'IAI39'),]
SMS35_rps <- rps[which(rps$strain == 'SMS-3-5'),]
UMN026_rps <- rps[which(rps$strain == 'UMN026'),]
CE10_rps <- rps[which(rps$strain == 'CE10'),]
D042_rps <- rps[which(rps$strain == '042'),]
#E
TW14359_rps <- rps[which(rps$strain == 'TW14359'),]
Sakai_rps <- rps[which(rps$strain == 'Sakai'),]
EDL933_rps <- rps[which(rps$strain == 'EDL933'),]

# SANIITY
#A
MG1655_rps
BW2952_rps
REL606_rps
#B1
APEC078_rps
IAI1_rps
E11368_rps
#B2
S88_rps
UTI89_rps
E2348_rps
#D
IAI39_rps
SMS35_rps
UMN026_rps
CE10_rps
D042_rps
#E
TW14359_rps
Sakai_rps
EDL933_rps

# rpl 5 phylogroups
rps_list <- list(MG1655_rps,
                 APEC078_rps,
                 S88_rps,
                 UMN026_rps,
                 TW14359_rps)
names(rps_list) <-  c('MG1655',
                      'APEC078',
                      'S88',
                      'UMN026',
                      'TW14359')

rps_list

# add another column conntaining the end point
# using tibble library
rps_list <- lapply(rps_list, function(x) add_column(x, end = x$start + 30000, .after = 3))

# drop the 1st column in all
rps_list <- lapply(rps_list, function(x) x[names(x) != 'strain'])

#sanity
rps_list

#                 GenoplotR and dna_seg objects
#_________________________________________________________________

MGrps_dna <- dna_seg(rps_list$MG1655)
APrps_dna <- dna_seg(rps_list$APEC078)
S8rps_dna <- dna_seg(rps_list$S88)
UMrps_dna <- dna_seg(rps_list$UMN026)
TWrps_dna <- dna_seg(rps_list$TW14359)

#sanity
MGrps_dna
APrps_dna
S8rps_dna
UMrps_dna
TWrps_dna

#dna_seg list
list_dna_rps <- list(
  MGrps_dna,
  APrps_dna,
  S8rps_dna,
  UMrps_dna,
  TWrps_dna
)

# use lapply to sort list_dna based on name
# this will allow the colours to apply to all properly
list_dna_rps <- lapply(list_dna_rps, function(y){
  y[order(y$name),]
})

#sanity
list_dna_rps

# name
names(list_dna_rps) <- c(
  'MG1655',
  'APEC078',
  'S88',
  'UMN026',
  'TW14359')


#sanity
list_dna_rps


# different colours for each name (tRNA gene)
# random colour generator
col_rps = sample(colours(), 17)
col_rps #sanity

# manually change the fill to the randomly generated colours
# can change colours later if they are not pretty enough.. do not worry about them
list_dna_rps$MG1655[c(5,6)] <- col_rps
list_dna_rps$APEC078[c(5,6)] <- col_rps
list_dna_rps$S88[c(5,6)] <- col_rps
list_dna_rps$UMN026[c(5,6)] <- col_rps
list_dna_rps$TW14359[c(5,6)] <- col_rps

#final sanity before creating the plot
list_dna_rps

#               Specify plot 
#________________________________________________________

# middle positions
mid1rps <- middle(list_dna_rps[[1]])
mid2rps <- middle(list_dna_rps[[2]])
mid3rps <- middle(list_dna_rps[[3]])
mid4rps <- middle(list_dna_rps[[4]])
mid5rps <- middle(list_dna_rps[[5]])

# annotations
annot1rps <- annotation(x1 = mid1rps, text = list_dna_rps[[1]]$name, rot = 90, col = col_rps)
annot2rps <- annotation(x1 = mid2rps, text = list_dna_rps[[2]]$name, rot = 90, col = col_rps)
annot3rps <- annotation(x1 = mid3rps, text = list_dna_rps[[3]]$name, rot = 90, col = col_rps)
annot4rps <- annotation(x1 = mid4rps, text = list_dna_rps[[4]]$name, rot = 90, col = col_rps)
annot5rps <- annotation(x1 = mid5rps, text = list_dna_rps[[5]]$name, rot = 90, col = col_rps)
annot_listrps <- list(annot1rps, annot2rps, annot3rps, annot4rps, annot5rps)

# plot
plot_gene_map(dna_segs = list_dna_rps, 
              comparisons = NULL,
              annotations = annot_listrps,
              main = 'Ribosomal Protein (rps) Positions in E.coli Phylogroups',
              annotation_height = 3,
              dna_seg_scale=FALSE,
              dna_seg_label_cex=0.8,
              scale = FALSE,
              minimum_gap_size = 0.1, 
              gene_type = 'arrows',
              arrow_head_len = 100000)
