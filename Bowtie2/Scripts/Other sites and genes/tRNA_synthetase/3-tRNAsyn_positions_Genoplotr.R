# 6. Parse through genoplotR to map their positions then compare
# 7. PROFIT

library(Biostrings)
library(seqinr)
library(genoPlotR)
library(tibble)

# set working directory before proceeding
file_in_dir <- file.choose() #choose a file in the desired dir
dir <- dirname(file_in_dir) #get the directory of where that file is
setwd(dir) #set the dir
getwd() #sanity

# read in the csv from BT2 and samtools
df <- read.csv('tRNA_Ecoli.csv')
df #sanity

#             Clean the tRNA synthetase csv 
#           Cols = rname[2], qname[3], pos[4], strand[6]
#                 (Strain)    (Gene)      
#___________________________________________________________________________

seqs <- df[c(2,3,4,6)]
colnames(seqs) <- c('strain','name', 'start', 'strand')
seqs #sanity
toString(seqs$strain)

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
D042 <- seqs[which(seqs$strain == '42'),]
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
list_strains <- lapply(list_strains, function(x) add_column(x, end = x$start + 50000, .after = 3))

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

# different colours for each name (tRNA gene)
# random colour generator
col = sample(colours(), 23)
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
              main = 'tRNA Synthetase Positions in E.coli Phylogroups',
              annotation_height = 3,
              dna_seg_scale=FALSE,
              dna_seg_label_cex=0.8,
              scale = FALSE,
              arrow_head_len = 600000,
              minimum_gap_size = 0.1)










