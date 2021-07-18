# purpose: find and align the terY sites in all 17 E.coli genomes to see variability
# date created: 18/05/21
# used with other scripts: No

library(Biostrings)
library(msa)
library(IRanges)
library(tidyverse)

#_________________________________________________________________________________ Set terY sequence and load in E.coli
Ecoli <- readDNAStringSet(file.choose())


terY <- DNAString('TATGGGTACGTTGTAATTAGGGA')


# match
match_pos <- vmatchPattern(terY, Ecoli, max.mismatch = 5)
match_min <- vmatchPattern(reverseComplement(terY), Ecoli, max.mismatch = 5)

# place the IRanges into unlisted IRange objects depending on their mismatch value
mismatch_pos0 <- unlist(match_pos)
mismatch_min0 <- unlist(match_min)

# only 5 mismatches gets any further hits
mismatch_pos5 <- unlist(match_pos)
mismatch_min5 <- unlist(match_min)

# sanity
mismatch_pos0
mismatch_min0

mismatch_pos5
mismatch_min5

reduce(c(mismatch_min0, mismatch_pos5, mismatch_min5))

#___________________________________________________________________________________ Send to df
#                                                                                     And specify strands
df0 <- data.frame(mismatch_min0)
df0$strand <- '-'
df0



temp1 <- data.frame(mismatch_pos5)
temp2 <- data.frame(mismatch_min5)
temp1$strand <- '+'
temp2$strand <- '-'
df5 <- data.frame(rbind(temp1, temp2))
df5

# merge into one df, dropping duplicated strains 
# then send to csv for the final trigrob figure where we add the terY site
merged_df <- rbind(df0, df5)
merged_df <- temp[!duplicated(temp$names),]

# genome order

gen_order <- c('MG1655',
'BW2952',
'REL606',
'APECO78',      
'IAI1',
'11368',
'S88',
'UTI89',
'E2348/69',
'IAI39',
'SMS-3-5',
'CE10',
'042',
'TW14359',
'Sakai',
'EDL933')

# reorder the names to gen_order so when we export as csv its in order already
rownames(merged_df) <- NULL

merged_df <- merged_df[c(1,2,3,4,5,6,7,8,9,10,14,15,16,11,12,13),]

merged_df # sanity

setwd(dirname(file.choose()))
write.csv(merged_df, 'terY_Bioc_Ecoli.csv', row.names = FALSE)

#___________________________________________________________________________________ Subseq the ter sequences 

# not as easy as it looks as the df names and Ecoli set do not line up
MG_terY <- subseq(Ecoli$MG1655, start = df0[which(df0$names=='MG1655'),]$start, width = 23)
MG_terY

# SANITY
#A
MG1655  <- subseq(Ecoli$MG1655, start = df0[which(df0$names=='MG1655'),]$start, width = 23)
BW2952  <- subseq(Ecoli$BW2952, start = df0[which(df0$names=='BW2952'),]$start, width = 23)
REL606  <- subseq(Ecoli$REL606, start = df0[which(df0$names=='REL606'),]$start, width = 23)
#B1
APEC078 <- subseq(Ecoli$APECO78, start = df0[which(df0$names=='APECO78'),]$start, width = 23)
IAI1    <- subseq(Ecoli$IAI1, start = df0[which(df0$names=='IAI1'),]$start, width = 23)
E11368  <- subseq(Ecoli$`11368`, start = df0[which(df0$names=='11368'),]$start, width = 23)
#B2
S88   <- subseq(Ecoli$S88, start = df0[which(df0$names=='S88'),]$start, width = 23)
UTI89 <- subseq(Ecoli$UTI89, start = df0[which(df0$names=='UTI89'),]$start, width = 23)
E2348 <- subseq(Ecoli$`E2348/69 `, start = df0[which(df0$names=='E2348/69 '),]$start, width = 23)
#D
IAI39   <- subseq(Ecoli$IAI39, start = df5[which(df5$names=='IAI39'),]$start, width = 23)
SMS35   <- subseq(Ecoli$`SMS-3-5`, start = df5[which(df5$names=='SMS-3-5'),]$start, width = 23)
UMN026  # No terY
CE10    <- subseq(Ecoli$CE10, start = df5[which(df5$names=='CE10'),]$start, width = 23)
D042    <- subseq(Ecoli$`042`, start = df0[which(df0$names=='042'),]$start, width = 23)
#E
TW14359 <- subseq(Ecoli$TW14359, start = df5[which(df5$names=='TW14359'),]$start, width = 23)
Sakai   <- subseq(Ecoli$Sakai, start = df5[which(df5$names=='Sakai'),]$start, width = 23)
EDL933  <- subseq(Ecoli$EDL933, start = df5[which(df5$names=='EDL933'),]$start, width = 23)



# sanity for terY with 0 mismatches
MG1655
BW2952
REL606
APEC078
IAI1
E11368
S88
UTI89
E2348
D042

# sanity for terY with 5 mismatches
TW14359
Sakai
EDL933
IAI39
SMS35
CE10



#___________________________________________________________________________________________ set names for all the sequences
# 0 mismatches
set0 <- DNAStringSet(list(MG1655,
BW2952,
REL606,
APEC078,
IAI1,
E11368,
S88,
UTI89,
E2348,
D042)
)

# need the template strand sequence to align properly later
set0 <- reverseComplement(set0)
set0
names(set0) <- c('MG1655',
                 'BW2952',
                 'REL606',
                 'APEC078',
                 'IAI1',
                 'E11368',
                 'S88',
                 'UTI89',
                 'E2348',
                 '042')

set0 # sanity


# need to split up  set5 into plus and min matches to reverse complement the minus strand matches
set5_min_temp <- DNAStringSet(list(IAI39,
                              SMS35,
                              CE10))

set5_pos_temp <- DNAStringSet(list(TW14359,
                                   Sakai,
                                   EDL933))

set5_min <- reverseComplement(set5_min_temp)

set5 <- DNAStringSet(c(set5_min, set5_pos_temp))

set5

names(set5) <- c('IAI39',
                 'SMS-3-5',
                 'CE10',
                 'TW14359',
                 'Sakai',
                 'EDL993')

set5 # sanity


# combine sets into one that can be aligned with MSA
all_terY <- DNAStringSet(c(set0, set5))

all_terY


#________________________________________________________________________________ Align with MSA

terY_alignment <- msa(all_terY, order = 'input', method = 'Muscle')

# sanity
terY_alignment

# set working directory to where pseudo_ter is located in BT2
setwd(dirname(file.choose()))

# send to PDF
msaPrettyPrint(terY_alignment,file = 'terY_alignment_Ecoli_V3.pdf', output = c('pdf'),
               shadingMode = 'identical', verbose = TRUE,
               shadingColors = 'grays', shadingModeArg=NA,
               consensusThreshold=100,
               paperWidth=5, paperHeight=4, margins=c(0.3, 0.3))







