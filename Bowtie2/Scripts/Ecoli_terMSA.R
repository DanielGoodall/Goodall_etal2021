# Multiple sequence alignment from BLAST csv - updated  method

library(Biostrings)
library(DECIPHER)
library(msa)
library(seqinr)

# read in the CSV files from ter_clean program:
#A
MG1655 <- read.csv('C:\\Users\\Danie\\Documents\\R\\termination\\bowtie2\\Ecoli_matched_ter_sequences_csv\\MG1655.csv')
BW2952 <- read.csv('C:\\Users\\Danie\\Documents\\R\\termination\\bowtie2\\Ecoli_matched_ter_sequences_csv\\BW2952.csv')
REL606 <- read.csv('C:\\Users\\Danie\\Documents\\R\\termination\\bowtie2\\Ecoli_matched_ter_sequences_csv\\REL606.csv')

#B1
APEC078 <- read.csv('C:\\Users\\Danie\\Documents\\R\\termination\\bowtie2\\Ecoli_matched_ter_sequences_csv\\APEC078.csv')
IAI1    <- read.csv('C:\\Users\\Danie\\Documents\\R\\termination\\bowtie2\\Ecoli_matched_ter_sequences_csv\\IAI1.csv')
E11368  <- read.csv('C:\\Users\\Danie\\Documents\\R\\termination\\bowtie2\\Ecoli_matched_ter_sequences_csv\\E11368.csv')

#B2
S88   <- read.csv('C:\\Users\\Danie\\Documents\\R\\termination\\bowtie2\\Ecoli_matched_ter_sequences_csv\\S88.csv')
UTI89 <- read.csv('C:\\Users\\Danie\\Documents\\R\\termination\\bowtie2\\Ecoli_matched_ter_sequences_csv\\UTI89.csv')
E2348 <- read.csv('C:\\Users\\Danie\\Documents\\R\\termination\\bowtie2\\Ecoli_matched_ter_sequences_csv\\E2348.csv')

#D
IAI39  <- read.csv('C:\\Users\\Danie\\Documents\\R\\termination\\bowtie2\\Ecoli_matched_ter_sequences_csv\\IAI39.csv')
SMS35  <- read.csv('C:\\Users\\Danie\\Documents\\R\\termination\\bowtie2\\Ecoli_matched_ter_sequences_csv\\SMS35.csv')
UMN026 <- read.csv('C:\\Users\\Danie\\Documents\\R\\termination\\bowtie2\\Ecoli_matched_ter_sequences_csv\\UMN026.csv')
CE10   <- read.csv('C:\\Users\\Danie\\Documents\\R\\termination\\bowtie2\\Ecoli_matched_ter_sequences_csv\\CE10.csv')
D042   <- read.csv('C:\\Users\\Danie\\Documents\\R\\termination\\bowtie2\\Ecoli_matched_ter_sequences_csv\\D042.csv')

#E
TW14359 <- read.csv('C:\\Users\\Danie\\Documents\\R\\termination\\bowtie2\\Ecoli_matched_ter_sequences_csv\\TW14359.csv')
Sakai   <- read.csv('C:\\Users\\Danie\\Documents\\R\\termination\\bowtie2\\Ecoli_matched_ter_sequences_csv\\Sakai.csv')
EDL933  <- read.csv('C:\\Users\\Danie\\Documents\\R\\termination\\bowtie2\\Ecoli_matched_ter_sequences_csv\\EDL933.csv')


# Saniity
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

###################################
# CONCATENATE INTO ONE DF
###################################
# the do.call functions allows existing functions to be used but with more variables

ecoli <- do.call('rbind',list(MG1655,
BW2952,
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
EDL933))

#sanity
ecoli

####################################
# Put into 10 objects terA-J
####################################
terA <- ecoli[ecoli$qname=='terA',]
terB <- ecoli[ecoli$qname=='terB',]
terC <- ecoli[ecoli$qname=='terC',]
terD <- ecoli[ecoli$qname=='terD',]
terE <- ecoli[ecoli$qname=='terE',]
terF <- ecoli[ecoli$qname=='terF',]
terG <- ecoli[ecoli$qname=='terG',]
terH <- ecoli[ecoli$qname=='terH',]
terI <- ecoli[ecoli$qname=='terI',]
terJ <- ecoli[ecoli$qname=='terJ',]

# RESET INDEX
rownames(terA) <- NULL
rownames(terB) <- NULL
rownames(terC) <- NULL
rownames(terD) <- NULL
rownames(terE) <- NULL
rownames(terF) <- NULL
rownames(terG) <- NULL
rownames(terH) <- NULL
rownames(terI) <- NULL
rownames(terJ) <- NULL

# SANITY
terA
terB
terC
terD
terE
terF
terG
terH
terI
terJ

# NAMES OF THE GENOMES
gen_names <- toString(terA$rname)
gen_names



#################################
#           ALIGNMENT
################################

# DNASTRINGSET
terA_dna <- DNAStringSet(terA$seq)
names(terA_dna) <- terA$rname

terB_dna <- DNAStringSet(terB$seq)
names(terB_dna) <- terB$rname

terC_dna <- DNAStringSet(terC$seq)
names(terC_dna) <- terC$rname

terD_dna <- DNAStringSet(terD$seq)
names(terD_dna) <- terD$rname

terE_dna <- DNAStringSet(terE$seq)
names(terE_dna) <- terE$rname

terF_dna <- DNAStringSet(terF$seq)
names(terF_dna) <- terF$rname

terG_dna <- DNAStringSet(terG$seq)
names(terG_dna) <- terG$rname

terH_dna <- DNAStringSet(terH$seq)
names(terH_dna) <- terH$rname

terI_dna <- DNAStringSet(terI$seq)
names(terI_dna) <- terI$rname

terJ_dna <- DNAStringSet(terJ$seq)
names(terJ_dna) <- terJ$rname


# SANITY
terA_dna
terB_dna
terC_dna
terD_dna
terE_dna
terF_dna
terG_dna
terH_dna
terI_dna
terJ_dna

############################################### 
#                   terA

# ALIGN
terA_al <- AlignSeqs(terA_dna)

# Browse / View
BrowseSeqs(terA_al, highlight = 0)

# MSA
msa_terA <- msa(terA_al,  order = 'input', method = 'Muscle')
msa_terA

# MSAprettyprint
# from  http://www.bioinf.jku.at/software/msa/Example.R
msaPrettyPrint(x=msa_terA, file="Ecoli_terA_Alignment.pdf",
               shadingMode="identical",
               shadingColors="grays",
               logoColors="chemical",
               showLogoScale="right",
               showConsensus=c("bottom"),
               consensusColors=c("ColdHot"),
               showLegend=FALSE,
               askForOverwrite=FALSE, verbose = TRUE)


############################################### 
#                   terB

# ALIGN
terB_al <- AlignSeqs(terB_dna)

# Browse / View
BrowseSeqs(terB_al, highlight = 0)

# MSA
msa_terB <- msa(terB_al,  order = 'input')
msa_terB

# MSAprettyprint
# from  http://www.bioinf.jku.at/software/msa/Example.R
msaPrettyPrint(x=msa_terB, file="Ecoli_terB_Alignment.pdf",
               shadingMode="functional",
               shadingModeArg="structure",
               shadingColors="grays",
               logoColors="chemical",
               showLogoScale="right",
               showConsensus=c("bottom"),
               consensusColors=c("ColdHot"),
               showLegend=FALSE,
               askForOverwrite=FALSE, verbose = TRUE)

############################################### 
#                   terC

# ALIGN
terC_al <- AlignSeqs(terC_dna)

# Browse / View
BrowseSeqs(terC_al, highlight = 0)

# MSA
msa_terC <- msa(terC_al,  order = 'input')
msa_terC

# MSAprettyprint
# from  http://www.bioinf.jku.at/software/msa/Example.R
msaPrettyPrint(x=msa_terC, file="Ecoli_terC_Alignment.pdf",
               shadingMode="functional",
               shadingModeArg="structure",
               shadingColors="grays",
               logoColors="chemical",
               showLogoScale="right",
               showConsensus=c("bottom"),
               consensusColors=c("ColdHot"),
               showLegend=FALSE,
               askForOverwrite=FALSE, verbose = TRUE)


############################################### 
#                   terD

# ALIGN
terD_al <- AlignSeqs(terD_dna)

# Browse / View
BrowseSeqs(terD_al, highlight = 0)

# MSA
msa_terD <- msa(terD_al,  order = 'input')
msa_terD

# MSAprettyprint
# from  http://www.bioinf.jku.at/software/msa/Example.R
msaPrettyPrint(x=msa_terD, file="Ecoli_terD_Alignment.pdf",
               shadingMode="functional",
               shadingModeArg="structure",
               shadingColors="grays",
               logoColors="chemical",
               showLogoScale="right",
               showConsensus=c("bottom"),
               consensusColors=c("ColdHot"),
               showLegend=FALSE,
               askForOverwrite=FALSE, verbose = TRUE)


############################################### 
#                   terE

# ALIGN
terE_al <- AlignSeqs(terE_dna)

# Browse / View
BrowseSeqs(terE_al, highlight = 0)

# MSA
msa_terE <- msa(terE_al,  order = 'input')
msa_terE

# MSAprettyprint
# from  http://www.bioinf.jku.at/software/msa/Example.R
msaPrettyPrint(x=msa_terE, file="Ecoli_terE_Alignment.pdf",
               shadingMode="functional",
               shadingModeArg="structure",
               shadingColors="grays",
               logoColors="chemical",
               showLogoScale="right",
               showConsensus=c("bottom"),
               consensusColors=c("ColdHot"),
               showLegend=FALSE,
               askForOverwrite=FALSE, verbose = TRUE)


############################################### 
#                   terF

# ALIGN
terF_al <- AlignSeqs(terF_dna)

# Browse / View
BrowseSeqs(terF_al, highlight = 0)

# MSA
msa_terF <- msa(terF_al,  order = 'input')
msa_terF

# MSAprettyprint
# from  http://www.bioinf.jku.at/software/msa/Example.R
msaPrettyPrint(x=msa_terF, file="Ecoli_terF_Alignment.pdf",
               shadingMode="functional",
               shadingModeArg="structure",
               shadingColors="grays",
               logoColors="chemical",
               showLogoScale="right",
               showConsensus=c("bottom"),
               consensusColors=c("ColdHot"),
               showLegend=FALSE,
               askForOverwrite=FALSE, verbose = TRUE)


############################################### 
#                   terG

# ALIGN
terG_al <- AlignSeqs(terG_dna)

# Browse / View
BrowseSeqs(terG_al, highlight = 0)

# MSA
msa_terG <- msa(terG_al,  order = 'input')
msa_terG

# MSAprettyprint
# from  http://www.bioinf.jku.at/software/msa/Example.R
msaPrettyPrint(x=msa_terG, file="Ecoli_terG_Alignment.pdf",
               shadingMode="functional",
               shadingModeArg="structure",
               shadingColors="grays",
               logoColors="chemical",
               showLogoScale="right",
               showConsensus=c("bottom"),
               consensusColors=c("ColdHot"),
               showLegend=FALSE,
               askForOverwrite=FALSE, verbose = TRUE)


############################################### 
#                   terH

# ALIGN
terH_al <- AlignSeqs(terH_dna)

# Browse / View
BrowseSeqs(terH_al, highlight = 0)

# MSA
msa_terH <- msa(terH_al,  order = 'input')
msa_terH

# MSAprettyprint
# from  http://www.bioinf.jku.at/software/msa/Example.R
msaPrettyPrint(x=msa_terH, file="Ecoli_terH_Alignment.pdf",
               shadingMode="functional",
               shadingModeArg="structure",
               shadingColors="grays",
               logoColors="chemical",
               showLogoScale="right",
               showConsensus=c("bottom"),
               consensusColors=c("ColdHot"),
               showLegend=FALSE,
               askForOverwrite=FALSE, verbose = TRUE)

############################################### 
#                   terI

# ALIGN
terI_al <- AlignSeqs(terI_dna)

# Browse / View
BrowseSeqs(terI_al, highlight = 0)

# MSA
msa_terI <- msa(terI_al,  order = 'input')
msa_terI

# MSAprettyprint
# from  http://www.bioinf.jku.at/software/msa/Example.R
msaPrettyPrint(x=msa_terI, file="Ecoli_terI_Alignment.pdf",
               shadingMode="functional",
               shadingModeArg="structure",
               shadingColors="grays",
               logoColors="chemical",
               showLogoScale="right",
               showConsensus=c("bottom"),
               consensusColors=c("ColdHot"),
               showLegend=FALSE,
               askForOverwrite=FALSE, verbose = TRUE)

############################################### 
#                   terJ

# ALIGN
terJ_al <- AlignSeqs(terJ_dna)

# Browse / View
BrowseSeqs(terJ_al, highlight = 0)

# MSA
msa_terJ <- msa(terJ_al,  order = 'input')
msa_terJ

# MSAprettyprint
# from  http://www.bioinf.jku.at/software/msa/Example.R
msaPrettyPrint(x=msa_terJ, file="Ecoli_terJ_Alignment.pdf",
               shadingMode="functional",
               shadingModeArg="structure",
               shadingColors="grays",
               logoColors="chemical",
               showLogoScale="right",
               showConsensus=c("bottom"),
               consensusColors=c("ColdHot"),
               showLegend=FALSE,
               askForOverwrite=FALSE, verbose = TRUE)





#################################################################################################################################
#
#                                                             5 PHYLOGROUP MSA
#                                                               NEXT PROGRAM
#
#################################################################################################################################

















