# Multiple sequence alignment from BLAST csv - updated  method

library(Biostrings)
library(DECIPHER)
library(msa)
library(seqinr)
library(dplyr)

# E.coli genomes
genomes <- readDNAStringSet(file.choose())
genomes

# read in the CSV files from ter_clean program:
#A
MG1655 <- read.csv('MG1655.csv')
BW2952 <- read.csv('BW2952.csv')
REL606 <- read.csv('REL606.csv')

#B1
APEC078 <- read.csv('APEC078.csv')
IAI1    <- read.csv('IAI1.csv')
E11368  <- read.csv('E11368.csv')

#B2
S88   <- read.csv('S88.csv')
UTI89 <- read.csv('UTI89.csv')
E2348 <- read.csv('E2348.csv')

#D
IAI39  <- read.csv('IAI39.csv')
SMS35  <- read.csv('SMS35.csv')
UMN026 <- read.csv('UMN026.csv')
CE10   <- read.csv('CE10.csv')
D042   <- read.csv('D042.csv')

#E
TW14359 <- read.csv('TW14359.csv')
Sakai   <- read.csv('Sakai.csv')
EDL933  <- read.csv('EDL933.csv')


# Use subseq() to increase window of matched ter sites 25nt either side of ter site
#A
MG1655 <- MG1655 %>% mutate(pos = pos-25, qwidth = pos+48)
BW2952 <- BW2952 %>% mutate(pos = pos-25, qwidth = pos+48)
REL606 <- REL606 %>% mutate(pos = pos-25, qwidth = pos+48)

#B1
APEC078 <- APEC078 %>% mutate(pos = pos-25, qwidth = pos+48)
IAI1 <- IAI1 %>% mutate(pos = pos-25, qwidth = pos+48)
E11368 <- E11368 %>% mutate(pos = pos-25, qwidth = pos+48)

#B2
S88 <- S88 %>% mutate(pos = pos-25, qwidth = pos+48)
UTI89 <- UTI89 %>% mutate(pos = pos-25, qwidth = pos+48)
E2348 <- E2348 %>% mutate(pos = pos-25, qwidth = pos+48)

#D
IAI39 <- IAI39 %>% mutate(pos = pos-25, qwidth = pos+48)
SMS35 <- SMS35 %>% mutate(pos = pos-25, qwidth = pos+48)
UMN026 <- UMN026 %>% mutate(pos = pos-25, qwidth = pos+48)
CE10 <- CE10 %>% mutate(pos = pos-25, qwidth = pos+48)
D042 <- D042 %>% mutate(pos = pos-25, qwidth = pos+48)

#E
TW14359 <- TW14359 %>% mutate(pos = pos-25, qwidth = pos+48)
Sakai <- Sakai %>% mutate(pos = pos-25, qwidth = pos+48)
EDL933 <- EDL933 %>% mutate(pos = pos-25, qwidth = pos+48)






# Update the strain dataframes with new seqs
#________________________________________________________ MG1655
MG1655_list <- list()
plusMG <- MG1655 %>% filter(strand == '+')
minMG <- MG1655 %>% filter(strand == '-')
plusMG
minMG
for (number in plusMG$pos) {
  seq_plus <- list(subseq(genomes$MG1655, start = number, width = 73))
  MG1655_list <- append(MG1655_list, seq_plus)
}
for (number in minMG$pos) {
  seq_min <- list(reverseComplement(subseq(genomes$MG1655, start = number, width = 73)))
  MG1655_list <- append(MG1655_list,seq_min)
}
MG1655_dnaset <- DNAStringSet(MG1655_list)
names(MG1655_dnaset) <- c(plusMG$qname, minMG$qname)

tempMG <- data.frame(seq = MG1655_dnaset, name = c(plusMG$qname, minMG$qname))
tempMG<- tempMG %>% arrange(name)
tempMG
MG1655$seq <- tempMG$seq
MG1655

#________________________________________________________ BW2952
BW2952_list <- list()
plusBW <- BW2952 %>% filter(strand == '+')
minBW <- BW2952 %>% filter(strand == '-')
plusBW
minBW
for (number in plusBW$pos) {
  seq_plus <- list(subseq(genomes$BW2952, start = number, width = 73))
  BW2952_list <- append(BW2952_list, seq_plus)
}
for (number in minBW$pos) {
  seq_min <- list(reverseComplement(subseq(genomes$BW2952, start = number, width = 73)))
  BW2952_list <- append(BW2952_list, seq_min)
}
BW2952_dnaset <- DNAStringSet(BW2952_list)
names(BW2952_dnaset) <- c(plusBW$qname, minBW$qname)
BW2952_dnaset
tempBW <- data.frame(seq = BW2952_dnaset, name = c(plusBW$qname, minBW$qname))
tempBW<- tempBW %>% arrange(name)
tempBW
BW2952$seq <- tempBW$seq
BW2952

#________________________________________________________ REL606
REL606_list <- list()
plusREL <- REL606 %>% filter(strand == '+')
minREL <- REL606 %>% filter(strand == '-')
plusREL
minREL
for (number in plusREL$pos) {
  seq_plus <- list(subseq(genomes$REL606, start = number, width = 73))
  REL606_list <- append(REL606_list, seq_plus)
}
for (number in minREL$pos) {
  seq_min <- list(reverseComplement(subseq(genomes$REL606, start = number, width = 73)))
  REL606_list <- append(REL606_list, seq_min)
}
REL606_dnaset <- DNAStringSet(REL606_list)
names(REL606_dnaset) <- c(plusREL$qname, minREL$qname)

tempREL <- data.frame(seq = REL606_dnaset, name = c(plusREL$qname, minREL$qname))
tempREL<- tempREL %>% arrange(name)
tempREL
REL606$seq <- tempREL$seq
REL606

#________________________________________________________ APECO78
APEC078_list <- list()
plusAP <- APEC078 %>% filter(strand == '+')
minAP <- APEC078 %>% filter(strand == '-')
plusAP
minAP
for (number in plusAP$pos) {
  seq_plus <- list(subseq(genomes$APECO78, start = number, width = 73))
  APEC078_list <- append(APEC078_list, seq_plus)
}
for (number in minAP$pos) {
  seq_min <- list(reverseComplement(subseq(genomes$APECO78, start = number, width = 73)))
  APEC078_list <- append(APEC078_list, seq_min)
}
APEC078_dnaset <- DNAStringSet(APEC078_list)
names(MG1655_dnaset) <- c(plusAP$qname, minAP$qname)

tempAP <- data.frame(seq = MG1655_dnaset, name = c(plusAP$qname, minAP$qname))
tempAP<- tempAP %>% arrange(name)
tempAP
APEC078$seq <- tempAP$seq
APEC078

#________________________________________________________ IAI1
IAI1_list <- list()
plusIAI <- IAI1 %>% filter(strand == '+')
minIAI <- IAI1 %>% filter(strand == '-')
plusIAI
minIAI
for (number in plusIAI$pos) {
  seq_plus <- list(subseq(genomes$IAI1, start = number, width = 73))
  IAI1_list <- append(IAI1_list, seq_plus)
}
for (number in minIAI$pos) {
  seq_min <- list(reverseComplement(subseq(genomes$IAI1, start = number, width = 73)))
  IAI1_list <- append(IAI1_list, seq_min)
}
IAI1_dnaset <- DNAStringSet(IAI1_list)
names(MG1655_dnaset) <- c(plusIAI$qname, minIAI$qname)

tempIAI <- data.frame(seq = IAI1_dnaset, name = c(plusIAI$qname, minIAI$qname))
tempIAI<- tempIAI %>% arrange(name)
tempIAI
IAI1$seq <- tempIAI$seq
IAI1

#________________________________________________________ 11368
E11368_list <- list()
plusE11 <- E11368 %>% filter(strand == '+')
minE11 <- E11368 %>% filter(strand == '-')
plusE11
minE11
for (number in plusE11$pos) {
  seq_plus <- list(subseq(genomes$`11368`, start = number, width = 73))
  E11368_list <- append(E11368_list, seq_plus)
}
for (number in minE11$pos) {
  seq_min <- list(reverseComplement(subseq(genomes$`11368`, start = number, width = 73)))
  E11368_list <- append(E11368_list, seq_min)
}
E11368_dnaset <- DNAStringSet(E11368_list)
names(MG1655_dnaset) <- c(plusE11$qname, minE11$qname)

tempE11 <- data.frame(seq = E11368_dnaset, name = c(plusE11$qname, minE11$qname))
tempE11<- tempE11 %>% arrange(name)
tempE11
E11368$seq <- tempE11$seq
E11368

#________________________________________________________ S88
S88_list <- list()
plusS88 <- S88 %>% filter(strand == '+')
minS88 <- S88 %>% filter(strand == '-')
plusS88
minS88
for (number in plusS88$pos) {
  seq_plus <- list(subseq(genomes$S88, start = number, width = 73))
  S88_list <- append(S88_list, seq_plus)
}
for (number in minS88$pos) {
  seq_min <- list(reverseComplement(subseq(genomes$S88, start = number, width = 73)))
  S88_list <- append(S88_list, seq_min)
}
S88_dnaset <- DNAStringSet(S88_list)
names(S88_dnaset) <- c(plusS88$qname, minS88$qname)

tempS88 <- data.frame(seq = S88_dnaset, name = c(plusS88$qname, minS88$qname))
tempS88 <- tempS88 %>% arrange(name)
tempS88
S88$seq <- tempS88$seq
S88

#________________________________________________________ UTI89
UTI89_list <- list()
plusS8889 <- UTI89 %>% filter(strand == '+')
minUTI89 <- UTI89 %>% filter(strand == '-')
plusUTI89
minUTI89
for (number in plusUTI89$pos) {
  seq_plus <- list(subseq(genomes$UTI89, start = number, width = 73))
  UTI89_list <- append(UTI89_list, seq_plus)
}
for (number in minUTI89$pos) {
  seq_min <- list(reverseComplement(subseq(genomes$UTI89, start = number, width = 73)))
  UTI89_list <- append(UTI89_list, seq_min)
}
UTI89_dnaset <- DNAStringSet(UTI89_list)
names(UTI89_dnaset) <- c(plusUTI89$qname, minUTI89$qname)

tempUTI <- data.frame(seq = UTI89_dnaset, name = c(plusUTI89$qname, minUTI89$qname))
tempUTI <- tempUTI %>% arrange(name)
tempUTI
UTI89$seq <- tempUTI$seq
UTI89

#________________________________________________________ E2348
E2348_list <- list()
plusE2348 <- E2348 %>% filter(strand == '+')
minE2348 <- E2348 %>% filter(strand == '-')
plusE2348
minE2348
for (number in plusE2348$pos) {
  seq_plus <- list(subseq(genomes$E2348, start = number, width = 73))
  E2348_list <- append(E2348_list, seq_plus)
}
for (number in minE2348$pos) {
  seq_min <- list(reverseComplement(subseq(genomes$E2348, start = number, width = 73)))
  E2348_list <- append(E2348_list, seq_min)
}
E2348_dnaset <- DNAStringSet(E2348_list)
names(E2348_dnaset) <- c(plusE2348$qname, minE2348$qname)

tempE23 <- data.frame(seq = E2348_dnaset, name = c(plusE2348$qname, minE2348$qname))
tempE23<- tempE23 %>% arrange(name)
tempE23
E2348$seq <- tempE23$seq
E2348

#________________________________________________________ IAI39
IAI39_list <- list()
plusIAI39 <- IAI39 %>% filter(strand == '+')
minIAI39 <- IAI39 %>% filter(strand == '-')
plusIAI39
minIAI39
for (number in plusIAI39$pos) {
  seq_plus <- list(subseq(genomes$IAI39, start = number, width = 73))
  IAI39_list <- append(IAI39_list, seq_plus)
}
for (number in minIAI39$pos) {
  seq_min <- list(reverseComplement(subseq(genomes$IAI39, start = number, width = 73)))
  IAI39_list <- append(IAI39_list, seq_min)
}
IAI39_dnaset <- DNAStringSet(IAI39_list)
names(IAI39_dnaset) <- c(plusIAI39$qname, minIAI39$qname)

tempIAI39 <- data.frame(seq = IAI39_dnaset, name = c(plusIAI39$qname, minIAI39$qname))
tempIAI39<- tempIAI39 %>% arrange(name)
tempIAI39
IAI39$seq <- tempIAI39$seq
IAI39

#________________________________________________________ SMS-3-5
SMS35_list <- list()
plusSMS35 <- SMS35 %>% filter(strand == '+')
minSMS35 <- SMS35 %>% filter(strand == '-')
plusSMS35
minSMS35
for (number in plusSMS35$pos) {
  seq_plus <- list(subseq(genomes$`SMS-3-5`, start = number, width = 73))
  SMS35_list <- append(SMS35_list, seq_plus)
}
for (number in minSMS35$pos) {
  seq_min <- list(reverseComplement(subseq(genomes$`SMS-3-5`, start = number, width = 73)))
  SMS35_list <- append(SMS35_list, seq_min)
}
SMS35_dnaset <- DNAStringSet(SMS35_list)
names(SMS35_dnaset) <- c(plusSMS35$qname, minSMS35$qname)

tempSMS <- data.frame(seq = SMS35_dnaset, name = c(plusSMS35$qname, minSMS35$qname))
tempSMS<- tempSMS %>% arrange(name)
tempSMS
SMS35$seq <- tempSMS$seq
SMS35

#________________________________________________________ UMN026
UMN026_list <- list()
plusUMN026 <- UMN026 %>% filter(strand == '+')
minUMN026 <- UMN026 %>% filter(strand == '-')
plusUMN026
minUMN026
for (number in plusUMN026$pos) {
  seq_plus <- list(subseq(genomes$UMN026, start = number, width = 73))
  UMN026_list <- append(UMN026_list, seq_plus)
}
for (number in minUMN026$pos) {
  seq_min <- list(reverseComplement(subseq(genomes$UMN026, start = number, width = 73)))
  UMN026_list <- append(UMN026_list, seq_min)
}
UMN026_dnaset <- DNAStringSet(UMN026_list)
names(UMN026_dnaset) <- c(plusUMN026$qname, minUMN026$qname)

tempUMN <- data.frame(seq = UMN026_dnaset, name = c(plusUMN026$qname, minUMN026$qname))
tempUMN<- tempUMN %>% arrange(name)
tempUMN
UMN026$seq <- tempUMN$seq
UMN026

#________________________________________________________ CE10
CE10_list <- list()
plusCE10 <- CE10 %>% filter(strand == '+')
minCE10 <- CE10 %>% filter(strand == '-')
plusCE10
minCE10
for (number in plusCE10$pos) {
  seq_plus <- list(subseq(genomes$CE10, start = number, width = 73))
  CE10_list <- append(CE10_list, seq_plus)
}
for (number in minCE10$pos) {
  seq_min <- list(reverseComplement(subseq(genomes$CE10, start = number, width = 73)))
  CE10_list <- append(CE10_list, seq_min)
}
CE10_dnaset <- DNAStringSet(CE10_list)
names(CE10_dnaset) <- c(plusCE10$qname, minCE10$qname)

tempCE10 <- data.frame(seq = CE10_dnaset, name = c(plusCE10$qname, minCE10$qname))
tempCE10<- tempCE10 %>% arrange(name)
tempCE10
CE10$seq <- tempCE10$seq
CE10

#________________________________________________________ 042
D042_list <- list()
plus042 <- D042 %>% filter(strand == '+')
min042 <- D042 %>% filter(strand == '-')
plus042
min042
for (number in plus042$pos) {
  seq_plus <- list(subseq(genomes$`042`, start = number, width = 73))
  D042_list <- append(D042_list, seq_plus)
}
for (number in min042$pos) {
  seq_min <- list(reverseComplement(subseq(genomes$`042`, start = number, width = 73)))
  D042_list <- append(D042_list, seq_min)
}
D042_dnaset <- DNAStringSet(D042_list)
names(D042_dnaset) <- c(plus042$qname, min042$qname)

temp042 <- data.frame(seq = D042_dnaset, name = c(plus042$qname, min042$qname))
temp042<- temp042 %>% arrange(name)
temp042
D042$seq <- temp042$seq
D042

#________________________________________________________ TW14359
TW14359_list <- list()
plusTW <- TW14359 %>% filter(strand == '+')
minTW <- TW14359 %>% filter(strand == '-')
plusTW
minTW
for (number in plusTW$pos) {
  seq_plus <- list(subseq(genomes$TW14359, start = number, width = 73))
  TW14359_list <- append(TW14359_list, seq_plus)
}
for (number in minTW$pos) {
  seq_min <- list(reverseComplement(subseq(genomes$TW14359, start = number, width = 73)))
  TW14359_list <- append(TW14359_list, seq_min)
}
TW14359_dnaset <- DNAStringSet(TW14359_list)
names(TW14359_dnaset) <- c(plusTW$qname, minTW$qname)

tempTW <- data.frame(seq = TW14359_dnaset, name = c(plusTW$qname, minTW$qname))
tempTW<- tempTW %>% arrange(name)
tempTW
TW14359$seq <- tempTW$seq
TW14359

#________________________________________________________ Sakai
Sakai_list <- list()
plusSakai <- Sakai %>% filter(strand == '+')
minSakai <- Sakai %>% filter(strand == '-')
plusSakai
minSakai
for (number in plusSakai$pos) {
  seq_plus <- list(subseq(genomes$Sakai, start = number, width = 73))
  Sakai_list <- append(Sakai_list, seq_plus)
}
for (number in minSakai$pos) {
  seq_min <- list(reverseComplement(subseq(genomes$Sakai, start = number, width = 73)))
  Sakai_list <- append(Sakai_list, seq_min)
}
Sakai_dnaset <- DNAStringSet(Sakai_list)
names(Sakai_dnaset) <- c(plusSakai$qname, minSakai$qname)

tempSakai <- data.frame(seq = Sakai_dnaset, name = c(plusSakai$qname, minSakai$qname))
tempSakai<- tempSakai %>% arrange(name)
tempSakai
Sakai$seq <- tempSakai$seq
Sakai

#________________________________________________________ EDL933
EDL933_list <- list()
plusEDL <- EDL933 %>% filter(strand == '+')
minEDL <- EDL933 %>% filter(strand == '-')
plusEDL
minEDL
for (number in plusEDL$pos) {
  seq_plus <- list(subseq(genomes$EDL933, start = number, width = 73))
  EDL933_list <- append(EDL933_list, seq_plus)
}
for (number in minEDL$pos) {
  seq_min <- list(reverseComplement(subseq(genomes$EDL933, start = number, width = 73)))
  EDL933_list <- append(EDL933_list, seq_min)
}
EDL933_dnaset <- DNAStringSet(EDL933_list)
names(EDL933_dnaset) <- c(plusEDL$qname, minEDL$qname)

tempEDL <- data.frame(seq = EDL933_dnaset, name = c(plusEDL$qname, minEDL$qname))
tempEDL<- tempEDL %>% arrange(name)
tempEDL
EDL933$seq <- tempEDL$seq
EDL933


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
terA$seq
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
terA_al <- DNAMultipleAlignment(terA_dna)

# Browse / View
#BrowseSeqs(terF_al, highlight = 0)

# MSA
msa_terA <- msa(terA_al,  order = 'input', method = 'Muscle', gapOpening = NULL, gapExtension = NULL)
msa_terA

# MSAprettyprint
#setwd(dirname(file.choose()))
# y=c(start=25, end=48)
msaPrettyPrint(x=x, file="Ecoli_terA_Alignment.pdf",
               shadingMode="identical",
               shadingColors="grays",
               logoColors="rasmol",
               showLogoScale="right",
               showConsensus=c("bottom"),
               consensusColors=c("ColdHot"),
               showLegend=FALSE,
               askForOverwrite=FALSE, verbose = TRUE,
               paperWidth = 8, paperHeight = 5)


  ############################################### 
#                   terB

# ALIGN
terB_al <- DNAMultipleAlignment(terB_dna)



# MSA
msa_terB <- msa(terB_al,  order = 'input')
msa_terB

# MSAprettyprint
# from  http://www.bioinf.jku.at/software/msa/Example.R
msaPrettyPrint(x=terB_al, file="Ecoli_terB_Alignment.pdf",
               shadingMode="identical",
               shadingColors="grays",
               logoColors="rasmol",
               showLogoScale="right",
               showConsensus=c("bottom"),
               consensusColors=c("ColdHot"),
               showLegend=FALSE,
               askForOverwrite=FALSE, verbose = TRUE,
               paperWidth = 8, paperHeight = 5)

############################################### 
#                   terC

# ALIGN
terC_al <- DNAMultipleAlignment(terC_dna)



# MSA
msa_terC <- msa(terC_al,  order = 'input')
msa_terC

# MSAprettyprint
# from  http://www.bioinf.jku.at/software/msa/Example.R
msaPrettyPrint(x=terC_al, file="Ecoli_terC_Alignment.pdf",
               shadingMode="identical",
               shadingColors="grays",
               logoColors="rasmol",
               showLogoScale="right",
               showConsensus=c("bottom"),
               consensusColors=c("ColdHot"),
               showLegend=FALSE,
               askForOverwrite=FALSE, verbose = TRUE,
               paperWidth = 8, paperHeight = 5)


############################################### 
#                   terD

# ALIGN
terD_al <- DNAMultipleAlignment(terD_dna)



# MSA
msa_terD <- msa(terD_al,  order = 'input')
msa_terD

# MSAprettyprint
# from  http://www.bioinf.jku.at/software/msa/Example.R
msaPrettyPrint(x=terD_al, file="Ecoli_terD_Alignment.pdf",
               shadingMode="identical",
               shadingColors="grays",
               logoColors="rasmol",
               showLogoScale="right",
               showConsensus=c("bottom"),
               consensusColors=c("ColdHot"),
               showLegend=FALSE,
               askForOverwrite=FALSE, verbose = TRUE,
               paperWidth = 8, paperHeight = 5)


############################################### 
#                   terE

# ALIGN
terE_al <- DNAMultipleAlignment(terE_dna)



# MSA
msa_terE <- msa(terE_al,  order = 'input')
msa_terE

# MSAprettyprint
# from  http://www.bioinf.jku.at/software/msa/Example.R
msaPrettyPrint(x=terE_al, file="Ecoli_terE_Alignment.pdf",
               shadingMode="identical",
               shadingColors="grays",
               logoColors="rasmol",
               showLogoScale="right",
               showConsensus=c("bottom"),
               consensusColors=c("ColdHot"),
               showLegend=FALSE,
               askForOverwrite=FALSE, verbose = TRUE,
               paperWidth = 8, paperHeight = 5)


############################################### 
#                   terF

# ALIGN
terF_al <- DNAMultipleAlignment(terF_dna)



# MSA
msa_terF <- msa(terF_al,  order = 'input')
msa_terF

# MSAprettyprint
# from  http://www.bioinf.jku.at/software/msa/Example.R
msaPrettyPrint(x=terF_al, file="Ecoli_terF_Alignment.pdf",
               shadingMode="identical",
               shadingColors="grays",
               logoColors="rasmol",
               showLogoScale="right",
               showConsensus=c("bottom"),
               consensusColors=c("ColdHot"),
               showLegend=FALSE,
               askForOverwrite=FALSE, verbose = TRUE,
               paperWidth = 8, paperHeight = 5)


############################################### 
#                   terG

# ALIGN
terG_al <- DNAMultipleAlignment(terG_dna)



# MSA
msa_terG <- msa(terG_al,  order = 'input')
msa_terG

# MSAprettyprint
# from  http://www.bioinf.jku.at/software/msa/Example.R
msaPrettyPrint(x=terG_al, file="Ecoli_terG_Alignment.pdf",
               shadingMode="identical",
               shadingColors="grays",
               logoColors="rasmol",
               showLogoScale="right",
               showConsensus=c("bottom"),
               consensusColors=c("ColdHot"),
               showLegend=FALSE,
               askForOverwrite=FALSE, verbose = TRUE,
               paperWidth = 8, paperHeight = 5)


############################################### 
#                   terH

# ALIGN
terH_al <- DNAMultipleAlignment(terH_dna)



# MSA
msa_terH <- msa(terH_al,  order = 'input')
msa_terH

# MSAprettyprint
# from  http://www.bioinf.jku.at/software/msa/Example.R
msaPrettyPrint(x=terH_al, file="Ecoli_terH_Alignment.pdf",
               shadingMode="identical",
               shadingColors="grays",
               logoColors="rasmol",
               showLogoScale="right",
               showConsensus=c("bottom"),
               consensusColors=c("ColdHot"),
               showLegend=FALSE,
               askForOverwrite=FALSE, verbose = TRUE,
               paperWidth = 8, paperHeight = 5)

############################################### 
#                   terI

# ALIGN
terI_al <- DNAMultipleAlignment(terI_dna)



# MSA
msa_terI <- msa(terI_al,  order = 'input')
msa_terI

# MSAprettyprint
# from  http://www.bioinf.jku.at/software/msa/Example.R
msaPrettyPrint(x=terI_al, file="Ecoli_terI_Alignment.pdf",
               shadingMode="identical",
               shadingColors="grays",
               logoColors="rasmol",
               showLogoScale="right",
               showConsensus=c("bottom"),
               consensusColors=c("ColdHot"),
               showLegend=FALSE,
               askForOverwrite=FALSE, verbose = TRUE,
               paperWidth = 8, paperHeight = 5)

############################################### 
#                   terJ

# ALIGN
terJ_al <- DNAMultipleAlignment(terJ_dna)



# MSA
msa_terJ <- msa(terJ_al,  order = 'input')
msa_terJ

# MSAprettyprint
# from  http://www.bioinf.jku.at/software/msa/Example.R
msaPrettyPrint(x=terJ_al, file="Ecoli_terJ_Alignment.pdf",
               shadingMode="identical",
               shadingColors="grays",
               logoColors="rasmol",
               showLogoScale="right",
               showConsensus=c("bottom"),
               consensusColors=c("ColdHot"),
               showLegend=FALSE,
               askForOverwrite=FALSE, verbose = TRUE,
               paperWidth = 8, paperHeight = 5)





#################################################################################################################################
#
#                                                             Fin
#
#################################################################################################################################

















