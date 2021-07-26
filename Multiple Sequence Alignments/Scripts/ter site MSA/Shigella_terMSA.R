
library(Biostrings)
library(DECIPHER)
library(msa)
library(seqinr)
library(msa)

# read in the csv from program 2 
flex_2a <- read.csv('flex_2a.csv')
flex_Shi06 <- read.csv('flex_Shi06.csv')
flex_8401 <- read.csv('flex_8401.csv')
flex_2002017 <- read.csv('flex_2002017.csv')
flex_2003036 <- read.csv('flex_2003036.csv')
boy_3038 <- read.csv('boy_3038.csv')
boy_600080 <- read.csv('boy_600080.csv')
boy_600690 <- read.csv('boy_600690.csv')
boy_227 <- read.csv('boy_227.csv')
sd197 <- read.csv('sd197.csv')
ss046 <-read.csv('ss046.csv')

########
#UPDATE 15/01/2021 (ish)
#remove terD if it is not a true terD 
########
flex_2a      <- subset(flex_2a[-4,])
flex_Shi06   <- subset(flex_Shi06[-4,])
flex_8401    <- subset(flex_8401[-4,])
flex_2002017 <- subset(flex_2002017[-4,])
flex_2003036 <- subset(flex_2003036[-4,])
boy_600080   <- subset(boy_600080[-4,])



#SANITY
flex_2a
flex_Shi06
flex_8401
flex_2002017
flex_2003036
boy_3038
boy_600080
boy_600690
boy_227
sd197
ss046

############################################################################################################
# add in terD found with vmatchPattern program for sb600690
terD_600690 <- DNAStringSet(Views(DNAString(sb69_str), start = unlist(startD), end = unlist(endD))[3])
terD <- data.frame(terD_600690)
boy_600690[4,] <- c('boy_600690', 'terD', '2866624', 23, '+', terD)
rownames(boy_600690) <- NULL
boy_600690
############################################################################################################


###################################
# CONCATENATE INTO ONE DF
###################################
# the do.call functions allows existing functions to be used but with more variables

shigella <- do.call('rbind',list(flex_2a,
                              flex_Shi06,
                              flex_8401,
                              flex_2002017,
                              flex_2003036,
                              boy_3038,
                              boy_600080,
                              boy_600690,
                              boy_227,
                              sd197,
                              ss046))

#sanity
shigella

####################################
# Put into 10 objects terA-J
####################################
terA <- shigella[shigella$qname=='terA',]
terB <- shigella[shigella$qname=='terB',]
terC <- shigella[shigella$qname=='terC',]
terD <- shigella[shigella$qname=='terD',]

#reset index
rownames(terA) <- NULL
rownames(terB) <- NULL
rownames(terC) <- NULL
rownames(terD) <- NULL

#sanity
terA
terB
terC
terD


#################################
#
#           ALIGNMENT
#
#################################



############################################### 
#                   terA

# ALIGN then name
terA_al <- AlignSeqs(DNAStringSet(terA$seq))
names(terA_al) <- terA$rname

# Browse / View
BrowseSeqs(terA_al, highlight = 0)

# MSA
msa_terA <- msa(terA_al,  order = 'input', method = 'Muscle')
msa_terA

# MSAprettyprint
# from  http://www.bioinf.jku.at/software/msa/Example.R
msaPrettyPrint(x=msa_terA, file="Shigella_terA_Alignment.pdf",
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
terB_al <- AlignSeqs(DNAStringSet(terB$seq))
names(terB_al) <- terB$rname
terB_al
# Browse / View
BrowseSeqs(terB_al, highlight = 0)

# MSA
msa_terB <- msa(terB_al,  order = 'input', method = 'Muscle')
msa_terB

# MSAprettyprint
# from  http://www.bioinf.jku.at/software/msa/Example.R
msaPrettyPrint(x=msa_terB, file="Shigella_terB_Alignment.pdf",
               shadingMode="identical",
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
terC_al <- AlignSeqs(DNAStringSet(terC$seq))
names(terC_al) <- terC$rname
terC_al
# Browse / View
BrowseSeqs(terC_al, highlight = 0)

# MSA
msa_terC <- msa(terC_al,  order = 'input', method = 'Muscle')
msa_terC

# MSAprettyprint
# from  http://www.bioinf.jku.at/software/msa/Example.R
msaPrettyPrint(x=msa_terC, file="Shigella_terC_Alignment.pdf",
               shadingMode="identical",
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
terD_al <- AlignSeqs(DNAStringSet(terD$seq))
names(terD_al) <- terD$rname
terD_al

# Browse / View
BrowseSeqs(terD_al, highlight = 0)

# MSA
msa_terD <- msa(terD_al,  order = 'input', method = 'Muscle')
msa_terD

# MSAprettyprint
# from  http://www.bioinf.jku.at/software/msa/Example.R
msaPrettyPrint(x=msa_terD, file="Shigella_terD_Alignment.pdf",
               shadingMode="identical",
               shadingColors="grays",
               logoColors="chemical",
               showLogoScale="right",
               showConsensus=c("bottom"),
               consensusColors=c("ColdHot"),
               showLegend=FALSE,
               askForOverwrite=FALSE, verbose = TRUE)


# CHEKC TERD AS THERE ARE MORE THAN HALF THE SITES THAT DON'T HAVE G6