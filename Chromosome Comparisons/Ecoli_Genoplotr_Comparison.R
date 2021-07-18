# Using genoplotr from Bowtie2 data not BLAST to compare the differences

library(Biostrings)
library(DECIPHER)
library(msa)
library(seqinr)
library(genoPlotR)

# read in the CSV from find_ter_BT2_Samtools.R
df <- read.csv('C:\\Users\\Danie\\Documents\\R\\termination\\bowtie2\\Ecoli_geneome_ter_bam.csv', header = TRUE)

#sanity
df

# read in the CSV files from SAM_ter_tocsv.R program:
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




# We need 4 columns for genoplotr to work
# name, start, end, strand
# drop the rname column and seq

#A
MG1655 <- subset(MG1655, select = -c(rname, seq))
BW2952 <- subset(BW2952, select = -c(rname, seq))
REL606 <- subset(REL606, select = -c(rname, seq))

#B1
APEC078 <- subset(APEC078, select = -c(rname, seq))
IAI1    <- subset(IAI1, select = -c(rname, seq))
E11368  <- subset(E11368, select = -c(rname, seq))

#B2
S88   <- subset(S88, select = -c(rname, seq))
UTI89 <- subset(UTI89, select = -c(rname, seq))
E2348 <- subset(E2348, select = -c(rname, seq))

#D
IAI39   <- subset(IAI39, select = -c(rname, seq))
SMS35   <- subset(SMS35, select = -c(rname, seq))
UMN026  <- subset(UMN026, select = -c(rname, seq))
CE10    <- subset(CE10, select = -c(rname, seq))
D042    <- subset(D042, select = -c(rname, seq))

#E
TW14359 <- subset(TW14359, select = -c(rname, seq))
Sakai   <- subset(Sakai, select = -c(rname, seq))
EDL933  <- subset(EDL933, select = -c(rname, seq))


# turn qwidth into end by adding qwidth to pos
#A
MG1655$qwidth <- MG1655$pos + MG1655$qwidth
BW2952$qwidth <- BW2952$pos + BW2952$qwidth
REL606$qwidth <- REL606$pos + REL606$qwidth

#B1
APEC078$qwidth <- APEC078$pos + APEC078$qwidth
IAI1$qwidth <- IAI1$pos + IAI1$qwidth
E11368$qwidth <- E11368$pos + E11368$qwidth

#B2
S88$qwidth <- S88$pos + S88$qwidth
UTI89$qwidth <- UTI89$pos + UTI89$qwidth
E2348$qwidth <- E2348$pos + E2348$qwidth

#D
IAI39$qwidth <- IAI39$pos + IAI39$qwidth
SMS35$qwidth <- SMS35$pos + SMS35$qwidth
UMN026$qwidth <- UMN026$pos + UMN026$qwidth
CE10$qwidth <- CE10$pos + CE10$qwidth
D042$qwidth <- D042$pos + D042$qwidth

#E
TW14359$qwidth <- TW14359$pos + TW14359$qwidth
Sakai$qwidth <- Sakai$pos + Sakai$qwidth
EDL933$qwidth <- EDL933$pos + EDL933$qwidth




# reset the column names
#A
colnames(MG1655) <- c('name', 'start', 'end', 'strand')
colnames(BW2952) <- c('name', 'start', 'end', 'strand')
colnames(REL606) <- c('name', 'start', 'end', 'strand')

#B1
colnames(APEC078) <- c('name', 'start', 'end', 'strand')
colnames(IAI1) <- c('name', 'start', 'end', 'strand')
colnames(E11368) <- c('name', 'start', 'end', 'strand')

#B2
colnames(S88) <- c('name', 'start', 'end', 'strand')
colnames(UTI89) <- c('name', 'start', 'end', 'strand')
colnames(E2348) <- c('name', 'start', 'end', 'strand')

#D
colnames(IAI39) <- c('name', 'start', 'end', 'strand')
colnames(SMS35) <- c('name', 'start', 'end', 'strand')
colnames(UMN026) <- c('name', 'start', 'end', 'strand')
colnames(CE10) <- c('name', 'start', 'end', 'strand')
colnames(D042) <- c('name', 'start', 'end', 'strand')

#E
colnames(TW14359) <- c('name', 'start', 'end', 'strand')
colnames(Sakai) <- c('name', 'start', 'end', 'strand')
colnames(EDL933) <- c('name', 'start', 'end', 'strand')


# SANITY
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

###########################################################
#     now use Genoplotr package andf dna_seg objects
###########################################################

# Turn into dna_seg
dna1 <- dna_seg(MG1655)
dna2 <- dna_seg(BW2952)
dna3 <- dna_seg(REL606)

dna4 <- dna_seg(APEC078)
dna5 <- dna_seg(IAI1)
dna6 <- dna_seg(E11368)

dna7 <- dna_seg(S88)
dna8 <- dna_seg(UTI89)
dna9 <- dna_seg(E2348)

dna10 <- dna_seg(IAI39)
dna11 <- dna_seg(SMS35)
dna12 <- dna_seg(UMN026)
dna13 <- dna_seg(CE10)
dna14 <- dna_seg(D042)

dna15 <- dna_seg(TW14359)
dna16 <- dna_seg(Sakai)
dna17 <- dna_seg(EDL933)

cols <- c('red', 'blue', 'orange', 'green', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey')

# SPECIFY THE COL COLOUR
dna1$col <- cols
dna2$col <- cols
dna3$col <- cols

dna4$col <- cols
dna5$col <- cols
dna6$col <- cols

dna7$col <- cols
dna8$col <- cols
dna9$col <- cols

dna10$col <- cols
dna11$col <- cols
dna12$col <- cols
dna13$col <- cols
dna14$col <- cols

dna15$col <- cols
dna16$col <- cols
dna17$col <- cols

# turn into a massive list
dna_segs <-list(dna1,
                dna2,
                dna3,
                dna4,
                dna5,
                dna6,
                dna7,
                dna8,
                dna9,
                dna10,
                dna11,
                dna12,
                dna13,
                dna14,
                dna15,
                dna16,
                dna17)

# sanity
dna1
dna2
dna3

dna4
dna5
dna6

dna7
dna8
dna9

dna10
dna11
dna12
dna13
dna14

dna15
dna16
dna17



# give names which will be used as the plot label
names(dna_segs) <- c('MG1655','BW2952','REL606','APEC078','IAI1','11368',
                     'S88','UTI89','E2348/69','IAI39','SMS35','UMN026','CE10',
                     '042','TW14359','Sakai','EDL933')

# SANITY
dna_segs

## Calculating middle positions
mid_pos1 <- middle(dna_segs[[1]])
mid_pos2 <- middle(dna_segs[[2]])
mid_pos3 <- middle(dna_segs[[3]])
mid_pos4 <- middle(dna_segs[[4]])
mid_pos5 <- middle(dna_segs[[5]])
mid_pos6 <- middle(dna_segs[[6]])
mid_pos7 <- middle(dna_segs[[7]])
mid_pos8 <- middle(dna_segs[[8]])
mid_pos9 <- middle(dna_segs[[9]])
mid_pos10 <- middle(dna_segs[[10]])
mid_pos11 <- middle(dna_segs[[11]])
mid_pos12 <- middle(dna_segs[[12]])
mid_pos13 <- middle(dna_segs[[13]])
mid_pos14 <- middle(dna_segs[[14]])
mid_pos15 <- middle(dna_segs[[15]])
mid_pos16 <- middle(dna_segs[[16]])
mid_pos17 <- middle(dna_segs[[17]])

# Annotations
annot1 <- annotation(x1 = mid_pos1, text = dna_segs[[1]]$name, rot = 90, col = 'black')
annot2 <- annotation(x1 = mid_pos2, text = dna_segs[[2]]$name, rot = 90, col = 'black')
annot3 <- annotation(x1 = mid_pos3, text = dna_segs[[3]]$name, rot = 90, col = 'black')
annot4 <- annotation(x1 = mid_pos4, text = dna_segs[[4]]$name, rot = 90, col = 'black')
annot5 <- annotation(x1 = mid_pos5, text = dna_segs[[5]]$name, rot = 90, col = 'black')
annot6 <- annotation(x1 = mid_pos6, text = dna_segs[[6]]$name, rot = 90, col = 'black')
annot7 <- annotation(x1 = mid_pos7, text = dna_segs[[7]]$name, rot = 90, col = 'black')
annot8 <- annotation(x1 = mid_pos8, text = dna_segs[[8]]$name, rot = 90, col = 'black')
annot9 <- annotation(x1 = mid_pos9, text = dna_segs[[9]]$name, rot = 90, col = 'black')
annot10 <- annotation(x1 = mid_pos10, text = dna_segs[[10]]$name, rot = 90, col = 'black')
annot11 <- annotation(x1 = mid_pos11, text = dna_segs[[11]]$name, rot = 90, col = 'black')
annot12 <- annotation(x1 = mid_pos12, text = dna_segs[[12]]$name, rot = 90, col = 'black')
annot13 <- annotation(x1 = mid_pos13, text = dna_segs[[13]]$name, rot = 90, col = 'black')
annot14 <- annotation(x1 = mid_pos14, text = dna_segs[[14]]$name, rot = 90, col = 'black')
annot15 <- annotation(x1 = mid_pos15, text = dna_segs[[15]]$name, rot = 90, col = 'black')
annot16 <- annotation(x1 = mid_pos16, text = dna_segs[[16]]$name, rot = 90, col = 'black')
annot17 <- annotation(x1 = mid_pos17, text = dna_segs[[17]]$name, rot = 90, col = 'black')

#list of annotation objects
annots <- list(annot1, annot2, annot3, annot4, annot5, annot6, annot7, annot8, annot9, annot10, annot11, annot12, annot13, annot14, annot15, annot16, annot17)

# plot
plot_gene_map(dna_segs=dna_segs, comparisons=NULL,
              annotations = annots, annotation_height = 3, annotation_cex = 0.6,
              main = 'E.coli Ter Locations Determined by BOWTIE2',
              dna_seg_scale=FALSE, gene_type = 'side_blocks',  dna_seg_label_cex=0.8)


###################################################################################################################
#
#       FIGURE IS THE SAME AS FROM BLAST PROGRAM,  NO ISSUE HERE CAN USE THE BLAST PROGRAM IN M&M
#
###################################################################################################################




