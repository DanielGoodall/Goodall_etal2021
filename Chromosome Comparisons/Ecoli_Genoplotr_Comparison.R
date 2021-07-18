# updating the figure after paper rejection
## make the ter sites more visible
## delete the comparison lines
## colour code the ter labels with the arrow fill colour

# Using genoplotr from Bpwtie2 data not BLAST to compare the differences

library(Biostrings)
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


###############################
#         ADD DIF
###############################
dif <- read.csv('C:\\Users\\Danie\\Documents\\R\\termination\\BLAST_results\\Ecoli\\Ecoli_dif_blastn.csv', header = FALSE)
dif <- dif[-c(4,6,7)]
colnames(dif) <- c('name', 'start', 'end', 'strand')
dif$strand[dif$strand == 'minus'] <- '-'
dif$strand[dif$strand == 'plus'] <- '+'
dif


# add dif to the pre-made ter dfs
MG1655 <- rbind(MG1655, dif[dif$name=='MG1655',])
BW2952 <- rbind(BW2952, dif[dif$name=='BW2952',])
REL606 <- rbind(REL606, dif[dif$name=='REL606',])
APEC078 <- rbind(APEC078, dif[dif$name=='APECO78',])
IAI1 <- rbind(IAI1, dif[dif$name=='IAI1',])
E11368 <- rbind(E11368, dif[dif$name=='11368',])
S88 <- rbind(S88, dif[dif$name=='S88',])
UTI89 <- rbind(UTI89, dif[dif$name=='UTI89',])
E2348 <- rbind(E2348, dif[dif$name=='E2348/69',])
IAI39 <- rbind(IAI39, dif[dif$name=='IAI39',])
SMS35 <- rbind(SMS35, dif[dif$name=='SMS-3-5',])
UMN026 <- rbind(UMN026, dif[dif$name=='UMN026',])
CE10 <- rbind(CE10, dif[dif$name=='CE10',])
D042 <- rbind(D042, dif[dif$name=='42',])
TW14359 <- rbind(TW14359, dif[dif$name=='TW14359',])
Sakai <- rbind(Sakai, dif[dif$name=='Sakai',])
EDL933 <- rbind(EDL933, dif[dif$name=='EDL933',])


# RENAME [11] AS dif
MG1655$name[11] <- 'dif'
BW2952$name[11] <- 'dif'
REL606$name[11] <- 'dif'
APEC078$name[11] <- 'dif'
IAI1$name[11] <- 'dif'
E11368$name[11] <- 'dif'
S88$name[11] <- 'dif'
UTI89$name[11] <- 'dif'
E2348$name[11] <- 'dif'
IAI39$name[11] <- 'dif'
SMS35$name[11] <- 'dif'
UMN026$name[11] <- 'dif'
CE10$name[11] <- 'dif'
D042$name[11] <- 'dif'
TW14359$name[11] <- 'dif'
Sakai$name[11] <- 'dif'
EDL933$name[11] <- 'dif'

# reset indexes
#A
rownames(MG1655) <- NULL
rownames(BW2952) <- NULL
rownames(REL606) <- NULL

#B1
rownames(APEC078) <- NULL
rownames(IAI1) <- NULL
rownames(E11368) <- NULL

#B2
rownames(S88) <- NULL
rownames(UTI89) <- NULL
rownames(E2348) <- NULL

#D
rownames(IAI39) <- NULL
rownames(SMS35) <- NULL
rownames(UMN026) <- NULL
rownames(CE10) <- NULL
rownames(D042) <- NULL

#E
rownames(TW14359) <- NULL
rownames(Sakai) <- NULL
rownames(EDL933) <- NULL

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

cols <- c('red', 'blue', 'orange', 'darkgreen', 'purple', 'grey', 'grey', 'grey', 'grey', 'grey', 'black')

# SPECIFY THE COL COLOUR
dna1$col <- 'black'
dna2$col <- 'black'
dna3$col <- 'black'

dna4$col <- 'black'
dna5$col <- 'black'
dna6$col <- 'black'

dna7$col <- 'black'
dna8$col <- 'black'
dna9$col <- 'black'

dna10$col <- 'black'
dna11$col <- 'black'
dna12$col <- 'black'
dna13$col <- 'black'
dna14$col <- 'black'

dna15$col <- 'black'
dna16$col <- 'black'
dna17$col <- 'black'

# SPECIFY THE FILL COLOUR
dna1$fill <- cols
dna2$fill <- cols
dna3$fill <- cols

dna4$fill <- cols
dna5$fill <- cols
dna6$fill <- cols

dna7$fill <- cols
dna8$fill <- cols
dna9$fill <- cols

dna10$fill <- cols
dna11$fill <- cols
dna12$fill <- cols
dna13$fill <- cols
dna14$fill <- cols

dna15$fill <- cols
dna16$fill <- cols
dna17$fill <- cols

# INCREASING END POS BY 10000 NT TO ENLARGE ARROWS
dna1$end[1:10] <- c(dna1$start[1:10] + 50000)
dna2$end[1:10] <- c(dna2$start[1:10] + 50000)
dna3$end[1:10] <- c(dna3$start[1:10] + 50000)

dna4$end[1:10] <- c(dna4$start[1:10] + 50000)
dna5$end[1:10] <- c(dna5$start[1:10] + 50000)
dna6$end[1:10] <- c(dna6$start[1:10] + 50000)

dna7$end[1:10] <- c(dna7$start[1:10] + 50000)
dna8$end[1:10] <- c(dna8$start[1:10] + 50000)
dna9$end[1:10] <- c(dna9$start[1:10] + 50000)

dna10$end[1:10] <- c(dna10$start[1:10] + 50000)
dna11$end[1:10] <- c(dna11$start[1:10] + 50000)
dna12$end[1:10] <- c(dna12$start[1:10] + 50000)
dna13$end[1:10] <- c(dna13$start[1:10] + 50000)
dna14$end[1:10] <- c(dna14$start[1:10] + 50000)

dna15$end[1:10] <- c(dna15$start[1:10] + 50000)
dna16$end[1:10] <- c(dna16$start[1:10] + 50000)
dna17$end[1:10] <- c(dna17$start[1:10] + 50000)

# change dif to be side blocks
dna1$gene_type[11] <- 'side_blocks'
dna2$gene_type[11] <- 'side_blocks'
dna3$gene_type[11] <- 'side_blocks'
dna4$gene_type[11] <- 'side_blocks'
dna5$gene_type[11] <- 'side_blocks'
dna6$gene_type[11] <- 'side_blocks'
dna7$gene_type[11] <- 'side_blocks'
dna8$gene_type[11] <- 'side_blocks'
dna9$gene_type[11] <- 'side_blocks'
dna10$gene_type[11] <- 'side_blocks'
dna11$gene_type[11] <- 'side_blocks'
dna12$gene_type[11] <- 'side_blocks'
dna13$gene_type[11] <- 'side_blocks'
dna14$gene_type[11] <- 'side_blocks'
dna15$gene_type[11] <- 'side_blocks'
dna16$gene_type[11] <- 'side_blocks'
dna17$gene_type[11] <- 'side_blocks'

# specify dif parameters to stop terC overlap EDL933 and SMS35
dna17$start[11] <- c(dna17$start[3] + 50700)
dna17$end[11] <- c(dna17$start[3] + 50750)
dna17

dna11$start[11] <- c(dna11$start[3] + 50700)
dna11$end[11] <- c(dna11$start[3] + 50750)
dna11

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
names <- c('MG1655','BW2952','REL606','APEC078','IAI1','11368',
              'S88','UTI89','E2348/69','IAI39','SMS35','UMN026','CE10',
              '042','TW14359','Sakai','EDL933')
names(dna_segs) <- names

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
annot1 <- annotation(x1 = mid_pos1, text = dna_segs[[1]]$name, rot = 90, col = cols)
annot2 <- annotation(x1 = mid_pos2, text = dna_segs[[2]]$name, rot = 90, col = cols)
annot3 <- annotation(x1 = mid_pos3, text = dna_segs[[3]]$name, rot = 90, col = cols)
annot4 <- annotation(x1 = mid_pos4, text = dna_segs[[4]]$name, rot = 90, col = cols)
annot5 <- annotation(x1 = mid_pos5, text = dna_segs[[5]]$name, rot = 90, col = cols)
annot6 <- annotation(x1 = mid_pos6, text = dna_segs[[6]]$name, rot = 90, col = cols)
annot7 <- annotation(x1 = mid_pos7, text = dna_segs[[7]]$name, rot = 90, col = cols)
annot8 <- annotation(x1 = mid_pos8, text = dna_segs[[8]]$name, rot = 90, col = cols)
annot9 <- annotation(x1 = mid_pos9, text = dna_segs[[9]]$name, rot = 90, col = cols)
annot10 <- annotation(x1 = mid_pos10, text = dna_segs[[10]]$name, rot = 90, col = cols)
annot11 <- annotation(x1 = mid_pos11, text = dna_segs[[11]]$name, rot = 90, col = cols)
annot12 <- annotation(x1 = mid_pos12, text = dna_segs[[12]]$name, rot = 90, col = cols)
annot13 <- annotation(x1 = mid_pos13, text = dna_segs[[13]]$name, rot = 90, col = cols)
annot14 <- annotation(x1 = mid_pos14, text = dna_segs[[14]]$name, rot = 90, col = cols)
annot15 <- annotation(x1 = mid_pos15, text = dna_segs[[15]]$name, rot = 90, col = cols)
annot16 <- annotation(x1 = mid_pos16, text = dna_segs[[16]]$name, rot = 90, col = cols)
annot17 <- annotation(x1 = mid_pos17, text = dna_segs[[17]]$name, rot = 90, col = cols)

#list of annotation objects
annots <- list(annot1, annot2, annot3, annot4, annot5, annot6, annot7, annot8, annot9, annot10, annot11, annot12, annot13, annot14, annot15, annot16, annot17)

# plot
plot_gene_map(dna_segs=dna_segs, comparisons=NULL,
              annotations = annots, annotation_height = 3, annotation_cex = 0.6,
              main = 'E.coli Ter Locations Determined by BOWTIE2',
              dna_seg_scale=FALSE,  dna_seg_label_cex=0.8, scale = FALSE,
              arrow_head_len = 30000, gene_type = 'side_blocks')


###################################################################################################################
#
#                                                   REVERSE PLOT
#
###################################################################################################################


# reverse each dna_seg object
dna1_r <- reverse(dna1)
dna2_r <- reverse(dna2)
dna3_r <- reverse(dna3)

dna4_r <- reverse(dna4)
dna5_r <- reverse(dna5)
dna6_r <- reverse(dna6)

dna7_r <- reverse(dna7)
dna8_r <- reverse(dna8)
dna9_r <- reverse(dna9)

dna10_r <- reverse(dna10)
dna11_r <- reverse(dna11)
dna12_r <- reverse(dna12)
dna13_r <- reverse(dna13)
dna14_r <- reverse(dna14)

dna15_r <- reverse(dna15)
dna16_r <- reverse(dna16)
dna17_r <- reverse(dna17)


# sanity
dna1_r
dna2_r
dna3_r

dna4_r
dna5_r
dna6_r

dna7_r
dna8_r
dna9_r

dna10_r
dna11_r
dna12_r
dna13_r
dna14_r

dna15_r
dna16_r
dna17_r

# reset the strand column back to the original dna_seg 
## this will show the ter site in correct polarity
dna1_r$strand <- dna1$strand
dna2_r$strand <- dna2$strand
dna3_r$strand <- dna3$strand

dna4_r$strand <- dna4$strand
dna5_r$strand <- dna5$strand
dna6_r$strand <- dna6$strand

dna7_r$strand <- dna7$strand
dna8_r$strand <- dna8$strand
dna9_r$strand <- dna9$strand

dna10_r$strand <- dna10$strand
dna11_r$strand <- dna11$strand
dna12_r$strand <- dna12$strand
dna13_r$strand <- dna13$strand
dna14_r$strand <- dna14$strand

dna15_r$strand <- dna15$strand
dna16_r$strand <- dna16$strand
dna17_r$strand <- dna17$strand







#################################################

#       TURN ALL GENOMES INTO plusGrob or MinGrob

#################################################
###########################################
#     Using triangle shape (grob)
###########################################

## TEMPLATE
## Functions returning grobs.
## Creates a triangle for ter sites
triangleGrob <- function(gene, ...) {
  x <- c(gene$start, (gene$start+gene$end)/2, gene$end)
  y1 <- 0.5 + 0.5*gene$strand
  y <- c(y1, 0.5, y1)
  polygonGrob(x, y,gp=gpar(fill=gene$fill, col=gene$col, lty=gene$lty,
                           lwd=gene$lwd), default.units="native")
}

# template set...
# set the genetype in each df to triangleGrob
dna1_r$gene_type[1:10] <- 'triangleGrob'
dna4_r$gene_type[1:10] <- 'triangleGrob'
dna7_r$gene_type[1:10] <- 'triangleGrob'
dna12_r$gene_type[1:10] <- 'triangleGrob'
dna15_r$gene_type[1:10] <- 'triangleGrob'

##################################################################################
## Creates a permissible triangle shape grob for ter sites
# plus grobs
plusGrob <- function(gene, ...) {
  plus_x <- c((gene$start+gene$end)/2, (gene$start+gene$end)/2, gene$end+5000)
  plus_y <- c(1.1, -0.1, 0.5)
  polygonGrob(plus_x, plus_y,gp=gpar(fill=gene$fill, col=gene$col, lty=gene$lty,
                                     lwd=gene$lwd), default.units="native")
}

# min grob
minGrob <- function(gene, ...) {
  min_x <- c(gene$start-5000, (gene$start+gene$end)/2, (gene$start+gene$end)/2)
  min_y <- c(0.5, -0.1, 1.1)
  polygonGrob(min_x, min_y,gp=gpar(fill=gene$fill, col=gene$col, lty=gene$lty,
                                   lwd=gene$lwd), default.units="native")
}
###################################################################################



# Turn gene_type into plusGrob if + strand
dna1_r$gene_type[1:10][dna1_r$strand[1:10]==1] <- 'plusGrob'
dna2_r$gene_type[1:10][dna2_r$strand[1:10]==1] <- 'plusGrob'
dna3_r$gene_type[1:10][dna3_r$strand[1:10]==1] <- 'plusGrob'

dna4_r$gene_type[1:10][dna4_r$strand[1:10]==1] <- 'plusGrob'
dna5_r$gene_type[1:10][dna5_r$strand[1:10]==1] <- 'plusGrob'
dna6_r$gene_type[1:10][dna6_r$strand[1:10]==1] <- 'plusGrob'

dna7_r$gene_type[1:10][dna7_r$strand[1:10]==1] <- 'plusGrob'
dna8_r$gene_type[1:10][dna8_r$strand[1:10]==1] <- 'plusGrob'
dna9_r$gene_type[1:10][dna9_r$strand[1:10]==1] <- 'plusGrob'

dna10_r$gene_type[1:10][dna10_r$strand[1:10]==1] <- 'plusGrob'
dna11_r$gene_type[1:10][dna11_r$strand[1:10]==1] <- 'plusGrob'
dna12_r$gene_type[1:10][dna12_r$strand[1:10]==1] <- 'plusGrob'
dna13_r$gene_type[1:10][dna13_r$strand[1:10]==1] <- 'plusGrob'
dna14_r$gene_type[1:10][dna14_r$strand[1:10]==1] <- 'plusGrob'

dna15_r$gene_type[1:10][dna15_r$strand[1:10]==1] <- 'plusGrob'
dna16_r$gene_type[1:10][dna16_r$strand[1:10]==1] <- 'plusGrob'
dna17_r$gene_type[1:10][dna17_r$strand[1:10]==1] <- 'plusGrob'


# Turn gene_type into minGrob if - strand
dna1_r$gene_type[1:10][dna1_r$strand[1:10]==-1] <- 'minGrob'
dna2_r$gene_type[1:10][dna2_r$strand[1:10]==-1] <- 'minGrob'
dna3_r$gene_type[1:10][dna3_r$strand[1:10]==-1] <- 'minGrob'

dna4_r$gene_type[1:10][dna4_r$strand[1:10]==-1] <- 'minGrob'
dna5_r$gene_type[1:10][dna5_r$strand[1:10]==-1] <- 'minGrob'
dna6_r$gene_type[1:10][dna6_r$strand[1:10]==-1] <- 'minGrob'

dna7_r$gene_type[1:10][dna7_r$strand[1:10]==-1] <- 'minGrob'
dna8_r$gene_type[1:10][dna8_r$strand[1:10]==-1] <- 'minGrob'
dna9_r$gene_type[1:10][dna9_r$strand[1:10]==-1] <- 'minGrob'

dna10_r$gene_type[1:10][dna10_r$strand[1:10]==-1] <- 'minGrob'
dna11_r$gene_type[1:10][dna11_r$strand[1:10]==-1] <- 'minGrob'
dna12_r$gene_type[1:10][dna12_r$strand[1:10]==-1] <- 'minGrob'
dna13_r$gene_type[1:10][dna13_r$strand[1:10]==-1] <- 'minGrob'
dna14_r$gene_type[1:10][dna14_r$strand[1:10]==-1] <- 'minGrob'

dna15_r$gene_type[1:10][dna15_r$strand[1:10]==-1] <- 'minGrob'
dna16_r$gene_type[1:10][dna16_r$strand[1:10]==-1] <- 'minGrob'
dna17_r$gene_type[1:10][dna17_r$strand[1:10]==-1] <- 'minGrob'




# turn into list
dna_segs_r <- list(dna1_r,
               dna2_r,
               dna3_r,
               dna4_r,
               dna5_r,
               dna6_r,
               dna7_r,
               dna8_r,
               dna9_r,
               dna10_r,
               dna11_r,
               dna12_r,
               dna13_r,
               dna14_r,
               dna15_r,
               dna16_r,
               dna17_r)

names(dna_segs_r) <- names
dna_segs_r



## Calculating middle positions
mid_pos1_r <- middle(dna_segs_r[[1]])
mid_pos2_r <- middle(dna_segs_r[[2]])
mid_pos3_r <- middle(dna_segs_r[[3]])
mid_pos4_r <- middle(dna_segs_r[[4]])
mid_pos5_r <- middle(dna_segs_r[[5]])
mid_pos6_r <- middle(dna_segs_r[[6]])
mid_pos7_r <- middle(dna_segs_r[[7]])
mid_pos8_r <- middle(dna_segs_r[[8]])
mid_pos9_r <- middle(dna_segs_r[[9]])
mid_pos10_r <- middle(dna_segs_r[[10]])
mid_pos11_r <- middle(dna_segs_r[[11]])
mid_pos12_r <- middle(dna_segs_r[[12]])
mid_pos13_r <- middle(dna_segs_r[[13]])
mid_pos14_r <- middle(dna_segs_r[[14]])
mid_pos15_r <- middle(dna_segs_r[[15]])
mid_pos16_r <- middle(dna_segs_r[[16]])
mid_pos17_r <- middle(dna_segs_r[[17]])

# Annotations
annot1_r <- annotation(x1 = mid_pos1_r, text = dna_segs_r[[1]]$name, rot = 90, col = cols)
annot2_r <- annotation(x1 = mid_pos2_r, text = dna_segs_r[[2]]$name, rot = 90, col = cols)
annot3_r <- annotation(x1 = mid_pos3_r, text = dna_segs_r[[3]]$name, rot = 90, col = cols)
annot4_r <- annotation(x1 = mid_pos4_r, text = dna_segs_r[[4]]$name, rot = 90, col = cols)
annot5_r <- annotation(x1 = mid_pos5_r, text = dna_segs_r[[5]]$name, rot = 90, col = cols)
annot6_r <- annotation(x1 = mid_pos6_r, text = dna_segs_r[[6]]$name, rot = 90, col = cols)
annot7_r <- annotation(x1 = mid_pos7_r, text = dna_segs_r[[7]]$name, rot = 90, col = cols)
annot8_r <- annotation(x1 = mid_pos8_r, text = dna_segs_r[[8]]$name, rot = 90, col = cols)
annot9_r <- annotation(x1 = mid_pos9_r, text = dna_segs_r[[9]]$name, rot = 90, col = cols)
annot10_r <- annotation(x1 = mid_pos10_r, text = dna_segs_r[[10]]$name, rot = 90, col = cols)
annot11_r <- annotation(x1 = mid_pos11_r, text = dna_segs_r[[11]]$name, rot = 90, col = cols)
annot12_r <- annotation(x1 = mid_pos12_r, text = dna_segs_r[[12]]$name, rot = 90, col = cols)
annot13_r <- annotation(x1 = mid_pos13_r, text = dna_segs_r[[13]]$name, rot = 90, col = cols)
annot14_r <- annotation(x1 = mid_pos14_r, text = dna_segs_r[[14]]$name, rot = 90, col = cols)
annot15_r <- annotation(x1 = mid_pos15_r, text = dna_segs_r[[15]]$name, rot = 90, col = cols)
annot16_r <- annotation(x1 = mid_pos16_r, text = dna_segs_r[[16]]$name, rot = 90, col = cols)
annot17_r <- annotation(x1 = mid_pos17_r, text = dna_segs_r[[17]]$name, rot = 90, col = cols)

#list of annotation objects
annots_r <- list(annot1_r, 
                 annot2_r, 
                 annot3_r, 
                 annot4_r, 
                 annot5_r, 
                 annot6_r, 
                 annot7_r, 
                 annot8_r, 
                 annot9_r, 
                 annot10_r, 
                 annot11_r, 
                 annot12_r, 
                 annot13_r, 
                 annot14_r, 
                 annot15_r, 
                 annot16_r, 
                 annot17_r)





# plot all E.coli genomes figure with plus/minGrob
plot_gene_map(dna_segs=dna_segs_r, comparisons=NULL,
              annotations = annots_r, annotation_height = 3, annotation_cex = 0.6,
              main = 'E.coli Ter Locations Showing Permissive Directionality',
              dna_seg_scale=FALSE,  dna_seg_label_cex=0.8, scale = FALSE)




#################################################

#    5 PHYLOGROUP FIGURE

#################################################




# set the genetype in each df to plusGrob if + strand
dna1_r$gene_type[1:10][dna1_r$strand[1:10]==1] <- 'plusGrob'
dna4_r$gene_type[1:10][dna4_r$strand[1:10]==1] <- 'plusGrob'
dna7_r$gene_type[1:10][dna7_r$strand[1:10]==1] <- 'plusGrob'
dna12_r$gene_type[1:10][dna12_r$strand[1:10]==1] <- 'plusGrob'
dna15_r$gene_type[1:10][dna15_r$strand[1:10]==1] <- 'plusGrob'

# set the genetype in each df to minGrob if - strand
dna1_r$gene_type[1:10][dna1_r$strand[1:10]==-1] <- 'minGrob'
dna4_r$gene_type[1:10][dna4_r$strand[1:10]==-1] <- 'minGrob'
dna7_r$gene_type[1:10][dna7_r$strand[1:10]==-1] <- 'minGrob'
dna12_r$gene_type[1:10][dna12_r$strand[1:10]==-1] <- 'minGrob'
dna15_r$gene_type[1:10][dna15_r$strand[1:10]==-1] <- 'minGrob'




# plot 5 Phylogroup figure
plot_gene_map(dna_segs=list(dna1_r, dna4_r,dna7_r,dna12_r,dna15_r), comparisons=NULL,
              annotations = list(annot1_r, annot4_r, annot7_r,annot12_r,annot15_r), annotation_height = 3, annotation_cex = 0.9,
              main = 'Phylogroup Analysis of E.coli ter sites',
              dna_seg_scale=FALSE,  dna_seg_label_cex=1, scale = FALSE, gene_type = NULL,
              dna_seg_labels = c('MG1655 \nGroup A','APEC078 \nGroup B1','S88 \nGroup B2','UNM026 \n Group D','TW14359 \nGroup E'))



