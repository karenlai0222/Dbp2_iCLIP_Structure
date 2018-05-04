
##########
##### This is to do the metagene analysis based on the script from Siwen
##### The goal is to find out the distribution of Dbp2-binding sites on all Dbp2-binding mRNAs
#########
##### Notes: the number of bins for each region is porportional to their median length
### This script focuses on sample Dbp2-2
### The binding sites used here are re-filtered in Dec2018 by Nadia

library(Rsamtools)
library(GenomicRanges)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
library(GenomicAlignments)
sc.bs3 <- BSgenome.Scerevisiae.UCSC.sacCer3
yeast.length <- seqlengths(sc.bs3)

##### Read in Dbp2-binding mRNAs
dbp2.mRNA <- read.delim("~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/010418_Data_from_Nadia/091317_Dbp2_iCLIP_intersected_filtered_sites_count_mat_for_kept_mRNAs.txt", quote = "", as.is = T)

##### Read in the site file (this script only focuses on Dbp2-iCLIP replicate 1)
D22.all.sites <- read.delim("~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/010418_Data_from_Nadia/091317_iCLIP_Dbp22_intersected_reads_xlinking_sites_filtered.txt", quote = "", as.is = T)
D22.site.gr <- GRanges(seqnames = Rle(D22.all.sites$chr), ranges = IRanges(start = D22.all.sites$D2_xlinking_site, end = D22.all.sites$D2_xlinking_site), strand = Rle(D22.all.sites$strand), id = D22.all.sites$id, seqlengths = yeast.length)

##### Get the region information of these mRNAs
yn.utr.k <- read.delim("~/Desktop/Research_projects_2015-2016/Data_general_use/071016_Yassour_Nagalakshmi_outermost_UTR.txt", quote = "", as.is = T)
dbp2.mRNA.info <- yn.utr.k[which(yn.utr.k$id%in%rownames(dbp2.mRNA)), ]
dim(dbp2.mRNA.info)
### Separate them by strandness
dbp2.mRNA.info.pos <- dbp2.mRNA.info[which(dbp2.mRNA.info$strand=="+"), ]
dbp2.mRNA.info.neg <- dbp2.mRNA.info[which(dbp2.mRNA.info$strand=="-"), ]

##### Make GRanges of three genomic regions (5'UTR, ORF, 3'UTR) to estimate the width
utr5.gr <- GRanges(seqnames = Rle(dbp2.mRNA.info$chromosome), ranges = IRanges(start = dbp2.mRNA.info$utr5.start, end = dbp2.mRNA.info$utr5.end), strand = Rle(dbp2.mRNA.info$strand), id = dbp2.mRNA.info$id, seqlengths = yeast.length)
orf.gr <- GRanges(seqnames = Rle(dbp2.mRNA.info$chromosome), ranges = IRanges(start = dbp2.mRNA.info$orf.start, end = dbp2.mRNA.info$orf.end), strand = Rle(dbp2.mRNA.info$strand), id = dbp2.mRNA.info$id, seqlengths = yeast.length)
utr3.gr <- GRanges(seqnames = Rle(dbp2.mRNA.info$chromosome), ranges = IRanges(start = dbp2.mRNA.info$utr3.start, end = dbp2.mRNA.info$utr3.end), strand = Rle(dbp2.mRNA.info$strand), id = dbp2.mRNA.info$id, seqlengths = yeast.length)
### Get the width of each region
utr5.width <- width(utr5.gr)
summary(utr5.width) ### Median: 124
### Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
### 1.0    62.0   124.0   165.9   218.0   796.0 
orf.width <- width(orf.gr)
summary(orf.width) ### Median: 1368
### Min. 1st Qu.  Median   Mean   3rd Qu.  Max. 
###  78     855    1368    1718    2142   14730 
utr3.width <- width(utr3.gr)
summary(utr3.width) ### Median: 151
### Min. 1st Qu.  Median   Mean   3rd Qu.  Max.
### 6.0   104.0   151.0   188.9   233.0  1423.0 
#***# Based on the length, 5'UTR has 5 bins, ORF has 56 bins, 3'UTR has 6 bins
#***# Each bin has about 25 nt

#----------# This part does the meta-analysis in 10:80:10
##### Make the matrix for all three genomic regions (5'UTR, ORF, 3'UTR)
### 5'UTR
utr5.shift.mat.pos <- matrix(nrow = length(dbp2.mRNA.info.pos$id), ncol = 10)
colnames(utr5.shift.mat.pos) <- 1:10
rownames(utr5.shift.mat.pos) <- dbp2.mRNA.info.pos$id
for (i in 1:10){
    utr5.shift.mat.pos[,i] <- (dbp2.mRNA.info.pos$utr5.end - dbp2.mRNA.info.pos$utr5.start + 1)*0.1*i
}
utr5.pos.gr <- GRanges(seqnames = Rle(dbp2.mRNA.info.pos$chr), ranges = IRanges(start = dbp2.mRNA.info.pos$utr5.start-((dbp2.mRNA.info.pos$utr5.end - dbp2.mRNA.info.pos$utr5.start + 1)*0.1), end = dbp2.mRNA.info.pos$utr5.start-1), strand = Rle(dbp2.mRNA.info.pos$strand), id = dbp2.mRNA.info.pos$id, seqlengths = yeast.length)  

utr5.shift.mat.neg <- matrix(nrow = length(dbp2.mRNA.info.neg$id), ncol = 10)
colnames(utr5.shift.mat.neg) <- 1:10
rownames(utr5.shift.mat.neg) <- dbp2.mRNA.info.neg$id
for (i in 1:10){
    utr5.shift.mat.neg[,i] <- (dbp2.mRNA.info.neg$utr5.end - dbp2.mRNA.info.neg$utr5.start + 1)*0.1*i
}
utr5.neg.gr <- GRanges(seqnames = Rle(dbp2.mRNA.info.neg$chr), ranges = IRanges(end = dbp2.mRNA.info.neg$utr5.end+((dbp2.mRNA.info.neg$utr5.end - dbp2.mRNA.info.neg$utr5.start + 1)*0.1), start = dbp2.mRNA.info.neg$utr5.end+1), strand = Rle(dbp2.mRNA.info.neg$strand), id = dbp2.mRNA.info.neg$id, seqlengths = yeast.length)  
## Make the dataframes (countoverlaps) for iCLIP Dbp2-1, 5'UTR
D22.utr5.pos.count.mat <- matrix(nrow=length(utr5.pos.gr), ncol=10)
colnames(D22.utr5.pos.count.mat) <- 1:10
rownames(D22.utr5.pos.count.mat) <- mcols(utr5.pos.gr)$id
for(i in 1:10){
    temp.shift <- shift(utr5.pos.gr, shift=utr5.shift.mat.pos[,i])
    temp.count <- countOverlaps(temp.shift, D22.site.gr)
    D22.utr5.pos.count.mat[, i] <- temp.count
}

D22.utr5.neg.count.mat <- matrix(nrow=length(utr5.neg.gr), ncol=10)
colnames(D22.utr5.neg.count.mat) <- 1:10
rownames(D22.utr5.neg.count.mat) <- mcols(utr5.neg.gr)$id
for(i in 1:10){
    temp.shift <- shift(utr5.neg.gr, shift= -utr5.shift.mat.neg[,i])
    temp.count <- countOverlaps(temp.shift, D22.site.gr)
    D22.utr5.neg.count.mat[, i] <- temp.count
}

D22.utr5.count.mat <- rbind(D22.utr5.pos.count.mat, D22.utr5.neg.count.mat)

### ORF
orf.shift.mat.pos <- matrix(nrow = length(dbp2.mRNA.info.pos$id), ncol = 80)
colnames(orf.shift.mat.pos) <- 1:80
for (i in 1:80){
    orf.shift.mat.pos[,i] <- (dbp2.mRNA.info.pos$orf.end - dbp2.mRNA.info.pos$orf.start+1)*(1/80)*i
}
orf.pos.gr <- GRanges(seqnames = Rle(dbp2.mRNA.info.pos$chr), ranges = IRanges(start = dbp2.mRNA.info.pos$orf.start-((dbp2.mRNA.info.pos$orf.end - dbp2.mRNA.info.pos$orf.start+1)*(1/80)), end = dbp2.mRNA.info.pos$orf.start-1), strand = Rle(dbp2.mRNA.info.pos$strand), id = dbp2.mRNA.info.pos$id, seqlengths = yeast.length)  

orf.shift.mat.neg <- matrix(nrow = length(dbp2.mRNA.info.neg$id), ncol = 80)
colnames(orf.shift.mat.neg) <- 1:80
for (i in 1:80){
    orf.shift.mat.neg[,i] <- (dbp2.mRNA.info.neg$orf.end - dbp2.mRNA.info.neg$orf.start+1)*(1/80)*i
}
orf.neg.gr <- GRanges(seqnames = Rle(dbp2.mRNA.info.neg$chr), ranges = IRanges(end = dbp2.mRNA.info.neg$orf.end+((dbp2.mRNA.info.neg$orf.end - dbp2.mRNA.info.neg$orf.start+1)*(1/80)), start = dbp2.mRNA.info.neg$orf.end+1), strand = Rle(dbp2.mRNA.info.neg$strand), id = dbp2.mRNA.info.neg$id, seqlengths = yeast.length)  
## Make the dataframes (countoverlaps) for iCLIP Dbp2-1, ORF
D22.orf.pos.count.mat <- matrix(nrow=length(orf.pos.gr), ncol=80)
colnames(D22.orf.pos.count.mat) <- 1:80
for(i in 1:80){
    temp.shift <- shift(orf.pos.gr, shift=orf.shift.mat.pos[,i])
    temp.count <- countOverlaps(temp.shift, D22.site.gr)
    D22.orf.pos.count.mat[, i] <- temp.count
}

D22.orf.neg.count.mat <- matrix(nrow=length(orf.neg.gr), ncol=80)
colnames(D22.orf.neg.count.mat) <- 1:80
for(i in 1:80){
    temp.shift <- shift(orf.neg.gr, shift= -orf.shift.mat.neg[,i])
    temp.count <- countOverlaps(temp.shift, D22.site.gr)
    D22.orf.neg.count.mat[, i] <- temp.count
}

D22.orf.count.mat <- rbind(D22.orf.pos.count.mat, D22.orf.neg.count.mat)

### 3' UTR
utr3.shift.mat.pos <- matrix(nrow = length(dbp2.mRNA.info.pos$id), ncol = 10)
colnames(utr3.shift.mat.pos) <- 1:10
for (i in 1:10){
    utr3.shift.mat.pos[,i] <- (dbp2.mRNA.info.pos$utr3.end - dbp2.mRNA.info.pos$utr3.start + 1)*0.1*i
}
utr3.pos.gr <- GRanges(seqnames = Rle(dbp2.mRNA.info.pos$chr), ranges = IRanges(start = dbp2.mRNA.info.pos$utr3.start-((dbp2.mRNA.info.pos$utr3.end - dbp2.mRNA.info.pos$utr3.start +1)*(1/10)), end = dbp2.mRNA.info.pos$utr3.start-1), strand = Rle(dbp2.mRNA.info.pos$strand), id = dbp2.mRNA.info.pos$id, seqlengths = yeast.length)  

utr3.shift.mat.neg <- matrix(nrow = length(dbp2.mRNA.info.neg$id), ncol = 10)
colnames(utr3.shift.mat.neg) <- 1:10
for (i in 1:10){
    utr3.shift.mat.neg[,i] <- (dbp2.mRNA.info.neg$utr3.end - dbp2.mRNA.info.neg$utr3.start + 1)*0.1*i
}
utr3.neg.gr <- GRanges(seqnames = Rle(dbp2.mRNA.info.neg$chr), ranges = IRanges(end = dbp2.mRNA.info.neg$utr3.end+((dbp2.mRNA.info.neg$utr3.end - dbp2.mRNA.info.neg$utr3.start +1)*(1/10)), start = dbp2.mRNA.info.neg$utr3.end+1), strand = Rle(dbp2.mRNA.info.neg$strand), id = dbp2.mRNA.info.neg$id, seqlengths = yeast.length)  
## Make the dataframes (countoverlaps) for iCLIP Dbp2-1, 3'UTR
D22.utr3.pos.count.mat <- matrix(nrow=length(utr3.pos.gr), ncol=10)
colnames(D22.utr3.pos.count.mat) <- 1:10
for(i in 1:10){
    temp.shift <- shift(utr3.pos.gr, shift=utr3.shift.mat.pos[,i])
    temp.count <- countOverlaps(temp.shift, D22.site.gr)
    D22.utr3.pos.count.mat[, i] <- temp.count
}

D22.utr3.neg.count.mat <- matrix(nrow=length(utr3.neg.gr), ncol=10)
colnames(D22.utr3.neg.count.mat) <- 1:10
for(i in 1:10){
    temp.shift <- shift(utr3.neg.gr, shift= -utr3.shift.mat.neg[,i])
    temp.count <- countOverlaps(temp.shift, D22.site.gr)
    D22.utr3.neg.count.mat[, i] <- temp.count
}

D22.utr3.count.mat <- rbind(D22.utr3.pos.count.mat, D22.utr3.neg.count.mat)

##### Get the probability table
D22.utr5.prop.table <- prop.table(D22.utr5.count.mat, 1) ### 1 means calculating the proportion of the whole row
D22.orf.prop.table <- prop.table(D22.orf.count.mat, 1)
D22.utr3.prop.table <- prop.table(D22.utr3.count.mat, 1)
# ### If any cell is NaN, put 0
# ## 5'UTR
# D22.utr5.prop.table.mod <- D22.utr5.prop.table
# utr5.nan <- is.nan(D22.utr5.prop.table.mod)
# D22.utr5.prop.table.mod[utr5.nan] <- 0
# ## ORF
# D22.orf.prop.table.mod <- D22.orf.prop.table
# orf.nan <- is.nan(D22.orf.prop.table.mod)
# D22.orf.prop.table.mod[orf.nan] <- 0
# ## 3'UTR
# D22.utr3.prop.table.mod <- D22.utr3.prop.table
# utr3.nan <- is.nan(D22.utr3.prop.table.mod)
# D22.utr3.prop.table.mod[utr3.nan] <- 0

##### Plot the 3 regions individually
plot(apply(D22.utr5.prop.table, 2, mean, na.rm=T), type="l", lwd=2, col="red", ylim=c(0, 0.5), main = "iCLIP_D22, 5'UTR")
# plot(apply(D22.utr5.prop.table.mod, 2, mean, na.rm=T), type="l", lwd=2, col="red", ylim=c(0, 0.5), main = "iCLIP_D22, 5'UTR")
# plot(apply(D22.utr5.prop.table, 2, median, na.rm=T), type="l", lwd=2, col="red", ylim=c(0, 0.5), main = "iCLIP_D22, 5'UTR")
### Using median, the trend is a little sharper
plot(apply(D22.utr3.prop.table, 2, mean, na.rm=T), type="l", lwd=2, col="red", ylim=c(0, 0.5), main = "iCLIP_D22, 3'UTR")
### The medians for 3'UTRs are zero! Don't use median for the downstream analysis
plot(apply(D22.orf.prop.table, 2, mean, na.rm=T), type="l", lwd=2, col="red", ylim=c(0, 0.5), main = "iCLIP_D22, ORF")

##### Now combine the 3 regions
D22.whole.mat <- cbind(D22.utr5.count.mat, D22.orf.count.mat)
D22.whole.mat <- cbind(D22.whole.mat, D22.utr3.count.mat)
rownames(D22.whole.mat) <- c(mcols(utr5.pos.gr)$id, mcols(utr5.neg.gr)$id)
D22.whole.prop.table <- prop.table(D22.whole.mat, 1)
### Output the count table and the probability table
write.table(D22.whole.mat, file = "~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/010418_Data_from_Nadia/010518_iCLIP_Dbp22_all_filtered_sites_in_bound_mRNAs_count_mat_100bins.txt", quote = F, sep = "\t")
write.table(D22.whole.prop.table, file = "~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/010418_Data_from_Nadia/010518_iCLIP_Dbp22_all_filtered_sites_in_bound_mRNA_prop_table_100bins.txt", quote = F, sep = "\t")
### Plot the whole transcript
setwd("~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/091417_Dbp2_meta-analysis_after_revised_filtering/")
tiff("091417_Dbp2_iCLIP_all_sites_metagene_in_bound_mRNA_prop_100bins_res600.tiff", height = 4, width = 5, units = 'in', res = 600)
par(lwd = 2)
plot(0.5:99.5, apply(D22.whole.prop.table, 2, mean, na.rm=T), cex.lab = 1.2, cex.axis = 1.1, type="l", lwd=4, col="blue3", yaxt="n", ylim = c(0, 0.03), xaxt = "n", xlab = "Relative transcript positions", ylab = "Dbp2 occupancy", main = "iCLIP_D22, Dbp2-binding mRNAs, all sites")
axis(1, at = seq(0, 100, by = 1), las = 2, lwd = 2)
axis(2, at = seq(0, 0.08, by = 0.01), lwd = 2)
abline(v=10, lty=2, lwd = 3)
abline(v=90, lty=2, lwd = 3)
dev.off()

#################################################################################################################

#----------# This part does the meta-analysis in 5:55:6
##### Make the matrix for all three genomic regions (5'UTR, ORF, 3'UTR)
### 5'UTR
utr5.shift.mat.pos <- matrix(nrow = length(dbp2.mRNA.info.pos$id), ncol = 5)
colnames(utr5.shift.mat.pos) <- 1:5
rownames(utr5.shift.mat.pos) <- dbp2.mRNA.info.pos$id
for (i in 1:5){
    utr5.shift.mat.pos[,i] <- (dbp2.mRNA.info.pos$utr5.end - dbp2.mRNA.info.pos$utr5.start + 1)*(1/5)*i
}
utr5.pos.gr <- GRanges(seqnames = Rle(dbp2.mRNA.info.pos$chr), ranges = IRanges(start = dbp2.mRNA.info.pos$utr5.start-((dbp2.mRNA.info.pos$utr5.end - dbp2.mRNA.info.pos$utr5.start + 1)*(1/5)), end = dbp2.mRNA.info.pos$utr5.start-1), strand = Rle(dbp2.mRNA.info.pos$strand), id = dbp2.mRNA.info.pos$id, seqlengths = yeast.length)  

utr5.shift.mat.neg <- matrix(nrow = length(dbp2.mRNA.info.neg$id), ncol = 5)
colnames(utr5.shift.mat.neg) <- 1:5
rownames(utr5.shift.mat.neg) <- dbp2.mRNA.info.neg$id
for (i in 1:5){
    utr5.shift.mat.neg[,i] <- (dbp2.mRNA.info.neg$utr5.end - dbp2.mRNA.info.neg$utr5.start + 1)*(1/5)*i
}
utr5.neg.gr <- GRanges(seqnames = Rle(dbp2.mRNA.info.neg$chr), ranges = IRanges(end = dbp2.mRNA.info.neg$utr5.end+((dbp2.mRNA.info.neg$utr5.end - dbp2.mRNA.info.neg$utr5.start + 1)*(1/5)), start = dbp2.mRNA.info.neg$utr5.end+1), strand = Rle(dbp2.mRNA.info.neg$strand), id = dbp2.mRNA.info.neg$id, seqlengths = yeast.length)  
## Make the dataframes (countoverlaps) for iCLIP Dbp2-1, 5'UTR
D22.utr5.pos.count.mat <- matrix(nrow=length(utr5.pos.gr), ncol=5)
colnames(D22.utr5.pos.count.mat) <- 1:5
rownames(D22.utr5.pos.count.mat) <- mcols(utr5.pos.gr)$id
for(i in 1:5){
    temp.shift <- shift(utr5.pos.gr, shift=utr5.shift.mat.pos[,i])
    temp.count <- countOverlaps(temp.shift, D22.site.gr)
    D22.utr5.pos.count.mat[, i] <- temp.count
}

D22.utr5.neg.count.mat <- matrix(nrow=length(utr5.neg.gr), ncol=5)
colnames(D22.utr5.neg.count.mat) <- 1:5
rownames(D22.utr5.neg.count.mat) <- mcols(utr5.neg.gr)$id
for(i in 1:5){
    temp.shift <- shift(utr5.neg.gr, shift= -utr5.shift.mat.neg[,i])
    temp.count <- countOverlaps(temp.shift, D22.site.gr)
    D22.utr5.neg.count.mat[, i] <- temp.count
}

D22.utr5.count.mat <- rbind(D22.utr5.pos.count.mat, D22.utr5.neg.count.mat)

### ORF
orf.shift.mat.pos <- matrix(nrow = length(dbp2.mRNA.info.pos$id), ncol = 55)
colnames(orf.shift.mat.pos) <- 1:55
for (i in 1:55){
    orf.shift.mat.pos[,i] <- (dbp2.mRNA.info.pos$orf.end - dbp2.mRNA.info.pos$orf.start+1)*(1/55)*i
}
orf.pos.gr <- GRanges(seqnames = Rle(dbp2.mRNA.info.pos$chr), ranges = IRanges(start = dbp2.mRNA.info.pos$orf.start-((dbp2.mRNA.info.pos$orf.end - dbp2.mRNA.info.pos$orf.start+1)*(1/55)), end = dbp2.mRNA.info.pos$orf.start-1), strand = Rle(dbp2.mRNA.info.pos$strand), id = dbp2.mRNA.info.pos$id, seqlengths = yeast.length)  

orf.shift.mat.neg <- matrix(nrow = length(dbp2.mRNA.info.neg$id), ncol = 55)
colnames(orf.shift.mat.neg) <- 1:55
for (i in 1:55){
    orf.shift.mat.neg[,i] <- (dbp2.mRNA.info.neg$orf.end - dbp2.mRNA.info.neg$orf.start+1)*(1/55)*i
}
orf.neg.gr <- GRanges(seqnames = Rle(dbp2.mRNA.info.neg$chr), ranges = IRanges(end = dbp2.mRNA.info.neg$orf.end+((dbp2.mRNA.info.neg$orf.end - dbp2.mRNA.info.neg$orf.start+1)*(1/55)), start = dbp2.mRNA.info.neg$orf.end+1), strand = Rle(dbp2.mRNA.info.neg$strand), id = dbp2.mRNA.info.neg$id, seqlengths = yeast.length)  
## Make the dataframes (countoverlaps) for iCLIP Dbp2-1, ORF
D22.orf.pos.count.mat <- matrix(nrow=length(orf.pos.gr), ncol=55)
colnames(D22.orf.pos.count.mat) <- 1:55
for(i in 1:55){
    temp.shift <- shift(orf.pos.gr, shift=orf.shift.mat.pos[,i])
    temp.count <- countOverlaps(temp.shift, D22.site.gr)
    D22.orf.pos.count.mat[, i] <- temp.count
}

D22.orf.neg.count.mat <- matrix(nrow=length(orf.neg.gr), ncol=55)
colnames(D22.orf.neg.count.mat) <- 1:55
for(i in 1:55){
    temp.shift <- shift(orf.neg.gr, shift= -orf.shift.mat.neg[,i])
    temp.count <- countOverlaps(temp.shift, D22.site.gr)
    D22.orf.neg.count.mat[, i] <- temp.count
}

D22.orf.count.mat <- rbind(D22.orf.pos.count.mat, D22.orf.neg.count.mat)

### 3' UTR
utr3.shift.mat.pos <- matrix(nrow = length(dbp2.mRNA.info.pos$id), ncol = 6)
colnames(utr3.shift.mat.pos) <- 1:6
for (i in 1:6){
    utr3.shift.mat.pos[,i] <- (dbp2.mRNA.info.pos$utr3.end - dbp2.mRNA.info.pos$utr3.start + 1)*(1/6)*i
}
utr3.pos.gr <- GRanges(seqnames = Rle(dbp2.mRNA.info.pos$chr), ranges = IRanges(start = dbp2.mRNA.info.pos$utr3.start-((dbp2.mRNA.info.pos$utr3.end - dbp2.mRNA.info.pos$utr3.start +1)*(1/6)), end = dbp2.mRNA.info.pos$utr3.start-1), strand = Rle(dbp2.mRNA.info.pos$strand), id = dbp2.mRNA.info.pos$id, seqlengths = yeast.length)  

utr3.shift.mat.neg <- matrix(nrow = length(dbp2.mRNA.info.neg$id), ncol = 6)
colnames(utr3.shift.mat.neg) <- 1:6
for (i in 1:6){
    utr3.shift.mat.neg[,i] <- (dbp2.mRNA.info.neg$utr3.end - dbp2.mRNA.info.neg$utr3.start + 1)*(1/6)*i
}
utr3.neg.gr <- GRanges(seqnames = Rle(dbp2.mRNA.info.neg$chr), ranges = IRanges(end = dbp2.mRNA.info.neg$utr3.end+((dbp2.mRNA.info.neg$utr3.end - dbp2.mRNA.info.neg$utr3.start +1)*(1/6)), start = dbp2.mRNA.info.neg$utr3.end+1), strand = Rle(dbp2.mRNA.info.neg$strand), id = dbp2.mRNA.info.neg$id, seqlengths = yeast.length)  
## Make the dataframes (countoverlaps) for iCLIP Dbp2-1, 3'UTR
D22.utr3.pos.count.mat <- matrix(nrow=length(utr3.pos.gr), ncol=6)
colnames(D22.utr3.pos.count.mat) <- 1:6
for(i in 1:6){
    temp.shift <- shift(utr3.pos.gr, shift=utr3.shift.mat.pos[,i])
    temp.count <- countOverlaps(temp.shift, D22.site.gr)
    D22.utr3.pos.count.mat[, i] <- temp.count
}

D22.utr3.neg.count.mat <- matrix(nrow=length(utr3.neg.gr), ncol=6)
colnames(D22.utr3.neg.count.mat) <- 1:6
for(i in 1:6){
    temp.shift <- shift(utr3.neg.gr, shift= -utr3.shift.mat.neg[,i])
    temp.count <- countOverlaps(temp.shift, D22.site.gr)
    D22.utr3.neg.count.mat[, i] <- temp.count
}

D22.utr3.count.mat <- rbind(D22.utr3.pos.count.mat, D22.utr3.neg.count.mat)

##### Get the probability table
D22.utr5.prop.table <- prop.table(D22.utr5.count.mat, 1) ### 1 means calculating the proportion of the whole row
D22.orf.prop.table <- prop.table(D22.orf.count.mat, 1)
D22.utr3.prop.table <- prop.table(D22.utr3.count.mat, 1)
# ### If any cell is NaN, put 0
# ## 5'UTR
# D22.utr5.prop.table.mod <- D22.utr5.prop.table
# utr5.nan <- is.nan(D22.utr5.prop.table.mod)
# D22.utr5.prop.table.mod[utr5.nan] <- 0
# ## ORF
# D22.orf.prop.table.mod <- D22.orf.prop.table
# orf.nan <- is.nan(D22.orf.prop.table.mod)
# D22.orf.prop.table.mod[orf.nan] <- 0
# ## 3'UTR
# D22.utr3.prop.table.mod <- D22.utr3.prop.table
# utr3.nan <- is.nan(D22.utr3.prop.table.mod)
# D22.utr3.prop.table.mod[utr3.nan] <- 0

##### Plot the 3 regions individually
plot(apply(D22.utr5.prop.table, 2, mean, na.rm=T), type="l", lwd=2, col="red", ylim=c(0, 0.5), main = "iCLIP_D22, 5'UTR")
# plot(apply(D22.utr5.prop.table.mod, 2, mean, na.rm=T), type="l", lwd=2, col="red", ylim=c(0, 0.5), main = "iCLIP_D22, 5'UTR")
# plot(apply(D22.utr5.prop.table, 2, median, na.rm=T), type="l", lwd=2, col="red", ylim=c(0, 0.5), main = "iCLIP_D22, 5'UTR")
### Using median, the trend is a little sharper
plot(apply(D22.utr3.prop.table, 2, mean, na.rm=T), type="l", lwd=2, col="red", ylim=c(0, 0.5), main = "iCLIP_D22, 3'UTR")
### The medians for 3'UTRs are zero! Don't use median for the downstream analysis
plot(apply(D22.orf.prop.table, 2, mean, na.rm=T), type="l", lwd=2, col="red", ylim=c(0, 0.5), main = "iCLIP_D22, ORF")

##### Now combine the 3 regions
D22.whole.mat <- cbind(D22.utr5.count.mat, D22.orf.count.mat)
D22.whole.mat <- cbind(D22.whole.mat, D22.utr3.count.mat)
rownames(D22.whole.mat) <- c(mcols(utr5.pos.gr)$id, mcols(utr5.neg.gr)$id)
D22.whole.prop.table <- prop.table(D22.whole.mat, 1)
### Output the count table and the probability table
write.table(D22.whole.mat, file = "~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/091417_Dbp2_meta-analysis_after_revised_filtering/091417_iCLIP_Dbp22_all_filtered_sites_in_bound_mRNAs_count_mat_66bins.txt", quote = F, sep = "\t")
write.table(D22.whole.prop.table, file = "~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/091417_Dbp2_meta-analysis_after_revised_filtering/091417_iCLIP_Dbp22_all_filtered_sites_in_bound_mRNA_prop_table_66bins", quote = F, sep = "\t")
### Plot the whole transcript
setwd("~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/091417_Dbp2_meta-analysis_after_revised_filtering/")
tiff("091417_Dbp22_iCLIP_all_sites_metagene_in_bound_mRNA_prop_66bins_res600.tiff", height = 4, width = 5, units = 'in', res = 600)
par(lwd = 2)
plot(0.5:65.5, apply(D22.whole.prop.table, 2, mean, na.rm=T), cex.lab = 1.2, cex.axis = 1.1, type="l", lwd=4, col="blue3", yaxt="n", ylim = c(0, 0.04), xaxt = "n", xlab = "Relative transcript positions", ylab = "Dbp2 occupancy", main = "iCLIP_D22, Dbp2-binding mRNAs, all sites")
axis(1, at = seq(0, 66, by = 1), las = 2, lwd = 2)
axis(2, at = seq(0, 0.08, by = 0.01), lwd = 2)
abline(v=5, lty=2, lwd = 3)
abline(v=60, lty=2, lwd = 3)
dev.off()




