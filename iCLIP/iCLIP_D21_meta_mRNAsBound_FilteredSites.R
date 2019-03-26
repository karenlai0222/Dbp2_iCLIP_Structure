
##########
##### This is to do the metagene analysis
##### The goal is to find out the distribution of Dbp2-binding sites on all Dbp2-binding mRNAs
#########
### This script focuses on sample Dbp2-1
### This script is run on cluster because the files are big

library(Rsamtools)
library(GenomicRanges)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
library(GenomicAlignments)
sc.bs3 <- BSgenome.Scerevisiae.UCSC.sacCer3
yeast.length <- seqlengths(sc.bs3)

##### Read in Dbp2-binding mRNAs
dbp2.mRNA <- read.delim("060718_Dbp2iCLIP_SitesFromIntersectedReads_countMat_for_kept_mRNAs.txt",
                        quote = "", as.is = T)

##### Read in the site file (this script only focuses on Dbp2-iCLIP replicate 1)
D21.all.sites <- read.delim("060618_iCLIP_Dbp21_intersectedReads_allxlinkingSites.txt",
                            quote = "", as.is = T)
D21.site.gr <- GRanges(seqnames = Rle(D21.all.sites$chr),
                       ranges = IRanges(start = D21.all.sites$D1_xlinking_site, end = D21.all.sites$D1_xlinking_site),
                       strand = Rle(D21.all.sites$strand), id = D21.all.sites$id,
                       seqlengths = yeast.length)

##### Get the region information of these mRNAs
yn.utr.k <- read.delim("071016_Yassour_Nagalakshmi_outermost_UTR.txt", quote = "", as.is = T)
dbp2.mRNA.info <- yn.utr.k[which(yn.utr.k$id%in%rownames(dbp2.mRNA)), ]
dim(dbp2.mRNA.info)
### Separate them by strandness
dbp2.mRNA.info.pos <- dbp2.mRNA.info[which(dbp2.mRNA.info$strand=="+"), ]
dbp2.mRNA.info.neg <- dbp2.mRNA.info[which(dbp2.mRNA.info$strand=="-"), ]

##### Make GRanges of three genomic regions (5'UTR, ORF, 3'UTR)
utr5.gr <- GRanges(seqnames = Rle(dbp2.mRNA.info$chromosome),
                   ranges = IRanges(start = dbp2.mRNA.info$utr5.start, end = dbp2.mRNA.info$utr5.end),
                   strand = Rle(dbp2.mRNA.info$strand), id = dbp2.mRNA.info$id,
                   seqlengths = yeast.length)
orf.gr <- GRanges(seqnames = Rle(dbp2.mRNA.info$chromosome),
                  ranges = IRanges(start = dbp2.mRNA.info$orf.start, end = dbp2.mRNA.info$orf.end),
                  strand = Rle(dbp2.mRNA.info$strand), id = dbp2.mRNA.info$id,
                  seqlengths = yeast.length)
utr3.gr <- GRanges(seqnames = Rle(dbp2.mRNA.info$chromosome),
                   ranges = IRanges(start = dbp2.mRNA.info$utr3.start, end = dbp2.mRNA.info$utr3.end),
                   strand = Rle(dbp2.mRNA.info$strand), id = dbp2.mRNA.info$id,
                   seqlengths = yeast.length)

#----------# This part does the meta-analysis in 10:80:10
##### Make the matrix for all three genomic regions (5'UTR, ORF, 3'UTR)
### 5'UTR
utr5.shift.mat.pos <- matrix(nrow = length(dbp2.mRNA.info.pos$id), ncol = 10)
colnames(utr5.shift.mat.pos) <- 1:10
rownames(utr5.shift.mat.pos) <- dbp2.mRNA.info.pos$id
for (i in 1:10){
    utr5.shift.mat.pos[,i] <- (dbp2.mRNA.info.pos$utr5.end - dbp2.mRNA.info.pos$utr5.start + 1)*0.1*i
}
utr5.pos.gr <- GRanges(seqnames = Rle(dbp2.mRNA.info.pos$chr),
                       ranges = IRanges(start = dbp2.mRNA.info.pos$utr5.start-((dbp2.mRNA.info.pos$utr5.end - dbp2.mRNA.info.pos$utr5.start + 1)*0.1),
                                        end = dbp2.mRNA.info.pos$utr5.start-1),
                       strand = Rle(dbp2.mRNA.info.pos$strand),
                       id = dbp2.mRNA.info.pos$id, seqlengths = yeast.length)  

utr5.shift.mat.neg <- matrix(nrow = length(dbp2.mRNA.info.neg$id), ncol = 10)
colnames(utr5.shift.mat.neg) <- 1:10
rownames(utr5.shift.mat.neg) <- dbp2.mRNA.info.neg$id
for (i in 1:10){
    utr5.shift.mat.neg[,i] <- (dbp2.mRNA.info.neg$utr5.end - dbp2.mRNA.info.neg$utr5.start + 1)*0.1*i
}
utr5.neg.gr <- GRanges(seqnames = Rle(dbp2.mRNA.info.neg$chr),
                       ranges = IRanges(end = dbp2.mRNA.info.neg$utr5.end+((dbp2.mRNA.info.neg$utr5.end - dbp2.mRNA.info.neg$utr5.start + 1)*0.1),
                                        start = dbp2.mRNA.info.neg$utr5.end+1),
                       strand = Rle(dbp2.mRNA.info.neg$strand), id = dbp2.mRNA.info.neg$id,
                       seqlengths = yeast.length)  
## Make the dataframes (countoverlaps) for iCLIP Dbp2-1, 5'UTR
D21.utr5.pos.count.mat <- matrix(nrow=length(utr5.pos.gr), ncol=10)
colnames(D21.utr5.pos.count.mat) <- 1:10
rownames(D21.utr5.pos.count.mat) <- mcols(utr5.pos.gr)$id
for(i in 1:10){
    temp.shift <- shift(utr5.pos.gr, shift=utr5.shift.mat.pos[,i])
    temp.count <- countOverlaps(temp.shift, D21.site.gr)
    D21.utr5.pos.count.mat[, i] <- temp.count
}

D21.utr5.neg.count.mat <- matrix(nrow=length(utr5.neg.gr), ncol=10)
colnames(D21.utr5.neg.count.mat) <- 1:10
rownames(D21.utr5.neg.count.mat) <- mcols(utr5.neg.gr)$id
for(i in 1:10){
    temp.shift <- shift(utr5.neg.gr, shift= -utr5.shift.mat.neg[,i])
    temp.count <- countOverlaps(temp.shift, D21.site.gr)
    D21.utr5.neg.count.mat[, i] <- temp.count
}

D21.utr5.count.mat <- rbind(D21.utr5.pos.count.mat, D21.utr5.neg.count.mat)

### ORF
orf.shift.mat.pos <- matrix(nrow = length(dbp2.mRNA.info.pos$id), ncol = 80)
colnames(orf.shift.mat.pos) <- 1:80
for (i in 1:80){
    orf.shift.mat.pos[,i] <- (dbp2.mRNA.info.pos$orf.end - dbp2.mRNA.info.pos$orf.start+1)*(1/80)*i
}
orf.pos.gr <- GRanges(seqnames = Rle(dbp2.mRNA.info.pos$chr),
                      ranges = IRanges(start = dbp2.mRNA.info.pos$orf.start-((dbp2.mRNA.info.pos$orf.end - dbp2.mRNA.info.pos$orf.start+1)*(1/80)),
                                       end = dbp2.mRNA.info.pos$orf.start-1),
                      strand = Rle(dbp2.mRNA.info.pos$strand), id = dbp2.mRNA.info.pos$id,
                      seqlengths = yeast.length)  

orf.shift.mat.neg <- matrix(nrow = length(dbp2.mRNA.info.neg$id), ncol = 80)
colnames(orf.shift.mat.neg) <- 1:80
for (i in 1:80){
    orf.shift.mat.neg[,i] <- (dbp2.mRNA.info.neg$orf.end - dbp2.mRNA.info.neg$orf.start+1)*(1/80)*i
}
orf.neg.gr <- GRanges(seqnames = Rle(dbp2.mRNA.info.neg$chr),
                      ranges = IRanges(end = dbp2.mRNA.info.neg$orf.end+((dbp2.mRNA.info.neg$orf.end - dbp2.mRNA.info.neg$orf.start+1)*(1/80)),
                                       start = dbp2.mRNA.info.neg$orf.end+1),
                      strand = Rle(dbp2.mRNA.info.neg$strand), id = dbp2.mRNA.info.neg$id,
                      seqlengths = yeast.length)  
## Make the dataframes (countoverlaps) for iCLIP Dbp2-1, ORF
D21.orf.pos.count.mat <- matrix(nrow=length(orf.pos.gr), ncol=80)
colnames(D21.orf.pos.count.mat) <- 1:80
for(i in 1:80){
    temp.shift <- shift(orf.pos.gr, shift=orf.shift.mat.pos[,i])
    temp.count <- countOverlaps(temp.shift, D21.site.gr)
    D21.orf.pos.count.mat[, i] <- temp.count
}

D21.orf.neg.count.mat <- matrix(nrow=length(orf.neg.gr), ncol=80)
colnames(D21.orf.neg.count.mat) <- 1:80
for(i in 1:80){
    temp.shift <- shift(orf.neg.gr, shift= -orf.shift.mat.neg[,i])
    temp.count <- countOverlaps(temp.shift, D21.site.gr)
    D21.orf.neg.count.mat[, i] <- temp.count
}

D21.orf.count.mat <- rbind(D21.orf.pos.count.mat, D21.orf.neg.count.mat)

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
D21.utr3.pos.count.mat <- matrix(nrow=length(utr3.pos.gr), ncol=10)
colnames(D21.utr3.pos.count.mat) <- 1:10
for(i in 1:10){
    temp.shift <- shift(utr3.pos.gr, shift=utr3.shift.mat.pos[,i])
    temp.count <- countOverlaps(temp.shift, D21.site.gr)
    D21.utr3.pos.count.mat[, i] <- temp.count
}

D21.utr3.neg.count.mat <- matrix(nrow=length(utr3.neg.gr), ncol=10)
colnames(D21.utr3.neg.count.mat) <- 1:10
for(i in 1:10){
    temp.shift <- shift(utr3.neg.gr, shift= -utr3.shift.mat.neg[,i])
    temp.count <- countOverlaps(temp.shift, D21.site.gr)
    D21.utr3.neg.count.mat[, i] <- temp.count
}

D21.utr3.count.mat <- rbind(D21.utr3.pos.count.mat, D21.utr3.neg.count.mat)

##### Get the probability table
D21.utr5.prop.table <- prop.table(D21.utr5.count.mat, 1) ### 1 means calculating the proportion of the whole row
D21.orf.prop.table <- prop.table(D21.orf.count.mat, 1)
D21.utr3.prop.table <- prop.table(D21.utr3.count.mat, 1)

##### Plot the 3 regions individually
plot(apply(D21.utr5.prop.table, 2, mean, na.rm=T), type="l", lwd=2, col="red", ylim=c(0, 0.5), main = "iCLIP_D21, 5'UTR")
### Using median, the trend is a little sharper
plot(apply(D21.utr3.prop.table, 2, mean, na.rm=T), type="l", lwd=2, col="red", ylim=c(0, 0.5), main = "iCLIP_D21, 3'UTR")
### The medians for 3'UTRs are zero! Don't use median for the downstream analysis
plot(apply(D21.orf.prop.table, 2, mean, na.rm=T), type="l", lwd=2, col="red", ylim=c(0, 0.5), main = "iCLIP_D21, ORF")

##### Now combine the 3 regions
D21.whole.mat <- cbind(D21.utr5.count.mat, D21.orf.count.mat)
D21.whole.mat <- cbind(D21.whole.mat, D21.utr3.count.mat)
rownames(D21.whole.mat) <- c(mcols(utr5.pos.gr)$id, mcols(utr5.neg.gr)$id)
D21.whole.prop.table <- prop.table(D21.whole.mat, 1)
### Output the count table and the probability table
write.table(D21.whole.mat, file = "060718_iCLIP_Dbp21_FilteredSites_inBoundmRNAs_countMat100bins.txt",
            quote = F, sep = "\t")
write.table(D21.whole.prop.table,
            file = "060718_iCLIP_Dbp21_all_FilteredSites_inBoundmRNA_propTable100bins.txt",
            quote = F, sep = "\t")
### Plot the whole transcript
D21.whole.prop.table <- read.delim("060718_iCLIP_Dbp21_all_FilteredSites_inBoundmRNA_propTable100bins.txt",
                                   as.is = T)
par(lwd = 2)
plot(0.5:99.5, apply(D21.whole.prop.table, 2, mean, na.rm=T), cex.lab = 1.2, cex.axis = 1.1,
     type="l", lwd=4, col="blue3", yaxt="n", ylim = c(0, 0.03), xaxt = "n",
     xlab = "Relative transcript positions", ylab = "Dbp2 occupancy",
     main = "iCLIP_D21, Dbp2-binding mRNAs, all sites")
axis(1, at = seq(0, 100, by = 1), las = 2, lwd = 2)
axis(2, at = seq(0, 0.08, by = 0.01), lwd = 2)
abline(v=10, lty=2, lwd = 3)
abline(v=90, lty=2, lwd = 3)
dev.off()




