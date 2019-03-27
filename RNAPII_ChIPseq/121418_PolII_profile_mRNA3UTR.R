
library(Rsamtools)
library(GenomicRanges)
library(GenomicAlignments)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
sc.bs3 <- BSgenome.Scerevisiae.UCSC.sacCer3
yeast.length <- seqlengths(sc.bs3)

#----------# 
##### Read the resulted matix
mRNA3UTRend.mat <- read.delim("121718_PolII_mRNAs_3UTRend_matForPlot.txt",
                              comment.char = "#", as.is = T)
wt.mRNAends.mat <- mRNA3UTRend.mat[, grep("WT", colnames(mRNA3UTRend.mat))]
dbp2.mRNAends.mat <- mRNA3UTRend.mat[, grep("dbp2", colnames(mRNA3UTRend.mat))]

### Plotting
wt.mRNAends.mean <- apply(wt.mRNAends.mat, 2, mean, na.rm = T)
wt.mRNAends.mean.norm <- (wt.mRNAends.mean/sum(wt.mRNAends.mean))*100
dbp2.mRNAends.mean <- apply(dbp2.mRNAends.mat, 2, mean, na.rm = T)
dbp2.mRNAends.mean.norm <- (dbp2.mRNAends.mean/sum(dbp2.mRNAends.mean))*100
x = seq(-490, 500, by = 10)
tiff("121718_PolII_nonOverlap_mRNA3UTR_withColors.tiff", height = 4, width = 5, units = 'in', res = 600)
par(lwd = 1, ps = 12, cex = 1)
plot(x = x, y = wt.mRNAends.mean.norm, type = "l", col = "skyblue", lwd = 2, xlim = c(-458,470),
     ylim = c(0.7, 1.2), xaxt = "n", yaxt = "n",
     xlab = "distance to mRNA 3' UTR ends (nt)", ylab = "Normalized Pol II occupancy",
     main = "RNA Pol II profile around the end of mRNA 3' UTRs")
polygon(c(-480, x, 510), c(0,wt.mRNAends.mean.norm,0), col = adjustcolor("skyblue", alpha.f = 0.4),
        border = NA)
lines(x = x, y = dbp2.mRNAends.mean.norm, type = "l", col = "pink", lwd = 2)
polygon(c(-480, x, 510), c(0,dbp2.mRNAends.mean.norm,0), col = adjustcolor("pink", alpha.f = 0.4),
        border = NA)
axis(1, at = seq(-600, 600, by = 200), las = 1, lwd = 2)
axis(2, at = seq(0, 1.4, by = 0.2), las = 2, lwd = 2)
abline(v=0, lty=2, lwd = 1)
legend("topright", legend = c("1", "2"), cex=1.2, text.font=8,
       y.intersp = 1.5, box.lty=0, bg = "transparent", fill = c("skyblue", "pink"))
dev.off()

##### Read the resulted matix of read-through mRNAs only
rt.mat <- read.delim("121718_PolII_readthrough_mRNAs_3UTRend_matForPlot.txt",
                     comment.char = "#", as.is = T)
wt.rt.mat <- rt.mat[, grep("WT", colnames(rt.mat))]
dbp2.rt.mat <- rt.mat[, grep("dbp2", colnames(rt.mat))]

### Plotting
wt.rt.mean <- apply(wt.rt.mat, 2, mean, na.rm = T)
wt.rt.mean.norm <- (wt.rt.mean/sum(wt.rt.mean))*100
dbp2.rt.mean <- apply(dbp2.rt.mat, 2, mean, na.rm = T)
dbp2.rt.mean.norm <- (dbp2.rt.mean/sum(dbp2.rt.mean))*100
tiff("121718_PolII_readthrough_mRNA3UTR_withColors.tiff", height = 4, width = 5, units = 'in', res = 600)
par(lwd = 1, ps = 12, cex = 1)
plot(x = x, y = wt.rt.mean.norm, type = "l", col = "skyblue", lwd = 2, xlim = c(-458,470),
     ylim = c(0.6, 1.2), xaxt = "n", yaxt = "n",
     xlab = "distance to readthrough mRNA 3' ends", ylab = "Normalized Pol II occupancy",
     main = "RNA Pol II profile around the end of readthrough mRNAs")
polygon(c(-480, x, 510), c(0,wt.rt.mean.norm,0), col = adjustcolor("skyblue", alpha.f = 0.4), border = NA)
lines(x = x, y = dbp2.rt.mean.norm, type = "l", lwd = 2, col = "pink")
polygon(c(-480, x, 510), c(0,dbp2.rt.mean.norm,0), col = adjustcolor("pink", alpha.f = 0.4),
        border = NA)
axis(1, at = seq(-600, 600, by = 200), las = 1, lwd = 2)
axis(2, at = seq(0, 2, by = 0.2), las = 2, lwd = 2)
abline(v=0, lty=2, lwd = 1)
legend("topright", legend = c("1", "2"), cex=1.2, text.font=8,
       y.intersp = 1.5, box.lty=0, bg = "transparent", fill = c("skyblue", "pink"))
dev.off()
