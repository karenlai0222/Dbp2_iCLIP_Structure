
library(Rsamtools)
library(GenomicRanges)
library(GenomicAlignments)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
sc.bs3 <- BSgenome.Scerevisiae.UCSC.sacCer3
yeast.length <- seqlengths(sc.bs3)

#----------# 
##### Read the resulted matix (use the region center as a reference point)
mRNA3UTRend.mat <- read.delim("~/Desktop/Research_2018to19/Code_outputs/121718_PolII_mRNAs_3UTRend_matForPlot.txt",
                              comment.char = "#", as.is = T)
# rownames(sno.term.mat) <- sno.term.bed$name
wt.mRNAends.mat <- mRNA3UTRend.mat[, grep("WT", colnames(mRNA3UTRend.mat))]
dbp2.mRNAends.mat <- mRNA3UTRend.mat[, grep("dbp2", colnames(mRNA3UTRend.mat))]
dbp2.mRNAends.scaled <- dbp2.mRNAends.mat*(max(wt.mRNAends.mat, na.rm = T)/max(dbp2.mRNAends.mat, na.rm = T))

max(wt.mRNAends.mat, na.rm = T) ### 17.74
max(dbp2.mRNAends.mat, na.rm = T) ### 19.61
### Plotting
wt.mRNAends.mean <- apply(wt.mRNAends.mat, 2, mean, na.rm = T)
dbp2.mRNAends.scaled.mean <- apply(dbp2.mRNAends.scaled, 2, mean, na.rm = T)
x = seq(-490, 500, by = 10)
setwd("~/Desktop/Research_2018to19/Code_outputs/")
tiff("121718_PolII_nonOverlap_mRNA3UTR.tiff", height = 4, width = 5, units = 'in', res = 600)
par(lwd = 1, ps = 12, cex = 1)
plot(x = x, y = wt.mRNAends.mean, type = "l", lty = "longdash", lwd = 2, xlim = c(-458,470),
     ylim = c(0.9, 1.4), xaxt = "n", yaxt = "n",
     xlab = "distance to mRNA 3' UTR ends (nt)", ylab = "Normalized Pol II occupancy",
     main = "RNA Pol II profile around the end of mRNA 3' UTRs")
# polygon(c(-480, x, 510), c(0,wt.sno.term.mean,0), col = adjustcolor("skyblue", alpha.f = 0.4), border = NA)
lines(x = x, y = dbp2.mRNAends.scaled.mean, type = "l", lty = "solid", lwd = 2)
# polygon(c(-480, x, 510), c(0,dbp2.sno.term.scaled.mean,0), col = adjustcolor("pink", alpha.f = 0.4), border = NA)
axis(1, at = seq(-600, 600, by = 200), las = 1, lwd = 2)
axis(2, at = seq(0, 1.4, by = 0.2), las = 2, lwd = 2)
abline(v=0, lty=2, lwd = 1.5)
legend("topright", legend = c("1", "2"), lty = c("longdash", "solid"), lwd = 2, cex=1.2,
       text.font=8, y.intersp = 1.5, box.lty=0, bg = "transparent")
dev.off()
### With colors
setwd("~/Desktop/Research_2018to19/Code_outputs/")
tiff("121718_PolII_nonOverlap_mRNA3UTR_withColors.tiff", height = 4, width = 5, units = 'in', res = 600)
par(lwd = 1, ps = 12, cex = 1)
plot(x = x, y = wt.mRNAends.mean, type = "l", col = "skyblue", lwd = 2, xlim = c(-458,470),
     ylim = c(0.9, 1.4), xaxt = "n", yaxt = "n",
     xlab = "distance to mRNA 3' UTR ends (nt)", ylab = "Normalized Pol II occupancy",
     main = "RNA Pol II profile around the end of mRNA 3' UTRs")
polygon(c(-480, x, 510), c(0,wt.mRNAends.mean,0), col = adjustcolor("skyblue", alpha.f = 0.4),
        border = NA)
lines(x = x, y = dbp2.mRNAends.scaled.mean, type = "l", col = "pink", lwd = 2)
polygon(c(-480, x, 510), c(0,dbp2.mRNAends.scaled.mean,0), col = adjustcolor("pink", alpha.f = 0.4),
        border = NA)
axis(1, at = seq(-600, 600, by = 200), las = 1, lwd = 2)
axis(2, at = seq(0, 1.4, by = 0.2), las = 2, lwd = 2)
abline(v=0, lty=2, lwd = 1)
legend("topright", legend = c("1", "2"), cex=1.2, text.font=8,
       y.intersp = 1.5, box.lty=0, bg = "transparent", fill = c("skyblue", "pink"))
dev.off()

##### Read the resulted matix of read-through mRNAs only
rt.mat <- read.delim("~/Desktop/Research_2018to19/Code_outputs/121718_PolII_readthrough_mRNAs_3UTRend_matForPlot.txt",
                     comment.char = "#", as.is = T)
# rownames(sno.termS.mat) <- sno.term.bed$name
wt.rt.mat <- rt.mat[, grep("WT", colnames(rt.mat))]
dbp2.rt.mat <- rt.mat[, grep("dbp2", colnames(rt.mat))]
dbp2.rt.scaled <- dbp2.rt.mat*(max(wt.rt.mat, na.rm = T)/max(dbp2.rt.mat, na.rm = T))

# max(wt.sno.term.mat) ### 19.5
# max(dbp2.sno.term.mat) ### 15.97
### Plotting
wt.rt.mean <- apply(wt.rt.mat, 2, mean, na.rm = T)
dbp2.rt.mean <- apply(dbp2.rt.mat, 2, mean, na.rm = T)
dbp2.rt.scaled.mean <- apply(dbp2.rt.scaled, 2, mean, na.rm = T)
x = seq(-490, 500, by = 10)
setwd("~/Desktop/Research_2018to19/Code_outputs/")
tiff("121718_PolII_readthrough_mRNA3UTR.tiff", height = 4, width = 5, units = 'in', res = 600)
par(lwd = 1, ps = 12, cex = 1)
plot(x = x, y = wt.rt.mean, type = "l", lty = "longdash", lwd = 2, xlim = c(-458,470),
     ylim = c(1, 1.75), xaxt = "n", yaxt = "n",
     xlab = "distance to readthrough mRNA 3' ends", ylab = "Normalized Pol II occupancy",
     main = "RNA Pol II profile around the end of readthrough mRNAs")
# polygon(c(-480, x, 510), c(0,wt.sno.termS.mean,0), col = adjustcolor("skyblue", alpha.f = 0.4), border = NA)
# lines(x = x, y = dbp2.rt.mean,  type = "l", lwd = 2)
lines(x = x, y = dbp2.rt.scaled.mean, type = "l", lwd = 2)
# polygon(c(-480, x, 510), c(0,dbp2.sno.termS.scaled.mean,0), col = adjustcolor("pink", alpha.f = 0.4), border = NA)
axis(1, at = seq(-600, 600, by = 200), las = 1, lwd = 2)
axis(2, at = seq(0, 2, by = 0.2), las = 2, lwd = 2)
abline(v=0, lty=2, lwd = 1)
legend("topright", legend = c("1", "2"), lty = c("longdash", "solid"), lwd = 2, cex=1.2,
       text.font=8, y.intersp = 1.6, box.lty=0, bg = "transparent")
dev.off()
### With colors
setwd("~/Desktop/Research_2018to19/Code_outputs/")
tiff("121718_PolII_readthrough_mRNA3UTR_withColors.tiff", height = 4, width = 5, units = 'in', res = 600)
par(lwd = 1, ps = 12, cex = 1)
plot(x = x, y = wt.rt.mean, type = "l", col = "skyblue", lwd = 2, xlim = c(-458,470),
     ylim = c(1, 1.75), xaxt = "n", yaxt = "n",
     xlab = "distance to readthrough mRNA 3' ends", ylab = "Normalized Pol II occupancy",
     main = "RNA Pol II profile around the end of readthrough mRNAs")
polygon(c(-480, x, 510), c(0,wt.rt.mean,0), col = adjustcolor("skyblue", alpha.f = 0.4), border = NA)
# lines(x = x, y = dbp2.rt.mean,  type = "l", lwd = 2)
lines(x = x, y = dbp2.rt.scaled.mean, type = "l", lwd = 2, col = "pink")
polygon(c(-480, x, 510), c(0,dbp2.rt.scaled.mean,0), col = adjustcolor("pink", alpha.f = 0.4),
        border = NA)
axis(1, at = seq(-600, 600, by = 200), las = 1, lwd = 2)
axis(2, at = seq(0, 2, by = 0.2), las = 2, lwd = 2)
abline(v=0, lty=2, lwd = 1)
legend("topright", legend = c("1", "2"), cex=1.2, text.font=8,
       y.intersp = 1.5, box.lty=0, bg = "transparent", fill = c("skyblue", "pink"))
dev.off()
