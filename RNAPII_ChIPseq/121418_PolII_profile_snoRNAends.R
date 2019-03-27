
library(Rsamtools)
library(GenomicRanges)
library(GenomicAlignments)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
sc.bs3 <- BSgenome.Scerevisiae.UCSC.sacCer3
yeast.length <- seqlengths(sc.bs3)

##### Read annotation and data in
snoRNAs <- read.delim("~/Desktop/Research_projects_2015-2016/Data_general_use/snoRNA.txt", as.is = T)
plot.mat <- read.delim("~/Desktop/Research_2018to19/Code_outputs/121418_PolII_snoRNAends_matForPlot.txt",
                       as.is = T, comment.char = "#")
rownames(plot.mat) <- snoRNAs$id
wt.plot.mat <- plot.mat[, grep("WT", colnames(plot.mat))]
dbp2.plot.mat <- plot.mat[, grep("dbp2", colnames(plot.mat))]
dbp2.plot.mat.scaled <- dbp2.plot.mat*(max(wt.plot.mat)/max(dbp2.plot.mat))

max(wt.plot.mat) ### 19.43
max(dbp2.plot.mat) ### 15.95

wt.snoRNA.mean <- apply(wt.plot.mat, 2, mean)
dbp2.snoRNA.scaled.mean <- apply(dbp2.plot.mat.scaled, 2, mean)
x = seq(-490, 500, by = 10)
plot(x = x, y = wt.snoRNA.mean, type = "l", col = "skyblue", lwd = 2, xlim = c(-460,470),
     ylim = c(3, 9), xaxt = "n", yaxt = "n",
     xlab = "distance to mature snoRNA ends (nt)", ylab = "Normalized Pol II occupancy",
     main = "RNA Pol II profile around the end of mature snoRNAs")
polygon(c(-480, x, 510), c(0,wt.snoRNA.mean,0), col = adjustcolor("skyblue", alpha.f = 0.4), border = NA)
lines(x = x, y = dbp2.snoRNA.scaled.mean, type = "l", col = "pink", lwd = 2)
polygon(c(-480, x, 510), c(0,dbp2.snoRNA.scaled.mean,0), col = adjustcolor("pink", alpha.f = 0.4), border = NA)
axis(1, at = seq(-400, 400, by = 200), las = 1, lwd = 2)
axis(2, at = seq(0, 9, by = 1), las = 2, lwd = 2)
abline(v=0, lty=2, lwd = 1.5)
legend("topright", legend=c("WT", "dbp2∆"), col=c("skyblue", "pink"), cex=1.2,
       text.font=8, y.intersp = 1.5, box.lty=0, bg = "transparent", fill = c("skyblue", "pink"))

#----------# Use termination sites determined by Schaughency et al. PLoS Genet. 2014
##### Create a bed file for it
sno.term <- read.delim("~/Desktop/Research_projects_2015-2016/Data_general_use/snoRNA termination site from Schaughency 2014 PLOSGenet.txt",
                       as.is = T)
sno.term.gr <- GRanges(seqnames = Rle(sno.term$chromosome),
                       ranges = IRanges(start = sno.term$termination_site-1, end = sno.term$termination_site),
                       strand = Rle(sno.term$strand), id = sno.term$id, seqlengths = yeast.length)
sno.term.gr.broad <- resize(sno.term.gr, width = 50, fix = "center")
sno.term.bed <- data.frame(chrom = seqnames(sno.term.gr.broad), start = start(sno.term.gr.broad),
                           end = end(sno.term.gr.broad), name = mcols(sno.term.gr.broad)$id,
                           score = 0, strand = strand(sno.term.gr.broad))
### Output this file
write.table(sno.term.bed, file = "~/Desktop/Research_2018to19/Code_outputs/121418_49snoRNAs_termination_region_50nt.bed",
            quote = F, sep = "\t", row.names = F, col.names = F)

##### Read the resulted matix (use the region center as a reference point)
sno.term.mat <- read.delim("~/Desktop/Research_2018to19/Code_outputs/121418_PolII_49snoRNAs_termRegion_matForPlot.txt",
                           comment.char = "#", as.is = T)
rownames(sno.term.mat) <- sno.term.bed$name
wt.sno.term.mat <- sno.term.mat[, grep("WT", colnames(sno.term.mat))]
dbp2.sno.term.mat <- sno.term.mat[, grep("dbp2", colnames(sno.term.mat))]
dbp2.sno.term.scaled <- dbp2.sno.term.mat*(max(wt.sno.term.mat)/max(dbp2.sno.term.mat))

max(wt.sno.term.mat) ### 19.5
max(dbp2.sno.term.mat) ### 15.97
### Plotting
wt.sno.term.mean <- apply(wt.sno.term.mat, 2, mean)
dbp2.sno.term.scaled.mean <- apply(dbp2.sno.term.scaled, 2, mean)
x = seq(-490, 500, by = 10)
plot(x = x, y = wt.sno.term.mean, type = "l", col = "skyblue", lwd = 2, xlim = c(-458,470),
     ylim = c(1, 9), xaxt = "n", yaxt = "n",
     xlab = "distance to mature snoRNA ends (nt)", ylab = "Normalized Pol II occupancy",
     main = "RNA Pol II profile around the end of mature snoRNAs")
polygon(c(-480, x, 510), c(0,wt.sno.term.mean,0), col = adjustcolor("skyblue", alpha.f = 0.4), border = NA)
lines(x = x, y = dbp2.sno.term.scaled.mean, type = "l", col = "pink", lwd = 2)
polygon(c(-480, x, 510), c(0,dbp2.sno.term.scaled.mean,0), col = adjustcolor("pink", alpha.f = 0.4), border = NA)
axis(1, at = seq(-600, 600, by = 200), las = 1, lwd = 2)
axis(2, at = seq(0, 9, by = 2), las = 2, lwd = 2)
abline(v=0, lty=2, lwd = 1.5)
legend("topright", legend=c("WT", "dbp2∆"), col=c("skyblue", "pink"), cex=1.2,
       text.font=8, y.intersp = 1.5, box.lty=0, bg = "transparent", fill = c("skyblue", "pink"))

##### Read the resulted matix (use the region start as a reference point)
sno.termS.mat <- read.delim("~/Desktop/Research_2018to19/Code_outputs/121418_PolII_49snoRNAs_termStart_matForPlot.txt",
                            comment.char = "#", as.is = T)
rownames(sno.termS.mat) <- sno.term.bed$name
wt.sno.termS.mat <- sno.termS.mat[, grep("WT", colnames(sno.termS.mat))]
dbp2.sno.termS.mat <- sno.termS.mat[, grep("dbp2", colnames(sno.termS.mat))]
dbp2.sno.termS.scaled <- dbp2.sno.termS.mat*(max(wt.sno.termS.mat)/max(dbp2.sno.termS.mat))

# max(wt.sno.term.mat) ### 19.5
# max(dbp2.sno.term.mat) ### 15.97
### Plotting
wt.sno.termS.mean <- apply(wt.sno.termS.mat, 2, mean)
dbp2.sno.termS.scaled.mean <- apply(dbp2.sno.termS.scaled, 2, mean)
x = seq(-490, 500, by = 10)
setwd("~/Desktop/Research_2018to19/Code_outputs/")
tiff("121718_PolII_49snoRNAterm.tiff", height = 4, width = 5, units = 'in', res = 600)
par(lwd = 1, ps = 12, cex = 1)
plot(x = x, y = wt.sno.termS.mean, type = "l", col = "skyblue", lwd = 2, xlim = c(-458,470),
     ylim = c(1, 9), xaxt = "n", yaxt = "n",
     xlab = "distance to mature snoRNA ends (nt)", ylab = "Normalized Pol II occupancy",
     main = "RNA Pol II profile around the end of mature snoRNAs")
polygon(c(-480, x, 510), c(0,wt.sno.termS.mean,0), col = adjustcolor("skyblue", alpha.f = 0.4), border = NA)
lines(x = x, y = dbp2.sno.termS.scaled.mean, type = "l", col = "pink", lwd = 2)
polygon(c(-480, x, 510), c(0,dbp2.sno.termS.scaled.mean,0), col = adjustcolor("pink", alpha.f = 0.4), border = NA)
axis(1, at = seq(-600, 600, by = 200), las = 1, lwd = 2)
axis(2, at = seq(0, 10, by = 2), las = 2, lwd = 2)
abline(v=0, lty=2, lwd = 1.5)
legend("topright", legend=c("WT", "dbp2∆"), col=c("skyblue", "pink"), cex=1.2,
       text.font=8, y.intersp = 1.5, box.lty=0, bg = "transparent", fill = c("skyblue", "pink"))
dev.off()

setwd("~/Desktop/Research_2018to19/Code_outputs/")
tiff("121718_PolII_49snoRNAterm_noColor.tiff", height = 4, width = 5, units = 'in', res = 600)
par(lwd = 1, ps = 12, cex = 1)
plot(x = x, y = wt.sno.termS.mean, type = "l", lwd = 2, xlim = c(-458,470), lty = "twodash",
     ylim = c(1, 9), xaxt = "n", yaxt = "n",
     xlab = "distance to mature snoRNA ends (nt)", ylab = "Normalized Pol II occupancy",
     main = "RNA Pol II profile around the end of mature snoRNAs")
# polygon(c(-480, x, 510), c(0,wt.sno.termS.mean,0), col = adjustcolor("skyblue", alpha.f = 0.4), border = NA)
lines(x = x, y = dbp2.sno.termS.scaled.mean, type = "l", lwd = 2)
# polygon(c(-480, x, 510), c(0,dbp2.sno.termS.scaled.mean,0), col = adjustcolor("pink", alpha.f = 0.4), border = NA)
axis(1, at = seq(-600, 600, by = 200), las = 1, lwd = 2)
axis(2, at = seq(0, 10, by = 2), las = 2, lwd = 2)
abline(v=0, lty=2, lwd = 1.5)
legend("topright", legend=c("WT", "dbp2∆"), cex=1.2, lty = c("twodash", "solid"), lwd = 2,
       text.font=8, y.intersp = 1.5, box.lty=0, bg = "transparent")
dev.off()
