
library(Rsamtools)
library(GenomicRanges)
library(GenomicAlignments)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
sc.bs3 <- BSgenome.Scerevisiae.UCSC.sacCer3
yeast.length <- seqlengths(sc.bs3)

#----------# Use termination sites determined by Schaughency et al. PLoS Genet. 2014
##### Create a bed file for it
sno.term <- read.delim("snoRNA termination site from Schaughency 2014 PLOSGenet.txt",
                       as.is = T)
sno.term.gr <- GRanges(seqnames = Rle(sno.term$chromosome),
                       ranges = IRanges(start = sno.term$termination_site-1, end = sno.term$termination_site),
                       strand = Rle(sno.term$strand), id = sno.term$id, seqlengths = yeast.length)
sno.term.gr.broad <- resize(sno.term.gr, width = 50, fix = "center")
sno.term.bed <- data.frame(chrom = seqnames(sno.term.gr.broad), start = start(sno.term.gr.broad),
                           end = end(sno.term.gr.broad), name = mcols(sno.term.gr.broad)$id,
                           score = 0, strand = strand(sno.term.gr.broad))
### Output this file
write.table(sno.term.bed, file = "121418_49snoRNAs_termination_region_50nt.bed",
            quote = F, sep = "\t", row.names = F, col.names = F)


##### Read the resulted matix (use the region start as a reference point)
sno.termS.mat <- read.delim("121418_PolII_49snoRNAs_termStart_matForPlot.txt",
                            comment.char = "#", as.is = T)
rownames(sno.termS.mat) <- sno.term.bed$name
wt.sno.termS.mat <- sno.termS.mat[, grep("WT", colnames(sno.termS.mat))]
dbp2.sno.termS.mat <- sno.termS.mat[, grep("dbp2", colnames(sno.termS.mat))]

### Plotting
wt.sno.termS.mean <- apply(wt.sno.termS.mat, 2, mean)
wt.sno.termS.mean.norm <- (wt.sno.termS.mean/sum(wt.sno.termS.mean))*100
dbp2.sno.termS.mean <- apply(dbp2.sno.termS.mat, 2, mean)
dbp2.sno.termS.mean.norm <- (dbp2.sno.termS.mean/sum(dbp2.sno.termS.mean))*100
x = seq(-490, 500, by = 10)
tiff("121718_PolII_49snoRNAterm.tiff", height = 4, width = 5, units = 'in', res = 600)
par(lwd = 1, ps = 12, cex = 1)
plot(x = x, y = wt.sno.termS.mean.norm, type = "l", col = "skyblue", lwd = 2, xlim = c(-458,470),
     ylim = c(0.05, 1.6), xaxt = "n", yaxt = "n",
     xlab = "distance to mature snoRNA ends (nt)", ylab = "Normalized Pol II occupancy",
     main = "RNA Pol II profile around the end of mature snoRNAs")
polygon(c(-480, x, 510), c(0,wt.sno.termS.mean.norm,0), col = adjustcolor("skyblue", alpha.f = 0.4), border = NA)
lines(x = x, y = dbp2.sno.termS.mean.norm, type = "l", col = "pink", lwd = 2)
polygon(c(-480, x, 510), c(0,dbp2.sno.termS.mean.norm,0), col = adjustcolor("pink", alpha.f = 0.4), border = NA)
axis(1, at = seq(-600, 600, by = 200), las = 1, lwd = 2)
axis(2, at = seq(0, 2, by = 0.5), las = 2, lwd = 2)
abline(v=0, lty=2, lwd = 1.5)
legend("topright", legend=c("WT", "dbp2∆"), col=c("skyblue", "pink"), cex=1.2,
       text.font=8, y.intersp = 1.5, box.lty=0, bg = "transparent", fill = c("skyblue", "pink"))
dev.off()
