
##########
##### This is to normalize the counts in each replicate for the metagene analysis to get a figure out of the three replicates
##########
##### After normalizing to the library size, the binding counts are also normalized to the RPKM of each transcript

##### Read the count tables in
d1.counts <- read.delim("060718_iCLIP_Dbp21_FilteredSites_inBoundmRNAs_countMat100bins.txt",
                        as.is = T)
d2.counts <- read.delim("060718_iCLIP_Dbp22_FilteredSites_inBoundmRNAs_countMat100bins.txt",
                        as.is = T)
d3.counts <- read.delim("060718_iCLIP_Dbp23_FilteredSites_inBoundmRNAs_countMat100bins.txt",
                        as.is = T)

##### Find out the "library size": the total count in each sample
sum(d1.counts) ### 427802
sum(d2.counts) ### 407930
sum(d3.counts) ### 455859

##### Divide the counts by "library size", per million
d1.counts.norm <- d1.counts/0.427802
d2.counts.norm <- d2.counts/0.407930
d3.counts.norm <- d3.counts/0.455859

##### Read the RPKM information in
rpkm.info <- read.delim("052617_S288C_gene_RPKM_from_structure-seq_control_no1153.txt",
                        as.is = T)

##### Normalize the count in each bin to the expression level of the transcript
for (i in 1:length(rownames(d1.counts.norm))) {
    d1.counts.norm[i,] <- d1.counts.norm[i,]/(rpkm.info$WT_RPKM_mean[which(rownames(rpkm.info)%in%rownames(d1.counts.norm)[[i]])])
}

for (i in 1:length(rownames(d2.counts.norm))) {
    d2.counts.norm[i,] <- d2.counts.norm[i,]/(rpkm.info$WT_RPKM_mean[which(rownames(rpkm.info)%in%rownames(d2.counts.norm)[[i]])])
}

for (i in 1:length(rownames(d3.counts.norm))) {
    d3.counts.norm[i,] <- d3.counts.norm[i,]/(rpkm.info$WT_RPKM_mean[which(rownames(rpkm.info)%in%rownames(d3.counts.norm)[[i]])])
}

##### Calculate the probability
### Convert tables into array
d1.counts.norm <- as.matrix(d1.counts.norm)
d2.counts.norm <- as.matrix(d2.counts.norm)
d3.counts.norm <- as.matrix(d3.counts.norm)
### prop.table
d1.counts.norm.prop <- prop.table(d1.counts.norm, 1)
d2.counts.norm.prop <- prop.table(d2.counts.norm, 1)
d3.counts.norm.prop <- prop.table(d3.counts.norm, 1)

##### Take average of the three replicates
avg.prop.table <- (d1.counts.norm.prop + d2.counts.norm.prop + d3.counts.norm.prop)/3

##### Plot Dbp2 occupancy over the whole transcript
tiff("060718_Dbp2iCLIP_FilteredSites_meta_inBoundmRNA_prop100bins3rep_norm_toRPKM_res600.tiff",
     height = 4, width = 5, units = 'in', res = 600)
par(ps = 12, cex = 1, lwd = 2)
plot(0.5:99.5, apply(avg.prop.table, 2, mean, na.rm=T), cex.lab = 1.2, cex.axis = 1.1, type="l",
     lwd=2, col="blue3", yaxt="n", ylim = c(0, 0.025), xaxt = "n", xlab = "Relative transcript positions",
     ylab = "Dbp2 occupancy", main = "3 replicates averaged, Dbp2-binding mRNAs, normalized to RPKM")
axis(1, at = seq(0, 100, by = 1), las = 2, lwd = 2)
axis(2, at = seq(0, 0.08, by = 0.01), lwd = 2)
abline(v=10, lty=2, lwd = 1.5)
abline(v=90, lty=2, lwd = 1.5)
dev.off()
### Normalize the value to 100%
avg.prop.table.avg <- apply(avg.prop.table, 2, mean, na.rm=T)
range(avg.prop.table.avg) ### 0.003584631 to 0.022862476
sum(avg.prop.table.avg)
tiff("060718_Dbp2iCLIP_FilteredSites_meta_inBoundmRNA_prop100bins3rep_norm_toRPKM_percent_res600.tiff",
     height = 4, width = 5, units = 'in', res = 600)
par(ps = 12, cex = 1, lwd = 2)
plot(0.5:99.5, apply(avg.prop.table, 2, mean, na.rm=T)*100, cex.lab = 1.2, cex.axis = 1.1, type="l",
     lwd=2, col="blue3", yaxt="n", ylim = c(0, 2.5), xaxt = "n", xlab = "Relative transcript positions",
     ylab = "Dbp2 occupancy", main = "3 replicates averaged, Dbp2-binding mRNAs, normalized to RPKM")
axis(1, at = seq(0, 100, by = 1), las = 2, lwd = 2)
axis(2, at = seq(0, 2.5, by = 0.5), las = 2, lwd = 2)
abline(v=10, lty=2, lwd = 1.5)
abline(v=90, lty=2, lwd = 1.5)
dev.off()

