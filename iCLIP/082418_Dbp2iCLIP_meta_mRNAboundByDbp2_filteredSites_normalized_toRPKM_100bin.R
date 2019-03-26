
##########
##### This is to normalize the counts in each replicate for the metagene analysis to get a figure out of the three replicates
##########
##### After normalizing to the library size, the binding counts are also normalized to the RPKM of each transcript
##### In this transcript, I seperate the intron-containing genes out

##### Read the count tables in
d1.counts <- read.delim("~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/060618_Dbp2iCLIPdata_reprocessing/060718_iCLIP_Dbp21_FilteredSites_inBoundmRNAs_countMat100bins.txt",
                        as.is = T)
d2.counts <- read.delim("~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/060618_Dbp2iCLIPdata_reprocessing/060718_iCLIP_Dbp22_FilteredSites_inBoundmRNAs_countMat100bins.txt",
                        as.is = T)
d3.counts <- read.delim("~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/060618_Dbp2iCLIPdata_reprocessing/060718_iCLIP_Dbp23_FilteredSites_inBoundmRNAs_countMat100bins.txt",
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
rpkm.info <- read.delim("~/Desktop/Research_projects_2015-2016/DMS/Structure-seq_data/052417_ORF_DE_analysis_without_1153/052617_S288C_gene_RPKM_from_structure-seq_control_no1153.txt",
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

### Only keep those without introns
### Read in the intron-containing genes that are expressed based on the revised DE analysis
intron.gene.exp <- read.delim("~/Desktop/Research_projects_2015-2016/DMS/Structure-seq_data/071717_Nov_STAR_alignment_splicing_analysis_no1153/071717_all_expressed_intron_containing_genes_in_SC.txt",
                              header = F, quote = "", as.is = T)
d1.counts.norm.nosp <- d1.counts.norm[-which(rownames(d1.counts.norm)%in%intron.gene.exp$V1), ]
dim(d1.counts.norm.nosp) ### 1876 entries
d2.counts.norm.nosp <- d2.counts.norm[-which(rownames(d2.counts.norm)%in%intron.gene.exp$V1), ]
dim(d2.counts.norm.nosp) ### 1876 entries
d3.counts.norm.nosp <- d3.counts.norm[-which(rownames(d3.counts.norm)%in%intron.gene.exp$V1), ]
dim(d3.counts.norm.nosp) ### 1876 entries

##### Calculate the probability
### Convert tables into array
d1.counts.norm.nosp <- as.matrix(d1.counts.norm.nosp)
d2.counts.norm.nosp <- as.matrix(d2.counts.norm.nosp)
d3.counts.norm.nosp <- as.matrix(d3.counts.norm.nosp)
### prop.table
d1.counts.norm.prop <- prop.table(d1.counts.norm.nosp, 1)
d2.counts.norm.prop <- prop.table(d2.counts.norm.nosp, 1)
d3.counts.norm.prop <- prop.table(d3.counts.norm.nosp, 1)

##### Take average of the three replicates
avg.prop.table.nosp <- (d1.counts.norm.prop + d2.counts.norm.prop + d3.counts.norm.prop)/3

##### Plot Dbp2 occupancy over the whole transcript
setwd("~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/060618_Dbp2iCLIPdata_reprocessing/")
tiff("060718_Dbp2iCLIP_FilteredSites_meta_inBoundmRNA_prop100bins3rep_norm_toRPKM_res600.tiff",
     height = 4, width = 5, units = 'in', res = 600)
par(ps = 12, cex = 1, lwd = 2)
plot(0.5:99.5, apply(avg.prop.table.nosp, 2, mean, na.rm=T), cex.lab = 1.2, cex.axis = 1.1, type="l",
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
setwd("~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/060618_Dbp2iCLIPdata_reprocessing/")
tiff("060718_Dbp2iCLIP_FilteredSites_meta_inBoundmRNA_prop100bins3rep_norm_toRPKM_percent_res600.tiff",
     height = 4, width = 5, units = 'in', res = 600)
par(ps = 12, cex = 1, lwd = 2)
plot(0.5:99.5, apply(avg.prop.table, 2, mean, na.rm=T)*100, cex.lab = 1.2, cex.axis = 1.1, type="l",
     lwd=2, col="blue3", yaxt="n", ylim = c(0, 2.5), xaxt = "n", xlab = "Relative transcript positions",
     ylab = "Dbp2 occupancy", main = "3 replicates averaged, Dbp2-binding mRNAs, normalized to RPKM")
axis(1, at = seq(0, 100, by = 1), las = 2, lwd = 2)
axis(2, at = seq(0, 0.8, by = 0.5), las = 2, lwd = 2)
abline(v=10, lty=2, lwd = 1.5)
abline(v=90, lty=2, lwd = 1.5)
dev.off()

####################################################################################################
#----------# The following is to compare the binding pattern in readthrough v.s. non-readthrough transcripts
##### Do the analyis for readthrough genes in dbp2-null (also bound)
### import the gene list with a readthrough defect in dbp2-null (new list produced without 1153 replicate)
rt1 <- read.delim(file = "~/Desktop/Research_projects_2015-2016/DMS/Structure-seq_data/101717_Nov_STAR_alignment_readthrough_analysis_YN0809data_no1153_nonOverlap_v3/101717_all_readthrough_mRNAs_v3.txt", header = F, quote = "", as.is = T)
### get the probability table for genes that are bound and readthrough in dbp2-null
avg.prop.table.nosp.rt1 <- avg.prop.table.nosp[which(rownames(avg.prop.table.nosp)%in%rt1$V1),]
### plot Dbp2 occupancy over these transcripts
plot(0.5:99.5,apply(avg.prop.table.nosp.rt1, 2, mean, na.rm=T), cex.lab = 1.2, cex.axis = 1.1,
     type="l", lwd=4, col="blue3", yaxt="n", ylim = c(0, 0.025), xaxt = "n",
     xlab = "Relative transcript positions", ylab = "Dbp2 occupancy",
     main = "3 replicates normalized, Dbp2-binding and readthrough mRNAs, all sites")
axis(1, at = seq(0, 100, by = 1), las = 2, lwd = 2)
axis(2, at = seq(0, 0.08, by = 0.01), lwd = 2)
abline(v=10, lty=2, lwd = 1.5)
abline(v=90, lty=2, lwd = 1.5)

##### Do the analyis for bound mRNAs but no readthrough in dbp2-null
### get the probability table for genes that are bound but not readthrough in dbp2-null
avg.prop.table.nosp.nort1 <- avg.prop.table.nosp[-which(rownames(avg.prop.table.nosp)%in%rt1$V1),]
### plot Dbp2 occupancy over these transcripts
plot(0.5:99.5,apply(avg.prop.table.nosp.nort1, 2, mean, na.rm=T), cex.lab = 1.2, cex.axis = 1.1, type="l",
     lwd=4, col="blue3", yaxt="n", ylim = c(0, 0.025), xaxt = "n", xlab = "Relative transcript positions",
     ylab = "Dbp2 occupancy", main = "3 replicates normalized, Dbp2-binding without readthrough mRNAs, all sites")
axis(1, at = seq(1, 68, by = 1), las = 2, lwd = 2)
axis(2, at = seq(0, 0.08, by = 0.01), lwd = 2)
abline(v=6, lty=2, lwd = 3)
abline(v=62, lty=2, lwd = 3)

# ##### Try to overlay the two plots
# setwd("~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/101817_Dbp2_refiltered_binding_vs_revised_readthrough_no1153/")
# tiff("101817_Dbp2_iCLIP_all_filtered_sites_rt_vs_NOrt_metagene_3rep_norm_toRPKM_res600.tiff", height = 4, width = 5, units = 'in', res = 600)
# par(lwd = 2)
# plot(0.5:65.5, apply(avg.prop.table.nort1, 2, mean, na.rm=T), cex.lab = 1.2, cex.axis = 1.1, type="l", lwd=2, col="blue3", yaxt="n", ylim = c(0, 0.04), xaxt = "n", xlab = "Relative transcript positions", ylab = "Dbp2 occupancy", main = "3 replicates averaged, Dbp2-binding mRNAs, rt v.s. NOrt , normalized to RPKM")
# lines(0.5:65.5, apply(avg.prop.table.rt1, 2, mean, na.rm=T), cex.lab = 1.2, cex.axis = 1.1, lty = 1, lwd=2, col="firebrick", yaxt="n", ylim = c(0, 0.04), xaxt = "n", xlab = "Relative transcript positions", ylab = "Dbp2 occupancy")
# axis(1, at = seq(0, 66, by = 1), las = 2, lwd = 2)
# axis(2, at = seq(0, 0.1, by = 0.01), lwd = 2)
# abline(v=5, lty=2, lwd = 2)
# abline(v=60, lty=2, lwd = 2)
# legend("top", legend=c("Readthrough", "No readthrough"), col=c("firebrick", "blue3"), lty=c(1,1), lwd = 2, cex=1.2, text.font=8, y.intersp = 1.5, box.lty=0, bg = "transparent")
# dev.off()

#----------# The following is to compare the binding pattern in up, down, not changed mRNAs
##### keep all bound mRNAs
### Read in the intron-containing genes that are expressed based on the revised DE analysis
dim(d1.counts.norm) ### 2045 entries
dim(d2.counts.norm) ### 2045 entries
dim(d3.counts.norm) ### 2045 entries

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
avg.prop.table<- (d1.counts.norm.prop + d2.counts.norm.prop + d3.counts.norm.prop)/3

##### Group them based on the expression level change
edgr.sum <- read.delim("~/Desktop/Research_projects_2015-2016/DMS/Structure-seq_data/052417_ORF_DE_analysis_without_1153/052617_Nov_STAR_aligned_ORF_no1153_edgeR_summary.txt",
                       as.is = T)
up <- rownames(edgr.sum)[which(edgr.sum$change == 1)]
down <- rownames(edgr.sum)[which(edgr.sum$change == -1)]
no.change <- rownames(edgr.sum)[which(edgr.sum$change == 0)]

avg.prop.table.up <- avg.prop.table[which(rownames(avg.prop.table)%in%up),]
avg.prop.table.down <- avg.prop.table[which(rownames(avg.prop.table)%in%down),]
avg.prop.table.no <- avg.prop.table[which(rownames(avg.prop.table)%in%no.change),]
###
write.table(rownames(avg.prop.table.up),
            file = "~/Desktop/Research_2018to19/Code_outputs/021419_Dbp2bound_up_mRNAsList.txt",
            quote = F, sep = "\t", row.names = F, col.names = F)

##### Plot
plot(0.5:99.5,apply(avg.prop.table.no, 2, mean, na.rm=T), cex.lab = 1.2, cex.axis = 1.1, type="l",
     lwd=4, col="black", yaxt="n", ylim = c(0, 0.028), xaxt = "n", xlab = "Relative transcript positions",
     ylab = "Dbp2 occupancy", main = "Dbp2-binding, no change ")
axis(1, at = seq(0, 100, by = 10), las = 1, lwd = 2)
axis(2, at = seq(0, 0.08, by = 0.005), lwd = 2, las = 1)
abline(v=10, lty=2, lwd = 1.5)
abline(v=90, lty=2, lwd = 1.5)
plot(0.5:99.5,apply(avg.prop.table.up, 2, mean, na.rm=T), cex.lab = 1.2, cex.axis = 1.1, type="l",
     lwd=4, col="black", yaxt="n", ylim = c(0, 0.028), xaxt = "n", xlab = "Relative transcript positions",
     ylab = "Dbp2 occupancy", main = "Dbp2-binding, up")
axis(1, at = seq(0, 100, by = 10), las = 1, lwd = 2)
axis(2, at = seq(0, 0.08, by = 0.005), lwd = 2, las = 1)
abline(v=10, lty=2, lwd = 1.5)
abline(v=90, lty=2, lwd = 1.5)
plot(0.5:99.5,apply(avg.prop.table.down, 2, mean, na.rm=T), cex.lab = 1.2, cex.axis = 1.1, type="l",
     lwd=4, col="black", yaxt="n", ylim = c(0, 0.028), xaxt = "n", xlab = "Relative transcript positions",
     ylab = "Dbp2 occupancy", main = "Dbp2-binding, down")
axis(1, at = seq(0, 100, by = 10), las = 1, lwd = 2)
axis(2, at = seq(0, 0.08, by = 0.005), lwd = 2, las = 1)
abline(v=10, lty=2, lwd = 1.5)
abline(v=90, lty=2, lwd = 1.5)

#####
setwd("~/Desktop/Research_2018to19/Code_outputs/")
tiff("022219_Dbp2_meta_in_bound_mRNA_3levels_100bins_res600.tiff", height = 4, width = 5,
     units = 'in', res = 600)
par(lwd = 2)
data <- data.frame(x = c(0.5:99.5), y = apply(avg.prop.table.no, 2, mean, na.rm=T))
loess_fit <- loess(y ~ x, data, span = 0.2)
plot(data$x, predict(loess_fit), col = "black", lwd = 2, xaxt = "n", yaxt = "n", ylab = "", xlab = "",
      type = "l", ylim = c(0, 0.027))
data <- data.frame(x = c(0.5:99.5), y = apply(avg.prop.table.down, 2, mean, na.rm=T))
loess_fit <- loess(y ~ x, data, span = 0.2)
lines(data$x, predict(loess_fit), col = "blue1", lwd = 2, xaxt = "n", yaxt = "n", ylab = "", xlab = "",
     type = "l", ylim = c(0, 0.027))
data <- data.frame(x = c(0.5:99.5), y = apply(avg.prop.table.up, 2, mean, na.rm=T))
loess_fit <- loess(y ~ x, data, span = 0.2)
lines(data$x, predict(loess_fit), col = "red", lwd = 2, xaxt = "n", yaxt = "n", ylab = "", xlab = "",
     type = "l", ylim = c(0, 0.027))
axis(1, at = seq(0, 100, by = 10), las = 1, lwd = 2)
axis(2, at = seq(0, 0.08, by = 0.005), lwd = 2, las = 1)
abline(v=10, lty=2, lwd = 1.5)
abline(v=90, lty=2, lwd = 1.5)
legend("top", legend=c("No change", "Up", "Down"), col=c("black", "red", "blue1"), lty=c(1,1),
       lwd = 2, cex=1, text.font=8, y.intersp = 1.5, box.lty=0, bg = "transparent")
dev.off()

