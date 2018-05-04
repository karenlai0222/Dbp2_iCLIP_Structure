
##########
##### This is to normalize the counts in each replicate for the metagene analysis to get a figure out of the three replicates
##########
##### After normalizing to the library size, the binding counts are also normalized to the RPKM of each transcript

##### Read the count tables in
d1.counts <- read.delim("~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/010418_Data_from_Nadia/010518_iCLIP_Dbp21_all_filtered_sites_in_bound_mRNAs_count_mat_100bins.txt", as.is = T)
d2.counts <- read.delim("~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/010418_Data_from_Nadia/010518_iCLIP_Dbp22_all_filtered_sites_in_bound_mRNAs_count_mat_100bins.txt", as.is = T)
d3.counts <- read.delim("~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/010418_Data_from_Nadia/010518_iCLIP_Dbp23_all_filtered_sites_in_bound_mRNAs_count_mat_100bins.txt", as.is = T)

##### Find out the "library size": the total count in each sample
sum(d1.counts) ### 310008
sum(d2.counts) ### 309953
sum(d3.counts) ### 326715

##### Divide the counts by "library size", per million
d1.counts.norm <- d1.counts/0.310008
d2.counts.norm <- d2.counts/0.309953
d3.counts.norm <- d3.counts/0.326715

##### Read the RPKM information in
rpkm.info <- read.delim("~/Desktop/Research_projects_2015-2016/DMS/Structure-seq_data/052417_ORF_DE_analysis_without_1153/052617_S288C_gene_RPKM_from_structure-seq_control_no1153.txt", as.is = T)

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
setwd("~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/010418_Data_from_Nadia/")
tiff("010518_Dbp2_iCLIP_all_sites_metagene_in_bound_mRNA_prop_100bins_3rep_norm_toRPKM_res600.tiff", height = 4, width = 5, units = 'in', res = 600)
par(ps = 12, cex = 1, lwd = 2)
plot(0.5:99.5, apply(avg.prop.table, 2, mean, na.rm=T), cex.lab = 1.2, cex.axis = 1.1, type="l", lwd=2, col="blue3", yaxt="n", ylim = c(0, 0.025), xaxt = "n", xlab = "Relative transcript positions", ylab = "Dbp2 occupancy", main = "3 replicates averaged, Dbp2-binding mRNAs, normalized to RPKM")
axis(1, at = seq(0, 100, by = 1), las = 2, lwd = 2)
axis(2, at = seq(0, 0.08, by = 0.01), lwd = 2)
abline(v=10, lty=2, lwd = 1.5)
abline(v=90, lty=2, lwd = 1.5)
dev.off()
### Normalize the value to 100%
avg.prop.table.avg <- apply(avg.prop.table, 2, mean, na.rm=T)
range(avg.prop.table.avg) ### 0.003256246 to 0.022916678
sum(avg.prop.table.avg)
setwd("~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/010418_Data_from_Nadia/")
tiff("010918_Dbp2_iCLIP_all_sites_metagene_in_bound_mRNA_prop_100bins_3rep_norm_toRPKM_percent_res600.tiff", height = 4, width = 5, units = 'in', res = 600)
par(ps = 12, cex = 1, lwd = 2)
plot(0.5:99.5, apply(avg.prop.table, 2, mean, na.rm=T)*100, cex.lab = 1.2, cex.axis = 1.1, type="l", lwd=2, col="blue3", yaxt="n", ylim = c(0, 2.5), xaxt = "n", xlab = "Relative transcript positions", ylab = "Dbp2 occupancy", main = "3 replicates averaged, Dbp2-binding mRNAs, normalized to RPKM")
axis(1, at = seq(0, 100, by = 1), las = 2, lwd = 2)
axis(2, at = seq(0, 2.5, by = 0.5), las = 2, lwd = 2)
abline(v=10, lty=2, lwd = 1.5)
abline(v=90, lty=2, lwd = 1.5)
dev.off()

####################################################################################################
##### The folowing parts have not been done for the new data
#----------# The following is to compare the binding pattern in readthrough v.s. non-readthrough transcripts
##### Do the analyis for readthrough genes in dbp2-null (also bound)
### import the gene list with a readthrough defect in dbp2-null (new list produced without 1153 replicate)
rt1 <- read.delim(file = "~/Desktop/Research_projects_2015-2016/DMS/Structure-seq_data/101717_Nov_STAR_alignment_readthrough_analysis_YN0809data_no1153_nonOverlap_v3/101717_all_readthrough_mRNAs_v3.txt", header = F, quote = "", as.is = T)
### get the probability table for genes that are bound and readthrough in dbp2-null
avg.prop.table.rt1 <- avg.prop.table[which(rownames(avg.prop.table)%in%rt1$V1),]
### plot Dbp2 occupancy over these transcripts
plot(0.5:65.5,apply(avg.prop.table.rt1, 2, mean, na.rm=T), cex.lab = 1.2, cex.axis = 1.1, type="l", lwd=4, col="blue3", yaxt="n", ylim = c(0, 0.06), xaxt = "n", xlab = "Relative transcript positions", ylab = "Dbp2 occupancy", main = "3 replicates normalized, Dbp2-binding and readthrough mRNAs, all sites")
axis(1, at = seq(1, 68, by = 1), las = 2, lwd = 2)
axis(2, at = seq(0, 0.08, by = 0.01), lwd = 2)
abline(v=6, lty=2, lwd = 3)
abline(v=62, lty=2, lwd = 3)

##### Do the analyis for bound mRNAs but no readthrough in dbp2-null
### get the probability table for genes that are bound but not readthrough in dbp2-null
avg.prop.table.nort1 <- avg.prop.table[-which(rownames(avg.prop.table)%in%rt1$V1),]
### plot Dbp2 occupancy over these transcripts
plot(0.5:65.5,apply(avg.prop.table.nort1, 2, mean, na.rm=T), cex.lab = 1.2, cex.axis = 1.1, type="l", lwd=4, col="blue3", yaxt="n", ylim = c(0, 0.06), xaxt = "n", xlab = "Relative transcript positions", ylab = "Dbp2 occupancy", main = "3 replicates normalized, Dbp2-binding without readthrough mRNAs, all sites")
axis(1, at = seq(1, 68, by = 1), las = 2, lwd = 2)
axis(2, at = seq(0, 0.08, by = 0.01), lwd = 2)
abline(v=6, lty=2, lwd = 3)
abline(v=62, lty=2, lwd = 3)

##### Try to overlay the two plots
setwd("~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/101817_Dbp2_refiltered_binding_vs_revised_readthrough_no1153/")
tiff("101817_Dbp2_iCLIP_all_filtered_sites_rt_vs_NOrt_metagene_3rep_norm_toRPKM_res600.tiff", height = 4, width = 5, units = 'in', res = 600)
par(lwd = 2)
plot(0.5:65.5, apply(avg.prop.table.nort1, 2, mean, na.rm=T), cex.lab = 1.2, cex.axis = 1.1, type="l", lwd=2, col="blue3", yaxt="n", ylim = c(0, 0.04), xaxt = "n", xlab = "Relative transcript positions", ylab = "Dbp2 occupancy", main = "3 replicates averaged, Dbp2-binding mRNAs, rt v.s. NOrt , normalized to RPKM")
lines(0.5:65.5, apply(avg.prop.table.rt1, 2, mean, na.rm=T), cex.lab = 1.2, cex.axis = 1.1, lty = 1, lwd=2, col="firebrick", yaxt="n", ylim = c(0, 0.04), xaxt = "n", xlab = "Relative transcript positions", ylab = "Dbp2 occupancy")
axis(1, at = seq(0, 66, by = 1), las = 2, lwd = 2)
axis(2, at = seq(0, 0.1, by = 0.01), lwd = 2)
abline(v=5, lty=2, lwd = 2)
abline(v=60, lty=2, lwd = 2)
legend("top", legend=c("Readthrough", "No readthrough"), col=c("firebrick", "blue3"), lty=c(1,1), lwd = 2, cex=1.2, text.font=8, y.intersp = 1.5, box.lty=0, bg = "transparent")
dev.off()

# setwd("~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/080717_Dbp2_iCLIP_metagene_analysis_after_bug_fixed/")
# tiff("082117_Dbp2_iCLIP_all_xlinking_sites_rt_vs_NOrt_noSet1targets_metagene_3rep_norm_toRPKM_res600.tiff", height = 4, width = 5, units = 'in', res = 600)
# par(lwd = 2)
# plot(apply(avg.prop.table.nort1, 2, mean, na.rm=T), cex.lab = 1.2, cex.axis = 1.1, type="l", lwd=2, col="blue3", yaxt="n", ylim = c(0, 0.04), xaxt = "n", xlab = "Relative transcript positions", ylab = "Dbp2 occupancy", main = "3 replicates averaged, Dbp2-binding mRNAs, rt v.s. NOrt , normalized to RPKM")
# lines(apply(avg.prop.table.rt1, 2, mean, na.rm=T), cex.lab = 1.2, cex.axis = 1.1, lty = 1, lwd=2, col="firebrick", yaxt="n", ylim = c(0, 0.04), xaxt = "n", xlab = "Relative transcript positions", ylab = "Dbp2 occupancy")
# axis(1, at = seq(1, 68, by = 1), las = 2, lwd = 2)
# axis(2, at = seq(0, 0.08, by = 0.01), lwd = 2)
# abline(v=6, lty=2, lwd = 2)
# abline(v=62, lty=2, lwd = 2)
# legend("top", legend=c("Readthrough", "No readthrough"), col=c("firebrick", "blue3"), lty=c(1,1), lwd = 2, cex=1.2, text.font=8, y.intersp = 1.5, box.lty=0, bg = "transparent")
# dev.off()

##### Find the difference between the two groups
avg.prop.table.rt1.mean <- apply(avg.prop.table.rt1, 2, mean, na.rm=T)
avg.prop.table.nort1.mean <- apply(avg.prop.table.nort1, 2, mean, na.rm=T)
dif <- avg.prop.table.rt1.mean - avg.prop.table.nort1.mean
plot(0.5:66.5, dif, cex.lab = 1.2, cex.axis = 1.1, type="l", lwd=2, col="blue3", yaxt="n", ylim = c(-0.012, 0.01), xaxt = "n", xlab = "Relative transcript positions", ylab = "Dbp2 occupancy", main = "3 replicates averaged, Dbp2-binding mRNAs, rt v.s. NOrt , normalized to RPKM")
abline(h=0, lty=2, lwd = 2)
axis(1, at = seq(0, 67, by = 1), las = 2, lwd = 2)
axis(2, at = seq(-0.02, 0.02, by = 0.002), lwd = 2)
abline(v=5, lty=2, lwd = 3)
abline(v=60, lty=2, lwd = 3)

##### Plot the trend line for the difference values
require(lattice)
x <- c(0.5:66.5)
y <- dif

data <- data.frame(x = x, y = y)

plot(0.5:66.5, dif, cex.lab = 1.2, cex.axis = 1.1, type = "l", lwd=2, col="black", yaxt="n", ylim = c(-0.012, 0.008), xaxt = "n", xlab = "Relative transcript positions", ylab = "Dbp2 occupancy", main = "Difference in Dbp2-binding, rt v.s. NOrt mRNAs")
lines(predict(loess(y~x, data = data)), lwd=3, col="red")
abline(h=0, lty=2, lwd = 2, col = "gray40")
axis(1, at = seq(0, 67, by = 1), las = 2, lwd = 2)
axis(2, at = seq(-0.02, 0.02, by = 0.002), lwd = 2)
abline(v=5, lty=2, lwd = 2, col = "gray40")
abline(v=61, lty=2, lwd = 2, col = "gray40")
### Output in high resolution
setwd("~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/080717_Dbp2_iCLIP_metagene_analysis_after_bug_fixed/")
tiff("090817_Dbp2_iCLIP_all_xlinking_sites_rt_vs_NOrt_noSet1targets_difference_res600.tiff", height = 4, width = 5, units = 'in', res = 600)
par(lwd = 2)
plot(0.5:66.5, dif, cex.lab = 1.2, cex.axis = 1.1, type = "l", lwd=2, col="black", yaxt="n", ylim = c(-0.012, 0.01), xaxt = "n", xlab = "Relative transcript positions", ylab = "Difference in Dbp2 occupancy", main = "Difference in Dbp2-binding, rt v.s. NOrt mRNAs")
lines(predict(loess(y~x, data = data, span = 0.5)), lwd=2, col="red")
abline(h=0, lty=2, lwd = 2, col = "gray40")
axis(1, at = seq(0, 67, by = 1), las = 2, lwd = 2)
axis(2, at = seq(-0.02, 0.02, by = 0.005), lwd = 2)
abline(v=5, lty=2, lwd = 2, col = "gray40")
abline(v=61, lty=2, lwd = 2, col = "gray40")
dev.off()


