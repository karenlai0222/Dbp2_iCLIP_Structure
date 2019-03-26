
##########
##### This is to assign Dbp2-binding sites (reads) to different RNA classes
##########
##### This is run on Snyder
setwd("/scratch/snyder/l/lai64/060618_Dbp2iCLIPdata_reprocessing/Dbp2iCLIP_sitesFromReproducedReads")

library(Rsamtools)
library(GenomicRanges)
library(GenomicAlignments)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
sc.bs3 <- BSgenome.Scerevisiae.UCSC.sacCer3
yeast.length <- seqlengths(sc.bs3)

##### Read Dbp2 sites in
d1.filt <- read.delim("060618_iCLIP_Dbp21_intersectedReads_allxlinkingSites.txt", as.is = T)
d2.filt <- read.delim("060618_iCLIP_Dbp22_intersectedReads_allxlinkingSites.txt", as.is = T)
d3.filt <- read.delim("060618_iCLIP_Dbp23_intersectedReads_allxlinkingSites.txt", as.is = T)

##### Read R64-2-1 in
r64 <- read.delim("../R64-2-1_yeast_id_name.txt", as.is = T)
### Only keep the types that we are interested in
r64.keep <- r64[-which(r64$type=="blocked_reading_frame"|r64$type=="pseudogene"|r64$type=="transposable_element_gene"),]

# #----------# First try to incorporate UTR regions and lncRNAs into the annotation to count the "reads"
##### Read UTR information in
yn.utr.s <- read.delim("../062216_Yassour_Nagalakshmi_outermost_UTR.txt", as.is = T)
yn.utr.s <- yn.utr.s[order(yn.utr.s$chromosome, yn.utr.s$transcript.start),]

##### Put the transcript start and end as the start and end of genes
r64.mod <- r64.keep
for (i in 1:length(r64.mod$id)) {
    if (r64.mod$id[i]%in%yn.utr.s$id){
        index <- which(yn.utr.s$id%in%r64.mod$id[i])
        r64.mod$start[i] <- yn.utr.s$transcript.start[index]
        r64.mod$end[i] <- yn.utr.s$transcript.end[index]
    }
}

##### also import lncRNA information
lncrna <- read.delim("../010417_all_lncRNA_merged.txt", as.is = T)

##### Make GRanges of the annotations
r64.gr <- GRanges(seqnames = Rle(r64.mod$chromosome),
                  ranges = IRanges(start = r64.mod$start, end = r64.mod$end),
                  strand = Rle(r64.mod$strand), id = r64.mod$id, type = r64.mod$type,
                  seqlengths = yeast.length)
r64.gr <- trim(r64.gr)
lncrna.gr <- GRanges(seqnames = Rle(lncrna$chromosome),
                     ranges = IRanges(start = lncrna$start, end = lncrna$end),
                     strand = Rle(lncrna$strand), id = lncrna$ID, type = "lncRNA",
                     seqlengths = yeast.length)
all.gr <- c(r64.gr, lncrna.gr)

# ##### Make GRanges of "reads"
# d1.reads.gr <- GRanges(seqnames = Rle(d1.filt$chr), ranges = IRanges(start = d1.filt$read_start, end = d1.filt$read_end), strand = Rle(d1.filt$strand), id = d1.filt$id, seqlengths = yeast.length)
# d2.reads.gr <- GRanges(seqnames = Rle(d2.filt$chr), ranges = IRanges(start = d2.filt$read_start, end = d2.filt$read_end), strand = Rle(d2.filt$strand), id = d2.filt$id, seqlengths = yeast.length)
# d3.reads.gr <- GRanges(seqnames = Rle(d3.filt$chr), ranges = IRanges(start = d3.filt$read_start, end = d3.filt$read_end), strand = Rle(d3.filt$strand), id = d3.filt$id, seqlengths = yeast.length)

# ##### Use summarizeOverlap to count the "reads"
# ### Generate an empty matrix first
# reads.ct <- data.frame(row.names = mcols(all.gr)$id, Dbp21 = NA, Dbp22 = NA, Dbp23 = NA, type = mcols(all.gr)$type)
# # reads.ct <- matrix(nrow = length(all.gr), ncol = 4)
# # rownames(reads.ct) <- mcols(all.gr)$id
# # colnames(reads.ct) <- c("Dbp2-1", "Dbp2-2", "Dbp2-3", "type")
# # reads.ct[,4] <- mcols(all.gr)$type
# ### start counting
# reads.ct[,1] <- assays(summarizeOverlaps(all.gr, d1.reads.gr, ignore.strand = F, mode = "IntersectionNotEmpty"))$counts
# reads.ct[,2] <- assays(summarizeOverlaps(all.gr, d2.reads.gr, ignore.strand = F, mode = "IntersectionNotEmpty"))$counts
# reads.ct[,3] <- assays(summarizeOverlaps(all.gr, d3.reads.gr, ignore.strand = F, mode = "IntersectionNotEmpty"))$counts
# colnames(reads.ct) <- c("Dbp2-1", "Dbp2-2", "Dbp2-3", "type")
# 
# colSums(reads.ct[,1:3])
# ###  Dbp21  Dbp22  Dbp23 
# ### 385643 441733 452372 
# 
# ##### Calculate the distribution
# ### Make a table to summarize
# types <- as.character(unique(reads.ct$type))
# length(types)
# sum.table <- data.frame(row.names = types, Dbp21 = rep(0, length(types)), Dbp22 = rep(0, length(types)), Dbp23 = rep(0, length(types)), Dbp21_percent = rep(0, length(types)), Dbp22_percent = rep(0, length(types)), Dbp23_percent = rep(0, length(types)))
# ### Count mRNAs
# for (i in 1:3){
#     sum.table[1,i] <- sum(reads.ct[which(reads.ct$type=="gene"),i])
# }
# for (i in 1:3) {
#     sum.table[1,i+3] <- sum.table[1,i]/colSums(reads.ct[,i])*100
# }
# ### Count ncRNA_gene
# for (i in 1:3){
#     sum.table[2,i] <- sum(reads.ct[which(reads.ct$type=="ncRNA_gene"),i])
# }
# for (i in 1:3) {
#     sum.table[2,i+3] <- sum.table[2,i]/colSums(reads.ct[,i])*100
# }
# ### Count rRNA_gene
# for (i in 1:3){
#     sum.table[3,i] <- sum(reads.ct[which(reads.ct$type=="rRNA_gene"),i])
# }
# for (i in 1:3) {
#     sum.table[3,i+3] <- sum.table[3,i]/colSums(reads.ct[,i])*100
# }
# ### Count snoRNA_gene
# for (i in 1:3){
#     sum.table[4,i] <- sum(reads.ct[which(reads.ct$type=="snoRNA_gene"),i])
# }
# for (i in 1:3) {
#     sum.table[4,i+3] <- sum.table[4,i]/colSums(reads.ct[,i])*100
# }
# ### Count snRNA_gene
# for (i in 1:3){
#     sum.table[5,i] <- sum(reads.ct[which(reads.ct$type=="snRNA_gene"),i])
# }
# for (i in 1:3) {
#     sum.table[5,i+3] <- sum.table[5,i]/colSums(reads.ct[,i])*100
# }
# ### Count telomerase_RNA_gene
# for (i in 1:3){
#     sum.table[6,i] <- sum(reads.ct[which(reads.ct$type=="telomerase_RNA_gene"),i])
# }
# for (i in 1:3) {
#     sum.table[6,i+3] <- sum.table[6,i]/colSums(reads.ct[,i])*100
# }
# ### Count tRNA_gene
# for (i in 1:3){
#     sum.table[7,i] <- sum(reads.ct[which(reads.ct$type=="tRNA_gene"),i])
# }
# for (i in 1:3) {
#     sum.table[7,i+3] <- sum.table[7,i]/colSums(reads.ct[,i])*100
# }
# ### Count lncRNA 
# for (i in 1:3){
#     sum.table[8,i] <- sum(reads.ct[which(reads.ct$type=="lncRNA"),i])
# }
# for (i in 1:3) {
#     sum.table[8,i+3] <- sum.table[8,i]/colSums(reads.ct[,i])*100
# }
# ### Sum ncRNA_gene, telomerase_RNA_gene, and lncRNA together for lncRNA_total
# lncRNA_total <- data.frame(row.names = "lncRNA_total", Dbp21 = 0, Dbp22 = 0, Dbp23 = 0, Dbp21_percent = 0, Dbp22_percent = 0, Dbp23_percent = 0)
# for (i in 1:3) {
#     lncRNA_total[,i] <- sum(sum.table[c(2,6,8),i])
# }
# for (i in 1:3) {
#     lncRNA_total[,i+3] <- lncRNA_total[1,i]/colSums(reads.ct[,i])*100
# }
# sum.table <- rbind(sum.table, lncRNA_total)
# ### Count the average
# for (i in 1:length(sum.table$Dbp21)) {
#     sum.table$average_percent[i] <- rowMeans(sum.table[i,4:6])
# }
# 
# ##### Output the table
# write.table(reads.ct, "~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/091317_Dbp2_xlinking_sites_filtering/091317_Dbp2_iCLIP_intersected_filtered_reads_count_mat.txt", quote = F, sep = "\t")
# write.table(sum.table, "~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/091317_Dbp2_xlinking_sites_filtering/091317_Dbp2_iCLIP_intersected_filtered_reads_distribution_summary.txt", quote = F, sep = "\t")


#----------# Count the "sites"
##### Make GRanges of "sites"
d1.sites.gr <- GRanges(seqnames = Rle(d1.filt$chr), ranges = IRanges(start = d1.filt$D1_xlinking_site, end = d1.filt$D1_xlinking_site), strand = Rle(d1.filt$strand), id = d1.filt$id, seqlengths = yeast.length)
d2.sites.gr <- GRanges(seqnames = Rle(d2.filt$chr), ranges = IRanges(start = d2.filt$D2_xlinking_site, end = d2.filt$D2_xlinking_site), strand = Rle(d2.filt$strand), id = d2.filt$id, seqlengths = yeast.length)
d3.sites.gr <- GRanges(seqnames = Rle(d3.filt$chr), ranges = IRanges(start = d3.filt$D3_xlinking_site, end = d3.filt$D3_xlinking_site), strand = Rle(d3.filt$strand), id = d3.filt$id, seqlengths = yeast.length)

##### Use summarizeOverlap to count the "sites"
### Generate an empty matrix first
sites.ct <- data.frame(row.names = mcols(all.gr)$id, Dbp21 = NA, Dbp22 = NA, Dbp23 = NA, type = mcols(all.gr)$type)
# sites.ct <- matrix(nrow = length(all.gr), ncol = 4)
# rownames(sites.ct) <- mcols(all.gr)$id
colnames(sites.ct) <- c("Dbp2-1", "Dbp2-2", "Dbp2-3", "type")
# sites.ct[,4] <- mcols(all.gr)$type
### start counting
sites.ct[,1] <- assays(summarizeOverlaps(all.gr, d1.sites.gr, ignore.strand = F, mode = "IntersectionNotEmpty"))$counts
sites.ct[,2] <- assays(summarizeOverlaps(all.gr, d2.sites.gr, ignore.strand = F, mode = "IntersectionNotEmpty"))$counts
sites.ct[,3] <- assays(summarizeOverlaps(all.gr, d3.sites.gr, ignore.strand = F, mode = "IntersectionNotEmpty"))$counts

colSums(sites.ct[,1:3]) ### This method resutls in less counts (likely because some sites are in the overlapping region of multiple transcripts)
### Dbp2-1  Dbp2-2  Dbp2-3 
### 882199  904585 1006183 

##### Calculate the distribution
### Make a table to summarize
types <- as.character(unique(sites.ct$type))
length(types)
sites.sum.table <- data.frame(row.names = types, Dbp21 = rep(0, length(types)), Dbp22 = rep(0, length(types)), Dbp23 = rep(0, length(types)), Dbp21_percent = rep(0, length(types)), Dbp22_percent = rep(0, length(types)), Dbp23_percent = rep(0, length(types)))
### Count mRNAs
for (i in 1:3){
    sites.sum.table[1,i] <- sum(sites.ct[which(sites.ct$type=="gene"),i])
}
for (i in 1:3) {
    sites.sum.table[1,i+3] <- sites.sum.table[1,i]/colSums(sites.ct[,i])*100
}
### Count ncRNA_gene
for (i in 1:3){
    sites.sum.table[2,i] <- sum(sites.ct[which(sites.ct$type=="ncRNA_gene"),i])
}
for (i in 1:3) {
    sites.sum.table[2,i+3] <- sites.sum.table[2,i]/colSums(sites.ct[,i])*100
}
### Count rRNA_gene
for (i in 1:3){
    sites.sum.table[3,i] <- sum(sites.ct[which(sites.ct$type=="rRNA_gene"),i])
}
for (i in 1:3) {
    sites.sum.table[3,i+3] <- sites.sum.table[3,i]/colSums(sites.ct[,i])*100
}
### Count snoRNA_gene
for (i in 1:3){
    sites.sum.table[4,i] <- sum(sites.ct[which(sites.ct$type=="snoRNA_gene"),i])
}
for (i in 1:3) {
    sites.sum.table[4,i+3] <- sites.sum.table[4,i]/colSums(sites.ct[,i])*100
}
### Count snRNA_gene
for (i in 1:3){
    sites.sum.table[5,i] <- sum(sites.ct[which(sites.ct$type=="snRNA_gene"),i])
}
for (i in 1:3) {
    sites.sum.table[5,i+3] <- sites.sum.table[5,i]/colSums(sites.ct[,i])*100
}
### Count telomerase_RNA_gene
for (i in 1:3){
    sites.sum.table[6,i] <- sum(sites.ct[which(sites.ct$type=="telomerase_RNA_gene"),i])
}
for (i in 1:3) {
    sites.sum.table[6,i+3] <- sites.sum.table[6,i]/colSums(sites.ct[,i])*100
}
### Count tRNA_gene
for (i in 1:3){
    sites.sum.table[7,i] <- sum(sites.ct[which(sites.ct$type=="tRNA_gene"),i])
}
for (i in 1:3) {
    sites.sum.table[7,i+3] <- sites.sum.table[7,i]/colSums(sites.ct[,i])*100
}
### Count lncRNA 
for (i in 1:3){
    sites.sum.table[8,i] <- sum(sites.ct[which(sites.ct$type=="lncRNA"),i])
}
for (i in 1:3) {
    sites.sum.table[8,i+3] <- sites.sum.table[8,i]/colSums(sites.ct[,i])*100
}
### Sum ncRNA_gene, telomerase_RNA_gene, and lncRNA together for lncRNA_total
lncRNA_total <- data.frame(row.names = "lncRNA_total", Dbp21 = 0, Dbp22 = 0, Dbp23 = 0, Dbp21_percent = 0, Dbp22_percent = 0, Dbp23_percent = 0)
for (i in 1:3) {
    lncRNA_total[,i] <- sum(sites.sum.table[c(2,6,8),i])
}
for (i in 1:3) {
    lncRNA_total[,i+3] <- lncRNA_total[1,i]/colSums(sites.ct[,i])*100
}
sites.sum.table <- rbind(sites.sum.table, lncRNA_total)
### Count the average
for (i in 1:length(sites.sum.table$Dbp21)) {
    sites.sum.table$average_percent[i] <- rowMeans(sites.sum.table[i,4:6])
}

##### Output the tables
write.table(sites.ct, "../Dbp2iCLIP_sitesFromReproducedReads/060618_Dbp2iCLIP_SitesFromIntersectedReads_countMatrix.txt",
            quote = F, sep = "\t")
write.table(sites.sum.table, "../Dbp2iCLIP_sitesFromReproducedReads/060618_Dbp2iCLIP_SitesFromIntersectedReads_distributionSummary.txt",
            quote = F, sep = "\t")
### Proceed to the following analysis using the result from this method

#----------# Only keep RNAs with > 5 counts in each library for the following analysis
library(edgeR)
##### This part can be done on my computer
##### Read the files in
sites.ct <- read.delim("~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/060618_Dbp2iCLIPdata_reprocessing/060618_Dbp2iCLIP_SitesFromIntersectedReads_countMatrix.txt",
                       as.is = T)

##### Only keep RNAs with > 5 counts in ALL three samples
### Convert the table into a matrix for calculation
sites.ct.mat <- as.matrix(sites.ct[,1:3])
### Filtering
gene.keep <- rowSums(sites.ct.mat > 5) >= 3
sites.ct.mat.filt <- sites.ct.mat[gene.keep,] ### 2391 transcripts left

##### Find out the number of mRNAs bound
sites.ct.filt <- sites.ct[which(rownames(sites.ct)%in%rownames(sites.ct.mat.filt)),]
type.tb <- table(sites.ct.filt$type)
### Output this table
write.table(type.tb,
            "~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/060618_Dbp2iCLIPdata_reprocessing/060618_Dbp2iCLIP_SitesFromIntersectedReads_RNAtype_number.txt",
            quote = F, sep = "\t", row.names = F)

##### Calculate the distribution
### Make a table to summarize
types <- as.character(unique(sites.ct.filt$type))
length(types)
sites.sum.table <- data.frame(row.names = types, Dbp21 = rep(0, length(types)),
                              Dbp22 = rep(0, length(types)), Dbp23 = rep(0, length(types)),
                              Dbp21_percent = rep(0, length(types)),
                              Dbp22_percent = rep(0, length(types)),
                              Dbp23_percent = rep(0, length(types)))

### Count mRNAs
for (i in 1:3){
    sites.sum.table[1,i] <- sum(sites.ct.filt[which(sites.ct.filt$type=="gene"),i])
}
for (i in 1:3) {
    sites.sum.table[1,i+3] <- sites.sum.table[1,i]/sum(sites.ct.filt[,i])*100
}
### Count ncRNA_gene
for (i in 1:3){
    sites.sum.table[2,i] <- sum(sites.ct.filt[which(sites.ct.filt$type=="ncRNA_gene"),i])
}
for (i in 1:3) {
    sites.sum.table[2,i+3] <- sites.sum.table[2,i]/sum(sites.ct.filt[,i])*100
}
### Count rRNA_gene
for (i in 1:3){
    sites.sum.table[3,i] <- sum(sites.ct.filt[which(sites.ct.filt$type=="rRNA_gene"),i])
}
for (i in 1:3) {
    sites.sum.table[3,i+3] <- sites.sum.table[3,i]/sum(sites.ct.filt[,i])*100
}
### Count snoRNA_gene
for (i in 1:3){
    sites.sum.table[4,i] <- sum(sites.ct.filt[which(sites.ct.filt$type=="snoRNA_gene"),i])
}
for (i in 1:3) {
    sites.sum.table[4,i+3] <- sites.sum.table[4,i]/sum(sites.ct.filt[,i])*100
}
### Count snRNA_gene
for (i in 1:3){
    sites.sum.table[5,i] <- sum(sites.ct.filt[which(sites.ct.filt$type=="snRNA_gene"),i])
}
for (i in 1:3) {
    sites.sum.table[5,i+3] <- sites.sum.table[5,i]/sum(sites.ct.filt[,i])*100
}
### Count telomerase_RNA_gene
for (i in 1:3){
    sites.sum.table[6,i] <- sum(sites.ct.filt[which(sites.ct.filt$type=="telomerase_RNA_gene"),i])
}
for (i in 1:3) {
    sites.sum.table[6,i+3] <- sites.sum.table[6,i]/sum(sites.ct.filt[,i])*100
}
### Count tRNA_gene
for (i in 1:3){
    sites.sum.table[7,i] <- sum(sites.ct.filt[which(sites.ct.filt$type=="tRNA_gene"),i])
}
for (i in 1:3) {
    sites.sum.table[7,i+3] <- sites.sum.table[7,i]/sum(sites.ct.filt[,i])*100
}
### Count lncRNA 
for (i in 1:3){
    sites.sum.table[8,i] <- sum(sites.ct.filt[which(sites.ct.filt$type=="lncRNA"),i])
}
for (i in 1:3) {
    sites.sum.table[8,i+3] <- sites.sum.table[8,i]/sum(sites.ct.filt[,i])*100
}
### Sum ncRNA_gene, telomerase_RNA_gene, and lncRNA together for lncRNA_total
lncRNA_total <- data.frame(row.names = "lncRNA_total", Dbp21 = 0, Dbp22 = 0, Dbp23 = 0,
                           Dbp21_percent = 0, Dbp22_percent = 0, Dbp23_percent = 0)
for (i in 1:3) {
    lncRNA_total[,i] <- sum(sites.sum.table[c(2,6,8),i])
}
for (i in 1:3) {
    lncRNA_total[,i+3] <- lncRNA_total[1,i]/sum(sites.ct.filt[,i])*100
}
sites.sum.table <- rbind(sites.sum.table, lncRNA_total)
### Count the average
for (i in 1:length(sites.sum.table$Dbp21)) {
    sites.sum.table$average_percent[i] <- rowMeans(sites.sum.table[i,4:6])
}

##### Output the tables
write.table(sites.ct.filt,
            "~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/060618_Dbp2iCLIPdata_reprocessing/060718_Dbp2iCLIP_SitesFromIntersectedReads_countMatrixFiltered.txt",
            quote = F, sep = "\t")
write.table(sites.sum.table, "~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/060618_Dbp2iCLIPdata_reprocessing/060718_Dbp2iCLIP_SitesFromIntersectedReads_distributionSummaryFiltered.txt",
            quote = F, sep = "\t")

##### Get the list of mRNAs that are bound
mrnas.left <- rownames(sites.ct.filt)[which(sites.ct.filt$type=="gene")]
length(mrnas.left) ### 2057 mRNAs are identified as binding targets
### Only keep those that are expressed based on structure-seq controls
exp.genes <- read.delim("~/Desktop/Research_projects_2015-2016/DMS/Structure-seq_data/052417_ORF_DE_analysis_without_1153/052617_gene_list_passed_cpm_filter.txt",
                        header = F, as.is = T)
mrnas.left.exp <- mrnas.left[which(mrnas.left%in%exp.genes$V1)]
length(mrnas.left.exp) ### 2045

### Output the list and counts
sites.ct.mat.filt.mrnas <- sites.ct.mat.filt[which(rownames(sites.ct.mat.filt)%in%mrnas.left.exp),]
write.table(sites.ct.mat.filt.mrnas,
            "~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/060618_Dbp2iCLIPdata_reprocessing/060718_Dbp2iCLIP_SitesFromIntersectedReads_countMat_for_kept_mRNAs.txt",
            quote = F, sep = "\t")


# snoRNA.left <- rownames(sites.ct.filt)[which(sites.ct.filt$type=="snoRNA_gene")]
# length(snoRNA.left) ### 64 snoRNAs are targets

