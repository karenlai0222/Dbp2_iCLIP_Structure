
##########
##### This is to partition filtered Dbp2 binding sites in mRNA targets into different regions
##########
### The binding sites used here are re-filtered in Dec2018 by Nadia

library(Rsamtools)
library(GenomicRanges)
library(GenomicAlignments)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
sc.bs3 <- BSgenome.Scerevisiae.UCSC.sacCer3
yeast.length <- seqlengths(sc.bs3)

#----------# This part is to find out what transcripts are bound in the 3' UTR and the 50 nt at the 3' end of ORF
##### Read in the filtered sites in bound mRNAs
d1 <- read.delim("~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/010418_Data_from_Nadia/010518_Dbp21_iCLIP_filtered_sites_in_target_mRNA_transcripts.txt", as.is = T)
d2 <- read.delim("~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/010418_Data_from_Nadia/010518_Dbp22_iCLIP_filtered_sites_in_target_mRNA_transcripts.txt", as.is = T)
d3 <- read.delim("~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/010418_Data_from_Nadia/010518_Dbp23_iCLIP_filtered_sites_in_target_mRNA_transcripts.txt", as.is = T)

##### Make GRanges for these sites
d1.gr <- GRanges(seqnames = Rle(d1$chr), ranges = IRanges(start = d1$D1_xlinking_site, end = d1$D1_xlinking_site), strand = Rle(d1$strand), id = d1$gene_id, seqlengths = yeast.length)
d2.gr <- GRanges(seqnames = Rle(d2$chr), ranges = IRanges(start = d2$D2_xlinking_site, end = d2$D2_xlinking_site), strand = Rle(d2$strand), id = d2$gene_id, seqlengths = yeast.length)
d3.gr <- GRanges(seqnames = Rle(d3$chr), ranges = IRanges(start = d3$D3_xlinking_site, end = d3$D3_xlinking_site), strand = Rle(d3$strand), id = d3$gene_id, seqlengths = yeast.length)

##### Read the UTR information in 
yn.utr.k <- read.delim("~/Desktop/Research_projects_2015-2016/Data_general_use/071016_Yassour_Nagalakshmi_outermost_UTR.txt", as.is = T)
### Only keep those that are bound
bound.mrna <- read.delim("~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/010418_Data_from_Nadia/091317_Dbp2_iCLIP_intersected_filtered_sites_count_mat_for_kept_mRNAs.txt", as.is = T)
yn.utr.k.keep <- yn.utr.k[which(yn.utr.k$id%in%rownames(bound.mrna)),]

##### Make GRanges for the 3' UTRs and the last 50 nt of ORFs
utr3.gr <- GRanges(seqnames = Rle(yn.utr.k.keep$chromosome), ranges = IRanges(start = yn.utr.k.keep$utr3.start, yn.utr.k.keep$utr3.end), strand = Rle(yn.utr.k.keep$strand), id = yn.utr.k.keep$id, seqlengths = yeast.length)
orf50nt.gr <- flank(utr3.gr, 50) ### This would create a range that includes 50 nt upstream of the 3' UTR ranges

##### Count binding sites in these two regions for each replicate
count.mat.list <- list(d1.mat <- matrix(nrow = length(utr3.gr), ncol = 2), d2.mat <- matrix(nrow = length(utr3.gr), ncol = 2), d3.mat <- matrix(nrow = length(utr3.gr), ncol = 2))
samples <- list(d1.gr, d2.gr, d3.gr)

for (i in 1:length(count.mat.list)){
    count.mat.list[[i]][,1] <- assays(summarizeOverlaps(orf50nt.gr, samples[[i]], ignore.strand = F, mode = "IntersectionNotEmpty"))$counts
    count.mat.list[[i]][,2] <- assays(summarizeOverlaps(utr3.gr, samples[[i]], ignore.strand = F, mode = "IntersectionNotEmpty"))$counts
    colnames(count.mat.list[[i]]) <- c("orf50nt", "3'utr")
    rownames(count.mat.list[[i]]) <- mcols(utr3.gr)$id
}

##### Only keep those that have counts in either range in all three replicates
count.mat.list.keep <- list()
for (i in 1:length(count.mat.list)) {
    count.mat.list.keep[[i]] <- count.mat.list[[i]][which(count.mat.list[[i]][,1]+count.mat.list[[i]][,2]!=0),]
}

##### Get the list of transcripts (only keep those that appear in all replicates)
d1.3prime.mrna <- rownames(count.mat.list.keep[[1]])
length(d1.3prime.mrna) ### 1086
d2.3prime.mrna <- rownames(count.mat.list.keep[[2]])
length(d2.3prime.mrna) ### 918
d3.3prime.mrna <- rownames(count.mat.list.keep[[3]])
length(d3.3prime.mrna) ### 1013
### Get the intersection
binding.3prime <- intersect(d1.3prime.mrna, d2.3prime.mrna)
length(binding.3prime) ### 893
binding.3prime <- intersect(binding.3prime, d3.3prime.mrna)
length(binding.3prime) ### 837
### What about the union?
binding.3prime.union <- union(d1.3prime.mrna, d2.3prime.mrna)
length(binding.3prime.union) ### 1111
binding.3prime.union <- union(binding.3prime.union, d3.3prime.mrna)
length(binding.3prime.union) ### 1140

##### Output the list (right now, use the union list)
write.table(binding.3prime.union, file = "~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/010418_Data_from_Nadia/010518_Dbp2_iCLIP_mRNAs_bound_at_3prime_end.txt", quote = F, row.names = F, col.names = F)
