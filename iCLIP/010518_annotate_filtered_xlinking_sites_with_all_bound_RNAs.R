
##########
##### This is to annotate the filtered xlinking sites from all intersected reads with mRNA transcript regions (including UTRs)
##########
##### Xlinking sites are derived from all intersected reads in each replicate
### Each replicate is analyzed separately because I want to get the genes that are targeted by all 3 replicates and filter out those are not reproduced in all 3 samples
### The binding sites used here are re-filtered in Dec2018 by Nadia

library(Rsamtools)
library(GenomicRanges)
library(GenomicAlignments)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
sc.bs3 <- BSgenome.Scerevisiae.UCSC.sacCer3
yeast.length <- seqlengths(sc.bs3)

##### Read in the targeted mRNA from all intersected transcripts
mrna.targets <- read.delim("~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/010418_Data_from_Nadia/091317_Dbp2_iCLIP_intersected_filtered_sites_count_mat_for_kept_mRNAs.txt", as.is = T)

##### Read in all the crosslinking sites in each replicate
D1.sites <- read.delim("~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/010418_Data_from_Nadia/091317_iCLIP_Dbp21_intersected_reads_xlinking_sites_filtered.txt", as.is = T)
D2.sites <- read.delim("~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/010418_Data_from_Nadia/091317_iCLIP_Dbp22_intersected_reads_xlinking_sites_filtered.txt", as.is = T)
D3.sites <- read.delim("~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/010418_Data_from_Nadia/091317_iCLIP_Dbp23_intersected_reads_xlinking_sites_filtered.txt", as.is = T)

##### Make GRanges of xlinking sites
D1.sites.gr <- GRanges(seqnames = Rle(D1.sites$chr), ranges = IRanges(start = D1.sites$D1_xlinking_site, end = D1.sites$D1_xlinking_site), strand = Rle(D1.sites$strand), seqlengths = yeast.length)
D2.sites.gr <- GRanges(seqnames = Rle(D2.sites$chr), ranges = IRanges(start = D2.sites$D2_xlinking_site, end = D2.sites$D2_xlinking_site), strand = Rle(D2.sites$strand), seqlengths = yeast.length)
D3.sites.gr <- GRanges(seqnames = Rle(D3.sites$chr), ranges = IRanges(start = D3.sites$D3_xlinking_site, end = D3.sites$D3_xlinking_site), strand = Rle(D3.sites$strand), seqlengths = yeast.length)

##### Load the annotation file with UTR information
yn.utr.s <- read.delim("~/Desktop/Research_projects_2015-2016/Data_general_use/062216_Yassour_Nagalakshmi_outermost_UTR.txt", quote = "", as.is = T)
rna.gr <- GRanges(seqnames = Rle(yn.utr.s$chromosome), ranges = IRanges(start = yn.utr.s$transcript.start, end = yn.utr.s$transcript.end), strand = Rle(yn.utr.s$strand), id = yn.utr.s$id, seqlengths = yeast.length)
rna.gr <-trim(rna.gr)

##### Annotate the sites
### Dbp2-1
D1.sites.index <- findOverlaps(rna.gr, D1.sites.gr, ignore.strand = F)
D1.sites.index.df <- data.frame(rna_index = queryHits(D1.sites.index), D1_index = subjectHits(D1.sites.index))
D1.sites.index.df$chr <- (as.character(seqnames(rna.gr)))[D1.sites.index.df$rna_index]
D1.sites.rna.range <- ranges(rna.gr)[D1.sites.index.df$rna_index]
D1.sites.rna.range.df <- as.data.frame(D1.sites.rna.range, as.is = T)
D1.sites.index.df$transcript.start <- D1.sites.rna.range.df$start
D1.sites.index.df$transcript.end <- D1.sites.rna.range.df$end
D1.sites.index.df$strand <- (as.character(strand(rna.gr)))[D1.sites.index.df$rna_index]
D1.sites.index.df$gene_id <- (mcols(rna.gr)$id)[D1.sites.index.df$rna_index]
D1.sites.range <- ranges(D1.sites.gr)[D1.sites.index.df$D1_index]
D1.sites.range.df <- as.data.frame(D1.sites.range, as.is=T)
D1.sites.index.df$D1_xlinking_site <- D1.sites.range.df$start
# D1.sites.index.df$counts <- (mcols(D1.sites.gr)$count)[D1.sites.index.df$D1_index]
dim(D1.sites.index.df) ### 340984
length(unique(D1.sites.index.df$D1_index)) ### 330402, this indicates some sites are counted twice because they overlap with more than one transcripts
## Need to decide which transcript those sites belong to manually. Keep the transcript if the site is closer to the transcript's ORF
yn.utr <- read.delim("~/Desktop/Research_projects_2015-2016/Data_general_use/071016_Yassour_Nagalakshmi_outermost_UTR.txt", quote = "", as.is = T)
yn.utr.keep <- yn.utr[, c(1,3,4,5,6,7)]
yn.utr.keep2<- yn.utr.keep[,c(2,5,6)]
## Start assigning
D1.dup <- D1.sites.index.df[duplicated(D1.sites.index.df$D1_index), ]
D1.uni <- D1.sites.index.df[-which(D1.sites.index.df$D1_index%in%D1.dup$D1_index), ]
dim(D1.uni) ### 321448 sites
D1.dup.sites <- D1.sites.index.df[D1.sites.index.df$D1_index%in%D1.dup$D1_index, ]
D1.dup.order <- D1.dup.sites[order(D1.dup.sites[ ,2]), ]
D1.dup.order.orf <- merge(D1.dup.order, yn.utr.keep2, by.x = c("gene_id"), by.y = c("id"), sort = F)
D1.dup.order.orf <- D1.dup.order.orf[order(D1.dup.order.orf[ ,3]), ]
a <- D1.dup.order.orf
dim(a) ### 19536 entries

a.keep.fun <- function(x){
    if (a$D1_index[x]==a$D1_index[x-1]){
        xd <- distance(IRanges(start = a$D1_xlinking_site[x]-1, end = a$D1_xlinking_site[x]), IRanges(start = a$orf.start[x], end = a$orf.end[x]))
        x1d <- distance(IRanges(start = a$D1_xlinking_site[x-1]-1, end = a$D1_xlinking_site[x-1]), IRanges(start = a$orf.start[x-1], end = a$orf.end[x-1]))
        if (xd>x1d){y <- a[x-1, ]} else{
            if (xd<x1d){
                y <- a[x, ]
            } 
        }
    }
}

D1.dup.keep <- lapply(2:length(a$gene_id), a.keep.fun)

y <- D1.dup.keep[[1]]
for (i in 2:length(D1.dup.keep)){
    y <- rbind(y, D1.dup.keep[[i]])
}

D1.dup.keep.df <- unique(y)
dim(D1.dup.keep.df) ### 7841 entries after cleaning up 
D1.dup.keep.df <- D1.dup.keep.df[,c(2,3,4,5,6,7,1,8)]
D1.final <- rbind(D1.uni, D1.dup.keep.df)
dim(D1.final) ### 329289 sites are mapped to mRNA transcripts
## Get the genelist for D1
D1.sites.genelist <- unique(D1.final$gene_id)
length(D1.sites.genelist) ### 3771
## Only keep sites/genes that match the target mRNAs (derived in Sep.2017)
D1.final.filt <- D1.final[which(D1.final$gene_id%in%rownames(mrna.targets)),]
D1.final.filt.genes <- unique(D1.final.filt$gene_id)
length(D1.final.filt.genes) ### 2036 transcripts are really xlinked in Dbp2-1 replicate
## Output the table
write.table(D1.final.filt, file = "~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/010418_Data_from_Nadia/010518_Dbp21_iCLIP_filtered_sites_in_target_mRNA_transcripts.txt", quote = F, sep = "\t", row.names = F)

### D2
D2.sites.index <- findOverlaps(rna.gr, D2.sites.gr, ignore.strand = F)
D2.sites.index.df <- data.frame(rna_index = queryHits(D2.sites.index), D2_index = subjectHits(D2.sites.index))
D2.sites.index.df$chr <- (as.character(seqnames(rna.gr)))[D2.sites.index.df$rna_index]
D2.sites.rna.range <- ranges(rna.gr)[D2.sites.index.df$rna_index]
D2.sites.rna.range.df <- as.data.frame(D2.sites.rna.range, as.is = T)
D2.sites.index.df$transcript.start <- D2.sites.rna.range.df$start
D2.sites.index.df$transcript.end <- D2.sites.rna.range.df$end
D2.sites.index.df$strand <- (as.character(strand(rna.gr)))[D2.sites.index.df$rna_index]
D2.sites.index.df$gene_id <- (mcols(rna.gr)$id)[D2.sites.index.df$rna_index]
D2.sites.range <- ranges(D2.sites.gr)[D2.sites.index.df$D2_index]
D2.sites.range.df <- as.data.frame(D2.sites.range, as.is=T)
D2.sites.index.df$D2_xlinking_sites <- D2.sites.range.df$end
D2.sites.index.df$counts <- (mcols(D2.sites.gr)$count)[D2.sites.index.df$D2_index]
dim(D2.sites.index.df) ### 348583
length(unique(D2.sites.index.df$D2_index)) ### 337598, this indicates some sites are counted twice because they overlap with more than one transcripts
## Need to decide which transcript those sites belong to manually. Keep the transcript if the site is closer to the transcript's ORF
## Start assigning
D2.dup <- D2.sites.index.df[duplicated(D2.sites.index.df$D2_index), ]
D2.uni <- D2.sites.index.df[-which(D2.sites.index.df$D2_index%in%D2.dup$D2_index), ]
dim(D2.uni) ### 328708 sites
D2.dup.sites <- D2.sites.index.df[D2.sites.index.df$D2_index%in%D2.dup$D2_index, ]
D2.dup.order <- D2.dup.sites[order(D2.dup.sites[ ,2]), ]
D2.dup.order.orf <- merge(D2.dup.order, yn.utr.keep2, by.x = c("gene_id"), by.y = c("id"), sort = F)
D2.dup.order.orf <- D2.dup.order.orf[order(D2.dup.order.orf[ ,3]), ]
a <- D2.dup.order.orf
dim(a) ### 19875 entries

a.keep.fun <- function(x){
    if (a$D2_index[x]==a$D2_index[x-1]){
        xd <- distance(IRanges(start = a$D2_xlinking_site[x]-1, end = a$D2_xlinking_site[x]), IRanges(start = a$orf.start[x], end = a$orf.end[x]))
        x1d <- distance(IRanges(start = a$D2_xlinking_site[x-1]-1, end = a$D2_xlinking_site[x-1]), IRanges(start = a$orf.start[x-1], end = a$orf.end[x-1]))
        if (xd>x1d){y <- a[x-1, ]} else{
            if (xd<x1d){y <- a[x, ]} 
        }
    }
}

D2.dup.keep <- lapply(2:length(a$gene_id), a.keep.fun)

y <- D2.dup.keep[[1]]
for (i in 2:length(D2.dup.keep)){
    y <- rbind(y, D2.dup.keep[[i]])
}

D2.dup.keep.df <- unique(y)
dim(D2.dup.keep.df) ### 7320 entries after cleaning up 
D2.dup.keep.df <- D2.dup.keep.df[,c(2,3,4,5,6,7,1,8)]
D2.final <- rbind(D2.uni, D2.dup.keep.df)
dim(D2.final) ### 336028 sites are mapped to mRNA transcripts
## Get the genelist for D2
D2.sites.genelist <- unique(D2.final$gene_id)
length(D2.sites.genelist) ### 3622
## Only keep sites/genes that match the target mRNAs (by intersected reads after filtering)
D2.final.filt <- D2.final[which(D2.final$gene_id%in%rownames(mrna.targets)),]
D2.final.filt.genes <- unique(D2.final.filt$gene_id)
length(D2.final.filt.genes) ### 2036 transcripts are really xlinked in Dbp2-3 replicate
## Output the table
write.table(D2.final.filt, file = "~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/010418_Data_from_Nadia/010518_Dbp22_iCLIP_filtered_sites_in_target_mRNA_transcripts.txt", quote = F, sep = "\t", row.names = F)

### D3
D3.sites.index <- findOverlaps(rna.gr, D3.sites.gr, ignore.strand = F)
D3.sites.index.df <- data.frame(rna_index = queryHits(D3.sites.index), D3_index = subjectHits(D3.sites.index))
D3.sites.index.df$chr <- (as.character(seqnames(rna.gr)))[D3.sites.index.df$rna_index]
D3.sites.rna.range <- ranges(rna.gr)[D3.sites.index.df$rna_index]
D3.sites.rna.range.df <- as.data.frame(D3.sites.rna.range, as.is = T)
D3.sites.index.df$transcript.start <- D3.sites.rna.range.df$start
D3.sites.index.df$transcript.end <- D3.sites.rna.range.df$end
D3.sites.index.df$strand <- (as.character(strand(rna.gr)))[D3.sites.index.df$rna_index]
D3.sites.index.df$gene_id <- (mcols(rna.gr)$id)[D3.sites.index.df$rna_index]
D3.sites.range <- ranges(D3.sites.gr)[D3.sites.index.df$D3_index]
D3.sites.range.df <- as.data.frame(D3.sites.range, as.is=T)
D3.sites.index.df$D3_xlinking_sites <- D3.sites.range.df$end
D3.sites.index.df$counts <- (mcols(D3.sites.gr)$count)[D3.sites.index.df$D3_index]
dim(D3.sites.index.df) ### 365254
length(unique(D3.sites.index.df$D3_index)) ### 351321, this indicates some sites are counted twice because they overlap with more than one transcripts
## Need to decide which transcript those sites belong to manually. Keep the transcript if the site is closer to the transcript's ORF
## Start assigning
D3.dup <- D3.sites.index.df[duplicated(D3.sites.index.df$D3_index), ]
D3.uni <- D3.sites.index.df[-which(D3.sites.index.df$D3_index%in%D3.dup$D3_index), ]
dim(D3.uni) ### 340682 sites
D3.dup.sites <- D3.sites.index.df[D3.sites.index.df$D3_index%in%D3.dup$D3_index, ]
D3.dup.order <- D3.dup.sites[order(D3.dup.sites[ ,2]), ]
D3.dup.order.orf <- merge(D3.dup.order, yn.utr.keep2, by.x = c("gene_id"), by.y = c("id"), sort = F)
D3.dup.order.orf <- D3.dup.order.orf[order(D3.dup.order.orf[ ,3]), ]
a <- D3.dup.order.orf
dim(a) ### 24572 entries

a.keep.fun <- function(x){
    if (a$D3_index[x]==a$D3_index[x-1]){
        xd <- distance(IRanges(start = a$D3_xlinking_site[x]-1, end = a$D3_xlinking_site[x]), IRanges(start = a$orf.start[x], end = a$orf.end[x]))
        x1d <- distance(IRanges(start = a$D3_xlinking_site[x-1]-1, end = a$D3_xlinking_site[x-1]), IRanges(start = a$orf.start[x-1], end = a$orf.end[x-1]))
        if (xd>x1d){y <- a[x-1, ]} else{
            if (xd<x1d){
                y <- a[x, ]
            } 
        }
    }
}

D3.dup.keep <- lapply(2:length(a$gene_id), a.keep.fun)

y <- D3.dup.keep[[1]]
for (i in 2:length(D3.dup.keep)){
    y <- rbind(y, D3.dup.keep[[i]])
}

D3.dup.keep.df <- unique(y)
dim(D3.dup.keep.df) ### 8890 entries after cleaning up 
D3.dup.keep.df <- D3.dup.keep.df[,c(2,3,4,5,6,7,1,8)]
D3.final <- rbind(D3.uni, D3.dup.keep.df)
dim(D3.final) ### 349572 sites are mapped to mRNA transcripts
## Get the genelist for D3
D3.sites.genelist <- unique(D3.final$gene_id)
length(D3.sites.genelist) ### 3729
## Only keep sites/genes that match the target mRNAs (by intersected reads after filtering)
D3.final.filt <- D3.final[which(D3.final$gene_id%in%rownames(mrna.targets)),]
D3.final.filt.genes <- unique(D3.final.filt$gene_id)
length(D3.final.filt.genes) ### 2036 transcripts are really xlinked in Dbp2-3 replicate
## Output the table
# Record of xlinking sites in all mRNA transcripts
write.table(D3.final.filt, file = "~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/010418_Data_from_Nadia/010518_Dbp23_iCLIP_filtered_sites_in_target_mRNA_transcripts.txt", quote = F, sep = "\t", row.names = F)

# ##### Get the genes that appear in all the three erplicates
# inter12.genes <- intersect(D1.final.filt.genes, D2.final.filt.genes)
# inter123.genes <- intersect(inter12.genes, D3.final.filt.genes)
# length(inter123.genes) ### 2058 transcripts
# ### Output the genes that are found in all three Dbp2-iCLIP replicates
# write.table(inter123.genes, file = "~/Desktop/Research projects 2015-2016/iCLIP/results_from_scripts/011717_Dbp2_iCLIP_reads_refiltering/011817_iCLIP_Dbp2_xlinking_sites_in_bound_mRNA_transcripts/011817_iCLIP_Dbp2_binding_mRNA_transcript_list.txt", quote = F, sep = "\t", row.names = F, col.names = F)

# ##### Get the information of xlinking sites in these transcripts in each replicate
# D1.final.filt2 <- D1.final.filt[which(D1.final.filt$gene_id%in%inter123.genes), ]
# D1.final.filt2 <- D1.final.filt2[,c(3,7,6,4,5,8,9)]
# D2.final.filt2 <- D2.final.filt[which(D2.final.filt$gene_id%in%inter123.genes),]
# D2.final.filt2 <- D2.final.filt2[,c(3,7,6,4,5,8,9)]
# D3.final.filt2 <- D3.final.filt[which(D3.final.filt$gene_id%in%inter123.genes),]
# D3.final.filt2 <- D3.final.filt2[,c(3,7,6,4,5,8,9)]
# ### Output the tables
# write.table(D1.final.filt2, file = "~/Desktop/Research projects 2015-2016/iCLIP/results_from_scripts/011717_Dbp2_iCLIP_reads_refiltering/011817_iCLIP_Dbp2_xlinking_sites_in_bound_mRNA_transcripts/011817_iCLIP_Dbp21_xlinking_sites_in_bound_mRNA_transcript.txt", quote = F, sep = "\t", row.names = F)
# write.table(D2.final.filt2, file = "~/Desktop/Research projects 2015-2016/iCLIP/results_from_scripts/011717_Dbp2_iCLIP_reads_refiltering/011817_iCLIP_Dbp2_xlinking_sites_in_bound_mRNA_transcripts/011817_iCLIP_Dbp22_xlinking_sites_in_bound_mRNA_transcript.txt", quote = F, sep = "\t", row.names = F)
# write.table(D3.final.filt2, file = "~/Desktop/Research projects 2015-2016/iCLIP/results_from_scripts/011717_Dbp2_iCLIP_reads_refiltering/011817_iCLIP_Dbp2_xlinking_sites_in_bound_mRNA_transcripts/011817_iCLIP_Dbp23_xlinking_sites_in_bound_mRNA_transcript.txt", quote = F, sep = "\t", row.names = F)

