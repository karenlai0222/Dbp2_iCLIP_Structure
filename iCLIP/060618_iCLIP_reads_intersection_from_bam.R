##########
##### This is to find the reads that are overlapping in all 3 replicates, based on the BAM files
##########
##### The purpose of using the original bam file is to take spliced/unspliced reads into consideration
### Also, the UTR regions are taken into account

##### This is run on cluster interactively

library(Rsamtools)
library(GenomicRanges)
library(GenomicAlignments)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
sc.bs3 <- BSgenome.Scerevisiae.UCSC.sacCer3
yeast.length <- seqlengths(sc.bs3)

##### Generate GRanges for protein coding genes (SGD R64-2-1)
yn.utr <- read.delim("CDS_UTR_coordinates.txt", as.is = T, header = T)
yn.utr.gr <- GRanges(seqnames = Rle(yn.utr$chromosome),
                     ranges = IRanges(start = yn.utr$transcript.start, end = yn.utr$transcript.end),
                     strand = Rle(yn.utr$strand), id = yn.utr$id, seqlengths = yeast.length)
yn.utr.gr <- trim(yn.utr.gr)

##### Also generate Granges for lncRNAs and all possible RNAs
all.rna <- read.delim("R64-2-1_yeast_id_name.txt", as.is = T)
all.rna.gr <- GRanges(seqnames = Rle(all.rna$chromosome),
                      ranges = IRanges(start = all.rna$start, end = all.rna$end),
                      strand = Rle(all.rna$strand), id = all.rna$id, seqlengths = yeast.length)
lncrna <- read.delim("All_lncRNA_merged.txt", as.is = T)
lncrna.gr <- GRanges(seqnames = Rle(lncrna$chromosome),
                     ranges = IRanges(start = lncrna$start, end = lncrna$end),
                     strand = Rle(lncrna$strand), id = lncrna$ID, seqlengths = yeast.length)

##### Get bam file names
bam.files <- c("Dbp2_1.sorted.twoMatches.bam", "Dbp2_2.sorted.twoMatches.bam", "Dbp2_3.sorted.twoMatches.bam")

##### parameters for scanBam
### Not filtering by mapq >= 30 (there may be multi-mappers)
### Make GRange for scanBamParam
yeast.gr <- GRanges(seqnames = Rle(names(yeast.length)),
                    ranges = IRanges(start = rep(1, length(yeast.length)), end = yeast.length),
                    strand = rep("*", length(yeast.length)), seqlengths = yeast.length)
my.m.fmr.flag <- scanBamFlag(isUnmappedQuery = F)
my.param <- ScanBamParam(which = yeast.gr, flag = my.m.fmr.flag,
                         what = c("flag", "strand", "mapq", "cigar", "qname"))

##### read bam files
# good.mapq <- 30
good.flags <- c(0, 16) 
### If the read is mapped to multiple places, it'll have several alignments,
### but only one 'primary alignment' will be kept using this flag filter

sample.names <- c("D21", "D22", "D23")
### D21
D21.bam.ga <- readGAlignments(bam.files[1], param = my.param) ### 19961599 alignments
## only use the reads that pass through the flag score
D21.good.reads <- which(mcols(D21.bam.ga)$flag %in% good.flags)
D21.bam.ga <- D21.bam.ga[D21.good.reads] ### still the same number of alignments; the filtering is not neccessary (Nadia might have done it) 
### D22
D22.bam.ga <- readGAlignments(bam.files[2], param = my.param) ### 17845843 alignments
### D23
D23.bam.ga <- readGAlignments(bam.files[3], param = my.param) ### 25253224 reads

##### Filter out reads that cover abnormal length on the chromosome (too long)
### Find out the reasonable lengths of transcripts
yn.utr.range <- ranges(yn.utr.gr)
yn.utr.range.df <- as.data.frame(yn.utr.range)
summary(yn.utr.range.df$width) ### Min: 189; Max: 15190

all.rna.range <- ranges(all.rna.gr)
all.rna.range.df <- as.data.frame(all.rna.range)
summary(all.rna.range.df$width) ### Min: 51; Max: 14730

lncrna.range <- ranges(lncrna.gr)
lncrna.range.df <- as.data.frame(lncrna.range)
summary(lncrna.range.df$width) ### Min: 33; Max: 8784

### Filter out reads that cover a region that is too long (>15200 nt)
D21.bam.ga.filt <- D21.bam.ga[which(width(D21.bam.ga)<=15200)] ### 19931680 alignments left
summary(width(D21.bam.ga.filt)) ### Min: 11, Max: 15200
D22.bam.ga.filt <- D22.bam.ga[which(width(D22.bam.ga)<=15200)] ### 17780376 alignments left
summary(width(D22.bam.ga.filt)) ### Min: 11, Max: 15200
D23.bam.ga.filt <- D23.bam.ga[which(width(D23.bam.ga)<=15200)] ### 25187733 alignments left
summary(width(D23.bam.ga.filt)) ### Min: 11, Max: 15160

##### convert bam.ga into gr objects
D21.filt.gr <- GRanges(seqnames = seqnames(D21.bam.ga.filt), ranges = ranges(D21.bam.ga.filt),
                       strand = strand(D21.bam.ga.filt), id = as.character(mcols(D21.bam.ga.filt)$qname),
                       seqlengths = yeast.length)
D22.filt.gr <- GRanges(seqnames = seqnames(D22.bam.ga.filt), ranges = ranges(D22.bam.ga.filt),
                       strand = strand(D22.bam.ga.filt), id = as.character(mcols(D22.bam.ga.filt)$qname),
                       seqlengths = yeast.length)
D23.filt.gr <- GRanges(seqnames = seqnames(D23.bam.ga.filt), ranges = ranges(D23.bam.ga.filt),
                       strand = strand(D23.bam.ga.filt), id = as.character(mcols(D23.bam.ga.filt)$qname),
                       seqlengths = yeast.length)

##### Find the intersection between D21 and D22 samples, and use these ranges to filter D23 reads
D21.D22.inter <- intersect(D21.filt.gr, D22.filt.gr, ignore.strand = F)
### Only keep D23 reads that intersect with D21.D22.inter
D3.inter12.index <- findOverlaps(D21.D22.inter, D23.filt.gr, ignore.strand = F)
D3.inter12.index.df <- data.frame(inter12_index = queryHits(D3.inter12.index),
                                  D3_index = subjectHits(D3.inter12.index))
D3.inter12.index.df$id<- (mcols(D23.filt.gr)$id)[D3.inter12.index.df$D3_index]
D3.inter12.index.range <- ranges(D23.filt.gr)[D3.inter12.index.df$D3_index]
D3.inter12.index.range.df <- as.data.frame(D3.inter12.index.range)
D3.inter12.index.df$read_start <- D3.inter12.index.range.df$start
D3.inter12.index.df$read_end <- D3.inter12.index.range.df$end
D3.inter12.index.df$strand <- (as.character(strand(D23.filt.gr)))[D3.inter12.index.df$D3_index]
D3.inter12.index.df$chr <- (as.character(seqnames(D23.filt.gr)))[D3.inter12.index.df$D3_index]
## re-format the table and don't let a read appear multiple times
D3.inter12.index.df2 <- D3.inter12.index.df[,c(3,7,6,4,5)] ### 24963456 entries
D3.inter12.index.df2 <- unique(D3.inter12.index.df2)
dim(D3.inter12.index.df2) ### 24920578 entries
D3.inter12.reads <- unique(D3.inter12.index.df$id)
length(D3.inter12.reads) ### 24920578 reads kept (originally in D23.filt.gr: 25187733, less than 1% was not kept)
## Output the reads
write.table(D3.inter12.index.df2,
            file = "/scratch/snyder/l/lai64/060618_Dbp2iCLIPdata_reprocessing/Filtered_Dbp2iCLIPreads/060618_Dbp23_iCLIPreads_reproduced.txt",
            quote = F, sep = "\t", row.names = F)

##### Find the intersection between D21 and D23 samples, and use these ranges to filter D22 reads
D21.D23.inter <- intersect(D21.filt.gr, D23.filt.gr, ignore.strand = F)
### Only keep D22 reads that intersect with D21.D23.inter
D2.inter13.index <- findOverlaps(D21.D23.inter, D22.filt.gr, ignore.strand = F)
D2.inter13.index.df <- data.frame(inter13_index = queryHits(D2.inter13.index), D2_index = subjectHits(D2.inter13.index))
D2.inter13.index.df$id<- (mcols(D22.filt.gr)$id)[D2.inter13.index.df$D2_index]
D2.inter13.index.range <- ranges(D22.filt.gr)[D2.inter13.index.df$D2_index]
D2.inter13.index.range.df <- as.data.frame(D2.inter13.index.range)
D2.inter13.index.df$read_start <- D2.inter13.index.range.df$start
D2.inter13.index.df$read_end <- D2.inter13.index.range.df$end
D2.inter13.index.df$strand <- (as.character(strand(D22.filt.gr)))[D2.inter13.index.df$D2_index]
D2.inter13.index.df$chr <- (as.character(seqnames(D22.filt.gr)))[D2.inter13.index.df$D2_index]
## re-format the table and don't let a read appear multiple times
D2.inter13.index.df2 <- D2.inter13.index.df[,c(3,7,6,4,5)] ### 17704641 entries
D2.inter13.index.df2 <- unique(D2.inter13.index.df2)
dim(D2.inter13.index.df2) ### 17663576 entries
D2.inter13.reads <- unique(D2.inter13.index.df$id)
length(D2.inter13.reads) ### 17663576 (originally in D22.filt.gr: 17780376)
## Output the reads
write.table(D2.inter13.index.df2,
            file = "/scratch/snyder/l/lai64/060618_Dbp2iCLIPdata_reprocessing/Filtered_Dbp2iCLIPreads/060618_Dbp22_iCLIPreads_reproduced.txt",
            quote = F, sep = "\t", row.names = F)

##### Find the intersection between D22 and D23 samples, and use these ranges to filter D21 reads
D22.D23.inter <- intersect(D22.filt.gr, D23.filt.gr, ignore.strand = F)
### Only keep D21 reads that intersect with D22.D23.inter
D1.inter23.index <- findOverlaps(D22.D23.inter, D21.filt.gr, ignore.strand = F)
D1.inter23.index.df <- data.frame(inter23_index = queryHits(D1.inter23.index), D1_index = subjectHits(D1.inter23.index))
D1.inter23.index.df$id<- (mcols(D21.filt.gr)$id)[D1.inter23.index.df$D1_index]
D1.inter23.index.range <- ranges(D21.filt.gr)[D1.inter23.index.df$D1_index]
D1.inter23.index.range.df <- as.data.frame(D1.inter23.index.range)
D1.inter23.index.df$read_start <- D1.inter23.index.range.df$start
D1.inter23.index.df$read_end <- D1.inter23.index.range.df$end
D1.inter23.index.df$strand <- (as.character(strand(D21.filt.gr)))[D1.inter23.index.df$D1_index]
D1.inter23.index.df$chr <- (as.character(seqnames(D21.filt.gr)))[D1.inter23.index.df$D1_index]
## re-format the table and don't let a read appear multiple times
D1.inter23.index.df2 <- D1.inter23.index.df[,c(3,7,6,4,5)] ### 19710582 entries
D1.inter23.index.df2 <- unique(D1.inter23.index.df2)
dim(D1.inter23.index.df2) ### 19676903 entries
D1.inter23.reads <- unique(D1.inter23.index.df$id)
length(D1.inter23.reads) ### 19676903 (originally in D21.filt.gr: 19931680)
## Output the reads
write.table(D1.inter23.index.df2,
            file = "/scratch/snyder/l/lai64/060618_Dbp2iCLIPdata_reprocessing/Filtered_Dbp2iCLIPreads/060618_Dbp21_iCLIPreads_reproduced.txt",
            quote = F, sep = "\t", row.names = F)



