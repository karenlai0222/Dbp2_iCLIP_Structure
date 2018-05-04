#Find non-overlapping transcripts

library(seqinr)
library(Biostrings)
library(rmngb)


#-------------------------
#Analyze all annotated genes
genomic_seq <- read.fasta("S288C_reference_sequence_R64-2-1_20150113.fsa", forceDNAtolower=F)
genomic_names <- c()
for (i in 1:length(genomic_seq)) {
  ch = strsplit(strsplit(attr(genomic_seq[[i]], "Annot"),split=" ")[[1]][6],split="=")[[1]][2]
  genomic_names <- c(genomic_names, substr(ch, 1, nchar(ch)-1))
}
genomic_names <- replace(genomic_names, which(genomic_names=="mitochondrion"), "M")

names = c("I", "II", "III", "IV", "V", "VI", "VII",
          "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV",
          "XVI", "M")
strand = c("+", "-")

chr_len <- c(230218, 813184, 316620, 1531933, 576874, 270161, 1090940, 562643, 439888, 745751, 666816,
             1078177, 924431, 784333, 1091291, 948066, 85779)

all_strands_names = apply(expand.grid(names, strand), 1, paste0, collapse="")

genes <- read.table("Yassour_Nagalakshmi_outermost_UTR_part1.txt", header= T)

#     chromosome type      id strand orf.start orf.end utr3.start utr3.end utr5.start utr5.end
#6204      chrXV gene YOR245C      -    794076  795332     794076   794075     795333   795631
genes[6204, 8] =  genes[6204, 7]

all_strands <- list()

for (i in 1:34) {
  all_strands[[all_strands_names[i]]] <- data.frame(chromosome=c(), type=c(), id=c(), strand=c(), orf.start=c(),
                                                    orf.end=c(), utr3.start= c(), utr3.end=c(), utr5.start=c(),
                                                    utr5.end= c())
}

for (i in 1: nrow(genes)) {
  curr <- genes[i, ]
  curr_strand <- paste0(strsplit(as.character(curr$chromosome), split="chr")[[1]][2], curr$strand)
  all_strands[[curr_strand]] <- rbind(all_strands[[curr_strand]], curr)
}
#No genes on - strand of chrM

all_strands <- all_strands[1:33]

for (i in 1:33) {
  all_strands[[i]] <- cbind(all_strands[[i]], data.frame(r_start= pmin(all_strands[[i]][, 5], all_strands[[i]][, 6], 
                                                                       all_strands[[i]][, 7], all_strands[[i]][, 8], 
                                                                       all_strands[[i]][, 9], all_strands[[i]][, 10]),
                                                         r_end= pmax(all_strands[[i]][, 5], all_strands[[i]][, 6], 
                                                                     all_strands[[i]][, 7], all_strands[[i]][, 8], 
                                                                     all_strands[[i]][, 9], all_strands[[i]][, 10])))
}

all_strands <- lapply(all_strands, function(x) if (nrow(x) != 0) x[order(x$r_start), ])

isolation_threshold  <- 100
for (i in 1:33) {
  curr <- all_strands[[i]]
  n = nrow(curr)
  overlap =  logical(n)
  isolated = logical(n)
  dist_from_previous = numeric(n)
  dist_from_next = numeric(n)
  
  if (curr$r_end[1] >= curr$r_start[2]) overlap[1] <- TRUE
  if (curr$r_start[2] > curr$r_end[1] + isolation_threshold) isolated[1] <- TRUE
  dist_from_previous[1] =NA
  dist_from_next[1] = curr$r_start[2] - curr$r_end[1]
  
  for (j in 2:(n-1)) {
    if (curr$r_end[j] >= curr$r_start[j+1] | curr$r_start[j] <= curr$r_end[j-1]) overlap[j] <- TRUE
    if (curr$r_start[j+1] > curr$r_end[j] + isolation_threshold & 
        curr$r_end[j-1] < curr$r_start[j] - isolation_threshold) isolated[j] <- TRUE
    dist_from_previous[j] = curr$r_start[j] - curr$r_end[j-1]
    dist_from_next[j] = curr$r_start[j+1] - curr$r_end[j]
  }
  
  if (curr$r_start[n] <= curr$r_end[n-1]) overlap[n] <- TRUE
  if (curr$r_end[n-1] < curr$r_start[n] - isolation_threshold) isolated[n] <- TRUE
  dist_from_previous[n] = curr$r_start[n] - curr$r_end[n-1]
  dist_from_next[n] = NA
  
  overlap_with = numeric(n)
  
  j=0
  olap= F
  while (j < n) {
    j =j+1
    curr_index_start = j
    curr_overlap=0
    
    if (j!=n) {
      if (dist_from_next[j] <0) {
        curr_overlap= curr_overlap +1
        olap= T
      } 
    } else olap=F
    
    while (olap== T) {
      j = j+1
      
      if (j!=n) {
        if (dist_from_previous[j] <0 & dist_from_next[j] <0) {
          curr_overlap= curr_overlap +1
          olap= T
        } else if (dist_from_previous[j] <0 & dist_from_next[j] >=0) {
          curr_overlap= curr_overlap +1
          olap= F
        } 
      } else {
        if (dist_from_previous[j] <0) {
          curr_overlap= curr_overlap +1
          olap= F
        }
      }
      
      overlap_with[curr_index_start:j] <- curr_overlap
    }
    
  }
  
  all_strands[[i]] <- cbind(curr, data.frame(overlap= overlap, overlap_with = overlap_with, isolated= isolated,
                                             d_from_pre= dist_from_previous, d_from_next= dist_from_next,
                                             gene_length= curr$r_end- curr$r_start))
}



lapply(all_strands, function(x) print(c(sum(x$overlap), nrow(x))))
lapply(all_strands, function(x) max(x$overlap_with))

sum(mapply(function(x) nrow(x) - sum(x$overlap), all_strands))

#---------------
#Find coordinates of non-overlapping genes

gene_non_overlapping <- mapply(function(x) {
  x[which(!x$overlap), c("id", "overlap", "isolated")]
  }, all_strands, SIMPLIFY = F)

gene_non_overlapping = do.call("rbind", gene_non_overlapping)

strand <- character(nrow(gene_non_overlapping))
chr <- character(nrow(gene_non_overlapping))
gene_start <- numeric(nrow(gene_non_overlapping))
gene_end <- numeric(nrow(gene_non_overlapping))
for (i in 1:nrow(gene_non_overlapping)) {
  curr <- genes[which(as.character(genes$id) == as.character(gene_non_overlapping$id[i])), ]
  chr[i] <- strsplit(as.character(curr$chromosome), split="chr")[[1]][2]
  strand[i] <- as.character(curr$strand)
  curr_strand <- all_strands[[paste0(chr[i], strand[i])]]
  curr_row <- which(as.character(curr_strand$id) == as.character(gene_non_overlapping$id[i]))
  if (strand[i] == "-") {
    gene_start[i] <- max(curr_strand[curr_row ,"r_start"] - min(curr_strand[curr_row, "d_from_pre"]-1,
                                                                1000, na.rm = T), 1)
    gene_end[i] <- curr_strand[curr_row, "r_end"] 
  } else {
    gene_start[i] <- curr_strand[curr_row ,"r_start"]
    gene_end[i] <- min(curr_strand[curr_row, "r_end"] + min(curr_strand[curr_row, "d_from_next"]-1,
                                                            1000, na.rm = T), chr_len[which(names== chr[i])])
  }
  
}

gene_non_overlapping <- cbind(gene_non_overlapping, data.frame(chr= chr, strand= strand, start= gene_start, end= gene_end))
#Deleting this line: chromosome type      id strand orf.start orf.end utr3.start utr3.end utr5.start utr5.end r_start r_end
# 973      chrVI gene YFL068W      +        53     535        536      670        -82       52     -82   670
# overlap overlap_with isolated d_from_pre d_from_next gene_length
# 973   FALSE            0    FALSE         NA          31         752

#For following line, UTR annotation goes beyond chrM length. Retaining the gene for now but may delete it later.
# chromosome type    id strand orf.start orf.end utr3.start utr3.end utr5.start utr5.end
# 28       chrM gene Q0297      +     85554   85709      85710    85844      85419    85553

gene_non_overlapping <- gene_non_overlapping[-which(as.character(gene_non_overlapping$id) == "YFL068W"), ]

for (i in names(all_strands)) {
  curr_chr <- substr(i, 1, nchar(i)-1)
  if (curr_chr == "Mito") curr_chr = "M"
  curr_strand <- substr(i, nchar(i), nchar(i))
  curr_header <- paste0("chr", curr_chr, "_", if (curr_strand == "+") "sense" else "antisense")
  write(curr_header, file="candidates_non_overlapping.txt", append=T, sep="\t")
  write("\n", file="candidates_non_overlapping.txt", append=T, sep="\t")
  curr_table <- gene_non_overlapping[which(as.character(gene_non_overlapping$chr) == curr_chr & 
                                             as.character(gene_non_overlapping$strand) == curr_strand), c(1, 6,7)]
  for (RNA in 1:nrow(curr_table)) {
    DNAseq = DNAString(paste0(rmAttr(genomic_seq[[which(genomic_names==curr_chr)]], 
                                     except=NULL)[curr_table[RNA, 2]:curr_table[RNA, 3]], collapse=""))
    if (curr_strand == "+") {
      RNAseq <- RNAString(DNAseq)
    } else {
      RNAseq <- paste0(RNAString(reverseComplement(DNAseq)))
    }
    
    write(paste0(">", curr_table[RNA, 1]), file= "candidates_non_overlapping_seq.fa", append=T, sep="\t")
    write(paste0(RNAseq, "\n"), file= "candidates_non_overlapping_seq.fa", append=T, sep="\t")
  }
  
  write.table(curr_table, file="candidates_non_overlapping.txt", append= T, 
              quote=F, sep="\t", row.names=F, col.names= F)
  write("\n", file="candidates_non_overlapping.txt", append=T, sep="\t")
}

#--------------
#Dbp2_binding_transcripts

Dbp2_binding = read.table("011917_iCLIP_Dbp2_binding_mRNA_transcript_list_after_removing_Set1_targets.txt")$V1
#---------------
#Bound in regions of interest

splicing_candidates = read.table("013017_Dbp2_binding_regions_in_target_mRNA_3SpliceSite.txt", 
                                 header= T)$Transcript_id
splicing_candidates = as.character(unique(splicing_candidates))

#TSS = Translation start site
TSS_candidates = read.table("013017_Dbp2_binding_regions_in_target_mRNA_5TranslationalStart.txt",
                            header= T)$Transcript_id
TSS_candidates = as.character(unique(TSS_candidates))

UTR_candidates = read.table("012717_Dbp2_binding_regions_in_target_mRNA_3utr.txt", header= T)$Transcript_id
UTR_candidates = unique(as.character(UTR_candidates))

Dbp2_binding_interesting = Reduce(union, list(splicing_candidates, TSS_candidates, UTR_candidates))

#---------------
#Genes with defects

rt_defect = read.table("Dbp2_bound_readthrough_transcripts.txt")$V1
rt_defect = as.character(unique(rt_defect))

splicing_defect = read.table("013117_splicing_defect_analysis_IAI1.5_genelist_revised.txt")$V1
splicing_defect = as.character(splicing_defect)

all_defective = union(rt_defect, splicing_defect)

sum(all_defective %in% Dbp2_binding)
#---------------
#Intron containing genes
intron_containing = read.table("010417_R64-2-1_intron_chromosome_coordinates.txt", header=T)$ID
intron_containing = unique(as.character(intron_containing))

#---------------
sum(Dbp2_binding %in% as.character(gene_non_overlapping$id))
length(as.character(gene_non_overlapping$id))
#--------------
#All genes with any alternative evidence of structural alteration
candidates = union(Dbp2_binding, all_defective)