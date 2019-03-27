#reactivity calculation

library(seqinr)
library(rmngb)
library(parallel)

file_names = c("110716_1152DMS.cnts", "110716_1152YPD.cnts",
               "110716_1154DMS.cnts", "110716_1154YPD.cnts", "110716_WT2DMS.cnts", "110416_WT2YPD.cnts",
               "110716_WT3DMS.cnts", "110716_WT3YPD.cnts", "110716_WT4DMS.cnts", "110716_WT4YPD.cnts")
conds = c("p1_mt", "m1_mt", "p3_mt", "m3_mt", "p1_wt", "m1_wt", "p2_wt", "m2_wt", "p3_wt", "m3_wt")

seq <- read.fasta("candidates_non_overlapping_seq.fa", forceDNAtolower=F)
closeAllConnections()

counts <- list()
for (RNA in names(seq)) {
  if (RNA == "Q0297") next #See Gene_annotation_analyses.R. Inconsistently annotated UTR.
  counts[[RNA]] <- data.frame(seq= rmAttr(seq[[RNA]], except=NULL))
}

for (i in 1:length(file_names)) {
  assign(conds[i], readLines(file(file_names[i], open="r")))
  closeAllConnections()
}

for (i in 1:length(conds)) {
  curr <- get(conds[i])
  curr_list <- list()
  for (j in 1:(length(curr)/4)) {
    curr_name = curr[4*j - 3]
    if (curr_name == "Q0297") next #See Gene_annotation_analyses.R. Inconsistently annotated UTR.
    counts[[curr_name]][paste0(conds[i], "_counts")] = as.numeric(strsplit(curr[4*j - 2], split="\t")[[1]])
    counts[[curr_name]][paste0(conds[i], "_cov")] = as.numeric(strsplit(curr[4*j - 1], split="\t")[[1]])
  }
}

#Note that the counts information for minus strand transcripts needs to be reversed. 
#-------------
#Get genes start and end positions
genes <- read.table("Yassour_Nagalakshmi_outermost_UTR_part1.txt", header= T)
genes <- genes[which(as.character(genes$id) %in% names(counts)), ]
#----------------
#Invert counts sequence for minus strand transcripts
#extract_counts.py generates them such that the counts corresponding to 5' end 
#are to the end of counts sequence.

counts <- mcmapply(function(x, y) {
  if ( as.character(genes$strand[which(as.character(genes$id) == y)]) == "-"
  ) {
    data.frame(seq= as.character(x[, 1]), apply(x[, -1], 2, rev))
  } else x
}, counts, names(counts), SIMPLIFY = F, mc.cores=parallel::detectCores()-1)

#-------------
#Get reactivities for all RNAs.
rates <- mclapply(counts, function(x) {
  curr <- x[, paste0(conds, "_counts")]/x[, paste0(conds, "_cov")]
  colnames(curr) <- conds
  curr$seq <- x$seq
  states <- paste0(substr(conds, 4, 5)[c(1:5)*2], substr(conds, 2, 2)[c(1:5)*2])
  for (i in states) {
    rep <- substr(i, 3, 3)
    strain <- substr(i, 1, 2)
    curr[, i] <- (curr[, paste0("p", rep, "_", strain)] - curr[, paste0("m", rep, "_", strain)])/(1-curr[, paste0("m", rep, "_", strain)])
    curr[, i][curr[, i] < 0] <- 0
  }
  curr}, mc.cores=parallel::detectCores()-1)


#--------------
#Add annotation information to counts object.
for (i in 1:length(counts)) {
  curr_name <- names(counts)[i]
  
  curr_ann_info <- genes[which(as.character(genes$id) == curr_name), ]
  curr_utr5 <- curr_ann_info$utr5.end - curr_ann_info$utr5.start +1
  curr_cds <- curr_ann_info$orf.end - curr_ann_info$orf.start +1
  curr_utr3 <- curr_ann_info$utr3.end - curr_ann_info$utr3.start +1
  
  curr_unannotated <- nrow(counts[[i]]) - (curr_utr5 + curr_cds + curr_utr3)
  counts[[i]]$ann <- c(rep("5\' UTR", curr_utr5),
                       rep("CDS", curr_cds),
                       rep("3\' UTR", curr_utr3),
                       rep("Unannotated", curr_unannotated))
}

#--------------
#Add intron annotation information to counts object.

intron <- read.table("010417_R64-2-1_intron_chromosome_coordinates.txt", header=T)

for (i in 1:length(counts)) {
  curr_name <- names(counts)[i]
  
  curr_ann_info <- genes[which(as.character(genes$id) == curr_name), ]
  curr_utr5 <- curr_ann_info$utr5.end - curr_ann_info$utr5.start +1
  curr_cds <- curr_ann_info$orf.end - curr_ann_info$orf.start +1
  curr_utr3 <- curr_ann_info$utr3.end - curr_ann_info$utr3.start +1
  
  curr_unannotated <- nrow(counts[[i]]) - (curr_utr5 + curr_cds + curr_utr3)
  counts[[i]]$ann <- c(rep("5\' UTR", curr_utr5),
                       rep("CDS", curr_cds),
                       rep("3\' UTR", curr_utr3),
                       rep("Unannotated", curr_unannotated))
  
  curr_intron <- intron[which(as.character(intron$ID) == curr_name), ]
  if (nrow(curr_intron)) {
    for (intr in 1:nrow(curr_intron)) {
      if (curr_ann_info$strand == "+") {
        curr_intron_start = curr_intron[intr, "start"] - min(curr_ann_info[, -(1:4)])+1
        curr_intron_end = curr_intron[intr, "end"] - min(curr_ann_info[, -(1:4)])+1
      } else {
        curr_intron_start = max(curr_ann_info[, -(1:4)]) - curr_intron[intr, "end"] +1
        curr_intron_end = max(curr_ann_info[, -(1:4)]) - curr_intron[intr, "start"] +1
      }
      
      counts[[i]]$ann[curr_intron_start:curr_intron_end] = "Intron"
    }
  }
  
}


#-----------------------
#2-8% normalization
states <- c(paste0("wt", 1:3), paste0("mt", c(1,3)))
reac_norm <- mcmapply(function(x, y) {
  curr_valid_sites <- as.character(y$ann) %in% c("5\' UTR", "CDS", "3\' UTR") & as.character(y$seq) %in% c("A", "C")
  apply(x[, states], 2, 
        function(z) {
          
          curr_reac <- z[curr_valid_sites]
          
          sorted <- curr_reac[order(curr_reac)]
          if (any(is.na(sorted))) {
            normalize.range <- c(round((min(which(is.na(sorted)))-1) * .9), round((min(which(is.na(sorted)))-1) * .98))
          } else {
            normalize.range <- c(round(length(sorted) * .9), round(length(sorted)* .98))
          }
          normalizer <- mean(sorted[normalize.range[1]:normalize.range[2]])
          
          z/normalizer
        })
}, rates, counts, SIMPLIFY= FALSE, mc.cores=parallel::detectCores()-1)


#-----------------------
#Some coverage statistics
coverages_all <- mcmapply(function(y) {
  plus_cov_columns = grepl("_cov", names(y)) & grepl("p", names(y))
  AC_rows = (as.character(y$seq) %in% c("A", "C")) & (as.character(y$ann) %in% c("3\' UTR", "5\' UTR", "CDS"))
  cov_AC = y[AC_rows, plus_cov_columns]
  mean(as.matrix(cov_AC), na.rm=T)
}, counts, mc.cores=parallel::detectCores()-1)

coverages_utr <- mcmapply(function(y) {
  plus_cov_columns = grepl("_cov", names(y)) & grepl("p", names(y))
  AC_rows = (as.character(y$seq) %in% c("A", "C")) & (as.character(y$ann) %in% c("3\' UTR"))
  cov_AC = y[AC_rows, plus_cov_columns]
  mean(as.matrix(cov_AC), na.rm=T)
}, counts, mc.cores=parallel::detectCores()-1)