#reactivity calculation

library(seqinr)
library(rmngb)
library(parallel)

file_names_snoRNA = c("110716_1152DMS_snoRNA.cnts", "110716_1152YPD_snoRNA.cnts",
               "110716_1154DMS_snoRNA.cnts", "110716_1154YPD_snoRNA.cnts", 
               "110716_WT2DMS_snoRNA.cnts", "110416_WT2YPD_snoRNA.cnts",
               "110716_WT3DMS_snoRNA.cnts", "110716_WT3YPD_snoRNA.cnts",
               "110716_WT4DMS_snoRNA.cnts", "110716_WT4YPD_snoRNA.cnts")
conds = c("p1_mt", "m1_mt", "p3_mt", "m3_mt", "p1_wt", "m1_wt", 
          "p2_wt", "m2_wt", "p3_wt", "m3_wt")

seq_snoRNA <- read.fasta("snoRNA.fa", forceDNAtolower=F)
closeAllConnections()

counts_snoRNA <- list()
for (RNA in names(seq_snoRNA)) {
  counts_snoRNA[[RNA]] <- data.frame(seq= rmAttr(seq_snoRNA[[RNA]], 
                                                 except=NULL))
}

for (i in 1:length(file_names_snoRNA)) {
  assign(conds[i], readLines(file(file_names_snoRNA[i], open="r")))
  closeAllConnections()
}



#To confirm, check YDR363W-A reactivities match previous calculations.
# ydr363wa = data.frame(n = 1:1700)
# YDR363W-A counts from new and old scripts match.
for (i in 1:length(conds)) {
  curr <- get(conds[i])
  curr_list <- list()
  for (j in 1:(length(curr)/4)) {
    curr_name = curr[4*j - 3]
    if (curr_name == "YDR363W-A") {
      # ydr363wa[paste0(conds[i], "_counts")] = as.numeric(strsplit(curr[4*j - 2], split="\t")[[1]])
      # ydr363wa[paste0(conds[i], "_cov")] = as.numeric(strsplit(curr[4*j - 1], split="\t")[[1]])
      next
    }
    counts_snoRNA[[curr_name]][paste0(conds[i], "_counts")] = as.numeric(strsplit(curr[4*j - 2], split="\t")[[1]])
    counts_snoRNA[[curr_name]][paste0(conds[i], "_cov")] = as.numeric(strsplit(curr[4*j - 1], split="\t")[[1]])
  }
}

#Note that the counts information for minus strand transcripts needs to be reversed. 
#-------------
#Get genes start and end positions
genes_snoRNA <- read.table("snoRNA.txt", header= T)
genes_snoRNA <- genes_snoRNA[which(as.character(genes_snoRNA$name) %in% names(counts_snoRNA)), ]
#----------------
#Invert counts sequence for minus strand transcripts
#extract_counts.py generates them such that the counts corresponding to 5' end 
#are to the end of counts sequence.

counts_snoRNA <- mapply(function(x, y) {
  if ( as.character(genes_snoRNA$strand[which(as.character(genes_snoRNA$name) == y)]) == "-"
  ) {
    data.frame(seq= as.character(x[, 1]), apply(x[, -1], 2, rev))
  } else x
}, counts_snoRNA, names(counts_snoRNA), SIMPLIFY = F)

#-------------
#Get reactivities for all RNAs.
rates_snoRNA <- lapply(counts_snoRNA, function(x) {
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
  curr})



#-----------------------
#2-8% normalization
states <- c(paste0("wt", 1:3), paste0("mt", c(1,3)))
reac_norm_snoRNA <- mcmapply(function(x, y) {
  curr_valid_sites <- as.character(y$seq) %in% c("A", "C")
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
}, rates_snoRNA, counts_snoRNA, SIMPLIFY= FALSE)


#-----------------------
#Some coverage statistics
coverages_all_snoRNA <- mapply(function(y) {
  plus_cov_columns = grepl("_cov", names(y)) & grepl("p", names(y))
  AC_rows = (as.character(y$seq) %in% c("A", "C"))
  cov_AC = y[AC_rows, plus_cov_columns]
  mean(as.matrix(cov_AC), na.rm=T)
}, counts_snoRNA)

save.image("snoRNA_reactivity.RData")
