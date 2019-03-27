load("snoRNA_reactivity.RData")

filenames <- c("wt1.reac", "wt2.reac", "wt3.reac", "mt1.reac", "mt3.reac")

for (r in 1:5) {
  sink(filenames[r])
  for (i in 1:length(reac_norm_snoRNA)) {
    curr_rna = names(counts)[i]
    curr_cov = coverages_all[curr_rna]
    curr_reac = reac_norm[[curr_rna]][, strsplit(filenames[r], split= ".reac")[[1]][1]]
    curr_reac = paste0(curr_reac, collapse= '\t')
    
    cat(paste(curr_rna, curr_cov, sep= '\t'))
    cat('\n')
    cat(curr_reac)
    cat('\n\n')
    
  }
  sink()
}
