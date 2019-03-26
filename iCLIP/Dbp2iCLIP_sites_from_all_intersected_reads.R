
##########
##### This is to get the crosslinking sites from all the intersected reads in all replicates
##########
##### This is run on cluster

##### Get read information
### Dbp2-1
D1.over.reads <- read.delim("060618_Dbp21_iCLIPreads_reproduced.txt", quote = "", as.is = T)
## Get crosslinking sites
D1.reads.watson <- D1.over.reads[D1.over.reads$strand=="+", ]
D1.reads.watson.sites <- unlist(lapply(D1.reads.watson$read_start, function(x){x-1}))
D1.reads.watson$D1_xlinking_site <- D1.reads.watson.sites
D1.reads.crick <- D1.over.reads[D1.over.reads$strand=="-", ]
D1.reads.crick.sites <- unlist(lapply(D1.reads.crick$read_end, function(x){x+1}))
D1.reads.crick$D1_xlinking_site <- D1.reads.crick.sites
D1.over.reads.sites <- rbind(D1.reads.watson, D1.reads.crick)
## Output the table
write.table(D1.over.reads.sites,
            file = "060618_iCLIP_Dbp21_intersectedReads_allxlinkingSites.txt",
            quote = F, sep = "\t", row.names = F)

### Dbp2-2
D2.over.reads <- read.delim("060618_Dbp22_iCLIPreads_reproduced.txt", quote = "", as.is = T)
## Get crosslinking sites
D2.reads.watson <- D2.over.reads[D2.over.reads$strand=="+", ]
D2.reads.watson.sites <- unlist(lapply(D2.reads.watson$read_start, function(x){x-1}))
D2.reads.watson$D2_xlinking_site <- D2.reads.watson.sites
D2.reads.crick <- D2.over.reads[D2.over.reads$strand=="-", ]
D2.reads.crick.sites <- unlist(lapply(D2.reads.crick$read_end, function(x){x+1}))
D2.reads.crick$D2_xlinking_site <- D2.reads.crick.sites
D2.over.reads.sites <- rbind(D2.reads.watson, D2.reads.crick)
## Output the table
write.table(D2.over.reads.sites,
            file = "060618_iCLIP_Dbp22_intersectedReads_allxlinkingSites.txt",
            quote = F, sep = "\t", row.names = F)

### Dbp2-3
D3.over.reads <- read.delim("060618_Dbp23_iCLIPreads_reproduced.txt", quote = "", as.is = T)
## Get crosslinking sites
D3.reads.watson <- D3.over.reads[D3.over.reads$strand=="+", ]
D3.reads.watson.sites <- unlist(lapply(D3.reads.watson$read_start, function(x){x-1}))
D3.reads.watson$D3_xlinking_site <- D3.reads.watson.sites
D3.reads.crick <- D3.over.reads[D3.over.reads$strand=="-", ]
D3.reads.crick.sites <- unlist(lapply(D3.reads.crick$read_end, function(x){x+1}))
D3.reads.crick$D3_xlinking_site <- D3.reads.crick.sites
D3.over.reads.sites <- rbind(D3.reads.watson, D3.reads.crick)
## Output the table
write.table(D3.over.reads.sites,
            file = "060618_iCLIP_Dbp23_intersectedReads_allxlinkingSites.txt",
            quote = F, sep = "\t", row.names = F)

##########
##### This part is to get the site counts for each replicate
##########

##### Get the "identification" and the frequency
### Dbp2-1
D1.chr <- D1.over.reads.sites$chr
D1.site <- as.character(D1.over.reads.sites$D1_xlinking_site)
D1.strand <- D1.over.reads.sites$strand
D1.position <- as.factor(paste(D1.chr, D1.site, D1.strand, sep = " "))
D1.position.tb <- as.data.frame(table(D1.position))
## Re-format the dataframe of positions and counts
D1.position.list <- as.character(D1.position.tb$D1.position)
D1.position.sp <- strsplit(D1.position.list, split = " ")
D1.pos.chr <- unlist(lapply(1:length(D1.position.sp), function(x){D1.position.sp[[x]][1]}))
D1.pos.site <- unlist(lapply(1:length(D1.position.sp), function(x){D1.position.sp[[x]][2]}))
D1.pos.strand <- unlist(lapply(1:length(D1.position.sp), function(x){D1.position.sp[[x]][3]}))
D1.pos.count.df <- data.frame(chr = D1.pos.chr, site = D1.pos.site, strand = D1.pos.strand, count = D1.position.tb$Freq)
dim(D1.pos.count.df) ### 136836 unique sites
## Output the table
write.table(D1.pos.count.df,
            file = "060618_iCLIP_Dbp21_intersectedReads_allxlinkingSiteCounts.txt",
            quote = F, sep = "\t", row.names = F)

### Dbp2-2
D2.chr <- D2.over.reads.sites$chr
D2.site <- as.character(D2.over.reads.sites$D2_xlinking_site)
D2.strand <- D2.over.reads.sites$strand
D2.position <- as.factor(paste(D2.chr, D2.site, D2.strand, sep = " "))
D2.position.tb <- as.data.frame(table(D2.position))
## Re-format the dataframe of positions and counts
D2.position.list <- as.character(D2.position.tb$D2.position)
D2.position.sp <- strsplit(D2.position.list, split = " ")
D2.pos.chr <- unlist(lapply(1:length(D2.position.sp), function(x){D2.position.sp[[x]][1]}))
D2.pos.site <- unlist(lapply(1:length(D2.position.sp), function(x){D2.position.sp[[x]][2]}))
D2.pos.strand <- unlist(lapply(1:length(D2.position.sp), function(x){D2.position.sp[[x]][3]}))
D2.pos.count.df <- data.frame(chr = D2.pos.chr, site = D2.pos.site, strand = D2.pos.strand, count = D2.position.tb$Freq)
dim(D2.pos.count.df) ### 59836 unique sites
## Output the table
write.table(D2.pos.count.df,
            file = "060618_iCLIP_Dbp22_intersectedReads_allxlinkingSiteCounts.txt",
            quote = F, sep = "\t", row.names = F)

### Dbp2-3
D3.chr <- D3.over.reads.sites$chr
D3.site <- as.character(D3.over.reads.sites$D3_xlinking_site)
D3.strand <- D3.over.reads.sites$strand
D3.position <- as.factor(paste(D3.chr, D3.site, D3.strand, sep = " "))
D3.position.tb <- as.data.frame(table(D3.position))
## Re-format the dataframe of positions and counts
D3.position.list <- as.character(D3.position.tb$D3.position)
D3.position.sp <- strsplit(D3.position.list, split = " ")
D3.pos.chr <- unlist(lapply(1:length(D3.position.sp), function(x){D3.position.sp[[x]][1]}))
D3.pos.site <- unlist(lapply(1:length(D3.position.sp), function(x){D3.position.sp[[x]][2]}))
D3.pos.strand <- unlist(lapply(1:length(D3.position.sp), function(x){D3.position.sp[[x]][3]}))
D3.pos.count.df <- data.frame(chr = D3.pos.chr, site = D3.pos.site, strand = D3.pos.strand, count = D3.position.tb$Freq)
dim(D3.pos.count.df) ### 109168 unique sites
## Output the table
write.table(D3.pos.count.df,
            file = "060618_iCLIP_Dbp23_intersectedReads_allxlinkingSiteCounts.txt",
            quote = F, sep = "\t", row.names = F)

##########
##### This is to merge binding sites from the three replicates
##########

##### Get the information about unique sites in each replicates
D1.unique.sites <- D1.pos.count.df[,c(1,2,3)]
D2.unique.sites <- D2.pos.count.df[,c(1,2,3)]
D3.unique.sites <- D3.pos.count.df[,c(1,2,3)]

##### combine all unique sites
all.unique.sites <- rbind(D1.unique.sites, D2.unique.sites)
all.unique.sites <- rbind(all.unique.sites, D3.unique.sites)
dim(all.unique.sites) ### 256681
### Don't let a site appear multiple times
all.unique.sites2 <- unique(all.unique.sites)
dim(all.unique.sites2) ### 233264 sites in total
### Order the sites by positions
all.unique.sites2.order <- all.unique.sites2[order(all.unique.sites2$chr, all.unique.sites2$site), ]
### Output the table
write.table(all.unique.sites2,
            file = "060618_Dbp2iCLIP_intersectedReads_allxlinkingSites_combined_unique.txt",
            quote = F, row.names = F)


