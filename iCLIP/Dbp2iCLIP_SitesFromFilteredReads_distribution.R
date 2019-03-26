
##########
##### This is to get the site distribution for all RNA types (Making a pie chart)
##########

##### Read the summary in
sum.table <- read.delim("060718_Dbp2iCLIP_SitesFromIntersectedReads_distributionSummaryFiltered.txt",
                        as.is = T)

##### Make the pie chart
### reorder the types of Dbp2-binding RNAs
sum.table.order <- sum.table[c(1,4,7,5,9,3),]
### Plot
tiff("060718_Dbp2iCLIP_SitesFromIntersectedReads_average_RNAtypes_pie_res600.tiff", height = 4,
     width = 4, units = 'in', res = 600)
par(lwd = 1.5)
# pie(sum.table.order$average_percent, labels = c("mRNA", "snoRNA", "tRNA", "snRNA", "lncRNA", "rRNA"), cex = 0.8)
pie(sum.table.order$average_percent, labels = c("", "", "", "", "", ""), cex = 0.8)
dev.off()
