
##########
##### This is to get the site distribution for all RNA types
##########

##### Read the summary in
sum.table <- read.delim("~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/060618_Dbp2iCLIPdata_reprocessing/060718_Dbp2iCLIP_SitesFromIntersectedReads_distributionSummaryFiltered.txt",
                        as.is = T)

##### Make the pie chart
### reorder the types of Dbp2-binding RNAs
sum.table.order <- sum.table[c(1,4,7,5,9,3),]
### Plot
setwd("~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/060618_Dbp2iCLIPdata_reprocessing/")
tiff("060718_Dbp2iCLIP_SitesFromIntersectedReads_average_RNAtypes_pie_res600.tiff", height = 4,
     width = 4, units = 'in', res = 600)
par(lwd = 1.5)
# pie(sum.table.order$average_percent, labels = c("mRNA", "snoRNA", "tRNA", "snRNA", "lncRNA", "rRNA"), cex = 0.8)
pie(sum.table.order$average_percent, labels = c("", "", "", "", "", ""), cex = 0.8)
dev.off()

# ##### Make the pie chart
# ##### Types of Dbp2-binding RNAs (from average of 3 replicates)
# # colors = c("red", "yellow", "green", "violet", "orange", "blue", "pink", "cyan") 
# pie(c(44.2,38.6,4.6,3.7,6.9,2.1), labels = c("mRNA", "snoRNA", "lncRNA", "snRNA", "tRNA", "rRNA"), cex = 1.5)
# ### Try to output a high resolution figure
# setwd("~/Desktop/Research_projects_2015-2016/iCLIP/results_from_scripts/051817_Dbp2_iCLIP_reads_counting_after_map_to_RDN37RDN5_using_STAR/")
# tiff("071217_iCLIP_Dbp2_average_RNAtypes_pie_chart_res600.tiff", height = 4, width = 4, units = 'in', res = 600)
# par(lwd = 1)
# pie(c(44.2,38.6,4.6,3.7,6.9,2.1), labels = c("", "", "", "", "", ""))
# dev.off()
