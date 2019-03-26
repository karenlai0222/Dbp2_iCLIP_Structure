# Genome-wide Discovery of DEAD-box RNA Helicase Targets Reveals RNA Structural Remodeling in Transcription Termination

## Workflow for iCLIP
1. Removing Solexa adaptors using Trimmomatic (v0.36)
2. De-mutiplexing reads (use the script here: https://github.com/qczhang/icSHAPE/blob/master/scripts/splitFastq.pl)
3. Removing PCR duplicates based on built-in random barcodes (use the script here: https://github.com/qczhang/icSHAPE/blob/master/scripts/readCollapse.pl)
4. Trimming the barcode sequence from the 5' end of retained forward reads using cutadapt (v1.9.1)
5. Mapping processed reads to the S288C reference genome (R64-2-1, from Saccharomyces Genome Database) using STAR (v2.5.2b)
6. Reads mapped to one or two sites were kept in the SAM/BAM file (grep 'NH:i:[1-2]' in unix)
7. For each replicate, reads that did not overlap with any read in the other two replicates were discarded (iCLIP_reads_intersection_from_bam.R)
8. The nucleotide position before the start of each read was extracted from the forward reads as the crosslinking site in each replicate (Dbp2iCLIP_sites_from_all_intersected_reads.R)
9. Assigning Dbp2 binding sites to different RNA classes (Dbp2_sites_toRNAclasses.R)
..-Transcripts that had less than 5 counts in each library were filtered from the analysis.
..-Only transcripts that were identified in all three replicates were regarded as binding targets.
..-A pie chart can be made based on the data from this step (Dbp2iCLIP_SitesFromFilteredReads_distribution.R)
10. 
