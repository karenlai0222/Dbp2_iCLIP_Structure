# Genome-wide Discovery of DEAD-box RNA Helicase Targets Reveals RNA Structural Remodeling in Transcription Termination

### The raw data used in these analysis can be found at Gene Expression Omnibus (GEO): GSE106479
### For detailed sequence information, please check the associated publication:
Genome-Wide Discovery of DEAD-Box RNA Helicase Targets Reveals RNA Structural Remodeling in Transcription Termination.
Yu-Hsuan Lai, Krishna Choudhary, Sara C. Cloutier, Zheng Xing, Sharon Aviran* and Elizabeth J. Tran*.
GENETICS 2019; https://doi.org/10.1534/genetics.119.302058

## Workflow for iCLIP-seq data analysis
1. Removing Solexa adaptors using Trimmomatic (v0.36)
2. De-mutiplexing reads (use the script here: https://github.com/qczhang/icSHAPE/blob/master/scripts/splitFastq.pl)
3. Removing PCR duplicates based on built-in random barcodes (use the script here: https://github.com/qczhang/icSHAPE/blob/master/scripts/readCollapse.pl)
4. Trimming the barcode sequence from the 5' end of retained forward reads using cutadapt (v1.9.1)
5. Mapping processed reads to the S288C reference genome (R64-2-1, from Saccharomyces Genome Database) using STAR (v2.5.2b)
6. Reads mapped to one or two sites were kept in the SAM/BAM file (grep 'NH:i:[1-2]' in unix)
7. For each replicate, reads that did not overlap with any read in the other two replicates were discarded (iCLIP_reads_intersection_from_bam.R)
8. The nucleotide position before the start of each read was extracted from the forward reads as the crosslinking site in each replicate (Dbp2iCLIP_sites_from_all_intersected_reads.R)
9. Assigning Dbp2 binding sites to different RNA classes (Dbp2_sites_toRNAclasses.R)  
   - Transcripts that had less than 5 counts in each library were filtered from the analysis.  
   - Only transcripts that were identified in all three replicates were regarded as binding targets.  
   - A pie chart can be made based on the data from this step (Dbp2iCLIP_SitesFromFilteredReads_distribution.R)
10. Meta-analysis of Dbp2 binding sites  
    - Doing three replicates separately (iCLIP_D21_meta_mRNAsBound_FilteredSites.R; iCLIP_D22_meta_mRNAsBound_FilteredSites.R; iCLIP_D23_meta_mRNAsBound_FilteredSites.R)  
    - Combining three replicates after normalization (Dbp2iCLIP_meta_mRNAboundByDbp2_filteredSites_normalized_toRPKM_100bin.R)  

## Workflow for Structure-seq data analysis
1. Removing adaptor sequences using Trimmomatic (v0.36)
2. Trimming random trimers from the 5' end of forward reads using cutadapt (v1.9.1)
3. Mapping processed reads to the S288C reference genome (R64-2-1, from Saccharomyces Genome Database) using STAR (v2.5.2b)  
   - Only uniquely mapped reads (MAPQ = 255 after STAR alignment) were kept for the subsequent analysis  
   - Protein-coding genes (mRNAs) overlapping with at least one other gene on the same strand are not included in the following analysis (Gene_annotation_analysis.R)
4. Get detection counts for each nucleotide  
   - mRNAs: extract_counts_mRNAs.py
   - snoRNAs: extract_counts_snoRNAs.py
5. Calaulate reactivities  
   - mRNAs: Reactivity_analysis_mRNAs.R; generate_reactivity_files_mRNAs.R
   - snoRNAs: Reactivity_analysis_snoRNAs.R; generate_reactivity_files_snoRNAs.R	

## Workflow for RNAPII ChIP-seq data analysis
1. Removing adaptor sequences using Trimmomatic (v0.36)
2. Mapping reads to the S288C reference genome (R64-2-1, from Saccharomyces Genome Database) using Bowtie 2 (v2.3.3.1)
3. Determining peaks using MACS2 (v2.1.2)  
   - Using input samples as a control to find peaks (121118_MACS_myChIP_repMerged.sh)  
   - Deriving the fold enrichment as the noramlized signals shown in figures (121118_MACS_myChIP_foldEnrichment_merged.sh)
4. Analyzing the overall RNAPII occupancy around termination sites of snoRNAs or 3' ends of mRNAs using deepTools (v3.1.1)  
   - snoRNAs (121418_PolIIProfile_49snoRNAtermStart.sh)  
   - mRNAs (121718_PolIIProfile_mRNA3UTRend.sh; 121718_PolIIProfile_mRNA3UTRend_readthroughOnly.sh)
   - Making figures in R (121418_PolII_profile_snoRNAends.R; 121418_PolII_profile_mRNA3UTR.R)
   
