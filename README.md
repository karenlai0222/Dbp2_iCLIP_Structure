# Genome-wide Discovery of DEAD-box RNA Helicase Targets Reveals RNA Structural Remodeling in Transcription Termination

## Workflow for iCLIP
1. Removing Solexa adaptors using Trimmomatic (v0.36)
2. De-mutiplexing reads (use the script here: https://github.com/qczhang/icSHAPE/blob/master/scripts/splitFastq.pl)
3. Removing PCR duplicates based on built-in random barcodes (use the script here: https://github.com/qczhang/icSHAPE/blob/master/scripts/readCollapse.pl)
4. Trimming the barcode sequence from the 5' end of retained forward reads using cutadapt (v1.9.1)
5. Mapping processed reads to the S288C reference genome (R64-2-1, from Saccharomyces Genome Database) using STAR (v2.5.2b)
6. Reads mapped to one or two sites were kept in the SAM/BAM file. (grep 'NH:i:[1-2]' in unix)
7. 
