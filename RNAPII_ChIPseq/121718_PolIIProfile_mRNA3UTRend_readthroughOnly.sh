
#!/bin/sh -l
#PBS -q biochem
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=4,naccesspolicy=shared
#PBS -N 121718_PolIIProfile_mRNA3UTRend_readthroughOnly

cd $PBS_O_WORKDIR
module load deeptools/3.1.1


computeMatrix reference-point -S 121118_WT_merged_FE.bw 121118_dbp2_merged_FE.bw -R 121718_mRNA3UTR_Readthrough.bed --outFileName 121718_PolII_readthrough_mRNAs_3UTRend_mat.gz --outFileNameMatrix 121718_PolII_readthrough_mRNAs_3UTRend_matForPlot.txt --afterRegionStartLength 500 --referencePoint TES -p max
