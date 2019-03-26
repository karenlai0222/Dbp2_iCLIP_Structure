
#!/bin/sh -l
#PBS -q biochem
#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=6,naccesspolicy=shared
#PBS -N 121118_MACS_myChIP_foldEnrichment_merged

cd $PBS_O_WORKDIR
module load macs/2.1.2

macs2 bdgcmp -t 121118_WT_merged_treat_pileup.bdg -c 121118_WT_merged_control_lambda.bdg -o 121118_WT_merged_FE.bdg -m FE
macs2 bdgcmp -t 121118_WT_merged_treat_pileup.bdg -c 121118_WT_merged_control_lambda.bdg -o 121118_WT_merged_logLR.bdg -m logLR -p 0.00001

macs2 bdgcmp -t 121118_dbp2_merged_treat_pileup.bdg -c 121118_dbp2_merged_control_lambda.bdg -o 121118_dbp2_merged_FE.bdg -m FE
macs2 bdgcmp -t 121118_dbp2_merged_treat_pileup.bdg -c 121118_dbp2_merged_control_lambda.bdg -o 121118_dbp2_merged_logLR.bdg -m logLR -p 0.00001
