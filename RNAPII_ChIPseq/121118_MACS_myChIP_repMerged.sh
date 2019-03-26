
#!/bin/sh -l
#PBS -q biochem
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=6,naccesspolicy=shared
#PBS -N 121118_MACS_myChIP_repMerged

cd $PBS_O_WORKDIR
module load macs/2.1.2

macs2 callpeak -t 120918_WT1_IP_sorted.bam 120918_WT2_IP_sorted.bam 120918_WT3_IP_sorted.bam -c 120918_WT1_input_sorted.bam 120918_WT2_input_sorted.bam 120918_WT3_input_sorted.bam -n 121118_WT_merged --outdir ../121118_MACS_myChIP_repMerged -f BAMPE -g 1.2e7 -B --SPMR

macs2 callpeak -t 120918_dbp201_IP_sorted.bam 120918_dbp202_IP_sorted.bam 120918_dbp203_IP_sorted.bam -c 120918_dbp201_input_sorted.bam 120918_dbp202_input_sorted.bam 120918_dbp203_input_sorted.bam -n 121118_dbp2_merged --outdir ../121118_MACS_myChIP_repMerged -f BAMPE -g 1.2e7 -B --SPMR



