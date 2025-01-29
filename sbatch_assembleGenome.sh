#!/bin/bash
#BATCH --partition=short               # choose from debug, express, or short
#SBATCH --job-name=assembleGenome
#SBATCH --time=04:00:00                 # the code pieces should run in far less than 4 hours
#SBATCH -N 1                            # nodes requested
#SBATCH -n 1                            # task per node requested
#SBATCH --output="batch-%x-%j.output"   # where to direct standard output; will be batch-jobname-jobID.output

echo "Starting our analysis $(date)"

ORGANISM="Escherichia coli"  # in future, we will define this as part of a config file
SRR_ID=SRR24007554  # in future, we will define this as part of a config file

# Create the results folder
mkdir -p results/
mkdir -p results/logs

echo "$ORGANISM SRR reads to process: $SRR_ID"

echo "Loading our BINF6308 Anaconda environment."
module load anaconda3/2021.11
source activate BINF-12-2021

echo "Downloading $SRR_ID reads $(date)"
bash /scratch/parikh.khu/fp_f2022_binf6308_sec3-khu03parikh/scripts/getNGS.sh 1>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-getNGS.log 2>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-getNGS.err

echo "Trimming $SRR_ID reads $(date)"
bash /scratch/parikh.khu/fp_f2022_binf6308_sec3-khu03parikh/scripts/trim.sh 1>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-trim.log 2>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-trim.err

echo "Assembling genome from trimmed $SRR_ID reads $(date)"
bash /scratch/parikh.khu/fp_f2022_binf6308_sec3-khu03parikh/scripts/runSpades.sh 1>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-runSpades.log 2>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-runSpades.err

echo "Analyzing genome assembly $(date)"
bash /scratch/parikh.khu/fp_f2022_binf6308_sec3-khu03parikh/scripts/runQuast.sh 1>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-runQuast.log 2>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-runQuast.err

echo "Assembly and analysis complete $(date)"
