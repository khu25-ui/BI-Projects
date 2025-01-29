# Prgram Name
Final Project

# Date
April 7, 2023

# Description
In this final project we work with genome assembly from DNA-seq data, Transcriptome assembly from RNA-seq data, Gene annotation and Protein function prediction.
In this final project the organism we are working with is Escherichia coli. The data for this organism was taken from NCBI site(NCBI Sequence Reads ArchiveLinks to an external site..) The filters used to find the data for DNA sequence were, Source: DNA, Library Layout: paired, Platform: Illumina, Strategy: Genome, File Type: fastq.
For RNA sequence the filters used were the same except for the source : RNA. For the DNA sequence we choose one SRR number which has bases more that 100 million.
For the RNA sequence we choose upto two SRR numbers which are around 200 million bases.

The main objectives for the final project are :
1. Perform a preliminary genome assembly from short-read DNA-seq.
2. Assemble a transcriptome from RNA-seq data.
3. Annotate genes within the genome.
4. Predict protein function for those genes.

# Scripts used in the module
Steps in module 6

1. getNGS.sh - In this script we take the SRR number from the NCBI site. The SRR number is SRR24007554. We will use the SRR to retrieve sequence data, then use one of
the bacterial genome assemblers to assemble thereads into a genome. For the getNGS.sh we use fasterq-dump. The command-line utility to retrieve sequence data in FASTQ  format from the SRR is fasterq-dump.

2. trim.sh - In this script we will be performing quality trimming with trimmomatic. The parameters used in this script are:
trimmomatic - Trimmomatic is a Java-based quality trimmer that uses a sliding window to determine where quality scores have dropped below a specified threshold. In addition to trimming based on quality scores, Trimmomatic also removes any adapter sequences from the reads.

PE is the first of the parameters we provide. PE indicates that we have paired-end reads.

threads indicates how many server threads to use for this job.

phred33 indicates the quality encoding method used for the reads. The spaces and backslashes (\) at the ends of lines are critical. The backslash allows commands to span multiple lines, so if you don't put a space before the backslash, it's like having no spaces between parameters. There has to be a space between the parameters.
The next parameters are the left and the right read files. we mention their path in the scripts folder.

The next four parameters are the output files for paired and unpaired output. The unpaired output files contain reads where the mate was deleted due to overall low quality. Since you always have to have the same number of reads in the left file and right file for paired-end reads, if one read is deleted entirely due to low quality, its mate has to be deleted as well. The reads whose mates have been deleted end up in the unpaired files.

HEADCROP indicates the number of bases to remove from the beginning, regardless of quality.

ILLUMINACLIP specifies a file of adapter sequences and the number of mismatches allowed in an adapter match.
LEADING and TRAILING determine the minimum quality for trimming the start and end of reads. Since, we are not dropping any leading or lagging bases, the leading and lagging bases are not specified in this script.

SLIDINGWINDOW indicates the sliding window size and the minimum average quality for the bases in that window.

MINLEN specifies the minimum length for a read to be kept.

3. runSpades.sh - In this script we will use the SPAdes assembler, so run spades.py without any parameters to see the SPAdes help menu. Based on the output of the
spades.py we will write a shell scrip to assemble the genome of the given Organism (Escherechia Coli), using the quality trimmed reads in data/trimmed/Paired files.

4. runQuast.sh - quast.py to run a basic analysis report for your genome and determine the N50 for your assembly.

5. sbatch_assembleGenome.sh - This shell script calls all the other scripts. We use the echo function to print the results of all the scripts. This script is performed in the anaconda environment of the discovery.
```
#!/bin/bash
#SBATCH --partition=short               # choose from debug, express, or short
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
```
Steps used in module 7.
1. AipBuild.sh -  In the script we build an index of the reference genome using the gmap_build command. In this script we provide the file paths that contain the data
about the Miseq of the organism Escherechia Coli. GMAP is used to align RNAseq datasets to genome. GSNAP will use this database to perform the alignment of the RNA-Seq reads. The GMAP database is indexed and optimized to allow much faster alignment than possible if the alignment were directly against the FASTA file. -D indicates the
directory in which to build the database. -d indicates the name of the database.

2. getNGSRNA.sh - In this script we take two SRR number from the NCBI site. The SRR numbers are SRR21849116 SRR21849117 . We will use the SRR to retrieve RNA sequence  data, then use one of the bacterial genome assemblers to assemble thereads into a genome. For the getNGSRNA.sh we use fasterq-dump. The command-line utility to retrieve sequence data in FASTQ  format from the SRR is fasterq-dump. 

3. trimAll.sh - This script is used for quality trimming. Trimmomatic is used as the tool. Trimmomatic is recommended by one of the most popular transcriptome assemblers (Trinity). Trimmomatic is a Java-based quality trimmer that uses a sliding window to determine where quality scores have dropped below a specified threshold. In addition to trimming based on quality scores, Trimmomatic also removes any adapter sequences from the reads. Sometimes during library preparation, extra copies of adapters get attached to the beginning or end of the cDNA fragments. These adapters are what sequencers use to immobilize cDNA fragments on the flow cell, but you don't want them treated as actual sequence data. If included in your sequence data, they confuse assemblers and read aligners.

Steps used in module 8.
1. trinityDeNovo.sh - This script works with de novo assemblies. Trinity requires two comma-separated lists of files for the left and right reads. That is leftReads and rightReads. We can find the output in results/trinity_de_novo.

2. analyzeTrinityDeNovo.sh - This shell script is created to check N50 and other stats for the assembly. This will be the same as analyzeTrinity.sh, but the directory and filename for the output will be results/trinity_de_novo and Trinity.fasta. This script runs TrinityStats.pl with your assembled transcriptome as input.

3. sbatch_trinity.sh. In this script all the scripts are called and results are printed using the echo command. We get the log files and trinity_de_novo. 
``` 
#!/usr/bin/bash
#SBATCH --partition=short               # choose from debug, express, or short
#SBATCH --job-name=trinity
#SBATCH --time=20:00:00                 # the code pieces should run in far less than 1 hour
#SBATCH -N 1                            # nodes requested
#SBATCH -n 4                            # task per node requested
#SBATCH --mem=10Gb
#SBATCH --exclusive
#SBATCH --output="batch-%x-%j.output"   # where to direct standard output; will be batch-jobname-jobID.output
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<parikh.khu>@northeastern.edu # Update to your user name!

# Usage: sbatch sbatch_transcriptome.sh


echo "Starting our analysis $(date)"

echo "Loading our BINF6308 Anaconda environment, which includes Trinity."
module load anaconda3/2021.11
source activate BINF-12-2021
echo "Loading samtools."
module load samtools/1.10


echo "Make directory for log files"
mkdir -p results/logs/


echo "Starting De Novo Assembly $(date)"
echo "Assemble the De Novo Transcriptome $(date)"
bash scripts/trinityDeNovo.sh 1>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-trinityDeNovo.log 2>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-trinityDeNovo.err

echo "Analyze the De Novo Transcriptome $(date)"
bash scripts/analyzeTrinityDeNovo.sh 1>results/$SLURM_JOB_NAME-$SLURM_JOB_ID-trinity_de_novo_stats.txt 2>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-analyzeTrinityDeNovo.err

echo "De Novo Assembly complete $(date)"

echo "Assemblies complete $(date)"
```
Steps in module 9
1. longOrfs_args.sh - This script we use Transdecoder.LongOrfs, they find the longest open reading frames and translates them to amino acid (protein) sequences.

2. blastPep_args.sh. - This script aligns proteins from orf to swissport. We can see from the output there is one hit per search because of our -max_target_seqs 1 flag.

3. pfamScan_args.sh - This script uses hmmscan, uses a Hidden Markov Model (HMM) to find protein domains to guide the prediction process.

4. predictProteins_args.sh - This script uses TransDecoder.Predict, takes in the open reading frames, the BLAST output, and the domain information to refine the protein predictions and produce a protein fasta file (pep). 

5. alignPredicted_args.sh. - In this script we use blastp, aligns the long ORFs to SwissProt to identify similar proteins that could guide the prediction process. 
To get a tabular output we use outfmt.

6. sbatch_transdecoder.sh. - In this script we call all the other scripts and print the result using the echo command. We get output files such as predictedproteins.
```
#!/usr/bin/bash
#SBATCH --partition=short               # choose from debug, express, or short
#SBATCH --job-name=transdecoder
#SBATCH --time=10:00:00
#SBATCH -N 1                            # nodes requested
#SBATCH -n 4                            # task per node requested
#SBATCH --mem=10Gb
#SBATCH --exclusive
#SBATCH --output="batch-%x-%j.output"   # where to direct standard output; will be batch-jobname-jobID.output
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<parikh.khu>@northeastern.edu

# Usage: sbatch sbatch_transdecoder.sh
# Assumes input data is in /home/$USER/AiptasiaRNASeq/data/

echo "Starting our analysis $(date)"
echo

# define key constants
TRANSCRIPTOME=data/trinity_de_novo/Trinity.fasta
SWISSPROT_DB=/work/courses/BINF6308/inputFiles/blastDB/swissprot
TRANSDECODER_DIR=results/trinity_de_novo.transdecoder_dir
LONGEST_ORFS=$TRANSDECODER_DIR/longest_orfs.pep
OUTFMT=results/blastPep_args.outfmt6
DOMTBLOUT=results/pfam.domtblout
PFAMA_PATH=/work/courses/BINF6308/inputFiles/SampleDataFiles/Pfam-A.hmm
PREDICTED_PROTEIN_PATH=results/predictedProteins
FINAL_PROTEINS=$PREDICTED_PROTEIN_PATH/*transdecoder.pep

# record these key constants to our batch*.output file by echoing them:
echo "Key parameters"
echo "TRANSCRIPTOME: $TRANSCRIPTOME"
echo "SWISSPROT_DB: $SWISSPROT_DB"
echo "TRANSDECODER_DIR: $TRANSDECODER_DIR"
echo "LONGEST_ORFS: $LONGEST_ORFS"
echo "OUTFMT: $OUTFMT"
echo "DOMTBLOUT: $DOMTBLOUT"
echo "PFAMA_PATH: $PFAMA_PATH"
echo "PREDICTED_PROTEIN_PATH: $PREDICTED_PROTEIN_PATH"
echo "FINAL_PROTEINS: $FINAL_PROTEINS"
echo
echo
echo "Loading our BINF6308 Anaconda environment."
module load anaconda3/2021.11
source activate BINF-12-2021

echo "Copy the data for trinity_de_novo"
cp -r /scratch/parikh.khu/fp_f2022_binf6308_sec3-khu03parikh/results/trinity_de_novo/ /scratch/parikh.khu/fp_f2022_binf6308_sec3-khu03parikh/data/

echo "Starting ORF prediction pipeline $(date)"
echo "Identify longORFs with TransDecoder.LongOrfs on $TRANSCRIPTOME $(date)"
bash scripts/longOrfs_args.sh $TRANSCRIPTOME $TRANSDECODER_DIR \
  1>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-longOrfs_args.log \
  2>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-longOrfs_args.err

echo "BLASTp of longest_orfs.pep against SwissProt BLAST DB at $SWISSPROT_DB $(date)"
bash scripts/blastPep_args.sh $LONGEST_ORFS $SWISSPROT_DB \
  1>$OUTFMT \
  2>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-blastPep_args.err

echo "Create pfamScan with hmmscan using the Pfam-A.hmm file found $PFAMA_PATH $(date)"
bash scripts/pfamScan_args.sh $DOMTBLOUT $PFAMA_PATH $LONGEST_ORFS \
  1>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-pfamScan_args.log \
  2>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-pfamScan_args.err
echo "Predict protiens with TransDecoder.Predict $(date)"
bash scripts/predictProteins_args.sh $TRANSCRIPTOME $TRANSDECODER_DIR $DOMTBLOUT $OUTFMT \
  1>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-predictProteins_args.log \
  2>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-predictProteins_args.err
echo "Copy TransDecoder.Predict outputs to $PREDICTED_PROTEIN_PATH"
mkdir -p $PREDICTED_PROTEIN_PATH
mv *transdecoder* $PREDICTED_PROTEIN_PATH

echo "Align predicted proteins to SwissProt DB $(date)"
bash scripts/alignPredicted_args.sh $FINAL_PROTEINS $SWISSPROT_DB \
  1>results/alignPredicted.txt \
  2>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-alignPredicted_args.err

echo "ORF prediction pipeline complete $(date)"



echo "Analysis complete $(date)"

```
# How to call the program
1. We run the sbatch_assembleGenome.sh script to get all the results. The program is run using ``` sbatch sbatch_assembleGenome.sh```.  It takes around 28 minutes to 
run this proram and get the results.

2. We run the sbatch_trinity.sh, using ```sbatch sbatch_trinity.sh```. It take around 8 hours to run this code.

3. We run the sbatch_transdecoder.sh using ```sbatch sbatch_transdecoder.sh```. It takes around 4 hours to run this code.

# Results/Output
1. sbatch_assembleGenome.sh
![Screenshot (161)](https://user-images.githubusercontent.com/99061467/232373279-17d33662-42b1-4572-a794-52b6f328aa98.png)


2. sbatch_trinity.sh
![Screenshot (162)](https://user-images.githubusercontent.com/99061467/232373467-b86242f9-6d4f-42b9-b6a9-0886a4e845d0.png)


3. sbatch_transdecoder.sh
![Screenshot (163)](https://user-images.githubusercontent.com/99061467/232373505-7bae21f6-9ea7-4db3-89fd-68327f28d66b.png)
