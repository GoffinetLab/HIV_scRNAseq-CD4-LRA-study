#!/bin/bash

#SBATCH --job-name=cr_c_hg38

# Set the file to write the stdout and stderr to (if -e is not set; -o or --output).
#SBATCH --output=logs/%x-%j.log

# Set the number of cores (-n or --ntasks).
#SBATCH --ntasks=16

# Set the memory per CPU. Units can be given in T|G|M|K.
#SBATCH --mem-per-cpu=20G

# Set the partition to be used (-p or --partition).
#SBATCH --partition=medium
 
# Set the expected running time of your job (-t or --time).
# Formats are MM:SS, HH:MM:SS, Days-HH, Days-HH:MM, Days-HH:MM:SS
#SBATCH --time=20:00:00


#1 = name for run or name to give sample, e.g. 001, 002 etc.

## For running CR v9.0.0 with GrCj38-2024 only reference (10x pre-built version)

# get sample number and append to P1831_ to complete path
# e.g. CG_JK_ + 001
num=$1
sample="CG_JK_"
sample+=$num
id+=GRCh38-2024_"$num"

echo "----------------------------------------------------------------------------"
echo This is the processing of
#echo $1
echo "Sample number"
echo $1
echo "ID assigned"
echo $id
echo "----------------------------------------------------------------------------"



cellranger count --id=$id \
--fastqs="/data/cephfs-2/unmirrored/projects/goffinet-hiv-cd4-lra/fastq/" \
--transcriptome="/data/cephfs-1/work/groups/goffinet/users/postmusd_c/cellranger/refdata-gex-GRCh38-2024-A/" \
--sample=$sample \
--create-bam false \
  --jobmode=slurm \
  --maxjobs=100 \
  --jobinterval=1000