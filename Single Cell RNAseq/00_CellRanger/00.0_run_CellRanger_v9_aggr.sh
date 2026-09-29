#!/bin/bash

#SBATCH --job-name=cr_aggr

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


# change to new directory with aggregated data
cd "/data/cephfs-2/unmirrored/projects/goffinet-hiv-cd4-lra/cd4lra.cr/aggr/"


cellranger aggr --id=aggHIVpos --csv=241211_CRaggr_csv.csv