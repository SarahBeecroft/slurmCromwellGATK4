#!/bin/bash -l
#SBATCH --job-name=gsutil
#SBATCH --account=pawsey0001
#SBATCH --nodes=1
#SBATCH --time=05:00:00
#SBATCH --partition=copyq
#SBATCH --ntasks=16
#SBATCH --export=NONE

conda activate google_cloud

mkdir -p /scratch/pawsey0001/sbeecroft/hg38/

gsutil du -sh gs://genomics-public-data/references/hg38/v0
gsutil -m cp -r gs://genomics-public-data/references/hg38/v0 /scratch/pawsey0001/sbeecroft/hg38/
