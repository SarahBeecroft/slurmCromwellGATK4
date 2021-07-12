#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH -c 2
#SBATCH --job-name=cromwell2
#SBATCH --partition=longq
#SBATCH --account=pawsey0001
#SBATCH --mem=8000
#SBATCH --time=96:00:00
#SBATCH --export=NONE

module load java
conda activate gatk4_pipeline

java -Dconfig.file=/scratch/pawsey0001/sbeecroft/cromwell2/slurm.conf -jar /scratch/pawsey0001/sbeecroft/tool/cromwell-58.jar run /scratch/pawsey0001/sbeecroft/cromwell2/Multisample_Fastq_to_Gvcf_GATK4.wdl \
   -i /scratch/pawsey0001/sbeecroft/cromwell2/Multisample_Fastq_to_Gvcf_GATK4_inputs_hg38.json \
   -o /scratch/pawsey0001/sbeecroft/cromwell2/cromwell.options
