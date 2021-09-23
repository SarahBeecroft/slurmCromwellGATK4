# SLURM HPC Cromwell implementation of GATK4 germline variant calling pipeline
See the [GATK](https://gatk.broadinstitute.org/hc/en-us) website for more information on this toolset 
## Assumptions
- Using hg38 human reference genome build
- Running using HPC/SLURM scheduling. This repo was specifically tested on Pawsey Zeus machine, primarily running in the `/scratch` partition. 
- Starting from short-read Illumina paired-end fastq files as input

### Dependencies
The following versions have been tested and work, but GATK and Cromwell are regularly updated and so one must consider whether they would like to use newer versions of these tools. 
- BWA/0.7.15
- GATK v4.0.6.0
- SAMtools/1.5
- picard/2.9
- Python/2.7
- Cromwell v61

## Quick start guide
### Installing and preparing environment for GATK4 with Cromwell

1. Clone repository
```
git clone https://github.com/SarahBeecroft/slurmCromwellGATK4.git
cd slurmCromwellGATK4
chmod +x *.sh
```

2. Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) if you haven’t already. This is best placed in your `/group` directory to avoid filling your small `/home` directory, or being purged is placed in the `/scratch` directory.

3. Create Conda environment using the supplied conda environment file

```
conda env create --file gatk4_pipeline.yml
```

3. Download the necessary .jar files
    - The Cromwell workfow orchestration engine can be downloaded from https://github.com/broadinstitute/cromwell/releases/ 
    - GATK can be downloaded from https://github.com/broadinstitute/gatk/releases. Unzip the file with `unzip` 
    - Picard can be downloaded from https://github.com/broadinstitute/picard/releases/


4. If you do not have the resource bundle files already, these need to be downloaded. In future they will be cached on Pawsey systems. The bundle data should be download from the [Google Cloud bucket](https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0;tab=objects?_ga=2.98248159.1769807612.1582055494-233304531.1578854612&pli=1&prefix=&forceOnObjectsSortingFiltering=false) and not from the FTP site, which is missing various files. Refer to this handy [blog post](https://davetang.org/muse/2020/02/21/using-google-cloud-sdk-to-download-gatk-resource-bundle-files/) on how to download the resource files using Google Cloud SDK. There is a Slurm script (download_bundle.slurm) that can be used to download all hg38 files from the Google Cloud bucket. The files were downloaded in /scratch/pawsey0001/sbeecroft/hg38/v0, which needs to be moved before the data becomes purged after 30 days. Note that Homo_sapiens_assembly38.dbsnp138.vcf.gz was from the FTP bundle as this file could not be downloaded using the Conda version of Google Cloud SDK.

Note that the `hg38_wgs_scattered_calling_intervals.txt` will need to be to generated using the following:

```
cd <your_resource_dir>
find `pwd` -name "scattered.interval_list" -print | sort > hg38_wgs_scattered_calling_intervals.txt
```

These files are required for Multisample_Fastq_to_Gvcf_GATK4.

```
Homo_sapiens_assembly38.dict
Homo_sapiens_assembly38.fasta
Homo_sapiens_assembly38.fasta.fai
Homo_sapiens_assembly38.fasta.64.alt
Homo_sapiens_assembly38.fasta.64.amb
Homo_sapiens_assembly38.fasta.64.ann
Homo_sapiens_assembly38.fasta.64.bwt
Homo_sapiens_assembly38.fasta.64.pac
Homo_sapiens_assembly38.fasta.64.sa
Homo_sapiens_assembly38.fasta.amb
Homo_sapiens_assembly38.fasta.ann
Homo_sapiens_assembly38.fasta.bwt
Homo_sapiens_assembly38.fasta.pac
Homo_sapiens_assembly38.fasta.sa
Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
Homo_sapiens_assembly38.dbsnp138.vcf
Homo_sapiens_assembly38.dbsnp138.vcf.idx
Homo_sapiens_assembly38.known_indels.vcf.gz
Homo_sapiens_assembly38.known_indels.vcf.gz.tbi
```

These files are required for Multisample_jointgt_GATK4.

```
wgs_evaluation_regions.hg38.interval_list
hg38.custom_100Mb.intervals
Homo_sapiens_assembly38.dbsnp138.vcf
Homo_sapiens_assembly38.dbsnp138.vcf.idx
1000G_phase1.snps.high_confidence.hg38.vcf.gz
1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
1000G_omni2.5.hg38.vcf.gz
1000G_omni2.5.hg38.vcf.gz.tbi
Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz
Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi
hapmap_3.3.hg38.vcf.gz
hapmap_3.3.hg38.vcf.gz.tbi
```


5. Set up the config files. Files that you need to edit with the correct paths to your data/jar files or other specific configurations are:
    - `Multisample_Fastq_to_Gvcf_GATK4_inputs_hg38.json`
    - `Multisample_jointgt_GATK4_inputs_hg38.json`
        - both json files will need the correct paths to your reference file locations, and the file specifying your inputs i.e. `samples.txt` or `gvcfs.txt`
    - `samples.txt`
    - `gvcfs.txt`
        - These are the sample input files (tab seperated)
        - The format for samples.txt is sampleID, sampleID_readgroup, path_to_fastq_R1_file, path_to_fastq_R2_file,
        - The format for gvcfs.txt is sample ID, gvcf, gvcf .tbi index file
        - Examples are included in this repo
        - NOTE: Having tabs, not spaces, is vital for parsing the file. Visual studio code tends to introduce spaces, so if you are having issues, check the file with another text editor such as sublime. 
    - `launch_cromwell.sh`
    - `launch_jointgt.sh`
        - These are the scripts which launch the pipeline. 
        - `launch_cromwell.sh` launches the fastq to gvcf stage
        - `launch_jointgt.sh` launched the gvcf joint genotyping to cohort vcf step. This is perfomed when you have run all samples through the fastq to gvcf stage.
        - Check the paths and parameters make sense for your machine
    - `slurm.conf`
        - the main options here relate to the job scheduler. If you are running on Zeus at Pawsey, you should not need to alter these parameters.
    - `cromwell.options`
        - `cromwell.options` requires editing to provide the directory where you would like the final workflow outputs to be written
    - `Multisample_Fastq_to_Gvcf_GATK4.wdl`
    - `ruddle_fastq_to_gvcf_single_sample_gatk4.wdl`
        - The paths to your jar files will need to be updated
        - The path to your conda `activate` binary will need to be updated (e.g. `/group/projectID/userID/miniconda/bin/activate`)

6. Launch the job using `sbatch launch_cromwell.sh`. When that has completed successfully, you can launch the second stage of the pipeline (joint calling) with `sbatch launch_jointgt.sh`.

### Overview of the steps in `Multisample_Fastq_to_Gvcf_GATK4.wdl`
This part of the pipeline takes short-read, Illumina paired-end fastq files as the input. The outputs generated are sorted, duplicate marked bam files and their indices, duplicate metric information, and a GVCF file for each sample. The GVCF files are used as input for the second part of the pipeline (joint genotyping).

```
FastqToUbam
GetBwaVersion
SamToFastqAndBwaMem
MergeBamAlignment
SortAndFixTags
MarkDuplicates
CreateSequenceGroupingTSV
BaseRecalibrator
GatherBqsrReports
ApplyBQSR
GatherBamFiles
HaplotypeCaller
MergeGVCFs
```

### Overview of the steps in `Multisample_jointgt_GATK4.wdl`
This part of the pipeline takes GVCF files (one per sample), and performs joint genotyping across all of the provided samples. This means that old previously generated GVCFs can be joint-called with new GVCFs whenever you need to add new samples. The key output from this is a joint-genotyped, cohort-wide VCF file.

```
GetNumberOfSamples
ImportGVCFs
GenotypeGVCFs
HardFilterAndMakeSitesOnlyVcf
IndelsVariantRecalibrator
SNPsVariantRecalibratorCreateModel
SNPsVariantRecalibrator
GatherTranches
ApplyRecalibration
GatherVcfs
CollectVariantCallingMetrics
GatherMetrics
DynamicallyCombineIntervals
```
