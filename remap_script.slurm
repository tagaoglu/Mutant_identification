#!/bin/bash

#SBATCH --mail-user=tugba.agaoglu@unifr.ch
#SBATCH --mail-type=ALL
#SBATCH --job-name="Align_RGI_dedup"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=72:00:00
#SBATCH --mem=32G
#SBATCH --partition=pcourse80
#SBATCH --output=Align_RGI_dedup-%j.out
#SBATCH --error=Align_RGI_dedup-%j.error

#Part 1: trim DNA reads [trimmomatic], align them on the reference genome [bwa] -->.sam and add READ GROUP INFORMATION (@RG) with [bwa]

#define attributes. Do not introduce extensions (e.g. .fastq.gz or.fa) in the sbatch command
REF=$1
READS=$2
SAMPLE=$3
FLOWCELL_ID=$4
FLOWCELL_LANE=$5
THREADS=8

#load required module
export PATH=/software/bin:$PATH;
module use /software/module/
module add vital-it
module add UHTS/Analysis/trimmomatic/0.36
module add UHTS/Analysis/samtools/1.4
module add UHTS/Aligner/bwa/0.7.17
module add UHTS/Analysis/picard-tools/2.9.0

#create a general working directory for _RG.bam files (with @RG) and go there
mkdir -p /data/users/${USER}/1_RG_BAM_dedup
cd /data/users/${USER}/1_RG_BAM_dedup
mkdir RG_BAM_dedup_${READS}
cd RG_BAM_dedup_${READS}

#link REF and SAMPLE/READS locally
ln -s /data/users/tagaoglu/BC7107_20/reference/${REF}.fa ${REF}.fa
ln -s /data/users/tagaoglu/BC7107_20/reads/${SAMPLE}*R1.fastq.gz ${SAMPLE}.R1.fastq.gz
ln -s /data/users/tagaoglu/BC7107_20/reads/${SAMPLE}*R2.fastq.gz ${SAMPLE}.R2.fastq.gz

#with trimmotatic: warning adjust MINLEN read length to sequencing parameters
trimmomatic PE -threads ${THREADS} -phred33 ${SAMPLE}.R1.fastq.gz ${SAMPLE}.R2.fastq.gz ${READS}_R1trim.fastq.gz ${READS}_R1unpaired.fastq.gz ${READS}_R2trim.fastq.gz ${READS}_R2unpaired.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:90

#index REF for bwa
bwa index ${REF}.fa

#map READS onto REF with bwa and add READ GROUP INFORMATION at the same time
#'@RG\tID:'${FLOWCELL_ID}.${FLOWCELL_LANE}'\tSM:'${SAMPLE}'\tPL:illumina\tLB:'lib_${SAMPLE}'\tPU:'${FLOWCELL_ID}.${FLOWCELL_LANE}.${SAMPLE}'
#READ GROUP NAME (ID)={FLOWCELL_ID}.{FLOWCELL_LANE}
#SAMPLE_NAME (SM)={SAMPLE}
#LIBRARY_INFORMATION=lib_${SAMPLE} since I will not have the name of most libraries
#PLATFORM_UNIT (PU)=${FLOWCELL_ID}.${FLOWCELL_LANE}.${SAMPLE}
#PLATFORM (PL)=ILLUMINA
bwa mem -t ${THREADS} -M -R '@RG\tID:'${FLOWCELL_ID}.${FLOWCELL_LANE}'\tSM:'${SAMPLE}'\tPL:illumina\tLB:'lib_${SAMPLE}'\tPU:'${FLOWCELL_ID}.${FLOWCELL_LANE}.${SAMPLE}'' ${REF}.fa ${READS}_R1trim.fastq.gz ${READS}_R2trim.fastq.gz > ${READS}_RG.sam


#Part 2: sort and convert _RG.sam into _RG_sorted.bam (Picard), markduplicates (Picard)


#convert sam to bam with Picard tools
#create a mandatory tpm directory in the current directory (not in /tmp) for memory reasons
mkdir tmp
picard-tools SortSam INPUT=${READS}_RG.sam OUTPUT=${READS}_RG_sorted.bam SORT_ORDER=coordinate TMP_DIR='pwd'/tmp

#mark duplicates with Picard tools
picard-tools MarkDuplicates INPUT=${READS}_RG_sorted.bam OUTPUT=${READS}_RG_sorted_dedup.bam METRICS_FILE=${READS}_metrics.txt TMP_DIR='pwd'/tmp

