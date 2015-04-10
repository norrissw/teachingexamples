#!/bin/bash


######################################

#$ -cwd
#$ -e ./
#$ -o ./
#$ -S /bin/bash
#$ -N pipeline

######################################
FILE1=$1
#1796_1.fastq
FILE2=$2
#1796_2.fastq
REF=$3
#hg19.fa
SAMBAM=$4
#1796

source /home/norrissw/bin/source_BNFO620.sh
source /usr/global/jdk7/java_env.sh
java -version
env

#1. Use FastQC to generate statistics
#fastqc ${FILE1} 
#fastqc ${FILE2}

#Align fastq files with indexed genome
#/usr/global/blp/bin/bwa aln ${REF} ${FILE1} >${FILE1}aln_.sai
#/usr/global/blp/bin/bwa aln ${REF} ${FILE2} >${FILE2}aln_.sai

#Transform .sai to .SAM format  use bwa sampe command for paired end reads
#/usr/global/blp/bin/bwa sampe ${REF} ${FILE1}aln_.sai ${FILE2}aln_.sai ${FILE1} ${FILE2} >${SAMBAM}.sam

#Use samtools to transform SAM to BAM
#/usr/global/blp/bin/samtools view -bS ${SAMBAM}.sam  -o ${SAMBAM}.bam

# Using mem now instead of aln and sampe
/usr/global/blp/bin/bwa mem ${REF} ${FILE1} ${FILE2} > ${FILE1}_merged.sam
#Sort BAM file
/usr/global/blp/bin/samtools sort ${SAMBAM}.bam ${SAMBAM}.sorted

#Add Read Groups
#java -Xmx2g -jar /usr/global/blp/picard-tools-1.95/AddOrReplaceReadGroups.jar INPUT=${SAMBAM}.bam OUTPUT=${SAMBAM}.sorted.bam SORT_ORDER=coordinate RGLB=8 RGPL=Illumina RGPU=1 RGSM=${SAMBAM}

#Remove duplicate reads
#java -jar /usr/global/blp/picard-tools-1.95/MarkDuplicates.jar REMOVE_DUPLICATES=true M=${SAMBAM}.duplicate_metrics.out I=${SAMBAM}.sorted.bam O=${SAMBAM}.no_dups.bam

#Index BAM File - Picard BuildBamIndex
#java -jar /usr/global/blp/picard-tools-1.95/BuildBamIndex.jar INPUT=${SAMBAM}.no_dups.bam

#Stats of  BAM File - picard BamIndexStats
#java -jar /usr/global/blp/picard-tools-1.95/BamIndexStats.jar INPUT=${SAMBAM}.no_dups.bam

#validate BAM file
#java -jar /usr/global/blp/picard-tools-1.95/ValidateSamFile.jar INPUT=${SAMBAM}.no_dups.bam

#generate the sequence dictionary
#java -jar /usr/global/blp/picard-tools-1.95/CreateSequenceDictionary.jar REFERENCE=${REFERENCE} OUTPUT=/gpfs_fs/home/bnfo620/norrissw/data/hg19.dict

#generate the fasta file index
#/usr/global/blp/bin/samtools faidx ${REFERENCE}

#Realignment - GenomeAnalysisTk.jar
#java -jar /usr/global/blp/GenomeAnalysisTK-3.1.1/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${REFERENCE} -o ${SAMBAM}.bam.list -I ${SAMBAM}.no_dups.bam
#java -jar /usr/global/blp/GenomeAnalysisTK-3.1.1/GenomeAnalysisTK.jar -I ${SAMBAM}.no_dups.bam -R /gpfs_fs/home/bnfo620/norrissw/data/hg19.fa -T IndelRealigner -targetIntervals ${SAMBAM}.bam.list -o ${SAMBAM}.realigned.bam
#java -jar /usr/global/blp/picard-tools-1.95/FixMateInformation.jar INPUT=${SAMBAM}.realigned.bam OUTPUT=${SAMBAM}.realigned_fixed.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true

#Base recalibration
#java -Xmx4g -jar GenomeAnalysisTK.jar -T BaseRecalibrator -I my_reads.bam -R resources/Homo_sapiens_assembly18.fasta \
 #  -knownSites bundle/hg18/dbsnp_132.hg18.vcf \
 #  -knownSites another/optional/setOfSitesToMask.vcf \
 #  -o recal_data.table



#Base Quality Score Recalibration
#java -jar /usr/global/blp/GenomeAnalysisTK-3.1.1/GenomeAnalysisTK.jar -T PrintReads -R ${REFERENCE} -I ${SAMBAM}.realigned_fixed.bam -BQSR recalibration_report.grp -o ${REFERENCE}.final.bam

 
