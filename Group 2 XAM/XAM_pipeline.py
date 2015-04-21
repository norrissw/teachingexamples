#!/usr/bin/env python
# Written by: Shaun Norris
# VCU BNFO 620 - Bioinformatics practicum
# Group 2 - XAM Pipeline

### Imports ###
import subprocess,sys,re
### Global Variables/Constants ###
try:
	F1 = sys.argv[1]
	F2 = sys.argv[2]
	REF = sys.argv[3]
	SB = sys.argv[4]
	SNAME = sys.argv[5]
except:
	print "Not enough arguments provided... Usage: XAM_pipeline.py file1.fastq file2.fastq reference.fa outfilename samplename"
	sys.exit(1)
bwapath="/usr/global/blp/bin/bwa"
sampath="/usr/global/blp/bin/samtools"
java="java -Xmx6g -jar"
picardpath="/usr/global/blp/picard-tools-1.95/"
GATK="/usr/global/blp/GenomeAnalysisTK-3.1.1/GenomeAnalysisTK.jar"
### Functions ###
def QC():
	try:
		print "Running FastQC..."
		subprocess.Popen(['fastqc',F1,F2]).wait()
	except:
		print "FastQC Failed..."
	#pass
def align():
#	try:
#	 /usr/global/blp/bin/bwa mem -t 10 -M ./ucsc.hg19.fasta -R "@RG\tID:TEST_TEST_TEST\tLB:TEST\tSM:TEST\tPL:ILLUMINA" 1796_1.fastq 1796_2.fastq > 1796_norrissw_merged.bam
	RG = '"@RG\\tID:%s_%s_%s\\tLB:%s\\tSM:%s\\tPL:ILLUMINA"' % (SNAME,SNAME,SNAME,SNAME,SNAME)
	merged = SNAME + '_merged.sam'
	print "Running bwa alignment..."
	print bwapath,'mem','-t','10','-M',REF,'-R',RG,F1,F2,'>',merged
	subprocess.Popen((bwapath,'mem','-t','10','-M',REF,'-R',RG,F1,F2,'>',merged)).wait()
#	except:
#		print "Alignment Failed..."
#		sys.exit()
def sam2bam():
	print "Converting SAM to BAM..."
	subprocess.Popen([sampath,'view','-bS','%s_merged.sam','-o','%s.bam' % (SNAME,SNAME)]).wait()
def sort():
	print "Sorting BAM File..."
	subprocess.Popen([sampath,'sort','%s.bam','%s_sorted.bam' % (SNAME,SNAME)]).wait()
def readgroups():
	print "Changing Read Groups..."
	picardjar = picardpath + "AddOrReplaceReadGroups.jar"
	subprocess.Popen([java,picardjar,'INPUT=%s_sorted.bam','OUTPUT=%s_sorted_grouped.bam','SORT_ORDER=coordinate','RGLB=8','RGPL=Illumina','RGPU=1','RGSM=%s' % (SNAME,SNAME)]).wait()
def rmv_dups():
	picardjar = picardpath + "MarkDuplicates.jar"
	print "Removing Duplicate Reads..."
	subprocess.Popen([java,picardjar,'REMOVE_DUPLICATES=true','M=%s_dup_mets.out','I=%s_sorted_grouped.bam','O=%s.nodup.bam' % (SNAME,SNAME,SNAME)]).wait()
def index():
	print "Creating Index..."
	picardjar = picardpath + "BuildBamIndex.jar"
	subprocess.Popen([java,picardjar,'INPUT=%s.nodup.bam' % (SNAME)]).wait()
def index_stats():
	print "Running Index stats..."
	picardjar = picardpath + "BamIndexStats.jar"
	subprocess.Popen([java,picardjar,'INPUT=%s.nodup.bam','>','%s.BamIDXstats.txt' % (SNAME,SNAME)]).wait()
def validate():
	print "Validating BAM File..."
	picardjar = picardpath + "ValidateSamFile.jar"
	subprocess.Popen([java,picardjar,'INPUT=%s.nodup.bam' % (SNAME)]).wait()
def seq_dict():
	print "Generating the Sequence Dictionary..."
	picardjar = picardpath + "CreateSequenceDictionary.jar"
	subprocess.Popen([java,picardjar,'REFERENCE=%s','OUTPUT=%s.dict' % (REF,REF)]).wait()
def fasta_idx():
	print "Creating index..."
	subprocess.Popen([sampath,'faidx',REF]).wait()
def reorder():
	print "Reordering BAM file..."
	picardjar = picardpath + "ReorderSam.jar"
	subprocess.Popen([java,picardjar,'I=%s.nodup.bam','O=%s.nodup_reorder.bam','R=%s' % (SNAME,SNAME,REF)]).wait()
def realign():
	print "Running Realignment..."
	picardjar = picardpath + "FixMateInformation.jar"
	subprocess.Popen([java,GATK,'-T','RealignerTargetCreator','-R',REF,'-o','%s.bam.list','-I','%s.nodup_reorder.bam' %(SNAME,SNAME)]).wait()
	subprocess.Popen([java,GATK,'-T','IndelRealigner','-targetIntervals','%s.bam.list','-I','%s.nodup_reorder.bam','-R',REF,'-o','%s.realigned.bam' % (SNAME,SNAME,SNAME)]).wait()
	subprocess.Popen([java,picardjar,'INPUT=%s.realigned.bam','OUTPUT=%s.realigned_fixmate.bam','SO=coordinate','VALIDATION_STRINGENCY=LENIENT','CREATE_INDEX=true' % (SNAME)]).wait()
def recalibrate():
	print "Running base recalibration..."
	subprocess.Popen([java,GATK,'-T','BaseRecalibrator','-I','%s.realigned_fixmate.bam','-R',REF,'-o','%s.recal.table' % (SNAME)]).wait()
def rescore():
	print "Running base quality score recalibration..."
	subprocess.Popen([java,GATK,'-T','PrintReads','-R',REF,'-I','%s.realigned_fixmate.bam','-BQSR','recalibration_report.grp','-o','%s_rescored.bam' % (SNAME,SNAME)]).wait()
	subprocess.Popen([java,GATK,'-T','IndelGenotyperV2','-R',REF,'-I','%s.realigned_fixmate.bam','-O','%s_indels.txt','--verbose','-o','%s_indel_stats.txt' % (SNAME,SNAME,SNAME)]).wait()
def snps_indels():
	print "Running SNP and InDels analysis"
	subprocess.Popen([java,GATK,'-T','UnifiedGenotyper','-R',REF,'-I','%s.realigned_fixmate.bam','-varout','%s.geli.calls','-vf','GELI','-stand_call_conf','30.0','-stand_emit_conf','10.0','-pl','SOLEXA' % (SNAME,SNAME)]).wait()	

#QC()
align()
sam2bam()
sort()
readgroups()
rmv_dups()
index()
index_stats()
validate()
#seq_dict() #Only needs to run once
#fasta_idx() #Only needs to run once
reorder()
realign()
recalibrate()
rescore()
snps_indels()
