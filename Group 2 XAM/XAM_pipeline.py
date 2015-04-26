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
java="java"
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
	RG = '@RG\\tID:%s_%s_%s\\tLB:%s\\tSM:%s\\tPL:ILLUMINA' % (SNAME,SNAME,SNAME,SNAME,SNAME)
	merged = SNAME + '_aligned.sam'
	mergedfile = open(merged,'w')
	print "Running bwa alignment..."
	#print bwapath,'mem','-t','10',REF,'-R',RG,F1,F2,'>',merged
	memalign = subprocess.Popen(['bwa','mem','-t','20','-M',REF,'-R',RG,F1,F2], stdout=subprocess.PIPE)#,stdin=subprocess.PIPE,stderr=subprocess.PIPE).wait() #maybe use .call?
	memout = memalign.communicate()[0]
#	print memout
	mergedfile.write(memout)
	return merged
#	except:
#		print "Alignment Failed..."
#		sys.exit()
def sam2bam(merged):
	print "Converting SAM to BAM..."
	BAM = merged.strip('.sam') + '.bam'
	subprocess.Popen([sampath,'view','-bS',merged,'-o',BAM]).wait()
	return BAM
def sort(BAM):
	print "Sorting BAM File..."
	SORTED = BAM.strip('.bam') + '_sorted'
	subprocess.Popen([sampath,'sort',BAM,SORTED]).wait()
	return SORTED
def readgroups(SORTED):
	# Not used anymore #
	print "Changing Read Groups..."
	picardjar = picardpath + "AddOrReplaceReadGroups.jar"
	subprocess.Popen([java,picardjar,'INPUT=%s_sorted.bam','OUTPUT=%s_sorted_grouped.bam','SORT_ORDER=coordinate','RGLB=8','RGPL=Illumina','RGPU=1','RGSM=%s' % (SNAME,SNAME)]).wait()
def rmv_dups(SORTED):
	picardjar = picardpath + "MarkDuplicates.jar"
	print "Removing Duplicate Reads..."
	MOUT = 'M=' + SNAME + '_dup_mets.out'
	IN = 'I=' + SORTED + '.bam'
	OUT = 'O=' + SNAME +'.nodup.bam'
	subprocess.Popen([java,'-Xmx6g','-jar',picardjar,'REMOVE_DUPLICATES=true',MOUT,IN,OUT]).wait()
	return (SNAME + '.nodup.bam')
def index(NODUPS):
	INPUT = 'I=' + SNAME + '.nodup.bam'
	print "Creating Index..."
	picardjar = picardpath + "BuildBamIndex.jar"
	subprocess.Popen([java,'-Xmx6g','-jar',picardjar,INPUT]).wait()
def index_stats(IDX):
	print "Running Index stats..."
	picardjar = picardpath + "BamIndexStats.jar"
	INPUT = 'INPUT=' + IDX
	OUTPUT = SNAME + '.bam_IDXstats.txt'
	outfile = open(OUTPUT,'w')
	idx_p = subprocess.Popen([java,'-Xmx6g','-jar',picardjar,INPUT], stdout=subprocess.PIPE)
	idxout = idx_p.communicate()[0]
	outfile.write(idxout)
def validate(NODUPS):
	print "Validating BAM File..."
	picardjar = picardpath + "ValidateSamFile.jar"
	INPUT = 'INPUT=' + NODUPS
	subprocess.Popen([java,'-Xmx6g','-jar',picardjar,INPUT]).wait()
def seq_dict():
	print "Generating the Sequence Dictionary..."
	picardjar = picardpath + "CreateSequenceDictionary.jar"
	REFERENCE = 'REFERENCE=%s' % REF
	OUTPUT = 'OUTPUT=%s.dict' % REF
	subprocess.Popen([java,'-Xmx6g','-jar',picardjar,REFERENCE,OUTPUT]).wait()
def fasta_idx():
	print "Creating index..."
	subprocess.Popen([sampath,'faidx',REF]).wait()
def reorder(NODUPS):
	print "Reordering BAM file..."
	picardjar = picardpath + "ReorderSam.jar"
	INPUT = 'I=' + NODUPS
	OUTPUT = 'O=%s.nodup_reorder.bam' % SNAME
	REFERENCE = 'R=' + REF
	subprocess.Popen([java,'-Xmx6g','-jar',picardjar,INPUT,OUTPUT,REFERENCE]).wait()
	return ('%s.nodup_reorder.bam' % SNAME)
def realign(REORDER):
	print "Running Realignment..."
	picardjar = picardpath + "FixMateInformation.jar"
	OUTLIST = '%s.bam.list' % SNAME
	REALIGNED = '%s.realigned.bam' % SNAME
	print java,'-Xmx6g','-jar',GATK,'-T','RealignerTargetCreator','-R',REF,'-o',OUTLIST,'-I',REORDER
	subprocess.Popen([java,'-Xmx6g','-jar',GATK,'-T','RealignerTargetCreator','-R',REF,'-o',OUTLIST,'-I',REORDER]).wait()
	subprocess.Popen([java,'-Xmx6g','-jar',GATK,'-T','IndelRealigner','-targetIntervals',OUTLIST,'-I',REORDER,'-R',REF,'-o',REALIGNED]).wait()
	INPUT = "INPUT=%s" % REALIGNED
	OUTPUT = "OUTPUT=%s.realigned_fixmate.bam" % SNAME
	subprocess.Popen([java,picardjar,INPUT,OUTPUT,'SO=coordinate','VALIDATION_STRINGENCY=LENIENT','CREATE_INDEX=true']).wait()
def recalibrate():
	print "Running base recalibration..."
	subprocess.Popen([java,'-Xmx6g','-jar',GATK,'-T','BaseRecalibrator','-I','%s.realigned_fixmate.bam','-R',REF,'-o','%s.recal.table' % (SNAME)]).wait()
def rescore():
	print "Running base quality score recalibration..."
	subprocess.Popen([java,'-Xmx6g','-jar',GATK,'-T','PrintReads','-R',REF,'-I','%s.realigned_fixmate.bam','-BQSR','recalibration_report.grp','-o','%s_rescored.bam' % (SNAME,SNAME)]).wait()
	subprocess.Popen([java,'-Xmx6g','-jar',GATK,'-T','IndelGenotyperV2','-R',REF,'-I','%s.realigned_fixmate.bam','-O','%s_indels.txt','--verbose','-o','%s_indel_stats.txt' % (SNAME,SNAME,SNAME)]).wait()
def snps_indels():
	print "Running SNP and InDels analysis"
	subprocess.Popen([java,'-Xmx6g','-jar',GATK,'-T','UnifiedGenotyper','-R',REF,'-I','%s.realigned_fixmate.bam','-varout','%s.geli.calls','-vf','GELI','-stand_call_conf','30.0','-stand_emit_conf','10.0','-pl','SOLEXA' % (SNAME,SNAME)]).wait()	

#QC()
#merged = align() COMMENTED OUTFOR TESTING
#merged = 'norrissw_796_merged.sam' #for testing only
#BAM = sam2bam(merged)
#SORTED = sort(BAM)
#readgroups(SORTED) #SHOULD NO LONGER BE NECESSARY
#SORTED = 'norrissw_796_merged_sorted' #for testing only
#NODUPS = rmv_dups(SORTED)
#print NODUPS, 'norrissw_796.nodup.bam' #for testing only
#index(NODUPS)
#index_stats(NODUPS)
#validate(NODUPS)
#seq_dict() #Only needs to run once
#fasta_idx() #Only needs to run once
NODUPS = 'norrissw_796.nodup.bam'
REORDER = reorder(NODUPS)
REALIGNED = realign(REORDER)
recalibrate()
rescore()
snps_indels()
