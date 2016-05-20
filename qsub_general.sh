#! /bin/bash
######################################################
#  execute script in current directory
#$ -cwd
#  want any .e/.o stuff to show up here too
#$ -e ./
#$ -o ./
#  shell for qsub to use:
#$ -S /bin/bash
#  name for the job; used by qstat
#$ -N rosetta_model 
#$ -pe smp 10
# specify to avoid getting godel199 (if not needed) ... VNK (120608) 
#$ -q all.q

######################################################

# specify all parameters to your program in double quotes
# Example
# qsub  qsub_general.sh fullpath_of_the_program "param1 param2 param3 param4 ..."
export PATH=$PATH:/usr/global/blp/bin
echo $1
echo $2 $3
$1 $2 $3

echo "Finished $1 ...."
