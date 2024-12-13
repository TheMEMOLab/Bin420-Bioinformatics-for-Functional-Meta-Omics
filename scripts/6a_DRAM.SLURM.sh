#!/bin/bash
#########################################################################
#	SLURM scrip for running DRAM annotator on SAGA cluster using all MAGs
#########################################################################

###############SLURM SCRIPT###################################

## Job name:
#SBATCH --job-name=DRAM
#
## Wall time limit:
#SBATCH --time=90:00:00
###Account
#SBATCH --account=nn9864k
## Other parameters:
#SBATCH --nodes 1
#SBATCH --cpus-per-task 16
#SBATCH --mem=80G
#SBATCH --gres=localscratch:50G
#SBATCH --partition=bigmem
#SBATCH --out slurm-%x-%A.out

###########################################################

##########Variables

dir=$1 ##directory with fasta files
ext=$2 #Extension of the fasta files e.g. fasta
OUTDIR=$3 #Outputdirectory
RSYNC='rsync -aLhv --no-perms --no-owner --no-group'

##Activate conda environments ## Arturo

module --quiet purge  # Reset the modules to the system default
module load Miniconda3/23.10.0-1


##Activate conda environments

eval "$(conda shell.bash hook)"

conda activate /cluster/projects/nn9987k/.share/conda_environments/DRAM/ 

####Do some work:########

## For debuggin
echo "Hello" $USER
echo "my submit directory is:"
echo $SLURM_SUBMIT_DIR
echo "this is the job:"
echo $SLURM_JOB_ID
echo "I am running on:"
echo $SLURM_NODELIST
echo "I am running with:"
echo $SLURM_CPUS_ON_NODE "cpus"
echo "Today is:"
date


## Copying data to local node for faster computation

cd $LOCALSCRATCH
echo "copying files to" $LOCALSCRATCH

#Copy the MAGs to the $LOCALSCRATCH for local computation

echo "copying MAGs to" $LOCALSCRATCH

time $RSYNC $dir/ ./MAGs


##################DRAM##############################

echo "DRAM for annotation at"
date +%d\ %b\ %T

time DRAM.py annotate \
-i MAGs'/*.'$ext \
-o dram.annotation.dir \
--min_contig_size 500 \
--threads $SLURM_CPUS_ON_NODE

echo "Distilling ..."
date +%d\ %b\ %T  

time DRAM.py distill \
-i dram.annotation.dir/annotations.tsv \
-o dram.genome_summaries.dir \
--trna_path dram.annotation.dir/trnas.tsv \
--rrna_path dram.annotation.dir/rrnas.tsv
echo "DRAM finished at"
date +%d\ %b\ %T

##Generate a master directory to keep all the results together 

mkdir DRAM.Results.dir
mv dram.annotation.dir DRAM.Results.dir
mv dram.genome_summaries.dir DRAM.Results.dir

###########Moving results to $SLURM_SUBMIT_DIR partition or anywhere the main script was submitted############

echo "moving results to" $OUTDIR/

cd $LOCALSCRATCH

time $RSYNC DRAM.Results.dir $OUTDIR/

echo "DRAM results are in: " $OUTDIR/DRAM.Results.dir


