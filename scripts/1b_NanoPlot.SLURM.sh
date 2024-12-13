#!/bin/bash

##############SLURM SCRIPT###################################

## Job name:
#SBATCH --job-name=NanoPlot
#
## Wall time limit:
#SBATCH --time=04:00:00
###Account
#SBATCH --account=nn9987k
## Other parameters:
#SBATCH --nodes 1
#SBATCH --cpus-per-task 12
#SBATCH --mem=26G
#SBATCH --gres=localscratch:100G
#SBATCH --partition=normal,bigmem,hugemem
#SBATCH --output=slurm-%x_%j.out
#########################################	



#Variables
RSYNC='rsync -aPLhv --no-perms --no-owner --no-group'
rawreads=$1
chopped=$2
input=$3
OUTDIR=$4

##Main script

##Activate conda environments ## Arturo

module --quiet purge  # Reset the modules to the system default
module load Anaconda3/2022.10


##Activate conda environments

export PS1=\$
source ${EBROOTANACONDA3}/etc/profile.d/conda.sh
conda deactivate &>/dev/null
conda activate /cluster/projects/nn9987k/.share/conda_environments/NANOPAKQC/
echo "I am workung with this" $CONDA_PREFIX

###Do some work:########

## For debuggin
echo "Hello" $USER
echo "my submit directory is:"
echo $SLURM_SUBMIT_DIR
echo "this is the job:"
echo $SLURM_JOB_ID_\_$SLURM_ARRAY_TASK_ID
echo "I am running on:"
echo $SLURM_NODELIST
echo "I am running with:"
echo $SLURM_CPUS_ON_NODE "cpus"
echo "Today is:"
date

## Copying data to local node for faster computation

cd $LOCALSCRATCH

echo "copying files to" $LOCALSCRATCH

echo "Copy raw fq file"

time $RSYNC $rawreads .
RAW=$(ls -1|grep fq.gz)

echo "Nanoplot on RawReads"

echo "Performing NanoPlot"

date +%d\ %b\ %T


time NanoPlot \
        -t $SLURM_CPUS_ON_NODE \
        --fastq $input.fq.gz \
        --N50 \
        --loglength \
        -o $input.Nanoplot.dir
###Creating a directory in $OUDIR for results
time $RSYNC $input.Nanoplot.dir $OUTDIR/$input\_QC/
rm -rf *

###

time $RSYNC $chopped .

CHOP=$(ls -1|grep fq.gz)

echo "Nanoplot on RawReads"

echo "Performing NanoPlot"

date +%d\ %b\ %T

echo "NanoPlot on Choppered reads..."

time NanoPlot \
        -t $SLURM_CPUS_ON_NODE \
        --fastq $CHOP \
        --N50 \
        --loglength \
        -o $input.chopper.Nanoplot.dir

###
echo "Moving files to $OUTDIR"

###Creating a directory in $OUDIR for results

time $RSYNC $input.chopper.Nanoplot.dir $OUTDIR/$input\_QC/


#######
echo "I've done"
date
