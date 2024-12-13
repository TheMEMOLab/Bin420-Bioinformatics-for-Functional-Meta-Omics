#!/bin/bash

##############SLURM SCRIPT###################################

## Job name:
#SBATCH --job-name=Chopper
#
## Wall time limit:
#SBATCH --time=24:00:00
###Account
#SBATCH --account=nn9987k
## Other parameters:
#SBATCH --nodes 1
#SBATCH --cpus-per-task 12
#SBATCH --gres=localscratch:50G
#SBATCH --partition=normal
#SBATCH --output=slurm-%x_%j.out
#########################################	



#Variables
RSYNC='rsync -aLhv --no-perms --no-owner --no-group'
input=$1
INDIR=$2
OUTDIR=$3

##Main script

##Activate conda environments ## Arturo

module --quiet purge  # Reset the modules to the system default
module load Miniconda3/23.10.0-1


##Activate conda environments

eval "$(conda shell.bash hook)"
conda activate /cluster/projects/nn9987k/.share/conda_environments/NANOPAKQC/
echo "I am workung with this" $CONDA_PREFIX

###Do some work:########

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

echo "Copy fq file"

time $RSYNC $INDIR/$input.fq.gz .

echo "Nanoplot on RawReads"

echo "Performing NanoPlot"

date +%d\ %b\ %T


time NanoPlot \
        -t $SLURM_CPUS_ON_NODE \
        --fastq $input.fq.gz \
        --N50 \
        --loglength \
        -o $input.Nanoplot.dir


###

echo "Decompress ..."

time gzip -d $input.fq.gz

###
echo "Starting QC cleanning with chooper ..."
date +%d\ %b\ %T

mkdir $OUTDIR/$input\_Chopper/

time cat $input.fq | chopper \
--threads $SLURM_CPUS_ON_NODE \
-q 10 \
-l 1000  > $OUTDIR/$input\_Chopper/$input.chopper.fq

rm $input.fq
echo "NanoPlot on Choppered reads..."

$RSYNC $OUTDIR/$input\_Chopper/$input.chopper.fq .

time parallel ::: "NanoPlot \
        -t 6 \
        --fastq $input.chopper.fq \
        --N50 \
        --loglength \
        -o $input.chopper.Nanoplot.dir" \
        "pigz -p 6 $OUTDIR/$input\_Chopper/$input.chopper.fq"

###
echo "Moving files to $OUTDIR"

###Creating a directory in $OUDIR for results

time $RSYNC $input.chopper.Nanoplot.dir $OUTDIR/$input\_QC/
time $RSYNC $input.Nanoplot.dir $OUTDIR/$input\_QC/


#######
echo "I've done"
date
