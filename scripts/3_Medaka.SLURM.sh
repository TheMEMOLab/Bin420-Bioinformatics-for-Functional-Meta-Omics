#!/bin/bash

##############SLURM SCRIPT###################################

## Job name:
#SBATCH --job-name=MedakaPolishing
#
## Wall time limit:
#SBATCH --time=2:00:00
###Account
#SBATCH --account=nn9987k
## Other parameters:
#SBATCH --nodes 1
#SBATCH --cpus-per-task 16
#SBATCH --gres=localscratch:50G
#SBATCH --output=slurm-%x_%j.out
#########################################

#Variables
RSYNC='rsync -aLhv --no-perms --no-owner --no-group'
input=$1
READIR=$2
ASSDIR=$3
OUTDIR=$4

##Main script

##Activate conda environments ## Arturo

module --quiet purge  # Reset the modules to the system default
module load Miniconda3/23.10.0-1

##Activate conda environments

eval "$(conda shell.bash hook)"
conda activate /cluster/projects/nn10039k/shared/condaenvironments/MEDAKA
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

time $RSYNC $READIR/$input.chopper.fq.gz .

echo "Copy assembly"

time $RSYNC $ASSDIR/$input.flye.outdir/assembly.fasta ./$input.assembly.fasta

##MEdaking

echo "Starting Medaka..."
date +%d\ %b\ %T

time medaka_consensus \
-i $input.chopper.fq.gz \
-d $input.assembly.fasta \
-o $input.medaka.dir \
-t $SLURM_CPUS_ON_NODE

echo "Cleaning and changing names..."

cd $input.medaka.dir
echo "I am on:"
pwd
mv consensus.fasta $input.toto

###Cleaning

ls -1|grep -v toto|\
while read -r line; 
    do
    rm -r $line;
done

mv $input.toto $input.medaka.consensus.fasta

##moving resutls

echo "Rsync results to $OUTDIR"
cd $LOCALSCRATCH
$RSYNC $input.medaka.dir $OUTDIR/

###
echo "I've done"
date

