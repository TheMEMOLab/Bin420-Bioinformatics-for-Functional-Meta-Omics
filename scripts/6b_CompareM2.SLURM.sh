#!/bin/bash
#########################################################################

###############SLURM SCRIPT###################################

## Job name:
#SBATCH --job-name=COMPAREM2
#
## Wall time limit:
#SBATCH --time=90:00:00
###Account
#SBATCH --account=nn9864k
## Other parameters:
#SBATCH --nodes 1
#SBATCH --cpus-per-task 16
#SBATCH --mem=80G
#SBATCH --gres=localscratch:200G
#SBATCH --partition=bigmem
#SBATCH --out slurm-%x-%A.out

###########################################################

##########Variables

dir=$1 ##directory with fasta files
OUTDIR=$2 #Outputdirectory
RSYNC='rsync -aLhv --no-perms --no-owner --no-group'


##Activate conda environments ## Arturo

module --quiet purge  # Reset the modules to the system default
mmodule load Miniconda3/23.10.0-1


##Activate conda environments

eval "$(conda shell.bash hook)"

conda activate /cluster/projects/nn9987k/.share/conda_environments/COMPAREM2/ 

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


##################COMPAREM2##############################
########Configuration of COMPAREM###############

export COMPAREM2_DATABASES="/cluster/projects/nn9987k/.share/db/COMPAREM2"
export APPTAINER_TMPDIR=$LOCALSCRATCH
export APPTAINER_CACHEDIR=$LOCALSCRATCH

echo "This configuration is running:"

echo "DB" $COMPAREM2_DATABASES
echo "APTTMPDID" $APPTAINER_TMPDIR
echo "APTCACHE" $APPTAINER_CACHEDIR

echo "Starting COMPAREM2"
date +%b\ %d\ %T

time comparem2 \
--cores $SLURM_CPUS_ON_NODE \
--use-singularity \
--singularity-prefix $(pwd) \
--config input_genomes="$(pwd)/MAGs/*.fasta" output_directory="CompareM.out.dir" \
--until assembly_stats checkm2 prokka gtdbtk

###########Moving results to $SLURM_SUBMIT_DIR partition or anywhere the main script was submitted############

echo "moving results to" $OUTDIR/

cd $LOCALSCRATCH

time $RSYNC CompareM.out.dir $OUTDIR/

echo "COMPAREM results are in: " $OUTDIR/CompareM.out.dir


