#!/bin/bash

###################################
## Job name:
#SBATCH --job-name=CheckM2
#
## Wall time limit:
#SBATCH --time=24:00:00
###Account
#SBATCH --account=nn9864k
## Other parameters:
#SBATCH --cpus-per-task 16
#SBATCH --mem=80G
#SBATCH --gres=localscratch:150G
#SBATCH --partition=normal,bigmem,hugemem
#SBATCH --out slurm-%x_%j.out
#######################################

###Variables###
input=$1 #input string
indir=$2 #Input directory
ext=$3 #extension of fasta file e.g fasta
outdir=$4 #
RSYNC='rsync -aLhv --no-perms --no-owner --no-group'
##Activate conda environments ## Arturo

module --quiet purge  # Reset the modules to the system default
module load Anaconda3/2022.10


##Activate conda environments

export PS1=\$
source ${EBROOTANACONDA3}/etc/profile.d/conda.sh
conda deactivate &>/dev/null
conda activate /cluster/projects/nn9987k/.share/conda_environments/CHEKM2
echo "I am workung with this" $CONDA_PREFIX

####Do some work:########

echo "Hello" $USER
echo "my submit directory is:"
echo $SLURM_SUBMIT_DIR
echo "this is the job:"
echo $SLURM_JOB_ID
echo "I am running on:"
echo $SLURM_NODELIST
echo "I am running with:"
echo $SLURM_CPUS_ON_NODE "cpus"
echo "I am working with this enviroment loaded"
echo $CONDA_PREFIX
echo "Today is:"
date

## Copying data to local node for faster computation

cd $LOCALSCRATCH

workdir=$(pwd)
echo "My working directory is" $workdir
echo "copying MAGs files ..."
time $RSYNC $indir/*.$ext .

mkdir $input.MAGs
mv *.$ext $input.MAGs/

#####CHECKM2

echo "Start checkm2"

date +%b\ %d\ %T

checkm2 predict \
--threads $SLURM_CPUS_ON_NODE \
--input $input.MAGs \
-x $ext \
--output-directory $input.MAGs.checkm2.dir \
--database_path /cluster/projects/nn9987k/.share/db/CheckM2_database/uniref100.KO.1.dmnd

ls

time $RSYNC $input.MAGs.checkm2.dir $outdir/

echo "results are in: " $outdir/$input.MAGs.checkm2.dir

###
echo "I've done at"
date