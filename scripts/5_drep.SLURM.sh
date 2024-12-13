#!/bin/bash

###################################
## Job name:
#SBATCH --job-name=dRep
#
## Wall time limit:
#SBATCH --time=24:00:00
###Account
#SBATCH --account=nn9864k
## Other parameters:
#SBATCH --cpus-per-task 16
#SBATCH --mem=120G
#SBATCH --gres=localscratch:150G
#SBATCH --partition=bigmem
#SBATCH --out slurm-%x_%j.out
#######################################


## Set up job environment:
#set -o errexit  # Exit the script on any error
#set -o nounset  # Treat any unset variables as an error

###Variables###
input=$1 #input string
indir=$2 #Input directory
ext=$3 #extension of fasta file e.g fasta
comp=$4 #completeness score 
con=$5 #Contamination score
outdir=$6 #output directory

##Activate conda environments ## Arturo

module --quiet purge  # Reset the modules to the system default
module load Miniconda3/23.10.0-1

##Activate conda environments

export PS1=\$
eval "$(/cluster/software/Miniconda3/23.10.0-1/bin/conda shell.bash hook)"
conda activate /cluster/projects/nn9987k/.share/conda_environments/DEREPLICATION
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
time rsync -aL $indir/*.$ext .

##Check for complete contigs and unbined files and remove it from the list

unbin=$(ls -1 | grep "unbin")
contigs=$(ls -1 | grep "contigs.fasta")
lowDepth=$(ls -1 | grep "lowDepth.fasta")
tooshort=$(ls -1 | grep "tooShort.fasta")

# Array of file variables
files=("$contigs" "$unbin" "$lowDepth" "$tooshort")

# Loop through files and check existence
for file in "${files[@]}"; do
    if [[ -f "$file" ]]; then
        echo "File $file found, removing it..."
        rm "$file"
    fi
done

###

echo "Number of genomes to drep:"
ls -1|wc -l

mkdir $input.MAGs
mv *.$ext $input.MAGs/

#####dREP dereplicate pipeline################
echo "start dREP at"
date +%d\ %b\ %T

time dRep \
	dereplicate \
	$input.DREP.$comp.$con.out \
	-g $input.MAGs/*.$ext \
	-p $SLURM_CPUS_ON_NODE \
	-comp $comp \
	-con $con

##Copy files to the $SLURM_SUBMIT_DIR
cd $workdir

time rsync -aP $input.DREP.$comp.$con.out $outdir

echo "results are in: " $outdir/$input.DREP.$comp.$con.out

###
echo "I've done at"
date
