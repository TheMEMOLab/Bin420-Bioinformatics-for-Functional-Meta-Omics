#!/bin/bash

##############SLURM SCRIPT###################################

## Job name:
#SBATCH --job-name=Binning
#
## Wall time limit:
#SBATCH --time=24:00:00
###Account
#SBATCH --account=nn9987k
## Other parameters:
#SBATCH --nodes 1
#SBATCH --cpus-per-task 12
#SBATCH --gres=localscratch:100G
#SBATCH --output=slurm-%x_%j.out
#########################################	

#Variables
RSYNC='rsync -aLhv --no-perms --no-owner --no-group'
input=$1
READIR=$2
ASSDIR=$3
OUTDIR=$4

##Activate conda environments ## Arturo

##Activate conda environments ## Arturo

module --quiet purge  # Reset the modules to the system default
module load Miniconda3/23.10.0-1


##Activate conda environments

eval "$(conda shell.bash hook)"
conda activate /cluster/projects/nn9987k/.share/conda_environments/MetaG_Assembly_And_Binning/
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

##createing a directory for Metabat
mkdir $input.Binning.dir
cd $input.Binning.dir

echo "Copy fq file"

time $RSYNC $READIR/$input.chopper.fq.gz .

echo "Copy assembly"

time $RSYNC $ASSDIR/$input.medaka.consensus.fasta ./$input.assembly.fasta

##Aling

echo "Start minimap2 "
date +%d\ %b\ %T

time minimap2 \
-ax map-ont \
-t $SLURM_CPUS_ON_NODE \
$input.assembly.fasta \
$input.chopper.fq.gz > $input.assembly.sam

#Get the bam and the files

echo "Samtools view.."
 
time samtools view \
-@ $SLURM_CPUS_ON_NODE \
-bS $input.assembly.sam > $input.assembly.bam

rm -r $input.assembly.sam
rm -r $input.chopper.fq.gz

echo "samtools sort"

time samtools sort \
-@ $SLURM_CPUS_ON_NODE $input.assembly.bam > $input.assembly.sorted.bam

rm -r $input.assembly.bam

###Binning

echo "Calcaulating deepth..."

time jgi_summarize_bam_contig_depths \
        --outputDepth $input.depth.txt \
        $input.assembly.sorted.bam

echo "Modifying $input.depth.txt to a MaxBin2 format:"

cut -f1,3 $input.depth.txt | tail -n+2 > $input.depth_maxbin.txt  

##Using Parallel to run metaba2 and Maxbin2 at the same time with 8 CPUs each...

cpu=$(($SLURM_CPUS_ON_NODE/2))
echo "Running parallel with $cpu cpus"

time parallel -j2 ::: "metabat2 \
        -i $input.assembly.fasta \
        -a $input.depth.txt \
        -m 1500 \
        --seed 100 \
        -t $cpu \
        --unbinned \
        -o $input.Metabat2" \
        "run_MaxBin.pl \
        -contig $input.assembly.fasta \
        -out $input.MaxBin.out \
        -abund $input.depth_maxbin.txt \
        -thread $cpu"

#Changing suffix fa to fasta useful for further analysis

for i in *.fa;
        do
        a=$(basename $i .fa);
        mv $i $a.fasta;
done

#Remove the original assemlby

rm $input.assembly.fasta

###
echo "Moving Binnig files to $OUTDIR"

cd $LOCALSCRATCH

time $RSYNC $input.Binning.dir $OUTDIR/

#######
echo "I've done"
date
