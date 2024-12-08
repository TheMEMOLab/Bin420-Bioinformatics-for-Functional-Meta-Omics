# Recovering MAGs from Metagenomic data



## 1. Quality control of raw reads.

Let's create a set of directories to work on the FRAM HPC:

```bash
cd /cluster/projects/nn9987k
mkdir $USER
cd $USER
mkdir scripts && mkdir data && mkdir -p results/MetaG
tree
```

This should be the structure of the directory:

```
.
├── data
├── results
│   └── MetaG
└── scripts

4 directories, 0 files

```

Create a softlink of the data:

```bash
ln -s /cluster/projects/nn9987k/.data/MetaG/ONT/rawreads/D01T6_T.fq.gz data/
tree
```

Now that we have all the data we can use a combination of tools to perform a QC and cleanning the reads

Use the SLURM script: ```1_chopper.SLURM.sh``` to clean the reads using [Chopper](https://github.com/wdecoster/chopper) and 
[Nanoplot](https://github.com/wdecoster/NanoPlot) to visualize the QC stats:


This is the template:

```bash
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
RSYNC='rsync -aPLhv --no-perms --no-owner --no-group'
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

time cat $input.fq | chopper \
--threads $SLURM_CPUS_ON_NODE \
-q 10 \
-l 1000  > $input.chopper.fq

echo "Compressing"

pigz -p $SLURM_CPUS_ON_NODE $input.chopper.fq

echo "NanoPlot on Choppered reads..."

time NanoPlot \
        -t $SLURM_CPUS_ON_NODE \
        --fastq $input.chopper.fq.gz \
        --N50 \
        --loglength \
        -o $input.chopper.Nanoplot.dir

###
echo "Moving files to $OUTDIR"

###Creating a directory in $OUDIR for results

time $RSYNC $input.chopper.Nanoplot.dir $OUTDIR/$input\_QC/
time $RSYNC $input.Nanoplot.dir $OUTDIR/$input\_QC/

time $RSYNC $input.chopper.fq.gz $OUTDIR/$input\_Chopper/

#######
echo "I've done"
date

```

A copy of this script is available in the ```/cluster/projects/nn9987k/.scripts/1_chopper.SLURM.sh```

Let's submit the job:

The script requires 3 arguments:
    -Input name
    -Input direcory (Absolut path)
    -Output direcory (Absolut path)

```bash

sbatch /cluster/projects/nn9987k/.scripts/1_chopper.SLURM.sh D01T6_T /cluster/projects/nn9987k/auve/data /cluster/projects/nn9987k/auve/results/MetaG/D01T6_T_Chopper && mkdir -p /cluster/projects/nn9987k/auve/results/MetaG/D01T6_T_Chopper

```

```
Submitted batch job 6011498
```

Monitor the job:

```bash
squeue -u $USER
```

```
JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
           6011498    normal  Chopper     auve  R       1:18      1 c7-10

```

The slurm-Chopper_<JOBID>.out should be generated:

```bash
tree
```

```
.
├── data
│   └── D01T6_T.fq.gz -> /cluster/projects/nn9987k/.data/MetaG/ONT/rawreads/D01T6_T.fq.gz
├── results
│   └── MetaG
│       └── D01T6_T_Chopper
├── scripts
└── slurm-Chopper_6011498.out
```


5 directories, 2 files

The final direcory looks like this:


Let's inspect the Stats of NanoPlot:



## 2. Assembly reads with MetaFlye

To extend the ONT reads we will use [Flye]() assembler with the ```--meta``` flag:

This is the template and the arguments are:

    - Input Name
    - Iniput dir
    -Output dir



```bash

#!/bin/bash

##############SLURM SCRIPT###################################

## Job name:
#SBATCH --job-name=MetaFly
#
## Wall time limit:
#SBATCH --time=72:00:00
###Account
#SBATCH --account=nn9987k
## Other parameters:
#SBATCH --nodes 1
#SBATCH --cpus-per-task 16
#SBATCH --mem=90G
#SBATCH --gres=localscratch:250G
#SBATCH --partition=bigmem
#SBATCH --output=slurm-%x_%j.out
#########################################	

###Basic usage help for this script#######

print_usage() {
        echo "Usage: sbatch $0 input indir outputdir"
}

if [ $# -lt 3 ]
        then
                print_usage
                exit 1
        fi


###############Main SCRIPT####################

##Variables###

input=$1
INDIR=$2
outdir=$3
RSYNC='rsync -aLhv --no-perms --no-owner --no-group'



##Activate conda environments ## Arturo

module --quiet purge  # Reset the modules to the system default
module load Miniconda3/23.10.0-1

##Activate conda environments

export PS1=\$
eval "$(/cluster/software/Miniconda3/23.10.0-1/bin/conda shell.bash hook)"
conda deactivate &>/dev/null

conda activate /cluster/projects/nn9987k/.share/conda_environments/MetaG_Assembly_And_Binning

echo "I'm working with this CONDAENV"
echo $CONDA_PREFIX

###Do some work:########

## For debuggin
echo "Hello" $USER
echo "my submit directory is:"
echo $SLURM_SUBMIT_DIR
echo "this is the job:"
echo $SLURM_JOB_ID\_$SLURM_ARRAY_TASK_ID
echo "I am running on:"
echo $SLURM_NODELIST
echo "I am running with:"
echo $SLURM_CPUS_ON_NODE "cpus"
echo "Today is:"
date

## Copying data to local node for faster computation

cd $LOCALSCRATCH

echo "copying Reads to" $LOCALSCRATCH

$RSYNC $INDIR/$input.*gz .

####Assembly#######################

echo "Starting assembly by Flye...."
date +%d\ %b\ %T

time flye \
--nano-raw $input.*gz \
--meta \
--out-dir $input.flye.outdir \
-t $SLURM_CPUS_ON_NODE

echo "Final results are in: "$outdir

$RSYNC $input.flye.outdir $outdir/

####removing tmp dir. Remember to do this for not filling the HDD in the node!!!!###

echo "I've done at"
date

```

Running using sbatch:


```bash
sbatch /cluster/projects/nn9987k/.scripts/2_flye.SLURM.chr.sh D01T6_T /cluster/projects/nn9987k/.auve/results/MetaG/D01T6_T_Chopper/ /cluster/projects/nn9987k/auve/results
```

Unfortunatelly for the BIN240 Course Sigma2 has only assigned a copuple of nodes in FRAM computer, so it is most likely the Job never runs. But you can copy the results of this assembly by:

```bash
rsync -aLhv /cluster/projects/nn9987k/.results/MetaG/D01T6_T.flye.outdir/assembly* /cluster/projects/nn9987k/$USER/results/MetaG/D01T6_T.flye.outdir/
```
[!WARNING]
!NB: Remember to kill the Fly job by ```scancel <JOBID>```

