# Recovering MAGs from Metagenomic data

<img src="https://github.com/TheMEMOLab/Bin420-Bioinformatics-for-Functional-Meta-Omics/blob/main/img/Assemlby.webp" height="400">

## 1. Quality control of raw reads.

Let's create a set of directories to work on the FRAM HPC:

```bash
cd /cluster/projects/nn9987k/$USER
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

### Running Chopper and NanoPlot.

Now that we have all the data we can use a combination of tools to perform a QC and cleanning the reads

Use the SLURM script: ```1_chopper.SLURM.sh``` to clean the reads using [Chopper](https://github.com/wdecoster/chopper) and 
[Nanoplot](https://github.com/wdecoster/NanoPlot) to visualize the QC stats:

<details>
<summary>This is the template:</summary>

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
</details>

A copy of this script is available in the ```/cluster/projects/nn9987k/.scripts/1_chopper.SLURM.sh```

Let's submit the job:

The script requires 3 arguments:
    -Input name
    -Input direcory (Absolut path)
    -Output direcory (Absolut path)

```bash

sbatch /cluster/projects/nn9987k/.scripts/1_chopper.SLURM.sh D01T6_T /cluster/projects/nn9987k/$USER/data /cluster/projects/nn9987k/$USER/results/MetaG/D01T6_T_Chopper && mkdir -p /cluster/projects/nn9987k/$USER/results/MetaG/D01T6_T_Chopper

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

Let's take a look on these files:

> [!Important]
> As we will start working with files let's ask for an interactive session in FRAM:

```bash
srun \
--account=nn9987k \
--gres=localscratch:20G \
--cpus-per-task 4 \
--nodes 1 \
--time=02:00:00 \
--pty bash \
-i
```

```
srun: job 6019143 queued and waiting for resources
srun: job 6019143 has been allocated resources
[auve@c5-1.FARM: ~]$
```
Now we are loged into a computing node. Let's work:

```bash
 cd /cluster/projects/nn9987k/$USER/results/MetaG/D01T6_T_Chopper/
cd D01T6_T_QC
```

The pipeline created 2 folders here. Use tree to display the content:


<details>


<summary>tree</summary>

```
/cluster/projects/nn9987k/auve/results/MetaG/D01T6_T_Chopper/D01T6_T_QC
├── D01T6_T.chopper.Nanoplot.dir
│   ├── LengthvsQualityScatterPlot_dot.html
│   ├── LengthvsQualityScatterPlot_dot.png
│   ├── LengthvsQualityScatterPlot_kde.html
│   ├── LengthvsQualityScatterPlot_kde.png
│   ├── LengthvsQualityScatterPlot_loglength_dot.html
│   ├── LengthvsQualityScatterPlot_loglength_dot.png
│   ├── LengthvsQualityScatterPlot_loglength_kde.html
│   ├── LengthvsQualityScatterPlot_loglength_kde.png
│   ├── NanoPlot_20241206_1637.log
│   ├── NanoPlot-report.html
│   ├── NanoStats.txt
│   ├── Non_weightedHistogramReadlength.html
│   ├── Non_weightedHistogramReadlength.png
│   ├── Non_weightedLogTransformed_HistogramReadlength.html
│   ├── Non_weightedLogTransformed_HistogramReadlength.png
│   ├── WeightedHistogramReadlength.html
│   ├── WeightedHistogramReadlength.png
│   ├── WeightedLogTransformed_HistogramReadlength.html
│   ├── WeightedLogTransformed_HistogramReadlength.png
│   ├── Yield_By_Length.html
│   └── Yield_By_Length.png
└── D01T6_T.Nanoplot.dir
    ├── LengthvsQualityScatterPlot_dot.html
    ├── LengthvsQualityScatterPlot_dot.png
    ├── LengthvsQualityScatterPlot_kde.html
    ├── LengthvsQualityScatterPlot_kde.png
    ├── LengthvsQualityScatterPlot_loglength_dot.html
    ├── LengthvsQualityScatterPlot_loglength_dot.png
    ├── LengthvsQualityScatterPlot_loglength_kde.html
    ├── LengthvsQualityScatterPlot_loglength_kde.png
    ├── NanoPlot_20241206_1613.log
    ├── NanoPlot-report.html
    ├── NanoStats.txt
    ├── Non_weightedHistogramReadlength.html
    ├── Non_weightedHistogramReadlength.png
    ├── Non_weightedLogTransformed_HistogramReadlength.html
    ├── Non_weightedLogTransformed_HistogramReadlength.png
    ├── WeightedHistogramReadlength.html
    ├── WeightedHistogramReadlength.png
    ├── WeightedLogTransformed_HistogramReadlength.html
    ├── WeightedLogTransformed_HistogramReadlength.png
    ├── Yield_By_Length.html
    └── Yield_By_Length.png

2 directories, 42 files
```

</details>


There are a bunch of files but we really need the stats to compare, for example

```bash
head -8 D01T6_T.Nanoplot.dir/NanoStats.txt
```

```
General summary:
Mean read length:                  5,242.1
Mean read quality:                    17.6
Median read length:                4,404.0
Median read quality:                  18.9
Number of reads:               2,332,038.0
Read length N50:                   5,894.0
STDEV read length:                 2,930.4

```

We can use R to plot these and easily. We can do this:

```bash
cp D01T6_T.Nanoplot.dir/NanoStats.txt $LOCALSCRATCH/D01T6_T_Raw.txt && cp D01T6_T.chopper.Nanoplot.dir/NanoStats.txt $LOCALSCRATCH/D01T6_T_choper.txt
tree $LOCALSCRATCH
```
```
/localscratch/6019143
├── D01T6_T_choper.txt
└── D01T6_T_Raw.txt

0 directories, 2 files

```

R is intalled as conda environment so let's call Miniconda3 and R_env

```bash
module load Miniconda3/23.10.0-1
eval "$(conda shell.bash hook)"
conda activate /cluster/projects/nn9987k/.share/conda_environments/R_env/
```

<details>

And then we can run the following:

<summary> R script</summary>

```R
# Load required libraries
library(tidyverse)

# Function to read and process data from a file
read_summary <- function(file, label) {
  read_lines(file) %>%
    # Split each line into name and value
    str_split_fixed(":\\s+", 2) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    # Name the columns
    rename(Metric = V1, Value = V2) %>%
    # Convert Value column to numeric
    mutate(Value = as.numeric(gsub(",", "", Value)),
           File = label) # Add file label
}

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Ensure two arguments are provided
if (length(args) != 2) {
  stop("Usage: Rscript plot_metrics.R <file_A> <file_B>")
}

# File paths from arguments
file_a <- args[1]
file_b <- args[2]

# Extract file labels (names without .txt)
label_a <- tools::file_path_sans_ext(basename(file_a))
label_b <- tools::file_path_sans_ext(basename(file_b))

# Read data from files with labels
data_a <- read_summary(file_a, label_a)
data_b <- read_summary(file_b, label_b)

# Combine data from both files
combined_data <- bind_rows(data_a, data_b)

# Filter and scale the required metrics
filtered_data <- combined_data %>%
  filter(Metric %in% c("Mean read length", "Mean read quality",
                       "Median read length", "Number of reads", "Total bases")) %>%
  # Scale metrics
  mutate(Value = case_when(
           Metric == "Mean read length" ~ Value / 1e3,
           Metric == "Median read length" ~ Value / 1e3,
           Metric == "Total bases" ~ Value / 1e9,
           Metric == "Number of reads" ~ Value / 1e6,
           TRUE ~ Value
         ),
         Metric = case_when(
           Metric == "Mean read length" ~ "Mean read length (Thousands)",
           Metric == "Median read length" ~ "Median read length (Thousands)",
           Metric == "Total bases" ~ "Total bases (Billions)",
           Metric == "Number of reads" ~ "Number of reads (Millions)",
           TRUE ~ Metric
         ))

# Check unique File values to debug issues
print(unique(filtered_data$File))

# Plot the data
plot <- ggplot(filtered_data, aes(x = Metric, y = Value, fill = File)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Comparison of Selected Metrics",
       x = "Metric",
       y = "Value",
       fill = "File") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # Use File names directly for colors
  scale_fill_manual(values = setNames(c("#1f78b4", "#33a02c"), c(label_a, label_b)))

# Save the plot to a file
output_file <- "comparison_NanoStats.pdf"
ggsave(output_file, plot, width = 10, height = 6)
cat("Plot saved to", output_file, "\n")


```

</details>

```bash
Rscript /cluster/projects/nn9987k/.scripts/plotNanoStats.r $LOCALSCRATCH/D01T6_T_Raw.txt $LOCALSCRATCH/D01T6_T_choper.txt
```

This will produce a pdf file:

```
Plot saved to comparison_NanoStats.pdf
```
![NanoStats](https://github.com/TheMEMOLab/Bin420-Bioinformatics-for-Functional-Meta-Omics/blob/main/img/NanoStats.PNG)

> [!Note]
> To copy data from FRAM to our laptop we can use SCP:

```bash
scp auve@fram.sigma2.no:/cluster/projects/nn9987k/.auve/results/MetaG/D01T6_T_Chopper/D01T6_T_QC/comparison_NanoStats.pdf .
```

> Rember to change ```auve``` to your user name in FRAM

**What can we say about this plot?**

> [!Important]
> Remember to finish your interactive session by ```exit```

## 2. Assembly reads with MetaFlye

To extend the ONT reads we will use [Flye]() assembler with the ```--meta``` flag:

<details>
<summary>This is the template</summary>

The arguments are:

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

</details>

### Running MetaFlye

Running using sbatch:


```bash
sbatch /cluster/projects/nn9987k/.scripts/2_flye.SLURM.chr.sh D01T6_T /cluster/projects/nn9987k/$USER/results/MetaG/D01T6_T_Chopper/ /cluster/projects/nn9987k/$USER/results
```
> [!NOTE]
> Unfortunatelly for the BIN240 Course Sigma2 has only assigned a copuple of nodes in the FRAM computer, so it is most likely the Job never runs/finish.

But you can copy the results of this assembly by:

```bash
rsync -aLhv /cluster/projects/nn9987k/.results/MetaG/D01T6_T.flye.outdir/assembly* /cluster/projects/nn9987k/$USER/results/MetaG/D01T6_T.flye.outdir/
```

> [!WARNING] 
> !NB: Remember to kill the Fly job by scancel ```<JOBID>``` .

Let's take a look on these files:

> [!Important]
> As we will start working with files let's ask for an interactive session in FRAM:

```bash
srun \
--account=nn9987k \
--gres=localscratch:20G \
--cpus-per-task 4 \
--nodes 1 \
--time=02:00:00 \
--pty bash \
-i
```

```
srun: job 6019143 queued and waiting for resources
srun: job 6019143 has been allocated resources
[auve@c5-1.FARM: ~]$
```
Now we are loged into a computing node.

Display the content of the Flye assembly folder:

```bash
ls /cluster/projects/nn9987k/auve/results/MetaG/D01T6_T.flye.outdir/
```

```
assembly.fasta  assembly_graph.gfa  assembly_graph.gv  assembly_info.txt
```

We can use ```assembly-stats``` tool to check for the main stats in the assembly:

```bash
module load Miniconda3/23.10.0-1
eval "$(conda shell.bash hook)"
conda activate /cluster/projects/nn9987k/.share/conda_environments/MetaG_Assembly_And_Binning/
assembly-stats /cluster/projects/nn9987k/auve/results/MetaG/D01T6_T.flye.outdir/assembly.fasta
```

```
stats for /cluster/projects/nn9987k/auve/results/MetaG/D01T6_T.flye.outdir/assembly.fasta
sum = 821745693, n = 45282, ave = 18147.29, largest = 1393968
N50 = 25699, n = 8036
N60 = 19600, n = 11713
N70 = 15449, n = 16442
N80 = 11905, n = 22498
N90 = 8290, n = 30706
N100 = 63, n = 45282
N_count = 0
Gaps = 0
```
> [!Important]
> Remember to finish your interactive session by ```exit```

## 3. Polishing.

A basic model of how polishing works is that the polisher stacks all relevant reads on top of the genome and decides for each position whether the present nucleotide letter is the best representative for that position, or not. There are several sources of variation that make draft assemblies polishable. The main sources are multi-strain variation from closely related species as well as incorporation of sequencing errors during the sequencing process. Ideally, assemblers would be perfect, and we wouldn't have to perform polishing. But because of some noise or artefacts that are present in our data, we might make our genomes more truthful to their biological origin by performing these polishing steps.

Genome polishing is reminiscent of generating a consensus genome. Consensus genome creation is a term used in reference mapping. This is why you may incidentally see the term consensus being used in the tools that we're gonna run.

[medaka](https://github.com/nanoporetech/medaka) is a tool to create consensus sequences and variant calls from nanopore sequencing data. This task is performed using neural networks applied a pileup of individual sequencing reads against a reference sequence, mostly commonly either a draft assembly or a database reference sequence. It provides state-of-the-art results outperforming sequence-graph based methods and signal-based methods, whilst also being faster.

**Features**

    -Requires only basecalled data. (.fasta or .fastq)
    -Improved accuracy over graph-based methods (e.g. Racon).
    -50X faster than Nanopolish (and can run on GPUs).
    -Includes extras for implementing and training bespoke correction networks.
    -Works on Linux and MacOS.
    -Open source (Oxford Nanopore Technologies PLC. Public License Version 1.0)

### Running Medaka.

Medaka needs two parameters to run:

- Fasta file of the assembly.
- Fastq files used for the assembly.

The main syntax would be something like this:

```bash
medaka_consensus \
-i $input.chopper.fq.gz \
-d $input.assembly.fasta \
-o $input.medaka.dir \
-t $SLURM_CPUS_ON_NODE

```
<details>

<summary> The following SLURM script can be used to perform the Medaka polishing </summary>

```bash
#!/bin/bash

##############SLURM SCRIPT###################################

## Job name:
#SBATCH --job-name=MedakaPolishing
#
## Wall time limit:
#SBATCH --time=24:00:00
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



```

</details>

We can submit it by:

```bash
cd /cluster/projects/nn9987k/$USER
sbatch /cluster/projects/nn9987k/.scripts/3_Medaka.SLURM.sh D01T6_T /cluster/projects/nn9987k/$USER/results/MetaG/D01T6_T_Chopper/D01T6_T_Chopper /cluster/projects/nn9987k/$USER/results/MetaG /cluster/projects/nn9987k/$USER/results/MetaG/D01T6_T.MEDAKA.dir && mkdir /cluster/projects/nn9987k/$USER/results/MetaG/D01T6_T.MEDAKA.dir
```

> [!NOTE]
> Unfortunatelly for the BIN240 Course Sigma2 has only assigned a copuple of nodes in the FRAM computer, so it is most likely the Job never runs/finish.


Same as the fly command, we have already prepared the data for the course and you can copy the Medaka results by:

> [!WARNING] 
> !NB: Remember to kill the Fly job by scancel ```<JOBID>``` Before copying the data .

```bash
rsync -avLh /cluster/projects/nn9987k/.results/MetaG/D01T6_T.MEDAKA.dir/D01T6_T.medaka.dir/D01T6_T.medaka.consensus.fasta /cluster/projects/nn9987k/$USER/results/MetaG/D01T6_T.MEDAKA.dir/
tree /cluster/projects/nn9987k/$USER/results/MetaG/D01T6_T.MEDAKA.dir/
```

```
/cluster/projects/nn9987k/auve/results/MetaG/D01T6_T.MEDAKA.dir/
└── D01T6_T.medaka.consensus.fasta

0 directories, 1 file
```

### Comparing Assemblies before and after polishing:

The best way to compare the assemblies is to perform ```assembly-stats``` on both and then compare the results.

Let's do it.

> [!Important]
> As we will start working with files let's ask for an interactive session in FRAM:

```bash
srun \
--account=nn9987k \
--gres=localscratch:20G \
--cpus-per-task 4 \
--nodes 1 \
--time=02:00:00 \
--pty bash \
-i
```

Now we can run ```assemlby-stats``` on both files like this:

```bash
module load Miniconda3/23.10.0-1
conda activate /cluster/projects/nn9987k/.share/conda_environments/MetaG_Assembly_And_Binning/
FLYE="/cluster/projects/nn9987k/$USER/results/MetaG/D01T6_T.flye.outdir/assembly.fasta"
MEDAKA="/cluster/projects/nn9987k/$USER/results/MetaG/D01T6_T.MEDAKA.dir/D01T6_T.medaka.consensus.fasta"
assembly-stats $FLYE $MEDAKA
```

```
stats for /cluster/projects/nn9987k/auve/results/MetaG/D01T6_T.flye.outdir/assembly.fasta
sum = 821745693, n = 45282, ave = 18147.29, largest = 1393968
N50 = 25699, n = 8036
N60 = 19600, n = 11713
N70 = 15449, n = 16442
N80 = 11905, n = 22498
N90 = 8290, n = 30706
N100 = 63, n = 45282
N_count = 0
Gaps = 0
-------------------------------------------------------------------------------
stats for /cluster/projects/nn9987k/auve/results/MetaG/D01T6_T.MEDAKA.dir/D01T6_T.medaka.consensus.fasta
sum = 818830768, n = 45282, ave = 18082.92, largest = 1393760
N50 = 25629, n = 8025
N60 = 19526, n = 11700
N70 = 15391, n = 16430
N80 = 11860, n = 22487
N90 = 8258, n = 30698
N100 = 63, n = 45282
N_count = 0
Gaps = 0

```

These results are good but not very comparable ```assembly-stats``` can perform an output as a table using the ```-t``` flag:

```
assembly-stats -t $FLYE $MEDAKA

```

We can save this into a file:

```bash
assembly-stats -t $FLYE $MEDAKA > /cluster/projects/nn9987k/$USER/results/MetaG/Flye.Medaka.stats.tsv
cat !$
```

```
cat /cluster/projects/nn9987k/$USER/results/Flye.Medaka.stats.tsv
filename        total_length    number  mean_length     longest shortest        N_count Gaps    N50     N50n    N70     N70n    N90     N90n
/cluster/projects/nn9987k/auve/results/MetaG/D01T6_T.flye.outdir/assembly.fasta 821745693       45282   18147.29        1393968 63      0       0       25699   8036    15449   16442   829030706
/cluster/projects/nn9987k/auve/results/MetaG/D01T6_T.MEDAKA.dir/D01T6_T.medaka.consensus.fasta  818830768       45282   18082.92        1393760 63      0       0       25629   8025    15391       16430   8258    30698
```

We can plot this result table:

```bash
conda activate /cluster/projects/nn9987k/.share/conda_environments/R_env/
cd /cluster/projects/nn9987k/$USER/results/MetaG/
Rscript /cluster/projects/nn9987k/.scripts/assemblyStats.r /cluster/projects/nn9987k/$USER/results/MetaG/Flye.Medaka.stats.tsv
```

This producess the plot:

```
Plot saved as AssemblyStats.pdf
```

![AssemblyStats](https://github.com/TheMEMOLab/Bin420-Bioinformatics-for-Functional-Meta-Omics/blob/main/img/Assemblystats.PNG)

**What can we say about these results?**

> [!Important]
> Remember to finish your interactive session by ```exit```

## 4. Binning.

So far we created an assembly containing all contigs (continuous sequences) from each of the organisms in the microbial community that we sequenced from the cow rumen.

Presently, the assembly consists of thousands of contigs, each coming from a single species. By grouping together the contigs from each species present in our sample, we can create what is referred to as a MAG, short for Metagenome Assembled Genome.

Popular binning algorithms like the ones used in Metabat2 utilize contig depth as a tell tale to link the individual contigs together that come from the same species. This is done by mapping the original reads onto the assembly and then counting the read depth of each contig. The smart thing here is that contigs coming from the same species will have similar depth. Another vital contig statistic that binners use is the GC-content. Each species has its own intrinsic GC-content, and by grouping contigs further on GC-content -in this case by counting the tetranucleotide frequency- we might get good representatives for distinct species in our sample. If our bins live up to our requirements, we can refer to them as MAGs.