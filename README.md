# RNA-seq and RRBS-seq snakemake pipelines

[Sanakemake](https://bitbucket.org/snakemake/snakemake/src/master/) is a workflow management system which allows to create analysis pipelines.
Snakemake pipelines were written for analyzing RNA-seq and RRBS-seq data for paired-end and single-end reads 
Workflows were created to be able to run on a cluster with [SLURM](https://slurm.schedmd.com/quickstart.html) management system.

## Installation

```
git clone https://github.com/alesniewska/Methyl.git
```
## Configuration

Before running the pipeline, change the paths to the directories in the right place in the file Snakefile.

For example:

```bash
#---------------change to match your configutation-------------------
DIR = "/home/plgrid/plgkraszewska/RNA_PE/"
SAMPLE_PATH = DIR +"sample_rna/"
GENOME_PATH = DIR + "genome/GRCh37_75.fa"
GENOME_ANNOTATION_PATH = DIR + "genome_ann/GRCh37.75.gtf"
SAMPLES, = glob_wildcards(SAMPLE_PATH + "{sample}_1.fastq")
#--------------------------------------------------------------------
```
You need to load the appropriate modules in the file snakemake.sbatch. 

```bash
# load all modules
module load python/3.7.3
module load trimgalore/0.6.0
module load fastqc/0.11.4
module load star/2.5.3a
module load samtools/1.2
module load subread/1.5.2
```
Make sure your cluster contains the required modules. If your cluster doesn't contain all needed modules you have to install missing packages and set path in right place in Snakefile rules.
All system requirements have been saved in the requirements.txt file in the appropriate directories.

## Running snakemake pipeline

Choose the directory that contains all required Python scripts the configured snakemake.sbatch and Snakefile files run:
```
sbatch snakemake.sbatch
```

## RNA-seq-snamemake pipeline schema

![alt text](https://github.com/alesniewska/Methyl/blob/master/images/rna_shema.jpg?raw=true)

## RRBS-seq-snamemake pipeline schema

![alt text](https://github.com/alesniewska/Methyl/blob/master/images/rrbs-schema.jpg?raw=true)

## Make table count

After running pipeline you get result files which contain Geneid, chromosome, start, end, strand, length, and numer of mapping reads (RNA-seq) or percentage of methylated CpGs in each island (RRBS-seq).
If you want to get a csv file, which contains only  Geneid and numer of mapping reads or percentage of methylated CpGs you can use the table_count.py

```
python table_count.py <path to *.cnt files> <output files>
```


# CpGislandsFind

You can use CpGislandsFind python script to find CpGislands in any DNA sequence in FASTA format. 

## Running CpGislandsFind

```
CpGislandFind.py <path to FASTA file> <optional arguments>
```
## Arguments

* Fasta file with DNA sequence
* minimal CpGisland length; default 500
* minimal Obs to Exp value; default 0.6
* minimal sum of C percentage and G percentage in island; default 50


## How does CpGislandsFind work?

1. Find distance between CpG in sequence and calculate median of distances
2. If distance between CpG is shorter than calculated median distance add CpG to pre-island
3. Select pre-island which is biggest or equal than minimal CpGisland length 
4. If sum of C% and G% is biggest than minimal value and Obs/Exp is biggest than minimal value write CpGisland

![alt text](https://github.com/alesniewska/Methyl/blob/master/images/obs_exp.jpg?raw=true)














