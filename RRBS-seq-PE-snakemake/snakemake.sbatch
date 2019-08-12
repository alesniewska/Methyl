#!/bin/bash

#SBATCH -p standard 
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=2-00:00:00

# load all modules
module load python/3.7.3
module load bedtools/2.21.0
module load bedtools/2.21.0/intersectBed
module load bedtools/2.21.0/groupBy
module load trimgalore/0.6.0
module load fastqc
module load bismark/0.21.0
module load bismark/0.21.0/bismark_genome_preparation
module load bismark/0.21.0/bismark_methylation_extractor
module load samtools
module load bowtie2

# activate virtual environment
python -m venv venv
source venv/bin/activate

# install missing packages
pip install pandas numpy snakemake biopython

# run a snakemake pipe
snakemake --unlock
snakemake --latency-wait 90 --cluster "sbatch -p standard --mem=200000 -N 1 -n 4 --time=2-00:00:00" --jobs 250 --jobname "{rulename}.{jobid}"

