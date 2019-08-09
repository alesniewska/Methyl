#!/bin/bash

#SBATCH -p standard 
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=2-00:00:00

module load bedtools/2.21.0
module load bedtools/2.21.0/intersectBed
module load bedtools/2.21.0/groupBy
module load python/3.7.3
module load trimgalore/3.7.3
module load fastqc
module load bismark/0.21.0
module load bismark/0.21.0/bismark_genome_preparation
module load bismark/0.21.0/bismark_methylation_extractor
module load bowtie2

python -m venv venv
source venv/bin/activate
pip install pandas numpy biopython snakemake

snakemake --unlock
snakemake --latency-wait 180 --cluster "sbatch -p standard --mem=200000 -N 1 -n 4 --time=2-00:00:00" --jobs 20 --jobname "{rulename}.{jobid}"
