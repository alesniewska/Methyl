#!/bin/bash

#SBATCH -p standard
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=24:00:00

module load python/3.7.3
module load trimgalore/0.6.0
module load fastqc
module load star
module load samtools
module load subread

python -m venv venv
source venv/bin/activate
pip install pandas numpy snakemake

snakemake --unlock
snakemake --latency-wait 60 --cluster "sbatch -p standard --mem=200000 -N 1 -n 4 --time=24:00:00 " --jobs 150 --jobname "{rulename}.{jobid}"
