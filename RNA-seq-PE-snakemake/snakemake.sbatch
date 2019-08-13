#!/bin/bash

#SBATCH -p standard
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=24:00:00

# load all modules
module load python/3.7.3
module load trimgalore/0.6.0
module load fastqc/0.11.4
module load star/2.5.3a
module load samtools/1.2
module load subread/1.5.2

# activate virtual environment
python -m venv venv
source venv/bin/activate

# install missing packages
pip install pandas numpy snakemake

# run a snakemake pipe
snakemake --unlock
snakemake --latency-wait 60 --cluster "sbatch -p standard --mem=200000 -N 1 -n 4 --time=24:00:00 " --jobs 20 --jobname "{rulename}.{jobid}"
