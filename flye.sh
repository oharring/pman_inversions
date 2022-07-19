#!/bin/bash
#SBATCH -n 12                # Number of cores
#SBATCH -N 1
#SBATCH -t 30-0:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p bigmem          # Partition to submit to
#SBATCH --mem=499G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o BK_flye.out  
#SBATCH -e BK_flye.err  
#SBATCH -J BK_flye.submit

#--pacbio-raw: indicates uncorrected pacbio reads. Options for corrected reads.
#-o: path where assembly will be output
#--threads: number of threads used for assembly
#--asm-coverage is coverage used for disjointig assembly
#--genome-size needed for asm-coverage
flye --pacbio-raw BK.fastq.gz -o output_path/BK --threads 12 --asm-coverage 40 --genome-size 2.7g
