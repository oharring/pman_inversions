#!/bin/bash
#SBATCH -n 1                # Number of cores
#SBATCH -N 1
#SBATCH -t 0-12:00
#SBATCH -p hoekstra,shared   # Partition to submit to
#SBATCH --mem=80000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o XXX_dotplotly.out  
#SBATCH -e XXX_dotplotly.err  
#SBATCH -J XXX_dotplotly.submit

POP=$1
CONTIG=$2
WORDSIZE=$3
POP2=$4
CONTIG2=$5

echo $POP
echo $CONTIG
echo $WORDSIZE
echo $POP2
echo $CONTIG2

nucmer --maxmatch --nosimplify -l $WORDSIZE -c 100 $POP-$CONTIG.fasta $POP2-$CONTIG2.fasta -p $POP-$CONTIG-$POP2-$CONTIG2.nucmer

show-coords -c $POP-$CONTIG-$POP2-$CONTIG2.nucmer.delta > $POP-$CONTIG-$POP2-$CONTIG2.nucmer.delta.coords

