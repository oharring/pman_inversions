#mask repeats with repeatmasker

module load RepeatMasker/4.0.5-fasrc05
module load samtools
sbatch -J mask.submit -e mask.err -o mask.out -n 1 -N 1 --mem 10000 -t 0-24:00 -p hoekstra,shared --wrap="RepeatMasker --species rodentia --xsmall breakpoints_all.fasta -dir breakpoints/repeatmasker"

#run sedef
module load GCC/8.2.0-2.31.1
module load parallel/20180522-fasrc01
module load samtools

./sedef.sh -o sedef_all -j 10 breakpoints/repeatmasker/breakpoints_all.fasta.masked
