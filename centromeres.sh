#module load GCC/8.2.0-2.31.1 OpenMPI/3.1.3 BLAST+/2.9.0
#module load blast/2.2.29+-fasrc03

#downloaded Pman1 satellite sequence (pman1_sat.fasta) associated with centromeres (from Smalec et al 2019), NCBI GenBank accession number KX555281
#note that I only use one satellite sequence since all P. maniculatus centromere satellite sequences have sequence identity > 95%

#created database from BW flye assembly for blasting
makeblastdb -in flye_assembly/BW/assembly.fasta -dbtype nucl -title BW_flye_db -out centromere/data/blast/BW_flye_db

#blast satellite sequences against BW flye assembly
blastn -query pman1_sat.fasta -db centromere/data/blast/BW_flye_db -out centromere/data/blast/pman1_sat_BW_flye.csv -outfmt 7

#repeat for each genome
