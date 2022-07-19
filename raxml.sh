#select samples from VCF and perform initial thinning

bcftools view -m2 -M2 -i 'F_MISSING<0.05' -s SAMPLE_LIST -r CHROM CHROM_SNPs.vcf.gz | vcftools --vcf - --thin 1000 --recode --recode-INFO-all --out CHROM_SNPs_thinned.vcf

#merge vcfs across all chromosomes

bcftools concat chr1_SNPs_thinned.recode.vcf chr2_SNPs_thinned.recode.vcf chr3_SNPs_thinned.recode.vcf chr4_SNPs_thinned.recode.vcf chr5_SNPs_thinned.recode.vcf chr6_SNPs_thinned.recode.vcf chr7_SNPs_thinned.recode.vcf chr8_SNPs_thinned.recode.vcf chr9_SNPs_thinned.recode.vcf chr10_SNPs_thinned.recode.vcf chr11_SNPs_thinned.recode.vcf chr12_SNPs_thinned.recode.vcf chr13_SNPs_thinned.recode.vcf chr14_SNPs_thinned.recode.vcf chr15_SNPs_thinned.recode.vcf chr16_SNPs_thinned.recode.vcf chr17_SNPs_thinned.recode.vcf chr18_SNPs_thinned.recode.vcf chr19_SNPs_thinned.recode.vcf chr20_SNPs_thinned.recode.vcf chr21_SNPs_thinned.recode.vcf chr22_SNPs_thinned.recode.vcf chr23_SNPs_thinned.recode.vcf -o all_chr.vcf.recode.vcf

#thin to 1 SNP per 100kb

vcftools --vcf all_chr.vcf.recode.vcf --thin 100000 --recode --recode-INFO-all --out all_chr_thinned_100kb.vcf

#convert vcf to PHYLIP (alignment) matrix
#https://github.com/edgardomortiz/vcf2phylip

module load python/3.6.3-fasrc01
python vcf2phylip.py -i all_chr_thinned_100kb.vcf.recode.vcf

#remove invariant sites
#https://github.com/btmartin721/raxml_ascbias

python3 ascbias.py -p all_chr_thinned_100kb.min4.recode.min4.phy -o all_chr_thinned_100kb.min4_no_invariants.phy

#make tree with RAxML

sbatch -J bootstrap100.submit -e bootstrap100.err -o bootstrap100.out -n 1 -N 1 --mem 10000 -t 0-48:00 -p hoekstra,shared,commons --wrap="raxml/standard-RAxML-master/raxmlHPC -m ASC_GTRCAT --asc-corr=lewis -f a -p 12489 -x 12489 -# 100 -s all_chr_thinned_100kb.min4_no_invariants.phy -n RAxML_wholegenome_thinned_100kb_100bootstrap"


