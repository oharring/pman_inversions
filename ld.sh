bcftools view -s SAMPLE_LIST -r CHROM INPUT_VCF.vcf.gz | vcftools --vcf - --maf 0.05 --thin 1000 -c --geno-r2 --max-missing-count 0 | perl emerald2windowldcounts_osh.pl | gzip > OUTPUT.ld.txt.gz
