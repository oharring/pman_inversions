#!/usr/bin/env Rscript

#run local PCA with lostruct

argv<-commandArgs(TRUE)
chrom<-argv[1]
pop<-argv[2]
bcf_header<-argv[3]
output_header<-argv[4]

library(data.table)
library(lostruct)

bcf<-paste(paste(bcf_header,chrom,'SNPs',pop,sep='_'),'.bcf',sep='')

snps<-vcf_windower(bcf,size=1e5,type='bp')
pcs<-eigen_windows(snps,k=2)
pcdist<-pc_dist(pcs,npc=2)

#remove nas
nas <- is.na(pcdist[,1])

#MDS 
fit<-cmdscale(pcdist[!nas,!nas],eig=TRUE,k=2)

all<-fit$points

output_fl<-paste(output_header,'mds_',chrom,'_',pop,'_wins100kb.rds',sep='')
saveRDS(all,file=output_fl)

output_fl2<-paste(output_header,'pcdist_',chrom,'_',pop,'_wins100kb.rds',sep='')
saveRDS(pcdist,file=output_fl2)

output_fl3<-paste(output_header,'pcs_',chrom,'_',pop,'_wins100kb.rds',sep='')
saveRDS(pcs,file=output_fl3)
