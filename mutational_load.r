#!/usr/bin/env Rscript

#Record non-syn and syn mutations and nucleotide diversity using PopGenome

library(PopGenome)

argv<-commandArgs(TRUE)
chr<-argv[1]
size<-argv[2]
size2<-argv[3]
inv<-argv[4]
haplo<-argv[5]
vcf_name<-argv[6]

#Define populations
a<-read.csv('inversion_genotypes.csv')
refs<-as.character(a2$sample_ID[a2[,inv]=='0/0'])
invs<-as.character(a2$sample_ID[a2[,inv]=='1/1'])
all<-as.character(a2$sample_ID)

print(refs)
print(invs)
print(all)

#Read in data
vcf<-paste(vcf_name,'_',chr,'_SNPs.vcf.gz',sep='')
gff<-paste('gff/',chr,'.gff',sep='')
fasta<-paste('fasta/',chr,'.fasta',sep='')

#perform analyses in 5-Mb windows to reduce runtime and memory requirements
if(haplo=='inv'){
smpls<-invs
}
if(haplo=='ref'){
smpls<-refs
}

start<-size
end<-size2
print(start)
print(end)
s0<-readVCF(vcf,numcols=10000,tid=chr,from=start,to=end,gffpath=gff,include.unknown=F,samplenames=smpls)  #use include.unknown=F to exclude positions with unknown/missing data

s<-s0

print('indivs of s')
print(get.individuals(s))

#verify syn-nonsyn sites
s<-set.synnonsyn(s,ref.chr=fasta)

#to calculate nuc. diversity for each site, use sliding window=1 SNP
s.slide<-sliding.window.transform(s,1,1,type=1) #sliding window by 1 SNP, type=1 uses biallelic snps only

#calculate nucleotide diversity (pi) for nonsyn sites
s.slide<-diversity.stats(s.slide,subsites='nonsyn')
pN<-s.slide@nuc.diversity.within

#calculate nucleotide diversity (pi) for syn sites
s.slide<-diversity.stats(s.slide,subsites='syn')
pS<-s.slide@nuc.diversity.within

#Save pi results
df<-data.frame(pN=pN,pS=pS,pos=s@region.data@biallelic.sites[[1]])
filename = paste('piN_piS_',inv,'_',haplo,'_',start,'.csv',sep='')
write.csv(df,file=filename)