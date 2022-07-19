#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scipy
import allel; print(allel.__version__)
import numpy as np
import os
import zarr
import numcodecs
import sys
import pandas as pd
import csv
import random


#Set-up                                                                                                                                                                             
 
#Set directory paths                                                                                                                                                                

base_path = '' #set base path for files                                                                                                                                             
genome_path = '' #set path for genome                                                                                                                                               
refseq = os.path.join(genome_path,'Pman2.1.3.chromosomes-unplaced_scaffolds.fasta')
popmap_file = 'sample_pop_map.txt'
ouput_path = '' #set path for outputs                                                                                                                                               
 

#Set up input chromosomes                                                                                                                                                           
 
chromosomes = dict()
with open(base_path+'population_genomics/scikit_allel/vcfs/Pman2.1.3.chromosomes-unplaced_scaffolds.chrom.sizes') as csvfile:
    f = csv.reader(csvfile,delimiter='\t')
    for row in f:
            chromosomes[row[0]] = row[1]

chromosomes.keys()


#Convert vcfs to zarr                                                                                                                                                               
#set zarr directory                                                                                                                                                                 
zarr_dir = base_path + 'scikit_allel/zarr_vcfs'

#set vcf directory                                                                                                                                                                  
in_vcf_dir = base_path + 'vcfs/'

#convert vcf to zarr for each chromosome into same zarr directory                                                                                                                   
 
for chrom in chromosomes.keys():
    print("processing:",chrom)
    in_vcf = in_vcf_dir+'hf.'+chrom+'biallelic_snps.vcf.gz'
    allel.vcf_to_zarr(in_vcf, zarr_dir, fields='*', overwrite=False,
                      group='{}'.format(chrom), region='{}'.format(chrom))

callset = zarr.open_group(zarr_dir, mode='r')

#Population and sample data                                                                                                                                                         
#make dataframe of samples and indices                                                                                                                                              
 
samples = pd.DataFrame(callset['chr1']['samples'], columns=['sample_ID'])

#add population information from popmap file                                                                                                                                        
pop_table = pd.read_csv(popmap_file, sep="\t", header=None, names=['sample_ID','pop'])
pop_map = pd.merge(samples,pop_table,how='outer').reset_index()

# create dictionary for sample_ID names per population 
inds_per_pop = {pop:list(samples) for pop,samples in pop_map.groupby('pop')['sample_ID']}

# create dictionary for sample_ID indices per population
idx_per_pop = {pop:list(index) for pop,index in pop_map.groupby('pop')['index']}

set(pop_map['pop'])

#Allele counts

# set up dictionaries for positions, genotypes and population allele counts
positions = {}
genotypes = {}
acs = {}
# make dictionaries 
for chrom in list(chromosomes.keys()):
    print(chrom)
    positions[chrom] = allel.SortedIndex(callset[chrom]['variants']['POS'])
    genotypes[chrom] = allel.GenotypeChunkedArray(callset[chrom]['calldata']['GT'])
    # calculate allele counts for each pop
    #acs[chrom] = genotypes[chrom].count_alleles_subpops(idx_per_pop2)

#Choose homozygous samples for inversion and standard haplotypes

#Read in genotypes file

genos = base_path+'inversion_genotypes.csv'
invs = pd.read_csv(genos)

#Record which populations to use for each inversion
invs_dict = {'chr15.0':['SW','BK'],'chr14.0':['SW','BK'],
             'chr13.0':['SW','BK'],'chr19.0':['SW','BK'],'chr21.0':['BW'],
             'chr22.0':['SW','BK'],'chr6.0':['BK','SW'],'chr7.0':['BK'],'chr7.1':['PO'],'chr7.2':['SW','BK'],
             'chr7.3':['BK','BW'],'chr3.0':['SW','BK'],'chr9.0':['SW','BK'],'chr9.1':['BW','SW'],
             'chr18.0':['SW','BK'],'chr23.0':['PO'],'chr20.0':['BK','SW'],'chr15.1':['BW'],'chr15.2':['BK','SW'],
             'chr10.0':['BW','NUB'],'chr11.0':['BW','NUB']}

samples_to_use={}
for inv in invs_dict.keys():
    samples_to_use[inv]={}
    samples_to_use[inv]['hom0']=[]
    samples_to_use[inv]['hom1']=[]
    hom0=invs[invs[inv]=='0/0']
    for pop in invs_dict[inv]:
        inds = [x for x in hom0['sample_ID'][hom0['pop']==pop]]
        samples_to_use[inv]['hom0']+=inds
    hom1=invs[invs[inv]=='1/1']
    for pop in invs_dict[inv]:
        inds = [x for x in hom1['sample_ID'][hom1['pop']==pop]]
        samples_to_use[inv]['hom1']+=inds


#Save chosen samples
inv_list = []
hom0_list = []
hom1_list = []
for inv in samples_to_use.keys():
    inv_list.append(inv)
    hom0_list.append(','.join(samples_to_use[inv]['hom0']))
    hom1_list.append(','.join(samples_to_use[inv]['hom1']))


#Allele counts by inversion genotype

def compute_acs(uniq_ID, chrom):
        
    #get sample IDs for hom0 and hom1 individuals
    hom0 = samples_to_use[uniq_ID]['hom0']
    hom1 = samples_to_use[uniq_ID]['hom1']
        
    idx_0 = [mm for mm in pop_map['index'][pop_map['sample_ID'].isin(hom0)]]
    idx_1 = [mm for mm in pop_map['index'][pop_map['sample_ID'].isin(hom1)]]
        
    idx_per_geno = {'hom0':idx_0, 'hom1':idx_1}
    
    #return allele counts
    acs[uniq_ID] = genotypes[chrom].count_alleles_subpops(idx_per_geno)

    #return sample IDs
    return({'hom0':hom0,'hom1':hom1})



#Fst

def fst(pop1, pop2, chrom, uniq_ID, winsize=10000, 
             save_fst=True, return_fst=True):
    
    start = 0
    stop = int(chromosomes[chrom])
    # get positions and allele count for haplotype
    variants = positions[chrom]
    variants_subset = np.logical_and(variants > start, variants < stop) #subset variants to relevant region
    
    # create combined allele count array and create filter
    acu = allel.AlleleCountsArray(acs[uniq_ID][pop1].compress(variants_subset, axis=0)[:] + acs[uniq_ID][pop2].compress(variants_subset, axis=0)[:])
    flt = acu.is_segregating()
    print('retaining', np.count_nonzero(flt), 'SNPs')
    
    # apply filter to positions and allele count arrays
    variants2 = variants.compress(variants_subset)
    pos = variants2.compress(flt)
    ac1_0 = acs[uniq_ID][pop1].compress(variants_subset, axis=0)
    ac2_0 = acs[uniq_ID][pop2].compress(variants_subset, axis=0)
    ac1 = ac1_0.compress(flt, axis=0)[:, :2]
    ac2 = ac2_0.compress(flt, axis=0)[:, :2]

    # run windowed fst
    fst, windows, counts = allel.windowed_hudson_fst(pos, ac1, ac2, size=winsize)
    y = fst
    # record mid-block position
    x = windows[:,0:1]+(winsize/2)
    
    # record fst into dataframe
    if save_fst==True:
        d = {'window_position': [dd[0] for dd in x.tolist()], 'fst': y.tolist(), 'chr': [chrom]*len(y)}
        df = pd.DataFrame(data=d)
        df.to_csv(plots_path+'fst_{}_'.format(uniq_ID)+pop1+'_'+pop2+'.csv',index=False)
        
    # return list of fst values for all windows
    if return_fst==True:
        return(y.tolist())


#Fst for each inversion

#Compute allele counts

acs = {}
for inv in invs_dict.keys():
    print(inv)
    chrm = inv.split('.')[0]
    aa=compute_acs(inv,chrm)


#Compute Fst
for i in invs_dict.keys():
    chrm = i.split('.')[0]
    fst('hom0', 'hom1', chrm, i, winsize=10000, 
            save_fst=True, return_fst=False)


