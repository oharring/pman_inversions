#!/usr/bin/env python
# coding: utf-8

import scipy
import allel; print(allel.__version__)
import numpy as np
import os
import zarr
import numcodecs
import sys
import matplotlib.cm as cm
import pandas as pd
import csv


#For each inversion, use informative populations to perform overall PCA and then project all samples onto top PC axes
#Loading samples onto top PC axes facilitates genotyping individual samples

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

#fill in PacBio samples pop in pop map dataframe
#note that PacBio samples are labeled as SAMPLE with number, corresponding to order in which merged
pop_map[pop_map['pop'].isna()] #view rows with missing pop
pop_map.at[pop_map.index[pop_map['sample_ID']=='SAMPLE'],'pop']='SW0'
pop_map.at[pop_map.index[pop_map['sample_ID']=='3:SAMPLE'],'pop']='BK0'
pop_map.at[pop_map.index[pop_map['sample_ID']=='4:SAMPLE'],'pop']='BW0'
pop_map.at[pop_map.index[pop_map['sample_ID']=='5:SAMPLE'],'pop']='PO0'
pop_map.at[pop_map.index[pop_map['sample_ID']=='6:SAMPLE'],'pop']='NUB0'

#create dictionary for sample_ID names per population 
inds_per_pop = {pop:list(samples) for pop,samples in pop_map.groupby('pop')['sample_ID']}

#create dictionary for sample_ID indices per population
idx_per_pop = {pop:list(index) for pop,index in pop_map.groupby('pop')['index']}

set(pop_map['pop'])

#Genotypes

#set up dictionaries for positions and genotypes 
positions = {}
genotypes = {}

#make dictionaries
for chrom in list(chromosomes.keys()):
    print(chrom)
    positions[chrom] = allel.SortedIndex(callset[chrom]['variants']['POS'])
    genotypes[chrom] = allel.GenotypeChunkedArray(callset[chrom]['calldata']['GT'])

#Calculate heterozygosity per sample

def percent_heterozygosity(chrom, start=0, stop=int(chromosomes[chrom]), pops='all'):
   
   loc = positions[chrom].locate_range(start, stop)
   
   # if all pops are specified, keep indices and pull out all genotypes in region
   if pops == 'all':
       new_idx_per_pop = idx_per_pop
       indices_subset = np.array([sample for sublist in idx_per_pop.values() for sample in sublist])
       genos = genotypes[chrom][loc]
   
   # if list of pops specified
   elif type(pops) is list:
       # make sure that list is sorted alphabetically (same order as full list)
       pops.sort()
       # pull out indices for subset and take genotypes
       idx_per_pop_subset = {pop:idx for (pop,idx) in idx_per_pop.items() if pop in pops}
       indices_subset = np.array([sample for sublist in idx_per_pop_subset.values() for sample in sublist])
       genos = genotypes[chrom][loc].take(indices_subset,axis=1)
       
       # create new pop_map based on number of samples and assuming alphabetical order
       new_idx_per_pop = {}
       last_idx = 0

       for pop in pops:
           new_idx_per_pop[pop] = list(range(last_idx, len(idx_per_pop_subset[pop]) + last_idx))
           last_idx += len(idx_per_pop_subset[pop])
       
   else:
       raise Exception("ERROR: pops should be 'all' or a list of pop_IDs")
   
   positions2 = positions[chrom].tolist()
   
   # subset to region of interest
   positions2 = [x for x in positions2 if x<stop and x>=start] 

   # calculate heterozygosity (percent heterozygous)
   pc_het = genos.count_het(axis=0)[:] * 100 / (stop - start)
   
   # save results
   d = {'original_sample_index': indices_subset, 'percent_het': pc_het}
   df = pd.DataFrame(data=d)
   pop_map2 = pop_map[pop_map['index'].isin(indices_subset.tolist())]
   df2=df.merge(pop_map2, left_on='original_sample_index', right_on='index')
   return(df2)


#PCA

def plot_PCA(chrom,start=0, stop=int(chromosomes[chrom]), pops='all', geno=0.8, thinning=50,
            pcX=1, pcY=2, scree_plot=True, plot_het=True, score_df=False, save_plots=False):

    loc = positions[chrom].locate_range(start, stop)

    #RUN PCA for populations with inversion                                                                                                                                                     
    #populations for performing PCA are listed as pops                                                                                                                                         
    # make sure that list is sorted alphabetically (same order as full list)                                                                                                                   
 
    pops.sort()

    # pull out indices for subset and slice / take genotypes                                                                                                                                   
 
    idx_per_pop_subset = {pop:idx for (pop,idx) in idx_per_pop.items() if pop in pops}
    indices_subset = np.array([sample for sublist in idx_per_pop_subset.values() for sample in sublist])
    genos = genotypes[chrom][loc].take(indices_subset,axis=1)

    new_idx_per_pop = {}
    last_idx = 0
    for pop in pops:
        new_idx_per_pop[pop] = list(range(last_idx, len(idx_per_pop_subset[pop]) + last_idx))
        last_idx += len(idx_per_pop_subset[pop])

    # get allele counts for relevant pops                                                                                                                                                       
    acs_all = genos.count_alleles()

    # define min. number of chromosomes / non-missing data                                                                                                                                      
    min_chroms = int(len(indices_subset) * 2 * geno)

    # filter and thin 
    pca_selection = (acs_all.max_allele() == 1) & (acs_all[:, :2].min(axis=1) >= 2) & (acs_all.sum(axis=1) >= min_chroms)
    indices = np.nonzero(pca_selection)[0]
    thinned = indices[::thinning]
    print("Retained SNPs:",len(thinned))

    # apply filter and get genotypes                                                                                                                                                            
    genotypes_pca = genos.take(thinned, axis=0)
    gn = genotypes_pca.to_n_alt()[:]

    # do PCA                                                                                                                                                                                    
    coords, model = allel.pca(gn, n_components=10, copy=True, scaler='patterson', ploidy=2)

    #PROJECT ALL SAMPLES ONTO PC axes                                                                                                                                                           

    new_idx_per_pop2 = idx_per_pop
    indices_subset2 = np.array([sample for sublist in idx_per_pop.values() for sample in sublist])
    genos2 = genotypes[chrom][loc].take(indices_subset2,axis=1)

    # apply filter and get genotypes for ALL samples                                                                                                                                            
    genotypes_pca2 = genos2.take(thinned, axis=0)
    gn2 = genotypes_pca2.to_n_alt()[:]

    # project onto PCA using transform                                                                                                                                                          
    coords2 = model.transform(gn2)

    # record PC1,PC2 scores and heterozygosity into dataframe                                                                                                                                   
    if score_df==True:
        heterozyg = percent_heterozygosity(chrom,start,stop,'all')
        x = coords2[:, pcX-1]
        y = coords2[:, pcY-1]
        pc1_pve = [model.explained_variance_ratio_[pcX-1]*100]*len(x)
        pc2_pve = [model.explained_variance_ratio_[pcY-1]*100]*len(y)
        d = {'original_sample_index': indices_subset2, 'pc1_score': x, 'pc2_score': y, 'pc1_pve': pc1_pve, 'pc2_pve':pc2_pve}
        df = pd.DataFrame(data=d)
        df2=df.merge(heterozyg, left_on='original_sample_index', right_on='index')
        df2.to_csv(plots_path+'projected_PC1_score_{}_'.format(chrom)+str(start)+'_'+str(stop)+'_'+'_'.join(pops)+'.csv',index=False)



#Record PCA and heterozygosity for each inversion

#read in inversions
haps = pd.read_csv(base_path+'inversion_intervals.csv', sep=",")

#record PCA and heterozygosity

for i in [x for x in haps.index.values]:
    chroms = haps['chr'][i]
    starts = int(haps['start'][i])
    ends = int(haps['end'][i])
    pop = haps['pop'][i] #pre-define which populations to use for defining PCs per inversion
    popps = pop.split('_')
    thin=1
    PCA(chroms,start=starts, stop=ends, pops=popps, geno=0.8, thinning=thin,
            pcX=1, pcY=2, score_df=True)
        

