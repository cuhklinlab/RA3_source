#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import numpy as np
import scanpy as sc
import os
import time
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import adjusted_mutual_info_score
from sklearn.metrics.cluster import homogeneity_score
from sklearn.metrics.cluster import normalized_mutual_info_score
from scipy import stats
# import rpy2.robjects as robjects
# from rpy2.robjects import pandas2ri
import hdf5storage
import re

workdir = './clustering/'
path_clusters = os.path.join(workdir,'clusters/')
path_metrics = os.path.join(workdir,'metrics/')
path_label = './data/' # this path contains count matrix and label for the data
path_mat = './output_RA3/' # this path contains RA3 outputs (.mat file format)


def getNClusters(adata,n_cluster,range_min=0,range_max=3,max_steps=20):
    this_step = 0
    this_min = float(range_min)
    this_max = float(range_max)
    while this_step < max_steps:
        print('step ' + str(this_step))
        this_resolution = this_min + ((this_max-this_min)/2)
        sc.tl.louvain(adata,resolution=this_resolution)
        this_clusters = adata.obs['louvain'].nunique()

        print('got ' + str(this_clusters) + ' at resolution ' + str(this_resolution))

        if this_clusters > n_cluster:
            this_max = this_resolution
        elif this_clusters < n_cluster:
            this_min = this_resolution
        else:
            return(this_resolution, adata)
        this_step += 1

    print('Cannot find the number of clusters')
    print('Clustering solution from last iteration is used:' + str(this_clusters) + ' at resolution ' + str(this_resolution))


df_metrics = pd.DataFrame(columns=["DR methods", "Clustering methods", "peak selection k", "Dimension of latent space", "Number of peaks", "Number of cells", "ARI", "AMI", "Homogeneity", "NMI", "Clustering path", "Dataset"])
datanames = ['donor_BM0828']
# cluster2 = ['GMvsHL','GMvsHek']
for i in range(len(datanames)):
    dataname = datanames[i]
    files = [x for x in os.listdir(path_mat) if re.match('Our_'+dataname+'_peak.*',x)]
    label = pd.read_csv(os.path.join(path_label,dataname+'_label.txt'),sep='\t',header=None)
    # if dataname in cluster2:
        # num_clusters=2
    # else:
        # num_clusters = len(np.unique(label))
    num_clusters = len(np.unique(label))
    for file in files:
        method = file.split('.')[0]
        print(method)
        res_mat = hdf5storage.loadmat(os.path.join(path_mat,file))
        H_hat = res_mat['H_hat']
        K1 = res_mat['K1'][0][0]
        if  res_mat['K2_left'][0][0]== 0:
            fm_mat = H_hat[:K1, :].T
        else:
            fm_mat = H_hat[:K1,:].T
            fm_mat = np.column_stack((fm_mat, res_mat['H2_trun'].T))

        for x in method.split('_'):
            if re.match(r'.*peak.*', x):
                peak = x[4:]
            if re.match(r'.*dim.*', x):
                dim = x[3:]
        df_metrics.loc[method,'Clustering methods']='Louvain clustering'
        df_metrics.loc[method,'Number of peaks']=res_mat['p'][0][0]
        df_metrics.loc[method,'Number of cells']=res_mat['n'][0][0]
        
        DRmethod = "RA3"
        df_metrics.loc[method, 'DR methods']= DRmethod
        df_metrics.loc[method, 'peak selection k'] = int(peak)/100
        df_metrics.loc[method, 'Dimension of latent space'] = dim
        df_metrics.loc[method, 'Clustering path'] = os.path.join(path_clusters, method + '_clusters.tsv')
        df_metrics.loc[method, 'Dataset'] = dataname

        adata = sc.AnnData(fm_mat)
        adata.obs['label'] = label.T.values
        # Louvain clustering
        np.random.seed(2020)
        sc.pp.neighbors(adata, n_neighbors=15,use_rep='X')
        #sc.tl.louvain(adata)
        getNClusters(adata,n_cluster=num_clusters)
        
        # adjusted rand index
        ari_louvain = adjusted_rand_score(adata.obs['label'], adata.obs['louvain'])
        # adjusted mutual information
        ami_louvain = adjusted_mutual_info_score(adata.obs['label'], adata.obs['louvain'],average_method='arithmetic')
        # homogeneity
        homo_louvain = homogeneity_score(adata.obs['label'], adata.obs['louvain'])
        # normalized mutual infromation
        nmi_louvain = normalized_mutual_info_score(adata.obs['label'], adata.obs['louvain'],average_method='arithmetic')
        # assign values
        df_metrics.loc[method,['ARI','AMI','Homogeneity','NMI']] = [ari_louvain,ami_louvain,homo_louvain,nmi_louvain]
        adata.obs[['label','louvain']].to_csv(os.path.join(path_clusters, method + '_clusters.tsv'),sep='\t')


df_metrics.to_csv(os.path.join(path_metrics,'clustering_scores_label_known.csv'), index=False)
