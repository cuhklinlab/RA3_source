import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from sklearn.metrics import confusion_matrix
import os

from scale.plot import plot_confusion_matrix, plot_embedding, plot_heatmap
from scale.utils import read_labels, reassign_cluster_with_ref, pairwise_pearson
from scale.specifity import *

%load_ext autoreload
%autoreload 2
%matplotlib inline

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

from scipy.optimize import linear_sum_assignment as linear_assignment
from sklearn.preprocessing import LabelEncoder

def reassign_cluster_with_ref(Y_pred, Y):
    """
    Reassign cluster to reference labels
    Inputs:
        Y_pred: predict y classes
        Y: true y classes
    Return:
        f1_score: clustering f1 score
        y_pred: reassignment index predict y classes
        indices: classes assignment
    """
    def reassign_cluster(y_pred, index):
        y_ = np.zeros_like(y_pred)
        for i, j in zip(index[0], index[1]):
            y_[np.where(y_pred==i)] = j
        return y_
    assert Y_pred.size == Y.size
    D = max(Y_pred.max(), Y.max())+1
    w = np.zeros((D,D), dtype=np.int64)
    for i in range(Y_pred.size):
        w[Y_pred[i], Y[i]] += 1
    ind = linear_assignment(w.max() - w)

    return reassign_cluster(Y_pred, ind)



mainDRdim = '10'
ourDRdim = '5'
peakR = 3
multi_donor = True
for data_name in ['MPP_LMPP_CLP','donor_BM0828','ALL_blood','forebrain_half','tissue4_forebrain']:
    for DRmethod in ["scABC","SCALE","Scasat","cisTopic","Cusanovich2018","SnapATAC","RA3"]:
        if DRmethod=="RA3":
            DRdim = ourDRdim
            if multi_donor==True and data_name in ['ALL_blood','MPP_LMPP_CLP']:
                file_name = '%s_%s_peak%03d_dim%s_multidonor'%(DRmethod,data_name,peakR*100,DRdim)
            else:
                file_name = '%s_%s_peak%03d_dim%s'%(DRmethod,data_name,peakR*100,DRdim)
            label_file_name = '/home/zhangwenyu/dataForComparison_clustering/clusters_our/%s_clusters.tsv'%(file_name)
        else:
            if DRmethod=="scABC":
                DRdim = '999'
            else:
                DRdim = mainDRdim
            file_name = '%s_%s_peak%03d_dim%s'%(DRmethod,data_name,peakR*100,DRdim)
            label_file_name = '/home/zhangwenyu/dataForComparison_clustering/clusters/%s_clusters.tsv'%(file_name)
        if os.path.exists(label_file_name):
            cluster_assign = pd.read_csv(label_file_name, sep='\t', index_col=0)
        else:
            print('===[DOES NOT EXIST]=== %s'%label_file_name)
            continue
        ref = cluster_assign['label'].values
        try:
            pred = cluster_assign['louvain'].values
        except:
            pred = cluster_assign['scABC'].values
            
        if data_name == 'tissue4_forebrain':
            ref[np.where(ref == "Astrocytes")] = 'AC'
            ref[np.where(ref == "Ex. neurons CPN")] = 'EX CPN'
            ref[np.where(ref == "Ex. neurons CThPN")] = 'EX CThPN'
            ref[np.where(ref == "Ex. neurons SCPN")] = 'EX SCPN'
            ref[np.where(ref == "Inhibitory neurons")] = 'IN'
            ref[np.where(ref == "Microglia")] = 'MG'
            ref[np.where(ref == "Oligodendrocytes")] = 'OC'
        

        encode = LabelEncoder()
        ref_label = encode.fit_transform(ref.squeeze())
        ref_class = encode.classes_

        encode = LabelEncoder()
        pred_label = encode.fit_transform(pred.squeeze())
        pred_class = encode.classes_

        reassign_pred = reassign_cluster_with_ref(pred_label, ref_label)

        cm = confusion_matrix(reassign_pred, ref_label, labels=np.unique(ref_label))

        file_name = '%s-%s'%(data_name,DRmethod)
        plot_confusion_matrix(cm, ref_class, np.unique(ref_label)+1, title=file_name, normalize=False, 
                              figsize=(5,5),show_cbar=False, save='./results/plot/1112/confusion_matrix/'+file_name+'.pdf')
