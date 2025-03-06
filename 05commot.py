import os
import gc
import ot
import pickle
import anndata
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse
from scipy.stats import spearmanr, pearsonr
from scipy.spatial import distance_matrix
import matplotlib.pyplot as plt

import commot as ct


out_dir = "/public/workspace/stu21230110/SPA_result/05commot"
annotation_dir = "/public/workspace/stu21230110/SPA_result/03cell2loc"
sample = ["C1","C2","C3","P1","P2","P3"]
data = "/public/workspace/stu21230110/SPA"
species = 'human'

name = sample[5]

data = data + '/' + name

annotation = annotation_dir + '/' + name + '/' + name + '.spatial.deconvolution.csv'

sample_name = name

input_sp = data

outdir = out_dir + '/' + name



adata = sc.read_visium(input_sp,library_id = name)

adata.var_names_make_unique()

sc.pp.normalize_total(adata, inplace=True)

sc.pp.log1p(adata)

adata_dis500 = adata.copy()

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

adata = adata[:, adata.var.highly_variable]

sc.tl.pca(adata, svd_solver='arpack')

#python3.11
#pip install scikit-learn==1.2.1
#pip install threadpoolctl==3.1.0
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

sc.tl.umap(adata)

sc.tl.leiden(adata, resolution=0.4)

df_cellchat = ct.pp.ligand_receptor_database(species=species, signaling_type='Secreted Signaling', database='CellChat')

df_cellchat_filtered = ct.pp.filter_lr_database(df_cellchat, adata_dis500, min_cell_pct=0.05)

ct.tl.spatial_communication(adata_dis500,
    database_name='cellchat', df_ligrec=df_cellchat_filtered, dis_thr=500, heteromeric=True, pathway_sum=True)

#import matplotlib as mpl
#default_backend = mpl.get_backend()
#mpl.use(default_backend)

for pathway in list(set(df_cellchat_filtered.iloc[:,2])):
    ct.tl.communication_direction(adata_dis500, database_name='cellchat', pathway_name=pathway, k=5)
    ct.pl.plot_cell_communication(adata_dis500, database_name='cellchat', pathway_name=pathway, plot_method='grid', background_legend=True,scale=0.0000014, ndsize=8, grid_density=0.4, summary='sender', background='image', clustering='leiden', cmap='Alphabet',normalize_v = True, normalize_v_quantile=0.995)
    plt.savefig(outdir + '/' + name + '.' + pathway + '.signal.arrow.grid.spatial.pdf',bbox_inches = 'tight')
    ct.pl.plot_cell_communication(adata_dis500, database_name='cellchat', pathway_name=pathway, plot_method='stream', background_legend=True,scale=0.00008, ndsize=8, grid_density=0.4, summary='sender', background='image', clustering='leiden', cmap='Alphabet',normalize_v = True, normalize_v_quantile=0.995)
    plt.savefig(outdir + '/' + name + '.' + pathway + '.signal.arrow.stream.spatial.png',bbox_inches = 'tight')

adata.write(outdir + '/' + name + '.spatial.h5ad')