#! /PROJ/development/zhaoyunfei/miniconda/miniconda/envs/scrna/bin/python
### 20220325
### zhaoyunfei
### https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html
### https://github.com/BayraktarLab/cell2location


import os
import sys
import gc
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import random
import cell2location
import scvi
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
import seaborn as sns
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel
from cell2location.plt import plot_spatial
from cell2location import run_colocation
from scipy.sparse import csr_matrix
import warnings
warnings.filterwarnings('ignore')
import igraph
import scipy

results_folder = './03cell2loc/'
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'
outdir = "/public/workspace/stu21230110/SPA_result/03cell2loc/"

if os.path.exists(outdir) is not True:os.mkdir(outdir)

os.system('cd %s'%(outdir))

sc.settings.verbosity = 3

sc.settings.set_figure_params(dpi=200, facecolor='white')

if sc_format == '10x':

	sc_adata = sc.read_10x_mtx(sc_input,var_names='gene_symbols',cache=True)

elif sc_format == 'h5':

	sc_adata=sc.read_10x_h5(sc_input,genome=None, gex_only=True)

elif sc_format == 'h5ad':

	sc_adata=sc.read("scRNA.h5ad")

else :

	sc_adata=sc.read_loom(sc_input,sparse = True)
	
anno=pd.read_csv("celltype.csv",index_col=1)

sc_adata.var.index = sc_adata.var['_index']

sc_adata = sc_adata[anno.index,:]

sc_adata.obs['celltype']=anno['V2']

sc_adata.var_names_make_unique()

#sc_adata.var['SYMBOL'] = sc_adata.var.index

#sc_adata.var_names = sc_adata.var['gene_ids']
###cell filter
#Barcode = pd.read_csv(celllabels,index_col = 0)

#sc_adata = sc_adata[Barcode.index,:]

sc_adata.obs['cluster'] = sc_adata.obs['celltype'].astype('category')##cluster因子型

sc_adata = sc_adata[~sc_adata.obs['cluster'].isna(), :]##去掉没有注释的细胞

try :

	sc_adata.obs['orig.ident']

except KeyError :

	print("KeyError:'sample' is not exist in sc_adata.obs,plesae set this information")

	exit(0)

####filter gene
sc_adata.var['MT_gene'] = [gene.startswith('MT-') for gene in sc_adata.var.index]

set(sc_adata.var['MT_gene'])

sc_adata.obsm['MT'] = sc_adata[:, sc_adata.var['MT_gene'].values].X.toarray()

sc_adata = sc_adata[:, ~sc_adata.var['MT_gene'].values]

sc.pp.filter_genes(sc_adata, min_cells=5)

# calculate the mean of each gene across non-zero cells
sc_adata.var['n_cells'] = (sc_adata.X.toarray() > 0).sum(0)

sc_adata.var['nonz_mean'] = sc_adata.X.toarray().sum(0) / sc_adata.var['n_cells']

nonz_mean_cutoff = np.log10(2) # cut off for expression in non-zero cells

# select genes based on mean expression in non-zero cells
sc_adata = sc_adata[:,np.array(np.log10(sc_adata.var['nonz_mean']) > nonz_mean_cutoff)]

sc_adata = sc_adata.copy()
# prepare anndata for the regression model

TF_CPP_MIN_LOG_LEVEL=0

if len(set(sc_adata.obs['orig.ident'])) > 1:

	cell2location.models.RegressionModel.setup_anndata(adata=sc_adata,
                        # cell type, covariate used for constructing signatures
                        labels_key='cluster',
                        batch_key='orig.ident'
                       )

else:

	cell2location.models.RegressionModel.setup_anndata(adata=sc_adata,labels_key='cluster')


mod = RegressionModel(sc_adata)

# Use all data for training (validation not implemented yet, train_size=1)
mod.train(max_epochs=500,train_size=1, lr=0.002)
###1:11:00##

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
sc_adata = mod.export_posterior(
    sc_adata, sample_kwargs={'num_samples': 1000}
)

sc_adata.__dict__['_raw'].__dict__['_var'] = sc_adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'genes'})
sc_adata.write(outdir + '/' + 'train.scRNA.h5ad')



#####Spatial data

sp_name = ["C1","C2","C3","P1","P2","P3"]

name = sp_name[0]

#sc_adata=sc.read(outdir + "/train.scRNA.h5ad")

outdir = "/public/workspace/stu21230110/SPA_result/03cell2loc/" + name

#adata_vis = sc.read(outdir + "/" + name + ".spatial.scRNA.h5ad")

adata_vis = sc.read_visium(path="/public/workspace/stu21230110/SPA/" + name + "/", library_id = name)

adata_vis.obs['sample'] = name

adata_vis.var_names_make_unique()

#adata_vis.var['SYMBOL'] = adata_vis.var.index

#adata_vis.var.index = adata_vis.var['gene_ids']

adata_vis.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_vis.var.index]

# remove MT genes for spatial mapping (keeping their counts in the object)
adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()

adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]
####spatial mapping

# export estimated expression in each cluster
inf_aver = sc_adata.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                for i in sc_adata.uns['mod']['factor_names']]].copy()



inf_aver.columns = sc_adata.uns['mod']['factor_names']

intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)

adata_vis = adata_vis[:, intersect].copy()

inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis)

# create and train the model
mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=30,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection (using default here):
    detection_alpha=200
)

mod.train(max_epochs=3000,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1)

##1:23:00
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs}
)

adata_vis.var.drop(columns='gene_ids', inplace=True)

####Visualising cell abundance in spatial coordinates
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']

adata_vis.obs[adata_vis.uns['mod']['factor_names']].to_csv(outdir + '/' + name + '.spatial.deconvolution.csv')

adata_vis.obsm['means_cell_abundance_w_sf'] = scipy.sparse.csr_matrix(adata_vis.obsm['means_cell_abundance_w_sf'])
adata_vis.obsm['q05_cell_abundance_w_sf'] = scipy.sparse.csr_matrix(adata_vis.obsm['q05_cell_abundance_w_sf'])
adata_vis.obsm['q95_cell_abundance_w_sf'] = scipy.sparse.csr_matrix(adata_vis.obsm['q95_cell_abundance_w_sf'])
adata_vis.obsm['stds_cell_abundance_w_sf'] = scipy.sparse.csr_matrix(adata_vis.obsm['stds_cell_abundance_w_sf'])
adata_vis.obsm['means_cell_abundance_w_sf'].columns = ['scvi_1','scvi_2','scvi_3','scvi_4','scvi_5','scvi_6','scvi_7','scvi_8','scvi_9','scvi_10','scvi_11']

adata_vis.write(outdir + '/' + name + '.spatial.scRNA.h5ad')

#name = sp_name[0]
#outdir = "/public/workspace/stu21230110/SPA_result/03cell2loc/" + name
#adata_vis = sc.read(outdir + "/" + name + ".spatial.scRNA.h5ad")

# plot in spatial coordinates
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
	
	for celltype in list(set(sc_adata.obs['cluster'])):
	
		sc.pl.spatial(adata_vis, cmap='bwr',##bwr(蓝红) hsv(变化明显) PuRd（红）
            # show first 8 cell types
            color=celltype,
            size=1.3,
            img_key='hires',
            # limit color scale at 99.2% quantile of cell abundance
            vmin=0.05, vmax='p99.2')
		plt.savefig(outdir + '/' + name + '_' + celltype + '_spatial.pdf',bbox_inches = 'tight')


	
clust_labels = random.sample(list(set(sc_adata.obs['cluster'])),3)

clust_labels=['Stromal_cells','Macrophage','Tissue_stem_cells']
clust_col = ['' + str(i) for i in clust_labels]
with mpl.rc_context({'figure.figsize': (17, 15)}):
	fig = plot_spatial(
	adata=adata_vis,
# labels to show on a plot
	color=clust_col, labels=clust_labels,
	show_img=True,
# 'fast' (white background) or 'dark_background'
	style='fast',
# limit color scale at 99.2% quantile of cell abundance
	max_color_quantile=0.992,
# size of locations (adjust depending on figure size)
	circle_diameter=6,
	colorbar_position='right')
plt.savefig(outdir + '/' + name + '_' + '_'.join(clust_col) + '_spatial.pdf',bbox_inches = 'tight')


####Downstream analysis
#####NMF 

# compute KNN using the cell2location output stored in adata.obsm
sc.pp.neighbors(adata_vis, use_rep='q05_cell_abundance_w_sf',
                n_neighbors = 15)

# Cluster spots into regions using scanpy
sc.tl.leiden(adata_vis, resolution=0.5)

# add region as categorical variable
adata_vis.obs["region_cluster"] = adata_vis.obs["leiden"].astype("category")

# compute UMAP using KNN graph based on the cell2location output
sc.tl.umap(adata_vis, min_dist = 0.3, spread = 1)

# plot in spatial coordinates
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [7, 5]}):
	sc.pl.spatial(adata_vis, color=['region_cluster'],
                  size=1.3, img_key='hires', alpha=0.8)
	
plt.savefig(outdir + '/' + name + '_cluster_spatial.pdf',bbox_inches = 'tight')

####Identifying cellular compartments / tissue zones using matrix factorisation (NMF)
res_dict, adata_vis = run_colocation(
    adata_vis,
    model_name='CoLocatedGroupsSklearnNMF',
    train_args={
      'n_fact': np.arange(11, 13), # IMPORTANT: use a wider range of the number of factors (5-30)
      'sample_name_col': 'sample', # columns in adata_vis.obs that identifies sample
      'n_restarts': 3 # number of training restarts
    },
    export_args={'path': f'{outdir}/NMF/'}
)

####cell_type_fractions_heatmap/: a dot plot of the estimated NMF weights of cell types (rows) across NMF components (columns)

####cell_type_fractions_mean/: the data used for dot plot

####factor_markers/: tables listing top 10 cell types most speficic to each NMF factor

####models/: saved NMF models

####predictive_accuracy/: 2D histogram plot showing how well NMF explains cell2location output

####spatial/: NMF weights across locatinos in spatial coordinates

####location_factors_mean/: the data used for the plot in spatial coordiantes

####stability_plots/: stability of NMF weights between training restarts

# Here we plot the NMF weights (Same as saved to `cell_type_fractions_heatmap`)
res_dict['n_fact12']['mod'].plot_cell_type_loadings()

plt.savefig(outdir + '/' + name + '_NMF_celltype.pdf',bbox_inches = 'tight')

adata_vis.write(outdir + '/' + name + '.Plot.spatial.scRNA.h5ad')
