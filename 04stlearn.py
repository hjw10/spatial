import os
import numpy as np
import stlearn as st
import pandas as pd
from matplotlib import pyplot as plt
import scanpy as sc
import random
from collections import Counter
import matplotlib
matplotlib.use('Agg')

out_dir = "/public/workspace/stu21230110/SPA_result/04stlearn"
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

sc.settings.verbosity = 1

sc.settings.set_figure_params(dpi=200,facecolor='white')

os.chdir(outdir)

adata = st.Read10X(input_sp)

adata.obs['sample'] = sample_name

adata.var_names_make_unique()
#st.add.image(adata=adata, imgpath=img,library_id=sample_name,visium=True)

st.pp.filter_genes(adata,min_cells=5)

st.pp.normalize_total(adata,target_sum=1e4) # NOTE: no log1p

adata.raw = adata

anno = pd.read_csv(annotation,sep=',',index_col = 0)

adata = adata[anno.index,:] ###filter Barcode

# 提取每一行的最大值及其对应的列名
#max_values = anno.max(axis=1)
#max_columns = anno.idxmax(axis=1)

# 创建新的数据表
#anno_max = pd.DataFrame({
#    'max_score': max_values,
#    'max_column': max_columns
#}, index=anno.index)

anno['max'] = anno.max(axis=1)

anno.loc[:,'prediction.score.max'] = anno.loc[:,'max']

anno.loc[:,'predicted.id'] = anno.idxmax(axis=1).to_frame(name='Label')['Label']

labels = anno.loc[:,'predicted.id'].values.astype(str)

spot_mixtures = anno.drop(['predicted.id','prediction.score.max'],axis=1)

spot_mixtures.columns = [col.replace('prediction.score.', '') for col in spot_mixtures.columns]

adata.obs['cell_type'] = labels

adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

adata.uns['cell_type'] = spot_mixtures

st.pl.cluster_plot(adata,use_label="cell_type",dpi=300, title='Cell type', fname= outdir + '/' + sample_name + '.annotation.pdf',size = 4,crop = False)


###Running the Ligand-Receptor Analysis

lrs = st.tl.cci.load_lrs(['connectomeDB2020_lit'], species=species)

###spot+immediate neighbours

st.tl.cci.run(adata, lrs,
                  min_spots = 5, #Filter out any LR pairs with no scores for less than min_spots
                  distance=None, # None defaults to spot+immediate neighbours; distance=0 for within-spot mode
                  n_pairs=2000, # Number of random pairs to generate; low as example, recommend ~10,000
                  n_cpus=8, # Number of CPUs for parallel. If None, detects & use all available.
                  )


#adata.write(outdir + '/' + sample_name + '.spatial.h5ad')

lr_info = adata.uns['lr_summary']

#lr_info.to_csv(outdir + '/' + sample_name + '.lr_info.spatial.csv')

####P-value adjustment

st.tl.cci.adj_pvals(adata, correct_axis='spot',pval_adj_cutoff=0.01, adj_method='fdr_bh')

####Visualise

st.pl.lr_summary(adata, n_top=500)

plt.savefig(outdir + '/' + sample_name + '.LR_Rank_500.pdf',dpi = 300)

st.pl.lr_summary(adata, n_top=50, figsize=(10,7))

plt.savefig(outdir + '/' + sample_name + '.LR_Rank_50.pdf',dpi = 300)

####A key aspect of the LR analysis is to control for LR expression level and frequency when calling significant hotspots.

####Hence, our diagnostic plots should show next to no correlation between the hotspots of the LR pairs and the expression level and frequency of expression.

st.pl.lr_diagnostics(adata, figsize=(10,4.5))

plt.savefig(outdir + '/' + sample_name + '.Diagnostic.pdf',dpi = 300)

####Left plot: Relationship between LR expression level (non-zero spots average median expression of genes in the LR pair) and the ranking of the LR.

####Right plot: Relationship between LR expression frequency (average proportion of zero spots for each gene in the LR pair) and the ranking of the LR.

st.pl.lr_n_spots(adata, n_top=50, figsize=(11, 8),max_text=100)

plt.savefig(outdir + '/' + sample_name + '.LR.sig.top50.pdf',dpi = 300)

st.pl.lr_n_spots(adata, n_top=500, figsize=(11, 5),max_text=100)

plt.savefig(outdir + '/' + sample_name + '.LR.sig.top500.pdf',dpi = 300)

####LR Statistics Visualisations(top 10)

adata.uns['lr_summary']['interaction_score'] = sum(np.asarray(adata.obsm['lr_scores']))

lr_summary = pd.DataFrame(adata.uns['lr_summary'])
lr_summary.to_csv(outdir + '/' + sample_name + '.lr_summary.spatial.csv')

best_lr = pd.DataFrame(adata.uns['lr_summary']).sort_values(['n_spots_sig','interaction_score'],ascending=[False,False]).index.values[:20]

lr_interaction_score = pd.DataFrame(adata.obsm['lr_scores'],index = adata.obs.index,columns = adata.uns['lr_summary'].index)

lr_interaction_score.to_csv(outdir + '/' + sample_name + '.LR.interaction.score.csv')

lr_interaction_pval = pd.DataFrame(adata.obsm['p_adjs'],index = adata.obs.index,columns = adata.uns['lr_summary'].index)

lr_interaction_pval.to_csv(outdir + '/' + sample_name + '.LR.interaction.padj.csv')

for lr in best_lr:
	stats = ['lr_scores','lr_sig_scores', 'p_vals', 'p_adjs', '-log10(p_adjs)']
	fig, axes = plt.subplots(ncols=len(stats), figsize=(20,10))
	for i, stat in enumerate(stats):
		st.pl.lr_result_plot(adata, use_result=stat, use_lr=lr, show_color_bar=False, ax=axes[i])
		axes[i].set_title(f'{lr} {stat}')
		plt.savefig(outdir + '/' + sample_name + '.LR.%s.pdf'%(lr),dpi = 300)
		try:
			st.pl.lr_plot(adata, lr, inner_size_prop=1, outer_mode='binary', pt_scale=10,use_label=None,show_image=True,sig_spots=True)
			plt.savefig(outdir + '/' + sample_name + '.LR.%s.interaction.pdf'%(lr),dpi = 300) ###The receptor is in green, the ligand is in red. The inner-point is the receptor, the outter point is the ligand.
			st.pl.lr_plot(adata,lr,inner_size_prop=0.04, middle_size_prop=.07, outer_size_prop=.4,outer_mode='continuous',
pt_scale=60,use_label=None, show_image=True,sig_spots=True)
			plt.savefig(outdir + '/' + sample_name + '.LR.%s.coexpression.pdf'%(lr),dpi = 300)
		except:
			print("%s is not a lr pair ~~~"%(lr))
####Predicting significant CCIs


st.pl.gene_plot(adata, gene_symbols="COL1A1", contour=True,cell_alpha=0.5,figsize=(10, 6))
plt.savefig(outdir + '/' + sample_name + ".COL1A1.pdf",dpi = 300)


# Running the counting of co-occurence of cell types and LR expression hotspots #
st.tl.cci.run_cci(adata, 'cell_type', # Spot cell information either in data.obs or data.uns
                  min_spots=3, # Minimum number of spots for LR to be tested.
                  spot_mixtures=True, # If True will use the label transfer scores,
                                      # so spots can have multiple cell types if score>cell_prop_cutoff
                  cell_prop_cutoff=0.1, # Spot considered to have cell type if score>0.1
                  sig_spots=True, # Only consider neighbourhoods of spots which had significant LR scores.
                  n_perms=1000 # Permutations of cell information to get background, recommend ~1000
                 )

####Diagnostic plot: check interaction and cell type frequency correlation

adata.write(outdir + '/' + sample_name + '.spatial.interaction.h5ad')

adata = sc.read(outdir + '/' + sample_name + '.spatial.interaction.h5ad')

st.pl.cci_check(adata, 'cell_type',figsize=(16, 14))

plt.savefig(outdir + '/' + sample_name + '.CCI.LR.interaction.pdf',dpi = 300)

####CCI network

# Visualising the no. of interactions between cell types across all LR pairs #

pos_1 = st.pl.ccinet_plot(adata, 'cell_type', return_pos=True)

plt.savefig(outdir + '/' + sample_name + '.cell-cell.interaction.pdf',dpi = 300)

lrs = adata.uns['lr_summary'].index.values[0:3]

for lr in lrs[0:3]:
    
	st.pl.ccinet_plot(adata, 'cell_type', lr, min_counts=2 , figsize=(10,7.5), pos=pos_1)

	plt.savefig(outdir + '/' + sample_name + '.cell-cell.interaction.%s.pdf'%(lr),dpi = 300)


####CCI chord-plot

st.pl.lr_chord_plot(adata, 'cell_type',title = 'interaction strength')

plt.savefig(outdir + '/' + sample_name + '.cell-cell.interaction.circle.pdf',dpi = 300)

for lr in lrs:

	st.pl.lr_chord_plot(adata, 'cell_type', lr,title = 'interaction strength')
	
	plt.savefig(outdir + '/' + sample_name + '.cell-cell.interaction.%s.circle.pdf'%(lr),dpi = 300)

####Heatmap Visualisations

st.pl.lr_cci_map(adata, 'cell_type', lrs=lrs, min_total=100, figsize=(20,8))

plt.savefig(outdir + '/' + sample_name + '.cell-cell.interaction.heatmap.pdf',dpi = 300)

st.pl.cci_map(adata, 'cell_type',figsize = (20,8))

plt.savefig(outdir + '/' + sample_name + '.cell.to.cell.interaction.heatmap.pdf',dpi = 300)

lrs = adata.uns['lr_summary'].index.values[0:5]

for lr in lrs[0:5]:

	st.pl.cci_map(adata, 'cell_type', lr)

	plt.savefig(outdir + '/' + sample_name + '.cell.to.cell.interaction.%s.heatmap.pdf'%(lr),dpi = 300)




#adata.write(outdir + '/' + sample_name + 'spatial.h5ad')
