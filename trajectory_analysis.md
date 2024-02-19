## Trajectory analysis for monocyte and macrophage populations (related to figure3)

```python
import palantir
import scanpy as sc
import numpy as np
import os

# Plotting 
import matplotlib
import matplotlib.pyplot as plt
# Inline plotting
%matplotlib inline
# Reset random seed
np.random.seed(5)
import pandas as pd
import pickle

monomac_metadata = pd.read_csv("../palantir/monomac_metadata.csv", index_col = 0)
monomac_count = pd.read_csv("../palantir/monomac_count.csv", index_col = 0) # this matrix is normalized gene expression matrix with hvgs

data_matrix = monomac_count
metadata = monomac_metadata
pcadims = [150]
tsne_perplexity = [700]
num_wp = 2000
start_cell = "CTGCGAGTCTGGAAGG-S1" # start cell was chosen randomly in the monocyte cluster.
filename = "monomac"
savedir = "../palantir/monomac_new/"

i = pcadims
pca_projections, _ = palantir.utils.run_pca(data_matrix, n_components=i, use_hvg=False)
pca_projections.shape
pca_projections.index = data_matrix.index
# diffusion maps
dm_res = palantir.utils.run_diffusion_maps(pca_projections , n_components=10)
ms_data = palantir.utils.determine_multiscale_space(dm_res)
imp_df = palantir.utils.run_magic_imputation(data_matrix, dm_res)
# Running palantir
with open(savedir + filename + '_ms_data_' + str(i) + '.pickle', 'wb') as f:
    pickle.dump(ms_data, f, pickle.HIGHEST_PROTOCOL)
    
pr_res = palantir.core.run_palantir(ms_data, start_cell, num_waypoints=num_wp)
with open(savedir + filename + '_pr_res_' + str(i) + '.pickle', 'wb') as f:
    pickle.dump(pr_res, f, pickle.HIGHEST_PROTOCOL)
        
pr_res.branch_probs.to_csv(savedir + filename + '_branch_probs_' + str(i) + '.csv')
pr_res.pseudotime.to_csv(savedir + filename + '_pseudotime_' + str(i) + '.csv')
    
j = tsne_perplexity
tsne = palantir.utils.run_tsne(ms_data, perplexity=j)
tsne.to_csv(savedir + filename + '_tsne_' + str(i) + '_' + str(j) + '.csv')
        
palantir.plot.highlight_cells_on_tsne(tsne, start_cell)
plt.savefig(savedir + filename + '_pca_' + str(i) + '_' + str(j) + '_tsne_start_cell.png')

palantir.plot.plot_palantir_results(pr_res, tsne)
plt.savefig(savedir + filename + '_pca_' + str(i) + '_' + str(j) + '_pr_res.png')
        
palantir.plot.plot_cell_clusters(tsne, metadata['condition'])
plt.savefig(savedir + filename + '_pca_' + str(i) + '_' + str(j) + '_tsne_condition.png')
        
palantir.plot.plot_cell_clusters(tsne, metadata['new.celltype'])
plt.savefig(savedir + filename + '_pca_' + str(i) + '_' + str(j) + '_tsne_celltype.png')
```

## Trajectory analysis for adipose lineage (related to figure6)

```python
Adipo_metadata = pd.read_csv("../palantir/Adipo_metadata.csv", index_col = 0)
Adipo_count = pd.read_csv("../palantir/Adipo_count.csv", index_col = 0) # this matrix is normalized gene expression matrix with hvgs

data_matrix = Adipo_count
metadata = Adipo_metadata
pcadims = [50,100,150]
tsne_perplexity = [300,500,700,800,900,1000]
num_wp = 5000
start_cell = "GTAGAAACACCAATTG-S6" # expressing high Dpp4 and Cd34
filename = "APC_Adipocyte"
savedir = "../palantir/APC_Adipocyte_final/"

i = pcadims
pca_projections, _ = palantir.utils.run_pca(data_matrix, n_components=i, use_hvg=False)
pca_projections.shape
pca_projections.index = data_matrix.index
# diffusion maps
dm_res = palantir.utils.run_diffusion_maps(pca_projections , n_components=15)
ms_data = palantir.utils.determine_multiscale_space(dm_res)
imp_df = palantir.utils.run_magic_imputation(data_matrix, dm_res)
# Running palantir
with open(savedir + filename + '_ms_data_' + str(i) + '.pickle', 'wb') as f:
    pickle.dump(ms_data, f, pickle.HIGHEST_PROTOCOL)

pr_res = palantir.core.run_palantir(ms_data, start_cell, num_waypoints=num_wp)
with open(savedir + filename + '_pr_res_' + str(i) + '.pickle', 'wb') as f:
    pickle.dump(pr_res, f, pickle.HIGHEST_PROTOCOL)
        
pr_res.branch_probs.to_csv(savedir + filename + '_branch_probs_' + str(i) + '.csv')
pr_res.pseudotime.to_csv(savedir + filename + '_pseudotime_' + str(i) + '.csv')

j = tsne_perplexity
tsne = palantir.utils.run_tsne(ms_data, perplexity=j)
tsne.to_csv(savedir + filename + '_tsne_' + str(i) + '_' + str(j) + '.csv')
        
palantir.plot.highlight_cells_on_tsne(tsne, start_cell)
plt.savefig(savedir + filename + '_pca_' + str(i) + '_' + str(j) + '_tsne_start_cell.png')

palantir.plot.plot_palantir_results(pr_res, tsne)
plt.savefig(savedir + filename + '_pca_' + str(i) + '_' + str(j) + '_pr_res.png')
        
palantir.plot.plot_cell_clusters(tsne, metadata['condition'])
plt.savefig(savedir + filename + '_pca_' + str(i) + '_' + str(j) + '_tsne_condition.png')
```

## Trajectory analysis for monocyte and macrophage populations (related to figureS4)

```R
library(monocle3)

# generate cds
monomac$celltype <- as.character(monomac$detailed.celltype)

cell_metadata = monomac@meta.data
gene_metadata = data.frame(gene_short_name = rownames(monomac), row.names = rownames(monomac))

cds <- new_cell_data_set(monomac@assays$originalexp@data, 
                         cell_metadata = cell_metadata, 
                         gene_metadata = gene_metadata)
# pre-process the data
cds <- preprocess_cds(cds, "PCA", num_dim = 20, norm_method = "none", use_genes = VariableFeatures(monomac))
plot_pc_variance_explained(cds)

# Reduce dimensionality and visualize the cells
cds <- reduce_dimension(cds, preprocess_method = 'PCA') # UMAP
cds = cluster_cells(cds, resolution = 0.0001)
plot_cells(cds, color_cells_by = "cluster")
cds = learn_graph(cds, learn_graph_control = list(prune_graph = TRUE))

plot_cells(cds, color_cells_by = "cluster", 
           label_groups_by_cluster = F, 
           label_leaves = F, 
           label_branch_points = F, group_label_size = 5)

cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "celltype",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE, trajectory_graph_color='black',cell_stroke = 0
           )+ scale_color_manual(values = c('#ffb703','#606c38','#2d8fd5','#ff4d6d'))+
  theme(axis.title = element_blank(),legend.position = 'none', axis.text = element_blank(),
        axis.line.x = element_blank(),axis.ticks = element_blank(),axis.line.y = element_blank())




