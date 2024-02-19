The co-expressing gene modules were obtained using hdWGCNA R pakcage (related to figure4 and S5)

```R
library(Seurat)
library(hdWGCNA)
library(patchwork)

# subcluster WT_HFD condition in monocyte and macrophage populations for analysis
subclustering_func <- function(seurat, idents=NULL, cell=NULL, PCA, resolution, invert = F, fdr.threshold = 0.05, hvg.number = NULL, doublet.score = FALSE){
  set.seed(123)
  subset <- subset(seurat, idents = idents, cells = cell, invert = invert)
  subset <- DietSeurat(subset)
  SCEset <- as.SingleCellExperiment(subset)
  
  clusters <- quickCluster(SCEset)
  SCEset <- computeSumFactors(SCEset, clusters = clusters)
  print(summary(sizeFactors(SCEset)))
  SCEset <- logNormCounts(SCEset, pseudo_count = 1)
  
  if(doublet.score == TRUE){
    library(scDblFinder)
    SCEset <- scDblFinder(SCEset, samples = 'batch')
  }
  
  dec <- modelGeneVar(SCEset)
  hvg.norm <- getTopHVGs(dec, fdr.threshold = fdr.threshold, n = hvg.number)
  print(paste0("The number of HVG is ",length(hvg.norm)))
  
  
  name <- as.Seurat(SCEset)
  name@assays$originalexp@var.features <- hvg.norm
  name <- ScaleData(name, features = rownames(name))
  name@reductions$PCA <- NULL
  name@reductions$UMAP <- NULL
  
  PCA=PCA
  name<- RunPCA(name, npcs = PCA)
  name<- FindNeighbors(name, dims= 1:PCA)
  name<- FindClusters(name, resolution = resolution)
  name<- RunUMAP(name, dims = 1:PCA, seed.use = 42, n.neighbors = 30)
  return(name)
}

monomac_WTHFD <- subclustering_func(monomac, cell = colnames(monomac)[monomac$condition == 'WT_HFD'], PCA = 10, resolution = 0.8)

set.seed(1234)

monomac_WTHFD <- SetupForWGCNA(
  monomac_WTHFD,
  gene_select = "fraction", # the gene selection approach
  wgcna_name = "fraction_gene", # the name of the hdWGCNA experiment
  fraction = 0.05
)

# construct metacells  in each group
monomac_WTHFD <- MetacellsByGroups(
  seurat_obj = monomac_WTHFD,
  group.by = c("celltype", "condition"), # specify the columns in seurat_obj@meta.data to group by
  k = 25, # nearest-neighbors parameter
  max_shared = 15, # maximum number of shared cells between two metacells
  ident.group = 'celltype' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
monomac_WTHFD <- NormalizeMetacells(monomac_WTHFD)

# optional : process the metacell seurat object
metacell_monomac <- GetMetacellObject(monomac_WTHFD)

# select soft-power thershold
# Test different soft powers:
monomac_WTHFD <- TestSoftPowers(
  monomac_WTHFD,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(monomac_WTHFD)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)

power_table <- GetPowerTable(monomac_WTHFD)

# construct co-expression network:
monomac_WTHFD <- ConstructNetwork(
  monomac_WTHFD, soft_power=3,
  setDatExpr=FALSE, overwrite_tom=TRUE
)

PlotDendrogram(monomac_WTHFD, main='monomac hdWGCNA Dendrogram') # FigureS5A was generated from this code

# Module Eigengenes and connectivity 
# compute harmonized module eigengenes
monomac_WTHFD <- ScaleData(monomac_WTHFD, features=VariableFeatures(monomac_WTHFD))

# compute all MEs in the full single-cell dataset
monomac_WTHFD <- ModuleEigengenes(
 monomac_WTHFD
)

# module eigengenes:
MEs <- GetMEs(monomac_WTHFD, harmonized=FALSE)

# Compute module connectivity
# compute eigengene-based connectivity (kME):
monomac_WTHFD <- ModuleConnectivity(
  monomac_WTHFD
)

# Getting the module assignment table
# get the module assignment table:
modules <- GetModules(monomac_WTHFD)

# make a featureplot of MEs for each module
plot_list <- ModuleFeaturePlot(
  monomac_WTHFD,
  features='MEs', # plot the MEs
  order=TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=2) # Figure 4B and S5C were generated from this code

# plot module correlagram
ModuleCorrelogram(hdWGCNA_monomac_subset) # Figure S5B was generated from this code
```
