
After removing putative doublets and low quality cells (as described in manuscript), remaining cells were subclustered and annotated as below (related to figure3C)

```R
library(scater)
library(scran)
library(harmony)
library(Seurat)
library(SingleCellExperiment)

# subclustering function
subclustering_harmony <- function(seurat, condition=condition, idents=NULL, cell=NULL, PCA, resolution, invert = F, fdr.threshold = 0.05, hvg.number = NULL, doublet.score = FALSE){
  subset <- subset(seurat, idents=idents, cells = cell, invert = invert)
  subset <- DietSeurat(subset)
  SCEset <- as.SingleCellExperiment(subset)
  
  clusters <- quickCluster(SCEset)
  SCEset <- computeSumFactors(SCEset, clusters = clusters)
  SCEset <- logNormCounts(SCEset, pseudo_count = 1)
  
  if(doublet.score == TRUE){
    library(scDblFinder)
    SCEset <- scDblFinder(SCEset, samples="batch")
  }
  
  set.seed(123)
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
  name<- RunUMAP(name, dims = 1:PCA, seed.use = 42, n.neighbors = 30)
  name<- RunHarmony(name, condition, plot_convergence = TRUE)
  name<- FindNeighbors(name, reduction = "harmony", dims = 1:PCA)
  name<- FindClusters(name, resolution = resolution)
  name<- RunUMAP(name, reduction = "harmony", dims = 1:PCA)
  return(name)
}

# subcluster remaining cells without doublets and low quality cells
v3_seurat <- subclustering_harmony(seurat_v2, condition = 'batch', cell = c(colnames(seurat_v2)[seurat_v2$seurat_clusters %in% c(12,27,29,30)],colnames(c26)[c26$seurat_clusters %in% c(0,2,3)]), invert=T,PCA=15, resolution=0.8)

# subcluster for immune cell annotation
v3_immune <- subclustering_harmony(v3_seurat, condition = 'batch', idents = c(7,8,9,11), PCA=15, resolution = 1.0)
v3_immune$celltype <- 'NA'
v3_immune$celltype[v3_immune$seurat_clusters %in% c(20)] <- 'DC'
v3_immune$celltype[v3_immune$seurat_clusters %in% c(19)] <- 'B'
v3_immune$celltype[v3_immune$seurat_clusters %in% c(16)] <- 'T'
v3_immune$celltype[v3_immune$seurat_clusters %in% c(15,2,13,14)] <- 'Monocyte'
v3_immune$celltype[v3_immune$seurat_clusters %in% c(17)] <- 'Mac.Prg4'
v3_immune$celltype[v3_immune$seurat_clusters %in% c(0,6,9,12,21)] <- 'Mac.Lyve1'
v3_immune$celltype[v3_immune$seurat_clusters %in% c(1,3,4,5,7,8,10,11,18)] <- 'Mac.Trem2'

v3_seurat$celltype <- 'Adipocyte'
v3_seurat$celltype[v3_seurat$seurat_clusters %in% c(7,8,9,11)] <- v3_immune$celltype
v3_seurat$celltype[v3_seurat$seurat_clusters %in% c(1,12)] <- 'MC'
v3_seurat$celltype[v3_seurat$seurat_clusters %in% c(13)] <- 'LEC'
v3_seurat$celltype[v3_seurat$seurat_clusters %in% c(10,15)] <- 'VEC'
v3_seurat$celltype[v3_seurat$seurat_clusters %in% c(2,6)] <- 'APC'

# subcluster for APC and SMC annotation
APC <- subclustering_harmony(v3_seurat, condition = 'batch',idents = c(2,6), PCA=10, resolution = 0.5)
APC$celltype <- 'APC'
APC$celltype[APC$seurat_clusters %in% c(8)] <- 'SMC'

v3_seurat$celltype[v3_seurat$celltype == 'APC'] <- APC$celltype
```


