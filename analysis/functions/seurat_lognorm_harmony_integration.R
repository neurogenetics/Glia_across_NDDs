library(Seurat)
library(scCustomize)
library(harmony)

lognorm_harmonyintegrate <- function(seurat_obj, n_varfeats = 2000, vars_to_regress, harmony_batches,
                                     n_dims = NULL, cluster_res){
  seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = n_varfeats)
  seurat_obj <- ScaleData(seurat_obj, vars.to.regress = vars_to_regress)
  seurat_obj <- RunPCA(seurat_obj)
  
  seurat_obj <- RunHarmony(object = seurat_obj, reduction = "pca", group.by.vars = harmony_batches,
                           reduction.save = 'harmony', plot_convergence = T, lambda = NULL)
  print(ElbowPlot(seurat_obj, reduction = "harmony", ndims = 50))
  
  if (is.null(n_dims)) {
    n_dims <- as.numeric(readline(prompt = "Enter the number of dimensions to use for neighbors/UMAP: "))
  }
  
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:n_dims)
  seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:n_dims, reduction.name = "umap_harmony")
  seurat_obj <- FindClusters(seurat_obj, resolution = cluster_res)
  
}