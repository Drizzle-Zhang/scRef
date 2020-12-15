# import python package: sklearn.metrics
library(reticulate)
use_python('/home/zy/tools/anaconda3/bin/python3', required = T)
# py_config()
py_module_available('sklearn')
metrics <- import('sklearn.metrics')

source('/home/zy/my_git/scRef/main/scRef.v19.R')

############# regard sc-counts data as reference
path.input <- '/home/zy/scRef/summary/'
path.output <- '/home/zy/scRef/Benchmark/mouse_brain/'
dataset <- 'Tasic'
file.data.unlabeled <- paste0(path.input, dataset, '_exp_sc_mat.txt')
file.label.unlabeled <- paste0(path.input, dataset, '_exp_sc_mat_cluster_merged.txt')
OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
exp_Tasic <- OUT$data.filter
label_Tasic <- OUT$label.filter
ref.labels <-label_Tasic[,1]
ref.mtx <- exp_Tasic
ref.dataset <- 'Tasic'

############### import unlabeled data
############### Habib
path.input <- '/home/zy/scRef/summary/'
path.output <- '/home/zy/scRef/Benchmark/mouse_brain/'
dataset <- 'Habib'
file.data.unlabeled <- paste0(path.input, dataset, '_exp_sc_mat.txt')
file.label.unlabeled <- paste0(path.input, dataset, '_exp_sc_mat_cluster_original.txt')
OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
exp_Habib <- OUT$data.filter
label_Habib <- OUT$label.filter
exp_sc_mat <- exp_Habib
label_sc <- label_Habib

# run methods
#############################################
### scRef
source('/home/zy/my_git/scRef/main/scRef.v20.R')
setwd('~/my_git/scRef')
res.scRef <- SCREF(exp_sc_mat, ref.mtx, ref.labels,
                      type_ref = 'sc-counts', use.RUVseq = T, 
                      cluster.speed = T, cluster.cell = 5,
                      min_cell = 10, CPU = 8)
path.res <- '/home/zy/scRef/figure/mouse_brain'
file.res <- paste0(path.res, 'list_result_', ref.dataset, '_', dataset, '_scRef.Rdata')
saveRDS(res.scRef, file = file.res)
res.scRef <- readRDS(file.res)

res.scRef <- annotate.UnassignedCell(res.scRef, exp_sc_mat, atlas = 'MCA', CPU = 10)

### original plot
library(Seurat)
# data preparing
seurat.unlabeled <- CreateSeuratObject(counts = exp_sc_mat)
seurat.unlabeled <- NormalizeData(seurat.unlabeled, normalization.method = "LogNormalize", 
                                  scale.factor = 10000)
seurat.unlabeled <- FindVariableFeatures(seurat.unlabeled, selection.method = "vst", nfeatures = 2000)
# VariableFeatures(seurat.unlabeled)
# VariableFeaturePlot(seurat.unlabeled)
# seurat.unlabeled[["percent.mt"]] <- PercentageFeatureSet(seurat.unlabeled, pattern = "^mt-")
# seurat.unlabeled <- ScaleData(seurat.unlabeled, vars.to.regress = "percent.mt")
# seurat.unlabeled <- ScaleData(seurat.unlabeled, features = all.genes)
seurat.unlabeled <- ScaleData(seurat.unlabeled)

# add label
use.cells <- dimnames(seurat.unlabeled@assays$RNA@counts)[[2]]
names(label_sc) <- 'cell_type'
# label_sc <- cbind(label_sc, seurat.unlabeled@reductions$umap@cell.embeddings)
seurat.unlabeled@meta.data$original.label <- label_sc[use.cells, 'cell_type']
pred.scRef <- res.scRef$final.out
# pred.scRef <- cbind(pred.scRef, seurat.unlabeled@reductions$umap@cell.embeddings)
seurat.unlabeled@meta.data$scRef.tag <- pred.scRef$scRef.tag
pred.unassign <- res.scRef$pred.new[use.cells,]
# pred.unassign <- cbind(pred.unassign, seurat.unlabeled@reductions$umap@cell.embeddings)
seurat.unlabeled@meta.data$new.tag <- pred.unassign$scRef.tag

# PCA
seurat.unlabeled <- RunPCA(seurat.unlabeled, npcs = 100, verbose = F)

# cluster
# seurat.unlabeled <- FindNeighbors(seurat.unlabeled, reduction = "pca", dims = 1:75, nn.eps = 0.5)
# seurat.unlabeled <- FindClusters(seurat.unlabeled, resolution = 3, n.start = 20)

# UMAP
seurat.unlabeled <- RunUMAP(seurat.unlabeled, dims = 1:100, n.neighbors = 50)

library(ggplot2)
# figure1: ture label
plot.umap <- 
    DimPlot(seurat.unlabeled, reduction = "umap", label = T, repel = T, group.by = 'original.label') + 
    xlim(-16, 14) +
    theme_bw() + 
    theme(axis.text = element_text(size = 9),
          panel.grid = element_blank(),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 11))
ggsave(filename = paste0('cluster_', ref.dataset, '_', dataset, '.png'), 
       path = path.res, plot = plot.umap,
       units = 'cm', height = 18, width = 24)

# figure2: cluster label
# DimPlot(seurat.unlabeled, reduction = "umap", label = T, group.by = 'seurat_clusters')
# figure3: scRef plus label
plot.umap.scRef <- 
    DimPlot(seurat.unlabeled, reduction = "umap", label = T, repel = T, group.by = 'scRef.tag') + 
    scale_color_manual(values = c('#F8766D', '#D39200', '#00BA38', '#00B9E3', '#619CFF', '#DB72FB', 'gray'),
                       breaks = c('Astrocyte', 'Endothelial Cell', 'Microglia', 'Neuron', 
                                  'Oligodendrocyte', 'Oligodendrocyte Precursor Cell', 'Unassigned')) + 
    theme_bw() + 
    theme(axis.text = element_text(size = 9),
          panel.grid = element_blank(),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 11))
ggsave(filename = paste0('cluster_scRef_', ref.dataset, '_', dataset, '.png'), 
       path = path.res, plot = plot.umap.scRef,
       units = 'cm', height = 18, width = 26)

# figure3: scRef and annotate unassigned
plot.umap.scRef.unassign <- 
    DimPlot(seurat.unlabeled, reduction = "umap", label = T, repel = T, group.by = 'new.tag') + 
    scale_color_manual(values = c('#F8766D', '#D39200', '#00BA38', '#00B9E3', '#619CFF', '#DB72FB', 
                                  '#93AA00', '#00C19F', '#A52A2A', 'gray'),
                       breaks = c('Astrocyte', 'Endothelial Cell', 'Microglia', 'Neuron', 
                                  'Oligodendrocyte', 'Oligodendrocyte Precursor Cell', 
                                  'Hypothalamic ependymal cell', 'Oligodendrocyte precursor cell',
                                  'Astroglial cell', 'Unassigned')) + 
    xlim(-18, 15) +
    theme_bw() + 
    theme(axis.text = element_text(size = 9),
          panel.grid = element_blank(),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 11))
ggsave(filename = paste0('cluster_scRef_unassign_', ref.dataset, '_', dataset, '.png'), 
       path = path.res, plot = plot.umap.scRef.unassign,
       units = 'cm', height = 18, width = 26)


