library(reticulate)
use_python('/home/zy/tools/anaconda3/bin/python3', required = T)
py_config()

# reference
path.output <- '/mdshare/zy/scRef/Benchmark/mouse_brain/'
list.Ref <- readRDS(paste0(path.output, dataset, '.Rdata'))
exp_Tasic <- list.Ref$mat_exp
label_Tasic <- list.Ref$label
ref.labels <-label_Tasic[, 1]
ref.mtx <- exp_Tasic
ref.dataset <- 'Tasic'

# sample test dataset
dataset <- 'Campbell'
OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
# exp_sc_mat <- OUT$mat_exp
# label_sc <- OUT$label
exp_Habib <- OUT$mat_exp
label_Habib <- OUT$label
set.seed(1234)
cells.sample <- sample(rownames(label_Habib),2000)
exp_sc_mat <- exp_Habib[,cells.sample]
label_sc <- label_Habib[cells.sample,]
list.demo <- list()
list.demo$mat_exp <- exp_sc_mat
list.demo$label <- label_sc
saveRDS(list.demo, file = '/home/zy/my_git/scMAGIC_scripts/data/Campbell_2k.Rdata')

# load data
path.output <- '/mdshare/zy/scRef/Benchmark/mouse_brain/'
ref.dataset <- 'Tasic'
list.Ref <- readRDS(paste0(path.output, ref.dataset, '.Rdata'))
ref.mtx <- list.Ref$mat_exp
ref.labels <-list.Ref$label[, 1]

list.demo <- readRDS('/local/zy/my_git/scMAGIC_scripts/data/Campbell_2k.Rdata')
exp_sc_mat <- list.demo$mat_exp
label_sc <-list.demo$label

source('/local/zy/my_git/scRef/main/scRef.v21.R')
setwd('~/my_git/scRef')
output.scMAGIC <- SCREF(exp_sc_mat, ref.mtx, ref.labels, CPU = 4)
print(output.scMAGIC$run.time)
pred.scMAGIC <- output.scMAGIC$final.out$scRef.tag
table(label_sc, pred.scMAGIC)

library(Seurat)
seurat.unlabeled <- CreateSeuratObject(counts = exp_sc_mat)
seurat.unlabeled <- NormalizeData(seurat.unlabeled)
seurat.unlabeled <- FindVariableFeatures(seurat.unlabeled, nfeatures = 2000)
seurat.unlabeled <- ScaleData(seurat.unlabeled)
seurat.unlabeled <- RunPCA(seurat.unlabeled)
seurat.unlabeled <- RunUMAP(seurat.unlabeled, dims = 1:50)
seurat.unlabeled@meta.data$original.label <- label_sc
seurat.unlabeled@meta.data$pred.tag <- pred.scMAGIC


library(ggplot2)
path.res <- "/mdshare/node9/zy/scRef/demo"
# origin label
plot.umap <-
    DimPlot(seurat.unlabeled, reduction = "umap",
            label = T, repel = T, group.by = 'original.label') +
    labs(title = 'True labels') + theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = 'TrueLabels.png',
       path = path.res, plot = plot.umap,
       units = 'cm', height = 15, width = 20)
# pred label
plot.umap.pred <-
    DimPlot(seurat.unlabeled, reduction = "umap",
            label = T, repel = T, group.by = 'pred.tag') +
    labs(title = 'Prediction labels') + theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = 'PredLabels.png',
       path = path.res, plot = plot.umap.pred,
       units = 'cm', height = 15, width = 22)


library(ggplot2)
plot.umap <-
    DimPlot(seurat.unlabeled, reduction = "umap",
            label = T, repel = T, label.size = 2.5,
            group.by = 'original.label') +
    scale_color_manual(values = c('#24B700', '#E18A00', '#BE9C00', '#00BE70', '#24B700', '#8CAB00',
                                  '#00C1AB', '#00BBDA', '#00ACFC', '#8B93FF', '#D575FE'),
                       breaks = c('Astrocytes', 'Endothelial cells', 'Ependymocytes', 'Mural cells',
                                  'Neurons', 'Oligodendrocytes', 'OPC', 'Pars tuberalis',
                                  'PVMs & Microglia', 'Tanycytes', 'VLMCs')) +
    theme_bw() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_text(size = 8, face = 'bold'),
          legend.text = element_text(size = 6),
          legend.position = 'bottom') +
    guides(color = guide_legend(ncol = 3,  byrow = TRUE, reverse = F,
                                override.aes = list(size=3),
                                keywidth = 0.1, keyheight = 0.1, default.unit = 'cm'))
ggsave(filename = paste0('cluster_', ref.dataset, '_', dataset, '.png'),
       path = path.res, plot = plot.umap,
       units = 'cm', height = 10.5, width = 8.5)


plot.umap.pred <-
    DimPlot(seurat.unlabeled, reduction = "umap",
            label = T, repel = T, label.size = 2.5,
            group.by = 'pred.tag') +
    scale_color_manual(values = c('#24B700', '#E18A00', '#00ACFC', '#24B700', '#8CAB00', '#00C1AB', 'gray'),
                       breaks = c('Astrocyte', 'Endothelial Cell', 'Microglia', 'Neuron',
                                  'Oligodendrocyte', 'Oligodendrocyte Precursor Cell', 'Unassigned')) +
    theme_bw() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_text(size = 8, face = 'bold'),
          legend.text = element_text(size = 6),
          legend.position = 'bottom') +
    guides(color = guide_legend(ncol = 3,  byrow = TRUE, reverse = F,
                                override.aes = list(size=3),
                                keywidth = 0.1, keyheight = 0.1, default.unit = 'cm'))
ggsave(filename = paste0('cluster_scRef_', ref.dataset, '_', dataset, '.png'),
       path = path.res, plot = plot.umap.scRef,
       units = 'cm', height = 10, width = 8.5)
