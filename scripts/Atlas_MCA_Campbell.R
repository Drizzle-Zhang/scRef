# import python package: sklearn.metrics
library(reticulate)
use_python('/home/zy/tools/anaconda3/bin/python3', required = T)
# py_config()
py_module_available('sklearn')
metrics <- import('sklearn.metrics')

# function of data preparation
prepare.data <- function(file.data.unlabeled, file.label.unlabeled, 
                         del.label = c('miss')) {
    library(stringr)
    data.unlabeled <- read.delim(file.data.unlabeled, row.names=1)
    data.unlabeled <- floor(data.unlabeled)
    names(data.unlabeled) <- str_replace_all(names(data.unlabeled), '_', '.')
    names(data.unlabeled) <- str_replace_all(names(data.unlabeled), '-', '.')
    # read label file
    file.label.unlabeled <- file.label.unlabeled
    label.unlabeled <- read.delim(file.label.unlabeled, row.names=1)
    row.names(label.unlabeled) <- str_replace_all(row.names(label.unlabeled), '_', '.')
    row.names(label.unlabeled) <- str_replace_all(row.names(label.unlabeled), '-', '.')
    col.name1 <- names(data.unlabeled)[1]
    if (substring(col.name1, 1, 1) == 'X') {
        row.names(label.unlabeled) <- paste0('X', row.names(label.unlabeled))
    }
    # filter data
    use.cols <- row.names(label.unlabeled)[!label.unlabeled[,1] %in% del.label]
    data.filter <- data.unlabeled[,use.cols]
    label.filter <- data.frame(label.unlabeled[use.cols,], row.names = use.cols)
    
    OUT <- list()
    OUT$data.filter <- data.filter
    OUT$label.filter <- label.filter
    return(OUT)
    
}

path.input <- '/home/zy/scRef/sc_data/'
path.output <- '/home/zy/scRef/atlas_anno'
dataset <- 'Campbell'
file.data.unlabeled <- paste0(path.input, dataset, '_exp_sc_mat.txt')
file.label.unlabeled <- paste0(path.input, dataset, '_exp_sc_mat_cluster_original.txt')
# OUT <- prepare.data(file.data.unlabeled, file.label.unlabeled, del.label = c('miss'))
# saveRDS(OUT, file = paste0(path.output, dataset, '.Rdata'))
OUT <- readRDS(paste0('/home/zy/scRef/Benchmark/mouse_brain/', dataset, '.Rdata'))
exp_Habib <- OUT$mat_exp
label_Habib <- OUT$label
exp_sc_mat <- exp_Habib
label_sc <- label_Habib

source('/home/zy/my_git/scRef/main/scRef.v20.R')
setwd('~/my_git/scRef')
df.atlas <- .imoprt_outgroup('MCA', normalization = F)
# df.atlas <- df.atlas[, colnames(df.atlas) != 'Pan gabaergic']

result.scref <- SCREF(exp_sc_mat, df.atlas,
                      type_ref = 'sum-counts', use.RUVseq = F, 
                      method1 = 'kendall',
                      cluster.speed = T, cluster.cell = 5,
                      min_cell = 3, CPU = 8)
pred.scRef <- result.scref$final.out$scRef.tag

true.tags <- label_sc$label.unlabeled.use.cols...
table(true.tags, pred.scRef)
# df.tags <- result.scref$combine.out
# df.view <- merge(label_sc, df.tags, by = 'row.names')
# View(df.view)

library(ggplot2)
path.res <- '/home/zy/scRef/figure/atlas_anno'

# heatmap
method <- 'scMAGIC'
mytable <- table(true.tags, pred.scRef)
mydata <- data.frame(stringsAsFactors = F)
table.true <- table(true.tags)
for (label1 in rownames(mytable)) {
    row.sum <- table.true[label1]
    for (label2 in c(colnames(df.atlas), "Unassigned")) {
        if (label2 %in% colnames(mytable)) {
            mydata <- rbind(mydata, data.frame(origin = label1, annotation = label2, 
                                               count = mytable[label1, label2], 
                                               prop = mytable[label1, label2]/row.sum))
        } else {
            mydata <- rbind(mydata, data.frame(origin = label1, annotation = label2, 
                                               count = 0, prop = 0))
        }
    }
}
mydata$origin <- factor(mydata$origin, levels = c(rownames(mytable)))
ref.cells <- c("Unassigned", 
               "Astrocyte", "Astroglial cell", 
               "Vascular endothelial cell", 
               "Ovarian vascular surface endothelium cell",
               "Hypothalamic ependymal cell", 
               "Dopaminergic neurons", "Ganglion cell", "Granule neuron", "Hippocampus neuron",
               "Interstitial macrophage", "Microglia",
               "Oligodendrocyte precursor cell")
mydata$annotation <- factor(mydata$annotation, 
                            levels = c(ref.cells, setdiff(colnames(df.atlas), ref.cells)))

plot.heatmap <- 
    ggplot(data = mydata, aes(x = origin, y = annotation)) + 
    geom_tile(aes(fill = prop)) + 
    scale_fill_gradient2(low = "#000000", high = "#FFFF00", mid = "#32CD32", midpoint = 0.5) + 
    labs(fill = 'Proportion') + 
    theme_bw() +
    theme(
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    ) 
    # geom_text(aes(label = round(prop, 2)), family = "Arial", size = 2.5)
ggsave(filename = paste0('heatmap_allMCA_', dataset, '_', method, '.png'), 
       path = path.res, plot = plot.heatmap,
       units = 'cm', height = 60, width = 20)




