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

# evaluation
simple.evaluation <- function(true.tag, )

source('/home/zy/my_git/scRef/main/scRef.v10.R')

############# regard counts data as reference
path.input <- '/home/zy/scRef/summary/'
path.output <- '/home/zy/scRef/Benchmark/mouse_brain/'
dataset <- 'Tasic'
file.data.unlabeled <- paste0(path.input, dataset, '_exp_sc_mat.txt')
file.label.unlabeled <- paste0(path.input, dataset, '_exp_sc_mat_cluster_merged.txt')
# OUT <- prepare.data(file.data.unlabeled, file.label.unlabeled, del.label = c('Unclassified'))
# saveRDS(OUT, file = paste0(path.output, dataset, '.Rdata'))
OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
data.filter <- OUT$data.filter
label.filter <- OUT$label.filter
label.in <- data.frame(cell_id = row.names(label.filter), tag = label.filter$label.unlabeled.use.cols...)
exp.Tasic.sum <- .generate_ref(data.filter, label.in, M='SUM')
exp_ref_mat <- exp.Tasic.sum

############### import unlabeled data
############### Habib
dataset <- 'Habib'
file.data.unlabeled <- paste0(path.input, dataset, '_exp_sc_mat.txt')
file.label.unlabeled <- paste0(path.input, dataset, '_exp_sc_mat_cluster_original.txt')
# OUT <- prepare.data(file.data.unlabeled, file.label.unlabeled, del.label = c('miss'))
# saveRDS(OUT, file = paste0(path.output, dataset, '.Rdata'))
OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
exp_sc_mat <- OUT$data.filter
label.origin <- OUT$label.filter

ref.names <- colnames(exp_ref_mat)
# list of cell names
all.cell <- unique(label.origin[,1])
sc.name <- c("Astrocyte", "EndothelialCells",
             "microglia", "Neurons", "Oligodend", "OPC")
unknow.cell <- setdiff(all.cell, sc.name)
df.cell.names <- data.frame(ref.name = ref.names, sc.name = sc.name, idx = 1:length(sc.name))


#################
### scRef
source('/home/zy/my_git/scRef/main/scRef.v10.R')
setwd('~/my_git/scRef')
result.scref <- SCREF(exp_sc_mat, exp_ref_mat, type_ref = 'count', 
                      cluster.speed = T, cluster.cell = 10,
                      min_cell = 10, CPU = 8)

library(ggplot2)
library(mclust)
sub.astrocyte <- df.tags[df.tags$scRef.tag == 'Astrocyte', ]
sub.astrocyte <- df.tags1[df.tags1$scRef.tag == 'Astrocyte', ]
ggplot(sub.astrocyte, aes(x = log10Pval)) + geom_density()
model.astrocyte <- densityMclust(sub.astrocyte$log10Pval, G=5)
summary(model.astrocyte, parameters = T)
sub.astrocyte$cluster <- model$classification

sub.neuron <- df.tags[df.tags$scRef.tag.12 == 'Neuron', ]
sub.neuron <- df.tags1[df.tags1$scRef.tag == 'Neuron', ]
ggplot(sub.neuron, aes(x = log10Pval)) + geom_density()
model.neuron <- densityMclust(sub.neuron$log10Pval, G=8)
summary(model.neuron, parameters = T)

sub.Endo <- df.tags1[df.tags1$scRef.tag == 'Endothelial Cell', ]
ggplot(sub.Endo, aes(x = log10Pval)) + geom_density()
model.Endo <- densityMclust(sub.Endo$log10Pval)
summary(model.Endo, parameters = T)

sub.Microglia <- df.tags[df.tags$scRef.tag == 'Microglia', ]
ggplot(sub.Microglia, aes(x = log10Pval)) + geom_density()
model.Microglia <- densityMclust(sub.Microglia$log10Pval)
summary(model.Microglia, parameters = T)

sub.Oligo <- df.tags[df.tags$scRef.tag == 'Oligodendrocyte', ]
ggplot(sub.Oligo, aes(x = log10Pval)) + geom_density()
model.Oligo <- densityMclust(sub.Oligo$log10Pval)
summary(model.Oligo, parameters = T)

sub.OPC <- df.tags1[df.tags1$scRef.tag == 'Oligodendrocyte Precursor Cell', ]
ggplot(sub.OPC, aes(x = log10Pval)) + geom_density()
model.OPC <- densityMclust(sub.OPC$log10Pval, G = 5)
summary(model.OPC, parameters = T)

meta.tag <- merge(result.scref$final.out, label.origin, by = 'row.names')
row.names(meta.tag) <- meta.tag$Row.names
meta.tag$Row.names <- NULL
names(meta.tag) <- c('scRef.tag', 'log10Pval', 'ori.tag')

### evaluation
true.tag = meta.tag$ori.tag
scRef.tag = meta.tag$scRef.tag

# uniform tags
for (j in 1:dim(df.cell.names)[1]) {
    scRef.tag[scRef.tag == df.cell.names[j, 'ref.name']] <- 
        df.cell.names[j, 'sc.name']
}
meta.tag$scRef.tag <- scRef.tag

# default cutoff
true.tag[true.tag %in% unknow.cell] <- 'Unassigned'
our.tag <- meta.tag$scRef.tag
metrics$f1_score(true.tag, our.tag, average = 'weighted')
metrics$f1_score(true.tag, our.tag, average = 'macro')
metrics$accuracy_score(true.tag, our.tag)

