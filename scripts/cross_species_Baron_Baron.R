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
simple.evaluation <- function(true.tag, scRef.tag, df.ref.names, df.sc.names) {
    # uniform tags
    for (j in 1:dim(df.ref.names)[1]) {
        scRef.tag[scRef.tag == df.ref.names[j, 'ref.name']] <- df.ref.names[j, 'name']
    }
    scRef.tag[!(scRef.tag %in% df.ref.names$name)] <- 'Unassigned'
    for (j in 1:dim(df.sc.names)[1]) {
        true.tag[true.tag == df.sc.names[j, 'sc.name']] <- df.sc.names[j, 'name']
    }
    
    # true.labels <- setdiff(unique(true.tag), 'Unassigned')
    true.labels <- unique(true.tag)
    our.tag <- scRef.tag
    weighted_macro_f1 <- metrics$f1_score(true.tag, our.tag, average = 'weighted', labels = true.labels)
    macro_f1 <- metrics$f1_score(true.tag, our.tag, average = 'macro', labels = true.labels)
    accuracy <- metrics$accuracy_score(true.tag, our.tag)
    # rm unassigned in tags
    true.tag.rm <- true.tag[our.tag != 'Unassigned']
    our.tag.rm <- our.tag[our.tag != 'Unassigned']
    accuracy.rm.unassigned <- metrics$accuracy_score(true.tag.rm, our.tag.rm)
    
    f1 <- c()
    for (label in true.labels) {
        tmp.true.tag <- true.tag
        tmp.our.tag <- our.tag
        tmp.true.tag[tmp.true.tag != label] <- '0'
        tmp.our.tag[tmp.our.tag != label] <- '0'
        sub.f1 <- metrics$f1_score(tmp.true.tag, tmp.our.tag, average = 'binary', pos_label = label)
        f1 <- c(f1, sub.f1)
    }
    names(f1) <- true.labels
    
    our.labels <- setdiff(unique(our.tag), 'Unassigned')
    precision <- c()
    for (label in our.labels) {
        tmp.true.tag <- true.tag.rm
        tmp.our.tag <- our.tag.rm
        tmp.true.tag[tmp.true.tag != label] <- '0'
        tmp.our.tag[tmp.our.tag != label] <- '0'
        sub.precision <- metrics$precision_score(tmp.true.tag, tmp.our.tag, average = 'binary', pos_label = label)
        precision <- c(precision, sub.precision)
        
    }
    names(precision) <- our.labels
    mean.precision <- mean(precision)
    
    out <- list()
    out$weighted_macro_f1 <- weighted_macro_f1
    out$macro_f1 <- macro_f1
    out$accuracy <- accuracy
    out$f1 <- f1
    out$accuracy.rm.unassigned <- accuracy.rm.unassigned
    out$precision.rm.unassigned <- precision
    out$mean.precision.rm.unassigned <- mean.precision
    out$conf <- table(true.tag, our.tag)
    
    return(out)
    
}

source('/home/zy/my_git/scRef/main/scRef.v20.R')

############# regard sc-counts data as reference
path.input <- '/home/zy/scRef/sc_data/'
path.output <- '/home/zy/scRef/Benchmark/cross_species/'
dataset <- 'BaronM'
file.data.unlabeled <- paste0(path.input, dataset, '/cell_exp.txt')
file.label.unlabeled <- paste0(path.input, dataset, '/cell_meta.txt')
OUT <- prepare.data(file.data.unlabeled, file.label.unlabeled, del.label = c('Unclassified'))
saveRDS(OUT, file = paste0(path.output, dataset, '.Rdata'))
OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
ref.mtx <- OUT$data.filter
ref.labels <- OUT$label.filter$label.unlabeled.use.cols...
ref.dataset <- 'BaronM'

############### import unlabeled data
############### BaronH
path.input <- '/home/zy/scRef/sc_data/'
path.output <- '/home/zy/scRef/Benchmark/cross_species/'
dataset <- 'BaronH'
file.data.unlabeled <- paste0(path.input, dataset, '/cell_exp.txt')
file.label.unlabeled <- paste0(path.input, dataset, '/cell_meta.txt')
OUT <- prepare.data(file.data.unlabeled, file.label.unlabeled, del.label = c('miss'))
saveRDS(OUT, file = paste0(path.output, dataset, '.Rdata'))
OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
exp_sc_mat <- OUT$data.filter
label_sc <- OUT$label.filter

ref.names <- unique(ref.labels)
# list of cell names
all.cell <- unique(label_sc[,1])
uniform.names <- c("Neuron", "Endothelial Cell", "Astrocyte", "Microglia/PVM", 
                   "Oligo/OPC", "Oligo/OPC")
df.ref.names <- data.frame(ref.name = ref.names, name = uniform.names)
uniform.names <- c("Neuron", "Neuron", "Neuron", "Neuron", 
                   "Neuron", "Neuron", "Neuron", "Neuron", 
                   "Oligo/OPC", "Neuron", "Neuron", "Neuron", 
                   "Neuron", "Endothelial Cell", "Neuron", "Astrocyte", 
                   "Neuron", "Unassigned", "Neuron", "Microglia/PVM", 
                   "Neuron", "Unassigned", "Unassigned", "Neuron",  "Neuron")
df.sc.names <- data.frame(sc.name = all.cell, name = uniform.names)

# run methods
#############################################
### scRef
source('/home/zy/my_git/scRef/main/scRef.v20.R')
exp_sc_mat <- transform.HomoloGene(exp_sc_mat)
setwd('~/my_git/scRef')
result.scref <- SCREF(exp_sc_mat, ref.mtx, ref.labels,
                      type_ref = 'sc-counts', use.RUVseq = T, 
                      cluster.speed = F, cluster.cell = 3,
                      # GMM.num_cluster = 6, GMM.neg_cutoff = 2, GMM.floor_cutoff = 3, GMM.ceiling_cutoff = 20, 
                      min_cell = 3, CPU = 10)
pred.scRef <- result.scref$final.out$scRef.tag

true.tags <- label_sc$label.unlabeled.use.cols...
table(true.tags, pred.scRef)
