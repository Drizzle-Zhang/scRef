# import python package: sklearn.metrics
library(reticulate)
library(ggplot2)
library(ggalluvial)
library(networkD3)
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
simple.evaluation <- function(true.tag, scRef.tag, df.cell.names) {
    # uniform tags
    for (j in 1:dim(df.cell.names)[1]) {
        scRef.tag[scRef.tag == df.cell.names[j, 'ref.name']] <- 
            df.cell.names[j, 'sc.name']
    }
    
    true.tag[true.tag %in% unknow.cell] <- 'Unassigned'
    true.labels <- unique(true.tag)
    our.tag <- scRef.tag
    weighted_macro_f1 <- metrics$f1_score(true.tag, our.tag, average = 'weighted', labels = true.labels)
    macro_f1 <- metrics$f1_score(true.tag, our.tag, average = 'macro', labels = true.labels)
    accuracy <- metrics$accuracy_score(true.tag, our.tag)
    
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
    
    out <- list()
    out$weighted_macro_f1 <- weighted_macro_f1
    out$macro_f1 <- macro_f1
    out$accuracy <- accuracy
    out$f1 <- f1
    out$conf <- table(true.tag, our.tag)
    
    return(out)
    
}

source('/home/zy/my_git/scRef/main/scRef.v12.R')

############# regard sc-counts data as reference
library(stringr)
file.mtx <- '/home/disk/scRef/MouseAtlas_SingleCell_Han2018/MCA_by_tissue/Brain/Count_all_batch.txt'
df.mtx <- read.delim(file.mtx, stringsAsFactors = F, row.names = 1)
file.cellid <- '/home/disk/scRef/MouseAtlas_SingleCell_Han2018/MCA_by_tissue/Brain/CellAssignments_all_batch.txt'
df.cellid <- read.delim(file.cellid, stringsAsFactors = F, row.names = 1)
df.labels <- df.cellid[colnames(df.mtx),]
ref.labels <- df.labels$CellType
ref.mtx <- df.mtx
ref.dataset <- 'MCA'

############### import unlabeled data
############### Habib
path.input <- '/home/zy/scRef/summary/'
path.output <- '/home/zy/scRef/Benchmark/mouse_brain/'
dataset <- 'Habib'
file.data.unlabeled <- paste0(path.input, dataset, '_exp_sc_mat.txt')
file.label.unlabeled <- paste0(path.input, dataset, '_exp_sc_mat_cluster_original.txt')
# OUT <- prepare.data(file.data.unlabeled, file.label.unlabeled, del.label = c('miss'))
# saveRDS(OUT, file = paste0(path.output, dataset, '.Rdata'))
OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
exp_Habib <- OUT$data.filter
label_Habib <- OUT$label.filter
exp_sc_mat <- exp_Habib
label_sc <- label_Habib

ref.names <- unique(ref.labels)
# list of cell names
all.cell <- unique(label_sc[,1])
sc.name <- c("Oligodend", "microglia", "Astrocyte", "Neurons", "Unassigned",
             "Unassigned", "Unassigned", "OPC", "Unassigned", "Unassigned", "Ependymocytes")
unknow.cell <- setdiff(all.cell, sc.name)
df.cell.names <- data.frame(ref.name = ref.names, sc.name = sc.name, idx = 1:length(sc.name))

path.output <- '/home/zy/scRef/figure/single_two_round/'
true.tags <- label_sc$label.unlabeled.use.cols...

# run methods
#############################################
### scRef two-round
source('/home/zy/my_git/scRef/main/scRef.v19.R')
setwd('~/my_git/scRef')
result.scref <- SCREF(exp_sc_mat, ref.mtx, ref.labels,
                      type_ref = 'sc-counts', use.RUVseq = T, 
                      cluster.speed = T, cluster.cell = 5,
                      min_cell = 10, CPU = 8)
pred.scRef <- result.scref$final.out$scRef.tag
saveRDS(pred.scRef, file = paste0(path.output, 'two_round_scRef.Rdata'))

rda.scRef <- paste0(path.output, 'two_round_scRef.Rdata')
pred.scRef <- readRDS(rda.scRef)
# res.scRef <- simple.evaluation(true.tags, pred.scRef, df.ref.names, df.sc.names)

# rename
df.ref.names <- data.frame(ref.name = ref.names, 
                           name = c('Oligodendrocyte', "Microglia", "Astrocyte", "Neuron", "Macrophage", 
                                    "Granulocyte", "Oligodendrocyte precursor cell", "Schwann cell", 
                                    "Astroglial cell", "Hypothalamic ependymal cell"))
df.sc.names <- data.frame(sc.name = all.cell, 
                          name = c("Neuron", "ParsTuber", "MuralCells", "Hypothalamic ependymal cell", 
                                   'Oligodendrocyte', "Tanycyte", "Endothelial Cell", "Fibroblast",
                                   "Astrocyte", "Microglia", "Oligodendrocyte precursor cell"))
for (j in 1:dim(df.ref.names)[1]) {
    pred.scRef[pred.scRef == df.ref.names[j, 'ref.name']] <- df.ref.names[j, 'name']
}
for (j in 1:dim(df.sc.names)[1]) {
    true.tags[true.tags == df.sc.names[j, 'sc.name']] <- df.sc.names[j, 'name']
}

mytable.two <- table(true.tags, pred.scRef)
mydata.two <- data.frame(stringsAsFactors = F)
for (label1 in rownames(mytable.two)) {
    for (label2 in colnames(mytable.two)) {
        mydata.two <- rbind(mydata.two, data.frame(origin = label1, annotation = label2, 
                                           count = mytable.two[label1, label2]))
    }
}
table.tags <- table(label_sc$label.unlabeled.use.cols...)
table.annos <- table(pred.scRef)
mydata.two$origin <- factor(mydata.two$origin, levels = names(table.tags)[order(table.tags, decreasing = T)])
mydata.two$annotation <- factor(mydata.two$annotation, levels = names(table.annos)[order(table.annos, decreasing = T)])
names <- sort(unique(c(levels(mydata.two$origin), levels(mydata.two$annotation))))
names <- data.frame(name=c(names))
sankeyNetwork(Links = mydata.two, Nodes = names, 
              Source = 'origin', Target = 'annotation', Value = 'count')

#########
# plot.1 <- ggplot(data = mydata,
#        aes(axis1 = origin, axis2 = annotation, y = count)) +
#     scale_x_discrete(limits = c("Origin labels", "scRef annotations"), expand = c(.01, .05)) +
#     geom_alluvium(aes(fill = origin)) +
#     geom_stratum(width = 0.25) + 
#     geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) + 
#     labs(fill = 'Origin labels', y = '') + 
#     theme(
#         legend.text = element_text(size = 12),
#         legend.title = element_text(size = 15),
#         panel.grid = element_blank(),
#         panel.background = element_rect(fill = 'transparent'),
#         axis.text.x = element_text(size = 15, vjust = 15),
#         panel.border = element_blank(),
#         axis.ticks = element_blank(),
#         axis.line = element_blank(),
#         axis.text.y = element_blank()
#     )
# pathout <- '/home/zy/scRef/figure'
# ggsave(filename = 'two-round.png', 
#        path = pathout, plot = plot.1,
#        units = 'cm', height = 40, width = 30)
####################

### scRef single-round
source('/home/zy/my_git/scRef/main/scRef.v12.R')
setwd('~/my_git/scRef')
result.scref.single <- SCREF(exp_sc_mat, ref.mtx, ref.labels,
                      type_ref = 'sc-counts', use.RUVseq = T, single_round = T, 
                      method1 = 'kendall',
                      cluster.speed = T, cluster.cell = 10,
                      min_cell = 10, CPU = 8)
pred.scRef.single <- result.scref.single$final.out$scRef.tag
mytable <- table(label_sc$label.unlabeled.use.cols..., pred.scRef.single)
mydata <- data.frame(stringsAsFactors = F)
for (label1 in rownames(mytable)) {
    for (label2 in colnames(mytable)) {
        mydata <- rbind(mydata, data.frame(origin = label1, annotation = label2, 
                                           count = mytable[label1, label2]))
    }
}
table.tags <- table(label_sc$label.unlabeled.use.cols...)
table.annos <- table(pred.scRef)
mydata$origin <- factor(mydata$origin, levels = names(table.tags)[order(table.tags, decreasing = T)])
mydata$annotation <- factor(mydata$annotation, levels = names(table.annos)[order(table.annos, decreasing = T)])



#########################
plot.2 <- ggplot(data = mydata,
                 aes(axis1 = origin, axis2 = annotation, y = count)) +
    scale_x_discrete(limits = c("Origin labels", "scRef annotations"), expand = c(.01, .05)) +
    geom_alluvium(aes(fill = origin)) +
    geom_stratum(width = 0.25) + 
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) + 
    labs(fill = 'Origin labels', y = '') + 
    theme(
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = 'transparent'),
        axis.text.x = element_text(size = 15, vjust = 15),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.y = element_blank()
    )
pathout <- '/home/zy/scRef/figure'
ggsave(filename = 'single-round.png', 
       path = pathout, plot = plot.2,
       units = 'cm', height = 40, width = 30)
##########################


