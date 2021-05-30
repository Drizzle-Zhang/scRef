library(reticulate)
use_python('/home/zy/tools/anaconda3/bin/python3', required = T)
py_config()

# reference
path.output <- '/home/zy/scRef/Benchmark/mouse_brain/'
OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
exp_Tasic <- OUT$mat_exp
label_Tasic <- OUT$label
ref.labels <-label_Tasic[, 1]
ref.mtx <- exp_Tasic
ref.dataset <- 'Tasic'

# sample test dataset
dataset <- 'Campbell'
OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
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
path.output <- '/home/zy/scRef/Benchmark/mouse_brain/'
ref.dataset <- 'Tasic'
OUT <- readRDS(paste0(path.output, ref.dataset, '.Rdata'))
ref.mtx <- OUT$mat_exp
ref.labels <-OUT$label[, 1]

list.demo <- readRDS('/home/zy/my_git/scMAGIC_scripts/data/Campbell_2k.Rdata')
exp_sc_mat <- list.demo$mat_exp
label_sc <-list.demo$label

source('/home/zy/my_git/scRef/main/scRef.v21.R')
setwd('~/my_git/scRef')
result.scref <- SCREF(exp_sc_mat, ref.mtx, ref.labels, CPU = 4)
print(result.scref$run.time)
pred.scMAGIC <- result.scref$final.out$scRef.tag
table(label_sc, pred.scMAGIC)


