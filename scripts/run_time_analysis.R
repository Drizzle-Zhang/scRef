############### Tasic2018
path.output <- '/home/zy/scRef/Benchmark/mouse_brain/'
dataset <- 'Tasic2018'
OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
ref_mat <- OUT$mat_exp
label_ref <- OUT$label

############### Campbell
path.output <- '/home/zy/scRef/Benchmark/mouse_brain/'
dataset <- 'Campbell'
OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
target_mat <- OUT$mat_exp
label_sc <- OUT$label


#### ref 2000 | target 2000
n.ref <- 2000
n.target <- 2000
set.seed(1234)

cells.ref <- sample(rownames(ref_mat),2000)
cells.target <- sample(rownames(target_mat),2000)

exp_sc_mat <- exp_Habib[,cells.sample]
label_sc <- label_Habib[cells.sample,]
