setwd('/home/drizzle_zhang/my_git/scRef/Benchmark/cross_validation/train4_test1')
source('./Cross_Validation.R')
source('./method_functions.R')
source('./evaluate.R')

path.input <- '/home/drizzle_zhang/scRef/'
path.output <- '/home/drizzle_zhang/scRef/cross_validation/train4_test1/'

# generate cross validation dataset
LabelsPath <- paste0(path.input, 'summary/Zeisel_exp_sc_mat_cluster_original.txt')
OutputDir <- path.output
# Cross_Validation(LabelsPath, OutputDir)

DataPath <- paste0(path.input, 'summary/Zeisel_exp_sc_mat.txt')
LabelsPath <- paste0(path.input, 'summary/Zeisel_exp_sc_mat_cluster_original.txt')
CV_RDataPath <- paste0(path.output, 'CV_folds.RData')

# SingleR
run_SingleR(DataPath,LabelsPath,CV_RDataPath,OutputDir)

# scID
run_scmap(DataPath,LabelsPath,CV_RDataPath,OutputDir)

# CHETAH
run_CHETAH(DataPath,LabelsPath,CV_RDataPath,OutputDir)

# scPred
run_scPred(DataPath,LabelsPath,CV_RDataPath,OutputDir)

# sciBet
run_sciBet(DataPath,LabelsPath,CV_RDataPath,OutputDir)

# scRef


