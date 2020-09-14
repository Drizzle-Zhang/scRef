setwd('/home/zy/my_git/scRef/Benchmark/cross_validation/train4_test1')
source('./Cross_Validation.R')
source('./method_functions.R')
source('./evaluate.R')

path.input <- '/home/zy/scRef/'
path.output <- '/home/zy/scRef/cross_validation/train4_test1/'

# generate cross validation dataset
LabelsPath <- paste0(path.input, 'summary/Zeisel_exp_sc_mat_cluster_original.txt')
OutputDir <- path.output
# Cross_Validation(LabelsPath, OutputDir)

DataPath <- paste0(path.input, 'summary/Zeisel_exp_sc_mat.txt')
LabelsPath <- paste0(path.input, 'summary/Zeisel_exp_sc_mat_cluster_original.txt')
CV_RDataPath <- paste0(path.output, 'CV_folds.RData')

# SingleR
run_SingleR(DataPath,LabelsPath,CV_RDataPath,OutputDir)
TrueLabelsPath <- paste0(OutputDir, 'SingleR_True_Labels.csv')
PredLabelsPath <- paste0(OutputDir, 'SingleR_Pred_Labels.csv')
res.SingleR <- evaluate(TrueLabelsPath, PredLabelsPath)

# scmap
run_scmap(DataPath,LabelsPath,CV_RDataPath,OutputDir)
TrueLabelsPath <- paste0(OutputDir, 'scmapcell_True_Labels.csv')
PredLabelsPath <- paste0(OutputDir, 'scmapcell_Pred_Labels.csv')
res.scmapcell <- evaluate(TrueLabelsPath, PredLabelsPath)
PredLabelsPath <- paste0(OutputDir, 'scmapcluster_Pred_Labels.csv')
res.scmapcluster <- evaluate(TrueLabelsPath, PredLabelsPath)

# CHETAH
run_CHETAH(DataPath,LabelsPath,CV_RDataPath,OutputDir)
TrueLabelsPath <- paste0(OutputDir, 'CHETAH_True_Labels.csv')
PredLabelsPath <- paste0(OutputDir, 'CHETAH_Pred_Labels.csv')
res.CHETAH <- evaluate(TrueLabelsPath, PredLabelsPath)

# scPred
run_scPred(DataPath,LabelsPath,CV_RDataPath,OutputDir)
TrueLabelsPath <- paste0(OutputDir, 'scPred_True_Labels.csv')
PredLabelsPath <- paste0(OutputDir, 'scPred_Pred_Labels.csv')
res.scPred <- evaluate(TrueLabelsPath, PredLabelsPath)

# sciBet
run_sciBet(DataPath,LabelsPath,CV_RDataPath,OutputDir)
TrueLabelsPath <- paste0(OutputDir, 'sciBet_True_Labels.csv')
PredLabelsPath <- paste0(OutputDir, 'sciBet_Pred_Labels.csv')
res.sciBet <- evaluate(TrueLabelsPath, PredLabelsPath)

# scRef
run_scRef(DataPath,LabelsPath,CV_RDataPath,OutputDir)
TrueLabelsPath <- paste0(OutputDir, 'scRef_True_Labels.csv')
PredLabelsPath <- paste0(OutputDir, 'scRef_Pred_Labels_cell.csv')
# PredLabelsPath <- paste0(OutputDir, 'scRef_Pred_Labels.csv')
res.scRef <- evaluate(TrueLabelsPath, PredLabelsPath)

# singleCellNet
run_singleCellNet(DataPath,LabelsPath,CV_RDataPath,OutputDir)
TrueLabelsPath <- paste0(OutputDir, 'singleCellNet_True_Labels.csv')
PredLabelsPath <- paste0(OutputDir, 'singleCellNet_Pred_Labels.csv')
res.singleCellNet <- evaluate(TrueLabelsPath, PredLabelsPath)

# CaSTLe


