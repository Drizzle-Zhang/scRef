Cross_Validation <- function(LabelsPath, OutputDir, del.cells, col_Index = 1) {
  "
  Cross_Validation
  Function returns train and test indices for 5 folds stratified across unique cell populations,
  also filter out cell populations with less than 10 cells.
  It return a 'CV_folds.RData' file which then used as input to classifiers wrappers.
  
  Parameters
  ----------
  LabelsPath : Cell population annotations file path (.tsv).
  col_Index : column index (integer) defining which level of annotation to use,
  in case of multiple cell type annotations (default is 1)
  OutputDir : Output directory defining the path of the exported file.
  "
  
  Labels <- as.matrix(read.delim(LabelsPath, row.names = 1, header = T))
  Labels <- as.vector(Labels[,col_Index])
  
  Removed_classes <- !(table(Labels) > 10)
  Cells_to_Keep <- !(is.element(Labels,names(Removed_classes)[Removed_classes]))
  Labels <- Labels[Cells_to_Keep]
  
  # Getting training and testing Folds
  library(rBayesianOptimization)
  n_folds = 5
  Folds <- KFold(Labels,nfolds = n_folds, stratified = TRUE)
  Test_Folds <- c(n_folds:1)
  Train_Idx <- list()
  Test_Idx <- list()
  for (i in c(1:length(Folds))) {
    Temp_Folds <- Folds
    Temp_Folds[Test_Folds[i]] <- NULL
    sub_train <- unlist(Folds[Test_Folds[i]])
    sub_train <- list(sub_train[!(Labels[sub_train] %in% del.cells)])
    Train_Idx[i] <- sub_train
    Test_Idx[i] <- list(unlist(Temp_Folds))
    # Train_Idx[i] <- list(unlist(Temp_Folds))
    # Test_Idx[i] <- Folds[Test_Folds[i]]
  }
  remove(Temp_Folds,i,Folds)
  # setwd(OutputDir)
  file.out <- paste0(OutputDir, 'CV_folds.RData')
  save(n_folds,Train_Idx,Test_Idx,col_Index,Cells_to_Keep,file = file.out)
}