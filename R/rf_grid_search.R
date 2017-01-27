source("correlation_features.R")
source("sample_data.R")
source("utils.R")
library(DMwR)
library(dplyr)
library(ggplot2)
library(pROC)

# getModelInfo("rf")


# TODO
# In stratified k-fold cross-validation, the folds are 
# selected so that the mean response value is approximately 
# equal in all the folds. In the case of a dichotomous classification, 
# this means that each fold contains roughly the same proportions of the two types of class labels.
# https://gist.github.com/multidis/8093372

# http://blog.kaggle.com/2016/12/27/a-kagglers-guide-to-model-stacking-in-practice/
# Let’s make a sad attempt at solving this classification problem using a K-Nearest Neighbors model. In order to select the best value for K, we’ll use 5-fold Cross-Validation combined with Grid Search where K=(1, 2, … 30). In pseudo code:

# 1. Partition the training data into five equal size folds. Call these test folds.
# 2. For K = 1, 2, … 10
#   1. For each test fold
#     1. Combine the other four folds to be used as a training fold
#     2. Fit a K-Nearest Neighbors model on the training fold (using the current value of K)
#     3. Make predictions on the test fold and measure the resulting accuracy rate of the predictions
#   2. Calculate the average accuracy rate from the five test fold predictions
# 3. Keep the K value with the best average CV accuracy rate
ensemble_modeling.kfoldid <- function(window_size = 30, quick=F, n_folds=5) {
  # trainset <- ensemble_modeling.kfoldid(window_size = 300, quick=F, n_folds=5)
  # trainset <- ensemble_modeling.kfoldid(window_size = 120, quick=F, n_folds=5)
  # trainset <- ensemble_modeling.kfoldid(window_size = 60, quick=F, n_folds=5)
  # trainset <- ensemble_modeling.kfoldid(window_size = 30, quick=F, n_folds=5)
  # trainset <- ensemble_modeling.kfoldid(window_size = 10, quick=F, n_folds=5)
  # trainset <- ensemble_modeling.kfoldid(window_size = 5, quick=F, n_folds=5)
  
  # trainset <- ensemble_modeling.kfoldid(window_size = 30, quick=T, n_folds=5)
  # trainset <- ensemble_modeling.kfoldid(window_size = 10, quick=T, n_folds=5)
  # trainset <- ensemble_modeling.kfoldid(window_size = 5, quick=T, n_folds=5)
  set.seed(1234)
  features_filename = sprintf("../data/features/corr_fft_basicstats.20161202/train_1_window_%s_secs_correlation_and_fft.testing.txt", 
                              window_size)
  trainset <- load_window_features(output_filename=features_filename)
  
  # remove rows that are all None
  trainset <- trainset[rowSums(is.na(trainset)) == 0,]
  trainset$target <- factor(trainset$target, levels=c(0,1), labels=c('interictal', 'preictal'))
  
  if (quick){
    # Sample for testing!!! build models faster...
    trainset <- sample_data(trainset, n_neg_samples=1000, n_pos_samples=500)
    ntree=10
  }
  folds <- createFolds(y=trainset$target, k=n_folds, list=TRUE,returnTrain=FALSE)
  print(sapply(folds, length))
  
  # add row ids 
  trainset$rowid <- seq(nrow(trainset))
  trainset$fold <- NA
  
  # add folds
  for (i in seq(length(folds))){
    vn_fold <- folds[[i]]
    print(sprintf("Fold: %s, size=%s", i, length(vn_fold)))
    trainset[vn_fold,]$fold <- i
  }
  
  outfile <- sprintf("train_1_window_%s_secs_correlation_and_fft.kfold.quick_%s.txt", window_size, quick)
  print(sprintf("Writing kfold data to : %s", outfile))
  write.csv(trainset, outfile, row.names = F)
  trainset
}

rf.kfold.grid_search_parallel <- function(quick=T, n_folds=5, cores=4){
  # See : getModelInfo("rf")

  set.seed(1234)
  
  if (quick) {
    v_window_size <- c(10, 30)
    v_rf_ntree <- c(5, 10)
    v_rf_mtry <- c(2,74)
    v_rf_metrics <- c("Kappa")
    
  } else {
    # v_window_size <- c(5, 10, 30, 60, 120, 300)
    v_window_size <- c(30, 60, 120, 300)
    v_rf_ntree <- c(2, 10, 50, 100, 300)
    v_rf_mtry <- c(2,38,74,147)
    v_rf_metrics <- c("Accuracy", "Kappa")
  }
  
  
  n_parameters <- length(v_rf_metrics) * length(v_rf_ntree) * length(v_rf_mtry)
  tasks <- NULL
  
  i <- 1
  for (metric in v_rf_metrics) {
    for (ntree in v_rf_ntree) { 
      for (mtry in v_rf_mtry) {
        for (window_size in v_window_size){
          task <- list(window_size=window_size, quick=quick, n_folds=n_folds, metric=metric, ntree=ntree, mtry=mtry)
          tasks <- c(tasks, list(task)) 
        }
      }
    }
  }
  
  print(sprintf("Number of tasks: %s", length(tasks)))
  
  runtime <- system.time({
    all_stats_df <- mclapply(tasks, rf.kfold.grid_search, mc.cores = cores, mc.silent = FALSE)
  })[3]
  
  print(sprintf("runtime : %s", runtime))
  print(sprintf("length : %s", length(all_stats_df)))
  all_stats_df <- do.call("rbind", all_stats_df)

  names(all_stats_df) <- make.names(names(all_stats_df))
  print("Writing rf_grid_search.csv")
  write.csv(all_stats_df, "rf_grid_search.csv", row.names = F)
  
  # sort by F1, recall
  ord <- order(all_stats_df$F1, all_stats_df$recall, decreasing = T)
  all_stats_df <- all_stats_df[ord,]
  all_stats_df
  
}


rf.kfold.grid_search <- function(task){
  # See: getModelInfo("rf")

  set.seed(1234)
  
  # Expects to be called by rf.kfold.grid_search_parallel
  # task has the parameters below:
  
  #
  window_size <- task$window_size
  n_folds <- task$n_folds
  quick <- task$quick
  metric <- task$metric
  ntree <- task$ntree
  mtry <- task$mtry
  
  # 
  kfold_trainfile <- sprintf("train_1_window_%s_secs_correlation_and_fft.kfold.quick_%s.txt", window_size, quick)
  trainset <- read.csv(kfold_trainfile, header=T, stringsAsFactors = F)
  print(table(trainset$fold))
  v_all_folds = seq(n_folds)
  
  print(sprintf("GRID Params: metric: %s, ntree: %s, mtry: %s, window_size: %s", metric, ntree, mtry, window_size))
  all_stats_df <- data.frame()
  
  for (n_fold in seq(n_folds)){
    v_training_folds <- setdiff(v_all_folds, c(n_fold))
    v_test_fold <- c(n_fold)
    print(sprintf("test_fold: %s, training_folds: %s", paste(v_training_folds, collapse=","), paste(v_test_fold, collapse = ",")))
  
    # test on this fold
    df_testing <- filter(trainset, fold %in% v_test_fold)
  
    # Combine the other four folds to be used as a training fold
    df_training <- filter(trainset, fold %in% v_training_folds)
  
    print(sprintf("train_fold length: %s, test_fold_length: %s", dim(df_training)[1], dim(df_testing)[1]))
  
    # cleanup cleanup rows and smote
    df_training$target <- factor(df_training$target, levels=c('interictal', 'preictal'))
    df_training <- subset(df_training, select=-c(id, window_id, segnum, n_dropout_rows, fold, rowid))
    smote_train <- SMOTE(target ~ ., data=df_training)
  
    # train a random forest model
    fit_rf <- train(target ~ ., data=smote_train, method="rf",
                    ntree=ntree,
                    tuneGrid=expand.grid(mtry=c(mtry)),
                    metric=metric)
  
    # calc stats
    prediction_rf <- predict(fit_rf, df_testing)
    cm = as.matrix(table("Actual"=df_testing$target, "Predicted"=prediction_rf))
    n = sum(cm) # number of instances
    nc = nrow(cm) # number of classes
    diag = diag(cm) # number of correctly classified instances per class
    rowsums = apply(cm, 1, sum) # number of instances per class
    colsums = apply(cm, 2, sum) # number of predictions per class
    p = rowsums / n # distribution of instances over the actual classes
    q = colsums / n # distribution of instances over the predicted classes
  
    # accuracy
    accuracy = sum(diag) / n
  
    # precision/recall/f1
    precision = diag / colsums
    recall = diag / rowsums
    f1 = 2 * precision * recall / (precision + recall)
  
  
    df_stats <- data.frame(precision=precision[2],
                           recall=recall[2],
                           F1=f1[2])
    all_stats_df <- rbind(all_stats_df, df_stats)
  }
  
  # get the average precision, recall, F1
  all_stats_df <- data.frame(precision=mean(all_stats_df$precision),
                             recall=mean(all_stats_df$recall),
                             F1=mean(all_stats_df$F1),
                             ntree=ntree,
                             mtry=mtry,
                             metric=metric,
                             window_size=window_size)
  
  print(all_stats_df)

  names(all_stats_df) <- make.names(names(all_stats_df))
  all_stats_df
  
}

rf.train_full <- function(window_size=30, metric="Accuracy", ntree=50, mtry=74, quick=T) {
  # train the RF model using the parameters found during the grid search (see show_grid_search_results())
  #     precision    recall        F1 ntree mtry   metric window_size
  # 41  0.2722213 0.4697898 0.3443397    50   74 Accuracy          30
  
  # train the RF model using the parameters found during the grid search
  set.seed(1234)
  
  # load the datums with kfold IDs
  kfold_trainfile <- sprintf("train_1_window_%s_secs_correlation_and_fft.kfold.quick_%s.txt", window_size, quick)
  trainset <- read.csv(kfold_trainfile, header=T, stringsAsFactors = F)
  print(table(trainset$fold))
  
  # cleanup cleanup rows
  trainset$target <- factor(trainset$target, levels=c('interictal', 'preictal'))
  trainset <- subset(trainset, select=-c(id, window_id, segnum, n_dropout_rows, fold, rowid))
  
  # Train with all the DATA with SMOTE
  print("SMOTE all the datums")
  smote_trainset <- SMOTE(target ~ ., data=trainset)
  print("Train with all the DATA with SMOTE")
  fit_rf <- train(target ~ .,
                  data=smote_trainset,
                  method="rf",
                  ntree=ntree,
                  metric=metric,
                  tuneGrid=expand.grid(mtry=c(mtry)))
  
  # get performance on training datums. hopefully not overfit. :-)
  prediction_rf <- predict(fit_rf, trainset)
  cm = as.matrix(table("Actual"=trainset$target, "Predicted"=prediction_rf))
  df_stats <- calc_perf_stats(cm)
  
  # store the model
  save_model_filename=sprintf("../data/models/train_1_window_%s_secs_correlation_and_fft.quick_%s.random_forest.rds", window_size, quick)
  print(sprintf("Saving RDS: %s", save_model_filename))
  saveRDS(fit_rf, save_model_filename)
  
  # 
  df_stats
}

create_meta_trainset <- function(window_size=30, metric="Accuracy", ntree=50, mtry=74, n_folds=5, quick=T) {
  # For each test fold
  #   train a model on 4 other folds, test on left out fold
  #   store predictions in trainset_meta
  # Store train_meta
  #
  # Parameters: 
  # train the RF model using the parameters found during the grid search (see show_grid_search_results())
  #     precision    recall        F1 ntree mtry   metric window_size
  # 41  0.2722213 0.4697898 0.3443397    50   74 Accuracy          30
  set.seed(1234)
  
  # load the datums with kfold IDs
  kfold_trainfile <- sprintf("train_1_window_%s_secs_correlation_and_fft.kfold.quick_%s.txt", window_size, quick)
  trainset <- read.csv(kfold_trainfile, header=T, stringsAsFactors = F)
  print(table(trainset$fold))
  v_all_folds = seq(n_folds)
  
  # cleanup cleanup rows
  trainset$target <- factor(trainset$target, levels=c('interictal', 'preictal'))
  trainset_meta <- subset(trainset, select=c(id, fold, rowid, target))
  trainset_meta <- mutate(trainset_meta, M1=NA)
  
  for (n_fold in seq(n_folds)){
    v_training_folds <- setdiff(v_all_folds, c(n_fold))
    v_test_fold <- c(n_fold)
    print(sprintf("test_fold: %s, training_folds: %s", paste(v_training_folds, collapse=","), paste(v_test_fold, collapse = ",")))
    
    # test on this fold
    df_testing <- filter(trainset, fold %in% v_test_fold)
    
    # Combine the other four folds to be used as a training fold
    df_training <- filter(trainset, fold %in% v_training_folds)
    
    print(sprintf("train_fold length: %s, test_fold_length: %s", dim(df_training)[1], dim(df_testing)[1]))
    
    # cleanup cleanup rows and smote
    df_training$target <- factor(df_training$target, levels=c('interictal', 'preictal'))
    df_training <- subset(df_training, select=-c(id, window_id, segnum, n_dropout_rows, fold, rowid))
    smote_train <- SMOTE(target ~ ., data=df_training)
    
    # train a random forest model
    fit_rf <- train(target ~ .,
                    data=smote_train,
                    method="rf",
                    ntree=ntree,
                    metric=metric,
                    tuneGrid=expand.grid(mtry=c(mtry)))
    
    # calc stats
    df_testing$M1 <- predict(fit_rf, df_testing)
    cm = as.matrix(table("Actual"=df_testing$target, "Predicted"=df_testing$M1))
    df_stats <- calc_perf_stats(cm)
    print(df_stats)
    
    # 
    trainset_meta <- merge(trainset_meta, subset(df_testing, select=c(rowid, M1)), by="rowid", all.x=T)
    
    # merge the M1 field, x = trainset_meta, y = df_testing. 
    # Use df_testing.M1 if trainset_meta.M1 is null
    trainset_meta <- mutate(trainset_meta, M1=ifelse(is.na(M1.x), M1.y, M1.x))
    # remove merge cols (.x, .y)
    trainset_meta <- subset(trainset_meta, select=c(id, fold, rowid, target, M1))

    print("trainset_meta folds with predictions")
    print(table(trainset_meta[!is.na(trainset_meta$M1),]$fold))
  }
  
  # change the labels back to the factor strings
  trainset_meta$M1 <- factor(trainset_meta$M1, levels=c(1,2), labels=c('interictal', 'preictal'))
  
  outfile <- sprintf("train_meta_1_MODEL_1_window_%s_secs_correlation_and_fft.kfold.quick_%s.txt", window_size, quick)
  print(sprintf("Writing train_meta data to : %s", outfile))
  write.csv(trainset_meta, outfile, row.names = F)
  
  trainset_meta
  
}


show_grid_search_results.rf <- function(filename="../data/models/rf_grid_search.csv") {
  # > head(show_grid_search_results.rf())
  # precision    recall        F1 ntree mtry   metric window_size
  # 69  0.2797997 0.4558048 0.3466332   300   38 Accuracy          30
  # 149 0.2797997 0.4558048 0.3466332   300   38    Kappa          30
  # 41  0.2722213 0.4697898 0.3443397    50   74 Accuracy          30
  # 121 0.2722213 0.4697898 0.3443397    50   74    Kappa          30
  # 74  0.2706729 0.4539553 0.3389850   300   74 Accuracy          60
  # 154 0.2706729 0.4539553 0.3389850   300   74    Kappa          60
  
  
  
  df_stats <- read.csv(filename, header=T, stringsAsFactors = F)
  ord <- order(df_stats$F1, df_stats$recall, decreasing = T)
  df_stats[ord,]
  
}