source("correlation_features.R")
source("sample_data.R")
source("utils.R")
library(DMwR)
library(dplyr)
library(ggplot2)
library(pROC)

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
  if (quick) {
    v_window_size <- c(10, 30)
    v_rf_ntree <- c(5, 10)
    v_rf_mtry <- c(2,74)
    v_rf_metrics <- c("Kappa")
    
  } else {
    v_window_size <- c(5, 10, 30, 60, 120, 300)
    v_rf_ntree <- c(2, 10, 30, 50, 101, 300, 500, 1000)
    v_rf_mtry <- c(2,38,74,110,300)
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
  
  print(sprintf("GRID Params: metric: %s, ntree: %s, mtry: %s", metric, ntree, mtry))
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

  names(all_stats_df) <- make.names(names(all_stats_df))
  all_stats_df
  
  
  
  # 
  # 
  # 
  # # convert to data.frame
  # df_grid_parameters <- mutate(df_grid_parameters, rf_ntree=as.numeric(rf_ntree), rf_mtry=as.numeric(rf_mtry))
  # all_stats_df <- data.frame()
  # 
  # # for each rf parameter
  # for (i in seq(n_parameters)){
  #     # for each test fold
  #     for (n_fold in seq(n_folds)){
  #       print(sprintf("parameters: %s", paste(df_grid_parameters[i,], collapse =",")))
  #       metric <- df_grid_parameters[i,]$rf_metric
  #       ntree <- df_grid_parameters[i,]$rf_ntree
  #       mtry <- df_grid_parameters[i,]$rf_mtry
  # 
  #       print(sprintf("GRID Params: metric: %s, ntree: %s, mtry: %s", metric, ntree, mtry))
  #       v_training_folds <- setdiff(v_all_folds, c(n_fold))
  #       v_test_fold <- c(n_fold)
  #       print(sprintf("test_fold: %s, training_folds: %s", paste(v_training_folds, collapse=","), paste(v_test_fold, collapse = ",")))
  #       
  #       # test on this fold
  #       df_testing <- filter(trainset, fold %in% v_test_fold)
  #       
  #       # Combine the other four folds to be used as a training fold
  #       df_training <- filter(trainset, fold %in% v_training_folds)
  # 
  #       print(sprintf("train_fold length: %s, test_fold_length: %s", dim(df_training)[1], dim(df_testing)[1]))
  #       
  #       # cleanup cleanup rows and smote
  #       df_training$target <- factor(df_training$target, levels=c('interictal', 'preictal'))
  #       df_training <- subset(df_training, select=-c(id, window_id, segnum, n_dropout_rows, fold, rowid))
  #       smote_train <- SMOTE(target ~ ., data=df_training)
  # 
  #       # train a random forest model
  #       fit_rf <- train(target ~ ., data=smote_train, method="rf", 
  #                       ntree=ntree,
  #                       tuneGrid=expand.grid(mtry=c(mtry)),
  #                       metric=metric)
  #       
  #       # calc stats
  #       prediction_rf <- predict(fit_rf, df_testing)
  #       cm = as.matrix(table("Actual"=df_testing$target, "Predicted"=prediction_rf))
  #       n = sum(cm) # number of instances
  #       nc = nrow(cm) # number of classes
  #       diag = diag(cm) # number of correctly classified instances per class 
  #       rowsums = apply(cm, 1, sum) # number of instances per class
  #       colsums = apply(cm, 2, sum) # number of predictions per class
  #       p = rowsums / n # distribution of instances over the actual classes
  #       q = colsums / n # distribution of instances over the predicted classes
  #       
  #       # accuracy
  #       accuracy = sum(diag) / n
  #       
  #       # precision/recall/f1
  #       precision = diag / colsums 
  #       recall = diag / rowsums 
  #       f1 = 2 * precision * recall / (precision + recall) 
  #       
  # 
  #       df_stats <- data.frame(precision=precision[2], 
  #                              recall=recall[2], 
  #                              F1=f1[2],
  #                              ntree=ntree,
  #                              mtry=mtry,
  #                              metric=metric,
  #                              test_fold=n_fold) 
  #       
  #       all_stats_df <- rbind(all_stats_df, df_stats)
  # 
  #     }
  # }
  # 
  # 
  # names(all_stats_df) <- make.names(names(all_stats_df))
  # print("Writing rf_grid_search.csv")
  # write.csv(all_stats_df, "rf_grid_search.csv", row.names = F)
  # ord <- order(all_stats_df$F1, all_stats_df$recall, decreasing = T)
  # all_stats_df
}

# ensemble_modeling.1.rf.grid_search <- function(window_size = 30, quick=T) {
#   # TODO: Grid Search
#   # http://machinelearningmastery.com/tuning-machine-learning-models-using-the-caret-r-package/
#   
#   # getModelInfo("rf")
#   # attributes(fit_rf)
# 
#   features_filename = sprintf("../data/features/corr_fft_basicstats.20161202/train_1_window_%s_secs_correlation_and_fft.testing.txt", 
#                               window_size)
#   train_1_preds_filename = sprintf("../data/features/train_1_window_%s_secs_correlation_and_fft.predictions.txt", 
#                                    window_size)
#   save_model_filename=sprintf("../data/models/train_1_window_%s_secs_correlation_and_fft.random_forest.rds", window_size)
#   
#   set.seed(1234)
#   
#   print(sprintf("Loading: %s", features_filename))
#   trainset <- load_window_features(output_filename=features_filename)
#   
#   # remove rows that are all None
#   trainset <- trainset[rowSums(is.na(trainset)) == 0,]
#   trainset$target <- factor(trainset$target, levels=c(0,1), labels=c('interictal', 'preictal'))
#   
#   if (quick){
#     # Sample for testing!!! build models faster...
#     trainset <- sample_data(trainset, n_neg_samples=1000, n_pos_samples=500)
#     ntree=10
#   } else {
#     ntree = 100
#   }
#   
#   # get the metadata cols, 
#   trainset.2 <- subset(trainset, select=c(id, window_id, segnum, n_dropout_rows, target))
#   # remove metadata cols
#   trainset <- subset(trainset, select=-c(id, window_id, segnum, n_dropout_rows))
#   
#   # train/test split
#   inTrain = createDataPartition(trainset$target, p = 3/4)[[1]]
#   training = trainset[ inTrain,]
#   testing = trainset[-inTrain,]
#   
#   # smote
#   smote_train <- SMOTE(target ~ ., data=training)
#   print("smote info")
#   print(table(smote_train$target))
#   
#   print("Training RF using test/train split")
#   all_stats_df <- data.frame()
#   # ntree=c(5,10,50,100,500)
#   
#   for (metric in c("Accuracy", "Kappa")){
#     for (ntree in c(5, 10)) {
#       
#       # mtry
#       # Number of variables randomly sampled as candidates at each split. 
#       # Note that the default values are different for 
#       # classification (sqrt(p) where p is number of variables in x) 
#       # and regression (p/3)
#       
#       # mtry=c(2,38,74,110,300)
#       grid <- expand.grid(mtry=c(2,110))
#       
#       fit_rf <- train(target ~ ., data=smote_train, method="rf", 
#                       ntree=ntree,
#                       tuneGrid=grid,
#                       metric=metric,
#                       trControl = trainControl(allowParallel=T, method="cv", number=4))
#       prediction_rf <- predict(fit_rf, testing)
#       cm = as.matrix(table("Actual"=testing$target, "Predicted"=prediction_rf))
#       n = sum(cm) # number of instances
#       nc = nrow(cm) # number of classes
#       diag = diag(cm) # number of correctly classified instances per class 
#       rowsums = apply(cm, 1, sum) # number of instances per class
#       colsums = apply(cm, 2, sum) # number of predictions per class
#       p = rowsums / n # distribution of instances over the actual classes
#       q = colsums / n # distribution of instances over the predicted classes
#       
#       # accuracy
#       accuracy = sum(diag) / n
# 
#       # precision/recall/f1
#       precision = diag / colsums 
#       recall = diag / rowsums 
#       f1 = 2 * precision * recall / (precision + recall) 
#       
#       print("Predictions RF")
#       confuse_rf <- confusionMatrix(prediction_rf, testing$target, positive='preictal')
#       # df_stats <- as.data.frame(t(confuse_rf$byClass))
#       
#       print(t(confuse_rf$byClass))
# 
#       
#       # df_stats$ntree <- ntree
#       # df_stats$mtry <- fit_rf$bestTune$mtry
#       # df_stats$metric <- fit_rf$metric
#       df_stats <- data.frame(precision=precision[2], 
#                              recall=recall[2], 
#                              F1=f1[2],
#                              ntree=ntree,
#                              mtry=fit_rf$bestTune$mtry,
#                              metric=fit_rf$metric) 
#       
#       all_stats_df <- rbind(all_stats_df, df_stats)
#       
#       # print(confuse_rf)
#       # print(fit_rf)
#       # print(summary(fit_rf))
#       
#   
#     } # ntree
#   } # metric
#   
#   names(all_stats_df) <- make.names(names(all_stats_df))
#   print("Writing rf_grid_search.csv")
#   write.csv(all_stats_df, "rf_grid_search.csv", row.names = F)
#   
#   ord <- order(all_stats_df$F1, all_stats_df$recall, decreasing = T)
#   # all_stats_df[ord, c("model", "Recall", "sample_type", "F1")]
#   
#   all_stats_df
#   
#   # fit_rf <- fit_rf$finalModel
#   # 
#   # # Train with all the DATA with SMOTE
#   # # print("SMOTE all the datums")
#   # # smote_trainset <- SMOTE(target ~ ., data=trainset)
#   # # print("Train with all the DATA with SMOTE")
#   # # fit_rf <- train(target ~ ., 
#   # #                 data=smote_trainset, 
#   # #                 method="rf", 
#   # #                 ntree=ntree,
#   # #                 trControl = trainControl(allowParallel=T, method="cv", number=4))
#   # 
#   # trainset.2$prediction <- predict(fit_rf, trainset)
#   # confuse_rf <- confusionMatrix(trainset.2$prediction, trainset$target, positive='preictal')
#   # print(confuse_rf)
#   # 
#   # # trainset <- subset(trainset, select=-c(id, window_id, segnum, fft_phase_entropy, n_dropout_rows))
#   # write.csv(trainset.2, train_1_preds_filename, row.names = F)
#   # print(sprintf("Wrote : %s", train_1_preds_filename))
#   # 
#   # print(sprintf("Saving RDS: %s", save_model_filename))
#   # saveRDS(fit_rf, save_model_filename)
#   # train_1_preds_filename
#   
# }