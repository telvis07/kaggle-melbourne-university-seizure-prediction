source("correlation_features.R")
source("sample_data.R")
source("utils.R")
library(DMwR)
library(dplyr)
library(ggplot2)
library(pROC)

# getModelInfo("glm")

group_meta_trainset_by_id <- function(df_meta) {
  # group training by id (filename)
  df_meta.byid <- summarize(group_by(df_meta, id), preictal_count=sum(M1=='preictal'), 
                                interictal_count=sum(M1=='interictal'),
                                total=n())
  df_meta.byid$target <- sapply(df_meta.byid$id, get_target_from_id)
  df_meta.byid$target <- factor(df_meta.byid$target, levels=c(0,1), labels=c('interictal', 'preictal'))
  
  df_meta.byid 
}


glm.kfold.grid_search_parallel <- function(quick=T, n_folds=5, cores=4){
  set.seed(1234)
  
  v_window_size <- c(30)
  v_weights <- c(1, 0.75, 0.5, 0.25, 0.12, 0.06, 0.03, 0.01)
    
  tasks <- NULL
  
  i <- 1
  for (weights in v_weights) {
      print(sprintf("weights: %s",weights))
      for (window_size in v_window_size){
        task <- list(window_size=window_size, quick=quick, n_folds=n_folds, weight_preictal=1, weight_interictal=weights)
        tasks <- c(tasks, list(task)) 
      }
  }
  
  print(sprintf("Number of tasks: %s", length(tasks)))
  
  runtime <- system.time({
    all_stats_df <- mclapply(tasks, glm.kfold.grid_search, mc.cores = cores, mc.silent = FALSE)
  })[3]
  
  print(sprintf("runtime : %s", runtime))
  print(sprintf("length : %s", length(all_stats_df)))
  all_stats_df <- do.call("rbind", all_stats_df)
  print(all_stats_df)
  
  names(all_stats_df) <- make.names(names(all_stats_df))
  print("Writing glm_grid_search.csv")
  write.csv(all_stats_df, "glm_grid_search.csv", row.names = F)
  
  # sort by F1, recall
  ord <- order(all_stats_df$F1, all_stats_df$recall, decreasing = T)
  all_stats_df <- all_stats_df[ord,]
  all_stats_df
}


glm.kfold.grid_search <- function(task){
  # task <- list(window_size=30, n_folds=5, quick=T, weight_preictal=1, weight_interictal=0.5)
  set.seed(1234)

  #
  window_size <- task$window_size
  n_folds <- task$n_folds
  quick <- task$quick
  weights <- task$weights
  weight_interictal = task$weight_interictal
  weight_preictal = task$weight_preictal

  print(sprintf("[glm.kfold.grid_search] task: %s",task))
  print(sprintf("[glm.kfold.grid_search] weights$weight_interictal: %s, weights$weight_preictal: %s",weight_interictal, weight_preictal))
  
  # load the META training data
  kfold_trainfile <- sprintf("train_meta_1_MODEL_1_window_%s_secs_correlation_and_fft.kfold.quick_%s.txt", window_size, quick)
  trainset <- read.csv(kfold_trainfile, header=T, stringsAsFactors = F)
  print(table(trainset$fold))
  v_all_folds = seq(n_folds)
  
  print(sprintf("[glm.kfold.grid_search] GRID Params: weights: %s, window_size=%s", weights, window_size))
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
    
    
    df_testing.byid <- group_meta_trainset_by_id(df_testing)
    print(sprintf("df_testing.byid nrow: %s", nrow(df_testing.byid)))
    
    df_training.byid <- group_meta_trainset_by_id(df_training)
    print(sprintf("df_training.byid nrow: %s", nrow(df_training.byid)))
    
    #################################
    
    train_weights <- sapply(df_training.byid$target, function(x){
      ifelse(x=='preictal', weight_preictal, weight_interictal)
    })


    # train a glm model
    fit_glm <- train(target ~ preictal_count + interictal_count,
                     data= df_training.byid,
                     preProcess=c("center","scale"),
                     method="glm", family="binomial", weights=train_weights)
    
    # calc stats
    predictions <- predict(fit_glm, df_testing.byid)
    cm = as.matrix(table("Actual"=df_testing.byid$target, "Predicted"=predictions))
    df_stats <- calc_perf_stats(cm)
    print(df_stats)

    all_stats_df <- rbind(all_stats_df, df_stats)
  }
  

  
  # # get the average precision, recall, F1
  all_stats_df <- data.frame(precision=mean(all_stats_df$precision),
                             recall=mean(all_stats_df$recall),
                             F1=mean(all_stats_df$F1),
                             weight_preictal=weight_preictal,
                             weight_interictal=weight_interictal,
                             window_size=window_size)
  # 
  print(all_stats_df)
  names(all_stats_df) <- make.names(names(all_stats_df))
  all_stats_df
}

glm.train_full <- function(window_size=30, quick=T, weight_interictal=0.25, weight_preictal=1.0) {
  # precision     recall         F1 weight_preictal weight_interictal window_size
  # 4 0.3019874 0.40059894 0.34261546               1              0.25          30
  set.seed(1234)
  
  
  # load the META training data
  kfold_trainfile <- sprintf("train_meta_1_MODEL_1_window_%s_secs_correlation_and_fft.kfold.quick_%s.txt", window_size, quick)
  trainset <- read.csv(kfold_trainfile, header=T, stringsAsFactors = F)
  
  # group trainset by id (filename)
  trainset.byid <- group_meta_trainset_by_id(trainset)
  print(sprintf("trainset.byid nrow: %s", nrow(trainset.byid)))
  
  # train the GLM model using the parameters found during the grid search
  train_weights <- sapply(trainset.byid$target, function(x){
    ifelse(x=='preictal', weight_preictal, weight_interictal)
  })
  
  
  # train a glm model
  fit_glm <- train(target ~ preictal_count + interictal_count,
                   data= trainset.byid,
                   preProcess=c("center","scale"),
                   method="glm", family="binomial", weights=train_weights)
  
  # get performance on training datums. hopefully not overfit. :-)
  predictions <- predict(fit_glm, trainset.byid)
  cm = as.matrix(table("Actual"=trainset.byid$target, "Predicted"=predictions))
  df_stats <- calc_perf_stats(cm)
  
  # store the model
  save_model_filename=sprintf("../data/models/train_1_window_%s_secs_correlation_and_fft.quick_%s.glm_ensemble.rds", window_size, quick)
  print(sprintf("Saving RDS: %s", save_model_filename))
  saveRDS(fit_glm, save_model_filename)
  
  # 
  df_stats
}

show_grid_search_results.glm <- function(filename="../data/models/glm_grid_search.csv") {
  # > head(show_grid_search_results.glm())
  # precision     recall         F1 weight_preictal weight_interictal window_size
  # 4 0.3019874 0.40059894 0.34261546               1              0.25          30
  # 5 0.1342071 0.95258955 0.23510198               1              0.12          30
  # 6 0.1098255 1.00000000 0.19791437               1              0.06          30
  # 7 0.1094934 1.00000000 0.19737521               1              0.03          30
  # 8 0.1094934 1.00000000 0.19737521               1              0.01          30
  # 3 0.2944444 0.02598943 0.04743746               1              0.50          30
  
  
  df_stats <- read.csv(filename, header=T, stringsAsFactors = F)
  ord <- order(df_stats$F1, df_stats$recall, decreasing = T)
  df_stats[ord,]
  
}



