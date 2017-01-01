# generate features for various window sizes
source("correlation_features.R")
source("sample_data.R")
library(DMwR)
library(dplyr)


get_target_from_id <- function(filename) {
  s_base_filename <- basename(filename)
  v_filename_parts <- strsplit(s_base_filename, "_")[[1]]
  s_target <- v_filename_parts[3]
  s_target <- gsub(".mat", "", s_target)
  n_seg_num <- as.numeric(v_filename_parts[2])
  as.numeric(s_target)
}

ensemble_modeling.1.rf.train <- function(window_size = 30, quick=T) {
  # TODO: Grid Search
  # http://machinelearningmastery.com/tuning-machine-learning-models-using-the-caret-r-package/
  
  features_filename = sprintf("../data/features/corr_fft_basicstats.20161202/train_1_window_%s_secs_correlation_and_fft.testing.txt", 
                            window_size)
  train_1_preds_filename = sprintf("../data/features/train_1_window_%s_secs_correlation_and_fft.predictions.txt", 
                              window_size)
  save_model_filename="../data/models/train_1_correlation_and_fft.random_forest.rds"
  
  set.seed(1234)
  ntree=10
  
  print(sprintf("Loading: %s", features_filename))
  trainset <- load_window_features(output_filename=features_filename)
  
  # remove rows that are all None
  trainset <- trainset[rowSums(is.na(trainset)) == 0,]
  trainset$target <- factor(trainset$target, levels=c(0,1), labels=c('interictal', 'preictal'))
  
  if (quick){
    # Sample for testing!!! build models faster...
    trainset <- sample_data(trainset, n_neg_samples=1000, n_pos_samples=500)
  }
  
  # get the metadata cols, 
  trainset.2 <- subset(trainset, select=c(id, window_id, segnum, n_dropout_rows, target))
  # remove metadata cols
  trainset <- subset(trainset, select=-c(id, window_id, segnum, n_dropout_rows))
  
  # train/test split
  inTrain = createDataPartition(trainset$target, p = 3/4)[[1]]
  training = trainset[ inTrain,]
  testing = trainset[-inTrain,]
  
  # smote
  smote_train <- SMOTE(target ~ ., data=training)
  print("smote info")
  print(table(smote_train$target))
  
  print("Training RF")
  fit_rf <- train(target ~ ., data=smote_train, method="rf", 
                  ntree=ntree,
                  trControl = trainControl(allowParallel=T, method="cv", number=4))
  prediction_rf <- predict(fit_rf, testing)
  print("Predictions RF")
  confuse_rf <- confusionMatrix(prediction_rf, testing$target, positive='preictal')
  print(confuse_rf)
  
  # Train with all the DATA with SMOTE
  smote_trainset <- SMOTE(target ~ ., data=trainset)
  fit_rf <- train(target ~ ., 
                  data=smote_trainset, 
                  method="rf", 
                  ntree=ntree,
                  trControl = trainControl(allowParallel=T, method="cv", number=4))
  
  trainset.2$prediction <- predict(fit_rf, trainset)
  confuse_rf <- confusionMatrix(trainset.2$prediction, trainset$target, positive='preictal')
  print(confuse_rf)
  
  # trainset <- subset(trainset, select=-c(id, window_id, segnum, fft_phase_entropy, n_dropout_rows))
  write.csv(trainset.2, train_1_preds_filename, row.names = F)
  print(sprintf("Wrote : %s", train_1_preds_filename))
  
  print(sprintf("Saving RDS: %s", save_model_filename))
  saveRDS(fit_rf, save_model_filename)
  train_1_preds_filename

  # Train with 
}


ensemble.2.glm <- function(window_size = 30) {
  train_1_preds_filename = sprintf("../data/features/train_1_window_%s_secs_correlation_and_fft.predictions.txt", 
                             window_size)
  train_2_filename = sprintf("../data/features/train_1_window_%s_secs_correlation_and_fft.train_2.txt", 
                             window_size)
  save_model_filename="../data/models/train_1_correlation_and_fft.glm_2.rds"
  trainset.2 <- read.csv(train_1_preds_filename, stringsAsFactors = F)
  
  # group by 
  trainset.2 <- summarize(group_by(trainset.2, id), preictal_count=sum(prediction=='preictal'), interictal_count=sum(prediction=='interictal'))
  trainset.2$target <- sapply(trainset.2$id, get_target_from_id)
  trainset.2$target <- factor(trainset.2$target, levels=c(0,1), labels=c('interictal', 'preictal'))
  
  # write features with IDs
  print(sprintf("Wrote to : %s", train_2_filename))
  write.csv(trainset.2, train_2_filename, row.names=F)
  head(trainset.2)
  
  trainset.2 <- subset(trainset.2, select=-c(id))
  inTrain = createDataPartition(trainset.2$target, p = 3/4)[[1]]
  training = trainset.2[ inTrain,]
  testing = trainset.2[-inTrain,]
  print(sprintf("training dim: %s", dim(training)[1] ))
  print(sprintf("testing dim: %s", dim(testing)[1] ))
  
  print (table(training$target))
  print (table(testing$target))
  # 
  
  
  print("Training Logistic Regression")
  fit_glm <- train(target ~ .,
                   data=training,
                   method="glm", family="binomial", trControl = trainControl(allowParallel=T, method="cv", number=4))
  prediction_glm <- predict(fit_glm, testing)
  print("Predictions GLM")
  confuse_glm <- confusionMatrix(prediction_glm, testing$target, positive='preictal')
  print(confuse_glm)

  # TODO: train a model on all the training data...  
  # print(sprintf("Saving RDS: %s", save_model_filename))
  # saveRDS(fit_glm, save_model_filename)
  
}

