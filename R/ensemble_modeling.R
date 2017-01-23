# generate features for various window sizes
source("correlation_features.R")
source("sample_data.R")
source("utils.R")
library(DMwR)
library(dplyr)
library(ggplot2)
library(pROC)




ensemble_modeling.1.rf.train <- function(window_size = 30, quick=T) {
  # TODO: Grid Search
  # http://machinelearningmastery.com/tuning-machine-learning-models-using-the-caret-r-package/
  features_filename = sprintf("../data/features/corr_fft_basicstats.20161202/train_1_window_%s_secs_correlation_and_fft.testing.txt", 
                            window_size)
  train_1_preds_filename = sprintf("../data/features/train_1_window_%s_secs_correlation_and_fft.predictions.txt", 
                              window_size)
  save_model_filename=sprintf("../data/models/train_1_window_%s_secs_correlation_and_fft.random_forest.rds", window_size)
  
  set.seed(1234)
  
  print(sprintf("Loading: %s", features_filename))
  trainset <- load_window_features(output_filename=features_filename)
  
  # remove rows that are all None
  trainset <- trainset[rowSums(is.na(trainset)) == 0,]
  trainset$target <- factor(trainset$target, levels=c(0,1), labels=c('interictal', 'preictal'))
  
  if (quick){
    # Sample for testing!!! build models faster...
    trainset <- sample_data(trainset, n_neg_samples=1000, n_pos_samples=500)
    ntree=10
  } else {
    ntree = 100
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
  
  print("Training RF using test/train split")
  grid <- expand.grid(mtry=c(2,38,74,110,300), ntree=c(5,10,50,100,500))
  fit_rf <- train(target ~ ., data=smote_train, method="rf", 
                  # ntree=ntree,
                  # tuneLength=10,
                  trControl = trainControl(allowParallel=T, method="cv", number=4))
  prediction_rf <- predict(fit_rf, testing)
  print("Predictions RF")
  confuse_rf <- confusionMatrix(prediction_rf, testing$target, positive='preictal')
  print(confuse_rf)
  print(fit_rf)
  print(summary(fit_rf))
  
  fit_rf <- fit_rf$finalModel
  
  # Train with all the DATA with SMOTE
  # print("SMOTE all the datums")
  # smote_trainset <- SMOTE(target ~ ., data=trainset)
  # print("Train with all the DATA with SMOTE")
  # fit_rf <- train(target ~ ., 
  #                 data=smote_trainset, 
  #                 method="rf", 
  #                 ntree=ntree,
  #                 trControl = trainControl(allowParallel=T, method="cv", number=4))
  
  trainset.2$prediction <- predict(fit_rf, trainset)
  confuse_rf <- confusionMatrix(trainset.2$prediction, trainset$target, positive='preictal')
  print(confuse_rf)
  
  # trainset <- subset(trainset, select=-c(id, window_id, segnum, fft_phase_entropy, n_dropout_rows))
  write.csv(trainset.2, train_1_preds_filename, row.names = F)
  print(sprintf("Wrote : %s", train_1_preds_filename))
  
  print(sprintf("Saving RDS: %s", save_model_filename))
  saveRDS(fit_rf, save_model_filename)
  train_1_preds_filename

}


ensemble.2.glm <- function(window_size = 30, patient_num = 1) {
  set.seed(1234)
  train_1_preds_filename = sprintf("../data/ensemble/train_1_window_%s_secs_correlation_and_fft.predictions.txt", 
                              window_size)
  train_2_filename = sprintf("../data/features/train_1_window_%s_secs_correlation_and_fft.train_2.txt", 
                             window_size)
  save_model_filename=sprintf("../data/models/train_1_window_%s_secs_correlation_and_fft.glm_2.rds", window_size)
  trainset.2 <- read.csv(train_1_preds_filename, stringsAsFactors = F)
  
  
  # group by 
  trainset.2 <- summarize(group_by(trainset.2, id), preictal_count=sum(prediction=='preictal'), 
                          interictal_count=sum(prediction=='interictal'),
                          total=n())
  trainset.2$target <- sapply(trainset.2$id, get_target_from_id)
  trainset.2$target <- factor(trainset.2$target, levels=c(0,1), labels=c('interictal', 'preictal'))
  print(sprintf("trainset nrow: %s", nrow(trainset.2)))
  
  # write features with IDs
  print(sprintf("Wrote to : %s", train_2_filename))
  write.csv(trainset.2, train_2_filename, row.names=F)
  head(trainset.2)
  
  trainset.2 <- subset(trainset.2, select=-c(id))
  # trainset.2 <- filter(trainset.2, total>1)
  
  # test/train split
  inTrain = createDataPartition(trainset.2$target, p = 3/4)[[1]]
  training = trainset.2[ inTrain,]
  testing = trainset.2[-inTrain,]
  print(sprintf("training dim: %s", dim(training)[1] ))
  print(sprintf("testing dim: %s", dim(testing)[1] ))
  
  print (table(training$target))
  print (table(testing$target))
  
  # 
  print("Training Logistic Regression")
  # cctrl1 <- trainControl(method = "cv", number = 3, returnResamp = "all",
  #                        classProbs = TRUE,
  #                        summaryFunction = twoClassSummary)
  # fit_glm <- train(target ~ preictal_count + interictal_count,
  #                  data=training,
  #                  method="glm", family="binomial",
  #                  metric = "ROC",
  #                  preProc = c("center", "scale"),
  #                  trControl = cctrl1)
  fit_glm <- train(target ~ preictal_count + interictal_count,
                   data=training,
                   method="glm", family="binomial", trControl = trainControl(allowParallel=T, method="cv", number=4))
  
  
  prediction_glm <- predict(fit_glm, testing)
  print("Predictions GLM")
  confuse_glm <- confusionMatrix(prediction_glm, testing$target, positive='preictal')
  print(confuse_glm)
  
  
  # # quick plot of 2nd classifier features
  # testing$prediction <- prediction_glm
  # p <- ggplot(data=testing, aes(x=preictal_count, y=interictal_count, color=prediction)) +
  #   geom_point(position = position_jitter(w = 0.3, h = 0.3))

  
  
  # print("Train on ALL the Data")
  # fit_glm <- train(target ~ preictal_count + interictal_count,
  #                  data=trainset.2,
  #                  method="glm", family="binomial",
  #                  metric = "ROC", 
  #                  preProc = c("center", "scale"),
  #                  trControl = cctrl1)
  # fit_glm <- fit_glm$finalModel
  
  
  trainset.2$prediction <- predict(fit_glm, trainset.2)
  p <- ggplot(data=trainset.2, aes(x=preictal_count, y=interictal_count, color=prediction)) +
    geom_point(position = position_jitter(w = 0.3, h = 0.3)) +
    labs(title="GLM class training labels")
  print(p)
  dev.copy(png, width = 960, height = 960, units = "px", sprintf("train_%s_window_%s_preds.png", patient_num, window_size))
  dev.off()
  
  
  print(p)
  
  print("final model stats")
  print(summary(fit_glm))
  

  # TODO: train a model on all the training data...  
  # print(sprintf("Saving RDS: %s", save_model_filename))
  # saveRDS(fit_glm, save_model_filename)
  
  ################################
  ## TESTING
  patient_num <- 1
  rf_model_filename <- sprintf("../data/ensemble/train_%s_window_%s_secs_correlation_and_fft.random_forest.rds", 
                               patient_num,
                               window_size)
  
  set.seed(1234)
  
  print("Loading models")
  fit_rf <- readRDS(rf_model_filename)

  testing_filename = sprintf("../data/features/corr_fft_basicstats.20161202/test_%s_new_window_%s_secs_correlation_and_fft.testing.txt", 
                             patient_num,
                             window_size)
  print(sprintf("loading TEST data: %s", testing_filename)) 
  
  testset <- load_window_features(output_filename=testing_filename)
  testset <- testset[rowSums(is.na(testset)) == 0,]
  
  # store metacols for the 2nd classifier
  testset.2 <- subset(testset, select=c(id, window_id, segnum, n_dropout_rows))
  # remove metadata cols
  testset <- subset(testset, select=-c(id, window_id, segnum, n_dropout_rows, target))
  
  # Performing predictions with 1st classifier
  print("Performing predictions with 1st classifier (RF)")
  testset.2$prediction <- predict(fit_rf, testset)
  print(table(testset.2$prediction))
  
  # Generate features for 2nd classifier using predictions from the 1st classifier
  print("Generate features for 2nd classifier")
  testset.2 <- summarize(group_by(testset.2, id), preictal_count=sum(prediction=='preictal'), 
                         interictal_count=sum(prediction=='interictal'),
                         total=n())
  head(testset.2)
  print(sprintf("testset.2 nrow: %s", nrow(testset.2)))
  
  print("Performing predictions with 2nd classifier : (GLM)")
  testset.2$prediction <- predict(fit_glm, testset.2)
  
  # transform factor label to (0,1)
  testset.2 <- mutate(testset.2, target=ifelse(prediction=='preictal', 1, 0))
  
  print("verify target and prediction counts are the same...")
  print(table(testset.2$target))
  print(table(testset.2$prediction))
  
  test_2_filename = sprintf("../data/features/test_%s_new_window_%s_secs_correlation_and_fft.mod_2.txt", 
                            patient_num,
                            window_size)
  
  # quick plot of 2nd classifier features
  p <- ggplot(data=testset.2, aes(x=preictal_count, y=interictal_count, color=prediction)) +
    geom_point(position = position_jitter(w = 0.3, h = 0.3))
  labs(title="GLM class submission labels")

  print(p)
  dev.copy(png, width = 960, height = 960, units = "px", sprintf("test_%s_pred_window_%s.png", patient_num, window_size))
  dev.off()
  
  
  # write features with IDs
  print(sprintf("Wrote to : %s", test_2_filename))
  write.csv(testset.2, test_2_filename, row.names=F)
  head(testset.2)
  
  # create submission
  submission_df <- subset(testset.2, select=c(id, target))
  submission_filename = sprintf("test_%s_new_window_%s_secs_predictions.csv",
                                patient_num,
                                window_size)
  
  print(sprintf("Writing predictions to csv: %s", submission_filename))
  write.csv(submission_df, submission_filename, row.names = F)
  
}

ensemble.2.plot <- function(window_size = 30) {
  train_1_preds_filename = sprintf("../data/ensemble/train_1_window_%s_secs_correlation_and_fft.predictions.txt", 
                                   window_size)
  # train_2_filename = sprintf("../data/features/train_1_window_%s_secs_correlation_and_fft.train_2.txt", 
  #                            window_size)
  # save_model_filename="../data/models/train_1_correlation_and_fft.glm_2.rds"
  trainset.2 <- read.csv(train_1_preds_filename, stringsAsFactors = F)
  
  # group by 
  trainset.2 <- summarize(group_by(trainset.2, id), 
                          preictal_count=sum(prediction=='preictal'), 
                          interictal_count=sum(prediction=='interictal'),
                          total=n())
  # trainset.2 <- filter(trainset.2, total>10)
  
  trainset.2$target <- sapply(trainset.2$id, get_target_from_id)
  trainset.2$target <- factor(trainset.2$target, levels=c(0,1), labels=c('interictal', 'preictal'))
  
  p <- ggplot(data=trainset.2, aes(x=preictal_count, y=interictal_count, color=target)) +
    geom_point(position = position_jitter(w = 0.3, h = 0.3))
  print(p)
  # trainset.2
}

