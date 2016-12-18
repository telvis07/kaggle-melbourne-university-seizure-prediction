# generate features for various window sizes
source("correlation_features.R")
source("sample_data.R")
library(DMwR)

run_buncha_models <- function (window_size = 30) {
  output_filename = sprintf("../data/features/corr_fft_basicstats.20161202/train_2_window_%s_secs_correlation_and_fft.testing.txt", 
                            window_size)
  trainset <- load_window_features(output_filename=output_filename)
  trainset <- trainset[rowSums(is.na(trainset)) == 0,]
  trainset <- subset(trainset, select=-c(id, window_id, segnum, n_dropout_rows))
  trainset$target <- factor(trainset$target, levels=c(0,1), labels=c('interictal', 'preictal'))
  print(names(trainset))
  
  # nosampling
  print("nosampling info")
  print(table(trainset$target))
  buncha_models.scale(trainset = trainset, 
                      save_stats_filename="patient_2_buncha_model_stats_quick_FALSE_nosampling_SCALE.csv")
  
  # downsample
  down_train <- downSample(x=subset(trainset, select=-c(target)), y=trainset$target, yname="target")
  print("downsample info")
  print(table(down_train$target))
  buncha_models.scale(trainset = down_train, 
                      save_stats_filename="patient_2_buncha_model_stats_quick_FALSE_downsample_CARET_SCALE.csv")
  
  # upsample
  up_train <- upSample(x=subset(trainset, select=-c(target)), y=trainset$target, yname="target")
  print("upsample info")
  print(table(up_train$target))
  buncha_models.scale(trainset = up_train, 
                      save_stats_filename="patient_2_buncha_model_stats_quick_FALSE_upsample_CARET_SCALE.csv")
  
  # smote
  smote_train <- SMOTE(target ~ ., data=trainset)
  print("smote info")
  print(table(smote_train$target))
  buncha_models.scale(trainset = smote_train, 
                      save_stats_filename="patient_2_buncha_model_stats_quick_FALSE_smote_CARET_SCALE.csv")

}

run_buncha_models.manual_downsample <- function (window_size = 30) {
  output_filename = sprintf("../data/features/corr_fft_basicstats.20161202/train_1_window_%s_secs_correlation_and_fft.testing.txt", 
                            window_size)
  trainset <- load_window_features(output_filename=output_filename)
  trainset <- trainset[rowSums(is.na(trainset)) == 0,]
  trainset <- subset(trainset, select=-c(id, window_id, segnum, n_dropout_rows))
  trainset$target <- factor(trainset$target, levels=c(0,1), labels=c('interictal', 'preictal'))
  buncha_models.scale(trainset = trainset, 
                downsample_negclass=10000,
                save_stats_filename="buncha_model_stats_quick_FALSE_downsample_10000_SCALE.csv")
}

buncha_models <- function(trainset, seed=1234, quick=FALSE, downsample_negclass=0, save_stats_filename="buncha_model_stats.csv") {
  set.seed(seed)
  
  if (quick) {
    # downsample, so we can run quicker
    trainset <- sample_data(trainset, n_neg_samples=1000, n_pos_samples=500)
  } else if (downsample_negclass){
    trainset <- sample_data(trainset, 
                            n_neg_samples=downsample_negclass, 
                            n_pos_samples=0)
  }
  
  # remove columns that we don't train on
  # trainset <- subset(trainset, select=-c(id, window_id, segnum, fft_phase_entropy, n_dropout_rows))
  
  print(sprintf("after removing columns: %s, %s", dim(trainset)[1],
                dim(trainset)[2]))
  
  
  inTrain = createDataPartition(trainset$target, p = 3/4)[[1]]
  training = trainset[ inTrain,]
  testing = trainset[-inTrain,]
  
  print(sprintf("training dim: %s", dim(training)[1] ))
  print(sprintf("testing dim: %s", dim(testing)[1] ))
  
  print (table(training$target))
  print (table(testing$target))
  
  
  print("Training RF")
  fit_rf <- train(target ~ ., data=training, method="rf", trControl = trainControl(allowParallel=T, method="cv", number=4))
  prediction_rf <- predict(fit_rf, testing)
  print("Predictions RF")
  confuse_rf <- confusionMatrix(prediction_rf, testing$target, positive='preictal')
  print(confuse_rf)

  print("Training GBM")
  fit_gbm <- train(target ~ ., data=training, method="gbm", verbose=FALSE, trControl = trainControl(allowParallel=T, method="cv", number=4))
  prediction_gbm <- predict(fit_gbm, testing)
  print("Predictions GBM")
  confuse_gbm <- confusionMatrix(prediction_gbm, testing$target, positive='preictal')
  print(confuse_gbm)
  
  print("Training LDA")
  fit_lda <- train(target ~ ., data=training, method="lda", trControl = trainControl(allowParallel=T, method="cv", number=4))
  prediction_lda <- predict(fit_lda, testing)
  print("Predictions LDA")
  confuse_lda <- confusionMatrix(prediction_lda, testing$target, positive='preictal')
  print(confuse_lda)
  
  print("Training Logistic Regression")
  fit_glm <- train(target ~ .,
                  data=training,
                  method="glm", family="binomial", trControl = trainControl(allowParallel=T, method="cv", number=4))
  prediction_glm <- predict(fit_glm, testing)
  print("Predictions GLM")
  confuse_glm <- confusionMatrix(prediction_glm, testing$target, positive='preictal')
  print(confuse_glm)
  
  # stacked using rf
  # fit a model that combines predictors
  predDF <- data.frame(prediction_rf, prediction_gbm, prediction_lda, target=testing$target)
  
  print("Training ensemble : RF, GBM, LDA")
  combModFit <- train(target ~ ., method="rf", data=predDF, trControl = trainControl(allowParallel=T, method="cv", number=4))
  combPred <- predict(combModFit, predDF)
  print("Predictions Combined: RF, GBM, LDA")
  confuse_comb <- confusionMatrix(combPred, testing$target, positive='preictal')
  print
  
  # Stacked Accuracy: 0.80 is better than random forests and lda and the same as boosting.
  
  #####
  df_rf <- as.data.frame(t(confuse_rf$byClass))
  df_rf$model <- "rf"
  names(df_rf) <- make.names(names(df_rf))
  
  df_gbm <- as.data.frame(t(confuse_gbm$byClass))
  df_gbm$model <- "gbm"
  names(df_gbm) <- make.names(names(df_gbm))
  
  df_lda <- as.data.frame(t(confuse_lda$byClass))
  df_lda$model <- "lda"
  names(df_lda) <- make.names(names(df_lda))
  
  df_glm <- as.data.frame(t(confuse_glm$byClass))
  df_glm$model <- "glm"
  names(df_glm) <- make.names(names(df_glm))
  
  df_ens <- as.data.frame(t(confuse_comb$byClass))
  df_ens$model <- "ensemble"
  names(df_ens) <- make.names(names(df_ens))
  
  all_stats_df <- rbind(df_rf, df_gbm, df_lda, df_glm, df_ens)
  write.csv(all_stats_df, save_stats_filename, row.names = F)
  
  all_stats_df
}


buncha_models.scale <- function(trainset, seed=1234, quick=FALSE, downsample_negclass=0, save_stats_filename="buncha_model_stats.csv") {
  set.seed(seed)
  
  # if (quick) {
  #   # downsample, so we can run quicker
  #   trainset <- sample_data(trainset, n_neg_samples=1000, n_pos_samples=500)
  # } else if (downsample_negclass){
  #   trainset <- sample_data(trainset, 
  #                           n_neg_samples=downsample_negclass, 
  #                           n_pos_samples=0)
  # } 
  
  # remove columns that we don't train on
  # trainset <- subset(trainset, select=-c(id, window_id, segnum, fft_phase_entropy, n_dropout_rows))
  
  print(sprintf("after removing columns: %s, %s", dim(trainset)[1],
                dim(trainset)[2]))
  
  
  inTrain = createDataPartition(trainset$target, p = 3/4)[[1]]
  training = trainset[ inTrain,]
  testing = trainset[-inTrain,]
  
  print(sprintf("training dim: %s", dim(training)[1] ))
  print(sprintf("testing dim: %s", dim(testing)[1] ))
  
  print (table(training$target))
  print (table(testing$target))
  
  
  print("Training RF")
  fit_rf <- train(target ~ ., data=training, method="rf", trControl = trainControl(allowParallel=T, method="cv", number=4))
  prediction_rf <- predict(fit_rf, testing)
  print("Predictions RF")
  confuse_rf <- confusionMatrix(prediction_rf, testing$target, positive='preictal')
  print(confuse_rf)
  
  print("Training GBM")
  fit_gbm <- train(target ~ ., data=training, method="gbm", preProcess=c("center","scale"), verbose=FALSE, trControl = trainControl(allowParallel=T, method="cv", number=4))
  prediction_gbm <- predict(fit_gbm, testing)
  print("Predictions GBM")
  confuse_gbm <- confusionMatrix(prediction_gbm, testing$target, positive='preictal')
  print(confuse_gbm)
  
  print("Training LDA")
  fit_lda <- train(target ~ ., data=training, method="lda", preProcess=c("center","scale"), trControl = trainControl(allowParallel=T, method="cv", number=4))
  prediction_lda <- predict(fit_lda, testing)
  print("Predictions LDA")
  confuse_lda <- confusionMatrix(prediction_lda, testing$target, positive='preictal')
  print(confuse_lda)
  
  print("Training Logistic Regression")
  fit_glm <- train(target ~ .,
                   data=training,
                   preProcess=c("center","scale"),
                   method="glm", family="binomial", trControl = trainControl(allowParallel=T, method="cv", number=4))
  prediction_glm <- predict(fit_glm, testing)
  print("Predictions GLM")
  confuse_glm <- confusionMatrix(prediction_glm, testing$target, positive='preictal')
  print(confuse_glm)
  
  # stacked using rf
  # fit a model that combines predictors
  predDF <- data.frame(prediction_rf, prediction_gbm, prediction_lda, target=testing$target)
  
  print("Training ensemble : RF, GBM, LDA")
  combModFit <- train(target ~ ., method="rf", data=predDF, preProcess=c("center","scale"), trControl = trainControl(allowParallel=T, method="cv", number=4))
  combPred <- predict(combModFit, predDF)
  print("Predictions Combined: RF, GBM, LDA")
  confuse_comb <- confusionMatrix(combPred, testing$target, positive='preictal')
  print
  
  # Stacked Accuracy: 0.80 is better than random forests and lda and the same as boosting.
  
  #####
  df_rf <- as.data.frame(t(confuse_rf$byClass))
  df_rf$model <- "rf"
  names(df_rf) <- make.names(names(df_rf))
  
  df_gbm <- as.data.frame(t(confuse_gbm$byClass))
  df_gbm$model <- "gbm"
  names(df_gbm) <- make.names(names(df_gbm))
  
  df_lda <- as.data.frame(t(confuse_lda$byClass))
  df_lda$model <- "lda"
  names(df_lda) <- make.names(names(df_lda))
  
  df_glm <- as.data.frame(t(confuse_glm$byClass))
  df_glm$model <- "glm"
  names(df_glm) <- make.names(names(df_glm))
  
  df_ens <- as.data.frame(t(confuse_comb$byClass))
  df_ens$model <- "ensemble"
  names(df_ens) <- make.names(names(df_ens))
  
  all_stats_df <- rbind(df_rf, df_gbm, df_lda, df_glm, df_ens)
  write.csv(all_stats_df, save_stats_filename, row.names = F)
  
  all_stats_df
}