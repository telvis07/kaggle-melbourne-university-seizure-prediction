


buncha_models <- function(trainset, seed=1234, quick=FALSE, save_model_filename="../data/models/train_1_glm_correlation_and_fft.rds") {
  set.seed(seed)
  
  if (quick) {
    # downsample, so we can run quicker
    trainset <- sample_data(trainset)
    ntree = 10
    n_rf_number = 5
  } else {
    ntree = 100
    n_rf_number = 10
  }
  
  # remove columns that we don't train on
  trainset <- subset(trainset, select=-c(id, window_id, segnum, fft_phase_entropy, n_dropout_rows))
  
  print(sprintf("after removing columns: %s, %s", dim(trainset)[1],
                dim(trainset)[2]))
  
  
  inTrain = createDataPartition(trainset$target, p = 3/4)[[1]]
  training = trainset[ inTrain,]
  testing = trainset[-inTrain,]
  
  print(sprintf("training dim: %s", dim(training)[1] ))
  print(sprintf("testing dim: %s", dim(testing)[1] ))
  
  print (table(training$target))
  print (table(testing$target))
  
  
  # TODO : Fix me
  
  set.seed(62433)
  print("Training RF")
  fit_rf <- train(training$diagnosis ~ ., data=training, method="rf")
  prediction_rf <- predict(fit_rf, testing)
  print("Predictions RF")
  print(confusionMatrix(prediction_rf, testing$diagnosis))
  # Accuracy : 0.7805 
  
  print("Training GBM")
  fit_gbm <- train(training$diagnosis ~ ., data=training, method="gbm", verbose=FALSE)
  prediction_gbm <- predict(fit_gbm, testing)
  print("Predictions GBM")
  print(confusionMatrix(prediction_gbm, testing$diagnosis))
  # Accuracy : 0.8049 
  
  print("Training LDA")
  fit_lda <- train(training$diagnosis ~ ., data=training, method="lda")
  prediction_lda <- predict(fit_lda, testing)
  print("Predictions LDA")
  print(confusionMatrix(prediction_lda, testing$diagnosis))
  # Accuracy : 0.7683
  
  # stacked using rf
  # fit a model that combines predictors
  predDF <- data.frame(prediction_rf, prediction_gbm, prediction_lda, diagnosis=testing$diagnosis)
  
  print("Training ensemble of mod1 and mod2: GAM")
  combModFit <- train(diagnosis ~ ., method="rf", data=predDF)
  combPred <- predict(combModFit, predDF)
  print(confusionMatrix(combPred, testing$diagnosis))
  # Accuracy : 0.8171 
  
  # Stacked Accuracy: 0.80 is better than random forests and lda and the same as boosting.
  
  
  
  #####
  
  
}