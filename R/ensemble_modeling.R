# generate features for various window sizes
source("correlation_features.R")
source("glm_grid_search.R")
source("utils.R")
library(DMwR)
library(dplyr)
library(ggplot2)

predict_kaggle_rf_glm <- function(window_size=30, patient_num=1, quick=T){
  
  # load the RF model
  rf_model_filename = sprintf("../data/models/train_1_window_%s_secs_correlation_and_fft.quick_%s.random_forest.rds", window_size, quick)
  print(sprintf("Loading RF model: %s", rf_model_filename))
  fit_rf <- readRDS(rf_model_filename)
  
  # lod the GLM model
  glm_model_filename <- sprintf("../data/models/train_1_window_%s_secs_correlation_and_fft.quick_%s.glm_ensemble.rds", window_size, quick)
  print(sprintf("Loading GLM model: %s", glm_model_filename))
  fit_glm <- readRDS(glm_model_filename)
  
  # load the test data
  testing_filename = sprintf("../data/features/corr_fft_basicstats.20161202/test_%s_new_window_%s_secs_correlation_and_fft.testing.txt", 
                             patient_num,
                             window_size)
  print(sprintf("loading TEST data: %s", testing_filename))
  
  # load the test data
  testset <- load_window_features(output_filename=testing_filename)
  
  # get all the filenames
  all_test_filenames = unique(testset$id)
  
  # remove rows with all NA for features
  testset <- testset[rowSums(is.na(testset)) == 0,]
  testset$target <- factor(testset$target, levels=c('interictal', 'preictal'))
  
  # create meta trainset
  testset_meta <- subset(testset, select=c(id))
  testset_meta <- mutate(testset_meta, M1=NA)
  
  # predict on 30 second windows using the RF model
  testset_meta$M1 <- predict(fit_rf, testset)
  
  print("Window Predictions using RF")
  print(table(testset_meta$M1))
  # print(table(testset_meta$target))
  
  # Now group testset predictions by id (filename)
  # returns id, preictal_count, interictal_count
  testset_meta.byid <- group_meta_trainset_by_id(testset_meta)
  print(sprintf("testset.byid nrow: %s", nrow(testset_meta.byid)))
  
  # now predict using GLM(preictal_count, interictal_count)
  testset_meta.byid$target <- predict(fit_glm, testset_meta.byid)
  
  # now let's get a quick plot of our predictions
  p <- ggplot(data=testset_meta.byid, aes(x=preictal_count, y=interictal_count, color=target)) +
    geom_point(position = position_jitter(w = 0.3, h = 0.3))
  labs(title="GLM class submission labels")
  save_plot_filename = sprintf("train_%s_window_%s_quick_%s_preds.png", patient_num, window_size, quick)
  save_plot(p, save_plot_filename)
  
  table(testset_meta.byid$target)
  print(sprintf("number of filenames in raw data: %s",length(all_test_filenames)))
  print(sprintf("number of filenames in the predictions: %s", dim(testset_meta.byid)[1]))
  
  # unfortunately, some segment files had all NA features for all windows
  # so just guess 'interictal' for missing files
  submission_df <- merge(data.frame(id=all_test_filenames), 
                         subset(testset_meta.byid, select=c(id, target)), 
                         by="id", all.x=T)
  na_targets <- is.na(submission_df$target)
  submission_df[na_targets,]$target = "interictal"
  
  # now change 'interical', 'preictal' factor strings to '0', '1'
  submission_df$target = factor(submission_df$target, levels=c('interictal', 'preictal'), labels=c(0, 1))
  
  # finally, let's create a file with submissions
  submission_filename = sprintf("test_%s_new_window_%s_secs_quick_%s_predictions.csv",
                                patient_num,
                                window_size,
                                quick)
  
  print(sprintf("Writing predictions to csv: %s", submission_filename))
  write.csv(submission_df, submission_filename, row.names = F)
  
  submission_df
}

combine_submission_files <- function(dirname="../data/submissions"){
  filenames = c(
    "../data/submissions/test_1_new_window_30_secs_quick_FALSE_predictions.csv",
    "../data/submissions/test_2_new_window_30_secs_quick_FALSE_predictions.csv",
    "../data/submissions/test_3_new_window_30_secs_quick_FALSE_predictions.csv"
  )
  
  submission_df <- data.frame()
  
  for (filename in filenames){
    df <- read.csv(filename, header = T, stringsAsFactors = F)
    submission_df <- rbind(submission_df, df)
  }
  
  # finally, let's create a file with submissions
  submission_filename = "test_ALL_new_window_30_secs_predictions.csv"
  
  names(submission_df) <- c("File", "Class")
  
  
  print(sprintf("Writing predictions to csv: %s", submission_filename))
  write.csv(submission_df, submission_filename, row.names = F)
  
  
  submission_df
  
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

