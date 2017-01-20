library(caret)
library(randomForest)
library(rpart)
library(DMwR)
library(gbm)
library(dplyr)
library(ggplot2)
source("correlation_features.R")
source("utils.R")


run_submission_predictions <- function(window_size=30, patient_num=1){
  
  rf_model_filename <- sprintf("../data/ensemble/train_%s_window_%s_secs_correlation_and_fft.random_forest.rds", 
                               patient_num,
                               window_size)
  glm_model_filename <- sprintf("../data/ensemble/train_%s_window_%s_secs_correlation_and_fft.glm_2.rds", 
                                patient_num,
                                window_size)
  
  set.seed(1234)
  
  print("Loading models")
  fit_rf <- readRDS(rf_model_filename)
  fit_glm <- readRDS(glm_model_filename)
  
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
  print(p)
  
  
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

# [telvis@yoyoma data]$ ls test_1_new | wc -l
# 216
# [telvis@yoyoma data]$ ls test_2_new | wc -l
# 1002
# [telvis@yoyoma data]$ ls test_3_new | wc -l
# 690
