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

