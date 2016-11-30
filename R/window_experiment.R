# generate features for various window sizes

print("generating features")
for (window_size in c(1, 5, 10, 30, 60, 120, 300)) {
  output_filename=sprintf("../data/features/train_1_window_%s_secs_correlation_and_fft.testing.txt", window_size)
  print(sprintf("Generating : %s", output_filename))
  df <- process_windows_parallel(cores=4, 
                                 inputdir="../data/train_1/",
                                 output_filename=output_filename,
                                 secs_per_window=window_size)
  
  # gonna reload from disk using load_window_features()
  rm(df)
  
  # build a model and test
  b_quick = F
  
  trainset <- load_window_features(output_filename=output_filename)
  
  # remove rows with NULLs
  trainset <- trainset[rowSums(is.na(trainset)) == 0,]
  
  # convert to 'target' to a factor variable
  trainset$target <- factor(trainset$target, levels=c(0,1),
                            labels=c('interictal', 'preictal'))
  
  
  # generate models with imputed data  
  save_model_filename=sprintf("../data/models/train_1_window_RF_%s_secs_imputed_TRUE_correlation_and_fft.rds", window_size)
  print(sprintf("Model : %s", save_model_filename))
  train_rf_all_features(trainset = trainset, quick=b_quick, save_model_filename=save_model_filename)
  
  # generate model without imputed data
  save_model_filename=sprintf("../data/models/train_1_window_RF_%s_secs_imputed_FALSE_correlation_and_fft.rds", window_size)
  print(sprintf("Model : %s", save_model_filename))
  
  # remove windows with dropout rows
  trainset <- filter(trainset, n_dropout_rows==0)
  train_rf_all_features(trainset = trainset, quick=b_quick, save_model_filename=save_model_filename)
}



