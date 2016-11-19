source("correlation_features.R")
trainset <- load_window_features()

# remove rows with NULLs
train_rf_all_features(trainset = trainset)