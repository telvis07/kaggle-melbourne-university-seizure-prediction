source("correlation_features.R")
trainset <- load_window_features()

# remove rows with NULLs
trainset <- trainset[rowSums(is.na(trainset)) == 0,]
train_rf_all_features(trainset = trainset)