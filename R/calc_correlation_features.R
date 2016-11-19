source("correlation_features.R")
source("sample_data.R")
trainset <- load_window_features()

# # remove rows with NULLs
trainset <- trainset[rowSums(is.na(trainset)) == 0,]

# convert to 'target' to a factor variable
trainset$target <- factor(trainset$target, levels=c(0,1), 
                          labels=c('interictal', 'preictal'))

# sample datums
# trainset <- sample_data(trainset)

# remove rows with NULLs
train_rf_all_features(trainset = trainset)