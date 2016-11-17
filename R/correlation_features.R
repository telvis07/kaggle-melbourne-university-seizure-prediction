#matlab extractor
rm(list=ls())
setwd("/Users/telvis/work/kaggle-melbourne-university-seizure-prediction/R")

library(ggplot2)
library(entropy)


feat_corr_eigen <- function(window_all_channels, verbose=T) {
  # Calculate correlation matrix and its eigenvalues (b/w channels)
  
  m_channel_corrs <- abs(cor(window_all_channels, method="pearson"))
  if (verbose) {
    print(m_channel_corrs)
  }
  m_channel_corrs[!is.finite(m_channel_corrs)] <- 0
  diag(m_channel_corrs) <- 0
  
  if (verbose){
    print(sprintf("Shape correlation matrix: %s", as.character(dim(m_channel_corrs)), 
                  collapse = ","))
  }
  
  l_eigen = eigen(m_channel_corrs, symmetric = T, only.values = T)
  if (verbose){
    print(l_eigen$values)
    print(sprintf("Eigenvalues: %s, %s", 
                  paste(length(l_eigen$values), collapse = ","), 
                  class(l_eigen$values)))
  }
  l_eigen$values
}

feat_fft_sums <- function(window_all_channels, n_samples_per_window=4000, verbose=T) {
  # transpose the eeg data : N x 16 -> 16 x N
  # calc fft
  # 
  m_chan_t <- t(window_all_channels)
  # if (dim(m_chan)[2] < n_samples_per_window){
  #   n_num_cols <- dim(m_chan)[2]
  #   n_num_rows <- dim(m_chan)[1]
  #   m_chan_t <- cbind(m_chan_t, matrix(0, n_num_rows, n_samples_per_window-n_num_cols))
  # }
  
  print(sprintf("dim(m_chan_t) : %s", paste(dim(m_chan_t), " ")))
  
  m_fft <- fft(m_chan_t)
  v_freq_magnitude_sums <- colSums(abs(m_fft))
  v_freq_magnitude_sums
  
  # library(entropy)
  # m_fft <- fft(t(m))
  # m_mags <- t(abs(m_fft))
  # 
  # > entropy.empirical(m_mags, unit="log2")
  # [1] 1.961536
  # > -sum(m_mags * log2(m_mags))
  # [1] NaN
  
}

feat_fft_means <- function(window_all_channels, n_samples_per_window=4000, verbose=T) {
  # transpose the eeg data : N x 16 -> 16 x N
  # calc fft
  # 
  m_chan_t <- t(window_all_channels)
  # if (dim(m_chan)[2] < n_samples_per_window){
  #   n_num_cols <- dim(m_chan)[2]
  #   n_num_rows <- dim(m_chan)[1]
  #   m_chan_t <- cbind(m_chan_t, matrix(0, n_num_rows, n_samples_per_window-n_num_cols))
  # }
  
  print(sprintf("dim(m_chan_t) : %s", paste(dim(m_chan_t), " ")))
  
  m_fft <- fft(m_chan_t)
  v_freq_magnitude_means <- apply((abs(m_fft)), 2, mean)
  v_freq_magnitude_means
  
  # library(entropy)
  # m_fft <- fft(t(m))
  # m_mags <- t(abs(m_fft))
  # 
  # > entropy.empirical(m_mags, unit="log2")
  # [1] 1.961536
  # > -sum(m_mags * log2(m_mags))
  # [1] NaN
  
}


feat_fft_mag_entropy <- function(window_all_channels, n_samples_per_window=4000, verbose=T) {
  # transpose the eeg data : N x 16 -> 16 x N
  m_chan_t <- t(window_all_channels)

  # sanity check the dimensions   
  print(sprintf("dim(m_chan_t) : %s", paste(dim(m_chan_t), " ")))
  
  # calc fft for channel data
  m_fft <- fft(m_chan_t)
  
  # calculate the magnitude fft values
  m_mags <- abs(m_fft)
  
  # return the 
  entropy.empirical(m_mags, unit="log2")
}

feat_all_channels_entropy <- function(window_all_channels, n_samples_per_window=4000, verbose=T) {
  # return the entropy of all channels
  entropy.empirical(as.matrix(window_all_channels), unit="log2")
}


make_feature_names <- function() {
  # generate the column names for the features
  v_names <- c("id", "target", "window_id", "segnum",
               paste("chan_eigenval_", 1:16, sep=""), 
               "fft_mag_entropy",
               "all_chan_entropy")
  
  v_names
}


process_windows_corr_fft <- function(data, id, target, n_seg_num, secs_per_window=10, verbose=T){
  # (1) eeg data
  # (2) nSamplesSegment: total number of time samples (number of rows in the data field).
  # (3) iEEGsamplingRate: data sampling rate, i.e. the number of data samples 
  #   representing 1 second of EEG data. 
  # (4) channelIndices: an array of the electrode indexes corresponding to 
  #     the columns in the data field.
  
  df_eeg = as.data.frame(data[[1]][[1]])
  
  # fs - 400Hz (samples/second) 
  iEEGsamplingRate <- data$dataStruct[[2]][1]
  # n - 240000
  nSamplesSegment <- data$dataStruct[[3]][1]
  # 16 channels
  channelIndices <- data$dataStruct[[4]]
  n_num_channels <- length(channelIndices)
  
  # samples / (samples/second) = seconds
  # seconds / (seconds/window) = number of windows
  n_num_windows <- (nSamplesSegment/iEEGsamplingRate) / secs_per_window
  
  # nsamples for a 10 second window
  # (samples/second) * (seconds/window) = samples/window 
  n_samples_per_window <- iEEGsamplingRate * secs_per_window
  n_window_id <- 0
  
  v_feature_names <- make_feature_names()
  trainset = matrix(,n_num_windows,length(v_feature_names))
  trainset = as.data.frame(trainset)
  
  # window_id, min, max, mean, median, 
  for (n_offset in seq(1, nSamplesSegment, n_samples_per_window)) {
    n_window_id <- n_window_id + 1
    
    # beginning and end of window
    n_offset_end <- n_offset+n_samples_per_window - 1
    
    # check for data drop-out, where values across all channels are 0
    window_all_channels <- df_eeg[n_offset:n_offset_end, ]
    if (! all(dim(window_all_channels) == c(n_samples_per_window, n_num_channels))){
      stop(sprintf("[process_windows] %s: Bad Dimensions: %s", id, paste(dim(window_all_channels), collapse="x ")))                            
    }
    
    li_absum_by_row <- as.numeric(apply(window_all_channels, 
                                        1,
                                        function (x) sum(abs(x))))
    n_dropout_rows <- sum(li_absum_by_row == 0)
    
    if (n_dropout_rows > 0) {
      print(sprintf("[process_windows] %s: n_dropout_rows = %s", id, n_dropout_rows))
    }
    
    v_dropout_rows <- (li_absum_by_row == FALSE)
    
    if (n_dropout_rows > 3500) {
      print(sprintf("[process_windows] - too many dropout rows in window: %s", n_window_id))
      # insert row with all features as NA
      trainset[n_window_id,] = c(id,
                                 target,
                                 n_window_id,
                                 n_seg_num,
                                 rep(NA, length(v_feature_names) - 4))
      next
    } 
    
    if( n_dropout_rows > 0) {
      ## Impute the drop out rows with the row-mean
      print(sprintf("Imputing file %s window %s", id, n_window_id))
      
      # set the dropout rows to NA
      window_all_channels[v_dropout_rows, ] <- NA
      
      # R fills matrices by columns, so transpose the matrix
      # to fill missing values with mean for value for the channel
      m_t <- t(window_all_channels)
      v_impute_vals <- apply(m_t, 1, mean, na.rm=T)
      m_t[,v_dropout_rows] <- v_impute_vals
      
      # now transpose it again, return to the normal shape
      window_all_channels <- t(m_t)
      # cleanup m_t
      rm(m_t)
    }
    
    print(dim(window_all_channels))
    if (dim(window_all_channels)[1] != n_samples_per_window) {
      stop(sprintf("[process_windows] - invalid window length. Expected: %s, Found: %s",
                   n_samples_per_window,
                   dim(window_all_channels)[1]))
    }
    
    # correlation and eigenvalues
    v_eigen_values <- feat_corr_eigen(window_all_channels = window_all_channels)
    n_fft_mag_entropy <- feat_fft_mag_entropy(window_all_channels = window_all_channels)
    n_all_chan_entropy <- feat_all_channels_entropy(window_all_channels = window_all_channels)
    
    trainset[n_window_id,] = c(id,
                               target,
                               n_window_id,
                               n_seg_num,
                               v_eigen_values,
                               n_fft_mag_entropy,
                               n_all_chan_entropy)
    
  }
  colnames(trainset) <- v_feature_names
  trainset
  
}

process_file_windows_single <- function(filename) {
  # trainset <- process_file_windows_single("../data/train_1/1_999_0.mat")
  
  # parse the 'preictal, interical' label
  data = readMat(filename)
  s_base_filename <- basename(filename)
  v_filename_parts <- strsplit(s_base_filename, "_")[[1]]
  s_target <- v_filename_parts[3]
  s_target <- gsub(".mat", "", s_target)
  n_seg_num <- as.numeric(v_filename_parts[2])
  
  # calculate window features
  # trainset <- process_windows(data = data, id=s_base_filename, target=s_target, n_seg_num=n_seg_num)
  trainset <- process_windows_corr_fft(data = data, id=s_base_filename, target=s_target, n_seg_num=n_seg_num)
  trainset
}

process_windows_parallel <- function(inputdir="../data/train_1.small/",
                                     output_filename="../data/features/train_1_correlation_and_fft.txt",
                                     cores=4) {
  # df <- process_windows_parallel(cores=8, inputdir = "../data/train_1/")
  
  filenames <- list.files(inputdir, pattern="*.mat", full.names=TRUE)
  runtime <- system.time({
    trainset <- mclapply(filenames, process_file_windows_single, mc.cores = cores)
  })[3]
  print(sprintf("runtime: %s", runtime))
  print(sprintf("length: %s", length(trainset)))
  df <- do.call("rbind", trainset)
  
  v_names <- make_feature_names()
  
  print(sprintf("Dimensions: %s", paste(dim(df), collapse="x ")))
  colnames(df) <- v_names
  
  # sort by target then segment number
  ord = order(df$target, df$segnum, decreasing = F)
  df <- df[ord,]
  
  # write to file
  print(sprintf("Writing features to : %s", output_filename))
  write.table(df, output_filename,
              quote=F,row.names=F,col.names=T,sep="\t")
  df
}


#############################################
# modeling
############################################
load_window_features <- function(output_filename="../data/features/train_1_correlation_and_fft.txt") {
  df <- read.table(output_filename, header=T, stringsAsFactors = F)
  df$segnum <- as.numeric(df$segnum)
  df$target <- as.numeric(df$target)
  
  # sort by target then segment number
  ord = order(df$target, df$segnum, decreasing = F)
  df <- df[ord,]
  df
}

train_rf_all_features <- function(trainset, seed=1234, save_model_filename="../data/models/train_1_correlation_and_fft.rds") {
  set.seed(seed)
  
  inTrain = createDataPartition(trainset$target, p = 3/4)[[1]]
  training = trainset[ inTrain,]
  testing = trainset[-inTrain,] 
  

  modFit <- train(target ~ ., 
                  data=training, 
                  method="rf", ntree=100,
                  trControl = trainControl(method="cv"), number=10)
  
  # testing
  prediction_rf <- predict(modFit, testing)
  print("Predictions RF : testing")
  print(confusionMatrix(prediction_rf, testing$target))
  
  # training
  prediction_rf <- predict(modFit, training)
  print("Predictions RF:training")
  print(confusionMatrix(prediction_rf, training$classe))
  
  print("Saving RDS")
  saveRDS(modFit, save_model_filename)
  
  # return model with all features
  modFit
}


#############################################
# Exploration
#############################################


get_first_windows_from_file <- function(filename="../data/train_1/1_999_0.mat"){
  # Get the first window from a file.
  
  # plot eignenvectors
  
  secs_per_window=10
  # filename = "../data/train_1/1_999_0.mat"
  
  # parse the 'preictal, interical' label
  data = readMat(filename)
  s_base_filename <- basename(filename)
  v_filename_parts <- strsplit(s_base_filename, "_")[[1]]
  s_target <- v_filename_parts[3]
  s_target <- gsub(".mat", "", s_target)
  n_seg_num <- as.numeric(v_filename_parts[2])
  
  # 
  # trainset <- compute_corr_for_windows(data = data, id=s_base_filename, target=s_target, n_seg_num=n_seg_num)
  df_eeg = as.data.frame(data[[1]][[1]])
  
  # fs - 400Hz (samples/second) 
  iEEGsamplingRate <- data$dataStruct[[2]][1]
  # n - 240000
  nSamplesSegment <- data$dataStruct[[3]][1]
  # 16 channels
  channelIndices <- data$dataStruct[[4]]
  n_num_channels <- length(channelIndices)
  # samples / (samples/second) = seconds
  # seconds / (seconds/window) = number of windows
  n_num_windows <- (nSamplesSegment/iEEGsamplingRate) / secs_per_window
  
  # nsamples for a 10 second window
  # (samples/second) * (seconds/window) = samples/window 
  n_samples_per_window <- iEEGsamplingRate * secs_per_window
  
  
  n_offset <- 1
  n_offset_end <- n_offset+n_samples_per_window - 1
  window_all_channels <- df_eeg[n_offset:n_offset_end, ]
  li_absum_by_row <- as.numeric(apply(window_all_channels, 
                                      1,
                                      function (x) sum(abs(x))))
  n_dropout_rows <- sum(li_absum_by_row == 0)
  
  print(sprintf("Num dropout rows : %s", n_dropout_rows))
  
  window_all_channels
}

plot_feat_corr_eigen <- function() {

  window_all_channels_0 <- get_windows_from_file(filename = "../data/train_1/1_999_0.mat")
  
  v_eigen_values <- feat_corr_eigen(window_all_channels = window_all_channels_0)
  
  plt1 <- ggplot(data.frame(eigenvalues=v_eigen_values, 
                           num=seq(length(v_eigen_values))), 
                aes(x=num, y=eigenvalues)) +
    geom_bar(stat="identity") +
    labs(title="cor eigenvalues (interical)")
  
  
  window_all_channels_1 <- get_windows_from_file(filename = "../data/train_1/1_100_1.mat")
  
  v_eigen_values <- feat_corr_eigen(window_all_channels = window_all_channels_1)
  
  plt2 <- ggplot(data.frame(eigenvalues=v_eigen_values, 
                           num=seq(length(v_eigen_values))), 
                aes(x=num, y=eigenvalues)) +
    geom_bar(stat="identity") +
    labs(title="cor eignenvalues (preictal)")
  
  grid.arrange(plt1, plt2)
}

plot_feat_fft_means <- function() {
  
  window_all_channels_0 <- get_windows_from_file(filename = "../data/train_1/1_999_0.mat")
  
  v_feats.0 <- feat_fft_means(window_all_channels = window_all_channels_0)
  
  plt1 <- ggplot(data.frame(feature=v_feats.0, 
                            num=seq(length(v_feats.0))), 
                 aes(x=num, y=v_feats.0)) +
    geom_bar(stat="identity") +
    labs(title="fft means (interical)")
  
  
  window_all_channels_1 <- get_windows_from_file(filename = "../data/train_1/1_100_1.mat")
  
  v_feats.1 <- feat_fft_means(window_all_channels = window_all_channels_1)
  
  plt2 <- ggplot(data.frame(feature=v_feats.1, 
                            num=seq(length(v_feats.1))), 
                 aes(x=num, y=v_feats.1)) +
    geom_bar(stat="identity") +
    labs(title="fft means (preictal)")
  
  grid.arrange(plt1, plt2)
}


plot_feat_fft_entropy <- function() {
  
  trainset <- calc_corr_for_windows()
  
  window_all_channels_0 <- get_windows_from_file(filename = "../data/train_1/1_999_0.mat")
  
  v_feats.0 <- feat_fft_entropy(window_all_channels = window_all_channels_0)
  print(v_feats.0)
  
  # plt1 <- ggplot(data.frame(feature=v_feats.0, 
  #                           num=seq(length(v_feats.0))), 
  #                aes(x=num, y=v_feats.0)) +
  #   geom_bar(stat="identity") +
  #   labs(title="fft means (interical)")
  
  
  window_all_channels_1 <- get_windows_from_file(filename = "../data/train_1/1_100_1.mat")
  
  v_feats.1 <- feat_fft_entropy(window_all_channels = window_all_channels_1)
  print(v_feats.1)
  
  # plt2 <- ggplot(data.frame(feature=v_feats.1, 
  #                           num=seq(length(v_feats.1))), 
  #                aes(x=num, y=v_feats.1)) +
  #   geom_bar(stat="identity") +
  #   labs(title="fft means (preictal)")
  # 
  # grid.arrange(plt1, plt2)
}

