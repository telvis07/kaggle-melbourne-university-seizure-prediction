#matlab extractor

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

library(entropy)

feat_fft_entropy <- function(window_all_channels, n_samples_per_window=4000, verbose=T) {
  # transpose the eeg data : N x 16 -> 16 x N
  m_chan_t <- t(window_all_channels)

  # sanity check the dimensions   
  print(sprintf("dim(m_chan_t) : %s", paste(dim(m_chan_t), " ")))
  
  m_fft <- fft(m_chan_t)
  m_mags <- abs(m_fft)

  # m_fft <- fft(t(m))
  # m_mags <- t(abs(m_fft))
  # 
  # > entropy.empirical(m_mags, unit="log2")
  # [1] 1.961536
  # > -sum(m_mags * log2(m_mags))
  # [1] NaN
  entropy.empirical(m_mags, unit="log2")
}


calc_corr_for_windows <- function(data, id, target, n_seg_num, secs_per_window=10, verbose=T){
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
      # TODO:
      # trainset[n_window_id,] = c(id,
      #                            target,
      #                            n_window_id,
      #                            n_seg_num,
      #                            rep(NA, n_num_channels),
      #                            rep(NA, n_num_channels),
      #                            rep(NA, n_num_channels),
      #                            rep(NA, n_num_channels),
      #                            rep(NA, n_num_channels),
      #                            rep(NA, n_num_channels))
      next
      
    } 
    
    if( n_dropout_rows > 0) {
      # remove dropout rows
      # window_all_channels <- window_all_channels[! v_dropout_rows,]
      windows_all_channels[v_dropout_rows, ] <- NA
      
      # R fills matrices by columns, so transpose the matrix
      # to fill missing values with mean for value for the channel
      m_t <- t(windows_all_channels)
      v_impute_vals <- apply(m_t, 1, mean, na.rm=T)
      m_t[,v_dropout_rows] <- v_impute_vals
      
      # now transpose it again, return to the normal shape
      windows_all_channels <- t(m_t)
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
    
  }
  
}

library(kernlab); data(spam); library(ggplot2)

spam_correlated_predictors <- function() {
  inTrain <- createDataPartition(y=spam$type, p=0.75, list=FALSE)
  training <- spam[inTrain,]
  testing <- spam[-inTrain,]
  # col 58 is the 'type' col
  M <- abs(cor(training[,-58], method="pearson"))
  diag(M) <- 0
  print(which(M > 0.8, arr.ind = T))
  
  # output shows
  # row col
  # num415  34  32
  # direct  40  32
  # num857  32  34
  # direct  40  34
  # num857  32  40
  # num415  34  40
  
  # > names(spam[c(32,34)])
  # [1] "num857" "num415"
  
  # show column names for cols 34 and 42
  print(names(spam[c(32,34)]))
  qplot(spam[,34], spam[,32])
}
