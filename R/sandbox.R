run_single <- function (secs_per_window=10){
  secs_per_window=10
  filename = "../data/train_1/1_999_0.mat"
  
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
  dim(window_all_channels)
  
  v_lvl = c(0.1, 4, 8, 12, 30, 70, 180)
  v_lseg = round(nSamplesSegment/iEEGsamplingRate*v_lvl)
  print(v_lseg)
}