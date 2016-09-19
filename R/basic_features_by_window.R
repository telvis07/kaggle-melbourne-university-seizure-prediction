# extract features for each window in each file.
library("parallel")

process_windows <- function(data, id, target, secs_per_window=10) {
  # (2) nSamplesSegment: total number of time samples (number of rows in the data field).
  # (3) iEEGsamplingRate: data sampling rate, i.e. the number of data samples 
  #   representing 1 second of EEG data. 
  # (4) channelIndices: an array of the electrode indexes corresponding to 
  #     the columns in the data field.
  df_eeg = as.data.frame(data[[1]][[1]])
  iEEGsamplingRate <- data$dataStruct[[2]][1]
  nSamplesSegment <- data$dataStruct[[3]][1]
  channelIndices <- data$dataStruct[[4]]
  
  i_num_windows <- (nSamplesSegment/iEEGsamplingRate) / secs_per_window
  # samples for a 10 second window
  i_samples_per_window <- iEEGsamplingRate * secs_per_window
  i_window_id <- 1
  
  # data (16 channels, 6 stats) + (id, target, window_id)
  trainset = matrix(,i_num_windows,(16*6+3))
  trainset = as.data.frame(trainset)
  
  # window_id, min, max, mean, median, 
  for (i_offset in seq(1, nSamplesSegment, i_samples_per_window)) {
    v_mean <- c()
    v_median <- c()
    v_sem <- c()
    v_abssum <- c()
    v_min <- c()
    v_max <- c()
    
    # TODO: check for data drop-out, where values across all channels 
    for (i_channel in channelIndices) {
      print(sprintf("window_id: %s, channel_id: %s", i_window_id, i_channel))
      print(sprintf("i_offset: %s, i_offset:end : %s", i_offset, 
                    i_offset+i_samples_per_window))
      i_offset_end <- i_offset+i_samples_per_window
      
      # 10-sec window for channel
      window <- df_eeg[i_offset:i_offset_end, i_channel]
      
      # window stats
      v_min <- c(v_min, min(window))
      v_max <- c(v_max, max(window))
      v_mean <- c(v_mean, mean(window))
      v_median <- c(v_median, median(window))
      v_sem <- c(v_sem, sd(window)/length(window))
      v_abssum <- c(v_abssum, sum(abs(window)))
    }
    trainset[i_window_id,] = c(id,
                              target,
                              i_window_id,
                              v_min,
                              v_max,
                              v_mean,
                              v_median,
                              v_sem,
                              v_abssum)
    
    i_window_id <- i_window_id + 1
  }

  trainset
}


process_file_windows_single <- function(filename) {
  # output
  # trainset = matrix(,1,(16*6+2))
  # trainset = as.data.frame(trainset)
  
  # parse the 'preictal, interical' label

  data = readMat(filename)
  s_base_filename <- basename(filename)
  s_target <- strsplit(fn, "_")[[1]][3]
  s_target <- gsub(".mat", "", s_target)
  
  # calculate window features
  trainset <- process_windows(data = data, id=s_base_filename, target=s_target)
  
  trainset
}

process_windows_parallel <- function(inputdir="../data/train_1",
                             output_filename="../data/features/train_1_windows_set.txt",
                             cores=4) {
  filenames <- list.files(inputdir, pattern="*.mat", full.names=TRUE)
  runtime <- system.time({
    trainset <- mclapply(filenames, process_file_windows_single, mc.cores = cores)
  })[3]
  print(sprintf("runtime: %s", runtime))
  print(sprintf("length: %s", length(trainset)))
  df <- do.call("rbind", trainset)
  
  v_names <- c("id", "target", "window_id",
               paste("min_", 1:16, sep=""), 
               paste("max_", 1:16, sep=""), 
               paste("mean_", 1:16, sep=""), 
               paste("median_", 1:16, sep=""), 
               paste("sem_", 1:16, sep=""), 
               paste("abssum_", 1:16, sep=""))
  print(sprintf("Dimensions: %s", paste(dim(df), collapse="x ")))
  colnames(df) <- v_names
  print(sprintf("Writing features to : %s", output_filename))
  write.table(df, output_filename,
              quote=F,row.names=F,col.names=T,sep="\t")
  df
}


