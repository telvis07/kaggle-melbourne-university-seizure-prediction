# extract features for each window in each file.
library("parallel")
library("R.matlab")
library("ggplot2")
library("gridExtra")
library("dplyr")

process_windows <- function(data, id, target, n_seg_num, secs_per_window=10) {
  # (2) nSamplesSegment: total number of time samples (number of rows in the data field).
  # (3) iEEGsamplingRate: data sampling rate, i.e. the number of data samples 
  #   representing 1 second of EEG data. 
  # (4) channelIndices: an array of the electrode indexes corresponding to 
  #     the columns in the data field.
  df_eeg = as.data.frame(data[[1]][[1]])
  iEEGsamplingRate <- data$dataStruct[[2]][1]
  nSamplesSegment <- data$dataStruct[[3]][1]
  channelIndices <- data$dataStruct[[4]]
  n_num_channels <- length(channelIndices)
  
  n_num_windows <- (nSamplesSegment/iEEGsamplingRate) / secs_per_window
  
  # nsamples for a 10 second window
  n_samples_per_window <- iEEGsamplingRate * secs_per_window
  n_window_id <- 0
  
  # data (16 channels, 6 stats) + (id, target, window_id, segnum)
  trainset = matrix(,n_num_windows,(16*6+4))
  trainset = as.data.frame(trainset)
  
  # window_id, min, max, mean, median, 
  for (n_offset in seq(1, nSamplesSegment, n_samples_per_window)) {
    n_window_id <- n_window_id + 1
    
    v_mean <- c()
    v_median <- c()
    v_sem <- c()
    v_abssum <- c()
    v_min <- c()
    v_max <- c()
    
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
      trainset[n_window_id,] = c(id,
                                 target,
                                 n_window_id,
                                 n_seg_num,
                                 rep(NA, n_num_channels),
                                 rep(NA, n_num_channels),
                                 rep(NA, n_num_channels),
                                 rep(NA, n_num_channels),
                                 rep(NA, n_num_channels),
                                 rep(NA, n_num_channels))
      next
      
    } 
    
    for (n_channel in channelIndices) {
      # 10-sec window for channel
      window <- df_eeg[n_offset:n_offset_end, n_channel]
      
      if( n_dropout_rows > 0) {
        # remove dropout rows
        window <- window[! v_dropout_rows]
      }
      
      if (length(window) != (n_samples_per_window - n_dropout_rows)) {
        print(sprintf("window_id: %s, channel_id: %s", n_window_id, n_channel))
        print(sprintf("n_offset: %s, n_offset:end : %s", n_offset, 
                      n_offset+n_samples_per_window))
        stop(sprintf("[process_windows] - invalid window length. Expected: %s, Found: %s",
                     n_samples_per_window - n_dropout_rows,
                     length(window)))
      }
      
      # window stats
      v_min <- c(v_min, min(window))
      v_max <- c(v_max, max(window))
      v_mean <- c(v_mean, mean(window))
      v_median <- c(v_median, median(window))
      v_sem <- c(v_sem, sd(window)/length(window))
      v_abssum <- c(v_abssum, sum(abs(window)))
    }
    trainset[n_window_id,] = c(id,
                               target,
                               n_window_id,
                               n_seg_num,
                               v_min,
                               v_max,
                               v_mean,
                               v_median,
                               v_sem,
                               v_abssum)
    
  }

  trainset
}

process_file_windows_single <- function(filename) {
  # parse the 'preictal, interical' label
  data = readMat(filename)
  s_base_filename <- basename(filename)
  v_filename_parts <- strsplit(s_base_filename, "_")[[1]]
  s_target <- v_filename_parts[3]
  s_target <- gsub(".mat", "", s_target)
  n_seg_num <- as.numeric(v_filename_parts[2])
  
  # calculate window features
  trainset <- process_windows(data = data, id=s_base_filename, target=s_target, n_seg_num=n_seg_num)
  
  trainset
}

process_windows_parallel <- function(inputdir="../data/train_1.small/",
                             output_filename="../data/features/train_1_windows_set.txt",
                             cores=4) {
  filenames <- list.files(inputdir, pattern="*.mat", full.names=TRUE)
  runtime <- system.time({
    trainset <- mclapply(filenames, process_file_windows_single, mc.cores = cores)
  })[3]
  print(sprintf("runtime: %s", runtime))
  print(sprintf("length: %s", length(trainset)))
  df <- do.call("rbind", trainset)
  
  v_names <- c("id", "target", "window_id", "segnum",
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

load_window_features <- function(output_filename="../data/features/train_1_windows_set.txt") {
  df <- read.table(output_filename, header=T, stringsAsFactors = F)
  df$segnum <- as.numeric(df$segnum)
  df$target <- as.numeric(df$target)

  # sort by target then segment number
  ord = order(df$target, df$segnum, decreasing = F)
  df <- df[ord,]
  df
}


plot_feature_for_channel <- function(df, offset, channel=1, label="min") {
  n_colnum <- offset + channel
  n_mean <- mean(df[,n_colnum], na.rm = T)
  n_std <- sd(df[,n_colnum], na.rm = T)
  datums <- data.frame(channel_data=df[,n_colnum],
                       times=seq(length(v_channel_data)),
                       target=df$target)

  # datums <- transform(datums, channel_data=(channel_data-n_mean)/n_std)
  obj <- ggplot(datums, aes(x=times, y=channel_data, color=target)) + 
    geom_line() + 
    ylab(sprintf("%s-%s", label, channel)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position="none")
  obj
}

plot_feature_min_channels <- function(df, num_channels=16, output_png="min_eeg_channels.png") {
  offset=4
  l_plots <- list()
  
  for (n_channel in seq(num_channels)) {
    plt <- plot_feature_for_channel(df, offset = offset, channel = n_channel, label="min")
    l_plots <- c(l_plots, list(plt))
  }
  
  print(length(l_plots))
  do.call(grid.arrange, c(l_plots, top="MIN per iEEG channel ('interical', 'preictal')"))
  dev.copy(png, width = 960, height = 960, units = "px", output_png)
  dev.off()
  l_plots
}

plot_feature_max_channels <- function(df, num_channels=16, output_png="max_eeg_channels.png") {
  offset=20
  l_plots <- list()
  
  for (n_channel in seq(num_channels)) {
    plt <- plot_feature_for_channel(df, offset = offset, channel = n_channel, label = "max" )
    l_plots <- c(l_plots, list(plt))
  }
  
  print(length(l_plots))
  do.call(grid.arrange, c(l_plots, top="MAX per iEEG channel ('interical', 'preictal')"))
  dev.copy(png, width = 960, height = 960, units = "px", output_png)
  dev.off()
  l_plots
}

plot_feature_mean_channels <- function(df, num_channels=16, output_png="mean_eeg_channels.png") {
  offset=36
  l_plots <- list()
  
  for (n_channel in seq(num_channels)) {
    plt <- plot_feature_for_channel(df, offset = offset, channel = n_channel, label="mean" )
    l_plots <- c(l_plots, list(plt))
  }
  
  print(length(l_plots))
  do.call(grid.arrange, c(l_plots, top="MEAN per iEEG channel ('interical', 'preictal')"))
  dev.copy(png, width = 960, height = 960, units = "px", output_png)
  dev.off()
  l_plots
}

