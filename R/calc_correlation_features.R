source("correlation_features.R")

b_quick = F

if (b_quick) {
  filename="../data/train_1/1_999_0.mat"
  df_1 <- process_file_windows_single(filename)
  
  filename="../data/train_1/1_1_1.mat"
  df_2 <- process_file_windows_single(filename)
  
  df <- rbind(df_1, df_2)
} else {
  window_size = 60
  output_filename=sprintf("../data/features/train_1_window_%s_secs_correlation_and_fft.testing.txt", window_size)
  df <- process_windows_parallel(cores=4, 
                                 inputdir="../data/train_1.small/",
                                 output_filename=output_filename,
                                 secs_per_window=window_size)
}