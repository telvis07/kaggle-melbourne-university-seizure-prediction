source("correlation_features.R")

b_quick = T

if (b_quick) {
  filename="../data/train_1/1_999_0.mat"
  df <- process_file_windows_single(filename)
  
  filename="../data/train_1/1_1_1.mat"
  df <- process_file_windows_single(filename)
} else {
  output_filename="../data/features/train_1_correlation_and_fft.txt"
  df <- process_windows_parallel(cores=4, 
                                 inputdir = "../data/train_1/",
                                 output_filename=output_filename)
}