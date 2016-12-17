source("correlation_features.R")


run_features <- function(){
  b_quick = T
  
  if (b_quick) {
    filename="../data/train_1/1_999_0.mat"
    df_1 <- process_file_windows_single(filename)
    
    filename="../data/train_1/1_1_1.mat"
    df_2 <- process_file_windows_single(filename)
    
    df <- rbind(df_1, df_2)
  } else {
    window_size = 30
    output_filename=sprintf("../data/features/train_2_window_%s_secs_correlation_and_fft.testing.txt", window_size)
    df <- process_windows_parallel(cores=4, 
                                   inputdir="../data/train_2/",
                                   output_filename=output_filename,
                                   secs_per_window=window_size)
  }
}

run_features.train <- function(){
  
  for (patient_id in c(2,3)){
    print(sprintf("Generating features for patient_id : %s", patient_id))
    window_size = 30
    output_filename=sprintf("../data/features/train_%s_window_%s_secs_correlation_and_fft.testing.txt", patient_id, window_size)
    print(sprintf("Generating : %s", output_filename))
    df <- process_windows_parallel(cores=4, 
                                   inputdir=sprintf("../data/train_%s/", patient_id),
                                   output_filename=output_filename,
                                   secs_per_window=window_size)
  }
}


run_features.test <- function(){
  
  for (patient_id in c(1,2,3)){
    print(sprintf("Generating features for patient_id : %s", patient_id))
    window_size = 30
    output_filename=sprintf("../data/features/test_%s_new_window_%s_secs_correlation_and_fft.testing.txt", patient_id, window_size)
    print(sprintf("Generating : %s", output_filename))
    df <- process_windows_parallel(cores=4, 
                                   inputdir=sprintf("../data/test_%s_new/", patient_id),
                                   output_filename=output_filename,
                                   secs_per_window=window_size)
  }
}