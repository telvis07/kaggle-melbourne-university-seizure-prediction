library(dplyr)

generate_training_label_file <- function(inputdir="../data/train_1.small/", outputcsv="train_1_small_labels.csv") {
  filenames <- list.files(inputdir, pattern="*.mat", full.names=TRUE)
  
  # rows = number of filenames. cols = (id, segnum, target)
  m_datums <- matrix(NA, length(filenames), 3)
  for (i in seq(length(filenames))){
    filename = filenames[i]
    s_base_filename <- basename(filename)
    v_filename_parts <- strsplit(s_base_filename, "_")[[1]]
    n_seg_num <- v_filename_parts[2]
    s_target <- v_filename_parts[3]
    
    s_target <- gsub(".mat", "", s_target)
    print(c(s_base_filename, n_seg_num, s_target))
    m_datums[i,] <- c(s_base_filename, n_seg_num, as.numeric(s_target))
  }
  
  print(head(m_datums))
  
  df_datums <- as.data.frame(m_datums, stringsAsFactors = F)
  names(df_datums) <- c("id", "segnum", "target")
  df_datums <- mutate(df_datums, segnum=as.numeric(segnum), target=as.numeric(target))
  # df_datums$target <- factor(df_datums$target, levels=c(0,1), labels=c('interictal', 'preictal'))
  ord <- order(df_datums$target, df_datums$segnum)
  df_datums <- df_datums[ord,]
  
  write.csv(df_datums, outputcsv, row.names = F)
  print(sprintf("Wrote : %s", outputcsv))
  df_datums
}