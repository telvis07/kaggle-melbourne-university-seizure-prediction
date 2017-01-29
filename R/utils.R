get_target_from_id <- function(filename) {
  s_base_filename <- basename(filename)
  v_filename_parts <- strsplit(s_base_filename, "_")[[1]]
  s_target <- v_filename_parts[3]
  s_target <- gsub(".mat", "", s_target)
  n_seg_num <- as.numeric(v_filename_parts[2])
  as.numeric(s_target)
}

calc_perf_stats <- function(cm){
  # cm : confusion matrix
  
  # calc stats
  n = sum(cm) # number of instances
  nc = nrow(cm) # number of classes
  diag = diag(cm) # number of correctly classified instances per class
  rowsums = apply(cm, 1, sum) # number of instances per class
  colsums = apply(cm, 2, sum) # number of predictions per class
  p = rowsums / n # distribution of instances over the actual classes
  q = colsums / n # distribution of instances over the predicted classes
  
  # accuracy
  accuracy = sum(diag) / n
  
  # precision/recall/f1
  precision = diag / colsums
  recall = diag / rowsums
  f1 = 2 * precision * recall / (precision + recall)
  
  # add to data.frame 
  df_stats <- data.frame(precision=precision[2],
                         recall=recall[2],
                         F1=f1[2])
  
  df_stats
}

save_plot <- function(obj, filename="plots/default.png") {
  png(filename, width = 960, height = 960, units = "px")
  print(obj)
  dev.off()
  print(sprintf("Saved %s", filename))
}