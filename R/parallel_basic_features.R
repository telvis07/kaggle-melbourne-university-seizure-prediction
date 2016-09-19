
library("parallel")


process_file_single <- function(filename) {
  # output
  # trainset = matrix(,1,(16*6+2))
  # trainset = as.data.frame(trainset)
  
  # parse the 'preictal, interical' label
  targets = strsplit(filename,"/")
  name = strsplit(targets[[1]][4],"_")
  target = strsplit(name[[1]][3],".mat")
  data = readMat(filenames[i])
  df = as.data.frame(data[[1]][[1]])
  
  # summary stats
  sum_list = as.numeric(apply(df, 2, sum))
  mean_list = as.numeric(apply(df, 2, mean))
  median_list = as.numeric(apply(df, 2, median))
  sem_list = as.numeric(apply(df, 2, function (x) {sd(x)/length(x)}))
  abssum_list = as.numeric(apply(df, 2, function (x) {sum(abs(x))}))
  min_list = as.numeric(apply(df, 2, min))
  max_list = as.numeric(apply(df, 2, max)) 
  
  print(target)
  trainset <- c(targets[[1]][4],
                target,
                min_list,
                max_list,
                mean_list,
                median_list,
                sem_list,
                abssum_list)
  trainset
}

process_parallel <- function(inputdir="../data/train_1",
                             output_filename="../data/features/train_1_set.txt",
                             cores=4) {
  filenames <- list.files(inputdir, pattern="*.mat", full.names=TRUE)
  runtime <- system.time({
    trainset <- mclapply(filenames, process_file_single, mc.cores = cores)
  })[3]
  print(sprintf("runtime: %s", runtime))
  print(trainset[[1]][1])
  print(trainset[[1]][2])
  df <- data.frame(matrix(unlist(trainset), 
                          nrow=length(trainset), byrow=T), 
                   stringsAsFactors = F)
  
  v_names <- c("id", "target", 
               paste("min_", 1:16, sep=""), 
               paste("max_", 1:16, sep=""), 
               paste("mean_", 1:16, sep=""), 
               paste("median_", 1:16, sep=""), 
               paste("sem_", 1:16, sep=""), 
               paste("abssum_", 1:16, sep=""))
  print(sprintf("Dimensions: %s", dim(df)))
  colnames(df) <- v_names
  print(sprintf("Writing features to : %s", output_filename))
  write.table(df, output_filename,
              quote=F,row.names=F,col.names=T,sep="\t")
  df
}


