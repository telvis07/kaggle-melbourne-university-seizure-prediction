# From: https://www.kaggle.com/c/melbourne-university-seizure-prediction/forums/t/23332/start-with-r-generate-basic-features/133838

# use R.matlab to read mat files and 
# extract min, max, mean, median, sem and 
# sum of abs value for 16 channel.

#Here is the example for train_1. 
# Use these 96 features and xgboost model, AUC = 0.56

library("R.matlab")

generate_basic_features <- function(){
  filenames <- list.files("../data/train_1", pattern="*.mat", full.names=TRUE)
  
  # 16 channels
  # times (6 summary statistics + 2 cols(id, target))
  trainset = matrix(,length(filenames),(16*6+2))
  
  trainset = as.data.frame(trainset)
  
  colnames(trainset) = c("id",
                         "target",
                         paste("min_",rep(1:16),sep=""),
                         paste("max_",rep(1:16),sep=""),
                         paste("mean_",rep(1:16),sep=""),
                         paste("median_",rep(1:16),sep=""),
                         paste("sem_",rep(1:16),sep=""),
                         paste("abssum_",rep(1:16),sep=""))
  
  for (i in 1: length(filenames)){
    
    targets = strsplit(filenames[i],"/")
    
    name = strsplit(targets[[1]][4],"_")
    
    target = strsplit(name[[1]][3],".mat")
    
    data = readMat(filenames[i])
    
    df = as.data.frame(data[[1]][[1]])
    
    sum_list = as.numeric(apply(df, 2, sum))
    
    mean_list = as.numeric(apply(df, 2, mean))
    
    median_list = as.numeric(apply(df, 2, median))
    
    sem_list = as.numeric(apply(df, 2, function (x) {sd(x)/length(x)}))
    
    abssum_list = as.numeric(apply(df, 2, function (x) {sum(abs(x))}))
    
    min_list = as.numeric(apply(df, 2, min))
    
    max_list = as.numeric(apply(df, 2, max)) 
    
    trainset[i,] = c(targets[[1]][4],target,min_list,max_list,mean_list,median_list,sem_list,abssum_list)
  }
  
  write.table(trainset,"train_1_set.txt",quote=F,row.names=F,col.names=T,sep="\t")
}