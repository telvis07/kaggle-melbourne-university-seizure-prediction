library(dplyr)

sample_data <- function(trainset, n_neg_samples=0, n_pos_samples=0) {
  # make sure we don't subsample rows with all NULLs
  trainset <- trainset[rowSums(is.na(trainset)) == 0,]
  
  # separate pos and negative class
  neg_samples <- trainset[trainset$target == "interictal",]
  pos_samples <- trainset[trainset$target == "preictal",]
  print(sprintf("before downsample: neg dim: %s", dim(neg_samples)[1] ))
  print(sprintf("before downsample: pos dim: %s", dim(pos_samples)[1] ))
  
  # downsample negative class (interictal)
  if ((n_neg_samples != 0) && (n_neg_samples < dim(neg_samples)[1])) {
    sample_inds <- sample(dim(neg_samples)[1], n_neg_samples, replace=FALSE)
    neg_samples <- neg_samples[sample_inds,]
  }
  
  # downsample positive class (preictal)
  if((n_pos_samples != 0) && (n_pos_samples < dim(pos_samples)[1])){
    sample_inds <- sample(dim(pos_samples)[1], n_pos_samples, replace=FALSE)
    pos_samples <- pos_samples[sample_inds,]
  }
  
  print(sprintf("after downsample: neg dim: %s", dim(neg_samples)[1] ))
  print(sprintf("after downsample: pos dim: %s", dim(pos_samples)[1] ))
  
  trainset <- rbind(neg_samples, pos_samples)
  
  print(sprintf("downsampled trainset: %s rows, %s cols", dim(trainset)[1],
                dim(trainset)[2]))
  print(table(trainset$target))
  
  trainset
}