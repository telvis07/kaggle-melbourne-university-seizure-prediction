
# > all_stats_df[ord, c("model", "Recall", "sample_type", "F1")]
# model      Recall sample_type          F1
# 6        rf 0.453947368       smote 0.354983923
# 12      gbm 0.707236842    upsample 0.332946187
# 7       gbm 0.442434211       smote 0.331280788
# 16       rf 0.787828947  downsample 0.328082192
# 17      gbm 0.756578947  downsample 0.320334262
# 13      lda 0.717105263    upsample 0.279397629
# 14      glm 0.703947368    upsample 0.277471637
# 8       lda 0.592105263       smote 0.277349769
# 9       glm 0.577302632       smote 0.275943396
# 18      lda 0.707236842  downsample 0.269592476
# 19      glm 0.685855263  downsample 0.265098538
# 3       lda 0.031250000    nosample 0.056886228
# 11       rf 0.021381579    upsample 0.039513678
# 1        rf 0.011513158    nosample 0.022400000
# 2       gbm 0.009868421    nosample 0.019386107
# 4       glm 0.009868421    nosample 0.018957346
# 5  ensemble 0.003289474    nosample 0.006557377
# 10 ensemble 0.000000000       smote          NA
# 15 ensemble 0.000000000    upsample          NA
# 20 ensemble 0.000000000  downsample          NA

analyze_buncha_models_stats <- function() {
  nosample_stats_df <- read.csv("../stats.sampling.patient_1/patient_1_buncha_model_stats_quick_FALSE_nosampling_SCALE.csv")
  nosample_stats_df$sample_type <- "nosample"
  
  smote_stats_df <- read.csv("../stats.sampling.patient_1/patient_1_buncha_model_stats_quick_FALSE_smote_CARET_SCALE.csv")
  smote_stats_df$sample_type <- "smote"
  
  upsample_stats_df <- read.csv("../stats.sampling.patient_1/patient_1_buncha_model_stats_quick_FALSE_upsample_CARET_SCALE.csv")
  upsample_stats_df$sample_type <- "upsample"
  
  downsample_stats_df <- read.csv("../stats.sampling.patient_1/patient_1_buncha_model_stats_quick_FALSE_downsample_CARET_SCALE.csv")
  downsample_stats_df$sample_type <- "downsample"
  
  all_stats_df <- rbind(nosample_stats_df, smote_stats_df, upsample_stats_df, downsample_stats_df)
  
  ord <- order(all_stats_df$F1, all_stats_df$Recall, decreasing = T)
  all_stats_df[ord, c("model", "Recall", "sample_type", "F1")]
}

