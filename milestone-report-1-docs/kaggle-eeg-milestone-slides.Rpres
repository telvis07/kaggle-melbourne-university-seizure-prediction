Kaggle Seizure Prediction Competition
========================================================
author: Telvis Calhoun
date: February 3, 2017
autosize: true

Kaggle Competition Goal
========================================================

Detect seizure (`preictal`) or non-seizure (`interical`) segments of intracranial electroencephalography (iEEG) data. [See Kaggle EEG Competition page for more details](https://www.kaggle.com/c/melbourne-university-seizure-prediction).

My Approach: 

- Extract basic stats and FFT features for non-overlapping 30-second iEEG windows
- Detect signal drop out and impute missing data with mean for each feature per window
- Predict seizure and non-seizure segments using a stacked model.


Overview of Stacked Model
========================================================

![Kaggle-EEG-Model-Stacking](images/Kaggle-EEG-Model-Stacking.png)

Model Features
========================================================

Random Forest
* Min, Max, Mean, Median, Std Err, Absolute Sum per EEG channel per window
* Kurtosis, Skewness per EEG channel per EEG Channel per window
* Entropy of FFT magnitude and phases
* Entropy, Eignenvalues of EEG channels

***
GLM

* Number of windows per file predicted `preictal` by RF
* Number of windows per file predicted `interictal` by RF

Model Parameters
========================================================

[RandomForest](https://cran.r-project.org/web/packages/randomForest/index.html)
* Number of Trees (ntree) - `50`
* Number of variables randomly sampled as candidates at each split (mtry) - `74`
* Metric - `Accuracy`
* Number training folds - `5`
* Use SMOTE due to `preictal`/`interictal` class imbalance

***

[GLM](https://stat.ethz.ch/R-manual/R-patched/library/stats/html/glm.html)
* Weight for `Preictal` samples: `1.0`
* Weight for `Interictal` samples : `0.25`
* Number training folds - `5`

GLM Predictions
========================================================

* Each point represents a 10-minute EEG segment file
* X-axis is the number of 30-second windows predicted as `preictal` by the RF
* Y-axis is the number of 30-second windows predicted as `interictal` by the RF
* The closer to an axis, the more certain the predicion

***

![GLM Predictions](images/train_1_window_30_quick_FALSE_preds.png)


Final Thoughts
========================================================

* This is my first Kaggle competition. I acheived my goal of making a competition submission. [See my profile](https://www.kaggle.com/telvis).
* I submitted after the deadline but my submission would have ranked 391 of 2440 submissions. [Screenshot](./images/kaggle_eeg_submission_capture.png) 
* The [Kaggle model stacking tutorial](http://blog.kaggle.com/2016/12/27/a-kagglers-guide-to-model-stacking-in-practice/) helped me understand cross-fold validation with stacked models.
* [Deep's Kernel](https://www.kaggle.com/deepcnn/melbourne-university-seizure-prediction/feature-extractor-matlab2python-translated) and [Tony Reina's Kernel](https://www.kaggle.com/treina/melbourne-university-seizure-prediction/feature-extractor-matlab2python-translated) helped me understand EEG features. 

