# Soil_classification
Use of prediction modeling in classification of soil from multiple environments.
Four models are tested for classification and predictabiliy.
Models tested are multinomial model using 
  i) all the soil features hereby called as base 
  ii) first two principal components hereby called as PC 
  iii) radom forest model abbreviated as RF
  iv) self-organizing map abbreviated as SOM
  
:There are two .R files
  i) MET_soil_models_initiation.R  - has scripts for initial quality checking and estimation and prediction evaluation of the dataset 
  ii) MET_soil_models_prediciton.R  - has scripts for testing the predictability of the models under multiple scenarios.
  
:The required packages are 
  library(ggplot2)
  library(reshape2)
  library(psych) # for correlation matrix visualization
  library(devtools) # required for ggbiplot downloading
  library(ggbiplot) # for biplot
  library(nnet) # for multinomial logistic regression
  library(caret) # for confustion matrix
  # library(pROC) # for ROC curve - can use for multi-class as well
  library(randomForest) # for RandomForest model
  library(kohonen) # for self organizing maps
  library(factoextra) # to extract and visualise results of multivariate analysis
  library(cvms) # to plot confusion matrix from tibble
  library(cluster) # for k-medoid clustering from unsupervised RF
  library(gtools) # for combinations of variables 
  library(dplyr) # for grouping and applying functions across groups
