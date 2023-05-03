### Multi-environmental soil predictions/extrapolations ###
### Ani A. Elas 	February 2023 ###


# prediction modeling on reduced number of features, training data size, integrated segments, and cross-checking of clusters


# required libraries

library(ggplot2)
library(reshape2)
library(psych) # for correlation matrix visualization
# install.packages("devtools", type = "win.binary") # if devtools are not loaded properly

library(devtools) # required for ggbiplot downloading
install_github("vqv/ggbiplot")
library(ggbiplot) # for biplot
library(nnet) # for multinomial logistic regression
library(caret) # for confustion matrix
# library(pROC) # for ROC curve - can use for multi-class as well
library(randomForest) # for RandomForest model
library(kohonen) # for self organizing maps
library(factoextra) # to extract and visualise results of multivariate analysis
# library(pmclust) # various algorithms for clustering
# library(tibble) # to create tibble
library(cvms) # to plot confusion matrix from tibble
library(cluster) # for k-medoid clustering from unsupervised RF

library(gtools) # for combinations of variables 
library(dplyr) # for grouping and applying functions across groups



setwd() # set working directory

data1 <- read.csv("data1.csv", h=T, sep=",") # data contains segments, plots, blocks, and locations

data1b <- na.omit(data1) # influential points for features are removed


### grouping of segments, plots, blocks and CV the cluster
 ## check with all the four models each cluster

dataC_base <- read.csv("MET_soil_cluster_base.csv", h=T, sep= ",") # data with cluster label from base model
dataC_PC <- read.csv("MET_soil_cluster_PC.csv", h=T, sep= ",") # data with cluster label from PC model
dataC_RF <- read.csv("MET_soil_cluster_RF.csv", h=T, sep= ",") # data with cluster label from RF model
dataC_SOM <- read.csv("MET_soil_cluster_SOM.csv", h=T, sep= ",") # data with cluster label from SOM model


dataC <- cbind(data1b, dataC_base[,"cluster"], dataC_PC[,"cluster"],
              dataC_RF[,"cluster"], dataC_SOM[,"cluster"])

colnames(dataC)[23:26] <- c("ClusterB", "ClusterP", "ClusterR", "ClusterS")

# merging multiple clusters
dataC$Plot3 <- dataC$Plot 
dataC$Plot3[dataC$Plot3 == '2'] <- '1'
dataC$Plot3[dataC$Plot3 == '3'] <- '2'
dataC$Plot3[dataC$Plot3 == '4'] <- '2'
dataC$Plot3[dataC$Plot3 == '5'] <- '3'
dataC$Plot3[dataC$Plot3 == '6'] <- '3'

dataC$Plot4 <- dataC$Plot
dataC$Plot4[dataC$Plot4 == '2'] <- '1'
dataC$Plot4[dataC$Plot4 == '3'] <- '1'
dataC$Plot4[dataC$Plot4 == '4'] <- '2'
dataC$Plot4[dataC$Plot4 == '5'] <- '2'
dataC$Plot4[dataC$Plot4 == '6'] <- '2'




dataC_plot <- dataC %>% group_by(Plot, Block, Location) %>%
            summarise(across(c(Moisture, pH, Salinity, Sand, Clay, Silt, Carbon,
                              Phosphorous, Nitrogen, Sulphur, Sodium, Calcium, 
                              Potassium, ClusterB, ClusterP, ClusterR, ClusterS), mean),
            .groups = 'drop') %>%
            as.data.frame()
# rounding the cluster to nearest whole number (for example, 2 out of 3 segments belong to same group)
dataC_plot$ClusterB <- round(dataC_plot$ClusterB, 0)
dataC_plot$ClusterP <- round(dataC_plot$ClusterP, 0)
dataC_plot$ClusterR <- round(dataC_plot$ClusterR, 0)
dataC_plot$ClusterS <- round(dataC_plot$ClusterS, 0)



dataC_plot3 <- dataC %>% group_by(Plot3, Block, Location) %>%
            summarise(across(c(Moisture, pH, Salinity, Sand, Clay, Silt, Carbon,
                              Phosphorous, Nitrogen, Sulphur, Sodium, Calcium, 
                              Potassium, ClusterB, ClusterP, ClusterR, ClusterS), mean),
            .groups = 'drop') %>%
            as.data.frame()

dataC_plot3$ClusterB <- as.factor(round(dataC_plot3$ClusterB, 0))
dataC_plot3$ClusterP <- as.factor(round(dataC_plot3$ClusterP, 0))
dataC_plot3$ClusterR <- as.factor(round(dataC_plot3$ClusterR, 0))
dataC_plot3$ClusterS <- as.factor(round(dataC_plot3$ClusterS, 0))


dataC_plot4 <- dataC %>% group_by(Plot4, Block, Location) %>%
            summarise(across(c(Moisture, pH, Salinity, Sand, Clay, Silt, Carbon,
                              Phosphorous, Nitrogen, Sulphur, Sodium, Calcium, 
                              Potassium, ClusterB, ClusterP, ClusterR, ClusterS), mean),
            .groups = 'drop') %>%
            as.data.frame()
dataC_plot4$ClusterB <- as.factor(round(dataC_plot4$ClusterB, 0))
dataC_plot4$ClusterP <- as.factor(round(dataC_plot4$ClusterP, 0))
dataC_plot4$ClusterR <- as.factor(round(dataC_plot4$ClusterR, 0))
dataC_plot4$ClusterS <- as.factor(round(dataC_plot4$ClusterS, 0))


dataC_segment <- dataC[,c(1,4,5, 9:11, 13:26)]
dataC_segment$ClusterB <- as.factor(dataC_segment$ClusterB)
dataC_segment$ClusterP <- as.factor(dataC_segment$ClusterP)
dataC_segment$ClusterR <- as.factor(dataC_segment$ClusterR)
dataC_segment$ClusterS <- as.factor(dataC_segment$ClusterS)



~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




## applying models on merged clusters ###
# modeling done on variation in plot area 

# two functions are written below
  # Soil_sim(...) is using all the variables 
  # Soil_sim_small(...) is using only the five most important variables based on models Na, Ca, K, C, P, and Sand




# SOM grid
data <- dataC_plot # for example. The clusters should be factors in the data
  # data should have c(1:3) as c(Plot, Block, Location) and c(17:20) as clusters
    # exception is for the original data with segment included 
dataP.grid <- somgrid(xdim = 10, ydim = 10, topo = "hexagonal", toroidal = T, neighbourhood.fct =  "gaussian") # for full data    
dataP.grid <- somgrid(xdim = 6, ydim = 6, topo = "hexagonal", toroidal = T, neighbourhood.fct =  "gaussian") # for by plot data
dataP.grid <- somgrid(xdim = 4, ydim = 4, topo = "hexagonal", toroidal = T, neighbourhood.fct =  "gaussian") # for by three plots per block data
dataP.grid <- somgrid(xdim = 2, ydim = 2, topo = "hexagonal", toroidal = T, neighbourhood.fct =  "gaussian") # for by two plots per block data

train_test <- c(0.8, 0.2) # training proportions as c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95)

# imp: change the end file beginning names for saving with the correct name.

Soil_sim <- function(data, train_test, SOM_grid){


cfm1 <- cfm2 <- cfm3 <- cfm4 <- cfm5 <- cfm6 <- cfm7 <- cfm8 <- cfm9 <- cfm10 <- cfm11 <- cfm12 <- list()
 cfm13 <- cfm14 <- cfm15 <- cfm16 <- list()
# cfm <- list() # can be used which will be recycled through all 16 analysis

AccuracyB <- AccuracyP <-  AccuracyR <- AccuracyS  <- AccuracyE <- AccuracyEP <- AccuracyER <- AccuracyES <- matrix(NA, 100, 4) 
colnames(AccuracyB) <- colnames(AccuracyP) <- colnames(AccuracyR) <- colnames(AccuracyS) <- c("Base", "PC", "RF", "SOM")
colnames(AccuracyE) <- colnames(AccuracyEP) <- colnames(AccuracyER) <- colnames(AccuracyES) <- c("Base", "PC", "RF", "SOM")

RuntimeB <- RuntimeP <- RuntimeR <- RuntimeS <- matrix(NA, 100, 4)
colnames(RuntimeB) <- colnames(RuntimeP) <- colnames(RuntimeR) <- colnames(RuntimeS) <- c("Base", "PC", "RF", "SOM")


for (i in 1:100){
  ind <- sample(2, nrow(data),
          replace = TRUE,
          prob = train_test)
  training <- data[ind==1,]
  testing <- data[ind==2,]

  # standardizing the data
  train1 <- as.matrix(scale(training[,c(4:16)]))
  test1 <- scale(testing[,c(4: 16)], center = attr(train1, "scaled:center")) 

  train2 <- cbind(training[,c(1:3,17:20)], train1)
  test2 <- cbind(testing[,c(1:3,17:20)], test1)

  ### checking the validity of ClusterB with other models

  ## base model
  start_time <- Sys.time()
  base.tr <- multinom(ClusterB ~ Moisture + pH + Salinity + Sand + Clay + Silt + Carbon +
                       Phosphorous +  Nitrogen + Sulphur + Sodium + Calcium + Potassium, train2, trace = FALSE)
  # if (class(base.tr) != "try-error"){
  p.tst1 <- predict(base.tr, test2)
  cfm1[[i]] <- table(Predict = p.tst1, Actual = test2$ClusterB)
  a <- try(tr(cfm1[[i]])/sum(cfm1[[i]]))
  	if (class(a) != "try-error"){
  		AccuracyB[i, 1] <- a
  		end_time <- Sys.time()
  		RuntimeB[i, 1] <- end_time - start_time
  	} else {
  		AccuracyE[i, 1] <- 1
  	}
  # } else {
  	# AccuracyE[i, 1] <- 1
  # }
  



  ## PC model
  start_time <- Sys.time()

  pc.tr <- prcomp(train2[,c(8:20)]) 
  # obtain PC values for all the rows and dependent variable to the data
  trg <- predict(pc.tr, train2) # predict() in PCA: predict the projection of new rows with PCA.
  trg <- as.data.frame(cbind(trg, train2$ClusterB)) # adding the response variable to data
  colnames(trg)[14] <- "ClusterB"
  tst <- predict(pc.tr, test2)
  tst <- as.data.frame(cbind(tst, test2$ClusterB))
  colnames(tst)[14] <- "ClusterB"

  model1 <- multinom(ClusterB ~ PC1+PC2, data = trg)

  # if (class(model1) != "try-error"){ 
  # Confusion matrix and misclassification error - testing data
  p.tst2 <- predict(model1, tst)
  cfm2[[i]] <- table(Predict = p.tst2, Actual = tst$ClusterB) # confusion matrix
  a <- try(tr(cfm2[[i]])/sum(cfm2[[i]]))
  	if (class(a) != "try-error"){
  		AccuracyB[i, 2] <- a
  		end_time <- Sys.time()
  		RuntimeB[i, 2] <- end_time - start_time
  	} else {
  		AccuracyE[i, 2] <- 1
  	}
  # } else {
  	# AccuracyE[i, 2] <- 1
  # }

  
  	


  ## RF model
  start_time <- Sys.time()
  rf.tr_cl <- try(randomForest(x = train2[,c(8:20)], 
            y = as.factor(train2$ClusterB), ntree = 1000))

if (class(rf.tr_cl) != "try-error") {
  # Predicting the Test set results
  trg = predict(rf.tr_cl, newdata = test2[,c(8: 20)])
  
  # Confusion Matrix
  cfm3[[i]] = table(Predict = trg, Actual = test2$ClusterB)
  a <- try(tr(cfm3[[i]])/sum(cfm3[[i]]))
  	if (class(a) != "try-error"){
  		AccuracyB[i, 3] <- a
  		end_time <- Sys.time()
  		RuntimeB[i, 3] <- end_time - start_time
  	} else {
  		AccuracyE[i, 3] <- 1
  	}
} else {
	# AccuracyE[i, 3] <- 1
}




  ## SOM model
  start_time <- Sys.time()  
  testXY <- list(independent = test1, dependent = factor(testing$ClusterB))

  train_som <- try(xyf(train1, classvec2classmat(factor(training$ClusterB)), SOM_grid, rlen = 1000, radius = 2.5, 
                      keep.data = TRUE, dist.fcts = "euclidean"))

  if (class(train_som) != "try-error"){
  # testing the model
  test_cluster <- predict(train_som, newdata = as.matrix(test1), whatmap = 1)
  cfm4[[i]] <- try(table(Predict = test_cluster$predictions[[2]], Actual = testing$ClusterB))
  if (class(cfm4[[i]]) != "try-error"){
  	a <- try(tr(cfm4[[i]])/sum(cfm4[[i]]))
  	if (class(a) != "try-error"){
  		AccuracyB[i, 4] <- a
  		end_time <- Sys.time()
  		RuntimeB[i, 4] <- end_time - start_time
  	} else {
  		AccuracyE[i, 4] <- 1
  	}
  } else {
  	AccuracyE[i, 4] <- 1
  }

  } else {
  	AccuracyE[i, 4] <- 1
  }
  
  	



### checking the validity of ClusterP with other models

  ## base model
  start_time <- Sys.time()
  base.tr <- multinom(ClusterP ~ Moisture + pH + Salinity + Sand + Clay + Silt + Carbon +
                       Phosphorous +  Nitrogen + Sulphur + Sodium + Calcium + Potassium, train2, trace = FALSE)

# if(class(base.tr) != "try-error"){
  p.tst3 <- predict(base.tr, test2)
  cfm5[[i]] <- table(Predict = p.tst3, Actual = test2$ClusterP)
  a <- try(tr(cfm5[[i]])/sum(cfm5[[i]]))
  	if (class(a) != "try-error"){
  		AccuracyP[i, 1] <- a
  		end_time <- Sys.time()
  		RuntimeP[i, 1] <- end_time - start_time
  	} else {
  		AccuracyEP[i, 1] <- 1
  	}
# } else {
   # AccuracyEP[i, 1] <- 1
# }
   


  ## PC model
  start_time <- Sys.time()
  pc.tr <- prcomp(train2[,c(8:20)])

  # obtain PC values for all the rows and dependent variable to the data
  trg <- predict(pc.tr, train2) # predict() in PCA: predict the projection of new rows with PCA.
  trg <- as.data.frame(cbind(trg, train2$ClusterP)) # adding the response variable to data
  colnames(trg)[14] <- "ClusterP"
  tst <- predict(pc.tr, test2)
  tst <- as.data.frame(cbind(tst, test2$ClusterP))
  colnames(tst)[14] <- "ClusterP"

  model1 <- multinom(ClusterP ~ PC1+PC2, data = trg)

  # if (class(model1) != "try-error"){
  p.tst4 <- predict(model1, tst)
  # Confusion matrix and misclassification error - testing data
  cfm6[[i]] <- table(Predict = p.tst4, Actual = tst$ClusterP) # confusion matrix
  a <- try(tr(cfm6[[i]])/sum(cfm6[[i]]))
  	if (class(a) != "try-error"){
  		AccuracyP[i, 2] <- a
  		end_time <- Sys.time()
  		RuntimeP[i, 2] <- end_time - start_time
  	} else {
  		AccuracyEP[i, 2] <- 1
  	}	
  # } else {
  	AccuracyEP[i, 2] <- 1
  # }

 
 


  ## RF model
  start_time <- Sys.time()
  rf.tr_cl <- try(randomForest(x = train2[,c(8:20)], 
            y = as.factor(train2$ClusterP), ntree = 1000))

  if (class(rf.tr_cl) != "try-error") {
  # Predicting the Test set results
  trg = predict(rf.tr_cl, newdata = test2[,c(8: 20)])
  
  # Confusion Matrix
  cfm7[[i]] = table(Predict = trg, Actual = test2$ClusterP)
  a <- try(tr(cfm7[[i]])/sum(cfm7[[i]]))
  	if (class(a) != "try-error"){
  		AccuracyP[i, 3] <- a
  		end_time <- Sys.time()
  		RuntimeP[i, 3] <- end_time - start_time
  	} else {
  		AccuracyEP[i, 3] <- 1
  	}
  } else {
  	AccuracyEP[i, 3] <- 1
  }
  


  ## SOM model
  start_time <- Sys.time()  
  testXY <- list(independent = test1, dependent = factor(testing$ClusterP))

  train_som <- try(xyf(train1, classvec2classmat(factor(training$ClusterP)), SOM_grid, rlen = 1000, radius = 2.5, 
                      keep.data = TRUE, dist.fcts = "euclidean"))

 if (class(train_som) != "try-error") {
  # testing the model
  test_cluster <- predict(train_som, newdata = as.matrix(test1), whatmap = 1)
  cfm8[[i]] <- try(table(Predict = test_cluster$predictions[[2]], Actual = testing$ClusterP))
  if (class(cfm8[[i]]) != "try-error"){
  	a <- try(tr(cfm8[[i]])/sum(cfm8[[i]]))
  	if (class(a) != "try-error"){
  		AccuracyP[i, 4] <- a
  		end_time <- Sys.time()
  		RuntimeP[i, 4] <- end_time - start_time
  	} else {
  		AccuracyEP[i, 4] <- 1
  	}
  } else {
  	AccuracyEP[i, 4] <- 1
  }
 } else {
  	AccuracyEP[i, 4] <- 1
  }

  
  


### checking the validity of ClusterR with other models

  ## base model
  start_time <- Sys.time()
  base.tr <- multinom(ClusterR ~ Moisture + pH + Salinity + Sand + Clay + Silt + Carbon +
                       Phosphorous +  Nitrogen + Sulphur + Sodium + Calcium + Potassium, train2, trace = FALSE)

  # if (class(base.tr) != "try-error") {
  p.tst5 <- predict(base.tr, test2)
  cfm9[[i]] <- table(Predict = p.tst5, Actual = test2$ClusterR)
  a <- try(tr(cfm9[[i]])/sum(cfm9[[i]]))
  	if (class(a) != "try-error"){
  		AccuracyR[i, 1] <- a
  		end_time <- Sys.time()
  		RuntimeR[i, 1] <- end_time - start_time
  	} else {
  		AccuracyER[i, 1] <- 1
  	}
  # } else {
  	# AccuracyER[i, 1] <- 1
  # }
 



  ## PC model
  start_time <- Sys.time()
  pc.tr <- prcomp(train2[,c(8:20)])

  # obtain PC values for all the rows and dependent variable to the data
  trg <- predict(pc.tr, train2) # predict() in PCA: predict the projection of new rows with PCA.
  trg <- as.data.frame(cbind(trg, train2$ClusterR)) # adding the response variable to data
  colnames(trg)[14] <- "ClusterR"
  tst <- predict(pc.tr, test2)
  tst <- as.data.frame(cbind(tst, test2$ClusterR))
  colnames(tst)[14] <- "ClusterR"

  model1 <- multinom(ClusterR ~ PC1+PC2, data = trg)

  # if (class(model1) != "try-error") {

  # Confusion matrix and misclassification error - testing data
  p.tst6 <- predict(model1, tst)
  cfm10[[i]] <- table(Predict = p.tst6, Actual = tst$ClusterR) # confusion matrix
  a <- try(tr(cfm10[[i]])/sum(cfm10[[i]]))
  	if (class(a) != "try-error"){
  		AccuracyR[i, 2] <- a
  		end_time <- Sys.time()
  		RuntimeR[i, 2] <- end_time - start_time
  	} else {
  		AccuracyER[i, 2] <- 1
  	}
  # } else {
  	# AccuracyER[i, 2] <- 1
  # }

  	


  ## RF model
  start_time <- Sys.time()
  rf.tr_cl <- try(randomForest(x = train2[,c(8:20)], 
            y = as.factor(train2$ClusterR), ntree = 1000))

  if (class(rf.tr_cl) != "try-error") {
  # Predicting the Test set results
  trg = predict(rf.tr_cl, newdata = test2[,c(8: 20)])
  
  # Confusion Matrix
  cfm11[[i]] = table(Predict = trg, Actual = test2$ClusterR)
  a <- try(tr(cfm11[[i]])/sum(cfm11[[i]]))
  	if (class(a) != "try-error"){
  		AccuracyR[i, 3] <- a
  		end_time <- Sys.time()
  		RuntimeR[i, 3] <- end_time - start_time
  	} else {
  		AccuracyER[i, 3] <- 1
  	}
  } else {
  	AccuracyER[i, 3] <- 1
  }

 


  ## SOM model
  start_time <- Sys.time()  
  testXY <- list(independent = test1, dependent = factor(testing$ClusterR))

  train_som <- try(xyf(train1, classvec2classmat(factor(training$ClusterR)), SOM_grid, rlen = 1000, radius = 2.5, 
                      keep.data = TRUE, dist.fcts = "euclidean"))

  if (class(train_som) != "try-error") {
  # testing the model
  test_cluster <- predict(train_som, newdata = as.matrix(test1), whatmap = 1)
  cfm12[[i]] <- try(table(Predict = test_cluster$predictions[[2]], Actual = testing$ClusterR))
  if (class(cfm12[[i]]) != "try-error"){
  	a <- try(tr(cfm12[[i]])/sum(cfm12[[i]]))
  	if (class(a) != "try-error"){
  		AccuracyR[i, 4] <- a
  		end_time <- Sys.time()
  		RuntimeR[i, 4] <- end_time - start_time
  	} else {
  		AccuracyER[i, 4] <- 1
  	}
  } else {
  	AccuracyER[i, 4] <- 1
  }

  } else {
  	AccuracyER[i, 4] <- 1
  }


	


### checking the validity of ClusterS with other models

  ## base model
  start_time <- Sys.time()
  base.tr <- multinom(ClusterS ~ Moisture + pH + Salinity + Sand + Clay + Silt + Carbon +
                       Phosphorous +  Nitrogen + Sulphur + Sodium + Calcium + Potassium, train2, trace = FALSE)

  # if (class(base.tr) != "try-error") {
  p.tst7 <- predict(base.tr, test2)
  cfm13[[i]] <- table(Predict = p.tst7, Actual = test2$ClusterS)
  a <- try(tr(cfm13[[i]])/sum(cfm13[[i]]))
  	if (class(a) != "try-error"){
  		AccuracyS[i, 1] <- a
  		end_time <- Sys.time()
  		RuntimeS[i, 1] <- end_time - start_time
  	} else {
  		AccuracyES[i, 1] <- 1
  	}
  # } else {
  	# AccuracyES[i, 1] <- 1
  # }
 


  ## PC model
  start_time <- Sys.time()
  pc.tr <- prcomp(train2[,c(8:20)])

  # obtain PC values for all the rows and dependent variable to the data
  trg <- predict(pc.tr, train2) # predict() in PCA: predict the projection of new rows with PCA.
  trg <- as.data.frame(cbind(trg, train2$ClusterS)) # adding the response variable to data
  colnames(trg)[14] <- "ClusterS"
  tst <- predict(pc.tr, test2)
  tst <- as.data.frame(cbind(tst, test2$ClusterS))
  colnames(tst)[14] <- "ClusterS"

  model1 <- multinom(ClusterS ~ PC1+PC2, data = trg)

  # if (class(model1) != "try-error"){
  # Confusion matrix and misclassification error - testing data
  p.tst8 <- predict(model1, tst)
  cfm14[[i]] <- table(Predict = p.tst8, Actual = tst$ClusterS) # confusion matrix
  a <- try(tr(cfm14[[i]])/sum(cfm14[[i]]))
  	if (class(a) != "try-error"){
  		AccuracyS[i, 2] <- a
  		end_time <- Sys.time()
  		RuntimeS[i, 2] <- end_time - start_time
  	} else {
  		AccuracyES[i, 2] <- 1
  	}
  # } else {
  	# AccuracyES[i, 2] <- 1
  # }


  	


  ## RF model
  start_time <- Sys.time()
  rf.tr_cl <- try(randomForest(x = train2[,c(8:20)], 
            y = as.factor(train2$ClusterS), ntree = 1000))

  if (class(rf.tr_cl) != "try-error"){
  # Predicting the Test set results
  trg = predict(rf.tr_cl, newdata = test2[,c(8: 20)])
  
  # Confusion Matrix
  cfm15[[i]] = table(Predict = trg, Actual = test2$ClusterS)
  a <- try(tr(cfm15[[i]])/sum(cfm15[[i]]))
  	if (class(a) != "try-error"){
  		AccuracyS[i, 3] <- a
  		end_time <- Sys.time()
  		RuntimeS[i, 3] <- end_time - start_time
  	} else {
  		AccuracyES[i, 3] <- 1
  	}
  } else {
  	AccuracyES[i, 3] <- 1
  }





  ## SOM model
  start_time <- Sys.time()  
  testXY <- list(independent = test1, dependent = factor(testing$ClusterS))

  train_som <- try(xyf(train1, classvec2classmat(factor(training$ClusterS)), SOM_grid, rlen = 1000, radius = 2.5, 
                      keep.data = TRUE, dist.fcts = "euclidean"))

  if (class(train_som) != "try-error") {

  # testing the model
  test_cluster <- predict(train_som, newdata = as.matrix(test1), whatmap = 1)
  cfm16[[i]] <- try(table(Predict = test_cluster$predictions[[2]], Actual = testing$ClusterS))
  if (class(cfm16[[i]]) != "try-error"){
  	a <- try(tr(cfm16[[i]])/sum(cfm16[[i]]))
  	if (class(a) != "try-error"){
  		AccuracyS[i, 4] <- a
  		end_time <- Sys.time()
  		RuntimeS[i, 4] <- end_time - start_time
  	} else {
  		AccuracyES[i, 4] <- 1
  	}
  	} else {
  		AccuracyES[i, 4] <- 1
  	}

  } else {
  	AccuracyES[i, 4] <- 1
  }
  
} # end of for - prediction of clusters 


# note: change the names, for example, 'Two plots' to the specific name of the data for easy identification
write.csv(AccuracyB, paste0("Two_plots_", train_test[1], "_Accuracy_base.csv"))
write.csv(AccuracyP, paste0("Two_plots_", train_test[1], "_Accuracy_PC.csv"))
write.csv(AccuracyR, paste0("Two_plots_", train_test[1], "_Accuracy_RF.csv"))
write.csv(AccuracyS, paste0("Two_plots_", train_test[1], "_Accuracy_SOM.csv"))


# Accuracy_base <- apply(AccuracyB, 2, mean, na.rm = TRUE)
# Accuracy_PC <- apply(AccuracyP, 2, mean, na.rm = TRUE)
# Accuracy_RF <- apply(AccuracyR, 2, mean, na.rm = TRUE)
# Accuracy_SOM <- apply(AccuracyS, 2, mean, na.rm = TRUE)

base_miss <- apply(AccuracyE, 2, sum, na.rm = TRUE)
PC_miss <- apply(AccuracyEP, 2, sum, na.rm = TRUE)
RF_miss <- apply(AccuracyER, 2, sum, na.rm = TRUE)
SOM_miss <- apply(AccuracyES, 2, sum, na.rm = TRUE)

# write.csv(AccuracyE, paste0("Three_plots_", train_test[1], "_missed_base.csv"))
# write.csv(AccuracyEP, paste0("Three_plots_", train_test[1], "_missed_PC.csv"))
# write.csv(AccuracyER, paste0("Three_plots_", train_test[1], "_missed_RF.csv"))
# write.csv(AccuracyES, paste0("Three_plots_", train_test[1], "_missed_SOM.csv"))

print("end")

Missed_iter <- rbind(base_miss, PC_miss, RF_miss, SOM_miss)
write.csv(Missed_iter, paste0("Two_plots_", train_test[1], "_missed_iterations.csv"))

Runtime_base <- apply(RuntimeB, 2, mean, na.rm = TRUE)
Runtime_PC <- apply(RuntimeP, 2, mean, na.rm = TRUE)
Runtime_RF <- apply(RuntimeR, 2, mean, na.rm = TRUE)
Runtime_SOM <- apply(RuntimeS, 2, mean, na.rm = TRUE)

Runtime <- rbind(Runtime_base, Runtime_PC, Runtime_RF, Runtime_SOM)
write.csv(Runtime, paste0("Two_plots_", train_test[1], "_runtime.csv"))

} # end of simulation function

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~







# names of datasets 
dataC_plot # by plot

dataC_plot3 # by three plots per block

dataC_plot4 # by two plots per block

dataC_segment # full data



dataP.grid <- somgrid(xdim = 10, ydim = 10, topo = "hexagonal", toroidal = T, neighbourhood.fct =  "gaussian") # for full data    
dataP.grid <- somgrid(xdim = 6, ydim = 6, topo = "hexagonal", toroidal = T, neighbourhood.fct =  "gaussian") # for by plot data
dataP.grid <- somgrid(xdim = 4, ydim = 4, topo = "hexagonal", toroidal = T, neighbourhood.fct =  "gaussian") # for by three plots per block data
dataP.grid <- somgrid(xdim = 2, ydim = 2, topo = "hexagonal", toroidal = T, neighbourhood.fct =  "gaussian") # for by two plo



Soil_sim_small <- function(data, train_test, SOM_grid){


cfm1 <- cfm2 <- cfm3 <- cfm4 <- cfm5 <- cfm6 <- cfm7 <- cfm8 <- cfm9 <- cfm10 <- cfm11 <- cfm12 <- list()
 cfm13 <- cfm14 <- cfm15 <- cfm16 <- list()
# cfm <- list() # can be used which will be recycled through all 16 analysis

AccuracyB <- AccuracyP <-  AccuracyR <- AccuracyS  <- AccuracyE <- AccuracyEP <- AccuracyER <- AccuracyES <- matrix(NA, 100, 4) 
colnames(AccuracyB) <- colnames(AccuracyP) <- colnames(AccuracyR) <- colnames(AccuracyS) <- c("Base", "PC", "RF", "SOM")
colnames(AccuracyE) <- colnames(AccuracyEP) <- colnames(AccuracyER) <- colnames(AccuracyES) <- c("Base", "PC", "RF", "SOM")

RuntimeB <- RuntimeP <- RuntimeR <- RuntimeS <- matrix(NA, 100, 4)
colnames(RuntimeB) <- colnames(RuntimeP) <- colnames(RuntimeR) <- colnames(RuntimeS) <- c("Base", "PC", "RF", "SOM")


for (i in 1:100){
  ind <- sample(2, nrow(data),
          replace = TRUE,
          prob = train_test)
  training <- data[ind==1,]
  testing <- data[ind==2,]

  # standardizing the data
  train1 <- as.matrix(scale(training[,c(7, 10:11, 14:16)]))
  test1 <- scale(testing[,c(7, 10:11, 14:16)], center = attr(train1, "scaled:center")) 

  train2 <- cbind(training[,c(1:3,17:20)], train1)
  test2 <- cbind(testing[,c(1:3,17:20)], test1)

  ### checking the validity of ClusterB with other models

  ## base model
  start_time <- Sys.time()
  base.tr <- multinom(ClusterB ~ Sand + Carbon + Phosphorous + Sodium + Calcium + Potassium, train2, trace = FALSE)
 
  # if (class(base.tr) != "try-error"){
  p.tst1 <- predict(base.tr, test2)
  cfm1[[i]] <- table(Predict = p.tst1, Actual = test2$ClusterB)
  a <- try(tr(cfm1[[i]])/sum(cfm1[[i]]))
  	if (class(a) != "try-error"){
  		AccuracyB[i, 1] <- a
  		end_time <- Sys.time()
  		RuntimeB[i, 1] <- end_time - start_time
  	} else {
  		AccuracyE[i, 1] <- 1
  	}
  # } else {
  	# AccuracyE[i, 1] <- 1
  # }
  



  ## PC model
  start_time <- Sys.time()

  pc.tr <- prcomp(train2[,c(8:13)]) 
  # obtain PC values for all the rows and dependent variable to the data
  trg <- predict(pc.tr, train2) # predict() in PCA: predict the projection of new rows with PCA.
  trg <- as.data.frame(cbind(trg, train2$ClusterB)) # adding the response variable to data
  colnames(trg)[7] <- "ClusterB"
  tst <- predict(pc.tr, test2)
  tst <- as.data.frame(cbind(tst, test2$ClusterB))
  colnames(tst)[7] <- "ClusterB"

  model1 <- multinom(ClusterB ~ PC1+PC2, data = trg, trace = FALSE)

  # if (class(model1) != "try-error"){ 
  # Confusion matrix and misclassification error - testing data
  p.tst2 <- predict(model1, tst)
  cfm2[[i]] <- table(Predict = p.tst2, Actual = tst$ClusterB) # confusion matrix
  a <- try(tr(cfm2[[i]])/sum(cfm2[[i]]))
  	if (class(a) != "try-error"){
  		AccuracyB[i, 2] <- a
  		end_time <- Sys.time()
  		RuntimeB[i, 2] <- end_time - start_time
  	} else {
  		AccuracyE[i, 2] <- 1
  	}
  # } else {
  	# AccuracyE[i, 2] <- 1
  # }

  
  	


  ## RF model
  start_time <- Sys.time()
  rf.tr_cl <- try(randomForest(x = train2[,c(8:13)], 
            y = as.factor(train2$ClusterB), ntree = 1000))

if (class(rf.tr_cl) != "try-error") {
  # Predicting the Test set results
  trg = predict(rf.tr_cl, newdata = test2[,c(8: 13)])
  
  # Confusion Matrix
  cfm3[[i]] = table(Predict = trg, Actual = test2$ClusterB)
  a <- try(tr(cfm3[[i]])/sum(cfm3[[i]]))
  	if (class(a) != "try-error"){
  		AccuracyB[i, 3] <- a
  		end_time <- Sys.time()
  		RuntimeB[i, 3] <- end_time - start_time
  	} else {
  		AccuracyE[i, 3] <- 1
  	}
} else {
	# AccuracyE[i, 3] <- 1
}




  ## SOM model
  start_time <- Sys.time()  
  testXY <- list(independent = test1, dependent = factor(testing$ClusterB))

  train_som <- try(xyf(train1, classvec2classmat(factor(training$ClusterB)), SOM_grid, rlen = 1000, radius = 2.5, 
                      keep.data = TRUE, dist.fcts = "euclidean"))

  if (class(train_som) != "try-error"){
  # testing the model
  test_cluster <- predict(train_som, newdata = as.matrix(test1), whatmap = 1)
  cfm4[[i]] <- try(table(Predict = test_cluster$predictions[[2]], Actual = testing$ClusterB))
  if (class(cfm4[[i]]) != "try-error"){
  	a <- try(tr(cfm4[[i]])/sum(cfm4[[i]]))
  	if (class(a) != "try-error"){
  		AccuracyB[i, 4] <- a
  		end_time <- Sys.time()
  		RuntimeB[i, 4] <- end_time - start_time
  	} else {
  		AccuracyE[i, 4] <- 1
  	}
  } else {
  	AccuracyE[i, 4] <- 1
  }

  } else {
  	AccuracyE[i, 4] <- 1
  }
  
  	



### checking the validity of ClusterP with other models

  ## base model
  start_time <- Sys.time()
  base.tr <- multinom(ClusterP ~ Sand + Carbon + Phosphorous + Sodium + Calcium + Potassium, train2, trace = FALSE)
                       
# if(class(base.tr) != "try-error"){
  p.tst3 <- predict(base.tr, test2)
  cfm5[[i]] <- table(Predict = p.tst3, Actual = test2$ClusterP)
  a <- try(tr(cfm5[[i]])/sum(cfm5[[i]]))
  	if (class(a) != "try-error"){
  		AccuracyP[i, 1] <- a
  		end_time <- Sys.time()
  		RuntimeP[i, 1] <- end_time - start_time
  	} else {
  		AccuracyEP[i, 1] <- 1
  	}
# } else {
   # AccuracyEP[i, 1] <- 1
# }
   


  ## PC model
  start_time <- Sys.time()
  pc.tr <- prcomp(train2[,c(8:13)])

  # obtain PC values for all the rows and dependent variable to the data
  trg <- predict(pc.tr, train2) # predict() in PCA: predict the projection of new rows with PCA.
  trg <- as.data.frame(cbind(trg, train2$ClusterP)) # adding the response variable to data
  colnames(trg)[7] <- "ClusterP"
  tst <- predict(pc.tr, test2)
  tst <- as.data.frame(cbind(tst, test2$ClusterP))
  colnames(tst)[7] <- "ClusterP"

  model1 <- multinom(ClusterP ~ PC1+PC2, data = trg, trace = FALSE)

  # if (class(model1) != "try-error"){
  p.tst4 <- predict(model1, tst)
  # Confusion matrix and misclassification error - testing data
  cfm6[[i]] <- table(Predict = p.tst4, Actual = tst$ClusterP) # confusion matrix
  a <- try(tr(cfm6[[i]])/sum(cfm6[[i]]))
  	if (class(a) != "try-error"){
  		AccuracyP[i, 2] <- a
  		end_time <- Sys.time()
  		RuntimeP[i, 2] <- end_time - start_time
  	} else {
  		AccuracyEP[i, 2] <- 1
  	}	
  # } else {
  	AccuracyEP[i, 2] <- 1
  # }

 
 


  ## RF model
  start_time <- Sys.time()
  rf.tr_cl <- try(randomForest(x = train2[,c(8:13)], 
            y = as.factor(train2$ClusterP), ntree = 1000))

  if (class(rf.tr_cl) != "try-error") {
  # Predicting the Test set results
  trg = predict(rf.tr_cl, newdata = test2[,c(8: 13)])
  
  # Confusion Matrix
  cfm7[[i]] = table(Predict = trg, Actual = test2$ClusterP)
  a <- try(tr(cfm7[[i]])/sum(cfm7[[i]]))
  	if (class(a) != "try-error"){
  		AccuracyP[i, 3] <- a
  		end_time <- Sys.time()
  		RuntimeP[i, 3] <- end_time - start_time
  	} else {
  		AccuracyEP[i, 3] <- 1
  	}
  } else {
  	AccuracyEP[i, 3] <- 1
  }
  


  ## SOM model
  start_time <- Sys.time()  
  testXY <- list(independent = test1, dependent = factor(testing$ClusterP))

  train_som <- try(xyf(train1, classvec2classmat(factor(training$ClusterP)), SOM_grid, rlen = 1000, radius = 2.5, 
                      keep.data = TRUE, dist.fcts = "euclidean"))

 if (class(train_som) != "try-error") {
  # testing the model
  test_cluster <- predict(train_som, newdata = as.matrix(test1), whatmap = 1)
  cfm8[[i]] <- try(table(Predict = test_cluster$predictions[[2]], Actual = testing$ClusterP))
  if (class(cfm8[[i]]) != "try-error"){
  	a <- try(tr(cfm8[[i]])/sum(cfm8[[i]]))
  	if (class(a) != "try-error"){
  		AccuracyP[i, 4] <- a
  		end_time <- Sys.time()
  		RuntimeP[i, 4] <- end_time - start_time
  	} else {
  		AccuracyEP[i, 4] <- 1
  	}
  } else {
  	AccuracyEP[i, 4] <- 1
  }
 } else {
  	AccuracyEP[i, 4] <- 1
  }

  
  


### checking the validity of ClusterR with other models

  ## base model
  start_time <- Sys.time()
  base.tr <- multinom(ClusterR ~ Sand + Carbon + Phosphorous + Sodium + Calcium + Potassium, train2, trace = FALSE)
                     

  # if (class(base.tr) != "try-error") {
  p.tst5 <- predict(base.tr, test2)
  cfm9[[i]] <- table(Predict = p.tst5, Actual = test2$ClusterR)
  a <- try(tr(cfm9[[i]])/sum(cfm9[[i]]))
  	if (class(a) != "try-error"){
  		AccuracyR[i, 1] <- a
  		end_time <- Sys.time()
  		RuntimeR[i, 1] <- end_time - start_time
  	} else {
  		AccuracyER[i, 1] <- 1
  	}
  # } else {
  	# AccuracyER[i, 1] <- 1
  # }
 



  ## PC model
  start_time <- Sys.time()
  pc.tr <- prcomp(train2[,c(8:13)])

  # obtain PC values for all the rows and dependent variable to the data
  trg <- predict(pc.tr, train2) # predict() in PCA: predict the projection of new rows with PCA.
  trg <- as.data.frame(cbind(trg, train2$ClusterR)) # adding the response variable to data
  colnames(trg)[7] <- "ClusterR"
  tst <- predict(pc.tr, test2)
  tst <- as.data.frame(cbind(tst, test2$ClusterR))
  colnames(tst)[7] <- "ClusterR"

  model1 <- multinom(ClusterR ~ PC1+PC2, data = trg, trace = FALSE)

  # if (class(model1) != "try-error") {

  # Confusion matrix and misclassification error - testing data
  p.tst6 <- predict(model1, tst)
  cfm10[[i]] <- table(Predict = p.tst6, Actual = tst$ClusterR) # confusion matrix
  a <- try(tr(cfm10[[i]])/sum(cfm10[[i]]))
  	if (class(a) != "try-error"){
  		AccuracyR[i, 2] <- a
  		end_time <- Sys.time()
  		RuntimeR[i, 2] <- end_time - start_time
  	} else {
  		AccuracyER[i, 2] <- 1
  	}
  # } else {
  	# AccuracyER[i, 2] <- 1
  # }

  	


  ## RF model
  start_time <- Sys.time()
  rf.tr_cl <- try(randomForest(x = train2[,c(8:13)], 
            y = as.factor(train2$ClusterR), ntree = 1000))

  if (class(rf.tr_cl) != "try-error") {
  # Predicting the Test set results
  trg = predict(rf.tr_cl, newdata = test2[,c(8: 13)])
  
  # Confusion Matrix
  cfm11[[i]] = table(Predict = trg, Actual = test2$ClusterR)
  a <- try(tr(cfm11[[i]])/sum(cfm11[[i]]))
  	if (class(a) != "try-error"){
  		AccuracyR[i, 3] <- a
  		end_time <- Sys.time()
  		RuntimeR[i, 3] <- end_time - start_time
  	} else {
  		AccuracyER[i, 3] <- 1
  	}
  } else {
  	AccuracyER[i, 3] <- 1
  }

 


  ## SOM model
  start_time <- Sys.time()  
  testXY <- list(independent = test1, dependent = factor(testing$ClusterR))

  train_som <- try(xyf(train1, classvec2classmat(factor(training$ClusterR)), SOM_grid, rlen = 1000, radius = 2.5, 
                      keep.data = TRUE, dist.fcts = "euclidean"))

  if (class(train_som) != "try-error") {
  # testing the model
  test_cluster <- predict(train_som, newdata = as.matrix(test1), whatmap = 1)
  cfm12[[i]] <- try(table(Predict = test_cluster$predictions[[2]], Actual = testing$ClusterR))
  if (class(cfm12[[i]]) != "try-error"){
  	a <- try(tr(cfm12[[i]])/sum(cfm12[[i]]))
  	if (class(a) != "try-error"){
  		AccuracyR[i, 4] <- a
  		end_time <- Sys.time()
  		RuntimeR[i, 4] <- end_time - start_time
  	} else {
  		AccuracyER[i, 4] <- 1
  	}
  } else {
  	AccuracyER[i, 4] <- 1
  }

  } else {
  	AccuracyER[i, 4] <- 1
  }





### checking the validity of ClusterS with other models

  ## base model
  start_time <- Sys.time()
  base.tr <- multinom(ClusterS ~ Sand + Carbon + Phosphorous + Sodium + Calcium + Potassium, train2, trace = FALSE)

  # if (class(base.tr) != "try-error") {
  p.tst7 <- predict(base.tr, test2)
  cfm13[[i]] <- table(Predict = p.tst7, Actual = test2$ClusterS)
  a <- try(tr(cfm13[[i]])/sum(cfm13[[i]]))
  	if (class(a) != "try-error"){
  		AccuracyS[i, 1] <- a
  		end_time <- Sys.time()
  		RuntimeS[i, 1] <- end_time - start_time
  	} else {
  		AccuracyES[i, 1] <- 1
  	}
  # } else {
  	# AccuracyES[i, 1] <- 1
  # }
 


  ## PC model
  start_time <- Sys.time()
  pc.tr <- prcomp(train2[,c(8:13)])

  # obtain PC values for all the rows and dependent variable to the data
  trg <- predict(pc.tr, train2) # predict() in PCA: predict the projection of new rows with PCA.
  trg <- as.data.frame(cbind(trg, train2$ClusterS)) # adding the response variable to data
  colnames(trg)[7] <- "ClusterS"
  tst <- predict(pc.tr, test2)
  tst <- as.data.frame(cbind(tst, test2$ClusterS))
  colnames(tst)[7] <- "ClusterS"

  model1 <- multinom(ClusterS ~ PC1+PC2, data = trg, trace = FALSE)

  # if (class(model1) != "try-error"){
  # Confusion matrix and misclassification error - testing data
  p.tst8 <- predict(model1, tst)
  cfm14[[i]] <- table(Predict = p.tst8, Actual = tst$ClusterS) # confusion matrix
  a <- try(tr(cfm14[[i]])/sum(cfm14[[i]]))
  	if (class(a) != "try-error"){
  		AccuracyS[i, 2] <- a
  		end_time <- Sys.time()
  		RuntimeS[i, 2] <- end_time - start_time
  	} else {
  		AccuracyES[i, 2] <- 1
  	}
  # } else {
  	# AccuracyES[i, 2] <- 1
  # }


  	


  ## RF model
  start_time <- Sys.time()
  rf.tr_cl <- try(randomForest(x = train2[,c(8:13)], 
            y = as.factor(train2$ClusterS), ntree = 1000))

  if (class(rf.tr_cl) != "try-error"){
  # Predicting the Test set results
  trg = predict(rf.tr_cl, newdata = test2[,c(8: 13)])
  
  # Confusion Matrix
  cfm15[[i]] = table(Predict = trg, Actual = test2$ClusterS)
  a <- try(tr(cfm15[[i]])/sum(cfm15[[i]]))
  	if (class(a) != "try-error"){
  		AccuracyS[i, 3] <- a
  		end_time <- Sys.time()
  		RuntimeS[i, 3] <- end_time - start_time
  	} else {
  		AccuracyES[i, 3] <- 1
  	}
  } else {
  	AccuracyES[i, 3] <- 1
  }





  ## SOM model
  start_time <- Sys.time()  
  testXY <- list(independent = test1, dependent = factor(testing$ClusterS))

  train_som <- try(xyf(train1, classvec2classmat(factor(training$ClusterS)), SOM_grid, rlen = 1000, radius = 2.5, 
                      keep.data = TRUE, dist.fcts = "euclidean"))

  if (class(train_som) != "try-error") {

  # testing the model
  test_cluster <- predict(train_som, newdata = as.matrix(test1), whatmap = 1)
  cfm16[[i]] <- try(table(Predict = test_cluster$predictions[[2]], Actual = testing$ClusterS))
  if (class(cfm16[[i]]) != "try-error"){
  	a <- try(tr(cfm16[[i]])/sum(cfm16[[i]]))
  	if (class(a) != "try-error"){
  		AccuracyS[i, 4] <- a
  		end_time <- Sys.time()
  		RuntimeS[i, 4] <- end_time - start_time
  	} else {
  		AccuracyES[i, 4] <- 1
  	}
  	} else {
  		AccuracyES[i, 4] <- 1
  	}

  } else {
  	AccuracyES[i, 4] <- 1
  }

	
  
} # end of for - prediction of clusters 



write.csv(AccuracyB, paste0("Two_plots_", train_test[1], "_Accuracy_base_reduced.csv"))
write.csv(AccuracyP, paste0("Two_plots_", train_test[1], "_Accuracy_PC_reduced.csv"))
write.csv(AccuracyR, paste0("Two_plots_", train_test[1], "_Accuracy_RF_reduced.csv"))
write.csv(AccuracyS, paste0("Two_plots_", train_test[1], "_Accuracy_SOM_reduced.csv"))


base_miss <- apply(AccuracyE, 2, sum, na.rm = TRUE)
PC_miss <- apply(AccuracyEP, 2, sum, na.rm = TRUE)
RF_miss <- apply(AccuracyER, 2, sum, na.rm = TRUE)
SOM_miss <- apply(AccuracyES, 2, sum, na.rm = TRUE)

print("end")

Missed_iter <- rbind(base_miss, PC_miss, RF_miss, SOM_miss)
write.csv(Missed_iter, paste0("Two_plots_", train_test[1], "_missed_iterations_reduced.csv"))

Runtime_base <- apply(RuntimeB, 2, mean, na.rm = TRUE)
Runtime_PC <- apply(RuntimeP, 2, mean, na.rm = TRUE)
Runtime_RF <- apply(RuntimeR, 2, mean, na.rm = TRUE)
Runtime_SOM <- apply(RuntimeS, 2, mean, na.rm = TRUE)

Runtime <- rbind(Runtime_base, Runtime_PC, Runtime_RF, Runtime_SOM)
write.csv(Runtime, paste0("Two_plots_", train_test[1], "_runtime_reduced.csv"))

} # end of simulation function



~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~






### visualization ###

# example

Tp <- read.csv("Three_plots_2_Accuracy_SOM_reduced.csv") 
# Tp <- melt(Tp[,2:5])
Tp <- melt(Tp[,2:4]) # SOM is all NAs
colnames(Tp) <- c("Model", "Accuracy")

Tp$Cluster <- "SOM"
Tp$Area <- "9x"
Tp$Train.prop <- "0.95"

write.csv(Tp, "Three_plots_0.95_SOM_melt_reduced.csv")



a <- read.csv("Two_plot_0.5_base_reduced_melt.csv")
b <- read.csv("Two_plot_0.5_PC_reduced_melt.csv")
d <- read.csv("Two_plot_0.5_RF_reduced_melt.csv")
e <- read.csv("Two_plot_0.5_SOM_reduced_melt.csv")

g <- read.csv("Two_plot_0.6_base_reduced_melt.csv")
h <- read.csv("Two_plot_0.6_PC_reduced_melt.csv")
i <- read.csv("Two_plot_0.6_RF_reduced_melt.csv")
j <- read.csv("Two_plot_0.6_SOM_reduced_melt.csv")

k <- read.csv("Two_plot_0.7_base_reduced_melt.csv")
l <- read.csv("Two_plot_0.7_PC_reduced_melt.csv")
m <- read.csv("Two_plot_0.7_RF_reduced_melt.csv")
n <- read.csv("Two_plot_0.7_SOM_reduced_melt.csv")

o <- read.csv("Two_plot_0.8_base_reduced_melt.csv")
p <- read.csv("Two_plot_0.8_PC_reduced_melt.csv")
q <- read.csv("Two_plot_0.8_RF_reduced_melt.csv")
r <- read.csv("Two_plot_0.8_SOM_reduced_melt.csv")

# s <- read.csv("Two_plot_0.85_base_reduced_melt.csv")
# t <- read.csv("Two_plot_0.85_PC_reduced_melt.csv")
# u <- read.csv("Two_plot_0.85_RF_reduced_melt.csv")
# v <- read.csv("Two_plot_0.85_SOM_reduced_melt.csv")

s <- read.csv("Two_plot_0.9_base_melt_reduced.csv")
t <- read.csv("Two_plot_0.9_PC_melt_reduced.csv")
u <- read.csv("Two_plot_0.9_RF_melt_reduced.csv")
v <- read.csv("Two_plot_0.9_SOM_melt_reduced.csv")

w <- read.csv("Two_plot_0.95_base_melt_reduced.csv")
x <- read.csv("Two_plot_0.95_PC_melt_reduced.csv")
y <- read.csv("Two_plot_0.95_RF_melt_reduced.csv")
z <- read.csv("Two_plot_0.95_SOM_melt_reduced.csv")

Plot_red <- rbind(a, b, d, e, g, h, i, j, k, l, m, n, 
	o, p, q, r, s, t, u, v, w, x, y, z)

Plot_red <- rbind(a, b, d, e, g, h, i, j, k, l, m, n, 
	o, p, q, r)

Plot_red <- rbind(s, t, u, v, w, x, y, z)


write.csv(Plot_red, "Two_plot_reduced_melt.csv")

P1 <- read.csv("Plot_reduced_melt.csv")
P2 <- read.csv("Two_plot_reduced_melt.csv")
P3 <- read.csv("Three_plot_reduced_melt.csv")
P4 <- read.csv("Segment_reduced_melt.csv")

P1234 <- rbind(P1, P2, P3, P4)



P1234avg <- P1234 %>% group_by(Train.prop, Area, Cluster, Model) %>%
            summarise(across(c(Accuracy), mean, na.rm = TRUE),
            .groups = 'drop') %>%
            as.data.frame()

# P14avg <- P14 %>% group_by(Train.prop, Area, Cluster, Model) %>%
#             summarise(across(c(Accuracy), mean, na.rm = TRUE),
#             .groups = 'drop') %>%
#             as.data.frame()

# P14avg$Train.prop <- as.factor(P14avg$Train.prop)
P1234avg$Train.prop <- as.factor(P1234avg$Train.prop)           


ggplot(P1234avg, aes(x = Train.prop, y = Accuracy)) +
	geom_point(aes(color = Model)) + 
	facet_grid(Cluster ~ Area) + 
	ylab("Classification accuracy") +
    xlab("Proportion of training population") + 
    theme_bw() +
    theme(axis.text.x = element_text(hjust = 1, angle =45), 
    axis.title = element_text(size =14, face = "bold"), plot.title = element_text(face = "bold")
        )
