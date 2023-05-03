### Multi-environmental soil data modeling###
### Ani A. Elas 	December 2022 ###

# soil data are evaluated using different supervised models, clustered using unsupervised models, and predictabiltiy tested 


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

library(gtools) # for combinations of variables for simulation
library(dplyr) # for grouping and applying functions across groups






setwd () 

data1 <- read.csv("XXX.csv", h=T, sep=",") # provide the name of the dataset

data3 <- data1[data1$Location == "Pune",]

# initial visualization example
ggplot(data3, aes(x=Plot, y=Segment, z= Potassium)) +
	geom_contour(aes(colour = after_stat(level)), size = 1) +
	theme_bw() +
	xlab("Plot") + 
	ylab("Segment") +
	facet_wrap(~Block) +
	theme(axis.title = element_text(size =14, face = "bold"), plot.title = element_text(face = "bold")) +
	  labs(fill="") +
	  # ggtitle("Observed electrical conductivity(μS/cm) from blocks A-D \n in Nagpur field")
	  # ggtitle("Gravimetric moisture content (%) from \n blocks A-D in Delhi field")
    ggtitle("P-xiii")
	  # ggtitle("Observed pH from blocks A-D in Nagpur field")

# multiple histograms
	ggplot(data = melt(data1), aes(x=value))+
  # geom_histogram(color="black", fill="white", scales = "free")+
  geom_histogram(scales = "free")+
  facet_grid(variable ~ .) +
  ggtitle("Histogram of soil variables") +
  theme_bw()


# correlation graph
  # find multicollinearity and choose variables 

  data1$Location <- as.factor(data1$Location)

  pairs.panels(data1[,-c(1:8)],
             gap = 0,
             bg = c( "black", "blue", "green", "red")[data1$Location],
             pch = 21 + as.numeric(data1$Location),
             main="Correlation and variability of soil features among locations", cex.main = 0.9)

  legend("topright", inset = c(-0.01, 1.2), horiz=TRUE, bty = "n",
	 fill = unique(data1$Location),
   legend = c( levels(data1$Location)))


  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# multinomial regression with normalized data - base model

data1b <- na.omit(data1)
data1a <- cbind(data1b[,1:8], scale(data1b[,-c(1:8, 12)]))


base <- multinom(Location ~ Moisture + pH + Salinity + Sand + Clay + Silt + Carbon +
                       Phosphorous +  Nitrogen + Sulphur + Sodium + Calcium + Potassium, data1a)
summary(base)


# prediction using base - all variables

# Location as the response variable
cfm <- list()
Accuracy <- vector()
for (i in 1:100){

# partitioning the data 
ind <- sample(2, nrow(data1a), replace = TRUE, prob = c(0.7, 0.3))
training <- data1a[ind==1,]
testing <- data1a[ind==2,]

# training model
base.tr <- multinom(Location ~ Moisture + pH + Salinity + Sand + Clay + Silt + Carbon +
                       Phosphorous +  Nitrogen + Sulphur + Sodium + Calcium + Potassium, training)


# Confusion matrix and misclassification error - testing data
p.tst <- predict(base.tr, testing)
cfm[[i]] <- table(Predict = p.tst, Actual = testing$Location) # confusion matrix
Accuracy[i] <- tr(cfm[[i]])/sum(cfm[[i]])
} # end of for - 70:30 with response variable as Location

 
cfm2 <- round(Reduce("+", cfm)/length(cfm),0) 
cfm2 <- Reduce("+", cfm)/length(cfm) 
mat <- melt(cfm2)
colnames(mat) <- c("pred", "Actual", "n")

plot_confusion_matrix(mat, target_col = "Actual", prediction_col = "pred", 
              counts_col = "n") +
ggtitle("Average confusion matrix obtained from 100 replications 
  of 70:30 training: test cross-validation with Location  
  as the response variable in multinomial base model") +

  ylab("Prediction") +
  xlab("Actual") +
  theme_bw()

mean(Accuracy)
  hist(Accuracy, main = "Accuracy in classification from 100 replications
             of 70:30 training: test cross-validation with Location  
             as the response variable in multinomial base model" )
  text(0.98, 30, "mean accuracy = 0.99")



# clustering using all variables - base model
set.seed(1)
# Clustering
fviz_nbclust(data1a[c(9:21)], cluster::pam, method = "wss") # determine the number of clusters


set.seed(1)
# do k-medoid clustering using PC1 and 2- kk means - medoids 
pam.base <- pam(data1a[c(9:21)], 5)


# add cluster to dataset
data_cluster <- data.frame(data1a, cluster = pam.base$cluster)

data_cluster$cluster <- as.factor(data_cluster$cluster)

ggplot(data_cluster, aes(x = Block, y = Plot, fill = cluster)) +
  facet_grid(. ~ Location)  + 
  geom_tile() +
  scale_fill_discrete() +
  ggtitle("Clustering of plot segments in various locations \n based on principal components") +
  ylab("Plot") +
  xlab("Blocks") +
  theme_bw()


# base with cluster - CV
cfm <- list()
cfm_miss <- list()
Accuracy <- vector()


for (i in 1:100){
  ind <- sample(2, nrow(data_cluster),
          replace = TRUE,
          prob = c(0.7, 0.3))
  training <- data_cluster[ind==1,]
  testing <- data_cluster[ind==2,]

  # training the model
  base.tr <- multinom(cluster ~ Moisture + pH + Salinity + Sand + Clay + Silt + Carbon +
                       Phosphorous +  Nitrogen + Sulphur + Sodium + Calcium + Potassium, training)

  p.tst <- predict(base.tr, testing)
  cfm[[i]] <- table(Predict = p.tst, Actual = testing$cluster)
  Accuracy[i] <- tr(cfm[[i]])/sum(cfm[[i]])
  
} # end of for - prediction of clusters using PC1 and 2

cfm2 <- round(Reduce("+", cfm)/length(cfm),0) 
mat <- melt(cfm2)
colnames(mat) <- c("pred", "Actual", "n")

plot_confusion_matrix(mat, target_col = "Actual", prediction_col = "pred", 
              counts_col = "n") +
ggtitle("Average confusion matrix obtained from 100 replications 
  of 70:30 training: test cross-validation with pc based cluster  
  as the response variable in multinomial base model") +
  ylab("Prediction") +
  xlab("Actual") +
  theme_bw()

mean(Accuracy)
  hist(Accuracy, main = "Accuracy in classification from 100 replications
             of 70:30 training: test cross-validation with pc based 10  
             clusters as response variable in multinomial model" )
  text(0.965, 30, "mean accuracy = 0.99")


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





# PCA - works in R 4.2.

pc <- prcomp(na.omit(data1[,-c(1:8, 12)]),
             center = TRUE,
            scale. = TRUE)
attributes(pc)

print(pc)

summary(pc)

# orthogonality of PCs
pairs.panels(pc$x,
             gap=0,
             bg = c("red", "yellow", "blue")[data1$Location],
             pch=21)

# Scree plot
ggscreeplot(pc) + 
	ggtitle ("Scree plot from PCA") +
	theme_bw() + 
	theme(plot.title = element_text(hjust = 0.5)) 
	  

# biplot
g <- ggbiplot(pc,
              obs.scale = 1,
              var.scale = 1,
              groups = na.omit(data1)$Location,
              ellipse = TRUE,
              # circle = TRUE,
              ellipse.prob = 0.68)
g <- g + scale_color_discrete(name = '') + theme_bw()
g <- g + ggtitle("Classification of soil based on features") +
  theme(plot.title = element_text(hjust = 0.5))
g <- g + theme(legend.direction = 'horizontal',
               legend.position = 'bottom')
print(g)

pc_predict <- predict(pc, data1)
pc_data1 <- data.frame(pc_predict, data1[1])
model_nm <- multinom(Location~PC1+PC2, data = pc_data1)
summary(model_nm)

# prediction using PCA - all variables

# Location as the response variable
cfm <- list()
Accuracy <- vector()
for (i in 1:100){

# partitioning the data 
ind <- sample(2, nrow(data1),
              replace = TRUE,
              prob = c(0.7, 0.3))
training <- data1[ind==1,]
testing <- data1[ind==2,]

# training model
pc.tr <- prcomp(na.omit(training[,-c(1:8, 12)]),
             center = TRUE,
            scale. = TRUE)

# obtain PC values for all the rows and dependent variable to the data
trg <- predict(pc.tr, training) # predict() in PCA: predict the projection of new rows with PCA.
trg <- data.frame(trg, training[1]) # adding the response variable to data
tst <- predict(pc.tr, testing)
tst <- data.frame(tst, testing[1])

# Because, dependent variable has 4 level, use multinomial logistic regression.
	# trg$Location <- relevel(trg$Location, ref = "Delhi")
model1 <- multinom(Location~PC1+PC2, data = trg)


# Confusion matrix and misclassification error - testing data
p.tst <- predict(model1, tst)
cfm[[i]] <- table(Predict = p.tst, Actual = tst$Location) # confusion matrix
Accuracy[i] <- tr(cfm[[i]])/sum(cfm[[i]])
} # end of for - 70:30 with response variable as Location


 
cfm2 <- round(Reduce("+", cfm)/length(cfm),0) 
cfm2 <- Reduce("+", cfm)/length(cfm) 
mat <- melt(cfm2)
colnames(mat) <- c("pred", "Actual", "n")

plot_confusion_matrix(mat, target_col = "Actual", prediction_col = "pred", 
	            counts_col = "n") +
ggtitle("Average confusion matrix obtained from 100 replications 
	of 70:30 training: test cross-validation with Location  
	as the response variable in multinomial model") +
  ylab("Prediction") +
  xlab("Actual") +
  theme_bw()

mean(Accuracy)
  hist(Accuracy, main = "Accuracy in classification from 100 replications
  	         of 70:30 training: test cross-validation with Location  
	           as the response variable in multinomial model" )
  text(0.852, 20, "mean accuracy = 0.92")

# Clustering based on PC1 and 2, then CV with response as cluster - all variables
set.seed(2)
# Clustering
fviz_nbclust(na.omit(pc_data1[c(1,2)]), cluster::pam, method = "wss") # determine the number of clusters


set.seed(3)
# do k-medoid clustering using PC1 and 2- kk means - medoids 
pam.mn <- pam(na.omit(pc_data1[c(1,2)]), 4)



# add cluster to dataset
data_cluster <- data.frame(na.omit(data1), cluster = pam.mn$cluster)

data_cluster$cluster <- as.factor(data_cluster$cluster)

ggplot(data_cluster, aes(x = Block, y = Plot, fill = cluster)) +
  facet_grid(. ~ Location)  + 
  geom_tile() +
  scale_fill_discrete() +
  ggtitle("Clustering of plot segments in various locations \n based on principal components") +
  ylab("Plot") +
  xlab("Blocks") +
  theme_bw()


cfm <- list()
cfm_miss <- list()
Accuracy <- vector()
pc_cluster <- data.frame(na.omit(pc_data1), cluster = pam.mn$cluster)
model2 <- multinom(cluster ~ PC1+PC2, data = pc_cluster)

for (i in 1:100){
	ind <- sample(2, nrow(pc_cluster),
          replace = TRUE,
          prob = c(0.7, 0.3))
	training <- pc_cluster[ind==1,]
	testing <- pc_cluster[ind==2,]

	# training the model
	train1 <- multinom(cluster ~ PC1+PC2, data = training)
	p.tst <- predict(train1, testing)
  cfm[[i]] <- table(Predict = p.tst, Actual = testing$cluster)
  Accuracy[i] <- tr(cfm[[i]])/sum(cfm[[i]])
	
} # end of for - prediction of clusters using PC1 and 2

cfm2 <- round(Reduce("+", cfm)/length(cfm),0) 
mat <- melt(cfm2)
colnames(mat) <- c("pred", "Actual", "n")

plot_confusion_matrix(mat, target_col = "Actual", prediction_col = "pred", 
	            counts_col = "n") +
ggtitle("Average confusion matrix obtained from 100 replications 
	of 70:30 training: test cross-validation with pc based cluster  
	as the response variable in multinomial model") +
  ylab("Prediction") +
  xlab("Actual") +
  theme_bw()

mean(Accuracy)
  hist(Accuracy, main = "Accuracy in classification from 100 replications
  	         of 70:30 training: test cross-validation with pc based 10  
	           clusters as response variable in multinomial model" )
  text(0.94, 25, "mean accuracy = 0.98")

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Random Forest (RF)


data1a <- na.omit(data1)
# Fitting Random Forest to the whole dataset
set.seed(4)  # Setting seed
model_rf = randomForest(x = data1a[,-c(1:8, 12)], y = as.factor(data1a$Location), ntree = 1000) 
#variables of importance
	importance(model_rf)
	varImpPlot(model_rf, main = "Soil variables ranked by importance based on \n Random forest model across locations ")

# fitting RF on data per location
	data2 <- data1[data1$Location == "Haveri",]

	model_rf_2 = randomForest(x = data2[,-c(1:9, 12, 15, 22)], y = as.factor(data2$Block), ntree = 500) 

  #variables of importance
		varImpPlot(model_rf_2, main = "Soil variables ranked by importance \n based on Random forest model in Haveri")


# Fitting Random Forest to the reduced variable dataset
set.seed(14)  # Setting seed
model_rf_2 = randomForest(x = data1[,c(10, 13, 16, 20, 21)], y = as.factor(data1$Location), ntree = 1000) 
#variables of importance
	importance(model_rf_2)
	varImpPlot(model_rf_2, main = "Soil variables ranked by importance based on \n reduced variable Random forest model across locations ")


# prediction classificaton - whole data - prediction of location

cfm <- list()
Accuracy <- vector()
for (i in 1:100){
	ind <- sample(2, nrow(data1a),
              replace = TRUE,
              prob = c(0.7, 0.3))
	training <- data1a[ind==1,]
	testing <- data1a[ind==2,]

	# Fitting RF model for prediction
  rf.tr <- randomForest(x = training[,-c(1:8, 12)], 
  	        y = as.factor(training$Location), ntree = 1000) 

	# Predicting the Test set results
	trg = predict(rf.tr, newdata = testing[,-c(1:8, 12)])
  
  # Confusion Matrix
  cfm[[i]] = table(Predict = trg, Actual = testing$Location)
  Accuracy[i] <- tr(cfm[[i]])/sum(cfm[[i]])
}

cfm2 <- round(Reduce("+", cfm)/length(cfm),0) 
mat <- melt(cfm2)
colnames(mat) <- c("pred", "Actual", "n")

plot_confusion_matrix(mat, target_col = "Actual", prediction_col = "pred", 
	            counts_col = "n") +
ggtitle("Average confusion matrix obtained from 100 replications 
	of 70:30 training: test cross-validation with location   
	as the response variable in Random Forest model") +
  ylab("Prediction") +
  xlab("Actual") +
  theme_bw()

mean(Accuracy)
  hist(Accuracy, main = "Accuracy in classification from 100 replications
  	         of 70:30 training: test cross-validation with location 
	           as response variable in Random Forest model" )
  text(0.98, 60, "mean accuracy = 1")


# Fitting Random Forest - unsupervised

set.seed(5)  # Setting seed
rf_un = randomForest(x = data1a[,-c(1:8, 12)], ntree = 1000, proximity = TRUE) 

# proximity matrix
prox <- rf_un$proximity

distance <- sqrt(1-prox)

# clustering using distance matrix from RF
fviz_nbclust(distance, pam, method = "wss") # determine the number of clusters

set.seed(6)
# do k-medoid clustering using the proximity matrix- kk means - medoids 
pam.rf <- pam(distance, 5)


# add cluster to dataset
data_cluster <- data.frame(data1a, cluster = pam.rf$cluster)

data_cluster$cluster <- as.factor(data_cluster$cluster)

ggplot(data_cluster, aes(x = Block, y = Plot, fill = cluster)) +
  facet_grid(. ~ Location)  + 
  geom_tile() +
  scale_fill_discrete() +
  ggtitle("Clustering of plot segments in various locations \n based on similarity matrix from unsupervised RF model") +
  ylab("Plot") +
  xlab("Blocks") +
  theme_bw()

# supervised RF based on cluster
set.seed(7)  # Setting seed
model_rf_cl <-  randomForest(x = data_cluster[,-c(1:8, 12, 23)], y = as.factor(data_cluster$cluster), 
	ntree = 1000) 
 # maxnodes - play with it to see if accuracy changes.
 # mtry - p number of variables at each split
#variables of importance
	importance(model_rf_cl)
	varImpPlot(model_rf_cl, main = "Soil variables ranked by importance based on \n unsupervised Random forest model across clusters ")




# supervised RF based on cluster
set.seed(27)  # Setting seed
model_rf_cl_2 <-  randomForest(x = data_cluster[,-c(1:9, 12, 15, 22, 23)], y = as.factor(data_cluster$cluster), ntree = 1000) 
#variables of importance
	importance(model_rf_cl_2)
	varImpPlot(model_rf_cl_2, main = "Soil variables ranked by importance based on \n unsupervised Random forest model across locations ")


# prediction based on cluster

cfm <- list()
cfm_miss <- list()
missed <- rep(0,5)
seq.col <- seq(1:5)
Accuracy <- vector()

for (i in 1:100){
	ind <- sample(2, nrow(data_cluster),
          replace = TRUE,
          prob = c(0.7, 0.3))
	training <- data_cluster[ind==1,]
	testing <- data_cluster[ind==2,]

	# Fitting RF model for prediction
  rf.tr_cl <- randomForest(x = training[,-c(1:8, 12, 23)], 
  	        y = as.factor(training$cluster), ntree = 1000) 

	# Predicting the Test set results
	trg = predict(rf.tr_cl, newdata = testing[,-c(1:8, 12, 23)])
  
  # Confusion Matrix
  cfm[[i]] = table(Predict = trg, Actual = testing$cluster)
  Accuracy[i] <- tr(cfm[[i]])/sum(cfm[[i]])
	
} # end of for - prediction of clusters 


cfm2 <- round(Reduce("+", cfm)/length(cfm),0) 
mat <- melt(cfm2)
colnames(mat) <- c("pred", "Actual", "n")

plot_confusion_matrix(mat, target_col = "Actual", prediction_col = "pred", 
	            counts_col = "n") +
ggtitle("Average confusion matrix obtained from 100 replications 
	of 70:30 training: test cross-validation with unsupervised RF based 
	cluster as the response variable in supervised RF model") +
  ylab("Prediction") +
  xlab("Actual") +
  theme_bw()

  mean(Accuracy)
  hist(Accuracy, main = "Accuracy in classification from 100 replications
  	         of 70:30 training: test cross-validation with unsupervised RF based 
	           cluster as the response variable in supervised RF model" )
  text(0.94, 20, "mean accuracy = 0.97")

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



## Self-organizing maps (SOM) - provides plot segment comparison inter and intra comparision


# unsupervised SOM 

data1$Location <- as.factor(data1$Location)
data1b <- na.omit(data1)

# data set that is scaled and converted to be a matrix cause kohonen function accept numeric matrix
data1a <- as.matrix(scale(na.omit(data1[,-c(1:8, 12)])))

# SOM grid
data1.grid <- somgrid(xdim = 10, ydim = 10, topo = "hexagonal", toroidal = T, neighbourhood.fct =  "gaussian")

# SOM model
set.seed(8)
model_som <- som(data1a, data1.grid, rlen = 1000, radius = 2.5, keep.data = TRUE, dist.fcts = "euclidean")
# model_som <- som(data1a, data1.grid, rlen = 800, keep.data = TRUE, dist.fcts = "euclidean", mode = "online") # since gaussian function is specified, no radius required
plot(model_som, type = "changes") # increase the interations (rlen) if line is not plateauing



# classification of each nodes by codes plot
plot(model_som, type = "mapping", pchs = 19, shape = "round")

plot(model_som, type = "codes", main = "SOM codes plot for soil features", palette.name = rainbow)

# plot the average distance of each SOM unit to its neighbors
plot(model_som, type = "dist.neighbours") 

# find the distance between codes or objects
# This function is used in the calculation of the U matrix in function
		# plot.kohonen using the type = "dist.neighbours" argument.
# Value: An object of class dist, which can be directly fed into (e.g.) a hierarchical clustering.
object.distances(kohobj, type = c("data", "codes"), whatmap)

### see if cluster labels can be derived drom plot('dist.neighbours', or from object.distances) 
 ## if not, use this distance matrix in kmeans clustring.


# Heatmap - visualisation of the distribution of a single variable across the map
p <- ncol(data1a)
heatmap.som <- function(model){
  for (i in 1:p) {
   plot(model, type = "property", property = getCodes(model)[,i], 
        main = colnames(getCodes(model))[i]) 
  }
}
heatmap.som(model_som)

plot(model_som, type = "property", property = getCodes(model_som)[,1], 
        main = colnames(getCodes(model_som))[1]) 

# cluster the observations based on SOM
fviz_nbclust(model_som$codes[[1]], pam, method = "wss") # determine the number of clusters

# pam clustering
set.seed(9)
pam.som <- pam(model_som$codes[[1]], 5) 

plot(model_som, type = "codes", bgcol = rainbow(5)[pam.som$cluster], main = "Cluster map based on SOM of soil features")
add.cluster.boundaries(model_som, pam.som$cluster)

vec <- rep(1,10)
names(vec) <- c(1:10)
barplot(vec, col = rainbow(10), yaxt='n', ylab="" , xlab = "Cluster number", main = "Color legend of clusters")

# add cluster to dataset
data1_cluster <- data.frame(data1b, cluster = pam.som$cluster[model_som$unit.classif])

data1_cluster$cluster <- as.factor(data1_cluster$cluster)

ggplot(data1_cluster, aes(x = Block, y = Plot, fill = cluster)) +
  facet_grid(. ~ Location)  + 
  geom_tile() +
  scale_fill_discrete() +
  ggtitle("Clustering of plot segments in various locations \n based on unsupervised SOM") +
  ylab("Plot") +
  xlab("Blocks") +
  theme_bw()



# predictability of SOM- test if the new labels 'cluster' can be predicted correctly. 
  # the labels are derived from clustering based on soil features. So, if labels can be predicted
   # correctly, it is a test of the predictability of SOM in classifying the soil based on their features. 
cfm <- list()
cfm_miss <- list()
missed <- rep(0,10)
seq.col <- seq(1:10)
Accuracy <- vector()
for (i in 1:100){
	ind <- sample(2, nrow(data1_cluster),
          replace = TRUE,
          prob = c(0.7, 0.3))
	training <- data1_cluster[ind==1,]
	testing <- data1_cluster[ind==2,]

	# set.seed(123) # guarantee reproduiable results
	# training the model
	train1 <- as.matrix(scale(training[,-c(1:8, 12, 23)]))
	# test1<- scale(testing[,-c(1:9, 12, 15, 22, 23)], center = attr(train1, "scaled:center"), scale = attr(train1, "scaled:scale")) # scaling of test data using training attributes
	test1 <- scale(testing[,-c(1:8, 12, 23)], center = attr(train1, "scaled:center")) # scaling of test data using training attributes


	testXY <- list(independent = test1, dependent = factor(testing$cluster))

	train_som <- xyf(train1, classvec2classmat(factor(training$cluster)), data1.grid, rlen = 1000, radius = 2.5, 
											keep.data = TRUE, dist.fcts = "euclidean")
	plot(train_som, type = "changes") 

	# testing the model
	# test_cluster <- predict(train_som, newdata = data.matrix(testXY), whatmap =  c("independent", "dependent"))
	test_cluster <- predict(train_som, newdata = as.matrix(test1), whatmap = 1)
	# cfm[[i]] <- table(Predict = test_cluster$predictions[[2]], Actual = testing$cluster, useNA = 'always')
	cfm[[i]] <- table(Predict = test_cluster$predictions[[2]], Actual = testing$cluster)
	Accuracy[i] <- tr(cfm[[i]])/sum(cfm[[i]])
} # end of for


cfm2 <- round(Reduce("+", cfm)/length(cfm),0) 
mat <- melt(cfm2)
colnames(mat) <- c("pred", "Actual", "n")

plot_confusion_matrix(mat, target_col = "Actual", prediction_col = "pred", 
	            counts_col = "n") +
ggtitle("Average confusion matrix obtained from 100 replications  
	       of 70:30 training: test cross-validationwith unsupervised SOM
         based 10 clusters as the response variable in supervised SOM") +
  ylab("Prediction") +
  xlab("Actual") +
  theme_bw()

  mean(Accuracy)
  hist(Accuracy, main = "Accuracy in classification from 100 replications
  	         of 70:30 training: test cross-validation with unsupervised SOM based 
	           10 clusters as the response variable in supervised SOM model" )
  text(0.92, 20, "mean accuracy = 0.98")




# supervised SOM - using Location as the response variable

data1a <- as.matrix(scale(na.omit(data1[,-c(1:8, 12)])))  
data1b <- na.omit(data1)

set.seed(10)
model2_som <- xyf(data1a, classvec2classmat(factor(data1b$Location)), data1.grid, rlen = 1000, radius = 2.5, 
											keep.data = TRUE, dist.fcts = "euclidean")
	plot(model2_som, type = "changes") 


# classification of each nodes by codes plot
plot(model2_som, type = "mapping", pchs = 19, shape = "round")

plot(model2_som, type = "codes", main = "SOM codes plot for soil features", palette.name = rainbow)

plot(model2_som, type = "dist.neighbours") # darker colored nodes - closer vector inputs


# predictability of supervised SOM with Location as the response variable
cfm <- list()
cfm_miss <- list()
missed <- rep(0,13)
seq.col <- seq(1:13)
Accuracy <- vector()

for (i in 1:100){
	ind <- sample(2, nrow(data1b),
          replace = TRUE,
          prob = c(0.7, 0.3))
	training <- data1b[ind==1,]
	testing <- data1b[ind==2,]

	# set.seed(123) # guarantee reproduiable results
	# training the model
	train1 <- as.matrix(scale(training[,-c(1:8, 12)]))
	# test1<- scale(testing[,-c(1:9, 12, 15, 22, 23)], center = attr(train1, "scaled:center"), scale = attr(train1, "scaled:scale")) # scaling of test data using training attributes
	test1 <- scale(testing[,-c(1:8, 12)], center = attr(train1, "scaled:center")) # scaling of test data using training attributes


	testXY <- list(independent = test1, dependent = factor(testing$Location))

	train_som <- xyf(train1, classvec2classmat(factor(training$Location)), data1.grid, rlen = 1000, radius = 2.5, 
											keep.data = TRUE, dist.fcts = "euclidean")
	plot(train_som, type = "changes") 

	# testing the model
	# test_cluster <- predict(train_som, newdata = data.matrix(testXY), whatmap =  c("independent", "dependent"))
	test_pred <- predict(train_som, newdata = as.matrix(test1), whatmap = 1)
	# cfm[[i]] <- table(Predict = test_cluster$predictions[[2]], Actual = testing$cluster, useNA = 'always')
	cfm[[i]] <- table(Predict = test_pred$predictions[[2]], Actual = testing$Location)
	Accuracy[i] <- tr(cfm[[i]])/sum(cfm[[i]])
} # end of for


cfm2 <- round(Reduce("+", cfm)/length(cfm),0) 
mat <- melt(cfm2)
colnames(mat) <- c("pred", "Actual", "n")

plot_confusion_matrix(mat, target_col = "Actual", prediction_col = "pred", 
	            counts_col = "n") +
ggtitle("Average confusion matrix obtained from 100 replications  
	       of 70:30 training: test cross-validation with supervised SOM
         using Locations as the response variable") +
  ylab("Prediction") +
  xlab("Actual") +
  theme_bw()

  mean(Accuracy)
  hist(Accuracy, main = "Accuracy in classification from 100 replications
  	         of 70:30 training: test cross-validation with supervised SOM 
	           using Locations as the response variable" )
  text(0.99, 60, "mean accuracy = 1")


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

