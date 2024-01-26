library(data.table)
library(tidyverse)
library(e1071)
library(pROC)



orgroup <- row.names(dtmat)
dtmat_group <- sapply(strsplit(row.names(dtmat), split=':', fixed=TRUE), function(x) (x[2]))
dtmat_group <- as.factor(dtmat_group)
dtmat_num_group <- gsub("control", "1", dtmat_group)
dtmat_num_group <- gsub("case", "2", dtmat_group)
dtmat_num_group <- as.factor(dtmat_group)
wts <- 100 / table(dtmat_group)
n_train <- nrow(dtmat)
alltruegroup <- c()
allpredgroup <- c()
allpredval <- c()

#run leave-one-out svm
for (x in 1:n_train) {
  traindt <- dtmat[-x,]
  traingroup <- dtmat_group[-x]
  testdt <- t(dtmat[x,])
  #train.pca <- prcomp(traindt)
  #print(dim(train.pca$x))
  #project.test <- predict(train.pca, testdt)
  testgroup <- dtmat_group[x]
  svm.model <- svm(traindt, traingroup, class.weights = wts, kernel="linear", cost=1, scale=FALSE,decision.values=TRUE)
  svm.pred <- predict(svm.model,testdt,decision.values=TRUE)
  alltruegroup <- c(alltruegroup,testgroup)
  allpredgroup <- c(allpredgroup,attributes(svm.pred)$decision.values)
  allpredval <- c(allpredval,svm.pred)
  predres <- data.frame(orgroup[x],svm.pred)
  print(predres)
}
svmperfm <- roc(response=alltruegroup,predictor=allpredgroup,ci=TRUE)
print(svmperfm$auc)
print(svmperfm$ci)
runacc <- table(alltruegroup,allpredval)
print(runacc)
acclist <- (runacc[1,1]+runacc[2,2])/sum(runacc)
print(acclist)
