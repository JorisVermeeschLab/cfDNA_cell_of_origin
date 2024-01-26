library(data.table)
library(tidyverse)
library(e1071)
library(pROC)

args <- commandArgs(TRUE)

output_dir <- args[1]
myseed <- args[2]
cost <- args[3]


#pivot wider combined data file
combined <- fread(paste0(output_dir, "/combined_data.txt"), header=T)
combined$group <- paste0(combined$sample, ":", combined$status)
combined.wider <- combined %>% select(-sample, -status) %>% pivot_wider(names_from = cell_type, values_from = rank) %>% as.data.frame()
rownames(combined.wider) <- combined.wider$group
dtmat <- combined.wider %>% select(-group) %>% data.matrix()

set.seed(myseed)

#set variables
orgroup <- row.names(dtmat)
dtmat_group <- sapply(strsplit(row.names(dtmat), split=':', fixed=TRUE), function(x) (x[2]))
dtmat_group <- as.factor(dtmat_group)

#get weights
wts <- 100 / table(dtmat_group)
n_train <- nrow(dtmat)
                      
alltruegroup <- c()
allpredgroup <- c()
allpredval <- c()
predres <- c()
                      
#run leave-one-out svm
for (x in 1:n_train) {
  traindt <- dtmat[-x,]
  traingroup <- dtmat_group[-x]
  testdt <- t(dtmat[x,])
  #train.pca <- prcomp(traindt)
  #print(dim(train.pca$x))
  #project.test <- predict(train.pca, testdt)
  testgroup <- dtmat_group[x]
  svm.model <- svm(traindt, traingroup, class.weights = wts, kernel="linear", cost=cost, scale=FALSE, decision.values=TRUE)
  svm.pred <- predict(svm.model,testdt,decision.values=TRUE)
  alltruegroup <- c(alltruegroup,testgroup)
  allpredgroup <- c(allpredgroup,attributes(svm.pred)$decision.values)
  allpredval <- c(allpredval,svm.pred)
  predres <- rbind(predres, data.frame(orgroup[x],attributes(svm.pred)$decision.values))
  print(data.frame(orgroup[x],svm.pred))
}
write.table(predres, paste0(output_dir, "svm_loo_decision_values.txt", sep="\t", row.names=F, col.names=T, quote=F)
svmperfm <- roc(response=alltruegroup,predictor=allpredgroup,ci=TRUE)
print(svmperfm$auc)
print(svmperfm$ci)
runacc <- table(alltruegroup,allpredval)
print(runacc)
acclist <- (runacc[1,1]+runacc[2,2])/sum(runacc)
print(paste0("Accuracy: ", acclist))

pdf(paste0(output_dir, "/svm_loo_roc.pdf"))
plot(svmperfm,print.auc=TRUE)
dev.off()
