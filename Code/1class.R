rm(list=ls())

library(splines)
library(parallel)
library(survival)
library(caret)
library(mlbench)
library(gbm)
library(corrplot)
library(pROC)
library(FSelector)
library(qtlcharts)
library(lattice)
library(energy)
library(RWeka)

## preprocessing the data (TODO REDUCE REDUNDANCY IN THE DATA!!!!!!!!!!!!!!!)
clinical <- read.delim("Train_clinical(1).txt", header =TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
call <- read.delim("Train_call(1).txt", header =TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
call = t(as.data.frame(call))
call = call[5:length(call[,1]),]
rownames(clinical) <- clinical[,1] 
combine <- merge(clinical, call, by="row.names")
row.names(combine)<-combine$Row.name
combine$Row.names<-NULL
combine$Sample<-NULL

## write out to file WEKA
write.arff(combine, file = "combine.arff")

## Feature selection by importance ()
fil=filterVarImp(combine[,2:ncol(combine)],combine[,1],nonpara = FALSE)

############### FEATURE SELECTION NOT DONE WELL, PROCEDING ANYWAY

# simple feature selection
fill = sort(rowMeans(fil), decreasing = T)
features <- names(fill[1:30])

#################### SETTING UP THE SCHEME ###############################

### DATASETS: GENERATE CROSSVALIDATION (TRAIN-TEST) SETS AND FINAL TEST SETS ###
crossval = 5
loops = 5

classes = c('HER2+',"HR+",'Triple Neg')


# Make 3 sets all with two classes (one class against not this one class)
twoclasses = function(x) {
  
  HER2_VS_rest = x
  HR = x
  Triple = x
  
  HER2$Subgroup[HER2_VS_rest$Subgroup == classes[2]] <- "NOHER2+"
  HER2$Subgroup[HER2_VS_rest$Subgroup == classes[3]] <- "NOHER2+"
  
  HR$Subgroup[HR$Subgroup == classes[1]] <- "NOHR+"
  HR$Subgroup[HR$Subgroup == classes[3]] <- "NOHR+"
  
  Triple$Subgroup[Triple$Subgroup == classes[1]] <- "NOTriple Neg"
  Triple$Subgroup[Triple$Subgroup == classes[2]] <- "NOTriple Neg"
  
  return(list(HER2,HR,Triple))
}

  
  


HER2plus <-combine[combine$Subgroup=="HER2+",]
HRplus <- combine[combine$Subgroup=="HR+",]
TrNeg <- combine[combine$Subgroup=="Triple Neg",]

cross <- function(x,n) split(x, factor(sort(rank(x)%%n)))

# randomize in the groups
HE <- sample(1:length(HER2plus[,1]),replace = F)
HR <- sample(1:length(HRplus[,1]),replace = F)
Tr <- sample(1:length(TrNeg[,1]),replace= F)

# divide in n testsets
total_HE = cross(HE,loops)
total_HR = cross(HR,loops)
total_Tr = cross(Tr,loops)

sets <- vector("list", loops)

for (i in 1:loops) {
  # extracting the testset
  test <- rbind(HER2plus[total_HE[[i]],],HRplus[total_HR[[i]],],TrNeg[total_Tr[[i]],])
  
  # extracting the trainset 
  # indices 
  HE_train = unlist(setdiff(total_HE,total_HE[i]))
  HR_train = unlist(setdiff(total_HR,total_HR[i]))
  Tr_train = unlist(setdiff(total_Tr,total_Tr[i]))
  
  # dividing in m crossvalidations
  HE_t <- cross(HE_train,crossval)
  HR_t <- cross(HR_train,crossval)
  Tr_t <- cross(Tr_train,crossval)
  
  # saving all sets combinations
  cross_sets <- vector("list",crossval)  
  for (j in 1:crossval) {
    test_cross <- rbind(HER2plus[HE_t[[j]],],HRplus[HR_t[[j]],],TrNeg[Tr_t[[j]],])
    HE_tr = unlist(setdiff(HE_t,HE_t[j]))
    HR_tr = unlist(setdiff(HE_t,HE_t[j]))
    Tr_tr = unlist(setdiff(HE_t,HE_t[j]))
    train_cross = rbind(HER2plus[HE_tr,],HRplus[HR_tr,],TrNeg[Tr_tr,]) 
    dim(train_cross)
    cross_sets[[j]]<- list(train_cross,test_cross)
  } 
  sets[[i]]<-list(cross_sets,test)
}

### TRAINING MODELS ###

# making a list with n different models
Nmodel = 3
models <- vector("list",Nmodel)
for (i in 1:Nmodel) {
  models[[i]] <- vector("list",loops)
  for (j in 1:loops)
    models[[i]][[j]] <- vector("list",crossval)
}

for (i in 1:loops) {
  for (j in 1:crossval) {
    ## build the three models
    
    # nnet
    features[1:(5+j)]
    model = paste(features, collapse=' + ')
    model = paste("Subgroup ~ ", model, collapse ='')
    model
    fitControl <- trainControl(## 10-fold CV
      method = "repeatedcv",
      number = 10,
      ## repeated ten times
      repeats = 10)
    models[[1]][[i]][[j]] <- train(eval(parse(text=model)), data = sets[[i]][[1]][[j]][[1]], method = "nnet", trControl = fitControl, verbose = FALSE)
    
    # gbm
    features[1:(5+j)]
    model = paste(features, collapse=' + ')
    model = paste("Subgroup ~ ", model, collapse ='')
    models[[2]][[i]][[j]] <- train(eval(parse(text=model)), data = sets[[i]][[1]][[j]][[1]], method = "gbm", verbose = FALSE)
    
    #svm
    fitControl <- trainControl(## 10-fold CV
      method = "repeatedcv",
      number = 10,
      ## repeated ten times
      repeats = 10)
    features[1:(5+j)]
    model = paste(features, collapse=' + ')
    model = paste("Subgroup ~ ", model, collapse ='')
    models[[3]][[i]][[j]] <- train(eval(parse(text=model)), data = sets[[i]][[1]][[j]][[1]], method = "svmRadial", trControl = fitControl, verbose = FALSE)
  }    
}

### COMPARE THE MODELS NOW ###

for (k in 1:Nmodel) {
  for (i in 1:loops) {
    for (j in 1:crossval) {
      #varimp <- varImp(models[[k]][[i]][[j]], varImp.train=FALSE)
      pred = predict(models[[k]][[i]][[j]],data=sets[[i]][[1]][[j]][[2]],na.omit(sets[[i]][[1]][[j]][[2]]))
      print(confusionMatrix(sets[[i]][[1]][[j]][[2]][,1], pred)[2:3])
      if (1) {
        x = as.numeric(sets[[i]][[1]][[j]][[2]][,1])
        y = as.numeric(pred)
        roc = multiclass.roc(x,y)
        print(auc(roc))
        rs <- roc[['rocs']]
        plot.roc(rs[[1]])
        print("k")
        print(k)
        print("i")
        print(i)
        print("j")
        print(j)
        sapply(2:length(rs),function(h) lines.roc(rs[[h]],col=h))
      }
    }
  }
}

models[[1]][[1]][[1]]

rs[[2]]

