############################################################################################
###
### functions used in main.R
###
############################################################################################

load_depend = function() {
  
  source("code/model_tuning.R")
  source("code/feature_selection.R")
  
  # Loading all the dependicies
  library(splines); library(parallel); library(survival); library(caret); library(mlbench)
  library(gbm); library(corrplot); library(pROC); library(FSelector); library(qtlcharts)
  library(lattice); library(energy); library(RWeka); library(obliqueRF); library(stepPlr);
  library(dplyr); library(LiblineaR); library(kernlab)
  
  #write.arff(combine, file = "combine.arff")
  
}

# load unlabeled data
Load_Unlabeled_Data <- function(input){
  
  Validation_call <- read.delim(input, header =TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
  head(Validation_call[,1:10])
  Validation_call = t(as.data.frame(Validation_call))
  Validation_call = Validation_call[5:length(Validation_call[,1]),]
  
}

# load labeled data
Load_labeled_Data = function() {
  
  clinical <- read.delim("Data/Train_clinical.txt", header =TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
  call <- read.delim("Data/Train_call.txt", header =TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
  call = t(as.data.frame(call))
  call = call[5:length(call[,1]),]
  rownames(clinical) <- clinical[,1] 
  combine <- merge(clinical, call, by="row.names")
  row.names(combine)<-combine$Row.name
  combine$Row.names<-NULL
  combine$Sample<-NULL
  levels(combine$Subgroup)[1]<-'HER2Plus'
  levels(combine$Subgroup)[2]<-'HRPlus'
  levels(combine$Subgroup)[3]<-'TripleNeg'
  
  return(combine)
  
}

# make train and test sets for double loop cross validation (workswith 3 classes)
double_crossval_datasets = function(combine, crossval, loops) {
  
  # group data by class
  HER2plus <-combine[combine$Subgroup=="HER2Plus",]
  HRplus <- combine[combine$Subgroup=="HRPlus",]
  TrNeg <- combine[combine$Subgroup=="TripleNeg",]
  
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
  
}

# make train and test sets for single loop cross validation (workswith 3 classes)
crossval_sets = function(x, crossval) {
  
  # group data by class
  HER2plus <-x[x$Subgroup=="HER2Plus",]
  HRplus <- x[x$Subgroup=="HRPlus",]
  TrNeg <- x[x$Subgroup=="TripleNeg",]
  
  cross <- function(y,n) split(y, factor(sort(rank(y)%%n)))
  
  # randomize in the groups
  HE <- sample(1:length(HER2plus[,1]),replace = F)
  HR <- sample(1:length(HRplus[,1]),replace = F)
  Tr <- sample(1:length(TrNeg[,1]),replace= F)
  
  # divide in n testsets
  total_HE = cross(HE,crossval)
  total_HR = cross(HR,crossval)
  total_Tr = cross(Tr,crossval)

  sets <- vector("list", crossval)
  
  if (crossval != 1) {
    
    for (i in 1:crossval) {
      # extracting the testset
      test <- rbind(HER2plus[total_HE[[i]],],HRplus[total_HR[[i]],],TrNeg[total_Tr[[i]],])
      
      # extracting the trainset 
      HE_t = unlist(setdiff(total_HE,total_HE[i]))
      HR_t = unlist(setdiff(total_HR,total_HR[i]))
      Tr_t = unlist(setdiff(total_Tr,total_Tr[i]))
      
      train = rbind(HER2plus[HE_t,],HRplus[HR_t,],TrNeg[Tr_t,])
      
      #saving the sets
      sets[[i]]<-list(train,test)
    }
    
   } else {
     
    sample_size = floor(0.75 * nrow(HER2plus))
    train_ind <- sample(seq_len(nrow(HER2plus)), size = sample_size)
    train_HER2 <- HER2plus[train_ind, ]
    test_HER2 <- HER2plus[-train_ind, ]
    
    sample_size = floor(0.75 * nrow(HRplus))
    train_ind <- sample(seq_len(nrow(HRplus)), size = sample_size)
    train_HR <- HRplus[train_ind, ]
    test_HR <- HRplus[-train_ind, ]
    
    sample_size = floor(0.75 * nrow(TrNeg))
    train_ind <- sample(seq_len(nrow(TrNeg)), size = sample_size)
    train_Tr <- TrNeg[train_ind, ]
    test_Tr <- TrNeg[-train_ind, ]

    train <- rbind(train_HER2,train_HR,train_Tr)
    test <- rbind(test_HER2,test_HR,test_Tr)
    sets[[1]]<-list(train,test)
   }
  
  return(sets)
  
}

# make train and test sets for single loop cross validation (workswith 2 classes)
crossval_sets2 = function(x, crossval) {
  
  # group data by class
  c1 <-x[x$Subgroup==levels(x$Subgroup)[1],]
  c2 <- x[x$Subgroup==levels(x$Subgroup)[2],]

  cross <- function(y,n) split(y, factor(sort(rank(y)%%n)))
  
  # randomize in the groups
  c1ran <- sample(1:length(c1[,1]),replace = F)
  c2ran <- sample(1:length(c2[,1]),replace = F)

  # divide in n testsets
  total_c1ran = cross(c1ran,crossval)
  total_c2ran = cross(c2ran,crossval)

  sets <- vector("list", crossval)

  if (crossval != 1) {
    
    for (i in 1:crossval) {
      # extracting the testset
      test <- rbind(c1[total_c1ran[[i]],],c2[total_c2ran[[i]],])
      
      # extracting the trainset 
      c1ran_t = unlist(setdiff(total_c1ran,total_c1ran[[i]]))
      c2ran_t = unlist(setdiff(total_c2ran,total_c2ran[[i]]))

      train = rbind(c1[c1ran_t,],c2[c2ran_t,])
      
      #saving the sets
      sets[[i]]<-list(train,test)
    }
    
  } else {
    
    sample_size = floor(0.75 * nrow(c1))
    train_ind <- sample(seq_len(nrow(c1)), size = sample_size)
    train_c1 <- c1[train_ind, ]
    test_c1 <- c1[-train_ind, ]
    
    sample_size = floor(0.75 * nrow(c2))
    train_ind <- sample(seq_len(nrow(c2)), size = sample_size)
    train_c2 <- c2[train_ind, ]
    test_c2 <- c2[-train_ind, ]

    train <- rbind(train_c1,train_c2)
    test <- rbind(test_c1,test_c2)
    sets[[1]]<-list(train,test)
  }
  
  return(sets)
  
}

# Make 3 sets all with two classes
twoclasses = function(x) {
  
  HER2 = x
  levels(HER2$Subgroup)[2]<-"NOHER2Plus"
  levels(HER2$Subgroup)[3]<-"NOHER2Plus"
  
  HR = x
  levels(HR$Subgroup)[1]<-'NOHRPlus'
  levels(HR$Subgroup)[3]<-'NOHRPlus'
  
  Triple = x
  levels(Triple$Subgroup)[1]<-'NOTripleNeg'
  levels(Triple$Subgroup)[2]<-'NOTripleNeg'
  sets <- list(HER2,HR,Triple)
  
  return(sets)
  
}

#make the two classes datasets with classes of same proportions
make_equal_sets = function(x,class) {
  
  classes = levels(x$Subgroup)
  # check which is the longest
  # add duplicates to the lowerst group (double)
  class1 <-x[x$Subgroup==classes[1],]
  class2 <- x[x$Subgroup==classes[2],]
  if (length(class1[,1]) > length(class2[,2])) {
    result = rbind(class1,class2,class2)
  } else {
    result = rbind(class1,class1,class2)
  }
  
  return(result)
  
}

train_all_models = function(x, features){
  
  gbmFit = gbm_tuning(x,features) 
  svmFit = svm_tuning(x,features)
  nnetFit = nnet_tuning(x,features)
  mrFit = mr_tuning(x,features)
  rfFit = rf_tuning(x,features)
  resamps <- resamples(list(GBM = gbmFit, SVM = svmFit, NNET = nnetFit, MR = mrFit, RF = rfFit))
  
  return(list(gbmFit,svmFit,nnetFit,mrFit,rfFit,resamps))
  
}

train_all_models_with_feature_selection = function(train_set,features) {
  
  gbmFit = mr_tuning(train_set,features)
  gbmFit_features = feature_var_imp(gbmFit, 10)
  gbmFit = gbm_tuning(train_set,features)
  
  svmFit = svm_tuning(train_set,features)
  
  nnetFit = nnet_tuning(train_set,features)
  
  mrFit = mr_tuning(train_set,features)
  mrFit_features = feature_var_imp(mrFit, 10)
  mrFit = mr_tuning(train_set,mrFit_features) 
  
  rfFit = rf_tuning(train_set,features)
  rfFit_features = feature_var_imp(rfFit, 10)
  rfFit = rf_tuning(train_set,rfFit_features)
  
  resamps <- resamples(list(GBM = gbmFit, SVM = svmFit, NNET = nnetFit, MR = mrFit, RF = rfFit))
  
  return(list(gbmFit,svmFit,nnetFit,mrFit,rfFit,resamps))
  
}

results_prediction = function(pred) {
  
  gbm = c(sum(pred$gbmFit == pred$real), sum(pred$gbmFit[pred$real == 1] == pred$real[pred$real == 1]),
          sum(pred$gbmFit[pred$real == 2] == pred$real[pred$real == 2]),
          sum(pred$gbmFit[pred$real == 3] == pred$real[pred$real == 3]))
  
  svm = c(sum(pred$svmFit == pred$real), sum(pred$svmFit[pred$real == 1] == pred$real[pred$real == 1]),
          sum(pred$svmFit[pred$real == 2] == pred$real[pred$real == 2]),
          sum(pred$svmFit[pred$real == 3] == pred$real[pred$real == 3]))
  
  nnet = c(sum(pred$nnetFit == pred$real), sum(pred$nnetFit[pred$real == 1] == pred$real[pred$real == 1]),
         sum(pred$nnetFit[pred$real == 2] == pred$real[pred$real == 2]),
         sum(pred$nnetFit[pred$real == 3] == pred$real[pred$real == 3]))
  
  mr = c(sum(pred$mrFit == pred$real), sum(pred$mrFit[pred$real == 1] == pred$real[pred$real == 1]),
         sum(pred$mrFit[pred$real == 2] == pred$real[pred$real == 2]),
         sum(pred$mrFit[pred$real == 3] == pred$real[pred$real == 3]))
  
  rf = c(sum(pred$rfFit == pred$real), sum(pred$rfFit[pred$real == 1] == pred$real[pred$real == 1]),
         sum(pred$rfFit[pred$real == 2] == pred$real[pred$real == 2]),
         sum(pred$rfFit[pred$real == 3] == pred$real[pred$real == 3]))
  total = c(length(pred$real), sum(pred$real == 1), sum(pred$real == 2), sum(pred$real == 3))
  results = data.frame(gbm, svm, nnet, mr, rf, total)
  
  return(results)
  
}

save_models = function(models, models_trained) {

  models[[length(models)+1]]  <- list(gbmFit = models_trained[[1]],
                               svmFit = models_trained[[2]],
                               nnetFit = models_trained[[3]],
                               mrFit = models_trained[[4]],
                               rfFit = models_trained[[5]])
  
  return(models)
  
}

save_pred =function(pred, models_trained, test_set) {

  temp = cbind(predict(models_trained[[1]],test_set),
               predict(models_trained[[2]],test_set),
               predict(models_trained[[3]],test_set),
               predict(models_trained[[4]],test_set),
               predict(models_trained[[5]],test_set),
               test_set$Subgroup)
  for (i in 1:length(temp[,1])) {
    pred[nrow(pred)+1,] <- temp[i,]
  }
  return(pred)
  
}

save_accur = function(accur, resample) {

  accur[nrow(accur)+1,] <- c(summary(resample)$statistics$Accuracy[1,4],
                             summary(resample)$statistics$Accuracy[2,4],
                             summary(resample)$statistics$Accuracy[3,4],
                             summary(resample)$statistics$Accuracy[4,4],
                             summary(resample)$statistics$Accuracy[5,4])
  
  return(accur)
  
}

empty_accur = function() {
  
  accur <- data.frame(x= numeric(0), y= numeric(0), z = numeric(0),
                      x= numeric(0), y= numeric(0))
  colnames(accur) <- c("gbmFit","svmFit","nnetFit","mrFit","rfFit")
  
  return(accur)
  
}


empty_pred = function() {
  
  pred <- data.frame(x= numeric(0), y= numeric(0), z = numeric(0),
                     x= numeric(0), y= numeric(0), z = numeric(0))
  colnames(pred) <- c("gbmFit","svmFit","nnetFit","mrFit","rfFit","real")
  
  return(pred)
  
}

empty_models = function() {
  
  models <- vector("list", 5)
  
  return(models)
  
}

plot_results = function(accur, results) {
  
  par(mfrow=c(2,2))
  main=1
  boxplot(accur,xlab='models', ylab='accuracy')
  boxplot(t(accur),xlab='cross validations', ylab='accuracy')
  boxplot(t(results[,1:5]),xlab='average predictions', ylab='correct')
  coll = rainbow(5)
  plot(0,0, xlim = c(0,10), ylim = c(0.7,1))
  for ( i in seq(1,length( accur ),1) ) lines(accur[,i],ylab=names(accur[i]),type="l", col = coll[i])
  par(mfrow=c(1,1))
  legend("bottomright",legend=colnames(accur),col=rainbow(5),pch=1)
}
