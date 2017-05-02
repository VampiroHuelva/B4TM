# In this script we do a comples analysis of the data with a simple feature selection and adding a vector with a vector selection used in weka. Algorithm used will be on the report.
# We use this script to first select the number of features is optimus for each algorithm, and then select the best algortihm
rm(list = ls())
#Preprocesing data

clinical <- read.delim("Train_clinical.txt", header =TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
call <- read.delim("Train_call.txt", header =TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
head(clinical[,1:2])
head(call[,1:10])
call = t(as.data.frame(call))
call = call[5:length(call[,1]),]
head(call[,1:10])
rownames(clinical) <- clinical[,1] 
head(clinical[,1:2])
combine <- merge(clinical, call, by="row.names")
row.names(combine)<-combine$Row.name
combine$Row.names<-NULL
combine$Sample<-NULL
head(combine[,1:10])



## write out to file WEKA
write.arff(combine, file = "combine.arff")

## Feature selection by importance ()
fil=filterVarImp(combine[,2:ncol(combine)],combine[,1],nonpara = FALSE)

############### FEATURE SELECTION NOT DONE WELL, PROCEDING ANYWAY

# simple feature selection
fill = sort(rowMeans(fil), decreasing = T)
features <- names(fill[1:30])

# Advance feature selction using WEKA. Algorithm used describe in the report

CorrelationAttribute <- c("V2185","V2214","V2215","V2211","V1657","V1679","V675","V1668","V674","V1680","V2207","V2224","V2224","V1673","V673","V1669","V1665","V2220","V693","V696","V855","V1678","V2212","V1664","V1670","V1671","V2225","V1674","V853","V1675","V850","V1663","V672","V680","V312","V673","V696","V773","V1570","V1657","V2185","V2593","V2663","V2733","V2185")

J48_previous <- c("V2733", "V1647","V1647","V1564","V624","V624","V1679","V1679","V2548","V148","V2548","V227","V2743")
filtering <- c("V312","V673","V696","V773","V1570","V1657","V2185","V25932","V2663","V2733")
Huge_combination <- c(J48_previous,CorrelationAttribute[1:25],filtering, features[1:20])
Huge_Unic <- unique(Huge_combination)

#################### SETTING UP THE SCHEME ###############################

### DATASETS: GENERATE CROSSVALIDATION (TRAIN-TEST) SETS AND FINAL TEST SETS ###
crossval = 5
loops = 5

# group data by class
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


# lets train nnet, ft, svmRadial,J48 and gbm
# making a list with n different models

####  Getting scores of the models to choose the number of features to use and the best model
Nmodel = 3
models <- vector("list",Nmodel)
for (i in 1:Nmodel) {
  models[[i]] <- vector("list",loops)
  for (j in 1:loops)
    models[[i]][[j]] <- vector("list",crossval)
}

for (i in 1:loops) {
  for (j in 1:crossval) {
    
    
    
    features =Huge_Unic[1:10+i]
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
    features =Huge_Unic[1:15+i]
    #[1:(5+j)]
    model = paste(features, collapse=' + ')
    model = paste("Subgroup ~ ", model, collapse ='')
    models[[2]][[i]][[j]] <- train(eval(parse(text=model)), data = sets[[i]][[1]][[j]][[1]], method = "nnet", verbose = FALSE)
    
    #svm
    features =Huge_Unic[1:20+i]
    fitControl <- trainControl(## 10-fold CV
      method = "repeatedcv",
      number = 10,
      ## repeated ten times
      repeats = 10)
    #features[1:(5+j)]
    model = paste(features, collapse=' + ')
    model = paste("Subgroup ~ ", model, collapse ='')
    models[[3]][[i]][[j]] <- train(eval(parse(text=model)), data = sets[[i]][[1]][[j]][[1]], method = "nnet", trControl = fitControl, verbose = FALSE)
  }    
}
### COMPARE THE MODELS NOW ###

for (k in 1:Nmodel) {
  for (i in 1:loops) {
    r = c()#r[i]=
    for (j in 1:crossval) {
      #varimp <- varImp(models[[k]][[i]][[j]], varImp.train=FALSE)
      pred = predict(models[[k]][[i]][[j]],data=sets[[i]][[1]][[j]][[2]],na.omit(sets[[i]][[1]][[j]][[2]]))
      #print(confusionMatrix(sets[[i]][[1]][[j]][[2]][,1], pred)[2:3])
      if (1) {
        x = as.numeric(sets[[i]][[1]][[j]][[2]][,1])
        y = as.numeric(pred)
        roc = multiclass.roc(x,y)
        r <- append(r,auc(roc))
        rs <- roc[['rocs']]
        #plot.roc(rs[[1]])
        #sapply(2:length(rs),function(h) lines.roc(rs[[h]],col=h))
        
        
      }
#scoring of the dirent models selected with the number of features selected      
    }
    print(mean(r))
    print(k)
  }
}