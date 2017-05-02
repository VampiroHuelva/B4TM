### DATASETS: GENERATE CROSSVALIDATION SETS ###
setwd("C:/Program Files/R/R-3.3.3/bin/x64")

## preprocessing the data (TODO REDUCE REDUNDANCY IN THE DATA!!!!!!!!!!!!!!!)
clinical <- read.delim("Train_clinical.txt", header =TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
call <- read.delim("Train_call.txt", header =TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
call = t(as.data.frame(call))
call = call[5:length(call[,1]),]
rownames(clinical) <- clinical[,1] 
combine <- merge(clinical, call, by="row.names")
row.names(combine)<-combine$Row.name
combine$Row.names<-NULL
combine$Sample<-NULL

Load_Unlabeled_Data <- function(input){

  Validation_call <- read.delim(input, header =TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
  head(Validation_call[,1:10])
  Validation_call = t(as.data.frame(Validation_call))
  Validation_call = Validation_call[5:length(Validation_call[,1]),]
  #head(Validation_call[,1:10])
}


## Feature selection CHANGE THE FEATURES

features <- c("V2185","V2214","V2215","V2211","V1657","V1679","V675","V1668","V674","V1680","V2207","V2224","V2224","V1673","V673","V1669","V1665","V2220","V693","V696","V855","V1678","V2212","V1664","V1670","V1671","V2225","V1674","V853","V1675","V850","V1663","V672","V680","V312","V673","V696","V773","V1570","V1657","V2185","V2593","V2663","V2733","V2185")
CorrelationAttribute <- c("V2185","V2214","V2215","V2211","V1657","V1679","V675","V1668","V674","V1680","V2207","V2224","V2224","V1673","V673","V1669","V1665","V2220","V693","V696","V855","V1678","V2212","V1664","V1670","V1671","V2225","V1674","V853","V1675","V850","V1663","V672","V680","V312","V673","V696","V773","V1570","V1657","V2185","V2593","V2663","V2733","V2185")
J48_previous <- c("V2733", "V1647","V1647","V1564","V624","V624","V1679","V1679","V2548","V148","V2548","V227","V2743")
filtering <- c("V312","V673","V696","V773","V1570","V1657","V2185","V25932","V2663","V2733")
Huge_combination <- c(J48_previous,CorrelationAttribute[1:25],filtering, features[1:20])
Huge_Unic <- unique(Huge_combination)

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
crossval = 5
total_HE = cross(HE,crossval)
total_HR = cross(HR,crossval)
total_Tr = cross(Tr,crossval)

sets <- vector("list", crossval)

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

### TRAINING MODELS ###

# making a list with n different models
models <- vector("list",crossval)

for (i in 1:crossval) {
  features = Huge_Unic[1:20]
  model = paste(features, collapse=' + ')
  model = paste("Subgroup ~ ", model, collapse ='')
  models[[i]] <- train(eval(parse(text=model)), data = sets[[i]][[1]], method = "nnet", verbose = FALSE)
}


"### COMPARE THE MODELS NOW ###
roc_average= c()
for (i in 1:crossval) {
  pred = predict(models[[i]],sets[[i]][[2]])
  print(confusionMatrix(sets[[i]][[2]][,1], pred)[2:3])
  if (1) {
    x = as.numeric(sets[[i]][[2]][,1])
    y = as.numeric(pred)
    roc = multiclass.roc(x,y)
    print(auc(roc))
    rs <- roc[['rocs']]
    plot.roc(rs[[1]])
    sapply(2:length(rs),function(h) lines.roc(rs[[h]],col=h))
    roc_average = append(roc_average,auc(roc))
  }
}
print(roc_average)
print(mean(roc_average))"
#################################### TRAIN FINAL MODEL ###########################################
#NNET_Model <- function(){
fitControl <- trainControl(## 10-fold CV
method = "repeatedcv",
number = 10,
  ## repeated ten times
repeats = 10)
features = Huge_Unic[1:12]
model = paste(features, collapse=' + ')
model = paste("Subgroup ~ ", model, collapse ='')
MODEL <- train(eval(parse(text=model)), data = combine, method = "nnet", trControl = fitControl, verbose = FALSE)
#}
save(MODEL, file= "model.pkl")
#################################### DO THE FINAL PREDICTIONS ####################################

#final <- read.delim("Validation_call.txt", header =TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
#final = t(as.data.frame(final))
Model <- NNET_Model(final_model)
Unlabelled_data <- Load_Unlabeled_Data("unlabelled_sample.txt")
pred = predict(MODEL,Unlabelled_data)
pred
