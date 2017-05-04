## functions used in B4TM

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
  return(combine)
}

# simple feature selection
filter_Var_selection = function() {
  fill=filterVarImp(combine[,2:ncol(combine)],combine[,1],nonpara = FALSE)
  fill = sort(rowMeans(fill), decreasing = T)
  features <- names(fill[1:30])
  return(features)
}

# make train and test sets for double loop cross validation
double_crossval_datasets = function(combine, crossval, loops) {
  
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
}

# make train and test sets for double loop cross validation (werkt alleen met 3 classes)
crossval_sets = function(combine, crossval) {
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
  return(sets)
}

# Make 3 sets all with two classes
twoclasses = function(x) {
  HER2 = x
  levels(HER2$Subgroup)[2]<-"NOHER2+"
  levels(HER2$Subgroup)[3]<-"NOHER2+"
  
  HR = x
  levels(HR$Subgroup)[1]<-'NOHR+'
  levels(HR$Subgroup)[3]<-'NOHR+'
  
  Triple = x
  levels(Triple$Subgroup)[1]<-'NOTriple Neg'
  levels(Triple$Subgroup)[2]<-'NOTriple Neg'
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