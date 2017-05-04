#######################################################################################
########################################################################################
###
### Script to select, make and train a model for the classification of different
### cancer types. First the data is read in to perform a data inspection.... TODO FINISH  
### .....
### 
### Dependecies: misc.R
### input: ... datasets
### output: .... models
###
###
########################################################################################
########################################################################################

# clear environment and load source file
rm(list = ls())
source("code/misc.R")

# Loading all the dependicies
library(splines); library(parallel); library(survival); library(caret); library(mlbench)
library(gbm); library(corrplot); library(pROC); library(FSelector); library(qtlcharts)
library(lattice); library(energy); library(RWeka); 

# get the train data and write to WEKA file
combine = Load_labeled_Data()
#write.arff(combine, file = "combine.arff")


################################## FEATURE SELECTION ####################################


# Feature selection by importance ()
features = filter_Var_selection()

# Advance feature selction using WEKA. Algorithm used describe in the report
CorrelationAttribute <- c("V2185","V2214","V2215","V2211","V1657","V1679",
                          "V675","V1668","V674","V1680","V2207","V2224","V2224","V1673",
                          "V673","V1669","V1665","V2220","V693","V696","V855","V1678",
                          "V2212","V1664","V1670","V1671","V2225","V1674","V853","V1675",
                          "V850","V1663","V672","V680","V312","V673","V696","V773","V1570",
                          "V1657","V2185","V2593","V2663","V2733","V2185")
J48_previous <- c("V2733", "V1647","V1647","V1564","V624","V624","V1679","V1679","V2548",
                  "V148","V2548","V227","V2743")
filtering <- c("V312","V673","V696","V773","V1570","V1657","V2185","V25932","V2663","V2733")
Huge_combination <- c(J48_previous, CorrelationAttribute[1:25], filtering, features[1:20])
Huge_Unic <- unique(Huge_combination)


############################### GENERATING TEST/TRAIN SETS ###############################


# create simple cross validation sets (works only with 3 classes)
crossval = 10
simple_cross_sets = crossval_sets(combine, crossval) 

# create double cross validation sets
dcrosssets = double_crossval_datasets(combine, crossval=5, loops = 5)

# make from 3 classes sets 3 sets of two classes
two_class_data = twoclasses(combine)

# make the 2 classes data equal in classes sample numbers
HER2 = make_equal_sets(two_class_data[[1]])
HR = make_equal_sets(two_class_data[[2]])
Triple = make_equal_sets(two_class_data[[3]])

############################## TRAINING MODELS ############################################

## model selection (store the settings for models to train with)##
models = c("nnet", "rf", "gbm", "glm")

# make list for storing the trained models
trained_models <- vector("list",4)

## TRAIN: here i try the double loop crossval feature selection

# DO THE FOLOWWING FOR ALL FEATURE SELECTION METHODS
for (i in 1:length(models)) {
  # TODO specify for each model certain properties (number of features for instance)
  
  # make a list for the predictions on the test sets and for the real values in the test set
  pred_test = c()
  real_values = c()
  
  # feature selection on each crossval train and do test and store predictions for roc plot
  for (j in 1:crossval) {
    
    # do feature selectionA (TODO choose feature selection function)
    features = Huge_Unic[1:12]
    
    # set spefic settings for model
    fitControl <- trainControl(## 10-fold CV
      method = "repeatedcv",
      number = 10,
      ## repeated ten times
      repeats = 10) 
    
    # makeformula 
    formula_model = paste(features, collapse=' + ')
    formula_model = paste("Subgroup ~ ", model, collapse ='')
    
    model <- train(eval(parse(text=model)), data = trainset[i], method = models[i], trControl = fitControl, verbose = FALSE)
    # save model
    
    # add real and predicted values
    real = append(real, testset)      
    pred = append(pred(predict(model,testset[j])))
  }
  
  # save for each crossval the model, predictions and real classes
  trained_models[[i]] <- list(model, real, pred)
  
}


##################### NOW we can compare the models with different feature selection ###################

# plot rocplots of all models
for (i in 1:length(trained_models)) {
    r = c()
    print(confusionMatrix(real[i], pred[i])[2:3])
    x = as.numeric(sets[[i]][[1]][[j]][[2]][,1])
    y = as.numeric(pred)
    roc = multiclass.roc(x,y)
    r <- append(r,auc(roc))
    rs <- roc[['rocs']]
    plot.roc(rs[[1]])
    sapply(2:length(rs),function(h) lines.roc(rs[[h]],col=h))
}

####################### NOW we should train 4 models with best features selected on cross val ###########

# compare the 4 models with each other and do the final predictions (TODO)

#################################### DO THE FINAL PREDICTIONS ####################################

save(MODEL, file= "model.pkl")
Model <- NNET_Model(final_model)
Unlabelled_data <- Load_Unlabeled_Data("Data/unlabelled_sample.txt")
pred = predict(MODEL,Unlabelled_data)
pred





# lets train nnet, ft, svmRadial,J48 and gbm
# making a list with n different models

# ####  Getting scores of the models to choose the number of features to use and the best model
# Nmodel = 3
# models <- vector("list",Nmodel)
# for (i in 1:Nmodel) {
#   models[[i]] <- vector("list",loops)
#   for (j in 1:loops)
#     models[[i]][[j]] <- vector("list",crossval)
# }
# 
# for (i in 1:loops) {
#   for (j in 1:crossval) {
#      features =Huge_Unic[1:10+i]
#     model = paste(features, collapse=' + ')
#     model = paste("Subgroup ~ ", model, collapse ='')
#     model
#     fitControl <- trainControl(## 10-fold CV
#       method = "repeatedcv",
#       number = 10,
#       ## repeated ten times
#       repeats = 10)
#     models[[1]][[i]][[j]] <- train(eval(parse(text=model)), data = sets[[i]][[1]][[j]][[1]], method = "nnet", trControl = fitControl, verbose = FALSE)
#     
#     # gbm
#     features =Huge_Unic[1:15+i]
#     #[1:(5+j)]
#     model = paste(features, collapse=' + ')
#     model = paste("Subgroup ~ ", model, collapse ='')
#     models[[2]][[i]][[j]] <- train(eval(parse(text=model)), data = sets[[i]][[1]][[j]][[1]], method = "nnet", verbose = FALSE)
#     
#     #svm
#     features =Huge_Unic[1:20+i]
#     fitControl <- trainControl(## 10-fold CV
#       method = "repeatedcv",
#       number = 10,
#       ## repeated ten times
#       repeats = 10)
#     #features[1:(5+j)]
#     model = paste(features, collapse=' + ')
#     model = paste("Subgroup ~ ", model, collapse ='')
#     models[[3]][[i]][[j]] <- train(eval(parse(text=model)), data = sets[[i]][[1]][[j]][[1]], method = "nnet", trControl = fitControl, verbose = FALSE)
#   }    
# }
# ### COMPARE THE MODELS NOW ###
# 
# for (k in 1:Nmodel) {
#   for (i in 1:loops) {
#     r = c()#r[i]=
#     for (j in 1:crossval) {
#       #varimp <- varImp(models[[k]][[i]][[j]], varImp.train=FALSE)
#       pred = predict(models[[k]][[i]][[j]],data=sets[[i]][[1]][[j]][[2]],na.omit(sets[[i]][[1]][[j]][[2]]))
#       #print(confusionMatrix(sets[[i]][[1]][[j]][[2]][,1], pred)[2:3])
#       if (1) {
#         x = as.numeric(sets[[i]][[1]][[j]][[2]][,1])
#         y = as.numeric(pred)
#         roc = multiclass.roc(x,y)
#         r <- append(r,auc(roc))
#         rs <- roc[['rocs']]
#         #plot.roc(rs[[1]])
#         #sapply(2:length(rs),function(h) lines.roc(rs[[h]],col=h))
#         
#         
#       }
# #scoring of the dirent models selected with the number of features selected      
#     }
#     print(mean(r))
#     print(k)
#   }
# }