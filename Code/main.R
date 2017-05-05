#######################################################################################
########################################################################################
###
### Script to select, make and train a model for the classification of different
### cancer types. First the data is read in to perform a data inspection.... TODO FINISH  
### .....
### 
### Dependecies: misc.R and model_tuning.R
### input: Train_call.txt, Train_clinical.txt andunlabelled_sample.txt
### output: .... models
###
###
########################################################################################
########################################################################################

# clear environment and load source files
rm(list = ls())
source("code/misc.R")
source("code/model_tuning.R")
source("code/feature_selection.R")

# Loading all the dependicies
library(splines); library(parallel); library(survival); library(caret); library(mlbench)
library(gbm); library(corrplot); library(pROC); library(FSelector); library(qtlcharts)
library(lattice); library(energy); library(RWeka); library(obliqueRF); library(stepPlr);
library(dplyr); library(LiblineaR); library(kernlab)

# get the train data and write to WEKA file
combine = Load_labeled_Data()
#write.arff(combine, file = "combine.arff")


################################## FEATURE SELECTION ####################################


# Feature selection by importance ()
features = filter_Var_selection()[1:10]

# use rfe featureselection (not returningfeatures still)
features_rfe = feature_selection_rfe(combine)

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


############################## TRAINING MODELS ON TWO CLASSES ##############################


# see in source file model_tuning.R for link to documentation
result_all_data = train_all_models(combine,features)

# see what happenson two class sets (HER2, HR and Triple)
features_HER2 = filter_Var_selection(HER2,number_of_features=30)
result_HER2 = train_all_models(HER2, features_HER2)

features_HR = filter_Var_selection(HR,number_of_features=30)
result_HR = train_all_models(HR, features_HR)

features_Triple = filter_Var_selection(Triple,number_of_features=30)
results_Triple = train_all_models(Triple, features_Triple)

summary(result_HER2)[3]$statistics$Accuracy
summary(result_HR)[3]$statistics$Accuracy
summary(results_Triple)[3]$statistics$Accuracy


############################## TRAINING MODELS ON THREE CLASSES ##############################


real = c()
pred_gbmFit = c()
pred_svmFit = c()
pred_nnetFit = c()
pred_mrFit = c()
pred_rfFit = c()
real_values = c()

# feature selection on each crossval train and do test and store predictions for roc plot
for (j in 1:crossval) {
  
  # data and feature selection
  train_set = simple_cross_sets[[j]][[1]]
  test_set = simple_cross_sets[[j]][[2]]
  #features = filter_Var_selection(train_set,10)
  features = CorrelationAttribute
  
  # train the models and store on the trainset (feature selection before on whole set)
  models_trained = get_trained_models(train_set,features) 
  gbmFit = models_trained[[1]]
  svmFit = models_trained[[2]]
  nnetFit = models_trained[[3]]
  mrFit = models_trained[[4]]
  rfFit = models_trained[[5]]
  
  # add real and predicted values
  real = append(real, test_set$Subgroup)      
  pred_gbmFit = append(pred_gbmFit, predict(gbmFit,test_set))
  pred_svmFit = append(pred_svmFit, predict(svmFit,test_set))
  pred_nnetFit = append(pred_nnetFit, predict(nnetFit,test_set))
  pred_mrFit = append(pred_mrFit, predict(mrFit,test_set))
  pred_rfFit = append(pred_rfFit, predict(rfFit,test_set))
}


####################################### USING THE ABOVE RESULTS ##################################

## TODO MAKE SOME MORE STATISTICAL ESTIMATES POSSIBLY BY GETTING RESULTS FOR EACH CROSS VALLIDATION
## SO DO 10 times EVERY CROSS VALIDATION: summary(FITMODEL)[3]$statistics$Accuracy 
## USE THIS DATA TO MAKE A CONFIDENCE INTERVAL OR REPORT FOR EACH ALL THE CONVIDENCE INTERVALS

sum(real == pred_nnetFit)
sum(real == pred_mrFit)
sum(real == pred_svmFit)
sum(real == pred_rfFit)
sum(real == pred_gbmFit)
real==2

class_cancer = 1
sum(real==class_cancer)
sum(real[real==class_cancer] == pred_nnetFit[real==class_cancer])
sum(real[real==class_cancer] == pred_mrFit[real==class_cancer])
sum(real[real==class_cancer] == pred_svmFit[real==class_cancer])
sum(real[real==class_cancer] == pred_rfFit[real==class_cancer])
sum(real[real==class_cancer] == pred_gbmFit[real==class_cancer])

sum(real=class_cancer)
class_cancer = 2
sum(real==class_cancer)
sum(real[real==class_cancer] == pred_nnetFit[real==class_cancer])
sum(real[real==class_cancer] == pred_mrFit[real==class_cancer])
sum(real[real==class_cancer] == pred_svmFit[real==class_cancer])
sum(real[real==class_cancer] == pred_rfFit[real==class_cancer])
sum(real[real==class_cancer] == pred_gbmFit[real==class_cancer])

class_cancer = 3
sum(real==class_cancer)
sum(real[real==class_cancer] == pred_nnetFit[real==class_cancer])
sum(real[real==class_cancer] == pred_mrFit[real==class_cancer])
sum(real[real==class_cancer] == pred_svmFit[real==class_cancer])
sum(real[real==class_cancer] == pred_rfFit[real==class_cancer])
sum(real[real==class_cancer] == pred_gbmFit[real==class_cancer])


################################### EXPERIMENTAL STUFF ###########################################
################################### EXPERIMENTAL STUFF ###########################################
################################### EXPERIMENTAL STUFF ###########################################
################################### EXPERIMENTAL STUFF ###########################################
################################### EXPERIMENTAL STUFF ###########################################
################################### EXPERIMENTAL STUFF ###########################################

################# TRAINING MODELS ON THREE CLASSES FEATURE SELECTIONON TRAIN DATA ################


features = CorrelationAttribute
real = c()
pred_gbmFit = c()
pred_svmFit = c()
pred_nnetFit = c()
pred_mrFit = c()
pred_rfFit = c()
real_values = c()

# feature selection on each crossval train and do test and store predictions for roc plot
for (j in 1:crossval) {
  
  # get data
  train_set = simple_cross_sets[[j]][[1]]
  test_set = simple_cross_sets[[j]][[2]]

  # train the models and store on the trainset (feature selection on train set)
  gbmFit = gbm_tuning(train_set,features)
  gbmFit_features = feature_var_imp(gbmFit,10)
  gbmFit = gbm_tuning(train_set,gbmFit_features)
  
  svmFit = svm_tuning(train_set,features)
  svmFit_features = feature_var_imp(svmFit,10)
  svmFit = svm_tuning(train_set,svmFit_features)

  nnetFit = nnet_tuning(train_set,features)
  nnetFit_features = feature_var_imp(nnetFit,10)
  nnetFit = nnet_tuning(train_set,nnetFit_features)
 
  mrFit = mr_tuning(train_set,features)
  mrFit_features = feature_var_imp(mrFit,10)
  mrFit = mr_tuning(train_set,mrFit_features) 
  
  rfFit = rf_tuning(train_set,features)
  rfFit_features = feature_var_imp(rfFit,10)
  rfFit = rf_tuning(train_set,rfFit_features) 
  
  # add real and predicted values
  real = append(real, test_set$Subgroup)      
  pred_gbmFit = append(pred_gbmFit, predict(gbmFit,test_set))
  pred_svmFit = append(pred_svmFit, predict(svmFit,test_set))
  pred_nnetFit = append(pred_nnetFit, predict(nnetFit,test_set))
  pred_mrFit = append(pred_mrFit, predict(mrFit,test_set))
  pred_rfFit = append(pred_rfFit, predict(rfFit,test_set))
}


##################### NOW we can compare the models with different feature selection ###################


# further stuff doesnt work still

# plot rocplots of all models

# real = c(simple_cross_sets[[1]][[2]]$Subgroup,simple_cross_sets[[2]][[2]]$Subgroup,
#          simple_cross_sets[[3]][[2]]$Subgroup,simple_cross_sets[[4]][[2]]$Subgroup,
#          simple_cross_sets[[5]][[2]]$Subgroup,simple_cross_sets[[6]][[2]]$Subgroup,
#          simple_cross_sets[[7]][[2]]$Subgroup,simple_cross_sets[[8]][[2]]$Subgroup,
#          simple_cross_sets[[9]][[2]]$Subgroup,simple_cross_sets[[10]][[2]]$Subgroup)

multiclass.roc(as.numeric(real), as.numeric(pred_mrFit))

multiclass.roc(real, pred_mrFit)[1]plot_roc(real, pred_svmFit)
multiclass.roc(real, pred_mrFit)[1]plot_roc(real, pred_nnetFit)

plot_roc(as.numeric(real), as.numeric(pred_mrFit))

plot_roc = function(real, pred) {
  print(confusionMatrix(real, pred)[2:3])
  test_all = as.numeric(real)
  pred_num = as.numeric(pred)
  roc = multiclass.roc(pred_num,test_all)
  auc = auc(roc)
  rs <- roc[['rocs']]
  plot.roc(rs[[1]])
  sapply(2:length(rs),function(h) lines.roc(rs[[h]],col=h))
  print(auc)
}


####################### NOW we should train 4 models with best features selected on cross val ###########

# compare the 4 models with each other and do the final predictions (TODO)

# EXAMPLE OFTWO MODELS (gbmFit3 settings used above)

gbmFit3 <- train(Subgroup ~ ., data = combine, 
                 method = "gbm", 
                 trControl = fitControl, 
                 verbose = FALSE, 
                 tuneGrid = gbmGrid,
                 ## Specify which metric to optimize
                 metric = "ROC")

features = Huge_Unic[1:12]
formula_model = paste(features, collapse=' + ')
formula_model = paste("Subgroup ~ ", formula_model, collapse ='')

svmFit <- train(eval(parse(text=formula_model)), data = combine, 
                method = "svmRadial", 
                trControl = fitControl, 
                preProc = c("center", "scale"),
                tuneLength = 8,
                metric = "ROC")
svmFit
resamps <- resamples(list(GBM = gbmFit3,SVM = svmFit))
                     
summary(resamps)

#################################### DO THE FINAL PREDICTIONS ####################################

save(MODEL, file= "model.pkl")
Model <- NNET_Model(final_model)
Unlabelled_data <- Load_Unlabeled_Data("Data/unlabelled_sample.txt")
pred = predict(MODEL,Unlabelled_data)
pred
