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
features = filter_Var_selection(combine,30)

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
acc_gbmFit = c()
acc_svmFit = c()
acc_nnetFit = c()
acc_mrFit = c()
acc_rfFit = c()
acc_values = c()

# feature selection on each crossval train and do test and store predictions for roc plot
for (j in 1:crossval) {
  
  # data and feature selection
  train_set = simple_cross_sets[[j]][[1]]
  test_set = simple_cross_sets[[j]][[2]]
  #features = filter_Var_selection(train_set,10)
  features = CorrelationAttribute
  
  # train the models and store on the trainset (feature selection before on whole set)
  models_trained = train_all_models(train_set,features) 
  gbmFit = models_trained[[1]]
  svmFit = models_trained[[2]]
  nnetFit = models_trained[[3]]
  mrFit = models_trained[[4]]
  rfFit = models_trained[[5]]
  resample = models_trained[[6]]

  acc_gbmFit = append(acc_gbmFit,summary(resample)$statistics$Accuracy[1,4])
  acc_svmFit = append(acc_svmFit,summary(resample)$statistics$Accuracy[2,4])
  acc_nnetFit = append(acc_nnetFit,summary(resample)$statistics$Accuracy[3,4])
  acc_mrFit = append(acc_mrFit,summary(resample)$statistics$Accuracy[4,4])
  acc_rfFit = append(acc_rfFit,summary(resample)$statistics$Accuracy[5,4])

  # add real and predicted values
  real = append(real, test_set$Subgroup)      
  pred_gbmFit = append(pred_gbmFit, predict(gbmFit,test_set))
  pred_svmFit = append(pred_svmFit, predict(svmFit,test_set))
  pred_nnetFit = append(pred_nnetFit, predict(nnetFit,test_set))
  pred_mrFit = append(pred_mrFit, predict(mrFit,test_set))
  pred_rfFit = append(pred_rfFit, predict(rfFit,test_set))
}


################# TRAINING MODELS ON THREE CLASSES FEATURE SELECTION ON TRAIN DATA ################


### select first all features from the feature selecting methods combined
features = CorrelationAttribute

real = c()
pred_gbmFit = c()
pred_svmFit = c()
pred_nnetFit = c()
pred_mrFit = c()
pred_rfFit = c()
real_values = c()
acc_gbmFit = c()
acc_svmFit = c()
acc_nnetFit = c()
acc_mrFit = c()
acc_rfFit = c()
acc_values = c()

# feature selection on each crossval train and do test and store predictions for roc plot
for (j in 1:crossval) {
  
  # get data
  train_set = simple_cross_sets[[j]][[1]]
  test_set = simple_cross_sets[[j]][[2]]

  # train the models (feature selection on train set)
  print("gbm")
  # do some feature selection that is woring???
  gbmFit = gbm_tuning(train_set,features)

  # do some feature selection that is woring???
  print("svm")
  svmFit = svm_tuning(train_set,features)

  # do some feature selection that is woring???
  print("nnet")
  nnetFit = nnet_tuning(train_set,features)

  print("mr")
  mrFit = mr_tuning(train_set,features)
  mrFit_features = feature_var_imp(mrFit,10)
  mrFit = mr_tuning(train_set,mrFit_features) 
  
  print("rf")
  rfFit = rf_tuning(train_set,features)
  rfFit_features = feature_var_imp(rfFit,10)
  rfFit = rf_tuning(train_set,rfFit_features) 
  
  # saving all the accuracies of the models
  resamps <- resamples(list(GBM = gbmFit, SVM = svmFit, NNET = nnetFit, MR = mrFit, RF = rfFit))
  acc_gbmFit = append(acc_gbmFit,summary(resample)$statistics$Accuracy[1,4])
  acc_svmFit = append(acc_svmFit,summary(resample)$statistics$Accuracy[2,4])
  acc_nnetFit = append(acc_nnetFit,summary(resample)$statistics$Accuracy[3,4])
  acc_mrFit = append(acc_mrFit,summary(resample)$statistics$Accuracy[4,4])
  acc_rfFit = append(acc_rfFit,summary(resample)$statistics$Accuracy[5,4])

    # add real and predicted values
  real = append(real, test_set$Subgroup)      
  pred_gbmFit = append(pred_gbmFit, predict(gbmFit,test_set))
  pred_svmFit = append(pred_svmFit, predict(svmFit,test_set))
  pred_nnetFit = append(pred_nnetFit, predict(nnetFit,test_set))
  pred_mrFit = append(pred_mrFit, predict(mrFit,test_set))
  pred_rfFit = append(pred_rfFit, predict(rfFit,test_set))
}


##################### NOW we can compare the models with different feature selection ###################

acc_values = cbind(acc_gbmFit,acc_svmFit,acc_nnetFit,acc_mrFit,acc_rfFit)
 
total_HER2 = sum(real==1)
total_HR = sum(real==2)
total_TRIPLE = sum(real==3)

gbm = results_prediction(real, pred_gbmFit)
svm = results_prediction(real, pred_svmFit)
nnet = results_prediction(real, pred_nnetFit)
mr = results_prediction(real, pred_mrFit)
rf = results_prediction(real, pred_rfFit)

results = cbind(gbm,svm,nnet,mr,rf)

final = rbind(acc_values,results)
final



######################################## TRAIN THE FINAL MODEL ON ALL DATA #############################



#################################### DO THE FINAL PREDICTIONS ##########################################


save(MODEL, file= "model.pkl")
Model <- NNET_Model(final_model)
Unlabelled_data <- Load_Unlabeled_Data("Data/unlabelled_sample.txt")
pred = predict(MODEL,Unlabelled_data)
pred
