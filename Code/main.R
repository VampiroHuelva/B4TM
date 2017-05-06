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
source("code/misc.R")
load_depend()

# get the train data and write to WEKA file
combine = Load_labeled_Data()

################################## FEATURE SELECTION ####################################

# load: features, features_rfe, CorrelationAttribute, J48_previous. Huge_combination and
# Huge_Unic
load_features()

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


# train models for two class sets (HER2 and NOHER), (HR and NOHR) and (Triple and NOTriple)
features_HER2 = filter_Var_selection(HER2,number_of_features=30)
result_HER2 = train_all_models(HER2, features_HER2)

features_HR = filter_Var_selection(HR,number_of_features=30)
result_HR = train_all_models(HR, features_HR)

features_Triple = filter_Var_selection(Triple,number_of_features=30)
results_Triple = train_all_models(Triple, features_Triple)

summary(result_HER2)[3]$statistics$Accuracy
summary(result_HR)[3]$statistics$Accuracy
summary(results_Triple)[3]$statistics$Accuracy

### todo: analyse which features matter for which class

############################## TRAINING MODELS ON THREE CLASSES ##############################


accur = empty_accur()
pred = empty_pred()
models = empty_models()

# feature selection before training
for (j in 1:crossval) {
  
  # data and feature selection
  train_set = simple_cross_sets[[j]][[1]]
  test_set = simple_cross_sets[[j]][[2]]
  #features = filter_Var_selection(train_set,10)
  features = CorrelationAttribute
  
  # train the models and store on the trainset (feature selection before on whole set)
  models_trained = train_all_models(train_set,features)
  
  # saving all the models
  models = save_models(models, models_trained)
  
  # saving all the accuracies of the models
  accur = save_accur(accur, models_trained[[6]])
  
  # save predictions
  pred = save_pred(pred,models_trained, test_set) 
 
}


################# TRAINING MODELS ON THREE CLASSES FEATURE SELECTION ON TRAIN DATA ################


### select first all features from the feature selecting methods combined
features = CorrelationAttribute
accur = empty_accur()
pred = empty_pred()
models = empty_models()

# feature selection on each crossval
for (j in 1:crossval) {
  
  # get data
  train_set = simple_cross_sets[[j]][[1]]
  test_set = simple_cross_sets[[j]][[2]]

  # train the models (feature selection on train set)
  models_trained = train_all_models_with_feature_selection(train_set,features) 

  # saving all the models
  models = save_models(models, models_trained)
    
  # saving all the accuracies of the models
  accur = save_accur(accur, models_trained[[6]])
  
  # save predictions
  pred = save_pred(pred,models_trained, test_set)

}


##################### NOW we can compare the models with different feature selection ###################

acc_values = cbind(acc_gbmFit,acc_svmFit,acc_nnetFit,acc_mrFit,acc_rfFit)

accur
pred


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
unlist(final[1:10,1])



plot(acc_gbmFit)
plot(acc_gbmFit)
plot(acc_gbmFit)
plot(acc_gbmFit)
plot(acc_gbmFit)
plot(acc_values[,5])
######################################## TRAIN THE FINAL MODEL ON ALL DATA #############################



#################################### DO THE FINAL PREDICTIONS ##########################################


save(MODEL, file= "model.pkl")
Model <- NNET_Model(final_model)
Unlabelled_data <- Load_Unlabeled_Data("Data/unlabelled_sample.txt")
pred = predict(MODEL,Unlabelled_data)
pred
