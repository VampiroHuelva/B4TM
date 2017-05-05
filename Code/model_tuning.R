############################################################################################
###
### Script with all the models used in main.R
###
############################################################################################


# gbm https://topepo.github.io/caret/model-training-and-tuning.html#model-training-and-parameter-tuning
# https://cran.r-project.org/web/packages/gbm/gbm.pdf
# 
gbm_tuning = function(x,features) {
  gbmGrid <-  expand.grid(interaction.depth = c(1, 5, 9), 
                          n.trees = (1:30)*50, 
                          shrinkage = 0.1,
                          n.minobsinnode = 20)
  
  fitControl <- trainControl(method = "repeatedcv",
                             number = 10,
                             repeats = 10,
                             ## Estimate class probabilities
                             classProbs = TRUE)

  formula_model = paste(features, collapse=' + ')
  formula_model = paste("Subgroup ~ ", formula_model, collapse ='')
  
  gbmFit <- train(eval(parse(text=formula_model)), data = x, 
                   method = "gbm", 
                   trControl = fitControl, 
                   verbose = FALSE, 
                   tuneGrid = gbmGrid)
  return(gbmFit)
}

# svm
svm_tuning = function(x,features) {
  
  fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10,
                           ## Estimate class probabilities
                           classProbs = TRUE)
  
  formula_model = paste(features, collapse=' + ')
  formula_model = paste("Subgroup ~ ", formula_model, collapse ='')

  svmFit <- train(eval(parse(text=formula_model)), data = x, 
                method = "svmRadial", 
                trControl = fitControl, 
                preProc = c("center", "scale"),
                tuneLength = 8)
  return(svmFit)
}

# nnet Neural Networks with Feature Extraction
# Tuning parameters:
# size (#Hidden Units)
# decay (Weight Decay)
nnet_tuning = function(x,features) {
  
  fitControl <- trainControl(method = "repeatedcv",
                             number = 10,
                             repeats = 10,
                             ## Estimate class probabilities
                             classProbs = TRUE)
  
  formula_model = paste(features, collapse=' + ')
  formula_model = paste("Subgroup ~ ", formula_model, collapse ='')
  
  nnetFIT <- train(eval(parse(text=formula_model)), data = x, 
                  method = "pcaNNet", 
                  trControl = fitControl, 
                  preProc = c("center", "scale"))
  return(nnetFIT)
}

# Penalized Multinomial Regression https://topepo.github.io/caret/train-models-by-tag.html#l2-regularization
# Tuning parameters:
# decay (Weight Decay)
mr_tuning = function(x,features) {
  
  fitControl <- trainControl(method = "repeatedcv",
                             number = 10,
                             repeats = 10,
                             ## Estimate class probabilities
                             classProbs = TRUE)
  
  formula_model = paste(features, collapse=' + ')
  formula_model = paste("Subgroup ~ ", formula_model, collapse ='')
  
  mrFIT <- train(eval(parse(text=formula_model)), data = x, 
                 method = "multinom", 
                 trControl = fitControl)
  return(mrFIT)
}

# rt  
# No need to perform feature selection with random forest, it will do so automatically,
# ignoring features that make bad split points.
rf_tuning = function(x,features) {
  
  fitControl <- trainControl(method = "repeatedcv",
                             number = 10,
                             repeats = 10,
                             ## Estimate class probabilities
                             classProbs = TRUE)
  
  formula_model = paste(features, collapse=' + ')
  formula_model = paste("Subgroup ~ ", formula_model, collapse ='')
  
  rfFIT <- train(eval(parse(text=formula_model)), data = x, 
                 method = "rf", 
                 trControl = fitControl)
  return(rfFIT)
}

