############################################################################################
###
### Script with all the models used in main.R
### dependicies: source("code/misc.R")
############################################################################################

# loading dependicies
source("code/misc.R")

# simple feature selection
filter_Var_selection = function(x,number_of_features) {
  
  fill=filterVarImp(x[,2:ncol(x)],x[,1],nonpara = FALSE)
  fill = sort(rowMeans(fill), decreasing = T)
  features <- names(fill[1:number_of_features])
  
  return(features)
  
}

feature_var_imp = function(model,number_of_features) {
  
  col_index <- varImp(model, scale = FALSE)$importance %>% 
    mutate(names=row.names(.)) %>%
    arrange(-Overall)
  imp_names <- col_index$names[1:number_of_features]
  
  return(imp_names)
  
}

feature_selection_rfe = function(x) {

  sets = crossval_sets(x,1)[[1]]
  train = sets[[1]]
  test = sets[[2]]
  
  ldaProfile <- rfe(train[,2:ncol(train)], train[,1],
                    sizes = c(1:10, 15, 30),
                    rfeControl = rfeControl(functions = ldaFuncs, method = "cv"))
  plot(ldaProfile, type = c("o", "g"))
  
  postResample(predict(ldaProfile, test[,2:ncol(test)]), test[,1])
  
  return("this should be a feature list")
  
}

filter_Var_selection = function(x,number_of_features) {
  
  fill=filterVarImp(x[,2:ncol(x)],x[,1],nonpara = FALSE)
  fill = sort(rowMeans(fill), decreasing = T)
  features <- names(fill[1:number_of_features])
  
  return(features)
  
}



