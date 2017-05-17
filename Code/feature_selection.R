############################################################################################
###
### Script with all the models used in main.R
### dependicies: source("code/misc.R")
############################################################################################

# loading dependicies
source("code/misc.R")


# Feature selection by importance ()
features = filter_Var_selection(combine,30)

# use rfe featureselection (not returningfeatures still)
#features_rfe = feature_selection_rfe(combine)

# Advance feature selction using WEKA. Algorithm used describe in the report
CorrelationAttribute <- c("V2185","V2214","V2215","V2211","V1657","V1679",
                          "V675","V1668","V674","V1680","V2207","V2224","V2224","V1673",
                          "V673","V1669","V1665","V2220","V693","V696","V855","V1678",
                          "V2212","V1664","V1670","V1671","V2225","V1674","V853","V1675",
                          "V850","V1663","V672","V680","V312","V673","V696","V773","V1570",
                          "V1657","V2185","V2593","V2663","V2733","V2185")
J48_previous <- c("V2733", "V1647","V1647","V1564","V624","V624","V1679","V1679","V2548",
                  "V148","V2548","V227","V2743")
filtering <- c("V312","V673","V696","V773","V1570","V1657","V2185","V2593","V2663","V2733")
Huge_combination <- c(J48_previous, CorrelationAttribute, filtering, features)
Huge_Unic <- unique(Huge_combination)


##### Combining and arranging classses ####

### HR~TRIPLE
HRTR_f <- c("V59","V68","V313","V386","V624","V625","V692","V696", "V736","V773","V842","V843", "V1567","V1568","V1640","V1647","V1657","V1658","V1679","V1688","V1839","V1974","V2018","V2026","V2548",
            "V2026","V2548", "V2593","V2733",  "V2734", "V2735","V2740","V2749", "V2752","V2808")

# J48 in selected features. 10 fold cross-validation. 0.76 ac. 
HRTRJ48_f <- c("V2733","V1647","V1567","V842","V2548","V68")
HRTRPEAR_f <- c("V1679", "V1657", "V1673", "V1680", "V1669", "V1668", "V1678", "V1665", "V675")
HRTR_features = unique(c(HRTR_f,HRTRJ48_f, HRTRPEAR_f))

### HER+~TRIPLE-    1 acuarracy
HERTR_f <- c("V17","V33","V55","V485","V731","V905","V1001","V1671","V2026","V2105","V2114","V2126","V2167","V2185","V2220","V2248","V2609","V2654","V2774")
HERTRJ48_f <- c("V2185")
HERTR_features = unique(c(HERTR_f, HERTRJ48_f))

###  HER+~HR       1 acuarracy
HERHR_f <- c("V118","V475","V673","V853","V1092","V1574","V2169","V2184","V2185","V2214","V2277","V2380","V2663")
HERHRJ48_f <- c("V2185")
HERHRJ48_features = c(HERHR_f, HERHRJ48_f)

### HER+HR~TRIPLE-
# Best accuaracy with all the feature than with the feature selection. 0.88
HERHR_TRIPLE_J48 <- c("V2185","V2026","V1657","V736","V1688","V306","V736","V59","V487","V2762")



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



