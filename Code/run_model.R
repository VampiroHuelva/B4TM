# Author: Chao (Cico) Zhang
# Date: 31 Mar 2017
# Usage: Rscript run_model.R -i unlabelled_sample.txt -m model.pkl -o output.txt
# If you are using python, please use the Python script template instead.
# Set up R error handling to go to stderr
options(show.error.messages=F, error=function(){cat(geterrmessage(),file=stderr());q("no",1,F)})

# Since the nature of the model, the accuaraty I sent the deadline can be slightly different to the model built to handly in run_model.R. Each run it change
# This script have be made to be loadin the windows 10 home terminal.
# Import required libraries
# You might need to load other packages here. This are the libreries requiered for all the code, scoring of the models... For this script are not neccesary
suppressPackageStartupMessages({
  library('getopt')
  library('caret')
  library(splines)
  library(parallel)
  library(survival)
  library(caret)
  library(mlbench)
  library(gbm)
  library(corrplot)
  library(pROC)
  library(FSelector)
  library(qtlcharts)
  library(lattice)
  library(energy)
  library(RWeka)
})

# Take in trailing command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Get options using the spec as defined by the enclosed list
# Read the options from the default: commandArgs(TRUE)
option_specification <- matrix(c(
  'input', 'i', 2, 'character',
  'model', 'm', 2, 'character',
  'output', 'o', 2, 'character'
), byrow=TRUE, ncol=4);

# Parse options
options <- getopt(option_specification);

# Start your coding

# suggested steps
# Run the script with the function for loading the data
source("misc.R")
# Step 1: load unlabeled data with Load_Unlabeled_Data() function from Building_Final_mode.R script
Unlabelled_data <- Load_Unlabeled_Data(options$input)

# Step 2: apply the model to the input file (options$input) to do the prediction
load(options$model)
Prediction <- predict(MODEL,Unlabelled_data)
Prediction
# Step 3: write the prediction into the desinated output file (options$output)
Pred_as_charact <- as.character(Prediction)
output <- file(options$output)
writeLines(Pred_as_charact,output)
close(output)

# End your coding
message ("Done!")






