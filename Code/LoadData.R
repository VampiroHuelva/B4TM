# Preprocessing the data

Load_Unlabeled_Data <- function(input){
  
  Validation_call <- read.delim(input, header =TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
  head(Validation_call[,1:10])
  Validation_call = t(as.data.frame(Validation_call))
  Validation_call = Validation_call[5:length(Validation_call[,1]),]
  
}