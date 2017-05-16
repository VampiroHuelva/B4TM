library(ggplot2)
library("cluster")
library("factoextra")
setwd("/home/bjorn/Documents/school/TM/B4TM/Code")

# Processing the data
cgh_data <- read.delim("Train_call.txt", header =TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
drop_columns <- c("Chromosome","Start","End","Nclone")
cgh_data <- cgh_data[ , !(names(cgh_data) %in% drop_columns)]
cgh_data <- t(as.matrix(sapply(cgh_data, as.numeric)))
clinical <- read.delim("Train_clinical.txt", header =TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
rownames(cgh_data) <- paste(rownames(cgh_data), clinical[[2]], sep=' ')
# rownames(cgh_data) <- cgh_data[,1] 
# rownames(cgh_data) <- clinical[[2]]

# cgh_data <- cbind(id=rownames(cgh_data), cgh_data)
test <- cgh_data
test <- scale(test)
test[is.nan(test)] <- 0

desc_stats <- data.frame(
  Min = apply(test, 2, min), # minimum
  Med = apply(test, 2, median), # median
  Mean = apply(test, 2, mean), # mean
  SD = apply(test, 2, sd), # Standard deviation
  Max = apply(test, 2, max) # Maximum
)


# Clustering
res.dist <- get_dist(test, stand = FALSE, method = "pearson")
fviz_dist(res.dist, 
          gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

km.res <- kmeans(test, 10, nstart = 1)
plot(test, col = km.res$cluster, pch = 19, frame = FALSE,
     main = "K-means with k = 3")
points(km.res$centers, col = 1:2, pch = 8, cex = 3)

fviz_nbclust(test, kmeans, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2)

fviz_cluster(km.res, data = test)

# 2. Compute dissimilarity matrix
d <- dist(test, method = "euclidean")
# Hierarchical clustering using Ward's method
res.hc <- hclust(d, method = "ward.D2" )
# Cut tree into 4 groups
grp <- cutree(res.hc, k = 10)
# Visualize
plot(res.hc, cex = 0.6) # plot tree
rect.hclust(res.hc, k = 3, border = 2:5) # add rectangle

# Compute hierarchical clustering and cut into 4 clusters
res <- hcut(test, k = 3, stand = TRUE)
# Visualize
fviz_dend(res, rect = TRUE, cex = 0.5,
          k_colors = c("#00AFBB","#2E9FDF", "#E7B800"))

