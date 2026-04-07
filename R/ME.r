library(ape);

# Reading distances
fileMatrix<- read.table("R/matrixME.txt", sep = ",", header = FALSE);
dist <- as.matrix(fileMatrix);
dist_obj <- as.dist(dist);
taxa <- readLines("Data/names.txt");
attr(dist_obj, "Labels") <- taxa

# Building tree
METree <- fastme.bal(dist_obj);
# Writing tree on file
write.tree(METree, file = "R/TreeME.tree");
