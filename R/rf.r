library(TreeDist);
library(ape);

tree1<-read.tree("R/tree1.tree");
tree2<-read.tree("R/tree2.tree");

distRF = JaccardRobinsonFoulds(tree1, tree2, normalize=TRUE);
write(distRF, "R/dist.txt");
