library(ape)
library(phangorn)

setwd("../results/")

tre = read.tree(file="07-supertree.tre")
plot(tre)

rtre = root(tre,outgroup="H_vulgare_HVens23", resolve.root=TRUE)
plot(rtre)

gene_trees <- read.tree("04-all_gene_trees.tre")
st_parsimony<-superTree(gene_trees)
st_parsimony<-root(st,"H_vulgare_HVens23",resolve.root = T)
plot(st_parsimony)
