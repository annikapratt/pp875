library(ape)
library(phangorn)
library(phytools)
library(ggplot2)

getwd() #Check the working directory. we want to be in the results/RAxML/10Mb-concatenation folder
#setwd("PathTo/results/RAxML/10Mb-concatenation") #replace PathTo with the correct path and run if not in the correct folder

tree_files <-list.files(pattern="\\.raxml.bestTree$") #List all .bestTree files. $ ensures the end of the name

trees<- list() # list with all the trees
class(trees)<- "multiPhylo" #make it a multiphylo object for ease of use with other 

i<-1
for(tree_file in tree_files){ ##go thru each file and read the tree
  trees[[i]]<- read.tree(tree_file)
  i<-i+1
}

#re-reroot all our gene trees by the respective outgroup
for(i in 1:length(trees)){
  trees[[i]]<- root(trees[[i]],
                         outgroup = "H_vulgare_HVens23",
                         resolve.root=TRUE)
  trees[[i]]<-chronos(trees[[i]]) ## make ultrametric for nicer densitree
}

st<-superTree(trees)
st<-root(st,"H_vulgare_HVens23",resolve.root = T)
plot(st)

densiTree(trees,consensus=st,scaleX=T,type='cladogram', alpha=0.1)

trees2 <- read.nexus("../../../data/Wheat_Relative_History_Data_Glemin_et_al/Densitree_OneCopyGenes-modified.nex")

#re-reroot all our gene trees by the respective outgroup
for(i in 1:length(trees2)){
  trees2[[i]]<- root(trees2[[i]],
                         outgroup = "H_vulgare_HVens23",
                         resolve.root=TRUE)
  trees2[[i]]<-chronos(trees2[[i]]) ## make ultrametric for nicer densitree
}

st2<-superTree(trees2)
st2<-root(st2,"H_vulgare_HVens23",resolve.root = T)
plot(st2)

densiTree(trees2,consensus=st2,scaleX=T,type='cladogram', alpha=0.1)
