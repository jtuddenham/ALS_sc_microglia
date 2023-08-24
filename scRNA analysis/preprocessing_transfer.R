library(Seurat)
library(Matrix)
library(dplyr)
library(readxl)
library(SeuratWrappers)

source("clustering_diffex_functions.R")

##initial processing and merging of data##
allfiles = list.files("../raw_data/")
sampdata=read_excel("../intermediate_data/Table S1.xlsx", sheet = 1)
featuredata = read_excel("../intermediate_data/Table S1.xlsx", sheet = 3)

#set up sample lists for aggregation#
allfiles_hash = allfiles[allfiles %in% featuredata$`Sequencing ID`]
allfiles_nohash = allfiles[!(allfiles %in% allfiles_hash)]

#make full data objects#
dat_barcodes <- aggfilter_nonhash(allfiles_nohash, base_directory = "../raw_data/", umimin = 500, umimax = 10000, RP_filt = T, MT_filt = T)
alldat = dat_barcodes[[1]]; barcodes = dat_barcodes[[2]]
hash_counts <- aggfilter_hash(allfiles_hash, base_directory = "../raw_data/", umimin = 500, umimax = 10000, 
                              HTO_filter = 0, feature_sheet = "../intermediate_data/Table S1.xlsx",
                              RP_filt = T, MT_filt = T)
allhash = hash_counts[[1]]; counts = hash_counts[[2]]

#generate the combined seurat object with added metadata and finalized hash assignments#
old_new_seurat <- makeseurat_metadata(alldat = alldat, allhash = allhash, sampdata = sampdata, featuredata = featuredata)

#remove non-ALS samples hashed in - should probably do this beforehand anyways#
old_new_seurat$Label[old_new_seurat$Label == "Hashtag7" & old_new_seurat$orig.ident == "PM080"] = "negative"
old_new_seurat$Label[old_new_seurat$Label == "Hashtag6" | old_new_seurat$Label == "Hashtag7" & old_new_seurat$orig.ident == "PM077"] = "negative"
old_new_seurat = SetIdent(object = old_new_seurat, value = old_new_seurat$Label)
old_new_seurat <- subset(old_new_seurat, idents = c("negative", "uncertain", "doublet"), invert = TRUE)

###start with all cells - no integration here. This object is just a concatenation of the old data matrix and the new data (with no annotation)####
old_new_seurat=SCTransform(old_new_seurat)

####pairwise classification#####
###RF classification, old microglia clusters (1-9) only###
library(randomForest)
datmat=old_new_seurat@assays$SCT
allpreds_pair=list()
for (typ1 in 1:8) {
  for (typ2 in (typ1+1):9) {
    nam=paste0(typ1,"-",typ2)
    allpreds_pair[[nam]]=list()
    for (seedval in 1:20) {
      set.seed(seedval)
      trainset=c()
      for (ii in c(typ1,typ2)) {
        allset=which(old_new_seurat@meta.data$old_clust==ii)
        if (length(allset)>200) {
          trainset=c(trainset,sample(allset,200))
        } else {
          trainset=c(trainset,sample(allset,200,replace=T))
        }
      }
      trainclass=old_new_seurat@meta.data$old_clust[trainset]
      ###all pairwise genes###
      keepgen1=c()
      for (ii in c(typ1,typ2)) {
        inds1=which(trainclass==ii)
        inds2=setdiff(1:length(trainclass),inds1)
        diffgen=apply(datmat[,trainset],1,function(x){return(wilcox.test(x[inds1],x[inds2])$p.value)})
        tempgen=which(diffgen<0.01)
        keepgen1=unique(c(keepgen1,tempgen))
        print(c(ii,length(tempgen),length(keepgen1)))
      }
      #diffgen=apply(datmat[,trainset],1,function(x){return(kruskal.test(x,trainclass)$p.value)})
      #keepgen=which(diffgen<0.05)
      rfmod=randomForest(t(as.matrix(datmat[keepgen1,trainset])),y=as.factor(trainclass))
      rfpred=predict(rfmod,t(as.matrix(datmat[keepgen1,])),type = "prob")
      allpreds_pair[[nam]][[seedval]]=rfpred
      print(c(nam,seedval))
    }
  }
  save(allpreds_pair,file="rf_20fold_prediction_allcells_pairwise.rda")
}
load("rf_20fold_prediction_allcells_pairwise.rda")

###Identify where query cells are classified into, across the 20 runs###
maxval=matrix(0,nrow=nrow(allpreds_pair[[1]][[1]]),ncol=length(allpreds_pair))
for (ii in 1:length(allpreds_pair)) {
  allmat=matrix(0,nrow=nrow(allpreds_pair[[ii]][[1]]),ncol=length(allpreds_pair[[ii]]))
  for (jj in 1:length(allpreds_pair[[ii]])) {
    allmat[,jj]=apply(allpreds_pair[[ii]][[jj]],1,which.max)
  }
  maxval[,ii]=apply(allmat,1,function(x){ttt=table(x);return(names(ttt)[which.max(ttt)])})
  pairs=strsplit(names(allpreds_pair)[ii],"-")[[1]]
  maxval[,ii]=pairs[as.numeric(maxval[,ii])]
  print(ii)
}
colnames(maxval)=names(allpreds_pair)
save(maxval,file="rf_20fold_prediction_allcells_pairwise_maxval.rda")
load("rf_20fold_prediction_allcells_pairwise_maxval.rda")

maxtab=apply(maxval,1,function(x){ttt=table(x);return(names(ttt)[which.max(ttt)])})
maxtabval=apply(maxval,1,function(x){ttt=table(x);return(max(ttt))})
maxtab[which(old_new_seurat$old_clust %in% 1:9)]=old_new_seurat@meta.data$old_clust[oldvals]
write.csv(maxtab,file="cluster_assignments_mapped_microglia_only.csv")

####Integration with Harmony for visualization purposes, and to assign non-microglial cells
library(harmony)
old_new_seurat=SCTransform(old_new_seurat)
old_new_seurat=RunPCA(old_new_seurat)
old_new_seurat@meta.data$new_clust=maxtab2
thetaval=2
dimval=30
old_new_seurat=RunHarmony(old_new_seurat,assay.use = "SCT",theta=thetaval,group.by.vars="orig.ident")
old_new_seurat=RunUMAP(old_new_seurat,reduction="harmony",dims=1:dimval)
old_new_seurat=RunTSNE(old_new_seurat,reduction="harmony",dims=1:dimval)
old_new_seurat=FindNeighbors(old_new_seurat,reduction="harmony",dims=1:dimval)
old_new_seurat=FindClusters(old_new_seurat,resolution=0.3)

old_new_seurat@meta.data$old_clust=as.numeric(as.character(old_new_seurat@meta.data$old_clust))
old_new_seurat=NormalizeData(old_new_seurat)

###clusters: 7= Mono, 5 = T-cell, 11 = RBC, 16 = B-cell, 6+8+15 = Astro/Neuron
filtered_cluster=old_new_seurat@meta.data$new_clust
filtered_cluster[old_new_seurat@meta.data$seurat_clusters==7]=10
filtered_cluster[old_new_seurat@meta.data$seurat_clusters==5]=11
filtered_cluster[old_new_seurat@meta.data$seurat_clusters==16]=12
filtered_cluster[old_new_seurat@meta.data$seurat_clusters %in% c(6,8,15)]=13
filtered_cluster[old_new_seurat@meta.data$seurat_clusters==11]=14
filtered_cluster[old_new_seurat@meta.data$seurat_clusters==10]=15
filtered_cluster[old_new_seurat@meta.data$seurat_clusters==13]=16
filtered_cluster[old_new_seurat@meta.data$seurat_clusters==9 & old_new_seurat@meta.data$new_clust!=6]=5
old_new_seurat@meta.data$filtered_cluster=as.numeric(as.character(filtered_cluster))
save(old_new_seurat,file="old_new_seurat_mapped_clusters.rda")

##build and plot phylotree of microglial clusters##
#use only the microglial clusters for this#
old_new_seurat <- subset(old_new_seurat, idents = c(1:9))

old_new_seurat <- BuildClusterTree(old_new_seurat,
                                   dims = 1:30, 
                                   reorder = F, 
                                   reorder.numeric = F, 
                                   verbose = T)
PlotClusterTree(old_new_seurat)