library(Seurat)
library(Matrix)
library(dplyr)
library(readxl)
library(SeuratWrappers)
library(DropletUtils)
library(clusterProfiler)
library(ReactomePA)

###load in the data with labels transferred###
load("old_new_seurat_mapped_clusters.rda")
old_new_seurat$filtered_cluster = factor(old_new_seurat$filtered_cluster, levels = c(1:16))
old_new_seurat = SetIdent(old_new_seurat, value = old_new_seurat$filtered_cluster)

###set assay to RNA and normalize###
DefaultAssay(old_new_seurat) = "RNA"
old_new_seurat <- NormalizeData(old_new_seurat, assay = "RNA")

###add in a new set of metadata###
donor_region = paste(old_new_seurat$autopsy_id, old_new_seurat$region, sep = "_")
old_new_seurat = AddMetaData(object = old_new_seurat, metadata = donor_region, col.name = "donor_region")
disease_region = paste(old_new_seurat$diagnosis, old_new_seurat$region, sep = "_")
old_new_seurat = AddMetaData(object = old_new_seurat, metadata = disease_region, col.name = "disease_region")
ALS_status = old_new_seurat$diagnosis
ALS_status[grep("ALS", ALS_status, invert = T)] = "non-ALS"
old_new_seurat = AddMetaData(object = old_new_seurat, metadata = ALS_status, col.name = "ALS_status")

#set a colorscheme for UMAP visualization#
colorscheme_hex = c("#49589F","#91CADB","#7A589C","#913B35","#F2EA57","#E5A745","#3D643D","#86BA5A","#CE4236","#6183B1","#EBC2CB","#A75E9C", "grey","black", "turquoise","coral2")
names(colorscheme_hex) = c(1:16)

##downsample v3 data##
v3_only = subset(old_new_seurat, cells = names(old_new_seurat$tech)[old_new_seurat$tech == "v3"])
v2_only = subset(old_new_seurat, cells = names(old_new_seurat$tech)[old_new_seurat$tech == "v2"])
count_mat = v3_only[["RNA"]]@counts
count_mat = downsampleMatrix(count_mat, prop = 0.5)
v3_only = CreateSeuratObject(count_mat, meta.data = v3_only@meta.data)
v2_only = CreateSeuratObject(v2_only[["RNA"]]@counts, meta.data = v2_only@meta.data)
old_new_seurat = merge(v3_only, v2_only)
old_new_seurat = NormalizeData(old_new_seurat)
old_new_seurat$filtered_cluster = factor(old_new_seurat$filtered_cluster, levels = c(1:16))
old_new_seurat = SetIdent(old_new_seurat, value = old_new_seurat$filtered_cluster)

####All cluster marker genes####
pairwise_markers=list()
for (ii in 1:16) {
  for (jj in 1:16) {
    nam=paste0(ii,"-",jj)
    marklist=FindMarkers(old_new_seurat,ident.1=ii,ident.2=jj,group.by = "filtered_cluster",metod="MAST")
    pairwise_markers[[nam]]=marklist
    print(c(ii,jj))
  }
}
save(pairwise_markers,file="pairwise_cluster_markers.rda")

####compile pairwise marker genes and ALS vs nonALS genes###
load("pairwise_cluster_markers.rda")
names(pairwise_markers)
###for each cluster, find genes that are upregulated versus at least one other cluster, and not in others###
load("old_new_seurat_mapped_clusters.rda")
###microglia only, no 15/16###
fullmat=c()
for (ii in c(1:9)) {
  allmat=data.frame(gene=rownames(old_new_seurat),uptype=ii,downtypes='n',num_downtypes=0,sum_log2foldchange=0,row.names = rownames(old_new_seurat),stringsAsFactors = F)
  keepnam=grep(paste0("^",ii,"-"),names(pairwise_markers),val=T)
  for (jj in keepnam) {
    dtype=strsplit(jj,"-")[[1]][2]
    if (dtype!=ii & dtype %in% c(1:9)) {
      mat1=pairwise_markers[[jj]]
      mat1=mat1[mat1$p_val_adj<0.05,]
      ###remove genes that are downregulated###
      remgenes=rownames(mat1)[mat1$avg_log2FC<0]
      allmat=allmat[setdiff(rownames(allmat),remgenes),]
      mat1=mat1[which(mat1$avg_log2FC>0 & rownames(mat1) %in% rownames(allmat)),]
      keeprows=match(rownames(mat1),rownames(allmat))
      allmat$downtypes[keeprows]=paste0(allmat$downtypes[keeprows],";",dtype)
      allmat$num_downtypes[keeprows]=allmat$num_downtypes[keeprows]+1
      allmat$sum_log2foldchange[keeprows]=allmat$sum_log2foldchange[keeprows]+mat1$avg_log2FC
    }
  }
  allmat2=allmat[allmat$num_downtypes>0,]
  allmat2$downtypes=gsub("n;","",allmat2$downtypes)
  allmat2=allmat2[order(-allmat2$num_downtypes),]
  fracvals=apply(old_new_seurat@assays$RNA[rownames(allmat2),which(old_new_seurat$filtered_cluster==ii)]>0,1,mean)
  allmat2$on_fraction=fracvals
  fullmat=rbind(fullmat,allmat2)
}
write.csv(fullmat,file="pairwise_upregulated_genes_microgliaclusters1-9_20210802.csv")

###same as above, but downregulated###
fullmat=c()
for (ii in c(1:9)) {
  allmat=data.frame(gene=rownames(old_new_seurat),downtype=ii,uptypes='n',num_uptypes=0,sum_log2foldchange=0,row.names = rownames(old_new_seurat),stringsAsFactors = F)
  keepnam=grep(paste0("^",ii,"-"),names(pairwise_markers),val=T)
  for (jj in keepnam) {
    dtype=strsplit(jj,"-")[[1]][2]
    if (dtype!=ii & dtype %in% c(1:9)) {
      mat1=pairwise_markers[[jj]]
      mat1=mat1[mat1$p_val_adj<0.05,]
      ###remove genes that are upregulated###
      remgenes=rownames(mat1)[mat1$avg_log2FC>0]
      allmat=allmat[setdiff(rownames(allmat),remgenes),]
      mat1=mat1[which(mat1$avg_log2FC<0 & rownames(mat1) %in% rownames(allmat)),]
      keeprows=match(rownames(mat1),rownames(allmat))
      allmat$uptypes[keeprows]=paste0(allmat$uptypes[keeprows],";",dtype)
      allmat$num_uptypes[keeprows]=allmat$num_uptypes[keeprows]+1
      allmat$sum_log2foldchange[keeprows]=allmat$sum_log2foldchange[keeprows]+mat1$avg_log2FC
    }
  }
  allmat2=allmat[allmat$num_uptypes>0,]
  allmat2$uptypes=gsub("n;","",allmat2$uptypes)
  allmat2=allmat2[order(-allmat2$num_uptypes),]
  fracvals=apply(old_new_seurat@assays$RNA[rownames(allmat2),which(old_new_seurat$filtered_cluster==ii)]>0,1,mean)
  allmat2$on_fraction=fracvals
  fullmat=rbind(fullmat,allmat2)
}
write.csv(fullmat,file="pairwise_downregulated_genes_microgliaclusters1-9_20210802.csv")

#compute ALS DEGs on dsamp data#
allde=list()
for (ii in c(1:16)) {
  de1=FindMarkers(old_new_seurat,ident.1=colnames(old_new_seurat)[which(old_new_seurat$filtered_cluster==ii & old_new_seurat$diagnosis=="ALS")],
                  ident.2=colnames(old_new_seurat)[which(old_new_seurat$filtered_cluster==ii & old_new_seurat$diagnosis!="ALS")],test.use="MAST")
  allde[[paste0("c",ii)]]=de1
  print(ii)
}
save(allde,file="dsamp_within_cluster_als_vs_nonals.rda")

###als vs non-als cluster data###
load("dsamp_within_cluster_als_vs_nonals.rda")
for (ii in 1:15) {
  mat1=allde[[paste0("c",ii)]]
  write.csv(mat1,paste0("dsamp_als_vs_nonals_",ii,".csv"))
}


#pathway annotation#
#load in the pairwise diffex lists#
uglia_out_stub = "degene_lists/"
#####visualize DE genes#######################
markers = read.csv(paste(uglia_out_stub, "pairwise_upregulated_genes_microgliaclusters1-9_20210802.csv", sep = ""))
markers = subset(markers, uptype %in% c(1:9))
markers_down = read.csv(paste(uglia_out_stub, "pairwise_downregulated_genes_microgliaclusters1-9_20210802.csv", sep = ""))
markers_down = subset(markers_down, downtype %in% c(1:9))
genes = rownames(old_new_seurat)

#now do reactome analysis#
entrez_names = mapIds(org.Hs.eg.db, genes, "ENTREZID", "SYMBOL")
FDR_thresh = 0.05
total_num = 50

#run a grouped clusterprofiler run#
genes_grouped_up = list()
genes_grouped_down = list()
for(i in unique(markers_down$downtype)){
  #up#
  module_df = data.frame(Entrez = entrez_names, group = 20)
  module_df$group[genes %in% markers$gene[markers$uptype == unique(markers$uptype)[i]][1:total_num]] = i
  module_df = module_df[!(is.na(module_df$Entrez)),]
  genes_grouped_up[[i]] = module_df$Entrez[module_df$group == i]
  #down#
  module_df = data.frame(Entrez = entrez_names, group = 20)
  module_df$group[genes %in% markers_down$gene[markers_down$downtype == unique(markers_down$downtype)[i]][1:total_num]] = i
  module_df = module_df[!(is.na(module_df$Entrez)),]
  genes_grouped_down[[i]] = module_df$Entrez[module_df$group == i]
}
names(genes_grouped_up) = c(1:9); names(genes_grouped_down) = c(1:9)

#now go ahead and do reactome analysis on what is up/down in ALS per cluster#
file_loc = "degene_lists/"

ALS_genes_grouped_up = list()
ALS_genes_grouped_down = list()
#now do it in the aggregate#
for(i in 1:9){
  markers = read.csv(paste0(file_loc, "dsamp_als_vs_nonals_", i, ".csv"))
  markers = markers[order(markers$avg_log2FC),]
  downreg_markers = markers$X[1:total_num]
  upreg_markers = markers$X[(dim(markers)[1]-total_num):dim(markers)[1]]
  
  #up#
  module_df = data.frame(Entrez = entrez_names, group = 20)
  module_df$group[genes %in% upreg_markers] = i
  module_df = module_df[!(is.na(module_df$Entrez)),]
  ALS_genes_grouped_up[[i]] = module_df$Entrez[module_df$group == i]
  #down#
  module_df = data.frame(Entrez = entrez_names, group = 20)
  module_df$group[genes %in% markers_down$gene[markers_down$downtype == unique(markers_down$downtype)[i]][1:total_num]] = i
  module_df$group[genes %in% downreg_markers] = i
  ALS_genes_grouped_down[[i]] = module_df$Entrez[module_df$group == i]
}
names(ALS_genes_grouped_up) = c(1:9); names(ALS_genes_grouped_down) = c(1:9)


#monocle3 analysis#
###run monocle pipeline with seurat clusters###
##convert to monocle object##
cds <- as.cell_data_set(old_new_seurat)
##monocle INSISTS that you do the clustering in monocle, although we'll try to hackily use the seurat IDs anyways##
cds <- cluster_cells(cds, resolution = 1e-5, num_iter = 5)
plot_cells(cds, show_trajectory_graph = FALSE)
cds@clusters$UMAP$clusters = Idents(old_new_seurat)
cds@clusters$UMAP$cluster_result$optim_res$membership = Idents(old_new_seurat)

##learn the pseudotime graph and overlay it onto our structure##
cds <- learn_graph(cds)
#choose the root#
max.avp <- which.max(unlist(FetchData(old_new_seurat, "AVP")))
max.avp <- colnames(old_new_seurat)[max.avp]
cds <- order_cells(cds, root_cells = max.avp)
cds <- order_cells(cds) #manually choose


###region/disease enrichment with Kruskal-Wallis###
DimPlot(old_new_seurat, group.by = "ALS_status")
##now, let's generate some nicer figures from the clusters differentially represented in ALS disease regions##
split.by = "donor_region"
conditions = "ALS"
##table for donor region versus the seurat identities##
donor_cluster_region_enrichment = prop.table(table(Idents(old_new_seurat), old_new_seurat@meta.data[,split.by]), margin = 2)
classifications_donor_region = data.frame(unique(old_new_seurat$donor_region)); colnames(classifications_donor_region) = "donor_region"
samp_locs = match(classifications_donor_region$donor_region, old_new_seurat$donor_region)
classifications_donor_region$region = old_new_seurat$region[samp_locs]
classifications_donor_region$diagnosis = old_new_seurat$diagnosis[samp_locs]
classifications_donor_region$gender = old_new_seurat$gender[samp_locs]
classifications_donor_region$ALS = old_new_seurat$ALS_status[samp_locs]
#in order#
classifications_donor_region = classifications_donor_region[order(classifications_donor_region$donor_region),]

#Evaluating ALS status vs. not#
for(j in conditions){
  kw.df = data.frame(matrix(ncol = 5, nrow = 0))
  for(i in 1:dim(donor_cluster_region_enrichment)[[1]]){
    #run the KW test#
    kw = kruskal.test(donor_cluster_region_enrichment[i,], classifications_donor_region[,j])
    kw_row = data.frame(i, j, kw[1], kw[2], kw[3])
    kw.df = rbind(kw.df, kw_row)
    temp = data.frame(cbind(donor_cluster_region_enrichment[i,], classifications_donor_region[,j]))
    colnames(temp) = c(paste("Cluster_proportion_", i, sep = ""), j)
    
    #generate ggplot figure
    ggplot(temp, aes(x=temp[[2]], y = as.numeric(temp[[1]]), fill = temp[[2]])) + geom_boxplotMod(alpha=0.9, color = "black") +  
      ylab("Proportion") + scale_fill_npg() + theme_classic() +  ggtitle(paste("Cluster", i)) +
      theme(axis.text.x = element_text(size = 40, color = "black"),axis.text.y = element_text(size = 30), panel.border = element_rect(colour = "black", fill=NA, size=1.25),
            axis.title.x = element_blank(), axis.title.y = element_text(size = 40), legend.position = "none", title = element_text(hjust = 0.5, size = 45)) 
    ggsave(paste("KW",j, i, ".png", sep = "_"), units = "mm", width = 325, height = 350, dpi = 300) 
  }
  kw.df = setNames(kw.df, nm =c("Cluster", "metadata", "chi-squared", "df", "p_val"))
  kw.df$p_val_BH = p.adjust(kw.df$p_val, "BH")
  kw.df$p_val_holm = p.adjust(kw.df$p_val, "holm")
  write.csv(kw.df, file = paste("KW", j, ".csv", sep ="_"))
}

##Evaluating region within ALS samples##
conditions = "region"

##looking exclusively in ALS##
classifications_donor_region = classifications_donor_region[classifications_donor_region$ALS == "ALS",]
donor_cluster_region_enrichment = donor_cluster_region_enrichment[,colnames(donor_cluster_region_enrichment) %in% classifications_donor_region$donor_region]

#statistics and figures for your condition of choice#
for(j in conditions){
  kw.df = data.frame(matrix(ncol = 5, nrow = 0))
  for(i in 1:dim(donor_cluster_region_enrichment)[[1]]){
    #run the KW test#
    kw = kruskal.test(donor_cluster_region_enrichment[i,], classifications_donor_region[,j])
    kw_row = data.frame(i, j, kw[1], kw[2], kw[3])
    kw.df = rbind(kw.df, kw_row)
    temp = data.frame(cbind(donor_cluster_region_enrichment[i,], classifications_donor_region[,j]))
    colnames(temp) = c(paste("Cluster_proportion_", i, sep = ""), j)
    
    #generate ggplot figure
    ggplot(temp, aes(x=temp[[2]], y = as.numeric(temp[[1]]), fill = temp[[2]])) + geom_boxplotMod(alpha=0.9, color = "black") +  
      ylab("Proportion") + scale_fill_npg() + theme_classic() +  ggtitle(paste("Cluster", i)) +
      theme(axis.text.x = element_text(size = 40, color = "black"),axis.text.y = element_text(size = 30), panel.border = element_rect(colour = "black", fill=NA, size=1.25),
            axis.title.x = element_blank(), axis.title.y = element_text(size = 40), legend.position = "none", title = element_text(hjust = 0.5, size = 45)) 
    ggsave(paste("KW_ALS_only",j, i, ".png", sep = "_"), units = "mm", width = 325, height = 350, dpi = 300) 
  }
  kw.df = setNames(kw.df, nm =c("Cluster", "metadata", "chi-squared", "df", "p_val"))
  kw.df$p_val_BH = p.adjust(kw.df$p_val, "BH")
  kw.df$p_val_holm = p.adjust(kw.df$p_val, "holm")
  write.csv(kw.df, file = paste("KW_ALS_only", j, ".csv", sep ="_"))
}