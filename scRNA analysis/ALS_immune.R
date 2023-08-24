library(Matrix)
library(Seurat)
library(dplyr)
library(readxl)
library(SeuratWrappers)
library(batchelor)
library(pheatmap)
library(patchwork)
library(ggplot2)
library(DESeq2)
library(ggfortify)
library(ggrepel)
library(corrplot)
library(DESeq2)
library(MAST)
library(edgeR)
library(ReactomePA)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(DropletUtils)
library(schex)
library(viridis)
library(ggsci)
library(igraph)
library(harmony)
library(data.table)

####setup locations and files to read in. Please ensure that you are using R 4.1.0+ for downstream analysis ####
source("clustering_diffex_functions.R")
source("ChooseR/helper_functions.R")
source("ChooseR/pipeline.R")
#functions#
chooseR <- function(obj, npcs = 20, resolutions = c(0.3,0.6,0.9,1.2,1.5,2,4,6), assay = "SCT", reduction = "pca", results_path = "results/"){
  # Run pipeline
  for (res in resolutions) {
    message(paste0("Clustering ", res, "..."))
    message("\tFinding ground truth...")
    
    # "Truths" will be stored at glue::glue("{reduction}.{assay}_res.{res}")
    obj <- find_clusters(
      obj,
      reduction = reduction,
      assay = assay,
      resolution = res,
      npcs = npcs
    )
    clusters <- obj[[glue::glue("{reduction}.{assay}_res.{res}")]]
    
    # Now perform iterative, sub-sampled clusters
    results <- multiple_cluster(
      obj,
      n = 50,
      size = 0.8,
      npcs = npcs,
      res = res,
      reduction = reduction,
      assay = assay
    )
    
    # Now calculate the co-clustering frequencies
    message(paste0("Tallying ", res, "..."))
    # This is the more time efficient vectorisation
    # However, it exhausts vector memory for (nearly) all datasets
    # matches <- purrr::map(columns, find_matches, df = results)
    # matches <- purrr::reduce(matches, `+`)
    columns <- colnames(dplyr::select(results, -cell))
    mtchs <- matrix(0, nrow = dim(results)[1], ncol = dim(results)[1])
    i <- 1 # Counter
    for (col in columns) {
      message(paste0("\tRound ", i, "..."))
      mtchs <- Reduce("+", list(
        mtchs,
        find_matches(col, df = results)
      ))
      i <- i + 1
    }
    
    message(paste0("Scoring ", res, "..."))
    mtchs <- dplyr::mutate_all(
      dplyr::as_tibble(mtchs),
      function(x) dplyr::if_else(Re(x) > 0, percent_match(x), 0)
    )
    
    # Now calculate silhouette scores
    message(paste0("Silhouette ", res, "..."))
    sil <- cluster::silhouette(
      x = as.numeric(as.character(unlist(clusters))),
      dmatrix = (1 - as.matrix(mtchs))
    )
    saveRDS(sil, paste0(results_path, "silhouette_", res, ".rds"))
    
    # Finally, calculate grouped metrics
    message(paste0("Grouping ", res, "..."))
    grp <- group_scores(mtchs, unlist(clusters))
    saveRDS(grp, paste0(results_path, "frequency_grouped_", res, ".rds"))
    sil <- group_sil(sil, res)
    saveRDS(sil, paste0(results_path, "silhouette_grouped_", res, ".rds"))
  }
  return(obj)
}
chooseR_viz <- function(obj, npcs = 20, resolutions = c(0.3,0.6,0.9,1.2,1.5,2,4,6), assay = "SCT", reduction = "pca", results_path = "results/"){
  scores <- purrr::map(
    paste0(results_path, "silhouette_grouped_", resolutions, ".rds"),
    readRDS
  )
  scores <- dplyr::bind_rows(scores) %>%
    dplyr::group_by(res) %>%
    dplyr::mutate("n_clusters" = dplyr::n()) %>%
    dplyr::ungroup()
  meds <- scores %>%
    dplyr::group_by(res) %>%
    dplyr::summarise(
      "boot" = list(boot_median(avg_sil)),
      "n_clusters" = mean(n_clusters)
    ) %>%
    tidyr::unnest_wider(boot)
  
  writexl::write_xlsx(meds, paste0(results_path, "median_ci.xlsx"))
  
  # Find thresholds
  threshold <- max(meds$low_med)
  choice <- as.character(
    meds %>%
      dplyr::filter(med >= threshold) %>%
      dplyr::arrange(n_clusters) %>%
      tail(n = 1) %>%
      dplyr::pull(res)
  )
  
  # And plot!
  ggplot(meds, aes(factor(res), med)) +
    geom_crossbar(
      aes(ymin = low_med, ymax = high_med),
      fill = "grey",
      size = 0.25
    ) +
    geom_hline(aes(yintercept = threshold), colour = "blue") +
    geom_vline(aes(xintercept = choice), colour = "red") +
    geom_jitter(
      data = scores,
      aes(factor(res), avg_sil),
      size = 0.35,
      width = 0.15
    ) +
    scale_x_discrete("Resolution") +
    scale_y_continuous(
      "Silhouette Score",
      expand = c(0, 0),
      limits = c(-1, 1),
      breaks = seq(-1, 1, 0.25),
      oob = scales::squish
    ) +
    cowplot::theme_minimal_hgrid() +
    theme(
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 7),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black"),
      axis.ticks = element_line(colour = "black"),
    )
  
  ggsave(
    filename = paste0(results_path, "silhouette_distribution_plot.png"),
    dpi = 300,
    height = 3.5,
    width = 3.5,
    units = "in"
  )
  
  # Finally, a dot plot of silhouette scores to help identify less robust clusters
  # The initial pipe is to order the clusters by silhouette score
  scores %>%
    dplyr::filter(res == choice) %>%
    dplyr::arrange(dplyr::desc(avg_sil)) %>%
    dplyr::mutate_at("cluster", ordered, levels = .$cluster) %>%
    ggplot(aes(factor(cluster), avg_sil)) +
    geom_point() +
    scale_x_discrete("Cluster") +
    scale_y_continuous(
      "Silhouette Score",
      expand = c(0, 0),
      limits = c(-1, 1),
      breaks = seq(-1, 1, 0.25),
      oob = scales::squish
    ) +
    cowplot::theme_minimal_grid() +
    theme(
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 7),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black"),
      axis.ticks = element_line(colour = "black"),
    )
  
  ggsave(
    filename = paste0(results_path, "silhouette_point_plot_", choice, ".png"),
    dpi = 300,
    height = 3.5,
    width = 3.5,
    units = "in"
  )
  return(choice)
}
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

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

##Immune cell subclustering from the data##
#first check to see which are immune and which are other#
immune_seurat = subset(old_new_seurat, idents = c(10:12))
other_seurat = subset(old_new_seurat, idents = c(13))
rm(old_new_seurat); gc()

#get and subcluster the immune cells#
immune_seurat=SCTransform(immune_seurat)
immune_seurat=RunPCA(immune_seurat)
thetaval=2
dimval=20
immune_seurat=RunHarmony(immune_seurat,assay.use = "SCT",theta=thetaval,group.by.vars="orig.ident")
immune_seurat=RunUMAP(immune_seurat,reduction="harmony",dims=1:dimval)
immune_seurat=FindNeighbors(immune_seurat,reduction="harmony",dims=1:dimval)
immune_seurat=FindClusters(immune_seurat,resolution=c(0.2, 0.3, 0.5))

#run ChooseR to identify the correct clustering parameter#
immune_seurat_chooseR <- chooseR(obj = immune_seurat, npcs = 20, resolutions = c(0.25, 0.35, 0.5, 0.6, 0.75, 0.9,1,1.25,1.5,2,3,4,6), assay = "SCT", reduction = "harmony", results_path = "2021_alt/Analysis_Outs/ChooseR/")
choice <- chooseR_viz(obj = immune_seurat_chooseR, npcs = 20, resolutions = c(0.25, 0.35, 0.5, 0.6, 0.75, 0.9,1,1.25,1.5,2,3,4,6), assay = "SCT", reduction = "harmony", results_path = "2021_alt/Analysis_Outs/ChooseR/")
immune_seurat_chooseR <- SetIdent(object = immune_seurat_chooseR, value = paste0("harmony.SCT_res.", choice))
DefaultAssay(immune_seurat_chooseR) = "RNA"
immune_seurat_chooseR = NormalizeData(immune_seurat_chooseR)
save(immune_seurat_chooseR,file="immune_cell_seurat.rda")

#restart after loading
load("immune_cell_seurat.rda")
new.ids = 1:length(unique(immune_seurat_chooseR[[paste0("harmony.SCT_res.", choice)]])[,1])
names(new.ids) <- 0:(length(unique(immune_seurat_chooseR[[paste0("harmony.SCT_res.", choice)]])[,1])-1)
immune_seurat_chooseR <- RenameIdents(immune_seurat_chooseR, new.ids)
immune_seurat_chooseR <- subset(immune_seurat_chooseR, idents = c(1:19)) #removing the <10cell clusters

#add the final identities back in#
immune_seurat_chooseR$final_ident = Idents(immune_seurat_chooseR)

#do a KW enrichment for initial evaluation#
KW_enrichment(immune_seurat_chooseR, split.by = "donor_region", conditions = c("gender", "ALS"), filestub = "2021_alt/Analysis_Outs/KW/immune/2022_KW_", label = F)

###immune cluster marker genes###
pairwise_markers=list()
for (ii in 1:19) {
  for (jj in 1:18) {
    nam=paste0(ii,"-",jj)
    marklist=FindMarkers(immune_seurat_chooseR,ident.1=ii,ident.2=jj,group.by = "final_ident",method="MAST")
    pairwise_markers[[nam]]=marklist
    print(c(ii,jj))
  }
}
save(pairwise_markers,file="immune_pairwise.rda")


####compile pairwise marker genes and ALS vs nonALS genes###
load("immune_pairwise.rda")
names(pairwise_markers)
###for each cluster, find genes that are upregulated versus at least one other cluster, and not in others###
###all clusters###
fullmat=c()
for (ii in 1:19) {
  allmat=data.frame(gene=rownames(immune_seurat_chooseR),uptype=ii,downtypes='n',num_downtypes=0,sum_log2foldchange=0,row.names = rownames(immune_seurat_chooseR),stringsAsFactors = F)
  keepnam=grep(paste0("^",ii,"-"),names(pairwise_markers),val=T)
  for (jj in keepnam) {
    dtype=strsplit(jj,"-")[[1]][2]
    if (dtype!=ii) {
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
  fracvals=apply(immune_seurat_chooseR@assays$RNA[rownames(allmat2),which(immune_seurat_chooseR$final_ident==ii)]>0,1,mean)
  allmat2$on_fraction=fracvals
  fullmat=rbind(fullmat,allmat2)
}
write.csv(fullmat,file="/pairwise_upregulated_genes_immuneclusters.csv")

###same as above, but downregulated###
load("immune_pairwise.rda")
fullmat=c()
for (ii in 1:19) {
  allmat=data.frame(gene=rownames(immune_seurat_chooseR),downtype=ii,uptypes='n',num_uptypes=0,sum_log2foldchange=0,row.names = rownames(immune_seurat_chooseR),stringsAsFactors = F)
  keepnam=grep(paste0("^",ii,"-"),names(pairwise_markers),val=T)
  for (jj in keepnam) {
    dtype=strsplit(jj,"-")[[1]][2]
    if (dtype!=ii) {
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
  fracvals=apply(immune_seurat_chooseR@assays$RNA[rownames(allmat2),which(immune_seurat_chooseR$final_ident==ii)]>0,1,mean)
  allmat2$on_fraction=fracvals
  fullmat=rbind(fullmat,allmat2)
}
write.csv(fullmat,file="pairwise_downregulated_genes_immuneclusters.csv")

#now annotate markers in the aggregate#
genes = rownames(immune_seurat_chooseR)
entrez_names = mapIds(org.Hs.eg.db, genes, "ENTREZID", "SYMBOL")
FDR_thresh = 0.05
total_num = 50
for(i in 1:19){
  markers = read.csv("2021_alt/Analysis_Outs/immune_DEGs/pairwise_upregulated_genes_immuneclusters_07252022.csv")
  markers <- markers[markers$uptype == i,]
  upreg_markers = markers$gene[1:total_num]
  
  #up#
  module_df = data.frame(Entrez = entrez_names, group = 20)
  module_df$group[genes %in% upreg_markers] = i
  module_df = module_df[!(is.na(module_df$Entrez)),]
  ALS_genes_grouped_up[[i]] = module_df$Entrez[module_df$group == i]
}
names(ALS_genes_grouped_up) = c(1:19)


#azimuth mapping#
#restart after loading
load("immune_cell_seurat.rda")
new.ids = 1:length(unique(immune_seurat_chooseR[[paste0("harmony.SCT_res.", choice)]])[,1])
names(new.ids) <- 0:(length(unique(immune_seurat_chooseR[[paste0("harmony.SCT_res.", choice)]])[,1])-1)
immune_seurat_chooseR <- RenameIdents(immune_seurat_chooseR, new.ids)
#projected UMAP
projected.umap <- readRDS('azimuth/azimuth_umap.Rds')
immune_seurat_chooseR <- immune_seurat_chooseR[, Cells(projected.umap)]
immune_seurat_chooseR[['umap.proj']] <- projected.umap
#imputed ADT assay#
imputed.assay <- readRDS('azimuth/azimuth_impADT.Rds')
immune_seurat_chooseR <- immune_seurat_chooseR[, Cells(imputed.assay)]
immune_seurat_chooseR[['impADT']] <- imputed.assay
#cell predictions#
predictions <- read.delim('azimuth/azimuth_pred.tsv', row.names = 1)
immune_seurat_chooseR <- AddMetaData(
  object = immune_seurat_chooseR,
  metadata = predictions)
#set to primary predicted cell types#
immune_seurat_chooseR <- SetIdent(immune_seurat_chooseR, value = "predicted.celltype.l2")

#KW enrichment for azimuth predictions
KW_enrichment(immune_seurat_chooseR, split.by = "donor_region", conditions = c("gender", "ALS"), filestub = "2021_alt/Analysis_Outs/KW/immune/azimuth/2022_KW_", label = F)
