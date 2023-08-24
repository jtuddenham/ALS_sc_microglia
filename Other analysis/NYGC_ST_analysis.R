library(Matrix)
library(readxl)
library(readxl)
library(patchwork)
library(ggplot2)
library(ggsci)
library(ggfortify)
library(ggrepel)
library(tidyverse)
library(reshape2)

stmeta = read_xlsx("Metadata_05132019_SM.xlsx") ###remove/change this?###
# stmeta=read.table("metadata/human_sample_names_sra.tsv",as.is=T,skip=1)
allhuman=list()
for (ii in 1:nrow(stmeta)) {
  temp=read.table(paste0("count_matrices/",stmeta$V1[ii],"_stdata_aligned_counts_IDs.txt"),as.is=T,check.names=F,header=T,row.names=NULL)
  temp=temp[grep("ambiguous",temp[,1],invert=T),]
  temp=temp[which(!duplicated(temp[,1])),]
  rownames(temp)=temp[,1]
  temp=as.matrix(temp[,-1])
  allhuman[[stmeta$V1[ii]]]=Matrix(temp)
}
save(allhuman,file="allhuman_data.rda")

#first, reconstruct the data list from the new posterior mean data#
posterior_mean_human_ST = list()
file_names = list.files("posterior_means/lambda_means_human/")
for(name in file_names){
  data = read_csv(paste0("posterior_means/lambda_means_human/", name))
  sparse_mat_dat = Matrix(data = as.matrix(data[,-1]), sparse = T)
  row.names(sparse_mat_dat) = data[,1][[1]]
  posterior_mean_human_ST[[strsplit(name, "[.]")[[1]][1]]] = sparse_mat_dat
}
names(posterior_mean_human_ST) = sapply(names(posterior_mean_human_ST), FUN = function(x){return(substr(x, start = 1, stop = str_length(x)-8))})
save(posterior_mean_human_ST, file = "allhuman_posteriors.rda")

#now use the predicted spot depth#
#knowing that the median spot depth is 1203, calculate spot depth for each spot across all genes#
load("allhuman_data.rda")
allspot_depth = c()
for(samp in allhuman){
  allspot_depth = c(allspot_depth, colSums(samp))
}

#now, using this and the posterior means to estimate the predicted counts: pred_counts = lambda * spot_depth/median_spot_depth#
load("allhuman_posteriors.rda")
allhuman_pred_counts = list()
for(i in 1:length(posterior_mean_human_ST)){
  samp = posterior_mean_human_ST[[i]]
  pred_count_dat = sapply(X = colnames(samp), function(x){
    return((samp[,x] * allspot_depth[x])/1203)
  })
  allhuman_pred_counts[[names(posterior_mean_human_ST)[[i]]]] = pred_count_dat
}
save(allhuman_pred_counts, file = "allhuman_predicted_counts.rda")

#now, conduct a paired t-test between DH and VH acoss all genes#
#function requires two "region" matrices with identical rownames and column names#
#no log-transform upstream- all norm done here#
paired_ttest_diffex_linlog <- function(region1, region2, factor = 10^6){
  if( any(dim(region1) != dim(region2)) | any(rownames(region1) != rownames(region2)) | any(colnames(region1) != colnames(region2)) ){
    stop("Dissimilar matrices. Check dimensions, rownames, and colnames.")
  }
  final_vals = as.data.frame(matrix(nrow = dim(region1)[2], ncol = 7))
  for(gene in 1:dim(region1)[2]){
    t1 = sum(region1[,gene]); t2 =  sum(region2[,gene])
    if (t1 + t2 > 0){
      ttest = t.test(x= log2(region1[,gene]*factor+1), y= log2(region2[,gene]*factor+1), paired = T, var.equal = F, conf.level = 0.95)
      tmp = data.frame(t(unlist(ttest)[1:5]))
      if (t1==0) {
        lfc=Inf
      } else if(t2 == 0){
        lfc=-Inf
      } else{
        lfc = log2(mean(region2[,gene])/mean(region1[,gene]))
      }
      tmp = cbind(tmp, lfc, T)
      final_vals[gene,] = tmp
      colnames(final_vals) = colnames(tmp)
    } else{
      final_vals[gene,7] =  F
    }
  }
  final_vals$p_val_BH = p.adjust(final_vals$p.value, "BH")
  final_vals$gene = colnames(region1)
  return(final_vals)
}


###now let's try to rank based on genes ####
load("allhuman_posteriors.rda")
load("allhuman_predicted_counts.rda")

meta1=read.table("metadata/human_sample_names_sra.tsv",as.is=T,skip = 1,sep="\t")

#subset the samples for the high-quality samples used in analysis#
retain_list = c("L8CN13_C1", "L8CN13_D2", "L8CN13_E1", "L8CN13_E2", "L8CN14_C2", "L8CN14_D1", "L8CN14_D2", "L8CN14_E1", "L8CN14_E2", "L8CN155_C1", "L8CN156_D2", "L8CN156_E1", "L8CN6_C1", "L8CN7_E2", 
                "L8CN153_C1", "L8CN154_C1", "L8CN154_C2", "L8CN154_D1", "L8CN154_E1", "L8CN151_D2", "L8CN152_C1", "L8CN152_C2", "L8CN159_D1", "L8CN159_E2", "L8CN160_C1", "L8CN160_C2")
posterior_mean_human_ST = posterior_mean_human_ST[names(posterior_mean_human_ST) %in% retain_list]
allhuman_pred_counts = allhuman_pred_counts[names(allhuman_pred_counts) %in% retain_list]

###get signatures for each cluster###
uglia_genes=read_xls("ALS cluster signatures_v2.xls")
uglia_genes = as.data.frame(as.table(as.matrix(uglia_genes))) ##make it long###
uglia_genes = uglia_genes[,-1]
colnames(uglia_genes) = c("cluster", "gene")
####read in annotations###
anno_files = list.files(".")[grep("annotations", list.files("."))]
sampnames = sapply(strsplit(anno_files, split = "[.]"), FUN = "[", 1)
human_names = names(posterior_mean_human_ST)
tempmeta = read.table(anno_files[1], header = T)
regions = colnames(tempmeta[,5:17])
#focus the analysis on dorsal versus ventral horns#
regions = regions[regions %in% c("Vent_Horn", "Dors_Horn")]

###find the aggregate unique gene list-should be the same for both###
gene_num = list()
max = 0
max_loc = 0
genes = c()
for(samp in 1:length(posterior_mean_human_ST)){
  if(dim(posterior_mean_human_ST[[samp]])[[1]] > max){
    max = dim(posterior_mean_human_ST[[samp]])[[1]]
    max_loc = samp
  }
  gene_num[[samp]] = dim(posterior_mean_human_ST[[samp]])[[1]]
  genes = union(genes, rownames(posterior_mean_human_ST[[samp]]))
}

#done with predicted counts#
#lump regions together so you can directly compare overall region as a pseudobulk lump versus other region#
regional_matrices = list()
sum_matrices = list()
median_matrices = list()
sd_matrices = list()
for(region in regions){
  regmat = Matrix(ncol = length(genes), nrow = 0)
  summat = Matrix(ncol = length(genes), nrow = 0)
  medianmat = Matrix(ncol = length(genes), nrow = 0)
  sdmat = Matrix(ncol = length(genes), nrow = 0)
  colnames(regmat) = genes
  for(samp in 1:length(allhuman_pred_counts)) {
    startmat=as.matrix(t(allhuman_pred_counts[[samp]]))
    ###dividing by sum of each row, but not log-normalizing###
    startmat=sweep(startmat,1,rowSums(startmat),"/")
    startnam = human_names[samp]
    startmeta = read.table(anno_files[grep(startnam, sampnames)], header = T)
    startmeta$combined_pos = paste(startmeta$xPos, startmeta$yPos, sep = "_")
    startmeta$location = colnames(startmeta[,5:17])[max.col(startmeta[,5:17],ties.method="first")]
    if(sum(rownames(startmat) %in% startmeta$combined_pos[startmeta$location == region]) == 1){
      reg_mat = t(Matrix(startmat[rownames(startmat) %in% startmeta$combined_pos[startmeta$location == region],]))
      rownames(reg_mat) = rownames(startmat)[rownames(startmat) %in% startmeta$combined_pos[startmeta$location == region]]
    } else{
      reg_mat = startmat[rownames(startmat) %in% startmeta$combined_pos[startmeta$location == region],]
    }
    ###sort out the genes that are missing in each sample so we can jam all these matrices together##
    missing_genes = colnames(regmat)[!(colnames(regmat) %in% colnames(reg_mat))]
    missing_mat = Matrix(0, nrow = dim(reg_mat)[[1]], ncol = length(missing_genes), sparse = F)
    colnames(missing_mat) = missing_genes
    reg_mat = cbind(reg_mat, missing_mat)
    reg_mat = reg_mat[,genes]
    sum_mat = apply(reg_mat, MARGIN = 2, sum)
    med_mat = apply(reg_mat, MARGIN = 2, median)
    sd_mat = apply(reg_mat, MARGIN = 2, sd)
    medianmat = rbind(medianmat, med_mat)
    sdmat = rbind(sdmat, sd_mat)
    summat = rbind(summat, sum_mat)
    regmat = rbind(regmat, reg_mat)
  }
  colnames(summat) = genes
  rownames(summat) = names(allhuman_pred_counts)
  colnames(medianmat) = genes; colnames(sdmat) = genes
  rownames(medianmat) = names(allhuman_pred_counts); rownames(sdmat) = names(allhuman_pred_counts)
  median_matrices[[region]] = medianmat
  sd_matrices[[region]] = sdmat
  regional_matrices[[region]] = regmat
  sum_matrices[[region]] = summat
}

#now retrieve t-test results for all genes#
DH_VH_ttest_results = paired_ttest_diffex_linlog(region1 = sum_matrices$Vent_Horn, region2 = sum_matrices$Dors_Horn, factor = 10^6)
write.csv(x = DH_VH_ttest_results, file = "DH_VH_pairedTtest.csv")