
Sys.setenv(RETICULATE_PYTHON = "/.conda/envs/leiden/bin/python")
Sys.setenv(RETICULATE_PYTHON = "/.conda/envs/leiden/bin/python")

library(reticulate)
library(Seurat)
library(scalop)
library(leiden)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(scales)
library(reshape2)
library(reshape2)
library(scales)
library(NMF)
library(MetBrewer)
library(colorspace)
library(tibble)
library(dplyr)
library(data.table)
library(stringr)
library(readr)
library(Matrix)
library(bigmemory)
library(doMC)
library(patchwork)

###### Per sample Leiden clustering##############################

samples_names <- (read.delim("general/GBM_samples.txt", header = FALSE, sep = "\t"))$V1

# set parameters 
complexity_filter <- 1000
mitochondrial_filter <- 20
genes_filter <- 7000
n_dim <- 20
dim_filter <- 5^(-15)
res_param <- 1
sig_th <- .005
mp_num <- 13
distinct16_pal<-c("#11A579","#F2B701","#66C5CC","#80BA5A","#F6CF71","#7F3C8D","#CF1C90","#3969AC","#f97b72","#E73F74","#4b4b8f","#ACA4E2","#31C53F","#B95FBB","#D4D915","#28CECA")
dim2use_list <- c(8,20,19,9,12,14,16,10,10,14,8,10,12,16,7,5,14,16,4,12,6,6,8,12,13,10) # set by jackstraw

leiden_clustering <- sapply(c(1:length(samples_names)), function(i){
  print(samples_names[i])
  # load spatial data
  unfilt_obj<-Load10X_Spatial(data.dir = paste("general/GBM_data/",samples_names[i],"/outs", sep = ""),
                              filename = "filtered_feature_bc_matrix.h5",
                              assay = "Spatial",
                              slice = "detected_tissue_image.jpg",
                              filter.matrix = TRUE,
                              to.upper = FALSE
  )
  
  # filtering
  exp_obj <- subset(unfilt_obj, subset = nCount_Spatial > complexity_filter) # filter out spots with less than #complexity_filter UMIs
  exp_obj[["percent.mt"]] <- PercentageFeatureSet(exp_obj, pattern = "^MT-") # filter spots with high percent (mitochondrial_filter) mitochondrial genes
  exp_obj <- subset(exp_obj, subset = percent.mt<mitochondrial_filter)
  
  # normalization and centering
  exp_obj <- NormalizeData(exp_obj, assay = "Spatial", normalization.method="LogNormalize", scale.factor=10000)
  exp_obj <- FindVariableFeatures(exp_obj, selection.method = "vst", nfeatures = genes_filter) # use only #genes_filter most variable genes
  exp_obj <- ScaleData(exp_obj)
  
  # PCA
  exp_obj <- RunPCA(exp_obj, features = VariableFeatures(object = exp_obj))
  #PCA dims to use determined by Jackstraw
  #exp_obj <- JackStraw(exp_obj, num.replicate = 100)
  #exp_obj <- ScoreJackStraw(exp_obj, dims = 1:n_dim)
  #js_scores <- exp_obj@reductions$pca@jackstraw$overall.p.values
  #dim2use <- max(js_scores[,"PC"][js_scores[,"Score"] < dim_filter])
  dim2use <- dim2use_list[i]
  
  # clustering
  exp_obj <- FindNeighbors(exp_obj, dims = 1:dim2use)
  exp_obj <- FindClusters(exp_obj, algorithm=4, resolution = res_param) # leiden based clustering

  spatial_sample_programs <- FindAllMarkers(exp_obj, only.pos = TRUE, return.thresh = sig_th)
  
  genes_num <- table(spatial_sample_programs$cluster)
  spots_clusters <- exp_obj@meta.data[,"seurat_clusters"]
  names(spots_clusters) <- row.names(exp_obj@meta.data)
  cluster_num <- table(exp_obj@meta.data$seurat_clusters)
  
  return(list(spatial_sample_programs,genes_num,spots_clusters,cluster_num)) # output: #genes per cluster, #spots per cluster, clusters genes sig, spots assignment to cluster
  
})


# clusters quality 
genes_num <- c()
spots_num <- c()
sapply(c(1:length(samples_names)), function(i){
  genes_num <<- c(genes_num,as.numeric(leiden_clustering[2,i][[1]]))
  spots_num <<- c(spots_num,as.numeric(leiden_clustering[4,i][[1]]))
})
hist(genes_num, breaks = 30)
hist(spots_num, breaks = 30)

# create genes sig table for save
genes_sig <- leiden_clustering[1,1][[1]]
genes_sig$sample <- rep(samples_names[1],dim(genes_sig)[1])

sapply(c(2:length(samples_names)),function(i){
  genes_sig_temp <- leiden_clustering[1,i][[1]]
  genes_sig_temp$sample <- rep(samples_names[i],dim(genes_sig_temp)[1])
  genes_sig <<- rbind(genes_sig,genes_sig_temp)
})

# create clusters assignment table for save
spots_assign <- as.data.frame(leiden_clustering[3,1][[1]])
spots_assign$sample <- rep(samples_names[1],dim(spots_assign)[1])
row.names(spots_assign) <- paste(row.names(spots_assign), "_", samples_names[1], sep = "")
colnames(spots_assign) <- c("cluster","sample")

sapply(c(2:length(samples_names)),function(i){
  spots_assign_temp <- as.data.frame(leiden_clustering[3,i][[1]])
  spots_assign_temp$sample <- rep(samples_names[i],dim(spots_assign_temp)[1])
  row.names(spots_assign_temp) <- paste(row.names(spots_assign_temp), "_", samples_names[i], sep = "")
  colnames(spots_assign_temp) <- c("cluster","sample")
  spots_assign <<- rbind(spots_assign,spots_assign_temp)
})

# set gene signatures list in the right format 
genes_list <- lapply(c(1:length(samples_names)), function(i){
  sample_table <- leiden_clustering[1,i][[1]]
  sample_genes <- lapply(c(1:length(unique(sample_table$cluster))),function(j){
    cluster_table <- sample_table[sample_table$cluster == j,]
    genes <- cluster_table$gene[order(cluster_table$avg_log2FC, decreasing = TRUE)][1:min(50,nrow(cluster_table))]
    return(genes)
  })
  return(sample_genes)
})
genes_list <- unlist(genes_list, recursive = FALSE)
names(genes_list) <- unlist(sapply(c(1:length(samples_names)),function(i){
  paste(samples_names[i],c(1:length(leiden_clustering[4,i][[1]])), sep = "_")
}))

samples_list <- unlist(sapply(c(1:length(samples_names)),function(i){
  rep(samples_names[i],length(leiden_clustering[4,i][[1]]))
}))
#saveRDS(genes_list,"/home/labs/tirosh/spatial_glioma/results/leiden_MPs/all_leiden_programs_v2.rds")
# generate an ordered jaccard matrix and clusters 
jac_mat = Jaccard(genes_list)
jac_matord = scalop::hca_reorder(jac_mat)
jac_plot <- melt(jac_matord)

prog_ordered <- colnames(jac_matord)
clusters <- scalop::hca_groups(jac_matord,
                               cor.method="none",
                               k=mp_num,
                               min.size=4,
                               max.size=0.5)

clusterCut <- sapply(names(clusters), function(i){
  new_vec <- rep(i,length(clusters[[i]]))
  names(new_vec) <- clusters[[i]]
  return(new_vec)
})

names(clusterCut) <- NULL
clusterCut <- unlist(clusterCut)

# plot jaccard 

jac_plot$value[jac_plot$value >= 0.35] <- 0.35
p1 <- ggplot(data = jac_plot, aes(x = Var1, y = Var2, fill = value)) + geom_raster() + 
  scale_fill_gradientn(colours = colorspace::sequential_hcl(9, palette="YlOrRd", rev=T), limits = c(0,0.35), name = "Jaccard\nIndex") +
  theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),  axis.line = element_blank(), legend.text.align = 0.5) +
  scale_x_discrete(name = "\nPrograms", labels = element_blank(), breaks = element_blank()) + scale_y_discrete(name = "\nPrograms", labels = element_blank(), breaks = element_blank())


metaprog_df <- data.frame(row.names = prog_ordered,
                          cluster = as.factor(clusterCut[prog_ordered]),
                          samples = stringr::str_replace(prog_ordered, "_.*",""),
                          cells = prog_ordered)

samples_annotation_plot <- ggplot(metaprog_df, aes(x = factor(cells, levels = cells), y = 1)) +
  geom_raster(aes(fill = samples)) + theme(legend.position = "none",axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), axis.text = element_blank(), axis.line = element_blank(), axis.text.x = element_blank()) +
  labs(x = "", y = "samples", fill = "")
metaprog_annotation_plot <- ggplot(metaprog_df, aes(x = factor(cells, levels = cells), y = 1)) +
  geom_raster(aes(fill = cluster)) + theme(legend.position = "none", axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), axis.text = element_blank(), axis.line = element_blank(), axis.text.x = element_blank()) +
  labs(x = "", y = "", fill = "")
MP_plot <- egg::ggarrange(p1, metaprog_annotation_plot, samples_annotation_plot, ncol = 1, nrow = 3, heights = c(40, 2, 1))


#Define meta-programs from program clusters, Sort by frequency across programs in cluster (>= 3)
cut_off <- table(metaprog_df$cluster)*0.25
cut_off[cut_off < 3] <- 3

mp_freqs = sapply(clusters, function(k) sort(table(unlist(genes_list[k])),decreasing = T),simplify = F)
metaprograms = sapply(names(mp_freqs), function(c) head(names(mp_freqs[[c]])[mp_freqs[[c]] >= as.numeric(cut_off[c])], 50), simplify = F)
names(metaprograms)<-c("vascular","neuron","hypoxia","hypoxia.immune","oligo","AC.immune","low.quality","NPC2","macrophage","MES.stromal","OPC.AC","metabolism","MES1.RGL")


metaprog_df$samp_cluster <- str_split(metaprog_df$cells, "_", simplify = T)[,2]
metaprog_df$anno_cluster <- metaprog_df$cluster
class(metaprog_df$anno_cluster)
levels(metaprog_df$anno_cluster) <- c("vascular","MES.stromal","OPC.AC","metabolism","MES1.RGL","neuron","hypoxia","hypoxia.immune","oligo","AC.immune","low.quality","NPC2","macrophage")

sapply(c(1:26), function(i){
  print(samples_names[i])
  spots_gen_assign <- as.data.frame(leiden_clustering[3,i][[1]])
  spots_gen_assign$generalized <- apply(spots_gen_assign,1,function(r){
    return(metaprog_df$anno_cluster[metaprog_df$samples == samples_names[i] & metaprog_df$samp_cluster == r])
  })
  spots_gen_assign$barcodes <- row.names(spots_gen_assign)
  spots_gen_assign <- spots_gen_assign[,c("barcodes", "generalized")]
})





###### NMF (Previously run by server. Possible to run per sample below)---------------------------------------------------------------------

#!/usr/bin/env Rscript
#args = commandArgs(trailingOnly = TRUE)


sname <- "MGH258"# enter sample name here 
rank_lb <- 2
rank_ub <- 11

m <- readRDS(paste0("MP/NMF/NMF_mats_GBM/", sname, ".rds"))
m <- as.matrix(m)
res <- NMF::nmf(x = m, rank = rank_lb:rank_ub, nrun = 5, method = "snmf/r", .opt = list(debug=F, parallel=F, shared.memory=F, verbose=T))


###### Generation of NMF + Leiden GBM spatial metaprograms#########

#

#robust_nmf_programs function

robust_nmf_programs <- function(nmf_programs, intra_min = 30, intra_max = 10, inter_filter=T, inter_min = 10) {
  
  # Select NMF programs based on the minimum overlap with other NMF programs from the same sample
  intra_intersect <- lapply(nmf_programs, function(z) apply(z, 2, function(x) apply(z, 2, function(y) length(intersect(x,y))))) 
  intra_intersect_max <- lapply(intra_intersect, function(x) apply(x, 2, function(y) sort(y, decreasing = T)[2]))             
  nmf_sel <- lapply(names(nmf_programs), function(x) nmf_programs[[x]][,intra_intersect_max[[x]]>=intra_min]) 
  names(nmf_sel) <- names(nmf_programs)
  
  # Select NMF programs based on i) the maximum overlap with other NMF programs from the same sample and
  # ii) the minimum overlap with programs from another sample
  nmf_sel_unlist <- do.call(cbind, nmf_sel)
  inter_intersect <- apply(nmf_sel_unlist , 2, function(x) apply(nmf_sel_unlist , 2, function(y) length(intersect(x,y)))) ## calculating intersection between all programs
  
  final_filter <- NULL 
  for(i in names(nmf_sel)) {
    a <- inter_intersect[grep(i, colnames(inter_intersect), invert = T),grep(i, colnames(inter_intersect))]
    b <- sort(apply(a, 2, max), decreasing = T) # for each sample, ranks programs based on their maximum overlap with programs of other samples
    if(inter_filter==T) b <- b[b>=inter_min] # selects programs with a maximum intersection of at least 10
    if(length(b) > 1) {
      c <- names(b[1]) 
      for(y in 2:length(b)) {
        if(max(inter_intersect[c,names(b[y])]) <= intra_max) c <- c(c,names(b[y])) # selects programs iteratively from top-down. Only selects programs that have a intersection smaller than 10 with a previously selected programs
      }
      final_filter <- c(final_filter, c)
    } else {
      final_filter <- c(final_filter, names(b))
    }
  }
  return(final_filter)                                                      
}


# Custom color palette

custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))

## Create list of NMF matrics where each sample is an entry
path <- "MP/NMF/out_dir/"
sample_ls <- list.files(path)

## Create list of NMF matrics where each sample is an entry
prog_genes_ls <- list()
for(i in seq_along(sample_ls)) {
  nmf_obj <- readRDS(paste(path, sample_ls[[i]], sep = "/"))
  samp_name <- stringr::str_split(sample_ls[[i]], pattern = "_")[[1]][[1]]
  nmf_mats <- c()
  for(j in names(nmf_obj$fit)) {
    get_basis <- basis(nmf_obj$fit[[j]])
    colnames(get_basis)  <- paste0(samp_name, ".", j, ".", 1:j)
    nmf_mats <- cbind(nmf_mats, get_basis)
  }
  prog_genes_ls[[i]] <- nmf_mats
  names(prog_genes_ls)[i] <- paste0(samp_name)
  rm(nmf_obj, nmf_mats, get_basis)
}
Genes_nmf_w_basis <- do.call(c, list(prog_genes_ls))

# Find robust NMFs
# get gene programs (top 50 genes by NMF score)
nmf_programs_sig <- lapply(prog_genes_ls, function(x) apply(x, 2, function(y) names(sort(y, decreasing = T))[1:50]))

# for each sample, select robust NMF programs (i.e. observed using different ranks in the same sample), remove redundancy due to multiple ranks, and apply a filter based on the similarity to programs from other samples. 
nmf_filter_all <- robust_nmf_programs(nmf_programs_sig, intra_min = 35, intra_max = 10, inter_filter=T, inter_min = 10)
nmf_programs_sig <- lapply(nmf_programs_sig, function(x) x[, is.element(colnames(x), nmf_filter_all),drop=F])
nmf_programs_sig <- do.call(cbind, nmf_programs_sig)

#leiden clusters
all_leiden_programs<-readRDS("MP/all_leiden_programs_v2.rds")
all_leiden_programs <- do.call(cbind, all_leiden_programs)
nmf_programs_sig<-cbind(nmf_programs_sig,all_leiden_programs)

# calculate similarity between programs
nmf_intersect <- apply(nmf_programs_sig , 2, function(x) apply(nmf_programs_sig , 2, function(y) length(intersect(x,y)))) 

# hierarchical clustering of the similarity matrix 
nmf_intersect_hc_all <- hclust(as.dist(50-nmf_intersect), method="average") 
nmf_intersect_hc_all <- reorder(as.dendrogram(nmf_intersect_hc_all), colMeans(nmf_intersect))
nmf_intersect         <- nmf_intersect[order.dendrogram(nmf_intersect_hc_all), order.dendrogram(nmf_intersect_hc_all)]

nmf_intersect<-readRDS("MP/NMF/nmf_intersect_124.RDS")
nmf_programs_sig<-readRDS("MP/NMF/nmf_programs_sig_124.RDS")


### use a clustering approach that updates MPs in each iteration 

nmf_intersect_KEEP    <- nmf_intersect
nmf_programs_sig_KEEP <- nmf_programs_sig


### Parameters (later change to function form)v1-keep!

Min_intersect_initial <- 12  # the minimal intersection cutoff for defining the Founder NMF program of a cluster
Min_intersect_cluster <- 12 # the minimal intersection cuttof for adding a new NMF to the forming cluster 
Min_group_size        <- 4     # the minimal group size to consider for defining the Founder_NMF of a MP 

Sorted_intersection       <-  sort(apply(nmf_intersect , 2, function(x) (length(which(x>=Min_intersect_initial))-1)  ) , decreasing = TRUE)

Cluster_list <- list()   ### Every entry contains the NMFs of a chosec cluster
k <- 1
Curr_cluster <- c()
MP_list      <- list()

while (Sorted_intersection[1]>Min_group_size) {   
  
  Curr_cluster <- c(Curr_cluster , names(Sorted_intersection[1]))
  
  ### intersection between all remaining NMFs and Genes in MP 
  Genes_MP                   <- nmf_programs_sig[,names(Sorted_intersection[1])] # initial genes are those in the first NMF. Genes_MP always has only 50 genes consisting of the current MP
  nmf_programs_sig           <- nmf_programs_sig[,-match(names(Sorted_intersection[1]) , colnames(nmf_programs_sig))]  # remove selected NMF
  Intersection_with_Genes_MP <- sort(apply(nmf_programs_sig, 2, function(x) length(intersect(Genes_MP,x))) , decreasing = TRUE) # intersection between all other NMFs and Genes_MP  
  NMF_history                <- Genes_MP  # has all genes in all NMFs in the current cluster, for newly defining Genes_MP after adding a new NMF 
  
  ### Create gene list - composed of intersecting genes in descending order. Update Curr_cluster each time.
  
  while ( Intersection_with_Genes_MP[1] >= Min_intersect_cluster) {   ### Define current cluster 
    
    Curr_cluster  <- c(Curr_cluster , names(Intersection_with_Genes_MP)[1])
    
    Genes_MP_temp   <- sort(table(c(NMF_history , nmf_programs_sig[,names(Intersection_with_Genes_MP)[1]])), decreasing = TRUE)   ## Genes_MP is newly defined each time according to all NMFs in the current cluster 
    Genes_at_border <- Genes_MP_temp[which(Genes_MP_temp == Genes_MP_temp[50])]   ### genes with overlap equal to the 50th gene
    
    if (length(Genes_at_border)>1){
      ### Sort last genes in Genes_at_border according to maximal NMF gene scores
      ### Run over all NMF programs in Curr_cluster and extract NMF scores for each gene
      Genes_curr_NMF_score <- c()
      for (i in Curr_cluster) {
        curr_study           <- strsplit(i, "[.]")[[1]][[1]]
        Q                    <- Genes_nmf_w_basis[[curr_study]][ match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]])))[!is.na(match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]]))))]   ,i] 
        #names(Q)             <- names(Genes_at_border[!is.na(match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]]))))])  ### sometimes when adding genes the names do not appear 
        Genes_curr_NMF_score <- c(Genes_curr_NMF_score,  Q )
      }
      Genes_curr_NMF_score_sort <- sort(Genes_curr_NMF_score , decreasing = TRUE)
      Genes_curr_NMF_score_sort <- Genes_curr_NMF_score_sort[unique(names(Genes_curr_NMF_score_sort))]   ## take only the maximal score of each gene - which is the first entry after sorting
      
      Genes_MP_temp <- c(names(Genes_MP_temp[which(Genes_MP_temp > Genes_MP_temp[50])]) , names(Genes_curr_NMF_score_sort))
      
    } else {
      Genes_MP_temp <- names(Genes_MP_temp)[1:50] 
    }
    
    NMF_history   <- c(NMF_history , nmf_programs_sig[,names(Intersection_with_Genes_MP)[1]]) 
    Genes_MP <- Genes_MP_temp[1:50]
    
    nmf_programs_sig      <- nmf_programs_sig[,-match(names(Intersection_with_Genes_MP)[1] , colnames(nmf_programs_sig))]  # remove selected NMF
    
    Intersection_with_Genes_MP <- sort(apply(nmf_programs_sig, 2, function(x) length(intersect(Genes_MP,x))) , decreasing = TRUE) # intersection between all other NMFs and Genes_MP  
    
  }
  
  Cluster_list[[paste0("Cluster_",k)]] <- Curr_cluster
  MP_list[[paste0("MP_",k)]]           <- Genes_MP
  k <- k+1
  
  
  nmf_intersect             <- nmf_intersect[-match(Curr_cluster,rownames(nmf_intersect) ) , -match(Curr_cluster,colnames(nmf_intersect) ) ]  # remove current chosen cluster
  
  Sorted_intersection       <-  sort(apply(nmf_intersect , 2, function(x) (length(which(x>=Min_intersect_initial))-1)  ) , decreasing = TRUE)   ### Sort intersection of remaining NMFs not included in any of the previous clusters
  
  Curr_cluster <- c()
  print(dim(nmf_intersect)[2])
}

#### *****  Sort Jaccard similarity plot according to new clusters:

inds_sorted <- c()

for (j in 1:length(Cluster_list)){
  
  inds_sorted <- c(inds_sorted , match(Cluster_list[[j]] , colnames(nmf_intersect_KEEP)))
  
}
inds_new <- c(inds_sorted   ,   which(is.na( match(1:dim(nmf_intersect_KEEP)[2],inds_sorted)))) ### combine inds in clusters with inds without clusters


# plot re-ordered similarity matrix heatmap     
nmf_intersect_meltI_NEW <- reshape2::melt(nmf_intersect_KEEP[inds_new,inds_new]) 

custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))


ggplot(data = nmf_intersect_meltI_NEW, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
  geom_tile() + 
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +                                
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))

ggplot(data = nmf_intersect_meltI_NEW, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
  geom_tile() + 
  scale_color_gradient2(limits=c(2,30), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 16, oob=squish, name="Similarity\n(Jaccard index)") +                                
  scale_fill_gradient2(limits=c(2,30), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 16, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))


MP <-  do.call(cbind, MP_list)

names(MP_list)=c("Neuron","Vasc","MES.Hyp","Mac","OPC.AC","Oligo","LQ.Chromatin.reg","MES","Prolif.Metab","MES.Ast","Reactive.Ast","NPC","Inflammatory.Mac")


###### RECLUSTERING/HEIRARCHICAL CLUSTERING, EXTENDED MPs, AND SUBCLUSTER PROGRAMS ################
programs_ls<-readRDS("MP/NMF/programs_ls_124.rds") #this is "nmf_programs_sig"
Cluster_list<-readRDS("MP/NMF/Cluster_list_124.rds")
str(Cluster_list)
MP_list<-readRDS("MP/NMF/combined_gbm_metaprograms_raw_124.rds")


#reclustering OPC.AC into 2 clusters
opc.ac <- Cluster_list[[5]]

mat2list <- function(x) {
  stopifnot(is.matrix(x))
  lapply(seq_len(ncol(x)), function(i) x[, i])
}
program_ls<-mat2list(programs_ls)
names(program_ls)<-colnames(programs_ls)
opc.ac_mat <- program_ls[names(program_ls) %in% opc.ac]

jac_mat <- scalop::jaccard(opc.ac_mat)
hc <- hclust(dist(1 - jac_mat), method = "average")
jac_plot <- melt(jac_mat[hc$order, hc$order])
#summary(jac_plot$value)

p1 <- ggplot(data = jac_plot, aes(x = Var1, y = Var2, fill = value)) + geom_raster() + 
  scale_fill_gradient2(limits = c(0, 0.4), low = rev(magma(323, begin = 0.15)), mid = "antiquewhite", high = rev(magma(323, begin = 0.18)), midpoint = 0 , oob = squish, name = "Pearson\nCorrelation") + 
  theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),  axis.line = element_blank(),axis.title = element_text(size=24), legend.text.align = 0.5) + 
  scale_x_discrete(name = "\nPrograms", breaks = colnames(jac_mat)[hc$order][seq(5, ncol(jac_mat), by = 5)], labels = seq(5, ncol(jac_mat), by = 5)) + 
  scale_y_discrete(name = "\nPrograms", breaks = colnames(jac_mat)[hc$order][seq(5, ncol(jac_mat), by = 5)], labels = seq(5, ncol(jac_mat), by = 5))  

prog_ordered <- colnames(jac_mat[hc$order, hc$order])
plot(hc)
rect.hclust(hc , k = 2, border = 2:6)
clusterCut <- stats::cutree(tree = hc, k = 2)

metaprog_df <- data.frame(row.names = prog_ordered)
metaprog_df$cluster <- (clusterCut[prog_ordered])
metaprog_df$cells <- prog_ordered
metaprog_df$cluster <- as.factor(metaprog_df$cluster)

metaprog_annotation_plot <- ggplot(metaprog_df, aes(x = factor(cells, levels = cells), y = 1)) + 
  geom_raster(aes(fill = cluster)) + theme(legend.position = "bottom", axis.ticks = element_blank(),axis.title = element_text(size=24), panel.background = element_rect(fill = "white"), axis.text = element_blank(), axis.line = element_blank(), axis.text.x = element_blank()) + 
  labs(x = "", y = "", fill = "") 

egg::ggarrange(p1, metaprog_annotation_plot, ncol = 1, nrow = 2, heights = c(40, 2))

clust_list_ac_opc <- list(names(clusterCut[clusterCut == 1]), names(clusterCut[clusterCut == 2])) # use this ordering to reorder cluster 5 in main MP heatmap
mp_freqs <- sapply(clust_list_ac_opc, function(k) sort(table(unlist(opc.ac_mat[k])), decreasing = T), simplify = F)
opc.ac_metaprograms <- sapply(mp_freqs, function(tab) head(names(tab)[tab >= 3], 50), simplify = F)
names(opc.ac_metaprograms)<-c("AC","OPC")


#dividing low quality from chromatin.reg cluster
chromatin <- Cluster_list[[7]]
mat2list <- function(x) {
  stopifnot(is.matrix(x))
  lapply(seq_len(ncol(x)), function(i) x[, i])
}
program_ls<-mat2list(programs_ls)
names(program_ls)<-colnames(programs_ls)
chromatin_mat <- program_ls[names(program_ls) %in% chromatin]

jac_mat <- scalop::jaccard(chromatin_mat)
hc <- hclust(dist(1 - jac_mat), method = "average")
jac_plot <- melt(jac_mat[hc$order, hc$order])
#summary(jac_plot$value)

p1 <- ggplot(data = jac_plot, aes(x = Var1, y = Var2, fill = value)) + geom_raster() + 
  scale_fill_gradient2(limits = c(0, 0.4), low = rev(magma(323, begin = 0.15)), mid = "antiquewhite", high = rev(magma(323, begin = 0.18)), midpoint = 0 , oob = squish, name = "Pearson\nCorrelation") + 
  theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),  axis.line = element_blank(),axis.title = element_text(size=24), legend.text.align = 0.5) + 
  scale_x_discrete(name = "\nPrograms", breaks = colnames(jac_mat)[hc$order][seq(5, ncol(jac_mat), by = 5)], labels = seq(5, ncol(jac_mat), by = 5)) + 
  scale_y_discrete(name = "\nPrograms", breaks = colnames(jac_mat)[hc$order][seq(5, ncol(jac_mat), by = 5)], labels = seq(5, ncol(jac_mat), by = 5))  


prog_ordered <- colnames(jac_mat[hc$order, hc$order])
plot(hc)
rect.hclust(hc , k = 3, border = 2:6)
clusterCut <- stats::cutree(tree = hc, k = 3)

metaprog_df <- data.frame(row.names = prog_ordered)
metaprog_df$cluster <- (clusterCut[prog_ordered])
metaprog_df$cells <- prog_ordered
metaprog_df$cluster <- as.factor(metaprog_df$cluster)

metaprog_annotation_plot <- ggplot(metaprog_df, aes(x = factor(cells, levels = cells), y = 1)) + 
  geom_raster(aes(fill = cluster)) + theme(legend.position = "bottom", axis.ticks = element_blank(),panel.background = element_rect(fill = "white"), axis.text = element_blank(), axis.line = element_blank(), axis.text.x = element_blank(),axis.title = element_text(size=16)) + 
  labs(x = "", y = "", fill = "") 

egg::ggarrange(p1, metaprog_annotation_plot, ncol = 1, nrow = 2, heights = c(40, 2))

clust_list <- list(names(clusterCut[clusterCut == 1]), names(clusterCut[clusterCut == 2]),names(clusterCut[clusterCut == 3]))
mp_freqs <- sapply(clust_list, function(k) sort(table(unlist(chromatin_mat[k])), decreasing = T), simplify = F)
chromatin_metaprograms <- sapply(mp_freqs, function(tab) head(names(tab)[tab >= 3], 50), simplify = F)
names(chromatin_metaprograms)<-c("LQ1","chromatin_clean","LQ2")
chromatin_metaprogram_clean<-chromatin_metaprograms$chromatin_clean


clust_list2 <- list(names(clusterCut[clusterCut == 1]), names(clusterCut[clusterCut == 3]),names(clusterCut[clusterCut == 2]))
#this is the order to use for the final MP heatmap
clust_list_chromatin <- c(names(clusterCut[clusterCut == 1]), names(clusterCut[clusterCut == 2]),names(clusterCut[clusterCut == 3]))

###### MPs after OP/AC split and LQ/chromatin.reg split#######
MP_list2<-MP_list[-c(5,7)]

MP_list2[["chromatin.reg"]] <- chromatin_metaprogram_clean
MP_list2[["OPC"]] <- opc.ac_metaprograms$OPC
MP_list2[["AC"]] <- opc.ac_metaprograms$AC


###### SPOT ASSIGNMENTS TO CLEANED METAPROGRAMS##########

# generate normalized exp matrices
file_paths <- list.files(path ="general/exp_mats_GBM", pattern = "\\.rds", full.names = TRUE)
sample_ls <-  gsub(pattern = "\\.rds$", replacement = "", x = basename(file_paths))
sample_ls->samples_names
sample_ls->samples
per_sample_mat <- lapply(file_paths, readRDS)

for (i in seq_along(per_sample_mat)){ 
  m <- as.matrix(per_sample_mat[[i]])
  m <- m[-grep("^MT-|^RPL|^RPS", rownames(m)), ]
  if(min(colSums(m)) == 0){m <- m[, colSums(m) != 0]}
  scaling_factor <- 1000000/colSums(m)
  m_CPM <- sweep(m, MARGIN = 2, STATS = scaling_factor, FUN = "*")
  m_loged <- log2(1 + (m_CPM/10))
  
  # removing genes with zero variance across all cells
  var_filter <- apply(m_loged, 1, var)
  m_proc <- m_loged[var_filter != 0, ]
  # filtering out lowly expressed genes
  exp_genes <- rownames(m_proc)[(rowMeans(m_proc) > 0.4)]
  m_proc <- m_proc[exp_genes, ]
  
  # output to a list of gene expression profiles (GEP)
  per_sample_mat[[i]] <- m_proc
  names(per_sample_mat)[i] <- sample_ls[[i]]
  rm(m,m_loged, var_filter, exp_genes, m_proc)
}

#generate a list with the filtered exp matrix for each sample, i.e. m_proc<-(per_sample_mat[[1]])
score_mat <- lapply(c(1:length(per_sample_mat)), function(i){
  m_proc<-per_sample_mat[[i]]
  metaprograms_gene_list <- readRDS("MP/clean_spatial_gbm_metaprograms_124.rds")
  signatures <- scalop::sigScores(m_proc, metaprograms_gene_list, expr.center = TRUE, conserved.genes = 0.5)
  
  spot_scores <- data.frame(spot_names = rownames(signatures))
  spot_scores$Neuron <- signatures$Neuron
  spot_scores$Vasc<- signatures$Vasc
  spot_scores$MES.Hyp <- signatures$MES.Hyp
  spot_scores$Mac <- signatures$Mac
  spot_scores$Oligo <- signatures$Oligo
  spot_scores$MES <- signatures$MES
  spot_scores$Prolif.Metab <- signatures$Prolif.Metab
  spot_scores$MES.Ast <- signatures$MES.Ast
  spot_scores$Reactive.Ast <- signatures$Reactive.Ast
  spot_scores$NPC <- signatures$NPC
  spot_scores$Inflammatory.Mac <- signatures$Inflammatory.Mac
  spot_scores$chromatin.reg <- signatures$chromatin.reg
  spot_scores$OPC <- signatures$OPC
  spot_scores$AC <- signatures$AC
  return(spot_scores)
})

names(score_mat) <- sample_ls

for (i in seq_along(score_mat)){
  score_df <- as.data.frame(score_mat[[i]])
  score_df <-column_to_rownames(score_df, 'spot_names')
  maxcol_meta<-maxcol_strict(score_df)
  maxcol_meta<-stack(maxcol_meta)
  setnames(maxcol_meta,2,"spot_type_meta_new")
  setnames(maxcol_meta,1,"SpotID")
}
