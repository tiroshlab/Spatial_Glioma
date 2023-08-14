library(tidyverse)
library(scalop)
library(patchwork)
library(colorspace)
library(ggpubr)
library(Matrix)
library(biomaRt)

# Run functions first

# load data ---------------------------------------------------------------


# Load Data (counts mat exported from Seurat post-spot filtering:

samples_names <- (read.delim("CNA/cna_samples.txt", header = FALSE))$V1
sample_ls <- (read.delim("general/GBM_samples.txt", header = FALSE))$V1

norm_brain_ref1<-readRDS("CNA/normal_brain_ref/UKF256_C.rds")
colnames(norm_brain_ref1) <- paste(colnames(norm_brain_ref1), "_ref1", sep="") 

norm_brain_ref2<-readRDS("CNA/normal_brain_ref/UKF265_C.rds")
colnames(norm_brain_ref2) <- paste(colnames(norm_brain_ref2), "_ref2", sep="") 


samples_regions <- read.delim("CNA/samples_regions.txt", header = TRUE)


# merged samples cna --------------------------------------------------------

samples_num <- c(15:18)
merge_cna <- mclapply(samples_num, merge_cna_md_fx, mc.cores = 4)

merge_cna_md_fx <- function(i){
  m1 <- readRDS(paste("general/exp_mats_GBM/", samples_names[i], "counts.rds", sep = ""))
  m <- cbind(m1,norm_brain_ref1,norm_brain_ref2)
  
  scaling_factor <- 1000000/colSums(m) #tumor mat processing
  m_CPM <- sweep(m, MARGIN = 2, STATS = scaling_factor, FUN = "*")
  m_loged <- log2(1 + (m_CPM/10))
  var_filter <- apply(m_loged, 1, var) # removing genes with zero variance across all cells
  m_proc <- m_loged[var_filter != 0, ]
  exp_genes <- rownames(m_proc)[(rowMeans(m_proc) > 0.4)] # filtering out lowly expressed genes
  m_proc <- m_proc[exp_genes, ]
  
  # create normal ref
  ref <-list(colnames(norm_brain_ref1), colnames(norm_brain_ref2))
  query <- colnames(m_proc)[!is.element(colnames(m_proc), unname(unlist(ref)))]
  
  # create metadata table 
  
  
  metadata <- tibble(CellID = colnames(m_proc), 
                     org_ref = ifelse(colnames(m_proc) %in% unname(unlist(ref)), "reference", "not reference"))
  
  metadata$sample <- "control"
  metadata$region <- "control"
  
  if (str_detect(samples_names[i], "merge")) {
    sapply(c(1:length(samples_regions$sample[samples_regions$cna_samples_name == samples_names[i]])), function(s){
      s_name <- samples_regions$sample[samples_regions$cna_samples_name == samples_names[i]][s]
      s_region <- samples_regions$region[samples_regions$cna_samples_name == samples_names[i]][s]
      metadata$sample <<- ifelse(str_detect(metadata$CellID, s_name), s_name, metadata$sample)
      metadata$region <<- ifelse(str_detect(metadata$CellID, s_name), s_region, metadata$region)
    })
  } else {
    metadata$sample <- ifelse(metadata$org_ref == "not reference", samples_regions$sample[samples_regions$cna_samples_name == samples_names[i]], metadata$sample)
    metadata$region <- ifelse(metadata$org_ref == "not reference", samples_regions$region[samples_regions$cna_samples_name == samples_names[i]], metadata$region)
    
  }
  
  
  # calc CNA
  hg19<- readRDS("CNA/hg38.rds")
  cna_score_mat <- calc_cna(matrix = m_proc, query = query, ref = ref, genome = hg19, range = c(-3,3), window = 150, noise = 0.15, isLog = TRUE, per_chr = TRUE, scale = NULL, top_genes = NULL, verbose = TRUE)
  sig_and_cor <- cna_sig_cor(cna_score_mat, epi_cells = query, ref_cells = ref, cna_sig = "abs", top_region = 1/3, top_cells = 1/3)
  metadata$CNAsig_org <- sig_and_cor$cna_signal[match(metadata$CellID, names(sig_and_cor$cna_signal))]
  metadata$CNAcor_org <- sig_and_cor$cna_correlation[match(metadata$CellID, names(sig_and_cor$cna_correlation))]
  
  # update reference with combined external ref + tumor defined by cna sig / corr thresholds
  metadata$new_ref <- ifelse((metadata$CNAcor_org <= 0.25 & metadata$CNAsig_org <= 0.13) | (metadata$org_ref == "reference"), "reference", "query")
  
  ref <- metadata$CellID[metadata$new_ref == "reference"]
  query <- colnames(m_proc)[!is.element(colnames(m_proc), unname(unlist(ref)))]
  
  # calc CNA with updated ref
  cna_score_mat <- calc_cna(matrix = m_proc, query = query, ref = ref, genome = hg19, range = c(-3,3), window = 150, noise = 0.15, isLog = TRUE, per_chr = TRUE, scale = NULL, top_genes = NULL, verbose = TRUE)
  sig_and_cor <- cna_sig_cor(cna_score_mat, epi_cells = query, ref_cells = ref, cna_sig = "abs", top_region = 1/3, top_cells = 1/3)
  metadata$CNAsig <- sig_and_cor$cna_signal[match(metadata$CellID, names(sig_and_cor$cna_signal))]
  metadata$CNAcor <- sig_and_cor$cna_correlation[match(metadata$CellID, names(sig_and_cor$cna_correlation))]
  
  metadata$CNAtot <- metadata$CNAcor + (metadata$CNAsig * max(metadata$CNAcor) / max(metadata$CNAsig))
  
  fin_cna <- cna_score_mat[,metadata$CellID[metadata$org_ref == "not reference"]]
  md <- metadata[metadata$org_ref == "not reference",]
  
  return(list(md,fin_cna))
}


# all cna + metadata --------------------------------------------------------

samples_num <- c(1:26)
all_cna <- mclapply(samples_num, all_cna_fx, mc.cores = 26)

all_cna_fx <- function(i){
  m1 <- readRDS(paste("general/exp_mats_GBM/", samples_names[i], "counts.rds", sep = ""))
  m <- cbind(m1,norm_brain_ref1,norm_brain_ref2)
  
  scaling_factor <- 1000000/colSums(m) #tumor mat processing
  m_CPM <- sweep(m, MARGIN = 2, STATS = scaling_factor, FUN = "*")
  m_loged <- log2(1 + (m_CPM/10))
  var_filter <- apply(m_loged, 1, var) # removing genes with zero variance across all cells
  m_proc <- m_loged[var_filter != 0, ]
  exp_genes <- rownames(m_proc)[(rowMeans(m_proc) > 0.4)] # filtering out lowly expressed genes
  m_proc <- m_proc[exp_genes, ]
  
  ref <-list(colnames(norm_brain_ref1), colnames(norm_brain_ref2))
  query <- colnames(m_proc)[!is.element(colnames(m_proc), unname(unlist(ref)))]
  
  # create metadata table 
  metadata <- tibble(CellID = colnames(m_proc), 
                     org_ref = ifelse(colnames(m_proc) %in% unname(unlist(ref)), "reference", "not reference"))
  
  metadata$sample <- "control"
  
  metadata$sample <- ifelse(metadata$org_ref == "not reference", samples_names[i], metadata$sample)
  
  
  # calc CNA
  hg19<- readRDS("CNA/hg38.rds")
  cna_score_mat <- calc_cna(matrix = m_proc, query = query, ref = ref, genome = hg19, range = c(-3,3), window = 150, noise = 0.15, isLog = TRUE, per_chr = TRUE, scale = NULL, top_genes = NULL, verbose = TRUE)
  sig_and_cor <- cna_sig_cor(cna_score_mat, epi_cells = query, ref_cells = ref, cna_sig = "abs", top_region = 1/3, top_cells = 1/3)
  metadata$CNAsig_org <- sig_and_cor$cna_signal[match(metadata$CellID, names(sig_and_cor$cna_signal))]
  metadata$CNAcor_org <- sig_and_cor$cna_correlation[match(metadata$CellID, names(sig_and_cor$cna_correlation))]
  
  # update reference with combined external ref + tumor defined by cna sig / corr thresholds
  metadata$new_ref <- ifelse((metadata$CNAcor_org <= 0.25 & metadata$CNAsig_org <= 0.13) | (metadata$org_ref == "reference"), "reference", "query")
  
  ref <- metadata$CellID[metadata$new_ref == "reference"]
  query <- colnames(m_proc)[!is.element(colnames(m_proc), unname(unlist(ref)))]
  
  # calc CNA with updated ref
  cna_score_mat <- calc_cna(matrix = m_proc, query = query, ref = ref, genome = hg19, range = c(-3,3), window = 150, noise = 0.15, isLog = TRUE, per_chr = TRUE, scale = NULL, top_genes = NULL, verbose = TRUE)
  sig_and_cor <- cna_sig_cor(cna_score_mat, epi_cells = query, ref_cells = ref, cna_sig = "abs", top_region = 1/3, top_cells = 1/3)
  metadata$CNAsig <- sig_and_cor$cna_signal[match(metadata$CellID, names(sig_and_cor$cna_signal))]
  metadata$CNAcor <- sig_and_cor$cna_correlation[match(metadata$CellID, names(sig_and_cor$cna_correlation))]
  
  metadata$CNAtot <- metadata$CNAcor + (metadata$CNAsig * max(na.omit(metadata$CNAcor)) / max(na.omit(metadata$CNAsig)))
  
  fin_cna <- cna_score_mat[,metadata$CellID[metadata$org_ref == "not reference"]]
  md <- metadata[metadata$org_ref == "not reference",]
  return(list(md,fin_cna))
  
}


# plot combined CNA ----------------------------------------------------------------


c1 <- as.data.frame(all_cna[[1]][[2]])
c2 <- as.data.frame(all_cna[[2]][[2]])

genes_cna <- intersect(row.names(c1),row.names(c2))

sapply(c(3:26), function(i){
  ci <- as.data.frame(all_cna[[i]][[2]])
  genes_cna <<- intersect(genes_cna, row.names(ci))
})

set.seed(50)
samp_spots <- sapply(c(1:26), function(i){
  m <- all_cna[[i]][[1]]
  return(sample(m$CellID[m$sample == samples_names[i]],700, replace = F))
})

cna_comb <- c1[genes_cna,as.character(samp_spots[,1])]
colnames(cna_comb) <- paste(samples_names[1], colnames(cna_comb), sep = "_")

sapply(c(2:26), function(i){
  ci_temp <- as.data.frame(all_cna[[i]][[2]])
  ci <- ci_temp[genes_cna,as.character(samp_spots[,i])]
  colnames(ci) <- paste(samples_names[i], colnames(ci), sep = "_")
  cna_comb <<- cbind(cna_comb, ci)
})


spots_list <- lapply(c(1:26), function(i){
  change_spots <- paste(samples_names[i], samp_spots[,i], sep = "_")
  return(change_spots)
})

names(spots_list) <- samples_names

find_chr <- genome_break(row.names(cna_comb))
sig7 <- apply(cna_comb[1772:2006,],2,mean)
sig10 <- apply(cna_comb[2363:2538,],2,mean)

spots_ord <- sapply(spots_list, function(sp){
  return(mean(sig7[sp])+abs(mean(sig10[sp])))
})
names(spots_ord) <- names(spots_list)

spots_list <- spots_list[names(sort(spots_ord))]

cna_pl <- infercna::cnaPlot(cna_comb, limit=c(-0.7,0.7), order.cells = spots_list, ratio = NULL)   # plot cna - important that the ratio is set to NULL!

cna_pl$p + theme(legend.position = "right")

old_spots_list <- spots_list[rev(c("ZH1019T1","UKF259","UKF255","UKF248","ZH1007nec","UKF251","MGH258","UKF243","ZH1019inf","UKF275","UKF260","UKF313","ZH916T1",
                                   "ZH1007inf","UKF334","UKF269","ZH8811Bbulk","ZH8811Abulk","ZH916bulk","UKF296","ZH881T1","UKF304","ZH8812bulk","UKF266","ZH881inf","ZH916inf"))]


# purity and binning ------------------------------------------------------------------


all_tot <- sapply(c(1:14), function(i){
  md <- all_cna[[i]][[1]]
  return(md$CNAtot)
})
names(all_tot) <- samples_names[1:14]

all_tot_merge <- sapply(c(1:4), function(i){
  md <- merge_cna[[i]][[1]]
  return(md$CNAtot)
})
names(all_tot_merge) <- samples_names[15:18]

hist(c(unlist(all_tot), unlist(all_tot_merge)), breaks = 150)

all_tot_split_tmp <- sapply(names(all_tot_merge), function(m){
  tmp <- all_tot_merge[[m]]
  split_samps <- samples_regions$sample[samples_regions$cna_samples_name == m]
  split_scores <- sapply(split_samps, function(s){
    tmp_split <- tmp[grepl(s,names(tmp))]
    names(tmp_split) <- sapply(strsplit(names(tmp_split),"_"),function(ns){ns[2]}) 
    return(tmp_split)
  })
  return(split_scores)
})

all_tot_split <- unlist(all_tot_split_tmp, recursive = F)
names(all_tot_split) <- sapply(strsplit(names(all_tot_split),"\\."),function(ns){ns[2]})

purity_score <- c(all_tot,all_tot_split)
purity_score <- purity_score[sample_ls]

min_vl <- min(unlist(purity_score))
max_vl <- max(unlist(purity_score))

purity_score_scaled <- sapply(purity_score, function(p){
  p_sc <- (p - min_vl)/(max_vl-min_vl)
  return(p_sc)
})


#  malignancy binning ------------------------------------

hist(c(purity_score_scaled[[15]], purity_score_scaled[[23]]), breaks = 100)
purity_th <- quantile(c(purity_score_scaled[[15]], purity_score_scaled[[23]]))

malignancy_bin <- sapply(purity_score_scaled, function(p){
  bin <- ifelse(p < purity_th[2] , "non_malignant", 
                ifelse(p >= purity_th[2] & p < purity_th[3], "mix_low", 
                       ifelse(p >= purity_th[3]  & p < purity_th[4],"mix_high","malignant")))
  return(bin)
})


# zh1019 ------------------------------------------------------------------

# plot samples separately 

cna_pl <- infercna::cnaPlot(all_cna[[23]][[2]], limit=c(-0.7,0.7), order.cells = NULL, ratio = NULL)   # plot cna - important that the ratio is set to NULL!
cna_pl$p + theme(legend.position = "left")

cna_pl <- infercna::cnaPlot(all_cna[[15]][[2]], limit=c(-0.7,0.7), order.cells = NULL, ratio = NULL)   # plot cna - important that the ratio is set to NULL!
cna_pl$p + theme(legend.position = "left")

# plot samples together

zh1019inf_bin <- malignancy_bin[[15]]
zh1019t1_bin <- malignancy_bin[[23]]
names(zh1019inf_bin) <- paste("ZH1019inf_", names(zh1019inf_bin), sep = "")
names(zh1019t1_bin) <- paste("ZH1019T1_", names(zh1019t1_bin), sep = "")
zh1019_bin <- c(zh1019inf_bin, zh1019t1_bin)

zh1019inf_purity <- purity_score_scaled[[15]]
zh1019t1_purity <- purity_score_scaled[[23]]
names(zh1019inf_purity) <- paste("ZH1019inf_", names(zh1019inf_purity), sep = "")
names(zh1019t1_purity) <- paste("ZH1019T1_", names(zh1019t1_purity), sep = "")
zh1019_purity <- c(zh1019inf_purity, zh1019t1_purity)
spots_by_purity <- names(sort(zh1019_purity))

zh1019_md <- merge_cna[[1]][[1]]
zh1019_cna <- merge_cna[[1]][[2]]
ord_cna <- zh1019_cna[,spots_by_purity]

least_pure_sig <- apply(ord_cna[,1:259], 1, mean)
ord_cna_clean <- ord_cna - least_pure_sig

cna_pl <- infercna::cnaPlot(ord_cna_clean, limit=c(-0.7,0.7), order.cells = F, ratio = NULL)   # plot cna - important that the ratio is set to NULL!
anno_plot <- cna_pl$p + theme(legend.position = "right")

samp_col<- c(ZH1019inf="palegreen1", ZH1019T1 ="darkgreen")
bin_col<- c(malignant="#7c1d6f", mix_high = "#C77CFF", mix_low="#F2B701", non_malignant="#fcde9c")
annotate_sample <- ggbar(zh1019_md$sample[match(zh1019_md$CellID,colnames(ord_cna_clean))], dir = "v",legend_title = "sample", cols = samp_col) + theme(legend.direction = "vertical")
annotate_purity <- ggbar(as.character(zh1019_bin[colnames(ord_cna_clean)]), dir = "v", legend_title = "purity binning", cols = bin_col) + theme(legend.direction = "vertical")

anno_plot + annotate_sample + theme(legend.position = "none") + annotate_purity + theme(legend.position = "none") + plot_layout(ncol = 4, widths = c(1, 0.05, 0.05, 0.05), guides = "collect")

# extract the legends to paste them later in power-point :)
lej1 <- cowplot::get_legend(annotate_sample)
lej2 <- cowplot::get_legend(annotate_purity)
grid.newpage()
grid.draw(lej1)
grid.newpage()
grid.draw(lej2)


######## Functions ---------------------------------------------------------------
# CNA Score Calculation ---------------------------------------------------


#### calc_cna - calculate CNA scores per cell/spot given a gene x cell matrix, query cells and reference cell list/vector
# NOTE: matrix should NOT be row(gene)-centered
calc_cna <- function(matrix,
                     query,
                     ref,
                     top_genes = NULL,
                     window = 100,
                     range = c(-3, 3),
                     per_chr = TRUE,
                     scale = 0.05,
                     noise = NULL,
                     isLog = FALSE,
                     genome = NULL,
                     min_ref_leng = 10,
                     verbose = FALSE){
  
  if(all(round(range(rowMeans(matrix)), 3) == 0)) {
    stop(print("Matrix is row-centered. Please provide non-centered data."))
  }
  if(is.list(ref)){
    # Remove small references
    ref_leng <- lapply(ref, length)
    ref <- ref[ref_leng >= min_ref_leng]
    # Prepare CNA matrix to work on
    cna_mat <- matrix[, c(query, unlist(ref))]
  }
  else if(!is.list(ref)){
    if(length(ref) < min_ref_leng){
      stop(print("There are not enough reference cells!"))
    }
    cna_mat <- matrix[, c(query, ref)]
  }
  
  # Order chromosomes
  if (verbose) message("Ordering the genes by their genomic position.")
  order_genes <- genome_sort(rownames(cna_mat), genome = genome)
  genes <- order_genes$cna_genes
  chr_breaks <- order_genes$chr_breaks
  
  # Optional list of top expressed genes
  if(!is.null(top_genes)){
    if(verbose) message("Filtering the expression matrix to include only top ", length(top_genes), "genes...")
    if(isTRUE(isLog)){
      cna_mat <- un_log(cna_mat)}
    cna_mat <- cna_mat[genes, ]
    cna_mat <- apply(cna_mat, 1, mean)
    cna_mat <- names(tail(sort(cna_mat), n = top_genes))
    order_genes <- genome_sort(cna_mat, genome = genome)
    genes <- order_genes$cna_genes
    chr_breaks <- order_genes$chr_breaks
  }  
  
  # Reorder
  ordered_mat <- cna_mat[genes, ]
  # Log before first centering
  if(isFALSE(isLog)){
    ordered_mat <- log_norm(ordered_mat)}
  # First row centering step
  if (verbose) message("Performing mean-centering of the genes.")
  avg <- apply(ordered_mat, 1, mean)
  ordered_mat <- sweep(ordered_mat, 1, avg)
  # Set 3 and -3 as extreme values (as set by the argument "range")
  if (verbose) message("Restricting expression matrix values to between ", range[[1]], " and ", range[[2]], ".")
  ordered_mat <- apply(ordered_mat, 2, function(x) pmax(x, range[1]))
  ordered_mat <- apply(ordered_mat, 2, function(x) pmin(x, range[2]))
  # Unlog to CPM/TPM again for moving average
  ordered_mat <- un_log(ordered_mat)
  
  # Calculate moving average per chromosome by window of 100 (as set by the argument "window")
  if(isTRUE(per_chr)){
    if (verbose) message("Calculating rolling means with a window size of ", window, " genes, on each chromosome in turn.")
    num <- seq(1:(length(chr_breaks) -1))
    perchr <- lapply(num, function(y){
      if(y == length(num)){
        end <- nrow(ordered_mat)
      }
      if(y != length(num)){
        end <- chr_breaks[y + 1] - 1
      }
      chr <- ordered_mat[chr_breaks[y]:end, ]
      chr_mat <- apply(chr, 2, function(x) caTools::runmean(x, k = window, endrule = "mean"))
    })
    calc_mat <- do.call(rbind, perchr)  
  }
  # Calculate moving average for all genes as one chromosome
  if(isFALSE(per_chr)){
    if (verbose) message("Calculating rolling means with a window size of ", window, " genes, on all genes as one chromosome.")
    calc_mat <- apply(ordered_mat, 2, function(x) caTools::runmean(x, k = window, endrule = "mean"))
  }
  
  # Log before second centering
  if (verbose) message("Converting CNA score values to log(2) space.")
  calc_mat <- log_norm(calc_mat)
  # Substract median per cell
  if (verbose) message("Performing median-centering of the cells.")
  cell_med <- apply(calc_mat, 2, median)
  calc_mat <- sweep(calc_mat, 2, cell_med)
  # Unlog to CPM/TPM again for reference removal
  calc_mat <- un_log(calc_mat)
  
  # Create max/min values per gene from reference cells
  if(is.list(ref)){
    ref_leng <- seq(1, length(ref))
    mean_mat <- lapply(ref_leng, function(x){
      idx_ref <- ref[[x]]
      m1 <- apply(calc_mat[, idx_ref], 1, mean)
      m1
    })
    ref_mat <- do.call(cbind, mean_mat)
    # Log references
    ref_mat <- log_norm(ref_mat)
    ref_max <- apply(ref_mat, 1, max)
    ref_min <- apply(ref_mat, 1, min)
  }
  else if(!is.list(ref)){
    ref_mat <- apply(calc_mat[, ref], 1, mean)
    ref_mat <- log_norm(ref_mat)
    ref_max <- ref_mat
    ref_min <- ref_mat
  }
  
  # Expand reference boundaries by scaling percentage
  if(!is.null(scale)){
    rmax <- ref_max + scale * abs(ref_max)
    rmin <- ref_min - scale * abs(ref_min)
  }
  # Or expand by fixed noise factor
  if(!is.null(noise)){
    rmax <- ref_max + noise
    rmin <- ref_min - noise 
  }
  
  # Log CNA matrix
  calc_mat <- log_norm(calc_mat)
  # Centre by reference
  if (verbose) message("Correcting CNA profiles using CNA values from reference cells.")
  score_mat <- ifelse(calc_mat > rmax, calc_mat - rmax,
                      ifelse(calc_mat < rmin, calc_mat - rmin, 0))
  rownames(score_mat) <- rownames(ordered_mat)
  
  if (verbose) message("Done!")
  return(score_mat)
}



# CNA Matrix Subclones ----------------------------------------------------


#### Dived the CNA matrix to subclones by reducing feature dimension and clustering with Louvain:
spatial_subclones <- function(cna_matrix,
                              epi_cells,
                              separate = c("arm", "chr"),
                              genome = NULL,
                              genecut = 10,
                              top_region = 1/3,
                              top_method = c("abs", "sd"),
                              reduction_dims = 30,
                              cluster_method = c("louvain", "dbscan"),
                              cluster_k = 10,
                              dbs_minpts = NULL,
                              dbs_ptscale = NULL,
                              diffcut = 10){
  
  ## Subset by chromosome
  breaks <- genome_break(genes = rownames(cna_matrix), separate = separate,
                         genecut = genecut, genome = genome)
  break_idx <- breaks$break_idx
  labels <- breaks$breaks_df$labels
  dispose <- breaks$dispose
  if(length(dispose) != 0) {
    matrix <- cna_matrix[-dispose, ]
  } else {matrix <- cna_matrix}
  break_num <- unique(break_idx)
  
  ## Set working matrix
  matrix <- matrix[, epi_cells] 
  names(break_idx) <- rownames(matrix)
  
  ## Select the most CNA-rich regions
  if(top_method == "abs"){
    abs_vals <- apply(matrix, 1, function(x) mean(abs(x)))
    mat <- matrix[abs_vals > quantile(abs_vals, probs = 1 - top_region), ]
  }
  ## or alternatively - Select the most variant gene regions  
  if(top_method == "sd"){
    stand_dev <- apply(matrix, 1, sd)
    mat <- matrix[stand_dev > quantile(stand_dev, probs = 1 - top_region), ]
  }
  transpose_mat <- as.matrix(t(mat))
  
  # Dimred
  pca <- prcomp(transpose_mat, center = FALSE, scale. = FALSE)
  dim_red_mat <- pca$x[, 1:reduction_dims]
  rownames(dim_red_mat) <- colnames(matrix)
  
  ## Clustering methods:
  if(cluster_method == "louvain"){
    clusters <- cluster_coord(dim_red_mat, method = "louvain", louvain = cluster_k)
  }
  if(cluster_method == "dbscan"){
    if(is.null(dbs_minpts) & is.null(dbs_ptscale)) {minpts <- log(nrow(trans_mat3))}
    if(!is.null(dbs_ptscale)) {minpts <- dbs_ptscale * log(nrow(trans_mat3))}
    if(!is.null(dbs_minpts)) {minpts <- dbs_minpts}
    clusters <- cluster_coord(dim_red_mat, method = "dbscan", min_pts = minpts, diffcut = diffcut)
    
    clusters <- dbscan.cluster(trans_mat3, min_pts = minpts, diffcut = diffcut)
  }
  
  names(clusters) <- colnames(matrix)
  return(clusters)
}



# Arrange Genes by Chromosomes --------------------------------------------


#### genome_sort - with gene list as input, reorder genes by chromosome  position and get breaks for chromosomes and arms
genome_sort <- function(genes, genome = NULL, attributes = c("hgnc_symbol", "start_position", "chromosome_name", "band"), downURL = "https://uswest.ensembl.org"){
  
  # Choose which species to use and server to download from
  if(is.null(genome)){
    mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = downURL, port = 80)
    #Get attributes for gene list
    cna_attr <- biomaRt::getBM(attributes = attributes, filters = "hgnc_symbol", values = genes, mart = mart, uniqueRows = TRUE)
  }
  if(!is.null(genome)){
    cna_attr <- as.data.frame(genome)
    rownames(cna_attr) <- cna_attr[, 1]
    cna_attr <- cna_attr[rownames(cna_attr) %in% genes, ]
    colnames(cna_attr)[1] <- "hgnc_symbol"
  }
  
  # Remove genes not mapped to proper chromosomes
  chr_names <- c(seq(1:22), "X", "Y")
  genes2rm <- setdiff(unique(cna_attr$chromosome_name), chr_names)
  rmgenes <- lapply(genes2rm, function(x){
    gene_drop <- cna_attr[grep(x, cna_attr$chromosome_name), ]
  })
  gene2rm <- do.call(rbind, rmgenes)
  cna_attr <- cna_attr[!is.element(rownames(cna_attr), rownames(gene2rm)), ]
  # Remove doubles
  cna_attr <- cna_attr[!duplicated(cna_attr$hgnc_symbol), ]
  # Remove NA's
  cna_attr <- na.omit(cna_attr)
  
  # Change X and Y chromosomes to numbers for ordering
  cna_attr$chromosome_name <- gsub("X", 23, cna_attr$chromosome_name)
  cna_attr$chromosome_name <- gsub("Y", 24, cna_attr$chromosome_name)
  cna_attr$chromosome_name <- as.numeric(cna_attr$chromosome_name)
  # Order
  cna_attr <- cna_attr[order(cna_attr$chromosome_name, cna_attr$start_position, decreasing = FALSE), ]
  # Chromosome number as sorting vector
  chr_names <- as.character(unique(cna_attr$chromosome_name))
  # Chromosome length
  chr_length <- lapply(chr_names, function(x){
    l <- nrow(subset(cna_attr, cna_attr$chromosome_name == x))
  })
  chr_length <- do.call(rbind, chr_length)
  # Break vector at each chromosome end
  chr_breaks <- cumsum(chr_length)
  # Set all genes as p or q
  cna_attr$band <- gsub("[[:digit:]]", "", cna_attr$band)
  cna_attr$band <- gsub("p.", "p", cna_attr$band)
  cna_attr$band <- gsub("q.", "q", cna_attr$band)
  # Chromosome arm lenghts
  arm_breaks <- integer(length = length(chr_breaks))
  for(i in seq_along(chr_breaks)){
    starting <- chr_breaks[i]
    ending <- nrow(cna_attr[cna_attr$chromosome_name == i & cna_attr$band == "q", ])
    gap_length <- starting - ending
    arm_breaks[i] <- round(gap_length, digits = 0)
  }
  
  # Breaks and labels for arms and full chromosomes
  full_breaks <- sort(c(1, chr_breaks, arm_breaks))
  # Add p and q labels
  q2q <- paste(seq(1, length(chr_breaks)), "q", sep = "")
  p2p <- paste(seq(1, length(chr_breaks)), "p", sep = "")
  full_labels <- sort(c(seq(1, length(chr_breaks)), seq(1, length(chr_breaks))))
  full_labels[seq(2, length(full_labels), 2)] <- q2q
  full_labels[seq(1, length(full_labels), 2)] <- p2p
  # Empty label at end of Y chromosome
  full_labels <- c(full_labels, " ")
  # Name X and Y chromosomes
  full_labels[45:48] <- c("Xp", "Xq", "Yp", "Yq")
  chr_breaks <- c(1, chr_breaks)
  chr_names[23:25] <- c("X", "Y", " ")
  
  output <- list(cna_genes = cna_attr$hgnc_symbol,
                 chr_breaks = chr_breaks,
                 chr_names = chr_names,
                 arm_breaks = arm_breaks,
                 full_breaks = full_breaks,
                 full_labels = full_labels,
                 all = cna_attr)
  return(output)
}


#### Index and label genes by chromosomal location - chromosome number and arm

genome_break <- function(genes, genecut = 10, separate = "chr", genome = NULL){
  
  chr_sort <- genome_sort(genes, genome = genome)
  genes <- chr_sort$cna_genes
  
  if(separate == "chr") {breaks = chr_sort$chr_breaks
  labels = chr_sort$chr_names}
  if(separate == "arm") {breaks = chr_sort$full_breaks
  labels = chr_sort$full_labels}
  if(length(breaks) != length(labels)) {labels = labels[1:length(breaks)]}
  
  # Breaks as factor
  break_idx <- character(length = length(genes))
  for(i in 1:(length(breaks) - 1)){
    if(i == (length(breaks) - 1)) {end = length(genes)}
    if(i != (length(breaks) - 1)) {end = breaks[i + 1]}
    chr_leng <- end - breaks[i]
    if(chr_leng < genecut & chr_leng > 0 & i != 1) {break_idx[breaks[i] + 1:end] = NA}
    if(chr_leng < genecut & chr_leng > 0 & i == 1) {break_idx[breaks[i]:end] = NA}
    if(chr_leng == 0 & i != 1) {break_idx[breaks[i] + 1] = NA}
    if(chr_leng == 0 & i == 1) {break_idx[breaks[i]] = NA}
    if(chr_leng > genecut & i != 1) {break_idx[breaks[i] + 1:end] = i}
    if(chr_leng > genecut & i == 1) {break_idx[breaks[i]:end] = i}
  }
  
  # Remove genes on chromosomes with too few genes
  dispose <- which(is.na(break_idx))
  break_idx <- na.omit(break_idx)
  break_idx <- sort(as.numeric(break_idx))
  genes <- chr_sort$cna_genes[-dispose]
  names(break_idx) <- genes
  
  # Keep only labels with chromosomes that have enough genes
  breaks_df <- cbind.data.frame(breaks, labels)
  rownames(breaks_df) <- seq(1, nrow(breaks_df))
  breaks_df <- breaks_df[unique(break_idx), ]
  
  out <- list(break_idx = break_idx, breaks_df = breaks_df, dispose = dispose)
  return(out)
}




# Utility functions -------------------------------------------------------


#### log and un_log functions:
log_norm <- function(x, centre = FALSE){
  logged <- log2(1 + (x/10))
  if(isTRUE(centre)){
    avg <- rowMeans(logged)
    logged <- sweep(logged, 1, avg)
  }
  return(logged)
}

un_log <- function(x){
  tpmed <- 10 * (2^x - 1)
  return(tpmed)
}


#### Clustering function and methods
## General cluster coordinate function
cluster_coord <- function(matrix, method = c("dbscan", "louvain"), louvain = 25, ...){
  if(method == "dbscan"){
    clusters <- cluster_dbscan(matrix, ...)
  }
  if(method == "louvain"){
    clusters <- cluster_louvain(matrix, k = louvain)
  }
  
  rename_clusts <- function(x){paste("cluster", x, sep = "")}
  clusters <- rename_clusts(clusters)
  names(clusters) <- rownames(matrix)
  return(clusters)
}

## Louvain clustering
cluster_louvain <- function(data, k) {
  knn <- FNN::get.knn(as.matrix(data), k = k)
  knn <- data.frame(from = rep(1:nrow(knn$nn.index), k), to = as.vector(knn$nn.index), weight = 1 / (1 + as.vector(knn$nn.dist)))
  simple_igraph <- igraph::simplify(igraph::graph_from_data_frame(knn, directed = FALSE))
  louvain_clusts <- igraph::cluster_louvain(simple_igraph)
  clusters <- igraph::membership(louvain_clusts)
  return (clusters)
}

## DBScan clustering
cluster_dbscan <- function(matrix, min_pts = log(nrow(matrix)), probs = 3/4, diff_cut = 10){
  dist <- dbscan::kNNdist(matrix, min_pts)
  dist_cutoff <- quantile(sort(dist), probs = probs)
  dist <- dist[dist > dist_cutoff]
  dist <- dist[order(dist)]
  dist <- dist / max(dist)
  der_dist <- diff(dist) / (1 / length(dist))
  knee <- dist[length(der_dist) - length(der_dist[der_dist > diff_cut])] + dist_cutoff
  db_run <- dbscan::dbscan(matrix, eps = knee, minPts = min_pts)
  clusters <- as.factor(db_run$cluster)
  return(clusters)
}


#Dor Simkin
#11:18 AM (0 minutes ago)
#to me

#### Calculate CNAsignal & CNAcorrelation scores

cna_sig_cor <- function(cna_matrix,
                        
                        epi_cells,
                        
                        ref_cells,
                        
                        top_region = 1/3,
                        
                        top_cells = 1/3,
                        
                        cna_sig = "abs"){
  
  
  
  if(is.list(ref_cells)){ref_cells <- unlist(ref_cells)}
  ref_cells <- ref_cells[ref_cells %in% colnames(cna_matrix)]
  
  if(!is.null(intersect(ref_cells, epi_cells))){
    
    rmcells <- intersect(ref_cells, epi_cells)
    
    ref_cells <- setdiff(ref_cells, rmcells)
    
  } 
  
  cna_matrix <- cna_matrix[, c(epi_cells, ref_cells)]
  
  all_cells <- colnames(cna_matrix)
  
  epi_matrix <- cna_matrix[, epi_cells]
  
  
  
  # Absolute CNA value per gene for all epithelial cells
  
  absvals <- apply(epi_matrix, 1, function(x) mean(abs(x)))
  
  # Top one-third of genes with highest absolute value (as set by the argument "top_region")
  
  top_abs <- ifelse(absvals > quantile(absvals, probs = 1 - top_region), absvals, 0)
  
  # Define matrices where to calculate cell averages - only the regions with most CNA
  
  score_cna <- cna_matrix[top_abs > 0, ]
  
  # CNA score for each tumor cell across relevant region
  
  if(cna_sig == "abs"){cna_score <- apply(score_cna, 2, function(x) mean(abs(x)))}
  
  if(cna_sig == "sqrt"){cna_score <- apply(score_cna, 2, function(x) mean(sqrt(x^2)))}
  
  if(cna_sig == "square"){cna_score <- apply(score_cna, 2, function(x) mean(x^2))}
  
  # Calculate correlation only with epithelial cells with strongest signal (defined by a cutoff set by the argument "top_cells")
  
  cellcut <- quantile(cna_score, probs = 1 - top_cells)
  
  epi_cna <- epi_matrix[top_abs > 0, cna_score[epi_cells] > cellcut]
  
  # Correlation vector
  
  cor_vec <- cor(score_cna, epi_cna)
  
  cna_cor <- apply(cor_vec, 1, mean)
  
  
  
  out <- list(cna_signal = cna_score,
              
              cna_correlation = cna_cor)
  
  return(out)
  
}

## plotting function


ggbar = function(colVar,
                 legend_title = NULL,
                 obs=NULL,
                 dir=c('h','v'),
                 cols = c('blue','red','orange','green','magenta')) {
  
  L = list(i = obs,col = colVar)
  lens = sapply(L, length, simplify = F)
  lens = lens[lens!=0]
  lens = unlist(lens)
  
  stopifnot(length(lens)>0 & length(unique(lens))==1)
  
  len = unique(unlist(lens))
  if (is.null(obs)) obs = 1:len
  if (is.null(colVar)) colVar = rep('', len)
  
  dir = match.arg(dir)
  d = data.frame(id = obs,
                 colVar=colVar,
                 stringsAsFactors=F)
  
  d = d %>% dplyr::mutate(id = factor(id, levels = unique(id)))
  
  if (dir=='h') {G = ggplot(d, aes(y=1,x=id,fill=colVar))}
  else {G = ggplot(d, aes(y=id,x=1,fill=colVar))}
  
  
  G + geom_tile() +
    scale_fill_manual(values=cols, name = legend_title, guide = guide_legend(override.aes = list(size = 10))) +
    theme_void() +
    theme(legend.position='top',
          plot.margin=margin(0,0,0,0,'cm')) 
}

