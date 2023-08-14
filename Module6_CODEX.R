# CODEX related code

###################################################################
#Load packages, create color palettes and factor levels, load data#
###################################################################

# Load packages
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(RColorBrewer)
library(ggvoronoi)
library(ggraph)
library(scales)
library(ggpubr)
library(tidyverse)

`%ni%` <- Negate(`%in%`)

# factor levels and color palettes for plotting
an_cols_cell_type <-c(
  #cancer
  "OPC" = "#F97B72",
  "MES" = "#F6CF71",
  "MES_hyp"="#F2B701",
  "AC"="#E68310",
  "NPC" = "#7F3C8D",
  "chromatin_reg" = "#3969AC",
  #normal
  "neuron"="#4B4B8F",
  "oligo"="#B95FBB",
  "astrocyte" = "#D4D915",
  #endothelial
  "vascular"="#CF1C90",
  #immune
  "macrophage" = "#80BA5A",
  "T_cell" = "red",
  "B_cell" = "#39FF14"
)

an_cols_cell_subtype <-c(
  #cancer
  "OPC" = "#F97B72",
  "MES" = "#F6CF71",
  "MES_hyp"="#F2B701",
  "AC"="#E68310",
  "NPC" = "#7F3C8D",
  "chromatin_reg" = "#3969AC",
  #normal
  "neuron"="#4B4B8F",
  "oligo"="#B95FBB",
  "ast1" = "#D4D915",
  "ast2" = "#fffdaf",
  #endothelial
  "vasc1" = "#CF1C90", 
  "vasc2" = "pink",
  #immune
  "T_cell_CD4" = "red",
  "T_cell_CD8" = "red",
  "B_cell" = "#39FF14",
  "mac1" = "#80BA5A",
  "mac2" = "#1AA579",
  "mac3" = "#d3f8d3",
  #other
  "artifact" = "grey50",
  "unknown" = "black"
)

an_cols_ivygap <- c(
  "necrosis" = "red",
  "cellular_tumor" = "#80BA5A",
  "MVP" = "#CF1C90",
  "infiltrating" = "#4B4B8F",
  "pseudopalisading" = "#F2B701",
  "unspecified" = "grey50"
)


levels_cell_type <-
  c(
    "B_cell",
    "T_cell",
    "macrophage",
    "vascular",
    "astrocyte",
    "neuron",
    "oligo",
    "chromatin_reg",
    "NPC",
    "AC",
    "MES_hyp",
    "MES",
    "OPC"
  )

levels_samples_structured <- c(
  "ZH1007_INF",
  "ZH1019_T1_A",
  "ZH916_T1_B",
  "ZH811_T1_v6",
  "ZH1007_NEC_B",
  "ZH916_2_B",
  "ZH811_1_B",
  "ZH1041_T1_B"
)

levels_samples_disorganised <-
  c("ZH916_INF", 
    "ZH811_INF_v6", 
    "MGH258_v6", 
    "ZH1019_INF"
    )


#Import log-transformed and normalized matrix and cell_table
nmat <- readRDS("Codex/mat_norm.RDS")
cells <- readRDS("Codex/cell_table.RDS")
marker <- readRDS("Codex/marker.RDS")

##########################################
# Calculate cell densities per cell type##
##########################################

# remove artifact and unknown cells
cells_clean <- cells %>% filter(cell_type %ni% c("unknown", "artifact"))

# create spatial experiment object
spe_clean <- SpatialExperiment::SpatialExperiment(assays = list(logcounts = nmat[, which(cells$cell_type %ni% c("unknown", "artifact"))]),
                                                  sample_id = cells_clean$sample,
                                                  spatialCoords = as.matrix(cells_clean %>% select(c(centroid_x, centroid_y))),
                                                  rowData = marker,
                                                  colData = list(cell_name = cells_clean$cell_name,
                                                                 cell_type = as.factor(cells_clean$cell_type),
                                                                 ivygap = as.factor(cells_clean$ivygap),
                                                                 region = as.factor(cells_clean$region),
                                                                 centroid_x = cells_clean$centroid_x,
                                                                 centroid_y = cells_clean$centroid_y,
                                                                 img_id = cells_clean$sample,
                                                                 all = rep("all", dim(nmat[, which(cells$cell_type %ni% c("unknown", "artifact"))])[2]))
)

# Create spatial graph by centroid expansion by r = 27.5um (equal to d = 55um)
spe_clean <- imcRtools::buildSpatialGraph(spe_clean, 
                                          img_id = "sample_id", 
                                          type = "expansion", 
                                          threshold = 27.5, 
                                          coords = c("centroid_x", "centroid_y"),
                                          name = "expansion_graph")


# Aggreagate neighbors by absolute counts (not fraction) to asses average density per cell type
spe <- imcRtools::aggregateNeighbors(spe_clean, 
                                     colPairName = "expansion_graph", 
                                     aggregate_by = "metadata",
                                     count_by = "cell_type",
                                     name = "nhood_abs_mat", 
                                     proportions = F)

nhood_abs <- SingleCellExperiment::colData(spe)[["nhood_abs_mat"]] %>% 
  as_tibble()

nhood_abs <- nhood_abs %>%
  mutate(count = rowSums(nhood_abs)) %>%
  mutate(cell_type = cells_clean$cell_type) %>%
  mutate(ivygap = cells_clean$ivygap)


# Plot boxplot per cell type
ggplot(nhood_abs, aes(x=fct_reorder(cell_type, count, .desc = F), y=count, fill=cell_type)) +
  stat_boxplot(geom = "errorbar", color = "black") +
  geom_boxplot(color = "black", outlier.shape = NA, notch = T) +
  scale_fill_manual(values = an_cols_cell_type) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position = "none",
        axis.text = element_text(size = 20, color = "black"), 
        axis.title = element_text(size = 25)) +
  labs(x=NULL, y="cell count per spot-radius") +
  ylim(0,30)


#####################################
#Script for co-localisation analysis#
#####################################

# 1. Create function: Test neighborhood

Test_nhood <- function(sample = "ZH1019_INF", 
                       radius = 27.5, 
                       cell_type_column = "cell_type", 
                       fraction_coherance = 0.8,
                       iter = 10,
                       workers = 10) {
  
  # Observed part  
  cells_clean <- cells %>% 
    filter(cell_type %ni% c("unknown", "artifact")) %>% 
    filter(sample %in% c({{sample}}))
  
  spe_clean <- SpatialExperiment::SpatialExperiment(
    sample_id = cells_clean$sample,
    spatialCoords = as.matrix(cells_clean %>% select(c(centroid_x, centroid_y))),
    rowData = marker,
    colData = list(cell_name = cells_clean$cell_name,
                   cell_type = as.factor(cells_clean[[cell_type_column]]),
                   centroid_x = cells_clean$centroid_x,
                   centroid_y = cells_clean$centroid_y,
                   img_id = cells_clean$sample,
                   all = rep("all", dim(cells_clean)[1]))
  )
  
  spe_clean <- imcRtools::buildSpatialGraph(spe_clean, 
                                            img_id = "sample_id", 
                                            type = "expansion", 
                                            threshold = radius, 
                                            coords = c("centroid_x", "centroid_y"),
                                            name = "expansion_graph")
  
  
  spe_clean <- imcRtools::aggregateNeighbors(spe_clean, 
                                             colPairName = "expansion_graph", 
                                             aggregate_by = "metadata", 
                                             count_by = "cell_type",
                                             name = "nhood_mat", 
                                             proportions = F)
  
  
  nhood_mat <- cbind((SingleCellExperiment::colData(spe_clean)[["nhood_mat"]]), 
                     from_cell_type = cells_clean[[cell_type_column]]) %>% as_tibble
  
  
  nhood_mat_num <- nhood_mat %>% select_if(is.numeric)
  
  row_sums <- nhood_mat_num %>% rowSums()
  fraction_coherance_mat <- nhood_mat_num/row_sums
  drop_rows <- apply(fraction_coherance_mat, 1, function(x) any(x > fraction_coherance))
  nhood_norm_coherant <- nhood_mat[!drop_rows,]
  nhood_norm_coherant <- nhood_norm_coherant %>% drop_na()
  
  ## Summarise interactions
  
  nhood_sum <- nhood_norm_coherant %>% 
    group_by(from_cell_type) %>% 
    summarise(across(where(is.numeric), sum))
  
  ## Calc rowSums and subtract self-pairs
  nhood_rowSums <- nhood_sum %>% 
    select_if(is.numeric) %>% 
    rowSums()
  
  ## Pivot_long and extract self-pairs
  self_pair_count <- nhood_sum %>% 
    pivot_longer(!from_cell_type, names_to = "to_cell_type", values_to = "count") %>% 
    filter(from_cell_type == to_cell_type) %>% 
    pull(count)
  
  ## Normalize nhood_sum by difference rowSums and self_pair_count
  
  nhood_norm <- sweep(nhood_sum %>% select_if(is.numeric), 1, (nhood_rowSums-self_pair_count), "/") %>% 
    cbind(from_cell_type = nhood_sum$from_cell_type, .) %>% as_tibble()
  
  nhood_norm_obs <- nhood_norm %>% 
    pivot_longer(!from_cell_type, names_to = "to_cell_type", values_to = paste0("obs_count_" ,{{sample}}))
  
  # Shuffeled part 
  
  Tmp <- function(sample. = sample, 
                  radius. = radius, 
                  cell_type_column. = cell_type_column, 
                  fraction_coherance. = fraction_coherance,
                  iter. = iter,
                  workers. = workers) 
  {
    
    shuffled_counts <- BiocParallel::bplapply(1:iter., 
                                              BPPARAM = BiocParallel::MulticoreParam(workers = workers., progressbar = T),
                                              function(i) {
                                                
                                                cells_clean <- cells %>% 
                                                  filter(cell_type %ni% c("unknown", "artifact")) %>% 
                                                  filter(sample %in% c(sample.))
                                                
                                                cells_clean[[cell_type_column.]] <- sample(cells_clean[[cell_type_column.]])
                                                
                                                spe_clean <- SpatialExperiment::SpatialExperiment(
                                                  sample_id = cells_clean$sample,
                                                  spatialCoords = as.matrix(cells_clean %>% select(c(centroid_x, centroid_y))),
                                                  rowData = marker,
                                                  colData = list(cell_name = cells_clean$cell_name,
                                                                 cell_type = as.factor(cells_clean[[cell_type_column.]]),
                                                                 centroid_x = cells_clean$centroid_x,
                                                                 centroid_y = cells_clean$centroid_y,
                                                                 img_id = cells_clean$sample,
                                                                 all = rep("all", dim(cells_clean)[1]))
                                                )
                                                
                                                spe_clean <- imcRtools::buildSpatialGraph(spe_clean, 
                                                                                          img_id = "sample_id", 
                                                                                          type = "expansion", 
                                                                                          threshold = radius., 
                                                                                          coords = c("centroid_x", "centroid_y"),
                                                                                          name = "expansion_graph")
                                                
                                                
                                                spe_clean <- imcRtools::aggregateNeighbors(spe_clean, 
                                                                                           colPairName = "expansion_graph", 
                                                                                           aggregate_by = "metadata", 
                                                                                           count_by = "cell_type",
                                                                                           name = "nhood_mat", 
                                                                                           proportions = F)
                                                
                                                
                                                nhood_mat <- cbind((SingleCellExperiment::colData(spe_clean)[["nhood_mat"]]), 
                                                                   from_cell_type = cells_clean[[cell_type_column.]]) %>% as_tibble
                                                
                                                
                                                nhood_mat_num <- nhood_mat %>% select_if(is.numeric)
                                                
                                                row_sums <- nhood_mat_num %>% rowSums()
                                                fraction_coherance_mat <- nhood_mat_num/row_sums
                                                drop_rows <- apply(fraction_coherance_mat, 1, function(x) any(x > fraction_coherance.))
                                                nhood_norm_coherant <- nhood_mat[!drop_rows,]
                                                nhood_norm_coherant <- nhood_norm_coherant %>% drop_na()
                                                
                                                ## Summarise interactions
                                                
                                                nhood_sum <- nhood_norm_coherant %>% 
                                                  group_by(from_cell_type) %>% 
                                                  summarise(across(where(is.numeric), sum))
                                                
                                                ## Calc rowSums and subtract self-pairs
                                                nhood_rowSums <- nhood_sum %>% 
                                                  select_if(is.numeric) %>% 
                                                  rowSums()
                                                
                                                ## Pivot_long and extract self-pairs
                                                self_pair_count <- nhood_sum %>% 
                                                  pivot_longer(!from_cell_type, names_to = "to_cell_type", values_to = "count") %>% 
                                                  filter(from_cell_type == to_cell_type) %>% 
                                                  pull(count)
                                                
                                                ## Normalize nhood_sum by difference rowSums and self_pair_count
                                                
                                                nhood_norm <- sweep(nhood_sum %>% select_if(is.numeric), 1, (nhood_rowSums-self_pair_count), "/") %>% 
                                                  cbind(from_cell_type = nhood_sum$from_cell_type, .) %>% as_tibble()
                                                
                                                nhood_norm_shuff <- nhood_norm %>% 
                                                  pivot_longer(!from_cell_type, names_to = "to_cell_type", values_to = paste0("shuff_", i)) #%>% 
                                                #select(paste0("shuff_", i))
                                                
                                                return(nhood_norm_shuff)
                                              }
    )
    
  }
  
  list_nhood_norm_shuff <- Tmp()
  
  nhood_norm_shuff <- left_join(nhood_norm_obs, reduce(list_nhood_norm_shuff, left_join, by = c("from_cell_type", "to_cell_type"))) %>% as_tibble()
  
  nhood_scaled <- nhood_norm_shuff %>% select_if(is.numeric) %>% t %>%  scale(center = T, scale = T) %>% t
  
  nhood_scaled <- cbind(nhood_norm_shuff %>% select(from_cell_type, to_cell_type), nhood_scaled) %>% as_tibble
  
  
  return(nhood_scaled)
  
}


# 2. Run function for 1 visium spot resolution (r=27.5):

list_nhood <- lapply(cells$sample %>% unique, function(x) Test_nhood(sample = x, 
                                                                     iter = 500, 
                                                                     workers = 500, 
                                                                     cell_type_column = "cell_type", 
                                                                     radius = 27.5))

names(list_nhood) <- cells$sample %>% unique

# 3. Summarize and plot the results

Summarise_nhood <- function(sample) {
  
  tmp <- sample
  
  tmp$z_score <- tmp[[3]]
  tmp$perm_min <- apply(tmp %>% select_if(str_detect(colnames(.), pattern = "shuff")), 1, min)
  tmp$perm_max <- apply(tmp %>% select_if(str_detect(colnames(.), pattern = "shuff")), 1, max)
  tmp$perm_mean <- apply(tmp %>% select_if(str_detect(colnames(.), pattern = "shuff")), 1, mean)
  tmp$perm_median <- apply(tmp %>% select_if(str_detect(colnames(.), pattern = "shuff")), 1, median)
  tmp$perm_sd <- apply(tmp %>% select_if(str_detect(colnames(.), pattern = "shuff")), 1, sd)
  
  tmp$count_larger <- rowSums(tmp %>% select_if(str_detect(colnames(.), pattern = "shuff")) >  tmp[[3]])
  tmp$count_smaller <- rowSums(tmp %>% select_if(str_detect(colnames(.), pattern = "shuff")) <  tmp[[3]])
  
  tmp$interaction_type <- if_else(tmp[[3]] >0, true = "attraction", false = "avoidance")
  tmp$p_value <- if_else(tmp[[3]]>0, true = 1/tmp$count_smaller, false = 1/tmp$count_larger)
  tmp$significant <- if_else(tmp$p_value < 0.05, true = TRUE, false = FALSE)
  tmp$pair <- paste0(tmp$from_cell_type, "_", tmp$to_cell_type)
  
  tmp <- tmp %>% select(from_cell_type, 
                        to_cell_type, 
                        z_score, 
                        interaction_type, 
                        p_value, 
                        significant, 
                        perm_mean, 
                        perm_median, 
                        perm_min, 
                        perm_max, 
                        perm_sd, 
                        count_larger, 
                        count_smaller, 
                        pair)
  
  return(tmp)
  
}


Plot_summary <- function(list_nhood = list_nhood, cell_type_column = cell_type) {
  
  tmp <- BiocParallel::bplapply(list_nhood, 
                                FUN = function(x) Summarise_nhood(sample = x), 
                                BPPARAM = BiocParallel::MulticoreParam(workers = 6, progressbar = T))
  
  nhood_summary <- do.call(rbind.data.frame, tmp)
  nhood_summary <- nhood_summary %>% select(from_cell_type, to_cell_type, interaction_type, z_score)
  nhood_summary <- nhood_summary %>% 
    group_by(from_cell_type, to_cell_type, interaction_type) %>% 
    summarise(mean = mean(z_score), n = n()) %>% 
    group_by(from_cell_type, to_cell_type) %>% 
    arrange(-n) %>% 
    slice(1) %>% 
    ungroup() %>%
    complete(from_cell_type = cells %>% filter({{cell_type_column}} %ni% c("artifact", "unknown")) %>% pull(cell_type) %>% unique(), 
             to_cell_type= cells %>% filter({{cell_type_column}} %ni% c("artifact", "unknown")) %>% pull(cell_type) %>% unique(), 
             fill = list(mean = NA, n = 0))
  
  nhood_summary$percent_samples <- (nhood_summary$n/12)
  
  return(nhood_summary)
}


levels_layers <- c("neuron","NPC","oligo", "astrocyte", "OPC", "AC","vascular", "T_cell", "B_cell", "macrophage","MES","MES_hyp","chromatin_reg")

tmp1 <- Plot_summary(list_nhood = list_nhood[levels_samples_structured])

tmp2 <- Plot_summary(list_nhood = list_nhood[levels_samples_disorganised])


ggpubr::ggarrange(
  
  
  ggplot(tmp1 %>% filter(from_cell_type != to_cell_type), aes(factor(to_cell_type, levels = levels_layers), factor(from_cell_type, levels = levels_layers), col=mean, size=percent_samples)) +
    geom_tile(col = "black", size=0, fill = "white") +
    geom_point(shape = 19) + 
    scale_color_gradient2(low = "white", high = "#b2182b", midpoint=0, na.value = "white") +
    scale_size(range=c(9,16),limits = c(0.1,1)) +
    theme_classic() +
    guides(size = guide_legend(override.aes = list(color = "black"))) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1,size=24)) + 
    theme(axis.text.y =element_text(size=24)) +
    labs(x=NULL, y=NULL, size = "significant in \nfraction of samples") +
    ggtitle("CODEX structured - spot 1")
  
  ,
  
  ggplot(tmp2 %>% filter(from_cell_type != to_cell_type), aes(factor(to_cell_type, levels = levels_layers), factor(from_cell_type, levels = levels_layers), col=mean, size=percent_samples)) +
    geom_tile(col = "black", size=0, fill = "white") +
    geom_point(shape = 19) + 
    scale_color_gradient2(low = "white", high = "#b2182b", midpoint=0, na.value = "white") +
    scale_size(range=c(9,16),limits = c(0.1,1)) +
    theme_classic() +
    guides(size = guide_legend(override.aes = list(color = "black"))) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1,size=24)) + 
    theme(axis.text.y =element_text(size=24)) +
    labs(x=NULL, y=NULL, size = "significant in \nfraction of samples") +
    ggtitle("CODEX disorganised - spot 1")
)
