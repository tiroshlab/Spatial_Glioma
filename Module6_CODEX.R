# CODEX related code
# Load libraries, custom functions and levels  -----------------------------------------
## Load libraries
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(RColorBrewer)
library(ggraph)
library(scales)
library(ggpubr)
library(tidyverse)
#library(ggvoronoi)

## Load custom functions and levels
Count_cell_type <- function(cell_type = cell_type_figure) {
  cells %>% count({{cell_type}}) %>% print(n=Inf)
}

`%ni%` <- Negate(`%in%`)

Create_aligned_pseudospots <- function(sample_codex="ZH916_INF", sample_visium="ZH916inf") {
  
  tmp_codex <- cells_aligned %>% 
    filter(sample %in% c({{sample_codex}})) %>% 
    filter(cell_type_tidy %ni% c("T_cell", "B_cell", "excluded"))
  
  # Create aligned pseudospots codex
  spot_distance <- scalef %>% 
    filter(sample == {{sample_visium}}) %>% 
    pull(scalef_hires)
  
  spot_diameter <- spot_distance*0.55
  
  ## Import spot grid
  spots <- visium_positions %>% 
    filter(sample == {{sample_visium}}) 
  
  ## Calculate distances between tmp and spot centers
  distances <- spots %>%
    expand_grid(cell_name = tmp_codex$cell_name) %>%
    left_join(tmp_codex, by = "cell_name") %>%
    mutate(distance = sqrt((spot_x - aligned_x)^2 + (spot_y - aligned_y)^2))
  
  ## Assign tmp_codex to spots based on distances
  assigned_spots <- distances %>%
    filter(distance <= spot_diameter / 2) %>%
    group_by(cell_name) %>%
    arrange(distance) %>%
    slice_head(n = 1)
  
  rm(distances)
  gc()
  
  ## Append spot names to the original tibble
  tmp_final <- tmp_codex %>%
    left_join(assigned_spots %>% dplyr::select(cell_name, spot_id), by = "cell_name") %>%
    mutate(spot_id = ifelse(is.na(spot_id), NA, spot_id))
  
  ## Calculate and plot dominant cell_type per spot
  tmp_dom <- tmp_final %>%
    filter(cell_type_tidy %ni% c("excluded")) %>%
    group_by(spot_id) %>% 
    count(cell_type_tidy) %>% 
    arrange(spot_id, desc(n)) %>% 
    mutate(n_total = sum(n)) %>% 
    mutate(fraction = n/sum(n)) %>% 
    filter(row_number()==1) %>% 
    left_join(spots %>% select(spot_x, spot_y, spot_id), by = c("spot_id")) %>%
    ungroup() %>% 
    distinct(.)
  
  return(tmp_dom)
}

Run_DEPs <- function(sample) {
  nmat_sample <- nmat[, cells %>% 
                        filter(sample %in% c({{sample}})) %>% 
                        filter(cell_type_figure %ni% c("excluded")) %>%
                        pull(cell_name)
  ]
  
  cells_sample <- cells %>% 
    filter(sample %in% c({{sample}})) %>%  
    filter(cell_type_figure %ni% c("excluded"))
  
  DEPs_hm_sample <-
    tibble(dim(nmat_sample)[1]:dim(nmat_sample)[1]) # Create df in which the DEG for each cluster will be inserted
  
  for (i in sort(unique(cells_sample$cell_type_figure))) {
    # i iterates all clusters
    top_DEP_value <-
      (rowMeans(nmat_sample[, which(i == cells_sample$cell_type_figure)]) - rowMeans(nmat_sample[, which(i != cells_sample$cell_type_figure)]))
    DEPs_hm_sample <-
      cbind(DEPs_hm_sample, top_DEP_value)   # select top gene names
  }
  
  DEPs_hm_sample <-
    DEPs_hm_sample[, -1] # Delete redundant first column
  colnames(DEPs_hm_sample) <-
    paste0(sort(unique(cells_sample$cell_type_figure)), "_", {
      {
        sample
      }
    })
  return(DEPs_hm_sample)
}

Run_DEPs_hm <- function(sample, cell_types, marker) {
  nmat <- nmat[{{marker}} ,cells %>% 
                 filter(sample %in% {{sample}}) %>% 
                 filter(cell_type_figure %in% {{cell_types}}) %>% 
                 pull(cell_name) ]
  
  cells <- cells %>%
    filter(sample %in% {{sample}}) %>%
    filter(cell_type_figure %in% {{cell_types}})
  
  
  DEPs_hm <- tibble(dim(nmat)[1]:dim(nmat)[1]) # Create df in which the DEG for each cluster will be inserted
  
  for (i in sort(unique(cells$cell_type_figure))) {    # i iterates all clusters
    top_DEP_value <- (rowMeans(nmat[,which(i == cells$cell_type_figure) ]) - rowMeans(nmat[,which(i != cells$cell_type_figure)]))
    DEPs_hm <- cbind(DEPs_hm, top_DEP_value)   # select top gene names
  }
  
  DEPs_hm <- DEPs_hm[,-1] # Delete redundant first column
  colnames(DEPs_hm) <- rep( x = paste(sort(unique(cells$cell_type_figure))))
  
  return(DEPs_hm)
}

outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

`%ni%` <- Negate(`%in%`)

Plot_tiles_9 <- function(sample = "ZH916T1") {
  
  mp_cons <- tibble(mp_cons = vspots$mp_cons %>% unique()) %>% 
    drop_na()
  
  vis <- spots_visium_aligned %>% 
    filter(spot_id_sample %in% intersect(spots_codex_aligned$spot_id_sample, spots_visium_aligned$spot_id_sample)) %>% 
    filter(sample == {{sample}}) %>% 
    select(spot_id, mp_cons, spot_x, spot_y)
  
  cod <- spots_codex_aligned %>% 
    filter(spot_id_sample %in% intersect(spots_codex_aligned$spot_id_sample, spots_visium_aligned$spot_id_sample)) %>% 
    filter(sample == {{sample}}) %>% 
    mutate(mp_cons = mp_cons) %>%
    select(spot_id, mp_cons, spot_x, spot_y)
  
  # Calculate the range of x and y coordinates
  x_range <- range(vis$spot_x)
  y_range <- range(vis$spot_y)
  
  # Define the size of the tiles (one-third of the range)
  tile_size_x <- (x_range[2] - x_range[1]) / 3
  tile_size_y <- (y_range[2] - y_range[1]) / 3
  
  # Initialize a list to store the 3x3 tiles
  tiles_vis <- vector("list", length = 9)
  
  # Loop through each tile
  for (i in 1:3) {
    for (j in 1:3) {
      # Define the boundaries of the current tile
      x_min <- x_range[1] + (i - 1) * tile_size_x
      x_max <- x_range[1] + i * tile_size_x
      y_min <- y_range[1] + (j - 1) * tile_size_y
      y_max <- y_range[1] + j * tile_size_y
    }
  }
  
  ggpubr::ggarrange(
    ggplot(vis, aes(x=spot_x, y=spot_y, color = mp_cons)) + 
      geom_point(size = 1, alpha = 1, shape = 16) + 
      guides(colour = guide_legend(override.aes = list(size=2, alpha = 1))) +
      scale_colour_manual(values = an_cols_mp) +
      labs(color = "mp_cons", x = "", y = "") +
      theme_classic() + 
      geom_vline(xintercept = seq(x_range[1], x_range[2], length.out = 4)) +
      geom_hline(yintercept = seq(y_range[1], y_range[2], length.out = 4)) +
      scale_y_reverse() +
      guides(color="none") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ,
    ggplot(cod, aes(x=spot_x, y=spot_y, color=mp_cons)) +
      geom_point(size = 1, alpha=1, shape = 16) +
      scale_colour_manual(values = an_cols_mp) +
      guides(colour = guide_legend(override.aes = list(size=2, alpha = 1))) +
      labs(color = "mp_cons", x = "", y = "") +
      geom_vline(xintercept = seq(x_range[1], x_range[2], length.out = 4)) +
      geom_hline(yintercept = seq(y_range[1], y_range[2], length.out = 4)) +
      theme_classic() +
      scale_y_reverse() +
      guides(color="none") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ,
    nrow = 2
  )
}

Get_prop_9 <- function(sample="ZH1019inf", min_spots_per_tile = 50) {
  
  mp_cons <- tibble(mp_cons = vspots$mp_cons %>% unique()) %>% 
    drop_na()
  
  vis <- spots_visium_aligned %>% 
    filter(spot_id_sample %in% intersect(spots_codex_aligned$spot_id_sample, spots_visium_aligned$spot_id_sample)) %>% 
    filter(sample == {{sample}}) %>% 
    select(spot_id, mp_cons, spot_x, spot_y)
  
  cod <- spots_codex_aligned %>% 
    filter(spot_id_sample %in% intersect(spots_codex_aligned$spot_id_sample, spots_visium_aligned$spot_id_sample)) %>% 
    filter(sample == {{sample}}) %>% 
    mutate(mp_cons = mp_cons) %>%
    select(spot_id, mp_cons, spot_x, spot_y)
  
  ## For vis
  # Calculate the range of x and y coordinates
  x_range <- range(vis$spot_x)
  y_range <- range(vis$spot_y)
  
  # Define the size of the tiles (one-third of the range)
  tile_size_x <- (x_range[2] - x_range[1]) / 3
  tile_size_y <- (y_range[2] - y_range[1]) / 3
  
  # Initialize a list to store the 3x3 tiles
  tiles_vis <- vector("list", length = 9)
  
  # Loop through each tile
  for (i in 1:3) {
    for (j in 1:3) {
      # Define the boundaries of the current tile
      x_min <- x_range[1] + (i - 1) * tile_size_x
      x_max <- x_range[1] + i * tile_size_x
      y_min <- y_range[1] + (j - 1) * tile_size_y
      y_max <- y_range[1] + j * tile_size_y
      
      # Extract cells within the current tile
      current_tile <- vis[
        vis$spot_x >= x_min & vis$spot_x <= x_max &
          vis$spot_y >= y_min & vis$spot_y <= y_max, ]
      
      # Store the current tile in the list
      tiles_vis[[(i - 1) * 3 + j]] <- current_tile
    }
  }
  
  names(tiles_vis) <- paste0("vis_", 1:9)
  
  ## For cod
  # Calculate the range of x and y coordinates
  x_range <- range(cod$spot_x)
  y_range <- range(cod$spot_y)
  
  # Define the size of the tiles (one-third of the range)
  tile_size_x <- (x_range[2] - x_range[1]) / 3
  tile_size_y <- (y_range[2] - y_range[1]) / 3
  
  # Initialize a list to store the 3x3 tiles
  tiles_cod <- vector("list", length = 9)
  
  # Loop through each tile
  for (i in 1:3) {
    for (j in 1:3) {
      # Define the boundaries of the current tile
      x_min <- x_range[1] + (i - 1) * tile_size_x
      x_max <- x_range[1] + i * tile_size_x
      y_min <- y_range[1] + (j - 1) * tile_size_y
      y_max <- y_range[1] + j * tile_size_y
      
      # Extract cells within the current tile
      current_tile <- cod[
        cod$spot_x >= x_min & cod$spot_x <= x_max &
          cod$spot_y >= y_min & cod$spot_y <= y_max, ]
      
      # Store the current tile in the list
      tiles_cod[[(i - 1) * 3 + j]] <- current_tile
    }
  }
  
  names(tiles_cod) <- paste0("cod_", 1:9)
  
  list_prop <- lapply(1:9, function(x) {
    x <- left_join(
      tiles_vis[[x]] %>% 
        count(mp_cons) %>% 
        mutate(prop_vis = n/sum(n))
      ,
      tiles_cod[[x]] %>% 
        count(mp_cons) %>% 
        mutate(prop_cod = n/sum(n))
      ,
      by = "mp_cons") %>% 
      right_join(., mp_cons, by = "mp_cons") %>% 
      replace_na(list(prop_cod = 0, prop_vis = 0, n.x = 0, n.y = 0))
    
  })
  
  
  # Only use quadrants that cover > 150 spots
  filtered_list_prop <- lapply(list_prop, function(df) {
    if(sum(df$n.x) > min_spots_per_tile) {
      return(df)
    }
  })
  
  prop_final <- do.call(rbind, filtered_list_prop)
  
  return(prop_final)
}

## Min distance to cell type with distribution

Min_dist_to_nth_cell_type_distribution <- function(cell_type = "Vasc", 
                                                   sample = "ZH916_INF", 
                                                   n = 1,
                                                   cell_type_column = "cell_type_figure") {
  # Create distance matrix
  cells_clean <- cells %>% filter(cell_type_figure %ni% c("excluded"))
  
  dist <- cells_clean %>% 
    filter(sample == {{sample}}) %>% 
    select(c(centroid_x,centroid_y)) %>% 
    dist() %>% 
    as.matrix()
  
  # Find the index of the closest spot of cell_type "endothelial
  ix_endothelial <- cells_clean %>% 
    filter(sample == {{sample}}) %>% 
    pull({{cell_type_column}}) == {{cell_type}}
  
  
  # Calculate min distance for each spot to the closest spot of cell type {{cell_type}}
  min_dist_to_endothelial <- apply(dist, 1, function(x) sort(x[ix_endothelial])[n])
  
  
  
  # Calculate average distance per cell type to closest endothelial cell
  cells_clean <- cells_clean %>% 
    filter(sample == {{sample}}) %>% 
    select(cell_name, {{cell_type_column}}, sample) %>% 
    cbind(min_dist_to_endothelial) %>% 
    group_by(!!sym(cell_type_column)) #%>% 
  # summarise(mean_dist = mean(min_dist_to_endothelial), 
  #           median_dist = median(min_dist_to_endothelial)) %>% 
  # arrange(desc(mean_dist))
  
  rm(dist)
  gc()
  
  return(cells_clean)
}

## Test neighborhood with shuffling

Test_nhood <- function(sample = "ZH1019_INF", 
                       radius = 27.5, 
                       cell_type_figure_column = "cell_type_figure", 
                       fraction_coherance = 0.8,
                       iter = 10,
                       workers = 10) {
  
  # Observed part  
  cells_clean <- cells %>% 
    filter(cell_type_figure %ni% c("excluded")) %>% 
    filter(sample %in% c({{sample}}))
  
  spe_clean <- SpatialExperiment::SpatialExperiment(
    sample_id = cells_clean$sample,
    spatialCoords = as.matrix(cells_clean %>% select(c(centroid_x, centroid_y))),
    rowData = marker,
    colData = list(cell_name = cells_clean$cell_name,
                   cell_type_figure = as.factor(cells_clean[[cell_type_figure_column]]),
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
                                             count_by = "cell_type_figure",
                                             name = "nhood_mat", 
                                             proportions = F)
  
  
  nhood_mat <- cbind((SingleCellExperiment::colData(spe_clean)[["nhood_mat"]]), 
                     from_cell_type_figure = cells_clean[[cell_type_figure_column]]) %>% as_tibble
  
  
  nhood_mat_num <- nhood_mat %>% select_if(is.numeric)
  
  row_sums <- nhood_mat_num %>% rowSums()
  fraction_coherance_mat <- nhood_mat_num/row_sums
  drop_rows <- apply(fraction_coherance_mat, 1, function(x) any(x > fraction_coherance))
  nhood_norm_coherant <- nhood_mat[!drop_rows,]
  nhood_norm_coherant <- nhood_norm_coherant %>% drop_na()
  
  ## Summarise interactions
  
  nhood_sum <- nhood_norm_coherant %>% 
    group_by(from_cell_type_figure) %>% 
    summarise(across(where(is.numeric), sum))
  
  ## Calc rowSums and subtract self-pairs
  nhood_rowSums <- nhood_sum %>% 
    select_if(is.numeric) %>% 
    rowSums()
  
  ## Pivot_long and extract self-pairs
  self_pair_count <- nhood_sum %>% 
    pivot_longer(!from_cell_type_figure, names_to = "to_cell_type_figure", values_to = "count") %>% 
    filter(from_cell_type_figure == to_cell_type_figure) %>% 
    pull(count)
  
  ## Normalize nhood_sum by difference rowSums and self_pair_count
  
  nhood_norm <- sweep(nhood_sum %>% select_if(is.numeric), 1, (nhood_rowSums-self_pair_count), "/") %>% 
    cbind(from_cell_type_figure = nhood_sum$from_cell_type_figure, .) %>% as_tibble()
  
  nhood_norm_obs <- nhood_norm %>% 
    pivot_longer(!from_cell_type_figure, names_to = "to_cell_type_figure", values_to = paste0("obs_count_" ,{{sample}}))
  
  # Shuffeled part 
  
  Tmp <- function(sample. = sample, 
                  radius. = radius, 
                  cell_type_figure_column. = cell_type_figure_column, 
                  fraction_coherance. = fraction_coherance,
                  iter. = iter,
                  workers. = workers) 
  {
    
    shuffled_counts <- BiocParallel::bplapply(1:iter., 
                                              BPPARAM = BiocParallel::MulticoreParam(workers = workers., progressbar = T),
                                              function(i) {
                                                
                                                cells_clean <- cells %>% 
                                                  filter(cell_type_figure %ni% c("excluded")) %>% 
                                                  filter(sample %in% c(sample.))
                                                
                                                cells_clean[[cell_type_figure_column.]] <- sample(cells_clean[[cell_type_figure_column.]])
                                                
                                                spe_clean <- SpatialExperiment::SpatialExperiment(
                                                  sample_id = cells_clean$sample,
                                                  spatialCoords = as.matrix(cells_clean %>% select(c(centroid_x, centroid_y))),
                                                  rowData = marker,
                                                  colData = list(cell_name = cells_clean$cell_name,
                                                                 cell_type_figure = as.factor(cells_clean[[cell_type_figure_column.]]),
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
                                                                                           count_by = "cell_type_figure",
                                                                                           name = "nhood_mat", 
                                                                                           proportions = F)
                                                
                                                
                                                nhood_mat <- cbind((SingleCellExperiment::colData(spe_clean)[["nhood_mat"]]), 
                                                                   from_cell_type_figure = cells_clean[[cell_type_figure_column.]]) %>% as_tibble
                                                
                                                
                                                nhood_mat_num <- nhood_mat %>% select_if(is.numeric)
                                                
                                                row_sums <- nhood_mat_num %>% rowSums()
                                                fraction_coherance_mat <- nhood_mat_num/row_sums
                                                drop_rows <- apply(fraction_coherance_mat, 1, function(x) any(x > fraction_coherance.))
                                                nhood_norm_coherant <- nhood_mat[!drop_rows,]
                                                nhood_norm_coherant <- nhood_norm_coherant %>% drop_na()
                                                
                                                ## Summarise interactions
                                                
                                                nhood_sum <- nhood_norm_coherant %>% 
                                                  group_by(from_cell_type_figure) %>% 
                                                  summarise(across(where(is.numeric), sum))
                                                
                                                ## Calc rowSums and subtract self-pairs
                                                nhood_rowSums <- nhood_sum %>% 
                                                  select_if(is.numeric) %>% 
                                                  rowSums()
                                                
                                                ## Pivot_long and extract self-pairs
                                                self_pair_count <- nhood_sum %>% 
                                                  pivot_longer(!from_cell_type_figure, names_to = "to_cell_type_figure", values_to = "count") %>% 
                                                  filter(from_cell_type_figure == to_cell_type_figure) %>% 
                                                  pull(count)
                                                
                                                ## Normalize nhood_sum by difference rowSums and self_pair_count
                                                
                                                nhood_norm <- sweep(nhood_sum %>% select_if(is.numeric), 1, (nhood_rowSums-self_pair_count), "/") %>% 
                                                  cbind(from_cell_type_figure = nhood_sum$from_cell_type_figure, .) %>% as_tibble()
                                                
                                                nhood_norm_shuff <- nhood_norm %>% 
                                                  pivot_longer(!from_cell_type_figure, names_to = "to_cell_type_figure", values_to = paste0("shuff_", i)) #%>% 
                                                #select(paste0("shuff_", i))
                                                
                                                return(nhood_norm_shuff)
                                              }
    )
    
  }
  
  list_nhood_norm_shuff <- Tmp()
  
  nhood_norm_shuff <-
    left_join(nhood_norm_obs,
              reduce(
                list_nhood_norm_shuff,
                left_join,
                by = c("from_cell_type_figure", "to_cell_type_figure")
              )) %>% as_tibble()
  
  nhood_scaled <-
    nhood_norm_shuff %>% select_if(is.numeric) %>% t %>%  scale(center = T, scale = T) %>% t
  
  nhood_scaled <-
    cbind(
      nhood_norm_shuff %>% select(from_cell_type_figure, to_cell_type_figure),
      nhood_scaled
    ) %>% as_tibble
  
  
  return(nhood_scaled)
  
}

## Summarize nhood 

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
  tmp$pair <- paste0(tmp$from_cell_type_figure, "_", tmp$to_cell_type_figure)
  
  tmp <- tmp %>% select(from_cell_type_figure, 
                        to_cell_type_figure, 
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


## Plot summary

Plot_summary <- function(list_nhood = list_nhood, cell_type_figure_column) {
  
  tmp <- BiocParallel::bplapply(list_nhood, 
                                FUN = function(x) Summarise_nhood(sample = x), 
                                BPPARAM = BiocParallel::MulticoreParam(workers = 6, progressbar = T))
  
  nhood_summary <- do.call(rbind.data.frame, tmp)
  nhood_summary <- nhood_summary %>% select(from_cell_type_figure, to_cell_type_figure, interaction_type, z_score)
  nhood_summary <- nhood_summary %>% 
    group_by(from_cell_type_figure, to_cell_type_figure, interaction_type) %>% 
    summarise(mean = mean(z_score), n = n()) %>% 
    group_by(from_cell_type_figure, to_cell_type_figure) %>% 
    arrange(-n) %>% 
    slice(1) %>% 
    ungroup() %>%
    complete(from_cell_type_figure = cells %>% filter({{cell_type_figure_column}} %ni% c("excluded")) %>% pull({{cell_type_figure_column}}) %>% unique(), 
             to_cell_type_figure= cells %>% filter({{cell_type_figure_column}} %ni% c("excluded")) %>% pull({{cell_type_figure_column}}) %>% unique(), 
             fill = list(mean = NA, n = 0))
  
  nhood_summary$percent_samples <- (nhood_summary$n/12)
  
  return(nhood_summary)
}


levels_samples_structured <- c(
  "ZH916_T1_B", 
  "ZH1007_NEC_B",
  "ZH916_2_B",
  "ZH1019_T1_A",
  "ZH1041_T1_B",
  "ZH881_T1_v6",
  "ZH881_INF_v6",
  "ZH881_1_B",
  "ZH1007_INF"
)

levels_samples_disorganised <- c(
  "MGH258_v6",
  "ZH1019_INF",
  "ZH916_INF")

levels_samples <-
  c(
    "ZH916_INF",
    "ZH881_INF_v6",
    "ZH1007_INF",
    "ZH1019_INF",
    "ZH916_T1_B",
    "ZH881_T1_v6",
    "ZH1007_NEC_B",
    "ZH1019_T1_A",
    "ZH916_2_B",
    "ZH881_1_B",
    "ZH1041_T1_B",
    "MGH258_v6"
  )

key_sample <-
  tibble(
    codex = c(
      "ZH916_INF",
      "ZH881_INF_v6",
      "ZH1007_INF",
      "ZH1019_INF",
      "ZH916_T1_B",
      "ZH881_T1_v6",
      "ZH1007_NEC_B",
      "ZH1019_T1_A",
      "ZH881_1_B",
      "MGH258_v6"
    ),
    visium = c(
      "ZH916inf",
      "ZH881inf",
      "ZH1007inf",
      "ZH1019inf",
      "ZH916T1",
      "ZH881T1",
      "ZH1007nec",
      "ZH1019T1",
      "ZH8811Bbulk",
      "MGH258"
    )
  )

# Import data -------------------------------------------------------------
##Import log-transformed and normalized matrix and cell_table
nmat <- readRDS("Codex/mat_norm.RDS")
cells <- readRDS("Codex/cells_table.RDS")
marker <- readRDS("Codex/marker.RDS")

vspots <- readRDS("Codex/visium_spots.RDS")
visium_positions <- readRDS("Codex/visium_spot_positions_all.RDS")
spots_visium_aligned <- readRDS("Codex/spots_aligned_visium.RDS")
spots_codex_aligned <- readRDS("Codex/pseudospots_aligned_codex.RDS")
cells_aligned <- readRDS("Codex/cells_aligend_codex.RDS")
scalef <- readRDS("Codex/scale_factors_hires.RDS")

spot_dom <- readRDS("Codex/CODEX_pseudospot_dominant.RDS")
spot_comp <- readRDS("Codex/CODEX_pseudospot_composition.RDS")


# Calculate cell densities per cell type -----------------------------------
## Create SpatialExperiment (spe)
cells_clean <- cells %>% filter(cell_type_figure %ni% c("excluded"))

spe <- SpatialExperiment::SpatialExperiment(assays = list(logcounts = nmat[, which(cells$cell_type_figure %ni% c("excluded"))]),
                                            sample_id = cells_clean$sample,
                                            spatialCoords = as.matrix(cells_clean %>% select(c(centroid_x, centroid_y))),
                                            rowData = marker,
                                            colData = list(cell_name = cells_clean$cell_name,
                                                           cell_type_figure = as.factor(cells_clean$cell_type_figure),
                                                           ivygap = as.factor(cells_clean$ivygap),
                                                           region = as.factor(cells_clean$region),
                                                           centroid_x = cells_clean$centroid_x,
                                                           centroid_y = cells_clean$centroid_y,
                                                           img_id = cells_clean$sample,
                                                           all = rep("all", dim(nmat[, which(cells$cell_type_figure %ni% c("excluded"))])[2]))
)

## Create spatial graph by centroid expansion and add to spe (colPair)
spe <- imcRtools::buildSpatialGraph(spe, 
                                    img_id = "sample_id", 
                                    type = "expansion", 
                                    threshold = 27.5, 
                                    coords = c("centroid_x", "centroid_y"),
                                    name = "expansion_graph")

SingleCellExperiment::colPairNames(spe) #shows constructed graph
colPair(spe, "expansion_graph") # prints summary of constructed graph
colPair(spe, "expansion_graph") %>% class #

# The graph is stored in form of a SelfHits object in colPair(object, name). This object can be regarded as an edgelist and coerced to an igraph object via:

igraph::graph_from_edgelist(as.matrix(colPair(spe, "expansion_graph")))

## Aggregate neighbors and create neighborhood matrix
spe <- imcRtools::aggregateNeighbors(spe, 
                                     colPairName = "expansion_graph", 
                                     aggregate_by = "metadata", 
                                     count_by = "cell_type_figure",
                                     name = "nhood_mat")

as_tibble(SingleCellExperiment::colData(spe)[["nhood_mat"]]) %>% dim
as_tibble(SingleCellExperiment::colData(spe)[["nhood_mat"]])

nhood <- as_tibble(SingleCellExperiment::colData(spe)[["nhood_mat"]])

## Aggregate neighbors for cell density plots
spe <- imcRtools::aggregateNeighbors(spe, 
                                     colPairName = "expansion_graph", 
                                     aggregate_by = "metadata",
                                     count_by = "cell_type_figure",
                                     name = "nhood_abs_mat", 
                                     proportions = F)

nhood_abs <- SingleCellExperiment::colData(spe)[["nhood_abs_mat"]] %>% as_tibble()
nhood_abs <- nhood_abs %>% 
  mutate(count = rowSums(nhood_abs)) %>% 
  mutate(cell_type_figure = cells_clean$cell_type_figure) %>% 
  mutate(ivygap = cells_clean$ivygap)


## Add malignant column for later plots
nhood_abs <- nhood_abs %>% 
  mutate(malignant = if_else(cell_type_figure %in% c("Chromatin-Reg","MES-Hyp", "MES", "AC", "OPC", "NPC"), true = "malignant", false = cell_type_figure)) %>% 
  mutate(malignant = if_else(cell_type_figure %in% c("Mac", "Inflammatory-Mac", "T-cell", "B-cell", "Vasc"), true = "immune", false = malignant)) %>% 
  mutate(malignant = if_else(cell_type_figure %in% c("Reactive-Ast", "Oligo", "Neuron"), true = "normal", false = malignant))

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


# Co-localisation analysis ------------------------------------------------
## 1. See custom function: Test neighborhood
## 2. Run function for 1 visium spot resolution (r=27.5) for all samples:

list_nhood <- lapply(cells$sample %>% unique(), function(x) Test_nhood(sample = x, 
                                                                       iter = 500, 
                                                                       workers = 100, 
                                                                       cell_type_figure_column = "cell_type_tidy", 
                                                                       radius = 27.5))

names(list_nhood) <- cells$sample %>% unique()

## 3. Summarize and plot the results

levels_layers <-
  c(
    "Neuron",
    "NPC",
    "Oligo",
    "Reactive_Ast",
    "OPC",
    "AC",
    "Mac",
    "Vasc",
    "T_cell",
    "B_cell",
    "MES",
    "Inflammatory_Mac",
    "MES_Hyp",
    "Chromatin_Reg"
  )


tmp1 <- Plot_summary(list_nhood = list_nhood[levels_samples_structured], cell_type_figure_column = "cell_type_tidy")

tmp2 <- Plot_summary(list_nhood = list_nhood[levels_samples_disorganised], cell_type_figure_column = "cell_type_tidy")


#plot 10*22
ggpubr::ggarrange(
  
  
  ggplot(tmp1 %>% 
           filter(from_cell_type_figure != to_cell_type_figure) %>% 
           drop_na(),
         aes(factor(to_cell_type_figure, levels = levels_layers), factor(from_cell_type_figure, levels = levels_layers), col=mean, size=percent_samples)) +
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
  
  ggplot(tmp2 %>% 
           filter(from_cell_type_figure != to_cell_type_figure) %>% 
           drop_na(),
         aes(factor(to_cell_type_figure, levels = levels_layers), factor(from_cell_type_figure, levels = levels_layers), col=mean, size=percent_samples)) +
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



## Structured only
ggplot(tmp1 %>% 
         filter(from_cell_type_figure != to_cell_type_figure) %>% 
         drop_na(),
       aes(factor(to_cell_type_figure, levels = levels_layers), factor(from_cell_type_figure, levels = levels_layers), col=mean, size=percent_samples)) +
  geom_tile(col = "white", size=0, fill = "white") +
  geom_point(shape = 19) + 
  scale_color_gradient2(
    low = "white",
    high = "#b2182b",
    midpoint = 0,
    na.value = "white",
    limit = c(0, 20),
    oob = scales::squish, 
    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  scale_size(range=c(5,12),limits = c(0.1,1)) +
  theme_classic() +
  theme(axis.text.x = element_text(
    angle = 40,
    vjust = 1,
    hjust = 1,
    size = 20)) + 
  theme(axis.text.y = element_text(size = 20)) + 
  theme(legend.text = element_text(size = 20)) + 
  theme(legend.title = element_text(size = 20)) +
  labs(x=NULL, y=NULL, size = "significant in \nfraction of samples")

## 4. Calculate mean of each reciprocal pair for network graph

pair_struct_mean <- tmp1 %>% 
  rowwise() %>% 
  mutate(pair = paste0(sort(c(from_cell_type_figure, to_cell_type_figure)), collapse = ".")) %>% 
  arrange(pair) %>% 
  filter(from_cell_type_figure != to_cell_type_figure) %>% 
  filter(from_cell_type_figure %ni% c("excluded")) %>% 
  filter(to_cell_type_figure %ni% c("excluded")) 

pair_struct_mean <- pair_struct_mean %>% 
  group_by(pair) %>% 
  summarise(pair_mean = mean(mean),
            pair_n = mean(n),
            pair_percent_samples = mean(percent_samples))

pair_struct_mean <- pair_struct_mean %>% separate(col = pair, into = c("from", "to"), sep = "\\.", remove = F)

write_tsv(pair_struct_mean %>% select(-pair), file = "pair_struct_mean.tsv")


## 5. Calculate max of each pair for network graph (To not loose assymetric pairs like MES-Hyp + MES)
pair_struct_max <- tmp1 %>% 
  rowwise() %>% 
  mutate(pair = paste0(sort(c(from_cell_type_figure, to_cell_type_figure)), collapse = ".")) %>% 
  arrange(pair) %>% 
  filter(from_cell_type_figure != to_cell_type_figure) %>% 
  filter(from_cell_type_figure %ni% c("excluded")) %>% 
  filter(to_cell_type_figure %ni% c("excluded")) 

pair_struct_max <- pair_struct_max %>% 
  arrange(pair, desc(mean)) %>% 
  group_by(pair) %>% 
  filter(row_number() == 1) %>% 
  ungroup() %>% 
  arrange(desc(mean))


# Create aligned pseudospots ----------------------------------------------
## Create list of df with per sample CODEX pseudospots
list_spots_aligned <- map2(.x = key_sample$codex,
                           .y = key_sample$visium,
                           .f = ~Create_aligned_pseudospots(sample_codex = .x, sample_visium = .y))

names(list_spots_aligned) <- key_sample$visium

## Flatten list do df
list_spots_aligned <- map(names(list_spots_aligned), 
                          function(x) list_spots_aligned[[x]] %>% mutate(sample = x)
)

spots_codex_aligned <- do.call(rbind, list_spots_aligned) %>% 
  mutate(spot_id_sample = paste0(sample, "_", spot_id)) %>% 
  rename(mp_cons = cell_type_tidy) %>% 
  drop_na()

## Filter aligned visium spots
spots_visium_aligned <- vspots %>% 
  filter(sample %in% key_sample$visium) %>% 
  mutate(spot_id_sample = paste0(sample, "_", spot_id))
