library(patchwork)
library(parallel)

#!!! Run Functions section first

# load data ---------------------------------------------------------------

sample_ls <- (read.delim("general/GBM_samples.txt", header = FALSE))$V1

gen_clusters <- as.character(unique(unlist(sapply(c(1:length(sample_ls)), function(i){
  mp_assign <- readRDS(paste("MP/mp_assign_124/", sample_ls[i], ".rds", sep = ""))
  return(unique(mp_assign$spot_type_meta_new))
}))))



# calculate spatial coherence (Figures 4A-C) --------

rand_num <- 100

sapply(c(1:length(sample_ls)), function(i){
  
  print(sample_ls[i])
  
  # load data
  spots_positions <- read.csv(paste("general/GBM_data/", sample_ls[i] , "/outs/spatial/tissue_positions_list.csv", sep = ""), header = FALSE, stringsAsFactors = FALSE)
  row.names(spots_positions) <- spots_positions$V1
  
  spots_clusters <- readRDS(paste("MP/mp_assign_124/", sample_ls[i], ".rds", sep = ""))
  spots_clusters <- na.omit(spots_clusters)
  colnames(spots_clusters) <- c("barcodes", "spot_type")
  row.names(spots_clusters)<- spots_clusters$barcodes  
  
  # abundance 
  
  programs_comp <- sample_programs_composition(spots_clusters,gen_clusters)
  
  # neighbors tables   
  
  neighbors_table <- neighbors_table_func(spots_positions,spots_clusters)
  
  rand_neighbors_table <- lapply(c(1:rand_num), function(i){
    new_pos <- sample(spots_positions$V1, length(spots_positions$V1), replace = FALSE)
    pos_table <- spots_positions
    pos_table$V1 <- new_pos
    pos_table$V2 <- spots_positions[new_pos, "V2"]
    
    neighbors_table <- neighbors_table_func(pos_table,spots_clusters)
    return(neighbors_table)
  })
  
  
  # spatial coherence
  
  programs_spatial_score <- sapply(sort(gen_clusters), function(cluster){
    if (!(cluster %in% spots_clusters$spot_type)) {
      prog_score <- NaN
    } else {
      program_neighbors_table = neighbors_table[row.names(neighbors_table) %in% spots_clusters$barcodes[spots_clusters$spot_type == cluster],]
      obs <- obs_program_spatial_score(program_neighbors_table, cluster)
      one_spatial_score <- one_val(dim(program_neighbors_table)[1])
      zero_spatial_score <- zero_val(rand_neighbors_table, spots_clusters, cluster)
      if (obs>one_spatial_score){obs <- one_spatial_score}
      if (obs<zero_spatial_score){obs <- zero_spatial_score}
      
      prog_score <- (obs - zero_spatial_score)/(one_spatial_score - zero_spatial_score)
    }
    return(prog_score)
  })
  
  saveRDS(programs_spatial_score, "save_path")
})

# plot 

sapply(c(1:length(sample_ls)), function(i){
  
  print(sample_ls[i])
  
  # load data
  spots_clusters <- readRDS(paste("MP/mp_assign/mp_assign_124/", sample_ls[i], ".rds", sep = ""))
  spots_clusters <- na.omit(spots_clusters)
  colnames(spots_clusters) <- c("barcodes", "spot_type")
  row.names(spots_clusters)<- spots_clusters$barcodes  
  
  # abundance 
  
  programs_comp <- sample_programs_composition(spots_clusters,gen_clusters)
  saveRDS(programs_comp, paste("save_path_abund", sample_ls[i],"_nodes_abund.rds", sep = ""))
})

all_comp <- list.files("save_path_abund")

compositions <- sapply(all_comp, function(comp){
  samp_comp <- readRDS(paste("save_path_abund", comp, sep = ""))
  return(samp_comp)
})

colnames(compositions) <- as.character(sapply(colnames(compositions), function(x){substr(x, 1, nchar(x)-16)}))

all_spatial_score <- list.files("save_path")

spatial_score <- sapply(all_spatial_score, function(s_score){
  samp_scores <- readRDS(paste("save_path", s_score, sep = ""))
  return(samp_scores)
})


spatial_score <- as.data.frame(spatial_score)
colnames(spatial_score) <- as.character(sapply(colnames(spatial_score), function(x){substr(x,1, nchar(x)-18)}))

new_score <- ifelse(compositions < 0.05, 0,1) * spatial_score # filtered <0.05
new_score[new_score == 0] <- NA

new_score$metaprogram <- row.names(spatial_score)
spatial_score_filt <- new_score



mal_mp <- c("AC", "NPC", "OPC", "MES", "chromatin.reg", "AC.MES", "metabolism", "inflammatory.resp", "macrophage", "vascular", "oligo", "reactive.astrocyte")
prog_order <- rev(spatial_score_filt$metaprogram[order(apply(spatial_score_filt[,c(1:length(sample_ls))],1, function(x){mean(na.omit(x[x!=0]))}))])
samples_order <- rev(colnames(spatial_score_filt)[c(1:length(sample_ls))][order(apply(spatial_score_filt[spatial_score_filt$metaprogram %in% c("AC", "NPC", "OPC", "MES", "chromatin.reg", "AC.MES", "inflammatory.resp", "macrophage", "vascular","oligo", "reactive.astrocyte"),c(1:length(sample_ls))],2, function(x){mean(na.omit(x[x!=0]))}))])

spatial_long <- melt(setDT(spatial_score_filt), id.vars = "metaprogram", variable.name = "sample", value.name = "spatial_score")
spatial_long <- spatial_long[spatial_long$spatial_score != 0,]


spatial_long$sample <- factor(spatial_long$sample, levels = samples_order)
spatial_long$metaprogram <- factor(spatial_long$metaprogram, levels = prog_order)

mp_spatial_mean <- as.data.frame(apply(new_score[,c(1:26)],1, function(x){mean(na.omit(x))}))
colnames(mp_spatial_mean) <- c("mean_spatial")
mp_spatial_mean$metaprogram <- new_score$metaprogram
mp_spatial_mean$sd <- apply(new_score[,c(1:26)],1, function(x){sd(na.omit(x))})
mp_ord_tmp <- apply(new_score[,c(1:26)],1, function(x){mean(na.omit(x))})
names(mp_ord_tmp) <- new_score$metaprogram
mp_ord <- names(sort(mp_ord_tmp, decreasing = T))

mp_spatial_mean$metaprogram <- factor(mp_spatial_mean$metaprogram, levels = mp_ord)

samp_pal<-met.brewer("Veronese",26)
meta_pal<-c("hypoxia"="#F2B701", "MES"="#F6CF71", "AC.MES"="#66C5CC","AC"="#E68310","OPC"="#f97b72", "NPC"="#7F3C8D","chromatin.reg"="#3969AC", "metabolism"="#EF4868", "neuron"="#4b4b8f","oligo"="#B95FBB", "reactive.astrocyte"="#D4D915", "macrophage"="#80BA5A", "inflammatory.resp"="#11A579","vascular"="#CF1C90")

# plot was not not used in paper
p1_mp <- ggplot(spatial_long, aes(x=factor(metaprogram), y=spatial_score, fill = sample)) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(size= 12, angle = 45, vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=14, face="bold"),
        plot.title = element_text(size=14, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) + 
  ylab("spatial cohernece score") +
  xlab("metaprogram") +
  ggtitle("Spatial Coherence by Cancer Metaprogram") +
  scale_fill_manual(values = samp_pal)

# fig 4C
p2_mp <- ggplot(data=mp_spatial_mean, aes(x=metaprogram, y=mean_spatial, group=1)) + 
  geom_line(color = "grey", size = 2)+
  geom_point() +
  geom_errorbar(aes(ymin=mean_spatial-sd, ymax=mean_spatial+sd), width=.2,
                position=position_dodge(0.05)) + theme_minimal() +
  theme(axis.text.x = element_text(size= 12, angle = 45, vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=14, face="bold"),
        plot.title = element_text(size=14, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) + 
  ylab("mean spatial cohernece score") 

#fig 4B
p1_samp <- ggplot(spatial_long[spatial_long$metaprogram %in% c("AC", "NPC", "OPC", "MES", "chromatin.reg", "AC.MES", "inflammatory.resp", "macrophage", "vascular","oligo", "reactive.astrocyte")], aes(x=factor(sample), y=spatial_score, fill = metaprogram)) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(size= 12, angle = 45, vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=14, face="bold"),
        plot.title = element_text(size=14, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) + 
  ylab("spatial cohernece score") +
  xlab("sample") +
  ggtitle("Spatial Coherence by samples") +
  scale_fill_manual(values = meta_pal)

#[ plot was not used in paper]
spatial_mean <- as.data.frame(apply(new_score[c(1:3,5:7,10,12,14),c(1:26)],2, function(x){mean(na.omit(x))})) # claculate mean spatial coherence excluding hypoxia and metabolism and normal brain (oligo, neuron and reactive.ast) 
colnames(spatial_mean) <- c("mean_spatial")
spatial_mean$sample <- row.names(spatial_mean)
spatial_mean$sd <- apply(new_score[c(1:3,5:7,10,12,14),c(1:26)],2, function(x){sd(na.omit(x))})

samples_ord <- names(sort(apply(new_score[c(1:3,5:7,10,12,14),c(1:26)],2, function(x){mean(na.omit(x))}), decreasing = T))

spatial_mean$sample <- factor(spatial_mean$sample, levels = samples_ord)

p2_samp <- ggplot(data=spatial_mean[spatial_mean$sample,], aes(x=sample, y=mean_spatial, group=1)) +
  geom_line(color = "grey", size = 2)+
  geom_point() + 
  geom_errorbar(aes(ymin=mean_spatial-sd, ymax=mean_spatial+sd), width=.2,
                position=position_dodge(0.05)) + theme_minimal() +
  theme(axis.text.x = element_text(size= 12, angle = 45, vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=14, face="bold"),
        plot.title = element_text(size=14, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) + 
  ylab("mean spatial cohernece score") 



# calculate spatial coherence by window to define structured vs organized spots (Figures 4D-E) (Alt. skip and use saved results in next section) ----------------------------------------------------------

win_size <- 5 # run for c(5,8,11) 
rand_num <- 100


samples_num <- c(1:length(sample_ls))
all_scores <- mclapply(samples_num, all_scores_fx, mc.cores = 26)

all_scores_fx <- function(i){
  
  #print(sample_ls[i])
  
  # load data
  spots_positions <- read.csv(paste("general/GBM_data/", sample_ls[i] , "/outs/spatial/tissue_positions_list.csv", sep = ""), header = FALSE, stringsAsFactors = FALSE)
  row.names(spots_positions) <- spots_positions$V1
  
  spots_clusters <- readRDS(paste("MP/mp_assign_124/", sample_ls[i], ".rds", sep = ""))
  spots_clusters <- na.omit(spots_clusters)
  colnames(spots_clusters) <- c("barcodes", "spot_type")
  row.names(spots_clusters)<- spots_clusters$barcodes  
  
  
  # neighbors tables   
  all_win_neighbors_table <- neighbors_table_funcV2(spots_positions,spots_clusters)
  neighbors_table <- neighbors_table_func(spots_positions,spots_clusters)
  
  # spatial coherence by window 
  
  all_spots <- spots_clusters$barcodes
  
  spatial_score_win <- sapply(all_spots, function(spot){
    win_spots <- c(spot)
    sapply(c(1:win_size), function(j){
      win_spots <<- unique(c(win_spots,unique(na.omit(as.character(all_win_neighbors_table[win_spots,])))))
      win_spots <<- win_spots[win_spots != "NaN"]
    })
    
    win_rand_neighbors_table <- lapply(c(1:rand_num), function(k){
      win_spots_positions <- spots_positions[spots_positions$V1 %in% win_spots,]
      new_pos <- sample(win_spots_positions$V1, length(win_spots_positions$V1), replace = FALSE)
      pos_table <- win_spots_positions
      row.names(pos_table) <- new_pos 
      pos_table$V1 <- new_pos
      pos_table$V2 <- win_spots_positions[new_pos, "V2"]
      
      neighbors_table <- win_prox_neighbors_table_func(pos_table,spots_clusters[row.names(pos_table),])
      return(neighbors_table)
    }) 
    
    win_abund <- sample_programs_composition(spots_clusters[spots_clusters$barcodes %in% win_spots,], gen_clusters) 
    
    programs_spatial_score <- lapply(sort(gen_clusters), function(cluster){
      if (!(cluster %in% spots_clusters$spot_type[spots_clusters$barcodes %in% win_spots]) | 
          length(spots_clusters$barcodes[spots_clusters$spot_type == cluster & spots_clusters$barcodes %in% win_spots]) < 10) {
        prog_score <- NaN
        return(prog_score)
      } else {
        program_neighbors_table = neighbors_table[row.names(neighbors_table) %in% spots_clusters$barcodes[spots_clusters$spot_type == cluster & spots_clusters$barcodes %in% win_spots],]
        obs <- obs_program_spatial_score(program_neighbors_table, cluster)
        one_spatial_score <- one_val(dim(program_neighbors_table)[1])
        zero_spatial_score <- zero_val(win_rand_neighbors_table, spots_clusters[spots_clusters$barcodes %in% win_spots,], cluster)
        
        if (win_abund[cluster] >= 0.6) {
          if (obs>one_spatial_score){obs <- one_spatial_score}
          return(obs/one_spatial_score)
        } else {
          if (obs>one_spatial_score){obs <- one_spatial_score}
          if (obs<zero_spatial_score){obs <- zero_spatial_score}
          
          prog_score <- (obs - zero_spatial_score)/(one_spatial_score - zero_spatial_score)
          return(prog_score)
        }
      }
    })
    return(mean(na.omit(unlist(programs_spatial_score))))
  })
  return(spatial_score_win)
}


# define dis-organized zones ------------------------------------------------------------

all_scores <- readRDS("Spatial_coh_zones/sc_windows/spatial_win5v3.rds")
hist(unlist(all_scores), breaks = 100, main = "radius 5: spatial coherence socres", xlab = "spatial coherence score")
th1 <- as.numeric(quantile(na.omit(unlist(all_scores)),probs = (seq(0,1,0.1)))[5])

all_zones5 <- sapply(c(1:length(sample_ls)), function(i){
  set_zones <- ifelse(all_scores[[i]] <= th1, "dis", "other")
  print(table(set_zones))
  return(set_zones)
})

all_scores <- readRDS("Spatial_coh_zones/sc_windows/spatial_win8v3.rds")
hist(unlist(all_scores)[unlist(all_scores) > 0 & unlist(all_scores) < 1], breaks = 100, main = "radius 8: spatial coherence socres", xlab = "spatial coherence score")
th1 <- as.numeric(quantile(na.omit(unlist(all_scores)),probs = (seq(0,1,0.1)))[5])

all_zones10 <- sapply(c(1:length(sample_ls)), function(i){
  set_zones <- ifelse(all_scores[[i]] <= th1, "dis", "other")
  print(table(set_zones))
  return(set_zones)
})

all_scores <- readRDS("Spatial_coh_zones/sc_windows/spatial_win11v3.rds")
hist(unlist(all_scores)[unlist(all_scores) > 0 & unlist(all_scores) < 1], breaks = 100, main = "radius 11: spatial coherence socres", xlab = "spatial coherence score")
th1 <- as.numeric(quantile(na.omit(unlist(all_scores)),probs = (seq(0,1,0.1)))[5])

all_zones15 <- sapply(c(1:length(sample_ls)), function(i){
  set_zones <- ifelse(all_scores[[i]] <= th1, "dis", "other")
  print(table(set_zones))
  return(set_zones)
})



zones_intersect <- sapply(c(1:26), function(x){
  inter1 <- intersect(names(all_zones5[[x]])[all_zones5[[x]] == "dis"],names(all_zones10[[x]])[all_zones10[[x]] == "dis"])
  inter2 <- intersect(inter1,names(all_zones15[[x]])[all_zones15[[x]] == "dis"])
  return(inter2)
})


all_zones <- sapply(c(1:26),function(i){
  set_zones <- ifelse(names(all_scores[[i]]) %in% zones_intersect[[i]],"dis","other")
  names(set_zones) <- names(all_scores[[i]])
  return(set_zones)
})

win_size <- 4

all_smooth <- sapply(c(1:26),function(i){
  print(sample_ls[i])
  spots_positions <- read.csv(paste("general/GBM_data/", sample_ls[i] , "/outs/spatial/tissue_positions_list.csv", sep = ""), header = FALSE, stringsAsFactors = FALSE)
  row.names(spots_positions) <- spots_positions$V1
  
  spots_clusters <- readRDS(paste("MP/mp_assign_124/", sample_ls[i], ".rds", sep = ""))
  spots_clusters <- na.omit(spots_clusters)
  colnames(spots_clusters) <- c("barcodes", "spot_type")
  row.names(spots_clusters)<- spots_clusters$barcodes
  
  neighbors_table <- prox_neighbors_table_func(spots_positions,spots_clusters)
  
  all_spots <- spots_clusters$barcodes
  
  smoothing <- sapply(all_spots, function(spot){
    win_spots <- c(spot)
    sapply(c(1:win_size), function(i){
      win_spots <<- unique(c(win_spots,unique(na.omit(as.character(neighbors_table[win_spots,])))))
    })
    win_zones <- all_zones[[i]][names(all_zones[[i]]) %in% win_spots]
    if (length(win_zones[win_zones == "dis"])/length(win_zones) >= 0.5) {return("dis")}
    else {return("other")}
  })
  return(smoothing)
})


# define struct zones -----------------------------------------------------

all_scores <- readRDS("Spatial_coh_zones/sc_windows/spatial_win5v3.rds")
hist(unlist(all_scores), breaks = 100)
th1 <- as.numeric(quantile(na.omit(unlist(all_scores)),probs = (seq(0,1,0.1)))[5])

all_zones5 <- sapply(c(1:length(sample_ls)), function(i){
  set_zones <- ifelse(all_scores[[i]] >= th1, "struct","other")
  print(table(set_zones))
  return(set_zones)
})

all_scores <- readRDS("Spatial_coh_zones/sc_windows/spatial_win8v3.rds")
hist(unlist(all_scores)[unlist(all_scores) > 0 & unlist(all_scores) < 1], breaks = 100)
th1 <- as.numeric(quantile(na.omit(unlist(all_scores)),probs = (seq(0,1,0.1)))[5])

all_zones10 <- sapply(c(1:length(sample_ls)), function(i){
  set_zones <- ifelse(all_scores[[i]] >= th1, "struct", "other")
  print(table(set_zones))
  return(set_zones)
})

all_scores <- readRDS("Spatial_coh_zones/sc_windows/spatial_win11v3.rds")
hist(unlist(all_scores)[unlist(all_scores) > 0 & unlist(all_scores) < 1], breaks = 100)
th1 <- as.numeric(quantile(na.omit(unlist(all_scores)),probs = (seq(0,1,0.1)))[5])

all_zones15 <- sapply(c(1:length(sample_ls)), function(i){
  set_zones <- ifelse(all_scores[[i]] >= th1, "struct", "other")
  print(table(set_zones))
  return(set_zones)
})



zones_intersect <- sapply(c(1:26), function(x){
  inter1 <- intersect(names(all_zones5[[x]])[all_zones5[[x]] == "struct"],names(all_zones10[[x]])[all_zones10[[x]] == "struct"])
  inter2 <- intersect(inter1,names(all_zones15[[x]])[all_zones15[[x]] == "struct"])
  return(inter2)
})


all_zones <- sapply(c(1:26),function(i){
  set_zones <- ifelse(names(all_scores[[i]]) %in% zones_intersect[[i]],"struct","other")
  names(set_zones) <- names(all_scores[[i]])
  return(set_zones)
})

win_size <- 4
all_smooth_struct <- sapply(c(1:26),function(i){
  print(sample_ls[i])
  spots_positions <- read.csv(paste("general/GBM_data/", sample_ls[i] , "/outs/spatial/tissue_positions_list.csv", sep = ""), header = FALSE, stringsAsFactors = FALSE)
  row.names(spots_positions) <- spots_positions$V1
  
  spots_clusters <- readRDS(paste("MP/mp_assign_124/", sample_ls[i], ".rds", sep = ""))
  spots_clusters <- na.omit(spots_clusters)
  colnames(spots_clusters) <- c("barcodes", "spot_type")
  row.names(spots_clusters)<- spots_clusters$barcodes
  
  neighbors_table <- prox_neighbors_table_func(spots_positions,spots_clusters)
  
  all_spots <- spots_clusters$barcodes
  
  smoothing <- sapply(all_spots, function(spot){
    win_spots <- c(spot)
    sapply(c(1:win_size), function(i){
      win_spots <<- unique(c(win_spots,unique(na.omit(as.character(neighbors_table[win_spots,])))))
    })
    win_zones <- all_zones[[i]][names(all_zones[[i]]) %in% win_spots]
    if (length(win_zones[win_zones == "struct"])/length(win_zones) >= 0.5) {return("struct")}
    else {return("other")}
  })
  return(smoothing)
})


# final zones -------------------------------------------------------------

all_zones <- sapply(c(1:26),function(i){
  sapply(names(all_smooth[[i]]), function(s){
    if (all_smooth[[i]][s] == "dis" & all_smooth_struct[[i]][s] == "struct") {
      return("Intermediate")
    } else if (all_smooth[[i]][s] == "dis") {
      return("dis")
    } else if (all_smooth_struct[[i]][s] == "struct") {
      return("struct")
    } else {
      return("Intermediate")
    }
  })
})

names(all_zones) <- sample_ls

is_small <- sapply(all_zones, function(x){
  tx <- table(x)
  
  old_clusters <- names(tx)
  add_clusters <- c("dis","Intermediate", "struct")[!c("dis","Intermediate", "struct") %in% old_clusters]
  if (length(add_clusters) > 0) {
    tx <- c(as.numeric(tx),0)
    names(tx) <- c(old_clusters, add_clusters)
  }
  tx <- tx[sort(names(tx))]
  tx10 <- sum(tx)*0.12
  print(tx)
  return(c(tx[1]<tx10, tx[2]<tx10,tx[3]<tx10))
})
colnames(is_small) <- sample_ls

dis2use <- c(1:26)[!is_small[1,]]
sturct2use <- c(1:26)[!is_small[3,]] 

spatial_score_zones <- lapply(c(1:length(sample_ls)), function(i){
  # load data
  print(sample_ls[i])
  spots_positions <- read.csv(paste("general/GBM_data/", sample_ls[i] , "/outs/spatial/tissue_positions_list.csv", sep = ""), header = FALSE, stringsAsFactors = FALSE)
  row.names(spots_positions) <- spots_positions$V1
  
  spots_clusters <- readRDS(paste("MP/mp_assign_124/", sample_ls[i], ".rds", sep = ""))
  spots_clusters <- na.omit(spots_clusters)
  colnames(spots_clusters) <- c("barcodes", "spot_type")
  row.names(spots_clusters)<- spots_clusters$barcodes  
  zone_score <- lapply(c("dis", "struct"), function(z){
    
    zspots <- na.omit(names(all_zones[[i]])[all_zones[[i]] == z]) # checkk!!!!!!!!!!!
    if (length(zspots) == 0) {
      return(NA)
    } else {
      # abundance 
      programs_comp <- sample_programs_composition(spots_clusters[spots_clusters$barcodes %in% zspots,],gen_clusters)
      
      # neighbors tables   
      
      neighbors_table <- neighbors_table_func(spots_positions,spots_clusters)
      zone_neighbors_table <- neighbors_table[zspots,] 
      
      win_rand_neighbors_table <- lapply(c(1:rand_num), function(i){
        zspots_positions <- spots_positions[spots_positions$V1 %in% zspots,]
        new_pos <- sample(zspots_positions$V1, length(zspots_positions$V1), replace = FALSE)
        pos_table <- zspots_positions
        row.names(pos_table) <- new_pos 
        pos_table$V1 <- new_pos
        pos_table$V2 <- zspots_positions[new_pos, "V2"]
        
        neighbors_table <- win_prox_neighbors_table_func(pos_table,spots_clusters[spots_clusters$barcodes %in% zspots,])
        return(neighbors_table)
      }) 
      
      # spatial coherence
      
      programs_spatial_score <- sapply(sort(gen_clusters), function(cluster){
        if (!(cluster %in% spots_clusters$spot_type[spots_clusters$barcodes %in% zspots]) | 
            length(spots_clusters$barcodes[spots_clusters$spot_type == cluster & spots_clusters$barcodes %in% zspots]) < 3) {
          prog_score <- NaN
        } else {
          program_neighbors_table = neighbors_table[row.names(neighbors_table) %in% spots_clusters$barcodes[spots_clusters$spot_type == cluster & spots_clusters$barcodes %in% zspots],]
          obs <- obs_program_spatial_score(program_neighbors_table, cluster)
          one_spatial_score <- one_val(dim(program_neighbors_table)[1])
          zero_spatial_score <- zero_val(win_rand_neighbors_table, spots_clusters[spots_clusters$barcodes %in% zspots,], cluster)
          if (obs>one_spatial_score){obs <- one_spatial_score}
          if (obs<zero_spatial_score){obs <- zero_spatial_score}
          
          prog_score <- (obs - zero_spatial_score)/(one_spatial_score - zero_spatial_score)
        }
        return(prog_score)
      })    
      return(list(programs_comp,programs_spatial_score))
    }
  })
  names(zone_score) <- c("dis", "struct")
  return(zone_score)
})


# plot Relationship between spot purity and spatial coherence -------------

purity_score_scaled <- readRDS("CNA/mal_lev.rds")

all_dis_purity <- sapply(c(1:26), function(i){
  dis_purity <- purity_score_scaled[[i]][all_zones[[i]] == "dis"]
  return(dis_purity)
})

hist(unlist(all_dis_purity))

all_st_purity <- sapply(c(1:26), function(i){
  st_purity <- purity_score_scaled[[i]][all_zones[[i]] == "struct"]
  return(st_purity)
})

hist(unlist(all_st_purity))


zones_class_fin <- sapply(c(1:26), function(i){
  print(sample_ls[[i]])
  
  spots_clusters <- readRDS(paste("MP/mp_assign_124/", sample_ls[[i]], ".rds", sep = ""))
  spots_clusters <- na.omit(spots_clusters)
  colnames(spots_clusters) <- c("barcodes", "spot_type")
  row.names(spots_clusters)<- spots_clusters$barcodes  
  
  v2 <- sapply(c(1:length(all_zones[[i]])), function(s){
    if(all_zones[[i]][s] == "dis" & purity_score_scaled[[i]][names(all_zones[[i]][s])] < 0.5) {
      return("dis_norm")
    } else if (all_zones[[i]][s] == "dis" & purity_score_scaled[[i]][names(all_zones[[i]][s])] >= 0.5) {
      return("dis_mal")
    } else if (all_zones[[i]][s] == "struct" & purity_score_scaled[[i]][names(all_zones[[i]][s])] < 0.4) {
      return("st_norm")
    } else if (all_zones[[i]][s] == "struct" & purity_score_scaled[[i]][names(all_zones[[i]][s])] >= 0.4) {
      return("st_mal")
    } else {
      return(all_zones[[i]][s])
    }
  })
  names(v2) <- names(all_zones[[i]])
  
  return(v2)
})
names(zones_class_fin) <- sample_ls


zone.colors <- c(structured_mal = "#cf5e4e",structured_norm = "#a82203", disorganized_mal = "#208cc0", disorganized_norm = "#003967",intermediate ="#f1af3a")
sapply(c(1:26), function(i){
  print(sample_ls[[i]])
  spots_positions <- read.csv(paste("general/GBM_data/", sample_ls[[i]] , "/outs/spatial/tissue_positions_list.csv", sep = ""), header = FALSE, stringsAsFactors = FALSE)
  row.names(spots_positions) <- spots_positions$V1
  
  spots_clusters <- readRDS(paste("MP/mp_assign_124/", sample_ls[[i]], ".rds", sep = ""))
  spots_clusters <- na.omit(spots_clusters)
  colnames(spots_clusters) <- c("barcodes", "spot_type")
  row.names(spots_clusters)<- spots_clusters$barcodes  
  
  spots_filt = spots_positions[spots_positions$V1 %in% spots_clusters$barcodes,]
  row2plot1 <- spots_clusters$barcodes[spots_clusters$barcodes %in% names(zones_class_fin[[i]])[zones_class_fin[[i]] == "st_mal"]]
  row2plot2 <- spots_clusters$barcodes[spots_clusters$barcodes %in% names(zones_class_fin[[i]])[zones_class_fin[[i]] == "st_norm"]]
  row2plot3 <- spots_clusters$barcodes[spots_clusters$barcodes %in% names(zones_class_fin[[i]])[zones_class_fin[[i]] == "dis_mal"]]
  row2plot4 <- spots_clusters$barcodes[spots_clusters$barcodes %in% names(zones_class_fin[[i]])[zones_class_fin[[i]] == "dis_norm"]]
  spots_filt$plot <- factor(ifelse(spots_filt$V1 %in% row2plot1, "structured_mal", 
                                   ifelse(spots_filt$V1 %in% row2plot2, "structured_norm", 
                                          ifelse(spots_filt$V1 %in% row2plot3, "disorganized_mal",
                                                 ifelse(spots_filt$V1 %in% row2plot4, "disorganized_norm","intermediate")))), levels = c("disorganized_mal","disorganized_norm", "structured_mal", "structured_norm","intermediate"))
  spots_filt$V3_ops = -(spots_filt$V3)
  
  
  gg = ggplot(spots_filt, aes(x=V4, y=V3_ops)) + 
    geom_point(aes(col=plot), size=2) + 
    labs(title=paste(sample_ls[[i]], "spatial coherence zones"), y="pos y", x="pos x") +
    scale_color_manual(values = zone.colors, name = "zone") + 
    theme_void()
  print(gg)
  
})

zones_cat <- c("dis_mal","dis_norm","Intermediate","st_mal","st_norm")
zone_abund <- sapply(zones_class_fin, function(z){
  spots_zone <- as.data.frame(z)
  spots_zone <- na.omit(spots_zone)
  spots_zone$barcodes <- row.names(spots_zone)
  colnames(spots_zone) <- c("spot_type","barcodes")
  
  # abundance 
  
  zones_comp <- sample_programs_composition(spots_zone,zones_cat)
  return(zones_comp)
})



######## Functions  --------------------------------------------------------------


sample_programs_composition <- function(spots_clusters, gen_clusters){
  composition <- table(spots_clusters$spot_type)
  old_clusters <- names(composition)
  add_clusters <- gen_clusters[!gen_clusters %in% old_clusters]
  sapply(add_clusters,function(clust){
    composition <<- c(composition, clust = 0)
  })
  
  names(composition) <- c(old_clusters, add_clusters)
  final_composition <- composition[sort(names(composition))]/sum(composition)
  return(final_composition)
}

obs_program_spatial_score <- function(program_neighbors, cluster){
  cluster_neighbors_bin <- ifelse(program_neighbors == cluster, 1, 0)
  if(is.null(dim(program_neighbors))){
    cluster_neighbors_sum <- sum(cluster_neighbors_bin)
  } else {
    cluster_neighbors_sum <- apply(cluster_neighbors_bin,1,function(rx){sum(na.omit(rx))})
  }
  obs <- mean(cluster_neighbors_sum)
  return(obs)
}

one_val <- function(spots_num){
  a <- sqrt((4*spots_num)/(6*sqrt(3)))
  oneval <- (6*spots_num-12*a-6)/spots_num
  return(oneval)
}

zero_val <- function(rand_table, spots_clusters, cluster){
  all_zeroval <- sapply(rand_table, function(neighbors_rand_table){
    program_rand_neighbors_table = neighbors_rand_table[row.names(neighbors_rand_table) %in% spots_clusters$barcodes[spots_clusters$spot_type == cluster],]
    rand_obs <- obs_program_spatial_score(program_rand_neighbors_table, cluster)
    return(rand_obs)
  })
  zeroval <- mean(all_zeroval)
  return(zeroval)
}


calc_adj_mat <- function(neighbored_cluster, spots_clusters, cluster_neighbors_table){
  if (!(neighbored_cluster %in% spots_clusters$spot_type)) {
    return(0)
  } else {
    cluster_neighbors_bin <- ifelse(cluster_neighbors_table == neighbored_cluster, 1, 0)
    
    if (is.null(dim(cluster_neighbors_bin))){
      cluster_neighbors_sum <- sum(na.omit(cluster_neighbors_bin))
    } else {
      cluster_neighbors_sum <- sum(apply(cluster_neighbors_bin,1,function(x){sum(na.omit(x))}))
    }
    return(cluster_neighbors_sum)
  }
}


prog_connectivity_score <- function(program_neighbors, cluster){
  cluster_neighbors_bin <- ifelse(program_neighbors == cluster, 0, 1)
  if(is.null(dim(program_neighbors))){
    prog_connect <- 1
  } else {
    prog_connect <- length(which(apply(cluster_neighbors_bin,1,function(x){sum(na.omit(x))}) >0))
  }
  return(prog_connect)
}

neighbors_table_func <- function(spots_positions,spots_clusters){
  neighbors_table <- sapply(spots_clusters$barcodes, function(spot){
    spots_row = spots_positions[spots_positions$V1 == spot, 3]
    spots_col = spots_positions[spots_positions$V1 == spot, 4]
    
    if (spots_col == 0 | spots_row == 0) {
      c1 = NaN
    } else {
      n1 = spots_positions$V1[spots_positions$V3 == spots_row - 1 & spots_positions$V4 == spots_col - 1]
      if (spots_positions$V2[spots_positions$V1 == n1] == 0 | !(n1 %in% spots_clusters$barcodes)){
        c1 = NaN
      } else {
        c1 = as.character(spots_clusters$spot_type[spots_clusters$barcodes == n1])
      }
    }
    
    if (spots_col == 127 | spots_row == 0) {
      c2 = NaN
    } else {
      n2 = spots_positions$V1[spots_positions$V3 == spots_row - 1 & spots_positions$V4 == spots_col + 1]
      if (spots_positions$V2[spots_positions$V1 == n2] == 0 | !(n2 %in% spots_clusters$barcodes)){
        c2 = NaN
      } else {
        c2 = as.character(spots_clusters$spot_type[spots_clusters$barcodes == n2])
      }
    }
    
    if (spots_col == 0 | spots_col == 1) {
      c3 = NaN
    } else {
      n3 = spots_positions$V1[spots_positions$V3 == spots_row & spots_positions$V4 == spots_col - 2]
      if (spots_positions$V2[spots_positions$V1 == n3] == 0 | !(n3 %in% spots_clusters$barcodes)){
        c3 = NaN
      } else {
        c3 = as.character(spots_clusters$spot_type[spots_clusters$barcodes == n3])
      }
    }
    
    if (spots_col == 126 | spots_col == 127) {
      c4 = NaN
    } else {
      n4 = spots_positions$V1[spots_positions$V3 == spots_row & spots_positions$V4 == spots_col + 2]
      if (spots_positions$V2[spots_positions$V1 == n4] == 0 | !(n4 %in% spots_clusters$barcodes)){
        c4 = NaN
      } else {
        c4 = as.character(spots_clusters$spot_type[spots_clusters$barcodes == n4])
      }
    }
    
    if (spots_col == 0 | spots_row == 77) {
      c5 = NaN
    } else {
      n5 = spots_positions$V1[spots_positions$V3 == spots_row + 1 & spots_positions$V4 == spots_col - 1]
      if (spots_positions$V2[spots_positions$V1 == n5] == 0 | !(n5 %in% spots_clusters$barcodes)){
        c5 = NaN
      } else {
        c5 = as.character(spots_clusters$spot_type[spots_clusters$barcodes == n5])
      }
    }
    
    if (spots_col == 127 | spots_row == 77) {
      c6 = NaN
    } else {
      n6 = spots_positions$V1[spots_positions$V3 == spots_row + 1 & spots_positions$V4 == spots_col + 1]
      if (spots_positions$V2[spots_positions$V1 == n6] == 0 | !(n6 %in% spots_clusters$barcodes)){
        c6 = NaN
      } else {
        c6 = as.character(spots_clusters$spot_type[spots_clusters$barcodes == n6])
      }
    }
    
    
    return(c(c1,c2,c3,c4,c5,c6))
    
  })
  
  neighbors_table = t(neighbors_table)
  row.names(neighbors_table) = spots_clusters$barcodes
  
  return(neighbors_table)
}

neighbors_table_funcV2 <- function(spots_positions,spots_clusters){
  neighbors_table <- sapply(spots_clusters$barcodes, function(spot){
    spots_row = spots_positions[spots_positions$V1 == spot, 3]
    spots_col = spots_positions[spots_positions$V1 == spot, 4]
    
    if (spots_col == 0 | spots_row == 0) {
      c1 = NaN
    } else {
      n1 = spots_positions$V1[spots_positions$V3 == spots_row - 1 & spots_positions$V4 == spots_col - 1]
      if (spots_positions$V2[spots_positions$V1 == n1] == 0 | !(n1 %in% spots_clusters$barcodes)){
        c1 = NaN
      } else {
        c1 = as.character(n1)
      }
    }
    
    if (spots_col == 127 | spots_row == 0) {
      c2 = NaN
    } else {
      n2 = spots_positions$V1[spots_positions$V3 == spots_row - 1 & spots_positions$V4 == spots_col + 1]
      if (spots_positions$V2[spots_positions$V1 == n2] == 0 | !(n2 %in% spots_clusters$barcodes)){
        c2 = NaN
      } else {
        c2 = as.character(n2)
      }
    }
    
    if (spots_col == 0 | spots_col == 1) {
      c3 = NaN
    } else {
      n3 = spots_positions$V1[spots_positions$V3 == spots_row & spots_positions$V4 == spots_col - 2]
      if (spots_positions$V2[spots_positions$V1 == n3] == 0 | !(n3 %in% spots_clusters$barcodes)){
        c3 = NaN
      } else {
        c3 = as.character(n3)
      }
    }
    
    if (spots_col == 126 | spots_col == 127) {
      c4 = NaN
    } else {
      n4 = spots_positions$V1[spots_positions$V3 == spots_row & spots_positions$V4 == spots_col + 2]
      if (spots_positions$V2[spots_positions$V1 == n4] == 0 | !(n4 %in% spots_clusters$barcodes)){
        c4 = NaN
      } else {
        c4 = as.character(n4)
      }
    }
    
    if (spots_col == 0 | spots_row == 77) {
      c5 = NaN
    } else {
      n5 = spots_positions$V1[spots_positions$V3 == spots_row + 1 & spots_positions$V4 == spots_col - 1]
      if (spots_positions$V2[spots_positions$V1 == n5] == 0 | !(n5 %in% spots_clusters$barcodes)){
        c5 = NaN
      } else {
        c5 = as.character(n5)
      }
    }
    
    if (spots_col == 127 | spots_row == 77) {
      c6 = NaN
    } else {
      n6 = spots_positions$V1[spots_positions$V3 == spots_row + 1 & spots_positions$V4 == spots_col + 1]
      if (spots_positions$V2[spots_positions$V1 == n6] == 0 | !(n6 %in% spots_clusters$barcodes)){
        c6 = NaN
      } else {
        c6 = as.character(n6)
      }
    }
    
    
    return(c(c1,c2,c3,c4,c5,c6))
    
  })
  
  neighbors_table = t(neighbors_table)
  row.names(neighbors_table) = spots_clusters$barcodes
  
  return(neighbors_table)
}


win_prox_neighbors_table_func <- function(spots_positions,spots_clusters){
  neighbors_table <- sapply(spots_clusters$barcodes, function(spot){
    spots_row = spots_positions[spots_positions$V1 == spot, 3]
    spots_col = spots_positions[spots_positions$V1 == spot, 4]
    
    if (spots_col == 0 | spots_row == 0) {
      n1 = NA
    } else {
      n1_temp = spots_positions$V1[spots_positions$V3 == spots_row - 1 & spots_positions$V4 == spots_col - 1]
      if (length(n1_temp) == 0) {
        n1 = NA
      } else if (spots_positions$V2[spots_positions$V1 == n1_temp] == 0 | !(n1_temp %in% spots_clusters$barcodes)){
        n1 = NA
      } else {
        n1 = as.character(spots_clusters$spot_type[spots_clusters$barcodes == n1_temp])
      }
    }
    
    if (spots_col == 127 | spots_row == 0) {
      n2 = NA
    } else {
      n2_temp = spots_positions$V1[spots_positions$V3 == spots_row - 1 & spots_positions$V4 == spots_col + 1]
      if (length(n2_temp) == 0) {
        n2 = NA
      } else if (spots_positions$V2[spots_positions$V1 == n2_temp] == 0 | !(n2_temp %in% spots_clusters$barcodes)){
        n2 = NA
      } else {
        n2 = as.character(spots_clusters$spot_type[spots_clusters$barcodes == n2_temp])
      }
    }
    
    if (spots_col == 0 | spots_col == 1) {
      n3 = NA
    } else {
      n3_temp = spots_positions$V1[spots_positions$V3 == spots_row & spots_positions$V4 == spots_col - 2]
      if (length(n3_temp) == 0) {
        n3 = NA
      } else if (spots_positions$V2[spots_positions$V1 == n3_temp] == 0 | !(n3_temp %in% spots_clusters$barcodes)){
        n3 = NA
      } else {
        n3 = as.character(spots_clusters$spot_type[spots_clusters$barcodes == n3_temp])
      }
    }
    
    if (spots_col == 126 | spots_col == 127) {
      n4 = NA
    } else {
      n4_temp = spots_positions$V1[spots_positions$V3 == spots_row & spots_positions$V4 == spots_col + 2]
      if (length(n4_temp) == 0) {
        n4 = NA
      } else if (spots_positions$V2[spots_positions$V1 == n4_temp] == 0 | !(n4_temp %in% spots_clusters$barcodes)){
        n4 = NA
      } else {
        n4 = as.character(spots_clusters$spot_type[spots_clusters$barcodes == n4_temp])
      }
    }
    
    if (spots_col == 0 | spots_row == 77) {
      n5 = NA
    } else {
      n5_temp = spots_positions$V1[spots_positions$V3 == spots_row + 1 & spots_positions$V4 == spots_col - 1]
      if (length(n5_temp) == 0) {
        n5 = NA
      } else if (spots_positions$V2[spots_positions$V1 == n5_temp] == 0 | !(n5_temp %in% spots_clusters$barcodes)){
        n5 = NA
      } else {
        n5 = as.character(spots_clusters$spot_type[spots_clusters$barcodes == n5_temp])
      }
    }
    
    if (spots_col == 127 | spots_row == 77) {
      n6 = NA
    } else {
      n6_temp = spots_positions$V1[spots_positions$V3 == spots_row + 1 & spots_positions$V4 == spots_col + 1]
      if (length(n6_temp) == 0) {
        n6 = NA
      } else if (spots_positions$V2[spots_positions$V1 == n6_temp] == 0 | !(n6_temp %in% spots_clusters$barcodes)){
        n6 = NA
      } else {
        n6 = as.character(spots_clusters$spot_type[spots_clusters$barcodes == n6_temp])
      }
    }
    
    
    return(c(n1,n2,n3,n4,n5,n6))
    
  })
  
  neighbors_table = t(neighbors_table)
  row.names(neighbors_table) = spots_clusters$barcodes
  
  return(neighbors_table)
}

prox_neighbors_table_func <- function(spots_positions,spots_clusters){
  neighbors_table <- sapply(spots_clusters$barcodes, function(spot){
    spots_row = spots_positions[spots_positions$V1 == spot, 3]
    spots_col = spots_positions[spots_positions$V1 == spot, 4]
    
    if (spots_col == 0 | spots_row == 0) {
      n1 = NA
    } else {
      n1_temp = spots_positions$V1[spots_positions$V3 == spots_row - 1 & spots_positions$V4 == spots_col - 1]
      if (spots_positions$V2[spots_positions$V1 == n1_temp] == 0 | !(n1_temp %in% spots_clusters$barcodes)){
        n1 = NA
      } else {
        n1 = n1_temp
      }
    }
    
    if (spots_col == 127 | spots_row == 0) {
      n2 = NA
    } else {
      n2_temp = spots_positions$V1[spots_positions$V3 == spots_row - 1 & spots_positions$V4 == spots_col + 1]
      if (spots_positions$V2[spots_positions$V1 == n2_temp] == 0 | !(n2_temp %in% spots_clusters$barcodes)){
        n2 = NA
      } else {
        n2 = n2_temp
      }
    }
    
    if (spots_col == 0 | spots_col == 1) {
      n3 = NA
    } else {
      n3_temp = spots_positions$V1[spots_positions$V3 == spots_row & spots_positions$V4 == spots_col - 2]
      if (spots_positions$V2[spots_positions$V1 == n3_temp] == 0 | !(n3_temp %in% spots_clusters$barcodes)){
        n3 = NA
      } else {
        n3 = n3_temp
      }
    }
    
    if (spots_col == 126 | spots_col == 127) {
      n4 = NA
    } else {
      n4_temp = spots_positions$V1[spots_positions$V3 == spots_row & spots_positions$V4 == spots_col + 2]
      if (spots_positions$V2[spots_positions$V1 == n4_temp] == 0 | !(n4_temp %in% spots_clusters$barcodes)){
        n4 = NA
      } else {
        n4 = n4_temp
      }
    }
    
    if (spots_col == 0 | spots_row == 77) {
      n5 = NA
    } else {
      n5_temp = spots_positions$V1[spots_positions$V3 == spots_row + 1 & spots_positions$V4 == spots_col - 1]
      if (spots_positions$V2[spots_positions$V1 == n5_temp] == 0 | !(n5_temp %in% spots_clusters$barcodes)){
        n5 = NA
      } else {
        n5 = n5_temp
      }
    }
    
    if (spots_col == 127 | spots_row == 77) {
      n6 = NA
    } else {
      n6_temp = spots_positions$V1[spots_positions$V3 == spots_row + 1 & spots_positions$V4 == spots_col + 1]
      if (spots_positions$V2[spots_positions$V1 == n6_temp] == 0 | !(n6_temp %in% spots_clusters$barcodes)){
        n6 = NA
      } else {
        n6 = n6_temp
      }
    }
    
    
    return(c(n1,n2,n3,n4,n5,n6))
    
  })
  
  neighbors_table = t(neighbors_table)
  row.names(neighbors_table) = spots_clusters$barcodes
  
  return(neighbors_table)
}
