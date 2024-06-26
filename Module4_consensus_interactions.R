library(ggvenn)

# Load data ---------------------------------------------------------------


file_list <- list.files("MP/mp_assign_124/")

gen_clusters <- as.character(unique(unlist(sapply(c(1:26), function(i){
  mp_assign <- readRDS(paste("MP/mp_assign_124/", file_list[i], sep = ""))
  return(unique(mp_assign$spot_type_meta_new))
}))))

pairs <- combn(sort(gen_clusters),2)
pairs_names <- apply(pairs, 2, function(x){return(paste(x[1],x[2], sep = " "))})

st_samp <-c("UKF243","UKF248","UKF255","UKF259","UKF260","UKF266","UKF269","UKF275","UKF304","UKF313","UKF334","ZH1019inf","ZH881T1","ZH916bulk","ZH916inf","ZH916T1","ZH1007nec","ZH1019T1","ZH8811Abulk","ZH8811Bbulk")
dis_samp <-c("MGH258","UKF251","UKF259","UKF260","UKF275","UKF296","UKF313","ZH1019inf","ZH881inf","ZH916inf","ZH1007inf","ZH1019T1", "ZH8812bulk")


# structured ----------------------------------------------------------------------

coloc_st <- readRDS("Summary_tables/colocal_structured_summary.rds")
conn_st <- readRDS("Summary_tables/conn_st.rds")
prox_st <- readRDS("Summary_tables/prox_st.rds")
(sum(coloc_st$enriched)+sum(conn_st$enr)+sum(prox_st[,5]$enr))/(length(pairs_names)*length(st_samp)*3)

coloc_st$mean_scaled <- ifelse(coloc_st$mean_enrichment < 1,
                               ((coloc_st$mean_enrichment - min(coloc_st$mean_enrichment))/(max(coloc_st$mean_enrichment)-min(coloc_st$mean_enrichment))) - 1,
                               ((coloc_st$mean_enrichment - min(coloc_st$mean_enrichment))/(max(coloc_st$mean_enrichment)-min(coloc_st$mean_enrichment)))) 
conn_st$mean_scaled <- ifelse(conn_st$mean_score < 0,
                              ((conn_st$mean_score - min(na.omit(conn_st$mean_score)))/(max(na.omit(conn_st$mean_score))-min(na.omit(conn_st$mean_score)))) - 1,
                              ((conn_st$mean_score - min(na.omit(conn_st$mean_score)))/(max(na.omit(conn_st$mean_score))-min(na.omit(conn_st$mean_score)))))


summary_coloc <- data.frame(pair = coloc_st$pair,
                            analysis = rep("coloc", 91),
                            prop = coloc_st$proportion_sig,
                            mean_scaled = coloc_st$mean_scaled)

new_pair <- as.character(sapply(summary_coloc$pair, function(p1){
  splt <- strsplit(p1, " ")[[1]]
  new_name <- sapply(row.names(conn_st),function(p2){
    splt2 <- strsplit(p2, " ")[[1]]
    if(splt[1] %in% splt2 & splt[2] %in% splt2) {
      return(p2)
    } else
      return(NA)
  })
  return(na.omit(new_name))
}))

summary_coloc$pair <- new_pair


summary_adj <- data.frame(pair = row.names(conn_st),
                          analysis = rep("adj", 91),
                          prop = conn_st$prop,
                          mean_scaled = conn_st$mean_scaled)

summary_prox5 <- data.frame(pair = row.names(conn_st),
                            analysis = rep("prox5", 91),
                            prop = prox_st[,5]$prop,
                            mean_scaled = prox_st[,5]$mean_score)

summary_prox5$mean_scaled <- ifelse(summary_prox5$mean_scaled < 0,
                                    ((summary_prox5$mean_scaled - min(na.omit(summary_prox5$mean_scaled)))/(max(na.omit(summary_prox5$mean_scaled))-min(na.omit(summary_prox5$mean_scaled)))) - 1,
                                    ((summary_prox5$mean_scaled - min(na.omit(summary_prox5$mean_scaled)))/(max(na.omit(summary_prox5$mean_scaled))-min(na.omit(summary_prox5$mean_scaled)))))


summary_prox8 <- data.frame(pair = row.names(conn_st),
                            analysis = rep("prox8", 91),
                            prop = prox_st[,8]$prop,
                            mean_scaled = prox_st[,8]$mean_score)

summary_prox8$mean_scaled <- ifelse(summary_prox8$mean_scaled < 0,
                                    ((summary_prox8$mean_scaled - min(na.omit(summary_prox8$mean_scaled)))/(max(na.omit(summary_prox8$mean_scaled))-min(na.omit(summary_prox8$mean_scaled)))) - 1,
                                    ((summary_prox8$mean_scaled - min(na.omit(summary_prox8$mean_scaled)))/(max(na.omit(summary_prox8$mean_scaled))-min(na.omit(summary_prox8$mean_scaled)))))

summary_prox15 <- data.frame(pair = row.names(conn_st),
                             analysis = rep("prox15", 91),
                             prop = prox_st[,15]$prop,
                             mean_scaled = prox_st[,15]$mean_score)

summary_prox15$mean_scaled <- ifelse(summary_prox15$mean_scaled < 0,
                                     ((summary_prox15$mean_scaled - min(na.omit(summary_prox15$mean_scaled)))/(max(na.omit(summary_prox15$mean_scaled))-min(na.omit(summary_prox15$mean_scaled)))) - 1,
                                     ((summary_prox15$mean_scaled - min(na.omit(summary_prox15$mean_scaled)))/(max(na.omit(summary_prox15$mean_scaled))-min(na.omit(summary_prox15$mean_scaled)))))


summary_df <- rbind(summary_coloc,summary_adj,summary_prox5,summary_prox8,summary_prox15)

pairs_mean <- sapply(unique(summary_df$pair),function(p){
  p_df <- summary_df[summary_df$pair == p,]
  return(mean(na.omit(p_df$mean_scaled)))
})

summary_df$pair_mean <- sapply(summary_df$pair, function(p){
  return(pairs_mean[p])
})

summary_df$pattern <- ifelse(summary_df$pair_mean > 0.35, "connected",
                             ifelse(summary_df$pair_mean < -0.35, "dis connected","mix"))
st_summary <- summary_df


# dis-organized  ---------------------------------------------------------------------


coloc_dis <- readRDS("Summary_tables/colocal_summary_dis.rds")
conn_dis <- readRDS("Summary_tables/conn_dis.rds")
prox_dis <- readRDS("Summary_tables/prox_dis.rds")
round((sum(coloc_dis$enriched)+sum(conn_dis$enr)+sum(prox_dis[,5]$enr))/(length(pairs_names)*length(dis_samp)*3),2)

coloc_dis$mean_scaled <- ifelse(coloc_dis$mean_enrichment < 1,
                                ((coloc_dis$mean_enrichment - min(coloc_dis$mean_enrichment))/(max(coloc_dis$mean_enrichment)-min(coloc_dis$mean_enrichment))) - 1,
                                ((coloc_dis$mean_enrichment - min(coloc_dis$mean_enrichment))/(max(coloc_dis$mean_enrichment)-min(coloc_dis$mean_enrichment)))) 
conn_dis$mean_scaled <- ifelse(conn_dis$mean_score < 0,
                               ((conn_dis$mean_score - min(na.omit(conn_dis$mean_score)))/(max(na.omit(conn_dis$mean_score))-min(na.omit(conn_dis$mean_score)))) - 1,
                               ((conn_dis$mean_score - min(na.omit(conn_dis$mean_score)))/(max(na.omit(conn_dis$mean_score))-min(na.omit(conn_dis$mean_score)))))


summary_coloc <- data.frame(pair = coloc_dis$pairs,
                            analysis = rep("coloc", 91),
                            prop = coloc_dis$prop_sig,
                            mean_scaled = coloc_dis$mean_scaled)

new_pair <- as.character(sapply(summary_coloc$pair, function(p1){
  splt <- strsplit(p1, " ")[[1]]
  new_name <- sapply(row.names(conn_dis),function(p2){
    splt2 <- strsplit(p2, " ")[[1]]
    if(splt[1] %in% splt2 & splt[2] %in% splt2) {
      return(p2)
    } else
      return(NA)
  })
  return(na.omit(new_name))
}))

summary_coloc$pair <- new_pair


summary_adj <- data.frame(pair = row.names(conn_dis),
                          analysis = rep("adj", 91),
                          prop = conn_dis$prop,
                          mean_scaled = conn_dis$mean_scaled)

summary_prox5 <- data.frame(pair = row.names(conn_dis),
                            analysis = rep("prox5", 91),
                            prop = prox_dis[,5]$prop,
                            mean_scaled = prox_dis[,5]$mean_score)

summary_prox5$mean_scaled <- ifelse(summary_prox5$mean_scaled < 0,
                                    ((summary_prox5$mean_scaled - min(na.omit(summary_prox5$mean_scaled)))/(max(na.omit(summary_prox5$mean_scaled))-min(na.omit(summary_prox5$mean_scaled)))) - 1,
                                    ((summary_prox5$mean_scaled - min(na.omit(summary_prox5$mean_scaled)))/(max(na.omit(summary_prox5$mean_scaled))-min(na.omit(summary_prox5$mean_scaled)))))


summary_prox8 <- data.frame(pair = row.names(conn_dis),
                            analysis = rep("prox8", 91),
                            prop = prox_dis[,8]$prop,
                            mean_scaled = prox_dis[,8]$mean_score)

summary_prox8$mean_scaled <- ifelse(summary_prox8$mean_scaled < 0,
                                    ((summary_prox8$mean_scaled - min(na.omit(summary_prox8$mean_scaled)))/(max(na.omit(summary_prox8$mean_scaled))-min(na.omit(summary_prox8$mean_scaled)))) - 1,
                                    ((summary_prox8$mean_scaled - min(na.omit(summary_prox8$mean_scaled)))/(max(na.omit(summary_prox8$mean_scaled))-min(na.omit(summary_prox8$mean_scaled)))))

summary_prox15 <- data.frame(pair = row.names(conn_dis),
                             analysis = rep("prox15", 91),
                             prop = prox_dis[,15]$prop,
                             mean_scaled = prox_dis[,15]$mean_score)

summary_prox15$mean_scaled <- ifelse(summary_prox15$mean_scaled < 0,
                                     ((summary_prox15$mean_scaled - min(na.omit(summary_prox15$mean_scaled)))/(max(na.omit(summary_prox15$mean_scaled))-min(na.omit(summary_prox15$mean_scaled)))) - 1,
                                     ((summary_prox15$mean_scaled - min(na.omit(summary_prox15$mean_scaled)))/(max(na.omit(summary_prox15$mean_scaled))-min(na.omit(summary_prox15$mean_scaled)))))


summary_df <- rbind(summary_coloc,summary_adj,summary_prox5,summary_prox8,summary_prox15)

pairs_mean <- sapply(unique(summary_df$pair),function(p){
  p_df <- summary_df[summary_df$pair == p,]
  return(mean(na.omit(p_df$mean_scaled)))
})

summary_df$pair_mean <- sapply(summary_df$pair, function(p){
  return(pairs_mean[p])
})

summary_df$pattern <- ifelse(summary_df$pair_mean > 0.35, "connected",
                             ifelse(summary_df$pair_mean < -0.35, "dis connected","mix"))

dis_summary <- summary_df

# st vs dis strongest couplings ---------------------------------------------------------------


dim(st_summary[st_summary$pattern == "connected" & st_summary$analysis %in% c("coloc", "adj", "prox5") & st_summary$mean_scaled > 0.35 & st_summary$prop >= 0.2,])
t_st <- table(as.character(st_summary$pair[st_summary$pattern == "connected" & st_summary$analysis %in% c("coloc", "adj") & st_summary$mean_scaled > 0.35 & st_summary$prop >= 0.2]))
colocOadj_st <- names(t_st)
colocPadj_st <- names(t_st[t_st > 1])
t_st2 <- table(as.character(st_summary$pair[st_summary$pattern == "connected" & st_summary$analysis %in% c("prox5","prox8","prox15") & st_summary$mean_scaled > 0.35 & st_summary$prop >= 0.2]))
prox_st <- names(t_st2)
st_pairs <- unique(c(colocPadj_st, intersect(colocOadj_st,prox_st)))


dim(dis_summary[dis_summary$pattern == "connected" & dis_summary$analysis %in% c("coloc", "adj", "prox5") & dis_summary$mean_scaled > 0.35 & dis_summary$prop >= 0.2,])
t_dis <- table(as.character(dis_summary$pair[dis_summary$pattern == "connected" & dis_summary$analysis %in% c("coloc", "adj") & dis_summary$mean_scaled > 0.35 & dis_summary$prop >= 0.2]))
colocOadj_dis <- names(t_dis)
colocPadj_dis <- names(t_dis[t_dis > 1])
t_dis2 <- table(as.character(dis_summary$pair[dis_summary$pattern == "connected" & dis_summary$analysis %in% c("prox5","prox8","prox15") & dis_summary$mean_scaled > 0.35 & dis_summary$prop >= 0.2]))
prox_dis <- names(t_dis2)
dis_pairs <- unique(c(colocPadj_dis, intersect(colocOadj_dis,prox_dis)))




# st vs dis strong coupling by scale  -------------------------------------


prox_st <- readRDS("Summary_tables/prox_st.rds")
summary_prox10 <- data.frame(pair = pairs_names,
                             analysis = rep("prox10", 91),
                             prop = prox_st[,10]$prop,
                             mean_scaled = prox_st[,10]$mean_score)

summary_prox10$mean_scaled <- ifelse(summary_prox10$mean_scaled < 0,
                                     ((summary_prox10$mean_scaled - min(na.omit(summary_prox10$mean_scaled)))/(max(na.omit(summary_prox10$mean_scaled))-min(na.omit(summary_prox10$mean_scaled)))) - 1,
                                     ((summary_prox10$mean_scaled - min(na.omit(summary_prox10$mean_scaled)))/(max(na.omit(summary_prox10$mean_scaled))-min(na.omit(summary_prox10$mean_scaled)))))


summary_prox12 <- data.frame(pair = pairs_names,
                             analysis = rep("prox12", 91),
                             prop = prox_st[,12]$prop,
                             mean_scaled = prox_st[,12]$mean_score)

summary_prox12$mean_scaled <- ifelse(summary_prox12$mean_scaled < 0,
                                     ((summary_prox12$mean_scaled - min(na.omit(summary_prox12$mean_scaled)))/(max(na.omit(summary_prox12$mean_scaled))-min(na.omit(summary_prox12$mean_scaled)))) - 1,
                                     ((summary_prox12$mean_scaled - min(na.omit(summary_prox12$mean_scaled)))/(max(na.omit(summary_prox12$mean_scaled))-min(na.omit(summary_prox12$mean_scaled)))))

st_coloc  <- length(st_summary$pair[st_summary$analysis == "coloc" & st_summary$mean_scaled > 0.35 & st_summary$prop >= 0.2])
st_adj  <- length(st_summary$pair[st_summary$analysis == "adj" & st_summary$mean_scaled > 0.35 & st_summary$prop >= 0.2])
st_p8  <- length(st_summary$pair[st_summary$analysis == "prox8" & st_summary$mean_scaled > 0.35 & st_summary$prop >= 0.2])
st_p10  <- length(summary_prox10$pair[summary_prox10$analysis == "prox10" & summary_prox10$mean_scaled > 0.35 & summary_prox10$prop >= 0.2])
st_p12  <- length(summary_prox12$pair[summary_prox12$analysis == "prox12" & summary_prox12$mean_scaled > 0.35 & summary_prox12$prop >= 0.2])
st_p15  <- length(st_summary$pair[st_summary$analysis == "prox15" & st_summary$mean_scaled > 0.35 & st_summary$prop >= 0.2])


prox_dis <- readRDS("Summary_tables/prox_dis.rds")
summary_prox10 <- data.frame(pair = pairs_names,
                             analysis = rep("prox10", 91),
                             prop = prox_dis[,10]$prop,
                             mean_scaled = prox_dis[,10]$mean_score)

summary_prox10$mean_scaled <- ifelse(summary_prox10$mean_scaled < 0,
                                     ((summary_prox10$mean_scaled - min(na.omit(summary_prox10$mean_scaled)))/(max(na.omit(summary_prox10$mean_scaled))-min(na.omit(summary_prox10$mean_scaled)))) - 1,
                                     ((summary_prox10$mean_scaled - min(na.omit(summary_prox10$mean_scaled)))/(max(na.omit(summary_prox10$mean_scaled))-min(na.omit(summary_prox10$mean_scaled)))))


summary_prox12 <- data.frame(pair = pairs_names,
                             analysis = rep("prox12", 91),
                             prop = prox_dis[,12]$prop,
                             mean_scaled = prox_dis[,12]$mean_score)

summary_prox12$mean_scaled <- ifelse(summary_prox12$mean_scaled < 0,
                                     ((summary_prox12$mean_scaled - min(na.omit(summary_prox12$mean_scaled)))/(max(na.omit(summary_prox12$mean_scaled))-min(na.omit(summary_prox12$mean_scaled)))) - 1,
                                     ((summary_prox12$mean_scaled - min(na.omit(summary_prox12$mean_scaled)))/(max(na.omit(summary_prox12$mean_scaled))-min(na.omit(summary_prox12$mean_scaled)))))


dis_coloc  <- length(dis_summary$pair[dis_summary$analysis == "coloc" & dis_summary$mean_scaled > 0.35 & dis_summary$prop >= 0.2])
dis_adj  <- length(dis_summary$pair[dis_summary$analysis == "adj" & dis_summary$mean_scaled > 0.35 & dis_summary$prop >= 0.2])
dis_p8  <- length(dis_summary$pair[dis_summary$analysis == "prox8" & dis_summary$mean_scaled > 0.35 & dis_summary$prop >= 0.2])
dis_p10  <- length(summary_prox10$pair[summary_prox10$analysis == "prox10" & summary_prox10$mean_scaled > 0.35 & summary_prox10$prop >= 0.2])
dis_p12  <- length(summary_prox12$pair[summary_prox12$analysis == "prox12" & summary_prox12$mean_scaled > 0.35 & summary_prox12$prop >= 0.2])
dis_p15  <- length(dis_summary$pair[dis_summary$analysis == "prox15" & dis_summary$mean_scaled > 0.35 & dis_summary$prop >= 0.2])

coupling_by_scale <- data.frame(num_coup = c(st_coloc,st_adj,st_p8,
                                             dis_coloc,dis_adj,dis_p8),
                                region = c(rep("st",3),rep("dis",3)),
                                analysis = rep(c("coloc","adj","prox8"),2))


coupling_by_scale$analysis <- factor(coupling_by_scale$analysis, levels = c("coloc","adj","prox8"))

ggplot(coupling_by_scale, aes(x=analysis, y=num_coup, group=region)) +
  geom_line(aes(color=region))+
  geom_point(aes(color=region))

st_p8_pairs  <- unique(st_summary$pair[st_summary$analysis %in% c("prox5","prox8","prox15") & st_summary$mean_scaled > 0.35 & st_summary$prop >= 0.2])
dis_p8_pairs  <- unique(dis_summary$pair[dis_summary$analysis %in% c("prox5","prox8","prox15") & dis_summary$mean_scaled > 0.35 & dis_summary$prop >= 0.2])


x <- list(
  struct = st_p8_pairs, 
  disorganised = dis_p8_pairs
)

ggvenn(
  x, 
  fill_color = c("white", "white"),
  stroke_size = 0.5, set_name_size = 4
)

stdis_int <- intersect(st_p8_pairs,dis_p8_pairs)
only_st <- st_p8_pairs[!(st_p8_pairs %in% stdis_int)]
only_dis <- dis_p8_pairs[!(dis_p8_pairs %in% stdis_int)]



