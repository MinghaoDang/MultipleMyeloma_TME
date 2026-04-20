source('global_config.R')

library('NMF')
library('igraph')
library('ggraph')

###################### Figure 4A NMF (k=5) heatmap ######################

meta <- readRDS("MM.TME.235samples.TME.meta.data.rds") %>%
  add_count(NewID, name = "n_cells") %>%
  filter(n_cells >= 200) %>%
  count(NewID, cell.status, name = "Freq") %>%
  complete(NewID, cell.status, fill = list(Freq = 0)) %>%
  group_by(NewID) %>%
  mutate(
    count = sum(Freq),
    percentage = Freq / count
  )

ratio <- meta %>%
  dplyr::select(NewID,cell.status,percentage) %>%
  tidyr::pivot_wider(
    names_from = cell.status,
    values_from = percentage,
    values_fill = 0
  ) %>%
  tibble::column_to_rownames("NewID") %>%
  as.matrix()

z_ratio <- scale(ratio) %>% as.data.frame() %>% t()

anno_use <- readRDS("sample.level.info.rds") %>%
  as.data.frame() %>%
  dplyr::select(
    NewID,
    Diagnosis,
    Primary_Cytogenetics,
    Gender,
    Race,
    Age
  ) %>%
  mutate(
    Age_group = case_when(
      Age <= 50 ~ "<=50",
      Age <= 60 ~ "<=60",
      Age <= 70 ~ "<=70",
      Age <= 80 ~ "<=80",
      TRUE ~ ">80"),
    
    Primary_Cytogenetics = case_when(
      Diagnosis == "nBM" ~ "Normal",
      Primary_Cytogenetics == "HY;t(11;14)" ~ "t(11;14)",
      Primary_Cytogenetics == "HY;t(14;16)" ~ "t(14;16)",
      Primary_Cytogenetics == "HY;t(4;14)" ~ "t(4;14)",
      Primary_Cytogenetics %in% c("NA", "Other") ~ "other/undetermined",
      TRUE ~ Primary_Cytogenetics)
  ) %>%
  filter(NewID %in% colnames(z_ratio)) %>%
  tibble::column_to_rownames("NewID") %>%
  dplyr::select(-Age)

anno_color<-list('Diagnosis'=c('nBM'='#4dac26','MGUS'='#1c75bc','SMM'='#fbb040','NDMM'='#d01c8b', 'RRMM'='#8856a7'),
                 'Primary_Cytogenetics'=c('Normal'='#00a651','t(4;14)'='#8dd3c7','t(6;14)'='#6a3d9a','t(11;14)'='#bebada','t(14;16)'='#80b1d3','t(14;20)'='#fccde5','HY'='#ff7f00','other/undetermined'='#878787'),
                 'Gender'=c('Male'='#9ecae1','Female'='#fa9fb5'),
                 'Race'=c('Asian'='Yellow', 'Black or African American'='Black', 'Other'='grey', 'White or Caucasian'='white'),
                 'Age_group'=c('<=50'='#f7fcfd','<=60'='#bfd3e6','<=70'='#8c96c6','<=80'='#88419d','>80'='#6e016b'))

nmf.rank<-readRDS('NMF.res_k5.rds')

B<-basis(nmf.rank)

group_map <- c(`1` = 5, `2` = 4, `3` = 1, `4` = 3, `5` = 2)

sample_group <- data.frame(
  SampleID = names(predict(nmf.rank, what = "samples")),
  group = predict(nmf.rank, what = "samples")
) %>%
  mutate(
    group = factor(unname(group_map[as.character(group)]), levels = 1:5),
    diagnosis = factor(gsub("[0-9]+", "", SampleID), levels = c("nBM", "MGUS", "SMM", "NDMM", "RRMM"))
  )

sample_order <- with(sample_group, order(group, diagnosis))

cell_group <- data.frame(
  SampleID = names(predict(nmf.rank, what = "features")),
  group = predict(nmf.rank, what = "features"),
  value = apply(B, 1, max)
) %>%
  mutate(
    group = factor(unname(group_map[as.character(group)]), levels = 1:5)
  )

cell_order <- with(cell_group, order(group, -value))


z_ratio_2<-z_ratio[cell_order,sample_order]

cell_group_split<-cell_group$group[cell_order]
sample_group_split<-sample_group$group[sample_order]

top_ann = HeatmapAnnotation( df=anno_use[colnames(z_ratio_2),],
                             which='column',
                             annotation_name_side = "right",
                             col=anno_color,
                             show_legend = T
)

plot_matrix<-z_ratio_2/4
plot_matrix[plot_matrix>0.5]<- 0.5
plot_matrix[plot_matrix< -0.5]<- -0.5

pdf('sample_NMF_k5.pdf',height=15,width=25)
ht_list<-Heatmap(plot_matrix, name = "ratio",
                 top_annotation = top_ann,
                 row_split = cell_group_split,
                 column_split = sample_group_split,
                 row_gap = unit(2, "mm"),
                 column_gap = unit(2, "mm"),
                 cluster_rows = F,
                 cluster_columns = F,
                 #column_order = order(factor(info.matrix$group, levels = c(1:8))),
                 row_names_gp = grid::gpar(fontsize = 10),
                 column_names_gp = grid::gpar(fontsize = 5))
draw(ht_list, annotation_legend_side = "right")
dev.off()


###################### Figure 4B NMF (k=5) spearman correlation between cell status in each cell status group ######################

meta <- readRDS("MM.TME.235samples.TME.meta.data.rds") %>%
  add_count(NewID, name = "n_cells") %>%
  filter(n_cells >= 200)


celltype <- meta %>%
  distinct(cell.type.level2, cell.status) %>%
  left_join(
    meta %>% count(cell.status, name = "Freq"),
    by = "cell.status"
  ) %>%
  mutate(log2Freq = log2(Freq))

df <- table(meta$NewID,meta$cell.status)
ratio <- as.data.frame(df / rowSums(df))
colnames(ratio) <- c("sampleID","cell.status","Freq")


row_clusters<-readRDS('NMF.k5.cell_group.rds')

rownames(row_clusters)<-row_clusters$SampleID
row_clusters$SampleID<-NULL
row_clusters$value<-NULL

size_range <- c(min(celltype$log2Freq),max(celltype$log2Freq))
weight_range <- c(0.5,1)

gp.list<-NULL
cell.use.list<-NULL
for(i in unique(row_clusters$group)){
  
  cell.type.use<-rownames(subset(row_clusters,group==i))
  ratio.use<-subset(ratio,cell.status%in%cell.type.use)
  
  ratio.use <- dcast(ratio.use, sampleID ~ ratio.use$cell.status, value.var = "Freq")
  rownames(ratio.use)<-ratio.use$sampleID
  ratio.use<-ratio.use[,-1]
  
  col_cor_mat <- cor(ratio.use, method = "spearman")
  
  cell.use.list[[length(cell.use.list)+1]]<-rownames(col_cor_mat)[rowSums(col_cor_mat>0.5)>1]
  
  col_cor_mat[lower.tri(col_cor_mat, diag = TRUE)] <- 0
  col_cor_mat[abs(col_cor_mat) < 0.5] <- 0
  
  edge_list <- which(col_cor_mat != 0, arr.ind = TRUE)
  edges <- data.frame(
    from = rownames(col_cor_mat)[edge_list[, 1]],
    to = colnames(col_cor_mat)[edge_list[, 2]],
    weight = col_cor_mat[edge_list]
  )
  
  # Optional: node metadata for coloring
  node_df <- data.frame(
    nodes = unique(c(edges$from,edges$to))
  )
  node_df<-left_join(node_df,celltype,by=c('nodes'='cell.status'))
  
  # Create igraph object
  graph <- graph_from_data_frame(d = edges, vertices = node_df, directed = FALSE)
  
  
  plotx<-ggraph(graph, layout = "circle") +
    geom_edge_link(aes(width = abs(weight)), color = "grey60", alpha = 0.6) +
    geom_node_point(aes(fill = cell.type.level2, size=log2Freq), shape = 21, color = "black") +
    geom_node_text(aes(label = name), repel = TRUE, size = 3, fontface = "bold") +
    scale_fill_manual(values = cell.type.level2.color) +
    scale_edge_width(limits = weight_range,range = c(0.5, 3)) +
    scale_size(limits = size_range,range = c(10, 40))+
    theme_void() +
    ggtitle("Correlation Network (Spearman)") +
    theme(legend.position = "none")
  
  gp.list[[length(gp.list)+1]]<-plotx
  
}
names(cell.use.list)<-c(1:5)
saveRDS(cell.use.list,'cell.use.in.4B.list.rds')

do.call(ggpubr::ggarrange,c(gp.list,ncol = 5,nrow = 1)) -> combined.gp
pdf('k5.cell.status.correlation.graph.in.each.group.pdf',height=8,width = 49)
print(combined.gp)
dev.off()

###################### Figure 4C, S6D module profile and compositions ######################

res_k<-readRDS('NMF.res_k5.rds')

W <- basis(res_k)
colnames(W)<-c('E5','E4','E1','E3','E2')

cell.order<-readRDS('NMF.k5.cell_group.rds')
W<-W[cell.order$SampleID,c('E1','E2','E3','E4','E5')]
#W[W>1]<-1
W<-t(scale(t(W)))


pdf('Cell_Subset.Composition.Profiles.of.Ecotypes.pdf',height=8,width=8)
ht_list<-Heatmap(W, name = "Scaled Contribution",
                 #top_annotation = top_ann,
                 #row_split = cell.order,
                 #column_split = sample.group,
                 row_gap = unit(2, "mm"),
                 column_gap = unit(2, "mm"),
                 cluster_rows = F,
                 cluster_columns = F,
                 #column_order = order(factor(info.matrix$group, levels = c(1:8))),
                 row_names_gp = grid::gpar(fontsize = 10),
                 column_names_gp = grid::gpar(fontsize = 10))
draw(ht_list, annotation_legend_side = "right")
dev.off()



H <- coef(res_k)
rownames(H)<-c('E5','E4','E1','E3','E2')
sample.order<-readRDS('NMF.k5.sample_group.rds')
H<-H[c('E1','E2','E3','E4','E5'),sample.order$SampleID]
H<-scale(H)

info<-readRDS('sample.level.info.rds')
anno<-info[,c('NewID','Ecotype')]
anno<-anno[!is.na(anno$Ecotype),]

rownames(anno)<-anno$NewID
anno$NewID<-NULL

anno_color<-list('Ecotype'=c('E1'='#a6cee3','E2'='#b2df8a','E3'='#fb9a99','E4'='#fdbf6f', 'E5'='#cab2d6'))

top_ann = HeatmapAnnotation( df=anno[colnames(H),,drop=F],
                             which='column',
                             annotation_name_side = "right",
                             col=anno_color,
                             show_legend = T
)


pdf('signature.contribution.per.patient.heatmap.pdf',height=5,width=15)
ht_list<-Heatmap(H, name = "Scaled Contribution",
                 top_annotation = top_ann,
                 #row_split = cell.order,
                 #column_split = sample.group,
                 row_gap = unit(2, "mm"),
                 column_gap = unit(2, "mm"),
                 cluster_rows = F,
                 cluster_columns = F,
                 #column_order = order(factor(info.matrix$group, levels = c(1:8))),
                 row_names_gp = grid::gpar(fontsize = 10),
                 column_names_gp = grid::gpar(fontsize = 10))
draw(ht_list, annotation_legend_side = "right")
dev.off()

H <- coef(res_k)
rownames(H)<-c('E5','E4','E1','E3','E2')
df <- melt(H)
colnames(df) <- c("Signature","Sample","Value")

df<-left_join(df,info[,c('NewID','Ecotype')],by=c('Sample'='NewID'))
df$Signature<-factor(df$Signature,levels=c('E1','E2','E3','E4','E5'))
df$Sample<-factor(df$Sample,levels=sample.order$SampleID)

plotx<-ggplot(df, aes(x = Sample, y = Value, fill = Signature)) +
  geom_bar(stat="identity", position="fill") +
  scale_fill_manual(values=c('E1'='#a6cee3','E2'='#b2df8a','E3'='#fb9a99','E4'='#fdbf6f', 'E5'='#cab2d6'))+
  facet_grid(.~Ecotype,scales = 'free',space = 'free') +
  ylab("Contribution proportion") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(angle = 90,size = 5,vjust=0.2,hjust=0.95),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        legend.title=element_text(size=0), 
        legend.text=element_text(size=12),
        legend.key = element_rect(colour = NA, fill = NA),
        strip.text.x = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5,size=20))

ggsave('patient.level.signature.contribution.pdf',useDingbats=F,plotx,width = 15,height=5)

###################### Figure 4D cell status proportion compare between ecotype ######################

meta <- readRDS("MM.TME.235samples.TME.meta.data.rds") %>%
  add_count(NewID, name = "n_cells") %>%
  filter(n_cells >= 200) %>%
  count(NewID, cell.status, name = "Freq") %>%
  complete(NewID, cell.status, fill = list(Freq = 0)) %>%
  group_by(NewID) %>%
  mutate(
    count = sum(Freq),
    percentage = Freq / count
  )

sample_group<-readRDS('NMF.k5.sample_group.rds')

data<-left_join(meta,sample_group,by=c('NewID'='SampleID'))

data.use<-subset(data,cell.status%in%c('HSC','CD4T_Tn','CD8T_Teff_FGFBP2','CD14CD16_Mono','CD16_Mono_JUND','CD4T_Th17') & diagnosis!='nBM')

plotx<-ggplot(data.use, aes(x = Ecotype, y = percentage)) +
  geom_boxplot(alpha = 0.8, show.legend = T,width = 0.75, aes(fill=Ecotype), outlier.shape=NA) + 
  geom_jitter(position=position_jitter(width=0.1, height=0), size=0.5) +
  stat_compare_means(method = "kruskal.test") + # Add p-values
  scale_fill_manual(values=c('E1'='#a6cee3','E2'='#b2df8a','E3'='#fb9a99','E4'='#fdbf6f', 'E5'='#cab2d6'))+
  facet_grid(rows = vars(cell.status),cols = vars(diagnosis),scales='free') +
  scale_y_continuous(expand = expansion(mult = c(0.15)))+
  theme(strip.text.x = element_text(size=10, color="black", face="bold"),
        axis.text.x = element_text(angle = 90,size = 10,vjust=0.2,hjust=0.95,colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        text = element_text(size=12)) +
  labs(x="", y="% (of TME cells)") 
ggsave('cell.status.in.ecotype.compare.pdf',plotx,width=8,height=12,limitsize = FALSE)

###################### Figure S6A choose k ######################

meta <- readRDS("MM.TME.235samples.TME.meta.data.rds") %>%
  add_count(NewID, name = "n_cells") %>%
  filter(n_cells >= 200) %>%
  count(NewID, cell.status, name = "Freq") %>%
  complete(NewID, cell.status, fill = list(Freq = 0)) %>%
  group_by(NewID) %>%
  mutate(
    count = sum(Freq),
    percentage = Freq / count
  )

ratio <- meta %>%
  dplyr::select(NewID,cell.status,percentage) %>%
  tidyr::pivot_wider(
    names_from = cell.status,
    values_from = percentage,
    values_fill = 0
  ) %>%
  tibble::column_to_rownames("NewID") %>%
  as.matrix()

scale_ratio <- apply(ratio, MARGIN = 2, function(x) (x-min(x))/(max(x)-min(x)))
scale_ratio <- as.data.frame(scale_ratio) %>% t()

ranks <- 2:10
estim.coad <- nmf(scale_ratio, ranks, nrun=500,method = "lee")
saveRDS(estim.coad,'estim.coad.rds')

pdf('NMF_rank_survey.pdf',height=6,width=10)
plot(estim.coad)
dev.off()

###################### Figure S6B consensusmap ######################
k_vec<-c(2:10)
seed = 123

nmf_list <- lapply(k_vec, function(k){
  cat("Running k =", k, "\n")
  nmf(scale_ratio, k, nrun = 100, seed = 123, method = "lee")
})

names(nmf_list) <- as.character(k_vec)

saveRDS(nmf_list,'NMF.res_k2-10.rds')

pdf("k2-10.consensusmap.pdf", width = 6.5, height = 6)
for(i in c(1:9)){

  res_k<-nmf_list[[i]]
  consensusmap(res_k, main = paste0("Consensus Clustering (k=",i+1,")"))
}
dev.off()

###################### Figure S6C co-occurrences of cell subsets in NMF ######################

meta <- readRDS("MM.TME.235samples.TME.meta.data.rds") %>%
  add_count(NewID, name = "n_cells") %>%
  filter(n_cells >= 200)

all_samples<-unique(meta$NewID)

module_merge <- list()

for(i in 1:200){
    #remove 20% samples randomly
    random.samples <- sample(all_samples, length(all_samples)*0.8)
    
    meta.use <- meta[meta$NewID %in% random.samples,]
    
    meta.use <- meta.use %>%
      count(NewID, cell.status, name = "Freq") %>%
      complete(NewID, cell.status, fill = list(Freq = 0)) %>%
      group_by(NewID) %>%
      mutate(
        count = sum(Freq),
        percentage = Freq / count
      )
    
    ratio <- meta.use %>%
      dplyr::select(NewID,cell.status,percentage) %>%
      tidyr::pivot_wider(
        names_from = cell.status,
        values_from = percentage,
        values_fill = 0
      ) %>%
      tibble::column_to_rownames("NewID") %>%
      as.matrix()
    
    scale_ratio <- apply(ratio, MARGIN = 2, function(x) (x-min(x))/(max(x)-min(x)))
    scale_ratio <- as.data.frame(scale_ratio) %>% t()
    
    seed = 123
    for(rk in 2:10){
      nmf.rank5 <- nmf(scale_ratio, 
                       rank = rk, 
                       nrun = 100,
                       seed = seed,
                       method = "lee")
      
      aa<-predict(nmf.rank5,what='features')
      bb<-split(names(aa),aa)
      names(bb)<-NULL
      module_merge<-c(module_merge,bb)
      
    }
    
}

statuses <- unique(meta$cell.status)

Mat <- outer(
  statuses, statuses,
  Vectorize(function(x, y) {
    sum(sapply(module_merge, function(m) x %in% m & y %in% m))
  })
)

rownames(Mat) <- statuses
colnames(Mat) <- statuses

#saveRDS(Mat,'robustness_of_NMF_Mat.rds')

custom_magma <- c(colorRampPalette(c('white', rev(magma(50, begin = 0.15))[1]))(10), rev(magma(50, begin = 0.18)))

cols <- rev(magma(256, begin = 0.5, end = 0.95))

pdf('robustness_of_NMF.pdf',height = 10,width = 12)
pheatmap(as.matrix(Mat), cluster_cols=T, cluster_rows=T,
         clustering_distance_rows="euclidean",color=custom_magma,
         fontsize=10,treeheight_row=0,treeheight_col=30,
         #cellheight = 7,cellwidth = 7,
         show_rownames=T,clustering_method = "ward.D2",
         show_colnames=F, border_color = NA)
dev.off()


###################### Figure S7A NMF heatmap by diagnosis ######################

meta <- readRDS("MM.TME.235samples.TME.meta.data.rds") %>%
  add_count(NewID, name = "n_cells") %>%
  filter(n_cells >= 200)

info<-readRDS("sample.level.info.rds")
info<-info[!is.na(info$Ecotype),]
rownames(info)<-info$NewID

anno_color<-list('Diagnosis'=c('nBM'='#4dac26','MGUS'='#1c75bc','SMM'='#fbb040','NDMM'='#d01c8b', 'RRMM'='#8856a7'),
                 'Ecotype'=c('E1'='#a6cee3','E2'='#b2df8a','E3'='#fb9a99','E4'='#fdbf6f', 'E5'='#cab2d6'))

for(i in c('MGUS','SMM','NDMM','RRMM')){
  
  meta.use<-subset(meta,Diagnosis==i)
  
  scale_ratio <- meta.use %>%
    count(NewID, cell.status, name = "n") %>%
    group_by(NewID) %>%
    mutate(Freq = n / sum(n)) %>%
    ungroup() %>%
    select(NewID, cell.status, Freq) %>%
    pivot_wider(
      names_from = cell.status,
      values_from = Freq,
      values_fill = 0
    ) %>%
    column_to_rownames("NewID") %>%
    as.matrix() %>%
    apply(2, \(x) (x - min(x)) / (max(x) - min(x))) %>%
    t() %>%
    as.data.frame()
  
  ranks <- 2:10
  estim.coad <- nmf(scale_ratio, ranks, nrun=500) # ,method = "lee"
  saveRDS(estim.coad,paste0(i,'.estim.coad.rds'))
  
  
  nmf_list <- lapply(ranks, function(k){
    cat("Running k =", k, "\n")
    nmf(scale_ratio, k, nrun = 100, seed = 123, method = "lee")
  })
  names(nmf_list) <- as.character(ranks)
  saveRDS(nmf_list,paste0(i,'.res_k2-10.rds'))
  
}

pdf('sample_by_diagnosis_NMF_k5.pdf',height=15,width=25)
for(i in c('MGUS','SMM','NDMM','RRMM')){
  
  meta.use<-subset(meta,Diagnosis==i)
  
  info.use<-subset(info,Diagnosis==i)
  
  df <- table(meta.use$NewID,meta.use$cell.status)
  ratio <- as.data.frame(df / rowSums(df))
  colnames(ratio) <- c("sampleID","cell.type","Freq")
  
  ratio <- dcast(ratio, sampleID ~ ratio$cell.type, value.var = "Freq")
  
  rownames(ratio) <- ratio$sampleID
  
  ratio <- ratio[,-1]
  ratio[is.na(ratio)] <- 0
  
  z_ratio <- scale(ratio)
  z_ratio <- as.data.frame(z_ratio)
  z_ratio <- t(z_ratio)
  
  nmf_list<-readRDS(paste0(i,'.res_k2-10.rds'))
  
  for(j in 1:length(nmf_list)){
    
    nmf.rank<-nmf_list[[j]]
    
    B<-basis(nmf.rank)
    
    sample.group <- predict(nmf.rank,what='samples')
    sample.group.2<-data.frame('SampleID'=names(sample.group),'group'=sample.group)
    sample.order<-order(sample.group.2[[2]])
    
    cell.group <- predict(nmf.rank,what='features')
    cell.group.2<-data.frame('SampleID'=names(cell.group),'group'=cell.group,'value'=apply(B,1,max))
    cell.group.2<-cell.group.2[cell.group.2$value>quantile(cell.group.2$value,0.25),]
    z_ratio_2<-z_ratio[cell.group.2$SampleID,]
    
    cell.order<-order(cell.group.2[[2]],-cell.group.2[[3]])
    
    z_ratio_2<-z_ratio_2[cell.order,sample.order]
    
    cell.group<-cell.group.2$group[cell.order]
    sample.group<-sample.group.2$group[sample.order]
    
    
    top_ann = HeatmapAnnotation( df=info.use[colnames(z_ratio_2),c('Diagnosis','Ecotype')],
                                 which='column',
                                 annotation_name_side = "right",
                                 col=anno_color,
                                 show_legend = T)
    
    plot_matrix<-z_ratio_2/4
    plot_matrix[plot_matrix>0.5]<- 0.5
    plot_matrix[plot_matrix< -0.5]<- -0.5
    
    ht_list<-Heatmap(plot_matrix, name = "ratio",
                     top_annotation = top_ann,
                     row_split = cell.group,
                     column_split = sample.group,
                     row_gap = unit(2, "mm"),
                     column_gap = unit(2, "mm"),
                     cluster_rows = F,
                     cluster_columns = F,
                     column_title=paste0(i,' k=',j+1),
                     #column_order = order(factor(info.matrix$group, levels = c(1:8))),
                     row_names_gp = grid::gpar(fontsize = 10),
                     column_names_gp = grid::gpar(fontsize = 5))
    draw(ht_list, annotation_legend_side = "right")
    
  }
  
}
dev.off()

###################### Figure S7B 53 normal + this study NMF heatmap (k=5) ######################

meta <- readRDS("MM.TME.235samples.TME.meta.data.rds") %>%
  dplyr::select(NewID, cell.status, Diagnosis)

public_nBM<-readRDS("53public_nBM.meta.data.rds") %>%
  select(orig.ident, predicted.id, Diagnosis) %>%
  rename(
    NewID = orig.ident,
    cell.status = predicted.id
  )

meta<-rbind(meta,public_nBM)

meta<- meta %>%
  add_count(NewID, name = "n_sample") %>%
  add_count(cell.status, name = "n_status") %>%
  filter(
    n_sample >= 200,
    n_status >= 50
  ) %>%
  select(-n_sample, -n_status)

ratio <- meta %>%
  count(NewID, cell.status, name = "n") %>%
  group_by(NewID) %>%
  mutate(Freq = n / sum(n)) %>%
  ungroup() %>%
  select(-n) %>%
  pivot_wider(
    names_from = cell.status,
    values_from = Freq,
    values_fill = 0
  ) %>%
  column_to_rownames("NewID") %>%
  as.data.frame()

z_ratio <- scale(ratio) %>% as.data.frame() %>% t()

info.use<-unique(meta[,c('NewID','Diagnosis')])
rownames(info.use)<-info.use$NewID
info.use$NewID<-NULL


nmf_list<-readRDS('MM.TME.235samples.and.53public_nBM.NMF.res_k2-10.rds')
nmf.rank<-nmf_list[[4]]
B<-basis(nmf.rank)

group_map <- c(`1` = 3, `2` = 4, `3` = 2, `4` = 1, `5` = 5)

sample_group <- data.frame(
  SampleID = names(predict(nmf.rank, what = "samples")),
  group = predict(nmf.rank, what = "samples")
) %>%
  mutate(group = unname(group_map[as.character(group)])) %>%
  left_join(
    meta %>% distinct(NewID, Diagnosis),
    by = c("SampleID" = "NewID")
  ) %>%
  mutate(
    Diagnosis = factor(Diagnosis, levels = c("nBM", "MGUS", "SMM", "NDMM", "RRMM")),
    group = factor(group, levels = 1:5)
  )

sample_order <- with(sample_group, order(group, Diagnosis))

cell_group <- data.frame(
  SampleID = names(predict(nmf.rank, what = "features")),
  group = predict(nmf.rank, what = "features"),
  value = apply(B, 1, max)
) %>%
  mutate(
    group = factor(unname(group_map[as.character(group)]), levels = 1:5)
  )

cell_order <- with(cell_group, order(group, -value))

z_ratio_2 <- z_ratio[cell_group$SampleID[cell_order], sample_group$SampleID[sample_order]]

cell_group_split <- cell_group$group[cell_order]
sample_group_split <- sample_group$group[sample_order]

top_ann = HeatmapAnnotation( df=info.use[colnames(z_ratio_2),,drop=F],
                             which='column',
                             annotation_name_side = "right",
                             col=list('Diagnosis'=c('nBM'='#4dac26','MGUS'='#1c75bc','SMM'='#fbb040','NDMM'='#d01c8b', 'RRMM'='#8856a7')),
                             show_legend = T)

plot_matrix<-z_ratio_2/4
plot_matrix[plot_matrix>0.5]<- 0.5
plot_matrix[plot_matrix< -0.5]<- -0.5

pdf('MM.TME.235samples.and.53public_nBM.k5.heatmap.pdf',height=15,width=25)

ht_list<-Heatmap(plot_matrix, name = "ratio",
                 top_annotation = top_ann,
                 row_split = cell_group_split,
                 column_split = sample_group_split,
                 row_gap = unit(2, "mm"),
                 column_gap = unit(2, "mm"),
                 cluster_rows = F,
                 cluster_columns = F,
                 row_names_gp = grid::gpar(fontsize = 10),
                 column_names_gp = grid::gpar(fontsize = 5))
draw(ht_list, annotation_legend_side = "right")

dev.off()


colnames(sample_group)[2]<-'Ecotype'

data_plot <- sample_group %>%
  count(Ecotype, Diagnosis, name = "Freq") %>%
  group_by(Diagnosis) %>%
  mutate(
    count = sum(Freq),
    percentage = round(Freq * 100 / count, 2),
    label_text = paste0(round(percentage, 1), "%")
  ) %>%
  ungroup()

plotx<-ggplot(data_plot, aes(x = 2, y = percentage, fill = Ecotype)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +  # makes the hole
  facet_grid(. ~ Diagnosis) +
  geom_text(aes(label = label_text),
            position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(values = c(
    '1'='#a6cee3','2'='#b2df8a','3'='#fb9a99','4'='#fdbf6f', '5'='#cab2d6'
  )) +
  theme_void() +
  theme(legend.position = "bottom")

ggsave('normal.ecotype.piechart.pdf',useDingbats=F,plotx,width = 18,height=6)

###################### Figure S8A-C 4 independent datasets umap, cell count, stacked barplot ######################

# # label transfer
# TME.obj.CD8T<-readRDS('TME.4datasets.CD8T.scs.rds')
# ref.CD8T.obj<-readRDS('MM.TME.235samples.CD8T.filter.scs.rds')
# CD8T.anchors <- FindTransferAnchors(reference = ref.CD8T.obj, query = TME.obj.CD8T, dims = 1:30,
#                                     reference.reduction = "pca")
# predictions <- TransferData(anchorset = CD8T.anchors, refdata = ref.CD8T.obj$cell.status, dims = 1:30)
# TME.obj.CD8T <- AddMetaData(TME.obj.CD8T, metadata = predictions)

TME.obj<-readRDS('TME.4datasets.filter.scs.rds')

cell.type.level2.color.2<-c('Mature_B'='#1f78b4',
                          'CD8T'='#33a02c','CD4T'='#936BD8','NK'='#a6cee3',
                          'CD14_Mono'='#ff7f00','CD14CD16_Mono'='#b2df8a','CD16_Mono'='#7D6294','Macrophage'='#E1D1AF',
                          'cDC'='#D947E4','pDC'='#dc9ae9',
                          'Platelet'='#6a3d9a','Erythrocytes'='#fc4e2a',
                          'Progenitors'='#f4a900',
                          'Stroma'='#b15928')

p1<-DimPlot(object = TME.obj, reduction = "umap",pt.size = 0.5,group.by = 'project',label = T,label.size = 6)+
  scale_color_manual(values=c('GSE161801'='#beaed4','GSE223060'='#fdc086','GSE271107'='#7fc97f','GSE124310'='#386cb0'))+
  theme_umap()

p2<-DimPlot(object = TME.obj, reduction = "umap",pt.size = 0.5,group.by = 'cell.type.level2',label = T,label.size = 6)+
  scale_color_manual(values=cell.type.level2.color.2)+
  theme_umap()

do.call(ggarrange,c(list(p1,p2),ncol = 1,nrow = 1)) -> combined.gp
pdf(paste0('TME.4datasets.umap.',Sys.Date(),'.pdf'),height=8,width = 8)
print(combined.gp)
dev.off()


data=TME.obj[[c('project','cell.type.level2')]]
#data<-subset(data,cell.type.level2!='Erythrocytes')

data1=as.data.frame(table(data))
data2 = data1 %>% group_by(project) %>% summarise(count=sum(Freq))
data3<-inner_join(data1,data2,by='project')
data3$percentage<- data3$Freq*100/data3$count;


plotx1<-ggplot(data3, aes(x=project,fill=cell.type.level2,y=percentage)) +
  geom_bar(stat="identity",position = 'stack') +
  scale_fill_manual(values = cell.type.level2.color) +
  #scale_x_continuous(breaks=c(0,25,50,75,100)) +
  labs(x='',y = '% (of TME cells)') +
  #facet_grid(.~Diagnosis,scales = 'free',space = 'free') +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(angle = 90,size = 5,vjust=0.2,hjust=0.95),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        legend.title=element_text(size=0), 
        legend.text=element_text(size=12),
        legend.key = element_rect(colour = NA, fill = NA),
        strip.text.x = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5,size=20))

plotx2<-ggplot(subset(data3,cell.type.level2=='Progenitors'), aes(x=project,y=count,fill=project)) +
  geom_bar(stat="identity",position = 'stack') +
  scale_fill_manual(values=c('GSE161801'='#beaed4','GSE223060'='#fdc086','GSE271107'='#7fc97f','GSE124310'='#386cb0')) +
  #scale_x_continuous(breaks=c(0,25,50,75,100)) +
  labs(x='',y = '% (of TME cells)') +
  #facet_grid(.~Diagnosis,scales = 'free',space = 'free') +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(angle = 90,size = 5,vjust=0.2,hjust=0.95),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        legend.title=element_text(size=0), 
        legend.text=element_text(size=12),
        legend.key = element_rect(colour = NA, fill = NA),
        strip.text.x = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5,size=20))

plotx<-plot_grid(
  plotx2, plotx1,
  align = "vbht", axis = "rb",
  nrow = 2, rel_heights = c(1,2)
)

ggsave('TME.5datasets.composition.stackbarplot.pdf',useDingbats=F,plotx,width = 6,height=10)

###################### Figure S8D-E NMF curve, heatmap ######################

meta<- readRDS("TME.4datasets.meta.data.rds") %>%
  add_count(SampleID, name = "n_sample") %>%
  add_count(cell.status, name = "n_status") %>%
  filter(
    n_sample >= 200,
    n_status >= 50
  ) %>%
  select(-n_sample, -n_status) %>%
  count(SampleID, cell.status, name = "Freq") %>%
  complete(SampleID, cell.status, fill = list(Freq = 0)) %>%
  group_by(SampleID) %>%
  mutate(
    count = sum(Freq),
    percentage = Freq / count
  )

ratio <- meta %>%
  dplyr::select(SampleID,cell.status,percentage) %>%
  tidyr::pivot_wider(
    names_from = cell.status,
    values_from = percentage,
    values_fill = 0
  ) %>%
  tibble::column_to_rownames("SampleID") %>%
  as.matrix()

scale_ratio <- apply(ratio, MARGIN = 2, function(x) (x-min(x))/(max(x)-min(x)))
scale_ratio <- as.data.frame(scale_ratio) %>% t()

ranks <- 2:10
estim.coad <- nmf(scale_ratio, ranks, nrun=500,method = "lee")
#saveRDS(estim.coad,'TME.4datasets.estim.coad.rds')

pdf('TME.4datasets.NMF_rank_survey.pdf',height=6,width=10)
plotx<-plot(estim.coad)
print(plotx)
dev.off()

# seed = 123
# nmf.rank <- nmf(scale_ratio, rank=5, nrun = 100, seed = seed)

nmf.rank<-readRDS('TME.4datasets.NMF.res_k5.rds')

#z_score
z_ratio <- scale(ratio) %>% as.data.frame() %>% t()

anno<-readRDS('4public_datasets.annotation.rds')

rownames(anno)<-anno$SampleID
anno$SampleID<-NULL

anno_color<-list('Diagnosis'=c('nBM'='#4dac26','MGUS'='#1c75bc','SMM'='#fbb040','NDMM'='#d01c8b', 'RRMM'='#8856a7'),
                 'project'=c('GSE161801'='#beaed4','GSE223060'='#fdc086','GSE271107'='#ffff99','GSE124310'='#386cb0'),
                 'SampleType'=c('CD138_neg'='#b2df8a',
                                'WBM'='#fb9a99',
                                'CD138_neg/WBM Mix'='#cab2d6'),
                 'Ecotype'=c('E1'='#a6cee3','E2'='#b2df8a','E3'='#fb9a99','E4'='#fdbf6f', 'E5'='#cab2d6'),
                 'cohort'=c('GSE223060_WashU_cohort1'='#7fc97f',
                            'GSE223060_WashU_cohort2'='#e31a1c',
                            'GSE223060_MMRF'='#a65628',
                            'GSE161801'='#beaed4','GSE223060'='#fdc086','GSE271107'='#ffff99','GSE124310'='#386cb0'))


group_map <- c(`1` = 3, `2` = 1, `3` = 5, `4` = 4, `5` = 2)

B<-basis(nmf.rank)

sample_group <- data.frame(
  SampleID = names(predict(nmf.rank, what = "samples")),
  group = predict(nmf.rank, what = "samples")
) %>%
  mutate(
    group = factor(unname(group_map[as.character(group)]), levels = 1:5)
)
sample_group$diagnosis<-anno[rownames(sample_group),'Diagnosis']

sample_order <- with(sample_group, order(group, diagnosis))

cell_group <- data.frame(
  SampleID = names(predict(nmf.rank, what = "features")),
  group = predict(nmf.rank, what = "features"),
  value = apply(B, 1, max)
) %>%
  mutate(
    group = factor(unname(group_map[as.character(group)]), levels = 1:5)
  )

cell_order <- with(cell_group, order(group, -value))

z_ratio_2<-z_ratio[cell_group$SampleID[cell_order],sample_group$SampleID[sample_order]]

cell_group_split<-cell_group$group[cell_order]
sample_group_split<-sample_group$group[sample_order]

top_ann = HeatmapAnnotation( df=anno[colnames(z_ratio_2),],
                             which='column',
                             annotation_name_side = "right",
                             col=anno_color,
                             show_legend = T
)

plot_matrix<-z_ratio_2/4
plot_matrix[plot_matrix>0.5]<- 0.5
plot_matrix[plot_matrix< -0.5]<- -0.5


pdf('4datasets_NMF_t12_k5.pdf',height=15,width=25)
ht_list<-Heatmap(plot_matrix, name = "ratio",
                 top_annotation = top_ann,
                 row_split = cell_group_split,
                 column_split = sample_group_split,
                 row_gap = unit(2, "mm"),
                 column_gap = unit(2, "mm"),
                 cluster_rows = F,
                 cluster_columns = F,
                 #column_order = order(factor(info.matrix$group, levels = c(1:8))),
                 row_names_gp = grid::gpar(fontsize = 10),
                 column_names_gp = grid::gpar(fontsize = 5))
draw(ht_list, annotation_legend_side = "right")
dev.off()


