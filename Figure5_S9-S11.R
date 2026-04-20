source('global_config.R')

library('CellChat')

###################### Figure 5A, S9B cellchat heatmap ######################

cellchat_subset.list<-readRDS('cellchat_subset.list.by.Ecotype.rds')

cellchat <- mergeCellChat(cellchat_subset.list, add.names = paste0('E',1:5))

group.cellType <- c('B', 
                    'T/NK', 'T/NK', 'T/NK', 'T/NK', 'T/NK', 
                    'Tumor', 
                    'Progenitors', 
                    'Myeloid', 'Myeloid', 'Myeloid', 'Myeloid', 'Myeloid', 'Myeloid', 'Myeloid') # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- c('Mature_B', 'CD8T', 'CD4T', 'NK', 'gdT', 'MAIT', 
                           'cPC', 'Progenitors', 'CD14_Mono', 'CD14CD16_Mono', 'CD16_Mono', 'Macrophage', 'Neutrophil', 'cDC', 'pDC')


pdf('CXCL.in.E1.pdf',height = 5,width = 5)
netVisual_aggregate(cellchat_subset.list[[1]], signaling = c('CXCL'), layout = "chord",group = group.cellType,
                    color.use=c('Mature_B'='#1f78b4',
                                'CD8T'='#33a02c','CD4T'='#936BD8','NK'='#a6cee3','gdT'='#FFA8BB','MAIT'='#84DF67',
                                'Stroma'='#b15928', 'cPC'='#525252', 'Platelet'='#6a3d9a', 'Progenitors'='#f4a900',
                                'CD14_Mono'='#ff7f00','CD14CD16_Mono'='#b2df8a','CD16_Mono'='#7D6294','Macrophage'='#E1D1AF','Neutrophil'='#e31a1c','cDC'='#D947E4',
                                'pDC'='#dc9ae9'))
netVisual_chord_gene(cellchat_subset.list[[1]], signaling = c("CXCL"),legend.pos.x = 8,
                     color.use=c('Mature_B'='#1f78b4',
                                 'CD8T'='#33a02c','CD4T'='#936BD8','NK'='#a6cee3','gdT'='#FFA8BB','MAIT'='#84DF67',
                                 'Stroma'='#b15928', 'cPC'='#525252', 'Platelet'='#6a3d9a', 'Progenitors'='#f4a900',
                                 'CD14_Mono'='#ff7f00','CD14CD16_Mono'='#b2df8a','CD16_Mono'='#7D6294','Macrophage'='#E1D1AF','Neutrophil'='#e31a1c','cDC'='#D947E4',
                                 'pDC'='#dc9ae9'))

dev.off()


pdf('CCL.in.E3.pdf',height = 5,width = 5)
netVisual_aggregate(cellchat_subset.list[[3]], signaling = c('CCL'), layout = "chord",group = group.cellType,
                    color.use=c('Mature_B'='#1f78b4',
                                'CD8T'='#33a02c','CD4T'='#936BD8','NK'='#a6cee3','gdT'='#FFA8BB','MAIT'='#84DF67',
                                'Stroma'='#b15928', 'cPC'='#525252', 'Platelet'='#6a3d9a', 'Progenitors'='#f4a900',
                                'CD14_Mono'='#ff7f00','CD14CD16_Mono'='#b2df8a','CD16_Mono'='#7D6294','Macrophage'='#E1D1AF','Neutrophil'='#e31a1c','cDC'='#D947E4',
                                'pDC'='#dc9ae9'))
netVisual_chord_gene(cellchat_subset.list[[3]], signaling = c("CCL"),legend.pos.x = 8,
                     color.use=c('Mature_B'='#1f78b4',
                                 'CD8T'='#33a02c','CD4T'='#936BD8','NK'='#a6cee3','gdT'='#FFA8BB','MAIT'='#84DF67',
                                 'Stroma'='#b15928', 'cPC'='#525252', 'Platelet'='#6a3d9a', 'Progenitors'='#f4a900',
                                 'CD14_Mono'='#ff7f00','CD14CD16_Mono'='#b2df8a','CD16_Mono'='#7D6294','Macrophage'='#E1D1AF','Neutrophil'='#e31a1c','cDC'='#D947E4',
                                 'pDC'='#dc9ae9'))
dev.off()

pdf('IL1.in.E4.pdf',height = 5,width = 5)
netVisual_aggregate(cellchat_subset.list[[4]], signaling = c('IL1'), layout = "chord",group = group.cellType,
                    color.use=c('Mature_B'='#1f78b4',
                                'CD8T'='#33a02c','CD4T'='#936BD8','NK'='#a6cee3','gdT'='#FFA8BB','MAIT'='#84DF67',
                                'Stroma'='#b15928', 'cPC'='#525252', 'Platelet'='#6a3d9a', 'Progenitors'='#f4a900',
                                'CD14_Mono'='#ff7f00','CD14CD16_Mono'='#b2df8a','CD16_Mono'='#7D6294','Macrophage'='#E1D1AF','Neutrophil'='#e31a1c','cDC'='#D947E4',
                                'pDC'='#dc9ae9'))
netVisual_chord_gene(cellchat_subset.list[[4]], signaling = c("IL1"),legend.pos.x = 8,
                     color.use=c('Mature_B'='#1f78b4',
                                 'CD8T'='#33a02c','CD4T'='#936BD8','NK'='#a6cee3','gdT'='#FFA8BB','MAIT'='#84DF67',
                                 'Stroma'='#b15928', 'cPC'='#525252', 'Platelet'='#6a3d9a', 'Progenitors'='#f4a900',
                                 'CD14_Mono'='#ff7f00','CD14CD16_Mono'='#b2df8a','CD16_Mono'='#7D6294','Macrophage'='#E1D1AF','Neutrophil'='#e31a1c','cDC'='#D947E4',
                                 'pDC'='#dc9ae9'))
dev.off()

pdf('IFN-II.in.E5.pdf',height = 5,width = 5)
netVisual_aggregate(cellchat_subset.list[[5]], signaling = c('IFN-II'), layout = "chord",group = group.cellType,
                    color.use=c('Mature_B'='#1f78b4',
                                'CD8T'='#33a02c','CD4T'='#936BD8','NK'='#a6cee3','gdT'='#FFA8BB','MAIT'='#84DF67',
                                'Stroma'='#b15928', 'cPC'='#525252', 'Platelet'='#6a3d9a', 'Progenitors'='#f4a900',
                                'CD14_Mono'='#ff7f00','CD14CD16_Mono'='#b2df8a','CD16_Mono'='#7D6294','Macrophage'='#E1D1AF','Neutrophil'='#e31a1c','cDC'='#D947E4',
                                'pDC'='#dc9ae9'))
netVisual_chord_gene(cellchat_subset.list[[5]], signaling = c("IFN-II"),legend.pos.x = 8,
                     color.use=c('Mature_B'='#1f78b4',
                                 'CD8T'='#33a02c','CD4T'='#936BD8','NK'='#a6cee3','gdT'='#FFA8BB','MAIT'='#84DF67',
                                 'Stroma'='#b15928', 'cPC'='#525252', 'Platelet'='#6a3d9a', 'Progenitors'='#f4a900',
                                 'CD14_Mono'='#ff7f00','CD14CD16_Mono'='#b2df8a','CD16_Mono'='#7D6294','Macrophage'='#E1D1AF','Neutrophil'='#e31a1c','cDC'='#D947E4',
                                 'pDC'='#dc9ae9'))
dev.off()

pdf('TGFb.in.E5.pdf',height = 5,width = 5)
netVisual_aggregate(cellchat_subset.list[[5]], signaling = c('TGFb'), layout = "chord",group = group.cellType,
                    color.use=c('Mature_B'='#1f78b4',
                                'CD8T'='#33a02c','CD4T'='#936BD8','NK'='#a6cee3','gdT'='#FFA8BB','MAIT'='#84DF67',
                                'Stroma'='#b15928', 'cPC'='#525252', 'Platelet'='#6a3d9a', 'Progenitors'='#f4a900',
                                'CD14_Mono'='#ff7f00','CD14CD16_Mono'='#b2df8a','CD16_Mono'='#7D6294','Macrophage'='#E1D1AF','Neutrophil'='#e31a1c','cDC'='#D947E4',
                                'pDC'='#dc9ae9'))
netVisual_chord_gene(cellchat_subset.list[[5]], signaling = c("TGFb"),legend.pos.x = 8,
                     color.use=c('Mature_B'='#1f78b4',
                                 'CD8T'='#33a02c','CD4T'='#936BD8','NK'='#a6cee3','gdT'='#FFA8BB','MAIT'='#84DF67',
                                 'Stroma'='#b15928', 'cPC'='#525252', 'Platelet'='#6a3d9a', 'Progenitors'='#f4a900',
                                 'CD14_Mono'='#ff7f00','CD14CD16_Mono'='#b2df8a','CD16_Mono'='#7D6294','Macrophage'='#E1D1AF','Neutrophil'='#e31a1c','cDC'='#D947E4',
                                 'pDC'='#dc9ae9'))
dev.off()

pdf('TNF.in.E5.pdf',height = 5,width = 5)
netVisual_aggregate(cellchat_subset.list[[5]], signaling = c('TNF'), layout = "chord",group = group.cellType,
                    color.use=c('Mature_B'='#1f78b4',
                                'CD8T'='#33a02c','CD4T'='#936BD8','NK'='#a6cee3','gdT'='#FFA8BB','MAIT'='#84DF67',
                                'Stroma'='#b15928', 'cPC'='#525252', 'Platelet'='#6a3d9a', 'Progenitors'='#f4a900',
                                'CD14_Mono'='#ff7f00','CD14CD16_Mono'='#b2df8a','CD16_Mono'='#7D6294','Macrophage'='#E1D1AF','Neutrophil'='#e31a1c','cDC'='#D947E4',
                                'pDC'='#dc9ae9'))
netVisual_chord_gene(cellchat_subset.list[[5]], signaling = c("TNF"),legend.pos.x = 8,
                     color.use=c('Mature_B'='#1f78b4',
                                 'CD8T'='#33a02c','CD4T'='#936BD8','NK'='#a6cee3','gdT'='#FFA8BB','MAIT'='#84DF67',
                                 'Stroma'='#b15928', 'cPC'='#525252', 'Platelet'='#6a3d9a', 'Progenitors'='#f4a900',
                                 'CD14_Mono'='#ff7f00','CD14CD16_Mono'='#b2df8a','CD16_Mono'='#7D6294','Macrophage'='#E1D1AF','Neutrophil'='#e31a1c','cDC'='#D947E4',
                                 'pDC'='#dc9ae9'))
dev.off()

###################### Figure S9A cellchat heatmap ######################

cellchat_subset.list<-readRDS('cellchat_subset.list.by.Ecotype.rds')

cellchat <- mergeCellChat(cellchat_subset.list, add.names = paste0('E',1:5))

aa<-cellchat@LR

LR_pairs<-unique(c(aa$E1$LRsig$interaction_name,
                   aa$E2$LRsig$interaction_name,
                   aa$E3$LRsig$interaction_name,
                   aa$E4$LRsig$interaction_name,
                   aa$E5$LRsig$interaction_name))

gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", #pairLR=LR_pairs,
               slot.name = "netP",
               stacked = T, 
               color.use = c('#a6cee3','#b2df8a','#fb9a99','#fdbf6f', '#cab2d6'),
               return.data=T,
               do.stat = TRUE,comparison = c(1:5))

data_plot <- gg1$signaling.contribution %>%
  reshape2::dcast(name ~ group, value.var = "contribution") %>%
  tibble::column_to_rownames("name") %>%
  as.matrix() %>%
  t() %>%
  scale() %>%
  t()

data_plot<-data_plot[,paste0('E',1:5)]

highlight <- c(
  # E1 (Progenitor-enriched)
  "ICAM", "LIGHT", "CXCL", "IGF", "IGFBP", "PECAM1", "PECAM2", "VCAM", "CEACAM","CD39",
  "CD200", "BTLA",
  # E2 (Naïve/Memory/Treg-enriched)
  "MK","SELL",
  # E3 (Effector T cell-enriched)
  "JAM", "COMPLEMENT","EGF","ncWNT",
  # Lows to note in caption
  "ICOS", "TNF", "IFN-II",
  # E4 (Monocyte-enriched)
  "BAFF", "APRIL", "IL1", "CCL", "NECTIN",
  # Lows to note in caption
  "CD86", "TRAIL",
  # Myeloid checkpoints
  "SIRP", "VISTA",
  # E5 (Stressed/Neutrophil/Myeloma-influenced)
  "THBS", "PLAU","TGFb","CD96",
  # Lows to note in caption
  "MHC-I"
)

idx <- which(rownames(data_plot) %in% highlight)

ha <- rowAnnotation(
  mark = anno_mark(at = idx, labels = rownames(data_plot)[idx],
                   labels_gp = gpar(fontsize = 9))
)

pdf('scaled.contribution.between.ecotypes.heatmap.pdf',height = 6,width = 5)
ComplexHeatmap::Heatmap(data_plot, 
                        name="contribution",
                        row_names_gp = gpar(fontsize = 9),
                        show_column_names = T,
                        show_row_names = F,
                        column_names_gp = gpar(fontsize = 9),
                        right_annotation = ha,
                        #bottom_annotation = regulon_ann_v,
                        cluster_columns = F,
                        cluster_rows = T)
dev.off()

###################### Figure S9C CCL, CXCL, IL1, IFNG, TGFB, TNF ligand/receptor expression in each ecotype ######################

TME.seu.obj<-readRDS('MM.TME.235samples.TME.no_Eryth.filter.scs.rds')

ligand<-c('CXCL12','CXCL8','CXCL16',
          'CCL3','CCL5',
          'IL1B','IL18',
          'IFNG',
          'TGFB1',
          'TNF')
receptor<-c('CXCR4','CXCR2','CXCR6',
            'CCR1',
            'IL1R2','IL18R1','IL18RAP',
            'IFNGR1','IFNGR2',
            'TGFBR1','TGFBR2','ACVR1B',
            'TNFRSF1A','TNFRSF1B')

p<-DotPlot(TME.seu.obj, features = ligand,col.max = 1,col.min=-1,group.by = 'cell.status') 
data<-p$data[,c('id','features.plot','pct.exp','avg.exp.scaled')]

row_clusters<-readRDS('NMF.k5.cell_group.rds')

data<-left_join(data,row_clusters[,c('SampleID','group')],by=c('id'='SampleID'))
data<-data[!is.na(data$group),]
#data$pct.exp[data$pct.exp>60]<-60

cell.use.list<-readRDS('cell.use.in.S9C.list.rds')

data<-data[data$id%in%unlist(cell.use.list),]

plotx<-ggplot(data, aes(y = features.plot,x = id)) +
  geom_point(aes(fill = avg.exp.scaled,size =pct.exp),color='black',shape=21) +
  scale_fill_gradientn(colours=c('#f7fcfd','#f7fcfd','#e0ecf4','#bfd3e6','#9ebcda','#8c96c6','#8c6bb1','#88419d','#810f7c','#4d004b')) +
  scale_size(range = c(0,6.5), limits = c(0, 100), breaks = c(0,20,40,60,80,100)) +
  facet_grid(cols=vars(group),scales='free',space='free')+
  labs(x='',y='')+
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "grey90",size=0.2),
        strip.background = element_blank(),strip.text = element_text(angle = 90),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle = 90,vjust=0.2,hjust=0.95) )

ggsave('Ligand.expression.in.each.ecotype.pdf',useDingbats=F,plotx,width = 15,height=4.5)


p<-DotPlot(TME.seu.obj, features = receptor,col.max = 1,col.min=-1,group.by = 'cell.status') 
data<-p$data[,c('id','features.plot','pct.exp','avg.exp.scaled')]

data<-left_join(data,row_clusters[,c('SampleID','group')],by=c('id'='SampleID'))
data<-data[!is.na(data$group),]

data<-data[data$id%in%unlist(cell.use.list),]

plotx<-ggplot(data, aes(y = features.plot,x = id)) + 
  geom_point(aes(fill = avg.exp.scaled,size =pct.exp),color='black',shape=21) +
  scale_fill_gradientn(colours=c('#f7fcfd','#f7fcfd','#e0ecf4','#bfd3e6','#9ebcda','#8c96c6','#8c6bb1','#88419d','#810f7c','#4d004b')) +
  scale_size(range = c(0,6.5), limits = c(0, 100), breaks = c(0,20,40,60,80,100)) + 
  facet_grid(cols=vars(group),scales='free',space='free')+
  labs(x='',y='')+
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "grey90",size=0.2),
        strip.background = element_blank(),strip.text = element_text(angle = 90),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle = 90,vjust=0.2,hjust=0.95) )

ggsave('Receptor.expression.in.each.ecotype.pdf',useDingbats=F,plotx,width = 15,height=5.5)

###################### Figure 5B, S10A IREA heatmap ######################

enrichmentscore<-c()

for( i in c('CD8T','CD4T','Macrophage','NK','CD14_Mono','CD16_Mono','Neutrophil','Bcell','cDC2','cDC1','pDC','LAMP3_DC')){
  
  data<-read.table(paste0(i,'.immune.dictionary.txt'),header = T,sep = '\t')
  enrichmentscore<-rbind(enrichmentscore,data)
  
}

saveRDS(enrichmentscore,'merge.enrichmentscore.immune.dictionary.rds')

enrichmentscore<-readRDS('merge.enrichmentscore.immune.dictionary.rds')

enrichmentscore$logfdr[enrichmentscore$enrichmentscore<0]<- -enrichmentscore$logfdr[enrichmentscore$enrichmentscore<0]

aa<-reshape2::dcast(enrichmentscore,cytokine~cell.status,value.var = 'logfdr')
rownames(aa)<-aa$cytokine
aa$cytokine<-NULL

aa[aa > -log10(0.05)]<-3
aa[aa < -log10(0.05)]<- 0

aa<-aa[apply(aa,1,function(x){sum(x> -log10(0.05))})>=1,]
aa<-aa[,apply(aa,2,function(x){sum(x> -log10(0.05))})>=1]

aa<-t(aa)

row_clusters<-readRDS('NMF.k5.cell_group.rds')

anno.use<-row_clusters[rownames(aa),]

right_ann <- rowAnnotation(
  group = anno_simple(anno.use[,c('group'),drop=F], 
                      col = c('1'='#a6cee3','2'='#b2df8a','3'='#fb9a99','4'='#fdbf6f', '5'='#cab2d6')),
  bar = anno_barplot(rowSums(aa>0), 
                     border = FALSE, 
                     gp = gpar(fill = "steelblue")),
  show_annotation_name = TRUE,
  annotation_name_side = "top"
)

bar_data <- colSums(aa>0)
top_ann <- HeatmapAnnotation(
  bar = anno_barplot(bar_data,
                     which = "column", 
                     gp = gpar(fill = "steelblue"), 
                     border = FALSE),
  annotation_name_side = "right"
)


pdf('merge.IREA.logfdr.heatmap.pdf',height=10,width=15)
ht_list<-Heatmap(aa, name = "signed log10p",
                 right_annotation = right_ann,
                 top_annotation = top_ann,
                 col=colorRampPalette(c("white", "firebrick3"))(100),
                 #row_split = cell.group,
                 #column_split = sample.group,
                 row_gap = unit(2, "mm"),
                 column_gap = unit(2, "mm"),
                 cluster_rows = T,
                 cluster_columns = T,
                 #column_order = order(factor(info.matrix$group, levels = c(1:8))),
                 row_names_gp = grid::gpar(fontsize = 10),
                 column_names_gp = grid::gpar(fontsize = 5))
draw(ht_list, annotation_legend_side = "right")
dev.off()


cell.use<-readRDS('cell.use.in.S9C.list.rds')
aa<-reshape2::dcast(enrichmentscore,cytokine~cell.status,value.var = 'logfdr')
rownames(aa)<-aa$cytokine
aa$cytokine<-NULL


pdf('IREA.logfdr.heatmap.pdf')

for(i in 1:5){
  aa.use<-aa[,intersect(colnames(aa),cell.use[[i]])]
  aa.use<-aa.use[apply(aa.use,1,function(x){sum(x> -log10(0.05))})>=2,]
  aa.use<-aa.use[,apply(aa.use,2,function(x){sum(x> -log10(0.05))})>=1]
  
  aa.use[aa.use>3]<-3
  aa.use[aa.use< -3]<- -3
  
  plotx<-pheatmap::pheatmap(aa.use[names(sort(apply(aa.use,1,function(x){sum(x> -log10(0.05))}),decreasing = T)),],#main='Jaccard index',x
                            color=colorRampPalette(c("navy", "white", "firebrick3"))(100),#colorRampPalette(c('blue','white','red'))(100),
                            scale = 'none',
                            show_colnames =T,
                            #angle_col = 0,
                            show_rownames =T,
                            cluster_rows = F,
                            cluster_cols = T
  )
  
  print(plotx)
}
dev.off()

###################### Figure S10B SCENIC Jaccard index heatmap ######################

library('SCopeLoomR')
library('SCENIC')

# data<-c()
# rss_all<-c()
# for(celltype in c('CD4T','CD8T','NK','Bcell','CD14_Mono','CD16_Mono','Macrophage','Neutrophil','DC')){ # 'CD4T',
#   
#   loom<-open_loom(paste0('../SCENIC/',celltype,'.filter.SCENIC.loom'))
#   
#   cellInfo <- get_cell_annotation(loom)
#   regulonsAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
#   
#   
#   regulonsAUC <- regulonsAUC[onlyNonDuplicatedExtended(rownames(regulonsAUC)),]
#   
#   rss <- calcRSS(AUC=getAUC(regulonsAUC), cellAnnotation=cellInfo[colnames(regulonsAUC), 'CellType'])
#   
#   pdf(paste0(celltype,'.SCENIC.plotRSS.pdf'),height = 5,width = 5.5)
#   for(i in unique(cellInfo$CellType)){
#     p<-plotRSS_oneSet(rss, setName = i, n=5)
#     print(p)
#   }
#   dev.off()
#   
#   #rss_all<-bind_cols(rss_all,rss)
#   rss_all<-merge(rss_all, rss, by = "row.names", all = TRUE)
#   rownames(rss_all) <- rss_all$Row.names
#   rss_all <- rss_all[, -1]
#   
# }
# 
# saveRDS(rss_all,'cell.status.SCENIC.rss.rds')
# 
# rss_all$TF<-rownames(rss_all)
# data<-melt(rss_all)
# data<-data[!is.na(data$value),]
# colnames(data)<-c('TF','cell.type','rss')
# data<-data[,c('cell.type','TF','rss')]
# data <- data %>% group_by(cell.type) %>% top_n(n = 50, wt = rss)
# 
# saveRDS(data,'cell.status.SCENIC.top50.TFs.rds')

# data1<-readRDS('cell.status.SCENIC.top50.TFs.rds')
# data<-left_join(data,data1[,c('cell.type','TF','cell.type.level2','group')],by=c('cell.type'='cell.type','TF'='TF'))
# data<-as.data.frame(data)

data<-readRDS('cell.status.SCENIC.top50.TFs.rds')  

celltype<-unique(data$cell.type)

jaccard_index<-matrix(0,length(celltype),length(celltype))

for(i in 1:nrow(jaccard_index)){
  for(j in 1:ncol(jaccard_index)){
    
    I<-intersect(data$TF[data$cell.type==celltype[i]],data$TF[data$cell.type==celltype[j]])
    U<-union(data$TF[data$cell.type==celltype[i]],data$TF[data$cell.type==celltype[j]])
    jaccard_index[i,j]<-length(I)/50 #length(U)
    
  }
}

rownames(jaccard_index)<-celltype
colnames(jaccard_index)<-celltype

jaccard_index[jaccard_index==1]<-NA

jaccard_index[jaccard_index>0.5]<-0.5


row_clusters<-readRDS('NMF.k5.cell_group.rds')

anno.use<-row_clusters[rownames(jaccard_index),]
top_ann = HeatmapAnnotation( df=anno.use[,c('group'),drop=F],
                             which='column',
                             annotation_name_side = "right",
                             col=list('group'=c('1'='#a6cee3','2'='#b2df8a','3'='#fb9a99','4'='#fdbf6f', '5'='#cab2d6')),
                             show_legend = T
)


pdf('SCENIC.JaccardIndex.pdf',height=10,width=10)
ht_list<-Heatmap(jaccard_index, name = "Jaccard index",
                 top_annotation = top_ann,
                 col=colorRampPalette(c('#f7fcfd','#f7fcfd','#f7fcfd','#f7fcfd','#e0ecf4','#bfd3e6','#9ebcda','#8c96c6','#8c6bb1','#88419d','#810f7c','#4d004b'))(100),
                 #row_split = cell.group,
                 #column_split = sample.group,
                 row_gap = unit(2, "mm"),
                 column_gap = unit(2, "mm"),
                 cluster_rows = T,
                 cluster_columns = T,
                 #column_order = order(factor(info.matrix$group, levels = c(1:8))),
                 row_names_gp = grid::gpar(fontsize = 10),
                 column_names_gp = grid::gpar(fontsize = 5))
draw(ht_list, annotation_legend_side = "right")
dev.off()

###################### Figure 5C, 5D, S10C, S10D SCENIC Jaccard index heatmap by ecotype ######################
data<-readRDS('cell.status.SCENIC.top50.TFs.rds')  

celltype<-unique(data$cell.type)

jaccard_index<-matrix(0,length(celltype),length(celltype))

for(i in 1:nrow(jaccard_index)){
  for(j in 1:ncol(jaccard_index)){
    
    I<-intersect(data$TF[data$cell.type==celltype[i]],data$TF[data$cell.type==celltype[j]])
    U<-union(data$TF[data$cell.type==celltype[i]],data$TF[data$cell.type==celltype[j]])
    jaccard_index[i,j]<-length(I)/50 #length(U)
    
  }
}

rownames(jaccard_index)<-celltype
colnames(jaccard_index)<-celltype

jaccard_index[jaccard_index==1]<-NA

jaccard_index[jaccard_index>0.5]<-0.5

row_clusters<-readRDS('NMF.k5.cell_group.rds')
# data<-left_join(data,row_clusters[,c('SampleID','group')],by=c('cell.type'='SampleID'))
row_clusters.use<-subset(row_clusters,SampleID%in%unique(data$cell.type))


cell.use.list<-readRDS('cell.use.in.S9C.list.rds')


pdf('ecotype.overlapped.TF.pdf',width=10)
for(i in c(1:5)){
  
  cell.use<-cell.use.list[[i]]
  cell.use<-cell.use[cell.use%in%rownames(jaccard_index)]
  
  jaccard_index.use<-jaccard_index[cell.use,cell.use]
  plotx<-pheatmap::pheatmap(jaccard_index.use,color = colorRampPalette(c('#f7fcfd','#e0ecf4','#bfd3e6','#9ebcda','#8c96c6','#8c6bb1','#88419d','#810f7c','#4d004b'))(100))
  print(plotx)
  
  aa<-data[data$cell.type%in%cell.use,]
  aa<-aa[aa$TF%in%names(table(aa$TF))[table(aa$TF)>2],]
  
  
  gene_counts <- aa %>%
    group_by(TF, cell.type.level2) %>%
    summarise(count = n(), .groups = "drop")
  
  gene_counts<-gene_counts[gene_counts$TF%in%names(table(gene_counts$TF))[table(gene_counts$TF)>1],] # only keep TFs show up in more than 1 cell type
  aa<-aa[aa$TF%in%gene_counts$TF,]
  
  bb<-sort(table(aa$TF),decreasing = T)
  
  gene_counts<-gene_counts[gene_counts$TF%in%names(bb),]
  gene_counts$TF<-factor(gene_counts$TF,levels=names(bb))
  # Plot
  plotx<-ggplot(gene_counts, aes(x = TF, y = count, fill = cell.type.level2)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c('CD4T'='#936BD8','CD8T'='#33a02c','NK'='#a6cee3','Bcell'='#1f78b4',
                                 'CD14_Mono'='#ff7f00','CD14CD16_Mono'='#b2df8a','CD16_Mono'='#7D6294','Mac'='#E1D1AF','cDC'='#D947E4',
                                 'Neutro'='#e31a1c','pDC'='#dc9ae9') )+
    theme_minimal() +
    labs(x = "Gene", y = "Count", fill = "Cell Type") +
    theme(axis.text.x = element_text(size = 12,angle=90,vjust=0.2,hjust=0.5),
          strip.text.x = element_text(size=12, color="black", face="bold"),
          text = element_text(size=15),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),axis.line = element_line(colour = "black"))
  
  print(plotx)
  
}
dev.off()

###################### Figure 5E top 50 genes overlap with deg ######################
marker_list<-readRDS('MM.TME.235samples.TME.no_Eryth.Spectra.modified_marker_list.rds')

marker_list.use<-marker_list[1:12,]

subset.deg<-readRDS('cell_subset.deg.rds')

deg_all<-c()

for(celltype in names(subset.deg)){
  
  deg_top50 <- subset.deg[[celltype]] %>% filter(p_val_adj<0.01) %>% 
    group_by(cluster) %>% 
    top_n(n = 50, wt = dplyr::desc(p_val_adj)) %>% 
    top_n(n = 50, wt = avg_log2FC) %>% 
    as.data.frame()
  
  deg_all<-rbind(deg_all,deg_top50)
}

celltype<-unique(deg_all$cluster)

jaccard_index<-matrix(0,length(celltype),nrow(marker_list.use))

for(i in 1:nrow(jaccard_index)){
  for(j in 1:ncol(jaccard_index)){
    
    I<-intersect(deg_all$gene[deg_all$cluster==celltype[i]],marker_list.use[j,])
    U<-union(deg_all$gene[deg_all$cluster==celltype[i]],marker_list.use[j,])
    jaccard_index[i,j]<-length(I)/50#length(U)
    
  }
}

rownames(jaccard_index)<-celltype
colnames(jaccard_index)<-paste0('gene_module_',0:11)

jaccard_index[jaccard_index>0.5]<-0.5

plotx<-pheatmap::pheatmap(t(jaccard_index),main='Jaccard index',
                          color=colorRampPalette(c('white','red'))(100),
                          show_colnames =T,
                          #angle_col = 0,
                          show_rownames =T,
                          display_numbers=T,fontsize_number=3,
                          cluster_rows = T,
                          cluster_cols = T
)
ggsave('gene.module.overlap.with.deg.JaccardIndex.pdf',plotx,width=10,height=4)

###################### Figure 5F spectra score umap in each cell type ######################

spectra_score<-readRDS('MM.TME.235samples.TME.no_Eryth.Spectra.cell_scores.rds')

gp.list<-NULL
for(i in c('CD4T','CD8T','NK','Bcell','CD14_Mono','CD16_Mono','Macrophage','DC')){ # 
  
  seu.obj<-readRDS(paste0('MM.TME.235samples.',i,'.filter.scs.rds'))
  
  meta<-FetchData(seu.obj,vars=c('orig.ident'))
  meta$barcode<-rownames(meta)
  meta<-left_join(meta,spectra_score,by=c('barcode'='barcode'))
  
  seu.obj$GM1<-meta$gene_module_1
  seu.obj$GM3<-meta$gene_module_3
  seu.obj$GM5<-meta$gene_module_5
  
  seu.obj$GM1[seu.obj$GM1>quantile(seu.obj$GM1,0.99,na.rm=T)]<-quantile(seu.obj$GM1,0.99,na.rm=T)
  seu.obj$GM3[seu.obj$GM3>quantile(seu.obj$GM3,0.99,na.rm=T)]<-quantile(seu.obj$GM3,0.99,na.rm=T)
  seu.obj$GM5[seu.obj$GM5>quantile(seu.obj$GM5,0.99,na.rm=T)]<-quantile(seu.obj$GM5,0.99,na.rm=T)
  
  p2<-FeaturePlot(seu.obj,features = c('GM1','GM3','GM5'),
                  pt.size=1,cols=c("white", "#67000d"),order = F,ncol = 3, raster = T)&
    #scale_color_gradientn(colors=pals::jet()[c(1:11,15:25)],na.value = "#e6e6e6")&
    scale_color_gradientn(colors=c('#fff5f0','#67000d'),na.value = "#e6e6e6")&
    theme_umap()
  
  gp.list[[length(gp.list)+1]]<-p2
  
}

do.call(ggarrange,c(gp.list,ncol = 1,nrow = 1)) -> combined.gp
pdf('gene_module.umap.pdf',height=6,width = 18)
print(combined.gp)
dev.off()

###################### Figure 5G TME spectra score by ecotype heatmap ######################

spectra_score <- readRDS("MM.TME.235samples.TME.no_Eryth.Spectra.cell_scores.rds")

meta <- readRDS("MM.TME.235samples.TME.meta.data.rds") %>%
  left_join(spectra_score, by = "barcode")

sample_group <- readRDS("NMF.k5.sample_group.rds") %>%
  dplyr::select(SampleID, Ecotype)

data <- meta %>%
  dplyr::select(NewID, cell.type, gene_module_1, gene_module_3, gene_module_5) %>%
  inner_join(sample_group, by = c("NewID" = "SampleID")) %>%
  pivot_longer(
    cols = starts_with("gene_module_"),
    names_to = "variable",
    values_to = "value"
  )

data_global <- data %>%
  filter(cell.type %in% c(
    "CD4T", "CD8T", "gdT", "MAIT", "NK", "Mature_B",
    "CD14_Mono", "CD14CD16_Mono", "CD16_Mono",
    "Macrophage", "cDC", "pDC"
  ))


data_global <- data_global %>%
  group_by(cell.type, variable) %>%
  mutate(
    value = pmin(value, quantile(value, 0.95, na.rm = TRUE)),
    value = as.numeric(scale(value))
  ) %>%
  ungroup() %>%
  group_by(Ecotype, cell.type, variable) %>%
  summarise(
    scaled_mean_value = mean(value),
    .groups = "drop"
  ) %>%
  mutate(
    scaled_mean_value = pmax(pmin(scaled_mean_value, 0.5), -0.5),
    cell.type = factor(
      cell.type,
      levels = c(
        "CD4T", "CD8T", "gdT", "MAIT", "NK", "Mature_B",
        "CD14_Mono", "CD14CD16_Mono", "CD16_Mono",
        "Macrophage", "cDC", "pDC"
      )
    ),
    variable = factor(variable, levels = rev(unique(variable)))
  )


plotx<-ggplot(data_global, aes(x = Ecotype, y = variable, fill = scaled_mean_value)) +
  geom_tile(color = "black") +
  scale_fill_gradientn(colours=rev(c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061')),na.value = 'white') +
  facet_grid(cols=vars(cell.type))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+
  coord_fixed()

ggsave('global_gene_module_by_ecotype_heatmap.pdf',plotx,height=4,width=20)

###################### Figure S11A TME spectra score by cell type heatmap ######################

spectra_score<-readRDS('MM.TME.235samples.TME.no_Eryth.Spectra.cell_scores.rds')

meta <- readRDS("MM.TME.235samples.TME.meta.data.rds")

meta<-FetchData(pbmc.f.TME,vars=c('cell.type'))

meta<-left_join(meta,spectra_score,by=c('barcode'='barcode'))
meta$barcode<-NULL
meta<-meta[,c(17,20:176)]


data <- meta %>% group_by(cell.type) %>% summarise(across(everything(), \(x) mean(x, na.rm = TRUE))) %>% as.data.frame()
rownames(data)<-data$cell.type
data$cell.type<-NULL
data<-t(data)
data[data>0.01]<-0.01

plotx<-pheatmap::pheatmap(data,
                          #color=c('blue','white','yellow'),
                          show_colnames =T,
                          #display_numbers=T,
                          number_color = 'black',
                          fontsize_number = 10,
                          #angle_col = 0,
                          show_rownames =T,
                          cluster_rows = F,
                          cluster_cols = F
)

ggsave('TME.cell.type.spectra_score.heatmap.pdf',plotx,width=10,height=10)

###################### Figure S11B top 50 genes per gene module ######################
gene_score<-readRDS('MM.TME.235samples.TME.no_Eryth.Spectra.gene_scores.rds')

gene_score<-gene_score[1:12,]
aa<-reshape2::melt(gene_score)

top50 <- aa %>% group_by(module) %>% top_n(n = 50, wt = value)

gp.list<-c()

for(i in paste0('gene_module_',0:11)){
  
  data.use<-subset(top50,module==i)
  data.use<-data.use[order(data.use$value,decreasing = T),]
  data.use$variable<-factor(data.use$variable,levels=rev(data.use$variable))
  
  plotx<-ggplot(data.use, aes(x=value,y=variable)) +
    geom_col(fill='#69b3a2') +
    #scale_fill_manual(values = color1) +
    #scale_y_continuous(breaks=c(0,25,50,75,100)) +
    labs(x='',y = '',title = i) +
    theme(axis.title.x = element_text(size = 16),
          axis.text.x = element_text(angle = 90,size = 15,vjust=0.2,hjust=0.95),
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          legend.title=element_text(size=0), 
          legend.text=element_text(size=12),
          legend.key = element_rect(colour = NA, fill = NA),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5,size=20))
  
  gp.list[[length(gp.list)+1]]<-plotx
  
}

do.call(ggarrange,c(gp.list,ncol = 6,nrow = 2)) -> combined.gp
pdf('gene_module.top50.gene_score.pdf',height=20,width = 20)
print(combined.gp)
dev.off()


###################### Figure S11C gene module top 50 genes GO ######################

marker_list<-readRDS('MM.TME.235samples.TME.no_Eryth.Spectra.modified_marker_list.rds')

marker_list.use<-marker_list[1:12,]

library('clusterProfiler')
library('org.Hs.eg.db')
library('GOSemSim')
library('enrichplot')

pdf('compGO.CC.emapplot.pdf',height=15,width=15)

for(i in c(2,4,6)){
  
  set.seed(20250703)
  res<-enrichGO(marker_list.use[i,],
                OrgDb=org.Hs.eg.db,
                keyType = 'SYMBOL',
                pvalueCutoff  = 0.01,
                qvalueCutoff=0.05,
                pAdjustMethod = "BH", 
                ont = 'BP',
                readable = F)
  
  aa<-pairwise_termsim(res)
  
  p<-emapplot(aa)+
    scale_fill_gradientn(colors = c("#67001f","#b2182b","#d6604d","#f4a582","#fddbc7","white","#d1e5f0","#92c5de","#4393c3","#2166ac","#053061"))+
    labs(title=i-1)+
    theme(plot.title = element_text(hjust = 0.5,size=18,face='bold'))
  
  print(p)
  
}

dev.off()

###################### Figure S11D scatter plot of correlation between GM5 and GM1, GM3 ######################

library(tidyverse)

# 1. read data

df<-readRDS('spectra_GM_score.tumor_burden.per.sample.rds')

stage_colors <- c(
  "MGUS" = "#1c75bc",
  "SMM"  = "#fbb040",
  "NDMM" = "#d01c8b",
  "RRMM" = "#8856a7"
)

df <- df %>%
  filter(
    Diagnosis %in% c("MGUS", "SMM", "NDMM", "RRMM"),
    variable %in% c("gene_module_1", "gene_module_3", "gene_module_5")
  ) %>%
  mutate(
    Diagnosis = factor(Diagnosis, levels = c("MGUS", "SMM", "NDMM", "RRMM"))
  )


# 2. wide format

df2 <- df %>%
  dplyr::select(NewID, Diagnosis, variable, mean_value) %>%
  pivot_wider(names_from = variable, values_from = mean_value)


# 3. long format for plotting

plot_df <- df2 %>%
  pivot_longer(
    cols = c(gene_module_1, gene_module_3),
    names_to = "module",
    values_to = "x_value"
  ) %>%
  mutate(
    module = factor(module, levels = c("gene_module_1", "gene_module_3"))
  )


# 4. helper functions

get_cor <- function(x, y) {
  d <- stats::na.omit(data.frame(x = x, y = y))
  if (nrow(d) < 5) {
    return(list(r = NA_real_, p = NA_real_))
  }
  ct <- suppressWarnings(stats::cor.test(d$x, d$y, method = "spearman"))
  list(r = unname(ct$estimate), p = ct$p.value)
}

fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "NA", format(round(x, digits), nsmall = digits))
}

fmt_p <- function(x) {
  ifelse(is.na(x), "NA", formatC(x, format = "e", digits = 2))
}


# 5. make label table

overall_stats <- plot_df %>%
  group_by(module) %>%
  summarise(
    stat = list(get_cor(x_value, gene_module_5)),
    .groups = "drop"
  ) %>%
  mutate(
    r = purrr::map_dbl(stat, "r"),
    p = purrr::map_dbl(stat, "p"),
    line = paste0("Overall: r=", fmt_num(r), ", p=", fmt_p(p))
  ) %>%
  dplyr::select(module, line)

diag_stats <- plot_df %>%
  group_by(module, Diagnosis) %>%
  summarise(
    stat = list(get_cor(x_value, gene_module_5)),
    .groups = "drop"
  ) %>%
  mutate(
    r = purrr::map_dbl(stat, "r"),
    p = purrr::map_dbl(stat, "p"),
    line = paste0(as.character(Diagnosis), ": r=", fmt_num(r), ", p=", fmt_p(p))
  ) %>%
  dplyr::select(module, Diagnosis, line)

label_df <- bind_rows(
  overall_stats %>% mutate(order = 1),
  diag_stats %>%
    mutate(
      order = dplyr::case_when(
        Diagnosis == "MGUS" ~ 2,
        Diagnosis == "SMM"  ~ 3,
        Diagnosis == "NDMM" ~ 4,
        Diagnosis == "RRMM" ~ 5
      )
    ) %>%
    dplyr::select(module, line, order)
) %>%
  arrange(module, order)

# panel positions
panel_pos <- plot_df %>%
  group_by(module) %>%
  summarise(
    x_min = min(x_value, na.rm = TRUE),
    x_max = max(x_value, na.rm = TRUE),
    y_min = min(gene_module_5, na.rm = TRUE),
    y_max = max(gene_module_5, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    x = x_min + 0.03 * (x_max - x_min),
    y_top = y_max - 0.03 * (y_max - y_min),
    step = 0.08 * (y_max - y_min)
  )

label_df <- label_df %>%
  group_by(module) %>%
  mutate(line_id = row_number()) %>%
  ungroup() %>%
  left_join(panel_pos, by = "module") %>%
  mutate(
    y = y_top - (line_id - 1) * step
  )


# 6. plot

p <- ggplot(plot_df, aes(x = x_value, y = gene_module_5, color = Diagnosis)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
  facet_wrap(~ module, scales = "free", ncol = 2) +
  geom_text(
    data = label_df,
    aes(x = x, y = y, label = line),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.4
  ) +
  scale_color_manual(values = stage_colors) +
  theme_classic(base_size = 13) +
  labs(
    x = "Gene module score",
    y = "Gene module 5",
    color = "Diagnosis"
  ) +
  theme(
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.5, "cm")
  )

ggsave("GM5_vs_GM1_GM3_colored_by_stage_one_figure.pdf", p, width = 11, height = 5)

