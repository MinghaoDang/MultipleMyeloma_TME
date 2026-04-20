source('global_config.R')

###################### Figure 1B Clinical summary figure ######################
data<-read_excel('TableS1.xlsx',sheet = 1)%>%as.data.frame()
data$Diagnosis<-factor(data$Diagnosis,levels=c('nBM','MGUS','SMM','NDMM','RRMM'))
data$Primary_Cytogenetics<-factor(data$Primary_Cytogenetics,levels=c('HY',
                                                                     'HY;t(4;14)','t(4;14)',
                                                                     't(6;14)',
                                                                     'HY;t(11;14)','t(11;14)',
                                                                     'HY;t(14;16)','t(14;16)',
                                                                     't(14;20)',
                                                                     'Other','NA'))
data$`Other`<-ifelse(data$Primary_Cytogenetics%in%c('Other'),'Pos','Neg')
data$`NA`<-ifelse(data$Primary_Cytogenetics%in%c('NA'),'Pos','Neg')

d1<-data[,c('NewID','Hyperdiploidy','t(4;14)','t(6;14)','t(11;14)','t(14;16)','t(14;20)','Other','NA','1q gain/CKS1B Pos','monosomy 13','17p loss/p53 Del')]


rownames(d1)<-d1$NewID
d1$NewID<-NULL

Anno<-data[,c('NewID','Diagnosis','Gender','Race','Primary_Cytogenetics','R-ISS')]
rownames(Anno)<-Anno$NewID
Anno$NewID<-NULL 

annotation_col<-Anno#[pp_order,]
colnames(annotation_col)[1]<-'Diagnosis'

annotation_col<-annotation_col[order(annotation_col$Diagnosis,annotation_col$Primary_Cytogenetics),]

anno_color<-list('Diagnosis'=c('nBM'='#4dac26','MGUS'='#1c75bc','SMM'='#fbb040','NDMM'='#d01c8b', 'RRMM'='#8856a7'),
                 'Gender'=c('Female'='#d6604d','Male'='#92c5de'),
                 'Race'=c('White or Caucasian'='#fee0d2','Black or African American'='#525252','Other'='#a65628'),
                 'R-ISS'=c('I'='#fee8c8','II'='#fc8d59','III'='#7f0000'))


ha<-HeatmapAnnotation(df=annotation_col[,c('Diagnosis','Gender','Race','R-ISS')],
                      col=anno_color,
                      show_annotation_name = T,
                      annotation_height = 0.02)

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = '#F0F0F0', col = NA))
  },
  `Pos` = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h, gp = gpar(fill = "#34A02C", col = NA))
  },
  `Neg` = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h, gp = gpar(fill = 'grey80', col = NA))
  },
  `NA` = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h, gp = gpar(fill = 'grey80', col = NA))
  }
)

d2<-t(d1)[,rownames(annotation_col)]

pdf('patient_overview.pdf',width = 20,height = 4);
pp<-oncoPrint(d2,
              alter_fun = alter_fun, #col = col,
              row_names_gp = gpar(fontsize = 8),
              top_annotation= ha,
              show_pct = F,
              row_names_side = "left",
              column_split = annotation_col$Diagnosis,
              row_split = c(rep('P',8),rep('S',3)),
              show_column_names=F,
              column_title = "Summary of data"
)
#row_order<-rownames(aa)[pp@row_order]
print(pp)
dev.off();

###################### Figure 1D area plot & p value ######################

meta <- readRDS("MM.TME.235samples.TME.meta.data.rds") %>%
  dplyr::select(NewID, cell.type.level2, Diagnosis) %>%
  add_count(NewID, name = "n_cells") %>%
  filter(n_cells >= 200) %>%
  dplyr::select(-n_cells)

# diagnosis-level composition
diagnosis_comp <- meta %>%
  count(Diagnosis, cell.type.level2, name = "Freq") %>%
  group_by(Diagnosis) %>%
  mutate(
    total = sum(Freq),
    percentage = Freq * 100 / total,
    Diagnosis_level = as.numeric(Diagnosis)
  ) %>%
  ungroup()


plotx<-ggplot(diagnosis_comp, aes(x = Diagnosis_level, y = percentage, fill = cell.type.level2)) +
  geom_area() +
  scale_fill_manual(values = cell.type.level2.color) +
  theme_minimal(base_size = 14) +
  labs(y = "% among all TME cells", x = NULL, fill = NULL)
ggsave('TME.diagnosis.area.plot.pdf',height = 5,width=5,plotx)

# sample-level composition
sample_comp <- meta %>%
  count(NewID, cell.type.level2, name = "Freq") %>%
  complete(NewID, cell.type.level2, fill = list(Freq = 0)) %>%
  left_join(meta %>% distinct(NewID, Diagnosis), by = "NewID") %>%
  group_by(NewID) %>%
  mutate(
    total = sum(Freq),
    percentage = Freq * 100 / total
  ) %>%
  ungroup()

# Kruskal-Wallis test by cell type
kruskal_results <- sample_comp %>%
  group_by(cell.type.level2) %>%
  summarise(
    p_value = kruskal.test(percentage ~ Diagnosis)$p.value
  )

###################### Figure 1E, S2G TME cell type level2 compare ######################

meta <- readRDS("MM.TME.235samples.TME.meta.data.rds") %>%
  dplyr::select(NewID, cell.type.level2, Diagnosis) %>%
  add_count(NewID, name = "n_cells") %>%
  filter(n_cells >= 200) %>%
  dplyr::select(-n_cells)

sample_comp <- meta %>%
  count(NewID, cell.type.level2, name = "Freq") %>%
  complete(NewID, cell.type.level2, fill = list(Freq = 0)) %>%
  left_join(meta %>% distinct(NewID, Diagnosis), by = "NewID") %>%
  group_by(NewID) %>%
  mutate(
    total = sum(Freq),
    percentage = Freq * 100 / total
  ) %>%
  ungroup() %>%
  as.data.frame()


plotx<-signif_plot(sample_comp,
                   x='Diagnosis',
                   y='percentage',
                   split='cell.type.level2',
                   signif.cutoff=0.05,
                   ncol=17,
                   fill_by='Diagnosis',
                   fill_color=c('nBM'='#4dac26','MGUS'='#1c75bc','SMM'='#fbb040','NDMM'='#d01c8b', 'RRMM'='#8856a7'),
                   ylab='% (of TME cells)')
ggsave('TME.cell.type.level2.compare2.pdf',plotx,width=45,height=5)


###################### Figure 1F cell subset UMAP ######################

pdf(paste0('cell.subsets.umap.',Sys.Date(),'.pdf'),height=8.5,width = 8)

for(celltype in c('CD8T','CD4T','NK','Macrophage','Neutrophil','DC','CD16_Mono','CD14_Mono','Bcell', 'Progenitors')){ # 
  
  seu.obj<-readRDS(paste0('MM.TME.235samples.',celltype,'.filter.scs.rds'))
  
  plotx<-DimPlot(object = seu.obj, reduction = "umap",pt.size = 1.5,group.by = 'cell.status',label = T,label.size = 5, repel = T, raster = T)+
    scale_color_manual(values=cell.status.color)+
    theme_umap() +
    theme(legend.position="none")
  
  print(plotx)
  
  
}

dev.off()

###################### Figure S1A ALL cell type UMAP ######################

MM.seu.obj<-readRDS('MM.TME.235samples.filter.scs.rds')

p1<-DimPlot(object = MM.seu.obj, reduction = "umap",pt.size = 0.1,group.by = c('NewID','seurat_clusters','Diagnosis'),label = F,label.size = 4, raster = T)&
  theme(legend.position = "none")&
  theme_umap()

ggsave(paste0('ALL.umap.',Sys.Date(),'.pdf'),p1,height=5,width=15)

p2<-FeaturePlot(MM.seu.obj,features = c('SDC1','CD38','TNFRSF17'),
                pt.size=0.1,max.cutoff=3,cols=c("lightgrey", "red"),order = F,ncol = 3)&
  #theme(legend.position = "none")+
  theme_umap()

ggsave(paste0('ALL.SDC1_CD38_BCMA.umap.',Sys.Date(),'.pdf'),p2,height=4.75,width=15)

p3<-FeaturePlot(MM.seu.obj,features = c('BCR_clonotype_proportion'),
                pt.size=0.1,cols=c("lightgrey", "red"),order = F,ncol = 3)&
  #theme(legend.position = "none")+
  theme_umap()

ggsave(paste0('ALL.BCR_clonotype_proportion.umap.',Sys.Date(),'.pdf'),p3,height=4.75,width=15)

p4<-FeaturePlot(MM.seu.obj,features = c('TCR_clonotype_proportion'),
                pt.size=0.1,max.cutoff=0.01,cols=c("lightgrey", "red"),order = F,ncol = 3)&
  #theme(legend.position = "none")+
  theme_umap()

ggsave(paste0('ALL.TCR_clonotype_proportion.umap.',Sys.Date(),'.pdf'),p4,height=4.75,width=15)


###################### Figure S1B tumor proportion compare ######################

data<-readRDS('sample.level.info.rds')

data<-subset(data,Non_Eryth.cell.count>=200)


P<-ggplot(data, aes(x=Diagnosis, y=Tumor.pct) ) + 
  geom_boxplot(show.legend = T,width = 0.75, aes(fill=Diagnosis), outlier.shape=NA) + 
  geom_jitter(position=position_jitter(width=0.1, height=0), size=1.5) +
  ggsignif::geom_signif(test="wilcox.test", comparisons = list(c('nBM','MGUS'),c('MGUS','SMM'),c('SMM','NDMM'),c('NDMM','RRMM')),#combn(levels(factor(data$Diagnosis)),2, simplify = F),
                        map_signif_level = F,   # c("***"=0.001, "**"=0.01, "*"=0.05)
                        vjust=-0.15,
                        textsize=6,
                        size=0.75,
                        step_increase = 0.15)+
  scale_fill_manual(values=c('nBM'='#4dac26','MGUS'='#1c75bc','SMM'='#fbb040','NDMM'='#d01c8b', 'RRMM'='#8856a7'))+
  scale_y_continuous(expand = expansion(mult = c(0.15)))+
  theme(strip.text.x = element_text(size=14, color="black", face="bold"),
        axis.text.x = element_text(angle = 90,size = 14,vjust=0.2,hjust=0.95,colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        text = element_text(size=20)) +
  labs(x="", y="% (of ALL cells)") 
ggsave('Tumor.pct.pdf',P,width=5,height=6)

###################### Figure S1C, S1E TME before batch correction UMAP (this study) ######################

MM_TME.seu.obj<-readRDS('MM.TME.235samples.TME.no_Eryth.filter.scs.rds')


p1<-DimPlot(object = MM_TME.seu.obj, reduction = "umap",pt.size = 0.1,group.by = 'cell.type.level2',label = T,label.size = 6)+
  scale_color_manual(values=cell.type.level2.color)+
  theme(legend.position = "none")+
  theme_umap()

ggsave('TME.umap.before_batch_correction.pdf',useDingbats=F,p1,width = 8,height=8)

p2<-DimPlot(object = MM_TME.seu.obj, reduction = "umap",pt.size = 0.5,group.by = 'ProcolNumber',split.by = 'ProcolNumber',label = F,label.size = 6)&
  #scale_color_manual(values=cell.type.level2.color)+
  theme(legend.position = "none")+
  theme_umap()

ggsave('TME.umap.before_batch_correction.split.by.study.pdf',useDingbats=F,p2,width = 25,height=6)

p3<-DimPlot(object = MM_TME.seu.obj, reduction = "harmony_umap",pt.size = 0.5,group.by = 'ProcolNumber',split.by = 'ProcolNumber',label = F,label.size = 6)&
  #scale_color_manual(values=cell.type.level2.color)+
  theme(legend.position = "none")+
  theme_umap()

ggsave('TME.umap.after_batch_correction.split.by.study.pdf',useDingbats=F,p3,width = 25,height=6)

###################### Figure S1D, S1F LISI ######################

library('lisi')

MM_TME.seu.obj<-readRDS('MM.TME.235samples.TME.no_Eryth.filter.scs.rds')

aa=Embeddings(object = MM_TME.seu.obj, reduction = "umap");

res <- compute_lisi(aa, MM_TME.seu.obj@meta.data, c('Batch','cell.type.level2'),perplexity = 150)

plotx1<-ggplot(res, aes(x=Batch)) +
  labs(title='iLISI before batch correction')+
  geom_density()

plotx2<-ggplot(res, aes(x=cell.type.level2)) +
  labs(title='cLISI before batch correction')+
  geom_density();

bb=Embeddings(object = MM_TME.seu.obj, reduction = "harmony_umap");

res2 <- compute_lisi(bb, MM_TME.seu.obj@meta.data, c('Batch','cell.type.level2'),perplexity = 150)

plotx3<-ggplot(res2, aes(x=Batch)) +
  labs(title='iLISI after batch correction')+
  geom_density()

plotx4<-ggplot(res2, aes(x=cell.type.level2)) +
  labs(title='cLISI after batch correction')+
  geom_density();

do.call(ggarrange,c(list(plotx1,plotx2,plotx3,plotx4),ncol = 2,nrow=2)) -> combined.gp

pdf(paste0('LISI_batch_correction_evaluation.pdf'),width = 12,height=6)
print(combined.gp)
dev.off()


###################### Figure S2A TME cell type level2 DEG bubble plot ######################
MM_TME.seu.obj<-readRDS('MM.TME.235samples.TME.no_Eryth.filter.scs.rds')

load('TME.cell.type.level2.deg.RData')
top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

Progenitors<-c('PRTN3','AZU1','MPO','ELANE','SOX4')
CD4T<-c('LTB','IL7R','CD3E','TCF7','LEF1')
CD8T<-c('CD8A','CD8B','CCL5','GZMH','IL32')
gdT<-c("TRDV2","TRGV9","TRDC","NKG7","DUSP2")
MAIT<-c("KLRB1","TRAV1-2","SLC4A10","TNFAIP3",'NCR3')
NK<-c('GNLY','PRF1','GZMB','SPON2','FGFBP2')
Mature_B<-c('CD79A','MS4A1','IGHM','BANK1','TCL1A')
CD14_Mono<-c('CD14','VCAN','LYZ','FCN1','IFI30')
CD14CD16_Mono<-c('IFITM3','APOBEC3A','AIF1','COTL1','CST3')
CD16_Mono<-c('FCGR3A','LST1','CDKN1C','MS4A7','SERPINA1')
Macrophage<-c('C1QB','C1QA','APOE','C1QC','HMOX1')
cDC<-c('FCER1A','CD1C','HLA-DQA1','HLA-DRB1','HLA-DPA1')
Neutrophil<-c('LCN2','CAMP','S100A12','MMP9','S100A8')
pDC<-c('PLD4','ITM2C','LILRA4','JCHAIN','IRF7')
Platelet<-c('PPBP','NRGN','TUBB1','PF4','CAVIN2')
#Erythrocyte<-c('HBA1','HBA2','HBD','HBM')
Stroma<-c('FABP4','CXCL12','CCL14','NNMT','SELENOP')

gene<-c(Progenitors,CD4T,CD8T,gdT,MAIT,NK,Mature_B,CD14_Mono,CD14CD16_Mono,CD16_Mono,Macrophage,cDC,Neutrophil,pDC,Platelet,Stroma)

# gene<-unique(top10$gene)

p<-DotPlot(MM_TME.seu.obj, features = gene,col.max = 1.5,col.min=-1.5,group.by = 'cell.type.level2') 
data<-p$data[,c('id','features.plot','pct.exp','avg.exp.scaled')]

data$gene_function='Mature_B'
data$gene_function[which(data$features.plot%in%CD14_Mono)]='CD14_Mono'
data$gene_function[which(data$features.plot%in%CD14CD16_Mono)]='CD14CD16_Mono'
data$gene_function[which(data$features.plot%in%CD16_Mono)]='CD16_Mono'
data$gene_function[which(data$features.plot%in%CD4T)]='CD4T'
data$gene_function[which(data$features.plot%in%CD8T)]='CD8T'
data$gene_function[which(data$features.plot%in%gdT)]='gdT'
data$gene_function[which(data$features.plot%in%MAIT)]='MAIT'
data$gene_function[which(data$features.plot%in%cDC)]='cDC'
data$gene_function[which(data$features.plot%in%Macrophage)]='Macrophage'
data$gene_function[which(data$features.plot%in%Neutrophil)]='Neutrophil'
data$gene_function[which(data$features.plot%in%NK)]='NK'
data$gene_function[which(data$features.plot%in%pDC)]='pDC'
data$gene_function[which(data$features.plot%in%Platelet)]='Platelet'
data$gene_function[which(data$features.plot%in%Progenitors)]='Progenitors'
#data$gene_function[which(data$features.plot%in%Erythrocyte)]='Erythrocyte'
data$gene_function[which(data$features.plot%in%Stroma)]='Stroma'

data$gene_function<-factor(data$gene_function,levels=c('Progenitors','CD4T','CD8T','gdT','MAIT','NK','Mature_B','CD14_Mono','CD14CD16_Mono','CD16_Mono','Macrophage','cDC','Neutrophil','pDC','Platelet','Stroma'))

plotx<-ggplot(data, aes(x = features.plot,y = id)) +        ## global aes
  geom_point(aes(fill = avg.exp.scaled,size =pct.exp),color='black',shape=21)  +    ## geom_point for circle illusion
  #scale_fill_gradientn(colours=rev(color),limits=c(0,max(data$avg.exp)))+       ## color of the corresponding aes
  scale_fill_gradientn(colours=c('#f7fcfd','#f7fcfd','#e0ecf4','#bfd3e6','#9ebcda','#8c96c6','#8c6bb1','#88419d','#810f7c','#4d004b')) +
  scale_size(range = c(0,6), limits = c(0, 100), breaks = c(0,20,40,60,80,100)) +             ## to tune the size of circles
  facet_grid(cols=vars(gene_function),scales='free',space='free')+
  labs(x='',y='')+
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "grey90",size=0.2), #panel.grid.minor = element_blank(),
        strip.background = element_blank(),strip.text = element_text(angle = 90),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle = 90,vjust=0.2,hjust=0.95) )

ggsave(paste0('TME.cell.type.level2.top4deg.BubblePlot.',Sys.Date(),'.pdf'),plotx,height=5.7,width=18)

###################### Figure S2B-E neutrophils ######################

seu.obj<-readRDS(paste0('MM.TME.235samples.','Neutrophil','.filter.scs.rds'))

gene<-c('S100A8','S100A9','CTSD','LYZ','LGALS3','FCGR3B','CSF3R',
        'CEACAM8',#'AZU1','ELANE','MPO',
        'IL1B','CXCL8','TNFAIP3','NFKBIA','IFITM1',
        'ARG1','OLR1','S100A12',
        'CXCR2','ITGAM','SELL',
        'MMP8','MMP9','LCN2','OLFM4','CAMP')

p<-DotPlot(seu.obj, features = gene,group.by = 'cell.status',col.max = 1.5,col.min=-1.5) 
data4<-p$data[,c('id','features.plot','pct.exp','avg.exp.scaled')]

data4$id<-factor(data4$id,levels=rev(levels(data4$id)))
plotx<-ggplot(data4, aes(x = features.plot,y = id)) +        ## global aes
  geom_point(aes(fill = avg.exp.scaled,size =pct.exp),color='black',shape=21)  +    ## geom_point for circle illusion
  scale_fill_gradientn(colours=c('#f7fcfd','#f7fcfd','#e0ecf4','#bfd3e6','#9ebcda','#8c96c6','#8c6bb1','#88419d','#810f7c','#4d004b')) + 
  scale_size(range = c(0,5), limits = c(0, 100), breaks = c(0,20,40,60,80,100)) +             ## to tune the size of circles
  #facet_grid(cols=vars(gene_function),scales='free',space='free')+
  labs(x='',y='')+
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "grey90",size=0.2), #panel.grid.minor = element_blank(),
        strip.background = element_blank(),strip.text = element_text(angle = 90),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle = 90,vjust=0.2,hjust=0.95) )
ggsave('Neutrophil.subset.gene.bubbleplot.pdf',plotx,height=3,width=6.5)



MM_TME.seu.obj<-readRDS('MM.TME.235samples.TME.no_Eryth.filter.scs.rds')

p<-DotPlot(MM_TME.seu.obj, features = gene,group.by = 'cell.type.level2',col.max = 1.5,col.min=-1.5)

data4<-p$data[,c('id','features.plot','pct.exp','avg.exp.scaled')]

data4$id<-factor(data4$id,levels=rev(levels(data4$id)))
plotx<-ggplot(data4, aes(x = features.plot,y = id)) +        ## global aes
  geom_point(aes(fill = avg.exp.scaled,size =pct.exp),color='black',shape=21)  +    ## geom_point for circle illusion
  scale_fill_gradientn(colours=c('#f7fcfd','#f7fcfd','#e0ecf4','#bfd3e6','#9ebcda','#8c96c6','#8c6bb1','#88419d','#810f7c','#4d004b')) + 
  scale_size(range = c(0,5), limits = c(0, 100), breaks = c(0,20,40,60,80,100)) +             ## to tune the size of circles
  #facet_grid(cols=vars(gene_function),scales='free',space='free')+
  labs(x='',y='')+
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "grey90",size=0.2), #panel.grid.minor = element_blank(),
        strip.background = element_blank(),strip.text = element_text(angle = 90),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle = 90,vjust=0.2,hjust=0.95) )
ggsave('Neutrophil.gene.bubbleplot.pdf',plotx,height=4,width=7)


plotx<-VlnPlot(MM_TME.seu.obj,features = 'nFeature_RNA',group.by = 'cell.type.level2',pt.size = 0,adjust=3,ncol = 1)&
  scale_fill_manual(values=cell.type.level2.color)&
  stat_summary(fun.y=mean,geom='point',position = position_dodge(width = 0.9))&
  theme(axis.text.x = element_text(size=12, color="black",angle = 90,vjust = 0.5),
        #axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),#axis.ticks=element_blank(),
        #strip.text.x = element_text(size=12, color="black", face="bold"),
        #strip.background = element_blank(),strip.text = element_text(angle = 0),
        #strip.text.y = element_text(angle = 0),
        text = element_text(size=15),
        #plot.title = element_blank(),#element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA, size=1))
ggsave(paste0('TME.cell.type.level2.nFeature.vlnplot.',Sys.Date(),'.pdf'),plotx,height=5,width=8)

deg_blood2025<-read_excel('blood_bld-2025-028963-mmc4.xlsx',sheet = 'DGE myeloid cells') %>% as.data.frame()
deg_blood2025_top50 <- deg_blood2025 %>% filter(p_val_adj<0.05) %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
deg_blood2025_top50 <- deg_blood2025_top50[deg_blood2025_top50$`cluster name`%notin%c('S100A8+CD14+MAC','HLA-DR+CD16+mac'),]

deg_TME<-read_excel('TableS2.xlsx',sheet = 'Neutrophil') %>% as.data.frame()


celltype1<-unique(deg_TME$cluster)
celltype2<-unique(deg_blood2025_top50$`cluster name`)

jaccard_index<-matrix(0,length(celltype1),length(celltype2))

for(i in 1:nrow(jaccard_index)){
  for(j in 1:ncol(jaccard_index)){
    
    I<-intersect(deg_TME$gene[deg_TME$cluster==celltype1[i]],deg_blood2025_top50$gene[deg_blood2025_top50$`cluster name`==celltype2[j]])
    U<-union(deg_TME$gene[deg_TME$cluster==celltype1[i]],deg_blood2025_top50$gene[deg_blood2025_top50$`cluster name`==celltype2[j]])
    jaccard_index[i,j]<-length(I)/50#length(U)
    
  }
}

rownames(jaccard_index)<-celltype1
colnames(jaccard_index)<-celltype2

jaccard_index[jaccard_index>0.3]<-0.3
#jaccard_index[jaccard_index<0.1]<-0

plotx<-pheatmap::pheatmap(t(jaccard_index),main='Jaccard index',
                          color=colorRampPalette(c('white','red'))(100),
                          show_colnames =T,
                          #angle_col = 0,
                          show_rownames =T,
                          #display_numbers=T,fontsize_number=3,
                          cluster_rows = T,
                          cluster_cols = T
)
ggsave('neutrophil.deg.overlap.with.public.deg.JaccardIndex.pdf',plotx,width=6,height=5)

###################### Figure S1D TME samples cell composition stackbarplot ######################
meta <- readRDS("MM.TME.235samples.TME.meta.data.rds") %>%
  dplyr::select(NewID, cell.type.level2, Diagnosis) %>%
  add_count(NewID, name = "n_cells") %>%
  filter(n_cells >= 200) %>%
  dplyr::select(-n_cells)

sample_comp <- meta %>%
  count(NewID, cell.type.level2, name = "Freq") %>%
  complete(NewID, cell.type.level2, fill = list(Freq = 0)) %>%
  left_join(meta %>% distinct(NewID, Diagnosis), by = "NewID") %>%
  group_by(NewID) %>%
  mutate(
    total = sum(Freq),
    percentage = Freq * 100 / total
  ) %>%
  ungroup() %>%
  as.data.frame()

Progenitor<-subset(sample_comp,cell.type.level2=='Progenitors')
Progenitor<-Progenitor%>%group_by(Diagnosis)%>%arrange(desc(percentage),.by_group=T)

sample_comp$NewID<-factor(sample_comp$NewID,levels=Progenitor$NewID)


plotx1<-ggplot(sample_comp, aes(x=NewID,fill=cell.type.level2,y=percentage)) +
  geom_bar(stat="identity",position = 'stack') +
  scale_fill_manual(values = cell.type.level2.color) +
  #scale_x_continuous(breaks=c(0,25,50,75,100)) +
  labs(x='',y = '% (of TME cells)') +
  facet_grid(.~Diagnosis,scales = 'free',space = 'free') +
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

plotx2<-ggplot(subset(sample_comp,cell.type.level2=='Progenitors'), aes(x=NewID,y=total)) +
  geom_bar(stat="identity",position = 'stack') +
  #scale_fill_manual(values = cell.type.level2.color) +
  #scale_x_continuous(breaks=c(0,25,50,75,100)) +
  labs(x='',y = '% (of TME cells)') +
  facet_grid(.~Diagnosis,scales = 'free',space = 'free') +
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

ggsave('TME.cell.type.level2.composition.stackbarplot.pdf',useDingbats=F,plotx,width = 20,height=10)
