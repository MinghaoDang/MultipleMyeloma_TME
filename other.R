library('CellChat')
library('patchwork')
library('ggplot2')
library('readxl')
library('ggpubr')
library('grid')
library('gridExtra')
library('ComplexHeatmap')
library('circlize')

options(stringsAsFactors = FALSE)

###################### CellChat ######################

subset_cellchat_by_celltype <- function(cellchat, celltypes_to_keep) {
  # Get the cell names for the selected cell types
  keep.cells <- colnames(cellchat@data)[cellchat@idents %in% celltypes_to_keep]
  
  meta <- cellchat@meta
  keep_cells <- rownames(meta)[meta$cell.type %in% celltypes_to_keep]
  data_input <- cellchat@data[, keep_cells]
  meta_subset <- meta[keep_cells, , drop = FALSE]
  
  cellchat.sub <- createCellChat(object = data_input, meta = meta_subset, group.by = "cell.type")
  cellchat.sub@DB <- CellChatDB.use
  
  cellchat.sub@idents<-factor(cellchat.sub@idents,levels=celltypes_to_keep)
  
  cellchat.sub <- subsetData(cellchat.sub)
  cellchat.sub <- identifyOverExpressedGenes(cellchat.sub)
  cellchat.sub <- identifyOverExpressedInteractions(cellchat.sub)
  cellchat.sub <- computeCommunProb(cellchat.sub, type = "truncatedMean") # 
  cellchat.sub <- filterCommunication(cellchat.sub)
  cellchat.sub <- computeCommunProbPathway(cellchat.sub)
  cellchat.sub <- aggregateNet(cellchat.sub)
  cellchat.sub <- netAnalysis_computeCentrality(cellchat.sub, slot.name = "netP")
  
  return(cellchat.sub)
}

CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB)

options(future.globals.maxSize = 20 * 1024^3)
future::plan("multisession", workers = 4) # do parallel


cellchat.list.by.Ecotype<-NULL
for(i in paste0('E',1:5)){
  
  print(i)

  seu.obj<-readRDS(paste0(i,'.seu.obj.rds'))
  seu.obj$samples<-seu.obj$NewID
  cellchat <- createCellChat(object = seu.obj, group.by = "cell.type.level2", assay = "RNA")
  
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 100)
  
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  
  cellchat.list.by.Ecotype[[length(cellchat.list.by.Ecotype)+1]]<-cellchat
  
}

names(cellchat.list.by.Ecotype)<-paste0('E',1:5)
saveRDS(cellchat.list.by.Ecotype,'cellchat.list.by.Ecotype.rds')

selected_types <- c('Mature_B','CD8T','CD4T','NK','gdT','MAIT','cPC','Progenitors',
                    'CD14_Mono','CD14CD16_Mono','CD16_Mono','Macrophage','Neutrophil',
                    'cDC','pDC')


cellchat_subset.list <- lapply(cellchat.list.by.Ecotype, function(obj) {
  subset_cellchat_by_celltype(obj, selected_types)
})

saveRDS(cellchat_subset.list,'cellchat_subset.list.by.Ecotype.rds')

###################### Spectra ######################
convert_mouse_to_human <- function(gene_list) { 
  output = c()
  mouse_human_genes = read.csv("https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
  
  for(gene in gene_list) {
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name == "mouse, laboratory"))[['DB.Class.Key']]
    if( !identical(class_key, integer(0)) ) {
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
      for(human_gene in human_genes) {
        output = rbind(c(gene, human_gene), output)
      }
    }
  }
  return (output)
}

convert_human_to_mouse <- function(gene_list) {
  output = c()
  mouse_human_genes = read.csv("https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
  
  for(gene in gene_list) {
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name == "human"))[['DB.Class.Key']]
    if( !identical(class_key, integer(0)) ) {
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="mouse, laboratory"))[,"Symbol"]
      for(human_gene in human_genes) {
        output = rbind(c(gene, human_gene), output)
      }
    }
  }
  return (output)
}

#### CD8T ####
load('MM.TME.235samples.CD8T.cluster.markers.RData')
data <- cluster.markers %>% filter(p_val_adj<0.01) %>% 
  group_by(cluster) %>% 
  top_n(n = 50, wt = dplyr::desc(p_val_adj)) %>% 
  top_n(n = 50, wt = avg_log2FC) %>% 
  as.data.frame()

ref_cd8<-readRDS('./Immune_Dictionary_ref/ref_data_T_cell_CD8.RDS')
cytokine<-unique(ref_cd8$sample)[unique(ref_cd8$sample)!='PBS']

enrichmentscore<-c()
for(i in unique(data$cluster)){
  
  d1<-subset(data,cluster==i)
  aa<-convert_human_to_mouse(d1$gene)
  bb<-FetchData(ref_cd8,vars = aa[,2],slot = 'data')
  cc<-rowSums(bb)
  dd<-data.frame('group'=ref_cd8$sample,'score'=cc)
  
  out.all<-c()
  for(j in cytokine){
    pval<-wilcox.test(dd$score[dd$group==j],dd$score[dd$group=='PBS'])
    es<-sum(dd$score[dd$group==j])-sum(dd$score[dd$group=='PBS'])
    out<-data.frame('cytokine'=j,'enrichmentscore'=es,'pval'=pval$p.value)
    out.all<-rbind(out.all,out)
  }
  
  out.all$padj<-p.adjust(out.all$pval,method = 'fdr')
  out.all$logfdr<- -log10(out.all$padj)
  out.all<-out.all[order(out.all$enrichmentscore,decreasing = T),]
  out.all$cytokine<-factor(out.all$cytokine,levels=rev(out.all$cytokine))
  out.all$cell.status=i
  
  enrichmentscore<-rbind(enrichmentscore,out.all)
  
}
write.table(enrichmentscore,file = 'CD8T.immune.dictionary.txt',sep = '\t',quote = F,col.names = T,row.names = F)

#### NK ####
load('MM.TME.235samples.NK.cluster.markers.RData')
data <- cluster.markers %>% filter(p_val_adj<0.01) %>% 
  group_by(cluster) %>% 
  top_n(n = 50, wt = dplyr::desc(p_val_adj)) %>% 
  top_n(n = 50, wt = avg_log2FC) %>% 
  as.data.frame()

ref_cd8<-readRDS('./Immune_Dictionary_ref/ref_data_NK_cell.RDS')
cytokine<-unique(ref_cd8$sample)[unique(ref_cd8$sample)!='PBS']

enrichmentscore<-c()
for(i in unique(data$cluster)){
  
  d1<-subset(data,cluster==i)
  aa<-convert_human_to_mouse(d1$gene)
  bb<-FetchData(ref_cd8,vars = aa[,2],slot = 'data')
  cc<-rowSums(bb)
  dd<-data.frame('group'=ref_cd8$sample,'score'=cc)
  
  out.all<-c()
  for(j in cytokine){
    pval<-wilcox.test(dd$score[dd$group==j],dd$score[dd$group=='PBS'])
    es<-sum(dd$score[dd$group==j])-sum(dd$score[dd$group=='PBS'])
    out<-data.frame('cytokine'=j,'enrichmentscore'=es,'pval'=pval$p.value)
    out.all<-rbind(out.all,out)
  }
  
  out.all$padj<-p.adjust(out.all$pval,method = 'fdr')
  out.all$logfdr<- -log10(out.all$padj)
  out.all<-out.all[order(out.all$enrichmentscore,decreasing = T),]
  out.all$cytokine<-factor(out.all$cytokine,levels=rev(out.all$cytokine))
  out.all$cell.status=i
  
  enrichmentscore<-rbind(enrichmentscore,out.all)
  
}
write.table(enrichmentscore,file = 'NK.immune.dictionary.txt',sep = '\t',quote = F,col.names = T,row.names = F)

#### Macrophage ####
load('MM.TME.235samples.Macrophage.cluster.markers.RData')
data <- cluster.markers %>% filter(p_val_adj<0.01) %>% 
  group_by(cluster) %>% 
  top_n(n = 50, wt = dplyr::desc(p_val_adj)) %>% 
  top_n(n = 50, wt = avg_log2FC) %>% 
  as.data.frame()

ref_cd8<-readRDS('./Immune_Dictionary_ref/ref_data_Macrophage.RDS')
cytokine<-unique(ref_cd8$sample)[unique(ref_cd8$sample)!='PBS']

enrichmentscore<-c()

for(i in unique(data$cluster)){
  
  d1<-subset(data,cluster==i)
  aa<-convert_human_to_mouse(d1$gene)
  bb<-FetchData(ref_cd8,vars = aa[,2],slot = 'data')
  cc<-rowSums(bb)
  dd<-data.frame('group'=ref_cd8$sample,'score'=cc)
  
  out.all<-c()
  for(j in cytokine){
    pval<-wilcox.test(dd$score[dd$group==j],dd$score[dd$group=='PBS'])
    es<-sum(dd$score[dd$group==j])-sum(dd$score[dd$group=='PBS'])
    out<-data.frame('cytokine'=j,'enrichmentscore'=es,'pval'=pval$p.value)
    out.all<-rbind(out.all,out)
  }
  
  out.all$padj<-p.adjust(out.all$pval,method = 'fdr')
  out.all$logfdr<- -log10(out.all$padj)
  out.all<-out.all[order(out.all$enrichmentscore,decreasing = T),]
  out.all$cytokine<-factor(out.all$cytokine,levels=rev(out.all$cytokine))
  out.all$cell.status=i
  
  enrichmentscore<-rbind(enrichmentscore,out.all)
  
}
write.table(enrichmentscore,file = 'Macrophage.immune.dictionary.txt',sep = '\t',quote = F,col.names = T,row.names = F)

#### CD4T ####
load('MM.TME.235samples.CD4T.cluster.markers.RData')
data <- cluster.markers %>% filter(p_val_adj<0.01) %>% 
  group_by(cluster) %>% 
  top_n(n = 50, wt = dplyr::desc(p_val_adj)) %>% 
  top_n(n = 50, wt = avg_log2FC) %>% 
  as.data.frame()

ref_cd8<-readRDS('./Immune_Dictionary_ref/ref_data_T_cell_CD4.RDS')
cytokine<-unique(ref_cd8$sample)[unique(ref_cd8$sample)!='PBS']


enrichmentscore<-c()
for(i in unique(data$cluster)){
  
  d1<-subset(data,cluster==i)
  aa<-convert_human_to_mouse(d1$gene)
  bb<-FetchData(ref_cd8,vars = aa[,2],slot = 'data')
  cc<-rowSums(bb)
  dd<-data.frame('group'=ref_cd8$sample,'score'=cc)
  
  out.all<-c()
  for(j in cytokine){
    pval<-wilcox.test(dd$score[dd$group==j],dd$score[dd$group=='PBS'])
    es<-sum(dd$score[dd$group==j])-sum(dd$score[dd$group=='PBS'])
    out<-data.frame('cytokine'=j,'enrichmentscore'=es,'pval'=pval$p.value)
    out.all<-rbind(out.all,out)
  }
  
  out.all$padj<-p.adjust(out.all$pval,method = 'fdr')
  out.all$logfdr<- -log10(out.all$padj)
  out.all<-out.all[order(out.all$enrichmentscore,decreasing = T),]
  out.all$cytokine<-factor(out.all$cytokine,levels=rev(out.all$cytokine))
  out.all$cell.status=i
  
  enrichmentscore<-rbind(enrichmentscore,out.all)
  
}
write.table(enrichmentscore,file = 'CD4T.immune.dictionary.txt',sep = '\t',quote = F,col.names = T,row.names = F)

#### Bcell ####
load('MM.TME.235samples.Bcell.cluster.markers.RData')
data <- cluster.markers %>% filter(p_val_adj<0.01) %>% 
  group_by(cluster) %>% 
  top_n(n = 50, wt = dplyr::desc(p_val_adj)) %>% 
  top_n(n = 50, wt = avg_log2FC) %>% 
  as.data.frame()

ref_b<-readRDS('./Immune_Dictionary_ref/ref_data_B_cell.RDS')
cytokine<-unique(ref_b$sample)[unique(ref_b$sample)!='PBS']

enrichmentscore<-c()

for(i in unique(data$cluster)){
  
  d1<-subset(data,cluster==i)
  aa<-convert_human_to_mouse(d1$gene)
  bb<-FetchData(ref_b,vars = aa[,2],slot = 'data')
  cc<-rowSums(bb)
  dd<-data.frame('group'=ref_b$sample,'score'=cc)
  
  out.all<-c()
  for(j in cytokine){
    pval<-wilcox.test(dd$score[dd$group==j],dd$score[dd$group=='PBS'])
    es<-sum(dd$score[dd$group==j])-sum(dd$score[dd$group=='PBS'])
    out<-data.frame('cytokine'=j,'enrichmentscore'=es,'pval'=pval$p.value)
    out.all<-rbind(out.all,out)
  }
  
  out.all$padj<-p.adjust(out.all$pval,method = 'fdr')
  out.all$logfdr<- -log10(out.all$padj)
  out.all<-out.all[order(out.all$enrichmentscore,decreasing = T),]
  out.all$cytokine<-factor(out.all$cytokine,levels=rev(out.all$cytokine))
  out.all$cell.status=i
  
  enrichmentscore<-rbind(enrichmentscore,out.all)
  
}
write.table(enrichmentscore,file = 'Bcell.immune.dictionary.txt',sep = '\t',quote = F,col.names = T,row.names = F)


#### Neutrophil ####
load('MM.TME.235samples.Neutrophil.cluster.markers.RData')
data <- cluster.markers %>% filter(p_val_adj<0.01) %>% 
  group_by(cluster) %>% 
  top_n(n = 50, wt = dplyr::desc(p_val_adj)) %>% 
  top_n(n = 50, wt = avg_log2FC) %>% 
  as.data.frame()

ref_b<-readRDS('./Immune_Dictionary_ref/ref_data_Neutrophil.RDS')
cytokine<-unique(ref_b$sample)[unique(ref_b$sample)!='PBS']

enrichmentscore<-c()

for(i in unique(data$cluster)){
  
  d1<-subset(data,cluster==i)
  aa<-convert_human_to_mouse(d1$gene)
  bb<-FetchData(ref_b,vars = aa[,2],slot = 'data')
  cc<-rowSums(bb)
  dd<-data.frame('group'=ref_b$sample,'score'=cc)
  
  out.all<-c()
  for(j in cytokine){
    pval<-wilcox.test(dd$score[dd$group==j],dd$score[dd$group=='PBS'])
    es<-sum(dd$score[dd$group==j])-sum(dd$score[dd$group=='PBS'])
    out<-data.frame('cytokine'=j,'enrichmentscore'=es,'pval'=pval$p.value)
    out.all<-rbind(out.all,out)
  }
  
  out.all$padj<-p.adjust(out.all$pval,method = 'fdr')
  out.all$logfdr<- -log10(out.all$padj)
  out.all<-out.all[order(out.all$enrichmentscore,decreasing = T),]
  out.all$cytokine<-factor(out.all$cytokine,levels=rev(out.all$cytokine))
  out.all$cell.status=i
  
  enrichmentscore<-rbind(enrichmentscore,out.all)
  
}
write.table(enrichmentscore,file = 'Neutrophil.immune.dictionary.txt',sep = '\t',quote = F,col.names = T,row.names = F)


#### CD14_Mono ####
load('MM.TME.235samples.CD14_Mono.cluster.markers.RData')
data <- cluster.markers %>% filter(p_val_adj<0.01) %>% 
  group_by(cluster) %>% 
  top_n(n = 50, wt = dplyr::desc(p_val_adj)) %>% 
  top_n(n = 50, wt = avg_log2FC) %>% 
  as.data.frame()

ref_b<-readRDS('./Immune_Dictionary_ref/ref_data_Monocyte.RDS')
cytokine<-unique(ref_b$sample)[unique(ref_b$sample)!='PBS']

enrichmentscore<-c()

for(i in unique(data$cluster)){
  
  d1<-subset(data,cluster==i)
  aa<-convert_human_to_mouse(d1$gene)
  bb<-FetchData(ref_b,vars = aa[,2],slot = 'data')
  cc<-rowSums(bb)
  dd<-data.frame('group'=ref_b$sample,'score'=cc)
  
  out.all<-c()
  for(j in cytokine){
    pval<-wilcox.test(dd$score[dd$group==j],dd$score[dd$group=='PBS'])
    es<-sum(dd$score[dd$group==j])-sum(dd$score[dd$group=='PBS'])
    out<-data.frame('cytokine'=j,'enrichmentscore'=es,'pval'=pval$p.value)
    out.all<-rbind(out.all,out)
  }
  
  out.all$padj<-p.adjust(out.all$pval,method = 'fdr')
  out.all$logfdr<- -log10(out.all$padj)
  out.all<-out.all[order(out.all$enrichmentscore,decreasing = T),]
  out.all$cytokine<-factor(out.all$cytokine,levels=rev(out.all$cytokine))
  out.all$cell.status=i
  
  enrichmentscore<-rbind(enrichmentscore,out.all)
  
}
write.table(enrichmentscore,file = 'CD14_Mono.immune.dictionary.txt',sep = '\t',quote = F,col.names = T,row.names = F)


#### CD16_Mono ####
load('MM.TME.235samples.CD16_Mono.cluster.markers.RData')
data <- cluster.markers %>% filter(p_val_adj<0.01) %>% 
  group_by(cluster) %>% 
  top_n(n = 50, wt = dplyr::desc(p_val_adj)) %>% 
  top_n(n = 50, wt = avg_log2FC) %>% 
  as.data.frame()

ref_b<-readRDS('./Immune_Dictionary_ref/ref_data_Monocyte.RDS')
cytokine<-unique(ref_b$sample)[unique(ref_b$sample)!='PBS']

enrichmentscore<-c()

for(i in unique(data$cluster)){
  
  d1<-subset(data,cluster==i)
  aa<-convert_human_to_mouse(d1$gene)
  bb<-FetchData(ref_b,vars = aa[,2],slot = 'data')
  cc<-rowSums(bb)
  dd<-data.frame('group'=ref_b$sample,'score'=cc)
  
  out.all<-c()
  for(j in cytokine){
    pval<-wilcox.test(dd$score[dd$group==j],dd$score[dd$group=='PBS'])
    es<-sum(dd$score[dd$group==j])-sum(dd$score[dd$group=='PBS'])
    out<-data.frame('cytokine'=j,'enrichmentscore'=es,'pval'=pval$p.value)
    out.all<-rbind(out.all,out)
  }
  
  out.all$padj<-p.adjust(out.all$pval,method = 'fdr')
  out.all$logfdr<- -log10(out.all$padj)
  out.all<-out.all[order(out.all$enrichmentscore,decreasing = T),]
  out.all$cytokine<-factor(out.all$cytokine,levels=rev(out.all$cytokine))
  out.all$cell.status=i
  
  enrichmentscore<-rbind(enrichmentscore,out.all)
  
}
write.table(enrichmentscore,file = 'CD16_Mono.immune.dictionary.txt',sep = '\t',quote = F,col.names = T,row.names = F)


#### cDC2 ####
load('MM.TME.235samples.DC.cluster.markers.RData')
data <- cluster.markers %>% filter(p_val_adj<0.01) %>% 
  group_by(cluster) %>% 
  top_n(n = 50, wt = dplyr::desc(p_val_adj)) %>% 
  top_n(n = 50, wt = avg_log2FC) %>% 
  as.data.frame()

ref_b<-readRDS('./Immune_Dictionary_ref/ref_data_cDC2.RDS')
cytokine<-unique(ref_b$sample)[unique(ref_b$sample)!='PBS']

enrichmentscore<-c()

for(i in c('cDC2_MNDA','cDC2_VEGFA')){
  
  d1<-subset(data,cluster==i)
  aa<-convert_human_to_mouse(d1$gene)
  bb<-FetchData(ref_b,vars = aa[,2],slot = 'data')
  cc<-rowSums(bb)
  dd<-data.frame('group'=ref_b$sample,'score'=cc)
  
  out.all<-c()
  for(j in cytokine){
    pval<-wilcox.test(dd$score[dd$group==j],dd$score[dd$group=='PBS'])
    es<-sum(dd$score[dd$group==j])-sum(dd$score[dd$group=='PBS'])
    out<-data.frame('cytokine'=j,'enrichmentscore'=es,'pval'=pval$p.value)
    out.all<-rbind(out.all,out)
  }
  
  out.all$padj<-p.adjust(out.all$pval,method = 'fdr')
  out.all$logfdr<- -log10(out.all$padj)
  out.all<-out.all[order(out.all$enrichmentscore,decreasing = T),]
  out.all$cytokine<-factor(out.all$cytokine,levels=rev(out.all$cytokine))
  out.all$cell.status=i
  
  enrichmentscore<-rbind(enrichmentscore,out.all)
  
}
write.table(enrichmentscore,file = 'cDC2.immune.dictionary.txt',sep = '\t',quote = F,col.names = T,row.names = F)


#### cDC1 ####
load('MM.TME.235samples.DC.cluster.markers.RData')
data <- cluster.markers %>% filter(p_val_adj<0.01) %>% 
  group_by(cluster) %>% 
  top_n(n = 50, wt = dplyr::desc(p_val_adj)) %>% 
  top_n(n = 50, wt = avg_log2FC) %>% 
  as.data.frame()

ref_b<-readRDS('./Immune_Dictionary_ref/ref_data_cDC1.RDS')
cytokine<-unique(ref_b$sample)[unique(ref_b$sample)!='PBS']

enrichmentscore<-c()

for(i in c('cDC1')){
  
  d1<-subset(data,cluster==i)
  aa<-convert_human_to_mouse(d1$gene)
  bb<-FetchData(ref_b,vars = aa[,2],slot = 'data')
  cc<-rowSums(bb)
  dd<-data.frame('group'=ref_b$sample,'score'=cc)
  
  out.all<-c()
  for(j in cytokine){
    pval<-wilcox.test(dd$score[dd$group==j],dd$score[dd$group=='PBS'])
    es<-sum(dd$score[dd$group==j])-sum(dd$score[dd$group=='PBS'])
    out<-data.frame('cytokine'=j,'enrichmentscore'=es,'pval'=pval$p.value)
    out.all<-rbind(out.all,out)
  }
  
  out.all$padj<-p.adjust(out.all$pval,method = 'fdr')
  out.all$logfdr<- -log10(out.all$padj)
  out.all<-out.all[order(out.all$enrichmentscore,decreasing = T),]
  out.all$cytokine<-factor(out.all$cytokine,levels=rev(out.all$cytokine))
  out.all$cell.status=i
  
  enrichmentscore<-rbind(enrichmentscore,out.all)
  
}
write.table(enrichmentscore,file = 'cDC1.immune.dictionary.txt',sep = '\t',quote = F,col.names = T,row.names = F)


#### pDC ####
load('MM.TME.235samples.DC.cluster.markers.RData')
data <- cluster.markers %>% filter(p_val_adj<0.01) %>% 
  group_by(cluster) %>% 
  top_n(n = 50, wt = dplyr::desc(p_val_adj)) %>% 
  top_n(n = 50, wt = avg_log2FC) %>% 
  as.data.frame()

ref_b<-readRDS('./Immune_Dictionary_ref/ref_data_pDC.RDS')
cytokine<-unique(ref_b$sample)[unique(ref_b$sample)!='PBS']

enrichmentscore<-c()

for(i in c('pDC')){
  
  d1<-subset(data,cluster==i)
  aa<-convert_human_to_mouse(d1$gene)
  bb<-FetchData(ref_b,vars = aa[,2],slot = 'data')
  cc<-rowSums(bb)
  dd<-data.frame('group'=ref_b$sample,'score'=cc)
  
  out.all<-c()
  for(j in cytokine){
    pval<-wilcox.test(dd$score[dd$group==j],dd$score[dd$group=='PBS'])
    es<-sum(dd$score[dd$group==j])-sum(dd$score[dd$group=='PBS'])
    out<-data.frame('cytokine'=j,'enrichmentscore'=es,'pval'=pval$p.value)
    out.all<-rbind(out.all,out)
  }
  
  out.all$padj<-p.adjust(out.all$pval,method = 'fdr')
  out.all$logfdr<- -log10(out.all$padj)
  out.all<-out.all[order(out.all$enrichmentscore,decreasing = T),]
  out.all$cytokine<-factor(out.all$cytokine,levels=rev(out.all$cytokine))
  out.all$cell.status=i
  
  enrichmentscore<-rbind(enrichmentscore,out.all)
  
}
write.table(enrichmentscore,file = 'pDC.immune.dictionary.txt',sep = '\t',quote = F,col.names = T,row.names = F)


#### LAMP3_DC ####
load('MM.TME.235samples.DC.cluster.markers.RData')
data <- cluster.markers %>% filter(p_val_adj<0.01) %>% 
  group_by(cluster) %>% 
  top_n(n = 50, wt = dplyr::desc(p_val_adj)) %>% 
  top_n(n = 50, wt = avg_log2FC) %>% 
  as.data.frame()

ref_b<-readRDS('./Immune_Dictionary_ref/ref_data_MigDC.RDS')
cytokine<-unique(ref_b$sample)[unique(ref_b$sample)!='PBS']

enrichmentscore<-c()

for(i in c('LAMP3_DC')){
  
  d1<-subset(data,cluster==i)
  aa<-convert_human_to_mouse(d1$gene)
  bb<-FetchData(ref_b,vars = aa[,2],slot = 'data')
  cc<-rowSums(bb)
  dd<-data.frame('group'=ref_b$sample,'score'=cc)
  
  out.all<-c()
  for(j in cytokine){
    pval<-wilcox.test(dd$score[dd$group==j],dd$score[dd$group=='PBS'])
    es<-sum(dd$score[dd$group==j])-sum(dd$score[dd$group=='PBS'])
    out<-data.frame('cytokine'=j,'enrichmentscore'=es,'pval'=pval$p.value)
    out.all<-rbind(out.all,out)
  }
  
  out.all$padj<-p.adjust(out.all$pval,method = 'fdr')
  out.all$logfdr<- -log10(out.all$padj)
  out.all<-out.all[order(out.all$enrichmentscore,decreasing = T),]
  out.all$cytokine<-factor(out.all$cytokine,levels=rev(out.all$cytokine))
  out.all$cell.status=i
  
  enrichmentscore<-rbind(enrichmentscore,out.all)
  
}
write.table(enrichmentscore,file = 'LAMP3_DC.immune.dictionary.txt',sep = '\t',quote = F,col.names = T,row.names = F)


###################### MMRF immune atlas label transfer ######################

MMRF_obj <- readRDS('COMBINED_VALIDATION_MMRF_ImmuneAtlas_Full_Censored_Metadata.rds')

MMRF_obj.use <- subset(
  MMRF_obj,
  lineage_group %notin% c('LQ','P') &
    VJ_INTERVAL == 'Baseline' &
    cohort == 'discovery'
)

MMRF_obj.use <- subset(
  MMRF_obj.use,
  !cellID_short %in% c(
    'EB','EB_MKi67_1','EB_MKi67_2',
    'Ery.0','Ery.1','Ery.2','Ery.3','Ery.4','Ery.5','Ery.6','Ery.8','Ery.9'
  )
)

saveRDS(MMRF_obj.use,'MMRF.263baseline.scs.rds')

MM_TME <- readRDS('MM.TME.235samples.TME.no_Eryth.filter.scs.rds')

process_seurat <- function(obj) {
  obj %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA()
}

run_transfer <- function(query, ref, ref_label, outfile) {
  anchors <- FindTransferAnchors(
    reference = ref,
    query = query,
    dims = 1:30,
    reference.reduction = "pca"
  )
  pred <- TransferData(
    anchorset = anchors,
    refdata = ref[[ref_label]],
    dims = 1:30
  )
  saveRDS(pred, outfile)
}

lineage_list <- list(
  CD8T = list(
    subset = expression(lineage_group == 'CD8' & cellID_short != 'MAIT'),
    ref = readRDS('MM.TME.235samples.CD8T.filter.scs.rds')
  ),
  CD4T = list(
    subset = expression(lineage_group == 'CD4'),
    ref = readRDS('MM.TME.235samples.CD4T.filter.scs.rds')
  ),
  B = list(
    subset = expression(lineage_group == 'B'),
    ref = subset(MM_TME, cell.status %in% c(
      'Bcell_Fos/Jun','Bcell_IgMhi_TransB','Bcell_Isg',
      'Bcell_Mem','Bcell_Naive','Pre/Immature_B',
      'Pro_B','Prolif_Pre_B','Prolif_Pro_B'
    ))
  ),
  Nk = list(
    subset = expression(lineage_group == 'Nk'),
    ref = readRDS('MM.TME.235samples.NK.filter.scs.rds')
  ),
  M = list(
    subset = expression(cellID_short %in% c(
      'CD14+Mono_CTSS','CD14+Mono_hypo','CD14+Mono_IFN',
      'CD14+Mono_pro-inflam','CD14+Mono_S100A','CD16+Mono',
      'cDC1','cDC2','M2_Macro','Macro/Mono','pDCs_b'
    )),
    ref = subset(MM_TME, cell.type.level2 %in% c(
      'CD14_Mono','CD14CD16_Mono','CD16_Mono','Macrophage','cDC'
    ))
  ),
  Other = list(
    subset = expression(cellID_short %in% c(
      'Fibro','HSC','pDCs_a','GMP','Granulocyte','Neutrophil_RPS/RPL',
      'MAIT','MastC','MegaK','Neutrohpil_ARG1'
    )),
    ref = subset(MM_TME, cell.status %in% c(
      'BaEoMa','Endothelial','Fibroblast','gdT','GMP','HSC','MAIT',
      'Neutro_ABCA13','Neutro_CAMP','Neutro_CXCL8','Neutro_Isg',
      'Neutro_MMP9','Neutro_NEAT1','Neutro_TENM1','pDC','Platelet',
      'Pro_E','Pro_Mk','Prolif_CD4T','Prolif_CD8T','Prolif_NK','Prolif_pDC'
    ))
  )
)

for (nm in names(lineage_list)) {
  
  cat("Processing:", nm, "\n")
  
  query <- subset(MMRF_obj.use, eval(lineage_list[[nm]]$subset))
  query <- process_seurat(query)
  
  ref <- lineage_list[[nm]]$ref
  ref <- process_seurat(ref)
  
  run_transfer(
    query,
    ref,
    ref_label = "cell.status",
    outfile = paste0("MMRF.", nm, ".predictions.rds")
  )
}

pred_files <- list.files(pattern = "MMRF\\..*\\.predictions.rds")

predictions_ALL <- lapply(pred_files, readRDS) %>%
  bind_rows() %>%
  tibble::rownames_to_column("Barcode") %>%
  select(Barcode, predicted.id)

meta <- MMRF_obj.use@meta.data %>%
  left_join(predictions_ALL, by = c("cellname" = "Barcode"))

MMRF_obj.use$predicted.id <- meta$predicted.id

saveRDS(meta, "MMRF.263baseline.meta.data.rds")
saveRDS(MMRF_obj.use, "MMRF.263baseline.scs.rds")
