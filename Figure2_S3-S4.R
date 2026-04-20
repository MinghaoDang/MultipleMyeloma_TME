source('global_config.R')

###################### Figure 2A lymphoid cell status compare boxplot ######################

meta<-readRDS('MM.TME.235samples.TME.meta.data.rds')

for(celltype in c('CD8T','CD4T','NK','Mature_B')){
  
  meta.sub<-subset(meta,cell.type.level2==celltype)
  
  meta.sub <- meta.sub %>% add_count(NewID, name = "n_cells") %>%
    filter(n_cells > 50) %>%
    dplyr::select(-n_cells)
  
  sample_comp <- meta.sub %>%
    count(NewID, cell.status, name = "Freq") %>%
    complete(NewID, cell.status, fill = list(Freq = 0)) %>%
    left_join(meta.sub %>% distinct(NewID, Diagnosis), by = "NewID") %>%
    group_by(NewID) %>%
    mutate(
      total = sum(Freq),
      percentage = Freq * 100 / total
    ) %>%
    ungroup() %>%
    as.data.frame()
  
  library('gridExtra')
  library('grid')
  library('purrr')
  library('ggsignif')
  
  plotx<-signif_plot(sample_comp,
                     x='Diagnosis',
                     y='percentage',
                     split='cell.status',
                     signif.cutoff=0.05,
                     ncol=11,
                     fill_by='Diagnosis',
                     fill_color=c('nBM'='#4dac26','MGUS'='#1c75bc','SMM'='#fbb040','NDMM'='#d01c8b', 'RRMM'='#8856a7'),
                     ylab=paste0('% (of ',celltype,' cells)'))
  
  ggsave(paste0(celltype,'.cell.status.compare.pdf'),plotx,width=30,height=5)
  
}

###################### Figure 2B CD8 TCR stackplot ######################
CD8T.obj<-readRDS('MM.TME.235samples.CD8T.filter.scs.rds')

data<-CD8T.obj@meta.data[,c('NewID','Diagnosis','cell.status','TCR_clonotype_freq')]
data<-data[!is.na(data$TCR_clonotype_freq),]
data$group=ifelse(data$TCR_clonotype_freq==1,'=1',
                  ifelse(data$TCR_clonotype_freq<=5,'2-5',
                         ifelse(data$TCR_clonotype_freq<=10,'>5',
                                ifelse(data$TCR_clonotype_freq<=20,'>10',
                                       ifelse(data$TCR_clonotype_freq<=50,'>20',
                                              ifelse(data$TCR_clonotype_freq<=100,'>50','>100'))))))

data1=as.data.frame(table(data[,c('cell.status','group')]))
data2 = data1 %>% group_by(cell.status) %>% summarise(count=sum(Freq))
data3<-inner_join(data1,data2,by='cell.status')
data3$percentage<- data3$Freq*100/data3$count;

data3$group<-factor(data3$group,levels=c('>100','>50','>20','>10','>5','2-5','=1'))

color=c('#fef0d9','#fdd49e','#fdbb84','#fc8d59','#ef6548','#d7301f','#990000')
names(color)=c('=1','2-5','>5','>10','>20','>50','>100')

plotx<-ggplot(data3, aes(x=cell.status,fill=group,y=percentage)) +
  geom_bar(stat="identity",position = 'stack') +
  scale_fill_manual(values = color) +
  scale_y_continuous(breaks=c(0,25,50,75,100)) +
  labs(x='',y = '% (of CD8T cells)') +
  theme_vlnplot()

ggsave('CD8T.TCR.diversity.stackplot.pdf',useDingbats=F,plotx,width = 6,height=6)

###################### Figure 2C CD8 monocle2 ######################
library('monocle')

CD8T.obj<-readRDS('MM.TME.235samples.CD8T.filter.scs.rds')

set.seed(12345)
obj = CD8T.obj[,sort(sample(1:ncol(CD8T.obj),10000))]
expr_matrix <- Seurat::GetAssayData(obj, slot = 'counts')

pd <- data.frame(data = obj@meta.data)
fd <- data.frame(gene_short_name = row.names(expr_matrix), row.names = row.names(expr_matrix))

pd <- new("AnnotatedDataFrame", data = pd)
fd <- new("AnnotatedDataFrame", data = fd)

cds <- newCellDataSet(expr_matrix, phenoData = pd, featureData = fd)

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- detectGenes(cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 100))

cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- orderCells(cds)

plotx<-plot_cell_trajectory(cds, color_by = "data.cell.status",cell_size = 2,color='white')+
  scale_color_manual(values=CD8T_color)+
  theme_umap()
ggsave(paste0('CD8T.downsampling10000.monocle2.color.by.cell.status.pdf'),plotx,height=8,width=7)  


plotx<-plot_cell_trajectory(cds, color_by = "Pseudotime",cell_size = 2,color='white')+
  scale_color_viridis(option = "C")+
  theme_umap()
ggsave(paste0('CD8T.downsampling10000.monocle2.color.by.Pseudotime.pdf'),plotx,height=8,width=7)  


data<-pData(cds)[,c('State','data.cell.status')]

data1 = as.data.frame(table(data))
data2 = data1 %>% group_by(State) %>% summarise(count=sum(Freq))
data3 <- inner_join(data1,data2,by='State')
data3$percentage<- data3$Freq*100/data3$count;

plotx<-ggplot(data3, aes(x = '', y = percentage, fill = data.cell.status)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  scale_fill_manual(values=CD8T_color)+
  facet_wrap(~State,ncol = 7)+
  theme_void()
ggsave('CD8T.downsampling10000.monocle2.cell_type.piechart.pdf',useDingbats=F,plotx,width = 21,height=3)

###################### Figure 2C CD8 CytoTRACE ######################

Sys.setenv(RETICULATE_PYTHON="/bin/python3")
library('CytoTRACE')

dir.create('./CytoTRACE')

CD8T.obj<-readRDS('MM.TME.235samples.CD8T.filter.scs.rds')

umap_cord<-FetchData(CD8T.obj,vars=c('UMAP_1','UMAP_2'))
pheno<-as.character(CD8T.obj$cell.status)
names(pheno)<-names(CD8T.obj$cell.status)
expr<-as.matrix(GetAssayData(CD8T.obj,slot='data'))
results <- CytoTRACE(expr,ncores = 4)
CD8T.obj$CytoTRACE<-results$CytoTRACE[colnames(CD8T.obj)]

plotCytoTRACE(results,phenotype=pheno,emb=umap_cord,outputDir = './CytoTRACE/')

aa<-read.table('./CytoTRACE/CytoTRACE_plot_table.txt',sep='\t',header = T,row.names = 1)
aa$Phenotype<-factor(aa$Phenotype,levels=c('CD8T_Tn','CD8T_Tcm','CD8T_Trm','CD8T_Tem','CD8T_Teff_MHCII','CD8T_Teff_FGFBP2','CD8T_Teff_NKlike','CD8T_Tisg','CD8T_Teff_TGFB1','CD8T_Teff_CCL4','CD8T_Teff_Fos/Jun'))

CD8T_color<-c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#FF7F00","#E31A1C",
              "#FDBF6F","#CAB2D6","#6A3D9A","#FFFF99")
names(CD8T_color)<-c('CD8T_Tn','CD8T_Tcm','CD8T_Tem','CD8T_Trm','CD8T_Teff_CCL4','CD8T_Teff_MHCII','CD8T_Teff_Fos/Jun','CD8T_Teff_TGFB1','CD8T_Teff_FGFBP2','CD8T_Teff_NKlike','CD8T_Tisg')

plotx<-ggplot(aa, aes(x=Phenotype, y=CytoTRACE) ) + 
  geom_violin(adjust=0.25,bw=0.3,scale='width',aes(fill=factor(Phenotype))) + 
  stat_summary(fun.y=median,geom='point',
               position = position_dodge(width = 0.9)) +
  scale_fill_manual(values=CD8T_color)+
  theme(axis.text.x = element_text(size = 12,angle=90,vjust=0.2,hjust=0.95),
        strip.text.x = element_text(size=12, color="black", face="bold"),
        text = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA, size=1), legend.position = "none") +
  labs(x="",y='Predited ordering by CytoTRACE')

ggsave(paste0('CD8T.CytoTRACE.violin.pdf'),plotx,height=8,width=8)


###################### Figure 2D TCR diversity across disease spectrum ######################

library('tidyr')
library('emmeans')
library('patchwork')

diagnosis_colors <- c(
  "SMM"  = "#fbb040",
  "NDMM" = "#d01c8b",
  "RRMM" = "#8856a7")

dat <- readRDS("combined.TCR.rds")

calc_shannon <- function(x) {
  p <- table(x) / length(x)
  -sum(p * log(p))
}


# 1. preprocess

dat2 <- dat %>%
  filter(!is.na(CTaa), CTaa != "", CTaa != "NA") %>%
  filter(!Diagnosis %in% c("nBM", "MGUS")) %>%
  mutate(Diagnosis = factor(Diagnosis, levels = c("SMM", "NDMM", "RRMM")))

sample_summary <- dat2 %>%
  group_by(sample, Diagnosis) %>%
  summarise(
    Tcell_count = n(),
    Shannon = calc_shannon(CTaa),
    n_clonotype = n_distinct(CTaa),
    .groups = "drop"
  ) %>%
  filter(Tcell_count >= 100)

dat2 <- dat2 %>%
  filter(sample %in% sample_summary$sample)

sample_summary <- dat2 %>%
  group_by(sample, Diagnosis) %>%
  summarise(
    Tcell_count = n(),
    Shannon = calc_shannon(CTaa),
    n_clonotype = n_distinct(CTaa),
    .groups = "drop"
  ) %>%
  mutate(
    log10_Tcell_count = log10(Tcell_count),
    norm_Shannon = Shannon / log(Tcell_count),
    clonality = 1 - norm_Shannon
  )

# pairwise comparisons
my_comparisons <- list(
  c("SMM", "NDMM"),
  c("SMM", "RRMM"),
  c("NDMM", "RRMM")
)


# 2. Panel A: T-cell count vs Shannon

cor_res <- cor.test(sample_summary$Tcell_count, sample_summary$Shannon, method = "pearson")

cor_label <- paste0(
  "Pearson r = ", round(unname(cor_res$estimate), 3),
  "\nP = ", signif(cor_res$p.value, 3)
)

pA <- ggplot(sample_summary, aes(x = Tcell_count, y = Shannon, color = Diagnosis)) +
  geom_point(size = 2.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, color = "black") +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = cor_label,
    hjust = 1.1, vjust = 1.5,
    size = 4
  ) +
  scale_color_manual(values = diagnosis_colors) +
  theme_classic() +
  labs(x = "T-cell count", y = "Shannon entropy")


# 3. Panel B: Original Shannon

ymax_B <- max(sample_summary$Shannon, na.rm = TRUE)

pB <- ggplot(sample_summary, aes(x = Diagnosis, y = Shannon, color = Diagnosis)) +
  geom_boxplot(outlier.shape = NA, width = 0.65) +
  geom_jitter(width = 0.15, size = 2.2) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.format",
    hide.ns = FALSE
  ) +
  scale_color_manual(values = diagnosis_colors) +
  theme_classic() +
  labs(x = NULL, y = "Shannon entropy") +
  coord_cartesian(ylim = c(min(sample_summary$Shannon), ymax_B * 1.25))


# 4. Downsampled Shannon

depth <- 100
nrep <- 100
set.seed(123)

downsample_once <- function(df, depth) {
  df %>%
    group_by(sample, Diagnosis) %>%
    slice_sample(n = depth) %>%
    summarise(
      Shannon_down = calc_shannon(CTaa),
      .groups = "drop"
    )
}

ds_res <- map_dfr(seq_len(nrep), function(i) {
  downsample_once(dat2, depth) %>%
    mutate(rep = i)
})

ds_summary <- ds_res %>%
  group_by(sample, Diagnosis) %>%
  summarise(
    Shannon_down_mean = mean(Shannon_down),
    Shannon_down_sd = sd(Shannon_down),
    .groups = "drop"
  )

ymax_C <- max(ds_summary$Shannon_down_mean, na.rm = TRUE)

pC <- ggplot(ds_summary, aes(x = Diagnosis, y = Shannon_down_mean, color = Diagnosis)) +
  geom_boxplot(outlier.shape = NA, width = 0.65) +
  geom_jitter(width = 0.15, size = 2.2) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.format",
    hide.ns = FALSE
  ) +
  scale_color_manual(values = diagnosis_colors) +
  theme_classic() +
  labs(x = NULL, y = "Downsampled Shannon entropy") +
  coord_cartesian(ylim = c(min(ds_summary$Shannon_down_mean), ymax_C * 1.25))


# 5. Panel D: Normalized Shannon

ymax_D <- max(sample_summary$norm_Shannon, na.rm = TRUE)

pD <- ggplot(sample_summary, aes(x = Diagnosis, y = norm_Shannon, color = Diagnosis)) +
  geom_boxplot(outlier.shape = NA, width = 0.65) +
  geom_jitter(width = 0.15, size = 2.2) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.format",
    hide.ns = FALSE
  ) +
  scale_color_manual(values = diagnosis_colors) +
  theme_classic() +
  labs(x = NULL, y = "Normalized Shannon entropy") +
  coord_cartesian(ylim = c(min(sample_summary$norm_Shannon), ymax_D * 1.25))


# 6. Adjust for T-cell count

fit <- lm(Shannon ~ log10_Tcell_count + Diagnosis, data = sample_summary)

library(emmeans)

emm <- emmeans(fit, ~ Diagnosis)
emm_df <- as.data.frame(emm)

pairwise_res <- pairs(emm)
pairwise_df <- as.data.frame(pairwise_res)
pairwise_df <- pairwise_df %>%
  separate(contrast, into = c("group1", "group2"), sep = " - ")

pE <- ggplot(emm_df,
             aes(x = Diagnosis, y = emmean, color = Diagnosis)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.15) +
  scale_color_manual(values = diagnosis_colors) +
  theme_classic() +
  labs(y = "Adjusted Shannon (EMM)", x = NULL) +
  stat_pvalue_manual(
    pairwise_df,
    label = "p.value",
    xmin = "group1",
    xmax = "group2",
    y.position = max(emm_df$emmean) * c(1.05, 1.12, 1.19)
  )


# 7. Panel F: Rarefaction curves
#   depth from 50 to 500 by 50
#   only keep samples with >=500 cells

depth_seq <- seq(50, 500, by = 50)

dat_rare <- dat2 %>%
  group_by(sample, Diagnosis) %>%
  filter(n() >= 500) %>%
  ungroup()

rarefy_one_sample <- function(df, depth_seq, nrep = 50) {
  map_dfr(depth_seq, function(d) {
    map_dfr(seq_len(nrep), function(i) {
      tibble(
        depth = d,
        richness = df %>%
          slice_sample(n = d) %>%
          summarise(richness = n_distinct(CTaa)) %>%
          pull(richness)
      )
    }) %>%
      summarise(
        depth = d,
        richness_mean = mean(richness),
        richness_sd = sd(richness)
      )
  })
}

rare_df <- dat_rare %>%
  group_split(sample, Diagnosis) %>%
  map_dfr(function(x) {
    sample_id <- unique(x$sample)
    dx <- unique(x$Diagnosis)
    rarefy_one_sample(x, depth_seq, nrep = 50) %>%
      mutate(sample = sample_id, Diagnosis = dx)
  })

# Group-level mean ± SE across samples
rare_diag_df <- rare_df %>%
  group_by(Diagnosis, depth) %>%
  summarise(
    mean_richness = mean(richness_mean),
    sd_richness = sd(richness_mean),
    n_sample = n(),
    se_richness = sd_richness / sqrt(n_sample),
    .groups = "drop"
  )

pF <- ggplot(rare_diag_df, aes(x = depth, y = mean_richness, color = Diagnosis, fill = Diagnosis)) +
  geom_ribbon(
    aes(ymin = mean_richness - se_richness, ymax = mean_richness + se_richness),
    alpha = 0.2,
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = diagnosis_colors) +
  scale_fill_manual(values = diagnosis_colors) +
  scale_x_continuous(breaks = depth_seq) +
  theme_classic() +
  labs(
    x = "Number of T cells sampled",
    y = "Observed clonotype richness"
  )


# 8. combine panels

final_plot <- (pA | pB | pC) / (pD | pE | pF)

ggsave("TCR.control.for.Tcell_number.pdf", final_plot, width = 12, height = 6)

###################### Figure 2E Myeloid signature score radar plot ######################

library('fmsb')

data<-readRDS('Myeloid.signature.score.rds')

max_min<-data.frame(
  M1_score1 = c(max(data$median[data$variable=='M1_score1']), min(data$median[data$variable=='M1_score1'])), 
  M2_score2 = c(max(data$median[data$variable=='M2_score2']), min(data$median[data$variable=='M2_score2'])), 
  Angio_score3 = c(max(data$median[data$variable=='Angio_score3']), min(data$median[data$variable=='Angio_score3'])),
  phago_score4 = c(max(data$median[data$variable=='phago_score4']), min(data$median[data$variable=='phago_score4'])), 
  MDSC_score5 = c(max(data$median[data$variable=='MDSC_score5']), min(data$median[data$variable=='MDSC_score5']))
)
rownames(max_min) <- c("Max", "Min")

create_beautiful_radarchart <- function(data, color = "#00AFBB", 
                                        vlabels = colnames(data), vlcex = 0.7,
                                        caxislabels = NULL, title = NULL, ...){
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.25), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 'dashed', cglwd = 0.8,
    # Customize the axis
    #axislabcol = "grey", 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}

pdf('Myeloid.signature.score.radar.pdf')

for(i in unique(data$cell.type)){
  
  d1<-subset(data,cell.type==i)
  d1<-reshape2::dcast(d1,Diagnosis~variable)
  rownames(d1)<-d1$Diagnosis
  d1$Diagnosis<-NULL
  
  df <- rbind(max_min, d1)
  
  create_beautiful_radarchart(df,color=c('#4dac26','#1c75bc','#fbb040','#d01c8b', '#8856a7'), title=i)
  
}

dev.off()

###################### Figure 2F Myeloid cell status compare boxplot ######################
meta<-readRDS('MM.TME.235samples.TME.meta.data.rds')

meta<-subset(meta,cell.type%in%c('CD14_Mono','CD14CD16_Mono','CD16_Mono','cDC',
                                          'Macrophage','Neutrophil','pDC'))

meta <- meta %>% add_count(NewID, name = "n_cells") %>%
  filter(n_cells >= 100) %>%
  dplyr::select(-n_cells)

sample_comp <- meta %>%
  count(NewID, cell.status, name = "Freq") %>%
  complete(NewID, cell.status, fill = list(Freq = 0)) %>%
  left_join(meta %>% distinct(NewID, Diagnosis), by = "NewID") %>%
  group_by(NewID) %>%
  mutate(
    total = sum(Freq),
    percentage = Freq * 100 / total
  ) %>%
  ungroup() %>%
  as.data.frame()

library('gridExtra')
library('grid')
library('purrr')
library('ggsignif')

plotx<-signif_plot(sample_comp,
                   x='Diagnosis',
                   y='percentage',
                   split='cell.status',
                   signif.cutoff=0.05,
                   ncol=11,
                   fill_by='Diagnosis',
                   fill_color=c('nBM'='#4dac26','MGUS'='#1c75bc','SMM'='#fbb040','NDMM'='#d01c8b', 'RRMM'='#8856a7'),
                   ylab='% (of Myeloid cells)')
ggsave('Myeloid.cell.status.compare2.pdf',plotx,width=30,height=13)

###################### Figure 2G vlnplot ######################

CD14_Mono.seu.obj<-readRDS('MM.TME.235samples.CD14_Mono.filter.scs.rds')

aa<-FetchData(CD14_Mono.seu.obj,vars = c('Diagnosis','HLA-DRA','CD74'))

bb<-reshape2::melt(aa)
bb$value<-bb$value + rnorm(n = length(x =bb$value)) / 100000

P<-ggplot(bb,aes(x=Diagnosis,y=value,fill=Diagnosis)) +
  geom_violin() + # adjust=1.5,bw=0.3,scale='area',color='white' scale='width',
  stat_summary(fun.y=median,geom='point',
               position = position_dodge(width = 0.9)) +
  ggsignif::geom_signif(test="wilcox.test", comparisons = list(c('nBM','MGUS'),c('MGUS','SMM'),c('SMM','NDMM'),c('NDMM','RRMM')),
                        map_signif_level = F,   # c("***"=0.001, "**"=0.01, "*"=0.05)
                        vjust=-0.15,
                        textsize=6,
                        size=0.75,
                        step_increase = 0.15)+
  scale_fill_manual(values=c('nBM'='#4dac26','MGUS'='#1c75bc','SMM'='#fbb040','NDMM'='#d01c8b', 'RRMM'='#8856a7'))+
  facet_grid(cols=vars(variable),scales = 'free',space='free_x') +
  theme(strip.text.x = element_text(size=20, color="black", face="bold"),
        axis.text.x = element_text(angle = 90,size = 15,vjust=0.5,hjust=0.95), element_blank(), #
        text = element_text(size=20),
        #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA)) +
  labs(x='', y="Score")

ggsave('CD14_Mono.MHCII.gene.vlnplot.pdf',useDingbats=F,P,width=20,height=10)



###################### Figure 2H, S4B Myeloid subsets genes bubbleplot ######################

TME<-readRDS('MM.TME.235samples.TME.no_Eryth.filter.scs.rds')

Myeloid<-subset(TME,cell.type%in%c('CD14_Mono','CD14CD16_Mono','CD16_Mono','Macrophage','cDC','pDC',
                                          'Neutrophil'))

Idents(Myeloid)<-factor(Myeloid$cell.status)

M1<-c('IL23','TNF','CXCL9','CXCL10','CXCL11','CD86','IL1A','IL1B',
      'IL6','CCL5','IRF5','IRF1', 'CD40','IDO1','KYNU','CCR7' )
M2<-c('IL4R','CCL4','CCL13','CCL20','CCL17', 
      'CCL18','CCL22','CCL24','LYVE1','VEGFA','VEGFB','VEGFC', 
      'VEGFD','EGF','CTSA','CTSB','CTSC','CTSD','TGFB1','TGFB2', 
      'TGFB3','MMP14','MMP19','MMP9','CLEC7A','WNT7B','FASL', 
      'TNFSF12','TNFSF8','CD276','VTCN1','MSR1','FN1','IRF4' )
Angiogenesis<-c('CCND2','CCNE1','CD44','CXCR4','E2F3','EDN1','EZH2','FGF18',
                'FGFR1','FYN','HEY1','ITGAV','JAG1','JAG2','MMP9','NOTCH1', 
                'PDGFA','PTK2','SPP1','STC1','TNFAIP6','TYMP','VAV2','VCAN','VEGFA')
Phagocytosis<-c('MRC1','CD163','MERTK','C1QB' )
MDSC<-unique(c('STAT1','STAT3','STAT6','NFKB1','S100A8','S100A9','ANXA1',
               'LYZ','CXCL1','CXCL2','CXCL8','LILRA3','TREM1',
               'PTGS2','ARG1','ARG2','WFDC17','OLR1','CD84','NOS2','IL10'))
checkpoint.gene<-c('SIRPA','VSIR','LILRB1','LILRB2','SIGLEC10','TREM2','CD200R1','HAVCR2','PDCD1LG2','CD274','CSF1R','AXL','PILRA','NECTIN2','PVR','BTLA','ADORA2A')

gene<-unique(c(M1,M2,Angiogenesis,Phagocytosis,MDSC,checkpoint.gene))

p<-DotPlot(Myeloid, features = gene,col.max = 1.5,col.min=-1.5) 
data<-p$data[,c('id','features.plot','pct.exp','avg.exp.scaled')]
gene.sel<-data %>% group_by(features.plot) %>% 
  top_n(n=1, wt = pct.exp) %>% 
  filter(pct.exp>=20) %>% 
  pull(features.plot) %>%
  as.character()
data<-data[data$features.plot%in%gene.sel,]


data$gene_function='M1'
data$gene_function[which(data$features.plot%in%M2)]='M2'
data$gene_function[which(data$features.plot%in%Angiogenesis)]='Angio'
data$gene_function[which(data$features.plot%in%Phagocytosis)]='Phago'
data$gene_function[which(data$features.plot%in%MDSC)]='MDSC'
data$gene_function[which(data$features.plot%in%checkpoint.gene)]='checkpoint.gene'

data$gene_function<-factor(data$gene_function,levels=c('M1','M2','Angio','Phago','MDSC','checkpoint.gene'))

data$cell.type<-as.character(data$id)
data$cell.type[grepl('DC',data$id)]<-'DC'
data$cell.type[grepl('CD14_Mono',data$id)]<-'CD14_Mono'
data$cell.type[grepl('CD14CD16_Mono',data$id)]<-'CD14CD16_Mono'
data$cell.type[grepl('^CD16_Mono',data$id)]<-'CD16_Mono'
data$cell.type[grepl('Mac',data$id)]<-'Macro'
data$cell.type[grepl('Neutro',data$id)]<-'Neutro'


plotx<-ggplot(data, aes(x = features.plot,y = id)) +        ## global aes
  geom_point(aes(fill = avg.exp.scaled,size =pct.exp),color='black',shape=21)  +    ## geom_point for circle illusion
  #scale_fill_gradientn(colours=rev(color),limits=c(0,max(data$avg.exp)))+       ## color of the corresponding aes
  scale_fill_gradientn(colours=c('#f7fcfd','#f7fcfd','#e0ecf4','#bfd3e6','#9ebcda','#8c96c6','#8c6bb1','#88419d','#810f7c','#4d004b')) + 
  scale_size(range = c(0,5), limits = c(0, 100), breaks = c(0,20,40,60,80,100)) +             ## to tune the size of circles
  facet_grid(cols=vars(gene_function),rows=vars(cell.type),scales='free',space='free')+
  labs(x='',y='')+
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "grey90",size=0.2), #panel.grid.minor = element_blank(),
        strip.background = element_blank(),strip.text.x = element_text(angle = 0),
        strip.text.y = element_text(angle = 0),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle = 90,vjust=0.2,hjust=0.95) )
ggsave(paste0('Myeloid.functional.gene.BubblePlot.',Sys.Date(),'.pdf'),plotx,height=6.5,width=16)


###################### Figure S3A lymphoid cell status diagnosis fold change ######################
meta<-readRDS('MM.TME.235samples.TME.meta.data.rds')

diagnosis_levels <- c("nBM", "MGUS", "SMM", "NDMM", "RRMM")
cell.status.fc<-c()

for(celltype in c('CD8T','CD4T','NK','Mature_B')){ 
  
  meta.sub<-subset(meta,cell.type.level2==celltype)
  
  meta.sub=meta.sub[,c('NewID','cell.status','Diagnosis')]
  
  if(nrow(meta.sub)>10000){
    meta.sub<-subset(meta.sub,NewID%in%names(which(table(meta.sub$NewID)>50)))
  } else{meta.sub<-subset(meta.sub,NewID%in%names(which(table(meta.sub$NewID)>20)))}
  
  
  cell_status_fc_current <- meta.sub %>%
    count(NewID, cell.status, name = "Freq") %>%
    group_by(NewID) %>%
    mutate(
      total = sum(Freq),
      percentage = 100 * Freq / total
    ) %>%
    ungroup() %>%
    left_join(
      meta.sub %>% distinct(NewID, Diagnosis),
      by = "NewID"
    ) %>%
    group_by(Diagnosis, cell.status) %>%
    summarise(
      mean = mean(percentage),
      .groups = "drop"
    ) %>%
    mutate(Diagnosis = factor(Diagnosis, levels = diagnosis_levels)) %>%
    reshape2::dcast(cell.status ~ Diagnosis, value.var = "mean") %>%
    mutate(
      across(all_of(diagnosis_levels[-1]), ~ log2(.x / .data[[diagnosis_levels[1]]]))
    ) %>%
    mutate(cell.type = celltype)
  
  cell.status.fc <- bind_rows(cell.status.fc, cell_status_fc_current)
  
}

cell.status.fc$nBM<-0

aa<-reshape2::melt(cell.status.fc)

colnames(aa)[c(3,4)]<-c('Diagnosis','Log2FC')


last_points <- aa %>%
  group_by(cell.status) %>%
  filter(Diagnosis == "RRMM")

aa$cell.type<-factor(aa$cell.type,levels=c("CD8T","CD4T","NK","Mature_B"))
p<-ggplot(aa, aes(x=Diagnosis, y=Log2FC) ) + 
  #geom_bar(stat="identity") + 
  geom_line(aes(group = cell.status, color = cell.status)) +
  geom_point(aes(color=cell.status),size=1) +
  geom_text(data = last_points, aes(label = cell.status), hjust = -0.1, size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "red")+
  facet_wrap(.~cell.type,ncol=5,scales='free') +
  scale_color_manual(values=cell.status.color) +
  theme_clean() +
  theme(legend.position="none")+
  coord_cartesian(clip = "off")+
  labs(x='', y="Log2FC") 
ggsave(paste0('Lymphoid.cell.status.compare.dotline.',Sys.Date(),'.pdf'),p,width=15,height=5)


###################### Figure S3B CD8 signature score ######################
CD8T.obj<-readRDS('MM.TME.235samples.CD8T.filter.scs.rds')

gene.list<-read_excel('41591_2023_2371_MOESM3_ESM.xlsx',sheet='Table S4',skip = 1,col_names = T) # Chu et al., Nature Medicine 2023
gene.list<-as.list(gene.list)

gene.list<-lapply(gene.list, function(x) x[!is.na(x)])


CD8T.obj<-AddModuleScore(CD8T.obj,gene.list,name=names(gene.list))

gp.list<-NULL
for(i in colnames(CD8T.obj@meta.data)[24:42]){
  
  p<-FeaturePlot(object = CD8T.obj, reduction = "umap",features=i,pt.size = 0.75,label = F,label.size = 4,
                 max.cutoff = quantile(CD8T.obj@meta.data[,i],0.95),
                 min.cutoff = quantile(CD8T.obj@meta.data[,i],0.05),
                 raster = T,order=F) +
    scale_color_gradientn(colors=pals::jet()[c(1:11,15:25)],na.value = "#e6e6e6")+
    theme_umap()
  gp.list[[length(gp.list)+1]]<-p
}

do.call(ggpubr::ggarrange,c(gp.list,ncol = 1,nrow = 1)) -> combined.gp
pdf(paste0('CD8T.functional.score.umap.',Sys.Date(),'.pdf'),height=8,width = 8)
print(combined.gp)
dev.off()

###################### Figure S3C CD8 HSPA1A, HSPA1B expression ######################
CD8T.obj<-readRDS('MM.TME.235samples.CD8T.filter.scs.rds')

Idents(CD8T.obj)<-CD8T.obj$cell.status
p<-VlnPlot(CD8T.obj,features = c('HSPA1A','HSPA1B'),pt.size = 0,ncol = 2)&
  stat_summary(fun.y=median,geom='point',
               position = position_dodge(width = 0.9)) &
  theme(axis.text.x = element_text(angle = 90),
        #axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),#axis.ticks=element_blank(),
        #strip.text.x = element_text(size=12, color="black", face="bold"),
        strip.background = element_blank(),strip.text = element_text(angle = 0),
        #strip.text.y = element_text(angle = 0),
        text = element_text(size=15),
        #plot.title = element_blank(),#element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA, size=1))
ggsave('CD8T.HSPA1A_HSPA1B.vlnplot.pdf',useDingbats=F,p,width = 6,height=5)

###################### Figure S3D ER stress related genes bubble plot in CD8 ######################
CD8T.obj<-readRDS('MM.TME.235samples.CD8T.filter.scs.rds')

p<-DotPlot(CD8T.obj, features = c('ATF6','ERN1','EIF2AK3','HSPA5','MBTPS1','MBTPS2','XBP1','TRAF2','NFKB1','MAPK8','JUN','FOS','EIF2A','ATF4','EIF4EBP1','DDIT3','PPP1R15A','HIF1A','SLC2A1','LDHA','GFPT1'),group.by = 'cell.status',col.min = -1, col.max = 1) # ,'SCD','FASN','ACLY','ACACA','HMGCR','HMGCS1'
data4<-p$data[,c('id','features.plot','pct.exp','avg.exp.scaled')]

data4$group<-'Intracellular Signal Transduction'
data4$group[data4$features.plot%in%c('ATF6','ERN1','EIF2AK3','HSPA5')]<-'Sensors in ER'
data4$group[data4$features.plot%in%c('HIF1A','SLC2A1','LDHA','GFPT1')]<-'Response to hypoxia'

data4$group<-factor(data4$group,levels=c('Sensors in ER','Intracellular Signal Transduction','Response to hypoxia'))

plotx<-ggplot(data4, aes(x = features.plot,y = id)) +        ## global aes
  geom_point(aes(fill = avg.exp.scaled,size =pct.exp),color='black',shape=21)  +    ## geom_point for circle illusion
  scale_fill_gradientn(colours=c('#f7fcfd','#f7fcfd','#e0ecf4','#bfd3e6','#9ebcda','#8c96c6','#8c6bb1','#88419d','#810f7c','#4d004b')) +       ## color of the corresponding aes
  #scale_fill_viridis()+
  #scale_size(range = c(0,6), limits = c(0, 100), breaks = c(0,20,40,60,80,100)) +             ## to tune the size of circles
  facet_grid(cols=vars(group),scales='free',space='free')+
  labs(x='',y='')+
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "grey90",size=0.2), #panel.grid.minor = element_blank(),
        strip.background = element_blank(),strip.text.y = element_text(angle = 0),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle = 90,vjust=0.2,hjust=0.95) )

ggsave(paste0('CD8T.cell.status.Stress_genes.BubblePlot.',Sys.Date(),'.pdf'),plotx,height=3.5,width=7)


###################### Figure S3E CD8 GOBP enrichment ######################
CD8T.obj<-readRDS('MM.TME.235samples.CD8T.filter.scs.rds')

library('clusterProfiler')
library('msigdbr')
library('org.Hs.eg.db')
library('GOSemSim')
library('enrichplot')

GOBP_sig <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = 'GO:BP') %>% 
  dplyr::select(gs_name, gene_symbol)

m_t2g<-GOBP_sig

cluster.markers<-FindMarkers(CD8T.obj,min.pct = 0.25,ident.1 = c('CD8T_Teff_Fos/Jun','CD8T_Teff_TGFB1'),logfc.threshold = 0,only.pos =F)
cluster.markers<-cluster.markers%>%arrange(desc(avg_log2FC))

geneList<-cluster.markers$avg_log2FC
names(geneList)<-rownames(cluster.markers)
geneList<-geneList[!is.na(geneList)]

em<-GSEA(geneList,TERM2GENE=m_t2g,minGSSize = 5,eps=0, pvalueCutoff = 1)

i='CD8T_Tstr'

out<-as.data.frame(em)
out$group<-ifelse(out$NES>0,paste0('Enriched in ',i),paste0('Low in ',i))
out$`-log10(FDR)`<- -log10(out$qvalues)
out$group<-factor(out$group,levels=c(paste0('Enriched in ',i),paste0('Low in ',i)))
out<-out[out$`-log10(FDR)`> -log10(0.01),]

out$`-log10(FDR)`[out$group==paste0('Low in ',i)]<- -out$`-log10(FDR)`[out$group==paste0('Low in ',i)]

out<-out%>%group_by(group)%>%arrange(desc(`-log10(FDR)`), .by_group = T)
out$ID<-gsub('GOBP_','',out$ID)
out$ID<-tolower(out$ID)
out$ID<-factor(out$ID,levels=rev(out$ID))

p1<-ggplot(out, aes(`-log10(FDR)`,ID)) +
  geom_col(aes(fill=group)) +
  #scale_fill_manual(values = c(paste0('Enriched in ',i)='red',paste0('Low in ',i)='blue'))+
  labs(x="-log10(FDR)", y="")+
  theme(axis.text.x = element_text(size = 12,angle=0,vjust=0.2,hjust=0.5),
        strip.text.x = element_text(size=12, color="black", face="bold"),
        text = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))

ggsave(paste0(i,'.vs.other_CD8T','.GOBP.GSEA.pdf'),p1,width = 20,height = 12)

###################### Figure S3F CD8 TCR overlap ######################
CD8T.obj<-readRDS('MM.TME.235samples.CD8T.filter.scs.rds')

CD8T_meta<-CD8T.obj@meta.data
tcr_overlap<-matrix(0,length(levels(CD8T_meta$cell.status)),length(levels(CD8T_meta$cell.status)))

for(i in 1:nrow(tcr_overlap)){
  for(j in 1:ncol(tcr_overlap)){
    if(i==j){tcr_overlap[i,j]<-0}else{
      tcr_overlap[i,j]<-length(intersect(CD8T_meta$TCR_clonotype_raw[CD8T_meta$cell.status==levels(CD8T_meta$cell.status)[i]],CD8T_meta$TCR_clonotype_raw[CD8T_meta$cell.status==levels(CD8T_meta$cell.status)[j]]))/length(union(CD8T_meta$TCR_clonotype_raw[CD8T_meta$cell.status==levels(CD8T_meta$cell.status)[i]],CD8T_meta$TCR_clonotype_raw[CD8T_meta$cell.status==levels(CD8T_meta$cell.status)[j]]))
    }
  }
}

rownames(tcr_overlap)<-levels(CD8T_meta$cell.status)
colnames(tcr_overlap)<-levels(CD8T_meta$cell.status)

#tcr_overlap[tcr_overlap>0.2]<-0.2

plotx<-pheatmap::pheatmap(t(tcr_overlap),main='TCR overlap',
                          color=colorRampPalette(rev(c('#d53e4f',
                                                       '#f46d43',
                                                       '#fdae61',
                                                       '#fee08b',
                                                       '#ffffbf',
                                                       '#e6f598',
                                                       '#abdda4',
                                                       '#66c2a5',
                                                       '#3288bd')))(100),
                          show_colnames =T,
                          #display_numbers=T,
                          number_color = 'black',
                          fontsize_number = 10,
                          #angle_col = 0,
                          show_rownames =T,
                          cluster_rows = F,
                          cluster_cols = F
)
ggsave('CD8T.subcluster.TCR.overlap.pdf',plotx,width=8,height=8)

###################### Figure S3G NK cells checkpoints ######################
NK.obj<-readRDS('MM.TME.235samples.NK.filter.scs.rds')
Idents(NK.obj)<-factor(NK.obj$cell.status,levels=c('NK_GZMK_SELL','NK_GZMK_TOX2','NK_GZMK','NK_CD16high','NK_CD16high_MHCII','NK_JUND','NK_Tisg'))
NK.obj$cell.status<-Idents(NK.obj)

NK.checkpoint.gene=read_excel('~/data/gene_sets/NK-receptor-checkpoint.xlsx',sheet="Sheet1")

p<-DotPlot(NK.obj,features = c(NK.checkpoint.gene$Symbol,'ITGAL','TMIGD2','CD2','KLRC1','IL12RB1','IL12RB2','IL18R1','IL21R','IL15RA','IL2RB','TGFBR1','TGFBR2','ACVR1','ADORA2A','NCAM1','GZMB','NKG7','PRF1'),group.by = 'cell.status') # LFA1, CD28H, CD2, NKG2A
data<-p$data

data<-left_join(data,NK.checkpoint.gene[,c('Symbol','Function')],by=c('features.plot'='Symbol'))
data$Function[data$features.plot=='CD244']<-'Activating' # 2B4, https://www.nature.com/articles/s41568-020-0272-z/figures/1
data$Function[data$features.plot=='CD226']<-'Activating' # DNAM1, https://www.nature.com/articles/s41568-020-0272-z/figures/1
data$Function[data$features.plot%in%c('ITGAL','TMIGD2','CD2')]<-'Activating'
data$Function[data$features.plot%in%c('NCAM1')]<-'Activating'
data$Function[data$features.plot%in%c('KLRC1')]<-'Inhibitory'
data$Function[data$features.plot%in%c('IL12RB1','IL12RB2','IL18R1','IL21R','IL15RA','IL2RB')]<-'Stimulatory_Receptors'
data$Function[data$features.plot%in%c('TGFBR1','TGFBR2','ACVR1','ADORA2A')]<-'Suppressive_Receptors'
data$Function[data$features.plot%in%c('GZMB','NKG7','PRF1')]<-'Cytotoxic'

data$id<-factor(data$id,levels=rev(levels(NK.obj)))

plotx<-ggplot(data, aes(x = features.plot,y = id)) +        ## global aes
  geom_point(aes(fill = avg.exp.scaled,size =pct.exp),color='black',shape=21)  +    ## geom_point for circle illusion
  #scale_fill_gradientn(colours=rev(color),limits=c(0,max(data$avg.exp)))+       ## color of the corresponding aes
  scale_fill_gradientn(colours=c('#f7fcfd','#f7fcfd','#e0ecf4','#bfd3e6','#9ebcda','#8c96c6','#8c6bb1','#88419d','#810f7c','#4d004b')) + 
  scale_size(range = c(0,6), limits = c(0, 100), breaks = c(0,20,40,60,80,100)) +             ## to tune the size of circles
  labs(x='',y='')+
  facet_grid(cols=vars(Function),scales='free',space='free')+
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "grey90",size=0.2), #panel.grid.minor = element_blank(),
        strip.background = element_blank(),strip.text = element_text(angle = 45),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle = 90,vjust=0.2,hjust=0.95) )
ggsave('NK.checkpoint.gene.bubbleplot.pdf',plotx,height=3.5,width=15)

###################### Figure S3H IGSF8 expression in myeloma cell ######################

tumor <- readRDS('MM.TME.235samples.Tumor.filter.scs.rds')

tumor_sub <- tumor %>%
  subset(cell.type == "cPC") %>%
  subset(NewID %in% names(which(table(NewID) >= 50)))

Idents(tumor_sub)<-tumor_sub$NewID

meta <- FetchData(tumor_sub, vars = c("NewID", "Diagnosis")) %>%
  distinct()

expr_df <- AverageExpression(
  tumor_sub,
  features = "IGSF8"
)$RNA %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("NewID")

colnames(expr_df)[2] <- "expression"

plot_df <- expr_df %>%
  left_join(meta, by = "NewID") %>%
  mutate(
    Diagnosis = factor(Diagnosis, levels = c('MGUS','SMM','NDMM','RRMM'))
  )

P<-ggplot(plot_df, aes(x=Diagnosis, y=expression) ) + 
  geom_boxplot(show.legend = T,width = 0.75, aes(fill=Diagnosis), outlier.shape=NA) + 
  geom_jitter(position=position_jitter(width=0.1, height=0), size=1.5) +
  ggsignif::geom_signif(test="wilcox.test", comparisons = combn(levels(factor(cc$Diagnosis)),2, simplify = F),
                        map_signif_level = F,   # c("***"=0.001, "**"=0.01, "*"=0.05)       
                        vjust=-0.15,
                        textsize=6,
                        size=0.75,
                        step_increase = 0.15)+
  scale_fill_manual(values=c('MGUS'='#1c75bc','SMM'='#fbb040','NDMM'='#d01c8b', 'RRMM'='#8856a7'))+
  #facet_wrap(.~Var1,ncol=7,scales='free') +
  scale_y_continuous(expand = expansion(mult = c(0.15)))+
  theme(strip.text.x = element_text(size=14, color="black", face="bold"),
        axis.text.x = element_text(angle = 90,size = 14,vjust=0.2,hjust=0.95,colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        text = element_text(size=20)) +
  labs(x="", y="Expression") 
ggsave('Tumor.IGSF8.compare.pdf',P,width=4.5,height=4.5)


###################### Figure S4A Myeloid subsets signature score heatmap ######################
TME<-readRDS('MM.TME.235samples.TME.no_Eryth.filter.scs.rds')

Myeloid<-subset(TME,cell.type%in%c('CD14_Mono','CD14CD16_Mono','CD16_Mono','Macrophage','cDC','pDC',
                                          'Neutrophil'))

Idents(Myeloid)<-factor(Myeloid$cell.status)

M1<-c('IL23','TNF','CXCL9','CXCL10','CXCL11','CD86','IL1A','IL1B',
      'IL6','CCL5','IRF5','IRF1', 'CD40','IDO1','KYNU','CCR7' )
M2<-c('IL4R','CCL4','CCL13','CCL20','CCL17', 
      'CCL18','CCL22','CCL24','LYVE1','VEGFA','VEGFB','VEGFC', 
      'VEGFD','EGF','CTSA','CTSB','CTSC','CTSD','TGFB1','TGFB2', 
      'TGFB3','MMP14','MMP19','MMP9','CLEC7A','WNT7B','FASL', 
      'TNFSF12','TNFSF8','CD276','VTCN1','MSR1','FN1','IRF4' )
Angiogenesis<-c('CCND2','CCNE1','CD44','CXCR4','E2F3','EDN1','EZH2','FGF18',
                'FGFR1','FYN','HEY1','ITGAV','JAG1','JAG2','MMP9','NOTCH1', 
                'PDGFA','PTK2','SPP1','STC1','TNFAIP6','TYMP','VAV2','VCAN','VEGFA')
Phagocytosis<-c('MRC1','CD163','MERTK','C1QB' )
MDSC<-unique(c('STAT1','STAT3','STAT6','NFKB1','S100A8','S100A9','ANXA1',
               'LYZ','CXCL1','CXCL2','CXCL8','LILRA3','TREM1',
               'PTGS2','ARG1','ARG2','WFDC17','OLR1','CD84','NOS2','IL10'))


Myeloid<-AddModuleScore(Myeloid,features = list(M1,M2,Angiogenesis,Phagocytosis,MDSC),name=c('M1_score','M2_score','Angio_score','phago_score','MDSC_score') )

aa<-FetchData(Myeloid,vars = c('cell.status','M1_score1','M2_score2','Angio_score3','phago_score4','MDSC_score5'))

bb<-reshape2::melt(aa)

bb$value<-bb$value + rnorm(n = length(x =bb$value)) / 100000

bb$cell.type<-as.character(bb$cell.status)
bb$cell.type[grepl('DC',bb$cell.status)]<-'DC'
bb$cell.type[grepl('CD14_Mono',bb$cell.status)]<-'CD14_Mono'
bb$cell.type[grepl('CD14CD16_Mono',bb$cell.status)]<-'CD14CD16_Mono'
bb$cell.type[grepl('^CD16_Mono',bb$cell.status)]<-'CD16_Mono'
bb$cell.type[grepl('Mac',bb$cell.status)]<-'Macro'
bb$cell.type[grepl('Neutro',bb$cell.status)]<-'Neutro'

cc<-bb%>%group_by(cell.status,variable)%>%summarise(median=median(value))%>%as.data.frame()
dd<-reshape2::dcast(cc,cell.status~variable)
rownames(dd)<-dd$cell.status
dd<-dd[,-1]

dd<-scale(dd)
dd[dd>1.5]<- 1.5
dd[dd< -1.5]<- -1.5

library(pheatmap);
pdf('Myeloid.signature.score.heatmap.pdf',height = 6,width = 4)
pheatmap::pheatmap(dd,fontsize_col = 10,
                   #scale ='row',
                   border_color = NA,
                   treeheight_row =15,
                   treeheight_col =30,
                   cluster_rows = F,
                   show_rownames = T,
                   show_colnames = T,
                   cluster_cols = T
);
dev.off()

